using Plots
using Measures
using Statistics
using Printf
using Oceananigans

using CairoMakie
using ArgParse

ENV["GKSwstype"] = "nul" # if on remote HPC

###########-------- SIMULATION PARAMETERS ----------------#############

# curved topo parameters
const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best


function parse_commandline()
  s = ArgParseSettings()
  @add_arg_table s begin
      "paramset"
          help = "sets which parameters to use"
          default = "U100N100Lz100g100"
      "--resScale"
     	  help = "scale of resolution from default dh = 4 and dz=2"
	        arg_type = Float64
	        default = 1.0  
     "path"
          help = "pathname to save data under"
          default = ""
  end
  return parse_args(s)
end

args=parse_commandline()

@info("command line args:")
for (arg,val) in args
  @info("   $arg => $val")
end

path_name = args["path"]
setname = args["paramset"]

resS = args["resScale"]
@info "Loading in parameters..."

include("../parameters.jl")
pm = getproperty(SimParams(), Symbol(setname))

dzr = pm.dz * resS
dhr = pm.dh * resS

pm = merge(pm, (; dzr=dzr, dhr=dhr, Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/dzr),
		            nx = round(Int, pm.Lx/dhr),
                m = -π/pm.Lzˢ,
                l = sqrt(((π/pm.Lzˢ)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                Tf = 2*π/pm.f, 
                Tσ = 2*π/pm.σ))

# other params for setting up the grid
z_start = pm.Lz - pm.Lzˢ

Sp_extra = ifelse(z_start>0, 250.0, 0.0)
Ly = pm.Lyˢ+Sp_extra
ny = round(Int,Ly/pm.dhr)
slope_end = pm.Lzˢ/pm.Tanα

pm = merge(pm, (;Ly=Ly,ny=ny, slope_end=slope_end, Sp_extra=Sp_extra))

# if slope is in different spot than usual, need to move the curved part too!
const zSlopeSameˢ = -pm.Tanαˢ * ySlopeSameˢ
ySlopeSame = zSlopeSameˢ / -pm.Tanα
ΔySlopeSame = ySlopeSameˢ - ySlopeSame

@inline heaviside(X) = ifelse(X <0, 0.0, 1.0)
# exponential gaussian for curved corner
@inline expcurve(y, ystar, smo) = -pm.Lzˢ + pm.Lzˢ * exp(-(y-ystar)^2/(2*smo^2))
# planar slope line
@inline linslope(y) = -pm.Tanα*y
# combining the 2 with heaviside split at ySlopeSame
@inline curvedslope(y) = linslope(y) + (-linslope(y) + expcurve(y, gausT_center-ΔySlopeSame, gausT_width)) * heaviside(y-ySlopeSame)

@info "Pulling Data..."
path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"
name_prefix = "IntWave_mp_" * setname
#name_prefix = "IntWave_" * setname * @sprintf("_R%0.0f", resS*100)
filepath = path_name * name_prefix * ".jld2"

v_timeseries = FieldTimeSeries(filepath,"v");
cs_timeseries = FieldTimeSeries(filepath,"Cs");
cg_timeseries = FieldTimeSeries(filepath,"Cg");
b_timeseries = FieldTimeSeries(filepath,"b");
ε_timeseries = FieldTimeSeries(filepath, "ϵ");
N2_timeseries = FieldTimeSeries(filepath, "N2");

xc, yc, zc = nodes(b_timeseries) #CCC
xv, yv, zv = nodes(v_timeseries) #CFC
xn, yn, zn = nodes(N2_timeseries) #CFC

land = curvedslope.(yc) 

kwargs = Dict(:linewidth => 0,
:colorbar => :true,
:xlims => (0, 4000),
:ylims => (-500, 0),
:size => (1150, 1200))

delta = pm.U₀/pm.Ñ

big_title = @sprintf("Internal Wave Breaking, U₀=%0.3f, N=%0.2f×10⁻³, δ=%0.1f", pm.U₀, 10^3*pm.Ñ, delta)

xlocat = 19

N21 = interior(N2_timeseries)[xlocat, :,:, 1] 

anim = @animate for (i, t) in enumerate(b_timeseries.times)
        @info "Plotting frame $i of $(length(b_timeseries.times))..."
    
        v = v_timeseries[i]
        b = b_timeseries[i]
        cg = cg_timeseries[i]
        cs = cs_timeseries[i]
        ε  = ε_timeseries[i]
        N2  = N2_timeseries[i]

        vi = interior(v)[xlocat, :,:]
        bi = interior(b)[xlocat, :,:]
        cgi = interior(cg)[xlocat, :,:]
        csi = interior(cs)[xlocat, :,:]
        εi = interior(ε)[xlocat, :,:]
        N2i = interior(N2)[xlocat, :,:] .- N21

        #MOMENTUM PLOTS
                 
        v_plot = heatmap(yv, zv, vi'; color = :balance, ylabel = "z",
                        guidefontsize = 14, titlefont=14, tickfont = 8, xticks=false,
                        clims = (-pm.U₀, pm.U₀), left_margin = 10.0mm, right_margin = 5.0mm, 
                        kwargs...)
                 plot!(yc, land, linecolor=:black, lw = 1, legend = false)
    
        # DYE
        cs_plot = heatmap(yc, zc, log10.(clamp.(csi, 10^(-8), 1))'; color = :thermal,
                        guidefontsize = 14, titlefont=14, tickfont = 8, ylabel = "z", xticks=false,
                        clims = (-4, 0), right_margin = 5.0mm, kwargs...)
                plot!(yc, land, linecolor=:black, lw = 1, legend = false)
    
        cg_plot = heatmap(yc, zc, log10.(clamp.(cgi, 10^(-8), 1))'; color = :thermal,
                guidefontsize = 14, titlefont=14, tickfont = 8, ylabel = "z", xlabel = "y",
                clims = (-4, 0), right_margin = 5.0mm, kwargs...)
        plot!(yc, land, linecolor=:black, lw = 1, legend = false)

        #last column
        b_plot = heatmap(yc, zc, bi'; color = :thermal, yticks=false, xticks=false,
                        guidefontsize = 14, titlefont=14, tickfont = 8,  
                        clims = (-0.007, 0), right_margin = 10.0mm,
                        kwargs...)
                    plot!(yc, land, linecolor=:black, lw = 1, legend = false)
                Plots.contour!(yc, zc, bi'; color=:black, lw=2)

        n_plot = heatmap(yn, zn, (N2i .* 1e5)'; color = :balance, yticks=false, xticks=false,
                guidefontsize = 14, titlefont=14, tickfont = 8, 
                clims = (-1.5, 1.5), right_margin = 10.0mm,
                kwargs...)
            plot!(yc, land, linecolor=:black, lw = 1, legend = false)
        
        e_plot = heatmap(yc, zc, log10.(clamp.(εi, 10^(-10), 1))'; color = :thermal, xlabel = "y", yticks=false,
            guidefontsize = 14, titlefont=14, tickfont = 8, 
            clims = (-9, -5), right_margin = 10.0mm,
            kwargs...)
        plot!(yc, land, linecolor=:black, lw = 1, legend = false)

        y = 1:5
    
        v_title = @sprintf("v velocity at tσ/2π = %.1f", t/pm.Tσ)
        cs_title = @sprintf("log₁₀Cs at tf/2π = %.1f", t/pm.Tf)
        cg_title = @sprintf("log₁₀Cg at t = %.1f hrs", t/3600)
        b_title = @sprintf("buoyancy at t = %.1f hrs", t/3600)
        n_title = @sprintf("N²'×10⁻⁵ at t = %.1f hrs", t/3600)
        e_title = @sprintf("ε at t = %.1f hrs", t/3600)
        
        BigT = Plots.scatter(y, marker=0, markeralpha=0, annotations=(3, y[3],
                            Plots.text(big_title)), axis=nothing, grid=false, leg=false,
                            foreground_color_subplot=colorant"white")
    
        data_plots = plot(v_plot, b_plot, cs_plot, n_plot, cg_plot, e_plot,
                    title = [v_title b_title cs_title n_title cg_title e_title],
                    layout = (3, 2))
        fin = Plots.plot(BigT, data_plots, layout=grid(2, 1, heights=[0.05,0.95]))
end

apath = path_name * "Analysis/"    
cont_name = apath * "Contours_" * name_prefix
mp4(anim, cont_name * ".mp4", fps = 8)


##############
#   CAIRO PLOT MOVIE
##############
n = Observable(1)
xlocat = 19

v = @lift interior(v_timeseries[$n],xlocat, :, :,);
b = @lift interior(b_timeseries[$n], xlocat, :, :);
cg = @lift log10.(clamp.(interior(cg_timeseries[$n],xlocat, :, :), 1e-8,1));
cs = @lift log10.(clamp.(interior(cs_timeseries[$n], xlocat, :, :), 1e-8,1));
ε  = @lift log10.(clamp.(interior(ε_timeseries[$n],xlocat, :, :), 1e-10,1));
N2  = @lift interior(N2_timeseries[$n], xlocat, :, :);

title = @lift "Internal Wave Breaking, t = " * string(round(b_timeseries.times[$n]/3600, digits=2)) * " hrs, Tσ = "*string(round(b_timeseries.times[$n]/pm.Tσ, digits=2))

f = Figure(resolution = (1150, 1200),fontsize=26) 

ga = f[1, 1] = GridLayout()
axv = Axis(ga[1, 1], ylabel = "z [m]")
axb = Axis(ga[1, 3])

axcg = Axis(ga[2, 1], ylabel = "z [m]")
axcs = Axis(ga[2, 3])

axe = Axis(ga[3, 1], ylabel = "z [m]", xlabel = "y [m]")
axn = Axis(ga[3, 3], xlabel = "y [m]")

axe.xticks = 500:1000:1500
axn.xticks = 500:1000:1500

axv.yticks = [-250, 0]
axcg.yticks = [-250, 0]
axe.yticks = [-250, 0]

limits!(axv, 0, 2000, -450, 0)
limits!(axcg, 0, 2000, -450, 0)
limits!(axcs, 0, 2000, -450, 0)
limits!(axb, 0, 2000, -450, 0)
limits!(axe, 0, 2000, -450, 0)
limits!(axn, 0, 2000, -450, 0)

hidedecorations!(axb)
hidedecorations!(axcs)
hidexdecorations!(axv)
hidexdecorations!(axcg)
hideydecorations!(axn)

global hmv = heatmap!(axv, yv, zv, v, colormap = :balance, colorrange = (-pm.U₀, pm.U₀))
lines!(axv, yc, land, color=:black, lw = 4)

global hmb = heatmap!(axb, yc, zc, b, colormap= :thermal, colorrange = (-0.007, 0),)
contour!(axb, yc, zc, b, color = :black, lw = 6, levels = -0.006:0.001:0, alpha = 0.5)
lines!(axb, yc, land, color=:black, lw = 4)

global hmcg = heatmap!(axcg, yc, zc, cg, colormap = :thermal, colorrange = (-4, 0))
lines!(axcg, yc, land, color=:black, lw = 4)

global hmcs = heatmap!(axcs, yc, zc, cs, colormap = :thermal, colorrange = (-4, 0))
lines!(axcs, yc, land, color=:black, lw = 4)

global hme = heatmap!(axe, yc, zc, ε, colormap = :thermal, colorrange = (-9, -5))
lines!(axe, yc, land, color=:black, lw = 4)

global hmn= heatmap!(axn, yn, zn, N2, colormap = :balance, colorrange = (-1.5e-5, 1.5e-5))
lines!(axn, yc, land, color=:black, lw = 4)


cb1 = Colorbar(ga[1,2], hmv, ticks = (-0.2:.1:0.2), size =35)
cb2 = Colorbar(ga[1,4], hmb, ticks = (-5e-3:2e-3:-1e-3, ["-0.005", "-0.003", "-0.001"] ), size = 35)
cb3 = Colorbar(ga[2,2], hmcg, ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ), size = 35)
cb4 = Colorbar(ga[2,4], hmcs, ticks = (-4:1:-1, ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"] ), size = 35)
cb5 = Colorbar(ga[3,2], hme, ticks = (-9:2:-4, ["10⁻⁹", "10⁻⁷", "10⁻⁵"] ), size = 35)
cb6 = Colorbar(ga[3,4], hmn, ticks = (-1e-5:1e-5:1e-5, ["-10⁻⁵", "0", "10⁻⁵"] ), size = 35)

colsize!(ga, 2, Relative(0.04))
colsize!(ga, 4, Relative(0.04))


colgap!(ga, 15)
rowgap!(ga, 5)

savename = "contours_m_" * setname
apath  = path_name * "Analysis/"

frames = 1:length(b_timeseries.times)
record(f, apath * savename * ".mp4", frames, framerate=8) do j
    msg = string("Plotting frame ", j, " of ", frames[end])
    print(msg * " \r")
    n[] = j
end
