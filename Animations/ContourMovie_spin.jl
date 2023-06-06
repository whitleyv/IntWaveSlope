using Plots
using Measures
using Statistics
using Printf
using Oceananigans

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

name_prefix1 = "IntWave_spinup_" * setname
name_prefix2 = "IntWave_postspin_" * setname 
filepath1 = path_name * name_prefix1 * ".jld2"
filepath2 = path_name * name_prefix2 * ".jld2"


v_timeseries1 = FieldTimeSeries(filepath1,"v");
b_timeseries1 = FieldTimeSeries(filepath1,"b");
e_timeseries1 = FieldTimeSeries(filepath1,"ϵ");
N2_timeseries1 = FieldTimeSeries(filepath1, "N2");

cs_timeseries = FieldTimeSeries(filepath2,"Cs");
cg_timeseries = FieldTimeSeries(filepath2,"Cg");
b_timeseries2 = FieldTimeSeries(filepath2,"b");
v_timeseries2 = FieldTimeSeries(filepath2,"v");
e_timeseries2 = FieldTimeSeries(filepath2, "ϵ");
N2_timeseries2 = FieldTimeSeries(filepath2, "N2");

xc, yc, zc = nodes(b_timeseries1) #CCC
xv, yv, zv = nodes(v_timeseries1) #CFC
xn, yn, zn = nodes(N2_timeseries1) #CFC

land = curvedslope.(yc) 

kwargs = Dict(:linewidth => 0,
:colorbar => :true,
:xlims => (0, 4000),
:ylims => (-500, 0),
:size => (1150, 1200))

delta = pm.U₀/pm.Ñ

big_title = @sprintf("Internal Wave Breaking, U₀=%0.3f, N=%0.2f×10⁻³, δ=%0.1f", pm.U₀, 10^3*pm.Ñ, delta)

xlocat = 19

N2ini = interior(N2_timeseries1)[xlocat, :,:, 1] 

spinlength = length(b_timeseries1.times)
fulltimes = vcat(b_timeseries1.times,b_timeseries2.times)
cini = zeros(size(N2ini))

anim = @animate for (i, t) in enumerate(fulltimes)
        @info "Plotting frame $i of $(length(fulltimes))..."
    
        if i <= spinlength
            vi = interior(v_timeseries1[i])[xlocat,:,:]
            bi = interior(b_timeseries1[i])[xlocat,:,:]
            ei  = interior(e_timeseries1[i])[xlocat,:,:]
            N2i  = interior(N2_timeseries1[i])[xlocat,:,:] .- N2ini
            cgi = cini
            csi = cini
        else
            vi = interior(v_timeseries2[i - spinlength])[xlocat,:,:]
            bi = interior(b_timeseries2[i - spinlength])[xlocat,:,:]
            ei  = interior(e_timeseries2[i - spinlength])[xlocat,:,:]
            N2i  = interior(N2_timeseries2[i - spinlength])[xlocat,:,:] .- N2ini
            cgi  = interior(cg_timeseries[i - spinlength])[xlocat,:,:]
            csi  = interior(cs_timeseries[i - spinlength])[xlocat,:,:]

        end


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
        
        e_plot = heatmap(yc, zc, log10.(clamp.(ei, 10^(-10), 1))'; color = :thermal, xlabel = "y", yticks=false,
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


