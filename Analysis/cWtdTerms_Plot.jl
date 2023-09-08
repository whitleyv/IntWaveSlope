using Statistics
using Printf
using Oceananigans
using JLD2
using CairoMakie

#ENV["GKSwstype"] = "nul" # if on remote HPC

#path_name = "/glade/scratch/whitleyv/NewAdvection/Parameters/VaryU03C/wPE/"

apath = "Analysis/Plots/"
dpath =  "Data/"

include("parameters.jl")

setname = "U300N100Lz100g100"

# how many times were saved?
#zlength = 250
ylength = 1000 # roughly out to y = 4000 (this should be enough, at least for .2 ??)
xlength = 38
tlength = 161 

const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best

pm = getproperty(SimParams(), Symbol(setname))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
m = -π/pm.Lz,
l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
Tf = 2*π/pm.f, 
Tσ = 2*π/pm.σ))

zSlopeSameˢ = -pm.Tanαˢ * ySlopeSameˢ
ySlopeSame = zSlopeSameˢ / -pm.Tanα
ΔySlopeSame = ySlopeSameˢ - ySlopeSame
    
@inline heaviside(X) = ifelse(X <0, 0.0, 1.0)
# exponential gaussian for curved corner
@inline expcurve(y, ystar, smo) = -pm.Lzˢ + pm.Lzˢ * exp(-(y-ystar)^2/(2*smo^2))
# planar slope line
@inline linslope(y) = -pm.Tanα*y
# combining the 2 with heaviside split at ySlopeSame
@inline curvedslope(y) = linslope(y) + (-linslope(y) + expcurve(y, gausT_center-ΔySlopeSame, gausT_width)) * heaviside(y-ySlopeSame)

filepath = dpath * "cWtd_Gauss_terms.jld2"
filepath2 = dpath * "cWtd_Gauss_terms2.jld2"

f1 = jldopen(filepath)
f2 = jldopen(filepath2)

fkeys = keys(f1)
fkeys2 = keys(f2)

##### data onslope gaussian t = 0
∫b∂Cg_∂t_Cg_sum = f1[fkeys[1]];
∫Cg∂b_∂t_Cgsum = f1[fkeys[2]];
∂b̄_∂t_Cg = f1[fkeys[3]];
b̄_Cg = f1[fkeys[10]];
dt_times = f1[fkeys[13]];
##### data offslope gaussian t = 0
∫b∂Cgr_∂t_Cgr_sum =  f1[fkeys[4]];
∫Cgr∂b_∂t_Cgrsum =  f1[fkeys[5]];
∂b̄_∂t_Cgr =  f1[fkeys[6]];
b̄_Cgr = f1[fkeys[11]];
##### data onslope gaussian t = 4 Ts
∫b∂Cgl_∂t_Cgl_sum =  f1[fkeys[7]];
∫Cgl∂b_∂t_Cglsum = f1[fkeys[8]];
∂b̄_∂t_Cgl =   f1[fkeys[9]];
b̄_Cgl = f1[fkeys[12]];
dt_times_spin = f1[fkeys[14]];
##### data onslope gaussian t = 5 Ts
∫b∂Cgl_∂t_Cgl2_sum =  f2[fkeys2[1]];
∫Cgl∂b_∂t_Cgl2sum = f2[fkeys2[2]];
∂b̄_∂t_Cgl2 =  f2[fkeys2[3]];
b̄_Cgl2 = f2[fkeys2[4]];
dt_times_spin2 = f2[fkeys2[5]];

######################
#                  dtb PLOT FOR ALL GAUSSIAN
######################
f = Figure(resolution = (1000, 700), fontsize=26)
ga = f[1, 1] = GridLayout() 
ax1 = Axis(ga[1, 1], ylabel = "∂b̄/∂t [ms⁻³]", xlabel = "Wave Periods [Tσ]")
ax1.xticks =  2:2:10
#ax1.yticks =  (-5e-8:5e-8:5e-8, ["-5×10⁻⁸", "0", "5×10⁻⁸"])

limits!(ax1, 0, 7, -1e-7,   1e-7)

dbr = lines!(ax1, dt_times, ∂b̄_∂t_Cgr, color = :gray32, linewidth = 5)
db = lines!(ax1, dt_times, ∂b̄_∂t_Cg, color = :firebrick2, linewidth = 5)
dbl = lines!(ax1, dt_times_spin, ∂b̄_∂t_Cgl, color = :dodgerblue2, linewidth = 5)
dbl2 = lines!(ax1, dt_times_spin2, ∂b̄_∂t_Cgl2, color = :darkgreen, linewidth = 5)
axislegend(ax1, [db, dbl, dbl2, dbr], ["On Slope, t = 0", "On Slope, t = 4 Tσ",  "On Slope, t = 5 Tσ", "Off Slope, t = 0"],
position = :lt)

savename = "Gauss_cwtd_dbdtonly_" * setname

save(apath * savename * ".png", f)

######################
#                  b CENTROID PLOT FOR ALL GAUSSIAN
######################
f = Figure(resolution = (1000, 700), fontsize=26)
ga = f[1, 1] = GridLayout() 
ax1 = Axis(ga[1, 1], ylabel = "Δb̄ [ms⁻²]", xlabel = "Wave Periods [Tσ]")
ax1.xticks =  2:2:10
#ax1.yticks =  (-5e-8:5e-8:5e-8, ["-5×10⁻⁸", "0", "5×10⁻⁸"])

limits!(ax1, 0, 11, -2.3e-4, 2.3e-4)

dbr = lines!(ax1, dt_times, b̄_Cgr[2:end] .- b̄_Cgr[1] , color = :gray32, linewidth = 5)
db = lines!(ax1, dt_times, b̄_Cg[2:end] .- b̄_Cg[1] , color = :firebrick2, linewidth = 5)
dbl = lines!(ax1, dt_times_spin, b̄_Cgl[2:end] .- b̄_Cgl[1] , color = :dodgerblue2, linewidth = 5)
dbl2 = lines!(ax1, dt_times_spin2, b̄_Cgl2[2:end] .- b̄_Cgl2[1] , color = :darkgreen, linewidth = 5)

axislegend(ax1, [db, dbl, dbl2, dbr], ["On Slope, t = 0", "On Slope, t = 4 Tσ", "On Slope, t = 5 Tσ", "Off Slope, t = 0"],
position = :lt)

savename = "Gauss_cwtd_bbaronly_" * setname

save(apath * savename * ".png", f)

######################
#                  b CENTROID PLOT FOR ALL GAUSSIAN at actual location
######################
f = Figure(resolution = (1000, 700), fontsize=26)
ga = f[1, 1] = GridLayout() 
ax1 = Axis(ga[1, 1], ylabel = "Δb̄ [ms⁻²]", xlabel = "Wave Periods [Tσ]")
ax1.xticks =  2:2:10
#ax1.yticks =  (-5e-8:5e-8:5e-8, ["-5×10⁻⁸", "0", "5×10⁻⁸"])

#limits!(ax1, 0, 11, -3.66e-3, -2.9e-3)

tims = (0:600:(pm.Tσ*11))./pm.Tσ
db = lines!(ax1, tims, b̄_Cg , color = :firebrick2, linewidth = 5)
dbl = lines!(ax1, tims[59:end], b̄_Cgl , color = :dodgerblue2, linewidth = 5)
dbl2 = lines!(ax1, tims[73:end], b̄_Cgl2, color = :darkgreen, linewidth = 5)

axislegend(ax1, [db, dbl, dbl2], ["On Slope, t = 0", "On Slope, t = 4 Tσ", "On Slope, t = 5 Tσ"],
position = :lb)

savename = "Gauss_cwtd_bbaronly_true_" * setname

save(apath * savename * ".png", f)

######################
#                  dtb COMPONENTS PLOT FOR ALL GAUSSIAN
######################

f = Figure(resolution = (1000, 1000), fontsize=26)
ga = f[1, 1] = GridLayout() 
ax1 = Axis(ga[1, 1], ylabel = "∫ b∂c/∂t dV [ms⁻³]")

ax2 = Axis(ga[2, 1], ylabel = "∫ c∂b/∂t dV [ms⁻³]", xlabel = "Wave Periods [Tσ]")
ax2.xticks =  2:2:10

limits!(ax1, 0, 7, -6e-7, 2e-6)
limits!(ax2, 0, 7, -1.5e-6, 6e-7)

hidexdecorations!(ax1, grid = false)

dbr = lines!(ax1, dt_times, ∫b∂Cgr_∂t_Cgr_sum, color = :gray32, linewidth = 5)
db = lines!(ax1, dt_times, ∫b∂Cg_∂t_Cg_sum, color = :firebrick2, linewidth = 5)
dbl = lines!(ax1, dt_times_spin, ∫b∂Cgl_∂t_Cgl_sum, color = :dodgerblue2, linewidth = 5)

lines!(ax2, dt_times, ∫Cgr∂b_∂t_Cgrsum, color = :gray32, linewidth = 5)
lines!(ax2, dt_times, ∫Cg∂b_∂t_Cgsum, color = :firebrick2, linewidth = 5)
lines!(ax2, dt_times_spin, ∫Cgl∂b_∂t_Cglsum, color = :dodgerblue2, linewidth = 5)

axislegend(ax1, [db, dbl, dbr], ["On Slope, t = 0", "On Slope, t = 4 Tσ", "Off Slope, t = 0"],
position = :rt)

savename = "Gauss_cwtd_dbdtterms_" * setname

save(apath * savename * ".png", f)

######################
#                  dtb MAGNITUDE COMPONENTS PLOT FOR SLOPE GAUSSIAN ON SAME PLOT
######################
f = Figure(resolution = (1000, 700), fontsize=26)
ga = f[1, 1] = GridLayout() 
ax1 = Axis(ga[1, 1], ylabel = "|∂b̄/∂t| [ms⁻³]", xlabel = "Wave Periods [Tσ]")
ax1.xticks =  2:2:10

limits!(ax1, 0, 7, -6e-7, 2e-6)

db = lines!(ax1, dt_times, abs.(∫b∂Cg_∂t_Cg_sum), color = :firebrick2, linewidth = 5)
db2 = lines!(ax1, dt_times, abs.(∫Cg∂b_∂t_Cgsum), color = :firebrick2, linewidth = 5, linestyle = :dot)

dbl = lines!(ax1, dt_times_spin, abs.(∫b∂Cgl_∂t_Cgl_sum), color = :dodgerblue2, linewidth = 5)
dbl2 = lines!(ax1, dt_times_spin,  abs.(∫Cgl∂b_∂t_Cglsum), color = :dodgerblue2, linewidth = 5, linestyle = :dot)

axislegend(ax1, [db, dbl], ["t = 0", "t = 4 Tσ"],
position = :rt)

savename = "Gauss_cwtd_dbdttermstogether_" * setname

save(apath * savename * ".png", f)

######################
#                  dtb SUM of COMPONENTS PLOT FOR SLOPE GAUSSIAN ON SAME PLOT
######################
f = Figure(resolution = (1000, 700), fontsize=26)
ga = f[1, 1] = GridLayout() 
ax1 = Axis(ga[1, 1], ylabel = "∫ b∂c/∂t dV (-), -∫ c∂b/∂t dV (⋅⋅⋅) [ms⁻³]", xlabel = "Wave Periods [Tσ]")
ax1.xticks =  2:2:10

limits!(ax1, 0, 7, -1e-6, 1e-6)

db = lines!(ax1, dt_times, ∫b∂Cg_∂t_Cg_sum .+ ∫Cg∂b_∂t_Cgsum, color = :firebrick2, linewidth = 5)
db2 = lines!(ax1, dt_times, ∂b̄_∂t_Cg, color = (:firebrick2, 0.5), linewidth = 5)

dbl = lines!(ax1, dt_times_spin, ∫b∂Cgl_∂t_Cgl_sum .+ ∫Cgl∂b_∂t_Cglsum, color = :dodgerblue2, linewidth = 5)
dbl2 = lines!(ax1, dt_times_spin, ∂b̄_∂t_Cgl, color = (:dodgerblue2, 0.5), linewidth = 5)

axislegend(ax1, [db, db2, dbl, dbl2], ["t = 0, ∑", "t = 0, ∂ₜb̄", "t = 4 Tσ,  ∑", "t = 4 Tσ, ∂ₜb̄"],
position = :rt)

savename = "Gauss_cwtd_dbdttermssum_" * setname

save(apath * savename * ".png", f)

######################
#                  GAUSSIAN, RAISED, and LATE RELEASE WTD VELOCITY TERMS PLOT
#######################

f = Figure(resolution = (1200, 1000), fontsize=26)
ga = f[1, 1] = GridLayout() 

ax1 = Axis(ga[1, 1], ylabel = "∂b̄/∂t [ms⁻³]")
ax2 = Axis(ga[2, 1], ylabel = "∂b̄/∂t [ms⁻³]")
ax3 = Axis(ga[3, 1], ylabel = "∂b̄/∂t [ms⁻³]", xlabel = "Wave Periods [Tσ]")
ax3.xticks =  2:2:10

ax1.yticks =  (-5e-7:5e-7:5e-7, ["-5×10⁻⁷", "0", "5×10⁻⁷"])
ax2.yticks =   (-5e-7:5e-7:5e-7, ["-5×10⁻⁷", "0", "5×10⁻⁷"])
ax3.yticks =  (-5e-7:5e-7:5e-7, ["-5×10⁻⁷", "0", "5×10⁻⁷"])

limits!(ax1, 0, 11, -1e-6,   1e-6)
limits!(ax2, 0, 11, -1e-6,   1e-6)
limits!(ax3, 0, 11, -1e-6,   1e-6)

hidexdecorations!(ax1)
hidexdecorations!(ax2)

dbb = lines!(ax1, dt_times, ∂b̄_∂t_Cg, color = :gray32, linewidth = 5)
cdb = lines!(ax1, dt_times, ∫Cg∂b_∂t_Cgsum, color = :firebrick2, linewidth = 5)
bdc = lines!(ax1, dt_times, ∫b∂Cg_∂t_Cg_sum, color = :dodgerblue2, linewidth = 5)

lines!(ax2, dt_times_spin, ∂b̄_∂t_Cgl, color = :gray32, linewidth = 5)
lines!(ax2, dt_times_spin, ∫Cgl∂b_∂t_Cglsum, color = :firebrick2, linewidth = 5)
lines!(ax2, dt_times_spin, ∫b∂Cgl_∂t_Cgl_sum, color = :dodgerblue2, linewidth = 5)

lines!(ax3, dt_times, ∂b̄_∂t_Cgr, color = :gray32, linewidth = 5)
lines!(ax3, dt_times, ∫Cgr∂b_∂t_Cgrsum, color = :firebrick2, linewidth = 5)
lines!(ax3, dt_times, ∫b∂Cgr_∂t_Cgr_sum, color = :dodgerblue2, linewidth = 5)

Legend(ga[1, 1, Top()], [dbb, cdb, bdc], ["∂b̄/∂t", "∫c ∂b/∂t dV", "∫b ∂c/∂t dV"],
                tellheight = false, tellwidth = false,
                margin = (10, 10, 10, 30), framevisible = false, patchlabelgap = 7,
                halign = :center, valign = :bottom, orientation = :horizontal)

Label(ga[1, 1, Top()],"On Slope, Release t = 0 Tσ", valign = :top,
                font = :bold, fontsize = 25,
                padding = (0, 0, 10, 0))
                
Label(ga[2, 1, Top()],"On Slope, Release t = 4 Tσ", valign = :top,
                font = :bold, fontsize = 25,
                padding = (0, 0, 10, 0))
    
Label(ga[3, 1, Top()], "Off Slope, Release t = 0 Tσ", valign = :top,
                font = :bold, fontsize = 25,
                padding = (0, 0, 10, 0))

savename = "Gauss_cwtd_dbdt_" * setname

save(apath * savename * ".png", f)
