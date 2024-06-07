using JLD2
using CairoMakie
using Printf

dpath = "../Data/"
apath = "Plots/"

fname= "DCwtd_Gauss.jld2"
filepath = dpath * fname

wtd_file = jldopen(filepath, "r+")

fkeys = keys(wtd_file)

Cstrings = ["Cg", "Cgr", "Cgs"]
wl = 134
tl = 161

Wₜ_Wavg = zeros(3,wl)
Kₜ_Wavg = zeros(3,wl)
∂ₜb̄_Wavg = zeros(3,wl)
∇b̄_Wavg = zeros(3,wl+1)

Wₜ = zeros(3,tl-1)
∂ₜb̄ = zeros(3,tl-1)
∇b̄ = zeros(3,tl)

for (i, cs) in enumerate(Cstrings)
    Wₜ_Wavg[i, :] = wtd_file["Wₜ_Wavg_" * cs]
    Kₜ_Wavg[i, :] = wtd_file["Kₜ_Wavg_" * cs]
    ∂ₜb̄_Wavg[i, :] = wtd_file["∂ₜb̄_Wavg_" * cs]
    ∇b̄_Wavg[i, :] = wtd_file["∇b̄_Wavg_" * cs]

    Wₜ[i, :] = wtd_file["Wₜ_" * cs]
    ∂ₜb̄[i, :] = wtd_file["∂ₜb̄_" * cs]
    ∇b̄[i, :] = wtd_file["∇b̄_" * cs]
end

f = Figure(resolution = (900, 800), fontsize=26)
ga = f[1, 1] = GridLayout() # vhat
#gb = f[1, 2] = GridLayout()

ax1 = Axis(ga[1, 1], ylabel = "K [m²s⁻¹]",
 title = "Tracer Weighted Buoyancy Diffusivity") 
#ax1.yticks = (0:4e-4:1.2e-3 , ["0", "4×10⁻⁴", "8×10⁻⁴", "1.2×10⁻³"])
hidexdecorations!(ax1)
#limits!(ax1, 1, 10, -3e-4, 1.3e-3)

ax2 = Axis(ga[2, 1], ylabel = "W [ms⁻¹]", xlabel = "Tσ",
title =  "Tracer Weighted Buoyancy Velocity") 
#ax2.yticks = (-8e-5:4e-5:8e-5, ["-8×10⁻⁵", "-4×10⁻⁵", "0", "-4×10⁻⁵", "8×10⁻⁵"])
ax2.xticks = 2:2:10
#limits!(ax2, 1, 10, -1.2e-4, 1.2e-4)

yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2])
ax1.yticklabelspace = yspace
ax2.yticklabelspace = yspace

rowgap!(ga, 10)

wavs = (15 .+ 1)./15 : 1 ./15: 10 .- 1 ./15
labels = ["Center-Slope", "Center-Interior", "Up-SLope"]
linestyle = [:solid, :dot, :dash]
for i = 1:3
    lines!(ax1, wavs, Kₜ_Wavg[i,:], linewidth = 5, color = :dodgerblue2, 
            linestyle = linestyle[i], label = labels[i])

    lines!(ax2, wavs, Wₜ_Wavg[i,:], linewidth = 5, color = :firebrick2,linestyle = linestyle[i] )
end

leg = Legend(ga[1, 2], ax1)

#colsize!(f.layout, 2, Relative(0.05))

savename = "Paper_Cwtd_GaussKW_wavg" 

save(apath * savename * ".png", f)


f = Figure(resolution = (1000, 700), fontsize=26)
ga = f[1, 1] = GridLayout() # vhat

ax1 = Axis(ga[1, 1], ylabel = "", 
 title = "Diapycnal Velocity, W [ms⁻¹]") 
ax1.yticks = (-2e-4:1e-4:2e-4 , ["-2×10⁻⁴", "-1×10⁻⁴", "0", "1×10⁻⁴", "2×10⁻⁴"])
hidexdecorations!(ax1)
limits!(ax1, 1, 10, -2.3e-4, 2.3e-4)

ax2 = Axis(ga[2, 1], ylabel = "", xlabel = "Tσ",
title =  "Buoyancy Velocity, ∂ₜb̄ [ms⁻³]") 
#ax2.yticks = (-1e-8:5e-9:1e-8, ["-1×10⁻⁸", "-5×10⁻⁹", "0", "5×10⁻⁹", "1×10⁻⁸"])
ax2.xticks = 2:2:10
limits!(ax2, 1, 10, -1.2e-8, 1.2e-8)

ax3 = Axis(ga[2, 2], ylabel = "", xlabel = "Tσ",
title =  "Buoyancy Gradient, |∇b̄| [s⁻²]") 
ax3.yticks = (3e-5:1e-5:5e-5, ["3×10⁻⁵", "4×10⁻⁵", "5×10⁻⁵"])
ax3.xticks = 2:2:10
limits!(ax3, 1, 10, 2e-5, 6e-5)

rowgap!(ga, 10)

wavs = (15 .+ 1)./15 : 1 ./15: 10 .- 1 ./15
wavs2 = (15)./15 : 1 ./15: 10 .- 1 ./15
labels = ["On Slope, z = -250 m", "Off Slope, z = -136 m", "On Slope, z = -150 m"]
colors = [:dodgerblue2, :firebrick2, :gray30]
linestyle = [:solid, :solid, :solid]
for i = 3
    lines!(ax1, wavs, Wₜ_Wavg[i,:], linewidth = 5, color = colors[i], 
            linestyle = linestyle[i], label = labels[i])
    lines!(ax2, wavs, ∂ₜb̄_Wavg[i,:], linewidth = 5, color = colors[i],linestyle = linestyle[i] )
    lines!(ax3, wavs2, ∇b̄_Wavg[i,:], linewidth = 5, color = colors[i],linestyle = linestyle[i] )
end

for i = 1:2
    lines!(ax1, wavs, Wₜ_Wavg[i,:], linewidth = 5, color = colors[i], 
            linestyle = linestyle[i], label = labels[i])
    lines!(ax2, wavs, ∂ₜb̄_Wavg[i,:], linewidth = 5, color = colors[i],linestyle = linestyle[i] )
    lines!(ax3, wavs2, ∇b̄_Wavg[i,:], linewidth = 5, color = colors[i],linestyle = linestyle[i] )
end

leg = Legend(ga[1, 2], ax1, "Initial Centroid Location")
colsize!(ga, 2, Relative(0.5))
#leg.tellheight = true

Label(ga[1, 1:2, Top()], "Tracer-Weighted Buoyancy Velocity", valign = :bottom,
    font = :bold,
    padding = (0, 0, 40, 0))


#yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2, ax3])
#ax1.yticklabelspace = yspace
#ax2.yticklabelspace = yspace
#ax3.yticklabelspace = yspace

savename = "Paper_Cwtd_all3GaussW_wavg" 

save(apath * savename * ".png", f)


f = Figure(resolution = (900, 800), fontsize=26)
ga = f[1, 1] = GridLayout() # vhat

ax1 = Axis(ga[1, 1], ylabel = "W [ms⁻¹]", 
 title = "Tracer Weighted Diapycnal Diffusivity") 
#ax1.yticks = (0:4e-4:1.2e-3 , ["0", "4×10⁻⁴", "8×10⁻⁴", "1.2×10⁻³"])
hidexdecorations!(ax1)
limits!(ax1, 6, 10, -1e-4, 1e-4)

ax2 = Axis(ga[2, 1], ylabel = "∂ₜb̄ [ms⁻³]", 
title =  "Tracer Weighted Buoyancy Velocity") 
#ax2.yticks = (-8e-5:4e-5:8e-5, ["-8×10⁻⁵", "-4×10⁻⁵", "0", "-4×10⁻⁵", "8×10⁻⁵"])
limits!(ax2, 6, 10, -3e-9, 3e-9)

ax3 = Axis(ga[3, 1], ylabel = "|∇b̄| [s⁻²]", xlabel = "Tσ",
title =  "Tracer Weighted Buoyancy Gradient") 
#ax2.yticks = (-8e-5:4e-5:8e-5, ["-8×10⁻⁵", "-4×10⁻⁵", "0", "-4×10⁻⁵", "8×10⁻⁵"])
ax3.xticks = 2:2:10
limits!(ax3, 6, 10, 3e-5, 4e-5)

rowgap!(ga, 10)

wavs = (15 .+ 1)./15 : 1 ./15: 10 .- 1 ./15
wavs2 = (15)./15 : 1 ./15: 10 .- 1 ./15
labels = ["Center-Slope", "Center-Interior", "Up-SLope"]
colors = [:dodgerblue2, :firebrick2, :black]
linestyle = [:solid, :solid, :solid]
for i = 1:3
    lines!(ax1, wavs, Wₜ_Wavg[i,:], linewidth = 5, color = colors[i], 
            linestyle = linestyle[i], label = labels[i])
    lines!(ax2, wavs, ∂ₜb̄_Wavg[i,:], linewidth = 5, color = colors[i],linestyle = linestyle[i] )
    lines!(ax3, wavs2, ∇b̄_Wavg[i,:], linewidth = 5, color = colors[i],linestyle = linestyle[i] )
end

leg = Legend(ga[1, 2], ax1, title = "Initial Gaussian Centroid Location" )

#colsize!(f.layout, 2, Relative(0.05))

#yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2, ax3])
#ax1.yticklabelspace = yspace
#ax2.yticklabelspace = yspace
#ax3.yticklabelspace = yspace

savename = "Paper_Cwtd_GaussWLtd_wavg" 

save(apath * savename * ".png", f)


f = Figure(resolution = (900, 800), fontsize=26)
ga = f[1, 1] = GridLayout() # vhat
#gb = f[1, 2] = GridLayout()

ax1 = Axis(ga[1, 1], ylabel = "K [m²s⁻¹]",
 title = "Tracer Weighted Buoyancy Diffusivity") 
ax1.yticks = (0:0.002:0.006 , ["0", "2×10⁻³", "4×10⁻³", "6×10⁻³"])
hidexdecorations!(ax1)
limits!(ax1, 1, 10, 0.0, 8e-3)

ax2 = Axis(ga[2, 1], ylabel = "W [ms⁻¹]", xlabel = "Tσ",
title =  "Tracer Weighted Buoyancy Velocity") 
ax2.yticks = (-5e-5:5e-5:5e-5, ["-5×10⁻⁵", "0", "5×10⁻⁵"])
ax2.xticks = 2:2:10
limits!(ax2, 1, 10, -1.0e-4, 1.0e-4)

yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2])
ax1.yticklabelspace = yspace
ax2.yticklabelspace = yspace

rowgap!(ga, 10)

wavs = (15 .+ 1)./15 : 1 ./15: 10 .- 1 ./15

#labels = ["Center-Slope", "Up-SLope", "Center-Interior"]
i = 1   
lines!(ax1, wavs, Kₜ_Wavg[i,:], linewidth = 5, color = :dodgerblue2, )
lines!(ax2, wavs, Wₜ_Wavg[i,:], linewidth = 5, color = :firebrick2,)


#colsize!(f.layout, 2, Relative(0.05))

savename = "Paper_Cwtd_CentGaussKW_wavg" 

save(apath * savename * ".png", f)


f = Figure(resolution = (900, 800), fontsize=26)
ga = f[1, 1] = GridLayout() # vhat
#gb = f[1, 2] = GridLayout()

ax1 = Axis(ga[1, 1], ylabel = "K [m²s⁻¹]",
 title = "Tracer Weighted Buoyancy Diffusivity") 
ax1.yticks = (0:0.002:0.008 , ["0", "2×10⁻³", "4×10⁻³", "6×10⁻³", "8×10⁻³"])
hidexdecorations!(ax1)
limits!(ax1, 1, 10, 0.0, 1e-2)

ax2 = Axis(ga[2, 1], ylabel = "W [ms⁻¹]", xlabel = "Tσ",
title =  "Tracer Weighted Buoyancy Velocity") 
ax2.yticks = (-5e-5:5e-5:5e-5, ["-5×10⁻⁵", "0", "5×10⁻⁵"])
ax2.xticks = 2:2:10
limits!(ax2, 1, 10, -1.0e-4, 1.0e-4)

yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2])
ax1.yticklabelspace = yspace
ax2.yticklabelspace = yspace

rowgap!(ga, 10)

wavs = (15 .+ 1)./15 : 1 ./15: 10 .- 1 ./15

#labels = ["Center-Slope", "Up-SLope", "Center-Interior"]
i = 2    
lines!(ax1, wavs, Kₜ_Wavg[i,:], linewidth = 5, color = :dodgerblue2, )
lines!(ax2, wavs, Wₜ_Wavg[i,:], linewidth = 5, color = :firebrick2,)


#colsize!(f.layout, 2, Relative(0.05))

savename = "Paper_Cwtd_RaisedGaussKW_wavg" 

save(apath * savename * ".png", f)


f = Figure(resolution = (900, 800), fontsize=26)
ga = f[1, 1] = GridLayout() # vhat
#gb = f[1, 2] = GridLayout()

ax1 = Axis(ga[1, 1], ylabel = "K [m²s⁻¹]",
 title = "Tracer Weighted Buoyancy Diffusivity") 
ax1.yticks = (0.0:0.004:0.012 , ["0", "4×10⁻³","8×10⁻³", "1.2×10⁻²"])
hidexdecorations!(ax1)
limits!(ax1, 1, 10, -4e-3, 1.5e-2)

ax2 = Axis(ga[2, 1], ylabel = "W [ms⁻¹]", xlabel = "Tσ",
title =  "Tracer Weighted Buoyancy Velocity") 
ax2.yticks = (-1.5e-4:5e-5:1.5e-4, ["-1.5×10⁻⁴", "-1×10⁻⁴", "-5×10⁻⁵", "0", "5×10⁻⁵", "1×10⁻⁴", "1.5×10⁻⁴"])
ax2.xticks = 2:2:10
limits!(ax2, 1, 10, -2e-4, 2e-4)

yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2,])
ax1.yticklabelspace = yspace
ax2.yticklabelspace = yspace

rowgap!(ga, 10)

wavs = (15 .+ 1)./15 : 1 ./15: 10 .- 1 ./15

#labels = ["Center-Slope", "Up-SLope", "Center-Interior"]
i = 3  
lines!(ax1, wavs, Kₜ_Wavg[i,:], linewidth = 5, color = :dodgerblue2, )
lines!(ax2, wavs, Wₜ_Wavg[i,:], linewidth = 5, color = :firebrick2,)


#colsize!(f.layout, 2, Relative(0.05))

savename = "Paper_Cwtd_SlidUpGaussKW_wavg" 

save(apath * savename * ".png", f)



f = Figure(resolution = (900, 800), fontsize=26)
ga = f[1, 1] = GridLayout() # vhat

ax1 = Axis(ga[1, 1], ylabel = "W [ms⁻¹]", 
 title = "Tracer Weighted Diapycnal Diffusivity") 
#ax1.yticks = (0:4e-4:1.2e-3 , ["0", "4×10⁻⁴", "8×10⁻⁴", "1.2×10⁻³"])
hidexdecorations!(ax1)
#limits!(ax1, 1, 10, -3e-4, 1.3e-3)

ax2 = Axis(ga[2, 1], ylabel = "∂ₜb̄ [ms⁻³]", 
title =  "Tracer Weighted Buoyancy Velocity") 
hidexdecorations!(ax2)
#ax2.yticks = (-8e-5:4e-5:8e-5, ["-8×10⁻⁵", "-4×10⁻⁵", "0", "-4×10⁻⁵", "8×10⁻⁵"])
#limits!(ax2, 1, 10, -1.2e-4, 1.2e-4)

ax3 = Axis(ga[3, 1], ylabel = "|∇b̄| [s⁻²]", xlabel = "Tσ",
title =  "Tracer Weighted Buoyancy Gradient") 
#ax2.yticks = (-8e-5:4e-5:8e-5, ["-8×10⁻⁵", "-4×10⁻⁵", "0", "-4×10⁻⁵", "8×10⁻⁵"])
ax3.xticks = 2:2:10
#limits!(ax2, 1, 10, -1.2e-4, 1.2e-4)

rowgap!(ga, 10)

wavs1 = (600:600:160*600 )./ pm.Tσ
wavs2 = (0:600:160*600) ./ pm.Tσ

labels = ["Center-Slope", "Center-Interior", "Up-SLope"]
colors = [:dodgerblue2, :firebrick2, :black]
linestyle = [:solid, :solid, :solid]
for i = 1:2
    lines!(ax1, wavs1, Wₜ[i,:], linewidth = 5, color = colors[i], 
            linestyle = linestyle[i], label = labels[i])
    lines!(ax2, wavs1, ∂ₜb̄[i,:], linewidth = 5, color = colors[i],linestyle = linestyle[i] )
    lines!(ax3, wavs2, ∇b̄[i,:], linewidth = 5, color = colors[i],linestyle = linestyle[i] )
end

leg = Legend(ga[1, 2], ax1)

#colsize!(f.layout, 2, Relative(0.05))

#yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2, ax3])
#ax1.yticklabelspace = yspace
#ax2.yticklabelspace = yspace
#ax3.yticklabelspace = yspace

savename = "Paper_Cwtd_GaussW" 

save(apath * savename * ".png", f)



f = Figure(resolution = (500, 800), fontsize=26)
ga = f[1, 1] = GridLayout() # vhat

ax1 = Axis(ga[1, 1], ylabel = "W [ms⁻¹]", 
 title = "Tracer Weighted Diapycnal Velocity") 
#ax1.yticks = (0:4e-4:1.2e-3 , ["0", "4×10⁻⁴", "8×10⁻⁴", "1.2×10⁻³"])
hidexdecorations!(ax1)
#limits!(ax1, 1, 10, -3e-4, 1.3e-3)

ax2 = Axis(ga[2, 1], ylabel = "∂ₜb̄ [ms⁻³]", 
title =  "Tracer Weighted Buoyancy Velocity") 
hidexdecorations!(ax2)
#ax2.yticks = (-8e-5:4e-5:8e-5, ["-8×10⁻⁵", "-4×10⁻⁵", "0", "-4×10⁻⁵", "8×10⁻⁵"])
#limits!(ax2, 1, 10, -1.2e-4, 1.2e-4)

ax3 = Axis(ga[3, 1], ylabel = "|∇b̄| [s⁻²]", xlabel = "Tσ",
title =  "Tracer Weighted Buoyancy Gradient") 
#ax2.yticks = (-8e-5:4e-5:8e-5, ["-8×10⁻⁵", "-4×10⁻⁵", "0", "-4×10⁻⁵", "8×10⁻⁵"])
ax3.xticks = 2:2:10
#limits!(ax2, 1, 10, -1.2e-4, 1.2e-4)

rowgap!(ga, 10)

wavs1 = (600:600:160*600 )./ pm.Tσ
wavs2 = (0:600:160*600) ./ pm.Tσ

labels = ["Center-Slope", "Center-Interior", "Up-SLope"]
colors = [:dodgerblue2, :firebrick2, :black]
linestyle = [:solid, :solid, :solid]
for i = 3
    lines!(ax1, wavs1, Wₜ[i,:], linewidth = 5, color = colors[i], 
            linestyle = linestyle[i], label = labels[i])
    lines!(ax2, wavs1, ∂ₜb̄[i,:], linewidth = 5, color = colors[i],linestyle = linestyle[i] )
    lines!(ax3, wavs2, ∇b̄[i,:], linewidth = 5, color = colors[i],linestyle = linestyle[i] )
end

leg = Legend(ga[1, 2], ax1)

#colsize!(f.layout, 2, Relative(0.05))

#yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2, ax3])
#ax1.yticklabelspace = yspace
#ax2.yticklabelspace = yspace
#ax3.yticklabelspace = yspace

savename = "Paper_Cwtd_GaussW" 

save(apath * savename * ".png", f)
