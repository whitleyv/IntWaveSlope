using Statistics
using Printf
using JLD2
using CairoMakie
using GeometryBasics

dpath = "ReviewerAnalysis/ReviewerData/"
apath = "ReviewerAnalysis/ReviewerPlots/"

sn = "U350N100Lz100g100"

files_fluxes = dpath * "BFlux_Rotated_" * sn * ".jld2"

fluxfile = jldopen(files_fluxes, "r+")

vb_phasedepWavg = fluxfile["vb_phasedepWavg"]
wb_phasedepWavg = fluxfile["wb_phasedepWavg"]
#ŵ_Wavg = fluxfile["ŵ_Wavg"]
#v̂_Wavg = fluxfile["v̂_Wavg"]
phase_times = fluxfile["phase_times"]
boolZYX = fluxfile["boolZYX"]
yb = fluxfile["yb"]
zb = fluxfile["zb"]

ylength = 880 

include("../parameters.jl")
pm = getproperty(SimParams(), Symbol(sn))

pm = merge(pm, (; Tanθ = sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                Tanα = pm.γ * sqrt((pm.σ^2 - pm.f^2)/(pm.Ñ^2-pm.σ^2)),
                nz = round(Int,pm.Lz/2),
                m = π/pm.Lz,
                l = sqrt(((π/pm.Lz)^2 * (pm.f^2 - pm.σ^2)) / (pm.σ^2 - pm.Ñ^2)),
                Tf = 2*π/pm.f, 
                Tσ = 2*π/pm.σ))


α = atan(pm.Tanα)

# rotated derivatives
function y_deriv(A, h)
    step = 1/h
    step2 = 1/(h*2)
    da = zeros(size(A));
    da[:,1, :] = (A[:,2,:] .- A[:,1,:]) .* step
    da[:,end, :] = (A[:,end, :] .- A[:,end-1,:]) .* step
    da[:,2:end-1,:] = (A[:,3:end,:] .- A[:,1:end-2,:]) .* step2

    return da
end

function z_deriv(A, h)
    step = 1/h
    step2 = 1/(h*2)
    da = zeros(size(A));
    da[:,:,1] = (A[:,:,2] .- A[:,:,1]) .* step
    da[:,:,end] = (A[:,:,end] .- A[:,:, end-1]) .* step
    da[:,:,2:end-1] = (A[:,:,3:end] .- A[:,:, 1:end-2]) .* step2

    return da
end

Δy = (yb[2]-yb[1])
Δz = (zb[2]-zb[1])

∂ŷ_vb_phasedepWavg = y_deriv(vb_phasedepWavg, Δy) .* cos(α) .- z_deriv(vb_phasedepWavg, Δz) .* sin(α)
∂ẑ_wb_phasedepWavg = y_deriv(wb_phasedepWavg, Δy) .* sin(α) .+ z_deriv(wb_phasedepWavg, Δz) .* cos(α)

∂ŷ_vb_phasedepWavg_xavg = mean(∂ŷ_vb_phasedepWavg, dims=1)[1,:,:]
∂ẑ_wb_phasedepWavg_xavg = mean(∂ẑ_wb_phasedepWavg, dims=1)[1,:,:]

∇̂_ub_phasedepWavg_xavg = ∂ŷ_vb_phasedepWavg_xavg .+ ∂ẑ_wb_phasedepWavg_xavg

# curved topo parameters
const gausT_center = 895                                 # gaussian paramereters for curved "corner" of slope
const gausT_width = 180
const ySlopeSameˢ = 1332.22                           # point where planar and curved corner math up the best

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

land = curvedslope.(yb)
land_pdel = (curvedslope.(yb) .+ pm.U₀/pm.Ñ)[1:382]
δ = pm.U₀/pm.Ñ
zlength = length(zb)
ylength = 880

nabmaxp = 3e-7 

###[ v'b'p  ][ w'b'p  ][ -∇p        ]

phlims = (-nabmaxp, nabmaxp)
    
yb_cut = yb[1:ylength]
        
f1 = Figure(resolution = (1800, 500), fontsize=26)
        ga = f1[1, 1] = GridLayout()    
        gb = f1[1, 2] = GridLayout()
        gc = f1[1, 3] = GridLayout()
    
        gcb1 = f1[1, 4] = GridLayout()       
    
        axvb = Axis(ga[1, 1], ylabel = "z [m]", xlabel = "y [m]",) #aspect = DataAspect()) #vh
        axwb = Axis(gb[1, 1], xlabel = "y [m]") 
        axgrad = Axis(gc[1, 1], xlabel = "y [m]") #dz w'b'

        axvb.xticks = 500:500:1000
        axwb.xticks = 500:500:1000
        axgrad.xticks = 500:500:1000

        axvb.yticks = [-250, 0]
        
        limits!(axvb, 0, 1500, -500, 0)
        limits!(axwb, 0, 1500, -500, 0)
        limits!(axgrad, 0, 1500, -500, 0)
    
        hideydecorations!(axwb)
        hideydecorations!(axgrad)
    
        colsize!(f1.layout, 4, Relative(0.05))
    
        hv = heatmap!(axvb, yb_cut[2:end], zb[2:end], -1 .* ∂ŷ_vb_phasedepWavg_xavg, colormap = :balance, colorrange = phlims)
            band!(axvb, yb, land,-500, color=:black)
            lines!(axvb, yb[1:382], land_pdel, color=:black, linewidth = 4, linestyle = :dash)
            text!(axvb, Point.(50, -450), text = "- ∂ŷ(v̂̃b̃)", align = (:left, :center), color = :white, 
            font = :bold, fontsize = 26)
    
        hw = heatmap!(axwb, yb_cut, zb[2:end], -1 .* ∂ẑ_wb_phasedepWavg_xavg, colormap = :balance, colorrange = phlims)
            band!(axwb, yb, land,-500, color=:black)
            lines!(axwb, yb[1:382], land_pdel, color=:black, linewidth = 4, linestyle = :dash)
            text!(axwb, Point.(50, -450), text = "- ∂ẑ(ŵ̃b̃)", align = (:left, :center), color = :white, 
            font = :bold, fontsize = 26)

        hnab = heatmap!(axgrad, yb_cut[2:end], zb[2:end], -1 .* ∇̂_ub_phasedepWavg_xavg, colormap = :balance, colorrange = phlims)
            band!(axgrad, yb, land,-500, color=:black)
            lines!(axgrad, yb[1:382], land_pdel, color=:black, linewidth = 4, linestyle = :dash)
            text!(axgrad, Point.(50, -450), text = "- ∇̂⋅(u⃗̂̃b̃)", align = (:left, :center), color = :white, 
            font = :bold, fontsize = 26)

        ycut = 250/pm.Tanα:1:250/pm.Tanα + 10
        ycut2 = (250/pm.Tanα + 30):1:(250/pm.Tanα + 40)

        arrows!(axvb, [250/pm.Tanα], [-250], [1.5*δ], [-1.5*δ*pm.Tanα], arrowsize = 17, linewidth = 7,  
            arrowcolor = :gray30, linecolor = :gray30) # normal
        arrows!(axvb, [250/pm.Tanα], [-250], [1.5*δ*pm.Tanα], [1.5*δ], arrowsize = 17, linewidth = 7,  
            arrowcolor = :gray30, linecolor = :gray30) # tangential
        #lines!(axvb,  yb, yb./pm.Tanα .- 250/pm.Tanα^2 .- 250)
        #lines!(axvb,  yb, yb./pm.Tanα .- 250/pm.Tanα^2 .- 250 .- 30/pm.Tanα .+ 30*pm.Tanα)
        #lines!(axvb,  yb, -1 .* yb*pm.Tanα )
        #lines!(axvb,  yb, -1 .* yb*pm.Tanα .+ 60)

        text!(axvb, [250/pm.Tanα + 70], [-250+ 110], text="ẑ", align = (:left, :center), 
            color = :gray30, font = :bold, fontsize = 35)
        text!(axvb, [250/pm.Tanα + 10], [-250 - 50], text="ŷ", align = (:left, :center), 
            color = :gray30, font = :bold, fontsize = 35)

        ySL = 30
        zSL = 30
                  #corner left,  
        #squarey = [250/pm.Tanα, 250/pm.Tanα + ySL, 250/pm.Tanα + ySL + zSL, 250/pm.Tanα + zSL]
        #squarez = [-250, -250-ySL*pm.Tanα, -250+ySL/pm.Tanα + (zSL*pm.Tanα), -250 + zSL*pm.Tanα]
        #scatter!(axvb,Point2f.(squarey,squarez))
        #poly!(axvb, Polygon(Point2f.(squarey,squarez)), color=:white, strokecolor = :gray30, linewidth = 5)
    
        # create colorbars the size of the whole data set
        cb1 = Colorbar(gcb1[1,1], hv, ticks = (-nabmaxp/2:nabmaxp/2:nabmaxp/2), size =35, label = "Wave Terms [ms⁻³]")

        colgap!(ga, 15)
        colgap!(gb, 15)
        colgap!(gc, 15)
    
save(apath * "BFluxTrip_Rotated_" * sn * ".png", f1, px_per_unit = 2)


