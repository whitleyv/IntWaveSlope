using CUDA: has_cuda_gpu

using Oceananigans
using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
using Oceananigans.Grids: Center, Face
using Oceananigans.Operators
import Oceananigans.TurbulenceClosures: viscosity, diffusivity
using Oceananigans.Fields: @compute 
using Oceanostics


function get_outputs_tuple(model; Fulld=false)

    # Output: primitive fields + computations
    b = model.tracers.b
    νₑ = model.diffusivity_fields.νₑ
    
    κ = diffusivity(model.closure, model.diffusivity_fields, Val(:b))
    # eddy diffusivity field
    κₑ = κ.a
    # scalar
    κₗ = κ.b
 
    outputs = merge(model.velocities, model.tracers)

    if Fulld
        ε = Field(KineticEnergyDissipationRate(model))

        dBdz = Field(@at (Center, Center, Center) ∂z(b))
        dBdy = Field(@at (Center, Center, Center) ∂x(b))
        dBdx = Field(@at (Center, Center, Center) ∂y(b))

        ∇b² = Field(@at (Center, Center, Center) dBdx^2 + dBdy^2 + dBdz^2)

        outputs = merge(outputs, (; ϵ=ε, N2=dBdz, νₑ=νₑ, κₑ=κₑ, ∇b²=∇b²,))
    end

    return outputs

end