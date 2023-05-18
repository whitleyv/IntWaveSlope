using CUDA: has_cuda_gpu

using Oceananigans
using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
using Oceananigans.Grids: Center, Face
using Oceananigans.Operators

function get_outputs_tuple(model; Fulld=false, )

    # Output: primitive fields + computations
    u, v, w, b, c = merge(model.velocities, model.tracers)
    νₑ = model.diffusivity_fields.νₑ
    κₑ = model.diffusivity_fields.κₑ.b

    _, yb, zb = nodes(model.tracers.b)

    outputs = merge(model.velocities, model.tracers)

    ∂b∂x = Field(∂x(b))
    ∂b∂y = Field(∂y(b))
    ∂b∂z = Field(∂z(b))

    ∇b² = Field(@at (Center, Center, Center) ∂b∂x^2 + ∂b∂y^2 + ∂b∂z^2)

    (xlength, ylength, zlength) = size(b)

    Ygrid = reshape(repeat(yb, zlength), ylength, zlength)
    YXgrid = reshape(repeat(yb, zlength*xlength), ylength, zlength, xlength)
    SlopeGridYX = curvedslope.(YXgrid)
    SlopeGridY = curvedslope.(Ygrid)
    boolZYX = (SlopeGridYX .+ 2) .<= zb' # all the values greater than slope
    non_boolZY = (SlopeGridY .+ 2) .> zb' # all the values in the slope

    first_allowed_y_index = sum(non_boolZY, dims=1) .+ 1
    first_z_index_in_domain = sum(first_allowed_y_index .> ylength) + 1

    boolZ_perm = permutedims(boolZYX, [2, 1,3])

    permdimsb = permutedims(interior(b), [3,2,1])

    sortb = sort(Array(permdimsb[boolZ_perm]))

    background_APEb = zeros(xlength, ylength, zlength)

    running_index_total = 1
    # to reshape you need to take into account the slope:
    for k = first_z_index_in_domain:zlength
        first_yidx = first_allowed_y_index[k]
        yL = length(first_yidx:ylength)
        num_indices_in_plane = yL*xlength
        new_total = running_index_total + num_indices_in_plane
        background_APEb[:,first_yidx:ylength,k] = reshape(sortb[running_index_total:new_total-1], xlength,yL)
        running_index_total = new_total
    end
    

    if Fulld
        ## ertel PV
        dWdy = Field(∂y(w))
        dVdz = Field(∂z(v))

        dUdz = Field(∂z(u))
        dWdx = Field(∂x(w))

        dVdx = Field(@at (Face, Center, Center) ∂x(v)) # failed (FFC)
        dUdy = Field(@at (Center, Face, Center) ∂y(u)) # failed (FFC)

        ## pseudo dissipation rate
        dUdx = Field(∂x(u)) # ccc
        dVdy = Field(∂y(v)) # ccc
        dWdz = Field(∂z(w)) # ccc

        ddx² = Field(dUdx^2 + dVdx^2 + dWdx^2)  
        ddy² = Field(dUdy^2 + dVdy^2 + dWdy^2)
        ddz² = Field(dUdz^2 + dVdz^2 + dWdz^2) 

        dsum = Field(ddx² + ddy² + ddz²)
        ϵ = Field(νₑ * dsum)

        outputs = merge(outputs, (; ϵ=ϵ, νₑ=νₑ, κₑ=κₑ, N2=∂b∂z, Eb=background_APEb))

    else
        outputs = merge(outputs, (; κₑ=κₑ, N2=∂b∂z, Eb=background_APEb))
    end

    return outputs

end


