using Statistics
using Printf
using CUDA: has_cuda_gpu
using ArgParse

using Oceananigans
using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
using Oceananigans.Grids: Center, Face
using Oceananigans.Operators


arch = has_cuda_gpu() ? GPU() : CPU()

underlying_grid = RectilinearGrid(arch; size = (3, 5, 10), 
                                            x = (0, 1), y = (0, 1),
                                            z = (-1, 0),
                                            halo = (4, 4, 4),
                                            topology = (Periodic, Bounded, Bounded))

@inline curvedslope(y) = -2.0 * y

@inline is_immersed(x, y, z) = z < curvedslope(y)

immersed_grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBoundary(is_immersed))

start_time = time_ns() # time to start the clock for how long it's been running

f = 1.0

kwargs_model = (closure = SmagorinskyLilly(),
              advection = WENO(),
                tracers = (:b, :c),
                coriolis = FPlane(f=f),
                buoyancy = BuoyancyTracer(),
            timestepper = :RungeKutta3)

model = NonhydrostaticModel(grid = immersed_grid, ; kwargs_model...)

b₀(x, y, z) = 0.003^2 * z  

set!(model, u=0, v=0, w=0, b=b₀, c=1)
path_name = "/Users/vicwhit/Desktop/Research/Code/IntWaveSlope/"
@show data_path = path_name * "PEtestrun"

Δt = 0.5
simulation = Simulation(model, Δt = Δt, stop_time = 2.0)

progress_message(sim) =
        @printf("i: %04d, t: %s, Δt: %s, wall time: %s\n",
        sim.model.clock.iteration, prettytime(sim.model.clock.time),
        prettytime(sim.Δt), prettytime((time_ns() - start_time) * 1e-9))

simulation.callbacks[:progress] = Callback(progress_message, TimeInterval(0.5))


# Output: primitive fields + computations
u, v, w, b, c = merge(model.velocities, model.tracers)
νₑ = model.diffusivity_fields.νₑ
κₑ = model.diffusivity_fields.κₑ.b

outputs = merge(model.velocities, model.tracers)

∂b∂x = Field(∂x(b))
∂b∂y = Field(∂y(b))
∂b∂z = Field(∂z(b))

∇b² = Field(@at (Center, Center, Center) ∂b∂x^2 + ∂b∂y^2 + ∂b∂z^2)


function flattenedsort(A, dim_order::Union{Tuple, AbstractVector}, boolZ::BitArray)
    return sort(Array(permutedims(A, dim_order)[boolZ]))
end

function reshape_around_slope(sortb::AbstractVector, grid, k1, j1)
    background_APEb = zeros(xlength, ylength, zlength)

    running_index_total = 1

    # to reshape you need to take into account the slope:
    for k = k1:grid.Nz
        first_yidx = j1[k]
        yL = length(first_yidx:grid.Ny)
        num_indices_in_plane = yL*grid.Nx
        new_total = running_index_total + num_indices_in_plane
        background_APEb[:,first_yidx:grid.Ny,k] = reshape(sortb[running_index_total:new_total-1], grid.Nx,yL)
        running_index_total = new_total
    end

    return background_APEb
end

function sort_b(model; )

    b = model.tracers.b
    _, yb, zb = nodes(model.tracers.b)
    grid = model.grid

    # creating a grid of y values
    Ygrid = reshape(repeat(yb, grid.Nz), grid.Ny, grid.Nz)
    YXgrid = reshape(repeat(yb, grid.Nz*grid.Nx), grid.Ny, grid.Nz, grid.Nx)

    # creating a grid of the slope value associated with y
    SlopeGridYX = curvedslope.(YXgrid)
    SlopeGridY = curvedslope.(Ygrid)

    # creating a boolean grid of values above the topography
    boolZYX = (SlopeGridYX ) .<= zb' # all the values greater than slope
    non_boolZY = (SlopeGridY) .> zb' # all the values in the slope

    # first y index to the right of the slope for each z level
    first_allowed_y_index = sum(non_boolZY, dims=1) .+ 1

    #first z index where the allowable y index is in the first_z_index_in_domain
    # since tanh function bottom vals don't come into play ...
    first_z_index_in_domain = sum(first_allowed_y_index .> ylength) + 1

    # permuting boolean to correct order
    boolZ_perm = permutedims(boolZYX, [2, 1, 3])

    # creating a 1D array of sorted buoyancy values above the slope
    sorted_B = flattenedsort(interior(b), [3,2,1], boolZ_perm)

    # reshaping 1D array into region above slope 
    return reshape_around_slope(sorted_B, grid, first_z_index_in_domain, first_allowed_y_index)

end

E_b = (mod)->sort_b(mod; )

outputs = merge(outputs, (; κₑ=κₑ, N2=∂b∂z, Eb=E_b, ∇b²=∇b²))

simulation.output_writers[:fields] = JLD2OutputWriter(model, outputs;
                                                          schedule = TimeInterval(1.0),
                                                          filename = data_path,
                                                          overwrite_existing = true)


run!(simulation)

