module AC274

export neighbor, has_neighbor, DG1D, DG2D, CG2D, solve, generateMesh, Mesh1D,
    plotSolution, invchi, chi, 𝜒⁻¹, 𝜒, plotMesh, evaluate, evaluate_ref, interpolate,
    chi, invchi, vertices, mesh, nbf, 𝖓, cid, do_quad_ref, ∂, rewriteMesh, calciᴰ, ℳ,
    𝒩, createDualMesh, apply∂D!, applyA∂D!, drawSupport, generateMatrices, stiffness,
    loadV

include("simpletypes.jl")
include("mesh.jl")
include("simulation.jl")
include("quadrature.jl")
include("neighbors.jl")
include("nodalbasis.jl")
include("RHS.jl")
include("solve.jl")
include("plot1d.jl")
include("plot2d.jl")
include("cg.jl")


end