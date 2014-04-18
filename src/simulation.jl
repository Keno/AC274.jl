abstract Simulation
abstract DG <: Simulation


# 1 Dimensional Discontinuous Galerkin Method

type DG1D <: DG
    p::Int64 # The order of ϕ
    K::Int64 # The number of cells to use
    a::Float64 # [a,b] is the interval on which to cumpute
    b::Float64 #
    fflux::Function # The numerical flux function
    Δt::Float64 # Time-step
    nt::Int64 # Number of time steps
    C::Float64 # Diffusion constant
    q₀
    mesh::Any
    boundary::Function
    hackinextraflux::Bool
    useMUSCL::Bool
    DG1D(p, K, a::Float64, b::Float64, fflux::Function, Δt::Float64,
        nt::Int64, C::Float64, q₀; mesh = nothing, boundary = t->zero(typeof(q₀(0))), hackinextraflux=false, useMUSCL = false) =
        new(p,K,a,b,fflux,Δt,nt,C,q₀,mesh,boundary,hackinextraflux, useMUSCL)
end

porder(p::DG1D) = p.p
nbf(p::DG1D) = p.p+1

zerop(p::DG1D) = zero(Float64)
numcells(p::DG1D) = p.K

# 2 Dimensional Discontinuous Galerkin Method

type DG2D <: DG
    mesh::Meshes.Mesh{Vertex2}
    porder::Int64
    fflux::Function # The numerical flux function
    Δt::Float64 # Time-step
    nt::Int64 # Number of time steps
    C::Float64 # Diffusion constant
    q₀
    boundary::Function
end

numcells(p::DG2D) = length(p.mesh.faces)
porder(p::DG2D) = p.porder

zerop(p::DG2D) = zero(Vertex2)

# sum_{k=0}^p (k + n - 1 \choose n - 1) for n = 2
nbf(p::DG2D) = div((porder(p)+1)*(porder(p)+2),2)
