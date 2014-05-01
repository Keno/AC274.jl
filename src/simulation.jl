abstract Simulation
abstract Galerkin <: Simulation
abstract DG <: Galerkin
abstract CG <: Galerkin

import Base: show

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
        nt::Int64, C::Float64, q₀; mesh = nothing, boundary = (args...)->zero(typeof(q₀(0))), hackinextraflux=false, useMUSCL = false) =
        new(p,K,a,b,fflux,Δt,nt,C,q₀,mesh,boundary,hackinextraflux, useMUSCL)
end

function show(io::IO, p::DG1D)
    println(io,"1 Dimensional Discontinuous Galerkin Simulation:")
    println(io,"    Ω = [",p.a,",",p.b,"]")
    println(io," mesh = 1D Uniform (K = ",p.K,")")
    println(io,"    p = ",p.p)
    println(io,"   Δt = ",p.Δt)
    println(io,"   nt = ",p.nt)
end

porder(p::DG1D) = p.p
𝖓(p::DG1D) = p.p+1

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

type CG2D{porder} <: DG
    mesh::Meshes.Mesh{Vertex2}
    dualmesh::Meshes.Mesh{Vertex2}
end


const Galerkin2D = Union(CG2D,DG2D)

function show(io::IO, p::Galerkin2D)
    println(io,"2D Galerkin Simulation:")
    println(io," mesh = 2D Triangular")
    println(io,"    p = ",porder(p))
end

numcells(p::DG2D) = length(p.mesh.faces)
porder(p::DG2D) = p.porder
porder{pp}(::CG2D{pp}) = pp

zerop(p::Galerkin2D) = zero(Vertex2)

# sum_{k=0}^p (k + n - 1 \choose n - 1) for n = 2
𝖓(p::Galerkin2D) = div((porder(p)+1)*(porder(p)+2),2)

# Utilities
const nbf = 𝖓
mesh(p::Galerkin) = p.mesh
