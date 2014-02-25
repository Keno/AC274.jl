# AC 274 Assignment 1

# Note, this is a quick-and-dirty implementation.
# It's lacking abstraction and performance and it is simply meant 
# to be quick to write and easy to understand. 

using Polynomial
using Gadfly

import Base: getindex, start, next, done, length


# Mesh generation

abstract Mesh
abstract Cell

immutable Cell1D
    left::Float64
    right::Float64
    coord::Int64
end

start(x::Cell1D) = start((:l,:r))
next(x::Cell1D,i) = next((:l,:r),i)
done(x::Cell1D,i) = done((:l,:r),i)

immutable Mesh1D <: Mesh
    elements::Base.FloatRange{Float64}
    isperiodic::Bool
end

# This terminology is slightly misleading as we probably should consider the
# mesh a finite cochain complex rather than the coordinate neighborhoods of a
# smooth manifold (in particular they are not overlapping), but I don't have 
# defer that to a later time.
atlas(m::Mesh1D) = (k,x)->((m.elements[k]+(1/2)*step(m.elements))-x)/step(m.elements)

invatlas(m::Mesh1D) = (k,x)->( m.elements[k] + (1/2)*(1+x)*(step(m.elements)) )

function has_neighbor(m::Mesh1D, cell, face)
    m.isperiodic && return true
    if cell.coord == 1 && face == :l
        return false
    elseif cell.coord == length(m.elements)-1 && face == :r
        return false
    else 
        return true
    end 
end

function neighbor(m::Mesh1D, cell, face)
    @assert has_neighbor(m, cell, face)
    if face == :l
        if m.isperiodic && cell.coord == 1
            return m[length(m.elements)-1]
        else
            return m[cell.coord - 1]
        end
    else
        @assert face == :r
        if m.isperiodic && cell.coord == length(m.elements)-1
            return m[1]
        else
            return m[cell.coord + 1]
        end
    end
end

length(x::Mesh1D) = length(x.elements) - 1
start(x::Mesh1D) = 1
next(x::Mesh1D,i) = (i,i+1)
done(x::Mesh1D,i) = i > length(x)

generateMesh(::Type{Mesh1D},a,b,K; periodic=false) =  Mesh1D(a:((b-a)/K):b, periodic)

getindex(m::Mesh1D, cell::Int64) = Cell1D(m.elements[cell], m.elements[cell+1], cell)
getindex(m::Mesh1D, cell::Int64, face) = m[cell][face]
function getindex(c::Cell1D, face) 
    @assert face == :r || face == :l
    (face == :r) ? c.right : c.left
end

# Interpolation basis

# Chebyshev nodes
nodes(p) = [cos((2i-1)/(2p)*pi) for i = 1:p]
#nodes(p) = [-1,-1/3,1/3,1]

phi(nodes::Vector) = (x = Poly([0.0,1.0]); [Poly([1.0])*prod([k == i ? 1 : (x-nodes[i])/(nodes[k]-nodes[i]) for i=1:length(nodes)]) for k=1:length(nodes)])
phi(p::Int64) = phi(nodes(p+1))


integrate(poly::Poly, a, b) = -(map(x->fapply(Polynomial.integrate(poly),x),(b,a))...)
integrate(x::Number, a, b) = x*(b-a)

# In general, Gauss-Kronrod quadrature of order 7
integrate(x::Function, a, b) = quadgk(x,a,b;order=7)[1]

# Driver code

abstract Simulation


type DG1D <: Simulation
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
    inflow::Function
    DG1D(p, K, a::Float64, b::Float64, fflux::Function, Δt::Float64,
        nt::Int64, C::Float64, q₀; mesh = nothing, inflow = t->0) = 
        new(p,K,a,b,fflux,Δt,nt,C,q₀,mesh,inflow)
end

immutable Coeffs{N} 
    coeffs::Vector{Float64} #NTuple{N,Float64}
end

order(x::Coeffs) = length(x.coeffs) #N
order(x::Vector) = length(x)

fapply(f::Poly,x) = papply(f,x)
fapply(vec::Vector,x) = map(f->fapply(f,x),vec)
fapply(f,x) = f(x)
fapply(x::Number,y) = x

+{N}(x::Coeffs{N}, y::Coeffs{N}) = Coeffs{N}(x.coeffs+y.coeffs)
+{N}(x::Coeffs{N},y::Vector{Float64}) = Coeffs{N}(x.coeffs+y)
-{N}(x::Coeffs{N},y::Vector{Float64}) = Coeffs{N}(x.coeffs-y)
*{N}(s::Number, x::Coeffs{N}) = Coeffs{N}(s*x.coeffs)
.*{N}(s::Number, x::Coeffs{N}) = Coeffs{N}(s.*x.coeffs)
.^{N}(s::Coeffs{N}, y) = Coeffs{N}(s.coeffs.^y)

for f in (:getindex, :start, :done, :next, :length)
    @eval ($f)(x::Coeffs,args...) = ($f)(x.coeffs,args...)
end

import Base: zero, one, zeros, ones, conj

zero{N}(::Type{Coeffs{N}}) = Coeffs{N}(zeros(Float64,N))
one{N}(::Type{Coeffs{N}}) = Coeffs{N}(ones(Float64,N))
zeros{N}(::Type{Coeffs{N}},n) = [zero(Coeffs{N}) for _=1:n]
ones{N}(::Type{Coeffs{N}},n) = [one(Coeffs{N}) for _=1:n]
conj(x::Coeffs) = x

function evaluate{N}(coeffs::Coeffs{N},basis)
    order(basis) == order(basis) || error("Basis and vector must agree")
    sum([a*f for (f,a) in zip(coeffs,basis)])
end

# ElemJacobian
elemJ(x,k) = (1//2)*(x.mesh[k,:r] - x.mesh[k,:l])

# MassMatrix

⊗(arg::Union(Poly,Number),args::Union(Poly,Number)...) = *(arg,args...)
⊗(f::Function,fs::Function...) = x->*(f(x),[f(x) for f in fs]...)

Mlocal(p::DG1D,k,basis) = Float64[integrate(⊗(basis[i],basis[j],elemJ(p,k)),-1,1) for i=1:(p.p+1), j=1:(p.p+1) ]
M(p::DG1D) = (basis = phi(p.p); Diagonal([Mlocal(p,k,basis) for k in p.mesh]))

# RHS
coord(face) = face == :l ? -1 : 1
oppcoord(face) = -coord(face)

n⁻(face) = face == :l ? -1 : 1

function RHS(p,k,Q,t)
    basis = phi(p.p)
    dphi = derivative(basis)

    Fk = p.fflux(Q[k])
    vRHS = Float64[ sum([ Fk[j]*integrate(⊗(basis[j],dphi[i]),-1,1) for j = 1:(p.p+1) ]) for i = 1:(p.p+1)]

    cell = p.mesh[k]
    for face in cell

        q⁻ = fapply(evaluate(Q[k],basis),coord(face))

        if has_neighbor(p.mesh, cell, face)

            
            q⁺ = fapply(evaluate(Q[neighbor(p.mesh, cell, face).coord],basis),oppcoord(face))
            normal = n⁻(face)

            # f﹡ϕ_i(x_r)n-
            vRHS[:] -= 
                n⁻(face) * fapply(basis,coord(face)) * (
                    (1//2)*(p.fflux(q⁻) + p.fflux(q⁺)) + # average term
                    (1//2)*p.C*( q⁻ - q⁺ )*n⁻(face) # Jump term
                )
        else 
            vRHS[:] -= n⁻(face) * fapply(basis,coord(face)) * (1//2)*p.fflux(q⁻)
            
            # boundary conditions
            # Outflow BC
            if face == :r
                vRHS[:] -= n⁻(face) * fapply(basis,coord(face)) * (1//2)*p.fflux(q⁻)
            elseif face == :l
                vRHS[:] -= n⁻(face) * fapply(basis,coord(face)) * p.fflux(p.inflow(t))
            end
        end
    end
    vRHS
end

function globalSolve(p::DG1D,MMatrix,Q,t)
    k = Array(Coeffs{p.p+1},p.K)
    for j = 1:p.K
        k[j] = Coeffs{p.p+1}(MMatrix[j]\RHS(p,j,Q,t));
    end
    k
end

function _solve(p, ℚ)
    basis = phi(p.p)
    # ℚ[1] is the inital condition
    M = [Mlocal(p,k,basis) for k = 1:p.K]
    for i in 2:p.nt
        k1 = globalSolve(p,M,ℚ[i-1,:],i*p.Δt)
        k2 = globalSolve(p,M,ℚ[i-1,:] + (0.5*p.Δt*k1)',i*p.Δt)
        k3 = globalSolve(p,M,ℚ[i-1,:] + (0.5*p.Δt*k2)',i*p.Δt)
        k4 = globalSolve(p,M,ℚ[i-1,:] + (p.Δt*k3)',i*p.Δt)
        ℚ[i,:] = ℚ[i-1,:] + p.Δt/6 * (k1+2*k2+2*k3+k4)';
    end

    ℚ
end

# Nodal interpolation
interpolate(f, nodes) = Coeffs{length(nodes)}(reverse(map(f,nodes)))

function solve(p::DG1D)
    if p.mesh === nothing
        p.mesh = generateMesh(Mesh1D,p.a,p.b,p.K; periodic=false)
    end

    # Save coefficients for all timesteps (nt) 
    # for all K cells
    ℚ = Array(Coeffs{p.p+1},p.nt,p.K)

    if isa(p.q₀,Vector)
        ℚ[1,:] = p.q₀ #[zero(Coeffs{p.p+1}) for i = 1:p.K]
    elseif isa(p.q₀,Function)
        # Use nodal interpolation
        𝜒⁻¹ = invatlas(p.mesh)
        for i = 1:p.K
            ℚ[1,i] = interpolate(p.q₀,map(n->𝜒⁻¹(i,n),nodes(p.p+1)))
        end
    elseif isa(p.q₀, Number)
        # Uniform initial conditions
        ℚ[1,:] = p.q₀*ones(Coeffs{p.p+1},p.K)
    else 
        @assert "Unknown option"
    end

    _solve(p,ℚ)
end

function plotSolution(p::DG1D,Q)
    𝜒 = atlas(p.mesh)
    basis = phi(p.p)
    plot(function (x)
        k = 1+div(x-p.a,(p.b-p.a)/p.K)
        fapply(evaluate(Q[k],basis),(𝜒(k,x)))
    end,p.a,p.b)
end

## End of library code
