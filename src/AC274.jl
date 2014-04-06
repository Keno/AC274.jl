# AC 274 Assignment 1

module AC274

# Note, this is a quick-and-dirty implementation.
# It's lacking abstraction and performance and it is simply meant
# to be quick to write and easy to understand.

using Polynomials
using Gadfly

import Base: getindex, start, next, done, length, size, eltype, promote_rule


export neighbor, has_neighbor, DG1D, solve, generateMesh, Mesh1D, 
    plotSolution, invchi, chi, 𝜒⁻¹, 𝜒

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

𝜒(m::Mesh1D,k,x) = 2*(x-m.elements[k])/step(m.elements) - 1
𝜒⁻¹(m::Mesh1D,k,x) = m.elements[k] + (1/2)*(1+x)*(step(m.elements))
𝜒(m::Mesh1D,c::Cell1D,x) = 𝜒(m,c.coord,x)
𝜒⁻¹(m::Mesh1D,c::Cell1D,x) = 𝜒⁻¹(m,c.coord,x)


const invchi = 𝜒⁻¹
const chi = 𝜒

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
next(x::Mesh1D,i) = (x[i],i+1)
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
nodes(p) = reverse([cos((2i-1)/(2p)*pi) for i = 1:p])
#nodes(p) = [-1,-1/3,1/3,1]

phi(nodes::Vector) = (x = Poly([0.0,1.0]); [Poly([1.0])*prod([k == i ? 1 : (x-nodes[i])/(nodes[k]-nodes[i]) for i=1:length(nodes)]) for k=1:length(nodes)])
phi(p::Int64) = phi(nodes(p+1))

integrate(poly::Poly, a, b) = (pp = polyint(poly); fapply(pp,b) - fapply(pp,a))
integrate(x::Number, a, b) = x*(b-a)

const gp = Base.gauss(Float64, 7)
function integrate(f::Function,a, b)
    @assert a == -1 && b == 1
    result = 0.0
    for (x,w) in zip(gp...)
        result += w*f(x)
    end
    result
end

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
    boundary::Function
    hackinextraflux::Bool
    useMUSCL::Bool
    DG1D(p, K, a::Float64, b::Float64, fflux::Function, Δt::Float64,
        nt::Int64, C::Float64, q₀; mesh = nothing, boundary = t->zero(typeof(q₀(0))), hackinextraflux=false, useMUSCL = false) =
        new(p,K,a,b,fflux,Δt,nt,C,q₀,mesh,boundary,hackinextraflux, useMUSCL)
end

immutable Coeffs{N,T}
    coeffs::Vector{T} #NTuple{N,Float64}
end
Coeffs{T}(N,x::Vector{T}) = (Coeffs{N,T}(x))

order(x::Coeffs) = length(x.coeffs) #N
order(x::Vector) = length(x)
eltype{N,T}(c::Coeffs{N,T}) = T
eltype{N,T}(::Type{Coeffs{N,T}}) = T

fapply(f::Poly,x) = polyval(f,x)
fapply(vec::Vector,x) = [fapply(f,x) for f in vec]
fapply(f,x) = f(x)
fapply(x::Number,y) = x

promote_rule{N,T}(::Type{Coeffs{N,T}},::Type{Coeffs{N,T}}) = Coeffs{N,T}
promote_rule{N,T,S<:Number}(::Type{Coeffs{N,T}},::Type{S}) = Coeffs{N,T}

+{N,T}(x::Coeffs{N,T}, y::Coeffs{N,T}) = Coeffs{N,T}(x.coeffs+y.coeffs)
+{N,T}(x::Coeffs{N,T},y::Vector{Float64}) = Coeffs{N,T}(x.coeffs+y)
-{N,T}(x::Coeffs{N,T},y::Vector{Float64}) = Coeffs{N,T}(x.coeffs-y)
*{N,T}(s::Number, x::Coeffs{N,T}) = Coeffs{N,T}(s*x.coeffs)
.*{N,T}(s, x::Coeffs{N,T}) = Coeffs{N,T}(s.*x.coeffs)
.^{N,T}(s::Coeffs{N,T}, y) = Coeffs{N,T}(s.coeffs.^y)
size(x::Coeffs) = size(x.coeffs)

for f in (:getindex, :start, :done, :next, :length)
    @eval ($f)(x::Coeffs,args...) = ($f)(x.coeffs,args...)
end

import Base: zero, one, zeros, ones, conj

zero{N,T}(::Type{Coeffs{N,T}}) = Coeffs{N,T}(zeros(Float64,N))
one{N,T}(::Type{Coeffs{N,T}}) = Coeffs{N,T}(ones(Float64,N))
zeros{N,T}(::Type{Coeffs{N,T}},n) = [zero(Coeffs{N,T}) for _=1:n]
ones{N,T}(::Type{Coeffs{N,T}},n) = [one(Coeffs{N,T}) for _=1:n]
conj(x::Coeffs) = x

function evaluate{N,T}(coeffs::Coeffs{N,T},basis)
    order(coeffs) == order(basis) || error("Basis and vector must agree")
    sum([a*f for (f,a) in zip(coeffs,basis)])
end

function evaluate{N,T}(coeffs::Coeffs{N,T},basis,p)
    order(coeffs) == order(basis) || error("Basis and vector must agree")
    result = 0.0
    for i = 1:length(basis)
        result += coeffs[i]*fapply(basis[i],p)
    end
    result
end


# ElemJacobian
elemJ(x,k) = (1//2)*(x.mesh[k,:r] - x.mesh[k,:l])

# MassMatrix

⊗(arg::Union(Poly,Number),args::Union(Poly,Number)...) = *(arg,args...)
⊗(f::Function,fs::Function...) = x->*(f(x),[f(x) for f in fs]...)
⊖(f,g) = x->fapply(f,x)-fapply(g,x)
∘(f,g) = x->fapply(f,fapply(g,x))

Mlocal(p::DG1D,k,basis) = Float64[integrate(⊗(basis[i],basis[j],elemJ(p,k)),-1,1) for i=1:(p.p+1), j=1:(p.p+1) ]
M(p::DG1D) = (basis = phi(p.p); Diagonal([Mlocal(p,k,basis) for k in p.mesh]))

# RHS
coord(face) = face == :l ? -1 : 1
oppcoord(face) = -coord(face)

n⁻(face) = face == :l ? -1 : 1

immutable Cache
    phidphi::Array{Float64,2}
    basis::Array{Poly{Float64},1}
    dphi::Array{Poly{Float64},1}
    ns::Array{Float64,1}
end

function generateMatrices(p)
    ns = nodes(p.p+1)
    basis = phi(ns)
    dphi = polyder(basis)

    phidphi = [integrate(⊗(basis[j],dphi[i]),-1,1) for j = 1:(p.p+1), i = 1:(p.p+1)]
    Cache(phidphi,basis,dphi,ns)
end

function RHS(p,k,Q,t; cache=generateMatrices(p))
    basis, dphi, ns = cache.basis, cache.dphi, cache.ns

    Qk = Q[k]
    Fk = typeof(Qk)([p.fflux(e,[𝜒⁻¹(p.mesh,k,x) for x in ns],t)::typeof(e) for e in Qk])
    vRHS = zeros(eltype(Fk),p.p+1)
    for i = 1:p.p+1, j = 1:p.p+1
        vRHS[i] += Fk[j]*cache.phidphi[j,i]
    end

    cell = p.mesh[k]
    for face in cell

        x = cell[face]

        q⁻ = evaluate(Qk,basis,coord(face))

        if has_neighbor(p.mesh, cell, face)

            q⁺ =  evaluate(Q[neighbor(p.mesh, cell, face).coord],basis,oppcoord(face))

            # f﹡ϕ_i(x_r)n-
            f1 = n⁻(face) * (
                (1/2)*(p.fflux(q⁻,x,t)::typeof(q⁻) + p.fflux(q⁺,x,t)::typeof(q⁺)) + # average term
                (1/2)*p.C*( q⁻ - q⁺ )*n⁻(face) # Jump term
            )
            for i=1:p.p+1
                vRHS[i] -= f1*fapply(basis[i],coord(face))
            end

            if p.hackinextraflux
            # Hack in extra flux
            vRHS[:] -= -integrate(x->(y = 𝜒(p.Mesh1D,k,x);fapply(evaluate(Qk,basis),x)*fapply(basis,x)*
                (2*t*(y^2-1)-2*y*(y^2+1)^2)/(2t*y + (y^2+1)^2)^2),-1,1)
            end
        else
            vRHS[:] -= n⁻(face) * fapply(basis,coord(face)) * (1//2)*p.fflux(q⁻,x,t)::typeof(q⁻)

            # boundary conditions
            vRHS[:] -= n⁻(face) * fapply(basis,coord(face)) * (1/2)*p.fflux(p.boundary(q⁻,t,face),x,t)::typeof(q⁻)
        end
    end
    vRHS
end

function globalSolve(p::DG1D,MMatrix,Q,t; cache=generateMatrices(p))
    k = Array(Coeffs{p.p+1,eltype(eltype(Q))},p.K)
    for j = 1:p.K
        r = Coeffs(p.p+1,MMatrix[j]\RHS(p,j,Q,t; cache=cache));
        k[j] = r
    end
    k
end

function _solve(p, ℚ; cache = generateMatrices(p))
    basis = phi(p.p)
    # ℚ[1] is the inital condition
    M = [factorize(Mlocal(p,k,basis)) for k = 1:p.K]
    for i in 2:p.nt
        k1 = globalSolve(p,M,ℚ[i-1,:],i*p.Δt; cache=cache)
        k2 = globalSolve(p,M,ℚ[i-1,:] + (0.5*p.Δt*k1)',i*p.Δt; cache=cache)
        k3 = globalSolve(p,M,ℚ[i-1,:] + (0.5*p.Δt*k2)',i*p.Δt; cache=cache)
        k4 = globalSolve(p,M,ℚ[i-1,:] + (p.Δt*k3)',i*p.Δt; cache=cache)
        ℚ[i,:] = ℚ[i-1,:] + p.Δt/6 * (k1+2*k2+2*k3+k4)';
    
        # post-processing
        if p.useMUSCL
            @assert p.p == 1 # For now
            for c in p.mesh
                Qk = ℚ[i,c.coord]
                ql = evaluate(Qk,basis,coord(:l))
                qr = evaluate(Qk,basis,coord(:r))

                h = (c.right-c.left)

                slope = (qr-ql)/h
                qkbar = (ql+qr)/2

                s = sign(slope)
                minabsa = abs(slope)

                for face in (:l,:r)
                    if has_neighbor(p.mesh,c, face)
                        N = neighbor(p.mesh, c, face)
                        Qk = ℚ[i,N.coord]
                        ql = evaluate(Qk,basis,coord(:l))
                        qr = evaluate(Qk,basis,coord(:r))
                        qNbar = (ql+qr)/2

                        eeslope = n⁻(face)*(qNbar-qkbar)/h

                        minabsa = min(minabsa,abs(eeslope))

                        s += sign(eeslope)
                    end
                end

                s/=3
                # Hack in knowledge about WaveState - will fix after
                # assignment is due

                local q1c, q2c
                if abs(s).q1 != 1
                    q1c = [qkbar.q1 for i=1:(p.p+1)]
                else
                    slope = s*minabsa
                    q1c = [(qkbar+(node)*(h/2)*slope).q1 for node in cache.ns]
                end

                if abs(s).q2 != 1
                    q2c = [qkbar.q2 for i=1:(p.p+1)]
                else
                    slope = s*minabsa
                    q2c = [(qkbar+(node)*(h/2)*slope).q2 for node in cache.ns]
                end

                @show (s, q1c, q2c)

                cs = typeof(s)[typeof(s)(q1c[i],q2c[i]) for i = 1:length(q1c)]
                ℚ[i,c.coord] = Coeffs(p.p+1,cs)
            end
        end
    end
    ℚ
end

# Nodal interpolation
interpolate(f, nodes) = (Coeffs(length(nodes),map(f,nodes)))

function solve(p::DG1D)
    if p.mesh === nothing
        p.mesh = generateMesh(Mesh1D,p.a,p.b,p.K; periodic=false)
    end

    cache = generateMatrices(p)


    # Save coefficients for all timesteps (nt)
    # for all K cells
    

    if isa(p.q₀,Vector)
        T = eltype(p.q₀)
        ℚ = Array(Coeffs{p.p+1,T},p.nt,p.K)
        ℚ[1,:] = p.q₀ #[zero(Coeffs{p.p+1}) for i = 1:p.K]
    elseif isa(p.q₀,Function)
        T = typeof(p.q₀(0))
        ℚ = Array(Coeffs{p.p+1,T},p.nt,p.K)
        # Use nodal interpolation
        for i = 1:p.K
            ℚ[1,i] = interpolate(p.q₀,map(n->𝜒⁻¹(p.mesh,i,n),nodes(p.p+1)))
        end
    elseif isa(p.q₀, Number)
        # Uniform initial conditions
        T = typeof(p.q₀)
        ℚ = Array(Coeffs{p.p+1,T},p.nt,p.K)
 
        ℚ[1,:] = p.q₀*ones(Coeffs{p.p+1,T},p.K)
    else
        @assert "Unknown option"
    end

    _solve(p,ℚ; cache=cache)
end

function pf(p,mesh,a,b,Q)
    basis = phi(p+1)
    function (x)
        k = min(1+div(x-a,(b-a)/K),K)
        poly = evaluate(Q[k],basis,𝜒(mesh,k,x))
    end
end

function pf(p,Q)
    basis = phi(p.p)
    function (x)
        k = min(1+div(x-p.a,(p.b-p.a)/p.K),p.K)
        poly = evaluate(Q[k],basis,𝜒(p.mesh,k,x))
    end
end


plotSolution(p::DG1D,Q) = plot(pf(p,Q),p.a,p.b)

end

#include("2d.jl")

## End of library code
