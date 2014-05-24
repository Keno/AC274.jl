# This file is ugly. Don't look at it. Should be replaced by fixed-size arrays/tuples

using Polynomials

import Base: getindex, start, done, next, length, size, eltype, 
             promote_rule, zero, one, zeros, ones, conj, copy

using ImmutableArrays
using Meshes
import Meshes: Vertex2

# Coefficient vector for a polynomial basis
immutable Coeffs{N,T}
    coeffs::Vector{T} #NTuple{N,Float64}
end
Coeffs{T}(N,x::Vector{T}) = (Coeffs{N,T}(x))

copy{N,T}(c::Coeffs{N,T}) = Coeffs{N,T}(copy(c.coeffs))

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

zero{N,T}(::Type{Coeffs{N,T}}) = Coeffs{N,T}(zeros(Float64,N))
one{N,T}(::Type{Coeffs{N,T}}) = Coeffs{N,T}(ones(Float64,N))
zeros{N,T}(::Type{Coeffs{N,T}},n) = [zero(Coeffs{N,T}) for _=1:n]
ones{N,T}(::Type{Coeffs{N,T}},n) = [one(Coeffs{N,T}) for _=1:n]
conj(x::Coeffs) = x

# Vertex2 (from Meshes)

zero(::Type{Vertex2}) = Vertex2(zero(Vector2{Float64}))

import Base: *, -, +, dot, ./, \

for f in (:(-),:(+))
    @eval ($f)(a::Vertex2,b::Vertex2) = Vertex2(($f)(a.coords,b.coords))
end
dot(a::Vertex2,b::Vertex2) = dot(a.coords,b.coords)
*(a::Float64,b::Vertex2) = Vertex2(a*b.coords.e1,a*b.coords.e2)
+(a::Vector2{Float64},b::Vertex2) = Vertex2(a+b.coords)
-(a::Vector2{Float64},b::Vertex2) = Vertex2(a-b.coords)
+(b::Vertex2,a::Vector2{Float64}) = a+b
./(a::Vertex2,b::Number) = Vertex2(a.coords./b)
\(a,b::Vertex2) = (r = a\[b.coords.e1,b.coords.e2]; Vertex2(r[1],r[2]))
\{T}(a::Array{T,2},b::Vector2{T}) = (r = a\[b.e1,b.e2]; Vector2{T}(r[1],r[2]))
*(a::Array{Float64,2}, b::Vertex2) = (r=a*[b.coords.e1,b.coords.e2]; Vertex2(r[1],r[2]))

# Manually inlined, because this function is hot in laxcf
dot(a::Vector2,b::Vertex2) = a.e1*b.coords.e1 + a.e2*b.coords.e2

function (*){T,S}(A::StridedMatrix{T}, x::Vector2{S})
    Vector2{S}((A*[x[1];x[2]])...)
end
