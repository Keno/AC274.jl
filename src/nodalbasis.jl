### Implementation of a nodal polynomial basis for elements in 1 and 2 dimensions
#
# Polynomials are written in coordinates on reference elements.
###


# 1D nodal basis

using Polynomials

# Chebyshev nodes
nodes1(p) = reverse([cos((2i-1)/(2p)*pi) for i = 1:p])
nodes(p::DG1D) = nodes1(nbf(p))
#nodes(p) = [-1,-1/3,1/3,1]

phi1(nodes::Vector) = (x = Poly([0.0,1.0]); Poly{Float64}[Poly([1.0])*prod([k == i ? 1 : (x-nodes[i])/(nodes[k]-nodes[i]) for i=1:length(nodes)]) for k=1:length(nodes)])
phi(p::DG1D) = phi1(nodes(p))
dphi(p::DG1D) = polyder(phi(p))

# 2D nodal basis

#
# Consider the folling sketch of a reference element. We label
# reference coordinates as ξ and η (not x and y as those refer to real
# coordinates in the mesh)
#
#
#               η
#               ▲
#               |
#             1 -      v3
#               |      |\
#               |      |  \
#               |      |    \
#               |      |      \
#               |      |        \
#               |      |          \
#             0 -      |____________\
#               |      v1           v2
#               |
#               |
#               +------|-------------|----▶ ξ
#                      0             1
#
#
#
# Note the choice of counterclockwise orientation when numbering
# The vertices.
#

## p = 0 ##
#
# One node and (1/3,1/3) with constant basis functions
#
###########

p0ϕ1(ξ,η) = 1.
∇p0ϕ1(ξ,η) = Vertex2(0.,0.)
const P0 = [p0ϕ1]
const DP0 = [∇p0ϕ1]
const N0 = [Vertex2(1//3,1//3)]

## p = 1 ##
#
# Three nodes at the vertices of the reference element with linear
# basis functions
#
###########

p1ϕ1(ξ,η) = 1 - ξ - η
p1ϕ2(ξ,η) = ξ
p1ϕ3(ξ,η) = η

∇p1ϕ1(ξ,η) = [-1,-1]
∇p1ϕ2(ξ,η) = [1,0]
∇p1ϕ3(ξ,η) = [0,1]

const P1 = [p1ϕ1,p1ϕ2,p1ϕ3]
const DP1 = [∇p1ϕ1,∇p1ϕ2,∇p1ϕ3]
const N1 = [Vertex2(0,0),Vertex2(1,0),Vertex2(0,1)]

## p = 2 ##
#
# Six nodes at the vertices + midpoints of the reference element with
# quadratic functions
#
###########

p2ϕ1(ξ,η) = 2.0*(0.5 - ξ - η)*(1.0 - ξ - η)
p2ϕ2(ξ,η) = 2.0ξ*(ξ - 0.5)
p2ϕ3(ξ,η) = 2.0η*(η - 0.5)
p2ϕ4(ξ,η) = 4.0ξ*(1. - ξ - η)
p2ϕ5(ξ,η) = 4.0η*ξ
p2ϕ6(ξ,η) = 4.0η*(1. - ξ - η)

# TODO: Make Vector2?
∇p2ϕ1(ξ,η) = Vertex2(4ξ+4η-3,4ξ+4η-3)
∇p2ϕ2(ξ,η) = Vertex2(4ξ - 1,0)
∇p2ϕ3(ξ,η) = Vertex2(0,4η - 1)
∇p2ϕ4(ξ,η) = Vertex2(-8ξ-4η+4,-4ξ)
∇p2ϕ5(ξ,η) = Vertex2(4η,4ξ)
∇p2ϕ6(ξ,η) = Vertex2(-4η,-8η-4ξ+4)

const P2 = [p2ϕ1,p2ϕ2,p2ϕ3,p2ϕ4,p2ϕ5,p2ϕ6]
const DP2 = [∇p2ϕ1,∇p2ϕ2,∇p2ϕ3,∇p2ϕ4,∇p2ϕ5,∇p2ϕ6]
const N2 = [Vertex2(0,0),Vertex2(1,0),Vertex2(0,1),Vertex2(0.5,0),Vertex2(0.5,0.5),Vertex2(0,0.5)]

### End of P=2 Basis functions

for f in [P0,P1,P2,DP0,DP1,DP2]
    @eval $(f.env.name)(v::Vertex2) = $(f.env.name)(v.coords.e1,v.coords.e2)
end

getPhi(p) = (@assert 0 <= p <= 2; p==0 ? P0 : p == 1 ? P1 : P2)
getDPhi(p) = (@assert 0 <= p <= 2; p==0 ? DP0 : p == 1 ? DP1 : DP2)
nodes2(p) = (@assert 0 <= p <= 2; p==0 ? N0 : p == 1 ? N1 : N2)
phi(p::DG2D) = getPhi(porder(p))
dphi(p::DG2D) = getDPhi(porder(p))
nodes(p::DG2D) = nodes2(porder(p))

# Evaluating coefficients on a basis
function evaluate_ref{N,T}(coeffs::Coeffs{N,T},basis)
    order(coeffs) == order(basis) || error("Basis and vector must agree")
    sum([a*f for (f,a) in zip(coeffs,basis)])
end

order(x::Array) = length(x)

function evaluate_ref(coeffs,basis,p)
    order(coeffs) == order(basis) || error("Basis and vector must agree (got $basis and $coeffs)")
    result = coeffs[1]*fapply(basis[1],p)
    for i = 2:length(basis)
        result += coeffs[i]*fapply(basis[i],p)
    end
    result
end

# Below are the hottest functions in the system. We unroll them manually 
# For maximum performce. This is usually a bad idea, but because they are 
# so important, we really have no other choice.

macro unroll_eval(call, coefff, basis)
    coefff = eval(coefff)
    funcs = eval(basis)
    adds = Expr(:block)
    for i = 2:length(funcs)
        addend = coefff(i,:($(Base.function_name(funcs[i]))(p)))
        push!(adds.args,:(@inbounds result += $addend))
    end
    Expr(:function,esc(call),quote
        @inbounds result = $(coefff(1,:($(Base.function_name(funcs[1]))(p))))
        $adds
        result
    end)
end

evalrefcoefff = (k,e)->:(coeffs[$k]*$e)
@unroll_eval evaluate_ref2dp0(coeffs,p) evalrefcoefff P0
@unroll_eval evaluate_ref2dp1(coeffs,p) evalrefcoefff P1
@unroll_eval evaluate_ref2dp2(coeffs,p) evalrefcoefff P2

cidevalrefcoefff = (k,e)->:(coeffs[1,cid,$k]*$e)
@unroll_eval evaluate_ref2dp0(coeffs,cid,p) cidevalrefcoefff P0
@unroll_eval evaluate_ref2dp1(coeffs,cid,p) cidevalrefcoefff P1
@unroll_eval evaluate_ref2dp2(coeffs,cid,p) cidevalrefcoefff P2

function evaluate_ref2d{T}(p,coeffs::Array{T},a,b) 
    po = porder(p)
    if po == 0
        return evaluate_ref2dp0(coeffs,a,b)
    elseif po == 1
        return evaluate_ref2dp1(coeffs,a,b)
    elseif po == 2
        return evaluate_ref2dp2(coeffs,a,b)
    else 
        error("Not specialized for this case")
    end
end
function evaluate_ref2d{T}(p,coeffs::Array{T},a) 
    po = porder(p)
    if po == 0
        return evaluate_ref2dp0(coeffs,a)
    elseif po == 1
        return evaluate_ref2dp1(coeffs,a)
    elseif po == 2
        return evaluate_ref2dp2(coeffs,a)
    else 
        error("Not specialized for this case")
    end
end
evaluate_ref2d(p,coeffs::Array{Coeffs},point) = evaluate_ref2d(p,coeffs.coeffs,point)


macro unroll_map(call, basis)
    funcs = eval(basis)
    cats = Expr(:vcat)
    for i = 1:length(funcs)
        push!(cats.args,:($(Base.function_name(funcs[i]))(p)))
    end
    Expr(:function,esc(call),quote
        $cats
    end)
end

@unroll_map evaluate_basis2dp0(p) P0
@unroll_map evaluate_basis2dp1(p) P1
@unroll_map evaluate_basis2dp2(p) P2

function evaluate_basis2d(p,point)
    po = porder(p)
    if po == 0
        return evaluate_basis2dp0(point)
    elseif po == 1
        return evaluate_basis2dp1(point)
    elseif po == 2
        return evaluate_basis2dp2(point)
    else 
        error("Not specialized for this case")
    end
end

evalbasiscoefff = (k,e)->:(vRHS[$k] += mul*$e)
@unroll_eval evaluate_basis2dp0!(vRHS,mul,p) evalbasiscoefff P0
@unroll_eval evaluate_basis2dp1!(vRHS,mul,p) evalbasiscoefff P1
@unroll_eval evaluate_basis2dp2!(vRHS,mul,p) evalbasiscoefff P2

function evaluate_basis2d!(vRHS,p,mul,point)
    po = porder(p)
    if po == 0
        return evaluate_basis2dp0!(vRHS,mul,point)
    elseif po == 1
        return evaluate_basis2dp1!(vRHS,mul,point)
    elseif po == 2
        return evaluate_basis2dp2!(vRHS,mul,point)
    else 
        error("Not specialized for this case")
    end
end

const evaluate = evaluate_ref
