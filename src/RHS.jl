# MassMatrix

⊗(arg::Union(Poly,Number),args::Union(Poly,Number)...) = *(arg,args...)
⊗(f::Function,fs::Function...) = x->*(f(x),[f(x) for f in fs]...)
⊖(f,g) = x->fapply(f,x)-fapply(g,x)
∘(f,g) = x->fapply(f,fapply(g,x))

Mlocal(p::DG,c,basis) = elemJ(p,c)*Float64[do_quad_ref(⊗(basis[i],basis[j]),p) for i=1:nbf(p), j=1:nbf(p) ]
M(p::DG) = (basis = phi(p); Diagonal([Mlocal(p,c,basis) for c in p.mesh]))

# RHS 2D 

type Cache2D
    phidphi::Array{Array{Vertex2,2},1}
    basis::Array{Function,1}
    dphi::Array{Function,1}
    ns::Array{Vertex2,1}
    neighbors::Array{Int64,2}
    elemJ::Array{Float64,1}
    Δ1quadpoints::(Array{Float64,1},Array{Float64,1})
    realspacenodes::Array{Array{Vertex2,1},1}
end

function generateMatrices(p::DG2D) 
    pphi, pdphi = phi(p),dphi(p)
    phidphi = [
        Vertex2[ do_quad_ref(x->(pphi[i](x)*inv(Ak(c)')*pdphi[j](x)),p) for i = 1:nbf(p), j = 1:nbf(p)]
        for c in p.mesh]
    ns = nodes(p)
    realspacenodes = [[𝜒⁻¹(p.mesh,c,point) for point in ns] for c in p.mesh]
    Cache2D(phidphi, pphi, pdphi, nodes(p), computeNeighbors(p.mesh), [elemJ(p,c) for c in p.mesh],
        Base.gauss(Float64,max(1,2*porder(p))),realspacenodes)
end

export n⃗, midp

function n⃗(f::Edge) 
    t⃗ = f.p1-f.p2
    t⃗ ./= norm(t⃗)
    Vector2{Float64}(t⃗.coords.e2,-t⃗.coords.e1)
end

n⁻(f::Edge) = n⃗(f::Edge)

midp(f::Edge) = (1/2)*(f.p1+f.p2)

# RHS

immutable Cache{S,T}
    phidphi::Array{S,2}
    basis::Array{T,1}
    dphi::Array{T,1}
    ns::Array{S,1}
end
Cache{S,T}(phidphi::Array{S,2},phi::Array{T,1},dphi::Array{T,1},ns::Array{S,1}) =
    Cache{S,T}(phidphi,phi,dphi,ns)

function nbd(p::DG1D)
    ns = nodes(p)
    basis = phi(p)
    db = dphi(basis)
    (ns,basis,db)
end

nbd(p) = (nodes(p),phi(p),dphi(p))

function generateMatrices(p::DG1D)
    ns, phi, dphi = nbd(p)
    phidphi = [ do_quad_ref(⊗(phi[i],dphi[j]),p) for j = 1:nbf(p), i = 1:nbf(p)]
    Cache(phidphi,phi,dphi,ns)
end

fflux(p,q,cell,face,point,t) = p.fflux(q,cell,face,point,t)::typeof(zerop(p))
#fflux(q,cell,face,t) = Main.fflux(q,cell,face,t)



#::typeof(zerop(p))
function laxcf(p,q⁻,q⁺,cell,face,point,t)
    a = n⁻(face)
    b = fflux(p,q⁻,cell,face,point,t) + 
        fflux(p,q⁺,cell,face,point,t) # average term
    return (1/2)*(a[1]*b[1]+a[2]*b[2] + p.C*( q⁻ - q⁺ )) # Jump term
end

function computeFlux!(p::DG1D,cell,face,Q,t; cache = nothing)
    oppedge = oppface(p.mesh, cell, face)
    Qn = Q[cid(oppedge)]

    q⁻ = evaluate(Qk,basis,face)
    q⁺ = evaluate(Qn,basis,oppedge)

    # f﹡ϕ_i(x_r)n-
    f1 = laxcf(p,q⁻,q⁺,cell,face,t)

    for i=1:nbf(p)
        vRHS[i] -= f1*fapply(basis[i],coord(face))
    end

    if isa(p,DG1D) p.hackinextraflux
    # Hack in extra flux
    vRHS[:] -= -integrate(x->(y = 𝜒(p.Mesh1D,cell,x);fapply(evaluate(Qk,basis),x)*fapply(basis,x)*
        (2*t*(y^2-1)-2*y*(y^2+1)^2)/(2t*y + (y^2+1)^2)^2),-1,1)
    end
end

oppcoord(edge,oppedge,point⁻) = 𝜒(oppedge,𝜒⁻¹(edge,point⁻))

function computeFlux!(vRHS::Array,p::DG2D,cell::Cell2D,face::Edge,Q,t,cache)
    basis, dphi, ns = cache.basis, cache.dphi, cache.ns

    oppedge = oppface(p.mesh, cell, face, cache.neighbors)

    gx, gw = cache.Δ1quadpoints
    # Manually inline quadrature
    factor = sqrt(dot(face.p2-face.p1,face.p2-face.p1))/2
    for i in 1:length(gw)
        point⁻, w = gx[i], gw[i]
        point⁺ = oppcoord(face,oppedge,point⁻)

        point = 𝜒⁻¹(face,point⁻)

        cpoint⁻ = ∂(cell,face,point⁻)
        q⁻ = evaluate_ref2d(p,Q,cid(cell),cpoint⁻)
        q⁺ = evaluate_ref2d(p,Q,cid(oppedge),∂(p.mesh[cid(oppedge)],oppedge,point⁺))

        # f﹡ϕ_i(x_r)n-
        f = laxcf(p,q⁻,q⁺,cell,face,point,t)
        #vRHS[:] -= factor*w*f*fapply(basis,cpoint⁻)
        b = evaluate_basis2d!(vRHS,p,-f*factor*w,cpoint⁻)
    end
end

boundaryflux(p,q⁻,cell,face,point,t) = p.boundaryflux(q⁻,cell,face,point,t)

function computeBoundary!(vRHS::Array,p::DG2D,cell::Cell2D,face::Edge,Q,t,cache)
    Qk = Q[1,cid(cell),:]

    gx, gw = cache.Δ1quadpoints
    # Manually inline quadrature
    factor = sqrt(dot(face.p2-face.p1,face.p2-face.p1))/2
    for i in 1:length(gw)
        point⁻, w = gx[i], gw[i]

        point = 𝜒⁻¹(face,point⁻)

        cpoint⁻ = ∂(cell,face,point⁻)
        q⁻ = evaluate_ref2d(p,Qk,cpoint⁻)

        # f﹡ϕ_i(x_r)n-
        f = boundaryflux(p,q⁻,cell,face,point,t)

        #vRHS[:] -= factor*w*f*fapply(basis,cpoint⁻)
        b = evaluate_basis2d!(vRHS,p,-f*factor*w,cpoint⁻)
    end
end


function RHS{T}(p,cell,Q::Array{T,3},t,cache=generateMatrices(p))
    basis, dphi, ns = cache.basis, cache.dphi, cache.ns

    vRHS = zeros(eltype(Q),nbf(p))

    for face in cell
        if has_neighbor(p.mesh, cell, face, cache.neighbors)
            computeFlux!(vRHS,p,cell,face,Q,t,cache)
        else
            computeBoundary!(vRHS,p,cell,face,Q,t,cache)

            #vRHS[:] = 0
            #vRHS[:] -= n⁻(face) * fapply(basis,coord(face)) * (1//2)*p.fflux(q⁻,cell,face,t)::typeof(q⁻)

            # boundary conditions
            #vRHS[:] -= n⁻(face) * fapply(basis,coord(face)) * (1/2)*p.fflux(p.boundary(q⁻,t,face),cell,face,t)::typeof(q⁻)
        end
    end

    for i = 1:nbf(p)
        Fki = fflux(p,Q[1,cid(cell),i],cell,0,cache.realspacenodes[cid(cell)][i],t)
        #@show Fkj
        for j = 1:nbf(p)
            vRHS[j] += cache.elemJ[cid(cell)]*dot(Fki,cache.phidphi[cid(cell)][i,j])
        end
    end

    vRHS
end