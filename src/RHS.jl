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
    nodes::Array{Vertex2,1}
    neighbors::Array{Int64,2}
    Δ1quadpoints::(Array{Float64,1},Array{Float64,1})
    realspacenodes::Array{Array{Vertex2,1},1}
end

function generateMatrices(p::DG2D) 
    pphi, pdphi = phi(p),dphi(p)
    phidphi = [
        Vertex2[ -elemJ(p,c)*do_quad_ref(x->(pphi[i](x)*inv(Ak(c)')*pdphi[j](x)),p) for i = 1:nbf(p), j = 1:nbf(p)]
        for c in p.mesh]
    ns = nodes(p)
    realspacenodes = [[𝜒⁻¹(p.mesh,c,point) for point in ns] for c in p.mesh]
    Cache2D(phidphi, pphi, pdphi, nodes(p), computeNeighbors(p.mesh),
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
    nodes::Array{S,1}
    realspacenodes::Array{Array{Float64,1},1}
end
Cache{S,T}(phidphi::Array{S,2},phi::Array{T,1},dphi::Array{T,1},ns::Array{S,1}) =
    Cache{S,T}(phidphi,phi,dphi,ns)

function nbd(p::DG1D)
    ns = nodes(p)
    basis = phi(p)
    db = polyder(basis)
    (ns,basis,db)
end

nbd(p) = (nodes(p),phi(p),dphi(p))

function generateMatrices(p::DG1D)
    ns, phi, dphi = nbd(p)
    phidphi = [ do_quad_ref(⊗(phi[i],dphi[j]),p) for i = 1:nbf(p), j = 1:nbf(p)]
    realspacenodes = [Float64[𝜒⁻¹(p.mesh,c,point) for point in ns] for c in p.mesh]
    Cache(phidphi,phi,dphi,ns,realspacenodes)
end

fflux(p,q,cell,face,point,t) = p.fflux(q,cell,face,point,t)
#fflux(q,cell,face,t) = Main.fflux(q,cell,face,t)



#::typeof(zerop(p))
dot(a::Union(Vector2,Vertex2),b::Union(Vector2,Vertex2)) = a[1]*b[1]+a[2]*b[2]
function laxcf(p,q⁻,q⁺,cell,face,point,t)
    a = n⁻(face)
    b = fflux(p,q⁻,cell,face,point,t) + 
        fflux(p,q⁺,cell,face,point,t) # average term
    return (1/2)*(dot(a,b) + p.C*( q⁻ - q⁺ )) # Jump term
end

function fluxpostp!(vRHS,p,cell,face,Qk,t,cache)
    basis = cache.basis
    if isa(p,DG1D) && p.hackinextraflux
        # Hack in extra flux
        vRHS[:] -= -integrate(x->(y = 𝜒(p.mesh,cell,x); evaluate_ref(Qk,basis,x)*fapply(basis,x)*
            (2*t*(y^2-1)-2*y*(y^2+1)^2)/(2t*y + (y^2+1)^2)^2),-1,1)
    end
end

function computeFlux!(vRHS::Array,p::DG1D,cell::Cell1D,face,Q,t,cache)
    basis, dphi, ns = cache.basis, cache.dphi, cache.nodes

    x = coord(face)

    oppedge = oppface(p.mesh, cell, face)
    Qk = Q[1,cid(cell),:]
    Qn = Q[1,cid(neighbor(p.mesh,cell,face)),:]

    q⁻ = evaluate_ref(Qk,basis,x)
    q⁺ = evaluate_ref(Qn,basis,coord(oppedge))

    # f﹡ϕ_i(x_r)n-
    f1 = laxcf(p,q⁻,q⁺,cell,face,cell[face],t)

    for i=1:nbf(p)
        vRHS[i] -= f1*fapply(basis[i],x)
    end

    fluxpostp!(vRHS,p,cell,face,Qk,t,cache)
end

function computeBoundary!(vRHS::Array,p::DG1D,cell::Cell1D,face,Q,t,cache)
    basis, dphi, ns = cache.basis, cache.dphi, cache.nodes

    x = coord(face)
    Qk = Q[1,cid(cell),:]
    q⁻ = evaluate(Qk,basis,x)
    vRHS[:] -= n⁻(face) * fapply(basis,x) * (1//2)*fflux(p,q⁻,cell,face,x,t)::eltype(vRHS)
            
    # boundary conditions
    vRHS[:] -= n⁻(face) * fapply(basis,coord(face)) * (1/2)*fflux(p,p.boundary(q⁻,t,face),cell,face,x,t)::eltype(vRHS)
end

oppcoord(edge,oppedge,point⁻) = 𝜒(oppedge,𝜒⁻¹(edge,point⁻))

function computeFlux!(vRHS::Array,p::DG2D,cell::Cell2D,face::Edge,Q,t,cache)
    basis, dphi, ns = cache.basis, cache.dphi, cache.nodes

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
    gx, gw = cache.Δ1quadpoints
    # Manually inline quadrature
    factor = sqrt(dot(face.p2-face.p1,face.p2-face.p1))/2
    for i in 1:length(gw)
        point⁻, w = gx[i], gw[i]

        point = 𝜒⁻¹(face,point⁻)

        cpoint⁻ = ∂(cell,face,point⁻)
        q⁻ = evaluate_ref2d(p,Q,cid(cell),cpoint⁻)

        f = boundaryflux(p,q⁻,cell,face,point,t)

        #vRHS[:] -= factor*w*f*fapply(basis,cpoint⁻)
        b = evaluate_basis2d!(vRHS,p,-f*factor*w,cpoint⁻)
    end
end

has_neighbor(p::Mesh1D, cell, face, cache) = has_neighbor(p,cell,face)
has_neighbor(p::Meshes.Mesh, cell::Cell2D, face::Edge, cache::Cache2D) = has_neighbor(p,cell,face,cache.neighbors)


phidphi(cache::Cache2D,cell,i,j) = cache.phidphi[cid(cell)][i,j]
phidphi(cache::Cache,cell,i,j) = cache.phidphi[i,j]

function RHS{T}(p,cell::Cell,Q::Array{T,3},t,cache=generateMatrices(p))
    basis, dphi, ns = cache.basis, cache.dphi, cache.nodes

    vRHS = zeros(eltype(Q),nbf(p))

    for i = 1:nbf(p)
        Fki = fflux(p,Q[1,cid(cell),i],cell,0,cache.realspacenodes[cid(cell)][i],t)
        for j = 1:nbf(p)
            vRHS[j] += dot(phidphi(cache,cell,i,j),Fki)
        end
    end

    for face in cell
        if has_neighbor(p.mesh, cell, face, cache)
            computeFlux!(vRHS,p,cell,face,Q,t,cache)
        else
            computeBoundary!(vRHS,p,cell,face,Q,t,cache)
        end
    end

    vRHS
end