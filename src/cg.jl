function localstiffness(p::CG2D{1},c,m::Mesh,cache)
    D = [(c.p3[1] - c.p2[1]) (c.p1[1] - c.p3[1]) (c.p2[1] - c.p1[1]);
          (c.p3[2] - c.p2[2]) (c.p1[2] - c.p3[2]) (c.p2[2] - c.p1[2])]
    D'D/(4*cache.elemJ[cid(c)])
end

function localstiffness(p::CG2D,c,m::Mesh,cache)
    pdphi = AC274.dphi(p)
    Al = Array(Float64,𝖓(p),𝖓(p))
    invA = inv(AC274.Ak(c)')
    for i = 1:𝖓(p), j = 1:𝖓(p)
        Al[i,j] = 
        cache.elemJ[cid(c)]*AC274.do_quad_ref(x->dot(invA*pdphi[i](x),invA*pdphi[j](x)),p)
    end
    Al
end

𝒩(p::CG2D{1}) = length(mesh(p).vertices)
𝒩(p::CG2D{2}) = length(mesh(p).vertices) + length(p.dualmesh.vertices)

# Assemble stiffness matrix
function stiffness{pp}(p::CG2D{pp},is∂D,cache=AC274.generateMatrices(p))
    m = mesh(p)

    A = spzeros(𝒩(p),𝒩(p))
    for c in mesh(p)
        Al = localstiffness(p,c,m,cache)
        for i = 1:𝖓(p), j = 1:𝖓(p)
            A[ℳ(p,c,pp,i),ℳ(p,c,pp,j)] += Al[i,j]
        end
    end
    applyA∂D!(A,calciᴰ(p,is∂D,cache.neighbors))
    A
end

function applyA∂D!(A,iᴰ)
    for i in iᴰ
        A[i,1:end] = spzeros(1,size(A,2))
        A[i,i] = 1.0
    end
    A
end
applyA∂D!(A,is∂D,cache,m::Mesh) = applyA∂D!(A,calciᴰ(m,is∂D,cache.neighbors))


function loadV(p::CG2D, is∂D, f, gN, gD, cache=AC274.generateMatrices(p))
    iᴰ = calciᴰ(p,is∂D,cache.neighbors)
    F = loadVbase(p, f, cache)
    apply∂N!(p, F, is∂D, gN, cache)
    apply∂D!(mesh(p), F, gD, iᴰ)
    F
end

function loadVbase(p::CG2D, f, cache)
    basis, dphi, ns = cache.basis, cache.dphi, cache.nodes
    m = mesh(p)
    F = zeros(length(m.vertices))
    
    for c in mesh(p)
        for i = 1:nbf(p)
            F[ℳ(p,c,i)] += cache.elemJ[cid(c)]*do_quad_ref(x->(f(x)*basis[i](x)),p)
        end
    end

    F
end

function apply∂N!(p, F, is∂D, gN, cache)
    m = mesh(p)
    for c in m
        for face in c
            if !has_neighbor(m, c, face, cache) && !is∂D(face.p1)
                vRHS = zeros(Float64,nbf(p))
                gx, gw = cache.Δ1quadpoints
                # Manually inline quadrature
                factor = sqrt(dot(face.p2-face.p1,face.p2-face.p1))/2
                for i in 1:length(gw)
                    point⁻, w = gx[i], gw[i]
                    cpoint⁻ = ∂(c,face,point⁻)
                    AC274.evaluate_basis2d!(vRHS,p,-gN(face.p1)*factor*w,cpoint⁻)
                end
                for i = 1:nbf(p)
                    F[ℳ(p,c,i)] += vRHS[i]
                end
            end
        end    
    end
end

function apply∂D!(m, F, gD, iᴰ)
    idxs = Int64[i for i in iᴰ]
    F[idxs] = [gD(x) for x in m.vertices[idxs]]
end

# Calculate indicies for dirichlet nodes indicated by a given filter function is∂D
function calciᴰ{pp}(p::CG2D{pp},is∂D,neighbors = computeNeighbors(mesh(p)))
    iᴰ = Set{Int64}()
    for c in mesh(p)
        for edge in c
            if !has_neighbor(mesh(p), c, edge, neighbors) && is∂D(edge.p1) && is∂D(edge.p2)
                c = mesh(p)[cid(edge)]
                id = AC274.edgeid(c,edge)
                # Add both to the dirichlet boundary
                push!(iᴰ,ℳ(p,c,pp,id))
                push!(iᴰ,ℳ(p,c,pp,mod1(id+1,3)))
                pp == 2 && push!(iᴰ,ℳ(p,c,pp,3+id))
            end
        end
    end
    iᴰ  
end