# mesh: Basic data structures for meshes

abstract Cell

import Meshes: AbstractMesh, Mesh

import Base: start, next, done, length, 
             getindex, norm, reverse

### Notation used in this file
#
# Note: If you are just trying to read the code, there is no need to read this section
#       it is included for preciseness when annotating code below with comments.
#
#
# The mesh for a domain Ω ⊂ R^k is denoted 𝕄(Ω) and Δⁿ(Ω) denote the set
# n(=dim(Ω)) simplicies on Ω.
#
# A mesh is the datum a set of vertices V(𝕄(Ω)) ⊂ Ω and cells (aka faces) C(𝕄(Ω)) ⊂ Δⁿ(Ω) togerther 
# with an ismorphism cid: C(𝕄(Ω)) → (𝕁 ⊂ ℕ) that assigns a unique index to every face of the 
# mesh. In the current implmentation cid⁻¹ = getindex(𝕄,⋅).
#
# We let E: Δⁿ(Ω) → ⊕ⁿ⁺¹ Δⁿ⁻¹(Ω) denote the set of faces of the cells (edges of the mesh) 
###

#
# Note currently we only support embeddings where k = n i.e. the dimension of the 
# domain is equal to the dimension of the ambient space (no 2d simulation on the surface of a sphere, etc.)
#

### Functionality provided
#
# Iteration:
#   - Over meshes 𝕄 gives iteration over F(𝕄)
#   - Over cells C gives iteration over E(C)
# 
# Indexing (getindex): 
#   - Into meshes 𝕄 gives the cell index (i.e. cid⁻¹ = getindex(𝕄,⋅))
#   - On cells C, depends on n, but the numbering will always depend on the numbering of 
#   - the corresponding reference element 𝜒(C) rather than the orientation of C itself
####

# 1D Meshes
#
# Meshes where dim(Ω) = 1 (i.e. at the moment also dim(k) = 1)
#
# C(𝕄(Ω)) ⊂ Δ¹ = Line segments
# E(C) ⊂ Δ⁰ = Points
#

# A cell in a 1D Mesh 
immutable Cell1D <: Cell
    left::Float64
    right::Float64
    coord::Int64
end

cid(c::Cell1D) = c.coord 
start(x::Cell1D) = start((:l,:r))
next(x::Cell1D,i) = next((:l,:r),i)
done(x::Cell1D,i) = done((:l,:r),i)

# The actual 1D mesh
immutable Mesh1D <: AbstractMesh
    elements::Base.FloatRange{Float64}
    isperiodic::Bool
end

generateMesh(::Type{Mesh1D},a,b,K; periodic=false) =  Mesh1D(a:((b-a)/K):b, periodic)

# Number of cells in the mesh
length(x::Mesh1D) = length(x.elements) - 1

# Iterate over cells in the mesh
start(x::Mesh1D) = 1
next(x::Mesh1D,i) = (x[i],i+1)
done(x::Mesh1D,i) = i > length(x)

𝜒(m::Mesh1D,k,x) = 2*(x-m.elements[k])/step(m.elements) - 1
𝜒⁻¹(m::Mesh1D,k,x) = m.elements[k] + (1/2)*(1+x)*(step(m.elements))
𝜒(m::Mesh1D,c::Cell1D,x) = 𝜒(m,c.coord,x)
𝜒⁻¹(m::Mesh1D,c::Cell1D,x) = 𝜒⁻¹(m,c.coord,x)

getindex(m::Mesh1D, cell::Int64) = Cell1D(m.elements[cell], m.elements[cell+1], cell)
getindex(m::Mesh1D, cell::Int64, face) = m[cell][face]
function getindex(c::Cell1D, face)
    @assert face == :r || face == :l
    (face == :r) ? c.right : c.left
end

coord(face) = face == :l ? -1 : 1
oppcoord(face) = -coord(face)

oppface(cell::Cell1D,face) = face == :l ? :r : :l

n⁻(face) = face == :l ? -1 : 1


# 2D Meshes
#
# Meshes where dim(Ω) = 2 (i.e. at the moment also dim(k) = 2)
#
# C(𝕄(Ω)) ⊂ Δ² = Triangles
# E(C) ⊂ Δ¹ = Line segments
#
# Note that in two dimensions the orientation becomes important. The orientation 
# is implicit by the order of the vertices of the cell. In particular, we assume
# that vertices are arranged in a clockwise manner. Note that this automatically 
# induces a notion of outward and inward normal on the edge, implicit by the
# ordering of the vertices. 
#
# We adopt the convention of ordering the verices of a cell counter-clockwise. The
# induced orientation of the Edge is as follows then: Imagine standing at p1 and 
# looking towards p2. The outward normal is on the right.
#
# For 2 dimensional meshes we make use of the meshes pacakage (hence no 
# of a 2D mesh type here). However, we still need to declare the appropriate
# types for C(𝕄(Ω)) and E(C) in Ω coordinates.
#

immutable Cell2D
    cid::Int64
    p1::Vertex2
    p2::Vertex2
    p3::Vertex2
end

cid(c::Cell2D) = c.cid

#
# Edges  
#

immutable Edge
    cid::Int64
    p1::Vertex2
    p2::Vertex2
end

cid(c::Edge) = c.cid

#
# Functions to go between edges and their numbering.
#
# Edges are numbered {1,2,3} in a clockwise notion from vertex 0 of the cell. 
# Please also note the disucussion about orientation at the beginning of this
# section since this is where that convention is enforced.
#
function edgeid(x::Cell2D,edge; throw = true)
    if edge.p1 == x.p1 && edge.p2 == x.p2 
        return 1
    elseif edge.p1 == x.p2 && edge.p2 == x.p3
        return 2
    elseif edge.p1 == x.p3 && edge.p2 == x.p1
        return 3
    else
        throw && error("Not an edge of this cell")
        return 0
    end
end
edgeid(m::Meshes.Mesh,edge) = edgeid(m[edge.cid],edge)

function getindex(x::Cell2D,i)
    Edge(x.cid,(if i == 1
        (x.p1,x.p2)
    elseif i == 2
        (x.p2,x.p3)
    elseif i == 3
        (x.p3,x.p1)
    else
        error("Invalid Edge index!")
    end)...)
end

# Iterate over edges over a cell
const edges = (1,2,3)
start(x::Cell2D) = start(edges)
next(x::Cell2D,i) = ((j,state) = next(edges,i); (x[j],state))
done(x::Cell2D,i) = done(edges,i)

# Reverse the orientation of an edge. 
# In particular note that for two adjacent triangles, the two 
# overlapping edges of the adjacent triangles will have opposite orientation.
reverse(e::Edge) = Edge(e.cid,e.p2,e.p1)

# Computes ||⋅||_2 of the edge
norm(e::Edge) = norm(e.p2-e.p1)

#### Iteration over meshes/indexing into meshes

# Create a cell from a mesh face (in the Meshes.Face sense)
cell2d(x,f,id) = Cell2D(id,x.vertices[f.v1],x.vertices[f.v2],x.vertices[f.v3])

length(x::Meshes.Mesh) = length(x.faces)
start(x::Meshes.Mesh) = start(x.faces)
function next(x::Meshes.Mesh,state)
    oldstate = state
    f, state = next(x.faces,state)
    (cell2d(x,f,oldstate),state)
end
done(x::Meshes.Mesh,state) = done(x.faces, state)
getindex(x::Meshes.Mesh,i) = cell2d(x,x.faces[i],i)

Ak(c) = [ (-c.p1[1] + c.p2[1]) (-c.p1[1] + c.p3[1])
           (-c.p1[2] + c.p2[2]) (-c.p1[2] + c.p3[2]) ]
Bk(c) = [ c.p1[1]; c.p1[2] ]

# Should fix multiplication. No time now tough
function 𝜒(c::Cell2D,p)
    iA = inv(Ak(c))
    Vertex2(([ iA[1,1]*p.coords[1] + iA[1,2]*p.coords[2];
       iA[2,1]*p.coords[1] + iA[2,2]*p.coords[2]] - 
       iA*Bk(c))...)
end
𝜒⁻¹(c::Cell2D,p) = p1ϕ1(p)*c.p1 + p1ϕ2(p)*c.p2 + p1ϕ3(p)*c.p3
𝜒⁻¹(mesh::Meshes.Mesh,c::Cell2D,p) = 𝜒⁻¹(c,p)

function 𝜒(c::Edge,p)
    s = c.p2-c.p1
    2*dot(p-c.p1,s)/dot(s,s) - 1
end
𝜒⁻¹(c::Edge,p) = c.p1 + (1/2)*(1+p)*(c.p2-c.p1)

# Mapping from a point on the edge to the point on the correspoding reference elemt
# in essence ∂(cell,edge,p) = 𝜒(cell,𝜒⁻¹(edge,p))
function ∂(cell,edge,p)
    id = edgeid(cell,edge)
    ret = id == 1 ? ((p+1)/2)*Vertex2(1.0,0.0) : 
    id == 2 ? Vertex2(1.0,0.0) + ((p+1)/2)*Vertex2(-1.0,1.0) :
    id == 3 ? Vertex2(0.0,1.0)-((p+1)/2)*Vertex2(0.0,1.0) : error("")
    ret
end

#### Other stuff

to2d(m) = Mesh{Vertex2}([(@assert v[3]==0; Vertex2(v[1],v[2])) for v in m.vertices], copy(m.faces))

export to2d, ∂Ω, ∇I

function ∂Ω(m;ns = computeNeighbors(m))
    edges = AC274_2D.Edge[]
    for c in m
        for edge in c
            if !has_neighbor(m, edge; ns=ns)
                push!(edges, edge)
            end
        end
    end
    edges
end

function ∇I(c,coeffs,p; DP=getDPhi(2))
    #@show AC274_2D.chi(c,p)
    Base.sum([a*dphi(AC274.chi(c,p)) for (a,dphi) in zip(coeffs,DP)])
end

# Aliases for those who don't like unicode
const invchi = 𝜒⁻¹
const chi = 𝜒
