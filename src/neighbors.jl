# Neighbors in 1D

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

# Neighbors in 2D

has_neighbor(m::Meshes.Mesh, edge::Edge; ns = computeNeighbors(m)) = has_neighbor(m,edge,ns)
has_neighbor(m::Meshes.Mesh, cell::Cell2D, edgeid::Int64; ns = computeNeighbors(m)) = has_neighbor(m,cell,edgeid,ns)
function has_neighbor(m::Meshes.Mesh, edge::Edge, ns)
    ns[edge.cid, edgeid(m,edge)] != -1
end

function has_neighbor(m::Meshes.Mesh, cell::Cell2D, edgeid::Int64, ns)
    ns[cell.cid, edgeid] != -1
end

has_neighbor(m::Meshes.Mesh, cell::Cell2D, edge::Edge, ns) = has_neighbor(m,edge, ns)
has_neighbor(m::Meshes.Mesh, cell::Cell2D, edge::Edge; ns = computeNeighbors(m)) = 
    has_neighbor(m,edge, ns)

function neighbor(m::Meshes.Mesh, cell::Cell2D, edgeid; ns = computeNeighbors(m))
    m[ns[cell.cid, edgeid]]
end

neighbor(m::Meshes.Mesh, cell::Cell2D, edge::Edge; ns = computeNeighbors(m)) =
    neighbor(m,cell,edge,ns)

function neighbor(m::Meshes.Mesh, cell::Cell2D, edge::Edge, ns)
    m[ns[cell.cid, edgeid(cell,edge)]]
end

function neighbor(m::Meshes.Mesh, edge::Edge, ns = computeNeighbors(m))
    m[ns[cid(edge), edgeid(m[cid(edge)],edge)]]
end

function oppface(m::Mesh, cell::Cell2D, face, ns) 
    n = neighbor(m, cell,face, ns)
    n[edgeid(n,reverse(face))]
end
oppface(m::Mesh, cell::Cell2D, face;  ns = computeNeighbors(m)) = oppface(m,cell,face,ns)

canon(a,b) = (a<b) ? (a,b) : (b,a)

function computeNeighbors(m)
    fs = ((Int64,Int64)=>(Int64,Int64))[]
    ns = Array(Int64,length(m.faces),3)
    fill!(ns,-1)
    for i in 1:length(m.faces)
        f = m.faces[i]
        for (nn,a,b) in ((1,:v1,:v2),(2,:v2,:v3),(3,:v3,:v1))
            c = canon(f.(a),f.(b))
            if haskey(fs,c)
                (j,num) = fs[c]
                ns[i,nn] = j
                ns[j,num] = i
            else
                fs[c] = (i,nn)
            end
        end
    end
    ns
end

# Tests if the vertex p is contained in cell x
function facingface(x::Cell2D,p::Vertex2)
    v0 = x.p2 - x.p1
    v1 = x.p3 - x.p1
    v2 = p - x.p1

    dot00 = dot(v0, v0)
    dot01 = dot(v0, v1)
    dot02 = dot(v0, v2)
    dot11 = dot(v1, v1)
    dot12 = dot(v1, v2)

    invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
    λ1 = u = (dot11 * dot02 - dot01 * dot12) * invDenom
    λ2 = v = (dot00 * dot12 - dot01 * dot02) * invDenom
    λ3 = 1-u-v

    pt = clamp(λ1,0.0,1.0)*x.p1 + clamp(λ2,0.0,1.0)*x.p2 + clamp(λ3,0.0,1.0)*x.p3

    dist = norm(p-pt)

    if (u >= 0) && (v >= 0) && (u + v < 1)
        return (0,0.0)
    elseif (u <= 0) && (v <= 0)
        return ((abs(u) > abs(v) ? 3 : 2),dist)
    elseif (u <= 0)
        return (3,dist)
    elseif (v <= 0)
        return (2,dist)
    else
        return (1,dist)
    end
end

contains(x::Cell2D,p::Vertex2) = facingface(x,p)[1] == 0