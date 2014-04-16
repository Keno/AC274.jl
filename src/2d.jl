#
# Polynomial Basis
# (For now use the ones provided in the lecture notes.
# Might come up with a better one (Chebyshev?) later
#

# Seperate module to test out infrastructure for now
module AC274_2D

using Meshes
using ImmutableArrays

import AC274: DG

import Meshes: Vertex2

type DG2D <: DG
    mesh::Meshes.Mesh{Vertex2}
    fflux::Function # The numerical flux function
    Δt::Float64 # Time-step
    nt::Int64 # Number of time steps
    C::Float64 # Diffusion constant
    q₀
    boundary::Function
end

immutable Cell2D
    cid::Int64
    p1::Vertex2
    p2::Vertex2
    p3::Vertex2
end

immutable Edge
    cid::Int64
    p1::Vertex2
    p2::Vertex2
end

import Base: start, next, done, length, getindex, norm

norm(e::Edge) = norm(e.p2-e.p1)

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

const edges = (:e1,:e2,:e3)
start(x::Cell2D) = start(edges)
next(x::Cell2D,i) = ((j,state) = next(edges,i); (x[j],state))
done(x::Cell2D,i) = done(edges,i)
function getindex(x::Cell2D,i)
    Edge(x.cid,(if i == :e1
        (x.p1,x.p2)
    elseif i == :e2
        (x.p2,x.p3)
    elseif i == :e3
        (x.p3,x.p1)
    end)...)
end

function edgeid(x::Cell2D,edge)
    if edge.p1 == x.p1 && edge.p2 == x.p2
        return 1
    elseif edge.p1 == x.p2 && edge.p2 == x.p3
        return 2
    elseif edge.p1 == x.p3 && edge.p2 == x.p1
        return 3
    else
        error("Not an edge of this cell")
    end
end
edgeid(m::Meshes.Mesh,edge) = edgeid(m[edge.cid],edge)

import AC274: neighbor, has_neighbor

function has_neighbor(m::Meshes.Mesh, edge::Edge; ns = computeNeighbors(m))
    ns[edge.cid, edgeid(m,edge)] != -1
end

function has_neighbor(m::Meshes.Mesh, cell::Cell2D, edgeid; ns = computeNeighbors(m))
    ns[cell.cid, edgeid] != -1
end

function neighbor(m::Meshes.Mesh, cell, face; ns = computeNeighbors(m))
    m[ns[cell.cid, face]]
end

Ak(c) = [ (-c.p1[1] + c.p2[1]) (-c.p1[1] + c.p3[1])
           (-c.p1[2] + c.p2[2]) (-c.p1[2] + c.p3[2]) ]
Bk(c) = [ c.p1[1]; c.p1[2] ]

##
# 
#
##

# Should fix multiplication. No time now tough
function 𝜒(c::Cell2D,p)
    iA = inv(Ak(c))
    Vertex2(([ iA[1,1]*p.coords[1] + iA[1,2]*p.coords[2];
       iA[2,1]*p.coords[1] + iA[2,2]*p.coords[2]] - 
       iA*Bk(c))...)
end
𝜒⁻¹(c::Cell2D,p) = p1ϕ1(p)*c.p1 + p1ϕ2(p)*c.p2 + p1ϕ3(p)*c.p3

function 𝜒(c::Edge,p)
    s = c.p2-c.p1
    proj = (dot(p-c.p1)/dot(s,s))*s
    2*norm(proj)/norm(s) - 1
end
𝜒⁻¹(c::Edge,p) = c.p1 + (1/2)*(1+p)*(c.p2-c.p1)

const chi = 𝜒
const invchi = 𝜒⁻¹

p0ϕ1(ξ,η) = 1
const P0 = [p0ϕ1]
const N0 = [Vertex2(1//3,1//3)]

p1ϕ1(ξ,η) = 1 - ξ - η
p1ϕ2(ξ,η) = ξ
p1ϕ3(ξ,η) = η

∇p1ϕ1(ξ,η) = [-1,-1]
∇p1ϕ2(ξ,η) = [1,0]
∇p1ϕ3(ξ,η) = [0,1]

const P1 = [p1ϕ1,p1ϕ2,p1ϕ3]
const DP1 = [∇p1ϕ1,∇p1ϕ2,∇p1ϕ3]
const N1 = [Vertex2(0,0),Vertex2(1,0),Vertex2(0,1)]

p2ϕ1(ξ,η) = 2*(0.5 - ξ - η)*(1 - ξ - η)
p2ϕ2(ξ,η) = 2ξ*(ξ - 0.5)
p2ϕ3(ξ,η) = 2η*(η - 0.5)
p2ϕ4(ξ,η) = 4ξ*(1 - ξ - η)
p2ϕ5(ξ,η) = 4η*ξ
p2ϕ6(ξ,η) = 4η*(1 - ξ - η)

∇p2ϕ1(ξ,η) = [4ξ+4η-3,4ξ+4η-3]
∇p2ϕ2(ξ,η) = [4ξ - 1,0]
∇p2ϕ3(ξ,η) = [0,4η - 1]
∇p2ϕ4(ξ,η) = [-8ξ-4η+4,-4ξ]
∇p2ϕ5(ξ,η) = [4η,4ξ]
∇p2ϕ6(ξ,η) = [-4η,-8η-4ξ+4]

const P2 = [p2ϕ1,p2ϕ2,p2ϕ3,p2ϕ4,p2ϕ5,p2ϕ6]
const DP2 = [∇p2ϕ1,∇p2ϕ2,∇p2ϕ3,∇p2ϕ4,∇p2ϕ5,∇p2ϕ6]
const N2 = [Vertex2(0,0),Vertex2(1,0),Vertex2(0,1),Vertex2(0.5,0),Vertex2(0.5,0.5),Vertex2(0,0.5)]

for f in [P0,P1,P2,DP1,DP2]
    @eval $(f.env.name)(v::Vertex2) = $(f.env.name)(v[1],v[2])
end

getPhi(p) = (@assert 0 <= p <= 2; p==0 ? P0 : p == 1 ? P1 : P2)
getDPhi(p) = (@assert 0 <= p <= 2; p==0 ? DP0 : p == 1 ? DP1 : DP2)
nodes(p) = (@assert 0 <= p <= 2; p==0 ? N0 : p == 1 ? N1 : N2)

# Get quadrature weights (from lecture notes)
function quadp(p) 
    if p == 1
        return (
            [Vertex2(1//3,1//3)],[0.5])
    elseif p == 2
        return (
            [Vertex2(2//3,1//6),Vertex2(1//6,2//3),
             Vertex2(1//6,1//6)],
             [1//6,1//6,1//6])
    elseif p == 3
        return (
            [Vertex2(0.1550510257, 0.1785587282),
             Vertex2(0.6449489742, 0.0750311102),
             Vertex2(0.1550510257, 0.6663902460),
             Vertex2(0.6449489742, 0.2800199154)],
             [0.1590206908,0.0909793091,
              0.1590206908,0.0909793091])
    elseif p == 4
        return (
            [Vertex2(0.4459484909, 0.4459484909),
             Vertex2(0.1081030181, 0.4459484909),
             Vertex2(0.4459484909, 0.1081030181),
             Vertex2(0.0915762135, 0.0915762135),
             Vertex2(0.8168475729, 0.0915762135),
             Vertex2(0.0915762135, 0.8168475729)],
             [0.1116907948,0.1116907948,0.1116907948,
              0.0549758718,0.0549758718,0.0549758718])
    end
end

function do_quad(f,p)
    result = 0.0
    for (x,w) in zip(quadp(p)...)
        result += w*f(x)
    end
    result
end

do_quad(c::Cell2D, f, p) = -det(Ak(c))*do_quad(x->f(𝜒⁻¹(c,x)),p)

function do_quad(m::Meshes.Mesh, edge::Edge, f, porder; points=Base.gauss(Float64,porder))
    gx, gw = points
    result = 0.0
    for i in 1:length(gw)
        x,w = gx[i],gw[i]
        point = 𝜒⁻¹(edge,x)
        #@show p
        result += w*f(point)
    end
    (norm(edge)/2)*result
end

export n⃗, midp

function n⃗(f::Edge) 
    t⃗ = f.p1-f.p2
    t⃗ /= norm(t⃗)
    Vector2{Float64}(t⃗[2],-t⃗[1])
end

midp(f::Edge) = (1/2)*(f.p1+f.p2)

# Plotting things

using Cairo

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

minx(c) = min(c.p1[1],c.p2[1],c.p3[1])
maxx(c) = max(c.p1[1],c.p2[1],c.p3[1])
miny(c) = min(c.p1[2],c.p2[2],c.p3[2])
maxy(c) = max(c.p1[2],c.p2[2],c.p3[2])

import Base: *, -, +, dot, ./

for f in (:(-),:(+),:dot)
    @eval ($f)(a::Vertex2,b::Vertex2) = ($f)(a.coords,b.coords)
end
*(a::Float64,b::Vertex2) = Vertex2(a*b.coords)
+(a::Vector2{Float64},b::Vertex2) = Vertex2(a+b.coords)
-(a::Vector2{Float64},b::Vertex2) = Vertex2(a-b.coords)
+(b::Vertex2,a::Vector2{Float64}) = a+b
./(a::Vertex2,b::Number) = Vertex2(a.coords./b)

function transformf(m,w,h)
    xs = map(v->v[1],m.vertices)
    ys = map(v->v[2],m.vertices)
    mminx, mmaxx = extrema(xs)
    mminy, mmaxy = extrema(ys)

    function transform(v)
        typeof(v)(tocr(v[1],w,mminx,mmaxx),tocr(v[2],h,mminy,mmaxy))
    end
end

export vertices

vertices(c::Cell2D) = [c.p1,c.p2,c.p3]

function drawVertices(m::Meshes.Mesh{Vertex2},vs::Vector{Vertex2},w,h; colors = nothing)
    if colors != nothing 
        @assert length(colors) == length(vs)
    end

    c,cr = setupDrawing(w,h)

    transform = transformf(m,w,h)

    set_source_rgb(cr,0.0,0.0,0.0);    # white

    for i in 1:length(vs)
        v = vs[i]
        colors != nothing && set_source(cr,colors[i])
        vʹ = transform(v)
        circle(cr,vʹ[1],vʹ[2],2.0)
        fill(cr)
    end

    c
end

function setupDrawing(w,h)
    c = CairoRGBSurface(w,h);
    cr = CairoContext(c);

    save(cr);
    set_source_rgb(cr,1.0,1.0,1.0);    # white
    rectangle(cr,0.0,0.0,w,h); # background
    fill(cr);
    restore(cr);

    scale(cr,1.0,-1.0)
    translate(cr,0.0,-h)

    (c,cr)
end

function _drawMesh(cr::CairoContext, m::Meshes.Mesh{Vertex2}, w, h; transform = transformf(m,w,h))
    tvs = map(transform,m.vertices)

    set_source_rgb(cr,0.0,0.0,0.0);    # black
    set_line_width (cr, 1.0);
    for f in m.faces
        s = tvs[f.v1][1:2]
        move_to(cr,s...)
        line_to(cr,tvs[f.v2][1:2]...)
        line_to(cr,tvs[f.v3][1:2]...)
        line_to(cr,s...)
    end    
    stroke(cr);
end
function drawMesh(m::Meshes.Mesh{Vertex2}, w, h)
    c,cr = setupDrawing(w,h)
    _drawMesh(cr,m,w,h)
    c
end

function drawNormals(m::Meshes.Mesh{Vertex2},edges::Vector{Edge},w,h; inward=false)
    c,cr = setupDrawing(w,h)

    transform = transformf(m,w,h)

    _drawMesh(cr,m,w,h,transform=transform)

    set_source_rgb(cr,1.0,0.0,0.0)
    for e in edges
        m = e.p1 + (e.p2-e.p1)/2
        n = inward ? -n⃗(e)/20 : n⃗(e)/20
        endp = m+n
        endp = transform(Vertex2(endp[1],endp[2]))
        m = transform(m)
        move_to(cr,m[1],m[2])
        line_to(cr,endp[1],endp[2])
    end
    stroke(cr)

    c
end

using Winston

#colorbar
function colorbar(dmin, dmax; orientation="horizontal", colormap=Winston._current_colormap, kvs...)
 
    if orientation == "vertical"
        p=FramedPlot(aspect_ratio=10.0)
        setattr(p.x, draw_ticks=false)
        setattr(p.y1, draw_ticks=false)
        setattr(p.x1, draw_ticklabels=false)
        setattr(p.y1, draw_ticklabels=false)
        setattr(p.y2, draw_ticklabels=true) 
 
        xr=(1,2)
        yr=(dmin,dmax)
 
        y=linspace(dmin, dmax, 256)*1.0
        data=[y y]
    elseif orientation == "horizontal"
        p=FramedPlot(aspect_ratio=0.1)
        setattr(p.y, draw_ticks=false)
        setattr(p.x1, draw_ticks=false)
        setattr(p.y1, draw_ticklabels=false)
        setattr(p.x1, draw_ticklabels=false)
        setattr(p.x2, draw_ticklabels=true) 
 
        yr=(1,2)
        xr=(dmin, dmax)
 
        x=linspace(dmin,dmax,256)*1.0
        data=[x', x']
    end
 
    setattr(p, :xrange, xr)
    setattr(p, :yrange, yr)
 
    setattr(p; kvs...)
 
    clims = (minimum(data),maximum(data))
 
    img = Winston.data2rgb(data, clims, colormap)
    add(p, Image(xr, yr, img))
    p
end

function default_colorize(data)
    colors = Winston.data2rgb(data, (minimum(data),maximum(data)+1), Winston._current_colormap)
    for x=1:size(data,1), y=1:size(data,2)
        if isnan(data[x,y])
            colors[x,y] = 0xD3D3D3
        end
    end

    #display(colorbar(minimum(data),maximum(data)))

    colors
end

immutable BBox 
    mminx::Float64
    mmaxx::Float64
    mminy::Float64
    mmaxy::Float64
end

tocr(x,s,mmin,mmax) = s*(x-mmin)/(mmax-mmin)
tom(x,s,mmin,mmax) = mmin+(x/s)*(mmax-mmin)

tocrx(x,s,bbox::BBox) = tocr(x,s,bbox.mminx,bbox.mmaxx)
tocry(x,s,bbox::BBox) = tocr(x,s,bbox.mminy,bbox.mmaxy)

tomx(x,s,bbox::BBox) = tom(x,s,bbox.mminx,bbox.mmaxx)
tomy(x,s,bbox::BBox) = tom(x,s,bbox.mminy,bbox.mmaxy)

tom(c::Cell2D,x,y,w,h,bbox::BBox) = Vertex2(tomx(x,w,bbox),tomy(y,h,bbox))
tocr(c::Cell2D,x,y,w,h,bbox::BBox)= Vertex2(tocrx(x,w,bbox),tocry(y,h,bbox))

function computebbox(m)
    xs = map(v->v[1],m.vertices)
    ys = map(v->v[2],m.vertices)
    mminx, mmaxx = extrema(xs)
    mminy, mmaxy = extrema(ys)
    BBox(mminx,mmaxx,mminy,mmaxy)
end

function pixelwise!(data,m::Meshes.Mesh{Vertex2},w,h,func; bbox=computebbox(m), naval=NaN)
    ns = computeNeighbors(m)
    for i in 1:length(m.faces)
        f = m.faces[i]
        # Find all the pixels whose center is in this triangle. 
        # There's probably a smarter way to do this, but this is fine
        # for now.
        #@show i
        cell = cell2d(m,f,i)
        for x=floor(tocrx(minx(cell),w,bbox)):ceil(tocrx(maxx(cell),w,bbox)),
            y=floor(tocry(miny(cell),h,bbox)):ceil(tocrx(maxy(cell),h,bbox))
            p = Vertex2(tomx(x,w,bbox),tomy(y,h,bbox))
            ff, dist = facingface(cell,p)
            if ff == 0
                data[clamp(x,1,w),clamp(y,1,h)] = func(cell,p)
            else 
                #@show dist
                mindist = Inf
                for face in 1:3
                    has_neighbor(m,cell,face;ns=ns) || continue
                    nn = neighbor(m,cell,face;ns=ns)
                    _, ndist = facingface(nn,p)
                    mindist = min(ndist,mindist)
                    #@show mindist
                    contains(nn,p) && assert(mindist == 0.)
                    mindist == 0. && break
                end
                if mindist >= dist-w*h*eps() && data[clamp(x,1,w),clamp(y,1,h)] === naval
                    data[clamp(x,1,w),clamp(y,1,h)] = func(cell,p)
                end
            end
        end
    end
end

function plotMesh(m::Meshes.Mesh{Vertex2},w,h; func=nothing, colorize=identity, drawMesh=true)
    c,cr = setupDrawing(w,h)

    if func !== nothing
        #typeof(func(first(m),Vertex2(0,0)))
        data = Array(Float64,w,h)
        fill!(data,NaN)
        pixelwise!(data,m,w,h,func)
        Cairo.image(cr,CairoImageSurface(colorize(data),Cairo.FORMAT_RGB24,flipxy=false),0,0,w,h)
    end

    drawMesh && _drawMesh(cr,m,w,h)

    c
end

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
    Base.sum([a*dphi(AC274_2D.chi(c,p)) for (a,dphi) in zip(coeffs,DP)])
end

function drawArrow(cr, from, to; arrow_length = 15.0, arrow_degrees = 0.5)
    angle = atan2(to[2] - from[2], to[1] - from[1]) + pi;
    p1 = to + Vertex2(arrow_length * cos(angle - arrow_degrees),
                arrow_length * sin(angle - arrow_degrees))
    p2 = to + Vertex2(arrow_length * cos(angle + arrow_degrees),
                arrow_length * sin(angle + arrow_degrees))
    set_source_rgb(cr,1.0,0.0,0.0)
    move_to(cr,from[1],from[2])
    line_to(cr,to[1],to[2])
    line_to(cr,p1[1],p1[2])
    move_to(cr,to[1],to[2])
    line_to(cr,p2[1],p2[2])
    stroke(cr)
end

norm(v::Vertex2) = norm(v.coords)

function draw_vector_field(m,Q,w,h;nump=10)
    c,cr = setupDrawing(w,h)

    #_drawMesh(cr,m,w,h)

    bbox = computebbox(m)

    xps = map(round,linspace(0,w,nump))
    yps = map(round,linspace(0,h,nump))

    data = zeros(Int64,w,h)
    pixelwise!(data,m,w,h,(cell,p)->cell.cid; bbox=bbox)

    for x in xps, y in yps
        cid = data[clamp(x,1,w),clamp(y,1,h)]
        if cid==0
            continue
        end
        p = tom(m[cid],x,y,w,h,bbox)

        gf = ∇I(m[cid],Q[cid],p)

        #v = tocr(m[cid],gf[1],gf[2],w,h,bbox)
        v = (Ak(m[cid])')\gf
        #v = gf
        v ./= norm(v)
        v = 20.0*v
        #
        drawArrow(cr,Vertex2(x,y),Vertex2(v[1]+x,v[2]+y); arrow_length = 5)
    end

    c
end

##
# f = a f_1 + b f_2 + c f_3
# g = f(chi(p))
# g' = f'(chi(p))*chi'(p)
## 

function drawref(w,h,m,cell,coeffs,func; colorize=identity, nump=10)
    c,cr = setupDrawing(w,h)
    data = Array(Float64,w,h)
    fill!(data,NaN)
    bbox = computebbox(m)
    for x=1:w, y=1:h
        point = cell.p1 + (x/w)*(cell.p2-cell.p1) + (y/h)*(cell.p3-cell.p1)
        if contains(cell,point)
            data[clamp(x,1,w),clamp(y,1,h)] = func(cell,point)
        end
    end

    Cairo.image(cr,CairoImageSurface(colorize(data),Cairo.FORMAT_RGB24,flipxy=false),0,0,w,h)

    xps = map(round,linspace(0,w,nump))
    yps = map(round,linspace(0,h,nump))

    for x in xps, y in yps
        point = cell.p1 + (x/w)*(cell.p2-cell.p1) + (y/h)*(cell.p3-cell.p1)
        if contains(cell,point)
            gf = ∇I(cell,coeffs,point)

            #v = tocr(m[cid],gf[1],gf[2],w,h,bbox)
            #v = 𝜒⁻¹(m[cid],Vertex2(gf[1],gf[2]))-m[cid].p1
            v = gf
            v ./= norm(v)
            v = 20.0*v
            #
            drawArrow(cr,Vertex2(x,y),Vertex2(v[1]+x,v[2]+y); arrow_length = 5)
        end
    end

    c
end

end