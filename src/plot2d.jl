using Cairo

minx(c) = min(c.p1[1],c.p2[1],c.p3[1])
maxx(c) = max(c.p1[1],c.p2[1],c.p3[1])
miny(c) = min(c.p1[2],c.p2[2],c.p3[2])
maxy(c) = max(c.p1[2],c.p2[2],c.p3[2])


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

function drawVertices(m::Meshes.Mesh{Vertex2},vs::Vector{Vertex2},w,h; colors = nothing, draw_mesh = false)
    if colors != nothing 
        @assert length(colors) == length(vs)
    end

    c,cr = setupDrawing(w,h)

    transform = transformf(m,w,h)

    set_source_rgb(cr,0.0,0.0,0.0);    # white

    draw_mesh && _drawMesh(cr,m,w,h,transform=transform)

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
        m = e.p1 + (e.p2-e.p1)./2.
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
        p=FramedPlot()
        setattr(p.x, draw_ticks=false)
        setattr(p.y1, draw_ticks=false)
        setattr(p.x1, draw_ticklabels=false)
        setattr(p.y1, draw_ticklabels=false)
        setattr(p.y2, draw_ticklabels=true) 
 
        xr=(1,2)
        yr=(dmin,dmax)
 
        setattr(p, :xrange, xr)
        setattr(p.y1, range = yr)
        setattr(p.y2, range = yr)

        y=linspace(dmin, dmax, 256)*1.0
        data=[y y]
    elseif orientation == "horizontal"
        p=FramedPlot()
        setattr(p.y, draw_ticks=false)
        setattr(p.x1, draw_ticks=false)
        setattr(p.y1, draw_ticklabels=false)
        setattr(p.x1, draw_ticklabels=false)
        setattr(p.x2, draw_ticklabels=true) 
 
        yr=(1,2)
        xr=(dmin, dmax)

        setattr(p, :yrange, yr)
        setattr(p.x1, range = xr)
        setattr(p.x2, range = xr)        
 
        x=linspace(dmin,dmax,256)*1.0
        data=[x', x']
    end
 
    #setattr(p; kvs...)
 
    clims = (minimum(data),maximum(data))
 
    img = Winston.data2rgb(data, clims, colormap)
    add(p, Image(xr, yr, img))
    p
end

function default_colorize(data)
    colors = Winston.data2rgb(data, (minimum(data),maximum(data)), Winston._current_colormap)
    for x=1:size(data,1), y=1:size(data,2)
        if isnan(data[x,y])
            colors[x,y] = 0xD3D3D3
        end
    end

    c = colorbar(maximum(data),minimum(data), orientation = "vertical")

    colors, c
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
            y=floor(tocry(miny(cell),h,bbox)):ceil(tocry(maxy(cell),h,bbox))
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
    c,cr = setupDrawing(w+(func !== nothing ? 10+ceil(h/10) : 0),h)

    if func !== nothing
        #typeof(func(first(m),Vertex2(0,0)))
        data = Array(Float64,w,h)
        fill!(data,NaN)
        pixelwise!(data,m,w,h,func)
        colors, colorbar = colorize(data)
        ccr = Winston.CairoRenderer(c)
        Winston.compose(colorbar, ccr, Base.Graphics.BoundingBox(w+10,w+10+ceil(h/10),0,h))
        Cairo.image(cr,CairoImageSurface(colors,Cairo.FORMAT_RGB24,flipxy=false),0,0,w,h)
    end

    drawMesh && _drawMesh(cr,m,w,h)

    c
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

function plotSolution{N,T}(p::DG2D, Q::Array{Coeffs{N,T},2}, t)
    coeffs = Q[t,:]
    plotMesh(p.mesh,256,256; func=function (c,point)
        evaluate(coeffs[c.cid],AC274.getPhi(porder(p)),chi(c,point))
    end,
    colorize=AC274.default_colorize,
    drawMesh=false)
end

function plotSolution{T}(p::DG2D, QQ::Array{T,3}, t; tfunc = identity, w =256, h=256)
    coeffs = QQ[t,:,:]
    plotMesh(p.mesh,w,h; func=function (c,point)
        evaluate(tfunc(coeffs[1,c.cid,:]),AC274.getPhi(porder(p)),chi(c,point))
    end,
    colorize=AC274.default_colorize,
    drawMesh=false)
end