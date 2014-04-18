# 1D quadrature

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

do_quad_ref(p::DG1D,f::Function) = integrate(f,-1,1)

# ElemJacobian
elemJ(x::DG1D,k) = (1//2)*(x.mesh[k,:r] - x.mesh[k,:l])

# 2D quadrature

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
    else 
        error("Unimplemented")
    end
end

function do_quad(f,porder)
    result = zero(typeof(f(zero(Vertex2))))
    for (x,w) in zip(quadp(porder)...)
        result += w*f(x)
    end
    result
end

do_quad_ref(f::Function,p::DG2D) = do_quad(f,clamp(porder(p)*2,1,4))
do_quad_ref(f::Function,p::DG2D,edge::Edge) = do_quad_ref(p.mesh,edge,f,porder(p)*2)

do_quad(c::Cell2D, f, p) = -det(Ak(c))*do_quad(x->f(𝜒⁻¹(c,x)),p)
elemJ(p::DG2D,c) = -det(Ak(c))

#point = 𝜒⁻¹(edge,x)
#@show p

macro do_quad_ref2d_Δ1(fargs, var, block)
    esc(quote
        args = $fargs
        p = args[1]
        edge = args[2]
        points = args[3]
        gx, gw = points
        # Manual loop peeling
        x,w = gx[1], gw[1]
        $var = x
        result = w*begin
            $block
        end
        for i in 2:length(gw)
            x,w = gx[i],gw[i]
            $var = x
            result += w*begin
                $block
            end
        end
        (norm(edge)/2)*result   
    end)
end

function do_quad_ref(m::Meshes.Mesh, edge::Edge, f, porder; points=Base.gauss(Float64,porder))
    gx, gw = points
    result = 0.0
    for i in 1:length(gw)
        x,w = gx[i],gw[i]
        result += w*f(x)
    end
    (norm(edge)/2)*result
end