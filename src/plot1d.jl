using Gadfly

function pf(porder,mesh,a,b,Q)
    basis = phi(porder+1)
    function (x)
        k = min(1+div(x-a,(b-a)/K),K)
        poly = evaluate(Q[1,k,:],basis,𝜒(mesh,k,x))
    end
end

function pf(p::DG1D,Q)
    basis = phi(p)
    function (x)
        k = min(1+div(x-p.a,(p.b-p.a)/p.K),p.K)
        poly = evaluate(Q[1,k,:],basis,𝜒(p.mesh,k,x))
    end
end


plotSolution(p::DG1D,Q) = Gadfly.plot(pf(p,Q),p.a,p.b)