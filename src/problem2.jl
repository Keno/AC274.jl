# Do Gauss quadrature by hand
# The routine that computes the weights are at
# https://github.com/JuliaLang/julia/blob/master/base/quadgk.jl
# note that the method is essentially the same as the one used
# in the Matlab example form class, but takes advantage of 
# the special sparsity pattern and allows these points to be 
# computed in arbitrary spaces
(gx,gw) = Base.gauss(Float64, 7)
do_quadrature(f,a,b) = sum([ w*f(x) for (w,x) in zip(gw,gx) ]) 
g(x) = cos(π*x)
data = Array(Array{Float64,1},4)
for p=1:4
    println("p: $p")
    # Pick the polynomial basis
    cheb = AC274.nodes(p+1) 
    basis = AC274.phi(cheb)
    error = zeros(Float64,5)
    i=1
    for K = 10:10:50
        println("K: $K")
        # Regular Mesh of K cells
        mesh = generateMesh(Mesh1D,0.0,2.0,K)

        for cell in mesh
            # Interpolate at the nodes
            ℐg = AC274.evaluate(AC274.Coeffs{p+1}(map(n->g(𝜒⁻¹(mesh,cell,n)),cheb)),basis)
            
            # Do the quadrature on the current node
            
            error[i] += 0.5*(cell.right-cell.left)*AC274.integrate(y->(AC274.fapply(g,𝜒⁻¹(mesh,cell,y))-AC274.fapply(ℐg,y))^2,-1,1)
        end
        error[i] = sqrt(error[i])
        i += 1
    end
    data[p] = error
    #xs = 2./[10:10:50]
    #display(MIME("text/html"),plot(layer(h->h^(p+1),minimum(xs),maximum(xs)),
    #    layer(x=xs,y=error, Geom.point),Scale.x_log10,Scale.y_log10))
end