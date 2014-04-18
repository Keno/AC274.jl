import NumericExtensions: add!, multiply!

function globalSolve(p::DG,MMatrix,Q,t,cache)
    k = Array(Coeffs{nbf(p),eltype(eltype(Q))},numcells(p))
    for c in p.mesh
        r = Coeffs(nbf(p),MMatrix[cid(c)]\RHS(p,c,Q,t,cache));
        k[cid(c)] = r
    end
    k
end

function _solve(p::DG, ℚ; cache = generateMatrices(p))
    basis = phi(p)
    # ℚ[1] is the inital condition
    M = [factorize(Mlocal(p,c,basis)) for c in p.mesh]
    for i in 2:p.nt
        k1 = globalSolve(p,M,ℚ[i-1,:],i*p.Δt,cache)
        k2 = globalSolve(p,M,ℚ[i-1,:] + (0.5*p.Δt*k1)',i*p.Δt,cache)
        k3 = globalSolve(p,M,ℚ[i-1,:] + (0.5*p.Δt*k2)',i*p.Δt,cache)
        k4 = globalSolve(p,M,ℚ[i-1,:] + (p.Δt*k3)',i*p.Δt,cache)
        for j = 1:numcells(p)
            # This is probably stupid, but 
            # oh well
            # ℚ[i,:] = ℚ[i-1,:] + p.Δt/6 * (k1+2*k2+2*k3+k4)';
            ℚ[i,j] = copy(ℚ[i-1,j])
            multiply!(k1[j].coeffs,p.Δt/6)
            add!(ℚ[i,j].coeffs,k1[j].coeffs)
            multiply!(k2[j].coeffs,2*p.Δt/6)
            add!(ℚ[i,j].coeffs,k2[j].coeffs)
            multiply!(k3[j].coeffs,2*p.Δt/6)
            add!(ℚ[i,j].coeffs,k3[j].coeffs)
            multiply!(k4[j].coeffs,p.Δt/6)
            add!(ℚ[i,j].coeffs,k4[j].coeffs)
        end
    
        # post-processing
        if isa(p,DG1D) && p.useMUSCL
            @assert p.p == 1 && isa(p,DG1D) # For now
            for c in p.mesh
                Qk = ℚ[i,c.coord]
                ql = evaluate(Qk,basis,coord(:l))
                qr = evaluate(Qk,basis,coord(:r))

                h = (c.right-c.left)

                slope = (qr-ql)/h
                qkbar = (ql+qr)/2

                s = sign(slope)
                minabsa = abs(slope)

                for face in (:l,:r)
                    if has_neighbor(p.mesh,c, face)
                        N = neighbor(p.mesh, c, face)
                        Qk = ℚ[i,N.coord]
                        ql = evaluate(Qk,basis,coord(:l))
                        qr = evaluate(Qk,basis,coord(:r))
                        qNbar = (ql+qr)/2

                        eeslope = n⁻(face)*(qNbar-qkbar)/h

                        minabsa = min(minabsa,abs(eeslope))

                        s += sign(eeslope)
                    end
                end

                s/=3
                # Hack in knowledge about WaveState - will fix after
                # assignment is due

                local q1c, q2c
                if abs(s).q1 != 1
                    q1c = [qkbar.q1 for i=1:(p.p+1)]
                else
                    slope = s*minabsa
                    q1c = [(qkbar+(node)*(h/2)*slope).q1 for node in cache.ns]
                end

                if abs(s).q2 != 1
                    q2c = [qkbar.q2 for i=1:(p.p+1)]
                else
                    slope = s*minabsa
                    q2c = [(qkbar+(node)*(h/2)*slope).q2 for node in cache.ns]
                end

                #@show (s, q1c, q2c)

                cs = typeof(s)[typeof(s)(q1c[i],q2c[i]) for i = 1:length(q1c)]
                ℚ[i,c.coord] = Coeffs(p.p+1,cs)
            end
        end
    end
    ℚ
end

# Nodal interpolation
interpolate(f, nodes) = (Coeffs(length(nodes),map(f,nodes)))

function solve(p::DG)
    if p.mesh === nothing
        if isa(p,DG1D)
            p.mesh = generateMesh(Mesh1D,p.a,p.b,p.K; periodic=false)
        else
            error("Can only generate 1D uniform meshes. Please specify a mesh!")
        end
    end

    cache = generateMatrices(p)

    # Save coefficients for all timesteps (nt)
    # for all K cells
    if isa(p.q₀,Vector)
        T = eltype(p.q₀)
        ℚ = Array(Coeffs{p.p+1,T},p.nt,numcells(p))
        ℚ[1,:] = p.q₀ #[zero(Coeffs{p.p+1}) for i = 1:p.K]
    elseif isa(p.q₀,Function)
        T = typeof(p.q₀(zerop(p)))
        ℚ = Array(Coeffs{nbf(p),T},p.nt,numcells(p))
        # Use nodal interpolation
        ns = nodes(p)
        for c in p.mesh
            ℚ[1,cid(c)] = interpolate(p.q₀,map(n->𝜒⁻¹(p.mesh,c,n),ns))
        end
    elseif isa(p.q₀, Number)
        # Uniform initial conditions
        T = typeof(p.q₀)
        ℚ = Array(Coeffs{p.p+1,T},p.nt,numcells(p))
 
        ℚ[1,:] = p.q₀*ones(Coeffs{p.p+1,T},numcells(p))
    else
        @assert "Unknown option"
    end

    _solve(p,ℚ; cache=cache)
end