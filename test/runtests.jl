using Base.Test

mesh = generateMesh(Mesh1D, 0, 5, 50)

ns = nodes(4)
𝜒 = atlas(mesh)
𝜒⁻¹ = invatlas(mesh)
k = 25
@test_approx_eq map(x->𝜒(k,𝜒⁻¹(k,x)),ns) ns

# 2D tests

for p,n in zip(AC274_2D.P1,AC274_2D.N1)
    @test fapply(p,n) == 1
end

for p,n in zip(AC274_2D.P2,AC274_2D.N2)
    @test fapply(p,n) == 1
end

for p = (1,2,3,4)
    @test AC274_2D.do_quad(x->1,p) == 0.5
end