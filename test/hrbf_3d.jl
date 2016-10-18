using SpatialFields
using Base.Test
using Iterators

points = SVector{3, Float64}[[1; 0; 0], [0; 1; 0], [-1; 0; 0], [0; -1; 0], [0; 0; 1], [0; 0; -1]]
normals = SVector{3, Float64}[[1; -1; 0], [0; 1; 0], [-1; 0; 0], [0; -1; 0], [1; 0; 1], [0; 0; -1]]


field = HermiteRadialField(points, normals)

X = linspace(-2, 2)
Y = linspace(-2, 2)
Z = linspace(-2, 2)

for (i, point) in enumerate(points)
    @test isapprox(field(point), 0, atol=1e-6)
    @test isapprox(SpatialFields.gradient(field)(point), normals[i], atol=1e-6)
end
