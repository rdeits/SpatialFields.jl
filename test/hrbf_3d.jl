using SpatialFields
using Base.Test
using Iterators

function hrbf_3d()
	points = SVector{3, Float64}[[1; 0; 0], [0; 1; 0], [-1; 0; 0], [0; -1; 0], [0; 0; 1], [0; 0; -1]]
	normals = SVector{3, Float64}[[1; -1; 0], [0; 1; 0], [-1; 0; 0], [0; -1; 0], [1; 0; 1], [0; 0; -1]]


	field = HermiteRadialField(points, normals)
	@time field = HermiteRadialField(points, normals)

	X = linspace(-2, 2)
	Y = linspace(-2, 2)
	Z = linspace(-2, 2)

	field([X[1], Y[1], Z[1]])
	v = [X[1], Y[1], Z[1]]
	@time for j = 1:1e5; field(v); end
	@time C = [field([x,y,z]) for (x,y,z) in product(X, Y, Z)];
end

hrbf_3d()
