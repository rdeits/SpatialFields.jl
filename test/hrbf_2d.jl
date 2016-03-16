using SpatialFields
using Base.Test

function hrbf_2d()
	points = Point{2, Float64}[[1; 0], [0; 1], [-1; 0], [0; -1]]
	# points = [1. 0; 0 1; -1 0; 0 -1]'
	normals = Normal{2, Float64}[[1; 1], [0; 1], [-1; 1], [0; -1]]
	# normals = [1. 1; 0 1; -1 1; 0 -1]'

	field = HermiteRadialField(points, normals)

	X = linspace(-2, 2)
	Y = linspace(-2, 2)
	Z = zeros(length(X), length(Y))

	@time for i = 1:length(X)
	    for j = 1:length(Y)
	        Z[j,i] = evaluate(field, [X[i], Y[j]])
	    end
	end

	for i in 1:size(points, 2)
		@test isapprox(evaluate(field, points[i]), 0, atol=1e-6)
		eps = 1e-4
		nudged_point = points[i] + eps * Point(normals[i])
		@test isapprox(dot(grad(field, points[i]), Point(normals[i]) * eps), evaluate(field, nudged_point), atol=1e-6)
	end
end

hrbf_2d()
