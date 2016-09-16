using SpatialFields
using Base.Test
# using PyPlot

function hrbf_2d_5th_power()
	points = SVector{2, Float64}[[1; 0], [0; 1], [-1; 0], [0;-1]]
	normals = SVector{2, Float64}[[1; 1], [0; 1], [-1; 1], [0; -1]]
	num_points = size(points, 2)

	field = HermiteRadialField(points, normals, x -> x^5)

	X = linspace(-2, 2)
	Y = linspace(-2, 2)
	Z = zeros(length(X), length(Y))

	@time for i = 1:length(X)
	    for j = 1:length(Y)
	        Z[j,i] = field([X[i], Y[j]])
	    end
	end

	for i in 1:size(points, 2)
		@test isapprox(field(points[i]), 0, atol=1e-6)
		eps = 1e-4
		nudged_point = points[i] + eps * SVector(normals[i])
		@test isapprox(dot(SpatialFields.gradient(field)(points[i]), SVector(normals[i])) * eps, field(nudged_point), atol=1e-6)
	end

	# clf()
	# hold(true)
	# PyPlot.contour(X, Y, Z, [-0.1, 0.0, 0.1])
	# for i = 1:num_points
	#     PyPlot.quiver(points[:,i]..., normals[:,i]...)
	# end
	# axis("equal")
	# show()
end

hrbf_2d_5th_power()
