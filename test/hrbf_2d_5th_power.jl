using SpatialFields
using Base.Test
# using PyPlot

function hrbf_2d_5th_power()
	points = [1. 0; 0 1; -1 0; 0 -1]'
	normals = [1. 1; 0 1; -1 1; 0 -1]'
	num_points = size(points, 2)

	field = HermiteRadialField(points, normals, SpatialFields.TwiceDifferentiableFunction(x -> x^5))

	X = linspace(-2, 2)
	Y = linspace(-2, 2)
	Z = zeros(length(X), length(Y))

	@time for i = 1:length(X)
	    for j = 1:length(Y)
	        Z[j,i] = evaluate(field, [X[i], Y[j]])
	    end
	end

	for i in 1:size(points, 2)
		@test isapprox(evaluate(field, points[:,i]), 0, atol=1e-6)
		eps = 1e-4
		nudged_point = points[:,i] + eps * normals[:,i]
		@test isapprox((grad(field, points[:,i])' * normals[:,i] * eps)[1], evaluate(field, nudged_point), atol=1e-6)
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