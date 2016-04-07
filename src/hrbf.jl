type HermiteRadialField{N, T, PhiType} <: ScalarField
	alphas::Vector{T}
	betas::Vector{Point{N, T}}
	points::Vector{Point{N, T}}
	phi::PhiType
end

function HermiteRadialField{Dimension, T, PhiType <: BaseTwiceDifferentiableFunction}(points::Vector{Point{Dimension, T}}, normals::Vector{Normal{Dimension, T}}, phi_function::PhiType=XCubed())
	@assert length(points) == length(normals)

	num_points = length(points)

	A = Array{T}(num_points * (Dimension + 1), num_points * (Dimension + 1))
	b = Array{T}(num_points * (Dimension + 1))
	# u = Array{T}(Dimension)
	# v = Array{T}(Dimension)

	for point_index = 1:num_points
		row = (point_index - 1) * (1 + Dimension) + 1
		for k = 1:num_points
			col = (k - 1) * (1 + Dimension) + 1
			u = points[point_index] - points[k]
			# for i = 1:Dimension
			# 	u[i] = points[i, point_index] - points[i, k]
			# end
			# u = points[:,point_index] - points[:,k]
			n = norm(u)
			if n <= 1e-9
				A[row + (0:Dimension), col + (0:Dimension)] = 0
			else
				# f = phi_function.f(n)
				# df = phi_function.df(n)
				# ddf = phi_function.ddf(n)
				f = phi(phi_function, n)
				df = dphi(phi_function, n)
				ddf = ddphi(phi_function, n)
				df_over_n = df / n
				v = df_over_n * u
				# for i = 1:Dimension
				# 	v[i] = df_over_n * u[i]
				# end

				A[row, col] = f
				for i = 1:Dimension
					A[row, col + i] = v[i]
					A[row + i, col] = v[i]
				end
				scaling = (ddf - df_over_n) / (n^2)
				for i = 1:Dimension
					for j = 1:Dimension
						A[row + i, col + j] = scaling * u[i] * u[j]
					end
				end
				# A[row + (1:Dimension), col + (1:Dimension)] = (ddf - df_over_n) / (n^2) * (u * u')
				for i = 1:Dimension
					A[row + i, col + i] += df_over_n
				end
			end
		end

		b[row] = 0
		for i = 1:Dimension
			b[row + i] = normals[point_index][i]
		end
	end

	y = A \ b
	y = reshape(y, Dimension + 1, num_points)
	alphas = vec(y[1,:])
	betas = [Point{Dimension, T}(y[2:end, i]) for i in 1:num_points]
	HermiteRadialField{Dimension, T, PhiType}(alphas, betas, points, phi_function)
end

function evaluate{Dimension, T, PhiType}(field::HermiteRadialField{Dimension, T, PhiType}, x::Point{Dimension, T})
	value::T = zero(T)
	num_points = length(field.points)
	@assert length(x) == Dimension
	@assert length(field.alphas) == num_points
	@assert length(field.betas) == length(field.points)
	@inbounds for i = 1:length(field.points)
		u = x - field.points[i]
		n = norm(u)
		if n > 0
			z::T = dot(field.betas[i], u)
			# p = field.phi.f(n)
			# dp = field.phi.df(n)
			p = phi(field.phi, n)
			dp = dphi(field.phi, n)
			value += field.alphas[i] * p
			value += dp / n * z
		end
	end
	value
end

evaluate{Dimension, T}(field::HermiteRadialField{Dimension, T}, x) = evaluate(field, convert(Point{Dimension, T}, x))

function grad{Dimension, T, PhiType}(field::HermiteRadialField{Dimension, T, PhiType}, x::Point{Dimension, T})
	num_points = length(field.points)
	@assert length(x) == Dimension
	@assert length(field.alphas) == num_points
	@assert length(field.betas) == length(field.points)
	g = Point{Dimension, T}(0)
	for i = 1:num_points
		u = x - field.points[i]
		n = norm(u)

		if n > 1e-5
			uhat = u ./ n
			# df = field.phi.df(n)
			# ddf = field.phi.ddf(n)
			df = dphi(field.phi, n)
			ddf = ddphi(field.phi, n)
			alpha_df = field.alphas[i] * df
			beta_uhat = dot(field.betas[i], uhat)

			g += alpha_df .* uhat + beta_uhat * (ddf * uhat - u * df / n^2) + field.betas[i] * df / n
		end
	end
	g
end

grad{Dimension, T}(field::HermiteRadialField{Dimension, T}, x) = grad(field, convert(Point{Dimension, T}, x))

function grad{Dimension, T}(field::HermiteRadialField{Dimension, T})
	FunctionalVectorField{Dimension, T}(x -> grad(field, x))
end

function bounds(field::HermiteRadialField)
	lb = Vec(minimum(field.points))
	ub = Vec(maximum(field.points))
	HyperRectangle(lb, ub - lb)
end

