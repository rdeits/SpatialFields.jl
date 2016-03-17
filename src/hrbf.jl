
abstract BaseTwiceDifferentiableFunction

type TwiceDifferentiableFunction <: BaseTwiceDifferentiableFunction
	f::Function
	df::Function
	ddf::Function
end
phi(func::TwiceDifferentiableFunction, x) = func.f(x)
dphi(func::TwiceDifferentiableFunction, x) = func.df(x)
ddphi(func::TwiceDifferentiableFunction, x) = func.ddf(x)

function TwiceDifferentiableFunction(f::Function)
	TwiceDifferentiableFunction(f,
		x -> ForwardDiff.derivative(f, x),
		x -> ForwardDiff.derivative(y -> ForwardDiff.derivative(f, y), x))
end

type PhiXCubed <: BaseTwiceDifferentiableFunction
end
phi(::PhiXCubed, x) = x^3
dphi(::PhiXCubed, x) = 3x^2
ddphi(::PhiXCubed, x) = 6x

type HermiteRadialField{N, T} <: ScalarField
	alphas::Vector{T}
	betas::Vector{Point{N, T}}
	points::Vector{Point{N, T}}
	phi::BaseTwiceDifferentiableFunction
end

function HermiteRadialField{Dimension, T}(points::Vector{Point{Dimension, T}}, normals::Vector{Normal{Dimension, T}}, phi_function::BaseTwiceDifferentiableFunction=PhiXCubed())
	@assert length(points) == length(normals)

	if any(n -> any(isnan, n), normals)
		error("Should not get nans as normals anymore")
		valid_points = !map(n -> any(isnan, n), normals)
		return HermiteRadialField(points[valid_points], normals[valid_points], phi_function)
	end
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
			if n == 0
				A[row + (0:Dimension), col + (0:Dimension)] = 0
			else
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
	if any(isnan(y))
		@show A
		@show b
		@show y
		error("got nans after linear system solve")
	end
	alphas = vec(y[1,:])
	betas = [Point{Dimension, T}(y[2:end, i]) for i in 1:num_points]
	# betas = y[2:end,:]
	HermiteRadialField{Dimension, T}(alphas, betas, points, phi_function)
end

function evaluate{Dimension, T}(field::HermiteRadialField{Dimension, T}, x::Point{Dimension, T})
	value::T = zero(T)
	num_points = length(field.points)
	@assert length(x) == Dimension
	@assert length(field.alphas) == num_points
	@assert length(field.betas) == length(field.points)
	for i = 1:length(field.points)
		u = x - field.points[i]
		# for j = 1:dimension
		# 	u[j] = x[j] - field.points[j, i]
		# end
		n = norm(u)
		if n > 0
			value += field.alphas[i] * phi(field.phi, n) + dphi(field.phi, n) / n * dot(field.betas[i], u)
		end
	end
	value
end

evaluate{Dimension, T}(field::HermiteRadialField{Dimension, T}, x) = evaluate(field, convert(Point{Dimension, T}, x))

function grad{Dimension, T}(field::HermiteRadialField{Dimension, T}, x::Point{Dimension, T})
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
