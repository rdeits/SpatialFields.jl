
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

type HermiteRadialField{T} <: ScalarField
	alphas::Vector{T}
	betas::Array{T, 2}
	points::Array{T, 2}
	phi::BaseTwiceDifferentiableFunction
end

function HermiteRadialField{T}(points::Array{T, 2}, normals::Array{T, 2}, phi_function::BaseTwiceDifferentiableFunction=PhiXCubed())
	@assert size(points) == size(normals)
	dimension = size(points, 1)
	num_points = size(points, 2)

	A = Array{T}(num_points * (dimension + 1), num_points * (dimension + 1))
	b = Array{T}(num_points * (dimension + 1))

	for point_index = 1:num_points
		row = (point_index - 1) * (1 + dimension) + 1
		for k = 1:num_points
			col = (k - 1) * (1 + dimension) + 1
			u = vec(points[:,point_index] - points[:,k])
			n = norm(u)
			if n == 0
				A[row + (0:dimension), col + (0:dimension)] = 0
			else
				f = phi(phi_function, n)
				df = dphi(phi_function, n)
				ddf = ddphi(phi_function, n)
				df_over_n = df / n
				v = df_over_n .* u

				A[row, col] = f
				A[row, col + (1:dimension)] = v
				A[row + (1:dimension), col] = v
				A[row + (1:dimension), col + (1:dimension)] = (ddf - df_over_n) / (n^2) * (u * u')
				for i = 1:dimension
					A[row + i, col + i] += df_over_n
				end
			end
		end

		b[row] = 0
		b[row + (1:dimension)] = normals[:,point_index]
	end

	y = A \ b
	y = reshape(y, dimension + 1, num_points)
	alphas = y[1,:]
	betas = y[2:end,:]
	HermiteRadialField{T}(vec(alphas), betas, points, phi_function)
end

function evaluate{T}(field::HermiteRadialField{T}, x::Vector{T})
	value::T = zero(T)
	dimension = size(field.points, 1)
	u = Array{T}(dimension)
	for i = 1:size(field.points, 2)
		for j = 1:dimension
			u[j] = x[j] - field.points[j, i]
		end
		n = norm(u)
		if n > 0
			value += field.alphas[i] * phi(field.phi, n) + dphi(field.phi, n) / n * (field.betas[:,i]' * u)[1]
		end
	end
	value
end

function grad{T}(field::HermiteRadialField{T}, x::Vector)
	dimension = size(field.points, 1)
	num_points = size(field.points, 2)
	g = zeros(T, dimension)
	for i = 1:num_points
		u = x - field.points[:,i]
		n = norm(u)

		if n > 1e-5
			uhat = u ./ n
			df = dphi(field.phi, n)
			ddf = ddphi(field.phi, n)
			alpha_df = field.alphas[i] * df
			beta_uhat = (field.betas[:,i]' * uhat)[1]

			g += alpha_df .* uhat + beta_uhat * (ddf * uhat - u * df / n^2) + field.betas[:,i] * df / n
		end
	end
	g
end

function grad{T}(field::HermiteRadialField{T})
	FunctionalVectorField(x -> grad(field, x))
end
