type HermiteRadialField{N, T, F <: Function} <: AbstractScalarField{N}
	alphas::Vector{T}
	betas::Vector{SVector{N, T}}
	points::Vector{SVector{N, T}}
	phi::F
end

function HermiteRadialField{Dimension, T, F <: Function}(points::Vector{SVector{Dimension, T}}, normals::Vector{SVector{Dimension, T}}, phi::F=XCubed())
	@assert length(points) == length(normals)

	num_points = length(points)

	A = Array{T}(num_points * (Dimension + 1), num_points * (Dimension + 1))
	b = Array{T}(num_points * (Dimension + 1))
	# u = Array{T}(Dimension)
	# v = Array{T}(Dimension)
    dphi = derivative(phi)
    ddphi = derivative(dphi)

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
				f = phi(n)
                df = dphi(n)
                ddf = ddphi(n)
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
	betas = [SVector{Dimension, T}(y[2:end, i]) for i in 1:num_points]
	HermiteRadialField{Dimension, T, F}(alphas, betas, points, phi)
end


function (field::HermiteRadialField{Dimension, T}){Dimension, T}(x::AbstractVector{T})
	value::T = zero(T)
	num_points = length(field.points)
	@assert length(x) == Dimension
	@assert length(field.alphas) == num_points
	@assert length(field.betas) == length(field.points)
    dphi = derivative(field.phi)
	@inbounds for i = eachindex(field.points)
		u = x - field.points[i]
		n = norm(u)
		if n > 0
			z = sum(field.betas[i] .* u)
			p = field.phi(n)
			dp = dphi(n)
			value += field.alphas[i] * p
			value += dp / n * z
		end
	end
	value
end

function gradient{Dimension, T, PhiType}(field::HermiteRadialField{Dimension, T, PhiType}, x::AbstractVector{T})
	num_points = length(field.points)
	@assert length(x) == Dimension
	@assert length(field.alphas) == num_points
	@assert length(field.betas) == num_points
	g = zeros(SVector{Dimension, T})
    dphi = derivative(field.phi)
    ddphi = derivative(dphi)
	for i = eachindex(field.points)
		u = x - field.points[i]
		n = norm(u)

		if n > 1e-5
			uhat = u ./ n
            df = dphi(n)
            ddf = ddphi(n)
			alpha_df = field.alphas[i] * df
			beta_uhat = sum(field.betas[i] .* uhat)

			g += alpha_df .* uhat + beta_uhat * (ddf * uhat - u * df / n^2) + field.betas[i] * df / n
		end
	end
	g
end

gradient(field::HermiteRadialField) = x -> gradient(field, x)

# function bounds(field::HermiteRadialField)
# 	lb = Vec(minimum(field.points))
# 	ub = Vec(maximum(field.points))
# 	HyperRectangle(lb, ub - lb)
# end
