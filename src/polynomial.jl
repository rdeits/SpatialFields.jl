type PolynomialScalarField{T} <: ScalarField
	polynomial::MultiPoly.MPoly{T}
end
convert{T}(::Type{ScalarField}, poly::MultiPoly.MPoly{T}) = PolynomialScalarField(poly)

type PolynomialVectorField{T} <: VectorField
	partials::Vector{MultiPoly.MPoly{T}}
end
convert{T}(::Type{VectorField}, partials::Vector{MultiPoly.MPoly{T}}) = PolynomialVectorField(partials)

function evaluate{T}(polynomial::MultiPoly.MPoly{T}, x)
	@assert length(x) == length(polynomial.vars)
	result::T = zero(T)
	for (powers, coeff) in polynomial.terms
		term::T = one(T)
		for i = 1:length(powers)
			term *= x[i] ^ powers[i]
		end
		result += term * coeff
	end
	result
end

function evaluate{T}(field::PolynomialScalarField{T}, x)
	evaluate(field.polynomial, x)
end

function grad{T}(field::PolynomialScalarField{T})
	diffs = [MultiPoly.diff(field.polynomial, v) for v in field.polynomial.vars]
	PolynomialVectorField(diffs)
end

evaluate{T}(field::PolynomialVectorField{T}, x) = [evaluate(p, x) for p in field.partials]

function convert{N, T}(::Type{Array{T, 2}}, points::Vector{Point{N, T}})
	A = Array{T}(N, length(points))
	for i = 1:N
		for j = 1:length(points)
			A[i,j] = points[j][i]
		end
	end
	A
end

function linear_fit{T}(coordinates::Array{T, 2}, data::Vector{T})
	dimension = size(coordinates, 1)
	num_points = size(coordinates, 2)
	A = hcat(coordinates', ones(num_points))
	v = A \ data
	coeffs = OrderedDict(zeros(dimension) => v[end])
	for j = 1:dimension
		powers = zeros(Int, dimension)
		powers[j] = 1
		coeffs[powers] = v[j]
	end
	return MultiPoly.MPoly{T}(coeffs, [:x, :y, :z][1:dimension])
end

function linear_fit{N, T}(coordinates::Vector{Point{N, T}}, data::Vector{T})
	linear_fit(convert(Array{T, 2}, coordinates), data)
end

function linear_fit{N, T}(coordinates::Vector{Point{N, T}}, data::Vector{Point{N, T}})
	A = convert(Array{T, 2}, coordinates)
	MultiPoly.MPoly{T}[linear_fit(A, [d[i] for d in data]) for i in 1:N]
end
