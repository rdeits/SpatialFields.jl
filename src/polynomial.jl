type PolynomialScalarField{T} <: ScalarField
	polynomial::MultiPoly.MPoly{T}
end

type PolynomialVectorField{T} <: VectorField
	partials::Vector{MultiPoly.MPoly{T}}
end

function evaluate{T}(field::PolynomialScalarField{T}, x)
	@assert length(x) == length(field.polynomial.vars)
	MultiPoly.evaluate(field.polynomial, x...)
end

function grad{T}(field::PolynomialScalarField{T})
	@assert length(x) == length(field.polynomial.vars)

	diffs = [MultiPoly.diff(field.polynomial, v) for v in vars]
	PolynomialVectorField(diffs)
end

evaluate{T}(field::PolynomialVectorField{T}, x) = [MultiPoly.evaluate(p, x) for p in field.partials]

