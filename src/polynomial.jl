type PolynomialScalarField{T} <: ScalarField
	polynomial::MultiPoly.MPoly{T}
end
convert{T}(::Type{ScalarField}, poly::MultiPoly.MPoly{T}) = PolynomialScalarField(poly)

type PolynomialVectorField{T} <: VectorField
	partials::Vector{MultiPoly.MPoly{T}}
end
convert{T}(::Type{VectorField}, partials::Vector{MultiPoly.MPoly{T}}) = PolynomialVectorField(partials)

function evaluate{T}(field::PolynomialScalarField{T}, x)
	@assert length(x) == length(field.polynomial.vars)
	MultiPoly.evaluate(field.polynomial, x...)
end

function grad{T}(field::PolynomialScalarField{T})
	diffs = [MultiPoly.diff(field.polynomial, v) for v in field.polynomial.vars]
	PolynomialVectorField(diffs)
end

evaluate{T}(field::PolynomialVectorField{T}, x) = [MultiPoly.evaluate(p, x...) for p in field.partials]

