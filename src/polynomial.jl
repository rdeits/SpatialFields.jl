type PolynomialScalarField{T} <: ScalarField
	polynomial::MultiPoly.MPoly{T}
end
convert{T}(::Type{ScalarField}, poly::MultiPoly.MPoly{T}) = PolynomialScalarField(poly)

type PolynomialVectorField{T} <: VectorField
	partials::Vector{MultiPoly.MPoly{T}}
end
convert{T}(::Type{VectorField}, partials::Vector{MultiPoly.MPoly{T}}) = PolynomialVectorField(partials)

function evaluate{T}(polynomial::MultiPoly.MPoly{T}, x::Vector{T})
	@assert length(x) == length(polynomial.vars)
	result::T = zero(T)
	for (powers, coeff) in polynomial.terms
		term::T = zero(T)
		for i = 1:length(powers)
			term += x[i] ^ powers[i]
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

