VERSION >= v"0.4" && __precompile__()

module SpatialFields

using GeometryTypes
import ForwardDiff
import MultiPoly
import Base: convert

export Point,
	grad,
	evaluate,
	ScalarField,
	VectorField,
	HermiteRadialField,
	FunctionalVectorField,
	PolynomialScalarField,
	PolynomialVectorField

abstract ScalarField{N, T}
abstract VectorField{N, T}

type FunctionalVectorField{N, T} <: VectorField{N, T}
	f::Function
end
evaluate(field::FunctionalVectorField, x) = field.f(x)

# Default gradient implementation for any scalar field, using ForwardDiff for
# automatic differentiation. This is likely to be slower than a custom
# gradient implementation for a particular type, but it's a useful fallback to
# have.
function grad{N, T}(field::ScalarField{N, T})
	FunctionalVectorField{N, T}(ForwardDiff.gradient(x -> evaluate(field, x)))
end

# Shortcut for evaluating the gradient without returning a vector field
function grad(field::ScalarField, x)
	evaluate(grad(field), x)
end

include("hrbf.jl")
include("polynomial.jl")

end # module
