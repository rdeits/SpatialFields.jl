VERSION >= v"0.4" && __precompile__()

module SpatialFields

using GeometryTypes
import ForwardDiff
import MultiPoly
import Base: convert
import DataStructures: OrderedDict

export Point,
	Normal,
	grad,
	evaluate,
	bounds,
	ScalarField,
	VectorField,
	HermiteRadialField,
	InterpolatingSurface,
	FunctionalVectorField,
	PolynomialScalarField,
	PolynomialVectorField

abstract ScalarField{N, T}
abstract VectorField{N, T}


type FunctionalVectorField{N, T} <: VectorField{N, T}
	f::Function
end
evaluate{N, T}(field::FunctionalVectorField{N, T}, x) = field.f(x)

# Default gradient implementation for any scalar field, using ForwardDiff for
# automatic differentiation. This is likely to be slower than a custom
# gradient implementation for a particular type, but it's a useful fallback to
# have.
type AutoDiffVectorField{N, T} <: VectorField{N, T}
	grad_function
end
evaluate{N, T}(field::AutoDiffVectorField{N, T}, x::Vector{T}) = field.grad_function(x)
evaluate{N, T}(field::AutoDiffVectorField{N, T}, x) = evaluate(field, convert(Vector{T}, x))
function grad{N, T}(field::ScalarField{N, T})
	AutoDiffVectorField{N, T}(ForwardDiff.gradient(x -> evaluate(field, x)))
end

# Shortcut for evaluating the gradient without returning a vector field
function grad(field::ScalarField, x)
	evaluate(grad(field), x)
end

include("radial_functions.jl")
include("hrbf.jl")
include("interpolating.jl")
include("polynomial.jl")

end # module
