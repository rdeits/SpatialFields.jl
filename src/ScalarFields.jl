module ScalarFields

import ForwardDiff
import MultiPoly

export grad, 
	evaluate,
	HermiteRadialField,
	ScalarField,
	VectorField

abstract ScalarField
abstract VectorField

type FunctionalVectorField <: VectorField
	f::Function
end
evaluate(field::FunctionalVectorField, x) = field.f(x)

# Default gradient implementation for any scalar field, using ForwardDiff for
# automatic differentiation. This is likely to be slower than a custom
# gradient implementation for a particular type, but it's a useful fallback to
# have.
function grad(field::ScalarField)
	FunctionalVectorField(ForwardDiff.gradient(x -> evaluate(field, x)))
end

# Shortcut for evaluating the gradient without returning a vector field
function grad(field::ScalarField, x)
	evaluate(grad(field), x)
end

include("hrbf.jl")

end # module
