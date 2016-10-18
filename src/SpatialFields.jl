__precompile__()

module SpatialFields

import ForwardDiff
import MultiPoly
import Base: convert
import DataStructures: OrderedDict
import StaticArrays: SVector, MVector

export InterpolatingSurface,
    HermiteRadialField,
    PolynomialScalarField,
    PolynomialVectorField

derivative(f::Function) = x -> ForwardDiff.derivative(f, x)
gradient(f::Function) = x -> ForwardDiff.gradient(f, x)

abstract AbstractScalarField{N} <: Function
abstract AbstractVectorField{N} <: Function

type ScalarField{N, F <: Function} <: AbstractScalarField{N}
    func::F
end
ScalarField{F <: Function}(N::Integer, f::F) = ScalarField{N, F}(f)
(field::ScalarField{N}){N}(x::AbstractVector) = field.func(x)

type VectorField{N, F <: Function} <: AbstractVectorField{N}
    func::F
end
VectorField{F <: Function}(N::Integer, f::F) = VectorField{N, F}(f)
(field::VectorField{N}){N}(x::AbstractVector) = field.func(x)

gradient{N}(field::AbstractScalarField{N}) = VectorField(N, x -> ForwardDiff.gradient(field, x))

include("radial_functions.jl")
include("hrbf.jl")
include("interpolating.jl")
include("polynomial.jl")

end # module
