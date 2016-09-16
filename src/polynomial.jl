function evaluate{T1, T2}(polynomial::MultiPoly.MPoly{T1}, x::AbstractVector{T2})
    @assert length(x) == length(polynomial.vars)
    result = zero(promote_type(T1, T2))
    for (powers, coeff) in polynomial.terms
        term = one(T2)
        for i = 1:length(powers)
            if powers[i] != 0
                term *= x[i] ^ powers[i]
            end
        end
        result += term * coeff
    end
    result
end

type PolynomialScalarField{N, T} <: AbstractScalarField{N}
	polynomial::MultiPoly.MPoly{T}
end
PolynomialScalarField{T}(poly::MultiPoly.MPoly{T}) = PolynomialScalarField{length(poly.vars), T}(poly)
(field::PolynomialScalarField)(x::AbstractVector) = evaluate(field.polynomial, x)

type PolynomialVectorField{N, T} <: AbstractVectorField{N}
	partials::SVector{N, MultiPoly.MPoly{T}}
end
(field::PolynomialVectorField)(x::AbstractVector) = map(p -> evaluate(p, x), field.partials)

gradient{N, T}(field::PolynomialScalarField{N, T}) = PolynomialVectorField{N, T}(map(v -> MultiPoly.diff(field.polynomial, v), field.polynomial.vars))

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

destructure{T,N}(A::Array{T,N}) = reinterpret(eltype(T), A, (size(T)..., size(A)...))

function linear_fit{N, T}(coordinates::Vector{SVector{N, T}}, data::Vector{T})
    linear_fit(destructure(coordinates), data)
end

function linear_fit{N, T}(coordinates::Vector{SVector{N, T}}, data::Vector{SVector{N, T}})
    A = destructure(coordinates)
    MultiPoly.MPoly{T}[linear_fit(A, [d[i] for d in data]) for i in 1:N]
end
