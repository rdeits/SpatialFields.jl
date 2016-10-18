type InterpolatingSurface{N, T, F <: Function} <: AbstractScalarField{N}
    weights::Vector{T}
    points::Vector{SVector{N, T}}
    offset::MultiPoly.MPoly{T}
    phi::F
end

function InterpolatingSurface{Dimension, T, Tv}(points::Vector{SVector{Dimension, T}}, values::AbstractVector{Tv},
    phi::Function=XCubed(),
    normalize=false)
    @assert length(points) == length(values)
    @assert Dimension <= 3

    num_points = length(points)
    A = zeros(T, num_points + Dimension + 1, num_points + Dimension + 1)
    b = vcat(values, zeros(Tv, Dimension + 1))

    for i = 1:num_points
        A[i, end-Dimension] = 1.0
        for j = 1:Dimension
            A[i, end-Dimension+j] = points[i][j]
        end
        A[num_points+1, i] = 1.0
        for j = 1:Dimension
            A[end-Dimension+j,i] = points[i][j]
        end
        for j = 1:num_points
            if i != j
                v = points[i] - points[j]
                n = norm(v)
                if n >= 1e-9
                    f = phi(n)
                    A[i,j] = f
                end
            end
        end
    end

    y = A \ b

    weights = y[1:num_points]
    p = y[num_points+1:end]
    if normalize
        # Normalize the surface's SDF values to try to ensure that they
        # represent real distances
        averaging_factor = (1 + (Dimension - 1) / 2) / Dimension
        average_slope = averaging_factor * 3 * sum(weights .* (p -> sum(p .^ 2)).(points))
        weights ./= average_slope
        p ./= average_slope
    end

    terms = OrderedDict(zeros(Int, Dimension) => p[1])
    for i = 1:Dimension
        powers = zeros(Int, Dimension)
        powers[i] = 1
        terms[powers] = p[i+1]
    end
    offset = MultiPoly.MPoly{T}(terms, [:x, :y, :z][1:Dimension])

    InterpolatingSurface(weights, points, offset, phi)
end

Base.norm{T <: Real}(x::SVector{2, T}) = sqrt(x[1]^2 + x[2]^2)
Base.norm{T <: Real}(x::SVector{3, T}) = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
Base.norm{T <: Real}(x::MVector{2, T}) = sqrt(x[1]^2 + x[2]^2)
Base.norm{T <: Real}(x::MVector{3, T}) = sqrt(x[1]^2 + x[2]^2 + x[3]^2)

function (surface::InterpolatingSurface{N, Ts}){N, Ts, Tx}(x::AbstractVector{Tx})::promote_type(Ts, Tx)
    typealias T promote_type(Ts, Tx)
    result::T = zero(T)
    for i = eachindex(surface.points)
        n = norm(x - surface.points[i])
        if n >= 1e-9
            result += surface.weights[i] * surface.phi(n)
        end
    end
    result + evaluate(surface.offset, x)
end
