type InterpolatingSurface{N, T, PhiType} <: ScalarField{N, T}
    weights::Vector{T}
    points::Vector{Point{N, T}}
    offset::MultiPoly.MPoly{T}
    phi::PhiType
end

function InterpolatingSurface{Dimension, T}(points::Vector{Point{Dimension, T}}, values, phi_function=XCubed())
    @assert length(points) == length(values)
    @assert Dimension <= 3

    num_points = length(points)
    A = zeros(T, num_points + Dimension + 1, num_points + Dimension + 1)
    b = vcat(values, zeros(Dimension + 1))

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
                    f = phi(phi_function, n)
                    A[i,j] = f
                end
            end
        end
    end

    y = A \ b

    if any(isnan, y)
        @show A
        @show b
        @show y
        @show values
        error("NaN after linear solve")
    end
    weights = y[1:num_points]
    p = y[num_points+1:end]

    terms = OrderedDict(zeros(Int, Dimension) => p[1])
    for i = 1:Dimension
        powers = zeros(Int, Dimension)
        powers[i] = 1
        terms[powers] = p[i+1]
    end
    offset = MultiPoly.MPoly{T}(terms, [:x, :y, :z][1:Dimension])
    InterpolatingSurface(weights, points, offset, phi_function)
end

function evaluate{Dimension, FieldType, InputType}(surface::InterpolatingSurface{Dimension, FieldType}, x::Point{Dimension, InputType})
    result = zero(promote_type(FieldType, InputType))
    num_points = length(surface.points)
    for i = 1:num_points
        n = norm(x - surface.points[i])
        if n >= 1e-9
            result += surface.weights[i] * phi(surface.phi, n)
        end
    end
    
    result += evaluate(surface.offset, x)
end

evaluate{Dimension, T}(field::InterpolatingSurface{Dimension}, x::Vector{T}) = evaluate(field, convert(Point{Dimension, T}, x))


    


