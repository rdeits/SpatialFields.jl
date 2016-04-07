using SpatialFields
using Base.Test

function test_interpolating_2d()
    edge_points = Point{2, Float64}[[0; 0], [1; 0], [1; 1], [0; 1]]
    interior_points = Point{2, Float64}[[0.5; 0.5]]

    points = vcat(edge_points, interior_points)
    values = vcat([0.0 for i in edge_points], [-1.0 for i in interior_points])
    surface = InterpolatingSurface(points, values)

    for point in edge_points
        @test isapprox(evaluate(surface, point), 0.0, atol=1e-6)
    end

    for point in interior_points
        @test isapprox(evaluate(surface, point), -1.0, atol=1e-6)
    end

    g = grad(surface)
    p = Point(1.0, 0.0)
    @test isapprox(grad(surface, p), evaluate(g, p))
end

test_interpolating_2d()
