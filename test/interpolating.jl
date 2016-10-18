using SpatialFields
using Base.Test
using StaticArrays

@testset "2d" begin
    edge_points = SVector{2, Float64}[[0; 0], [1; 0], [1; 1], [0; 1]]
    interior_points = SVector{2, Float64}[[0.5; 0.5]]

    points = vcat(edge_points, interior_points)
    values = vcat([0.0 for i in edge_points], [-1.0 for i in interior_points])
    surface = InterpolatingSurface(points, values)

    for point in edge_points
        @test isapprox(surface(point), 0.0, atol=1e-6)
    end

    for point in interior_points
        @test isapprox(surface(point), -1.0, atol=1e-6)
    end
end

@testset "3d" begin
    edge_points = SVector{3, Float64}[]
    for z = 0:1
        for x = 0:1
            for y = 0:1
                push!(edge_points, [x; y; z])
            end
        end
    end
    interior_points = SVector{3, Float64}[[0.5; 0.5; 0.5]]

    points = vcat(edge_points, interior_points)
    values = vcat([0.0 for i in edge_points], [-1.0 for i in interior_points])
    surface = InterpolatingSurface(points, values)

    for point in edge_points
        @test isapprox(surface(point), 0.0, atol=1e-6)
    end

    for point in interior_points
        @test isapprox(surface(point), -1.0, atol=1e-6)
    end
end

@testset "normalized" begin
    edge_points = SVector{2, Float64}[[0; 0], [1; 0], [1; 1], [0; 1]]
    interior_points = SVector{2, Float64}[[0.5; 0.5]]

    points = vcat(edge_points, interior_points)
    values = vcat([0.0 for i in edge_points], [-1.0 for i in interior_points])
    surface = InterpolatingSurface(points, values, SpatialFields.XCubed(), true)

    @show SpatialFields.gradient(surface)([10, 0.5])
    @test isapprox(SpatialFields.gradient(surface)([10, 0.5]), [1, 0], atol=1e-3)
    @test isapprox(SpatialFields.gradient(surface)([-10, 0.5]), [-1, 0], atol=1e-3)
    @test isapprox(SpatialFields.gradient(surface)([0.5, 10]), [0, 1], atol=1e-3)
    @test isapprox(SpatialFields.gradient(surface)([0.5, -10]), [0, -1], atol=1e-3)
end
