using SpatialFields
using BenchmarkTools
import StaticArrays: SVector

function evaluation_speed()
    println("Evaluation speed, hrbf:")
    points = [rand(SVector{3, Float64}) for i in 1:100]
    normals = [rand(SVector{3, Float64}) for i in 1:100];
    field = HermiteRadialField(points, normals)
    println(@benchmark ($field)(x) setup=(x=zeros(SVector{3, Float64})))

    println("Evaluation speed, interpolating:")
    edge_points = [rand(SVector{3, Float64}) for i in 1:50]
    interior_points = [rand(SVector{3, Float64}) for i in 1:50]
    points = vcat(edge_points, interior_points)
    values = vcat([0.0 for i in edge_points], [-1.0 for i in interior_points])
    surface = InterpolatingSurface(points, values)
    println(@benchmark ($surface)(x) setup=(x=zeros(SVector{3, Float64})))
end

evaluation_speed()
