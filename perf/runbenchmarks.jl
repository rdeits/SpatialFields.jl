using SpatialFields
using BenchmarkTools
import StaticArrays: SVector

suite = BenchmarkGroup()


let
    suite["hrbf"] = BenchmarkGroup(["hrbf"])
    points = [rand(SVector{3, Float64}) for i in 1:100]
    normals = [rand(SVector{3, Float64}) for i in 1:100];
    suite["hrbf"]["construct"] = @benchmarkable field = HermiteRadialField($points, $normals)

    field = HermiteRadialField(points, normals)
    x = zeros(SVector{3, Float64})
    suite["hrbf"]["evaluate"] = @benchmarkable ($field)($x)
end

let
    suite["interpolating"] = BenchmarkGroup(["interpolating"])
    edge_points = [rand(SVector{3, Float64}) for i in 1:50]
    interior_points = [rand(SVector{3, Float64}) for i in 1:50]
    points = vcat(edge_points, interior_points)
    values = vcat([0.0 for i in edge_points], [-1.0 for i in interior_points])

    suite["interpolating"]["construct"] = @benchmarkable surface = InterpolatingSurface($points, $values)
    surface = InterpolatingSurface(points, values)
    x = zeros(SVector{3, Float64})
    suite["interpolating"]["evaluate_svector"] = @benchmarkable ($surface)($x)

    x = zeros(3)
    suite["interpolating"]["evaluate_vector"] = @benchmarkable ($surface)($x)
end

tune!(suite)
results = run(suite, verbose=false)
showall(results)
println()
