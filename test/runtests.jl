using SpatialFields
using Base.Test

@testset "interpolating" begin
    include("interpolating.jl")
end

@testset "hrbf_2d" begin
    include("hrbf_2d.jl")
end

@testset "hrbf_3d" begin
    include("hrbf_3d.jl")
end

@testset "hrbf_2d_5th_power" begin
    include("hrbf_2d_5th_power.jl")
end

@testset "polynomial" begin
    include("polynomial.jl")
end

@testset "linear fit" begin
    include("linear_fit.jl")
end

@testset "benchmarks" begin
    include("../perf/runbenchmarks.jl")
end
