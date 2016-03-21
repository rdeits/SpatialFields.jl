using SpatialFields
using Base.Test

import Base: isapprox

isapprox{T1, D1, N1, T2, D2, N2}(v1::FixedSizeArrays.FixedArray{T1, D1, N1}, v2::FixedSizeArrays.FixedArray{T2, D2, N2}, args...; kwargs...) = begin
    isapprox(convert(Array{T1, D1}, v1), convert(Array{T2, D2}, v2), args...; kwargs...)
end

isapprox{T, D, N}(v1::FixedSizeArrays.FixedArray{T, D, N}, args...; kwargs...) = begin
    isapprox(convert(Array{T, D}, v1), args...; kwargs...)
end

isapprox{T, D, N}(v1, v2::FixedSizeArrays.FixedArray{T, D, N}, args...; kwargs...) = begin
    isapprox(v1, convert(Array{T, D}, v2), args...; kwargs...)
end

include("hrbf_2d.jl")
include("hrbf_3d.jl")
include("hrbf_2d_5th_power.jl")
include("polynomial.jl")
include("linear_fit.jl")
include("evaluation_speed.jl")
