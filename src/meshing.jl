function convert{T <: HomogenousMesh}(::Type{T}, field::ScalarField{3})
    surface_mesh(field, bounds(field), 0.0)
end

function surface_mesh(field::ScalarField{3}, bounds::HyperRectangle, iso_level=0.0, sampling_rate=20)
    w = widths(bounds)

    lb = minimum(bounds) - w / 2.0
    window = HyperRectangle(Vec(lb), 2 * w)
    sdf = SignedDistanceField(x -> evaluate(field, x) - iso_level, window, maximum(w) / sampling_rate)
    mesh = HomogenousMesh(sdf, 0.0)

    rescaled_points = Point{3,Float64}[Vec(v-1) ./ (Vec(size(sdf))-1) .* (2 * w) + lb for v in vertices(mesh)]
    return HomogenousMesh(rescaled_points, mesh.faces)
end
