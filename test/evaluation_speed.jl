using SpatialFields
# using ProfileView

function evaluation_speed()
    points = Point{3, Float64}[rand(3) for i in 1:100]
    normals = Normal{3, Float64}[rand(3) for i in 1:100];
    x = Point{3, Float64}(0)
    field = HermiteRadialField(points, normals);

    println("Evaluation speed:")
    # Profile.clear()
    for i = 1:1e4; evaluate(field, x); end
    @time for i = 1:1e4; evaluate(field, x); end
    # ProfileView.view()
end

evaluation_speed()
