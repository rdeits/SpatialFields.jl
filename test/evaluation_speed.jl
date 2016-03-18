using SpatialFields
using ProfileView

function evaluation_speed()
    points = SpatialFields.Point{3, Float64}[rand(3) for i in 1:100]
    normals = SpatialFields.Normal{3, Float64}[rand(3) for i in 1:100];
    x = SpatialFields.Point{3, Float64}(0)
    field = SpatialFields.HermiteRadialField(points, normals);

    println("Evaluation speed:")
    Profile.clear()
    for i = 1:1e4; SpatialFields.evaluate(field, x); end
    @time @profile for i = 1:1e4; SpatialFields.evaluate(field, x); end
    ProfileView.view()
end

evaluation_speed()
