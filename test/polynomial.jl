using SpatialFields
using Base.Test
import MultiPoly

function test_polynomial()
  x, y = MultiPoly.generators(MultiPoly.MPoly{Float64}, :x, :y)
  poly = x^2 + 3x*y + y + y^2 + 10

  field = PolynomialScalarField(poly)

  @test isapprox(field([3., 7.]), MultiPoly.evaluate(poly, 3., 7.))

  vec_field = SpatialFields.gradient(field)
  x = [100, 200]
  @test isapprox(vec_field(x),
                 Float64[MultiPoly.evaluate(MultiPoly.diff(poly, v),
                                     x[1], x[2]) for v in [:x, :y]])
end

test_polynomial()
