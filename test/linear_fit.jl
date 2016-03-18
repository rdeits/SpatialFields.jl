using SpatialFields
using Base.Test
import MultiPoly
using GeometryTypes

function test_linear_fit()
    x, y = MultiPoly.generators(MultiPoly.MPoly{Float64}, :x, :y)

    true_polynomial = 3 + 4x - 7y
    @show true_polynomial

    points = Point{2, Float64}[[5; -2], [-1; 10], [-1; 11]]
    data = Float64[MultiPoly.evaluate(true_polynomial, p...) for p in points]
    @show data

    fit_polynomial = SpatialFields.linear_fit(points, data)
    @show fit_polynomial

    for p in points
        @test isapprox(MultiPoly.evaluate(fit_polynomial, p...), MultiPoly.evaluate(true_polynomial, p...))
    end
    @test isapprox(fit_polynomial[0,0], true_polynomial[0,0])
    @test isapprox(fit_polynomial[1,0], true_polynomial[1,0])
    @test isapprox(fit_polynomial[0,1], true_polynomial[0,1])
end

test_linear_fit()

function test_vector_linear_fit()
    x, y = MultiPoly.generators(MultiPoly.MPoly{Float64}, :x, :y)

    true_polynomials = [3 + 4x - 7y; -1 - 8x + 2y]
    points = Point{2, Float64}[[5; -2], [-1;10], [-1; 11]]
    data = Point{2, Float64}[[MultiPoly.evaluate(poly, p...) for poly in true_polynomials] for p in points]

    fit_polynomials = SpatialFields.linear_fit(points, data)

    for i = 1:length(true_polynomials)
        for p in points
            @test isapprox(MultiPoly.evaluate(fit_polynomials[i], p...), MultiPoly.evaluate(true_polynomials[i], p...))
        end
        @test isapprox(fit_polynomials[i][0,0], true_polynomials[i][0,0])
        @test isapprox(fit_polynomials[i][1,0], true_polynomials[i][1,0])
        @test isapprox(fit_polynomials[i][0,1], true_polynomials[i][0,1])
    end
end

test_vector_linear_fit()
