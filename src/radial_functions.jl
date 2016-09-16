immutable XCubed <: Function end
@inline (::XCubed)(x) = x^3
derivative(::XCubed) = dXCubed()

immutable dXCubed <: Function end
@inline (::dXCubed)(x) = 3x^2
derivative(::dXCubed) = ddXCubed()

immutable ddXCubed <: Function end
@inline (::ddXCubed)(x) = 6x

immutable XSquaredLogX <: Function end
@inline (::XSquaredLogX)(x) = x^2 * log(x)
derivative(::XSquaredLogX) = dXSquaredLogX()

immutable dXSquaredLogX <: Function end
@inline (::dXSquaredLogX)(x) = 2x*log(x) + x
derivative(::dXSquaredLogX) = ddXSquaredLogX()

immutable ddXSquaredLogX <: Function end
@inline (::ddXSquaredLogX)(x) = 2*log(x) + 3
