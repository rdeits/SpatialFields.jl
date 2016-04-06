abstract BaseTwiceDifferentiableFunction

type TwiceDifferentiableFunction <: BaseTwiceDifferentiableFunction
	f
	df
	ddf
end
phi(func::TwiceDifferentiableFunction, x) = func.f(x)
dphi(func::TwiceDifferentiableFunction, x) = func.df(x)
ddphi(func::TwiceDifferentiableFunction, x) = func.ddf(x)

function TwiceDifferentiableFunction(f::Function)
	TwiceDifferentiableFunction(f,
		x -> ForwardDiff.derivative(f, x),
		x -> ForwardDiff.derivative(y -> ForwardDiff.derivative(f, y), x))
end

immutable XCubed <: BaseTwiceDifferentiableFunction
end

phi(::XCubed, x) = x^3
dphi(::XCubed, x) = 3x^2
ddphi(::XCubed, x) = 6x
