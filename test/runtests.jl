using FiniteDifferenceFormula
using Test

@testset "FiniteDifferenceFormula.jl" begin
	@test sum(compute(1, 0:1)[1]) == 0
	@test sum(compute(3, -2:5)[1]) == 0
	@test compute(3, -2:5)[2] > 0
end
