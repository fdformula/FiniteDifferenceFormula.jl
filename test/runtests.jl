using FiniteDifferenceFormula
using Test

@testset "FiniteDifferenceFormula.jl" begin
	@test sum(compute(1, 0:1)[3]) == 0
	@test sum(compute(3, -2:5)[3]) == 0
	@test compute(3, -2:5)[4] > 0
end
