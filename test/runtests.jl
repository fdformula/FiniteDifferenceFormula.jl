using FiniteDifferenceFormula
using Test

@testset "FiniteDifferenceFormula.jl" begin
	# confirmed
	@test sum(computecoefs(1, 0:1)[1]) == 0
	@test sum(computecoefs(1, 0:2)[1]) == 0
	@test sum(computecoefs(3, -2:5)[1]) == 0
	@test computecoefs(3, -2:5)[2] > 0
end
