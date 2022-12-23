using FiniteDifferenceFormula
using Test

@testset "FiniteDifferenceFormula.jl" begin
	# confirmed
	computecoefs(1, 0:1, true)
	computecoefs(1, 0:2, true)
	computecoefs(1, -1:0, true)
	computecoefs(1, -2:0, true)
	computecoefs(1, -1:1, true)
	computecoefs(1, -2:2, true)
	computecoefs(2, 0:2, true)
	computecoefs(2, 0:3, true)
	computecoefs(2, -2:0, true)
	computecoefs(2, -3:0, true)
	computecoefs(2, -1:1, true)
	computecoefs(2, -2:2, true)
	computecoefs(3, 0:3, true)
	computecoefs(3, 0:4, true)
	computecoefs(3, -3:0, true)
	computecoefs(3, -4:0, true)
	computecoefs(3, -2:2, true)
	computecoefs(3, -3:3, true)
	computecoefs(4, 0:4, true)
	computecoefs(4, 0:5, true)
	computecoefs(4, -4:0, true)
	computecoefs(4, -5:0, true)
	computecoefs(4, -2:2, true)
	computecoefs(4, -3:3, true)
	computecoefs(2, -3:3, true)
	computecoefs(2, -4:4, true)

	# for fun or exploring
	computecoefs(1, -2:5, true)
	computecoefs(2, -4:1)
	formula()

	# Doesn't work bc f^(4)(x[i]) term survives? More points are needed!
	computecoefs(5, -3:1, true)
	computecoefs(5, -2:2, true)
	computecoefs(5, 0:4, true)

	computecoefs(5, 0:5, true)
	computecoefs(1, 2:5, true)
	computecoefs(1, -5:-2, true)
	computecoefs(1, -6:6, true)
	computecoefs(2, -5:5, true)
	computecoefs(3, -6:6, true)
	computecoefs(8, -5:5, true)
	computecoefs(10, -5:5)
	formula()
	computecoefs(2, -3:8, true)
end
