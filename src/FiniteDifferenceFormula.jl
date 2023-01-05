module FiniteDifferenceFormula

#
# This Julia code shows how to derive n-point finite difference formulas for
# the 1st, 2nd, ..., derivatives by using Taylor series expansions of a function
# at evenly spaced points.
#
# David Wang, dwang at liberty dot edu, on 12/20/2022
#

export computecoefs, formula, decimalplaces, taylor

max_num_of_taylor_terms = 30 # number of the 1st terms of a Taylor series expansion
                             # variable, depending on input to computecoefs

decimal_places = 15          # use it to print Julia function for a formula
                             # call decimalplaces(n) to reset it

####################################################################################

using Printf
using RowEchelon             # provide rref

mutable struct FDData
    n
    points
    k
    m
    coefs
end

data = FDData                # share computing results between functions
computedq::Bool = false      # make sure computecoefs(n, points) is called first
formula_status = 0           # a formula may not be valid

# a vector of the coefficients of Taylor series expansion of the linear combination:
# k[1]*f(x[i+start]) + k[2]*f(x[i+start+1]) + ... + k[stop-start+1]*f(x[i+stop])
lcombination_coefs = Array{Any}(undef, max_num_of_taylor_terms)

# This function returns the first 'max_num_of_taylor_terms' of Taylor series of
# f(x[i+1]) centered at x=x[i] in a vector with f(x[i]), f'(x[i]), ..., removed. The
# purpose is to obtain Taylor series expansions for f(x[i±k]) = f(x[i]±kh]) which
# are used to derive the m-point finite difference formulas for the first, second,
# ..., order derivatives at x[i].
#
#       f(x[i+1]) = f(x[i]) + 1/1! f'(x[i]) h + 1/2! f''(x[i]) h^2 + ...
#
# where h = x[i+1] - x[i].
#
# Usage:
#   f(x[i+k]) - call taylor_coefs(k)
#   f(x[i-k]) - call taylor_coefs(-k)
#   where k = 0, 1, 2, 3, ...
function taylor_coefs(h)
    # the simplest implementation:
    #   return [1 // convert(BigInt,factorial(big(n))) * h^n for n in 0 : max_num_of_taylor_terms - 1]
    # but the following code shows better time performance
    result = Matrix{Rational{BigInt}}(undef, 1, max_num_of_taylor_terms)
    factorial::BigInt = 1
    for n in 1 : max_num_of_taylor_terms
        N = n - 1                        # order of a derivative in Taylor series
        if N > 0; factorial *= N; end    # 0! = 1
        result[n] = 1 // factorial * h^N
    end
    return result
end

function decimalplaces(n)
    global decimal_places
    if isinteger(n) && n >= 2
        decimal_places = n
    else
        error("decimalplaces(n): n must be integer greater than 1")
    end
end

# convert a coefficient to a readable string
function c2s(c, first_termq = false, floating = false)
    s = ""
    if c < 0
        s = first_termq ? "-" : s = " - "
    elseif !first_termq
        s = " + "
    end

    c = abs(c)
    if isinteger(c)
        if c != 1
            s *= string(round(BigInt, c), " ")
        end
    elseif floating
        global decimal_places
        fmt = Printf.Format("%.$(decimal_places)f ")
        s *= Printf.format(fmt, convert(Float64, c))
    else
        s *= string(numerator(c), "/", string(denominator(c)), " ");
    end

    return s
end  # c2s

# convert f(x[i+k]) to a readable string
function f2s(k)
    s = "f(x[i"
    if k != 0
        if k > 0; s *= "+"; end
        s = "$s$k"
    end
    return "$s])"
end  # f2s

# print readable Taylor series expansion of f(x[i + j]) about x[i] (for teaching!)
function taylor(j::Int, num_of_nonzero_terms = 10)
    global max_num_of_taylor_terms
    if max_num_of_taylor_terms < num_of_nonzero_terms
        max_num_of_taylor_terms = num_of_nonzero_terms
    end
    coefs = taylor_coefs(j)
    println("\nf(x[i" * (j == 0 ? "" : (j > 0 ? "+$j" : "$j")) * "]) =")
    print_taylor(coefs, num_of_nonzero_terms)
end

# print readable Taylor series
#
# Input: An array that contains the coefficients of the first
#        "max_num_of_taylor_terms" of Taylor series expansion of a function
function print_taylor(coefs, num_of_nonzero_terms = max_num_of_taylor_terms)
    first_termq = true
    for n in 0 : max_num_of_taylor_terms - 1
        N = n + 1
        if coefs[N] == 0; continue; end

        print(c2s(coefs[N], first_termq))
        if abs(coefs[N]) != 1; print("* "); end
        if n <= 3
            print("f$("'" ^ n)(x[i])")
        else
            print("f^($n)(x[i])")
        end
        if n >= 1
            print(" * h")
            if n > 1; print("^$n"); end
        end
        first_termq = false

        num_of_nonzero_terms -= 1
        if num_of_nonzero_terms == 0; break; end
    end
    println(" + ...\n")
end  # print_taylor

#
# Algorithm
# ---------
# General generator of finite difference formulas for the n-th order derivatives
# of f(x) at x[i] using points = start:stop.
#
# It uses the linear combination of f(x[i+j]), j in start:stop, to eliminate f(x[i]), f'(x[i]), ...,
# so that the first term of the Taylor series expansion of the linear combination is f^(n)(x[i]):
#
#    k[1]*f(x[i+start]) + k[2]*f(x[i+start+1]) + ... + k[stop-start+1]*f(x[i+stop])
#        = m*f^(n)(x[i]) + ..., m > 0, n > 0
#
# It is this equation that gives the formula for computing f^(n)(x[i]).
#
# Values of the cofficients k[:] and m will be determined, and the first few terms of the remainder
# (...) will be listed for estimating the truncation error of a formula.
#
# Julia's Rational type and related arithmetic operations, especially a // b, together with the
# RowEchelon package, are a perfect fit for obtaining an "exact" formula (in the sense of exact
# solutions).
#
# Input
# -----
#             n: the n-th order derivative to be found
#        points: in the format of a range, start:stop. (See examples below.)
# printformulaq: print the computed formula?
#
#    points  |  The points to be used
#    --------+----------------------------------------------
#     0:2    |  x[i], x[i+1], x[i+2]
#    -2:2    |  x[i-2], x[i-1], x[i], x[i+1], x[i+2]
#    -3:2    |  x[i-3], x[i-2], x[i-1], x[i], x[i+1], x[i+2]
#
# Output
# ------
# The function returns a tuple, ([k[1], k[2], ..., k[stop-start+1]], m).
#
function computecoefs(n::Int, points::UnitRange{Int}, printformulaq::Bool = false)
    # setup a linear system Ax = B first
    if n < 1; @error "Wrong order of derivatives."; end

    len = length(points)
    global max_num_of_taylor_terms = max(len, n) + 8
    global lcombination_coefs = Array{Any}(undef, max_num_of_taylor_terms)

    # setup the coefficients of Taylor series expansions of f(x) at each of the involved points
    coefs = Array{Any}(undef, len)
    for i in eachindex(points)
        coefs[i] = taylor_coefs(points[i])
    end

    # We find a linear combination of f(x[i+start]), f(x[i+start+1]) ... f(x[i+stop]),
    #
    # k[1]*f(x[i+start]) + k[2]*f(x[i+start+1]) + k[3]*f(x[i+start+2]) + k[stop-start+1]*f(x[i+stop])i
    #  = 0*f(x[i]) + 0*f'(x[i]) + ... + 0*f^(n-1)(x[i]) + m*f^(n)(x[i]) + ..., m != 0
    #
    # so that it must eliminate f(x[i]), f'(x[i]), ..., f^(n-1)(x[i]); if possible, it also
    # eliminates f^(n+1)(x[i]), f^(n+2)(x[i]), ....
    #
    # For example, to eliminate f(x[i]), we have
    #
    #    k[1]*coefs[1][1] + k[2]*coefs[2][1] + ... + k[stop-start+1]*coefs[stop-start+1][1] = 0 
    #
    # and to eliminate f'(x[i]), we have
    #
    #    k[1]*coefs[1][2] + k[2]*coefs[2][2] + ... + k[stop-start+1]*coefs[stop-start+1][2] = 0
    #
    # Therefore, a linear system is detemined by the following equations
    #
    #    k[1]*coefs[1][j] + k[2]*coefs[2][j] + ... + k[stop-start+1]*coefs[stop-start+1][j] = 0 ..... (1)
    #
    # where j = 1, 2, ..., n.
    #
    A = Matrix{Rational{BigInt}}(undef, len, len)
    row = 1
    for order = 0 : len - 1
        if order == n; continue; end

        # eliminating f^(order)(x[i])
        order_index = order + 1
        for i = eachindex(points)
            A[row, i] = coefs[i][order_index]
        end
        row += 1
    end

    # The homogeneous linear system (1) has no nontrivial solution or has inifinitely many nontrivial
    # solutions. It is understandable that it may not have a nontrivial solution. But why inifinitely
    # many solutions? Because if k[:] is a solution, λ k[:] is also a solution, where λ is any
    # nonzero real constant, i.e., nontrivial solutions are parallel to each other. So, in the case
    # that there are infinitely many nontrivial solutions, if we know one entry k[which] is nonzero
    # and let it be a nonzero constant (say, 1), then, the solution k[:] is unique.
    #
    # Purpose: determine 'which' so that k[which] != 0. It is not necessarily 1.
    #
    # Note: By commonsense, we assume that, except x[i], the closer a point is near x[i], the larger
    # its weight is as usual/expected, and the closest points to x[i] are x[i ± 1] if available.
    which = 1
    if -points.start == points.stop     # central
        which = points.stop             # x[i - 1]
    elseif points.start == 0            # forward
        which = 2                       # x[i + 1]
    elseif points.stop == 0             # backward
        which = len - 1                 # x[i - 1]
    else                                # mixed, e.g., -5:3, -4:-2, 2:5
        if points.start < 0 && points.stop > 0
            which = -points.start       # x[i - 1]
        elseif points.start < 0 && points.stop < 0
            which = len                 # the point closest to x[i]
        #else points.start > 0 && points.stop > 0
        #    which = 1                  # the point closest to x[i]
        end
    end

    A[len, 1 : len] .= 0;               # let k[which] = 1
    A[len, which] = 1

    B = zeros(Rational{BigInt}, len, 1)
    B[len] = 1                          # so that k[which] = 1

    # solve Ax = B for x, i.e., k[:]
    k = rref([A B])[:, len + 1]         # package RowEchelon provides rref
    k = k // gcd(k)

    # Taylor series expansion of the linear combination
    # k[1]*f(x[i+start]) + k[2]*f(x[i+start+1]) + ... + k[stop-start+1]*f(x[i+stop])
    lcombination_coefs = k[1] * coefs[1]
    for i = eachindex(points)
        if i == 1 || k[i] == 0; continue; end
        lcombination_coefs += k[i] * coefs[i]
    end

    # find the first nonzero term, v1.0.3
    m = lcombination_coefs[n + 1]
    for i = 1 : n
        if lcombination_coefs[i] != 0
            m = lcombination_coefs[i]
            break
        end
    end

    # "normalize" k[:] and m so that m is a positive integer
    if m < 0; k *= -1; lcombination_coefs *= -1; end
    m = lcombination_coefs[n + 1]
    x = round(BigInt, m)
    if x == m; m = x; end

    results = Matrix{Any}(undef, 1, len)
    for i in eachindex(points)
        x = round(BigInt, k[i])
        results[i] = x == k[i] ? x : k[i]
    end

    # save the results in a global variable for other functions
    global data, computedq
    data = FDData(n, points, k, m, coefs)
    computedq = true

    test_formula_validity()

    if printformulaq; formula(); end

    return (results, m)
end   # computecoefs

# print the linear combination
# k[1]*f(x[i+start]) + k[2]*f(x[i+start+1]) + ... + k[stop-start+1]*f(x[i+stop])
function print_lcombination(data::FDData, decimal = false)
    firstq = true
    for i = eachindex(data.points)
        if data.k[i] == 0; continue; end
        times = abs(data.k[i]) == 1 ? "" : "* "
        print(c2s(data.k[i], firstq, decimal), times, f2s(data.points[i]))
        firstq = false
    end
end

function test_formula_validity()
    # to find f^(n)(x[i]) by obtaining
    #
    #    k[1]*f(x[i+start]) + k[2]*f(x[i+start+1]) + ... + k[stop-start+1]*f(x[i+stop])
    #        = m*f^(n)(x[i]) + ..., m > 0
    #
    # the most important step is to know if f(x[i]), f'(x[i]), ..., f^(n-1)(x[i]) are all eliminated, i.e.,
    #    k[1]*coefs[1][j] + k[2]*coefs[2][j] + ... + k[stop-start+1]*coefs[stop-start+1][j] = 0
    # where j = 1:n
    global data, formula_status, computedq
    if !computedq
        @error "Please call computecoefs(n, points) first!"
    end

    n = data.n
    k = data.k
    coefs = data.coefs
    points = data.points
    len = length(points)

    formula_status = 0
    for i = 1 : n
        if lcombination_coefs[i] != 0
            println("-" ^ 81)
            println("***** Error:: n=$n, $(points.start):$(points.stop) :: i = $i, k[1]*coefs[1][$i] + k[2]*coefs[2][$i] + ... + k[$len]*coefs[$len][$i] != 0")
            if len <= n
                th = n == 1 ? "st" : (n == 2 ? "nd" : (n == 3 ? "rd" : "th"))
                println("The input $(points) is invalid because at least $(n + 1) points are needed for the $n$th derivative.\n")
            else
                println("A formula might not exist.\n")
            end
            return
        end
    end

    m = lcombination_coefs[n + 1]
    if m == 0
        println("-" ^ 81)
        println("***** Error:: n=$n, $(points.start):$(points.stop) :: m = 0, formula can't be found!")
        return
    end

    if sum(k) != 0   # sum of coefficients must be 0
        println("\n***** Warning:: n=$n, $(points.start):$(points.stop) :: sum(k[:]) != 0")
        formula_status += 1
    end

    # must coefficients of central formulas be symmetrical?
    if abs(points.start) == points.stop  # add && false to false to skip the test
        j = length(k)
        for i in 1 : round(Int64, length(k)/2)
            if abs(k[i]) != abs(k[j])
                println("\n***** Warning:: n=$n, $(points.start):$(points.stop) :: k[$i] != k[$j]")
                formula_status += 1
                break
            end
            j -= 1
        end
    end

    if formula_status == 0
        formula_status = 100    # perfect
    end
end  # test_formula_validity

function print_formula(data::FDData, bigO="", decimal = false)
    if bigO == ""    # printing Julia function
        th = data.n == 1 ? "st" : (data.n == 2 ? "nd" : (data.n == 3 ? "rd" : "th"))
        s = "mixed"
        if -data.points.start == data.points.stop
            s = "central"
        elseif data.points.start == 0
            s = "forward"
        elseif data.points.stop == 0
            s = "backward"
        end

        n = 0  # how many points are involved?
        for i in eachindex(data.points)
            if data.k[i] == 0; continue; end
            n += 1
        end

        print("f$(data.n)$(th)deriv$(n)pt$(s)(f, x, i, h)")
    else
        if data.n <= 3
            print("f" * "'"^(data.n))
        else
            print("f^($(data.n))")
        end
        print("(x[i])")
    end
    print(" = ( ")
    print_lcombination(data, decimal)
    print(" ) / ")

    if data.m > 1; print("($(data.m) * "); end
    print("h")
    if data.n > 1; print("^$(data.n)"); end
    if data.m > 1; print(")"); end
    if bigO != ""; print(" + $bigO"); end
    println("\n")
end  # print_formula

# print readable formula and other computing results
# use data stored in global variable data
function formula()
    global data, computedq, formula_status

    if !computedq
        @error "Please call computecoefs(n, points) first!"
    end

    if formula_status !=0
        println("-" ^ 81)
        print("The following formula ")
        if formula_status == 100
            print("passed all tests: sum of coefs being zero")
            if abs(data.points.start) == data.points.stop
                print(", symmetry of coefs about x[i]")
            end
            println(", etc.\n")
        else
            println("may still be valid, though it didn't pass tests like sum of coefs being zero.\n")
        end
    end

    # print Taylor series expansion of the linear combination:
    # k[1]*f(x[i+start]) + k[2]*f(x[i+start+1]) + ... + k[stop-start+1]*f(x[i+stop])
    println("The computing result:\n")
    print_lcombination(data)
    print(" =\n")
    print_taylor(lcombination_coefs, 5); # print at most 5 nonzero terms

    if formula_status != 0
        println("The exact formula:\n")
        # print the very formula with big-O notation of the trucation error
        #
        # find x of O(h^x)
        x = max_num_of_taylor_terms;
        for i in data.n + 2 : length(lcombination_coefs)
            if lcombination_coefs[i] != 0; x = i; break; end
        end
        x -= data.n + 1
        bigO = "O(h"
        if x > 1; bigO *= "^$x"; end
        bigO *= ")"

        print_formula(data, bigO)
        # print the formula in another format
        data1 = FDData
        if data.m > 1
            data1 = FDData(data.n, data.points, data.k // data.m, 1, data.coefs)
            print("Or\n\n")
            print_formula(data1, bigO)
        end

        print("Julia function:\n\n")
        print_formula(data)                # exact
        # print the function in decimal format
        if data.m > 1
            print("Or\n\n")
            print_formula(data1)           # exact
            print("Or\n\n")
            print_formula(data1, "", true) # decimal
        end
    end
end  # formula

end # module
