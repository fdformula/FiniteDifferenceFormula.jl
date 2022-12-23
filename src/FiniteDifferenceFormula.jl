module FiniteDifferenceFormula

#
# This Julia code shows how to derive n-point finite difference formulas for
# the 1st, 2nd, ..., derivatives by using Taylor series expansions of a function
# at evenly spaced points.
#
# David Wang, dwang at liberty dot edu, on 12/20/2022
#

export computecoefs, formula

max_num_of_taylor_terms = 30 # number of the 1st terms of a Taylor series expansion

using RowEchelon             # provides rref used in the code

mutable struct FDData
    n
    points
    k
    m
    coefs
end

data = FDData                # share computing results between functions
computedq::Bool = false      # make sure computecoefs(n, points) is called first

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
    return [1 // convert(BigInt,factorial(big(n))) * h^n for n in 0 : max_num_of_taylor_terms - 1]
end

# convert a coefficient to a readable string
function c2s(c, first_termq = false)
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
            print(" h")
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
    # setup a linear system Ax = B
    if n < 1
        @error "Wrong order of derivatives."
    end

    # setup the coefficients of Taylor series expansions of f(x) at each of the involved points
    len = points.stop - points.start + 1
    coefs = Array{Any}(undef, len)
    for i in eachindex(points)
        coefs[i] = taylor_coefs(points[i])
    end

    # The linear combination of f(x[i+start]), f(x[i+start+1]) ... f(x[i+stop])
    #
    # k[1]*f(x[i+start]) + k[2]*f(x[i+start+1]) + k[3]*f(x[i+start+2]) + k[stop-start+1]*f(x[i+stop])i
    #  = m*f^(n)(x[i]) + ..., m != 0
    #
    # shall eliminate f(x[i]), f'(x[i]), ..., f^(n-1)(x[i]); if needed, it shall eliminate also
    # f^(n+1)(x[i]), ..., on the right-hand side of all involved Taylor series expansions.
    #
    # For example, to eliminate f(x[i]) on the righ-hand side, we have
    #
    #    k[1]*coefs[1][1] + k[2]*coefs[2][1] + ... + k[stop-start+1]*coefs[stop-start+1][1] = 0
    #
    # and to eliminate f'(x[i]), we have
    #
    #    k[1]*coefs[1][2] + k[2]*coefs[2][2] + ... + k[stop-start+1]*coefs[stop-start+1][2] = 0
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
    A[len, 1] = 1; # let k[1] = 1
    A[len, 2 : len] .= 0

    B = Array{Rational{BigInt}}(undef, len, 1)
    fill!(B, 0)
    B[len] = 1       # so that k[1] = 1

    k = rref([A B])[:, len + 1]   # package RowEchelon provides rref
    k = k // gcd(k)
    m::Any = 0;
    order_index = n + 1
    for i in eachindex(points)
        m += k[i] * coefs[i][order_index]
    end

    if m < 0; k *= -1; m *= -1; end
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

    if printformulaq; formula(); end

    return (results, m)
end   # computecoefs

# print the linear combination
# k[1]*f(x[i+start]) + k[2]*f(x[i+start+1]) + ... + k[stop-start+1]*f(x[i+stop])
function print_lcombination(data::FDData)
    for i = eachindex(data.points)
        times = abs(data.k[i]) == 1 ? "" : "* "
        if i == 1
            print(c2s(data.k[i], true), times, f2s(data.points[i]))   # assume k[1] != 0
        else
            if data.k[i] != 0
                print(c2s(data.k[i]), times, f2s(data.points[i]))
            end
        end
    end
end

function print_formula(data::FDData, bigO="")
    if bigO == ""    # it is to print Julia formula
        th = data.n == 1 ? "st" : (data.n == 2 ? "nd" : "th")
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

        print("f_$(n)_points_$(s)_$(data.n)$(th)_derivative(f, x, i, h)")
    else
        if data.n <= 3
            print("f" * "'"^(data.n))
        else
            print("f^($(data.n))")
        end
        print("(x[i])")
    end
    print(" = ( ")
    print_lcombination(data)
    print(" ) / ")

    if data.m > 1; print("($(data.m) * "); end
    print("h")
    if data.n > 1; print("^$(data.n)"); end
    if data.m > 1; print(")"); end
    if bigO != ""; print(" + $bigO"); end
    print("\n\n")
end

# print readable formula and other computing results
# use data stored in global variable data
function formula()
    global data, computedq

    if !computedq
        @error "Please call computecoefs(n, points) first!"
    end

    # print Taylor series expansion of the linear combination:
    # k[1]*f(x[i+start]) + k[2]*f(x[i+start+1]) + ... + k[stop-start+1]*f(x[i+stop])
    println("****** This shows the computing result ******")
    print_lcombination(data)

    # Taylor series expansion of the linear combination
    sum_coefs = data.k[1] * data.coefs[1]
    for i = eachindex(data.points)
        if i > 1 && data.k[i] != 0
            sum_coefs += data.k[i] * data.coefs[i]
        end
    end

    print(" =\n")
    print_taylor(sum_coefs, 5); # print at most the number of nonzero terms

    # invalid input because no sufficient points are provided?
    valid_inputq = true
    for i in 1 : data.n
        if sum_coefs[i] != 0
            valid_inputq = false;
            break;
        end
    end
    if valid_inputq # why not verify validity of input in the beginning? for teaching's purpose
        println("****** This is the exact formula ******")
        # print the very formula with big-O notation of the trucation error
        #
        # find x of O(h^x)
        x = max_num_of_taylor_terms;
        for i in data.n + 2 : length(sum_coefs)
            if sum_coefs[i] != 0; x = i; break; end
        end
        x -= data.n + 1
        bigO = "O(h"
        if x > 1; bigO *= "^$x"; end
        bigO *= ")"

        print_formula(data, bigO)

        # print the formula in another format
        if data.m > 1
            data1 = FDData(data.n, data.points, data.k // data.m, 1, data.coefs)
            print("Or\n\n")
            print_formula(data1, bigO)
        end
        print("Julia function:\n\n")
        print_formula(data)
    else
        th = data.n == 1 ? "st" : (data.n == 2 ? "nd" : "th")
        print("!!!!! The input $(data.points) is invalid because at least $(data.n + 1) points are required for the $(data.n)$th derivative.\n\n")
    end
end  # formula

end # module
