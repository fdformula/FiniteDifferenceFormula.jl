module FiniteDifferenceFormula

#
# This Julia module generates n-point finite difference formulas for the 1st, 2nd, ...,
# derivatives by using Taylor series expansions of a function at evenly spaced points.
# It also allows users to examine and test formulas right away in Julia REPL. It surely
# helps when we teach/learn numerical computing, especially, the finite difference
# method.
#
# David Wang, dwang at liberty dot edu, on 12/20/2022
#

export computecoefs, formula, truncationerror
export decimalplaces, activatejuliafunction, taylor

_max_num_of_taylor_terms = 30 # number of the 1st terms of a Taylor series expansion
                              # variable, depending on input to computecoefs

_decimal_places = 16          # use it to print Julia function for a formula
                              # call decimalplaces(n) to reset it

####################################################################################

using Printf
using RowEchelon              # provide rref

mutable struct _FDData
    n
    points
    k
    m
    coefs
end

_data                        = _FDData    # share computing results between functions
_computedq::Bool             = false      # make sure computecoefs(n, points) is called first
_formula_status              = 0          # a formula may not be available
                                          # values? see function _test_formula_validity()

# a vector of the coefficients of Taylor series expansion of the linear combination:
# k[1]*f(x[i+start]) + k[2]*f(x[i+start+1]) + ... + k[stop-start+1]*f(x[i+stop])
_lcombination_coefs          = Array{Any}(undef, _max_num_of_taylor_terms)

_range_inputq                = false
_range_input::UnitRange{Int} = 0:0        # if computecoefs receives a range, save it

_julia_exact_func_expr       = ""         # string of 1st exact Julia function for f^(n)(x[i])
_julia_exact_func_expr1      = ""         # string of 2nd exact Julia function for f^(n)(x[i])
_julia_decimal_func_expr     = ""         # string of decimal Julia function for f^(n)(x[i])
_julia_func_basename         = ""

_bigO                        = ""         # truncation error of a finidifference formula

# This function returns the first '_max_num_of_taylor_terms' of Taylor series of
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
#   f(x[i+k]) - call _taylor_coefs(k)
#   f(x[i-k]) - call _taylor_coefs(-k)
#   where k = 0, 1, 2, 3, ...
function _taylor_coefs(h)
    # the simplest implementation:
    #   return [1 // convert(BigInt,factorial(big(n))) * h^n for n in 0 : _max_num_of_taylor_terms - 1]
    # but the following code shows better time performance
    result = Matrix{Rational{BigInt}}(undef, 1, _max_num_of_taylor_terms)
    factorial::BigInt = 1
    for n in 1 : _max_num_of_taylor_terms
        N = n - 1                        # order of a derivative in Taylor series
        if N > 0; factorial *= N; end    # 0! = 1
        result[n] = 1 // factorial * h^N
    end
    return result
end  # _taylor_coefs

function decimalplaces()
    global _decimal_places
    return _decimal_places
end  # decimalplaces

function decimalplaces(n)
    global _decimal_places, _computedq
    if isinteger(n) && n > 0
        _decimal_places = n
        if _computedq
            println("Please call 'formula' to generate (or 'activatejuliafunction' to generate and activate) a Julia function for the newly computed formula, using the new decimal places.")
        else
            println("You may start your work by calling 'computecoefs'.")
        end
    else
        println("decimalplaces(n): n must be a positive integer.")
    end

    return
end  # decimalplaces

# convert a coefficient to a readable string
function _c2s(c, first_termq = false, floatingq = false)
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
    elseif floatingq
        global _decimal_places
        fmt = Printf.Format("%.$(_decimal_places)f ")
        s *= Printf.format(fmt, convert(Float64, c))
    else
        s *= string(numerator(c), "/", string(denominator(c)), " ");
    end

    return s
end  # _c2s

# convert f(x[i+k]) to a readable string
function _f2s(k)
    s = "f(x[i"
    if k != 0
        if k > 0; s *= "+"; end
        s = "$s$k"
    end
    return "$s])"
end  # _f2s

# print readable Taylor series expansion of f(x[i + j]) about x[i] (for teaching!)
function taylor(j::Int, num_of_nonzero_terms = 10)
    global _max_num_of_taylor_terms
    if _max_num_of_taylor_terms < num_of_nonzero_terms
        _max_num_of_taylor_terms = num_of_nonzero_terms
    end
    coefs = _taylor_coefs(j)
    println("\nf(x[i" * (j == 0 ? "" : (j > 0 ? "+$j" : "$j")) * "]) =")
    _print_taylor(coefs, num_of_nonzero_terms)
    return
end  # taylor

# print readable Taylor series
#
# Input: An array that contains the coefficients of the first
#        "_max_num_of_taylor_terms" of Taylor series expansion of a function
function _print_taylor(coefs, num_of_nonzero_terms = _max_num_of_taylor_terms)
    first_termq = true
    for n in 0 : _max_num_of_taylor_terms - 1
        N = n + 1
        if coefs[N] == 0; continue; end

        print(_c2s(coefs[N], first_termq))
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
    return
end  # _print_taylor

function computecoefs(n::Int, points::UnitRange{Int}, printformulaq::Bool = false)
    if n < 1
        println("Wrong order of derivatives :: $n. It must be a positive integer.")
        return
    end

    global _range_inputq = true
    global _range_input  = points
    return _computecoefs(n, collect(points), printformulaq)
end

# computecoefs(2, [1 2 3 -1])
function computecoefs(n::Int, points::Matrix{Int}, printformulaq::Bool = false)
    if n < 1
        println("Wrong order of derivatives :: $n. It must be a positive integer.")
        return
    end

    if length(points) <= 1
        println("Invalid input :: $points - more points are needed.")
        return
    end

    m, = size(points)
    if m > 1 points = points'; end      # a column vector
    points = sort(unique(points))
    if length(points) == 1
         println("Invalid input :: points = $(points')")
         return
    end

    println("You input: $points")
    # is the input 'points' actually a range?
    if length(points[1]:points[end]) == length(points)
        return computecoefs(n, points[1]:points[end], printformulaq)
    else
        global _range_inputq = false
        global _range_input  = 0:0
        return _computecoefs(n, points, printformulaq)
    end
end  # computecoefs

# computecoefs(2, [1, 2, 3, -1])
function computecoefs(n::Int, points::Array{Int}, printformulaq::Bool = false)
    return computecoefs(n, hcat(points), printformulaq)
end

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
function _computecoefs(n::Int, points::Vector{Int}, printformulaq::Bool = false)
    global _computedq = false
    global _formula_status = 0
    global _bigO = ""

    # for teaching's purpose, we don't do so
    # if length(points) <= n
    #     pts = _range_inputq ? "$(_range_input)" : "$(points')"
    #     th = n == 1 ? "st" : (n == 2 ? "nd" : (n == 3 ? "rd" : "th"))
    #     println("$pts is invalid b/c at least $(n + 1) points are needed for the $n$th derivative.")
    #     return
    # end

    # setup a linear system Ax = B first
    len = length(points)
    global _max_num_of_taylor_terms = max(len, n) + 8
    global _lcombination_coefs = Array{Any}(undef, _max_num_of_taylor_terms)

    # setup the coefficients of Taylor series expansions of f(x) at each of the involved points
    coefs = Array{Any}(undef, len)
    for i in 1 : len
        coefs[i] = _taylor_coefs(points[i])
    end

    # We find a linear combination of f(x[i+start]), f(x[i+start+1]) ... f(x[i+stop]),
    #
    # k[1]*f(x[i+start]) + k[2]*f(x[i+start+1]) + k[3]*f(x[i+start+2]) + k[stop-start+1]*f(x[i+stop])i
    #  = 0*f(x[i]) + 0*f'(x[i]) + ... + 0*f^(n-1)(x[i]) + m*f^(n)(x[i]) + ..., m != 0
    #
    # so that it must eliminate f(x[i]), f'(x[i]), ..., f^(n-1)(x[i]); given more points are provided,
    # it also eliminates f^(n+1)(x[i]), f^(n+2)(x[i]), ....
    #
    # For example, to eliminate f(x[i]), we have
    #
    #    k[1]*coefs[1][1] + k[2]*coefs[2][1] + ... + k[len]*coefs[len][1] = 0
    #
    # and to eliminate f'(x[i]), we have
    #
    #    k[1]*coefs[1][2] + k[2]*coefs[2][2] + ... + k[len]*coefs[len][2] = 0
    #
    # Therefore, a linear system is detemined by the following equations
    #
    #    k[1]*coefs[1][j] + k[2]*coefs[2][j] + ... + k[len]*coefs[len][j] = 0 ..... (1)
    #
    # where j = 1, 2, ..., len but not n.
    #
    A = Matrix{Rational{BigInt}}(undef, len, len)
    row = 1
    for order in 0 : len - 1             # corresponding to f(x[i]), f'(x[i]), f''(x[i]), ...
        if order == n; continue; end     # skip f^(n)(x[i])

        # eliminating f^(order)(x[i])
        # A[:,j] stores coefs of k[1]*coefs[1][j] + k[2]*coefs[2][j] + ... + k[len]*coefs[len][j] = 0
        for j = 1 : len
            A[row, j] = coefs[j][order + 1]
        end
        row += 1
    end

    # The homogeneous linear system (1) has no nontrivial solution or has inifinitely many nontrivial
    # solutions. It is understandable that it may not have a nontrivial solution. But why inifinitely
    # many nontrivial solutions? It is because, if k[:] is a nontrivial solution, α k[:] is also a
    # nontrivial solution, where α is any nonzero real constant, i.e., all nontrivial solutions (a
    # subspace spanned by this k[:]) are parallel to each other. Therefore, in the case that there
    # are infinitely many nontrivial solutions, if we know one entry k[which] is nonzero and let it
    # be a nonzero constant (say, 1), then, a nontrivial solution k[:] is unique/found.
    #
    # Beware, there may be multiple nontrivial solutions,
    #    k[:] = [k1, k2, ...], k[:] = [K1, K2, ...], ..., or k[:] = [κ1, κ2, ...]
    # of which no two are parallel to each other. However, each of these solutions must satisfy the
    # condition that both k[1] and k[end] can't be zero. Why? If, say, k[1] = 0, in other words, a
    # formula only uses/depends on (at most) x[i+start+1], x[i+start+2], ..., x[i+stop], why should
    # we say it is a formula that uses/depends on x[i+start] (and x[i+start+1], ..., x[i+stop])?
    # Therefore, we can assume that k[1] != 0.
    #
    A[len, 2 : len] .= 0                  # setup an equation so that k[1] = 1
    A[len, 1]        = 1

    B = zeros(Rational{BigInt}, len, 1)
    B[len] = 1                            # so that k[1] = 1

    # solve Ax = B for x, i.e., k[:]
    k = rref([A B])[:, len + 1]           # package RowEchelon provides rref
    k = k // gcd(k)

    # Taylor series expansion of the linear combination
    # k[1]*f(x[i+start]) + k[2]*f(x[i+start+1]) + ... + k[stop-start+1]*f(x[i+stop])
    _lcombination_coefs = k[1] * coefs[1] # let Julia determine the type
    for i in 2 : len
        if k[i] == 0; continue; end
        _lcombination_coefs += k[i] * coefs[i]
    end

    # find the first nonzero term, v1.0.3
    m = _lcombination_coefs[n + 1]
    for i in 1 : n
        if _lcombination_coefs[i] != 0
            m = _lcombination_coefs[i]
            break
        end
    end

    # "normalize" k[:] and m so that m is a positive integer
    if m < 0; k *= -1; _lcombination_coefs *= -1; end
    m = _lcombination_coefs[n + 1]
    x = round(BigInt, m)
    if x == m; m = x; end  # so that m = 5 rather than 5//1

    # save the results in a global variable for other functions
    global _data = _FDData(n, points, k, m, coefs)
    _computedq = true

    _test_formula_validity()

    if printformulaq; formula(); end

    return (coefs, m)
end   # _computecoefs

# return a string of the linear combination
# k[1]*f(x[i+start]) + k[2]*f(x[i+start+1]) + ... + k[stop-start+1]*f(x[i+stop])
#
# To test in Julia if a formula is valid, we must convert the data to BigInt or BigBloat
# becasue, say, for computecoefs(2, -12:12), without the conversion, the formula will
# give obviously wrong answers because of overflow and rounding errors. julia_REPL_funcq
# is here for testing in Julia REPL only. For ordinary definition of a function, call
# formula().
#
function _lcombination_expr(data::_FDData, decimalq = false, julia_REPL_funcq = false)
    firstq = true
    s = ""
    for i = eachindex(data.points)
        if data.k[i] == 0; continue; end
        times = abs(data.k[i]) == 1 ? "" : "* "
        if julia_REPL_funcq && firstq
            c2s = _c2s(data.k[i], true, decimalq)
            if c2s == ""
                c2s = "1"
            elseif c2s == "-"
                c2s = "-1"
            end
            s *= "big(" * c2s * ")"     # convert one value to BigInt/BigFloat
        else
            s *= _c2s(data.k[i], firstq, decimalq)
        end
        s *= times * _f2s(data.points[i])
        firstq = false
    end
    return s
end  # _lcombination_expr

# check if the newly computed formula is valid. results are saved in the global variable
# _formula_status:
#    100 - Perfect, even satifiying some commonly seen "rules", such as sum of coefficients = 0,
#          symmetry of cofficients about x[i] in a central formula
#   -100 - No formula can't be found
#   > 0  - Mathematically, the formula is still valid but may not satisfy some commonly seen "rules"
#          such as sum of coefficients = 0 and symmetry of coefficients about x[i]
#
function _test_formula_validity()
    # to find f^(n)(x[i]) by obtaining
    #
    #    k[1]*f(x[i+start]) + k[2]*f(x[i+start+1]) + ... + k[stop-start+1]*f(x[i+stop])
    #        = m*f^(n)(x[i]) + ..., m > 0
    #
    # the most important step is to know if f(x[i]), f'(x[i]), ..., f^(n-1)(x[i]) are all eliminated, i.e.,
    #    k[1]*coefs[1][j] + k[2]*coefs[2][j] + ... + k[stop-start+1]*coefs[stop-start+1][j] = 0
    # where j = 1:n
    global _data, _range_inputq, _range_input

    n = _data.n
    k = _data.k
    coefs = _data.coefs
    points = _data.points
    len = length(points)

    input_points = _range_inputq ? _range_input : points'

    # Is there any equation in system (1) that is not satisfied?
    has_solutionq = true
    global _formula_status = 0
    for i = 1 : n
        if _lcombination_coefs[i] != 0
            println("***** Error:: $n, $input_points :: i = $i, k[1]*coefs[1][$i] + k[2]*coefs[2][$i] + ... + k[$len]*coefs[$len][$i] != 0")
            has_solutionq = false
            break
        end
    end

    # Is m = 0 ?
    if has_solutionq
        m = _lcombination_coefs[n + 1]
        if m == 0
            println("-" ^ 81)
            println("\n***** Error:: $n, $input_points :: m = 0, formula can't be found.")
            has_solutionq = false
        end
    end

    # there could be a solution that doesn't use the whole range of input (points)
    if has_solutionq
        start, stop = 1, len    # actual left and right endpoints
        while start <= len && k[start] == 0; start += 1; end
        while stop >= 1 && k[stop] == 0; stop -= 1; end
        if start >= stop
            has_solutionq = false
        elseif start > 1 || stop < len
            s = _range_inputq ? (points[start] : points[stop]) : points[start : stop]
            println("Error:: $n, $input_points :: A formula can't be found.\n")
            println("Warning:: $n, $s might be your input for which a formula is found.\n")
        end
    end

    if !has_solutionq
        if len <= n
            th = n == 1 ? "st" : (n == 2 ? "nd" : (n == 3 ? "rd" : "th"))
            println("Error:: $n, $input_points :: Invalid input because at least $(n + 1) points are needed for the $n$th derivative.\n")
        else
            println("Error:: $n, $input_points :: A formula can't be found.")
        end
        return
    end

    if sum(k) != 0   # sum of coefficients must be 0
        println("\n***** Warning:: $n, $input_points :: sum(k[:]) != 0")
        _formula_status += 1
    end

    # must coefficients of central formulas be symmetrical?
    if _range_inputq && abs(_range_input.start) == _range_input.stop  # add && false to skip the test
        j = length(k)
        for i in 1 : round(Int64, length(k)/2)
            if abs(k[i]) != abs(k[j])
                println("\n***** Warning:: $n, $input_points :: k[$i] != k[$j]")
                _formula_status += 1
                break
            end
            j -= 1
        end
    end

    if _formula_status == 0
        _formula_status = 100    # perfect
    end

    return
end  # _test_formula_validity

function _denominator_expr(data::_FDData)
    s  = data.m > 1 ? "($(data.m) * " : ""
    s *= "h"
    if data.n > 1; s *= "^$(data.n)"; end
    if data.m > 1; s *= ")"; end
    return s
end

# print and return the function for the newly computed finite difference formula
function _julia_func_expr(data::_FDData, decimalq = false, julia_REPL_funcq = false)
    global _range_inputq, _range_input

    s = ""
    if _range_inputq
        if -_range_input.start == _range_input.stop
            s = "central"
        elseif _range_input.start == 0
            s = "forward"
        elseif _range_input.stop == 0
            s = "backward"
        end
    end

    n = 0    # how many points are actually involved?
    for i in eachindex(data.points)
        if data.k[i] == 0; continue; end
        n += 1
    end

    th = data.n == 1 ? "st" : (data.n == 2 ? "nd" : (data.n == 3 ? "rd" : "th"))
    global _julia_func_basename = "f$(data.n)$(th)deriv$(n)pt$(s)"
    fexpr  = "(f, x, i, h) = "
    if julia_REPL_funcq; fexpr *= "Float64("; end
    fexpr *= "( "
    fexpr *= _lcombination_expr(data, decimalq, julia_REPL_funcq)
    fexpr *= " ) / " * _denominator_expr(data)
    if julia_REPL_funcq; fexpr *= ")"; end

    return fexpr
end  # _julia_func_expr

# print the formula with big-O notation for the newly computed formula
function _print_bigo_formula(data::_FDData, bigO)
    print(data.n <= 3 ? "f" * "'"^(data.n) : "f^($(data.n))")
    print("(x[i]) = ( ", _lcombination_expr(data, false))
    print(" ) / ", _denominator_expr(data), " + $bigO\n\n")
end  # _print_bigo_formula

# create and print readable formula and other computing results
# using data stored in global variable _data
#
# No valid formula can be found? Still dump the computing result for teaching.
function formula()
    global _data, _computedq, _formula_status, _range_inputq, _range_input, _bigO

    if !_computedq
        println("Please call computecoefs(n, points) first!")
        return
    end

    if _formula_status > 0
        println("-" ^ 80)
        print("The following formula ")
        if _formula_status == 100
            print("passed all tests: sum of coefs being zero")
            if _range_inputq && abs(_range_input.start) == _range_input.stop
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
    print(_lcombination_expr(_data), " =\n")
    _print_taylor(_lcombination_coefs, 5);       # print at most 5 nonzero terms

    if _formula_status > 0
        println("The exact formula:\n")
        # find x in the big-O notation O(h^x)
        x = _max_num_of_taylor_terms
        for i in _data.n + 2 : length(_lcombination_coefs)
            if _lcombination_coefs[i] != 0; x = i; break; end
        end
        x -= _data.n + 1
        _bigO = "O(h"
        if x > 1; _bigO *= "^$x"; end
        _bigO *= ")"

        _print_bigo_formula(_data, _bigO)
        data1 = _FDData
        if _data.m > 1           # print in another format
            data1 = _FDData(_data.n, _data.points, _data.k // _data.m, 1, _data.coefs)
            print("Or\n\n")
            _print_bigo_formula(data1, _bigO)
        end

        println("Julia function:\n")
        global _julia_exact_func_expr = _julia_func_expr(_data)
        print(_julia_func_basename, "e", _julia_exact_func_expr, "\n\n")
        if _data.m > 1           # other formats
            global _julia_exact_func_expr1  = _julia_func_expr(data1)
            print("Or\n\n", _julia_func_basename, "e1", _julia_exact_func_expr1, "\n\nOr\n\n")
            global _julia_decimal_func_expr = _julia_func_expr(data1, true)
            print(_julia_func_basename, "d", _julia_decimal_func_expr, "\n\n")
        end  # m = 1? No decimal formula
    else
        _formula_status = -100 # no formula can't be generated
    end

    return
end  # formula

# print _bigO, the estimation of trucation error in the big-O notation
function truncationerror()
    global _bigO, _computedq, _formula_status
    if _computedq
        if _bigO == ""; formula(); end

        if _formula_status == -100
            println("No valid formula is available. Please call 'computecoefs' again.")
        else
            println(_bigO)
        end
    else
        println("Please call 'computecoefs' first.")
    end
    return
end

# activate function(s) for the newly computed finite difference formula, allowing
# immediate evaluation of the formula in present Julia REPL session.
function activatejuliafunction()
    global _computedq, _formula_status, _data

    if !_computedq
        println("Please run 'computecoefs' first.")
        return
    end

    global _julia_exact_func_expr  = ""
    global _julia_exact_func_expr1  = ""
    global _julia_decimal_func_expr = ""

    # generate Julia function for current REPL session
    if _formula_status > 0
        global _julia_exact_func_expr = _julia_func_expr(_data, false, true)
        if _data.m > 1           # print in other formats
            data1 = _FDData(_data.n, _data.points, _data.k // _data.m, 1, _data.coefs)
            global _julia_exact_func_expr1  = _julia_func_expr(data1, false, true)
            global _julia_decimal_func_expr = _julia_func_expr(data1, true, true)
        end  # m = 1? No decimal formula
    else
        println("No valid formula is for activation. Please run 'computecoefs' with different input again.")
        return
    end

    eval(Meta.parse(_julia_func_basename * "e" * _julia_exact_func_expr))
    if _julia_decimal_func_expr != ""
        eval(Meta.parse(_julia_func_basename * "e1" * _julia_exact_func_expr1))
        eval(Meta.parse(_julia_func_basename * "d" * _julia_decimal_func_expr))
    end

    println("The following function(s) are available temporiarily in the FiniteDifferenceFormula module. Usage:\n")
    println("  import FiniteDifferenceFormula as fd\n")
    println("  f, x, i, h = sin, 0:0.01:10, 501, 0.01")
    eval(Meta.parse("f, x, i, h = sin, 0:0.01:10, 501, 0.01"))

    exact = sin(_data.n * pi /2 + 5)
    apprx = eval(Meta.parse("$(_julia_func_basename)e(f, x, i, h)"))
    println("  fd.$(_julia_func_basename)e(f, x, i, h)    # output: ", apprx)
    if _julia_decimal_func_expr != ""
        apprx = eval(Meta.parse("$(_julia_func_basename)e1(f, x, i, h)"))
        println("  fd.$(_julia_func_basename)e1(f, x, i, h)   # output: ", apprx)

        approx = eval(Meta.parse("$(_julia_func_basename)d(f, x, i, h)"))
        println("  fd.$(_julia_func_basename)d(f, x, i, h)    # output: ", apprx)
    end
    len = length("fd.$(_julia_func_basename)")
    println(" "^(len + 19) * "# cp:     ", exact)
    println("\nCall fd.formula() to view the very definition.")

    return  # stop Julia from returning something users never expect
end  # activatejuliafunction

end # module
