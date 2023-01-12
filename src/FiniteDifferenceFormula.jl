module FiniteDifferenceFormula

#
# This Julia module generates n-point finite difference formulas for the 1st,
# 2nd, ..., derivatives by using Taylor series expansions of a function at
# evenly spaced points. It also allows users to examine and test formulas right
# away in the current Julia REPL. It surely helps when we teach/learn numerical
# computing, especially, the finite difference method.
#
# David Wang, dwang at liberty dot edu, on 12/20/2022
#

export computecoefs, formula, truncationerror

# for teaching/learning/exploring
export decimalplaces, activatejuliafunction, taylor, printtaylor
export _set_default_max_num_of_taylor_terms

_default_max_num_of_taylor_terms = 30
# it should be a constant. usually, users never need to change it.

_decimal_places = 16          # use it to print Julia function for a formula
                              # call decimalplaces(n) to reset it

###############################################################################

# of the 1st terms of a Taylor series expansion
_max_num_of_taylor_terms =  _default_max_num_of_taylor_terms
                              # variable, depending on input to computecoefs

using Printf

mutable struct _FDData
    n
    points
    k
    m
    coefs
    max_num_of_taylor_terms
end

_data                        = _FDData    # share results between functions
_computedq::Bool             = false      # make sure computecoefs is called first
_formula_status              = 0          # a formula may not be available
                                          # values? see _test_formula_validity()

# a vector of the coefficients of Taylor series expansion of the linear combination:
# k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ... + k[len]*f(x[i+points[len]])
_lcombination_coefs          = Array{Any}

_range_inputq                = false
_range_input::UnitRange{Int} = 0:0        # computecoefs receives a range? save it

_julia_exact_func_expr       = ""         # 1st exact Julia function for f^(n)(x[i])
_julia_exact_func_expr1      = ""         # 2nd exact Julia function for f^(n)(x[i])
_julia_decimal_func_expr     = ""         # decimal Julia function for f^(n)(x[i])
_julia_func_basename         = ""

_bigO                        = ""         # truncation error of a formula

# This function returns the first '_max_num_of_taylor_terms' of Taylor series
# of f(x[i+1]) centered at x=x[i] in a vector with f(x[i]), f'(x[i]), ...,
# removed. The purpose is to obtain Taylor series expansions for f(x[i±k]) =
# f(x[i]±kh]) which are used to derive the m-point finite difference formulas
# for the first, second, ..., order derivatives at x[i].
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
    global _max_num_of_taylor_terms
    result = Matrix{Rational{BigInt}}(undef, 1, _max_num_of_taylor_terms)
    factorial::BigInt = 1
    for n in 1 : _max_num_of_taylor_terms
        N = n - 1          # order of derivatives in Taylor series
        if N > 0; factorial *= N; end
        result[n] = 1 // factorial * h^N
    end
    return result
end  # _taylor_coefs

# convert a coefficient to a readable string
function _c2s(c, first_termq = false, decimalq = false)
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
    elseif decimalq
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

# print readable Taylor series
#
# Input: An array that contains the coefficients of the first
#        "_max_num_of_taylor_terms" of Taylor series expansion of a function
function _print_taylor(coefs, num_of_nonzero_terms = 0) # 0? see the 1st statement
    global _max_num_of_taylor_terms
    if num_of_nonzero_terms == 0 # for translating to other programming language
        num_of_nonzero_terms = _max_num_of_taylor_terms
    end

    first_termq = true
    for n in 0 : length(coefs) - 1
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

# before v1.0.7, we thankfully used the function 'rref' provided by the package
# RowEchelon v0.2.1. since it has been removed from the base, it is safer for
# this package to have its own implementaion of the function. anyway, it is
# just a few lines of code.
function _rref(A::Matrix{Rational{BigInt}})  # change A to reduced row echelon form
    row, col = size(A)
    for i in 1 : min(row, col)
        # choose the largest entry in A[i:end, i] as the pivot
        lv, li = abs(A[i ,i]), i  # the largest value and its index
        for j = i + 1 : row
            absv = abs(A[j, i])
            if absv > lv; lv, li = absv, j; end
        end
        if lv == 0; continue; end
        if li != i                # interchange two rows
            A[i, i : end], A[li, i : end] = A[li, i : end], A[i, i : end]
        end

        A[i, i : end] /= A[i, i]  # on the main diagonal are 0's and 1's

        # now, use the pivot A[i, i] to eliminate entries above and below it
        for j = 1 : row
            if j == i || A[j, i] == 0; continue; end
            A[j, i : end] -= A[j, i] * A[i, i : end]
        end
    end
    return A
end

#
# Algorithm
# ---------
# General generator of finite difference formulas for the n-th order derivatives
# of f(x) at x[i] using points in a soredt list, e.g., [-2 -1 0 1 3 5 6 7].
#
# It uses the linear combination of f(x[i+j]), j in points, to eliminate
# f(x[i]), f'(x[i]), ..., so that the first term of the Taylor series expansion
# of the linear combination is f^(n)(x[i]):
#
#  k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ... + k[len]*f(x[i+points[len]])
#     = m*f^(n)(x[i]) + ..., m > 0, len = length(points)
#
# It is this equation that gives the formula for computing f^(n)(x[i]).
#
# Values of the cofficients k[:] and m will be determined, and the first few
# terms of the remainder will be listed for calculating the truncation error of
# a formula.
#
# Julia's Rational type and related arithmetic operations, especially a // b,
# are a perfect fit for obtaining an "exact" formula.
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
# The function returns a tuple, ([k[1], k[2], ..., k[len]], m).
#
function _computecoefs(n::Int, points::Vector{Int}, printformulaq::Bool = false)
    global _max_num_of_taylor_terms, _default_max_num_of_taylor_terms
    global _lcombination_coefs

    global _computedq = false
    global _formula_status = 0
    global _bigO = ""

    # for teaching's purpose, we don't do so
    # if length(points) <= n
    #     pts = _range_inputq ? "$(_range_input)" : "$(points')"
    #     th = n == 1 ? "st" : (n == 2 ? "nd" : (n == 3 ? "rd" : "th"))
    #     println("$pts is invalid because at least $(n + 1) points are ",
    #             "needed for the $n$th derivative.")
    #     return
    # end

    # setup a linear system Ax = B first
    len = length(points)
    _max_num_of_taylor_terms =
        max(_default_max_num_of_taylor_terms, max(len, n) + 8)
    _lcombination_coefs = Array{Any}(undef, _max_num_of_taylor_terms)

    # setup the coefficients of Taylor series expansions of f(x) at each of the
    # involved points
    coefs = Array{Any}(undef, len)
    for i in 1 : len
        coefs[i] = _taylor_coefs(points[i])
    end

    # We find a linear combination of
    # f(x[i+points[1]]), f(x[i+points[2]]), ..., f(x[i+points[len]]),
    #
    # k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ... + k[len]*f(x[i+points[len]])
    #  = 0*f(x[i]) + 0*f'(x[i]) + ... + 0*f^(n-1)(x[i]) + m*f^(n)(x[i]) + ..., m != 0
    #
    # so that it must eliminate f(x[i]), f'(x[i]), ..., f^(n-1)(x[i]); given
    # more points, it also eliminates f^(n+1)(x[i]), f^(n+2)(x[i]), ....
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
    #  k[1]*coefs[1][j] + k[2]*coefs[2][j] + ... + k[len]*coefs[len][j] = 0 ... (1)
    #
    # where j = 1, 2, ..., len but not n.
    #
    A = Matrix{Rational{BigInt}}(undef, len, len)
    row = 1
    for order in 0 : len - 1  # correspond to f(x[i]), f'(x[i]), f''(x[i]), ...
        if order == n; continue; end     # skip f^(n)(x[i])

        # eliminating f^(order)(x[i])
        # A[:,j] stores coefs of
        #  k[1]*coefs[1][j] + k[2]*coefs[2][j] + ... + k[len]*coefs[len][j] = 0
        for j = 1 : len
            A[row, j] = coefs[j][order + 1]
        end
        row += 1
    end

    # The homogeneous linear system (1) has no nontrivial solution or has
    # inifinitely many nontrivial solutions. It is understandable that it may
    # not have a nontrivial solution. But why inifinitely many nontrivial
    # solutions? It is because, if k[:] is a nontrivial solution, α k[:] is
    # also a nontrivial solution, where α is any nonzero real constant, i.e.,
    # all nontrivial solutions (a subspace spanned by this k[:]) are parallel
    # to each other. Therefore, in the case that there are infinitely many
    # nontrivial solutions, if we know one entry k[which] is nonzero and let
    # it be a nonzero constant (say, 1), then, a nontrivial solution k[:] is
    # uniquely determined.
    #
    # Beware, there may be multiple nontrivial solutions,
    #  k[:] = [k1, k2, ...], k[:] = [K1, K2, ...], ..., or k[:] = [κ1, κ2, ...]
    # of which no two are parallel to each other. However, each of these
    # solutions must satisfy the condition that both k[1] and k[end] can't be
    # zero. Why? If, say, k[1] = 0, in other words, a formula only uses/depends
    # on (at most) x[i+points[2]], x[i+points[3]], ..., x[i+points[len]], why
    # should we say it is a formula that uses/depends on x[i+points[1]] (and
    # x[i+points[2]], ..., x[i+points[len]])? Therefore, we can assume that
    # k[1] != 0.
    #
    A[len, 2 : len] .= 0                   # setup an equation so that k[1] = 1
    A[len, 1]        = 1

    B = zeros(Rational{BigInt}, len, 1)
    B[len] = 1                             # so that k[1] = 1

    # solve Ax = B for x, i.e., k[:]
    k = _rref([A B])[:, len + 1]
    k = k // gcd(k)

    # Taylor series expansion of the linear combination
    # k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ... + k[len]*f(x[i+points[len]])
    _lcombination_coefs = k[1] * coefs[1]  # let Julia determine the type
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
    global _data = _FDData(n, points, k, m, coefs, _max_num_of_taylor_terms)
    _computedq = true

    _test_formula_validity()

    if printformulaq; formula(); end

    return (coefs, m)
end   # _computecoefs

# return a string of the linear combination
# k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ... + k[len]*f(x[i+points[len]])
#
# To test in Julia if a formula is valid, we must convert the data to BigInt or
# BigBloat becasue, say, for computecoefs(2, -12:12), without the conversion,
# the formula will give obviously wrong answers because of overflow and
# rounding errors. julia_REPL_funcq is here for testing in Julia REPL only.
# For ordinary definition of a function, call formula().
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

# check if the newly computed formula is valid. results are saved in the global
# variable _formula_status:
#  100 - Perfect, even satifiying some commonly seen "rules", such as the sum of
#        coefficients = 0, symmetry of cofficients about x[i] in a central formula
# -100 - No formula can't be found
# > 0  - Mathematically, the formula is still valid but may not satisfy some
#        commonly seen "rules" such as the sum of coefficients = 0 and symmetry
#        of coefficients about x[i]
#
function _test_formula_validity()
    # to find f^(n)(x[i]) by obtaining
    #
    # k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ... + k[len]*f(x[i+points[len]])
    #   = m*f^(n)(x[i]) + ..., m > 0
    #
    # the most important step is to know if f(x[i]), f'(x[i]), ..., f^(n-1)(x[i])
    # are all eliminated, i.e.,
    #    k[1]*coefs[1][j] + k[2]*coefs[2][j] + ... + k[len]*coefs[len][j] = 0
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
            println("***** Error:: $n, $input_points :: i = $i, k[1]*coefs[1][$i]",
                    " + k[2]*coefs[2][$i] + ... + k[$len]*coefs[$len][$i] != 0")
            has_solutionq = false
            break
        end
    end

    # Is m = 0 ?
    if has_solutionq
        m = _lcombination_coefs[n + 1]
        if m == 0
            println("-" ^ 81)
            println("\n***** Error:: $n, $input_points :: m = 0, ",
                    "formula can't be found.")
            has_solutionq = false
        end
    end

    # there could be a solution that doesn't use the whole range of input (points)
    formula_for_inputq = true
    if has_solutionq
        start, stop = 1, len    # actual left and right endpoints
        while start <= len && k[start] == 0; start += 1; end
        while stop >= 1 && k[stop] == 0; stop -= 1; end
        if start >= stop
            has_solutionq = false
        elseif start > 1 || stop < len
            s = _range_inputq ? (points[start] : points[stop]) : points[start : stop]
            println("***** Error:: $n, $input_points :: can't find a formula.\n")
            println("***** Warning:: $n, $s might be your input for which a ",
                    "formula is found.\n")
            formula_for_inputq = false
        end
    end

    if !has_solutionq
        if len <= n
            th = n == 1 ? "st" : (n == 2 ? "nd" : (n == 3 ? "rd" : "th"))
            println("***** Error:: $n, $input_points :: Invalid input because",
                    " at least $(n + 1) points are needed for the $n$th ",
                    "derivative.\n")
        else
            println("***** Error:: $n, $input_points :: can't find a formula.")
        end
        return
    end

    if sum(k) != 0   # sum of coefficients must be 0
        println("***** Warning:: $n, $input_points :: sum(k[:]) != 0")
        _formula_status += 1
    end

    # are coefficients of central formulas symmetrical about x[i]?
    if formula_for_inputq && _range_inputq && abs(_range_input.start) == _range_input.stop
        j = length(k)
        for i in 1 : round(Int64, length(k)/2)
            if abs(k[i]) != abs(k[j])
                println("***** Warning:: $n, $input_points :: k[$i] != k[$j]")
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
    global _data, _computedq, _formula_status
    global _range_inputq, _range_input, _bigO, _julia_func_basename

    if !_computedq
        println("Please call computecoefs(n, points) first!")
        return
    end

    # it can be changed by talor and/or printtaylor; change it back
    global _max_num_of_taylor_terms = _data.max_num_of_taylor_terms

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
            println("may still be valid, though it didn't pass tests like sum ",
                    "of the coefficients being zero.\n")
        end
    end

    # print Taylor series expansion of the linear combination:
    # k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ... + k[len]*f(x[i+points[len]])
    println("The computing result:\n")
    print(_lcombination_expr(_data), " =\n")
    _print_taylor(_lcombination_coefs, 5);       # print at most 5 nonzero terms

    if _formula_status > 0
        println("The exact formula:\n")
        # find x in the big-O notation O(h^x)
        x = _data.max_num_of_taylor_terms
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
            data1 = _FDData(_data.n, _data.points, _data.k // _data.m, 1,
                            _data.coefs, _data.max_num_of_taylor_terms)
            print("Or\n\n")
            _print_bigo_formula(data1, _bigO)
        end

        println("Julia function:\n")
        global _julia_exact_func_expr = _julia_func_expr(_data)
        print(_julia_func_basename, "e", _julia_exact_func_expr, "\n\n")
        if _data.m > 1           # other formats
            global _julia_exact_func_expr1  = _julia_func_expr(data1)
            print("Or\n\n", _julia_func_basename, "e1",
                  _julia_exact_func_expr1, "\n\nOr\n\n")
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
            println("No valid formula is available. Please call ",
                    "'computecoefs' again.")
        else
            println(_bigO)
        end
    else
        println("Please call 'computecoefs' first.")
    end
    return
end

############################ for teaching/learning ############################

# show current decimal places
function decimalplaces()
    global _decimal_places
    return _decimal_places
end  # decimalplaces

# set decimal places to n
function decimalplaces(n)
    global _decimal_places, _computedq
    if isinteger(n) && n > 0
        _decimal_places = n
        if _computedq
            println("Please call 'formula' to generate (or ",
                    "'activatejuliafunction' to generate and activate) a ",
                    "Julia function for the newly computed formula, using ",
                    "the new decimal places.")
        else
            println("You may start your work by calling 'computecoefs'.")
        end
    else
        println("decimalplaces(n): n must be a positive integer.")
    end

    return
end  # decimalplaces

function _printexampleresult(suffix, exact)
    global _julia_func_basename

    apprx = eval(Meta.parse("$(_julia_func_basename)$(suffix)(f, x, i, h)"))
    relerr = abs((apprx - exact) / exact) * 100
    print("  fd.$(_julia_func_basename)$(suffix)(f, x, i, h)  ",
          suffix == "e1" ? "" : " ", "# result: ")
    @printf("%.16f, ", apprx)
    print("relative error = "); @printf("%.8f", relerr); println("%")
end

# activate function(s) for the newly computed finite difference formula,
# allowing immediate evaluation of the formula in Julia REPL
function activatejuliafunction()
    global _computedq, _formula_status, _data, _julia_func_basename

    if !_computedq
        println("Please run 'computecoefs' first.")
        return
    end

    # it can be changed by talor and/or printtaylor; change it back
    global _max_num_of_taylor_terms = _data.max_num_of_taylor_terms

    global _julia_exact_func_expr   = ""
    global _julia_exact_func_expr1  = ""
    global _julia_decimal_func_expr = ""

    # generate Julia function for current REPL session
    if _formula_status > 0
        global _julia_exact_func_expr = _julia_func_expr(_data, false, true)
        if _data.m > 1           # print in other formats
            data1 = _FDData(_data.n, _data.points, _data.k // _data.m, 1,
                            _data.coefs, _data.max_num_of_taylor_terms)
            _julia_exact_func_expr1  = _julia_func_expr(data1, false, true)
            _julia_decimal_func_expr = _julia_func_expr(data1, true, true)
        end  # m = 1? No decimal formula
    else
        print("No valid formula is for activation. Please run 'computecoefs' ",
              "with different input again.")
        return
    end

    eval(Meta.parse(_julia_func_basename * "e" * _julia_exact_func_expr))
    if _julia_decimal_func_expr != ""
        eval(Meta.parse(_julia_func_basename * "e1" * _julia_exact_func_expr1))
        eval(Meta.parse(_julia_func_basename * "d" * _julia_decimal_func_expr))
    end

    println("The following function(s) are available temporiarily in the ",
            "FiniteDifferenceFormula module. Usage:\n")
    println("  import FiniteDifferenceFormula as fd\n")
    println("  f, x, i, h = sin, 0:0.01:10, 501, 0.01")

    fmt = Printf.Format("%.8f")      # format of relative error

    # sine is taken as the example b/c sin^(n)(x) = sin(n π/2 + x), simply
    exact = sin(_data.n * pi /2 + 5) # x[501] = 5
    eval(Meta.parse("f, x, i, h = sin, 0:0.01:10, 501, 0.01"))

    _printexampleresult("e", exact)
    if _julia_decimal_func_expr != ""
        _printexampleresult("e1", exact)
        _printexampleresult("d", exact)
    end
    len = length("fd.$(_julia_func_basename)")
    print(" "^(len + 18) * "# cp:     ")
    @printf("%.16f", exact)
    println("\n\nCall fd.formula() to view the very definition.")

    return  # stop Julia from returning something users never expect
end  # activatejuliafunction

# calculate the coefficients of Taylor series of f(x[i + j]) about x[i]
# for teaching/learning!
function taylor(j::Int)
    global _max_num_of_taylor_terms, _default_max_num_of_taylor_terms
    # computecoefs could have changed the default value; change it back
    _max_num_of_taylor_terms = _default_max_num_of_taylor_terms
    return _taylor_coefs(j)
end  # taylor

# print readable Taylor series expansion of f(x[i + j]) about x[i]
# for teaching/learning!
function printtaylor(j::Int, num_of_nonzero_terms::Int = 10)
    global _max_num_of_taylor_terms, _default_max_num_of_taylor_terms
    # computecoefs could have changed the default value; change it back
    _max_num_of_taylor_terms = _default_max_num_of_taylor_terms
    if num_of_nonzero_terms > _max_num_of_taylor_terms
        println("can't print $num_of_nonzero_terms terms of the ",
                "Taylor series.\n")
        num_of_nonzero_terms = _max_num_of_taylor_terms
    end
    coefs = taylor(j)
    println("\nf(x[i" * (j == 0 ? "" : (j > 0 ? "+$j" : "$j")) * "]) =")
    _print_taylor(coefs, num_of_nonzero_terms)
    return
end  # printtaylor

# print readable Taylor series of a function/expression about x[i]. e.g.,
# fd.printtaylor(2*fd.taylor(0) - 5*fd.taylor(1) + 4*fd.taylor(2))
function printtaylor(coefs::Matrix{Rational{BigInt}},
                     num_of_nonzero_terms::Int = 10)
    global _max_num_of_taylor_terms, _default_max_num_of_taylor_terms
    # computecoefs could have changed the default value; change it back
    _max_num_of_taylor_terms = _default_max_num_of_taylor_terms
    if num_of_nonzero_terms > _max_num_of_taylor_terms
        println("can't print $num_of_nonzero_terms terms of the Taylor series.\n")
        num_of_nonzero_terms = _max_num_of_taylor_terms
    end
    _print_taylor(coefs, num_of_nonzero_terms)
    return
end  # printtaylor

# for explorers/researchers. users never need to know its existence
function _set_default_max_num_of_taylor_terms(n::Int)
    global _default_max_num_of_taylor_terms
    if n < 0
        println("Wrong input, $n. A positive integer is expected.")
    elseif n < 10
        println("$n? It doesn't have to be so small. The default value ",
                "is still $_default_max_num_of_taylor_terms.")
    else
        _default_max_num_of_taylor_terms = n
    end
    return
end

end # module
