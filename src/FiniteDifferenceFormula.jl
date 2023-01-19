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

# Warning: users should not call/access a function/variable starting with "_" !

using Printf

############################# EXPORTED FUNCTIONS ##############################
export compute, formula, truncationerror, verifyformula

# for teaching/learning/exploring
export decimalplaces, activatejuliafunction, taylor, printtaylor
export _set_default_max_num_of_taylor_terms

######################### BEGIN OF GLOBAL VARIABLES ###########################
_decimal_places::Int = 16     # use it to print Julia function for a formula
                              # call decimalplaces(n) to reset it

mutable struct _FDData
    n; points; k; m; coefs                 # on one line? separated by ;
end

_data                            = _FDData # share results between functions
_computedq::Bool                 = false   # make sure compute is called first
_formula_status::Int             = 0       # a formula may not be available
                                           # values? see _test_formula_validity()

# a vector of the coefficients of Taylor series of the linear combination:
# k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ... + k[len]*f(x[i+points[len]])
_lcombination_coefs              = Array{Any}

_range_inputq::Bool              = false
_range_input::UnitRange{Int}     = 0:0     # compute receives a range? save it

_julia_exact_func_expr::String   = ""      # 1st exact Julia function for f^(n)(x[i])
_julia_exact_func_expr1::String  = ""      # 2nd exact Julia function for f^(n)(x[i])
_julia_decimal_func_expr::String = ""      # decimal Julia function for f^(n)(x[i])
_julia_func_basename::String     = ""

_bigO::String                    = ""      # truncation error of a formula
_bigO_exp::Int                   = -1      # the value of n as in O(h^n)
########################### END OF GLOBAL VARIABLES ###########################

# for future coders/maintainers of this package:
# to compute a new formula, this function must be called first.
function _initialization()
    global _data                    = nothing
    global _computedq               = false
    global _formula_status          = 0
    global _lcombination_coefs      = nothing

    global _range_inputq            = false
    global _range_input             = 0:0

    global _julia_exact_func_expr   = ""
    global _julia_exact_func_expr1  = ""
    global _julia_decimal_func_expr = ""
    global _julia_func_basename     = ""

    global _bigO                    = ""
    global _bigO_exp                = -1
end  # _initialization

# This function returns the first 'max_num_of_terms' terms of Taylor series of
# f(x[i+1]) centered at x=x[i] in a vector with f(x[i]), f'(x[i]), ..., removed.
# The purpose is to obtain Taylor series expansions for f(x[i±k]) = f(x[i]±kh])
# which are used to derive the m-point finite difference formulas for the first,
# second, ..., order derivatives at x[i].
#
#       f(x[i+1]) = f(x[i]) + 1/1! f'(x[i]) h + 1/2! f''(x[i]) h^2 + ...
#
# where h = x[i+1] - x[i].
#
# Usage:
#   for f(x[i+k]), call _taylor_coefs(k, ...), k = ±0, ±1, ±2, ...
function _taylor_coefs(h, max_num_of_terms)
    result = Matrix{Rational{BigInt}}(undef, 1, max_num_of_terms)
    factorial::BigInt = 1
    for n in 1 : max_num_of_terms
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
# Input: An array that contains the coefficients of the first terms of Taylor
#        series expansion of a function
function _print_taylor(coefs, num_of_nonzero_terms = 10)
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

function _dashline(); return "-" ^ 105; end

function compute(n, points, printformulaq = false)
    if !isinteger(n) || n < 1
        println("Invalid order of derivatives, $n. A positive integer ",
                "is expected.")
        return
    end
    n = round(Int, n)  # 4.0 --> 4

    if (typeof(points) <: Tuple) || !(typeof(points[1]) <: Integer)
        println("Invalid input, $points. A list of integers like [-2, 1, ...]",
                " is expected.")
        return
    end
    if typeof(printformulaq) != Bool
        println("Invalid input, $printformulaq. A value, false or true, is ",
                "expected.")
        return
    end

    oldlen = length(points)
    points = sort(unique(collect(points)))
    len = length(points)
    if len < 2
        println("Invalid input, $points. A list of two or more different ",
                "points is expected.")
        return
    end

    _initialization()
    if len != oldlen
        input_points = _format_of_points(points)
        print(_dashline(), "\nYour input is converted to ($n, $input_points")
        if printformulaq; print(", true"); end
        println(").\n", _dashline())
    else
        _format_of_points(points)
    end
    return _compute(n, points, printformulaq)
end  # compute

# before v1.0.7, we thankfully used the function 'rref' provided by the package
# RowEchelon v0.2.1. since it has been removed from the base, it is safer for
# this package to have its own implementaion of the function. anyway, it is
# just a few lines of code.
#
# customized code from RowEchelon v0.2.1, also vectorized, simplified 1/10/23
function _rref(A::Matrix{Rational{BigInt}})
    nr, nc = size(A)
    i = j = 1                     # A[i, j], pivot entries
    while i <= nr && j <= nc
        (m, mi) = findmax(abs.(A[i : nr, j]))
        if m != 0
            mi += i - 1
            if mi != i
                A[i, j : nc], A[mi, j : nc] = A[mi, j : nc], A[i, j : nc]
            end
            J = j + 1               #
            A[i, J : nc] /= A[i, j] # A[i, j : nc] /= A[i, j] depends on how
            A[i, j] = 1             # Julia implements inside; not a safe way
            for r = 1 : nr
                if r == i || A[r, j] == 0; continue; end
                A[r, J : nc] -= A[r, j] * A[i, J : nc] # same as above
                A[r, j] = 0                            #
            end
            i += 1
        end
        j += 1
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
#     = m*f^(n)(x[i]) + ..., m > 0, len = length(points)     .............. (1)
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
# The function returns a tuple, (n, points, k[:], m]).
# where, m, k? refer to eq (1); n, points? refer to _compute()
#
function _compute(n::Int, points::Vector{Int}, printformulaq::Bool = false)
    global _lcombination_coefs, _formula_status, _range_inputq, _range_input

    # for teaching's purpose, we don't do so
    # if length(points) <= n
    #     pts = _range_inputq ? "$(_range_input)" : "$(points')"
    #     println("$pts is invalid because at least $(n + 1) points are ",
    #             "needed for the $(_nth(n)) derivative.")
    #     return
    # end

    # setup a linear system Ax = B first
    len = length(points)
    max_num_of_terms = max(len, n) + 8

    # setup the coefficients of Taylor series expansions of f(x) at each of the
    # involved points
    coefs = Array{Any}(undef, len)
    for i in 1 : len
        coefs[i] = _taylor_coefs(points[i], max_num_of_terms)
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
    #  k[1]*coefs[1][j] + k[2]*coefs[2][j] + ... + k[len]*coefs[len][j] = 0 ... (2)
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

    # The homogeneous linear system (2) has no nontrivial solution or has
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

    k = k // gcd(k)                        # change each element to an integer
    # the following code does the same     # for translating to other languages
    ## for i in 1 : len
    ##    if k[i] != round(BigInt, k[i]); k *= denominator(k[i]); end
    ## end

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
    if x == m; m = x; end                  # already integer; don't show like 5//1

    # save the results in a global variable for other functions
    global _data = _FDData(n, points, k, m, coefs)
    global _computedq = true

    _test_formula_validity()

    if printformulaq; formula(); end

    if _formula_status >= 0
        return (n, _range_inputq ? _range_input : points, round.(BigInt, k), m)
    else
        return nothing
    end
end   # _compute

# return a string of the linear combination
# k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ... + k[len]*f(x[i+points[len]])
#
# To test in Julia if a formula is valid, we must convert the data to BigInt or
# BigBloat becasue, say, for compute(2, -12:12), without the conversion,
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
            s *= "BigFloat(" * c2s * ")"
        else
            s *= _c2s(data.k[i], firstq, decimalq)
        end
        s *= times * _f2s(data.points[i])
        firstq = false
    end
    return s
end  # _lcombination_expr

# return string of 1st, 2nd, 3rd, 4th ...
function _nth(n::Int)
    th::String = n == 1 ? "st" : (n == 2 ? "nd" : (n == 3 ? "rd" : "th"))
    return "$n$th"
end

# check if the newly computed formula is valid. results are saved in the global
# variable _formula_status:
#  100 - perfect, even satifiying some commonly seen "rules", such as the sum of
#        coefficients = 0, symmetry of cofficients about x[i] in a central formula
# -100 - no formula can't be found
# -200 - same as -100, but do not try to find a formula if
#        activatejuliafunction(n, point, k, m) fails
#  250 - same as -100, but used for communication btw 2 'activatejuliafunction's
# > 0  - mathematically, the formula is still valid but may not satisfy some
#        commonly seen "rules" such as the sum of coefficients = 0 and symmetry
#        of coefficients about x[i]
#
# return m as in equation (1) for 'activatejuliafunction'
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
    global _data, _lcombination_coefs, _range_inputq, _range_input

    n = _data.n
    k = _data.k
    coefs = _data.coefs
    points = _data.points
    len = length(points)

    input_points = _range_inputq ? _range_input : points'

    # Is there any equation in system (2) that is not satisfied?
    has_solutionq = true
    global _formula_status = 0
    for i = 1 : n
        if _lcombination_coefs[i] != 0
            println("***** Error: $n, $input_points : i = $i, k[1]*coefs[1][$i]",
                    " + k[2]*coefs[2][$i] + ... + k[$len]*coefs[$len][$i] != 0")
            has_solutionq = false
            break
        end
    end

    m = _lcombination_coefs[n + 1]
    if m == 0; has_solutionq = false; end

    # there could be a solution that doesn't use the whole range of input (points)
    formula_for_inputq = true
    if has_solutionq
        start, stop = 1, len    # actual left and right endpoints
        while start <= len && k[start] == 0; start += 1; end
        while stop >= 1 && k[stop] == 0; stop -= 1; end
        #if start >= stop  # can't be true b/c m != 0
        #    has_solutionq = false
        #    # what's the error?
        #else
        if start > 1 || stop < len
            s = _range_inputq ?
                (points[start] : points[stop]) : points[start : stop]
            println("***** Warning: $n, $s might be your input for which a ",
                    "formula is found.\n")
            formula_for_inputq = false
        end
    end

    if !has_solutionq
        if len <= n
            println("***** Error: $n, $input_points : Invalid input because",
                    " at least $(n + 1) points are needed for the $(_nth(n)) ",
                    "derivative.\n")
            _formula_status = -200
            return m
        end
        println("\n***** Error: $n, $input_points : can't find a formula.\n")
        _formula_status = -100
        return m
    end

    if sum(k) != 0   # sum of coefficients must be 0
        println("***** Warning: $n, $input_points : sum(k[:]) != 0")
        _formula_status += 1
    end

    # are coefficients of central formulas symmetrical about x[i]?
    if formula_for_inputq && _range_inputq && abs(_range_input.start) == _range_input.stop
        j = length(k)
        for i in 1 : round(Int64, length(k)/2)
            if abs(k[i]) != abs(k[j])
                println("***** Warning: $n, $input_points : k[$i] != k[$j]")
                _formula_status += 1
                break
            end
            j -= 1
        end
    end

    if _formula_status == 0
        _formula_status = 100    # perfect
    end

    # now, determine the big-O notation - what is x in O(h^x)?
    x = length(_lcombination_coefs)
    for i in n + 2 : x            # skip f^(n)(x[i])
        if _lcombination_coefs[i] != 0; x = i; break; end
    end
    x -= n + 1
    global _bigO = "O(h"
    if x > 1; _bigO *= "^$x"; end
    _bigO *= ")"
    global _bigO_exp = x

    return m
end  # _test_formula_validity

function _denominator_expr(data::_FDData, julia_REPL_funcq::Bool = false)
    ms = ""
    if typeof(data.m) <: Rational
        ms = "$((data.m).num)/$((data.m).den)"
    else
        ms = "$(data.m)"
    end
    s  = data.m != 1 ? "($ms * " : ""
    if julia_REPL_funcq  # v1.1.1
        s *= "BigFloat(h)"    # v1.1.1, doesn't matter to Julia 1.8.x
    else
        s *= "h"
    end
    if data.n != 1; s *= "^$(data.n)"; end
    if data.m != 1; s *= ")"; end
    return s
end  # _denominator_expr

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

    global _julia_func_basename = "f$(_nth(data.n))deriv$(n)pt$(s)"
    fexpr  = "(f, x, i, h) = "
    if julia_REPL_funcq; fexpr *= "Float64("; end
    fexpr *= "( "
    fexpr *= _lcombination_expr(data, decimalq, julia_REPL_funcq)
    fexpr *= " ) / " * _denominator_expr(data, julia_REPL_funcq)
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
    global _data, _computedq, _formula_status, _lcombination_coefs
    global _range_inputq, _range_input, _bigO, _julia_func_basename

    if !_computedq
        println("Please call compute(n, points) first!")
        return
    end

    if _formula_status > 0
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
    println("Computing result:\n")
    print(_lcombination_expr(_data), " =\n")
    _print_taylor(_lcombination_coefs, 5);   # print at most 5 nonzero terms

    if _formula_status > 0
        println("The exact formula:\n")
        _print_bigo_formula(_data, _bigO)
        data1 = _FDData
        if _data.m != 1                 # print in another format
            data1 = _FDData(_data.n, _data.points, _data.k // _data.m, 1,
                            _data.coefs)
            print("Or\n\n")
            _print_bigo_formula(data1, _bigO)
        end

        println("Julia function:\n")
        global _julia_exact_func_expr = _julia_func_expr(_data)
        print(_julia_func_basename, _data.m != 1 ? "e" : "",
              _julia_exact_func_expr, "\n\n")
        if _data.m != 1                 # other formats
            global _julia_exact_func_expr1  = _julia_func_expr(data1)
            print("Or\n\n", _julia_func_basename, "e1",
                  _julia_exact_func_expr1, "\n\nOr\n\n")
            global _julia_decimal_func_expr = _julia_func_expr(data1, true)
            print(_julia_func_basename, "d", _julia_decimal_func_expr, "\n\n")
        end
    #else
    #    _formula_status = -100         # no formula can't be generated
    end

    return
end  # formula

# print _bigO, the estimation of trucation error in the big-O notation
#
# output:
#    (-1, "")      - There is no valid formula
#    (n, "O(h^n)") - There is a valid formula
function truncationerror()
    global _bigO, _bigO_exp, _computedq, _formula_status
    if _computedq
        if _formula_status <= -100
            println("No valid formula is available.")
            return (-1, "")
        else
            return (_bigO_exp, _bigO)
        end
    end
    println("Please call 'compute' first.")
    return (-1, "")
end  # truncationerror

######################## for teaching/learning/exploring #######################

_max_num_of_taylor_terms = 30
# it should be a constant. usually, users never need to change it.

# show current decimal places
function decimalplaces()
    global _decimal_places
    return _decimal_places
end  # decimalplaces

# set decimal places to n
function decimalplaces(n)
    global _decimal_places, _computedq
    if isinteger(n) && n >= 0
        _decimal_places = round(Int, n) # 10.0 --> 10
        if _computedq
            println("Please call 'formula' to generate (or ",
                    "'activatejuliafunction' to generate and activate) a ",
                    "Julia function for the newly computed formula, using ",
                    "the new decimal places.")
        else
            println("You may start your work by calling 'compute'.")
        end
    else
        println("decimalplaces(n): a nonnegative integer is expected.")
    end

    return _decimal_places
end  # decimalplaces

function _format_of_points(points)
    global _range_inputq, _range_input
    if length(points) == length(points[1] : points[end])
        global _range_inputq = true
        global _range_input  = points[1] : points[end]
        return _range_input
    end
    return points
end  # _format_of_points

function _printexampleresult(suffix, exact)
    global _julia_func_basename

    apprx = eval(Meta.parse("$(_julia_func_basename)$(suffix)(f, x, i, h)"))
    relerr = abs((apprx - exact) / exact) * 100
    print("  fd.$(_julia_func_basename)$(suffix)(f, x, i, h)  ",
          suffix == "e1" ? "" : " ", "# result: ")
    @printf("%.16f, ", apprx)
    print("relative error = "); @printf("%.8f", relerr); println("%")
end  # _printexampleresult

# the name is self-explanatory. it is exactly the same as the function
# activatejuliafunction(n, points, k, m)
function verifyformula(n, points, k, m)
    return activatejuliafunction(n, points, k, m)
end

# if you have data from textbooks or other sources, you may use this function
# to verify if it is right, activate related Julia function(s), evaluate and
# see the computiong results.
#
# it seemed to be very useful when I tried to port this package to Python
# (3.11.1, the newest version as of 1/13/2023). the effort failed in hours
# becasue Python's SymPy might not be able to handle very large integers.
# for example,
#   compute(3,[0,1,2,3,6,8,9,10,11,12,13,14,15,16,17,18,19]) .............. (3)
# though it could handle
#   compute(3,[0,1,2,3,6,8,9,10,11,12,13,14,15,16,17,18])
#
# 1/19/23: removed SymPy, include 'from fractions import Fraction', treated a
# "matrix" as a list of lists for which a special 'rref' was coded. a big
# progress! however, still failed to do, say, compute(1, -30:30).
#
# thus, the functionality of Julia's BigInt and Rational is simply amazing!!!
#
# while Python's output for command (3) was so different, I wanted to load it
# to some function here to test and evaluate. here it is.
#
function activatejuliafunction(n, points, k, m)
    if !isinteger(n) || n <= 0
        println("Error: invalid first argument, $n. A positive integer is ",
                "expected.")
        return
    end
    n = round(Int, n)  # 4.0 --> 4
    if length(points) == 1
        println("Error: invalid input, points = $points. A list of two or ",
                "more points is expected.")
        return
    end
    if !(typeof(points[1]) <: Integer)
        println("Error: invalid input, $points. Integers are expected.")
        return
    end

    # tuple input can cause trouble - Julia 1.8.5: x = (1, 2//3), y = collect(x)
    # typeof(y[1]) != typeof(y[2]) !!!
    if typeof(points) <: Tuple
        println("Invalid input, $points. A list like [1, 2, ...] is expected")
        return
    end
    if typeof(k) <: Tuple
        println("Invalid input, $k. A list like [1, 2, ...] is expected")
        return
    end

    oldlen = length(points)
    points = sort(unique(collect(points)))
    len = length(points)

    rewrittenq = oldlen != len
    # don't do so for teaching
    #if n <= len
    #    println("Error: at least $(n+1) points are needed for the $(_nth(n))",
    #            " derivative.")
    #    return
    #end
    if len != length(k)
        println("Error: The number of points != the number of coefficients.");
        return
    end

    _initialization()      # needed b/c it's like computing a new formula
    input_points = _format_of_points(points)

    # "normalize" input so that m > 0, and m is integer
    if m < 0; k *= -1; m *= -1; rewrittenq = true; end
    if isinteger(m)
        m = round(Int, m)  # 5.0 --> 5
    else
        if !(typeof(m) <: Rational)
            m = rationalize(convert(Float64, m))
        end
        k         *= m.den
        m          = m.num
        rewrittenq = true
    end
    if m == 0
        if rewrittenq; println("You input: $n, $input_points, $k, $m."); end
        println("Error: invalid input, the last argument m = 0. ",
                "It can't be zero.")
        return
    end

    # "normalize" k[:] so that each element is integer
    if !((typeof(k[1]) <: Integer) || (typeof(k[1]) <: Rational))
        k = rationalize.(convert.(Float64, k))
    end
    if typeof(k[1]) <: Rational
        for i in 1 : len
            if k[i].den == 1; continue; end
            m         *= k[i].den
            k         *= k[i].den
            rewrittenq = true
        end
        k = round.(BigInt, k)
    end

    if rewrittenq
        # print k[:] nicely
        ks = "[$(k[1])"
        for i in 2 : len; ks *= ", $(k[i])"; end
        ks *= "]"
        println(_dashline(), "\nYour input is converted to ($n, $input_points, ",
                "$ks, $m).\n", _dashline())
    end

    # setup the coefficients of Taylor series expansions of f(x) at each of the
    # involved points
    max_num_of_terms = max(len, n) + 8
    coefs = Array{Any}(undef, len)
    for i in 1 : len
        coefs[i] = _taylor_coefs(points[i], max_num_of_terms)
    end

    # Taylor series of the linear combination
    # k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ... + k[len]*f(x[i+points[len]])
    global _lcombination_coefs = k[1] * coefs[1]  # let Julia determine the type
    for i in 2 : len
        if k[i] == 0; continue; end
        _lcombination_coefs += k[i] * coefs[i]
    end

    global _data = _FDData(n, points, k, m, coefs)
    global _formula_status
    M = _test_formula_validity()
    find_oneq::Bool = _formula_status == -100

    # perhaps, the coefficients k[:] is correct, but m is not
    if _formula_status == 100 && M != m
        x = round(BigInt, M)
        ms = "$M"
        if M == x
            ms = "$x"              # don't print 5//1
        elseif typeof(M) <: Rational
            ms = "$(M.num)/$(M.den)"
        end
        println("***** Error: The last argument m = $m is incorrect. ",
                "It should be $ms.\n")
        find_oneq = true
    end

    if find_oneq                   # the input can't be a formula, but still
        _formula_status = 250      # force to activate Julia function
        # 250, reserved for communication w/ another 'activatejuliafunction'
    end
    global _computedq   = true     # assume a formula has been computed

    # force to activate even for invalid input formula
    b200::Bool = _formula_status == -200
    if b200; _formula_status = 100; end;   # assume it is valid
    activatejuliafunction(true)
    if b200; _formula_status = -200; end;  # change it back

    if find_oneq                   # use the basic input to generate a formula
        _formula_status = -100
        println(_dashline(), "\nFinding a formula using the points....")
        result = _compute(n, points)
        if _formula_status >= 0
            println("Call fd.formula() to view the results and fd.activatejulia",
                    "function() to activate the new Julia function(s).")
            return result
        end
    end
    return
end  # activatejuliafunction

# activate function(s) for the newly computed finite difference formula,
# allowing immediate evaluation of the formula in Julia REPL
function activatejuliafunction(external_dataq = false)
    global _computedq, _formula_status, _data
    global _julia_func_basename, _julia_exact_func_expr
    global _julia_exact_func_expr1, _julia_decimal_func_expr

    if !external_dataq && !_computedq
        println("Please run 'compute' first.")
        return
    end

    # generate Julia function in/for current REPL session
    count = 1
    if _formula_status > 0
        _julia_exact_func_expr = _julia_func_expr(_data, false, true)
        if _data.m != 1           # print in other formats
            data1 = _FDData(_data.n, _data.points, _data.k // _data.m, 1,
                            _data.coefs)
            _julia_exact_func_expr1  = _julia_func_expr(data1, false, true)
            _julia_decimal_func_expr = _julia_func_expr(data1, true, true)
            count = 3
        end  # m = 1? no decimal formula
    else
        print("No valid formula is for activation.")
        return
    end

    eval(Meta.parse(_julia_func_basename * (count == 1 ? "" : "e") * _julia_exact_func_expr))
    if count == 3
        eval(Meta.parse(_julia_func_basename * "e1" * _julia_exact_func_expr1))
        eval(Meta.parse(_julia_func_basename * "d" * _julia_decimal_func_expr))
    end

    println("The following function$(count == 3 ? "s are" : " is") available ",
            "temporarily in the FiniteDifferenceFormula module. Usage:\n")
    println("  import FiniteDifferenceFormula as fd\n")
    println("  f, x, i, h = sin, 0:0.01:10, 501, 0.01")

    fmt = Printf.Format("%.8f")      # format of relative error

    # sine is taken as the example b/c sin^(n)(x) = sin(n π/2 + x), simply
    exact = sin(_data.n * pi /2 + 5) # x[501] = 5
    eval(Meta.parse("f, x, i, h = sin, 0:0.01:10, 501, 0.01"))

    _printexampleresult(count == 1 ? "" : "e", exact)
    if count == 3
        _printexampleresult("e1", exact)
        _printexampleresult("d", exact)
    end
    len = length("fd.$(_julia_func_basename)")
    if count == 3; len += 1; end
    print(" "^(len + 17) * "# cp:     ")
    @printf("%.16f\n", exact)

    # 250, a sepcial value for communication w/ another 'activatejuliafunction'
    if !(external_dataq && _formula_status == 250)
        println("\nCall fd.formula() to view the very definition.")
    end
    return
end  # activatejuliafunction

# calculate the coefficients of Taylor series of f(x[i + j]) about x[i]
# for teaching/learning!
function taylor(j::Int)
    global _max_num_of_taylor_terms
    return _taylor_coefs(j, _max_num_of_taylor_terms)
end  # taylor

# print readable Taylor series expansion of f(x[i + j]) about x[i]
# for teaching/learning!
function printtaylor(j::Int, num_of_nonzero_terms::Int = 10)
    if num_of_nonzero_terms < 0
        println("$num_of_nonzero_terms? It is expected to be an positive integer.\n")
        return;
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
    if num_of_nonzero_terms < 0
        println("$num_of_nonzero_terms? It is expected to be an positive integer.\n")
        return;
    end
    _print_taylor(coefs, num_of_nonzero_terms)
    return
end  # printtaylor

# for explorers/researchers. users never need to know its existence
function _set_default_max_num_of_taylor_terms(n::Int)
    global _max_num_of_taylor_terms
    if n < 0
        println("Wrong input, $n. A positive integer is expected.")
    elseif n < 10
        println("$n? It doesn't have to be so small. The default value ",
                "is still $_max_num_of_taylor_terms.")
    else
        _max_num_of_taylor_terms = round(Int, n)
    end
    return _max_num_of_taylor_terms
end  # _set_default_max_num_of_taylor_terms

function _set_default_max_num_of_taylor_terms()
    return _max_num_of_taylor_terms
end

end # module
