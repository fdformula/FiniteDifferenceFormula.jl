#
# This Julia module, FiniteDifferenceFormula, generates n-point finite difference
# formulas for the 1st, 2nd, ..., derivatives by using Taylor series expansions of
# a function at evenly spaced points. It also allows users to examine and test
# formulas right away in current Julia REPL session. It surely helps when we
# teach/learn numerical computing, especially, the finite difference method.
#
# David Wang, dwang at liberty dot edu, on 12/20/2022
#
# Warning: users should not call/access a function/variable starting with "_".
#
module FiniteDifferenceFormula

using Printf

# v1.3.3 for parallel computing using multiple threads in _rref!
# Start Julia: julia -t auto
using Base.Threads

############################# EXPORTED FUNCTIONS ##############################
export compute, find, findforward, findbackward, formula, formulas
export truncationerror, verifyformula, activatejuliafunction
export decimalplaces, taylorcoefs, tcoefs, taylor

######################### BEGIN OF GLOBAL VARIABLES ###########################
_NUM_OF_EXTRA_TAYLOR_TERMS::Int     = 8       # for examining truncation error
_decimal_places::Int                = 16      # for generating Julia function(s)
                                              # call decimalplaces(n) to reset it

mutable struct _FDData
    n; points; k; m; coefs                    # on one line? separated by ;
end

_data                               = _FDData # share results between functions
_computedq::Bool                    = false   # make sure compute() is called first
_formula_status::Int                = 0       # a formula may not be available
                                              # values? see _test_formula_validity()

# a vector of the coefficients of Taylor series of the linear combination:
# k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ... + k[len]*f(x[i+points[len]])
_lcombination_coefs                 = Array{Any}

_range_inputq::Bool                 = false
_range_input::UnitRange{Int}        = 0:0     # compute receives a range? save it

_julia_exact_func_expr::String      = ""      # 1st exact Julia function for f^(n)(x[i])
_julia_exact_func_expr1::String     = ""      # 2nd exact Julia function for f^(n)(x[i])
_julia_decimal_func_expr::String    = ""      # decimal Julia function for f^(n)(x[i])
_julia_func_basename::String        = ""

_bigO_exp::Int                      = Threads.nthreads() # used temporarily
const _nthreads::Int                = _bigO_exp > 1 ? _bigO_exp - 1 : 1 # v.1.3.4

_bigO::String                       = ""      # truncation error of a formula
_bigO_exp::Int                      = -1      # the value of n as in O(h^n)
########################### END OF GLOBAL VARIABLES ###########################

# for future coders/maintainers of this package:
# to compute a new formula, this function must be called first.
function _reset() # 1.3.1, renamed from _initialization()
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

    if @isdefined(x) && length(x) > 201 #v1.3.2
        x = []
    end
end  # _reset

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
function _taylor_coefs(h, max_num_of_terms = 30)
    result = Matrix{Rational{BigInt}}(undef, 1, max_num_of_terms)
    factorial::BigInt = 1
    for n in 1 : max_num_of_terms
        N = n - 1          # order of derivatives in Taylor series
        if N > 0; factorial *= N; end
        result[n] = 1 // factorial * BigInt(h)^N
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
    println(" + ...")
    return
end  # _print_taylor

function _dashline(n = 105); return "-" ^ n; end

function _validate_input(n, points, printformulaq = false)
    if !isinteger(n) || n < 1
        println("Invalid order of derivatives, $n. A positive integer ",
                "is expected.")
        return []
    end
    n = round(Int, n)  # 4.0 --> 4

    len = length(points)
    if len == 0 ||  # v1.1.5, invalid input: 10:9
       (typeof(points) <: Tuple) ||
       !(typeof(points[1]) <: Int64)  # v1.3.1 from Integer
        println("Invalid input, $points. A list of integers like -1:2 or ",
                "[-1, 0, 1, 2] is expected.")
        return []
    end
    if typeof(printformulaq) != Bool
        println("Invalid input, $printformulaq. A value, false or true, is ",
                "expected.")
        return []
    end

    # v1.3.1, handling exceptions
    oldpoints = [] # define the variable
    try
        oldpoints = collect(points) # v1.2.8
        points = sort(unique(oldpoints))
    catch OutOfMemoryError
        println("Memory allocation error: _validate_input.")
        return []
    end
    len = length(points)
    if len < 2
        println("Invalid input, $points. A list of two or more different ",
                "points is expected.")
        return []
    end

    _reset()
    if oldpoints != points
        input_points = _format_of_points(points)
        print(_dashline(), "\nYour input is converted to ($n, $input_points")
        if printformulaq; print(", true"); end
        println(").\n", _dashline())
    else
        _format_of_points(points)   # no change; set _range_input[q]
    end

    return points
end  # _validate_input

"""
```compute```(n, points, printformulaq = false)

Compute a formula for the nth order derivative using the given points.
```
            n: the n-th order derivative to be found
       points: in the format of a range, start : stop, or a vector
printformulaq: print the computed formula or not
```
|   points     |   The points/nodes to be used                  |
|   ---------- | ---------------------------------------------- |
|    0:2       |   x[i], x[i+1], x[i+2]                         |
|   -2:2       |   x[i-2], x[i-1], x[i], x[i+1], x[i+2]         |
|   -3:2       |   x[i-3], x[i-2], x[i-1], x[i], x[i+1], x[i+2] |
|   [1, 1, -1] |   x[i-1], x[i+1]                               |
|   [1 0 1 -1] |   x[i-1], x[i], x[i+1]                         |

A vector can be like [1, 0, 2] or [1 0 2]. It will be rearranged so
that elements are ordered from lowest to highest with duplicate ones
removed.

Examples
====

```julia-repl
julia> import FiniteDifferenceFormula as fd
julia> fd.compute(2, [0 1 2 3])
julia> fd.compute(2, 0:3)
julia> fd.compute(2, [-5, -2, 1, 2, 4], true)
```
"""
function compute(n, points, printformulaq = false)
    points = _validate_input(n, points, printformulaq)
    if points == []; return nothing; end
    return _compute(n, points, printformulaq)
end  # compute

"""
```find```(n, points, printformulaq = false)

Compute a formula for the nth order derivative using the given points.

For the input, ```n``` and ```points``` (See [```compute```]), there may not be
formulas which use the two end points like -2 and 3 in -2 : 3 or [-2, -1, 0, 1,
2, 3]. In this case, ```find``` tries to find a formula by shrinking the
range to, first -1 : 3, then, -2 : 2, then -1 : 2, and so on, until a formula
is found or no formulas can be found at all.

See also [```compute```], [```findbackward```], and [```findforward```].

Examples
=====
```
import FiniteDifferenceFormula as fd
fd.find(2, -10:9)
```
"""
function find(n, points, printformulaq = false)
    global _range_input, _range_inputq
    points = _validate_input(n, points, printformulaq)
    if points == []; return nothing; end
    result = _compute(n, points, printformulaq)
    while result == nothing && length(points) > n + 1   # failed
        if _range_inputq; _range_input = points[2] : points[end]; end
        result = _compute(n, points[2 : end], printformulaq)

        if result == nothing
            if _range_inputq; _range_input = points[1] : points[end - 1]; end
            result = _compute(n, points[1 : end - 1], printformulaq)

            if result == nothing
                pop!(points); popfirst!(points)  # points = points[2 : end - 1]
                if _range_inputq; _range_input = points[1] : points[end]; end
            end
        end
    end
    return result
end

function _findforward(n, points, printformulaq = false, forwardq::Bool = true)
    global _range_input, _range_inputq
    points = _validate_input(n, points, printformulaq)
    if points == []; return nothing; end
    result = _compute(n, points, printformulaq)
    while result == nothing && length(points) > n + 1   # failed
        #points = forwardq ? points[2 : end] : points[1 : end - 1]
        if forwardq; popfirst!(points); else pop!(points); end
        if _range_inputq; _range_input = points[1] : points[end]; end
        result = _compute(n, points, printformulaq)
    end
    return result
end

"""
```findforward```(n, points, printformulaq = false)

Compute a formula for the nth order derivative using the given points.

For the input, ```n``` and ```points``` (See [```compute```]), there may not be
formulas which use the two end points like -2 and 3 in -2 : 3 or [-2, -1, 0,
1, 2, 3]. In this case, ```findforward``` tries to find a formula by
shrinking the range from the left endpoint to, first -1 : 3, then, 0 : 3, then
1 : 3, and so on, until a formula is found or no formulas can be found at all.

See also [```compute```], [```find```], and [```findbackward```].

Examples
=====
```
import FiniteDifferenceFormula as fd
fd.findforward(2, -10:9)
```
"""
function findforward(n, points, printformulaq = false)
    return _findforward(n, points, printformulaq, true)
end

"""
```findbackward```(n, points, printformulaq = false)

Compute a formula for the nth order derivative using the given points.

For the input, ```n``` and ```points``` (See [```compute```]), there may not be
formulas which use the two end points like -2 and 3 in -2 : 3 or [-2, -1, 0,
1, 2, 3]. In this case, ```findbackward``` tries to find a formula by
shrinking the range from the right endpoint to, first -2 : 2, then, -2 : 1, then
-2 : 0, and so on, until a formula is found or no formulas can be found at all.

See also [```compute```], [```find```], and [```findforward```].

Examples
=====
```
import FiniteDifferenceFormula as fd
fd.compute(3,-100:50)
# output: ***** Warning: 3, -100:22 might be your input for which a formula is found.
fd.findbackward(3,-99:50)
# output: (3, -99:23, ...)
```
"""
function findbackward(n, points, printformulaq = false)
    return _findforward(n, points, printformulaq, false)
end

# "evenly" assign at least smallest_chunk_size elements to each thread/task
#
function partition(chunks, smallest_chunk_size, nthreads = _nthreads)
    list = Vector()
    len = length(chunks)
    d, r = divrem(len, nthreads)

    if d < smallest_chunk_size
        nthreads = div(len, smallest_chunk_size)
        d, r = divrem(len, nthreads == 0 ? 1 : nthreads)
    end

    d -= 1
    start = chunks.start
    while start <= chunks.stop
        stop = start + d
        if r > 0
            stop += 1
            r -= 1
        end
        push!(list, start : stop)
        start = stop + 1
    end
    return list
end

# v1.2.9, reimplemented Gauss-Jordan Elimination and optimized on 10/4/2023 for
# the purpose of this project. at the end, # A is "virtually" an identity matrix
# (but not actually). the time performance has a 35% increase.
#
# input:  A and b as in Ax = b
# output: b is the solution x
#
# assume: A is square & invertible; A = [ 1 0 0 ... 0; x x x ...]
#
# v1.3.3, use Threads.@threads, Threads.@spawn, or Folds.map to improve time performance
# dramatically for a large matrix A.
#
# Parallelization through Threads.@threads is easier to implement, but proper use of Threads.@spawn
# is better. Both are better than Folds.map.
#
# v1.3.4 replaces the key @threads with @spawn
function _rref!(A::Matrix{Rational{BigInt}}, b::Matrix{Rational{BigInt}})
    nr, nc = size(A);
    i = 1
    while i < nr
        j = i + 1
        # make a[i, i] the pivotal entry
        if i != 1                    # A[1, 1] = 1 is already the pivotal entry
            # find the largest entry on A[i:end, i]
            (_, mi) = findmax(abs.(A[i : nr, i])) # It seems a little faster than Folds.findmax
            mi += i - 1                           # and code using @threads for (nr = 12800)!!!

            if mi != i               # interchange the two rows
                A[i, i : nc], A[mi, i : nc] = A[mi, i : nc], A[i, i : nc]
                b[i], b[mi] = b[mi], b[i]
            end
            A[i, j : nc] /= A[i, i]
            b[i] /= A[i, i]
            # A[i, i] = 1            # unnecessary
        end

        # v1.3.4, rewritten, 1.37X speedup when A is large
        chunks = partition(j : nc, 5)
        for r = j : nr
            Ari = A[r, i]
            if Ari != 0
                b[r] -= Ari * b[i]

                tasks = map(chunks) do chunk
                    Threads.@spawn begin
                        local lk = ReentrantLock()
                        lock(lk)              # though not needed theretically, it does matter!
                        try
                            A[r, chunk] -= Ari * A[i, chunk]  # for k in chunk; A[r, k] -= Ari * A[i, k]; end
                        finally
                            unlock(lk)
                        end
                    end
                end
                wait.(tasks)
                # A[r, i] = 0        # unnecessary
            end
        end
        i = j
    end
    b[nr] /= A[nr, nr]
    # A[nr, nr] = 1                  # unnecessary

    # eliminate entries above A[i, i]
    for i in nr : -1 : 2      # can't be parallelized !
        chunk = 1 : i - 1
        b[chunk] -= A[chunk, i] * b[i]
    end
end # _rref!

#
# Algorithm
# ---------
# General generator of finite difference formulas for the n-th order derivatives
# of f(x) at x[i] using points in a sorted list, e.g., [-2 -1 0 1 3 5 6 7].
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
# Values of the coefficients k[:] and m will be determined, and the first few
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
    global _NUM_OF_EXTRA_TAYLOR_TERMS

    # for teaching's purpose, we don't do so
    # if length(points) <= n
    #     pts = _range_inputq ? "$(_range_input)" : "$(points')"
    #     println("$pts is invalid because at least $(n + 1) points are ",
    #             "needed for the $(_nth(n)) derivative.")
    #     return nothing
    # end

    # setup a linear system Ax = k first
    len = length(points)
    max_num_of_terms = max(len, n) + _NUM_OF_EXTRA_TAYLOR_TERMS

    # setup the coefficients of Taylor series expansions of f(x) at each of the
    # involved points
    # v1.3.1, handling exceptions
    coefs = []  # define the variable
    try
        coefs = Array{Any}(undef, len)
    catch OutOfMemoryError
        println("Memory allocation error: _compute #1.")
        _reset()
        return nothing
    end

    ##f0(x) = _taylor_coefs(x, max_num_of_terms) # Folds v0.2.8 failed
    ##coefs = Folds.map(f0, points)
    @threads for i in 1 : len
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
    # v1.3.1, handling exceptions
    A = [] # define the variables
    k = []
    try
        A = Matrix{Rational{BigInt}}(undef, len, len)
        k = zeros(Rational{BigInt}, 1, len) # v1.3.3 for output: no k' is needed
    catch OutOfMemoryError
        println("Memory allocation error: _compute #2.")
        _reset()
        return nothing
    end
    A[1, 2 : len] .= 0                   # setup an equation so that k[1] = 1
    A[1, 1]        = 1
    k[1]           = 1

    row = 2
    for order in 0 : len - 1  # correspond to f(x[i]), f'(x[i]), f''(x[i]), ...
        if order == n; continue; end     # skip f^(n)(x[i])

        # eliminating f^(order)(x[i])
        # A[:,j] stores coefs of
        #  k[1]*coefs[1][j] + k[2]*coefs[2][j] + ... + k[len]*coefs[len][j] = 0
        tmp = order + 1
        @threads for j = 1 : len
            A[row, j] = coefs[j][tmp]
        end
        if row == len; break; end
        row += 1
    end

    # The homogeneous linear system (2) has no nontrivial solution or has
    # infinitely many nontrivial solutions. It is understandable that it may
    # not have a nontrivial solution. But why infinitely many nontrivial
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

    # solve Ax = k for x, i.e., k[:]
    _rref!(A, k); A = []                 # output: k is the solution
    tmp = gcd(k)                         # change each element to an integer
    @threads for i in 1 : len
        k[i] //= tmp
    end
    ## k //= gcd(k)
    ## #the following code does the same  # for translating to other languages
    ## for i in 1 : len
    ##    if k[i] != round(BigInt, k[i]); k *= denominator(k[i]); end
    ## end

    # Taylor series expansion of the linear combination
    # k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ... + k[len]*f(x[i+points[len]])
    _lcombination_coefs = k[1] * coefs[1]  # let Julia determine the type
    for i in 2 : len
        if k[i] == 0; continue; end
        @threads for j in 1 : max_num_of_terms
            _lcombination_coefs[j] += k[i] * coefs[i][j]
        end
    end
    ##function f2(x, y)
    ##    f(z) = x * z
    ##    return Folds.map(f, y)
    ##end
    ##_lcombination_coefs = Folds.sum(Folds.map(f2, k[1:len], coefs[1:len]))

    m = _lcombination_coefs[n + 1]
    # let _test_formula_validity() take care if m = 0 or m is the (n+1)th term.

    # "normalize" k[:] and m so that m is a positive integer
    if m < 0
        @threads for i in 1 : len
            k[i] *= -1
        end
        @threads for i in 1 : max_num_of_terms
            _lcombination_coefs[i] *= -1
        end
        m = _lcombination_coefs[n + 1]
    end

    x = round(BigInt, m)
    if x == m; m = x; end  # already integer. don't show like 5//1

    # save the results in a global variable for other functions
    global _data = _FDData(n, points, k, m, coefs)
    global _computedq = true

    _test_formula_validity()

    if printformulaq; formula(); end

    if _formula_status >= 0
        return (n, _range_inputq ? _range_input : points, round.(BigInt, k), m)
    else
        _reset()
        return nothing
    end
end  # _compute

"""
```loadcomputingresults```(results)

Input: 'results' is a tuple, (n, points, k[:], m). See compute(...).

Load computing results from the output of compute(...). After this
command, formula(), activatepythonfunction(), truncationerror(0, etc.,
are available. It allows users to work on saved computing results (say,
in a textfile). For example, when it takes hours to compute/find a
formula, users may run commands like the following one from OS terminal

julia -e "import FiniteDifferenceFormula as fd; println(fd.compute(1, range(-200, 201)))" > data.txt

and then mannually load data from data.txt to this function later.

See also [compute].
"""
function loadcomputingresults(results)
    global _lcombination_coefs, _data, _computedq, _NUM_OF_EXTRA_TAYLOR_TERMS

    if !(typeof(results) <: Tuple) || length(results) != 4
        println("Invalid input. A tuple of the form (n, points, k, m) is expected.")
        println("It is the output of compute(...).")
        return
    end

    _reset()
    n, points, k, m, = results
    points = collect(points)
    len = length(points)
    if len == length(k) + 1
        println("\nYour input might be from Python FiniteDifferenceFormula.")
        pop!(points)
        len -= 1
    end
    _format_of_points(points)
    k = Rational.(k)

    max_num_of_terms = max(len, n) + _NUM_OF_EXTRA_TAYLOR_TERMS

    # setup the coefficients of Taylor series expansions of f(x) at each of
    # the involved points
    coefs = Array{Any}(undef, len)
    for i in 1 : len
        coefs[i] = _taylor_coefs(points[i], max_num_of_terms)
    end

    # Taylor series of the linear combination
    # k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ...
    _lcombination_coefs = k[1] * coefs[1]  # let Julia determine the type
    for i in 2 : len
        if k[i] == 0; continue; end
        _lcombination_coefs += k[i] * coefs[i]
    end

    _data = _FDData(n, points, k, m, coefs)
    _computedq = true

    _test_formula_validity(true)
    return
end  # of loadcomputingresults

# return a string of the linear combination
# k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ... + k[len]*f(x[i+points[len]])
#
# To test in Julia if a formula is valid, we must convert the data to BigInt or
# BigBloat because, say, for compute(2, -12:12), without the conversion,
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
                # times = "*"   # -2x == -2 * x in Julia;  an error in other languages
                c2s = "1"
            elseif c2s == "-"
                # times = "*"   # same as above
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
#        coefficients = 0, symmetry of coefficients about x[i] in a central formula
# -100 - no formula can't be found
# -200 - same as -100, but do not try to find a formula if
#        activatejuliafunction(n, point, k, m) fails
#  250 - same as -100, but used for communication btw 2 'activatejuliafunction's
# > 0  - mathematically, the formula is still valid but may not satisfy some
#        commonly seen "rules" such as the sum of coefficients = 0 and symmetry
#        of coefficients about x[i]
#
# return m as in equation (1) for 'activatejuliafunction'
function _test_formula_validity(verifyingq::Bool = false)
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
    k = _data.k           # just convenient names for the data; not copies!
    coefs = _data.coefs   #
    points = _data.points #
    len = length(points)

    input_points = _range_inputq ? _range_input : points'

    # Is there any equation in system (2) that is not satisfied?
    has_solutionq = true
    global _formula_status = 0
    for i = 1 : n
        if _lcombination_coefs[i] != 0
            x = round(Integer, _lcombination_coefs[i])   # v1.1.7
            if x != _lcombination_coefs[i]
                x = Float64(_lcombination_coefs[i])
            end
            fnxi = "f" * (i <= 4 ? ("'" ^ (i - 1)) : "^($(i - 1))") * "(x[i])"

            print("***** Error: $n, $input_points")
            if verifyingq; print(", $k"); end
            println(" : k[1]*coefs[1][$i]",
                    " + k[2]*coefs[2][$i] + ... + k[$len]*coefs[$len][$i] ",
                    "= $x != 0, i.e., $fnxi can't be eliminated as indicated ",
                    "in the following computing result:")
            println(_dashline())
            print(_lcombination_expr(_data), " =\n    ")
            _print_taylor(_lcombination_coefs, 5)   # print at most 5 nonzero terms
            println(_dashline())
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
            if _range_inputq
                s = _range_input = points[start] : points[stop]  # v1.1.5
            else
                s = points[start : stop]
            end
            println("***** Warning: $n, $s might be your input for which a ",
                    "formula is found.\n")
            formula_for_inputq = false
        end
    end

    if !has_solutionq
        if len <= n
            println("***** Error: $n, $input_points : Invalid input. ",
                    "At least $(n + 1) points are needed for the $(_nth(n)) ",
                    "derivative.")
            _formula_status = -200
            return m
        end
        print("\n***** Error: $n, $input_points")
        if verifyingq; print(", $k"); end
        println(": can't find a formula.\n")
        _formula_status = -100
        return m
    end

    if sum(k) != 0   # sum of coefficients must be 0
        println("***** Warning: $n, $input_points : sum(k[:]) != 0")
        _formula_status += 1
    end

    # are coefficients of central formulas symmetrical about x[i]?
    if formula_for_inputq && _range_inputq && abs(_range_input.start) == _range_input.stop
        j = kLen = length(k)
        for i in 1 : round(Int64, kLen / 2)
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

    n = _num_of_used_points()

    global _julia_func_basename = "fd$(_nth(data.n))deriv$(n)pt$(s)"
    fexpr  = "(f, x, i, h) = "
    if julia_REPL_funcq; fexpr *= "Float64( "; end   # convert the final result
    fexpr *= "( "
    fexpr *= _lcombination_expr(data, decimalq, julia_REPL_funcq)
    fexpr *= " ) / " * _denominator_expr(data, julia_REPL_funcq)
    if julia_REPL_funcq; fexpr *= " )"; end
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
"""
```formula```()

Generate and list:

1. ```k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ... + k[len]*f(x[i+points[len]])
       = m*f^(n)(x[i]) + ..., m > 0```

1. The formula for f^(n)(x[i]), including estimation of accuracy in the
   big-O notation.

1. Julia function(s) for f^(n)(x[i]).

Calling ```compute(n, points, true)``` is the same as calling
```compute(n, points)``` and then ```formula()```.

Even if no formula can be found, it still lists the computing results from which
we can see why. For example, after ```compute(2,1:2)```, try ```formula()```.
"""
function formula()
    global _data, _computedq, _formula_status, _lcombination_coefs
    global _range_inputq, _range_input, _bigO, _julia_func_basename

    if !_computedq
        println("Please call 'compute', 'find', 'findbackward', or 'findforward' first!")
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
    print(_lcombination_expr(_data), " =\n    ")
    _print_taylor(_lcombination_coefs, 5)   # print at most 5 nonzero terms

    if _formula_status > 0
        println("\nThe exact formula:\n")
        _print_bigo_formula(_data, _bigO)
        data1 = _FDData
        if _data.m != 1                 # print in another format
            # _data.k // _data.m ==> rationalize.(convert.(Float64, _data.k // _data.m))
            # no differnce
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
    else
        println("\nNo formula can be found.") # v1.1.7
    end

    return
end  # formula

# print _bigO, the estimation of truncation error in the big-O notation
#
# output:
#    (-1, "")      - There is no valid formula
#    (n, "O(h^n)") - There is a valid formula
"""
```truncationerror```()

Show the truncation error of the last computed formula in the big_O notation.

Examples
====
```
import FiniteDifferenceFormula as fd
fd.compute(2,-3:3)
fd.truncationerror()
fd.find(3,[-2, 1, 2, 5, 7, 15])
fd.truncationerror()
```
"""
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
    println("Please call 'compute', 'find', 'findbackward', or 'findforward' first!")
    return (-1, "")
end  # truncationerror

######################## for teaching/learning/exploring #######################

# show current decimal places
"""
```decimalplaces```() or ```decimalplaces```(n)

Show present decimal places for generating Julia function(s) of computed
formulas, if no argument is passed to ```decimalplaces```. Otherwise, set
the decimal places to n.

Examples
====
```
import FiniteDifferenceFormula as fd
fd.compute(2,-3:3)
fd.formula()  # by default, use 16 decimal places to generate a Julia function
fd.decimalplaces(4)
fd.formula()  # now, use 4 decimal places to generate a Julia function
```
"""
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
            println("You may start your work by calling 'compute', 'find', 'findbackward', or 'findforward'.")
        end
    else
        println("decimalplaces(n): a nonnegative integer is expected.")
    end

    return _decimal_places
end  # decimalplaces

function _format_of_points(points)
    global _range_inputq, _range_input
    if length(points) == length(points[1] : points[end])
        _range_inputq = true
        _range_input  = points[1] : points[end]
        return _range_input
    end
    return points
end  # _format_of_points

function _printexampleresult(suffix, exact)
    global _julia_func_basename
    #1.3.1, handling exceptions
    apprx = 0.0    # define the variable
    try
        apprx = eval(Meta.parse("$(_julia_func_basename)$(suffix)(f, x, i, h)"))
    catch BoundsError
        println("BoundsError: _printexampleresult.")
    return -1  # failure
    end
    relerr = abs((apprx - exact) / exact) * 100
    print("  fd.$(_julia_func_basename)$(suffix)(f, x, i, h)  ",
          suffix == "e1" ? "" : " ", "# result: ")
    @printf("%.16f, ", apprx)
    print("relative error = "); @printf("%.8f", relerr); println("%")
    return 0       # success
end  # _printexampleresult

# the name is self-explanatory. it is exactly the same as the function
# activatejuliafunction(n, points, k, m)
"""
```verifyformula```(n, points, k, m = 1) or
```activatejuliafunction```(n, points, k, m = 1)

Verify if a formula is valid. If it is valid, generate and activate its Julia
function(s). If not, try to find a formula for the derivative using the points.

```
     n: the n-th order derivative
points: in the format of a range, start : stop, or a vector
     k: a list of the coefficients in a formula
     m: the coefficient in the denominator of a formula
```

Examples
====
```julia-repl
import FiniteDifferenceFormula as fd
fd.verifyformula(1,[-1,2],[-3,4],5)  # f'(x[i]) = (-3f(x[i-1])+4f(x[i+2]))/(5h)?
fd.verifyformula(2, -3:3, [2,-27,270,-490,270,-27,2], 18)
fd.activatejuliafunction(2, -3:3, [1/90,-3/20,3/2,-49/18,3/2,-3/20,1/90])
fd.verifyformula(2, -3:3, [1//90,-3//20,3//2,-49//18,3//2,-3//20,1//90])
```
"""
function verifyformula(n, points, k, m = 1)
    return activatejuliafunction(n, points, k, m)
end

# if you have data from textbooks or other sources, you may use this function
# to verify if it is right, activate related Julia function(s), evaluate and
# see the computiong results.
"""
```verifyformula(n, points, k, m = 1)``` or
```activatejuliafunction(n, points, k, m = 1)```

Verify if a formula is valid. If it is valid, generate and activate its Julia
function(s). If not, try to find a formula for the derivative using the points.

```
     n: the n-th order derivative
points: in the format of a range, start : stop, or a list
     k: a list of the coefficients in a formula
     m: the coefficient in the denominator of a formula
```

Examples
====
```julia-repl
import FiniteDifferenceFormula as fd
fd.verifyformula(1,[-1,2],[-3,4],5)  # f'(x[i]) = (-3f(x[i-1])+4f(x[i+2]))/(5h)?
fd.verifyformula(2, -3:3, [2,-27,270,-490,270,-27,2], 18)
fd.activatejuliafunction(2, -3:3, [1/90,-3/20,3/2,-49/18,3/2,-3/20,1/90])
fd.verifyformula(2, -3:3, [1//90,-3//20,3//2,-49//18,3//2,-3//20,1//90])
```
"""
function activatejuliafunction(n, points, k, m = 1)
    global _NUM_OF_EXTRA_TAYLOR_TERMS
    if !isinteger(n) || n <= 0
        println("Error: invalid first argument, $n. A positive integer is ",
                "expected.")
        return nothing
    end
    n = round(Int, n)  # 4.0 --> 4
    if length(points) == 1
        println("Error: invalid input, points = $points. A list of two or ",
                "more points is expected.")
        return nothing
    end
    if !(typeof(points[1]) <: Integer)
        println("Error: invalid input, $points. Integers are expected.")
        return nothing
    end

    # tuple input can cause trouble - Julia 1.8.5: x = (1, 2//3), y = collect(x)
    # typeof(y[1]) != typeof(y[2]) !!!
    if typeof(points) <: Tuple
        println("Invalid input, $points. A list like [1, 2, ...] is expected")
        return nothing
    end
    if typeof(k) <: Tuple
        println("Invalid input, $k. A list like [1, 2, ...] is expected")
        return nothing
    end

    len = length(points)
    if len != length(unique(collect(points))) # v1.2.7, removed sort(...)
        println("Error: Invalid input - $points contains duplicate points.")
        return nothing
    end
    points = collect(points)

    # don't do so for teaching
    #if n <= len
    #    println("Error: at least $(n+1) points are needed for the $(_nth(n))",
    #            " derivative.")
    #    return nothing
    #end
    if len != length(k)
        println("Error: The number of points != the number of coefficients.");
        return nothing
    end

    _reset()      # needed b/c it's like computing a new formula
    input_points = _format_of_points(points)

    rewrittenq = false
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
        return nothing
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
    max_num_of_terms = max(len, n) + _NUM_OF_EXTRA_TAYLOR_TERMS
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
    M = _test_formula_validity(true)
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

    if find_oneq
        # v1.2.7
        println(">>>>> Your input doesn't define a valid formula, but it is ",
                "still activated for your examination.\n")
        _formula_status = 250      # force to activate Julia function
        # 250, reserved for communication w/ another 'activatejuliafunction'
    else # v1.2.7
        println(">>>>> Your input defines a valid formula.\n")
    end
    global _computedq   = true     # assume a formula has been computed

    # force to activate even for invalid input formula
    b200::Bool = _formula_status == -200
    if b200; _formula_status = 100; end;   # assume it is valid
    activatejuliafunction(true)
    if b200; _formula_status = -200; end;  # change it back

    if find_oneq                   # use the basic input to generate a formula
        _formula_status = -100
        println(_dashline(), "\nFinding a formula using the points....\n")
        result = _compute(n, points)
        if _formula_status >= 0
            println("Call fd.formula() to view the results and fd.activatejulia",
                    "function() to activate the new Julia function(s).")
            return result
        end
    end
    return nothing
end  # activatejuliafunction

# activate function(s) for the newly computed finite difference formula,
# allowing immediate evaluation of the formula in Julia REPL
function activatejuliafunction(external_dataq = false)
    global _computedq, _formula_status, _data
    global _julia_func_basename, _julia_exact_func_expr
    global _julia_exact_func_expr1, _julia_decimal_func_expr

    if !(external_dataq || _computedq)
        println("Please call 'compute', 'find', 'findbackward', or 'findforward' first!")
        return nothing
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
        return nothing
    end

    eval(Meta.parse(_julia_func_basename * (count == 1 ? "" : "e") * _julia_exact_func_expr))
    if count == 3
        eval(Meta.parse(_julia_func_basename * "e1" * _julia_exact_func_expr1))
        eval(Meta.parse(_julia_func_basename * "d" * _julia_decimal_func_expr))
    end

    println("The following function$(count == 3 ? "s are" : " is") available ",
            "temporarily in the FiniteDifferenceFormula module. Usage:\n")
    println("  import FiniteDifferenceFormula as fd\n")
    # v1.3.0, data points are determined according to input points rather than
    # f, x, i, h = sin, 0:0.01:10, 501, 0.01
    h = 0.01
    center = max(abs(_data.points[1]), abs(_data.points[end]))
    if center < 99; center = 99; end
    tmp = (center * 2 + 1) * h
    center += 1
    stop = ceil(Int, tmp)
    example = "f, x, i, h = sin, 0:$h:$stop, $center, $h"
    # v1.3.1, handling exceptions
    # x = [] # can't define the variable (Julia 1.9.3) ! eval(Meta.parse(str)) must be very special
    try
        eval(Meta.parse(example))
    catch OutOfMemoryError
        x = []
        println("Memory allocation error: activatejuliafunction #1.")
        return nothing
    end

    #v1.3.2
    exact = 0.0 # define a local variable
    try
        # xi = (center - 1) * h   # or x[center]
        println("  $example  # xi = ", Printf.format(Printf.Format("%.2f"), x[center]))
        # sine is taken as the example b/c sin^(n)(x) = sin(n π/2 + x), simply
        exact = sin(_data.n * pi /2 + x[center])
    catch BoundsError
        x = []
        println("Memory allocation error: activatejuliafunction #2.")
        return nothing
    end

    if _printexampleresult(count == 1 ? "" : "e", exact) == 0 && count == 3
        if _printexampleresult("e1", exact) == 0
            _printexampleresult("d", exact)
        end
    end
    len = length("fd.$(_julia_func_basename)")
    if count == 3; len += 1; end
    print(" "^(len + 17) * "# cp:     ")
    @printf("%.16f\n", exact)

    # 250, a sepcial value for communication w/ another 'activatejuliafunction'
    if !(external_dataq && _formula_status == 250)
        println("\nCall fd.formula() to view the very definition.")
    end
    return nothing
end  # activatejuliafunction

# calculate the coefficients of Taylor series of f(x[i + j]) about x[i]
# for teaching/learning!
"""
```taylorcoefs```(j, n = 10))

Compute and return coefficients of the first n terms of the Taylor series of
f(x[i + j]) = f(x[i] + jh) about x[i], where h is the increment in x.

Examples
====
```
import FiniteDifferenceFormula as fd
fd.taylorcoefs(-2)
fd.taylorcoefs(5, 4)
```
"""
function taylorcoefs(j::Int, n::Int = 10)
    if n < 1
        println("n = $n? It is expected to be an positive integer.")
        return nothing
    end
    return _taylor_coefs(j, n)
end  # taylorcoefs

"""
```tcoefs```(j, n = 10))

Same as taylorcoefs(j, n).
"""
function tcoefs(j::Int, n::Int = 10)
    return taylorcoefs(j, n)
end  # tcoefs

# print readable Taylor series expansion of f(x[i + j]) about x[i]
# for teaching/learning!
"""
```taylor```()
  - Print the first few nonzero terms of the Taylor series of the linear
    combination k[0]f(x[i+j0]) + k[1]f(x[i+j1]) + ... for the newly
    computed formula (even if failed).

```taylor```(j, n = 10)
  - Print the 1st n terms of Taylor series of f(x[i+j]) about x[i].

```taylor```(coefs, n = 10), or
  - Print the 1st n terms of Taylor series with coefficients in 'coefs'

```taylor```(points, k, n::Int = 10)
  - Prints the 1st n nonzero terms of the Taylor series of the linear
    combination:  k[0]f(x[i+points[0]]) + k[1]f(x[i+points[1]]) + ...

The last two provide also another way to verify if a formula is mathematically
valid or not.

See also [```verifyformula```], [```activatejuliafunction```], and
[```taylorcoefs```].

Examples
====
```
import FiniteDifferenceFormula as fd
fd.compute(1, [0, 1, 5, 8])
fd.taylor()

fd.taylor(2)

n = 50
coefs = 2*fd.tcoefs(0, n) - 6*fd.tcoefs(1, n) + 4*fd.tcoefs(2, n)
fd.taylor(coefs, n) # this n can be any positive integer

fd.taylor(-fd.tcoefs(0) + 3*fd.tcoefs(1) - 3*fd.tcoefs(2) + fd.tcoefs(3))
fd.taylor(0:3, [-1, 3, -3, 1], 6)
```
"""
function taylor()   #v0.6.4
    global _data, _computedq, _lcombination_coefs
    # print the Taylor series of the linear combination of
    # k[1]f(x[i+points[1]]) + k[2]f(x[i+points[2]]) + ...
    if _computedq
        print(_lcombination_expr(_data), " =\n    ",)
        _print_taylor(_lcombination_coefs, 5)
    else
        println("Please call 'compute', 'find', 'findbackward', or",
                "'findforward' first!")
    end
    return
end  # taylor()

function taylor(j::Int, n::Int = 10)
    if n < 1
        println("n = $n? It is expected to be an positive integer.")
        return
    end
    coefs = taylorcoefs(j, n)
    print("f(x[i" * (j == 0 ? "" : (j > 0 ? "+$j" : "$j")) * "]) = ")
    _print_taylor(coefs, n)
    return
end  # taylor

# print readable Taylor series of a function/expression about x[i]. e.g.,
# fd.taylor(2*fd.taylorcoefs(0) - 5*fd.taylorcoefs(1) + 4*fd.taylorcoefs(2))
function taylor(coefs, n::Int = 10)
    if n < 1
        println("n = $n? It is expected to be an positive integer.")
        return
    end
    coefs = collect(coefs)
    _print_taylor(coefs, n)
    return
end  # taylor

# input: points and k[:] are as in the linear combination:
# k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ...
function taylor(points, k, n::Int = 10)
    global _NUM_OF_EXTRA_TAYLOR_TERMS
    if n < 1
        println("n = $n? It is expected to be an positive integer.")
        return
    end
    oldpoints = collect(points)
    points = sort(unique(oldpoints))
    len = length(points)
    if oldpoints != points #v1.2.8
        println(_dashline())
        println("Your input: points = $points.")
        println(_dashline())
    end
    oldpoints = []
    if len != length(k)
       println("Error: invalid input. The sizes of points and k are not the same.")
       return
    end
    k = collect(k)
    max_num_of_terms = max(n, len, 30) + _NUM_OF_EXTRA_TAYLOR_TERMS
    coefs = k[1] * _taylor_coefs(points[1], max_num_of_terms)
    for i in 2 : len
        coefs += k[i] * _taylor_coefs(points[i], max_num_of_terms)
    end
    _print_taylor(coefs, n)
    return
end  # taylor

# return the number of points actually used in a formula
function _num_of_used_points()
    global _data
    n = 0
    for i in eachindex(_data.points)
        if _data.k[i] == 0; continue; end
        n += 1
    end
    return n
end  # _num_of_used_points

"""
```formulas```(orders = 1:3, min_num_of_points = 2, max_num_of_points = 5)

By default, the function prints all forward, backward, and central finite
difference formulas for the 1st, 2nd, and 3rd derivatives, using 2 to 5 points.

Examples
====

```julia-repl
# The following examples show all forward, backward, and central finite
# difference formulas for the specified derivatives, using 4 to 11 points.
julia> import FiniteDifferenceFormula as fd
julia> fd.formulas(2:5, 4, 11)       # the 2nd, 3rd, .., 5th derivatives
juliq> fd.formulas([2, 4, 7], 4, 11) # the 2nd, 4th, and 7th derivatives
julia> fd.formulas(3, 4, 11)         # the 3rd derivative
```
"""
function formulas(orders = 1:3,
                  min_num_of_points::Int = 2,
                  max_num_of_points::Int = 5)
    global _data, _bigO
    if  !(typeof((collect(orders))[1]) <: Integer)
        println("Error: Invalid input, orders = $orders. ",
                "It must be a positive integer or a list of positive integers.")
        return
    end
    if min_num_of_points < 2
        println("Error: Invalid input, min_num_of_points = $min_num_of_points. ",
                "It must be at least 2.")
        return
    end
    if  max_num_of_points < min_num_of_points
        println("Error: Invalid input, max_num_of_points = $max_num_of_points. ",
                "It must be at least $min_num_of_points.")
        return
    end
    for i in eachindex(orders)
        if orders[i] < 1
            println("Error: Invalid input, orders = $orders. ",
                    "Positive integers are expected.")
            return
        end
    end

    oldorders = collect(orders)
    orders = sort(unique(oldorders))
    if oldorders != orders
        println(_dashline())
        println("Your input: formulas($orders, $min_num_of_points, ",
                "$max_num_of_points)")
        println(_dashline())
    end
    oldorders = []

    for n in orders
        # forward schemes
        start = max(n + 1, min_num_of_points)
        for num_of_points in start : max_num_of_points
            compute(n, 0 : num_of_points - 1)
            if _formula_status > 0
                println(_num_of_used_points(),
                        "-point forward finite difference formula:")
                _print_bigo_formula(_data, _bigO)
            end
        end

        # backward schemes
        for num_of_points in start : max_num_of_points
            compute(n, 1 - num_of_points : 0)
            if _formula_status > 0
                println(_num_of_used_points(),
                        "-point backward finite difference formula:")
                _print_bigo_formula(_data, _bigO)
            end
        end

        # central schemes
        start = floor(Int, max(n, min_num_of_points) / 2)
        stop  = ceil(Int, max_num_of_points / 2)
        for num_of_points in start : stop
            len = 2 * num_of_points + 1
            if n >= len; continue; end
            compute(n, -num_of_points : num_of_points)
            if _formula_status > 0
                x = _num_of_used_points()
                if x <= max_num_of_points
                    println(x, "-point central finite difference formula:")
                    _print_bigo_formula(_data, _bigO)
                end
            #else
            #    _reset() # v1.2.6
            end
        end
    end

    return
end  # formulas

end # module
