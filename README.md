# FiniteDifferenceFormula

This Julia package generates finite difference formulas for n-th order derivatives by
using Taylor series expansions of a function at evenly spaced points. It also gives the
truncation error of a formula in the big-O notation. You can use it to generate new
formulas in addition to verification of known formulas. See comments in the source code
for the algorithm.

It is for fun when you teach/learn numerical computing, especially the finite difference
method. You can explore the distribution, symmetry, and beauty in the coefficients in the
formulas and make some related conjectures and even give proofs to the conjectures.

There is a natural constraint, i.e., the number of points to be used in a formula can't
be too large. It's due to

1) the constraint of the largest integer computers can express and process.
2) that we want "exact" formulas, i.e., coefficients in the formulas are exact, otherwise,
   extra truncation errors occur already in the formulas.
3) that, to explore the properties of coefficients in the formulas, we had better have
   "exact" formulas.
4) that practical applications usually don't allow too many points to be used for
   computing derivatives of a function at a single point because of the constraints of
   computing power.

To run the code, you need the Julia programming language (https://julialang.org/), a
wonderful and amazing computing platform.

             *     *     *     *     *     *     *     *     *     *

The package exports three functions:

1. computecoefs
2. formula
3. set_decimal_places:          15 by default

function computecoefs(n::Int64, points::UnitRange{Int64}, printformulaq::Bool = false)
--------------------------------------------------------------------------------------
Input:

                   n: the n-th order derivative to be found
              points: in the format of a range, start:stop
       printformulaq: print the computed formula or not

Output:

The function returns a tuple, ([k[1], k[2], ..., k[stop-start+1]], m) where k[:] and m
are described belows. With the information, you may generate formulas for any
programming language of your choice.

The algorithm uses the linear combination of f(x[i+j]), j in start:stop, to eliminate
f(x[i]), f'(x[i]), ..., so that the first term of the Taylor series expansion of the
linear combination is f^(n)(x[i]):

    k[1]*f(x[i+start]) + k[2]*f(x[i+start+1]) + ... + k[stop-start+1]*f(x[i+stop])
        = m*f^(n)(x[i]) + ..., m > 0

It is this equation that gives the formula for computing f^(n)(x[i]) and the truncation
error in the big-O notation as well.

function formula()
------------------
The function presents

1) k[1]*f(x[i+start]) + k[2]*f(x[i+start+1]) + ... + k[stop-start+1]*f(x[i+stop])
       = m*f^(n)(x[i]) + ..., m > 0

2) The formula for finding f^(n)(x[i]), including estimation of accuracy in the big-O
   notation.

3) Julia function for f^(n)(x[i]).

             *     *     *     *     *     *     *     *     *     *

function set_decimal_places(n)
------------------------------
The function sets the decimal places to be n for printing Julia function for a formula.

Examples
========

using FiniteDifferenceFormula

set_decimal_places(6)        # n = 15 by default (without calling the function)

computecoefs(1, 0:2, true)   # find and print 3-point forward finite difference formula for f'(x[i])

computecoefs(2, -3:0, true)  # find and print 4-point backward finite difference formula for f''(x[i])

computecoefs(3, -9:9, true)  # find and print 19-point central finite difference formula for f'''(x[i])

computecoefs(4,-160:260)     # find 421-point finite difference formula for f''''(x[i])     # for fun
                             # it takes quite a while to finish the computation

formula()                    # print the formula computed last time you called computecoefs(...)
