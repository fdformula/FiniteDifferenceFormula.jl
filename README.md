# FiniteDifferenceFormula

This Julia package provides a general finite difference formula generator and a tool
for teaching/learning the finite difference method. It generates finite difference
formulas for derivatives of various orders by using Taylor series expansions of a
function at evenly spaced points. It also gives the truncation error of a formula
in the big-O notation. You can use it to generate new formulas in addition to
verification of known ones. See documentation in the source code for the algorithm.

You may play with this package when teaching/learning numerical computing, especially
the finite difference method, and explore the distribution, symmetry, and beauty in
the coefficients of the formulas. By changing decimal places, we can also see how
rounding errors affect a result.

Beware, though formulas are mathematically correct, they may not be numerically useful.
This is true especially when we derive formulas for a derivative of higher order. For
example, run computecoefs(9,-5:5), provided by this package, to generate a 10-point
central formula for the 9-th derivative. The formula is mathematically correct, but it
can hardly be put into use for numerical computing without, if possible, rewriting it
in a special way. Similarly, the more points are used, the more precise a formula
is mathematically. However, due to rounding errors, this may not be true numerically.

To run the code, you need the Julia programming language (https://julialang.org/), a
wonderful and amazing computing platform.

## How to install this package

In Julia REPL,

1. using Pkg
1. Pkg.add("FiniteDifferenceFormula")

## The package exports seven functions

- computecoefs
- decimalplaces
- formula
- activatejuliafunction
- truncationerror
- taylor
- printtaylor

### function computecoefs(n, points, printformulaq = false)

#### Input

```
            n: to find a formula for the n-th order derivative
       points: it can be in the format of a range, start : stop, or a vector
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

#### Output

The function returns a tuple, ([k[1], k[2], ..., k[stop-start+1]], m) where k[:] and m
are described below. With the information, you may generate formulas for any
programming language of your choice.

The algorithm uses the linear combination of f(x[i+j]), j ∈ points, a given list of points,
to eliminate f(x[i]), f'(x[i]), f''(x[i]), ..., so that the first nonzero term of the Taylor
series of the linear combination is f^(n)(x[i]):

```Julia
    k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ... + k[len]*f(x[i+points[len]])
        = m*f^(n)(x[i]) + ..., m > 0
```

where len = length(points). It is this equation that gives the formula for computing f^(n)(x[i])
and the truncation error in the big-O notation as well.

### function decimalplaces( ) or decimalplaces(n)

Without argument, the function returns current decimal places. With argument n, it sets the
decimal places to be n for generating Julia function for a formula. Without/before calling
the function, 16 decimal places are used by default.

### function formula( )

The function generates and lists

1. k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ... + k[len]*f(x[i+points[len]])
       = m*f^(n)(x[i]) + ..., m > 0

1. The formula for finding f^(n)(x[i]), including estimation of accuracy in the big-O
   notation.

1. Julia function for f^(n)(x[i]).

### function activatejuliafunction( )

Call this function to activate the Julia function(s) based on the newly computed
finite difference formula. E.g.,

```Julia
f1stderiv2ptcentrale(f, x, i, h)  = ( -f(x[i-1]) + f(x[i+1]) ) / (2 * h)
f1stderiv2ptcentrale1(f, x, i, h) = ( -1/2 * f(x[i-1]) + 1/2 * f(x[i+1]) ) / h
f1stderiv2ptcentrald(f, x, i, h)  = ( -0.5000 * f(x[i-1]) + 0.5000 * f(x[i+1]) ) / h
```
The suffixes 'e' and 'd' stand for 'exact' and 'decimal', respectively.

After activating the function(s), you can evaluate right away in the present Julia REPL
session.

```Julia
FiniteDifferenceFormula.f1stderiv2ptcentrale(sin, 0:0.01:pi, 3, 0.01)
```

### function truncationerror( )

The function shows the truncation error of the newly computed finite difference formula.

-----

The following functions are provided for teaching/learning the finite difference method.

### function taylor(j)

The function returns the coefficients of the first 30 terms of the Taylor series of f(x[i+j])
about x[i].

### function printtaylor(j, n = 10)

The function prints the first n terms of the Taylor series of f(x[i+j]) about x[i].

### function printtaylor(coefficients_of_taylor_series, n = 10)

The function prints the first n nonzero terms of a Taylor series of which the coefficents are
provided.

## Examples

```Julia
import FiniteDifferenceFormula as fd

fd.decimalplaces(6)               # use 6 decimal places to generate Julia functions of computed formulas

fd.computecoefs(1, 0:2, true)     # find, generate, and print "3"-point forward formula for f'(x[i])

fd.computecoefs(2, -3:0, true)    # find, generate, and print "4"-point backward formula for f''(x[i])

fd.computecoefs(3, -9:9)          # find "19"-point central formula for f'''(x[i])

fd.computecoefs(2, [-3 -2 1 2 7]) # find formula for f''(x[i]) using points x[i+j], j = -3, -2, 1, 2, and 7

fd.computecoefs(1,-230:230)       # find "461"-point central formula for f'(x[i]). does it exist? run the code!

fd.formula()                      # generate and print the formula computed last time you called computecoefs(...)

fd.truncationerror()              # print the truncation error of the newly computed formula

fd.printtaylor(-2, 5)             # print the first 5 terms of the Taylor series of f(x[i-2]) about x[i]

coefs = 2*fd.taylor(0) - 5*fd.taylor(1) + 4*fd.taylor(2) - fd.taylor(3);
fd.printtaylor(coefs, 7)          # print the 1st 7 nonzero terms of the Taylor series of
                                  # 2f(x[i]) - 5f(x[i+1]) + 4f(x[i+2]) - f(x[i+3])

fd.activatejuliafunction()        # activate Julia function(s) of the newly computed formula in present REPL session
```
