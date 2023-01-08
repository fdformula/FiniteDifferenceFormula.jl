# FiniteDifferenceFormula

This Julia package provides a finite difference formula generator. It generates
finite difference formulas for n-th order derivatives by using Taylor series expansions
of a function at evenly spaced points. It also gives the truncation error of a formula
in the big-O notation. You can use it to generate new formulas in addition to
verification of known ones. See documentation in the source code for the algorithm.

You may play with this package when teaching/learning numerical computing, especially
the finite difference method, and explore the distribution, symmetry, and beauty in
the coefficients of the formulas.

Beware, there is a natural constraint, i.e., the number of points and the order
of derivatives to be in a formula can't be too large, due to

1. the constraint of the largest integer computers can express and process.
1. that we want "exact" formulas, i.e., coefficients in the formulas are integers
   or rational numbers, otherwise, extra truncation errors occur already in the formulas.
1. that, to explore the properties of coefficients in the formulas, we had better have
   "exact" formulas.
1. that practical applications usually don't allow too many points to be used for
   computing derivatives of a function at a single point because of the constraints of
   computing power, not to mention that more involved operations tend to cause more
   rounding errors.

To run the code, you need the Julia programming language (https://julialang.org/), a
wonderful and amazing computing platform.

## The package exports six functions

- computecoefs
- decimalplaces
- formula
- activatejuliafunction
- truncationerror
- taylor

### function computecoefs(n, points, printformulaq = false)

#### Input

```
            n: to find a formula for the n-th order derivative
       points: it can be in the format of a range, start : stop, or
               a vector
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

The algorithm uses the linear combination of f(x[i+j]), j = start : stop (or a specifid
list), to eliminate f(x[i]), f'(x[i]), ..., so that the first nonzero term of the Taylor
series of the linear combination is f^(n)(x[i]):

```Julia
    k[1]*f(x[i+start]) + k[2]*f(x[i+start+1]) + ... + k[stop-start+1]*f(x[i+stop])
        = m*f^(n)(x[i]) + ..., m > 0
```

It is this equation that gives the formula for computing f^(n)(x[i]) and the truncation
error in the big-O notation as well.

### function decimalplaces() or decimalplaces(n)

Without argument, the function returns current decimal places. With arugment n, it sets the
decimal places to be n for generating Julia function for a formula. Without/before calling
the function, 16 decimal places are used by default.

### function formula()

The function generates and lists

1. k[1]*f(x[i+start]) + k[2]*f(x[i+start+1]) + ... + k[stop-start+1]*f(x[i+stop])
       = m*f^(n)(x[i]) + ..., m > 0

1. The formula for finding f^(n)(x[i]), including estimation of accuracy in the big-O
   notation.

1. Julia function for f^(n)(x[i]).

### function activatejuliafunction()

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

### function truncationerror()

The function shows the truncation error of the newly computed finite difference formula.

### function taylor(j, n = 10)

The function prints the first n terms of Taylor series expansion of f(x[i+j]) about x[i].
It is simply for teaching/learning the finite difference method.

## Examples

```Julia
import FiniteDifferenceFormula as fd

fd.decimalplaces(6)               # use 6 decimal places to generate Julia functions of computed formulas

fd.computecoefs(1, 0:2, true)     # find, generate, and print "3"-point forward formula for f'(x[i])

fd.computecoefs(2, -3:0, true)    # find, generate, and print "4"-point backward formula for f''(x[i])

fd.computecoefs(3, -9:9)          # find "19"-point central formula for f'''(x[i])

fd.computecoefs(2, [-3 -2 1 2 7]) # find 5-point formula for f''(x[i])

fd.computecoefs(3,-100:122)       # find "223"-point formula for f'''(x[i])

fd.formula()                      # generate and print the formula computed last time you called computecoefs(...)

fd.truncationerror()              # print the truncation error of the newly computed formula

fd.taylor(-2)                     # print the first 10 terms of the Taylor series of f(x[i-2]) about x[i]

fd.taylor(2, 7)                   # print the first 7 terms of the Taylor series of f(x[i+2]) about x[i]

fd.activatejuliafunction()        # activate Julia function(s) of the newly computed formula in present REPL session
```
