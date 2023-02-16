# FiniteDifferenceFormula

This Julia package provides a general finite difference formula generator and a tool
for teaching/learning the finite difference method. It generates finite difference
formulas for derivatives of various orders by using Taylor series expansions of a
function at evenly spaced points. It also gives the truncation error of a formula
in the big-O notation. We can use it to generate new formulas in addition to
verification of known ones. By changing decimal places, we can also see how rounding
errors may affect a result.

Beware, though formulas are mathematically correct, they may not be numerically useful.
This is true especially when we derive formulas for a derivative of higher order. For
example, run compute(9,-5:5), provided by this package, to generate a 10-point
central formula for the 9-th derivative. The formula is mathematically correct, but it
can hardly be put into use for numerical computing without, if possible, rewriting it
in a special way. Similarly, the more points are used, the more precise a formula
is mathematically. However, due to rounding errors, this may not be true numerically.

To run the code, you need the Julia programming language (https://julialang.org/).

Note: This package has been ported to Python, https://github.com/Winux2k/FiniteDifferenceFormula.py.

## How to install the package

In Julia REPL, execute the following two commands in order.

1. import Pkg
1. Pkg.add("FiniteDifferenceFormula")

## The package exports the following functions

```activatejuliafunction```, ```compute```, ```decimalplaces```, ```find```, ```findbackward```,
```findforward```, ```formula```, ```formulas```, ```loadcomputingresults```, ```taylor```,
```taylorcoefs```, ```tcoefs```, ```truncationerror```, ```verifyformula```

### functions, ```compute```, ```find```, ```findforward```, and ```findbackward```

All take the same arguments (n, points, printformulaq = false).

#### Input

```
            n: the n-th order derivative to be found
       points: in the format of a range, start : stop, or a vector
printformulaq: print the computed formula or not
```

|   points     |   The points/nodes to be used                  |
|   ---------- | ---------------------------------------------- |
|   -3:2       |   x[i-3], x[i-2], x[i-1], x[i], x[i+1], x[i+2] |
|   [1 0 1 -1] |   x[i-1], x[i], x[i+1]                         |

A vector can be like [1, 0, 2] or [1 0 2]. It will be rearranged so that elements are ordered
from lowest to highest with duplicate ones removed.

#### Output

Each function returns a tuple, (n, points, k[:], m), where n, points, k[:] and m are described below.
With the information, you may generate functions for any programming language of your choice.

While 'compute' may fail to find a formula using the points, others try to find one, if possible,
by using fewer points in different ways. (See the docstring of each function by, say,
```?fd.find``` after ```import FiniteDifferenceFormula as fd```, in Julia REPL.)

The algorithm uses the linear combination of f(x[i+j]) = f(x[i] + jh), where h is the increment
in x and j ∈ points, to eliminate f(x[i]), f'(x[i]), f''(x[i]), ..., so that the first nonzero
term of the Taylor series of the linear combination is f^(n)(x[i]):

```Julia
k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ... + k[len]*f(x[i+points[len]]) = m*f^(n)(x[i]) + ..., m > 0
```

where len = length(points). It is this equation that gives the formula for computing f^(n)(x[i])
and the truncation error in the big-O notation as well.

### function ```loadcomputingresults(results)```

The function loads results, a tuple of the form (n, points, k, m), returned by ```compute```.
For example, it may take hours to compute/find formulas invloving hundreds of points. In this
case, we can save the results in a text file and come back later to work on the results
with ```activatejuliafunction```, ```formula```, ```truncationerror```, and so on.

### function ```formula()```

The function generates and lists

1. k[1]\*f(x[i+points[1]]) + k[2]\*f(x[i+points[2]]) + ... + k[len]\*f(x[i+points[len]])
       = m\*f^(n)(x[i]) + ..., m > 0

1. The formula for f^(n)(x[i]), including estimation of accuracy in the big-O notation.

1. Julia function(s) for f^(n)(x[i]).

### function ```truncationerror()```

The function returns a tuple, (n, "O(h^n)"), the truncation error of the newly computed finite
difference formula in the big-O notation.

### function ```decimalplaces()``` or ```decimalplaces(n)```

Without an argument, the function returns current decimal places. With argument n, it sets the
decimal places to be n for generating Julia function(s) for formulas if n is a nonnegative
integer. It returns the (new) default decimal places. Without/before calling the function, 16
decimal places are used by default.

This function can only affect Julia functions with the suffix "d" such as fd1stderiv2ptcentrald.
See function activatejuliafunction().

### function ```activatejuliafunction()```

Call this function to activate the Julia function(s) for the newly computed finite
difference formula. For example, after compute(1, -1:1) and decimalplaces(4), it activates the
following Julia functions.

```Julia
fd1stderiv2ptcentrale(f, x, i, h)  = ( -f(x[i-1]) + f(x[i+1]) ) / (2 * h)
fd1stderiv2ptcentrale1(f, x, i, h) = ( -1/2 * f(x[i-1]) + 1/2 * f(x[i+1]) ) / h
fd1stderiv2ptcentrald(f, x, i, h)  = ( -0.5000 * f(x[i-1]) + 0.5000 * f(x[i+1]) ) / h
```
The suffixes 'e' and 'd' stand for 'exact' and 'decimal', respectively. No suffix? It is "exact".
After activating the function(s), we can evaluate right away in the present Julia REPL session. For example,

```Julia
FiniteDifferenceFormula.fd1stderiv2ptcentrale(sin, 0:0.01:pi, 3, 0.01)
```
Below is the output of activatejuliafunction(). It gives us the first chance to examine the usability
of the computed or tested formula.

```Julia
import FiniteDifferenceFormula as fd
f, x, i, h = sin, 0:0.01:10, 501, 0.01
fd.fd1stderiv2ptcentrale(f, x, i, h)   # result: 0.2836574577837647, relative error = 0.00166666%
fd.fd1stderiv2ptcentrale1(f, x, i, h)  # result: 0.2836574577837647, relative error = 0.00166666%
fd.fd1stderiv2ptcentrald(f, x, i, h)   # result: 0.2836574577837647, relative error = 0.00166666%
                                       # cp:     0.2836621854632262
```

### function ```activatejuliafunction(n, points, k, m)``` or ```verifyformula(n, points, k, m)```

They are the same. Each allows users to load a formula from some source to test and see if it is correct.
If it is valid, its truncation error in the big-O notation can be determined. Furthermore, if the input
data is not for a valid formula, it tries also to find one, if possible, using n and points.

Here, n is the order of a derivative, points are a list, k is a list of the corresponding
coefficients of a formula, and m is the coefficient of the term f^(n)(x[i]) in the linear
combination of f(x[i+j]), where j ∈ points. In general, m is the coefficient of h^n in the
denominator of a formula. For example,

```Julia
import FiniteDifferenceFormula as fd
fd.activatejuliafunction(2, [-1 0 2 3 6], [12 21 2 -3 -9], -12)
fd.truncationerror()
fd.activatejuliafunction(4, 0:4, [2//5 -8//5 12//5 -8//3 2//5], 5)
fd.activatejuliafunction(4, [0, 1, 2, 3, 4], [2/5 -8/5 12/5 -8/3 2/5], 5)
fd.activatejuliafunction(2, [-1 2 0 2 3 6], [1.257 21.16 2.01 -3.123 -9.5], -12)
```

### function ```taylorcoefs(j, n = 10)``` or ```tcoefs(j, n = 10)```

The function returns the coefficients of the first n terms of the Taylor series of f(x[i+j])
about x[i].

### function ```taylor(j, n = 10)```

The function prints the first n terms of the Taylor series of f(x[i+j]) about x[i].

### function ```taylor(coefs, n = 10)``` or ```taylor(points, k, n = 10)```

The function prints the first n nonzero terms of a Taylor series of which the coefficients are
provided in ```coefs``` or given through ```points``` and ```k[:]``` as in the linear combination
```k[1]*f(x[i+points[1]]) + k[2]*f(x[i+points[2]]) + ...``` It provides also another way
to verify if a formula is correct.

### function ```formulas(highest_order = 3, max_num_of_points = 5)```

By default, the function prints all forward, backward, and central finite difference formulas for
the 1st, 2nd, and 3rd derivatives, using at most 5 points.

## Examples

```Julia
import FiniteDifferenceFormula as fd
fd.compute(1, 0:2, true)             # find, generate, and print "3"-point forward formula for f'(x[i])
fd.compute(2, -3:0, true)            # find, generate, and print "4"-point backward formula for f''(x[i])
fd.compute(3, -9:9)                  # find "19"-point central formula for f'''(x[i])
fd.decimalplaces(6)                  # use 6 decimal places to generate Julia functions of computed formulas
fd.compute(2, [-3 -2 1 2 7])         # find formula for f''(x[i]) using points x[i+j], j = -3, -2, 1, 2, and 7
fd.compute(1,-230:230)               # find "461"-point central formula for f'(x[i]). it may take hours!
fd.formula()                         # generate and print the formula computed last time you called compute(...)
fd.truncationerror()                 # print and return the truncation error of the newly computed formula
fd.taylor(-2, 5)                     # print the first 5 terms of the Taylor series of f(x[i-2]) about x[i]
coefs = 2*fd.tcoefs(0) - 5*fd.tcoefs(1) + 4*fd.tcoefs(2) - fd.tcoefs(5);
fd.taylor(coefs, 7)                  # print the 1st 7 nonzero terms of the Taylor series of 2f(x[i]) - 5f(x[i+1]) + 4f(x[i+2]) - f(x[i+5])
fd.taylor([0, 1, 2, 5], [2, -5, 4, -1], 7)
fd.activatejuliafunction()           # activate Julia function(s) of the newly computed formula in present REPL session
fd.verifyformula(1, 2:3, [-4, 5], 6) # verify if f'(x[i]) = (-4f(x[i+2] + 5f(x[i+3)) / (6h) is a valid formula
fd.formulas(5, 10)                   # print all forward, backword, and central formulas for the 1st, 2nd, ..., 5th derivatives, using at most 10 points
```
