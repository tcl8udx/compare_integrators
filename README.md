# Compare different quadrature rules for integration

There are two examples provided for calculating the weights and abscissas for gaussian quadrature rules, try:

```
make
./gqconstants
```

or

```
python gqconstants.py
```

You can also use the C++ example as a guide to build your own executable

There is no need to look at rules >~25 for Gaussian quadrature.  And you can also stop at ~ 1000 divisions for the trapezoidal and Simpson's rules.  If you run much longer you'll see the numerical errors bevome visible for the trapezoidal, but hyou'll need to think about how to code efficiently or the running time may be very long.





Student (Tristen Lowrey) comments:

Landau questions:
1  Done
2  Done in terminal
3  Done
4  Slopes:
        Trapezoid: -2.037732627533854
        Simpson's: -4.07540553234387
        Gaussian:  -32.5508520285975

Workflow questions:
1  You can break them pretty easily with a singularity. Trap. and Simp. don't even run if there's a singularity and Gauss. never converges. See
   BadErrors_singularity.png
   You can also break it by attempting to integrate a function with many roots, since the methods rely on Taylor expanding and expressing as polynomials. See BadErrors.png
2  See attached .png's
3  See 1 and attached .png's
