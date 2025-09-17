from decimal import Decimal, getcontext
from typing import List
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def trapinteg(f, a: float, b: float, npoints: int) -> Decimal:
    """
    Integrate function f from a to b using trapeziod summation
    """
    x = np.linspace(a, b, npoints)
    y = f(x)
    h = (b - a) / (npoints - 1) # length of the interval
    w = np.ones(npoints) * h # default weights
    w[0] = w[-1] = h/2

    return np.dot(w, y)  # integral of function f

def simpinteg(f, a: float, b: float, npoints: int) -> Decimal:
    """Integrate function from a to b using Simpson's rule"""
    
    if npoints % 2 == 0:
        raise ValueError("Third argument (npoints) must be odd for Simpson's rule")
    
    x = np.linspace(a, b, npoints)
    y = f(x)
    h = (b - a) / (npoints - 1) # length of the interval
    w = np.ones(npoints) * 2 # set baseline to 2
    w[1:-1:2] = 4 # ignoring endpoints, change every second weight to 4
    w[0] = w[-1] = 1 # correct endpoints to 1

    return (h/3)*np.dot(w, y)

def gqinteg(f, a: float, b: float, npoints: int) -> Decimal:
    """Integrate function from a to b using Gaussian quadrature"""
    x, w = np.polynomial.legendre.leggauss(npoints)
    t = 0.5*(x + 1)*(b - a) + a
    w = 0.5*(b - a)*w
    return np.dot(w, f(t))

def invexp(t: float):
    return np.exp(-t)

def slope(table, column, pt1, pt2):

    err1 = table[column].iloc[pt1]
    err2 = table[column].iloc[pt2]
    N1 = table["N"].iloc[pt1]
    N2 = table["N"].iloc[pt2]

    if np.isnan(err1) or np.isnan(err2):
        return np.nan
    else: 
        return (np.log10(err2) - np.log10(err1))/(np.log10(N2) - np.log10(N1))

    
# Example usage and testing
if __name__ == "__main__":

    real = 1 - np.exp(-1)
    a = 0.0
    b = 1.0
    f = invexp

    alldata = []

    for N in range(3, 1000):
        try:
            tpzd_val = trapinteg(f, a, b, N)
            tpzd_err = abs((tpzd_val - real)/real)
        except Exception:
            tpzd_val, tpzd_err = np.nan, np.nan
        
        if N % 2 == 1:
            simp_val = simpinteg(f, a, b, N)
            simp_err = abs((simp_val - real)/real)
        else:
            simp_val, simp_err = np.nan, np.nan

        if N <= 50:
            gsqd_val = gqinteg(f, a, b, N)
            gsqd_err = abs((gsqd_val - real)/real)

        alldata.append([N, tpzd_val, tpzd_err, simp_val, simp_err, gsqd_val, gsqd_err])

    table = pd.DataFrame(alldata, columns=["N","Trap. Method", "Trap. Rel. Error", 
                                           "Simp. Method", "Simp. Rel. Error", "Gauss. Quad.", 
                                           "Gauss. Quad. Rel. Error"])
    simp_mask = ~np.isnan(table["Simp. Rel. Error"]) # creates an array the same length of "table" 
    # and reads the entries of "Simp. Rel. Error". If entry is a number, writes it as True

    # make a plot of the errors

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    ax.plot(table["N"], table["Trap. Rel. Error"], '-', label="Trapezoid")
    ax.plot(table["N"][simp_mask], table["Simp. Rel. Error"][simp_mask], '-', label="Simpson's")
    ax.plot(table["N"], table["Gauss. Quad. Rel. Error"], '-', label="Gaussian Quadrature")
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.ylim(bottom=5*10e-18)
    ax.set_xlabel("Number of Divisions $N$", fontsize=15)
    ax.set_ylabel(r'|$\varepsilon_{rel}$|', fontsize=15)
    ax.legend(fontsize=9)

    plt.savefig("Errors.png", dpi=300)

    # Print out the table of data
    print(f"\n\n'Real' Value = {real}\n\n")
    print(table)

    # Now let's calculate the slope
    # Based on Landau, for a log log plog, the slope is the order of the number of divisions

    slope_tpzd = slope(table, "Trap. Rel. Error", 50, 52) # for N = 53 to N= 55
    slope_simp = slope(table, "Simp. Rel. Error", 50, 52) # for N = 53 to N = 55
    slope_gsqd = slope(table, "Gauss. Quad. Rel. Error", 1, 2) # for N = 4 to N = 5

    print("\n\nPower Law Dependence: \u03B5 ~ C·N^\u03B1")
    print("Where \u03B1 is the slope:\n")
    print(f"Trapezoid:  \u03B1 ~ ln(\u03B5)/ln(N) ~ {slope_tpzd}")
    print(f"Simpson's:  \u03B1 ~ ln(\u03B5)/ln(N) ~ {slope_simp}")
    print(f"Gaussian:   \u03B1 ~ ln(\u03B5)/ln(N) ~ {slope_gsqd}\n")







