import sys
from sage.all import *
import cProfile
import math



def main():
    # example
    rhs_vals = get_rhs_vals(10000, 100000) # 8 hours
    rhs_precs = get_rhs_vals(10000, 1000, get_precs=True) # 5 minutes
    rhs_derivs = rhs_vals_to_deriv(rhs_vals, 100) # visible numerical artifacts for delta_i < 100

    p1 = points(rhs_vals.items(), dpi=300, size=1)
    x = var('x')
    p2 = plot((x/496.0)**(5.0/6), (x,0,496), dpi=300, color='black')
    show(p1 + p2) # Fig. 1.4
    
    p3 = points(rhs_derivs.items(), dpi=300, size=1, ymax=0.01)
    show(p3) # Fig. 1.5

    from matplotlib.ticker import FormatStrFormatter
    precs_plt = points(rhs_precs.items(), dpi=300, size=1)
    precs_plt.SHOW_OPTIONS['tick_formatter']=(None, FormatStrFormatter('%.1E')) # otherwise scientific notation exponent is dropped; known bug
    show(precs_plt) # Errors are of size ~10^-6

    to_ret = [rhs_vals, rhs_precs, rhs_derivs, p1, p2, p3, precs_plt]
    return to_ret


def integrand(beta, ell):
    # ell is lambda; can't use lambda in python
    C3 = -(27.0/4 * beta**2 + ell/64.0)
    C = RDF(C3).nth_root(3)
    upper = max([1, C])
    lower = max([-1, C])
    to_ret = 1.0/4 * (upper - lower)
    return to_ret


def integral_beta(ell):
    func = lambda beta: integrand(beta, ell)
    intval, interr = numerical_integral(func, -1, 1)
    to_ret = (intval, interr)
    return to_ret


def rho(p,m):
    p_part = ZZ(m).p_primary_part(p)
    n = log(p_part, p)
    if p == 3:
        if n == 0:
            rhoval = 1
        else:
            rhoval = 0    
    elif p == 2:
        if n == 0:
            rhoval = 1.0/2
        elif n in [1,2]:
            rhoval = 1.0/4
        else:
            rhoval = 0
    else:
        if n == 0:
            rhoval = 1 - 1.0/p**2
        elif n == 1:
            rhoval = 1.0/p**2 - 1.0/p**3
        elif n == 2:
            rhoval = 1.0/p**3 - 1.0/p**4
        elif n == 3:
            rhoval = 1.0/p**4 - 2.0/p**5 + 1.0/p**6
        elif n == 4:
            rhoval = 2.0/p**5 - 3.0/p**6 + 1.0/p**7
        elif n == 5:
            rhoval = 2.0/p**6 - 4.0/p**7 + 2.0/p**8
        elif n == 6:
            rhoval = 3.0/p**7 - 5.0/p**8 + 2.0/p**9
        elif n == 7:
            rhoval = 3.0/p**8 - 5.0/p**9 + 2.0/p**10
        elif n == 8:
            rhoval = 3.0/p**9 - 5.0/p**10 + 2.0/p**11
        else:
            rhoval = 2.0/p**(n+1) - 4.0/p**(n+2) + 2.0/p**(n+3)
    return rhoval


def get_rho_prod(m):
    rho_prod = rho(2,m) * rho(3,m)
    if rho_prod != 0:
        for p in ZZ(m).prime_divisors():
            if p >= 5:
                rho_prod *= rho(p,m) / (1 - p**(-2.0))
    return rho_prod


def m_term(m, ell, get_precs=False):
    rho_prod = get_rho_prod(m)
    if rho_prod != 0:
        zeta_ratio = zeta(10.0) / zeta(2.0)
        zeta_euler_facs = (1 - 2.0**(-10)) * (1 - 3.0**(-10)) / ((1 - 2.0**(-2)) * (1 - 3.0**(-2)))
        if abs(m*ell) < 496:
            F_val_plus, F_val_plus_err = integral_beta(m*ell)
            F_val_minus, F_val_minus_err = integral_beta(-m*ell)
        else:
            F_val_plus, F_val_plus_err = 1, 0
            F_val_minus, F_val_minus_err = 0, 0
        if not get_precs:
            termval = zeta_ratio * zeta_euler_facs * rho_prod * (F_val_plus - F_val_minus)
        else:
            termval = zeta_ratio * zeta_euler_facs * rho_prod * (F_val_plus_err + F_val_minus_err)
    else:
        termval = 0
    return termval


def rhs(ell, m_max, get_precs=False):
    sumval = 0
    for m in range(1,m_max):
        sumval += m_term(m, ell, get_precs=get_precs)
    return sumval


def get_rhs_vals(m_max, num_ells=1000, rhs_vals=None, get_precs=False):
    if rhs_vals is None:
        rhs_vals = {}
    delta_ell = 496.0/(num_ells-1)
    for i in range(num_ells):
        ell = delta_ell*i
        if ell not in rhs_vals:
            rhs_vals[ell] = rhs(ell, m_max, get_precs=get_precs)
    return rhs_vals


def rhs_vals_to_deriv(rhs_vals, delta_i):
    rhs_deriv_vals = {}
    rhs_items = sorted(list(rhs_vals.items()))
    rhs_deriv_vals = {}
    for i in range(delta_i, len(rhs_items)-delta_i):
        k0,v0 = rhs_items[i - delta_i]
        k1,v1 = rhs_items[i]
        k2,v2 = rhs_items[i + delta_i]
        deriv_val = 0.5*(v1 - v0)/(k1 - k0) + 0.5*(v2 - v1)/(k2 - k1)
        rhs_deriv_vals[k1] = deriv_val
    return rhs_deriv_vals



        
#This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    if '-profile' in sys.argv:
        cProfile.run('main()')
    else:
        main()
