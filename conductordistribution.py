import sys
from sage.all import *
import cProfile
import math
import random
import os

import reductiontables


YOUNGFAMILY = True#False

WRITE_FILENAME_PREFIX = 'rhs_vals_youngfamily'#'rhs_vals'


DATA_TYPE = 'diff'
REQUIRE_GLOBALLY_MINIMAL = True
LOCAL_TYPE_DATA_2_FILENAME = 'short_weierstrass_models_2_type.sobj'
LOCAL_TYPE_DATA_3_FILENAME = 'short_weierstrass_models_3_type_2.sobj'



def example():
    rhs_vals = get_rhs_vals(num_ells=10**5, m_max=10**4) # 8 hours
    rhs_precs = get_rhs_vals(num_ells=1000, m_max=10**4, get_precs=True) # 5 minutes
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


def main():
    # Make sure YOUNGFAMILY is set correctly
    num_instances = int(sys.argv[1])
    instance_index = int(sys.argv[2])
    num_ells = int(sys.argv[3])
    if len(sys.argv) > 4:
        m_max = int(sys.argv[4])
    else:
        m_max = None

    write_filename_prefix = WRITE_FILENAME_PREFIX
    write_filename_suffix = '_%s_%s_%s' % (str(num_instances), str(instance_index), str(num_ells))
    if m_max:
        write_filename_suffix += '_%s' % str(m_max)
    write_filename = write_filename_prefix + write_filename_suffix

    if m_max:
        parallelize_single(write_filename, num_instances, instance_index, num_ells, m_max=m_max)
    else:
        parallelize_single(write_filename, num_instances, instance_index, num_ells)
    return


def parallelize(num_instances, num_ells, m_max=None):
    test_get_index_range(num_instances, num_ells)
    for instance_index in range(num_instances):
        nohupsuffix = '_' + WRITE_FILENAME_PREFIX + '_%s_%s_%s' % (str(num_instances), str(instance_index), str(num_ells))
        if m_max is not None:
            nohupsuffix += '_%s' % str(m_max)
        command = ''
        command += 'nohup sage conductordistribution.py'
        command += ' %s %s %s' % (str(num_instances), str(instance_index), str(num_ells))
        if m_max is not None:
            command += ' %s' % str(m_max)
        command += ' > nohup' + nohupsuffix + '.out &'
        os.system(command) # apparently it's better to use subprocess.run
    return


def parallelize_single(write_filename, num_instances, instance_index, num_ells, m_max=None):
    index_range = get_index_range(num_instances, instance_index, num_ells)
    rhs_vals = get_rhs_vals(num_ells=num_ells, index_range=index_range, m_max=m_max)
    write_rhs_vals(write_filename, rhs_vals)
    return


def write_rhs_vals(write_filename, rhs_vals):
    save(rhs_vals, write_filename)
    return


def read_rhs_vals(filename_prefix, filename_suffix, data_path='./'):
    rhs_vals = {}
    for fname in os.listdir(data_path):
        if (fname[:len(filename_prefix)] == filename_prefix) and (fname[-len(filename_suffix):] == filename_suffix):
            dat = load(fname)
            rhs_vals.update(dat)
    return rhs_vals


def get_index_range(num_instances, instance_index, num_ells):
    full_range = list(range(num_ells))
    delta_index = int(ceil(ZZ(num_ells)/num_instances))
    assert(delta_index > 0)
    start = delta_index * instance_index
    end = delta_index * (instance_index+1)
    index_range = full_range[start:end]
    return index_range


def test_get_index_range(num_instances, num_ells):
    all_ranges = {i:set(get_index_range(num_instances,i,num_ells)) for i in range(num_instances)}
    full_range = set(range(num_ells))
    range_union = set()
    for ir in all_ranges.values():
        range_union = range_union.union(ir)
    if full_range != range_union and len(range_union) < 302:
        print(range_union)
    assert(full_range == range_union)
    for i in all_ranges:
        for j in all_ranges:
            if i != j:
                assert(len(all_ranges[i].intersection(all_ranges[j])) == 0)
    return


###

def get_rhs_vals(num_ells=1000, m_max=None, index_range=None, rhs_vals=None, get_precs=False, youngfamily=YOUNGFAMILY, cached=None):
    #if cached is None:
    #    cached = {}
    if rhs_vals is None:
        rhs_vals = {}
    delta_ell = ZZ(496)/(num_ells-1)
    if index_range is None:
        index_range = list(range(num_ells))
    for i in index_range:
        ell = delta_ell*i
        if ell not in rhs_vals:
            rhs_vals[ell] = rhs(ell, m_max=m_max, get_precs=get_precs, youngfamily=youngfamily, cached=cached)
    return rhs_vals


def rhs(ell, m_max=None, get_precs=False, youngfamily=YOUNGFAMILY, cached=None):
    m_max_was_None = m_max is None
    if m_max is None:
        if ell == 0:
            return 0
        # zeta(10)/zeta(2) sum_{m=1}^\infty (F(m*ell) - F(-m*ell)) prod_p rho(p,m)/(1 - p^(-2))
        # = 1 + zeta(10)/zeta(2) sum_{1 <= m <= 496/ell} (F(m*ell) - F(-m*ell) - 1) prod_p rho(p,m)/(1 - p^(-2))
        m_max = ceil(496.0/ell) + 1 # the above also holds for m <= 496/ell' for any ell' < ell
    sumval = 0
    for m in range(1,m_max+1):
        termval = m_term(m, ell, get_precs=get_precs, youngfamily=youngfamily, cached=cached, shift_by_1=m_max_was_None)
        sumval += termval
    if m_max_was_None:
        sumval += 1
    return sumval


def m_term(m, ell, get_precs=False, youngfamily=YOUNGFAMILY, cached=None, shift_by_1=False):
    rho_prod = get_rho_prod(m, youngfamily=youngfamily)
    if rho_prod != 0:
        zeta_ratio = zeta(10.0) / zeta(2.0)
        if youngfamily:
            zeta_euler_facs = (1 - 2.0**(-10)) * (1 - 3.0**(-10)) / ((1 - 2.0**(-2)) * (1 - 3.0**(-2)))
        else:
            zeta_euler_facs = 1
        if abs(m*ell) < 496:
            F_val_plus, F_val_plus_err = integral_beta(m*ell, cached=cached)
            F_val_minus, F_val_minus_err = integral_beta(-m*ell, cached=cached)
        else:
            F_val_plus, F_val_plus_err = 1, 0
            F_val_minus, F_val_minus_err = 0, 0
        if shift_by_1:
            F_val_plus -= 1
        if not get_precs:
            termval = zeta_ratio * zeta_euler_facs * rho_prod * (F_val_plus - F_val_minus)
        else:
            termval = zeta_ratio * zeta_euler_facs * rho_prod * (F_val_plus_err + F_val_minus_err)
    else:
        termval = 0
    return termval


def integrand(beta, ell):
    # ell is lambda; can't use lambda in python
    C3 = -(27.0/4 * beta**2 + ell/64.0)
    C = RDF(C3).nth_root(3)
    upper = max([1, C])
    lower = max([-1, C])
    to_ret = 1.0/4 * (upper - lower)
    return to_ret


def integral_beta(ell, cached=None):
    if (cached is None) or (ell not in cached):
        ellf = float(ell)
        func = lambda beta: integrand(beta, ellf)
        intval, interr = numerical_integral(func, -1, 1)
        to_ret = (intval, interr)
        if cached is not None:
            cached[ell] = to_ret
    else:
        to_ret = cached[ell]    
    return to_ret


def rho(p,m, youngfamily=YOUNGFAMILY):
    n = ZZ(m).valuation(p)
    if p == 3:
        if youngfamily:
            if n == 0:
                rhoval = 1
            else:
                rhoval = 0
        else:
            if n < 13:
                vals3 = [8.0/3**2,
                         2.0/3**3,
                         2.0/3**4,
                         0,
                         2.0/3**5,
                         2.0/3**7,
                         10.0/3**8,
                         10.0/3**9,
                         10.0/3**10,
                         4.0/3**11,
                         4.0/3**12,
                         4.0/3**13,
                         148.0/3**14,
                         ]
                rhoval = vals3[n]
            else:
                rhoval = 40.0/3**(n+2)
    elif p == 2:
        if youngfamily:
            if n == 0:
                rhoval = 1.0/2
            elif n in [1,2]:
                rhoval = 1.0/4
            else:
                rhoval = 0
        else:
            if n < 13:
                vals2 = [1.0/2,
                         1.0/2**2,
                         1.0/2**3,
                         0,
                         1.0/2**4,
                         1.0/2**6,
                         3.0/2**7,
                         3.0/2**7,
                         3.0/2**8,
                         3.0/2**9,
                         1.0/2**10,
                         1.0/2**11,
                         1.0/2**12,
                         21.0/2**13,
                         ]
                rhoval = vals2[n]
            else:
                rhoval = 5.0/2**(n+1)
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


def get_rho_prod(m, youngfamily=YOUNGFAMILY):
    rho_prod = rho(2,m,youngfamily=youngfamily) * rho(3,m,youngfamily=youngfamily)
    if rho_prod != 0:
        for p in ZZ(m).prime_divisors():
            if p >= 5:
                rho_prod *= rho(p,m) / (1 - p**(-2.0))
        if not youngfamily:
            rho_prod /= (1 - 2.0**(-2))
            rho_prod /= (1 - 3.0**(-2))
    return rho_prod


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



###


def is_good(a, b, r=None, t=None, require_globally_minimal=REQUIRE_GLOBALLY_MINIMAL):
    if r is not None:
        if a%6 != r:
            return False
    if t is not None:
        if b%6 != t:
            return False
    if 4*a**3 + 27*b**2 == 0:
        return False
    if require_globally_minimal:
        if a != 0:
            afac = list(factor(a))
            for (p,m) in afac:
                if m >= 4:
                    if b%p**6 == 0:
                        return False
        else:
            bfac = list(factor(b))
            for (p,m) in bfac:
                if m >= 6:
                    return False
    return True


def make_default_toupdate():
    return {(r,t):{} for r in range(6) for t in range(6)}


def get_reduction_types(kmin, kmax, to_update=None):
    if to_update is None:
        to_update = make_default_toupdate()
    for k in range(kmin, kmax):
        amin = 6**(2*k)
        amax = 6**(2*(k+1))
        if amin == 1:
            amin = 0
        bmin = 6**(3*k)
        bmax = 6**(3*(k+1))
        if bmin == 1:
            bmin = 0
        for a in range(amin,amax):
            if a%2 == 0: #otherwise already done
                for b in range(bmin,bmax):
                    if b%3 == 0: #otherwise already done
                        if is_good(a,b):
                            r = a%6
                            t = b%6
                            ld = extract_local_data(a,b)
                            to_update[(r,t)][(a,b)] = ld
    return to_update


def extract_local_data(a,b):
    E = EllipticCurve([a,b])
    ld2 = extract_local_data_p(E,2)
    ld3 = extract_local_data_p(E,3)
    return (ld2, ld3)


def extract_local_data_p(E, p, data_type=DATA_TYPE):
    ld = E.local_data(p)
    if data_type == 'diff':
        Nv = ld.conductor_valuation()
        Dv = ld.discriminant_valuation()
        diff = Dv - Nv
        to_ret = diff
    else:
        to_ret = str(ld.kodaira_symbol())
    return to_ret

            
def get_reduction_types_p(p, kmin, kmax, data_type=DATA_TYPE, require_globally_minimal=REQUIRE_GLOBALLY_MINIMAL, to_update=None):
    if to_update is None:
        to_update = {}
    for k in range(kmin, kmax):
        amin = p**(2*k)
        amax = p**(2*(k+1))
        if amin == 1:
            amin = 0
        bmin = p**(3*k)
        bmax = p**(3*(k+1))
        if bmin == 1:
            bmin = 0
        for a in range(amin,amax):
            if (p != 2) or (a%2 == 0): #otherwise already done
                for b in range(bmin,bmax):
                    if (p != 3) or (b%3 == 0): #otherwise already done
                        if is_good(a, b, require_globally_minimal=require_globally_minimal):
                            E = EllipticCurve([a,b])
                            ld = extract_local_data_p(E, p, data_type=data_type)
                            try:
                                to_update[ld].append((a,b))
                            except KeyError:
                                to_update[ld] = []
                                to_update[ld].append((a,b))
    return to_update
    

def load_local_data(p, filename=None):
    if filename is not None:
        data = load(filename)
    else:
        if p == 2:
            data = load(LOCAL_TYPE_DATA_2_FILENAME)
        elif p == 3:
            data = load(LOCAL_TYPE_DATA_3_FILENAME)
        else:
            raise NotImplementedError
    return data


def get_reversed_val(v, min_len=None, give_str=False, p=2):
    v = ZZ(v)
    #vb = v.binary()
    vb = v.str(p)
    if (min_len is not None) and (len(vb) < min_len):
        vb = '0'*(min_len - len(vb)) + vb
    vbr = vb[::-1]
    vr = ZZ(vbr, base=p)
    if give_str:
        to_ret = vbr
    else:
        to_ret = vr
    return to_ret


def get_vred(v, modulus):
    vred = set()
    for ab in v:
        vred.add((ab[0]%modulus, ab[1]%modulus))
    return vred


def get_region_plot(ab_vals, min_len, p=2, pt_scale=1, region_only=False, polygon_plot=False, **kwargs):
    pts = []
    for ab in ab_vals:
        ar = get_reversed_val(ab[0], min_len, p=p)
        br = get_reversed_val(ab[1], min_len, p=p)
        pts.append((ar,br))
    if not region_only:
        if polygon_plot:
            poly_plot = squareify(pts, pt_scale=pt_scale, **kwargs)
        else:
            poly_plot = points(pts, **kwargs)
        to_ret = [pts, poly_plot]
    else:
        to_ret = pts
    return to_ret


def get_regions(local_data, min_len=None, num_titles=1, **pts_kwargs):
    p = 2
    num_abs = sum([len(tmp) for tmp in local_data.values()])
    if min_len is None:
        min_len = ceil(log(num_abs, p))
    modulus = p**min_len
    region_list = []
    text_list = []
    for i,(k,v) in enumerate(local_data.items()):
        #color = hue((1 - sqrt(5.0))/2 * i)
        #color = hue(i/float(len(local_data)))
        #shift = ((1 + sqrt(5.0))/2 * i)%1
        shift = i%2
        color = hue(i/float(len(local_data)), v = shift*0.2 + 0.8, s = (1-shift)*0.2 + 0.8)
        if k[0][:2] in ['I' + str(tmp) for tmp in range(10)]:
            if k[0][-1] == '*':
                #k0str = '\\mathrm{' + k[0][0] + '}^*' + '_{' + k[0][1:-1] + '}' #+ '^*'
                k0str = '\\mathrm{' + k[0][0] + '}' + '_{' + k[0][1:-1] + '}' + '^*'
            else:
                k0str = '\\mathrm{' + k[0][0] + '}' + '_{' + k[0][1:] + '}'
        else:
            if k[0][-1] == '*':
                k0str = '\\mathrm{' + k[0][:-1] + '}^*'
            else:
                k0str = '\\mathrm{' + k[0] + '}'
        region_title = '$' + k0str + '(%i)' % k[1] + '$'
        vred = set()
        for ab in v:
            vred.add((ab[0]%modulus, ab[1]%modulus))
        #region, region_plot = get_region_plot(vred, min_len, dpi=200, wireframe='black', fill='white', xmin=0, xmax=modulus, ymin=0, ymax=modulus)
        region, region_plot = get_region_plot(vred, min_len, p=p, dpi=200, xmin=0, xmax=modulus, ymin=0, ymax=modulus, color=color, **pts_kwargs)
        #region_text = text(str(k), region.centroid(), color='black', xmin=0, xmax=modulus, ymin=0, ymax=modulus)
        region_list.append(region_plot)
        region.sort()
        for j in range(num_titles):
            #region_text = text(str(k), region[random.randint(0,len(region)-1)], color='black', xmin=0, xmax=modulus, ymin=0, ymax=modulus)
            #index = int(j*len(region)/num_titles)
            #index = random.randint(0,len(region)-1)
            #index = int(j*len(region)/num_titles) + random.randint(-int(len(region)/(3.0*num_titles)), int(len(region)/(3.0*num_titles)))
            index = int( ((j+1) * (1+sqrt(5.0))/2 * len(region)) ) % len(region)
            region_text = text(region_title, region[index], color='black', xmin=0, xmax=modulus, ymin=0, ymax=modulus, fontsize='small')
            text_list.append(region_text)
    return region_list + text_list



###

def test_type_locations(p=3, run_asserts=True):
    if p == 2:
        type_data = reductiontables.get_type_data_2()
        min_len = 10
    elif p == 3:
        type_data = reductiontables.get_type_data_3()
        min_len = 8
    else:
        raise NotImplementedError

    modulus = p**min_len

    for td in type_data:
        unchecked_locations = set(td.locations)
        while unchecked_locations:
            loc = unchecked_locations.pop()
            a,b = loc
            
            if a in ZZ:
                al = [a]
            else:
                al = [int(a), int(a)+1]
            if b in ZZ:
                bl = [b]
            else:
                bl = [int(b), int(b)+1]
            
            if (len(al) > 1) or (len(bl) > 1):
                for atmp in al:
                    for btmp in bl:
                        unchecked_locations.add((atmp, btmp))
            else:
                ar = int(loc[0])
                astr = ZZ(ar).str(p)
                astr = '0'*(min_len - len(astr)) + astr
                astr_r = astr[::-1]
                a = ZZ(astr_r, base=p)

                br = int(loc[1])
                bstr = ZZ(br).str(p)
                bstr = '0'*(min_len - len(bstr)) + bstr
                bstr_r = bstr[::-1]
                b = ZZ(bstr_r, base=p)

                E = EllipticCurve([a,b])

                N = E.conductor()
                Np = ZZ(N).p_primary_part(p)

                D = E.discriminant()
                Dp = ZZ(D).p_primary_part(p)
                #Dm = E.minimal_discriminant_ideal().gens()[0]
                #Dmp = ZZ(Dm).p_primary_part(p)

                ldat = E.local_data(p)

                Nv = ldat.conductor_valuation()
                N_ldat = p**Nv

                Dv = ldat.discriminant_valuation()
                D_ldat = p**Dv
                if not E.is_p_minimal(p):
                    D_ldat *= p**12

                kodaira_symbol = str(ldat.kodaira_symbol())

                N_manual_good = Np == td.conductor
                N_ldat_good = N_ldat == td.conductor
                D_manual_good = Dp == td.discriminant
                D_ldat_good = D_ldat == td.discriminant
                ratio_good = Dp/Np == td.ratio
                kodaira_symbol_good = kodaira_symbol == td.symbol
                td_good = all([N_manual_good, N_ldat_good, D_manual_good, D_ldat_good, ratio_good, kodaira_symbol_good])

                error_message = ''
                try:
                    error_message += 'loc = ' + str(loc) + ' = (' + ZZ(loc[0]).str(p) + ', ' + ZZ(loc[1]).str(p) + ')'
                except TypeError:
                    error_message += 'loc = ' + str(loc)
                    td_good = False
                error_message += '\n'
                error_message += 'a, b = %i, %i' % (a,b)
                error_message += '\n'
                error_message += 'N: %i,  Np: %i,  td.conductor: %i' % (N, Np, td.conductor)
                error_message += '\n'
                #error_message += 'D: %i,  Dm: %i,  Dp: %i,  td.discriminant: %i' % (D, Dm, Dp, td.discriminant)
                error_message += 'D: %i,  Dp: %i,  td.discriminant: %i' % (D, Dp, td.discriminant)
                error_message += '\n'
                error_message += 'td.ratio: ' + str(td.ratio)
                error_message += '\n'
                error_message += 'kodaira_symbol:  %s,  td.symbol:  %s' % (kodaira_symbol, td.symbol)
                error_message += '\n'
                error_message += str(ldat)

                if not td_good:
                    print(error_message)

                if run_asserts:
                    assert Np == td.conductor, 'loc: ' + str(loc) + ', a,b = %i,%i, N: %i, Np: %i, td.conductor: %i' % (a, b, N, Np, td.conductor)
                    assert Dp == td.discriminant, 'loc: ' + str(loc) + ', a,b = %i,%i, D: %i, Dp: %i, td.discriminant: %i' % (a, b, D, Dp, td.discriminant)
                    assert Dp/Np == td.ratio, 'loc: ' + str(loc) + ', a,b = %i,%i, Np: %i, td.conductor: %i, Dp: %i, td.discriminant: %i' % (a, b, Np, td.conductor, Dp, td.discriminant)
                    assert N_ldat == td.conductor
                    assert D_ldat == td.discriminant
                    assert D_ldat/N_ldat == td.ratio
                    assert kodaira_symbol == td.symbol

                    assert N_manual_good
                    assert N_ldat_good
                    assert D_manual_good
                    assert D_ldat_good
                    assert ratio_good, 'Dp/Np: %i,  td.ratio: %i' % (Dp/Np, td.ratio)
                    assert kodaira_symbol_good
                    assert td_good

    return


###


def get_region_pts_from_file(type_name, p=3):
    filename = type_name_to_region_pts_filename(type_name, p=p)
    region_pts = read_region(filename)
    return region_pts


def get_squares_from_file(type_name, p=3):
    filename = type_name_to_squares_filename(type_name, p=p)
    squares = read_squares(filename)
    return squares


def read_region(filename):
    #region = load(filename)
    region = []
    with open(filename,'r') as f:
        for line in f:
            line = line.strip()
            line = line.replace('(','')
            line = line.replace(')','')
            pt0, pt1 = line.split(',')
            pt0 = int(pt0)
            pt1 = int(pt1)
            region.append((pt0, pt1))
    return region


def read_squares(filename):
    squares = {}
    with open(filename,'r') as f:
        for line in f:
            line = line.strip()
            x, y, size = line.split(',')
            x = float(x)
            y = float(y)
            size = float(size)
            squares[(x,y)] = size
    return squares


def write_region(region, filename):
    #save(region, filename)
    region = sorted(region)
    with open(filename,'w') as f:
        first_pt = True
        for pt in region:
            if first_pt:
                first_pt = False
            else:
                tmp = f.write('\n')
            ptstr = str(pt)
            ptstr.replace(' ','')
            tmp = f.write(ptstr)
    return


def write_squares(square_list, filename):
    with open(filename,'w') as f:
        first_line = True
        for xy,size in sorted(list(square_list.items())):
            if first_line:
                first_line = False
            else:
                tmp = f.write('\n')
            tmp = f.write(str(xy[0]) + ',' + str(xy[1]) + ',' + str(size))
    return


def type_name_to_region_pts_filename(type_name, p=3):
    type_name = type_name.replace(' ','_')
    type_name = type_name.replace('^','pow')
    type_name = type_name.replace('*','star')
    filename = 'region_points_' + str(type_name)
    if p != 3:
        filename += '_' + str(p)
    return filename


def type_name_to_squares_filename(type_name, p=3):
    type_name = type_name.replace(' ','_')
    type_name = type_name.replace('^','pow')
    type_name = type_name.replace('*','star')
    filename = 'region_points_' + str(type_name) + '_squares'
    if p != 3:
        filename += '_' + str(p)
    return filename


def make_region_pts_files(local_data=None, squares_too=True, overwrite=False, p=3):
    if p == 3:
        if local_data is None:
            local_data = load_local_data(3)
        min_len = 8
        modulus = 3**min_len
    elif p == 2:
        if local_data is None:
            local_data = load_local_data(2)
        min_len = 10
        modulus = 2**min_len
    else:
        raise NotImplementedError

    if p == 3:
        type_data = reductiontables.get_type_data_3(with_divisible=True)
    elif p == 2:
        type_data = reductiontables.get_type_data_2(with_divisible=True)
    td_div = type_data[-1]
    assert(td_div.name == 'divisible')
    type_data = type_data[:-1]
    type_data = [td for td in type_data if td.modulus <= modulus]
    type_data = [td for td in type_data if td.dict_key() in local_data]        
    
    # p^4|a and p^6 region setup
    divisible_region = set()
    if p == 3:
        for a in range(int(modulus * ZZ(3)**(-4))):
            for b in range(int(modulus * ZZ(3)**(-6))):
                divisible_region.add((a,b))
    elif p == 2:
        for a in range(int(modulus/8)):
            for b in range(int(modulus/8)):
                divisible_region.add((a,b))

    # Create region pts for each type
    for i, td in enumerate(type_data):
        print(td.name)
        k = td.dict_key()
        v = local_data[k]
        vred = get_vred(v, modulus)
        region = get_region_plot(vred, min_len, p=p, region_only=True)
        divisible_region.difference_update(set(region))

        if p == 3:
            if td.name == 'I0':
                region = [pt for pt in region if pt[0] >= 2*modulus/3]
            elif td.name == '3^12 I0':
                region = [pt for pt in region if pt[0] < 2*modulus/3]
        
        filename = type_name_to_region_pts_filename(td.name, p=p)
        existing_filenames = os.listdir('./')
        if overwrite or (filename not in existing_filenames):
            write_region(region, filename)

        if squares_too:
            squares = points_to_squares(region)
            squares_filename = type_name_to_squares_filename(td.name, p=p)
            existing_filenames = os.listdir('./')
            if overwrite or (squares_filename not in existing_filenames):
                write_squares(squares, squares_filename)

    divisible_filename = type_name_to_region_pts_filename('divisible', p=p)
    write_region(divisible_region, divisible_filename)
    if squares_too:
        squares = points_to_squares(divisible_region)
        squares_filename = type_name_to_squares_filename('divisible', p=p)
        existing_filenames = os.listdir('./')
        if overwrite or (squares_filename not in existing_filenames):
            write_squares(squares, squares_filename)
    return


def get_region_figure_2(region_pts_from_file=True, local_data=None, dpi=None, text_dpi=None, figsize=20, min_len=10, pt_scale=1, axes_width_scale_factor=1.3, random_text_loc=False, polygon_plot=True):
    if region_pts_from_file:
        assert(min_len == 10)
    if (not region_pts_from_file) and (local_data is None):
        local_data = load_local_data(2)
    modulus = 2**min_len

    region_plot_list = []
    text_list = []
    type_data = reductiontables.get_type_data_2(with_divisible=True)
    td_div = type_data[-1]
    assert(td_div.name == 'divisible')
    type_data = type_data[:-1]
    type_data = [td for td in type_data if td.modulus <= modulus]
    if region_pts_from_file:
        type_data = [td for td in type_data if type_name_to_region_pts_filename(td.name, p=2) in os.listdir('./')]
    else:
        type_data = [td for td in type_data if td.dict_key() in local_data]
    
    # ['$\\mathrm{II}(4)$', '$\\mathrm{II}(6)$', '$\\mathrm{II}(7)$', '$\\mathrm{I}_0$', '$\\mathrm{I}_1$', '$\\mathrm{I}_2$', '$\\mathrm{I}_3$', '$\\mathrm{III}(3)$', '$\\mathrm{III}(5)$', '$\\mathrm{III}(7)$', '$\\mathrm{III}(8)$', '$\\mathrm{IV}$', '$\\mathrm{I}^{\\!*}_0\\!(4)$', '$\\mathrm{I}^{\\!*}_0\\!(5)$', '$\\mathrm{I}^{\\!*}_0\\!(6)$', '$\\mathrm{I}^{\\!*}_1$', '$\\mathrm{I}^{\\!*}_2\\!(4)$', '$\\mathrm{I}^{\\!*}_2\\!(6)$', '$\\mathrm{I}^{\\!*}_2\\!(7)$', '$\\mathrm{I}^{\\!*}_3\\!(4)$', '$\\mathrm{I}^{\\!*}_3\\!(5)$', '$\\mathrm{I}^{\\!*}_4\\!(4)$', '$\\mathrm{I}^{\\!*}_4\\!(6)$', '$\\mathrm{I}^{\\!*}_5\\!(4)$', '$\\mathrm{I}^{\\!*}_5\\!(6)$', '$\\mathrm{I}^{\\!*}_6\\!(4)$', '$\\mathrm{I}^{\\!*}_6\\!(6)$', '$\\mathrm{I}^{\\!*}_7\\!(4)$', '$\\mathrm{I}^{\\!*}_7\\!(6)$', '$\\mathrm{I}^{\\!*}_8\\!(6)$', '$\\mathrm{IV}^{\\!*}$', '$\\mathrm{III}^{\\!*}\\!(3)$', '$\\mathrm{III}^{\\!*}\\!(5)$', '$\\mathrm{III}^{\\!*}\\!(7)$', '$\\mathrm{III}^{\\!*}\\!(8)$', '$\\mathrm{II}^{\\!*}\\!(3)$', '$\\mathrm{II}^{\\!*}\\!(4)$', '$\\mathrm{II}^{\\!*}\\!(6)$']
    tdnames = ['II(4)', 'II(6)', 'II(7)', 'I0', 'I1', 'I2', 'I3', 'III(3)', 'III(5)', 'III(7)', 'III(8)', 'IV', 'I0*(4)', 'I0*(5)', 'I0*(6)', 'I1*', 'I2*(4)', 'I2*(6)', 'I2*(7)', 'I3*(4)', 'I3*(5)', 'I4*(4)', 'I4*(6)', 'I5*(4)', 'I5*(6)', 'I6*(4)', 'I6*(6)', 'I7*(4)', 'I7*(6)', 'I8*(6)', 'IV*', 'III*(3)', 'III*(5)', 'III*(7)', 'III*(8)', 'II*(3)', 'II*(4)', 'II*(6)']
    sorted_type_data = []
    for name in tdnames:
        for td in type_data:
            tdname_tmp = td.name
            tdname_tmp = tdname_tmp.replace('2^12 ','')
            if tdname_tmp == name:
                sorted_type_data.append(td)
    type_data = list(sorted_type_data)
    #print([td.name for td in type_data])
    # ['II(4)', 'II(6)', 'II(7)', '2^12 I0', '2^12 I1', '2^12 I2', '2^12 I3', 'III(3)', 'III(5)', 'III(7)', 'III(8)', 'IV', 'I0*(4)', 'I0*(5)', 'I0*(6)', 'I1*', 'I2*(4)', 'I2*(6)', 'I2*(7)', 'I3*(4)', 'I3*(5)', 'I4*(4)', 'I4*(6)', 'I5*(4)', 'I5*(6)', 'I6*(4)', 'I6*(6)', 'I7*(4)', 'I7*(6)', 'I8*(6)', 'IV*', 'III*(3)', 'III*(5)', 'III*(7)', 'III*(8)', 'II*(3)', 'II*(4)', 'II*(6)']

    # Figure options
    if dpi is None:
        dpi = 400
    if text_dpi is None:
        text_dpi = dpi
    
    tick_spacing = int(modulus/8)
    tick_list = [tick_spacing*k for k in range(int(modulus/tick_spacing))]
    tick_list.append(tick_spacing*int(modulus/tick_spacing) - 1)    
    tick_names = [r'$' + r'\dots\!' + get_reversed_val(tmp, min_len=min_len, give_str=True)[5:] + '$' for tmp in tick_list]
    tick_names_a = [r'$\,\;' + tmp[1:] + '      ' for tmp in tick_names]
    tick_list = [tmp*pt_scale for tmp in tick_list] # do this after setting names
    a_label = r'$a$'
    b_label = '\n' + r'$\quad\quad\quad\quad b \ (2\text{-adic})$'
    axes_labels = [a_label, b_label]
    axes_fontsize_scale_factor = figsize/20.0 * 2#3.25 #* 2.0
    for td in type_data:
        td.font_scale_factor *= figsize/20.0

    # p^4|a and p^6 region setup
    if not region_pts_from_file:
        divisible_region = set()
        for a in range(int(modulus/8)):
            for b in range(int(modulus/8)):
                divisible_region.add((a,b))
    else:
        if polygon_plot:
            divisible_region = get_squares_from_file('divisible', p=2)
        else:
            divisible_region = get_region_pts_from_file('divisible', p=2)
    
    # Create region plots for each type
    for i, td in enumerate(type_data):        
        region_title = td.latex_name
        if hasattr(td, 'color'):
            color = td.color
        else:
            shift = i%2
            hval = i*(1 + sqrt(5.0))/2 % 1
            vval = shift*0.2 + 0.8
            sval = (1-shift)*0.2 + 0.8
            
            if td.symbol == 'I3':
                hval -= 0.18
            if td.symbol == 'I7*' and td.conductor == 2**6:
                hval -= 0.14
                vval = 1
                sval = 0.8
            if td.symbol == 'I4*' and td.conductor == 2**6:
                hval -= 0.02
                vval = 1
                sval = 1
            if td.symbol == 'III' and td.conductor == 2**5:
                hval -= 0.22
                vval = 1
                sval = 1
            if td.symbol == 'II' and td.conductor == 2**7:
                hval -= 0.08
                vval = 1
                sval = 1
            if td.symbol == 'III*' and td.conductor == 2**3:
                hval -= 0.04
                vval = 0.8
                sval = 1
            if td.symbol == 'IV*':
                hval -= 0.04
                vval = 1
                sval = 1
            if td.symbol == 'I2*' and td.conductor == 2**6:
                hval -= 0.02
                vval = 0.8
                sval = 1
            if td.symbol == 'I0*' and td.conductor == 2**4:
                hval -= 0.08
                vval = 0.8
            if td.symbol == 'I1*':
                hval = 0.85
                vval = 1
                sval = 0.6
            if td.symbol == 'II*' and td.conductor == 2**6:
                hval -= 0.05
            if td.symbol == 'III' and td.conductor == 2**8:
                hval -= 0.02
            if td.symbol == 'I3*' and td.conductor == 2**5:
                sval = 0.6
            
            if 0.64 < hval%1 and hval%1 < 0.72 and sval > 0.99:
                sval = 0.87
            color = hue(hval, v=vval, s=sval)
        
        if not region_pts_from_file:
            v = local_data[td.dict_key()]
            vred = get_vred(v, modulus)
            region, region_plot = get_region_plot(vred, min_len, p=2, polygon_plot=polygon_plot, pt_scale=pt_scale,
                                              dpi=dpi, xmin=0, xmax=modulus*pt_scale, ymin=0, ymax=modulus*pt_scale, color=color, aspect_ratio=1,
                                              axes_labels=axes_labels, ticks=[tick_list, tick_list], tick_formatter=[tick_names_a, tick_names])
            divisible_region.difference_update(set(region))
        else:
            if False:
                region = get_region_pts_from_file(td.name, p=2)
                region_plot = squareify(region, pt_scale=pt_scale, p=2,
                                                    dpi=dpi, xmin=0, xmax=modulus*pt_scale, ymin=0, ymax=modulus*pt_scale, color=color, aspect_ratio=1,
                                                    axes_labels=axes_labels, ticks=[tick_list, tick_list], tick_formatter=[tick_names_a, tick_names])
            else:
                squares = get_squares_from_file(td.name, p=2)
                polygon_list = squares_to_polygons(squares, pt_scale=pt_scale,
                                                   dpi=dpi, xmin=0, xmax=modulus*pt_scale, ymin=0, ymax=modulus*pt_scale, color=color, aspect_ratio=1,
                                                   axes_labels=axes_labels, ticks=[tick_list, tick_list], tick_formatter=[tick_names_a, tick_names])
                region_plot = sum(polygon_list)

        region_plot_list.append([region_plot, td])
        
        for loc in td.locations:
            if color == 'black':
                text_color = 'white'
            else:
                text_color = 'black'
            
            text_loc = vector(loc)
            if False: #td.modulus < modulus:
                text_loc = vector(loc) - vector((ZZ(1)/2, ZZ(1)/2))
            if False: #'2^12' in td.name:
                text_loc = vector(loc) + vector((ZZ(1)/4, 0))
            if td.modulus >= modulus:
                text_loc = vector(loc) + vector((ZZ(1)/2, ZZ(1)/2))
            
            if td.symbol == 'I0' and loc[0] < modulus/2:
                fontsize = 2*td.font_size()
            elif td.name == 'I4*(6)':
                fontsize = td.font_size() #* 2
            elif td.name == 'I4*(4)':
                fontsize = td.font_size()# / 2
            else:
                fontsize = td.font_size()
            
            if fontsize * figsize > 10:
                region_text = text(region_title, text_loc, color=text_color, xmin=0, xmax=modulus*pt_scale, ymin=0, ymax=modulus*pt_scale, fontsize=fontsize, dpi=text_dpi,
                                   axes_labels=axes_labels, ticks=[tick_list, tick_list], tick_formatter=[tick_names_a, tick_names])
                text_list.append(region_text)
    
    if False:
        divisible_region_plot = squareify(divisible_region, pt_scale=pt_scale, p=2,
                                            dpi=dpi, xmin=0, xmax=modulus*pt_scale, ymin=0, ymax=modulus*pt_scale, color=td_div.color, aspect_ratio=1,
                                            axes_labels=axes_labels, ticks=[tick_list, tick_list], tick_formatter=[tick_names_a, tick_names])
    else:
        polygon_list = squares_to_polygons(divisible_region, pt_scale=pt_scale,
                                           dpi=dpi, xmin=0, xmax=modulus*pt_scale, ymin=0, ymax=modulus*pt_scale, color=td_div.color, aspect_ratio=1,
                                           axes_labels=axes_labels, ticks=[tick_list, tick_list], tick_formatter=[tick_names_a, tick_names])
        divisible_region_plot = sum(polygon_list)
    region_plot_list.append([divisible_region_plot, td_div])
    
    for loc in td_div.locations:
        text_loc = vector(loc) - vector((ZZ(1)/2, ZZ(1)/2))
        text_loc *= modulus/ZZ(2)**10 * pt_scale
        fontsize = 1.3 * td_div.font_size() * figsize/20.0
        if True: #fontsize * figsize > 10:
            region_text = text(td_div.latex_name, text_loc, color='black', xmin=0, xmax=modulus*pt_scale, ymin=0, ymax=modulus*pt_scale, fontsize=fontsize, dpi=text_dpi,
                               axes_labels=axes_labels, ticks=[tick_list, tick_list], tick_formatter=[tick_names_a, tick_names])
            text_list.append(region_text)

    region_plot_list.sort(key = lambda tup: -tup[1].modulus)
    region_plot_list = [tmp[0] for tmp in region_plot_list]
    
    plt = sum(region_plot_list + text_list)
    
    plt.axes_width(plt.axes_width()*axes_width_scale_factor)
    plt.fontsize(plt.fontsize()*axes_fontsize_scale_factor)
    
    show(plt, figsize=figsize)
    return plt


def get_region_figure_3(region_pts_from_file=True, local_data=None, dpi=None, text_dpi=None, figsize=40, min_len=8, pt_scale=1, axes_width_scale_factor=1, random_text_loc=False, polygon_plot=True):
    if region_pts_from_file:
        assert(min_len == 8)
    if (not region_pts_from_file) and (local_data is None):
        local_data = load_local_data(3)
    modulus = 3**min_len
    
    region_plot_list = []
    text_list = []
    type_data = reductiontables.get_type_data_3(with_divisible=True)
    td_div = type_data[-1]
    assert(td_div.name == 'divisible')
    type_data = type_data[:-1]
    type_data = [td for td in type_data if td.modulus <= modulus]
    if region_pts_from_file:
        type_data = [td for td in type_data if type_name_to_region_pts_filename(td.name) in os.listdir('./')]
    else:
        type_data = [td for td in type_data if td.dict_key() in local_data]

    # Figure options
    if dpi is None:
        dpi = 400
    if text_dpi is None:
        text_dpi = dpi
    
    tick_spacing = int(modulus/9)
    tick_list = [tick_spacing*k for k in range(int(modulus/tick_spacing))]
    tick_list.append(tick_spacing*int(modulus/tick_spacing) - 1)
    tick_names = [r'$' + r'\dots\!' + get_reversed_val(tmp, min_len=min_len, give_str=True, p=3)[4:] + '$' for tmp in tick_list]
    tick_names_a = [r'$\,\;' + tmp[1:] + '      ' for tmp in tick_names]
    tick_list = [tmp*pt_scale for tmp in tick_list] # do this after setting names
    a_label = r'$a$'
    b_label = '\n' + r'$\quad\quad\quad\quad b \ (3\text{-adic})$'
    axes_labels = [a_label, b_label]
    axes_fontsize_scale_factor = figsize/40.0 * 6.5#* 4.0
    for td in type_data:
        td.font_scale_factor *= figsize/40.0
        
    # p^4|a and p^6 region setup
    if not region_pts_from_file:
        divisible_region = set()
        for a in range(int(modulus * ZZ(3)**(-4))):
            for b in range(int(modulus * ZZ(3)**(-6))):
                divisible_region.add((a,b))
    else:
        if polygon_plot:
            divisible_region = get_squares_from_file('divisible')
        else:
            divisible_region = get_region_pts_from_file('divisible')

    # Create region plots for each type
    for i, td in enumerate(type_data):        
        region_title = td.latex_name
        if hasattr(td, 'color'):
            color = td.color
        else:
            shift = i%2
            hval = i*(1 + sqrt(5.0))/2 % 1
            vval = shift*0.2 + 0.8
            sval = (1-shift)*0.2 + 0.8            
            if 0.64 < hval%1 and hval%1 < 0.72 and sval > 0.99:
                sval = 0.87
            if td.name == '3^12 I1':
                hval -= 0.13
            color = hue(hval, v=vval, s=sval)
        
        if not region_pts_from_file:
            v = local_data[td.dict_key()]
            vred = get_vred(v, modulus)
            region, region_plot = get_region_plot(vred, min_len, p=3, polygon_plot=polygon_plot, pt_scale=pt_scale,
                                              dpi=dpi, xmin=0, xmax=modulus*pt_scale, ymin=0, ymax=modulus*pt_scale, color=color, aspect_ratio=1,
                                              axes_labels=axes_labels, ticks=[tick_list, tick_list], tick_formatter=[tick_names_a, tick_names])
            divisible_region.difference_update(set(region))
        else:
            if False:
                region = get_region_pts_from_file(td.name)
                region_plot = squareify(region, pt_scale=pt_scale,
                                                    dpi=dpi, xmin=0, xmax=modulus*pt_scale, ymin=0, ymax=modulus*pt_scale, color=color, aspect_ratio=1,
                                                    axes_labels=axes_labels, ticks=[tick_list, tick_list], tick_formatter=[tick_names_a, tick_names])
            else:
                squares = get_squares_from_file(td.name)
                polygon_list = squares_to_polygons(squares, pt_scale=pt_scale,
                                                   dpi=dpi, xmin=0, xmax=modulus*pt_scale, ymin=0, ymax=modulus*pt_scale, color=color, aspect_ratio=1,
                                                   axes_labels=axes_labels, ticks=[tick_list, tick_list], tick_formatter=[tick_names_a, tick_names])
                region_plot = sum(polygon_list)

        region_plot_list.append([region_plot, td])

        if random_text_loc:
            region.sort()
            num_titles = 20
            loc_list = []
            for j in range(num_titles):
                index = random.randint(0,len(region)-1)
                loc_list.append(region[index])
        else:
            loc_list = td.locations
        
        for loc in loc_list:
            try:
                text_loc = vector(loc)
                if ('(' in td.latex_name) or ('_' in td.latex_name) or (r'\mid' in td.latex_name):
                    text_loc += vector((ZZ(1)/2, ZZ(1)/2))
                if td.name == '3^12 I1':
                    text_loc += vector((1, 0))
                text_loc *= pt_scale
            except TypeError:
                print(td.name)
                print(loc)
                raise
            if not random_text_loc:
                text_loc *= (modulus/ZZ(3)**8)
                        
            if color == 'black':
                text_color = 'white'
            else:
                text_color = 'black'
            
            #if td.symbol == 'I0' and loc[0] > modulus/3:
            #    fontsize = 3**6*td.font_size()
            if (td.symbol == 'I0*') and (loc[1] < modulus/9):
                if figsize > 20:
                    fontsize = 3 * td.font_size() / td.font_scale_factor
                else:
                    fontsize = [tdtmp.font_size() for tdtmp in type_data if tdtmp.name == 'IV(5)'][0]
            elif (td.symbol in ['II*']) and (figsize <= 10):
                fontsize = td.font_size() / 2
            elif (figsize <= 7) and (td.symbol in ['III*', 'IV*']):
                fontsize = td.font_size() * 0.8
            else:
                fontsize = td.font_size()
            #fontsize *= pt_scale

            if fontsize * figsize > 2:
                region_text = text(region_title, text_loc, color=text_color, xmin=0, xmax=modulus*pt_scale, ymin=0, ymax=modulus*pt_scale, fontsize=fontsize, dpi=text_dpi,
                                   axes_labels=axes_labels, ticks=[tick_list, tick_list], tick_formatter=[tick_names_a, tick_names])
                text_list.append(region_text)

    if False:
        divisible_region_plot = squareify(divisible_region, pt_scale=pt_scale,
                                            dpi=dpi, xmin=0, xmax=modulus*pt_scale, ymin=0, ymax=modulus*pt_scale, color=td_div.color, aspect_ratio=1,
                                            axes_labels=axes_labels, ticks=[tick_list, tick_list], tick_formatter=[tick_names_a, tick_names])
    else:
        polygon_list = squares_to_polygons(divisible_region, pt_scale=pt_scale,
                                           dpi=dpi, xmin=0, xmax=modulus*pt_scale, ymin=0, ymax=modulus*pt_scale, color=td_div.color, aspect_ratio=1,
                                           axes_labels=axes_labels, ticks=[tick_list, tick_list], tick_formatter=[tick_names_a, tick_names])
        divisible_region_plot = sum(polygon_list)
    region_plot_list.append([divisible_region_plot, td_div])
    
    for loc in td_div.locations:
        text_loc = vector(loc) + vector((ZZ(1)/2, ZZ(1)/2))
        text_loc *= modulus/ZZ(3)**8 * pt_scale
        if figsize > 10:
            region_text = text(td_div.latex_name, text_loc, color='black', xmin=0, xmax=modulus*pt_scale, ymin=0, ymax=modulus*pt_scale, fontsize=td_div.font_size(), dpi=text_dpi,
                               axes_labels=axes_labels, ticks=[tick_list, tick_list], tick_formatter=[tick_names_a, tick_names])
            text_list.append(region_text)
    
    region_plot_list.sort(key = lambda tup: -tup[1].modulus)
    region_plot_list = [tmp[0] for tmp in region_plot_list]
    
    plt = sum(region_plot_list + text_list)
    
    plt.axes_width(plt.axes_width()*axes_width_scale_factor)
    plt.fontsize(plt.fontsize()*axes_fontsize_scale_factor)
    
    show(plt, figsize=figsize)
    return plt


###

def squareify(point_list, p=3, pt_scale=1, **kwargs):
    square_list = points_to_squares(point_list, p=p)
    polygon_list = squares_to_polygons(square_list, pt_scale=pt_scale, **kwargs)
    return sum(polygon_list)

    
def points_to_squares(point_list, p=3):
    square_list = {coords:1 for coords in point_list}
    max_counter = 100
    counter = 0
    changed = True
    while changed and (counter < max_counter):
        square_list, changed = merge_squares(square_list, p=p)
        counter += 1
    if changed:
        print('Warning!')
    return square_list


def merge_squares(square_list, sort_list=True, p=3):
    changed = False
    square_list_coords = list(square_list.keys())
    if sort_list:
        square_list_coords.sort()
    for (x,y) in square_list_coords:
        if (x,y) in square_list:
            size = square_list[(x,y)]
            good = True
            good_squares = []
            for i in range(p):
                for j in range(p):
                    new_pt = (x + i*size, y + j*size)
                    present = new_pt in square_list
                    if not present:
                        good = False
                        break
                    right_size = square_list[new_pt] == size
                    if not right_size:
                        good = False
                        break
                    good_squares.append(new_pt)
                if not good:
                    break
            if good:
                new_centre = sum([vector(pt) for pt in good_squares]) / ZZ(len(good_squares))
                new_centre = (int(new_centre[0]), int(new_centre[1]))
                new_size = size * p
                for pt in good_squares:
                    del square_list[pt]
                square_list[new_centre] = new_size
                changed = True
    return [square_list, changed]


def squares_to_polygons(squares, pt_scale=1, **kwargs):
    # squares should be dict {centre:side_length}
    polygon_list = []
    for (x,y), size in squares.items():
        xm = x - size/2.0 + 1/2.0
        xp = x + size/2.0 + 1/2.0
        ym = y - size/2.0 + 1/2.0
        yp = y + size/2.0 + 1/2.0
        xm *= pt_scale
        xp *= pt_scale
        ym *= pt_scale
        yp *= pt_scale
        poly = polygon([(xm,ym), (xm,yp), (xp,yp), (xp,ym)], **kwargs)
        polygon_list.append(poly)
    to_ret = polygon_list
    #to_ret = sum(polygon_list)
    return to_ret


def plot_squares(square_list, figsize=20):
    squares_to_polys = {}
    for (x,y), size in square_list.items():
        poly = squares_to_polygons({(x,y):size}, dpi=400, aspect_ratio=1)[0]
        squares_to_polys[((x,y),size)] = poly
    fig_list = []
    centre_list = []
    for k,v in squares_to_polys.items():
        centre_list.append(k[0])
        centre_pt = points([k[0]], dpi=400, color='black', aspect_ratio=1)
        centre_text = text(str(k[0]), (k[0][0], k[0][1]), dpi=800, color='black', fontsize=2, aspect_ratio=1)
        fig_list.append(v)
        #fig_list.append(centre_pt)
        fig_list.append(centre_text)
    show(sum(fig_list), figsize=figsize)
    if len(centre_list) < 100:
        print(str(centre_list))
    return


            

        
#This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    if '-profile' in sys.argv:
        cProfile.run('main()')
    else:
        main()
