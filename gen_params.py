from sage.all import log, ceil, floor, is_prime, is_prime_power, is_even, is_odd, euler_phi, multiplicative_order, radical, mod, ZZ, factor
from symbolic_variables import *
from ring_params import sym_ring_param
from util import *

def gen_order_1_prime(f, log_q = 128):
    # Input: conductor f, (optional) bit length log_q
    # Output: prime q of bit length in {log_q-1, log_q} such that multiplicative order of q modulo f is 1
    k = floor(2**log_q/f) + 1
    weight = log_q
    q = 2**log_q
    while q > 2**(log_q-1):
        k -= 1
        q = f*k + 1 # ensure that q = -1 mod f
        if q.is_prime():
            return q

def gen_order_2_prime(f, log_q = 64):
    # Input: conductor f, (optional) bit length log_q
    # Output: prime q of bit length in {log_q-1, log_q} such that multiplicative order of q modulo f is 2
    k = floor(2**log_q/f) + 1
    weight = log_q
    q = 2**log_q
    while q > 2**(log_q-1):
        k -= 1
        q = f*k - 1 # ensure that q = -1 mod f
        if q.is_prime():
            return q



def gen_ring_param(sym_mode = False, f = 2, log_q = 64, q = 0):
    # Input: conductor f, (optional) bit length log_q, (optional) modulus q
    # Output: 
    #   - ring degree phi 
    #   - inertia degree e of q 
    #   - Size C, expansion factors gamma_2 and gamma_inf, inverse expansion factors theta_2 and theta_inf of subtractive set
    # TODO: Costs are probably wrong. Fill in the correct costs.
    if sym_mode:
        return sym_ring_param
    if q == 0:
        if log_q < 100:
          q = gen_order_2_prime(f, log_q = log_q)
        else:
          q = gen_order_1_prime(f, log_q = log_q)
    assert is_prime(q)
    assert not q.divides(f)
    log_q = readable_log(q)
    phi = euler_phi(f)
    Rq_size = phi * ceil(log(q,2))
    e = multiplicative_order(mod(q, f))
    ring_exp_inf = phi
    f_hat = ZZ(f/2) if is_even(f) and is_odd(ZZ(f/2)) and f > 2 else f
    if is_prime_power(f_hat):
        p = radical(f_hat)
        C = p
        gamma_2 = phi
        theta_2 = phi
        gamma_inf = phi
        theta_inf = phi
    else:
        F = factor(f_hat)
        k = len(F)
        f_max = max([F[i][0]**F[i][1] for i in range(len(F))])
        C = f_hat/f_max
        gamma_2 = 1
        theta_2 = phi
        gamma_inf = 2
        theta_inf = phi
    if is_even(f):
      fhat = f/2
    else:
      fhat = f
    return {
        "f" : f,
        "fhat": fhat,
        "phi" : phi, 
        "q" : q, 
        "log_q": log_q,
        "Rq_size" : Rq_size,
        "e" : e, 
        "C" : C, 
        "gamma_2" : gamma_2, 
        "theta_2" : theta_2, 
        "gamma_inf" : gamma_inf, 
        "theta_inf" : theta_inf, 
        "ring_exp_inf" : ring_exp_inf
    }

# Estimating SIS hardness
def findMSISdelta(phi, n, log_beta_sis, log_q):
    # Function for estimating the MSIS hardness given parameters:
    # a (n x m) matrix in \Rq along with the solution bound B. It returns the
    # root Hermite factor \delta. We use the methodology presented by
    # [GamNgu08, Mic08] and hence m is irrelevant.
    if log_beta_sis >= log_q:                  # Check if the norm is above q
        return 2	
    log_delta = log_beta_sis**2 / (4*n*phi*log_q)
    return 2 ** log_delta

# Find commitment module rank such that SIS is hard in canonical-2 norm
def find_commit_module_rank(phi,log_beta_sis,log_q,rhf=1.0044):
    if log_beta_sis >= log_q:                  # Check if the norm is above q
        return 0
    n = ceil(log_beta_sis**2 / (4 * phi * log_q * log(rhf,2)))
    return n

sym_sis_param = {
    "n" : v_ntop,
    "log_beta_sis" : log(v_beta_sis)
}

def gen_sis_param(sym_mode = False, ring_param = sym_ring_param, log_beta_sis = 60,rhf=1.0044):
    if sym_mode:
        return sym_sis_param
    if log_beta_sis >= ring_param["log_q"]:                  # Check if the norm is above q
        return 0
    n = ceil(log_beta_sis**2 / (4 * ring_param["phi"] * ring_param["log_q"] * log(rhf,2)))
    return {
        "n" : n,
        "log_beta_sis" : log_beta_sis
    }

def gen_sis_param_aggresive(sym_mode = False, ring_param = sym_ring_param, log_beta_sis = 60,rhf=1.0044):
    if sym_mode:
        return sym_sis_param
    if log_beta_sis >= ring_param["log_q"]:                  # Check if the norm is above q
        return 0
    n = ceil((log_beta_sis - log(ring_param["phi"], 2))**2 / (4 * ring_param["phi"] * ring_param["log_q"] * log(rhf,2)))
    return {
        "n" : n,
        "log_beta_sis" : log_beta_sis
    }
    