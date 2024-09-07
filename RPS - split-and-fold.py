import math
from re import T
import random
from sympy.ntheory.factor_ import totient 
import sympy
import numpy as np

def Findfmax(f):
    # Function for finding the maximum power-power divisor of f.
    L = list(sympy.factorint(f).items())
    f_max = 0
    for t in range(len(L)):
        if L[t][0] ** L[t][1] > f_max:
            f_max = L[t][0] ** L[t][1]
    return f_max

def findMSISdelta(B, n, m, deg, logq):
    # Function for estimating the MSIS hardness given parameters:
    # a (n x m) matrix in \Rq along with the solution bound B. It returns the
    # root Hermite factor \delta. We use the methodology presented by
    # [GamNgu08, Mic08] and hence m is irrelevant.
    if B >= 2 ** logq:                  # Check if the norm is above q
        return 2
    logB = math.log(B, 2)		
    logdelta = logB**2 / (4*n*deg*logq)
    return 2 ** logdelta

# Main function to compute norm
def SplitAndFoldWithNorm(f, ell, m, n, n_out, n1, n_batch, d, mu, logq, e, beta, exp_factor, inv_exp_factor,size_of_subtrset): 

    deg_f = totient(f)
    # Send the witness in the clear
    if m == 1:
        return [deg_f * math.log(8.6 * beta, 2),0]

    b_ip = 2 * beta / (math.sqrt(m) * deg_f ** (3/2))
    ell_0 = math.ceil( math.log(2 * beta **2 + 1, b_ip))
    b = math.ceil( (2 * beta + 1) ** (1 / ell)) 
    beta1 = exp_factor * (ell_0 + 1) * beta
    beta2 = exp_factor * ell * d * 1/2 * math.sqrt(m) * deg_f ** (3/2) * b

    beta2_prime = beta2                        # Due to the exact norm proof

    beta1_prime = 4 * math.sqrt(d) * inv_exp_factor * beta * beta2_prime
    beta0_prime = 8 * inv_exp_factor * beta1_prime

    if findMSISdelta(beta0_prime, 1 , m * deg_f, deg_f, logq) > 1.0044:
        return [1,1]

    # Norm proof size
    Norm_BDecomp = n_out * (ell-1) * deg_f * logq 
    Norm_Evaluations = 3  * (ell + 1)
    Norm_BatchRows = 0

    # Split and fold proof size
    new_n_out = n1 + n_batch
    SP_BDecomp = new_n_out * (ell-1) * deg_f * logq 
    SP_Split = ell * new_n_out * (d-1) * deg_f * logq
    SP_Fold = 0

    current_round_proof_size = (Norm_BDecomp + Norm_Evaluations + Norm_BatchRows) + (SP_BDecomp + SP_Split + SP_Fold)
    n = n + 3
    n_out = new_n_out
    m = m /d 

    # Knowledge error
    q = 2 ** logq
    current_KE = (2 * ell + 1)/size_of_subtrset + totient(f) * (4 * (ell + 1) + 2 * m) / (e * q ** e)

    beta = beta2

    P = SplitAndFoldWithNorm(f, ell, m, n, n_out, n1, n_batch, d, mu - 1, logq, e, beta, exp_factor, inv_exp_factor,size_of_subtrset)
    P[0] = P[0] + current_round_proof_size
    P[1] = P[1] + current_KE
    return P

# Security parameters
lmbda = 128                                     
LogROqueries = 64


### Setting global parameters
logq = 100                                       # Proof system modulus

# Set up the ring 
f = 27720
L = sympy.primefactors(f)
max_prime = L[len(L)-1] 
f_max = Findfmax(f)
phi = totient(f)


# Witness and statement dimensions
logm = 12
beta = 2 ** (logm / 2)                          # Starting norm


n = 1                                           # Number of tensor relations        
n_out = n                                       # Number of unstructured linear relations
n_batch = 1                                     # Batching



# Subtractive sets
if f == f_max:                                      # Prime power case
    exp_factor = max_prime
    inv_exp_factor = f / (2 * math.sqrt(2))
    size_of_subtrset = max_prime
else:                                               # Composite case
    exp_factor = 1
    inv_exp_factor = f / (4 * math.sqrt(2))
    size_of_subtrset = f / f_max

Params = []


for logd in sympy.divisors(logm):
    if logd == logm:
        break
    mu = logm / logd             
    m = 2 ** logm
    d = 2 ** logd
    for ell in range(2, logq):
    
        #Starting parameters
        n1 = 1
        n_batch = 1
        n = 1
        e = 1                       # very naive bounds, but it doesn't have impact on soundness
    
        P =  SplitAndFoldWithNorm(f, ell, m, n, n_out, n1, n_batch, d, mu, logq, e, beta, exp_factor, inv_exp_factor,size_of_subtrset)
        KE = P[1]
        if  P[1] < 1:
            repetition = math.floor(- lmbda / math.log(KE,2))
            Params.append([repetition * P[0] / 2 ** 23, d, ell, f, repetition, KE , size_of_subtrset, phi])


# Define the set of valid parameters


Params = np.array(Params)
Params = Params[np.argsort(Params[:,0])]
                
print("Log of the witness length: ",  math.log(phi * m,2))
print("Proof size in MB, followed by parameters: ", Params[0])
print("Commitment size: ", round(phi*logq / 2 ** 13, 2))
