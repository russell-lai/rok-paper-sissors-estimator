from dataclasses import *
from typing import List
from sage.all import euler_phi, Expression, var, mod, ceil, is_prime_power, sqrt, radical
from lattice_lib import *
from lattice_lib.util import *

@dataclass
class SubtractiveSet:
    """
    Data class for subtractive sets.
    
    Examples:
    
    sage: SubtractiveSet.gen_klno24_cyclotomic(9)
    SubtractiveSet(cardinality=3, gamma_2=3, gamma_inf=6, theta_2=9/4*sqrt(2), theta_inf=6)
    
    sage: SubtractiveSet.gen_klno24_cyclotomic(60) 
    SubtractiveSet(cardinality=12, gamma_2=1, gamma_inf=1024/pi^3, theta_2=15/2*sqrt(2), theta_inf=1024/pi^3)
    """
    cardinality: int = 2    # cardinality of the subtractive set C
    gamma_2: float = 1      # forward expansion factor of C in canonical ell_2-norm
    gamma_inf: float = 1    # forward expansion factor of C in coefficient ell_inf-norm
    theta_2: float = 1      # inverse expansion factor of C in canonical ell_2-norm
    theta_inf: float = 1    # inverse expansion factor of C in coefficient ell_inf-norm
    
    def gen_klno24_cyclotomic(f: int):
        """ 
        Use the subtractive set constructions for cyclotomic fields reported in KLNO24. 
        """
        if mod(f,4) == 2:
            raise Exception("Conductor f cannot be congruent to 2 modulo 4.")
        if is_prime_power(f):
            if f <= 4:
                raise Exception("Conductor f <= 4 is not supported.")
            phi = euler_phi(f)
            p = radical(f)
            return SubtractiveSet(cardinality = p, gamma_2 = p, gamma_inf = phi, theta_2 = f/(2*sqrt(2)), theta_inf = phi)
        else:
            phi = euler_phi(f)
            fmax = max_prime_power_divisor(f) 
            c_rad = (4/pi)**len(f.prime_divisors())
            return SubtractiveSet(cardinality = f/fmax, gamma_2 = 1, gamma_inf = c_rad * phi, theta_2 = f/(4*sqrt(2)), theta_inf = c_rad * phi)

@dataclass
class RingParam:
    """
    Data class for ring parameters.
    
    Example: 
    sage: RingParam(f=60,log_q=32)
    RingParam(f=60, phi=16, log_q=32, e=3, ring_exp_inf=16)
    """
    f: int = 1                                  # conductor
    phi: int = field(init=False)                # degree
    log_q: int = 64                             # log of modulus q (assumed prime)
    e: int | None = None                        # inertia degree of q, i.e. R_q splits into fields of cardinality q**e 
    C: InitVar[SubtractiveSet | None] = None    # subtractive set parameters
    ring_exp_inf: float = field(init=False)     # ring expansion factor in coefficient ell_inf-norm

    def __post_init__(self, C):
        if mod(self.f,4) == 2:
            raise Exception("Conductor f cannot be congruent to 2 modulo 4.")
        self.phi = euler_phi(self.f)
        self.ring_exp_inf = euler_phi(self.f)
        if self.C == None:
            self.C = SubtractiveSet()
        if self.e == None:
            self.e = ceil(80/self.log_q)

    def size_Rq(self):
        return self.phi * self.log_q

@dataclass
class Relation:
    """Data class for relations."""
    ncom: int = 1                   # module rank of commitment, with extractability guarantee
    nrel_uncompressed: int = 1      # number of uncompressed relations, no extractability guarantee  
    nrel_compressed: int = 1        # number of compressed relations, no extractability guarantee 
    wdim: int = 1                   # dimension of each witness 
    wbdl: int = 1                   # bundle size of witnesses 
    log_beta_wit_2: float = 1       # log of canonical ell_2-norm bound of witness
    log_beta_wit_inf: float = 1     # log of coefficient ell_inf-norm bound of witness
    log_beta_ext_2: float = 1       # log of canonical ell_2-norm bound of witness
    log_beta_ext_inf: float = 1     # log of coefficient ell_inf-norm bound of witness

# The identity function
var('x')
id_expr = x
id_f = id_expr.function(x)

@dataclass
class RoK:
    """Data class for reductions of knowledge."""
    name: str       # name of RoK
    f_nrel_uncompressed: Expression = id_f      # function for new nrel_uncompressed
    f_nrel_compressed: Expression = id_f        # function for new nrel_compressed
    f_wdim: Expression = id_f                   # function for new wdim
    f_wbdl: Expression = id_f                   # function for new wbdl
    f_log_beta_wit_2: Expression = id_f         # function for new log_beta_wit_2
    f_log_beta_wit_inf: Expression = id_f       # function for new log_beta_wit_inf
    f_log_beta_ext_2: Expression = id_f         # function for new log_beta_ext_2
    f_log_beta_ext_inf: Expression = id_f       # function for new log_beta_ext_inf

    # TODO: communication cost
    # TODO: soundness cost

pi_noop     = RoK(name="noop")      # the no-operation RoK for testing purposes
pi_bdecomp  = RoK(name="bdecomp")   # TODO: set parameters
pi_split    = RoK(name="split")     # TODO: set parameters
pi_fold     = RoK(name="fold")      # TODO: set parameters
pi_batch    = RoK(name="batch")     # TODO: set parameters
pi_norm     = RoK(name="norm")      # TODO: set parameters
pi_ip       = RoK(name="ip")        # TODO: set parameters
pi_aut      = RoK(name="aut")       # TODO: set parameters
pi_finish   = RoK(name="finish")    # TODO: set parameters

class Protocol:
    # A history is a list of states
    # A state is a tuple consisting of 
    # - a relation, 
    # - the name of the protocol used to reach this state, 
    # - the added communication cost, 
    # - the accumulated communication cost, 
    # - the added soundness error, 
    # - the accumulated soundness error, 
    # - the ell_2-canonical norm of the extracted witness (only set if pi_finish is run), 
    # - the ell_inf-coefficient norm of the extracted witness (only set if pi_finish is run). 
    
    def execute(pi: RoK, rel: Relation):
        """
        Execute RoK pi on relation rel.

        Example:
        sage: rel = Relation()
        sage: new_rel = execute(pi_noop, rel)
        """
        
        # TODO: Raise exception if any protocol is run after pi_finish.
        
        rel.nrel_uncompressed = pi.f_nrel_uncompressed(rel.nrel_uncompressed)
        rel.nrel_compressed = pi.f_nrel_compressed(rel.nrel_compressed)
        rel.wdim = pi.f_wdim(rel.wdim)
        rel.wbdl = pi.f_wbdl(rel.wbdl)
        rel.log_beta_wit_2 = pi.f_log_beta_wit_2(rel.log_beta_wit_2)
        rel.log_beta_wit_inf = pi.f_log_beta_wit_inf(rel.log_beta_wit_inf)
        rel.log_beta_ext_2 = pi.f_log_beta_ext_2(rel.log_beta_ext_2)
        rel.log_beta_ext_inf = pi.f_log_beta_ext_inf(rel.log_beta_ext_inf)
        
        # TODO: Keep track of communication cost. 
        # TODO: Keep track of soundness cost. 
        
        return rel      # TODO: return communication and soundness costs