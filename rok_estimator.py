from dataclasses import *
from typing import List, Tuple
from sage.all import euler_phi, Expression, var, mod, ceil, floor, is_prime_power, sqrt, radical
from lattice_lib import *
from lattice_lib.util import *
import importlib
estimator = importlib.import_module("lattice-estimator.estimator")
import warnings
from itertools import chain

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def _kb(v):
    """
    Convert bits to kilobytes.
    """
    return round(float(v / 8.0 / 1024.0), 1)

# bytes pretty-printing
UNITS_MAPPING = [
    # (1<<53, ' PB'),
    # (1<<43, ' TB'),
    # (1<<33, ' GB'),
    (1<<23, ' MB'),
    (1<<13, ' KB'),
    (1<<3, ' B'),
]

def pretty_size(bits, units=UNITS_MAPPING):
    """Get human-readable file sizes.
    simplified version of https://pypi.python.org/pypi/hurry.filesize/
    """
    for factor, suffix in units:
        if bits >= factor:
            break
    # amount = int(bits / factor)
    amount = n(bits/factor, digits=4)

    if isinstance(suffix, tuple):
        singular, multiple = suffix
        if amount == 1:
            suffix = singular
        else:
            suffix = multiple
    return str(amount) + suffix

# The following bounds are from https://eprint.iacr.org/2024/1972.pdf Corollary 1.
def bound_canon_2_from_coeff_inf(ring, beta_inf, dim = 1):
    return sqrt(ring.fhat * ring.phi * dim) * beta_inf 

def bound_log_canon_2_from_log_coeff_inf(ring, log_beta_inf, dim = 1):
    return log(sqrt(ring.fhat * ring.phi * dim), 2) + log_beta_inf 

def bound_coeff_inf_from_canon_2(beta_2):
    return beta_2

def bound_log_coeff_inf_from_log_canon_2(log_beta_2):
    return log_beta_2

@dataclass
class SubtractiveSet:
    """
    Data class for subtractive sets.
    
    Examples:
    
    sage: SubtractiveSet.gen_klno24_cyclotomic(27)
    SubtractiveSet(cardinality=3, gamma_2=3, theta_2=27/4*sqrt(2), gamma_inf=18, theta_inf=18)
    
    sage: SubtractiveSet.gen_klno24_cyclotomic(60) 
    SubtractiveSet(cardinality=12, gamma_2=1, theta_2=15/2*sqrt(2), gamma_inf=1024/pi^3, theta_inf=1024/pi^3)
    """
    cardinality: int = 2    # cardinality of the subtractive set C
    gamma_2: float = 1      # forward expansion factor of C in canonical ell_2-norm
    theta_2: float = 1      # inverse expansion factor of C in canonical ell_2-norm
    gamma_inf: float = 1    # forward expansion factor of C in coefficient ell_inf-norm
    theta_inf: float = 1    # inverse expansion factor of C in coefficient ell_inf-norm
    
    def gen_klno24_cyclotomic(f: int):
        """ 
        Use the subtractive set constructions for cyclotomic fields reported in KLNO24. 
        """
        if mod(f,4) == 2:
            raise Exception("Conductor f cannot be congruent to 2 modulo 4.")
        # Hard-code values for empirically verified rings
        if f == 24:
            return SubtractiveSet(cardinality = 3, gamma_2 = 1, theta_2 = f/(4*sqrt(2)), gamma_inf = 2, theta_inf = 7)
        if f == 60:
            return SubtractiveSet(cardinality = 12, gamma_2 = 1, theta_2 = f/(4*sqrt(2)), gamma_inf = 4, theta_inf = 17)
        # General case
        if is_prime_power(f):
            if f <= 4:
                raise Exception("Conductor f <= 4 is not supported.")
            phi = euler_phi(f)
            p = radical(f)
            return SubtractiveSet(cardinality = p, gamma_2 = p, theta_2 = f/(2*sqrt(2)), gamma_inf = phi, theta_inf = phi)
        else:
            phi = euler_phi(f)
            fmax = max_prime_power_divisor(f) 
            c_rad = (4/pi)**len(f.prime_divisors())
            return SubtractiveSet(cardinality = f/fmax, gamma_2 = 1, theta_2 = f/(4*sqrt(2)), gamma_inf = c_rad * phi, theta_inf = c_rad * phi)
        
    def show(self):
        print(f'Subtractive set parameters:')
        print(f'    cardinality: {self.cardinality}')
        print(f'    forward ell_2 expansion factor gamma_2: 2^{ceil(log(self.gamma_2,2))}')
        print(f'    inverse ell_2 expansion factor theta_2: 2^{ceil(log(self.theta_2,2))}')
        print(f'    forward ell_inf expansion factor gamma_inf: 2^{ceil(log(self.gamma_inf,2))}')
        print(f'    inverse ell_inf expansion factor theta_inf: 2^{ceil(log(self.theta_inf,2))}')

@dataclass
class RingParam:
    """
    Data class for ring parameters.
    
    Example: 
    sage: RingParam(f=60,log_q=32)
    RingParam(f=60, phi=16, log_q=32, e=3, ring_exp_inf=16)
    """
    f: int = 1                                  # conductor
    fhat: int = 1                               # fhat = f if f is odd and f/2 if f is even
    phi: int = field(init=False)                # degree
    log_beta_sis_2: int = 32                    # log of ell_2-norm bound beta_sis_2 for SIS problem
    log_beta_sis_inf: int = 32                  # log of ell_inf-norm bound beta_sis_inf for SIS problem TODO: Currently ignored
    log_q: int = 64                             # log of modulus q (assumed prime)
    residue_deg: int | None = None              # inertia degree of q, i.e. R_q splits into fields of cardinality q**residue_deg 
    C: InitVar[SubtractiveSet | None] = None    # subtractive set parameters
    ring_exp_inf: float = field(init=False)     # ring expansion factor in coefficient ell_inf-norm
    n_sis: int | None = None                    # module rank of SIS instance
    secpar_target: int = 128                    # target SIS security
    secpar_result: int | None = None            # resulting SIS security
    kappa_target: int = 80                      # target Schwartz-Zippel security
    kappa_result: int | None = None             # resulting Schwartz-Zippel security
    kappa_hedge: int = 3                       # hedge against running the RoKs 2**kappa_hedge times

    def __post_init__(self, C):
        if mod(self.f,4) == 2:
            raise Exception("Conductor f cannot be congruent to 2 modulo 4.")
        if self.log_beta_sis_2 >= self.log_q:
            raise Exception("Norm bound beta_sis must be smaller than modulus q.")
        self.phi = euler_phi(self.f)
        self.ring_exp_inf = euler_phi(self.f)
        if self.C == None:
            self.C = SubtractiveSet.gen_klno24_cyclotomic(self.f)       # Use the subtractive set construction for cyclotomic fields reported in KLNO24
        if self.residue_deg == None:
            self.residue_deg = ceil((self.kappa_target + self.kappa_hedge)/self.log_q)       # Aim for kappa_target bits of soundness for Schwartz-Zippel. Adding kappa_hedge bits to hedge against running the RoKs 2**kappa_hedge times
            self.kappa_result = self.log_q * self.residue_deg
        if is_even(self.f):
            self.fhat = self.f/2
        else:
            self.fhat = self.f
        # Estimate vSIS security
        # Heuristics 1: vSIS is as hard as SIS
        # Heuristics 2: SIS over modules in canonical ell_2-norm is as hard as SIS over ZZ in ell_2-norm
        # TODO: Also take ell_inf-norm into account
        if self.n_sis == None:    
            for n_sis in range(1,2048):
                sis = estimator.SIS.Parameters(self.phi*n_sis, 2**self.log_q, 2**self.log_beta_sis_2)  
                with HiddenPrints():
                    costs = estimator.SIS.estimate(sis) # BUG: seems that estimator.SIS.estimate returns +Infinity when security is too low
                sec = min(cost["rop"] for cost in costs.values())
                if sec < oo and sec >= 2**self.secpar_target:    
                    self.n_sis = n_sis
                    self.secpar_result = floor(log(sec,2))
                    break
            if self.n_sis == None:
                raise Exception("All SIS module rank considered are too small for the target security level.")
        else:   
            sis = estimator.SIS.Parameters(self.phi*self.n_sis, 2**self.log_q, 2**self.log_beta_sis_2)
            with HiddenPrints():
                costs = estimator.SIS.estimate(sis)
            sec = min(cost["rop"] for cost in costs.values())
            if sec < 2**self.secpar_target:
                warnings.warn("Specified module rank for SIS is too small for the target security level.")
                # raise Exception("Specified module rank for SIS is too small for the target security level.")

    def size_Rq(self):
        return self.phi * self.log_q
    
    def show(self):
        print("Ring parameters:")
        print(f"    conductor f: {self.f}, degree phi: {self.phi}, modulus q: 2^{self.log_q}, beta_sis_2: 2^{self.log_beta_sis_2}")
        print(f"    SIS module rank n_sis: {self.n_sis}, target SIS security: {self.secpar_target}, resulting SIS security: {self.secpar_result}")
        print(f"    residue degree: {self.residue_deg}, target Schwartz-Zippel security: {self.kappa_target}, resulting Schwartz-Zippel security: {self.kappa_result}")
        print(f"    |R_q| = {pretty_size(self.size_Rq())}, |R_q^(n_sis)| = {pretty_size(self.size_Rq()*self.n_sis)}")
        print(f' ')
        self.C.show()

@dataclass
class Relation:
    """
    Data class for relations. A Relation object contains methods modelling reductions of knowledge (RoK) each of which returns a reduced relation of the RoK. 
    
    Example:
    sage: rel = Relation(ring=RingParam(f=60,n_sis=2),wdim=60)
    sage: rel_bdecomp = rel.pi_bdecomp(ell=2)
    sage: rel_bdecomp.show()
    
    Example when self if a Relation object, useful when defining new RoKs:
    sage: rel_params = {
            "ring": self.ring,
            "trivial": self.trivial,
            "n_compress": self.n_compress,
            "n_commit": self.n_commit,
            "n_rel": self.n_rel,
            "wdim": self.wdim,
            "rep": self.rep,
            "log_beta_wit_2": self.log_beta_wit_2,
            "log_beta_wit_inf": self.log_beta_wit_inf
        }
    sage: replace(self, **rel_params)
    """
    ring: RingParam = field(repr=False)             # ring parameters
    trivial : bool = False                          # True if the relation is the "True" relation
    n_compress: int | None = None                   # number of compressed relations, including both commitment and non-commitment relations
    n_commit: int | None = None                     # module rank of commitment, will be overwritten by self.ring.n_sis.
    n_rel: int = 0                                  # number of (non-commitment) relations  
    wdim: int = 1                                   # dimension of each witness 
    rep: int = 1                                    # bundle size of witnesses 
    log_beta_wit_2: float | None = None             # log of canonical ell_2-norm bound of witness
    log_beta_wit_inf: float | None = None           # log of coefficient ell_inf-norm bound of witness
    log_beta_ext_2: float | None = None             # log of canonical ell_2-norm bound of extracted witness, only computed after extraction
    log_beta_ext_inf: float | None = None           # log of coefficient ell_inf-norm bound of extracted witness, only computed after extraction
    comm : int = 0                                  # communication cost spent to arrive at this relation
    acc_comm: int = 0                               # accumulated communication cost
    snd_err : int = 0                               # soundness error spent to arrive at this relation
    acc_snd_err: int = 0                            # accumulated soundness error
    log_beta_ext_2_func     : function = lambda x : x   # function mapping old canonical ell_2-norm of extracted witness to new one 
    log_beta_ext_inf_func   : function = lambda x : x   # function mapping old coefficient ell_inf-norm of extracted witness to new one
    
    def __post_init__(self):
        self.n_commit = self.ring.n_sis
        if self.n_compress == None:
            self.n_compress = self.ring.n_sis + self.n_rel # Do not assume compression at the beginning
        if self.n_compress < self.n_commit:
            raise Exception("Number of compressed relations cannot be smaller than number of commitment relations.")
        if self.log_beta_wit_2 == None and self.log_beta_wit_inf == None:
            raise Exception("Relation must have either a canonical ell_2-norm bound or a coefficient ell_inf-norm bound.")
        if self.log_beta_wit_2 == None:
            # self.log_beta_wit_2 = ceil(log(sqrt(self.wdim * self.rep * self.ring.phi * self.ring.fhat) * 2**self.log_beta_wit_inf,2)) # Measured in Frobenius norm
            self.log_beta_wit_2 = ceil(log(sqrt(self.wdim * self.ring.phi * self.ring.fhat) * 2**self.log_beta_wit_inf,2)) # Measured in max ell_2-norm over all columns
        if self.log_beta_wit_inf == None:
            self.log_beta_wit_inf = self.log_beta_wit_2 # If ell_inf-norm is not specified, trivially bound it by ell_2-norm 
        
    def wit_size(self):
        return self.wdim * self.rep * self.ring.phi * ceil(self.log_beta_wit_inf + 1)
    
    def show(self,label=None,brief=False):
        label_str = f'{label:8s}' if label else 'Relation'
        flag_log_beta_wit_2 = f'*' if self.log_beta_wit_2 + 1 > self.ring.log_beta_sis_2 else ' '                                   # NOTE: Underestimating security when log_beta_wit_2 is measured in Frobenius norm 
        flag_log_beta_ext_2 = f'*' if self.log_beta_ext_2 != None and self.log_beta_ext_2 + 1> self.ring.log_beta_sis_2 else ' '    # NOTE: Underestimating security when log_beta_ext_2 is measured in Frobenius norm 
        if self.snd_err == 0:
            log_snd_err = -oo
        else:
            log_snd_err = ceil(log(self.snd_err,2))
        if self.acc_snd_err == 0:
            log_acc_snd_err = -oo
        else:
            log_acc_snd_err = ceil(log(self.acc_snd_err,2))    
        if self.trivial:
            print(f'{label_str}: communication (growth | total) = ({pretty_size(self.comm):8s} | {pretty_size(self.acc_comm):8s}), soundness error (growth | total) = (2^{log_snd_err} | 2^{log_acc_snd_err})')
        elif brief:
            print(f'{label_str}: wdim = {self.wdim:8d}, rep = {self.rep:3d}, log_2-norm (real | extr) = ({ceil(self.log_beta_wit_2):3d}{flag_log_beta_wit_2}|{ceil(self.log_beta_ext_2):3d}{flag_log_beta_ext_2}), log_inf-norm (real | extr) = ({ceil(self.log_beta_wit_inf):3d} |{ceil(self.log_beta_ext_inf):3d} ), wit size = {pretty_size(self.wit_size()):8s}, communication (growth | total) = ({pretty_size(self.comm):8s} | {pretty_size(self.acc_comm):8s}), soundness error (growth | total) = (2^{log_snd_err} | 2^{log_acc_snd_err})')
        else:
            print(f'Relation:')
            print(f'    H * F * W = Y')
            print(f'Statement:')
            print(f'    H: n_compress x (n_commit + n_rel)')
            print(f'    F: (n_commit + n_rel) x wdim')
            print(f'    Y: n_compress x rep')
            print(f'Witness:')
            print(f'    W: wdim x rep')
            print(f'    ||sigma(W)||_2 <= 2^log_beta_wit_2 measured in max column canonical ell_2-norm')
            print(f'    ||psi(W)||_inf <= 2^log_beta_wit_inf')
            print(f'Parameters:')
            print(f'    n_compress = {self.n_compress}, n_commit = {self.n_commit}, n_rel = {self.n_rel}')
            print(f'    wdim = {self.wdim}, rep = {self.rep}')
            print(f'    log_2-norm (real | extr) = ({ceil(self.log_beta_wit_2)}{flag_log_beta_wit_2}|{ceil(self.log_beta_ext_2)}{flag_log_beta_ext_2}), log_inf-norm (real | extr) = ({ceil(self.log_beta_wit_inf)}|{ceil(self.log_beta_ext_inf)})')
            print(f'    wit size = {pretty_size(self.wit_size()):8s}')
            print(f' ')
            
    def execute(self, op, **kwargs):
        match op:
            case "noop":
                return self.pi_noop()
            case "finish":
                return self.pi_finish()
            case "bdecomp":
                return self.pi_bdecomp(**kwargs)
            case "split":
                return self.pi_split(**kwargs)
            case "fold":
                return self.pi_fold(**kwargs)
            case "batch":
                return self.pi_batch()
            case "norm":
                return self.pi_norm()
            case "ip":
                return self.pi_ip()
            case "aut":
                return self.pi_aut()
        
    def pi_noop(self):
        """
        Returns the relation resulting from the pi_noop RoK. 
        
        The RoK pi_noop does nothing, it reduces any relation to itself. It is there for testing purposes.
        """
        return deepcopy(self)
    
    def pi_finish(self):
        """
        Returns the relation resulting from the pi_finish RoK. 
        
        The pi_finish RoK reduces any relation to True.
        """
        comm = self.wit_size()
        rel_param = {
            # "ring": self.ring,
            "trivial": True,
            # "n_compress": self.n_compress,
            # "n_commit": self.n_commit,
            # "n_rel": self.n_rel,
            # "wdim": self.wdim,
            # "rep": self.rep,
            # "log_beta_wit_2": self.log_beta_wit_2,
            # "log_beta_wit_inf": self.log_beta_wit_inf
            "comm": comm,
            "snd_err": 0,
            "acc_comm" : self.acc_comm + comm,
            "log_beta_ext_2_func" : lambda x : self.log_beta_wit_2, # perfect extraction
            "log_beta_ext_inf_func" : lambda x : self.log_beta_wit_inf, # perfect extraction
        }
        return replace(self, **rel_param)
      
    def pi_bdecomp(self,base: int | None = None,ell: int | None = None):
        """
        Returns the relation resulting from the pi_bdecomp RoK. 
        
        Parameters: Specify either the base 'base' or the target number of chunks 'ell'. 
        
        The pi_bdecomp increases the bundle size 'rep' and decreases the norm. 
        """
        if base == None and ell == None:
            raise Exception("Parameters base and ell cannot be both undefined.")
        if base != None and ell != None:
            raise Exception("Parameters base and ell cannot be both defined.")
        if base != None:
            ell = ceil(log(2*(2**self.log_beta_wit_inf)+1,base))
        else:
            base = ceil( (2*(2**self.log_beta_wit_inf)+1)**(1/ell) ) 
            if base <= 1:
                raise Exception("The choice of ell is impossible. Parameter base must be greater than 1.")
        comm = self.ring.size_Rq() * (ell-1) * self.n_compress * self.rep
        rel_param = {
            # "ring": self.ring,
            # "trivial": self.trivial,
            # "n_compress": self.n_compress,
            # "n_commit": self.n_commit,
            # "n_rel": self.n_rel,
            # "wdim": self.wdim,
            "rep" : self.rep * ell,
            # "log_beta_wit_2" : bound_log_canon_2_from_log_coeff_inf(self.ring, log(floor(base / 2),2), dim=self.wdim * self.rep * ell), # measured in Frobenius norm
            "log_beta_wit_2" : bound_log_canon_2_from_log_coeff_inf(self.ring, log(floor(base / 2),2), dim=self.wdim), # measured in max ell_2-norm over all columns
            "log_beta_wit_inf" : log(floor(base / 2),2),
            "comm": comm,
            "snd_err": 0,
            "acc_comm": self.acc_comm + comm,
            "log_beta_ext_2_func" : lambda x : x + log((base**ell-1)/(base-1),2),
            "log_beta_ext_inf_func" : lambda x : x + log((base**ell-1)/(base-1),2),
        }        
        return replace(self, **rel_param)
    
    def pi_split(self,d: int):
        """     
        Returns the relation resulting from the pi_split RoK.
        
        Parameters: Target number of chunks 'd'.
        """
        if not d.divides(self.wdim):
            raise Exception("Cannot split the witness into d chunks. Parameter d must divide wdim.")
        comm = self.ring.size_Rq() * ((d-1) * self.n_commit + (d**2-1) * (self.n_compress-self.n_commit)) * self.rep
        snd_err = (d-1) / 2**(self.ring.log_q * self.ring.residue_deg)
        rel_params = {
            # "ring": self.ring,
            # "trivial": self.trivial,
            # "n_compress": self.n_compress,
            # "n_commit": self.n_commit,
            # "n_rel": self.n_rel,
            "wdim": ZZ(self.wdim / d),
            "rep": self.rep * d,
            # "log_beta_wit_2": self.log_beta_wit_2,
            # "log_beta_wit_inf": self.log_beta_wit_inf
            "comm": comm,
            "acc_comm": self.acc_comm + comm,
            "snd_err": snd_err,
            "acc_snd_err": self.acc_snd_err + snd_err,
        }        
        return replace(self, **rel_params)
    
    def pi_fold(self, repout: int | None = None):
        """
        Returns the relation resulting from the pi_fold RoK. 
        
        Fold witness matrix W and challenge matrix C into W' = W * C. 
        
        Parameters: Output bundle size 'repout'.
        """
        repin = self.rep
        if repout == None:
            # Ensure that repin / (self.ring.C.cardinality**repout) <= 2^-kappa_target  
            repout = ceil((self.ring.kappa_target + log(repin,2) + self.ring.kappa_hedge) / log(self.ring.C.cardinality,2)) # +kappa_hedge hedges against running the RoK 2**kappa_hedge times
        snd_err = repin / (self.ring.C.cardinality**repout)
        rel_params = {
            # "ring": self.ring,
            # "trivial": self.trivial,
            # "n_compress": self.n_compress,
            # "n_commit": self.n_commit,
            # "n_rel": self.n_rel,
            # "wdim": self.wdim,
            "rep": repout,
            # "log_beta_wit_2": log(sqrt(repout) * repin * self.ring.C.gamma_2,2) + self.log_beta_wit_2, # Measured in Frobenius norm
            "log_beta_wit_2": log(repin * self.ring.C.gamma_2,2) + self.log_beta_wit_2, # Measured in max ell_2-norm over all columns
            "log_beta_wit_inf": log(repin * self.ring.C.gamma_inf,2) + self.log_beta_wit_inf, 
            "snd_err": snd_err,
            "acc_snd_err": self.acc_snd_err + snd_err,
            "log_beta_ext_2_func" : lambda x : x + log(2 * sqrt(repin) * self.ring.C.theta_2,2),
            "log_beta_ext_inf_func" : lambda x : x + log(2 * self.ring.C.theta_inf,2),
        }
        return replace(self, **rel_params)
    
    def pi_batch(self):
        """
        Returns the relation resulting from the pi_batch RoK. 
        """
        # Batching only when n_compress > n_commit.
        # Otherwise, there is an unintuitive indirect cost for pi_batch: It increases n_compress by 1 even if n_compress = n_commit (i.e. n_rel = 0), and in the the communication cost of pi_split n_compress is multiplied by d**2 instead of d.
        if self.n_compress > self.n_commit: 
            snd_err = self.rep * self.n_rel / (2**(self.ring.log_q * self.ring.residue_deg))
            rel_params = {
                # "ring": self.ring,
                # "trivial": self.trivial,
                "n_compress": self.n_commit + 1, # TODO: Allow batching into more than 1 row to allow smaller field size.
                # "n_commit": self.n_commit,
                # "n_rel": self.n_rel,
                # "wdim": self.wdim,
                # "rep": self.rep,
                # "log_beta_wit_2": self.log_beta_wit_2,
                # "log_beta_wit_inf": self.log_beta_wit_inf
                "snd_err": snd_err,
                "acc_snd_err": self.acc_snd_err + snd_err,
            }
            return replace(self, **rel_params)
        else:
            return deepcopy(self)
    
    def pi_norm(self):
        """
        Returns the relation resulting from the pi_norm RoK. 
        """
        base = 2 * 2**self.log_beta_wit_inf + 1  
        ell = ceil(log( self.ring.ring_exp_inf * self.wdim * 2**(self.log_beta_wit_inf * 2), base )) 
        # base = 2 * beta_wit_inf + 1
        # ell >= log( ring_exp_inf * wdim * beta_wit_inf^2, base ) so that base^ell >= ring_exp_inf * wdim * beta_wit_inf^2
        comm = self.ring.size_Rq() * (ell * (self.ring.n_sis + self.n_rel) + 3 * self.rep + 3 * ell)
        snd_err = ell / 2**(self.ring.log_q * self.ring.residue_deg)
        rel_params = {
            # "ring": self.ring,
            # "trivial": self.trivial,
            "n_compress": self.n_compress + 3,
            # "n_commit": self.n_commit,
            "n_rel": self.n_rel + 3,
            # "wdim": self.wdim,
            "rep": self.rep + ell,
            # "log_beta_wit_2": bound_log_canon_2_from_log_coeff_inf(self.ring,self.log_beta_wit_inf, dim=self.wdim * (self.rep + ell)), # Measured in Frobenius norm
            "log_beta_wit_2": max([self.log_beta_wit_2, bound_log_canon_2_from_log_coeff_inf(self.ring,self.log_beta_wit_inf, dim=self.wdim)]), # Measured in max ell_2-norm over all columns
            "log_beta_wit_inf": self.log_beta_wit_inf,
            "comm" : comm,
            "acc_comm" : self.acc_comm + comm,
            "snd_err": snd_err,
            "acc_snd_err": self.acc_snd_err + snd_err,
            "log_beta_ext_2_func" : lambda x : self.log_beta_wit_2, # perfect extraction
            "log_beta_ext_inf_func" : lambda x : bound_coeff_inf_from_canon_2(self.log_beta_wit_2), # extraction of ell_2-norm is perfect, then bound ell_inf-norm by norm conversion
        }
        return replace(self, **rel_params)
    
    def pi_ip(self): # TODO: Dummy protocol to be implemented
        """
        Returns the relation resulting from the pi_ip RoK. 
        """
        return deepcopy(self)

    
    def pi_aut(self): # TODO: Dummy protocol to be implemented
        """
        Returns the relation resulting from the pi_aut RoK. 
        """
        return deepcopy(self)

@dataclass
class Simulation:
    ring_params: dict = field(repr=False) # ring parameters
    rel_params: dict = field(repr=False) # relation parameters 
    
    ring: RingParam = field(repr=False,init=False)      # ring parameters
    trace : List[Tuple[str, Relation]] = field(repr=False,init=False) # execution trace
    max_log_beta_wit_2 : int = 0           # maximum log_beta_wit_2
    max_log_beta_wit_inf : int = 0         # maximum log_beta_wit_inf
    max_log_beta_ext_2 : int = 0           # maximum log_beta_ext_2
    max_log_beta_ext_inf : int = 0         # maximum log_beta_ext_inf
    error_log : List[str] = field(repr=False,init=False) # error log
    
    def __post_init__(self):
        self.ring = RingParam(**self.ring_params)
        rel = Relation(ring = self.ring, **self.rel_params)
        self.trace = [("init", rel)]
        self.max_log_beta_wit_2 = rel.log_beta_wit_2
        self.max_log_beta_wit_inf = rel.log_beta_wit_inf
        self.error_log = []
    
    def execute(self, ops):
        """
        Simulates the execution of a sequence of RoKs on a relation.
        """
        # Forward direction, a.k.a. "correctness direction"
        for op, params in ops:
            new_rel = self.trace[-1][1].execute(op, **params)
            self.trace += [(op, new_rel)]
            if new_rel.log_beta_wit_2 > self.max_log_beta_wit_2:
                self.max_log_beta_wit_2 = new_rel.log_beta_wit_2
            if new_rel.log_beta_wit_inf > self.max_log_beta_wit_inf:
                self.max_log_beta_wit_inf = new_rel.log_beta_wit_inf
            
    def extract(self):
        # Backward direction, a.k.a. "extraction direction"        
        if self.trace[-1][0] == "finish":
            for i in range(len(self.trace)-1):
                # The RoK is from `rel_src` to `rel_tgt`. The extracted norm function is stored in `rel_tgt`.
                rel_tgt = self.trace[-i-1][1]
                rel_src = self.trace[-i-2][1]
            
                rel_src.log_beta_ext_2 = rel_tgt.log_beta_ext_2_func(rel_tgt.log_beta_ext_2)   
                rel_src.log_beta_ext_inf = rel_tgt.log_beta_ext_inf_func(rel_tgt.log_beta_ext_inf) 
                
                if self.trace[-i-1][0] == "norm":
                    if rel_src.ring.log_q - 1 <= rel_src.log_beta_ext_inf * 2 + log(rel_src.ring.ring_exp_inf,2) + log(rel_src.wdim,2):
                        self.error_log += [f"Extraction failure for pi_norm: The norm of the square of the extracted witness is 2^{ceil(rel_src.log_beta_ext_inf * 2 + log(rel_src.ring.ring_exp_inf,2) + log(rel_src.wdim,2))} overflowing modulo q."]
                    
                # Check if any norm is overestimated. 
                ring = rel_src.ring
                wdim = rel_src.wdim
                rep = rel_src.rep
                if rel_src.log_beta_ext_2 > bound_log_canon_2_from_log_coeff_inf(ring, rel_src.log_beta_ext_inf, dim=wdim * rep):
                    rel_src.log_beta_ext_2 = bound_log_canon_2_from_log_coeff_inf(ring, rel_src.log_beta_ext_inf, dim=wdim * rep)
                    # print(f"{self.trace[-i-2][0]}: ell-2 norm is overestimated!")
                    
                if rel_src.log_beta_ext_inf > bound_log_coeff_inf_from_log_canon_2(rel_src.log_beta_ext_2):
                    rel_src.log_beta_ext_inf = bound_log_coeff_inf_from_log_canon_2(rel_src.log_beta_ext_2)
                    # print(f"{self.trace[-i-2][0]}: ell-inf norm is overestimated!")
                    
                # Record maximum log_beta_ext_2 and log_beta_ext_inf
                if rel_src.log_beta_ext_2 > self.max_log_beta_ext_2:
                    self.max_log_beta_ext_2 = rel_src.log_beta_ext_2
                if rel_src.log_beta_ext_inf > self.max_log_beta_ext_inf:
                    self.max_log_beta_ext_inf = rel_src.log_beta_ext_inf
    
    def show_trace(self,brief=False):
        print(f'Execution Trace:')
        for op, rel in self.trace:
            rel.show(label=op,brief=True)
        print(f' ')
        
    def show(self,brief=False):
        self.ring.show()
        print(f' ')
        
        self.trace[0][1].show(brief=brief)
        
        if not brief:
            self.show_trace(brief=brief)
        
        total_comm = self.trace[-1][1].acc_comm
        total_snd_err = self.trace[-1][1].acc_snd_err
        print(f'Total Cost: communication = {pretty_size(total_comm):8s}, soundness error = 2^{ceil(log(total_snd_err,2))}')
        flag_log_beta_wit_2 = f'*' if self.max_log_beta_wit_2 + 1 > self.ring.log_beta_sis_2 else ' '                                   # NOTE: Underestimating security when log_beta_wit_2 is measured in Frobenius norm 
        flag_log_beta_ext_2 = f'*' if self.max_log_beta_ext_2 != None and self.max_log_beta_ext_2 + 1> self.ring.log_beta_sis_2 else ' '    # NOTE: Underestimating security when log_beta_ext_2 is measured in Frobenius norm 
        print(f'Maximum log ell_2-norm (real | extr) = ({ceil(self.max_log_beta_wit_2):3d}{flag_log_beta_wit_2}|{ceil(self.max_log_beta_ext_2):3d}{flag_log_beta_ext_2}), log SIS norm bound = {self.ring.log_beta_sis_2}')
        for err in self.error_log:
            print(f'{err}')
        print(f' ')