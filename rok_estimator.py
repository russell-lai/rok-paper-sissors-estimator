from dataclasses import *
from typing import List
from sage.all import euler_phi, Expression, var, mod, ceil, is_prime_power, sqrt, radical
from lattice_lib import *
from lattice_lib.util import *
import importlib
estimator = importlib.import_module("lattice-estimator.estimator")

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

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
    log_betasis: int = 32                       # log of norm bound beta_sis for SIS problem
    log_q: int = 64                             # log of modulus q (assumed prime)
    e: int | None = None                        # inertia degree of q, i.e. R_q splits into fields of cardinality q**e 
    C: InitVar[SubtractiveSet | None] = None    # subtractive set parameters
    ring_exp_inf: float = field(init=False)     # ring expansion factor in coefficient ell_inf-norm
    nsis: int | None = None
    secpar: int | None = 128

    def __post_init__(self, C):
        if mod(self.f,4) == 2:
            raise Exception("Conductor f cannot be congruent to 2 modulo 4.")
        if self.log_betasis >= self.log_q:
            raise Exception("Norm bound beta_sis must be smaller than modulus q.")
        self.phi = euler_phi(self.f)
        self.ring_exp_inf = euler_phi(self.f)
        if self.C == None:
            self.C = SubtractiveSet.gen_klno24_cyclotomic(self.f)   # Use the subtractive set construction for cyclotomic fields reported in KLNO24
        if self.e == None:
            self.e = ceil(80/self.log_q)        # Aim for 80 bits of soundness for Schwartz-Zippel
        if is_even(self.f):
            self.fhat = self.f/2
        else:
            self.fhat = self.f
        for nsis in range(1, 500):
            sis = estimator.SIS.Parameters(self.phi*nsis, 2**self.log_q, 2**self.log_betasis)
            with HiddenPrints():
                costs = estimator.SIS.estimate(sis)
            sec = min(cost["rop"] for cost in costs.values())
            if sec > 2**self.secpar:    
                self.nsis = nsis
                self.secpar = floor(log(sec,2))
                break

    def size_Rq(self):
        return self.phi * self.log_q

@dataclass
class Cost:
    """
    Data class for costs of reductions of knowledge.
    
    Example:
    sage: Cost(log_beta_ext_2_exp=1,log_beta_ext_inf_exp=1,comm=0,snd=0)
    """
    log_beta_ext_2_exp : float = 1        # norm expansion factor for canonical ell_2-norm after extraction
    log_beta_ext_inf_exp : float = 1      # norm expansion factor for coefficient ell_inf-norm after extraction
    log_beta_wit_2 : float | None = None  # set canonical ell_2-norm of extracted witness to this value
    log_beta_wit_inf : float | None = None # set coefficient ell_inf-norm of extracted witness to this value
    comm : int = 0              # communication cost
    snd : int = 0               # soundness cost

@dataclass
class Relation:
    """
    Data class for relations.
    
    Example:
    sage: Relation(ntop=1,nbot=1,nout=1,wdim=1,rep=1,log_beta_wit_2=1,log_beta_wit_inf=1)
    """
    ring_params: RingParam = field(repr=False)      # ring parameters
    nout: int = 1                                   # number of compressed relations, including both commitment and non-commitment relations
    ntop: int = 1                                   # module rank of commitment
    nbot: int = 1                                   # number of non-commitment relations  
    wdim: int = 1                                   # dimension of each witness 
    rep: int = 1                                    # bundle size of witnesses 
    log_beta_wit_2: float = 0                       # log of canonical ell_2-norm bound of witness
    log_beta_wit_inf: float = 0                     # log of coefficient ell_inf-norm bound of witness
    
    def show(self):
        print(f'Relation: wdim = {self.wdim}, rep = {self.rep}, log_beta_wit_2 = {ceil(self.log_beta_wit_2)}, log_beta_wit_inf = {ceil(self.log_beta_wit_inf)}')
            
    def pi_noop(self):
        """
        Returns the relation resulting from the pi_noop RoK and its costs. 
        
        The RoK pi_noop does nothing. It is there for testing purposes.
        """
        return deepcopy(self), Cost()
    
    def pi_finish(self):
        """
        Returns the relation resulting from the pi_finish RoK and its costs. 
        """
        return deepcopy(self), Cost()
      
    def pi_bdecomp(self,base: int | None = None,ell: int | None = None):
        """
        Returns the relation resulting from the pi_bdecomp RoK and its costs. 
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
            
        rel = deepcopy(self)
        rel.rep = self.rep * ell
        rel.log_beta_wit_2      = log(sqrt(ell * self.rep * self.wdim * self.ring_params.fhat * self.ring_params.phi) * base / 2,2)
        rel.log_beta_wit_inf    = log(floor(base / 2),2)
        
        log_beta_ext_2_exp = log((base**ell-1)/(base-1),2)
        log_beta_ext_inf_exp = log((base**ell-1)/(base-1),2)
        comm = self.ring_params.size_Rq() * (ell-1) * self.nout * self.rep
        cost = Cost(log_beta_ext_2_exp=log_beta_ext_2_exp,log_beta_ext_inf_exp=log_beta_ext_inf_exp,comm=comm)
        
        return rel, cost
    
    def pi_split(self,d: int):
        """     
        Returns the relation resulting from the pi_split RoK and its costs.
        """
        if not d.divides(self.wdim):
            raise Exception("Cannot split the witness into d chunks. Parameter d must divide wdim.")
        
        rel = deepcopy(self)
        rel.wdim = self.wdim / d
        rel.rep = self.rep * d
        
        comm = self.ring_params.size_Rq() * ((d - 1) * self.ntop + (d**2 - 1) * (self.nout - self.ntop)) * self.rep
        cost = Cost(comm=comm)
        
        return rel, cost
    
    def pi_fold(self,repout: int):
        """
        Returns the relation resulting from the pi_fold RoK and its costs. 
        """
        rel = deepcopy(self)
        repin = self.rep
        rel.rep = repout
        rel.log_beta_wit_2 = log(sqrt(repout) * repin * self.ring_params.C.gamma_2,2) + self.log_beta_wit_2
        rel.log_beta_wit_inf = log(sqrt(repout) * repin * self.ring_params.C.gamma_inf,2) + self.log_beta_wit_inf
        
        log_beta_ext_2_exp = 2 * sqrt(repin) * self.ring_params.C.theta_2
        log_beta_ext_inf_exp = 2 * sqrt(repin) * self.ring_params.C.theta_inf
        snd = repin / (self.ring_params.C.cardinality**repout)
        cost = Cost(log_beta_ext_2_exp=log_beta_ext_2_exp,log_beta_ext_inf_exp=log_beta_ext_inf_exp,snd=snd)
        
        return rel, cost
    
    def pi_batch(self):
        """
        Returns the relation resulting from the pi_batch RoK and its costs. 
        """
        rel = deepcopy(self)
        if self.nout > self.ntop: 
            rel.nout = self.ntop + 1
        snd = self.rep * self.nbot / (2**(self.ring_params.log_q * self.ring_params.e))
        cost = Cost(snd=snd)
        return rel, cost
    
        # TODO: Michal claims that sometimes it is not worth it running pi_batch. This cost table does not reflect this observation.
    
    def pi_norm(self):
        """
        Returns the relation resulting from the pi_norm RoK and its costs. 
        """
        
        base = 2 * 2**self.log_beta_wit_inf + 1
        ell = ceil(log( self.ring_params.ring_exp_inf * self.wdim * 2**(self.log_beta_wit_inf * 2), base ))
        
        rel = deepcopy(self)
        rel.nout = self.nout + 3
        rel.nbot = self.nbot + 3
        rel.rep = self.rep + ell
        rel.log_beta_wit_2 = log(sqrt( 2**(self.log_beta_wit_2*2) +  ell * self.wdim * self.ring_params.fhat * 2**(self.log_beta_wit_inf*2) ),2) 
        
        # TODO: Add logic for resetting norm of extracted witness 
        comm = self.ring_params.size_Rq() * (ell * (self.ring_params.nsis + self.nbot) + 3 * self.rep + 3 * ell) 
        snd = 2 * self.wdim / (2**(self.ring_params.log_q * self.ring_params.e))
        log_beta_ext_2 = self.log_beta_wit_2
        log_beta_ext_inf = log(sqrt(self.ring_params.fhat * self.ring_params.phi * self.wdim * self.rep),2) + self.log_beta_wit_2
        cost = Cost(comm=comm,snd=snd)
        return rel, cost
    
    def pi_ip(self):
        """
        Returns the relation resulting from the pi_ip RoK and its costs. 
        """
        return deepcopy(self), Cost()
    
    def pi_aut(self):
        """
        Returns the relation resulting from the pi_aut RoK and its costs. 
        """
        return deepcopy(self), Cost()
    

# class Protocol:
#     # A history is a list of states
#     # A state is a tuple consisting of 
#     # - a relation, 
#     # - the name of the protocol used to reach this state, 
#     # - the added communication cost, 
#     # - the accumulated communication cost, 
#     # - the added soundness error, 
#     # - the accumulated soundness error, 
#     # - the ell_2-canonical norm of the extracted witness (only set if pi_finish is run), 
#     # - the ell_inf-coefficient norm of the extracted witness (only set if pi_finish is run). 
    
#     def execute(pi: RoK, rel: Relation):
#         """
#         Execute RoK pi on relation rel.

#         Example:
#         sage: rel = Relation()
#         sage: new_rel = execute(pi_noop, rel)
#         """
        
#         # TODO: Raise exception if any protocol is run after pi_finish.
        
#         rel.nbot = pi.f_nbot(rel.nbot)
#         rel.nout = pi.f_nout(rel.nout)
#         rel.wdim = pi.f_wdim(rel.wdim)
#         rel.rep = pi.f_rep(rel.rep)
#         rel.log_beta_wit_2 = pi.f_log_beta_wit_2(rel.log_beta_wit_2)
#         rel.log_beta_wit_inf = pi.f_log_beta_wit_inf(rel.log_beta_wit_inf)
#         rel.log_beta_ext_2_exp = pi.f_log_beta_ext_2_exp(rel.log_beta_ext_2_exp)
#         rel.log_beta_ext_inf_exp = pi.f_log_beta_ext_inf_exp(rel.log_beta_ext_inf_exp)
        
#         # TODO: Keep track of communication cost. 
#         # TODO: Keep track of soundness cost. 
        
#         return rel      # TODO: return communication and soundness costs