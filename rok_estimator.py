from dataclasses import *
from typing import List
from sage.all import euler_phi, Expression, var, mod, ceil, floor, is_prime_power, sqrt, radical
from lattice_lib import *
from lattice_lib.util import *
import importlib
estimator = importlib.import_module("lattice-estimator.estimator")
import warnings

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
    residue_deg: int | None = None              # inertia degree of q, i.e. R_q splits into fields of cardinality q**residue_deg 
    C: InitVar[SubtractiveSet | None] = None    # subtractive set parameters
    ring_exp_inf: float = field(init=False)     # ring expansion factor in coefficient ell_inf-norm
    nsis: int | None = None                     # module rank of SIS instance
    secpar: int | None = 128                    # target bit security

    def __post_init__(self, C):
        if mod(self.f,4) == 2:
            raise Exception("Conductor f cannot be congruent to 2 modulo 4.")
        if self.log_betasis >= self.log_q:
            raise Exception("Norm bound beta_sis must be smaller than modulus q.")
        self.phi = euler_phi(self.f)
        self.ring_exp_inf = euler_phi(self.f)
        if self.C == None:
            self.C = SubtractiveSet.gen_klno24_cyclotomic(self.f)   # Use the subtractive set construction for cyclotomic fields reported in KLNO24
        if self.residue_deg == None:
            self.residue_deg = ceil(80/self.log_q)        # Aim for 80 bits of soundness for Schwartz-Zippel
        if is_even(self.f):
            self.fhat = self.f/2
        else:
            self.fhat = self.f
        if self.nsis == None:    
            for nsis in range(1, 500):
                sis = estimator.SIS.Parameters(self.phi*nsis, 2**self.log_q, 2**self.log_betasis)
                with HiddenPrints():
                    costs = estimator.SIS.estimate(sis)
                sec = min(cost["rop"] for cost in costs.values())
                if sec >= 2**self.secpar:    
                    self.nsis = nsis
                    self.secpar = floor(log(sec,2))
                    break
        else:   
            sis = estimator.SIS.Parameters(self.phi*self.nsis, 2**self.log_q, 2**self.log_betasis)
            with HiddenPrints():
                costs = estimator.SIS.estimate(sis)
            sec = min(cost["rop"] for cost in costs.values())
            if sec < 2**self.secpar:
                warnings.warn("Specified module rank for SIS is too small for the target security level.")
                # raise Exception("Specified module rank for SIS is too small for the target security level.")

    def size_Rq(self):
        return self.phi * self.log_q

@dataclass
class Cost:
    """
    Data class for costs of reductions of knowledge.
    
    Example:
    sage: Cost(log_beta_ext_2_exp=1,log_beta_ext_inf_exp=1,comm=0,snd=0)
    """
    log_beta_ext_2_exp      : float = 1        # norm expansion factor for canonical ell_2-norm after extraction
    log_beta_ext_inf_exp    : float = 1      # norm expansion factor for coefficient ell_inf-norm after extraction
    log_beta_wit_2          : float | None = None  # set canonical ell_2-norm of extracted witness to this value
    log_beta_wit_inf        : float | None = None # set coefficient ell_inf-norm of extracted witness to this value
    comm : int = 0              # communication cost
    snd : int = 0               # soundness cost
    
    def show(self,label=None,brief=False):
        if self.snd == 0:
            log_snd = -oo
        else:
            log_snd = floor(log(self.snd,2))
        
        
        label_str = f'{label:8s}' if label else 'Cost'
        print(f'{label_str}: communication = {_kb(self.comm):6.2f} KB, soundness error = 2^{log_snd}') # TODO: show in KB
        if not brief:
            print(f' ')

@dataclass
class Relation:
    """
    Data class for relations. A Relation object contains methods modelling reductions of knowledge (RoK) each of which returns a reduced relation and the cost of the RoK. 
    
    Example:
    sage: rel = Relation(ring_params=RingParam(f=60,nsis=2),wdim=60)
    sage: rel_bdecomp, cost_bdecomp = rel.pi_bdecomp(ell=2)
    sage: rel_bdecomp.show()
    """
    ring_params: RingParam = field(repr=False)      # ring parameters
    trivial : bool = False                          # True if the relation is the "True" relation
    nout: int = 1                                   # number of compressed relations, including both commitment and non-commitment relations
    ntop: int = 1                                   # module rank of commitment
    nbot: int = 1                                   # number of non-commitment relations  
    wdim: int = 1                                   # dimension of each witness 
    rep: int = 1                                    # bundle size of witnesses 
    log_beta_wit_2: float = 0                       # log of canonical ell_2-norm bound of witness
    log_beta_wit_inf: float = 0                     # log of coefficient ell_inf-norm bound of witness
    log_beta_wit_2_extract: float = 0               # log of canonical ell_2-norm bound of extracted witness, only computed after extraction
    log_beta_wit_inf_extract: float = 0             # log of coefficient ell_inf-norm bound of extracted witness, only computed after extraction
    
    def show(self,label=None,brief=False):
        label_str = f'{label:8s}' if label else 'Relation'
        if self.trivial:
            print(f'{label_str}: True')
        elif brief:
            print(f'{label_str}: wdim = {self.wdim:6d}, rep = {self.rep:3d}, log_beta_wit_2 = {ceil(self.log_beta_wit_2):3d}, log_beta_wit_inf = {ceil(self.log_beta_wit_inf):3d}')
        else:
            print(f'Relation:')
            print(f'    H * F * W = Y')
            print(f'Statement:')
            print(f'    H: nout x (ntop + nbot)')
            print(f'    F: (ntop + nbot) x wdim')
            print(f'    Y: nout x rep')
            print(f'Witness:')
            print(f'    W: wdim x rep')
            print(f'    ||sigma(W)||_2 <= 2^log_beta_wit_2')
            print(f'    ||psi(W)||_inf <= 2^log_beta_wit_inf')
            print(f'Parameters:')
            print(f'    wdim = {self.wdim}, rep = {self.rep}, log_beta_wit_2 = {ceil(self.log_beta_wit_2)}, log_beta_wit_inf = {ceil(self.log_beta_wit_inf)}')
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
        Returns the relation resulting from the pi_noop RoK and its costs. 
        
        The RoK pi_noop does nothing, it reduces any relation to itself. It is there for testing purposes.
        """
        return deepcopy(self), Cost()
    
    def pi_finish(self):
        """
        Returns the relation resulting from the pi_finish RoK and its costs. 
        
        The pi_finish RoK reduces any relation to True.
        """
        rel = deepcopy(self)
        rel.trivial = True
        comm = self.ring_params.size_Rq() * self.wdim * self.rep # Overestimating. The actual communication is likely smaller because the norm of the witness is smaller than q/2. 
        log_beta_wit_2 = self.log_beta_wit_2
        log_beta_wit_inf = self.log_beta_wit_inf
        cost = Cost(log_beta_wit_2=log_beta_wit_2,log_beta_wit_inf=log_beta_wit_inf,comm=comm)
        return rel, cost
      
    def pi_bdecomp(self,base: int | None = None,ell: int | None = None):
        """
        Returns the relation resulting from the pi_bdecomp RoK and its costs. 
        
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
            
        rel = deepcopy(self)
        rel.rep = self.rep * ell
        rel.log_beta_wit_2      = log(sqrt(ell * self.rep * self.wdim * self.ring_params.fhat * self.ring_params.phi) * base / 2,2)
        rel.log_beta_wit_inf    = log(floor(base / 2),2)
        
        log_beta_ext_2_exp      = log((base**ell-1)/(base-1),2)
        log_beta_ext_inf_exp    = log((base**ell-1)/(base-1),2)
        comm = self.ring_params.size_Rq() * (ell-1) * self.nout * self.rep
        # snd = 0
        cost = Cost(log_beta_ext_2_exp=log_beta_ext_2_exp,log_beta_ext_inf_exp=log_beta_ext_inf_exp,comm=comm)
        
        return rel, cost
    
    def pi_split(self,d: int):
        """     
        Returns the relation resulting from the pi_split RoK and its costs.
        
        Parameters: Target number of chunks 'd'.
        """
        if not d.divides(self.wdim):
            raise Exception("Cannot split the witness into d chunks. Parameter d must divide wdim.")
        
        rel = deepcopy(self)
        rel.wdim    = ZZ(self.wdim / d)
        rel.rep     = self.rep * d
        
        comm    = self.ring_params.size_Rq() * ((d - 1) * self.ntop + (d**2 - 1) * (self.nout - self.ntop)) * self.rep
        snd     = (d-1) / 2**(self.ring_params.log_q * self.ring_params.residue_deg)
        cost = Cost(comm=comm,snd=snd)
        
        #TODO: Raise warning if 2 * extracted norm is greater than beta_sis
        
        return rel, cost
    
    def pi_fold(self,repout: int):
        """
        Returns the relation resulting from the pi_fold RoK and its costs. 
        
        Parameters: Output bundle size 'repout'.
        """
        rel = deepcopy(self)
        repin                   = self.rep
        rel.rep                 = repout
        rel.log_beta_wit_2      = log(sqrt(repout) * repin * self.ring_params.C.gamma_2,2) + self.log_beta_wit_2
        rel.log_beta_wit_inf    = log(sqrt(repout) * repin * self.ring_params.C.gamma_inf,2) + self.log_beta_wit_inf
        
        log_beta_ext_2_exp      = 2 * sqrt(repin) * self.ring_params.C.theta_2
        log_beta_ext_inf_exp    = 2 * sqrt(repin) * self.ring_params.C.theta_inf
        #comm = 0
        snd = repin / (self.ring_params.C.cardinality**repout)
        cost = Cost(log_beta_ext_2_exp=log_beta_ext_2_exp,log_beta_ext_inf_exp=log_beta_ext_inf_exp,snd=snd)
        
        return rel, cost
    
    def pi_batch(self):
        """
        Returns the relation resulting from the pi_batch RoK and its costs. 
        """
        rel = deepcopy(self)
        if self.nout > self.ntop: 
            rel.nout = self.ntop + 1 # TODO: Allow batching into more than 1 row to allow smaller field size.
            
        # comm = 0
        snd = self.rep * self.nbot / (2**(self.ring_params.log_q * self.ring_params.residue_deg))
        cost = Cost(snd=snd)
        return rel, cost
    
        # There is an indirect cost for pi_batch: It increases nout by 1, and in the the communication cost of pi_split nout is multiplied by d**2 instead of d.
    
    def pi_norm(self):
        """
        Returns the relation resulting from the pi_norm RoK and its costs. 
        """
        
        base = 2 * 2**self.log_beta_wit_inf + 1
        ell = ceil(log( self.ring_params.ring_exp_inf * self.wdim * 2**(self.log_beta_wit_inf * 2), base ))
        
        rel = deepcopy(self)
        rel.nout            = self.nout + 3
        rel.nbot            = self.nbot + 3
        rel.rep             = self.rep + ell
        rel.log_beta_wit_2  = log(sqrt( 2**(self.log_beta_wit_2*2) +  ell * self.wdim * self.ring_params.fhat * 2**(self.log_beta_wit_inf*2) ),2) 
        
        comm                = self.ring_params.size_Rq() * (ell * (self.ring_params.nsis + self.nbot) + 3 * self.rep + 3 * ell) 
        snd                 = 2 * self.wdim / (2**(self.ring_params.log_q * self.ring_params.residue_deg))
        log_beta_ext_2      = self.log_beta_wit_2
        log_beta_ext_inf    = log(sqrt(self.ring_params.fhat * self.ring_params.phi * self.wdim * self.rep),2) + self.log_beta_wit_2
        cost = Cost(log_beta_wit_2=log_beta_ext_2,log_beta_wit_inf=log_beta_ext_inf,comm=comm,snd=snd)
        return rel, cost
    
    def pi_ip(self): # TODO: Dummy protocol to be implemented
        """
        Returns the relation resulting from the pi_ip RoK and its costs. 
        """
        return deepcopy(self), Cost()

    
    def pi_aut(self): # TODO: Dummy protocol to be implemented
        """
        Returns the relation resulting from the pi_aut RoK and its costs. 
        """
        return deepcopy(self), Cost()
    
def simulate(rel, ops):
    """
    Simulates the execution of a sequence of RoKs on a relation.
    
    Example: 
    sage: ring_params = RingParam(f=60,log_betasis=32,log_q=64)
    sage: 
    sage: rep = 2**5
    sage: wdim = 2**15
    sage: log_beta_wit_inf = 0
    sage: log_beta_wit_2 = ceil(log(sqrt(wdim * ring_params.phi * ring_params.fhat) * 2**log_beta_wit_inf,2))
    sage: 
    sage: rel = Relation(ring_params=ring_params,wdim=wdim,rep=rep,log_beta_wit_inf=log_beta_wit_inf,log_beta_wit_2=log_beta_wit_2)
    sage: 
    sage: ell = 2
    sage: d = 4
    sage: 
    sage: opener = [("norm", {}), ("batch", {}), ("split", {"d":d}), ("fold", {"repout":rep})]
    sage: loop = [("bdecomp", {"ell":ell}), ("norm", {}), ("batch", {}), ("split", {"d":d}), ("fold", {"repout":rep})]
    sage: ops = opener + loop + loop + opener + loop + [("finish", {})]
    sage: 
    sage: trace, costs = simulate(rel, ops)
    """
    trace = [("init", rel)]
    costs = []

    for op, params in ops:
        new_rel, new_cost = trace[-1][1].execute(op, **params)
        trace += [(op, new_rel)]
        costs += [(op, new_cost)]
    
    return trace, costs
    
# class Protocol:
    # Variables: 
    # Soundness error budget
    # List of Relations # Potential upgrade: Maintain a tree of relations
    # List of ell_2-canonical norms of the extracted witnesses 
    # List of ell_inf-coefficient norms of the extracted witnesses
    # List of communication costs
    # List of soundness errors
    # Accumulated communication cost
    # Accumulated soundness error
    #
    # Methods:
    # - init creates a list containing a single starting relation
    # - execute executes a RoK on the last relation in the list, appends the resulting relation to the list, keep track of the costs. If pi_finish is run, simulate extraction of the witness.


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