from symbolic_variables import *
from ring_params import *
from gen_params import sym_sis_param
from util import *

class Step:
    def __init__(self, *, 
                 sym_mode = False, 
                 ring_param = sym_ring_param,
                 sis_param = sym_sis_param,
                 op = "init", 
                 rep = v_rep, 
                 wit_rdim = v_wit_rdim, 
                 ntot = v_ntot, 
                 nbot = v_nbot, 
                 beta_2 = v_beta_2, 
                 beta_inf = v_beta_inf
                ):    
        self.sym_mode = sym_mode            # if in sym_mode, instantiate parameters with symbolic parameters 
        self.op = op                        # the operation performed to result in the current parameters
        self.ring_param = ring_param        # ring parameters
        self.sis_param = sis_param          # ring parameters
        self.rep = rep                      # bundle size a.k.a. number of columns of witness
        self.wit_rdim = wit_rdim            # witness dimension over ring a.k.a. number of rows of witness
        self.wit_zdim = wit_rdim * self.ring_param["phi"]            # witness dimension over ZZ
        self.ntot = ntot
        self.nbot = nbot                    # number of bottom relations
        self.nout = sis_param["n"] + nbot 
        self.beta_wit_2 = beta_2            # norm of the current witness in canonical-2 norm
        self.beta_ext_2 = 0                 # norm of the extracted witness in canonical-2 norm (calculated after protocol.execute("finish") and protocol.extract())
        self.beta_wit_inf = beta_inf        # norm of the current witness in coefficient-inf norm
        self.beta_ext_inf = 0               # norm of the extracted witness in coefficient-inf norm (calculated after protocol.execute("finish") and protocol.extract())
        self.ext_expansion_2 = [1,0]        # expansion factor of norm of the extracted witness incurred by the current operation in canonical-2 norm
        self.ext_expansion_inf = [1,0]      # expansion factor of norm of the extracted witness incurred by the current operation in coefficient-inf norm
        self.snderr = 0                     # soundness error
        self.numtr = 1                      # number of transcripts needed for extraction if the current witness were to be sent in plain
        self.prover_comm = 0                # Communication sent by prover until now. Protocols used as reductions, so don't include sending the final witness (of size (ntop + nbot)*wit_rdim*rep * ceil(log(beta)) _bits_)
        self.verifier_comm = 0              # verifier communication

        self.update_readable()

    def update_readable(self):
        self.log_beta_wit_2 = readable_log(self.beta_wit_2, sym_mode = self.sym_mode)
        self.log_beta_ext_2 = readable_log(self.beta_ext_2, sym_mode = self.sym_mode)
        self.log_beta_wit_inf = readable_log(2*self.beta_wit_inf, sym_mode = self.sym_mode, do_floor = True) # Omitting +1 from 2 * beta + 1 for simplicity
        self.log_beta_ext_inf = readable_log(2*self.beta_ext_inf, sym_mode = self.sym_mode, do_floor = True) # Omitting +1 from 2 * beta + 1 for simplicity
        self.wit_zdim = self.wit_rdim * self.ring_param["phi"]
        self.wit_size = self.wit_zdim * self.rep * self.log_beta_wit_inf
        self.readable_prover_comm = self.prover_comm if self.sym_mode else pretty_size(self.prover_comm)
        self.readable_wit_size = self.wit_size if self.sym_mode else pretty_size(self.wit_size)
        self.nlog_snderr = readable_log(self.snderr, sym_mode = self.sym_mode, do_negate = True)
        self.canonicalize()
    
    def __repr__(self):
        return self.op
    
    def canonicalize(self):
        def simplify(x):
            my_canonicalize(x).factor() # FIXME: Not sure if ".factor()" is a good idea or not.
        self.wit_rdim = my_canonicalize(self.wit_rdim)
        self.rep = my_canonicalize(self.rep)
        self.ntot = my_canonicalize(self.ntot)
        self.nbot = my_canonicalize(self.nbot)
        self.nout = my_canonicalize(self.nout)
        self.beta_wit_2 = my_canonicalize(self.beta_wit_2)
        self.beta_ext_2 = my_canonicalize(self.beta_ext_2)
        self.beta_wit_inf = my_canonicalize(self.beta_wit_inf)
        self.beta_ext_inf = my_canonicalize(self.beta_ext_inf)
        self.log_beta_wit_2 = my_canonicalize(self.log_beta_wit_2)
        self.log_beta_ext_2 = my_canonicalize(self.log_beta_ext_2)
        self.log_beta_wit_inf = my_canonicalize(self.log_beta_wit_inf)
        self.log_beta_ext_inf = my_canonicalize(self.log_beta_ext_inf)
        self.snderr = my_canonicalize(self.snderr)
        self.numtr = my_canonicalize(self.numtr)
        self.prover_comm = my_canonicalize(self.prover_comm)
        self.verifier_comm = my_canonicalize(self.verifier_comm)
        self.readable_prover_comm = my_canonicalize(self.readable_prover_comm)
        self.readable_wit_size = my_canonicalize(self.readable_wit_size)
        self.nlog_snderr = my_canonicalize(self.nlog_snderr)