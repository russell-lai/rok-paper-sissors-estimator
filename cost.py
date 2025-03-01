from sage.all import log, sqrt, ceil, floor
from symbolic_variables import *

class Cost:
    def __init__(self, par, op, **kwargs):
        # On input a set of parameters and an operation in {bdecomp, split, fold, batch, norm, finish} create an object storing the cost of this operation. 
        if op == "bdecomp":
            if "base" in kwargs and "ell" not in kwargs:
                base = kwargs["base"]
                ell = v_ell if par.sym_mode else ceil(log(2*par.beta_wit_inf+1,base)) 
            elif "ell" in kwargs and "base" not in kwargs:
                ell = kwargs["ell"]
                base = ceil( (2*par.beta_wit_inf+1)**(1/ell) ) 
                assert base > 1
            elif "base" not in kwargs and "ell" not in kwargs:
                base = v_base
                ell = v_ell
            else:
                base = kwargs["base"]
                ell = kwargs["ell"]
            
            dim_reduction = 1 # new wit_rdim is wit_rdim/dim_reduction
            cost_ntot = [1,0] # new ntot is cost_ntot[0] * ntot + cost_ntot[1]
            cost_nbot = [1,0] # new nbot is cost_nbot[0] * nbot + cost_nbot[1]
            cost_bun = [ell,0] # new rep is cost_bun[0] * rep + cost_bun[1]
            # cadidate_0 = sqrt(ell * par.rep * par.wit_zdim * par.ring_param["fhat"]) * base / 2
            # cadidate_1 = 
            cost_wit_2 = [0, sqrt(ell * par.rep * par.wit_zdim * par.ring_param["fhat"]) * base / 2] # new beta_wit_2 is cost_wit_2[0] * beta_wit_2 + cost_wit_2[1]
            cost_ext_2 = [2*base**(ell-1),0] if par.sym_mode else [(base**ell-1)/(base-1),0] # Using simplified upper bound in symbolic mode
            cost_wit_inf = [0,base/2] if par.sym_mode else [0, floor(base / 2)]
            cost_ext_inf = [2*base**(ell-1),0] if par.sym_mode else [(base**ell-1)/(base-1),0] # Using simplified upper bound in symbolic mode
            cost_comm = par.ring_param["Rq_size"] * (ell-1) * par.nout * par.rep # new prover_comm is prover_comm + cost_comm
            cost_snd = 0

        if op == "split":
            tdim = kwargs["tdim"] if "tdim" in kwargs else v_tdim
            
            dim_reduction = tdim # new wit_rdim is wit_rdim/dim_reduction
            cost_ntot = [1,0] # new ntot is cost_ntot[0] * ntot + cost_ntot[1]
            cost_nbot = [1,0] # new nbot is cost_nbot[0] * nbot + cost_nbot[1]
            cost_bun = [tdim,0] # new rep is cost_bun[0] * rep + cost_bun[1]
            cost_wit_2 = [1,0] # new beta_wit_2 is cost_wit_2[0] * beta_wit_2 + cost_wit_2[1]
            cost_ext_2 = [sqrt(tdim),0] # ext_expansion_2 = cost_ext_2
            cost_wit_inf = [1,0]
            cost_ext_inf = [1,0]
            cost_comm = par.ring_param["Rq_size"] * ((tdim - 1) * (par.nout -  par.nbot) + (tdim * tdim - 1) * par.nbot) * par.rep # new prover_comm is prover_comm + cost_comm
            cost_snd = 0

        if op == "fold":
            repout = kwargs["repout"] if "repout" in kwargs else v_repout
            
            dim_reduction = 1 # new wit_rdim is wit_rdim/dim_reduction
            cost_ntot = [1,0] # new ntot is cost_ntot[0] * ntot + cost_ntot[1]
            cost_nbot = [1,0] # new nbot is cost_nbot[0] * nbot + cost_nbot[1]
            cost_bun = [0,repout] # new rep is cost_bun[0] * rep + cost_bun[1]
            cost_wit_2 = [sqrt(repout) * par.rep * par.ring_param["gamma_2"], 0] # new beta_wit_2 is cost_wit_2[0] * beta_wit_2 + cost_wit_2[1]
            cost_ext_2 = [2 * sqrt(par.rep) * par.ring_param["theta_2"], 0] # ext_expansion_2 = cost_ext_2
            cost_wit_inf = [par.ring_param["gamma_inf"], 0]
            cost_ext_inf = [2 * par.ring_param["theta_inf"], 0]
            # cost_comm = par.ring_param["Rq_size"] * repout * par.rep # new prover_comm is prover_comm + cost_comm
            cost_comm = 0 # new prover_comm is prover_comm + cost_comm
            cost_snd = par.rep / ( par.ring_param["C"]**repout ) 

        if op == "batch":
            dim_reduction = 1 # new wit_rdim is wit_rdim/dim_reduction
            cost_ntot = [1,0] # new ntot is cost_ntot[0] * ntot + cost_ntot[1]
            cost_nbot = [0,1] # new nbot is cost_nbot[0] * nbot + cost_nbot[1]
            cost_bun = [1,0] # new rep is cost_bun[0] * rep + cost_bun[1]
            cost_wit_2 = [1,0] # new beta_wit_2 is cost_wit_2[0] * beta_wit_2 + cost_wit_2[1]
            cost_ext_2 = [1,0] # ext_expansion_2 = cost_ext_2
            cost_wit_inf = [1,0]
            cost_ext_inf = [1,0]
            cost_comm = 0 # new prover_comm is prover_comm + cost_comm
            cost_snd = par.rep * par.nbot / (par.ring_param["q"]**par.ring_param["e"])

        if op == "norm":
            if par.sym_mode:
                ell = v_ell
                base = v_base
            else:
                base_multiplier = kwargs["base_multiplier"] if "base_multiplier" in kwargs else 2
                base = base_multiplier * par.beta_wit_inf + 1
                ell = ceil(log( par.ring_param["ring_exp_inf"] * par.wit_rdim * par.beta_wit_inf**2, base))
            

            # new_beta_wit_2 = sqrt( ell * par.wit_rdim * par.phi**2 ) * base/2
            new_beta_wit_2 = sqrt( par.beta_wit_2**2 +  ell * par.wit_zdim * par.ring_param["fhat"] * par.beta_wit_inf**2 ) 
            
            dim_reduction = 1 # new wit_rdim is wit_rdim/dim_reduction
            cost_ntot = [1,3] # new ntot is cost_ntot[0] * ntot + cost_ntot[1]
            cost_nbot = [1,3] # new nbot is cost_nbot[0] * nbot + cost_nbot[1]
            cost_bun = [1,ell] # new rep is cost_bun[0] * rep + cost_bun[1]
            # cost_wit_2 = [sqrt(2),0] # new beta_wit_2 is cost_wit_2[0] * beta_wit_2 + cost_wit_2[1]
            cost_wit_2 = [0, new_beta_wit_2] # new beta_wit_2 is cost_wit_2[0] * beta_wit_2 + cost_wit_2[1]
            cost_ext_2 = [0, par.beta_wit_2 / new_beta_wit_2] # ext_expansion_2 = cost_ext_2
            cost_wit_inf = [0, par.beta_wit_inf]
            cost_ext_inf = [0, par.beta_wit_2 / par.beta_wit_inf]
            cost_comm = par.ring_param["Rq_size"] * (ell * (par.sis_param["n"] + par.nbot) # Y'
                                                + 3 * par.rep + 3 * ell) # Y_E, Y'_E # new prover_comm is prover_comm + cost_comm
            cost_snd = 2 * par.wit_rdim / (par.ring_param["q"]**par.ring_param["e"])

        if op == "norm-hp":
            if par.sym_mode:
                ell = v_ell
                base = v_base
            else:
                base_multiplier = kwargs["base_multiplier"] if "base_multiplier" in kwargs else 2
                base = base_multiplier * par.beta_wit_inf + 1
                ell = ceil(log( par.ring_param["ring_exp_inf"] * par.wit_rdim * par.beta_wit_inf**2, base ))
            

            # new_beta_wit_2 = sqrt( ell * par.wit_rdim * par.phi**2 ) * base/2
            # new_beta_wit_2 = sqrt( par.beta_wit_2**2 +  ell * par.wit_zdim * par.ring_param["fhat"] * par.beta_wit_inf**2 ) 
            
            dim_reduction = 1 # new wit_rdim is wit_rdim/dim_reduction
            cost_ntot = [1,1] # new ntot is cost_ntot[0] * ntot + cost_ntot[1]
            cost_nbot = [1,1] # new nbot is cost_nbot[0] * nbot + cost_nbot[1]
            cost_bun = [1,ell] # new rep is cost_bun[0] * rep + cost_bun[1]
            # cost_wit_2 = [sqrt(2),0] # new beta_wit_2 is cost_wit_2[0] * beta_wit_2 + cost_wit_2[1]
            cost_wit_2 = [1, 0] # new beta_wit_2 is cost_wit_2[0] * beta_wit_2 + cost_wit_2[1]
            cost_ext_2 = [0, 1] # ext_expansion_2 = cost_ext_2
            cost_wit_inf = [1, 0]
            cost_ext_inf = [0, 0, 1]
            cost_comm = par.ring_param["Rq_size"] * 2*par.rep*ceil(log(par.wit_rdim,2)) # Y_E, Y'_E # new prover_comm is prover_comm + cost_comm
            cost_snd = 2 * par.wit_rdim / (par.ring_param["q"]**par.ring_param["e"]) # TODO


        if op == "norm-hp-tensor":
            if par.sym_mode:
                ell = v_ell
                base = v_base
            else:
                base_multiplier = kwargs["base_multiplier"] if "base_multiplier" in kwargs else 2
                base = base_multiplier * par.beta_wit_inf + 1
                ell = ceil(log( par.ring_param["ring_exp_inf"] * par.wit_rdim * par.beta_wit_inf**2, base ))
            

            # new_beta_wit_2 = sqrt( ell * par.wit_rdim * par.phi**2 ) * base/2
            # new_beta_wit_2 = sqrt( par.beta_wit_2**2 +  ell * par.wit_zdim * par.ring_param["fhat"] * par.beta_wit_inf**2 ) 
            
            dim_reduction = 1 # new wit_rdim is wit_rdim/dim_reduction
            cost_ntot = [1,1] # new ntot is cost_ntot[0] * ntot + cost_ntot[1]
            cost_nbot = [1,1] # new nbot is cost_nbot[0] * nbot + cost_nbot[1]
            cost_bun = [1,ell] # new rep is cost_bun[0] * rep + cost_bun[1]
            # cost_wit_2 = [sqrt(2),0] # new beta_wit_2 is cost_wit_2[0] * beta_wit_2 + cost_wit_2[1]
            cost_wit_2 = [1, 0] # new beta_wit_2 is cost_wit_2[0] * beta_wit_2 + cost_wit_2[1]
            cost_ext_2 = [0, 1] # ext_expansion_2 = cost_ext_2
            cost_wit_inf = [1, 0]
            cost_ext_inf = [0, 0, 1]
            cost_comm = par.ring_param["log_q"] * 2*par.rep*ceil(log(par.wit_rdim * par.ring_param["phi"],2)) # Y_E, Y'_E # new prover_comm is prover_comm + cost_comm
            cost_snd = 2 * par.wit_rdim / (par.ring_param["q"]**par.ring_param["e"]) # TODO

        if op == "finish":
            dim_reduction = 1 # new wit_rdim is wit_rdim/dim_reduction
            cost_ntot = [1,0] # new ntot is cost_ntot[0] * ntot + cost_ntot[1]
            cost_nbot = [0,0] # new nbot is cost_nbot[0] * nbot + cost_nbot[1]
            cost_bun = [0,0] # new rep is cost_bun[0] * rep + cost_bun[1]
            cost_wit_2 = [1,0] # new beta_wit_2 is cost_wit_2[0] * beta_wit_2 + cost_wit_2[1]
            cost_ext_2 = [0,1] # ext_expansion_2 = cost_ext_2
            cost_wit_inf = [1,0]
            cost_ext_inf = [0,1]
            cost_comm = par.wit_size # new prover_comm is prover_comm + cost_comm
            cost_snd = 0

        self.op = op
        self.wit_rdim = dim_reduction
        self.ntot = cost_ntot
        self.nbot = cost_nbot
        self.bun = cost_bun
        self.wit_2 = cost_wit_2
        self.ext_2 = cost_ext_2
        self.wit_inf = cost_wit_inf
        self.ext_inf = cost_ext_inf
        self.comm = cost_comm
        self.snderr = cost_snd

    def apply(self, par):
        # Apply cost to the input set of parameters
        cost = self
        par.op = cost.op
        par.wit_rdim = par.wit_rdim/cost.wit_rdim
        par.ntot = cost.ntot[0] * par.ntot + cost.ntot[1]
        par.nbot = cost.nbot[0] * par.nbot + cost.nbot[1]
        par.nout = par.sis_param["n"] + par.nbot
        par.rep = cost.bun[0] * par.rep + cost.bun[1]
        par.beta_wit_2 = cost.wit_2[0] * par.beta_wit_2 + cost.wit_2[1] if par.sym_mode else (cost.wit_2[0] * par.beta_wit_2 + cost.wit_2[1]).n()
        par.ext_expansion_2 = cost.ext_2
        par.beta_wit_inf = cost.wit_inf[0] * par.beta_wit_inf + cost.wit_inf[1] if par.sym_mode else (cost.wit_inf[0] * par.beta_wit_inf + cost.wit_inf[1]).n()
        par.ext_expansion_inf = cost.ext_inf
        par.prover_comm += cost.comm
        par.snderr += cost.snderr
        par.update_readable()