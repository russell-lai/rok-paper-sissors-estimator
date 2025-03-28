from sage.all import copy, show, table
from ring_params import *
from gen_params import sym_sis_param
from cost import Cost
from step import Step
from util import *

class Protocol:
    def __init__(
        self, *, 
        sym_mode = False, 
        ring_param = sym_ring_param,
        sis_param = sym_sis_param,
        rep = v_rep, 
        wit_rdim = v_wit_rdim, 
        nbot = 0, 
        beta_2 = v_beta_2, 
        beta_inf = v_beta_inf
    ):   
        self.sym_mode = sym_mode 
        self.ring_param = ring_param
        self.sis_param = sis_param
        
        self.history = [Step(
            sym_mode = sym_mode,
            ring_param = self.ring_param, 
            sis_param = self.sis_param, 
            op = "init",
            rep = rep, 
            wit_rdim = wit_rdim, 
            nbot = nbot, 
            beta_2 = beta_2, 
            beta_inf = beta_inf
        )]
        self.prover_comm = self.history[0].prover_comm
        self.wit_size = self.history[0].wit_size

    def execute(self, op, **kwargs):
        last_step = self.history[len(self.history)-1]
        next_step = copy(last_step)
        cost = Cost(next_step, op, **kwargs)
        cost.apply(next_step)
        self.history += [next_step]
        self.prover_comm = self.history[len(self.history)-1].prover_comm
        self.wit_size = self.history[len(self.history)-1].wit_size

    def extract(self):
        history = self.history
        l = len(history)
        assert history[l-1].op == "finish"
        history[l-1].beta_ext_2 = history[l-1].beta_wit_2
        history[l-1].beta_ext_inf = history[l-1].beta_wit_inf
        history[l-1].update_readable()
        for i in range(l-1):
            j = l-i-1
            history[j-1].beta_ext_2 = history[j-1].ext_expansion_2[0] * history[j].beta_ext_2 + history[j-1].ext_expansion_2[1] * history[j-1].beta_wit_2
            history[j-1].beta_ext_inf = history[j-1].ext_expansion_inf[0] * history[j].beta_ext_inf + history[j-1].ext_expansion_inf[1] * history[j-1].beta_wit_inf
            # note that ext_expansion_2 is of the form either [e,0] or [0,1].
            # In the first case, the extracted norm is e times that of the next round witness.
            # In the second case, the extracted norm is reset to beta_wit_2
            
            history[j-1].update_readable()
            if 2 * history[j-1].beta_ext_2 > 2**history[j-1].sis_param["log_beta_sis"]:
                print("Warning: SIS norm bound exceeded.")
    
    def print(self):
        print("Global parameters:")
        l = [
            ["conductor", v_f, self.ring_param["f"]],
            ["ring degree", v_phi, self.ring_param["phi"]],
            ["subtractive set size", v_C, self.ring_param["C"]],
            ["log modulus", log(v_q), self.ring_param["log_q"]],
            ["log sis norm", log(v_beta_sis), self.sis_param["log_beta_sis"]],
            ["Rq size", v_ringelq, pretty_size(self.ring_param["Rq_size"])]
        ]
        show(table(list(map(list, zip(*l)))))
        print("")
        print("Protocol history:")
        history = self.history
        l = [
            ["operation", " "] + [step.op for step in history],
            ["#relations", (v_ntop, v_nbot)] + [(step.sis_param["n"], step.nbot) for step in history],
            ["witness dimension", v_wit_rdim] + [step.wit_rdim for step in history],
            ["bundle size", v_rep] + [step.rep for step in history],
            ["log witness norm can-2", log(v_beta_wit_2)] + [step.log_beta_wit_2 for step in history],
            ["log extration norm can-2", log(v_beta_ext_2)] + [step.log_beta_ext_2 for step in history],
            ["log witness norm coeff-inf", log(2*v_beta_wit_inf)] + [step.log_beta_wit_inf for step in history], # Omitting +1 from 2 * beta + 1 for simplicity
            ["log extration norm coeff-inf", log(2*v_beta_ext_inf)] + [step.log_beta_ext_inf for step in history], # Omitting +1 from 2 * beta + 1 for simplicity
            ["prover comm", v_pi] + [step.readable_prover_comm for step in history],
            ["witness size", v_w] + [step.readable_wit_size for step in history],
            ["soundness", -log(v_snderr)] + [step.nlog_snderr for step in history],
        ]
        show(table(list(map(list, zip(*l)))))