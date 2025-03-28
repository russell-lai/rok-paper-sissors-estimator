from symbolic_variables import *
from sage.all import log

## Symbolic

sym_ring_param = {
    "f" : v_f,
    "fhat" : v_fhat,
    "phi" : v_phi, 
    "q" : v_q, 
    "log_q": log(v_q),
    "Rq_size" : 1, ## We count communication in number of R_q elements
    "e" : v_e, 
    "C" : v_C, 
    "gamma_2" : v_gamma_2, 
    "theta_2" : v_theta_2, 
    "gamma_inf" : v_gamma_inf, 
    "theta_inf" : v_theta_inf, 
    "ring_exp_inf" : v_ring_exp_inf
}

## The following are ring_param precomputed using gen_ring_param(f = ...)

## Powers of 2

ring_param_2048_1024 = {
  'f': 2048,
  'fhat': 1024,
  'phi': 1024,
  'q': 18446744069683058687,
  'log_q': 64,
  'Rq_size': 65536,
  'e': 2,
  'C': 2,
  'gamma_2': 1024,
  'theta_2': 1024,
  'gamma_inf': 1024,
  'theta_inf': 1024,
  'ring_exp_inf': 1024
}

ring_param_1024_512 = {
  'f': 1024,
  'fhat': 512,
 'phi': 512,
 'q': 18446744071864058879,
 'log_q': 64,
 'Rq_size': 32768,
 'e': 2,
 'C': 2,
 'gamma_2': 512,
 'theta_2': 512,
 'gamma_inf': 512,
 'theta_inf': 512,
 'ring_exp_inf': 512}

ring_param_64_32 = {
  'f': 64,
  'fhat': 32,
 'phi': 32,
 'q': 18446744073525002303,
 'log_q': 64,
 'Rq_size': 2048,
 'e': 2,
 'C': 2,
 'gamma_2': 32,
 'theta_2': 32,
 'gamma_inf': 32,
 'theta_inf': 32,
 'ring_exp_inf': 32}

## Primes

ring_param_31_30 = {
  'f': 31,
  'fhat': 31,
 'phi': 30,
 'q': 18446744073664069763,
 'log_q': 64,
 'Rq_size': 1920,
 'e': 2,
 'C': 31,
 'gamma_2': 30,
 'theta_2': 30,
 'gamma_inf': 30,
 'theta_inf': 30,
 'ring_exp_inf': 30}

## Power-Smooth

ring_param_9240_1920 = {'f': 9240,
 'fhat': 4620,
 'phi': 1920,
 'q': 18446744073709505399,
 'log_q': 64,
 'Rq_size': 122880,
 'e': 2,
 'C': 840,
 'gamma_2': 1,
 'theta_2': 1920,
 'gamma_inf': 2,
 'theta_inf': 1920,
 'ring_exp_inf': 1920}

ring_param_7560_1728 = {'f': 7560,
 'fhat': 3780,
 'phi': 1728,
 'q': 18446744073709541519,
 'log_q': 64,
 'Rq_size': 110592,
 'e': 2,
 'C': 280,
 'gamma_2': 1,
 'theta_2': 1728,
 'gamma_inf': 2,
 'theta_inf': 1728,
 'ring_exp_inf': 1728}


ring_param_4620_960 = {'f': 4620,
 'fhat': 2310,
 'phi': 960,
 'q': 18446744073709505399,
 'log_q': 64,
 'Rq_size': 61440,
 'e': 2,
 'C': 420,
 'gamma_2': 1,
 'theta_2': 960,
 'gamma_inf': 2,
 'theta_inf': 960,
 'ring_exp_inf': 960}

ring_param_2520_576 = {
  'f': 2520,
  'fhat': 1260,
 'phi': 576,
 'q': 18446744071570719239,
 'log_q': 64,
 'Rq_size': 36864,
 'e': 2,
 'C': 360,
 'gamma_2': 1,
 'theta_2': 576,
 'gamma_inf': 2,
 'theta_inf': 576,
 'ring_exp_inf': 576}



ring_param_8160_2048 = {'f': 8160,
 'fhat': 4080,
 'phi': 2048,
 'q': 18446744073709551359,
 'log_q': 64,
 'Rq_size': 131072,
 'e': 2,
 'C': 480,
 'gamma_2': 1,
 'theta_2': 2048,
 'gamma_inf': 2,
 'theta_inf': 2048,
 'ring_exp_inf': 2048}
 
ring_param_4080_1024 = {'f': 4080,
 'fhat': 2040,
 'phi': 1024,
 'q': 18446744073709551359,
 'log_q': 64,
 'Rq_size': 65536,
 'e': 2,
 'C': 240,
 'gamma_2': 1,
 'theta_2': 1024,
 'gamma_inf': 2,
 'theta_inf': 1024,
 'ring_exp_inf': 1024}

ring_param_2040_512 = {'f': 2040,
 'fhat': 1020,
 'phi': 512,
 'q': 18446744073709551359,
 'log_q': 64,
 'Rq_size': 32768,
 'e': 2,
 'C': 120,
 'gamma_2': 1,
 'theta_2': 512,
 'gamma_inf': 2,
 'theta_inf': 512,
 'ring_exp_inf': 512}

ring_param_2310_480 = {'f': 2310,
 'fhat': 1155,
 'phi': 480,
 'q': 18446744073709540049,
 'log_q': 64,
 'Rq_size': 30720,
 'e': 2,
 'C': 105,
 'gamma_2': 1,
 'theta_2': 480,
 'gamma_inf': 2,
 'theta_inf': 480,
 'ring_exp_inf': 480}

ring_param_72_24 = {
  'f': 72,
  'fhat': 36,
 'phi': 24,
 'q': 18446744073609183239,
 'log_q': 64,
 'Rq_size': 1536,
 'e': 2,
 'C': 8,
 'gamma_2': 1,
 'theta_2': 24,
 'gamma_inf': 2,
 'theta_inf': 24,
 'ring_exp_inf': 24}

ring_param_60_16 = {
  'f': 60,
  'fhat': 30,
 'phi': 16,
 'q': 18446744073586343939,
 'log_q': 64,
 'Rq_size': 1024,
 'e': 2,
 'C': 12,
 'gamma_2': 1,
 'theta_2': 16,
 'gamma_inf': 2,
 'theta_inf': 16,
 'ring_exp_inf': 16}

ring_param_24_8 = {'f': 24,
 'fhat': 12,
 'phi': 8,
 'q': 18446744073709551359,
 'log_q': 64,
 'Rq_size': 512,
 'e': 2,
 'C': 3,
 'gamma_2': 1,
 'theta_2': 8,
 'gamma_inf': 2,
 'theta_inf': 8,
 'ring_exp_inf': 8}

ring_param_12_4 = {'f': 12,
 'fhat': 6,
 'phi': 4,
 'q': 18446744073709551359,
 'log_q': 64,
 'Rq_size': 256,
 'e': 2,
 'C': 3,
 'gamma_2': 1,
 'theta_2': 4,
 'gamma_inf': 2,
 'theta_inf': 4,
 'ring_exp_inf': 4}

ring_param_6_2 = {'f': 6,
 'fhat': 3,
 'phi': 2,
 'q': 18446744073709551557,
 'log_q': 64,
 'Rq_size': 128,
 'e': 2,
 'C': 3,
 'gamma_2': 2,
 'theta_2': 2,
 'gamma_inf': 2,
 'theta_inf': 2,
 'ring_exp_inf': 2}

## Integers

int_param = {
  'f': 2,
  'fhat': 1,
 'phi': 1,
 'q': 18446744073707521283,
 'log_q': 64,
 'Rq_size': 64,
 'e': 1,
 'C': 2,
 'gamma_2': 1,
 'theta_2': 1,
 'gamma_inf': 1,
 'theta_inf': 1,
 'ring_exp_inf': 1}