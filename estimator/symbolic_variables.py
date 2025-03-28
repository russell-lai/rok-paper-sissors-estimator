# Setup variables. Assert domain positive, so that $sqrt(r) r = sqrt(r^3)$ for simplifications
# See: https://doc.sagemath.org/html/en/reference/calculus/sage/symbolic/expression.html#sage.symbolic.expression.Expression.is_negative

from sage.all import var

var("v_rep", latex_name="r", domain='positive')
var("v_repin", latex_name="r_{\\mathrm{in}}", domain='positive')
var("v_repout", latex_name="r_{\\mathrm{out}}", domain='positive')
var("v_wit_rdim", latex_name="m", domain='positive')
var("v_tdim", latex_name="d", domain='positive')
var("v_ntot", latex_name="n", domain='positive')
var("v_ntop", latex_name="n'", domain='positive')
var("v_nbot", latex_name="n''", domain='positive')
var("v_nout", latex_name="n^{\\mathrm{out}}", domain='positive')

var("v_beta", latex_name="\\beta", domain='positive')
var("v_beta_2", latex_name="\\beta_2", domain='positive')
var("v_beta_inf", latex_name="\\beta_\\infty", domain='positive')
var("v_beta_sis", latex_name="\\beta_{\\mathsf{sis}}", domain='positive')
var("v_beta_wit_2", latex_name="\\beta_{\\mathsf{wit},2}", domain='positive')
var("v_beta_ext_2", latex_name="\\beta_{\\mathsf{ext},2}", domain='positive')
var("v_beta_wit_inf", latex_name="\\beta_{\\mathsf{wit},\\infty}", domain='positive')
var("v_beta_ext_inf", latex_name="\\beta_{\\mathsf{ext},\\infty}", domain='positive')


var("v_snderr", latex_name="\\kappa", domain='positive')
var("v_numtr", latex_name="\\#\\mathrm{tr}", domain='positive')

var("v_ell", latex_name="\\ell", domain='positive')
var("v_base", latex_name="b", domain='positive')
var("v_ellip", latex_name="\\ell_{\\mathsf{ip}}", domain='positive')
var("v_baseip", latex_name="b_{\\mathsf{ip}}", domain='positive')

var("v_f", latex_name="\\mathfrak{f}", domain='positive')
var("v_fhat", latex_name="\\hat{\\mathfrak{f}}", domain='positive')
var("v_phi", latex_name="\\varphi", domain='positive')
var("v_q", latex_name="q", domain='positive')
var("v_e", latex_name="e", domain='positive')

var("v_gamma_2", latex_name="\\gamma_2", domain='positive')
var("v_theta_2", latex_name="\\theta_2", domain='positive')
var("v_gamma_inf", latex_name="\\gamma_inf", domain='positive')
var("v_theta_inf", latex_name="\\theta_inf", domain='positive')
var("v_ring_exp_inf", latex_name="\\gamma_{\\mathcal{R}, \\infty}", domain='positive')
var("v_C", latex_name="|\\mathcal{C}|", domain='positive')


var("v_w", latex_name="|w|", domain='positive')
var("v_pi", latex_name="|\\pi|", domain='positive')

var("v_ringelq", latex_name="\\mathcal{R}_{q}", domain='positive') # Ring element (size) placeholder