# To get the positive (or a canonical root) we need to be explicit in Sage.
# For this we can canoncialize radicals. (Otherwise, sqrt(x*y) and sqrt(x) * sqrt(y) are NOT equal, and such simplifications will not be made.)
# See: https://doc.sagemath.org/html/en/reference/calculus/sage/symbolic/expression.html#sage.symbolic.expression.Expression.canonicalize_radical

from sage.all import log, Infinity, ceil, floor

def has_method(obj, method):
    return callable(getattr(obj, method, None))

def my_canonicalize(x):
    if has_method(x, "canonicalize_radical"):
        x = x.canonicalize_radical()
    if has_method(x, "simplify_full"): 
        x = x.simplify_full()
    if has_method(x, "expand"):
        x = x.expand()
    if has_method(x, "simplify_log"):
        x = x.simplify_log('all')
    return x

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
    amount = int(bits / factor)

    if isinstance(suffix, tuple):
        singular, multiple = suffix
        if amount == 1:
            suffix = singular
        else:
            suffix = multiple
    return str(amount) + suffix


def readable_log(x, sym_mode = False, do_negate = False, do_floor = False):
    if x == 0:
        return Infinity if do_negate else 0
    elif sym_mode:
        return -log(x) if do_negate else log(x)
    elif do_floor:
        return floor(-log(x,2)) if do_negate else floor(log(x,2))
    else:
        return ceil(-log(x,2)) if do_negate else ceil(log(x,2))