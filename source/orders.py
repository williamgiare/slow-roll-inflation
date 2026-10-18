#========
# WG Edit
# Resolve any two of n, p and W without silently accepting inconsistent orders.
#========
from dataclasses import dataclass



class UnsupportedOrderError(ValueError):
    pass


def nonnegative_integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer.")
    return value


def checked_sr_order(p):
    nonnegative_integer(p, "p")
    return p


@dataclass(frozen=True)
class ConsistencyOrder:
    n: int
    p: int
    W: int

    @classmethod
    def resolve(cls, *, n=None, p=None, W=None):
        if sum(value is not None for value in (n, p, W)) < 2:
            raise ValueError("Supply any two of n, p and W; W=n+1+p.")
        for name, value in (("n", n), ("p", p), ("W", W)):
            if value is not None:
                nonnegative_integer(value, name)
        if n is None:
            n = W - p - 1
        elif p is None:
            p = W - n - 1
        elif W is None:
            W = n + 1 + p
        if n < 0 or p < 0 or W < 1 or W != n + 1 + p:
            raise ValueError("Incompatible orders: require n>=0, p>=0 and W=n+1+p.")
        checked_sr_order(p)
        return cls(n, p, W)

# WG Edit: validate the two algebra engines before cached dispatch.
def checked_symbolic(symbolic):
    if type(symbolic) is not bool:
        raise ValueError('symbolic must be a boolean.')
    return symbolic


def checked_method(method):
    if method not in ("iterative", "closed"):
        raise ValueError("Method must be 'iterative' or 'closed'.")
    return method
