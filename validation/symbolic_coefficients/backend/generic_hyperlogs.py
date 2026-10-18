#========
# WG Edit
# Exact ordered Mellin masters from a triangular differential system in rates.
# Coefficients are finite rational combinations of Goncharov polylogarithms.
# This is a research implementation, not a public coefficient backend.
#========
from collections import Counter, defaultdict
from functools import lru_cache
from fractions import Fraction as Q
from itertools import combinations_with_replacement, permutations, product
from math import comb, factorial
import heapq

import mpmath as mp
import sympy as sp
from sympy.polys.fields import field

from mixed_hyperlogs import gamma_jet, rational_jet, jet_mul, t

RF, rt = field('t', sp.QQ)


@lru_cache(None)
def zero_reduce(parameters, rates):
    """Raise negative integer exponents to zero, preserving meromorphic poles.

    A vanishing final rate is integrated out in its convergence half-plane.
    Completely scaleless masters are zero in Mellin analytic regularization.
    """
    rank = len(rates)
    if not any(rates):
        return {}
    if rates[-1] == 0:
        a, v = parameters[-1]
        merged = parameters[:-2]+((parameters[-2][0]+a, parameters[-2][1]+v),)
        return {key:-value/(a+v*rt) for key,value in zero_reduce(merged,rates[:-1]).items()}
    if rank == 1 or all(a == 0 for a,v in parameters):
        return {(parameters,rates):RF.one}
    i = next(j for j,(a,v) in enumerate(parameters) if a < 0)
    a,v = parameters[i]
    denominator = a+v*rt
    out = defaultdict(lambda:RF.zero)
    if rates[i]:
        shifted = parameters[:i]+((a+1,v),)+parameters[i+1:]
        for key,value in zero_reduce(shifted,rates).items():
            out[key] += rates[i]*value/denominator
    for j,sign in ((i,1),(i-1,-1)):
        if not 0 <= j < rank-1:
            continue
        left,right = parameters[j:j+2]
        p = parameters[:j]+((left[0]+right[0],left[1]+right[1]),)+parameters[j+2:]
        r = rates[:j]+(rates[j]+rates[j+1],)+rates[j+2:]
        for key,value in zero_reduce(p,r).items():
            out[key] += sign*value/denominator
    return {key:value for key,value in out.items() if value}


@lru_cache(None)
def core_words(slopes, rates, order):
    """[t^order] of the normalized zero-exponent ordered master.

    The state j has the last m-j+1 parameters and rates merged. Its
    differential equation is in epsilon form, with constant residue matrices.
    Iterating this polynomial expansion only appends standard G letters.
    """
    rank = len(slopes)
    if order == 0:
        return {():1}
    # All states are evaluated on the same path: nonfinal original rates
    # are multiplied by s, while the original final rate remains fixed.
    # The helper below retains that original final rate for merged states.
    return _path_words(slopes,rates,rates[-1],order)[-1]


@lru_cache(None)
def _path_words(slopes,rates,base,order):
    rank = len(slopes)
    if order == 0:
        return tuple({():1} for _ in range(rank))
    previous = _path_words(slopes,rates,base,order-1)
    out = []
    total = sum(slopes)
    for j in range(1,rank+1):
        prefix = sum(slopes[:j-1])
        suffix = total-prefix
        rate = sum(rates[j-1:])-base
        value = defaultdict(int)
        # d f_j = t[-prefix ds/s - suffix ds/(s-a_j)]f_j
        #         +t prefix[ds/s-ds/(s-a_j)]f_(j-1).
        for word,c in previous[j-1].items():
            if prefix: value[(Q(0),)+word] -= prefix*c
            if rate: value[(-Q(base,rate),)+word] -= suffix*c
        if prefix:
            for word,c in previous[j-2].items():
                value[(Q(0),)+word] += prefix*c
                if rate: value[(-Q(base,rate),)+word] -= prefix*c
        out.append({word:c for word,c in value.items() if c})
    return tuple(out)


class GPL(sp.Function):
    """Standard Goncharov G(word;1), with tangential endpoint regularization."""
    nargs = 1

    @classmethod
    def eval(cls,word):
        if not word: return sp.S.One
        if len(set(word)) == 1:
            a = word[0]
            if a in (0,1): return sp.S.Zero
            return sp.log(1-1/a)**len(word)/sp.factorial(len(word))
        if tuple(word[:-1]) == (0,)*(len(word)-1):
            return -sp.polylog(len(word),1/word[-1])

    def _eval_evalf(self,prec):
        digits = int(prec*0.30103)+10
        value = gpl_value(tuple(Q(int(a.p),int(a.q)) for a in self.args[0]),digits)
        return sp.Float(mp.nstr(value,digits),precision=prec)


def word_expression(words):
    return sp.Add(*(sp.Integer(value)*GPL(sp.Tuple(*(sp.Rational(a.numerator,a.denominator)
                                                     for a in word)))
                    for word,value in words.items()))


@lru_cache(None)
def core_expression(slopes,rates,order):
    return word_expression(core_words(slopes,rates,order))


@lru_cache(None)
def gpl_series(word, z, digits):
    """Convergent power-log series; z lies strictly inside all nonzero letters."""
    with mp.workdps(digits+15):
        x = mp.mpf(z.numerator)/z.denominator
        radius = min((abs(a) for a in word if a),default=Q(1))
        ratio = x/(mp.mpf(radius.numerator)/radius.denominator)
        if ratio >= 1: raise ValueError('Power series requested outside convergence disk.')
        size = max(30,int((digits+15)*mp.log(10)/(-mp.log(ratio)))+10)
        size = 64*((size+63)//64)
        logx = mp.log(x)
        return +sum(x**n*sum(c*logx**j for j,c in enumerate(p))
                    for n,p in enumerate(_power_log_series(word,size,digits)))


@lru_cache(2048)
def _power_log_series(word,size,digits):
    """Share common word suffixes between distinct polylogarithm constants."""
    with mp.workdps(digits+15):
        def primitive(poly,n):
            if n == 0: return [mp.mpf(0)]+[v/(j+1) for j,v in enumerate(poly)]
            out = [mp.mpf(0)]*len(poly)
            for j,v in enumerate(poly):
                for d in range(j+1):
                    out[j-d] += v*(-1)**d*(factorial(j)//factorial(j-d))/mp.mpf(n)**(d+1)
            return out
        if not word: return tuple((mp.mpf(n==0),) for n in range(size))
        tail = _power_log_series(word[1:],size,digits)
        if not word[0]: return tuple(tuple(primitive(p,n)) for n,p in enumerate(tail))
        a = mp.mpf(word[0].numerator)/word[0].denominator
        accum = []
        out = [(mp.mpf(0),)]
        for n in range(1,size):
            p = tail[n-1]
            if len(p)>len(accum): accum.extend([mp.mpf(0)]*(len(p)-len(accum)))
            for j in range(len(accum)):
                accum[j] = -(p[j] if j < len(p) else 0)/a+accum[j]/a
            out.append(tuple(primitive(accum,n)))
        return tuple(out)


def evaluate(expression,digits=30):
    """Evaluate all GPL constants at one shared working precision."""
    from coefficients import C,C_VALUE
    working=digits+15
    values={}
    for atom in sorted(expression.atoms(GPL),key=lambda a:(len(a.args[0]),str(a))):
        word=tuple(Q(int(a.p),int(a.q)) for a in atom.args[0])
        values[atom]=sp.Float(mp.nstr(gpl_value(word,working),working),working)
    return expression.xreplace(values).subs(C,C_VALUE).evalf(digits)


@lru_cache(None)
def gpl_value(word,digits=50):
    if not word: return mp.mpf(1)
    if any(0 < a < 1 for a in word):
        raise ValueError('An interior singular letter requires a contour prescription.')
    radius = min((abs(a) for a in word if a),default=Q(1))
    split = min(Q(1,2),radius/2)
    with mp.workdps(digits+10):
        return +sum((-1)**k*gpl_series(tuple(1-a for a in word[:k][::-1]),1-split,digits+5)
                     *gpl_series(word[k:],split,digits+5) for k in range(len(word)+1))


@lru_cache(None)
def core_jet(parameters,rates,maximum,phase_integer,phase_slope):
    """Restore Gamma factors and the original common outgoing scale (-2i)."""
    if len(rates)==1:
        a,v = parameters[0]
        extra = sp.Integer(rates[0])**(-a)
        # The nonunit real scale contributes exp[-v*t*log(rate)].
        gamma = gamma_jet(((a,v,1),),phase_integer,phase_slope,maximum,extra)
        exponential = {n:(-v*sp.log(rates[0]))**n/sp.factorial(n)
                       for n in range(maximum-min(gamma)+1)}
        return jet_mul(gamma,exponential,maximum)
    slopes = tuple(v for a,v in parameters)
    total = sum(slopes)
    prefix_factor = sp.prod(sum(slopes[:j]) for j in range(1,len(slopes)))
    gamma = gamma_jet(((0,total,1),),phase_integer,phase_slope,maximum,
                     1/(prefix_factor*t**(len(slopes)-1)))
    base = rates[-1]
    needed = maximum-min(gamma)
    cores = {n:sp.expand(sum((-total*sp.log(base))**j/sp.factorial(j)
                            *core_expression(slopes,rates,n-j) for j in range(n+1)))
             for n in range(needed+1)}
    return jet_mul(gamma,cores,maximum)


def ordered_jet(parameters,rates,maximum=2,phase_integer=0,phase_slope=0):
    out = defaultdict(lambda:sp.S.Zero)
    for (p,r),factor in zero_reduce(tuple(parameters),tuple(rates)).items():
        pole = len(r)
        prefactor = rational_jet(factor.as_expr(),maximum+pole)
        core = core_jet(p,r,maximum-min(prefactor),phase_integer,phase_slope)
        for n,value in jet_mul(prefactor,core,maximum).items(): out[n] += value
    return {n:sp.expand(value) for n,value in out.items()}


@lru_cache(None)
def cycle_patterns(rank):
    result = Counter()
    for order in permutations(range(rank)):
        positions = {label:i for i,label in enumerate(order)}
        result[tuple(sum(positions[v]>positions[label]
                         for v in ((label-1)%rank,(label+1)%rank)) for label in order)] += 1
    return result


@lru_cache(None)
def vertex_weights(rank):
    """Combine all cycle orderings before reducing their Mellin masters."""
    weights = defaultdict(Q)
    for degrees,multiplicity in cycle_patterns(rank).items():
        for conjugates in product(*(range(d+1) for d in degrees)):
            weight = Q(multiplicity,2**rank)
            for d,c in zip(degrees,conjugates): weight *= comb(d,c)
            weights[conjugates] += weight
    return dict(weights)


@lru_cache(None)
def trace_reduction(slopes):
    """Symmetric source trace, combining prefactors before expansion."""
    rank = len(slopes)
    reduced = defaultdict(lambda:RF.zero)
    assignments = Counter(permutations(slopes))
    for assigned,count in assignments.items():
        for conjugates,base_weight in vertex_weights(rank).items():
            weight = base_weight*Q(count,factorial(rank))
            rates = tuple(1-c for c in conjugates)
            # y^2, |y|^2 and conjugate(y)^2; coefficient of x^-p is i^p times these integers.
            polys = [(1,2,1),(1,0,-1),(1,-2,1)]
            for powers in product(range(3),repeat=rank):
                factor = weight
                for c,p in zip(conjugates,powers): factor *= polys[c][p]
                if not factor: continue
                integer_sum = -rank-sum(powers)
                factor *= (-1)**rank*2**(-integer_sum)
                params = tuple((-1-p,v) for p,v in zip(powers,assigned))
                for key,value in zero_reduce(params,rates).items():
                    reduced[key] += RF(factor.numerator)/factor.denominator*value
    return {key:value for key,value in reduced.items() if value}


def divide_linear_jet(value,a,v,maximum):
    if not a: return {n-1:c/Q(v) for n,c in value.items() if n-1<=maximum}
    out=defaultdict(Q)
    for n,c in value.items():
        term=c/Q(a)
        for power in range(n,maximum+1):
            out[power]+=term
            term*=-Q(v,a)
    return {n:c for n,c in out.items() if c}


@lru_cache(None)
def zero_reduce_jet(parameters,rates,maximum):
    """The same reduction with exact truncated rational Laurent polynomials."""
    rank=len(rates)
    if not any(rates): return {}
    if rates[-1]==0:
        a,v=parameters[-1]
        p=parameters[:-2]+((parameters[-2][0]+a,parameters[-2][1]+v),)
        return {key:{n:-c for n,c in divide_linear_jet(value,a,v,maximum).items()}
                for key,value in zero_reduce_jet(p,rates[:-1],maximum+(a==0)).items()}
    if rank==1 or all(a==0 for a,v in parameters):
        return {(parameters,rates):{0:Q(1)}}
    i=next(j for j,(a,v) in enumerate(parameters) if a<0)
    a,v=parameters[i]
    out=defaultdict(lambda:defaultdict(Q))
    branches=[]
    if rates[i]:
        branches.append((parameters[:i]+((a+1,v),)+parameters[i+1:],rates,rates[i]))
    for j,sign in ((i,1),(i-1,-1)):
        if 0<=j<rank-1:
            left,right=parameters[j:j+2]
            p=parameters[:j]+((left[0]+right[0],left[1]+right[1]),)+parameters[j+2:]
            r=rates[:j]+(rates[j]+rates[j+1],)+rates[j+2:]
            branches.append((p,r,sign))
    for p,r,multiplier in branches:
        for key,value in zero_reduce_jet(p,r,maximum).items():
            for n,c in divide_linear_jet(value,a,v,maximum).items(): out[key][n]+=multiplier*c
    return {key:{n:c for n,c in value.items() if c} for key,value in out.items()
            if any(value.values())}


@lru_cache(None)
def trace_reduction_jet(slopes,maximum):
    rank=len(slopes)
    reduced=defaultdict(lambda:defaultdict(Q))
    assignments=Counter(permutations(slopes))
    for assigned,count in assignments.items():
        for conjugates,base_weight in vertex_weights(rank).items():
            weight=base_weight*Q(count,factorial(rank))
            rates=tuple(1-c for c in conjugates)
            polys=((1,2,1),(1,0,-1),(1,-2,1))
            for powers in product(range(3),repeat=rank):
                factor=weight*(-1)**rank*2**(rank+sum(powers))
                for c,p in zip(conjugates,powers): factor*=polys[c][p]
                if not factor: continue
                params=tuple((-1-p,v) for p,v in zip(powers,assigned))
                for key,value in zero_reduce_jet(params,rates,maximum).items():
                    for n,c in value.items(): reduced[key][n]+=factor*c
    return {key:{n:sp.Rational(c.numerator,c.denominator) for n,c in value.items() if c}
            for key,value in reduced.items() if any(value.values())}


@lru_cache(None)
def trace_reduction_dynamic(slopes,degree):
    """Combine equal integrals before each integration-by-parts reduction.

    The rank and the negative integer degree strictly decrease along every
    branch. Processing these levels in descending order combines all parent
    contributions before reducing a state. Exact Laurent arithmetic retains
    degree+rank terms, sufficient for every remaining auxiliary pole.
    """
    rank=len(slopes)
    buckets={}
    queue=[]
    finals=defaultdict(lambda:defaultdict(Q))
    def add_state(parameters,rates,value,multiplier=1):
        if not any(rates): return
        r=len(rates)
        defect=-sum(a for a,v in parameters)
        level=(-r,-defect)
        if level not in buckets:
            buckets[level]={}
            heapq.heappush(queue,level)
        state=(parameters,rates)
        target=buckets[level].setdefault(state,defaultdict(Q))
        for n,c in value.items():
            if n<=degree+r: target[n]+=multiplier*c
    assignments=Counter(permutations(slopes))
    polys=((1,2,1),(1,0,-1),(1,-2,1))
    for assigned,count in assignments.items():
        for conjugates,base_weight in vertex_weights(rank).items():
            weight=base_weight*Q(count,factorial(rank))
            rates=tuple(1-c for c in conjugates)
            for powers in product(range(3),repeat=rank):
                factor=weight*(-1)**rank*2**(rank+sum(powers))
                for c,p in zip(conjugates,powers): factor*=polys[c][p]
                if factor:
                    add_state(tuple((-1-p,v) for p,v in zip(powers,assigned)),rates,{0:factor})
    while queue:
        level=heapq.heappop(queue)
        states=buckets.pop(level)
        for (parameters,rates),value in states.items():
            value={n:c for n,c in value.items() if c}
            if not value: continue
            r=len(rates)
            if rates[-1]==0:
                a,v=parameters[-1]
                p=parameters[:-2]+((parameters[-2][0]+a,parameters[-2][1]+v),)
                add_state(p,rates[:-1],divide_linear_jet(value,a,v,degree+r-1),-1)
            elif r==1 or all(a==0 for a,v in parameters):
                target=finals[(parameters,rates)]
                for n,c in value.items(): target[n]+=c
            else:
                i=next(j for j,(a,v) in enumerate(parameters) if a<0)
                a,v=parameters[i]
                divided=divide_linear_jet(value,a,v,degree+r)
                if rates[i]:
                    add_state(parameters[:i]+((a+1,v),)+parameters[i+1:],rates,divided,rates[i])
                for j,sign in ((i,1),(i-1,-1)):
                    if 0<=j<r-1:
                        left,right=parameters[j:j+2]
                        p=parameters[:j]+((left[0]+right[0],left[1]+right[1]),)+parameters[j+2:]
                        b=rates[:j]+(rates[j]+rates[j+1],)+rates[j+2:]
                        add_state(p,b,divided,sign)
    return {key:{n:sp.Rational(c.numerator,c.denominator) for n,c in value.items() if c}
            for key,value in finals.items() if any(value.values())}


def pump_generator(slopes):
    """Symmetric rational pump generator from its finite source subsets."""
    rank = len(slopes)
    values = {}
    for mask in range(1,1 << rank):
        weight = sum(v for j,v in enumerate(slopes) if mask & (1 << j))
        if mask & (mask-1) == 0:
            values[mask] = RF.one/(weight*rt-3)
        else:
            subtotal = RF.zero
            part = (mask-1)&mask
            while part:
                subtotal += values[part]*values[mask^part]
                part = (part-1)&mask
            values[mask] = -subtotal/(weight*rt-3)
    return (2*values[(1 << rank)-1]/(factorial(rank)*sum(slopes)*rt)).as_expr()


@lru_cache(None)
def directional_generator(slopes,maximum=2):
    """Exact chi_m(v1*t,...,vm*t), including its pump normalization."""
    rank = len(slopes)
    if any(type(v) is not int or v <= 0 for v in slopes):
        raise ValueError('Auxiliary directions must contain positive integers.')
    if type(maximum) is not int or maximum < 0:
        raise ValueError('The logarithmic degree must be a nonnegative integer.')
    if rank < 2: raise ValueError('Use the established one-insertion generator for rank one.')
    reduced = trace_reduction_dynamic(slopes,maximum)
    print(f'Rank {rank}: {len(reduced)} zero-exponent master terms.',flush=True)
    out = defaultdict(lambda:sp.S.Zero)
    for (p,r),prefactor in reduced.items():
        if not prefactor: continue
        core = core_jet(p,r,maximum-min(prefactor),0,sum(slopes))
        for n,value in jet_mul(prefactor,core,maximum).items(): out[n] += value
    pump = rational_jet(pump_generator(slopes),maximum)
    result = {}
    for n in range(-rank,maximum+1):
        value = sp.expand(out.get(n,0))
        real = sp.expand((value+value.xreplace({sp.I:-sp.I}))/2)
        result[n] = sp.expand(-sp.Rational(2,rank)*real+pump.get(n,0))
    return result


def diagonal_generator(rank,maximum=2):
    return directional_generator((1,)*rank,maximum)


def patterns(rank,degree,minimum=0):
    if rank == 1:
        if degree >= minimum: yield (degree,)
    else:
        for first in range(minimum,degree//rank+1):
            for tail in patterns(rank-1,degree-first,first): yield (first,)+tail


@lru_cache(None)
def family(rank,degree):
    """Exact mixed moments at arbitrary insertion rank and logarithmic degree.

    Known one-varying families supply only that independent part of the
    symmetric homogeneous polynomial. All mixed coefficients are generated
    analytically from the general Mellin master system, without numerical data.
    """
    from coefficients import moment
    if type(rank) is not int or rank < 1 or type(degree) is not int or degree < 0:
        raise ValueError('A positive insertion rank and a nonnegative logarithmic degree are required.')
    indices = tuple(patterns(rank,degree))
    known = {key:moment(key) for key in indices if rank<=2 or sum(j>0 for j in key)<=1}
    unknown = [key for key in indices if key not in known]
    if not unknown: return known
    def monomial(key,values):
        return sum(sp.prod(sp.Integer(v)**j for v,j in zip(values,exponents))
                   for exponents in set(permutations(key)))/sp.prod(sp.factorial(j) for j in key)
    rows,directions = [],[]
    for rest in combinations_with_replacement(range(1,degree+2),rank-1):
        direction = (1,)+rest
        row = [monomial(key,direction) for key in unknown]
        if sp.Matrix(rows+[row]).rank()>len(rows):
            rows.append(row);directions.append(direction)
        if len(rows)==len(unknown): break
    if len(rows)!=len(unknown): raise RuntimeError('Insufficient independent interpolation directions.')
    rhs=[]
    for direction in directions:
        jet=directional_generator(direction,degree)
        for n in range(-rank,0):
            if jet[n] != 0:
                raise ArithmeticError(f'Uncancelled auxiliary pole at rank {rank}, power {n}.')
        rhs.append(jet[degree]-sum(value*monomial(key,direction) for key,value in known.items()))
    solution = sp.Matrix(rows).inv()*sp.Matrix(rhs)
    return known | {key:sp.expand(value) for key,value in zip(unknown,solution)}
