#========
# WG Edit
# Exact incomplete-beta expansions for mixed Mellin insertions.
# Harmonic polylogarithms use kernels 1/x and 1/(x-1), with upper limit 1/2.
# No numerical coefficient recognition is used.
#========
from collections import defaultdict, Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product, permutations
from math import comb, factorial

import mpmath as mp
import sympy as sp
from sympy.polys.fields import field

t = sp.Symbol('t')


@lru_cache(None)
def shuffle(a, b):
    if not a:
        return {b:1}
    if not b:
        return {a:1}
    out = defaultdict(int)
    for w,c in shuffle(a[1:],b).items():
        out[(a[0],)+w] += c
    for w,c in shuffle(a,b[1:]).items():
        out[(b[0],)+w] += c
    return dict(out)


@lru_cache(None)
def primitive(k, word):
    """Integral from zero of x^k G_word(x), for nonnegative integer k."""
    if not word:
        return {(k+1,()):Q(1,k+1)}
    a, tail = word[0], word[1:]
    out = defaultdict(Q)
    out[k+1,word] += Q(1,k+1)
    for degree in ((k,) if a == 0 else range(k+1)):
        for key,value in primitive(degree,tail).items():
            out[key] -= value/(k+1)
    if a == 1:
        out[0,(1,)+tail] -= Q(1,k+1)
    return {key:value for key,value in out.items() if value}


def log_words(order, a, b):
    # exp(t(a log x+b log(1-x))) coefficient: shuffle factors cancel order!.
    return {word:Q(a)**word.count(0)*Q(b)**word.count(1)
            for word in product((0,1),repeat=order)}


def polynomial(a, b):
    return {a-1+j:Q((-1)**j*comb(b-1,j)) for j in range(b)}


def integrate_terms(terms):
    out = defaultdict(Q)
    for (k,word),value in terms.items():
        if not value:
            continue
        for key,c in primitive(k,word).items():
            out[key] += value*c
    return {key:value for key,value in out.items() if value}


def at_half(terms):
    out = defaultdict(Q)
    for (k,word),value in terms.items():
        out[word] += value*Q(1,2)**k
    return {word:value for word,value in out.items() if value}


@lru_cache(None)
def beta_core(a, b, va, vb, degree):
    terms = {(k,word):value*c for k,value in polynomial(a,b).items()
             for word,c in log_words(degree,va,vb).items() if c}
    return at_half(integrate_terms(terms))


@lru_cache(None)
def double_core(a,b,c,d,va,vb,vc,vd,degree):
    terms = defaultdict(Q)
    for inner_degree in range(degree+1):
        inner = integrate_terms({(k,w):v*c0 for k,v in polynomial(a,d).items()
                                 for w,c0 in log_words(inner_degree,va,vd).items() if c0})
        outer_logs = log_words(degree-inner_degree,vb,vc)
        for (k,word),value in inner.items():
            for power,c0 in polynomial(b,c).items():
                for word2,c1 in outer_logs.items():
                    if not c1:
                        continue
                    for joined,multiplicity in shuffle(word,word2).items():
                        terms[k+power,joined] += value*c0*c1*multiplicity
    return at_half(integrate_terms(terms))


class HPL(sp.Function):
    """Goncharov G word in the alphabet {0,1}, evaluated at 1/2."""
    nargs = 1

    @classmethod
    def eval(cls, word):
        if not word:
            return sp.S.One
        if len(set(word)) == 1:
            return (-sp.log(2))**len(word)/sp.factorial(len(word))
        if tuple(word) == (0,1):
            return -sp.pi**2/12+sp.log(2)**2/2
        if tuple(word) == (0,0,1):
            return -sp.Rational(7,8)*sp.zeta(3)+sp.pi**2*sp.log(2)/12-sp.log(2)**3/6
        if tuple(word) == (0,1,1):
            return sp.zeta(3)/8-sp.log(2)**3/6
        if tuple(word[:-1]) == (0,)*(len(word)-1) and word[-1] == 1:
            return -sp.polylog(len(word),sp.Rational(1,2))

    def _eval_evalf(self, prec):
        digits = int(prec*0.30103)+10
        value = hpl_value(tuple(map(int,self.args[0])),digits)
        return sp.Float(mp.nstr(value,digits),precision=prec)


def expression(words):
    return sp.Add(*(sp.Rational(v.numerator,v.denominator)*HPL(sp.Tuple(*w))
                    for w,v in words.items()))


def lyndon_factors(word):
    factors = []
    i = 0
    while i < len(word):
        j,k = i+1,i
        while j < len(word) and word[k] <= word[j]:
            k = i if word[k] < word[j] else k+1
            j += 1
        size = j-k
        while i <= k:
            factors.append(word[i:i+size])
            i += size
    return factors


@lru_cache(None)
def shuffle_reduction(maximum):
    """Reduce words to a standard Lyndon basis using exact shuffle identities."""
    result = {():sp.S.One}
    for size in range(1,maximum+1):
        for word in product((0,1),repeat=size):
            factors = lyndon_factors(word)
            if len(factors) == 1:
                result[word] = HPL(sp.Tuple(*word))
                continue
            combined = {():1}
            for factor in factors:
                updated = defaultdict(int)
                for w,m in combined.items():
                    for v,n in shuffle(w,factor).items(): updated[v] += m*n
                combined = updated
            value = sp.prod(result[f] for f in factors)
            value -= sum(n*result[w] for w,n in combined.items() if w != word)
            result[word] = sp.expand(value/combined[word])
    return result


def reduce_words(value):
    words = value.atoms(HPL)
    if not words: return value
    replacements = shuffle_reduction(max(len(w.args[0]) for w in words))
    return sp.expand(value.xreplace({w:replacements[tuple(map(int,w.args[0]))] for w in words}))


def compositions(n):
    if n == 0:
        yield ()
    for k in range(1,n+1):
        for tail in compositions(n-k): yield (k,)+tail


def index_word(index):
    return tuple(letter for k in index for letter in (0,)*(k-1)+(1,))


def word_index(word):
    out = []
    zeros = 0
    for letter in word:
        if letter == 0:
            zeros += 1
        else:
            out.append(zeros+1)
            zeros = 0
    assert zeros == 0
    return tuple(out)


@lru_cache(None)
def stuffle(a,b):
    if not a: return {b:1}
    if not b: return {a:1}
    out = defaultdict(int)
    for w,c in stuffle(a[1:],b).items(): out[(a[0],)+w] += c
    for w,c in stuffle(a,b[1:]).items(): out[(b[0],)+w] += c
    for w,c in stuffle(a[1:],b[1:]).items(): out[(a[0]+b[0],)+w] += c
    return dict(out)


@lru_cache(None)
def mzv_reduction(maximum):
    """Derive low-weight MZVs from shuffle, series product and duality identities."""
    result = {():sp.S.One}
    for weight in range(2,maximum+1):
        indices = [a for a in compositions(weight) if a[0]>=2]
        symbols = {a:sp.Symbol('Z_'+'_'.join(map(str,a))) for a in indices}
        equations = [symbols[(weight,)]-sp.zeta(weight)]
        for a in indices:
            dual = word_index(tuple(1-v for v in index_word(a)[::-1]))
            equations.append(symbols[a]-symbols[dual])
        for leftweight in range(1,weight):
            left = [(1,)] if leftweight == 1 else [a for a in compositions(leftweight) if a[0]>=2]
            right = [b for b in compositions(weight-leftweight) if b[0]>=2]
            for a in left:
                for b in right:
                    series = stuffle(a,b)
                    integral = defaultdict(int)
                    for w,c in shuffle(index_word(a),index_word(b)).items(): integral[word_index(w)] += c
                    # The single divergent leading-one symbol cancels in the difference.
                    diff = {key:series.get(key,0)-integral.get(key,0) for key in series.keys()|integral.keys()}
                    assert all(k[0]>=2 or not v for k,v in diff.items())
                    equations.append(sum(v*symbols[k] for k,v in diff.items() if v))
                    if a != (1,):
                        equations.append(sum(v*symbols[k] for k,v in series.items())-result[a]*result[b])
        solution = sp.linsolve(equations,list(symbols.values()))
        assert solution is not sp.EmptySet
        values = next(iter(solution))
        result.update({a:sp.expand(v) for a,v in zip(indices,values)})
    return result


@lru_cache(None)
def half_argument_reduction(maximum):
    """Reflect the two paths 0->1/2->1 and reduce the resulting shuffle identities."""
    mzv = mzv_reduction(maximum)
    replacements = {}
    for weight in range(2,maximum+1):
        equations = []
        for index in compositions(weight):
            if index[0]<2: continue
            word = index_word(index)
            expr = 0
            for k in range(weight+1):
                reflected = tuple(1-v for v in word[:k][::-1])
                expr += (-1)**k*HPL(sp.Tuple(*reflected))*HPL(sp.Tuple(*word[k:]))
            expr -= (-1)**len(index)*mzv[index]
            equations.append(sp.expand(reduce_words(expr).xreplace(replacements)))
        atoms = sorted(set().union(*(e.atoms(HPL) for e in equations)),key=str)
        if not atoms:
            assert all(e == 0 for e in equations),equations
            continue
        solution = sp.linsolve(equations,atoms)
        assert solution is not sp.EmptySet,weight
        values = next(iter(solution))
        replacements.update({a:sp.expand(v) for a,v in zip(atoms,values) if a != v})
    return replacements


def simplify_half(value):
    value = reduce_words(value)
    atoms = value.atoms(HPL)
    if not atoms: return value
    substitutions = half_argument_reduction(max(len(a.args[0]) for a in atoms))
    return sp.expand(value.xreplace(substitutions))


@lru_cache(None)
def hpl_value(word,digits=50):
    with mp.workdps(digits+10):
        size = 4*digits+30
        def integrate_polynomial(p,n):
            if n == 0:
                return [mp.mpf(0)]+[value/(j+1) for j,value in enumerate(p)]
            out = [mp.mpf(0)]*len(p)
            for j,value in enumerate(p):
                for d in range(j+1):
                    out[j-d] += value*(-1)**d*(factorial(j)//factorial(j-d))/mp.mpf(n)**(d+1)
            return out
        @lru_cache(None)
        def series(w):
            if not w:
                return tuple((mp.mpf(1) if n==0 else mp.mpf(0),) for n in range(size))
            tail = series(w[1:])
            if w[0] == 0:
                return tuple(tuple(integrate_polynomial(p,n)) for n,p in enumerate(tail))
            out = [(mp.mpf(0),)]
            accum = []
            for n in range(1,size):
                p = tail[n-1]
                if len(p)>len(accum): accum.extend([mp.mpf(0)]*(len(p)-len(accum)))
                for j,value in enumerate(p): accum[j] -= value
                out.append(tuple(integrate_polynomial(accum,n)))
            return tuple(out)
        L = -mp.log(2)
        return +sum(mp.mpf(2)**(-n)*sum(c*L**j for j,c in enumerate(p))
                    for n,p in enumerate(series(word)))


# Parameter tuples (integer part, slope in t) make every recurrence exact.
def add(a,b): return (a[0]+b[0],a[1]+b[1])
def sub(a,b): return (a[0]-b[0],a[1]-b[1])
def shift(a,n=1): return (a[0]+n,a[1])
def parameter(a): return a[0]+a[1]*t


@lru_cache(None)
def beta_reduce(a,b):
    """B_(1/2)(a,b) as convergent beta cores and elementary boundary terms."""
    if a[0] < 1:
        terms = {key:parameter(add(a,b))/parameter(a)*value
                 for key,value in beta_reduce(shift(a),b).items()}
        key = ('boundary',add(a,b))
        terms[key] = terms.get(key,0)+1/parameter(a)
    elif b[0] < 1:
        terms = {key:parameter(add(a,b))/parameter(b)*value
                 for key,value in beta_reduce(a,shift(b)).items()}
        key = ('boundary',add(a,b))
        terms[key] = terms.get(key,0)-1/parameter(b)
    else:
        return {('beta',a,b):sp.S.One}
    return {key:sp.cancel(value) for key,value in terms.items() if value != 0}


@lru_cache(None)
def double_reduce(a,b,c,d):
    """Integral u^(b-1)(1-u)^(c-1) B_u(a,d), through convergent cores."""
    if a[0] < 1:
        terms = {key:parameter(add(a,d))/parameter(a)*v
                 for key,v in double_reduce(shift(a),b,c,d).items()}
        for key,v in beta_reduce(add(a,b),add(c,d)).items():
            terms[key] = terms.get(key,0)+v/parameter(a)
    elif b[0] < 1 or c[0] < 1:
        move_b = b[0] < 1
        divisor = parameter(b if move_b else c)
        sign = 1 if move_b else -1
        args = (a,shift(b),c,d) if move_b else (a,b,shift(c),d)
        terms = {key:parameter(add(b,c))/divisor*v for key,v in double_reduce(*args).items()}
        # Boundary value 2^(-b-c) B_(1/2)(a,d).
        for key,v in beta_reduce(a,d).items():
            boundary_key = ('times_boundary',add(b,c),key)
            terms[boundary_key] = terms.get(boundary_key,0)+sign*v/divisor
        for key,v in beta_reduce(add(a,b),add(c,d)).items():
            terms[key] = terms.get(key,0)-sign*v/divisor
    else:
        return {('double',a,b,c,d):sp.S.One}
    return {key:sp.cancel(value) for key,value in terms.items() if value != 0}


@lru_cache(None)
def core_coefficient(key,degree):
    if degree < 0: return sp.S.Zero
    kind = key[0]
    if kind == 'boundary':
        a = key[1]
        return sp.Rational(2)**(-a[0])*(-a[1]*sp.log(2))**degree/sp.factorial(degree)
    if kind == 'times_boundary':
        return sum(core_coefficient(('boundary',key[1]),j)*core_coefficient(key[2],degree-j)
                   for j in range(degree+1))
    if kind == 'beta':
        a,b = key[1:]
        return expression(beta_core(a[0],b[0],a[1],b[1],degree))
    if kind == 'double':
        a,b,c,d = key[1:]
        return expression(double_core(a[0],b[0],c[0],d[0],a[1],b[1],c[1],d[1],degree))
    raise ValueError(key)


def integral_jet(terms,maximum,minimum=-3):
    result = {j:sp.S.Zero for j in range(minimum,maximum+1)}
    for key,rational in terms.items():
        series = sp.series(rational,t,0,maximum+1).removeO().expand()
        coefficients = {j:series.coeff(t,j) for j in range(minimum,maximum+1)}
        for n in result:
            result[n] += sum(value*core_coefficient(key,n-j)
                             for j,value in coefficients.items() if value and j <= n)
    return {j:sp.expand(v) for j,v in result.items()}


def rational_jet(value, maximum):
    num,den = map(lambda e:sp.Poly(e,t),sp.fraction(sp.cancel(value)))
    if num.is_zero: return {}
    nv = min(k[0] for k in num.as_dict())
    dv = min(k[0] for k in den.as_dict())
    first = nv-dv
    out = {}
    for n in range(maximum-first+1):
        out[first+n] = (num.nth(nv+n)-sum(den.nth(dv+j)*out[first+n-j]
                                              for j in range(1,n+1)))/den.nth(dv)
    return {k:sp.expand(v) for k,v in out.items() if v}


def jet_mul(a,b,maximum,minimum=None):
    if not a or not b: return {}
    if minimum is None: minimum = min(a)+min(b)
    out = defaultdict(lambda:sp.S.Zero)
    for i,x in a.items():
        for j,y in b.items():
            if minimum <= i+j <= maximum:
                out[i+j] += x*y
    return {k:sp.expand(v) for k,v in out.items() if v != 0}


@lru_cache(None)
def gamma_jet(factors, phase_power, phase_slope, maximum, extra=sp.S.One):
    """Exact Laurent jet of products of Gamma(n+v*t) and (-2i)^(b-w*t)."""
    from .fast import C
    rational = extra*(-2*sp.I)**phase_power
    slopes = defaultdict(lambda:sp.S.Zero)
    for n,v,sign in factors:
        if n >= 1:
            rational *= sp.prod(k+v*t for k in range(1,n))**sign
        else:
            rational *= sp.prod(k+v*t for k in range(n,1))**(-sign)
    r = rational_jet(rational,maximum)
    if not r: return {}
    needed = maximum-min(r)
    for n,v,sign in factors:
        for k in range(1,needed+1):
            slopes[k] += sign*sp.Integer(v)**k
    L = {1:-(C-sp.log(2)+2)*slopes[1]-(sp.log(2)-sp.I*sp.pi/2)*phase_slope}
    L.update({k:(-1)**k*sp.zeta(k)*slopes[k]/k for k in range(2,needed+1)})
    e = {0:sp.S.One}
    for n in range(1,needed+1):
        e[n] = sp.expand(sum(k*L[k]*e[n-k] for k in range(1,n+1))/n)
    return jet_mul(r,e,maximum)


@lru_cache(None)
def directional_three(slopes=(1,1,1), maximum=2):
    """Exact Laurent/Taylor coefficients of chi_3(v1*t,v2*t,v3*t)."""
    from .fast import C
    maxpole = 3
    work = maximum+maxpole
    B = lambda n,v:gamma_jet(((n-1,v,1),),1-n,v,work,
                             (n+v*t+1)/(n+v*t-3))
    def J(na,va,nb,vb):
        result = defaultdict(lambda:sp.S.Zero)
        for p,q in product(range(3),repeat=2):
            a,b = na-1-p,nb-1-q
            value = gamma_jet(((a+b,va+vb,1),(a,va,1),(1-a-b,-va-vb,1),(1-b,-vb,-1)),
                              -a-b,va+vb,work)
            factor = comb(2,p)*comb(2,q)*(-sp.I)**p*sp.I**q
            for degree,c in value.items(): result[degree] += factor*c
        return {k:sp.expand(v) for k,v in result.items()}
    # Combine source assignments and shifted beta representations before expansion.
    slope_sum = sum(slopes)
    assignments = Counter(permutations(slopes))
    rational_field,_ = field('t',sp.QQ_I)
    reduced = defaultdict(lambda:rational_field.zero)
    for (va,vb,vc),multiplicity in assignments.items():
        for p,q,r in product(range(3),repeat=3):
            total = 3+p+q+r
            factor = (multiplicity*comb(2,p)*comb(2,q)*comb(2,r)*(-sp.I)**p*sp.I**(q+r)
                      *(-2*sp.I)**total/sp.prod(slope_sum*t-k for k in range(total+1)))
            terms = double_reduce((-1-p,va),(-1-q,vb),(-1-r,vc),(1+total,-slope_sum))
            factor_exact = rational_field.from_expr(factor)
            for key,c in terms.items():
                reduced[key] += factor_exact*rational_field.from_expr(c)
    mode_T = defaultdict(lambda:sp.S.Zero)
    for key,c in reduced.items():
        prefactor = gamma_jet(((1,slope_sum,1),),0,slope_sum,maximum,c.as_expr())
        for i,value in prefactor.items():
            for d in range(maximum-i+1):
                mode_T[i+d] += value*core_coefficient(key,d)
    print(f'Analytic three-insertion beta reduction: {len(reduced)} regular cores.',flush=True)
    F = defaultdict(lambda:sp.S.Zero,mode_T)
    for (va,vb,vc),multiplicity in assignments.items():
        for p in (0,2):
            value = jet_mul(jet_mul(B(-1-p,va+vb),B(0,vc),work),
                            rational_jet(1/(va*t-1-p),work),maximum)
            for k,c in value.items(): F[k] += multiplicity*c
        for q in (0,2):
            first,second = J(0,va,-1-q,vb+vc),J(-1-q,va+vb,0,vc)
            difference = {k:first.get(k,0)-second.get(k,0) for k in first.keys()|second.keys()}
            value = jet_mul(difference,rational_jet(1/(vb*t-1-q),work),maximum)
            for k,c in value.items(): F[k] += multiplicity*c
        for p,q in product((0,2),repeat=2):
            value = jet_mul(B(-2-p-q,slope_sum),
                            rational_jet(2/((va*t-1-p)*((va+vb)*t-2-p-q)),work),maximum)
            for k,c in value.items(): F[k] += multiplicity*c
    raw = jet_mul(jet_mul(B(0,slopes[0]),B(0,slopes[1]),work),B(0,slopes[2]),maximum)
    pump = rational_jet(-4*sum(1/(3-(slope_sum-v)*t) for v in slopes)/
                         (3*slope_sum*t*(3-slope_sum*t)*sp.prod(3-v*t for v in slopes)),maximum)
    out = {}
    for k in range(-3,maximum+1):
        value = sp.expand(raw.get(k,0)+F.get(k,0))
        imaginary = sp.expand((value-value.xreplace({sp.I:-sp.I}))/(2*sp.I))
        out[k] = sp.expand(imaginary/12+pump.get(k,0))
    return out


def diagonal_three(maximum=2):
    return directional_three((1,1,1),maximum)


@lru_cache(None)
def three_family(total_log_degree):
    """All symmetric three-insertion moments at a chosen logarithmic degree.

    Exact directional Taylor coefficients determine a symmetric homogeneous
    polynomial through a finite rational linear system. No lower-rank mixed
    moment is fitted or used as an unknown input.
    """
    from .fast import moment
    degree = total_log_degree
    patterns = [(a,b,degree-a-b) for a in range(degree//3+1)
                for b in range(a,(degree-a)//2+1)]
    known = {key:moment(key) for key in patterns if sum(j>0 for j in key)<=1}
    unknown = [key for key in patterns if key not in known]
    if not unknown: return known
    def monomial(key,values):
        return sum(sp.prod(sp.Integer(v)**j for v,j in zip(values,exponents))
                   for exponents in set(permutations(key)))/sp.prod(sp.factorial(j) for j in key)
    rows,directions = [],[]
    for b in range(1,degree+2):
        for a in range(1,b+1):
            direction = (1,a,b)
            row = [monomial(key,direction) for key in unknown]
            if sp.Matrix(rows+[row]).rank()>len(rows):
                rows.append(row);directions.append(direction)
            if len(rows)==len(unknown): break
        if len(rows)==len(unknown): break
    assert len(rows)==len(unknown)
    rhs = []
    for direction in directions:
        jet = directional_three(direction,degree)
        for power in range(-3,0):
            assert simplify_half(jet[power]) == 0, (direction,power)
        value = jet[degree]-sum(v*monomial(key,direction) for key,v in known.items())
        rhs.append(reduce_words(value))
    solution = sp.Matrix(rows).inv()*sp.Matrix(rhs)
    return known | {key:simplify_half(value) if degree<=2 else reduce_words(value)
                    for key,value in zip(unknown,solution)}
