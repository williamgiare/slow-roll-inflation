#========
# WG Edit
# Independent high-order oracle from linear mode and linear growing-pump equations.
# Ordinary source amplitudes are propagated; their logarithm is taken only at
# the endpoint. This module does not use SRI's Riccati/pump coefficient generator
# or weighted-polynomial arithmetic.
#========
from fractions import Fraction as F
from math import factorial, prod
import mpmath as mp


def add(*polys):
    out=[F(0)]*max((len(p) for p in polys),default=0)
    for p in polys:
        for j,c in enumerate(p):out[j]+=c
    while out and out[-1]==0:out.pop()
    return out


def scale(poly,v):return add([c*v for c in poly])
def derivative(poly):return [F(j)*poly[j] for j in range(1,len(poly))]
def multiply(a,b):
    if not a or not b:return []
    out=[F(0)]*(len(a)+len(b)-1)
    for j,c in enumerate(a):
        for k,d in enumerate(b):out[j+k]+=c*d
    return add(out)


def inverse_derivative(poly,n):
    terms=[];power=1
    while poly:
        terms.append(scale(poly,-F(1,n**power)))
        poly=derivative(poly);power+=1
    return add(*terms)


def independent_basis(order):
    def counts(j,budget,prefix):
        if j==order:
            if any(prefix):yield prefix
            return
        for count in range(budget//(j+1)+1):
            yield from counts(j+1,budget-count*(j+1),(*prefix,count))
    zero=(0,)*order
    return (zero,*sorted(counts(0,order,()),key=lambda ex:(sum((j+1)*c for j,c in enumerate(ex)),ex)))


class LinearModeOracle:
    def __init__(self,order,initial_terms=18):
        self.order=order;self.basis=independent_basis(order)
        self.index={ex:k for k,ex in enumerate(self.basis)}
        self.weights=[sum((j+1)*c for j,c in enumerate(ex)) for ex in self.basis]
        self.size=len(self.basis)
        self.pairs=[[] for _ in self.basis]
        for k,ex in enumerate(self.basis[1:],1):
            for a,ea in enumerate(self.basis[1:],1):
                eb=tuple(x-y for x,y in zip(ex,ea))
                if min(eb)>=0 and any(eb) and eb in self.index:
                    self.pairs[k].append((a,self.index[eb]))
        self.predecessors=[]
        for ex in self.basis:
            self.predecessors.append([(j,self.index[tuple(v-(i==j) for i,v in enumerate(ex))])
                                      for j,c in enumerate(ex) if c])
        # Linear pump Z''-3 Z'=g Z, Z(0)=1, polynomial growing branch.
        Z=[[F(1)]];inverse=[[F(1)]];s=[[]]
        for k in range(1,self.size):
            rhs=add(*([[F(0)]*j+Z[b] for j,b in self.predecessors[k]]))
            first=inverse_derivative(rhs,3)
            Z.append([F(0),*[c/F(j+1) for j,c in enumerate(first)]])
            inverse.append(scale(add(Z[k],*(multiply(Z[a],inverse[b]) for a,b in self.pairs[k])),-1))
            s.append(add(derivative(Z[k]),*(multiply(derivative(Z[a]),inverse[b]) for a,b in self.pairs[k])))
        self.Z,self.inverse,self.s=Z,inverse,s
        # Independent inverse-t series for the raw linear mode ratio R=y/y0.
        # R''=(1+2t+2/(1+t))*R'+g R; primes here mean d/d ln(t).
        zero=[[] for _ in self.basis];raw=[]
        for n in range(1,initial_terms+1):
            row=[]
            for k in range(self.size):
                if n==1:
                    rhs=add(*([[F(0)]*j+[F(-1,2)] for j,b in self.predecessors[k] if b==0]))
                else:
                    prev=raw[-1][k]
                    first=add(derivative(prev),scale(prev,-(n-1)))
                    rhs=scale(add(derivative(first),scale(first,-(n-1)),scale(first,-1)),F(1,2))
                    rhs=add(rhs,*(scale([F(0)]*j+raw[-1][b],F(-1,2)) for j,b in self.predecessors[k]))
                    rhs=add(rhs,*(scale(add(derivative(raw[j-1][k]),scale(raw[j-1][k],-j)),
                        -(-1)**(n-2-j)) for j in range(1,n-1)))
                row.append(inverse_derivative(rhs,n))
            raw.append(row)
        self.raw=raw
        self.multiplicities=[factorial(sum(ex))//prod(factorial(c) for c in ex) for ex in self.basis]

    def integrate(self,upper=150,lower=1e-11,step=.125,stages=10,digits=50):
        with mp.workdps(digits):
            def value(poly,L):
                out=mp.mpc(0)
                for c in poly[::-1]:out=out*L+mp.mpf(c.numerator)/c.denominator
                return out
            # Independently integrate cardinal polynomials for a Gauss tableau.
            nodes,weights=mp.gauss_quadrature(stages,'legendre')
            nodes=[(v+1)/2 for v in nodes];weights=[v/2 for v in weights]
            A=mp.matrix(stages)
            for j in range(stages):
                poly=[mp.mpf(1)];den=mp.mpf(1)
                for k in range(stages):
                    if j==k:continue
                    den*=nodes[j]-nodes[k]
                    out=[mp.mpf(0)]*(len(poly)+1)
                    for i,c in enumerate(poly):out[i]-=nodes[k]*c;out[i+1]+=c
                    poly=out
                for i in range(stages):A[i,j]=sum(c*nodes[i]**(k+1)/(k+1) for k,c in enumerate(poly))/den
            l=mp.log(upper);stop=mp.log(lower);L=l+mp.j*mp.pi/2
            R=[mp.mpc(k==0) for k in range(self.size)];DR=[mp.mpc(0)]*self.size
            for n,row in enumerate(self.raw,1):
                for k,poly in enumerate(row):
                    R[k]+=value(poly,L)/mp.mpf(upper)**n
                    DR[k]+=value(add(derivative(poly),scale(poly,-n)),L)/mp.mpf(upper)**n
            V=[value(p,L) for p in self.inverse];S=[value(p,L) for p in self.s]
            amp=[R[k]+V[k]+sum(R[a]*V[b] for a,b in self.pairs[k]) if k else mp.mpc(1)
                 for k in range(self.size)]
            velocity=[DR[k]+sum(DR[a]*V[b] for a,b in self.pairs[k]) for k in range(self.size)]
            velocity=[velocity[k]-S[k]-sum(S[a]*amp[b] for a,b in self.pairs[k]) for k in range(self.size)]
            while l>stop:
                h=-min(mp.mpf(str(step)),l-stop)
                ts=[mp.exp(l+h*c) for c in nodes]
                diagonal=[1+2*t+2/(1+t) for t in ts]
                force=[2*t*t/(1+t) for t in ts]
                ss=[[value(p,l+h*c+mp.j*mp.pi/2) for p in self.s] for c in nodes]
                inv=(mp.eye(stages)-h*A*mp.diag(diagonal))**-1
                stage_amp=[[mp.mpc(1)]*stages]+[[mp.mpc(0)]*stages for _ in range(self.size-1)]
                stage_velocity=[[mp.mpc(0)]*stages for _ in self.basis]
                for k in range(1,self.size):
                    forcing=[force[i]*ss[i][k]+sum(-2*ss[i][a]*stage_velocity[b][i]
                        +force[i]*ss[i][a]*stage_amp[b][i] for a,b in self.pairs[k]) for i in range(stages)]
                    rhs=mp.matrix([velocity[k]+h*sum(A[i,j]*forcing[j] for j in range(stages)) for i in range(stages)])
                    v=list(inv*rhs);stage_velocity[k]=v
                    stage_amp[k]=[amp[k]+h*sum(A[i,j]*v[j] for j in range(stages)) for i in range(stages)]
                    amp[k]+=h*sum(weights[i]*v[i] for i in range(stages))
                    velocity[k]+=h*sum(weights[i]*(diagonal[i]*v[i]+forcing[i]) for i in range(stages))
                l+=h
            logarithm=[mp.mpc(0)]
            for k in range(1,self.size):
                logarithm.append(amp[k]-sum(amp[a]*self.weights[b]*logarithm[b]
                    for a,b in self.pairs[k])/self.weights[k])
            return {tuple(j for j,c in enumerate(ex) for _ in range(c)):
                float(2*mp.re(logarithm[k])/self.multiplicities[k])
                for k,ex in enumerate(self.basis[1:],1)}
#========
