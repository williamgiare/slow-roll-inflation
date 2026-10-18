#========
# WG Edit
# Canonical physical mode oracle with an analytic scale factor. No primordial
# slow-roll series is used in this oracle. Full linear scalar/tensor modes are
# integrated at finite parameters, with exact physical growing-pump normalization.
#========
import mpmath as mp
from .linear_modes import add,derivative


def gauss(mp,stages):
    x,w=mp.gauss_quadrature(stages,'legendre');c=[(a+1)/2 for a in x];w=[a/2 for a in w]
    A=mp.matrix(stages)
    for j in range(stages):
        poly=[mp.mpf(1)];den=mp.mpf(1)
        for k in range(stages):
            if k==j:continue
            out=[mp.mpf(0)]*(len(poly)+1)
            for i,v in enumerate(poly):out[i]-=c[k]*v;out[i+1]+=v
            poly=out;den*=c[j]-c[k]
        for i in range(stages):A[i,j]=sum(v*c[i]**(k+1)/(k+1) for k,v in enumerate(poly))/den
    return c,w,A


def physical_modes(e,b,*,upper=80,lower='1e-13',step='.16',stages=10,digits=50,initial_terms=18,power_terms=70):
    # a(x)=x^-1 exp[-e (x^b-1)/b], epsilon=v(1+v+b)/(1+v)^2,
    # v=e*x^b. Positive epsilon on the real contour reconstructs a canonical
    # potential via dphi/dN=sqrt(2 epsilon), V=H^2(3-epsilon).
    with mp.workdps(digits):
        e,b=mp.mpf(str(e)),mp.mpf(str(b));lower=mp.mpf(lower);hstar=1+e
        c,w,A=gauss(mp,stages)
        def epsilon(v):return v*(1+v+b)/(1+v)**2
        epsstar=epsilon(e)
        output={}
        for field in ('tensor','scalar'):
            if field=='tensor':
                source=[mp.mpf(0)]*(power_terms+1);source[1]=3-b;source[2]=1
            else:
                u=[mp.mpf(0)]*(power_terms+1);u[0]=-1+b/2
                for k in range(1,power_terms+1):u[k]=b/2*(-1)**(k-1)*((1+b)**(-k)-2)-(1 if k==1 else 0)
                source=[b*k*u[k]+sum(u[j]*u[k-j] for j in range(k+1))-u[k]-(2 if k==0 else 0)
                        for k in range(power_terms+1)]
            def pump(l):
                v=e*mp.exp(b*(l+mp.j*mp.pi/2))
                s=-v
                F=-e/b*(mp.exp(b*(l+mp.j*mp.pi/2))-1) if b else -e*(l+mp.j*mp.pi/2)
                if field=='scalar':
                    s+=b/2*(1+v/(1+v+b)-2*v/(1+v))
                    F+=mp.log(epsilon(v)/epsstar)/2
                return s,F,v
            # Raw linear-mode initial data, expanded in inverse t and analytic v.
            series=[]
            for n in range(1,initial_terms+1):
                row=[]
                for k in range(power_terms+1):
                    if n==1:rhs=-source[k]/2
                    else:
                        D=b*k-(n-1)
                        rhs=(D*D-D)*series[-1][k]/2
                        rhs-=sum(source[j]*series[-1][k-j] for j in range(k+1))/2
                        rhs-=sum((-1)**(n-2-j)*(b*k-j)*series[j-1][k] for j in range(1,n-1))
                    denominator=b*k-n
                    if abs(denominator)<mp.mpf('1e-35'):
                        raise ValueError('Initial-data power resonance: choose nonresonant b or use log-power data.')
                    row.append(rhs/denominator)
                series.append(row)
            def evaluate(poly,v):
                r=mp.mpc(0)
                for a in poly[::-1]:r=r*v+a
                return r
            l=mp.log(upper);stop=mp.log(lower);s,F,v=pump(l)
            R=mp.mpc(1);DR=mp.mpc(0)
            for n,row in enumerate(series,1):
                R+=evaluate(row,v)/mp.mpf(upper)**n
                DR+=evaluate([(b*k-n)*a for k,a in enumerate(row)],v)/mp.mpf(upper)**n
            amplitude=R*mp.exp(-F);velocity=(DR-s*R)*mp.exp(-F)
            while l>stop:
                dl=-min(mp.mpf(step),l-stop)
                ls=[l+dl*a for a in c];ts=[mp.exp(a) for a in ls];ss=[pump(a)[0] for a in ls]
                f=[1+2*t+2/(1+t)-2*s for t,s in zip(ts,ss)]
                j=[2*t*t/(1+t)*s for t,s in zip(ts,ss)]
                inv=mp.eye(stages)-dl*A*mp.diag(f)-dl*dl*A*mp.diag(j)*A
                rhs=mp.matrix([velocity+dl*sum(A[i,k]*j[k]*amplitude for k in range(stages)) for i in range(stages)])
                vv=list(mp.lu_solve(inv,rhs))
                aa=[amplitude+dl*sum(A[i,k]*vv[k] for k in range(stages)) for i in range(stages)]
                amplitude+=dl*sum(w[i]*vv[i] for i in range(stages))
                velocity+=dl*sum(w[i]*(f[i]*vv[i]+j[i]*aa[i]) for i in range(stages))
                l+=dl
            # On the BD contour the coefficient phase changes, not its modulus.
            output[field]=mp.nstr(abs(amplitude)**2/hstar**2, digits-5)
        return {'amplitudes':output,'epsilon_star':mp.nstr(epsstar,30),'h_star':mp.nstr(hstar,30),
                'settings':dict(upper=upper,lower=str(lower),step=step,stages=stages,digits=digits,
                                initial_terms=initial_terms,power_terms=power_terms)}


def local_flow(e,b,depth,digits=60):
    # Numerical Taylor algebra of the analytic background; no SRI flow dictionary.
    with mp.workdps(digits):
        e,b=mp.mpf(str(e)),mp.mpf(str(b));count=depth+2
        def h(l):return 1+e*mp.exp(b*l)
        def epsilon(l):
            v=e*mp.exp(b*l)
            return v*(1+v+b)/(1+v)**2
        def divide(a,den):
            q=[]
            for j in range(min(len(a),len(den))):
                q.append((a[j]-sum(den[k]*q[j-k] for k in range(1,j+1)))/den[0])
            return q
        hjet=mp.taylor(h,0,count)
        def DN(poly):
            d=[(j+1)*poly[j+1] for j in range(len(poly)-1)]
            return [-v for v in divide(d,hjet)]
        ej=mp.taylor(epsilon,0,count);eps=[ej[0]]
        current=ej
        for _ in range(depth):
            current=divide(DN(current),current)
            eps.append(current[0])
        return [mp.nstr(v,digits-5) for v in eps]
#========
