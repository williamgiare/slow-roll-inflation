# ========
# WG Edit
# Published comparison targets only; never imported by the production engine.
# Auclair & Ringeval, arXiv:2205.12608v2, Eqs. (44), (54), at k=k_*.
# ========
import sympy as sp

REFERENCE = 'https://arxiv.org/html/2205.12608v2'

def amplitudes(p):
    if p not in (0, 1, 2, 3):
        raise ValueError('The published comparison covers amplitude orders 0..3 only.')
    e1,e2,e3,e4 = sp.symbols('epsilon_1:5')
    C = sp.EulerGamma + sp.log(2) - 2
    P,Z = sp.pi**2,sp.zeta(3)
    common1 = -2*(C+1)*e1
    common2 = (P/2+2*C**2+2*C-3)*e1**2
    common3 = -(4*C**3+3*(P-8)*C+14*Z-16)*e1**3/3
    tensor = [sp.Integer(1),common1,
      common2 + (P-12*C**2-24*C-24)*e1*e2/12,
      common3 + (24*C**3+13*P+2*(5*P-36)*C+36*C**2-96)*e1**2*e2/12
      -(4*C**3-P-(P-24)*C+12*C**2+8*Z+8)*(e1*e2**2+e1*e2*e3)/12]
    scalar = [sp.Integer(1),common1-C*e2,
      common2+(7*P/12+C**2-C-6)*e1*e2+(P+4*C**2-8)*e2**2/8+(P-12*C**2)*e2*e3/24,
      common3-(4*C**3+3*(P-8)*C+14*Z-16)*e2**3/24
      +(13*P-8*(P-9)*C+36*C**2-84*Z)*e1**2*e2/12
      -(8*C**3-15*P+6*(P-4)*C-12*C**2+100*Z+16)*e1*e2**2/24
      +(P*C-4*C**3-8*Z+16)*(e2*e3**2+e2*e3*e4)/24
      +(12*C**3+(5*P-48)*C)*e2**2*e3/24
      +(8*C**3+P+6*(P-12)*C-12*C**2-8*Z-8)*e1*e2*e3/12]
    return {field:sp.expand(sum(rows[:p+1])) for field,rows in [('scalar',scalar),('tensor',tensor)]}
# ========
