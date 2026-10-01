"""Choose exact decimal-rational polynomial trials; eigensolver results are diagnostics."""
from pathlib import Path
import mpmath as mp,numpy as np,json
from scipy.linalg import eigh_tridiagonal
mp.mp.dps=120
out=Path(__file__).resolve().parent
def choose(n,K):
    c=10*mp.pi;par=n%2;deg=[2*j+par for j in range(K)]
    def a(l):return mp.mpf(l)/mp.sqrt(4*l*l-1) if l else mp.mpf(0)
    d=[l*(l+1)+c*c*(a(l)**2+a(l+1)**2) for l in deg]
    e=[c*c*a(l+1)*a(l+2) for l in deg[:-1]]
    eig,V=eigh_tridiagonal(np.array([float(t) for t in d]),np.array([float(t) for t in e]),select='i',select_range=(n//2,n//2))
    shift=mp.mpf(float(eig[0]))+mp.mpf('1e-7')
    v=[mp.mpf(float(t)) for t in V[:,0]]
    piv=[d[0]-shift];lo=[]
    for j in range(1,K):
        lo.append(e[j-1]/piv[-1]);piv.append(d[j]-shift-lo[-1]*e[j-1])
    for _ in range(16):
        rhs=[v[0]]
        for j in range(1,K):rhs.append(v[j]-lo[j-1]*rhs[-1])
        v=[mp.mpf(0)]*K;v[-1]=rhs[-1]/piv[-1]
        for j in range(K-2,-1,-1):v[j]=(rhs[j]-e[j]*v[j+1])/piv[j]
        norm=mp.sqrt(mp.fsum(t*t for t in v));v=[t/norm for t in v]
    coeff=[v[j]*mp.sqrt(mp.mpf(2*l+1)/2) for j,l in enumerate(deg)]
    if par==0:
        center=mp.fsum(coeff[j]*mp.legendre(l,0) for j,l in enumerate(deg));F=mp.sqrt(2)*v[0]/center
    else:
        # P_(2j+1)'(0)=(-1)^j(2j+1) binomial(2j,j)/4^j.
        center=mp.fsum(coeff[j]*(-1)**j*(2*j+1)*mp.binomial(2*j,j)/4**j for j in range(K))
        F=c*mp.sqrt(mp.mpf(2)/3)*v[0]/center
    leak=1-c/(2*mp.pi)*F*F
    return deg,coeff,leak
trials=[];diagnostic=[]
for n in range(5):
    deg,coeff,leak=choose(n,64)
    _,_,large=choose(n,88)
    difference=abs(leak-large)/large
    assert leak>0 and large>0 and difference<mp.mpf('1e-45')
    trials.append({'n':n,'degrees':deg,'coefficients':[mp.nstr(t,100) for t in coeff],'basis':'unnormalised Legendre P_l; coefficient decimal strings are exact rationals'})
    diagnostic.append({'n':n,'K1':64,'K2':88,'leak_K1':str(leak),'leak_K2':str(large),'relative_difference':str(difference),'certified':False})
    print(json.dumps(diagnostic[-1]),flush=True)
(out/'trials.json').write_text(json.dumps({'c':'10*pi','K':64,'digits':100,'trials':trials},indent=2)+'\n')
(out/'diagnostics.json').write_text(json.dumps(diagnostic,indent=2)+'\n')
print('COMPLETE: fixed exact-rational polynomial trials written; diagnostic values are not certificates.',flush=True)
