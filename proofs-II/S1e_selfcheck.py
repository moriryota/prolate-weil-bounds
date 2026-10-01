"""Own numerical self-checks. Arb certifies auxiliary constants only.
Galerkin eigenfunctions, eigenvalues and sampled variations are diagnostics.
Run from the project root with .venv/bin/python; no imported project code.
"""
import json, math, time
from pathlib import Path
from fractions import Fraction as Q
import mpmath as mp
import numpy as np
from flint import arb, ctx

ROOT=Path(__file__).resolve().parent
mp.mp.dps=110
ctx.prec=320
eps=Q(16*math.factorial(14),30**10)
lb=(1-5*eps)/(1-eps)
print('EXACT epsilon upper =',eps,'float =',float(eps),flush=True)
print('EXACT lambda lower =',lb,'float =',float(lb),flush=True)
assert eps<Q(1,400) and lb>Q(1,2)
bracket=Q(1)+Q(7,60)+Q(35,4*30**2)+Q(105,8*30**3)+Q(105,16*30**4)
assert bracket<2
print('tail recurrence bracket at R^2=30 =',bracket,'< 2',flush=True)
pi=arb.pi()
assert pi>3 and pi<4
print('Arb pi =',pi, '; certified 3 < pi < 4',flush=True)
print('Arb 16 exp(-30) 30^4 =',16*arb(-30).exp()*30**4,flush=True)

def pswf(c,K,parity):
    ks=list(range(parity,2*K+parity,2))
    def A(k): return mp.mpf(k)/mp.sqrt(4*k*k-1) if k else mp.mpf(0)
    mat=mp.matrix(K)
    for i,k in enumerate(ks):
        mat[i,i]=k*(k+1)+c*c*(A(k+1)**2+A(k)**2)
        if i+1<K: mat[i,i+1]=mat[i+1,i]=c*c*A(k+1)*A(k+2)
    E,V=mp.eigsy(mat)
    out={}
    for idx in range(3 if parity==0 else 2):
        n=2*idx+parity
        co=[V[i,idx] for i in range(K)]
        p0=mp.mpf(0); dp0=mp.mpf(0)
        for k,v in zip(ks,co):
            norm=mp.sqrt(mp.mpf(2*k+1)/2)
            if k%2==0: p0+=v*norm*(-1)**(k//2)*mp.binomial(k,k//2)/(2**k)
            else: dp0+=v*norm*(-1)**((k-1)//2)*(k+1)*mp.binomial(k,(k-1)//2)/(2**k)
        eigF=(mp.sqrt(2)*co[0]/p0 if not parity else 1j*c*mp.sqrt(mp.mpf(2)/3)*co[0]/dp0)
        lam=c/(2*mp.pi)*abs(eigF)**2
        if (p0 if not parity else dp0)<0: co=[-v for v in co]
        full=np.zeros(2*K+parity)
        for k,v in zip(ks,co): full[k]=float(v*mp.sqrt(mp.mpf(2*k+1)/2))
        out[n]={'coef':co,'ks':ks,'npcoef':full,'eig_SL':E[idx],'lambda':lam}
    return out

def variation(mu,states,count):
    lam=math.sqrt(mu); a=math.log(lam)
    I0=math.sqrt(lam)*float(states[0]['coef'][0])*math.sqrt(2)
    I4=math.sqrt(lam)*float(states[4]['coef'][0])*math.sqrt(2)
    H=math.hypot(I4,I0)
    coeff=(I4*states[0]['npcoef']-I0*states[4]['npcoef'])/math.sqrt(lam)
    N=math.ceil(mu)-1
    breaks=sorted(set([-a,a]+[a-math.log(n) for n in range(2,N+1)]))
    total=0.; previous=0.
    jumps=0.; smooth=0.
    for lo,hi in zip(breaks[:-1],breaks[1:]):
        y=np.linspace(lo,hi,count)
        u=np.exp(y); active=math.floor(lam/math.exp((lo+hi)/2))
        val=np.zeros(count)
        for n in range(1,active+1): val+=np.polynomial.legendre.legval(np.clip(n*u/lam,-1,1),coeff)
        val*=np.exp(y/2)
        jumps+=abs(val[0]-previous)
        smooth+=float(np.abs(np.diff(val)).sum())
        previous=val[-1]
    jumps+=abs(previous)
    total=smooth+jumps
    bound=9*mu**3*H
    sharper=(6*math.sqrt(2)*mu**1.5+4*math.pi*math.sqrt(2/3)*mu**2.5)*H
    return {'grid_points_per_piece':count,'Hnorm':H,'smooth_variation':smooth,'jump_variation':jumps,
      'full_zero_extension_TV':total,'bound_9mu3':bound,'ratio_bound_TV':bound/total,
      'sharper_bound':sharper,'ratio_sharper_TV':sharper/total}

result={'epsilon_exact':str(eps),'lambda_lower_exact':str(lb),'cases':[],'certification_scope':'Arb constants only; Galerkin and grids are numerical diagnostics'}
for mu in (5,7,11):
    started=time.time();c=2*mp.pi*mu
    states={}
    K1,K2=(96,112) if mu==11 else (56,68)
    for parity in (0,1): states.update(pswf(c,K1,parity))
    # Separate larger truncation to diagnose numerical truncation sensitivity.
    larger={}
    for parity in (0,1): larger.update(pswf(c,K2,parity))
    x=np.linspace(-1,1,40001)
    rows=[]
    for n in range(5):
        st=larger[n];eig=st['lambda'];co=st['npcoef']
        vmax=float(np.max(np.abs(np.polynomial.legendre.legval(x,co))))
        dmax=float(np.max(np.abs(np.polynomial.legendre.legval(x,np.polynomial.legendre.legder(co)))))
        b0=float(mp.sqrt(c/(mp.pi*eig)))
        b1=float(c**mp.mpf('1.5')/mp.sqrt(3*mp.pi*eig))
        deficit=1-eig
        assert 0<eig<1, ('nonphysical Galerkin diagnostic',mu,n,mp.nstr(deficit,20))
        row={'n':n,'lambda':mp.nstr(eig,72),'one_minus_lambda':mp.nstr(deficit,32),
          'lambda_lower':float(lb),'ratio_lower_lambda':float(mp.mpf(lb.numerator)/lb.denominator/eig),
          'ratio_half_lambda':float(mp.mpf('0.5')/eig),
          'lambda_truncation_difference':mp.nstr(abs(eig-states[n]['lambda']),8),
          'sup_psi_grid':vmax,'sup_psi_bound':b0,'ratio_psi_bound_grid':b0/vmax,
          'sup_derivative_grid':dmax,'derivative_bound_improved':b1,'ratio_derivative_bound_grid':b1/dmax,
          'derivative_bound_requested':math.sqrt(3)*b1}
        rows.append(row)
        print('mu',mu,'n',n,json.dumps(row),flush=True)
    va=variation(mu,larger,5001); vb=variation(mu,larger,10001)
    print('mu',mu,'VARIATION coarse',json.dumps(va),flush=True)
    print('mu',mu,'VARIATION fine',json.dumps(vb),flush=True)
    result['cases'].append({'mu':mu,'truncations':[K1,K2],'eigenfunctions':rows,'variation_coarse':va,'variation_fine':vb,'seconds':time.time()-started})
    (ROOT/'self_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print('FINISHED numerical diagnostics; these are not certified PSWF enclosures',flush=True)
