"""0784: mu12 source generator copied from 0778; midpoint diagnostics for a Liu-style bounded comparison, not a certificate.

All source quadratures use Arb for rounding. Frequency and pole quadrature
remainders are unpaid. No Liu source code or numerical arrays are used.
"""
from pathlib import Path
from fractions import Fraction
import sys,json,time,resource,hashlib,gc
from diagnostic_support import tick,quality,save_matrix,diagnose
import mpmath as mp
from flint import arb,acb,arb_mat,ctx

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'python'))
from weil_spectral import gauss_legendre,legendre_values,WeilSpectral
OUT=Path(__file__).parent
N,bits,q=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]);tag=sys.argv[4]
assert N in (640,832,1008) and 1280<=bits<=2560 and q in (224,256)
assert not (OUT/(tag+'_source_manifest.json')).exists(), 'Do not overwrite earlier evidence'
ctx.prec=bits;ctx.threads=1;mp.mp.prec=bits
t0=time.monotonic(); L=arb(12).log()/2; Omega=416; gamma=arb(65)/16
Ctail=(arb(Omega)/(2*arb.pi())).log()-1/arb(Omega)-gamma
assert Ctail>arb(1)/8
c=Omega*L
delta=(-2*(c+1)).exp()/(4*(c+1)*(c+4))
ehi=[(c+1)/(arb.pi()*c*c),(3*c*c+6*c+2)/(arb.pi()*c*c*c)]
u=[int(float((Ctail/e).lower())) for e in ehi]
assert u==[207,68]
assert all(u[j]*ehi[j]<Ctail for j in (0,1))
prime={2:2,3:3,4:2,5:5,7:7,8:2,9:3,11:11}
shifts={n:arb(n).log() for n in prime}
weights={n:arb(p).log()/arb(n).sqrt() for n,p in prime.items()}
norm=[(arb(2*j+1)/(2*L)).sqrt() for j in range(N)]
events=sorted({Fraction(1),Fraction(12),*[max(Fraction(12,n*n),Fraction(n*n,12)) for n in prime]})
ab=lambda r:arb(r.numerator)/r.denominator
ids=[list(range(p,N,2)) for p in (0,1)]
AA=[arb_mat(len(i),len(i)) for i in ids];BB=[arb_mat(len(i),len(i)) for i in ids]
xs,ws=gauss_legendre(N+2)
for a,b in zip(events,events[1:]):
    mid=(a+b)/2
    active=[(n,1) for n in prime if mid<Fraction(12,n*n)]+[(n,-1) for n in prime if mid>Fraction(n*n,12)]
    lo,hi=ab(a).log()/2,ab(b).log()/2;half=(hi-lo)/2;center=(hi+lo)/2
    phis=[];images=[];wts=[]
    for xi,wi in zip(xs,ws):
        x=center+half*xi
        phi=[v*norm[j] for j,v in enumerate(legendre_values(x/L,N))]
        im=[gamma*v for v in phi]
        for n,s in active:
            vals=legendre_values((x+s*shifts[n])/L,N)
            for j in range(N): im[j]-=weights[n]*norm[j]*vals[j]
        phis.append(phi);images.append(im);wts.append(2*half*wi)
    for p in (0,1):
        I=ids[p];P=arb_mat([[r[j] for j in I] for r in phis]);F=arb_mat([[r[j] for j in I] for r in images])
        WF=arb_mat([[w*r[j] for j in I] for w,r in zip(wts,images)])
        AA[p]+=P.transpose()*WF;BB[p]+=F.transpose()*WF
    tick('prime cell complete')
print('prime matrices',time.monotonic()-t0,flush=True)
del phis,images,wts,P,F,WF;gc.collect()

band=[arb_mat(len(i),len(i)) for i in ids];KB=[arb_mat(len(i),len(i)) for i in ids]
tx,tw=gauss_legendre(q);panels=[arb(0),arb(1)/2]+[arb(2)**k for k in range(0,9)]+[arb(Omega)]
for lo,hi in zip(panels,panels[1:]):
    half=(hi-lo)/2;center=(hi+lo)/2;rows=[];ww=[];wa=[]
    for xi,wi in zip(tx,tw):
        t=center+half*xi;z=L*t
        common=(arb.pi()/(2*z)).sqrt()
        values=[(2*L*(2*j+1)).sqrt()*((-1)**(j//2))*common*z.bessel_j(arb(j)+arb(1)/2) for j in range(N)]
        arch=acb(arb(1)/4,t/2).digamma().real-arb.pi().log()-gamma
        rows.append(values);ww.append(half*wi/arb.pi());wa.append(ww[-1]*arch)
    for p in (0,1):
        I=ids[p];R=arb_mat([[r[j] for j in I] for r in rows])
        WR=arb_mat([[w*r[j] for j in I] for w,r in zip(ww,rows)])
        AR=arb_mat([[w*r[j] for j in I] for w,r in zip(wa,rows)])
        band[p]+=R.transpose()*WR;KB[p]+=R.transpose()*AR
    print('frequency panel',float(lo),float(hi),time.monotonic()-t0,flush=True)
    tick('frequency panel complete')

engine=WeilSpectral(L,N,prec=bits)
MP,_=engine.l2_and_poles();hconst=arb(1)/4;hconst=hconst.digamma()-arb.pi().log()

manifest={'mu':12,'Omega':Omega,'gamma':str(gamma),'u':u,'delta':str(delta),'N':N,'bits':bits,'q':q,'status':'Rounding enclosures only; quadrature remainder unpaid','files':[], 'quality':[], 'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
blocks=[]

for parity in (0,1):
    ix=ids[parity];d=len(ix)
    h=arb_mat([[((1-delta) if j==0 else 0)-band[parity][j,0]] for j in range(d)])
    pole=arb_mat([[MP[i,j]-(hconst if i==j else 0) for j in ix] for i in ix])
    U=u[parity]*(h*h.transpose())
    J=KB[parity]+pole+U
    block=dict(A=AA[parity],B=BB[parity],J=J,band=band[parity],pole=pole,U=U)
    for name,M in block.items():
        check=quality(M);manifest['quality'].append({'parity':parity,'matrix':name,**check})
        assert arb(check['max_radius_row_sum'])<arb('1e-70'), (parity,name,check)
        manifest['files'].append({'parity':parity,'matrix':name,**save_matrix(OUT/f'{tag}_p{parity}_{name}.json.gz',M)})
        tick('saved '+name)
    blocks.append(J)
manifest['assembly_seconds']=time.monotonic()-t0
(OUT/(tag+'_source_manifest.json')).write_text(json.dumps(manifest,indent=2))
results=[]
for parity in (0,1):
    results.append(diagnose(AA[parity],BB[parity],blocks[parity],parity,OUT,tag))
result={'status':'MU12 RAW MIDPOINT DIAGNOSTIC; block errors NOT paid; frequency/pole quadrature remainder unpaid','mu':12,'Omega':Omega,'gamma':str(gamma),'u':u,'delta':str(delta),'N':N,'bits':bits,'q':q,'results':results,'seconds':time.monotonic()-t0,'maxrss_bytes_macos':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(OUT/(tag+'_result.json')).write_text(json.dumps(result,indent=2))
print('complete',result['seconds'],result['maxrss_bytes_macos'],flush=True)
