"""0792: low-memory source generator, one parity per process. Midpoint diagnostics only, not a certificate.

Same matrices as 0784/generate.py (A=E*ME, B=E*M^2E, J=E*(K+U)E, band, pole, U), with these changes:
  * one parity per process (argument), so only d=N/2 columns are ever held;
  * quadrature points are accumulated in chunks of CHUNK rows (no N x N row lists; default 256,
    since Arb products with inner dimension 64 cost about 4x more in total, 0792 README);
  * Bessel values in the frequency panels are evaluated for the requested parity only;
  * the pole block is built per parity from the moment vectors c, s (2cc-2ss, without adding
    and subtracting h0 as 0784 did through WeilSpectral.l2_and_poles);
  * band, pole, U, K_B are saved and freed before the diagnostic;
  * the approximate eigenvector step is optional (--eig).
Parameters for mu=12 are those of 0784; for mu=15 see the case table below (0791/PLAN.md, 0792 README). Frequency and pole quadrature remainders are unpaid. Rounding enclosures only.
"""
from pathlib import Path
from fractions import Fraction
import sys,json,time,hashlib,gc,argparse
import mpmath as mp
from flint import arb,acb,arb_mat,ctx
from matrix_utils import tick,quality,save_matrix,midpoint_ldl,peak_rss_bytes

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'python'))
from weil_spectral import gauss_legendre,legendre_values

ap=argparse.ArgumentParser()
ap.add_argument('case',choices=('mu12','mu15','mu15a'));ap.add_argument('N',type=int);ap.add_argument('bits',type=int)
ap.add_argument('q',type=int);ap.add_argument('parity',type=int,choices=(0,1));ap.add_argument('tag')
ap.add_argument('--chunk',type=int,default=256);ap.add_argument('--threads',type=int,default=1)
ap.add_argument('--eig',action='store_true');ap.add_argument('--out',default=None)
a=ap.parse_args()
N,bits,q,parity,tag=a.N,a.bits,a.q,a.parity,a.tag
OUT=Path(a.out) if a.out else Path(__file__).parent/'runs'
OUT.mkdir(exist_ok=True)
assert N%2==0 and N>=64 and 512<=bits<=8192 and 64<=q<=512
assert not (OUT/f'{tag}_p{parity}_source_manifest.json').exists(), 'Do not overwrite earlier evidence'
ctx.prec=bits;ctx.threads=a.threads;mp.mp.prec=bits

if a.case=='mu12':
    mu=12;Omega=416;gamma=arb(65)/16;expect_u=[207,68]
    prime={2:2,3:3,4:2,5:5,7:7,8:2,9:3,11:11}
    panel_ends=[arb(0),arb(1)/2]+[arb(2)**k for k in range(0,9)]+[arb(Omega)]   # as 0784
else:
    # mu15: s=66/16 with the weight 1+(2/3)(x/L)^2 (0792 scan), gamma=s+1/2, Omega=736.
    # mu15a: the 0778 refined_gate candidate s=69/16, gamma=77/16, Omega=896.
    mu=15;Omega,gamma,expect_u={'mu15':(736,arb(74)/16,[428,142]),'mu15a':(896,arb(77)/16,[557,185])}[a.case]
    prime={2:2,3:3,4:2,5:5,7:7,8:2,9:3,11:11,13:13}
    # 0778 NUMERICS_AND_PLAN: [0,1/2],[1/2,1],...,[32,64], then panels of length at most 64 up to Omega.
    ends=list(range(128,Omega,64))+[Omega]
    panel_ends=[arb(0),arb(1)/2]+[arb(2)**k for k in range(0,7)]+[arb(e) for e in ends]
t0=time.monotonic();L=arb(mu).log()/2
Ctail=(arb(Omega)/(2*arb.pi())).log()-1/arb(Omega)-gamma
assert Ctail>arb(1)/8
c=Omega*L
delta=(-2*(c+1)).exp()/(4*(c+1)*(c+4))
ehi=[(c+1)/(arb.pi()*c*c),(3*c*c+6*c+2)/(arb.pi()*c*c*c)]
u=[int(float((Ctail/e).lower())) for e in ehi]
assert u==expect_u,(u,expect_u)
assert all(u[j]*ehi[j]<Ctail for j in (0,1))
shifts={n:arb(n).log() for n in prime}
weights={n:arb(p).log()/arb(n).sqrt() for n,p in prime.items()}
norm=[(arb(2*j+1)/(2*L)).sqrt() for j in range(N)]
events=sorted({Fraction(1),Fraction(mu),*[max(Fraction(mu,n*n),Fraction(n*n,mu)) for n in prime]})
ab=lambda r:arb(r.numerator)/r.denominator
I=list(range(parity,N,2));d=len(I);CH=a.chunk

# Prime part: A=E*ME and B=E*M^2E, cell by cell, CH quadrature points at a time.
A=arb_mat(d,d);B=arb_mat(d,d)
xs,ws=gauss_legendre(N+2)
for lo_e,hi_e in zip(events,events[1:]):
    mid=(lo_e+hi_e)/2
    active=[(n,1) for n in prime if mid<Fraction(mu,n*n)]+[(n,-1) for n in prime if mid>Fraction(n*n,mu)]
    lo,hi=ab(lo_e).log()/2,ab(hi_e).log()/2;half=(hi-lo)/2;center=(hi+lo)/2
    for k0 in range(0,len(xs),CH):
        phis=[];images=[];wts=[]
        for xi,wi in zip(xs[k0:k0+CH],ws[k0:k0+CH]):
            x=center+half*xi
            P=legendre_values(x/L,N)
            phi=[P[j]*norm[j] for j in I]
            im=[gamma*v for v in phi]
            for n,s in active:
                vals=legendre_values((x+s*shifts[n])/L,N)
                for t,j in enumerate(I): im[t]-=weights[n]*norm[j]*vals[j]
            phis.append(phi);images.append(im);wts.append(2*half*wi)
        Pm=arb_mat(phis);Fm=arb_mat(images)
        WF=arb_mat([[w*v for v in r] for w,r in zip(wts,images)])
        A+=Pm.transpose()*WF;B+=Fm.transpose()*WF
        del phis,images,wts,Pm,Fm,WF
    tick('prime cell complete')
print('prime matrices',time.monotonic()-t0,flush=True)
gc.collect()

# Frequency part: band and K_B on the panels, Bessel values for this parity only.
band=arb_mat(d,d);KB=arb_mat(d,d)
tx,tw=gauss_legendre(q)
for lo,hi in zip(panel_ends,panel_ends[1:]):
    half=(hi-lo)/2;center=(hi+lo)/2
    for k0 in range(0,q,CH):
        rows=[];ww=[];wa=[]
        for xi,wi in zip(tx[k0:k0+CH],tw[k0:k0+CH]):
            t=center+half*xi;z=L*t
            common=(arb.pi()/(2*z)).sqrt()
            rows.append([(2*L*(2*j+1)).sqrt()*((-1)**(j//2))*common*z.bessel_j(arb(j)+arb(1)/2) for j in I])
            arch=acb(arb(1)/4,t/2).digamma().real-arb.pi().log()-gamma
            ww.append(half*wi/arb.pi());wa.append(ww[-1]*arch)
        R=arb_mat(rows)
        band+=R.transpose()*arb_mat([[w*v for v in r] for w,r in zip(ww,rows)])
        KB+=R.transpose()*arb_mat([[w*v for v in r] for w,r in zip(wa,rows)])
        del rows,R
    print('frequency panel',float(lo.mid()),float(hi.mid()),time.monotonic()-t0,flush=True)
    tick('frequency panel complete')

# Pole part: moment vectors as in WeilSpectral.l2_and_poles (N+40 Gauss points), this parity only.
px,pw=gauss_legendre(N+40);cv=[arb(0)]*d;sv=[arb(0)]*d
for xi,w in zip(px,pw):
    x=L*xi;P=legendre_values(xi,N);ch,sh=(x/2).cosh(),(x/2).sinh()
    for t,j in enumerate(I):
        v=P[j]*norm[j]*w*L;cv[t]+=v*ch;sv[t]+=v*sh
pole=arb_mat([[2*cv[i]*cv[j]-2*sv[i]*sv[j] for j in range(d)] for i in range(d)])
tick('pole complete')

h=arb_mat([[((1-delta) if j==0 else 0)-band[j,0]] for j in range(d)])
U=u[parity]*(h*h.transpose())
J=KB+pole+U
manifest={'case':a.case,'mu':mu,'Omega':Omega,'gamma':str(gamma),'u':u,'delta':str(delta),'N':N,'bits':bits,'q':q,'parity':parity,
          'chunk':CH,'threads':a.threads,'panel_ends':[str(p) for p in panel_ends],
          'status':'Rounding enclosures only; quadrature remainder unpaid','files':[],'quality':[],
          'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
for name,M in (('band',band),('pole',pole),('U',U),('A',A),('B',B),('J',J)):
    check=quality(M);manifest['quality'].append({'parity':parity,'matrix':name,**check})
    assert arb(check['max_radius_row_sum'])<arb('1e-70'),(parity,name,check)
    manifest['files'].append({'parity':parity,'matrix':name,**save_matrix(OUT/f'{tag}_p{parity}_{name}.json.gz',M)})
    tick('saved '+name)
del band,KB,pole,U,h;gc.collect()
manifest['assembly_seconds']=time.monotonic()-t0
(OUT/f'{tag}_p{parity}_source_manifest.json').write_text(json.dumps(manifest,indent=2))

# Diagnostic: same quantities as 0784 diagnostic_support.diagnose, freeing intermediates.
start=time.monotonic();Id=arb_mat(d,d)
for i in range(d):Id[i,i]=1
A=(A+A.transpose())/2;B=(B+B.transpose())/2;J=(J+J.transpose())/2;m=arb(1)/2
C=B-m*A;tick('inverse C start');X=C.solve(A-m*Id);del C;gc.collect()
W=(Id-(A-m*Id)*X)/m;W=(W+W.transpose())/2;del X;gc.collect()
Wi=W.inv();Wi=(Wi+Wi.transpose())/2;del W;gc.collect()
T=Wi+J;T=(T+T.transpose())/2
qc=quality(T);assert arb(qc['max_radius_row_sum'])<arb('1e-70')
tref=save_matrix(OUT/f'{tag}_p{parity}_T.json.gz',T)
tick('raw T; start midpoint LDL')
lo,piv=midpoint_ldl(T)
negative=sum(v<0 for v in piv)
rec={'parity':parity,'dimension':d,'source_quality':qc,'T':tref,'negative_midpoint_pivots':negative,
     'positive_midpoint_pivots':sum(v>0 for v in piv),'block_errors_paid':False,'integral_remainders_paid':False,
     'minimum_midpoint_pivot':mp.nstr(min(piv),55)}
(OUT/f'{tag}_p{parity}_pivots.json').write_text(json.dumps([mp.nstr(v,170) for v in piv]))
if negative==0:
    lm=arb_mat([[arb(mp.nstr(v,170)) for v in row] for row in lo]);del lo;gc.collect()
    Z=lm.inv().transpose();del lm;gc.collect()
    for j in range(d):
        scale=arb(mp.nstr(piv[j],170)).sqrt()
        for i in range(d):Z[i,j]/=scale
    R=Z.transpose()*T*Z-Id
    eta=max(sum(abs(R[i,j]).upper() for j in range(d)) for i in range(d));del R;gc.collect()
    z2=sum(abs(Z[i,j]).upper()*abs(Z[i,j]).upper() for i in range(d) for j in range(d))
    assert z2.is_finite() and z2>0
    rec.update(congruence_residual=str(eta),congruence_lt_one=bool(eta<1),frobenius_squared=str(z2),
               surrogate_lower_bound=str((1-eta)/z2),positive_scope='rounded quadrature surrogate only')
    rec['Z']=save_matrix(OUT/f'{tag}_p{parity}_Z.json.gz',Z)
    del Z;gc.collect()
else:
    del lo;gc.collect()
del piv;gc.collect();tick('congruence done')
if a.eig:
    ctx.prec=512
    vals,vecs=T.mid().eig(right=True,algorithm='approx');k=min(range(d),key=lambda j:float(vals[j].real.mid()))
    col=[vecs[i,k] for i in range(d)];j=max(range(d),key=lambda i:float(abs(col[i]).mid()));rot=[v/col[j] for v in col]
    v=arb_mat([[z.real.mid()] for z in rot]);ctx.prec=bits;v=v/(v.transpose()*v)[0,0].sqrt()
    ray=lambda M:(v.transpose()*M*v)[0,0]
    rq=ray(T);res=T*v-rq*v
    rn=sum(abs(res[i,0]).upper()*abs(res[i,0]).upper() for i in range(d)).sqrt()
    rec.update(raw_Winv_plus_J_rayleigh=str(rq),approx_eigenvalue=str(vals[k]),eigenvector_residual_norm=str(rn),
               compression_on_same_vector=str(ray(A+J)),inverse_transfer_loss=str(ray(A-Wi)),D_on_same_vector=str(ray(B-A*A)))
    (OUT/f'{tag}_p{parity}_direction.json').write_text(json.dumps([str(v[i,0]) for i in range(d)]))
rec['diagnostic_seconds']=time.monotonic()-start
rec['total_seconds']=time.monotonic()-t0
rec['peak_rss_bytes']=peak_rss_bytes()
(OUT/f'{tag}_p{parity}_diagnostic.json').write_text(json.dumps(rec,indent=2))
print('complete parity',parity,'negative_pivots',negative,'min_pivot',rec['minimum_midpoint_pivot'][:20],
      'seconds',rec['total_seconds'],'peak_RSS_MB',round(rec['peak_rss_bytes']/1e6,1),flush=True)
