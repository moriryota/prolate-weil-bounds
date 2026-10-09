"""Independent selected source entries: analytic polynomial A/B and q=320 J.

No project/candidate modules imported. Ball radii for quadrature rounding and
analytic quadrature remainders are reported separately.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import gzip, json, hashlib, time, resource, sys
from flint import arb, acb, ctx

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
ctx.prec=int(sys.argv[1]) if len(sys.argv)>1 else 1280
assert ctx.prec in (640,1280)
ctx.threads=1
t0=time.monotonic()
L=arb(12).log()/2; O=arb(416); pi=arb.pi(); gamma=arb(65)/16
pp={2:2,3:3,4:2,5:5,7:7,8:2,9:3,11:11}
pairs={0:[(1,1),(5,3)],1:[(2,0)]}
degrees={p:sorted({p,*[2*k+p for ij in pairs[p] for k in ij]}) for p in (0,1)}
def ab(q):return arb(q.numerator)/q.denominator
def add(a,b):
 z=[a[k] if k<len(a) else 0 for k in range(max(len(a),len(b)))]
 for k,v in enumerate(b):z[k]+=v
 return z
def scale(a,s):return [s*v for v in a]
def mul(a,b):
 z=[arb(0)]*(len(a)+len(b)-1)
 for i,v in enumerate(a):
  for j,w in enumerate(b):z[i+j]+=v*w
 return z
def shifted(a,h):
 return [sum((a[k]*comb(k,j)*h**(k-j) for k in range(j,len(a))),arb(0)) for j in range(len(a))]
def integral(a,left,right):
 return L*sum((v*(right**(k+1)-left**(k+1))/(k+1) for k,v in enumerate(a)),arb(0))
# Legendre recurrence in exact rational coefficients, then normalized balls.
leg=[[F(1)],[F(0),F(1)]]
for k in range(1,10):
 leg.append(scale(add(scale([F(0)]+leg[k],F(2*k+1)),scale(leg[k-1],F(-k))),F(1,k+1)))
phi={j:scale([ab(v) for v in leg[j]],(arb(2*j+1)/(2*L)).sqrt()) for j in set(sum(degrees.values(),[]))}
events=sorted([arb(0),arb(1)]+[abs(1-arb(n).log()/L) for n in pp])
AB={(p,i,j,k):arb(0) for p in pairs for i,j in pairs[p] for k in ('A','B')}
for left,right in zip(events,events[1:]):
 mid=(left+right)/2
 active=[(sign*arb(n).log()/L,arb(base).log()/arb(n).sqrt())
         for n,base in pp.items() for sign in (-1,1)
         if -1<mid+sign*arb(n).log()/L<1]
 images={j:scale(phi[j],gamma) for j in phi}
 for j in phi:
  for h,w in active:images[j]=add(images[j],scale(shifted(phi[j],h),-w))
 for p in pairs:
  for i,j in pairs[p]:
   di,dj=2*i+p,2*j+p
   AB[p,i,j,'A']+=2*integral(mul(phi[di],images[dj]),left,right)
   AB[p,i,j,'B']+=2*integral(mul(images[di],images[dj]),left,right)
print('analytic prime integrals completed',time.monotonic()-t0,flush=True)
# Different quadrature order from source: q=320, exactness degree 639.
q=320
nodes=[arb.legendre_p_root(q,k,weight=True) for k in range(q)]
band={(p,i,j):arb(0) for p in pairs for i,j in pairs[p]}
kband=band.copy()
col={(p,j):arb(0) for p in degrees for j in degrees[p]}
panels=[arb(0),arb(1)/2]+[arb(2)**k for k in range(9)]+[O]
for lo,hi in zip(panels,panels[1:]):
 half=(hi-lo)/2;mid=(hi+lo)/2
 for xi,wi in nodes:
  t=mid+half*xi;z=L*t
  vals={j:(2*L*(2*j+1)).sqrt()*(-1)**(j//2)*(pi/(2*z)).sqrt()*z.bessel_j(arb(j)+arb(1)/2) for j in phi}
  g=acb(arb(1)/4,t/2).digamma().real-pi.log()-gamma
  wt=half*wi/pi
  for p in pairs:
   for i,j in pairs[p]:
    v=wt*vals[2*i+p]*vals[2*j+p]
    band[p,i,j]+=v;kband[p,i,j]+=v*g
   for j in degrees[p]:col[p,j]+=wt*vals[j]*vals[p]
 print('frequency panel',float(lo),float(hi),time.monotonic()-t0,flush=True)
 assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<8_500_000_000
# Exact antiderivative moments of exp(x/2), not pole Gauss quadrature.
mom=[]
for k in range(11):
 moment=2*(L**k*(L/2).exp()-(-L)**k*(-L/2).exp())
 if k:moment-=2*k*mom[-1]
 mom.append(moment)
polecoef={j:sum((coef*mom[k]/L**k for k,coef in enumerate(phi[j])),arb(0)) for j in phi}
d=(-2*(O*L+1)).exp()/(4*(O*L+1)*(O*L+4))
J={}
for p in pairs:
 for i,j in pairs[p]:
  di,dj=2*i+p,2*j+p
  hi=(1-d if i==0 else arb(0))-col[p,di]
  hj=(1-d if j==0 else arb(0))-col[p,dj]
  J[p,i,j]=kband[p,i,j]+2*(-1)**p*polecoef[di]*polecoef[dj]+(207 if p==0 else 68)*hi*hj
# Entry error; intentionally omit beneficial 1/pi. Old-source and own errors.
H=arb(5)/2*arb(150).exp()
eb_new=4*O*H*arb(2)**(-639)
eb_old=4*O*H*arb(2)**(-511)
pins=json.loads((ROOT/'notes/0785_gpt_A4_mu12_certificate/SOURCE_PINS.json').read_text())
records=[]
for p in pairs:
 for name in ('A','B','J'):
  info=pins['files'][f'p{p}_{name}']; path=ROOT/info['path']
  assert hashlib.sha256(path.read_bytes()).hexdigest()==info['sha256']
  with gzip.open(path,'rt') as f:source=json.load(f)
  for i,j in pairs[p]:
   ref=arb(source[i][j]);v=J[p,i,j] if name=='J' else AB[p,i,j,name]
   diff=v-ref
   if name=='J':
    u=207 if p==0 else 68
    tol=1100*(eb_old+eb_new)+u*(2*(eb_old+eb_new)+eb_old**2+eb_new**2)+arb(2)**-2000
    print('J comparison',p,i,j,'own',v,'source',ref,'difference',diff,'bound',tol,flush=True)
    assert abs(diff)<tol
   else:
    tol=arb(0);assert v.overlaps(ref)
   records.append({'parity':p,'matrix':name,'row':i,'column':j,
    'own_rounding_ball':str(v),'source_rounding_ball':str(ref),
    'difference_ball':str(diff),'ratio_ball':str(v/ref),
    'combined_quadrature_bound':str(tol),
    'difference_over_quad_bound':str(abs(diff)/tol) if name=='J' else None,
    'raw_balls_overlap':v.overlaps(ref),
    'status':'PASS'})
   print(p,name,i,j,'PASS',str(diff),flush=True)
data={'status':'PASS','precision_bits':ctx.prec,'frequency_q':q,
 'prime_method':'Exact polynomial antiderivatives on physical activation cells',
 'pole_method':'Exponential moment antiderivative recurrence',
 'frequency_method':'Independent entry integrals; Arb Gauss q=320 vs source q=256',
 'source_delta':'d(c), not target 2^-1515',
 'records':records,'new_band_entry_error':str(eb_new),'old_band_entry_error':str(eb_old),
 'seconds':time.monotonic()-t0,'maxrss_bytes_macos':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(OUT/'source_spotchecks.json').write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
print('COMPLETE',data['seconds'],data['maxrss_bytes_macos'])
