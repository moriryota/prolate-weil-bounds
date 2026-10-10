"""0805 Sol: selected mu15 source entries from different formulas.

A/B: polynomial antiderivatives (no Gauss); high pair (400,401) included.
J: q256 instead of q224; poles via modified Bessel I integral identity.
No candidate/project modules imported. Arb balls and analytic remainders kept separate.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb, factorial
from flint import arb, acb, arb_poly, ctx
import gzip, json, hashlib, time, resource, gc, threading, os

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
ctx.prec=8192;ctx.threads=1
started=time.monotonic()
def guard():
    while True:
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>8_500_000_000:
            print('REJECT memory threshold 8.5GB',flush=True);os._exit(3)
        time.sleep(.5)
threading.Thread(target=guard,daemon=True).start()
def ss(v):return v.str(40)
def aq(v):return arb(v.numerator)/v.denominator
L=arb(15).log()/2;gamma=arb(37)/8;pi=arb.pi()
primes={2:2,3:3,4:2,5:5,7:7,8:2,9:3,11:11,13:13}
ABpairs={0:[(7,4),(400,401)],1:[(3,0)]}
Jpairs={0:[(7,4),(400,401)],1:[(3,0),(400,401)]}
ABdegrees=sorted({2*k+p for p,pairs in ABpairs.items() for pair in pairs for k in pair})
def normalized_polynomial(n):
    coeff=[arb(0)]*(n+1)
    for k in range(n//2+1):
        coeff[n-2*k]=arb((-1)**k*comb(n,k)*comb(2*n-2*k,n))*arb(2)**-n
    return arb_poly(coeff)*(arb(2*n+1)/(2*L)).sqrt()
phi={n:normalized_polynomial(n) for n in ABdegrees}
facts=[factorial(k) for k in range(max(ABdegrees)+1)]
def translate(poly,h):
    # Taylor shift by one polynomial convolution: coefficient_j=sum a_k binom(k,j) h^(k-j).
    n=len(poly)-1
    rev=arb_poly([poly[n-k]*facts[n-k] for k in range(n+1)])
    powers=[];power=arb(1)
    for k in range(n+1):powers.append(power/facts[k]);power*=h
    product=rev*arb_poly(powers)
    return arb_poly([product[n-j]/facts[j] for j in range(n+1)])
translations={(n,q,sign):translate(phi[n],sign*arb(q).log()/L)
              for n in ABdegrees for q in primes for sign in (-1,1)}
# Derive exact activation q thresholds and physical endpoints. No candidate cells imported.
events=sorted({F(1),F(15),*[max(F(15,q*q),F(q*q,15)) for q in primes]})
AB={(p,i,j,name):arb(0) for p,pairs in ABpairs.items() for i,j in pairs for name in ('A','B')}
without13=AB.copy()
def integ(poly,a,b):
    primitive=poly.integral()
    return 2*L*(primitive(b)-primitive(a))  # positive half, parity product even
for leftq,rightq in zip(events,events[1:]):
    a=aq(leftq).log()/(2*L);b=aq(rightq).log()/(2*L)
    mid=(leftq+rightq)/2
    active=[(q,sign) for q in primes for sign in (-1,1)
            if (sign==1 and mid<F(15,q*q)) or (sign==-1 and mid>F(q*q,15))]
    full={n:gamma*phi[n] for n in ABdegrees};no13={n:gamma*phi[n] for n in ABdegrees}
    for q,sign in active:
        w=arb(primes[q]).log()/arb(q).sqrt()
        for n in ABdegrees:
            term=w*translations[n,q,sign]
            full[n]-=term
            if q!=13:no13[n]-=term
    for p,pairs in ABpairs.items():
        for i,j in pairs:
            n,r=2*i+p,2*j+p
            AB[p,i,j,'A']+=integ(phi[n]*full[r],a,b)
            AB[p,i,j,'B']+=integ(full[n]*full[r],a,b)
            without13[p,i,j,'A']+=integ(phi[n]*no13[r],a,b)
            without13[p,i,j,'B']+=integ(no13[n]*no13[r],a,b)
    print('PRIME CELL',str(leftq),str(rightq),time.monotonic()-started,flush=True)
del phi,translations,full,no13;gc.collect()
for name in ('A','B'):
    contribution=AB[0,7,4,name]-without13[0,7,4,name]
    assert not contribution.contains(0),name+' has no resolved prime13 contribution'
    print('PRIME13',name,ss(contribution),flush=True)

# Distinct Gauss order, scalar entries only; special-function precision differs from source.
ctx.prec=2304
L=arb(15).log()/2;pi=arb.pi();Omega=736;c=Omega*L
q=256
nodes=[arb.legendre_p_root(q,k,weight=True) for k in range(q)]
degrees={p:sorted({p,*[2*k+p for pair in pairs for k in pair]}) for p,pairs in Jpairs.items()}
all_degrees=sorted(set(degrees[0]+degrees[1]))
ends=[arb(0),arb(1)/2]+[arb(2)**k for k in range(7)]+[arb(k) for k in range(128,736,64)]+[arb(736)]
KB={(p,i,j):arb(0) for p,pairs in Jpairs.items() for i,j in pairs}
cols={(p,n):arb(0) for p in degrees for n in degrees[p]}
for lo,hi in zip(ends,ends[1:]):
    half=(hi-lo)/2;center=(hi+lo)/2
    for xi,wi in nodes:
        t=center+half*xi;z=L*t
        vals={n:(-1)**(n//2)*(2*L*(2*n+1)).sqrt()*(pi/(2*z)).sqrt()*z.bessel_j(arb(n)+arb(1)/2) for n in all_degrees}
        weight=half*wi/pi
        arch=acb(arb(1)/4,t/2).digamma().real-pi.log()-gamma
        for p,pairs in Jpairs.items():
            for i,j in pairs:KB[p,i,j]+=weight*arch*vals[2*i+p]*vals[2*j+p]
            for n in degrees[p]:cols[p,n]+=weight*vals[n]*vals[p]
    print('FREQUENCY PANEL',ss(lo),ss(hi),time.monotonic()-started,flush=True)
# Rodrigues + modified-Bessel integral gives the exact exponential moment, no pole quadrature.
z=L/2
polecoef={n:L*(arb(2*n+1)/(2*L)).sqrt()*(2*pi/z).sqrt()*z.bessel_i(arb(n)+arb(1)/2) for n in all_degrees}
delta_old=(-2*(c+1)).exp()/(4*(c+1)*(c+4))
J={}
for p,pairs in Jpairs.items():
    for i,j in pairs:
        n,r=2*i+p,2*j+p
        hi=(1-delta_old if i==0 else arb(0))-cols[p,n]
        hj=(1-delta_old if j==0 else arb(0))-cols[p,r]
        pole=2*(-1)**p*polecoef[n]*polecoef[r]
        J[p,i,j]=KB[p,i,j]+pole+(428 if p==0 else 142)*hi*hj

# Entire full-dimension quadrature bounds dominate selected-entry errors.
HB=2*L*(48*L).exp()
eb=832*4*736*HB*arb(2)**-511
own_quad=1100*eb+428*(2*eb+eb*eb)
assert own_quad<arb(2)**-383
old_quad=arb(2)**-319
tol=(own_quad+old_quad).upper()
pins=json.loads((ROOT/'notes/0803_gpt_A5_mu15_certificate/SOURCE_PINS.json').read_text())
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()
records=[]
for p in (0,1):
    for name in ('A','B','J'):
        ref=pins['files'][f'p{p}_{name}'];path=ROOT/ref['path']
        assert digest(path)==ref['sha256']
        with gzip.open(path,'rt') as f:source=json.load(f)
        assert len(source)==832 and all(len(row)==832 for row in source)
        pairs=Jpairs[p] if name=='J' else ABpairs[p]
        for i,j in pairs:
            given=arb(source[i][j]);own=J[p,i,j] if name=='J' else AB[p,i,j,name]
            difference=own-given
            if name=='J':assert abs(difference)<tol
            else:assert own.overlaps(given)
            row={'parity':p,'matrix':name,'row':i,'column':j,'degrees':[2*i+p,2*j+p],
                 'own':ss(own),'source':ss(given),'difference':ss(difference),'ratio':ss(own/given),
                 'raw_balls_overlap':own.overlaps(given),'analytic_difference_bound':ss(tol) if name=='J' else '0',
                 'difference_over_bound':ss(abs(difference)/tol) if name=='J' else None,'status':'PASS'}
            if name!='J':row['prime13_increment']=ss(AB[p,i,j,name]-without13[p,i,j,name])
            records.append(row)
            print('ENTRY PASS',p,name,i,j,'own',row['own'],'source',row['source'],'difference',row['difference'],'ratio',row['ratio'],flush=True)
        del source;gc.collect()
result={'status':'PASS_INDEPENDENT_SOURCE_SPOTCHECKS','reviewer':'sol',
 'prime_precision_bits':8192,'frequency_precision_bits':2304,'frequency_q':256,
 'prime_method':'Legendre coefficient formula, factorial-convolution Taylor shifts, polynomial antiderivatives',
 'pole_method':'Exact modified-Bessel-I exponential moment; cosh/sinh by parity',
 'source_delta':'d(736L), not target 2^-2901',
 'own_frequency_operator_error':ss(own_quad),'source_full_J_quadrature_error':'2^-319',
 'records':records,'seconds':time.monotonic()-started,'maxrss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(OUT/'source_spotchecks.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print('COMPLETE',len(records),'entries',result['seconds'],'seconds',result['maxrss_bytes'],'bytes RSS',flush=True)
