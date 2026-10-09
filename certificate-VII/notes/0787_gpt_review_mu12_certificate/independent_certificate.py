"""0787 reviewer certificate verification, written independently.

No candidate/prior codec, bounds, checker or proposal imports.
1664-bit Arb, exact dyadic parsing, max(1-norm,infinity-norm) bounds.
W is recomputed through a rearranged polynomial expression.
"""
from pathlib import Path
from fractions import Fraction as F
import gzip,json,hashlib,time,resource,threading,os
from flint import arb,arb_mat,ctx
ctx.prec=1664;ctx.threads=1
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
BASE=ROOT/'notes/0785_gpt_A4_mu12_certificate';DIM=416
tstart=time.monotonic()
def memory_guard():
 while True:
  if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>8_500_000_000:
   print('REJECT memory threshold 8.5GB',flush=True);os._exit(3)
  time.sleep(.5)
threading.Thread(target=memory_guard,daemon=True).start()
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
fixed=json.loads((BASE/'MATRIX_SHA256.json').read_text())
for p,v in fixed.items():assert digest(ROOT/p)==v,p
cfg=json.loads((BASE/'packet/manifest.json').read_text())
assert {k:cfg[k] for k in ['schema','N','Omega','mu','m','gamma','b','u','delta','tau']}=={
 'schema':'A4-0785-v1','N':832,'Omega':416,'mu':12,'m':'1/2','gamma':'65/16',
 'b':'61/8','u':[207,68],'delta':'2^-1515','tau':'2^-1518'}
assert digest(BASE/'SOURCE_PINS.json')==cfg['source_pins_sha256']
pins=json.loads((BASE/'SOURCE_PINS.json').read_text())
assert (pins['N'],pins['bits'],pins['q'])==(832,1920,256)
assert digest(ROOT/pins['manifest']['path'])==pins['manifest']['sha256']
hashes_before={p:digest(ROOT/p) for p in fixed}
def lift(q):
 q=F(q);return arb(q.numerator)/q.denominator
def load_exact(ref):
 path=BASE/'packet'/ref['file'];assert path.name==ref['file'] and digest(path)==ref['sha256']
 j=json.loads(gzip.decompress(path.read_bytes()));assert j['denominator_power']==512
 assert len(j['rows'])==DIM and all(len(r)==DIM for r in j['rows'])
 nums=[[int(s) for s in r] for r in j['rows']]
 assert all(abs(z).bit_length()<=ctx.prec for r in nums for z in r)
 mat=arb_mat([[arb(z)*arb(2)**-512 for z in r] for r in nums])
 assert all(mat[i,k].is_exact() and mat[i,k].is_finite() for i in range(DIM) for k in range(DIM))
 return mat
def load_source(p,k):
 ref=pins['files'][f'p{p}_{k}'];path=ROOT/ref['path'];assert digest(path)==ref['sha256']
 obj=json.loads(gzip.decompress(path.read_bytes()));assert len(obj)==DIM and all(len(r)==DIM for r in obj)
 mat=arb_mat([[arb(s) for s in r] for r in obj])
 assert all(mat[i,j].is_finite() for i in range(DIM) for j in range(DIM))
 return mat
def upper_norm(M):
 # max(||M||_1,||M||_inf) >= sqrt(||M||_1 ||M||_inf) >= ||M||_2.
 rows=[arb(0) for _ in range(DIM)];cols=[arb(0) for _ in range(DIM)]
 for i in range(DIM):
  for j in range(DIM):
   v=abs(M[i,j]).upper();rows[i]+=v;cols[j]+=v
 val=max(v.upper() for v in rows+cols)
 assert val.is_finite();return val
def witness_size(Z):return sum((abs(Z[i,j]).upper()**2 for i in range(DIM) for j in range(DIM)),arb(0)).upper()
I=arb_mat(DIM,DIM)
for i in range(DIM):I[i,i]=1
def positivity(M,Z,error,stage):
 residual=upper_norm(Z.transpose()*(M*Z)-I)
 size=witness_size(Z);assert size>0
 paid=(residual+error*size).upper();assert paid<1,stage
 bound=((1-paid)/size).lower()
 assert bound>0
 return {'residual':str(residual),'witness_size':str(size),'paid_residual':str(paid),'lower':str(bound)}
budgets={'A':F(1,2**400),'B':F(1,2**400),'J':F(1,2**250),'W':F(1,2**350),'D':F(1,2**380)}
assert {k:F(v) for k,v in cfg['errors'].items()}==budgets
ea,eb,ej,ew,ed=[lift(budgets[k]) for k in ['A','B','J','W','D']]
e,n,h=F(1,2**311),F(1,2**634),F(1,2**629);m=F(1,2)
alpha=lift(e+2*e*e/(m-n));beta=lift((e+n)/(m*m)+2*h*h/(m*m*(m-n)))
jquad=lift(F(1,2**263));jdelta=lift(F(1,2**1504))
references=json.loads((BASE/'verification.json').read_text())
assert [b['parity'] for b in cfg['blocks']]==[0,1]
results=[]
for block in cfg['blocks']:
 p=block['parity'];print('BEGIN parity',p,'seconds',time.monotonic()-tstart,flush=True)
 assert set(block['files'])=={'A','B','J','X','W0','Y','ZC','ZW','ZT'}
 mats={k:load_exact(ref) for k,ref in block['files'].items()}
 A,B,J,X,W,Y,ZC,ZW,ZT=[mats[k] for k in ['A','B','J','X','W0','Y','ZC','ZW','ZT']]
 for k in ['A','B','J','W0','Y']:
  assert all(mats[k][i,j]==mats[k][j,i] for i in range(DIM) for j in range(i)),k
 source_errors={}
 for k,budget in [('A',ea),('B',eb),('J',ej)]:
  rounding=upper_norm(mats[k]-load_source(p,k));total=rounding+(jquad+jdelta if k=='J' else 0)
  assert total<=budget,k
  source_errors[k]={'rounding':str(rounding),'quad':str(jquad if k=='J' else arb(0)),
                    'delta':str(jdelta if k=='J' else arb(0)),'total':str(total),'budget':str(budget),'ratio':str(total/budget)}
 Cgate=positivity(B-A/2,ZC,eb+ea/2,'C')
 print('SOURCE/C PASS parity',p,flush=True)
 xnorm,anorm=upper_norm(X),upper_norm(A)
 # 2I -(2A-I)X -X^T(2A-I)+X^T(2B-A)X.
 Fmat=2*A-I
 Wexpr=2*I-Fmat*X-X.transpose()*Fmat+X.transpose()*((2*B-A)*X)
 wround=upper_norm(Wexpr-W)
 wneed=(xnorm*xnorm+4*xnorm)*ea+2*xnorm*xnorm*eb+wround
 dneed=eb+2*anorm*ea+ea*ea
 assert wneed<=ew and dneed<=ed
 Wgate=positivity(W,ZW,ew,'W')
 print('W PASS parity',p,flush=True)
 invres=upper_norm(I-(W+ew*I)*Y);assert invres<1
 ey=upper_norm(Y)*invres/(1-invres)
 et=2*lift(F(61,8))**2*ew+ey+ej+beta*ed
 assert et<=arb(2)**-240
 T=Y+J-alpha*I-beta*(B-A*A)
 Tgate=positivity(T,ZT,et,'T')
 own=arb(Tgate['lower']);given=arb(references['results'][p]['T']['positive_lower_bound'])
 assert own>arb(2)**-180
 row={'parity':p,'source_errors':source_errors,'C':Cgate,'W':Wgate,
      'W_required':str(wneed),'W_budget':str(ew),'W_ratio':str(wneed/ew),
      'D_required':str(dneed),'D_budget':str(ed),'D_ratio':str(dneed/ed),
      'inverse_residual':str(invres),'epsilon_Y':str(ey),'epsilon_T':str(et),
      'epsilon_T_budget':str(arb(2)**-240),'T':Tgate,
      'comparison':{'own_lower':str(own),'candidate_lower':str(given),'ratio':str(own/given),'difference':str(own-given)}}
 results.append(row)
 print('PASS parity',p,'own lower',own,'candidate lower',given,'ratio',own/given,flush=True)
 del mats,A,B,J,X,W,Y,ZC,ZW,ZT,T,Fmat,Wexpr
for p,v in hashes_before.items():assert digest(ROOT/p)==v,p
result={'status':'PASS_INDEPENDENT_CERTIFICATE','precision_bits':ctx.prec,'norm_bound':'max(row sum,column sum)',
        'all_matrix_hashes_passed':len(fixed),'results':results,'seconds':time.monotonic()-tstart,
        'maxrss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'scope':'independent final arithmetic check, analytic scalar gates reviewed separately; not Lean'}
(OUT/'independent_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
for row in results:
 print('SUMMARY',json.dumps(row),flush=True)
print('COMPLETE',result['seconds'],'seconds',result['maxrss_bytes'],'bytes RSS',flush=True)
