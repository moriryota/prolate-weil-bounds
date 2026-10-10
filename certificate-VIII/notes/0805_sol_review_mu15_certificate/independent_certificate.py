"""0805 Sol: independent fixed-mu15 final arithmetic checker.

Adapted from our 0787 implementation; no candidate codec/bounds/checker/proposer
imports. One source matrix at a time. Norm bound max(1-norm,infinity-norm).
Z is untrusted input. Scalar analytic bounds were reviewed in 0802.
"""
from pathlib import Path
from fractions import Fraction as F
from flint import arb, arb_mat, ctx
import gzip, json, hashlib, time, resource, threading, os, gc

ctx.prec=1664
ctx.threads=1
OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
BASE=ROOT/'notes/0803_gpt_A5_mu15_certificate'
DIM=832
started=time.monotonic()
def memory_guard():
    while True:
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>8_500_000_000:
            print('REJECT memory threshold 8.5GB',flush=True)
            os._exit(3)
        time.sleep(.5)
threading.Thread(target=memory_guard,daemon=True).start()

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''): h.update(chunk)
    return h.hexdigest()
def lift(q):
    q=F(q)
    return arb(q.numerator)/q.denominator
def textnum(q): return q.str(40)
class Reject(Exception): pass
def need(test,reason):
    if not test: raise Reject(reason)

cfg=json.loads((BASE/'packet/manifest.json').read_text())
expected={'schema':'A5-0803-v1','N':1664,'Omega':736,'mu':15,'m':'1/2',
          'gamma':'37/8','b':'35/4','u':[428,142],'delta':'2^-2901','tau':'2^-2904'}
need({k:cfg[k] for k in expected}==expected,'E_CONTRACT')
pins=json.loads((BASE/'SOURCE_PINS.json').read_text())
need(sha(BASE/'SOURCE_PINS.json')==cfg['source_pins_sha256'],'E_SOURCE_PINS')
need((pins['N'],pins['bits'],pins['q'])==(1664,3328,224),'E_SOURCE_PARAMETERS')
need([v['parity'] for v in cfg['blocks']]==[0,1],'E_PARITIES')
ends=[arb(0),arb(1)/2]+[arb(2)**k for k in range(7)]+[arb(k) for k in range(128,736,64)]+[arb(736)]
checked_hashes={}
for ref in pins['manifests']:
    path=ROOT/ref['path'];need(sha(path)==ref['sha256'],'E_MANIFEST_HASH')
    checked_hashes[ref['path']]=ref['sha256']
    data=json.loads(path.read_text())
    need((data['mu'],data['Omega'],data['N'],data['bits'],data['q'])==(15,736,1664,3328,224),'E_SOURCE_PARAMETERS')
    need(data['u']==[428,142] and arb(data['gamma'])==arb(37)/8,'E_SOURCE_PARAMETERS')
    need(len(data['panel_ends'])==len(ends) and all(arb(x)==y for x,y in zip(data['panel_ends'],ends)),'E_PANELS')
for pname in ['SOURCE_CODE_PINS.json','REVIEW_0802_PINS.json']:
    for name,digest in json.loads((BASE/pname).read_text()).items():
        need(sha(ROOT/name)==digest,'E_DEPENDENCY_HASH')
        checked_hashes[name]=digest

def load_exact(ref):
    name=ref['file'];path=BASE/'packet'/name
    need(path.name==name and sha(path)==ref['sha256'],'E_PACKET_HASH')
    checked_hashes[str(path.relative_to(ROOT))]=ref['sha256']
    with gzip.open(path,'rt') as f: data=json.load(f)
    need(type(data['denominator_power']) is int and data['denominator_power']==512,'E_DENOMINATOR')
    rows=data['rows']
    need(type(rows) is list and len(rows)==DIM and all(type(r) is list and len(r)==DIM for r in rows),'E_SHAPE')
    mat=arb_mat(DIM,DIM)
    for i,row in enumerate(rows):
        for j,string in enumerate(row):
            need(type(string) is str,'E_INTEGER_TEXT')
            integer=int(string)
            need(str(integer)==string and abs(integer).bit_length()<=ctx.prec,'E_EXACT_INPUT')
            mat[i,j]=arb(integer)*arb(2)**-512
            need(mat[i,j].is_exact() and mat[i,j].is_finite(),'E_EXACT_INPUT')
    return mat
def load_source(p,name):
    ref=pins['files'][f'p{p}_{name}'];path=ROOT/ref['path']
    need(sha(path)==ref['sha256'],'E_SOURCE_HASH')
    checked_hashes[ref['path']]=ref['sha256']
    with gzip.open(path,'rt') as f:rows=json.load(f)
    need(len(rows)==DIM and all(len(r)==DIM for r in rows),'E_SOURCE_SHAPE')
    mat=arb_mat(DIM,DIM)
    for i,row in enumerate(rows):
        for j,value in enumerate(row):
            mat[i,j]=arb(value)
            need(mat[i,j].is_finite(),'E_SOURCE_FINITE')
    return mat
def upper_norm(M):
    rows=[arb(0) for _ in range(M.nrows())]
    cols=[arb(0) for _ in range(M.ncols())]
    for i in range(M.nrows()):
        for j in range(M.ncols()):
            value=abs(M[i,j]).upper()
            rows[i]+=value;cols[j]+=value
    value=max(v.upper() for v in rows+cols)
    need(value.is_finite(),'E_NORM_FINITE')
    return value
def witness_size(Z):
    return sum((abs(Z[i,j]).upper()**2 for i in range(DIM) for j in range(DIM)),arb(0)).upper()
I=arb_mat(DIM,DIM)
for i in range(DIM): I[i,i]=1
def positive(M,Z,error,stage):
    eta=upper_norm(Z.transpose()*(M*Z)-I)
    size=witness_size(Z);need(size>0,'E_WITNESS')
    paid=(eta+error*size).upper()
    need(paid<1,'E_POSITIVITY_'+stage)
    lower=((1-paid)/size).lower();need(lower>0,'E_POSITIVITY_'+stage)
    return {'residual':textnum(eta),'witness_size':textnum(size),
            'paid_residual':textnum(paid),'lower':textnum(lower)}
budgets={key:F(1,2**exponent) for key,exponent in [('A',400),('B',400),('J',310),('W',350),('D',380)]}
need({k:F(v) for k,v in cfg['errors'].items()}==budgets,'E_FIXED_ERRORS')
ea,eb,ej,ew,ed=[lift(budgets[k]) for k in ['A','B','J','W','D']]
e,n,h=F(1,2**706),F(1,2**1424),F(1,2**1419)
m=F(1,2);b=F(35,4)
alpha_exact=e+2*e*e/(m-n)
beta_exact=(e+n)/(m*m)+2*h*h/(m*m*(m-n))
need(alpha_exact+beta_exact*b*b<F(1,2**697),'E_ALPHA_BETA')
alpha,beta=lift(alpha_exact),lift(beta_exact)
qerror=lift(F(1,2**319));delta_error=lift(F(1,2**2889))
def source_gate(M,S,name):
    center=upper_norm(M-S)
    analytic=qerror+delta_error if name=='J' else arb(0)
    total=(center+analytic).upper()
    need(total<=lift(budgets[name]),'E_SOURCE_'+name)
    return {'rounding_and_center':textnum(center),'quadrature':textnum(qerror if name=='J' else arb(0)),
            'delta_change':textnum(delta_error if name=='J' else arb(0)),
            'total':textnum(total),'budget':textnum(lift(budgets[name])),
            'ratio':textnum(total/lift(budgets[name]))}

results=[];negative_tests=[]
for block in cfg['blocks']:
    p=block['parity'];print('BEGIN parity',p,'seconds',time.monotonic()-started,flush=True)
    need(set(block['files'])=={'A','B','J','X','W0','Y','ZC','ZW','ZT'},'E_MATRIX_NAMES')
    mats={}
    for key,ref in block['files'].items():
        mats[key]=load_exact(ref)
        print('LOAD exact',p,key,flush=True)
    A,B,J,X,W,Y,ZC,ZW,ZT=[mats[k] for k in ['A','B','J','X','W0','Y','ZC','ZW','ZT']]
    for name in ['A','B','J','W0','Y']:
        need(all(mats[name][i,j]==mats[name][j,i] for i in range(DIM) for j in range(i)),'E_SYMMETRY')
    src={}
    for name in ['A','B','J']:
        S=load_source(p,name)
        src[name]=source_gate(mats[name],S,name)
        print('SOURCE PASS',p,name,'ratio',src[name]['ratio'],flush=True)
        if p==0 and name=='J':
            # Isolate the arithmetic source-error gate, without changing any input file.
            original=J[7,7];J[7,7]=original+arb(2)**-200
            try:
                source_gate(J,S,'J')
            except Reject as err:
                need(str(err)=='E_SOURCE_J','E_MUTATION_REASON')
                negative_tests.append({'case':'J_7_7_plus_2^-200','expected':'E_SOURCE_J',
                    'actual':str(err),'passed':True,'changed_input_files':False})
            else:
                raise Reject('E_MUTATION_NOT_REJECTED')
            finally:
                J[7,7]=original
            print('NEGATIVE PASS: J[7,7]+2^-200 -> E_SOURCE_J',flush=True)
        del S;gc.collect()
    Ctest=positive(B-A/2,ZC,eb+ea/2,'C')
    xn,an=upper_norm(X),upper_norm(A)
    # Rearranged polynomial, independent of candidate Wexpr operation order.
    Fmat=2*A-I
    Wcalc=2*I-Fmat*X-X.transpose()*Fmat+X.transpose()*((2*B-A)*X)
    wr=upper_norm(Wcalc-W)
    wneed=(xn*xn+4*xn)*ea+2*xn*xn*eb+wr
    dneed=eb+2*an*ea+ea*ea
    need(wneed<=ew,'E_W_ENCLOSURE');need(dneed<=ed,'E_D_ENCLOSURE')
    Wtest=positive(W,ZW,ew,'W')
    print('C/W PASS parity',p,flush=True)
    invres=upper_norm(I-(W+ew*I)*Y);need(invres<1,'E_INVERSE_RESIDUAL')
    ey=upper_norm(Y)*invres/(1-invres)
    need(ey<=arb(2)**-310,'E_INVERSE_ERROR')
    inverse_source=2*lift(b)**2*ew
    et=inverse_source+ey+ej+beta*ed
    need(et<=arb(2)**-300,'E_ERROR_BUDGET')
    T=Y+J-alpha*I-beta*(B-A*A)
    Ttest=positive(T,ZT,et,'T')
    row={'parity':p,'source_errors':src,'C':Ctest,'W':Wtest,'T':Ttest,
         'W_arithmetic_error':textnum(wr),'X_norm':textnum(xn),'A_norm':textnum(an),
         'epsilon_W_required':textnum(wneed),'W_ratio':textnum(wneed/ew),
         'epsilon_D_required':textnum(dneed),'D_ratio':textnum(dneed/ed),
         'inverse_residual':textnum(invres),'epsilon_Y':textnum(ey),
         'epsilon_inverse_source':textnum(inverse_source),'epsilon_T':textnum(et)}
    results.append(row)
    print('PASS parity',p,'T lower',Ttest['lower'],'paid',Ttest['paid_residual'],flush=True)
    del mats,A,B,J,X,W,Y,ZC,ZW,ZT,T,Fmat,Wcalc;gc.collect()
for name,digest in checked_hashes.items():need(sha(ROOT/name)==digest,'E_END_HASH')
need(len(negative_tests)>=1,'E_NEGATIVE_TEST')
result={'status':'PASS_INDEPENDENT_CERTIFICATE','reviewer':'sol','precision_bits':ctx.prec,
 'norm_bound':'max(sum absolute column entries, sum absolute row entries)',
 'alpha_exact':str(alpha_exact),'beta_exact':str(beta_exact),'fixed_tau':'2^-2904',
 'hashes_rechecked':len(checked_hashes),'results':results,'negative_tests':negative_tests,
 'seconds':time.monotonic()-started,'maxrss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
 'scope':'independent final arithmetic plus previously reviewed analytic bounds; external Arb/FLINT, not Lean'}
(OUT/'independent_certificate.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print('COMPLETE',result['seconds'],'seconds',result['maxrss_bytes'],'bytes RSS',flush=True)
