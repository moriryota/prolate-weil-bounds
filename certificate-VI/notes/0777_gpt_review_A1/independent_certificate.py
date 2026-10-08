"""Reviewer implementation. No imports from proposal/checker/codec/bounds.
1536-bit Arb; max(row sum,column sum) replaces geometric norm bound.
Both parity blocks, all source errors, C/W positivity, inverse and final T.
"""
import gzip,json,hashlib,time,resource
from pathlib import Path
from fractions import Fraction as F
from math import factorial
from flint import arb,arb_mat,ctx
ctx.prec=1536
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]
cfg=json.loads((HERE/'packet/manifest.json').read_text())
pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
def checksum(path):return hashlib.sha256(path.read_bytes()).hexdigest()
assert checksum(HERE/'SOURCE_PINS.json')==cfg['source_pins_sha256']
assert checksum(ROOT/pins['manifest']['path'])==pins['manifest']['sha256']
def rational(x):
    x=F(x);return arb(x.numerator)/x.denominator
def matrix(ref):
    path=HERE/'packet'/ref['file'];assert checksum(path)==ref['sha256']
    data=json.loads(gzip.decompress(path.read_bytes())); power=data['denominator_power']
    assert len(data['rows'])==320 and all(len(row)==320 for row in data['rows'])
    rows=[[int(s) for s in row] for row in data['rows']]
    assert all(abs(x).bit_length()<=1536 for row in rows for x in row)
    return arb_mat([[arb(x)*arb(2)**(-power) for x in row] for row in rows])
def source(parity,name):
    ref=pins['files'][f'p{parity}_{name}'];path=ROOT/ref['path']
    assert checksum(path)==ref['sha256']
    return arb_mat([[arb(s) for s in row] for row in json.loads(gzip.decompress(path.read_bytes()))])
def norm(M):
    rows=[sum((abs(M[i,j]).upper() for j in range(320)),arb(0)).upper() for i in range(320)]
    cols=[sum((abs(M[i,j]).upper() for i in range(320)),arb(0)).upper() for j in range(320)]
    return max(rows+cols)
I=arb_mat([[int(i==j) for j in range(320)] for i in range(320)])
def gate(M,Z,error):
    eta=norm(Z.transpose()*M*Z-I)
    s=sum((abs(Z[i,j]).upper()**2 for i in range(320) for j in range(320)),arb(0)).upper()
    paid=(eta+error*s).upper();assert paid<1
    return {'eta':str(eta),'s':str(s),'paid':str(paid),'lower':str(((1-paid)/s).lower())}
eb=F(320*1024*3*10**48,2**383)
eq=F(320*1024*4000*10**48,2**383)+81*(2*eb+eb*eb)+320*F(3,5)**721/factorial(721)
r,p,m=F(1,2**190),F(1,2**440),F(1,2)
e=1500*r+p;n=896*r*r+p;h=12000*r*r+p
alpha=rational(e+2*e*e/(m-n));beta=rational((e+n)/(m*m)+2*h*h/(m*m*(m-n)))
ea,ebudget,ej,ew,ed=[rational(cfg['errors'][k]) for k in ['A','B','J','W','D']]
t0=time.monotonic();results=[]
for block in cfg['blocks']:
    parity=block['parity'];M={name:matrix(ref) for name,ref in block['files'].items()}
    A,B,J,X,W,Y,ZC,ZW,ZT=[M[k] for k in ['A','B','J','X','W0','Y','ZC','ZW','ZT']]
    for mat in [A,B,J,W,Y]:assert all(mat[i,j]==mat[j,i] for i in range(320) for j in range(i))
    se={}
    for name,budget in [('A',ea),('B',ebudget),('J',ej)]:
        se[name]=norm(M[name]-source(parity,name))+(rational(eq) if name=='J' else 0)
        assert se[name]<=budget
    cc=gate(B-A/2,ZC,ebudget+ea/2)
    xn,an=norm(X),norm(A)
    WX=X+X.transpose()-X.transpose()*A*X+2*(I-A*X-X.transpose()*A+X.transpose()*B*X)
    we=(xn*xn+4*xn)*ea+2*xn*xn*ebudget+norm(WX-W)
    de=ebudget+2*an*ea+ea*ea
    assert we<=ew and de<=ed
    wc=gate(W,ZW,ew)
    residual=norm(I-(W+ew*I)*Y);assert residual<1
    iy=norm(Y)*residual/(1-residual)
    et=2*rational(F(13,2))**2*ew+iy+ej+beta*ed
    tc=gate(Y+J-alpha*I-beta*(B-A*A),ZT,et)
    lower=arb(tc['lower']);assert lower>arb(2)**(-144)
    results.append({'parity':parity,'source_errors':{k:str(v) for k,v in se.items()},'C':cc,'W':wc,
                    'W_required':str(we),'D_required':str(de),'inverse_residual':str(residual),
                    'epsilon_Y':str(iy),'epsilon_T':str(et),'T':tc,
                    'T_lower_over_2^-144':str(lower/(arb(2)**(-144)))})
    print('parity',parity,'PASS',tc['lower'],flush=True)
out={'status':'PASS','precision_bits':1536,'norm':'max(1-norm,infinity-norm)',
     'results':results,'seconds':time.monotonic()-t0,'rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
assert out['rss_bytes']<8_000_000_000
(HERE/'independent_results.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
