"""Untrusted proposal. The independent checker does not import this module."""
import json,sys,time,resource
from pathlib import Path
import mpmath as mp
import numpy as np
from flint import arb,arb_mat,ctx
from codec import dump_dyadic,read_dyadic,read_source,sha
from bounds import analytic_bounds,fraction_text

OUT=Path(__file__).parent;ROOT=OUT.parents[1];ctx.prec=1280;mp.mp.prec=512
def tick(s):
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(s,'RSS_MB',round(rss/1e6,1),flush=True)
    if rss>8_000_000_000:raise MemoryError('8GB checkpoint stop')
def sym(M):return (M+M.transpose())/2
def identity():
    I=arb_mat(320,320)
    for i in range(320):I[i,i]=1
    return I
def ab(q):return arb(q.numerator)/q.denominator
def double_whitener(M):
    a=np.array([[float(M[i,j].mid()) for j in range(320)] for i in range(320)])
    l=np.linalg.cholesky((a+a.T)/2)
    z=np.linalg.solve(l.T,np.eye(320))
    return arb_mat([[arb(float(x)) for x in row] for row in z])
def ldl_whitener(M):
    a=[[mp.mpf(M[i,j].mid().str(170,radius=False)) for j in range(320)] for i in range(320)]
    l=[[mp.mpf(0)]*320 for _ in range(320)];ds=[]
    for j in range(320):
        l[j][j]=mp.mpf(1);weighted=[l[j][k]*ds[k] for k in range(j)]
        ds.append(a[j][j]-sum(l[j][k]*weighted[k] for k in range(j)))
        assert ds[-1]>0
        for i in range(j+1,320):
            l[i][j]=(a[i][j]-sum(l[i][k]*weighted[k] for k in range(j)))/ds[-1]
    L=arb_mat([[arb(mp.nstr(x,170)) for x in row] for row in l])
    Z=L.inv().transpose()
    for j in range(320):
        scale=arb(mp.nstr(ds[j],170)).sqrt()
        for i in range(320):Z[i,j]/=scale
    return Z

def main():
    dest=OUT/'packet';dest.mkdir(exist_ok=True)
    assert not (dest/'manifest.json').exists()
    pins=json.loads((OUT/'SOURCE_PINS.json').read_text())
    vals=analytic_bounds();I=identity();start=time.monotonic();blocks=[]
    eps={'A':'1/'+str(2**400),'B':'1/'+str(2**400),'J':'1/'+str(2**185),
         'W':'1/'+str(2**350),'D':'1/'+str(2**380)}
    manifest={'schema':'A1-0776-v1','N':640,'Omega':256,'source_pins_sha256':sha(OUT/'SOURCE_PINS.json'),
              'errors':eps,'blocks':blocks}
    for p in (0,1):
        refs={};M={}
        def put(name,x):
            refs[name]=dump_dyadic(dest/f'p{p}_{name}.json.gz',x)
            M[name]=read_dyadic(dest/refs[name]['file'])
            return M[name]
        for name in ['A','B','J']:
            put(name,sym(read_source(ROOT,pins,p,name)))
        A,B,J=M['A'],M['B'],M['J'];C=B-A/2
        X=put('X',C.solve(A-I/2))
        W=put('W0',sym(X+X.transpose()-X.transpose()*A*X+2*(I-A*X-X.transpose()*A+X.transpose()*B*X)))
        Wup=W+arb(2)**-350*I
        Y=put('Y',sym(Wup.inv()))
        T=Y+J-ab(vals['alpha'])*I-ab(vals['beta'])*(B-A*A)
        put('ZC',double_whitener(C));put('ZW',double_whitener(W))
        put('ZT',ldl_whitener(T))
        blocks.append({'parity':p,'files':refs});tick('proposed parity '+str(p))
    (dest/'manifest.json').write_text(json.dumps(manifest,indent=2))
    (OUT/'proposal_result.json').write_text(json.dumps({'seconds':time.monotonic()-start,
        'maxrss_bytes_macos':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'manifest_sha256':sha(dest/'manifest.json'),'scope':'UNTRUSTED_PROPOSAL'},indent=2))
    print('proposal complete',flush=True)
if __name__=='__main__':main()
