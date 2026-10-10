"""LDL/congruence midpoint diagnostics. No integral certificate is claimed.
0792: copied from 0784; the peak-RSS stop is read from RSS_STOP_BYTES (default 8e9, as before)."""
import json,time,resource,gzip,hashlib,gc,os,sys
RSS_STOP=int(os.environ.get('RSS_STOP_BYTES','8000000000'))
from pathlib import Path
import mpmath as mp
from flint import arb,arb_mat,ctx

def peak_rss_bytes():
    # ru_maxrss is in bytes on macOS and in kilobytes on Linux (the VM).
    r=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r if sys.platform=='darwin' else r*1024

def tick(label):
    rss=peak_rss_bytes()
    print(label,'peak_RSS_MB',round(rss/1e6,1),flush=True)
    if rss>RSS_STOP:raise MemoryError('peak-RSS stop at %d bytes'%RSS_STOP)

def quality(M):
    d=M.nrows();rad=arb(0);sym=True
    for i in range(d):
        row=arb(0)
        for j in range(d):
            assert M[i,j].is_finite()
            row+=M[i,j].rad()
            sym=sym and (M[i,j]-M[j,i]).contains(0)
        rad=max(rad,row)
    assert sym, 'source does not enclose a symmetric matrix'
    return {'max_radius_row_sum':str(rad),'finite':True,'symmetry_enclosed':True,'row_radius_below_1e_minus48':bool(rad<arb('1e-48'))}

def save_matrix(path,M):
    assert not path.exists()
    # Stream rows, rather than retaining six large JSON string arrays.
    with gzip.open(path,'wt') as f:
        f.write('[')
        for i in range(M.nrows()):
            if i:f.write(',')
            json.dump([str(M[i,j]) for j in range(M.ncols())],f)
        f.write(']')
    return {'path':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size}

def midpoint_ldl(M,prec=512):
    mp.mp.prec=prec;d=M.nrows()
    vals=[[mp.mpf(M[i,j].mid().str(170,radius=False)) for j in range(d)] for i in range(d)]
    lower=[[mp.mpf(0)]*d for _ in range(d)];diagonal=[]
    for j in range(d):
        lower[j][j]=mp.mpf(1)
        weighted=[lower[j][k]*diagonal[k] for k in range(j)]
        pivot=vals[j][j]-sum(lower[j][k]*weighted[k] for k in range(j))
        assert pivot!=0,'Zero midpoint LDL pivot; cannot infer inertia'
        diagonal.append(pivot)
        for i in range(j+1,d):
            lower[i][j]=(vals[i][j]-sum(lower[i][k]*weighted[k] for k in range(j)))/pivot
        if j%80==79:tick('LDL column '+str(j+1))
    return lower,diagonal

