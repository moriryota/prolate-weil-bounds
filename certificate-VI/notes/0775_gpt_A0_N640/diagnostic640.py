"""LDL/congruence midpoint diagnostics. No integral certificate is claimed."""
import json,time,resource,gzip,hashlib,gc
from pathlib import Path
import mpmath as mp
from flint import arb,arb_mat,ctx

def tick(label):
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(label,'peak_RSS_MB',round(rss/1e6,1),flush=True)
    if rss>8_000_000_000:raise MemoryError('0775 conservative 8GB peak-RSS stop')

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

def parameters():
    m=arb(1)/2;r=arb(2)**-190;p=arb(2)**-440
    e=1500*r+p;n=896*r*r+p;h=12000*r*r+p
    assert n<m
    alpha=e+2*e*e/(m-n);beta=(e+n)/(m*m)+2*h*h/(m*m*(m-n))
    assert alpha+beta*(arb(13)/2)**2<arb(2)**-170
    return alpha,beta,{'r':'2^-190','p':'2^-440','m':'1/2','mu0':'1/2','theta':1,'chi':1,'n':str(n),'n_lt_m':True,'alpha':str(alpha),'beta':str(beta),'total_scalar_budget_lt_2_pow_minus170':True}

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

def diagnose(A,B,J,parity,out,tag):
    start=time.monotonic();bits=ctx.prec;d=A.nrows();I=arb_mat(d,d)
    for i in range(d):I[i,i]=1
    # Symmetrization retains the exact symmetric source matrices.
    A=(A+A.transpose())/2;B=(B+B.transpose())/2;J=(J+J.transpose())/2
    m=arb(1)/2;D=B-A*A;C=B-m*A
    tick('inverse C start')
    X=C.solve(A-m*I)
    W=(I-(A-m*I)*X)/m;W=(W+W.transpose())/2
    Winv=W.inv();Winv=(Winv+Winv.transpose())/2
    alpha,beta,pars=parameters();T=Winv+J-alpha*I-beta*D;T=(T+T.transpose())/2
    q=quality(T);assert q['row_radius_below_1e_minus48'], 'T rounding too wide for meaningful target scale'
    tick('paid T formed')
    tref=save_matrix(out/f'{tag}_p{parity}_T.json.gz',T)
    lo,piv=midpoint_ldl(T)
    negative=sum(v<0 for v in piv)
    rec={'parity':parity,'dimension':d,'parameters':pars,'T_source_rounding_quality':q,
         'LDL_arithmetic_bits':512,'negative_midpoint_pivots':negative,'positive_midpoint_pivots':sum(v>0 for v in piv),
         'minimum_midpoint_pivot':mp.nstr(min(piv),60),'T_file':tref}
    (out/f'{tag}_p{parity}_pivots.json').write_text(json.dumps([mp.nstr(v,170) for v in piv]))
    if negative==0:
        # Use LDL to propose an inverse-square-root congruence, then check
        # its residual against the interval-rounded quadrature matrix T.
        lm=arb_mat([[arb(mp.nstr(v,170)) for v in row] for row in lo])
        inv=lm.inv().transpose()
        for j in range(d):
            scale=arb(mp.nstr(piv[j],170)).sqrt()
            for i in range(d):inv[i,j]/=scale
        whiten=inv.transpose()*T*inv
        eta=max(sum(abs(whiten[i,j]-(i==j)) for j in range(d)) for i in range(d))
        bf=sum(inv[i,j]**2 for i in range(d) for j in range(d))
        rec.update(congruence_residual_row_bound=str(eta),congruence_residual_lt_one=bool(eta<1),
                   frobenius_squared=str(bf),quadrature_surrogate_lower_bound=str((1-eta)/bf),
                   positive_scope='Interval-rounded quadrature surrogate, with analytic block errors included; frequency/pole integration remainder still unpaid')
        del inv,whiten,lm;gc.collect()
    # One approximate eigendecomposition (not repeated for six auxiliaries).
    # Fast FLINT arithmetic replaces Python mpmath full spectra.
    tick('approximate minimum direction start')
    ctx.prec=512
    TM=T.mid();eigs,vecs=TM.eig(right=True,algorithm='approx')
    k=min(range(d),key=lambda j:float(eigs[j].real.mid()))
    column=[vecs[i,k] for i in range(d)]
    jmax=max(range(d),key=lambda i:float(abs(column[i]).mid()))
    rot=[z/column[jmax] for z in column]
    imaginary=max(abs(z.imag) for z in rot)
    v=arb_mat([[z.real.mid()] for z in rot]);norm=(v.transpose()*v)[0,0].sqrt();v=v/norm
    ctx.prec=bits
    ray=lambda M:(v.transpose()*M*v)[0,0]
    comp=ray(A+J);loss=ray(A-Winv);bpay=beta*ray(D);value=ray(T)
    residual=TM*v-value*v
    resnorm=sum(residual[i,0]**2 for i in range(d)).sqrt()
    rec.update(approximate_eigenvalue=str(eigs[k]),eigenvector_imaginary_rotation_max=str(imaginary),
               minimum_direction_rayleigh=str(value),eigenvector_residual_norm=str(resnorm),
               decomposition={'compression_M_plus_V':str(comp),'inverse_transfer_loss':str(loss),'alpha_payment':str(alpha),'beta_D_payment':str(bpay),'sum_check':str(comp-loss-alpha-bpay-value)})
    (out/f'{tag}_p{parity}_direction.json').write_text(json.dumps([str(v[i,0]) for i in range(d)]))
    rec['seconds']=time.monotonic()-start
    (out/f'{tag}_p{parity}_diagnostic.json').write_text(json.dumps(rec,indent=2))
    tick('parity complete '+str(parity))
    return rec
