"""0784: raw T only; congruence certifies rounded surrogate, not source integrals."""
import json,time,gc
import mpmath as mp
from flint import arb,arb_mat,ctx
from matrix_utils import tick,quality,save_matrix,midpoint_ldl

def diagnose(A,B,J,parity,out,tag):
 start=time.monotonic();bits=ctx.prec;d=A.nrows();Id=arb_mat(d,d)
 for i in range(d):Id[i,i]=1
 A=(A+A.transpose())/2;B=(B+B.transpose())/2;J=(J+J.transpose())/2;m=arb(1)/2
 C=B-m*A;tick('inverse C start');X=C.solve(A-m*Id)
 W=(Id-(A-m*Id)*X)/m;W=(W+W.transpose())/2
 Wi=W.inv();Wi=(Wi+Wi.transpose())/2
 T=Wi+J;T=(T+T.transpose())/2
 qc=quality(T);assert arb(qc['max_radius_row_sum'])<arb('1e-70')
 tref=save_matrix(out/f'{tag}_p{parity}_T.json.gz',T)
 tick('raw T; start midpoint LDL')
 lo,piv=midpoint_ldl(T)
 negative=sum(v<0 for v in piv)
 rec={'parity':parity,'dimension':d,'source_quality':qc,'T':tref,'negative_midpoint_pivots':negative,'positive_midpoint_pivots':sum(v>0 for v in piv),'block_errors_paid':False,'integral_remainders_paid':False,'minimum_midpoint_pivot':mp.nstr(min(piv),55)}
 (out/f'{tag}_p{parity}_pivots.json').write_text(json.dumps([mp.nstr(v,170) for v in piv]))
 if negative==0:
  lm=arb_mat([[arb(mp.nstr(v,170)) for v in row] for row in lo]);Z=lm.inv().transpose()
  for j in range(d):
   scale=arb(mp.nstr(piv[j],170)).sqrt()
   for i in range(d):Z[i,j]/=scale
  R=Z.transpose()*T*Z-Id
  eta=max(sum(abs(R[i,j]).upper() for j in range(d)) for i in range(d))
  z2=sum(abs(Z[i,j]).upper()*abs(Z[i,j]).upper() for i in range(d) for j in range(d))
  assert z2.is_finite() and z2>0
  rec.update(congruence_residual=str(eta),congruence_lt_one=bool(eta<1),frobenius_squared=str(z2),surrogate_lower_bound=str((1-eta)/z2),positive_scope='rounded quadrature surrogate only')
  rec['Z']=save_matrix(out/f'{tag}_p{parity}_Z.json.gz',Z)
  del lm,Z,R;gc.collect()
 del lo,piv;gc.collect();tick('approx eig start')
 ctx.prec=512
 vals,vecs=T.mid().eig(right=True,algorithm='approx');k=min(range(d),key=lambda j:float(vals[j].real.mid()))
 col=[vecs[i,k] for i in range(d)];j=max(range(d),key=lambda i:float(abs(col[i]).mid()));rot=[v/col[j] for v in col]
 v=arb_mat([[z.real.mid()] for z in rot]);ctx.prec=bits;v=v/(v.transpose()*v)[0,0].sqrt()
 ray=lambda M:(v.transpose()*M*v)[0,0]
 rq=ray(T);res=T*v-rq*v
 rn=sum(abs(res[i,0]).upper()*abs(res[i,0]).upper() for i in range(d)).sqrt()
 rec.update(raw_Winv_plus_J_rayleigh=str(rq),approx_eigenvalue=str(vals[k]),eigenvector_residual_norm=str(rn),compression_on_same_vector=str(ray(A+J)),inverse_transfer_loss=str(ray(A-Wi)),D_on_same_vector=str(ray(B-A*A)),seconds=time.monotonic()-start)
 weights=[float(abs(v[i,0]).mid())**2 for i in range(d)]
 top=sorted(range(d),key=lambda j:weights[j],reverse=True)[:12]
 rec['mode_profile']={'top_modes':[[2*j+parity,weights[j]] for j in top],'mass_below_degree':{str(n):sum(w for j,w in enumerate(weights) if 2*j+parity<n) for n in (16,32,64,128,256,512)}}
 (out/f'{tag}_p{parity}_direction.json').write_text(json.dumps([str(v[i,0]) for i in range(d)]))
 (out/f'{tag}_p{parity}_diagnostic.json').write_text(json.dumps(rec,indent=2))
 print('diagnostic',parity,'rq',rq.str(18),'negative_pivots',negative,'seconds',rec['seconds'],flush=True)
 tick('parity complete');return rec
