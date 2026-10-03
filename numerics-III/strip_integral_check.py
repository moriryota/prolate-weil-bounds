"""Paper III, Lemma 6.1: numerical comparison of the strip integral with its bound (mpmath quadrature; not part of the proof).
Usage: python numerics-III/strip_integral_check.py"""
import mpmath as m,json
m.mp.dps=100
b=m.mpf('11.15');T0=m.mpf(3)*10**12
rows=[]
for logmu in [m.log(50),m.log(1000),m.log(10**6),m.log(10**12),m.mpf(1000),m.mpf(10)**6]:
 x=m.log(3)+max(m.log(T0),3*(logmu+m.log(2*m.pi)))
 r=1/(b*x);delta=logmu*r;upper=logmu*(1-r)
 cuts=[delta]+[m.mpf(t) for t in [1,3,10,30,100] if delta<t<upper]+[upper]
 actual=m.fsum(m.quad(lambda w:m.exp(-w)/w,[aa,bb]) for aa,bb in zip(cuts,cuts[1:]))/2
 e1=(m.e1(delta)-m.e1(upper))/2
 bound=(m.log(b*x/logmu)+m.exp(-1))/2
 assert actual<bound and abs(actual-e1)<m.mpf('1e-90')
 rows.append({'log_mu':m.nstr(logmu,30),'delta':m.nstr(delta,30),'integral_over_lambda_quadrature':m.nstr(actual,50),'integral_over_lambda_e1':m.nstr(e1,50),'quadrature_minus_e1':m.nstr(actual-e1,20),'bound_over_lambda':m.nstr(bound,50),'actual_over_bound':m.nstr(actual/bound,40)})
print(json.dumps({'digits':m.mp.dps,'interval_certified':False,'method':'Independent scalar quadrature; not a substitute for the analytic integral inequality.','rows':rows},indent=2))
