"""Rigorous finite-polynomial Rayleigh/min-max certificate via a sinc Taylor kernel."""
from pathlib import Path
from decimal import Decimal,ROUND_CEILING,localcontext
from flint import arb,ctx
import json,math,time
out=Path(__file__).resolve().parent
ctx.prec=512;pi=arb.pi();c=10*pi;N=160
raw=json.loads((out/'trials.json').read_text())['trials']
trial=[{'degrees':r['degrees'],'a':[arb(t) for t in r['coefficients']]} for r in raw]
def zero_matrix():return [[arb(0) for _ in range(5)] for _ in range(5)]
G=zero_matrix();Q=zero_matrix()
for i in range(5):
    for j in range(i,5):
        if i%2!=j%2:continue
        G[i][j]=G[j][i]=sum((2*ai*aj/(2*l+1) for l,ai,aj in zip(trial[i]['degrees'],trial[i]['a'],trial[j]['a'])),arb(0))
# Moments ∫x^r P_l(x) dx: 0 unless r>=l and parity agrees;
# M_(r,p)=2/(r+p+1), M_(r,l+2)=M_(r,l)*(r-l)/(r+l+3).
moments=[[arb(0) for _ in range(2*N+1)] for _ in range(5)]
scale=arb(1)
for r in range(2*N+1):
    if r:scale*=c/r
    parity=r%2;m=arb(2)/(r+parity+1)
    indices=list(range(parity,5,2))
    temp={i:arb(0) for i in indices}
    for j,l in enumerate(range(parity,min(r,trial[parity]['degrees'][-1])+1,2)):
        for i in indices:temp[i]+=trial[i]['a'][j]*m
        m*=arb(r-l)/(r+l+3)
    for i in indices:moments[i][r]=scale*temp[i]
print('Gram and scaled moments constructed; bits=',ctx.prec,'N=',N,flush=True)
for i in range(5):
    for j in range(i,5):
        if i%2!=j%2:continue
        accum=arb(0);par=i%2
        for k in range(N+1):
            r=2*k
            conv=sum((moments[i][l]*moments[j][r-l] for l in range(par,r+1,2)),arb(0))
            accum+=(-1)**(k+par)*conv/(2*k+1)
        Q[i][j]=Q[j][i]=10*accum
        print('Q_N entry',i,j,str(Q[i][j]),flush=True)
L=[[G[i][j]-Q[i][j] for j in range(5)] for i in range(5)]
ratio=(2*c)**2/((2*N+4)*(2*N+5));assert ratio<1
kernel_operator_error=20*(2*c)**(2*N+2)/math.factorial(2*N+3)
assert kernel_operator_error<arb('1e-80')
print('uniform operator-norm Taylor error <=',kernel_operator_error,flush=True)
def maximum(v):
    result=v[0]
    for t in v[1:]:result=(result+t+abs(result-t))/2
    return result
def round_up(ball,digits=8):
    with localcontext() as context:
        context.prec=160
        b=Decimal(ball.upper().str(130,radius=False))
        step=Decimal(1).scaleb(b.adjusted()-digits+1)
        cap=b.quantize(step,rounding=ROUND_CEILING)
        while not arb(str(cap))>=ball:cap+=step
        assert arb(str(cap))>=ball
        return str(cap)
results=[]
for n in range(5):
    eta=maximum([sum((abs(G[i][j]-(1 if i==j else 0)) for j in range(n+1)),arb(0)) for i in range(n+1)])
    lowG=1-eta;assert lowG>0
    upperL=maximum([L[i][i]+sum((abs(L[i][j]) for j in range(n+1) if j!=i),arb(0)) for i in range(n+1)])
    assert upperL>0
    bound=upperL/lowG+kernel_operator_error
    Un=round_up(bound)
    assert arb(Un)>=bound
    diagnostic=arb(json.loads((out/'diagnostics.json').read_text())[n]['leak_K2'])
    assert arb(Un)/diagnostic<2
    result={'n':n,'gram_eta':str(eta),'gram_lower':str(lowG),'gershgorin_L_upper':str(upperL),'rayleigh_leak_upper':str(bound),'U_n':Un,'diagnostic_leak':str(diagnostic),'U_over_diagnostic':str(arb(Un)/diagnostic)}
    results.append(result);print('CERTIFIED',json.dumps(result),flush=True)
(out/'certificate.json').write_text(json.dumps({'bits':ctx.prec,'N':N,'degree_max':127,'kernel_operator_error':str(kernel_operator_error),'G':[[str(t) for t in row] for row in G],'Q_N':[[str(t) for t in row] for row in Q],'L_N':[[str(t) for t in row] for row in L],'results':results,'method':'full (n+1)-dimensional polynomial trial space, sinc Taylor remainder, Gershgorin and Gram lower bound'},indent=2)+'\n')
print('COMPLETE: all five rigorous starting upper bounds and safe-rounding assertions passed.',flush=True)
