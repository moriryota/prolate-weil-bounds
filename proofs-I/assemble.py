"""Arb evaluation of the prescribed Theorem D assembly with certified U_n."""
from pathlib import Path
from decimal import Decimal, ROUND_CEILING, localcontext
from flint import arb, ctx
import sympy as sp
import math, json

out = Path(__file__).resolve().parent
ctx.prec = 512
pi = arb.pi()
c = 10*pi
u = sp.symbols('u')
certificate = json.loads((out/'certificate.json').read_text())

def round_up(ball, digits=8):
    with localcontext() as context:
        context.prec = 160
        endpoint = Decimal(ball.upper().str(130, radius=False))
        step = Decimal(1).scaleb(endpoint.adjusted()-digits+1)
        cap = endpoint.quantize(step, rounding=ROUND_CEILING)
        while not arb(str(cap)) >= ball:
            cap += step
        assert arb(str(cap)) >= ball
        return str(cap)

def tail(j):
    poly = sp.Poly(sp.expand(sp.hermite(j, u)**2), u)
    X = c.sqrt()
    assert c >= poly.degree()
    total = sum((abs(int(cf))*(-c).exp()*X**(m-1)
                 for (m,), cf in poly.terms()), arb(0))
    return 2*total/(2**j*math.factorial(j)*pi.sqrt())

results = []
for row in certificate['results']:
    n = row['n']
    start = arb(row['U_n'])
    assert start >= arb(row['rayleigh_leak_upper'])
    tau = sum((tail(j) for j in range(n%2, n+1, 2)), arb(0))
    a = (2*n+1)/(c*(1-tau))
    assert 0 < a < 1 and 0 < tau < 1
    beta = 30+(2*n+1)*c*tau/(2*(1-tau)) + arb(9)*(2*n+1)**2/(32*(1-tau)**2*(1-a))
    original_beta = c**2*(a/2-arb(2*n+1)/(2*c)+arb(9)*a**2/(32*(1-a)))+30
    # Algebraic identity; Arb overlap is a computation check, not its proof.
    assert (beta-original_beta).contains(0)
    factor = (2*c).exp()*c**(-(arb(n)+arb(1)/2))*(beta/c).exp()
    C = start*factor
    cap = round_up(C)
    assert arb(cap) >= C
    oldC = arb('7e-9')*factor
    comparison = 4*pi.sqrt()*8**n/math.factorial(n)
    result = {'n':n, 'U_n':row['U_n'], 'tau':str(tau), 'a':str(a), 'beta':str(beta),
              'C_n_ball':str(C), 'C_n_safe_upper':cap,
              'comparison_expression':str(comparison),
              'C_upper_over_comparison':str(arb(cap)/comparison),
              'old_C_n':str(oldC), 'old_start_over_U':str(arb('7e-9')/start)}
    results.append(result)
    print(json.dumps(result), flush=True)
(out/'constants.json').write_text(json.dumps({'bits':ctx.prec, 'c0':'10*pi',
    'assembly':'python/theoremD_all_n.py, start replaced by U_n',
    'comparison_status':'evaluation of given 4*sqrt(pi)*8^n/n!, not a literature verification',
    'results':results}, indent=2)+'\n')
print('COMPLETE: all five C_n upper-endpoint rounding assertions passed.', flush=True)
