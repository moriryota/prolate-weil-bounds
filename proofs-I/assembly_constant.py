from flint import arb,ctx
from pathlib import Path
import json
ctx.prec=256
p=arb.pi();c0=10*p;d=arb('1.370');h=arb('.162');m=arb('1.837');r=arb('1.2');a=arb('.7');z=d/c0
parts={'main':2*r,'cross1':8*h*d/a,'cross2':p*h*(2*m)*(1+2*r/c0**2),'quadratic':2*d*d*((2*z).exp()+z.exp())}
total=sum(parts.values(),arb(0))
for k,v in parts.items():print(k,v)
print('sum',total,'safe b=15')
assert total<15
Path(__file__).with_suffix('.json').write_text(json.dumps({'parts':{k:str(v) for k,v in parts.items()},'sum':str(total),'safe_b':15,'condition':'Volterra identities and regular solution (paper I, Sections 3-4)'},indent=2)+'\n')
