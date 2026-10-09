from pathlib import Path
import difflib,json,hashlib
r=Path.cwd();o=r/'notes/0784_gpt_A3_mu12';old=r/'notes/0778_gpt_A2_beyond_L76'
s=(old/'diagnostic_mu12.py').read_text()
s=s.replace('0775: N640 source generator adapted from 0773','0784: mu12 source generator copied from 0778')
s=s.replace('assert N in (128,256,384,640) and bits in (768,1024,1280) and q in (160,192,224)','assert N in (640,832,1008) and 1280<=bits<=2560 and q in (224,256)')
s=s.replace('ctx.prec=bits;mp.mp.prec=bits','ctx.prec=bits;ctx.threads=1;mp.mp.prec=bits')
s=s.replace('assert all(u[j]*ehi[j]<Ctail for j in (0,1))','assert u==[207,68]\nassert all(u[j]*ehi[j]<Ctail for j in (0,1))')
s=s.replace("assert check['row_radius_below_1e_minus48'], (parity,name,check)","assert arb(check['max_radius_row_sum'])<arb('1e-70'), (parity,name,check)")
(o/'generate.py').write_text(s)
(o/'generator.diff').write_text(''.join(difflib.unified_diff((old/'diagnostic_mu12.py').read_text().splitlines(True),s.splitlines(True),fromfile='0778/diagnostic_mu12.py',tofile='0784/generate.py')))
base=(r/'notes/0775_gpt_A0_N640/diagnostic640_v2.py').read_text()
base=base[:base.index('def parameters():')]+base[base.index('def midpoint_ldl'):base.index('def diagnose(')]
(o/'matrix_utils.py').write_text(base)
pins={str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [old/'diagnostic_mu12.py',old/'diagnostic_support.py',r/'notes/0775_gpt_A0_N640/diagnostic640_v2.py',r/'python/weil_spectral.py']}
(o/'REUSE_PINS.json').write_text(json.dumps(pins,indent=2))
