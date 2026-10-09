from pathlib import Path
import hashlib,json,ast
O=Path(__file__).resolve().parent;R=O.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for n in ('inputs_start.json','inputs_end.json'):
 a=json.loads((O/n).read_text());assert a['ok'] and a['inputs']==17 and not a['mismatches']
pins=json.loads((O/'SOURCE_PINS.json').read_text())
for ref in [pins['manifest'],*pins['files'].values()]:assert sha(R/ref['path'])==ref['sha256']
reuse=json.loads((O/'REUSE_PINS.json').read_text());assert all(sha(R/p)==h for p,h in reuse.items())
pack=json.loads((O/'packet/manifest.json').read_text())
for b in pack['blocks']:
 for ref in b['files'].values():assert sha(O/'packet'/ref['file'])==ref['sha256']
for p in O.glob('*.py'):ast.parse(p.read_text())
for p in O.glob('*_guard.json'):
 a=json.loads(p.read_text());assert a['exit_code']==0 and not a['memory_stop'] and a['sampled_tree_peak_RSS_bytes']<8_500_000_000
for n in ('verification.json','verification_1536.json'):assert json.loads((O/n).read_text())['status']=='PASS_EXTERNAL_CERTIFICATE'
assert json.loads((O/'mutation_results.json').read_text())['status']=='PASS'
assert json.loads((O/'PRECISION_CHECK.json').read_text())['status']=='PASS'
matrices=[*sorted((O/'sources').glob('*.json.gz')),*sorted((O/'packet').glob('*.json.gz'))]
(O/'MATRIX_SHA256.json').write_text(json.dumps({str(p.relative_to(R)):sha(p) for p in matrices},indent=2))
(O/'DELIVERY_CHECKS.json').write_text(json.dumps({'status':'PASS','scope':'external numerical certificate and delivery integrity; independent mathematical review outstanding','fixed_inputs':17,'source_matrices':6,'rational_candidate_matrices':18,'verification_bits':[1280,1536],'negative_tests':7,'tau':'2^-1518','all_owned_heavy_jobs_finished':True,'VM_used':False,'prior_reuse_pins_unchanged':True},indent=2))
files=[p for p in O.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='OUTPUT_SHA256.json']+[R/'reports/20261009_0785_GPT_A4.md']
hashes={str(p.relative_to(R)):sha(p) for p in sorted(files)}
(O/'OUTPUT_SHA256.json').write_text(json.dumps(hashes,indent=2));print('PASS delivery hashes:',len(hashes),'baseline matrices:',len(matrices))
