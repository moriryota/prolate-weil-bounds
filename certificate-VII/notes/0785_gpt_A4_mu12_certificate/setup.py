from pathlib import Path
import json,hashlib,shutil,difflib
R=Path(__file__).resolve().parents[2];O=Path(__file__).parent;old=R/'notes/0776_gpt_A1_certificate';src=R/'notes/0784_gpt_A3_mu12';tag='n832p1920q256'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
pinned=json.loads((src/'OUTPUT_SHA256.json').read_text());manifest=src/(tag+'_source_manifest.json');assert sha(manifest)==pinned[str(manifest.relative_to(R))]
pins={'N':832,'q':256,'bits':1920,'original_delta':'d(416 L)','target_delta':'2^-1515','manifest':{'path':str(manifest.relative_to(R)),'sha256':sha(manifest)},'files':{}}
(O/'sources').mkdir(exist_ok=True)
for p in (0,1):
 for name in ('A','B','J'):
  source=src/f'{tag}_p{p}_{name}.json.gz';assert sha(source)==pinned[str(source.relative_to(R))]
  dest=O/'sources'/source.name;assert not dest.exists();shutil.copyfile(source,dest)
  pins['files'][f'p{p}_{name}']={'path':str(dest.relative_to(R)),'sha256':sha(dest),'original_path':str(source.relative_to(R))}
(O/'SOURCE_PINS.json').write_text(json.dumps(pins,indent=2))
for name in ('codec.py','propose.py','check_certificate.py','test_mutations_v2.py'):
 s=(old/name).read_text().replace('320','416').replace('640','832')
 if name=='propose.py':
  s=s.replace('ctx.prec=1280','ctx.prec=1280;ctx.threads=1').replace("2**185","2**250").replace("'A1-0776-v1'","'A4-0785-v1'").replace("'Omega':256","'Omega':416")
  s=s.replace("'errors':eps,'blocks':blocks", "'mu':12,'m':'1/2','gamma':'65/16','b':'61/8','u':[207,68],'delta':'2^-1515','tau':'2^-1518','errors':eps,'blocks':blocks")
 if name=='check_certificate.py':
  s=s.replace('def verify(packet):','def verify(packet,bits=1280):').replace('ctx.prec=1280;start=', 'ctx.prec=bits;ctx.threads=1;start=')
  s=s.replace("cfg['schema']=='A1-0776-v1' and cfg['N']==832 and cfg['Omega']==256", "cfg['schema']=='A4-0785-v1' and cfg['N']==832 and cfg['Omega']==416 and cfg.get('mu')==12 and cfg.get('m')=='1/2' and cfg.get('gamma')=='65/16' and cfg.get('b')=='61/8' and cfg.get('u')==[207,68] and cfg.get('delta')=='2^-1515' and cfg.get('tau')=='2^-1518'")
  s=s.replace("(832,1280,192)","(832,1920,256)").replace("(arb(13)/2)*(arb(13)/2)","(arb(61)/8)*(arb(61)/8)")
  s=s.replace("arb('1e-46')","arb(2)**-240").replace('fixed 1e-46 budget','fixed 2^-240 budget')
  s=s.replace("'tau':'2^-49162'","'tau':'2^-1518','delta':'2^-1515','precision_bits':bits,'scalar_checks':scalar_report")
  s=s.replace('from bounds import analytic_bounds','from bounds import analytic_bounds, validate_scalars')
  s=s.replace("constants=analytic_bounds();alpha", "scalar_report=validate_scalars();constants=analytic_bounds();alpha")
  s=s.replace("parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()", "parser.add_argument('--output',type=Path,required=True);parser.add_argument('--bits',type=int,choices=[1024,1280,1536,1920],default=1280);args=parser.parse_args()")
  s=s.replace('verify(args.packet);code=0','verify(args.packet,args.bits);code=0')
 if name=='test_mutations_v2.py':
  s=s.replace("mutations_v2","mutations").replace("mutation_results_v2.json","mutation_results.json")
  s=s.replace("    # A real positivity", "    p,c=fixture('understated_J_error');c['errors']['J']='1/'+str(2**300)\n    results.append(run(p,c,'E_SOURCE_J'))\n    p,c=fixture('wrong_contract');c['tau']='2^-1517'\n    results.append(run(p,c,'E_CONTRACT'))\n    # A real positivity")
  s=s.replace('All 5 mutations','All 7 mutations')
 (O/name).write_text(s)
 (O/(name+'.diff')).write_text(''.join(difflib.unified_diff((old/name).read_text().splitlines(True),s.splitlines(True),fromfile='0776/'+name,tofile='0785/'+name)))
(O/'.gitignore').write_text('*.json.gz\n__pycache__/\n')
(O/'REUSE_PINS.json').write_text(json.dumps({str(p.relative_to(R)):sha(p) for p in [*(old/n for n in ['codec.py','propose.py','check_certificate.py','test_mutations_v2.py']),R/'notes/0778_gpt_A2_beyond_L76/scalar_gate.py',R/'python/weil_spectral.py']},indent=2))
print('copied six pinned sources and four programs')
