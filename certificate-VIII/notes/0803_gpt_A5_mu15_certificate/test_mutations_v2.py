"""Negative tests succeed only for the specified mathematical/integrity reason."""
import json,os,subprocess,sys,copy,time
from pathlib import Path
from flint import arb,ctx
from codec import read_dyadic,dump_dyadic,sha
import check_certificate as check
from bounds import analytic_bounds

OUT=Path(__file__).parent;BASE=OUT/'packet';ctx.prec=1280
def fixture(name):
    p=OUT/'mutations'/name;p.mkdir(parents=True,exist_ok=False)
    cfg=json.loads((BASE/'manifest.json').read_text())
    for b in cfg['blocks']:
        for ref in b['files'].values():os.link(BASE/ref['file'],p/ref['file'])
    return p,cfg
def replace(p,cfg,name,M):
    old=cfg['blocks'][0]['files'][name]
    path=p/old['file'];path.unlink() # Never write through the hard link to BASE.
    cfg['blocks'][0]['files'][name]=dump_dyadic(path,M)
def run(p,cfg,expected):
    (p/'manifest.json').write_text(json.dumps(cfg,indent=2))
    dest=p/'result.json'
    proc=subprocess.run([sys.executable,str(OUT/'check_certificate.py'),'--packet',str(p),'--output',str(dest)],capture_output=True,text=True)
    (p/'stdout.log').write_text(proc.stdout);(p/'stderr.log').write_text(proc.stderr)
    got=json.loads(dest.read_text()) if dest.exists() else {'status':'MISSING'}
    if not(proc.returncode==1 and got.get('status')=='REJECT' and got.get('reason')==expected):
        raise AssertionError((p.name,expected,proc.returncode,got))
    return {'case':p.name,'expected_reason':expected,'actual_reason':got['reason'],'exit_code':proc.returncode,'passed':True}
def main():
    start=time.monotonic();results=[]
    base_manifest=sha(BASE/'manifest.json')
    p,c=fixture('sign_reversed_J');replace(p,c,'J',-read_dyadic(BASE/c['blocks'][0]['files']['J']['file']))
    results.append(run(p,c,'E_SOURCE_J'))
    p,c=fixture('excess_error');c['errors']['J']='1/1'
    results.append(run(p,c,'E_ERROR_BUDGET'))
    p,c=fixture('B_replaced_by_A_squared');A=read_dyadic(BASE/c['blocks'][0]['files']['A']['file'])
    replace(p,c,'B',A*A);results.append(run(p,c,'E_SOURCE_B'))
    p,c=fixture('file_tampering')
    target=p/c['blocks'][0]['files']['A']['file'];raw=target.read_bytes();target.unlink();target.write_bytes(raw+b'x')
    results.append(run(p,c,'E_FILE_HASH'))
    p,c=fixture('understated_J_error');c['errors']['J']='1/'+str(2**400)
    results.append(run(p,c,'E_SOURCE_J'))
    p,c=fixture('wrong_contract');c['tau']='2^-2903'
    results.append(run(p,c,'E_CONTRACT'))
    # A real positivity rejection, in addition to the end-to-end provenance tests.
    c=json.loads((BASE/'manifest.json').read_text());refs=c['blocks'][0]['files']
    mats={k:read_dyadic(BASE/refs[k]['file']) for k in ['A','B','J','Y','ZT']}
    scalars=analytic_bounds();I=check.identity()
    T=mats['Y']+mats['J']-check.ab(scalars['alpha'])*I-check.ab(scalars['beta'])*(mats['B']-mats['A']*mats['A'])
    try:check.congruence(-T,mats['ZT'],arb(0),'T')
    except check.Rejected as e:
        assert e.code=='E_POSITIVITY_T'
        results.append({'case':'negative_T_same_congruence_gate','expected_reason':e.code,'actual_reason':e.code,'passed':True})
    else:raise AssertionError('negative T was accepted')
    assert sha(BASE/'manifest.json')==base_manifest
    for b in c['blocks']:
        for ref in b['files'].values():assert sha(BASE/ref['file'])==ref['sha256']
    out=OUT/'mutation_results.json';assert not out.exists()
    out.write_text(json.dumps({'status':'PASS','tests':results,'base_packet_unchanged':True,'seconds':time.monotonic()-start},indent=2))
    print('All 7 mutations rejected for their specified reasons')
if __name__=='__main__':main()
