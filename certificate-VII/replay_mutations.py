"""Replay the seven mutation tests of notes/0785_gpt_A4_mu12_certificate without touching the shipped files.

The original test script (whose hash is recorded in the provenance files) refuses to overwrite the shipped
mutation_results.json. This wrapper makes a disposable copy of this folder (hard links where possible, so it
costs no space), removes the shipped results from the copy only, runs the original script there, and checks
that every case is rejected for the same reason as in the shipped mutation_results.json.

    python replay_mutations.py            # from the certificate-VII folder; prints the copy's location
"""
import json,os,shutil,subprocess,sys,time
from pathlib import Path

HERE=Path(__file__).resolve().parent
NOTE=Path('notes/0785_gpt_A4_mu12_certificate')
SKIP={'.git','__pycache__','replay'}

def link_tree(src,dst):
    for root,dirs,files in os.walk(src):
        dirs[:]=[d for d in dirs if d not in SKIP]
        rel=Path(root).relative_to(src);(dst/rel).mkdir(parents=True,exist_ok=True)
        for f in files:
            try:os.link(Path(root)/f,dst/rel/f)
            except OSError:shutil.copy2(Path(root)/f,dst/rel/f)

def main():
    shipped=json.loads((HERE/NOTE/'mutation_results.json').read_text())
    assert not (HERE/NOTE/'mutations').exists(),'remove the mutations/ folder left by an earlier direct run first'
    copy=HERE/'replay'/time.strftime('mutations-%Y%m%d-%H%M%S');assert not copy.exists()
    link_tree(HERE,copy)
    (copy/NOTE/'mutation_results.json').unlink()   # removes the link in the copy only
    proc=subprocess.run([sys.executable,str(copy/NOTE/'test_mutations_v2.py')],cwd=copy)
    if proc.returncode!=0:sys.exit('FAIL: the mutation script exited with %d (copy kept at %s)'%(proc.returncode,copy))
    got=json.loads((copy/NOTE/'mutation_results.json').read_text())
    pairs=lambda r:[(t['case'],t['actual_reason']) for t in r['tests']]
    assert got['status']=='PASS' and got['base_packet_unchanged'] is True
    assert pairs(got)==pairs(shipped),(pairs(got),pairs(shipped))
    for case,reason in pairs(got):print(f'  {case}: rejected with {reason}')
    print('PASS: all 7 mutations rejected for the same reasons as in the shipped mutation_results.json')
    print('disposable copy:',copy.relative_to(HERE),'(may be deleted)')

if __name__=='__main__':main()
