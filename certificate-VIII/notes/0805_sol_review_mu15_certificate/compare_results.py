"""Compare the independently proved gates with pinned 1280/1536-bit records.
Numbers, differences and ratios are retained. Candidate records are not proof inputs.
"""
from pathlib import Path
from flint import arb,ctx
import json
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
BASE=ROOT/'notes/0803_gpt_A5_mu15_certificate'
ctx.prec=1664
own=json.loads((OUT/'independent_certificate.json').read_text())
spots=json.loads((OUT/'source_spotchecks.json').read_text())
assert own['status']=='PASS_INDEPENDENT_CERTIFICATE' and spots['status']=='PASS_INDEPENDENT_SOURCE_SPOTCHECKS'
assert own['reviewer']==spots['reviewer']=='sol'
refs={bits:json.loads((BASE/f'check{bits}_final.json').read_text()) for bits in (1280,1536)}
for bits,ref in refs.items():
    assert ref['status']=='PASS_EXTERNAL_CERTIFICATE' and ref['precision_bits']==bits
def ss(x):return x.str(40)
comparisons=[];precision_rows=[]
for p in (0,1):
    row=own['results'][p]
    for stage in ('C','W','T'):
        value=arb(row[stage]['lower']);assert value>0
        for bits,ref in refs.items():
            given=arb(ref['results'][p][stage]['positive_lower_bound'])
            diff=value-given;rel=abs(diff/given)
            assert rel<arb('1e-37')
            comparisons.append({'parity':p,'stage':stage,'reference_bits':bits,
                'independent_lower':ss(value),'reference_lower':ss(given),
                'difference':ss(diff),'ratio':ss(value/given),'relative_difference_upper':ss(rel.upper())})
        v1280=arb(refs[1280]['results'][p][stage]['positive_lower_bound'])
        v1536=arb(refs[1536]['results'][p][stage]['positive_lower_bound'])
        relative=abs((v1280-v1536)/v1536)
        assert relative<arb(2)**-1000
        precision_rows.append({'parity':p,'stage':stage,'lower1280':ss(v1280),'lower1536':ss(v1536),
            'difference':ss(v1280-v1536),'ratio':ss(v1280/v1536),'relative_difference_upper':ss(relative.upper())})
    for bits,ref in refs.items():
        eps=arb(row['epsilon_T']);given=arb(ref['results'][p]['epsilon_T'])
        assert abs((eps-given)/given)<arb('1e-37')
    for stage in ('C','W','T'):assert arb(row[stage]['paid_residual'])<1
    assert arb(row['epsilon_T'])<arb(2)**-300 and arb(row['epsilon_Y'])<arb(2)**-310

mutations=json.loads((BASE/'mutation_results.json').read_text())
expected={'sign_reversed_J':'E_SOURCE_J','excess_error':'E_ERROR_BUDGET',
 'B_replaced_by_A_squared':'E_SOURCE_B','file_tampering':'E_FILE_HASH',
 'understated_J_error':'E_SOURCE_J','wrong_contract':'E_CONTRACT',
 'negative_T_same_congruence_gate':'E_POSITIVITY_T'}
assert mutations['status']=='PASS' and len(mutations['tests'])==7 and mutations['base_packet_unchanged'] is True
for test in mutations['tests']:
    assert test['expected_reason']==test['actual_reason']==expected[test['case']] and test['passed'] is True
    if test['case']!='negative_T_same_congruence_gate':
        assert test['exit_code']==1
        result=json.loads((BASE/'mutations'/test['case']/'result.json').read_text())
        assert result['status']=='REJECT' and result['reason']==expected[test['case']]
assert own['negative_tests'][0]['passed'] is True
assert len(spots['records'])==10 and all(r['status']=='PASS' for r in spots['records'])
assert any(r['row']==400 and r['column']==401 for r in spots['records'])
result={'status':'PASS_COMPARISON_AND_NEGATIVE_RECORD_AUDIT','reviewer':'sol',
 'independent_vs_candidate':comparisons,'candidate_two_precisions':precision_rows,
 'candidate_negative_tests':'7 specified reasons checked, 6 result files read; negative T is an in-memory gate test',
 'independent_negative_test':own['negative_tests'],
 'certified_T_display_lower_bounds':['2.3979369226e-70','3.9458005898e-66']}
for p,bound in enumerate(result['certified_T_display_lower_bounds']):
    assert arb(own['results'][p]['T']['lower'])>arb(bound)
(OUT/'compare_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
for row in comparisons:
    if row['reference_bits']==1280:
        print('COMPARE',row['parity'],row['stage'],'own',row['independent_lower'],'reference',row['reference_lower'],'difference',row['difference'],'ratio',row['ratio'],flush=True)
for row in precision_rows:print('TWO PRECISIONS',row['parity'],row['stage'],'relative difference upper',row['relative_difference_upper'],flush=True)
print('COMPLETE: 12 independent/reference comparisons, 6 precision comparisons, 7 negative records, 1 own mutation, 10 source entries',flush=True)
