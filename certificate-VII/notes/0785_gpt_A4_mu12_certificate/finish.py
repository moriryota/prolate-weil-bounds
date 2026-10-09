from pathlib import Path
from flint import arb,ctx
import json,hashlib,sys,platform,importlib.metadata as md,difflib
ctx.prec=512;O=Path(__file__).resolve().parent;R=O.parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
a=json.loads((O/'verification.json').read_text());b=json.loads((O/'verification_1536.json').read_text());mut=json.loads((O/'mutation_results.json').read_text())
assert a['status']==b['status']=='PASS_EXTERNAL_CERTIFICATE' and a['tau']==b['tau']=='2^-1518'
assert mut['status']=='PASS' and len(mut['tests'])==7 and all(r['passed'] for r in mut['tests'])
comparisons=[]
for ra,rb in zip(a['results'],b['results']):
 lower=arb(ra['T']['positive_lower_bound']);lower2=arb(rb['T']['positive_lower_bound']);k=180 if ra['parity']==0 else 166
 assert lower>arb(2)**-k and lower2>arb(2)**-k
 rel=abs(lower-lower2)/lower;assert rel<arb('1e-80')
 assert arb(ra['epsilon_T'])<arb(2)**-249 and arb(rb['epsilon_T'])<arb(2)**-249
 comparisons.append({'parity':ra['parity'],'relative_lower_bound_difference':str(rel),'certified_dyadic_lower':'2^-'+str(k),'epsilon_T_less_than':'2^-249'})
(O/'PRECISION_CHECK.json').write_text(json.dumps({'status':'PASS','runs':[1280,1536],'results':comparisons},indent=2))
env={'python':sys.version,'platform':platform.platform(),'packages':{k:md.version(k) for k in ('python-flint','mpmath','numpy')},'VM_used':False}
(O/'ENVIRONMENT.json').write_text(json.dumps(env,indent=2))
old=R/'notes/0776_gpt_A1_certificate'
for name in ('codec.py','propose.py','check_certificate.py','test_mutations_v2.py'):
 (O/(name+'.diff')).write_text(''.join(difflib.unified_diff((old/name).read_text().splitlines(True),(O/name).read_text().splitlines(True),fromfile='0776/'+name,tofile='0785/'+name)))
fmt=lambda s:arb(s).str(14)
lines=['# 0785 GPT：μ12の完全な外部区間証明書','',
'2026-10-09。依頼0785を完了。**通常の解析論証と外部区間証明書による候補。追加成果の独立査読・Lean化は未実施。**','',
'## 結論','',
'固定窓 I=[−log(12)/2, log(12)/2] に台を持つ複素C²関数f（f=0も含む）について、','',
r'$$Q(f)\ge\tau\|f\|_2^2,\qquad\tau=2^{-1518}>0.$$','',
'を支える有理行列証明書を作成し、提案器をimportしない検査器が全誤差を支払って受理した。1280bitと1536bitで同じパケットを検査し、ともにPASS_EXTERNAL_CERTIFICATE。これはLeanや第三者の承認ではない。','',
'## 認証結果','',
'|パリティ|全支払い後の有限条件の下界（表示は丸め）|簡単な厳密下界|εT|','|---|---|---|---|']
for z,c in zip(a['results'],comparisons):
 lines.append(f"|{['偶','奇'][z['parity']]}|{fmt(z['T']['positive_lower_bound'])}|>{c['certified_dyadic_lower']}|{fmt(z['epsilon_T'])}|")
lines+=['','両方でεT<2^-249。最終のη+εT||Z||F²<1に加え、C_M=B−mAとWの正定値性、Xの源誤差、D=B−A²の誤差、Yの逆残差も検査済み。','',
'## 0784から埋めた残件','',
'- δを2^-1515に固定。d(416L)は2^-1515以上2^-1514未満、C>1/8なのでτ=δ/8。',
'- 0784の1920bit/q256の源6行列をSHA照合して複製。旧δ=d(c)から新δへの変更を2^-1504としてJに加算した。',
'- μ12の素数Schur評価、帯域・極・Uの求積上界、射影誤差を検査器が毎回再評価。μ10の係数の名前だけの置換はしていない。',
'- q256のJ求積上界2^-263、δ変更、丸めをまとめてεJ=2^-250で包含。有理A0/B0/J0/X/W0/Y/Zを2^512分母で固定。',
'- 有限条件から全Hの有界比較作用素へ、そこから試験関数域のWeil形式へつなぐ論証をTHEOREM.mdに記した。','',
'## 負例試験','',
'|変異|所定の拒否理由|','|---|---|']
for z in mut['tests']:lines.append(f"|{z['case']}|{z['actual_reason']}|")
lines+=['','CLIの負例は終了コード1、status=REJECT、期待reasonの完全一致を必要とした。予期しない例外・基盤障害のERROR（終了コード2）を合格扱いしない。負のTの合同検査は、同じ数学的判定関数に対する直接試験として区別して記録した。元パケットのSHA不変も検査した。','',
'## 資源と再現','',
'|工程|実測時間|ピークRSS|','|---|---:|---:|']
prop=json.loads((O/'proposal_result.json').read_text())
for title,r in [('有理候補の提案',prop),('1280bitの検査',a),('1536bitの検査',b)]:
 lines.append(f"|{title}|{r['seconds']:.2f}秒|{r['maxrss_bytes_macos']/1e9:.3f}GB|")
g=json.loads((O/'mutations_guard.json').read_text());lines.append(f"|7種の負例（監視器計測）|{g['seconds']:.2f}秒|{g['sampled_tree_peak_RSS_bytes']/1e9:.3f}GB（標本化したプロセスツリー）|")
lines+=['','重い実行は順次、監視停止8.5GB、VM不使用。源の再生成は不要だった。検査だけの再現コマンド：','',
'```sh',
'.venv/bin/python notes/0785_gpt_A4_mu12_certificate/run_guarded.py replay1280 notes/0785_gpt_A4_mu12_certificate/check_certificate.py --bits 1280 --output notes/0785_gpt_A4_mu12_certificate/replay1280.json',
'```','',
'必要データ・版・1536bit再現はREADME.md。既存の出力は上書きしない。大きい源・候補行列はGit対象外で、別途SHAとともに渡す。','',
'## 信頼の範囲と次の査読','',
'FLINT/Arb、Python、源生成コード、保存ball、古典的解析入力と通常の解析証明が信頼基盤に残る。検査器は提案器から独立した判定経路だが、別チームの実装ではなく、bounds/codecを共有する。SHAは同一性の証拠であり数学的正しさの代わりではない。',
'0783は再開後に完了しており、欠損補題とμ10系は0782・0783で受理済み。0784報告の停止中という記述を今回の時点では更新する。GENERAL_TAILは0782の1回の査読。今回のμ12の新しい接続の別担当査読はこれから。',
'読む順番：CONTRACT → THEOREM → ANALYTIC_BOUNDS → CHANGES_FROM_0776 → bounds/check/codec → verification・負例。特にδ摂動、ρ=2の複素楕円でのdigamma上界、ρ=3の射影誤差、GからWeil形式への接続を査読してほしい。',
'幅log12の固定窓であり、任意幅の正値性やRHを結論しない。文献上の優先権・新規性は今回調査していない。','',
'成果物：notes/0785_gpt_A4_mu12_certificate/。既存の0770〜0784、原稿、台帳、Gitは変更していない。入力確認と最終SHAはinputs_start/end.json、DELIVERY_CHECKS.json、OUTPUT_SHA256.jsonを参照。','']
(R/'reports/20261009_0785_GPT_A4.md').write_text('\n'.join(lines))
print('precision comparison and report complete')
