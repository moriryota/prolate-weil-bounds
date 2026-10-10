"""Independent final checker: no import of propose.py, no eigensolver.

Trust: pinned Arb source evaluations + the analytic lemmas in ANALYTIC_BOUNDS.md,
Python integer arithmetic and Arb/FLINT interval arithmetic. This is not Lean.
"""
import argparse,json,time,sys,resource
from pathlib import Path
from fractions import Fraction
from flint import arb,arb_mat,ctx
from codec import read_dyadic,read_source,sha
from bounds import analytic_bounds, validate_scalars

OUT=Path(__file__).parent;ROOT=OUT.parents[1]
class Rejected(Exception):
    def __init__(self,code,detail):
        self.code=code;self.detail=detail;super().__init__(code+': '+detail)
def need(ok,code,detail):
    if not ok:raise Rejected(code,detail)
def ab(q):return arb(q.numerator)/q.denominator
def identity():
    I=arb_mat(832,832)
    for i in range(832):I[i,i]=1
    return I
def norm_bound(M):
    """sqrt(||M||_1 ||M||_infty), also valid for non-Hermitian residuals."""
    rows=[sum(abs(M[i,j]).upper() for j in range(M.ncols())).upper() for i in range(M.nrows())]
    cols=[sum(abs(M[i,j]).upper() for i in range(M.nrows())).upper() for j in range(M.ncols())]
    v=(max(rows)*max(cols)).sqrt().upper()
    need(v.is_finite(),'E_NONFINITE','matrix norm')
    return v
def frob_squared(Z):
    return sum(abs(Z[i,j]).upper()*abs(Z[i,j]).upper()
               for i in range(Z.nrows()) for j in range(Z.ncols())).upper()
def symmetric(M,name):
    need(all(M[i,j]==M[j,i] for i in range(832) for j in range(i)),
         'E_SYMMETRY',name+' must be exactly symmetric')
def congruence(M,Z,perturbation,stage):
    I=identity()
    eta=norm_bound(Z.transpose()*M*Z-I)
    s=frob_squared(Z);need(s>0,'E_WITNESS',stage)
    total=eta+perturbation*s
    need(total<1,'E_POSITIVITY_'+stage,'congruence residual plus source error is not below 1')
    return {'residual':str(eta),'frob_squared_upper':str(s),'paid_residual':str(total),
            'positive_lower_bound':str(((1-total)/s).lower())}

def verify(packet,bits=1280):
    ctx.prec=bits;ctx.threads=1;start=time.monotonic();I=identity()
    cfg=json.loads((packet/'manifest.json').read_text())
    need(cfg['schema']=='A5-0803-v1' and cfg['N']==1664 and cfg['Omega']==736 and cfg.get('mu')==15 and cfg.get('m')=='1/2' and cfg.get('gamma')=='37/8' and cfg.get('b')=='35/4' and cfg.get('u')==[428,142] and cfg.get('delta')=='2^-2901' and cfg.get('tau')=='2^-2904','E_CONTRACT','fixed N and Omega')
    need(cfg['source_pins_sha256']==sha(OUT/'SOURCE_PINS.json'),'E_SOURCE_PINS','pin file changed')
    pins=json.loads((OUT/'SOURCE_PINS.json').read_text())
    for pinfile in ['SOURCE_CODE_PINS.json','REVIEW_0802_PINS.json']:
        for name,digest in json.loads((OUT/pinfile).read_text()).items():
            need(sha(ROOT/name)==digest,'E_DEPENDENCY_HASH',name)
    need((pins['N'],pins['bits'],pins['q'])==(1664,3328,224),'E_SOURCE_PINS','source configuration')
    for ref in pins['manifests']:
        need(sha(ROOT/ref['path'])==ref['sha256'],'E_SOURCE_HASH','source manifest')
        metadata=json.loads((ROOT/ref['path']).read_text())
        need((metadata['mu'],metadata['Omega'],metadata['N'],metadata['bits'],metadata['q'])==(15,736,1664,3328,224),'E_SOURCE_PINS','source metadata')
        need(metadata['parity'] in [0,1] and metadata['u']==[428,142] and arb(metadata['gamma'])==arb(37)/8,'E_SOURCE_PINS','source parity / constants')
        ends=[arb(0),arb(1)/2]+[arb(2)**k for k in range(7)]+[arb(k) for k in range(128,736,64)]+[arb(736)]
        need(len(metadata['panel_ends'])==len(ends) and all(arb(x)==y for x,y in zip(metadata['panel_ends'],ends)),'E_SOURCE_PINS','frequency panels')
    errors={k:ab(Fraction(v)) for k,v in cfg['errors'].items()}
    need(set(errors)=={'A','B','J','W','D'} and all(x>0 for x in errors.values()),'E_ERRORS','positive error bounds required')
    eA,eB,eJ,eW,eD=[errors[k] for k in ['A','B','J','W','D']]
    scalar_report=validate_scalars();constants=analytic_bounds();alpha=ab(constants['alpha']);beta=ab(constants['beta'])
    need([b['parity'] for b in cfg['blocks']]==[0,1],'E_CONTRACT','both parities, once each')
    results=[]
    print('scalar gates passed',flush=True)
    for block in cfg['blocks']:
        p=block['parity'];mat={}
        print('loading packet parity',p,flush=True)
        need(set(block['files'])=={'A','B','J','X','W0','Y','ZC','ZW','ZT'},'E_FORMAT','matrix names')
        for name,ref in block['files'].items():
            need(Path(ref['file']).name==ref['file'],'E_PATH','local basenames only')
            path=packet/ref['file']
            need(path.is_file() and sha(path)==ref['sha256'],'E_FILE_HASH',str(p)+' '+name)
            mat[name]=read_dyadic(path)
        A,B,J,X,W0,Y,ZC,ZW,ZT=[mat[k] for k in ['A','B','J','X','W0','Y','ZC','ZW','ZT']]
        for name in ['A','B','J','W0','Y']:symmetric(mat[name],name)
        source_errors={}
        for name,claimed in [('A',eA),('B',eB),('J',eJ)]:
            print('checking source',p,name,flush=True)
            S=read_source(ROOT,pins,p,name)
            measured=norm_bound(mat[name]-S)
            total=measured+(ab(constants['J_quad']) if name=='J' else 0)
            need(total<=claimed,'E_SOURCE_'+name,'candidate does not enclose pinned true source with the declared bound')
            source_errors[name]={'rounding_and_center':str(measured),'analytic_quadrature':str(ab(constants['J_quad']) if name=='J' else arb(0))}
        # C>0 for true A,B, not just the rational centre.
        print('checking C W inverse T',p,flush=True)
        C=B-A/2
        ctest=congruence(C,ZC,eB+eA/2,'C')
        xn=norm_bound(X);an=norm_bound(A)
        Wexpr=X+X.transpose()-X.transpose()*A*X+2*(I-A*X-X.transpose()*A+X.transpose()*B*X)
        wround=norm_bound(Wexpr-W0)
        werror=(xn*xn+4*xn)*eA+2*xn*xn*eB+wround
        need(werror<=eW,'E_W_ENCLOSURE','rational W0 does not contain true W(X)')
        derror=eB+2*an*eA+eA*eA
        need(derror<=eD,'E_D_ENCLOSURE','B-A^2 error too small')
        wtest=congruence(W0,ZW,eW,'W')
        # Wup dominates true W(X), which dominates G; Wup is positive.
        Wup=W0+eW*I
        residual=norm_bound(I-Wup*Y)
        need(residual<1,'E_INVERSE_RESIDUAL','Neumann residual not below 1')
        eY=norm_bound(Y)*residual/(1-residual)
        need(eY<=arb(2)**-310,'E_INVERSE_ERROR','inverse error exceeds 0802 cap')
        inverse_source_error=2*(arb(35)/4)*(arb(35)/4)*eW
        total_error=inverse_source_error+eY+eJ+beta*eD
        need(total_error<=arb(2)**-300,'E_ERROR_BUDGET','total T error exceeds fixed 2^-300 budget')
        T=Y+J-alpha*I-beta*(B-A*A)
        ttest=congruence(T,ZT,total_error,'T')
        results.append({'parity':p,'source_errors':source_errors,'C':ctest,'W':wtest,'T':ttest,
                        'epsilon_W_required':str(werror),'epsilon_D_required':str(derror),
                        'inverse_residual':str(residual),'epsilon_Y':str(eY),
                        'epsilon_inverse_source':str(inverse_source_error),'epsilon_T':str(total_error)})
        print('parity accepted',p,flush=True)
        rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        need(rss<8_000_000_000,'E_MEMORY','8GB checkpoint exceeded')
    return {'status':'PASS_EXTERNAL_CERTIFICATE','scope':'Numerical certificate plus documented analytic source bounds; not Lean or independent human review',
            'tau':'2^-2904','delta':'2^-2901','precision_bits':bits,'scalar_checks':scalar_report,'results':results,'seconds':time.monotonic()-start,
            'maxrss_bytes_macos':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--packet',type=Path,default=OUT/'packet')
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--bits',type=int,choices=[1024,1280,1536,1920],default=1280);args=parser.parse_args()
    if args.output.exists():raise RuntimeError('refuse to overwrite an earlier result')
    try:result=verify(args.packet,args.bits);code=0
    except Rejected as e:result={'status':'REJECT','reason':e.code,'detail':e.detail};code=1
    except Exception as e:result={'status':'ERROR','reason':'E_INFRASTRUCTURE','detail':repr(e)};code=2
    args.output.write_text(json.dumps(result,indent=2));print(json.dumps(result));sys.exit(code)
if __name__=='__main__':main()
