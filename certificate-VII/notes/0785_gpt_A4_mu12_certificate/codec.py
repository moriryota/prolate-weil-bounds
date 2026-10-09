"""Only exact dyadic I/O. No positivity decision or source-error mathematics."""
import gzip,hashlib,json
from pathlib import Path
from flint import arb,arb_mat,ctx

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def dump_dyadic(path,M,bits=512):
    path=Path(path);assert not path.exists()
    scale=arb(2)**bits;rows=[]
    for i in range(M.nrows()):
        row=[]
        for j in range(M.ncols()):
            assert M[i,j].is_finite()
            z=(M[i,j].mid()*scale).floor().unique_fmpz()
            assert z is not None
            row.append(str(z))
        rows.append(row)
    raw=json.dumps({'denominator_power':bits,'rows':rows},separators=(',',':')).encode()
    path.write_bytes(gzip.compress(raw,mtime=0))
    return {'file':path.name,'sha256':sha(path)}

def read_dyadic(path):
    data=json.loads(gzip.decompress(Path(path).read_bytes()))
    if type(data.get('denominator_power')) is not int or not 0<=data['denominator_power']<=2048:
        raise ValueError('invalid dyadic exponent')
    rows=data['rows']
    if len(rows)!=416 or any(len(r)!=416 for r in rows):
        raise ValueError('matrix shape must be 416 x 416')
    den=arb(2)**data['denominator_power']
    ints=[[int(s) for s in row] for row in rows]
    if any(abs(v).bit_length()>2048 for row in ints for v in row):
        raise ValueError('oversized integer')
    if ctx.prec<2048:
        # All admitted dyadics in this packet are exact at 1280 bits.
        if any(abs(v).bit_length()>ctx.prec for row in ints for v in row):
            raise ValueError('precision insufficient for exact dyadic input')
    return arb_mat([[arb(v)/den for v in row] for row in ints])

def read_source(root,pins,parity,name):
    ref=pins['files'][f'p{parity}_{name}'];path=Path(root)/ref['path']
    if sha(path)!=ref['sha256']:raise ValueError('pinned source hash mismatch')
    return arb_mat([[arb(s) for s in row] for row in json.loads(gzip.decompress(path.read_bytes()))])
