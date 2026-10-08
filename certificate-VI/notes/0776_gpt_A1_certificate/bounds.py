"""Exact rational scalars justified in ANALYTIC_BOUNDS.md, not numerical fits."""
from fractions import Fraction as F
from math import factorial

def analytic_bounds():
    # 10 panels have total length 256; Gauss q=192, polynomial degree 383.
    eb=F(320*1024*3*10**48,2**383)
    ek=F(320*1024*4000*10**48,2**383)
    dp=2*F(3,5)**721/factorial(721)
    ep=160*dp
    eu=81*(2*eb+eb*eb)
    assert dp<1 and ek+ep+eu<F(1,2**190)
    m=F(1,2);r=F(1,2**190);p=F(1,2**440)
    e=1500*r+p;n=896*r*r+p;h=12000*r*r+p
    assert n<m
    alpha=e+2*e*e/(m-n);beta=(e+n)/(m*m)+2*h*h/(m*m*(m-n))
    assert alpha+beta*F(13,2)**2<F(1,2**170)
    return dict(band=eb,frequency=ek,pole=ep,U=eu,J_quad=ek+ep+eu,
                alpha=alpha,beta=beta,n=n)

def fraction_text(q):
    return str(q.numerator)+'/'+str(q.denominator)
