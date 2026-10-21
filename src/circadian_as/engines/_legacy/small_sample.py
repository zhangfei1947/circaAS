"""BB-Pearson-df: null residual-df dispersion + common-rho nested mean test.

An empirical estimating-equation route, not Cox-Reid/REML. The null BB mean
profile is used in the Pearson equation Q(rho)=m-rank(Xnull). Bootstrap refits
all nuisance estimation. No coverage gate is changed or applied here.
"""
from __future__ import annotations
from functools import lru_cache
import numpy as np
from scipy.optimize import brentq, minimize
from scipy.special import expit,logit,gammaln,betaln,digamma
from scipy.stats import chi2

RHO_MAX=.98
COEFFICIENT_MAX=18.


def design(t,block=None,harmonic=True,levels=None):
    t=np.asarray(t,float);cols=[np.ones(len(t))]
    if block is not None:
        block=np.asarray(block).astype(str);levels=sorted(set(block)) if levels is None else levels
        if not set(block).issubset(levels):raise ValueError('Unknown prediction experiment')
        cols.extend((block==lev).astype(float) for lev in levels[1:])
    if harmonic:cols.extend([np.cos(2*np.pi*t/24),np.sin(2*np.pi*t/24)])
    return np.column_stack(cols)


def _logpdf(y,n,p,rho):
    p=np.clip(p,1e-10,1-1e-10);const=gammaln(n+1)-gammaln(y+1)-gammaln(n-y+1)
    if rho==0:return const+y*np.log(p)+(n-y)*np.log1p(-p)
    k=(1-rho)/rho
    return const+betaln(y+p*k,n-y+(1-p)*k)-betaln(p*k,(1-p)*k)


def _fixed_fit(y,n,X,offset,rho,start=None):
    initial=np.linalg.lstsq(X,logit((y+.5)/(n+1))-offset,rcond=None)[0] if start is None else np.asarray(start)
    const=np.sum(gammaln(n+1)-gammaln(y+1)-gammaln(n-y+1))
    def objective(beta):
        eta=X@beta+offset;p=np.clip(expit(eta),1e-10,1-1e-10)
        if rho==0:return float(np.sum(n*np.logaddexp(0,eta)-y*eta)-const),X.T@(n*p-y)
        k=(1-rho)/rho;a=p*k;b=(1-p)*k
        ll=_logpdf(y,n,p,rho)
        grad_p=k*(digamma(y+a)-digamma(a)-digamma(n-y+b)+digamma(b))
        return -float(ll.sum()),-X.T@(grad_p*p*(1-p))
    z=minimize(objective,np.clip(initial,-17,17),jac=True,method='L-BFGS-B',bounds=[(-18,18)]*X.shape[1],options={'maxiter':250,'ftol':1e-11,'gtol':1e-6,'maxls':30})
    return dict(beta=z.x,loglik=float(-z.fun),converged=bool(z.success),coefficient_boundary=bool(np.max(np.abs(z.x))>=17.99),optimizer_message=str(z.message))


def _null_profile(y,n,X0,offset,block,rho):
    # With event-constant effective-length offset, the null separates by experiment.
    # At fixed rho, BB loglik is strictly concave in mean p, so a bracketed score root is unique.
    if np.max(offset)-np.min(offset)>1e-12:
        return _fixed_fit(y,n,X0,offset,rho)
    groups=np.unique(block);probs=[]
    for g in groups:
        ix=block==g;yy=y[ix];nn=n[ix]
        if yy.sum()==0:p=1e-10
        elif yy.sum()==nn.sum():p=1-1e-10
        elif rho==0:p=float(yy.sum()/nn.sum())
        else:
            k=(1-rho)/rho
            def score(p):
                a=p*k;b=(1-p)*k
                return float(np.sum(digamma(yy+a)-digamma(a)-digamma(nn-yy+b)+digamma(b)))
            p=brentq(score,1e-10,1-1e-10,xtol=1e-12,rtol=1e-12)
        probs.append(p)
    eta=logit(probs)-offset[0];beta=np.r_[eta[0],eta[1:]-eta[0]]
    p=expit(X0@beta+offset)
    return dict(beta=beta,loglik=float(_logpdf(y,n,p,rho).sum()),converged=True,coefficient_boundary=bool(np.max(np.abs(beta))>=17.99),optimizer_message='independent strictly-concave mean-score roots')


def fit(y,n,t,offset,block=None):
    y,n,t=np.broadcast_arrays(np.asarray(y,float),np.asarray(n,float),np.asarray(t,float));offset=np.broadcast_to(np.asarray(offset,float),y.shape)
    bb=np.zeros(len(y),dtype=str) if block is None else np.asarray(block).astype(str)
    keep=np.isfinite(y+n+t+offset)&(n>0)&(y>=0)&(y<=n)
    y,n,t,offset,bb=[a[keep] for a in (y,n,t,offset,bb)]
    levels=sorted(set(bb));X0=design(t,bb,False,levels);X=design(t,bb,True,levels)
    df=len(y)-np.linalg.matrix_rank(X0)
    fail=dict(family='bb_pearson_df',p_asymptotic=1.,statistic=np.nan,beta=np.full(X.shape[1],np.nan),null_beta=np.full(X0.shape[1],np.nan),rho=np.nan,null_rho=np.nan,converged=False,loglik=np.nan,null_loglik=np.nan,block_levels=levels,n_samples=len(y),dispersion_residual_df=int(df),failure_reason='insufficient_design_or_splice_support')
    if len(y)<=X.shape[1]+1 or np.linalg.matrix_rank(X)<X.shape[1] or y.sum()==0 or y.sum()==n.sum():return fail
    cache={}
    def profile(rho):
        key=float(rho)
        if key not in cache:
            f=_null_profile(y,n,X0,offset,bb,key);p=expit(X0@f['beta']+offset)
            Q=float(np.sum((y-n*p)**2/np.maximum(n*p*(1-p)*(1+(n-1)*key),1e-15)))
            cache[key]=(Q-df,f,Q)
        return cache[key]
    q0,_,_=profile(0.)
    if q0<=0:rho=0.
    else:
        qmax,_,_=profile(RHO_MAX)
        if qmax>=0:return fail|dict(failure_reason='pearson_dispersion_no_root_below_upper_bound',pearson_q_at_zero=q0+df,pearson_q_at_upper=qmax+df)
        try:rho=float(brentq(lambda r:profile(r)[0],0.,RHO_MAX,xtol=2e-8,rtol=1e-8,maxiter=50))
        except (ValueError,RuntimeError) as e:return fail|dict(failure_reason='pearson_dispersion_root_failure:'+str(e))
    _,null,Q=profile(rho)
    nested=np.r_[null['beta'],0.,0.]
    alt=_fixed_fit(y,n,X,offset,rho,start=nested)
    raw=2*(alt['loglik']-null['loglik'])
    if (not alt['converged']) or raw < -1e-6:
        second=_fixed_fit(y,n,X,offset,rho)
        if second['converged'] and second['loglik']>alt['loglik']:alt=second
        raw=2*(alt['loglik']-null['loglik'])
    okay=bool(null['converged'] and alt['converged'] and not null['coefficient_boundary'] and not alt['coefficient_boundary'] and raw>=-1e-5)
    beta=alt['beta'];amp=np.hypot(beta[-2],beta[-1]);bases=np.r_[beta[0],beta[0]+beta[1:len(levels)]]
    return dict(family='bb_pearson_df',p_asymptotic=float(chi2.sf(max(0,raw),2)) if okay else 1.,statistic=float(max(0,raw)),raw_lr=float(raw),beta=beta,null_beta=null['beta'],rho=rho,null_rho=rho,converged=okay,loglik=alt['loglik'],null_loglik=null['loglik'],block_levels=levels,n_samples=len(y),dispersion_residual_df=int(df),pearson_q=Q,pearson_equation_error=float(Q-df),rho_boundary_zero=bool(rho==0),null_profile_evaluations=len(cache),coefficient_boundary=null['coefficient_boundary'] or alt['coefficient_boundary'],failure_reason='' if okay else 'optimizer_boundary_or_nesting_failure',peak_time=float(np.arctan2(beta[-1],beta[-2])%(2*np.pi)*24/(2*np.pi)),psi_peak_to_trough_blockmean=float(np.mean(expit(bases+amp)-expit(bases-amp))))


def predict_psi(t,fit_result,block=None,harmonic=True):
    return expit(design(t,block,harmonic,fit_result['block_levels'])@fit_result['beta' if harmonic else 'null_beta'])


def log_predictive(y,n,t,offset,fit_result,block=None,harmonic=True):
    p=expit(design(t,block,harmonic,fit_result['block_levels'])@fit_result['beta' if harmonic else 'null_beta']+offset)
    return _logpdf(np.asarray(y),np.asarray(n),p,fit_result['rho'])


def draw(rng,n,p,rho):
    if rho==0:return rng.binomial(np.asarray(n,int),p)
    k=(1-rho)/rho
    return rng.binomial(np.asarray(n,int),rng.beta(np.clip(p*k,1e-10,None),np.clip((1-p)*k,1e-10,None)))


def bootstrap(y,n,t,offset,block=None,n_boot=199,seed=20261019,observed=None):
    observed=fit(y,n,t,offset,block) if observed is None else observed
    n=np.asarray(n,int);t=np.asarray(t,float);offset=np.broadcast_to(np.asarray(offset,float),n.shape)
    if not observed['converged']:return dict(p_bootstrap=1.,bootstrap_valid=0,bootstrap_failures=0,bootstrap_failure_rate=np.nan,bootstrap_attempted=0,bootstrap_requested=n_boot,bootstrap_reliable=False)
    X=design(t,block,False,observed['block_levels']);p=expit(X@observed['null_beta']+offset);rng=np.random.default_rng(seed)
    stats=[];failures=0
    for _ in range(n_boot):
        yy=draw(rng,n,p,observed['rho'])
        try:
            r=fit(yy,n,t,offset,block)
            if r['converged']:stats.append(r['statistic'])
            else:failures+=1
        except (ValueError,RuntimeError,FloatingPointError):failures+=1
    exceed=sum(s>=observed['statistic'] for s in stats);upper=(1+exceed+failures)/(1+n_boot);reliable=failures/n_boot<=.05
    return dict(p_bootstrap=upper if reliable else 1.,p_bootstrap_upper=upper,p_bootstrap_valid_only=(1+exceed)/(1+len(stats)),bootstrap_valid=len(stats),bootstrap_failures=failures,bootstrap_failure_rate=failures/n_boot,bootstrap_attempted=n_boot,bootstrap_requested=n_boot,bootstrap_reliable=reliable)
