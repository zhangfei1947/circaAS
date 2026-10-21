"""Sample-level beta-binomial harmonic regression with conditional-depth bootstrap.

I | N follows BB(mean p, intraclass correlation rho). logit(p) is the
effective-length offset log(LI/LS), replicate intercepts and a fixed 24 h
harmonic. Dispersion is shared across samples in an event-condition group.
Zero-depth samples contribute no likelihood and are retained at zero during
bootstrap. The model conditions on observed event coverage, not gene TPM.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, logit, gammaln, betaln, digamma
from scipy.stats import chi2


def design(t, block=None, harmonic=True, levels=None):
    t = np.asarray(t, float)
    cols = [np.ones(len(t))]
    if block is not None:
        block = np.asarray(block).astype(str)
        levels = sorted(set(block)) if levels is None else list(levels)
        if not set(block).issubset(levels):
            raise ValueError('Prediction block not present in training data')
        cols += [(block == level).astype(float) for level in levels[1:]]
    if harmonic:
        cols += [np.cos(2*np.pi*t/24), np.sin(2*np.pi*t/24)]
    return np.column_stack(cols)


def _logpdf(y, n, p, rho):
    p = np.clip(p, 1e-10, 1-1e-10)
    const = gammaln(n+1)-gammaln(y+1)-gammaln(n-y+1)
    if rho == 0:
        return const + y*np.log(p)+(n-y)*np.log1p(-p)
    k = (1-rho)/rho
    return const+betaln(y+p*k, n-y+(1-p)*k)-betaln(p*k,(1-p)*k)


def _fit_design(y, n, X, offset, nested=None):
    b0 = np.clip(np.linalg.lstsq(X, logit((y+.5)/(n+1))-offset, rcond=None)[0], -12,12)
    def binobj(beta):
        eta = X@beta+offset
        return np.sum(n*np.logaddexp(0,eta)-y*eta), X.T@(n*expit(eta)-y)
    bf = minimize(binobj,b0,jac=True,method='L-BFGS-B',bounds=[(-18,18)]*len(b0),
                  options={'maxiter':300,'ftol':1e-11,'gtol':1e-6})
    candidates=[]
    if bf.success:
        candidates.append(dict(beta=bf.x,rho=0.,ll=float(_logpdf(y,n,expit(X@bf.x+offset),0).sum()),
                               converged=True,dispersion_boundary=True,optimizer='binomial_boundary'))
    def objective(theta):
        p = np.clip(expit(X@theta[:-1]+offset),1e-10,1-1e-10)
        k=np.exp(theta[-1]); a=p*k; b=(1-p)*k
        ll=_logpdf(y,n,p,1/(1+k))
        ga=digamma(y+a)-digamma(n+k)-digamma(a)+digamma(k)
        gb=digamma(n-y+b)-digamma(n+k)-digamma(b)+digamma(k)
        gradient=np.r_[X.T@(k*(ga-gb)*p*(1-p)),np.sum(k*(p*ga+(1-p)*gb))]
        return -float(ll.sum()),-gradient
    starts=[np.r_[bf.x,np.log(20.)],np.r_[bf.x,np.log(200.)]]
    if nested is not None:
        starts.insert(0,np.r_[nested['beta'],np.log(200. if nested['rho']==0 else (1-nested['rho'])/nested['rho'])])
    bounds=[(-18,18)]*len(b0)+[(np.log(.02),np.log(1e5))]
    for start in starts:
        ff=minimize(objective,start,jac=True,method='L-BFGS-B',bounds=bounds,
                    options={'maxiter':300,'ftol':1e-10,'gtol':1e-6})
        if ff.success and np.isfinite(ff.fun):
            candidates.append(dict(beta=ff.x[:-1],rho=float(1/(1+np.exp(ff.x[-1]))),ll=-float(ff.fun),
                converged=True,dispersion_boundary=bool(ff.x[-1] >= np.log(1e5)-1e-4),optimizer='L-BFGS-B'))
    if not candidates:
        return dict(beta=bf.x,rho=np.nan,ll=-np.inf,converged=False,dispersion_boundary=False,optimizer='failed')
    best=max(candidates,key=lambda r:r['ll'])
    best['coefficient_boundary']=bool(np.any(np.abs(best['beta'])>=17.99))
    best['dispersion_upper_boundary']=bool(best['rho']>=1/(1+.02001))
    best['converged'] &= not best['coefficient_boundary'] and not best['dispersion_upper_boundary']
    return best


def fit(y,n,t,offset,block=None):
    y,n,t=np.broadcast_arrays(np.asarray(y,float),np.asarray(n,float),np.asarray(t,float))
    offset=np.broadcast_to(np.asarray(offset,float),y.shape)
    valid=np.isfinite(y+n+t+offset)&(n>0)&(y>=0)&(y<=n)
    levels=[] if block is None else sorted(set(np.asarray(block).astype(str)[valid]))
    b=None if block is None else np.asarray(block).astype(str)[valid]
    yy,nn,tt,oo=[a[valid] for a in (y,n,t,offset)]
    X0=design(tt,b,False,levels); X=design(tt,b,True,levels)
    fail=dict(family='bb',p_asymptotic=1.,statistic=np.nan,beta=np.full(X.shape[1],np.nan),
        null_beta=np.full(X0.shape[1],np.nan),rho=np.nan,null_rho=np.nan,converged=False,
        loglik=np.nan,null_loglik=np.nan,block_levels=levels,n_samples=int(valid.sum()),
        n_zero_depth=int((n==0).sum()),failure_reason='insufficient_design_or_splice_support')
    if len(yy)<=X.shape[1]+1 or np.linalg.matrix_rank(X)<X.shape[1] or yy.sum()==0 or (nn-yy).sum()==0:
        return fail
    null=_fit_design(yy,nn,X0,oo)
    nested=dict(beta=np.r_[null['beta'],0.,0.],rho=null['rho']) if null['converged'] else None
    alt=_fit_design(yy,nn,X,oo,nested=nested)
    raw_lr=2*(alt['ll']-null['ll'])
    converged=bool(null['converged'] and alt['converged'] and np.isfinite(raw_lr) and raw_lr>=-1e-5)
    stat=max(0.,raw_lr) if np.isfinite(raw_lr) else np.nan
    beta=alt['beta']; amp=float(np.hypot(*beta[-2:])); nblocks=max(1,len(levels))
    bases=np.r_[beta[0], beta[0]+beta[1:nblocks]]
    return dict(family='bb',p_asymptotic=float(chi2.sf(stat,2)) if converged else 1.,statistic=stat,
        raw_lr=raw_lr,beta=beta,null_beta=null['beta'],rho=alt['rho'],null_rho=null['rho'],
        converged=converged,loglik=alt['ll'],null_loglik=null['ll'],block_levels=levels,
        n_samples=int(valid.sum()),n_zero_depth=int((n==0).sum()),
        peak_time=float(np.arctan2(beta[-1],beta[-2])%(2*np.pi)*24/(2*np.pi)),
        psi_peak_to_trough_blockmean=float(np.mean(expit(bases+amp)-expit(bases-amp))),
        dispersion_boundary=alt['dispersion_boundary'],null_dispersion_boundary=null['dispersion_boundary'],
        coefficient_boundary=alt.get('coefficient_boundary',False) or null.get('coefficient_boundary',False),
        failure_reason='' if converged else 'optimizer_boundary_or_nesting_failure')


def predict_psi(t,fit_result,block=None,harmonic=True):
    beta=fit_result['beta' if harmonic else 'null_beta']
    return expit(design(t,block,harmonic,fit_result['block_levels'])@beta)


def log_predictive(y,n,t,offset,fit_result,block=None,harmonic=True):
    beta=fit_result['beta' if harmonic else 'null_beta']
    eta=design(t,block,harmonic,fit_result['block_levels'])@beta+np.asarray(offset,float)
    return _logpdf(np.asarray(y,float),np.asarray(n,float),expit(eta),fit_result['rho' if harmonic else 'null_rho'])


def draw(rng,n,p,rho):
    if rho==0:
        return rng.binomial(np.asarray(n,int),p)
    k=(1-rho)/rho
    return rng.binomial(np.asarray(n,int),rng.beta(np.clip(p*k,1e-10,None),np.clip((1-p)*k,1e-10,None)))


def bootstrap(y,n,t,offset,block=None,n_boot=199,seed=173,observed=None,max_failure_fraction=.05):
    """Fit-specific parametric calibration, not an external model-comparison DGP.

    Failed replicates are excluded from valid statistics and explicitly counted.
    Upper-bound p treats all failed replicates as exceedances. Primary p uses that
    conservative bound; >5% failures or an invalid observed fit yields p=1.
    """
    y,n,t=np.broadcast_arrays(np.asarray(y,float),np.asarray(n,float),np.asarray(t,float))
    offset=np.broadcast_to(np.asarray(offset,float),y.shape)
    observed=fit(y,n,t,offset,block) if observed is None else observed
    base=dict(bootstrap_requested=int(n_boot),bootstrap_seed=int(seed))
    if not observed['converged']:
        return base|dict(p_bootstrap=1.,p_bootstrap_valid_only=np.nan,p_bootstrap_upper=1.,
            bootstrap_attempted=0,bootstrap_valid=0,bootstrap_failures=0,bootstrap_failure_rate=np.nan,
            bootstrap_exceed=0,bootstrap_reliable=False,bootstrap_skipped_reason='observed_fit_invalid')
    finite=np.isfinite(y+n+t+offset)&(n>=0)&(y>=0)&(y<=n)
    y,n,t,offset=[v[finite] for v in (y,n,t,offset)]
    block=None if block is None else np.asarray(block).astype(str)[finite]
    # A block represented only by n=0 samples has no fitted intercept: those rows
    # can be discarded because their generated y is deterministically zero.
    if block is not None:
        ok=np.isin(block,observed['block_levels'])
        y,n,t,offset,block=[v[ok] for v in (y,n,t,offset,block)]
    p0=expit(design(t,block,False,observed['block_levels'])@observed['null_beta']+offset)
    rng=np.random.default_rng(seed); vals=[]; failures=0
    for _ in range(n_boot):
        yy=draw(rng,n,p0,observed['null_rho'])
        try:
            rr=fit(yy,n,t,offset,block)
            if rr['converged'] and np.isfinite(rr['statistic']):
                vals.append(rr['statistic'])
            else:
                failures+=1
        except (ValueError,FloatingPointError,np.linalg.LinAlgError):
            failures+=1
    exceed=int(np.sum(np.asarray(vals)>=observed['statistic']-1e-9))
    upper=(1+exceed+failures)/(1+n_boot)
    reliable=bool(n_boot>0 and failures/max(n_boot,1)<=max_failure_fraction)
    return base|dict(p_bootstrap=float(upper) if reliable else 1.,
        p_bootstrap_valid_only=float((1+exceed)/(1+len(vals))) if vals else np.nan,
        p_bootstrap_upper=float(upper),bootstrap_attempted=int(n_boot),bootstrap_valid=len(vals),bootstrap_failures=failures,
        bootstrap_failure_rate=float(failures/max(n_boot,1)),bootstrap_exceed=exceed,bootstrap_reliable=reliable)
