"""Development candidate: bootstrap bias-correct null rho, original v2 BB LRT.

Pilot99 estimates null-dispersion bias at the null MLE. rho_bc is clipped
2*rho_hat-mean(rho_pilot). Only the bootstrap-generating null changes. Neither
the observed statistic nor final simulated refit statistic is bias-corrected.
This is approximate plug-in calibration, not a double bootstrap or exact test.
"""
from __future__ import annotations
import numpy as np
from scipy.special import expit
from . import small_sample

from . import bb_model as base
PILOT_B=99
MAX_FAILURE_FRACTION=.05
RHO_BC_UPPER=.98
design=base.design
predict_psi=base.predict_psi
log_predictive=base.log_predictive


def fit(y,n,t,offset,block=None):
    r=base.fit(y,n,t,offset,block)
    r['family']='bb_bias_corrected'
    r['inference']='Original v2 BB LRT; use pilot99 bias-corrected-null bootstrap calibration'
    return r


def bootstrap(y,n,t,offset,block=None,n_boot=199,seed=173,observed=None):
    observed=fit(y,n,t,offset,block) if observed is None else observed
    result=dict(bootstrap_requested=int(n_boot),bootstrap_seed=int(seed),pilot_requested=PILOT_B,pilot_attempted=0,pilot_valid=0,pilot_failures=0,pilot_failure_rate=np.nan,bootstrap_attempted=0,bootstrap_valid=0,bootstrap_failures=0,bootstrap_failure_rate=np.nan,bootstrap_reliable=False,p_bootstrap=1.,p_bootstrap_upper=1.,p_bootstrap_valid_only=np.nan)
    if not observed['converged']:return result|dict(bootstrap_skipped_reason='observed_fit_invalid')
    y,n,t=np.broadcast_arrays(np.asarray(y,float),np.asarray(n,float),np.asarray(t,float));offset=np.broadcast_to(np.asarray(offset,float),y.shape)
    block=None if block is None else np.asarray(block).astype(str)
    finite=np.isfinite(y+n+t+offset)&(n>=0)&(y>=0)&(y<=n)
    y,n,t,offset=[v[finite] for v in (y,n,t,offset)]
    if block is not None:
        block=block[finite]
        known=np.isin(block,observed['block_levels'])
        y,n,t,offset,block=[v[known] for v in (y,n,t,offset,block)]
    positive=n>0
    nn=n[positive];tt=t[positive];oo=offset[positive]
    bb=np.zeros(len(nn),dtype=str) if block is None else block[positive]
    levels=sorted(set(bb));X0=small_sample.design(tt,bb,False,levels)
    X0all=base.design(t,block,False,observed['block_levels'])
    original_null_p=expit(X0all@observed['null_beta']+offset)
    streams=np.random.SeedSequence(int(seed)).spawn(2)
    pilot_rng=np.random.default_rng(streams[0]);final_rng=np.random.default_rng(streams[1])
    pilot_rhos=[];pilot_failures=0
    for _ in range(PILOT_B):
        yy=base.draw(pilot_rng,n,original_null_p,observed['null_rho'])[positive]
        try:
            if yy.sum()==0 or yy.sum()==nn.sum():pilot_failures+=1;continue
            f0=base._fit_design(yy,nn,X0,oo)
            if f0['converged'] and np.isfinite(f0['rho']):pilot_rhos.append(f0['rho'])
            else:pilot_failures+=1
        except (ValueError,FloatingPointError,np.linalg.LinAlgError):pilot_failures+=1
    result.update(pilot_attempted=PILOT_B,pilot_valid=len(pilot_rhos),pilot_failures=pilot_failures,pilot_failure_rate=pilot_failures/PILOT_B,null_rho_original=float(observed['null_rho']),rng_streams='SeedSequence(seed).spawn(2): pilot child0, final child1')
    if pilot_failures/PILOT_B>MAX_FAILURE_FRACTION or not pilot_rhos:
        return result|dict(bootstrap_skipped_reason='pilot_failure_rate_above_5pct')
    pilot_mean=float(np.mean(pilot_rhos));raw=2*observed['null_rho']-pilot_mean;corrected=float(np.clip(raw,0,RHO_BC_UPPER))
    result.update(pilot_mean_null_rho=pilot_mean,pilot_sd_null_rho=float(np.std(pilot_rhos,ddof=1)),estimated_null_rho_bias=pilot_mean-observed['null_rho'],null_rho_corrected_raw=float(raw),null_rho_corrected=corrected,rho_correction_clipped_low=bool(raw<0),rho_correction_clipped_high=bool(raw>RHO_BC_UPPER))
    corrected_mean=small_sample._null_profile(y[positive],nn,X0,oo,bb,corrected)
    if not corrected_mean['converged'] or corrected_mean['coefficient_boundary']:
        return result|dict(bootstrap_skipped_reason='corrected_null_mean_invalid')
    result['corrected_null_beta']=corrected_mean['beta'].tolist()
    corrected_p=expit(X0all@corrected_mean['beta']+offset)
    failures=0;stats=[]
    for _ in range(n_boot):
        yy=base.draw(final_rng,n,corrected_p,corrected)
        try:r=base.fit(yy,n,t,offset,block)
        except (ValueError,FloatingPointError,np.linalg.LinAlgError):failures+=1;continue
        if r['converged'] and np.isfinite(r['statistic']):stats.append(r['statistic'])
        else:failures+=1
    exceed=int(np.sum(np.asarray(stats)>=observed['statistic']-1e-9));upper=(1+exceed+failures)/(n_boot+1);reliable=bool(n_boot>0 and failures/max(n_boot,1)<=MAX_FAILURE_FRACTION)
    result.update(bootstrap_attempted=int(n_boot),bootstrap_valid=len(stats),bootstrap_failures=failures,bootstrap_failure_rate=failures/max(n_boot,1),bootstrap_reliable=reliable,p_bootstrap=float(upper) if reliable else 1.,p_bootstrap_upper=float(upper),p_bootstrap_valid_only=float((1+exceed)/(1+len(stats))),bootstrap_exceed=exceed)
    return result
