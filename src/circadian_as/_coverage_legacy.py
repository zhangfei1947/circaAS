"""Coverage-only operational identifiability gates, fixed before outcome fitting.

These are design diagnostics, NOT guarantees of power or absence of rhythm.
No inclusion/skipping proportion enters the primary gate.
"""
import numpy as np

RULE={'period':24.,'min_positive_samples':10,'min_positive_times':6,
      'min_positive_samples_each_repeat':4,'min_total_count':100,
      'min_distinct_phases':4,'max_circular_gap_h':12.,
      'min_effective_samples':8.,'max_harmonic_condition':30.,
      'sample_depth_cap_for_design':50.}

def assess_coverage(n,t,block,period=24.):
    n=np.asarray(n,float);t=np.asarray(t,float);block=np.asarray(block).astype(str)
    assert n.shape==t.shape==block.shape and np.isfinite(n).all() and (n>=0).all()
    pos=n>0;times=np.unique(t);levels=np.unique(block);ph=np.unique(np.round(t[pos]%period,8))
    gaps=np.diff(np.r_[ph,ph[0]+period]) if len(ph) else np.array([period])
    w=np.minimum(n,RULE['sample_depth_cap_for_design']);neff=w.sum()**2/(w@w) if (w@w)>0 else 0.
    nuisance=np.column_stack([np.ones(len(t))]+[(block==l).astype(float) for l in levels[1:]])
    h=np.column_stack([np.cos(2*np.pi*t/period),np.sin(2*np.pi*t/period)])
    design=np.column_stack([nuisance,h]);rank=int(np.linalg.matrix_rank(design[pos])) if pos.any() else 0
    sw=np.sqrt(w);a=nuisance*sw[:,None];b=h*sw[:,None]
    residual=b-a@np.linalg.lstsq(a,b,rcond=None)[0]
    ev=np.linalg.eigvalsh(residual.T@residual)
    cond=float(ev[-1]/ev[0]) if ev[0]>1e-8 else float('inf')
    counts_per_rep=[int((pos&(block==l)).sum()) for l in levels]
    time_depth=np.array([n[t==x].sum() for x in times])
    reasons=[]
    if pos.sum()<10:reasons.append('fewer_than_10_positive_samples')
    if len(np.unique(t[pos]))<6:reasons.append('fewer_than_6_positive_absolute_times')
    if len(levels)<2 or min(counts_per_rep,default=0)<4:reasons.append('repeat_coverage_insufficient')
    if n.sum()<100:reasons.append('total_informative_count_below_100')
    if len(ph)<4:reasons.append('fewer_than_4_observed_24h_phases')
    if gaps.max()>12.+1e-8:reasons.append('circular_phase_gap_exceeds_12h')
    if neff<8:reasons.append('effective_sample_count_below_8')
    if rank<design.shape[1]:reasons.append('harmonic_design_rank_deficient')
    if cond>30:reasons.append('harmonic_design_ill_conditioned')
    return dict(sufficient=not reasons,reason=';'.join(reasons) or 'passes_operational_coverage_gate',
        n_samples=len(n),n_positive_samples=int(pos.sum()),n_positive_times=int((time_depth>0).sum()),
        n_zero_times=int((time_depth==0).sum()),n_low_times=int((time_depth<20).sum()),
        n_distinct_phases=len(ph),max_phase_gap=float(gaps.max()),total_count=int(n.sum()),
        min_positive_each_repeat=min(counts_per_rep,default=0),effective_samples=float(neff),
        harmonic_condition=cond,design_rank=rank)

def coverage_eligibility(times,N,reps):
    return assess_coverage(N,times,reps)
