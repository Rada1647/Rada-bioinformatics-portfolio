"""Independent reviewer-1 checks; does not execute the authors' self-tests.

Writes reviewer evidence only. Production inputs/results are read-only.
"""
from pathlib import Path
import sys, json, hashlib, itertools, math, datetime
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import methods as impl
import benchmark_extension as bench
import prepare_extension as prep

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def independent_resolvent(z, t):
    """Closed cosine eigenbasis of the reflecting path, not a band solve."""
    n = len(z)
    k = np.arange(n)
    u = np.cos(np.pi * (np.arange(n)[:, None] + .5) * k[None, :] / n)
    u[:, 0] = 1 / np.sqrt(n)
    u[:, 1:] *= np.sqrt(2 / n)
    return u @ ((u.T @ z) / (1 + t * np.sin(np.pi * k / (2*n))**2)[:, None])

def reflect(i, n):
    v = i % (2 * n)
    return v if v < n else 2*n - 1 - v

def independent_filter(z, name, params):
    if name == 'moving_average':
        radius = params['half_width']
        weights = np.ones(2*radius+1)/(2*radius+1)
    elif name == 'gaussian':
        sigma = params['sigma']
        radius = int(4*sigma + .5)
        weights = np.exp(-np.arange(-radius,radius+1)**2/(2*sigma*sigma))
        weights /= weights.sum()
    else:
        radius = params['width']//2
    out = np.zeros_like(z, dtype=float)
    for i in range(len(z)):
        x = z[[reflect(j,len(z)) for j in range(i-radius,i+radius+1)]]
        out[i] = np.sort(x,axis=0)[radius] if name == 'median' else weights @ x
    return out

def reference_predict(z, name, params, epsilon=None):
    z = z.astype(float)
    if name == 'raw': return z
    if name in ('moving_average','gaussian','median'):
        return independent_filter(z, name, params)
    h = np.clip(independent_resolvent(z,params['t']),0,1)
    h[:,np.all(z==0,axis=0)] = 0
    h[:,np.all(z==1,axis=0)] = 1
    if name == 'linear': return h
    if name == 'paper332': return 3*h*h-2*h*h*h
    if name == 'hard_linear': return (h>=.5).astype(float)
    if name == 'oracle_debiased_linear': return np.clip((h-epsilon)/(1-2*epsilon),0,1)
    a = params['alpha']
    return h**a/(h**a+(1-h)**a)

def independent_masks(x, f):
    rare = f <= .1
    short = np.zeros_like(x,dtype=bool)
    trans = np.zeros_like(x,dtype=bool)
    events=0
    for j in range(x.shape[1]):
        i=0
        while i<len(x):
            if not x[i,j]: i+=1; continue
            start=i
            while i<len(x) and x[i,j]: i+=1
            if start>0 and i<len(x) and i-start<=3:
                short[start:i,j]=True; events+=1
        for i in range(1,len(x)):
            if x[i,j] != x[i-1,j]: trans[i-1:i+1,j]=True
    return {'rare_features':rare,'rare_positive':x & rare[None,:],
            'short_positive':short,'transition':trans,'brief_event_count':events}

def independent_metrics(p,x,m):
    error=p.astype(float)-x.astype(float)
    yes=p>=.5
    tp=sum(int(a and b) for a,b in zip(yes.flat,x.flat))
    fp=sum(int(a and not b) for a,b in zip(yes.flat,x.flat))
    fn=sum(int(not a and b) for a,b in zip(yes.flat,x.flat))
    def div(a,b): return a/b if b else None
    def avg(a): return float(np.mean(a)) if a.size else None
    return {'brier':float(np.mean(error**2)),'mae':float(np.mean(abs(error))),
            'occupancy_mae':float(np.mean([abs(float(np.mean(p[:,j]))-float(np.mean(x[:,j]))) for j in range(x.shape[1])])),
            'occupancy_signed_bias':float(np.sum(error)/error.size),
            'precision':div(tp,tp+fp),'recall':div(tp,tp+fn),'f1':div(2*tp,2*tp+fp+fn),
            'rare_brier':avg(error[:,m['rare_features']]**2),
            'rare_positive_recall':avg(yes[m['rare_positive']]),
            'short_positive_recall':avg(yes[m['short_positive']]),
            'transition_brier':avg(error[m['transition']]**2)}

def assert_metrics(a,b,tol=1e-12):
    errors=[]
    for k in a:
        assert (a[k] is None)==(b[k] is None), (k,a[k],b[k])
        if a[k] is not None:
            errors.append(abs(a[k]-b[k])); assert errors[-1]<=tol,(k,a[k],b[k])
    return max(errors,default=0)

def independent_noise(x,xyz,pairs,s,words):
    rng=np.random.Generator(np.random.PCG64(np.random.SeedSequence(words)))
    e=s['level']
    if s['kind']=='bitflip': return x ^ (rng.uniform(size=x.shape)<e)
    if s['kind']=='correlated':
        # Consume a matrix in row-major order, and use explicit 0/1 transitions.
        uniform=rng.uniform(size=x.shape)
        err=uniform.copy()<e
        for i in range(1,len(x)):
            for j in range(x.shape[1]):
                err[i,j] = uniform[i,j] < (1-(1-e)*(1-s['rho']) if err[i-1,j] else e*(1-s['rho']))
        return x ^ err
    r=xyz.astype(float)+rng.normal(size=xyz.shape)*e
    return np.stack([np.sum((r[:,a]-r[:,b])**2,axis=1)<s['cutoff_A']**2 for a,b in pairs],axis=1)

def run():
    ev={}
    rng=np.random.default_rng(9127301)
    dense_error=0.; spectral_error=0.; median_error=0.; path_cases=0
    for n in range(1,8):
        z=np.asarray(list(itertools.product((0.,1.),repeat=n))).T
        a=np.eye(n)
        if n>1:
            a=np.zeros((n,n))
            for i in range(n-1): a[i,i+1]=a[i+1,i]=.25
            for i in range(n): a[i,i]=1-a[i].sum()
        assert np.array_equal(a,a.T) and np.all(a.sum(axis=1)==1)
        for t in (1,2,4,8,16,32,64):
            r=np.linalg.inv(np.eye(n)+t*(np.eye(n)-a))
            assert r.min()>-1e-14 and np.max(abs(r.sum(axis=1)-1))<1e-14
            h=impl.predict(z,'linear',{'t':t})
            dense_error=max(dense_error,float(np.max(abs(h-r@z))))
            spectral_error=max(spectral_error,float(np.max(abs(h-independent_resolvent(z,t)))))
            # Enumerate all endpoint triples for an expectation independent of cubic code.
            if n<=4:
                expected=np.zeros_like(z)
                for p,q,s in itertools.product(range(n),repeat=3):
                    majority=(z[p]+z[q]+z[s]>=2)
                    expected += (r[:,p]*r[:,q]*r[:,s])[:,None]*majority[None,:]
                y=impl.predict(z,'paper332',{'t':t})
                median_error=max(median_error,float(np.max(abs(y-expected))))
            path_cases+=z.shape[1]
    assert dense_error<2e-14 and spectral_error<2e-14 and median_error<2e-14
    ev['resolvent_and_expected_median']={'binary_path_parameter_cases':path_cases,'lengths':list(range(1,8)),
        't_grid':[1,2,4,8,16,32,64],'max_dense_matrix_error':dense_error,'max_cosine_spectrum_error':spectral_error,
        'max_triple_endpoint_enumeration_error':median_error,'endpoint_triple_enumeration_lengths':[1,2,3,4],
        'formula':'R=(I+t(I-A))^-1=(1/(1+t))*sum_{k>=0}(t/(1+t))^k A^k; mean geometric step count t. Three independent endpoint binary labels have majority probability 3h^2-2h^3.'}
    filter_error=0.; constant_error=0.; filter_cases=0; median_defects=[]
    for n in (1,2,3,5,12,20,71):
        z=np.column_stack([rng.integers(0,2,(n,4)),np.zeros(n),np.ones(n)]).astype(float)
        for name,grid in impl.candidates().items():
            for p in grid:
                actual=impl.predict(z,name,p,noise_level=.1)
                constant_error=max(constant_error,float(np.max(abs(actual[:,-2:]-z[:,-2:]))))
                assert actual.dtype==np.float64
                if name!='hard_linear':
                    case_error=float(np.max(abs(actual-reference_predict(z,name,p,.1))))
                    if case_error>2e-13 and name=='median':
                        median_defects.append({'n_frames':n,'parameters':p,'max_error':case_error,
                            'nonbinary_output_count':int(np.count_nonzero((actual!=0)&(actual!=1)))})
                    else:
                        filter_error=max(filter_error,case_error)
                filter_cases+=1
    assert filter_error<5e-11 and constant_error<1e-14,(filter_error,constant_error)
    grid=np.r_[np.linspace(0,1,10001),np.nextafter(.5,0),.5,np.nextafter(.5,1)]
    power_error=0.
    for a in (.5,.75,1,1.25,1.5,2,3,4,8):
        direct=grid**a/(grid**a+(1-grid)**a)
        power_error=max(power_error,float(np.max(abs(impl._power(grid,a)-direct))))
    assert power_error<1e-14
    ev['boundaries_constants_power']={'candidate_input_cases':filter_cases,'manual_half_sample_reflection':True,
        'max_reference_probability_error':filter_error,'max_constant_probability_error':constant_error,
        'max_direct_power_error':power_error,'hard_near_ties_compared_separately':True,'median_short_path_defects':median_defects,
        'cubic_threshold_disagreements_on_10004_values':int(np.count_nonzero((grid>=.5)!=(impl._cubic(grid)>=.5))),
        'classification_fact':'phi(h)-1/2=(h-1/2)*(3/2-2*(h-1/2)^2); sign preserved on [0,1] in exact arithmetic.'}
    actual_length_median=[]
    for n in (71,249,250,5001):
        z=np.column_stack([rng.integers(0,2,(n,7)),np.arange(n)%2,np.zeros(n),np.ones(n)]).astype(float)
        for w in (3,5,9,17,33,65,129):
            a=impl.predict(z,'median',{'width':w}); b=independent_filter(z,'median',{'width':w})
            assert np.array_equal(a,b),(n,w)
            actual_length_median.append({'n_frames':n,'width':w,'columns':z.shape[1],'exact_match':True})
    ev['boundaries_constants_power']['actual_length_manual_median_checks']=actual_length_median
    metric_err=0.; cases=0
    for n in (1,2,6,11):
        for _ in range(15):
            x=rng.integers(0,2,(n,5)).astype(bool)
            f=np.array([.05,.1,.10000001,.4,.8])
            p=rng.choice(np.r_[0,.2,np.nextafter(.5,0),.5,np.nextafter(.5,1),.9,1],size=x.shape)
            m=independent_masks(x,f); actual=bench.strata(x,f)
            for k in m: assert np.array_equal(m[k],actual[k]),k
            metric_err=max(metric_err,assert_metrics(independent_metrics(p,x,m),bench.metrics(p,x,actual)))
            cases+=1
    # Explicit zero-denominator conventions and strict 8 A contacts.
    for val in (False,True):
        x=np.full((4,2),val); f=np.array([.3,.6])
        metric_err=max(metric_err,assert_metrics(independent_metrics(x,x,independent_masks(x,f)),bench.metrics(x.astype(float),x,bench.strata(x,f))))
    xyz=np.zeros((3,5,3),np.float32); xyz[:,4,0]=[np.nextafter(np.float32(8),np.float32(0)),8,np.nextafter(np.float32(8),np.float32(9))]
    expected=np.array([[True],[False],[False]])
    assert np.array_equal(prep.contact_matrix(xyz,np.array([[0,4]]),8),expected)
    assert np.array_equal(bench.contact_map(xyz,np.array([[0,4]]),8),expected)
    assert np.array_equal(prep.all_pairs(6),np.array([[0,4],[0,5],[1,5]]))
    ev['metrics_and_contact_threshold']={'random_adversarial_cases':cases,'max_metric_error':metric_err,
        'boundary_censored_short_runs_checked':True,'overlapping_transition_union_checked':True,
        'rare_threshold_inclusive_point_1_checked':True,'empty_strata_null_checked':True,'strict_distance_cutoff_checked':True}
    protocol=json.loads((ROOT/'protocol/protocol_v1.1.json').read_text())
    for si,s in enumerate(protocol['scenarios']):
        x=rng.integers(0,2,(17,3)).astype(bool); coords=rng.normal(size=(17,6,3)); pairs=prep.all_pairs(6)
        words=[20261007,2,3,314,si,7]
        assert words==bench.seed_words(2,3,314,si,7)
        assert np.array_equal(independent_noise(x,coords,pairs,s,words),bench.corrupt(x,coords,pairs,s,words))
    ev['noise']={'all_8_scenarios_exact_array_match':True,'explicit_PCG64_SeedSequence':True,
        'stationary_correlated_proof':'With p01=e(1-rho), p10=(1-e)(1-rho), stationary mass p01/(p01+p10)=e and second eigenvalue 1-p01-p10=rho. Independent per-contact random inputs implement independent error chains.',
        'jitter_variance':'level is per-atom/component SD; relative coordinate noise is shared across contacts involving one atom, as intended.'}
    # Recompute selection and every reported aggregate from immutable stored rows.
    mp=ROOT/'data/prepared/villin/manifest.json'; manifest=json.loads(mp.read_text())
    runs={r['id']:r for r in manifest['runs']}
    summaries=[]; agg_error=0.; validation_error=0.; row_count=0; output_hashes={}
    sampled_error=0.; sampled_checks=[]; tie_effects=[]
    for si,s in enumerate(protocol['scenarios']):
        paths=[ROOT/'results/villin'/f"{s['name']}_{phase}.json" for phase in ('validation','test')]
        v,r=[json.loads(p.read_text()) for p in paths]
        output_hashes.update({str(p.relative_to(ROOT)):sha(p) for p in paths})
        for doc in (v,r):
            assert doc['protocol_sha256']==sha(ROOT/'protocol/protocol_v1.1.json')
            assert doc['prepared_manifest_sha256']==sha(mp)
            for fn,h in doc['code_sha256'].items(): assert sha(ROOT/fn)==h
        assert r['validation_file_sha256']==sha(paths[0])
        assert len(v['records'])==24*8 and len(r['per_run'])==48
        best={}
        for method,candidates in v['candidate_grids'].items():
            mean=np.mean([rec['candidate_brier'][method] for rec in v['records']],axis=0)
            validation_error=max(validation_error,float(np.max(abs(mean-v['mean_candidate_brier'][method]))))
            pick=min(range(len(mean)),key=lambda i:mean[i])
            assert candidates[pick]==v['selected_parameters'][method]
            best[method]=float(mean[pick])
        baseline_order=['raw','linear','hard_linear','moving_average','gaussian','median','power_full','power_matched']
        assert min(baseline_order,key=lambda k:best[k])==v['selected_baseline']
        assert min([k for k in baseline_order if k!='power_full'],key=lambda k:best[k])==v['selected_matched_baseline']
        for run in r['per_run']:
            assert len(run['noise_rows'])==16
            for ni,nrow in enumerate(run['noise_rows']):
                assert nrow['seed_words']==[20261007,2,3,run['run_index'],si,ni]
                row_count+=1
            for name in r['selected_parameters']:
                for metric,target in run['aggregate'][name].items():
                    vals=[row['methods'][name][metric] for row in run['noise_rows'] if row['methods'][name][metric] is not None]
                    assert len(vals)==run['aggregate_contributing_noise_count'][name][metric]
                    if vals: agg_error=max(agg_error,abs(float(np.mean(vals))-target))
                    else: assert target is None
        for name in r['selected_parameters']:
            for metric,target in r['aggregate'][name].items():
                vals=[run['aggregate'][name][metric] for run in r['per_run'] if run['aggregate'][name][metric] is not None]
                assert len(vals)==r['aggregate_contributing_run_count'][name][metric]
                if vals: agg_error=max(agg_error,abs(float(np.mean(vals))-target))
                else: assert target is None
        # Two selected trajectories: 250 and 71 frames; one independent noise draw each.
        for ri in (0,1):
            result_run=r['per_run'][ri]; run=runs[result_run['run_id']]
            with np.load(mp.parent/run['contacts'][str(s['cutoff_A'])]) as d:
                x=d['X']; pairs=d['pairs']; f=d['frequency']
            with np.load(mp.parent/run['coordinates_path']) as d: coords=d['xyz']
            words=[20261007,2,3,run['run_index'],si,0]
            z=independent_noise(x,coords,pairs,s,words)
            masks=independent_masks(x,f)
            for name,p in r['selected_parameters'].items():
                pred=reference_predict(z,name,p,s['level'])
                got=independent_metrics(pred,x,masks)
                saved=result_run['noise_rows'][0]['methods'][name]
                production=impl.predict(z,name,p,noise_level=s['level'])
                disagreement=(pred>=.5)!=(production>=.5)
                if np.any(disagreement):
                    before=got.copy()
                    if name=='hard_linear':
                        h1=independent_resolvent(z.astype(float),p['t'])
                        h2=impl.predict(z,'linear',{'t':p['t']})
                        tie_distance=max(float(np.max(abs(h1[disagreement]-.5))),float(np.max(abs(h2[disagreement]-.5))))
                    else:
                        tie_distance=max(float(np.max(abs(pred[disagreement]-.5))),float(np.max(abs(production[disagreement]-.5))))
                    assert tie_distance<1e-12,(s['name'],name,tie_distance)
                    pred=pred.copy();pred[disagreement]=production[disagreement]
                    got=independent_metrics(pred,x,masks)
                    tie_effects.append({'scenario':s['name'],'run':run['id'],'noise_index':0,'method':name,
                        'classification_disagreements':int(disagreement.sum()),'max_distance_from_half':tie_distance,
                        'metrics_before_tie_alignment':before,'saved_metrics':saved})
                sampled_error=max(sampled_error,assert_metrics(got,saved,2e-12))
                sampled_checks.append({'scenario':s['name'],'run':run['id'],'noise_index':0,'method':name})
        summaries.append({'scenario':s['name'],'selected_baseline':r['selected_baseline'],
            'paper332_brier':r['aggregate']['paper332']['brier'],
            'baseline_brier':r['aggregate'][r['selected_baseline']]['brier'],
            'same_t_threshold_disagreements':r['same_t_threshold_disagreements']})
    assert agg_error<1e-14 and validation_error<1e-14
    ev['villin_audit']={'scenarios':summaries,'test_noise_rows_checked':row_count,'max_stored_aggregate_error':agg_error,
        'max_validation_candidate_mean_error':validation_error,'independently_rebuilt_method_metric_rows':len(sampled_checks),
        'max_independent_rebuilt_metric_error':sampled_error,'sampled_reconstruction_scope':sampled_checks,
        'classification_roundoff_effects':tie_effects,
        'tie_policy_in_review':'Where independent spectral and production band solves produce opposite classifications only within 1e-12 of .5, report raw discrepancies then align those cells to production threshold convention before metric comparison.',
        'output_sha256':output_hashes}
    clean_path=ROOT/'results/villin/primary_clean.json'; clean=json.loads(clean_path.read_text())
    ev['villin_audit']['output_sha256'][str(clean_path.relative_to(ROOT))]=sha(clean_path)
    cleanerr=0.
    for name in clean['selected_parameters']:
        for metric,target in clean['aggregate'][name].items():
            vals=[r['methods'][name][metric] for r in clean['per_run'] if r['methods'][name][metric] is not None]
            assert len(vals)==clean['aggregate_contributing_run_count'][name][metric]
            if vals: cleanerr=max(cleanerr,abs(float(np.mean(vals))-target))
            else: assert target is None
    ev['villin_audit']['max_clean_aggregate_error']=cleanerr
    assert clean['oracle_epsilon']==.1
    ev['theorem_scope']={'brier_guarantee':False,'physical_coordinate_reconstruction_guarantee':False,
        'reason':'Squared whole-vector L1 distances in a metric Markov cotype inequality do not equal mean coordinate squared error against clean labels. The pipeline outputs marginal contact probabilities, with no Euclidean embedding or physical-force constraints.',
        'elementary_counterexample':'At h=.6 and clean label x=0, linear Brier=.36 whereas cubic=.648 gives Brier=.419904; cubic is not uniformly Brier-improving.',
        'paper_proof_review':'Not a line-by-line proof verification of the cited manuscript; review checks the implemented binary reduction and claim boundary.'}
    result={'reviewer':'Independent AI technical audit 1: mathematics and implementation','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'status':'initial_complete_pending_t4','production_code_implemented_by_reviewer':False,
        'scope':'protocol v1.1; methods.py; benchmark_extension.py; prepare_extension.py; all eight completed Villin validation/test outputs and primary clean aggregate; independent small-array mathematical checks plus 156 held-out method/noise metric reconstructions. No T4 outcome audit yet.',
        'reviewed_code_sha256':{fn:sha(ROOT/fn) for fn in ('protocol/protocol_v1.1.json','methods.py','benchmark_extension.py','prepare_extension.py')},
        'checker_sha256':sha(__file__),'verdict':'No critical discrepancy affecting reviewed benchmark outputs detected; off-protocol very-short-path SciPy median edge defect documented.',
        'findings':[{'severity':'scope_limit','id':'M1','finding':'Cotype theorem does not establish denoising superiority, Brier improvement, equilibrium populations, or physical coordinate validity. Existing protocol and code disclaimers are appropriate.'},
                    {'severity':'reporting_note','id':'M2','finding':'Short-positive recall is weighted by positive cells within complete sampled runs, not by event count; transition Brier uses a union of both endpoints. These match the clarified protocol.'},
                    {'severity':'reporting_note','id':'M3','finding':'Power alpha=3 is a prior-ADK selected control and full power grid has 63 settings versus seven for cubic. Both declared full-budget and matched-budget comparisons must remain visible.'},
                    {'severity':'limitation','id':'M4','finding':'Independent data-to-metric reconstruction samples two Villin trajectories and one noise draw per scenario, while stored aggregation/selection checks cover all completed rows. Full numerical reproduction is a separate reviewer task.'},
                    {'severity':'low_unused_domain_defect','id':'M5','finding':'SciPy 1.17.0 reflect median may return nonbinary/incorrect values for some sub-20-frame multicolumn arrays and very wide windows. Reproduced for n=2 with widths17/33/65/129. These paths violate the dataset minimum20 rule. Explicit manual reflection matches every declared median width at actual lengths71/250, length249 and anticipated T4 length5001. No demonstrated effect on saved benchmark results.'},
                    {'severity':'reporting_note','id':'M6','finding':'An independent cosine solver produces small threshold metric differences at values within about1e-14 of .5; raw disagreements are recorded. This is numerical tie sensitivity rather than an observed Brier or algorithm mismatch. Production cubic-versus-linear disagreement totals remain as saved.'}],
        'independent_checks':ev}
    out=ROOT/'reviews/01_math_initial.json'; out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'report':str(out),'resolvent_max_error':dense_error,'median_max_error':median_error,'villin_sample_metric_max_error':sampled_error,'sampled_method_rows':len(sampled_checks),'verdict':result['verdict']},indent=2))

if __name__=='__main__': run()
