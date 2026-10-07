"""Independent T4 mathematical review. No production/results files are modified.

--self-check validates the efficient spectral method while outcomes are pending.
The default final audit requires every prespecified T4 result file to exist.
"""
from pathlib import Path
import argparse, datetime, hashlib, importlib.util, json, sys
import numpy as np
from scipy.fft import dct, idct

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import methods

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def cosine_resolvent(z,t):
    n=len(z)
    transformed=dct(z.astype(float),type=2,axis=0,norm='ortho')
    transformed/=(1+t*np.sin(np.pi*np.arange(n)/(2*n))**2)[:,None]
    h=idct(transformed,type=2,axis=0,norm='ortho')
    h=np.clip(h,0,1)
    h[:,np.all(z==0,axis=0)]=0
    h[:,np.all(z==1,axis=0)]=1
    return h

def residual(z,h,t):
    v=h-z.astype(float)
    d=h[:-1]-h[1:]
    v[:-1]+=(t/4)*d
    v[1:]-=(t/4)*d
    return float(np.max(abs(v)))

def selfcheck():
    rng=np.random.default_rng(401772)
    dense=0.; residual_max=0.; difference=0.
    for n in (1,2,3,7,19,71):
        z=rng.integers(0,2,(n,13))
        a=np.eye(n)
        if n>1:
            for i in range(n-1):
                a[i,i]-=.25;a[i+1,i+1]-=.25
                a[i,i+1]=a[i+1,i]=.25
        for t in (1,2,4,8,16,32,64):
            h=cosine_resolvent(z,t)
            expected=np.linalg.solve(np.eye(n)+t*(np.eye(n)-a),z)
            dense=max(dense,float(np.max(abs(h-expected))))
            residual_max=max(residual_max,residual(z,h,t))
    z=rng.integers(0,2,(5001,32))
    for t in (1,2,4,8,16,32,64):
        h=cosine_resolvent(z,t)
        difference=max(difference,float(np.max(abs(h-methods.predict(z,'linear',{'t':t})))))
        residual_max=max(residual_max,residual(z,h,t))
    assert dense<1e-13 and residual_max<1e-12 and difference<1e-13
    import benchmark_extension as implementation
    protocol=json.loads((ROOT/'protocol/protocol_v1.1.json').read_text())
    coords=rng.normal(size=(131,10,3)).astype(np.float32)
    ii,jj=np.triu_indices(10,4);pairs=np.column_stack((ii,jj))
    x=rng.integers(0,2,(131,len(pairs))).astype(bool)
    columns=np.array([0,3,7,11,14,20])
    for index,scenario in enumerate(protocol['scenarios']):
        words=[20261007,2,1,3,index,0]
        z,score=stream_noise(x,coords,pairs,scenario,words,columns)
        full=implementation.corrupt(x,coords,pairs,scenario,words)
        assert np.array_equal(z,full[:,columns])
        assert score==float(np.mean(full!=x))
    return {'max_small_dense_error':dense,'max_tridiagonal_residual':residual_max,
            'max_5001_frame_spectral_vs_banded_error':difference,
            'length_5001_columns':32,'t_grid':[1,2,4,8,16,32,64],
            'streamed_noise_exact_matches_all_8_scenarios':True,'streamed_noise_test_shape':[131,len(pairs)]}

def stream_noise(x,coords,pairs,s,words,columns):
    """Preserve full-width random stream while storing only selected channels.

    Also recompute the full feature raw Brier numerator, independently of
    saved summaries and without retaining a full floating-point noise array.
    """
    rng=np.random.Generator(np.random.PCG64(np.random.SeedSequence(words)))
    n,c=x.shape; out=np.empty((n,len(columns)),bool); wrong=0
    epsilon=s['level']; previous=None
    for start in range(0,n,64):
        stop=min(n,start+64); truth=x[start:stop]
        if s['kind']=='bitflip':
            corruption=rng.uniform(size=(stop-start,c))<epsilon
            z=truth^corruption
        elif s['kind']=='correlated':
            uniform=rng.uniform(size=(stop-start,c)); corruption=np.empty_like(truth)
            for i in range(stop-start):
                if previous is None: previous=uniform[i]<epsilon
                else:
                    p=epsilon*(1-s['rho'])+previous.astype(float)*s['rho']
                    previous=uniform[i]<p
                corruption[i]=previous
            z=truth^corruption
        else:
            noisy=coords[start:stop].astype(float)+rng.standard_normal(size=coords[start:stop].shape)*epsilon
            z=np.empty_like(truth)
            for f in range(0,c,128):
                p=pairs[f:f+128]
                delta=noisy[:,p[:,0]]-noisy[:,p[:,1]]
                z[:,f:f+128]=np.sum(delta*delta,axis=2)<s['cutoff_A']**2
        wrong+=int(np.count_nonzero(z!=truth))
        out[start:stop]=z[:,columns]
    return out,wrong/(n*c)

def audit_stored(v,r,scenario_index,dataset_id,expected_valruns,expected_testruns):
    aggregate_error=0.; selection_error=0.; baseline_order=['raw','linear','hard_linear','moving_average','gaussian','median','power_full','power_matched']
    assert len(v['records'])==expected_valruns*8
    assert len(r['per_run'])==expected_testruns
    scores={}
    for name,grid in v['candidate_grids'].items():
        means=np.mean([row['candidate_brier'][name] for row in v['records']],axis=0)
        selection_error=max(selection_error,float(np.max(abs(means-v['mean_candidate_brier'][name]))))
        idx=min(range(len(grid)),key=lambda i:means[i])
        assert v['selected_parameters'][name]==grid[idx]
        assert r['selected_parameters'][name]==grid[idx]
        scores[name]=means[idx]
    assert min(baseline_order,key=lambda k:scores[k])==v['selected_baseline']==r['selected_baseline']
    assert min([k for k in baseline_order if k!='power_full'],key=lambda k:scores[k])==v['selected_matched_baseline']==r['selected_matched_baseline']
    for row in v['records']:
        assert row['seed_words']==[20261007,1,dataset_id,row['run_index'],scenario_index,row['noise_index']]
    for run in r['per_run']:
        assert len(run['noise_rows'])==16
        for row in run['noise_rows']:
            assert row['seed_words']==[20261007,2,dataset_id,run['run_index'],scenario_index,row['noise_index']]
            if row['same_t_threshold_disagreements']==0:
                for metric in ('precision','recall','f1','rare_positive_recall','short_positive_recall'):
                    assert row['methods']['paper332'][metric]==row['same_t_linear'][metric]
        for name in r['selected_parameters']:
            for key,target in run['aggregate'][name].items():
                vals=[row['methods'][name][key] for row in run['noise_rows'] if row['methods'][name][key] is not None]
                assert len(vals)==run['aggregate_contributing_noise_count'][name][key]
                if vals: aggregate_error=max(aggregate_error,abs(float(np.mean(vals))-target))
                else: assert target is None
    for name in r['selected_parameters']:
        for key,target in r['aggregate'][name].items():
            vals=[row['aggregate'][name][key] for row in r['per_run'] if row['aggregate'][name][key] is not None]
            assert len(vals)==r['aggregate_contributing_run_count'][name][key]
            if vals: aggregate_error=max(aggregate_error,abs(float(np.mean(vals))-target))
            else: assert target is None
    assert r['same_t_threshold_disagreements']==sum(run['same_t_threshold_disagreements'] for run in r['per_run'])
    for run in r['per_run']:
        assert run['same_t_threshold_disagreements']==sum(row['same_t_threshold_disagreements'] for row in run['noise_rows'])
    assert aggregate_error<1e-14 and selection_error<1e-14
    return {'max_candidate_selection_mean_error':selection_error,'max_aggregate_error':aggregate_error}

def final_audit():
    pre=selfcheck()
    protocol_path=ROOT/'protocol/protocol_v1.1.json'; protocol=json.loads(protocol_path.read_text())
    mp=ROOT/'data/prepared/t4/manifest.json'; manifest=json.loads(mp.read_text()); runs={r['id']:r for r in manifest['runs']}
    assert [r['id'] for r in manifest['runs']]==['md1us4','md1us5','md1us6','md1us7','md1us8']
    assert all(r['n_frames']==5001 and r['n_ca']==162 for r in manifest['runs'])
    assert manifest['protocol_sha256']==sha(protocol_path)
    scenarios=[]; inspected_hashes={}; numerical=[]; raw_error=0.
    for si,s in enumerate(protocol['scenarios']):
        paths=[ROOT/'results/t4'/f"{s['name']}_{phase}.json" for phase in ('validation','test')]
        v,r=[json.loads(p.read_text()) for p in paths]
        inspected_hashes.update({str(p.relative_to(ROOT)):sha(p) for p in paths})
        for doc in (v,r):
            assert doc['protocol_sha256']==sha(protocol_path) and doc['prepared_manifest_sha256']==sha(mp)
            for fn,h in doc['code_sha256'].items(): assert sha(ROOT/fn)==h
        assert r['validation_file_sha256']==sha(paths[0])
        agg=audit_stored(v,r,si,1,1,3)
        for runresult in r['per_run']:
            run=runs[runresult['run_id']]
            with np.load(mp.parent/run['contacts'][str(s['cutoff_A'])]) as data:
                x=data['X'];pairs=data['pairs']
            coords=None
            if s['kind']=='jitter':
                with np.load(mp.parent/run['coordinates_path']) as data: coords=data['xyz']
            columns=np.unique(np.linspace(0,x.shape[1]-1,min(32,x.shape[1]),dtype=int))
            # One full-width seeded draw per run/scenario; spectral check on 32 columns.
            z,raw_brier=stream_noise(x,coords,pairs,s,[20261007,2,1,run['run_index'],si,0],columns)
            delta=abs(raw_brier-runresult['noise_rows'][0]['methods']['raw']['brier'])
            raw_error=max(raw_error,delta); assert delta<1e-15
            params=r['selected_parameters']['paper332'];t=params['t']
            h=cosine_resolvent(z,t)
            band=methods.predict(z,'linear',{'t':t})
            y=3*h*h-2*h*h*h
            production=methods.predict(z,'paper332',params)
            error=float(np.max(abs(h-band))); cubic_error=float(np.max(abs(y-production)))
            res=residual(z,h,t)
            assert error<1e-13 and cubic_error<2e-13 and res<1e-12
            disagree=(y>=.5)!=(production>=.5)
            threshold_error=float(np.max(abs(y[disagree]-.5))) if np.any(disagree) else 0.
            assert threshold_error<1e-12
            same_t=int(np.count_nonzero((band>=.5)!=(production>=.5)))
            powers=[]
            for method in ('power_full','power_matched'):
                p=r['selected_parameters'][method];ph=cosine_resolvent(z,p['t']);alpha=p['alpha']
                expected=ph**alpha/(ph**alpha+(1-ph)**alpha)
                observed=methods.predict(z,method,p)
                err=float(np.max(abs(expected-observed)))
                # alpha<1 amplifies harmless endpoint solver roundoff.
                assert err<1e-7
                powers.append({'method':method,'parameters':p,'max_probability_error':err})
            numerical.append({'scenario':s['name'],'run_id':run['id'],'noise_index':0,
                'shape':list(x.shape),'contact_column_indices':columns.tolist(),'cubic_selected_t':t,
                'independent_full_raw_brier':raw_brier,'raw_brier_saved_difference':delta,
                'max_resolvent_error':error,'max_cubic_probability_error':cubic_error,
                'max_tridiagonal_residual':res,'spectral_vs_production_cubic_threshold_disagreements':int(disagree.sum()),
                'threshold_disagreement_max_distance_from_half':threshold_error,
                'production_same_t_linear_cubic_disagreements_on_sample':same_t,'power_comparisons':powers})
        scenarios.append({'scenario':s['name'],'selected_baseline':r['selected_baseline'],
            'selected_matched_baseline':r['selected_matched_baseline'],'selected_parameters':r['selected_parameters'],
            'paper332_brier':r['aggregate']['paper332']['brier'],'baseline_brier':r['aggregate'][r['selected_baseline']]['brier'],
            'same_t_linear_brier':r['same_t_linear_aggregate']['brier'],
            'same_t_threshold_disagreements':r['same_t_threshold_disagreements'],'aggregation_check':agg,
            'primary_success':r.get('primary_success'),'matched_budget_success_secondary':r.get('matched_budget_success_secondary')})
        print('reviewed T4 '+s['name'],flush=True)
    clean_path=ROOT/'results/t4/primary_clean.json';clean=json.loads(clean_path.read_text())
    inspected_hashes[str(clean_path.relative_to(ROOT))]=sha(clean_path)
    assert clean['selected_parameters']==scenarios[0]['selected_parameters'] and clean['oracle_epsilon']==.1
    clean_error=0.
    for name in clean['selected_parameters']:
        for k,target in clean['aggregate'][name].items():
            vals=[r['methods'][name][k] for r in clean['per_run'] if r['methods'][name][k] is not None]
            assert len(vals)==clean['aggregate_contributing_run_count'][name][k]
            if vals: clean_error=max(clean_error,abs(float(np.mean(vals))-target))
            else: assert target is None
    assert clean_error<1e-14
    initial=json.loads((ROOT/'reviews/01_math_initial.json').read_text())
    primary=json.loads((ROOT/'results/t4/flip10_cut8_test.json').read_text())
    selected=primary['selected_baseline'];b=primary['aggregate'][selected];p=primary['aggregate']['paper332']
    gain=(b['brier']-p['brier'])/b['brier']
    run_evidence=[]
    for run in primary['per_run']:
        rb=run['aggregate'][selected];rp=run['aggregate']['paper332']
        run_evidence.append({'run_id':run['run_id'],'brier_gain':rb['brier']-rp['brier'],
            'rare_recall_loss':rb['rare_positive_recall']-rp['rare_positive_recall'],
            'brief_recall_loss':rb['short_positive_recall']-rp['short_positive_recall']})
    rules=protocol['primary_success']
    decision_conditions={'relative_gain_at_least_minimum':gain>=rules['relative_mean_brier_gain_min'],
        'positive_gain_every_test_run':all(x['brier_gain']>0 for x in run_evidence),
        'rare_guardrail_every_test_run':all(x['rare_recall_loss']<=rules['max_rare_recall_loss_each_run'] for x in run_evidence),
        'brief_guardrail_every_test_run':all(x['brief_recall_loss']<=rules['max_brief_recall_loss_each_run'] for x in run_evidence)}
    decision='pass' if all(decision_conditions.values()) else 'fail'
    assert decision_conditions==primary['primary_success']['conditions'] and decision==primary['primary_success']['status']
    initial_checks=initial['independent_checks']
    initial_evidence={
        'resolvent_and_expected_median':initial_checks['resolvent_and_expected_median'],
        'villin_test_noise_rows_checked':initial_checks['villin_audit']['test_noise_rows_checked'],
        'villin_independently_rebuilt_method_metric_rows':initial_checks['villin_audit']['independently_rebuilt_method_metric_rows'],
        'villin_max_rebuilt_metric_error_after_disclosed_tie_alignment':initial_checks['villin_audit']['max_independent_rebuilt_metric_error'],
        'villin_max_stored_aggregate_error':initial_checks['villin_audit']['max_stored_aggregate_error'],
        'actual_length_manual_median_check_count':len(initial_checks['boundaries_constants_power']['actual_length_manual_median_checks']),
        'villin_scenarios':initial_checks['villin_audit']['scenarios'],
        'scope_and_limitations':'Initial report retains every inspected output hash, exact sampled reconstruction scope, raw threshold discrepancies and unused very-short-path dependency defect.'}
    final={'reviewer':'Independent AI technical audit 1: mathematics and implementation',
        'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'final_complete',
        'production_code_implemented_by_reviewer':False,'initial_review':'reviews/01_math_initial.json',
        'initial_review_sha256':sha(ROOT/'reviews/01_math_initial.json'),
        'reviewed_code_sha256':{str(p):sha(ROOT/p) for p in ('protocol/protocol_v1.1.json','methods.py','benchmark_extension.py','prepare_extension.py')},
        'checker_sha256':sha(__file__),'scope':'Extends initial all-scenario Villin audit with all eight T4 validation/test output selections and aggregates, clean aggregate, and independent full raw Brier plus 32-column mathematical checks for all three T4 test runs at noise0 in every scenario. Not a line-by-line proof review of the cotype manuscript and not an independent full data rerun.',
        'verdict':'No critical discrepancy affecting reviewed benchmark outputs detected. Initial unused very-short-path median dependency defect and floating-point tie caveats remain.',
        'preparation_hardening_review':{'artifact':'protocol/preparation_integrity_hardening.json','sha256':sha(ROOT/'protocol/preparation_integrity_hardening.json'),
            'conclusion':'Diff from preserved initial preparer adds only T4 source record, extraction completion/hash, geometry and metadata integrity checks. No contact definition, discovery selection, split, seed, denoiser, metric or statistical rule changed. Numerical benchmark and protocol hashes unchanged.'},
        'self_check':pre,'scenario_results':scenarios,'independent_t4_checks':numerical,
        'independent_primary_decision':{'selected_baseline':selected,'relative_brier_gain':gain,
            'conditions':decision_conditions,'decision':decision,'per_run':run_evidence},
        'initial_villin_and_mathematical_evidence':initial_evidence,
        'draft_report_claims_review':{'file':'RESULTS.txt','sha256_at_review':sha(ROOT/'RESULTS.txt'),
            'conclusion':'Mathematical formula, binary-only specialization, theorem scope and no physical-coordinate reconstruction claims are appropriate. Recommended explicit T4 same-t attribution and distinguishing independently tuned linear from linear at cubic-selected t.',
            't4_same_t_linear_brier':primary['same_t_linear_aggregate']['brier'],
            't4_same_t_cubic_relative_brier_gain':(primary['same_t_linear_aggregate']['brier']-p['brier'])/primary['same_t_linear_aggregate']['brier'],
            't4_same_t_threshold_disagreements':primary['same_t_threshold_disagreements'],
            't4_same_t_rare_and_brief_recall_identical':all(primary['same_t_linear_aggregate'][k]==p[k] for k in ('rare_positive_recall','short_positive_recall'))},
        'max_full_raw_brier_error':raw_error,'max_clean_aggregate_error':clean_error,'output_sha256':inspected_hashes,
        'retained_findings':initial['findings'],
        'claim_limit':'The binary cubic is the exact coordinatewise expected majority of three independent resolvent endpoints. Metric Markov cotype controls a different geometric inequality; neither it nor these artificial-noise experiments guarantee biological correctness, realizable coordinates or universal Brier improvement.'}
    path=ROOT/'reviews/01_math_final.json';path.write_text(json.dumps(final,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'report':str(path),'verdict':final['verdict'],'raw_brier_max_error':raw_error,
        'max_spectral_resolvent_error':max(x['max_resolvent_error'] for x in numerical)},indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--self-check',action='store_true');args=ap.parse_args()
    if args.self_check: print(json.dumps(selfcheck(),indent=2))
    else: final_audit()
