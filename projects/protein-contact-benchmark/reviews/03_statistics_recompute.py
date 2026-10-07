"""Independent reviewer aggregation audit. Does not import benchmark code."""
from pathlib import Path
import collections, datetime, hashlib, json, math, os, zipfile
import numpy as np
from scipy.stats import t

R = Path(__file__).resolve().parents[1]
P = json.loads((R/'protocol/protocol_v1.1.json').read_text())
C = collections.Counter()
FAIL = []
MAX = collections.defaultdict(float)

def load(p): return json.loads((R/p).read_text())
def sha(p): return hashlib.sha256((R/p).read_bytes()).hexdigest()
def eq(a,b,label,category='numeric'):
    if isinstance(a,dict) and isinstance(b,dict):
        if set(a)!=set(b): FAIL.append({'field':label,'key_mismatch':[sorted(a),sorted(b)]})
        for k in a.keys() & b.keys(): eq(a[k],b[k],label+'.'+str(k),category)
    elif isinstance(a,list) and isinstance(b,list):
        if len(a)!=len(b): FAIL.append({'field':label,'length_mismatch':[len(a),len(b)]})
        for i,(x,y) in enumerate(zip(a,b)): eq(x,y,label+f'[{i}]',category)
    elif isinstance(a,(int,float)) and not isinstance(a,bool) and isinstance(b,(int,float)) and not isinstance(b,bool):
        C[category+'_numeric_fields_checked']+=1
        delta=abs(a-b);MAX[category]=max(MAX[category],delta)
        if not math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12): FAIL.append({'field':label,'expected':a,'actual':b,'absolute_difference':delta})
    else:
        C[category+'_nonnumeric_fields_checked']+=1
        if a!=b: FAIL.append({'field':label,'expected':a,'actual':b})

def mean(xs):
    v=[x for x in xs if x is not None]
    return math.fsum(v)/len(v) if v else None
def meanrows(rows): return {k:mean([r[k] for r in rows]) for k in rows[0]}
def counts(rows): return {k:sum(r[k] is not None for r in rows) for k in rows[0]}
def finite_range(xs):
    values=[x for x in xs if x is not None]
    return [min(values),max(values)] if values else None
def independent_strata(X,f):
    rare=X[:,f<=0.1]
    trans=np.zeros(X.shape,bool)
    trans[:-1] |= X[:-1]!=X[1:]
    trans[1:] |= X[:-1]!=X[1:]
    brief=np.zeros_like(X);events=0
    for col in range(X.shape[1]):
        i=0
        while i<len(X):
            if not X[i,col]: i+=1;continue
            start=i
            while i<len(X) and X[i,col]:i+=1
            if start>0 and i<len(X) and i-start<=3:
                brief[start:i,col]=True;events+=1
    return {'rare_features':int((f<=0.1).sum()),'rare_positive':int(rare.sum()),'transition':int(trans.sum()),'short_positive':int(brief.sum()),'brief_event_count':events}

def composite(runs,agg,baseline):
    rel=(agg[baseline]['brier']-agg['paper332']['brier'])/agg[baseline]['brier'] if agg[baseline]['brier'] else None
    rr=[]
    for run in runs:
        b,p=run['aggregate'][baseline],run['aggregate']['paper332']
        rare=None if b['rare_positive_recall'] is None or p['rare_positive_recall'] is None else b['rare_positive_recall']-p['rare_positive_recall']
        brief=None if b['short_positive_recall'] is None or p['short_positive_recall'] is None else b['short_positive_recall']-p['short_positive_recall']
        rr.append({'run_id':run['run_id'],'brier_gain':b['brier']-p['brier'],'positive_brier_gain':b['brier']>p['brier'],'rare_recall_loss':rare,'brief_recall_loss':brief,'rare_guardrail_pass':None if rare is None else rare<=0.02,'brief_guardrail_pass':None if brief is None else brief<=0.02})
    conditions={'relative_gain_at_least_minimum':None if rel is None else rel>=0.02,'positive_gain_every_test_run':all(r['positive_brier_gain'] for r in rr),'rare_guardrail_every_test_run':None if any(r['rare_guardrail_pass'] is None for r in rr) else all(r['rare_guardrail_pass'] for r in rr),'brief_guardrail_every_test_run':None if any(r['brief_guardrail_pass'] is None for r in rr) else all(r['brief_guardrail_pass'] for r in rr)}
    status='not_evaluable' if None in conditions.values() else ('pass' if all(conditions.values()) else 'fail')
    return {'baseline':baseline,'relative_mean_brier_gain':rel,'per_run':rr,'conditions':conditions,'status':status}

def audit(dataset):
    M=load(f'data/prepared/{dataset}/manifest.json'); midx={r['id']:r for r in M['runs']}
    D=R/f'data/prepared/{dataset}'
    summary=load('summary.json')['datasets'].get(dataset)
    S=[];allintervals=[];allt=[];allv=[];stratacache={}
    order=['raw','linear','hard_linear','moving_average','gaussian','median','power_full','power_matched']
    for si,s in enumerate(P['scenarios']):
        n=s['name'];v=load(f'results/{dataset}/{n}_validation.json');z=load(f'results/{dataset}/{n}_test.json')
        grid={k:g for k,g in P['methods'].items() if k!='oracle_debiased_linear' or s['kind'] in ['bitflip','correlated']}
        eq(grid,v['candidate_grids'],n+'.grids','candidate_grid')
        eq({k:len(g) for k,g in grid.items()},v['candidate_counts'],n+'.candidate_counts','candidate_count')
        expected_ids=[r['id'] for r in M['runs'] if r['split']=='validation']
        eq(collections.Counter({i:8 for i in expected_ids}),dict(collections.Counter(r['run_id'] for r in v['records'])),n+'.validation_run_counts','run_counts')
        for rid in expected_ids:
            eq(list(range(8)),[r['noise_index'] for r in v['records'] if r['run_id']==rid],n+'.'+rid+'.validation_noise_indices','validation_noise_count')
        vmeans={};sel={};best={};ties={}
        for method,g in grid.items():
            vals=[[r['candidate_brier'][method][j] for r in v['records']] for j in range(len(g))]
            C['validation_candidate_scores_read']+=sum(map(len,vals))
            vmeans[method]=[mean(x) for x in vals]
            eq(vmeans[method],v['mean_candidate_brier'][method],n+'.means.'+method,'validation_mean')
            # Tie rule checked using archived arithmetic as exact ties may be rounding-sensitive.
            rawmeans=v['mean_candidate_brier'][method]
            idx=min(range(len(g)),key=lambda j:rawmeans[j]);sel[method]=g[idx]
            best[method]=rawmeans[idx];ties[method]=[j for j,x in enumerate(rawmeans) if x==rawmeans[idx]]
            independent_idx=min(range(len(g)),key=lambda j:vmeans[method][j])
            eq(idx,independent_idx,n+'.stable_selection.'+method,'selection_stability')
        eq(sel,v['selected_parameters'],n+'.selected','selected_parameters')
        eq(best,v['selected_validation_brier'],n+'.best','selected_validation')
        baseline=min(order,key=best.get);matched=min([m for m in order if m!='power_full'],key=best.get)
        eq(baseline,v['selected_baseline'],n+'.baseline','baseline_choice');eq(matched,v['selected_matched_baseline'],n+'.matched','baseline_choice')
        eq(sel,z['selected_parameters'],n+'.test_settings','locked_test_settings')
        eq(baseline,z['selected_baseline'],n+'.test_baseline','locked_test_settings');eq(matched,z['selected_matched_baseline'],n+'.test_matched','locked_test_settings')
        eq(sha(f'results/{dataset}/{n}_validation.json'),z['validation_file_sha256'],n+'.validation_hash','identity')
        for result in [v,z]:
            eq(sha('protocol/protocol_v1.1.json'),result['protocol_sha256'],n+'.protocol_hash','identity')
            eq(sha(f'data/prepared/{dataset}/manifest.json'),result['prepared_manifest_sha256'],n+'.manifest_hash','identity')
            eq(s,result['scenario'],n+'.scenario','identity')
        eq(True,P['created_utc']<v['completed_utc']<z['completed_utc'],n+'.freeze_validation_test_order','chronology')
        for row in v['records']:
            r=midx[row['run_id']];eq([20261007,1,M['dataset_id'],r['run_index'],si,row['noise_index']],row['seed_words'],n+'.validation_seed','seed')
            allv.append(tuple(row['seed_words']))
        runaggs=[];recomputed_runs=[];nulls=collections.Counter()
        eq([r['id'] for r in M['runs'] if r['split']=='test'],[r['run_id'] for r in z['per_run']],n+'.test_run_order','run_counts')
        for run in z['per_run']:
            rid=run['run_id'];rows=run['noise_rows'];eq(list(range(16)),[r['noise_index'] for r in rows],n+'.noise_indices','noise_count')
            for row in rows:
                eq([20261007,2,M['dataset_id'],midx[rid]['run_index'],si,row['noise_index']],row['seed_words'],n+'.test_seed','seed')
                allt.append(tuple(row['seed_words']))
            ag={method:meanrows([r['methods'][method] for r in rows]) for method in grid}
            eq(ag,run['aggregate'],n+'.'+rid+'.aggregate','per_run_aggregate')
            eq({method:counts([r['methods'][method] for r in rows]) for method in grid},run['aggregate_contributing_noise_count'],n+'.'+rid+'.contributing_noise','noise_contributor')
            eq(meanrows([r['same_t_linear'] for r in rows]),run['same_t_linear_aggregate'],n+'.'+rid+'.same_t_aggregate','same_t_per_run_aggregate')
            eq(counts([r['same_t_linear'] for r in rows]),run['same_t_linear_contributing_noise_count'],n+'.'+rid+'.same_t_contributing','same_t_noise_contributor')
            eq(sum(r['same_t_threshold_disagreements'] for r in rows),run['same_t_threshold_disagreements'],n+'.'+rid+'.threshold','threshold')
            cachekey=(rid,int(s['cutoff_A']))
            if cachekey not in stratacache:
                with np.load(D/midx[rid]['contacts'][str(s['cutoff_A'])]) as a: stratacache[cachekey]=independent_strata(a['X'],a['frequency'])
            eq(stratacache[cachekey],run['stratum_counts'],n+'.'+rid+'.strata','stratum_denominator')
            for method in grid:
                for metric,val in ag[method].items():
                    if val is None:nulls[metric]+=1
            for method in grid:
                if method=='paper332':continue
                a=[r['methods'][method]['brier'] for r in rows];b=[r['methods']['paper332']['brier'] for r in rows]
                dif=[x-y for x,y in zip(a,b)];dm=mean(dif);var=math.fsum((x-dm)**2 for x in dif)/(len(dif)-1)
                half=float(t.ppf(.975,len(dif)-1))*math.sqrt(var/len(dif))
                pi={'baseline':method,'n_noise_realizations':len(dif),'baseline_mean_brier':mean(a),'paper332_mean_brier':mean(b),'mean_brier_gain':dm,'relative_mean_brier_gain':dm/mean(a) if mean(a)>0 else None,'noise_only_paired_95pct_ci':[dm-half,dm+half],'positive_noise_differences':sum(x>0 for x in dif),'negative_noise_differences':sum(x<0 for x in dif)}
                eq(pi,{k:run['paired_comparisons'][method][k] for k in pi},n+'.'+rid+'.CI.'+method,'paired_interval')
                C['intervals_checked']+=1
                allintervals.append({'scenario':n,'run_id':rid,**pi})
            runaggs.append(ag);recomputed_runs.append({'run_id':rid,'aggregate':ag})
        agg={method:meanrows([a[method] for a in runaggs]) for method in grid}
        eq(agg,z['aggregate'],n+'.aggregate','overall_aggregate')
        eq({method:counts([a[method] for a in runaggs]) for method in grid},z['aggregate_contributing_run_count'],n+'.contributing_runs','run_contributor')
        eq(meanrows([r['same_t_linear_aggregate'] for r in z['per_run']]),z['same_t_linear_aggregate'],n+'.same_t_aggregate','same_t_overall')
        eq(counts([r['same_t_linear_aggregate'] for r in z['per_run']]),z['same_t_linear_contributing_run_count'],n+'.same_t_contributors','same_t_overall_count')
        eq(sum(r['same_t_threshold_disagreements'] for r in z['per_run']),z['same_t_threshold_disagreements'],n+'.threshold_total','threshold')
        c=agg['paper332']['brier'];bb=agg[baseline]['brier'];mb=agg[matched]['brier'];lin=agg['linear']['brier']
        item={'scenario':n,'validation_record_count':len(v['records']),'selected_parameters':sel,'candidate_minimum_exact_tie_indices':ties,'baseline':baseline,'baseline_tied_methods':[m for m in order if best[m]==best[baseline]],'matched_baseline':matched,'paper332_brier':c,'selected_baseline_brier':bb,'matched_baseline_brier':mb,'linear_brier':lin,'gain_vs_selected_baseline_percent':100*(bb-c)/bb,'gain_vs_matched_percent':100*(mb-c)/mb,'gain_vs_linear_percent':100*(lin-c)/lin,'null_run_method_metric_counts':dict(nulls),'test_runs':len(runaggs),'contributing_run_counts_unique':sorted(set(x for m in z['aggregate_contributing_run_count'].values() for x in m.values())),'contributing_noise_counts_unique':sorted(set(x for r in z['per_run'] for m in r['aggregate_contributing_noise_count'].values() for x in m.values()))}
        if s['primary']:
            for stored,base in [('stress_success_descriptive' if dataset=='villin' else 'primary_success',baseline),('matched_budget_success_secondary',matched)]:
                comp=composite(recomputed_runs,agg,base);eq(comp,{k:z[stored][k] for k in comp},n+'.'+stored,'composite')
                item[stored]={k:v for k,v in comp.items() if k!='per_run'}
                item[stored]['run_counts']={'positive_brier':sum(r['positive_brier_gain'] for r in comp['per_run']),'rare_guardrail_pass':sum(r['rare_guardrail_pass'] is True for r in comp['per_run']),'brief_guardrail_pass':sum(r['brief_guardrail_pass'] is True for r in comp['per_run'])}
                item[stored]['rare_loss_range']=finite_range(r['rare_recall_loss'] for r in comp['per_run'])
                item[stored]['brief_loss_range']=finite_range(r['brief_recall_loss'] for r in comp['per_run'])
                item[stored]['nonevaluable_rare_runs']=[r['run_id'] for r in comp['per_run'] if r['rare_guardrail_pass'] is None]
                item[stored]['nonevaluable_brief_runs']=[r['run_id'] for r in comp['per_run'] if r['brief_guardrail_pass'] is None]
            item['weighting_evidence']={'test_frames_distribution':dict(collections.Counter(r['shape'][0] for r in z['per_run'])),'equal_run_cubic_brier':c,'cell_weighted_cubic_brier_for_diagnostic_only':math.fsum(r['aggregate']['paper332']['brier']*math.prod(r['shape']) for r in z['per_run'])/sum(math.prod(r['shape']) for r in z['per_run']),'used_for_score':'equal run means; cell weighted value is diagnostic only'}
        if summary:
            ss=next(x for x in summary['scenarios'] if x['name']==n)
            expected={'name':n,'baseline':baseline,'matched_baseline':matched,'cubic_brier':c,'baseline_brier':bb,'linear_brier':lin,'gain_vs_baseline_percent':100*(bb-c)/bb,'gain_vs_linear_percent':100*(lin-c)/lin,'selected_parameters':sel,'aggregate':agg,'source_file_sha256':sha(f'results/{dataset}/{n}_test.json'),'same_t_threshold_disagreements':z['same_t_threshold_disagreements']}
            eq(expected,ss,n+'.summary','summary')
            if s['primary']:
                eq(agg,summary['primary']['aggregate'],n+'.primary_summary_aggregate','summary')
                eq(z['aggregate_contributing_run_count'],summary['primary']['contributing_runs'],n+'.summary_contributors','summary')
                eq(sel,summary['primary']['selected_parameters'],n+'.summary_selected_parameters','summary')
                eq(baseline,summary['primary']['baseline'],n+'.summary_baseline','summary')
                eq(matched,summary['primary']['matched_baseline'],n+'.summary_matched_baseline','summary')
                eq(z['stress_success_descriptive' if dataset=='villin' else 'primary_success'],summary['primary']['criterion'],n+'.summary_criterion','summary')
                eq(z['matched_budget_success_secondary'],summary['primary']['matched_criterion'],n+'.summary_matched_criterion','summary')
                for sr,tr in zip(summary['primary']['per_run'],z['per_run']):eq(sr,{k:tr[k] for k in sr},n+'.summary_per_run','summary')
        S.append(item)
    eq(len(allv),len(set(allv)),'validation_unique_seeds','seed_uniqueness');eq(len(allt),len(set(allt)),'test_unique_seeds','seed_uniqueness')
    eq(0,len(set(allv)&set(allt)),'validation_test_disjoint_seeds','seed_uniqueness')
    clean=load(f'results/{dataset}/primary_clean.json')
    eq(load(f'results/{dataset}/flip10_cut8_validation.json')['selected_parameters'],clean['selected_parameters'],'clean.primary_settings','locked_clean_settings')
    ca={method:meanrows([r['methods'][method] for r in clean['per_run']]) for method in clean['selected_parameters']}
    eq(ca,clean['aggregate'],'clean.aggregate','clean_aggregate')
    eq({method:counts([r['methods'][method] for r in clean['per_run']]) for method in clean['selected_parameters']},clean['aggregate_contributing_run_count'],'clean.contributors','clean_contributor')
    if summary:eq(ca,summary['clean'],'clean.summary','summary')
    return {'scenarios':S,'validation_noise_records':len(allv),'test_noise_records':len(allt),'clean_aggregate_brier':{m:a['brier'] for m,a in ca.items()}}

def audit_villin_selection_coverage():
    M=load('data/prepared/villin/manifest.json');TM=load('data/villin/trajectory_manifest.json')
    D=R/'data/prepared/villin';cov=load('data/prepared/villin/coverage.json');S=load('summary.json')['datasets']['villin']
    original=sorted(TM,key=lambda r:r['source_path']);fullidx={r['trajectory_id']:i for i,r in enumerate(original)}
    chosen=[];groups={};splits={'1':('discovery',24),'2':('validation',24),'3':('test',48)}
    for group,(split,cap) in splits.items():
        eligible=[r for r in original if str(r['group'])==group and r['eligibility_min20_frames'] and r['geometry_pass'] and r['n_ca']==35 and abs(r['dt_ps']-200)<.001]
        target=sorted(eligible,key=lambda r:hashlib.sha256(r['source_path'].encode()).hexdigest())[:cap]
        actual=[r for r in M['runs'] if r['split']==split]
        eq(sorted(r['trajectory_id'] for r in target),sorted(r['id'] for r in actual),'hash_caps.'+group,'hash_cap')
        for r in actual:
            eq(fullidx[r['id']],r['run_index'],'run_index.'+r['id'],'run_index');eq(group,str(r['group']),'split_group.'+r['id'],'split_group')
            eq(hashlib.sha256(r['source_path'].encode()).hexdigest(),r['path_selection_sha256'],'path_hash.'+r['id'],'path_hash')
        groups[group]={'split':split,'eligible':len(eligible),'selected':len(target),'unique_ancestral_roots_selected':len(set(r['ancestry_root_or_missing_parent'] for r in target))}
    xyz={r['id']:np.load(D/r['coordinates_path'])['xyz'].astype(np.float64) for r in M['runs']}
    pairs=np.array([(i,j) for i in range(35) for j in range(i+4,35)],dtype=np.int64)
    distances={rid:np.sqrt(((x[:,pairs[:,0]]-x[:,pairs[:,1]])**2).sum(axis=2)) for rid,x in xyz.items()}
    coverage_summary={}
    for cut in ['7','8','9']:
        allX={rid:d<float(cut) for rid,d in distances.items()};disc=np.concatenate([allX[r['id']] for r in M['runs'] if r['split']=='discovery'])
        sums=disc.sum(axis=0);keep=(sums>0)&(sums<len(disc));f=sums[keep]/len(disc)
        exp={'all_eligible_pair_count':len(pairs),'selected_pair_count':int(keep.sum()),'discovery_frame_count':len(disc),'discovery_constant_absent_pairs':int((sums==0).sum()),'discovery_constant_present_pairs':int((sums==len(disc)).sum())}
        eq(exp,{k:cov[cut][k] for k in exp},'coverage.'+cut,'coverage')
        for r in M['runs']:
            rid=r['id'];X=allX[rid];total=int(X.sum());absent=int(X[:,sums==0].sum());selected=int(X[:,keep].sum())
            rc={'split':r['split'],'all_positive_cells':total,'selected_positive_cells':selected,'discovery_absent_positive_cells':absent,'discovery_absent_positive_fraction':absent/total if total else None,'omitted_positive_cells':total-selected}
            eq(rc,cov[cut]['runs'][rid],'coverage.'+cut+'.'+rid,'coverage')
            with np.load(D/r['contacts'][cut]) as saved:
                eq(True,bool(np.array_equal(saved['X'],X[:,keep])),'contacts_from_coordinates.'+rid+'.'+cut,'feature_integrity')
                eq(True,bool(np.array_equal(saved['pairs'],pairs[keep])),'discovery_pairs.'+rid+'.'+cut,'feature_integrity')
                eq(True,bool(np.array_equal(saved['frequency'],f)),'discovery_frequency.'+rid+'.'+cut,'feature_integrity')
                C['contact_cells_recomputed_from_coordinates']+=saved['X'].size
        tests=[cov[cut]['runs'][r['id']] for r in M['runs'] if r['split']=='test']
        total=sum(r['all_positive_cells'] for r in tests);absent=sum(r['discovery_absent_positive_cells'] for r in tests)
        coverage_summary[cut]={**exp,'test_total_positive_cells':total,'test_discovery_absent_positive_cells':absent,'pooled_test_absent_positive_fraction':absent/total}
        if cut=='8':eq(absent/total,S['primary_omitted_discovery_absent_positive_fraction_pooled'],'summary.coverage','summary')
    return {'hash_selection':groups,'coverage':coverage_summary}

def audit_t4_selection_coverage():
    """Bounded-memory coordinate audit suitable for all T4 frame/pair cells."""
    M=load('data/prepared/t4/manifest.json');D=R/'data/prepared/t4';cov=load('data/prepared/t4/coverage.json')
    summary=load('summary.json')['datasets'].get('t4')
    expected=[{'id':'md1us'+str(i),'run_index':i-4,'split':'discovery' if i==4 else 'validation' if i==5 else 'test'} for i in range(4,9)]
    eq(expected,[{k:r[k] for k in expected[0]} for r in M['runs']],'t4.frozen_run_splits','t4_split')
    nca=M['runs'][0]['n_ca'];pairs=np.array([(i,j) for i in range(nca) for j in range(i+4,nca)],dtype=np.int64)
    def getxyz(r):
        with np.load(D/r['coordinates_path']) as data:
            xyz=data['xyz'].astype(np.float64);times=data['time_ps']
        eq(True,bool(np.isfinite(xyz).all()),'t4.finite.'+r['id'],'coordinate_time')
        eq(True,bool(np.allclose(np.diff(times),200,atol=.001,rtol=0)),'t4.cadence.'+r['id'],'coordinate_time')
        eq([len(xyz),nca,3],list(xyz.shape),'t4.shape.'+r['id'],'coordinate_time')
        return xyz
    sums={cut:np.zeros(len(pairs),np.int64) for cut in ['7','8','9']};frames=0
    for r in M['runs']:
        if r['split']!='discovery':continue
        xyz=getxyz(r);frames+=len(xyz)
        for start in range(0,len(xyz),64):
            block=xyz[start:start+64]
            dist=np.sqrt(((block[:,pairs[:,0]]-block[:,pairs[:,1]])**2).sum(axis=2))
            for cut in sums:sums[cut]+=(dist<float(cut)).sum(axis=0)
    keep={cut:(ss>0)&(ss<frames) for cut,ss in sums.items()}
    result={}
    for cut in sums:
        fields={'all_eligible_pair_count':len(pairs),'selected_pair_count':int(keep[cut].sum()),'discovery_frame_count':frames,'discovery_constant_absent_pairs':int((sums[cut]==0).sum()),'discovery_constant_present_pairs':int((sums[cut]==frames).sum())}
        eq(fields,{k:cov[cut][k] for k in fields},'t4.coverage.'+cut,'coverage')
        result[cut]=fields
        if summary:eq(fields['selected_pair_count'],summary['selected_features'][cut],'t4.summary_feature_count.'+cut,'summary')
    for r in M['runs']:
        xyz=getxyz(r);record={cut:{'all':0,'absent':0,'selected':0} for cut in sums};saved={}
        for cut in sums:
            with np.load(D/r['contacts'][cut]) as z:
                saved[cut]=z['X']
                eq(True,bool(np.array_equal(z['pairs'],pairs[keep[cut]])),'t4.pairs.'+r['id']+'.'+cut,'feature_integrity')
                eq(True,bool(np.array_equal(z['frequency'],sums[cut][keep[cut]]/frames)),'t4.frequency.'+r['id']+'.'+cut,'feature_integrity')
                C['contact_cells_recomputed_from_coordinates']+=z['X'].size
        for start in range(0,len(xyz),64):
            block=xyz[start:start+64];dist=np.sqrt(((block[:,pairs[:,0]]-block[:,pairs[:,1]])**2).sum(axis=2))
            for cut in sums:
                X=dist<float(cut);rec=record[cut]
                rec['all']+=int(X.sum());rec['absent']+=int(X[:,sums[cut]==0].sum());rec['selected']+=int(X[:,keep[cut]].sum())
                eq(True,bool(np.array_equal(saved[cut][start:start+64],X[:,keep[cut]])),'t4.contact_cells.'+r['id']+'.'+cut+f'.{start}','feature_integrity')
                C['all_eligible_reference_contact_cells_recomputed']+=X.size
        for cut in sums:
            rr=record[cut];total=rr['all'];absent=rr['absent'];selected=rr['selected']
            fields={'split':r['split'],'all_positive_cells':total,'selected_positive_cells':selected,'discovery_absent_positive_cells':absent,'discovery_absent_positive_fraction':absent/total if total else None,'omitted_positive_cells':total-selected}
            eq(fields,cov[cut]['runs'][r['id']],'t4.coverage.'+cut+'.'+r['id'],'coverage')
    for cut in sums:
        tests=[cov[cut]['runs'][r['id']] for r in M['runs'] if r['split']=='test']
        total=sum(r['all_positive_cells'] for r in tests);absent=sum(r['discovery_absent_positive_cells'] for r in tests)
        result[cut].update(test_total_positive_cells=total,test_discovery_absent_positive_cells=absent,pooled_test_absent_positive_fraction=absent/total if total else None)
        if summary and cut=='8':
            eq(absent/total if total else None,summary['primary_omitted_discovery_absent_positive_fraction_pooled'],'t4.summary_coverage','summary')
            eq({r['id']:cov[cut]['runs'][r['id']] for r in M['runs'] if r['split']=='test'},summary['primary_coverage_per_run'],'t4.summary_coverage_per_run','summary')
    return {'frozen_split_verified':expected,'coverage':result}

def freeze_receipt_evidence():
    receipt=load('protocol/library_save_receipts.json')
    observed=[]
    for row in receipt['results']:
        path=R/'protocol'/row['file_name']
        attrs={k:os.getxattr(path,k).decode() for k in os.listxattr(path) if k.startswith('user.library-')}
        eq(row['local_xattrs'],attrs,'library_xattrs.'+row['file_name'],'receipt_metadata')
        observed.append({'file_name':row['file_name'],'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'local_xattrs':attrs})
    with zipfile.ZipFile(R/'protocol/protein332_frozen_protocol.zip') as archive:
        h=hashlib.sha256(archive.read('protocol_v1.1.json')).hexdigest()
        eq(sha('protocol/protocol_v1.1.json'),h,'protocol_zip_content','frozen_zip')
    invocations=[load(str(x.relative_to(R))) for x in (R/'results/villin').glob('invocation_*.json')]
    first=min(x['started_utc'] for x in invocations)
    eq(True,P['created_utc']<first,'freeze_before_first_invocation','chronology')
    return {'receipt_copy':receipt,'observed_files':observed,'first_villin_invocation_utc':first,'protocol_created_utc':P['created_utc'],'interpretation':'Local archive bytes and persistent-save identifier xattrs agree. The receipt copy explicitly was written during analysis and lacks an independent server timestamp; original tool ordering is reported by parent-session provenance, not independently proved by local receipt. Not a public preregistration.'}

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--dataset',default='villin');ap.add_argument('--out',default='reviews/03_statistics_initial.json');args=ap.parse_args()
    if args.dataset=='all':
        a={'villin':audit('villin'),'t4':audit('t4')}
        a['villin'].update(audit_villin_selection_coverage());a['t4'].update(audit_t4_selection_coverage())
    else:
        a=audit(args.dataset)
        if args.dataset=='villin':a.update(audit_villin_selection_coverage())
        if args.dataset=='t4':a.update(audit_t4_selection_coverage())
    a['freeze_receipt_evidence']=freeze_receipt_evidence()
    report={'reviewer':'Independent AI technical audit 03: statistics and protocol','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':args.dataset+' completed statistics; T4 confirmation pending' if args.dataset=='villin' else args.dataset,'independence':'Recomputed from archived candidate scores and noise rows with standalone reviewer code; no import of production analysis modules. Recomputed discovery features, reference contacts, coverage, and strata from coordinates/contacts. Does not independently regenerate method predictions, which belongs to reproduction audit.','status':'pass_with_interpretation_limits' if not FAIL else 'issues_found','protocol_sha256':sha('protocol/protocol_v1.1.json'),'code_sha256':sha('benchmark_extension.py'),'machine_checks':dict(C),'maximum_absolute_numeric_differences':dict(MAX),'failures':FAIL,'evidence':a,'interpretation_limits':['All intervals use 16 paired artificial-noise replicates per fixed trajectory, df=15. Frames, contacts, trajectories or seeds are not asserted to be independent biological replicates.','Per-run 95% intervals are pointwise and unadjusted across scenarios, runs and comparators. They support descriptive conditional noise statements, not familywise significance or biological generalization.','Villin 48 test paths belong to one held-out adaptive archive group; within-group ancestry and adaptive selection imply dependence. Held-out groups have unverified campaign independence. Hash caps prevent outcome selection but do not restore independence or equilibrium sampling.','Validation-selected full power has 63 candidates while cubic and matched alpha=3 family have seven; both full-budget primary and matched secondary comparisons must be reported.','Prospective freeze is locally timestamped and hash-identified; an external Library receipt was not yet available to this reviewer. No claim of public preregistration is justified.'],'critical_findings':[]}
    report['interpretation_limits'][-1]='Frozen archive, current JSON hash and version-0 Library identifier xattrs agree. Receipt copy was created during analysis and has no independent server timestamp; original save-before-run ordering is supported by parent conversation provenance. This is a private prospective freeze, not a public registry preregistration.'
    report['summary_numeric_field_count']=sum(v for k,v in C.items() if k.endswith('_numeric_fields_checked'))
    report['findings']=[{'severity':'none','finding':'No discrepancy in candidate aggregation, tuning/tie-breaking, equal-run metric weighting, contributor/null handling, conditional CIs, effect arithmetic, composite decision, summary, hash-capped selection, feature selection or coverage.'},{'severity':'interpretation','finding':'Villin primary cubic has higher Brier than full/matched power in all 48 held-out paths. Its 35.895% gain over linear alone does not establish comparative success. One secondary jitter scenario has only 0.02194% gain over selected baseline and cannot replace the primary scenario.'},{'severity':'interpretation','finding':'The 48 selected test trajectories descend from only six named roots within one adaptive archive group. Even roots may be coupled by adaptive selection; report these as dependent stress paths, not 48 independent biological simulations.'}]
    (R/args.out).write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    ss=a['scenarios'] if args.dataset!='all' else a['t4']['scenarios']+a['villin']['scenarios']
    print(json.dumps({'status':report['status'],'failures':len(FAIL),'machine_checks':dict(C),'scenario_effects':[{k:s[k] for k in ['scenario','paper332_brier','selected_baseline_brier','gain_vs_selected_baseline_percent','gain_vs_linear_percent']} for s in ss]},indent=2))
