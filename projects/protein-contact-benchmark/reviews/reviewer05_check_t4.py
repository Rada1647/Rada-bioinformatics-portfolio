"""Independent final T4 audit: no imports from the benchmark implementation."""
from pathlib import Path
import collections,datetime,hashlib,json,math,platform,subprocess,sys
import numpy as np
from scipy.fft import dct,idct
R=Path(__file__).resolve().parents[1];P=R/'data/prepared/t4'
checks=collections.Counter();cache={}
def read(p):return json.loads(Path(p).read_text())
def check(ok,label):
 checks[label]+=1
 if not bool(ok):raise AssertionError(label)
def sha(p):
 p=Path(p).resolve()
 if p not in cache:
  h=hashlib.sha256()
  with p.open('rb') as f:
   for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
  cache[p]=h.hexdigest()
 return cache[p]
def stamp(s):return datetime.datetime.fromisoformat(s)
def progress(s):print(s,flush=True)
report={'reviewer':'05 independent numerical reproduction AI audit','checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'running','critical_findings':[]}
protocol=read(R/'protocol/protocol_v1.1.json');manifest=read(P/'manifest.json');psha=sha(R/'protocol/protocol_v1.1.json');msha=sha(P/'manifest.json')
check(manifest['protocol_sha256']==psha,'manifest protocol hash');check(manifest['code_sha256']==sha(R/'prepare_extension.py'),'preparation code hash')
check(psha==(R/'protocol/protocol_v1.1.sha256').read_text().split()[0],'frozen hash sidecar')
check(stamp(protocol['created_utc'])<stamp(manifest['prepared_utc']),'protocol before T4 preparation')
check([(r['id'],r['run_index'],r['split']) for r in manifest['runs']]==[('md1us4',0,'discovery'),('md1us5',1,'validation'),('md1us6',2,'test'),('md1us7',3,'test'),('md1us8',4,'test')],'frozen T4 trajectory split')
input_hashes={};coordinates={};source_checks=[]
for run in manifest['runs']:
 cp=P/run['coordinates_path'];check(sha(cp)==run['coordinates_sha256'],'prepared CA hash')
 check(sha(R/'data/t4'/(run['id']+'_ca_200ps.npz'))==sha(cp),'source/derived CA hash')
 extraction=read(R/'data/t4'/(run['id']+'_extraction.json'));source=read(R/'data/t4'/(run['id']+'.xtc.checksum.json'))
 check(extraction['output_sha256']==sha(cp),'extraction output hash')
 check(source['published_checksum']=='md5:'+source['md5'] and source['published_size']==source['verified_size'],'T4 recorded original checksum and size')
 source_checks.append({'run':run['id'],'source_xtc_sha256':source['sha256'],'source_xtc_md5':source['md5'],'derived_ca_sha256':sha(cp),'scope':'Verified retained source-checksum/extraction records against prepared arrays; this reviewer does not redownload original T4 source XTCs.'})
 input_hashes[run['id']]={'coordinates':sha(cp),'contacts':{}}
 for cutoff,p in run['contacts'].items():
  check(sha(P/p)==run['contacts_sha256'][cutoff],'prepared contact hash');input_hashes[run['id']]['contacts'][cutoff]=sha(P/p)
 with np.load(cp) as d:coordinates[run['id']]=(d['xyz'],d['time_ps'])
progress('T4 prepared arrays, source records, split and hashes verified.')
# Independent reconstruction: squared distances, blocks along frames.
allpairs=np.array([(i,j) for i in range(162) for j in range(i+4,162)],dtype=np.int64)
def contacts(xyz,pairs,cutoff):
 out=np.empty((len(xyz),len(pairs)),bool)
 for start in range(0,len(xyz),96):
  frame=np.asarray(xyz[start:start+96],dtype=float);delta=frame[:,pairs[:,0],:]-frame[:,pairs[:,1],:]
  out[start:start+96]=np.einsum('ijk,ijk->ij',delta,delta)<float(cutoff)**2
 return out
coverage=read(P/'coverage.json');reconstruction=[]
for cutoff in ['7','8','9']:
 discovery=contacts(coordinates['md1us4'][0],allpairs,cutoff);sums=discovery.sum(axis=0);total=len(discovery);keep=(sums>0)&(sums<total);freq=sums[keep]/total
 for run in manifest['runs']:
  recon=discovery if run['id']=='md1us4' else contacts(coordinates[run['id']][0],allpairs,cutoff)
  with np.load(P/run['contacts'][cutoff]) as d:
   check(np.array_equal(d['X'],recon[:,keep]),'T4 independently reconstructed contact array')
   check(np.array_equal(d['pairs'],allpairs[keep]),'T4 discovery feature selection')
   check(np.array_equal(d['frequency'],freq),'T4 discovery frequency')
   check(np.array_equal(d['time_ps'],coordinates[run['id']][1]),'T4 frame-time identity')
  cr=coverage[cutoff]['runs'][run['id']]
  check(cr['all_positive_cells']==int(recon.sum()),'T4 coverage all positives')
  check(cr['selected_positive_cells']==int(recon[:,keep].sum()),'T4 coverage selected positives')
  check(cr['discovery_absent_positive_cells']==int(recon[:,sums==0].sum()),'T4 coverage absent positives')
  if run['id']!='md1us4':del recon
 reconstruction.append({'cutoff':int(cutoff),'all_pairs':len(allpairs),'selected_pairs':int(keep.sum()),'discovery_frames':total,'contact_arrays_reconstructed':5})
 del discovery
progress('All 15 T4 contact arrays, discovery frequencies and coverage independently reconstructed.')
# Independent comparison is root-path restricted and checks scalar types/finiteness too.
allowed={'completed_utc','seconds','method_prediction_and_metrics_seconds','validation_file_sha256'}
excluded=collections.Counter();numeric=other=0;maxabs=0.
def compare(a,b,path):
 global numeric,other,maxabs
 check(type(a)==type(b),'paired field types')
 if isinstance(a,dict):
  check(set(a)==set(b),'paired keys')
  for k in a:
   if k in allowed:check(path.count('/')==0,'excluded metadata only at root');excluded[k]+=1
   else:compare(a[k],b[k],path+'/'+k)
 elif isinstance(a,list):
  check(len(a)==len(b),'paired lengths')
  for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+'/'+str(i))
 elif type(a) in (float,int):
  numeric+=1;check(math.isfinite(a) and math.isfinite(b),'numeric finiteness');maxabs=max(maxabs,abs(a-b));check(a==b,'exact numerical reproduction')
 else:other+=1;check(a==b,'exact scalar reproduction')
names=[s['name']+'_'+phase+'.json' for s in protocol['scenarios'] for phase in ['validation','test']]+['primary_clean.json']
summary=read(R/'reviews/fresh_t4_comparison.json');summaryfiles={x['file']:x for x in summary['files']};check(set(names)==set(summaryfiles),'full 17 expected output names')
check(set(summary['excluded_metadata_keys'])==allowed,'default comparator exclusions retained')
outputs={};chronology=[];seed_counts={};environments={};order=['raw','linear','hard_linear','moving_average','gaussian','median','power_full','power_matched']
for dirname in ['t4','reproduced_t4']:
 base=R/'results'/dirname;seeds=set();invocations=sorted(base.glob('invocation_*.json'));check(len(invocations)==8,'eight scenario invocations');environments[dirname]=[]
 for ip in invocations:
  inv=read(ip);environments[dirname].append(inv['environment']);check(stamp(protocol['created_utc'])<stamp(inv['started_utc']),'protocol predates invocation')
  for p,h in inv['inputs_and_code_sha256'].items():check(sha(p)==h,'invocation actual input/code hash')
 for si,s in enumerate(protocol['scenarios']):
  val=read(base/(s['name']+'_validation.json'));test=read(base/(s['name']+'_test.json'))
  check(test['validation_file_sha256']==sha(base/(s['name']+'_validation.json')),'within-run validation pointer hash')
  check(test['selected_parameters']==val['selected_parameters'],'saved validation choices used for test')
  inferred=stamp(test['completed_utc'])-datetime.timedelta(seconds=test['seconds']);check(stamp(val['completed_utc'])<=inferred,'selection before inferred test start')
  check((base/(s['name']+'_validation.json')).stat().st_mtime<=(base/(s['name']+'_test.json')).stat().st_mtime,'selection file before test file')
  chronology.append({'output':dirname,'scenario':s['name'],'validation_completed_utc':val['completed_utc'],'test_start_inferred_utc':inferred.isoformat(),'gap_seconds':(inferred-stamp(val['completed_utc'])).total_seconds()})
  for phase,kind,obj in [(1,'validation',val),(2,'test',test)]:
   outputs[(dirname,s['name'],kind)]=obj
   check(obj['protocol_sha256']==psha and obj['prepared_manifest_sha256']==msha,'output protocol/manifest hash')
   check(obj['prepared_inputs_sha256']==input_hashes,'output input hashes')
   check(obj['code_sha256']=={p:sha(R/p) for p in ['benchmark_extension.py','methods.py']},'output code hashes')
   check(obj['scenario']==s and obj['scenario_index']==si,'output scenario identity')
   records=obj['records'] if kind=='validation' else [dict(x,run_id=r['run_id'],run_index=r['run_index']) for r in obj['per_run'] for x in r['noise_rows']]
   expected={(r['id'],r['run_index'],n) for r in manifest['runs'] if r['split']==kind for n in range(8 if phase==1 else 16)}
   check({(x['run_id'],x['run_index'],x['noise_index']) for x in records}==expected,'full seed row coverage')
   check(len(records)==len(expected),'full seed row count')
   for row in records:
    words=(20261007,phase,1,row['run_index'],si,row['noise_index']);check(tuple(row['seed_words'])==words,'frozen seed derivation');check(words not in seeds,'no seed-key duplication');seeds.add(words)
  for method,grid in val['candidate_grids'].items():
   means=np.array([row['candidate_brier'][method] for row in val['records']]).mean(axis=0)
   check(np.max(np.abs(means-np.asarray(val['mean_candidate_brier'][method])))<2e-15,'validation candidate means')
   check(val['selected_parameters'][method]==grid[int(np.argmin(means))],'validation candidate argmin')
   check(abs(val['selected_validation_brier'][method]-float(means.min()))<2e-15,'validation minimum')
  check(val['selected_baseline']==min(order,key=lambda k:val['selected_validation_brier'][k]),'full baseline argmin')
  check(val['selected_matched_baseline']==min([k for k in order if k!='power_full'],key=lambda k:val['selected_validation_brier'][k]),'matched baseline argmin')
 seed_counts[dirname]=len(seeds)
 for name in names:
  check(sha(base/name)==summaryfiles[name]['original_sha256' if dirname=='t4' else 'reproduction_sha256'],'comparison file hashes')
for name in names:compare(read(R/'results/t4'/name),read(R/'results/reproduced_t4'/name),name)
check(numeric==summary['numeric_fields_compared'] and other==summary['other_scalar_fields_compared'],'comparison scalar counts')
progress(f'Full default comparison independently verified: {numeric} numeric, {other} other scalars, max difference {maxabs}.')
# Independent spectral solver using analytic path eigenvalues and orthonormal DCT.
# A eigenvalues are (1+cos(pi*k/n))/2; inverse resolvent weights follow directly.
def predict_spectral(z,settings):
 z=np.asarray(z,float);n=len(z);coeff=dct(z,type=2,axis=0,norm='ortho');lambdal=np.sin(np.pi*np.arange(n)/(2*n))**2
 solved={};preds={};residual=0.
 for method in ['raw','linear','paper332','power_full','power_matched']:
  if method=='raw':preds[method]=z;continue
  t=settings[method]['t']
  if t not in solved:
   h=idct(coeff/(1+t*lambdal[:,None]),type=2,axis=0,norm='ortho')
   rh=(1+t*.5)*h.copy();rh[1:]-=t*.25*h[:-1];rh[:-1]-=t*.25*h[1:];rh[0]-=t*.25*h[0];rh[-1]-=t*.25*h[-1]
   residual=max(residual,float(np.max(np.abs(rh-z))));del rh
   h[:,z.sum(axis=0)==0]=0;h[:,z.sum(axis=0)==n]=1;solved[t]=np.clip(h,0,1)
  h=solved[t]
  if method=='linear':preds[method]=h
  elif method=='paper332':preds[method]=3*h*h-2*h*h*h
  else:
   alpha=settings[method]['alpha'];num=h**alpha;preds[method]=num/(num+(1-h)**alpha)
 return preds,residual
def predict_thomas(z,settings):
 """Independent scalar-coefficient forward elimination/back substitution."""
 z=np.asarray(z,float);n=len(z);cache={};predictions={}
 for method in ['raw','linear','paper332','power_full','power_matched']:
  if method=='raw':predictions[method]=z;continue
  t=settings[method]['t']
  if t not in cache:
   if n==1:h=z.copy()
   else:
    diagonal=np.full(n,1+t/2.);diagonal[0]=diagonal[-1]=1+t/4.;off=-t/4.;h=z.copy()
    for i in range(1,n):
     factor=off/diagonal[i-1];diagonal[i]-=factor*off;h[i]-=factor*h[i-1]
    h[-1]/=diagonal[-1]
    for i in range(n-2,-1,-1):h[i]=(h[i]-off*h[i+1])/diagonal[i]
   h[:,z.sum(axis=0)==0]=0;h[:,z.sum(axis=0)==n]=1;cache[t]=np.clip(h,0,1)
  h=cache[t]
  if method=='linear':predictions[method]=h
  elif method=='paper332':predictions[method]=3*h*h-2*h*h*h
  else:
   alpha=settings[method]['alpha'];num=h**alpha;predictions[method]=num/(num+(1-h)**alpha)
 return predictions

maskcache={}
def masks_for(run,x,freq):
 if run in maskcache:return maskcache[run]
 rare=freq<=.1;short=np.zeros_like(x);transition=np.zeros_like(x);events=0;n=len(x)
 for col in range(x.shape[1]):
  i=0
  while i<n:
   if not x[i,col]:i+=1;continue
   start=i
   while i<n and x[i,col]:i+=1
   if start>0 and i<n and i-start<=3:short[start:i,col]=True;events+=1
 for i in range(1,n):
  change=x[i]!=x[i-1];transition[i,change]=True;transition[i-1,change]=True
 rarepos=x&rare
 counts={'rare_features':int(rare.sum()),'rare_positive':int(rarepos.sum()),'transition':int(transition.sum()),'short_positive':int(short.sum()),'brief_event_count':events}
 maskcache[run]=(rare,rarepos,short,transition,counts);return maskcache[run]
def metrics(pred,x,masks):
 rare,rarepos,short,transition,counts=masks;binary=pred>=.5;err=pred-x;sq=err**2;n=len(x)
 tp=int((binary&x).sum());fp=int((binary&~x).sum());fn=int((~binary&x).sum())
 def masked(a,m):return float(a[m].sum()/m.sum()) if m.any() else None
 return {'brier':float(sq.sum()/sq.size),'mae':float(np.abs(err).sum()/err.size),'occupancy_mae':float(np.abs((pred.sum(axis=0)-x.sum(axis=0))/n).mean()),'occupancy_signed_bias':float(err.sum()/err.size),'precision':tp/(tp+fp) if tp+fp else None,'recall':tp/(tp+fn) if tp+fn else None,'f1':2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None,'rare_brier':float(sq[:,rare].mean()) if rare.any() else None,'rare_positive_recall':masked(binary,rarepos),'short_positive_recall':masked(binary,short),'transition_brier':masked(sq,transition)}
def corrupt(x,xyz,pairs,s,words):
 rng=np.random.default_rng(np.random.SeedSequence(words));eps=s['level']
 if s['kind']=='bitflip':return x!=(rng.random(x.shape)<eps)
 if s['kind']=='jitter':return contacts(np.asarray(xyz,float)+rng.normal(scale=eps,size=xyz.shape),pairs,s['cutoff_A'])
 e=np.empty_like(x);e[0]=rng.random(x.shape[1])<eps
 for i in range(1,len(x)):
  p=np.full(x.shape[1],eps*(1-s['rho']));p[e[i-1]]=eps+(1-eps)*s['rho'];e[i]=rng.random(x.shape[1])<p
 return x!=e
boundary_findings=[];max_solver_probability_difference=0.;max_spectral_metric_difference=0.
def compare_algorithm_variants(spectral,thomas,x,masks,target,context):
 global max_solver_probability_difference,max_spectral_metric_difference
 classified={'precision','recall','f1','rare_positive_recall','short_positive_recall'}
 for method,pred in spectral.items():
  difference=float(np.max(np.abs(pred-thomas[method])));max_solver_probability_difference=max(max_solver_probability_difference,difference)
  check(difference<2e-12,'independent spectral/Thomas probability agreement')
  m=metrics(pred,x,masks);mdiff={}
  for metric,value in m.items():
   wanted=target[method][metric]
   if value is not None:
    delta=value-wanted;max_spectral_metric_difference=max(max_spectral_metric_difference,abs(delta))
    if abs(delta)>2e-12:
     check(metric in classified,'spectral differences restricted to classification metrics');mdiff[metric]={'spectral':value,'recorded':wanted,'signed_difference':delta}
  mismatch=(pred>=.5)!=(thomas[method]>=.5);count=int(mismatch.sum())
  if count:
   proximity=float(np.max(np.abs(pred[mismatch]-.5)));check(proximity<2e-12,'classification differences occur only near threshold')
   cells=np.argwhere(mismatch)
   boundary_findings.append({**context,'method':method,'classification_disagreements_spectral_vs_Thomas':count,'max_probability_difference':difference,'maximum_distance_to_threshold_at_disagreement':proximity,'spectral_metric_differences_vs_recorded':mdiff,'example_cells':[{'frame_index':int(i),'feature_index':int(j),'spectral':float(pred[i,j]),'Thomas':float(thomas[method][i,j])} for i,j in cells[:20]]})
  check(not mdiff or count>0,'spectral metric differences accounted for by threshold disagreements')

runmap={r['id']:r for r in manifest['runs']};spot=[];metric_n=0;maxdelta=0.;maxresidual=0.
def check_metrics(predictions,x,masks,target):
 global metric_n,maxdelta
 localdelta=0.
 for method,pred in predictions.items():
  actual=metrics(pred,x,masks)
  for metric,value in actual.items():
   wanted=target[method][metric];metric_n+=1;check((value is None)==(wanted is None),'independent metric null identity')
   if value is not None:
    diff=abs(value-wanted);maxdelta=max(maxdelta,diff);localdelta=max(localdelta,diff)
    check(diff<2e-12,'independent spectral metric tolerance')
 return localdelta
for scenario_name in ['flip10_cut8','jitter050_cut8','correlated10_cut8']:
 test=outputs[('t4',scenario_name,'test')];val=outputs[('t4',scenario_name,'validation')]
 for idx in ([0,1,2] if scenario_name=='flip10_cut8' else [0]):
  rr=test['per_run'][idx];run=runmap[rr['run_id']]
  with np.load(P/run['contacts']['8']) as d:x=d['X'];pairs=d['pairs'];freq=d['frequency']
  masks=masks_for(run['id'],x,freq);check(masks[-1]==rr['stratum_counts'],'independent T4 stratum counts')
  for ni in [0,15]:
   row=rr['noise_rows'][ni];z=corrupt(x,coordinates[run['id']][0],pairs,test['scenario'],row['seed_words']);spectral,res=predict_spectral(z,test['selected_parameters']);predictions=predict_thomas(z,test['selected_parameters']);compare_algorithm_variants(spectral,predictions,x,masks,row['methods'],{'scenario':scenario_name,'run_id':run['id'],'noise_index':ni});del spectral;maxresidual=max(maxresidual,res)
   delta=check_metrics(predictions,x,masks,row['methods']);spot.append({'scenario':scenario_name,'run_id':run['id'],'noise_index':ni,'shape':list(x.shape),'max_absolute_metric_difference':delta,'max_spectral_solve_residual':res});del z,predictions
  progress('Independent spectral metrics '+scenario_name+' '+run['id']+' verified.')
 vr=val['records'][0];run=runmap[vr['run_id']]
 with np.load(P/run['contacts']['8']) as d:x=d['X'];pairs=d['pairs']
 z=corrupt(x,coordinates[run['id']][0],pairs,val['scenario'],vr['seed_words'])
 for t in [1,2,4,8,16,32,64]:
  settings={'raw':{},'linear':{'t':t},'paper332':{'t':t},'power_full':{'t':t,'alpha':3},'power_matched':{'t':t,'alpha':3}}
  predictions,res=predict_spectral(z,settings);maxresidual=max(maxresidual,res)
  for method in ['linear','paper332']:
   index=val['candidate_grids'][method].index({'t':t});diff=abs(float(np.mean((predictions[method]-x)**2))-vr['candidate_brier'][method][index]);metric_n+=1;maxdelta=max(maxdelta,diff);check(diff<2e-12,'independent validation score tolerance')
  del predictions
 del z
clean=read(R/'results/t4/primary_clean.json')
for rr in clean['per_run']:
 run=runmap[rr['run_id']]
 with np.load(P/run['contacts']['8']) as d:x=d['X'];freq=d['frequency']
 masks=masks_for(run['id'],x,freq);check(masks[-1]==rr['stratum_counts'],'clean T4 strata')
 spectral,res=predict_spectral(x,clean['selected_parameters']);predictions=predict_thomas(x,clean['selected_parameters']);compare_algorithm_variants(spectral,predictions,x,masks,rr['methods'],{'scenario':'primary_clean','run_id':run['id']});del spectral;maxresidual=max(maxresidual,res);check_metrics(predictions,x,masks,rr['methods']);del predictions
progress(f'Independent T4 numerical reconstruction: {metric_n} scores/metrics; maximum metric difference {maxdelta:g}; solver residual {maxresidual:g}.')
packages=json.loads(subprocess.check_output([sys.executable,'-m','pip','list','--format=json'],text=True))
check({x['name']:x['version'] for x in packages}=={'numpy':'2.3.5','pip':'25.0.1','scipy':'1.17.0'},'fresh benchmark dependency package set')
check('include-system-site-packages = false' in Path(sys.prefix,'pyvenv.cfg').read_text(),'fresh venv excludes system packages')
report.update({'status':'pass_with_scope_limits','protocol_sha256':psha,'prepared_manifest_sha256':msha,'source_record_checks':source_checks,'all_contacts_reconstruction':reconstruction,'full_comparison':{'files':17,'numeric_fields':numeric,'other_scalar_fields':other,'max_absolute_difference':maxabs,'excluded_root_key_occurrences':dict(excluded),'comparator_mode':'default, without allow-regenerated-manifest'},'seeds':seed_counts,'chronology':chronology,'environments':environments,'independent_numerical_reconstruction':{'methods':['raw','linear','paper332','power_full','power_matched'],'test_cases':spot,'validation_candidate_scores':42,'clean_runs':3,'metric_values_total':metric_n,'absolute_tolerance':2e-12,'max_absolute_metric_difference':maxdelta,'max_spectral_solver_residual':maxresidual,'implementation':'Independent Thomas forward elimination/back substitution supplies reconstructed test/clean metrics; analytic path eigenvalues sin(pi*k/(2*n))^2 and orthonormal DCT-II independently cross-check continuous predictions and validation candidate scores. Separate corruption, squared-distance contacts, event scanner and metric formulas. Imports no benchmark code.','case_selection':'All 3 primary test runs, first and last noise draws; first test run for jitter050 and correlated10 first and last draws; first validation row for 7 linear and cubic candidates; all 3 clean test runs.'},'fresh_environment':{'python':sys.version,'executable':sys.executable,'venv_config':Path(sys.prefix,'pyvenv.cfg').read_text(),'packages':packages},'readme_followup':'Confirmed explicit .venv interpreter/dependency-install commands and final raw md1us8 retention disclosure now present.','limitations':['Fresh installation/replay covers benchmark NumPy/SciPy dependencies only, on the same host and Python interpreter family. No independent-host or freshly installed raw MD extraction-stack claim.','T4 source XTC-to-CA extraction is reviewed through records and arrays here, not independently redownloaded/re-extracted by reviewer05. A separate source-level reviewer may supply that evidence.','Local timestamps support chronology but do not establish public preregistration or tamper-proof proof of no previous outcomes; test-start timestamps are inferred from completion minus elapsed runtime.','Full numerical equality is same-implementation replay; independent Thomas and spectral reconstruction sample selected methods/cases. Spectral classification is rounding-sensitive near 0.5 and is not claimed to be bitwise identical.','AI technical audit is not human peer review or experimental biological validation.'],'checks_total':sum(checks.values()),'check_counts':dict(checks),'script_sha256':sha(__file__),'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
primary=outputs[('t4','flip10_cut8','test')];baseline=primary['selected_baseline'];bg=primary['aggregate'][baseline]['brier']-primary['aggregate']['paper332']['brier']
report['finite_precision_threshold_audit']={'finding':'Independent DCT spectral probabilities agree with Thomas to roundoff, but a few near-0.5 calls change due to floating-point evaluation. Thomas metrics reproduce recorded metrics at the declared 2e-12 tolerance; this does not establish classification invariance across solvers or platforms.','max_spectral_vs_Thomas_probability_difference':max_solver_probability_difference,'max_raw_spectral_metric_difference_vs_recorded':max_spectral_metric_difference,'affected_method_cases':boundary_findings,'initial_failed_strict_spectral_check_log':'reviews/05_reproduction_final_machinecheck_attempt1.log','primary_decision_sensitivity':{'recorded_status':primary['primary_success']['status'],'baseline':baseline,'aggregate_brier_gain':bg,'required_gain':0.02*primary['aggregate'][baseline]['brier'],'Brier_condition_remains_failed':bg<0,'explanation':'Primary composite requires Brier improvement. Its observed aggregate Brier gain is negative by orders of magnitude more than solver roundoff. Threshold rounding cannot make this continuous-metric condition pass; therefore the primary failure is unchanged regardless of these classification-boundary cells.'}}
villin=read(R/'reviews/05_reproduction_initial.json')
optional=read(R/'reviews/05_optional_comparator_review.json')
report['villin_prior_completed_review']={'review_file':'reviews/05_reproduction_initial.json','sha256':sha(R/'reviews/05_reproduction_initial.json'),'status':villin['status'],'full_comparison':villin['exact_full_comparison'],'independent_contact_arrays':villin['data_audit']['all_contact_arrays_reconstructed'],'independent_metric_values':villin['independent_numerical_reconstruction']['metric_and_candidate_values_compared'],'max_independent_metric_difference':villin['independent_numerical_reconstruction']['max_absolute_metric_difference'],'machinecheck_log_sha256':sha(R/'reviews/05_reproduction_machinecheck.log')}
report['optional_manifest_comparator_review']={'file':'reviews/05_optional_comparator_review.json','sha256':sha(R/'reviews/05_optional_comparator_review.json'),'status':optional['status'],'tested_cases':len(optional['tests']),'used_for_fresh_T4_comparison':False}
report['audited_code_sha256']={name:sha(R/name) for name in ['benchmark_extension.py','methods.py','prepare_extension.py','compare_reproduction.py','run_scenarios.py']}
report['spectral_preflight']={'file':'reviews/05_spectral_audit_preflight.json','sha256':sha(R/'reviews/05_spectral_audit_preflight.json'),'cases':28,'maximum_absolute_difference':read(R/'reviews/05_spectral_audit_preflight.json')['max_absolute_difference']}
report['overall_scope']='Full same-host fresh-benchmark-dependency replay of all eight scenarios and clean controls for both T4 and villin; independent contact reconstruction for all prepared arrays and independent algorithm/metric spot checks for both datasets.'
(R/'reviews/05_reproduction_final.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
progress('PASS: reviews/05_reproduction_final.json written.')
