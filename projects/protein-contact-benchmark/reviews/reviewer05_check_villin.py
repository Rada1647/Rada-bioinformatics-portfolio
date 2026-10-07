"""Independent read-only reproduction audit. Imports no project implementation."""
from pathlib import Path
import collections, datetime, hashlib, json, math, platform, subprocess, sys
import numpy as np

R=Path(__file__).resolve().parents[1]
P=R/'data/prepared/villin'
report={'reviewer':'05 independent numerical reproduction AI audit','checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'running','critical_findings':[], 'limitations':['Same host and same Python 3.12 interpreter family; no independent-host claim.','Fresh installation covers benchmark NumPy/SciPy dependencies only, not the raw MD extraction dependency stack.','Independent numeric reconstruction samples cases and methods; exact full output comparison repeats the production implementation.','Local file timestamps and embedded chronology are audit evidence, not a public preregistration or tamper-proof proof of no prior outcomes.','Raw XTC-to-CA extraction is not independently rerun in this review. T4 is pending.']}
checks=collections.Counter()
def check(condition,label):
 checks[label]+=1
 if not bool(condition):raise AssertionError(label)
def read(p):return json.loads(Path(p).read_text())
cache={}
def sha(p):
 p=Path(p).resolve()
 if p not in cache:
  h=hashlib.sha256()
  with p.open('rb') as f:
   for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
  cache[p]=h.hexdigest()
 return cache[p]
def stamp(v):return datetime.datetime.fromisoformat(v)
def progress(s):print(s,flush=True)
protocol=read(R/'protocol/protocol_v1.1.json'); manifest=read(P/'manifest.json')
psha=sha(R/'protocol/protocol_v1.1.json');msha=sha(P/'manifest.json')
check(psha==(R/'protocol/protocol_v1.1.sha256').read_text().split()[0],'protocol hash sidecar')
check(manifest['protocol_sha256']==psha,'manifest protocol identity')
check(manifest['code_sha256']==sha(R/'prepare_extension.py'),'preparation code hash')
check(protocol['supersedes_sha256']==sha(R/'protocol/protocol_v1.json'),'v1 ancestry hash')
check(stamp(protocol['created_utc'])<stamp(manifest['prepared_utc']),'protocol before input preparation')
report['protocol']={'sha256':psha,'prepared_manifest_sha256':msha,'created_utc':protocol['created_utc'],'prepared_utc':manifest['prepared_utc']}
progress('Protocol/manifest identities verified.')
# Re-hash all prepared inputs once, then use independent maps for invocation/output checks.
input_hashes={}
for run in manifest['runs']:
 cp=P/run['coordinates_path'];check(sha(cp)==run['coordinates_sha256'],'prepared CA hash')
 check(sha(cp)==run['source_metadata']['npz_sha256'],'prepared/source CA identity')
 check(sha(R/'data/villin'/run['source_metadata']['npz_path'])==sha(cp),'source CA copy hash')
 input_hashes[run['id']]={'coordinates':sha(cp),'contacts':{}}
 for cutoff,p in run['contacts'].items():
  check(sha(P/p)==run['contacts_sha256'][cutoff],'prepared contact hash')
  input_hashes[run['id']]['contacts'][cutoff]=sha(P/p)
full=sorted(read(R/'data/villin/trajectory_manifest.json'),key=lambda row:row['source_path'])
expected=[]
for group,split,cap in [('1','discovery',24),('2','validation',24),('3','test',48)]:
 eligible=[(i,row) for i,row in enumerate(full) if str(row['group'])==group and row['eligibility_min20_frames'] and row['geometry_pass'] and row['n_ca']==35 and abs(row['dt_ps']-200)<.001]
 chosen=sorted(eligible,key=lambda ir:hashlib.sha256(ir[1]['source_path'].encode()).hexdigest())[:cap]
 expected.extend((i,row['trajectory_id'],split) for i,row in chosen)
check(sorted(expected)==[(r['run_index'],r['id'],r['split']) for r in manifest['runs']],'stable metadata-only split and run indices')
# Independent distance arithmetic uses squared component sum rather than production norm.
pairs_all=np.array([(i,j) for i in range(35) for j in range(i+4,35)],dtype=np.int64)
def contacts(xyz,pairs,cutoff):
 delta=np.asarray(xyz,dtype=np.float64)[:,pairs[:,0],:]-np.asarray(xyz,dtype=np.float64)[:,pairs[:,1],:]
 return np.einsum('ijk,ijk->ij',delta,delta)<float(cutoff)**2
source_arrays={}
for r in manifest['runs']:
 with np.load(P/r['coordinates_path']) as d:source_arrays[r['id']]=(d['xyz'],d['time_ps'])
coverage=read(P/'coverage.json');counters={}
for cutoff in ['7','8','9']:
 total=sum(len(source_arrays[r['id']][0]) for r in manifest['runs'] if r['split']=='discovery')
 sums=sum(contacts(source_arrays[r['id']][0],pairs_all,cutoff).sum(axis=0) for r in manifest['runs'] if r['split']=='discovery')
 keep=(sums>0)&(sums<total);freq=sums[keep]/total
 for r in manifest['runs']:
  with np.load(P/r['contacts'][cutoff]) as d:
   check(np.array_equal(d['pairs'],pairs_all[keep]),'independent discovery feature pairs')
   check(np.array_equal(d['frequency'],freq),'independent discovery frequency')
   recon=contacts(source_arrays[r['id']][0],pairs_all,cutoff)
   check(np.array_equal(d['X'],recon[:,keep]),'independent CA to contacts')
   check(np.array_equal(d['time_ps'],source_arrays[r['id']][1]),'contact/CA time identity')
   cr=coverage[cutoff]['runs'][r['id']]
   check(cr['all_positive_cells']==int(recon.sum()),'coverage all positives')
   check(cr['selected_positive_cells']==int(recon[:,keep].sum()),'coverage selected positives')
   check(cr['discovery_absent_positive_cells']==int(recon[:,sums==0].sum()),'coverage absent positives')
 counters[cutoff]={'all_pairs':len(pairs_all),'selected_pairs':int(keep.sum()),'discovery_frames':total}
report['data_audit']={'prepared_runs':len(manifest['runs']),'prepared_files':384,'all_contact_arrays_reconstructed':288,'cutoffs':counters,'split_and_seed_run_indices_reconstructed_from_complete_source_manifest':True}
progress('384 prepared hashes and 288 contact arrays independently reconstructed; metadata-only splits verified.')
# Full archive checksums against retained source metadata, without re-extracting MD.
archive=R/'data/villin/protein_folding_datasets.zip';hmd5=hashlib.md5();hsha=hashlib.sha256()
with archive.open('rb') as f:
 for b in iter(lambda:f.read(16*1024**2),b''):hmd5.update(b);hsha.update(b)
q=read(R/'data/villin/qualification.json');fm=read(R/'data/villin/figshare_metadata.json')
filemeta=next(x for x in fm['files'] if x['name']==archive.name)
check(hmd5.hexdigest()==q['archive_md5']==filemeta['computed_md5']==filemeta['supplied_md5'],'full archive MD5')
check(hsha.hexdigest()==q['archive_sha256'],'full archive SHA256')
check(archive.stat().st_size==q['archive_bytes']==filemeta['size'],'full archive size')
report['archive']={'bytes':archive.stat().st_size,'md5':hmd5.hexdigest(),'sha256':hsha.hexdigest(),'source_metadata_scope':'Compared retained Figshare metadata and qualification record, no fresh remote download.'}
progress('Full retained villin source archive SHA256 and published-MD5 metadata match.')
# Path-restricted exclusion check (production comparator excludes keys recursively).
allowed={'completed_utc','seconds','method_prediction_and_metrics_seconds','validation_file_sha256'}
excluded_occurrences=collections.Counter();numeric=other=0;max_abs=0.;excluded_paths=[]
def compare(a,b,path):
 global numeric,other,max_abs
 check(type(a)==type(b),'paired field types')
 if isinstance(a,dict):
  check(set(a)==set(b),'paired object keys')
  for k in a:
   if k in allowed:
    check(path.count('/')==0,'excluded keys exist only at result root')
    excluded_occurrences[k]+=1;excluded_paths.append(path+'/'+k)
   else:compare(a[k],b[k],path+'/'+k)
 elif isinstance(a,list):
  check(len(a)==len(b),'paired array lengths')
  for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+'/'+str(i))
 elif type(a) in (int,float):
  check(math.isfinite(a) and math.isfinite(b),'numeric finiteness');numeric+=1
  max_abs=max(max_abs,abs(a-b));check(a==b,'exact numerical reproduction')
 else:other+=1;check(a==b,'exact other scalar reproduction')
names=[s['name']+'_'+phase+'.json' for s in protocol['scenarios'] for phase in ['validation','test']]+['primary_clean.json']
root_comparison=read(R/'reviews/fresh_villin_comparison.json');root_files={x['file']:x for x in root_comparison['files']}
check(set(root_files)==set(names),'comparator expected filenames')
outputs={};chronology=[];seed_sets={};selection_checks=0
order=['raw','linear','hard_linear','moving_average','gaussian','median','power_full','power_matched']
for dirname in ['villin','reproduced_villin']:
 out=R/'results'/dirname;seed_sets[dirname]=set()
 invocations=sorted(out.glob('invocation_*.json'))
 check(len(invocations)==8,'eight scenario invocations')
 for ip in invocations:
  inv=read(ip);check(stamp(protocol['created_utc'])<stamp(inv['started_utc']),'protocol before benchmark invocation')
  for p,h in inv['inputs_and_code_sha256'].items():check(sha(p)==h,'invocation input and code hash')
 report.setdefault('environments',{})[dirname]=[read(ip)['environment'] for ip in invocations]
 for si,s in enumerate(protocol['scenarios']):
  val=read(out/(s['name']+'_validation.json'));test=read(out/(s['name']+'_test.json'))
  check(test['validation_file_sha256']==sha(out/(s['name']+'_validation.json')),'excluded validation pointer integrity')
  check(val['selected_parameters']==test['selected_parameters'],'frozen validation choices used in test')
  inferred=stamp(test['completed_utc'])-datetime.timedelta(seconds=test['seconds'])
  check(stamp(val['completed_utc'])<=inferred,'selection before inferred test start')
  check((out/(s['name']+'_validation.json')).stat().st_mtime<=(out/(s['name']+'_test.json')).stat().st_mtime,'selection file before test file')
  chronology.append({'output':dirname,'scenario':s['name'],'validation_completed_utc':val['completed_utc'],'test_start_inferred_from_completed_minus_elapsed_utc':inferred.isoformat(),'gap_seconds':(inferred-stamp(val['completed_utc'])).total_seconds()})
  for kind,obj in [('validation',val),('test',test)]:
   outputs[(dirname,s['name'],kind)]=obj
   check(obj['protocol_sha256']==psha and obj['prepared_manifest_sha256']==msha,'output protocol and manifest identity')
   check(obj['code_sha256']=={p:sha(R/p) for p in ['benchmark_extension.py','methods.py']},'output code identity')
   check(obj['prepared_inputs_sha256']==input_hashes,'output prepared file identity')
   check(obj['scenario']==s and obj['scenario_index']==si,'output scenario identity')
   if kind=='validation':
    records=obj['records'];phase=1;count=8
    check(len(records)==24*8,'full validation row count')
   else:
    check(len(obj['per_run'])==48,'full test run count');phase=2;count=16
    records=[dict(x,run_id=r['run_id'],run_index=r['run_index']) for r in obj['per_run'] for x in r['noise_rows']]
   expected_runs=[r for r in manifest['runs'] if r['split']==kind.replace('validation','validation')]
   expected_rows={(r['id'],r['run_index'],n) for r in expected_runs for n in range(count)}
   check({(row['run_id'],row['run_index'],row['noise_index']) for row in records}==expected_rows,'full scenario seed row coverage')
   for row in records:
    words=(20261007,phase,3,row['run_index'],si,row['noise_index'])
    check(tuple(row['seed_words'])==words,'declared seed derivation')
    check(words not in seed_sets[dirname],'unique seed key across dataset phases and scenarios');seed_sets[dirname].add(words)
  for method,grid in val['candidate_grids'].items():
   scores=np.array([row['candidate_brier'][method] for row in val['records']])
   means=scores.mean(axis=0)
   check(np.max(np.abs(means-np.asarray(val['mean_candidate_brier'][method])))<2e-15,'candidate means from saved validation rows within summation roundoff')
   check(val['selected_parameters'][method]==grid[int(np.argmin(means))],'candidate argmin tie order')
   check(abs(val['selected_validation_brier'][method]-float(means.min()))<2e-15,'candidate best Brier within summation roundoff')
   selection_checks+=1
  check(val['selected_baseline']==min(order,key=lambda k:val['selected_validation_brier'][k]),'full baseline validation argmin')
  check(val['selected_matched_baseline']==min([k for k in order if k!='power_full'],key=lambda k:val['selected_validation_brier'][k]),'matched baseline validation argmin')
 for name in names:
  f=out/name
  root_h=root_files[name]['original_sha256' if dirname=='villin' else 'reproduction_sha256']
  check(sha(f)==root_h,'reported comparison file hashes')
for name in names:compare(read(R/'results/villin'/name),read(R/'results/reproduced_villin'/name),name)
check(numeric==root_comparison['numeric_fields_compared']==1105791,'reported numeric field count')
check(other==root_comparison['other_scalar_fields_compared']==15848,'reported other scalar field count')
report['exact_full_comparison']={'files':len(names),'numeric_fields':numeric,'other_scalar_fields':other,'max_absolute_difference':max_abs,'excluded_root_key_occurrences':dict(excluded_occurrences),'excluded_paths':excluded_paths,'justification':{'completed_utc':'Execution timestamp','seconds':'Elapsed runtime','method_prediction_and_metrics_seconds':'Per-method elapsed runtime','validation_file_sha256':'Differs because validation timestamps/runtime differ; independently checked as an exact within-run file hash, and validation scientific content independently compared.'}}
report['seeds_and_selection']={'unique_seed_keys_per_execution':{k:len(v) for k,v in seed_sets.items()},'validation_candidate_families_recomputed':selection_checks,'candidate_mean_tolerance':2e-15,'reviewer_check_correction':'Initial reviewer used exact equality for independent numpy mean versus sequential production mean. Raw one-column mean differs by <=8.33e-17 due to reduction order; corrected reviewer tolerance to 2e-15, without changing any benchmark output.','chronology':chronology,'scope':'No per-phase explicit test-start timestamp is recorded. Test start inferred from its completion time minus elapsed monotonic seconds; file mtimes also order correctly. Production control flow writes validation before entering test phase and before loading test contacts for evaluation; all input files are hashed at invocation start.'}
progress(f'Full reproduction exact: {numeric} numeric + {other} other scalar fields; seed coverage and validation selection confirmed.')
# Independent metrics: positive runs found with a while scanner, not production edges.
def own_metrics(pred,true,freq):
 true=true.astype(bool);calls=pred>=.5;error=pred-true;sq=error**2
 tp=int((calls&true).sum());fp=int((calls&~true).sum());fn=int((~calls&true).sum())
 rare=freq<=.1;rare_pos=true&rare
 short=np.zeros(true.shape,bool);n=len(true);events=0
 for c in range(true.shape[1]):
  i=0
  while i<n:
   if not true[i,c]:i+=1;continue
   start=i
   while i<n and true[i,c]:i+=1
   if start>0 and i<n and i-start<=3:short[start:i,c]=True;events+=1
 transition=np.zeros(true.shape,bool)
 for i in range(1,n):
  changed=true[i]!=true[i-1];transition[i,changed]=True;transition[i-1,changed]=True
 def mean(a,mask):return float(a[mask].sum()/mask.sum()) if mask.any() else None
 metrics={'brier':float(sq.sum()/sq.size),'mae':float(np.abs(error).sum()/error.size),'occupancy_mae':float(np.abs((pred.sum(axis=0)-true.sum(axis=0))/n).mean()),'occupancy_signed_bias':float(error.sum()/error.size),'precision':tp/(tp+fp) if tp+fp else None,'recall':tp/(tp+fn) if tp+fn else None,'f1':2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None,'rare_brier':float(sq[:,rare].sum()/sq[:,rare].size) if rare.any() else None,'rare_positive_recall':mean(calls,rare_pos),'short_positive_recall':mean(calls,short),'transition_brier':mean(sq,transition)}
 counts={'rare_features':int(rare.sum()),'rare_positive':int(rare_pos.sum()),'transition':int(transition.sum()),'short_positive':int(short.sum()),'brief_event_count':events}
 return metrics,counts
# Dense linear algebra from the explicit stochastic matrix, independent of banded production solver.
def dense_predictions(z,settings):
 z=z.astype(float);n=len(z);a=np.eye(n)
 if n>1:
  a=np.diag(np.full(n,.5))+np.diag(np.full(n-1,.25),1)+np.diag(np.full(n-1,.25),-1);a[0,0]=a[-1,-1]=.75
 hcache={};out={};residual=0.
 for method in ['raw','linear','paper332','power_full','power_matched']:
  if method=='raw':out[method]=z;continue
  t=settings[method]['t']
  if t not in hcache:
   matrix=np.eye(n)+t*(np.eye(n)-a);h=np.linalg.solve(matrix,z)
   residual=max(residual,float(np.max(np.abs(matrix@h-z))))
   h[:,z.sum(axis=0)==0]=0;h[:,z.sum(axis=0)==n]=1;hcache[t]=np.clip(h,0,1)
  h=hcache[t]
  if method=='linear':out[method]=h
  elif method=='paper332':out[method]=3*h*h-2*h*h*h
  else:
   alpha=settings[method]['alpha'];numer=h**alpha;out[method]=numer/(numer+(1-h)**alpha)
 return out,residual

def own_corruption(x,xyz,pairs,scenario,words):
 rng=np.random.default_rng(np.random.SeedSequence(words));eps=scenario['level']
 if scenario['kind']=='bitflip':return x!=(rng.random(x.shape)<eps)
 if scenario['kind']=='jitter':return contacts(xyz.astype(float)+rng.normal(scale=eps,size=xyz.shape),pairs,scenario['cutoff_A'])
 e=np.empty(x.shape,bool);e[0]=rng.random(x.shape[1])<eps
 for i in range(1,len(x)):
  p=np.full(x.shape[1],eps*(1-scenario['rho']));p[e[i-1]]=eps+(1-eps)*scenario['rho'];e[i]=rng.random(x.shape[1])<p
 return x!=e
spot=[];metric_n=0;max_metric_delta=0.;max_residual=0.
runmap={r['id']:r for r in manifest['runs']}
for scenario_name in ['flip10_cut8','jitter050_cut8','correlated10_cut8']:
 test=outputs[('villin',scenario_name,'test')];val=outputs[('villin',scenario_name,'validation')]
 for idx in [0,1,47]:
  rr=test['per_run'][idx];run=runmap[rr['run_id']]
  with np.load(P/run['contacts']['8']) as d:x=d['X'];pairs=d['pairs'];freq=d['frequency']
  xyz,_=source_arrays[run['id']]
  # First and last draws give a deterministic spread without looking at performance.
  for ni in [0,15]:
   row=rr['noise_rows'][ni];z=own_corruption(x,xyz,pairs,test['scenario'],row['seed_words']);preds,res=dense_predictions(z,test['selected_parameters']);max_residual=max(max_residual,res);delta=0.
   for method,pred in preds.items():
    actual,counts=own_metrics(pred,x,freq);check(counts==rr['stratum_counts'],'independent stratum counts')
    for metric,value in actual.items():
     target=row['methods'][method][metric];metric_n+=1
     check((value is None)==(target is None),'spot metric null identity')
     if value is not None:
      diff=abs(value-target);delta=max(delta,diff);max_metric_delta=max(max_metric_delta,diff)
      check(diff<2e-12,'spot independent metric reproduction tolerance')
   spot.append({'scenario':scenario_name,'run_id':run['id'],'noise_index':ni,'shape':list(x.shape),'methods':list(preds),'metrics_per_method':11,'max_absolute_metric_difference':delta,'max_dense_solve_residual':res})
 # Reconstruct first validation row for all seven cubic/linear candidates, plus raw.
 vr=val['records'][0];run=runmap[vr['run_id']]
 with np.load(P/run['contacts']['8']) as d:x=d['X'];pairs=d['pairs']
 z=own_corruption(x,source_arrays[run['id']][0],pairs,val['scenario'],vr['seed_words'])
 for t in [1,2,4,8,16,32,64]:
  settings={'raw':{},'linear':{'t':t},'paper332':{'t':t},'power_full':{'t':t,'alpha':3},'power_matched':{'t':t,'alpha':3}}
  preds,res=dense_predictions(z,settings);max_residual=max(max_residual,res)
  for method in ['linear','paper332']:
   j=val['candidate_grids'][method].index({'t':t});diff=abs(float(np.mean((preds[method]-x)**2))-vr['candidate_brier'][method][j]);check(diff<2e-12,'spot validation candidate score');metric_n+=1;max_metric_delta=max(max_metric_delta,diff)
# Clean controls on same first, short second, last selected test trajectories.
clean=read(R/'results/villin/primary_clean.json')
for idx in [0,1,47]:
 rr=clean['per_run'][idx];run=runmap[rr['run_id']]
 with np.load(P/run['contacts']['8']) as d:x=d['X'];freq=d['frequency']
 preds,res=dense_predictions(x,clean['selected_parameters']);max_residual=max(max_residual,res)
 for method,pred in preds.items():
  actual,counts=own_metrics(pred,x,freq);check(counts==rr['stratum_counts'],'clean stratum reconstruction')
  for metric,value in actual.items():
   target=rr['methods'][method][metric];metric_n+=1;check((value is None)==(target is None),'clean metric null identity')
   if value is not None:
    diff=abs(value-target);max_metric_delta=max(max_metric_delta,diff);check(diff<2e-12,'clean metric reconstruction')
report['independent_numerical_reconstruction']={'implementation':'This script imports no project module. Uses explicit dense path transition matrix and numpy.linalg.solve (production uses SciPy banded solve), polynomial cubic expression, direct positive-power ratio, separate positive-run scanner and metrics, separate noise/contacts routines.','methods':['raw','linear','paper332','power_full','power_matched'],'test_cases':spot,'validation_candidate_scores':42,'clean_control_runs':3,'metric_and_candidate_values_compared':metric_n,'absolute_tolerance':2e-12,'max_absolute_metric_difference':max_metric_delta,'max_dense_solve_residual':max_residual,'selection_rule':'First, second, last test trajectory by output order; first/last noise draws; first validation record. Chosen without reading case performance.'}
progress(f'Independent reconstruction: {metric_n} metric/candidate values; maximum difference {max_metric_delta:g}; dense residual {max_residual:g}.')
# Confirm actual interpreter package isolation, not just installation logs.
report['fresh_environment']={'python_executable':sys.executable,'python_version':platform.python_version(),'numpy_version':np.__version__,'venv_config':Path(sys.prefix,'pyvenv.cfg').read_text(),'packages':json.loads(subprocess.check_output([sys.executable,'-m','pip','list','--format=json'],text=True)),'install_log_sha256':sha(R/'reviews/fresh_environment_install.log')}
check('include-system-site-packages = false' in report['fresh_environment']['venv_config'],'fresh venv excludes system packages')
check({x['name']:x['version'] for x in report['fresh_environment']['packages']}=={'numpy':'2.3.5','pip':'25.0.1','scipy':'1.17.0'},'isolated benchmark dependency package set')
report['check_counts']=dict(checks);report['checks_total']=sum(checks.values());report['script_sha256']=sha(__file__);report['status']='pass_with_scope_limits';report['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
(R/'reviews/05_reproduction_initial.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
progress('PASS: reviews/05_reproduction_initial.json written; no production/protocol/results files modified.')
