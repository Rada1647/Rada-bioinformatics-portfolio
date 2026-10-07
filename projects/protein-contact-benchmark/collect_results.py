"""Summarize immutable benchmark JSON; never tune or alter benchmark results."""
from pathlib import Path
import argparse,datetime,json,hashlib,csv
R=Path(__file__).resolve().parent
LABELS={'raw':'Raw noisy contacts','linear':'Linear resolvent','paper332':'#332 cubic','oracle_debiased_linear':'Known-noise oracle','hard_linear':'Hard threshold','moving_average':'Moving average','gaussian':'Gaussian filter','median':'Temporal median','power_full':'Power sharpening (63)','power_matched':'Power alpha=3 (7)'}
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--allow-partial',action='store_true');a=ap.parse_args()
 p=load(R/'protocol/protocol_v1.1.json');summary={'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'protocol_sha256':digest(R/'protocol/protocol_v1.1.json'),'datasets':{},'ubiquitin_qualification':load(R/'data/ubiquitin/qualification.json')}
 long=[]
 for dataset in ['t4','villin']:
  d=R/'results'/dataset
  expected=[d/(s['name']+'_test.json') for s in p['scenarios']]
  if not all(x.exists() for x in expected):
   if a.allow_partial:continue
   raise SystemExit('Incomplete '+dataset+' results')
  scenarios=[]
  for f in expected:
   x=load(f);b=x['selected_baseline'];matched=x['selected_matched_baseline'];agg=x['aggregate'];c=agg['paper332']['brier']
   s={'name':x['scenario']['name'],'baseline':b,'matched_baseline':matched,'cubic_brier':c,'baseline_brier':agg[b]['brier'],'linear_brier':agg['linear']['brier'],'gain_vs_baseline_percent':100*(agg[b]['brier']-c)/agg[b]['brier'],'gain_vs_linear_percent':100*(agg['linear']['brier']-c)/agg['linear']['brier'],'selected_parameters':x['selected_parameters'],'aggregate':agg,'source_file_sha256':digest(f),'same_t_threshold_disagreements':x['same_t_threshold_disagreements']}
   scenarios.append(s)
   for method,m in agg.items():
    long.append({'dataset':dataset,'scenario':s['name'],'method':method,**m})
  x=load(expected[0]);manifest=load(R/'data/prepared'/dataset/'manifest.json');counts={split:sum(r['split']==split for r in manifest['runs']) for split in ['discovery','validation','test']}
  primary={'aggregate':x['aggregate'],'selected_parameters':x['selected_parameters'],'baseline':x['selected_baseline'],'matched_baseline':x['selected_matched_baseline'],'criterion':x.get('primary_success',x.get('stress_success_descriptive')),'matched_criterion':x['matched_budget_success_secondary'],'per_run':[{'run_id':r['run_id'],'shape':r['shape'],'aggregate':r['aggregate'],'stratum_counts':r['stratum_counts'],'paired_comparisons':r['paired_comparisons']} for r in x['per_run']],'same_t_threshold_disagreements':x['same_t_threshold_disagreements'],'same_t_linear_aggregate':x['same_t_linear_aggregate'],'contributing_runs':x['aggregate_contributing_run_count']}
  clean=load(d/'primary_clean.json')
  coverage=load(R/'data/prepared'/dataset/'coverage.json')
  selected_counts={cut:coverage[cut]['selected_pair_count'] for cut in ['7','8','9']}
  cov=[v for v in coverage['8']['runs'].values() if v['split']=='test']
  totalpos=sum(v['all_positive_cells'] for v in cov);abspos=sum(v['discovery_absent_positive_cells'] for v in cov)
  summary['datasets'][dataset]={'role':p['datasets'][dataset]['role'],'split_counts':counts,'n_ca':manifest['runs'][0]['n_ca'],'selected_features':selected_counts,'primary':primary,'clean':clean['aggregate'],'scenarios':scenarios,'primary_omitted_discovery_absent_positive_fraction_pooled':abspos/totalpos if totalpos else None,'primary_selected_positive_fraction_pooled':sum(v['selected_positive_cells'] for v in cov)/totalpos if totalpos else None,'primary_omitted_total_positive_fraction_pooled':sum(v['omitted_positive_cells'] for v in cov)/totalpos if totalpos else None,'primary_coverage_per_run':{k:v for k,v in coverage['8']['runs'].items() if v['split']=='test'},'better_than_linear_scenarios':sum(s['gain_vs_linear_percent']>0 for s in scenarios),'better_than_selected_baseline_scenarios':sum(s['gain_vs_baseline_percent']>0 for s in scenarios),'result_files':[str(f.relative_to(R)) for f in expected]}
 (R/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
 if long:
  with (R/'metrics.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(long[0]));w.writeheader();w.writerows(long)
 print(json.dumps({d:{'criterion':v['primary']['criterion']['status'],'primary':v['scenarios'][0]} for d,v in summary['datasets'].items()},indent=2))
if __name__=='__main__':main()
