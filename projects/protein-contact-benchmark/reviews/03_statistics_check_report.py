"""Check report statistical rows independently of build_report.py."""
from pathlib import Path
import collections, datetime, hashlib, json, math, re
from scipy.stats import t
R=Path(__file__).resolve().parents[1]
L={'raw':'Raw contacts','linear':'Linear resolvent','paper332':'#332 cubic','oracle_debiased_linear':'Known-noise oracle','hard_linear':'Hard threshold','moving_average':'Moving average','gaussian':'Gaussian filter','median':'Temporal median','power_full':'Power (63 settings)','power_matched':'Power alpha=3 (7)'}
def load(p):return json.loads((R/p).read_text())
def mean(v):
    x=[z for z in v if z is not None]
    return math.fsum(x)/len(x) if x else None
def f(v,n=6):return 'not evaluable' if v is None else f'{v:.{n}f}'
def pct(v):return 'not evaluable' if v is None else f'{100*v:.2f}%'
text=(R/'RESULTS.txt').read_text();fails=[];counts=collections.Counter();datasets={};coverage={}
def check(expected,nums,tag):
    counts['report_rows_checked']+=1;counts['report_numeric_fields_checked']+=nums
    if expected not in text:fails.append({'check':tag,'expected_text':expected})
for dataset,title in [('t4','T4'),('villin','Villin')]:
    primary=load(f'results/{dataset}/flip10_cut8_test.json');baseline=primary['selected_baseline'];matched=primary['selected_matched_baseline'];ag={}
    for m in primary['selected_parameters']:
        ag[m]={k:mean([mean([row['methods'][m][k] for row in run['noise_rows']]) for run in primary['per_run']]) for k in primary['aggregate'][m]}
        x=ag[m];check(' | '.join([L[m],f(x['brier']),pct(x['rare_positive_recall']),pct(x['short_positive_recall']),f(x['occupancy_mae'])]),4,dataset+'.primary.'+m)
    gain=100*(ag[baseline]['brier']-ag['paper332']['brier'])/ag[baseline]['brier']
    check(' | '.join([title+' primary' if dataset=='t4' else title+' stress',f(ag['paper332']['brier']),L[baseline],f(ag[baseline]['brier']),f'{gain:+.2f}%']),3,dataset+'.headline')
    for scenario in load('protocol/protocol_v1.1.json')['scenarios']:
        z=load(f"results/{dataset}/{scenario['name']}_test.json");base=z['selected_baseline'];a={m:mean([mean([row['methods'][m]['brier'] for row in run['noise_rows']]) for run in z['per_run']]) for m in ['paper332','linear',base]}
        g=100*(a[base]-a['paper332'])/a[base];gl=100*(a['linear']-a['paper332'])/a['linear']
        check(' | '.join([title,scenario['name'],f(a['paper332'],5),L[base],f'{g:+.2f}%',f'{gl:+.2f}%']),3,dataset+'.scenario.'+scenario['name'])
    clean=load(f'results/{dataset}/primary_clean.json')
    for m in ['raw','linear','paper332',baseline]:
        x={k:mean([r['methods'][m][k] for r in clean['per_run']]) for k in ['brier','short_positive_recall','occupancy_mae']}
        check(' | '.join([title,L[m],f(x['brier']),pct(x['short_positive_recall']),f(x['occupancy_mae'])]),3,dataset+'.clean.'+m)
    cv=load(f'data/prepared/{dataset}/coverage.json')['8'];held=[r for r in cv['runs'].values() if r['split']=='test'];total=sum(r['all_positive_cells'] for r in held);sel=sum(r['selected_positive_cells'] for r in held);absent=sum(r['discovery_absent_positive_cells'] for r in held)
    coverage[dataset]={'all_eligible_pairs':cv['all_eligible_pair_count'],'selected_pairs':cv['selected_pair_count'],'selected_positive_cells':sel,'all_positive_cells':total,'selected_positive_fraction':sel/total,'omitted_positive_fraction':1-sel/total,'discovery_absent_positive_fraction':absent/total}
    check(f'Held-out positive-cell fraction in discovery-constant-absent pairs: {100*absent/total:.4f}%',1,dataset+'.absent_coverage')
    check(f'Selected features contain {100*sel/total:.4f}% of all held-out positive cells; total omitted positive-cell fraction is {100*(1-sel/total):.4f}%',2,dataset+'.selected_and_total_omitted_coverage')
    if dataset=='t4':
        cis=[]
        for run in primary['per_run']:
            b={k:mean([r['methods'][baseline][k] for r in run['noise_rows']]) for k in ['brier','rare_positive_recall','short_positive_recall']};p={k:mean([r['methods']['paper332'][k] for r in run['noise_rows']]) for k in b}
            gain=b['brier']-p['brier'];rare=b['rare_positive_recall']-p['rare_positive_recall'];brief=b['short_positive_recall']-p['short_positive_recall']
            check(' | '.join([run['run_id'],f(gain),pct(rare),pct(brief),str(gain>0),str(rare<=.02),str(brief<=.02)]),3,'t4.gates.'+run['run_id'])
            diff=[r['methods'][baseline]['brier']-r['methods']['paper332']['brier'] for r in run['noise_rows']];dm=mean(diff);half=float(t.ppf(.975,len(diff)-1))*math.sqrt(math.fsum((x-dm)**2 for x in diff)/(len(diff)-1)/len(diff));lo,hi=dm-half,dm+half
            check(' | '.join([run['run_id'],f(dm),f(lo),f(hi)]),3,'t4.CI.'+run['run_id']);cis.append({'run_id':run['run_id'],'mean_brier_gain':dm,'conditional_noise_95pct_ci':[lo,hi]})
        gm=100*(ag[matched]['brier']-ag['paper332']['brier'])/ag[matched]['brier']
        check(f'Its validation-selected comparator is {L[matched]}, Brier {f(ag[matched]["brier"])}; cubic relative gain {gm:+.2f}%.',2,'t4.matched_budget')
        datasets[dataset]={'selected_baseline':baseline,'matched_baseline':matched,'conditional_intervals':cis}
for phrase in ['not human peer review','at least 2% lower mean Brier error','no more than two percentage points','Villin is reported separately and cannot determine the T4 primary decision','pointwise and unadjusted for multiple comparisons','not independent biological replicates','positive gain means lower error']:
    counts['scope_phrases_checked']+=1
    if phrase.lower() not in text.lower():fails.append({'missing_scope_phrase':phrase})
scenario_rows=re.findall(r'^(?:T4|Villin) \| (?:flip|jitter|correlated)\S+ \|',text,re.M)
if len(scenario_rows)!=16:fails.append({'expected_scenario_rows':16,'observed':len(scenario_rows)})
out={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'pass' if not fails else 'issues_found','results_txt_sha256':hashlib.sha256((R/'RESULTS.txt').read_bytes()).hexdigest(),'counts':dict(counts),'scenario_rows_count':len(scenario_rows),'failures':fails,'evidence':datasets,'positive_coverage':coverage,'interpretation':'All statistical report rows checked against independently recomputed means/contrasts from archived noise rows; 95% Student-t intervals are approximate pointwise conditional-noise intervals, not simultaneous or biological confidence statements.'}
(R/'reviews/03_statistics_report_checks.json').write_text(json.dumps(out,indent=2)+'\n')
report=load('reviews/03_statistics_final.json');report['scope']='Final completed T4 confirmation plus Villin adaptive stress set: all 16 dataset-scenarios, validation selections, source-derived coverage and report statistics';report['report_checks']=out
report['findings']=[x for x in report['findings'] if x.get('id') not in ['t4_primary','coverage_scope']]
report['findings'].append({'id':'t4_primary','severity':'interpretation','finding':'T4 primary cubic Brier 0.04036348805980258 versus selected full/matched power 0.03719128753665906: 8.5294184% worse. Every composite gate fails, with 0/3 positive Brier gains, 0/3 rare guardrail passes and 0/3 brief guardrail passes. Rare losses are 8.7996–10.8051 percentage points; brief losses are 11.3532–11.6204 points. All strata are evaluable.'})
report['findings'].append({'id':'coverage_scope','severity':'interpretation','resolution':'Report now explicitly states both selected-positive coverage and total omitted-positive fraction; numeric statements independently checked.','finding':'T4 primary scores cover 523/12,561 eligible pairs and 63.3229% of held-out positive cells. The 0.1983% positives in discovery-absent pairs are only part of omission: discovery-constant-present pairs are also excluded. Total omitted positive-cell fraction is 36.6771%. This is protocol-consistent feature-restricted performance, not all-contact accuracy.'})
if fails:report['status']='issues_found';report['failures'].extend(fails)
report['report_numeric_fields_checked']=counts['report_numeric_fields_checked']
(R/'reviews/03_statistics_final.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps({'status':out['status'],'failures':fails,'counts':dict(counts),'coverage':coverage},indent=2))
