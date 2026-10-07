"""Create a readable report and publication-style figure from fixed outputs."""
from pathlib import Path
import argparse,datetime,html,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R=Path(__file__).resolve().parent
LABEL={'raw':'Raw contacts','linear':'Linear resolvent','paper332':'#332 cubic','oracle_debiased_linear':'Known-noise oracle','hard_linear':'Hard threshold','moving_average':'Moving average','gaussian':'Gaussian filter','median':'Temporal median','power_full':'Power (63 settings)','power_matched':'Power alpha=3 (7)'}
def load(p):return json.loads(p.read_text())
def num(v,d=6):return 'not evaluable' if v is None else f'{v:.{d}f}'
def pct(v,d=2):return 'not evaluable' if v is None else f'{100*v:.{d}f}%'

def figure(summary):
 fig,axs=plt.subplots(1,2,figsize=(11,5.4),sharey=True)
 methods=list(LABEL)
 for ax,dataset in zip(axs,['t4','villin']):
  data=summary['datasets'][dataset];a=data['primary']['aggregate'];base=data['primary']['baseline']
  values=[a[m]['brier'] for m in methods]
  colors=['#24608c' if m=='paper332' else '#cb7628' if m==base else '#9aa8b3' for m in methods]
  bars=ax.barh(np.arange(len(methods)),values,color=colors,height=.70)
  for bar,m,v in zip(bars,methods,values):
   if m=='oracle_debiased_linear':bar.set_hatch('///');bar.set_edgecolor('#555555')
   ax.text(v+.001,np.where(np.array(methods)==m)[0][0],f'{v:.5f}',va='center',fontsize=8)
  ax.set_yticks(np.arange(len(methods)),[LABEL[m] for m in methods]);ax.set_ylim(len(methods)-.5,-.5);ax.set_xlim(0,.125)
  ax.set_xlabel('Mean Brier error (lower is better)',fontsize=10)
  ax.set_title('T4 lysozyme: primary' if dataset=='t4' else 'Villin: adaptive stress test',fontsize=11)
  ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',alpha=.15);ax.set_axisbelow(True)
 fig.suptitle('Prespecified 10% independent-flip scenario',fontsize=13)
 fig.text(.5,.015,'Equal weight per held-out trajectory. Orange: control selected on validation. Hatched: oracle supplied the true error rate.',ha='center',fontsize=8)
 fig.tight_layout(rect=[0,.055,1,.96]);fig.savefig(R/'primary_results.png',dpi=220);plt.close(fig)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--pdf',action='store_true');args=ap.parse_args()
 s=load(R/'summary.json')
 if set(s['datasets'])!={'t4','villin'}:raise SystemExit('Both complete benchmark datasets are required for the report')
 refs=load(R/'references.json');refnum={x['id']:i+1 for i,x in enumerate(refs)}
 def cite(*ids):return '['+', '.join(str(refnum[i]) for i in ids)+']'
 sections=[]
 def heading(x):sections.append(('heading',x))
 def para(x):sections.append(('paragraph',x))
 def table(headers,rows):sections.append(('table',(headers,rows)))
 def page():sections.append(('page',None))
 primary=s['datasets']['t4']['primary'];decision=primary['criterion']['status'];v=s['datasets']['villin']['primary']
 heading('Prospective benchmark of a theorem-linked protein-contact smoother')
 para('Computational research report | Frozen protocol 1.1 | AI-assisted analysis and five independent AI technical review roles. This is not human peer review.')
 para('The T4 primary practical-success decision is '+decision.upper()+'. The decision requires at least 2% lower mean Brier error than the validation-selected strongest prespecified nonoracle control, improvement on every test run, and no more than two percentage points of rare- or brief-contact recall lost on any run. These are engineering thresholds, not a theorem-derived guarantee.')
 rows=[]
 for name,title in [('t4','T4 primary'),('villin','Villin stress')]:
  d=s['datasets'][name];p=d['primary'];a=p['aggregate'];b=p['baseline'];gain=d['scenarios'][0]['gain_vs_baseline_percent']
  rows.append([title,num(a['paper332']['brier']),LABEL[b],num(a[b]['brier']),f'{gain:+.2f}%'])
 table(['Experiment','Cubic error','Selected control','Control error','Cubic gain'],rows)
 para('Positive gain means lower error from the cubic; negative gain means the selected control has lower error. Villin is reported separately and cannot determine the T4 primary decision. Ubiquitin was excluded by the predeclared provenance gate, so there is no ubiquitin denoising result.')
 sections.append(('figure','primary_results.png'))
 para('The method computes H = [I + t(I-A)]^-1 Z and Y = 3H^2 - 2H^3 for binary contact input Z on a reflecting temporal path. This is a faithful binary specialization of the expected median of three independent geometric-walk endpoints. The released geometric inequality does not guarantee lower Brier reconstruction error or improved protein physics. The cubic itself predates the released proof. '+cite('OpenAI2026','ChengWangXiang2026'))
 page();heading('T4 lysozyme: held-out simulation runs')
 d=s['datasets']['t4'];p=d['primary'];a=p['aggregate']
 para('The five deposited one-microsecond runs were downloaded completely and checked against their published MD5s. Run md1us4 supplied discovery features; md1us5 supplied validation settings; md1us6, md1us7 and md1us8 were held out for testing. We used 162 CA atoms at 200 ps intervals. The deposited sequence is reported exactly in the provenance audit and is not labelled wild type. The original study concerns side-chain NMR relaxation; its ABSURDer reweighting pipeline is not implemented here. '+cite('Kummerer2021','KummererData2020'))
 para('Relative to the current UniProt P00720 sequence, the deposited 162-residue construct has differences at positions 12, 54, 97 and 137 and lacks the terminal NL. This is a sequence comparison, not a claim that every difference was an engineered mutation. Periodic-boundary reconstruction uses the bonded backbone; all 2,430,486 CA coordinate values in the first run matched full-protein reconstruction exactly. UniProt P00720 was accessed on 7 October 2026; the reference bytes and their hash are retained. '+cite('UniProtP00720'))
 table(['Method','Brier','Rare recall','Brief recall','Occupancy error'],[[LABEL[m],num(a[m]['brier']),pct(a[m]['rare_positive_recall']),pct(a[m]['short_positive_recall']),num(a[m]['occupancy_mae'])] for m in LABEL])
 para('Recall is the fraction of reference-positive cells recovered at a 0.5 decision threshold. Brief recall concerns cells within complete positive runs of at most 600 ps; it is not an event-detection rate. Occupancy error compares each contact\'s average predicted and reference occurrence.')
 st=p['same_t_linear_aggregate']
 para('At the cubic-selected t='+str(p['selected_parameters']['paper332']['t'])+', the same-t linear comparator has Brier '+num(st['brier'])+' and the cubic has '+num(a['paper332']['brier'])+'. Their threshold disagreements total '+str(p['same_t_threshold_disagreements'])+'; rare and brief recall are identical at this same smoothing scale. The independently tuned linear method in the preceding table uses t='+str(p['selected_parameters']['linear']['t'])+', which explains its different recall.')
 rows=[]
 for r in p['criterion']['per_run']:
  rows.append([r['run_id'],num(r['brier_gain']),pct(r['rare_recall_loss']),pct(r['brief_recall_loss']),str(r['positive_brier_gain']),str(r['rare_guardrail_pass']),str(r['brief_guardrail_pass'])])
 table(['Test run','Brier gain','Rare loss','Brief loss','Gain >0','Rare gate','Brief gate'],rows)
 para('Losses in the preceding table are percentage-point differences, displayed as percentage values: e.g. 2.00% means two percentage points. Required gates use unrounded numbers. Primary settings: '+json.dumps(p['selected_parameters'],sort_keys=True)+'.')
 matched=p['matched_baseline'];mg=p['matched_criterion']['relative_mean_brier_gain']
 para('The matched-budget contrast fixes sharpening alpha=3 using prior adenylate-kinase development information and tunes seven t values, matching the cubic\'s seven candidates. Its validation-selected comparator is '+LABEL[matched]+', Brier '+num(a[matched]['brier'])+'; cubic relative gain '+pct(mg)+'. Its separately recorded practical criterion is '+p['matched_criterion']['status'].upper()+'.')
 table(['Test run','Brier gain vs selected control','Noise-only 95% lower','Noise-only 95% upper'],[[r['run_id'],num(r['paired_comparisons'][p['baseline']]['mean_brier_gain']),num(r['paired_comparisons'][p['baseline']]['noise_only_paired_95pct_ci'][0]),num(r['paired_comparisons'][p['baseline']]['noise_only_paired_95pct_ci'][1])] for r in p['per_run']])
 para('These paired intervals describe corruption-draw variability on each fixed reference trajectory. They are pointwise and unadjusted for multiple comparisons; they do not quantify uncertainty across proteins or establish population significance.')
 page();heading('Villin: a separate adaptive-sampling stress test')
 d=s['datasets']['villin'];p=d['primary'];a=p['aggregate']
 para('The complete checksum-verified archive contains 2,137 short trajectories across three archive groups. The frozen metadata-only selection chose 24 discovery trajectories from group 1, 24 validation from group 2 and 48 test trajectories from group 3, using SHA256 path ordering. Every path stays separate. These groups resemble adaptive campaigns, but their statistical independence and generation seeds are not certified by the source. The analyzed subset is 96 trajectories, not the entire approximately 107-microsecond archive. '+cite('Doerr2016','DoerrData2026','HTMDTutorial2026'))
 para('The topology has 35 CA atoms, with norleucine residues at source positions 65 and 70 and HSP at position 68. This is engineered HP35, not wild-type villin. CA selection includes noncanonical residues. The source tutorial reports 360 K; the original coordinates are already whole and pass adjacent-CA geometry checks. Crystal pH from the corroborating PDB entry is not used to infer simulation pH. '+cite('Kubelka2006','PDB2F4K'))
 table(['Method','Brier','Rare recall','Brief recall','Occupancy error'],[[LABEL[m],num(a[m]['brier']),pct(a[m]['rare_positive_recall']),pct(a[m]['short_positive_recall']),num(a[m]['occupancy_mae'])] for m in LABEL])
 para('The descriptive stress criterion is '+p['criterion']['status'].upper()+'. Cubic and linear predictions at the same t had '+str(p['same_t_threshold_disagreements'])+' threshold disagreements. Both primary tuned methods selected t='+str(p['selected_parameters']['paper332']['t'])+'; their equal classification recall is consistent with the cubic\'s monotonic threshold-preserving transformation. Lower Brier error alone can reflect sharper probability estimates without recovering additional binary contacts.')
 heading('Ubiquitin qualification: no temporal benchmark')
 para('The verified PROTHON archive contains native Q99 and five biased partially folded ensembles, each with 5,000 frames. The article reports five native runs and 40 ps sampling, but inspected archive/author code does not document original run boundaries or chronological assembly. MDAnalysis interprets the DCD header spacing as 0.04888821 ps. We therefore did not infer five consecutive 1,000-frame runs, repair the timing by assumption, or compute denoising scores. The archive remains useful for structural-ensemble work; our exclusion applies to this temporal experiment. '+cite('Aina2023','AinaData2023'))
 page();heading('All prespecified scenarios')
 rows=[]
 for name,title in [('t4','T4'),('villin','Villin')]:
  for x in s['datasets'][name]['scenarios']:
   rows.append([title,x['name'],num(x['cubic_brier'],5),LABEL[x['baseline']],f"{x['gain_vs_baseline_percent']:+.2f}%",f"{x['gain_vs_linear_percent']:+.2f}%"])
 table(['Dataset','Scenario','Cubic error','Selected control','Gain vs control','Gain vs independently tuned linear'],rows)
 para('The independent 10% flip / 8-Angstrom scenario is primary. Other rows are prespecified sensitivity analyses, not opportunities to redefine success. Coordinate-jitter settings are standard deviations in Angstrom per Cartesian component. Correlated errors have stationary prevalence 10% and lag-one correlation 0.8; the error indicator has 0-to-1 and 1-to-0 transition probabilities 0.02 and 0.18, respectively. All control selection used validation scores only.')
 heading('Damage on clean reference input')
 rows=[]
 for name,title in [('t4','T4'),('villin','Villin')]:
  d=s['datasets'][name]
  for m in ['raw','linear','paper332',d['primary']['baseline']]:
   x=d['clean'][m];rows.append([title,LABEL[m],num(x['brier']),pct(x['short_positive_recall']),num(x['occupancy_mae'])])
 table(['Dataset','Method','Clean Brier','Brief recall','Occupancy error'],rows)
 para('All clean-input controls use settings selected on noisy primary validation, without retuning. Unchanged clean input has zero Brier error. Smoothing can erase real changes in the reference trajectory even while reducing artificial noise. This does not establish whether those short contacts are experimentally functional.')
 page();heading('Methods and limitations')
 para('Nonlocal CA pairs have positional sequence separation at least four. Contacts use strict distance <8 Angstrom, with 7/9-Angstrom sensitivities. Distances use float64 arithmetic from stored float32 coordinates. Discovery-variable features are selected separately at each cutoff. Retrospective smoothers use both past and future observations inside each trajectory. No time interpolation, cross-trajectory filtering or test-set tuning occurs.')
 para('The path operator A stays in place with probability one-half at interior frames and moves to either neighboring frame with probability one-quarter. At each endpoint the missing-neighbor probability is added to staying in place. Every contact column and trajectory is processed separately. The output is a set of contact scores, not a reconstructed three-dimensional protein structure.')
 para('Rare contacts have discovery frequency at most 10%. Brief-contact recall counts positive cells within complete runs of at most three sampled frames, excluding runs cut off by a trajectory boundary. Constant discovery features are excluded from the scored set, including constant-present pairs. Coverage of newly appearing contacts is only one coverage diagnostic; it does not measure every form of omitted dynamics.')
 para('Prespecified controls include linear graph smoothing, hard thresholding, moving average, Gaussian filtering, temporal median filtering, a known-flip-rate oracle and power sharpening. These are explicitly defined generic controls, not full reproductions of cited software pipelines. Power sharpening uses h^alpha/[h^alpha+(1-h)^alpha], equivalent to binary probability sharpening with temperature 1/alpha. Brier error is mean single-binary squared error, without doubling for the complementary class. '+cite('SadhanalaWangTibshirani2016','Berthelot2019','Brier1950'))
 para('Contact-based time-series analysis has prior applications in molecular-dynamics event detection, including TimeScapes. That literature motivates the application context; TimeScapes itself is not a tested baseline here. '+cite('Wriggers2009'))
 for name,title in [('t4','T4'),('villin','Villin')]:
  d=s['datasets'][name]
  para(title+' primary selected features: '+str(d['selected_features']['8'])+'. Held-out positive-cell fraction in discovery-constant-absent pairs: '+pct(d['primary_omitted_discovery_absent_positive_fraction_pooled'],4)+ '. Selected features contain '+pct(d['primary_selected_positive_fraction_pooled'],4)+' of all held-out positive cells; total omitted positive-cell fraction is '+pct(d['primary_omitted_total_positive_fraction_pooled'],4)+', including discovery-constant-present pairs. These coverage diagnostics are pooled over cells and do not change the equally weighted trajectory-level primary scores.')
 para('Eight validation and sixteen test corruption draws are used per trajectory. Scores average draws within each run and then give each run equal weight. Noise-only paired 95% intervals are descriptive, pointwise, unadjusted for multiple comparisons, and conditional on fixed trajectories. Frames, contacts and noise draws are not independent biological replicates. The small protein panel, artificial corruption models, fixed 200 ps resolution, finite tuning ranges and discovery-feature selection limit external inference. Boundary-selected hyperparameters are reported as such; no grid expansion followed test inspection.')
 page();heading('Reproduction and independent checks')
 for name in ['villin','t4']:
  path=R/'reviews'/('fresh_'+name+'_comparison.json')
  if path.exists():
   c=load(path);para(name.capitalize()+': '+str(c['numeric_fields_compared'])+' scientific numeric fields compared across all 17 validation/test/clean JSON files; '+str(c['difference_count'])+' differences. Runtime/timestamps and validation-file hashes that include timestamps are explicitly excluded. Reproduction used newly installed benchmark dependencies in a fresh virtual environment on the same host, not an independent computer.')
  else:para(name.capitalize()+': full fresh-environment reproduction is pending; this draft is not final.')
 reviewfiles=sorted((R/'reviews').glob('*_final.json'))
 para('Final independent AI review records present: '+str(len(reviewfiles))+' of five. The review roles cover mathematics/implementation, data/provenance, statistics/protocol, citations/claims and numerical reproduction. Consult the individual reports for checks performed, limitations and resolved findings. AI technical audits do not constitute human peer review.')
 para('Independent spectral implementations agree on continuous predictions to numerical precision, but some scores extremely close to 0.5 give different binary classifications under floating-point rounding. In sampled independent-solver checks, the largest observed classified-metric shift was 0.0558 percentage points in clean-input rare-contact recall. The review records preserve those discrepancies. Exact same-algorithm reproduction is a separate claim from bitwise agreement between independent mathematical implementations.')
 para('From the package directory: create a Python 3.12 virtual environment; install requirements-benchmark.txt; run check_methods.py; then run benchmark_extension.py --dataset t4 --out results/local_t4 and the corresponding villin command. README.txt documents complete acquisition and reconstruction. Core software: NumPy 2.3.5, SciPy 1.17.0 and MDAnalysis 2.9.0; plotting uses Matplotlib. '+cite('MichaudAgrawal2011','Gowers2016','Harris2020','Virtanen2020','Hunter2007'))
 para('Villin and ubiquitin source records explicitly specify CC BY 4.0. The examined T4 record supplies no explicit dataset license. Public availability is not treated as a blanket redistribution grant. Original atomistic trajectories are not bundled. This task did not publish a GitHub repository.')
 para('Frozen protocol SHA256: '+s['protocol_sha256']+'. The private prospective save predates new-protein outcomes; it is not a public registry preregistration. Protocol v1.0 and its pre-outcome amendment to v1.1 are preserved.')
 page();heading('References and source attribution')
 reflines=[]
 for i,x in enumerate(refs,1):
  authors='; '.join(x['authors'])+('; et al.' if x.get('additional_authors') else '')
  line=f"[{i}] {authors} ({x['year']}). {x['title']}. "
  for k in ['journal','booktitle','series','volume','issue','pages']:
   if x.get(k):line+=str(x[k])+'; '
  line+=x.get('url','https://doi.org/'+x.get('doi',''))
  if x.get('accessed'):line+=' (accessed '+x['accessed']+')'
  reflines.append(line);sections.append(('reference',line))
 (R/'REFERENCES.txt').write_text('\n\n'.join(reflines)+'\n')
 text=[]
 for kind,obj in sections:
  if kind in ['heading','paragraph','reference']:text.append(obj)
  elif kind=='table':
   headers,rows=obj;text.append(' | '.join(headers)+'\n'+'\n'.join(' | '.join(str(v) for v in row) for row in rows))
  elif kind=='figure':text.append('Figure: '+obj)
 (R/'RESULTS.txt').write_text('\n\n'.join(text)+'\n')
 figure(s)
 if not args.pdf:
  print('Text report and figure generated; PDF not requested');return
 from reportlab.lib import colors
 from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
 from reportlab.lib.enums import TA_LEFT
 from reportlab.lib.pagesizes import A4
 from reportlab.pdfbase import pdfmetrics
 from reportlab.pdfbase.ttfonts import TTFont
 from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Image,PageBreak,KeepTogether
 font='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf';bold='/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'
 pdfmetrics.registerFont(TTFont('DV',font));pdfmetrics.registerFont(TTFont('DVB',bold))
 styles=getSampleStyleSheet();styles.add(ParagraphStyle(name='BodyR',fontName='DV',fontSize=8.9,leading=12.3,spaceAfter=8));styles.add(ParagraphStyle(name='HeadR',fontName='DVB',fontSize=14,leading=18,spaceBefore=8,spaceAfter=11));styles.add(ParagraphStyle(name='CellR',fontName='DV',fontSize=7.0,leading=9));styles.add(ParagraphStyle(name='RefR',fontName='DV',fontSize=7.5,leading=10.2,spaceAfter=8,wordWrap='CJK'))
 def P(x,style='BodyR'):return Paragraph(html.escape(str(x)),styles[style])
 width=A4[0]-88;story=[]
 for kind,obj in sections:
  if kind=='heading':story.append(P(obj,'HeadR'))
  elif kind=='paragraph':story.append(P(obj))
  elif kind=='reference':story.append(P(obj,'RefR'))
  elif kind=='page':story.append(PageBreak())
  elif kind=='figure':story.append(Image(str(R/obj),width=width,height=width*5.4/11));story.append(Spacer(1,8))
  elif kind=='table':
   headers,rows=obj;n=len(headers)
   weights=([1.3,1,1.6,1,1] if n==5 else [1]*n)
   if headers[0]=='Method':weights=[1.7,1,1,1,1.1]
   if n==6:weights=[.7,1.45,1,1.7,1.1,1.1]
   if n==7:weights=[1.05,1,.85,.85,.75,.85,.85]
   widths=[width*x/sum(weights) for x in weights]
   t=Table([[P(x,'CellR') for x in headers]]+[[P(x,'CellR') for x in row] for row in rows],colWidths=widths,repeatRows=1,hAlign='LEFT')
   t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e3edf4')),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#62849b')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f4f6f8')]),('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),3.5),('BOTTOMPADDING',(0,0),(-1,-1),3.5)]));story.extend([t,Spacer(1,10)])
 def footer(canvas,doc):
  canvas.setFont('DV',7);canvas.setFillColor(colors.HexColor('#60717d'));canvas.drawString(44,24,'Protein-contact benchmark | Protocol 1.1 | AI-assisted research artifact');canvas.drawRightString(A4[0]-44,24,str(doc.page))
 doc=SimpleDocTemplate(str(R/'benchmark_report.pdf'),pagesize=A4,rightMargin=44,leftMargin=44,topMargin=35,bottomMargin=38,title='Prospective protein-contact benchmark',author='AI-assisted computational analysis')
 doc.build(story,onFirstPage=footer,onLaterPages=footer);print('Created benchmark_report.pdf')
if __name__=='__main__':main()
