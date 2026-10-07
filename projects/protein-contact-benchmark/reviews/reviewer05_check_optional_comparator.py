"""Check optional manifest-hash exclusion on disposable fixtures, never results."""
from pathlib import Path
import datetime, hashlib, json, subprocess, sys, tempfile
R=Path(__file__).resolve().parents[1]
results=[]
with tempfile.TemporaryDirectory(prefix='protein332_comparator_review_') as tmp:
 base=Path(tmp);a=base/'original';b=base/'regenerated';a.mkdir();b.mkdir()
 names=[f'scenario{i}_{phase}.json' for i in range(8) for phase in ['validation','test']]+['primary_clean.json']
 fixture={'prepared_manifest_sha256':'first-manifest','prepared_inputs_sha256':{'run1':{'coordinates':'same-coordinate-hash','contacts':{'8':'same-contact-hash'}}},'protocol_sha256':'same-protocol','code_sha256':{'methods.py':'same-code'},'selected_parameters':{'paper332':{'t':8}},'brier':0.123,'completed_utc':'old-time','seconds':1.5,'validation_file_sha256':'old-validation-file'}
 def run(label,modifications,optional,expect_pass):
  changed=json.loads(json.dumps(fixture));changed.update(modifications)
  for name in names:
   (a/name).write_text(json.dumps(fixture));(b/name).write_text(json.dumps(changed))
  out=base/'comparison.json';args=[sys.executable,str(R/'compare_reproduction.py'),str(a),str(b),'--out',str(out)]
  if optional:args.append('--allow-regenerated-manifest')
  completed=subprocess.run(args,text=True,capture_output=True);report=json.loads(out.read_text())
  assert (completed.returncode==0)==expect_pass,(label,completed.stdout,completed.stderr)
  assert ('prepared_manifest_sha256' in report['excluded_metadata_keys'])==optional
  results.append({'case':label,'optional':optional,'expected_pass':expect_pass,'exit_code':completed.returncode,'difference_count':report['difference_count'],'excluded_metadata_keys':report['excluded_metadata_keys']})
 run('Default rejects changed preparation hash',{'prepared_manifest_sha256':'new-manifest'},False,False)
 run('Explicit optional mode allows changed preparation hash',{'prepared_manifest_sha256':'new-manifest'},True,True)
 run('Optional mode rejects changed prepared contact hash',{'prepared_manifest_sha256':'new-manifest','prepared_inputs_sha256':{'run1':{'coordinates':'same-coordinate-hash','contacts':{'8':'different-contact-hash'}}}},True,False)
 run('Optional mode rejects changed benchmark code hash',{'prepared_manifest_sha256':'new-manifest','code_sha256':{'methods.py':'different-code'}},True,False)
 run('Optional mode rejects changed protocol hash',{'prepared_manifest_sha256':'new-manifest','protocol_sha256':'different-protocol'},True,False)
 run('Optional mode rejects changed scientific metric',{'prepared_manifest_sha256':'new-manifest','brier':.124},True,False)
 run('Optional mode rejects changed selected smoothing parameter',{'prepared_manifest_sha256':'new-manifest','selected_parameters':{'paper332':{'t':16}}},True,False)
report={'reviewer':'05','checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'pass','comparator_sha256':hashlib.sha256((R/'compare_reproduction.py').read_bytes()).hexdigest(),'fixture_scope':'Disposable synthetic 17-output fixtures to test CLI exclusion behavior; no benchmark rerun or production output mutation.','tests':results,'readme_review':{'raw_reconstruction_paths':'Source directories renamed to reference; exact extractor+metadata copies recreate required relative source paths and preparation uses them correctly.','manifest_exception_reason':'Preparation timestamps and source-acquisition provenance change the manifest hash on independent extraction. Exact prepared-array file hashes and all scientific result fields remain required.','requested_documentation_corrections':['Install requirements-analysis.txt explicitly before extraction commands and use the same .venv interpreter throughout.','Clarify T4 script retains final md1us8 raw XTC for source-level review although md1us4–7 newly downloaded XTCs are removed.']},'limits':['Comparator itself is not a provenance attestation. Optional mode permits any preparation-manifest bytes to differ; users must inspect regenerated provenance.','Default production result comparison and independent reviewer hash checks remain stronger for the existing prepared inputs.','No claim of independently rerunning source extraction follows from these tests.']}
(R/'reviews/05_optional_comparator_review.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'status':'pass','cases':len(results),'comparator_sha256':report['comparator_sha256']}))
