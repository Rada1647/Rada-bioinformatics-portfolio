"""Compare all scientific fields, excluding explicitly listed run metadata."""
from pathlib import Path
import argparse,datetime,hashlib,json,math

EXCLUDED={'completed_utc','seconds','method_prediction_and_metrics_seconds','validation_file_sha256'}
def canonical(obj):
 if isinstance(obj,dict):return {k:canonical(v) for k,v in obj.items() if k not in EXCLUDED}
 if isinstance(obj,list):return [canonical(x) for x in obj]
 return obj
def main():
 p=argparse.ArgumentParser();p.add_argument('original',type=Path);p.add_argument('reproduction',type=Path);p.add_argument('--out',type=Path,required=True);p.add_argument('--allow-regenerated-manifest',action='store_true',help='Permit a different preparation provenance manifest; all prepared input file hashes and scientific fields must still match exactly');a=p.parse_args()
 if a.allow_regenerated_manifest:EXCLUDED.add('prepared_manifest_sha256')
 names=sorted([p.name for p in a.original.glob('*_validation.json')]+[p.name for p in a.original.glob('*_test.json')]+['primary_clean.json'])
 if len(names)!=17:raise RuntimeError(f'Expected8validation+8test+1clean, found{len(names)}')
 numeric=0;other=0;differences=[];maxabs=0.;files=[]
 def compare(x,y,path):
  nonlocal numeric,other,maxabs
  if isinstance(x,bool) or x is None or isinstance(x,str):
   other+=1
   if type(x)!=type(y) or x!=y:differences.append({'path':path,'original':x,'reproduction':y})
  elif isinstance(x,(int,float)) and isinstance(y,(int,float)):
   numeric+=1;delta=abs(x-y);maxabs=max(maxabs,delta)
   if delta!=0:differences.append({'path':path,'original':x,'reproduction':y,'absolute_difference':delta})
  elif isinstance(x,dict) and isinstance(y,dict):
   if x.keys()!=y.keys():differences.append({'path':path,'key_difference':sorted(set(x)^set(y))})
   for k in x.keys()&y.keys():compare(x[k],y[k],path+'/'+k)
  elif isinstance(x,list) and isinstance(y,list):
   if len(x)!=len(y):differences.append({'path':path,'lengths':[len(x),len(y)]})
   for i,(u,v) in enumerate(zip(x,y)):compare(u,v,path+'/'+str(i))
  else:differences.append({'path':path,'type_difference':[type(x).__name__,type(y).__name__]})
 for name in names:
  x,y=a.original/name,a.reproduction/name
  if not y.exists():raise RuntimeError('Missing reproduced file '+str(y))
  compare(canonical(json.loads(x.read_text())),canonical(json.loads(y.read_text())),name)
  files.append({'file':name,'original_sha256':hashlib.sha256(x.read_bytes()).hexdigest(),'reproduction_sha256':hashlib.sha256(y.read_bytes()).hexdigest()})
 result={'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'original':str(a.original),'reproduction':str(a.reproduction),'files':files,'excluded_metadata_keys':sorted(EXCLUDED),'numeric_fields_compared':numeric,'other_scalar_fields_compared':other,'difference_count':len(differences),'max_absolute_difference':maxabs,'differences':differences[:100],'status':'exact_match' if not differences else 'differences','scope':'All17 scientific output JSON files; invocation timestamps/paths and runtime metadata are not reproducibility targets. This comparison is numerical reproduction of the same algorithm, not an independent algorithm implementation.'}
 a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 print(json.dumps({k:result[k] for k in ['status','numeric_fields_compared','other_scalar_fields_compared','difference_count','max_absolute_difference']}))
 if differences:raise SystemExit(1)
if __name__=='__main__':main()
