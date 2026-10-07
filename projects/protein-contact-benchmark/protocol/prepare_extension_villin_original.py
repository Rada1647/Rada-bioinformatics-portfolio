"""Apply the frozen metadata-only splits, then discovery-only contact selection.

Source CA arrays have already passed trajectory/provenance qualification.
This script never changes a split based on structural motions or contacts.
"""
from pathlib import Path
import argparse,datetime,hashlib,json,shutil
import numpy as np

ROOT=Path(__file__).resolve().parent
def digest(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def save(p,obj):p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def all_pairs(n):
 i,j=np.triu_indices(n,4);return np.column_stack([i,j]).astype(np.int64)
def contact_matrix(xyz,pairs,cutoff):
 X=np.empty((len(xyz),len(pairs)),dtype=bool)
 for a in range(0,len(xyz),128):
  block=xyz[a:a+128].astype(np.float64)
  d=np.linalg.norm(block[:,pairs[:,0]]-block[:,pairs[:,1]],axis=2)
  X[a:a+128]=d<float(cutoff)
 return X
def source_runs(dataset):
 base=ROOT/'data'/dataset
 if dataset=='ubiquitin':
  q=json.loads((base/'qualification.json').read_text())
  if not q.get('eligible',False):raise RuntimeError('Ubiquitin failed temporal qualification; no benchmark may be prepared.')
  raise RuntimeError('Eligible ubiquitin requires a documented mapping adapter; no assumption is implemented.')
 if dataset=='t4':
  rows=[]
  for i,name in enumerate(['md1us4','md1us5','md1us6','md1us7','md1us8']):
   xyz=base/(name+'_ca_200ps.npz')
   if not xyz.exists():raise RuntimeError('All five verified T4 CA runs are required: missing '+str(xyz))
   check=json.loads((base/(name+'.xtc.checksum.json')).read_text())
   assert check['published_checksum']=='md5:'+check['md5']
   assert check['published_size']==check['verified_size']
   rows.append({'id':name,'run_index':i,'split':('discovery' if i==0 else 'validation' if i==1 else 'test'), 'source':xyz, 'source_metadata':{'raw_xtc':check,'coordinate_sha256':digest(xyz)}})
  return rows,{'qualification':'Five published whole-file MD5s checked; see individual extraction/provenance records.'}
 manifest=json.loads((base/'trajectory_manifest.json').read_text())
 # Enumeration of the complete canonical manifest determines RNG run_index.
 manifest=sorted(manifest,key=lambda r:r['source_path'])
 allrows=[];excluded=[]
 for i,r in enumerate(manifest):
  ok=r['eligibility_min20_frames'] and r['geometry_pass'] and r['n_ca']==35 and abs(r['dt_ps']-200)<1e-3
  if not ok:excluded.append({'id':r['trajectory_id'],'reason':'metadata qualification failure','metadata':r});continue
  allrows.append({'id':r['trajectory_id'],'run_index':i,'group':str(r['group']),'source':base/r['npz_path'],'source_path':r['source_path'],'path_selection_sha256':hashlib.sha256(r['source_path'].encode()).hexdigest(),'source_metadata':r})
 selected=[]
 for g,split,cap in [('1','discovery',24),('2','validation',24),('3','test',48)]:
  group=sorted([r for r in allrows if r['group']==g],key=lambda r:r['path_selection_sha256'])
  if not group:raise RuntimeError('Empty eligible villin archive group '+g)
  for r in group[:cap]:r['split']=split;selected.append(r)
 return sorted(selected,key=lambda r:r['run_index']),{'complete_trajectory_count':len(manifest),'eligible_count':len(allrows),'excluded':excluded,'caps':{'discovery':24,'validation':24,'test':48},'group_split':{'1':'discovery','2':'validation','3':'test'},'independence':'Adaptive stress set; archive-group statistical independence unverified.'}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--dataset',choices=['t4','ubiquitin','villin'],required=True);ap.add_argument('--protocol',type=Path,default=ROOT/'protocol/protocol_v1.1.json');args=ap.parse_args()
 protocol=json.loads(args.protocol.read_text());assert protocol['version']=='1.1'
 out=ROOT/'data/prepared'/args.dataset
 if out.exists():raise SystemExit('Refusing to overwrite prepared input directory '+str(out))
 rows,qualification=source_runs(args.dataset)
 out.mkdir(parents=True);(out/'coordinates').mkdir();(out/'contacts').mkdir()
 manifest={'dataset':args.dataset,'dataset_id':{'t4':1,'ubiquitin':2,'villin':3}[args.dataset],'protocol_sha256':digest(args.protocol),'protocol_version':protocol['version'],'prepared_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'code_sha256':digest(__file__),'qualification':qualification,'runs':[]}
 nca=None
 for r in rows:
  with np.load(r['source']) as d:
   xyz=d['xyz'];t=d['time_ps']
   if nca is None:nca=xyz.shape[1]
   assert xyz.dtype==np.float32 and xyz.ndim==3 and xyz.shape[1:]==(nca,3)
   assert np.isfinite(xyz).all() and len(t)==len(xyz) and np.allclose(np.diff(t),200,atol=1e-3,rtol=0)
   assert len(t)>=20
  if args.dataset=='villin':assert digest(r['source'])==r['source_metadata']['npz_sha256']
  target=out/'coordinates'/(r['id']+'.npz');shutil.copyfile(r['source'],target)
  m={k:v for k,v in r.items() if k!='source'};m['coordinates_path']=str(target.relative_to(out));m['coordinates_sha256']=digest(target);m['n_frames']=len(t);m['n_ca']=nca;m['contacts']={}
  manifest['runs'].append(m)
 pairs=all_pairs(nca)
 coverage={}
 for cut in (7,8,9):
  sums=np.zeros(len(pairs),np.int64);total=0
  for r in manifest['runs']:
   if r['split']!='discovery':continue
   with np.load(out/r['coordinates_path']) as d:X=contact_matrix(d['xyz'],pairs,cut)
   sums+=X.sum(axis=0);total+=len(X)
  keep=(sums>0)&(sums<total)
  if not keep.any():raise RuntimeError(f'No discovery-variable contacts at cutoff {cut}; not evaluable')
  frequency=sums[keep]/total
  coverage[str(cut)]={'all_eligible_pair_count':len(pairs),'selected_pair_count':int(keep.sum()),'discovery_frame_count':total,'discovery_constant_absent_pairs':int((sums==0).sum()),'discovery_constant_present_pairs':int((sums==total).sum()),'runs':{}}
  for r in manifest['runs']:
   with np.load(out/r['coordinates_path']) as d:
    X=contact_matrix(d['xyz'],pairs,cut);times=d['time_ps']
   allpositive=int(X.sum());absentpositive=int(X[:,sums==0].sum());selectedpositive=int(X[:,keep].sum())
   coverage[str(cut)]['runs'][r['id']]={'split':r['split'],'all_positive_cells':allpositive,'selected_positive_cells':selectedpositive,'discovery_absent_positive_cells':absentpositive,'discovery_absent_positive_fraction':absentpositive/allpositive if allpositive else None,'omitted_positive_cells':allpositive-selectedpositive}
   dest=out/'contacts'/(r['id']+f'_cut{cut}.npz')
   np.savez_compressed(dest,X=X[:,keep],pairs=pairs[keep],frequency=frequency,time_ps=times)
   r['contacts'][str(cut)]=str(dest.relative_to(out));r.setdefault('contacts_sha256',{})[str(cut)]=digest(dest)
  print(args.dataset,'cutoff',cut,'selected',int(keep.sum()),'of',len(pairs),'eligible pairs',flush=True)
 manifest['coverage_path']='coverage.json'
 save(out/'coverage.json',coverage);save(out/'manifest.json',manifest)
 print('Prepared',len(rows),'trajectories;',out/'manifest.json',flush=True)
if __name__=='__main__':main()
