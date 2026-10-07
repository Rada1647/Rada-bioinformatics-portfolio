#!/usr/bin/env python3
"""Retrieve published T4 runs, checksum, extract whole-protein CA at 200 ps.

Raw XTC files are ~10 GB each; after successful checksum/extraction only raw
files downloaded by this invocation are deleted. TPRs and provenance remain.
No contact features, benchmark results, or method metrics are calculated.

Public-release adaptation: completed records are skipped only after the
corresponding coordinate output passes its recorded size and SHA256 checks.
The scientific extraction operations are unchanged.
"""
from pathlib import Path
import argparse, concurrent.futures, datetime, hashlib, json, os, threading, time, urllib.error, urllib.request
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.mdamath import make_whole, triclinic_vectors

ROOT = Path(__file__).resolve().parent
DOWNLOAD_WORKERS=32
def stamp(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def writej(p, obj):
    tmp=p.with_suffix(p.suffix+'.tmp'); tmp.write_text(json.dumps(obj, indent=2)+'\n'); tmp.replace(p)
def filehash(p):
    hs={k:hashlib.new(k) for k in ['md5','sha256']}
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):
            for h in hs.values(): h.update(b)
    return {k:h.hexdigest() for k,h in hs.items()}
def ranged_download(rec,tmp):
    """Bounded range retrieval, then verify the complete original file."""
    total=rec['size']; checkpoint=tmp.with_suffix(tmp.suffix+'.ranges.json')
    if checkpoint.exists():
        state=json.loads(checkpoint.read_text()); done={tuple(x) for x in state['completed']};prefix=state['prefix_bytes']
    else:
        prefix=tmp.stat().st_size if tmp.exists() else 0;done=set()
        writej(checkpoint,{'prefix_bytes':prefix,'completed':[]})
    fd=os.open(tmp,os.O_CREAT|os.O_RDWR,0o644);os.ftruncate(fd,total)
    spans=[(i,min(i+64*1024**2,total)-1) for i in range(prefix,total,64*1024**2)]
    pending=[s for s in spans if s not in done]
    abort=threading.Event()
    def work(span):
        a,b=span
        url=f'https://zenodo.org/records/3989057/files/{rec["key"]}?download=1'
        for attempt in range(4):
            if abort.is_set(): raise RuntimeError('Download stopped after server rate limiting')
            try:
                req=urllib.request.Request(url,headers={'Range':f'bytes={a}-{b}','User-Agent':'protein-contact332-reproduction/1.0'})
                with urllib.request.urlopen(req,timeout=120) as r:
                    if r.status!=206 or r.headers.get('Content-Range')!=f'bytes {a}-{b}/{total}':raise RuntimeError('Server ignored or changed requested range')
                    pos=a
                    while pos<=b:
                        if abort.is_set(): raise RuntimeError('Download stopped after server rate limiting')
                        block=r.read(min(1024**2,b-pos+1))
                        if not block:raise RuntimeError('Truncated range')
                        wrote=os.pwrite(fd,block,pos)
                        if wrote!=len(block):raise RuntimeError('Short local write')
                        pos+=len(block)
                return span
            except Exception as exc:
                print(stamp(),'RANGE ERROR',rec['key'],a,b,'attempt',attempt+1,repr(exc),flush=True)
                if isinstance(exc,urllib.error.HTTPError) and exc.code in [429,503]:
                    abort.set();raise
                if abort.is_set():raise
                if attempt==3:raise
                time.sleep(2**attempt)
    start=time.time();last=start
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=DOWNLOAD_WORKERS) as pool:
            futures={pool.submit(work,span):span for span in pending}
            for future in concurrent.futures.as_completed(futures):
                span=future.result()
                done.add(span);writej(checkpoint,{'prefix_bytes':prefix,'completed':sorted(done)})
                if time.time()-last>=30:
                    n=prefix+sum(b-a+1 for a,b in done)
                    print(stamp(),rec['key'],f'{n/1e9:.3f}/{total/1e9:.3f} GB, {DOWNLOAD_WORKERS} ranges',flush=True);last=time.time()
        os.fsync(fd)
    finally:os.close(fd)
    hashes=filehash(tmp)
    if tmp.stat().st_size!=total or hashes['md5']!=rec['checksum'].split(':')[1]:raise RuntimeError(f'Complete-file checksum mismatch {hashes}')
    checkpoint.unlink()
    return hashes
def obtain(rec):
    p=ROOT/rec['key']; new=not p.exists()
    if new:
        tmp=p.with_suffix(p.suffix+'.part')
        if rec['size']>100000000:
            print(stamp(),'RANGED DOWNLOAD',rec['key'],rec['size'],flush=True)
            hashes=ranged_download(rec,tmp);tmp.replace(p)
            evidence={'name':p.name,'source_url':rec['links']['self'],'published_size':rec['size'],
                      'published_checksum':rec['checksum'],'verified_size':p.stat().st_size,
                      **hashes,'verified_utc':stamp(),'download_method':f'Up to {DOWNLOAD_WORKERS} concurrent HTTP ranges; complete original-file MD5 and SHA256'}
            writej(ROOT/(p.name+'.checksum.json'),evidence)
            return p,new,evidence
        if tmp.exists(): raise RuntimeError(f'Incomplete file exists: {tmp}; inspect before retry')
        start=time.time(); n=0; last=start
        print(stamp(),'DOWNLOAD',rec['key'],rec['size'],flush=True)
        hs={k:hashlib.new(k) for k in ['md5','sha256']}
        req=urllib.request.Request(rec['links']['self'],headers={'User-Agent':'protein-contact332-reproduction/1.0'})
        with urllib.request.urlopen(req,timeout=120) as r, tmp.open('wb') as out:
            while True:
                b=r.read(8*1024**2)
                if not b: break
                out.write(b); n+=len(b)
                for h in hs.values(): h.update(b)
                if time.time()-last>=30:
                    print(stamp(),rec['key'],f'{n/1e9:.3f}/{rec["size"]/1e9:.3f} GB',f'{n/(time.time()-start)/1e6:.1f} MB/s',flush=True);last=time.time()
        hashes={k:h.hexdigest() for k,h in hs.items()}
        if n!=rec['size'] or hashes['md5']!=rec['checksum'].split(':')[1]:
            raise RuntimeError(f'Checksum/size failure: {rec["key"]} {n} {hashes}')
        tmp.replace(p)
    else:
        hashes=filehash(p)
        if p.stat().st_size!=rec['size'] or hashes['md5']!=rec['checksum'].split(':')[1]:
            raise RuntimeError(f'Existing file fails checksum: {p}')
    evidence={'name':p.name,'source_url':rec['links']['self'],'published_size':rec['size'],
              'published_checksum':rec['checksum'],'verified_size':p.stat().st_size,
              **hashes,'verified_utc':stamp()}
    writej(ROOT/(p.name+'.checksum.json'),evidence)
    return p,new,evidence

def extract(run, tpr, xtc, sources):
    out=ROOT/(run+'_ca_200ps.npz')
    if out.exists(): raise RuntimeError(f'Refusing to overwrite {out}')
    whole=mda.Universe(str(tpr))
    protein=whole.select_atoms('protein')
    if not np.array_equal(protein.indices,np.arange(protein.n_atoms)):
        raise ValueError('Protein topology is not the first contiguous atom block')
    whole.load_new(np.zeros((1,whole.atoms.n_atoms,3),np.float32))
    u=mda.Merge(protein)
    u.load_new(str(xtc))
    if u.atoms.n_atoms!=protein.n_atoms: raise ValueError('Protein-only XTC atom count mismatch')
    ca=u.select_atoms('name CA')
    backbone=u.select_atoms('name N CA C')
    if len(ca)!=162 or not np.array_equal(ca.resids,np.arange(1,163)):
        raise ValueError('Unexpected T4 CA topology')
    if len(u.atoms.fragments)!=1: raise ValueError('Expected one covalent protein fragment')
    source_nframes=len(u.trajectory); dt=float(u.trajectory.dt)
    stride=round(200.0/dt)
    if not np.isclose(stride*dt,200.0): raise ValueError(f'Cannot exactly sample 200 ps from {dt}')
    frame_indices=np.arange(0,source_nframes,stride,dtype=np.int64)
    xyz=np.empty((len(frame_indices),len(ca),3),np.float32)
    times=np.empty(len(frame_indices),np.float64)
    box=np.empty((len(frame_indices),6),np.float32)
    max_ca_bond=0.; max_shift=0.
    bond_indices=u.atoms.bonds.to_indices()
    full_unwrap_calls=0;already_whole_frames=0
    print(stamp(),'EXTRACT START',run,'source_frames',source_nframes,'selected_frames',len(frame_indices),flush=True)
    for k, frame in enumerate(frame_indices):
        ts=u.trajectory[int(frame)]
        before=ca.positions.copy()
        bonds=u.atoms.positions[bond_indices[:,0]]-u.atoms.positions[bond_indices[:,1]]
        shortest_box_bound=np.linalg.svd(triclinic_vectors(ts.dimensions),compute_uv=False)[-1]
        if shortest_box_bound>6.0 and np.all(np.einsum('ij,ij->i',bonds,bonds)<9.0):
            already_whole_frames+=1
        else:
            make_whole(backbone,inplace=True);full_unwrap_calls+=1
        xyz[k]=ca.positions;times[k]=ts.time;box[k]=ts.dimensions
        max_ca_bond=max(max_ca_bond,float(np.linalg.norm(np.diff(xyz[k],axis=0),axis=1).max()))
        max_shift=max(max_shift,float(np.abs(before-xyz[k]).max()))
        if (k+1)%1000==0: print(stamp(),'EXTRACT PROGRESS',run,k+1,'/',len(frame_indices),flush=True)
    if not np.all(np.isfinite(xyz)) or not np.allclose(np.diff(times),200.0):
        raise ValueError('Invalid coordinates or sample times')
    if max_ca_bond>5.0: raise ValueError(f'Adjacent CA separation >5 A: {max_ca_bond}')
    np.savez_compressed(out, xyz=xyz, time_ps=times, ca_atom_indices=ca.indices.astype(np.int64),
                        resids=ca.resids.astype(np.int64),resnames=np.asarray(ca.resnames,dtype='U3'),
                        frame_indices=frame_indices,box_dimensions=box)
    u.trajectory.close()
    whole.trajectory.close()
    record={'run':run,'created_utc':stamp(),'source_record':'https://zenodo.org/records/3989057',
            'source_files':sources,'source_total_topology_atoms':whole.atoms.n_atoms,
            'source_xtc_protein_atoms':u.atoms.n_atoms,'n_ca':len(ca),'source_nframes':source_nframes,
            'source_dt_ps':dt,'sampling_stride':stride,'sample_dt_ps':200.,'sample_nframes':len(times),
            'first_time_ps':float(times[0]),'last_time_ps':float(times[-1]),
            'ca_atom_indices_0based':ca.indices.tolist(),'resids':ca.resids.tolist(),'resnames':ca.resnames.tolist(),
            'pbc':'Bonded MDAnalysis make_whole on connected N-CA-C backbone when needed. Fast path: all raw protein covalent bond lengths <3 A and smallest singular value of triclinic box matrix >6 A; otherwise call make_whole(backbone). Topology bonds retained from TPR; first N is the same reference as full protein.',
            'backbone_make_whole_frames':full_unwrap_calls,'already_whole_frames':already_whole_frames,
            'max_adjacent_ca_distance_angstrom':max_ca_bond,'max_unwrap_coordinate_shift_angstrom':max_shift,
            'units':{'xyz':'angstrom','time_ps':'picosecond','box_dimensions':'angstrom,degree'},
            'output':out.name,'output_size':out.stat().st_size,'output_sha256':filehash(out)['sha256'],
            'MDAnalysis_version':mda.__version__,'numpy_version':np.__version__,
            'license_note':'Zenodo record is open access but metadata contains no explicit data license. Do not assume a redistribution license.'}
    writej(ROOT/(run+'_extraction.json'),record)
    print(stamp(),'EXTRACTED',run,xyz.shape,'dt',dt,'max_CA_bond_A',max_ca_bond,flush=True)
    return record

def main():
    global DOWNLOAD_WORKERS
    ap=argparse.ArgumentParser();ap.add_argument('--runs',nargs='+',default=['md1us4','md1us5','md1us6','md1us7','md1us8']);ap.add_argument('--keep-raw',action='store_true');ap.add_argument('--workers',type=int,default=32);args=ap.parse_args()
    if not 1<=args.workers<=32:raise ValueError('workers must be between 1 and 32')
    DOWNLOAD_WORKERS=args.workers
    meta=json.loads((ROOT/'zenodo_record.json').read_text()); files={r['key']:r for r in meta['files']}
    for run in args.runs:
        completion=ROOT/(run+'_extraction.json')
        if completion.exists():
            saved=json.loads(completion.read_text())
            expected=ROOT/(run+'_ca_200ps.npz')
            if (saved.get('run')!=run or saved.get('output')!=expected.name
                or not expected.is_file() or expected.stat().st_size!=saved.get('output_size')
                or filehash(expected)['sha256']!=saved.get('output_sha256')):
                raise RuntimeError(f'Completion record does not verify existing coordinates: {completion}. Preserve this record for inspection and reconstruct in a fresh directory; no missing data will be silently skipped.')
            print('SKIP verified completed',run,flush=True);continue
        tpr,_,te=obtain(files[run+'.tpr'])
        xtc,new,xe=obtain(files[run+'.xtc'])
        record=extract(run,tpr,xtc,[te,xe])
        if new and not args.keep_raw and run!='md1us8':
            xtc.unlink()
            record['raw_xtc_removed_after_verified_extraction']=True
            writej(ROOT/(run+'_extraction.json'),record)
            print(stamp(),'REMOVED newly downloaded',xtc.name,flush=True)

if __name__=='__main__':main()
