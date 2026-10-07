"""Independent AI data audit; imports no production benchmark functions."""
from pathlib import Path
import collections, datetime, hashlib, json, re, struct, tempfile, warnings, zipfile
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import minimize_vectors

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'reviews'
def hash_file(p, algo='sha256'):
    h = hashlib.new(algo)
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''): h.update(b)
    return h.hexdigest()
def read(p): return json.loads((ROOT / p).read_text())

def main():
    warnings.filterwarnings('ignore', category=UserWarning)
    warnings.filterwarnings('ignore', category=DeprecationWarning)
    out = {'reviewer':'Independent AI technical auditor #2, data/provenance',
           'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
           'software':{'MDAnalysis':mda.__version__, 'numpy':np.__version__},
           'scope':'Initial Villin and ubiquitin data audit; T4 extension pending',
           'code_and_protocol_sha256':{str(p.relative_to(ROOT)):hash_file(p) for p in
              [ROOT/'prepare_extension.py',ROOT/'data/villin/qualify_extract.py',
               ROOT/'data/ubiquitin/inspect_archive.py',ROOT/'protocol/protocol_v1.1.json',Path(__file__)]}}
    manifest = read('data/prepared/villin/manifest.json')
    source = read('data/villin/trajectory_manifest.json')
    source_by_path = {x['source_path']:x for x in source}
    prepared_by_path = {x['source_path']:x for x in manifest['runs']}
    prepared = ROOT/'data/prepared/villin'
    v = {}
    archive_path = ROOT/'data/villin/protein_folding_datasets.zip'
    fig = read('data/villin/figshare_metadata.json')
    expected = next(x for x in fig['files'] if x['name']==archive_path.name)
    v['archive']={'bytes':archive_path.stat().st_size,'md5':hash_file(archive_path,'md5'),'sha256':hash_file(archive_path)}
    assert v['archive']['md5']==expected['computed_md5'] and v['archive']['bytes']==expected['size']
    assert v['archive']['sha256']==read('data/villin/qualification.json')['archive_sha256']
    assert manifest['protocol_sha256']==out['code_and_protocol_sha256']['protocol/protocol_v1.1.json']
    archived_prepare=ROOT/'protocol/prepare_extension_villin_original.py'
    if manifest['code_sha256']!=out['code_and_protocol_sha256']['prepare_extension.py']:
        assert archived_prepare.exists() and hash_file(archived_prepare)==manifest['code_sha256']
        out['preparation_manifest_code_note']='Villin manifest matches archived original preparation code; current preparation adds T4 integrity gates.'
        out['code_and_protocol_sha256']['protocol/prepare_extension_villin_original.py']=hash_file(archived_prepare)
    coords={}; topologies={}; raw_checks=[]
    with zipfile.ZipFile(archive_path) as z, tempfile.TemporaryDirectory(prefix='review2_villin_') as td:
        td=Path(td)
        names=sorted(n for n in z.namelist() if n.endswith('.xtc'))
        assert len(names)==len(set(names))==len(source)==2137
        assert set(names)==set(source_by_path)
        groups=collections.Counter(n.split('/')[1] for n in names)
        chosen=[]
        for group,split,cap in [('1','discovery',24),('2','validation',24),('3','test',48)]:
            eligible=[]
            for n in names:
                r=source_by_path[n]
                if n.split('/')[1]==group and r['eligibility_min20_frames'] and r['geometry_pass'] and r['n_ca']==35 and abs(r['dt_ps']-200)<1e-3:
                    eligible.append(n)
            ranked=sorted(eligible,key=lambda n:hashlib.sha256(n.encode()).hexdigest())[:cap]
            chosen.extend(ranked)
            assert set(ranked)==set(r['source_path'] for r in manifest['runs'] if r['split']==split)
        assert set(chosen)==set(prepared_by_path) and len(chosen)==96
        full_run_index={n:i for i,n in enumerate(names)}
        assert all(r['run_index']==full_run_index[n] for n,r in prepared_by_path.items())
        v['selection']={'archive_group_counts':dict(groups),'selected_split_counts':dict(collections.Counter(r['split'] for r in manifest['runs'])),
            'recomputed_from_archive_paths':True,'manifest_eligible_count':len(source),'selection_mismatches':0,
            'source_cadence_min_ps':min(r['source_dt_min_ps'] for r in source),
            'source_cadence_max_ps':max(r['source_dt_max_ps'] for r in source)}
        assert all(r['source_dt_min_ps']==100 and r['source_dt_max_ps']==100 for r in source)
        nodes={r['child_node']:r for r in source}
        roots=collections.Counter()
        for r in source:
            assert r['group']==r['source_path'].split('/')[1]
            seen=set();n=r['child_node']
            while nodes[n]['parent_node'] is not None:
                assert n not in seen;seen.add(n)
                child=nodes[n];parent=nodes[child['parent_node']]
                assert child['group']==parent['group'] and child['epoch']>parent['epoch']
                assert 0<=child['parent_frame']<parent['n_source_frames']
                n=child['parent_node']
            roots[n]+=1
        v['ancestry']={'root_count':len(roots),'parent_links':sum(r['parent_node'] is not None for r in source),
            'roots_per_group':dict(collections.Counter(k.split('/')[0] for k in roots)),
            'missing_cross_group_cyclic_or_out_of_bounds_links':0,
            'independence':'Not established by metadata; adaptive stress reporting remains necessary.'}
        for g in ('1','2','3'):
            topname=f'protein_folding_datasets/{g}/filtered/filtered.pdb'
            raw=z.read(topname);p=td/f'g{g}.pdb';p.write_bytes(raw)
            u=mda.Universe(str(p)); ca=u.select_atoms('name CA')
            t={'n_atoms':len(u.atoms),'n_ca':len(ca),'resids':ca.resids.tolist(), 'resnames':ca.resnames.tolist(),
               'atom_indices':ca.indices.tolist(),'chain_ids':sorted(set(ca.chainIDs)), 'segids':sorted(set(ca.segids)),
               'sha256':hashlib.sha256(raw).hexdigest()}
            assert t['n_atoms']==582 and t['n_ca']==35 and t['resids']==list(range(42,77))
            assert [i for i,n in zip(t['resids'],t['resnames']) if n=='NLE']==[65,70]
            assert [i for i,n in zip(t['resids'],t['resnames']) if n=='HSP']==[68]
            topologies[g]=t
        assert len(set(t['sha256'] for t in topologies.values()))==1
        v['topologies']=topologies
        for k,n in enumerate(sorted(chosen)):
            r=prepared_by_path[n];g=r['group'];raw=z.read(n);xp=td/'sample.xtc';xp.write_bytes(raw)
            assert hashlib.sha256(raw).hexdigest()==r['source_metadata']['source_xtc_sha256']
            u=mda.Universe(str(td/f'g{g}.pdb'),str(xp));ca=u.select_atoms('name CA')
            times=np.array([ts.time for ts in u.trajectory]); frames=np.arange(0,len(times),2)
            assert np.all(np.diff(times)==100)
            arr=np.stack([ca.positions.copy() for ts in u.trajectory[::2]])
            boxes=np.stack([ts.dimensions.copy() for ts in u.trajectory[::2]])
            cp=prepared/r['coordinates_path'];assert hash_file(cp)==r['coordinates_sha256']
            with np.load(cp) as d:
                err=float(np.max(np.abs(arr-d['xyz'])))
                assert err==0 and np.array_equal(d['time_ps'],times[frames])
                assert np.array_equal(d['source_frame'],frames)
                assert np.array_equal(d['ca_indices'],ca.indices) and np.array_equal(d['ca_resnames'],ca.resnames)
            assert np.all(np.diff(times[frames])==200) and np.isfinite(arr).all()
            step=np.diff(arr.astype('float64'),axis=1)
            adj=np.linalg.norm(step,axis=2)
            pbc_error=max(float(np.max(np.abs(minimize_vectors(st,bx)-st))) for st,bx in zip(step,boxes))
            assert adj.max()<5 and pbc_error<1e-4
            ur=mda.Universe(str(td/f'g{g}.pdb'),str(xp),convert_units=False)
            native=ur.select_atoms('name CA').positions.copy()
            uniterr=float(np.max(np.abs(native*10-arr[0])))
            assert uniterr==0
            coords[n]=arr
            raw_checks.append({'id':r['id'],'n_source_frames':len(times),'n_sampled_frames':len(frames),
                'xyz_max_abs_error_A':err,'time_exact':True,'angstrom_vs_native_nm_max_abs_error_A':uniterr,
                'adjacent_ca_min_A':float(adj.min()),'adjacent_ca_max_A':float(adj.max()),'adjacent_minimum_image_max_abs_correction_A':pbc_error})
            u.trajectory.close();ur.trajectory.close()
            if k%24==0:print('Independent raw trajectories checked:',k+1,flush=True)
    v['raw_reextraction_checks']=raw_checks
    v['raw_coordinate_values_checked']=sum(a.size for a in coords.values())
    # Pair enumeration and squared-distance threshold independently implemented.
    pairs=np.array([(i,j) for i in range(35) for j in range(i+4,35)],dtype=int)
    sq={n:((a[:,pairs[:,0]].astype('float64')-a[:,pairs[:,1]].astype('float64'))**2).sum(axis=2) for n,a in coords.items()}
    v['contacts']={}
    cov=read('data/prepared/villin/coverage.json')
    for cutoff in (7,8,9):
        discovery=[sq[r['source_path']]<cutoff**2 for r in manifest['runs'] if r['split']=='discovery']
        counts=sum(x.sum(axis=0) for x in discovery);total=sum(len(x) for x in discovery)
        mask=(counts!=0)&(counts!=total);freq=counts[mask]/total
        checked=0
        for n,r in prepared_by_path.items():
            cp=prepared/r['contacts'][str(cutoff)];assert hash_file(cp)==r['contacts_sha256'][str(cutoff)]
            ref=sq[n]<cutoff**2
            with np.load(cp) as d:
                assert np.array_equal(d['pairs'],pairs[mask]) and np.array_equal(d['frequency'],freq)
                assert np.array_equal(d['X'],ref[:,mask])
                checked+=d['X'].size
            cr=cov[str(cutoff)]['runs'][r['id']]
            assert cr['all_positive_cells']==int(ref.sum())
            assert cr['selected_positive_cells']==int(ref[:,mask].sum())
            assert cr['discovery_absent_positive_cells']==int(ref[:,counts==0].sum())
        v['contacts'][str(cutoff)]={'eligible_pair_count':len(pairs),'selected_pair_count':int(mask.sum()),
            'discovery_frame_count':total,'constant_absent_count':int((counts==0).sum()),'constant_present_count':int((counts==total).sum()),
            'contact_cells_checked':checked,'contact_pair_frequency_or_coverage_mismatches':0}
    out['villin']=v
    b=ROOT/'data/ubiquitin';up=b/'ubiquitin-md-generated-ensemble.zip'
    ub={'archive_bytes':up.stat().st_size,'archive_md5':hash_file(up,'md5'),'archive_sha256':hash_file(up)}
    uq=read('data/ubiquitin/qualification.json')
    assert ub['archive_md5']==uq['archive_md5'] and ub['archive_sha256']==uq['archive_sha256']
    zen=read('data/ubiquitin/zenodo_record.json')
    assert next(f['checksum'] for f in zen['files'] if f['key']==up.name)=='md5:'+ub['archive_md5']
    with zipfile.ZipFile(up) as z:
        ub['archive_members']=z.namelist();assert z.testzip() is None
        ub['all_member_checksums_match_extracted']=all(hashlib.sha256(z.read(n)).hexdigest()==hash_file(b/'archive'/n) for n in z.namelist() if not n.endswith('/'))
    u=mda.Universe(str(b/'archive/topology.pdb'),str(b/'archive/Q99.dcd'))
    ca=u.select_atoms('name CA')
    aa3='ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split()
    aa=dict(zip(aa3,'ARNDCQEGHILKMFPSTWYV'))
    ub['topology']={'n_atoms':len(u.atoms),'n_residues':len(u.residues),'n_ca':len(ca),'sequence':''.join(aa[n] for n in ca.resnames)}
    ub['q99']={'n_frames':len(u.trajectory),'reader_dt_ps':u.trajectory.dt,
               'first_time_ps':u.trajectory[0].time,'last_time_ps':u.trajectory[-1].time}
    with (b/'archive/Q99.dcd').open('rb') as f:
        raw=f.read(1024)
    endian='<' if struct.unpack('<i',raw[:4])[0]==84 else '>'
    ctl=struct.unpack(endian+'20i',raw[8:88])
    ub['q99']['header']={'nset':ctl[0],'istart':ctl[1],'nsavc':ctl[2],'delta_float':struct.unpack(endian+'f',raw[44:48])[0],
      'titles':[raw[100+i*80:180+i*80].split(b'\x00')[0].decode('ascii',errors='replace').strip()
                for i in range(struct.unpack(endian+'i',raw[96:100])[0])]}
    u.trajectory.close()
    out['ubiquitin']=ub
    (OUT/'02_data_numerical_checks.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'villin_checked_raw':len(raw_checks),'coordinates':v['raw_coordinate_values_checked'],
        'contacts':v['contacts'],'ubiquitin':ub},indent=2),flush=True)

if __name__=='__main__':main()
