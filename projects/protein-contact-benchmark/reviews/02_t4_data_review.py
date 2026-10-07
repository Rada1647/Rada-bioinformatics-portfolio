"""Independent reviewer #2 T4 source/prepared audit. No production imports."""
from pathlib import Path
import argparse, collections, datetime, hashlib, json, warnings
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.mdamath import make_whole

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/t4'
def js(p):return json.loads(p.read_text())
def hashes(p):
    h={a:hashlib.new(a) for a in ('md5','sha256')}
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):
            for x in h.values():x.update(b)
    return {a:x.hexdigest() for a,x in h.items()}
def metadata():
    meta=js(DATA/'zenodo_record.json');pub={f['key']:f for f in meta['files']}
    result={'reviewer':'Independent AI data/provenance audit #2',
      'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'code_sha256':hashes(Path(__file__))['sha256'],
      'acquisition_code_sha256':hashes(DATA/'acquire_extract.py')['sha256'],
      'preparation_code_sha256':hashes(ROOT/'prepare_extension.py')['sha256'],
      'protocol_sha256':hashes(ROOT/'protocol/protocol_v1.1.json')['sha256'],
      'topologies':{},'completed_outputs':{}}
    expected=None
    aa=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))
    for i in range(4,9):
        name=f'md1us{i}';tp=DATA/(name+'.tpr');hh=hashes(tp)
        assert hh['md5']==pub[tp.name]['checksum'].split(':')[1]
        assert tp.stat().st_size==pub[tp.name]['size']
        u=mda.Universe(str(tp));protein=u.select_atoms('protein');ca=protein.select_atoms('name CA')
        assert np.array_equal(protein.indices,np.arange(2612)) and len(u.atoms)==51932
        assert len(ca)==162 and np.array_equal(ca.resids,np.arange(1,163)) and len(protein.fragments)==1
        bonds=protein.bonds.to_indices();seq=''.join(aa[r] for r in ca.resnames)
        seq_audit=js(DATA/'sequence_audit.json');assert seq==seq_audit['observed_sequence']
        mapping=(ca.indices.tolist(),ca.resids.tolist(),ca.resnames.tolist(),protein.names.tolist(),bonds.tolist())
        if expected is None:expected=mapping
        assert mapping==expected
        bb=protein.select_atoms('name N CA C');bbset=set(bb.indices.tolist())
        links=collections.defaultdict(list)
        for a,b in bonds:
            if int(a) in bbset and int(b) in bbset:links[int(a)].append(int(b));links[int(b)].append(int(a))
        seen=set();stack=[int(bb.indices[0])]
        while stack:
            n=stack.pop()
            if n in seen:continue
            seen.add(n);stack.extend(links[n])
        assert seen==bbset and len(bb)==486
        result['topologies'][name]={'md5':hh['md5'],'sha256':hh['sha256'],'bytes':tp.stat().st_size,
          'total_atoms':len(u.atoms),'protein_atoms':len(protein),'ca_count':len(ca),
          'sequence':seq,'residue_start':int(ca.resids[0]),'residue_end':int(ca.resids[-1]),
          'connected_backbone_atoms':len(bb),'bonds_sha256':hashlib.sha256(bonds.tobytes()).hexdigest()}
        reference_sequence=''.join(line.strip() for line in (DATA/'reference_P00720.fasta').read_text().splitlines() if not line.startswith('>'))
        result['construct_reference_comparison']={'observed_sequence':seq,'observed_length':len(seq),
          'reference_sequence':reference_sequence,'reference_length':len(reference_sequence),
          'aligned_mismatches':[{'position':k+1,'reference':a,'observed':b} for k,(a,b) in enumerate(zip(reference_sequence,seq)) if a!=b],
          'reference_extra_C_terminal_residues':reference_sequence[len(seq):],
          'interpretation':'Exact deposited162-residue construct; aligned differences do not by themselves establish engineered versus strain/reference substitutions.'}
        ep=DATA/(name+'_extraction.json')
        if not ep.exists():continue
        e=js(ep);npz=DATA/e['output'];oh=hashes(npz)
        assert oh['sha256']==e['output_sha256'] and npz.stat().st_size==e['output_size']
        for evidence in e['source_files']:
            original=pub[evidence['name']]
            assert evidence['published_size']==evidence['verified_size']==original['size']
            assert evidence['published_checksum']=='md5:'+evidence['md5']==original['checksum']
            separate=js(DATA/(evidence['name']+'.checksum.json'))
            assert evidence==separate
        with np.load(npz) as d:
            assert d['xyz'].shape==(5001,162,3) and d['xyz'].dtype==np.float32 and np.isfinite(d['xyz']).all()
            assert np.array_equal(d['time_ps'],np.arange(5001,dtype=np.float64)*200)
            assert np.array_equal(d['frame_indices'],np.arange(5001,dtype=np.int64)*200)
            assert np.array_equal(d['ca_atom_indices'],ca.indices) and np.array_equal(d['resnames'],ca.resnames)
            assert np.array_equal(d['resids'],ca.resids)
            adj=np.sqrt(np.sum(np.diff(d['xyz'].astype(np.float64),axis=1)**2,axis=2))
            assert adj.max()<5
            assert e['source_nframes']==1000001 and e['source_dt_ps']==1 and e['sampling_stride']==200
            assert e['sample_nframes']==5001 and e['sample_dt_ps']==200
            assert e['already_whole_frames']+e['backbone_make_whole_frames']==5001
            result['completed_outputs'][name]={'sha256':oh['sha256'],'n_frames':5001,'n_ca':162,
              'times_exact_0_to_1000000_ps_step200':True,'min_adjacent_CA_A':float(adj.min()),'max_adjacent_CA_A':float(adj.max()),
              'all_CA_scalars_finite':True,'raw_published_md5_manifest_verified':True,
              'raw_rehash_performed':False,'coordinate_hash_matches_extraction_manifest':True}
    reference=DATA/'md1us4_full_make_whole_reference.npz'
    if reference.exists():
        with np.load(reference) as a,np.load(DATA/'md1us4_ca_200ps.npz') as b:
            assert set(a.files)==set(b.files)
            assert all(np.array_equal(a[k],b[k]) for k in a.files)
        result['md1us4_saved_full_vs_backbone_reference']={'all_fields_bitwise_equal':True,'scalar_count':5001*162*3,
          'scope':'Independent comparison of saved arrays; raw4 removed so not independent raw re-extraction.'}
    result['all_five_complete']=len(result['completed_outputs'])==5
    result['prepared_exists']=(ROOT/'data/prepared/t4/manifest.json').exists()
    return result

def raw8(result):
    name='md1us8';p=DATA/(name+'.xtc');npz=DATA/(name+'_ca_200ps.npz')
    assert p.exists() and npz.exists() and (DATA/(name+'_extraction.json')).exists()
    hh=hashes(p);e=js(DATA/(name+'.xtc.checksum.json'))
    assert hh['md5']==e['md5'] and hh['sha256']==e['sha256'] and p.stat().st_size==e['published_size']
    top=mda.Universe(str(DATA/(name+'.tpr')))
    top.load_new(np.zeros((1,len(top.atoms),3),np.float32))
    u=mda.Merge(top.select_atoms('protein'));u.load_new(str(p));ca=u.select_atoms('name CA')
    assert len(u.trajectory)==1000001
    raw_first_CA=ca.positions.copy()
    native_reader=mda.coordinates.XTC.XTCReader(str(p),convert_units=False)
    native_first_CA=native_reader[0].positions[ca.indices].copy()
    native_reader.close()
    unit_error=float(np.max(np.abs(native_first_CA*10-raw_first_CA)))
    assert unit_error==0
    # Decode every XTC header's magic/atom-count/time, using reader offsets only
    # for random-access positioning. This verifies native cadence without
    # decompressing a million protein coordinate frames.
    offsets=np.asarray(u.trajectory._xdr.offsets,dtype=np.int64)
    assert len(offsets)==1000001 and np.all(np.diff(offsets)>0)
    raw=np.memmap(p,mode='r',dtype=np.uint8)
    def header_values(byte_offset,dtype):
        b=raw[offsets[:,None]+byte_offset+np.arange(4)[None,:]].copy()
        return b.reshape(-1).view(dtype)
    assert np.all(header_values(0,'>i4')==1995)
    assert np.all(header_values(4,'>i4')==2612)
    alltimes=header_values(12,'>f4').astype(np.float64)
    assert np.array_equal(alltimes,np.arange(1000001,dtype=np.float64))
    sampled=np.unique(np.rint(np.linspace(0,5000,129)).astype(int))
    cells=0;maxerr=0.;rawshift=0.;contact_cells=0
    sample_distances={}
    pairs=np.array([(i,j) for i in range(162) for j in range(i+4,162)],dtype=int)
    with np.load(npz) as d:
        for k in sampled:
            ts=u.trajectory[int(k*200)];before=ca.positions.copy()
            make_whole(u.atoms,inplace=True)
            after=ca.positions.copy();diff=float(np.max(np.abs(after-d['xyz'][k])))
            assert diff==0 and ts.time==d['time_ps'][k]
            cells+=after.size;maxerr=max(maxerr,diff);rawshift=max(rawshift,float(np.max(np.abs(after-before))))
            q=after.astype(np.float64);sq=((q[pairs[:,0]]-q[pairs[:,1]])**2).sum(axis=1)
            sample_distances[int(k)]=sq
    result['raw_md1us8']={'hashes':hh,'bytes':p.stat().st_size,'raw_source_headers_checked':len(offsets),
      'native_source_cadence_exact_ps':1,'raw_source_time_range_ps':[float(alltimes[0]),float(alltimes[-1])],
      'native_nm_to_angstrom_max_error_A':unit_error,
      'independent_full_protein_make_whole_sample_frames':sampled.tolist(),'ca_scalars_checked':cells,
      'max_abs_CA_error_A':maxerr,'max_raw_to_whole_CA_shift_A':rawshift,
      'contact_sample_check':'Pending prepared contacts' }
    prep=ROOT/'data/prepared/t4';mp=prep/'manifest.json'
    if mp.exists():
        pr=next(r for r in js(mp)['runs'] if r['id']==name)
        for cutoff in (7,8,9):
            with np.load(prep/pr['contacts'][str(cutoff)]) as d:
                pairlookup={tuple(p):i for i,p in enumerate(pairs)}
                cols=np.array([pairlookup[tuple(p)] for p in d['pairs']])
                for k in sampled:
                    ref=sample_distances[int(k)][cols]<cutoff**2
                    assert np.array_equal(ref,d['X'][k]);contact_cells+=len(ref)
        result['raw_md1us8']['contact_sample_check']={'contact_cells_checked':contact_cells,'mismatches':0}
    del raw;u.trajectory.close();top.trajectory.close()

def prepared(result):
    out=ROOT/'data/prepared/t4';manifest=js(out/'manifest.json')
    assert result['all_five_complete']
    assert manifest['protocol_sha256']==result['protocol_sha256']
    expected=[(f'md1us{i}',i-4,'discovery' if i==4 else 'validation' if i==5 else 'test') for i in range(4,9)]
    assert [(r['id'],r['run_index'],r['split']) for r in manifest['runs']]==expected
    allpairs=np.array([(i,j) for i in range(162) for j in range(i+4,162)],dtype=np.int64)
    assert len(allpairs)==12561
    coverage=js(out/'coverage.json');summary={}
    for r in manifest['runs']:
        hh=hashes(out/r['coordinates_path'])['sha256']
        assert hh==r['coordinates_sha256']==result['completed_outputs'][r['id']]['sha256']
    for cut in (7,8,9):
        r=manifest['runs'][0]
        with np.load(out/r['coordinates_path']) as d:
            a=d['xyz'].astype(np.float64)
        counts=np.zeros(len(allpairs),np.int64)
        for k in range(0,len(a),32):
            x=a[k:k+32];sq=((x[:,allpairs[:,0]]-x[:,allpairs[:,1]])**2).sum(axis=2)
            counts+=(sq<cut*cut).sum(axis=0)
        mask=(counts>0)&(counts<len(a));freq=counts[mask]/len(a)
        checked=0
        for r in manifest['runs']:
            cp=out/r['contacts'][str(cut)];assert hashes(cp)['sha256']==r['contacts_sha256'][str(cut)]
            allpos=selpos=abspos=0
            with np.load(out/r['coordinates_path']) as d,np.load(cp) as c:
                assert np.array_equal(c['pairs'],allpairs[mask]) and np.array_equal(c['frequency'],freq)
                assert np.array_equal(c['time_ps'],d['time_ps'])
                for k in range(0,len(d['xyz']),32):
                    x=d['xyz'][k:k+32].astype(np.float64)
                    ref=(((x[:,allpairs[:,0]]-x[:,allpairs[:,1]])**2).sum(axis=2)<cut*cut)
                    assert np.array_equal(ref[:,mask],c['X'][k:k+32])
                    checked+=ref[:,mask].size;allpos+=int(ref.sum());selpos+=int(ref[:,mask].sum());abspos+=int(ref[:,counts==0].sum())
            cr=coverage[str(cut)]['runs'][r['id']]
            assert [cr['all_positive_cells'],cr['selected_positive_cells'],cr['discovery_absent_positive_cells']]==[allpos,selpos,abspos]
        summary[str(cut)]={'eligible_pairs':len(allpairs),'selected_pairs':int(mask.sum()),'discovery_frames':len(a),'contact_cells_checked':checked,'mismatches':0}
    result['prepared_contact_checks']=summary

def main():
    warnings.filterwarnings('ignore',category=UserWarning)
    ap=argparse.ArgumentParser();ap.add_argument('--raw8',action='store_true');ap.add_argument('--prepared',action='store_true');args=ap.parse_args()
    result=metadata()
    if args.prepared:prepared(result)
    if args.raw8:raw8(result)
    filename='02_t4_final_numerical_checks.json' if args.prepared and args.raw8 else '02_t4_progress_checks.json'
    (ROOT/'reviews'/filename).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'completed':list(result['completed_outputs']),'all_five_complete':result['all_five_complete'],
      'prepared':result['prepared_exists'],'raw8_checked':'raw_md1us8' in result,'output':filename},indent=2),flush=True)
if __name__=='__main__':main()
