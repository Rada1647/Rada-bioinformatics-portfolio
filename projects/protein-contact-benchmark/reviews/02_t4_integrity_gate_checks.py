"""Exercise T4 completion gates in isolated temporary trees; production unchanged."""
from pathlib import Path
import copy, hashlib, importlib.util, json, shutil, tempfile
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('review_prepare_gate',ROOT/'prepare_extension.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
old=(ROOT/'protocol/prepare_extension_villin_original.py').read_text()
new=(ROOT/'prepare_extension.py').read_text()
marker=" manifest=json.loads((base/'trajectory_manifest.json').read_text())"
assert old[old.index(marker):]==new[new.index(marker):]
assert old[:old.index(" if dataset=='t4':")]==new[:new.index(" if dataset=='t4':")]
results=[]
def trial(name,mutate):
    with tempfile.TemporaryDirectory(prefix='review2_t4_gate_') as td:
        td=Path(td);base=td/'data/t4';base.mkdir(parents=True)
        shutil.copyfile(ROOT/'data/t4/zenodo_record.json',base/'zenodo_record.json')
        for run in ('md1us4','md1us5'):
            for suffix in ('.xtc.checksum.json','_extraction.json'):
                shutil.copyfile(ROOT/'data/t4'/(run+suffix),base/(run+suffix))
            (base/(run+'_ca_200ps.npz')).symlink_to(ROOT/'data/t4'/(run+'_ca_200ps.npz'))
        mutate(base)
        mod.ROOT=td
        try:mod.source_runs('t4')
        except Exception as e:
            message=str(e)
            if name=='control_complete4_and5_then_missing6':assert isinstance(e,RuntimeError) and 'md1us6_ca' in message
            else:assert not (isinstance(e,RuntimeError) and 'md1us6_ca' in message),name
            results.append({'test':name,'rejected':True,'exception_type':type(e).__name__,'message':message})
        else:raise AssertionError('Incomplete or inconsistent tree was accepted: '+name)
        assert not (td/'data/prepared').exists()
def change_json(base,file,key,value):
    p=base/file;j=json.loads(p.read_text());j[key]=value;p.write_text(json.dumps(j))
trial('control_complete4_and5_then_missing6',lambda b:None)
trial('missing_extraction_manifest',lambda b:(b/'md1us4_extraction.json').unlink())
for key,value in [('run','md1us5'),('output','wrong.npz'),('output_sha256','0'*64),('output_size',1),
                  ('n_ca',161),('sample_nframes',5000),('sample_dt_ps',100),('max_adjacent_ca_distance_angstrom',6)]:
    trial('bad_extraction_'+key,lambda b,k=key,v=value:change_json(b,'md1us4_extraction.json',k,v))
trial('checksum_record_wrong_published_md5',lambda b:change_json(b,'md1us4.xtc.checksum.json','md5','0'*32))
trial('checksum_record_wrong_source_size',lambda b:change_json(b,'md1us4.xtc.checksum.json','verified_size',1))
trial('checksum_record_wrong_sha256',lambda b:change_json(b,'md1us4.xtc.checksum.json','sha256','0'*64))
def partial_only(b):
    (b/'md1us4_ca_200ps.npz').unlink()
    (b/'md1us4.xtc.part').write_bytes(b'partial-not-a-trajectory')
    (b/'md1us4.xtc.part.ranges.json').write_text('{"completed": []}')
trial('partial_raw_and_checkpoint_without_CA',partial_only)
def checksum_self_consistent_but_not_published(b):
    p=b/'md1us4.xtc.checksum.json';d=json.loads(p.read_text());d['md5']='0'*32;d['published_checksum']='md5:'+d['md5'];p.write_text(json.dumps(d))
trial('self_consistent_checksum_not_published',checksum_self_consistent_but_not_published)
def corrupt_CA(b):
    p=b/'md1us4_ca_200ps.npz';p.unlink();p.write_bytes(b'incomplete zip bytes')
trial('truncated_CA',corrupt_CA)
mod.ROOT=ROOT
report={'reviewer':'Independent AI data/provenance audit #2','code_sha256':hashlib.sha256((ROOT/'prepare_extension.py').read_bytes()).hexdigest(),
  'only_T4_source_gate_changed':True,'tests':results,'all_rejected_before_prepared_output':True,
  'scope':'Isolated temporary copies of metadata with symlinked read-only completed4/5 coordinates; no production file changed. Control reaches missing6 as expected. Complete5-run positive gate remains pending.'}
(ROOT/'reviews/02_t4_integrity_gate_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print('Passed',len(results),'T4 integrity rejection/control cases; non-T4 source text unchanged.')
