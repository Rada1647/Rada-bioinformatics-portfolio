"""Wait for complete original T4 outputs, then execute fresh-environment replay."""
from pathlib import Path
import datetime,json,os,subprocess,sys,time
R=Path(__file__).resolve().parents[1]
protocol=json.loads((R/'protocol/protocol_v1.1.json').read_text())
names=[s['name']+'_'+phase+'.json' for s in protocol['scenarios'] for phase in ['validation','test']]+['primary_clean.json']
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
print(stamp(),'WAITING for 17 complete original T4 JSONs',flush=True)
while True:
 complete=0
 for name in names:
  try:
   d=json.loads((R/'results/t4'/name).read_text())
   if d.get('status')=='complete' and d.get('dataset')=='t4':complete+=1
  except (OSError,json.JSONDecodeError):pass
 if complete==17:break
 print(stamp(),'ORIGINAL_PROGRESS',complete,'/',len(names),flush=True)
 time.sleep(30)
print(stamp(),'ORIGINAL_COMPLETE; REPRODUCTION_START',flush=True)
out=R/'results/reproduced_t4'
if out.exists() and any(out.iterdir()):raise RuntimeError('Refusing to overwrite any T4 reproduction outputs')
command=[sys.executable,str(R/'run_scenarios.py'),'--dataset','t4','--workers','4','--out',str(out)]
print('COMMAND',json.dumps(command),flush=True)
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
with (R/'reviews/t4_fresh_reproduction.log').open('x') as log:
 log.write(stamp()+' '+json.dumps(command)+'\n');log.flush()
 rc=subprocess.call(command,stdout=log,stderr=subprocess.STDOUT,cwd=R,env=env)
print(stamp(),'REPRODUCTION_EXIT',rc,flush=True)
if rc:raise SystemExit(rc)
command=[sys.executable,str(R/'compare_reproduction.py'),str(R/'results/t4'),str(out),'--out',str(R/'reviews/fresh_t4_comparison.json')]
rc=subprocess.call(command,cwd=R,env=env)
print(stamp(),'DEFAULT_COMPARISON_EXIT',rc,flush=True)
raise SystemExit(rc)
