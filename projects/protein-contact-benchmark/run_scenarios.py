"""Execute independent frozen scenarios concurrently without changing results."""
from pathlib import Path
import argparse,concurrent.futures,json,os,subprocess,sys
ROOT=Path(__file__).resolve().parent
def main():
 p=argparse.ArgumentParser();p.add_argument('--dataset',required=True,choices=['t4','villin']);p.add_argument('--workers',type=int,default=4);p.add_argument('--out',type=Path);a=p.parse_args()
 out=(a.out or ROOT/'results'/a.dataset).resolve();out.mkdir(parents=True,exist_ok=True)
 protocol=ROOT/'protocol/protocol_v1.1.json';scenarios=json.loads(protocol.read_text())['scenarios']
 env=dict(os.environ);env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
 def run(s):
  name=s['name'];command=[sys.executable,str(ROOT/'benchmark_extension.py'),'--dataset',a.dataset,'--protocol',str(protocol),'--scenario',name,'--phase','all','--out',str(out)]
  print('START',a.dataset,name,flush=True)
  with (out/(name+'.log')).open('x') as log:r=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=env)
  print('FINISH',a.dataset,name,'exit',r.returncode,flush=True)
  return name,r.returncode
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
  result=list(pool.map(run,scenarios))
 if any(code for _,code in result):raise SystemExit('Failed scenarios: '+str(result))
if __name__=='__main__':main()
