"""Monitor only this job's process tree, run no VM and stop at 8.5GB."""
import subprocess,time,json,sys,os
from pathlib import Path
O=Path(__file__).parent;tag=sys.argv[1];cmd=[sys.executable,*sys.argv[2:]];log=O/(tag+'.log');record=O/(tag+'_guard.json');assert not log.exists() and not record.exists()
start=time.monotonic();peak=0;stopped=False
with log.open('w') as f:
 p=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1'),start_new_session=True)
 while p.poll() is None:
  z=subprocess.run(['ps','-axo','pid=,ppid=,rss='],capture_output=True,text=True);rows=[list(map(int,line.split())) for line in z.stdout.splitlines() if len(line.split())==3];owned={p.pid}
  for _ in range(10):
   new=owned|{pid for pid,ppid,rss in rows if ppid in owned}
   if new==owned:break
   owned=new
  peak=max(peak,sum(rss*1024 for pid,ppid,rss in rows if pid in owned))
  if peak>8_500_000_000:
   stopped=True;os.killpg(p.pid,15)
   try:p.wait(timeout=10)
   except subprocess.TimeoutExpired:os.killpg(p.pid,9);p.wait()
   break
  time.sleep(2)
 code=p.wait()
r={'command':cmd,'exit_code':code,'seconds':time.monotonic()-start,'sampled_tree_peak_RSS_bytes':peak,'memory_stop':stopped};record.write_text(json.dumps(r,indent=2));print(json.dumps(r));sys.exit(1 if stopped else code)
