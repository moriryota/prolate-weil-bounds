"""Run one owned child, stop before 10 GB. No other process is modified."""
from pathlib import Path
import subprocess,sys,time,json,os,signal
out=Path(__file__).parent;N,bits,q,tag=sys.argv[1:];log=out/(tag+'.log');record=out/(tag+'_guard.json')
assert not log.exists() and not record.exists()
cmd=[sys.executable,str(out/'generate.py'),N,bits,q,tag]
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
start=time.monotonic();peak=0;stopped=False
with log.open('w') as f:
 p=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT,env=env)
 while p.poll() is None:
  z=subprocess.run(['ps','-o','rss=','-p',str(p.pid)],capture_output=True,text=True)
  if z.stdout.strip():peak=max(peak,int(z.stdout.strip())*1024)
  if peak>8_500_000_000:
   stopped=True;p.terminate()
   try:p.wait(timeout=10)
   except subprocess.TimeoutExpired:p.kill();p.wait()
   break
  time.sleep(2)
 code=p.wait()
record.write_text(json.dumps({'command':cmd,'exit_code':code,'seconds':time.monotonic()-start,'sampled_peak_RSS_bytes':peak,'memory_stop':stopped,'threshold_bytes':8500000000},indent=2))
print(record.read_text());sys.exit(1 if stopped else code)
