"""0792: run one generate_lowmem.py child, sample its RSS, stop it above GUARD_BYTES (default 8.5e9)."""
from pathlib import Path
import subprocess,sys,time,json,os
out=Path(__file__).parent/'runs';out.mkdir(exist_ok=True)
args=sys.argv[1:];tag=f"{args[5]}_p{args[4]}"
log=out/(tag+'.log');record=out/(tag+'_guard.json')
assert not log.exists() and not record.exists()
limit=int(os.environ.get('GUARD_BYTES','8500000000'))
cmd=[sys.executable,str(Path(__file__).parent/'generate_lowmem.py'),*args]
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
start=time.monotonic();peak=0;stopped=False
with log.open('w') as f:
 p=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT,env=env)
 while p.poll() is None:
  z=subprocess.run(['ps','-o','rss=','-p',str(p.pid)],capture_output=True,text=True)
  if z.stdout.strip():peak=max(peak,int(z.stdout.strip())*1024)
  if peak>limit:
   stopped=True;p.terminate()
   try:p.wait(timeout=10)
   except subprocess.TimeoutExpired:p.kill();p.wait()
   break
  time.sleep(2)
 code=p.wait()
record.write_text(json.dumps({'command':cmd[1:],'exit_code':code,'seconds':time.monotonic()-start,'sampled_peak_RSS_bytes':peak,'memory_stop':stopped,'threshold_bytes':limit},indent=2))
print(record.read_text());sys.exit(1 if stopped else code)
