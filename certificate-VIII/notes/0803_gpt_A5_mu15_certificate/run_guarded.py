"""Run only a local task script, record exit/RSS, never overwrite prior logs."""
import sys, subprocess, os, time, json, hashlib
from pathlib import Path
OUT=Path(__file__).resolve().parent
label,script,*args=sys.argv[1:]
assert Path(label).name==label
path=(OUT/script).resolve();assert path.parent==OUT and path.suffix=='.py'
log=OUT/(label+'.log');record=OUT/(label+'_guard.json')
assert not log.exists() and not record.exists()
limit=8_500_000_000;peak=0;stop=False;start=time.monotonic()
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
command=[sys.executable,str(path),*args]
with log.open('x') as f:
 p=subprocess.Popen(command,stdout=f,stderr=subprocess.STDOUT,env=env)
 while p.poll() is None:
  proc=subprocess.run(['ps','-o','rss=','-p',str(p.pid)],capture_output=True,text=True)
  if proc.stdout.strip():peak=max(peak,int(proc.stdout.strip())*1024)
  if peak>limit:
   stop=True;p.terminate()
   try:p.wait(timeout=10)
   except subprocess.TimeoutExpired:p.kill();p.wait()
   break
  time.sleep(2)
 code=p.wait()
r={'command':command,'exit_code':code,'memory_stop':stop,'sampled_peak_RSS_bytes':peak,'limit_bytes':limit,'seconds':time.monotonic()-start,'script_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
record.write_text(json.dumps(r,indent=2));print(json.dumps(r),flush=True)
sys.exit(1 if stop else code)
