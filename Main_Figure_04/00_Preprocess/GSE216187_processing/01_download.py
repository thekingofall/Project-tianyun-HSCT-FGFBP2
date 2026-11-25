from pathlib import Path
import json,time,subprocess,hashlib,gzip,csv,sys
root=Path(__file__).resolve().parents[1]
plan=json.loads((root/'resources/download_plan.json').read_text())
for attempt in range(1,7):
 pending=[x for x in plan if not (root/x['path']).is_file() or (root/x['path']).with_name(x['filename']+'.aria2').exists()]
 print('ROUND',attempt,'pending',len(pending),flush=True)
 if not pending:break
 p=root/'resources'/('aria2_round_%02d.txt'%attempt)
 with p.open('x') as f:
  for x in pending:f.write(x['url']+'\n  dir='+str(root/'data/raw')+'\n  out='+x['filename']+'\n')
 cmd=['aria2c','--continue=true','--allow-overwrite=false','--auto-file-renaming=false','--file-allocation=none','--max-concurrent-downloads=3','--split=1','--max-connection-per-server=1','--max-tries=8','--retry-wait=15','--connect-timeout=30','--timeout=60','--summary-interval=30','--input-file='+str(p)]
 with (root/'logs'/('aria2_round_%02d.log'%attempt)).open('x') as f:rc=subprocess.call(cmd,stdout=f,stderr=subprocess.STDOUT)
 if rc==0:break
 if rc:time.sleep(20)
else:raise RuntimeError('NCBI downloads incomplete; completed files preserved')
for x in plan:
 p=root/x['path'];assert p.exists() and not p.with_name(p.name+'.aria2').exists()
 sha=hashlib.sha256()
 with p.open('rb') as f:
  while True:
   b=f.read(1024*1024)
   if not b:break
   sha.update(b)
 with gzip.open(p,'rb') as f:
  while f.read(1024*1024):pass
 x['bytes']=p.stat().st_size;x['sha256']=sha.hexdigest()
 print('COMPLETE',x['filename'],x['bytes'],flush=True)
with (root/'data/raw/download_manifest.json').open('x') as f:json.dump({'status':'complete','files':plan},f,indent=2)
print('ALL_EIGHT_FILES_VALIDATED',flush=True)
