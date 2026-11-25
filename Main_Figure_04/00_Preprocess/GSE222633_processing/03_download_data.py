from pathlib import Path
import urllib.request,urllib.error,json,re,csv,hashlib,gzip,time,concurrent.futures,threading

root=Path(__file__).resolve().parents[1]
series=(root/'resources/GSE222633_series.soft.txt').read_text(encoding='utf-8')
urls=json.loads((root/'resources/geo_supplementary_urls.json').read_text(encoding='utf-8'))
lock=threading.Lock()

def event(z):
    z['utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
    with lock:
        with (root/'logs/download_events.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(z)+'\n')
        print(json.dumps(z),flush=True)

def get_text(gsm):
    path=root/'resources'/('GSM_'+gsm+'.soft.txt')
    url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+gsm+'&targ=self&form=text&view=full'
    with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 GEO-analysis'}),timeout=60) as r:b=r.read()
    assert b.startswith(b'^SAMPLE'),gsm
    with path.open('xb') as f:f.write(b)
    text=b.decode('utf-8')
    title=next(x.split(' = ',1)[1].strip() for x in text.splitlines() if x.startswith('!Sample_title = '))
    return {'GSM':gsm,'title':title,'metadata_path':str(path.relative_to(root))}

def download(source):
    name=source.rsplit('/',1)[-1];url=source.replace('ftp://','https://',1)
    for attempt in range(2):
        folder=root/'data/raw' if attempt==0 else root/'data/raw/retry_1'
        folder.mkdir(parents=True,exist_ok=True);path=folder/name
        try:
            assert not path.exists(),str(path)
            req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 GEO-analysis'})
            sha=hashlib.sha256();md5=hashlib.md5();size=0;start=time.time()
            with urllib.request.urlopen(req,timeout=120) as r:
                expected=int(r.headers['Content-Length']) if r.headers.get('Content-Length') else None
                assert r.status==200,(r.status,r.headers.get('Content-Range'))
                with path.open('xb') as f:
                    while True:
                        b=r.read(1024*1024)
                        if not b:break
                        f.write(b);sha.update(b);md5.update(b);size+=len(b)
            if expected is not None:assert size==expected,(name,size,expected)
            with path.open('rb') as f:magic=f.read(8)
            if name.endswith('.h5'):assert magic==b'\x89HDF\r\n\x1a\n',name
            elif name.endswith('.gz'):
                assert magic[:2]==b'\x1f\x8b',name
                with gzip.open(path,'rb') as f:
                    while f.read(1024*1024):pass
            item={'filename':name,'path':str(path.relative_to(root)),'url':url,'bytes':size,'http_expected_bytes':expected,'sha256':sha.hexdigest(),'md5':md5.hexdigest(),'seconds':round(time.time()-start,1),'status':'complete'}
            event(item.copy());return item
        except Exception as e:
            event({'filename':name,'attempt':attempt+1,'status':'failed_attempt','error':repr(e)})
            if attempt==1:raise

event({'status':'started','files':len(urls)})
gsms=re.findall(r'^!Series_sample_id = (GSM\d+)',series,re.M)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:metadata=list(pool.map(get_text,gsms))
with (root/'resources/GSM_sample_titles.tsv').open('x',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(metadata[0]),delimiter='\t');w.writeheader();w.writerows(metadata)
event({'status':'sample_metadata_complete','GSM_samples':len(metadata)})
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:files=list(pool.map(download,urls))
with (root/'data/raw/download_manifest.json').open('x',encoding='utf-8') as f:json.dump({'status':'complete','files':files,'bytes':sum(x['bytes'] for x in files)},f,indent=2)
samples=[]
for item in files:
    match=re.fullmatch(r'GSE222633_([DR])(\d+)(pre|post|d90|d180)PBMC_filtered_counts\.h5',item['filename'])
    if not match:continue
    role,pair,tp=match.groups();sample=role+pair+tp+'PBMC'
    tcr=next(x for x in files if x['filename']=='GSE222633_'+sample+'_filtered_contigs.csv.gz')
    samples.append({'sample':sample,'pair_id':pair,'role':'Donor' if role=='D' else 'Recipient','timepoint':{'pre':'Pre_GCSF','post':'Post_GCSF','d90':'D90','d180':'D180'}[tp],'RNA_file':item['path'],'TCR_file':tcr['path']})
assert len(samples)==17
with (root/'data/processed/sample_manifest.tsv').open('x',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(samples[0]),delimiter='\t');w.writeheader();w.writerows(samples)
event({'status':'complete','RNA_TCR_samples':len(samples),'files':len(files),'bytes':sum(x['bytes'] for x in files)})
