"""Download público, com hashes e licença, das bases de validação externa."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import time
import requests

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT/'data'/'external'
OUT = ROOT/'results'/'validation_2026_09_28'
PHYSIO = 'https://physionet.org/files/eegmat/1.0.0/'
MIRROR = 'https://physionet-open.s3.amazonaws.com/eegmat/1.0.0/'
OSF = {
    'readme.csv': 'https://osf.io/download/utvsj/',
    'eeg_G6.csv': 'https://osf.io/download/afc64/',
    'eeg_G40.csv': 'https://osf.io/download/x3t5u/',
    'relax_G6.csv': 'https://osf.io/download/hqbzv/',
    'relax_G40.csv': 'https://osf.io/download/4t3xv/',
}


def fetch(task):
    dataset, filename, url, expected = task
    path = DATA/dataset/filename
    path.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(3):
        try:
            if not path.exists():
                response = requests.get(url, timeout=(20,120))
                response.raise_for_status()
                if 'text/html' in response.headers.get('Content-Type',''):
                    raise ValueError(f'HTML inesperado: {url}')
                content = response.content
                digest = hashlib.sha256(content).hexdigest()
                if expected and digest != expected:
                    raise ValueError(f'Hash incorreto: {filename}')
                path.write_bytes(content)
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            if expected and digest != expected:
                raise ValueError(f'Arquivo existente não confere com SHA256SUMS: {filename}')
            return dict(dataset=dataset,file=str(path.relative_to(ROOT)),url=url,
                        bytes=path.stat().st_size,sha256=digest,upstream_hash_verified=bool(expected))
        except (requests.RequestException,ValueError):
            if attempt==2:
                raise
            time.sleep(attempt+1)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    checksum = requests.get(PHYSIO+'SHA256SUMS.txt',timeout=30)
    checksum.raise_for_status()
    hashes = {line.split()[1].lstrip('*'):line.split()[0] for line in checksum.text.splitlines() if line.strip()}
    tasks = [('eegmat', f'Subject{s:02d}_{condition}.edf',MIRROR+f'Subject{s:02d}_{condition}.edf',hashes[f'Subject{s:02d}_{condition}.edf'])
             for s in range(36) for condition in [1,2]]
    tasks += [('eegmat',name,PHYSIO+name,hashes.get(name)) for name in ['README.txt','RECORDS','subject-info.csv','SHA256SUMS.txt']]
    tasks += [('sudre2024',name,url,None) for name,url in OSF.items()]
    with ThreadPoolExecutor(max_workers=4) as pool:
        records = []
        for item in pool.map(fetch,tasks):
            records.append(item)
            print(f'{len(records)}/{len(tasks)} {item["file"]}',flush=True)
    manifest = dict(downloaded_utc=datetime.now(timezone.utc).isoformat(),
                    sources={'eegmat':'https://physionet.org/content/eegmat/1.0.0/',
                             'sudre2024':'https://osf.io/stzfd/'},files=records)
    (OUT/'download_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':
    main()
