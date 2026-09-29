"""Configuração única e proveniência dos quatro registros deste piloto."""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import platform
import subprocess
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data' / 'raw'
RESULTS = ROOT / 'results' / 'audit_2026_09_26'
SESSIONS = [('Base', '15-56-59'), ('BB', '16-04-54'), ('Pós 1', '16-36-02'), ('Pós 2', '16-44-15')]


def session_path(stamp):
    return RAW / f'OpenBCI-RAW-2023-11-28_{stamp}.txt'


def provenance(name, parameters):
    RESULTS.mkdir(parents=True, exist_ok=True)
    versions = {}
    for package in ('numpy', 'scipy', 'pandas', 'scikit-learn', 'torch', 'matplotlib'):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            pass
    paths = list(RAW.glob('*.txt')) + list((ROOT / 'src').rglob('*.py'))
    manifest = {'generated_utc': datetime.now(timezone.utc).isoformat(),
        'python': platform.python_version(), 'versions': versions, 'parameters': parameters,
        'git_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}}
    (RESULTS / f'{name}_manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False)+'\n')
