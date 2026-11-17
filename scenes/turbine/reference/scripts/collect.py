"""Fetch or verify the pinned public SWiFT subset. Standard library only.

Run from any directory. Existing files must match the lock; mismatches fail rather
than silently replacing evidence. Raw data lives in the repository's ignored local/.
"""
import argparse
import hashlib
import json
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parents[2]


def collect(offline=False, data_root=ROOT / 'data/raw'):
    manifest = json.loads((ROOT / 'reference/data-sources.json').read_text())
    total = 0
    for source in manifest['files']:
        path = data_root / source['path']
        if not path.exists():
            if offline:
                raise FileNotFoundError(path)
            path.parent.mkdir(parents=True, exist_ok=True)
            with urllib.request.urlopen(source['url'], timeout=90) as response:
                data = response.read()
            if (len(data) != source['bytes'] or
                    hashlib.sha256(data).hexdigest() != source['sha256']):
                raise ValueError(f'Download differs from pinned source: {source["url"]}')
            partial = path.with_suffix(path.suffix + '.partial')
            partial.write_bytes(data)
            partial.replace(path)
        data = path.read_bytes()
        if len(data) != source['bytes'] or hashlib.sha256(data).hexdigest() != source['sha256']:
            raise ValueError(f'Changed local evidence: {path}')
        total += len(data)
    print(f'PASS: {len(manifest["files"])} pinned source files, {total:,} bytes.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--data-root', type=Path, default=ROOT / 'data/raw')
    args = parser.parse_args()
    collect(args.offline, args.data_root)
