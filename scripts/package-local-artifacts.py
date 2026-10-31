#!/usr/bin/env python3
"""Stream inventoried local artifacts into release-sized tar parts (no credentials)."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
entries = json.loads((root / 'provenance/local-artifacts.json').read_text())
for entry in entries:
    p = Path(entry['path'])
    if not p.is_file() or p.stat().st_size != entry['bytes']:
        raise SystemExit(f'Missing or changed artifact: {p}')
args.output.mkdir(parents=True, exist_ok=False)
with subprocess.Popen(['split', '-b', '1900M', '-d', '-a', '3', '-', str(args.output / 'evtgen-local.tar.gz.part-')], stdin=subprocess.PIPE) as split:
    with tarfile.open(fileobj=split.stdin, mode='w|gz', compresslevel=1) as tar:
        for entry in entries:
            p = Path(entry['path'])
            tar.add(p, arcname=p.relative_to('/home/darthlucas/programs/test'))
    split.stdin.close()
    if split.wait() != 0:
        raise SystemExit('split failed')
with (args.output / 'SHA256SUMS').open('w') as sums:
    for p in sorted(args.output.glob('*.part-*')):
        with p.open('rb') as stream:
            digest = hashlib.file_digest(stream, 'sha256').hexdigest()
        sums.write(f'{digest}  {p.name}\n')
print(args.output)
