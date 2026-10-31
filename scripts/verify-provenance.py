#!/usr/bin/env python3
import hashlib
import json
from pathlib import Path
root = Path(__file__).resolve().parents[1]
entries = json.loads((root / 'provenance/copied-files.json').read_text())
for entry in entries:
    path = root / entry['destination']
    with path.open('rb') as stream:
        actual = hashlib.file_digest(stream, 'sha256').hexdigest()
    if actual != entry['sha256']:
        raise SystemExit(f'Checksum mismatch: {path}')
print(f'{len(entries)} historical files verified.')
