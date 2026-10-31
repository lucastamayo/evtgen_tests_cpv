#!/usr/bin/env python3
"""Finish the authorized historical-data publication and record remote checksums."""
import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GH = str(Path.home() / '.local/bin/gh')
REPO = 'lucastamayo/evtgen_tests_cpv'
TAG = 'historical-data-v1'
STATUS = ROOT / 'build/publication-status.json'


def run(*args):
    return subprocess.check_output(args, cwd=ROOT, text=True)


def api(path):
    return json.loads(run(GH, 'api', path))


def save_status(state, **details):
    STATUS.write_text(json.dumps({'state': state, **details}, indent=2) + '\n')
    print(state, details, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wait-for-pid', type=int)
    args = parser.parse_args()
    save_status('uploading')
    if args.wait_for_pid:
        proc = Path('/proc') / str(args.wait_for_pid)
        while proc.exists():
            try:
                if proc.joinpath('exe').resolve().name != 'gh':
                    break
            except OSError:
                break
            time.sleep(15)
    manifest_path = ROOT / 'provenance/release-assets.json'
    manifest = json.loads(manifest_path.read_text())
    folder = ROOT / 'release-assets/historical-data'
    expected = {part['name']: {'bytes': part['bytes']} for part in manifest['parts']}
    for line in manifest['checksums'].splitlines():
        digest, name = line.split()
        expected[name]['digest'] = 'sha256:' + digest
    checksum_file = folder / 'SHA256SUMS'
    expected['SHA256SUMS'] = {'bytes': checksum_file.stat().st_size,
                             'digest': 'sha256:' + hashlib.sha256(checksum_file.read_bytes()).hexdigest()}
    release = next(r for r in api(f'repos/{REPO}/releases') if r['tag_name'] == TAG)
    release_id = release['id']
    def valid(asset, target):
        return asset is not None and asset['state'] == 'uploaded' and asset['size'] == target['bytes'] and asset.get('digest') == target['digest']
    for attempt in range(3):
        assets = {a['name']: a for a in api(f'repos/{REPO}/releases/{release_id}/assets')}
        missing = [name for name, target in expected.items() if not valid(assets.get(name), target)]
        if not missing:
            break
        if not release['draft']:
            raise RuntimeError('Existing public release differs from expected assets; refusing replacement.')
        for name in missing:
            save_status('uploading', file=name, attempt=attempt + 1)
            result = subprocess.run([GH, 'release', 'upload', TAG, str(folder / name), '--repo', REPO, '--clobber'], cwd=ROOT)
            if result.returncode:
                time.sleep(10)
    assets = {a['name']: a for a in api(f'repos/{REPO}/releases/{release_id}/assets')}
    if any(not valid(assets.get(name), target) for name, target in expected.items()):
        raise RuntimeError('Release assets are incomplete or differ from local checksums.')
    save_status('verified', asset_count=len(expected))
    if release['draft']:
        run(GH, 'release', 'edit', TAG, '--repo', REPO, '--draft=false')
    published = api(f'repos/{REPO}/releases/{release_id}')
    if published['draft']:
        raise RuntimeError('Release is still a draft.')
    url = published['html_url']
    manifest['publication'] = 'Published; all five remote asset sizes and SHA256 digests verified.'
    manifest['release_url'] = url
    manifest['published_at'] = published['published_at']
    manifest['remote_assets'] = [{'name': assets[name]['name'], 'size': assets[name]['size'], 'digest': assets[name]['digest'], 'url': assets[name]['browser_download_url']} for name in expected]
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    doc = ROOT / 'docs/PUBLICAR.md'
    text = doc.read_text().replace('Esses arquivos estão locais; ainda não foram\nenviados ao GitHub. O comando abaixo permite recriar o pacote em outro diretório.', f'Os arquivos foram publicados na [release histórica]({url}), com tamanhos e\nSHA256 remotos conferidos. O comando abaixo permite recriar o pacote em outro diretório.')
    doc.write_text(text)
    readme = ROOT / 'README.md'
    text = readme.read_text()
    text += f'\n## Download dos dados históricos\n\nOs quatro arquivos de dados e `SHA256SUMS` estão na [release histórica]({url}).\nConsulte `docs/PUBLICAR.md` para conferir os checksums e recuperar os eventos.\n'
    readme.write_text(text)
    subprocess.run(['git', 'add', 'scripts/finish-publication.py', 'provenance/release-assets.json', 'docs/PUBLICAR.md', 'README.md'], cwd=ROOT, check=True)
    subprocess.run(['git', 'commit', '-q', '-m', 'Publish historical release data with verified remote checksums'], cwd=ROOT, check=True)
    subprocess.run(['git', 'push', 'origin', 'main'], cwd=ROOT, check=True)
    save_status('complete', repository=f'https://github.com/{REPO}', release=url, asset_count=len(expected))


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        save_status('failed', error=str(exc))
        raise
