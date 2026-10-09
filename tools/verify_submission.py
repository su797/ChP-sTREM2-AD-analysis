"""Verify the manuscript code snapshot, optionally against Additional file 3."""

import argparse
import hashlib
import json
from pathlib import Path
import sys
import zipfile


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--zip', type=Path, help='Original journal submission ZIP')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / 'submission_manifest.json').read_text(encoding='utf-8'))
    failures = []
    for entry in manifest['files']:
        path = root / manifest['repository_directory'] / entry['path']
        if not path.is_file():
            failures.append(f"Missing repository file: {entry['path']}")
        elif digest(path.read_bytes()) != entry['sha256']:
            failures.append(f"Repository hash mismatch: {entry['path']}")
    snapshot = root / manifest['repository_directory']
    actual = {p.relative_to(snapshot).as_posix() for p in snapshot.rglob('*') if p.is_file()}
    expected = {entry['path'] for entry in manifest['files']} | set(manifest['repository_only_files'])
    for path in sorted(actual - expected):
        failures.append(f'Unexpected snapshot file: {path}')
    if args.zip:
        if digest(args.zip.read_bytes()) != manifest['archive_sha256']:
            failures.append('Archive SHA-256 mismatch')
        with zipfile.ZipFile(args.zip) as archive:
            names = [info.filename for info in archive.infolist() if not info.is_dir()]
            expected_names = {manifest['archive_directory'] + '/' + entry['path'] for entry in manifest['files']}
            if len(names) != len(set(names)) or set(names) != expected_names:
                failures.append('Archive file list mismatch')
            for entry in manifest['files']:
                name = manifest['archive_directory'] + '/' + entry['path']
                if name not in names or digest(archive.read(name)) != entry['sha256']:
                    failures.append(f"Archive hash mismatch: {entry['path']}")
    if failures:
        print('\n'.join(failures), file=sys.stderr)
        return 1
    print(f"Verified {len(manifest['files'])} manuscript files" + (' and the original submission ZIP.' if args.zip else '.'))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
