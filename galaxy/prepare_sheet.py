import os
import re
import sys

import pandas as pd

sheet_path, sheet_ext, sheet_index, id_col, pairs_path, out_path = sys.argv[1:7]

if sheet_ext == 'xlsx':
    df = pd.read_excel(sheet_path, sheet_name=int(sheet_index), dtype=str)
elif sheet_ext == 'csv':
    df = pd.read_csv(sheet_path, dtype=str)
else:
    df = pd.read_csv(sheet_path, sep='\t', dtype=str)

df.columns = [str(c).strip() for c in df.columns]
if id_col not in df.columns:
    sys.exit('ERROR: ID column "%s" not found in the tracking sheet. Columns present: %s'
             % (id_col, ', '.join(df.columns)))

# Paths in a user's sheet point at their own filesystem, never at Galaxy's.
stale = [c for c in df.columns if c.lower() in ('forward', 'fwd', 'reverse', 'rev')]
if stale:
    df = df.drop(columns=stale)
    sys.stderr.write('Ignoring file path columns from the tracking sheet: %s\n' % ', '.join(stale))

df[id_col] = df[id_col].astype(str).str.strip()
df = df[df[id_col].ne('') & df[id_col].ne('nan')]

pairs = pd.read_csv(pairs_path, sep='\t', header=None, dtype=str,
                    names=['sample', 'forward', 'reverse'], keep_default_na=False)
pairs['sample'] = pairs['sample'].str.strip()

bad = [s for s in pairs['sample'] if not re.fullmatch(r'[A-Za-z0-9._-]+', s)]
if bad:
    sys.exit('ERROR: these collection element identifiers contain characters other than '
             'letters, digits, dot, underscore or hyphen: %s' % ', '.join(bad))

dupes = sorted(set(df[id_col][df[id_col].duplicated()]))
if dupes:
    sys.exit('ERROR: duplicate sample IDs in the tracking sheet: %s' % ', '.join(dupes))

missing = sorted(set(pairs['sample']) - set(df[id_col]))
if missing:
    sys.exit('ERROR: these read pairs have no row in the tracking sheet (the collection element '
             'identifier must equal the value in the "%s" column): %s' % (id_col, ', '.join(missing)))

unused = sorted(set(df[id_col]) - set(pairs['sample']))
if unused:
    sys.stderr.write('Skipping %d tracking sheet rows with no matching reads: %s\n'
                     % (len(unused), ', '.join(unused)))
df = df[df[id_col].isin(set(pairs['sample']))].copy()

# Later pipeline steps match files to samples by substring, so an ID that
# occurs inside another ID would pick up the wrong sample's files.
ids = sorted(df[id_col])
clashes = ['%s / %s' % (a, b) for a in ids for b in ids if a != b and a in b]
if clashes:
    sys.exit('ERROR: some sample IDs occur inside other IDs, which the pipeline cannot keep apart: '
             '%s. Rename the samples, for example by zero-padding (S01 rather than S1).'
             % '; '.join(clashes))

paths = pairs.set_index('sample')
df['forward'] = [os.path.abspath(paths.at[s, 'forward']) for s in df[id_col]]
df['reverse'] = [os.path.abspath(paths.at[s, 'reverse']) for s in df[id_col]]
df.to_csv(out_path, index=False)
sys.stderr.write('Tracking sheet prepared for %d samples\n' % len(df))
