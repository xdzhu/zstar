"""Create an isolated clean manuscript candidate without changing artwork."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import zipfile


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('article', type=Path)
parser.add_argument('output', type=Path)
args = parser.parse_args()
source = args.article/'zstar_CPC-full.tex'
text = source.read_text(encoding='utf-8')
figures = re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}', text)
assert len(figures) == 8 and len(set(figures)) == 8, figures
assert all(Path(f).name == f and f.startswith('Figure_') and f.endswith('.pdf') for f in figures)
args.output.mkdir(parents=True, exist_ok=False)
(args.output/'zstar_CPC.tex').write_text('\\def\\ZStarClean{1}\n'+text, encoding='utf-8')
for name in [*figures, 'zstar.bib', 'zstar-elsarticle-num.bst']:
    shutil.copy2(args.article/name, args.output/name)
shutil.copy2(args.article/'zstar_CPC-clean-20260904.bbl', args.output/'zstar_CPC.bbl')
(args.output/'BUILD.md').write_text(
    '# ZStar CPC candidate 0.3.0rc1\n\n'
    'Compile with `latexmk -pdf zstar_CPC.tex`. All figure paths are basenames.\n\n'
    '**Not yet submitted:** Figure 1 is retained unchanged at the author\'s request\n'
    'and still requires the final Unified-framework artwork. No other figure\n'
    'has been edited by this packaging step.\n', encoding='utf-8')
hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
          for p in sorted(args.output.iterdir()) if p.is_file()}
(args.output/'SHA256.json').write_text(json.dumps(hashes, indent=2)+'\n', encoding='utf-8')
archive_path = args.output.with_name(args.output.name + '.zip')
with zipfile.ZipFile(archive_path, 'w', zipfile.ZIP_DEFLATED) as archive:
    for path in sorted(args.output.iterdir()):
        archive.write(path, path.name)
print(archive_path)
