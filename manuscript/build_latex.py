#!/usr/bin/env python3
"""T1: main.qmd -> main.pdf + latex/main.tex (arXiv / JoE submission build).

House format (matches the reference paper.tex of 2026-09-16): 12pt letterpaper
article, 1in margins, \\setstretch{1.5} body, title \\thanks for acknowledgements
and funding, author \\thanks for affiliation and contact, Keywords line after the
abstract, colorlinks. Everything except the two \\thanks strings comes from the
qmd's `pdf` format block; the \\thanks are patched in here because they are
LaTeX-only and would leak into the HTML if put in the YAML.

Re-run after any edit to main.qmd. Requires quarto + pdflatex.
"""
import subprocess, shutil, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
LTX  = os.path.join(HERE, 'latex')

# --- the two footnotes on the title page -------------------------------------
# The title footnote gives the fixed reference to the public replication archive
# (repository + tag + the commit checked against this manuscript) and the funding.
# ARCHIVE_COMMIT is the seven-character SHA of the checked commit; the literal
# COMMIT7 is a deliberate placeholder that make_arxiv_zips.sh refuses to submit,
# so the paper can never go out naming a commit that does not exist.
ARCHIVE_REPO   = "https://github.com/sokubo/paper-panel-conditioning-replication"
ARCHIVE_TAG    = "paper-v1.0"
ARCHIVE_COMMIT = "f06283e"

TITLE_THANKS = (
    r"\thanks{Code reproducing every numerical result in this paper, including the "
    r"downstream identities of Section 6 and every row of Tables 1 and 2, is in the "
    r"replication archive at \url{%s} (fixed version: tag \texttt{%s}, commit "
    r"\texttt{%s}). No restricted data are used. This work was supported by JSPS "
    r"KAKENHI Grant Number 22K13525.}" % (ARCHIVE_REPO, ARCHIVE_TAG, ARCHIVE_COMMIT)
)
AUTHOR_THANKS = (
    r"\thanks{Department of Sociology, Toyo University, Tokyo, Japan. "
    r"Email: okubo080@toyo.jp. Website: sokubo.github.io.}"
)
DATE = 'September 23, 2026'
KEYWORDS = (r"\noindent\textbf{Keywords:} panel conditioning; identification; "
            r"age-period-cohort; two-way fixed effects; event study; refreshment samples")
# -----------------------------------------------------------------------------

# Refuse to build a submission PDF that names a commit which does not exist. A draft build for internal
# reading is allowed with ALLOW_PLACEHOLDER=1, and is stamped DRAFT on the title page so it cannot be
# mistaken for a release build (the build rejects placeholder identifiers).
DRAFT = False
if re.fullmatch(r'COMMIT[0-9]*', ARCHIVE_COMMIT) or not re.fullmatch(r'[0-9a-f]{7,40}', ARCHIVE_COMMIT):
    if os.environ.get('ALLOW_PLACEHOLDER') == '1':
        DRAFT = True
        print('WARNING: ARCHIVE_COMMIT is the placeholder %r; building a DRAFT (not for submission)' % ARCHIVE_COMMIT)
    else:
        sys.exit('ARCHIVE_COMMIT is the placeholder %r: publish the archive, set the checked commit here and in '
                 'main.qmd (Data and code availability), then rebuild. For an internal draft: ALLOW_PLACEHOLDER=1 python3 build_latex.py'
                 % ARCHIVE_COMMIT)
qmd = open(os.path.join(HERE, 'main.qmd')).read()
if not DRAFT and re.search(r'`COMMIT[0-9]*`', qmd):
    sys.exit('main.qmd still contains a COMMIT placeholder in the Data and code availability section')

subprocess.run(['quarto', 'render', 'main.qmd', '--to', 'pdf'], cwd=HERE, check=True)

tex_src = os.path.join(HERE, 'main.tex')          # keep-tex output
tex = open(tex_src).read()
if 'CSLReferences' not in tex:
    sys.exit('bibliography not resolved into the tex — check citeproc')

# global theorem numbering (Theorem 1, 2, ...) so that the PDF matches the HTML and the prose cross-references
tex = re.sub(r'(\\newtheorem\{(theorem|corollary|proposition|lemma|definition)\}\{[^}]*\})\[section\]', r'\1', tex)

# title \thanks: attach to the end of the \title{...} argument
m = re.search(r'\\title\{', tex)
if not m:
    sys.exit('no \\title in generated tex')
i, depth = m.end(), 1
while depth:
    if tex[i] == '{': depth += 1
    elif tex[i] == '}': depth -= 1
    i += 1
tex = tex[:i-1] + TITLE_THANKS + tex[i-1:]

# author \thanks
m = re.search(r'\\author\{([^}]*)\}', tex)
if not m:
    sys.exit('no \\author in generated tex')
tex = tex[:m.end()-1] + AUTHOR_THANKS + tex[m.end()-1:]

# visible Keywords line after the abstract (house format)
if r'\textbf{Keywords:}' not in tex:
    tex = tex.replace(r'\end{abstract}', '\\end{abstract}\n\n' + KEYWORDS + '\n', 1)

# date, spelled out as in the house format (a draft build is stamped as such)
tex = re.sub(r'\\date\{[^}]*\}', r'\\date{' + DATE + (r' --- DRAFT BUILD, archive commit not yet fixed' if DRAFT else '') + '}', tex, count=1)

os.makedirs(LTX, exist_ok=True)
open(os.path.join(LTX, 'main.tex'), 'w').write(tex)
os.remove(tex_src)
shutil.copy(os.path.join(HERE, 'references.bib'), os.path.join(LTX, 'references.bib'))
# citations are already resolved by citeproc, so a plain two-pass run suffices
for _ in range(2):
    r = subprocess.run(['xelatex', '-interaction=nonstopmode', 'main.tex'],
                       cwd=LTX, capture_output=True, text=True)
open(os.path.join(LTX, 'full.log'), 'w').write(r.stdout[-20000:])

log = open(os.path.join(LTX, 'main.log'), errors='ignore').read()
errs  = [l for l in log.splitlines() if l.startswith('!')]
undef = [l for l in log.splitlines() if 'Reference' in l and 'undefined' in l]
print('LaTeX errors:', len(errs), errs[:5])
print('undefined references:', len(undef), undef[:5])
shutil.copy(os.path.join(LTX, 'main.pdf'), os.path.join(HERE, 'main.pdf'))
print('pages:', subprocess.run(['pdfinfo', os.path.join(HERE, 'main.pdf')],
      capture_output=True, text=True).stdout.split('Pages:')[1].split()[0])
