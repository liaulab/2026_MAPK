"""Static site generator for the interactive figures.

Produces a plain HTML site suitable for GitHub Pages. No build step, no
framework, no server: the output is HTML, one CSS file, and the figure pages
themselves. Each figure is embedded with ``<iframe loading="lazy">``, so a
section page with two hundred figures costs nothing until the reader opens
that group, and the browser caches plotly.js once across every figure.

The only JavaScript is a few lines that filter the group list on a section
page. Turn it off and the page still works: every group is a plain anchor.
"""

from __future__ import annotations

import html
import json
from datetime import date
from pathlib import Path

# ---------------------------------------------------------------------------
# SECTION DEFINITIONS
# ---------------------------------------------------------------------------

SECTIONS = [
    {
        'slug': 'scatterplots-gof',
        'title': 'GOF screens',
        'blurb': 'Resistance screens: base-editing Z-score against amino-acid '
                 'position, one panel per gene, every inhibitor on the same axes. '
                 'Circles are ABE, diamonds CBE.',
        'group_label': 'Gene',
    },
    {
        'slug': 'scatterplots-lof',
        'title': 'LOF screens',
        'blurb': 'Dropout screens: base-editing Z-score against amino-acid '
                 'position, one panel per gene. Circles are ABE, diamonds CBE.',
        'group_label': 'Gene',
    },
    {
        'slug': 'scatterplots-meki',
        'title': 'MEKi screens',
        'blurb': 'MEK-inhibitor screens: base-editing Z-score against amino-acid '
                 'position, one panel per gene, every condition on the same axes. '
                 'Circles are ABE, diamonds CBE.',
        'group_label': 'Gene',
    },
    {
        'slug': 'pocket-pdb',
        'title': 'Pocket hits — co-crystal structures',
        'blurb': 'Resistance hits classified as inside or outside the drug-binding '
                 'pocket, using five inhibitor co-crystal structures.',
        'group_label': 'Structure',
    },
    {
        'slug': 'pocket-alphafold',
        'title': 'Pocket hits — AlphaFold models',
        'blurb': 'The same analysis extended to every paralog, using AlphaFold '
                 'models with the inhibitor aligned into the pocket.',
        'group_label': 'Structure',
    },
    {
        'slug': 'splice-heatmaps',
        'title': 'Splice-guide validation heatmaps',
        'blurb': 'Mean Z-score per gene and condition for the essential-splice-site '
                 'guides, i.e. the phenotype of ablating each node. Each panel is '
                 'drawn across the whole pathway and across the RAF paralogs alone.',
        'group_label': 'Gene set',
    },
    {
        'slug': 'boxplots',
        'title': 'Per-gene boxplots',
        'blurb': 'Z-score distribution of every guide targeting each gene, one panel '
                 'per screen, editor and condition, with the non-targeting and '
                 'essential-splice-site controls on the same axis. The last panel '
                 'brackets the two-sided Mann-Whitney U test of controls against '
                 'essential guides, within each editor.',
        'group_label': 'Screen',
    },
    {
        'slug': 'sankey',
        'title': 'Validated hits by inhibitor',
        'blurb': 'Where each inhibitor’s resistance hits fall in the pathway. '
                 'Grey links are on-target, pink downstream, blue upstream.',
        'group_label': 'Figure',
    },
    {
        'slug': 'trajectory',
        'title': 'GMM clustering and trajectory',
        'blurb': 'Guides clustered on their inhibitor-response profiles. Each '
                 'clustered heatmap is blocked by cluster and carries a gene bar '
                 'per guide; the validation screen additionally gets a '
                 'minimum-spanning-tree pseudotime through the mixture components.',
        'group_label': 'Screen',
    },
    {
        'slug': 'lof-vs-dms',
        'title': 'Hyperactivation versus deep mutational scanning',
        'blurb': 'Guides that hyperactivate the pathway, compared against published '
                 'DMS data for PTPN11, KRAS and MAPK1.',
        'group_label': 'Figure',
    },
    {
        'slug': 'cbioportal',
        'title': 'cBioPortal lollipops',
        'blurb': 'Pocket distance, validated screen hits and patient mutation counts '
                 'on one position axis.',
        'group_label': 'Gene',
    },
]

SECTION_BY_SLUG = {s['slug']: s for s in SECTIONS}

STYLESHEET = """\
/* MAPK interactive data pages. Plain CSS, no framework. */

:root {
  --ink: #1a1a1a;
  --ink-soft: #565656;
  --rule: #e2e2e2;
  --bg: #ffffff;
  --bg-soft: #f7f7f8;
  --accent: #7584E6;
  --accent-dark: #4c5bb5;
  --max: 1180px;
}

* { box-sizing: border-box; }

body {
  margin: 0;
  font: 16px/1.6 -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
  color: var(--ink);
  background: var(--bg);
}

header.site {
  border-bottom: 1px solid var(--rule);
  background: var(--bg-soft);
}
header.site .inner, main, footer.site .inner {
  max-width: var(--max);
  margin: 0 auto;
  padding: 0 24px;
}
header.site .inner { padding-top: 28px; padding-bottom: 28px; }

/* Title on the left, theme toggle on the right. */
.masthead { display: flex; align-items: flex-start; justify-content: space-between; gap: 16px; }
.masthead .lede { margin: 0; }
button.theme {
  font: inherit;
  font-size: 0.85rem;
  color: var(--ink-soft);
  background: var(--bg);
  border: 1px solid var(--rule);
  border-radius: 999px;
  padding: 5px 13px;
  cursor: pointer;
  flex: 0 0 auto;
}
button.theme:hover { color: var(--accent-dark); border-color: var(--accent-dark); }

h1 { font-size: 1.7rem; margin: 0 0 6px; letter-spacing: -0.01em; }
h2 { font-size: 1.2rem; margin: 40px 0 10px; letter-spacing: -0.01em; }
h3 { font-size: 1rem; margin: 0; font-weight: 600; }
p  { margin: 0 0 12px; }
.lede { color: var(--ink-soft); max-width: 68ch; }

a { color: var(--accent-dark); }
a:hover { color: var(--accent); }

nav.crumbs { font-size: 0.85rem; color: var(--ink-soft); margin: 18px 0 0; }
nav.crumbs a { text-decoration: none; }

main { padding-bottom: 64px; }

/* Section cards on the landing page, five across. */
.cards {
  display: grid;
  grid-template-columns: repeat(5, 1fr);
  gap: 14px;
  margin: 24px 0 8px;
  padding: 0;
  list-style: none;
}
.card {
  border: 1px solid var(--rule);
  border-radius: 8px;
  padding: 14px 15px;
  background: var(--bg);
  display: flex;
  flex-direction: column;
}
.card a { text-decoration: none; }
.card p { font-size: 0.85rem; color: var(--ink-soft); margin: 8px 0 0; }
.card .count { font-size: 0.8rem; color: var(--ink-soft); margin-top: auto; padding-top: 10px; display: block; }

/* Five across needs the room; step down rather than squeeze. */
@media (max-width: 1100px) { .cards { grid-template-columns: repeat(3, 1fr); } }
@media (max-width: 760px)  { .cards { grid-template-columns: repeat(2, 1fr); } }
@media (max-width: 520px)  { .cards { grid-template-columns: 1fr; } }

/* Filter box */
.filter {
  margin: 22px 0 8px;
  display: flex;
  gap: 10px;
  align-items: center;
  flex-wrap: wrap;
}
.filter input {
  font: inherit;
  padding: 7px 11px;
  border: 1px solid var(--rule);
  border-radius: 6px;
  min-width: 240px;
}
.filter .hint { font-size: 0.85rem; color: var(--ink-soft); }

/* Figure groups */
details.group {
  border: 1px solid var(--rule);
  border-radius: 8px;
  margin: 10px 0;
  background: var(--bg);
}
details.group > summary {
  cursor: pointer;
  padding: 12px 16px;
  font-weight: 600;
  list-style: none;
  display: flex;
  align-items: baseline;
  gap: 8px;
}
details.group > summary::-webkit-details-marker { display: none; }
details.group > summary::before {
  content: "\\25B8";
  color: var(--ink-soft);
  flex: 0 0 auto;
}
details.group[open] > summary::before { content: "\\25BE"; }
details.group > summary .label { flex: 1 1 auto; }
details.group > summary .n { font-weight: 400; color: var(--ink-soft); font-size: 0.85rem; flex: 0 0 auto; }
details.group > summary:hover { color: var(--accent-dark); }
details.group .body { padding: 4px 16px 18px; border-top: 1px solid var(--rule); }

.figure { margin: 18px 0 26px; }
.figure figcaption {
  font-size: 0.85rem;
  color: var(--ink-soft);
  margin-bottom: 6px;
  display: flex;
  justify-content: space-between;
  gap: 12px;
  flex-wrap: wrap;
}
.figure figcaption .name { color: var(--ink); font-weight: 600; }
.figure iframe {
  width: 100%;
  border: 1px solid var(--rule);
  border-radius: 6px;
  background: #fff;
  display: block;
}

table.summary { border-collapse: collapse; width: 100%; font-size: 0.9rem; margin: 12px 0 20px; }
table.summary th, table.summary td { text-align: left; padding: 7px 10px; border-bottom: 1px solid var(--rule); }
table.summary th { font-weight: 600; }
table.summary td.num { text-align: right; font-variant-numeric: tabular-nums; }

footer.site {
  border-top: 1px solid var(--rule);
  background: var(--bg-soft);
  font-size: 0.85rem;
  color: var(--ink-soft);
}
footer.site .inner { padding-top: 20px; padding-bottom: 28px; }

code { background: var(--bg-soft); padding: 1px 5px; border-radius: 4px; font-size: 0.9em; }

/* Dark palette. The toggle writes data-theme on <html>; with it unset the page
   follows the operating system, which is what it did before the toggle. Both
   directions have to be spelled out so an explicit choice beats the OS either
   way. Figure frames stay white in both themes: the plots inside them are
   drawn for a light background. */
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --ink: #ececec; --ink-soft: #a6a6a6; --rule: #333;
    --bg: #151517; --bg-soft: #1d1d20; --accent-dark: #9aa6ee;
  }
}
:root[data-theme="dark"] {
  --ink: #ececec; --ink-soft: #a6a6a6; --rule: #333;
  --bg: #151517; --bg-soft: #1d1d20; --accent-dark: #9aa6ee;
}
"""

# Runs in <head>, before anything paints, so a stored choice does not flash the
# other theme first.
THEME_BOOT = (
    "(function(){try{var t=localStorage.getItem('mapk-theme');"
    "if(t==='light'||t==='dark'){document.documentElement.dataset.theme=t;}}"
    "catch(e){}})();"
)

THEME_SCRIPT = """\
// Theme toggle, cycling auto -> light -> dark. The choice is remembered per
// browser. Without JS the button stays hidden and the page follows the OS.
(function () {
  var button = document.getElementById('theme');
  if (!button) return;
  var order = ['auto', 'light', 'dark'];

  function stored() {
    try { return localStorage.getItem('mapk-theme') || 'auto'; } catch (e) { return 'auto'; }
  }
  function apply(mode) {
    if (mode === 'auto') { delete document.documentElement.dataset.theme; }
    else { document.documentElement.dataset.theme = mode; }
    button.textContent = 'Theme: ' + mode;
    try { localStorage.setItem('mapk-theme', mode); } catch (e) {}
  }

  button.hidden = false;
  apply(stored());
  button.addEventListener('click', function () {
    apply(order[(order.indexOf(stored()) + 1) % order.length]);
  });
})();
"""

FILTER_SCRIPT = """\
// Filters the group list. Without JS every group is still present and open-able.
(function () {
  var box = document.getElementById('filter');
  if (!box) return;
  var groups = Array.prototype.slice.call(document.querySelectorAll('details.group'));
  box.addEventListener('input', function () {
    var q = box.value.trim().toLowerCase();
    groups.forEach(function (g) {
      var hit = !q || g.dataset.search.indexOf(q) !== -1;
      g.style.display = hit ? '' : 'none';
      if (q && hit) { g.open = true; }
    });
  });
})();
"""


# ---------------------------------------------------------------------------
# PAGE BUILDING
# ---------------------------------------------------------------------------


def _page(title: str, body: str, depth: int = 0, subtitle: str = '',
          script: str = '') -> str:
    """Wrap body content in the shared page chrome."""
    root = '../' * depth
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<link rel="stylesheet" href="{root}assets/style.css">
<script>{THEME_BOOT}</script>
</head>
<body>
<header class="site">
  <div class="inner">
    <div class="masthead">
      <div>
        <h1><a href="{root}index.html" style="text-decoration:none;color:inherit">MAPK base-editing screens</a></h1>
        <p class="lede">{subtitle}</p>
      </div>
      <button id="theme" class="theme" type="button" hidden
              aria-label="Switch between automatic, light and dark appearance">Theme</button>
    </div>
  </div>
</header>
<main>
{body}
</main>
<footer class="site">
  <div class="inner">
    <p>Interactive companion to the MAPK base-editing screen manuscript.
       Figures generated with <a href="https://github.com/liaulab/be-scan">be_scan</a>
       and Plotly. Built {date.today().isoformat()}.</p>
    <p>Every figure is also available as a print-ready SVG in the repository
       under <code>Outputs/</code>. Screen data is in
       <code>TableS1-ScreenData.xlsx</code>.</p>
  </div>
</footer>
<script>{THEME_SCRIPT}{script}</script>
</body>
</html>
"""


def _figure_block(entry: dict, depth: int, height: int) -> str:
    """One lazily-loaded figure with its caption.

    The frame height comes from the figure itself where known, so a short wide
    histogram and a tall lollipop each get the space they need.
    """
    root = '../' * depth
    height = int(entry.get('height') or height) + 16
    caption = (f'<span class="note">{html.escape(entry["caption"])}</span>'
               if entry.get('caption') else '')
    return f"""      <figure class="figure">
        <figcaption>
          <span class="name">{html.escape(entry['title'])}</span>
          {caption}
        </figcaption>
        <iframe src="{root}{entry['path']}" height="{height}" loading="lazy"
                title="{html.escape(entry['title'])}"></iframe>
      </figure>
"""


def _iframe_height(section_slug: str) -> int:
    """Fallback frame height, used only when a figure records no height."""
    return {
        'scatterplots-gof': 340,
        'scatterplots-lof': 340,
        'scatterplots-meki': 340,
        'pocket-pdb': 400,
        'pocket-alphafold': 400,
        'splice-heatmaps': 460,
        'sankey': 560,
        'trajectory': 620,
        'lof-vs-dms': 420,
        'cbioportal': 700,
    }.get(section_slug, 440)


def build_section_page(section: dict, entries: list[dict], depth: int = 0) -> str:
    """A section page: one collapsible group per gene or structure."""
    groups: dict[str, list[dict]] = {}
    for entry in entries:
        groups.setdefault(entry['group'], []).append(entry)

    height = _iframe_height(section['slug'])
    blocks = []
    for group, items in groups.items():
        search = ' '.join([group] + [i['title'] for i in items]).lower()
        figures = ''.join(_figure_block(i, depth, height) for i in items)
        # The first group opens on arrival so the page is not a wall of
        # closed rows; the rest stay shut so nothing else loads.
        is_first = not blocks
        blocks.append(f"""  <details class="group"{' open' if is_first else ''} data-search="{html.escape(search)}">
    <summary><span class="label">{html.escape(group)}</span><span class="n">{len(items)} figure{'s' if len(items) != 1 else ''}</span></summary>
    <div class="body">
{figures}    </div>
  </details>
""")

    body = f"""  <nav class="crumbs"><a href="{'../' * depth}index.html">Overview</a> &rsaquo; {html.escape(section['title'])}</nav>
  <h2>{html.escape(section['title'])}</h2>
  <p class="lede">{section['blurb']}</p>
  <div class="filter">
    <input id="filter" type="search" placeholder="Filter by {section['group_label'].lower()} or figure name"
           aria-label="Filter figures">
    <span class="hint">{len(entries)} figures in {len(groups)} group{'s' if len(groups) != 1 else ''} &middot; figures load when a group is opened</span>
  </div>
{''.join(blocks)}"""
    return _page(section['title'], body, depth=depth,
                 subtitle=section['title'], script=FILTER_SCRIPT)


def build_index(manifest: list[dict], summary_rows: list[dict] | None = None) -> str:
    """The landing page: what this is, section cards, dataset summary."""
    counts: dict[str, int] = {}
    for entry in manifest:
        counts[entry['section']] = counts.get(entry['section'], 0) + 1

    cards = []
    for section in SECTIONS:
        n = counts.get(section['slug'], 0)
        if not n:
            continue
        cards.append(f"""    <li class="card">
      <h3><a href="{section['slug']}.html">{html.escape(section['title'])}</a></h3>
      <p>{section['blurb']}</p>
      <span class="count">{n} figure{'s' if n != 1 else ''}</span>
    </li>
""")

    table = ''
    if summary_rows:
        rows = ''.join(
            f"      <tr><td>{html.escape(str(r['screen']))}</td>"
            f"<td>{html.escape(str(r['editor']))}</td>"
            f"<td class='num'>{r['guides']:,}</td>"
            f"<td class='num'>{r['genes']}</td>"
            f"<td class='num'>{r['conditions']}</td></tr>\n"
            for r in summary_rows)
        table = f"""  <h2>Screens</h2>
  <table class="summary">
    <thead><tr><th>Screen</th><th>Editor</th><th>Guides</th><th>Genes</th><th>Conditions</th></tr></thead>
    <tbody>
{rows}    </tbody>
  </table>
"""

    body = f"""  <h2>Explore the data</h2>
  <p class="lede">{len(manifest)} figures across {len(counts)} analyses.</p>
  <ul class="cards">
{''.join(cards)}  </ul>
{table}"""
    return _page('MAPK base-editing screens', body, depth=0,
                 subtitle='An interactive companion to the manuscript')


def build_site(manifest: list[dict] | str | Path, site_root: Path,
               summary_rows: list[dict] | None = None) -> dict:
    """Write the whole site. Returns a per-page count."""
    if isinstance(manifest, (str, Path)):
        manifest = json.loads(Path(manifest).read_text())

    site_root = Path(site_root)
    site_root.mkdir(parents=True, exist_ok=True)
    (site_root / 'assets').mkdir(exist_ok=True)
    (site_root / 'assets' / 'style.css').write_text(STYLESHEET)

    # GitHub Pages runs Jekyll by default, which skips files and directories
    # beginning with an underscore. This turns that off.
    (site_root / '.nojekyll').write_text('')

    written = {}
    by_section: dict[str, list[dict]] = {}
    for entry in manifest:
        by_section.setdefault(entry['section'], []).append(entry)

    for slug, entries in by_section.items():
        section = SECTION_BY_SLUG.get(
            slug, {'slug': slug, 'title': slug, 'blurb': '', 'group_label': 'Group'})
        page = build_section_page(section, entries)
        (site_root / f'{slug}.html').write_text(page)
        written[f'{slug}.html'] = len(entries)

    (site_root / 'index.html').write_text(build_index(manifest, summary_rows))
    written['index.html'] = len(manifest)
    return written
