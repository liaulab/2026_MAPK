"""Analysis code for the MAPK base-editing screen manuscript.

Plotting primitives come from ``be_scan`` (https://github.com/liaulab/be-scan).
This package holds only what is specific to this manuscript.

``analysis`` computes and ``charts`` / ``lollipop`` / ``trajectory`` draw for
print; ``interactive`` draws the same figures for the web and ``website``
assembles them into the published pages. Both notebooks import ``analysis``,
so the print and interactive figures are guaranteed to show the same numbers.
"""

from . import analysis
from . import boxplot
from . import build_inputs
from . import charts
from . import config
from . import data
from . import dms
from . import interactive
from . import lollipop
from . import sankey
from . import structures
from . import trajectory
from . import website

__all__ = [
    'analysis', 'boxplot', 'build_inputs', 'charts', 'config', 'data', 'dms',
    'interactive',
    'lollipop', 'sankey', 'structures', 'trajectory', 'website',
]
