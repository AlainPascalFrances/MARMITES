# -*- coding: utf-8 -*-
"""Every configuration key is READ by the model, or says why not.

The cookbook's Appendix B is the table: panel -> TOML key -> the code that
reads it. This keeps it true. A panel field the run ignores is worse than no
field -- the modeller believes they have set something -- and the audit that
produced Appendix B found 29 of them, including the whole of the State
variables panel and the Soil panel's soil parameters.

Two failures, both deliberate:

  * a NEW key nothing reads, and not listed below with a reason -- the next
    decorative field;
  * a key listed as NOT WIRED that the model now DOES read -- take it off
    the list and update Appendix B, or the table is lying the other way.
"""

import os
import re
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'app')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import marmites_config as mcfg                                # noqa: E402
from lib import schema                                        # noqa: E402

# Keys the run does NOT read, each with the reason. Mirrors Appendix B.
NOT_WIRED = {
    'obs.table': 'no reader',
    'obs.name_column': 'no reader',
    'obs.heads_prefix': 'obs_series uses its own default prefix',
    'obs.sm_prefix': 'no reader',
    'obs.aet_prefix': 'no reader',
    'obs.ro_prefix': 'no reader',
    'surface.meteo_zones': 'no reader in marmites_surface',
    'et.extwc_source': 'the build always takes extwc from thtr',
    'sfr.source': 'the network comes from the dataset CSV',
    'sfr.min_slope': 'no reader',
    'sfr.monotonic_bed': 'no reader',
    'lak.source': 'the run reads lak.geometry instead',
    'lak.polygons': 'no reader',
    'lak.surfdep': 'no reader',
    'lak.maxiter': 'no reader',
    'lak.stagechg': 'no reader',
    'postproc.hydro_year_start': 'no reader',
    'postproc.wb_unit': 'no reader',
    'postproc.obs_series': 'no reader',
    'postproc.sankey': 'no reader',
    'postproc.input_maps': 'no reader',
    'postproc.result_maps': 'no reader',
    'postproc.tick_trimester_years': 'no reader',
    'postproc.tick_semester_years': 'no reader',
    'paths.nwt_reference': 'defined only',
}
# Not model inputs at all, or asked ahead of the work package that uses them.
NOT_AN_INPUT = {
    'meta.config_version': 'file-format version, checked by validate()',
    'meta.description': 'provenance, copied into resolved_config.toml',
    'ui.execution': 'front-end', 'ui.runs_dir': 'front-end',
    'ui.port': 'front-end', 'ui.address': 'front-end',
    'ui.poll_secs': 'front-end',
}
AHEAD = ('crr.', 'pest.')                 # WP5, WP7


def _model_sources():
    out = {}
    for d, _dirs, files in os.walk(CODE):
        if (os.sep + 'app' in d or os.sep + 'tests' in d
                or os.sep + 'legacy' in d or '__pycache__' in d):
            continue
        for f in files:
            if f.endswith('.py') and f != 'marmites_config.py':
                p = os.path.join(d, f)
                out[p] = open(p, encoding='utf-8', errors='replace').read()
    p = os.path.join(CODE, 'tests', 'run_lamata_mf6.py')
    out[p] = open(p, encoding='utf-8', errors='replace').read()
    return out


SRC = _model_sources()
PROPS = SRC.get(os.path.join(CODE, 'ppMF6', 'marmites_props.py'), '')


def _is_read(dotted):
    """The same rules Appendix B used, and nothing looser."""
    parts = dotted.split('.')
    leaf, parent, sec_path = parts[-1], parts[-2], '.'.join(parts[:-1])
    for s in SRC.values():
        if re.search(r'\b%s\.%s\b' % (re.escape(parent), re.escape(leaf)), s):
            return True
        for a in set(re.findall(r'\b(\w+)\s*=\s*(?:self\.)?cfg\.%s\b(?!\.)'
                                % re.escape(sec_path), s)):
            if re.search(r'\b%s\.%s\b' % (a, re.escape(leaf)), s):
                return True
        if re.search(r"getattr\([^,()]*\b%s\s*,\s*'%s'"
                     % (re.escape(parent), re.escape(leaf)), s):
            return True
        if re.search(r'asdict\(\s*[\w.]*\b%s\s*\)' % re.escape(parent), s):
            return True
    # marmites_props reads the layer, UZF and boundary keys in loops over a
    # table of NAMES -- getattr(cfg.layers, name), pkg.<leaf> -- so the key
    # appears there as a quoted name or as pkg.<leaf>.
    if parts[0] in ('layers', 'uzf', 'ghb', 'drn'):
        if re.search(r"""['"]%s['"]""" % re.escape(leaf), PROPS) or \
                re.search(r'\bpkg\.%s\b' % re.escape(leaf), PROPS):
            return True
    if dotted == 'meta.model':
        return any('meta.model_name(' in s for s in SRC.values())
    return False


def _all_keys():
    cfg = mcfg.RunConfig.from_dict({})
    keys = []
    for sec in mcfg._SECTIONS:
        keys += [d for d, _v in schema.fields_of(cfg, sec)]
        keys += ['%s.%s' % (sec, n)
                 for n in getattr(type(getattr(cfg, sec)), '_ELEMENTS', {})]
    return keys


KEYS = _all_keys()


def test_every_key_is_read_or_says_why_not():
    silent = [k for k in KEYS
              if not _is_read(k) and k not in NOT_WIRED
              and k not in NOT_AN_INPUT and not k.startswith(AHEAD)]
    assert not silent, (
        'read by nothing, and not listed with a reason -- a panel field the '
        'run ignores: %s. Wire it, or add it to NOT_WIRED and to Appendix B '
        'of the cookbook.' % ', '.join(silent))


@pytest.mark.parametrize('key', sorted(NOT_WIRED))
def test_a_key_listed_as_not_wired_is_still_not_wired(key):
    """If the model now reads it, the list and Appendix B are wrong."""
    assert key in KEYS, '%s is no longer a configuration key' % key
    assert not _is_read(key), (
        '%s is READ now -- take it off NOT_WIRED and update Appendix B of '
        'the cookbook' % key)


def test_the_lists_do_not_overlap():
    assert not set(NOT_WIRED) & set(NOT_AN_INPUT)
