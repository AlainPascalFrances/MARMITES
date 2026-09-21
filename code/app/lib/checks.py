# -*- coding: utf-8 -*-
"""WP1d -- everything that can be said about a configuration before it runs.

DELIBERATELY STREAMLIT-FREE, like editor.py and runs.py, so the rules can be
tested in the model environment and so that ONE answer serves both places
that need it: the validation panel, which lists them, and the Launch button,
which refuses on them. Two implementations of "is this configuration ready"
would drift, and the one that drifts is always the one that guards the run.

THE DISTINCTION THAT MATTERS:

``error``    the run cannot proceed. A configuration the schema refuses, a
             MODFLOW library that is not there, a switch the panels and the
             file disagree about. Launch is blocked.
``warning``  the run will proceed and may not mean what was intended -- a
             truncated forcing series, saved state from another grid, a
             stress-period cap left on from a trial. Launch diverts to the
             validation panel so it is SEEN, and can then be taken anyway.
``info``     worth knowing, nothing to fix.

A check that cannot answer says nothing. Silence here means "not known", so
no check may report a problem it has not actually established -- a panel that
cries wolf is a panel that gets skipped, and this one guards the run.
"""

import os

__all__ = ['Check', 'ERROR', 'WARNING', 'INFO', 'collect', 'collect_default',
           'worst', 'count_by_level', 'CheckError']

ERROR = 'error'
WARNING = 'warning'
INFO = 'info'

_ORDER = {ERROR: 0, WARNING: 1, INFO: 2}

# The panels that are not in schema.PANELS, because they settle nothing.
# Numbered here and nowhere else: a panel is named in prose, but a Check has
# to carry the number so the validation page can link to it.
RUN_PANEL = 8


class CheckError(Exception):
    """A check itself could not run."""


class Check(object):
    """One thing to say about a configuration.

    ``level``   error / warning / info
    ``panel``   the panel number that answers it, or None
    ``title``   one line, the problem itself
    ``detail``  what to do about it; may be empty
    ``key``     the dotted configuration key, when there is exactly one
    """

    __slots__ = ('level', 'panel', 'title', 'detail', 'key')

    def __init__(self, level, title, panel=None, detail='', key=''):
        if level not in _ORDER:
            raise CheckError('unknown level %r' % level)
        self.level = level
        self.title = title
        self.panel = panel
        self.detail = detail
        self.key = key

    def __repr__(self):                      # pragma: no cover - debugging
        return '<Check %s %s: %s>' % (self.level, self.key or '-', self.title)

    def __eq__(self, other):
        return (isinstance(other, Check)
                and (self.level, self.title, self.panel, self.detail,
                     self.key)
                == (other.level, other.title, other.panel, other.detail,
                    other.key))


def worst(checks):
    """The most severe level present, or ''."""
    if not checks:
        return ''
    return sorted(checks, key=lambda c: _ORDER[c.level])[0].level


def count_by_level(checks):
    out = {ERROR: 0, WARNING: 0, INFO: 0}
    for c in checks:
        out[c.level] += 1
    return out


# --------------------------------------------------------------- the checks
# Each one takes the configuration and yields Checks. They are separate
# functions so a test can exercise one without standing up the others, and so
# a check that raises cannot take the panel down with it -- see collect.

def check_schema(cfg, **_):
    """What ``RunConfig.problems`` refuses, one entry each.

    The schema's own rules, which used to reach the modeller only as a
    ConfigError from a save -- one wall of text, at the moment they were
    trying to write something else.
    """
    for msg in cfg.problems():
        key = msg.split(' ')[0].split('=')[0].strip()
        if '.' not in key or ' ' in key:
            key = ''
        yield Check(ERROR, msg, panel=_panel_of(key), key=key)


def check_switches(cfg, unsaved=(), **_):
    """Panel toggles that disagree with the file (see panelui)."""
    for switch, live, saved in unsaved or ():
        yield Check(
            ERROR,
            '%s: the panel says %s, the file says %s'
            % (switch, 'ON' if live else 'off',
               'true' if saved else 'false'),
            panel=_panel_of(switch), key=switch,
            detail='The run reads the FILE. Save the configuration, or set '
                   'the switch back to what the file says.')


def check_libmf6(cfg, runner=None, **_):
    """Is the MODFLOW 6 library where the configuration says?

    Blank is a CHOICE -- build the input files and stop -- so it is info,
    not a problem. Anything else that cannot be resolved is an error,
    because it is discovered otherwise only after the whole model is built.
    """
    given = (cfg.paths.libmf6 or '').strip()
    if not given:
        yield Check(
            INFO, 'paths.libmf6 is blank: the run stops after writing the '
            'MODFLOW 6 input files', panel=RUN_PANEL, key='paths.libmf6',
            detail='Set it to "auto" for a coupled run.')
        return
    if runner is None:
        return
    try:
        lib = runner.check_libmf6(runner.resolve_libmf6(given))
    except Exception as exc:
        yield Check(ERROR, 'paths.libmf6: %s' % exc, panel=RUN_PANEL,
                    key='paths.libmf6')
    else:
        yield Check(INFO, 'MODFLOW 6 library: %s' % lib, panel=RUN_PANEL,
                    key='paths.libmf6')


def check_forcing(cfg, dataset_dir=None, surface=None, **_):
    """The daily forcing: present, and the shape its block count implies.

    With ``run.surface`` on the series are run OUTPUT and their absence is
    normal -- MMsurf is about to write them. With it off they must already
    exist, and a truncated one is the failure this catches: the committed
    inputZONRFe_veg_d.txt held 4869 values against a 1949-day record.
    """
    if surface is None or dataset_dir is None:
        return
    if cfg.run.surface:
        yield Check(INFO, 'MMsurf will write the daily forcing', panel=2,
                    key='run.surface')
        return
    try:
        spec = surface.forcing_spec(cfg, str(dataset_dir))
        ndays = surface.check_forcing(spec, must_exist=True,
                                      nper=cfg.run.nsp or None)
    except Exception as exc:
        yield Check(ERROR, 'the daily forcing cannot be used: %s' % exc,
                    panel=2, key='run.surface',
                    detail='Switch MMsurf on to produce it, or repair the '
                           'series named above.')
    else:
        yield Check(INFO, '%d day(s) of daily forcing, every file the shape '
                    'its block count implies' % ndays, panel=2,
                    key='run.surface')


def check_state_scope(cfg, workspace=None, **_):
    """Saved state is only valid for the grid and layer set that made it.

    A spin-up written on the structured grid says nothing about a voronoi
    mesh. The run stops on this; said here it is a warning, because the
    answer may be to clear the field rather than to regenerate the state.
    """
    for key, name in (('strt_heads', 'spinup.strt_heads'),
                      ('steady_means', 'spinup.steady_means')):
        named = (getattr(cfg.spinup, key, '') or '').strip()
        if not named:
            continue
        if workspace is None:
            continue
        side = os.path.join(str(workspace), named + '.scope.json')
        if not os.path.exists(side):
            yield Check(
                WARNING,
                '%s = %r has no scope sidecar' % (name, named), panel=4,
                key=name,
                detail='It was written before the state guard existed, so '
                       'the grid it belongs to is not recorded. On grid.kind '
                       '= %r the run falls back to starting the water table '
                       'from the land surface. Regenerate it with a spin-up '
                       'on this grid, or clear the field.' % cfg.grid_kind)


def check_run_scope(cfg, **_):
    """Settings that quietly shrink a run.

    ``run.nsp`` is the one that has cost real time: set for a trial, left
    on, and the result read as though it covered the record.
    """
    if cfg.run.nsp:
        yield Check(
            WARNING, 'run.nsp = %d: the run stops after %d stress period(s) '
            'instead of covering the record' % (cfg.run.nsp, cfg.run.nsp),
            panel=2, key='run.nsp',
            detail='Useful to try a change; misleading in a result. Set it '
                   'to 0 to cover the whole record.')
    if cfg.run.build_only:
        yield Check(
            WARNING, 'run.build_only is on: the MODFLOW files are written '
            'and nothing is run', panel=RUN_PANEL, key='run.build_only')
    if not cfg.run.model:
        yield Check(
            WARNING, 'run.model is off: MMsoil and MODFLOW 6 will not run',
            panel=4, key='run.model')


def check_dataset(cfg, dataset_dir=None, **_):
    """Does the case folder this configuration names actually exist?"""
    if dataset_dir is None:
        return
    if not os.path.isdir(str(dataset_dir)):
        yield Check(ERROR, 'the dataset folder does not exist: %s'
                    % dataset_dir, panel=0, key='paths.case')


CHECKS = (check_schema, check_switches, check_dataset, check_forcing,
          check_state_scope, check_run_scope, check_libmf6)


def collect_default(cfg, unsaved=()):
    """:func:`collect` with the pieces a real application has to hand.

    Here rather than on the panel so that the validation page and the Launch
    button ask the SAME question with the SAME arguments. Anything that
    cannot be imported is simply not consulted -- the checks that depend on
    it stay silent, which is the rule for the whole module.
    """
    import mm_paths
    import marmites_config as mcfg
    try:
        import marmites_surface as surface
    except Exception:                                   # noqa: BLE001
        surface = None
    try:
        from lib import runs as runlib
        runner = runlib.import_runner()
    except Exception:                                   # noqa: BLE001
        runner = None
    try:
        ws = mcfg.state_workspace(cfg, str(mm_paths.WS_ROOT))
    except Exception:                                   # noqa: BLE001
        ws = None
    return collect(cfg, dataset_dir=mm_paths.dataset_dir(cfg.paths.case),
                   workspace=ws, unsaved=unsaved, surface=surface,
                   runner=runner)


def collect(cfg, dataset_dir=None, workspace=None, unsaved=(), surface=None,
            runner=None):
    """Every check, worst first. Never raises.

    A check that blows up becomes an error naming itself, rather than taking
    the panel down: the panel's whole job is to be reachable when something
    is wrong.
    """
    out = []
    for fn in CHECKS:
        try:
            out.extend(fn(cfg, dataset_dir=dataset_dir, workspace=workspace,
                          unsaved=unsaved, surface=surface, runner=runner))
        except Exception as exc:                        # noqa: BLE001
            out.append(Check(ERROR, '%s could not run: %s'
                             % (fn.__name__, exc)))
    return sorted(out, key=lambda c: (_ORDER[c.level],
                                      c.panel if c.panel is not None else 99))


# Which panel answers a given key. Built from the schema's own panel table, so
# it cannot drift from the sidebar.
def _panel_of(key):
    if not key or '.' not in key:
        return None
    section = key.split('.')[0]
    try:
        from lib import schema
    except ImportError:                                 # pragma: no cover
        return None
    for num, _title, _icon, switch, sections, _blurb in schema.PANELS:
        if section in (sections or ()):
            return num
        if switch == key:
            return num
    return None
