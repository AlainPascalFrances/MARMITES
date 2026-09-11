# -*- coding: utf-8 -*-
"""WP1d -- editing a run configuration from the panels.

DELIBERATELY STREAMLIT-FREE. The panels collect widget values into a plain
dict and hand it here; everything that could actually be wrong -- a change
that would make the configuration invalid, a table row with a missing field,
a save that silently does nothing -- is decided and tested in this module.

The rule is WP0's, unchanged: a configuration is VALIDATED BEFORE IT IS
WRITTEN. A file on disk that the model would refuse is worse than an error in
the browser, because it is discovered hours later by a run.
"""

import copy
import dataclasses
import os

__all__ = ['EditError', 'changes', 'apply_changes', 'save', 'table_rows',
           'set_table', 'add_row', 'drop_row', 'row_defaults']


class EditError(Exception):
    """An edit cannot be applied to the configuration."""


def _get(cfg, dotted):
    obj = cfg
    for part in dotted.split('.'):
        if not hasattr(obj, part):
            raise EditError('no such key: %s' % dotted)
        obj = getattr(obj, part)
    return obj


def changes(cfg, edited):
    """``section.key=value`` strings for the entries that actually differ.

    Compared as STRINGS, the same way ``apply_overrides`` coerces them back,
    so a float that survives a round trip through a widget does not register
    as a change and nothing is written for a form no one touched.
    """
    out = []
    for dotted, val in sorted(edited.items()):
        try:
            cur = _get(cfg, dotted)
        except EditError:
            raise
        if isinstance(cur, bool):
            same = bool(val) == cur
        else:
            same = str(val) == str(cur)
        if not same:
            out.append('%s=%s' % (dotted, val))
    return out


def apply_changes(cfg, edited):
    """Apply the edits to a COPY, so a rejected form leaves the original alone.

    Returns ``(new_cfg, applied)``. Raises ``ConfigError`` from the schema if
    the result would be invalid -- which is the point: the panels show that
    message instead of writing the file.
    """
    todo = changes(cfg, edited)
    new = copy.deepcopy(cfg)
    if todo:
        new.apply_overrides(todo, echo=False)
    return new, todo


def save(cfg, path, edited=None):
    """Validate, then write. Returns (applied, config_hash).

    ``write_toml`` writes the RESOLVED configuration, so the file that comes
    back is what the model will actually see -- provenance, not a draft.
    """
    new, applied = apply_changes(cfg, edited or {})
    new.validate()
    d = os.path.dirname(os.path.abspath(path))
    if not os.path.isdir(d):
        raise EditError('no such folder: %s' % d)
    new.write_toml(path)
    return applied, new.config_hash()


# ---------------------------------------------------------------- tables

def _resolve_table(cfg, dotted):
    section, name = dotted.split('.', 1)
    sec = getattr(cfg, section, None)
    if sec is None or not hasattr(sec, name):
        raise EditError('no such table: %s' % dotted)
    elements = getattr(type(sec), '_ELEMENTS', {})
    if name not in elements:
        raise EditError('%s is not an array of tables' % dotted)
    return sec, name, elements[name]


def row_defaults(cfg, dotted):
    """A new row, with the dataclass defaults already in it."""
    _sec, _name, element = _resolve_table(cfg, dotted)
    return {f.name: getattr(element(), f.name)
            for f in dataclasses.fields(element)}


def table_rows(cfg, dotted):
    """The table as a list of plain dicts, ready for a data editor."""
    sec, name, element = _resolve_table(cfg, dotted)
    names = [f.name for f in dataclasses.fields(element)]
    return [{n: getattr(row, n) for n in names} for row in getattr(sec, name)]


def set_table(cfg, dotted, rows):
    """Replace a table from a list of dicts, on a COPY.

    Every row must carry every field: a data editor that drops a column would
    otherwise quietly reset it to the dataclass default, which for a
    vegetation type means silently swapping in FAO-56 grass.
    """
    sec, name, element = _resolve_table(cfg, dotted)
    fields_ = [f.name for f in dataclasses.fields(element)]
    new = copy.deepcopy(cfg)
    built = []
    for k, row in enumerate(rows):
        missing = [f for f in fields_ if f not in row]
        if missing:
            raise EditError('%s row %d is missing %s'
                            % (dotted, k + 1, ', '.join(missing)))
        extra = [f for f in row if f not in fields_]
        if extra:
            raise EditError('%s row %d has unknown field(s) %s -- valid: %s'
                            % (dotted, k + 1, ', '.join(sorted(extra)),
                               ', '.join(fields_)))
        built.append(element(**{f: row[f] for f in fields_}))
    setattr(getattr(new, dotted.split('.', 1)[0]), name, built)
    new.validate()
    return new


def add_row(cfg, dotted):
    rows = table_rows(cfg, dotted)
    rows.append(row_defaults(cfg, dotted))
    return set_table(cfg, dotted, rows)


def drop_row(cfg, dotted, index):
    rows = table_rows(cfg, dotted)
    if not (0 <= index < len(rows)):
        raise EditError('%s has no row %d' % (dotted, index + 1))
    del rows[index]
    return set_table(cfg, dotted, rows)
