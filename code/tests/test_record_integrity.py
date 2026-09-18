# -*- coding: utf-8 -*-
"""A corrupt cell stops the read; it does not delete the row.

La Mata's meteorological record and irrigation series both carried
``#VALUE!`` where the timestamp for 2013-08-16 13:00 belongs -- an Excel
error written into the file. numpy's ``loadtxt`` treats a leading ``#``
as a comment, so it stripped the line to nothing and DROPPED it, leaving
both series an hour short and saying so only through a warning about
``max_rows``. The measurements for that hour were intact all along.
"""

import importlib.util
import os
import re
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
READERS = (os.path.join(CODE, 'MARMITESsurf', 'startMARMITESsurface.py'),
           os.path.join(CODE, 'ppMF_FloPy', 'ppMODFLOW_flopy_v3.py'))


def _rows_of(path):
    """Lift _rows out of a module that cannot be imported standalone."""
    src = open(path, encoding='utf-8').read()
    block = re.search(r'EXCEL_ERRORS = .*?\n    return out\n', src, re.S)
    assert block, '%s has no _rows guard' % os.path.basename(path)
    ns = {'os': os}
    exec(block.group(0), ns)
    return ns['_rows'], ns['EXCEL_ERRORS']


@pytest.mark.parametrize('reader', READERS, ids=os.path.basename)
def test_an_excel_error_stops_the_read(reader, tmp_path):
    rows, errors = _rows_of(reader)
    good = tmp_path / 'clean.txt'
    good.write_text('date\tvalue\n2013-08-16\t1.0\n2013-08-17\t2.0\n',
                    encoding='utf-8')
    assert len(rows(str(good))) == 3

    bad = tmp_path / 'broken.txt'
    bad.write_text('date\tvalue\n2013-08-16\t1.0\n#VALUE!\t2.0\n',
                   encoding='utf-8')
    with pytest.raises(ValueError) as e:
        rows(str(bad))
    assert 'line 3' in str(e.value)
    assert '#VALUE!' in str(e.value)
    assert 'broken.txt' in str(e.value)


@pytest.mark.parametrize('reader', READERS, ids=os.path.basename)
def test_every_excel_error_is_caught_not_only_the_one_we_met(reader, tmp_path):
    rows, errors = _rows_of(reader)
    for token in errors:
        f = tmp_path / 'x.txt'
        f.write_text('h\n1.0\n%s\n' % token, encoding='utf-8')
        with pytest.raises(ValueError) as e:
            rows(str(f))
        assert token in str(e.value)


@pytest.mark.parametrize('reader', READERS, ids=os.path.basename)
def test_a_trailing_newline_is_still_just_a_trailing_newline(reader,
                                                             tmp_path):
    """The blank-line filter was right; it was only not the whole story."""
    rows, _e = _rows_of(reader)
    f = tmp_path / 'x.txt'
    f.write_text('h\n1.0\n2.0\n\n', encoding='utf-8')
    assert len(rows(str(f))) == 3


@pytest.mark.parametrize('reader', READERS, ids=os.path.basename)
def test_hash_no_longer_means_comment_in_a_data_file(reader):
    """Belt as well as braces: if a corrupt cell ever slips past the guard
    above, loadtxt must fail to parse it rather than make the row vanish."""
    src = open(reader, encoding='utf-8').read()
    # by LINE: a regex stopping at the first ')' closes on _rows(...) and
    # never sees the keywords that follow it
    reads = [ln for ln in src.splitlines() if 'np.loadtxt(_rows(' in ln]
    assert reads, 'no guarded data read in %s' % os.path.basename(reader)
    for ln in reads:
        assert 'comments = None' in ln or 'comments=None' in ln, (
            'a data read still lets "#" start a comment: %s' % ln.strip())
