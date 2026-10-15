# -*- coding: utf-8 -*-
"""Shared by every test file: ``lib`` is the APP's package, in any order.

The pages, panelui and many tests write ``from lib import schema``, i.e.
``code/app/lib`` imported as the top-level package ``lib``. On Windows
pywin32 puts ``site-packages\\win32`` on sys.path, and the ``lib`` folder in
it imports as a NAMESPACE package ``lib``. A test that reached ``from lib
import ...`` while code/app was not on sys.path -- checks._panel_of does, and
swallows the ImportError -- left that namespace in sys.modules, and every
later ``from lib import`` of the session failed, whatever it put on sys.path
first. test_app_lib.py run whole failed on ``from lib import panelui``
(2026-10-15); the full suite passed only because another file's module-level
``from lib import`` happened to be collected first.

So the app's ``lib`` is registered here, by its file, before any test file
is collected: sys.path is left alone, and every run, a single test or the
whole suite, imports the same package.
"""

import importlib.util
import os
import sys

_LIB = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'app',
                                    'lib'))

if getattr(sys.modules.get('lib'), '__file__', None) is None:
    _spec = importlib.util.spec_from_file_location(
        'lib', os.path.join(_LIB, '__init__.py'),
        submodule_search_locations=[_LIB])
    _mod = importlib.util.module_from_spec(_spec)
    sys.modules['lib'] = _mod
    _spec.loader.exec_module(_mod)
