# -*- coding: utf-8 -*-
"""Run one PEST tool for the front-end and record how it ended (WP7
front-end, phase 4; app/lib/pest_jobs.py starts it, detached).

    python pest_job.py <job folder> <tool.py> [tool arguments...]

The tool's output goes where this script's goes (the job's job.log); when
the tool returns, ``result.json`` in the job folder records its exit code
-- which another process cannot read back on Windows, so the page would
otherwise have to guess success from the log.
"""

import json
import os
import subprocess
import sys
import time


def main(argv):
    if len(argv) < 2:
        sys.exit(__doc__)
    job_dir, tool, args = argv[0], argv[1], argv[2:]
    t0 = time.time()
    print('job: %s %s' % (os.path.basename(tool), ' '.join(args)), flush=True)
    sys.stdout.flush()
    # the log handed on EXPLICITLY: a process started detached has no
    # console, and on Windows its child inherits no output handle unless
    # one is given -- the tool's output was lost (2026-10-20)
    rc = subprocess.call([sys.executable, '-u', tool] + list(args),
                         stdout=sys.stdout, stderr=subprocess.STDOUT)
    res = {'returncode': int(rc),
           'finished': time.strftime('%Y-%m-%d %H:%M:%S'),
           'minutes': round((time.time() - t0) / 60.0, 1)}
    tmp = os.path.join(job_dir, 'result.json.tmp')
    with open(tmp, 'w', encoding='utf-8') as fh:
        json.dump(res, fh, indent=1)
    os.replace(tmp, os.path.join(job_dir, 'result.json'))
    print('job ended: exit %d after %.1f min' % (rc, res['minutes']),
          flush=True)
    return rc


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
