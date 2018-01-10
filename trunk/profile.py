#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      alf
#
# Created:     13-05-2011
# Copyright:   (c) alf 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python
import sys
sys.path.append(r'E:\00code\MARMITES\trunk')

import cProfile
cProfile.run('import startMARMITES_v3', 'profile.tmp')

import pstats
p = pstats.Stats('profile.tmp')
p.sort_stats('cumulative').print_stats(50)