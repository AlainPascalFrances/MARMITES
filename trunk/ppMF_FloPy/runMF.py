#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      alf
#
# Created:     02-05-2011
# Copyright:   (c) alf 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python
import sys
sys.path.append(r'E:\00code\MARMITES\trunk\ppMF_FloPy')
import ppMODFLOW_flopy as ppMF

MF_ws = r'E:\00code\00ws\zz_TESTS\MARMITESv2_r13c6l2_REF\MF_ws'

SP_d, nrow, ncol, delr, delc, nlay, perlen, nper, top, hnoflo, hdry, ibound, laytyp, h_MF, cbc, top_array, inputFileMF_fn, lenuni = ppMF.ppMF(MF_ws, rch_input = 0.0001, rch_dft = 0.0001)

print "\nMODFLOW done!"

