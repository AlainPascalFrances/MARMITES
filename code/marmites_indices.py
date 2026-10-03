# -*- coding: utf-8 -*-
"""Shared definitions of the MARMITES HDF5 output-array indices (Phase 2).

These maps name the last-axis columns of the ``MM`` (per-cell fluxes) and
``MM_S`` (per-soil-layer fluxes) HDF5 datasets. They were previously a
literal block inside startMARMITES_v3.py and re-declared in tests; keeping
them here removes that duplication and gives one source of truth used by
the driver, the soil model callers, and the test suite.
"""

__all__ = ['INDEX_MM', 'INDEX_MM_SOIL']

# per-cell fluxes (MM dataset, last axis)
INDEX_MM = {
    'iP': 0, 'iPT': 1, 'iPE': 2, 'iPe': 3, 'iSsurf': 4, 'iRo': 5, 'iEXFg': 6,
    'iEow': 7, 'iMB': 8, 'iEi': 9, 'iEo': 10, 'iEg': 11, 'iTg': 12, 'idSsurf': 13,
    'iETg': 14, 'iETsoil': 15, 'iSsoil_pc': 16, 'idSsoil': 17, 'iperc': 18,
    'ihcorr': 19, 'idgwt': 20, 'iuzthick': 21, 'iI': 22, 'iMBsurf': 23,
    # WP2, appended so the legacy MODFLOW-NWT reference (24 columns) still
    # indexes as before: UZF's PET demand, its ACTUAL ET (from the MF6
    # budget), the total ET of the five sources, and the rejected
    # infiltration UZF returned to the soil column from below
    'iPETuzf': 24, 'iETuzf': 25, 'iETtot': 26, 'iRejInf': 27,
    # WP2 row 2: the open-water evaporation kept apart by package -- the
    # stream reaches (SFR SIMEVAP) and the ponds (LAK EVAP). iEow stays
    # their sum; how much a pond loses is a result in its own right.
    'iEow_sfr': 28, 'iEow_lak': 29,
    # what MF6's UZF ET routine removed BEYOND the demand, within its own
    # wave-merging tolerance: a numerical storage loss of UZF, not ET
    'iETuzf_num': 30,
    # WP5, the CRR cascade: the run-on a soil column receives from upslope,
    # the part of it the soil took (inside iI), and the runoff the cascade
    # evaporates -- the 1 - beta share and a sink's water -- booked at the
    # cell it left and counted in iETtot. All zero with [crr] off.
    'iRunon': 31, 'iEcrr': 32, 'iReinf': 33,
}

# per-soil-layer fluxes (MM_S dataset, last axis)
INDEX_MM_SOIL = {
    'iEsoil': 0, 'iTsoil': 1, 'iSsoil_pc_s': 2, 'iRsoil': 3, 'iExf': 4,
    'idSsoil_s': 5, 'iSsoil': 6, 'iSAT': 7, 'iMB_s': 8,
}
