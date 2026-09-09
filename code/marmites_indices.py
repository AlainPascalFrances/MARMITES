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
}

# per-soil-layer fluxes (MM_S dataset, last axis)
INDEX_MM_SOIL = {
    'iEsoil': 0, 'iTsoil': 1, 'iSsoil_pc_s': 2, 'iRsoil': 3, 'iExf': 4,
    'idSsoil_s': 5, 'iSsoil': 6, 'iSAT': 7, 'iMB_s': 8,
}
