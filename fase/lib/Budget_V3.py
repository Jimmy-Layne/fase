
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 13:55:21 2016

@author: James Layne, and Thijs
"""


import matplotlib.pyplot as plt
import numpy as np 

from scipy.interpolate import griddata
from copy import deepcopy
# local modules
from lib.flight_info import get_flight
import lib.structures as stct
import lib.computations as cmp
import lib.entrainment as ent


def volumeBudget(date):
    """This is the main function which commands the computation of the qt budget of the volume"""
    # Begin by getting flight information and physical constants based on the date provided
    flight = get_flight(date)
    constants, wp = cmp.compute_constants(flight['Wfile'])
    # Initialize dict structures
    phi = stct.get_phi()
    profiles = stct.get_profiles()
    budget = stct.get_budget()
    # Initialize the main dataframe
    df = cmp.init_cbfile(flight['Cfile'], constants)
    # Initialize the Indices of the flight legs
    idxb, idxe = cmp.get_index(flight['Ifile'], df.index)

    # Here is our empty height profile
    zi = np.linspace(0., flight['Ctop'], flight['Nz'])

    # Now we perform the linear interpolation of the conserved variables.
    for key,beg in idxb.items():
        end = idxe[key]
        for var,val in phi.items():
            phi[var][key] = griddata(df.alt[beg:end].values, df[var][beg:end].values, zi, method="nearest")
    
    # >>>Functions which edit By Reference
    # First call the function to compute the two vertical profile legs Tendency and advection
    cmp.compute_profiles(df, phi, profiles,idxb, idxe, constants['sample_freq'])
    # Now compute Surface Flux
    cmp.compute_sflux(df, budget, phi, idxb, idxe)
    # >>>> End reference edits

    # Create Entrainment Profile object
    entr = ent.EntrAnalysis(date=date, df=df)

    for key in phi.keys():
        # Entrainment budget term is given by the average flux
        budget['eflux1'][key] = entr.avg_flux[key]

    # Now integrate all of the profiles to create the budget terms for those flight legs
    for term,val in profiles.items():
        budget[term] = cmp.vertIntegrate(val, zi)
    
    # To compute tendency term we take a weighted average of all the legs used to measure tend
    dt = cmp.get_time_weights(df, idxb, idxe)

    for var in phi.keys():
        num=0.0
        den=0.0
        for k in dt.keys():
            num=num+dt[k]*budget[k][var]
            den=den+dt[k]
        budget['FullTend'][var]=num/den
        # The residual is the remaining quantity which closes the budget.
        budget['FTres'][var] = budget['FullTend'][var] - \
                               (budget['advx1'][var]+budget['advy1'][var]+budget['eflux1'][var]+budget['sflux1'][var])
    
    return phi, zi, profiles, budget




