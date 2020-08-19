''' This module contains all of the budget term computations and other mathematical analysis tools'''

import numpy as np 
import pandas as pd
from math import factorial
# Local modules
from lib.structures import get_rename, get_constants


def vertIntegrate(phi, z):
    '''This function integrates our vertical profiles numerically using either trapezoidal or simpson rule '''

    Nz = len(z)
    out = {}
    if Nz % 2 == 0:
        # Nz is even so use simpsons rule
        h = (z[1] - z[0]) / Nz
        for var, val in phi.items():
            f_ev = phi[var][2:-1:2]
            f_od = phi[var][1:-1:2]
            out[var] = 1.0 / 3.0 * h * (phi[var][0] + 4.0 * (np.sum(f_od)) + 2.0 * (np.sum(f_ev)) + phi[var][-1])
    else:
        # Nz is odd so use the trapezoidal rule.
        h = (z[1] - z[0]) / Nz
        for var, val in phi.items():
            out[var] = h * (.5 * phi[var][0] + sum(phi[var][1:-1]) + .5 * phi[var][-1])
    return (out)


def get_box_dim(dfm, idx_b, idx_e):
    '''This function computes the average length and width of the Volume for use in advection calculation'''
    # compute the box dimensions using the distance covered during the
    # advection leg of the flight and the True air speed
    time = dfm['time'][idx_b['advAMP1']:idx_e['advMPB1']].values
    tas = dfm['TAS (m/s)'][idx_b['advAMP1']:idx_e['advMPB1']].values
    AB1 = np.sum((time[1:] - time[:-1]) * tas[:-1])
    time = dfm['time'][idx_b['advBMP1']:idx_e['advMPC1']].values
    tas = dfm['TAS (m/s)'][idx_b['advBMP1']:idx_e['advMPC1']].values
    BC1 = np.sum((time[1:] - time[:-1]) * tas[:-1])
    time = dfm['time'][idx_b['advCMP1']:idx_e['advMPD1']].values
    tas = dfm['TAS (m/s)'][idx_b['advCMP1']:idx_e['advMPD1']].values
    CD1 = np.sum((time[1:] - time[:-1]) * tas[:-1])
    time = dfm['time'][idx_b['advDMP1']:idx_e['advMPA1']].values
    tas = dfm['TAS (m/s)'][idx_b['advDMP1']:idx_e['advMPA1']].values
    DA1 = np.sum((time[1:] - time[:-1]) * tas[:-1])

    time = dfm['time'][idx_b['advAMP2']:idx_e['advMPB2']].values
    tas = dfm['TAS (m/s)'][idx_b['advAMP2']:idx_e['advMPB2']].values
    AB2 = np.sum((time[1:] - time[:-1]) * tas[:-1])
    time = dfm['time'][idx_b['advBMP2']:idx_e['advMPC2']].values
    tas = dfm['TAS (m/s)'][idx_b['advBMP2']:idx_e['advMPC2']].values
    BC2 = np.sum((time[1:] - time[:-1]) * tas[:-1])
    time = dfm['time'][idx_b['advCMP2']:idx_e['advMPD2']].values
    tas = dfm['TAS (m/s)'][idx_b['advCMP2']:idx_e['advMPD2']].values
    CD2 = np.sum((time[1:] - time[:-1]) * tas[:-1])
    time = dfm['time'][idx_b['advDMP2']:idx_e['advMPA2']].values
    tas = dfm['TAS (m/s)'][idx_b['advDMP2']:idx_e['advMPA2']].values
    DA2 = np.sum((time[1:] - time[:-1]) * tas[:-1])

    # These are the average distances in the box
    dy = (AB1 + CD1 + AB2 + CD2) / 4.
    dx = (BC1 + DA1 + BC2 + DA2) / 4.

    return dx, dy


def get_time_weights(df, idxb, idxe):
    '''This function computes the total time for each of the flight legs used in tend. analysis for computation of the weighted average'''
    dt = {}
    dt['AMtend'] = df['time'][idxe['advAMP1']] - df['time'][idxb['advAMP1']]
    dt['MBtend'] = df['time'][idxe['advMPB1']] - df['time'][idxb['advMPB1']]
    dt['BMtend'] = df['time'][idxe['advBMP1']] - df['time'][idxb['advBMP1']]
    dt['MCtend'] = df['time'][idxe['advMPC1']] - df['time'][idxb['advMPC1']]
    dt['CMtend'] = df['time'][idxe['advCMP1']] - df['time'][idxb['advCMP1']]
    dt['MDtend'] = df['time'][idxe['advMPD1']] - df['time'][idxb['advMPD1']]
    dt['DMtend'] = df['time'][idxe['advDMP1']] - df['time'][idxb['advDMP1']]
    dt['MAtend'] = df['time'][idxe['advMPA1']] - df['time'][idxb['advMPA1']]
    dt['tend1'] = df['time'][idxe['tend1']] - df['time'][idxb['tend1']]

    return dt


def compute_sflux(dfm, bg, ph, idx_b, idx_e):
    '''This function computes the surface flux, and amends the budget term by reference.'''
    zsurf1 = dfm.alt[idx_b['sst1']:idx_e['sst1']].mean()
    zsurf2 = dfm.alt[idx_b['sst2']:idx_e['sst2']].mean()

    wsurf1 = dfm.w[idx_b['sst1']:idx_e['sst1']].values - dfm.w[idx_b['sst1']:idx_e['sst1']].mean()
    wsurf2 = dfm.w[idx_b['sst2']:idx_e['sst2']].values - dfm.w[idx_b['sst2']:idx_e['sst2']].mean()
    for var, val in ph.items():
        bg['sflux1'][var] = np.mean(
            wsurf1 * (dfm[var][idx_b['sst1']:idx_e['sst1']].values - dfm[var][idx_b['sst1']:idx_e['sst1']].mean()))
        bg['sflux2'][var] = np.mean(
            wsurf2 * (dfm[var][idx_b['sst2']:idx_e['sst2']].values - dfm[var][idx_b['sst2']:idx_e['sst2']].mean()))


def compute_profiles(df, phi, prf, idx_b, idx_e, freq):
    '''This function computes the profiles of advection and tendency'''
    # start by determining box dimensions and average phi over each leg
    dx, dy = get_box_dim(df, idx_b, idx_e)
    for var, val in phi.items():
        phi[var]['AB1'] = (phi[var]['advAMP1'] + phi[var]['advMPB1']) / 2
        phi[var]['BC1'] = (phi[var]['advBMP1'] + phi[var]['advMPC1']) / 2
        phi[var]['CD1'] = (phi[var]['advCMP1'] + phi[var]['advMPD1']) / 2
        phi[var]['DA1'] = (phi[var]['advDMP1'] + phi[var]['advMPA1']) / 2
        phi[var]['CTR1'] = (phi[var]['AB1'] + phi[var]['BC1'] + phi[var]['CD1'] + phi[var]['DA1']) / 4

        phi[var]['AB2'] = (phi[var]['advAMP2'] + phi[var]['advMPB2']) / 2
        phi[var]['BC2'] = (phi[var]['advBMP2'] + phi[var]['advMPC2']) / 2
        phi[var]['CD2'] = (phi[var]['advCMP2'] + phi[var]['advMPD2']) / 2
        phi[var]['DA2'] = (phi[var]['advDMP2'] + phi[var]['advMPA2']) / 2
        phi[var]['CTR2'] = (phi[var]['AB2'] + phi[var]['BC2'] + phi[var]['CD2'] + phi[var]['DA2']) / 4

    # Compute Advection terms
    for var, val in phi.items():
        prf['advx1'][var] = -phi['u']['CTR1'] * (phi[var]['AB1'] - phi[var]['CD1']) / dx
        prf['advy1'][var] = -phi['v']['CTR2'] * (phi[var]['BC1'] - phi[var]['DA1']) / dx
        prf['advx2'][var] = -phi['u']['CTR1'] * (phi[var]['AB2'] - phi[var]['CD2']) / dy
        prf['advy2'][var] = -phi['v']['CTR2'] * (phi[var]['BC2'] - phi[var]['DA2']) / dy
    
    # Compute tendency terms
    for var, val in phi.items():
        prf['tend1'][var] = 2. * (val['tend2'] - val['tend1']) / (
                    1 / freq * ((idx_b['tend2'] + idx_e['tend2']) - (idx_b['tend1'] + idx_e['tend1'])))
        prf['tend2'][var] = prf['tend1'][var] * 2. * (val['tend3'] - val['tend2']) / (
                    1 / freq * ((idx_b['tend3'] + idx_e['tend3']) - (idx_b['tend2'] + idx_e['tend2'])))
        prf['AMtend'][var] = 2. * (val['advAMP2'] - val['advAMP1']) / (
                    1 / freq * ((idx_b['advAMP2'] + idx_e['advAMP2']) - (idx_b['advAMP1'] + idx_e['advAMP1'])))
        prf['MBtend'][var] = 2. * (val['advMPB2'] - val['advMPB1']) / (
                    1 / freq * ((idx_b['advMPB2'] + idx_e['advMPB2']) - (idx_b['advMPB1'] + idx_e['advMPB1'])))
        prf['BMtend'][var] = 2. * (val['advBMP2'] - val['advBMP1']) / (
                    1 / freq * ((idx_b['advBMP2'] + idx_e['advBMP2']) - (idx_b['advBMP1'] + idx_e['advBMP1'])))
        prf['MCtend'][var] = 2. * (val['advMPC2'] - val['advMPC1']) / (
                    1 / freq * ((idx_b['advMPC2'] + idx_e['advMPC2']) - (idx_b['advMPC1'] + idx_e['advMPC1'])))
        prf['CMtend'][var] = 2. * (val['advCMP2'] - val['advCMP1']) / (
                    1 / freq * ((idx_b['advCMP2'] + idx_e['advCMP2']) - (idx_b['advCMP1'] + idx_e['advCMP1'])))
        prf['MDtend'][var] = 2. * (val['advMPD2'] - val['advMPD1']) / (
                    1 / freq * ((idx_b['advMPD2'] + idx_e['advMPD2']) - (idx_b['advMPD1'] + idx_e['advMPD1'])))
        prf['DMtend'][var] = 2. * (val['advDMP2'] - val['advDMP1']) / (
                    1 / freq * ((idx_b['advDMP2'] + idx_e['advDMP2']) - (idx_b['advDMP1'] + idx_e['advDMP1'])))
        prf['MAtend'][var] = 2. * (val['advMPA2'] - val['advMPA1']) / (
                    1 / freq * ((idx_b['advMPA2'] + idx_e['advMPA2']) - (idx_b['advMPA1'] + idx_e['advMPA1'])))

def compute_constants(Wfile):
    '''This function returns global physical constants for the flight including lat,lon of flight volume waypoints'''
    constants = get_constants()   
    # from this compute the orientation (heading) of the box
    # Set the Waypoints of the box 
    latwp={}
    lonwp={}
    
    WPlist=['A','B','C','D']
    tmplat=[]
    tmplon=[]

    with open(str(Wfile),"r") as f:
        for line in f.readlines():
            line=line.rstrip()
            lat,lon=line.split(",")
            tmplat.append(float(lat))
            tmplon.append(float(lon))
            
    #TODO: Use lamb functions for this 
    for i in range(4):
        key=WPlist[i]
        latwp[key]=tmplat[i]
        lonwp[key]=tmplon[i]
        
    constants['headAB'] = np.arctan(constants['kmpdeglat']*(latwp['B']-latwp['A'])/(constants['kmpdeglon']*(lonwp['B']-lonwp['A'])))+np.deg2rad(270.)
    wp = {'lat' : latwp, 'lon' : lonwp}
    return constants, wp


def init_cbfile(fname, constants):
    ''' This function initializes the cabin file for analysis'''
    df = pd.read_csv(fname, sep=",", index_col="Rec ID.",na_values=-9999)

    df.rename(columns=get_rename(), inplace=True)
    # Derived Quantities:
    # time in seconds
    df['time'] = pd.Series(df['UTC time (HH.hhhhh)'] * 3600.)
    # x-component of wind velocity positive from A to B
    df['u'] = pd.Series(df.ws * np.sin(constants['headAB'] - np.deg2rad(df.wdir)), index=df.index)
    # y-velocity
    df['v'] = pd.Series(df.ws * np.cos(constants['headAB'] - np.deg2rad(df.wdir)), index=df.index)
    # Liquid water
    df['ql'] = pd.Series(df['HPVM LWC'] / constants['rho_air'], index=df.index)
    # Total Water(g/kg)
    df['qt'] = pd.Series(df.ql + df.qv, index=df.index)
    # moist static energy
    # We'll break this into three terms to plot individually later
    T1 = pd.Series(constants['Cp'] * 1.0e3 * (df['T'] + 273.15), index=df.index)
    T2 = pd.Series(constants['g'] * df['alt'], index=df.index)
    T3 = pd.Series(constants['Lv'] * df['r'], index=df.index)
    df['se_t1'] = T1
    df['se_t2'] = T2
    df['se_t3'] = T3
    df['se'] = T1 + T2 + T3

    #Now we need to do a cleanup of any nan values in our derived quantities so that they dont fuck with the data
    flag_nans=['time', 'u', 'v','ql','qt','se_t1','se_t2','se_t3']
    df.dropna(how='any',subset=flag_nans,inplace=True)
    return df

def get_index(fname, dfindx):
    '''This function initializes the indexes which divide up the cabin file into its individual flight legs'''
    idxb = {}
    idxe = {}
    with open(str(fname), "r") as f:
        for line in f.readlines():
            line = line.rstrip()
            key, strt, end = line.split(",")
            idxb[key] = dfindx[dfindx.searchsorted(int(strt))]
            idxe[key] = dfindx[dfindx.searchsorted(int(end))]
    return idxb, idxe


def create_histogram(dim, data):
    '''This function creates a histogram for the events provided and includes a reference list of indexes sorted to bins

    This function creates a histogram of "events" with respect to some scale dataset. it expects a dimensions dict
    with at least dn or h_freq included, from which it will compute the remaining dimensions of the histogram based
    on the scale data. the data argument should be a pd series object containing the data to be sorted'''


    # compute dimensions of histogram
    h_dim = {
        'st': dim.get('st'),
        'ed': dim.get('ed'),
        'dn': dim.get('dn'),
        'h_freq': dim.get('h_freq')
    }

    if not h_dim['st']:
        h_dim['st'] = np.floor(np.min(data))
    if not h_dim['ed']:
        h_dim['ed'] = np.ceil(np.max(data))
    if h_dim['dn'] and not h_dim['h_freq']:
        h_dim['h_freq'] = 1.0/h_dim['dn']
    elif h_dim['h_freq'] and not h_dim['dn']:
        h_dim['dn'] = 1.0/h_dim['h_freq']
    else:
        print("Malformed dimension object!")
        raise ValueError

    h_dim['N'] = int(h_dim['h_freq']*(h_dim['ed'] - h_dim['st']))

    h_edges = [h_dim['st'] + x*h_dim['dn'] for x in range(h_dim['N'])]
    h_edges.append(h_dim['ed'])
    # These lists will hold the data from each bin 
    counts = []
    h_index= []
    for i in range(h_dim['N']):

        l_edge = h_edges[i]
        r_edge = h_edges[i+1]
        tmp_ind = []
        _counts = 0
        for index, entry in data.iteritems():


            if l_edge <= entry < r_edge:
                _counts += 1
                tmp_ind.append(index)
        counts.append(_counts)
        h_index.append(tmp_ind)

    histogram={
        'edges': h_edges,
        'counts':counts,
        'index': h_index,
        'dim': h_dim
    }

    return histogram

