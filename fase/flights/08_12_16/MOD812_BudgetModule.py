
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 13:55:21 2016

@author: Thijs
"""


import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd
import math
from scipy.interpolate import griddata

import pdb




def distance_on_earth(lat1, long1, lat2, long2):
     
    # Convert latitude and longitude to
    # spherical coordinates in radians.
    degrees_to_radians = math.pi/180.0
     
    # phi = 90 - latitude
    phi1 = (90.0 - lat1)*degrees_to_radians
    phi2 = (90.0 - lat2)*degrees_to_radians
     
    # theta = longitude
    theta1 = long1*degrees_to_radians
    theta2 = long2*degrees_to_radians
     
    # Compute spherical _ from spherical coordinates.
     
    # For two locations in spherical coordinates
    # (1, theta, phi) and (1, theta', phi')
    # cosine( arc length ) =
    # sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length
     
    cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) +
    math.cos(phi1)*math.cos(phi2))
    arc = math.acos( cos )
     
    # Remember to multiply arc by the radius of the earth
    # in your favorite set of units to get length.
    arc*=6.373e6
    return arc







def headingsplot(fname):
    #Commenting these out for now
    #flights = np.array([fname])
    #fname = flights[0]

    df = pd.read_csv(fname,sep = ",",index_col = "Rec ID.", na_values=[-9999])
    df.rename(columns={'Wind Speed (m/s)': 'ws','Wind dir(Deg)': 'wdir','Palt (m)':'alt','Heading (deg)':'heading','MR-h2o (g/Kg)':'qv','Vert. Wind (m/s)':'w','L5200 CO2 (V)':'co2','Thetae':'thetae'},inplace=True)

    pdb.set_trace()
    #This plots all the different plots on the same figure
    f, ax = plt.subplots(4,sharex=True)
    ax[0].plot(df.index,df.Lat.values)
    ax[1].plot(df.index,df.Long.values)
    ax[2].plot(df.index,df.alt.values)
    ax[3].plot(df.index,df.heading.values)
    
    #This plots them each on separate figures
#    fig1=plt.figure()
#    ax1=fig1.add_subplot(111)
#    ax1.plot(df.index,df.Lat.values)
#    plt.title("Latitude plot")
#    
#    fig2=plt.figure()
#    ax2=fig2.add_subplot(111)
#    ax2.plot(df.index,df.Long.values)
#    plt.title("Longitude plot")    
#
#    fig3=plt.figure()
#    ax3=fig3.add_subplot(111)
#    ax3.plot(df.index,df.alt.values)
#    plt.title("Altitude plot")
#    
#    fig4=plt.figure()
#    ax4=fig4.add_subplot(111)
#    ax4.plot(df.index,df.heading.values)
#    plt.title("Heading plot")
    
        
    
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.zeros(y.shape)
    for i in range (y.shape[0]):
        y_smooth[i,:] = np.convolve(y[i,:], box, mode='same')
    return y_smooth

def massBudget(Cfile,Ifile,Wfile,Nz=50,CT=800.):
    #Start by hard coding some constants
    #Air density, probably change this later  
    rho_air = 1.225 
    #km per degree for lat and lon
    kmpdeglat = 111.
    kmpdeglon = 90.
     
    ##################################################
    #This section sets the waypoints of the box
    ################################################
    
    latwp={}
    lonwp={}
    
    WPlist=['A','B','C','D']
    tmplat=[]
    tmplon=[]
    with open(Wfile,"r") as f:
        for line in f.readlines():
            line=line.rstrip()
            lat,lon=line.split(",")
            tmplat.append(float(lat))
            tmplon.append(float(lon))
    for i in range(4):
        key=WPlist[i]
        latwp[key]=tmplat[i]
        lonwp[key]=tmplon[i]
        
    #from this compute the heading
    headAB = np.arctan(kmpdeglat*(latwp['B']-latwp['A'])/(kmpdeglon*(lonwp['B']-lonwp['A'])))+np.deg2rad(270.)




    #########################################################################
    #This section reads and organizes the main dataframe 
    
    df = pd.read_csv(Cfile,sep = ",",index_col = "Rec ID.", na_values=[-9999])
    #Some of the headers are really inconveniently named. let's change that
    df.rename(columns={'Wind Speed (m/s)': 'ws','Wind dir(Deg)': 'wdir','Palt (m)':'alt','Heading (deg)':'heading','MR-h2o (g/Kg)':'qv','Vert. Wind (m/s)':'w','L5200 CO2 (V)':'co2','Thetae':'thetae'},inplace=True)
    #From the data file we'll compute a few easy quantities first
    #Time in seconds
    df['time'] = pd.Series(df['UTC time (HH.hhhhh)']*3600.)
    # x-velocity positive from A to B
    df['u'] = pd.Series(df.ws*np.sin(headAB-np.deg2rad(df.wdir)),index=df.index)
    #y-velocity
    df['v'] = pd.Series(df.ws*np.cos(headAB-np.deg2rad(df.wdir)),index=df.index)
    #Liquid water
    df['ql'] = pd.Series(df['HPVM LWC']/rho_air,index=df.index)
    #Total Water(g/kg?)
    df['qt'] = pd.Series(df.ql+df.qv,index=df.index)
    
    #
    #In this section we set the indices marking the start and end of the individual legs of the flight.
    #
    idxb={}
    idxe={}
    with open(Ifile,"r") as f:
        for line in f.readlines():
            line=line.rstrip()
            key,strt,end=line.split(",")
            idxb[key]=df.index[df.index.searchsorted(int(strt))]
            idxe[key]=df.index[df.index.searchsorted(int(end))]
    
    
    ######################################################################   
    #                       Data Processing
    ####################################################################
    

    #This is a height profile.
    zi = np.linspace(0.,CT,Nz)
    
    #compute the box dimensions using the distance covered during the 
    #advection leg of the flight and the True air speed
    time =  df['time'][idxb['advAMP1']:idxe['advMPB1']].values
    tas =   df['TAS (m/s)'][idxb['advAMP1']:idxe['advMPB1']].values
    AB1 = np.sum((time[1:]-time[:-1])*tas[:-1])
    time =  df['time'][idxb['advBMP1']:idxe['advMPC1']].values
    tas =   df['TAS (m/s)'][idxb['advBMP1']:idxe['advMPC1']].values
    BC1 = np.sum((time[1:]-time[:-1])*tas[:-1])
    time =  df['time'][idxb['advCMP1']:idxe['advMPD1']].values
    tas =   df['TAS (m/s)'][idxb['advCMP1']:idxe['advMPD1']].values
    CD1 = np.sum((time[1:]-time[:-1])*tas[:-1])
    time =  df['time'][idxb['advDMP1']:idxe['advMPA1']].values
    tas =   df['TAS (m/s)'][idxb['advDMP1']:idxe['advMPA1']].values
    DA1 = np.sum((time[1:]-time[:-1])*tas[:-1])
    
    time =  df['time'][idxb['advAMP2']:idxe['advMPB2']].values
    tas =   df['TAS (m/s)'][idxb['advAMP2']:idxe['advMPB2']].values
    AB2 = np.sum((time[1:]-time[:-1])*tas[:-1])
    time =  df['time'][idxb['advBMP2']:idxe['advMPC2']].values
    tas =   df['TAS (m/s)'][idxb['advBMP2']:idxe['advMPC2']].values
    BC2 = np.sum((time[1:]-time[:-1])*tas[:-1])
    time =  df['time'][idxb['advCMP2']:idxe['advMPD2']].values
    tas =   df['TAS (m/s)'][idxb['advCMP2']:idxe['advMPD2']].values
    CD2 = np.sum((time[1:]-time[:-1])*tas[:-1])
    time =  df['time'][idxb['advDMP2']:idxe['advMPA2']].values
    tas =   df['TAS (m/s)'][idxb['advDMP2']:idxe['advMPA2']].values
    DA2 = np.sum((time[1:]-time[:-1])*tas[:-1])
    
    #These are the average distances in the box
    dy = (AB1+CD1+AB2+CD2)/4.
    dx = (BC1+DA1+BC2+DA2)/4.
    
    #This takes count of the conserved variables that will be analyzed in the budget
    #phi holds the linearly interpolated data based on the profile zi.
    #What would be the impact of resolving zi smaller??
    phi={'u':{},'v':{},'w':{},'qv':{},'ql':{},'qt':{},'thetae':{},'co2':{}}
    profiles={'tend1':{},'tend2':{},'advx1':{},'advx2':{},'advy1':{},'advy2':{},'res1':{},'res2':{},'fdiv1':{},'fdiv2':{}}
    Budget={'tend1':{},'tend2':{},'advx1':{},'advx2':{},'advy1':{},'advy2':{},'res1':{},'res2':{},'fdiv1':{},'fdiv2':{},'eflux1':{},'eflux2':{},'sflux1':{},'sflux2':{}}
    for key,beg in idxb.items():
        end = idxe[key]
        for var,val in phi.items():
            phi[var][key] = griddata(df.alt[beg:end].values,df[var][beg:end].values,zi,method="nearest")

    
    #    u[key] = griddata(df.alt[beg:end].values,df.u[beg:end].values,zi)
    #    v[key] = griddata(df.alt[beg:end].values,df.v[beg:end].values,zi)
    #    w[key] = griddata(df.alt[beg:end].values,df.w[beg:end].values,zi)
    #    qv[key] = griddata(df.alt[beg:end].values,df.qv[beg:end].values,zi)
    #    ql[key] = griddata(df.alt[beg:end].values,df.ql[beg:end].values,zi)
    #    qt[key] = griddata(df.alt[beg:end].values,df.qt[beg:end].values,zi)
    #    thetae[key] = griddata(df.alt[beg:end].values,df.Thetae[beg:end].values,zi)
    #    co2[key] = griddata(df.alt[beg:end].values,df.co2[beg:end].values,zi)
    
    
    #Averages over the advection leg, used in finite difference computation
    for var,val in phi.items():    
        phi[var]['AB1'] = (phi[var]['advAMP1']+phi[var]['advMPB1'])/2
        phi[var]['BC1'] = (phi[var]['advBMP1']+phi[var]['advMPC1'])/2
        phi[var]['CD1'] = (phi[var]['advCMP1']+phi[var]['advMPD1'])/2
        phi[var]['DA1'] = (phi[var]['advDMP1']+phi[var]['advMPA1'])/2
        phi[var]['CTR1'] =(phi[var]['AB1']+phi[var]['BC1']+phi[var]['CD1']+phi[var]['DA1'])/4    
    
        phi[var]['AB2'] = (phi[var]['advAMP2']+phi[var]['advMPB2'])/2
        phi[var]['BC2'] = (phi[var]['advBMP2']+phi[var]['advMPC2'])/2
        phi[var]['CD2'] = (phi[var]['advCMP2']+phi[var]['advMPD2'])/2
        phi[var]['DA2'] = (phi[var]['advDMP2']+phi[var]['advMPA2'])/2
        phi[var]['CTR2'] = (phi[var]['AB2']+phi[var]['BC2']+phi[var]['CD2']+phi[var]['DA2'])/4    
    
    #The tendency term is the time rate of change of the desired quantity
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    ax1.plot(zi,phi['qt']['tend1'],'ro')
#    ax1.plot(zi,phi['qt']['tend2'],'bo')
#    plt.show()
#    
    #Tendency    
    freq=10.
    for var,val in phi.items():
        profiles['tend1'][var] = 2.*(val['tend2']-val['tend1'])/(1/freq*((idxb['tend2']+idxe['tend2'])-(idxb['tend1']+idxe['tend1'])))
        profiles['tend2'][var] = profiles['tend1'][var]*2.*(val['tend3']-val['tend2'])/(1/freq*((idxb['tend3']+idxe['tend3'])-(idxb['tend2']+idxe['tend2'])))




    
    #Advection
    for var,val in phi.items():
        profiles['advx1'][var]=-phi['u']['CTR1']*(phi[var]['AB1']-phi[var]['CD1'])/dx
        #profiles['advx2'][var]=-phi['u']['CTR2']*(phi[var]['BC1']-phi[var]['DA1'])/dx
        profiles['advy1'][var]=-phi['v']['CTR1']*(phi[var]['AB2']-phi[var]['CD2'])/dy
        #profiles['advy2'][var]=-phi['v']['CTR2']*(phi[var]['BC2']-phi[var]['DA2'])/dy


    
    #Entrainment
    ztop1 = df.alt[idxb['entrW1']:idxe['entrW1']].mean()
    #ztop2 = df.alt[idxb['entrW2']:idxe['entrW2']].mean()
    wentr1 = df.w[idxb['entrW1']:idxe['entrW1']].values - df.w[idxb['entrW1']:idxe['entrW1']].mean()
    #wentr2 = df.w[idxb['entrW2']:idxe['entrW2']].values - df.w[idxb['entrW2']:idxe['entrW2']].mean()
    for var,val in phi.items():
        Budget['eflux1'][var] = -np.mean(wentr1 * (df[var][idxb['entrW1']:idxe['entrW1']].values - df[var][idxb['entrW1']:idxe['entrW1']].mean()))
        #Budget['eflux2'][var] = -np.mean(wentr2 * (df[var][idxb['entrW2']:idxe['entrW2']].values - df[var][idxb['entrW2']:idxe['entrW2']].mean()))
    #wentr1 = df.w[idxb['entrW1']:idxe['entrW1']].mean()
    #wentr2 = df.w[idxb['entrW2']:idxe['entrW2']].mean()
    #for var,val in phi.items():
    #    entrflux1[var] = wentr1*(df[var][idxb['entrST1']:idxe['entrST1']].max()-df[var][idxb['entrST1']:idxe['entrST1']].min())
    #    entrflux2[var] = wentr2*(df[var][idxb['entrST2']:idxe['entrST2']].max()-df[var][idxb['entrST2']:idxe['entrST2']].min())




    
    #Surface
    zsurf1 = df.alt[idxb['sst1']:idxe['sst1']].mean()
    #zsurf2 = df.alt[idxb['sst2']:idxe['sst2']].mean()

    wsurf1 = df.w[idxb['sst1']:idxe['sst1']].values - df.w[idxb['sst1']:idxe['sst1']].mean()  
    #wsurf2 = df.w[idxb['sst2']:idxe['sst2']].values - df.w[idxb['sst2']:idxe['sst2']].mean() 
    for var, val in phi.items():
        Budget['sflux1'][var] = np.mean(wsurf1 * (df[var][idxb['sst1']:idxe['sst1']].values - df[var][idxb['sst1']:idxe['sst1']].mean()))
        #Budget['sflux2'][var] = np.mean(wsurf2 * (df[var][idxb['sst2']:idxe['sst2']].values - df[var][idxb['sst2']:idxe['sst2']].mean()))
        
    

    #Now integrate all of the profiles        
    for term, var in profiles.items():
        #For each term in the profiles dict, integrate each one of the variables 
        Budget[term]=vertIntegrate(var,zi)
    #Using the height integrated budget terms, we now compute the residual, thus closing the budget.
    for var in phi.keys():
        Budget['res1'][var]=Budget['tend1'][var] - (Budget['advx1'][var]+Budget['advy1'][var]+Budget['eflux1'][var]+Budget['sflux1'][var])
        #Budget['res2'][var]=Budget['tend2'][var] - (Budget['advx2'][var]+Budget['advy2'][var]+Budget['eflux2'][var]+Budget['sflux2'][var])



    
    
#    #Flux divergence
#    for var, val in phi.items():
#        grad = (Budget['eflux1'][var]-Budget['sflux1'][var])/(ztop1-zsurf1)
#        profiles['fdiv1'][var] = -np.ones(zi.shape)*grad
#        
#        grad = (Budget['eflux2'][var]-Budget['sflux2'][var])/(ztop2-zsurf2)
#        profiles['fdiv2'][var] = -np.ones(zi.shape)*grad    
#    
#    #Residual Terms
#    for var,val in phi.items():
#        profiles['res1'][var] = profiles['tend1'][var] - (profiles['advx1'][var]+profiles['advy1'][var]+profiles['fdiv1'][var])
#        profiles['res2'][var] = profiles['tend2'][var] - (profiles['advx2'][var]+profiles['advy2'][var]+profiles['fdiv2'][var])
#        
#        
#    
#
#    #Now we return all the terms of the box, as well as the 
#    #Vertical profile and integrated budget terms.
    
    
    
    
    return(phi,zi,profiles,Budget)
        
        

def vertIntegrate(phi, z):
    #first determine the method that needs to be used.
    Nz=len(z)
    out={}
    h=z[1]-z[0]
    if Nz%2==0:
        #Nz is even so use simpsons rule
        
        for var,val in phi.items():            
            out[var]=np.sum(phi[var])*h
#            f_ev=phi[var][2:-2:2]
#            f_od=phi[var][1:-1:2]
#            out[var]=1.0/3.0*h*(phi[var][0]+4.0*(np.sum(f_od))+2.0*(np.sum(f_ev))+phi[var][-1])
    else:
        for var,val in phi.items():
            out[var]=h*(.5*phi[var][0]+sum(phi[var][1:-1])+.5*phi[var][-1])        
    return(out)
    
    
    
