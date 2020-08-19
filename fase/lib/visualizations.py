'''These form the visualization routines for our work in the fase project.'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
#from mpl_toolkits.basemap import Basemap
from os import getcwd
import pandas as pd


BASE_DIR = getcwd()


def budget(budget,phi='qt',tend='FullTend',res='FTres',cloud_filter=False,date="",save=False,show=True):
    fig=plt.figure(figsize=(18,11))
    ax=fig.add_subplot(111)
    ax.grid()
    ax.set_axisbelow(True)
    plt.title("{} Budget".format(phi))
    #The first box is the tend leg only model
    ax.bar(0.5,budget['advx1'][phi],1.0,color='red',label='Advx')
    ax.bar(2.5,budget['advy1'][phi],1.0,color='orange',label='Advy')
    if cloud_filter:
        ax.bar(4.5,budget['f_eflux1'][phi],1.0,color='blue',label='Entr')
        ax.bar(8.5,budget['filter_FTres'][phi],1.0,color='grey',label='Resid')
    else:
        ax.bar(4.5,budget['eflux1'][phi],1.0,color='blue',label='Entr')
        ax.bar(8.5,budget[res][phi],1.0,color='grey',label='Resid')
    ax.bar(6.5,budget['sflux1'][phi],1.0,color='green',label='Surf')
    ax.bar(10.5,budget[tend][phi],1.0,color='yellow',label='Tend')
    ax.legend(loc='center left',bbox_to_anchor=(1,0.5))

    if save:
        figname=BASE_DIR+"../outputs/figures/"+date+"/{}_budget-".format(phi)+date+".png"
        plt.savefig(figname)
    if show:
        plt.show()


def term_comparison(budget,cloud_filter):
    fig=plt.figure(figsize=(18,5))
    ax=fig.add_subplot(111)
    ax.grid()
    ax.set_axisbelow(True)
    plt.title('Moist Static energy entr term comp')
    if cloud_filter:
        ax.bar(0.5,budget['f_eflux1']['se_t1'],1.0,color='red',label='Cp*T')
        ax.bar(2.5,budget['f_eflux1']['se_t2'],1.0,color='orange',label='g*z')
        ax.bar(4.5,budget['f_eflux1']['se_t3'],1.0,color='blue',label='Lv*r')
    else:
        ax.bar(0.5,budget['eflux1']['se_t1'],1.0,color='red',label='Cp*T')
        ax.bar(2.5,budget['eflux1']['se_t2'],1.0,color='orange',label='g*z')
        ax.bar(4.5,budget['eflux1']['se_t3'],1.0,color='blue',label='Lv*r')
    ax.legend(loc='center left',bbox_to_anchor=(1,0.5))
    plt.show()


def entr_scatter(w,phi):
    fig=plt.figure(figsize=(18,5))
    ax=fig.add_subplot(111)
    ax.grid()
    ax.set_axisbelow(True)
    plt.title('Scatter plot of flux-jump terms')
    ax.plot(w,phi,'go')
    ax.axhline(y=0,color='k')
    ax.axvline(x=0,color='k')
    plt.show()


def profiles(date,save=False,show=True):
    ''' This function plots the qt,ql, and theta profiles for the date provided'''
    
    # Load the profile arrays
    pr_fname = BASE_DIR+"flights/"+date+"/profile-"+date+".csv"
    cbn_fname = BASE_DIR+"flights/"+date+"/CABIN_10hz_"+date+".TXT"
    
    pr_index=np.loadtxt(pr_fname,delimiter=",")
    cbn_file=pd.read_csv(cbn_fname)
    
    rh_a=1.225    

    alt1=cbn_file['NovAtel Alt (m)'][int(pr_index[0,0]):int(pr_index[0,1])]
    alt2=cbn_file['NovAtel Alt (m)'][int(pr_index[1,0]):int(pr_index[1,1])]
    
    qv1=cbn_file['MR-h2o (g/Kg)'][int(pr_index[0,0]):int(pr_index[0,1])]
    qv2=cbn_file['MR-h2o (g/Kg)'][int(pr_index[1,0]):int(pr_index[1,1])]

    ql1=pd.Series(cbn_file['HPVM LWC'][int(pr_index[0,0]):int(pr_index[0,1])]/rh_a,index=alt1.index)
    ql2=pd.Series(cbn_file['HPVM LWC'][int(pr_index[1,0]):int(pr_index[1,1])]/rh_a,index=alt2.index)

    qt1=pd.Series(qv1+ql1,index=alt1.index)
    qt2=pd.Series(qv2+ql2,index=alt2.index)

    th1=cbn_file['Theta (K)'][int(pr_index[0,0]):int(pr_index[0,1])]
    th2=cbn_file['Theta (K)'][int(pr_index[1,0]):int(pr_index[1,1])]
    
    # Create the figure
    fig, ax = plt.subplots(nrows=2,ncols=3,figsize=(18,15))
    ax[0][0].plot(qt1,alt1,'bo',label=r'$q_t$')
    ax[0][1].plot(ql1,alt1,'go',label=r'$q_l$')
    ax[0][2].plot(th1,alt1,'ro',label=r'$\theta$')
    
    ax[1][0].plot(qt2,alt2,'bo',label=r'$q_t$')
    ax[1][1].plot(ql2,alt2,'go',label=r'$q_l$')
    ax[1][2].plot(th2,alt2,'ro',label=r'$\theta$')
    
    # Loop through every axis object and set the formatting
    for i in range(3):
        for j in range(2):
            ax[j][i].grid()
            #ax[j][i].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    ax[0][0].set_ylabel("Alt (m)")
    ax[1][0].set_ylabel("Alt (m)")
    # Titles have to be set separately
    ax[0][0].set_title(r"$q_t$ start")
    ax[1][0].set_title(r"$q_t$ end")
    ax[0][0].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[1][0].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[0][1].set_title(r"$q_l$ start")  
    ax[1][1].set_title(r"$q_l$ end")
    ax[0][1].xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax[1][1].xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax[1][2].set_title(r"$\theta$ end")
    ax[0][2].set_title(r"$\theta$ start")
    ax[0][2].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax[1][2].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    if show:
        plt.show()
    if save:
        figname = "/home/mrmisanthropy/Projects/fase/outputs/figures/"+date+"/profiles-"+date+".png"
        plt.savefig(figname)
    
    # Return each profile object.


def combined_figure(bg1,bg2,bg3,bg4,phi='qt',tend='FullTend',res='FTres',show=True,save=False):

    fig, ax = plt.subplots(nrows=2,ncols=2,sharey=True)

    ax[0][0].bar(0.5,bg1['advx1'][phi],1.0,color='red',label='Advx')
    ax[0][0].bar(2.5,bg1['advy1'][phi],1.0,color='orange',label='Advy')
    ax[0][0].bar(4.5,bg1['eflux1'][phi],1.0,color='blue',label='Entr')
    ax[0][0].bar(6.5,bg1['sflux1'][phi],1.0,color='green',label='Surf')
    ax[0][0].bar(8.5,bg1[res][phi],1.0,color='grey',label='Resid')
    ax[0][0].bar(10.5,bg1[tend][phi],1.0,color='yellow',label='Tend')
    ax[0][0].grid()
    ax[0][0].set_axisbelow(True)
    ax[0][0].set_title("July 22")

    ax[0][1].bar(0.5,bg2['advx1'][phi],1.0,color='red',label='Advx')
    ax[0][1].bar(2.5,bg2['advy1'][phi],1.0,color='orange',label='Advy')
    ax[0][1].bar(4.5,bg2['eflux1'][phi],1.0,color='blue',label='Entr')
    ax[0][1].bar(6.5,bg2['sflux1'][phi],1.0,color='green',label='Surf')
    ax[0][1].bar(8.5,bg2[res][phi],1.0,color='grey',label='Resid')
    ax[0][1].bar(10.5,bg2[tend][phi],1.0,color='yellow',label='Tend')
    ax[0][1].grid()
    ax[0][1].set_axisbelow(True)
    ax[0][1].set_title("July 26")

    ax[1][0].bar(0.5,bg3['advx1'][phi],1.0,color='red',label='Advx')
    ax[1][0].bar(2.5,bg3['advy1'][phi],1.0,color='orange',label='Advy')
    ax[1][0].bar(4.5,bg3['eflux1'][phi],1.0,color='blue',label='Entr')
    ax[1][0].bar(6.5,bg3['sflux1'][phi],1.0,color='green',label='Surf')
    ax[1][0].bar(8.5,bg3[res][phi],1.0,color='grey',label='Resid')
    ax[1][0].bar(10.5,bg3[tend][phi],1.0,color='yellow',label='Tend')
    ax[1][0].grid()
    ax[1][0].set_axisbelow(True)
    ax[1][0].set_title("August 8")

    ax[1][1].bar(0.5,bg4['advx1'][phi],1.0,color='red',label='Advx')
    ax[1][1].bar(2.5,bg4['advy1'][phi],1.0,color='orange',label='Advy')
    ax[1][1].bar(4.5,bg4['eflux1'][phi],1.0,color='blue',label='Entr')
    ax[1][1].bar(6.5,bg4['sflux1'][phi],1.0,color='green',label='Surf')
    ax[1][1].bar(8.5,bg4[res][phi],1.0,color='grey',label='Resid')
    ax[1][1].bar(10.5,bg4[tend][phi],1.0,color='yellow',label='Tend')
    ax[1][1].grid()
    ax[1][1].set_axisbelow(True)
    ax[1][1].set_title("August 10")

    # This is the legend, I cant make it fit
    # ax[0][0].legend(loc='center left',bbox_to_anchor=(1,0.5))
    if show:
        plt.show()
    if save:
        figname = BASE_DIR + "../outputs/figures/combined_budget.png"
        plt.savefig(figname)
    plt.close()

def coastline_figures(date,save=False,show=True):
    ''' This function is designed to generate figures of the coastline during flight days

        inputs should just be a string formatted as "DD_MM_YY"'''
    root=BASE_DIR
    wp_file=root + "WP_" + date +".csv"
    st_img = plt.imread(root + "/images/goes-" + date + "-startImage.png")
    ed_img = plt.imread(root + "/images/goes-" + date + "-endImage.png")
    wp_list = np.loadtxt(wp_file, dtype='float',delimiter=',')
    # These are cast as lists to make things easier later
    lat = wp_list[:, 0].tolist()
    lon = wp_list[:, 1]* -1
    lon = lon.tolist()
    # Determine box dimensions
    lat_min = min(lat)
    lon_min = min(lon)

    # Now begin the image
    fig = plt.figure(figsize=(12,6))

    ax1 = plt.subplot2grid((2,2),(0,0),rowspan=2)
    ax1.set_title("Flight plan")
    ax2 = plt.subplot2grid((2,2),(0,1))
    ax2.set_title("Start Conditions")
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax3 = plt.subplot2grid((2,2),(1,1))
    ax3.set_title("End Conditions")
    ax3.set_xticklabels([])
    ax3.set_yticklabels([])
    # First plot the basemap image of our latitude.
    mp = Basemap(width=3E5, height=3E5, projection='lcc',resolution='h', lat_0=lat_min, lon_0=lon_min,ax=ax1)
    mp.drawcoastlines()
    mp.drawlsmask()
    # Convert the waypoints to lat lon coordinates based on the current map
    # lons,lats = mp(lon,lat)

    # This is an optional line that can be used if plotting a box, it sets
    # a coordinate at the end which will close the shape in plotting
    lon.append(lon[0])
    lat.append(lat[0])
    lons,lats = mp(lon,lat)
    mp.plot(lons,lats,'r--')

    # Else this way will plot them as dots
    # mp.scatter(lats,lons,marker='o',color='r')

    # Next plot the two images of the start and end of the flight
    ax2.imshow(st_img)
    ax3.imshow(ed_img)
    if save:
        save_dir = BASE_DIR+"../outputs/figures/" + date + "/coastline-"+date+".png"
        plt.savefig(save_dir)
    if show:
        plt.show()
    plt.close()

