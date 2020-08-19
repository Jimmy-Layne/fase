# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 19:53:12 2016

@author: mrmisanthropy
"""
#This is the mass budgeting stuff we wrote
import MOD812_BudgetModule as bm
import argparse
import numpy as np
import matplotlib.pyplot as plt

import pdb

##################################################################
#           Functions
################################################################
def figureLabels(var,Nz):
    if var=='thetae':
        title=r"Integrated budget terms of $\theta_e$ with $N_z$={}".format(Nz)
        label=r"$\theta_e$ (K)"
        
    elif var=='u':
        title=r"Integrated budget terms of $u$ with $N_z$={}".format(Nz)
        label=r"$u$ ($m s^{-1}$)"
        
    elif var=='v':
        title=r"Integrated budget terms of $v$ with $N_z$={}".format(Nz)
        label=r"$v$ ($m s^{-1}$)"
        
    elif var=='w':
        title=r"Integrated budget terms of $w$ with $N_z$={}".format(Nz)
        label=r"$w$ ($m s^{-1}$)"
        
    elif var=='qt':
        title=r"Estimated Water budget terms"
        label=r"$q_t$ Budget Contribution ($g kg^{-1}s^{-1}$)"
        
    elif var=='ql':
        title=r"Integrated budget terms of $q_l$ with $N_z$={}".format(Nz)
        label=r"$q_l$ ($g kg^{-1}$)"
        
    elif var=='qv':
        title=r"Integrated budget terms of $q_v$ with $N_z$={}".format(Nz)
        label=r"$q_v$ ($g kg^{-1}$)"
        
    elif var=='co2':
        title=r"Integrated budget terms of $CO_2$ with $N_z$={}".format(Nz)
        label=r"$CO_2$ ($g kg^{-1}$)"
    return(title,label)



#First parse out the arguments and make the file names
parser=argparse.ArgumentParser()
parser.add_argument('dFlag',metavar='Cf',type=str,nargs='+',help="Date of files to be analyzed, formatted like: 'month_day_2016'")
#parser.add_argument('Cfile',metavar='Cf',type=str,nargs='+',help="Cabin File Name")
#parser.add_argument('Ifile',metavar='If',type=str,nargs='+',help="Name of Index file")
#parser.add_argument('Wfile',metavar='Wf',type=str,nargs='+',help="Name of Waypoint file")

parser.add_argument('vertRes',metavar='vR',type=int,nargs='+',help="Number of points to include in the vertical profile")
parser.add_argument('Ctop',metavar='ct',type=int,nargs='+',help='Average cloud top level')
args=parser.parse_args()

Date=args.dFlag[0]
Cfile="/home/mrmisanthropy/Projects/fase/flights/"+Date+"/CABIN_10hz_"+Date+".TXT"
Ifile="/home/mrmisanthropy/Projects/fase/flights/"+Date+"/"+Date+"_Index.csv"
Wfile="/home/mrmisanthropy/Projects/fase/flights/"+Date+"/WP_"+Date+".csv"


Nz=args.vertRes[0]
CT=args.Ctop[0]
phi,zi,profiles,budget=bm.massBudget(Cfile,Ifile,Wfile,Nz,CT)



#Now Create visual Outputs, each plot should cover one of the conserved variables and show the change in each of the terms with changing z
#Variables: co2,qt,qv,ql,thetae,u,v,w

#First create the save directory

savedir="/home/mrmisanthropy/Projects/fase/flights/"+Date+"/Outputs/z{}".format(Nz)


###########################################
#   This is some diagnostic stuff
LHS1=['advx1','advy1','sflux1','eflux1']
print('Budget results')
print('------------------')
print('box -1:')
print('------------------')
for var in phi.keys():
    acc=0.0   
    acc=budget['tend1'][var]-(budget['advx1'][var]+budget['advy1'][var]+budget['eflux1'][var]+budget['sflux1'][var]+budget['res1'][var])
    print("{}-Diff: {}".format(var,np.abs(acc)))



    
    
########################################################
###################
#Profiles plot
##################
for var,val in phi.items():
        figname=savedir+"/Profiles/"+var+"_z{}_".format(Nz)+Date+".png"
        f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(11,8))
        ax1.plot(np.zeros(zi.shape),zi,linewidth=1,color='k',linestyle=':')
        l1 = ax1.plot(profiles['tend1'][var],zi,linewidth=2,label='Tend.',color='k')
        ax1.plot(profiles['advx1'][var],zi,linewidth=2,label='U-Adv.',color='r')
        ax1.plot(profiles['advy1'][var],zi,linewidth=2,label='V-Adv.',color='b')
        ax1.grid()
        plt.ylim(0,300)
        ax1.set_title('Box 1')
        
#        ax2.plot(np.zeros(zi.shape),zi,linewidth=1,color='k',linestyle=':')
#        ax2.plot(profiles['tend2'][var],zi,linewidth=2,label='Tend.',color='k')
#        ax2.plot(profiles['advx2'][var],zi,linewidth=2,label='U-Adv.',color='r')
#        ax2.plot(profiles['advy2'][var],zi,linewidth=2,label='V-Adv.',color='b')
#        ax2.grid()
#        ax2.set_title('Box 2')
        ax1.locator_params(axis='x',nbins=2)
        #ax2.locator_params(axis='x',nbins=2)
    
        ax1.legend(loc='lower right',ncol=4)
        plt.suptitle("Vertical profiles of " + var + r" with $N_z$={}".format(Nz),fontsize=24)
        plt.savefig(figname)
        plt.close()

    ####################
    #Budget terms plot
    ####################
        
        figname=savedir+"/Budget_Terms/"+var+"_z{}_".format(Nz)+Date+".png"
        f, ax1 = plt.subplots(1, 1, sharey=True, figsize=(11,8))
        #Plot the first Box terms
        ax1.bar(0.5,budget['advx1'][var],1.0,color='red',label='Advx')
        ax1.bar(2.5,budget['advy1'][var],1.0,color='orange',label='Advy')
        ax1.bar(4.5,budget['eflux1'][var],1.0,color='blue',label='Entr')
        ax1.bar(6.5,budget['sflux1'][var],1.0,color='green',label='Surf')
        ax1.bar(8.5,budget['res1'][var],1.0,color='grey',label='Resid')
        ax1.bar(10.5,budget['tend1'][var],1.0,color='yellow',label='Tend')        
        box=ax1.get_position()
        ax1.set_position([box.x0,box.y0,box.width*0.8,box.height])
        ax1.legend(loc='center left',bbox_to_anchor=(1,0.5))
        ax1.grid()
        ax1.set_axisbelow(True)
#        #Plot the second box terms
#        ax2.bar(0.5,budget['advx2'][var],1.0,color='red',label='Advx')
#        ax2.bar(2.5,budget['advy2'][var],1.0,color='orange',label='Advy')
#        ax2.bar(4.5,budget['eflux2'][var],1.0,color='blue',label='Entr')
#        ax2.bar(6.5,budget['sflux2'][var],1.0,color='green',label='Surf')
#        ax2.bar(8.5,budget['res2'][var],1.0,color='grey',label='Resid')
#        ax2.bar(10.5,budget['tend2'][var],1.0,color='yellow',label='Tend')        
        #Figure Formatting
        title,label=figureLabels(var,Nz)
        plt.suptitle(title,fontsize=24)
        plt.ylabel(label,fontsize=16)
        plt.savefig(figname)
        plt.close()
        