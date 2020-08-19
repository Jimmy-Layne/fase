''' This module contains all of the file utilities for the fase analysis'''


import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd
import math
from lib.structures import get_rename, get_phi
from lib.computations import create_histogram

def distance_on_earth(lat1, long1, lat2, long2):
    '''This function converts lat and lon coordinates to distances on earth''' 
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
     
    cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) + math.cos(phi1)*math.cos(phi2))
    arc = math.acos(cos)
     
    # Remember to multiply arc by the radius of the earth
    # in your favorite set of units to get length.
    arc *= 6.373e6
    return arc


def flight_plot(fname,plot=False):
    '''This function plots the lat,lon,alt,and heading of the flight, from which we derive our index files'''

    df = pd.read_csv(fname,sep = ",",index_col = "Rec ID.", na_values=[-9999])
    df.rename(columns=get_rename(), inplace=True)


    lat = df.Lat.values
    lon = df.Long.values
    alt = df.alt.values
    hed = df.heading.values

    if plot:
        #This plots all the different plots on the same figure
        f, ax = plt.subplots(4,sharex=True)
        ax[0].plot(df.index,lat)
        ax[1].plot(df.index,lon)
        ax[2].plot(df.index,alt)
        ax[3].plot(df.index,hed)



        # optionally we can have them on separate plots
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
    return(lat,lon,alt,hed)

# >>>> Begin Cloudtop Determination functions
# Visualization of cloudtop guesses
def ct_visualizations(leg_data, grad, guess, leg):
    bar_loc = [x[3] for x in grad]
    gradients = [x[1] for x in grad]
    diffs = [x[2] for x in grad]
    fig, ax = plt.subplots(2, 2)
    # Potential energy profile
    ax[0][0].plot(leg_data['theta'], leg_data['alt'], 'go')
    ax[0][0].axhline(guess[3], color='red', linestyle="--")
    ax[0][0].axhline(guess[4], color='red', linestyle='--')
    plt.ylabel("Alt (m)")
    plt.xlabel(r"$\theta$ (K)")
    # Liquid water profile
    ax[0][1].plot(leg_data['ql'], leg_data['alt'], 'bo')
    ax[0][1].axhline(guess[3], color='red', linestyle="--")
    ax[0][1].axhline(guess[4], color='red', linestyle='--')
    plt.ylabel("Alt (m)")
    plt.xlabel("ql")
    # Bar plot of gradients
    ax[1][0].bar(bar_loc, gradients)
    plt.ylabel(r"$\nabla \theta$")
    # Bar plot of max differences
    ax[1][1].bar(bar_loc, diffs)
    plt.ylabel(r"$\Delta \theta$")

    plt.suptitle(r"visualization of leg: {}".format(leg))
    plt.show()

# Routine for making guesses at the correct bin
def make_guesses(data):
    # the three greatest metric terms will form the top three guesses, so we'll start by sorting the list
    while True:
        user_guess = input("Please enter the integer number that you think is closest to the true value")
        user_guess = int(user_guess)
        for i in range(len(data)):
            low = data[i][3]
            high = data[i][4]
            if low <= user_guess < high:
                return data[i]
        print("I couldn't find the bin containing {}".format(user_guess))

# Main Function governing determination of cloudtop
def cloud_top_determination(leg, name):
    '''A diagnostic function to determine cloudtop based on gradients in phi

    This function determines the cloud top buffer-region of each section of the entrainment leg as well as
    the flux jump in qt over the gap
    '''
    # Create a histogram object wrt altitude
    alt = leg['alt']
    dim = {'dn': 10}
    z_hist = create_histogram(dim, alt)

    # loop over the bins to determine gradients
    gradient = []
    for i in range(z_hist['dim']['N']):
        ind_list = z_hist['index'][i]
        events = leg[leg.index.isin(ind_list)]
        tmp_grad = np.average(np.abs(np.gradient(events['theta'])))
        tmp_diff = np.max(np.abs(events['theta'].diff()))
        # Take the sum of gradient and difference. to form our metric for decision making
        # Note, here we will scale the gradient term by 10 in order to make the two comparable
        tmp_metric = tmp_grad * 10 + tmp_diff
        grad = (tmp_metric, tmp_grad, tmp_diff, z_hist['edges'][i], z_hist['edges'][i+1], i)
        gradient.append(grad)
    # Now Make a guess at the place where the cloudtop most likely falls, based on the steepest gradients.
    # Sorting a list of tuples works based on the first element in the tuple, which in this case is our metric
    # Thus we take the top three of the sorted list.
    # Since our analysis has shown that cloud-top is better represented as a buffer layer between the STBL and the
    # Free troposphere, cloud top is analyzed over the top two bins, the last cloudy air bin and the first clear bin
    # here buffer_layer[0] is the
    gradient.sort()
    guesses = gradient[-3:]
    buffer_layer = []
    ct_found = False
    while not ct_found:

        for g in guesses:
            ct_visualizations(leg, gradient, g, name)
            # Now confirm the guess
            found = input("Does that guess work?")
            if found == 'y' or found == 'Y':
                indx = g[5]
                ct_height = (g[3], g[4], indx)
                ct_found = True
                bin = input("Is this the first clear bin(1), or the last cloud bin(2)")

                if bin == "1":
                    buffer_layer.append((z_hist['edges'][indx-1], z_hist['edges'][indx], indx -1))
                    buffer_layer.append(ct_height)
                elif bin == "2":
                    buffer_layer.append(ct_height)
                    buffer_layer.append((z_hist['edges'][indx+1], z_hist['edges'][indx + 2], indx+1))
                break

        if not ct_found:
            while True:
                guess = make_guesses(gradient)
                ct_visualizations(leg, gradient, guess, name)

                found = input("Does that guess work?")
                if found == 'y' or found == 'Y':
                    indx = guess[5]
                    ct_height = (guess[3], guess[4], guess[5])
                    ct_found = True
                    bin = input("Is this the first clear bin(1), or the last cloud bin(2)")
                    if bin == "1":
                        buffer_layer.append((z_hist['edges'][indx-1], z_hist['edges'][indx], indx))

                    elif bin == "2":
                        buffer_layer.append(ct_height)
                        buffer_layer.append((z_hist['edges'][indx + 1], z_hist['edges'][indx], indx))

                    break

        # now that we've determined the altitude of cloud top, we need to determine Delta phi over cloud top.
    phi_keys = get_phi()
    phi_keys = [x for x in phi_keys.keys()]
    ev_1 = leg[leg.index.isin(z_hist['index'][buffer_layer[0][2]])]
    ev_1 = ev_1[phi_keys]
    ev_2 = leg[leg.index.isin(z_hist['index'][buffer_layer[1][2]])]
    ev_2 = ev_2[phi_keys]
    Dphi={}

    for var in ev_1.keys():
        phi_1=np.average(ev_1[var])
        phi_2=np.average(ev_2[var])
        Dphi[var] = phi_2 - phi_1
    # phi_1 = np.average(ev_1['qt'])
    # phi_2 = np.average(ev_2['qt'])
    # Dphi = phi_2 - phi_1

    return Dphi, buffer_layer
