"""This impliments the class defining the Ent_analysis object"""
import numpy as np
import json
from pandas import DataFrame, read_csv, read_json, Series
from copy import deepcopy
from scipy.signal import savgol_filter
from scipy.integrate import simps, trapz
import matplotlib.pyplot as plt
# Local modules
from lib.flight_info import get_flight as getFlight
from lib.utilities import cloud_top_determination
from lib.computations import get_index, init_cbfile, vertIntegrate, create_histogram, compute_constants
import lib.structures as stct


class EntrAnalysis:

    help = "A class to hold the entrainment leg mapping and fitting, requires a date to analyze."

    def __init__(self, date, df=DataFrame(), const={}, ctfind=False):

        # ~~Initialize attributes
        self.date = date

        # Get all of the filepaths for the flight
        self.fl_info = getFlight(date)

        # Check to see if df or const need to be loaded
        if df.empty or not const:
            self.en_cbfile, self.const = self.cb_data_init()
        else:
            self.en_cbfile = df
            self.const = const

        # Begining and ending indexes from the entrainment index file
        self.ein_b, self.ein_e = get_index(self.fl_info['EIfile'],self.en_cbfile.index)
        # Obtain raw and smoothed phi from our indexes
        self.phi_raw, self.phi_smooth = self.get_conserved_quantities()
        # Check to see if CT determination needs to be run before entrainment is computed
        if ctfind:
            self.find_cloudtop()    
        # compute entrainment flux
        self.avg_flux, self.ent_flux, self.ent_data=self.compute_entrainment(ctfind=False)



    # METHODS
    def cb_data_init(self):
        constants, lat_lon = compute_constants(self.fl_info['Wfile'])
        # First Initialize the index and data files for the cabin file
        tmp = init_cbfile(self.fl_info['Cfile'], constants)
        return tmp, constants

    def get_conserved_quantities(self):
        '''This function determines phi along our entrainment leg, it returns raw and smoothed structures

        Smoothing is accomplished using savitzky-golay filtering. No interpolation is done in this step, so phi in the
        entrainment object differs from phi in the budgeting calculations, even though the structure is the same.'''
        phi_raw = stct.get_phi()
        phi_smooth = stct.get_phi()
        for key, bg in self.ein_b.items():
            ed = self.ein_e[key]
            for var, val in phi_raw.items():
                zi = self.en_cbfile['alt'][bg:ed].values
                cb_slice = self.en_cbfile[var][bg:ed].values
                phi_raw[var][key]=np.vstack((zi, cb_slice))
                phi_smooth[var][key] = savgol_filter(self.en_cbfile[var][bg:ed].values, 51, 3)
        return phi_raw, phi_smooth


    def compute_entrainment(self, ctfind = False):
        ''' This function computes the histogram and average w for use in our computation of the w'phi' calculation

        Here we compute the histogram that will allow us to average w over the whole bin for use in the RMS calculation
        bin frequency is set by h_freq usually 2Hz which corresponds to around _______ vertical change. '''
        
        # These three structures will be our outputs
        # This one is for full diagnostics (Note Eflux includes only qt here)
        ent_data = DataFrame(columns=['leg', 'eflux', 'dphi', 'wprime'])
        # This one holds the ent flux for each leg
        ent_flux = {p:[] for p in self.phi_raw.keys()}
        # This one is for the average flux, used in budgeting
        avg_flux = {p: None for p in self.phi_raw.keys()}

        # hist. frequency
        dt = self.const['hist_width']
        # hist. dimensions
        # averaged v_velocity data structure
        flight_legs = {leg: None for leg in self.ein_b.keys()}
        # this structure will hold the (w', zbar) mapping for each leg of the flight. it's primary use is for visuals
        wp_mapping = {}
        # Here we'll determine the cloud-top buffer zone and qt flux jump for each of the flight legs
        # Load the 
        ctop_fname = "/home/mrmisanthropy/Projects/fase/fase/flights/{}/CloudtopHeights_{}.json".format(self.date, self.date)
        try:
            cloud_top = read_json(ctop_fname, orient='index')
        except FileNotFoundError:
            print("Ct file not found, run ct determination?")
            raise FileNotFoundError

        # Here's a quick way check of our ct determination
        # for leg in flight_legs.keys():
        #     leg_data = self.en_cbfile[self.ein_b[leg]:self.ein_e[leg]]
        #     ct_entry = cloud_top.loc[leg, :]
        #     self.ct_check(leg_data, ct_entry, leg, save=True, show=False)

        for leg in flight_legs.keys():
            # Take a slice of the data for the leg
            leg_data = self.en_cbfile[self.ein_b[leg]:self.ein_e[leg]]
            
            # Create histogram of arrival times for averaging
            ts_data = leg_data['time']
            dim = {'dn': dt}
            time_series = create_histogram(dim, ts_data)
            leg_map = []
            for i in range(time_series['dim']['N']):
                ind_list = time_series['index'][i]
                events = leg_data[leg_data.index.isin(ind_list)]
                wbar = np.average(events['w'])
                zbar = np.average(events['alt'])
                w_prime_sq = np.average((events['w'] - wbar)**2)
                wprime = np.sqrt(w_prime_sq)
                leg_map.append((wprime, zbar))
            leg_map = DataFrame(leg_map, columns=["w_prime", "alt"])
            wp_mapping[leg] = leg_map
            # Now using our estimation of cloudtop, compute wprime at the altitude of the center of the buffer layer
            leg_ctop = cloud_top.loc[leg, :]
            z_top = leg_ctop['cld_bin'][1]
            wp = leg_map.iloc[(leg_map['alt'] - z_top).abs().argsort()[0]]
            
            # This bit computes only for qt, its used for diagnostics
            leg_eflux = {}
            leg_eflux['leg'] = leg
            leg_eflux['eflux'] = wp['w_prime'] * leg_ctop['Dphi']['qt']
            leg_eflux['dphi'] = leg_ctop['Dphi']['qt']
            leg_eflux['wprime'] = wp['w_prime']
            ent_data = ent_data.append(Series(leg_eflux), ignore_index=True)

            # Here we append the ent flux value for each phi variable
            for ph in self.phi_raw.keys():
                ent_flux[ph].append(wp['w_prime'] * leg_ctop['Dphi'][ph])
        
        # Here we average each ent flux        
        for ph in self.phi_raw.keys():
            avg_flux[ph] = np.average(ent_flux[ph])
        
        return avg_flux, ent_flux, ent_data
    

    def find_cloudtop(self):
        # Set the name for the cloudtop file
        ctop_fname = "/home/mrmisanthropy/Projects/fase/fase/flights/{}/CloudtopHeights_{}.json".format(self.date, self.date)
        flight_legs = {leg: None for leg in self.ein_b.keys()}
        out_data={}
        for leg in flight_legs:
            # Cloudtop file headers look like:
            # leg, Delta qt,Delta se cloudbin_low, cloudbin_high, clearbin_low, clearbin_high
            leg_data = self.en_cbfile[self.ein_b[leg]:self.ein_e[leg]]
            Dphi, buffer_layer = cloud_top_determination(leg_data, leg)
            out_data[leg] = {
            'Dphi': Dphi,
            'cld_bin':(buffer_layer[0][0],buffer_layer[0][1]),
            'clr_bin':(buffer_layer[1][0],buffer_layer[1][1])
            }

        with open(ctop_fname, 'w') as f:
            json.dump(out_data,f)

    # VISUALIZATIONS
    def plot_raw_vs_smooth(self):
        # This function is designed to compare smoothed profiles to raw ones
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for leg in self.phi_raw['qt'].keys():
            ax.plot(self.phi_raw['qt'][leg][1], self.phi_raw['qt'][leg][0])
            ax.plot(self.phi_smooth['qt'][leg][:], self.phi_raw['qt'][leg][0,:], '--')
            plt.xlabel("Total Water")
            plt.ylabel("Altitude")

    def plot_gradient(self):
        # This function is designed to plot gradient profiles
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for leg in self.grad_profiles['qt'].keys():
            ax.plot(self.grad_profiles['qt'][leg][:], self.phi_raw['qt'][leg][0])
            plt.xlabel("grad qt")
            plt.ylabel("Altitude")


    def flight_leg_map_plot(self, fl_cbfile,lvl_ind):
        """This function charts the lat lon of the level and st legs for testing"""
        st_ind = (min(self.ein_b.values()),max(self.ein_e.values()))
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(fl_cbfile['Lat'][lvl_ind[0]:lvl_ind[1]], fl_cbfile['Long'][lvl_ind[0]:lvl_ind[1]],'r--', label= "Level Leg")
        ax.plot(self.en_cbfile['Lat'][st_ind[0]:st_ind[1]], self.en_cbfile['Long'][st_ind[0]:st_ind[1]], 'b--', label= "Sawtooth leg")
        plt.xlabel('latitude')
        plt.ylabel('longitude')
        plt.title("Lat vs Lon plot for st vs lvl legs")
        plt.legend()
        plt.show()


