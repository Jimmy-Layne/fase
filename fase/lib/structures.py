''' This file contains all of the dictionary structures that are used in the computation of the Budget

    #   -phi holds the linearly interpolated data for each of our conserved quantities to analyze
    #   -profiles holds the interpolated height profiles for each flight leg
    #   -Budget holds the completed budget with both integrated profie terms and computed quantities
    #   -rename is the list of name changes that goes into the initial loading of the Cabin file
    #   -constants is a list of physical and mathematical constants used in the analysis

'''

# Rename
def get_rename():
    rename = {
            'Wind Speed (m/s)': 'ws',
            'Wind dir(Deg)': 'wdir',
            'Palt (m)':'alt',
            'Heading (deg)':'heading',
            'MR-h2o (g/Kg)':'qv',
            'Vert. Wind (m/s)':'w',
            'L5200 CO2 (V)':'co2',
            'Thetae':'thetae',
            'Theta (K)': 'theta',
            'SPHUM (g/Kg)': 'r',
            'Tamb (C)': 'T',
            'UTC time(HH.)': 'time'
            }

    return rename

# constants
def get_constants():
    constants = {
        'rho_air': 1.225,
        'Cp': 1.003,
        'Lv': 2260.,
        'g': 9.8,
        'kmpdeglat': 111.,
        'kmpdeglon': 90.,
        'sample_freq': 10.,
        'hist_width': 2.0
        }

    return constants

# Phi
def get_phi():
    phi = {
            'u':{},
            'v':{},
            'w':{},
            'qv':{},
            'ql':{},
            'qt':{},
            'theta':{},
            'co2':{},
            'se':{},
            'se_t1':{},
            'se_t2':{},
            'se_t3':{}
            }

    return phi

# Profiles
def get_profiles():
    profiles = {
            'tend1':{},
            'tend2':{},
            'advx1':{},
            'advx2':{},
            'advy1':{},
            'advy2':{},
            'AMtend':{},
            'MBtend':{},
            'BMtend':{},
            'MCtend':{},
            'CMtend':{},
            'MDtend':{},
            'DMtend':{},
            'MAtend':{}
            }

    return profiles

# Budget
def get_budget():
    Budget = {
             'tend1':{},
             'tend2':{},
             'advx1':{},
             'advx2':{},
             'advy1':{},
             'advy2':{},
             'res1':{},
             'res2':{},
             'eflux1':{},
             'eflux2':{},
             'sflux1':{},
             'sflux2':{},
             'AMtend':{},
             'MBtend':{},
             'BMtend':{},
             'MCtend':{},
             'CMtend':{},
             'MDtend':{},
             'DMtend':{},
             'MAtend':{},
             'Advtend':{},
             'Advres':{},
             'FullTend':{},
             'FTres':{},
             'f_eflux1':{},
             'filter_FTres':{}
             }
    
    return Budget

