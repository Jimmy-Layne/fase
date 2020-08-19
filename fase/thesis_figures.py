"""
This script runs our fase analysis for all dates included in the paper

This forms the master generator for all of our analysis figures. It runs for each of the
days in question, and saves those figures silently to the outputs directory.
"""

import lib.Budget_V3 as bm
import lib.visualizations as vis
import lib.flight_info as fd

##################################################################
#           Functions
################################################################

date = [
    fd.Date[0],
    fd.Date[1],
]
# Now Obtain budget data for each day.
phi_1,zi_1,prf_1,bg_1=bm.volumeBudget(date[0])
phi_2,zi_2,prf_2,bg_2=bm.volumeBudget(date[1])
#First parse out the arguments and make the file names


# Next generate all the figures

# Coastline figures (Can't figure out how to suppress show with basemap)
# basemap is broken curently
# vis.coastline_figures(date[0],save=True, show=False)
# vis.coastline_figures(date[1],save=True, show=False)
# vis.coastline_figures(date[2],save=True, show=False)

# Mass budget Figures
vis.budget(bg_1,date=date[0], save=True, show=False)
vis.budget(bg_2,date=date[1], save=True, show=False)

# Moist Static Energy budget Figures
vis.budget(bg_1,date=date[0], phi='se',save=True, show=False)
vis.budget(bg_2,date=date[1], phi='se',save=True, show=False)

# Potential Temperature budget figures
vis.budget(bg_1,date=date[0], phi='theta',save=True, show=False)
vis.budget(bg_2,date=date[1], phi='theta',save=True, show=False)

# Profiles
vis.profiles(date=date[0], save=True, show=False)
vis.profiles(date=date[1], save=True, show=False)


