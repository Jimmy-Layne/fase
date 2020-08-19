'near''This file holds all of the specific information for the flights including filepaths and cloudtop height'''
from pathlib import Path, PosixPath

cwd = Path.cwd()
flight_722={}
flight_726={}
flight_804={}
flight_810={}
Date=['07_22_16','07_26_16','08_04_16','08_10_16']


# flight on 07-22-2016 near Half Moon Bay

flight_722['Cfile'] = cwd/"flights"/Date[0]/PosixPath("CABIN_10hz_"+Date[0]+".TXT")
flight_722['Ifile'] = cwd/"flights"/Date[0]/PosixPath(Date[0]+"_Index.csv")
flight_722['Wfile'] = cwd/"flights"/Date[0]/PosixPath("WP_"+Date[0]+".csv")
flight_722['EIfile'] = cwd/"flights"/Date[0]/PosixPath("entrST-"+Date[0]+".csv")
flight_722['Ctop'] = 250
flight_722['Nz'] = 10

# flight on 07-26-2016 near Big Sur
flight_726['Cfile'] = cwd/"flights"/Date[1]/PosixPath("CABIN_10hz_"+Date[1]+".TXT")
flight_726['Ifile'] = cwd/"flights"/Date[1]/PosixPath(Date[1]+"_Index.csv")
flight_726['Wfile'] = cwd/"flights"/Date[1]/PosixPath("WP_"+Date[1]+".csv")
flight_726['EIfile'] = cwd/"flights"/Date[1]/PosixPath("entrST-"+Date[1]+".csv")
flight_726['Ctop'] = 250
flight_726['Nz'] = 10

# flight on 08-04-2016 near Marina
flight_726['Cfile'] = cwd/"flights"/Date[2]/PosixPath("CABIN_10hz_"+Date[2]+".TXT")
flight_726['Ifile'] = cwd/"flights"/Date[2]/PosixPath(Date[2]+"_Index.csv")
flight_726['Wfile'] = cwd/"flights"/Date[2]/PosixPath("WP_"+Date[2]+".csv")
flight_726['EIfile'] = cwd/"flights"/Date[2]/PosixPath("entrST-"+Date[2]+".csv")
flight_804['Ctop'] = 450
flight_804['Nz'] = 10

# flight on 08-10-2016 near Half Moon Bay
# Deprecated
flight_810['Cfile'] = cwd/Date[3]/"CABIN_10hz_"/Date[3]/".TXT"
flight_810['Ifile'] = cwd/Date[3]/""/Date[3]/"_Index.csv"
flight_810['Wfile'] = cwd/Date[3]/"WP_"/Date[3]/".csv"
flight_810['Ctop'] = 570
flight_810['Nz'] = 10

def get_flight(dt):
    if dt == Date[0]:
        return flight_722
    elif dt == Date[1]:
        return flight_726
    elif dt == Date[2]:
        return flight_804
    elif dt == Date[3]:
        return flight_810
    else:
        return {}
