'''
This script runs our cloud top determination routine for the entrainment object

Starts the ct determination routine for the entrainment flux computation. For the dates provided
it saves the user-determined cloud top regions in a file /flights/<date>/CloudtopHeights_<date>.json
this file must be present for the entrainment computation to be completed.
'''
import lib.entrainment as ent
Dates = ['07_26_16']

for dt in Dates:
	ent_obj = ent.EntrAnalysis(date=dt, ctfind=True)
