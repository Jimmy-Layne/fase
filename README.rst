*****
Fog And Stratocumulus Evolution experiment analysis tools.
*****
this package contains the analysis for the FASE experiment. It is designed to determine the most influential factors in stratocumulus lifetime over a finite volume of space. Data files are from flights of the CIRPAS twin-otter aircraft during July and August of 2016.

Installation
#########
In order to install first run:
::
	git clone <BIG URL HERE>
Then navigate to the correct package and run <make?> to install all required packages and setup.
Usage
#########
Tools from the fase analysis package can be used by calling:
::
	import fase
	fase.VolumeBudget('date')

Where date is a MM_DD_YY string which corresponds to a date in the fase/flights directory, for which there is valid analysis. The analysis from the FASE project can be run by executing the script at bin/FASE.py which will generate all figures associated with the analysis.