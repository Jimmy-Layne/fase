# Fog And Stratocumulus Evolution experiment analysis tools.
this package contains the analysis for the FASE experiment. It is designed to determine the most influential factors in stratocumulus lifetime over a volume of space using aircraft observation. Data files are from flights of the CIRPAS twin-otter aircraft during July and August of 2016.

### Installation
In order to install first run
```
    git clone https://github.com/Jimmy-Layne/fase.git
```

Then install all required packages from the requirements file using pip with:
pip install -r requirements.txt

### Usage
Analysis from the FASE project and all associated figures may be obtained by running the thesis_figures.py script. Determination of cloud top was accomplished with the ct_determination.py script.  The analysis tools may be used directly from python as well:
```python
    import fase
    phi, profiles, budget = fase.VolumeBudget('date')
```
Where date is a MM_DD_YY string which corresponds to a date in the fase/flights directory, for which there are valid data.
The VolumeBudget function will return conserved quantities, interpolated profiles, and budget terms in that order
