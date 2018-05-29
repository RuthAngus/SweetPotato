# SweetPotato

Inferring trends in exoplanet populations over time.

Steps:
=======

1. Check what the influence of metallicity is on planet occurrence rate.
2. Is there a trend between vertical action and temperature?
3. Is there a trend between vertical action and any other stellar parameter?
4. Does the cut on parallax signal to noise bias you somehow?
5.

Code:
=====


teff_occurrence_rate.ipynb: An implementation of planet occurrence rate as a
function of temperature.

detection_efficiency_functions.py: code that implements Chris Burke's
completeness model.

data_wrangling.py: Some functions for downloading stellar catalogues.

plotting.py: some plotting functions used to make marginalised 1 and 2D
histograms. No longer used.

teff_bv.py: converts teff to b-v and vice versa. Not currently used.

exopops/action_occurrence.ipynb: Exoplanet occurrence rate as a function of
vertical action.
