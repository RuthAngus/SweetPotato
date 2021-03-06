# SweetPotato

Inferring trends in exoplanet populations over time.

Steps:
=======

1. [x] Check what the influence of metallicity is on planet occurrence rate.
    -- (see figure below.)
2. [x] Is there a trend between vertical action and temperature?
    -- [x] yes, so we need to plot occurrence rate as a function of this too.
    -- [x] looks like only the hot stars show this effect, is that because they
       have a broader range of vertical actions? There just maybe are not
       enough cool stars.
3. Is there a trend between vertical action and any other stellar parameter?
4. Does the cut on parallax signal to noise bias you somehow?
5. Is it ok to use vertical actions and not vertical action dispersion.
6. Does binarity influence this?
7. Flags epsi and sepsi - astrometric excess noise and astrometric excess
   noise significance.
8. Triple check the Gaia data - are there some systematic effects that might
   be influencing this.

Gijs Mulders' Metallicity -- Occurrence rate plot looks like this -->

![alt text](https://github.com/ruthangus/SweetPotato/blob/master/Mulders.png)

So cutting on a Metallicity of -.1 to +.1 should remove most of this trend.


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

exopops/action_occurrence.ipynb: Exoplanet occurrence rate as a function of vertical action.

code/Checking_actions.ipynb: comparing my actions with Adrian's.
