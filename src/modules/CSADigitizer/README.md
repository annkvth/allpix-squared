# CSADigitizer
**Maintainer**: 
**Status**: 
**Input**: PixelCharge  
**Output**: PixelHit  

### Description
Digitization module which translates the collected charges into a digitized signal.


With the `output_plots` parameter activated, the module produces histograms of the charge distribution at the different stages of the simulation, i.e. before processing, with electronics noise, after threshold selection, and with ADC smearing applied.
A 2D-histogram of the actual pixel charge in electrons and the converted charge in ADC units is provided if ADC simulation is enabled by setting `adc_resolution` to a value different from zero.
In addition, the distribution of the actually applied threshold is provided as histogram.


### Parameters
* `gain_smearing` : Standard deviation of the Gaussian uncertainty in the gain factor. Defaults to 0.
* `threshold` : Threshold for considering the collected charge as a hit. Defaults to 600 electrons.
* `output_plots` : Enables output histograms to be be generated from the data in every step (slows down simulation considerably). Disabled by default.


### Usage
The default configuration is equal to the following:

```ini
[CSADigitizer]
electronics_noise = 110e
threshold = 600e
threshold_smearing = 30e
adc_smearing = 300e
```
