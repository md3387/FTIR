# FTIR
Reverse-engineered FTIR post-processing along with integrated absorbance functions\
Written to interface with data output from Thermo-Fisher Scientific Fourier Transform Infrared Spectrometers (FTIR), which typically utilize a proprietary software called "OMNIC." \

## Quick Function Descriptions
\
There are two primary functions:\
1) FTIR External Detector Data Converter\
2) Integrated Absorbance.\
\
There are multiple "External Detector Data Converter" files that are split between transmittance and absorbance, batch and single script processing.  I'm hoping to combine those into a single file with
a simple set of user-defined inputs to switch betweeen those options.  For now, they're similar, but separate functions.

## FTIR Data Processing Concepts
An FTIR is a common instrument in chemistry capable of scanning across nearly the entire infrared spectrum (~500-12000nm) using a michelson interferometer.  The details of it's operation won't be discussed here.  
Suffice to say the operator can choose any range of wavelengths to scan across.  A sample is then placed between the laser emission and a detector, and the transmission of light at those wavelengths is measured.
In a typical FTIR experiment, two scans are taken; a background and a sample. We are usually interested in either the transmittance or absorbance of the sample.  That is found using Beer's law:\
Transmittance = I/I_0\
Where I is the intensity of light transmitted through the sample, and I_0 was the original intensity of light that had been transmitted during the background scan. Put more simply:\
I=Sample\
I_0=Background\
Absorbance is just the -ln of transmittance:\
Absorbance=-ln(Transmittance)=-ln(I/I_0).\
Both transmittance and absorbance are wavelength dependent, so the two equations presented above are evaluated at every wavelength within the scan. Furthermore, the Michelson Interferometer doesn't output
either of these measurements directly.  It outputs what is called an "interferogram," which looks like a noisy signal with a big spike in the middle, called a "centerburst."  The problem with this signal is
that it's a signal intensity on the y-axis, and time on the x-axis. In the end, we are trying to create a plot with Absorbance in arbitrary units on the y-axis, and wavelength (nm or um) or wavenumber (cm-1) 
on the x-axis. This requires a couple more signals from the FTIR, and a fourier transform.  The other two signals needed are a "mirror trigger" and a "Hene Fringe."  The michelson interferometer scans the wavelengths 
using a moving mirror.  When the mirror reaches the extent of it's travel, the scan is over.  So the "mirror trigger" signal is a high/low signal that defines the boundaries of the scan.  A Helium Neon (HeNe) laser 
has a well-defined wavelength, and when passed throught the interferometer, results in a sine wave signal, the peaks of which represent a fairly precise measure of mirror travel distance.  So, the "HeNe Fringe" is a
sine wave that can be used to determine the distance the mirror travels relative to a stationary mirror (termed the optical path difference (OPD) per unit time during the scan.  
Therefore, three signals are required to produce an absorbance plot:  The "Mirror Trigger", the "HeNe Fringe", and
the "Interferogram." \ 
\
The FTIR has a built-in detector and a sample compartment.  We modified ours for sampling outside the sample compartment.  This allows us to do in-situe FITR sampling in some of our experimental setups. This also required an external
detector, which prevented us from using the proprietary OMNIC software.  The external detector signal (the "Interferogram") was then sent to an external oscilloscope. The 'Mirror Trigger" and "HeNe Fringe" had to be sent from the FTIR 
to that external Oscilloscope so that all three signals were on the same time base and recorded by the same instrument.  We then had to write a Matlab code that converts those three signals into either a transmittance or absorbance
plot, like OMNIC does in the standard setup.\
\
So, the postprocessing steps completed in FTIR "External Detector Data Converter" are:\
1) Read in and parse data to process only that between the mirror triggers\
2) Approximate the Optical Path Difference using the HeNe Fringe (Needs work)\
3) Appodization of the Interferogram  - bell-shaped weighting function that multiplies the centerburst by one, and multiplies data farthes from the centerburst by zero.\
4) Zero Filling - fourier transform works best with 2^n data points, so add zeros at the far edges until data length = 2^n.) (Needs to be implemented)\
5) Fast Fourier Transform (fft) -Converts interferogram from time domain into frequency domain, then separates real and imaginary parts\
6) Phase Correction\
\
All IR absorbance is a function of wavenumber, composition, temperature, and pressure.  It's possible, over narrow absorption bands, to remove some of the temperature dependence using integrated band intensity. Over a narrow band,
if one absorption feature gets stronger, another typically becomes weaker so the area under the curve remains constant. However, it's possible that rather than absorbing or transmitting, the sample may also scatter some light so that
wavelengths where there is no absorption may still show reduced transmittance. We call that baseline drift.  It may be necessary in some cases to correct the baseline back to zero if you can accurately determine the drift.\
\
So, the the postprocessing steps completed in FTIR "Integrated Absorbance" are:
1) Baseline Correction
2) Integrate Absorbance
   



