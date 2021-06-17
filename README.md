# DHPSFU
Script to find double-dots and extract molecular coordinates from DHPSF microscopy data, pre-processed by extracting single localisations.

_Authors: Aleks Ponjavic, University of Leeds; Aleksandra Jartseva, University of Cambridge_

To run this, you need both the main script DHPSFUAngle.m and the auxiliary function calibAngle.m in your working directory. Open DHPSFUAngle.m, edit the necessary parameters, and run it.

## Inputs - files:
- Data for analysis in the form of a text file containing a list of single localisations. To enhance batch processing, the code loops through all files with the correct extension in the supplied folder.
  - Coordinates in pixels.
  - The recommended fitting software is GDSC SMLM PeakFit.
- Calibration file, similarly formatted.
  - Must contain a single fiducial, imaged at equally-spaced z-coordinates in the desired range. E.g. it could have been imaged every 50 nm from -2 um to 2 um in z. There must be one frame per z-coordinate, i.e. there must be a one-step movement between each frame.

## Inputs - parameters:
- Pixel size (in nm).
- Precision cutoff, above which the Peak Fit localisations will be filtered out.
- The step size in z in the calibration sequence.
- Filtering parameters: the initial values for filtering out pairs of localisations that are too far and too close, and the allowed deviation in the distance and intensity ratio, compared to the calibration.

## Output
Output will be written into the specified folder in the ViSP .3d format: x y z intensity frame, tab-separated.

