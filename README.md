# dotf-simple-example
this script computes dOTF from two sets of psf frames, one of them being captured with a small pupil amplitude modification. 

The Python script accompanies a presentation of the differential optical transfer function used as a focal plane wavefront sensor (see: Codona, J. L. (2013). Differential optical transfer function wavefront sensing. Optical Engineering, 52(9), 097105-097105).

Requirements: Python 3, NumPy, Matplotlib, and OpenCV.

This presentation took place in french language on March 4, 2026, on the ASTRO-FR Discord server, a gathering place for members of the French-speaking amateur astronomy community.

Its purpose was to present a technique also used in professional laboratories and observatories, and to showcase the initial results for amateur telescope makers.

The Python script reads two sets of diffraction spot images.

For one of the two image series, the wave amplitude was modified by a local obstruction of the pupil edge, representing slightly less than one percent of the total transmitted flux.
The two image series are otherwise identical: the observed source is a point source with a known average wavelength, and the camera response is required to be linear for recording, which is done in the rawest possible format.
To average the effects of air turbulence and increase the signal-to-noise ratio, several thousand images per series are necessary to improve the quality of the wavefront estimation.

The full width at half maximum (FWHM) of the star in the images must be at least two pixels, the maximum must be less than the camera's saturation value and the minimum above 0 (no clipping).

Roughly, the script compute the optical transfer function of the point spread function of the two images sets, subtract the first to the second and extracts a phase signal close to that of the wavefront at the output of the optical instrument to be tested, then writes the later to an intermediate file with the extension .wft. This script can be loaded by an application allowing for more detailed analysis of the results (https://github.com/githubdoe/DFTFringe).

The script, being only a showcase, is adapted to the datasets shown during the presentation and will require adjustments for other data.

Possible uses:

- In the workshop, checking slightly deformed mirrors (aspherical or not) at the center of curvature.

- In the workshop, checking an instrument in autocollimation on a reference flat mirror, potentially evaluating chromatic aberration.

- In the sky, checking a high-resolution instrument, comparing configurations (different Barlow lenses, impact of filters or ADC).

- In the sky, characterizing collimation according to the direction pointed (instrumental flexures).
