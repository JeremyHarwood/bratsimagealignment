# BRATS Image Aligner (CASA Version)

Automatic pixel-based alignment tool for radio images for use in the creation of spectral ageing and spectral index maps.

If you have made use of this script, please cite Harwood, Vernstrom & Stroe 2019, MNRAS 491 803.

**Note this is the CASA region file version. If you are using DS9 style region files, there is a separate (easier to use) version available.**

This has primarily been created for use with the Broadband Radio Astronomy Tools (BRATS) but is a stand alone module and so can be used for any alignment purposes.

The module aligned images based the peak of Gaussians fitted to the specified regions. The standard alignment process is as follows:

- Create region file(s) around a bright point source as close to the target source as is practicable. The source should be visible in all images, or a staged approach should be used (see below).
- Set the required alignment parameters.
- Call the setup function.
- Call the align function.
- Ensure from the output the min, max, and mean offsets are to an acceptable level of accuracy.

Where it is not practicable for all images to be align on a single point source, multiple runs should be performed. Staged approach:
- Align an intial set of images in the standard manner.
- Either manually or via a coded loop align the second set of images to the **aligned set** using a secondary point source via one of the following options:
  - Option 1: Set "mode=0" and "reference_image" to the index of the previously aligned image.
  - Option 2: Manually obtain the Gaussian peak of the previously aligned image, and use "mode=1" with "reference_location" set to the peak value.
- Manually check that the target source has been aligned correctly.

Only bratsalign_casa.py is required to allow usage, which can then be imported in the
standard manner (#import bratsalign_casa).

See bia_casa_example.py and accompanying .fits files for example usage.

### Known issues
On certain setups when using certain version combinations of CASA and Python, it is possible that CASA will not run automatically due to a Python path conflict. In such cases, it is recommended to run the Python script created by bratsimagealigner manually. Note that the final alignment checking will not be done in such cases and should be checked by hand by the user on completion.




