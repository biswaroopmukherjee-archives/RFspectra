# RF spectra analysis

MATLAB code to plot local RF spectra across an atomic cloud in a hybrid trap. Use Functions/rfspectra.m to evaluate the spectra:

```MATLAB
[spec,clocks] = rfspectra(images,rf);
```

where `images` is a cell array of image filepaths and `rf` is a cell array of RF frequencies. You can supply an optional array of crop coordinates as a third argument `crop`. To test the program on the sample data from 2015-11-18, just run `rfspectra;` in the command line. The program crops the central part of each image, averages over the radial direction, and appends this profile to `spectra`, an array of RF freq vs axial position. The other output, `clocks` is an array of mean RF transition frequencies.

![Samples](https://raw.githubusercontent.com/biswaroopmukherjee/RFspectra/master/Figures/samples.png)

Atoms are transferred to state 3 with an RF pulse. Since k_F varies across the cloud, the RF transfer fraction depends on the position inside the cloud. Varying the RF frequency results in 'bands' of atoms in equipotential volumes transferred and visible in the arrival images (see above). 

The spectra represented as an RF vs axial position image:

![Spectra](https://raw.githubusercontent.com/biswaroopmukherjee/RFspectra/master/Figures/spectrum_out.png)

Horizontal slices through the image above are the RF spectra for a particular k_F. Below, we plot a few slices to show the clock shift:

![spectra2](https://raw.githubusercontent.com/biswaroopmukherjee/RFspectra/master/Figures/spectra.png)

The mean transition frequency varies as a function of position in the cloud. There are interesting features superimposed on an overall quadratic profile.

![Clock](https://raw.githubusercontent.com/biswaroopmukherjee/RFspectra/master/Figures/clock_out.png)

