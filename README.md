# RF spectra analysis

MATLAB code to plot local RF spectra across an atomic cloud in a hybrid trap. Use Functions/rfspectra.m to evaluate the spectra:

```MATLAB
[spec,clock] = rfspectra(images,rf,crop);
```

where images is a cell array of image filepaths and rf is a cell array of rf frequencies. Crop is an optional array of crop coordinates.

![Samples](https://raw.githubusercontent.com/biswaroopmukherjee/RFspectra/master/Figures/samples.png)

Atoms are transferred to state 3 with an RF pulse. Since k_F varies across the cloud, the RF transfer fraction depends on the position inside the cloud. Varying the RF frequency results in 'bands' of atoms in equipotential volumes transferred and visible in the arrival images (see above). 

![Spectra](https://raw.githubusercontent.com/biswaroopmukherjee/RFspectra/master/Figures/spectrum_out.png)

![Clock](https://raw.githubusercontent.com/biswaroopmukherjee/RFspectra/master/Figures/clock_out.png)