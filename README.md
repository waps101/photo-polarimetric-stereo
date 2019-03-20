# Linear Differential Constraints for Photo-polarimetric Height Estimation

This is a Matlab implementation of our ICCV 2017 paper "Linear Differential Constraints for Photo-polarimetric Height Estimation"

The main function that executes our shape estimation method is:

```matlab
[ height,mask ] = photoPolStereo( theta,Iun,phi,s,mask,Iun2,t,options );
```

The function can be called in a number of forms (corresponding to the different formulations in the paper). Use the options structure to select which constraints are included in the optimisation. At least two of: options.usephase, options.usepolratio and options.useintensityratio must be true. Depending on which constraints use use, you may not need all of the other parameters. Just pass in [] if a parameter is not needed. s and t are the two light source directions. (theta,Iun,phi) are derived from the polarisation image (theta is the zenith angle obtained by inverting the degree of polarisation equation). mask is a bindary foreground mask. Iun2 is the second unpolarised intensity image with light source direction t. Note that the mask may be modified if there are pixels where a finite difference gradient could not be estimated so you should use the modified one in any further processing. The returned variable height contains the estimated height map.

A demo script shows how to generate some data and then call the shape estimation code. This is run with:

```matlab
demoICCV;
```

## Multichannel polarisation image estimation

We also include code for our multichannel polarisation image estimation method. This is called with:

```matlab
[rho_est,phi_est,Iun_est]=AMPolarisationImage(polImages,polAng,mask);
```

Here, polImages is size [row,colum,nChannel,nImages] and polAng contains the polariser orientation angle and should have the same length as nImages. This code is for using a colour image with a single light source direction so nChannel should be 3.

## Citation

If you use this code in your research, please cite the following paper:

S. Tozza, W.A.P. Smith, D. Zhu, R. Ramamoorthi and E.R. Hancock. Linear Differential Constraints for Photo-polarimetric Height Estimation. In Proc. IEEE International Conference on Computer Vision (ICCV), pp. 2279-2287, 2017.

Bibtex:

    @inproceedings{tozza2017linear,
      title={Linear differential constraints for photo-polarimetric height estimation},
      author={Tozza, Silvia and Smith, William AP and Zhu, Dizhong and Ramamoorthi, Ravi and Hancock, Edwin R},
      booktitle={Proceedings of the IEEE International Conference on Computer Vision (ICCV)},
      pages={2279--2287},
      year={2017}
    }
