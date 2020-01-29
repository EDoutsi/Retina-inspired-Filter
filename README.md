# Retina-inspired-Filter
========================================================================

-----------COPYRIGHT NOTICE STARTS WITH THIS LINE------------

Copyright (c) 2018 Université Côte d'Azur and 4G-TECHNOLOGY All rights reserved.

Permission is hereby granted, without written agreement and without license or royalty fees, to use, copy,  modify, and distribute this code (the source files) and its documentation for any purpose, provided that the copyright notice in its entirety appear in all copies of this code, and the original source of this code, Université Côte d'Azur (http://univ-cotedazur.fr/en#.Xe5q6S2B0UE) and 4G-TECHNOLOGY (https://4g-technology.fr) are acknowledged in any publication that reports research using this code. The research is to be cited in the bibliography as:

E. Doutsi, L. Fillatre, M. Antonini and J. Gaulmin, "Retina-Inspired Filter," in IEEE Transactions on Image Processing, vol. 27, no. 7, pp. 3484-3499, July 2018. doi: 10.1109/TIP.2018.2812079.

-----------COPYRIGHT NOTICE ENDS WITH THIS LINE------------%

Author  : Effrosyni Doutsi
Version : 1.0

Kindly report any suggestions or corrections to doutsi_efrosini@hotmail.com

========================================================================

This is a demonstration of the Retina-inspired Filter (RIF) which is a spatiotemporal transform  with the shape of a Weighted Difference of Gaussian (WDoG).  The RIF approximates the Outer Plexiform Layer (OPL) of the retina which is the innest multi-layer tissue of the early visual system.

The algorithm is described in:

E. Doutsi, L. Fillatre, M. Antonini and J. Gaulmin, "Retina-Inspired Filter".

You can change this program as you like and use it anywhere, but please refer to its original source (cite our paper).

Running on Matlab the main.m function

1. Load the image, for example

  I = imread('testimage1.tiff'); 

2. Apply the retina-inspired transform:

  [RIF_I,fftRIF_I,fftRIF_Filter] = RIF_Transform(I)
  
3. Apply the inverse of the retina-inspired transform:

  [Iout,MSE] = RIF_Inverse_Transform(I,fftRIF_Filter,fftRIF_I,10000);

Dependencies: 

MATLAB files:  main.m, GaussianKernel.m, ComputingRs.m, ComputingRc.m, RIFFilter.m, RIF_Transform.m, RIF_Inverse_Transform.m, Matrix_reshape.m  (provided with release)

Image Files: testimage1.tiff, testimage2.tiff
