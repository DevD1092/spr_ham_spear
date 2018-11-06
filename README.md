# spr_ham_spear
As researchers, we firmly believe in the concept of reproducible research.

Thus, we release the codes of simulation, algorithm development, and testing online on GitHub for the research paper "A hamming distance and spearman-correlation based star identification algorithm" (https://ieeexplore.ieee.org/abstract/document/8374923).

This repository contains the codes for simulation and testing of the star identification algorithm based on hamming distance and spearman-correlation.

Please cite our paper (citation details below) properly, if you use the code for simulation, development or testing of a new star identification algorithm or comparison & benchmarking with the above-mentioned alogirhtm.

This code is only for academic and research purposes. Commerical use of this code is not permitted.

Citation of the paper "A hamming distance and spearman-correlation based star identification algorithm":
@article{samirbhai2018hamming, <br />
  title={A hamming distance and spearman-correlation based star identification algorithm},<br />
  author={Samirbhai, Mehta Deval and Chen, Shoushun and Low, Kay Soon},<br />
  journal={IEEE Transactions on Aerospace and Electronic Systems},<br />
  year={2018},<br />
  publisher={IEEE}<br />
}<br />

## Repository details

This repository contains four folders. The information about each of the folder is described as below.

#### 1. simulate - Codes for simulating the star images
  Convert_Axis_2_AttitudeMatrix.m -- For converting the ECI (Star position in the catalog) frame to the camera frame (Star sesnsor).<br />
  Find_neighbor_star_FOV.m -- For finding the number and position of the neighboring stars in a specified FOV from the center star.<br />
  PSF.m -- Point Spread Function simulation of the star amongst the pixels.<br />
  Plot_sky_images.m -- For simulating the star images at a specific RA & DEC angle along with a defined FOV (this function is used by the Testing technique eventually).<br />
  centroider.m -- Finding the centroid of the stars in the image.<br />
  SKY2000_Magnitude6_doublestars_0.12.txt -- Star catalog (adopted from SAO) containing the star ID and it's corresponding RA, DEC and Mv information. Stars having a relative magnitude threshold (Mv) of less than 6.0 are selected for making this catalog.<br />
  
#### 2. Testing - Testing as well as implementation code of the proposed star identification algorithm
  Proposed_technique.m -- Code for testing and implementation of the proposed technique in an ideal case scenario (i.e. without any positional deviation, false stars or magnitude uncertainty).
  Testing_false_stars.m -- Code for testing and implementation of the proposed technique when false stars are added to the simulated star images. The number of false stars to be added can be specified in this script.
  Testing_magnitude_uncertainty.m -- Code for testing and implementation of the proposed technique in the scenario of magnitude uncertainty in the star images. The value of the magnitude uncertainty can be specified in the script and corresponding low magnitude stars will be deleted in the simulated star image.
  Testing_positional_deviation.m -- Code for testing and implementation of the proposed technique in the when the star position is deviated from its original position. The range for the positional deviation to be introduced for each star can be specified in the script and positional deviation will be introduced by selecting a random value between the range for a star.
###### NOTE: The testing and implementation scripts utilize the simulation scripts as well as the input from the LUT and the SPD directories. So, please change the path of this input accordingly.

#### 3. SPD - Generating the SPD for the propsoed technique
  SPD_generate.m -- For generating the SPD for the propsoed technique. Specific parameters (such as the FOV, pixel size, Mv, bin_size, etc.) can be specified inside the script.
  SPD.txt -- SPD generated from the above script with a magnitude threshold (Mv) of 6.0. This SPD should be given as an input to the Testing scripts.
 
#### 4. LUT - Generating the LUT for the proposed technique
  LUT_generation.m -- For generating the LUT for the proposed technique. Specific parameters (such as the FOV, pixel size, Mv, bin_size, number of bins, etc.) can be spcecified inside the script.
  LUT.txt -- LUT generated from the above script with a magnitude threshold (Mv) of 6.0. This LUT should be given as an input to the Testing scripts.
