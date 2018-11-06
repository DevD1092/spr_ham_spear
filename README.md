# spr_ham_spear
As researchers, we firmly believe in the concept of reproducible research.

Thus, we release the codes of simulation, algorithm development, and testing online on GitHub for the research paper "A hamming distance and spearman-correlation based star identification algorithm" (https://ieeexplore.ieee.org/abstract/document/8374923).

This repository contains the codes for simulation and testing of the star identification algorithm based on hamming distance and spearman-correlation.

Please cite our paper (citation details below) properly, if you use the code for simulation, development or testing of a new star identification algorithm or comparison & benchmarking with the above-mentioned alogirhtm.

This code is only for academic and research purposes. Commerical use of this code is not permitted.

Citation of the paper "A hamming distance and spearman-correlation based star identification algorithm":
@article{samirbhai2018hamming,
  title={A hamming distance and spearman-correlation based star identification algorithm},
  author={Samirbhai, Mehta Deval and Chen, Shoushun and Low, Kay Soon},
  journal={IEEE Transactions on Aerospace and Electronic Systems},
  year={2018},
  publisher={IEEE}
}

# Repository details

This repository contains four folders. The information about each of the folder is described as below.

1. simulate - Codes for simulating the star images. Code files inside the directory.
  Convert_Axis_2_AttitudeMatrix.m -- For converting the ECI (Star position in the catalog) frame to the camera frame (Star sesnsor).
  Find_neighbor_star_FOV.m -- For finding the number and position of the neighboring stars in a specified FOV from the center star.
  PSF.m -- Point Spread Function simulation of the star amongst the pixels.
  Plot_sky_images.m -- For simulating the star images at a specific RA & DEC angle along with a defined FOV (this function is used by the Testing technique eventually).
  centroider.m -- Finding the centroid of the stars in the image.
  SKY2000_Magnitude6_doublestars_0.12.txt -- Star catalog (adopted from SAO) containing the star ID and it's corresponding RA, DEC and Mv information. Stars having a relative magnitude threshold (Mv) of less than 6.0 are selected for making this catalog.
  
 2. 
