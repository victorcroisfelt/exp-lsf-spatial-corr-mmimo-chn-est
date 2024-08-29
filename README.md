# Exponential Spatial Correlation with Large-Scale Fading Variations in Massive MIMO Channel Estimation

This is a research-oriented code package that is primarily intended to allow readers to replicate the results of the article mentioned below and also encourage and accelerate further research on this topic:

Victor Croisfelt Rodrigues, José Carlos Marinello Filho, and Taufik Abrão, "[Exponential Spatial Correlation with Large-Scale Fading Variations in Massive MIMO Channel Estimation](https://doi.org/10.1002/ett.3563)", Transactions on Emerging Telecommunications, 2019;e3563. 

The package is based on the Matlab language and can, in fact, reproduce all the numerical results and figures discussed in the article. To contextualize, in the sequel, we present the abstract of the article and other important information.

I hope this content helps in your reaseach and contributes to building the precepts behind open science. Remarkably, in order to boost the idea of open science and further drive the evolution of science, we also motivate you to share your published results to the public.

If you have any questions and if you have encountered any inconsistency, please do not hesitate to contact me via victorcroisfelt@gmail.com.

## Abstract
To provide the vast exploitation of the large number of antennas on massive multiple‐input–multiple‐output (M‐MIMO), it is crucial to know as accurately as possible the channel state information in the base station. This knowledge is canonically acquired through channel estimation procedures conducted after a pilot signaling phase, which adopts the widely accepted time‐division duplex scheme. However, the quality of channel estimation is very impacted either by pilot contamination or by spatial correlation of the channels. There are several models that strive to match the spatial correlation in M‐MIMO channels, the exponential correlation model being one of these. To observe how the channel estimation and pilot contamination are affected by this correlated fading model, this work proposes to investigate an M‐MIMO scenario applying the standard minimum mean square error channel estimation approach over uniform linear arrays and uniform planar arrays (ULAs and UPAs, respectively) of antennas. Moreover, the elements of the array are considered to contribute unequally on the communication, owing to large‐scale fading variations over the array. Thus, it was perceived that the spatially correlated channels generated by this combined model offer a reduction of pilot contamination, consequently the estimation quality is improved. The UPA acquired better results regarding pilot contamination since it has been demonstrated that this type of array generates stronger levels of spatial correlation than the ULA. In contrast to the favorable results in channel estimation, the channel hardening effect was impaired by the spatially correlated channels, where the UPA imposes the worst performance of this effect for the discussed model.

## Content
The codes provided herein can be used to simulate the Figs. 1 to 9 contained in the article; this is done by running the scripts that have "simulation" in their names, while those with "function" in their names are called by the main scripts. Further details about each file can be found inside them.

## Acknowledgments
This research was supported in part by the Conselho Nacional de Desenvolvimento Científico e Tecnológico (CNPq) of Brazil, under grants 304066/2015‐0 and 102050/2018‐0, and in part by the Universidade Estadual de Londrina (UEL), Brazil under the PIBITI grant of the public notice 02/2018.

## Citing this Repository and License
This code is subject to the GPLv3 license. If you use any part of this repository for research, please consider to cite our aforementioned work.

```bibtex
@article{https://doi.org/10.1002/ett.3563,
author = {Croisfelt Rodrigues, Victor and Marinello, José Carlos and Abrão, Taufik},
title = {Exponential spatial correlation with large-scale fading variations in massive MIMO channel estimation},
journal = {Transactions on Emerging Telecommunications Technologies},
volume = {30},
number = {5},
pages = {e3563},
doi = {https://doi.org/10.1002/ett.3563},
url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/ett.3563},
eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/ett.3563},
note = {e3563 ett.3563}
}
