OPEN PLATFORM FOR ULTRASOUND LOCALIZATION MICROSCOPY: PERFORMANCE ASSESSMENT OF LOCALIZATION ALGORITHMS

DATE 2020.04.03 - VERSION 1.0.0
AUTHROS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris  
Code Available under Creative Commons Non-Commercial 4.0
ACADEMIC REFERENCES TO BE CITED
Details of the code published in 2020 article by Heiles, Chavignon, Hingot and Couture.
Open Platform for Ultrasound Localization Microscopy: performance assessment of localization algorithms
General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018 

** ABSTRACT **
Ultrasound Localization Microscopy (ULM) is an ultrasound imaging technique that relies on the acoustic response of sub-wavelength 
ultrasound scatterers to map the microcirculation with an order of magnitude increase in resolution.Initially demonstrated in vitro, 
this technique has matured and sees implementation in vivo for vascular imaging of organs, and tumors in both animal models and humans. 
The performance of the localization algorithm greatly defines the quality of vascular mapping. We compiled and implemented a collection 
of ultrasound localization algorithms and devised three datasets in silico and in vivo to compare their performance through 18 metrics. 
We also present two novel algorithms designed to increase speed and performance. By openly providing a complete package to perform ULM 
with the algorithms, the datasets used, and the metrics, we aim to give researchers a tool to identify the optimal localization algorithm 
for their usage, benchmark their software and enhance the overall image quality in the field while uncovering its limits.


This article provides all materials and post-processing scripts and functions.

** MAIN SCRIPTS **
For each data sets (InSilicoPSF, InSilicoFlow, InVivoBrain), scripts are separated in two parts:
    - post processing script performing localization and tracking if needed
    - displaying scirpt analysing results from the localization
Scripts are located in 'paper_v4' folder

** PATH AND LOCATIONS **
In main scripts, two paths are required:
    - PARLA_addons_folder: the addons folder with custom functions
    - PARLA_data_folder: root path of data
    - SimppleTracker_folder: simpletracker folder , by Jean-Yves Tinevez

** MAIN ULM FUNCTION **
The two main function of superresolution process are
    - ULM_Superloc2D: localizing bubbles in a stack of images
    - ULM_tracking2D: pairing a list of bubbles into a set of tracks
    - Track2MatOut: generates a density image (matout) by accumulating tracks (https://github.com/tinevez/simpletracker)

** LOCALIZATION ALGORITHMS **
In this article with compared a list of 7 algorithms:
    - no_shift: no localization, bubble taken in the center of the max intensity pixel.
    - weighted average: intensity ponderation of the center.
    - interpolation: using cubic/lanczos/spline kernel.
    - gaussian fit: using a minimization algortihm and with a gaussian kernel.
    - radial symmertry: based on the minimlization of the intensity gradiant.
You can easily add you custom localization kernel in the function 'ULM_Superloc2D' (and 'ULM_Superloc2D_mesh') by adding a new 'LocMethod'.
Your custom functoin musts take in input a square image [5x5 or 7x7 pixels] centered on the bubble, and return the localization of 
the bubble in the pixel grid.
If the bubble is located at pixel [3 3], the function must return something like [3.1 3.2].









