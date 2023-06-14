### LOTUS : An accompanying toolbox for the article "OPEN PLATFORM FOR ULTRASOUND LOCALIZATION MICROSCOPY: PERFORMANCE ASSESSMENT OF LOCALIZATION ALGORITHMS
[![Generic badge](https://img.shields.io/badge/NBME-10.1038/s41551021008248-red.svg)](https://doi.org/10.1038/s41551-021-00824-8)
![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4343435.svg)](https://doi.org/10.5281/zenodo.4343435)

<p align="center">
<img src="https://github.com/AChavignon/PALA/blob/master/LOTUS/splashscreen.png" width="400">
</p>

###### DATE 2023.12.-VERSION 1.3
**AUTHORS: Baptiste Heiles, Arthur Chavignon, CNRS, Sorbonne Universite, INSERM.**
Directed by: Olivier Couture, Research Director CNRS, Sorbonne Universite, INSERM
Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris
Code Available under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (see https://creativecommons.org/licenses/by-nc-sa/4.0/)

Partly funded by the European Research Council under the European Union Horizon H2020 programme/ERC Consolidator grant agreement No 772786-ResolveStroke

#### ACADEMIC REFERENCES TO BE CITED
Details of the code in the article by Heiles, Chavignon, Hingot, Lopez, Teston and Couture.  
[*Performance benchmarking of microbubble-localization algorithms for ultrasound localization microscopy*, Nature Biomedical Engineering, 2022 (10.1038/s41551-021-00824-8)](https://www.nature.com/articles/s41551-021-00824-8).

General description of super-resolution in: Couture et al., [*Ultrasound localization microscopy and super-resolution: A state of the art*, IEEE UFFC 2018](https://doi.org/10.1109/TUFFC.2018.2850811).

#### CORRESPONDING AUTHORS
- Article: Baptiste Heiles (baptiste.heiles@gmail.com)
- Materials, scripts, and codes: Arthur Chavignon (arthur.chavignon.pro@gmail.com)
- Collaborations, rights and others: Olivier Couture (olivier.couture@sorbonne-universite.fr)

#### ABSTRACT
Ultrasound Localization Microscopy (ULM) is an ultrasound imaging technique that relies on the acoustic response of sub-wavelength ultrasound scatterers to map the microcirculation with an order of magnitude increase in resolution. Initially demonstrated in vitro, this technique has matured and sees implementation in vivo for vascular imaging of organs, and tumors in both animal models and humans. The performance of the localization algorithm greatly defines the quality of vascular mapping. We compiled and implemented a collection of ultrasound localization algorithms and devised three datasets in silico and in vivo to compare their performance through 18 metrics. We also present two novel algorithms designed to increase speed and performance. By openly providing a complete package to perform ULM with the algorithms, the datasets used, and the metrics, we aim to give researchers a tool to identify the optimal localization algorithm for their usage, benchmark their software and enhance the overall image quality in the field while uncovering its limits.

This article provides all materials and post-processing scripts and functions.

#### 1. MINIMUM REQUIREMENTS
LOTUS software integrates MATLAB Runtime R2021a (9.10) (MathWorks). We recommend using Microsoft Windows 10 (version 1803 or higher), with 4 GB RAM minimum, and 5 GB of free space.

#### 2. INSTALLATION

1) Launch the installer `LOTUS_v1_3.exe`
2) Follow the installers' instructions
3) Restart your computer
4) LOTUS should be now available for use in the start-up menu

#### 3. USE
Please refer to the [User Guide](https://github.com/AChavignon/PALA/blob/master/LOTUS/LOTUS_UserGuide.pdf) for a step by step explanation

#### 4. DISCLAIMER
THIS SOFTWARE IS PROVIDED BY THE AUTHORS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS AND CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
