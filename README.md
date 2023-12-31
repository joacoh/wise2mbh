<div id="hi" align="center">
  <h1>
    Welcome to the WISE2MBH repository!
    <img src="https://media.giphy.com/media/hvRJCLFzcasrR4ia7z/giphy.gif" width="30px"/>
  </h1>
</div>
<div id="header" align="center">
  <img src="logos/WISE2MBH positivo sin fondo.png" width="600"/>
</div>

---
### Instalation 

To install `wise2mbh-0.5.2` you will need to have `git` installed. If you don't have it, you can install it in **Linux** with the following command:

    sudo apt install git

To install `wise2mbh-0.5.2`, use the following command:

    pip install git+https://github.com/joacoh/wise2mbh.git

Pre-requisites are a Python version `>=3.8` and have `numpy-1.23.5`, `scipy-1.9.3`, `astropy-5.1.1`, `pandas-1.5.2` and `astroquery-0.4.6`

---
### Scripts and Tutorials

A solid pipeline has been developed for internal-use only on the ETHER sample and a modified version will be provided in the near future!
Right now, everybody can build a script with the provided functions and tutorials. Be aware of matrix sizes when using the MC approach.

Actual stage of public pipeline: **Ready**

Actual stage of tutorials: **Ready**

---

### About this project

- WISE2MBH is a simple an effective algorithm that uses infrared cataloged data from the Wide-field Infrared Survey Explorer (WISE) to estimate the mass of supermassive black holes (SMBH). It was presented at XVII RRLA-LARIM in Montevideo as a poster that can be seen [here](https://joacoh.github.io/talks/2023-11-29-talk) and a proceeding was submitted.

- It implements a Monte Carlo approach for error propagation, considering mean photometric errors from WISE magnitudes, errors in fits of scaling relations used and scatter of those relations, if available.

- Actual stage of publication: **Submitted to MNRAS**

- Actual stage of RRLA-LARIM proceeding: **Submitted to RevMexAA**

- Actual stage of final sample: **Awaiting acceptance in MNRAS**

- Actual stage of package: **stable-0.5.2**

---

### Acknowledgements

- I greatly appreciate the support from my collaborators: **Neil Nagar** (MSc thesis advisor), **Vicente Arratia** and **Thomas H. Jarrett**. Also to **Yuri Kovalev**, **Angelo Ricarte**, and **Dominic Pesce** for useful discussions during my visit at Black Hole Initiative in Harvard and to **Yuhan Yao** for providing comparisons to put in the paper. 
- We, as a team, acknowledge funding from **ANID Chile via Nucleo Milenio TITANs (Project NCN19-058)**, **Fondecyt Regular (Project 1221421)** and **Basal (Project FB210003)**. T.H.J. acknowledges support from the **National Research Foundation (South Africa)**.
- All lookup tables were produced by T.H.J and come originally from the publication: [A New Wide-field Infrared Survey Explorer Calibration of Stellar Mass](https://iopscience.iop.org/article/10.3847/1538-4357/acb68f/meta) in ApJ. The original TBL files can be found in the wise2mbh/kcorrections folder. If you just want to use that tables, please consider referencing that publication instead of WISE2MBH.
- Logos were made by **Benjamín Ramírez**, a very good friend of mine! You can contact him on [Instagram](https://www.instagram.com/iamtwentythreee/) or [Behance](https://www.behance.net/be23r/).