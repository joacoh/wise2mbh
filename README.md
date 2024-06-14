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

To install `wise2mbh-0.6` you will need to have `git` installed. If you don't have it, you can install it in **Linux** with the following command:

    sudo apt install git

To install `wise2mbh-0.6`, use the following command:

    pip install git+https://github.com/joacoh/wise2mbh.git

Pre-requisites are a Python version `>=3.8` and have `numpy-1.23.5`, `scipy-1.9.3`, `astropy-5.1.1`, `pandas-1.5.2` and `astroquery-0.4.6`

---
### Scripts and Tutorials

A solid pipeline has been developed for internal-use only on the ETHER sample and a modified version will be provided in the near future!
Right now, everybody can build a script with the provided functions and tutorials. Be aware of matrix sizes when using the MC approach.

Actual stage of public pipeline: **Ready**

Actual stage of tutorials: **Ready**

Actual stage of documentation: **WIP**

---

### About this project

- WISE2MBH is a simple an effective algorithm that uses infrared cataloged data from the Wide-field Infrared Survey Explorer (WISE) to estimate the mass of supermassive black holes (SMBH). It was presented at XVII RRLA-LARIM in Montevideo as a poster that can be seen [here](https://joacoh.github.io/talks/2023-11-29-talk) and a proceeding was submitted.

- It implements a Monte Carlo approach for error propagation, considering mean photometric errors from WISE magnitudes, errors in fits of scaling relations used and scatter of those relations, if available.

- Main Publication: **[MNRAS](https://doi.org/10.1093/mnras/stae1372)**

  If the results, discussion, or code is useful for your research, please use the following reference:

  ```
  @article{10.1093/mnras/stae1372,
      author = {Hernández-Yévenes, J and Nagar, N and Arratia, V and Jarrett, T H},
      title = "{WISE2MBH: A scaling-based algorithm for probing supermassive black hole masses through WISE catalogs}",
      journal = {Monthly Notices of the Royal Astronomical Society},
      pages = {stae1372},
      year = {2024},
      month = {06},
      issn = {0035-8711},
      doi = {10.1093/mnras/stae1372},
      url = {https://doi.org/10.1093/mnras/stae1372},
      eprint = {https://academic.oup.com/mnras/advance-article-pdf/doi/10.1093/mnras/stae1372/58213090/stae1372.pdf},
  }
  ```

- RRLA-LARIM Proceeding: **Accepted to RevMexAA**

- WIS2MBH final sample: **Available in MNRAS as supplementary material**

- WISE2MBH last version: **stable-0.6**

---

### Acknowledgements

- I greatly appreciate the support from my collaborators: **Neil Nagar** (MSc thesis advisor), **Vicente Arratia** and **Thomas H. Jarrett**. Also to **Yuri Kovalev**, **Angelo Ricarte**, **Dominic Pesce** and **Priyamvada Natarajan** for useful discussions during my visit at Black Hole Initiative in Harvard and to **Yuhan Yao** for providing comparisons to put in the paper. 
- We, as a team, acknowledge funding from **ANID Chile via Nucleo Milenio TITANs (Project NCN19-058)**, **Fondecyt Regular (Project 1221421)** and **Basal (Project FB210003)**. T.H.J. acknowledges support from the **National Research Foundation (South Africa)**.
- All lookup tables were produced by T.H.J and come originally from the publication: [A New Wide-field Infrared Survey Explorer Calibration of Stellar Mass](https://iopscience.iop.org/article/10.3847/1538-4357/acb68f/meta) in ApJ. The original TBL files can be found in the wise2mbh/kcorrections folder. If you just want to use that tables, please consider referencing that publication instead of WISE2MBH.
- Logos were made by **Benjamín Ramírez**, a very good friend of mine! You can contact him on [Instagram](https://www.instagram.com/iamtwentythreee/) or [Behance](https://www.behance.net/be23r/).