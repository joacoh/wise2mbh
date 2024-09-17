<div id="hi" align="center">
  <h1>
    Welcome to the WISE2MBH repository!
    <img src="https://media.giphy.com/media/hvRJCLFzcasrR4ia7z/giphy.gif" width="30px"/>
  </h1>
</div>
<div id="header" align="center">
  <img src="logos/WISE2MBH positivo sin fondo.png" width="600"/>
</div>

<div align="center">
    <a href="https://doi.org/10.1093/mnras/stae1372">
        <img src="https://img.shields.io/badge/DOI-10.1093%2Fmnras%2Fstae1372-blue" alt="!Paper">
    </a>
    <a href="https://doi.org/10.48550/arXiv.2405.18336">
        <img src="https://img.shields.io/badge/arXiv-arXiv%3A2405.18336-orange" alt="!arxiv">
    </a>
    <a href="#">
        <img src="https://img.shields.io/badge/version-0.7.3-white" alt="!version">
    </a>
</div>


---
### Instalation 

To install `wise2mbh-0.7.3` you will need to have `git` installed. If you don't have it, you can install it in **Linux** with the following command:

    sudo apt install git

To install `wise2mbh-0.7.3`, use the following command:

    pip install git+https://github.com/joacoh/wise2mbh.git

Pre-requisites are a Python version `>=3.10` and have `numpy-1.23`, `scipy-1.12.0`, `astropy-5.3.4`, `pandas-2.2.1` and `astroquery-0.4.6`, and for this reason **we highly recommend using Wise2mbh in a dedicated Python environment with exact versions**.

---
### Scripts and Tutorials

A solid pipeline has been developed for **internal-use only on the ETHER database** and a modified version is already available for public use.
Right now, everybody can build a script with the provided functions and tutorials. Be aware of matrix sizes when using the MC approach.

- Public pipeline: Available as **pipeline.py** on main

- Tutorials: Notebooks on main/notebooks

- Documentation: **WIP**

---

### About this project

- WISE2MBH is a simple an effective algorithm that uses infrared cataloged data from the Wide-field Infrared Survey Explorer (WISE) to estimate the mass of supermassive black holes (SMBH). It was presented at XVII RRLA-LARIM in Montevideo as a poster that can be seen [here](https://joacoh.github.io/talks/2023-11-29-talk) and a proceeding was submitted.

- It implements a Monte Carlo approach for error propagation, considering mean photometric errors from WISE magnitudes, errors in fits of scaling relations used and scatter of those relations, if available.

- Main Publication: **[MNRAS](https://doi.org/10.1093/mnras/stae1372)**

  If the results, discussion, and/or code are useful for your research, **please consider referencing us**:

  ```
  @ARTICLE{2024MNRAS.531.4503H,
        author = {{Hern{\'a}ndez-Y{\'e}venes}, J. and {Nagar}, N. and {Arratia}, V. and {Jarrett}, T.~H.},
          title = "{WISE2MBH: a scaling-based algorithm for probing supermassive black hole masses through WISE catalogues}",
        journal = {\mnras},
      keywords = {Astrophysics - Astrophysics of Galaxies},
          year = 2024,
          month = jul,
        volume = {531},
        number = {4},
          pages = {4503-4523},
            doi = {10.1093/mnras/stae1372},
  archivePrefix = {arXiv},
        eprint = {2405.18336},
  primaryClass = {astro-ph.GA},
        adsurl = {https://ui.adsabs.harvard.edu/abs/2024MNRAS.531.4503H},
        adsnote = {Provided by the SAO/NASA Astrophysics Data System}
  }
  ```

- RRLA-LARIM Proceeding: **Accepted to RevMexAA**

- WISE2MBH final sample: Available in MNRAS as **supplementary material** --- Soon on CDS/Vizier

- WISE2MBH last version: **stable-0.7.3**

---

### Acknowledgements

- I greatly appreciate the support from my collaborators: **Neil Nagar** (MSc thesis advisor), **Vicente Arratia** and **Thomas H. Jarrett**. Also to **Yuri Kovalev**, **Angelo Ricarte**, **Dominic Pesce** and **Priyamvada Natarajan** for useful discussions during my visit at Black Hole Initiative in Harvard and to **Yuhan Yao** for providing comparisons to put in the paper. 
- We, as a team, acknowledge funding from **ANID Chile via Nucleo Milenio TITANs (Project NCN19-058)**, **Fondecyt Regular (Project 1221421)** and **Basal (Project FB210003)**. T.H.J. acknowledges support from the **National Research Foundation (South Africa)**.
- All lookup tables were produced by T.H.J and come originally from the publication: [A New Wide-field Infrared Survey Explorer Calibration of Stellar Mass](https://iopscience.iop.org/article/10.3847/1538-4357/acb68f/meta) in ApJ. The original TBL files can be found in the wise2mbh/kcorrections folder. If you just want to use that tables, please consider referencing that publication instead of WISE2MBH.
- Logos were made by **Benjamín Ramírez**, a very good friend of mine! You can contact him on [Instagram](https://www.instagram.com/iamtwentythreee/) or [Behance](https://www.behance.net/be23r/).