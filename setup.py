from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

VERSION = '0.5.1'
DESCRIPTION = 'Galaxy properties and MBH estimates from WISE'

# Setting up
setup(
    name="wise2mbh",
    version=VERSION,
    author="joacoh (Joaquin Hernández Yévenes)",
    author_email="<jheryev@gmail.com>",
    description=DESCRIPTION,
    packages=find_packages(),
    url='https://github.com/joacoh/wise2mbh',
    package_data={
        'wise2mbh.kcorrections': ['*.tbl']},
    install_requires=['numpy','astropy','scipy','pandas','astroquery'],
    keywords=['python'],
    license="MIT",
    classifiers=[
        "Development Status :: 1 - Planning",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ],
    python_requires='>=3.8',
)