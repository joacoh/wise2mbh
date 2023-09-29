from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

VERSION = '0.0.1'
DESCRIPTION = 'Galaxy properties and MBH estimates from WISE'

# Setting up
setup(
    name="wise2mbh",
    version=VERSION,
    author="joacoh (Joaquín Hernández-Yévenes)",
    author_email="<jheryev@gmail.com>",
    description=DESCRIPTION,
    packages=find_packages(),
    install_requires=['numpy','astropy','scipy'],
    keywords=['python'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)