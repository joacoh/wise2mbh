from setuptools import setup
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

VERSION = '0.0.2'
DESCRIPTION = 'Galaxy properties and MBH estimates from WISE'

# Setting up
setup(
    name="wise2mbh",
    version=VERSION,
    author="joacoh (Joaquin Hernandez Yevenes)",
    author_email="<jheryev@gmail.com>",
    description=DESCRIPTION,
    packages=['wise2mbh'],
    install_requires=['numpy','astropy','scipy','pandas'],
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