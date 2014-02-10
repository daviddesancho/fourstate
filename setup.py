#!/usr/bin/env python

# Setup script for bestmsm package

import os
from setuptools import setup,find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
		name='FourState',
		version='0.1dev',
		description='committors of a four state model',
		url='http://github.com/daviddesancho/fourstate',
		author='David De Sancho',
		author_email='daviddesancho.at.gmail.com',
		license='GPL',
		packages=find_packages(),
		keywords= "pfold, kinetics",
		long_description=read('README.txt'),
		classifiers = ["""\
				Development Status :: 1 - Planning
				Operating System :: POSIX :: Linux
				Operating System :: MacOS
				Programming Language :: Python :: 2.7
				Topic :: Scientific/Engineering :: Bio-Informatics
				Topic :: Scientific/Engineering :: Chemistry
				"""]
		)
