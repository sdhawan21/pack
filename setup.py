import sys

from distutils.core import setup

with open('README.rst') as file:
	long_description = file.read()

setup{name:'pack',
version='0.1',
description='bolometric light curve analysis',
long_description=long_description,
author='Suhail Dhawan',
url='https://github.com/sdhawan21/pack',
packages=['pack'],


}
