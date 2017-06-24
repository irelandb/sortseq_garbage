import sys
import os
import re

from setuptools import setup
from setuptools import Extension



if (sys.version_info[0], sys.version_info[1]) != (2, 7):
    raise RuntimeError('sortseq is currently only compatible with Python 2.7.\nYou are using Python %d.%d' % (sys.version_info[0], sys.version_info[1]))


# main setup command
setup(
    name = 'sortseq', 
    description = 'Tools for analysis of Sort-Seq experiments.',
    version = '0.0.14',
    #long_description = readme,
    install_requires = [\
        ],
    platforms = 'Linux (and maybe also Mac OS X).',
    packages = ['sortseq_garbage'],
    package_dir = {'sortseq_garbage':'src'},
    download_url = '',
    scripts = [
            'scripts/sortseq_garbage'
            ],
    #package_data = {'sortseq':['*.txt']}, # template from weblogo version 3.4
)

