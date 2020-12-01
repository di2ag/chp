import os
import sys
import re
import io

from setuptools import find_packages
from setuptools import setup

__version__ = re.search(r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        io.open('_version.py', encoding='utf_8_sig').read()).group(1)

REQUIRED_PACKAGES = [
    'pybkb'
]

setup(
    name='chp',
    version=__version__,
    author='Chase Yakaboski',
    author_email='chase.th@dartmouth.edu',
    description='NCATS Connections Hypothesis Provider',
    packages=find_packages(),
    install_requires=REQUIRED_PACKAGES,
    python_requires='>=3.8',
    dependency_links=[
        'git+https://github.com/di2ag/PyBKB.git@master#egg=pybkb-1.0.0'
    ]
)

