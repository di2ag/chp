import os
import sys

from setuptools import find_packages
from setuptools import setup

REQUIRED_PACKAGES = [
    'pybkb'
]

setup(
    name='chp',
    version='1.0.0',
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

