#!/usr/bin/env python

from setuptools import setup

version = '1.0.0'

required = open('requirements.txt').read().split('\n')

setup(
    name='ancpy',
    version=version,
    description='A python package for ancient DNA',
    author='Joe Marcus, Arjun Biddanda, Chi-Chun Liu',
    author_email='jhmarcus@uchicago.edu, abiddanda@uchicago.edu, chichun@uchicago.edu',
    url='https://github.com/NovembreLab/ancpy',
    packages=['ancpy'],
    install_requires=required,
    long_description='See ' + 'https://github.com/NovembreLab/ancpy',
    license='MIT'
)
