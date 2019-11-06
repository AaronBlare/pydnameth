#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

with open('requirements_dev.txt', 'r') as f:
    requirements = list(f)

setup_requirements = [ ]

test_requirements = [ ]

packages = find_packages(include=['pydnameth'])
packages.extend('pydnameth.' + item for item in find_packages(where='pydnameth'))

setup(
    author="Aaron Blare",
    author_email='aaron.blare@mail.ru',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="DNA Methylation Analysis Package",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='pydnameth',
    name='pydnameth',
    packages=packages,
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/AaronBlare/pydnameth',
    version='0.2.4',
    zip_safe=False,
)
