from setuptools import setup, find_packages

setup(
    name='ScanningTools',
    version='1.0',
    packages=find_packages(exclude=['tests*']),
    license='MIT',
    description='Some astronomical utilities for studing sky scanning strategies',
    long_description=open('README.md').read(),
    install_requires=['numpy', 'healpy', 'astropy'],
    url='https://github.com/fincardona/ScanningTools.py.git',
    author='Federico Incardona',
    author_email='incardona.federico@gmail.com'
)
