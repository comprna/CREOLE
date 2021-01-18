from setuptools import setup, find_packages

setup (
       name='CREOLE',
       version='0.1',
       packages=find_packages(),

       # Declare your packages' dependencies here, for eg:
       install_requires=['foo>=3'],

       # Fill in these to make your Egg ready for upload to
       # PyPI
       author='Joel A. Indi',
       author_email='joeamaky@gmail.com',

       #summary = 'Just another Python package for the cheese shop',
       url='',
       license='',
       long_description='Reference-free identification of open reading frames and encoded proteins from nanopore transcriptomic long reads.',

       # could also include long_description, download_url, classifiers, etc.

  
       )