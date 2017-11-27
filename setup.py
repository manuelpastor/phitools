from setuptools import setup

setup(name='phitools',
      version='0.1',
      description='Collection of simple tools used for manipulating series of chemical compounds',
      url='https://github.com/phi-grib/phitools',
      author='Manuel Pastor and Elisabet Gregori',
      author_email='manuel.pastor@upf.edu',
      license='Apache License, Version 2.0',
      packages=['phitools', 'phitools.src'],
      zip_safe=False)