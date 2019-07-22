
from setuptools import setup, find_packages

VERSION = '2.0.0a2'

setup(name='primerserver2',
      version=VERSION,
      description="a high-throughput primer design and specificity-checking platform",
      long_description=__doc__,
      classifiers=[
          'Development Status :: 2 - Pre-Alpha',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Operating System :: Unix',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      keywords='primer bioinformatics PCR',
      author='Tao Zhu',
      author_email='taozhu@mail.bnu.edu.cn',
      url='https://github.com/billzt/PrimerServer2',
      license='MIT',
      packages=find_packages(),
      include_package_data=True,
      python_requires='>=3.6',
      install_requires=[
          'primer3-py'
      ],
      entry_points={
          'console_scripts': [
              'primertool = primerserver2.cmd.primertool:main'
          ]
      },
)