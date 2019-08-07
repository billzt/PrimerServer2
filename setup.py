
from setuptools import setup, find_packages

VERSION = '2.0.0b3'

setup(name='primerserver2',
      version=VERSION,
      description="a high-throughput primer design and specificity-checking platform",
      long_description=__doc__,
      classifiers=[
          'Development Status :: 4 - Beta',
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
          'primer3-py',
          'progressbar2',
          'flask'
      ],
      entry_points={
          'console_scripts': [
              'primertool = primerserver2.cmd.primertool:main',
              'primerserver-config = primerserver2.web.config:main'
          ]
      },
)