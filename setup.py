
from setuptools import setup, find_packages

from primerserver2.core import version

setup(name='primerserver2',
      version=version.get(),
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
      python_requires='>=3.6, <3.9',
      install_requires=[
          'primer3-py < 0.7',
          'progressbar2 < 3.56',
          'flask < 2.1 ',
          'python-dotenv < 0.20',
          'python-utils < 2.7',
          'six < 1.17'
      ],
      entry_points={
          'console_scripts': [
              'primertool = primerserver2.cmd.primertool:main',
              'primertool-junction = primerserver2.cmd.junctions:main',
              'primertool-isoform = primerserver2.cmd.isoforms:main',
              'primerserver-config = primerserver2.web.config:prepare'
          ]
      },
)
