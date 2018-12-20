from setuptools import setup


def readme():
    with open('README.rst') as fh:
        return fh.read()


setup(
    name = 'GridTools',
    version = '1.2',
    description = 'GRID-seq Analysis Toolkit',
    long_description = readme(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      keywords='GRID-seq, RNA–chromatin, RNA-DNA, interactome',
      url='http://github.com/biz007/gridtools',
      author='Bing Zhou',
      author_email='bzhou.sd@gmail.com',
      license='MIT',
      packages=['gridseq'],
      install_requires=[
          'markdown', 'numpy', 'pandas', 'pysam', 'h5py'
      ],
      test_suite='nose.collector',
      tests_require=['nose', 'nose-cover3'],
      entry_points={
          'console_scripts': ['GridTools=gridtools.main:main'],
      },
      include_package_data=True,
      zip_safe=False
)