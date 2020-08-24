from setuptools import setup
import re

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open('cgnsutilities/__init__.py').read(),
)[0]

setup(name='cgnsutilities',
      version=__version__,
      description="cgnsUtilities is a package to create, modify, and use CGNS meshes.",
      keywords='CGNS',
      author='',
      author_email='',
      url='https://github.com/mdolab/cgnsutilities',
      license='Apache 2.0',
      packages=[
          'cgnsutilities',
      ],
      package_data={
          'cgnsutilities': ['*.so']
      },
      install_requires=[
          'numpy>=1.16',
      ],
      classifiers=[
        "Operating System :: Linux",
        "Programming Language :: Python, Fortran"],
      scripts=['cgnsutilities/cgns_utils'],
      )
