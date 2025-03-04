from setuptools import setup
import re

__version__ = re.findall(r"""__version__ = ["']+([0-9\.]*)["']+""", open("cgnsutilities/__init__.py").read())[0]

setup(
    name="cgnsutilities",
    version=__version__,
    description="cgnsUtilities is a package to create, modify, and use CGNS meshes.",
    keywords="CGNS",
    author="",
    author_email="",
    url="https://github.com/mdolab/cgnsutilities",
    license="Apache 2.0",
    packages=["cgnsutilities"],
    package_data={"cgnsutilities": ["*.so"]},
    install_requires=["numpy>=1.21", "scipy>=1.7"],
    extras_require={
        "testing": ["mdolab-baseclasses>=1.3", "testflo", "parameterized"],
        "advanced": ["pyspline"],
    },
    classifiers=["Operating System :: Linux", "Programming Language :: Python, Fortran"],
    entry_points={"console_scripts": ["cgns_utils = cgnsutilities.cgns_utils:main"]},
)
