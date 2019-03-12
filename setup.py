import sys
import os
import subprocess
import re

########################################################################
########################################################################
# Import setuptools
# Use existing setuptools, otherwise try ez_setup.
try:
    import setuptools
except ImportError:
    # try to get via ez_setup
    # ez_setup did not work on all machines tested as
    # it uses curl with https protocol, which is not
    # enabled in ScientificLinux
    import ez_setup
    ez_setup.use_setuptools()

from setuptools import setup, find_packages

from distutils.version import LooseVersion
if LooseVersion(setuptools.__version__) < LooseVersion('1.1'):
    print(("Version detected:", LooseVersion(setuptools.__version__)))
    raise ImportError(
        "the CGAT code collection requires setuptools 1.1 higher")

########################################################################
########################################################################
IS_OSX = sys.platform == 'darwin'

########################################################################
########################################################################
# collect CGAT version
sys.path.insert(0, "scripts")
import version

version = version.__version__

###############################################################
###############################################################
# Check for external dependencies
#
# Not exhaustive, simply execute a representative tool from a toolkit.
external_dependencies = (
    ("wigToBigWig", "UCSC tools", 255),
    ("bedtools", "bedtools", 0),
    )

for tool, toolkit, expected in external_dependencies:
    try:
        # py3k
        from subprocess import DEVNULL
    except ImportError:
        DEVNULL = open(os.devnull, 'wb')

    try:
        retcode = subprocess.call(tool, shell=True,
                                  stdout=DEVNULL, stderr=DEVNULL)
    except OSError as msg:
        print(("WARNING: depency check for %s failed: %s" % (toolkit, msg)))

    # UCSC tools return 255 when called without arguments
    if retcode != expected:
        print(("WARNING: depency check for %s(%s) failed, error %i" %
               (toolkit, tool, retcode)))

major, minor1, minor2, s, tmp = sys.version_info
if major < 3:
    raise SystemExit("""CGAT requires Python 3 or later.""")

cgat_packages = find_packages()
cgat_package_dirs = {'cgatpipelines': 'cgatpipelines'}

# Classifiers
classifiers = """
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

setup(
    # package information
    name="cgatpipelines",
    version=version,
    description="cgat-flow : Next-generation sequencing pipelines",
    author="Andreas Heger",
    author_email="andreas.heger@gmail.com",
    license="MIT",
    platforms=["any"],
    keywords="computational genomics",
    long_description="cgatflow : Next-generation sequencing pipelines",
    classifiers=[_f for _f in classifiers.split("\n") if _f],
    url="https://github.com/cgat-developers/cgat-flow",
    # package contents
    packages=cgat_packages,
    package_dir=cgat_package_dirs,
    include_package_data=True,
    entry_points={
        "console_scripts": ["cgatflow = cgatpipelines.cgatflow:main"]
    },
    # extension modules
    ext_modules=[],
    # other options
    zip_safe=False,
    test_suite="tests",
)
