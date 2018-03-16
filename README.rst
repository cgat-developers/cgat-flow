
==========
CGAT-flow 
==========

CGAT-flow is a set of `ruffus <http://www.ruffus.org.uk/>`_ based Next generation sequencing (NGS)
pipelines that we have built using our `CGAT-core <https://github.com/cgat-developers/cgat-core>`_ workflow 
management system. For documentation please see `here <https://www.cgat.org/downloads/public/cgatpipelines/documentation>`_.

These workflows are in constant development and we encourage people to contribute to our code,
for more information please see documentation on the contributing section `here <>`_


=========================
Installation Instructions
=========================

For detailed installation instructions please refer to the `documentation <https://www.cgat.org/downloads/public/cgatpipelines/documentation/InstallingPipelines.html>`_

The preferred method to install the CGAT Pipelines is using the installation script, which uses `conda <https://conda.io/docs/>`_.

Here are the steps::

    # download installation script:
    curl -O https://raw.githubusercontent.com/CGATOxford/CGATPipelines/master/install-CGAT-tools.sh

    # see help:
    bash install-CGAT-tools.sh

    # install the development version (recommended, no production version yet):
    bash install-CGAT-tools.sh --devel --no-dashboard [--location </full/path/to/folder/without/trailing/slash>]

    # the code is downloaded in zip format by default. If you want to get a git clone, use:
   --git # for an HTTPS clone
   --git-ssh # for a SSH clone (you need to be a CGATOXford contributor on GitHub to do this)

   # the pipelines are intended to run on a cluster using the DRMAA API. If that's not your case, please use:
   --no-cluster

   # if you want to download and install IDEs like Spyder or RStudio with this installation, please use:
   --ide

   # once the installation is finished, enable the conda environment as requested by the installation script:
   source </full/path/to/folder/without/trailing/slash>/conda-install/bin/activate cgat-p

   # finally, please run the cgatflow command-line tool to check the installation:
   cgatflow --help

The installation script will put everything under the specified location. It needs 15 GB of disk space and it takes about 35 minutes to complete.
The aim of the script is to provide a portable installation that does not interfere with the existing software. As a result, you will have a
conda environment working with the CGAT Pipelines which can be enabled on demand according to your needs.
