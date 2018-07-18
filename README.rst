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
=======
.. image:: https://travis-ci.org/cgat-developers/cgat-flow.svg?branch=master
    :target: https://travis-ci.org/cgat-developers/cgat-flow

=========
CGAT Flow
=========

We have developed a set of ruffus_ based pipelines in comparative genomics
and NGS analysis.

We are working on improving the existing documentation and portability of the code
to release a set of production pipelines soon so please stay tuned.

We are currently testing a script to automate the installation with conda_. Feel
free to give it a go::

        # download installation script:
        curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/master/install-CGAT-tools.sh

        # see help:
        bash install-CGAT-tools.sh

        # install the development version (recommended, no production version yet):
        bash install-CGAT-tools.sh --devel [--location </full/path/to/folder/without/trailing/slash>]

        # the code is downloaded in zip format by default. If you want to get a git clone, use:
           --git # for an HTTPS clone
           --git-ssh # for a SSH clone (you need to be a CGATOXford contributor on GitHub to do this)

        # the pipelines are intended to run on a cluster using the DRMAA API. If that's not your case, please use:
           --no-cluster

        # if you want to download and install IDEs like Spyder or RStudio with this installation, please use:
           --ide

        # once the installation is finished, enable the conda environment as requested by the installation script
        # NB: you probably want to automate this by adding the instructions below to your .bashrc
        source </full/path/to/folder/without/trailing/slash>/conda-install/etc/profile.d/conda.sh
        conda activate base
        conda activate cgat-f

        # finally, please run the cgatflow command-line tool to check the installation:
        cgatflow --help

The installation script will put everything under the specified location. It needs
15 GB of disk space and it takes about 35 minutes to complete. The aim of the
script is to provide a portable installation that does not interfere with the existing
software. As a result, you will get a conda environment working with CGAT Flow
which can be enabled on demand according to your needs.

On top of the instructions above, please make sure that you configure the following
environment variables::

        # Access to the DRMAA library: https://en.wikipedia.org/wiki/DRMAA
        export DRMAA_LIBRARY_PATH=/<full-path>/libdrmaa.so

        # You can get this value from your configured environment:
        env | grep DRMAA_LIBRARY_PATH

        # or just look for the library:
        find <path-to-DRMS-install-folder> -name "*libdrmaa.so"

        # Also, make sure you have defined temporary folders
        # 1. Local to execution hosts with
        export TMPDIR=/tmp
        # 2. Shared to pipeline working directory
        export SHARED_TMPDIR=/<path-to-network-folder>/scratch

For questions, please open a new issue on
`GitHub
<https://github.com/cgat-developers/cgat-flow/issues>`_.

.. _ruffus: http://www.ruffus.org.uk
.. _conda: https://conda.io
