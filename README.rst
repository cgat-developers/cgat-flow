.. image:: https://travis-ci.org/cgat-developers/cgat-flow.svg?branch=master
    :target: https://travis-ci.org/cgat-developers/cgat-flow

=========
CGAT Flow
=========

We have developed a set of ruffus_ based pipelines in comparative genomics
and NGS analysis.

We are working on improving the existing documentation and portability of the code
to release a set of production pipelines soon so please stay tuned. However, to view our
current documentation please see `here <https://www.cgat.org/downloads/public/cgatpipelines/documentation/>`_

We are currently testing a script to automate the installation with conda_. Feel
free to give it a go::

    # download installation script:
    curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/master/install-devel.sh

    # see help:
    bash install-devel.sh -h

    ./install-devel.sh
    	 --install-repo
	 --install-pipeline-dependencies
	 --clone-from-repo
	 --location </full/path/to/folder/without/trailing/slash>

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

