#!/usr/bin/env bash

# Travis: repo is provided by travis.
#         build conda environment
#
# Jenkins: repo is provided by jenkins
#
# Local: code can be downloaded (--clone-from-repo)
# 

# References
# http://kvz.io/blog/2013/11/21/bash-best-practices/
# http://jvns.ca/blog/2017/03/26/bash-quirks/

# exit when a command fails
set -o errexit

# exit if any pipe commands fail
set -o pipefail

# exit when your script tries to use undeclared variables
#set -o nounset

# trace what gets executed
#set -o xtrace
#set -o errtrace

# Bash traps
# http://aplawrence.com/Basics/trapping_errors.html
# https://stelfox.net/blog/2013/11/fail-fast-in-bash-scripts/

SCRIPT_NAME="$0"
SCRIPT_PARAMS="$@"

error_handler() {
    echo
    echo " ########################################################## "
    echo
    echo " An error occurred in:"
    echo
    echo " - line number: ${1}"
    shift
    echo " - exit status: ${1}"
    shift
    echo " - command: ${@}"
    echo
    echo " The script will abort now. User input was: "
    echo
    echo " ${SCRIPT_NAME} ${SCRIPT_PARAMS}"
    echo
    echo " Please copy and paste this error and report it via Git Hub: "
    echo " https://github.com/cgat-developers/cgat-flow/issues "
    print_env_vars
    echo " ########################################################## "
}

trap 'error_handler ${LINENO} $? ${BASH_COMMAND}' ERR INT TERM

# log installation information
log() {
    echo "# install-devel.sh log | `hostname` | `date` | $1 "
}

# report error and exit
report_error() {
    echo
    echo $1
    echo
    echo "Aborting."
    echo
    exit 1
}

# detect CGAT installation
detect_cgat_installation() {

    if [[ -z "$CGAT_HOME" ]] ; then

	if [[ -d "$HOME/cgat-install/conda-install" ]] ; then
	    UNINSTALL_DIR="$HOME/cgat-install"
	fi

    else

	if [[ -d "$CGAT_HOME/conda-install" ]] ; then
	    UNINSTALL_DIR="$CGAT_HOME"
	fi

    fi

} # detect_cgat_installation


# configure environment variables 
# set: CGAT_HOME, CONDA_INSTALL_DIR, CONDA_INSTALL_TYPE_PIPELINES
get_cgat_env() {

    if [[ -z $CGAT_HOME  ]] ; then
	CGAT_HOME=$HOME/cgat-install
    fi

    if [[ -z $CGATFLOW_REPO ]] ; then
	CGATFLOW_REPO="$CGAT_HOME/cgat-flow"
    fi
    
    if [[ $INSTALL_DEVEL ]] || [[ $INSTALL_PIPELINE_DEPENDENCIES ]]; then
	CONDA_INSTALL_TYPE_PIPELINES="cgat-flow.yml"
	CONDA_INSTALL_TYPE_APPS="cgat-apps.yml"
	CONDA_INSTALL_TYPE_CORE="cgat-core.yml"
    elif [[ $RUN_TESTS ]] || [[ $INSTALL_UPDATE ]] ; then
	if [[ -d $CGAT_HOME/conda-install ]] ; then
	    AUX=`find $CGAT_HOME/conda-install/envs/cgat-* -maxdepth 0`
	    CONDA_INSTALL_TYPE_PIPELINES=`basename $AUX`
	else
	    echo " The location of the CGAT code was not found (function: get_cgat_env). "
	    echo " Please install it first or use --location option with full path to your installation. "
	    echo
	    exit 1
	fi
    else
	echo
	echo " Wrong installation type! "
	echo " Installation aborted. "
	echo
	exit 1
    fi # if install type

    # set installation folder
    CONDA_INSTALL_DIR=$CGAT_HOME/conda-install

    # set conda environment name
    [[ ${CONDA_INSTALL_ENV} ]] || CONDA_INSTALL_ENV="cgat-flow"

} # get_cgat_env


# check whether the 'cgat-flow' conda environment is enabled
# and try to enable it if not
is_env_enabled() {
    # disable error checking
    set +e

    # is conda available?
    CONDA_PATH=$(which conda)

    if [[ $? -ne 0 ]] ; then
        # conda is not available
        get_cgat_env
        source ${CONDA_INSTALL_DIR}/etc/profile.d/conda.sh
    fi

    # is conda available?
    CONDA_PATH=$(which conda)

    if [[ $? -eq 0 ]] ; then
        # conda is available
        ENV_PATH=$(dirname $(dirname $CONDA_PATH))
        stat ${ENV_PATH}/envs/cgat-flow >& /dev/null
    else
        # conda is not available
        report_error " Conda can't be found! "
    fi

    # enable error checking again
    set -e
}


# setup environment variables
setup_env_vars() {

    export CFLAGS=$CFLAGS" -I$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_ENV/include -L$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_ENV/lib"
    export CPATH=$CPATH" -I$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_ENV/include -L$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_ENV/lib"
    export C_INCLUDE_PATH=$C_INCLUDE_PATH:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_ENV/include
    export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_ENV/include
    export LIBRARY_PATH=$LIBRARY_PATH:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_ENV/lib
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_ENV/lib:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_ENV/lib/R/lib

} # setup_env_vars

# print related environment variables
print_env_vars() {

    echo
    echo " Debugging: "
    echo " CFLAGS: "$CFLAGS
    echo " CPATH: "$CPATH
    echo " C_INCLUDE_PATH: "$C_INCLUDE_PATH
    echo " CPLUS_INCLUDE_PATH: "$CPLUS_INCLUDE_PATH
    echo " LIBRARY_PATH: "$LIBRARY_PATH
    echo " LD_LIBRARY_PATH: "$LD_LIBRARY_PATH
    echo " CGAT_HOME: "$CGAT_HOME
    echo " CGATFLOW_DIR: "$CGATFLOW_DIR
    echo " CGATFLOW_REPO: "$CGATFLOW_REPO
    echo " CONDA_INSTALL_DIR: "$CONDA_INSTALL_DIR
    echo " CONDA_INSTALL_TYPE_CORE:"$CONDA_INSTALL_TYPE_CORE
    echo " CONDA_INSTALL_TYPE_APPS: "$CONDA_INSTALL_TYPE_APPS
    echo " CONDA_INSTALL_TYPE_PIPELINES: "$CONDA_INSTALL_TYPE_PIPELINES
    echo " CONDA_INSTALL_ENV: "$CONDA_INSTALL_ENV
    echo " PYTHONPATH: "$PYTHONPATH
    [[ ! $RUN_TESTS ]] && echo " CGATFLOW_BRANCH: "$CGATFLOW_BRANCH
    [[ ! $RUN_TESTS ]] && echo " CGATAPPS_BRANCH: "$CGATAPPS_BRANCH
    [[ ! $RUN_TESTS ]] && echo " CGATCORE_BRANCH: "$CGATCORE_BRANCH
    [[ ! $RUN_TESTS ]] && echo " RELEASE: "$RELEASE
    echo " CODE_DOWNLOAD_TYPE: "$CODE_DOWNLOAD_TYPE
    echo " CLUSTER: "$CLUSTER
    echo

} # print_env_vars

# Travis installations are running out of RAM
# with large conda installations. Issue has been submitted here:
# https://github.com/conda/conda/issues/1197
# While we wait for a response, we'll try to clean up the conda
# installation folder as much as possible
conda_cleanup() {
    conda clean --index-cache
    conda clean --lock
    conda clean --tarballs -y
    conda clean --packages -y
}


# proceed with conda installation
conda_install() {

    log "installing conda"

    detect_cgat_installation

    if [[ -n "$UNINSTALL_DIR" ]] ; then

	echo
	echo " An installation of the CGAT code was found in: $UNINSTALL_DIR"
	echo " Please use --location to install CGAT code in a different location "
	echo " or uninstall the current version before proceeding."
	echo
	echo " Installation is aborted."
	echo
	exit 1
    fi

    # get environment variables: CGAT_HOME, CONDA_INSTALL_DIR, CONDA_INSTALL_TYPE_PIPELINES
    get_cgat_env

    mkdir -p $CGAT_HOME
    cd $CGAT_HOME

    # select Miniconda bootstrap script depending on Operating System
    MINICONDA=

    if [[ `uname` == "Linux" ]] ; then

	# Conda 4.4 breaks everything again!
	# Conda 4.5 looks better
	MINICONDA="Miniconda3-latest-Linux-x86_64.sh"
	#MINICONDA="Miniconda3-4.3.31-Linux-x86_64.sh"

    elif [[ `uname` == "Darwin" ]] ; then

	# Conda 4.4 breaks everything again!
	# Conda 4.5 looks better
	MINICONDA="Miniconda3-latest-MacOSX-x86_64.sh"
	#MINICONDA="Miniconda3-4.3.31-MacOSX-x86_64.sh"

    else

	echo
	echo " Unsupported operating system detected. "
	echo
	echo " Aborting installation... "
	echo
	exit 1

    fi

    log "downloading miniconda"
    # download and install conda
    curl -O https://repo.continuum.io/miniconda/${MINICONDA}

    log "installing miniconda"

    bash ${MINICONDA} -b -p $CONDA_INSTALL_DIR
    source ${CONDA_INSTALL_DIR}/etc/profile.d/conda.sh
    hash -r

    # install cgat environment
    # Conda 4.4 breaks everything again!
    # Conda 4.5 looks better
    #conda install --quiet --yes 'conda=4.3.33'
    conda update --all --yes
    conda info -a

    log "installing conda CGAT environment"

    # Now using conda environment files:
    # https://conda.io/docs/using/envs.html#use-environment-from-file

    log "curl -o env-cgat-core.yml -O https://raw.githubusercontent.com/cgat-developers/cgat-core/${CGATCORE_BRANCH}/conda/environments/${CONDA_INSTALL_TYPE_CORE}"
    curl -o env-cgat-core.yml -O https://raw.githubusercontent.com/cgat-developers/cgat-core/${CGATCORE_BRANCH}/conda/environments/${CONDA_INSTALL_TYPE_CORE}

    log "curl -o env-cgat-apps.yml -O https://raw.githubusercontent.com/cgat-developers/cgat-apps/${CGATAPPS_BRANCH}/conda/environments/${CONDA_INSTALL_TYPE_APPS}"
    curl -o env-cgat-apps.yml -O https://raw.githubusercontent.com/cgat-developers/cgat-apps/${CGATAPPS_BRANCH}/conda/environments/${CONDA_INSTALL_TYPE_APPS}

    log "curl -o env-cgat-flow.yml -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${TRAVIS_BRANCH}/conda/environments/${CONDA_INSTALL_TYPE_PIPELINES}"
    curl -o env-cgat-flow.yml -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${TRAVIS_BRANCH}/conda/environments/${CONDA_INSTALL_TYPE_PIPELINES}

    [[ ${CLUSTER} -eq 0 ]] && sed -i'' -e '/drmaa/d' env-cgat-flow.yml

    log "creating environment with cgat-core dependencies"
    conda env create --name ${CONDA_INSTALL_ENV} --file env-cgat-core.yml

    log "updating with cgat-apps dependencies"
    conda env update --name ${CONDA_INSTALL_ENV} --file env-cgat-apps.yml

    log "updating with cgat-flow dependencies"
    conda env update --name ${CONDA_INSTALL_ENV} --file env-cgat-flow.yml

    # activate cgat environment
    conda activate ${CONDA_INSTALL_ENV}

    log "installing CGAT code into conda environment"
    DEV_RESULT=0

    # install extra deps
    install_extra_deps

    # Set up other environment variables
    #setup_env_vars

    # install cgat-core
    install_cgat_core

    # install cgat-apps
    install_cgat_apps

    # make sure you are in the CGAT_HOME folder
    cd $CGAT_HOME

    # download the code out of jenkins
    if [[ ${CLONE_FROM_REPO} -eq 1 ]] ; then

	if [[ $CODE_DOWNLOAD_TYPE -eq 0 ]] ; then
	    # get the latest version from Git Hub in zip format
	    curl -LOk https://github.com/cgat-developers/cgat-flow/archive/$CGATFLOW_BRANCH.zip
	    unzip $CGATFLOW_BRANCH.zip
	    rm $CGATFLOW_BRANCH.zip
	    if [[ ${RELEASE} ]] ; then
		NEW_NAME=`echo $CGATFLOW_BRANCH | sed 's/^v//g'`
		mv cgat-flow-$NEW_NAME/ $CGATFLOW_REPO
	    else
		mv cgat-flow-$CGATFLOW_BRANCH/ $CGATFLOW_REPO
	    fi
	elif [[ $CODE_DOWNLOAD_TYPE -eq 1 ]] ; then
	    # get latest version from Git Hub with git clone
	    git clone --branch=$CGATFLOW_BRANCH https://github.com/cgat-developers/cgat-flow.git $CGATFLOW_REPO
	elif [[ $CODE_DOWNLOAD_TYPE -eq 2 ]] ; then
	    # get latest version from Git Hub with git clone
	    git clone --branch=$CGATFLOW_BRANCH git@github.com:cgat-developers/cgat-flow.git $CGATFLOW_REPO
	else
	    report_error " Unknown download type for CGAT code... "
	fi

	# make sure you are in the CGAT_HOME/cgat-flow folder
	cd $CGATFLOW_REPO
    else
	log "using existing cgat-flow repo in $CGATFLOW_REPO"
	cd "$CGATFLOW_REPO"
    fi

    # Python preparation
    log "linking cgat-flow code into conda environment"
    sed -i'' -e '/REPO_REQUIREMENT/,/pass/d' setup.py
    sed -i'' -e '/# dependencies/,/dependency_links=dependency_links,/d' setup.py
    python setup.py develop

    if [[ $? -ne 0 ]] ; then
	echo
	echo " There was a problem doing: 'python setup.py develop' "
	echo " Installation did not finish properly. "
	echo 
	echo " Please submit this issue via Git Hub: "
	echo " https://github.com/cgat-developers/cgat-flow/issues "
	echo

	print_env_vars

    fi # if-$?

    # revert setup.py if downloaded with git
    [[ $CODE_DOWNLOAD_TYPE -ge 1 ]] && git checkout -- setup.py

    # environment pinning
    # python scripts/conda.py

    # check whether conda create went fine
    if [[ $DEV_RESULT -ne 0 ]] ; then
	echo
	echo " There was a problem installing the code with conda. "
	echo " Installation did not finish properly. "
	echo
	echo " Please submit this issue via Git Hub: "
	echo " https://github.com/cgat-developers/cgat-flow/issues "
	echo

	print_env_vars

    else
	clear
	echo 
	echo " The code successfully installed!"
	echo
	echo " To activate the CGAT environment type: "
	echo " $ source $CONDA_INSTALL_DIR/etc/profile.d/conda.sh"
	echo " $ conda activate base"
	echo " $ conda activate $CONDA_INSTALL_ENV"
	echo
	echo " To deactivate the environment, use:"
	echo " $ conda deactivate"
	echo
    fi # if-$ conda create

} # conda install


# install extra dependencies
install_extra_deps() {

    log "install extra deps"

    curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${TRAVIS_BRANCH}/conda/environments/pipelines-extra.yml

    conda env update --name ${CONDA_INSTALL_ENV} --file pipelines-extra.yml

}

# install dependencies for running the pipelines
install_pipeline_deps() {

    cd $CGAT_HOME

    get_cgat_env
    
    # activate cgat environment
    is_env_enabled
    
    log "install pipeline deps"

    curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${TRAVIS_BRANCH}/conda/environments/cgat-flow-pipelines.yml

    conda env update --name ${CONDA_INSTALL_ENV} --file cgat-flow-pipelines.yml

    install_extra_envs
}

# install Python 2 dependencies in a different conda environment
install_extra_envs() {

    log "install extra environments"

    curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${TRAVIS_BRANCH}/conda/environments/pipelines-macs2.yml

    conda env update --file pipelines-macs2.yml

    curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${TRAVIS_BRANCH}/conda/environments/pipelines-tophat2.yml

    conda env update --file pipelines-tophat2.yml

    curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${TRAVIS_BRANCH}/conda/environments/pipeline-peakcalling-sicer.yml

    conda env update --file pipeline-peakcalling-sicer.yml

    curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${TRAVIS_BRANCH}/conda/environments/pipelines-splicing.yml

    conda env update --file pipelines-splicing.yml

    curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${TRAVIS_BRANCH}/conda/environments/pipelines-sailfish.yml
    
    conda env update --file pipelines-sailfish.yml

    curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${TRAVIS_BRANCH}/conda/environments/pipelines-salmon.yml
    
    conda env update --file pipelines-salmon.yml
    
}


# helper function to install cgat-apps
install_cgat_apps() {

    log "install cgat apps"

    OLDWD=`pwd`
    cd $CGAT_HOME

    if [[ $CODE_DOWNLOAD_TYPE -eq 0 ]] ; then
	# get the latest version from Git Hub in zip format
	curl -LOk https://github.com/cgat-developers/cgat-apps/archive/$CGATAPPS_BRANCH.zip
	unzip $CGATAPPS_BRANCH.zip
	rm $CGATAPPS_BRANCH.zip
	if [[ ${RELEASE} ]] ; then
	    NEW_NAME=`echo $CGATAPPS_BRANCH | sed 's/^v//g'`
	    mv cgat-apps-$NEW_NAME/ cgat-apps/
	else
	    mv cgat-apps-$CGATAPPS_BRANCH/ cgat-apps/
	fi
    elif [[ $CODE_DOWNLOAD_TYPE -eq 1 ]] ; then
	# get latest version from Git Hub with git clone
	git clone --branch=$CGATAPPS_BRANCH https://github.com/cgat-developers/cgat-apps.git
    elif [[ $CODE_DOWNLOAD_TYPE -eq 2 ]] ; then
	# get latest version from Git Hub with git clone
	git clone --branch=$CGATAPPS_BRANCH git@github.com:cgat-developers/cgat-apps.git
    else
	report_error " Unknown download type for CGAT code... "
    fi

    cd cgat-apps/

    # remove install_requires (no longer required with conda package)
    sed -i'' -e '/REPO_REQUIREMENT/,/pass/d' setup.py
    sed -i'' -e '/# dependencies/,/dependency_links=dependency_links,/d' setup.py
    python setup.py develop

    if [[ $? -ne 0 ]] ; then
	echo
	echo " There was a problem doing: 'python setup.py develop' "
	echo " Installation did not finish properly. "
	echo
	echo " Please submit this issue via Git Hub: "
	echo " https://github.com/cgat-developers/cgat-flow/issues "
	echo
	print_env_vars

    fi # if-$?

    # revert setup.py if downloaded with git
    [[ $CODE_DOWNLOAD_TYPE -ge 1 ]] && git checkout -- setup.py

    # go back to old working directory
    cd $OLDWD

} # install_cgat_apps

# helper function to install cgat-core
install_cgat_core() {

    log "install cgat core"

    OLDWD=`pwd`
    cd $CGAT_HOME

    if [[ $CODE_DOWNLOAD_TYPE -eq 0 ]] ; then
	# get the latest version from Git Hub in zip format
	curl -LOk https://github.com/cgat-developers/cgat-core/archive/$CGATCORE_BRANCH.zip
	unzip $CGATCORE_BRANCH.zip
	rm $CGATCORE_BRANCH.zip
	if [[ ${RELEASE} ]] ; then
	    NEW_NAME=`echo $CGATCORE_BRANCH | sed 's/^v//g'`
	    mv cgat-core-$NEW_NAME/ cgat-core/
	else
	    mv cgat-core-$CGATCORE_BRANCH/ cgat-core/
	fi
    elif [[ $CODE_DOWNLOAD_TYPE -eq 1 ]] ; then
	# get latest version from Git Hub with git clone
	git clone --branch=$CGATCORE_BRANCH https://github.com/cgat-developers/cgat-core.git
    elif [[ $CODE_DOWNLOAD_TYPE -eq 2 ]] ; then
	# get latest version from Git Hub with git clone
	git clone --branch=$CGATCORE_BRANCH git@github.com:cgat-developers/cgat-core.git
    else
	report_error " Unknown download type for CGAT core... "
    fi

    cd cgat-core/

    # remove install_requires (no longer required with conda package)
    sed -i'' -e '/REPO_REQUIREMENT/,/pass/d' setup.py
    sed -i'' -e '/# dependencies/,/dependency_links=dependency_links,/d' setup.py
    python setup.py develop

    if [[ $? -ne 0 ]] ; then
	echo
	echo " There was a problem doing: 'python setup.py develop' "
	echo " Installation did not finish properly. "
	echo
	echo " Please submit this issue via Git Hub: "
	echo " https://github.com/cgat-developers/cgat-apps/issues "
	echo
	print_env_vars

    fi # if-$?

    # revert setup.py if downloaded with git
    [[ $CODE_DOWNLOAD_TYPE -ge 1 ]] && git checkout -- setup.py

    # go back to old working directory
    cd $OLDWD

} # install_cgat_core


# test code with conda install
conda_test() {

    log "starting conda_test"

    # get environment variables: CGAT_HOME, CONDA_INSTALL_DIR, CONDA_INSTALL_TYPE_PIPELINES
    get_cgat_env

    #setup_env_vars

    # setup environment and run tests
    if [[ $TRAVIS_INSTALL ]] || [[ $JENKINS_INSTALL ]] ; then

        # activate cgat environment
        is_env_enabled

	# show conda environment used for testing
	log "conda env export"
	conda env export

	# install cgat-core
	install_cgat_core

	# install cgat-apps
	install_cgat_apps

	# python preparation
	log "install CGAT code into conda environment"
	cd $CGAT_HOME

	# Python preparation
	sed -i'' -e 's/CGATScripts/scripts/g' setup.py
	sed -i'' -e '/REPO_REQUIREMENT/,/pass/d' setup.py
	sed -i'' -e '/# dependencies/,/dependency_links=dependency_links,/d' setup.py
	python setup.py develop

	log "starting tests"
	# run nosetests
	if [[ $TEST_ALL ]] ; then
	    log "test_style.py" && nosetests -v tests/test_style.py && \
		log "test_import.py" && nosetests -v tests/test_import.py && \
		echo -e "restrict:\n    manifest:\n" > tests/_test_commandline.yaml && \
		log "test_commandline" && nosetests -v tests/test_commandline.py && \
		log "test_scripts" && nosetests -v tests/test_scripts.py ;
	elif [[ $TEST_IMPORT ]] ; then
	    nosetests -v tests/test_import.py ;
	elif [[ $TEST_STYLE ]] ; then
	    nosetests -v tests/test_style.py ;
	elif [[ $TEST_CMDLINE ]] ; then
	    echo -e "restrict:\n    manifest:\n" > tests/_test_commandline.yaml
	    nosetests -v tests/test_commandline.py ;
	elif [[ $TEST_PRODUCTION_SCRIPTS  ]] ; then
	    echo -e "restrict:\n    manifest:\n" > tests/_test_scripts.yaml
	    nosetests -v tests/test_scripts.py ;
	else
	    nosetests -v tests/test_scripts.py ;
	fi

    else

	if [[ $CONDA_INSTALL_TYPE_PIPELINES ]] ; then

            # activate cgat environment
            is_env_enabled

            # show conda environment used for testing
            log "conda env export"
            conda env export

	    # make sure you are in the CGAT_HOME/cgat-flow folder
	    cd $CGATFLOW_REPO

	    OUTPUT_DIR=`pwd`

	    # run tests, note: OSX has no options (no -v, -o)
	    /usr/bin/time nosetests -v tests/test_import.py
	    if [[ $? -eq 0 ]] ; then
		echo
		echo " test_import.py passed successfully! "
		echo
	    else
		echo
		echo " test_import.py failed. Please see $OUTPUT_DIR/test_import.out file for detailed output. "
		echo

		print_env_vars

	    fi

	    /usr/bin/time nosetests -v tests/test_scripts.py
	    if [[ $? -eq 0 ]] ; then
		echo
		echo " test_scripts.py passed successfully! "
		echo
	    else
		echo
		echo " test_scripts.py failed. Please see $OUTPUT_DIR/test_scripts.out file for detailed output. "
		echo

		print_env_vars

	    fi
	    
	else
	    echo
	    echo " There was an error running the tests. "
	    echo " Execution aborted. "
	    echo

	    print_env_vars

	    exit 1
	fi

    fi # if travis or jenkins

} # conda_test


# update conda installation
conda_update() {

    # get environment variables: CGAT_HOME, CONDA_INSTALL_DIR, CONDA_INSTALL_TYPE_PIPELINES
    get_cgat_env

    # activate cgat environment
    is_env_enabled

    conda update --all

    if [[ ! $? -eq 0 ]] ; then

	echo
	echo " There was a problem updating the installation. "
	echo 
	echo " Please submit this issue via Git Hub: "
	echo " https://github.com/cgat-developers/cgat-flow/issues "
	echo 

    else 

	echo
	echo " All packages were succesfully updated. "
	echo 

    fi 

} # conda_update


# unistall CGAT code collection
uninstall() {

    detect_cgat_installation

    if [[ -z "$UNINSTALL_DIR" ]] ; then

	echo
	echo " The location of the CGAT code was not found. "
	echo " Please uninstall manually."
	echo
	exit 1
	
    else

	rm -rf $UNINSTALL_DIR
	if [[ $? -eq 0 ]] ; then
	    echo
	    echo " CGAT code successfully uninstalled."
	    echo 
	    exit 0
	else
	    echo
	    echo " There was a problem uninstalling the CGAT code."
	    echo " Please uninstall manually."
	    echo
	    exit 1
	fi
    fi

}


# test whether --git and --git-ssh download is doable
test_git() {
    git --version >& /dev/null || GIT_AVAIL=$?
    if [[ $GIT_AVAIL -ne 0 ]] ; then
	echo
	echo " Git is not available but --git or --git-ssh option was given."
	echo " Please rerun this script on a computer with git installed "
	echo " or try again without --git or --git-ssh"
	report_error " "
    fi
}


# test whether --git-ssh download is doable
test_git_ssh() {
    ssh-add -L >& /dev/null || SSH_KEYS_LOADED=$?
    if [[ $SSH_KEYS_LOADED -ne 0 ]] ; then
	echo
	echo " Please load your ssh keys for GitHub before proceeding!"
	echo
	echo " Try: "
	echo " 1. eval \$(ssh-agent)"
	echo " 2. ssh-add ~/.ssh/id_rsa # or the file where your private key is"
	report_error " and run this script again. "
    fi
}


# don't mix branch and release options together
test_mix_branch_release() {
    # don't mix branch and release options together
    if [[ $RELEASE ]] ; then
	if [[ "$CGATFLOW_BRANCH" != "master" ]] || [[ $CGATAPPS_BRANCH != "master" ]] ; then
            echo
            echo " You cannot mix git branches and releases for the installation."
            echo
            echo " Your input was: "$SCRIPT_PARAMS
            report_error " Please either use branches or releases but not both."
	fi
    fi
}


# test whether a branch exists in the cgat-core repository
# https://stackoverflow.com/questions/12199059/how-to-check-if-an-url-exists-with-the-shell-and-probably-curl
test_core_branch() {
    RELEASE_TEST=0
    curl --output /dev/null --silent --head --fail https://raw.githubusercontent.com/cgat-developers/cgat-core/${CGATCORE_BRANCH}/README.md || RELEASE_TEST=$?
    if [[ ${RELEASE_TEST} -ne 0 ]] ; then
	echo
	echo " The branch provided for cgat-core does not exist: ${CGATCORE_BRANCH}"
	echo
	echo " Please have a look at valid branches here: "
	echo " https://github.com/cgat-developers/cgat-core/branches"
	echo
	report_error " Please use a valid branch and try again."
    fi
}


# test whether a branch exists in the cgat-apps repository
# https://stackoverflow.com/questions/12199059/how-to-check-if-an-url-exists-with-the-shell-and-probably-curl
test_apps_branch() {
    RELEASE_TEST=0
    curl --output /dev/null --silent --head --fail https://raw.githubusercontent.com/cgat-developers/cgat-apps/${CGATAPPS_BRANCH}/README.rst || RELEASE_TEST=$?
    if [[ ${RELEASE_TEST} -ne 0 ]] ; then
	echo
	echo " The branch provided for cgat-apps does not exist: ${CGATAPPS_BRANCH}"
	echo
	echo " Please have a look at valid branches here: "
	echo " https://github.com/cgat-developers/cgat-apps/branches"
	echo
	report_error " Please use a valid branch and try again."
    fi
}


# test whether a release exists or not
# https://stackoverflow.com/questions/12199059/how-to-check-if-an-url-exists-with-the-shell-and-probably-curl
test_release() {
    RELEASE_PIPELINES=0
    # check pipelines
    curl --output /dev/null \
	 --silent --head --fail \
	 https://raw.githubusercontent.com/cgat-developers/cgat-flow/${RELEASE}/README.rst || RELEASE_PIPELINES=$?

    if [[ ${RELEASE_PIPELINES} -ne 0 ]] ; then
	echo
	echo " The release number provided for the pipelines does not exist: ${RELEASE}"
	echo
	echo " Please have a look at valid releases here: "
	echo " https://github.com/cgat-developers/cgat-flow/releases"
	echo
	echo " An example of valid release is: --release v0.4.0"
	report_error " Please use a valid release and try again."
    fi
    
    RELEASE_SCRIPTS=0
    # check scripts
    curl --output /dev/null \
	 --silent --head --fail \
	 https://raw.githubusercontent.com/cgat-developers/cgat-apps/${RELEASE}/README.rst || RELEASE_SCRIPTS=$?
    
    if [[ ${RELEASE_SCRIPTS} -ne 0 ]] ; then
	echo
	echo " The release number provided for the scripts does not exist: ${RELEASE}"
	echo
	echo " Please have a look at valid releases here: "
	echo " https://github.com/cgat-developers/cgat-apps/releases"
	echo
	echo " An example of valid release is: --release v0.4.0"
	report_error " Please use a valid release and try again."
    fi
}


# test whether a C/C++ compiler is available
test_compilers() {
    which gcc &> /dev/null || report_error " C compiler not found "
    which g++ &> /dev/null || report_error " C++ compiler not found "
}


# clean up environment
# deliberately use brute force
cleanup_env() {
    set +e
    conda deactivate >& /dev/null || true
    conda deactivate >& /dev/null || true
    unset -f conda || true
    unset PYTHONPATH || true
    # Next actions disabled. Please see:
    # https://github.com/cgat-developers/cgat-core/issues/44
    #module purge >& /dev/null || true
    #mymodule purge >& /dev/null || true
    set -e
}


# function to display help message
help_message() {
    cat << EOF
usage="$(basename "$0") cgatflow install script

where:
    -h/--help			show this help text. The following options are supported:

    --install-repo   		install cgatflow from repository. This will also install
    		     		conda and the conda environment.

    --install-pipeline-dependencies install dependencies for running the pipelines
                     		included in the cgat-flow repository.

    --clone-from-repo  		clone a fresh copy from git for the install. Please specify
                       		either --clone-from-repo or --use-repo.

    --use-repo PATH    		install from an existing repository in PATH, see above.

    --location PATH    		installation PATH for the conda environment. Must not exist.

    --cgatflow-branch BRANCH   	branch to use for cgatflow repository checkout. Default
                               	is master.

    --cgatapps-branch BRANCH   	branch to install from cgatapps. Default is master.

    --cgatcore-branch BRANCH   	branch to install from cgatapps. Default is master.

    --env-name NAME            	name of conda environment. Default is 'cgat-flow'

Examples:

For an install from a clone in the current directory into a new conda installation, type:

    ./install-devel.sh
    	 --install-repo
	 --install-pipeline-dependencies
	 --use-repo .
	 --location </full/path/to/folder/without/trailing/slash>

To install a specific branch of the code on GitHub from scratch:

    curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/master/install-devel.sh

    ./install-devel.sh
    	 --install-repo
	 --install-pipeline-dependencies
	 --clone-from-repo
	 --cgatflow-branch <branch>
	 --location </full/path/to/folder/without/trailing/slash>

Please submit any issues via Git Hub:

https://github.com/cgat-developers/cgat-flow/issues

EOF

} # help_message

# the script starts here

cleanup_env
test_compilers

if [[ $# -eq 0 ]] ; then

    help_message

fi

# conda installation type. If set to 1, cgat-core and cgat-apps
# will be installed from repository.
INSTALL_DEVEL=
# If set to 1, install all pipeline dependencies as well
INSTALL_PIPELINE_DEPENDENCIES=
# whether or not to clone from repo
CLONE_FROM_REPO=
# pre-existing directory of a cgat-flow repo
CGATFLOW_REPO=
# test current installation
RUN_TESTS=
# update current installation
INSTALL_UPDATE=
# uninstall CGAT code
UNINSTALL=
UNINSTALL_DIR=
# where to install CGAT code
# default value is in HOME
CGAT_HOME=$HOME/cgat-install
# how to download CGAT code:
# 0 = as zip (default)
# 1 = git clone with https
# 2 = git clone with ssh
CODE_DOWNLOAD_TYPE=0
# which github branch to use (default: master)
CGATFLOW_BRANCH="master"
CGATAPPS_BRANCH="master"
CGATCORE_BRANCH="master"
# type of installation
CONDA_INSTALL_TYPE_PIPELINES=
CONDA_INSTALL_TYPE_APPS=
CONDA_INSTALL_TYPE_CORE=
# rename conda environment
CONDA_INSTALL_ENV=
# Use cluster?
# 0 = no
# 1 = yes (default)
CLUSTER=1
# Install a released version?
RELEASE=

if [[ `uname` == "Darwin" ]] ; then
    READLINK=greadlink
else
    READLINK=readlink
fi

# parse input parameters
# https://stackoverflow.com/questions/402377/using-getopts-in-bash-shell-script-to-get-long-and-short-command-line-options
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash

while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
	
	-h)
	    help_message
	    exit 0
	    ;;
	    
	--help)
	    help_message
	    exit 0
	    ;;

	--install)
	    INSTALL_DEVEL=1
	    shift
	    ;;

	--install-repo)
	    INSTALL_DEVEL=1
	    shift
	    ;;

	--install-pipeline-dependencies)
	    INSTALL_PIPELINE_DEPENDENCIES=1
	    shift
	    ;;

	--clone-from-repo)
	    CLONE_FROM_REPO=1
	    shift # past argument
	    ;;

	--location)
            CGAT_HOME=$("$READLINK" --canonicalize "$2")
	    shift 2
	    ;;

	--use-repo)
	    CGATFLOW_REPO="$2"
	    shift 2
	    ;;
	
	--cgatflow-branch)
	    CGATFLOW_BRANCH="$2"
	    test_mix_branch_release
	    shift 2
	    ;;

	--cgatapps-branch)
	    CGATAPPS_BRANCH="$2"
	    test_mix_branch_release
	    shift 2
	    ;;

	--cgatcore-branch)
	    CGATCORE_BRANCH="$2"
	    test_core_branch
	    shift 2
	    ;;

	--env-name)
	    CONDA_INSTALL_ENV="$2"
	    shift 2
	    ;;

	--zip)
	    CODE_DOWNLOAD_TYPE=0
	    shift
	    ;;

	--git)
	    CODE_DOWNLOAD_TYPE=1
	    shift
	    test_git
	    ;;

	--git-ssh)
	    CODE_DOWNLOAD_TYPE=2
	    shift
	    test_git
	    test_git_ssh
	    ;;

	--test)
	    RUN_TESTS=1
	    shift
	    ;;

	--update)
	    INSTALL_UPDATE=1
	    shift
	    ;;

	--uninstall)
	    UNINSTALL=1
	    shift
	    ;;

	--no-cluster)
	    CLUSTER=0
	    shift
	    ;;

	--release)
	    RELEASE="$2"
	    test_mix_branch_release
	    test_release
	    CGATFLOW_BRANCH="$2"
	    CGATAPPS_BRANCH="$2"
	    shift 2
	    ;;

	*)
	    echo
	    echo
	    echo " Wrong input: ${SCRIPT_NAME} ${SCRIPT_PARAMS}"
	    echo
	    help_message
	    ;;

    esac
done


# sanity check: make sure one installation option is selected
if [[ -z $RUN_TESTS ]] && \
       [[ -z $INSTALL_DEVEL ]] && \
       [[ -z $INSTALL_PIPELINE_DEPENDENCIES ]] && \
       [[ -z $INSTALL_TRAVIS ]] ; then
    report_error " You need to select either --run-tests or --install-repo"
fi

# sanity check: make sure one installation option for repo is selected
if [[ ! -z $INSTALL_DEVEL ]] && \
       [[ -z $CGATFLOW_REPO ]] && \
       [[ -z $CLONE_FROM_REPO ]] ; then
    report_error " You need to select either --clone-from-repo or --use-repo XYZ "
fi

if [[ ! -z "$CGATFLOW_REPO" ]] ; then
    CGATFLOW_REPO=$("$READLINK" -f "$CGATFLOW_REPO")
    if [[ ! -e "$CGATFLOW_REPO/setup.py" ]] ; then
	report_error "No setup.py present in $CGATFLOW_REPO"
    fi
fi

# sanity check: make sure there is space available in the destination folder (20 GB) in 512-byte blocks
[[ -z ${TRAVIS_INSTALL} ]] && \
    mkdir -p ${CGAT_HOME} && \
    [[ `df -P ${CGAT_HOME} | awk '/\// {print $4}'` -lt 41943040 ]] && \
    report_error " Not enough disk space available on the installation folder: $CGAT_HOME"


[[ -z ${TRAVIS_BRANCH} ]] && TRAVIS_BRANCH=${CGATFLOW_BRANCH}

if [[ $INSTALL_DEVEL ]] ; then
    conda_install
fi

if [[ $INSTALL_PIPELINE_DEPENDENCIES -eq 1 ]]; then
    install_pipeline_deps
fi

if [[ $RUN_TESTS ]] ; then
    conda_test
fi

if [[ $INSTALL_UPDATE ]] ; then
    conda_update
fi

if [[ $UNINSTALL ]] ; then
    uninstall
fi
