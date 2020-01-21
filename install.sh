#!/usr/bin/env bash

# cgat-flow-conda: Integration test for pipelines running against the
# latest conda packages of cgat-core and cgat-apps. The idea is that this
# is the situation of a user just trying to run or extend our pipelines
# without getting into cgat-core or cgat-apps internals.

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
    echo "# install.sh log | `hostname` | `date` | $1 "
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

    # minimum pinning pipeline
    CONDA_INSTALL_TYPE_PIPELINES="cgat-flow.yml"

    if [[ ! $CONDA_INSTALL  ]] ; then
	# set installation folder
	CONDA_INSTALL_DIR=$CGAT_HOME/conda-install
    else
	CONDA_INSTALL_DIR=${CONDA_PREFIX}
    fi
    # set conda environment name
    [[ ${CONDA_INSTALL_ENV} ]] || CONDA_INSTALL_ENV="cgat-flow"

} # get_cgat_env


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
        # activate cgat-flow env
        conda activate ${CONDA_INSTALL_ENV}
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
    echo " CONDA_INSTALL_DIR: "$CONDA_INSTALL_DIR
    echo " CONDA_INSTALL_TYPE_CORE:"$CONDA_INSTALL_TYPE_CORE
    echo " CONDA_INSTALL_TYPE_APPS: "$CONDA_INSTALL_TYPE_APPS
    echo " CONDA_INSTALL_TYPE_PIPELINES: "$CONDA_INSTALL_TYPE_PIPELINES
    echo " CONDA_INSTALL_ENV: "$CONDA_INSTALL_ENV
    echo " PYTHONPATH: "$PYTHONPATH
    [[ ! $INSTALL_TEST ]] && echo " BRANCH: "$BRANCH
    [[ ! $INSTALL_TEST ]] && echo " APPS_BRANCH: "$APPS_BRANCH
    [[ ! $INSTALL_TEST ]] && echo " CORE_BRANCH: "$CORE_BRANCH
    [[ ! $INSTALL_TEST ]] && echo " RELEASE: "$RELEASE
    echo " CODE_DOWNLOAD_TYPE: "$CODE_DOWNLOAD_TYPE
    echo " INSTALL_IDE: "$INSTALL_IDE
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


# install and activate miniconda
miniconda_install() {

    log "installing miniconda"
    
    # get environment variables: CGAT_HOME, CONDA_INSTALL_DIR, CONDA_INSTALL_TYPE_PIPELINES
    get_cgat_env

    # select Miniconda bootstrap script depending on Operating System
    MINICONDA=

    if [[ `uname` == "Linux" ]] ; then

	# Conda 4.4 breaks everything again!
	# Conda 4.5 looks better
	MINICONDA="Miniconda3-latest-Linux-x86_64.sh"

    elif [[ `uname` == "Darwin" ]] ; then

	# Conda 4.4 breaks everything again!
	# Conda 4.5 looks better
	MINICONDA="Miniconda3-latest-MacOSX-x86_64.sh"

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
    # conda install --quiet --yes 'conda=4.3.33'
    # conda update --all --yes
    conda info -a
}


# proceed with conda installation
conda_install() {

    log "installing conda"

    detect_cgat_installation

    get_cgat_env

    mkdir -p $CGAT_HOME
    cd $CGAT_HOME

    # activate cgat environment
    is_env_enabled

    log "romoving old environment ${CONDA_INSTALL_ENV} if it exists"
    conda env remove -y -n ${CONDA_INSTALL_ENV} >& /dev/null || echo "Not removing environment ${CONDA_INSTALL_ENV} because it does not exist"

    log "installing cgat-core, cgat-apps and its dependencies"
    conda create -y --name ${CONDA_INSTALL_ENV} -c conda-forge -c bioconda cgatcore cgat-apps

    log "curl -o env-cgat-flow.yml -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${BRANCH}/conda/environments/${CONDA_INSTALL_TYPE_PIPELINES}"
    curl -o env-cgat-flow.yml -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${BRANCH}/conda/environments/${CONDA_INSTALL_TYPE_PIPELINES}
    
    [[ ${CLUSTER} -eq 0 ]] && sed -i'' -e '/drmaa/d' env-cgat-flow.yml

    conda env update --name ${CONDA_INSTALL_ENV} --file env-cgat-flow.yml

    log "curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${BRANCH}/conda/environments/cgat-flow-pipelines."
    curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${BRANCH}/conda/environments/cgat-flow-pipelines.yml

    conda env update --name ${CONDA_INSTALL_ENV} --file cgat-flow-pipelines.yml

    # install extra deps
    install_extra_deps

    # install Python 2 deps
    install_py2_deps
}

code_install() {
	
    detect_cgat_installation

    get_cgat_env

    # activate the environment
    is_env_enabled

    OLDPWD=`pwd`

    log "installing cgat-flow code into conda environment"
    if [[ ${CODE_INSTALL} || ${FULL_INSTALL} ]] ; then

	cd ${CGAT_HOME}

        if [[ $CODE_DOWNLOAD_TYPE -eq 0 ]] ; then
	    # get the latest version from Git Hub in zip format
	    curl -LOk https://github.com/cgat-developers/cgat-flow/archive/$BRANCH.zip
	    unzip $BRANCH.zip
	    rm $BRANCH.zip
	    if [[ ${RELEASE} ]] ; then
	        NEW_NAME=`echo $BRANCH | sed 's/^v//g'`
		mv cgat-flow-$NEW_NAME/ cgat-flow/
	    else
	        mv cgat-flow-$BRANCH/ cgat-flow/
            fi
        elif [[ $CODE_DOWNLOAD_TYPE -eq 1 ]] ; then
	    # get latest version from Git Hub with git clone
	    git clone --branch=$BRANCH https://github.com/cgat-developers/cgat-flow.git
        elif [[ $CODE_DOWNLOAD_TYPE -eq 2 ]] ; then
	    # get latest version from Git Hub with git clone
	    git clone --branch=$BRANCH git@github.com:cgat-developers/cgat-flow.git
        else
	    report_error " Unknown download type for CGAT code... "
	fi

        # make sure you are in the CGAT_HOME/cgat-flow folder
	cd ${CGAT_HOME}/cgat-flow

    else
	cd ${REPO_DIR}
    fi

    # Set up other environment variables
    setup_env_vars
    
    # Python preparation
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
	
    else
	clear
	echo 
	echo " The code was successfully installed!"
	echo
	echo " To activate the CGAT environment type: "
	echo " $ source $CONDA_INSTALL_DIR/etc/profile.d/conda.sh"
	echo " $ conda activate base"
	echo " $ conda activate $CONDA_INSTALL_ENV"
	[[ $INSTALL_PRODUCTION ]] && echo " cgatflow --help"
	echo
	echo " To deactivate the environment, use:"
	echo " $ conda deactivate"
	echo
    fi # if-$ conda create

    # revert setup.py if downloaded with git
    [[ $CODE_DOWNLOAD_TYPE -ge 1 ]] && git checkout -- setup.py

    conda env export > environment.yml

    # go back to old working directory
    cd $OLDPWD

} # code_install


# install extra dependencies
install_extra_deps() {

    log "installing extra deps"

    curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${BRANCH}/conda/environments/pipelines-extra.yml

    conda env update --quiet --name ${CONDA_INSTALL_ENV} --file pipelines-extra.yml

    if [[ ${INSTALL_IDE} -eq 1 ]] ; then
	curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${BRANCH}/conda/environments/pipelines-ide.yml
	conda env update --quiet --name ${CONDA_INSTALL_ENV} --file pipelines-ide.yml
    fi

}


# install Python 2 dependencies in a different conda environment
install_py2_deps() {

    log "installing Python 2 deps: macs2"
    curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${BRANCH}/conda/environments/pipelines-macs2.yml

    conda env update --quiet --file pipelines-macs2.yml

    log "installing Python 2 deps: tophat2"
    curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${BRANCH}/conda/environments/pipelines-tophat2.yml

    conda env update --quiet --file pipelines-tophat2.yml

    log "installing Python 2 deps: sicer"
    curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${BRANCH}/conda/environments/pipeline-peakcalling-sicer.yml

    conda env update --quiet --file pipeline-peakcalling-sicer.yml

    log "installing Python 2 deps: splicing"
    curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${BRANCH}/conda/environments/pipelines-splicing.yml

    conda env update --quiet --file pipelines-splicing.yml

}

code_test() {

    detect_cgat_installation

    get_cgat_env

    # activate the environment
    is_env_enabled

    OLDPWD=`pwd`

    log "starting Python tests"

    conda env list
    conda env export

    cd $CGAT_HOME/cgat-flow

    nosetests -v tests/test_import.py
    if [[ $? -eq 0 ]] ; then
        echo
	echo " test_imports.py passed sucessfully!"
	echo
    else
	echo
	echo " test_script.py failed. "
	echo
	print_env_vars
    fi

    nosetests -v tests/test_scripts.py
    if [[ $? -eq 0 ]] ; then
        echo
	echo " test_imports.py passed sucessfully!"
	echo
    else
	echo
	echo " test_script.py failed. "
	echo
	print_env_vars
    fi


    # go back to old working directory
    cd $OLDPWD

} # code_test


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
    # source deactivate >& /dev/null || true
    # source deactivate >& /dev/null || true
    # unset -f conda || true
    # unset PYTHONPATH || true

    # Next actions disabled. Please see:
    # https://github.com/cgat-developers/cgat-core/issues/44
    #module purge >& /dev/null || true
    #mymodule purge >& /dev/null || true
    set -e
}


# function to display help message
help_message() {
    echo
    echo " This script uses Conda to install cgat-flow. To do a full install, please type:"
    echo " ./install.sh --full-install --install-dir </full/path/to/folder/without/trailing/slash>"
    echo
    echo " You can test your installation with:"
    echo "./install.sh --test --install-dir </full/path/to/install/dir>"
    echo
    exit 1
} # help_message


cleanup_env
test_compilers

# jenkins testing
JENKINS_INSTALL=
# conda install - requires an activate environment
FULL_INSTALL=
# conda install - requires an activate environment
CONDA_INSTALL=
# conda install - requires an activate environment
CODE_INSTALL=
# Location of repository. By default location of install script
REPO_DIR="$(readlink -f "$(dirname "${BASH_SOURCE[0]}")")"
# test current installation
INSTALL_TEST=
# update current installation
INSTALL_UPDATE=
# whether to clone repository
CLONE_REPO=
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
BRANCH="master"
# type of installation
CONDA_INSTALL_TYPE_PIPELINES=
CONDA_INSTALL_TYPE_APPS=
CONDA_INSTALL_TYPE_CORE=
# rename conda environment
CONDA_INSTALL_ENV=
# install additional IDEs?
# 0 = no (default)
# 1 = yes
INSTALL_IDE=0
# Use cluster?
# 0 = no
# 1 = yes (default)
CLUSTER=1
# Install a released version?
RELEASE=

# parse input parameters
# https://stackoverflow.com/questions/402377/using-getopts-in-bash-shell-script-to-get-long-and-short-command-line-options
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash

# the script starts here

while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in

	--help)
	    help_message
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

	--clone)
	    CLONE_REPO=1
	    shift # past argument
	    ;;

	--full-install)
	    FULL_INSTALL=1
	    shift
	    ;;

	--conda-install)
	    CONDA_INSTALL=1
	    shift
	    ;;

	--code-install)
	    CODE_INSTALL=1
	    shift
	    ;;
	
	--devel)
	    INSTALL_DEVEL=1
	    shift
	    ;;

	--test)
	    INSTALL_TEST=1
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

	--install-dir)
	    CGAT_HOME="$2"
	    shift 2
	    ;;

	--pipelines-branch)
	    BRANCH="$2"
	    test_mix_branch_release
	    shift 2
	    ;;

	--scripts-branch)
	    APPS_BRANCH="$2"
	    test_mix_branch_release
	    shift 2
	    ;;

	--core-branch)
	    CORE_BRANCH="$2"
	    test_core_branch
	    shift 2
	    ;;

	--ide)
	    INSTALL_IDE=1
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
	    BRANCH="$2"
	    shift 2
	    ;;

	--env-name)
	    CONDA_INSTALL_ENV="$2"
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

if [[ !($CODE_INSTALL || $CONDA_INSTALL || $FULL_INSTALL || $INSTALL_TEST) ]] ; then
    help_message
fi

# sanity check: make sure there is space available in the destination folder (20 GB) in 512-byte blocks
[[ -z ${TRAVIS_INSTALL} ]] && \
    mkdir -p ${CGAT_HOME} && \
    [[ `df -P ${CGAT_HOME} | awk '/\// {print $4}'` -lt 41943040 ]] && \
    report_error " Not enough disk space available on the installation folder: "$CGAT_HOME

if [[ $FULL_INSTALL ]] ; then
    miniconda_install
fi
    
if [[ $CONDA_INSTALL || $FULL_INSTALL ]] ; then
    conda_install
fi

if [[ $CODE_INSTALL || $CONDA_INSTALL || $FULL_INSTALL ]] ; then
    code_install
fi

if [[ $INSTALL_TEST ]] ; then
    code_test
fi

if [[ $UNINSTALL ]] ; then
    uninstall
fi
