#!/usr/bin/env bash

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

# Bash traps
# http://aplawrence.com/Basics/trapping_errors.html
# https://stelfox.net/blog/2013/11/fail-fast-in-bash-scripts/

set -o errtrace

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
   echo "# install-CGAT-tools.sh log | `hostname` | `date` | $1 "
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

if [[ $TRAVIS_INSTALL ]] ; then

   CGAT_HOME=$TRAVIS_BUILD_DIR
   CONDA_INSTALL_TYPE_PIPELINES="pipelines-nosetests.yml"
   CONDA_INSTALL_TYPE_APPS="apps-nosetests.yml"
   CONDA_INSTALL_TYPE_CORE="core-production.yml"

elif [[ $JENKINS_INSTALL ]] ; then

   CGAT_HOME=$WORKSPACE
   CONDA_INSTALL_TYPE_PIPELINES="pipelines-devel.yml"
   CONDA_INSTALL_TYPE_APPS="apps-devel.yml"
   CONDA_INSTALL_TYPE_CORE="core-production.yml"

else

   if [[ -z $CGAT_HOME  ]] ; then
      CGAT_HOME=$HOME/cgat-install
   fi

   if [[ $INSTALL_PRODUCTION ]] ; then
      CONDA_INSTALL_TYPE_PIPELINES="pipelines-production.yml"
      CONDA_INSTALL_TYPE_APPS="apps-production.yml"
      CONDA_INSTALL_TYPE_CORE="core-production.yml"
   elif [[ $INSTALL_DEVEL ]] ; then
      CONDA_INSTALL_TYPE_PIPELINES="pipelines-devel.yml"
      CONDA_INSTALL_TYPE_APPS="apps-devel.yml"
      CONDA_INSTALL_TYPE_CORE="core-devel.yml"
   elif [[ $INSTALL_TEST ]] || [[ $INSTALL_UPDATE ]] ; then
      if [[ -d $CGAT_HOME/conda-install ]] ; then
         AUX=`find $CGAT_HOME/conda-install/envs/cgat-* -maxdepth 0`
         CONDA_INSTALL_TYPE_PIPELINES=`basename $AUX`
      else
         echo
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

fi # if travis install

# set installation folder
CONDA_INSTALL_DIR=$CGAT_HOME/conda-install

# set conda environment name
[[ ${CONDA_INSTALL_ENV} ]] || CONDA_INSTALL_ENV="cgat-f"

} # get_cgat_env


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
[[ ! $INSTALL_TEST ]] && echo " PIPELINES_BRANCH: "$PIPELINES_BRANCH
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
source ${CONDA_INSTALL_DIR}/bin/activate
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

[[ -z ${TRAVIS_BRANCH} ]] && TRAVIS_BRANCH=${PIPELINES_BRANCH}

curl -o env-core.yml -O https://raw.githubusercontent.com/cgat-developers/cgat-core/${CORE_BRANCH}/conda/environments/${CONDA_INSTALL_TYPE_CORE}

curl -o env-apps.yml -O https://raw.githubusercontent.com/cgat-developers/cgat-apps/${APPS_BRANCH}/conda/environments/${CONDA_INSTALL_TYPE_APPS}

curl -o env-pipelines.yml -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${TRAVIS_BRANCH}/conda/environments/${CONDA_INSTALL_TYPE_PIPELINES}

[[ ${CLUSTER} -eq 0 ]] && sed -i'' -e '/drmaa/d' env-pipelines.yml

conda env create --quiet --name ${CONDA_INSTALL_ENV} --file env-pipelines.yml
conda env update --quiet --name ${CONDA_INSTALL_ENV} --file env-apps.yml
conda env update --quiet --name ${CONDA_INSTALL_ENV} --file env-core.yml

# activate cgat environment
source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_ENV

log "installing CGAT code into conda environment"
# if installation is 'devel' (outside of travis), checkout latest version from github
if [[ -z ${TRAVIS_INSTALL} ]] ; then

   DEV_RESULT=0

   if [[ $INSTALL_DEVEL ]] || [[ $INSTALL_PRODUCTION ]] || [[ $JENKINS_INSTALL ]] ; then

      # install extra deps
      install_extra_deps

      # install Python 2 deps
      install_py2_deps

      # make sure you are in the CGAT_HOME folder
      cd $CGAT_HOME

      # download the code out of jenkins
      if [[ -z ${JENKINS_INSTALL} ]] ; then

         if [[ $CODE_DOWNLOAD_TYPE -eq 0 ]] ; then
            # get the latest version from Git Hub in zip format
            curl -LOk https://github.com/cgat-developers/cgat-flow/archive/$PIPELINES_BRANCH.zip
            unzip $PIPELINES_BRANCH.zip
            rm $PIPELINES_BRANCH.zip
            if [[ ${RELEASE} ]] ; then
               NEW_NAME=`echo $PIPELINES_BRANCH | sed 's/^v//g'`
               mv cgat-flow-$NEW_NAME/ cgat-flow/
            else
               mv cgat-flow-$PIPELINES_BRANCH/ cgat-flow/
            fi
         elif [[ $CODE_DOWNLOAD_TYPE -eq 1 ]] ; then
            # get latest version from Git Hub with git clone
            git clone --branch=$PIPELINES_BRANCH https://github.com/cgat-developers/cgat-flow.git
         elif [[ $CODE_DOWNLOAD_TYPE -eq 2 ]] ; then
            # get latest version from Git Hub with git clone
            git clone --branch=$PIPELINES_BRANCH git@github.com:cgat-developers/cgat-flow.git
         else
            report_error " Unknown download type for CGAT code... "
         fi

         # make sure you are in the CGAT_HOME/cgat-flow folder
         cd $CGAT_HOME/cgat-flow

      else

         # in jenkins the code is in CGAT_HOME/CGATPipelines
         cd $CGAT_HOME/CGATPipelines

      fi

      # Set up other environment variables
      setup_env_vars

      # brute force: modify console_scripts variable/entry point for cgat command
      sed -i'' -e 's/CGATScripts/scripts/g' setup.py

      # install cgat-core
      install_cgat_core

      # install cgat-apps
      install_cgat_apps

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

      fi # if-$?

      # revert setup.py if downloaded with git
      [[ $CODE_DOWNLOAD_TYPE -ge 1 ]] && git checkout -- setup.py

   fi # if INSTALL_DEVEL

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
      echo " $ source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_ENV"
      echo
      echo " To deactivate the environment, use:"
      echo " $ source deactivate"
      echo
   fi # if-$ conda create

fi # if travis install

} # conda install


# install extra dependencies
install_extra_deps() {

log "install extra deps"

curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${TRAVIS_BRANCH}/conda/environments/pipelines-extra.yml
curl -O https://raw.githubusercontent.com/cgat-developers/cgat-apps/${APPS_BRANCH}/conda/environments/apps-extra.yml

conda env update --quiet --name ${CONDA_INSTALL_ENV} --file pipelines-extra.yml
conda env update --quiet --name ${CONDA_INSTALL_ENV} --file apps-extra.yml

if [[ ${INSTALL_IDE} -eq 1 ]] ; then
   curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${TRAVIS_BRANCH}/conda/environments/pipelines-ide.yml
   conda env update --quiet --name ${CONDA_INSTALL_ENV} --file pipelines-ide.yml
fi

}


# install Python 2 dependencies in a different conda environment
install_py2_deps() {

log "install Python 2 deps"

curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${TRAVIS_BRANCH}/conda/environments/pipelines-macs2.yml

conda env update --quiet --file pipelines-macs2.yml

curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${TRAVIS_BRANCH}/conda/environments/pipelines-tophat2.yml

conda env update --quiet --file pipelines-tophat2.yml


curl -O https://raw.githubusercontent.com/cgat-developers/cgat-flow/${TRAVIS_BRANCH}/conda/environments/pipeline-peakcalling-sicer.yml

conda env update --quiet --file pipeline-peakcalling-sicer.yml

}


# helper function to install cgat-apps
install_cgat_apps() {

log "install cgat apps"

OLDWD=`pwd`
cd $CGAT_HOME

if [[ $JENKINS_INSTALL ]] ; then

   cd cgat/

else

   if [[ $CODE_DOWNLOAD_TYPE -eq 0 ]] ; then
      # get the latest version from Git Hub in zip format
      curl -LOk https://github.com/cgat-developers/cgat-apps/archive/$APPS_BRANCH.zip
      unzip $APPS_BRANCH.zip
      rm $APPS_BRANCH.zip
      if [[ ${RELEASE} ]] ; then
         NEW_NAME=`echo $APPS_BRANCH | sed 's/^v//g'`
         mv cgat-apps-$NEW_NAME/ cgat-apps/
      else
         mv cgat-apps-$APPS_BRANCH/ cgat-apps/
      fi
   elif [[ $CODE_DOWNLOAD_TYPE -eq 1 ]] ; then
      # get latest version from Git Hub with git clone
      git clone --branch=$APPS_BRANCH https://github.com/cgat-developers/cgat-apps.git
   elif [[ $CODE_DOWNLOAD_TYPE -eq 2 ]] ; then
      # get latest version from Git Hub with git clone
      git clone --branch=$APPS_BRANCH git@github.com:cgat-developers/cgat-apps.git
   else
      report_error " Unknown download type for CGAT code... "
   fi

   cd cgat-apps/

fi

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
   curl -LOk https://github.com/cgat-developers/cgat-core/archive/$CORE_BRANCH.zip
   unzip $CORE_BRANCH.zip
   rm $CORE_BRANCH.zip
   if [[ ${RELEASE} ]] ; then
      NEW_NAME=`echo $CORE_BRANCH | sed 's/^v//g'`
      mv cgat-core-$NEW_NAME/ cgat-core/
   else
      mv cgat-core-$CORE_BRANCH/ cgat-core/
   fi
elif [[ $CODE_DOWNLOAD_TYPE -eq 1 ]] ; then
   # get latest version from Git Hub with git clone
   git clone --branch=$CORE_BRANCH https://github.com/cgat-developers/cgat-core.git
elif [[ $CODE_DOWNLOAD_TYPE -eq 2 ]] ; then
   # get latest version from Git Hub with git clone
   git clone --branch=$CORE_BRANCH git@github.com:cgat-developers/cgat-core.git
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

setup_env_vars

# setup environment and run tests
if [[ $TRAVIS_INSTALL ]] || [[ $JENKINS_INSTALL ]] ; then

   # enable Conda env
   log "activating CGAT conda environment"
   source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_ENV

   # show conda environment used for testing
   conda env export

   # install cgat-core
   install_cgat_core

   # install cgat-apps
   install_cgat_apps

   # python preparation
   log "install CGAT code into conda environment"
   [[ $TRAVIS_INSTALL  ]] && cd $CGAT_HOME
   [[ $JENKINS_INSTALL ]] && cd $CGAT_HOME/CGATPipelines

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
      # prepare environment
      source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_ENV

      # make sure you are in the CGAT_HOME/cgat-flow folder
      cd $CGAT_HOME/cgat-flow

      # Python preparation
      sed -i'' -e 's/CGATScripts/scripts/g' setup.py
      sed -i'' -e '/REPO_REQUIREMENT/,/pass/d' setup.py
      sed -i'' -e '/# dependencies/,/dependency_links=dependency_links,/d' setup.py
      python setup.py develop
      OUTPUT_DIR=`pwd`

      # run tests
      /usr/bin/time -o test_import.time -v nosetests -v tests/test_import.py >& test_import.out
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

      /usr/bin/time -o test_scripts.time -v nosetests -v tests/test_scripts.py >& test_scripts.out
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

source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_ENV
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
      if [[ "$PIPELINES_BRANCH" != "master" ]] || [[ $APPS_BRANCH != "master" ]] ; then
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
   curl --output /dev/null --silent --head --fail https://raw.githubusercontent.com/cgat-developers/cgat-core/${CORE_BRANCH}/README.rst || RELEASE_TEST=$?
   if [[ ${RELEASE_TEST} -ne 0 ]] ; then
      echo
      echo " The branch provided for cgat-core does not exist: ${CORE_BRANCH}"
      echo
      echo " Please have a look at valid branches here: "
      echo " https://github.com/cgat-developers/cgat-core/branches"
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


# function to display help message
help_message() {
echo
echo " This script uses Conda to install cgat-flow. To proceed, please type:"
echo " ./install-CGAT-tools.sh --devel [--location </full/path/to/folder/without/trailing/slash>]"
echo
echo " The default install folder will be: $HOME/cgat-install"
echo
echo " It is also possible to install/test a specific branch of the code on GitHub:"
echo " ./install-CGAT-tools.sh --devel --pipelines-branch <branch> --scripts-branch <branch> --core-branch <branch>"
echo
echo " This will create an isolated Conda environment with both the pipelines and the scripts from:"
echo " https://github.com/cgat-developers/cgat-apps"
echo
echo " The default name of the newly created conda environment is cgat-f, but you can change it with:"
echo " --env-name name"
echo
echo " The code is downloaded in zip format by default. If you want to get a git clone, use:"
echo " --git # for an HTTPS clone"
echo " --git-ssh # for a SSH clone (you need to be a CGATOXford contributor on GitHub to do this)"
echo
echo " The pipelines are intended to run on a cluster using the DRMAA API. If that's not your case, please use:"
echo " --no-cluster"
echo
echo " If you want to download and install IDEs like Spyder or RStudio with this installation, please use:"
echo " --ide"
echo
echo " To test the installation:"
echo " ./install-CGAT-tools.sh --test [--location </full/path/to/folder/without/trailing/slash>]"
echo
echo " To update the Conda packages:"
echo " ./install-CGAT-tools.sh --update [--location </full/path/to/folder/without/trailing/slash>]"
echo
echo " To uninstall the CGAT code:"
echo " ./install-CGAT-tools.sh --uninstall [--location </full/path/to/folder/without/trailing/slash>]"
echo
echo " Please submit any issues via Git Hub:"
echo " https://github.com/cgat-developers/cgat-flow/issues"
echo
exit 1
} # help_message

# the script starts here

if [[ $# -eq 0 ]] ; then

   help_message

fi

# travis execution
TRAVIS_INSTALL=
# jenkins testing
JENKINS_INSTALL=
# conda installation type
INSTALL_PRODUCTION=
INSTALL_DEVEL=
# test current installation
INSTALL_TEST=
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
PIPELINES_BRANCH="master"
APPS_BRANCH="master"
CORE_BRANCH="master"
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

while [[ $# -gt 0 ]]
do
key="$1"

case $key in

    --help)
    help_message
    ;;

    --travis)
    TRAVIS_INSTALL=1
    shift # past argument
    ;;

    --jenkins)
    JENKINS_INSTALL=1
    shift # past argument
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

    --production)
    INSTALL_PRODUCTION=1
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

    --location)
    CGAT_HOME="$2"
    shift 2
    ;;

    --pipelines-branch)
    PIPELINES_BRANCH="$2"
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
    PIPELINES_BRANCH="$2"
    APPS_BRANCH="$2"
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

# sanity check 1: don't mix production and development installs
if [[ $INSTALL_PRODUCTION ]] && [[ $INSTALL_DEVEL ]] ; then

   report_error " Incorrect input arguments: mixing --production and --devel is not permitted. "

fi

# sanity check 2: make sure one installation option is selected
if [[ -z $INSTALL_TEST ]] && \
   [[ -z $INSTALL_PRODUCTION ]] && \
   [[ -z $INSTALL_DEVEL ]] && \
   [[ -z $TRAVIS_INSTALL ]] && \
   [[ -z $JENKINS_INSTALL ]] ; then

   report_error " You need to select either --devel or --production. "

fi

# sanity check 3: make sure there is space available in the destination folder (20 GB) in 512-byte blocks
[[ -z ${TRAVIS_INSTALL} ]] && \
mkdir -p ${CGAT_HOME} && \
[[ `df -P ${CGAT_HOME} | awk '/\// {print $4}'` -lt 41943040 ]] && \
   report_error " Not enought disk space available on the installation folder: "$CGAT_HOME

# perform actions according to the input parameters processed
if [[ $TRAVIS_INSTALL ]] || [[ $JENKINS_INSTALL ]] ; then

  conda_install
  conda_test

else 

  if [[ $OS_PKGS ]] ; then
     install_os_packages
  fi

  if [[ $INSTALL_PRODUCTION ]] || [[ $INSTALL_DEVEL ]] ; then
     conda_install
  fi

  if [[ $INSTALL_TEST ]] ; then
     conda_test
  fi

  if [[ $INSTALL_UPDATE ]] ; then
     conda_update
  fi

  if [[ $UNINSTALL ]] ; then
     uninstall
  fi

fi # if-variables

