#!/bin/bash -uxe

# Quick & dirty installer. This is not a replacement for install-CGAT-tools.sh.
#
# This script will create a cgat-devel environment in your conda
# installation. It assumes that:
#
# a. you have installed conda and activated it
# b. you have downloaded the three CGAT repositories cgat-core, cgat-apps and cgat-flow
# c. you have set the variable REPO_DIR to the location of the repositories

: ${REPO_DIR:=/ifs/devel/andreas/cgatnew}

PROJECTS="cgat-core cgat-apps cgat-flow"

conda create -y -n cgat-devel python=3.6
source activate cgat-devel

for x in $PROJECTS; do
	echo "installing dependencies for $x"
       fn="$REPO_DIR/$x/conda_requires.txt"
	conda install -y `cat "$fn" | grep -v "#" | xargs`;
done

echo "install miscellaneous non-conda dependencies"

pip install quicksect

echo "activating projects"

for x in $PROJECTS; do
	cd "$REPO_DIR/$x" && python setup.py develop
done

conda create -y -n macs2 python=2.7 macs2
conda create -y -n sicer python=2.6 sicer
