#!/usr/bin/env bash

# from https://conda.io/projects/conda/en/latest/user-guide/tasks/use-conda-with-travis-ci.html (March 1 2019)
if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    MINICONDA_OS=Linux;
elif [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    MINICONDA_OS=MacOSX;
fi
if [[ "$TRAVIS_PYTHON_VERSION" == "2.7" ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda2-latest-${MINICONDA_OS}-x86_64.sh -O miniconda.sh;
else
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-${MINICONDA_OS}-x86_64.sh -O miniconda.sh;
fi
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda
# Useful for debugging any issues with conda
conda info -a
