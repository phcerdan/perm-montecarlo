#!/usr/bin/env bash

# Assumes that python points to the python you want to use
# This should be managed from other script, or from CI/CD pipelines

script_dir=$(cd $(dirname $0) || exit 1; pwd)

echo ""
echo "$(python --version)"

python -m pip install -r $script_dir/requirements-deploy.txt
python -m pip install delocate


pushd ${script_dir}
python setup.py bdist_wheel --build-type Release -G Ninja -- \
  -DCMAKE_OSX_DEPLOYMENT_TARGET:STRING=10.9 \
  -DCMAKE_OSX_ARCHITECTURES:STRING=x86_64 \
  -DPERM_BUILD_TESTING:BOOL=OFF \
  -DPERM_WRAP_PYTHON:BOOL=ON \
  || exit 1
  # ${PYBIN}/python setup.py clean

export DEPENDENCIES_LD_LIBRARY_PATH=""

DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${DEPENDENCIES_LD_LIBRARY_PATH} delocate-listdeps $PWD/dist/*.whl # lists library dependencies
DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${DEPENDENCIES_LD_LIBRARY_PATH} delocate-wheel $PWD/dist/*.whl # copies library dependencies into wheel
popd
