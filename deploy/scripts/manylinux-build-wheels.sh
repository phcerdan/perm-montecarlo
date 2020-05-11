#!/usr/bin/env bash

# -----------------------------------------------------------------------
# These variables are set in common script:
#
ARCH=""
PYBINARIES=""

script_dir=$(cd $(dirname $0) || exit 1; pwd)
source "${script_dir}/manylinux-build-common.sh"
# -----------------------------------------------------------------------

deploy_dir=${script_dir}/..
pushd ${deploy_dir}
# Compile wheels re-using standalone project and archive cache
for PYBIN in "${PYBINARIES[@]}"; do
    PYTHON_EXECUTABLE=${PYBIN}/python
    PYTHON_INCLUDE_DIR=$( find -L ${PYBIN}/../include/ -name Python.h -exec dirname {} \; )
    PYTHON_INCLUDE_DIRS=${PYTHON_INCLUDE_DIR}

    echo ""
    echo "PYTHON_EXECUTABLE:${PYTHON_EXECUTABLE}"
    echo "PYTHON_INCLUDE_DIR:${PYTHON_INCLUDE_DIR}"

    # Remove when scikit-build includes cmake_target PR:
    # https://github.com/scikit-build/scikit-build/pull/477
    ${PYBIN}/python -m pip uninstall scikit-build -y
    ${PYBIN}/python -m pip install -r requirements-deploy.txt

    ${PYBIN}/python setup.py bdist_wheel --build-type Release -G Ninja -- \
      -DCMAKE_CXX_COMPILER_TARGET:STRING=$(uname -p)-linux-gnu \
      -DPERM_BUILD_TESTING:BOOL=OFF \
      -DPERM_WRAP_PYTHON:BOOL=ON \
      -DPYTHON_EXECUTABLE:FILEPATH=${PYTHON_EXECUTABLE} \
      -DPYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR} \
    || exit 1
    # ${PYBIN}/python setup.py clean
done

# auditwheel will bundle shared libraries in the wheel,
# but they have to be found first using LD_LIBRARY_PATH
export DEPENDENCIES_LD_LIBRARY_PATH=""
# This step will fixup the wheel switching from 'linux' to 'manylinux2014' tag and include third party libraries
for whl in dist/*linux_$(uname -p).whl; do
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${DEPENDENCIES_LD_LIBRARY_PATH} auditwheel repair ${whl} -w /work/dist/
    rm ${whl}
done
popd
