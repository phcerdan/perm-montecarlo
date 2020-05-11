from __future__ import print_function
from os import sys, path

try:
    from skbuild import setup
except ImportError:
    print('scikit-build is required to build from source.', file=sys.stderr)
    print('Please run:', file=sys.stderr)
    print('', file=sys.stderr)
    print('  python -m pip install scikit-build')
    sys.exit(1)

sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
from perm_montecarlo_version import get_versions

long_description= r'perm_montecarlo is an open-source, cross-platform c++ library ' \
                   'performing montecarlo simulations of self-avoiding walks.\n' \
                   'Includes rosenbluth sampling, and PERM (Prune and Enriched Rosenbluth Method. ' \
                   'It allows to control the amount of nearest neighbors, in 2D and 3D.'

setup(
    name='perm_montecarlo',
    version=get_versions()['package-version'],
    author='Pablo Hernandez-Cerdan',
    author_email='pablo.hernandez.cerdan@outlook.com',
    packages=['perm_montecarlo'],
    download_url=r'https://github.com/phcerdan/perm-montecarlo',
    description=r'PERM: Prune and Enrichment Rosenbluth Method for self-avoiding walks',
    long_description=long_description,
    classifiers=[
        "License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)",
        "Programming Language :: Python",
        "Programming Language :: C++",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Software Development :: Libraries",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Operating System :: MacOS"
        ],
    license='MPL2',
    keywords='PERM montecarlo SAW polymer simulation',
    url=r'https://github.com/phcerdan/perm-montecarlo',
    install_requires=[],
    cmake_source_dir='../',
    cmake_args=[
        '-DBUILD_SHARED_LIBS:BOOL=FALSE',
        '-DPERM_BUILD_TESTING:BOOL=FALSE',
        '-DPERM_WRAP_PYTHON:BOOL=TRUE',
    ],
    cmake_install_target="install-runtime"
    )

