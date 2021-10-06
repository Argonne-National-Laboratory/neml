# This file lightly adapted from the example CMake build provided in the pybind11 distribution
# located at https://github.com/pybind/cmake_example/blob/master/setup.py
# That file is released under an open source license, contained in the pybind11 subdirectory
# of this package.

import os
import re
import sys
import platform
import subprocess

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding = 'utf-8') as f:
  long_description = f.read()

class CMakeExtension(Extension):
  def __init__(self, name, sourcedir=''):
      Extension.__init__(self, name, sources=[])
      self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
  def run(self):
    try:
      out = subprocess.check_output(['cmake', '--version'])
    except OSError:
      raise RuntimeError("CMake must be installed to build the following extensions: " +
                          ", ".join(e.name for e in self.extensions))

    if platform.system() == "Windows":
      cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
      if cmake_version < '3.1.0':
        raise RuntimeError("CMake >= 3.1.0 is required on Windows")

    for ext in self.extensions:
        self.build_extension(ext)

  def build_extension(self, ext):
    extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
    # required for auto-detection of auxiliary "native" libs
    if not extdir.endswith(os.path.sep):
      extdir += os.path.sep

    cmake_args = [
                  '-DWRAP_PYTHON=ON',
                  '-DUSE_OPENMP=OFF',
                  '-DBUILD_UTILS=OFF',
                  '-DCMAKE_INSTALL_PREFIX=' + extdir,
                  '-DPYTHON_PACKAGE=ON',
                  '-DPYTHON_EXECUTABLE={}'.format(sys.executable)]

    cfg = 'Debug' if self.debug else 'Release'
    build_args = ['--config', cfg]
    cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]

    if platform.system() == "Windows":
      cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
      if sys.maxsize > 2**32:
          cmake_args += ['-A', 'x64']
      build_args += ['--', '/m']
    else:
      build_args += ['--', '-j2']

    env = os.environ.copy()
    env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                          self.distribution.get_version())
    if not os.path.exists(self.build_temp):
        os.makedirs(self.build_temp)
  
    subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
    subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)
    subprocess.check_call(['cmake', '--install', '.'], cwd = self.build_temp)

setup (
    # Name of the project
    name = 'neml',
    # Version
    version = '1.4.0',
    # One line-description
    description = "Nuclear Engineering Material model Library",
    # README
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    # Project webpage
    url='https://github.com/Argonne-National-Laboratory/neml',
    # Author
    author='Argonne National Laboratory',
    # email
    author_email = 'messner@anl.gov',
    # Various things for pypi
    classifiers=[
      'Intended Audience :: Science/Research',
      'License :: OSI Approved :: MIT License',
      'Programming Language :: C++',
      'Programming Language :: Python :: 3',
      'Operating System :: OS Independent'
      ],
    # Which version of python is needed
    python_requires='>=3.6',
    # Keywords to help find
    keywords='materials structures modeling',
    # Definitely not zip safe
    zip_safe=False,
    # Require the various modules
    ext_modules=[CMakeExtension('neml')],
    # Add the CMake builder
    cmdclass=dict(build_ext=CMakeBuild),
    # Get the python files
    packages=find_packages(),
    # Locate tests
    test_suite='nose.collector',
    tests_require=['nose'],
    # Python dependencies
    install_requires=[
      'numpy',
      'scipy',
      'matplotlib',
      'networkx'
      ],
)
