Installation
============

At this time, OpenMMVT has only been tested on Linux systems. Therefore, all
installation instructions are for Linux only.

Installation begins with cloning and installing the OpenMMVT frontent::

  git clone https://github.com/lvotapka/openmmvt.git
  cd openmmvt
  python setup.py install

Next, you must choose which MD engine you will use: either OpenMM or NAMD.
Each engine has their own advantages - OpenMM is faster on GPUs and is likely
to give slightly more accurate results. NAMD is optimized for distributed
computing systems, such as supercomputers or cluster which use large numbers
of CPU cores. If you are uncertain which to choose, OpenMM is a good default.

Install Conda
-------------

It is recommended, though not mandatory, that you install Conda along with 
OpenMMVT.

If you do not already have Conda, it can be easily installed by completing the
following steps:

Download Conda, run the script, and fill out the prompts::

  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  sh Miniconda3-latest-Linux-x86_64.sh

Make sure Conda is installed by running:

``which conda``

If you want you can create a conda environment, but you can also just install 
all packages straight to the base environment. Whenever installing or running
anything involving OpenMM or OpenMMVT, make sure that you have activated your 
environment by running ``conda activate``.

Install OpenMM and Plugin
-------------------------
If you desire to use OpenMM, you must install OpenMM from source (for now). 
Please see the official `OpenMM guide to installing from source <http://docs.openmm.org/latest/userguide/library.html#compiling-openmm-from-source-code>`_ 
or see the "OpenMM Installation from Source" section(s) at the bottom of this
page.

Make sure to take note of the CMAKE_INSTALL_PREFIX variable, which will be 
refered to as /path/to/openmm.

If you desire to use NAMD, then see the "Install NAMD" section below.

After pulling the repo and installing the package, install the OpenMMVT plugin.
Depending on where OpenMM is installed, some commands may or may not need
to be run with sudo.::

  cd openmmvt/plugin
  mkdir build
  cd build
  cmake -DCMAKE_INSTALL_PREFIX=/path/to/openmm -DOPENMM_DIR=/path/to/openmm ..
  make
  make install
  make PythonInstall
  make test
  cd ../..

It is recommended that you run tests

``python setup.py test``

Install NAMD
------------
If you desire to NAMD for your MD calculations, you should follow the `NAMD
installation instructions <https://www.ks.uiuc.edu/Research/namd/2.9/ug/node91.html>`_

NAMD is also often available on shared scientific computing resources such as
most supercomputers and clusters. Consult the resource manual or system
administrators to see if NAMD is installed on the available shared resource.

OpenMM Installation from Source on Local Machine
------------------------------------------------
If you desire to run OpenMMVT on your local machine with OpenMM, you'll need
to perform the following steps.

If you want to use a GPU to accelerate your OpenMM simulations (highly 
recommended) you must ensure that a recent version of CUDA is installed and
loaded. It is highly recommended that you contact your system administrator
about this, although if you have to do it by yourself, you should carefully read
and follow all instructions from 
`NVIDIA's CUDA toolkit installation instructions 
<https://developer.nvidia.com/cuda-toolkit>`_.

In order to use CUDA, you may also need to define the following environmental
variable by placing it in your .bashrc file: 
``export OPENMM_CUDA_COMPILER=/path/to/nvcc``

Obviously, you'll need to modify "/path/to/nvcc" with the actual path. The 
program "nvcc" will exist in your CUDA installation, and might be discoverable 
by typing ``which nvcc``.

Next, install the necessary programs and packages into Conda.

``conda install numpy scipy netcdf4 mpi4py swig``

Make sure 'git' is installed, if not already.

``conda install git``

Make sure 'ccmake' is installed

``which ccmake``

If nothing happens, you may need to ask your system administrator to install 
ccmake or install it yourself if you have sudo privileges:

``sudo apt-get install cmake-curses-gui``

Make sure 'doxygen' is installed.

``conda install -c conda-forge doxygen``

Install Cython:

``pip install --upgrade cython``

Clone OpenMM and cd into OpenMM directory, then perform necessary build steps.::

  git clone https://github.com/openmm/openmm.git
  cd openmm
  mkdir build
  cd build
  ccmake ..

The ccmake gui should come up. Press 'c' and then 't'

You should modify the following variables:

CMAKE_INSTALL_PREFIX: change to a local directory that exists (example: 
/home/USERNAME/bin/openmm). If such a directory doesn't exist, then make one.
You can also leave this variable at the default if you have sudo privileges
and don't mind installing OpenMM globally.

Check all the variables, and then type 'c' to configure. If there are any 
problems, it will let you know.

When the configuration is successful, type 'g' to generate. Then ccmake 
should close on its own.

If you are having trouble with assigning a variable, like 
CUDA_CUDA_LIBRARY-NOTFOUND, then run 'cmake' (instead of 'ccmake') and 
assign the missing variable using the -D argument:

For example:
``cmake -DCMAKE_LIBRARY_PATH=/usr/local/cuda/lib64/stubs ..``

Next, build, install, and test OpenMM::

  make
  make install
  make PythonInstall
  make test

If the PythonInstall step fails, then make sure you have installed cython

``pip install --upgrade cython``

Hopefully, with the final step, all the tests pass. If a few fail, then 
determine if those failures will be necessary for our calculations. If 
several or all fail, then you'll need to be sure that you fix whatever 
problem caused those failures. If the CUDA tests failed, then you either do
not have a working CUDA installation, or the proper environmental variables
such as OPENMM_CUDA_COMPILER have not been set.

Try to see if the python interface works. Inside a python shell, try:

from simtk import openmm

If you see no errors, then your OpenMM installation was probably successful.

OpenMM Installation from Source on Cluster or Shared Resource
-------------------------------------------------------------

Installation of OpenMM on a shared resource is almost identical to the
local installation of OpenMM as detailed in the previous section. However, the
shared resource is likely to have a number of specific features that will have
to be taken into account when installing OpenMM from source.

Some tips and advice:

* You should run all installation commands in an interactive node to avoid 
  clogging up the login nodes. Consider using the debug or development queue,
  if available.

* Use "wget" to obtain miniconda: ``wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh``

* If the cluster or shared resource has GPU computing capabilities, the 
  administrators have likely made CUDA available. You should consult the 
  resource's manual or reach out to the system administrators for how to 
  load or utilize CUDA.

* If 'ccmake' is not available, you can still use 'cmake' to install OpenMM,
  you just must provide any arguments using '-D'. For instance: 
  ``cmake -DCMAKE_INSTALL_PREFIX=/path/to/openmm -DCMAKE_LIBRARY_PATH=/path/to/cuda/lib64/stubs ..``
