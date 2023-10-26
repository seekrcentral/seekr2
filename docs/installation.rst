Installation
============

At this time, SEEKR2 has only been tested on Linux systems. Therefore, all
installation instructions are for Linux only.

First, you must choose which MD engine you will use: either OpenMM or NAMD.
Each engine has their own advantages - OpenMM is faster on GPUs and is likely
to give slightly more accurate results in SEEKR2. NAMD is optimized for 
distributed computing systems, such as supercomputers or clusters which use 
large numbers of CPU cores. If you are uncertain which to choose, OpenMM is 
a good default choice.

Install Conda
-------------

It is recommended, though not mandatory, that you install Conda to use with 
SEEKR2. Without Conda, all dependencies will need to be installed by hand.

If you do not already have Conda, it can be easily installed by completing the
following steps:

Download Conda, run the script, and fill out the prompts::

  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh

Make sure Conda is installed by running:

``which conda``

You will want to use Python 3.8, so you can see which version you are with
the command:

``python -V``

If it says any other version besides Python 3.8, then enter:

``conda install python=3.8``

If you want you can create a conda environment, 

``conda create --name SEEKR python=3.8``

but you can also just install all packages straight to the base environment
if you wish to. If using an environment, whenever you're installing or running 
anything involving OpenMM or SEEKR2, make sure that you have activated your 
environment by running ``conda activate SEEKR``.

Install OpenMM and Plugin with Conda
------------------------------------
This section describes the fasted and easiest way to get SEEKR2 working.

If you desire to use OpenMM, you must install OpenMM either from conda or from 
source. If you wish to install from source, see the "Installing OpenMM from
Source" sections below.

If you desire to use NAMD, then see the "Install NAMD" section below.

With Conda working, You may create and activate any environment you wish, 
or use the base environment. Install the SEEKR2 OpenMM Plugin:

``conda install -c conda-forge seekr2_openmm_plugin``

OpenMM will be installed automatically alongside the plugin.

One can test the installation by opening a Python terminal and typing:

``import seekr2plugin``

If you get an error such as "No module named seekr2plugin", you might
need to install with CUDA Toolkit version 10.2 and OpenMM 7.7:

``conda install -c conda-forge seekr2_openmm_plugin cudatoolkit=10.2 openmm=7.7``


Installation of SEEKR2 itself begins with cloning and installing the SEEKR2 
python API::

  git clone https://github.com/seekrcentral/seekr2.git
  cd seekr2
  python -m pip install .
  
  
Once OpenMM and the OpenMM SEEKR2 Plugin are installed, it is recommended that 
you run tests of SEEKR2. From within the "seekr2/" directory, run:

``pytest``

One or two tests may fail depending on whether NAMD2 and/or Browndye2 have been
installed, and can be safely ignored if those programs are not needed.

Additional continuous integration tests may be run from the Python scripts in
the seekr2/seekr2/continuous_integration/ directory if extra testing is
desired.

You should now be able to use SEEKR2.

Install NAMD (If not using OpenMM)
----------------------------------
If you desire not to use OpenMM, but rather to use NAMD for your MD 
calculations, you should follow the 
`NAMD installation instructions <https://www.ks.uiuc.edu/Research/namd/2.9/ug/node91.html>`_

NAMD is also often available on shared scientific computing resources such as
most supercomputers and clusters. Consult the resource manual or system
administrators to see if NAMD2 is installed on the available shared resource.

SEEKR2 itself must be installed by cloning and installing the SEEKR2 
python API::

  git clone https://github.com/seekrcentral/seekr2.git
  cd seekr2
  python -m pip install .

OpenMM and Plugin Installation from Source on Local Machine (If not using Conda to install OpenMM and Plugin)
-------------------------------------------------------------------------------------------------------------
Compiling OpenMM from source is tricky, but may be desirable if the Conda 
installation doesn't work, or if you wish to optimize OpenMM's performance.

Please see the official 
`OpenMM guide to installing from source <http://docs.openmm.org/latest/userguide/library.html#compiling-openmm-from-source-code>`_ 
for complete OpenMM installation instructions. 

If you want to use a GPU to accelerate your OpenMM simulations (highly 
recommended) you must ensure that a recent version of CUDA is installed and
loaded. 

If you need to install CUDA, it is highly recommended that you contact your 
system administrator about this, although if you have to install CUDA by 
yourself, you should carefully read and follow all instructions from 
`NVIDIA's CUDA toolkit installation instructions 
<https://developer.nvidia.com/cuda-toolkit>`_ or 
https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html
.

Many times, cuda is located in /usr/local/cuda::

  ls /usr/local/cuda
  
If CUDA is located here, then the OpenMM plugin should be able to automatically
detect it there. If CUDA is not in /usr/local/cuda, then you can also sometimes
find the CUDA compile 'nvcc' using 'which'. You can also see whether the 
CUDA_HOME environmental variable is defined::

  which nvcc
  echo $CUDA_HOME
  
If the commands didn't return a path to nvcc, or a value or CUDA_HOME, SEEKR2
is likely to have difficulty finding CUDA on it's own. You may have to take
more trouble to explicitly assign the necessary variables to the cmake or 
ccmake commands.  

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
ccmake or you can install it yourself if you have sudo privileges:

``sudo apt-get install cmake-curses-gui``

Make sure 'doxygen' is installed.

``conda install -c conda-forge doxygen``

Upgrade Cython:

``pip install --upgrade cython``

Clone OpenMM and cd into OpenMM directory, then perform necessary build steps.::

  git clone https://github.com/openmm/openmm.git
  cd openmm
  mkdir build
  cd build
  ccmake ..

The ccmake gui should come up. Press 'c' and then 't'

You should modify the following variables:

CMAKE_INSTALL_PREFIX: change to a local directory that exists (for example: 
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
  python -m openmm.testInstallation

If the PythonInstall step fails, then make sure you have upgraded cython

``pip install --upgrade cython``

Hopefully, with the final step, all the tests pass. If a few fail, then 
determine if those failures will be necessary for your calculations. If 
several or all fail, then you'll need to be sure that you fix whatever 
problem caused those failures. If the CUDA tests failed, then you either do
not have a working CUDA installation, or the proper environmental variables
such as OPENMM_CUDA_COMPILER have not been set.

You'll need to install the SEEKR2 Plugin on top of this version of OpenMM::

  cd ~ # or another directory of your choice
  git clone https://github.com/seekrcentral/seekr2_openmm_plugin.git
  cd seekr2_openmm_plugin/seekr2plugin
  mkdir build
  cd build
  ccmake ..
  
Now the ccmake gui should come up. Press 'c'.

You should modify the following variables:

* CMAKE_INSTALL_PREFIX and OPENMM_DIR: change to the directory that was
  CMAKE_INSTALL_PREFIX for the OpenMM installation above (example: 
  /home/USERNAME/bin/openmm).

* SEEKR2_BUILD_OPENCL_LIB should be set to OFF.

Close the GUI by pressing 'c' and then 'g'. Then make the plugin::
  
  make
  make install
  make PythonInstall
  make test # Optional

Installation of SEEKR2 itself begins with cloning and installing the SEEKR2 
python API::

  git clone https://github.com/seekrcentral/seekr2.git
  cd seekr2
  python -m pip install .

At this point, its a good idea to run the SEEKR2 tests. Navigate to where the 
"seekr2" git repository was cloned. From within the "seekr2/" directory, run:

``pytest``

If you get an error like "ImportError: libOpenMM.so.7.7: cannot open shared 
object file: No such file or directory", you will need to point your
LD_LIBRARY_PATH to the installed OpenMM library location::

  export LD_LIBRARY_PATH="/home/USERNAME/bin/openmm/lib:$LD_LIBRARY_PATH"
  export LD_LIBRARY_PATH="/home/USERNAME/bin/openmm/lib/plugins:$LD_LIBRARY_PATH"
  
Of course, change your path to be the actual location where CMAKE_INSTALL_PREFIX
was pointing.

OpenMM Installation from Source on Cluster or Shared Resource
-------------------------------------------------------------

A simple Conda installation on a Cluster or Supercomputer would probably
work just fine, but if you wish to install from source, this section provides
some helpful information to that end.

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
