#######################################################################
# Fast Auxiliary Space Preconditioners (FASP) 
#
################# User Defined Configuration Options #################
#
# 1. Copy this file to a file named "FASP.mk".
# 2. Edit "FASP.mk" to adjust options/settings for your system
#    following the directions below.
# 3. Type "make help" to see all build and configuration options.
########################################################################
#
# The default setting for build type for FASP is RELEASE. The RELEASE 
# build type by default has the "-O3". You may adjust the optimization
# compilation options according to your hardware and software setting.
# For example, on a macbook pro with Intel i7, best options could be
# "-Ofast -march=corei7 -mtune=corei7".
#
# If you want to work with build type DEBUG, then uncomment the next
# line (to include "-Wall -g")
#
# debug=yes
#
# In order to get debug output during run-time, you can uncomment the 
# following line (to include "-Wall -g -DDEBUG_MODE=3"):
# 
# debug=all
#
# The default setting for vebosity level for FASP is verbose=no. If you
# want to increase verbosity level, uncomment the next line:
#
# verbose=yes
#
# By default, FASP generates static libraries. If you need to generate 
# shared libs instead of static libs, uncomment the next line:
#
# shared=yes
#
# You may use multithread version after you enable OpenMP support. To
# setup the environment, you need 
#  >> export OMP_NUM_THREADS=4 (for bash)
#  >> setenv OMP_NUM_THREADS 4 (for tcsh)
# If you want to compile with OpenMP support, uncomment the next line:
#
# openmp=yes
#
# These user options can also be applied as make command line options.
# For example, to enforce the debug compiling options:
#
# make config debug=yes
#
#-------------------------------------------------------------------------
#
# By default, FASP uses the command-line Doxygen to generate a reference
# manual. If you want to use the GUI of Doxgen instead of command-line
# (if there is one installed on your system), uncomment the next line:
#
# doxywizard=yes
#
#-------------------------------------------------------------------------
# If you want to use UMFPACK (part of SuiteSparse package), uncomment the next 
# line (and read carefully the instructons below it):
# 
# umfpack=yes
#
# If you have installed SuiteSparse from source or for some other
# reason you want to specify the path to SuiteSparse libraries and
# header files, uncomment and edit the definition of "suitesparse_dir"
# below (and continue reading...)  
#
# suitesparse_dir="/path/to/SuiteSparse"
#
# IMPORTANT:
# This defines the path to the SuiteSparse library and include files.
# These are expected to be found in $(suitesparse_dir)/lib and
# $(suitesparse_dir)/include or in the system standard paths for libraries
# and header files. As a bare minimum, $(suitesparse_dir)/lib and
# $(suitesparse_dir)/include must contain the library and header files for
# AMD,UMFPACK, and SUITESPARSECONFIG. 
# -------------------------------------------------------------------------
# If you want to use SuperLU, uncomment the next line:
#
# superlu=yes
#
# If you want to specify the path to SuperLU, uncomment the next line
# and give the correct path to SuperLU here. For example:
#
# superlu_dir="/path/to/SuperLU"
#
#-------------------------------------------------------------------------
# If you want to use MUMPS, uncomment the next line:
#
 mumps=yes
#
# If you want to specify the path to MUMPS, uncomment the next line
# and give the correct path to MUMPS here. For example:
#
 mumps_dir="/opt/MUMPS_4.10.0.ifort"
# mumps_dir="/home/spring/spring/2020fasp/MUMPS_4.10.0.ifort"
#
#-------------------------------------------------------------------------
# If you want to use Intel MKL PARDISO, uncomment the next line:
#
# pardiso=yes
#
# If you want to specify the path to MKL, uncomment the next line
# and give the correct path to MKL here. For example:
#
  mkl_dir="/opt/intel/mkl/"
#
#  MKL_INCLUDE_DIRS: /opt/intel/composer_xe_2015.2.164/mkl/include
#  MKL_LIBRARY_DIRS: /opt/intel/composer_xe_2015.2.164/mkl/lib/intel64;/opt/intel/composer_xe_2015.2.164/compiler/lib/intel64
#-------------------------------------------------------------------------
# If you want to install the fasp library in another location, put the
# path below. It will install FASP library in $(prefix)/lib and the
# header files in $(prefix)/include
# prefix = /path/to/install
####################  User Defined Compiler Flags  #####################
#export CC=gcc
#cflags   = "-O3 -funroll-loops -funswitch-loops"
#cflags   = "-O3 -funroll-loops -funswitch-loops -march=corei7-avx -ffast-math"
#cxxflags = "-O3 -funroll-loops -funswitch-loops"
#cxxflags = "-O3 -funroll-loops -funswitch-loops -march=corei7-avx -ffast-math"
#fflags   = "-O3 -funroll-loops -funswitch-loops"
#fflags   = "-O3 -funroll-loops -funswitch-loops -march=corei7-avx -ffast-math"
