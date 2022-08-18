# build_f2py_diff.sh
# This script builds the Python interface for the fem_module.f90
# The differentiated modes are also included here

# Define variables
MOD=asa_module
AMOD=llt_module
SMOD=fem_module
FFLAGS="-O2 -fPIC -fdefault-real-8 -fbounds-check"
CFLAGS="-fPIC"

# Remove previous files
rm *.o
rm *.mod

# Make a copy of the file that maps Fortran types to Python types
cp f2py_f2cmap .f2py_f2cmap

# Compile Fortran codes
gfortran -c $FFLAGS ${AMOD}.f90 ${SMOD}.f90 ${MOD}.f90

# Create the f2py signature for the main module
# This generates the .pyf file
f2py3 ${MOD}.f90 -m ${MOD} -h ${MOD}.pyf --overwrite-signature

# Create the f2py module
# Add all *.o here
f2py3 -lgfortran --f90flags="$FFLAGS" -c -DF2PY_REPORT_ON_ARRAY_COPY=1 ${AMOD}.o ${SMOD}.o ${MOD}.o ${MOD}.pyf

# Clean-up
rm *.o
rm *.mod
rm .f2py_f2cmap
rm *.pyf
