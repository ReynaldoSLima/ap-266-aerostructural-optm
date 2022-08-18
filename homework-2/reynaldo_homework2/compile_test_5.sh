# Compile most used files from the ADFirstAidKit
gcc -c adStack.c
gfortran -c -fdefault-real-8 adBuffer.f

# Compile
gfortran -c -fdefault-real-8 llt_module_b.f90 llt_module_d.f90 test_5.F90
# Link
gfortran adStack.o adBuffer.o llt_module_b.o llt_module_d.o test_5.o -o test_5

# Execute
./test_5


