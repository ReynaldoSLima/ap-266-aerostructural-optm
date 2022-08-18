gfortran -c llt_module.f90 test_1.F90

gfortran llt_module.o test_1.o -o test_1

./test_1
