gfortran -c llt_module.f90 llt_module_d.f90 test_4.F90

gfortran llt_module.o llt_module_d.o test_4.o -o test_4

./test_4
