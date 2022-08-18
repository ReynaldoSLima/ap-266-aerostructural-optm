gfortran -c vector_module.F90 llt_module.F90 part_c.F90

gfortran vector_module.o llt_module.o part_c.o -o part_c

./part_c
