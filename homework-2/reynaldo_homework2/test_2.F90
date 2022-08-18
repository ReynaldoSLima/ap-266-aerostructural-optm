program main

	! MODULES 
	use llt_module, only : get_residuals, get_functions

	! DECLARATIONS
	implicit none
	integer :: ii
	real    :: uinf(3), X(3,4), res(3), Gama(3), alpha0(3), chords(3), cl0(3), cla(3), rho, Sref, CL, CD, L, D
	
	! DEFINING X
	X(1,1) = 0.0
	X(1,2) = 0.0
	X(1,3) = 0.0
	X(1,4) = 0.0
	X(2,1) = 0.0
	X(2,2) = 1.0
	X(2,3) = 2.0
	X(2,4) = 3.0
	X(3,1) = 0.0
	X(3,2) = 0.0
	X(3,3) = 0.0
	X(3,4) = 0.0

	! DEFINING UINF
	uinf(1) = 1.0
	uinf(2) = 0.0
	uinf(3) = 0.0	

	! DEFINING A BUNCH OF OTHER STUF
	rho = 1.0
	cla(1) = 6.283
	cla(2) = 6.283	
	cla(3) = 6.283
	cl0(1) = 0.0
	cl0(2) = 0.0
	cl0(3) = 0.0
	chords(1) = 0.3
	chords(2) = 0.3
	chords(3) = 0.3
	alpha0(1) = 0.0
	alpha0(2) = 0.0
	alpha0(3) = 0.0
	Gama(1) = 2.0
	Gama(2) = 2.0
	Gama(3) = 2.0
	call get_residuals(3, X, Gama, alpha0, chords, cl0, cla, uinf, rho,res)
	print *, 'res  = ', res
	call get_functions(3, X, Gama, alpha0, chords, cl0, cla, uinf, rho, Sref, CL, CD, L, D)
	print *, 'Sref = ', Sref
	print *, 'CL   = ', CL
	print *, 'CD   = ', CD
	print *, 'L    = ', L
	print *, 'D    = ', D

	! ALL OVER AGAIN, BUT NOT THE SAME
	! DEFINING X
	X(1,1) = 0.0
	X(1,2) = 0.0
	X(1,3) = 0.0
	X(1,4) = 0.0
	X(2,1) = 0.0
	X(2,2) = 1.0
	X(2,3) = 2.0
	X(2,4) = 3.0
	X(3,1) = 0.2
	X(3,2) = 0.0
	X(3,3) = 0.0
	X(3,4) = 0.2

	! DEFINING UINF
	uinf(1) = 1.0
	uinf(2) = 0.0
	uinf(3) = 0.0	

	! DEFINING A BUNCH OF OTHER STUF
	rho = 1.0
	cla(1) = 6.283
	cla(2) = 6.283	
	cla(3) = 6.283
	cl0(1) = 0.0
	cl0(2) = 0.0
	cl0(3) = 0.0
	chords(1) = 0.3
	chords(2) = 0.6
	chords(3) = 0.3
	alpha0(1) = 0.0
	alpha0(2) = 0.1
	alpha0(3) = 0.0
	Gama(1) = 1.0
	Gama(2) = 2.0
	Gama(3) = 1.0
	call get_residuals(3, X, Gama, alpha0, chords, cl0, cla, uinf, rho,res)
	print *, 'res  = ', res
	call get_functions(3, X, Gama, alpha0, chords, cl0, cla, uinf, rho, Sref, CL, CD, L, D)
	print *, 'Sref = ', Sref
	print *, 'CL   = ', CL
	print *, 'CD   = ', CD
	print *, 'L    = ', L
	print *, 'D    = ', D
	
end program main

