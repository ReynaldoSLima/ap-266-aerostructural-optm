program main

	! MODULES 
	use llt_module, only : get_residuals, get_functions, norm
	use LLT_MODULE_DIFF_D, only : get_residuals_d, get_functions_d

	! DECLARATIONS
	implicit none
	integer :: ii
	real    :: uinf(3), X(3,4), Xd(3,4), res(3), res1(3), Gama(3), Gamad(3)
	real    :: FD(3), AD(3), res2(3), deltaRes(3), h, alpha0(3), alpha0d(3)
	real    :: chords(3),chordsd(3), cl0(3), cla(3), rho, Sref, CL, CD, L, D
	real    :: deltaSref, CL2, CD2, L2, D2, Sref2
	real    :: srefd, cld, cdd, ld, dd
	real    :: delta
	
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
	! DEFINING XD
	Xd(1,1) = -0.1
	Xd(1,2) = 0.5
	Xd(1,3) = 0.2
	Xd(1,4) = 0.0
	Xd(2,1) = 1.0
	Xd(2,2) = -0.7
	Xd(2,3) = 0.1
	Xd(2,4) = -0.9
	Xd(3,1) = 0.3
	Xd(3,2) = 0.0
	Xd(3,3) = -0.6
	Xd(3,4) = 0.2

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
	
	! SEEDS!!!!!
	chordsd(1) = 1.0
	chordsd(2) = -0.2
	chordsd(3) = 0.7
	alpha0d(1) = 0.0
	alpha0d(2) = 0.1
	alpha0d(3) = -0.1
	Gamad(1) = -1.0
	Gamad(2) = 0.5
	Gamad(3) = 0.3
	
	! H VALUE
	h = 1E-4

	call get_residuals(3, X, Gama, alpha0, chords, cl0, cla, uinf, rho,res)
	call get_residuals(3, X+h*Xd, Gama+h*Gamad, alpha0 + h*alpha0d, chords+&
&	h*chordsd, cl0, cla, uinf, rho, res1)
	FD = (res1 - res)/h
	CALL GET_RESIDUALS_D(3, x, xd, gama, gamad, alpha0, alpha0d,&
&	chords, chordsd, cl0, cla, uinf, rho, res2, AD)
	deltaRes = FD-AD
	call norm(deltaRes, delta)
	print *, 'Delta : ', delta
	print *,

	! OUTPUT VALUES

	call get_functions(3, X, Gama, alpha0, chords, cl0, cla, uinf, rho, Sref, &
& 	CL, CD, L, D)
	call get_functions(3, X+h*Xd, Gama+h*Gamad, alpha0 + h*alpha0d, chords+&
&	h*chordsd, cl0, cla, uinf, rho, Sref2, CL2, CD2, L2, D2)
	call GET_FUNCTIONS_D(3, x, xd, gama, gamad, alpha0, alpha0d&
&   , chords, chordsd, cl0, cla, uinf, rho, sref, srefd, cl, cld, cd, &
&   cdd, l, ld, d, dd)
	print *, '--------------------------------------------------------------------------------'
	print *, 'resD  = ', FD, 'DF method'
	print *,
	
	print *, 'SrefD = ',-(sref-sref2)/h, 'DF method'
	
	print *, 'CLD   = ',-(cl-cl2)/h, 'DF method'
	
	print *, 'CDD   = ',-(cd-cd2)/h, 'DF method'
	
	print *, 'LD    = ',-(l-l2)/h, 'DF method'
		
	print *, 'DD    = ',-(d-d2)/h, 'DF method'
	print *,
	print *, '--------------------------------------------------------------------------------'
	print *, 'resD  = ', AD
	print *,
	print *, 'SrefD = ', srefd
	print *, 'CLD   = ', cld
	print *, 'CDD   = ', cdd
	print *, 'LD    = ', ld
	print *, 'DD    = ', dd 
	print *, '--------------------------------------------------------------------------------'
	
end program main

