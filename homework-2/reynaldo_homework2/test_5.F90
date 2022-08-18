program main

	! MODULES 
	use llt_module, only : norm
	use LLT_MODULE_DIFF_D, only : get_residuals_d, get_functions_d
	use LLT_MODULE_B, only : GET_RESIDUALS_B, GET_FUNCTIONS_B, GET_RESIDUALS

	! DECLARATIONS
	implicit none
	integer :: ii
	real    :: uinf(3), X(3,4), Xd(3,4), resd(3), res(3), Gama(3), Gamad(3)
	real    :: FD(3), AD(3), res2(3), deltaRes(3), alpha0(3), alpha0d(3)
	real    :: chords(3),chordsd(3), cl0(3), cla(3), rho, Sref, CL, CD, L, D
	real    :: delta, resb(3), xb(3,4), gamab(3), alpha0b(3), chordsb(3)
	real    :: srefd, cld, cdd, ld, dd, cdb, clb, srefb, db, lb
	real    :: dotprod, dotprod2
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
	
	! RESB VALUE
	resb(1) = 0.1
	resb(2) = -0.2
	resb(3) = 0.3
	 Srefb = 0.6
	CLb = -0.2
	CDb = 0.1
	Lb = -0.3
	Db = 0.7
	
	CALL GET_RESIDUALS_D(3, x, xd, gama, gamad, alpha0, alpha0d,&
&	chords, chordsd, cl0, cla, uinf, rho, res2, resd)
	
	print *, 'Give the derivative seeds of the outputs'
	print *, '----------------------------------------'
	print *, 'resd = ', resd	
	
	call GET_FUNCTIONS_D(3, x, xd, gama, gamad, alpha0, alpha0d&
&   , chords, chordsd, cl0, cla, uinf, rho, sref, srefd, cl, cld, cd, &
&   cdd, l, ld, d, dd)
	
	! DEFINING A BUNCH OF OTHER STUF
	chordsb(1) = 0.0
	chordsb(2) = 0.0
	chordsb(3) = 0.0
	alpha0b(1) = 0.0
	alpha0b(2) = 0.0
	alpha0b(3) = 0.0
	gamab(1) = 0.0
	gamab(2) = 0.0
	gamab(3) = 0.0
	Xb(1,1) = 0.0
	Xb(1,2) = 0.0
	Xb(1,3) = 0.0
	Xb(1,4) = 0.0
	Xb(2,1) = 0.0
	Xb(2,2) = 0.0
	Xb(2,3) = 0.0
	Xb(2,4) = 0.0
	Xb(3,1) = 0.0
	Xb(3,2) = 0.0
	Xb(3,3) = 0.0
	Xb(3,4) = 0.0
	res(1) = 0
	res(2) = 0
	res(3) = 0
	print *,
	print *, 'Give the derivative seeds of the outputs'
	print *, '----------------------------------------'
	print *, 'Srefd = ', Srefd	
	print *, 'CLd   = ', CLd
	print *, 'CDd   = ', CDd
	print *, 'Ld    = ', Ld
	print *, 'Dd    = ', Dd
	print *,
	call GET_FUNCTIONS_B(3, x, xb, gama, gamab, alpha0, alpha0b&
&   , chords, chordsb, cl0, cla, uinf, rho, sref, srefb, cl, clb, cd, &
&   cdb, l, lb, d, db)
	resb(1) = 0.1
	resb(2) = -0.2
	resb(3) = 0.3
	 Srefb = 0.6
	CLb = -0.2
	CDb = 0.1
	Lb = -0.3
	Db = 0.7
	print *, 'Give the derivative seeds of the inputs'
	print *, '----------------------------------------'
	print *, 'Xb      = ', Xb(1,:)
	print *, '          ', Xb(2,:)
	print *, '          ', Xb(3,:)
	print *, 'Gamab   = ', Gamab
	print *, 'alpha0b = ', alpha0b
	print *, 'Chordsb = ', chordsb
	print *,
	dotprod = sum(Xd*Xb) + sum(Gamad*Gamab) + sum(alpha0d*alpha0b) + sum(chordsd*chordsb) &
!&	- sum(resd*resb)
&	- Srefd*Srefb - CLd*CLb - CDd*CDb -Ld*Lb - Dd*Db
	print *, 'Write the value you obtained for dotprod'
	print *, '----------------------------------------'
	print *, dotprod
	print *, 
	! DEFINING A BUNCH OF OTHER STUF
	chordsb(1) = 0.0
	chordsb(2) = 0.0
	chordsb(3) = 0.0
	alpha0b(1) = 0.0
	alpha0b(2) = 0.0
	alpha0b(3) = 0.0
	gamab(1) = 0.0
	gamab(2) = 0.0
	gamab(3) = 0.0
	Xb(1,1) = 0.0
	Xb(1,2) = 0.0
	Xb(1,3) = 0.0
	Xb(1,4) = 0.0
	Xb(2,1) = 0.0
	Xb(2,2) = 0.0
	Xb(2,3) = 0.0
	Xb(2,4) = 0.0
	Xb(3,1) = 0.0
	Xb(3,2) = 0.0
	Xb(3,3) = 0.0
	Xb(3,4) = 0.0
	call get_residuals(3, X, Gama, alpha0, chords, cl0, cla, uinf, rho,res)

	call GET_RESIDUALS_B(3, x, xb, gama, gamab, alpha0, alpha0b&
&   , chords, chordsb, cl0, cla, uinf, rho, res, resb)
	print *,
	print *, 'Give the derivative seeds of the inputs'
	print *, '----------------------------------------'
	print *, 'Xb      = ', Xb(1,:)
	print *, '          ', Xb(2,:)
	print *, '          ', Xb(3,:)
	print *, 'Gamab   = ', Gamab
	print *, 'alpha0b = ', alpha0b
	print *, 'Chordsb = ', chordsb
	print *,
	resb(1) = 0.1
	resb(2) = -0.2
	resb(3) = 0.3
	dotprod2 = sum(Xd*Xb) + sum(Gamad*Gamab) + sum(alpha0d*alpha0b) + sum(chordsd*chordsb) &
&	- sum(resd*resb)
	print *, 'Write the value you obtained for dotprod'
	print *, '----------------------------------------'
	print *, dotprod2
	print *, 
end program main

