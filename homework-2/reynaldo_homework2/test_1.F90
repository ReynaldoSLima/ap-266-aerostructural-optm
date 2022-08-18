program main

	! MODULES 
	use llt_module, only : get_AIC_matrix

	! DECLARATIONS
	implicit none
	integer :: ii
	real    :: uinf(3), X(3,4), aix(3,3,3)
	
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

	call get_AIC_matrix(3, X, uinf, aix)
	do ii = 1, 3
		print *, 'aix: ', aix(3,ii,:)
	end do

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
	X(3,2) = 0.1
	X(3,3) = 0.2
	X(3,4) = 0.3
	
	print *, ' '
	
	call get_AIC_matrix(3, X, uinf, aix)
	do ii = 1, 3
		print *, aix(2,ii,:)
	end do

end program main

