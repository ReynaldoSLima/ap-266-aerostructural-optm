program main2

	! MODULES 
	use llt_module, only : induced_vel

	! DECLARATIONS
	integer :: nv, ii
	real :: uinf(3), x(3) 
	real, dimension(3,4) :: X_v
	real, dimension(3,3) :: ind_vel
	
	! Number of vortices
	nv = 3
	
	! INPUT
	uinf(1) = 1.0
	uinf(2) = 0.0
	uinf(3) = 0.0
	x(1) = 0.0
	x(2) = 0.5
	x(3) = 0.0
	X_v =reshape( (/0.0, 0.0, 0.0, 0.0, 1.0, 0.1, 0.0, 2.0, 0.2, 0.0, 3.0, 0.3 /), (/ 3, 4 /))
	
	! EXECUTION
	do ii = 1, nv
		call induced_vel(uinf, X_v(:, ii), X_v(:, ii+1), x, ind_vel(:, ii))
		print *, X_v(:, ii)
		print *, X_v(:, ii+1)
	end do	
	

	print *, 'Matrix of seq. linked horseshoes vortices: ', ind_vel
	

end program main2

