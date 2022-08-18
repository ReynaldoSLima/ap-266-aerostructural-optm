program main

	! MODULES 
	use llt_module, only : induced_vel

	! DECLARATIONS
	implicit none
	real :: uinf(3), x1(3), x2(3), x(3), ind_vel(3)
	x1(1)   = 0.0
	x1(2)   = 0.0
	x1(3)   = 0.0
	x2(1)   = 0.0
	x2(2)   = 1.0
	x2(3)   = 0.0
	x(1)    = 0.0
	x(2)    = 0.0
	x(3)    = 1.0
	uinf(1) = 1.0
	uinf(2) = 0.5
	uinf(3) = 0.0
	
	! EXECUTION TEST

	call induced_vel(uinf, x1, x2, x, ind_vel)
	print *, 'vel_ind, test case: ', ind_vel

	! EXECUTION CASE 1

	x2(1)   = 0.0
	x2(2)   = 2.0
	x2(3)   = 0.0
	x(1)    = 0.0
	x(2)    = 3.0
	x(3)    = 0.0
	uinf(1) = 1.0
	uinf(2) = 0.0
	uinf(3) = 0.0
	
	call induced_vel(uinf, x1, x2, x, ind_vel)
	print *, 'vel_ind, case 1: ', ind_vel

	! EXECUTION CASE 2

	x2(1)   = 0.0
	x2(2)   = 2.0
	x2(3)   = 0.0
	x(1)    = 0.0
	x(2)    = -1.0
	x(3)    = 0.0
	
	call induced_vel(uinf, x1, x2, x, ind_vel)
	print *, 'vel_ind, case 2: ', ind_vel

	! EXECUTION CASE 3

	x2(1)   = 0.5
	x2(2)   = 2.0
	x2(3)   = 0.5
	x(1)    = 0.0
	x(2)    = 1.0
	x(3)    = 0.0

	call induced_vel(uinf, x1, x2, x, ind_vel)
	print *, 'vel_ind, case 3: ', ind_vel

	! EXECUTION CASE 4

	x2(1)   = 0.0
	x2(2)   = 2.0
	x2(3)   = 0.0
	x(1)    = -0.1
	x(2)    = 3.0
	x(3)    = -0.5
	

	call induced_vel(uinf, x1, x2, x, ind_vel)
	print *, 'vel_ind, case 4: ', ind_vel
	
	! EXECUTION CASE 5

	x2(1)   = 0.0
	x2(2)   = 2.0
	x2(3)   = 0.0
	x(1)    = 0.0
	x(2)    = 1.0
	x(3)    = 0.0

	call induced_vel(uinf, x1, x2, x, ind_vel)
	print *, 'vel_ind, case 5: ', ind_vel


end program main

