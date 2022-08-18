module llt_module

contains
	subroutine induced_vel(uinf, x1, x2, x, gamma_ind)
	! This subroutine compute the velocity induced by
	! a horseshoe vortec with kinks at x1 and x2 at the
	! point x.
	
		! DECLARATIONS
		use vector_module, only: dot, cross, norm

		! Input variables
			! uinf -> Unitary vector along the free-stream dir
			! x1   -> Coordinates of the first kink
			! x2   -> Coordinates of the second kink
			! x    -> Coordinates where we want the vel

		real, intent(in) :: uinf(3), x1(3), x2(3), x(3)

		! Output variables
		real, intent(out) :: gamma_ind(3)

		! Declare local constant Pi
		REAL, PARAMETER :: Pi = 3.1415927

		! Working variables
		real :: prev(3), cross_aux(3), dot_aux, r1(3), r2(3), norm_r1, norm_r2		

		! EXECUTION
		r1 = x - x1
		r2 = x - x2
		
		call norm(3, r1, norm_r1)
		call norm(3, r2, norm_r2)
		gamma_ind = (/0, 0, 0 /)
		! Calculating first term within the expression
		call dot(3, uinf, r2, dot_aux)
		call cross(uinf, r2, cross_aux)
		IF ((norm_r2*(norm_r2-dot_aux)) /= 0) THEN
			gamma_ind = 1/(4*Pi)*(cross_aux/(norm_r2*(norm_r2-dot_aux)))
		END IF
		! Calculating second term within the expression
		call dot(3, r1, r2, dot_aux)
		call cross(r1, r2, cross_aux)
		prev = gamma_ind
		IF ((norm_r2*norm_r1*(norm_r2*norm_r1+dot_aux)) /= 0) THEN
			gamma_ind = prev + 1/(4*Pi)*(norm_r1 + norm_r2)*(cross_aux/(norm_r2*norm_r1*(norm_r2*norm_r1+dot_aux)))
		END IF		
		! Calculating second term within the expression
		call dot(3, uinf, r1, dot_aux)
		call cross(uinf, r1, cross_aux)
		prev = gamma_ind
		IF ((norm_r1*(norm_r1-dot_aux)) /= 0) THEN
			gamma_ind = prev - 1/(4*Pi)*(cross_aux/(norm_r1*(norm_r1-dot_aux)))
		END IF
	
	end subroutine induced_vel

	!===================================
	
end module llt_module

