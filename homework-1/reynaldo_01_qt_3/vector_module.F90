module vector_module

contains

	subroutine dot(nvec, a, b, dp)
	! Dot product for 3-dimensional vectors (it's general, actl)
		
		! DECLARATIONS
		implicit none

		! Input variables
		integer, intent(in) :: nvec
		real, intent(in) :: a(nvec), b(nvec)
		
		! Output variables
		real, intent(out) :: dp
		
		! Working variables
		integer :: ii

		! EXECUTION
		dp = 0
		do ii = 1, nvec
			dp = dp + a(ii)*b(ii)
		end do
		
	end subroutine dot

	!===================================

	subroutine cross(a, b, dp)
	! Cross product for 3-dimensional vectors (in this case, there is no such sense in doing a general cross product)
		
		! DECLARATIONS
		implicit none

		! Input variables
		real, intent(in) :: a(3), b(3)
		
		! Output variables
		real, intent(out) :: dp(3)

		! EXECUTION
		dp(1) = a(2)*b(3) - a(3)*b(2)
		dp(2) = a(3)*b(1) - a(1)*b(3)
		dp(3) = a(1)*b(2) - a(2)*b(1)

	end subroutine cross

	!===================================

	subroutine norm(nvec, a, dp)
	! Norm for 3-dimensional vectors (it's general, actl)
		
		! DECLARATIONS
		implicit none

		! Input variables
		integer, intent(in) :: nvec
		real, intent(in) :: a(nvec)
		
		! Output variables
		real, intent(out) :: dp
		
		! Working variables
		integer :: ii

		! EXECUTION
		dp = 0
		do ii = 1, nvec
			dp = dp + a(ii)*a(ii)
		end do
		
		dp = sqrt(dp)

	end subroutine norm

	!===================================

end module vector_module
