module llt_module

  ! Module that holds necessary subroutines for
  ! the generalized lifting line method.
  ! The formulation is based on:
  !   Modern Adaptation of Prandtl's Classic Lifting-Line Theory
  !   W. F. Phillips and D. O. Snyder
  !   Journal of Aircraft Vol. 37, No. 4, July-August 2000
  !
  ! Ney Secco 2019-03-24

  implicit none

  real, parameter :: zero = 0.0
  real, parameter :: pi = 3.1415926535897932384626
  real, parameter :: eps = epsilon(zero)

contains

  subroutine induced_vel(uinf, x1, x2, x_c, gamma_ind)
	! This subroutine compute the velocity induced by
	! a horseshoe vortec with kinks at x1 and x2 at the
	! point x.
	
	! DECLARATIONS

	real, intent(in) :: uinf(3), x1(3), x2(3), x_c(3)

	! Output variables
	real, intent(out) :: gamma_ind(3)

	! Declare local constant Pi
	REAL, PARAMETER :: Pi = 3.1415927

	! Working variables
	real :: prev(3), cross_aux(3), dot_aux, r1(3), r2(3), norm_r1, norm_r2		

	! EXECUTION
	r1 = x_c - x1
	r2 = x_c - x2
		
	call norm(r1, norm_r1)
	call norm(r2, norm_r2)
	gamma_ind = (/0, 0, 0 /)
	! Calculating first term within the expression
	call dot(uinf, r2, dot_aux)
	call cross(uinf, r2, cross_aux)
	IF ((norm_r2*(norm_r2-dot_aux)) > eps) THEN
		gamma_ind = 1/(4*Pi)*(cross_aux/(norm_r2*(norm_r2-dot_aux)))
	END IF
	! Calculating second term within the expression
	call dot(r1, r2, dot_aux)
	call cross(r1, r2, cross_aux)
	prev = gamma_ind
	IF ((norm_r2*norm_r1*(norm_r2*norm_r1+dot_aux)) > eps) THEN
		gamma_ind = prev + 1/(4*Pi)*(norm_r1 + norm_r2)*(cross_aux/(norm_r2*norm_r1*(norm_r2*norm_r1+dot_aux)))
	END IF		
	! Calculating second term within the expression
	call dot(uinf, r1, dot_aux)
	call cross(uinf, r1, cross_aux)
	prev = gamma_ind
	IF ((norm_r1*(norm_r1-dot_aux)) > eps) THEN
		gamma_ind = prev - 1/(4*Pi)*(cross_aux/(norm_r1*(norm_r1-dot_aux)))
	END IF
	
  end subroutine induced_vel

  !===================================
	

  !============================================================

  subroutine get_AIC_matrix(n_vort, X, Uinf, AIC)

    !   This subroutine computes the Aerodynamic Influence Coefficients (AIC)
    !   Input variables
	integer, intent(in) :: n_vort
	real, intent(in)    :: X(3,n_vort+1)
	real, intent(in)    :: Uinf(3)
    !   Auxiliar variables
	integer             :: ii, jj
	real                :: x_control(3)
	real                :: gamma_ind(3)
    !   Output variables
	real, intent(out)   :: AIC(3,n_vort,n_vort)
    
    !   Operations
	do ii = 1, n_vort
		x_control = (X(:,ii) + X(:,ii+1))/2
		do jj = 1, n_vort
			call induced_vel(Uinf, X(:,jj), X(:,jj+1), x_control, gamma_ind)
			AIC(:, jj, ii) = gamma_ind
		end do
	end do

	! Vind(1:3,ii) = AIC(1:3,:,ii)*Gamma(:)

  end subroutine get_AIC_matrix

  !============================================================

  subroutine get_geom_vectors(n_vort, X, alpha0, chords, Uai, Uni, si, areas)

    ! This subroutine computes some geometric properties of the
    ! aerodynamic strips.
    !
    ! INPUTS
    !
    ! n_vort: integer -> Number of horseshoe vortices.
    ! X: real(3,n_vort+1) -> Coordinates of the kinks of the horseshoe vortices.
    ! alpha0: real(n_vort) -> Array of local incidence angles.
    ! chords: real(n_vort) -> Array of local chords.
    !
    ! OUTPUTS
    !
    ! Uai: real(3,n_vort) -> Unitary vectors along local chords.
    ! Uni: real(3,n_vort) -> Unitary vectors along local normals.
    ! si: real(3,n_vort) -> Vectors along bound vortices. They are not unitary.
    ! areas: real(n_vort) -> Areas of each panel.

    implicit none

    ! Input variables
    integer, intent(in) :: n_vort
    real, intent(in) :: X(3,n_vort+1)
    real, intent(in) :: alpha0(n_vort), chords(n_vort)

    ! Output variables
    real, intent(out) :: Uai(3,n_vort), Uni(3,n_vort)
    real, intent(out) :: si(3,n_vort), areas(n_vort)

    ! Working variables
    integer :: ii
    real :: X1(3), X2(3), Usi(3), ni(3)
    real :: nim, sim

    ! EXECUTION

    ! Loop over all vortices
    do ii=1,n_vort

       ! Get kink coordinates of the current vortex
       X1 = X(:,ii)
       X2 = X(:,ii+1)

       ! Compute the vector along the local span
       si(:,ii) = X2 - X1
       call normalize(si(:,ii), Usi)

       ! Create a unitary vector along the chordline
       ! pointing to the trailing edge
       Uai(1,ii) = cos(alpha0(ii))
       Uai(2,ii) = zero
       Uai(3,ii) = -sin(alpha0(ii))

       ! Get the normal vector as the normal between the span
       ! and the chord vector
       call cross(Uai(:,ii), Usi, ni)
       call normalize(ni, Uni(:,ii))

       ! Compute the local area of the panel.
       ! The norm of the normal vector is the sine of the angle.
       call norm(ni, nim)
       call norm(si(:,ii), sim)
       areas(ii) = nim*chords(ii)*sim

    end do

  end subroutine get_geom_vectors

  !============================================================

  subroutine get_local_vels(n_vort, Gama, AIC, Vinf, Uai, Uni, Vlocal, alphaLocal)

    ! This subroutine computes the local velocity an angle of attack at each
    ! bound vortex for the given vorticity values.

    implicit none

    ! Input variables
    integer, intent(in)  :: n_vort
    real, intent(in) :: Gama(n_vort), AIC(3,n_vort,n_vort)
    real, intent(in) :: Vinf(3), Uai(3,n_vort), Uni(3,n_vort)

    ! Output variables
    real, intent(out) :: Vlocal(3,n_vort), alphaLocal(n_vort)

    ! Working variables
    integer :: ii,jj
    real :: num, den

    ! EXECUTION

    ! Do the matrix operation
    do ii=1,n_vort

       ! Start with the free-strem contribution
       Vlocal(:,ii) = Vinf

       ! Now add the velocities induced by each vortex
       do jj=1,n_vort
          Vlocal(:,ii) = Vlocal(:,ii) + AIC(:,jj,ii)*Gama(jj)
       end do

       ! Compute the local angle of attack
       call dot(Vlocal(:,ii), Uni(:,ii), num)
       call dot(Vlocal(:,ii), Uai(:,ii), den)
       alphaLocal(ii) = atan(num/den)

    end do

  end subroutine get_local_vels

  !============================================================

  subroutine compute_circulation_forces(n_vort, rho, Vlocal, si, Gama, Fcirc)

    ! This subroutine uses the circulation values to compute aerodynamic
    ! forces on every bound vortex.
    !
    ! INPUTS

    implicit none

    ! Input variables
    integer, intent(in) :: n_vort
    real, intent(in) :: rho, Vlocal(3,n_vort), si(3,n_vort)
    real, intent(in) :: Gama(n_vort)

    ! Output variables
    real, intent(out) :: Fcirc(3,n_vort)

    ! Working variables
    integer :: ii
    real :: li(3)

    ! EXECUTION

    do ii = 1,n_vort

       ! L=rho*V*gamma*l
       call cross(Vlocal(:,ii), si(:,ii), li)
       Fcirc(:,ii) = rho*Gama(ii)*li

    end do

  end subroutine compute_circulation_forces

  !============================================================

  subroutine compute_airfoil_forces(n_vort, rho, cl0, cla, areas, Vlocal, alphaLocal, si, Fairf)

    ! This subroutine compares the 2D forces and the forces given by
    ! the circulation.
    !
    ! INPUTS
    !
    ! cl0: real(n_vort) -> cl0 of each 2D section
    ! cla: real(n_vort) -> lift curve slope of each 2D section [1/rad]

    implicit none

    ! Input variables
    integer, intent(in) :: n_vort
    real, intent(in) :: rho, cl0(n_vort), cla(n_vort), areas(n_vort)
    real, intent(in) :: Vlocal(3,n_vort), alphaLocal(n_vort), si(3,n_vort)

    ! Output variables
    real, intent(out) :: Fairf(3,n_vort)

    ! Working variables
    real :: Fairfm, cl, li(3), Uli(3), Vlocalm
    integer :: ii

    ! EXECUTION

    ! Loop over all horseshoe vortices
    do ii=1,n_vort

       ! Compute the magnitude of the force given by the 2D section
       ! L = 0.5*Vlocal²*area*cl
       cl = cl0(ii) + alphaLocal(ii)*cla(ii)
       call norm(Vlocal(:,ii), Vlocalm)
       Fairfm = 0.5*rho*Vlocalm*Vlocalm*areas(ii)*cl

       ! Get normalized lift direction
       call cross(Vlocal(:,ii),si(:,ii), li)
       call normalize(li, Uli)

       ! Store the force
       Fairf(:,ii) = Fairfm*Uli

    end do

  end subroutine compute_airfoil_forces

  !============================================================

  subroutine get_residuals(n_vort, X, Gama, alpha0, chords, cl0, cla, Vinf, rho, res)

    ! This subroutine executes all steps of the LLT method to compute
    ! residuals for a given Gama.

    implicit none

    ! Input variables
    integer, intent(in) :: n_vort
    real, intent(in) :: X(3,n_vort+1)
    real, intent(in) :: Gama(n_vort), alpha0(n_vort)
    real, intent(in) :: chords(n_vort), cl0(n_vort)
    real, intent(in) :: cla(n_vort), Vinf(3), rho

    ! Output variables
    real, intent(out) :: res(n_vort)

    ! Working variables
    integer :: ii
    real :: Uinf(3), AIC(3,n_vort,n_vort)
    real :: Uai(3,n_vort), Uni(3,n_vort), si(3,n_vort)
    real :: areas(n_vort), Vlocal(3,n_vort), alphaLocal(n_vort)
    real :: Fcirc(3,n_vort), Fairf(3,n_vort), deltaF(3)

    ! EXECUTION

    ! Normalize free-stream velocity
    call normalize(Vinf, Uinf)

    ! Get AIC matrix
    call get_AIC_matrix(n_vort, X, Uinf, AIC)
    ! Get geometric vectors
    call get_geom_vectors(n_vort, X, alpha0, chords, Uai, Uni, si, areas)

    ! Get induced velocities
    call get_local_vels(n_vort, Gama, AIC, Vinf, Uai, Uni, Vlocal, alphaLocal)

    ! Get circulation forces
    call compute_circulation_forces(n_vort, rho, Vlocal, si, Gama, Fcirc)

    ! Get airfoil forces
    call compute_airfoil_forces(n_vort, rho, cl0, cla, areas, Vlocal, alphaLocal, si, Fairf)

    ! The residual is the difference between the forces for each horseshoe vortex
    do ii=1,n_vort

       ! The residual is the magnitude of the difference
       res(ii) = sum((Fcirc(:,ii) - Fairf(:,ii))**2)

    end do

  end subroutine get_residuals

  !============================================================

  subroutine get_functions(n_vort, X, Gama, alpha0, chords, cl0, cla, Vinf, rho, Sref, CL, CD, L, D)

    ! This function returns the total lift and drag of the body


    implicit none

    ! Input variables
    integer, intent(in) :: n_vort
    real, intent(in) :: X(3,n_vort+1)
    real, intent(in) :: Gama(n_vort), alpha0(n_vort)
    real, intent(in) :: chords(n_vort), cl0(n_vort)
    real, intent(in) :: cla(n_vort), Vinf(3), rho

    ! Output variables
    real, intent(out) :: Sref, CL, CD, L, D

    ! Working variables
    real :: Uinf(3), AIC(3,n_vort,n_vort)
    real :: Uai(3,n_vort), Uni(3,n_vort), si(3,n_vort)
    real :: areas(n_vort), Vlocal(3,n_vort), alphaLocal(n_vort)
    real :: Fcirc(3,n_vort), Fbody(3), Fbodym, Vinfm

    ! EXECUTION

    ! Normalize free-stream velocity
    call normalize(Vinf, Uinf)

    ! Get AIC matrix
    call get_AIC_matrix(n_vort, X, Uinf, AIC)

    ! Get geometric vectors
    call get_geom_vectors(n_vort, X, alpha0, chords, Uai, Uni, si, areas)

    ! Get induced velocities
    call get_local_vels(n_vort, Gama, AIC, Vinf, Uai, Uni, Vlocal, alphaLocal)

    ! Get forces on every bound vortex
    call compute_circulation_forces(n_vort, rho, Vlocal, si, Gama, Fcirc)

    ! Get total force in the system of coordinates of the body
    Fbody(1) = sum(Fcirc(1,:))
    Fbody(2) = sum(Fcirc(2,:))
    Fbody(3) = sum(Fcirc(3,:))

    ! The drag is aligned with the free-stream
    call dot(Fbody, Uinf, D)

    ! The rest is the lift force
    call norm(Fbody, Fbodym)
    L = sqrt(Fbodym*Fbodym - D*D)

    ! Use the total area as reference
    Sref = sum(areas)

    ! Make results non dimensional
    call norm(Vinf, Vinfm)
    CL = 2.0*L/Sref/Vinfm/Vinfm/rho
    CD = 2.0*D/Sref/Vinfm/Vinfm/rho

  end subroutine get_functions

  !============================================================

  subroutine norm(a, am)

    implicit none

    real, intent(in) :: a(3)
    real, intent(out) :: am

    am = sqrt(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))

  end subroutine norm

  !============================================================

  subroutine normalize(a, an)

    ! This subroutine normalizes a vector a to get the unitary vector an.

    implicit none

    real, intent(in) :: a(3)
    real, intent(out) :: an(3)
    real :: am

    call norm(a, am)
    an = a/am

  end subroutine normalize

  !============================================================

  subroutine dot(a, b, aDb)

    implicit none

    real, intent(in) :: a(3), b(3)
    real, intent(out) :: aDb

    aDb = sum(a*b)

  end subroutine dot

  !============================================================

  subroutine cross(a, b, aCb)

    implicit none

    real, intent(in) :: a(3), b(3)
    real, intent(out) :: aCb(3)

    aCb(1) = a(2)*b(3) - a(3)*b(2)
    aCb(2) = a(3)*b(1) - a(1)*b(3)
    aCb(3) = a(1)*b(2) - a(2)*b(1)

  end subroutine cross

  !============================================================

  subroutine tapenade_main(n_vort, X, Gama, alpha0, chords, cl0, cla, Vinf, rho, res, Sref, CL, CD, L, D)

    ! Dummy subroutine to allow differentiation of get_residuals and
    ! get_functions in a single Tapenade call.
    !
    ! INPUTS
    !
    ! n_vort: integer -> Number of horseshoe vortices.
    ! X: real(3,n_vort+1) -> Coordinates of the kinks of the horseshoe vortices.
    ! Gama: real(n_vort) -> Circulation intensity of each horseshoe vortex.
    ! alpha0: real(n_vort) -> Array of local incidence angles.
    ! chords: real(n_vort) -> Array of local chords.
    ! cl0: real(n_vort) -> cl0 of each 2D section
    ! cla: real(n_vort) -> lift curve slope of each 2D section [1/rad]
    ! Vinf: real(3) -> Free-stream velocity vector
    ! rho: real -> Air density
    !
    ! OUTPUTS
    !
    ! res: real(n_vort) -> Residuals of the LLT method
    ! Sref: real -> Area of all panels
    ! CL: real -> Lift coefficient (adimensionalized by Sref)
    ! CD: real -> Drag coefficient (adimensionalized by Sref)
    ! L: real -> Lift force
    ! D: real -> Drag force

    ! ADD YOUR CODE HERE

    ! Input variables
    integer, intent(in) :: n_vort
    real, intent(in) :: X(3,n_vort+1)
    real, intent(in) :: Gama(n_vort), alpha0(n_vort)
    real, intent(in) :: chords(n_vort), cl0(n_vort)
    real, intent(in) :: cla(n_vort), Vinf(3), rho

    ! Output variables
    real, intent(out) :: Sref, CL, CD, L, D, res(n_vort)

    call get_functions(n_vort, X, Gama, alpha0, chords, cl0, cla, Vinf, rho, Sref, CL, CD, L, D)
    call get_residuals(n_vort, X, Gama, alpha0, chords, cl0, cla, Vinf, rho, res)


  end subroutine tapenade_main

end module llt_module
