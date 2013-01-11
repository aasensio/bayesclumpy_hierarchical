module mcmc_class_hierarchical
use mcmc_class
use maths, only : alngam
implicit none

	type sed
		integer :: nFilters
		real(kind=8) :: noise
		real(kind=8), pointer :: time(:), delta(:), model(:), trend(:)
	end type sed
	
	type, EXTENDS (mcmc) :: mcmc_hierarchical

		integer :: nloops
		type(coronal_loop), pointer :: cloop(:)
	
 		contains
 			procedure :: evalPosterior => logPosteriorCoronalLoop
 			procedure :: initialValues => initialValuesCoronalLoop
	end type mcmc_hierarchical

	type polynomialCoef
		real(kind=8), pointer :: value(:)
	end type polynomialCoef
	
	real(kind=8), allocatable :: tauA_i(:), xi_i(:), lR_i(:), phi_i(:), A0_i(:), taud(:), P(:)
	
	type(polynomialCoef), pointer :: avec(:)

	contains

!------------------------------------------------------------------
! Set initial values for the parameters
!------------------------------------------------------------------
	subroutine initialValuesCoronalLoop(this)
	class(mcmc_hierarchical), intent(inout) :: this
	integer :: loop, i, j
	real(kind=8) :: value

		loop = 1
		
! tauA
		do i = 1, this%nloops
			this%pars(loop) = 300.d0 !400.d0 * randomu() + 100.d0
			loop = loop + 1
		enddo

! xi
		do i = 1, this%nloops
			this%pars(loop) = 15.d0
			loop = loop + 1
		enddo

! lR
		do i = 1, this%nloops
			this%pars(loop) = 0.2d0 !2.d0 * randomu()
			loop = loop + 1
		enddo

! phi
		do i = 1, this%nloops
			this%pars(loop) = 0.d0 !2.d0 * pi * randomu() - pi
			loop = loop + 1
		enddo

! A0
		do i = 1, this%nloops
			this%pars(loop) = 5.d0 !20.d0 * randomu()
			loop = loop + 1
		enddo

! alpha, beta, gamma, delta, eps, eta
		do i = 1, 4
			this%pars(loop) = 1.d0! randomu()
			loop = loop + 1
		enddo

		this%ls = 1.d0
		
	end subroutine initialValuesCoronalLoop

!------------------------------------------------------------------
! Evaluate the posterior. This function overrides the function in the parent
!------------------------------------------------------------------
	function logPosteriorCoronalLoop(this, outBounds) result (logP)
	class(mcmc_hierarchical), intent(in) :: this
	integer, intent(inout) :: outBounds
	real(kind=8) :: logP, alpha, beta, gamma, delta, eps, eta, nu
	integer :: i, j, loop, ierr

! The vector of trial parameters enters with the following order
! tauA(1..nloops)
! xi(1..nloops)
! lR(1..nloops)
! phi(1..nloops)
! A0(1..nloops)
! alpha, beta, gamma, delta
! a(1..nloops,1..ncoef_of_loop)

! Fill the vectors to simplify the notation and test for the boundaries
		outBounds = 0
		loop = 1
		do i = 1, this%nloops
			tauA_i(i) = this%trial(loop)
			loop = loop + 1
			if (tauA_i(i) < 1 .or. tauA_i(i) > 1500.d0) then
				outBounds = 1
				return
			endif
		enddo

		do i = 1, this%nloops
			xi_i(i) = this%trial(loop)
			loop = loop + 1
			if (xi_i(i) < 1 .or. xi_i(i) > 50.d0) then
				outBounds = 1
				return
			endif
		enddo

		do i = 1, this%nloops
			lR_i(i) = this%trial(loop)
			loop = loop + 1
			if (lR_i(i) < 0 .or. lR_i(i) > 2) then
				outBounds = 1
				return
			endif
		enddo

		do i = 1, this%nloops
			phi_i(i) = this%trial(loop)
			loop = loop + 1
			if (phi_i(i) < -pi .or. phi_i(i) > pi) then
				outBounds = 1
				return
			endif
		enddo

		do i = 1, this%nloops
			A0_i(i) = this%trial(loop)
			loop = loop + 1
			if (A0_i(i) < 0) then
				outBounds = 1
				return
			endif			
		enddo

! Hyperparameters
		alpha = this%trial(loop)
		loop = loop + 1
		if (alpha < 0) then
			outBounds = 1
			return
		endif
		
		beta = this%trial(loop)
		loop = loop + 1
		if (beta < 0) then
			outBounds = 1
			return
		endif
		
		gamma = this%trial(loop)
		loop = loop + 1
		if (gamma < 0) then
			outBounds = 1
			return
		endif
		
		delta = this%trial(loop)
		loop = loop + 1
		if (delta < 0) then
			outBounds = 1
			return
		endif
		
! 		eps = this%trial(loop)
! 		loop = loop + 1
! 		if (eps < 0 .or. eps > 0.02) then
! 			outBounds = 1
! 			return
! 		endif
! 		
! 		eta = this%trial(loop)
! 		loop = loop + 1
! 		if (eta < 0 .or. eta > 0.02) then
! 			outBounds = 1
! 			return
! 		endif

 		logP = 0.d0
 		
!-----------------
! LOG-PRIORS
!-----------------
		nu = 0.001d0
		
! tau_A ~ IG(gamma,delta)
		do i = 1, this%nloops
 			logP = logP - (gamma+1.d0)*log(tauA_i(i)) - delta/tauA_i(i) + gamma*log(delta) - alngam(gamma, ierr)
		enddo

! xi ~ IG(0.001,0.001)
		do i = 1, this%nloops
			logP = logP - (nu+1.d0) * log(xi_i(i)) - nu / xi_i(i)
		enddo

! l/R ~ Beta(alpha,beta,0,2)
		do i = 1, this%nloops
			logP = logP + (alpha-1.d0)*log(lR_i(i)) + (beta-1.d0)*log(2.d0-lR_i(i)) - (alpha+beta-1.d0)*log(2.d0) - &
				(alngam(alpha, ierr) + alngam(beta, ierr) - alngam(alpha+beta, ierr))
		enddo

! A0 ~ IG(0.001,0.001)
		do i = 1, this%nloops
			logP = logP - (nu+1.d0) * log(A0_i(i)) - nu / A0_i(i)
		enddo

! Priors for alpha, beta, gamma, delta, eps and eta
! ----- Beta prior
! alpha, beta -> shape parameters -> Uniform over [0,infinity]

! ----- IG priors over tau_A
! gamma -> shape parameter -> Uniform prior over [0,infinity]
! delta -> scale parameter -> Jeffreys prior (Gamma(nu,nu) with nu<<1)
 		logP = logP - (nu-1.d0) * log(delta) - nu * delta


! ----- IG priors over xi
! eps -> shape parameter -> Uniform prior over [0,infinity]
! eta -> scale parameter -> Jeffreys prior (Gamma(nu,nu) with nu<<1)
!  		logP = logP - (nu-1.d0) * log(eta) - nu * eta
 		

!-----------------
! DATA LOG-LIKELIHOOD
!-----------------
		P = tauA_i * sqrt(2.d0) * sqrt( (xi_i+1.d0) / xi_i )
		taud = P * 2.d0 / pi * (xi_i+1.d0) / (xi_i-1.d0) / lR_i
		
! Model prediction for the current value of the parameters
 		do i = 1, this%nloops

! First the oscillatory part
  			this%cloop(i)%model = this%cloop(i)%trend + A0_i(i) * sin(2.d0*pi/P(i) * this%cloop(i)%time + phi_i(i)) * exp(-this%cloop(i)%time / taud(i))

! And now compute the log-likelihood
 			logP = logP - 0.5d0*sum( (this%cloop(i)%model - this%cloop(i)%delta)**2 / (this%cloop(i)%noise**2) )
 		enddo

		return

	end function logPosteriorCoronalLoop

end module mcmc_class_hierarchical