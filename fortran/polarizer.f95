      ! =============================================================
      ! == Michael Schneider 2011/07/12                            ==
      ! == This work is licensed under the Creative Commons        ==
      ! == Attribution-ShareAlike 3.0 Unported License.            ==
      ! == http://creativecommons.org/licenses/by-sa/3.0/          ==
      ! =============================================================

      ! purpose
      ! =======
      ! Calculate total reflectivity and degree of polarization for a
      ! x-ray beam going through a set of four mirrors using
      ! transfer-matrix-method. The mirrors are made of an arbitrary
      ! number of layers with complex refractive indices. The first layer is
      ! ambient layer, which will be vacuum in almost all cases. The last layer
      ! in the list is assumed to be substrate. 
      ! Program is intended for use with f2py, e.g.:
      ! # f2py2 -c -m polarizer polarizer.f95
      ! which compiles on my machine without complaints to polarizer.so
      ! usage in python then is something like:
      ! >>> import polarizer
      ! >>> R, dop, phi = polarizer.polarizer(idx, thickness, angles, energy)
      ! idx (no unit): 
      ! ==============
      ! list of complex refractive indices of mirror, starting with
      ! ambient (e.g. vacuum: 1+0j) and including substrate
      !
      ! thickness (nm):
      ! ===============
      ! corresponding list of layer thicknesses, SAME LENGTH AS IDX(!)
      ! values for ambient and substrate need to be valid reals, but are
      ! arbitrary otherwise
      ! 
      ! angles (rad):
      ! =============
      ! list of input angles
      !
      ! energy (eV):
      ! ============
      ! photon energy
      !
      ! output is total reflectivity (=transmission), degree of circular
      ! polarization and rotation angle around long axis of polarizer
      ! for equal intensities in s/p polarization of whole polarizer
      ! (FOUR MIRRORS)

      ! {{{ worker function: loop through layers, angles
      SUBROUTINE polarizer(idx, thickn, n, angles, nA, energy,&
                           R_tot, dop, phi)

      ! {{{ declarations
      IMPLICIT NONE
      
      ! List of arguments
      ! =================
      ! INPUT 
      !  n       : number of layers
      !  nA      : number of incidence angles to calculate
      !  idx     : list of complex refractive indices (no unit)
      !  thickn  : corresponding list of layer thicknesses [m]
      !  energy  : photon energy of incident beam [eV]
      !  angles  : list of incidence angles [rad]
      INTEGER, INTENT(IN)                 :: n, nA
      COMPLEX, DIMENSION(n), INTENT(IN)   :: idx
      REAL, DIMENSION(n), INTENT(IN)      :: thickn
      REAL, INTENT(IN)                    :: energy
      REAL, DIMENSION(nA), INTENT(IN)     :: angles
       
      ! OUTPUT
      !  R_total : total balanced reflectivity (same intensity in s/p-pol) 
      !  dop     : degree of polarization (no unit)
      !  phi     : rotation angle of mirrorsurface against p-pol for
      !            balanced reflectivity [rad]
      !  all output variables are arrays corresponding to list of
      !  incidence angles
      REAL, DIMENSION(nA), INTENT(OUT)  :: R_tot, dop, phi

      ! List of variables
      ! =================
      INTEGER           :: lA, ln    ! loop index for angle/layer iteration
      COMPLEX, DIMENSION(2, 2) :: M_s, M_p  ! transfermatrix (s/p)
      COMPLEX :: lambda, angle          ! wavelength [m], actual angle
      COMPLEX :: r_s, r_p            ! fresnel reflection coefficients
      REAL :: abs_R_s, abs_R_p       ! absolute reflectivity
      !
      ! List of parameters
      ! ==================
      REAL, PARAMETER :: pi = 3.141592653589793
      REAL, PARAMETER :: h = 6.62606896E-34      ! planck constant (Js)
      REAL, PARAMETER :: c = 299792458.0         ! speed of light (m/s)
      REAL, PARAMETER :: e = 1.602176487E-19     ! electron charge (C)
      !
      ! List of functions
      ! =================
      REAL :: cplx_phase             ! phase of complex number [rad]
      COMPLEX :: refracted              ! snells law: refracted angle [rad]

      ! }}}

      ! convert energy[eV] to lambda[m]
      lambda = COMPLEX(h * c / (e * energy), 0)

      ! WRITE(*,*) 'angles:', angles * 180 / pi
      DO lA = 1, nA
        angle = COMPLEX(angles(lA), 0)
        ! generate ambient layer
        CALL M_amb(angle, idx(1), 's', M_s)
        CALL M_amb(angle, idx(1), 'p', M_p)

        ! iterate through actual mirror layers, except substrate and ambient
        DO ln = 2, n - 1
          ! we need refracted values of wavelength and angle.
          angle = refracted(angle, idx(ln-1), idx(ln))
          CALL M_lay(idx(ln), thickn(ln), angle,&
                     lambda / idx(ln), 's', M_s)
          CALL M_lay(idx(ln), thickn(ln), angle,&
                     lambda / idx(ln), 'p', M_p)
        END DO

        ! last step: substrate
        angle = refracted(angle, idx(n-1), idx(n))
        CALL M_sub(idx(n), angle, 's', M_s)
        CALL M_sub(idx(n), angle, 'p', M_p)

        ! calculate the output values from complete matrix
        r_s =  M_s(2, 1) / M_s(1, 1)
        r_p =  M_p(2, 1) / M_p(1, 1)
        ! WRITE(*,*) r_s, r_p
        abs_R_s = CABS(r_s)**2   ! reflectivity of single mirror
        abs_R_p = CABS(r_p)**2   ! reflectivity of single mirror

        phi(lA) = ATAN(abs_R_s / abs_R_p)
        R_tot(lA) = SQRT((SIN(phi(lA)) * abs_R_s)**2 + (COS(phi(lA)) * abs_R_p**2))
        dop(lA) = SIN(4 * (cplx_phase(r_p) - cplx_phase(r_s)))
      END DO

      END SUBROUTINE
      ! }}}

      ! {{{ transfer matrices
      ! there are different matrices used, dependent on the position of the
      ! layer in the mirror-stack: ambient (M_amb), middle layers (M_lay) and
      ! substrate (M_sub). All depend on polarization (pol = 's' ∨ 'p')
      
      ! {{{ ambient: M_amb
      SUBROUTINE M_amb(angle, idx, pol, M)
      IMPLICIT NONE

      COMPLEX, INTENT(IN)     :: angle          ! angle of incidence [rad]
      CHARACTER, INTENT(IN)   :: pol            ! polarization (s/p)
      COMPLEX, INTENT(IN)     :: idx            ! index of refraction
      COMPLEX, DIMENSION(2,2), INTENT(OUT) :: M ! matrix of ambient layer

      ! cumbersome, but more readable: single matrix elements
      IF (pol == 'p') THEN
          M(1, 1) = .5
          M(1, 2) = .5 * COS(angle) / idx
          M(2, 1) = -.5
          M(2, 2) = .5 * COS(angle) / idx
      ELSE
          M(1, 1) = .5
          M(1, 2) = .5 / (COS(angle) * idx)
          M(2, 1) = .5
          M(2, 2) = -.5 / (COS(angle) * idx)
      END IF
      ! WRITE(*,*) 'ambient:', M(1,1), M(1,2), M(2,1), M(2,2)

      END SUBROUTINE
      ! }}}

      ! {{{ middle layer: M_lay
      ! get existing transfermatrix (from ambient or previous layers)
      ! and multiply with new matrix...recursion ftw
      SUBROUTINE M_lay(idx, thickn, angle, lambda, pol, M)
      IMPLICIT NONE
      ! input
      ! =====
      REAL, INTENT(IN) :: thickn                ! ([rad], [m], [m])
      ! Need to make sure that incidence angle is already refracted!
      ! Thickness means physical thickness, not optical!
      COMPLEX, INTENT(IN) :: idx, angle, lambda ! index of refraction
      CHARACTER, INTENT(IN) :: pol              ! polarization (s/p)
      ! output
      ! ======
      COMPLEX, DIMENSION(2,2), INTENT(INOUT) :: M
      ! local variables
      ! ==============
      COMPLEX, DIMENSION(2,2) :: M_tmp   ! to make matrix product more readable
      COMPLEX                 :: beta    ! beta factor
      !
      ! parameters
      ! ==========
      REAL, PARAMETER :: pi = 3.141592653589793
      COMPLEX, PARAMETER :: ii = (0, 1)          ! imaginary i

      ! cumbersome, but more readable: single matrix elements. 
      ! I'll use a temporary matrix (M_tmp) that holds the single layer and
      ! multiply this with the input matrix (M). Think it's easier to
      ! read/understand/debug/maintain this way.

      ! beta is a phase factor in the fresnel-equations
      beta = 2 * pi * thickn * idx * COS(angle) / lambda
      ! WRITE(*,*) 'input:', pi, thickn, idx, angle, lambda

      IF (pol == 'p') THEN
          M_tmp(1, 1) = COS(beta)
          M_tmp(1, 2) = ii * COS(angle) * SIN(beta) / idx
          M_tmp(2, 1) = ii * idx * SIN(beta) / COS(angle)
          M_tmp(2, 2) = COS(beta)
      ELSE
          M_tmp(1, 1) = COS(beta)
          M_tmp(1, 2) = ii * SIN(beta) / (idx * COS(angle))
          M_tmp(2, 1) = ii * idx * SIN(beta) * COS(angle)
          M_tmp(2, 2) = COS(beta)
      END IF

      ! calculate the new compound transfermatrix
      M = MATMUL(M, M_tmp)
      ! WRITE(*,*) 'layer:', M_tmp(1,1), M_tmp(1,2), M_tmp(2,1), M_tmp(2,2)

      END SUBROUTINE
      ! }}}

      ! {{{ substrate: M_sub
      SUBROUTINE M_sub(idx, angle, pol, M)
      IMPLICIT NONE

      ! input/output
      ! ============
      COMPLEX, INTENT(IN) :: idx                  ! index of refraction
      COMPLEX, INTENT(IN) :: angle                ! angle of incidence [rad]
      CHARACTER, INTENT(IN) :: pol                ! polarization (s/p)
      COMPLEX, DIMENSION(2,2), INTENT(INOUT) :: M ! matrix of ambient layer

      ! local variables
      ! ===============
      COMPLEX, DIMENSION(2, 2) :: M_tmp ! same as in M_lay

      ! cumbersome, but more readable: single matrix elements
      IF (pol == 'p') THEN
          M_tmp(1, 1) = COS(angle) / idx
          M_tmp(1, 2) = 0
          M_tmp(2, 1) = 1
          M_tmp(2, 2) = 0
      ELSE
          M_tmp(1, 1) = 1 / (COS(angle) * idx)
          M_tmp(1, 2) = 0
          M_tmp(2, 1) = 1
          M_tmp(2, 2) = 0
      END IF

      ! glue it together...
      M = MATMUL(M, M_tmp)

      END SUBROUTINE
      ! }}}

      ! }}}

      ! {{{ helper functions
      ! {{{ phase of complex number
      FUNCTION cplx_phase(c)
      COMPLEX, INTENT(IN) :: c
      REAL :: cplx_phase
      REAL, PARAMETER :: pi = 3.141592653589793

      IF (REAL(c) > 0) THEN
          cplx_phase = ATAN(AIMAG(c) / REAL(c))
      ELSEIF ((REAL(c) < 0) .AND. (AIMAG(c) >= 0)) THEN
          cplx_phase = ATAN(AIMAG(c) / REAL(c)) + pi
      ELSEIF ((REAL(c) < 0) .AND. (AIMAG(c) < 0)) THEN
          cplx_phase = ATAN(AIMAG(c) / REAL(c)) - pi
      END IF

      END FUNCTION
      ! }}}

      ! {{{ snells law: refraction
      ! purpose: calculate refracted angle from input angle and indices
      ! of refraction from two adjacent layers
      FUNCTION refracted(a_in, n1, n2)
      COMPLEX, INTENT(IN) :: a_in, n1, n2 ! angle in [rad]
      COMPLEX :: refracted                ! [rad]
      refracted = ASIN(n1 / n2 * SIN(a_in))
      ! WRITE(*,*) 'refracted=', refracted
      END FUNCTION
      ! }}}
      ! }}}

      ! {{{ test driver
      ! simple test facility to test the polarizer submodules
!       PROGRAM call_polarizer
!       USE physicalconstants
! 
!       IMPLICIT NONE
! 
!       INTEGER                 :: i, n, nA     ! #layers, #angles
!       COMPLEX, DIMENSION(4)   :: idx          ! refractive indices
!       REAL,    DIMENSION(4)   :: thickn       ! layer thickness
!       REAL,    DIMENSION(10)   :: angles, R_tot, dop, phi  ! output
!       REAL                    :: energy       ! beam energy
! 
!       n = 4             ! 4 layers incl. ambient and substrate
!       nA = 10           ! 3 angles to calculate
!       idx(1) = (1., 0.)                      ! vacuum
!       idx(2) = (0.9149908945, -0.018499583)  ! B4C
!       idx(3) = (0.7841914, -0.12803553)      ! Mo
!       idx(4) = (0.9455383904, -0.0341861285) ! SiO2
!       thickn = [1., 3.E-9, 50.E-9, 1.]       ! layer thicknesses
!       angles = (/(pi / 4. + .1 * i * pi / 4, i=1, 10)/) ! [rad]
!       energy = 60.                           ! eV
! 
!       WRITE(*,*) "calling polarizer..."
!       CALL polarizer(idx, thickn, n, angles, nA, energy, R_tot, dop, phi)
!       WRITE(*,*) "Phi = ", 180. * phi / pi
!       WRITE(*,*) "DoP = ", dop
!       WRITE(*,*) "R = ", R_tot

!       END PROGRAM
      ! }}}

      ! vim: foldmethod=marker
