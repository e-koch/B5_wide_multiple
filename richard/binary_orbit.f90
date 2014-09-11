   SUBROUTINE binary_orbit(v1,v2,r1,r2,m1,m2,sma,per,e,Etot)
!
   USE phys_constants_mod
   IMPLICIT NONE
   DOUBLE PRECISION, intent(out) :: Etot                    ! Total energy of system
   DOUBLE PRECISION :: Ltot                    ! Total angular momentum of system
   DOUBLE PRECISION :: Lxsq,Lysq,Lzsq          ! Angular momentum components
   DOUBLE PRECISION :: magv                    ! velocity magnitude
   DOUBLE PRECISION :: magr                    ! separation magnitude
   DOUBLE PRECISION, INTENT(inout) :: v1(1:3)  ! velocity components, star 1
   DOUBLE PRECISION, INTENT(inout) :: v2(1:3)  ! velocity components, star 2
   DOUBLE PRECISION, INTENT(inout) :: r1(1:3)  ! separation components, star 1
   DOUBLE PRECISION, INTENT(inout) :: r2(1:3)  ! separation components, star 2
   DOUBLE PRECISION, INTENT(inout):: m1,m2     ! primary and secondary masses
   DOUBLE PRECISION, INTENT(out) :: sma        ! semi-major axis of system
   DOUBLE PRECISION, INTENT(out) :: per        ! period of system
   double PRECISION, INTENT(out) :: e          ! eccentricity of system
 
!
! =-=-= END OF DECLARATIONS =-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=    
   write(6,*) v1
   write(6,*) v2
   write(6,*) r1
   write(6,*) r2
   write(6,*) m1
   write(6,*) m2
!
! ************** Set up constants ***************
   twopi=8.*ATAN(1.)
! **************** End Constants ****************
!=============================================================================
!                          CONVERSION OF UNITS
!=============================================================================
! Convert masses from N-body units (Solar masses) into MKS units (kilograms)
!   m1=m1*msun
!   m2=m2*msun
!
! Convert separation components from AU/pc (N-body) to metres (MKS)
! From pc:
!   r1=r1*pc
!   r2=r2*pc
! From AU:
!   r1=r1*au
!   r2=r2*au  
!
! Convert velocity components from kms^-1 to ms^-1
!   v1=v1*1000.
!   v2=v2*1000.
!=============================================================================
!
! Determine the velocity magnitude
   magv=SQRT((v2(1) - v1(1))**2. + (v2(2) - v1(2))**2. + (v2(3) - v1(3))**2.)
!
! Determine the separation magnitude
   magr=SQRT((r2(1) - r1(1))**2. + (r2(2) - r1(2))**2. + (r2(3) - r1(3))**2.)
!   
   WRITE(6,*) 'magr = ',magr,'in au',magr/au
   write(6,*) 'Star 1 separation components: ',r1
   write(6,*) 'Star 2 separation components: ',r2
   WRITE(6,*) 'masses = ',m1,m2
   WRITE(6,*) 'magv = ',magv
   WRITE(6,*) 'Star 1 velocity components: ',v1 
   WRITE(6,*) 'Star 2 velocity components: ',v2
!
! Determine the energy of the binary system
   Etot=m1*m2*(((magv**2.)/(2.*(m1 + m2))) - G/magr)
   WRITE(6,*) 'Ebind = ',Etot

! ====== TOTAL ANGULAR MOMENTUM =====
   Lxsq=(((r1(2) - r2(2))**2.)*((v1(3) - v2(3))**2.) + ((r1(3) - r2(3))**2.)*((v1(2) - v2(2))**2.) &
&  - ((2.*(r1(2) - r2(2))*(r1(3) - r2(3))*(v1(2) - v2(2))*(v1(3) - v2(3))))) 
!
   Lysq=(((r1(3) - r2(3))**2.)*((v1(1) - v2(1))**2.) + ((r1(1) - r2(1))**2.)*((v1(3) - v2(3))**2.) &
&  - ((2.*(r1(1) - r2(1))*(r1(3) - r2(3))*(v1(1) - v2(1))*(v1(3) - v2(3))))) 
!
   Lzsq=(((r1(1) - r2(1))**2.)*((v1(2) - v2(2))**2.) + ((r1(2) - r2(2))**2.)*((v1(1) - v2(1))**2.) &
&  - ((2.*(r1(1) - r2(1))*(r1(2) - r2(2))*(v1(1) - v2(1))*(v1(2) - v2(2))))) 
!
   WRITE(6,*) 'Angular momentum components: ', Lxsq,Lysq,Lzsq
! Now multiply the "previous angular momentum" by the reduced mass of the system
   Ltot=(m1*m2/(m1 + m2))*SQRT(Lxsq + Lysq + Lzsq)   
   !WRITE(6,*) 'Ltot = ',Ltot
!
! Determine the eccentricity of the system
   if (1. + ((2.*(m1 + m2)*(Ltot**2.)*Etot)/((G**2.)*(m1**3.)*(m2**3.))) < 0.) then
      e=1.d-6
      write(6,*) 'eccentricity forced to 0.000001'
      !pause
   else
      e=SQRT(1. + ((2.*(m1 + m2)*(Ltot**2.)*Etot)/((G**2.)*(m1**3.)*(m2**3.))))
   end if
!   WRITE(6,*) 'eccentricity = ',e
!
! Determine the semi-major axis of the system
   sma=(G*m1*m2)/(2.*ABS(Etot))
   !WRITE(6,*) 'semi-major axis = ',sma
   !sma=(G*m1*m2)/(2.*ABS(Etot))!*(sqrt(1. - e**2.)))
   !wRITE(6,*) 'semi-minor axis = ',sma
!
! Determine the period of the system
   per=twopi*SQRT((sma**3.)/(G*(m1 + m2)))
!   write(6,*) 'sma + period',sma,per
   
!=============================================================================
!                          CONVERSION OF UNITS
!=============================================================================
! Convert masses back to N-body units (Solar Masses)
!   m1=m1/msun
!   m2=m2/msun
!
! Convert separation components to AU/pc
! In pc:
!  r1=r1/pc
!   r2=r2/pc
! In AU:
!   r1=r1/au
!   r2=r2/au
!
! Convert velocity components to kms^-1
!   v1=v1/1000.
!    v2=v2/1000.
!
! Convert semi-major axis of system into AU/pc
   sma=sma/au
!   sma=sma/pc 
!
! Convert period of the systems from seconds to days or years
!   per=per/day
   per=per/yr
!===========================================================================
!  

   RETURN 
   END SUBROUTINE binary_orbit
