   SUBROUTINE binding_energy(v1,v2,r1,r2,m1,m2,Ebind)

   USE phys_constants_mod
   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(in) :: v1(1:3)
   DOUBLE PRECISION, INTENT(in) :: v2(1:3)
   DOUBLE PRECISION, INTENT(in) :: r1(1:3)
   DOUBLE PRECISION, INTENT(in) :: r2(1:3)
   DOUBLE PRECISION, INTENT(in) :: m1,m2
   double precision :: ev1(1:3)
   double precision :: ev2(1:3)
   double precision :: er1(1:3)
   double precision :: er2(1:3)
   double precision :: em1,em2
   DOUBLE PRECISION :: magv
   DOUBLE PRECISION :: magr
   DOUBLE PRECISION, INTENT(out) :: Ebind
! 
! ======== END OF DECLARATIONS ===============

   ev1(1:3)=v1(1:3)
   ev2(1:3)=v2(1:3)
   er1(1:3)=r1(1:3)
   er2(1:3)=r2(1:3)
   em1=m1
   em2=m2
!   write(6,*) ev1,ev2,er1,er2,em1,em2
!=============================================================================
!                          CONVERSION OF UNITS
!=============================================================================
! Convert masses from N-body units (Solar masses) into MKS units (kilograms)
   em1=em1*msun
   em2=em2*msun
!
! Convert separation components from pc (N-body) to metres (MKS)
! From pc:
   er1=er1*pc
   er2=er2*pc
!
! Convert velocity components from kms^-1 to ms^-1
   ev1=ev1*1000.
   ev2=ev2*1000.
!
!=============================================================================
! Determine the velocity magnitude
   magv=SQRT((ev2(1) - ev1(1))**2. + (ev2(2) - ev1(2))**2. + (ev2(3) - ev1(3))**2.)
!
! Determine the separation magnitude
   magr=SQRT((er2(1) - er1(1))**2. + (er2(2) - er1(2))**2. + (er2(3) - er1(3))**2.)
!   WRITE(6,*) 'magr = ',magr
!   WRITE(6,*) 'Star 1 velocity components: ',ev1
!   WRITE(6,*) 'Star 2 velocity components: ',ev2
!   write(6,*) 'Star 1 separation components: ',er1
!   write(6,*) 'Star 2 separation components: ',er2
!   WRITE(6,*) 'masses = ',em1,em2
!   pause
!   WRITE(6,*) 'magv = ',magv
!   WRITE(6,*) 'Star 1 velocity components: ',ev1
!   WRITE(6,*) 'Star 2 velocity components: ',ev2


! Determine the energy of the binary system
   Ebind=em1*em2*(((magv**2.)/(2.*(em1 + em2))) - G/magr)
!   WRITE(6,*) 'Ebind = ',Ebind
   



   RETURN 
   END SUBROUTINE binding_energy
