   subroutine get_centre_of_mass(nstars,r,m,masscom)


   implicit none
   integer, intent(in) :: nstars
   double precision, intent(in) :: r(1:nstars,1:3)
   double precision, intent(in) :: m(1:nstars)
   integer :: i
   double precision :: totalmass
   double precision, intent(out) :: masscom(1:3)

   write(6,*) 'get centre of mass subroutine',nstars,r(nstars,1:3)

   masscom=0.d0
   totalmass=0.d0
   DO i=1,nstars
      masscom(1:2)=masscom(1:2) + m(i)*r(i,1:2)
      totalmass=totalmass + m(i)
   END DO

   masscom=masscom/totalmass

   write(6,*) 'mass weighted com',real(masscom)


   return
   end subroutine get_centre_of_mass
