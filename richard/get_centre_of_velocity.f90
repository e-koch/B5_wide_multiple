   subroutine get_centre_of_velocity(nstars,v,m,masscov)


   implicit none
   integer, intent(in) :: nstars
   double precision, intent(in) :: v(1:nstars,1:3)
   double precision, intent(in) :: m(1:nstars)
   integer :: i
   double precision :: totalmass
   double precision, intent(out) :: masscov(1:3)

   write(6,*) 'get centre of mass subroutine',nstars,v(nstars,1:3)



   masscov=0.d0
   totalmass=0.d0
   DO i=1,nstars
      masscov(1)=masscov(1) + m(i)*v(i,1)
      totalmass=totalmass + m(i)
   END DO

   masscov=masscov/totalmass

   write(6,*) 'mass weighted cov',real(masscov)

   return
   end subroutine get_centre_of_velocity
