   program multiple_system

   use phys_constants_mod
   implicit none
   character*15 :: name(1:4)
   double precision :: mass(1:4)
   double precision :: vel(1:4,1:3)
   double precision :: pos(1:4,1:3)    
   integer :: i,j,k
   integer :: nstars
   double precision :: masscom(1:3),velcom(1:3)
   double precision :: Ekin(1:4),Epot(1:4)
   double precision :: totEkin,totEpot
   double precision :: rkj
   double precision :: sma,e,per,Ebind
   double precision :: rcom(1:3),vcom(1:3),redmass
   double precision :: valt_crit
 

   write(6,*) 'find multiple system properties'
   
   name(1)='B5-IRS1'
   name(2)='Blob_1'
   name(3)='Blob_2'
   name(4)='Blob_3'

! masses (in Msun)
!   mass(1:4)=0.1    ! original values
   mass(1)=0.1
   mass(2)=0.1086
   mass(3)=0.078
   mass(4)=0.09
!   mass(1:4)=0.2
   
! velocities (line of sight, km/s)
!   vel(1,1:3)=10.215   ! original values
!   vel(2,1:3)=10.496   ! original values
!   vel(3,1:3)=10.235   ! original values
!   vel(4,1:3)=10.341   ! original values
   vel(1,1:3)=10.21   
   vel(2,1:3)=10.43   
   vel(3,1:3)=10.23   
   vel(4,1:3)=10.30   


! positions (au)
!   pos(1,1)=0.      ! original values
!   pos(1,2)=0.      ! original values
!   pos(2,1)=8157.   ! original values  
!   pos(2,2)=7778.   ! original values
!   pos(3,1)=-280.   ! original values
!   pos(3,2)=3559.   ! original values
!   pos(4,1)=-3967.  ! original values
!   pos(4,2)=-3908.  ! original values
   pos(1,1)=0.
   pos(1,2)=0.
   pos(2,1)=-8251.1
   pos(2,2)=7934.9
   pos(3,1)=248.9
   pos(3,2)=3309.9
   pos(4,1)=3873.9
   pos(4,2)=-3315.1


! no information on z-direction
   pos(1:4,3)=0. 
  

   nstars=4

   do i=1,nstars
      write(6,*) name(i),mass(i),vel(i,1),pos(i,1:3)
   end do

! put in SI units
   mass=mass*msun
   pos=pos*au
   vel=vel*1000.
   
   call get_centre_of_mass(nstars,pos,mass,masscom)
   call get_centre_of_velocity(nstars,vel,mass,velcom)

   Write(6,*) 'work out total energy of system'

   Ekin=0.d0
   Epot=0.d0
   totEkin=0.d0
   totEpot=0.d0
   
   do j=1,nstars
!      if (j==4) cycle
      Ekin(j)=0.5*mass(j)*(vel(j,1) - velcom(1))**2.
      do k=1,nstars
         rkj=0.d0
         if (k==j) cycle
         rkj=sqrt((pos(j,1) - pos(k,1))**2. + (pos(j,2) - pos(k,2))**2.)
         Epot(j)=Epot(j) + (G*mass(j)*mass(k))/rkj
      end do
      Epot(j)=Epot(j)*(-1.)
      totEkin=totEkin + Ekin(j)
      totEpot=totEpot + Epot(j)
      write(6,*) 'K + T ',Ekin(j) + Epot(j), 'K',Ekin(j),'T',Epot(j)
   end do

   write(6,*) 'virial ratio',(-2.)*totEkin/totEpot

   write(6,*) 'work out binary orbits'
   write(6,*) 'orbit of ',name(1),'and ',name(3)
   call binary_orbit(vel(1,1:3),vel(3,1:3),pos(1,1:3),pos(3,1:3),mass(1),mass(3),sma,per,e,Ebind)
   write(6,*) 'binary orbit',sma,per,e,Ebind
   write(6,*) 
   write(6,*) 'now work out "triple" properties'
   rcom(1:3)=(pos(1,1:3)*mass(1) + pos(3,1:3)*mass(3))/(mass(1) + mass(3))
   vcom(1:3)=(vel(1,1:3)*mass(1) + vel(3,1:3)*mass(3))/(mass(1) + mass(3))
   redmass=mass(1)*mass(3)/(mass(1) + mass(3))

   write(6,*) 'orbit of ',name(1),'+',name(3), 'AND',name(2)
   call binary_orbit(vel(2,1:3),vcom(1:3),pos(2,1:3),rcom(1:3),mass(2),mass(1)+mass(3),sma,per,e,Ebind)
   valt_crit=3.*((1. + mass(2)/(mass(1)+mass(3)))**(2./3.))*((1. - e)**(-1./6.))*((7./4. - 1./2. - 1.)**(1./3.)) ! i = 0 degrees
   write(6,*) 'binary orbit',sma,per,e,Ebind,valt_crit
   write(6,*) 'now using reduced mass'
   call binary_orbit(vel(2,1:3),vcom(1:3),pos(2,1:3),rcom(1:3),mass(2),redmass,sma,per,e,Ebind)
   valt_crit=3.*((1. + mass(2)/(redmass))**(2./3.))*((1. - e)**(-1./6.))*((7./4. - 1./2. - 1.)**(1./3.))
   write(6,*) 'binary orbit',sma,per,e,Ebind,valt_crit
   write(6,*)

   write(6,*) 'orbit of ',name(1),'+',name(3), 'AND',name(4)
   call binary_orbit(vel(4,1:3),vcom(1:3),pos(4,1:3),rcom(1:3),mass(4),mass(1)+mass(3),sma,per,e,Ebind)
   valt_crit=3.*((1. + mass(4)/(mass(1)+mass(3)))**(2./3.))*((1. - e)**(-1./6.))*((7./4. - 1./2. - 1.)**(1./3.)) ! i = 0 degrees 
   valt_crit=3.*((1. + mass(4)/(mass(1)+mass(3)))**(2./3.))*((1. - e)**(-1./6.))*((7./4. - 1./4. - 0.5)**(1./3.)) ! i = 45 degrees
   write(6,*) 'binary orbit',sma,per,e,Ebind,valt_crit
   write(6,*) 'now using reduced mass'
   call binary_orbit(vel(4,1:3),vcom(1:3),pos(4,1:3),rcom(1:3),mass(4),redmass,sma,per,e,Ebind)
   valt_crit=3.*((1. + mass(4)/(redmass))**(2./3.))*((1. - e)**(-1./6.))*((7./4. - 1./2. - 1.)**(1./3.))
   write(6,*) 'binary orbit',sma,per,e,Ebind,valt_crit
    

   masscom=masscom/au
   pos=pos/au   
!   write(6,*) masscom

   velcom=velcom/1000.
   vel=vel/1000.
!   write(6,*) velcom
 

   CALL pgbegin(0,'?',1,1)
   CALL pgslw(2)
   CALL pgsch(1.4)
   CALL pgsvp(0.2,0.85,0.15,.95)
   call pgswin(-10000.,10000.,-10000.,10000.)
   CALL pgbox('BCNTS',0.0,0,'BCNVTS',0.0,0)
   CALL pgmtxt('L',4.,0.5,0.5,'y/au')
   CALL pgmtxt('B',2.8,0.5,0.5,'x/au')
   CALL pgslw(5)
   call pgsch(3.)

     DO i=1,nstars
        call pgsci(i+1)
        write(6,*) real(pos(i,1:2))       
        CALL pgpoint(1,real(pos(i,1)),real(pos(i,2)),17)
     END DO
     call pgsci(1)
     call pgpoint(1,real(masscom(1)),real(masscom(2)),2)

   CALL pgend()

   


   end program multiple_system
