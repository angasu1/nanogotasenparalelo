Module phys_cons
  use types
  use global_variables   
  use math_subroutines
  implicit none


contains
   
        subroutine props_assignment(atnom,i)
         character(len=2),intent(in)::atnom
         integer(ik),intent(in)::i


         select case (trim(atnom))
         case('H')
          att(i)%atnum=1
          att(i)%ma=1.0079_rk*amu2au
         case('He')
          att(i)%atnum=2
          att(i)%ma=4.0026_rk*amu2au
         case('F')
          att(i)%atnum=9
          att(i)%ma=18.9984032_rk*amu2au
         case('Cl')
          att(i)%atnum=17
          att(i)%ma=35.453_rk*amu2au
         case('Br')
          att(i)%atnum=35
          att(i)%ma=79.904_rk*amu2au
         case('C')                      
          att(i)%atnum=6                
          att(i)%ma=12.011_rk*amu2au    
         case('N')                      
          att(i)%atnum=7                
          att(i)%ma=14.007_rk*amu2au 
         case default
           write(*,*) 'Atom ',atnom,'not in the list of mass assignments'
         end select

                
        end subroutine props_assignment

  subroutine cal_cm(xcmass,xx) 
    real(rk),intent(out)::xcmass(3)
    real(rk),intent(in)::xx(3,n_atxmol)
    integer(ik)::i

    xcmass=0._rk
    do i=1,n_atxmol
       xcmass(1)=xcmass(1)+xx(1,i)*att(i)%ma        
       xcmass(2)=xcmass(2)+xx(2,i)*att(i)%ma        
       xcmass(3)=xcmass(3)+xx(3,i)*att(i)%ma        
    enddo
    xcmass=xcmass/mtot

    return

  end subroutine cal_cm



subroutine rotcons_calculation
  real(rk)::rotcons(3)
  integer(ik)::i



 !diffusion constants in atomic units
  diffccm = 1._rk/(2._rk*mtot) 
  if (nhe.gt.0) then
  mHe=att(n_at+1)%ma
  diffcHe = 1._rk/(2._rk*mHe) 
  dmassred=mHe*mtot/(mhe+mtot)
  diffcHered = 1._rk/(2._rk*dmassred)
  endif 


  call rot_constants


  !A=10.403_rk/har2cm
  !B=A!OJO 

! A=A/har2cm
! B=B/har2cm
! C=C/har2cm



return
End Subroutine rotcons_calculation

subroutine rot_constants
real(rk)::rotcons(3)
real(rk)::inm(3),xcm(3),xx(3,n_atxmol)
integer(ik)::i
  xx=coor*ar2bo
  call cal_cm(xcm,xx)
  do i=1,n_at
  xx(:,i)=xx(:,i)-xcm        
  enddo
  call in_moment(inm,xx)
  do i=1,3
  rotcons(i)=1.0_rk/2._rk/inm(i)
  enddo
   A=rotcons(1)*factorrotx
   B=rotcons(2)*factorroty
   C=rotcons(3)*factorrotz
end subroutine rot_constants

subroutine in_moment(inm,xx)
real(rk),intent(out)::inm(3)
real(rk),intent(in)::xx(n_atxmol,3)
real(rk)::inten(3,3),v(3,3)
integer(ik)::j,i
 call in_tensor(inten,xx)
 if (moltyp.le.2) then
  inm(1)=max(inten(1,1),tol)
  inm(2)=max(inten(2,2),tol)
  inm(3)=max(inten(3,3),tol)
  return
 endif
!!Here we call the diagonalization of the inertia tensor       
 call jacobi(inten,3,3,inm,v,j)  
 call eigsrt(inm,v,3,3)  
!!write(*,*) 'eigenvectors'
!! do j=1,3
!!    write(*,*) (v(j,i),i=1,3)
!!    write(*,*)
!! enddo
!!read*
return
end subroutine in_moment

subroutine in_tensor(inten,xx)
real(rk),intent(out)::inten(3,3)
real(rk),intent(in)::xx(3,n_atxmol)
integer(ik)::i,j
inten=0._rk
do i=1,n_atxmol
 inten(1,1)=inten(1,1)+att(i)%ma*(xx(2,i)**2+xx(3,i)**2)       
 inten(2,2)=inten(2,2)+att(i)%ma*(xx(1,i)**2+xx(3,i)**2)       
 inten(3,3)=inten(3,3)+att(i)%ma*(xx(1,i)**2+xx(2,i)**2)       
 inten(1,2)=inten(1,2)-att(i)%ma*xx(1,i)*xx(2,i)       
 inten(1,3)=inten(1,3)-att(i)%ma*xx(1,i)*xx(3,i)       
 inten(2,3)=inten(2,3)-att(i)%ma*xx(2,i)*xx(3,i)       
enddo
do i=1,3
 do j=i+1,3
 inten(j,i)=inten(i,j)
 enddo
enddo

return
end subroutine in_tensor

  
  subroutine calculate_amm_par 
  integer(ik)::i
  real(rk):: r,z

!         call  cal_cm(xn,x_coor)

!         do i=1,n_at
!                 if (at_type(i).eq.'N') then 
!                 xn=(x_coor(i,:)-xn)*ar2bo
!                 endif
!         enddo
!         do i=1,n_at
!                 if (at_type(i).eq.'N') cycle 
!                 r=dsqrt(dot_product(x_coor(i,:),x_coor(i,:)))    
!                 z=x_coor(i,3)
!         umb_ang=pi/2.0_rk+abs(dasin(z/r))
!         exit
!         enddo




  end subroutine calculate_amm_par 
! Subroutine physical_constants_calc
!     real (kind=8)::dmassA,dmassB,dmassC
!     real (kind=8)::B_wnx,B_wny,B_wnz,grstate,enm1,enm2
!     character(len=15),parameter::gr='ground state =',et1='e(k=+-1) =',et2='e(k=0) ='

!      open (unit=300,file=trim(outdir)//"/e_teor.dat")
!     !Masses
!     dmassA =1.00794_rk * convmass!Hidrogeno 
!     dmassB =12.01115_rk*convmass!Carbono 
!     dmassC =14.00674_rk*convmass!Nitrogeno
!     dmassHe = 4.002602_rk*convmass*factormasshe!Helio  

!     if (molname.eq.'prd') Then
!     dmassMol = (5._rk*(dmassA + dmassB)+dmassC)*factormassmol  
!     B_wnx=0.1927_rk*factorrotx
!     B_wny=0.2159_rk*factorroty
!     B_wnz=0.1018_rk*factorrotz
!     else if (molname.eq.'ben') then
!     dmassMol = 6._rk*(dmassA + dmassB)
!     B_wnx=0.188_rk*factorrotx
!     B_wny=0.188_rk*factorroty
!     B_wnz=0.0938_rk*factorrotz

!     else 
!     stop 'recognized options for molname are "prd" for pyridine and "ben" for&
!            bencene'
!     endif
!     
!         
!     Bx=B_wnx/har2cm
!     By=B_wny/har2cm
!     Bz=B_wnz/har2cm

!     !Calculation of the theoretical energies
!     grstate=0_rk

!     enm1=grstate+b_wnx*2._rk-(b_wnx-b_wnz)
!     enm2=grstate+b_wnx*2._rk

!     write (300,*) 'Theoretical energies for the rigid rotor:'
!     write (300,*)
!     write (300,"('molname=',1X,A)") molname
!     write (300,*)
!     write (300,"(A,F8.4)") gr,grstate
!     write (300,*) 
!     write (300,"(A,F8.4)") et1,enm1
!     write (300,"(A,F8.4)") et2,enm2 
!    
!      !diffusion constants in atomic units
!      diffcHe = 1._rk/(2._rk*dmassHe) 
!      diffccm = 1._rk/(2._rk*dmassMol) 

!      close(300)

!      end Subroutine physical_constants_calc 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        end Module phys_cons
