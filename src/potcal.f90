    Module pot_calculation
      use types
      use global_variables   
      use pothx
      use potinter
      use math_subroutines
      use hcldimpot
      Implicit none

      Contains

      subroutine potcalc(y1,y2,y3,x,z,cthet,VV)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        !Subroutine that calculates the potential for different molecules
        !it gets the cartesian coordinates for the centerofmass-He distance in A. and
        !it returns the potential in cm-1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        real(rk),intent(in)::x(3),z(3),y1(3),y2(3),y3(3)
        real(rk),intent(out)::VV
        real(rk)::rvec(3),r,cthet,Vsum 
        integer(ik)::indx  
 
          
       Select case (trim(molname))
         case('pyridine') 
        !xxvar=(/xvar,yvar,zvar/)
         !v=PyridineHePot(xx)
         case('benzene')
         !call potential_benzene(x,y,z,v)
         case('ammonia') 
         !v=potw(x*ar2bo-xn(1),y*ar2bo-xn(2),z*ar2bo-xn(3),umb_ang)*har2cm
         !v=potw(x,y,z,umb_ang)  !OJO
         case('hcn') 
         rvec=z-x                       
         r=dsqrt(dot_product(rvec,rvec))
         cthet=dot_product(rvec,y3)/r   
         indx=4                         
         if(ptyp.eq.1)then              
           VV = V(r,cthet,indx)

         else if(ptyp.eq.2)then         
           !cthet=dacos(cthet)  
           call interpol(Vsum,r,cthet)  
           VV = Vsum                    
         endif                    
         case('hcl','hf','hfcl','hbr')
         rvec=z-x                       
         r=dsqrt(dot_product(rvec,rvec))
         cthet=dot_product(rvec,y3)/r ! si cthet es el cos(theta) entonces devo convertirlo a solo tener theta porque en mis listas por cada theta tengo valor
         indx=1                         
         if(ptyp.eq.1)then              
           VV = V(r,cthet,indx)         
         else if(ptyp.eq.2)then         
           !cthet=dacos(cthet)  
           call interpol(Vsum,r,cthet)  
           VV = Vsum                    
         endif  
       end select



      end subroutine potcalc

      subroutine potparam
       !Load the parameters for the corresponding potential

       Select case (trim(molname))
         case('pyridine') 
         !xxvar=(/xvar,yvar,zvar/)
         !v=PyridineHePot(xx)
         case('benzene')
         !call potential_benzene(x,y,z,v)
         case('ammonia') 
         !v=potw(x*ar2bo-xn(1),y*ar2bo-xn(2),z*ar2bo-xn(3),umb_ang)*har2cm
         !v=potw(x,y,z,umb_ang)  !OJO
         case('hf')
         call pars_assignment(1)
         case('hfcl')
         call pars_assignment(1)
         case('hcl')
         call pars_assignment(2)
         case('hbr')
         call pars_assignment(3)
         case('hcn')
         call pars_assignment(4)
       end select


      end subroutine potparam

      subroutine potdimparam
       !Load the parameters for the corresponding potential

       Select case (trim(molname))
         case('pyridine') 
         !xxvar=(/xvar,yvar,zvar/)
         !v=PyridineHePot(xx)
         case('benzene')
         !call potential_benzene(x,y,z,v)
         case('ammonia') 
         !v=potw(x*ar2bo-xn(1),y*ar2bo-xn(2),z*ar2bo-xn(3),umb_ang)*har2cm
         !v=potw(x,y,z,umb_ang)  !OJO
         case('hf')
         case('hfcl')
         call potkass
         case('hcl')
         call potkass
         case('hbr')
         case('hcn')
       end select


      end subroutine potdimparam

     subroutine potdimer(j,ri,rj,xcm1,xcm2,ctheta1c,ctheta2c,phic,potdim)
      integer(ik),intent(in)::j
      real(rk),intent(in)::ri(3),rj(3),xcm1(3),xcm2(3)
      real(rk)::potdim,ctheta1,ctheta2,phi
      real(rk)::cart_coor(3,2*nmon),error,rvec(3),rad
      integer(ik)::i
      real(rk)::ctheta1c,ctheta2c,phic,cart_coorc(3,2*nmon)
      real(rk)::ric(3),rjc(3),rvecc(3),radc
      real(rk)::xa,ya,za,xb,yb,zb,xcoma,ycoma,zcoma,xcomb,ycomb,zcomb
      real(rk)::tha2,thb2,pha2,phb2,ph


      !Subroutine to convert from dmc coordinates to internal and cartesian
      call internal_coord_conv(ri,rj,xcm1,xcm2,rvec,rad,ctheta1,ctheta2,phi,cart_coor)
      !Another method for converting coordinates
      ric=ri;rjc=rj;rvecc=rvec;radc=rad
      call calculated_angles(ric,rjc,rvecc,radc,ctheta1c,ctheta2c,phic,cart_coorc)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !David`s conversion of coordinates
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! xa=ri(1)+xcm1(1)
    ! ya=ri(2)+xcm1(2)
    ! za=ri(3)+xcm1(3)
    ! xb=rj(1)+xcm2(1)
    ! yb=rj(2)+xcm2(2)
    ! zb=rj(3)+xcm2(3)
    ! xcoma=xcm1(1)
    ! ycoma=xcm1(2)
    ! zcoma=xcm1(3)
    ! xcomb=xcm2(1)
    ! ycomb=xcm2(2)
    ! zcomb=xcm2(3)
    !    call angs(xa,ya,za,xb,yb,zb,xcoma,ycoma,zcoma,xcomb,ycomb,zcomb,j,&
    !   &tha2,thb2,pha2,phb2,ph)
    !  ctheta1c=dcos(tha2)
    !  ctheta2c=dcos(thb2)
    !  phic=ph

      !!!!!!!!!!!!Checking part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !if (abs(norm(ri)-1.0_rk).gt.1.0e-10_rk) stop 'norma1!!'
     !if (abs(norm(rj)-1.0_rk).gt.1.0e-10_rk) stop 'norma2!!'

     ! if  (abs(ctheta1-ctheta1c).gt.tol) stop 'theta1'
     ! if  (abs(ctheta2-ctheta2c).gt.tol) stop 'theta2'

     !if(abs(deg(phi)-deg(phic)).lt.tolb) then
     !  error=0.0_rk
     !else
     !  error=abs(deg(phi)-360.0_rk+deg(phic))
     !endif

     ! If (counter.lt.100) then
     !        counter=counter+1
     !        if (counter.eq.99) stop 'counter'
     !        write(234,1000) xcm1
     !        write(234,1000) xcm2
     !        write(234,1000) ri
     !        write(234,1000) rj
     !        write(234,1000) deg(dacos(ctheta1)),deg(dacos(ctheta2)),deg(phi),rvec,rad
     !        write(234,1000) deg(dacos(ctheta1c)),deg(dacos(ctheta2c)),deg(phic),rvecc,radc
     !        write(234,*)
!1000 format (*(f18.10,3x))
     ! endif

      !!!!!!!!!!!!!!!!End of checking!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       Select case (trim(molname))
         case('pyridine') 
         !xxvar=(/xvar,yvar,zvar/)G
         !v=PyridineHePot(xx)
         case('benzene')
         !call potential_benzene(x,y,z,v)
         case('ammonia') 
         !v=potw(x*ar2bo-xn(1),y*ar2bo-xn(2),z*ar2bo-xn(3),umb_ang)*har2cm
         !v=potw(x,y,z,umb_ang)  !OJO
         case('hf')
          call hf2pot(1.7374_rk,1.7374_rk,rad,deg(dacos(ctheta1c)),deg(dacos(ctheta2c)),deg(phic),potdim)
         !call hf2pot(1.7445d0, 1.7404d0, 5.144d0, 9.00d0, 64.14d0, 180.00d0,potdim)
         !write(*,*) 'atila',potdim
         case('hfcl')
         call pothcl2(rad*bo2ar,ctheta1c,ctheta2c,phic,potdim)
         case('hcl')
         call pothcl2(rad*bo2ar,ctheta1c,ctheta2c,phic,potdim)
        !call pothcl2(rad,ctheta1,ctheta2,phi,potdim)
           !write(*,*) 'rad',rad,ctheta1,ctheta2,phi,potdim
           !stop
         case('hbr')
         case('hcn')
       end select



            return
            end subroutine potdimer
       
           subroutine calculated_angles(ri,rj,rij,rad,ct1,ct2,ph,xx)
           real(rk),intent(in)::ri(3),rj(3),rij(3),rad
           real(rk)::ric(3),rjc(3),rijc(3)
           real(rk)::xx(3,2*nmon),ct1,ct2,ph,ph1,ph2,fac1,fac2
           real(rk)::phrot,throt

           fac1=-(att(1)%ma*requil/mtot)
           fac2=fac1+requil
           ric=ri
           rjc=rj
           rijc=rij

     !Rotation to the system
             phrot=datan2(-rijc(2),rijc(1))
             throt=dacos(rijc(3)/norm(rijc))
             throt=pi/2.0_rk-throt

            !write(*,*) 'rij ant',rij
             call rotz(rijc,-phrot)
             call rotz(ric,-phrot)
             call rotz(rjc,-phrot)
            !write(*,*) 'rij desp',rij
             call roty(rijc,throt)
             call roty(ric,throt)
             call roty(rjc,throt)
            !write(*,*) 'rij desp2',rij
            !stop



     !Calculation of cartesian coords
       xx(:,1)=fac2*ric-rijc/2.0_rk
       xx(:,2)=fac1*ric-rijc/2.0_rk
       xx(:,3)=fac2*rjc+rijc/2.0_rk
       xx(:,4)=fac1*rjc+rijc/2.0_rk


     !Angles for potential calculation
      ct1=dot_product(ric,rijc)/norm(rijc)
      ct2=dot_product(rjc,rijc)/norm(rijc)

      !phi, method 1
     !temp1=cross_unitary(ri,rij)
     !temp2=cross_unitary(rj,rij)
     !cphi=dot_product(temp1,temp2)
     !ph=dacos(cphi)

      !Phi method 2
      ph1=datan2(ric(2),ric(3))
      ph2=datan2(rjc(2),rjc(3))
      ph=abs(ph1-ph2)
                   



           return
           end subroutine calculated_angles

      subroutine internal_coord_conv(ri,rj,xcm1,xcm2,rvec,rad,ctheta1,ctheta2,phi,cart_coor)
      real(rk),intent(in)::ri(3),rj(3),xcm1(3),xcm2(3)
      real(rk),intent(out)::rvec(3),rad,ctheta1,ctheta2,phi,cart_coor(3,2*nmon)      
      real(rk)::fac1,fac2,temp1(3),temp2(3),cphi

      rvec=xcm1-xcm2  
      rad=norm(rvec)
      ctheta1=dot_product(rvec,ri)/rad
      ctheta2=dot_product(rvec,rj)/rad
      temp1=cross_unitary(ri,rvec)
      temp2=cross_unitary(rj,rvec)
      cphi=dot_product(temp1,temp2)
      phi=dacos(cphi)

      fac1=-(att(2)%ma*requil/mtot)
      fac2=fac1+requil

       cart_coor(:,1)=fac2*ri+xcm1
       cart_coor(:,2)=fac1*ri+xcm1
       cart_coor(:,3)=fac2*rj+xcm2
       cart_coor(:,4)=fac1*rj+xcm2

      end subroutine internal_coord_conv


      subroutine angs(xa,ya,za,xb,yb,zb,xcoma,ycoma,zcoma,xcomb,ycomb,zcomb,jj,&
        &tha2,thb2,pha2,phb2,ph)
	real(rk):: xA,yA,zA
	real(rk):: xcomA,ycomA,zcomA
	real(rk):: xB,yB,zB
	real(rk):: xcomB,ycomB,zcomB
        real(rk):: xa0,ya0,za0,xb0,yb0,zb0,xxa,yya,zza,xxb,yyb,zzb
        real(rk):: rr,xnewa,ynewa,znewa,xnewb,ynewb,znewb,xcmb,ycmb,zcmb
        real(rk):: xn1a,xn1b,yn1a,yn1b,zn1a,zn1b,xcm1b,ycm1b,zcm1b
        real(rk):: tha2,thb2,pha2,phb2,thet,phi,th,ph
        integer(ik)::j,jj

        j=jj

 xa0=0.d0
 ya0=0.d0
 za0=0.d0


 xb0=xcomB-xcomA
 yb0=ycomB-ycomA
 zb0=zcomB-zcomA

       xxa=xA-xcomA
               yya=yA-ycomA
               zza=zA-zcomA

               xxb=xB-xcomA
               yyb=yB-ycomA
               zzb=zB-zcomA



 rr=dsqrt((xa0-yb0)**2+(ya0-yb0)**2+(za0-zb0)**2)

 thet=dacos(zzb/rr)
 phi=datan2(yb0,xb0)



 xnewa=xxa*dcos(phi)+yya*dsin(phi)
 ynewa=-xxa*dsin(phi)+yya*dcos(phi)
 znewa=zza

 xnewb=xxb*dcos(phi)+yyb*dsin(phi)
 ynewb=-xxb*dsin(phi)+yyb*dcos(phi)
 znewb=zzb

 xcmb=xb0*dcos(phi)+yb0*dsin(phi)
 ycmb=-xb0*dsin(phi)+yb0*dcos(phi)
 zcmb=zb0



 rr=dsqrt(xcmb**2+ycmb**2+zcmb**2)

 th=dacos(zcmb/rr)
 

 xn1a=dcos(th)*xnewa-dsin(th)*znewa
 yn1a=ynewa
 zn1a=dsin(th)*xnewa+dcos(th)*znewa

 xn1b=dcos(th)*xnewb-dsin(th)*znewb
 yn1b=ynewb
 zn1b=dsin(th)*xnewb+dcos(th)*znewb

 xcm1b=dcos(th)*xcmb-dsin(th)*zcmb
 ycm1b=ycmb
 zcm1b=dsin(th)*xcmb+dcos(th)*zcmb

 
 tha2=dacos(zn1a)
 thb2=dacos(zn1b-zcm1b)

 pha2=datan2(yn1a,xn1a)
 phb2=datan2((yn1b-ycm1b),(xn1b-xcm1b))
 ph=pha2-phb2


 


return
end subroutine angs

           subroutine cartesian_coords_conv(ri,rj,rij,rad,xx,xxmon,zhe,rhe,thhe)
           real(rk),intent(in)::ri(3),rj(3),rij(3),rad
           real(rk)::xx(3,2*nmon+nhe),xxmon(3,nmon),zhe(3,nhe),xcm1(3),xcm2(3)
           real(rk)::phrot,throt,fac1,fac2,rhe(nmon,nhe),thhe(nmon,nhe)
           real(rk)::transl1(3),transl2(3)
           integer(ik)::i,j
           character(len=2)::atip(2*nmon+nhe)
           real(rk)::th1com(3),th2com(3),taucom(3),Rcom(3),rhecom(3,2,nhe),alhecom(3,2,nhe)
           real(rk)::rijcom(3),ricom(3),rjcom(3),xmoncom(3,2),zhecom(3,nhe)
           real(rk)::ph1,ph2,ph

          !write(*,*) 'cartesian'
           

           fac2=-(att(1)%ma*requil*ar2bo/mtot)
           fac1=fac2+requil*ar2bo


           if (nmon.eq.2) then
            !!!!!!!!!!!!!!!!!DEBUGGING
            atip(1:4)=(/'H ','Cl','H ','Cl'/)
            do j=1,nhe
                    atip(4+j)='He'
            enddo
           !xx(:,1)=(xxmon(:,1)+fac1*ri)*bo2ar
           !xx(:,2)=(xxmon(:,1)+fac2*ri)*bo2ar
           !xx(:,3)=(xxmon(:,2)+fac1*rj)*bo2ar
           !xx(:,4)=(xxmon(:,2)+fac2*rj)*bo2ar
           !xx(:,5)=zhe(:,1)*bo2ar
           !call write_xyz(238,5,atip,xx)


            !Calculos comprobacion 1
           !ricom=ri
           !rjcom=rj
           !xmoncom(:,1)=xx(:,1)*ar2bo-fac1*ricom
           !xmoncom(:,2)=xx(:,3)*ar2bo-fac1*rjcom
           !rijcom=xmoncom(:,2)-xmoncom(:,1)
           !rcom(1)=norm(rijcom)
           !zhecom=zhe(:,1)
           !th1com(1)=dacos(dot_product(ricom,rijcom)/rcom(1))*180.0_rk/pi
           !th2com(1)=dacos(dot_product(rjcom,rijcom)/rcom(1))*180.0_rk/pi
           !rhecom(1,1)=norm(zhecom-xmoncom(:,1))
           !rhecom(1,2)=norm(zhecom-xmoncom(:,2))
           !alhecom(1,1)=dacos(dot_product(zhecom-xmoncom(:,1),ricom)/rhecom(1,1))*180_rk/pi
           !alhecom(1,2)=dacos(dot_product(zhecom-xmoncom(:,2),rjcom)/rhecom(1,2))*180_rk/pi
           !write(239,1000) rcom(1),th1com(1),th2com(1),rhecom(1,1),rhecom(1,2),alhecom(1,1),alhecom(1,2)
1000     format(*(F18.10,3x))            


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






           !Translation of all the coordinates
           transl1=xxmon(:,1)
           xxmon(:,1)=xxmon(:,1)-transl1
           xxmon(:,2)=xxmon(:,2)-transl1
           do j=1,nhe
           zhe(:,j)=zhe(:,j)-transl1
           enddo
           !!!!!!!!!!!!!!DEbbuging
            xx(:,1)=(xxmon(:,1)+fac1*ri)*bo2ar
            xx(:,2)=(xxmon(:,1)+fac2*ri)*bo2ar
            xx(:,3)=(xxmon(:,2)+fac1*rj)*bo2ar
            xx(:,4)=(xxmon(:,2)+fac2*rj)*bo2ar
            do j=1,nhe
            xx(:,4+j)=zhe(:,j)*bo2ar
            enddo
          ! call write_xyz(238,5,atip,xx)

            !Calculos comprobacion 2
            ricom=ri
            rjcom=rj
            xmoncom(:,1)=xx(:,1)*ar2bo-fac1*ricom
            xmoncom(:,2)=xx(:,3)*ar2bo-fac1*rjcom
            rijcom=xmoncom(:,2)-xmoncom(:,1)
            rcom(2)=norm(rijcom)
            th1com(2)=dacos(dot_product(ricom,rijcom)/rcom(2))*180.0_rk/pi
            th2com(2)=dacos(dot_product(rjcom,rijcom)/rcom(2))*180.0_rk/pi
            do j=1,nhe
            zhecom(:,j)=zhe(:,j)
            rhecom(2,1,j)=norm(zhecom(:,j)-xmoncom(:,1))
            rhecom(2,2,j)=norm(zhecom(:,j)-xmoncom(:,2))
            alhecom(2,1,j)=dacos(dot_product(zhecom(:,j)-xmoncom(:,1),&
                           &ricom)/rhecom(2,1,j))*180_rk/pi
            alhecom(2,2,j)=dacos(dot_product(zhecom(:,j)-xmoncom(:,2),&
                           &rjcom)/rhecom(2,2,j))*180_rk/pi
            enddo
          ! write(239,1000) rcom(2),th1com(2),th2com(2),rhecom(2,1),rhecom(2,2),alhecom(2,1),alhecom(2,2)



           !Rotation to the system
            phrot=datan2(-rij(2),rij(1))
            throt=dacos(rij(3)/norm(rij))
            throt=pi/2.0_rk-throt

            !write(*,*) 'rij ant',rij
            !write(*,*) 'rij ant',xxmon(:,1)
            !write(*,*) 'rij ant',xxmon(:,2)
             call rotz(rij,-phrot)
             call rotz(ri,-phrot)
             call rotz(rj,-phrot)
             call rotz(xxmon(:,1),-phrot)
             call rotz(xxmon(:,2),-phrot)
             do j=1,nhe
             call rotz(zhe(:,j),-phrot)
             enddo
            !write(*,*) 'rij ant2',rij
            !write(*,*) 'rij ant2',xxmon(:,1)
            !write(*,*) 'rij ant2',xxmon(:,2)
             call roty(rij,throt)
             call roty(ri,throt)
             call roty(rj,throt)
             call roty(xxmon(:,1),throt)
             call roty(xxmon(:,2),throt)
             do j=1,nhe
             call roty(zhe(:,j),throt)
             enddo
            !write(*,*) 'rij ant3',rij
            !write(*,*) 'rij ant3',xxmon(:,1)
            !write(*,*) 'rij ant3',xxmon(:,2)
            !stop

           !do j=1,3
           !xx(j,1)=(xxmon(j,1)+fac1*ri(j))*bo2ar
           !xx(j,2)=(xxmon(j,1)+fac2*ri(j))*bo2ar
           !xx(j,3)=(xxmon(j,2)+fac1*rj(j))*bo2ar
           !xx(j,4)=(xxmon(j,2)+fac2*rj(j))*bo2ar
           !xx(j,5)=zhe(j,1)*bo2ar
           !enddo
           ph1=datan2(ri(2),ri(3))
           ph2=datan2(rj(2),rj(3))
           ph=abs(ph1-ph2)
            xx(:,1)=(xxmon(:,1)+fac1*ri)*bo2ar
            xx(:,2)=(xxmon(:,1)+fac2*ri)*bo2ar
            xx(:,3)=(xxmon(:,2)+fac1*rj)*bo2ar
            xx(:,4)=(xxmon(:,2)+fac2*rj)*bo2ar
            do j=1,nhe
            xx(:,4+j)=zhe(:,j)*bo2ar
            enddo
          ! call write_xyz(238,5,atip,xx)




            !Calculos de comprobacion 3
            ricom=ri
            rjcom=rj
            xmoncom(:,1)=xx(:,1)*ar2bo-fac1*ricom
            xmoncom(:,2)=xx(:,3)*ar2bo-fac1*rjcom
            rijcom=xmoncom(:,2)-xmoncom(:,1)
            rcom(3)=norm(rijcom)
            th1com(3)=dacos(dot_product(ricom,rijcom)/rcom(3))*180.0_rk/pi
            th2com(3)=dacos(dot_product(rjcom,rijcom)/rcom(3))*180.0_rk/pi
            do j=1,nhe
            zhecom(:,j)=zhe(:,j)
            rhecom(3,1,j)=norm(zhecom(:,j)-xmoncom(:,1))
            rhecom(3,2,j)=norm(zhecom(:,j)-xmoncom(:,2))
            alhecom(3,1,j)=dacos(dot_product(zhecom(:,j)-xmoncom(:,1),&
                           &ricom)/rhecom(3,1,j))*180_rk/pi
            alhecom(3,2,j)=dacos(dot_product(zhecom(:,j)-xmoncom(:,2),&
                           &rjcom)/rhecom(3,2,j))*180_rk/pi
            enddo

            write(250,1000) rcom(3),th1com(3),th2com(3),ph,&
                  &(rhecom(3,1,j),rhecom(3,2,j),alhecom(3,1,j),alhecom(3,2,j),j=1,nhe)


      !xcm1=-rij/2.0_rk
      !xcm2=rij/2.0_rk


       else

       xx(:,1)=fac2*ri
       xx(:,2)=fac1*ri
       xcm1=xxmon(:,1)

       do i=1,nhe
        xx(:,2+i)=zhe(:,i)-xcm1       
        rhe(1,i)=norm(xx(:,2+i))
        thhe(1,i)=dacos(dot_product(ri,xx(:,2+i))/rhe(1,i))
       enddo        
               
       endif
               

       return
       end subroutine cartesian_coords_conv


    end Module pot_calculation

