Module dmc_cal 
  use types
  use global_variables      
  use leg_integrator
  use input_values 
  use math_subroutines
  use phys_cons
  use pot_calculation 
  implicit none      
 Save
 integer,parameter::delta_N_max=10000
 integer,parameter::ndivbr=1000
 real(rk),allocatable:: y1(:,:,:),y2(:,:,:),y3(:,:,:)
 real(rk),allocatable:: zcor(:,:,:),xcor(:,:,:)
 real(rk),allocatable::Er(:),pot(:),eloc(:),erme_wn(:)
 real(rk),allocatable::signo(:),fonda(:,:),alphahe(:),betahe(:)
 real(rk),allocatable::probw(:)
 integer(rk),allocatable::ncell(:),ncount(:)
 real(rk)::ermedia,ermedia2
 real(rk),parameter::cc=20000d0,aa=0.5d0
 real(rk)::alpham,betam,gammam
 real(rk)::alphahe1,betahe1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!Programa que usa el metodo de DMC para calcular estados cuanticos! 
!y la funcion de onda correspondiente para el complejo He_NH3     !
!Variables:
!y1,y2,y3:vectores unitarios que definen las direcciones de los ejes
!principales de la molecula. En el caso de diatómicas, el vector y3 esta dirigido
!a lo largo del eje de la misma
!zcor: Coordenadas de los atomos de helio
!xcor: Coordenadas del centro de masa de la molecula
!!!!!!!!!!!!!!,xx,mass!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

Contains

Subroutine dmc_drv 
  real(rk)    :: media,sd
  integer(ik) :: icorrida,N,npaso,j,icontr
  real(rk)::rr,req,th1d,th2d,taud,energy,emin

   allocate(y1(3,nw+delta_N_max,nmon),y2(3,nw+delta_N_max,nmon),y3(3,nw+delta_N_max,nmon))
   allocate(zcor(3,nw+delta_N_max,nhe),xcor(3,nw+delta_N_max,nmon))
   allocate(Er(0:nstps),pot(nw+delta_N_max),eloc(nw+delta_N_max))
   allocate(signo(nw+delta_N_max),fonda(nw+delta_N_max,0:1))
   allocate(alphahe(nhe),betahe(nhe),erme_wn(nruns))
   allocate(probw(0:nw),ncount(nw),ncell(0:ndivbr))

         !open(unit=238,file='conf.xyz')!OJO
         !open(unit=239,file='comp.dat')!OJO

       idum=-1234
       if (.not.fixran) call seed_cal(idum)
       write(*,*) 'fixran',idum
       call rotcons_calculation         
       if (ptyp.eq.1) then              
       call potparam                    
       endif                            
       if (nmon.gt.1) call potdimparam  
       A=A*frotmol                      
       B=B*frotmol                      
       C=C*frotmol                      
       diffccm=diffccm*frothe           
                                        
       !write(*,*) "rotcons", A*har2cm,B*har2cm,diffccm*har2cm
       write(*,*) "rotcons", A*har2cm,B*har2cm,diffccm*har2cm
     !!!!!!!!!!!!!!!!!!!Constantes Articulo!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !HF
        A = 20.560_rk*cm2har 
        B = 20.560_rk*cm2har 
       !HCl
       !A = 10.403_rk*cm2har
       !B = 10.403_rk*cm2har
       !HBr
       !A = 8.473_rk*cm2har
       !B = 8.473_rk*cm2har
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  HCL DAVID
    !   A=10.406132287017446_rk*cm2har
    !   B=10.406132287017446_rk*cm2har
    !   diffccm=1.6732322591380111_rk*cm2har!OJO poniendo las mismas que David
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !HF DAVID  
    !A = 20.939_rk*cm2har
    !B = 20.939_rk*cm2har
    !diffccm =1.0_rk/(2.0_rk*36469.6_rk)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       call output_writing(1)


 main_loop:  do icorrida=1,nruns
     

   !initialization of values
         icontr=0               
         ermedia=0
         ermedia2=0d0 
         erme_wn(icorrida)=0d0
         er(0)=0d0 
         N=nw  

       
        call initial_conditions(icorrida) 
       !if (icorrida.eq.1) call dist_writing(N)

       !Ciclo sobre todos los pasos temporales 

steps:  do npaso=1,nstps

     !if (icorrida.eq.1) then
     !  if (npaso.ge.nstps-1000) then
     !  call dist_writing(N)
     !  endif
     !endif

       !Se llama a la subrutina que efectúa la difusión
        call diffusion(npaso,N)
       !Se calcula el potencial de cada replica
        pot=0d0
        call potential_calculation(npaso,N,icontr) 

        if (icorrida.eq.1) then    
           if (npaso.ge.nstps-1000) then
            !if (mod(npaso,100).eq.0) then
              call dist_writing(N)
            !endif
           endif 
        endif
!       !Se llama a la subrutina que efectúa los procesos de vida o muerte
        call branch(npaso,N,icorrida) 
        enddo steps


!       !Se escriben los valores de la energía media final en cada corrida
        write (300,*) erme_wn(icorrida)



enddo main_loop

        call write_emedia
        Close(100)
        close(250)
        Close(300)
        close(400)   

      


End Subroutine dmc_drv

      Subroutine initial_conditions(icorrida)
      real (rk)   :: rr,u,phiinicial,rot,rmin,rmax
      real (rk)   :: phi,theta,psi,rrot(3),rhe(3,nw+delta_N_max,nhe)
      integer :: jhe,jmon,j,icorrida,icont
 
    icont=0
    !Primero. Se asignan las posiciones del centro de masas de la molecula 
        xcor=0.0_rk
        if (moltyp.eq.2) then
                !If we have a dimer we start with the center
                !of masses separated by a distance 2*initd
                !initd is fixed in reading.f90
            xcor(1,:,1)=xcor(1,:,1)-initd
            xcor(1,:,2)=xcor(1,:,2)+initd
        endif
        

  !Segundo. Se asignan los vectores de orientacion de la molecula 
  
       do j=1,nw 

      rot=1.0d0 !poner a cero si se quiere iniciar alineado con el lab
        do jmon=1,nmon
         y1(:,j,jmon)=ei
         y2(:,j,jmon)=ej
         y3(:,j,jmon)=ek
         phiinicial = ran2(idum)*2.d0*pi*rot
        !phiinicial=0.0_rk!OJO
       !Se rota la molecula con respecto al ejex molecular un angulo phiinicial 
      call rotacion(y1(1,j,jmon),Y1(2,j,jmon),Y1(3,j,jmon),y3(1,j,jmon),y3(2,j,jmon),y3(3,j,jmon),phiinicial) 
      call rotacion(y1(1,j,jmon),Y1(2,j,jmon),Y1(3,j,jmon),y2(1,j,jmon),y2(2,j,jmon),y2(3,j,jmon),phiinicial) 
        
        phiinicial = ran2(idum)*2.d0*pi*rot
        !phiinicial=0.0_rk!OJO
       !se rota la molecula con respecto al ejey molecular un angulo phiinicial 
      call rotacion(y2(1,j,jmon),Y2(2,j,jmon),Y2(3,j,jmon),y3(1,j,jmon),y3(2,j,jmon),y3(3,j,jmon),phiinicial) 
      call rotacion(y2(1,j,jmon),Y2(2,j,jmon),Y2(3,j,jmon),y1(1,j,jmon),y1(2,j,jmon),y1(3,j,jmon),phiinicial) 
        

      if (moltyp.gt.2) then !Rotate for a non linear molecule
       phiinicial = ran2(idum)*2.d0*pi*rot 
       !se rota la molecula con respecto al ejez molecular un angulo phiinicial 
       call rotacion(y3(1,j,jmon),Y3(2,j,jmon),Y3(3,j,jmon),y1(1,j,jmon),y1(2,j,jmon),y1(3,j,jmon),phiinicial) 
       call rotacion(y3(1,j,jmon),Y3(2,j,jmon),Y3(3,j,jmon),y2(1,j,jmon),y2(2,j,jmon),y2(3,j,jmon),phiinicial)
      endif

      call euler(y1(:,j,jmon),y2(:,j,jmon),y3(:,j,jmon),phi,theta,psi)
         enddo
 !      
        
!       !Tercero. Se asignan los valores de las posiciones iniciales del helio aleatoriamente
!       !en una esfera de radio rr cerca del mínimo de potencial  
          do jhe=1,nhe
            rmin=6.0_rk
            rmax=15.0_rk
            if (simtyp.eq.2) then
            rmin=rfb;rmax=rfb
            endif       
            rr=rmin+ran2(idum)*(rmax-rmin)        
            u = ran2(idum)*2.d0-1.d0 
            phiinicial = ran2(idum)*2.d0*pi 
                             
            zcor(1,j,jhe) = rr*dsqrt(1-u**2)*dcos(phiinicial) 
            zcor(2,j,jhe) = rr*dsqrt(1-u**2)*dsin(phiinicial) 
            zcor(3,j,jhe) = rr*u  
      
       
          !call rotaciones(rhe(:,j,jhe),rrot,phi,theta,psi)
           
            if (icont.lt.10000) then 
             write (400,1000) rrot,rr,phi,theta,psi
             icont=icont+1
            endif
          enddo 
             enddo
            write(400,*) 'New run'

1000 format(7(F10.4))

        return
     end Subroutine initial_conditions
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!Esta subrutina hace dar un paso a todos los caminantes, rota la mo-!
!lécula,  manda a calcular el potencial y recalcula                 ! 
!sus posiciones                                                     ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     Subroutine diffusion(npaso,N) 
       real (rk),dimension(3)::Fq_x 
       real (rk),dimension(3,nhe)::Fq_z 
       real (rk)  ::phix,phiy,phiz,rvec_he(3)
       integer     :: npaso,j,i,jhe,jmon
       integer,intent(in)::N 
!        
       do j=1,N 
!        
        If (is) Then
        Call Fq(j,Fq_z,Fq_x)
        Else
        Fq_z=0d0
        Fq_x=0d0
        endif
        
        if (simtyp.eq.2) then
                do jhe=1,nhe 
       rvec_he=zcor(:,j,jhe)-xcor(:,j,1)        
       Phix=dsqrt(2d0*dtau*diffcHe)*gasdev(idum)!*0d0 
       Phiy=dsqrt(2d0*dtau*diffcHe)*gasdev(idum)!*0d0 
       Phiz=dsqrt(2d0*dtau*diffcHe)*gasdev(idum)!*0d0 

!       !1a-Se rota el helio con respecto al ejex de laboratorio un angulo phix 
       call rotacion(ei(1),ei(2),ei(3),rvec_he(1),rvec_he(2),rvec_he(3),phix) 
!        
!       !1b-Se rota el helio con respecto al ejey de laboratorio un angulo phiy 
       call rotacion(ej(1),ej(2),ej(3),rvec_he(1),rvec_he(2),rvec_he(3),phiy) 
!        
!       !1c-Se rota el helio con respecto al ejez de laboratorio un angulo phiz 
       call rotacion(ek(1),ek(2),ek(3),rvec_he(1),rvec_he(2),rvec_he(3),phiz) 
      !if (abs(rdis(xcor(:,j),zcor(:,j,jhe))-rfb).gt.tol) stop 'distancia'
             enddo

        else        

!       !primero, se mueve el centro de masa de la molecula y el atomo de helio
       !diffccm=1.6732_rk/har2cm
       do i=1,3 
       do jmon=1,nmon
       xcor(i,j,jmon)=xcor(i,j,jmon)+gasdev(idum)*dsqrt(dtau*2d0*diffccm)+Fq_x(i) 
       enddo
       do jhe=1,nhe
       zcor(i,j,jhe)=zcor(i,j,jhe)+gasdev(idum)*dsqrt(dtau*2d0*diffcHe)+Fq_z(i,jhe) 
       enddo
       enddo
       endif

       do jmon=1,nmon
!       !Ahora se calculan los angulos para rotar la molecula
          Phix=dsqrt(2d0*dtau*A)*gasdev(idum)!*0.0_rk 
          Phiy=dsqrt(2d0*dtau*B)*gasdev(idum)!*0.0_rk
          if (moltyp.gt.2) Phiz=dsqrt(2.0_rk*dtau*C)*gasdev(idum)!*0d0 
          
!          !Se llama a la subrutina que efectua la rotacion 
!           
!          !2a-Se rota la molecula con respecto al ejex molecular un angulo phix 
          call rotacion(y1(1,j,jmon),Y1(2,j,jmon),Y1(3,j,jmon),y3(1,j,jmon),y3(2,j,jmon),y3(3,j,jmon),phix) 
          call rotacion(y1(1,j,jmon),Y1(2,j,jmon),Y1(3,j,jmon),y2(1,j,jmon),y2(2,j,jmon),y2(3,j,jmon),phix) 
!           
!          !2b-se rota la molecula con respecto al ejey molecular un angulo phiy 
          call rotacion(y2(1,j,jmon),Y2(2,j,jmon),Y2(3,j,jmon),y3(1,j,jmon),y3(2,j,jmon),y3(3,j,jmon),phiy) 
          call rotacion(y2(1,j,jmon),Y2(2,j,jmon),Y2(3,j,jmon),y1(1,j,jmon),y1(2,j,jmon),y1(3,j,jmon),phiy) 
!           
!          !2c-se rota la molecula con respecto al ejez molecular un angulo phiz 
          if (moltyp.gt.2) then
          call rotacion(y3(1,j,jmon),Y3(2,j,jmon),Y3(3,j,jmon),y1(1,j,jmon),y1(2,j,jmon),y1(3,j,jmon),phiz) 
          call rotacion(y3(1,j,jmon),Y3(2,j,jmon),Y3(3,j,jmon),y2(1,j,jmon),y2(2,j,jmon),y2(3,j,jmon),phiz)
          endif
       enddo
      
   

       enddo


       return
       end Subroutine diffusion

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!Subrutina que calcula el potencial de cada una de las réplicas!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       Subroutine potential_calculation(npaso,N,icontr)
        real (rk) :: theta,psi,rhehe,vvhehe,v_molhe,sumpot      
        real (rk) :: phig,thetag,rg,xvec(3),rad,potdim      
        real (rk),dimension(3)::rrot,rrotb,rvec
        real(rk):: rhe(3,nw+delta_N_max,nhe)
        integer(ik) ::i,j,jhe,jmon,kmon,khe,icontr,npaso
        integer(ik),intent(in)::N 
        real(rk)::xcmass(3),cth,ctha,cthb,phi

        

     do j=1,N

          !Potential Molecule(s)-heliums
        sumpot=0.0_rk


        if (hmon) then
          do jhe=1,nhe
             do jmon=1,nmon
             call potcalc(y1(:,j,jmon),y2(:,j,jmon),y3(:,j,jmon),&
                     &xcor(:,j,jmon),zcor(:,j,jhe),cth,v_molhe)
             sumpot=sumpot+v_molhe
             enddo
          enddo
        endif

        pot(j)=sumpot*cm2har
  

        !if (.false.) then

          !Dimer Potential(in case there are more than one molecule)
          sumpot=0.0_rk
           if (mmon) then
             do jmon=1,nmon-1
                     do kmon=jmon+1,nmon
                             call potdimer(j,y3(:,j,jmon),y3(:,j,kmon)&
                                    &,xcor(:,j,jmon),xcor(:,j,kmon),ctha,cthb,phi,potdim)
                             sumpot=sumpot+potdim
                      enddo
              enddo
          endif

        pot(j)=pot(j)+sumpot*cm2har

       !A este potencial hay que sumarle todas las interacciones he-he

       sumpot=0.0_rk
          if (hhon) then
            do jhe=1,nhe-1
              do khe=jhe+1,nhe

          rhehe=dsqrt(dot_product(zcor(:,j,jhe)-zcor(:,j,khe),zcor(:,j,jhe)-zcor(:,j,khe)))
          call pothehe(vvhehe,rhehe)
          sumpot=sumpot+vvhehe

              enddo
            enddo
          endif

      pot(j)=pot(j)+sumpot

      !endif

      signo(j)=1.0_rk

      if (state.ne.0) then
       !if (nmon.eq.1) then
       !  ctha=1.0_rk
       !  do i=1,nhe 
       !    ctha=ctha*y3(3,j,i)
       !  enddo
       !endif
        ctha=y3(3,j,1)
        call wvfunct(ctha,cthb,phi,j,mod(npaso,2))
        signo(j)=fonda(j,1)*fonda(j,0)
      endif


     
     Enddo 
      
        return 
        End Subroutine potential_calculation 



        subroutine dist_writing(N)
          integer(ik),intent(in)::N      
          integer(ik)::jmon,kmon,j,k
          real(rk)::ri(3),rj(3),rij(3),rad,xx(3,2*nmon+nhe),avgxx(3,2*nmon+nhe)
          real(rk)::ct1,ct2,ph,xxmon(3,nmon),zhe(3,nhe),rhe(nmon,nhe),thhe(nmon,nhe)

      !Calculation of all the configurations
      
      avgxx=0.0_rk

      do j=1,N    
          
                !Warning:For the moment this is only implemented for a maximum of
                !two monomers

                   ri=y3(:,j,1)          
                   if (nmon.eq.2) then
                   rj=y3(:,j,2)          
                   rij=(xcor(:,j,2)-xcor(:,j,1))
                   rad=norm(rij)
                   endif
                   do k=1,nmon
                   xxmon(:,k)=xcor(:,j,k)
                   enddo
                   do k=1,nhe        
                   zhe(:,k)=zcor(:,j,k)
                   enddo

          call cartesian_coords_conv(ri,rj,rij,rad,xx,xxmon,zhe,rhe,thhe)


                
         !avgxx=avgxx+xx
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!debug
           ! write(250,1000) ri
           ! write(250,1000) rj
           ! write(250,1000) rij
           ! write(250,1000) rad
           ! write(250,1000) xxmon
           ! write(250,1000) zhe
           ! write(250,1000) rhe
           ! write(250,1000) thhe
            !call write_xyz(500,2*nmon+nhe,att%nom,xx*bo2ar) 
1000 format(*(F10.4,'&'))        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






      enddo
           !avgxx=avgxx/dfloat(N)
            call write_xyz(500,2*nmon+nhe,att%nom,xx)

          
           






                 return
                end Subroutine dist_writing
         
         
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        !Subrutina que calcula las particulas que mueren, las que  
        !sobreviven y las que nacen 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        Subroutine branch(npaso,N,icorrida) 
        real (rk) :: xnew(3,nw+delta_N_max,nmon),y1new(3,nw+delta_N_max,nmon),y2new(3,nw+delta_N_max,nmon)
        real (rk) ::y3new(3,nw+delta_N_max,nmon),znew(3,nw+delta_N_max,nhe),potnew(nw+delta_N_max)
        integer :: npaso,N,mbranch,nn,ikount,j,jhe,jmon,icorrida,kk 
        real (rk) :: W,alphar,erme_wn2,er_wn,erme,term,ekineticmol,ekinetiche,erme2
        real (rk) :: sumpotential,fondanew(nw+delta_N_max),aleat

       ikount=0 
!       !Este ciclo calcula mn para cada particula y suma el potencial total 
!       !para todos los caminantes vivos 
!         
       do j=1,N 


       If (is) then
       Call Ecinetica(j)
       Else
       Eloc(j)=0.0_rk
       Endif
       Eloc(j)=Eloc(j)+pot(j) 


       W=dexp(-(Eloc(j)-Er(npaso-1))*dtau) 

         
       aleat=ran2(idum)
       If(W+aleat.gt.3d0) Then 
       mbranch=3 
       Else 
       mbranch=Int(W+aleat) 
       Endif 


     if (state.ne.0) then
       If (signo(j).lt.0.0_rk) then
               mbranch=0
       else
              if (npaso.gt.1) then
               if (recross_corr(j)) mbranch=0
              endif
       endif 
      endif 


                                 
            ikount=ikount+mbranch


          !Si se excede el numero maximo de caminantes permitido, se para 
            If (ikount.gt.nw+delta_N_max) Then 
            write(*,*) ikount,npaso         
            stop 'Mas caminantes que los permitidos'  
            endif

           do nn=ikount-mbranch+1,ikount  
             potnew(nn)=pot(j) 

            kk=mod(npaso,2)
            fondanew(nn)=fonda(j,kk) 

           do jmon=1,nmon  
            xnew(:,nn,jmon) = xcor(:,j,jmon)
            y1new(:,nn,jmon) = y1(:,j,jmon)  
            y2new(:,nn,jmon) = y2(:,j,jmon)  
            y3new(:,nn,jmon) = y3(:,j,jmon)  
           enddo  
           

           do jhe=1,nhe  
            znew(:,nn,jhe) = zcor(:,j,jhe)  
           enddo
                                               
            enddo
                        
             Enddo 

!        !Este ciclo asigna los valores nuevos de todos las variables a los nuevos caminantes 
        do j=1,ikount 
!            
            pot(j)=potnew(j)
            fonda(j,kk)=fondanew(j)
            
          do jmon=1,nmon
             xcor(:,j,jmon) = xnew(:,j,jmon)  
          enddo
           

          do jhe=1,nhe  
            zcor(:,j,jhe) = znew(:,j,jhe)  
          enddo  
                                                  
          do jmon=1,nmon
             y1(:,j,jmon) = y1new(:,j,jmon)  
             y2(:,j,jmon) = y2new(:,j,jmon)  
             y3(:,j,jmon) = y3new(:,j,jmon)  
          enddo 

         enddo 
         
         
         
        !Aqui se calcula la nueva energia de referencia para el proximo paso 
          

           sumpotential=sum(pot)/dfloat(ikount)  



          N=ikount 
        
         term = dfloat(N)/dfloat(nw) 
         alphar = 1.0_rk/dtau 
         er(npaso) = sumpotential + alphar*(1.0_rk-term)         
         ermedia2=ermedia2+er(npaso)
         erme2=ermedia2/dfloat(npaso)
         erme_wn2=erme2*har2cm 


         if (npaso.gt.min(20000,nstps/2)) Then 
         ermedia=ermedia+er(npaso) 
         erme=ermedia/(dfloat(npaso)-dfloat(min(20000,nstps/2))) 
         erme_wn(icorrida)=erme*har2cm 
         endif 
         er_wn=er(npaso)*har2cm 
         if (mod(npaso,100).eq.0) Then 
         write(6,55) npaso,er_wn,erme_wn(icorrida),erme_wn2,ikount,icorrida 
         write(100,55) npaso,er_wn,erme_wn(icorrida),erme_wn2,ikount,icorrida 
         endif 

       55 format (I6,2x,3(F14.8,2x),I6,2x,I2)
        
          
        return 
        end Subroutine branch 

        Subroutine branch2(npaso,N,icorrida) 
        real (rk) :: xnew(3,nw+delta_N_max,nmon),y1new(3,nw+delta_N_max,nmon),y2new(3,nw+delta_N_max,nmon)
        real (rk) ::y3new(3,nw+delta_N_max,nmon),znew(3,nw+delta_N_max,nhe),potnew(nw+delta_N_max)
        integer :: npaso,N,mbranch,nn,ikount,i,j,jhe,jmon,icorrida,kk,k
        real (rk) :: sumw,W,alphar,erme_wn2,er_wn,erme,term,ekineticmol,ekinetiche,erme2,delta
        real (rk) :: sumpotential,fondanew(nw+delta_N_max),aleat

        k=mod(npaso,2)
         call copy_arrangments_global(xnew,xcor,znew,zcor,y1new,y1,y2new,y2,y3new,y3,potnew,pot,fondanew,fonda(:,k))

        !sumpotential=sum(pot)/N 
        !write(*,*) 'sum0',sumpotential*har2cm,N!OJO

       do j=1,N 
         If (is) then
         Call Ecinetica(j)
         Else 
         Eloc (j)=0.0_rk
         Endif
         Eloc(j)=Eloc(j)+pot(j) 
         probw(j)=dexp(-(Eloc(j)-Er(npaso-1))*dtau) 
       enddo

        probw=probw/sum(probw)
       !write(*,*) 'suma',sum(probw),N

         probw(0)=0.0_rk
        do j=1,N
         probw(j)=probw(j-1)+probw(j)
        enddo


        ncount=0
        ncount(ndivbr)=1
        ncell(0)=0
        ncell(ndivbr)=N
        kk=1
        delta=1.0_rk/ndivbr

        !Assignation of the correspondent interval for probw(j)
        do j=1,N-1
             if (probw(j).gt.dfloat(kk)/ndivbr) then
             ncell(kk)=j-1 
             kk=kk+1
             endif        
        enddo        
            if (kk.lt.ndivbr) ncell(kk)=N-1


         !  do j=1,ndivbr
         !   write(345,*) j,ncell(j)
         !  enddo

         !  do j=1,N
         !   write(345,*) j,probw(j),dexp(-(Eloc(j)-Er(npaso-1))*dtau),eloc(j)*dtau
         !  enddo

         !  stop

         !Generation of random numbers and assignation to a specific walker

         do j=1,N
            aleat=ran2(idum)     
            nn=int(aleat*ndivbr)+1     
            do i=ncell(nn-1),ncell(nn)
                   !if (i.eq.ncell(nn-1)) write(*,*) 'check',j,aleat,probw(ncell(nn-1)),probw(ncell(nn))
                   if (aleat.lt.probw(ncell(nn-1))) then
                         pot(j)=potnew(i)
                         fonda(j,k)=fondanew(i)
                      do jmon=1,nmon  
                        xcor(:,j,jmon) = xnew(:,i,jmon) 
                        y1(:,j,jmon)   = y1new(:,i,jmon)  
                        y2(:,j,jmon)   = y2new(:,i,jmon)  
                        y3(:,j,jmon)   = y3new(:,i,jmon)  
                      enddo
                      do jhe=1,nhe  
                        zcor(:,j,jhe)  = znew(:,i,jhe)  
                      enddo 
                      !write(345,*) j,i
                      exit
                   endif

                   if (aleat.gt.probw(ncell(nn))) then 
                         pot(j)=potnew(ncell(nn)+1)
                         fonda(j,k)=fondanew(ncell(nn)+1)
                      do jmon=1,nmon  
                        xcor(:,j,jmon) = xnew(:,ncell(nn)+1,jmon) 
                        y1(:,j,jmon)   = y1new(:,ncell(nn)+1,jmon)  
                        y2(:,j,jmon)   = y2new(:,ncell(nn)+1,jmon)  
                        y3(:,j,jmon)   = y3new(:,ncell(nn)+1,jmon)  
                      enddo
                      do jhe=1,nhe  
                        zcor(:,j,jhe)  = znew(:,ncell(nn)+1,jhe)  
                      enddo 
                      !write(345,*) j,i
                      exit
                   endif

                  if(aleat.lt.probw(i)) then
                         !write(*,*) 'Aqui',j,i
                         !if (j.ge.10) stop
                         !read*
                         pot(j)=potnew(i)
                         fonda(j,k)=fondanew(i)
                      do jmon=1,nmon  
                        xcor(:,j,jmon) = xnew(:,i,jmon) 
                        y1(:,j,jmon)   = y1new(:,i,jmon)  
                        y2(:,j,jmon)   = y2new(:,i,jmon)  
                        y3(:,j,jmon)   = y3new(:,i,jmon)  
                      enddo
                      do jhe=1,nhe  
                        zcor(:,j,jhe)  = znew(:,i,jhe)  
                      enddo 
                      !write(345,*) j,i
                      exit
                  endif
            enddo
         enddo



         sumpotential=sum(pot)/N 
        !write(*,*) 'sum',sumpotential*har2cm,N
         er(npaso) = sumpotential 
         ermedia2=ermedia2+er(npaso)
         erme2=ermedia2/dfloat(npaso)
         erme_wn2=erme2*har2cm 


         if (npaso.gt.min(20000,nstps/2)) Then 
         ermedia=ermedia+er(npaso) 
         erme=ermedia/(dfloat(npaso)-dfloat(min(20000,nstps/2))) 
         erme_wn(icorrida)=erme*har2cm 
         endif 
         er_wn=er(npaso)*har2cm 
         if (mod(npaso,100).eq.0) Then 
         if (mod(npaso,5000).eq.0) Then 
         write(6,55) npaso,er_wn,erme_wn(icorrida),erme_wn2,N,icorrida 
         endif
         write(100,55) npaso,er_wn,erme_wn(icorrida),erme_wn2,N,icorrida 
         endif 

       55 format (I6,2x,3(F14.8,2x),I6,2x,I2)
        
          
        return 
        end Subroutine branch2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   Subroutine write_emedia
   real (rk):: dato,media,sd
   integer:: j
   
       media=0.0_rk
       sd=0.0_rk
       
       do j=1,nruns
       media=media+erme_wn(j)
       enddo
       media=media/nruns
       do j=1,nruns
       sd=sd+(erme_wn(j)-media)**2
       enddo       
       sd=sqrt(sd/(nruns))
       
       write (600,1000) nruns,media,sd

1000  Format('After ',I3,'runs, the mean energy is',F10.4,'+/-',F10.4,' cm-1')

       return 

       End subroutine write_emedia
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     
        subroutine wvfunct(ctha,cthb,phi,j,ind)
        real (rk)::phi,ctha,cthb
        integer::ind,j
         
         if (state.eq.1) then
         !fonda(j,ind)=dcos(dacos(ctha)+dacos(cthb))
         fonda (j,ind) = ctha
         endif


         return
         
     
        end subroutine wvfunct
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        subroutine rotacion(xn,yn,zn,xh,yh,zh,phiang) 
        real (rk)::r,c,s,t,x,y,z,xn,yn,zn,r1,r2,r3,xh,yh,zh,phiang
        
        !  rotation from  
        !c     ` 
         
              r=dsqrt(xn**2+yn**2+zn**2) 
              
              c=dcos(phiang) 
              s=dsin(phiang) 
              t=1.d0-c 
               
        !c     unit vector along axis of rotation 
         
              x=xn/r 
              y=yn/r 
              z=zn/r 
               
        !c     point to be rotated 
         
              r1=xh 
              r2=yh 
              r3=zh 
         
              xh=c*r1+s*(z*r2-y*r3)+x*t*(x*r1+y*r2+z*r3) 
              yh=c*r2+s*(x*r3-z*r1)+y*t*(x*r1+y*r2+z*r3) 
              zh=c*r3+s*(y*r1-x*r2)+z*t*(x*r1+y*r2+z*r3) 
         
              return 
              end Subroutine rotacion 
         
         
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Esta subrutina rota un vector en el sistema rotado a un vector en el  labframe 
        !usando los ángulos de euler.   
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        Subroutine  rotaciones(y1,y2,phi,theta,psi)
          Real (rk),intent(in):: y1(3),phi,theta,psi
          Real (rk),intent(inout)::y2(3)
          Real (rk)::matb(3,3)
          Real (rk)::prueba1,prueba2
          
              !Comprobación. Se realiza la transformación para alinear los vectores
                  
                
                      !Matriz directa, al multiplicar las coordenadas de un vector en el sistema
                      !fijo, nos da las coordenadas del mismo en el sistema rotado.

                     matb=reshape((/cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi),&
                                   -sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi),&
                                    sin(theta)*sin(phi),&
                                    cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi),&
                                   -sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi),&
                                   -sin(theta)*cos(phi),&
                                    sin(theta)*sin(psi),&
                                    sin(theta)*cos(psi),&
                                    cos(theta)/),(/3,3/))
                    
               y2=matmul(matb,y1)
         
             
            return
            end subroutine rotaciones


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Subrutina que calcula el potencial he-he
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              subroutine pothehe(vv1,rhe) !he-he potential
              Real (rk) r,vv1,rm,rhe,x,vst1,vst2,fx,vst,econ,alpha,beta,aa,d,eps,c6,c8,c10
             
              aa=1.8443101d05
              alpha=10.43329537
              beta=-2.27965105
              rm=2.963
              d=1.4826
              eps=10.948
              econ=3.157865d05
              c6=1.36745214
              c8=0.42123807
              c10=0.17473318
              r=0.529177*rhe
              x=r/rm
              vst1=aa*dexp(-alpha*x+beta*x**2)
              vst2=c6/x**6+c8/x**8+c10/x**10
              if(x.lt.d) then
                 fx=dexp(-(d/x-1.d0)**2)
              else
                 fx=1.d0
              endif
              vst=vst1-fx*vst2
              vv1=eps*vst/econ
              return
              end subroutine pothehe 
              

           

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        !Subrutina que calcula la fuerza cuantica para molecula y atomos de helio!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE Fq (j,Fq_z,Fq_x)
               INTEGER :: j,k, kk, kkk
               REAL (rk) ::  R
               REAL (rk),  DIMENSION (3,nhe) :: FDER
               REAL (rk),  DIMENSION (3,nhe) :: Fq_z
               REAL (rk),  DIMENSION (3) :: Fq_x
               REAL (rk),  DIMENSION (3) :: Rcords
              

               DO k=1,nhe
                        Rcords=xcor(:,j,1)-zcor(:,j,k)
                        R=DSQRT(DOT_PRODUCT(Rcords,Rcords))

                        FDER(:,k)=(xcor(:,j,1)-zcor(:,j,k))*(5.d0*cc/R**7-aa/R)
               ENDDO

               DO k=1,nhe
                  Fq_z(:,k)=-dtau/mHe*FDER(:,k)
               ENDDO

              
                  DO kk=1,3
                  Fq_x(kk)=dtau/mtot*SUM(FDER(kk,:))
                  ENDDO   
               
               return
                END Subroutine  Fq
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !subrutina que calcula la parte de la energia cinetica para elocal. Posteriormente se suma el termino del potencial en      
        !el programa que llama a esta subrutina.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        Subroutine Ecinetica(j)
        real(rk) ::  R,kinetic,Rkkk,suma,ekineticmol,ekinetiche
        real(rk) ::  dd1,dd2
        real(rk),DIMENSION (3) :: Rcords,rkk
        real(rk),Dimension(nhe)::d2d
        integer::k,kk,kkk,j

        ekineticmol=0d0
        ekinetiche=0d0

       do k=1,nhe
           Rcords=xcor(:,j,1)-zcor(:,j,k) 
           R=DSQRT(DOT_PRODUCT(Rcords,Rcords))  
            
       do kk=1,3
        suma=0d0
        

        do kkk=1,nhe
            rkk=xcor(:,j,1)-zcor(:,j,kkk) 
            Rkkk=DSQRT(DOT_PRODUCT(rkk,rkk)) 
           suma=suma+dfx(rkkk,kk,kkk,j)
        enddo  

        dd2=d2fx(r,kk,k,j)
        dd1=dfx(r,kk,k,j)
        ekineticmol=ekineticmol+dd2+dd1*suma
        ekinetiche=ekinetiche+dd2+dd1**2 


       enddo
       enddo
          

       Eloc(j)=-1d0/(2d0*mHe)*ekinetiche-1d0/(2d0*mtot)*ekineticmol


        return 
        end subroutine Ecinetica
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Real (rk) Function df(rr)
        Real (rk)::rr

        df=5d0*cc/rr**6-aa

        endfunction df
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Real (rk) Function d2f(rr)
        Real (rk)::rr

        d2f=-30d0*cc/rr**7

        endfunction d2f


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Real (rk) Function dfx(rr,kder,jhe,j)
        Real (rk)::rr
        integer::kder,jhe,j

        dfx=df(rr)*(xcor(kder,j,1)-zcor(kder,j,jhe))/rr

        endfunction dfx

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Real (rk) Function d2fx(rr,kder,jhe,j)
           Real (rk)::rr,dif
          integer::kder,jhe,j

        dif=xcor(kder,j,1)-zcor(kder,j,jhe)
        d2fx=d2f(rr)*(dif**2)/rr**2+df(rr)*(1-dif**2/rr**2)/rr

        endfunction d2fx


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine euler(eir,ejr,ekr,phi,theta,psi)
         Real (rk):: matb(3,3)
         Real (rk):: eirr(3),ejrr(3),ekrr(3),cv(3)
         Real (rk), Intent(in)::eir(3),ejr(3),ekr(3)
         Real (rk),parameter::ei(3)=(/1d0,0d0,0d0/),ej(3)=(/0d0,1d0,0d0/),ek(3)=(/0d0,0d0,1d0/) 
         Real (rk)::phi,theta,psi,pru1,pru2,pru3,pru4,pru5,pru6,pi
         Real (rk)::phip,thetap,psip,valprueba1,valprueba2  
         Real (rk)::ter11,ter12,ter13,ter21,ter22,ter23,ter31,ter32,ter33,prudef 
         integer::ind1,ind2 


            pi=dacos(-1d0)

           
                !Vectores rotados           
                eirr=eir
                ejrr=ejr
                ekrr=ekr

                 !Calculo de los angulos de Euler

               
         
                 theta=dacos(ekrr(3))
                 if (theta.gt.1d-14) then
                 phi=datan2(ekrr(1),-ekrr(2))
                 if (phi.lt.0d0) phi=phi+2d0*pi
                 !print*,  'componentes Euler',eirr(3),e
                 psi=datan2(eirr(3),ejrr(3)) 
                 if (psi.lt.0d0) psi=psi+2d0*pi
                 else
                 phi=dacos(eirr(1))
                 psi=0d0
                 endif
              
            
                !Comprobación. Se realiza la transformación para alinear los vectores
                  
                
                      !Matriz directa, al multiplicar las coordenadas de un vector en el sistema
                      !fijo, nos da las coordenadas del mismo en el sistema rotado.

                     matb=reshape((/cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi),&
                                   -sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi),&
                                    sin(theta)*sin(phi),&
                                    cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi),&
                                   -sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi),&
                                   -sin(theta)*cos(phi),&
                                    sin(theta)*sin(psi),&
                                    sin(theta)*cos(psi),&
                                    cos(theta)/),(/3,3/))

                    
                   !Inversa de la matriz anterior, al multiplicar las coordenadas de un vector en el sistema
                   !rotado, nos da las coordenadas del mismo en el sistema fijo.



                  eirr=Matmul(Matb,eirr)
                  ejrr=Matmul(matb,ejrr)
                  ekrr=Matmul(matb,ekrr)

             
              
                 pru1=dot_product(ei-eirr,ei-eirr)
                 pru2=dot_product(ej-ejrr,ej-ejrr)
                 pru3=dot_product(ek-ekrr,ek-ekrr)
                 
                 prudef=pru1+pru2+pru3

                if (prudef.gt.1d-15) then
                print*,'mal algo en angulos de Euler'
                print*,'pruebas',pru1,pru2,pru3,prudef
                print*,'angulos',phi,theta,psi
                stop
                endif
                     

                  end Subroutine euler
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine copy_arrangments_global(x1,x2,z1,z2,y11,y12,y21,y22,y31,y32,pot1,pot2,fonda1,fonda2)

        real (rk) :: x1(3,nw+delta_N_max,nmon),y11(3,nw+delta_N_max,nmon),y21(3,nw+delta_N_max,nmon)
        real (rk) ::y31(3,nw+delta_N_max,nmon),z1(3,nw+delta_N_max,nhe),pot1(nw+delta_N_max)
        real (rk) :: fonda1(nw+delta_N_max)
        real (rk) :: x2(3,nw+delta_N_max,nmon),y12(3,nw+delta_N_max,nmon),y22(3,nw+delta_N_max,nmon)
        real (rk) ::y32(3,nw+delta_N_max,nmon),z2(3,nw+delta_N_max,nhe),pot2(nw+delta_N_max)
        real (rk) :: fonda2(nw+delta_N_max)

          x1=x2
          z1=z2
          y11=y12
          y21=y22
          y31=y32
          pot1=pot2
          fonda1=fonda2

           return
         end Subroutine copy_arrangments_global

       subroutine copy_arrangments_scalar(x1,x2,z1,z2,y11,y12,y21,y22,y31,y32,pot1,pot2,fonda1,fonda2)
        real (rk) :: x1,y11,y21,y31,z1,pot1,fonda1
        real (rk) :: x2,y12,y22,y32,z2,pot2,fonda2


          x1=x2
          z1=z2
          y11=y12
          y21=y22
          y31=y32
          pot1=pot2
          fonda1=fonda2

           return
         end Subroutine copy_arrangments_scalar


           Subroutine output_writing(ind)
           integer(ik),intent(in)::ind
           integer(ik)::i

           Select case (ind)
           case(1)
           write(600,1000) trim(molname)
           write(600,*) 
           write(600,*) 'Physical Properties:'
           write(600,*) 
           write(600,2000) diffccm*har2cm
           do i=1,nhe
           write(600,3000) diffcHe*har2cm
           enddo
           write(600,4000) A*har2cm!,B*har2cm!,C*har2cm
           case default
            write(*,*) 'molname ',trim(molname),' not recognized'
            stop
           end Select

1000   format(' Molecule name:',2x,A)
2000   format(' Molecule diffusion constant(cm-1):',2x,F18.10)
3000   format(' He diffusion constant(cm-1):',2x,F18.10)
4000   format(' Rotational constants(cm-1):',2x,F18.10)
5000   format(' Rotational. Constant(cm-1):',2x,F18.10)


           end subroutine output_writing        


           function recross_corr(j)
            integer(ik),intent(in)::j       
            logical::recross_corr
            real(rk)::prec,thnode

            recross_corr=.false.
            thnode=pi/2.0_rk
            prec=dexp(-2.0_rk*(0.5_rk/A)*abs(dacos(fonda(j,0))-thnode)*abs((dacos(fonda(j,1))-thnode))/dtau)


            if (prec.gt.ran2(idum)) recross_corr=.true.
           !recross_corr=.false.


             return
            end function recross_corr
                     




            end Module dmc_cal
