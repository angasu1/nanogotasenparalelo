c     
c     This program is a DMC code calculating the ground and 1st ex. state for a system of nHe's
c     and one molecule
c     
c     each type of a node has a recrossing correction
c
c     rotations are made in a molecular frame
c     
c     OCT  2006

      implicit real*8(a-h,o-z)

      real*8 xHe(10000,30),yHe(10000,30),zHe(10000,30)
      real*8 xnHe(10000,30),ynHe(10000,30),znHe(10000,30)

      real*8 xHel(10000),yHel(10000),zHel(10000)

      real*8 xcm(10000),ycm(10000),zcm(10000)
      real*8 xncm(10000),yncm(10000),zncm(10000)

      real*8 xA(10000),yA(10000),zA(10000),xray(10000)
      real*8 xnA(10000),ynA(10000),znA(10000)

      real*8 x_fdHe(10000,99), y_fdHe(10000,99), z_fdHe(10000,99) 
      real*8 x_fdMol(10000), y_fdMol(10000), z_fdMol(10000)

      real*8 x_sdMol(10000), y_sdMol(10000), z_sdMol(10000) 
      real*8 x_sdHe(10000,99), y_sdHe(10000,99), z_sdHe(10000,99) 

      real*8 x_axisX(10000),x_axisY(10000),x_axisZ(10000)
      real*8 y_axisX(10000),y_axisY(10000),y_axisZ(10000)

      real*8 x_axisXnew(10000),x_axisYnew(10000),x_axisZnew(10000)
      real*8 y_axisXnew(10000),y_axisYnew(10000),y_axisZnew(10000)

      real*8 v(500000),e(500000),e_wn(500000),avE(500000),aveavE(500000)
      real*8 elnew(10000),el(10000),afinal(99),NofHe(99)

      common/angles/ th(99),ph(99)
      common/cos_thetaHe/ cthHe(99)
      common/param/ c2,c2N,h,nHe
      common/nt/psi,rHe,ntype

********************************** SIMULATION INPUT ****************************

      wvfMol(r) = dexp(-cMol/r**5 - aMol*r)
      
      
      
      SDMol(crd1,crd2,r) = wvfMol(r) * (- 35.d0*cMol*(crd1-crd2)**2/r**9 
     /     + 5.d0*cMol/r**7 + aMol*(crd1-crd2)**2/r**3 - aMol/r )
     /     + wvfMol(r)*((crd1-crd2) * (5.d0*cMol/r**7 - aMol/r))**2
      
      FDMol(crd1,crd2,r) = dexp(-cMol/r**5-aMol*r)*(5.d0*cMol*
     /     (crd1-crd2)/(r**7)-aMol*(crd1 - crd2)/r) 

	
      idum=-693584
      write(*,*)idum
      neq = 1

      if (neq.eq.1) then
         INW =2000
      endif  
      write(*,*)'INW',INW      
      
      NOR = 1
      
      nimps = 0
      
      TSV = 10.d0
      write(*,*)'TSV',TSV   
      
      ntype = 0  !      ntype is d type of d node: 
      write(*,*)'ntype',ntype	
          
c     0 - groundstate
c     1 - cosbeta
c     2 - costeta
c     3 - general
c     4 - adiabatic

*******************************************CONSTANTS****************************
      
      pi = 2.d0*dasin(1.d0)
      cc = 2.99792458d10
      hp = 6.62608d-34
      amu = 1.66054d-27
      har2cm = 2.d0*13.605698d0/1.23985d-4
      aconvert = 0.5291772d0
      convmass = 1.6606d-27/9.1095d-31

      print*,har2cm
      pause

      dmassA =100d0*1.00794d0 * convmass
      dmassB =100d0*35.453d0*convmass
      dmassC = 0.d0 
      
      
      dmassMol =(dmassA + dmassB + dmassC)
      dmassHe = 100d0*4.002602d0 * convmass


      write(*,*)'redmass',dmassHe/convmass
      write(*,*)calc
      B_wn=10.664234d0*0.d0
      pause

      rAB=1.27d0

      rBond = rAB+rBC
	
      if (dmassC.eq.0.d0.and.B_wn.eq.0.d0)then
         
         redmassMol = dmassA * dmassB / (dmassMol)
         redmassMol_kg = (redmassMol / convmass) * amu 
	   rAB_m = rAB*1.d-10
         B_wn = hp / (8.d0*pi**2 * redmassMol_kg * rAB_m**2*cc)

      endif

      if (dmassC.gt.0.d0.and.B_wn.eq.0.d0)then

         dmassA_kg = amu*dmassA/convmass
         dmassC_kg = amu*dmassC/convmass
         rAB_m = rAB*1.d-10
         rBC_m = rBC*1.d-10
         dI = (dmassA_kg*rAB_m**2 + dmassC_kg*rBC_m**2)
         B_wn = hp/(8.d0 * pi**2 * dI * cc)

      endif

      B=B_wn/har2cm

      write(*,*)'B',B_wn   

      diffcHe = 1.d0/(2.d0*dmassHe)
      diffccm = 1.d0/(2.d0*dmassMol)

      er = 0.d0
      sumV = 0.d0
	
      if(nimps.eq.0) then
         write(*,*)'no_importance sampling'
      else
	 write(*,*)'importance sampling'
      endif

      do nn = 1,1
         
      nhe = nn
       

      NOTS=150000
      !NOTS = 150000+nHe*5000
      
      write(*,*)'#He',nHe

      aMol = .8d0
      c0 = 7.d3 * rBond
      cMol = c0 + c0 * (dfloat(nHe)/15.d0)

      print*,'cmol',cMol
      pause

      do l = 1,NOR !primer do
         
**************************INITIALAZING  POSITIONS OF  WALKERS *****************

         if (neq.eq.1) then!primer if
            
             if (nimps.eq.1)then!segundo if
               
              kount = 0
 
 5            kount = kount+1 
                          
               teta = 0.d0
               phi = 0.d0
               
               do j= 1,NHe!segundo do
                  
 10               continue              
                  
                  r = gasdev(idum)*2.d0+7.d0 
                  
                  if(r.gt.0.d0) then!tercer if
                     rr=r
                  else 
                     goto 10
                  endif!cierra tercer if
                  
                  phi  = 2.d0*pi*ran1(idum)
                  teta  = pi*ran1(idum)
                  
                  xHell = rr*dsin(teta)*dcos(phi)
                  yHell = rr**dsin(teta)*dsin(phi)
                  zHell = rr*dcos(teta)              
                  
                  rmax = (5.d0*cMol/aMol)**(1.d0/6.d0)
                  wvfMolmax = wvfMol(rmax)
                             
                  if(wvfMol(rr)/wvfMolmax.lt.ran2(idum))goto 10
                    
                  xHe(kount,j)=xHell
                  yHe(kount,j)=yHell
                  zHe(kount,j)=zHell
                  
               enddo!cierra segundo do
               
               xcm(kount) = 0.d0
               ycm(kount) = 0.d0
               zcm(kount) = 0.d0
                   
               xA(kount) = 0.d0
               yA(kount) = 0.d0
               zA(kount) = 1.d0   

              
               if(kount.lt.INW) go to 5
               
               nwlkrs = kount
               
            else!negativo segundo if (POR ACA!!)
               
               kount = 0
               i = 0
 15            i= i+1
               rr = 7.0d0
               dr = 0.5d0
               
               if(i.le.inw) then!cuarto if
                  
                  kount = kount+1 
                  
                  do j= 1,NHe!tercer do

                     u = ran2(idum)*2.d0-1.d0
                     teta = ran2(idum)*2.d0*pi
                     
                     xHe(kount,j) = rr*dsqrt(1-u**2)*dcos(teta)
                     yHe(kount,j) = rr*dsqrt(1-u**2)*dsin(teta)
                     zHe(kount,j) = rr*u  

		     rr = rr+dr
                  enddo  !cierra tercer do            

                  xcm(kount) = 0.d0
                  ycm(kount) = 0.d0
                  zcm(kount) = 0.d0               
                  
                  xA(kount) = 0.d0
                  yA(kount) = 0.d0	
                  zA(kount) = (-1.d0)**(i)

               endif!cierra cuarto if
               
               if(kount.lt.INW) go to 15
               
               nwlkrs = kount
               
            endif!cierra segundoif
         endif!cierra primer if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
   !no hacer caso a esta parte
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
       
         if(neq.eq.0) then

            kp = 0
 20       kp = kp+1  
            

               read(100+nhe,*,end=25)
     *              xcm(kp),ycm(kp),zcm(kp),xA(kp),yA(kp),zA(kp),
     *                x_axisX(kp),x_axisY(kp),x_axisZ(kp),
     *                y_axisX(kp),y_axisY(kp),y_axisZ(kp), 
     *              (xhe(kp,nhe1),yhe(kp,nhe1),zhe(kp,nhe1),
     *              nhe1=1,nhe)
            
            go to 20

 25         continue

            INW = kp-1
            nwlkrs = INW
            
         endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         nmax = nwlkrs
         n0 = nmax
         er = 0.d0
        
       !Comienza el ciclo en tau. Se hace caminar a los gueyes y se calcula su energia
          
          do k = 1,NOTS  !Inicia segundo do

            ikount = 0
            do i = 1,nmax!Inicia tercer do sobre todos los caminantes Para una unidad de tiempo.

               if (k.eq.1) then !Inicia primer if               
                  el_old = 0.d0                 
               else                 
                  el_old = el(i)             
               endif!termina primer if
               
               if (neq.eq.1)then!Inicia segundo if
                  
                  if (k.eq.1.) then!Inicia tercer if             
                    !Valores de los ejes para t=0 
                     x_axisXX = 1.d0
                     x_axisYY = 0.d0
                     x_axisZZ = 0.d0
                     
                     y_axisXX = 0.d0
                     y_axisYY = 1.d0
                     y_axisZZ = 0.d0
                     
                  else!negativo tercer if
                  !Esto es lo que se tiene para un tiempo distinto de 0   
                     x_axisXX = x_axisX(i)
                     x_axisYY = x_axisY(i)
                     x_axisZZ = x_axisZ(i)
                     
                     y_axisXX = y_axisX(i)
                     y_axisYY = y_axisY(i)
                     y_axisZZ = y_axisZ(i)
                     
                  endif!Termina tercer if
               endif!Termina segundo if
          
*************************************TRANSLATING He and ROTATING MOLECULE***********
               
               xA0 = xA(i) - xcm(i)
               yA0 = yA(i) - ycm(i)
               zA0 = zA(i) - zcm(i)
             
               if(ntype.eq.0) goto 30 ! for ground_state
******no hacer ni puto caso     ********
               rMol = dsqrt(xA0**2 + yA0**2 + zA0**2)
               th(1) = dacos(zA0/rMol)

               if(ntype.eq.1) goto 35 ! for cos_beta

               ph(1) = datan2(yA0,xA0)  

               do j = 1,NHe       
                  
                  xHe0 = xHe(i,j) - xcm(i)
                  yHe0 = yHe(i,j) - ycm(i)
                  zHe0 = zHe(i,j) - zcm(i)
                  
                  rHe = dsqrt(xHe0**2 + yHe0**2 + zHe0**2)

                  th(j+1) = dacos(zHe0/rHe)
                  ph(j+1) = datan2(yHe0,xHe0) 

               enddo 

 35            continue
               
               call node_type
      
               psiold=psi
********************************************************
 30            continue                   

               if(nimps.eq.0) then

                  do ii = 1,nHe

                     gaussxHe = dsqrt(2.d0*diffcHe*TSV) * gasdev(idum)
                     gaussyHe = dsqrt(2.d0*diffcHe*TSV) * gasdev(idum)
                     gausszHe = dsqrt(2.d0*diffcHe*TSV) * gasdev(idum)
                     
                     xHe(i,ii) = xHe(i,ii) + gaussxHe     
                     yHe(i,ii) = yHe(i,ii) + gaussyHe      
                     zHe(i,ii) = zHe(i,ii) + gausszHe  
  
                  enddo

               else

!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
           !!!!!!IMPORTANT SAMPLING!!!!!!!!!!!!       
                  do 40 ii = 1,nHe    
                      
                     
                     rMol = dsqrt((xcm(i) - xHe(i,ii))**2 + 
     *                            (ycm(i) - yHe(i,ii))**2 + 
     *                            (zcm(i) - zHe(i,ii))**2)

                     wfMol=wvfMol(rMol)
                     
                   x_fdMol(ii) = FDMol(xcm(i),xHe(i,ii),rMol)
     *                          / wfMol
                   y_fdMol(ii) = FDMol(ycm(i),yHe(i,ii),rMol)
     *                          / wfMol
                   z_fdMol(ii) = FDMol(zcm(i),zHe(i,ii),rMol)
     *                          / wfMol
          

 40             continue

                  FqxM = 0.d0
                  FqyM = 0.d0
                  FqzM = 0.d0
                          
                  do ii = 1,nHe
                                      
                     FqxM = FqxM + 2.d0*x_fdMol(ii)
                     FqyM = FqyM + 2.d0*y_fdMol(ii)
                     FqzM = FqzM + 2.d0*z_fdMol(ii)
                  
                     FqxHe = -2.d0*x_fdMol(ii)
                     FqyHe = -2.d0*y_fdMol(ii)
                     FqzHe = -2.d0*z_fdMol(ii)

                     driftxHe = TSV * diffcHe * FqxHe
                     driftyHe = TSV * diffcHe * FqyHe
                     driftzHe = TSV * diffcHe * FqzHe                       
                     
                     
                     gaussxHe = dsqrt(2.d0*diffcHe*TSV) * gasdev(idum)
                     gaussyHe = dsqrt(2.d0*diffcHe*TSV) * gasdev(idum)
                     gausszHe = dsqrt(2.d0*diffcHe*TSV) * gasdev(idum)
                     
                     xHe(i,ii) = xHe(i,ii) + gaussxHe + driftxHe       
                     yHe(i,ii) = yHe(i,ii) + gaussyHe + driftyHe       
                     zHe(i,ii) = zHe(i,ii) + gausszHe + driftzHe       
                     
                  enddo
               endif

               gaussxMol = dsqrt(2.d0*diffccm*TSV) * gasdev(idum)
               gaussyMol = dsqrt(2.d0*diffccm*TSV) * gasdev(idum)
               gausszMol = dsqrt(2.d0*diffccm*TSV) * gasdev(idum)


               if(nimps.eq.1) then
                  
                  driftxMol = TSV * diffccm * FqxM
                  driftyMol = TSV * diffccm * FqyM
                  driftzMol = TSV * diffccm * FqzM
                  
                  xtrans = gaussxMol + driftxMol
                  ytrans = gaussyMol + driftyMol
                  ztrans = gausszMol + driftzMol   
                  
               else
                  
                  xtrans = gaussxMol 
                  ytrans = gaussyMol 
                  ztrans = gausszMol  
                  
               endif
               
               xcm(i) = xcm(i) + xtrans           
               ycm(i) = ycm(i) + ytrans 
               zcm(i) = zcm(i) + ztrans 
               
               xA(i) = xA(i) + xtrans
               yA(i) = yA(i) + ytrans              
               zA(i) = zA(i) + ztrans 
               
               xA0 = xA(i) - xcm(i)
               yA0 = yA(i) - ycm(i)
               zA0 = zA(i) - zcm(i)
               
c               rotatational angles

               phiX = (dsqrt(2.d0 * TSV * B)) * gasdev(idum) 
               phiY = (dsqrt(2.d0 * TSV * B)) * gasdev(idum)
               phiZ = (dsqrt(2.d0 * TSV * B)) * gasdev(idum)             
              
c                  call rotate(xA0,yA0,zA0,   !lab frame rotation
c     *                 xnewA,ynewA,znewA,phiX,phiY,phiZ)                  

               h=1.d-2
               
c     function value before rotation
c!!!!!!!!!!!!!!!!Ni puto caso
               if(ntype.eq.0) goto 45 ! for ground_state
               
               rMol = dsqrt(xA0**2 + yA0**2 + zA0**2)        
               th(1) = dacos(zA0/rMol)

               if(ntype.eq.1) goto 50 ! for cos_beta

               ph(1) = datan2(yA0,xA0)  
               
               do j = 1,NHe       
                  
                  xHe0 = xHe(i,j) - xcm(i)
                  yHe0 = yHe(i,j) - ycm(i)
                  zHe0 = zHe(i,j) - zcm(i)
                  
                  rHe = dsqrt(xHe0**2 + yHe0**2 + zHe0**2)

                  th(j+1) = dacos(zHe0/rHe)
                  ph(j+1) = datan2(yHe0,xHe0) 

               enddo 

 50            continue

               call node_type
	
        
               psibefore = psi          
               
c     rotation around X by h
               
               xArecX = xA0
               yArecX = yA0
               zArecX = zA0
               
               call rot(x_axisXX,x_axisYY,x_axisZZ,
     *              xArecX,yArecX,zArecX,h)
               
               
               rMol = dsqrt(xArecX**2 + yArecX**2 + zArecX**2)
               th(1) = dacos(zArecX/rMol)
               ph(1) = datan2(yArecX,xArecX)
               
               call node_type
               
               psiX = psi 
               
               dphiX1 = (psiX-psibefore)/h
               
               
c     rotation around Y by h
               
               xArecY = xA0
               yArecY = yA0
               zArecY = zA0
               
               
               call rot(y_axisXX,y_axisYY,y_axisZZ,
     *              xArecY,yArecY,zArecY,h) 

               rMol = dsqrt(xArecX**2 + yArecX**2 + zArecX**2)
               th(1) = dacos(zArecX/rMol)
               ph(1) = datan2(yArecX,xArecX)

               call node_type
               
               psiY = psi 
               
               dphiY1 = (psiY-psibefore)/h

c     rotation around X by phi
c!!!!!!!!!!!!!!!!!!!!A partir de aqui sigue
 45            continue

               call rot(x_axisXX,x_axisYY,x_axisZZ,xA0,yA0,zA0,
     *              phiX) 
               
               
               
               call rot(x_axisXX,x_axisYY,x_axisZZ,y_axisXX,y_axisYY,
     *              y_axisZZ,phiX) 
               
c     rotation around Y by phi
               
               call rot(y_axisXX,y_axisYY,y_axisZZ,xA0,yA0,zA0,
     *              phiY)        
               
               call rot(y_axisXX,y_axisYY,y_axisZZ,x_axisXX,x_axisYY,
     *              x_axisZZ,phiY) 
               
               
c     function value after rotation
!!!!!!!!!!!!!!!!!!!Ni p c 
              if(ntype.eq.0) goto 55 ! for ground_state
               
               rMol = dsqrt(xA0**2 + yA0**2 + zA0**2)
               th(1) = dacos(zA0/rMol)
               ph(1) = datan2(yA0,xA0)
               
               call node_type
               
               psiafter = psi
                         
c     rotation around X by h
               
               xArecX = xA0
               yArecX = yA0
               zArecX = zA0
               
               call rot(x_axisXX,x_axisYY,x_axisZZ,
     *              xArecX,yArecX,zArecX,h)
               
               
               rMol = dsqrt(xArecX**2 + yArecX**2 + zArecX**2)
               th(1) = dacos(zArecX/rMol)
               ph(1) = datan2(yArecX,xArecX)

               call node_type
               
               psiX = psi 
               
               dphiX2 = (psiX-psiafter)/h      
               
c     rotation around Y by h
               
               xArecY = xA0
               yArecY = yA0
               zArecY = zA0
               
               
               call rot(y_axisXX,y_axisYY,y_axisZZ,
     *              xArecY,yArecY,zArecY,h) 
               
               rMol = dsqrt(xArecX**2 + yArecX**2 + zArecX**2)
               th(1) = dacos(zArecX/rMol)
               ph(1) = datan2(yArecX,xArecX)

               call node_type
               
               psiY = psi 
               
               dphiY2 = (psiY-psiafter)/h
               
               d1 = dabs(psibefore)/dsqrt(dphiX1**2+dphiY1**2)
               d2 = dabs(psiafter)/dsqrt(dphiX2**2+dphiY2**2)
               
               Precross = dexp(-d1*d2/(TSV*B))
              
               precross=-1.d0
               if (Precross.gt.ran1(idum))goto 60
               if (psiold*psiafter.lt.0.d0)goto 60
!!!!!!!!!!!!!!!continua aqui
 55            continue   
           
               x_axisX(i) = x_axisXX  
               x_axisY(i) = x_axisYY  
               x_axisZ(i) = x_axisZZ 
               
               
               y_axisX(i) = y_axisXX  
               y_axisY(i) = y_axisYY  
               y_axisZ(i) = y_axisZZ 


                  sumVheMol = 0.d0
       
                  do j = 1,nhe
                     
                     xHe0= xHe(i,j) - xcm(i)
                     yHe0= yHe(i,j) - ycm(i)
                     zHe0= zHe(i,j) - zcm(i)
                     
                     ab = xHe0*xA0 + yHe0*yA0 + zHe0*zA0
                     a2b2 = (xHe0**2 + yHe0**2 + zHe0**2)
     *                    * (xA0**2 + yA0**2 + zA0**2)
                     
                     rheMol = dsqrt(xHe0**2 + yHe0**2 + zHe0**2)

                     costh = ab/dsqrt(a2b2)

                     call pothemol(vheMol,rheMol,costh)                  

                     sumVheMol = sumVheMol+vheMol

                  enddo

                  if (NHe.gt.1)then 
     
                     sumVhehe = 0.d0
                     do m = 1,nhe
                     
                        do j = m+1,nhe

                           rhehe = dsqrt((xHe(i,m) - xHe(i,j))**2 +
     *                                   (yHe(i,m) - yHe(i,j))**2 +
     *                                   (zHe(i,m) - zHe(i,j))**2) 

                           call pothehe(vvhehe,rhehe)

                           sumVhehe = sumVhehe + vvhehe  
                        enddo
                        
                     enddo

                  endif 

                  vv = sumVheMol + sumVhehe
                  EE = vv

                  if(nimps.eq.1) then


                     do 65 kk = 1,nHe
                        
                        
                        rMol = dsqrt((xcm(i) - xHe(i,kk))**2 + 
     *                       (ycm(i) - yHe(i,kk))**2 + 
     *                       (zcm(i) - zHe(i,kk))**2)
                        
                        x_fdMol(kk) = fdMol(xHe(i,kk),xcm(i),rMol)
     *                       /wvfMol(rMol)
                        y_fdMol(kk) = fdMol(yHe(i,kk),ycm(i),rMol)
     *                       /wvfMol(rMol)
                        z_fdMol(kk) = fdMol(zHe(i,kk),zcm(i),rMol)
     *                       /wvfMol(rMol)
                        
                        x_sdMol(kk) = sdMol(xHe(i,kk),xcm(i),rMol)
     *                       /wvfMol(rMol)
                        y_sdMol(kk) = sdMol(yHe(i,kk),ycm(i),rMol)
     *                       /wvfMol(rMol)
                        z_sdMol(kk) = sdMol(zHe(i,kk),zcm(i),rMol)
     *                       /wvfMol(rMol)

 65                  continue
                     
c----------------->    interactions He-Mol  < --------------------          
                   
                     sum_sd_Mol = 0.d0
                  
                     do 70 j = 1,nHe 
   
                        sd_Mol = 0.d0
                        sd_Mol = x_sdMol(j) + y_sdMol(j) + z_sdMol(j)   

                     sum_sd_Mol = sum_sd_Mol + sd_Mol

 70               continue
                     
c-----------------------------------------------------------------
                     El_He = -1.d0/(2.d0*dmassHe) * sum_sd_Mol
         

c----------------->    interactions Mol-He  < -------------------- 
                     
                     sd_total = 0.d0
                     sd_Mol = 0.d0
                     sum_fd_Mol = 0.d0
                     fd_total = 0.d0
                     do 75 j = 1,nHe

                        fd_Mol = 0.d0
                        sd_Mol = (x_sdMol(j) + y_sdMol(j) + z_sdMol(j))
                        sd_total = sd_total + sd_Mol
                        
                        do 80 jj=1,nHe

                           if(jj.ne.j) then

                              fd_j_jj = (x_fdMol(jj) * x_fdMol(j) + 
     *                             y_fdMol(jj) * y_fdMol(j) + 
     *                             z_fdMol(jj) * z_fdMol(j))
                              
                              fd_Mol = fd_Mol + fd_j_jj 

                           endif
 80                   continue

                        fd_total = fd_total + fd_Mol

 75                  continue
                     
                      El_Mol = -1.d0 / (2.d0 * dmassMol) * 
     *                    (fd_total + sd_total + sum_genMol)

***************************************************************************************
                     
                     El_new = El_He + El_Mol + vv
                     
                     El(i) = El_new
                     
                     EE = (El_new + El_old)/2.d0            

                  endif

********************************************BRANCHING  TERM  ****************************               

                  xA(i)=xA0+xcm(i)
                  yA(i)=yA0+ycm(i)
                  zA(i)=zA0+zcm(i)

                  wx = dexp((ET - EE) * TSV) 
                  icomp = int(wx + ran2(idum))
                  mbranch = min(icomp,3)

                  
                  if(mbranch.eq.0) go to 60
*     
                  if(mbranch.eq.1) then
                     
                     ikount=ikount+1
                     v(ikount)=vv
                     elnew(ikount)=El(i)

                     x_axisXnew(ikount) = x_axisX(i) 
                     x_axisYnew(ikount) = x_axisY(i)  
                     x_axisZnew(ikount) = x_axisZ(i) 
                                          
                     y_axisXnew(ikount) = y_axisX(i)  
                     y_axisYnew(ikount) = y_axisY(i)  
                     y_axisZnew(ikount) = y_axisZ(i) 
                     
                     do j=1,nhe
                        xnHe(ikount,j)=xHe(i,j)
                        ynHe(ikount,j)=yHe(i,j)
                        znHe(ikount,j)=zHe(i,j)
                     enddo

                     xnA(ikount)=xA(i)
                     ynA(ikount)=yA(i)
                     znA(ikount)=zA(i)
                     
                     xncm(ikount)=xcm(i)
                     yncm(ikount)=ycm(i)
                     zncm(ikount)=zcm(i)
                  endif
*     
                  if(mbranch.eq.2) then
                     
                     ikount=ikount+2
                     
                     elnew(ikount-1)=El(i)
                     v(ikount-1)=vv
                     
                     x_axisXnew(ikount-1) = x_axisX(i) 
                     x_axisYnew(ikount-1) = x_axisY(i)  
                     x_axisZnew(ikount-1) = x_axisZ(i) 
                                          
                     y_axisXnew(ikount-1) = y_axisX(i)  
                     y_axisYnew(ikount-1) = y_axisY(i)  
                     y_axisZnew(ikount-1) = y_axisZ(i) 
                     
                     do j=1,nhe
                        xnHe(ikount-1,j)=xHe(i,j)
                        ynHe(ikount-1,j)=yHe(i,j)
                        znHe(ikount-1,j)=zHe(i,j)
                     enddo
                     
                     xnA(ikount-1)=xA(i)
                     ynA(ikount-1)=yA(i)
                     znA(ikount-1)=zA(i)
                     
                     xncm(ikount-1)=xcm(i)
                     yncm(ikount-1)=ycm(i)
                     zncm(ikount-1)=zcm(i)
                     
                     elnew(ikount)=El(i)
                     v(ikount)=vv

                     x_axisXnew(ikount) = x_axisX(i) 
                     x_axisYnew(ikount) = x_axisY(i)  
                     x_axisZnew(ikount) = x_axisZ(i) 
                                          
                     y_axisXnew(ikount) = y_axisX(i)  
                     y_axisYnew(ikount) = y_axisY(i)  
                     y_axisZnew(ikount) = y_axisZ(i) 

                     do j=1,nhe
                        xnHe(ikount,j)=xHe(i,j)
                        ynHe(ikount,j)=yHe(i,j)
                        znHe(ikount,j)=zHe(i,j)
                     enddo
                     
                     xnA(ikount)=xA(i)
                     ynA(ikount)=yA(i)
                     znA(ikount)=zA(i)
                     
                     xncm(ikount)=xcm(i)
                     yncm(ikount)=ycm(i)
                     zncm(ikount)=zcm(i)
                  endif
     
                  if(mbranch.eq.3) then
                     
                     ikount=ikount+3
                     
                     elnew(ikount)=El(i)
                     v(ikount)=vv

                     x_axisXnew(ikount) = x_axisX(i) 
                     x_axisYnew(ikount) = x_axisY(i)  
                     x_axisZnew(ikount) = x_axisZ(i) 
                                          
                     y_axisXnew(ikount) = y_axisX(i)  
                     y_axisYnew(ikount) = y_axisY(i)  
                     y_axisZnew(ikount) = y_axisZ(i) 

                     
                     do j=1,nhe
                        xnHe(ikount,j)=xHe(i,j)
                        ynHe(ikount,j)=yHe(i,j)
                        znHe(ikount,j)=zHe(i,j)
                     enddo

                     xnA(ikount)=xA(i)
                     ynA(ikount)=yA(i)
                     znA(ikount)=zA(i)
                     
                     xncm(ikount)=xcm(i)
                     yncm(ikount)=ycm(i)
                     zncm(ikount)=zcm(i)
                     
                     elnew(ikount-1)=El(i)
                     v(ikount-1)=vv

                     x_axisXnew(ikount-1) = x_axisX(i) 
                     x_axisYnew(ikount-1) = x_axisY(i)  
                     x_axisZnew(ikount-1) = x_axisZ(i) 
                                          
                     y_axisXnew(ikount-1) = y_axisX(i)  
                     y_axisYnew(ikount-1) = y_axisY(i)  
                     y_axisZnew(ikount-1) = y_axisZ(i) 

                     do j=1,nhe
                        xnHe(ikount-1,j)=xHe(i,j)
                        ynHe(ikount-1,j)=yHe(i,j)
                        znHe(ikount-1,j)=zHe(i,j)
                     enddo

                     xnA(ikount-1)=xA(i)
                     ynA(ikount-1)=yA(i)
                     znA(ikount-1)=zA(i)
                     
                     xncm(ikount-1)=xcm(i)
                     yncm(ikount-1)=ycm(i)
                     zncm(ikount-1)=zcm(i)
                     
                     elnew(ikount-2)=El(i)
                     v(ikount-2)=vv
                     
                     x_axisXnew(ikount-2) = x_axisX(i) 
                     x_axisYnew(ikount-2) = x_axisY(i)  
                     x_axisZnew(ikount-2) = x_axisZ(i) 
                                          
                     y_axisXnew(ikount-2) = y_axisX(i)  
                     y_axisYnew(ikount-2) = y_axisY(i)  
                     y_axisZnew(ikount-2) = y_axisZ(i) 

                     do j=1,nhe
                        xnHe(ikount-2,j)=xHe(i,j)
                        ynHe(ikount-2,j)=yHe(i,j)
                        znHe(ikount-2,j)=zHe(i,j)
                     enddo

                     xnA(ikount-2)=xA(i)
                     ynA(ikount-2)=yA(i)
                     znA(ikount-2)=zA(i)
                     
                     xncm(ikount-2)=xcm(i)
                     yncm(ikount-2)=ycm(i)
                     zncm(ikount-2)=zcm(i)
                     
                  endif
                  
 60               continue
               enddo
               
               nmax = ikount
               n1=nmax
               
               eltsum=0.d0
               vsum=0.d0
               
               do jj=1,nmax
                  
                  do j=1,nhe
                     xHe(jj,j)=xnHe(jj,j)
                     yHe(jj,j)=ynHe(jj,j)
                     zHe(jj,j)=znHe(jj,j)
                     
                     if(k.gt.NOTS-40000) then
                     
                     if(mod(k,10000).eq.0)then
                        write(500+j,*) xHe(jj,j),yHe(jj,j),zHe(jj,j)
                     endif
                     endif
                  enddo
                  
                  x_axisX(jj) = x_axisXnew(jj)
                  x_axisY(jj) = x_axisYnew(jj)
                  x_axisZ(jj) = x_axisZnew(jj)
                  
                  y_axisX(jj) = y_axisXnew(jj)
                  y_axisY(jj) = y_axisYnew(jj)
                  y_axisZ(jj) = y_axisZnew(jj)
                  
                  xA(jj)=xnA(jj)
                  yA(jj)=ynA(jj)
                  zA(jj)=znA(jj)
                  
                  xcm(jj)=xncm(jj)
                  ycm(jj)=yncm(jj)
                  zcm(jj)=zncm(jj)

                     if(k.gt.NOTS-40000) then
                    
                     if(mod(k,10000).eq.0)then
                     write(300,*) xA(jj),yA(jj),zA(jj)
                     write(400,*) xcm(jj),ycm(jj),zcm(jj)
                     endif

                     endif

                  vsum = vsum+v(jj)
                  el(jj) = Elnew(jj)
                  
                  eltsum = eltsum + el(jj)
                  
               enddo
                
*************************************************************************************************************
           
               if (k.eq.NOTS.and.neq.eq.1.and.ntype.eq.0) then
                  do kp=1,nmax
                     
                     write(100+nhe,*) xcm(kp),ycm(kp),zcm(kp),xA(kp),
     /                    yA(kp),zA(kp),
     /                    x_axisX(kp),x_axisY(kp),x_axisZ(kp),
     /                    y_axisX(kp),y_axisY(kp),y_axisZ(kp),     
     /                    (xhe(kp,nhe1),yhe(kp,nhe1),zhe(kp,nhe1),
     /                    nhe1=1,nhe)
                  enddo
               endif
               
               close(100+nhe)
               
**********************************************************CALCULATING ENERGY   ****************************************  
               
               term = dfloat(n1)/dfloat(n0)
               alpha = 1.d0/TSV
               
               if(nimps.eq.0)then
                  
                  er = vsum/dfloat(nmax) + alpha*(1.d0-term)

               else

                  er = eltsum/dfloat(nmax) + alpha*(1.d0-term)

               endif

               et = (et+er)/2.d0
               e_wn(k) = er*har2cm
               
               if(l.eq.NOR) then
                  write(20+nhe,*)k,e_wn(k)
               endif
              
*********************************************************RUNNING  AVERAGE    ******************************************            
                  sumE=0.d0
                  avEr=0.d0
          
                  if (k.gt.NOTS-50000) then 
                     
                     do jj = NOTS-50000,k
                        sumE = sumE + e_wn(jj)
                     enddo
                     
                     avEr = sumE
                     avEr = avEr/dfloat(k-(NOTS-50000))
                     
                  endif
                  
                  if (l.eq.NOR.and.k.gt.NOTS-50000) then
                     
                     write(40+nHe,*)k,avEr
                     
                  endif
                  
                  avE(k) = avEr

                  if(nHe.lt.6) then
                     
                     if(mod(k,100).eq.0)then
                        
                        write(6,85) l, k, Nhe, e_wn(k), ave(k),nmax                       
                        
 85              format(2X,I2,2X,I6,5X,'#ofHe = ',I2,5X,'ref.energy = ',
     *               F9.4,7X,'runave = ',F9.4,7X,'number of wlkrs= ',I5)
                     endif
                     
                  else
                     if (mod(k,10).eq.0)then
                        
                        write(6,90) l, k, Nhe, e_wn(k), ave(k),nmax
                        
 90              format(2X,I2,2X,I6,5X,'#ofHe = ',I2,5X,'ref.energy = ',
     *               F9.4,7X,'runave = ',F9.4,7X,'number of wlkrs= ',I5)
                        
                     endif
                  endif
                  
*********************************************************AVERAGING AVERAGES   *******************************************
                  
                  sumavE = 0.d0
                  aveavEr = 0.d0
                  
                  if (k.eq.NOTS) then
                     
                     do jjj = k-29999,k
                        sumavE = sumavE+avE(jjj)
                  enddo
                  
                  aveavEr = sumavE/3.d4
                  
               endif
              
            enddo
               
               write(60+nHe,95) nHe,aveavEr
 95            format(I2,1X,d20.10)
               
               write(11,96)nHe,aveavEr
 96              format(I2,1X,d20.10)
            
            if(NOR.lt.3)then
               
               write(6,*)'average energy ',aveavEr
               write(6,*)'energy/He=',aveaver/dfloat(nhe)
               
            endif
              
         enddo

         close(20+nHe)
         close(40+nHe)
         close(60+nHe)

         close(300)
         close(400)
         do j=1,nhe
            close(500+j)
          enddo

          jj=0
 800   jj=jj+1
      
      read(300,*,end=1000) xxA,yyA,zzA
      read(400,*,end=1000) xxcm,yycm,zzcm
      
      do j=1,nhe    
         
         read(500+j,*) xxHe,yyHe,zzHe
         
         xHel = xxHe-xxcm
         yHel = yyHe-yycm
         zHel = zzHe-zzcm

         xaMol=  xxA-xxcm
         yaMol = yyA-yycm
         zAMol = zzA-zzcm

         phiX = 0.d0
         phiY = 0.d0
         phiZ = -(pi/2.d0-datan2(yAMol,xAMol))

         call rotate(xHel,yHel,zHel,
     *        xHel1,yHel1,zHel1,phiX,phiY,phiZ)  

         r=dsqrt(xAmol**2+yAMol**2+zAMol**2)

         phiX = -dacos(zAmol/r)
         phiY = 0.d0
         phiZ = 0.d0
         
         call rotate(xHel1,yHel1,zHel1,
     *        xHel2,yHel2,zHel2,phiX,phiY,phiZ) 
      
      write(600+nHe,*)xHel2,yHel2,zHel2
      write(700+j,*)xHel2,yHel2,zHel2

      enddo
      
      go to 800
 1000  continue

********************************************  AVERAGED AVERAGE  OVER ALL RUNS  ********************************       
         aveavE(l) = aveavEr  
        
         if(NOR.gt.1)then
            i=0
 100     i=i+1
            READ (60+nHe,*,end=105)nHe,xray(i)
            go to 100
 105        continue
            

            nend = i-1
            sum = 0.d0
            do j = 1,nend
               sum = sum + xray(j)
            enddo             
            
            aven = sum/dfloat(nend) 
            
            sum1 = 0.d0
            do k = 1,nend
               sum1 = sum1 + (aven - xray(k))**2  
            enddo
            std=dsqrt(sum1/dfloat(nend-1))
            
            if (NOR.gt.2)then

               write(6,*)'average energy =' ,aven,' +/-',std
               write(6,*)'energy/He=',aven/dfloat(nhe)
               
               write(6,110)aven,std               
 110           format(d20.10,1X,d20.10)


               if (ntype.eq.0)then
                  write(15,115)aven,std               
 115              format(d20.10,1X,d20.10)
               else
                  write(16,120)aven,std               
 120              format(d20.10,1X,d20.10)
               endif

            endif
            
         endif 
         
         write(80+nHe,*)'#of He ',NHe,'INW =',INW,'TSV =',TSV
         write(80+nHe,*)'NOTS =',NOTS,'NOR =',NOR,'  B =',B_wn
         write(80+nHe,*) 

         if (nimps.eq.0)then
            write(80+nHe,*) 'no importance sampling'
            write(6,*) 'no importance sampling'
         else
            write(80+nHe,*) 'importance sampling'
            write(6,*) 'importance sampling'
         endif

         if (ntype.eq.0) then

            write (80+nHe,*) 'gr_state'
            write (6,*) 'gr_state'

         elseif (ntype.eq.1) then 
            
            write(80+nHe,*)'free_rotor'
            write(6,*)'free_rotor'
            
         elseif (ntype.eq.2) then 
            
            write(80+nHe,*)'cos_teta'
            write(6,*)'cos_teta'
            
         elseif (ntype.eq.3) then 
            
            write(80+nHe,*)'general_node'
            write(6,*)'general_node'

         elseif (ntype.eq.4) then 
            
            write(80+nHe,*)'adiabatic_node'
            write(6,*)'adiabatic_node'
            
         endif

            if (NOR.gt.2)then
               
               write(80+nHe,*)'final average =',aven, '+/-', std
               write(80+nHe,*)'energy/He=',aven/dfloat(nhe)
               
            else
               
            write(80+nHe,*)'final average =',aveaver
            write(80+nHe,*)'energy/He=',aveaver/dfloat(nhe)

         endif
         
         write(6,125) nhe,NOR
 125     format (I2,1x,'He-HX',3x,'NOR',1x,I2)
 
      enddo
      end
      
**************************************************SUBROUTINES  AND FUNCTIONS *****************************  
      subroutine pothehe(vv1,rhe) !he-he potential
      real*8 r,vv1,rm,rhe,x,vst1,vst2,fx,vst,econ,alpha,beta,
     /     aa,d,eps,c6,c8,c10
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
      end 
      

      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836)
      Parameter (NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END	
				
      FUNCTION gasdev(idum)
      INTEGER idum
      REAL*8 gasdev
      INTEGER iset
      REAL*8 fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
 1      v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1.0.or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      end

      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1

          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END

      subroutine rotate(xA,yA,zA,xnewA,ynewA,znewA,phiX,phiY,phiZ)
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension rb(3,3),rc(3,3),rd(3,3),a(3,3)

      rd(1,1) = 1.d0
      rd(1,2) = 0.d0
      rd(1,3) = 0.d0
      rd(2,1) = 0.d0
      rd(2,2) = dcos(phiX)
      rd(2,3) = dsin(phiX)
      rd(3,1) = 0.d0
      rd(3,2) = -dsin(phiX) 
      rd(3,3) = dcos(phiX)
      
      rc(1,1) = dcos(phiY)
      rc(1,2) = 0.d0
      rc(1,3) = -dsin(phiY)
      rc(2,1) = 0.d0
      rc(2,2) = 1.d0
      rc(2,3) = 0.d0
      rc(3,1) = dsin(phiY)
      rc(3,2) = 0.d0
      rc(3,3) = dcos(phiY)
      
      rb(1,1) = dcos(phiZ)
      rb(1,2) = dsin(phiZ)
      rb(1,3) = 0.d0
      rb(2,1) = -dsin(phiZ)
      rb(2,2) = dcos(phiZ)
      rb(2,3) = 0.d0
      rb(3,1) = 0.d0
      rb(3,2) = 0.d0
      rb(3,3) = 1.d0

c!se calcula el producto de las tres matrices de rotacion 
         
      a(1,1) = (rb(1,1)*rc(1,1) + rb(1,2)*rc(2,1) 
     *     + rb(1,3)*rc(3,1))*rd(1,1) 
     *     + (rb(1,1)*rc(1,2) + rb(1,2)*rc(2,2) 
     *     + rb(1,3)*rc(3,2))*rd(2,1) 
     *     + (rb(1,1)*rc(1,3) + rb(1,2)*rc(2,3) 
     *     + rb(1,3)*rc(3,3))*rd(3,1)
      
      a(1,2) =  (rb(1,1)*rc(1,1) + rb(1,2)*rc(2,1) 
     *     + rb(1,3)*rc(3,1))*rd(1,2) 
     *     + (rb(1,1)*rc(1,2) + rb(1,2)*rc(2,2) 
     *     + rb(1,3)*rc(3,2))*rd(2,2) 
     *     + (rb(1,1)*rc(1,3) + rb(1,2)*rc(2,3) 
     *     + rb(1,3)*rc(3,3))*rd(3,2)
      
      a(1,3) = (rb(1,1)*rc(1,1) + rb(1,2)*rc(2,1) 
     *     + rb(1,3)*rc(3,1))*rd(1,3) 
     *     + (rb(1,1)*rc(1,2) + rb(1,2)*rc(2,2) 
     *     + rb(1,3)*rc(3,2))*rd(2,3) 
     *     + (rb(1,1)*rc(1,3) + rb(1,2)*rc(2,3)
     *     + rb(1,3)*rc(3,3))*rd(3,3)
      
      
c     ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
      a(2,1) = (rb(2,1)*rc(1,1) + rb(2,2)*rc(2,1) 
     *     + rb(2,3)*rc(3,1))*rd(1,1) 
     *     + (rb(2,1)*rc(1,2) + rb(2,2)*rc(2,2) 
     *     + rb(2,3)*rc(3,2))*rd(2,1) 
     *     + (rb(2,1)*rc(1,3) + rb(2,2)*rc(2,3) 
     *     + rb(2,3)*rc(3,3))*rd(3,1)
      
      a(2,2) = (rb(2,1)*rc(1,1) + rb(2,2)*rc(2,1) 
     *     + rb(2,3)*rc(3,1))*rd(1,2) 
     *     + (rb(2,1)*rc(1,2) + rb(2,2)*rc(2,2) 
     *     + rb(2,3)*rc(3,2))*rd(2,2) 
     *     + (rb(2,1)*rc(1,3) + rb(2,2)*rc(2,3) 
     *     + rb(2,3)*rc(3,3))*rd(3,2)
      
      a(2,3) = (rb(2,1)*rc(1,1) + rb(2,2)*rc(2,1) 
     *     + rb(2,3)*rc(3,1))*rd(1,3) 
     *     + (rb(2,1)*rc(1,2) + rb(2,2)*rc(2,2) 
     *     + rb(2,3)*rc(3,2))*rd(2,3) 
     *     + (rb(2,1)*rc(1,3) + rb(2,2)*rc(2,3) 
     *     + rb(2,3)*rc(3,3))*rd(3,3)
      
      
c     ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^       
      a(3,1) = (rb(3,1)*rc(1,1) + rb(3,2)*rc(2,1) 
     *     + rb(3,3)*rc(3,1))*rd(1,1) 
     *     + (rb(3,1)*rc(1,2) + rb(3,2)*rc(2,2) 
     *     + rb(3,3)*rc(3,2))*rd(2,1) 
     *     + (rb(3,1)*rc(1,3) + rb(3,2)*rc(2,3) 
     *     + rb(3,3)*rc(3,3))*rd(3,1)
         
      a(3,2) = (rb(3,1)*rc(1,1) + rb(3,2)*rc(2,1) 
     *     + rb(3,3)*rc(3,1))*rd(1,2) 
     *     + (rb(3,1)*rc(1,2) + rb(3,2)*rc(2,2) 
     *     + rb(3,3)*rc(3,2))*rd(2,2) 
     *     + (rb(3,1)*rc(1,3) + rb(3,2)*rc(2,3) 
     *     + rb(3,3)*rc(3,3))*rd(3,2)
         
      a(3,3) = (rb(3,1)*rc(1,1) + rb(3,2)*rc(2,1) 
     *     + rb(3,3)*rc(3,1))*rd(1,3) 
     *     + (rb(3,1)*rc(1,2) + rb(3,2)*rc(2,2) 
     *     + rb(3,3)*rc(3,2))*rd(2,3) 
     *     + (rb(3,1)*rc(1,3) + rb(3,2)*rc(2,3) 
     *     + rb(3,3)*rc(3,3))*rd(3,3)
    
      xnewA = a(1,1)*xA + a(1,2)*yA + a(1,3)*zA 
      ynewA = a(2,1)*xA + a(2,2)*yA + a(2,3)*zA   
      znewA = a(3,1)*xA + a(3,2)*yA + a(3,3)*zA
      return
      end

      subroutine rot(xn,yn,zn,xh,yh,zh,phiang)
      implicit real*8(a-h,o-z)

c  rotation from 
c     `

      r=dsqrt(xn**2+yn**2+zn**2)
     
      c=dcos(phiang)
      s=dsin(phiang)
      t=1.d0-c
      
c     unit vector along axis of rotation

      x=xn/r
      y=yn/r
      z=zn/r
      
c     point to be rotated

      r1=xh
      r2=yh
      r3=zh

      xh=c*r1+s*(z*r2-y*r3)+x*t*(x*r1+y*r2+z*r3)
      yh=c*r2+s*(x*r3-z*r1)+y*t*(x*r1+y*r2+z*r3)
      zh=c*r3+s*(y*r1-x*r2)+z*t*(x*r1+y*r2+z*r3)

      return
      end

      subroutine node_type
      implicit real*8(a-h,o-z)
      common/nt/psi,rHe,ntype
      common/param/ c2,c2N,h,nHe 

      if (ntype.eq.1) then
         call cos_beta(psi)
      elseif (ntype.eq.2) then
         call cos_theta(psi)
      elseif (ntype.eq.3) then
         call general_node(psi)
      else
         call adiabatic_node(psi,rHe)
      endif

      return
      end


      subroutine adiabatic_node(psi,rHe)
      implicit real*8(a-h,o-z)
      common/angles/ th(99),ph(99)
      common/param/ c2,c2N,h,nHe   
      th1 = th(1)
      th2 = th(2)
      th3 = th(3)
      th4 = th(4)
      th5 = th(5)
      th6 = th(6)
      th7 = th(7)
      th8 = th(8)

      ph1 = ph(1)
      ph2 = ph(2)
      ph3 = ph(3)
      ph4 = ph(4)
      ph5 = ph(5)
      ph6 = ph(6)
      ph7 = ph(7)
      ph8 = ph(8)

      r=7.8d0 

	if (nhe.eq.1)then

      E=dexp(1.d0)

      


        psi1 = (0.334676599772923d0*Cos(th1))/
     -   E**(0.5d0*(-7.2d0 + r)**2) - 
     -  (0.08400207398708996d0*r*Cos(th1))/
     -   E**(0.5d0*(-7.2d0 + r)**2) + 
     -  (0.0058334773602145805d0*r**2*dcos(th1))/
     -   E**(0.5d0*(-7.2d0 + r)**2) - 
     -  (0.008025591660672151d0*dcos(th2))/
     -   E**(0.5d0*(-7.2d0 + r)**2) + 
     -  (1.0592246092641524d0*dcos(th1)*dcos(th2)**2)/
     -   E**(0.5d0*(-7.2d0 + r)**2) - 
     -  (0.2616488139211456d0*r*dcos(th1)*dcos(th2)**2)/
     -   E**(0.5d0*(-7.2d0 + r)**2) + 
     - (0.018170056522301775d0*r**2*dcos(th1)*dcos(th2)**2)/
     -   E**(0.5d0*(-7.2d0 + r)**2) + 
     -  (1.0594331024489496d0*dcos(ph1 - ph2)*dcos(th2)*
     -     dsin(th1)*dsin(th2))/E**(0.5d0*(-7.2d0 + r)**2) -
     -  (0.2617064633229274d0*r*dcos(ph1 - ph2)*dcos(th2)*
     -     dsin(th1)*dsin(th2))/E**(0.5d0*(-7.2d0 + r)**2) +
     -  (0.01817405995298107d0*r**2*dcos(ph1 -ph2)*dcos(th2)*
     -     dsin(th1)*dsin(th2))/E**(0.5d0*(-7.2d0 + r)**2)
	



	  PSI2=  0.04782817882442769d0*dcos(th1) - 
     -  0.002984155182973038d0*
     -   dcos(th1)**3 - 
     -  0.021696216415799233d0*dcos(th2) + 
     -  0.022084995404497928d0*
     -   dcos(th1)**2*dcos(th2) + 
     -  0.1765825621526455d0*dcos(th1)*
     -   dcos(th2)**2 + 
     -  0.008952465548919113d0*
     -   dcos(th1)**3*dcos(th2)**2 + 
     -  0.007120194266573668d0*
     -   dcos(th2)**3 - 
     -  0.021360582799721*dcos(th1)**2*
     -   dcos(th2)**3 + 
     -  0.014964801137924259d0*
     -   dcos(ph1 - ph2)*dcos(th1)*
     -   dSin(th1)*dSin(th2) + 
     -  0.17956671733561855d0*
     -   dcos(ph1 - ph2)*dcos(th2)*
     -   dSin(th1)*dSin(th2) + 
     -  0.011936620731892151d0*
     -   dcos(ph1 - ph2)*dcos(th1)**2*
     -   dcos(th2)*dSin(th1)*dSin(th2) - 
     -  0.028480777066294673*
     -   dcos(ph1 - ph2)*dcos(th1)*
     -   dcos(th2)**2*dSin(th1)*dSin(th2)
     -   + 0.002984155182973038d0*
     -   dcos(2*ph1 - 2*ph2)*dcos(th1)*
     -   dSin(th1)**2*dSin(th2)**2 - 
     -  0.007120194266573668d0*
     -   dcos(2*ph1 - 2*ph2)*dcos(th2)*
     -   dsin(th1)**2*dsin(th2)**2

	
  	psi1 = 0.04228304444914539d0*dcos(th1) - 
     -  0.010172018120532063*dcos(th2) + 
     -  0.18815263105180471d0*dCos(th1)*
     -   dCos(th2)**2 + 
     -  0.18815263105180471d0*
     -   dCos(ph1 - ph2)*dCos(th2)*
     -   dSin(th1)*dSin(th2)

        psi1 = 0.04150749431575186d0*dCos(th1) + 
     -  0.005371479329351467d0*
     -   dCos(th1)**3 - 
     -  0.019504527200276664d0*dCos(th2) + 
     -  0.02468955386670305d0*dCos(th1)**2*
     -   dCos(th2) + 
     -  0.19407670249461245d0*dCos(th1)*
     -   dCos(th2)**2 - 
     -  0.0161144379880544d0*dCos(th1)**3*
     -   Cos(th2)**2 + 
     -  0.008534683823302887d0*
     -   Cos(th2)**3 - 
     -  0.025604051469908663d0*
     -   dCos(th1)**2*dCos(th2)**3 + 
     -  0.01615487004340016d0*
     -   Cos(ph1 - ph2)*dCos(th1)*
     -   dSin(th1)*dSin(th2) + 
     -  0.188705223165261d0*
     -   dCos(ph1 - ph2)*dCos(th2)*
     -   dSin(th1)*dSin(th2) - 
     -  0.02148591731740587*
     -   dCos(ph1 - ph2)*dCos(th1)**2*
     -   dCos(th2)*dSin(th1)*dSin(th2) - 
     -  0.03413873529321155d0*
     -   dCos(ph1 - ph2)*dCos(th1)*
     -   dCos(th2)**2*dSin(th1)*dSin(th2)
     -   - 0.005371479329351467d0*
     -   dCos(2*ph1 - 2*ph2)*dCos(th1)*
     -   dSin(th1)**2*dSin(th2)**2 - 
     -  0.008534683823302887d0*
     -   dCos(2*ph1 - 2*ph2)*dCos(th2)*
     -   dSin(th1)**2*dSin(th2)**2

c morse now (2ND JUNe)

	psi=0.13783222385544802d0*Cos(th1) + 
     -  (0.03956749650218346d0*Cos(th1))/
     -   2.718281828459045d0**
     -    (1.56d0*(-6.3d0 + r)) - 
     -  (0.07913499300436692d0*Cos(th1))/
     -   2.718281828459045d0**
     -    (0.78d0*(-6.3 + r)) + 
     -  0.00030598753695909456d0*
     -   Cos(th2) + 
     -  (0.015943053333359672d0*Cos(th2))/
     -   2.718281828459045d0**
     -    (2.2d0*(-7.35 + r)) - 
     -  (0.031886106666719344d0*Cos(th2))/
     -   2.718281828459045d0**
     -    (1.1d0*(-7.35d0 + r)) + 
     -  (-0.001689999999999997d0 + 
     -     0.03839d0/
     -      2.718281828459045d0**
     -       (1.9d0*(-7.55d0 + r)) - 
     -     0.07678d0/
     -      2.718281828459045d0**
     -       (0.95d0*(-7.55d0 + r)))*
     -   (0.0010200000000000209d0 + 
     -     0.68318d0/
     -      2.718281828459045d0**
     -       (1.1d0*(-6.7d0 + r)) - 
     -     1.36636d0/
     -      2.718281828459045d0**
     -       (0.55d0*(-6.7d0 + r)))*
     -   (0.09746210015420952d0*
     -      Cos(th2) - 
     -     0.29238630046262853d0*
     -      Cos(th1)**2*Cos(th2) - 
     -     0.29238630046262853d0*
     -      Cos(ph1 - 1.*ph2)*Cos(th1)*
     -      Sin(th1)*Sin(th2))*
     -   (0.09746210015420952d0*
     -      Cos(th1) - 
     -     0.29238630046262853d0*Cos(th1)*
     -      Cos(th2)**2 - 
     -     0.29238630046262853d0*
     -      Cos(ph1 - 1.*ph2)*Cos(th2)*
     -      Sin(th1)*Sin(th2)) + 
     -  (-0.06483d0 - 
     -     0.03217d0/
     -      2.718281828459045d0**
     -       (2.2d0*(-7.5d0 + r)) + 
     -     0.06434d0/
     -      2.718281828459045d0**
     -       (1.1d0*(-7.5d0 + r)))*
     -   (0.17904931097838225d0*
     -      Cos(th2) - 
     -     0.5371479329351468d0*
     -      Cos(th1)**2*Cos(th2) - 
     -     0.2984155182973038d0*
     -      Cos(th2)**3 + 
     -     0.8952465548919113d0*
     -      Cos(th1)**2*Cos(th2)**3 - 
     -     0.238732414637843d0*
     -      Cos(ph1 - 1.*ph2)*Cos(th1)*
     -      Sin(th1)*Sin(th2) + 
     -     1.1936620731892151d0*
     -      Cos(ph1 - 1.*ph2)*Cos(th1)*
     -      Cos(th2)**2*Sin(th1)*
     -      Sin(th2) + 
     -     0.2984155182973038d0*
     -      Cos(2.*ph1 - 2.*ph2)*
     -      Cos(th2)*Sin(th1)**2*
     -      Sin(th2)**2)


      psi2=0.07529307994818765d0*dcos(th1) + 
     -  (0.06650131479822179d0*dcos(th1))/
     -   E**(1.*(-6.7d0 + r)) - 
     -  (0.13300262959644357d0*dcos(th1))/
     -   E**(0.5d0*(-6.7d0 + r)) + 
     -  (0.04011055546417393d0*dcos(th1))/
     -   E**(1.6d0*(-6.25 + r)) - 
     -  (0.08022111092834786d0*dcos(th1))/
     -   E**(0.8d0*(-6.25 + r)) + 
     -  0.03162009047467833d0*dcos(th2) + 
     -  (0.01581349104293555*dcos(th2))/
     -   E**(2.2d0*(-7.35d0 + r)) - 
     -  (0.0316269820858711d0*dcos(th2))/
     -   E**(1.1d0*(-7.35d0 + r)) - 
     -  0.39871257862586346d0*dcos(th1)*
     -   dcos(th2)**2 - 
     -  (0.19950394439466534d0*dcos(th1)*
     -     dcos(th2)**2)/
     -   E**(1.*(-6.7d0 + r)) + 
     -  (0.3990078887893307d0*dcos(th1)*
     -     dcos(th2)**2)/
     -   E**(0.5*(-6.7d0 + r)) - 
     -  0.39871257862586346d0*
     -   dcos(ph1 - ph2)*dcos(th2)*
     -   dSin(th1)*dSin(th2) - 
     -  (0.19950394439466534d0*
     -     dcos(ph1 - ph2)*dcos(th2)*
     -     dSin(th1)*dSin(th2))/
     -   E**(1.d0*(-6.7d0 + r)) + 
     -  (0.3990078887893307d0*
     -     dcos(ph1 - ph2)*dcos(th2)*
     -     dSin(th1)*dSin(th2))/
     -   E**(0.5d0*(-6.7d0 + r))



       
c     morse before	
	 psi2=0.07529307994818765d0*Cos(th1) + 
     -  (0.06650131479822179d0*Cos(th1))/
     -   E**(1.d0*(-6.7d0 + r)) - 
     -  (0.13300262959644357d0*Cos(th1))/
     -   E**(0.5d0*(-6.7d0 + r)) + 
     -  (0.04011055546417393d0*Cos(th1))/
     -   E**(1.6d0*(-6.25d0 + r)) - 
     -  (0.08022111092834786d0*Cos(th1))/
     -   E**(0.8d0*(-6.25d0 + r)) - 
     -  0.00024396303622414297d0*
     -   Cos(th2) - 
     -  (0.00011853571251568529d0*
     -     Cos(th2))/E**(2.1*(-7.7 + r))
     -    + (0.00023707142503137058d0*
     -     Cos(th2))/
     -   E**(1.05d0*(-7.7d0 + r)) - 
     -  0.39871257862586346d0*Cos(th1)*
     -   Cos(th2)**2 - 
     -  (0.19950394439466534d0*Cos(th1)*
     -     Cos(th2)**2)/
     -   E**(1.d0*(-6.7d0 + r)) + 
     -  (0.3990078887893307d0*Cos(th1)*
     -     Cos(th2)**2)/
     -   E**(0.5d0*(-6.7d0 + r)) - 
     -  0.39871257862586346d0*
     -   Cos(ph1 - ph2)*Cos(th2)*
     -   Sin(th1)*Sin(th2) - 
     -  (0.19950394439466534d0*
     -     Cos(ph1 - ph2)*Cos(th2)*
     -     Sin(th1)*Sin(th2))/
     -   E**(1.d0*(-6.7d0 + r)) + 
     -  (0.3990078887893307d0*
     -     Cos(ph1 - ph2)*Cos(th2)*
     -     Sin(th1)*Sin(th2))/
     -   E**(0.5d0*(-6.7d0 + r))

	endif
      return

      end

      subroutine aprox2_spacenode(psi)
      implicit real*8(a-b,d-h,o-z)
      common/angles/ th(99),ph(99)
      
      th1 = th(1)
      th2 = th(2)
      th3 = th(3)
      th4 = th(4)
      th5 = th(5)
      th6 = th(6)
      th7 = th(7)
      th8 = th(8)
      
      ph1 = ph(1)
      ph2 = ph(2)
      ph3 = ph(3)
      ph4 = ph(4)
      ph5 = ph(5)
      ph6 = ph(6)
      ph7 = ph(7)
      ph8 = ph(8)

      end 

      
      subroutine general_node(psi)
      implicit real*8(a-h,o-z)
      common/angles/ th(99),ph(99)
      common/param/ c2,c2N,h,nHe 

      c2N=dsqrt((1.d0-c2**2)/dfloat(nhe))
      cosheliums=0.d0

      do k=1,nhe  
         cosheliums=cosheliums+dcos(th(k+1))
      enddo

      psi=c2*dcos(th(1))+c2N*cosheliums
      end

       subroutine cos_beta(psi)
      implicit real*8(a-h,o-z)
      common/angles/ th(99),ph(99)
c     common/param/ c2,c2N,h,nHe 

      psi=dcos(th(1))

      end

      subroutine cos_theta(psi)
      implicit real*8(a-h,o-z)
      common/angles/ th(99),ph(99)
c      common/param/ c2,c2N,h,nHe 

      psi=dcos(th(2))

      end

      subroutine pothemol(V,r2,x)
      implicit none
      integer l
      double precision  r2, x
      double precision fa, fb, Vsh, Vas, bohr_A,
     &     cm_H, c(6:7,1:3), R,S,T,F,suma,V,r1
C     
      c(6,3)=-0.00313099d0
      c(6,1)=-0.00566357d0
      Bohr_A=0.52917725d0
      cm_H=4.55634d-6
C     calculamos f_6 e f_7      
      fa=1-dexp(S(x)*(r2*Bohr_A))*suma(r2,x)
c     print*,'fa',fa
      fb=1-dexp(S(x)*(r2*Bohr_A))*
     &     (suma(r2,x)+((dabs(S(x)*(r2*Bohr_A)))**7)/5040) 
c     print*,'fb=',fb
      Vsh=T(r2,x)*dexp(R(x)+S(x)*(r2*Bohr_A))
c     print*,'Vsh=',Vsh
      Vas=(fa*c(6,1))/((r2*Bohr_A)**6)
     &     +(((fa*c(6,3))/((r2*Bohr_A)**6)))*F(2,x)
c     print*,'Vas=',Vas
      V=(Vsh+Vas)*cm_H
c     hehfpot=(Vsh+Vas)*cm_H
      end
C     
C     P(n,x)
C     Polinomio de Legendre de orden n en x
C     
      double precision function F(l,x)
      implicit none
      integer l
      double precision x
      if (l .eq. 0) then
         F=1d0
      else if (l .eq. 1) then
         F=x
      else if (l .eq. 2) then
         F=0.5d0*(3d0*(x**2)-1d0)
      else if (l .eq. 3) then
         F=0.5d0*(5d0*(x**3)-3d0*x)
      else if (l .eq. 4) then
         F=(1d0/8d0)*(35d0*(x**4)-30d0*(x**2)+3d0)
      else if (l .eq. 5) then
         F=(1d0/8d0)*(63d0*(x**5)-70d0*(x**3)+15d0*x)
      else if (l .eq. 6) then
         F=(1d0/16d0)*(231d0*(x**6)-315d0*(x**4)+105d0*(x**2)-5d0)
      else if (l .eq. 7) then
         F=(1d0/16d0)*(429d0*(x**7)-693d0*(x**5)+315d0*(x**3)-35d0*x)
      else if (l .eq. 8) then
         F=(1d0/128d0)*(6435d0*(x**8d0)-12012d0*(x**6d0)
     &        +6930d0*(x**4d0)-1260d0*(x**2d0)+35d0)
c     print *,'F(',l,',',x,')=',F
      endif
      end
C     
C     Aqui calculamos e gardamos B(x), D(x) y G(x)
      double precision function R(x)
      implicit none
      integer y
      double precision  b(1:6),x,f
      R=0.d0
      b(1)=13.32010238d0
      b(2)=0.34286792d0
      b(3)=0.29816363d0
      b(4)=0.29425822d0
      b(5)=-0.00934684d0
      b(6)=-0.03540720d0
      do y = 1,6 
         R=R+b(y)*(F(y-1,x))
c     print*,'R(',x,')=',R
      enddo
      end
c     
      double precision function S(x)
      implicit none
      double precision d(1:6),x,f
      integer l 
      S=0d0
      d(1)=-2.10089330d0
      d(2)=0.01023381d0
      d(3)=-0.12817719d0
      d(4)=-0.01081582d0
      d(5)=0.02043170d0
      d(6)=0.00698361d0
      do l = 1,6
         S=S+D(l)*(F(l-1,x))
c     print*,'S(',x,')=',S
      enddo
      end
c     
      double precision function T(r2,x)
      implicit none
      double precision g(1:6,1:6),x,r2,bohr_A,f
      integer l 
      T=0.0d0
      g(1,1)=4.22760910d0
      g(1,2)=-0.22438673d0
      g(1,3)=1.58886495d0
      g(1,4)=0.31882907d0
      g(1,5)=0.39685627d0
      g(1,6)=0.17984821d0
      g(2,1)=-2.92257112d0
      g(2,2)=0.25162254d0
      g(2,3)=-1.21660466d0
      g(2,4)=-0.11149131d0
      g(2,5)=-0.26753161d0
      g(2,6)=-0.08211105d0
      g(3,1)=0.68212733d0
      g(3,2)=-0.08245550d0
      g(3,3)=0.31177458d0
      g(3,4)=0.00000000d0
      g(3,5)=0.05921384d0
      g(3,6)=0.00955425d0
      g(4,1)=-0.05598598d0
      g(4,2)=0.00891702d0
      g(4,3)=-0.02776035d0
      g(4,4)=0.00250968d0
      g(4,5)=-0.00434989d0
      g(4,6)=0.000000000d0
      Bohr_A= 0.52917725d0
      do l = 1,6
         T=T+(g(1,l)+g(2,l)*(r2*Bohr_A)+g(3,l)*((r2*Bohr_A)**2)+g(4,l)
     &        *((r2*Bohr_A)**3))*F(l-1,x)
c     print*,'T(',r2,',',x,')=',T
      enddo
      end
C     
C     Aqui calculamos suma
C     
      double precision function suma(r2,x)
      implicit none
      double precision r2,bohr_A,x,S,fact(1:7)
      integer k
      Bohr_A= 0.52917725d0
      suma=0d0
      fact(1)=1d0
      fact(2)=1d0
      fact(3)=2d0
      fact(4)=6d0
      fact(5)=24d0
      fact(6)=120d0
      fact(7)=720d0
      do k = 1,7
         suma=suma+((dabs(S(x)*(r2*Bohr_A)))**k)/fact(k)
C     print *, 'suma(',r2,',',x,')= ', suma
      enddo
      end
