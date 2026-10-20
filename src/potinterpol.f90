module potinter
!Module that calculates the potential He-HX with X=F,Cl,Br,CN (depending on the 
!indicator 1,2,3,4 in the subroutine pothemol 
!interpolación V_lam*P_lam =V(R,theta)       
   use types
   use global_variables 
   use input_values
   implicit none

   Contains
   subroutine interpol(Vsum,r,cthet)
   !interpolación, encontrar los vlam, plam y sacar el potencial, para nuevos r0 y
       real(rk),intent(out)            :: Vsum
       real(rk),intent(in)             :: r,cthet
       integer(ik)                     :: i, j
       real(rk)                        :: cth1, cth2, r1, r2 
       integer(ik)                     :: linear, lineath
       real(rk)                        :: pl, vl, sloper, slopeth
       
       !tamano del vector, dado por el # de terminos a utilizar 
       nColumna = ntermspot + 2 !ntermspot cuenta los terminos, lamda
     
        if (r.gt.rexit) then             
               vsum=0.0_rk              
               return                   
       endif 

       !encontar el num. de linea donde esta el valor
       linear = 1 + nint((1000._rk*(r-2._rk))/(rmaxpot-2._rk))
       !lineath = 1 + int((1000._rk*(cthet))/pi)
       lineath = 1 - nint(500._rk*(cthet-1._rk))

              
       cth1 = plam(lineath,1)            
       cth2 = plam(lineath+1,1)          
       r1 = vlam(linear,1)              
       r2 = vlam(linear+1,1)            
       
       ! Asegura estar en la linea correcta 
       if(abs(cthet-cth1)>(1._rk/1000._rk).or.abs(r-r1)>((rmaxpot-2._rk)/1000._rk))then
           write(*,*) "no estas en la línea correcta"
           write(*,*)  "angulos", lineath,cthet,cth1,cth2 
           write(*,*) "radio",linear,r,r1,r2
           stop                         
       end if  
       if(linear.gt.nLinea)stop "fuera del rango de theta"
       if(linear.gt.nLinea)stop "fuera del rango de r"
       
       !--------------------
       vsum = 0.0_rk
       slopeth = (cthet-cth1)/(cth2-cth1)
       sloper = (r-r1)/(r2-r1)

       !Encuentra los nuevos vlam, plam para calculo del potencial
       do j = 2, nColumna               
           pl = plam(lineath,j) + slopeth*(plam(lineath+1,j) - plam(lineath,j))
           vl = vlam(linear,j) + sloper*(vlam(linear+1,j) - vlam(linear,j))
           if(j ==  4)then
             vl = vl*1.2_rk
           end if 
           vsum = vsum + pl*vl
       end do !j = 2, nColumna              
                                        
   end subroutine interpol        

end module potinter
