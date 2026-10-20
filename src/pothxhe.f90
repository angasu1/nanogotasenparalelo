module pothx
!Module that calculates the potential He-HX with X=F,Cl,Br,CN (depending on the 
!indicator 1,2,3,4 in the subroutine pothemol        
   use types
   implicit none
   save
   integer(ik), parameter::l=5 
   real(rk)::r, theta 
   real(rk), dimension(0:l):: d, b 
   real(rk), dimension (2,2)::c 
   real(rk), dimension (0:3,0:l)::g 
   real(rk), dimension (0:l)::F 
   real(rk),parameter::bo2arvar=0.5291772_rk
   Integer(ik):: fact(0:7) 


    Contains

 subroutine pars_assignment(ind)
         integer(ik)::ind

  fact=(/1._rk,1._rk,2._rk,6._rk ,24._rk,120._rk,720._rk,5040._rk/)
   Select case (ind)
   case(1)
  d=(/ -2.53146116070142_rk,0.07421267974884_rk,-0.08276071628834_rk,-0.01498993894000_rk&
          &,0.04330076241291_rk,-0.03802078606627_rk/)
  b=(/ 13.35140347292370_rk,-0.09056257250910_rk,0.36222620495196_rk,0.10752025627937_rk,&
          &-0.10729093395687_rk,0.09867415096543_rk/)
   g(0,:)=(/4.57186377381132_rk,0.75361611194433_rk,0.72139454377371_rk,0.78271512913456_rk,&
           &0.29669542735975_rk,-0.02516881784584_rk/)
   g(1,:)=(/-3.73235058315367_rk,-0.48875737093263_rk,-0.66202397351030_rk,-0.61187512997162_rk,&
           &-0.19932910113773_rk,0.02744036687153_rk/)
   g(2,:)=(/1.04324447891595_rk,0.09402772497221_rk,0.20171778104736_rk,0.16229578225988_rk,&
           &0.03785716244561_rk,-0.00646574608438_rk/)
   g(3,:)=(/-0.10441433253613_rk,-0.00497315571833_rk,-0.02153063761904_rk,-0.01500658270722_rk,&
           &-0.00145052999200_rk,0.00000000000000_rk/)
   c(1,1)=-0.01401238000863_rk;c(1,2)=0.00793751597814_rk;c(2,1)=0._rk;c(2,2)=0._rk
     case (2)
  b=(/13.32010238_rk,0.34286792_rk,0.29816363_rk,0.29425822_rk,-0.00934684_rk,-0.03540720_rk/)
  d=(/-2.10089330_rk,0.01023381_rk,-0.12817719_rk,-0.01081582_rk,0.02043170_rk,0.00698361_rk/)
  g(0,:)=(/4.22760910_rk,-0.22438673_rk,1.58886495_rk,0.31882907_rk,0.39685627_rk,0.17984821_rk/)
  g(1,:)=(/-2.92257112_rk,0.25162254_rk,-1.21660466_rk,-0.11149131_rk,-0.26753161_rk,-0.08211105_rk/)
  g(2,:)=(/0.68212733_rk,-0.08245550_rk,0.31177458_rk,0.00000000_rk,0.05921384_rk,0.00955425_rk/)
  g(3,:)=(/-0.05598598_rk,0.00891702_rk,-0.02776035_rk,0.00250968_rk,-0.00434989_rk,0.000000000_rk/)
  c(1,1)=-0.00566357d0;c(1,2)=-0.00313099d0;c(2,1)=0.0_rk;c(2,2)=0.0_rk
    case(3)
 d=(/-2.039697201455_rk,-0.029313071712_rk,-0.093216954967_rk,-0.027635075041_rk&
        & ,-0.064847358666_rk,0.062929310818_rk/)
 b=(/13.389610799159_rk,0.524091825535_rk,0.851082514198_rk,0.280259558185_rk,&
         &0.275185409402_rk,-0.224208308256_rk/)
 g(0,:)=(/4.404601847646_rk,0.050846483162_rk,-1.270755306096_rk,0.912749039267_rk&
         &,0.557309473281_rk,-0.076662339377_rk/)
 g(1,:)=(/-2.879668513147_rk,0.031420907580_rk,0.746535010188_rk,-0.531871910461_rk,&
         &-0.345678532980_rk,0.090918816992_rk/)
 g(2,:)=(/0.637568549808_rk,-0.020484094623_rk,-0.148746579592_rk,0.103244893893_rk,&
         &0.075891728285_rk,-0.031802305582_rk/)
 g(3,:)=(/-0.049851078574_rk,0.002751241046_rk,0.010627371755_rk,-0.006629026827_rk,&
         &-0.006071430979_rk,0.003520180843_rk/)
 c(1,1)=-425.239879082696_rk;c(1,2)=-320.525651759542_rk
 c(2,1)=23.467575642035_rk;c(2,2)=29.052464149358_rk
     case(4)
 d=(/-1.69636850_rk,-0.18540393_rk,-0.19384131_rk,0.01302137_rk,0.00718626_rk,0.02570769_rk/)
 b=(/11.90598200_rk,1.33540950_rk,2.23321070_rk,0.31214304_rk,0.06320769_rk,-0.13334401_rk/)
 g(0,:)=(/0.12723862_rk,-0.36813262_rk,-0.18069923_rk,0.00010292_rk,0.24255687_rk,-0.00422193_rk/)
 g(1,:)=(/0.05485662_rk,0.08842465_rk,0.06785421_rk,-0.01495077_rk,-0.11373621_rk,0.02042454_rk/)
 g(2,:)=(/-0.01310749_rk,0.00037885_rk,-0.00098862_rk,0.00251898_rk,0.01929084_rk,-0.00553967_rk/)
 g(3,:)=(/0.00053377_rk,-0.00075759_rk,-0.00058862_rk,0.00001664_rk,-0.00111851_rk,0.00036288_rk/)
 c(1,1)=-19454.79100000_rk;c(1,2)=-9250.99370000_rk
 c(2,1)=-23297.75000000_rk;c(2,2)=-45021.96700000_rk

   end select

   end subroutine pars_assignment

      function V(r,x,indx)              
      real(rk),intent(in)::  r, x       
      integer(ik),intent(in)::indx      
      real(rk):: fn(2), Vsh, Vas, bohr_A,cm_H,V
      real(rk)::x2,x3,x4,x5,bb,dd,gg,xf,ra,r2,r6,r7
      integer(ik)::i                    
!                                       
      ra=r                              
      if (indx.le.3) ra=ra*bo2arvar     
      x2=x*x;x3=x2*x;x4=x2*x2;x5=x4*x   
      F=(/1._rk ,x ,0.5_rk*(3._rk*(x2)-1_rk) ,0.5_rk*(5._rk*(x3)-3._rk*x) ,&
      &(1._rk/8._rk)*(35._rk*(x4)-30._rk*(x2)+3._rk) ,(1._rk/8._rk)*(63._rk*&
      &(x5)-70._rk*(x3)+15._rk*x)/)     
                                        
      bb=bfun()                         
      dd=dfun()                         
      gg=gfun(ra)                       
      Vsh=gg*dexp(bb+dd*ra)             
      xf=dd*ra                          
      fn=f67(xf)                        
                                        
      Vas=0.0_rk                        
      r2=ra*ra;r6=r2*r2*r2;r7=r6*ra     
      Vas=fn(1)*(c(1,1)+c(1,2)*F(2))/R6 
      Vas=vas+fn(2)*(c(2,1)*F(1)+c(2,2)*F(3))/R7
      V=(Vsh+Vas)!*cm_H                 
!     hehfpot=(Vsh+Vas)*cm_H            
      if (indx.eq.4) v=v*219.474625_rk  
      end function V    

      function bfun()
        real(rk)::bfun
        integer(ik)::k

        bfun=0.0_rk
        do k=0,l
          bfun=bfun+b(k)*F(k)
        enddo

        return
       end function bfun

      function dfun()
        real(rk)::dfun
        integer(ik)::k

        dfun=0.0_rk
        do k=0,l
          dfun=dfun+d(k)*F(k)
        enddo

        return
       end function dfun

      function gfun(r)
        real(rk)::gfun,r,r2,r3
        integer(ik)::k

        gfun=0.0_rk
        r2=r*r;r3=r2*r
        do k=0,l
          gfun=gfun+(g(0,k)+g(1,k)*r+g(2,k)*r2+g(3,k)*r3)*f(k)
        enddo

        return
       end function gfun

       function f67(x)
        real(rk)::f67(2),x,abx,x2,x3,x4,x5,x6,x7,suma(2)
        integer(ik)::l

        abx=abs(x);x2=abx*abx;x3=abx*x2;x4=x2*x2;x5=x4*abx;x6=x4*x2;x7=x4*x3
        suma(1)=1.0_rk+abx+x2/fact(2)+x3/fact(3)+x4/fact(4)+x5/fact(5)+x6/fact(6)
        suma(2)=suma(1)+x7/fact(7)
        suma=suma*dexp(x)      
        f67=1.0_rk-suma      


        return
       end function f67



end module pothx
