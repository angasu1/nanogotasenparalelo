      block data
       implicit double precision (a-h,o-z)
      parameter (n=1)
      common/coordi/ x(n,3,4)
      common/energy/ e(n)
C approximate minimum structure of G|SO-3
      data x/1.726e0,0.243e0,0.e0,5.849e0,-1.595e0,0.e0,
     *      0.e0,0.e0,0.e0,5.159e0,0.e0,0.e0/
C approximate minimum structure of G|SC-2.9
C      data x/1.725e0,0.249e0,0.e0,5.840e0,-1.602e0,0.e0,
C     *      0.e0,0.e0,0.e0,5.167e0,0.e0,0.e0/
      end

      subroutine hf2pot(rm1var,rm2var,rabvar,th1var,th2var,tauvar,en)
       implicit double precision (a-h,o-z)


C
C potential energy surfaces for (HF)_2 
C coordinates x (common block /coordi/) --> 
C energy e/cm-1 (common block /energy/)
C see below
C
C PLEASE REPORT ANY PROBLEMS TO suhm@ir.phys.chem.ethz.ch 
C M. Suhm 1996/97
C
C checkpoints: (center-of-mass coordinates for SQSBDE and SN(A,B,C),
C               atomic coordinates for SC-2.9, SO-3)
C
C surface r1     r2     R      th1   th2  tau     :     e/cm-1
C
C SQSBDE  1.7445 1.7404 5.144  9.00 64.14 180.00  :  -1559.307
C +M91    2.     1.     4.    80.   70.   110.    : 726656.41
C 1991
C
C SNA+M   1.7408 1.7369 5.089  8.77 65.15 180.00  :  -1545.25
C 1992    2.     1.     4.    80.   70.   110.    : 692948.91
C
C SNB+M   1.7395 1.7356 5.142  8.69 65.00 180.00  :  -1544.21
C 1992    2.     1.     4.    80.   70.   110.    :  97393.86
C
C SNC+M   1.7396 1.7357 5.142  8.69 65.01 180.00  :  -1544.245
C 1992    2.     1.     4.    80.   70.   110.    :  97462.79
C
C SC-2.9  1.743  1.738  5.166  8.23 67.19 180.00  :  -1598.005
C +GPT    2.     1.     4.    80.   70.   110.    : 138468.86
C 1996
C
C SO-3    1.743  1.738  5.158  8.05 66.65 180.00  :  -1597.890
C +GPT    2.     1.     4.    80.   70.   110.    : 138519.45
C 1996
C
C M and M91 are Morse potentials, GPT is a generalized Poeschl Teller 
C potential for HF
C
C coordinates are  expressed in bohr(=a_0) and rad/deg 
C energy outputs are in cm-1
C
C length: 1ao = 52.9177249 pm
C energy: 1Eh = 219474.63067 cm-1
C m(H)/m(F)   = 0.053047888  = 1.007825/18.99840 (CRC Handbook)
C
C on 64bit computers, all 'implicit double precision' statements should
C be removed and all strings  ***d*0***, ***d*+*0*** and ***d*-*0*** 
C occuring in the source code should be replaced by ***e*0***,
C ***e*+0*** and ***e*-*0***, respectively
C (ignore the * above, they protect the text from the suggested change)
C
      !implicit double precision (a-h,o-z)
      parameter (n=1, ma=135) 
C n=number of simultaneous potential calls
C   if n is large, the time - consuming potential evaluation is vectorizable
C   to a very large extent, reducing the CPU time by typically one order
C   of magnitude
C ma=number of parameters
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      parameter (ipotty=0)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C potential flag:                      0=SQSBDE  1=SNA   2=SNB   3=SNC      C
C                                                4=SC-2.9-m 5=SC-2.9-gpt    C
C                                                6=SO-3-m 7=SO-3-gpt        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      parameter (icoor=2,iext=0)
C coordinate flag  0=cartesian  
C                  1=internal (c.o.m only for ipotty<4), angles in rad
C                  2=internal (c.o.m only for ipotty<4), angles in deg
C                    see JCP, 95 (1991), p30 for a definition
C iext: flag for inputs: if 1, input is via common block
C note that parameter nat in hf2p has to be used correctly
C (4 for dimer, 6 for trimer etc.)
      parameter (pi=3.1415926535898e0, phicon=0.01745329252e0)
      common/coordi/ x(n,3,4)
      common/intern/ rm1(n),rm2(n),rab(n),th1(n),th2(n),tau(n),
     1               r12(n),r14(n),r23(n),r34(n)
      common/energy/ e(n)
      common/helg/ dip0,dip1,alpa0,alpa1,alper0,alper1,
     *             rmix1,rmix2,rmix3,rmix4,
     *             gpt1,gpt2,gpt3,gpt4,gpt5
      dimension wellde(0:10)
      data wellde/1559.307e0,1545.25e0,
     1            1544.21e0,1544.242e0,7*0.e0/
C
C subroutine hf2pp establishes parameters and coordinate
C independent elements of the potential. This call has to be
C made only once for an arbitrary number of potential calls!
C
      call hf2pp(ipotty) 
C
      if(icoor.eq.0.and.iext.ne.1)then
      rewind(1)
      do 700 i=1,n
      read(1,*) x(i,1,1),x(i,2,1),x(i,3,1),
     1          x(i,1,2),x(i,2,2),x(i,3,2),
     2          x(i,1,3),x(i,2,3),x(i,3,3),
     3          x(i,1,4),x(i,2,4),x(i,3,4)
  700 continue
      elseif(icoor.ne.0.and.iext.ne.1)then
      rewind(1)
      do 701 i=1,n
      rm1(i)=rm1var
      rm2(i)=rm2var
      rab(i)=rabvar
      th1(i)=th1var
      th2(i)=th2var
      tau(i)=tauvar        
C     write(*,*) 'internal',rm1(i),rm2(i),rab(i),th1(i),th2(i),tau(i)        
C     read(1,*)  rm1(i),rm2(i),rab(i),th1(i),th2(i),tau(i)

        
  701 continue
      if(icoor.eq.2)then
      do 703 i=1,n
      th1(i)=th1(i)*phicon
      th2(i)=th2(i)*phicon
      tau(i)=tau(i)*phicon
  703 continue
      endif
C
C subroutine intcon converts internal to cartesian coordinates
C
      call intcon(ipotty)
      endif
      if(icoor.ne.0.and.iext.eq.1)then
      write(*,*) 'not implemented, would be inefficient'
      stop
      endif
C
C subroutine intcor calculates all necessary internal
C coordinates, namely r1,r2,r,th1,th2,tau,r12,r23,r13,r14
C and scales them for the potential, where appropriate.
C This is done for n sets of input coordinates at a time.
C
      call intcor(ipotty)
C
C subroutine hf2p calculates the potential energy e as a function
C of properly scaled internal coordinates r1...r14 for n
C sets of such coordinates.
C
      call hf2p(ipotty)
C
      if(mod(ipotty,10).le.3)then
      do 902 i=1,n
      e(i)=e(i)-wellde(ipotty)
  902 continue
      endif
C     write(*,*) 'energias',e(1)
      en=e(1)
      return
      end
C
C *** subroutine intcon ***  converts internal coordinates to cartesians
C
      subroutine intcon(ipotty)
       implicit double precision (a-h,o-z)
      !implicit double precision (a-h,o-z)
      parameter (n=1)
      parameter (zero=0.0e0, one=1.0e0, two=2.0e0, three=3.0e0) 
      parameter (ten=10.0e0, tenth=0.1e0, half=0.5e0)
      common/coordi/ x(n,3,4)
      common/energy/ e(n)
      common/intern/ rm1(n),rm2(n),rab(n),th1(n),th2(n),tau(n),
     1               r12(n),r14(n),r23(n),r34(n)
      romh=0.050375571e0
C mass fraction H/(F+H) for HF center of gravity : romh
      if(ipotty.gt.3) romh=zero
C atom-centered coordinates
      do 510 i=1,n

      x(i,1,3)=zero-rm1(i)*cos(th1(i))*romh
      x(i,2,3)=zero-rm1(i)*sin(th1(i))*romh
      x(i,3,3)=zero

      x(i,1,1)=zero+rm1(i)*cos(th1(i))*(one -romh)
      x(i,2,1)=zero+rm1(i)*sin(th1(i))*(one -romh)
      x(i,3,1)=zero

      x(i,1,4)=rab(i)-rm2(i)*cos(th2(i))*romh
      x(i,2,4)=zero  -rm2(i)*sin(th2(i))*cos(tau(i))*romh
      x(i,3,4)=zero  -rm2(i)*sin(th2(i))*sin(tau(i))*romh

      x(i,1,2)=rab(i)+rm2(i)*cos(th2(i))*(one -romh)
      x(i,2,2)=zero  +rm2(i)*sin(th2(i))*cos(tau(i))*(one -romh)
      x(i,3,2)=zero  +rm2(i)*sin(th2(i))*sin(tau(i))*(one -romh)

  510 continue
      return
      end
C
C *** subroutine intcor ***  computes internal coordinates
C
      subroutine intcor(ipotty)
       implicit double precision (a-h,o-z)
      parameter (n=1)
      parameter (zero=0.0e0, one=1.0e0, two=2.0e0, three=3.0e0) 
      parameter (ten=10.0e0, tenth=0.1e0, half=0.5e0)
      parameter (pi=3.1415926535898e0, phicon=0.01745329252e0)
      parameter (avf=0.999999e0)
C helps to avoid acos failures due to round off errors
C set to about 0.999999e0 for 32bit accuracy
C mass fraction H/(F+H) for HF center of gravity 
      common/coordi/ x(n,3,4)
      common/energy/ e(n)
      common/intern/ rm1(n),rm2(n),rab(n),th1(n),th2(n),tau(n),
     1               r12(n),r14(n),r23(n),r34(n)
      common/morcon/ r0,alpha,diss,scal,balo,coup,shi1,shi2
C kinter special for mod(ipotty,10)=4,5,6,7
      common/kinter/ ca(n),th1p(n),th2p(n),taup(n)
      dimension v1(n,3),v2(n,3),va(n,3),vb(n,3),dw(n)
      dimension vab(n,3),vsa(n,3),vsb(n,3),dra(n),drb(n)
      nat = 4
      nact1=1
      nact2=2
      romh=0.050375571e0
      if(ipotty.gt.3) romh=zero
      if(ipotty.gt.3)then
      shi1=zero
      shi2=zero
      endif
C build angles ...p from H-viewpoint
      if(mod(ipotty,10).gt.3)then
      do 410 i=1,n
      r14(i)=sqrt((x(i,1,nact1)-x(i,1,nat/2+nact2))**2+
     1            (x(i,2,nact1)-x(i,2,nat/2+nact2))**2+
     1            (x(i,3,nact1)-x(i,3,nat/2+nact2))**2)
      r23(i)=sqrt((x(i,1,nact2)-x(i,1,nat/2+nact1))**2+
     1            (x(i,2,nact2)-x(i,2,nat/2+nact1))**2+
     1            (x(i,3,nact2)-x(i,3,nat/2+nact1))**2)
      v1(i,1)=x(i,1,nat/2+nact1)-x(i,1,nact1)
      v1(i,2)=x(i,2,nat/2+nact1)-x(i,2,nact1)
      v1(i,3)=x(i,3,nat/2+nact1)-x(i,3,nact1)
      v2(i,1)=x(i,1,nat/2+nact2)-x(i,1,nact2)
      v2(i,2)=x(i,2,nat/2+nact2)-x(i,2,nact2)
      v2(i,3)=x(i,3,nat/2+nact2)-x(i,3,nact2)
      vab(i,1)=x(i,1,nact2)-x(i,1,nact1)
      vab(i,2)=x(i,2,nact2)-x(i,2,nact1)
      vab(i,3)=x(i,3,nact2)-x(i,3,nact1)
      r34(i)=sqrt(vab(i,1)**2+vab(i,2)**2+
     1            vab(i,3)**2)
      rm1(i)=sqrt(v1(i,1)**2+v1(i,2)**2+v1(i,3)**2)
      rm2(i)=sqrt(v2(i,1)**2+v2(i,2)**2+v2(i,3)**2)
      th1p(i)=acos(avf*((v1(i,1)*vab(i,1)+v1(i,2)*vab(i,2)
     1       +v1(i,3)*vab(i,3))/
     1       (r34(i)*rm1(i))))
      th2p(i)=acos(-avf*((v2(i,1)*vab(i,1)+v2(i,2)*vab(i,2)
     1       +v2(i,3)*vab(i,3))/
     1       (r34(i)*rm2(i))))
      th2p(i)=pi-th2p(i)
      vsa(i,1)=v1(i,3)*vab(i,2)-v1(i,2)*vab(i,3)
      vsa(i,2)=v1(i,1)*vab(i,3)-v1(i,3)*vab(i,1)
      vsa(i,3)=v1(i,2)*vab(i,1)-v1(i,1)*vab(i,2)
      vsb(i,1)=v2(i,3)*vab(i,2)-v2(i,2)*vab(i,3)
      vsb(i,2)=v2(i,1)*vab(i,3)-v2(i,3)*vab(i,1)
      vsb(i,3)=v2(i,2)*vab(i,1)-v2(i,1)*vab(i,2)
      dra(i)=sqrt(vsa(i,1)**2+vsa(i,2)**2+vsa(i,3)**2)
      drb(i)=sqrt(vsb(i,1)**2+vsb(i,2)**2+vsb(i,3)**2)
      dw(i)= v2(i,1)*(vab(i,2)*v1(i,3)-vab(i,3)*v1(i,2))+
     1       v2(i,2)*(vab(i,3)*v1(i,1)-vab(i,1)*v1(i,3))+
     2       v2(i,3)*(vab(i,1)*v1(i,2)-vab(i,2)*v1(i,1)) 
      if(dra(i).eq.zero.or.drb(i).eq.zero) then
      taup(i)=zero
      else
      taup(i)=acos(avf*((vsa(i,1)*vsb(i,1)+vsa(i,2)*vsb(i,2)
     1             +vsa(i,3)*vsb(i,3))/
     1             (dra(i)*drb(i))))
      endif
C build normal angles viewed from F
      v1(i,1)=x(i,1,nact1)-x(i,1,nat/2+nact1)
      v1(i,2)=x(i,2,nact1)-x(i,2,nat/2+nact1)
      v1(i,3)=x(i,3,nact1)-x(i,3,nat/2+nact1)
      v2(i,1)=x(i,1,nact2)-x(i,1,nat/2+nact2)
      v2(i,2)=x(i,2,nact2)-x(i,2,nat/2+nact2)
      v2(i,3)=x(i,3,nact2)-x(i,3,nat/2+nact2)
      vab(i,1)=x(i,1,nat/2+nact2)-x(i,1,nat/2+nact1)
      vab(i,2)=x(i,2,nat/2+nact2)-x(i,2,nat/2+nact1)
      vab(i,3)=x(i,3,nat/2+nact2)-x(i,3,nat/2+nact1)
      r34(i)=sqrt(vab(i,1)**2+vab(i,2)**2+
     1            vab(i,3)**2)
      rab(i)=r34(i)
      r12(i)=sqrt((x(i,1,nact1)-x(i,1,nact2))**2+
     1            (x(i,2,nact1)-x(i,2,nact2))**2+
     1            (x(i,3,nact1)-x(i,3,nact2))**2)
      rm1(i)=sqrt(v1(i,1)**2+v1(i,2)**2+v1(i,3)**2)
      rm2(i)=sqrt(v2(i,1)**2+v2(i,2)**2+v2(i,3)**2)
      ca(i)=(v1(i,1)*v2(i,1)+v1(i,2)*v2(i,2)+v1(i,3)*v2(i,3))/
     1 (rm1(i)*rm2(i))
      th1(i)=acos(avf*((v1(i,1)*vab(i,1)+v1(i,2)*vab(i,2)
     1       +v1(i,3)*vab(i,3))/
     1       (r34(i)*rm1(i))))
      th2(i)=acos(-avf*((v2(i,1)*vab(i,1)+v2(i,2)*vab(i,2)
     1       +v2(i,3)*vab(i,3))/
     1       (r34(i)*rm2(i))))
      th2(i)=pi-th2(i)
      vsa(i,1)=v1(i,3)*vab(i,2)-v1(i,2)*vab(i,3)
      vsa(i,2)=v1(i,1)*vab(i,3)-v1(i,3)*vab(i,1)
      vsa(i,3)=v1(i,2)*vab(i,1)-v1(i,1)*vab(i,2)
      vsb(i,1)=v2(i,3)*vab(i,2)-v2(i,2)*vab(i,3)
      vsb(i,2)=v2(i,1)*vab(i,3)-v2(i,3)*vab(i,1)
      vsb(i,3)=v2(i,2)*vab(i,1)-v2(i,1)*vab(i,2)
      dra(i)=sqrt(vsa(i,1)**2+vsa(i,2)**2+vsa(i,3)**2)
      drb(i)=sqrt(vsb(i,1)**2+vsb(i,2)**2+vsb(i,3)**2)
      dw(i)= v2(i,1)*(vab(i,2)*v1(i,3)-vab(i,3)*v1(i,2))+
     1       v2(i,2)*(vab(i,3)*v1(i,1)-vab(i,1)*v1(i,3))+
     2       v2(i,3)*(vab(i,1)*v1(i,2)-vab(i,2)*v1(i,1)) 
      if(dra(i).eq.zero.or.drb(i).eq.zero) then
      tau(i)=zero
      else
      tau(i)=acos(avf*((vsa(i,1)*vsb(i,1)+vsa(i,2)*vsb(i,2)
     1             +vsa(i,3)*vsb(i,3))/
     1             (dra(i)*drb(i))))
      endif
      tau(i)=sign(tau(i),dw(i))
 410  continue
      else
      do 210 i=1,n
      v1(i,1)=x(i,1,nact1)-x(i,1,nat/2+nact1)
      v1(i,2)=x(i,2,nact1)-x(i,2,nat/2+nact1)
      v1(i,3)=x(i,3,nact1)-x(i,3,nat/2+nact1)
      v2(i,1)=x(i,1,nact2)-x(i,1,nat/2+nact2)
      v2(i,2)=x(i,2,nact2)-x(i,2,nat/2+nact2)
      v2(i,3)=x(i,3,nact2)-x(i,3,nat/2+nact2)
      va(i,1)=x(i,1,nat/2+nact1)+romh*v1(i,1)
      va(i,2)=x(i,2,nat/2+nact1)+romh*v1(i,2)
      va(i,3)=x(i,3,nat/2+nact1)+romh*v1(i,3)
      vb(i,1)=x(i,1,nat/2+nact2)+romh*v2(i,1)
      vb(i,2)=x(i,2,nat/2+nact2)+romh*v2(i,2)
      vb(i,3)=x(i,3,nat/2+nact2)+romh*v2(i,3)
      vab(i,1)=vb(i,1)-va(i,1)
      vab(i,2)=vb(i,2)-va(i,2)
      vab(i,3)=vb(i,3)-va(i,3)
      rm1(i)=sqrt(v1(i,1)**2+v1(i,2)**2+v1(i,3)**2)
      rm2(i)=sqrt(v2(i,1)**2+v2(i,2)**2+v2(i,3)**2)
      rab(i)=sqrt(vab(i,1)**2+vab(i,2)**2+vab(i,3)**2)
      th1(i)=acos(avf*((v1(i,1)*vab(i,1)+v1(i,2)*vab(i,2)
     1       +v1(i,3)*vab(i,3))/
     1       (rab(i)*rm1(i))))
      th2(i)=acos(-avf*((v2(i,1)*vab(i,1)+v2(i,2)*vab(i,2)
     1       +v2(i,3)*vab(i,3))/
     1       (rab(i)*rm2(i))))
      th2(i)=pi-th2(i)
      vsa(i,1)=v1(i,3)*vab(i,2)-v1(i,2)*vab(i,3)
      vsa(i,2)=v1(i,1)*vab(i,3)-v1(i,3)*vab(i,1)
      vsa(i,3)=v1(i,2)*vab(i,1)-v1(i,1)*vab(i,2)
      vsb(i,1)=v2(i,3)*vab(i,2)-v2(i,2)*vab(i,3)
      vsb(i,2)=v2(i,1)*vab(i,3)-v2(i,3)*vab(i,1)
      vsb(i,3)=v2(i,2)*vab(i,1)-v2(i,1)*vab(i,2)
      dra(i)=sqrt(vsa(i,1)**2+vsa(i,2)**2+vsa(i,3)**2)
      drb(i)=sqrt(vsb(i,1)**2+vsb(i,2)**2+vsb(i,3)**2)
      tau(i)=acos(avf*((vsa(i,1)*vsb(i,1)+vsa(i,2)*vsb(i,2)
     1             +vsa(i,3)*vsb(i,3))/
     1             (dra(i)*drb(i))))
      dw(i)= v2(i,1)*(vab(i,2)*v1(i,3)-vab(i,3)*v1(i,2))+
     1       v2(i,2)*(vab(i,3)*v1(i,1)-vab(i,1)*v1(i,3))+
     2       v2(i,3)*(vab(i,1)*v1(i,2)-vab(i,2)*v1(i,1)) 
      tau(i)=sign(tau(i),dw(i))
  210 continue
      if(mod(ipotty,10).eq.0)then
      do 207 i=1,n
      r12(i)=sqrt((x(i,1,nact1)-x(i,1,nact2))**2+
     1            (x(i,2,nact1)-x(i,2,nact2))**2+
     2            (x(i,3,nact1)-x(i,3,nact2))**2)
      r14(i)=sqrt((x(i,1,nact1)-x(i,1,nat/2+nact2))**2+
     1            (x(i,2,nact1)-x(i,2,nat/2+nact2))**2+
     2            (x(i,3,nact1)-x(i,3,nat/2+nact2))**2)
      r23(i)=sqrt((x(i,1,nat/2+nact1)-x(i,1,nact2))**2+
     1            (x(i,2,nat/2+nact1)-x(i,2,nact2))**2+
     2            (x(i,3,nat/2+nact1)-x(i,3,nact2))**2)
      r34(i)=sqrt((x(i,1,nat/2+nact1)-x(i,1,nat/2+nact2))**2+
     1            (x(i,2,nat/2+nact1)-x(i,2,nat/2+nact2))**2+
     2            (x(i,3,nat/2+nact1)-x(i,3,nat/2+nact2))**2)
  207 continue
      else
      do 208 i=1,n
      dra(i)=shi1+shi2*((rm1(i)-r0)**2/(1.+(rm1(i)-r0)**2)
     1            +(rm2(i)-r0)**2/(1.+(rm2(i)-r0)**2))
      vab(i,1)=vab(i,1)*dra(i)
      vab(i,2)=vab(i,2)*dra(i)
      vab(i,3)=vab(i,3)*dra(i)
      r12(i)=sqrt((x(i,1,nact1)-x(i,1,nact2)-vab(i,1))**2+
     1            (x(i,2,nact1)-x(i,2,nact2)-vab(i,2))**2+
     2            (x(i,3,nact1)-x(i,3,nact2)-vab(i,3))**2)
      r14(i)=sqrt((x(i,1,nact1)-x(i,1,nat/2+nact2)-vab(i,1))**2+
     1            (x(i,2,nact1)-x(i,2,nat/2+nact2)-vab(i,2))**2+
     2            (x(i,3,nact1)-x(i,3,nat/2+nact2)-vab(i,3))**2)
      r23(i)=sqrt((x(i,1,nat/2+nact1)-x(i,1,nact2)-vab(i,1))**2+
     1            (x(i,2,nat/2+nact1)-x(i,2,nact2)-vab(i,2))**2+
     2            (x(i,3,nat/2+nact1)-x(i,3,nact2)-vab(i,3))**2)
      r34(i)=sqrt((x(i,1,nat/2+nact1)-x(i,1,nat/2+nact2)-vab(i,1))**2+
     1            (x(i,2,nat/2+nact1)-x(i,2,nat/2+nact2)-vab(i,2))**2+
     2            (x(i,3,nat/2+nact1)-x(i,3,nat/2+nact2)-vab(i,3))**2)
      rab(i)=rab(i)*(1.e0+dra(i))
  208 continue
      endif
      endif
      return
      end
  
C 
C *** subroutine hf2pp *** prepares computation of hf2 potential
C
      subroutine hf2pp(ipotty)
       implicit double precision (a-h,o-z)
      parameter(ma=135)
      parameter (zero=0.0e0, one=1.0e0, two=2.0e0, three=3.0e0) 
      parameter (ten=10.0e0, tenth=0.1e0, half=0.5e0)
      parameter (pi=3.1415926535898e0, phicon=0.01745329252e0)
      dimension fac(0:8)
      dimension fct(31),efct(31)
      common/morcon/ r0,alpha,diss,scal,balo,coup,shi1,shi2
      common/helg/ dip0,dip1,alpa0,alpa1,alper0,alper1,
     *             rmix1,rmix2,rmix3,rmix4,
     *             gpt1,gpt2,gpt3,gpt4,gpt5
      common/ang/ g(0:4,0:4,0:8,0:4),pf(0:4,0:4)
      common/par/ a(ma),ac(3,0:4,0:4,0:8)
C
C note: this routine supplies g(...)*pif instead of g(...)
C
      pif=half/sqrt(pi)

C
C define parameters for the surface
C
      do 750 j=1,ma
      a(j)=0.e0
  750 continue
      if(mod(ipotty,10).eq.0)then
       a(  1) =  0.1559307e04
       a(  2) = -0.2029e09
       a(  3) = -0.3703e010
       a(  4) = -0.1135e011
       a(  5) = -0.1736693095e07
       a(  6) = -0.4416487722e07
       a(  7) =  0.1601238951e08
       a(  8) = -0.6389444855e07
       a(  9) = -0.866867375551329e03
       a( 10) =  0.129789363470169e06
       a( 13) =  0.262218475665877e07
       a( 17) =  0.184254701843311e012
       a( 18) = -0.32e06
       a( 19) =  0.0332643894386025e0
       a( 20) = -0.00457688720688704e0
       a( 21) =  0.409730666039653e0
       a( 22) =  0.326235519432838e014
       a( 23) =  0.1e0
       a( 24) =  0.717544978942018e014
       a( 26) =  0.0503238787247912e0
       a( 27) =  0.146832057915281e01
       a( 28) =  0.338463155928653e04
       a( 31) = -0.160267101272529e06
       a( 36) = -0.162688946598104e05
       a( 39) =  0.280257610197247e04
       a( 46) =  0.286006394339599e04
       a( 52) =  0.509102022118840e03
       a( 68) = -0.346723417017878e04
       a( 73) =  0.230463285803830e04
       a( 84) = -0.644815083021031e03
       a(101) =  0.225309324112547e05
C       a(101) =  0.21489156e05
       a(102) =  0.163040722103571e04 
C       a(102) =  0.1895154e04 
C --trials with *3.0
       a(104) =  0.104290771467500e04
C       a(104) =  0.0768395e04
       a(105) = -0.366368726725276e03
C       a(105) = -0.395230e03
       a(106) =  0.242817244795167e04
C       a(106) =  0.2041094e04
       a(107) = -0.188237787794149e04
C       a(107) = -0.2422952e04
       a(110) = -0.326002341194829e03
C       a(110) = -0.297373e03
       a(113) =  0.116704136523402e03
C       a(113) =  0.039468e03
       a(115) =  0.572632667044664e03
C       a(115) =  0.720555e03
C --trials with *2.5
       r0=1.7374e0
       alpha=1.1953e0
       diss=47635.e0
       scal=one
       balo=zero
       coup=zero
       shi1=zero
       shi2=zero
      elseif(mod(ipotty,10).le.3.and.mod(ipotty,10).gt.0)then
C ipotty=3 directly, 1+2 indirectly
       a(  1) =  0.156935e04
       a(  2) = -0.2029e09
       a(  3) = -0.3703e010
       a(  4) = -0.1135e011
       a(  5) = -0.173e07
       a(  6) = -0.44e07
       a(  7) =  0.16e08
       a(  8) = -0.64e07
       a(  9) = -0.86687e03
       a( 10) =  0.129789e06
       a( 13) =  0.262219e07
       a( 17) =  0.184255e012
       a( 18) = -0.32e06
       a( 19) =  0.085e0
       a( 20) = -0.005e0
       a( 21) =  0.3e0
       a( 22) =  0.326236e014
       a( 23) =  0.05e0
       a( 24) =  0.71755e014
       a( 26) =  -0.03e0
       a( 27) =  0.155e0
       a( 31) = -0.160267e06
       a( 36) = -0.16269e05
       a( 39) =  0.280258e04
       a( 46) =  0.286006e04
       a( 52) =  0.509102e03
       a( 68) = -0.346723e04
       a( 73) =  0.230463e04
       a( 84) = -0.644815e03
       a(101) =  0.22531e05
       a(102) =  0.163041e04
       a(104) =  0.104291e04
       a(105) = -0.36637e03
       a(106) =  0.2800e04
       a(107) = -0.188238e04
       a(113) =  0.11670e03
       r0=1.7327e0
       alpha=1.16937e0
       diss=49652.7e0
       scal=0.984e0
       balo=9.5e0
       coup=50000.e0
       shi1=zero
       shi2=0.2e0
      elseif(mod(ipotty,10).gt.3.and.mod(ipotty,10).lt.6)then
C ipotty=4/5
CC&&&&&&&&&&&&&&&&&&&&&&&insert here program segment&&&&&&&&&&&&&&&&&&&&&&&
C GPT
C      a(  1)=   0.34075458e+06
C      a(  2)=  -0.16885361e+06
C...
C...      pf(4,4)=  0.42145971e-02
C GPT-MP2
C      a(  1)=   0.35571847e+06
C      a(  2)=  -0.18350144e+06
C...
C...      pf(4,4)=  0.42145971e-02
C GPT-MP2-BSSE
C      a(  1)=   0.35833432e+06
C      a(  2)=  -0.18632209e+06
C...
C...      pf(4,4)=  0.42145971e-02
C GPT-R12-BSSE
C      a(  1)=   0.35211964e+06
C      a(  2)=  -0.18013887e+06
C...
C...      pf(4,4)=  0.42145971e-02
C GPT-S-2.8
C      a(  1)=   0.3407545822e+06
C      a(  2)=  -0.1688536103e+06
C...
C...      pf(4,4)=  0.4214597071e-02
C GPT-SC-2.9
C      a(  1)=   0.3483216956e+06
C      a(  2)=  -0.1787395924e+06
C...
C...      pf(4,4)=  0.4214597071e-02
      a(  1)=   0.3483216956e+06
      a(  2)=  -0.1787395924e+06
      a(  3)=  -0.4098786351e+04
      a(  4)=   0.3703694359e+04
      a(  5)=   0.4822293282e+07
      a(  6)=  -0.4908056944e+06
      a(  7)=  -0.1265708233e+05
      a(  8)=   0.8509877308D+10
      a(  9)=   0.5388595589e+06
      a( 10)=  -0.3890129402e+06
      a( 11)=   0.3212599501e+05
      a( 12)=   0.2222942302e+05
      a( 13)=  -0.1152546844e+05
      a( 14)=   0.1293321969e+05
      a( 15)=   0.4875653619e+00
      a( 16)=   0.2975763926e+01
      a( 17)=  -0.1759533432e+05
      a( 18)=   0.8666478084e+03
      a( 19)=   0.5458168393e+00
      a( 20)=   0.8322852434e+04
      a( 21)=  -0.3936059465e+04
      a( 22)=  -0.2571017810e+07
      a( 23)=   0.4369201459e+06
      a( 24)=  -0.3751233102e+06
      a( 25)=  -0.9206230126e+01
      a( 26)=   0.8337489398e+05
      a( 27)=  -0.4682863351e+00
      a( 28)=   0.2470678922e+03
      a( 29)=  -0.3335775239e+05
      a( 30)=  -0.1201569315e+03
      a( 31)=   0.2152798832e+00
      a( 32)=   0.4994283912e+00
      a( 33)=   0.4319028236e+00
      a( 34)=   0.1521054195e+01
      a( 35)=   0.1075479991e+05
      a( 36)=   0.2138126034e+06
      a( 37)=   0.4560625013e+00
      a( 38)=  -0.3278489799e+05
      a( 39)=   0.2902356377e+01
      a( 40)=   0.2227647777e+04
      a( 41)=  -0.1402968123e+04
      a( 42)=   0.1562701790e+06
      a( 43)=   0.1429711730e+01
      a( 44)=   0.1251708032e+01
      a( 45)=  -0.2154473668e+01
      a( 46)=  -0.1976610764e+03
      a( 47)=  -0.1077811593e+05
      a( 48)=  -0.7958908853e+00
      a( 49)=   0.6529506884e+00
      a( 50)=  -0.9839368327e-01
      a( 51)=   0.1131649414e+01
      a( 52)=   0.1851459145e+01
      a( 53)=   0.2607012711e+00
      a( 54)=  -0.1895397639e+03
      a( 55)=  -0.1791664512e+04
      a( 56)=   0.1471881356e+01
      a( 57)=  -0.6604471722e+05
      a( 58)=   0.1296331624e+02
      a( 59)=   0.1146743598e+03
      a( 60)=   0.4371368533e+05
      a( 61)=   0.6827481531e+00
      a( 62)=   0.0000000000e+00
      a( 63)=   0.4114367187e+03
      a( 64)=  -0.2332580324e+05
      a( 65)=   0.0000000000e+00
      a( 66)=  -0.1789723462e+01
      a( 67)=  -0.1036106025e+01
      a( 68)=   0.3939911825e+00
      a( 69)=   0.1004744860e+06
      a( 70)=   0.0000000000e+00
      rmix1 =  0.2900000000e+01
      rmix2 = -0.1900000000e+01
      rmix3 =  0.0000000000e+00
      rmix4 =  0.0000000000e+00
      dip0  =  0.1847287981e+01
      dip1  =  0.9324918225e+00
      alpa0 =  0.6151261000e+01
      alpa1 =  0.5400010000e+01
      alper0=  0.4793851000e+01
      alper1=  0.1389360000e+01
      g(1,1,0,1)=  0.2000000000e+01
      g(1,1,2,1)= -0.1000000000e+01
      g(1,2,1,1)=  0.1732050808e+01
      g(1,2,3,1)= -0.1154700538e+01
      g(2,1,1,1)=  0.1732050808e+01
      g(2,1,3,1)= -0.1154700538e+01
      g(1,3,2,1)=  0.1632993162e+01
      g(1,3,4,1)= -0.1224744871e+01
      g(2,2,0,1)=  0.2000000000e+01
      g(2,2,0,2)=  0.2000000000e+01
      g(2,2,2,1)=  0.1000000000e+01
      g(2,2,2,2)= -0.2000000000e+01
      g(2,2,4,1)= -0.1333333333e+01
      g(2,2,4,2)=  0.3333333333e+00
      g(3,1,2,1)=  0.1632993162e+01
      g(3,1,4,1)= -0.1224744871e+01
      g(1,4,3,1)=  0.1581138830e+01
      g(1,4,5,1)= -0.1264911064e+01
      g(2,3,1,1)=  0.1885618083e+01
      g(2,3,1,2)=  0.1490711985e+01
      g(2,3,3,1)=  0.7071067812e+00
      g(2,3,3,2)= -0.2236067977e+01
      g(2,3,5,1)= -0.1414213562e+01
      g(2,3,5,2)=  0.4472135955e+00
      g(3,2,1,1)=  0.1885618083e+01
      g(3,2,1,2)=  0.1490711985e+01
      g(3,2,3,1)=  0.7071067812e+00
      g(3,2,3,2)= -0.2236067977e+01
      g(3,2,5,1)= -0.1414213562e+01
      g(3,2,5,2)=  0.4472135955e+00
      g(4,1,3,1)=  0.1581138830e+01
      g(4,1,5,1)= -0.1264911064e+01
      g(2,4,2,1)=  0.1825741858e+01
      g(2,4,2,2)=  0.1290994449e+01
      g(2,4,4,1)=  0.5477225575e+00
      g(2,4,4,2)= -0.2323790008e+01
      g(2,4,6,1)= -0.1460593487e+01
      g(2,4,6,2)=  0.5163977795e+00
      g(3,3,0,1)=  0.2000000000e+01
      g(3,3,0,2)=  0.2000000000e+01
      g(3,3,0,3)=  0.2000000000e+01
      g(3,3,2,1)=  0.1500000000e+01
      g(3,3,2,2)=  0.0000000000e+00
      g(3,3,2,3)= -0.2500000000e+01
      g(3,3,4,1)=  0.3333333333e+00
      g(3,3,4,2)= -0.2333333333e+01
      g(3,3,4,3)=  0.1000000000e+01
      g(3,3,6,1)= -0.1500000000e+01
      g(3,3,6,2)=  0.6000000000e+00
      g(3,3,6,3)= -0.1000000000e+00
      g(4,2,2,1)=  0.1825741858e+01
      g(4,2,2,2)=  0.1290994449e+01
      g(4,2,4,1)=  0.5477225575e+00
      g(4,2,4,2)= -0.2323790008e+01
      g(4,2,6,1)= -0.1460593487e+01
      g(4,2,6,2)=  0.5163977795e+00
      pf(0,0)=  0.2820947918e+00
      pf(1,0)=  0.4886025119e+00
      pf(1,1)= -0.3454941495e+00
      pf(2,0)=  0.6307831305e+00
      pf(2,1)= -0.2575161347e+00
      pf(2,2)=  0.1287580673e+00
      pf(3,0)=  0.7463526652e+00
      pf(3,1)= -0.2154534561e+00
      pf(3,2)=  0.6813236510e-01
      pf(3,3)= -0.2781492158e-01
      pf(4,0)=  0.8462843753e+00
      pf(4,1)= -0.1892349392e+00
      pf(4,2)=  0.4460310290e-01
      pf(4,3)= -0.1192068068e-01
      pf(4,4)=  0.4214597071e-02

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C GPT normal
       r0=1.7327e0
       alpha=1.16937e0
       diss=49652.7e0
       gpt1=0.5847305e0
       gpt2=1.28e0**2
       gpt3=953320.3e0
       gpt4=164766.75e0
       gpt5=49656.989942e0
C GPT MP2
C       r0=1.73916e0
C       alpha=1.16937e0
C       diss=49652.7e0
C       gpt1=0.581e0
C       gpt2=1.286e0**2
C       gpt3=949700.e0
C       gpt4=159771.e0
C       gpt5=49947.015e0
C GPT R12
C       r0=1.7370e0
C       alpha=1.16937e0
C       diss=49652.7e0
C       gpt1=0.5793e0
C       gpt2=1.28e0**2
C       gpt3=951322.5e0
C       gpt4=160317.e0
C       gpt5=50442.629e0
       scal=1.e0
       balo=9.5e0
       coup=50000.e0
       shi1=zero
       shi2=zero
      endif
      if(mod(ipotty,10).eq.1)then
       a(  1) =  0.1571983e04
       a( 19) =  0.145e0
       a( 23) =  0.047e0
       a( 26) =  -0.058e0
       a( 27) =  0.18e0
       scal=0.983e0
       balo=12.e0
       coup=zero
       shi1=0.01e0
       shi2=zero
      elseif(mod(ipotty,10).gt.5.and.mod(ipotty,10).lt.8)then
C ipotty=6/7 (SO-3)
CC&&&&&&&&&&&&&&&&&&&&&&&insert here program segment&&&&&&&&&&&&&&&&&&&&&&&
C GPT-SO-3
      a(  1)=   0.3579200000e+06
      a(  2)=  -0.1839030000e+06
      a(  3)=   0.0000000000e+00
      a(  4)=   0.3917400000e+04
      a(  5)=   0.5087280000e+07
      a(  6)=  -0.5807280000e+06
      a(  7)=  -0.1494530000e+05
      a(  8)=   0.7555000000D+10
      a(  9)=   0.5435040000e+06
      a( 10)=  -0.3938939000e+06
      a( 11)=   0.3562100000e+05
      a( 12)=   0.1934000000e+05
      a( 13)=  -0.1174580000e+05
      a( 14)=   0.1183500000e+05
      a( 15)=   0.5148500000e+00
      a( 16)=   0.1053000000e+02
      a( 17)=  -0.1794510000e+05
      a( 18)=   0.9666500000e+03
      a( 19)=   0.4839600000e+00
      a( 20)=   0.8115200000e+04
      a( 21)=  -0.4391950000e+04
      a( 22)=  -0.2596845000e+07
      a( 23)=   0.6719900000e+06
      a( 24)=  -0.3723690000e+06
      a( 25)=  -0.9280430000e+01
      a( 26)=   0.1058900000e+06
      a( 27)=  -0.1308200000e+02
      a( 28)=   0.2959500000e+03
      a( 29)=  -0.3241460000e+05
      a( 30)=  -0.1195130000e+03
      a( 31)=   0.1982800000e+00
      a( 32)=   0.5107000000e+00
      a( 33)=   0.5507900000e+00
      a( 34)=   0.1366000000e+01
      a( 35)=   0.1811600000e+05
      a( 36)=   0.1623500000e+06
      a( 37)=   0.3323000000e+00
      a( 38)=  -0.3803690000e+05
      a( 39)=   0.2234000000e+01
      a( 40)=   0.2316900000e+04
      a( 41)=  -0.1498590000e+04
      a( 42)=   0.1072600000e+06
      a( 43)=   0.1514200000e+01
      a( 44)=   0.9292900000e+00
      a( 45)=   0.0000000000e+00
      a( 46)=   0.0000000000e+00
      a( 47)=   0.0000000000e+00
      a( 48)=  -0.8305970000e+00
      a( 49)=   0.6727900000e+00
      a( 50)=  -0.1194400000e+00
      a( 51)=   0.1076200000e+01
      a( 52)=   0.1843400000e+01
      a( 53)=   0.2324900000e+00
      a( 54)=  -0.2225530000e+03
      a( 55)=  -0.2442710000e+04
      a( 56)=   0.1602200000e+01
      a( 57)=  -0.6928680000e+05
      a( 58)=   0.3168000000e+05
      a( 59)=   0.2960000000e+02
      a( 60)=   0.5040900000e+05
      a( 61)=   0.2104900000e+01
      a( 62)=   0.0000000000e+00
      a( 63)=   0.4401700000e+03
      a( 64)=  -0.2428640000e+05
      a( 65)=   0.0000000000e+00
      a( 66)=   0.0000000000e+00
      a( 67)=  -0.8560000000e+00
      a( 68)=   0.0000000000e+00
      a( 69)=   0.5626300000e+05
      a( 70)=   0.0000000000e+00
      rmix1 =  0.3000000000e+01
      rmix2 = -0.2000000000e+01
      rmix3 =  0.0000000000e+00
      rmix4 =  0.0000000000e+00
      dip0  =  0.1848868950e+01
      dip1  =  0.9347285625e+00
      alpa0 =  0.6145540000e+01
      alpa1 =  0.5394300000e+01
      alper0=  0.4788710000e+01
      alper1=  0.1389100000e+01
      g(1,1,0,1)=  0.2000000000e+01
      g(1,1,2,1)= -0.1000000000e+01
      g(1,2,1,1)=  0.1732050808e+01
      g(1,2,3,1)= -0.1154700538e+01
      g(2,1,1,1)=  0.1732050808e+01
      g(2,1,3,1)= -0.1154700538e+01
      g(1,3,2,1)=  0.1632993162e+01
      g(1,3,4,1)= -0.1224744871e+01
      g(2,2,0,1)=  0.2000000000e+01
      g(2,2,0,2)=  0.2000000000e+01
      g(2,2,2,1)=  0.1000000000e+01
      g(2,2,2,2)= -0.2000000000e+01
      g(2,2,4,1)= -0.1333333333e+01
      g(2,2,4,2)=  0.3333333333e+00
      g(3,1,2,1)=  0.1632993162e+01
      g(3,1,4,1)= -0.1224744871e+01
      g(1,4,3,1)=  0.1581138830e+01
      g(1,4,5,1)= -0.1264911064e+01
      g(2,3,1,1)=  0.1885618083e+01
      g(2,3,1,2)=  0.1490711985e+01
      g(2,3,3,1)=  0.7071067812e+00
      g(2,3,3,2)= -0.2236067977e+01
      g(2,3,5,1)= -0.1414213562e+01
      g(2,3,5,2)=  0.4472135955e+00
      g(3,2,1,1)=  0.1885618083e+01
      g(3,2,1,2)=  0.1490711985e+01
      g(3,2,3,1)=  0.7071067812e+00
      g(3,2,3,2)= -0.2236067977e+01
      g(3,2,5,1)= -0.1414213562e+01
      g(3,2,5,2)=  0.4472135955e+00
      g(4,1,3,1)=  0.1581138830e+01
      g(4,1,5,1)= -0.1264911064e+01
      g(2,4,2,1)=  0.1825741858e+01
      g(2,4,2,2)=  0.1290994449e+01
      g(2,4,4,1)=  0.5477225575e+00
      g(2,4,4,2)= -0.2323790008e+01
      g(2,4,6,1)= -0.1460593487e+01
      g(2,4,6,2)=  0.5163977795e+00
      g(3,3,0,1)=  0.2000000000e+01
      g(3,3,0,2)=  0.2000000000e+01
      g(3,3,0,3)=  0.2000000000e+01
      g(3,3,2,1)=  0.1500000000e+01
      g(3,3,2,2)=  0.0000000000e+00
      g(3,3,2,3)= -0.2500000000e+01
      g(3,3,4,1)=  0.3333333333e+00
      g(3,3,4,2)= -0.2333333333e+01
      g(3,3,4,3)=  0.1000000000e+01
      g(3,3,6,1)= -0.1500000000e+01
      g(3,3,6,2)=  0.6000000000e+00
      g(3,3,6,3)= -0.1000000000e+00
      g(4,2,2,1)=  0.1825741858e+01
      g(4,2,2,2)=  0.1290994449e+01
      g(4,2,4,1)=  0.5477225575e+00
      g(4,2,4,2)= -0.2323790008e+01
      g(4,2,6,1)= -0.1460593487e+01
      g(4,2,6,2)=  0.5163977795e+00
      pf(0,0)=  0.2820947918e+00
      pf(1,0)=  0.4886025119e+00
      pf(1,1)= -0.3454941495e+00
      pf(2,0)=  0.6307831305e+00
      pf(2,1)= -0.2575161347e+00
      pf(2,2)=  0.1287580673e+00
      pf(3,0)=  0.7463526652e+00
      pf(3,1)= -0.2154534561e+00
      pf(3,2)=  0.6813236510e-01
      pf(3,3)= -0.2781492158e-01
      pf(4,0)=  0.8462843753e+00
      pf(4,1)= -0.1892349392e+00
      pf(4,2)=  0.4460310290e-01
      pf(4,3)= -0.1192068068e-01
      pf(4,4)=  0.4214597071e-02
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C GPT normal
       r0=1.7327e0
       alpha=1.16937e0
       diss=49652.7e0
       gpt1=0.5847305e0
       gpt2=1.28e0**2
       gpt3=953320.3e0
       gpt4=164766.75e0
       gpt5=49656.989942e0
       scal=1.e0
       balo=9.5e0
       coup=50000.e0
       shi1=zero
       shi2=zero
      endif
      if(mod(ipotty,10).eq.1)then
       a(  1) =  0.1571983e04
       a( 19) =  0.145e0
       a( 23) =  0.047e0
       a( 26) =  -0.058e0
       a( 27) =  0.18e0
       scal=0.983e0
       balo=12.e0
       coup=zero
       shi1=0.01e0
       shi2=zero
      endif
      if(mod(ipotty,10).eq.2)then
       a(  1) =  0.156932e04
       a( 19) =  0.115e0
       a( 23) =  0.060e0
       a( 26) =  -0.045e0
       coup=zero
      endif
      if(mod(ipotty,10).lt.4)then
C
C 3j-parameters from general array of parameters
C
      icu=31
      do 211 j1=0,4,1
      do 212 j2=j1,4,1
      do 213 j3=(j2-j1),(j1+j2),2
      ac(1,j1,j2,j3)=a(icu)
      ac(2,j1,j2,j3)=a(icu+35)
      ac(3,j1,j2,j3)=a(icu+70)
      icu=icu+1
  213 continue
  212 continue
  211 continue
C
C prefactors for polynomials
C
      fac(0)=one
      do 111 i=1,8
      fac(i)=fac(i-1)*i
  111 continue
      do 112 l=0,4
      do 113 m=0,l
      pf(l,m)=(-1)**m*sqrt((2*l+1)*fac(l-m)/(4*pi*fac(l+m)))
  113 continue
  112 continue
C 
C 3j-symbols for potential energy surface
C (abc/f-f0)=g(a,b,c,f)
C Ulrich Schmitt, March 1987, June 1988
C Martin Suhm, May 1989, August 1989
C
      fct(1)=one
      efct(1)=zero
      do 114 i=1,30
      x=i
      j=i+1
      efct(j)=efct(i)
      fct(j)=fct(i)*x
    1 if(fct(j).ge.ten) then
        fct(j)=tenth*fct(j)
        efct(j)=efct(j)+one
        goto 1
      endif
  114 continue
C
C
      m2=0
      m3=0
      do 115 j1=0,8,2
      do 116 j2=0,8,2
      do 117 j3=iabs(j1-j2),iabs(j1+j2),4
      do 118 m1=0,min(j1,j2),2
      m2=-m1
      i1=max0(0,(j2-j3-m1)/2,(j1-j3+m2)/2)
      i2=min0((j1+j2-j3)/2,(j1-m1)/2,(j2+m2)/2)
      ia=i1
      ib=(j1+j2-j3)/2-i1
      ic=(j1-m1)/2-i1
      id=(j2+m2)/2-i1
      ie=(j3-j2+m1)/2+i1
      if=(j3-j1-m2)/2+i1
      sum=(-1)**i1/(fct(ia+1)*fct(ib+1)*fct(ic+1)*fct(id +1)*
     1 fct(ie +1)*fct(if+1))
      esum=-efct(ia+1)-efct(ib+1)-efct(ic+1)-efct(id +1)-
     1 efct(ie +1)-efct(if+1)
      nmax=i2-i1
      f=one
      if(nmax.gt.0) then
      do 119 i=nmax,1,-1
      x=(ib-i+1)*(ic-i+1)*(id - i+1)
      x=x/((ia+i)*(ie +i)*(if+i))
      f=one -x*f
  119 continue
      endif
      sum=sum*f
      f=fct((j1+j2-j3)/2+1)*fct((j1-j2+j3)/2+1)*
     1 fct((j2-j1+j3)/2+1)*fct((j1+m1)/2+1)*fct((j1-m1)/2+1)*
     2 fct((j2+m2)/2+1)*fct((j2-m2)/2+1)*fct((j3+m3)/2+1)*
     3 fct((j3-m3)/2+1)/fct((j1+j2+j3)/2+2)
      e=efct((j1+j2-j3)/2+1)+efct((j1-j2+j3)/2+1)+
     1 efct((j2-j1+j3)/2+1)+efct((j1+m1)/2+1)+efct((j1-m1)/2+1)+
     2 efct((j2+m2)/2+1)+efct((j2-m2)/2+1)+efct((j3+m3)/2+1)+
     3 efct((j3-m3)/2+1)-efct((j1+j2+j3)/2+2)
      cg3j=sum*sqrt(f)*ten**(esum+half*e)*(-1)**((j2-j1-m3)/2)
      g(j1/2,j2/2,j3/2,m1/2)=cg3j*pif
  118 continue
  117 continue
  116 continue
  115 continue
      endif
      return
      end
C
C *** subroutine hf2p *** computes the potential energy for n configurations
C
      subroutine hf2p(ipotty)
       implicit double precision (a-h,o-z)
      parameter(n=1,ma=135,nat=4)
C nat is the number of atoms (2*HF molecule number) for the total cluster
C this is needed in cofact to correct the one body term
      parameter (zero=0.0e0,one=1.0e0,two=2.0e0,three=3.0e0,oha=1.5e0)
      parameter (ten=10.0e0,tenth=0.1e0,half=0.5e0,third=1.e0/3.e0)
      parameter (pi=3.1415926535898e0,phicon=0.01745329252e0)
      dimension p1(n,0:4,0:4),p2(n,0:4,0:4),cs(n,4)
      dimension ar(n,0:26),ar4(n,0:4,0:4,0:8)
      dimension vp(n,0:26),vp4(n,0:4,0:4,0:8)
      dimension s(n),s1(n),s2(n),sq(n),s3(n)
      dimension c(n),c1(n),c2(n),cq(n),c3(n),c4(n),cl1(n),cl2(n)
      dimension c1p(n),c2p(n),ctp(n),s1p(n),s2p(n)
      dimension f014(n)
      dimension f023(n)
      dimension f0h(n)
      dimension f0f(n),f1f(n),f2f(n),f3f(n)
      dimension f4f(n),f5f(n),f6f(n),f7f(n)
      dimension rff3(n),rff6(n)
      dimension tdr1(n),tdr2(n)
      dimension qu1(n),qu2(n),adp1(n),adp2(n)
      dimension ap1(n),ap2(n),amu1(n),amu2(n)
      dimension dab1(n),dab2(n),dab3(n)
      dimension du(n),dv(n),dup(n),dvp(n),ds(n)
      dimension rr3(n),rr4(n),rr5(n),rr6(n),rr8(n),rr10(n)
      common/morcon/ r0,alpha,diss,scal,balo,coup,shi1,shi2
      common/energy/ e(n)
      common/intern/  r1(n), r2(n),  r(n),th1(n),th2(n),tau(n),
     1               r12(n),r14(n),r23(n),r34(n)
C kinter special for mod(ipotty,10)=4,5,6,7
      common/kinter/ ca(n),th1p(n),th2p(n),taup(n)
      common/helg/ dip0,dip1,alpa0,alpa1,alper0,alper1,
     *             rmix1,rmix2,rmix3,rmix4,
     *             gpt1,gpt2,gpt3,gpt4,gpt5
      common/ang/ g(0:4,0:4,0:8,0:4),pf(0:4,0:4)
      common/par/ a(ma),ac(3,0:4,0:4,0:8)
      pif=half/sqrt(pi)
CCCC
      if(mod(ipotty,10).gt.3)then
      if(mod(ipotty,10).eq.4.or.mod(ipotty,10).eq.5)then
      shdd=100.e0
      thirdo=third
      iso1=1
      iso2=0
      endif
      if(mod(ipotty,10).eq.6.or.mod(ipotty,10).eq.7)then
      shdd=200.e0
      thirdo=zero
      iso1=0
      iso2=1
      endif
C calculate legendre polynomials
      do 4100 id=1,n
      s1(id)=sin(th1(id))
      c1(id)=cos(th1(id))
      sq(id)=s1(id)*s1(id)
      s3(id)=sq(id)*s1(id)
      cq(id)=c1(id)*c1(id)
      c3(id)=cq(id)*c1(id)
      c4(id)=c3(id)*c1(id)
      p1(id,0,0)=pf(0,0)
      p1(id,1,0)=pf(1,0)*c1(id)
      p1(id,2,0)=pf(2,0)*0.5e0*(3*cq(id)-1)
      p1(id,3,0)=pf(3,0)*0.5e0*(5*c3(id)-3*c1(id))
      p1(id,4,0)=pf(4,0)*0.125e0*(35*c4(id)-30*cq(id)+3)
      p1(id,1,1)=pf(1,1)*s1(id)
      p1(id,2,1)=pf(2,1)*3*s1(id)*c1(id)
      p1(id,3,1)=pf(3,1)*1.5e0*s1(id)*(5*cq(id)-1)
      p1(id,4,1)=pf(4,1)*2.5e0*s1(id)*(7*c3(id)-3*c1(id))
      p1(id,2,2)=pf(2,2)*3*sq(id)
      p1(id,3,2)=pf(3,2)*15*sq(id)*c1(id)
      p1(id,4,2)=pf(4,2)*7.5e0*sq(id)*(7*cq(id)-1)
      p1(id,3,3)=pf(3,3)*15*s3(id)
      s2(id)=sin(th2(id))
      c2(id)=cos(th2(id))
      sq(id)=s2(id)*s2(id)
      s3(id)=sq(id)*s2(id)
      cq(id)=c2(id)*c2(id)
      c3(id)=cq(id)*c2(id)
      c4(id)=c3(id)*c2(id)
      p2(id,0,0)=pf(0,0)
      p2(id,1,0)=pf(1,0)*c2(id)
      p2(id,2,0)=pf(2,0)*0.5e0*(3*cq(id)-1)
      p2(id,3,0)=pf(3,0)*0.5e0*(5*c3(id)-3*c2(id))
      p2(id,4,0)=pf(4,0)*0.125e0*(35*c4(id)-30*cq(id)+3)
      p2(id,1,1)=pf(1,1)*s2(id)
      p2(id,2,1)=pf(2,1)*3*s2(id)*c2(id)
      p2(id,3,1)=pf(3,1)*1.5e0*s2(id)*(5*cq(id)-1)
      p2(id,4,1)=pf(4,1)*2.5e0*s2(id)*(7*c3(id)-3*c2(id))
      p2(id,2,2)=pf(2,2)*3*sq(id)
      p2(id,3,2)=pf(3,2)*15*sq(id)*c2(id)
      p2(id,4,2)=pf(4,2)*7.5e0*sq(id)*(7*cq(id)-1)
      p2(id,3,3)=pf(3,3)*15*s3(id)
      do 4110 j=1,3
      cs(id,j)=cos(j*tau(id))
 4110 continue
      c1p(id)=cos(th1p(id))
      c2p(id)=cos(th2p(id))
      ctp(id)=cos(taup(id))
      s1p(id)=sin(th1p(id))
      s2p(id)=sin(th2p(id))
      f014(id)=exp(-oha*r14(id))
      f023(id)=exp(-oha*r23(id))
      f0h(id)=exp(-oha*r12(id))
      f0f(id)=exp(-oha*r34(id))
      f1f(id)=f0f(id)*r34(id)
      f2f(id)=f1f(id)*r34(id)
      f3f(id)=f2f(id)*r34(id)
      f4f(id)=f3f(id)*r34(id)
      f5f(id)=f4f(id)*r34(id)
      f6f(id)=f5f(id)*r34(id)
      f7f(id)=f6f(id)*r34(id)
      rff3(id)=r34(id)**3
      rff6(id)=rff3(id)**2
      tdr1(id)=tanh(r1(id)-1.7327e0)
      tdr2(id)=tanh(r2(id)-1.7327e0)
c      amu1=1.045e0*r1(id)/((r1(id)/6.e0)**4+1)
c      if(r1(id).gt.2.e0) amu1=amu1*exp(-((r1(id)-2.e0)/2.e0)**2)
c      amu2=1.045e0*r2(id)/((r2(id)/6.e0)**4+1)
c      if(r2(id).gt.2.e0) amu2=amu2*exp(-((r2(id)-2.e0)/2.e0)**2)
C dipoles
      amu1(id)=dip0+iso2*dip1*tdr1(id)+iso1*dip0*
     1        (dip1*rff3(id)/(300+rff3(id)))*tdr1(id)
      amu2(id)=dip0+iso2*dip1*tdr2(id)+iso1*dip0*
     1        (dip1*rff3(id)/(300+rff3(id)))*tdr2(id)
c quadrupoles
      qu1(id)=2.35e0*(1+a(15)*tdr1(id))
      qu2(id)=2.35e0*(1+a(15)*tdr2(id))
c perpendicular polarizability
      ap1(id)=alper0+tdr1(id)*alper1
      ap2(id)=alper0+tdr2(id)*alper1
c parallel-perpendicular polarizability
      adp1(id)=(alpa0+tdr1(id)*alpa1-ap1(id))
      adp2(id)=(alpa0+tdr2(id)*alpa1-ap2(id))
      y=0.e0
c HH FF FH terms:
c --------------
      y=
     *  a(1)*f0h(id)*(1+a(50)*(tdr1(id)+tdr2(id)))
     *   +a(2)*f0h(id)*r12(id)**third
     *   +f0h(id)*(a(21)+a(4)*r12(id))*
     *   (((3*c2p(id)**2-1)*c1p(id)-2*s2p(id)*c2p(id)*s1p(id)*ctp(id))
     *   -((3*c1p(id)**2-1)*c2p(id)-2*s1p(id)*c1p(id)*s2p(id)*ctp(id)))
      y=y+(a(13)*f0h(id)*r12(id)**2+a(14)*f3f(id))*(c1p(id)*c2p(id)
     *   -half*s1p(id)*s2p(id)*ctp(id))
      y=y+a(5)*f0f(id)+a(8)*f0f(id)**3+a(6)*f1f(id)
      y=y+(f014(id)*(a(9)+r14(id)*(a(10)+a(11)*r14(id))))
     *  *((1+a(48)*tdr1(id)+a(53)*tdr2(id)+a(49)*tdr1(id)**2))
      y=y+(f023(id)*(a(9)+r23(id)*(a(10)+a(11)*r23(id))))
     *  *((1+a(48)*tdr2(id)+a(53)*tdr1(id)+a(49)*tdr2(id)**2))
c physical terms:
c --------------
      y=y+
c Dip-Dip:
     *   (-6.79437e04*amu1(id)*amu2(id)/(rff3(id)+shdd))
     *   *(c1(id)*c2(id)-half*s1(id)*s2(id)*cs(id,1))
c Dip-Q:
     *   +9.62963e04/(r34(id)**4+850-125*r34(id))*
     *    (amu1(id)*qu2(id)*( 
     *    (3*c2(id)**2-1)*c1(id)-2*s2(id)*c2(id)*s1(id)*cs(id,1))+
     *     qu1(id)*amu2(id)*(
     *   -(3*c1(id)**2-1)*c2(id)+2*s1(id)*c1(id)*s2(id)*cs(id,1)))
c Q-Q:
     *   +9.09868e04*qu1(id)*qu2(id)/(r34(id)**5+550)*
     *    (1-5*c1(id)**2-5*c2(id)**2-15*c1(id)**2*c2(id)**2+
     *     2*(4*c1(id)*c2(id)-s1(id)*s2(id)*cs(id,1))**2)
c polaris isotropic:
     *   -half*33972.e0/(rff6(id)+10**5)*
     *    (amu1(id)**2*ap2(id)*(3*c1(id)**2+1)
     *    +amu2(id)**2*ap1(id)*(3*c2(id)**2+1))
c polaris anisotropic:
     *   -half*33972.e0/(rff6(id)+10**5)*
     *    (amu1(id)**2*adp2(id)+amu2(id)**2*adp1(id))
     *    *(ca(id)+3*c2(id)*c1(id))**2
c dispersion:
     *   -5e06/(rff6(id)+10**5)*(1+(tdr1(id)+tdr2(id)))
      e(id)=y
c angular radial expansion terms:
c ------------------------------
      ar4(id,0,0,0)=zero
      ar4(id,2,2,0)=zero
      ar4(id,3,3,2)=zero
      ar4(id,3,3,6)=zero
C02
      ar4(id,0,1,1)= (a(22)*f1f(id))*(1+a(31)*tdr2(id))
     *          +(a(60)*f3f(id))*(1+a(37)*tdr2(id))
      ar4(id,1,0,1)=-(a(22)*f1f(id))*(1+a(31)*tdr1(id))
     *          -(a(60)*f3f(id))*(1+a(37)*tdr1(id))
C03
      ar4(id,0,2,2)= a(23)*f1f(id)*r34(id)**thirdo
     *          +a(42)*f0h(id)
      ar4(id,2,0,2)= a(23)*f1f(id)*r34(id)**thirdo
     *          +a(42)*f0h(id)
C04
      ar4(id,0,3,3)= a(24)*f1f(id)
      ar4(id,3,0,3)=-a(24)*f1f(id)
C05
      ar4(id,0,4,4)= a(35)*f2f(id)*r34(id)**thirdo
      ar4(id,4,0,4)= a(35)*f2f(id)*r34(id)**thirdo
C06
      ar4(id,1,1,0)= a(36)*f0h(id)+a(34)*f7f(id)+
     *           a(26)*f1f(id)*(1+a(44)*(tdr1(id)+tdr2(id)))
C07
      ar4(id,1,1,2)= a(55)*f5f(id)
C08
      ar4(id,1,2,1)= (a(38)*f2f(id)+a(46)*f4f(id))*
     *             (1+a(51)*tdr2(id)+a(19)*tdr1(id))
      ar4(id,2,1,1)=-(a(38)*f2f(id)+a(46)*f4f(id))*
     *             (1+a(51)*tdr1(id)+a(19)*tdr2(id))
C09
      ar4(id,1,2,3)= (a(54)*f5f(id))*(1+a(67)*tdr2(id)+a(61)*tdr1(id))
      ar4(id,2,1,3)=-(a(54)*f5f(id))*(1+a(67)*tdr1(id)+a(61)*tdr2(id))
C10
      ar4(id,1,3,2)= (a(40)*f3f(id))*(1+a(56)*tdr2(id))
      ar4(id,3,1,2)= (a(40)*f3f(id))*(1+a(56)*tdr1(id))
C11
      ar4(id,1,3,4)= (a(25)*f7f(id)+a(41)*f3f(id))
     *          *(1+a(43)*tdr2(id))
      ar4(id,3,1,4)= (a(25)*f7f(id)+a(41)*f3f(id))
     *          *(1+a(43)*tdr1(id))
C12
      ar4(id,1,4,3)= (a(7)*f1f(id))*(1+a(45)*tdr2(id))
      ar4(id,4,1,3)=-(a(7)*f1f(id))*(1+a(45)*tdr1(id))
C13
      ar4(id,1,4,5)= (a(28)*f4f(id))
      ar4(id,4,1,5)=-(a(28)*f4f(id))
C15
      ar4(id,2,2,2)= a(17)*f2f(id)*(1+a(32)*(tdr1(id)+tdr2(id)))
C16
      ar4(id,2,2,4)= a(63)*f5f(id)+a(64)*f3f(id)
C17
      ar4(id,2,3,1)= (a(57)*f0f(id))*(1+a(39)*tdr2(id))
      ar4(id,3,2,1)=-(a(57)*f0f(id))*(1+a(39)*tdr1(id))
C18
      ar4(id,2,3,3)= a(18)*f3f(id)
      ar4(id,3,2,3)=-a(18)*f3f(id)
C19
      ar4(id,2,3,5)= (a(20)*f4f(id)+a(29)*f3f(id)+a(30)*f6f(id))
     *           *(1+a(52)*tdr2(id)+a(33)*tdr1(id))
      ar4(id,3,2,5)=-(a(20)*f4f(id)+a(29)*f3f(id)+a(30)*f6f(id))
     *           *(1+a(52)*tdr1(id)+a(33)*tdr2(id))
C20
      ar4(id,2,4,2)= (a(66)*f3f(id)+a(68)*f4f(id))*(1+a(47)*tdr2(id))
      ar4(id,4,2,2)= (a(66)*f3f(id)+a(68)*f4f(id))*(1+a(47)*tdr1(id))
C21
      ar4(id,2,4,4)= a(27)*(iso1*f7f(id)+iso2*f5f(id))
      ar4(id,4,2,4)= a(27)*(iso1*f7f(id)+iso2*f5f(id))
C22
      ar4(id,2,4,6)= a(69)*f1f(id)+a(3)*f3f(id)+a(59)*
     *             f5f(id)*(1+a(16)*tdr2(id))
      ar4(id,4,2,6)= a(69)*f1f(id)+a(3)*f3f(id)+a(59)*
     *             f5f(id)*(1+a(16)*tdr1(id))
C23
      ar4(id,3,3,0)= a(12)*f0f(id)
C25
      ar4(id,3,3,4)= a(58)*(iso1*f5f(id)+iso2*f0f(id))
 4100 continue
      do 4102 j1j2=0,6
      do 4103 j1=max(0,j1j2-4),min(j1j2,4)
      j2=j1j2-j1
      do 4104 j3=abs(j2-j1),(j2+j1),2
      do 4101 i=1,n
      vp4(i,j1,j2,j3)=p1(i,j1,0)*p2(i,j2,0)
      do 4105 m=1,min(j1,j2),1
      vp4(i,j1,j2,j3)=vp4(i,j1,j2,j3)+g(j1,j2,j3,m)*
     *            cs(i,m)*p1(i,j1,m)*p2(i,j2,m)
 4105 continue
      e(i)=e(i)+vp4(i,j1,j2,j3)*ar4(i,j1,j2,j3)
 4101 continue
 4104 continue
 4103 continue
 4102 continue
      cofact=two/(nat-two)
      if(mod(ipotty,10).eq.4.or.mod(ipotty,10).eq.6)then
      do 4119 k=1,n
      e(k)=e(k)
     *         +cofact*diss*(one -exp(-alpha*(r1(k)-r0)))**2
     *         +cofact*diss*(one -exp(-alpha*(r2(k)-r0)))**2
 4119 continue
      elseif(mod(ipotty,10).eq.5.or.mod(ipotty,10).eq.7)then
      do 4118 k=1,n
      expa=exp(-2*gpt1*r1(k))
      bexpa=expa*gpt2
      e(k)=e(k)+
     * cofact*(expa*(-gpt3/(1+bexpa)**2+gpt4/(1-bexpa)**2)+gpt5)
      expa=exp(-2*gpt1*r2(k))
      bexpa=expa*gpt2
      e(k)=e(k)+
     * cofact*(expa*(-gpt3/(1+bexpa)**2+gpt4/(1-bexpa)**2)+gpt5)
 4118 continue
      endif

      else
CCCC
      do 210 j=1,2
      do 310 k=1,n
      cs(k,j)=cos(j*tau(k))
  310 continue
  210 continue
      do 311 k=1,n
      s(k)=sin(th1(k))
      c(k)=cos(th1(k))
      s2(k)=s(k)*s(k)
      c2(k)=c(k)*c(k)
      cl1(k)=c2(k)
      c3(k)=c2(k)*c(k)
      c4(k)=c3(k)*c(k)
      p1(k,0,0)=pf(0,0)
      p1(k,1,0)=pf(1,0)*c(k)
      p1(k,2,0)=pf(2,0)*0.5e0*(3*c2(k)-1.e0)
      p1(k,3,0)=pf(3,0)*0.5e0*(5*c3(k)-3*c(k))
      p1(k,4,0)=pf(4,0)*0.125e0*(35*c4(k)-30*c2(k)+3.e0)
      p1(k,1,1)=pf(1,1)*s(k)
      p1(k,2,1)=pf(2,1)*3*s(k)*c(k)
      p1(k,3,1)=pf(3,1)*1.5e0*s(k)*(5*c2(k)-1)
      p1(k,4,1)=pf(4,1)*2.5e0*s(k)*(7*c3(k)-3*c(k))
      p1(k,2,2)=pf(2,2)*3*s2(k)
      p1(k,3,2)=pf(3,2)*15*s2(k)*c(k)
      p1(k,4,2)=pf(4,2)*7.5e0*s2(k)*(7*c2(k)-1)
  311 continue
      do 312 k=1,n
      s(k)=sin(th2(k))
      c(k)=cos(th2(k))
      s2(k)=s(k)*s(k)
      c2(k)=c(k)*c(k)
      cl2(k)=c2(k)
      c3(k)=c2(k)*c(k)
      c4(k)=c3(k)*c(k)
      p2(k,0,0)=pf(0,0)
      p2(k,1,0)=pf(1,0)*c(k)
      p2(k,2,0)=pf(2,0)*0.5e0*(3*c2(k)-1.e0)
      p2(k,3,0)=pf(3,0)*0.5e0*(5*c3(k)-3*c(k))
      p2(k,4,0)=pf(4,0)*0.125e0*(35*c4(k)-30*c2(k)+3.e0)
      p2(k,1,1)=pf(1,1)*s(k)
      p2(k,2,1)=pf(2,1)*3*s(k)*c(k)
      p2(k,3,1)=pf(3,1)*1.5e0*s(k)*(5*c2(k)-1)
      p2(k,4,1)=pf(4,1)*2.5e0*s(k)*(7*c3(k)-3*c(k))
      p2(k,2,2)=pf(2,2)*3*s2(k)
      p2(k,3,2)=pf(3,2)*15*s2(k)*c(k)
      p2(k,4,2)=pf(4,2)*7.5e0*s2(k)*(7*c2(k)-1)
  312 continue
      do 410 k=1,n
      dab1(k)=exp(-0.20e0*(r(k)-3.e0)**2)
      dab2(k)=exp(-0.35e0*(r(k)-4.e0)**2)
      dab3(k)=exp(-0.80e0*(r(k)-4.9e0)**2)
      rr5(k)=one/r(k)**5
      rr4(k)=rr5(k)*r(k)
      rr3(k)=rr4(k)*r(k)
      rr6(k)=rr3(k)*rr3(k)
      rr8(k)=rr4(k)*rr4(k)
      rr10(k)=rr5(k)*rr5(k)
  410 continue
      do 411 k=1,n
      ar(k,0)=ac(1,0,0,0)*dab1(k)+ac(3,0,0,0)*dab3(k)
      ar(k,1)=ac(1,1,1,0)*dab1(k)+ac(3,1,1,0)*dab3(k)
      ar(k,2)=ac(1,1,2,3)*dab1(k)
      ar(k,4)=ac(1,2,2,4)*dab1(k)
      ar(k,5)=ac(1,2,4,6)*dab1(k)
      ar(k,7)=ac(2,0,2,2)*dab2(k)
      ar(k,9)=ac(2,1,2,1)*dab2(k)
      ar(k,11)=ac(2,2,3,5)*dab2(k)
      ar(k,13)=ac(3,0,1,1)*dab3(k)
      ar(k,15)=ac(3,0,3,3)*dab3(k)
      ar(k,17)=ac(3,0,4,4)*dab3(k)
      ar(k,19)=ac(3,1,1,2)*dab3(k)
      ar(k,24)=ac(3,1,3,2)*dab3(k)
      ar(k,20)=ac(3,1,4,5)*dab3(k)
      ar(k,26)=ac(3,2,2,2)*dab3(k)
  411 continue
      do 412 k=1,n
      ar(k,3)=-ar(k,2)
      ar(k,6)= ar(k,5)
      ar(k,8)= ar(k,7)
      ar(k,10)=-ar(k,9)
      ar(k,12)=-ar(k,11)
      ar(k,14)=-ar(k,13)
      ar(k,16)=-ar(k,15)
      ar(k,18)= ar(k,17)
      ar(k,25)= ar(k,24)
      ar(k,21)=-ar(k,20)
  412 continue
      do 910 k=1,n
      vp(k,0)    =(g(0,0,0,0)*p1(k,0,0)*p2(k,0,0)
     1           )*(+1)
      vp(k,1)    =(g(1,1,0,0)*p1(k,1,0)*p2(k,1,0)
     1          -2*g(1,1,0,1)*p1(k,1,1)*p2(k,1,1)*cs(k,1)
     2           )*(+1)
      vp(k,2)    =(g(1,2,3,0)*p1(k,1,0)*p2(k,2,0)
     1          -2*g(1,2,3,1)*p1(k,1,1)*p2(k,2,1)*cs(k,1)
     2           )*(-7)
      vp(k,3)    =(g(2,1,3,0)*p1(k,2,0)*p2(k,1,0)
     1          -2*g(2,1,3,1)*p1(k,2,1)*p2(k,1,1)*cs(k,1)
     2           )*(-7)
      vp(k,4)    =(g(2,2,4,0)*p1(k,2,0)*p2(k,2,0)
     1          -2*g(2,2,4,1)*p1(k,2,1)*p2(k,2,1)*cs(k,1)
     2          +2*g(2,2,4,2)*p1(k,2,2)*p2(k,2,2)*cs(k,2)
     3           )*(+9)
      vp(k,5)    =(g(2,4,6,0)*p1(k,2,0)*p2(k,4,0)
     1          -2*g(2,4,6,1)*p1(k,2,1)*p2(k,4,1)*cs(k,1)
     2          +2*g(2,4,6,2)*p1(k,2,2)*p2(k,4,2)*cs(k,2)
     3           )*(+13)
      vp(k,6)    =(g(4,2,6,0)*p1(k,4,0)*p2(k,2,0)
     1          -2*g(4,2,6,1)*p1(k,4,1)*p2(k,2,1)*cs(k,1)
     2          +2*g(4,2,6,2)*p1(k,4,2)*p2(k,2,2)*cs(k,2)
     3           )*(+13)
      vp(k,7)    =(g(0,2,2,0)*p1(k,0,0)*p2(k,2,0)
     1           )*(+5)
      vp(k,8)    =(g(2,0,2,0)*p1(k,2,0)*p2(k,0,0)
     1           )*(+5)
      vp(k,9)    =(g(1,2,1,0)*p1(k,1,0)*p2(k,2,0)
     1          -2*g(1,2,1,1)*p1(k,1,1)*p2(k,2,1)*cs(k,1)
     2           )*(-3)
      vp(k,10)   =(g(2,1,1,0)*p1(k,2,0)*p2(k,1,0)
     1          -2*g(2,1,1,1)*p1(k,2,1)*p2(k,1,1)*cs(k,1)
     2           )*(-3)
  910 continue
      do 911 k=1,n
      vp(k,11)  =(g(2,3,5,0)*p1(k,2,0)*p2(k,3,0)
     1          -2*g(2,3,5,1)*p1(k,2,1)*p2(k,3,1)*cs(k,1)
     2          +2*g(2,3,5,2)*p1(k,2,2)*p2(k,3,2)*cs(k,2)
     3           )*(-11)
      vp(k,12)  =(g(3,2,5,0)*p1(k,3,0)*p2(k,2,0)
     1          -2*g(3,2,5,1)*p1(k,3,1)*p2(k,2,1)*cs(k,1)
     2          +2*g(3,2,5,2)*p1(k,3,2)*p2(k,2,2)*cs(k,2)
     3           )*(-11)
      vp(k,13)   =(g(0,1,1,0)*p1(k,0,0)*p2(k,1,0)
     1           )*(-3)
      vp(k,14)   =(g(1,0,1,0)*p1(k,1,0)*p2(k,0,0)
     1           )*(-3)
      vp(k,15)   =(g(0,3,3,0)*p1(k,0,0)*p2(k,3,0)
     1           )*(-7)
      vp(k,16)   =(g(3,0,3,0)*p1(k,3,0)*p2(k,0,0)
     1           )*(-7)
      vp(k,17)   =(g(0,4,4,0)*p1(k,0,0)*p2(k,4,0)
     1           )*(+9)
      vp(k,18)   =(g(4,0,4,0)*p1(k,4,0)*p2(k,0,0)
     1           )*(+9)
      vp(k,19)   =(g(1,1,2,0)*p1(k,1,0)*p2(k,1,0)
     1          -2*g(1,1,2,1)*p1(k,1,1)*p2(k,1,1)*cs(k,1)
     2           )*(+5)
      vp(k,24)   =(g(1,3,2,0)*p1(k,1,0)*p2(k,3,0)
     1          -2*g(1,3,2,1)*p1(k,1,1)*p2(k,3,1)*cs(k,1)
     2           )*(+5)
      vp(k,25)   =(g(3,1,2,0)*p1(k,3,0)*p2(k,1,0)
     1          -2*g(3,1,2,1)*p1(k,3,1)*p2(k,1,1)*cs(k,1)
     2           )*(+5)
      vp(k,22)   =(g(1,3,4,0)*p1(k,1,0)*p2(k,3,0)
     1          -2*g(1,3,4,1)*p1(k,1,1)*p2(k,3,1)*cs(k,1)
     2           )*(+9)
      vp(k,23)   =(g(3,1,4,0)*p1(k,3,0)*p2(k,1,0)
     1          -2*g(3,1,4,1)*p1(k,3,1)*p2(k,1,1)*cs(k,1)
     2           )*(+9)
      vp(k,20)   =(g(1,4,5,0)*p1(k,1,0)*p2(k,4,0)
     1          -2*g(1,4,5,1)*p1(k,1,1)*p2(k,4,1)*cs(k,1)
     2           )*(-11)
      vp(k,21)   =(g(4,1,5,0)*p1(k,4,0)*p2(k,1,0)
     1          -2*g(4,1,5,1)*p1(k,4,1)*p2(k,1,1)*cs(k,1)
     2           )*(-11)
      vp(k,26)   =(g(2,2,2,0)*p1(k,2,0)*p2(k,2,0)
     1          -2*g(2,2,2,1)*p1(k,2,1)*p2(k,2,1)*cs(k,1)
     2          +2*g(2,2,2,2)*p1(k,2,2)*p2(k,2,2)*cs(k,2)
     3           )*(+5)
  911 continue
C
      do 321 k=1,n
      e(k)=0.0e0
  321 continue
      do 522 k=1,n
      e(k)=e(k)+vp(k,0)*ar(k,0)+vp(k,1)*ar(k,1)
     1         +vp(k,2)*ar(k,2)+vp(k,3)*ar(k,3)
     2         +vp(k,4)*ar(k,4)+vp(k,5)*ar(k,5)
     3         +vp(k,6)*ar(k,6)+vp(k,7)*ar(k,7)
     4         +vp(k,8)*ar(k,8)+vp(k,9)*ar(k,9)
     5         +vp(k,10)*ar(k,10)+vp(k,11)*ar(k,11)
     6         +vp(k,12)*ar(k,12)+vp(k,13)*ar(k,13)
     7         +vp(k,14)*ar(k,14)+vp(k,15)*ar(k,15)
     8         +vp(k,16)*ar(k,16)+vp(k,17)*ar(k,17)
     9         +vp(k,18)*ar(k,18)+vp(k,19)*ar(k,19)
     1         +vp(k,24)*ar(k,24)+vp(k,25)*ar(k,25)
     2         +vp(k,20)*ar(k,20)+vp(k,21)*ar(k,21)
     3         +vp(k,26)*ar(k,26)
  522 continue
C
      do 322 k=1,n
      e(k)=e(k)+a(1)
  322 continue
      do 323 k=1,n
      if(r(k).lt.7.2e0)then
      e(k)=e(k)+vp(k,0)*exp(-((7.2e0-r(k))/r(k))**2)*a(2)*rr6(k)
      else
      e(k)=e(k)+vp(k,0)*a(2)*rr6(k)
      endif
      if(r(k).lt.8.4e0)then
      e(k)=e(k)+vp(k,0)*exp(-((8.4e0-r(k))/r(k))**2)*a(3)*rr8(k)
      else
      e(k)=e(k)+vp(k,0)*a(3)*rr8(k)
      endif
      if(r(k).lt.9.6e0)then
      e(k)=e(k)+vp(k,0)*exp(-((9.6e0-r(k))/r(k))**2)*a(4)*rr10(k)
      else
      e(k)=e(k)+vp(k,0)*a(4)*rr10(k)
      endif
  323 continue
      do 324 k=1,n
      if(r(k).lt.9.5e0) then
      e(k)=e(k)+vp(k,19)*a(5)*rr3(k)*exp(-((9.5e0-r(k))/r(k))**2)
      else
      e(k)=e(k)+vp(k,19)*a(5)*rr3(k)
      endif
      if(r(k).lt.10.5e0) then
      e(k)=e(k)+(vp(k,3)-vp(k,2))*a(6)*rr4(k)*
     1     exp(-((10.5e0-r(k))/r(k))**2)
      else
      e(k)=e(k)+(vp(k,3)-vp(k,2))*a(6)*rr4(k)
      endif
      if(r(k).lt.11.5e0) then
      e(k)=e(k)+vp(k,4)*a(7)*rr5(k)*exp(-((11.5e0-r(k))/r(k))**2)
      e(k)=e(k)+(vp(k,23)+vp(k,22))*a(8)*rr5(k)*
     1     exp(-((11.5e0-r(k))/r(k))**2)
      else
      e(k)=e(k)+vp(k,4)*a(7)*rr5(k)
      e(k)=e(k)+(vp(k,23)+vp(k,22))*a(8)*rr5(k)
C was vp(k,23)+vp(k,23) until 1.2.95, only long range, negligible
      endif
  324 continue
      do 325 k=1,n
      e(k)=e(k)+a(9)*exp(-0.5e0*r12(k))
      e(k)=e(k)+a(10)*exp(-1.5e0*r12(k))
      e(k)=e(k)+a(13)*exp(-1.5e0*r34(k))
      e(k)=e(k)+a(17)*(exp(-10*r14(k))+exp(-10*r23(k)))
  325 continue
      do 326 k=1,n
      if(r(k).lt.4.5e0)then
      e(k)=e(k)+a(18)*rr6(k)*((3*(cl1(k))+1)+
     1               (3*(cl2(k))+1))*
     2                exp(-((4.5e0-r(k))/r(k))**2)
      else
      e(k)=e(k)+a(18)*rr6(k)*((3*(cl1(k))+1)+
     1               (3*(cl2(k))+1))
      endif
  326 continue
      do 337 k=1,n
      e(k)=e(k)+a(22)*exp(-4*r14(k))*exp(-4*r23(k))
      e(k)=e(k)+a(28)*(exp(-4*(r14(k)-r1(k))**2)
     1                +exp(-4*(r23(k)-r2(k))**2))
  337 continue
      if(mod(ipotty,10).gt.1)then
      do 338 k=1,n
      e(k)=e(k)+a(24)*(exp(-r14(k))+exp(-r23(k)))
     1               *(exp(-15*r0)+exp(-15*r0))
  338 continue
      else
      do 339 k=1,n
      e(k)=e(k)+a(24)*(exp(-r14(k))+exp(-r23(k)))
     1               *(exp(-15*r1(k))+exp(-15*r2(k)))
  339 continue
      endif
      if(mod(ipotty,10).gt.0)then
C preparation of damping term,
C barrier lowering, coupling and scaling of the intermolecular part
      do 328 k=1,n
      ds(k)=(one -exp(-half*(r23(k)-r14(k))**2))
      e(k)=e(k)
     1    -balo*r(k)**2*exp(-.04e0*r(k)**2-half*(r23(k)-r14(k))**2)
      e(k)=e(k)-coup*(exp(-r14(k))+exp(-r23(k)))*
     1     (r1(k)-r0)*(r2(k)-r0)/(one +(r1(k)-r0)**2+(r2(k)-r0)**2)
      e(k)=e(k)*scal
  328 continue
      endif

C
C monomer terms
C
C the following version avoids tanh-calls for vectorization reasons
C The identity used is 1-tanh x= 2/ (exp 2x + 1), if the argument 
C exceeds 30., it is replaced by 30. to avoid overflow
C
      if(mod(ipotty,10).gt.0)then
      do 327 k=1,n
      du(k)=min(a(27)*(r14(k)-three),30.e0)
      dup(k)=exp(2*du(k))+one
      du(k)=two/dup(k)
      dv(k)=min(a(21)*(r14(k)-three),30.e0)
      dvp(k)=exp(2*dv(k))+one
      dv(k)=two/dvp(k)
      e(k)=e(k)+diss*(one -a(19)*ds(k)*
     1         exp(-a(23)*(r14(k)-three)**2))*
     2         (one -exp(-alpha*(one -a(26)*
     3         ds(k)*
     4         du(k))*
     5         (r1(k)-r0*(one -a(20)*dv(k)))))**2
      du(k)=min(a(27)*(r23(k)-three),30.e0)
      dup(k)=exp(2*du(k))+one
      du(k)=two/dup(k)
      dv(k)=min(a(21)*(r23(k)-three),30.e0)
      dvp(k)=exp(2*dv(k))+one
      dv(k)=two/dvp(k)
      e(k)=e(k)+diss*(one -a(19)*ds(k)*
     1         exp(-a(23)*(r23(k)-three)**2))*
     2         (one -exp(-alpha*(one -a(26)*
     3         ds(k)*
     4         du(k))*
     5         (r2(k)-r0*(one -a(20)*dv(k)))))**2
  327 continue
C proper one - body implementation
      cofact=two/(real(nat)-two)
      do 819 k=1,n
      e(k)=e(k)
     *         -(one -cofact)*diss*(one -exp(-alpha*(r1(k)-r0)))**2
     *         -(one -cofact)*diss*(one -exp(-alpha*(r2(k)-r0)))**2
  819 continue
      else
      do 329 k=1,n
      du(k)=min(a(27)*(r14(k)-three),30.e0)
      dup(k)=exp(2*du(k))+one
      du(k)=two/dup(k)
      dv(k)=min(a(21)*(r14(k)-three),30.e0)
      dvp(k)=exp(2*dv(k))+one
      dv(k)=two/dvp(k)
      e(k)=e(k)+diss*(one -a(19)*exp(-a(23)*(r23(k)-three)**2))*
     1         (one -exp(-alpha*(one -a(26)*du(k))*
     2         (r1(k)-r0*(one -a(20)*dv(k)))))**2 
      du(k)=min(a(27)*(r23(k)-three),30.e0)
      dup(k)=exp(2*du(k))+one
      du(k)=two/dup(k)
      dv(k)=min(a(21)*(r23(k)-three),30.e0)
      dvp(k)=exp(2*dv(k))+one
      dv(k)=two/dvp(k)
      e(k)=e(k)+diss*(one -a(19)*exp(-a(23)*(r14(k)-three)**2))*
     1         (one -exp(-alpha*(one -a(26)*du(k))*
     2         (r2(k)-r0*(one -a(20)*dv(k)))))**2  
  329 continue
C proper one -body implementation
      cofact=two/(real(nat)-two)
      do 818 k=1,n
      e(k)=e(k)
     *         -(one -cofact)*diss*(one -exp(-alpha*(r1(k)-r0)))**2
     *         -(one -cofact)*diss*(one -exp(-alpha*(r2(k)-r0)))**2
  818 continue
      endif
      endif
      return
      end
