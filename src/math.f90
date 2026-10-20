    Module  math_subroutines
    use types
    implicit none
    real(rk),PRIVATE::tol=1.0e-14_rk



  Contains

              
!subroutines for solving the linear systems 
     SUBROUTINE ludcmp(a,n,np,indx,d)  
      integer(ik) n,np,indx(n),NMAX  
      REAL (rk):: d,a(np,np),TINY  
      PARAMETER (NMAX=500,TINY=1.0e-20_rk)  
      integer(ik) i,imax,j,k  
      REAL (rk):: aamax,dum,sum,vv(NMAX)  
      d=1.  
      do 12 i=1,n  
        aamax=0._rk  
        do 11 j=1,n  
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))  
11      continue  
        if (aamax.eq.0._rk) stop 'singular matrix in ludcmp'  
        vv(i)=1._rk/aamax  
12    continue  
      do 19 j=1,n  
        do 14 i=1,j-1  
          sum=a(i,j)  
          do 13 k=1,i-1  
            sum=sum-a(i,k)*a(k,j)  
13        continue  
          a(i,j)=sum  
14      continue  
        aamax=0._rk  
        do 16 i=j,n  
          sum=a(i,j)  
          do 15 k=1,j-1  
            sum=sum-a(i,k)*a(k,j)  
15        continue  
          a(i,j)=sum  
          dum=vv(i)*abs(sum)  
          if (dum.ge.aamax) then  
            imax=i  
            aamax=dum  
          endif  
16      continue  
        if (j.ne.imax)then  
          do 17 k=1,n  
            dum=a(imax,k)  
            a(imax,k)=a(j,k)  
            a(j,k)=dum  
17        continue  
          d=-d  
          vv(imax)=vv(j)  
        endif  
        indx(j)=imax  
        if(a(j,j).eq.0._rk)a(j,j)=TINY  
        if(j.ne.n)then  
          dum=1._rk/a(j,j)  
          do 18 i=j+1,n  
            a(i,j)=a(i,j)*dum  
18        continue  
        endif  
19    continue  
      return  
      END subroutine ludcmp  

      SUBROUTINE lubksb(a,n,np,indx,b)  
      integer(ik) n,np,indx(n)  
      REAL (rk):: a(np,np),b(n)  
      integer(ik) i,ii,j,ll  
      REAL (rk):: sum  
      ii=0  
      do 12 i=1,n  
        ll=indx(i)  
        sum=b(ll)  
        b(ll)=b(i)  
        if (ii.ne.0)then  
          do 11 j=ii,i-1  
            sum=sum-a(i,j)*b(j)  
11        continue  
        else if (sum.ne.0._rk) then  
          ii=i  
        endif  
        b(i)=sum  
12    continue  
      do 14 i=n,1,-1  
        sum=b(i)  
        do 13 j=i+1,n  
          sum=sum-a(i,j)*b(j)  
13      continue  
        b(i)=sum/a(i,i)  
14    continue  
      return  
      END subroutine lubksb 

  
      SUBROUTINE jacobi(a,n,np,d,v,nrot)  
      integer(ik) n,np,nrot,NMAX  
      REAL (rk)::a(np,np),d(np),v(np,np)  
      PARAMETER (NMAX=500)  
      integer(ik) i,ip,iq,j  
      REAL (rk)::c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)  
      do 12 ip=1,n  
        do 11 iq=1,n  
          v(ip,iq)=0._rk  
11      continue  
        v(ip,ip)=1._rk  
12    continue  
      do 13 ip=1,n  
        b(ip)=a(ip,ip)  
        d(ip)=b(ip)  
        z(ip)=0._rk  
13    continue  
      nrot=0  
      do 24 i=1,50  
        sm=0._rk  
        do 15 ip=1,n-1  
          do 14 iq=ip+1,n  
            sm=sm+abs(a(ip,iq))  
14        continue  
15      continue  
        if(sm.eq.0._rk)return  
        if(i.lt.4)then  
          tresh=0.2_rk*sm/n**2  
        else  
          tresh=0._rk  
        endif  
        do 22 ip=1,n-1  
          do 21 iq=ip+1,n  
            g=100._rk*abs(a(ip,iq))  
            if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then  
              a(ip,iq)=0._rk  
            else if(abs(a(ip,iq)).gt.tresh)then  
              h=d(iq)-d(ip)  
              if(abs(h)+g.eq.abs(h))then  
                t=a(ip,iq)/h  
              else  
                theta=0.5_rk*h/a(ip,iq)  
                t=1./(abs(theta)+sqrt(1.+theta**2))  
                if(theta.lt.0._rk)t=-t  
              endif  
              c=1._rk/sqrt(1._rk+t**2)  
              s=t*c  
              tau=s/(1._rk+c)  
              h=t*a(ip,iq)  
              z(ip)=z(ip)-h  
              z(iq)=z(iq)+h  
              d(ip)=d(ip)-h  
              d(iq)=d(iq)+h  
              a(ip,iq)=0._rk  
              do 16 j=1,ip-1  
                g=a(j,ip)  
                h=a(j,iq)  
                a(j,ip)=g-s*(h+g*tau)  
                a(j,iq)=h+s*(g-h*tau)  
16            continue  
              do 17 j=ip+1,iq-1  
                g=a(ip,j)  
                h=a(j,iq)  
                a(ip,j)=g-s*(h+g*tau)  
                a(j,iq)=h+s*(g-h*tau)  
17            continue  
              do 18 j=iq+1,n  
                g=a(ip,j)  
                h=a(iq,j)  
                a(ip,j)=g-s*(h+g*tau)  
                a(iq,j)=h+s*(g-h*tau)  
18            continue  
              do 19 j=1,n  
                g=v(j,ip)  
                h=v(j,iq)  
                v(j,ip)=g-s*(h+g*tau)  
                v(j,iq)=h+s*(g-h*tau)  
19            continue  
              nrot=nrot+1  
            endif  
21        continue  
22      continue  
        do 23 ip=1,n  
          b(ip)=b(ip)+z(ip)  
          d(ip)=b(ip)  
          z(ip)=0._rk  
23      continue  
24    continue  
      stop 'too many iterations in jacobi'  
      return  
      END SUBROUTINE jacobi      
             
      FUNCTION gasdev(idum)  
      integer,intent(inout)::idum        
      REAL(rk):: gasdev  
      integer(ik) iset  
      REAL(rk):: fac,gset,rsq,v1,v2
      SAVE iset,gset  
      DATA iset/0/  
      if (iset.eq.0) then  
1       v1=2._rk*ran1(idum)-1._rk  
        v2=2._rk*ran1(idum)-1._rk  
        rsq=v1**2+v2**2  
        if(rsq.ge.1._rk.or.rsq.eq.0._rk)goto 1  
        fac=sqrt(-2._rk*log(rsq)/rsq)  
        gset=v1*fac  
        gasdev=v2*fac  
        iset=1  
      else  
        gasdev=gset  
        iset=0  
      endif  
      return  
      END FUNCTION gasdev  

      FUNCTION ran1(idum)  
       integer,intent(inout)::idum        
      integer(ik) IA,IM,IQ,IR,NTAB,NDIV  
      REAL(rk):: ran1,AM,EPS,RNMX  
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,&  
     &NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7_rk,RNMX=1._rk-EPS)  
      integer(ik) j,k,iv(NTAB),iy  
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
      END FUNCTION ran1 

       FUNCTION ran2(idum) 
       integer,intent(inout)::idum        
       integer(ik) IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV 
       Real (rk) ran2,AM,EPS,RNMX 
       PARAMETER (IM1=2147483563,IM2=2147483399,AM=1._rk/IM1,IMM1=IM1-1,& 
      &IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,& 
      &NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7_rk,RNMX=1._rk-EPS) 
       integer(ik) idum2,j,k,iv(NTAB),iy 
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
       END FUNCTION ran2 
         


 SUBROUTINE eigsrt(d,v,n,np)  
              INTEGER :: n,np  
              REAL(rk):: d(np),v(np,np)  
              INTEGER(ik):: i,j,k  
              REAL(rk):: p   
              do 13 i=1,n-1  
                k=i  
                p=d(i)  
                do 11 j=i+1,n  
                  if(d(j).ge.p)then  
                    k=j  
                    p=d(j)  
                  endif  
        11      continue  
                if(k.ne.i)then  
                  d(k)=d(i)  
                  d(i)=p  
                  do 12 j=1,n  
                    p=v(j,i)  
                    v(j,i)=v(j,k)  
                    v(j,k)=p  
        12        continue  
                endif  
        13    continue  
              return  
              END  SUBROUTINE eigsrt

             subroutine write_xyz(nfile,nat,atyp,xx)
          integer(ik)::i,j,nfile,nat      
          real(rk)::xx(3,nat)
          character(len=2)::atyp(nat)

          write(nfile,*) nat
          write(nfile,*) 

         do i=1,nat
             write(nfile,1000) atyp(i),(xx(j,i),j=1,3)
         enddo

 1000    format(A,3(F20.6))
         return
        end subroutine write_xyz



        function per(x,b)
         real(rk)::x,b,per
         integer(ik)::nn

         nn=int(x/b)
         if (x.lt.0.0_rk) nn=nn-1
         per=x-nn*b

        return
        end function per        

       subroutine seed_cal(idum)
       integer(ik)   :: timeArray(3) 
       integer(ik),dimension(8)::value
       integer(ik)::idum

     !Seed calculation
     idum=-99
     call date_and_time(values=value)
     timeArray(1)=value(6)
     timeArray(2)=value(7)
     timeArray(3)=value(8)

       idum = -((timeArray(1)+1)*(timeArray(2)+1)*(timeArray(3)+1))
       idum = -INT(1.d3*ran1(idum)*timeArray(1)+1.d3*ran1(idum)&
              *timeArray(2)+1.d2*ran1(idum)*timeArray(3))

       return
       end subroutine seed_cal       

      subroutine sortbubble(nn,a,b)
      integer(ik) :: j  
      integer(ik), intent(in) :: nn
      real (rk), intent(inout) :: a(nn)
      integer (ik), intent(out) :: b(nn)
      real (rk) :: tmp,tmp2
      logical::swap

      do j=1,nn
        b(j)=j
      enddo      

      do410 : do
              swap=.false.
      do420 : do j=1,nn-1
      if(a(j) .gt. a(j+1)) then
          tmp=a(j)
          tmp2=b(j)
          a(j)=a(j+1)
          b(j)=b(j+1)
          a(j+1)=tmp
          b(j+1)=tmp2
          swap=.true.
      end if
      end do do420
          if (.not.swap) exit
      end do do410
!
      end subroutine sortbubble

      subroutine sortbubble_int(nn,a,b)
      integer(ik) :: j  
      integer(ik), intent(in) :: nn
      integer (ik), intent(inout) :: a(nn)
      integer (ik), intent(out) :: b(nn)
      integer (ik) :: tmp,tmp2
      logical::swap

      do j=1,nn
        b(j)=j
      enddo      

      do410 : do
              swap=.false.
      do420 : do j=1,nn-1
      if(a(j) .gt. a(j+1)) then
          tmp=a(j)
          tmp2=b(j)
          a(j)=a(j+1)
          b(j)=b(j+1)
          a(j+1)=tmp
          b(j+1)=tmp2
          swap=.true.
      end if
      end do do420
          if (.not.swap) exit
      end do do410

      return
      end subroutine sortbubble_int

      subroutine sortbubble_int1(nn,a)
      integer(ik) :: j  
      integer(ik), intent(in) :: nn
      integer (ik), intent(inout) :: a(nn)
      integer (ik) :: tmp
      logical::swap


      do410 : do
              swap=.false.
      do420 : do j=1,nn-1
      if(a(j) .gt. a(j+1)) then
          tmp=a(j)
          a(j)=a(j+1)
          a(j+1)=tmp
          swap=.true.
      end if
      end do do420
          if (.not.swap) exit
      end do do410
 
      return
      end subroutine sortbubble_int1


! subroutine listcheck(lista,nnei,diml,p,id)
! integer(ik),intent(in)::diml,lista(diml),nnei,p
! logical::id
! integer(ik)::i

!     id=.false.

!     do i=1,nnei
!        if (lista(i).eq.p) then
!          id=.true.
!          return
!        endif
!     enddo      
!  


!    return
!   end subroutine listcheck

    function rsq(x,y,box)
      real(rk),intent(in)::x(3),y(3),box(3)      
      real(rk)::rvec(3),rsq
      integer(ik)::k
    rvec=x-y

    if (maxval(box).gt.tol) then
       do k=1,3    
       rvec(k)=rvec(k)-box(k)*anint(rvec(k)/box(k)) 
       enddo 
    endif

    rsq=dot_product(rvec,rvec)
    return
    end function rsq

    function radius(x,y,box)
      real(rk)::x(3),y(3),box(3),radius            
      radius=dsqrt(rsq(x,y,box))
    return
    end function radius

    subroutine rsq_cal(x,y,box,rsq,rvec,ndev)
      real(rk),intent(in)::x(3),y(3),box(3)
      real(rk)::rvec(3),rsq
      integer(ik)::k,ndev(3)
    rvec=x-y

    if (maxval(box).gt.tol) then
       do k=1,3    
       ndev(k)=anint(rvec(k)/box(k)) 
       rvec(k)=rvec(k)-box(k)*ndev(k)
       enddo 
    endif

       rsq=dot_product(rvec,rvec)
    return
    end subroutine rsq_cal

    subroutine isot_expansion(xx,box,xi,bi,n,fac)
       integer(ik),intent(in)::n
       real(rk),intent(in)::xi(3,n),bi(3)
       real(rk),intent(inout)::xx(3,n),box(3)
       real(rk)::fac
       integer(ik)::i,j

           xx=xi*fac
           box=bi*fac 


       end subroutine isot_expansion

       function edev(e,e2,n)
        real(rk)::edev,e,e2,fact1
        integer(ik)::n

        e=e/n;e2=e2/n
        fact1=e2-e**2
        
        !if (fact1<0) then
        !write(*,*) 'Numero negativo','e2: ',e2,': e 2: ',e**2
        !endif
        edev=dsqrt(fact1)
        !edev=dsqrt(e2-e**2)

        return      
        end function edev

subroutine calculate_at_number(nn,a,box,volume,nhexver,nhexhor)
 !Given a number of atoms nn and a bond length a tries to calculate a graphene 
 !hexagonal sheet parameters of aproximately same length and width
 !It changes the number of atoms in order to complete it and it tries
real(rk),intent(out)::box(3),volume
integer(ik),intent(out)::nhexhor,nhexver
real(rk),intent(in)::a
real(rk)::distver,disthor
integer(ik)::natplus,natmin,k,nn

nhexhor=int(dsqrt(dfloat(nn)/dsqrt(3.0_rk)))
if (mod(nhexhor,2).ne.0) nhexhor=nhexhor+1
disthor=nhexhor*3.0_rk*a/2.0_rk

if (mod(nn,2*nhexhor).eq.0) then
nhexver=nn/(2*nhexhor)
distver=nhexver*dsqrt(3.0_rk)*a
else
natplus=nn;natmin=nn
do
natplus=natplus+1;natmin=natmin-1
if (mod(natplus,2*nhexhor).eq.0) then
nn=natplus
nhexver=nn/(2*nhexhor)
distver=nhexver*dsqrt(3.0_rk)*a
exit
endif
if (mod(natmin,2*nhexhor).eq.0) then
nn=natmin

nhexver=nn/(2*nhexhor)
distver=nhexver*dsqrt(3.0_rk)*a
exit
endif
enddo

endif

box(1)=disthor
box(2)=distver
box(3)=200.0_rk
volume=box(1)*box(2)*3.35_rk

if (Abs(distver/disthor-1.0_rk).gt.0.2) stop 'Couldn t get a square'

return
end subroutine calculate_at_number

 subroutine listcheck(lista,nnei,diml,p,id)
  integer(ik),intent(in)::diml,lista(diml),nnei,p
  logical::id
  integer(ik)::i

      id=.false.

      do i=1,nnei
         if (lista(i).eq.p) then
           id=.true.
           return
         endif
      enddo      
   


     return
    end subroutine listcheck

   subroutine write_vector(nfile,vect,n,m)
     integer(ik),intent(in)::nfile,n,m
     real(rk),intent(in)::vect(n,m)
     integer(ik)::i,j

     do i=1,n
      write(nfile,1000) (vect(i,j),j=1,m)
     enddo        

1000 format(100000(f25.10))
    end subroutine write_vector

      subroutine translationtozero(xx,bl,nn)
       integer(ik),intent(in)::nn
       real(rk)::minc(3)
       real(rk),intent(inout)::xx(3,nn)
       real(rk),intent(in)::bl(3)
       integer(ik)::i,j

       do j=1,3
       minc(j)=minval(xx(j,:))
       enddo

       do j=1,3
       xx(j,:)=xx(j,:)-minc(j)
       if (minval(xx(j,:)).lt.0.0_rk) stop 'lt 0 in translationtozero'
       if (maxval(xx(j,:)).gt.bl(j)) then
               write(*,*) bl(j),maxval(xx(j,:))
               stop 'gt bl(j) in translationtozero'
       endif
       enddo

       

      return
      end subroutine translationtozero

      subroutine calculate_gdr(nat,box)
      integer(ik),intent(in)::nat
      real(rk)::dens
      real(rk)::box(3)      
      integer(ik)::ierror,i,j,k,conf,ndiv,nn
      integer(ik),allocatable::A(:)
      real(rk)::r,rmax,gdr,deltar,pi,vol
      real(rk),allocatable::xx(:,:)
      character(len=2)::atyp      

      ndiv=250
      allocate(A(ndiv))
      allocate(xx(3,nat))
      A=0
      dens=0.0_rk

      rmax=minval(box)/2.0_rk
      pi=dacos(-1.0_rk)
      deltar=rmax/dfloat(ndiv)
      conf=0

      rewind(100)
      do
       read(100,*,IOSTAT=ierror) nn
       if (ierror.ne.0) exit
       if (nn.ne.nat) stop 'error in nat lecture for gdr calculation'
       conf=conf+1
       if (mod(conf,100).eq.0) write(*,*) 'configuration ',conf
       read(100,*) box
       vol=box(1)*box(2)*box(3)
       dens=dens+dfloat(nat)/vol

       do i=1,nat
        read(100,*) atyp,xx(:,i)
       enddo

       do i=1,nat-1
               do j=i+1,nat
                 r=rsq(xx(:,i),xx(:,j),box)
                 r=dsqrt(r)
                 if (r.ge.rmax) cycle
                 k=int(r*ndiv/rmax)+1
                 A(k)=A(k)+1
               enddo
       enddo
      enddo

      write(*,*) 'nconf=',conf
      dens=dens/conf

      do i=1,ndiv
         r=dfloat(i)*rmax/dfloat(ndiv)     
         vol=4.0_rk*pi*((r+deltar)**3-r**3)/3.0_rk
         gdr=2.0_rk*dfloat(A(i))/(nat*conf*vol*dens)
         write(600,1000) r,gdr
      enddo

1000 format(*(F18.10,2x))


      end subroutine calculate_gdr        

      subroutine at_interchange(nhexhor,nhexver,xc,n_at,box,stinval)
      integer,intent(in)::n_at,nhexhor,nhexver
      real(rk),intent(in)::box(3),stinval
      real(rk)::acop(3,n_at),xc(3,n_at),xt(3),rr
      integer(ik)::i,j
      integer(ik)::k,kk,kkk,idelta

      i=nhexhor*(nhexver)+nhexver
      j=i+2*nhexver
      rr=radius(xc(:,i),xc(:,j),box)
      if (rr.gt.2.0_rk*(1.0_rk+stinval)) j=i-2*nhexver
      if (abs(xc(2,i)-xc(2,j)).gt.tol) stop 'Bad choice of atoms in ts'
      rr=radius(xc(:,i),xc(:,j),box)
      if (rr.gt.2.0_rk*(1.0_rk+stinval)) stop 'Bad choice of atoms in ts big rr' 
      write(*,*) 'intercambio',i,j

       acop=xc
       xc(:,1)=xc(:,i)
       xc(:,2)=xc(:,j)

       if (xc(1,1).gt.xc(1,2)) then
         xt=xc(:,1)
         xc(:,1)=xc(:,2)
         xc(:,2)=xt
       endif      
       


       kk=0
       do k=3,n_at
        idelta=2
        if (k.ge.min(i,j)+2) idelta=1
        if (k.ge.max(i,j)+1) idelta=0
        xc(:,k)=acop(:,k-idelta)
       enddo


      return
      end subroutine at_interchange

      subroutine at_rotation(i,j,angrot,xc,n_at)
      integer,intent(in)::n_at
      real(rk)::xrot(3,2),xrotc(3,2),cm(3),angrot,xc(3,n_at)
      integer(ik),intent(in)::i,j
      integer(ik)::k

       xrot(:,1)=xc(:,i)
       xrot(:,2)=xc(:,j)
       cm=(xrot(:,1)+xrot(:,2))/2.0_rk
       
       do k=1,2
       xrot(k,:)=xrot(k,:)-cm(k) !Translation to the zero
       enddo

       xrotc=xrot

       do k=1,2
       xrot(1,k)=xrotc(1,k)*dcos(angrot)+xrotc(2,k)*dsin(angrot)  
       xrot(2,k)=-xrotc(1,k)*dsin(angrot)+xrotc(2,k)*dcos(angrot)  !Rotation
       enddo    

       do k=1,2
       xrot(k,:)=xrot(k,:)+cm(k) !Translation to the original position once rotated
       enddo
       xc(:,1)=xrot(:,1)
       xc(:,2)=xrot(:,2)


      return
      end subroutine at_rotation
      
      function cross(a,b)
      real (rk) ::a(3),b(3),cross(3)
       cross(1)=a(2)*b(3)-a(3)*b(2)
       cross(2)=a(3)*b(1)-a(1)*b(3)
       cross(3)=a(1)*b(2)-a(2)*b(1)
       return
      end function cross

    function cross_unitary(a,b)
    real (rk) ::a(3),b(3),cross_unitary(3),anorm
    cross_unitary=0._rk

    cross_unitary(1)=a(2)*b(3)-a(3)*b(2)
    cross_unitary(2)=a(3)*b(1)-a(1)*b(3)
    cross_unitary(3)=a(1)*b(2)-a(2)*b(1)
    anorm=norm(cross_unitary)
    if (anorm.lt.1.e-20_rk) stop 'problems calculating cross_unitary'
    cross_unitary=cross_unitary/anorm
    return
    end function cross_unitary

     function norm(a)
             real(rk)::a(3),norm

             norm=dsqrt(dot_product(a,a))


             end function norm

             function rad(a,b)
              real(rk)::rad,a(3),b(3),c(3)

               c=a-b
               rad=dsqrt(dot_product(c,c))

              end function rad


      function deg(x)
       real(rk)::pi,deg,x
        pi=dacos(-1.0_rk)

        deg=x*180.0_rk/pi

       end function deg       


   subroutine histogram(vec,ndim,ndiv,hist)
    !Given a vector of dimension ndim, it creates and histogram
    !taking ndivisions from the smallest to the biggest. The result
    !is stored on vector hist
    integer(ik)::ndim,ndiv,i,j,k
    real(rk)::vec(ndim),hist(2,ndiv),deltax,maxv,minv,delta

    maxv=maxval(vec)
    minv=minval(vec)
    delta=maxv-minv

    deltax=(maxv-minv)/ndiv
    hist=0.0_rk

    do i=1,ndiv
     hist(1,i)=minv+(i-0.5_rk)*deltax
    enddo        

    do i=1,ndim
            if (abs(maxv-vec(i)).lt.tol) then
            j=ndiv         
            else
            j=int(ndiv*((vec(i)-minv)/delta))+1
            if (j.gt.ndiv) stop 'error calculating histogram'
            endif
            hist(2,j)=hist(2,j)+1.0_rk
    enddo        


    return
   end subroutine histogram        

             subroutine rotz(rij,th)
               real(rk)::rij(3),th,Amat(3,3),rijt(3)
               integer(ik)::i,j

               amat(1,:)=(/dcos(th),dsin(th),0.0_rk/)
               amat(2,:)=(/-dsin(th),dcos(th),0.0_rk/)
               amat(3,:)=(/0.0_rk,0.0_rk,1.0_rk/)
               rijt=0.0_rk

               do i=1,3
                       do j=1,3
                         rijt(i)=rijt(i)+amat(i,j)*rij(j)
                       enddo
               enddo

                  rij=rijt      
                               

              return
              end subroutine rotz      

             subroutine roty(rij,th)
               real(rk)::rij(3),th,Amat(3,3),rijt(3)
               integer(ik)::i,j

               amat(1,:)=(/dcos(th),0.0_rk,dsin(th)/)
               amat(2,:)=(/0.0_rk,1.0_rk,0.0_rk/)
               amat(3,:)=(/-dsin(th),0.0_rk,dcos(th)/)
               rijt=0.0_rk

               do i=1,3
                       do j=1,3
                         rijt(i)=rijt(i)+amat(i,j)*rij(j)
                       enddo
               enddo

                  rij=rijt      
                               

              return
              end subroutine roty      
 
end module math_subroutines
