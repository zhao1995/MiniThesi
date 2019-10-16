      subroutine elmt23(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c-----------------------------------------------------------------------
c     quadrilateral linear IDKQ plate element                          
c
c     theory: 
c     1) C. Jeychandrabose, J. Kirkhope, L. Meekisho
c        An improved discrete Kirchhoff quadrilateral thin-plate bending
c        element, IJNME 1987; 24:635–654

c      2) S. Kapuria and S. D. Kulkarni
c         An improved discrete Kirchhoff quadrilateral element based
c         on third-order zigzag theory for static analysis of composite
c         and sandwich plates, IJNME 2007; 69:1948–1981(in the appendix) 
c
c     comments ww 
c         theory is for dofs w, w,_x, w,_y 
c         implementation for FEAP is w, w,_y, w,_x !
c         -> modification in SR BMATIDKQ
c
c     open
c         consistant load vector
c
c     Material parameters (input)                                      
c      1: E   = Young's modulus                                        
c      2: nu  = Poisson's ratio                                        
c      3: h   = thickness of shell                                     
c      4: q   = transverse load                                         
c
c     (c) WW UKA, FG TUD  4/2007
c
c-----------------------------------------------------------------------
      USE bdata      
      USE cdata
      USE eldata
      USE iofile
      USE pdata7
      USE prlod
      USE qload      
      USE strnam 
      implicit double precision (a-h,o-z)
      dimension xl(ndm,4),tl(*),ix(*),ul(ndf,4),s(nst,nst),p(nst),d(*),
     1    sg(4),tg(4),wg(4),gp(2),xjaci(2,2),gg(2,12),hg(2,12),db(3,3),
     2    sig(3),bb(3,12),shp(4)
      dimension h1(*),h2(*),h3(*)

c.... go to correct array processor
      go to(1,2,3,4,2,3,2,4,2,2,2,2,2,2,2,2,2,2,2,2,2,3) isw
c.... input material properties
1     if(ior.lt.0) write(*,1001)
      call dinput(d,4)
      write(iow,1002) (d(i),i=1,4)
      if(ior.lt.0) write(*,1002) (d(i),i=1,4)
      if(ndm.eq.3) ipla = 1
c.... description of stresses  
      strsus( 1) = '  MOMENT M_xx  '
      strsus( 2) = '  MOMENT M_xy  '
      strsus( 3) = '  MOMENT M_yy  '
      do is =4,25
        strsus(is) = '               '
      end do
2     return

3     call idkqval(xl,ndm)                   ! basic values
      l  = 2
      call pgauss(l,lint,sg,tg,wg) ! Gauss-Points
      call modidkq(db,d)                     ! elasticity matrix  
      do l = 1,lint
        xsi = sg(l)
        eta = tg(l)
        call jacoidkq(djacb,xjaci,xsi,eta,n) ! jacobian
        da = djacb*wg(l)
        call idkqn(xsi,eta,gg,hg,xjaci)      ! G,_i  H,_i
        call bmatidkq(bb,gg,hg)              ! B 
        call subidkq(bb,da,db,s,nst)         ! Int B^T*D*B*da
        call shpidkq(shp,xsi,eta)            ! N_i  
        call qloadidkq(qz,d,aqloa,numel,n,mqloa,propq,prop,isw) ! P
        do ino =1,4
          iq = (ino-1)*3+1
          p(iq) = p(iq) + shp(ino)*qz*da
        end do   
      end do
      call residkq(p,s,ul,nst)               ! G

      return

c.....stress resultants 
4     istv = -3
      call idkqval(xl,ndm)                   ! basic values
      l  = 2
      call pgauss(l,lint,sg,tg,wg)           ! Gauss-Points
      call modidkq(db,d)                     ! elasticity matrix  
      do l = 1,lint
        xsi = sg(l)
        eta = tg(l)
        call jacoidkq(djacb,xjaci,xsi,eta,n) ! jacobian
        da = djacb*wg(l)
        call idkqn(xsi,eta,gg,hg,xjaci)      ! G,_i  H,_i
        call bmatidkq(bb,gg,hg)              ! B 
        call stridkq(db,bb,ul,sig)           ! stresses 
        call shpidkq(shp,xsi,eta)            ! N_i
        if(isw.eq.4) then ! print M
          do 40 idm = 1,2                    ! GP-coordinates
            gp(idm) = 0.0
            do 40 ino = 1,4
   40       gp(idm) = gp(idm) + xl(idm,ino) * shp(ino) 
            mct = mct - 1
            if(mct.gt.0) go to 41
                         write(iow,2001) o,head
            if(ior.lt.0) write(*  ,2001) o,head
            mct = 50
   41                write(iow,2002) n,ma,(gp(i),i=1,2),(sig(i),i=1,3)
        if(ior.lt.0) write(*  ,2002) n,ma,(gp(i),i=1,2),(sig(i),i=1,3)
        else if(isw.eq.8) then ! plot M
          if(iplma(ma).eq.0)  return ! only if MATN
          call plotidkq(ix,strea,strea(1+numnp),shp,sig,da,numnp)
        end if  
      end do
      return
c
c.... formats
1001  format(' Input: E  nu  h  q',$)
1002  format(5x,'Material data IDKQ plate element:',/,
     + 5x,'elastic modulus............',   f15.4,/,
     + 5x,'poissons ratio................',f12.4,/,
     + 5x,'thickness.....................',f12.4,/,
     + 5x,'element load (transverse).....',f12.4)
2001  format(a1,20a4,//,2x,'E L E M E N T   S T R E S S E S',//,
     1 2x,'EL',1x,'MAT',1x,'1-COORD',1x,'2-COORD',
     2 2X,'   MX   ',3X,'   MY   ',3X,'   MXY   ',/)
2002  format(1x,i3,i3,2f8.3,1x,3e12.5)
      end
c
      subroutine bmatidkq(bb,gg,hg)
c-----------------------------------------------------------------------
c.....B-matrices for idkq-plate element
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  bb(3,12),gg(2,12),hg(2,12)

      do i = 1,12
        k = i
        if(i.eq.2.or.i.eq.5.or.i.eq.8.or.i.eq.11) k=i+1
        if(i.eq.3.or.i.eq.6.or.i.eq.9.or.i.eq.12) k=i-1
 
        bb(1,i) = gg(1,k)
        bb(2,i) = hg(2,k)
        bb(3,i) = gg(2,k) + hg(1,k)

c        bb(1,i) = gg(1,i)
c        bb(2,i) = hg(2,i)
c        bb(3,i) = gg(2,i) + hg(1,i)
      end do
      return
      end
c
      subroutine stridkq(db,bb,u,sig)
c-----------------------------------------------------------------------
c.....Strains and stresses for IDKQ-plate
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension bb(3,12),db(3,3),dbj(3,12),u(12),sig(3)
c.....dbj=d*bj
      call pzero(dbj,36)
      do 10 i = 1,3
        do 10 j = 1,12
          do 10 k = 1,3
   10     dbj(i,j) = dbj(i,j) + db(i,k) * bb(k,j)
c.... stresses
      call pzero(sig,3)
      do 20 i = 1,3
         do 20 j=1,12
           sig(i) = sig(i) + dbj(i,j)*u(j) 
   20      continue
      return
      end
c
      subroutine jacoidkq(djacb,xjaci,ss,tt,n)
c-----------------------------------------------------------------------
c.....Jacobi-Matrix 
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /idkq23/ x(4),y(4),b(4),c(4),p(4),q(4),r(4),s(4),t(4),dn
!$OMP THREADPRIVATE (/idkq23/)  
      dimension xjacm(2,2),xjaci(2,2)

c.....Jacobi-matrix xjacm
      xjacm(1,1)=0.25d0*(x(2)-x(1)+x(3)-x(4) + tt*(x(1)-x(2)+x(3)-x(4)))
      xjacm(1,2)=0.25d0*(y(2)-y(1)+y(3)-y(4) + tt*(y(1)-y(2)+y(3)-y(4)))
      xjacm(2,1)=0.25d0*(x(3)-x(2)+x(4)-x(1) + ss*(x(1)-x(2)+x(3)-x(4)))
      xjacm(2,2)=0.25d0*(y(3)-y(2)+y(4)-y(1) + ss*(y(1)-y(2)+y(3)-y(4)))

c.....Determinant  and Inverse of Jacobi-Matrix
      djacb=xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1)
      if(djacb) 1,1,2
    1   write(*,100) n
        stop
    2   xjaci(1,1)= xjacm(2,2)/djacb
        xjaci(2,2)= xjacm(1,1)/djacb
        xjaci(1,2)=-xjacm(1,2)/djacb
        xjaci(2,1)=-xjacm(2,1)/djacb
  100 format(' SR jacoidkq: zero or negative area for element ',i5)
      return
      end
c
      subroutine modidkq(db,d)
c-----------------------------------------------------------------------
c.....Elasticity matrix for plate 
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension db(3,3),d(*)
      e = d(1)
      v = d(2)
      h = d(3)
      c = e*h*h*h/(12.d0*(1.0d0-v*v))
      call pzero(db,9)
      db(1,1) = c
      db(2,2) = c
      db(1,2) = c*v
      db(2,1) = c*v
      db(3,3) = c*(1.0d0-v)*0.5d0
      return
      end
c
      subroutine subidkq(bb,da,db,s,nst)
c-----------------------------------------------------------------------
c.....K = Int B^T*D*B*da
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension bb(3,12),db(3,3),dbj(3,12),s(nst,nst)
c.....dbj=d*bj
      call pzero(dbj,36)
      do 10 i = 1,3
        do 10 j = 1,12
          do 10 k = 1,3
   10     dbj(i,j) = dbj(i,j) + db(i,k) * bb(k,j)
c.....bi T *(d*bj)
      do 20 i = 1,12
        do 20 j = 1,12
          do 20 k = 1,3
   20     s(i,j) = s(i,j) + bb(k,i) * dbj(k,j) * da
      return
      end
c
      subroutine residkq(p,s,u,nst)
c-----------------------------------------------------------------------
c.....G = p-K*v
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension p(nst),s(nst,nst),u(12)
      do i = 1,12
        do j = 1,12
          p(i) = p(i) - s(i,j) * u(j)
        end do
      end do
      return
       end

      subroutine shpidkq (shp,s,t)
c-----------------------------------------------------------------------
c.....Shape functions for linear isoparametric  2-d element
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension shp(4)
      st = s*t
      shp(1) = (1.d0-t-s+st)*0.25d0
      shp(2) = (1.d0-t+s-st)*0.25d0
      shp(3) = (1.d0+t+s+st)*0.25d0
      shp(4) = (1.d0+t-s-st)*0.25d0
      return
      end
c
      subroutine plotidkq(ix,dt,st,shp,sig,da,numnp)
c-----------------------------------------------------------------------
c.....Plot   mx(1)   mxy(2)  my(3) 
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension dt(numnp),st(numnp,*),shp(4),sig(3),ix(1)
        do 10  j = 1,4
           xsji = da*shp(j)
           ii = abs(ix(j))
           if(ii.eq.0) go to 10
             dt(ii) = dt(ii) + xsji
             st(ii,1) = st(ii,1) + sig(1)*xsji
             st(ii,2) = st(ii,2) + sig(3)*xsji
             st(ii,3) = st(ii,3) + sig(2)*xsji
10      continue
      return
      end
c      
      subroutine idkqval(xl,ndm)
c-----------------------------------------------------------------------
c.....idkq-values
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension xl(ndm,4)
      common /idkq23/ x(4),y(4),b(4),c(4),p(4),q(4),r(4),s(4),t(4),dn 
!$OMP THREADPRIVATE (/idkq23/)  

      do i = 1,4
        x(i) = xl(1,i) 
        y(i) = xl(2,i)
      end do   

      do i = 1,4
        k=i+1
        if(i.eq.4) k=1
        b(i) = y(i)-y(k)
        c(i) = x(k)-x(i)  
        dno  = b(i)*b(i)+c(i)*c(i)
        p(i) = 6.d0*c(i)/dno
        q(i) = 3.d0*b(i)*c(i)/dno
        r(i) = 3.d0*b(i)*b(i)/dno
        s(i) = 3.d0*c(i)*c(i)/dno
        t(i) = 6.d0*b(i)/dno
      end do  
      
      dn = x(1)*(y(2)-y(4)) 
     +   + x(2)*(y(3)-y(1)) 
     +   + x(3)*(y(4)-y(2)) 
     +   + x(4)*(y(1)-y(3))
c
      return
      end
c
      subroutine idkqn(ss,tt,gg,hg,xjaci)
c-----------------------------------------------------------------------
c..... local derivatives of shape functions N_i
c....  local derivatives of G, H
c.... global derivatives of G, H
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /idkq23/ x(4),y(4),b(4),c(4),p(4),q(4),r(4),s(4),t(4),dn
!$OMP THREADPRIVATE (/idkq23/)  
      real*8 n(2,9),ns,nt 
      dimension gg(2,12),hg(2,12),xjaci(2,2)
c.... constants
      d2 = 0.50d0
      d3 = 1.d0/3.d0
      d4 = 0.25d0
      d6 = 1.d0/6.d0
      d23= 2.d0/3.d0
      d32= 3.d0/2.d0
      d43= 4.d0/3.d0

c.....N_i,s
      n(1,1) =  d6*ss+d4*tt-d2*ss*tt-d4*tt*tt
      n(1,2) =  d6*ss-d4*tt-d2*ss*tt+d4*tt*tt
      n(1,3) =  d6*ss+d4*tt+d2*ss*tt+d4*tt*tt
      n(1,4) =  d6*ss-d4*tt+d2*ss*tt-d4*tt*tt

      n(1,5) = -d3*ss+ss*tt  
      n(1,6) =  d2+d23*ss-d2*tt*tt
      n(1,7) = -d3*ss-ss*tt  
      n(1,8) = -d2+d23*ss+d2*tt*tt          

      n(1,9) = -d43*ss
      
c.....N_i,t
      n(2,1) =  d4*ss+d6*tt-d4*ss*ss-d2*ss*tt
      n(2,2) = -d4*ss+d6*tt-d4*ss*ss+d2*ss*tt
      n(2,3) =  d4*ss+d6*tt+d4*ss*ss+d2*ss*tt
      n(2,4) = -d4*ss+d6*tt+d4*ss*ss-d2*ss*tt

      n(2,5) = -d2+d23*tt+d2*ss*ss
      n(2,6) = -d3*tt-ss*tt  
      n(2,7) =  d2+d23*tt-d2*ss*ss          
      n(2,8) = -d3*tt+ss*tt  

      n(2,9) = -d43*tt

c.... cartesian derivatives of N

      do i=1,9
        ns = n(1,i) 
        nt = n(2,i) 
        n(1,i) = xjaci(1,1)*ns + xjaci(1,2)*nt
        n(2,i) = xjaci(2,1)*ns + xjaci(2,2)*nt
      end do

c.....coordinate differences
      x12 = x(1)-x(2)
      x13 = x(1)-x(3)
      x14 = x(1)-x(4)

      x21 = x(2)-x(1)
      x23 = x(2)-x(3)
      x24 = x(2)-x(4)

      x31 = x(3)-x(1)
      x32 = x(3)-x(2)
      x34 = x(3)-x(4)

      x41 = x(4)-x(1)
      x42 = x(4)-x(2)
      x43 = x(4)-x(3)

      y12 = y(1)-y(2)
      y13 = y(1)-y(3)
      y14 = y(1)-y(4)

      y21 = y(2)-y(1)
      y23 = y(2)-y(3)
      y24 = y(2)-y(4)

      y31 = y(3)-y(1)
      y32 = y(3)-y(2)
      y34 = y(3)-y(4)

      y41 = y(4)-y(1)
      y42 = y(4)-y(2)
      y43 = y(4)-y(3)

c.... cartesian derivatives of G,H 

      do k=1,2

c.....G_i,x, G_i,y 

      gg(k, 1)= d4*p(1)*n(k,5)-d4*p(4)*n(k,8)+d32*y42*n(k,9)/dn
      gg(k, 4)=-d4*p(1)*n(k,5)+d4*p(2)*n(k,6)+d32*y13*n(k,9)/dn
      gg(k, 7)=-d4*p(2)*n(k,6)+d4*p(3)*n(k,7)+d32*y24*n(k,9)/dn
      gg(k,10)=-d4*p(3)*n(k,7)+d4*p(4)*n(k,8)+d32*y31*n(k,9)/dn

      gg(k, 3)=-d4*q(1)*n(k,5)-d4*q(4)*n(k,8)+d4*y42*(y21+y41)*n(k,9)/dn
      gg(k, 6)=-d4*q(1)*n(k,5)-d4*q(2)*n(k,6)+d4*y13*(y12+y32)*n(k,9)/dn
      gg(k, 9)=-d4*q(2)*n(k,6)-d4*q(3)*n(k,7)+d4*y24*(y23+y43)*n(k,9)/dn
      gg(k,12)=-d4*q(3)*n(k,7)-d4*q(4)*n(k,8)+d4*y31*(y14+y34)*n(k,9)/dn

      gg(k, 2)=-n(k,1)+d4*(s(1)-2.d0)*n(k,5)+d4*(s(4)-2.d0)*n(k,8)
     +                +d4*y42*(x21+x41)*n(k,9)/dn
      gg(k, 5)=-n(k,2)+d4*(s(1)-2.d0)*n(k,5)+d4*(s(2)-2.d0)*n(k,6)
     +                +d4*y13*(x12+x32)*n(k,9)/dn 
      gg(k, 8)=-n(k,3)+d4*(s(2)-2.d0)*n(k,6)+d4*(s(3)-2.d0)*n(k,7)
     +                +d4*y24*(x23+x43)*n(k,9)/dn
      gg(k,11)=-n(k,4)+d4*(s(3)-2.d0)*n(k,7)+d4*(s(4)-2.d0)*n(k,8)
     +                +d4*y31*(x14+x34)*n(k,9)/dn

c.....H_i,x   H_i,y

      hg(k, 1)=-d4*t(1)*n(k,5)+d4*t(4)*n(k,8)+d32*x24*n(k,9)/dn    ! wie Kapuria, Umstimmigkeit zu Original
      hg(k, 4)= d4*t(1)*n(k,5)-d4*t(2)*n(k,6)+d32*x31*n(k,9)/dn
      hg(k, 7)= d4*t(2)*n(k,6)-d4*t(3)*n(k,7)+d32*x42*n(k,9)/dn
      hg(k,10)= d4*t(3)*n(k,7)-d4*t(4)*n(k,8)+d32*x13*n(k,9)/dn

      hg(k, 2)=-d4*q(1)*n(k,5)-d4*q(4)*n(k,8)+d4*x24*(x21+x41)*n(k,9)/dn
      hg(k, 5)=-d4*q(1)*n(k,5)-d4*q(2)*n(k,6)+d4*x31*(x12+x32)*n(k,9)/dn
      hg(k, 8)=-d4*q(2)*n(k,6)-d4*q(3)*n(k,7)+d4*x42*(x23+x43)*n(k,9)/dn
      hg(k,11)=-d4*q(3)*n(k,7)-d4*q(4)*n(k,8)+d4*x13*(x14+x34)*n(k,9)/dn

      hg(k, 3)=-n(k,1)+d4*(r(1)-2.d0)*n(k,5)+d4*(r(4)-2.d0)*n(k,8)
     +                +d4*x24*(y21+y41)*n(k,9)/dn
      hg(k, 6)=-n(k,2)+d4*(r(1)-2.d0)*n(k,5)+d4*(r(2)-2.d0)*n(k,6)
     +                +d4*x31*(y12+y32)*n(k,9)/dn       
      hg(k, 9)=-n(k,3)+d4*(r(2)-2.d0)*n(k,6)+d4*(r(3)-2.d0)*n(k,7)
     +                +d4*x42*(y23+y43)*n(k,9)/dn
      hg(k,12)=-n(k,4)+d4*(r(3)-2.d0)*n(k,7)+d4*(r(4)-2.d0)*n(k,8)
     +                +d4*x13*(y14+y34)*n(k,9)/dn
 
      end do

      return
      end
c
      subroutine qloadidkq(qz,d,q,numel,n,mqloa,propq,prop,isw)
c----------------------------------------------------------
c.... set loads from macro qloa/mate
c----------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension q(numel,10),d(*)
      if(isw.eq.22) then
        qz=0.d0 
        if(mqloa.ne.1) qz = q(n,1)*propq 
      else
        qz = d(4)*prop
      end if
      return
      end
      