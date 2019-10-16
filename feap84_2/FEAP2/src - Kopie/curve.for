      subroutine curve ( x,y,z,m,nr,diff)
c-----------------------------------------------------------------------
c
c....   Purpose: calculation of defined curves
C                ityp 1  to 50  : standard curves     ( same in all feap version )
C                ityp 51 to 100 : user defined curves ( only in user version )
c
C     input :   x    = argument
C               y    = argument
C               z    = argument
C               nr   = curve number
c
C     output:   tol  = value  (only for 3-D)
C               y    = value  at x   (2-d)
C               z    = value  at x,y (3-d)
C               diff = 0 if x,y,z lies on curve (else diff<>0)
c
C     variables cpar = curve parametres
c
c-----------------------------------------------------------------------

      USE cdata
      USE iofile
      USE tdata
      implicit double precision (a-h,o-z)
      common /curvedat/cpar(20,8),nrt1(20,3),ic,nbe,nn3
      dimension xx(3),yy(3),xy(3)
      dimension m(*)
cpf   data pi /3.14159265359/
c
      do 100 i=1,20                   ! max 20 curves
        if( nrt1(i,1) .eq. nr ) then
          ityp = nrt1(i,2)             ! curve typ
          nn3  = nrt1(i,3)
         goto 200
        endif
100   continue
                   write(iow,2000) nr     ! error curve not input
      if(ior.lt.0) write(*,  2000) nr
      stop
200   continue
c
c-----------------------------------------------------------------------
c                  2 - D  Kurven
c-----------------------------------------------------------------------
C
      if ( ityp .eq. 1 ) then
c
c.... Gerade durch zwei Punkte (x1,y1) (x2,y2)
c
        x1 = cpar(nr,1)
        y1 = cpar(nr,2)
        x2 = cpar(nr,3)
        y2 = cpar(nr,4)
czr     a = (y1*x2-y2*x1) / (x2-x1)
czr     b = (y2-y1) / (x2-x1)
czr     y0 = a + b*x
czr>
        if (abs(x1-x2).gt.0.00001) then
          a = (y1*x2-y2*x1) / (x2-x1)
          b = (y2-y1) / (x2-x1)
cpf       y = a + b*x
          yc= a + b*x
        elseif ((abs(x-x2).lt.0.00001).or.(abs(x-x1).lt.0.00001)) then
cpf       y = z
          yc= y
        else
cpf       y = z + 100
          yc= 99999999.99
        endif
czr<
czr     z = b
cpf>
        diff=y-yc
cpf<
c
      elseif ( ityp .eq. 2 ) then
c
c.... Gerade durch einen Punkt (x1,y1) und Steigung b
c
        x1 = cpar(nr,1)
        y1 = cpar(nr,2)
        b  = cpar(nr,3)
        a  = y1 - b*x1
        yc = a + b*x
        z  = b
cpf>
        diff=y-yc
cpf<
c
      elseif ( ityp .eq. 3 ) then
c
C.... Halbkreis mit dem Radius r, Steuerparameter nx
c     und Mittelpunkt (x1,y1)
c
        r  = cpar(nr,1)
        xv = cpar(nr,2)
czr     if(abs(xv).ne.1.d0) xv=1.d0
cpf     if( .not.(abs(xv).lt.1.d0.or.abs(xv).gt.1.d0)) xv=1.d0
        x0 = cpar(nr,3)
        y0 = cpar(nr,4)
        if(x.gt.(x0+r) .or. x.lt.(x0-r)) then
          yc= 99999999.d0
        else
          yc= y0 + xv* dsqrt(r*r-(x-x0)*(x-x0))
cpf   diff = dsqrt((x-x0)**2+(y-y0)**2+(z-z0)**2) - r
        endif
        diff=yc-y
c
c-----------------------------------------------------------------------
c                   3 - D  Kurven
c-----------------------------------------------------------------------
C
      elseif ( ityp .eq. 32 ) then
c
c... Kugel mit dem Radius r und Mittelpunkt (x1,y1,z1)
c
        r  = cpar(nr,1)
        x0 = cpar(nr,2)
        y0 = cpar(nr,3)
        z0 = cpar(nr,4)
        diff = dsqrt((x-x0)**2+(y-y0)**2+(z-z0)**2) - r
c
      elseif ( ityp .eq. 31 ) then
c
C.... Ebene mit der Normalen x0 y0 z0 und der Konstanten c0
c         (x0*x+y0*y+z0*z+c0=0)
c
        x0 = cpar(nr,1)
        y0 = cpar(nr,2)
        z0 = cpar(nr,3)
        c0 = cpar(nr,4)
c
        diff = x*x0+y*y0+z*z0+c0
c
      elseif ( ityp .eq. 30 ) then
c
C.... Gerade mit den 2 Punkten x1,y1,z1 und x2,y2,z2,
c
        x1 = cpar(nr,1)
        y1 = cpar(nr,2)
        z1 = cpar(nr,3)
        x2 = cpar(nr,4)
        y2 = cpar(nr,5)
        z2 = cpar(nr,6)
        xx(1)=x2-x1
        xx(2)=y2-y1
        xx(3)=z2-z1
c        call vnorm (xx,xl)
        yy(1)=x-x1
        yy(2)=y-y1
        yy(3)=z-z1
        call vcross(xx,yy,xy)
        if (dabs(xy(1))+dabs(xy(2))+dabs(xy(3)).lt.1e-5) then
          diff=0.d0
        else
          diff=10.d0
        endif
c
      elseif ( ityp .eq. 50 ) then
c
c.... curve defined via points
       call ctyp50 (x,y,m(nbe),z,nn3)
C
c.... user defined curves
      elseif ( ityp.gt.50 .and. ityp.le.100 ) then
        call cvuser ( x,y,z,ityp,m,nr,diff )
      else
                     write(iow,2001) ityp
        if(ior.lt.0) write(*,2001) ityp
        stop
      endif
      return
2000  format('***ERROR*** curve number ',i3,' in not been input')
2001  format('***ERROR*** curve type ',i3,' in not defined')
      end

      subroutine curvin ( m )
c-----------------------------------------------------------------------
c
c.... Purpose:  this SR read the  curve-parameters from input file
c
c     cxy = array in blank common for curve data
c     nn1 = number of the curve
c     nn2 = type   of the curve
c     nn3 = number of points to read (for type 50)
c     cpar : max. 8 Parameter fÅr jede Kurve
c     nrt1 : Feld fÅr die Parameter nn1,nn2,nn3, fÅr jede Kurve
c
c-----------------------------------------------------------------------
      USE cdata
      USE iofile
      implicit double precision (a-h,o-z)
c
      logical ldummy
      common /curvedat/cpar(20,8),nrt1(20,3),ic,nbe,nn3
      dimension m(*),di(8)
c      call pzero ( cpar,20*8 )   !  in SR ini
c      call pzeroi( nrt1,20*3 )   !  in SR ini
10    continue
      ic = ic+1
      if(ic.gt.20) then                 ! max. number of curves is 20
                     write (iow,2100) 20
        if(ior.lt.0) write (*,  2100) 20
        stop
      endif
      call dinput ( di,3 )
      nn1 = di(1)
      nn2 = di(2)
      nn3 = di(3)
      if(nn1.eq.0) goto 501 ! return
czr                   write (iow,1020) nn1,nn2,nn3
      if(ior.lt.0) write (*,  1020) nn1,nn2,nn3
      nrt1(ic,1) = nn1
      nrt1(ic,2) = nn2
      nrt1(ic,3) = nn3
      if (nn2.gt.0  .and. nn2.le.49) then     ! standard curve
        call dinput ( di,8 )
        do i=1,8
          cpar(ic,i) = di(i)
        enddo
c.... print information
        if (nn2.eq. 1) write(iow,2001)
     +      (nrt1(ic,i),i=1,2),(cpar(ic,i),i=1,4)
        if (nn2.eq. 2) write(iow,2002)
     +      (nrt1(ic,i),i=1,2),(cpar(ic,i),i=1,3)
        if (nn2.eq. 3) write(iow,2003)
     +      (nrt1(ic,i),i=1,2),(cpar(ic,i),i=1,3)
        if (nn2.eq.30) write(iow,2030)
     +      (nrt1(ic,i),i=1,2),(cpar(ic,i),i=1,6)
        if (nn2.eq.31) write(iow,2031)
     +      (nrt1(ic,i),i=1,2),(cpar(ic,i),i=1,4)
        if (nn2.eq.32) write(iow,2032)
     +      (nrt1(ic,i),i=1,2),(cpar(ic,i),i=1,4)
        goto 10
      elseif ( nn2.eq.50 ) then  ! point curve (only one per data file)
        call pseta (nbe,2*nn3,ipr,ldummy,'CURVE')
        call pzero ( m(nbe),nn3*2 )
        call inp50 ( m(nbe),nn3 )
        goto 10
      elseif ( nn2.ge.51 .and. nn2.le.100 ) then   ! user  curves
        call dinput ( di,8 )
        do i=1,8
          cpar(ic,i) = di(i)
        enddo
c->     call to user info
        goto 10
      else
        write(iow,2000) nn2       ! error: curve is not defined
        write(*,2000)   nn2
        stop
      endif
c....
501   continue
      write(iow,*)
      do i=1,ic
        if (nrt1(i,3).eq.15) write(iow,3015) nrt1(i,1)
        if (nrt1(i,3).eq.16) write(iow,3016) nrt1(i,1)
        if (nrt1(i,3).eq.17) write(iow,3017) nrt1(i,1)
        if (nrt1(i,3).eq.18) write(iow,3018) nrt1(i,1)
        if (nrt1(i,3).ge.20 .and. nrt1(i,3).le.27 )
     +      write(iow,3020) nrt1(i,3)-20,nrt1(i,1)
        if (nrt1(i,3).ge.30 .and. nrt1(i,3).le.37 )
     +      write(iow,3030) nrt1(i,1),nrt1(i,3)-30,(cpar(i,il),il=7,8)
      enddo
      write(iow,*)
      return
czr1030  format(8g12.5)
1020  format(3x,'parameter for curve',3x,3(i5,2x))
2000  format('***ERROR in curve input***, specified curve-typ ',i3,
     +       ' is not defined')
2100  format('***ERROR in curve input***, max number of curves is ',i3,
     +       '   source code has to modified')
czr>
2001  format(' curve #',i2,' type',i2,
     +      '  * linear 2D-function (2 points): '/
     +      2x,2(f6.2,','),2x,2(f6.2,','))
2002  format(' curve #',i2,' type',i2,
     +      '  * linear 2D-function (1 point, incline): '/
     +      2x,3(f6.2,','),2x,1(f6.2,','))
2003  format(' curve #',i2,' type',i2,'  * half circle',
     +      '(radius, upper/lower (+1/-1) part, midpoint): '/
     +      2x,(f6.2,','),2x,(f6.2,','),2x,2(f6.2,','))
2030  format(' curve #',i2,' type',i2,
     +      '  * linear 3D-function (x1,y1,z1),(x2,y2,z2): '/
     +      2x,3(f6.2,','),2x,3(f6.2,','))
2031  format(' curve #',i2,' type ',i2,
     +      '  * plane function (x0*x + y0*y + z0*z + c0 = 0)'/
     +      2x,4(f6.2,','))
2032  format(' curve #',i2,' type ',i2,
     +      '  * sphere with radius r and midpoint (x,y,z)'/
     +      2x,4(f6.2,','))
3015  format(' curved surface will be generated along curve # ',i2)
3016  format(' curved boundary function ',
     +      'for plane surfaces along curve # ',i2)
3017  format(' curved boundary function ',
     +      'for curved surfaces along curve # ',i2)
3018  format(' linear boundary function ',
     +      'for curved surfaces along curve # ',i2)
3020  format(' b.c. # ',i1,' will be generated along curve #',i2)
3030  format(' Linear loading along curve #',i2,' in ',
     +      i1,'-dir with ordinates ',2(f6.2,','))
      end

      subroutine inp50 ( cxy,n )
c-----------------------------------------------------------------------
c
c.... Purpose: this SR is only needed for curve-typ 50
c....          read n points of data
c
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension cxy(2,*)
      write(iow,1090)
      do i=1,n                          ! read n points
        call dinput ( cxy(1,i),2 )
        write(iow,1100)   cxy(1,i),cxy(2,i)
      enddo
czr200   continue
czr1000  format(2f15.0)
1090  format('curve points')
1100  format(2(g15.5,3x))
      return
      end

      subroutine ctyp50 (x,y,cxy,deriv,n)
c-----------------------------------------------------------------------
c    Druck-Verschiebungskurve bei Silos
C    berechnet Druck y und Ableitung deriv
C   !!!!!!!!  nicht fÅr alle FÑlle anwendbar !!!!!!!!!!!!!!!!!
C
C     muss noch angepasst werden !!!!!!!!!!!!!!!!!!!!!
C
C     danger !!!!!  nicht verwenden  !!!!!!!!!!!!!!!!!!!!
C
c-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension cxy(2,*)
      do 100 i=1,n-1
        if( x.lt.cxy(1,i) .and. x.ge.cxy(1,i+1) ) then
           x1    = cxy(1,i)
           x2    = cxy(1,i+1)
           y1    = cxy(2,i)
           y2    = cxy(2,i+1)
           dx    = x2-x1
           dy    = y2-y1
           deriv = dy/dx
           y     = deriv*(x-x1)+y1
           return
        endif
100   continue
      deriv = 0.0d0
      y     = cxy(2,n)
      return
      end

