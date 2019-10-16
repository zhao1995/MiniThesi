      subroutine cvuser ( x,y,z,ityp,m,nr,diff )
c-----------------------------------------------------------------------
c
c....   Purpose: calculation of user defined curves
C                ityp 1  to 50  : standard curves     ( same in all feap version )
C                ityp 51 to 100 : user defined curves ( only in user version )
c
c                Knebel's curves
c
c      Input:   x    = argument
c               y    = argument 
c               nr   = curve number
c               ityp =
c               m    = m-array
c               nr
c 
c
c      Output:   y    = value at point x   (only for 2-d)
c                z    = value at point x,z (only for 3-d)
c                diff = 0 if x,y,z lies on curve (else diff<>0)
C                      ( only for implicit defined functions)
c-----------------------------------------------------------------------

      USE cdata
      USE iofile
      USE tdata
      implicit double precision (a-h,o-z)
      common /curvedat/ cpar(20,8),nrt1(20,3),ic,nbe,nn3
      dimension m(*)
      data pi /3.14159265359/
c
      if ( ityp.eq.51 ) then   ! ( war nr 10)

c....  - Konstante und Cosinus -
C
      s1 = cpar(nr,1) * pi/180.0
      s2 = cpar(nr,2) * pi/180.0
      a  = cpar(nr,3) 
      b  = cpar(nr,4)
      if(x .le. s1) then 
         y = a
      elseif ( x.gt.s1 .and. x.lt.s2 ) then
         a1 = 2.0* pi/(s2-s1)
         a2 = 0.5*(s1+s2)
         y  = 0.5*b*cos(a1*(x-a2)) + 0.5*b+a
      elseif ( x.ge.s2 .and. x.le.2.0*pi) then
         y = a
      else
        write(*,*)    ' Error in CURVE   Argument X does not fit'
        write(iow,*) ' Error in CURVE   Argument X does not fit'
        stop
      endif
C
      elseif ( ityp.eq.52 ) then   ! ( war nr 11 )
C
c....  - links Parabel, rechts horiz. Gerade -
C
C         y=ax^2 + po         po
C         def. einer Lastfunktion (Knebel) z=Ableitung
        p0 = cpar(nr,1)
        a  = cpar(nr,2)
        if( x .ge. 0.0 ) then
          y = p0
          z = 0.0d0
        else
          y = a*x*x + p0
          z = 2.0d0*a*x
        endif
      return
c
      elseif ( ityp.eq.53 ) then   ! ( war nr 12 )
C
c....  - links lineare Gerade, rechts horiz. Gerade -
C
        p0 = cpar(nr,1)         ! druck
        a  = cpar(nr,2)         ! bettung
        if( x .gt. 0.0d0 ) then
          y = p0
          z = 0.0d0
        else
          y = a*(-x) + p0 
          z = a
        endif
c        if(nn3.ne.0 .and. ttme .eq. 1.0d0) then
c        if(nn3.ne.0 .and. ttme .le. 1.0d0) then
        if(nn3.ne.0 .and. ttme .eq. 0.0d0) then
          y = p0           ! fuellen keine steifigkeit
          z = 0.0d0
        endif
      return
c
      elseif ( ityp.eq.54 ) then   ! ( war nr 13 )
C
c....  - links lineare Gerade, rechts lineare Gerade -
C        beidseitige gleiche Bettung  Z-links = Z-rechts
C
        p0 = cpar(nr,1)
        a  = cpar(nr,2)
        if( x .gt. 0.0d0 ) then
          y = p0 - a*x
          z = a
          xg = p0/tan(a)
          if(x.gt.xg) y=0.0d0
          if(x.gt.xg) z=0.0d0          
        else
          y = a*(-x) + p0
          z = a
        endif
c        if(nn3.ne.0 .and. ttme .eq. 1.0d0) then
c        if(nn3.ne.0 .and. ttme .le. 1.0d0) then
        if(nn3.ne.0 .and. ttme .le. 0.0d0) then
          y = p0
          z = 0.0d0
        endif
      return
c
      elseif ( ityp.eq.55 ) then   ! ( war nr 14 )
C
c....  - links lineare Gerade, rechts lineare Gerade -
C        beidseitige unterschiedliche Bettung
C        Z-links  =  zl
C        Z-rechts =  zr
C
        p0 = cpar(nr,1)
        zl = cpar(nr,2)
        zr = cpar(nr,3)
        if( x .gt. 0.0d0 ) then
          y = p0 - zr*x
          z = zr
          xg = p0/tan(zr)
          if(x.gt.xg) y=0.0d0
          if(x.gt.xg) z=0.0d0           
        else
          y = zl*(-x) + p0
          z = zl
        endif
c        if(nn3.ne.0 .and. ttme .eq. 1.0d0) then
c        if(nn3.ne.0 .and. ttme .le. 1.0d0) then
        if(nn3.ne.0 .and. ttme .le. 0.0d0) then
          y = p0
          z = 0.0d0
        endif
      return
c
      elseif ( ityp.eq.60 ) then   ! ( war nr 20 )
C---------------------------  teillast am vollzylinder
        p0   = cpar(nr,1)             ! basisdruck
        a    = 0.5*(p0-cpar(nr,2))    ! amplitude
        phi1 = cpar(nr,3)             ! anfangswinkel
        phi2 = cpar(nr,4)             ! endwinkel
        yy1  = cpar(nr,5)             ! anfangshîhe
        yy2  = cpar(nr,6)             ! endhîhe
         ic  = cpar(nr,7)             ! Auswahlparameter
C
      phi1 = phi1*pi/180.d0          ! in bogenma·
      phi2 = phi2*pi/180.d0          ! in bogenma·
      dphi = phi2-phi1
      dyy  = yy2-yy1
C
c     Umlaufswinkel  phi
      phi = datan(x/z) 
c
c     Umfangswelle am Vollzylinder : ic=3
      if(ic.eq.3) then
        if(y.gt.yy1.and.y.lt.yy2) then
        yty=y-yy1
        pphi=phi-phi1
        pa=-a*(cos(2.d0*pi*yty/dyy)-1.d0)
        p=p0-pa
        else
        pa=p0
        endif
       goto 123
       else
       endif
C
C     öberprÅfung,ob Punkt im Bereich konstanten Drucks liegt
C     teillast beim vollzylinder nur im ersten quadranten
      if(y.lt.yy1.or.y.gt.yy2.or.phi.lt.phi1.or.phi.gt.phi2) then
         p=p0
      else
      pa = 0.d0
      yty=y-yy1        
      pphi=phi-phi1
      if(ic.eq.0.and. x.ge.0.0.and. z.ge.0.0) then    ! beide Richtungen
        pa=a/2.d0*((cos(2.d0*pi*pphi/dphi)-1.d0)
     +             *(cos(2.d0*pi*yty/dyy)-1.d0))
      elseif(ic.eq.1) then                             ! Umfangsrichtung
        pa=-a*(cos(2.d0*pi*yty/dyy)-1.d0)
      elseif(ic.eq.2.and. x.ge.0.0.and. z.ge.0.0) then ! Laengsrichtung
        pa=-a*(cos(2.d0*pi*pphi/dphi)-1.d0)
      else
      endif
      p=p0-pa
      endif
123   continue
      z=max(p,0.d0)
      return
c
      elseif ( ityp.eq.61 ) then   ! ( war nr 21 )
        p0   = cpar(nr,1)             ! basisdruck
        fact = cpar(nr,2)             ! erhoehungsfaktor nach DIN
        phi1 = cpar(nr,3)             ! anfangswinkel
        phi2 = cpar(nr,4)             ! endwinkel
        yy1  = cpar(nr,5)             ! anfangshîhe
        yy2  = cpar(nr,6)             ! endhîhe
C
        phi1 = phi1*pi/180.d0         ! in bogenma·
        phi2 = phi2*pi/180.d0         ! in bogenma·
        dphi = phi2-phi1
        dyy  = yy2-yy1
C
c     Umlaufswinkel  phi
      phi = datan(x/z) 
c
C     öberprÅfung,ob Punkt im Bereich konstanten Drucks liegt
C     teillast beim vollzylinder nur im ersten quadranten
      if(y.lt.yy1.or.y.gt.yy2.or.phi.lt.phi1.or.phi.gt.phi2) then
        p=p0
      else
        p = p0*fact
      endif
      z=p
      return
c
      else
        write(iow,1000) ityp
        write(*  ,1000) ityp
        stop
      endif
      return
1000  format('***error*** user-curve type ',i3,' in not defined')
      end
