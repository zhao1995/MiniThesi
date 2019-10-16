c-----------------------------------------------------------------------
c
      subroutine colpal(ipal,nnc)
c-----------------------------------------------------------------------
c.... Purpose:  set color palette for ncc  colors
c
c     Input:
c       nnc    - number of colors
c
c     Output:
c       ipal     associated color palette
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      integer*4 ipal(14),ipalnc(14,14)
c.... set basis color table from blue to red for PHIGS/VGA
c.... set color tables for nnc .lt. 14
      data ipalnc
     +           /7, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     +            3,12, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     +            3, 7,12, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     +            3, 6,11,13, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     +            3, 5, 7, 9,13, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     +            3, 5, 7, 9,11,13, 1, 1, 1, 1, 1, 1, 1, 1,
     +            2, 4, 6, 7, 9,11,13, 1, 1, 1, 1, 1, 1, 1,
     +            3, 5, 6, 7, 9,10,11,13, 1, 1, 1, 1, 1, 1,
     +            2, 3, 5, 6, 7, 9,10,12,13, 1, 1, 1, 1, 1,
     +            2, 3, 4, 5, 6, 7, 9,11,12,13, 1, 1, 1, 1,
     +            2, 3, 4, 5, 6, 7, 9,10,11,12,13, 1, 1, 1,
     +            2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13, 1, 1,
     +            1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13, 1,
     +            1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14/
c
c.... define color table
      do i = 1,14
c....   define color number
        ii = ipalnc(i,nnc)
c....   find actual color
        ipal(i) = ii
      end do
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine dispola(x,b,ndm,ndf,numnp,ipola)
c-----------------------------------------------------------------------
c.... Purpose: modify nodal values of displacements cart/cylindrical
c
c     Input:
c       ipola        - number of colors
c       x(ndm,numnp) - ccordinates
c       b(ndf,numnp) - displacements
c
c     Output:
c       b(ndf,numnp) - displacements
c       for ipola = 0   x-y   displacements
c           ipola = 12  r-phi displacements  in 1-2 plane
c           ipola = 13  r-phi displacements  in 1 3 plane
c           ipola = 23  r-phi displacements  in 2-3 plane
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension x(ndm,*),b(ndf,*)
c
c.... transform to polar coordinates
      if(ipola.eq.12) then
        i=1
        k=2
      else if(ipola.eq.13) then
        i=1
        k=3
      else if(ipola.eq.23) then
        i=2
        k=3
      end if
      do n = 1,numnp
        xk = x(k,n)
        xi = x(i,n)
        if(xi.eq.0d0 .and. xk.eq.0.d0) then ! special case sphere at top
          sn = 0.d0
          cs = 0.d0
        else
          phi = datan2(xk,xi)
          sn = dsin(phi)
          cs = dcos(phi)
        end if
        ur = cs*b(i,n) + sn*b(k,n)
        ut =-sn*b(i,n) + cs*b(k,n)
        b(i,n) = ur
        b(k,n) = ut
        if(ndf.ge.6) then ! always for rotations
          ur = cs*b(i+3,n) + sn*b(k+3,n)
          ut =-sn*b(i+3,n) + cs*b(k+3,n)
          b(i+3,n) = ur
          b(k+3,n) = ut
        end if
      end do
      return
      end
c
      subroutine drawmess(text,icolo,ityp)
c-----------------------------------------------------------------------
c.... Purpose:  draw message
c
c     Input:
c       text           - string
c       icolo          - color
c       ityp = -3      - text    wait until click
c       ityp = -2      - text    wait some times
c       ityp = -1 or 0 - text    wait until click error
c       ityp =  1      - graphic wait until click error
c       ityp =  2      - graphic wait some times
c       ityp =  3      - graphic wait until click
c       ityp =  4      - open  graphic
c       ityp =  5      - close graphic
c       ityp =  6      - write always
c
c     Output:          - text
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE fdata
      USE iofile
      USE pdata2
      USE pdata8
      USE ptext
      implicit double precision(a-h,o-z)
      integer*2 icolo
      character*(*) text
c.....compare text, only 2 times same message
      if(text.eq.text1.and.text.eq.text2) then
        text2=text1
        text1=text
        return
      end if

      text2=text1
      text1=text

      if(ior.gt.0.and.pfr) then
          write(iow,*) text
      end if

      if(idev.lt.3) then
c....   HP (always the same)
        if(pfr) write(*,*) text
      end if
c.... WIN-INTEL/SALFORD
      if(idev.eq.3.or.idev.eq.4) then
c       call mess_win(text,ityp)
c       if(pfr) call mess_win(text,ityp)
        if(iplot.gt.0) call mess_win(text,ityp) ! not for batch FE2
      end if

c.... ALWAYS
      if(ityp.eq.6) then
        write(iow,*) text
        if(idev.lt.3) write(*,*) text
        if(idev.eq.3.or.idev.eq.4) call mess_win(text,0)
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine frame(x,ndm,numnp,isw,clip)
c-----------------------------------------------------------------------
c.... Purpose: compute scaling for plot area including clip
c
c     Input:
c       ipola        - number of colors
c       x(ndm,numnp) - ccordinates
c       b(ndf,numnp) - displacements
c       isw          -
c       clip         - clip
c
c     Output:        scale for plot area
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE iofile
      USE iwinio
      USE pclip
      USE pdata1
      USE pdata2
      USE pdata4
      USE pdatap
      USE pltran
      USE ppers
      implicit double precision (a-h,o-z)
c
c.... compute scaling for plot area
c
      logical clip,errv
      dimension x(ndm,*),xmn(3),xmx(3),xmino(3),xmaxo(3),xpl(4),xg(3),
     1          xcp(2,3),xzml(3,2)

      dimension nz(2)

      save xmino,xmaxo
      clip1 = .false.
      if(clip) clip1 = .true.
c.... determine coordinates from window
      if(nzm1.lt.0) then
cww     if(idev.eq.4) call clwopen('CLIP   Parameters',1,iwys-120,700,150,1,2)
        if(idev.eq.4) call clwopen('CLIP   Parameters',1,2)
        call pltzm(xpl)
        if(idev.eq.4) call clwclose(1,2,0)
c
c       Find closest node
        do i = 1,2
          nz(i) = 1
          x1    = xpl(2*i-1)
          y1    = xpl(2*i)
          xm = 1e6
          do n = 1,numnp
c.....      coordinates of node (cart or pers)
            xx1 = x(1,n)
            xx2 = x(2,n)
cww         if(ndm.eq.3) xx3 = x(3,n)
            if(ndm.ge.3) xx3 = x(3,n)
c....       screen coordinates of node
            s1 = scale*(xx1 + xx1 -sx(1)) + s0(1)
            s2 = scale*(xx2 + xx2 -sx(2)) + s0(2)
c....       if isometric recompute values
            if(iso) then
              xmul = 2.*max(dx(1),dx(2))/(dx(1)+dx(2))
              te   = (s0(1) + 0.5*(s1 - s2))*xmul
c....         again screen coordinates of node
              s2   = (0.2885*(s1 + s2))*xmul + scale*(xx3+xx3)
              s1   = te - 0.1
            end if
c....       distance (point-node in screen coordinates)
            ym = (s1-x1)**2+(s2-y1)**2
            if(ym.lt.xm) then
              xm    = ym
              nz(i) = n
            end if
100         continue
          end do
        end do
c.....  cart/pers coordinates of nearest nodes to clip points
        do i = 1,2
          xc(i,1) = x(1,nz(i))
          xc(i,2) = x(2,nz(i))
          if(ndm.ge.3) xc(i,3) = x(3,nz(i))
          xcp(i,1) = xc(i,1)
          xcp(i,2) = xc(i,2)
          if(ndm.ge.3) xcp(i,3) = xc(i,3)
        end do

        nzm1 = 0
c.....  Use macro like zoom,n1,n2
        nzm1 = nz(1)
        nzm2 = nz(2)
        clip = .false.

c.....  clear screen with next action
        iclear = 0
      end if
c.....ww??
      if(ndm.eq.1) then
        dx(2) = 0.
        sx(2) = 0.0
      end if
c
      ii = min(ndm,3)
      ij = min(ndm,2)
c.... define start values for min.and max. coordinate of input nodes (=max/min values!)

      if(isw.eq.1) then
        call pzero(xmin,3)
        call pzero(xmax,3)
        do 107 i = 1,ii
c......   values xmin
          xmin(i) = x(i,1)
          do 106 n = 1,numnp
            xmin(i) = max(xmin(i),x(i,n))
106       continue
c......   values xmax
          xmax(i) = x(i,1)
107     continue
      else
c....   choose old values
        do 108 i = 1,ii
          xmin(i) = xmino(i)
          xmax(i) = xmaxo(i)
108     continue
      end if
c.... find the minimum and maximum coordinate of all input nodes
      do 110 n = 1,numnp
        if(x(1,n).ne. -999.d0) then
          do 109 i = 1,ii
            xmin(i) = min(xmin(i),x(i,n))
            xmax(i) = max(xmax(i),x(i,n))
            if(mod(isw,2).eq.1) then
              xmino(i) = xmin(i)
              xmaxo(i) = xmax(i)
            end if
109       continue
        end if
110   continue
c.... plot region determination
      do 111 i = 1,ij
c....   complete region in window
        xw1 = xmin(i)
        xw2 = xmax(i)
c....   modify the window
        if (clip) then
c....     ... in case of clipping  x_1-x_2 plane
          xw1 = xc(1,i)
          xw2 = xc(2,i)
          if(kpers.eq.1) then
            xw1 = xcp(1,i)
            xw2 = xcp(2,i)
          end if
        else
c....     ... for  nzm1,nzm2  (zoom between nodes)
          if(nzm1.gt.0 .and. nzm1.le.numnp ) then
            xw1 = x(i,nzm1)
          end if
          if(nzm2.gt.0 .and. nzm2.le.numnp ) then
            xw2 = x(i,nzm2)
          end if
c....     ... for  nzm3  (zoom between coordinates)
          if(nzm3.gt.0) then
            call pmove(xzm,xzml,6)
c....       rotation
            call plxtrn(xzml,tra,vr,3,2)
c....       perspective projection
            if(kpers.eq.1) call perspj(xzml,xzml,3,3,2,errv)
            xw1 = xzml(i,1)
            xw2 = xzml(i,2)
          end if
        end if
c....   define window limits
        xmn(i) = min(xw1,xw2)
        xmx(i) = max(xw1,xw2)
        dx(i) = xmx(i) - xmn(i)
        sx(i) = xmx(i) + xmn(i)
111   continue
c.... rescale window
      if(dx(1).gt.dx(2)) then
        xmn(2) = (sx(2) - dx(1))/2.0
        xmx(2) = (sx(2) + dx(1))/2.0
      else
        xmn(1) = (sx(1) - dx(2))/2.0
        xmx(1) = (sx(1) + dx(2))/2.0
      end if
      do 112 i = 1,ij
        xmin(i) = max(xmin(i),xmn(i))
        xmax(i) = min(xmax(i),xmx(i))
112   continue
      scaleg = max(xmax(1)-xmin(1),xmax(2)-xmin(2))
      if(ii.gt.2) scaleg = max(scaleg,xmax(3)-xmin(3))
      do 113 i = 1,ij
        xmin(i) = xmin(i) - scaleg/100.
        xmax(i) = xmax(i) + scaleg/100.
113   continue
c.... default values
      scale  = max(xmax(1)-xmin(1),xmax(2)-xmin(2))
      if(scale.le.0.0d0) scale = scaleg/10
c.... reset values for deformed plotting
      do 114 i = 1,ii
        xcen = xmax(i)+xmin(i)
        xmax(i) = (xcen + 1.1*scale)/2.
        xmin(i) = (xcen - 1.1*scale)/2.
114   continue
      scaleg = 0.4/scale
      scale  = scaleg
      s0(1)  = 0.5
      s0(2)  = 0.5
c
      return
 2001 format(' Body occupies the space with:'/
     1  '                 X',13x,'Y',13x,'Z'/
     2  '   Minimum ',1p3e14.5/'   Maximum ',1p3e14.5/)
 2002 format(' Chosen body space :'/
     1  '                 X',13x,'Y'/
     2  '   Minimum ',1p2e14.5/'   Maximum ',1p2e14.5/)
 2003 format(' choose coordinates x_3 of point 1 and 2 :',$)
      end
c
c-----------------------------------------------------------------------
c
      subroutine maxcor(x,ndm,numnp)
c-----------------------------------------------------------------------
c
c.... Purpose: determine max window coordinates
c
c     Input:
c       x(ndm,numnp) - ccordinates
c
c     Output:        - vmin(3),vmax(3)
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata0
      implicit double precision (a-h,o-z)
      dimension x(ndm,*)
      ii = min(ndm,3)
      do 102 i = 1,ii
        vmin(i) = x(i,1)
        do 100 n = 1,numnp
            vmin(i) = max(vmin(i),x(i,n))
100     continue
        vmax(i) = vmin(i)
        do 101 n = 1,numnp
          if(x(1,n).ne.-999.d0) then
            vmin(i) = min(vmin(i),x(i,n))
          end if
101     continue
102   continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pppcolf(force)
c-----------------------------------------------------------------------
c
c.... Purpose: find color for plot of a force
c
c     Input:
c       force - number of force
c
c     Output:
c       icol  - color number
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
c
      USE forccont
      implicit double precision (a-h,o-z)
c
c.....find color
c     ipali,vc,nc
      icol = 1
      if(force.le. vc(1)) then
        icol=1
        goto 100
      end if
      do i = 2,nc
        v2 = vc(i)
        v1 = vc(i-1)
        if(force.gt.v1.and.force.le.v2) then
          icol = i
          goto 100
      end if
      end do
      if(force.gt. vc(nc)) then
       icol=nc+1
       goto 100
      end if
      if(nc.eq.0) then ! special case 1 color and F >F+eps
       icol=nc+1
       goto 100
      end if

      stop 'no force found'
100   icol1 = ipali(icol)
      call pppcol(16+icol1)
      return
      end
c
      subroutine perspe
c-----------------------------------------------------------------------
c
c.... Purpose: input of perspective parameters
c
c     Input:
c       kpers   - flag for perspective projection (1=perspective)
c       eold(3) - view point   def=[3v_x,2v_y,1.5v_z]_max
c       vold(3) - axial vector def=[0,0,1]
c
c     Output:
c
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
c
      USE iofile
      USE iwinio
      USE pback
      USE pdata0
      USE pdata2
      USE ppers
      USE vpoint
      implicit double precision (a-h,o-z)
      common /vpers/ v(3)
      dimension t(3,3)
c
    5 continue
CWW   if(idev.eq.4) call clwopen('Perspective  Parameters',1,iwys-260,500,300,1,2)
      if(idev.eq.4) call clwopen('Perspective  Parameters',1,2)

      eoldq = dsqrt(dot(eold,eold,3))

      if(eoldq.eq.0.0d0) then
        eold(1) = 3.0d0*vmax(1)
        eold(2) = 2.0d0*vmax(2)
        eold(3) = 1.5d0*vmax(3)
      end if

      if(ior.lt.0) write(*,2001) vmin,vmax,eold

      call dinput(e,3) ! read viewpoint

      eq = dsqrt(dot(e,e,3))

      if(eq.eq.0.0d0) then
        e(1) = eold(1)  ! use old vwp
        e(2) = eold(2)
        e(3) = eold(3)
      else
        eold(1) = e(1)  ! set old vwp
        eold(2) = e(2)
        eold(3) = e(3)
      end if

      tg(1)    = 0.5*(vmax(1) + vmin(1))
      tg(2)    = 0.5*(vmax(2) + vmin(2))
      tg(3)    = 0.5*(vmax(3) + vmin(3))

      tgold(1) = tg(1)
      tgold(2) = tg(2)
      tgold(3) = tg(3)

      if(ior.lt.0) write(*,2002) vold

      call dinput(v,3) ! read axial vector

      vq = dsqrt(dot(v,v,3))

      if(vq.eq.0.0d0) then
        v(1) = vold(1)  ! use old axial v
        v(2) = vold(2)
        v(3) = vold(3)
      else
        vold(1) = v(1)  ! set old axial v
        vold(2) = v(2)
        vold(3) = v(3)
      end if

      ff   = 1.d0
      fold = 1.d0
c
c.....projection matrix t
      q(1,1) = e(1) - tg(1)
      q(2,1) = e(2) - tg(2)
      q(3,1) = e(3) - tg(3)

      enorm = sqrt(q(1,1)*q(1,1)+q(2,1)*q(2,1)+q(3,1)*q(3,1))

      if(enorm.le.0.0d0) go to 901
      do i=1,3
        qq(i) = q(i,1)
      end do
      if(idev.eq.4) call clwclose(1,2,0)
c
      if(ibor.ge.0) call plview(qq,vmin,vmax)
c
      do 10 i=1,3
        q(i,3) = q(i,1) / enorm
   10 continue
      do 30 i=1,3
        do 20 j=1,3
          t(i,j) = - q(i,3) * q(j,3)
   20   continue
        t(i,i) = t(i,i) + 1.d0
   30 continue

c.....find projection of v
      call mult(t,v,q(1,2),3,3,1)
      vnorm = sqrt(q(1,2)*q(1,2)+q(2,2)*q(2,2)+q(3,2)*q(3,2))
      if(vnorm.le.0.0d0) go to 901
      do 40 i=1,3
        q(i,2) = q(i,2) / vnorm
   40 continue

c.....find u
      call vecp(q(1,2),q(1,3),q(1,1))

c.....find lambda
      call tran(q,3)
      call mult(q,t,xlbda,3,3,3)
      call tran(q,3)

c.... save viewpoint for smesh
      vwpt(1) = e(1)
      vwpt(2) = e(2)
      vwpt(3) = e(3)
      return
c
  901 if(ior.lt.0)
     +call drawmess(' Improper view has been specified',1,0)
      call errclr('PERSPE')
      goto 5

 2001 format(' Body occupies the space with:',/,
     1  '                 X',13x,'Y',13x,'Z',/,
     2  '   Minimum ',1p,3e14.5,/,'   Maximum ',1p,3e14.5,/,
     3  ' Enter coordinates of view point   (X,Y,Z).',/,
     4  ' Default/Last: X=',1pe9.2,', Y=',1pe9.2,', Z=',1pe9.2,/,' >',$)
 2002 format(' Enter comps. of vertical vector   (X,Y,Z).',/,
     1  ' Default/Last: X=',1pe9.2,', Y=',1pe9.2,', Z=',1pe9.2,/,' >',$)
      end
c
c-----------------------------------------------------------------------
c
      subroutine mult(a,b,c,nra,nca,ncb)
c-----------------------------------------------------------------------
c     multiply  c_ik = a_ij*b_jk
c-----------------------------------------------------------------------
      Implicit real*8(a-h,o-z)
      implicit integer(i-n)
      dimension a(nra,nca),b(nca,ncb),c(nra,ncb)
      do 200 i=1,nra
        do 200 j=1,ncb
          s=0.d0
          do 100 k=1,nca
            s=s+a(i,k)*b(k,j)
  100     continue
          c(i,j)=s
  200 continue
      return
      end
c
      subroutine tran(w,n)
c-----------------------------------------------------------------------
c.... transpose matrix                                                 |
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      implicit integer(i-n)
      dimension w(n,n)
      do 105 i=1,n
      ip1=i+1
      do 105 j=ip1,n
        temp=w(i,j)
        w(i,j)=w(j,i)
        w(j,i)=temp
  105 continue
      return
      end
c
      subroutine vecp (e1,e2,e3)
c-----------------------------------------------------------------------
c      vector product of two 3-d vectors                               |
c      input     e1(3),e2(3) = vectors to be multiplied                  |
c      output    e3(3)       = vector product ( e3 = e1xe2 )             |
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension e1(3),e2(3),e3(3)
      e3(1) = e1(2)*e2(3) - e1(3)*e2(2)
      e3(2) = e1(3)*e2(1) - e1(1)*e2(3)
      e3(3) = e1(1)*e2(2) - e1(2)*e2(1)
      return
      end
c
      subroutine plview(qq,vmin,vmax)
c-----------------------------------------------------------------------
c
c.... Purpose: plot viewpoint for perspectives
c
c     Input:
c         qq(3)
c       vmin(3)
c       vmax(3)
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension qq(3),vmin(3),vmax(3), bmax(3), view(3)
c
      call pppcol(5)
c.... vectors - 1
      call dplot( 0.9854650d0, 0.8825000d0, 3,0)
      call dplot( 1.1091850d0, 0.8825000d0, 2,0)
      call dplot( 1.1014525d0, 0.8902325d0, 2,0)
      call dplot( 1.1014525d0, 0.8747675d0, 2,0)
      call dplot( 1.1091850d0, 0.8825000d0, 2,0)
c.... vectors - 1
      call dplot( 1.1401150d0, 0.8825000d0, 3,0)
      call dplot( 1.2638350d0, 0.8825000d0, 2,0)
      call dplot( 1.2561025d0, 0.8902325d0, 2,0)
      call dplot( 1.2561025d0, 0.8747675d0, 2,0)
      call dplot( 1.2638350d0, 0.8825000d0, 2,0)
c.... vectors - 2
      call dplot( 1.0473250d0, 0.8306400d0, 3,0)
      call dplot( 1.0473250d0, 0.9443600d0, 2,0)
      call dplot( 1.0395925d0, 0.9366275d0, 2,0)
      call dplot( 1.0550575d0, 0.9366275d0, 2,0)
      call dplot( 1.0473250d0, 0.9443600d0, 2,0)
c.... vectors - 3
      call dplot( 1.2012750d0, 0.8306400d0, 3,0)
      call dplot( 1.2012750d0, 0.9443600d0, 2,0)
      call dplot( 1.1935425d0, 0.9366275d0, 2,0)
      call dplot( 1.2090075d0, 0.9366275d0, 2,0)
      call dplot( 1.2012750d0, 0.9443600d0, 2,0)
c... label = 1
c     call dplot( 1.0669175d0, 0.8925000d0, 3,0)
      call dplot( 1.06d0      , 0.8925000d0, 3,0)
      call plabl(1)
c     call dplot( 1.2269175d0, 0.8925000d0, 3,0)
      call dplot( 1.216d0    , 0.8925000d0, 3,0)
      call plabl(1)
c... label = 2
      call dplot( 1.0095925d0, 0.9471800d0, 3,0)
      call plabl(2)
c... label = 3
      call dplot( 1.1635425d0, 0.9471800d0, 3,0)
      call plabl(3)
c...
      qmax    = max(abs(qq(1)) , abs(qq(2)) , abs(qq(3)) )
      bmax(1) = (vmax(1) - vmin(1))/2.0d0
      bmax(2) = (vmax(2) - vmin(2))/2.0d0
      bmax(3) = (vmax(3) - vmin(3))/2.0d0
      absmax  = max(qmax , abs(bmax(1)) , abs(bmax(2)) , abs(bmax(3)) )
      rnorm   = 0.06186d0/absmax
      do 100 i = 1,3
        view(i) = rnorm*qq(i)
        bmax(i) = rnorm*bmax(i)
100   continue

c...  plot box of body limits - 12 -plane
      xx = 1.047325d0 - bmax(1)
      yy = 0.882500d0 - bmax(2)
      call ppbox( xx , yy, 2.d0*bmax(1), 2.d0*bmax(2), 3)

c...  plot box of body limits - 13 -plane
      xx = 1.201275d0 - bmax(1)
      call ppbox( xx , yy, 2.d0*bmax(1), 2.d0*bmax(2), 3)

c...  plot the eyes - 12 -plane
      call ppeye( view(1) + 1.047325d0 , view (2) + 0.8825000d0 )
c...  plot the eyes - 13 -plane
      call ppeye( view(1) + 1.201275d0 , view (3) + 0.8825000d0 )
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine ppeye(v1, v2)
c-----------------------------------------------------------------------
c
c.... Purpose: plot view point on screen
c
c     Input:
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      real*8 v1, v2
      call dplot( v1              , v2 + 0.0077325d0 , 3,0)
      call dplot( v1 - 0.015465d0 , v2               , 2,0)
      call dplot( v1              , v2 - 0.0077325d0 , 2,0)
      call dplot( v1 + 0.015465d0 , v2               , 2,0)
      call dplot( v1                , v2 + 0.0077325d0 , 2,0)
      call dplot( v1 - 0.0077325d0, v2               , 2,0)
      call dplot( v1                , v2 - 0.0077325d0 , 2,0)
      call dplot( v1 + 0.0077325d0, v2               , 2,0)
      call dplot( v1                , v2 + 0.0077325d0 , 2,0)
      end
c
c-----------------------------------------------------------------------
c
      subroutine perspj(x,xpr,ndm1,ndm2,numnp,errv)
c-----------------------------------------------------------------------
c
c.... Purpose: perspective projection of coordinates
c
c     Input:
c        x(ndm1,numnp) - global coordinates
c        numnp         - total number of nodal points
c
c     Output:
c       xpr(ndm2,numnp)- perspective projected coordinates
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE iofile
      USE ppers
      implicit double precision (a-h,o-z)
      logical errv
      dimension t1(3),x(ndm1,numnp),xpr(ndm2,numnp)
c
c.....loop over the data points and find projection
      do 500 n=1,numnp
        if(x(1,n).ne. -999.d0) then
          do 48 i=1,3
            t1(i) = x(i,n) - e(i)
   48     continue
          call mult(xlbda,t1,xpr(1,n),3,3,1)
c
          alpha = - ff*enorm/(t1(1)*q(1,3)+t1(2)*q(2,3)+t1(3)*q(3,3))
          if(alpha.lt.0.d0) then
            errv = .true.
            if(ior.lt.0) then
              write(*,2000)
              return
            else if(ior.ge.0) then
              write(iow,2000)
              stop 'SR perspj 2000'
            end if
          else
            errv = .false.
          end if
          do 50 i=1,3
            xpr(i,n) = alpha * xpr(i,n)
   50     continue
        end if
  500 continue
      return
2000  format(/,1x,' Point too close, choose another one!',/)
      end
c
c-----------------------------------------------------------------------
c
      subroutine pfeap(xl,yl,siz,icol)
c-----------------------------------------------------------------------
c
c.... Purpose: logo FEAP
c
c     Input:
c        xl,yl  - position
c        siz    - size
c        icol   - color
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata2
      implicit double precision (a-h,o-z)
      integer ifx(11),ify(11)
      integer iex(13),iey(13)
      integer iax(12),iay(12)
      integer ipx(15),ipy(15)
      integer ip1x(5),ip1y(5)
c     integer i9x(17),i9y(17)
c     integer i5x(16),i5y(16)
c     integer i2x(16),i2y(16)
c     integer i0x(15),i0y(15)

      integer ixl(8)


      data ifx/ 0,10,45,43,18,16,26,24,14,10, 0/
      data ify/ 0,50,50,40,40,30,30,20,20, 0, 0/

      data iex/ 0,10,45,43,18,16,26,24,14,12,37,35, 0/
      data iey/ 0,50,50,40,40,30,30,20,20,10,10, 0, 0/

      data iax/ 0,20,30,40,30,28,14,18,26,23,10, 0/
      data iay/ 0,50,50, 0, 0,10,10,20,20,33, 0, 0/

      data ipx/ 0,10,30,40,37,29,13,15,24,26,28,26,18,10,0/
      data ipy/ 0,50,50,40,24,15,15,25,25,28,38,40,40, 0,0/

      data ip1x/ 0,10,20,10, 0/
      data ip1y/ 0,50,50, 0, 0/

c     data i9x / 0,25,34,28,10, 4, 0, 5,19,20,10,13,23,19,17, 2, 0/
c     data i9y / 0, 0,45,50,50,45,25,20,20,27,27,40,40,20,10,10, 0/

c     data i5x / 0,20,25,29,22,16,18,33,35,10, 5,10,19,17, 2, 0/
c     data i5y / 0, 0, 5,23,30,30,40,40,50,50,25,20,20,10,10, 0/

c     data i2x / 0, 2,28,29,18,17, 7, 9,15,35,43,42,16,37,35, 0/
c     data i2y / 0,10,35,40,40,32,32,45,50,50,40,32,10,10, 0, 0/

c     data i0x / 2, 8,20,35,43,37,25,10,12,27,33,31,18,10, 2/
c     data i0y /10,40,50,50,40,10, 0, 0,10,10,39,40,40, 0,10/

      data ixl/25,65,105,150,195,207,0,0/      ! FEAP ||
c     data ixl/40,80,120,165,0,0,0,0/          ! FEAP
c     data ixl/25,65,105,150,195,218,0,0/      ! FEAP 95
c     data ixl/35,75,115,160,187,197,207,217/  ! FEAP 2000

      if(ibor.lt.0) return
      size = 200.0/siz
      call pppcol(icol)

c.... plot FEAP
c     f
      dixl = xl*size + ixl(1)
      call dplot((ifx(1)+dixl)/size,ify(1)/size+yl,ipgl,0)
      do i = 2,11
        call dplot((ifx(i)+dixl)/size,ify(i)/size+yl,2,0)
      end do
      if(ipgl.eq.1) call clpan

c     e
      dixl = xl*size + ixl(2)
      call dplot((iex(1)+dixl)/size,iey(1)/size+yl,ipgl,0)
      do i = 2,13
        call dplot((iex(i)+dixl)/size,iey(i)/size+yl,2,0)
      end do
      if(ipgl.eq.1) call clpan

c     a
      dixl = xl*size + ixl(3)
      call dplot((iax(1)+dixl)/size,iay(1)/size+yl,ipgl,0)
      do i = 2,12
        call dplot((iax(i)+dixl)/size,iay(i)/size+yl,2,0)
      end do
      if(ipgl.eq.1) call clpan

c     p
      dixl = xl*size + ixl(4)
      call dplot((ipx(1)+dixl)/size,ipy(1)/size+yl,ipgl,0)
      do i = 2,15
        call dplot((ipx(i)+dixl)/size,ipy(i)/size+yl,2,0)
      end do
      if(ipgl.eq.1) call clpan

c.... plot ||
      fac = 0.8
c     |1
      dixl = xl*size + ixl(5)
      call dplot((fac*ip1x(1)+dixl)/size,fac*ip1y(1)/size+yl,ipgl,0)
      do i = 2,5
        call dplot((fac*ip1x(i)+dixl)/size,fac*ip1y(i)/size+yl,2,0)
      end do
      if(ipgl.eq.1) call clpan

c     |2
      dixl = xl*size + ixl(6)
      call dplot((fac*ip1x(1)+dixl)/size,fac*ip1y(1)/size+yl,ipgl,0)
      do i = 2,5
        call dplot((fac*ip1x(i)+dixl)/size,fac*ip1y(i)/size+yl,2,0)
      end do
      if(ipgl.eq.1) call clpan

c.... plot 95
c     fac = 0.8

c     9
c     dixl = xl*size + ixl(5)
c     call dplot((fac*i9x(1)+dixl)/size,fac*i9y(1)/size+yl,ipgl,0)
c     do i = 2,17
c       call dplot((fac*i9x(i)+dixl)/size,fac*i9y(i)/size+yl,2,0)
c     end do
c     if(ipgl.eq.1) call clpan
c
c     5
c     dixl = xl*size + ixl(6)
c     call dplot((fac*i5x(1)+dixl)/size,fac*i5y(1)/size+yl,ipgl,0)
c     do i = 2,16
c       call dplot((fac*i5x(i)+dixl)/size,fac*i5y(i)/size+yl,2,0)
c     end do
c     if(ipgl.eq.1) call clpan
c
c.... plot 2000
c     fac= 0.3
c     2
c     dixl = xl*size + ixl(5)
c     call dplot((fac*i2x(1)+dixl)/size,fac*i2y(1)/size+yl,ipgl,0)
c     do i = 2,16
c       call dplot((fac*i2x(i)+dixl)/size,fac*i2y(i)/size+yl,2,0)
c     enddo
c     if(ipgl.eq.1) call clpan
c     0 3times
c     do k =1,3
c       dixl = xl*size + ixl(5+k)
c       call dplot((fac*i0x(1)+dixl)/size,fac*i0y(1)/size+yl,ipgl,0)
c       do i = 2,15
c         call dplot((fac*i0x(i)+dixl)/size,fac*i0y(i)/size+yl,2,0)
c       end do
c       if(ipgl.eq.1) call clpan
c     enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine plabl(n)
c-----------------------------------------------------------------------
c
c.... Purpose: Plot label
c
c     Input:
c        n     - string e.g. number
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE hpgl1
      USE pdata2
      USE pdatap
      USE pftn77
      USE plotter
      implicit double precision (a-h,o-z)
      character yy*7
      n1=iabs(n)
      if    (n1.gt. 9999) then
        write(yy,'(i6)') n
      else if(n1.gt. 999) then
        write(yy,'(i5)') n
      else if(n1.gt.  99) then
        write(yy,'(i4)') n
      else if(n1.gt.  9) then
        write(yy,'(i3)') n
      else
        write(yy,'(i2)') n
      end if
c
      if(idev.eq.1)  call gptx2(xp(3),5,yy)
      if(idev.eq.2)  call gtx(xpp(2),ypp(2),yy)
      if(idev.eq.3)  call draw_text(yy,ixa,iya,icc)
      if(idev.eq.4)  call draw_text(yy,ixa,iya,icc)
      if(nexte.gt.0.and.iprin.eq.1) then
        call hptext (xxp(2),yyp(2),yy,ihpgl)
      end if
      return
      end
c
      subroutine plbord(icl)
c-----------------------------------------------------------------------
c
c.... Purpose: Plot border around graphic window
c
c     Input:
c        icl     - color
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata2
      integer icl, icol
      if(ibor.lt.0) return
      icol = icl
      call pppcol(icol)
      call ppbox(0.0003d0, 0.0003d0, 1.279d0, 0.9697d0, 3)
c.... place box lines on figure
      call dplot(0.9700d0, 0.0003d0, 3,0)
      call dplot(0.9700d0, 0.9700d0, 2,0)
      call dplot(0.9700d0, 0.8250d0, 3,0)
      call dplot(1.2793d0, 0.8250d0, 2,0)
      call dplot(0.9700d0, 0.0960d0, 3,0)
      call dplot(1.2793d0, 0.0960d0, 2,0)
c.... place lines for perspective view
      call dplot(1.12465d0, 0.8250d0, 3,0)
      call dplot(1.12465d0, 0.9700d0, 2,0)
      end
c
c-----------------------------------------------------------------------
c
      subroutine pleigvt(omega,k1)
c-----------------------------------------------------------------------
c
c.... Purpose: Write Eigenvalue and No. of Eigenvalue on screen
c
c     Input:
c        omega   - valu
c        k1      - number of EV
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pftn77
      implicit double precision(a-h,o-z)
      character yy*16,yyy*15
      write(yyy,'(a12,i3)') 'Eigenvector ',k1
      x = 1.
      y = 0.76
      call drawtxt(1,x,y,1,1,15,yyy)
      write(yyy,'(a15)')  'Eigenvalue     '
      y = 0.73
      call drawtxt(1,x,y,1,1,15,yyy)
      write(yy,'(a4,e12.5)') 'w = ',omega
      y = 0.70
      call drawtxt(1,x,y,1,1,16,yy)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine plface(fe,id,ix,x,ndm,nen1,numnp,numel,iln,ct,nface,
     1                  icol)
c-----------------------------------------------------------------------
c
c.... Purpose: compute location of faces which face forward
c              for 8 node - brick
c
c     Input:
c
c     output:
c     il(4,i): node numbers of sides i of brick
c          i = 1 bottom         4---3    8---7
c          i = 2 front          |   |    |   |
c          i = 3 right          1---2    5---6
c          i = 4 back
c          i = 5 left           bottom     top
c          i = 6 top
c
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical fe(numel),flag
      integer   ix(nen1,numel), id(numnp,*), il(4,6)
      dimension x(ndm,numnp)
      data il/1,4,3,2, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 5,6,7,8/ ! ww numbering bottom?
      call pzerol(fe,.false.,numel)
      call pzeroi(id,10*numnp)

c.... compute location of boundary faces by eliminating repeated faces
      icm= 1
      do 140 n = 1,numel
        ma = ix(nen1,n)
        if(iplma(ma).eq.0)  goto 140
c
        flag = .false.

c....   find shell elements (look if node numbers at local nodes>4 exist)
        if(ix(il(4,2),n).eq.0) then
          fe(n) = .true.
          go to 140
        end if

        do 130 i = 1,6
          mix = min(ix(il(1,i),n),ix(il(2,i),n),
     1              ix(il(3,i),n),ix(il(4,i),n))
          do 100 j = 1,4
            if(mix.eq.ix(il(j,i),n)) go to 110
100       continue
110       k = mod(j+1,4) + 1
          k = ix(il(k,i),n)
c....     add to list if new
          ic = 1
120       if(id(mix,ic).eq.0) then
            id(mix,ic) = k
            flag = .true.
          else if(abs(id(mix,ic)).eq.k) then
            id(mix,ic) = -k
            flag = .false.
          else
            ic = ic + 1
            icm= ic
            if(icm.le.10) go to 120
          end if
130     continue

        fe(n) = flag
140   continue
c.... plot the faces which are unmatched
      if(icm.gt.10) then
        write(*,2000) icm
      end if

      nface = 0
      do 240 n = 1,numel
        ma = ix(nen1,n)
        if(iplma(ma).eq.0) goto 240

c....   find shell elements (look if node numbers at local nodes>4 exist)
        if(ix(il(4,2),n).eq.0) then
           nface = nface+1
           go to 240
        end if

c
        if(fe(n)) then
          fe(n) = .false.

          do 230 i = 1,6
            mix = min(ix(il(1,i),n),ix(il(2,i),n),
     1                ix(il(3,i),n),ix(il(4,i),n))
            do 200 j = 1,4
              if(mix.eq.ix(il(j,i),n)) go to 210
200         continue
210         k = mod(j+1,4) + 1
            k = ix(il(k,i),n)
c
            ic = 1
220         if(abs(id(mix,ic)).eq.k) then
              if(id(mix,ic).eq.k) then
c....           face is to be plotted if visible = plot mesh
                call pfacev(il(1,i),ix(1,n),x,ndm,iln,ct,fe(n),nface,
     1                      icol)
              end if
            else if(id(mix,ic).ne.0 .and. ic.le.9) then
              ic = ic + 1
              go to 220
            end if

230       continue
        end if
240   continue
c.... set line type to original
      call plline(iln)
      return
c
2000  format(' Plot failure: too many faces connected to nodes',/,
     + 'max = 10, actual: ',i5)
      end
c
c-----------------------------------------------------------------------
c
      subroutine pfacev(il,ix,x,ndm,iln,ct,fe,nface,icol)
c-----------------------------------------------------------------------
c
c.... Purpose: plot the visible face for 8 node - brick
c
c     Input:
c
c     output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical fe, visbl
      integer il(4),ix(8)
      real*8  x(ndm,*),xl(3,4)
      save ilnd
      data ilnd/1/
      do 110 j = 1,4
        i = ix(il(j))
        do 100 k = 1,3
          xl(k,j) = x(k,i)
100     continue
110   continue
      if(visbl(xl)) then
c....   plot the visible lines
        nface = nface + 1
        fe    = .true.
        call plline(ilnd)
c....      plot the face
cww        call plotl(xl(1,4),xl(2,4),xl(3,4),3)
cww        do 120 j = 1,4
cww          call plotl(xl(1,j),xl(2,j),xl(3,j),2)
cww120     continue
      else if(ct.ne.0.0d0) then
c....   plot the non-visible lines, only if v2.ne.0
        call plline(iln)
        call pppcol(8)
        call plotl(xl(1,4),xl(2,4),xl(3,4),3)
        do 130 j = 1,4
          call plotl(xl(1,j),xl(2,j),xl(3,j),2)
130     continue
        call pppcol(icol)
      end if
      end
c
c-----------------------------------------------------------------------
c
      logical function visbl(xl)
c-----------------------------------------------------------------------
c
c.... Purpose: Find a visible face for 8 node - brick
c
c     Input:
c
c     output:
c        visbl
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pltran
      USE ppers
      USE vpoint
      implicit real*8(a-h,o-z)
      real*8  xl(3,4),v1(3),v2(3),v3(3),vp(3)
      if(kpers.ne.0) then
        vp(1)=e(1)
        vp(2)=e(2)
        vp(3)=e(3)
      else ! mainly rot
        vp(1)=vwpt(1)
        vp(2)=vwpt(2)
        vp(3)=vwpt(3)
        call plxtrn(vp,tra,vr,3,1)
      end if
      do 100 k = 1,3
          v1(k) = xl(k,3) - xl(k,1)
          v2(k) = xl(k,4) - xl(k,2)
100   continue
      call vecp(v1,v2,v3)
      do 110 k = 1,3
          v1(k) = vp(k) - (xl(k,1)+xl(k,2)+xl(k,3)+xl(k,4))/4.d0
110   continue
      visbl = dot(v3,v1,3).gt.0.0d0
      end
c
c-----------------------------------------------------------------------
c
      subroutine plfacx(fe,id,ix,ixf,x,ndm,nen1,numnp,numel)
c-----------------------------------------------------------------------
c
c.... Purpose: compute location of faces which face forward
c     set field of nodes ixf at visible faces -> ix-field for hide!
c
c     Input:
c     ix(nen1,numel)  nodes on element 
c     output:
c     ixf(5,nnemax)   node on visible faces, faces treated as elements
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical fe(numel)
      integer   ix(nen1,numel), ixf(5,*), id(numnp,*), il(4,6)
      dimension x(ndm,numnp)
      data il/1,4,3,2, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 5,6,7,8/
      nf = 1

      do 240 n = 1,numel
        ma = ix(nen1,n)
        if(iplma(ma).eq.0) goto 240


c....   find shell elements (look if node numbers at local nodes>4 exist)
        if(ix(il(4,2),n).eq.0) then
c...      set node numbers of shell element
          do i = 1,4
            ixf(i,nf) = ix(i,n)
          end do
c....     set material number of shell element
          ixf(5,nf) = ix(nen1,n)
c
          nf=nf+1
          go to 240
        end if

        if(fe(n)) then
          do 230 i = 1,6
            mix = min(ix(il(1,i),n),ix(il(2,i),n),
     1                ix(il(3,i),n),ix(il(4,i),n))
            do 200 j = 1,4
              if(mix.eq.ix(il(j,i),n)) go to 210
200         continue
210         k = mod(j+1,4) + 1
            k = ix(il(k,i),n)
c
            ic = 1
220         if(abs(id(mix,ic)).eq.k) then
             if(id(mix,ic).eq.k) then
c....             face is to be plotted if visible
                  call pfacex(il(1,i),ix(1,n),ixf(1,nf),x,ndm,nen1,nf)
             end if
            else if(id(mix,ic).ne.0 .and. ic.le.9) then
             ic = ic + 1
             go to 220
            end if
230       continue
         end if
240   continue
      end
c
c-----------------------------------------------------------------------
c
      subroutine pfacex(il,ix,ixf,x,ndm,nen1,nf)
c-----------------------------------------------------------------------
c
c.... Purpose: plot the face
c
c     Input:
c
c     output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical visbl
      integer il(*),ix(*),ixf(*)
      real*8  x(ndm,*),xl(3,4)
      do 110 j = 1,4
        i = ix(il(j))
        do 100 k = 1,3
          xl(k,j) = x(k,i)
100     continue
110   continue
      if(visbl(xl)) then
c....   set nodes for the face
        do 120 j = 1,4
          ixf(j) = ix(il(j))
120     continue
c....   set a material number
        ixf(5) = ix(nen1)
        nf = nf + 1
      end if
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine pline(x,ie,ix,ic,numnp,numel,ndm,nen1,nen,nie,ct,isw)
c-----------------------------------------------------------------------
c
c.... Purpose: plot FE-mesh in optimized way
c
c     Input:
c       x(ndm,numnp)   - coordinates
c       ie(nie,numat)  - material
c       ix(nen1,numel) - nodes on elements
c       ic(numnp,*)    -
c       numnp          - max. number of nodes
c       numel          - max. number of elements
c       ndm            - dimesnion of problem
c       nen1           - nen+4
c       nen            -  node per element
c       nie            -
c       ct
c       isw
c
c     output:
c       FE-mesh
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE isogeo
      USE psize
      implicit double precision (a-h,o-z)
      logical ifl,iend,zoom,isw
      logical isNurbs ! isogeo
      dimension x(ndm,*),ie(nie,*),ix(nen1,*),ic(numnp,*),jplt(40),
     1          xl(3,40)
cwd   added include isogeo.h an additional xl array
      dimension xl_IGA(ndm,nen) 
      dimension dummyarray(nen)
c.... maximum connections to any node = nc
      data nc/10/
c.... initialize connection array
      do 100 i = 1,numnp
      do 100 j = 1,nc
        ic(i,j) = 0
100   continue
c.... loop through elements to set up list
      do 160 n = 1,numel
cwd   ploting isogeometric elements
        if (isNurbs()) then
c         creating local geometry array, needed for deformed mesh
          do 1008 i = 1,nen
            ii = ix(i,n)
            if(ii.gt.0) then
              do j = 1,ndm
                xl_IGA(j,i) = x(j,ii)
              end do
            end if
1008      continue
c....     get current patch and corresponding orders
          NURp = ngetNURnmpq(n,3,AInipa,AInmpq)
          NURq = ngetNURnmpq(n,4,AInipa,AInmpq)
c....     compute NURBS coordinates of recent element
          ni = ngetNURni(n,AIninc,AInien,AInipa,1)
          nj = ngetNURni(n,AIninc,AInien,AInipa,2)
          if (rNURknv1(ni).eq.rNURknv1(ni+1).or.
     +    rNURknv2(nj).eq.rNURknv2(nj+1)) goto 160
          do i = ni, ni+1
                Xi1 = rNURknv1(i)
                Xi2 = rNURknv2(nj)
c                call physicalCoordinates(Xi1,Xi2,ni,nj,x1,x2,x3,ndm) !kann nur unverformt rechnen
                call physCoor(Xi1,Xi2,ni,nj,x1,x2,x3,ndm,xl_IGA,
     +               NURp,NURq,AInkv1,AInkv2,3,dummyarray,dummyreal)
                call plotl(x1,x2,x3,3)
                Xi2 = rNURknv2(nj+1)
c                call physicalCoordinates(Xi1,Xi2,ni,nj,x1,x2,x3,ndm)! kann nur unverformt rechnen
                call physCoor(Xi1,Xi2,ni,nj,x1,x2,x3,ndm,xl_IGA,
     +               NURp,NURq,AInkv1,AInkv2,3,dummyarray,dummyreal)
                call plotl(x1,x2,x3,2)
          end do
          do j = nj, nj+1
                Xi1 = rNURknv1(ni)
                Xi2 = rNURknv2(j)
c                call physicalCoordinates(Xi1,Xi2,ni,nj,x1,x2,x3,ndm)  !  kann nur unverformt rechnen
                call physCoor(Xi1,Xi2,ni,nj,x1,x2,x3,ndm,xl_IGA,
     +               NURp,NURq,AInkv1,AInkv2,3,dummyarray,dummyreal)
                call plotl(x1,x2,x3,3)
                Xi1 = rNURknv1(ni+1)
c               call physicalCoordinates(Xi1,Xi2,ni,nj,x1,x2,x3,ndm)   !kann nur unverformt rechnen
                call physCoor(Xi1,Xi2,ni,nj,x1,x2,x3,ndm,xl_IGA,
     +               NURp,NURq,AInkv1,AInkv2,3,dummyarray,dummyreal)
                call plotl(x1,x2,x3,2)
          end do
          goto 160
        end if
cwd     end isogeometric elements
        ii = ix(nen1,n)
        if(iplma(ii).eq.0)  goto 160  ! only for matn
        jj = abs(ct)
        if(jj.eq.0) go to 110
        if(ii.ne.jj) go to 160
        if(ct.lt.0.0d0) call pppcol(jj)
110     call pltord(ie(nie-1,ii), iju,jplt)
c.... check window
        ij = 0
        do 120 i = 1,iju
        j = jplt(i)
        if(j.gt.0.and.j.le.nen) then
          jj = abs(ix(j,n))
          if(jj.gt.0) then
            ij       = ij + 1
            xl(1,ij) = x(1,jj)
            xl(2,ij) = x(2,jj)
            if(ndm.eq.3) then
              xl(3,ij) = x(3,jj)
            else
              xl(3,ij) = 0.0
            end if
          end if
        end if
120     continue
        if ( zoom (xl,3,ij) ) then
c.... look up element nodes
          ii = abs(ix(jplt(1),n))
          do 150 ij = 2,iju
            j = jplt(ij)
c           if((j.le.nen).and.(j.gt.0).and.(ix(j,n).ne.0)) then ! orginal
            if( (j.le.nen).and.(j.gt.0) ) then
           if ( ix(j,n).ne.0 ) then
              jj = abs(ix(j,n))
              if(jj.ne.ii) then
              n1 = min(ii,jj)
              n2 = max(ii,jj)
              do 130 k = 1,nc
                if(ic(n1,k).eq.0) then
                  ic(n1,k) = n2
                  go to 140
                else if(ic(n1,k).eq.n2) then
                  ic(n1,k) = -n2
                  go to 140
                end if
130           continue
140           ii = jj
              end if
             end if
            end if
150       continue
        end if
160   continue
c.... change signs to permit mesh plot
      if(isw) then
        do 250 n = 1,numnp
        do 250 i = 1,nc
          ic(n,i) = abs(ic(n,i))
250     continue
      end if
c.... plot outline of part with continuous lines
      x3 = 0.0
      do 304 ni = 1,numnp
        iend = .true.
        do 303 n = 1,numnp
          ifl = .true.
          n1  =  n
300       do 301 i = 1,nc
            if(ic(n1,i)) 301,303,302
301       continue
          go to 303
302       iend = .false.
          if(ndm.ge.3) x3 = x(3,n1)
          if(ifl) then
            call plotl(x(1,n1),x(2,n1),x3,3)
            ifl = .false.
          end if
          n2       =  ic(n1,i)
          ic(n1,i) = -n2
          if(ndm.ge.3) x3 = x(3,n2)
          call plotl(x(1,n2),x(2,n2),x3,2)
          n1 = n2
          go to 300
303     continue
        if(iend) go to 305
304   continue
305   return
      end
c
c-----------------------------------------------------------------------
c
      subroutine plfleg(idev,name,ival,nval,xmax)
c-----------------------------------------------------------------------
c
c.... Purpose: Plot force etc maximal values
c
c     Input:
c
c     output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      character yy*23, name*20
      dimension xmax(2,3)
      x = 1.
      if(idev.eq.4) x = 0.97
      y = 0.76
      write(yy,'(a20,a1)') name,'s'
      call drawtxt(1,x,y,1,1,21,yy)

c
      do i = 1,nval
        ii = ival-1+i
        y = 0.70-0.03*(i-1)*3
        write(yy,'(a20,a1,i1)') name,'_',ii
        call drawtxt(1,x,y,1,1,22,yy)
c
        y = 0.70-0.03*(3*i-2)
        write(yy,'(a6,e12.5)') ' Min ',xmax(1,i)
        call drawtxt(1,x,y,1,1,18,yy)
c
        y = 0.70-0.03*(3*i-1)
        write(yy,'(a6,e12.5)') ' Max ',xmax(2,i)
        call drawtxt(1,x,y,1,1,18,yy)
      end do
      end
c
c-----------------------------------------------------------------------
c
      subroutine plforc(idev)
c-----------------------------------------------------------------------
c
c.... Purpose: Plot text forces
c
c     Input:
c
c     output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE fornam
      USE pdata10
      USE pftn77
      implicit double precision(a-h,o-z)
      character yz*15,yy*19
      x = 1.
      y = 0.76
      if(forsus(nfp).eq.'               ') then
        call drawtxt(1,x,y,1,1,15,' Force / Moment')
      else
        call drawtxt(1,x,y,1,1,15,forsus(nfp))
      end if
      if(klayf.eq.0) klayf = 1
      write(yz,'(a12,i3)') '  Layer No. ',klayf
      x = 1.
      y = 0.72
      call drawtxt(1,x,y,1,1,16,yz)
      write(yy,'(a2,i2,a5,e10.4)') ' F',nfp,' Max ',xmaxf
      if(idev.eq.4) x = 0.97
      y = 0.64
      call drawtxt(1,x,y,1,1,19,yy)
      write(yy,'(a2,i2,a5,e10.4)') ' F',nfp,' Min ',xminf
      y = 0.61
      call drawtxt(1,x,y,1,1,19,yy)
      end
c
c-----------------------------------------------------------------------
c
      subroutine plogo(icol)
c-----------------------------------------------------------------------
c
c.... Purpose: Plot Logo in color icol
c
c     Input:
c         icol - color number
c
c     output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
        call pfeap(0.98d0,0.02d0,0.240d0,icol)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine plot2d(ie,ix,idis,b,x,xl,nie,ndm,nen,nen1,numel,n1,n2,
     +                  n3)
c-----------------------------------------------------------------------
c
c.... Purpose: plot two dimensional mesh
c
c     Input:
c
c     output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE eldata
      USE hdatam
      implicit double precision (a-h,o-z)
      logical fa,tr
      data fa,tr/.false.,.true./
      dimension ie(nie,*),ix(nen1,*),x(ndm,*),xl(ndm,*),b(*),idis(numel)
c.... initialize plot, draw border and add header
c...  loop over elements to draw mesh
      if(n3.eq.0) then
      if(n2.ne.0) call pppcol(n2)
        do 105 nn = 1,numel
          n = idis(nn) ! sort from hidden line
          ma = ix(nen1,n)
          if(n1.ne.0.and.ma.ne.n1) goto 105
          if(iplma(ma).eq.0)       goto 105
          if(n2.eq.0) call pppcol(ma+1)
          do 104 i = 1,nen
            ii = abs(ix(i,n))
            if(ii) 100,100,102
100         do 101 j = 1,ndm
             xl(j,i) = 0.0
101         continue
          go to 104
102       nel = i
            do 103 j = 1,ndm
              xl(j,i) = x(j,ii)
103         continue
104       continue
          call plot9(ie(nie-1,ma),ix(1,n),xl,ndm,nel,1)
c....     plot elements with only one node (stiffness)
          xsum = 0.d0
          do  i = 2,nen
            ii = abs(ix(i,n))
            if(ii.gt.0) then
              do j = 1,ndm
                xsum = xsum + x(j,ii)**2
              end do
            end if
          end do
          if(xsum.eq.0.d0) then
            ii = abs(ix(1,n))
            call pltel1n(x,ndm,ii)
          end if
105     continue
      else
        hflgu  = .false.
        h3flgu = .false.
        call formfe(b,x,x,x,fa,fa,fa,fa,13,1,numel,1)
      end if
      return
      end
c
      subroutine plot9(iel,ix,x,ndm,nel,ipc)
c-----------------------------------------------------------------------
c
c.... Purpose: plot single element
c
c     Input:
c
c     Output:
c
c     Comment
c     ipc necessary from pltpele: do not fill
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata7
      implicit double precision(a-h,o-z)
      logical zoom
      dimension ix(*),iplt(40),x(ndm,*)
      call pltord(iel, iu,iplt)
      v = 0.0
      if(ndm.eq.3) v = x(3,iplt(1))
      if(zoom(x,ndm,nel)) then
       ia = 1
       if(ipla.eq.2) ia = 3
       if(ipc.eq.2)  ia = 3
       call plotl(x(1,iplt(1)),x(2,iplt(1)),v,ia)
       do 100 i = 2,iu
        j = iplt(i)
        if((j.gt.nel).or.(ix(j).eq.0).or.(j.eq.0)) go to 100
        if(ndm.eq.3) v = x(3,j)
        call plotl(x(1,j),x(2,j),v,2)
100    continue
       if(ipc.eq.2) goto 200
         if(ia.eq.1) call clpan
200    continue
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine plot9s(ix,x,ndm,nel)
c-----------------------------------------------------------------------
c
c.... Purpose: plot forces for 2d/3d beam elements
c
c     Input:
c
c     Output:
c
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata2
      implicit double precision(a-h,o-z)
cww      logical zoom
      dimension ix(*),x(3,*)
      v = 0.0
      if(ndm.eq.3) v = x(3,1)
cww      if(zoom(x,ndm,nel)) then
       call plotl(x(1,1),x(2,1),v,ipgl)
       do 100 i = 1,nel
        j = ix(i+1)
        if(ndm.eq.3) v = x(3,j)
        call plotl(x(1,j),x(2,j),v,2)
100    continue
       if(ipgl.eq.1) call clpan
cww      end if
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine plot2dh(ie,ix,x,xl,nie,ndm,nen,nen1,numel,k2,k3,k4,
     1                   idis)
c----------------------------------------------------------------------
c
c.... Purpose: plot mesh/deformed mesh in hidden line technique
c
c     Input:
c
c     Output:
c
c     Comment:
c     plot elements from back to front
c     plot interior part of element in background color
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------
      USE eldata
      implicit double precision (a-h,o-z)
      logical zoom
      dimension ie(nie,*),ix(nen1,*),x(ndm,*),xl(ndm,*),idis(*)
cww   if(ndm.lt.3) return
c
c.... initialize plot, draw border and add header
c.... choose direction
         i1 =  1
         i2 =  numel
         inc = 1
      if(k4.eq.2) then
         i1 = numel
         i2 = 1
         inc = -1
      end if
c...  loop over elements to draw mesh
        do 105 n = i1,i2,inc
          new=idis(n)
          ma = ix(nen1,new)
          if(iplma(ma).eq.0) goto 105
          do 104 i = 1,nen
            ii = abs(ix(i,new))
            if(ii) 100,100,102
100         do 101 j = 1,ndm
             xl(j,i) = 0.0
101         continue
          go to 104
102       nel = i
            do 103 j = 1,ndm
              xl(j,i) = x(j,ii)
103         continue
104       continue
c.....  plot element
        if(zoom(xl,ndm,nen))
     +    call plot9h(ie(nie-1,ma),ix(1,new),xl,ndm,nel,k2,k3)
105     continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine plot9h(iel,ix,x,ndm,nel,k2,k3)
c-----------------------------------------------------------------------
c
c.... Purpose: plot elements, line in colo k2, interior part in colo k3
c
c     Input:
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata2
      implicit double precision(a-h,o-z)
      dimension ix(*),iplt(40),x(ndm,*)
      call pltord(iel, iu,iplt)
      v = 0.0
      if(ndm.eq.3) v = x(3,iplt(1))
c.... plot interior part
      icol = k3
      call pppcol(icol)
      call plotl(x(1,iplt(1)),x(2,iplt(1)),v,ipgl)
      do 100 i = 2,iu
        j = iplt(i)
      if((j.gt.nel).or.(j.eq.0)) go to 100
      if(ix(j).eq.0)             go to 100
        if(ndm.eq.3) v = x(3,j)
        call plotl(x(1,j),x(2,j),v,2)
100   continue
      if(ipgl.eq.1) call clpan
c.... plot element border
      call pppcol(k2)
      call plotl(x(1,iplt(1)),x(2,iplt(1)),v,3)
      do 200 i = 2,iu
        j = iplt(i)
      if((j.gt.nel).or.(j.eq.0)) go to 200
      if(ix(j).eq.0)             go to 200
        if(ndm.eq.3) v = x(3,j)
        call plotl(x(1,j),x(2,j),v,2)
200   continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine plotl(xx1,xx2,xx3,ipen)
c-----------------------------------------------------------------------
c
c.... Purpose: plot line
c
c     Input:
c       xx1  - coordinate x1
c       xx2  - coordinate x2
c       xx3  - coordinate x3
c       ipen - action to do, see dplot
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata1
      USE ppers
      implicit double precision (a-h,o-z)
      logical errv
      dimension xg(3)
c.... recover plot coordinates
      x1 = xx1
      x2 = xx2
      x3 = xx3
c.... perform perspective tranformation
      if(kpers.eq.1) then
        xg(1) = x1
        xg(2) = x2
        xg(3) = x3
        call perspj(xg,xg,3,3,1,errv)
        x1 = xg(1)
        x2 = xg(2)
        x3 = xg(3)
      end if
c.... compute the normal coordinates
      s1 = scale*(x1 + x1 - sx(1)) + s0(1)
      s2 = scale*(x2 + x2 - sx(2)) + s0(2)
c.... if isometric recompute values
      if(iso) then
         xmul = 2.*max(dx(1),dx(2))/(dx(1)+dx(2))
         te = (s0(1) + 0.5*(s1 - s2))*xmul
         s2 = (0.2885*(s1 + s2))*xmul + scale*(x3+x3)
         s1 = te - 0.1
      end if
      s1 = max(0.0d0,min(1.45d0,s1))
      s2 = max(0.0d0,min(1.00d0,s2))
      call dplot(s1,s2,ipen,1)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine plshowp(idev,ior,fplt,ii,iln,nsizt,cs,fact,xfac,rotang,
     1     deltax,deltay,kpers,ipola,msym,ipma,nummat,iso,defo,iopl)
c-----------------------------------------------------------------------
c
c.... Purpose: show plot informations
c
c     Input:
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE fornam
      USE strnam
      implicit double precision (a-h,o-z)
      logical iso,defo
      character*229 fplt
      dimension xfac(3),rotang(3),ipma(*)
      nc   = ipos(fplt,229)
      if(ior.lt.0) then
cww     if(idev.eq.4) call clwopen('Show Values PLOT',140,60,450,420,1,2)
        if(idev.eq.4) call clwopen('Show Values PLOT',1,2)
        if(idev.lt.4) write(*,2000)
        write(*,2001) fplt(1:nc),ii,iln,nsizt,cs,fact,xfac,rotang,
     +              deltax,deltay,ipola,msym
        ia=1
        ie=min(10,nummat)
        do i=1,nummat,10
          write(*,2002)(ipma(k),k=ia,ie)
          ia=ia+10
          ie=ie+10
          if(ie.gt.nummat) ie=nummat
        end do
        if(kpers.eq.0) then
          if(iso)      write(*,2005)
          if(.not.iso) write(*,2006)
        else
          write(*,2007)
        end if
        if(defo)      write(*,2016)
        if(.not.defo) write(*,2017)
c...    force names
        write(*,2018)
        do i=1,11
          if(forsus(i).ne.'               ') then
            write(*,2020) i, forsus(i)
          end if
        end do

c...    stress names
        write(*,2019)
        do i=1,25
          if(strsus(i).ne.'               ') then
            write(*,2020) i, strsus(i)
          end if
        end do

        if(idev.eq.4) then
          call clwclose(1,2,1)
        end if
      end if
      return
2000  format(3x,'P L O T    P A R A M E T E R S'//)
2001  format(  5x,'Plot file   = ',a/
     2       5x,'Frame No.   =',i3/5x,'Line Type   =',i3/
     2       5x,'Size of text=',i3/
     3       5x,'Disp. Scale =',e12.5/5x,'Fact  Scale =',e12.5/
     4       5x,'Scale x1    =',e12.5/5x,'Scale x2    =',e12.5/
     5       5x,'Scale x3    =',e12.5/5x,'Angl rot1   =',f12.3/
     6       5x,'Angl rot2   =',f12.3/5x,'Angl rot3   =',f12.3/
     7       5x,'Move 1-dir  =',f12.3/5x,'Move 2-dir  =',f12.3/
     8       5x,'Polar Coord.=',i4/   5x,'Symm.Cond.  =',i4)
2002  format(5x,'Mat.to plot =',10i5)
2005  format(5x,'Isometric view')
2006  format(5x,'Cartesian view')
2007  format(5x,'Perspective view')
2016  format(5x,'Plot on   deformed mesh')
2017  format(5x,'Plot on undeformed mesh')
2018  format(5x,'Names of Forces to plot')
2019  format(5x,'Names of Stresses to plot')
2020  format(5x,i3,a15)
      end
c
      subroutine pltangl(x,angl,ndm,numnp,n1)
c-----------------------------------------------------------------------
c
c.... Purpose: plot position of all nodes with angles
c
c     Input:
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata1
      USE pdata2
      implicit double precision (a-h,o-z)
      logical zoom
      dimension x(ndm,*),angl(numnp)
      dx1 = .002/scale
      x3 = 0.0
      do 100 n = 1,numnp
        angn = angl(n)
        if(angn.eq.0.d0) goto 100
        if(zoom(x(1,n),ndm,1)) then
          x1 = x(1,n)
          x2 = x(2,n)
          if(ndm.ge.3) x3 = x(3,n)
          call plotl(x1-dx1 , x2+dx1 , x3, 3)
          call plotl(x1-dx1 , x2+dx1 , x3, ipgl)
          call plotl(x1-dx1 , x2-dx1 , x3, 2)
          call plotl(x1+dx1 , x2-dx1 , x3, 2)
          call plotl(x1+dx1 , x2+dx1 , x3, 2)
          call plotl(x1-dx1 , x2+dx1 , x3, 2)
          if(ipgl.eq.1) call clpan
          call plotl(x1-12*dx1 , x2+2*dx1 , x3, 3)
          if(n1.ne.0) call plabl(n)
        end if
100   continue
      return
      end
c
c----------------------------------------------------------------------+
c
      subroutine pltangl1(x,angl,tra,vr,ndm,numnp,ipgl,xm)
c----------------------------------------------------------------------+
c
c.... Purpose: plot position of all nodes with angles and basis
c
c     Input:
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------+
      USE mdat2
      USE pdata1
      implicit double precision (a-h,o-z)
      logical zoom
      dimension x(ndm,*),angl(*),x0(3),x0z(3),tra(3,3),vr(3),xx(3,5),
     +          diri(3,3)
      pi=datan(1.d0)*4.d0
      xm = xm/scale/40.
c      if(ndm.lt.3) return
c...  dofs for angl
      ij1 = ia(1)
      ij2 = ia(2)
c.... basis at nodes with angl .ne.0
      do 100 ii = 1,numnp
        x0(1)  = x(1,ii)
        x0(2)  = x(2,ii)
        x0(3)  = 0.d0
        if(ndm.eq.3) x0(3)  = x(3,ii)
        x0z(1) = x0(1)
        x0z(2) = x0(2)
        x0z(3) = 0.d0
        if(ndm.eq.3) x0z(3) = x0(3)
        if(angl(ii).ne.0.0d0) then
c...      for nodes with angl
          call plxtrn(x0z,tra,vr,3,1)
          if(zoom(x0z,3,1)) then
c.....      get the basis vectors
            call pzero(diri,9)
            ang = angl(ii)*pi/180.d0
            cn  = cos(ang)
            sn  = sin(ang)
            if(ij1.eq.1.and.ij2.eq.2) then
              diri(1,1) =  cn
              diri(1,2) = -sn
              diri(2,1) =  sn
              diri(2,2) =  cn
              diri(3,3) =  1.d0
            else if(ij1.eq.1.and.ij2.eq.3) then
              diri(1,1) =  cn
              diri(1,3) = -sn
              diri(3,1) =  sn
              diri(3,3) =  cn
              diri(2,2) =  1.d0
            else if(ij1.eq.2.and.ij2.eq.3) then
              diri(2,2) =  cn
              diri(2,3) = -sn
              diri(3,2) =  sn
              diri(3,3) =  cn
              diri(1,1) =  1.d0
            end if
c.....      perspective projecton of axes
            do m = 1,ndm
              call pppcol(m)
              call pzero(xx,15)
              do n = 1,ndm
                fac1    = diri(n,m)*xm
                xx(n,1) = x0(n)
                xx(n,2) = xx(n,1) + fac1
                xx(n,5) = xx(n,2)
              enddo
              fac1 = diri(1,m)*xm
              fac2 = diri(2,m)*xm
              fac3 = diri(3,m)*xm
              xx(1,3) = xx(1,2) -.3*fac1 - .1*(fac2+fac3)
              xx(2,3) = xx(2,2) -.3*fac2 + .1*(fac1+fac3)
              xx(3,3) = xx(3,2) -.3*fac3 + .1*(fac1+fac2)
              xx(1,4) = xx(1,2) -.3*fac1 + .1*(fac2+fac3)
              xx(2,4) = xx(2,2) -.3*fac2 - .1*(fac1+fac3)
              xx(3,4) = xx(3,2) -.3*fac3 - .1*(fac1+fac2)
c.....        transform vector if necessary
              call plxtrn(xx,tra,vr,ndm,5)
c.....        plot the vector
              call plotl(xx(1,1),xx(2,1),xx(3,1),3)
              call plotl(xx(1,2),xx(2,2),xx(3,2),2)
              call plotl(xx(1,2),xx(2,2),xx(3,2),ipgl)
              call plotl(xx(1,3),xx(2,3),xx(3,3),2)
              call plotl(xx(1,4),xx(2,4),xx(3,4),2)
              call plotl(xx(1,2),xx(2,2),xx(3,2),2)
              if(ipgl.eq.1) call clpan
              call plotl(xx(1,2),xx(2,2),xx(3,2),3)
              call plabl(m)
            enddo
          end if
        end if
100   continue
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine pltaxs(tr,vr,xi,ndm,ct)
c----------------------------------------------------------------------
c
c.... Purpose: plot vectors for axes, including rot-macro
c
c     Input:
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------
      USE pdata2
      implicit double precision (a-h,o-z)
      dimension tr(3,3),xx(3,5),xi(3),vr(3),r(3)
cwd   limit ndm to 3 (as in IGA 4D projective coordinates are used),
cwd   ndm in original is replaced by ndm_local in whole subroutine
      ndm_local = min(ndm,3)
c.... compute plot location for axes and transform xi
      do 50 i = 1,ndm_local
          r(i) = tr(1,i)*xi(1)+tr(2,i)*xi(2)+tr(3,i)*xi(3)+vr(i)
50      continue
        xi(1) = r(1)
        xi(2) = r(2)
        xi(3) = 0.d0
        if (ndm_local.eq.3) xi(3) = r(3)
c.... perspective projecton of axes
      do 120 m = 1,ndm_local
      call pzero(xx,15)
      do 100 n = 1,ndm_local
      fac1 = tr(m,n)*ct
      xx(n,1) = xi(n)
      xx(n,2) = xx(n,1) + fac1
100   xx(n,5) = xx(n,2)
      fac1 = tr(m,1)*ct
      fac2 = tr(m,2)*ct
      fac3 = tr(m,3)*ct
      xx(1,3) = xx(1,2) -.3*fac1 - .1*(fac2+fac3)
      xx(2,3) = xx(2,2) -.3*fac2 + .1*(fac1+fac3)
      xx(3,3) = xx(3,2) -.3*fac3 + .1*(fac1+fac2)
      xx(1,4) = xx(1,2) -.3*fac1 + .1*(fac2+fac3)
      xx(2,4) = xx(2,2) -.3*fac2 - .1*(fac1+fac3)
      xx(3,4) = xx(3,2) -.3*fac3 - .1*(fac1+fac2)
c.... plot the vector
      call plotl(xx(1,1),xx(2,1),xx(3,1),3)
      call plotl(xx(1,2),xx(2,2),xx(3,2),2)
      call plotl(xx(1,2),xx(2,2),xx(3,2),ipgl)
      call plotl(xx(1,3),xx(2,3),xx(3,3),2)
      call plotl(xx(1,4),xx(2,4),xx(3,4),2)
      call plotl(xx(1,2),xx(2,2),xx(3,2),2)
      if(ipgl.eq.1) call clpan
c      call plotl(xx(1,1),xx(2,1),xx(3,1),3)
c      do 110 n = 2,5
c      call plotl(xx(1,n),xx(2,n),xx(3,n),2)
c110   continue
      call plotl(xx(1,2),xx(2,2),xx(3,2),3)
      call plabl(m)
120   continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltbou(id,x,angl,ix,nen1,ndm,ndf,numnp,k2,rk4,tra,vr)
c-----------------------------------------------------------------------
c
c.... Purpose: Plot 1 - 3-D boundary constraints
c
c     Input:
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE boun
      USE mdat2
      USE pdata1
      USE pdata7
      implicit double precision (a-h,o-z)
      logical zoom
      dimension id(ndf,*),x(ndm,*),angl(*),tra(3,3),vr(3),xz(3)
      integer ix(nen1,*)
      pi=datan(1.d0)*4.d0
      if(scale.eq.0.d0) then
        call drawmess(' Scale = 0 in PLTBOU',1,1)
        return
      end if
c...  test if b.c.occur
      do n = 1,numnp
        do idf = 1,ndf
            if (id(idf,n) .lt. 0) goto 99
        enddo
      enddo
      call drawmess(' No Boundary conditions available',1,1)
      return
c.... plot boundary restraints (lines = fixed)
99    dx1 = 0.015d0/scale*rk4
      do 100 n = 1,numnp
        if(iplmano(ix,n,nen1).eq.0)  goto 100  ! only if matn
c....   local boundary cond., rotate in plane ij1,ij2
        if(angl(n).ne.0.0d0) then
          ang = angl(n)*pi/180.d0
          cs  = cos(ang)
          sn  = sin(ang)
          ij1 = ia(1)
          ij2 = ia(2)
          if(ij1.eq.1 .and. ij2.eq.2) ij3 = 3
          if(ij1.eq.1 .and. ij2.eq.3) ij3 = 2
          if(ij1.eq.2 .and. ij2.eq.3) ij3 = 1
        else
          cs  = 1.d0
          sn  = 0.d0
          ij1 = 1
          ij2 = 2
          ij3 = 3
        end if
c....   rotation matrix
        call pzero(trb,9)
        trb(ij1,ij1) =  cs
        trb(ij2,ij2) =  cs
        trb(ij1,ij2) =  sn
        trb(ij2,ij1) = -sn
        trb(ij3,ij3) =  1.d0
        call pzero(vrb,3)
c
        xz(1) = x(1,n)
        xz(2) = x(2,n)
        xz(3) = 0.d0
        if(ndm.ge.3) xz(3) = x(3,n)
        call plxtrn(xz,tra,vr,ndm,1)
        if(zoom(xz,ndm,1)) then
          x1 = x(1,n)
          x2 = x(2,n)
          x3 = 0.d0
          if(ndm.ge.3) x3 = x(3,n)
c.....    Plate(u_3,phi_1,phi_2)
          if(ipla.eq.1) then
            if (id(1,n) .le. 0 .and. (k2.eq.1.or.k2.eq.0))
     +        call pltbou1(3,tra,vr,ndf,1)
            if (id(2,n) .le. 0 .and. (k2.eq.2.or.k2.eq.0))
     +        call pltbou1(4,tra,vr,ndf,1)
            if (id(3,n) .le. 0 .and. (k2.eq.3.or.k2.eq.0))
     +        call pltbou1(5,tra,vr,ndf,1)
c.....    2-D Beam, Axishell(u_1,u_2,phi_3)
          else if(ipla.eq.2) then
            if (id(1,n) .le. 0.and. (k2.eq.1.or.k2.eq.0))
     +        call pltbou1(1,tra,vr,ndf,1)
            if (id(2,n) .le. 0.and. (k2.eq.2.or.k2.eq.0))
     +        call pltbou1(2,tra,vr,ndf,1)
            if (id(3,n) .le. 0.and. (k2.eq.3.or.k2.eq.0))
     +        call pltbou1(7,tra,vr,ndf,1)
          else
c.....      other problems
c.....      1-D problem (u_1) + higher order elements
            if (id(1,n) .le. 0.and. (k2.eq.1.or.k2.eq.0))
     +          call pltbou1(1,tra,vr,ndf,1)
            if (ndm.ge.2. and. ndf.ge.2 ) then
c.....      2-D truss, plane stress (u_1,u_2) + higher order elements
              if (id(2,n) .le. 0.and. (k2.eq.2.or.k2.eq.0))
     +            call pltbou1(2,tra,vr,ndf,1)
            end if
c.....      3-D truss,3-D brick (u_1,u_2,u_3) + higher order elements
            if (ndm.ge.3. and. ndf.ge.3 ) then
              if (id(3,n) .le. 0.and. (k2.eq.3.or.k2.eq.0))
     +            call pltbou1(3,tra,vr,ndf,1)
            end if
c.....      3-D beam, shell(5),shell(6) (u_1,u_2,u_3,phi_1,phi_2,phi_3)
            if (ndm.ge.3. and. ndf.ge.4 ) then
              if (id(4,n) .le. 0.and. (k2.eq.4.or.k2.eq.0))
     +            call pltbou1(4,tra,vr,ndf,1)
            end if
c.....      3-D beam, shell(5), shell(6) (u_1,u_2,u_3,phi_1,phi_2,phi_3)
            if (ndm.ge.3. and. ndf.ge.5 ) then
              if (id(5,n) .le. 0.and.(k2.eq.5.or.k2.eq.0))
     +            call pltbou1(5,tra,vr,ndf,1)
            end if
c.....      3-D beam, shell(6) (u_1,u_2,u_3,phi_1,phi_2,phi_3)
            if (ndm.ge.3. and. ndf.ge.6 ) then
              if (id(6,n) .le. 0.and. (k2.eq.6.or.k2.eq.0))
     +            call pltbou1(6,tra,vr,ndf,1)
            end if
c.....      elements with more then 6 dofs
            if(ndf.gt.6) then
              iaddf = 0
              do 200 ii=7,ndf
                if(id(ii,n).le.0.and. (k2.eq.ii.or.k2.eq.0)) iaddf = 1
200           continue
              if(iaddf.eq.1) call pltbou1(8,tra,vr,ndf,1)
            end if
          end if
        end if
100   continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltbou1(isw,tra,vr,ndf,icol)
c-----------------------------------------------------------------------
c
c.... Purpose: Plot boundary vector
c
c     Input:
c       isw = 1-3 disp
c       isw = 4-6 rot
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE boun
      USE mdat2
      implicit double precision (a-h,o-z)
      dimension xx(3,6),tra(3,3),vr(3),tri(3,3),vri(3),symb(3,6)
      dx2= dx1*0.5d0
      call pzero(xx,18)
      call pzero(tri,9)
      call pzero(vri,3)
      irot = 0
c.... define symbol in 1-2 plane  1-3 in red, 4-6 in orange
      ddx = 0.0d0
      ic  = 2
      if(isw.ge.4 .and. isw.le.6) then
cww     ddx = dx1  ! plot symbol not directly at node
        ic  = 6    ! color
        if(itrot.eq.1) irot = 1  ! do not rotate for angles
      end if
      call pzero(symb,18)
      symb(2,1) =          -ddx
      symb(1,2) =  - dx2
      symb(2,2) =  - dx1   -ddx
      symb(1,3) =  + dx2
      symb(2,3) =  - dx1   -ddx
      symb(2,4) =  - dx1   -ddx
      symb(3,4) =  - dx2
      symb(2,5) =  - dx1   -ddx
      symb(3,5) =  + dx2
      symb(2,6) =  - dx1   -ddx
c
      goto(1,2,3,1,2,3,7,8), isw
c.... boundary x (rotate symbol -90 deg in 1-2 plane)
1     tri(1,2) = -1.d0
      tri(2,1) =  1.d0
      tri(3,3) =  1.d0
      call plxtrn(symb,tri,vri,3,6)
      if(irot.eq.0) call plxtrn(symb,trb,vrb,3,6)
      goto 10
c.... boundary y (do not rotate symbol )
2     if(irot.eq.0) call plxtrn(symb,trb,vrb,3,6)
      goto 10
c.... boundary z (rotate symbol +90 deg in 2-3 plane)
3     tri(1,1) =  1.d0
      tri(2,3) =  1.d0
      tri(3,2) = -1.d0
      call plxtrn(symb,tri,vri,3,6)
      if(irot.eq.0) call plxtrn(symb,trb,vrb,3,6)
      goto 10
c
c.... angle phi_3 in 2-D problems,no transformation necessary (orange)
7     if(icol.eq.1) call pppcol(8)
      call plotl(x1-dx2,x2-dx2,x3,3)
      call plotl(x1+dx2,x2-dx2,x3,2)
      call plotl(x1+dx2,x2+dx2,x3,2)
      call plotl(x1-dx2,x2+dx2,x3,2)
      call plotl(x1-dx2,x2-dx2,x3,2)
      return
c.... higher order boundarys (yellow)
8      xx(1,1) = x1 - dx1
       xx(2,1) = x2 - 0.2*dx2
       xx(3,1) = x3
       xx(1,2) = x1 + dx1
       xx(2,2) = x2 - 0.2*dx2
       xx(3,2) = x3
       call plxtrn(xx,tra,vr,3,2)
       if(icol.eq.1) call pppcol(7)
       call plotl(xx(1,1),xx(2,1),xx(3,1),3)
       call plotl(xx(1,2),xx(2,2),xx(3,2),2)
      return
c.... add the rotated symbol to node
10    continue
      do k = 1,6
        xx(1,k) = x1 + symb(1,k)
        xx(2,k) = x2 + symb(2,k)
        xx(3,k) = x3 + symb(3,k)
      enddo
c...  transform the boundary cond.
      call plxtrn(xx,tra,vr,3,6)
c.... plot the boundary cond.
      if(icol.eq.1) call pppcol(ic)
      call plotl(xx(1,1),xx(2,1),xx(3,1),3)
      call plotl(xx(1,2),xx(2,2),xx(3,2),2)
      call plotl(xx(1,3),xx(2,3),xx(3,3),2)
      call plotl(xx(1,1),xx(2,1),xx(3,1),2)
      call plotl(xx(1,4),xx(2,4),xx(3,4),2)
      call plotl(xx(1,5),xx(2,5),xx(3,5),2)
      call plotl(xx(1,1),xx(2,1),xx(3,1),2)
      call plotl(xx(1,6),xx(2,6),xx(3,6),2)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltcon(x,id,ie,ix,b,idis,nie,ndm,ndf,nen1,nxn,nne,ic,
     1                  mc,mmc,iev,k5,cinv)
c-----------------------------------------------------------------------
c
c.... Purpose: Plot contour values for displacements, stresses etc
c
c     Input:
c      x(ndm,numnp)   - nodal coordinates
c      id(ndf,numnp)  - equation numbers for each active dof
c      ie(nie,numat)  - assembly information for material set
c      ix(nen1,nnne) - element nodal connections of mesh
c      b(*)           - array to plot
c      idis
c      nie            - ndf+2
c      ndm            - spatial dimension of mesh
c      ndf            - number dof/node
c      nen1           - = nen+4 or nxd(Hide) 
c      nxn            - = nen   or 4(hide) 
c      nne            - number of elements to plot
c      ic             - number of eg. stress to plot, e.g. stress_1
c      mc             > 0 plot lines   (variable contour 1-14)
c                     <=0 filled plots (variable contour 1-13)
c      mmc            - = 1: stre, 2: disp, 3: velo, 4: acce, 5: eigv,
c                         6: resi
c      iev         - number of eigenvector to plot
c      k5            -    k5 .ne.0: plot element border in color k5
c      cinv          - =  1 colors from -=red ->+=blue set by colo,,cinv
c                      = -1 colors from -=blue->+=red
c
c     Output:
c
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE cdata
      USE contval
      USE errchk
      USE iofile
      USE isogeo
      USE pdata1
      USE pdata2
      USE pdata7
      USE pdatas
      USE plotter
      USE rpdata
      implicit double precision (a-h,o-z)
      character*1 y
      logical zoom,tvc(9,9),vflg,cont
cwd   added isogeo.h isNurbs and some arrays    
      logical isNurbs
      dimension xl_IGA(ndm,nen),v_IGA(nen)  
      dimension Xi1fani(9),Xi1fanip(9)
      dimension Xi2fanj(9),Xi2fanjp(9)
      dimension idis(*)
      integer*4 ipal(14),ipali(14)
      real*8  xl(3,40),x(ndm,*),b(*),v(40),vc(14),vt(4),xt(3,4)
      integer id(ndf,*),ie(nie,*),ix(nen1,*),iplt(40),ilc(4)
      save cont,vl,vu,nc,nnc,vc,xmx,ymx,xmn,ymn,vmn,vmx
c----------------------------------------------------------------------
c       
      iadd1 = 16
      vflg  = ipb.eq.0
      call pzerol ( tvc , .true. , 81 )
      dx1   = .024/scale
c.... symmetry: jump in case of repeated call
      if(icleas.eq.1) go to 11111
c
1     continue
15    continue
      if(icv.eq.1) then ! interactive
        if(ior.lt.0) write(*,2003)
        call dinput(contv,3)
      end if
      nci   = contv(1)
      nc    = min(14,nci)
      if(nc.le.0) nc = 14
      nnc = nc
      if(mc.le.0) nc = nc-1 ! filled plots
      vc(1) = contv(2)
      vc(2) = contv(3)
      if(errck) go to 15
      drm =  dabs(rmx-rmn)
      if( (drm.lt.1.e-5*dabs(rmx).or.drm.lt.1.e-10)
     + .and..not.(abs(mmc).eq.6) ) then
c....   no lines for nearly same values or zero values but not for residuum
        nnc = 1
        nc  = 0
        vc(1) = rmn
        vc(2) = rmx
        vmin  = rmn
        vmax  = rmx

      else if(vc(1).eq.vc(2)) then
c....   no input(default)
        vc(1) = rmn
        vc(2) = rmx
      else
c....   input (value from screen)
        if(nc.gt.1) vc(nc) = vc(2)
      end if
c
c.... generate intermediate values
      if(nc.ge.1) then !
c....   limit values
        vmin = rmn
        vmax = rmx
        vc0=vc(1)
        vcn=vc(2)
        if(vc0.lt.vmin) vmin = vc0
        if(vcn.gt.vmax) vmax = vcn
c....   intermediate values: nc colors -> nc-1 points, nc lines -> nc points
        do i = 1,nc
          vc(i) = vc0 + (vcn-vc0)*i/(nc+1)
        end do
      end if
      if(mc.gt.0) cont = .true.
      if(mc.le.0) cont = .false.
      if (icv.eq.1) then ! interactive
        if(ior.ge.0) write(iow,2000) (vc(i),i=1,nc)
        if(ior.lt.0) write(*  ,2000) (vc(i),i=1,nc)
c....   if interactive offer chance to change inputs
        if(ior.lt.0) then
          write(*,2002)
900       read (*,1000,err=901,end=902) y
          goto  903
901       call  errclr ('PLTCON')
          goto  900
902       call  endclr ('PLTCON',y)
903       if(y .eq. 'c' .or. y .eq. 'C') then
            if(idev.eq.4) call clwclose(1,2,0)
            return
          else if(y .ne. 'y' .and. y .ne. 'Y') then
            go to 1
          end if
          if(idev.eq.4) call clwclose(1,2,0)
        end if
      end if
c.... find max/min of plot variable
      j   = ic
      xmx = x(1,1)
      ymx = x(2,1)
      xmn = x(1,1)
      ymn = x(2,1)
      vmn = b(j)
      vmx = b(j)
cwd   changes for Nurbs
      if (isNurbs()) then
      do i = 1,numnp
        xmx = max(x(1,i),xmx)
        ymx = max(x(2,i),ymx)
        xmn = min(x(1,i),xmn)
        ymn = min(x(2,i),ymn)
      end do
      do i = 1,nne*9
        vmn = min(vmn,rNURstress(j,1))
        vmx = max(vmx,rNURstress(j,1))
        j   = j + ndf
      end do
      else
cwd   standard version
      do 17 i = 1,numnp
        xmx = max(x(1,i),xmx)
        ymx = max(x(2,i),ymx)
        xmn = min(x(1,i),xmn)
        ymn = min(x(2,i),ymn)
        vmn = min(vmn,b(j))
        vmx = max(vmx,b(j))
        j   = j + ndf
17    continue
      end if
cwd   end of changes
      delx  = (xmx-xmn)
      dely  = (ymx-ymn)
      deleps = 1.e-7    !cww  for values nearly zero?
      delx = max(delx,deleps)
      dely = max(dely,deleps)
      xmx = 8.2/delx
      ymx = 8.2/dely
      if(vmn.eq.vmx.and.(drm.ne.0.d0)) then
cww     if(ior.ge.0) write(iow,2001)
        if(ior.lt.0)
     +  call drawmess(' No plot - zero difference in values ',1,0)
        return
      end if
c.... open plot and loop through elements
11111 continue
c.... set colortable for nnc.le.14
      call colpal(ipal,nnc)
c.... set direction for colors
      ncol = nc
      if(mc.gt.0) ncol = nc-1
      call pltcol(ipal,ipali,ncol,mmc,cinv)
c
      if(icleas.eq.0) iclear = 0
      call plopen
      call pzero(xl,120)
c.... set colors for contour plots
      call setcol(2,0)
      if(icleas.eq.0) then
        if(cont) then
          call pltctx(vc,ic,nc,mmc,iev,ipali)
        else
          call pltftx(vc,-mc,nc,mmc,iev,ipali,vmin,vmax)
        end if
      end if
      ic = max(1,min(ic,ndf))
      maold = -1
c.... set for ps(prin) that no border is drawn
      ibps = 1

      do 250 n = 1,nne

c....   pickup new element number
        new = idis(n)
c....   plot only for specified material number
        if(iplma(ix(nen1,new)).eq.0)    goto 250
c....   get plot order for each element
        ma = ie(nie-1,ix(nen1,new))
        if(ma.ne.maold) then  ! then element can be used only with nel=const. WW
          call pltord(ma, iu,iplt)
          ns = iu
          do i = 1,iu-1  ! 102
            ii = iplt(i)
            if(ii.gt.nxn) then
              do j = i,iu-1 ! 101
                iplt(j) = iplt(j+1)
              end do ! 101
              ns = ns - 1
              iplt(iu) = 0
            end if
          end do ! 102
          iu = ns
          maold = ma
        end if
c....   set element values xl, v and set values of vl and vu
        vl = vmx
        vu = vmn
        ns = 0
cim        
        if(nxn.eq.16)then
          iu=16
          iplt(13)=13 !mid-nodes
          iplt(14)=14
          iplt(15)=15
          iplt(16)=16
        end if
cim
        do i = 1,iu ! 110
          ii = ix(iplt(i),new)
          if(ii.gt.0) then
            ns = ns + 1
            xl(1,ns) = x(1,ii)
            xl(2,ns) = x(2,ii)
            if(ndm.ge.3) xl(3,ns) = x(3,ii)
            j = ndf*(ii-1) + ic
            v(ns) = b(j)
            vl = min(vl,v(ns))
            vu = max(vu,v(ns))
          end if
        end do ! 110
        if ( zoom(xl,3,ns-1) ) then
c....     plot triangles using a center node, assign this node to nodes 1,4
c....     3 node: center = node 1
c....     4 node: center = sum 1:4/4 averaging
c....     8 node: center = sum 1:8/8 averaging, is this a good choice?? WW
c....     9 node: nodal values
          if(ns.gt.4) then ! center for 4/8/9/27 node elements
            nsi = 1
            if(nxn.eq.9.or.nxn.eq.27) then ! center for 9/27-node elements             
              ii=(ix(9,new))
              if(ii.gt.0) then
                xt(1,1) = x(1,ii)
                xt(2,1) = x(2,ii)
                if(ndm.ge.3) xt(3,1) = x(3,ii)
                j = ndf*(ii-1) + ic
                vt(1) = b(j)
              else
               go to 111
              end if ! ii
            else  !  center for 4-8 node elements
111           xt(1,1) = 0.d0
              xt(2,1) = 0.d0
              xt(3,1) = 0.d0
              vt(1)   = 0.d0
              xnn     = ns - 1
              do i = 1,xnn ! 120
                xt(1,1) = xt(1,1) + xl(1,i)
                xt(2,1) = xt(2,1) + xl(2,i)
                if(ndm.ge.3) xt(3,1) = xt(3,1) + xl(3,i)
                vt(1) = vt(1) + v(i)
              end do ! 120
              xt(1,1) = xt(1,1)/xnn
              xt(2,1) = xt(2,1)/xnn
              xt(3,1) = xt(3,1)/xnn
              vt(1)   = vt(1)/xnn
            end if ! nxn
          else !  center for 3 node element = node 1
            nsi     = 2
            xt(1,1) = xl(1,1)
            xt(2,1) = xl(2,1)
            xt(3,1) = xl(3,1)
            vt(1)   = v(1)
          end if ! ns
c....     assign center also to node 4
          xt(1,4) = xt(1,1)
          xt(2,4) = xt(2,1)
          xt(3,4) = xt(3,1)
          vt(4)   = vt(1)
cwd       for isogeometric: overload all xl values and load xt(1:3,1)
          if (isNurbs()) then
c         creating local geometry array, needed for deformed mesh
            do 1008 i = 1,nxn
                ii = ix(i,new)
                call evaluateoldix(AInoix,numel,nen,i,new,iiold)
                if(ii.gt.0) then
                 do j = 1,ndm
                     xl_IGA(j,i) = x(j,ii)
                 end do
                end if
1008        continue
c....       get current patch and corresponding orders
            NURp = ngetNURnmpq(n,3,AInipa,AInmpq)
            NURq = ngetNURnmpq(n,4,AInipa,AInmpq)
c....       compute NURBS coordinates of recent element
            ni = ngetNURni(n,AIninc,AInien,AInipa,1)
            nj = ngetNURni(n,AIninc,AInien,AInipa,2)
c....       check if element has zero measure
c            if (NURknv1(ni+1).eq.NURknv1(ni).or.NURknv2(nj+1)
c     +         .eq.NURknv2(nj)) goto 250
            Xi1 = 0.5d0*(rNURknv1(ni)+rNURknv1(ni+1))
            Xi2 = 0.5d0*(rNURknv2(nj)+rNURknv2(nj+1))
c            call physicalCoordinates(Xi1,Xi2,ni,nj,x1,x2,x3,ndm)
            call physCoor(Xi1,Xi2,ni,nj,x1,x2,x3,ndm,xl_IGA,
     +                    NURp,NURq,AInkv1,AInkv2,4,v_IGA,x4)
            xt(1,1)=x1;xt(2,1)=x2;xt(3,1)=x3;vt(1)=x4;
            vt(1) = rNURstress((n-1)*9+5,1)
c            assign center also to node 4
            xt(1,4) = xt(1,1)
            xt(2,4) = xt(2,1)
            xt(3,4) = xt(3,1)
            vt(4)   = vt(1)
c           calculate midside and corner points counter-clockwise
            Xi1fani = (/1.d0,0.5d0,0.d0,0.d0,0.d0,0.5d0,1.d0,1.d0,1.d0/)
            Xi1fanip= (/0.d0,0.5d0,1.d0,1.d0,1.d0,0.5d0,0.d0,0.d0,0.d0/)
            Xi2fanj = (/1.d0,1.d0,1.d0,0.5d0,0.d0,0.d0,0.d0,0.5d0,1.d0/)
            Xi2fanjp= (/0.d0,0.d0,0.d0,0.5d0,1.d0,1.d0,1.d0,0.5d0,0.d0/)
            do i = 1,9
              Xi1 = Xi1fani(i)*rNURknv1(ni)+Xi1fanip(i)*rNURknv1(ni+1)
              Xi2 = Xi2fanj(i)*rNURknv2(nj)+Xi2fanjp(i)*rNURknv2(nj+1)
              call physCoor(Xi1,Xi2,ni,nj,x1,x2,x3,ndm,xl_IGA,
     +                    NURp,NURq,AInkv1,AInkv2,4,v_IGA,x4)
              xl(1,i)=x1;xl(2,i)=x2;xl(3,i)=x3;v(i)=x4  
            end do
cwd         overloading value limits with NURBS values
            vl = vmx
            vu = vmn
            vl = min(vl,vt(1))
            vu = max(vu,vt(1))
            do int = 1,9
                vl = min(vl,v(int))
                vu = max(vu,v(int))
            end do
            ns=9
          end if
cwd       end of modification for isogeometric analysis
cim
c....     loop over subtriangles for 2d 16-node element
          jj=0
          if(nxn.eq.16)then
            jj=1
            jjj=1
            ns=5
            xt(3,1) = 0.d0
            xt(3,2) = 0.d0
            xt(3,3) = 0.d0
          endif  
455       if(nxn.eq.16) call pltcon16(jj,jjj,xl,v,xt,vt)
cim
c....     loop over subtriangles
          do ii = nsi,ns-nsi ! 240
            if(nxn.eq.16) goto 234
c....       set other points on the triangle 1=center, 4=center 
            xt(1,2) = xl(1,ii)
            xt(2,2) = xl(2,ii)
            xt(3,2) = xl(3,ii)
            vt(2)   = v(ii)
            xt(1,3) = xl(1,ii+1)
            xt(2,3) = xl(2,ii+1)
            xt(3,3) = xl(3,ii+1)
            vt(3)   = v(ii+1)
234         if(cont) then
c....         plot all contours which intersect this element
c....         define color table
              do nn = 1,nc ! 230
                vv = vc(nn)
                if(vv.ge.vl.and.vv.le.vu) then
                  call pppcol(iadd1+ipali(nn))
c....             loop over sides of the triangle to find plot points
                  j = 3
                  do i = 1,3 ! 220
                    if(vv.eq.vt(i)) then
                      x1 = xt(1,i)
                      y1 = xt(2,i)
                      z1 = xt(3,i)
                      call plotl(x1,y1,z1,j)
                      j = 2
                    else if((vt(i)-vv)*(vt(i+1)-vv).lt.0.0d0) then
                      s = (vv - vt(i))/(vt(i+1)-vt(i))
                      x1 = xt(1,i) + s*(xt(1,i+1) - xt(1,i))
                      y1 = xt(2,i) + s*(xt(2,i+1) - xt(2,i))
                      z1 = xt(3,i) + s*(xt(3,i+1) - xt(3,i))
                      call plotl(x1,y1,z1,j)
                      j = 2
                    end if ! vv
c....               add lables on contour lines
                    if(vflg.and.j.eq.2) then
                      ivc = (x1-xmn)*xmx + 1
                      ivc = max(1,min(9,ivc))
                      jvc = (y1-ymn)*ymx + 1
                      jvc = max(1,min(9,jvc))
                      if(tvc(ivc,jvc)) then
                        tvc(ivc,jvc) = .false.
                        call plotl(x1-dx1,y1,z1,3)
                        call plabl(nn)
                        call plotl(x1,y1,z1,3)
                      end if ! tvc
                    end if ! vflg
                  end do ! 220
                end if ! vv
              end do ! 230
            else
              call pltcor(3,ilc,vt,vc,nc)
              call pltefl(3,ilc,xt,vt,vc,nc,mmc,ipla,ipgl,idev,ipali)
            end if ! cont
          end do ! 240
cim
            if(nxn.eq.16)then
              if((jj.gt.0).and.(jj.le.9))then !for nxn=16 elements
                goto 455
              else
                goto 444
              end if
            end if
cim
        end if ! zoom
c....   plot element border in color k5
444     if(k5.gt.0) then
          call pppcol(k5)
          if ( zoom(xl,3,iu) ) then
            call plotl(xl(1,1),xl(2,1),xl(3,1),3)
            do ib = 2,iu
              call plotl(xl(1,ib),xl(2,ib),xl(3,ib),2)
            end do
          end if
        end if
250   continue

      ibps = 0
c.... set colors back
      call setcol(1,0)
      return
1000  format(a1)
2000  format('   ------ Contour Values for Plot ------'/
     +      (1x,7e11.4/1x,7e11.4))
cww2001  format(' ERROR: No plot - zero difference in values')
2002  format(' Input values correct? (y or n, c = cancel) > ',$)
2003  format(' Input No. of Lines/Colors, Min. Value, Max. Value for ',
     +       ' Plot'/3x,'>',$)
      end
c
      subroutine pltcon16(jj,jjj,xl,v,xt,vt)
c-----------------------------------------------------------------------
c
c.... Purpose: Plot contour values for 16 node 2D element
c
c     IM/MK BS KIT 11/14
c-----------------------------------------------------------------------
      integer jj,jjj
      real*8  xl(3,40),vt(4),xt(3,4),v(40)
      if(jj.eq.1)then
        xt(1,1) = (xl(1,1)+xl(1,2)+xl(1,13)+xl(1,12))/4
        xt(2,1) = (xl(2,1)+xl(2,2)+xl(2,13)+xl(2,12))/4
        vt(1) = (v(1)+v(2)+v(13)+v(12))/4
        if(jjj.eq.1)then
            xt(1,2) = xl(1,1)
            xt(2,2) = xl(2,1)
            vt(2) = v(1)
            xt(1,3) = xl(1,2)
            xt(2,3) = xl(2,2)
            vt(3) = v(2)
        end if
        if(jjj.eq.2)then
            xt(1,2) = xl(1,2)
            xt(2,2) = xl(2,2)
            vt(2) = v(2)
            xt(1,3) = xl(1,13)
            xt(2,3) = xl(2,13)
            vt(3) = v(13)
        end if
        if(jjj.eq.3)then
            xt(1,2) = xl(1,13)
            xt(2,2) = xl(2,13)
            vt(2) = v(13)
            xt(1,3) = xl(1,12)
            xt(2,3) = xl(2,12)
            vt(3) = v(12)
        end if
        if(jjj.eq.4)then
            xt(1,2) = xl(1,12)
            xt(2,2) = xl(2,12)
            vt(2) = v(12)
            xt(1,3) = xl(1,1)
            xt(2,3) = xl(2,1)
            vt(3) = v(1)
        end if
      else if(jj.eq.2)then
        xt(1,1) = (xl(1,2)+xl(1,3)+xl(1,14)+xl(1,13))/4
        xt(2,1) = (xl(2,2)+xl(2,3)+xl(2,14)+xl(2,13))/4
        vt(1) = (v(2)+v(3)+v(14)+v(13))/4
        if(jjj.eq.1)then
            xt(1,2) = xl(1,2)
            xt(2,2) = xl(2,2)
            vt(2) = v(2)
            xt(1,3) = xl(1,3)
            xt(2,3) = xl(2,3)
            vt(3) = v(3)
        end if
        if(jjj.eq.2)then
            xt(1,2) = xl(1,3)
            xt(2,2) = xl(2,3)
            vt(2) = v(3)
            xt(1,3) = xl(1,14)
            xt(2,3) = xl(2,14)
            vt(3) = v(14)
        end if
        if(jjj.eq.3)then
            xt(1,2) = xl(1,14)
            xt(2,2) = xl(2,14)
            vt(2) = v(14)
            xt(1,3) = xl(1,13)
            xt(2,3) = xl(2,13)
            vt(3) = v(13)
        end if
        if(jjj.eq.4)then
            xt(1,2) = xl(1,13)
            xt(2,2) = xl(2,13)
            vt(2) = v(13)
            xt(1,3) = xl(1,2)
            xt(2,3) = xl(2,2)
            vt(3) = v(2)
        end if
      else if(jj.eq.3)then
        xt(1,1) = (xl(1,3)+xl(1,4)+xl(1,5)+xl(1,14))/4
        xt(2,1) = (xl(2,3)+xl(2,4)+xl(2,5)+xl(2,14))/4
        vt(1) = (v(3)+v(4)+v(5)+v(14))/4
        if(jjj.eq.1)then
            xt(1,2) = xl(1,3)
            xt(2,2) = xl(2,3)
            vt(2) = v(3)
            xt(1,3) = xl(1,4)
            xt(2,3) = xl(2,4)
            vt(3) = v(4)
        end if
        if(jjj.eq.2)then
            xt(1,2) = xl(1,4)
            xt(2,2) = xl(2,4)
            vt(2) = v(4)
            xt(1,3) = xl(1,5)
            xt(2,3) = xl(2,5)
            vt(3) = v(5)
        end if
        if(jjj.eq.3)then
            xt(1,2) = xl(1,5)
            xt(2,2) = xl(2,5)
            vt(2) = v(5)
            xt(1,3) = xl(1,14)
            xt(2,3) = xl(2,14)
            vt(3) = v(14)
        end if
        if(jjj.eq.4)then
            xt(1,2) = xl(1,14)
            xt(2,2) = xl(2,14)
            vt(2) = v(14)
            xt(1,3) = xl(1,3)
            xt(2,3) = xl(2,3)
            vt(3) = v(3)
        end if
      else if(jj.eq.4)then
        xt(1,1) = (xl(1,12)+xl(1,13)+xl(1,16)+xl(1,11))/4
        xt(2,1) = (xl(2,12)+xl(2,13)+xl(2,16)+xl(2,11))/4
        vt(1) = (v(12)+v(13)+v(16)+v(11))/4
        if(jjj.eq.1)then
            xt(1,2) = xl(1,12)
            xt(2,2) = xl(2,12)
            vt(2) = v(12)
            xt(1,3) = xl(1,13)
            xt(2,3) = xl(2,13)
            vt(3) = v(13)
        end if
        if(jjj.eq.2)then
            xt(1,2) = xl(1,13)
            xt(2,2) = xl(2,13)
            vt(2) = v(13)
            xt(1,3) = xl(1,16)
            xt(2,3) = xl(2,16)
            vt(3) = v(16)
        end if
        if(jjj.eq.3)then
            xt(1,2) = xl(1,16)
            xt(2,2) = xl(2,16)
            vt(2) = v(16)
            xt(1,3) = xl(1,11)
            xt(2,3) = xl(2,11)
            vt(3) = v(11)
        end if
        if(jjj.eq.4)then
            xt(1,2) = xl(1,11)
            xt(2,2) = xl(2,11)
            vt(2) = v(11)
            xt(1,3) = xl(1,12)
            xt(2,3) = xl(2,12)
            vt(3) = v(12)
        end if
      else if(jj.eq.5)then
        xt(1,1) = (xl(1,13)+xl(1,14)+xl(1,15)+xl(1,16))/4
        xt(2,1) = (xl(2,13)+xl(2,14)+xl(2,15)+xl(2,16))/4
        vt(1) = (v(13)+v(14)+v(15)+v(16))/4
        if(jjj.eq.1)then
            xt(1,2) = xl(1,13)
            xt(2,2) = xl(2,13)
            vt(2) = v(13)
            xt(1,3) = xl(1,14)
            xt(2,3) = xl(2,14)
            vt(3) = v(14)
        end if
        if(jjj.eq.2)then
            xt(1,2) = xl(1,14)
            xt(2,2) = xl(2,14)
            vt(2) = v(14)
            xt(1,3) = xl(1,15)
            xt(2,3) = xl(2,15)
            vt(3) = v(15)
        end if
        if(jjj.eq.3)then
            xt(1,2) = xl(1,15)
            xt(2,2) = xl(2,15)
            vt(2) = v(15)
            xt(1,3) = xl(1,16)
            xt(2,3) = xl(2,16)
            vt(3) = v(16)
        end if
        if(jjj.eq.4)then
            xt(1,2) = xl(1,16)
            xt(2,2) = xl(2,16)
            vt(2) = v(16)
            xt(1,3) = xl(1,13)
            xt(2,3) = xl(2,13)
            vt(3) = v(13)
        end if
      else if(jj.eq.6)then
        xt(1,1) = (xl(1,14)+xl(1,5)+xl(1,6)+xl(1,15))/4
        xt(2,1) = (xl(2,14)+xl(2,5)+xl(2,6)+xl(2,15))/4
        vt(1) = (v(14)+v(5)+v(6)+v(15))/4
        if(jjj.eq.1)then
            xt(1,2) = xl(1,14)
            xt(2,2) = xl(2,14)
            vt(2) = v(14)
            xt(1,3) = xl(1,5)
            xt(2,3) = xl(2,5)
            vt(3) = v(5)
        end if
        if(jjj.eq.2)then
            xt(1,2) = xl(1,5)
            xt(2,2) = xl(2,5)
            vt(2) = v(5)
            xt(1,3) = xl(1,6)
            xt(2,3) = xl(2,6)
            vt(3) = v(6)
        end if
        if(jjj.eq.3)then
            xt(1,2) = xl(1,6)
            xt(2,2) = xl(2,6)
            vt(2) = v(6)
            xt(1,3) = xl(1,15)
            xt(2,3) = xl(2,15)
            vt(3) = v(15)
        end if
        if(jjj.eq.4)then
            xt(1,2) = xl(1,15)
            xt(2,2) = xl(2,15)
            vt(2) = v(15)
            xt(1,3) = xl(1,14)
            xt(2,3) = xl(2,14)
            vt(3) = v(14)
        end if
      else if(jj.eq.7)then
        xt(1,1) = (xl(1,11)+xl(1,16)+xl(1,9)+xl(1,10))/4
        xt(2,1) = (xl(2,11)+xl(2,16)+xl(2,9)+xl(2,10))/4
        vt(1) = (v(11)+v(16)+v(9)+v(10))/4
        if(jjj.eq.1)then
            xt(1,2) = xl(1,11)
            xt(2,2) = xl(2,11)
            vt(2) = v(11)
            xt(1,3) = xl(1,16)
            xt(2,3) = xl(2,16)
            vt(3) = v(16)
        end if
        if(jjj.eq.2)then
            xt(1,2) = xl(1,16)
            xt(2,2) = xl(2,16)
            vt(2) = v(16)
            xt(1,3) = xl(1,9)
            xt(2,3) = xl(2,9)
            vt(3) = v(9)
        end if
        if(jjj.eq.3)then
            xt(1,2) = xl(1,9)
            xt(2,2) = xl(2,9)
            vt(2) = v(9)
            xt(1,3) = xl(1,10)
            xt(2,3) = xl(2,10)
            vt(3) = v(10)
        end if
        if(jjj.eq.4)then
            xt(1,2) = xl(1,10)
            xt(2,2) = xl(2,10)
            vt(2) = v(10)
            xt(1,3) = xl(1,11)
            xt(2,3) = xl(2,11)
            vt(3) = v(11)
        end if
      else if(jj.eq.8)then
        xt(1,1) = (xl(1,16)+xl(1,15)+xl(1,8)+xl(1,9))/4
        xt(2,1) = (xl(2,16)+xl(2,15)+xl(2,8)+xl(2,9))/4
        vt(1) = (v(16)+v(15)+v(8)+v(9))/4
        if(jjj.eq.1)then
            xt(1,2) = xl(1,16)
            xt(2,2) = xl(2,16)
            vt(2) = v(16)
            xt(1,3) = xl(1,15)
            xt(2,3) = xl(2,15)
            vt(3) = v(15)
        end if
        if(jjj.eq.2)then
            xt(1,2) = xl(1,15)
            xt(2,2) = xl(2,15)
            vt(2) = v(15)
            xt(1,3) = xl(1,8)
            xt(2,3) = xl(2,8)
            vt(3) = v(8)
        end if
        if(jjj.eq.3)then
            xt(1,2) = xl(1,8)
            xt(2,2) = xl(2,8)
            vt(2) = v(8)
            xt(1,3) = xl(1,9)
            xt(2,3) = xl(2,9)
            vt(3) = v(9)
        end if
        if(jjj.eq.4)then
            xt(1,2) = xl(1,9)
            xt(2,2) = xl(2,9)
            vt(2) = v(9)
            xt(1,3) = xl(1,16)
            xt(2,3) = xl(2,16)
            vt(3) = v(16)
        end if
      else if(jj.eq.9)then
        xt(1,1) = (xl(1,15)+xl(1,6)+xl(1,7)+xl(1,8))/4
        xt(2,1) = (xl(2,15)+xl(2,6)+xl(2,7)+xl(2,8))/4
        vt(1) = (v(15)+v(6)+v(7)+v(8))/4
        if(jjj.eq.1)then
            xt(1,2) = xl(1,15)
            xt(2,2) = xl(2,15)
            vt(2) = v(15)
            xt(1,3) = xl(1,6)
            xt(2,3) = xl(2,6)
            vt(3) = v(6)
        end if
        if(jjj.eq.2)then
            xt(1,2) = xl(1,6)
            xt(2,2) = xl(2,6)
            vt(2) = v(6)
            xt(1,3) = xl(1,7)
            xt(2,3) = xl(2,7)
            vt(3) = v(7)
        end if
        if(jjj.eq.3)then
            xt(1,2) = xl(1,7)
            xt(2,2) = xl(2,7)
            vt(2) = v(7)
            xt(1,3) = xl(1,8)
            xt(2,3) = xl(2,8)
            vt(3) = v(8)
        end if
        if(jjj.eq.4)then
            xt(1,2) = xl(1,8)
            xt(2,2) = xl(2,8)
            vt(2) = v(8)
            xt(1,3) = xl(1,15)
            xt(2,3) = xl(2,15)
            vt(3) = v(15)
        end if
      end if
      jjj=jjj+1
      if(jjj.eq.5)then
        jj=jj+1
        jjj=1
      end if
c.....assign center also to node 4          
      xt(1,4) = xt(1,1)
      xt(2,4) = xt(2,1)
      xt(3,4) = xt(3,1)
      vt(4)   = vt(1)
      return
      end

c
c-----------------------------------------------------------------------
c
      subroutine pltconf(idev,xmaxf,xminf,k1,ityp)
c-----------------------------------------------------------------------
c
c.... Purpose: Plot contour values forces similar to pltcon
c
c     Input:
c       k1   number of force
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE contval
      USE errchk
      USE forccont
      USE iofile
      USE iwinio
      USE rpdata
      implicit double precision (a-h,o-z)
      character*1 y
      integer*4 ipal(14)
c
      rmx = xmaxf
      rmn = xminf
      if(xmaxf.eq.-1.e+38.or.xminf.eq.1.e+38) return
      if(icv.eq.1.and.idev.eq.4.and.ior.lt.0)
cww     + call clwopen('Values of    PROFILE',1,iwys-190,630,230,1,2)
     + call clwopen('Values of    PROFILE',1,2)
1     continue
      if(icv.eq.1) then ! interactive
        if(ior.lt.0) write(*,2003)
        call dinput(contv,3)
      end if
      nci   = contv(1)
      nc    = min(14,nci)
      if(nc.le.0) nc = 14
      nnc = nc
c.....filled plots
      nc = nc-1
c
      vc(1) = contv(2)
      vc(2) = contv(3)
      drm =  dabs(rmx-rmn)
      if( (drm.lt.1.e-5*dabs(rmx).or.drm.lt.1.e-10)) then
c....   no lines for nearly same values or zero values
        nnc = 1
        nc  = 0
        vc(1) = rmn
        vc(2) = rmx
        vmin  = rmn
        vmax  = rmx
      else if(vc(1).eq.vc(2)) then
c....   no input(default)
        vc(1) = rmn
        vc(2) = rmx
      else
c....   input (value from screen)
        if(nc.gt.1) vc(nc) = vc(2)
      end if
c.... generate intermediate values
      if(nc.ge.1) then
c....   limit values
        vmin = rmn
        vmax = rmx
        vc0=vc(1)
        vcn=vc(2)
        if(vc0.lt.vmin) vmin = vc0
        if(vcn.gt.vmax) vmax = vcn
c....   intermediate values: nc colors -> nc-1 points, nc lines -> nc points
        do i = 1,nc
          vc(i) = vc0 + (vcn-vc0)*i/(nc+1)
        end do
      end if
c
      if(nc.eq.1) vc(1)=0.d0  ! plot only in red(-) and blue(+) margin=0!
c
      if(icv.eq.1) then ! interactive
        if(ior.ge.0) write(iow,2000) (vc(i),i=1,nc)
        if(ior.lt.0) write(*  ,2000) (vc(i),i=1,nc)
        if(ior.lt.0) then
          write(*,2002)
900       read (*,1000,err=901,end=902) y
          goto  903
901       call  errclr ('PLTCONF')
          goto  900
902       call  endclr ('PLTCONF',y)
903       if(y .eq. 'c' .or. y .eq. 'C') then
            if(idev.eq.4) call clwclose(1,2,0)
            return
          else if(y .ne. 'y' .and. y .ne. 'Y') then
            go to 1
          end if
          if(idev.eq.4) call clwclose(1,2,0)
        end if
      end if
c
c.... mmc= no for disp,stre,...
      mmc=2
c.... color range
      cinv=-1
c.... set color palette
      call colpal(ipal,nnc)
c.... set direction for colors
      ncol = nc
      call pltcol(ipal,ipali,ncol,mmc,cinv)
c.... plot legend
      call pltftxf(vc,k1,nc,ipali,ityp,vmin,vmax)
c
      return
c
c.... formats
1000  format(a1)
2000  format('   ------ Contour Values for Plot ------'/
     +      (1x,7e11.4/1x,7e11.4))
2002  format(' Input values correct? (y or n, c = cancel) > ',$)
2003  format(' Input No. of Lines/Colors, Min. Value, Max. Value for ',
     +       ' Plot'/3x,'>',$)
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltcor(nel,ic,v,vc,nc)
c-----------------------------------------------------------------------
c
c.... Purpose: tag corners with contour number
c
c     Input:
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension ic(nel),v(nel),vc(nc)
      do 5 i = 1,nel
         ic(i) = 1
5     continue
      do 20 n = 1,nc
         do 10 i = 1,nel
            if(v(i).ge.vc(n)) then
              ic(i) = n + 1
            end if
10       continue
20    continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltcord(x,ndm,numnp,k1,k2,k3)
c-----------------------------------------------------------------------
c
c.... Purpose: plot coordinate axes with tics
c
c     Input:
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pclip
      USE pdata0
      USE pltran
      implicit double precision (a-h,o-z)
      character yy*80
      real*8 xp(3),xm(3),xb(3,8),xe(3),ddx(3),vminl(3),vmaxl(3)
      integer*4 ipl(15)
      data ipl /2,3,4,1,5,6,2,6,7,3,7,8,4,8,5/
      k4 = k2
      k2 = abs(k2)
      if(k1.eq.0) k1=5
      if(k2.eq.0) k2=5
      if(k3.eq.0) k3=5
      call pzero(xp,3)
      call pzero(xm,3)
      call plopen
c.... calculate extreme values
      call maxcor(x,ndm,numnp)
      do i = 1,3
        vminl(i) = vmin(i)
        vmaxl(i) = vmax(i)
      enddo
c.... compare with clip bounds 1=min,2=max!
      if(clip1) then
        do i = 1,3
          vminl(i) = max(vminl(i),xc(1,i))
          vmaxl(i) = min(vmaxl(i),xc(2,i))
        enddo
      end if
c.... calculate increment
      dpm = 0
      do i = 1,3
        dpm = max(dpm,0.075*(vmaxl(i)-vminl(i)))
      enddo
      ddpm = 0.5*dpm
c...  copy extreme values
      do i = 1,ndm
         xp(i) = vmaxl(i) + dpm
         xm(i) = vminl(i) - dpm
      enddo
c.... point 1
      xb(1,1)=xm(1)
      xb(2,1)=xm(2)
      xb(3,1)=xm(3)
c.... point 2
      xb(1,2)=xp(1)
      xb(2,2)=xm(2)
      xb(3,2)=xm(3)
c.... point 3
      xb(1,3)=xp(1)
      xb(2,3)=xp(2)
      xb(3,3)=xm(3)
c.... point 4
      xb(1,4)=xm(1)
      xb(2,4)=xp(2)
      xb(3,4)=xm(3)
c.... point 5
      xb(1,5)=xm(1)
      xb(2,5)=xm(2)
      xb(3,5)=xp(3)
c.... point 6
      xb(1,6)=xp(1)
      xb(2,6)=xm(2)
      xb(3,6)=xp(3)
c.... point 7
      xb(1,7)=xp(1)
      xb(2,7)=xp(2)
      xb(3,7)=xp(3)
c.... point 8
      xb(1,8)=xm(1)
      xb(2,8)=xp(2)
      xb(3,8)=xp(3)
c.... rotate values
      call plxtrn(xb,tra,vr,3,8)
c
c.... plot cubus
      if(k4.lt.0) then
        call pppcol(7)
        call plotl(xb(1,1),xb(2,1),xb(3,1),3)
        do i = 1,15
          call plotl(xb(1,ipl(i)),xb(2,ipl(i)),xb(3,ipl(i)),2)
        enddo
      end if
c
c.... plot marks on axes
      call pppcol(8)
c.... first axis
c.... make tics
      xh     = vminl(1)
      ddx(1) = (vmaxl(1)-xh)/k1
      do i = 1,k1+1
        xe(1) = xh+(i-1)*ddx(1)
        xe(2) = xm(2)
        xe(3) = vminl(3)
        call plxtrn(xe,tra,vr,3,1)
        call plotl(xe(1),xe(2),xe(3),3)
c
        xe(1) = xh+(i-1)*ddx(1)
        xe(2) = xm(2)+ddpm
        xe(3) = vminl(3)
        call plxtrn(xe,tra,vr,3,1)
        call plotl(xe(1),xe(2),xe(3),2)
      enddo
c.... make line
        call pppcol(7)
        xe(1) = xm(1)
        xe(2) = xm(2)
        xe(3) = vminl(3)
        call plxtrn(xe,tra,vr,3,1)
        call plotl(xe(1),xe(2),xe(3),3)
c
        xe(1) = xh
        xe(2) = xm(2)
        xe(3) = vminl(3)
        call plxtrn(xe,tra,vr,3,1)
        call plotl(xe(1),xe(2),xe(3),2)
c
        call pppcol(8)
        xe(1) = xh+k1*ddx(1)
        xe(2) = xm(2)
        xe(3) = vminl(3)
        call plxtrn(xe,tra,vr,3,1)
        call plotl(xe(1),xe(2),xe(3),2)
        call plabl(1)
c
c.... second axis
c.... make tics
      xh     = vminl(2)
      ddx(2) = (vmaxl(2)-xh)/k2
      do i = 1,k2+1
        xe(1) = xm(1)
        xe(2) = xh+(i-1)*ddx(2)
        xe(3) = vminl(3)
        call plxtrn(xe,tra,vr,3,1)
        call plotl(xe(1),xe(2),xe(3),3)
c
        xe(1) = xm(1)+ddpm
        xe(2) = xh+(i-1)*ddx(2)
        xe(3) = vminl(3)
        call plxtrn(xe,tra,vr,3,1)
        call plotl(xe(1),xe(2),xe(3),2)
      enddo
c.... make line
        call pppcol(7)
        xe(1) = xm(1)
        xe(2) = xm(2)
        xe(3) = vminl(3)
        call plxtrn(xe,tra,vr,3,1)
        call plotl(xe(1),xe(2),xe(3),3)
c
        xe(1) = xm(1)
        xe(2) = xh
        xe(3) = vminl(3)
        call plxtrn(xe,tra,vr,3,1)
        call plotl(xe(1),xe(2),xe(3),2)
c
        call pppcol(8)
        xe(1) = xm(1)
        xe(2) = xh+k2*ddx(2)
        xe(3) = vminl(3)
        call plxtrn(xe,tra,vr,3,1)
        call plotl(xe(1),xe(2),xe(3),2)
        call plabl(2)
c
c.... third axis
c.... make tics
      xh     = vminl(3)
      ddx(3) = (vmaxl(3)-xh)/k3
      do i = 1,k3+1
        xe(1) = xm(1)
        xe(2) = xm(2)
        xe(3) = xh+(i-1)*ddx(3)
        call plxtrn(xe,tra,vr,3,1)
        call plotl(xe(1),xe(2),xe(3),3)
c
        xe(1) = xm(1)+ddpm*0.7
        xe(2) = xm(2)+ddpm*0.7
        xe(3) = xh+(i-1)*ddx(3)
        call plxtrn(xe,tra,vr,3,1)
        call plotl(xe(1),xe(2),xe(3),2)
      enddo
c.... make line
        xe(1) = xm(1)
        xe(2) = xm(2)
        xe(3) = xh
        call plxtrn(xe,tra,vr,3,1)
        call plotl(xe(1),xe(2),xe(3),3)
c
        xe(1) = xm(1)
        xe(2) = xm(2)
        xe(3) = xh+k3*ddx(3)
        call plxtrn(xe,tra,vr,3,1)
        call plotl(xe(1),xe(2),xe(3),2)
        call plabl(3)
c
c....   plot extreme values on screen, correct only for c=cs=1!
        y0 = 0.76
        dy = 1./12.
        write(yy,'(15hGSCALE - VALUES)')
        call drawtxt(1,1.d0,y0,1,1,15,yy)
        do i = 1,3
          write(yy,2000) i, vminl(i)
          ycor = y0 - i*dy
          call drawtxt(1,1.d0,ycor,1,1,17,yy)
          write(yy,2001) i, vmaxl(i)
          ycor = y0 - (i+0.4)*dy
          call drawtxt(1,1.d0,ycor,1,1,17,yy)
        enddo
        y0 = y0-4.*dy
        dy = 1./24.
        do i = 1,3
          write(yy,2002) i,ddx(i)
          ycor = y0 - (i)*dy
          call drawtxt(1,1.d0,ycor,1,1,15,yy)
        enddo
      return
2000  format('x',i1,'min= ',e10.3)
2001  format('x',i1,'max= ',e10.3)
2002  format('ds',i1,'= ',e10.3)
      end
c
c----------------------------------------------------------------------+
c
      subroutine pltcol(ipal,ipali,nc,mmc,cinv)
c----------------------------------------------------------------------+
c
c.... Purpose: set color table
c
c     Input:
c      ipal    -
c      ipali   -
c      nc      - number of colors
c      mmc     - = 1: stre, 2: disp, 3: velo, 4: acce, 5: eigv, 6: resi
c      cinv    - =  1 colors from -=red ->+=blue  set by colo,,cinv
c                = -1 colors from -=blue->+=red
c       nc    number of colors
c
c     Output:
c       ipali  color table
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------+
      implicit double precision (a-h,o-z)
      integer*4  ipal(14),ipali(14)
c.... define color table sigma(-=red->+=blue), other (-=blue->+=red)
c.... table can be modified by colo,,cinv
      if(mmc.eq.1) ityp = +1*cinv
      if(mmc.ge.2) ityp = -1*cinv
      if(ityp.eq.1) then
          do j = 1,nc+1
            ipali(j) = ipal(j)
          enddo
      else if(ityp.eq.-1) then
          do j = 1,nc+1
            ipali(j) = ipal(nc+2-j)
          enddo
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltctx(vc,ic,nc,mmc,iev,ipali)
c-----------------------------------------------------------------------
c
c.... Purpose: plot legend for contour plots
c
c     Input:
c       vc(14)      - contour values
c       ic          - number of eg. stress to plot, e.g. stress_1
c       nc          - number of contours
c       mmc         - = 1: stre, 2: disp, 3: velo, 4: acce, 5: eigv,
c                         6: resi
c       iev         - number of eigenvector
c       ipali       - color palette
c
c     Output:
c       plot legend for contour plots
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata2
      USE pftn77
      USE plslay
      USE strnam
      implicit double precision (a-h,o-z)
      character strs(4)*13,yy*16,strph*15,eigv*4,disp*7,lin*2
cww      common /pdata9/  ifrm
      dimension vc(*)
      integer*4  ipali(14)
      data strs/' S T R E S S ' , ' DISPLACEMENT',
     +          '  VELOCITY   ' , ' ACCELERATION'/
      data eigv/'EIGV'/,disp/' DISPL '/,lin /'--'/
      iadd = 16
      ycor = 0.76
      xdv = 23.d0
cww      if(ifrm.ne.0) xdv = 40.d0
      dy = 1./xdv
      call pppcol(1)
      if(mmc.eq.1) then
c.....  stress (standard or user defined)
        if(strsus(ic).eq.'               ') then
           write(strph,'(a13,i2)' ) strs(mmc),ic
        else
           write(strph,'(a15)'    ) strsus(ic)
        end if
      else if(mmc.gt.1.and.mmc.le.4) then
c.....     displacement, velocity, acceleration
           write(strph,'(a13,i2)' ) strs(mmc),ic
      else if(mmc.eq.5) then
c.....     eigenvector
           write(strph,'(a4,i2,a7,i2)')eigv,iev,disp,ic
      else if(mmc.eq.6) then
c.....     residuum
           write(strph,'(a13,i2)' ) '   RESIDUUM  ',ic
      end if
      x = 1.
      y = ycor+0.06
      call drawtxt(1,x,y,1,1,15,strph)
      do 100 i = 1,nc
        x = 0.98
        y = ycor -i*dy
        write(yy, '(a2,i3,12x)' ) lin,i
        icol = iadd + ipali(i)
        call drawtxt(1,x,y,icol,1,16,yy)
        x = 1.06
        write(yy, '(1p,1e11.3)' ) vc(i)
        call drawtxt(1,x,y,1,1,16,yy)
100   continue
c.... draw layer number for stresses
      if(klay.ne.0.and.mmc.eq.1) then
        write(strph,'(6hLayer ,i2,5h Pos ,i2)') klay,mlay
        x = 1.
        y = 0.776
        call drawtxt(1,x,y,1,1,15,strph)
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltdam(gp,ndm,ko)
c-----------------------------------------------------------------------
c
c.... Purpose: plot state of damage at actual gauss-point/layer
c
c     Input:
c         g(3)      - Nodal coordinates of gauss point
c         ndm       - Dimension
c         ko        - 0=no damage  1= damage
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata1
      USE pdata2
      USE pltran
      implicit double precision (a-h,o-z)
      dimension gp(3),xp(3)
      dg = .002/scale
      xp(1) = gp(1)
      xp(2) = gp(2)
      z = 0.
      if(ndm.ge.3) z = gp(3)
c.... transform location of point
      call plxtrn(xp,tra,vr,ndm,1)
c.... plot damged point if ko = 1  (for composites)
      if(ko.eq.1) then
c....    plot point in x-1 x-2 plane
         call plotl(xp(1)-dg,xp(2)-dg,z,   3)
         call plotl(xp(1)-dg,xp(2)-dg,z,ipgl)
         call plotl(xp(1)+dg,xp(2)-dg,z+dg,2)
         call plotl(xp(1)+dg,xp(2)+dg,z,   2)
         call plotl(xp(1)-dg,xp(2)+dg,z-dg,2)
         call plotl(xp(1)-dg,xp(2)-dg,z,   2)
         if(ipgl.eq.1) call clpan
      end if
      return
      end
c
c----------------------------------------------------------------------+
c
      subroutine pltdraw(x,ndm,numnp,k1,k2,k3)
c----------------------------------------------------------------------+
c
c.... Purpose: plot line in color
c       (set  by colo) and linetype  (set  by line)
c       from node k2 to node k3  if k1=0 (on defm if k1<0)
c       interactively        if k1=1
c       left  button mouse:  each klick define define point
c....   right button mouse:  exit
c       with coordinates     if k1=2
c
c     Input:
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------+
      USE hptext1
      USE iwinio
      USE pdata2
      USE pftn77
      implicit double precision(a-h,o-z)
      integer*2 ihg,ivg,bstat
      real*8 xx(2,2),x(ndm,*),coor(3,2),vin(3)
      integer*4 status
      real*4 pos(3),area(6)
      integer data(20)
      if(k1.le.0) then
c....   nodes
        if(k2.gt.0.and.k2.le.numnp.and.k3.gt.0.and.k3.le.numnp) then
          x31 = 0.0d0
          x32 = 0.0d0
          if (ndm.eq.3) then
            x32 = x(3,k2)
            x33 = x(3,k3)
          end if
          call plotl(x(1,k2),x(2,k2),x32,3)
          call plotl(x(1,k3),x(2,k3),x33,2)
        end if
      else if(k1.eq.2) then
c....   ccordinates
cww     if(idev.eq.4) call clwopen('Coordinates for Line',1,iwys-170,500,200,1,2)
        if(idev.eq.4) call clwopen('Coordinates for Line',1,2)
        call pzero(coor,6)
        write(*,2001)
        call dinput(vin,3)
        coor(1,1) = vin(1)
        coor(2,1) = vin(2)
        coor(3,1) = vin(3)

        write(*,2002)
        call dinput(vin,3)
        coor(1,2) = vin(1)
        coor(2,2) = vin(2)
        coor(3,2) = vin(3)

        call plotl(coor(1,1),coor(2,1),coor(3,1),3)
        call plotl(coor(1,2),coor(2,2),coor(3,2),2)
        if(idev.eq.4) call clwclose(1,2,0)
      else if(k1.eq.1) then
c....   interactive
        if(idev.eq.1) then
c.....    inquire locator device state (15.42)
          call gpqlc(1,1,1,80,err,mode,esw,view,pos,echo,area,d1,data)
          write(*,2000)
c.....    initialize locator (new to define crosshair(type 2)) (9.30)
          call gpinlc(1,1,0,pos,2,area,0,data)
          write(*,2000)
c.....    request locator (get mouse position after click)
141       call gprqlc(1,1,status,view,pos)
          if(status.eq.0) then                    !no    button depressed
            goto 141
          else if(status.eq.1) then                !left  button depressed
c.....      nothing to do
          else if(status.eq.2) then                !right button depressed
            goto 161
          end if
          xx(1,ipt)  = (1.+ pos(1))*0.641
          xx(2,ipt)  = (1.+ pos(2))*0.641
          if(ipt.eq.1) then
            call dplot(xx(1,ipt),xx(2,ipt),3,0)
            ipt = 2
          else
            call dplot(xx(1,ipt),xx(2,ipt),2,0)
          end if
          goto 141
c.....    initialize locator (set back to standard values)
161       call gpinlc(1,1,0,pos,1,area,0,data)
        else if(idev.eq.2) then
          write(*,2000)
c.....    request locator (get mouse position after click)
142       call grqlc(1,1,status,it,xpos,ypos)
          if(status.eq.0) then                    !no    button depressed
            goto 142
          else if(status.eq.1) then                !left  button depressed
c.....      nothing to do
          else if(status.eq.2) then                !right button depressed
            return
          end if
          xx(1,ipt)  =     xpos*1.28
          xx(2,ipt)  =     ypos*1.28
          if(ipt.eq.1) then
            call dplot(xx(1,ipt),xx(2,ipt),3,0)
            ipt = 2
          else
            call dplot(xx(1,ipt),xx(2,ipt),2,0)
          end if
          goto 142
        else if(idev.eq.3) then

        else if(idev.eq.4) then
c....     set start position of cursor at midpoint
          ihg = 0
          ivg = 0
          ipt = 1
          write(*,2000)
143       call get_mouse_position(ihg,ivg,bstat)
          if(bstat.eq.0) then                     !no    button depressed
            goto 143
          else if(bstat.eq.1) then                 !left  button depressed
c.....      nothing to do
          else if(bstat.eq.2) then                 !right button depressed
            return
          end if
          xx(1,ipt)  =     ihg/(iwxgs*0.75)
          xx(2,ipt)  = 1.- ivg/(iwygs*1.00)
          if(ipt.eq.1) then
            call dplot(xx(1,ipt),xx(2,ipt),3,0)
            ipt = 2
          else
            call dplot(xx(1,ipt),xx(2,ipt),2,0)
          end if
          goto 143
        end if
      end if
      return
2000  format(' Each Left Mouse Button Click = Point,'
     +       ' Right Mouse Button Click = Exit')
2001  format(' Input Coordinates of Point 1'/3x,'>',$)
2002  format(' Input Coordinates of Point 2'/3x,'>',$)
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltefl(nel,ic,x,v,vc,nc,mmc,ipla,ipgl,idev,ipali)
c-----------------------------------------------------------------------
c
c.... Purpose: element fill for contour plots with variable no. of col (<=14)
c
c     Input:
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension ic(nel),xp(10),yp(10),zp(10),x(3,4),v(*),vc(*)
      integer*4 ipali(14)
      iadd = 16
      if(nc.gt.0) then
         vc(nc+1) = max(vc(nc),v(1))
      else
         vc(1) = v(1)
      end if
      vci = 0.0
      icx = ic(1)
      icn = ic(1)
      do 50 i = 1,nel
        vci      = max(vci,abs(v(i)))
        vc(nc+1) = max(vc(nc+1),v(i))
        icx      = max(ic(i),icx)
        icn      = min(ic(i),icn)
50    continue
      vc(nc+1)   = vc(nc+1)*1.001 + vci + 1.0d-8
      do 300 icol = icn,icx
        k = 0
        i = nel
        do 200 j = 1,nel
          if((ic(j).ge.icol .and. ic(i).le.icol)
     1                      .and. (ic(i).ne.ic(j))) then
            if(icol-1.ge.ic(i)) then
              s = (vc(icol-1)-v(i))/(v(j)-v(i))
              k = k + 1
              xp(k) = x(1,i) + (x(1,j)-x(1,i))*s
              yp(k) = x(2,i) + (x(2,j)-x(2,i))*s
              zp(k) = x(3,i) + (x(3,j)-x(3,i))*s
            end if
            s = (vc(icol)-v(i))/(v(j)-v(i))
            if(s.lt.1.0d0) then
              k = k + 1
              xp(k) = x(1,i) + (x(1,j)-x(1,i))*s
              yp(k) = x(2,i) + (x(2,j)-x(2,i))*s
              zp(k) = x(3,i) + (x(3,j)-x(3,i))*s
            end if
          else if((ic(i).ge.icol .and. ic(j).le.icol)
     1                          .and. (ic(i).ne.ic(j))) then
            s = (vc(icol)-v(i))/(v(j)-v(i))
            if(s.ge.0.0d0) then
              k = k + 1
              xp(k) = x(1,i) + (x(1,j)-x(1,i))*s
              yp(k) = x(2,i) + (x(2,j)-x(2,i))*s
              zp(k) = x(3,i) + (x(3,j)-x(3,i))*s
            end if
            if(icol-1.ge.ic(j)) then
              s = (vc(icol-1)-v(i))/(v(j)-v(i))
              k = k + 1
              xp(k) = x(1,i) + (x(1,j)-x(1,i))*s
              yp(k) = x(2,i) + (x(2,j)-x(2,i))*s
              zp(k) = x(3,i) + (x(3,j)-x(3,i))*s
            end if
          end if
          if(ic(j).eq.icol) then
            k = k + 1
            xp(k) = x(1,j)
            yp(k) = x(2,j)
            zp(k) = x(3,j)
          end if
          i = j
200     continue
c.... plot panel of this color
        ia=1
cww        if(ipla.eq.2.or.ipgl.eq.3) ia=3    ! necessary ww 29.11.2006??
        call plotl(xp(1),yp(1),zp(1),ia)
        do 150 j = 2,k
          call plotl(xp(j),yp(j),zp(j),2)
150     continue
        call pppcol(iadd+ipali(icol))
        if(ia.eq.1.or.ipgl.eq.1) call clpan
300   continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltelm(x,ie,ix,idis,scale,nie,ndm,nen,nen1,numel,k2)
c-----------------------------------------------------------------------
c
c.... Purpose: plot element labels
c
c     Input:
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical zoom
      dimension ie(nie,*),ix(nen1,*),x(ndm,*),xx(3),iplt(40),idis(*)
      if(scale.eq.0.d0) then
        call drawmess(' Scale = 0 in PLTELM',1,1)
        return
      end if
      dx1 = .005/scale
      do 100 n = 1,numel
        new=idis(n)
        if(k2.gt.0.and.k2.ne.new) goto 100
        ma = ix(nen1,new)
        if(iplma(ma).eq.0) goto 100
        xx(1) = 0.0
        xx(2) = 0.0
        xx(3) = 0.0
        jj = 0
        call pltord(ie(nie-1,ix(nen1,n)), nn,iplt)
        nn = max(1,nn-1)
        do 20 i = 1,nn
          j  = iplt(i)
          if(j.gt.0 .and. j.le.nen) then
            ii = ix(j,n)
            if(ii.gt.0) then
              jj = jj + 1
              xx(1) = xx(1) + x(1,ii)
              xx(2) = xx(2) + x(2,ii)
              if(ndm.ge.3) xx(3) = xx(3) + x(3,ii)
            end if
          end if
20      continue
c...... original version: do no plot at left corner do to zoom, 5=?
c        if(jj.gt.0) then
c          xx(1) = xx(1)/jj - dx1*5
c          xx(2) = xx(2)/jj - dx1
c          xx(3) = xx(3)/jj
c          if(zoom(xx(1),ndm,1)) then
c            call plotl(xx(1),xx(2),xx(3),3)
c            call plabl(n)
c          end if
c        end if
        if(jj.gt.0) then
          xx(1) = xx(1)/jj
          xx(2) = xx(2)/jj
          xx(3) = xx(3)/jj
          if(zoom(xx(1),ndm,1)) then
            xx(1) = xx(1) - dx1*5
            xx(2) = xx(2) - dx1
            call plotl(xx(1),xx(2),xx(3),3)
            call plabl(n)
          end if
        end if
100   continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltel1n(x,ndm,n)
c-----------------------------------------------------------------------
c
c.... Purpose: plot element with 1 node (=n)
c
c     Input:
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c.... plot  element                                    |
c-----------------------------------------------------------------------
      USE pdata1
      USE pdata2
      implicit double precision (a-h,o-z)
      dimension x(ndm,*)
c
      dx1 = .002/scale
      x3 = 0.0
      x1 = x(1,n)
      x2 = x(2,n)
      if(ndm.ge.3) x3 = x(3,n)
c
      call plotl(x1-dx1 , x2+dx1 , x3, 3)
      call plotl(x1-dx1 , x2+dx1 , x3, ipgl)
      call plotl(x1-dx1 , x2-dx1 , x3, 2)
      call plotl(x1+dx1 , x2-dx1 , x3, 2)
      call plotl(x1+dx1 , x2+dx1 , x3, 2)
      call plotl(x1-dx1 , x2+dx1 , x3, 2)
      if(ipgl.eq.1) call clpan
      call plotl(x1-12*dx1 , x2+2*dx1 , x3, 3)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine plterr(b,dr,ndf,numel)
c-----------------------------------------------------------------------
c
c.... Purpose: calculate and plot distribution of errors
c
c     Input:
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------

      USE errin1
      USE errin2
      USE errnam
      USE hdatam
      USE pdata2
      implicit double precision(a-h,o-z)
      dimension b(ndf,*), dr(*)
      logical fa
      character yy*21
      data  fa  /.false./
c.....error distribution  with respect to 'eval' percent of energy
      do i = 1,numerr
        e_om(i)  = 0.0
        e_bar(i) = eval/100.*sqrt(u_om(i)/numel)
        if(e_bar(i).eq.0.d0) e_bar(i)=1.d0
      end do
                    ypos = 0.8
      if(idev.eq.4) ypos = 0.81
      write(yy,1000)
      call drawtxt(1,1.0d0,ypos,1,1,20,yy)
      ival = eval
      xcor = 1.0
      ycor = ypos-0.03
      write(yy,1002) ival,e_name(iet(2))
      call drawtxt(1,xcor,ycor,1,1,20,yy)
c.....draw text
      xcor = 1.05
      do 10  i = 1,13
                      ypos = 0.68
        if(idev.eq.4) ypos = 0.70
        ycor = ypos - i*0.04
cww        write(yy,1001) i
        write(yy,1001) i*eval
        call drawtxt(1,xcor,ycor,1,1,20,yy)
10    continue

c.... set colors 
      call setcol(2,0)
c.....draw color boxes
      do 20  i = 1,14
        icol  = 31 -i
        xcor  = 1.00
        ycor  = 0.69 -i*0.04
c...    box
        call pppcol(icol)
        call ppbox(xcor,ycor, .025d0, .040d0,ipgl)
c...    border
        call pppcol(32)
        call ppbox(xcor,ycor, .025d0, .040d0,3)
20    continue
      hflgu  = .false.
      h3flgu = .false.
      call formfe(b,dr,dr,dr,fa,fa,fa,fa,9,1,numel,1)
c.... set colors back 
      call setcol(1,0)
1000  format('Error Analysis     ')
cww1001  format(i2,' Refinement       ')
1001  format(f5.0,'% Error        ')
1002  format(i2,'% ',a15,3x)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltflu(x,st,ndm,numnp,tra,vr,rk3)
c-----------------------------------------------------------------------
c
c.... Purpose: plot plot heat-flux vectors
c
c     Input:
c       x(ndm,numnp) - coordinates
c       st(numnp,2)  - array of fluxes to plot
c       ndm          - dimension
c       numnp        - no. of nodes
c       tr(3,3)      - rotation matrix
c       vr(3)        - coordinate shift
c       rk3          - scaling factor for lenght of fluxes
c
c     Output:
c
c     Comment:
c     heat-flux in 1-2 plane
c     storage  [ st(*,1) , st(*,2) ]
c     to be added: length of arrow tip = const. see pltfor
c     3D-plot  including rot-macros
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata1
      USE pdata4
      USE pdatas
      implicit double precision (a-h,o-z)
      dimension dd(2),xx(3,4),x(ndm,*),st(numnp,*),tra(3,3),vr(3),
     + txm(3),xcmax(2,3)
      character*20 name
c.... compute longest vector
      call pzero(xcmax,6)
      do i=1,ndm
        xcmax(1,i)= 1.e+10
        xcmax(2,i)=-1.e+10
      end do
      fm = 0.d0
      do 90 n = 1,numnp
        txm(1) = x(1,n)*fact
        txm(2) = x(2,n)*fact
        txm(3) = 0.0
        if(ndm.gt.2)txm(3) = x(3,n)*fact
c....   transform coordinates
        call plxtrn(txm,tra,vr,ndm,1)
        if(txm(1).lt.xmin(1).or.txm(1).gt.xmax(1)) go to 90
        if(txm(2).lt.xmin(2).or.txm(2).gt.xmax(2)) go to 90
        d  = st(n,1)**2 + st(n,2)**2
        do i = 1,ndm
          xcmax(1,i)= min(xcmax(1,i),st(n,i)) ! extremal values
          xcmax(2,i)= max(xcmax(2,i),st(n,i))
        end do
        fm = dmax1(fm,dsqrt(d))
90    continue
      fm = fm*scale*40.d0/rk3
c.... compute vector at each node
      do 120 n = 1,numnp
         x1 = x(1,n)*fact
         x2 = x(2,n)*fact
         txm(1) = x1
         txm(2) = x2
         txm(3) = 0.0
         if(ndm.gt.2) txm(3) = x(3,n)*fact
c....    transform coordinates
         call plxtrn(txm,tra,vr,ndm,1)
         if(txm(1).lt.xmin(1).or.txm(1).gt.xmax(1)) go to 120
         if(txm(2).lt.xmin(2).or.txm(2).gt.xmax(2)) go to 120
         dd(1) = st(n,1)/fm
         dd(2) = st(n,2)/fm
c..... modify due to quadrant for symmetry
         if(isym2(iadd+1,1,its).eq.1) dd(1) = - dd(1)
         if(isym2(iadd+1,2,its).eq.1) dd(2) = - dd(2)
c        if(isym2(iadd+1,3,its).eq.1) dd(3) = - dd(3)
c
         call pzero(xx,12)
         if(ndm.gt.2) then
           do 130 i = 1,4
              xx(3,i) = x(3,n)
130        continue
         end if
c....    compute plot locations
c....  center arrow on node
       ddx = -dd(1)*0.5d0
       ddy = -dd(2)*0.5d0
c
         xx(1,1) = x(1,n)                              + ddx
         xx(2,1) = x(2,n)                              + ddy
         xx(1,2) = xx(1,1) + dd(1)
         xx(2,2) = xx(2,1) + dd(2)
         xx(1,3) = xx(1,2) -.45d0*dd(1) - .15d0*dd(2)
         xx(2,3) = xx(2,2) -.45d0*dd(2) + .15d0*dd(1)
         xx(1,4) = xx(1,2) -.45d0*dd(1) + .15d0*dd(2)
         xx(2,4) = xx(2,2) -.45d0*dd(2) - .15d0*dd(1)
c....    transform the vector
         call plxtrn(xx,tra,vr,ndm,4)
c....    plot the vector
         call pppcol(2)
         call plotl(xx(1,1),xx(2,1),0.d0,3)
         do 110 i = 2,5
            j = i
            if(i.eq.5) j = 2
            call plotl(xx(1,j),xx(2,j),0.d0,2)
110      continue
120   continue
c...  plot extremal values
      name = '      FLUX Component'
      call plfleg(idev,name,1,ndm,xcmax)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltfor(x,f,angl,id,ndm,ndf,numnp,k2,rk3,tra,vr,is,
     +        ix,nen1,idtyp,isw)
c-----------------------------------------------------------------------
c
c.... Purpose: plot vectors for forces including transformation
c              for rot macro,  vector dof 1-3 or dof 4-6
c
c     Input:
c       isw=1 id>0
c       isw=2 id<0
c       idtyp=1 name = '  Nodal displacement'
c       idtyp=2 name = '      Nodal velocity'
c       idtyp=3 name = '  Nodal acceleration'
c       idtyp=4 name = '   Nodal eigenvector'
c       idtyp=5 name = '   Nodal eigenvector'
c       idtyp=6 name = '   Prescribed displ.'
c
c     Output:
c
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE iofile
      USE mdat2
      USE pdata1
      USE pdata2
      USE pdata7
      USE pdatas
      implicit double precision (a-h,o-z)
      logical vfl,zoom
      dimension dd(3),xx(3,4),x(ndm,*),f(ndf,*),angl(*),id(ndf,*),
     1          tra(3,3),vr(3),xz(3),xmax(2,3)
      character*20 name
      data ddx /0.6d0/, ddy /0.2d0/
      pi=datan(1.d0)*4.d0
      if(scale.eq.0.d0) then
        call drawmess(' Scale = 0 in PLTFOR',1,1)
        return
      end if
      dfm = 40.0d0
      ie = ndm
      if(ndf.eq.1) ie = 1 ! added ww
c.... compute longest vector
      do i=1,ndm
        xmax(1,i)= 1.e+10
        xmax(2,i)=-1.e+10
      end do
      fm = 0.d0
      do 90 n = 1,numnp
      if(iplmano(ix,n,nen1).eq.0)  goto 90  ! only if matn
        xz(1) = x(1,n)
        xz(2) = x(2,n)
        xz(3) = 0.d0
        if (ndm.eq.3) xz(3) = x(3,n)
        call plxtrn(xz,tra,vr,ndm,1)
        if(zoom(xz,ndm,1)) then
          d = 0.d0
          do i = 1,ie
            if(isw.eq.1) then
              if(id(i,n).gt.0) then
                xmax(1,i)= min(xmax(1,i),f(i,n)) ! extremal values
                xmax(2,i)= max(xmax(2,i),f(i,n))
                d = d + f(i,n)**2
              end if
            else if(isw.eq.2) then
              if(id(i,n).lt.0) then
                xmax(1,i)= min(xmax(1,i),f(i,n)) ! extremal values
                xmax(2,i)= max(xmax(2,i),f(i,n))
                d = d + f(i,n)**2
              end if
            end if
          end do
          fm = max(fm,d)
        end if
90    continue
      if(fm.le.0.0d0) then
        call drawmess(' No non-zero values acting on mesh',1,1)
      else if(fm.le.1.e-10) then
        return
      else
        fm = sqrt(fm)*scale*dfm/rk3
c....   compute vector at each node
        do 120 n = 1,numnp
          if(iplmano(ix,n,nen1).eq.0)  goto 120  ! only if matn
          xz(1) = x(1,n)
          xz(2) = x(2,n)
          xz(3) = 0.d0
          if (ndm.eq.3) xz(3) = x(3,n)
          call plxtrn(xz,tra,vr,ndm,1)
          if(zoom(xz,ndm,1)) then
            x1 = x(1,n)
            x2 = x(2,n)
            vfl = .false.
            call pzero(dd,3)
            if(ipla.eq.1) then
              if(isw.eq.1) then
                if((id(1,n).gt.0).and.(f(1,n).ne.0.0d0)) then
                  dd(3) = f(1,n)
                  vfl = .true.
                end if
              else if(isw.eq.2) then
                if((id(1,n).lt.0).and.(f(1,n).ne.0.0d0)) then
                  dd(3) = f(1,n)
                  vfl = .true.
                end if
              end if
            else
              ii = 0
              do 100 i = 1,ie
                ii    = ii+1
                dd(ii) = 0.d0
                if(isw.eq.1) then
                  if((id(i,n).gt.0).and.(f(i,n).ne.0.0d0)) then
                    dd(ii) = f(i,n)
                    vfl = .true.
                  end if
                else if(isw.eq.2) then
                  if((id(i,n).lt.0).and.(f(i,n).ne.0.0d0)) then
                    dd(ii) = f(i,n)
                    vfl = .true.
                  end if
                end if
100           continue
            end if
            if(vfl) then
              dd(1) = dd(1)/fm
              dd(2) = dd(2)/fm
              if(ndm.eq.3) dd(3) = dd(3)/fm
              if(angl(n).ne.0.0d0) then
                ang   = angl(n)*pi/180.d0
                cs    = cos(ang)
                sn    = sin(ang)
                ij1   = ia(1)
                ij2   = ia(2)
                tm      =  dd(ij1)*cs - dd(ij2)*sn
                dd(ij2) =  dd(ij1)*sn + dd(ij2)*cs
                dd(ij1) =  tm
              end if
c.....        modify due to quadrant for symmetry
              if(is.eq.1) then
                if(isym2(iadd+1,1,its).eq.1) dd(1) = - dd(1)
                if(isym2(iadd+1,2,its).eq.1) dd(2) = - dd(2)
                if(isym2(iadd+1,3,its).eq.1) dd(3) = - dd(3)
              end if
c....         do not multiply tip of vector
              dd1 = ddx/rk3
              dd2 = ddy/rk3
c
              x3 = 0.d0
              if(ndm.ge.3) then
                x3 = x(3,n)
                xx(3,1) = x3
                xx(3,2) = xx(3,1) + dd(3)
                xx(3,3) = xx(3,2) -dd1*dd(3) + dd2*(dd(1)+dd(2))
                xx(3,4) = xx(3,2) -dd1*dd(3) - dd2*(dd(1)+dd(2))
              end if
              xx(1,1) = x1
              xx(2,1) = x2
              xx(1,2) = xx(1,1) + dd(1)
              xx(2,2) = xx(2,1) + dd(2)
              xx(1,3) = xx(1,2) -dd1*dd(1) - dd2*(dd(2)+dd(3))
              xx(2,3) = xx(2,2) -dd1*dd(2) + dd2*(dd(1)+dd(3))
              xx(1,4) = xx(1,2) -dd1*dd(1) + dd2*(dd(2)+dd(3))
              xx(2,4) = xx(2,2) -dd1*dd(2) - dd2*(dd(1)+dd(3))
c....         transform  the vector
              call plxtrn(xx,tra,vr,ndm,4)
c....         plot the vector
              if(k2.eq.0) then
                dx1 = 0.d0
                dx2 = 0.d0
                dx3 = 0.d0
              else
                dx1 = xx(1,1) - xx(1,2)
                dx2 = xx(2,1) - xx(2,2)
                dx3 = xx(3,1) - xx(3,2)
              end if
              call plotl(xx(1,1)+dx1,xx(2,1)+dx2,xx(3,1)+dx3,3)
              call plotl(xx(1,2)+dx1,xx(2,2)+dx2,xx(3,2)+dx3,2)
              call plotl(xx(1,2)+dx1,xx(2,2)+dx2,xx(3,2)+dx3,ipgl)
              call plotl(xx(1,3)+dx1,xx(2,3)+dx2,xx(3,3)+dx3,2)
              call plotl(xx(1,4)+dx1,xx(2,4)+dx2,xx(3,4)+dx3,2)
              call plotl(xx(1,2)+dx1,xx(2,2)+dx2,xx(3,2)+dx3,2)
              if(ipgl.eq.1) call clpan
            end if
          end if
120     continue
      end if
c...  plot extremal values
      if(idtyp.eq.1) name = '  Nodal displacement'
      if(idtyp.eq.2) name = '      Nodal velocity'
      if(idtyp.eq.3) name = '  Nodal acceleration'
      if(idtyp.eq.4) name = '   Nodal eigenvector'
      if(idtyp.eq.5) name = '   Nodal eigenvector'
      if(idtyp.eq.6) name = '   Prescribed displ.'
      call plfleg(idev,name,1,ie,xmax)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltfor3(x,f,angl,id,ndm,ndf,numnp,k2,rk3,tra,vr,is,
     +        ix,nen1,idtyp)
c-----------------------------------------------------------------------
c
c.... Purpose: plot vectors for forces including transformation
c              for rot macro,  vector dof 1-3 or dof 4-6
c
c     Input:
c       isw=1 id>0
c       isw=2 id<0
c       idtyp=1 name = '  Nodal displacement'
c       idtyp=2 name = '      Nodal velocity'
c       idtyp=3 name = '  Nodal acceleration'
c       idtyp=4 name = '   Nodal eigenvector'
c       idtyp=5 name = '   Nodal eigenvector'
c       idtyp=6 name = '   Prescribed displ.'
c
c     Output:
c
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE iofile
      USE mdat2
      USE pdata1
      USE pdata2
      USE pdata7
      USE pdatas
      implicit double precision (a-h,o-z)
      logical vfl,zoom
      dimension dd(3),xx(3,4),x(ndm,*),f(ndf,*),angl(*),id(ndf,*),
     1          tra(3,3),vr(3),xz(3),xmax(2,3)
      character*20 name
      data ddx /0.6d0/, ddy /0.2d0/
      pi=datan(1.d0)*4.d0
      if(scale.eq.0.d0) then
        call drawmess(' Scale = 0 in PLTFOR3',1,1)
        return
      end if
      dfm = 40.0d0
      ie = ndm
      if(ndf.eq.1) ie = 1 ! added ww
c.... compute longest vector
      do i=1,ndm
        xmax(1,i)= 1.e+10
        xmax(2,i)=-1.e+10
      end do
      fm = 0.d0
      do 90 n = 1,numnp
      if(iplmano(ix,n,nen1).eq.0)  goto 90  ! only if matn
        xz(1) = x(1,n)
        xz(2) = x(2,n)
        xz(3) = 0.d0
        if (ndm.eq.3) xz(3) = x(3,n)
        call plxtrn(xz,tra,vr,ndm,1)
        if(zoom(xz,ndm,1)) then
          d = 0.d0
          do i = 1,ie
            if(id(i,n).le.0) then
              xmax(1,i)= min(xmax(1,i),f(i,n)) ! extremal values
              xmax(2,i)= max(xmax(2,i),f(i,n))
              d = d + f(i,n)**2
            end if
          end do
          fm = max(fm,d)
        end if
90    continue
      if(fm.le.0.0d0) then
        call drawmess(' No non-zero values acting on mesh',1,1)
      else if(fm.le.1.e-10) then
        return
      else
        fm = sqrt(fm)*scale*dfm/rk3
c....   compute vector at each node
        do 120 n = 1,numnp
          if(iplmano(ix,n,nen1).eq.0)  goto 120  ! only if matn
          xz(1) = x(1,n)
          xz(2) = x(2,n)
          xz(3) = 0.d0
          if (ndm.eq.3) xz(3) = x(3,n)
          call plxtrn(xz,tra,vr,ndm,1)
          if(zoom(xz,ndm,1)) then
            x1 = x(1,n)
            x2 = x(2,n)
            vfl = .false.
            call pzero(dd,3)
            if(ipla.eq.1) then
              if((id(1,n).le.0).and.(f(1,n).ne.0.0d0)) then
                dd(3) = f(1,n)
                vfl = .true.
              end if
            else
              ii = 0
              do 100 i = 1,ie
              ii    = ii+1
              dd(ii) = 0.d0
              if((id(i,n).le.0).and.(f(i,n).ne.0.0d0)) then
                dd(ii) = f(i,n)
                vfl = .true.
              end if
100           continue
            end if
            if(vfl) then
              dd(1) = dd(1)/fm
              dd(2) = dd(2)/fm
              if(ndm.eq.3) dd(3) = dd(3)/fm
              if(angl(n).ne.0.0d0) then
                ang   = angl(n)*pi/180.d0
                cs    = cos(ang)
                sn    = sin(ang)
                ij1   = ia(1)
                ij2   = ia(2)
                tm      =  dd(ij1)*cs - dd(ij2)*sn
                dd(ij2) =  dd(ij1)*sn + dd(ij2)*cs
                dd(ij1) =  tm
              end if
c.....        modify due to quadrant for symmetry
              if(is.eq.1) then
                if(isym2(iadd+1,1,its).eq.1) dd(1) = - dd(1)
                if(isym2(iadd+1,2,its).eq.1) dd(2) = - dd(2)
                if(isym2(iadd+1,3,its).eq.1) dd(3) = - dd(3)
              end if
c....         do not multiply tip of vector
              dd1 = ddx/rk3
              dd2 = ddy/rk3
c
              x3 = 0.d0
              if(ndm.ge.3) then
                x3 = x(3,n)
                xx(3,1) = x3
                xx(3,2) = xx(3,1) + dd(3)
                xx(3,3) = xx(3,2) -dd1*dd(3) + dd2*(dd(1)+dd(2))
                xx(3,4) = xx(3,2) -dd1*dd(3) - dd2*(dd(1)+dd(2))
              end if
              xx(1,1) = x1
              xx(2,1) = x2
              xx(1,2) = xx(1,1) + dd(1)
              xx(2,2) = xx(2,1) + dd(2)
              xx(1,3) = xx(1,2) -dd1*dd(1) - dd2*(dd(2)+dd(3))
              xx(2,3) = xx(2,2) -dd1*dd(2) + dd2*(dd(1)+dd(3))
              xx(1,4) = xx(1,2) -dd1*dd(1) + dd2*(dd(2)+dd(3))
              xx(2,4) = xx(2,2) -dd1*dd(2) - dd2*(dd(1)+dd(3))
c....         transform  the vector
              call plxtrn(xx,tra,vr,ndm,4)
c....         plot the vector
              if(k2.eq.0) then
                dx1 = 0.d0
                dx2 = 0.d0
                dx3 = 0.d0
              else
                dx1 = xx(1,1) - xx(1,2)
                dx2 = xx(2,1) - xx(2,2)
                dx3 = xx(3,1) - xx(3,2)
              end if
              call plotl(xx(1,1)+dx1,xx(2,1)+dx2,xx(3,1)+dx3,3)
              call plotl(xx(1,2)+dx1,xx(2,2)+dx2,xx(3,2)+dx3,2)
              call plotl(xx(1,2)+dx1,xx(2,2)+dx2,xx(3,2)+dx3,ipgl)
              call plotl(xx(1,3)+dx1,xx(2,3)+dx2,xx(3,3)+dx3,2)
              call plotl(xx(1,4)+dx1,xx(2,4)+dx2,xx(3,4)+dx3,2)
              call plotl(xx(1,2)+dx1,xx(2,2)+dx2,xx(3,2)+dx3,2)
              if(ipgl.eq.1) call clpan
            end if
          end if
120     continue
      end if
c...  plot extremal values
      if(idtyp.eq.1) name = '  Presc displacement'
      call plfleg(idev,name,1,ie,xmax)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltfor2(x,f,angl,id,ndmg,ndf,numnp,k2,rk3,tra,vr,is,k1,
     +        ix,nen1,idtyp)
c-----------------------------------------------------------------------
c
c.... Purpose: plot vectors for forces including transformation
c              for rot macro,  vector dof 1-3 or dof 4-6
c              separately for each dof
c     Input:
c        idtyp =1 load
c              =2 reac
c
c     Output:
c
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE iofile
      USE mdat2
      USE pdata1
      USE pdata2
      USE pdata7
      USE pdatas
      implicit double precision (a-h,o-z)
      logical vfl,zoom
      dimension dd(3),xx(3,7),x(ndmg,*),f(ndf,numnp),angl(*),id(ndf,*),
     1          tra(3,3),vr(3),xz(3),xmax(2,3)
      character*20 name
      data ddx /0.6d0/, ddy /0.2d0/
cwd   limit ndm to 3 as in IGA 4D-coordinates prevail
      ndm=min(ndmg,3)
      pi=datan(1.d0)*4.d0
      ictyp=1 ! typ color 1=load,2=moment
      if(scale.eq.0.d0) then
        call drawmess(' Scale = 0 in PLTFOR',1,1)
        return
      end if
      dfm = 40.0d0
c...  define dofs ia-ie for maximum calculation and dofs nca-nce to plot, nmax=no. of values to plot
      if((ndm.eq.2.and.ndf.eq.1)) then
c...    special cases(ndm,ndf): 2d-truss without transformation(2,1)
          ian = 1
          ie  = 1
          nca = 1
          nce = 1
          idof= 0
          nmax= 1
          ictyp=1
      else if(ipla.eq.2) then   ! (2,3) is also 2d-plate
c...    special cases(ndm,ndf): 2d-beam(2,3)
        if(k1.eq.0) then  ! all forces
          ian = 1
          ie  = 2
          nca = 1
          nce = 2
          idof= 0
          nmax= 2
          ictyp=1
        else if (k1.le.ndf) then  ! one value
cww          ian = 1 ! force
cww          ie  = 2
          ian = k1 ! force
          ie  = k1
          nca = k1
          nce = k1
          idof= 0
          nmax= 2
          ictyp=1
          if(k1.eq.3) then ! moment
            ian = 3
            ie  = 3
            idof= 2 !  not correct,  vector shows in 3-dir but this is not defined!
            nmax= 1
            ictyp=2
          end if
        else if (k1.eq.ndf+1) then  ! all forces
          ian = 1
          ie  = 2
          nca = 1
          nce = 2
          nmax= 2
          idof= 0
          ictyp=1
        else if (k1.eq.ndf+2) then  ! all moments
          ian = 3
          ie  = 3
          nca = 3
          nce = 3
          idof= 2 ! see above
          nmax= 1
          ictyp=2
        end if
      else if((ndm.eq.2.and.ndf.eq.3).or.ipla.eq.1) then
c...    special cases(ndm,ndf): 2d-plate(2,3),3d-plate(3,3)
        if(k1.eq.0) then  ! all forces
          ian = 1
          ie  = 1
          nca = 1
          nce = 1
          idof= 0 ! not correct 2D, vector must show in 3dir, but this is not defined!
          if(ipla.eq.1) idof= -2 ! 3D
          nmax= 1
          ictyp=1
        else if (k1.le.ndf) then  ! one value
cww          ian = 1 ! force
cww          ie  = 1
          ian = k1 ! force
          ie  = k1
          nca = k1
          nce = k1
          idof= 0  ! see above, 2D
          if(ipla.eq.1) idof= -2 ! 3D
          nmax= 1
          ictyp=1
          if(k1.eq.2.or.k1.eq.3) then ! moment
            ian = 2
            ie  = 3
            idof= 1
            nmax= 2
            ictyp=2
          end if
        else if (k1.eq.ndf+1) then  ! all forces
          ian = 1
          ie  = 1
          nca = 1
          nce = 1
          idof= 0  ! see above, 2D
          if(ipla.eq.1) idof= -2 ! 3D
          nmax= 1
          ictyp=1
        else if (k1.eq.ndf+2) then  ! all moments
          ian = 2
          ie  = 3
          nca = 2
          nce = 3
          idof= 1
          nmax= 2
          ictyp=2
        end if
      else
c...    standard cases(ndm,ndf): 2d-truss(2,2),3d-truss(3,3),plane stress(2,2),brick(3,3),
c                                3d-beam(3,6),shell(3,5),shell(3,6) and 3d-beam(3,7),shell(3,7)
        if(k1.eq.0) then  ! all forces
          ian = 1
          ie  = ndm
          nca = 1
          nce = ndm
          idof= 0
          nmax= ndm
          ictyp=1
        else if (k1.le.ndf) then  ! one value
cww          ian = 1 ! force  ... scale on all values
cww          ie  = ndm
          ian = k1 ! force    ... scale only on value k1
          ie  = k1
          nca = k1
          nce = k1
          idof= 0
          nmax= ndm
          ictyp=1
          if(k1.gt.ndm) then ! moment 4-5/6
            ian = 4
            ie  = min(ndf,6) ! 5 or 6
            idof= 3
            nmax= min(ndf,6)-3
            ictyp=2
          else if(k1.gt.6.and.ndf.gt.6) then ! higher order dof, like force
            ictyp=3
          end if
        else if (k1.eq.ndf+1) then  ! all forces
          ian = 1
          ie  = ndm
          nca = 1
          nce = ndm
          idof= 0
          nmax= ndm
          ictyp=1
        else if (k1.eq.ndf+2) then  ! all moments 4-5/6
          ian = 4
          ie  = min(ndf,6)
          nca = 4
          nce = min(ndf,6)
          idof= 3
          nmax= min(ndf,6)-3
          ictyp=2
        end if
      end if
c.... define color
      if(ictyp.eq.1) call pppcol(2)
      if(ictyp.eq.2) call pppcol(8)
      if(ictyp.eq.3) call pppcol(8)
c.... compute longest vector
      call pzero(xmax,6)
      do i=1,nmax
        xmax(1,i)= 1.e+10
        xmax(2,i)=-1.e+10
      end do
      fm = 0.
      do 90 n = 1,numnp
      if(iplmano(ix,n,nen1).eq.0)  goto 90  ! only if matn
        xz(1) = x(1,n)
        xz(2) = x(2,n)
        xz(3) = 0.d0
        if (ndm.eq.3) xz(3) = x(3,n)
        call plxtrn(xz,tra,vr,ndm,1)
        if(zoom(xz,ndm,1)) then
          d  = 0.d0
          if(idtyp.eq.1) then !load
            ii = 0
            do i = ian,ie
              if(id(i,n).gt.0) then
                ii = ii+1
                xmax(1,ii)= min(xmax(1,ii),f(i,n)) ! extremal values
                xmax(2,ii)= max(xmax(2,ii),f(i,n))
                d = d + f(i,n)**2
              end if
            end do
          else ! reac
            ii = 0
            do i = ian,ie
              ii = ii+1
              xmax(1,ii)= min(xmax(1,ii),f(i,n)) ! extremal values
              xmax(2,ii)= max(xmax(2,ii),f(i,n))
              d = d + f(i,n)**2
            end do
          end if
          fm = max(fm,d)
        end if
90    continue
      if(fm.le.0.0d0) then
        call drawmess(' No non-zero values acting on mesh',1,1)
      else if(fm.le.1.e-10) then
        return
      else
        fm = sqrt(fm)*scale*dfm/rk3
c....   compute vector components at each node
        do 120 n = 1,numnp ! loop over nodes
          if(iplmano(ix,n,nen1).eq.0)  goto 120  ! only if matn
c....     point
          xz(1) = x(1,n)
          xz(2) = x(2,n)
          xz(3) = 0.d0
          if (ndm.eq.3) xz(3) = x(3,n)
          call plxtrn(xz,tra,vr,ndm,1)
          if(zoom(xz,ndm,1)) then
            x1 = x(1,n)
            x2 = x(2,n)
c....       loop over components at node
            nci = 0
            do 121 nc = nca,nce
              vfl = .false.
              call pzero(dd,3)
              if(idtyp.eq.1) then ! find load
                if((id(nc,n).gt.0).and.(f(nc,n).ne.0.0d0)) then
                  vfl = .true.
                  dd(nc-idof) = f(nc,n)
                end if
              else ! find reac
                if(f(nc,n).ne.0.0d0) then
                  vfl = .true.
                  dd(nc-idof) = f(nc,n)
                end if
              end if
c....         plot
              if(vfl) then
                dd(1) = dd(1)/fm
                dd(2) = dd(2)/fm
                if(ndm.eq.3) dd(3) = dd(3)/fm
                if(angl(n).ne.0.0d0) then ! in case of angl
                  ang   = angl(n)*pi/180.d0
                  cs    = cos(ang)
                  sn    = sin(ang)
                  ij1   = ia(1)
                  ij2   = ia(2)
                  tm      =  dd(ij1)*cs - dd(ij2)*sn
                  dd(ij2) =  dd(ij1)*sn + dd(ij2)*cs
                  dd(ij1) =  tm
                end if
c.....          modify due to quadrant for symmetry
                if(is.eq.1) then
                  if(isym2(iadd+1,1,its).eq.1) dd(1) = - dd(1)
                  if(isym2(iadd+1,2,its).eq.1) dd(2) = - dd(2)
                  if(isym2(iadd+1,3,its).eq.1) dd(3) = - dd(3)
                end if
c....           do not multiply tip of vector
                dd1 = ddx/rk3
                dd2 = ddy/rk3
c
                x3 = 0.d0
                call pzero(xx,21)
                if(ndm.eq.3) then
                  x3 = x(3,n)
                  xx(3,1) = x3
                  xx(3,5) = xx(3,1) + dd(3)*1.5d0
                  xx(3,6) = xx(3,5) -dd1*dd(3) + dd2*(dd(1)+dd(2))
                  xx(3,7) = xx(3,5) -dd1*dd(3) - dd2*(dd(1)+dd(2))
                  xx(3,2) =(xx(3,6) + xx(3,7))*0.5d0
                  xx(3,3) = xx(3,2) -dd1*dd(3) + dd2*(dd(1)+dd(2))
                  xx(3,4) = xx(3,2) -dd1*dd(3) - dd2*(dd(1)+dd(2))
                end if
                xx(1,1) = x1
                xx(2,1) = x2
                xx(1,5) = xx(1,1) + dd(1)*1.5d0
                xx(2,5) = xx(2,1) + dd(2)*1.5d0
                xx(1,6) = xx(1,5) -dd1*dd(1) - dd2*(dd(2)+dd(3))
                xx(2,6) = xx(2,5) -dd1*dd(2) + dd2*(dd(1)+dd(3))
                xx(1,7) = xx(1,5) -dd1*dd(1) + dd2*(dd(2)+dd(3))
                xx(2,7) = xx(2,5) -dd1*dd(2) - dd2*(dd(1)+dd(3))
                xx(1,2) =(xx(1,6) + xx(1,7))*0.5d0
                xx(2,2) =(xx(2,6) + xx(2,7))*0.5d0
                xx(1,3) = xx(1,2) -dd1*dd(1) - dd2*(dd(2)+dd(3))
                xx(2,3) = xx(2,2) -dd1*dd(2) + dd2*(dd(1)+dd(3))
                xx(1,4) = xx(1,2) -dd1*dd(1) + dd2*(dd(2)+dd(3))
                xx(2,4) = xx(2,2) -dd1*dd(2) - dd2*(dd(1)+dd(3))
c
c....           transform  the vector
                call plxtrn(xx,tra,vr,ndm,7)
                if(k2.eq.0) then
                  dx1 = 0.d0
                  dx2 = 0.d0
                  dx3 = 0.d0
                else
                  dx1 = xx(1,1) - xx(1,5)
                  dx2 = xx(2,1) - xx(2,5)
                  dx3 = xx(3,1) - xx(3,5)
                end if
c....           plot the vector
c                     (3)   6
c          1 --------------(2)---5
c                     (4)   7
                call plotl(xx(1,1)+dx1,xx(2,1)+dx2,xx(3,1)+dx3,3) ! line
                call plotl(xx(1,5)+dx1,xx(2,5)+dx2,xx(3,5)+dx3,2)
c
                call plotl(xx(1,5)+dx1,xx(2,5)+dx2,xx(3,5)+dx3,ipgl) ! first tip
                call plotl(xx(1,6)+dx1,xx(2,6)+dx2,xx(3,6)+dx3,2)
                call plotl(xx(1,7)+dx1,xx(2,7)+dx2,xx(3,7)+dx3,2)
                call plotl(xx(1,5)+dx1,xx(2,5)+dx2,xx(3,5)+dx3,2)
                if(ipgl.eq.1) call clpan
c
                if(ictyp.eq.2) then
                 call plotl(xx(1,2)+dx1,xx(2,2)+dx2,xx(3,2)+dx3,ipgl) ! second tip
                 call plotl(xx(1,3)+dx1,xx(2,3)+dx2,xx(3,3)+dx3,2)
                 call plotl(xx(1,4)+dx1,xx(2,4)+dx2,xx(3,4)+dx3,2)
                 call plotl(xx(1,2)+dx1,xx(2,2)+dx2,xx(3,2)+dx3,2)
                 if(ipgl.eq.1) call clpan
                end if
              end if !vfl
121         continue
          end if !zoom
120       continue ! loop over nodes
      end if ! no load
c...  plot extremal values (disturbing in combination with stre etc.)
cww      if(idtyp.eq.1) name = '          Nodal Load'
cww      if(idtyp.eq.2) name = '      Nodal Reaction'
cww      if(idtyp.eq.3) name = '      Material Force'
cww   call plfleg(idev,name,ian,nmax,xmax)
      return
      end
c
c----------------------------------------------------------------------+
c
      subroutine pltfor1(x,ix,dir,xdir,knode,numnp,numel,nen,nen1,ndml,
     +                   xm,i1,i2,tra,vr,ipgl)
c----------------------------------------------------------------------+
c
c.... Purpose: plot vectors for basis, including rot-macro
c
c     Input:
c        x     coordinates
c        x0    position
c        dir   vector of basis
c        diri  local basis
c        xm    scaling factor
c
c     Output:
c
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------+
      USE pdata1
      implicit double precision (a-h,o-z)
      logical zoom
      dimension x(ndml,*), dir(10,knode,2),x0(3),diri(3,3),xdir(4),
     +          xx(3,5),tra(3,3),vr(3),x0z(3),ix(nen1,*)
        xm = xm/scale/40.
      idtyp = xdir(1)
cwd   restraining ndm locally to 3 for printing
      ndm=min(ndml,3)
      if(idtyp.eq.0.or.idtyp.gt.1) then
c....   basis at global nodes
        ia = 1
        ie = numnp
        if(i1.ge.1.and.i1.le.numnp) ia = i1
        if(i2.ge.1.and.i2.le.numnp) ie = i2
        do 100 ii = ia,ie
          if(iplmano(ix,n,nen1).eq.0)  goto 100  ! only if matn
          x0(1)  = x(1,ii)
          x0(2)  = x(2,ii)
          x0(3)  = x(3,ii)
          x0z(1) = x0(1)
          x0z(2) = x0(2)
          x0z(3) = x0(3)
            call plxtrn(x0z,tra,vr,3,1)
            if(zoom(x0z,3,1)) then
c.....      get the basis vectors
            call pzero(diri,9)
            do k = 1,3
              diri(k,1) = dir(k  ,ii,1)
              diri(k,2) = dir(k+3,ii,1)
              diri(k,3) = dir(k+6,ii,1)
            enddo
c.....      perspective projecton of axes
            do m = 1,ndm
              call pppcol(m)
              call pzero(xx,15)
              do n = 1,ndm
                fac1    = diri(n,m)*xm
                xx(n,1) = x0(n)
                xx(n,2) = xx(n,1) + fac1
                xx(n,5) = xx(n,2)
              enddo
              fac1 = diri(1,m)*xm
              fac2 = diri(2,m)*xm
              fac3 = diri(3,m)*xm
              xx(1,3) = xx(1,2) -.3*fac1 - .1*(fac2+fac3)
              xx(2,3) = xx(2,2) -.3*fac2 + .1*(fac1+fac3)
              xx(3,3) = xx(3,2) -.3*fac3 + .1*(fac1+fac2)
              xx(1,4) = xx(1,2) -.3*fac1 + .1*(fac2+fac3)
              xx(2,4) = xx(2,2) -.3*fac2 - .1*(fac1+fac3)
              xx(3,4) = xx(3,2) -.3*fac3 - .1*(fac1+fac2)
c.....        transform vector if necessary
                call plxtrn(xx,tra,vr,ndm,5)
c.....        plot the vector
              call plotl(xx(1,1),xx(2,1),xx(3,1),3)
              call plotl(xx(1,2),xx(2,2),xx(3,2),2)
              call plotl(xx(1,2),xx(2,2),xx(3,2),ipgl)
              call plotl(xx(1,3),xx(2,3),xx(3,3),2)
              call plotl(xx(1,4),xx(2,4),xx(3,4),2)
              call plotl(xx(1,2),xx(2,2),xx(3,2),2)
              if(ipgl.eq.1) call clpan
              call plotl(xx(1,2),xx(2,2),xx(3,2),3)
              call plabl(m)
            enddo
          end if
100     continue
      else if(idtyp.eq.1) then
c....   basis at element nodes
        ia = 1
        ie = numel
        if(i1.ge.1.and.i1.le.numel) ia = i1
        if(i2.ge.1.and.i2.le.numel) ie = i2
        do 200 ii = ia,ie
           if(iplma(ix(nen1,ii)).eq.0)  goto 200 ! only if MATN
         do kk = 1,nen
          jj = ix(kk,ii)
          x0(1)  = x(1,jj)
          x0(2)  = x(2,jj)
          x0(3)  = x(3,jj)
          x0z(1) = x0(1)
          x0z(2) = x0(2)
          x0z(3) = x0(3)
          call plxtrn(x0z,tra,vr,3,1)
          if(zoom(x0z,3,1)) then
c.......    get the basis vectors
            call pzero(diri,9)
            ll = (ii-1)*nen + kk
            do k = 1,3
              diri(k,1) = dir(k  ,ll,1)
              diri(k,2) = dir(k+3,ll,1)
              diri(k,3) = dir(k+6,ll,1)
            enddo
c.......    perspective projecton of axes
            do m = 1,ndm
              call pppcol(m)
              call pzero(xx,15)
              do n = 1,ndm
                fac1    = diri(n,m)*xm
                xx(n,1) = x0(n)
                xx(n,2) = xx(n,1) + fac1
                xx(n,5) = xx(n,2)
              enddo
              fac1 = diri(1,m)*xm
              fac2 = diri(2,m)*xm
              fac3 = diri(3,m)*xm
              xx(1,3) = xx(1,2) -.3*fac1 - .1*(fac2+fac3)
              xx(2,3) = xx(2,2) -.3*fac2 + .1*(fac1+fac3)
              xx(3,3) = xx(3,2) -.3*fac3 + .1*(fac1+fac2)
              xx(1,4) = xx(1,2) -.3*fac1 + .1*(fac2+fac3)
              xx(2,4) = xx(2,2) -.3*fac2 - .1*(fac1+fac3)
              xx(3,4) = xx(3,2) -.3*fac3 - .1*(fac1+fac2)
c.......      transform vector if necessary
              call plxtrn(xx,tra,vr,ndm,5)
c.......      plot the vector
              call plotl(xx(1,1),xx(2,1),xx(3,1),3)
              call plotl(xx(1,2),xx(2,2),xx(3,2),2)
              call plotl(xx(1,2),xx(2,2),xx(3,2),ipgl)
              call plotl(xx(1,3),xx(2,3),xx(3,3),2)
              call plotl(xx(1,4),xx(2,4),xx(3,4),2)
              call plotl(xx(1,2),xx(2,2),xx(3,2),2)
              if(ipgl.eq.1) call clpan
              call plotl(xx(1,2),xx(2,2),xx(3,2),3)
              call plabl(m)
            enddo
          end if
         enddo
200     continue
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltftx(vc,ic,nc,mmc,iev,ipali,vmin,vmax)
c-----------------------------------------------------------------------
c
c.... Purpose: set legend for filled plots
c
c     Input:
c       vc(14)      - contour values
c       ic          - number of eg. stress to plot, e.g. stress_1
c       nc          - number of contour values
c       mmc         - = 1: stre, 2: disp, 3: velo, 4: acce, 5: eigv,
c                         6: resi
c       iev         - number of eigenvector
c       ipali       - color palette
c       vmin,vmax   - min/max values of contour
c
c     Output:
c       plot legend for filled plots
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata2
      USE pftn77
      USE plslay
      USE rpdata
      USE strnam
      implicit double precision (a-h,o-z)
      character strs(4)*13,strph*15,eigv*4,disp*7
cww      common /pdata9/ ifrm
      dimension vc(14)
      dimension yfr(4)
      integer*4  ipali(14)
      data strs/' S T R E S S ' , ' DISPLACEMENT',
     +          '  VELOCITIY  ' , ' ACCELERATION'/
      data eigv/'EIGV'/,disp/' DISPL '/
cww      data yfr/0.90d0,0.70d0,0.50d0,0.30d0/
      iadd = 16
      xdv = 23.d0
cww      if(ifrm.ne.0) xdv = 40.d0
      dy = 1./xdv
      dyb = dy
      ycor  = 0.75
      ycor1 = 0.76
      if(idev.eq.4) ycor1 = 0.735
cww      if(ifrm.ne.0) then
cww        ycor  = yfr(ifrm)
cww        ycor1 = yfr(ifrm) - 0.01
cww      end if
c.... name
      if(mmc.eq.1) then
c.....  stress (standard or user defined)
        if(strsus(ic).eq.'               ') then
           write(strph,'(a13,i2)' ) strs(mmc),ic
        else
           write(strph,'(a15)'    ) strsus(ic)
        end if
      else if(mmc.gt.1.and.mmc.le.4) then
c.....  displacement, velocity, acceleration
        write(strph,'(a13,i2)' ) strs(mmc),ic
      else if(mmc.eq.5) then
c.....  eigenvector
        write(strph,'(a4,i2,a7,i2)')eigv,iev,disp,ic
      else if(mmc.eq.6) then
c.....  residuum
        write(strph,'(a13,i2)' ) '   RESIDUUM  ',ic
      end if
c
c.... minimum value
c
      x = 1.
      y = ycor+0.048
      if(idev.eq.4) y = ycor+0.06
      call drawtxt(1,x,y,1,1,15,strph)
      if(rmn.le.vmin) then
        write(strph, '(1p1e10.3,4h min)' ) rmn
      else
        write(strph, '(1p1e10.3)' ) vmin
      end if
      x = 1.04
      y = ycor
      call drawtxt(1,x,y,1,1,15,strph)
c.... box
      call pppcol(iadd+ipali(1))
      call ppbox(1.0d0,ycor1-dy,0.025d0,dyb,ipgl)
c.... border
      call pppcol(32)
      call ppbox(1.0d0,ycor1-dy,0.025d0,dyb,3)
c.... text
      if(nc.gt.0) then
        write(strph, '(1p1e10.3)') vc(1)
        y = ycor - dy
        call drawtxt(1,x,y,1,1,11,strph)
      end if
c
c.... interior values
c
      do 210 i = 2,nc
c....   box
        call pppcol(iadd+ipali(i))
        call ppbox(1.00d0,ycor1-i/xdv,0.025d0,dyb,ipgl)
c....   border
        call pppcol(32)
        call ppbox(1.00d0,ycor1-i/xdv,0.025d0,dyb,3)
c....   text
        write(strph, '(1p,1e10.3)') vc(i)
        y = ycor - i/xdv
        call drawtxt(1,x,y,1,1,11,strph)
210   continue
c
c.... maximum value
c
c.... box
      call pppcol(iadd+ipali(nc+1))
      call ppbox(1.00d0,ycor1-(nc+1)/xdv,0.025d0,dyb,ipgl)
c.... border
      call pppcol(32)
      call ppbox(1.00d0,ycor1-(nc+1)/xdv,0.025d0,dyb,3)
c.... text
      if(rmx.ge.vmax) then
        write(strph, '(1p1e10.3,4h max)' ) rmx
      else
        write(strph, '(1p1e10.3)' ) vmax
      end if
      y = ycor - (nc+1)/xdv
      call drawtxt(1,x,y,1,1,15,strph)

c.... draw layer number for stresses
      if(klay.ne.0.and.mmc.eq.1) then
        write(strph,'(6hLayer ,i2,5h Pos ,i2)') klay,mlay
        x = 1.
        y = 0.776
        call drawtxt(1,x,y,1,1,15,strph)
      end if

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltftxf(vc,ic,nc,ipali,ityp,vmin,vmax)
c-----------------------------------------------------------------------
c
c.... Purpose: set legend for filled plots for forces in pltconf
c              similar to pltftx
c
c     Input:
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE fornam
      USE pdata2
      USE pftn77
      USE plslay
      USE rpdata
      USE strnam
      implicit double precision (a-h,o-z)
      character strph*15
cww      common /pdata9/  ifrm
      dimension vc(14)
      dimension yfr(4)
      integer*4  ipali(14)
cww      data yfr/0.90d0,0.70d0,0.50d0,0.30d0/
      iadd = 16
      xdv = 23.d0
cww      if(ifrm.ne.0) xdv = 40.d0
      dy = 1./xdv
      dyb = dy
      ycor  = 0.75
      ycor1 = 0.76
      if(idev.eq.4) ycor1 = 0.735
cww      if(ifrm.ne.0) then
cww        ycor  = yfr(ifrm)
cww        ycor1 = yfr(ifrm) - 0.01
cww      end if
c.....name (user defined)
      if(ityp.eq.1) write(strph,'(a15)'    ) forsus(ic) ! from forc
      if(ityp.eq.2) write(strph,'(a15)'    ) strsus(ic) ! from str1
      x = 1.
      y = ycor+0.048
      if(idev.eq.4) y = ycor+0.06
      call drawtxt(1,x,y,1,1,15,strph)
      if(rmn.le.vmin) then
        write(strph, '(1p1e10.3,4h min)' ) rmn
      else
        write(strph, '(1p1e10.3)' ) vmin
      end if
      x = 1.04
      y = ycor
      call drawtxt(1,x,y,1,1,15,strph)
c...  box
      call pppcol(iadd+ipali(1))
      call ppbox(1.0d0,ycor1-dy,0.025d0,dyb,ipgl)
c...  border
      call pppcol(32)
      call ppbox(1.0d0,ycor1-dy,0.025d0,dyb,3)
      if(nc.gt.0) then
        write(strph, '(1p1e10.3)') vc(1)
        y = ycor - dy
        call drawtxt(1,x,y,1,1,11,strph)
      end if
      do 210 i = 2,nc
c...    box
        call pppcol(iadd+ipali(i))
        call ppbox(1.00d0,ycor1-i/xdv,0.025d0,dyb,ipgl)
c...    border
        call pppcol(32)
        call ppbox(1.00d0,ycor1-i/xdv,0.025d0,dyb,3)
        write(strph, '(1p,1e10.3)') vc(i)
        y = ycor - i/xdv
        call drawtxt(1,x,y,1,1,11,strph)
210   continue
c...  box
      call pppcol(iadd+ipali(nc+1))
      call ppbox(1.00d0,ycor1-(nc+1)/xdv,0.025d0,dyb,ipgl)
c...  border
      call pppcol(32)
      call ppbox(1.00d0,ycor1-(nc+1)/xdv,0.025d0,dyb,3)
      if(rmx.ge.vmax) then
        write(strph, '(1p1e10.3,4h max)' ) rmx
      else
        write(strph, '(1p1e10.3)' ) vmax
      end if
      y = ycor - (nc+1)/xdv
      call drawtxt(1,x,y,1,1,15,strph)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine plthid1(ix,x,xl,idis,dist,ndm,nen,nen1,numel,k1)
c-----------------------------------------------------------------------
c
c.... Purpose: sort elements to longest distance
c              midpoint of element - viepoint
c
c     Input:
c       vwpt: coordinates of viewpoint
c
c     Output:
c       idis(numel)
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE vpoint
      implicit double precision(a-h,o-z)
      dimension ix(nen1,*),x(ndm,*),xl(ndm,*),xc(3),dist(*),idis(*)
      if(k1.eq.0) then
c....   reset field idis to initial values 1-numel
        do n = 1,numel
          idis(n) = n
        enddo
      else if(k1.ne.0) then
c...    sort elements due to distance from viewpoint
        do 105 n = 1,numel
          call pzero(xc,3)
          np = nen           ! no. of nodes for bricks = 4 = only face
          if(nen.gt.4) np=4  ! could deleted
          npp = 0
          do 104 i = 1,np
            ii = abs(ix(i,n))
            if(ii.le.0) goto 104
            npp = npp+1
c....       add to midpoint
            do j=1,ndm
              xc(j) = xc(j)+x(j,ii)
            enddo
c
104       continue
          do j=1,ndm
             xc(j) = xc(j)/npp
          enddo
c
c....     calculate distance to viewpoint
          d1  = xc(1)-vwpt(1)
          d2  = xc(2)-vwpt(2)
          d3  = xc(3)-vwpt(3)
          dis = sqrt(d1*d1+d2*d2+d3*d3)
          dist(n) = dis
          idis(n) = n
105     continue
c....   recalculate field idis due to distance to viewpoint
        call plthid2(dist,idis,numel)
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine plthid2(dist,idis,numel)
c-----------------------------------------------------------------------
c
c.... Purpose: sort elements to longest distance from viewpoint
c
c     Input:
c
c     Output:
c
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension dist(*),idis(*)
c.... sort on first value
      do n = 1,numel-1
        do i = n+1,numel
          if(dist(i).gt.dist(n)) then
            sn = dist(n)
            dist(n) = dist(i)
            dist(i) = sn
            in = idis(n)
            idis(n) = idis(i)
            idis(i) = in
          end if
        enddo
      enddo
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine plthid3(idis,numel)
c-----------------------------------------------------------------------
c
c.... Purpose: initial field of element numbers-distance from viewpoint
c
c     Input:
c
c     Output:
c       idis
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension idis(*)
      do n = 1,numel
         idis(n) = n
      enddo
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine plthid4
c-----------------------------------------------------------------------
c
c.... Purpose: input of persp. parameters for hids macro in case of .not. pers
c
c     Input:
c
c     Output:
c
c     Comment:
c     can the viewpoint calc. automatically from the transformation ?
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE iofile
      USE iwinio
      USE pdata0
      USE pdata2
      USE ppers
      USE vpoint
      implicit double precision (a-h,o-z)

cww   if(idev.eq.4) call clwopen('Perspective  Parameters',1,iwys-260,500,300,1,2)
      if(idev.eq.4) call clwopen('Perspective  Parameters',1,2)

cww   if(ior.lt.0) write(*,2000)
cww   call dinput(vwpt,1)
cww   iused = vwpt(1)
cww   if(iused.eq.0) then
          vwptoq=dsqrt(dot(vwpto,vwpto,3))
          if(vwptoq.eq.0.0d0) then
            eq=dsqrt(dot(e,e,3))
            if(eq.ne.0.0d0) then ! from pers
              vwpto(1) = e(1)
              vwpto(2) = e(2)
              vwpto(3) = e(3)
            else ! no suggestion
              vwpto(1) = 3.0d0*vmax(1)
              vwpto(2) = 2.0d0*vmax(2)
              vwpto(3) = 1.5d0*vmax(3)
            end if
          end if
          if(ior.lt.0) write(*,2001) vmin,vmax,vwpto
          call dinput(vwpt,3)
          vwptq=dsqrt(dot(vwpt,vwpt,3))
          if(vwptq.eq.0.0d0) then
            vwpt(1)  = vwpto(1)
            vwpt(2)  = vwpto(2)
            vwpt(3)  = vwpto(3)
          else
            vwpto(1) = vwpt(1)
            vwpto(2) = vwpt(2)
            vwpto(3) = vwpt(3)
          end if
cww   else
cww       vwpt(1)  = vwpto(1)
cww       vwpt(2)  = vwpto(2)
cww       vwpt(3)  = vwpto(3)
cww   end if

      if(idev.eq.4) call clwclose(1,2,0)
      return
cww 2000 format(' Use old parameters? (0 = new or 1 = old) :',$)
 2001 format(' Body occupies the space with:',/,
     1  '                 X',13x,'Y',13x,'Z',/,
     2  '   Minimum ',1p,3e14.5,/,'   Maximum ',1p,3e14.5,/,
     3  ' Enter coordinates of view point   (X,Y,Z).',/,
     4  ' Default: X=',1pe9.2,', Y=',1pe9.2,', Z=',1pe9.2,/,'  >',$)
      end
c
c-----------------------------------------------------------------------
c
      subroutine plmaset(ipma,nummat,k1,k2,k3)
c-----------------------------------------------------------------------
c
c.... Purpose: set field for material numbers to plot
c
c     Input:
c     ipma(nummat)- array of materials
c     k1 < 0 : set all values to 0
c     k1 = 0 : set all values to i (plot for all numbers)
c     k1 > 0 : set values k1 to k2, inc = k3
c
c     Output:
c
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      integer ipma(nummat)
      if(k1.lt.0) then
        do i = 1,nummat
          ipma(i)=0
        end do
      else if(k1.eq.0) then
        do i = 1,nummat
          ipma(i)=i
        end do
      else if(k1.gt.0.and.k1.le.nummat) then
c....   set values for  k1 to k2, inc = k3
        k2 = min(max(0,k2),nummat)
        if(k2.ne.0) then
          if(k3.eq.0) k3=1
          do i = k1,k2,k3
             ipma(i)=i
          end do
c....   set value for  k1
        else if(k2.eq.0) then
          ipma(k1) = k1
        end if
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      function iplma(ma)
c-----------------------------------------------------------------------
c
c.... Purpose: test if material is to plot
c
c     Input:
c       ma
c
c     Output:
c       iplma(ma) = 0 do not plot
c       iplma(ma) = 1        plot
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE cdata
      USE pdata7
      integer*4 iplma,ma,i
      iplma = 0
      if(ma.eq.0) return
      do 10 i = 1,nummat
         if(ma.eq.aipma(i)) then
         iplma = 1
         goto 20
       end if
10    continue
20    return
      end
c
c-----------------------------------------------------------------------
c
      function imuse()
c-----------------------------------------------------------------------
c
c.... Purpose: test if material is to plot
c
c     Input:
c
c     Output:
c     imuse=0     all materials used
c     imuse=1 not all materials used
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE cdata
      imuse = 0
      do i = 1,nummat
        if(iplma(i).eq.0) then
          imuse = 1
          return
        end if
      enddo
      return
      end
c
c-------------------------------------------------------------
c
      function iplmano(ix,n,nen1)
c-------------------------------------------------------------
c
c.... Purpose: test if material is to plot
c
c     Input:
c
c     Output:
c     iplmano if n1     node of ipma isw = 1 -->        plot
c     iplmano if n1 not node of ipma isw = 0 --> do not plot
c
c     W. Wagner BS UKA 04/09
c-------------------------------------------------------------
      USE cdata
      USE pdata7
      integer ix(nen1,*)
c.... plot if all matn. used
      if(imuse().eq.0) then
         iplmano = 1
         return
      end if
c.... find element/node
      iplmano = 0
      do 200 j=1,numel            ! loop over elements
           if(iplma(ix(nen1,j)).eq.0) goto 200
         do k = 1,nen             ! loop over element-nodes
           if( ix(k,j) .eq. n) then
              iplmano = 1
              goto 300
           end if
         enddo
200   continue
300   continue
      return
      end

c
c-----------------------------------------------------------------------
      subroutine pltmate
c-----------------------------------------------------------------------
c
c.... Purpose: plot legend for PLOT,MATE
c
c     W. Wagner IBS KIT 04/15
c-----------------------------------------------------------------------
      USE cdata  
      USE pdata2
      implicit double precision(a-h,o-z)
      character yy*21
c.... headline 
      ypos  = 0.8
      write(yy,1000)
      call drawtxt(1,1.0d0,ypos,1,1,21,yy)

c.....draw color boxes and text
      nmat=nummat 
      do i = 1,nmat
c....   color 
        icol  = i+1
c....   position
        xcorb = 1.00
        xcort = 1.05
        ycorb = 0.77 -i*0.04
        ycort = 0.77 -(i-1)*0.04
c....   text
        write(yy,1001) i 
        call drawtxt(1,xcort,ycort,1,1,11,yy)
c...    box
        call pppcol(icol)
        call ppbox(xcorb,ycorb, .025d0, .040d0,ipgl)
c...    border
        call pppcol(32)
        call ppbox(xcorb,ycorb, .025d0, .040d0,3)
      end do
1000  format('Material distribution')
1001  format('Material ',i2)
      return
      end
c-----------------------------------------------------------------------
c
      subroutine pltlink1(x,ndm,numnp,iplk1,n2)
c-----------------------------------------------------------------------
c
c.... Purpose: plot location  of all linked nodes
c
c     Inputs:
c         x(ndm,*)      - Nodal coordinates of mesh
c         ndm           - Spatial dimension of mesh
c         numnp         - Number of nodes in mesh
c         iplk1(i)      - List of linked nodes to plot
c                         i<0 master node, |i| = node number master
c                         i>0 slave  node,  i  = node number master(!)
c
c         n2            - >0  plot node numbers
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata1
      USE pdata2
      implicit double precision (a-h,o-z)
      logical zoom
      dimension x(ndm,numnp),iplk1(numnp)
      dx1 = .002/scale
      x3 = 0.0d0
      do 10 n = 1,numnp
        nn = iplk1(n)
        if(nn.eq.0) go to 10        ! no link
        if(nn.lt.0) call pppcol(3)  ! master
        if(nn.gt.0) call pppcol(5)  ! slave
        if(zoom(x(1,n),ndm,1)) then
          x1 = x(1,n)
          x2 = x(2,n)
          if(ndm.ge.3) x3 = x(3,n)
          call plotl(x1-dx1 , x2+dx1 , x3, 3)
          call plotl(x1-dx1 , x2+dx1 , x3, ipgl)
          call plotl(x1-dx1 , x2-dx1 , x3, 2)
          call plotl(x1+dx1 , x2-dx1 , x3, 2)
          call plotl(x1+dx1 , x2+dx1 , x3, 2)
          call plotl(x1-dx1 , x2+dx1 , x3, 2)
          if(ipgl.eq.1) call clpan
          call plotl(x1-12*dx1 , x2+2*dx1 , x3, 3)
          if(n2.ne.0) call plabl(n)
        end if
10    continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltlink2(iplk2,x,angl,ndm,ndf,numnp,k1,rk3,tra,vr)
c-----------------------------------------------------------------------
c
c.... Purpose: Plot 1-3-D boundary constraints for linked nodes
c
c     Inputs:
c         iplk2(ndf,numnp)
c         iplks(j,i)     - List of linking conditions
c                          j=1 linked
c                          j=0 not linked
c
c         x(ndm,numnp)  - Nodal coordinates of mesh
c         ndm           - Spatial dimension of mesh
c         numnp         - Number of nodes in mesh
c         k1            - plot link for dof k1
c         k1=0          - plot link for all dofs
c         rk3           - scaling factor
c
c     Output:
c
c     W. Wagner BS KIT 12/10
c-----------------------------------------------------------------------
      USE boun
      USE mdat2
      USE pdata1
      USE pdata7
      implicit double precision (a-h,o-z)
      logical zoom
      dimension iplk2(ndf,*),x(ndm,*),angl(*),tra(3,3),vr(3),xz(3)
      integer idl(ndf)

      pi=datan(1.d0)*4.d0

c.... plot boundary restraints (lines = fixed)
      dx1 = 0.015d0/scale*rk3
      do 100 n = 1,numnp
c....   local boundary cond., rotate in plane ij1,ij2
        if(angl(n).ne.0.0d0) then
          ang = angl(n)*pi/180.d0
          cs  = cos(ang)
          sn  = sin(ang)
          ij1 = ia(1)
          ij2 = ia(2)
          if(ij1.eq.1 .and. ij2.eq.2) ij3 = 3
          if(ij1.eq.1 .and. ij2.eq.3) ij3 = 2
          if(ij1.eq.2 .and. ij2.eq.3) ij3 = 1
        else
          cs  = 1.d0
          sn  = 0.d0
          ij1 = 1
          ij2 = 2
          ij3 = 3
        end if
c....   rotation matrix
        call pzero(trb,9)
        trb(ij1,ij1) =  cs
        trb(ij2,ij2) =  cs
        trb(ij1,ij2) =  sn
        trb(ij2,ij1) = -sn
        trb(ij3,ij3) =  1.d0
        call pzero(vrb,3)
c
        xz(1) = x(1,n)
        xz(2) = x(2,n)
        xz(3) = 0.d0
        if(ndm.ge.3) xz(3) = x(3,n)
        call plxtrn(xz,tra,vr,ndm,1)

c....   find local link values

        if(k1.eq.0) then
          do i=1,ndf
            idl(i) = iplk2(i,n)
          end do
        else
          call pzeroi(idl,ndf)
          idl(k1) = iplk2(k1,n)
        end if

c....   plot local link values
        if(zoom(xz,ndm,1)) then
          x1 = x(1,n)
          x2 = x(2,n)
          x3 = 0.d0
          if(ndm.ge.3) x3 = x(3,n)
c.....    Plate(u_3,phi_1,phi_2)
          if(ipla.eq.1) then
            if (idl(1) .gt. 0)
     +        call pltbou1(3,tra,vr,ndf,0)
            if (idl(2) .gt. 0)
     +        call pltbou1(4,tra,vr,ndf,0)
            if (idl(3) .gt. 0)
     +        call pltbou1(5,tra,vr,ndf,0)
c.....    2-D Beam, Axishell(u_1,u_2,phi_3)
          else if(ipla.eq.2) then
            if (idl(1) .gt. 0)
     +        call pltbou1(1,tra,vr,ndf,0)
            if (idl(2) .gt. 0)
     +        call pltbou1(2,tra,vr,ndf,0)
            if (idl(3) .gt. 0)
     +        call pltbou1(7,tra,vr,ndf,0)
          else
c.....      other problems
c.....      1-D problem (u_1) + higher order elements
            if (idl(1) .gt. 0)
     +          call pltbou1(1,tra,vr,ndf,0)
            if (ndm.ge.2. and. ndf.ge.2 ) then
c.....      2-D truss, plane stress (u_1,u_2) + higher order elements
              if (idl(2) .gt. 0)
     +            call pltbou1(2,tra,vr,ndf,0)
            end if
c.....      3-D truss,3-D brick (u_1,u_2,u_3) + higher order elements
            if (ndm.ge.3. and. ndf.ge.3 ) then
              if (idl(3) .gt. 0)
     +            call pltbou1(3,tra,vr,ndf,0)
            end if
c.....      3-D beam, shell(5),shell(6) (u_1,u_2,u_3,phi_1,phi_2,phi_3)
            if (ndm.ge.3. and. ndf.ge.4 ) then
              if (idl(4) .gt. 0)
     +            call pltbou1(4,tra,vr,ndf,0)
            end if
c.....      3-D beam, shell(5), shell(6) (u_1,u_2,u_3,phi_1,phi_2,phi_3)
            if (ndm.ge.3. and. ndf.ge.5 ) then
              if (idl(5) .gt. 0)
     +            call pltbou1(5,tra,vr,ndf,0)
            end if
c.....      3-D beam, shell(6) (u_1,u_2,u_3,phi_1,phi_2,phi_3)
            if (ndm.ge.3. and. ndf.ge.6 ) then
              if (idl(6) .gt. 0)
     +            call pltbou1(6,tra,vr,ndf,0)
            end if
c.....      elements with more then 6 dofs
            if(ndf.gt.6) then
              iaddf = 0
              do 200 ii=7,ndf
                if(idl(ii).gt.0) iaddf = 1
200           continue
              if(iaddf.eq.1) call pltbou1(8,tra,vr,ndf,0)
            end if
          end if
        end if
100   continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltload(x,f,id,ndm,ndf,numnp,ityp)
c-----------------------------------------------------------------------
c
c.... Purpose: plot position of all nodes
c
c     Inputs:
c      x(ndm,numnp) -  coordinates
c      f(ndf,numnp) -  load vector
c      ityp = 1     - which are loaded
c      ityp = 3     - which have prescribed displacements .ne.0
c      ityp = 2     - which have values with respect to input
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata1
      USE pdata2
      implicit double precision (a-h,o-z)
      logical zoom
      dimension x(ndm,*),f(ndf,*),id(ndf,*)
      dx1 = .002/scale
      do n = 1,numnp
        x3 = 0.0
        fm = 0.d0
        if(ityp.eq.1) then ! only for free nodes
          do i = 1,ndf
            if(id(i,n).gt.0) fm = fm + f(i,n)**2
          end do
        else if(ityp.eq.3) then ! only for nodes with b.c.
          do i = 1,ndf
            if(id(i,n).le.0) fm = fm + f(i,n)**2
          end do
        else if(ityp.eq.2) then ! for all nodes
          do i = 1,ndf
            fm = fm + f(i,n)**2
          end do
        end if
        if(fm.gt.1.e-10) then
            if(zoom(x(1,n),ndm,1)) then
              x1 = x(1,n)
              x2 = x(2,n)
              if(ndm.ge.3) x3 = x(3,n)
              call plotl(x1-dx1 , x2+dx1 , x3, 3)
              call plotl(x1-dx1 , x2+dx1 , x3, ipgl)
              call plotl(x1-dx1 , x2-dx1 , x3, 2)
              call plotl(x1+dx1 , x2-dx1 , x3, 2)
              call plotl(x1+dx1 , x2+dx1 , x3, 2)
              call plotl(x1-dx1 , x2+dx1 , x3, 2)
              if(ipgl.eq.1) call clpan
c....       node number
cww           call plotl(x1-12*dx1 , x2+2*dx1 , x3, 3)
cww           call plabl(n)
            end if
          end if
      enddo
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltmain(x,st,ndm,numnp,tra,vr,k2,rk3,rk4,ipgl)
c-----------------------------------------------------------------------
c
c.... Purpose: plot main-stresses
c
c     Inputs:
c       x(ndm,numnp) -  coordinates
c       st(numnp,*)  -  stresses
c       k2=2: [ S_x=st(*,1), S_xy=st(*,2), S_y = st(*,3) ]
c       stresses are in 1-2 plane
c       sig(1) = S_x    sig(4) = S_1
c       sig(2) = S_xy   sig(5) = S_2
c       sig(3) = S_y    sig(6) = alpha
c
c       k2=3: [ S_xx=st(*,1), S_yy=st(*,2), S_zz = st(*,3) ]
c             [ S_xy=st(*,4), S_xz=st(*,2), S_yz = st(*,3) ]
c
c       rk3          - scaling factor for lenght of arrow
c       rk4          - scaling factor for lenght of tip
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata1
      USE pdata4
      USE pdatas
      implicit double precision (a-h,o-z)
      dimension dd(2),xx(3,9),x(ndm,*),st(numnp,*),sig(6),tra(3,3),
     1          vr(3),txm(3),smin(3),smax(3),sig3(6),sig0(3),z(3,3)
      pi=datan(1.d0)*4.d0
c...  typ 2=2D,3=3D
      if(k2.eq.0) k2 = 2
      k2 = max(2,min(k2,3))
c.... compute longest vector
      fm = 0.d0
      call pzero(smin,3)
      call pzero(smax,3)
c
      do 90 n = 1,numnp
        txm(1) = x(1,n)*fact
        txm(2) = x(2,n)*fact
        txm(3) = 0.0
c....   transform coordinates
        call plxtrn(txm,tra,vr,ndm,1)
cww ??       if(txm(1).lt.xmin(1).or.txm(1).gt.xmax(1)) go to 90
cww ??       if(txm(2).lt.xmin(2).or.txm(2).gt.xmax(2)) go to 90
cww        ...3 = ???
        if(k2.eq.2) then !2D
          sig(1) = st(n,1)
          sig(2) = st(n,2)
          sig(3) = st(n,3)
          call pstres(sig,sig(4),sig(5),sig(6))
          do 80 ii = 4,5
            d = dabs(sig(ii))
80        fm = dmax1(fm,d)
          smin(1) = dmin1(smin(1),sig(4))
          smax(1) = dmax1(smax(1),sig(4))
          smin(2) = dmin1(smin(2),sig(5))
          smax(2) = dmax1(smax(2),sig(5))
        else !3D
          sig3(1) = st(n,1)
          sig3(2) = st(n,2)
          sig3(3) = st(n,3)
          sig3(4) = st(n,4)
          sig3(5) = st(n,5)
          sig3(6) = st(n,6)
          call pstres1(sig3,sig0,z)
          do ii = 1,3
            d = dabs(sig0(ii))
            fm = dmax1(fm,d)
          enddo
          smin(1) = dmin1(smin(1),sig0(1))
          smax(1) = dmax1(smax(1),sig0(1))
          smin(2) = dmin1(smin(2),sig0(2))
          smax(2) = dmax1(smax(2),sig0(2))
          smin(3) = dmin1(smin(3),sig0(3))
          smax(3) = dmax1(smax(3),sig0(3))
        end if
90    continue
      fm = fm*scale*20.d0/rk3
c...  plot min/max on screen
      call pltmain1(smin,smax,k2)

c.... compute vector at each node
      do 120 n = 1,numnp
        x1 = x(1,n)*fact
        x2 = x(2,n)*fact
        txm(1) = x1
        txm(2) = x2
        txm(3) = 0.0
        if(ndm.gt.2) txm(3) = x(3,n)
c....   transform coordinates
        call plxtrn(txm,tra,vr,ndm,1)
cww ??       if(txm(1).lt.xmin(1).or.txm(1).gt.xmax(1)) go to 120
cww ??       if(txm(2).lt.xmin(2).or.txm(2).gt.xmax(2)) go to 120
c            3.. ????
        if(k2.eq.2) then !2D
          call pzero(xx,27)
          if(ndm.gt.2) then
            do 130 i = 1,9
              xx(3,i) = x(3,n)
130         continue
          end if
          sig(1) = st(n,1)
          sig(2) = st(n,2)
          sig(3) = st(n,3)
          call pstres(sig,sig(4),sig(5),sig(6))
          alfa = pi*sig(6)/180.d0 ! use only alpha
c
          do 118 ii = 4,5
            if(ii.eq.5) alfa = + alfa + pi/2.d0
            dd(1) = dabs(sig(ii))*dcos(alfa)/fm
            dd(2) = dabs(sig(ii))*dsin(alfa)/fm
c.....      modify due to quadrant for symmetry
            if(isym2(iadd+1,1,its).eq.1) dd(1) = - dd(1)
            if(isym2(iadd+1,2,its).eq.1) dd(2) = - dd(2)
c....       compute plot locations
            if(sig(ii).gt.0.d0) then
              xx(1,5) = x(1,n)  - dd(1)/2.d0
              xx(2,5) = x(2,n)  - dd(2)/2.d0
              xx(1,1) = xx(1,5) + dd(1)
              xx(2,1) = xx(2,5) + dd(2)
              xx(1,2) = xx(1,1) -.3d0*dd(1)*rk4 - .1d0*dd(2)*rk4
              xx(2,2) = xx(2,1) -.3d0*dd(2)*rk4 + .1d0*dd(1)*rk4
              xx(1,3) = xx(1,1) -.3d0*dd(1)*rk4 + .1d0*dd(2)*rk4
              xx(2,3) = xx(2,1) -.3d0*dd(2)*rk4 - .1d0*dd(1)*rk4
              xx(1,4) = xx(1,1)
              xx(2,4) = xx(2,1)
              xx(1,6) = xx(1,5) +.3d0*dd(1)*rk4 + .1d0*dd(2)*rk4
              xx(2,6) = xx(2,5) +.3d0*dd(2)*rk4 - .1d0*dd(1)*rk4
              xx(1,7) = xx(1,5) +.3d0*dd(1)*rk4 - .1d0*dd(2)*rk4
              xx(2,7) = xx(2,5) +.3d0*dd(2)*rk4 + .1d0*dd(1)*rk4
c....         transform the vector
              call plxtrn(xx,tra,vr,ndm,9)
c....         plot the vector
              call pppcol(3)
              call plotl(xx(1,1),xx(2,1),xx(3,1),3)
              call plotl(xx(1,1),xx(2,1),xx(3,1),ipgl)
              do 110 i = 2,8
                j = i
                if(i.eq.8) j = 5
                call plotl(xx(1,j),xx(2,j),xx(3,j),2)
110           continue
              if(ipgl.eq.1) call clpan
c
            else
              xx(1,6) = x(1,n)  - dd(1)/2.d0
              xx(2,6) = x(2,n)  - dd(2)/2.d0
              xx(1,1) = xx(1,6) + dd(1)
              xx(2,1) = xx(2,6) + dd(2)
              xx(1,2) = xx(1,1) - .1d0*dd(2)*rk4
              xx(2,2) = xx(2,1) + .1d0*dd(1)*rk4
              xx(1,3) = xx(1,1) - .3d0*dd(1)*rk4
              xx(2,3) = xx(2,1) - .3d0*dd(2)*rk4
              xx(1,4) = xx(1,1) + .1d0*dd(2)*rk4
              xx(2,4) = xx(2,1) - .1d0*dd(1)*rk4
              xx(1,5) = xx(1,1)
              xx(2,5) = xx(2,1)
              xx(1,7) = xx(1,6) + .1d0*dd(2)*rk4
              xx(2,7) = xx(2,6) - .1d0*dd(1)*rk4
              xx(1,8) = xx(1,6) + .3d0*dd(1)*rk4
              xx(2,8) = xx(2,6) + .3d0*dd(2)*rk4
              xx(1,9) = xx(1,6) - .1d0*dd(2)*rk4
              xx(2,9) = xx(2,6) + .1d0*dd(1)*rk4
c....         transform the vector
              call plxtrn(xx,tra,vr,ndm,9)
c....         plot the vector
              call pppcol(2)
              call plotl(xx(1,1),xx(2,1),xx(3,1),3)
              call plotl(xx(1,1),xx(2,1),xx(3,1),ipgl)
              do 112 i = 2,10
                j = i
                if(i.eq.10) j = 6
                call plotl(xx(1,j),xx(2,j),xx(3,j),2)
112           continue
              if(ipgl.eq.1) call clpan
            end if
118       continue
        else ! 3D-case
          sig3(1) = st(n,1)
          sig3(2) = st(n,2)
          sig3(3) = st(n,3)
          sig3(4) = st(n,4)
          sig3(5) = st(n,5)
          sig3(6) = st(n,6)
          call pstres1(sig3,sig0,z)
          do ii = 1,3
            call pzero(xx,27)
            dds = dabs(sig0(ii))/fm/2.0d0
            xx(1,1) = x(1,n)  - dds*z(1,ii) ! begin
            xx(2,1) = x(2,n)  - dds*z(2,ii)
            xx(3,1) = x(3,n)  - dds*z(3,ii)
            xx(1,2) = x(1,n)  + dds*z(1,ii) ! end
            xx(2,2) = x(2,n)  + dds*z(2,ii)
            xx(3,2) = x(3,n)  + dds*z(3,ii)
c....       transform the vector
            call plxtrn(xx,tra,vr,ndm,9)
c....       plot the vector
            call pppcol(2) ! blue
            if(sig0(ii).ge.0.0d0) call pppcol(3) ! red
            call plotl(xx(1,1),xx(2,1),xx(3,1),3)
            call plotl(xx(1,2),xx(2,2),xx(3,2),2)
          enddo
        end if
120   continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltmain1(smin,smax,k2)
c-----------------------------------------------------------------------
c
c.... Purpose: plot min/max for main-stresses
c
c     Inputs:
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      character*15 strph
      character*15 strs(3)
      dimension smin(3),smax(3)
      data strs/' STRESS   S_1  ',' STRESS   S_2  ',' STRESS   S_3  '/
      x = 1.d0
      y = 0.80d0
      write(strph,'(a15)' ) strs(1)
      call drawtxt(1,x,y,1,1,15,strph)
      y = 0.75d0
      write(strph, '(1p1e10.3,4h min)' ) smin(1)
      call drawtxt(1,x,y,1,1,15,strph)
      y = 0.70d0
      write(strph, '(1p1e10.3,4h max)' ) smax(1)
      call drawtxt(1,x,y,1,1,15,strph)
c
      y = 0.60d0
      write(strph,'(a15)' ) strs(2)
      call drawtxt(1,x,y,1,1,15,strph)
      y = 0.55d0
      write(strph, '(1p1e10.3,4h min)' ) smin(2)
      call drawtxt(1,x,y,1,1,15,strph)
      y = 0.50d0
      write(strph, '(1p1e10.3,4h max)' ) smax(2)
      call drawtxt(1,x,y,1,1,15,strph)
c
      if(k2.eq.3) then
        y = 0.40d0
        write(strph,'(a15)' ) strs(3)
        call drawtxt(1,x,y,1,1,15,strph)
        y = 0.35d0
        write(strph, '(1p1e10.3,4h min)' ) smin(3)
        call drawtxt(1,x,y,1,1,15,strph)
        y = 0.30d0
        write(strph, '(1p1e10.3,4h max)' ) smax(3)
        call drawtxt(1,x,y,1,1,15,strph)
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltmaxi(x,ndm,numnp)
c-----------------------------------------------------------------------
c
c.... Purpose: plot + print location of nodes with maximal values
c              of e.g. stress
c
c     Inputs:
c        x
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata1
      USE pdata2
      implicit double precision (a-h,o-z)
      logical zoom
      dimension x(ndm,*),nmax(2)
      common /rndata/ nmn,nmx
      nmax(1) = nmn
      nmax(2) = nmx
      write(*,1000) nmn
      write(*,1001) nmx
      dx1 = .002/scale
      x3 = 0.0
      do 100 nn = 1,2
        n = nmax(nn)
        if(zoom(x(1,n),ndm,1)) then
          x1 = x(1,n)
          x2 = x(2,n)
          if(ndm.ge.3) x3 = x(3,n)
          call plotl(x1-dx1 , x2+dx1 , x3, 3)
          call plotl(x1-dx1 , x2+dx1 , x3, ipgl)
          call plotl(x1-dx1 , x2-dx1 , x3, 2)
          call plotl(x1+dx1 , x2-dx1 , x3, 2)
          call plotl(x1+dx1 , x2+dx1 , x3, 2)
          call plotl(x1-dx1 , x2+dx1 , x3, 2)
          if(ipgl.eq.1) call clpan
          call plotl(x1-12*dx1 , x2+2*dx1 , x3, 3)
          call plabl(n)
        end if
100   continue
1000  format(' Minimum Value at node ',i8)
1001  format(' Maximum Value at node ',i8)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltnod(x,ix,nen1,ndm,numnp,n1,n2)
c-----------------------------------------------------------------------
c
c.... Purpose: plot position of all nodes
c
c     Inputs:
c        x
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata1
      USE pdata2
      implicit double precision (a-h,o-z)
      logical zoom
      dimension x(ndm,*)
      integer ix(nen1,*)
      if(scale.eq.0.d0) then
        call drawmess(' Scale = 0 in PLTNOD',1,1)
        return
      end if
      dx1 = .002/scale
      x3 = 0.0
      call pppcol(n2)
      do 100 n = 1,numnp
        if(n1.lt.0.and.abs(n1).ne.n) goto 100  ! only one node
      if(iplmano(ix,n,nen1).eq.0)  goto 100  ! only if matn
        if(zoom(x(1,n),ndm,1)) then
          x1 = x(1,n)
          x2 = x(2,n)
          if(ndm.ge.3) x3 = x(3,n)
          call plotl(x1-dx1 , x2+dx1 , x3, 3)
          call plotl(x1-dx1 , x2+dx1 , x3, ipgl)
          call plotl(x1-dx1 , x2-dx1 , x3, 2)
          call plotl(x1+dx1 , x2-dx1 , x3, 2)
          call plotl(x1+dx1 , x2+dx1 , x3, 2)
          call plotl(x1-dx1 , x2+dx1 , x3, 2)
          if(ipgl.eq.1) call clpan
          call plotl(x1-12*dx1 , x2+2*dx1 , x3, 3)
          if(n1.ne.0) then
            call pppcol(1) ! plot number in black
            call plabl(n)
            call pppcol(n2)
          end if
        end if
100   continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltnod1(x,ndm,numnp,n1)
c-----------------------------------------------------------------------
c
c.... Purpose: plot position of all nodes with neg D_ii
c
c     Inputs:
c        x
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE dii
      USE pdata1
      USE pdata2
      implicit double precision (a-h,o-z)
      logical zoom
      dimension x(ndm,*)
      dx1 = .002/scale
      x3 = 0.0
      if(ii(n1).eq.0) then
        call drawmess('No points available',6,1)
        return
      end if
      do 100 nn = 1,ii(n1)
        n = ndii(n1,nn)
        if(zoom(x(1,n),ndm,1)) then
          x1 = x(1,n)
          x2 = x(2,n)
          if(ndm.ge.3) x3 = x(3,n)
          call plotl(x1-dx1 , x2+dx1 , x3, 3)
          call plotl(x1-dx1 , x2+dx1 , x3, ipgl)
          call plotl(x1-dx1 , x2-dx1 , x3, 2)
          call plotl(x1+dx1 , x2-dx1 , x3, 2)
          call plotl(x1+dx1 , x2+dx1 , x3, 2)
          call plotl(x1-dx1 , x2+dx1 , x3, 2)
          if(ipgl.eq.1) call clpan
          call plotl(x1-12*dx1 , x2+2*dx1 , x3, 3)
          call plabl(n)
        end if
100   continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltnod2(x,ndm,numnp,n1)
c-----------------------------------------------------------------------
c
c.... Purpose: plot position of all nodes on intersections
c
c     Inputs:
c        x
c
c     Output:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE isecn
      USE pdata1
      USE pdata2
      implicit double precision (a-h,o-z)
      logical zoom
      dimension x(ndm,*)
      dx1 = .002/scale
      x3 = 0.0
      if(is.eq.0) then
        call drawmess('No points available',6,1)
        return
      end if
      do 100 nn = 1,is
        n = isecno(nn)
        if(zoom(x(1,n),ndm,1)) then
          x1 = x(1,n)
          x2 = x(2,n)
          if(ndm.ge.3) x3 = x(3,n)
          call plotl(x1-dx1 , x2+dx1 , x3, 3)
          call plotl(x1-dx1 , x2+dx1 , x3, ipgl)
          call plotl(x1-dx1 , x2-dx1 , x3, 2)
          call plotl(x1+dx1 , x2-dx1 , x3, 2)
          call plotl(x1+dx1 , x2+dx1 , x3, 2)
          call plotl(x1-dx1 , x2+dx1 , x3, 2)
          if(ipgl.eq.1) call clpan
          call plotl(x1-12*dx1 , x2+2*dx1 , x3, 3)
          if(n1.gt.0) call plabl(n)
        end if
100   continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltord(iel,iju,jplt)
c-----------------------------------------------------------------------
c
c.... Purpose: set up the plot order for special elements
c
c     Inputs:
c       iel  -  user element number
c
c     Output:
c       iju  -  number of nodes to plot element
c       jplt -  local node numbers to plot element
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata6
      USE pdatah
      dimension iplt(9),jplt(*)
c.... default orders for the 3-9 node 2-d elements
      data iplt/1,5,2,6,3,7,4,8,1/
c.... get the number of plot points to go around the element
      if(hide) then
          iju = 5
          do 50 ij = 1,iju-1
            jplt(ij) = ij
50        continue
          jplt(iju) = 1
      else
          iju = inord(iel)
          if(iju.eq.0) iju = 9
c....   set plot table
          do 100 ij = 1,iju
            if(inord(iel).eq.0) then
              jplt(ij) = iplt(ij)
            else
              jplt(ij) = ipord(ij,iel)
            end if
100       continue
      end if
      end
c
c----------------------------------------------------------------------+
c
      subroutine pltprf(jp,neq,lower)
c----------------------------------------------------------------------+
c
c     Purpose: Display plot of profile layout
c
c     Inputs:
c        jp(*)=ia   pointers or row pointers
c        ja(*)    - Column pointers
c        neq      - Number of equations
c        lower    - true=plot lower part
c
c     Outputs:
c        plot distribution of matrix A
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------+
      USE iscsr
      USE pdata1
      USE soltyp
      logical         lower,lower1
      integer         n,neq,jp(*)
      real*8          x0,y0, x,y,c, pfact,rfact,ptone,ptnin,onept,ds

      rfact = 0.5d0/(scaleg*fact)
      pfact = 0.8d0/dble(neq)*rfact
      x0    = 0.5d0*sx(1) - 0.5d0*rfact
      y0    = 0.5d0*sx(2) - 0.5d0*rfact
      ptone = 0.1d0*rfact
      onept = 1.0d0*rfact
      ptnin = 0.9d0*rfact
      ds = 0.4d0*rfact/neq
      if ( istyp.eq.0) then ! standard solver(0)

        if(neq.gt.100) then  ! lines , original
c         Upper part
          call pppcol(3)
          do n = 2,neq
            x = ptone + dble(n)*pfact
            y = onept - x + y0
            c = dble(jp(n) - jp(n-1))*pfact
            x = x + x0
            call plotl( x    , y    , 0.0d0, 3)
            call plotl( x    , y + c, 0.0d0, 2)
          end do
c         Lower part
          if(lower) then
            do n = 2,neq
              x = ptone + dble(n)*pfact
              y = onept - x + y0
              c = dble(jp(n) - jp(n-1))*pfact
              x = x + x0
              call plotl( x - c, y    , 0.0d0, 3)
              call plotl( x    , y    , 0.0d0, 2)
            end do
          end if

          else   ! with boxes for low resolution
          call pppcol(3)
c         1,1
          pix =  1 -0.5d0
          piy =  0 -0.5d0
          x = ptone  + pix*pfact + x0
          y = ptnin  + piy*pfact + y0
          call ppbox2(x,y,ds)

c         Upper part
          do n = 2,neq
            do k = -n+1,jp(n)-jp(n-1)-n+1
              pix =  n -0.5d0
              piy =  k -0.5d0
              x = ptone  + pix*pfact + x0
              y = ptnin  + piy*pfact + y0
              call ppbox2(x,y,ds)
            end do
          end do
c         Lower part
          if(lower) then
          do n = 2,neq
            do k = -n+1,jp(n)-jp(n-1)-n+1
              pix =  k -0.5d0
              piy =  n -0.5d0
              x = ptone  - pix*pfact + x0
              y = ptnin  - piy*pfact + y0
              call ppbox2(x,y,ds)
            end do
          end do
          end if
        end if

      else if (istyp.eq.1.or.istyp.eq.2) then
        call prof_csr(pfact,ptone,ptnin,x0,y0,ds,neq,jp,csrja,lower)


      else if (istyp.ge.3.and.istyp.le.8) then ! all other CSR solver
        if(isymcsr.eq.1) lower1=.false.
        if(isymcsr.eq.2) lower1=lower
        call prof_csr(pfact,ptone,ptnin,x0,y0,ds,neq,jp,csrja,lower1)

      end if

c.... for all solvers: Diagonal line and boundary lines
      call pppcol(32)
      call plotl( ptone + x0, ptnin + y0, 0.0d0, 3)
      call plotl( ptnin + x0, ptone + y0, 0.0d0, 2)
      call plotl( ptnin + x0, ptnin + y0, 0.0d0, 2)
      call plotl( ptone + x0, ptnin + y0, 0.0d0, 2)
      call plotl( ptone + x0, ptone + y0, 0.0d0, 2)
      call plotl( ptnin + x0, ptone + y0, 0.0d0, 2)

      end
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c
      subroutine pltrea(x,st,ndf,ndm,n1,n2,scfac,tra,vr,angl,ix,nen1)
c-----------------------------------------------------------------------
c
c     Purpose: plot reactions
c              storage  [ st(*,x) , st(*,y) , st(*,z)]
c
c     Inputs:
c
c
c     Outputs:
c
c     Remark: before using tang pltrea shows all nodal loads           |
c     set all values negative, only if tang is not used positive       |
c.... assumes that dofs occur only in the way that idm=idf             |
c.... thus only forces can be plotted ww 23.10.03                      |
c.... old version, routine will be not used, see pplotf>[reac]         |
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE fdata
      USE mdat2
      USE pdata1
      USE pdata2
      USE pdata7
      USE pdatas
      implicit double precision(a-h,o-z)
      dimension x(ndm,*),st(ndf,*),dd(3),xx(3,4),tra(3,3),vr(3),
     1          angl(*),xz(3),ix(nen1,*),xmax(2,3)
      logical zoom
      pi=datan(1.d0)*4.d0
c.... test if tang has been used
      fac = -1.d0
      if(fl(4)) fac= 1.d0
c.... multiply all terms of R by fac
c.... compute longest vector
      do i = 1,ndm
        xmax(1,i)= 1.e+10
        xmax(2,i)=-1.e+10
      enddo
      fm = 0.d0
      do 100 n = n1,n2
        if(iplmano(ix,n,nen1).eq.0)  goto 100  ! only if matn
        xz(1) = x(1,n)
        xz(2) = x(2,n)
        xz(3) = 0.d0
        if(ndm.ge.3) xz(3) = x(3,n)
        call plxtrn(xz,tra,vr,ndm,1)
        if(zoom(xz,ndm,1)) then
          if(ipla.eq.1) then
            d = st(1,n)**2
            xmax(1,1)=min(xmax(1,1),st(1,n)*fac)
            xmax(2,1)=max(xmax(2,1),st(1,n)*fac)
          else
            d = 0.d0
            do ii = 1,ndm
              xmax(1,ii)= min(xmax(1,ii),st(ii,n)*fac) ! extremal values
              xmax(2,ii)= max(xmax(2,ii),st(ii,n)*fac)
              d = d + st(ii,n)**2
            enddo
          end if
          fm = max(fm,sqrt(d))
        end if
100   continue
c.... for extremal values
      if(scfac.le.0.0d0) scfac = 1.0d0
      if(fm.le.0.0d0) then
        call drawmess(' No non-zero values acting on mesh',1,1)
      else
        fm = fm*scale*40.d0/scfac
c....   compute vector at each node
        do 120 n = n1,n2
          if(iplmano(ix,n,nen1).eq.0)  goto 120  ! only if matn
          xz(1) = x(1,n)
          xz(2) = x(2,n)
          xz(3) = 0.d0
          if (ndm.eq.3) xz(3) = x(3,n)
          call plxtrn(xz,tra,vr,ndm,1)
          if(zoom(xz,ndm,1)) then
            x3 = 0.0
            call pzero(xx,12)
            call pzero(dd,3)
            x1 = x(1,n)
            x2 = x(2,n)
            if(ndm.eq.3) x3 = x(3,n)
            if(ipla.eq.1) then
              dd(3) = fac*st(1,n)/fm
            else
                           dd(1) = fac*st(1,n)/fm
              if(ndf.gt.1) dd(2) = fac*st(2,n)/fm
              if(ndm.eq.3) dd(3) = fac*st(3,n)/fm
            end if
c....       transform in 1-2 plane
            if(angl(n).ne.0.0d0) then
              ang = angl(n)*pi/180.d0
              cs = cos(ang)
              sn = sin(ang)
              ij1   = ia(1)
              ij2   = ia(2)
              tm      =  dd(ij1)*cs - dd(ij2)*sn
              dd(ij2) =  dd(ij1)*sn + dd(ij2)*cs
              dd(ij1) =  tm
            end if
c.....      modify due to quadrant for symmetry
            if(isym2(iadd+1,1,its).eq.1) dd(1) = - dd(1)
            if(isym2(iadd+1,2,its).eq.1) dd(2) = - dd(2)
            if(isym2(iadd+1,3,its).eq.1) dd(3) = - dd(3)
c....       compute plot locations
            xx(1,1)=x(1,n)
            xx(2,1)=x(2,n)
            if(ndm.eq.3) xx(3,1)=x(3,n)
            xx(1,2)=xx(1,1) + dd(1)
            xx(2,2)=xx(2,1) + dd(2)
            if(ndm.eq.3)xx(3,2)=xx(3,1) + dd(3)
            xx(1,3)=xx(1,2) -.6*dd(1) - .2*(dd(2)+dd(3))
            xx(2,3)=xx(2,2) -.6*dd(2) + .2*(dd(1)+dd(3))
            if(ndm.eq.3)xx(3,3)=xx(3,2) -.6*dd(3) + .2*(dd(1)+dd(2))
            xx(1,4)=xx(1,2) -.6*dd(1) + .2*(dd(2)+dd(3))
            xx(2,4)=xx(2,2) -.6*dd(2) - .2*(dd(1)+dd(3))
            if(ndm.eq.3) xx(3,4)=xx(3,2) -.6*dd(3) - .2*(dd(1)+dd(2))
c....       transform the vector
            call plxtrn(xx,tra,vr,ndm,4)
c....       plot the vector
            call plotl(xx(1,1),xx(2,1),xx(3,1),3)
            call plotl(xx(1,2),xx(2,2),xx(3,2),2)
            call plotl(xx(1,2),xx(2,2),xx(3,2),ipgl)
            call plotl(xx(1,3),xx(2,3),xx(3,3),2)
            call plotl(xx(1,4),xx(2,4),xx(3,4),2)
            call plotl(xx(1,2),xx(2,2),xx(3,2),2)
            if(ipgl.eq.1) call clpan
          end if
120     continue
      end if
c...  plot extremal values
      nval = min(ndf,3)
      call plfleg(idev,'      Nodal Reaction',1,nval,xmax)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltrns(n,ct,tr,vr)
c-----------------------------------------------------------------------
c
c     Purpose: set rotations for isometric plots
c
c     Inputs:
c
c
c     Outputs:
c       tr(3,3) - rotation matrix
c       vr(3)   - coordinate shift
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension ct(3),tr(3,3),vr(3)
      pi=datan(1.d0)*4.d0
      ang = ct(1)*pi/180.d0
      cn  = cos(ang)
      sn  = sin(ang)
      j   = mod(n,3) + 1
      k   = mod(j,3) + 1
c.... construct the rotation matrix and coordinate shifts
      do 100 i = 1,3
        te      =  cn*tr(j,i) + sn*tr(k,i)
        tr(k,i) = -sn*tr(j,i) + cn*tr(k,i)
        tr(j,i) =  te
        te      =  cn*(vr(j) - ct(2)) - sn*(vr(k) - ct(3))
        vr(k)   =  sn*(vr(j) - ct(2)) + cn*(vr(k) - ct(3))
        vr(j)   =  te
100   continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltrot(tr,ndm,rotang,ityp)
c-----------------------------------------------------------------------
c
c     Purpose: plot  axes and values for rot-macro
c
c     Inputs:
c       ityp=1 at center
c       ityp=2 at left upper box
c
c     Outputs:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension tr(3,3),xx(3,5),xi(3),rotang(3)
      character*8 yy
      if(ndm.lt.3) return
      call pppcol(5)
c.... center of axis at midpoint
      if(ityp.eq.1) then
        ds = 0.25d0
        xi(1)=0.5d0
        xi(2)=0.5d0
        xi(3)=0.0d0
      else if(ityp.eq.2) then
        ds = 0.07d0
        xi(1)=1.047325d0
        xi(2)=0.897500d0
        xi(3)=0.0d0
      end if
c.... perspective projecton of axes
      do 120 m = 1,ndm
      call pzero(xx,15)
      do 100 n = 1,ndm
      fac1 = tr(m,n)*ds
      xx(n,1) = xi(n)
      xx(n,2) = xx(n,1) + fac1
100   xx(n,5) = xx(n,2)
      fac1 = tr(m,1)*ds
      fac2 = tr(m,2)*ds
      fac3 = tr(m,3)*ds
      xx(1,3) = xx(1,2) -.3*fac1 - .1*(fac2+fac3)
      xx(2,3) = xx(2,2) -.3*fac2 + .1*(fac1+fac3)
      xx(3,3) = xx(3,2) -.3*fac3 + .1*(fac1+fac2)
      xx(1,4) = xx(1,2) -.3*fac1 + .1*(fac2+fac3)
      xx(2,4) = xx(2,2) -.3*fac2 - .1*(fac1+fac3)
      xx(3,4) = xx(3,2) -.3*fac3 - .1*(fac1+fac2)
c.... plot the vector
      call dplot(xx(1,1),xx(2,1),3,0)
      call dplot(xx(1,2),xx(2,2),2,0)
      call dplot(xx(1,2),xx(2,2),1,0)
      call dplot(xx(1,3),xx(2,3),2,0)
      call dplot(xx(1,4),xx(2,4),2,0)
      call dplot(xx(1,2),xx(2,2),2,0)
      call clpan
c.... plot the label
      dx = 0.d0
      if(ityp.eq.2) dx = ds
      call dplot(xx(1,2)-dx,xx(2,2),3,0)
      call plabl(m)
120   continue
c.... plot the values
      x = 1.13
      y = 0.95
      write(yy,'(3ha1=,f5.0)') rotang(1)
      call drawtxt(1,x,y,1,1,8,yy)
      y = 0.91
      write(yy,'(3ha2=,f5.0)') rotang(2)
      call drawtxt(1,x,y,1,1,8,yy)
      y = 0.87
      write(yy,'(3ha3=,f5.0)') rotang(3)
      call drawtxt(1,x,y,1,1,8,yy)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltrsum(x,nodesrf,ndm,numnp,n1)
c-----------------------------------------------------------------------
c
c     Purpose: plot position of nodes for rsum
c
c     Inputs:
c
c     Outputs:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata1
      USE pdata2
      USE rsum
      implicit double precision (a-h,o-z)
      logical zoom
      dimension x(ndm,*),nodesrf(*)
      dx1 = .002/scale
      x3 = 0.0
      do 100 nn = 1,nfs1
        n = nodesrf(nn)
        if(zoom(x(1,n),ndm,1)) then
          x1 = x(1,n)
          x2 = x(2,n)
          if(ndm.ge.3) x3 = x(3,n)
          call plotl(x1-dx1 , x2+dx1 , x3, 3)
          call plotl(x1-dx1 , x2+dx1 , x3, ipgl)
          call plotl(x1-dx1 , x2-dx1 , x3, 2)
          call plotl(x1+dx1 , x2-dx1 , x3, 2)
          call plotl(x1+dx1 , x2+dx1 , x3, 2)
          call plotl(x1-dx1 , x2+dx1 , x3, 2)
          if(ipgl.eq.1) call clpan
          call plotl(x1-12*dx1 , x2+2*dx1 , x3, 3)
          if(n1.ne.0) call plabl(n)
        end if
100   continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltstr(dt,st,numnp)
c-----------------------------------------------------------------------
c
c     Purpose: calculate average stresses
c
c     Inputs:
c       st(numnp,*)  - stresses*weights
c       dt(numnp)    - weight at nodes
c
c     Outputs:
c       st(numnp,*)  - average stresses
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE strnam
      implicit double precision (a-h,o-z)
      dimension dt(numnp),st(numnp,*),sig(6)
      pi  = datan(1.0d0)*4.0d0
      do 100 ii = 1,numnp
        dh = dt(ii)
        if(dh.ne.0.0d0) then
c....     stresses 1,4
          do is=1,4
            st(ii,is) = st(ii,is)/dh
          end do

c....     stresses 5,7: normal or  principal stresses 2D
          if(istv.gt.0) then
            sig(1) = st(ii,1)
            sig(2) = st(ii,2)
            sig(3) = st(ii,3)
            call pstres(sig,sig(4),sig(5),sig(6))
            st(ii,5) = sig(4)
            st(ii,6) = sig(5)
            st(ii,7) = sig(6)
          else
            st(ii,5) = st(ii,5)/dh
            st(ii,6) = st(ii,6)/dh
            st(ii,7) = st(ii,7)/dh
          end if

c....     stresses 8,25
          do is=8,25
            st(ii,is) = st(ii,is)/dh
          end do

        end if
100   continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine plttie(x,ip,ndm,numnp,n1)
c-----------------------------------------------------------------------
c
c     Purpose: plot position of all tied nodes
c
c     Inputs:
c
c     Outputs:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata1
      USE pdata2
      implicit double precision (a-h,o-z)
      logical zoom
      dimension x(ndm,*),ip(numnp)
      dx1 = .002/scale
      x3 = 0.0
      do 100 n = 1,numnp
        if(n1.lt.0.and.abs(n1).ne.n) goto 100
        nn = ip(n)
        if(nn.eq.n) goto 100
        if(zoom(x(1,nn),ndm,1)) then
          x1 = x(1,nn)
          x2 = x(2,nn)
          if(ndm.ge.3) x3 = x(3,nn)
          call plotl(x1-dx1 , x2+dx1 , x3, 3)
          call plotl(x1-dx1 , x2+dx1 , x3, ipgl)
          call plotl(x1-dx1 , x2-dx1 , x3, 2)
          call plotl(x1+dx1 , x2-dx1 , x3, 2)
          call plotl(x1+dx1 , x2+dx1 , x3, 2)
          call plotl(x1-dx1 , x2+dx1 , x3, 2)
          if(ipgl.eq.1) call clpan
          call plotl(x1-12*dx1 , x2+2*dx1 , x3, 3)
          if(n1.ne.0) call plabl(n)
        end if
100   continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine plttit
c-----------------------------------------------------------------------
c
c     Purpose: plot title on screen
c
c     Inputs:
c
c     Outputs:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE bdata
      character str*76
      write(str,1000) (head(i),i=2,19)
      call drawtxt(1,0.01d0,0.03d0,1,1,76,str)
1000  format(19a4)
      end
c
c-----------------------------------------------------------------------
c
      subroutine plttxt(v3,nsizt)
c-----------------------------------------------------------------------
c
c     Purpose: Plot text on screen
c
c     Inputs:
c      v3 .ne. 0 gives 15 textlines in right window
c      v3:       line number (from top to bottom)
c
c     Outputs:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE iofile
      USE iwinio
      USE pdata2
      USE pdatap
      USE pftn77
      USE plotter
      implicit double precision (a-h,o-z)
      character*1 tx(80)
cww   character*11 strph
      integer*2 bstat,ih1,ih2,iyy,iyyy,ih3,ih4
      integer*4 status
      real*4 pos(3),xpos,ypos
      if(idev.eq.4) then
        iyy  = 0.9375 *iwxgs*0.75
        iyyy = 0.058455*iwygs !row 28
        ih1  = iwygs
        ih2  = iwygs*0.5
        ih3  = iwxgs*0.5*0.75
        ih4  = 0.9*iwygs
cww    call clwopen('TEXT  Input',1,iwys-120,600,150,1,2)
        call clwopen('TEXT  Input',1,2)
      end if
c.... request text inputs
      write(*,2030)
      read(*,1000) tx
c.... compute length of text record
      do 120 nn = 80,1,-1
        if(tx(nn).ne.' ') go to 130
120   continue
130   continue
      if(idev.eq.1) then
c.....  get mouse position  IBM
        if(v3.gt.0) then
          pos(1) = 0.54
          pos(2) = 0.27 - 0.07 * v3
        else
c.....    inquire locator device state (15.42)
          call gpqlc(1,1,1,80,err,mode,esw,view,pos,echo,area,d1,data)
c.....    initialize locator (new to define crosshair(type 2)) (9.30)
          call gpinlc(1,1,0,pos,2,area,0,data)
c.....    request locator (get mouse position after click)
          write(*,2004)
          call gprqlc(1,1,status,view,pos)
        end if
        xp(3) = pos(1)
        xp(4) = pos(2)
        if(nexte.gt.0) then
          xxp(2) =  0.64d0*(1.d0+pos(1))
          yyp(2) =  0.64d0*(1.d0+pos(2))
        end if
      else if(idev.eq.2) then
c.....  get mouse position  GKS
        if(v3.gt.0) then
          xpos = (1.01/1.28)
          ypos = (0.76-0.04*v3)/1.28
        else
c.....    request locator (get mouse position after click)
          write(*,2004)
          call grqlc(1,1,status,it,xpos,ypos)
        end if
        xpp(2) = xpos
        ypp(2) = ypos
        if(nexte.gt.0) then
          xxp(2) =  xpp(2)
          yyp(2) =  ypp(2)
        end if
      else if(idev.eq.3) then

      else if(idev.eq.4) then
        if(v3.gt.0) then
          ixa =  0.7656*iwxgs
          iya =  0.1875*iwygs + 0.04375*iwygs  * v3
        else
c.......  get mouse position  SALFORD
          write(*,2004)
140       call get_mouse_position(ixa,iya,bstat)
          if(bstat.eq.0) goto 140
          call get_mouse_position(ixa,iya,bstat)
        end if
        if(nexte.gt.0) then
          xxp(2) =        ixa/(iwxgs*0.75)
          yyp(2) = 1.d0 - iya/(iwygs*1.00)
        end if
      end if
c.....plot text at defined place
cww   call pltsiz(nsizt)
      call tplot(tx,nn)
cww   call pltsiz(1)
      if(idev.eq.4) call clwclose(1,2,0)
      return
c
1000  format(80a1)
2004  format(' Position mouse to location of first character',
     +       ' and press left button.')
2030    format('Enter text > ',$)
      end
c
c-----------------------------------------------------------------------
c
      subroutine plttraj(x,st,ndm,numnp,tra,vr,rk2,k3,k1)
c-----------------------------------------------------------------------
c
c     Purpose: plot trajectories of main-stresses
c
c     Inputs:
c       k1 = no 0f trajectory
c       k3=0,2: [ S_x=st(*,1), S_xy=st(*,2), S_y = st(*,3) ]
c       stresses are in 1-2 plane
c       sig(1) = S_x    sig(4) = S_1
c       sig(2) = S_xy   sig(5) = S_2
c       sig(3) = S_y    sig(6) = alpha
c
c       k3=3: [ S_xx=st(*,1), S_yy=st(*,2), S_zz = st(*,3) ]
c           : [ S_xy=st(*,4), S_xz=st(*,2), S_yz = st(*,3) ]
c
c
c     Outputs:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata1
      USE pdatas
      implicit double precision (a-h,o-z)
      dimension dd(2),xx(3,2),x(ndm,*),st(numnp,*),sig(6),tra(3,3),
     1          vr(3),sig3(6),sig0(3),z(3,3)
      pi=datan(1.d0)*4.d0
c...  typ 0,2=2D,3=3D
      if(k3.eq.0) k3 = 2
      k3 = max(2,min(k3,3))
      fm = .015/scale*rk2
c.... compute vector at each node
      do 120 n = 1,numnp
        if(k3.eq.2) then
c....     2D-case
          call pzero(xx,6)
          if(ndm.gt.2) then
            do i = 1,2
              xx(3,i) = x(3,n)
            end do
          end if
c....     calculate
          sig(1) = st(n,1)
          sig(2) = st(n,2)
          sig(3) = st(n,3)
          call pstres(sig,sig(4),sig(5),sig(6))
          alfa = pi*sig(6)/180.d0
          sig0(1) = sig(4)
          sig0(2) = sig(5)
c
c....     choose values to plot
          ia = 1
          ie = 2
          if(k1.eq.1) then
            ia=1
            ie=1
          else if(k1.eq.2) then
            ia=2
            ie=2
          end if
c....     plot
          do ii = ia,ie
            if(ii.eq.2) alfa = + alfa + pi/2.d0
            dd(1) = fm*dcos(alfa)
            dd(2) = fm*dsin(alfa)
c.....      modify due to quadrant for symmetry
            if(isym2(iadd+1,1,its).eq.1) dd(1) = - dd(1)
            if(isym2(iadd+1,2,its).eq.1) dd(2) = - dd(2)
c....       plot
            if(sig0(ii).ge.0.d0) then
              call pppcol(3) ! blue
            else
              call pppcol(2) ! red
            end if
            xx(1,2) = x(1,n)  - dd(1)/2.d0
            xx(2,2) = x(2,n)  - dd(2)/2.d0
            xx(1,1) = xx(1,2) + dd(1)
            xx(2,1) = xx(2,2) + dd(2)
c....       transform the vector
            call plxtrn(xx,tra,vr,ndm,2)
c....       plot the vector
            call plotl(xx(1,1),xx(2,1),xx(3,1),3)
            call plotl(xx(1,2),xx(2,2),xx(3,2),2)
          end do
        else ! 3D-case
c....     calculate
          sig3(1) = st(n,1)
          sig3(2) = st(n,2)
          sig3(3) = st(n,3)
          sig3(4) = st(n,4)
          sig3(5) = st(n,5)
          sig3(6) = st(n,6)
          call pstres1(sig3,sig0,z)
c....     choose values to plot
          ia = 1
          ie = 3
          if(k1.eq.1) then
            ia=1
            ie=1
          else if(k1.eq.2) then
            ia=2
            ie=2
          else if(k1.eq.3) then
            ia=3
            ie=3
          end if
c....     plot
          do ii = ia,ie
            call pzero(xx,6)
            dds = fm
            xx(1,1) = x(1,n)  - dds*z(1,ii) ! begin
            xx(2,1) = x(2,n)  - dds*z(2,ii)
            xx(3,1) = x(3,n)  - dds*z(3,ii)
            xx(1,2) = x(1,n)  + dds*z(1,ii) ! end
            xx(2,2) = x(2,n)  + dds*z(2,ii)
            xx(3,2) = x(3,n)  + dds*z(3,ii)
c....       transform the vector
            call plxtrn(xx,tra,vr,ndm,2)
c....       plot the vector
            if(sig0(ii).ge.0.0d0) then
              call pppcol(3) ! blue
            else
              call pppcol(2) ! red
            end if
            call plotl(xx(1,1),xx(2,1),xx(3,1),3)
            call plotl(xx(1,2),xx(2,2),xx(3,2),2)
          enddo
        end if
120   continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pltwip(k1)
c-----------------------------------------------------------------------
c
c     Purpose: wipe part of screen via plot in black
c
c     Inputs:
c       use mouse to define bounds
c       point 1 ix1,iy1  left top
c       point 2 ix2,iy2  right bottom
c       k1 = 1  between points
c       k1 = 2  outside points
c
c     Outputs:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE hptext1
      USE iwinio
      USE pdata2
      USE pftn77
      implicit double precision(a-h,o-z)
      integer*2 ihg,ivg,ix1,iy1,ix2,iy2,bstat
      real*8 x(2,11)
      integer*4 status
      real*4 pos(3),area(6),xpos,ypos
      integer data(20)
      call plopen
      if(idev.eq.1) then
c.....  inquire locator device state (15.42)
        call gpqlc(1,1,1,80,err,mode,esw,view,pos,echo,area,d1,data)
c....   get point 1
        write(*,2000)
c.....  initialize locator (new to define crosshair(type 2)) (9.30)
        call gpinlc(1,1,0,pos,2,area,0,data)
c.....  request locator (get mouse position after click)
        call gprqlc(1,1,status,view,pos)
        x(1,9) = (1.+ pos(1))*0.641
        x(2,9) = (1.+ pos(2))*0.641
c....   get point 2
        write(*,2001)
c.....  initialize locator (new to define rubberband(type 5))
        call gpinlc(1,1,0,pos,5,area,80,data)
c.....  request locator (get mouse position after click)
        call gprqlc(1,1,status,view,pos)
        x(1,7) = (1.+ pos(1))*0.641
        x(2,7) = (1.+ pos(2))*0.641
c.....  initialize locator (set back to standard values)
        call gpinlc(1,1,0,pos,1,area,0,data)
      else if(idev.eq.2) then
c....   get point 1
        write(*,2000)
c.....  request locator (get mouse position after click)
        call grqlc(1,1,status,it,xpos,ypos)
        x(1,9) = xpos*1.28
        x(2,9) = ypos*1.28
c....   get point 2
        write(*,2001)
c.....  initialize locator (new to define rubberband(type 5))
        call ginlc(1,1,0,xpos,ypos,5,0,1,0,1,80,data)
c.....  request locator (get mouse position after click)
        call grqlc(1,1,status,it,xpos,ypos)
        x(1,7) = xpos*1.28
        x(2,7) = ypos*1.28
c.....  initialize locator (set back to standard values)
        call ginlc(1,1,0,xpos,ypos,2,0,1,0,1,80,data)
      else if(idev.eq.3) then

      else if(idev.eq.4) then
CWW     call clwopen('WIPE   Parameters',1,iwys-120,640,150,1,2)
        call clwopen('WIPE   Parameters',1,2)
c....   set start position of cursor at midpoint
        ihg = 0
        ivg = 0
        write(*,2002)
c....   move cursor to point 1
        call set_graphics_selection(1)
140     call get_mouse_position(ihg,ivg,bstat)
        if(bstat.eq.0) goto 140 ! no button depressed
        ix1   = ihg
        iy1   = ivg
        x(1,9) =     ix1/(iwxgs*0.75)
        x(2,9) = 1.- iy1/(iwygs*1.00)
c....   move cursor to point 2
150     call get_mouse_position(ihg,ivg,bstat)
        if(bstat.eq.1) goto 150 ! left button depressed
        call get_mouse_position(ihg,ivg,bstat)
        call set_graphics_selection(0)
        ix2   = ihg
        iy2   = ivg
        x(1,7) =     ix2/(iwxgs*0.75)
        x(2,7) = 1.- iy2/(iwygs*1.00)
        call clwclose(1,2,0)
      end if
      call pppcol(32)
      if(k1.eq.1) then
c....   clear area between points
        if(idev.lt.4) then
          x0 = x(1,9)
          y0 = x(2,7)
          dx0= x(1,7)-x(1,9)
          dy0= x(2,9)-x(2,7)
          call ppbox(x0,y0,dx0,dy0,ipgl)
        else if(idev.eq.4) then
          call clear_screen_area(ix1,iy1,ix2,iy2,icc)
        end if
      else if(k1.eq.2) then
c....   clear area ouside points
c       4------------------ 3
c       |  9 ---------- 8   |
c       |  |            |   |
c       |  6(10) ------ 7   |
c       |                   |
c       1(5)(11)----------- 2
c
        x(1,1) = 0.004
        x(1,4) = x(1,1)
        x(1,5) = x(1,1)
        x(1,11)= x(1,1)
        x(1,6) = x(1,9)
        x(1,10)= x(1,9)
        x(1,8) = x(1,7)
        x(1,2) = 0.968
        x(1,3) = x(1,2)
        x(2,1) = x(1,1)
        x(2,2) = x(1,1)
        x(2,5) = x(1,1)
        x(2,11)= x(1,1)
        x(2,6) = x(2,7)
        x(2,10)= x(2,7)
        x(2,8) = x(2,9)
        x(2,3) = x(1,2)
        x(2,4) = x(2,3)
        call dplot(x(1,1),x(2,1),ipgl,0)
        do i = 2,10
          call dplot(x(1,i),x(2,i),2,0)
        enddo
        if(ipgl.eq.1) call clpan
c....   plot frame around clipped area
c       9 ---------- 8
c       |             |
c       6(10) ------ 7
c       call pppcol(1)
c       call dplot(x(1,6),x(2,6),3,0)
c       do i = 7,10
c         call dplot(x(1,i),x(2,i),2,0)
c       enddo
      end if
      return
2000  format(' mark left top point with left mouse button')
2001  format(' mark right bottom point with left mouse button')
2002  format(' Mark 1. wipe point with left mouse button,',
     +       ' hold button depressed and move ',//,
     +       ' mouse cursor to 2. wipe point.',
     +       ' Loose mouse button for chosen rectangle.')
      end
c
c-----------------------------------------------------------------------
c
      subroutine plxtrn(x,tr,vr,ndm,numnp)
c-----------------------------------------------------------------------
c
c     Purpose: transform the coordinates by the current transformation
c
c     Inputs:
c       x(ndm,numnp) - coordinates
c       tr(3,3)      - transformation matrix
c       vr(3)        - translation vector
c
c     Outputs:
c       x(ndm,numnp) - new coordinates
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension x(ndm,*),tr(3,3),vr(3),r(3)
c.... transform the coordinates by the current transformation
      if(ndm.lt.3) return
      do n = 1,numnp
        if(x(1,n).ne. -999.d0) then
          do i = 1,3
            r(i) = tr(1,i)*x(1,n)+tr(2,i)*x(2,n)+tr(3,i)*x(3,n)+vr(i)
          end do
          x(1,n) = r(1)
          x(2,n) = r(2)
          x(3,n) = r(3)
        end if
      end do
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine ppbox(x,y,dx,dy,is)
c-----------------------------------------------------------------------
c
c     Purpose: plot a contour or filled box
c
c     Inputs:
c       x
c       y
c       dx
c       dy
c       is = 1 - fill box
c       is = 3 - contour box
c
c     Outputs:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
c
      integer is
      real*8 x, y, dx, dy
c
      call dplot(x   ,y   ,is,0)
      call dplot(x+dx,y   ,2,0)
      call dplot(x+dx,y+dy,2,0)
      call dplot(x   ,y+dy,2,0)
      call dplot(x   ,y   ,2,0)
      if(is.eq.1) call clpan
      end
c
c-----------------------------------------------------------------------
c
      subroutine ppbox2(x,y,ds)
c-----------------------------------------------------------------------
c
c     Purpose: plot a square filled box
c
c     Inputs:
c       x
c       y
c       ds
c
c     Outputs:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
c
      real*8 x, y, ds
c
c.... fill
      call pppcol(3)
      call plotl(x-ds , y+ds , 0.d0, 3)
      call plotl(x-ds , y+ds , 0.d0, 1)
      call plotl(x-ds , y-ds , 0.d0, 2)
      call plotl(x+ds , y-ds , 0.d0, 2)
      call plotl(x+ds , y+ds , 0.d0, 2)
      call plotl(x-ds , y+ds , 0.d0, 2)
      call clpan

c.... outl
      call pppcol(2)
      call plotl(x-ds , y+ds , 0.d0, 3)
      call plotl(x-ds , y-ds , 0.d0, 2)
      call plotl(x+ds , y-ds , 0.d0, 2)
      call plotl(x+ds , y+ds , 0.d0, 2)
      call plotl(x-ds , y+ds , 0.d0, 2)

      end
c
c-----------------------------------------------------------------------
c
      subroutine pdefm(x,b,c,angl,ndm,ndf,numnp, dr,xfac)
c-----------------------------------------------------------------------
c
c     Purpose: compute the deformed mesh
c
c     Inputs:
c
c     Outputs:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE mdat2
      USE pdata7
      USE pdatas
      implicit double precision (a-h,o-z)
      logical isNurbs ! isogeo
      dimension x(ndm,*),b(ndf,*),angl(*),uu(9),dr(ndm,*),xfac(3)
      pi=datan(1.d0)*4.d0
cwd.. for IGA ndm inside the function has the be the geometrical value
      if (isNurbs()) then 
        ndm_geo = ndm - 1
      else 
        ndm_geo = ndm
      end if
c...  dofs for angl
      ij1 = ia(1)
      ij2 = ia(2)
      do 120 n = 1,numnp
        if(x(1,n).ne. -999.d0) then
          if(ipla.eq.1) then
c.....      plate elements (ndm=3!)
            uu(1)   = c*b(1,n)
            dr(1,n) = x(1,n) * xfac(1)
            dr(2,n) = x(2,n) * xfac(2)
            if(isym2(iadd+1,1,its).eq.1) uu(1)= -uu(1)
            dr(3,n) = (x(3,n) + uu(1))*xfac(3)
          else
            uu(1) = b(1,n)
            if(ndm_geo.gt.1) then
              if(angl(n).ne.0.0d0) then
                ang = angl(n)*pi/180.d0
                cn  = cos(ang)
                sn  = sin(ang)
              else
                cn  = 1.0
                sn  = 0.0
              end if
              do 100 i = 1,ndm_geo
                if(i.le.ndf) then
                  uu(i) = b(i,n)
                else
                  uu(i) = 0.0
                end if
 100          continue
              ut      = uu(ij1)*cn - uu(ij2)*sn
              uu(ij2) = uu(ij1)*sn + uu(ij2)*cn
              uu(ij1) = ut
            end if
c           do 110 i = 1,ndm
c               dr(i,n) = x(i,n) + c*uu(i)
c110        continue
cwd         changes for isogeometric version
            if (isNurbs()) then
                dr(1,n) = ( x(1,n) + c*uu(1) ) * xfac(1)
                dr(2,n) = ( x(2,n) + c*uu(2) ) * xfac(2)
                if(ndm_geo.ge.3) dr(3,n) = ( x(3,n)+c*uu(3))* xfac(3)
                if(ndm.eq.4) dr(4,n) = x(4,n)
            else
c               standard version
                dr(1,n) = ( x(1,n) + c*uu(1) ) * xfac(1)
                dr(2,n) = ( x(2,n) + c*uu(2) ) * xfac(2)
                if(ndm.eq.3) dr(3,n) = ( x(3,n) + c*uu(3) ) * xfac(3)
            end if
cwd         end of changes
          end if
        else
          dr(1,n) = -999.d0
        end if
120   continue
      return
      end
c
      subroutine prinstr(st,mpris,numnp)
c-----------------------------------------------------------------------
c
c     Purpose: calculate principal/equivalent stresses
c
c     Inputs:
c       st(numnp,npstr)  - stresses at nodes
c       mpris            - number of principal stress to plot
c       nprip            - position of stresses, default set in pmesh
c                          reset in element possible
c       nptyp            - similar
c
c     Outputs:
c       st(numnp,26)     - principal/eqivalent stress to plot
c
c     nptyp stresses                     |1 |2 |3      |4 |5 |6    |7   |8   |9    |10
c     1     sx,sxy,sy,                -> |s1|s2|alpha_s|  |  |     |sv2s|    |     |
c     2     nx,nxy,ny,                -> |n1|n2|alpha_n|  |  |     |sv2n|    |     |
c     3     mx,mxy,my                 -> |m1|m2|alpha_m|  |  |     |sv2m|    |     |
c     4     sx,sy,sz,sxy,sxz,syz      -> |s1|s2|s3     |  |  |     |    |    |     |sv3
c     5     nx,nxy,ny,mx,mxy,my,qx,qy -> |n1|n2|alpha_n|m1|m2|alpha|sv2n|sv2m|sv2.5|
c
c
c     W. Wagner BS KIT 06/10
c-----------------------------------------------------------------------
      USE prisdat
      USE strnam
      implicit double precision (a-h,o-z)
      dimension st(numnp,*),sig(8),z(3,3),sig0(3)

c.... 2D/3D name
      if(nptyp.le.3)  npris = 3
      if(nptyp.eq.4)  npris = 6
      if(nptyp.eq.5)  npris = 8

      strsus(26)= ' '

      if(nptyp.eq.1) then
        if(mpris.eq.1) strsus(26) = ' Stress S_1    '
        if(mpris.eq.2) strsus(26) = ' Stress S_2    '
        if(mpris.eq.3) strsus(26) = '   alpha_s1    '
        if(mpris.eq.7) strsus(26) = '  Sigma_v(s)   '
      else if(nptyp.eq.2) then
        if(mpris.eq.1) strsus(26) = '  Force n_1    '
        if(mpris.eq.2) strsus(26) = '  Force n_2    '
        if(mpris.eq.3) strsus(26) = '   alpha_n1    '
        if(mpris.eq.7) strsus(26) = '  Sigma_v2(n)  '
      else if(nptyp.eq.3) then
        if(mpris.eq.1) strsus(26) = '  Moment m_1   '
        if(mpris.eq.2) strsus(26) = '  Moment m_2   '
        if(mpris.eq.3) strsus(26) = '   alpha_m1    '
        if(mpris.eq.7) strsus(26) = '  Sigma_v2(m)  '
      elseif(nptyp.eq.4) then
        if(mpris.eq.1) strsus(26) = ' Stress S_1    '
        if(mpris.eq.2) strsus(26) = ' Stress S_2    '
        if(mpris.eq.3) strsus(26) = ' Stress S_3    '
        if(mpris.eq.10)strsus(26) = '  Sigma_v3(s)  '
      else if(nptyp.eq.5) then
        if(mpris.eq.1) strsus(26) = '  Force n_1    '
        if(mpris.eq.2) strsus(26) = '  Force n_2    '
        if(mpris.eq.3) strsus(26) = '   alpha_n1    '
        if(mpris.eq.4) strsus(26) = '  Moment m_1   '
        if(mpris.eq.5) strsus(26) = '  Moment m_2   '
        if(mpris.eq.6) strsus(26) = '   alpha_m1    '
        if(mpris.eq.7) strsus(26) = '  Sigma_v2(n)  '
        if(mpris.eq.8) strsus(26) = '  Sigma_v2(m)  '
        if(mpris.eq.9) strsus(26) = '  Sigma_v25(m) '
      end if


c.... loop over nodes
      do ii = 1,numnp
c....   set stresses
        do is = 1,npris
          sig(is) = st(ii,nprip(is))
        end do

        if(mpris.le.6) then
c....     principal stresses
          mmpris = mpris
          if(nptyp.eq.5) then
            if(mpris.le.3) call pstres(sig,sig0(1),sig0(2),sig0(3))! 5
            if(mpris.gt.3) then
              call pstres(sig(4),sig0(1),sig0(2),sig0(3))! 5
              mmpris = mpris-3 ! for correct storage
            end if
          else
            if(npris.eq.3) call pstres(sig,sig0(1),sig0(2),sig0(3))! 1,2,3
            if(npris.eq.6) call pstres1(sig,sig0,z) ! 4
          end if

c....     store
          st(ii,26) = sig0(mmpris)
        else if(mpris.gt.6) then
c....     eqivalent stresses
          if(mpris.eq.7) then
c....       eqivalent stress 2D s,n,m
            st(ii,26) = sigv2d(sig) 

          else if(mpris.eq.8) then
c....       eqivalent stress 2D m
            st(ii,26) = sigv2dm(sig) 
 
          else if(mpris.eq.9) then
c....       eqivalent stress 2.5D m,q
            st(ii,26) = sigv25d(sig) 

          else if(mpris.eq.10) then
c....       eqivalent stress 3D s
            st(ii,26) = sigv3d(sig) 

          end if
        end if
      end do
      return
      end
c
      double precision function sigv2d(sig)
c----------------------------------------------------------------------
c     Sigma_V2D  fo Sigma, N,M 
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension sig(*)
      sigv2d = sqrt(
     + sig(1)*sig(1)+sig(3)*sig(3)-sig(1)*sig(3)+ 3.d0*sig(2)*sig(2))
      return
      end
c
      double precision function sigv2dm(sig)
c----------------------------------------------------------------------
c     Sigma_V2D  for M 
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension sig(*)
      sigv2dm = sqrt(
     + sig(4)*sig(4)+sig(6)*sig(6)-sig(4)*sig(6)+ 3.d0*sig(5)*sig(5))
      return
      end
c
      double precision function sigv25d(sig)
c----------------------------------------------------------------------
c     Sigma_V2.5D  for M,Q 
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension sig(*)
      sigv25d = sqrt(
     + sig(4)*sig(4)+ sig(6)*sig(6) - sig(4)*sig(6)
     +              + 3.d0*sig(5)*sig(5) + 3.d0*sig(7)*sig(7)
     +              + 3.d0*sig(8)*sig(8))
      return
      end
c
      double precision function sigv3d(sig)
c----------------------------------------------------------------------
c     Sigma_V3D  fo Sigma
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension sig(*)
      sigv3d = sqrt(
     +         sig(1)*sig(1)+sig(2)*sig(2)+sig(3)*sig(3)
     +                      -sig(1)*sig(2)-sig(1)*sig(3)-sig(2)*sig(3)
     +       + 3.d0*sig(4)*sig(4)+3.d0*sig(5)*sig(5)+3.d0*sig(6)*sig(6))
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine psymm(a,nn,numnp,isymm,ndm,is,its)
c-----------------------------------------------------------------------
c
c     Purpose: compute coordinates + displ.(1--3) for adding quadrants
c              in case of symmetry
c
c     Inputs:
c       isym1: modify quadrant after quadrant
c       isym2: modify quadrant always with respect to 1. quadrant
c
c     Outputs:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension a(nn,numnp),isymm(8,3,3)
      if(is.eq.0) return
      mdm = 2
      if(ndm.eq.3) mdm = 3
      do 10 idm = 1,mdm
         if(isymm(is,idm,its).eq.1) then
            do 20 m = 1,numnp
               if(a(idm,m).ne.-999.d0) a(idm,m) = - a(idm,m)
20          continue
         end if
10    continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pwind(x,dr,ndm,ndf,numnp)
c-----------------------------------------------------------------------
c
c     Purpose:
c
c     Inputs:
c
c     Outputs:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension x(ndm,numnp),dr(ndf,numnp),xx(3,2)
      do 100 i = 1,ndm
        xx(i,1) = x(i,1)
        xx(i,2) = x(i,1)
100   continue
      do 110 n = 1,numnp
      do 110 i = 1,ndm
        xx(i,1) = min(xx(i,1),x(i,n),dr(i,n))
        xx(i,2) = max(xx(i,2),x(i,n),dr(i,n))
110   continue
      do 120 i = 1,ndm
        dr(i,1)     = xx(i,1)
        dr(i,numnp) = xx(i,2)
120   continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pzerol(fl,val,nn)
c-----------------------------------------------------------------------
c
c     Purpose: set a logical field fl to val
c
c     Inputs:
c
c     Outputs:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      logical fl(nn),val
      do 100 n = 1,nn
          fl(n) = val
100   continue
      end
c
c-----------------------------------------------------------------------
c
      subroutine rprint(dr,ix,x,ndm,numnp,ndf,nen1,nfl,idev)
c-----------------------------------------------------------------------
c
c     Purpose: compute the profile of values for stre,cont etc.
c              values are in dr(numnp)
c     Inputs:
c       dr(ndf,numnp)- array to plot is dr(1,n) n=1..numnp
c       ix(nen1,*)   - Element nodal connections of mesh
c       x(ndm,numnp) - coordinates
c         ndm        - Spatial dimension of mesh
c         numnp      - Number of nodes in mesh
c         ndf        - 1
c         nen1       - Dimension for ix array
c         nfl        - no of different values to plot
c         idev       - Graphic type
c
c     Outputs:       - rmx, rmn, pr
c
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE contval
      USE iofile
      USE iwinio
      USE pdata7
      USE rndata
      USE rpdata
      implicit double precision (a-h,o-z)
      dimension dr(ndf,*),x(ndm,numnp),pr(9)
      integer ix(nen1,*)
      data blank /-999.d0/
      save  im
      if(icv.eq.1.and.idev.eq.4.and.ior.lt.0)
cww     + call clwopen('Values of    PROFILE',1,iwys-190,630,230,1,2)
     + call clwopen('Values of    PROFILE',1,2)
      call pzero (pr,9)
      drv = 100./numnp

      if(imuse().eq.0) then
c....   plot in the range of all material numbers
        rmx = dr(1,1)
        rmn = dr(1,1)
        nmn = 1
        nmx = 1

        do 100 n = 1,numnp
          if(x(1,n).eq.blank) go to 100

          if(rmx.lt.dr(1,n)) then
            rmx = dr(1,n)
            nmx = n
          end if
          if(rmn.gt.dr(1,n)) then
            rmn = dr(1,n)
            nmn = n
          end if
100     continue

      else
c....   plot only in the range of specified material numbers
        rmx = 0.d0
        rmn = 0.d0
        nmn = 0
        nmx = 0
c....   min,max
        im = 0
        do 110 n = 1,numnp
          if(iplmano(ix,n,nen1).eq.0)  goto 110

          if(x(1,n).eq.blank) go to 110

          if(im.eq.0) then
            rmx = dr(1,n)
            rmn = dr(1,n)
            nmn = n
            nmx = n
            im = 1
          else
            if(rmx.lt.dr(1,n)) then
              rmx = dr(1,n)
              nmx = n
            end if
            if(rmn.gt.dr(1,n)) then
              rmn = dr(1,n)
              nmn = n
            end if
          end if
110     continue
      end if

      drm =  dabs(rmx-rmn)
      if(drm.lt.1.d-5*dabs(rmx).or.drm.lt.1.e-10) then
c....   nearly same values, no profile and plot, cont-->fill
          if(icv.eq.1) then
            if(ior.ge.0) write(iow,2001) rmn,rmx
            if(ior.lt.0) write(*  ,2001) rmn,rmx
          end if
cww          nfl = 0
      else
c....   profile
        do 200 n = 1,numnp
          if(x(1,n).eq.blank) go to 200
          rs = (dr(1,n) - rmn)/(rmx - rmn)
          do i = 1,9
            if(rs.ge.0.1d0*i) pr(i) = pr(i) + drv
          end do
200     continue
        if(icv.eq.1) then
          if(ior.ge.0) write(iow,2000) rmn,nmn,rmx,nmx,pr
          if(ior.lt.0) write(*  ,2000) rmn,nmn,rmx,nmx,pr
        end if
      end if
      return
2000  format(' Minimum is ',e15.5,' at node ',i6/,
     1       ' Maximum is ',e15.5,' at node ',i6,/,
     2  20x,'10%   20%   30%   40%   50%   60%   70%   80%   90%'/
     3       ' Profile above is:',10f6.1)
2001  format(' WARNING no profile, Minimum ',e15.5,' Maximum ',e15.5)
      end
c
c-----------------------------------------------------------------------
c
      subroutine rprint1(dr,ix,x,ndm,numnp,ndf,nfl,nen1,ii,idev)
c-----------------------------------------------------------------------
c
c     Purpose: compute the profile of values for stre,cont etc.
c              copy row ii values from dr(ii,numnp) to dr(numnp)
c
c     Inputs:
c       ii           - dof to plot in dr
c       other values - see SR rprint
c
c     Outputs:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension dr(ndf,*),x(ndm,*)
      integer ix(nen1,*)
      call rprint(dr(ii,1),ix,x,ndm,numnp,ndf,nen1,nfl,idev)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine setfor(f,f0,prop,nn, dr)
c-----------------------------------------------------------------------
c
c     Purpose: set forces to actual value
c
c     Inputs:
c
c     Outputs:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension f(*),f0(*),dr(*)
      do 100 n = 1,nn
         dr(n) = dr(n) + f(n)*prop + f0(n)
100   continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine tplot(tx,nn)
c-----------------------------------------------------------------------
c
c     Purpose: write text on screen
c
c     Inputs:
c
c     Outputs:
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata2
      USE pdatap
      USE pftn77
      USE plotter
      implicit double precision (a-h,o-z)
      character*1 tx(nn),textph*80
      write(textph,'(80a1)') (tx(i),i=1,nn)
      if(idev.eq.1)   call gptx2(xp(3),80,textph)
      if(idev.eq.2)   call gtx(xpp(2),ypp(2),textph)
      if(idev.eq.3)   call draw_text(textph,ixa,iya,icc)
      if(idev.eq.4)   call draw_text(textph,ixa,iya,icc)
      if(nexte.gt.0)  then
        call hptext (xxp(2),yyp(2),textph,ihpgl)
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      logical function zoom(xl,ndm,nel)
c-----------------------------------------------------------------------
c
c     Purpose: check if an element is in the window
c
c     Inputs:
c       xl(ndm,nel) - coordinates of nodes of element
c       ndm         - dimension
c       nel         - number of nodes on element
c
c     Outputs:
c       zoom        - true/false
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE pdata1
      USE pdata4
      USE ppers
      implicit double precision (a-h,o-z)
      logical errv
      dimension xl(ndm,*),xg(3)
      zoom = .false.
      if (fact .eq. 0.0d0) fact = 1.d0
      xm = xl(1,1)
      ym = xl(2,1)
      xx = xm
      yx = ym
      zm = 0.d0
      if(ndm.ge.3) zm = xl(3,1)
      zx = zm
      do 100 n = 1,nel
        if(xl(1,n).eq.-999.d0) return
        xm = min(xm,xl(1,n))
        xx = max(xx,xl(1,n))
        ym = min(ym,xl(2,n))
        yx = max(yx,xl(2,n))
        if(ndm.ge.3) then
          zm = min(zm,xl(3,n))
          zx = max(zx,xl(3,n))
        end if
100   continue
      if(kpers.eq.1) then
        xg(1) = xm
        xg(2) = ym
        xg(3) = zm
        call perspj(xg,xg,3,3,1,errv)
        xm    = xg(1)
        ym    = xg(2)
        xg(1) = xx
        xg(2) = yx
        xg(3) = zx
        call perspj(xg,xg,3,3,1,errv)
        xx    = xg(1)
        yx    = xg(2)
      end if
      zoom = ( (xm.ge.xmin(1)/fact) .and. (xx.le.xmax(1)/fact)
     1   .and. (ym.ge.xmin(2)/fact) .and. (yx.le.xmax(2)/fact) )
      return
      end
