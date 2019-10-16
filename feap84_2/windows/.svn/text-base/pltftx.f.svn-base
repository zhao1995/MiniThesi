c$Id:$
      subroutine pltftx(vc,ic,mc)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add eigenvector option                           09/11/2009
c          Add 'real' and 'imaginary' labels
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Place contour fill description on right hand side
c               of plot region.

c      Inputs:
c         vc(*)     - Values of contours to plotted
c         ic        - Component to plot
c         mc        - Plot type: 1 = stress;       2 = displacement;
c                                3 = velocity;     4 = acceleration;
c                                5 = prin. stress; 6 = streamline;
c                                7 = contact var.; 8 = constitutive.
c                                9 = eigenvector

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'complx.h'
      include  'elcapt.h'
      include  'iofile.h'
      include  'pdata1.h'
      include  'pdata2.h'
      include  'pdatxt.h'
      include  'pdatps.h'
      include  'pdatap.h'
      include  'plcapt.h'
      include  'prmptd.h'
      include  'psdat1.h'
      include  'rpdata.h'
      include  'sdata.h'

      logical   pcomp
      character yy*17, strs(9)*13,slab(7)*4,cplx(2)*13,blnk*17
      integer   ic, mc,mcmax, i,j,ii, ilnsv(2), iln(2)
      real*4    xph,yph
      real*8    xdv,dy,dinc,xleft,xright,xtext,xhead
      real*8    dcor,ycor,yphbot,yphtop,dx1,dx2,x1,x2,vc(*),yfr(4)

      save

      data      mcmax / 9 /
      integer   ipal(12)
      data      ipal  / 13, 4,12, 6,14, 3,10, 5,11, 8, 9, 2/
      integer   rpal(12)
      data      rpal  /  2, 9, 8,11, 5,10, 3,14, 6,12, 4,13/
      data strs/' S T R E S S ',' DISPLACEMENT','  VELOCITY',
     &          ' ACCELERATION',' Prin. Stress','  STREAMLINE ',
     &          ' CONTACT VAR.',' CONSTITUTIVE',' EIGENVECTOR '/
      data cplx/'   R E A L   ','  IMAGINARY  '/
      data blnk/'_________________'/
      data slab/'  1 ','  2 ','  3 ',' Ang',' I_1',' J_2',' J_3'/

c     Try some other y-positions for multiple contours

      data yfr/0.805d0,0.625d0,0.445d0,0.265d0/

c     Save line style and width

      ilnsv(1) = ilno(1)
      ilnsv(2) = ilno(2)
      iln  (1) = 0
      iln  (2) = 1
      call plline( iln )

c     PC device position parameters

      dtext  = 0.0600d0
      xhead  = 1.0000d0
      xleft  = 1.0400d0
      xright = 1.0800d0
      xtext  = 1.0900d0

c     Put box around caption region

      if(ifrm.eq.0) then
        call pppcol(-1,0)
        call ppbox(0.98d0,0.10d0,0.29d0,0.70d0,1)
        xdv  = 20.d0
        ycor = 0.75d0
      else
        xdv  = 70.d0
        ycor = yfr(ifrm)

      endif
      dy   = 1.d0/xdv
      xph  = 1./1.28
      yph  = ycor/1.28
      dcor = ycor - 1.35d0*dy
      dinc = 5.d0/(7.d0*xdv)

c     Draw color bars - with values

      do j = 1,12
        i = 13 - j
        if(psrevs) then
          ii = rpal(j)
        else
          ii = ipal(j)
        endif
        call pppcol(ii,2)
        yphbot = dcor - dinc*dble(i)
        yphtop = dcor - dinc*dble(i-1)
        call dplot( xleft,yphbot,1)
        call dplot(xright,yphbot,2)
        call dplot(xright,yphtop,2)
        call dplot( xleft,yphtop,2)
        call clpan
        call pppcol(1,1)
        if(j.gt.1) then
          yphbot = yphbot - 0.0075d0
          if(ifrm.eq.0 .or. j.eq.2 .or. j.eq.12) then
            write(yy, '(1p,1e9.2)' ) vc(j-1)
            call tplot(xtext,yphbot,yy,9,1)
          endif
        endif
      end do ! j

      if(ifrm.eq.0) then
        write(yy, '(1p,1e9.2)' ) rmx
        yphbot = dcor - 0.0075d0
        call tplot(xtext,yphbot,yy,9,1)
        write(yy, '(1p,1e9.2)' ) rmn
        yphbot = yphbot - 12.d0*dinc
        call tplot(xtext,yphbot,yy,9,1)
      endif

c     Real/Imaginary label

      if(ipc.gt.1) then
        call tplot(xhead,ycor+0.6d0*dy,blnk,13,1)
        if(cplxpl) then
          call tplot(xhead,ycor+0.6d0*dy,cplx(1),13,1)
        else
          call tplot(xhead,ycor+0.6d0*dy,cplx(2),13,1)
        endif
      endif

      dtext = 0.11d0
      call tplot(xhead,ycor,blnk,17,1)
      if(mc.eq.5) then
        if(ic.le.4) then
          write(yy,'(a13,a4)' ) strs(mc),slab(min(7,ic))
        elseif(ic.eq.5 .or. ic.eq.7) then
          write(yy,'(a13,a4)' ) '  Invariant  ',slab(min(7,ic))
        elseif(ic.eq.6) then
          write(yy,'(a17)' ) '  Mises  Stress  '
        endif
      elseif(mc.le.mcmax) then
        write(yy,'(a13,i3)' ) strs(mc),ic
        if(mc.eq.1) then
          if(.not.pcomp(ecapt(ic),'                 ',17)) then
            yy = ecapt(ic)
          endif
        endif
      else
        write(iow,'(a,i3)') ' Call to PLTFTX: MC =',mc
        call plstop()
      endif
      if(ncapt.gt.0) then
        yy       = ' '
        yy(2:16) = caption
      endif
      call tplot(xhead,ycor,yy,17,1)
      ncapt = 0

c     Draw box to contain color bars

      call dplot( xleft , dcor             , 3)
      call dplot( xright, dcor             , 2)
      call dplot( xright, dcor - 12.d0*dinc, 2)
      call dplot( xleft , dcor - 12.d0*dinc, 2)
      call dplot( xleft , dcor             , 2)

c     Add horizonal divider lines / ticks

      do i = 0,12
        call dplot(xleft        , dcor-dinc*dble(i),3)
        call dplot(xright+.005d0, dcor-dinc*dble(i),2)
      end do ! i

c     Place Minimum and Maximum values on plot

      if(ifrm.eq.0 .and. psmrk) then

c       Put min/max values on screen

        write(yy, '(6hMax = ,1p1e9.2)' ) psmx
        yphbot = dcor - 13.d0*dinc - 0.0075d0
        call tplot(1.02d0,yphbot,yy,15,1)

        write(yy, '(7hMin  = ,1p1e9.2)' ) psmn
        yphbot = dcor - 10.d0/xdv - 0.0075d0
        call tplot(1.02d0,yphbot,yy,16,1)

c       Place marker for view minimum

        if(psmmx) then
          dx1 = 0.008d0
          x1  = 1.00d0
          x2  = dcor - 10.d0/xdv - 0.0075d0

          call dplot(x1-dx1 , x2     , 3)
          call dplot(x1     , x2-dx1 , 2)
          call dplot(x1+dx1 , x2     , 2)
          call dplot(x1     , x2+dx1 , 2)
          call dplot(x1-dx1 , x2     , 2)

c         Place marker on plot for view minimum

          dx2 = .005/scale

          call plotl(xpsn(1)-dx2 , xpsn(2)     , xpsn(3), 3)
          call plotl(xpsn(1)     , xpsn(2)-dx2 , xpsn(3), 2)
          call plotl(xpsn(1)+dx2 , xpsn(2)     , xpsn(3), 2)
          call plotl(xpsn(1)     , xpsn(2)+dx2 , xpsn(3), 2)
          call plotl(xpsn(1)-dx2 , xpsn(2)     , xpsn(3), 2)

          call plotl(xpsn(1)     , xpsn(2)     , xpsn(3)-dx2, 3)
          call plotl(xpsn(1)     , xpsn(2)-dx2 , xpsn(3)    , 2)
          call plotl(xpsn(1)     , xpsn(2)     , xpsn(3)+dx2, 2)
          call plotl(xpsn(1)     , xpsn(2)+dx2 , xpsn(3)    , 2)
          call plotl(xpsn(1)     , xpsn(2)     , xpsn(3)-dx2, 2)

c         Place marker for view maximum

          x2  = dcor - 13.d0*dinc - 0.0075d0

          call dplot(x1-dx1 , x2+dx1 , 3)
          call dplot(x1+dx1 , x2-dx1 , 2)
          call dplot(x1-dx1 , x2-dx1 , 2)
          call dplot(x1+dx1 , x2+dx1 , 2)
          call dplot(x1-dx1 , x2+dx1 , 2)

c         Place marker on plot for view maximum

          call plotl(xpsx(1)-dx2 , xpsx(2)+dx2 , xpsx(3), 3)
          call plotl(xpsx(1)+dx2 , xpsx(2)-dx2 , xpsx(3), 2)
          call plotl(xpsx(1)-dx2 , xpsx(2)-dx2 , xpsx(3), 2)
          call plotl(xpsx(1)+dx2 , xpsx(2)+dx2 , xpsx(3), 2)
          call plotl(xpsx(1)-dx2 , xpsx(2)+dx2 , xpsx(3), 2)

          call plotl(xpsx(1)     , xpsx(2)+dx2 , xpsx(3)-dx2, 3)
          call plotl(xpsx(1)     , xpsx(2)-dx2 , xpsx(3)+dx2, 2)
          call plotl(xpsx(1)     , xpsx(2)-dx2 , xpsx(3)-dx2, 2)
          call plotl(xpsx(1)     , xpsx(2)+dx2 , xpsx(3)+dx2, 2)
          call plotl(xpsx(1)     , xpsx(2)+dx2 , xpsx(3)-dx2, 2)
        endif ! psmmx

      endif ! psmrk

c     Restore line style and width

      iln  (1) = ilnsv(1)
      iln  (2) = ilnsv(2)
      call plline( iln )

      end
