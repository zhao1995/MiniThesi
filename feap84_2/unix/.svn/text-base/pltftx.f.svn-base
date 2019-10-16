c$Id:$
      subroutine pltftx(vc,ic,mc)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add eigenvector option                           28/08/2007
c       2. Add real and imaginary labels                    09/11/2009
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
c                                9 = eigenvector; 10 = strain
c                               11 = history no.;

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
      character yy*17, strs(11)*13,slab(7)*4,cplx(2)*13,blnk*17
      integer   ic, mc,mcmax, i,j, ilnsv(2), iln(2)
      real*4    xph,yph
      real*8    xdv,dy,xleft,xright,xtext,xhead
      real*8    ycor,yphbot,yphtop, dx1,dx2, x1,x2, vc(*),yfr(4)

      save

      data      mcmax / 11 /
      integer   ipal(12)
      data      ipal  / 13, 4,12, 6,14, 3,10, 5,11, 8, 9, 2/
      integer   rpal(12)
      data      rpal  /  2, 9, 8,11, 5,10, 3,14, 6,12, 4,13/
      data strs/' S T R E S S ',' DISPLACEMENT','  VELOCITY',
     &          ' ACCELERATION',' Prin. Stress','  STREAMLINE ',
     &          ' CONTACT VAR.',' CONSTITUTIVE',' EIGENVECTOR ',
     &          ' S T R A I N ',' HISTORY No. '/
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

c     X11 device position parameters

      dtext  = 0.0000d0
      xhead  = 1.1200d0
      xleft  = 1.0400d0
      xright = 1.0800d0
      xtext  = 1.1600d0

      if(ifrm.eq.0) then
        xdv  = 20.d0
        ycor = 0.75d0
        call pppcol(-1,0)
        call ppbox(0.98d0,0.10d0,0.30d0,0.70d0,1)
      else
        xdv  = 70.d0
        ycor = yfr(ifrm)
      endif
      dy  = 1.d0/xdv
      xph = 1./1.28
      yph = real(ycor)/1.28

c     Draw color bars - with values

      do i = 1,12
        if(psrevs) then
          j = rpal(i)
        else
          j = ipal(i)
        endif
        call pppcol(j,2)
        yphbot = ycor - 1.35d0*dy - 5.d0/7.d0*i/xdv
        yphtop = ycor - 1.35d0*dy - 5.d0/7.d0*(i-1)/xdv
        call dplot( xleft,yphbot,1)
        call dplot(xright,yphbot,2)
        call dplot(xright,yphtop,2)
        call dplot( xleft,yphtop,2)
        call clpan
        call pppcol(1,1)
        if(i.lt.12) then
          yphbot = yphbot - 0.0075d0
          if(ifrm.eq.0 .or. i.eq.1 .or. i.eq.11) then
            write(yy, '(1p,1e9.2)' ) vc(i)
            call tplot(xtext,yphbot,yy,9,1)
          endif
        endif
      end do ! i

      if(ifrm.eq.0) then
        write(yy, '(1p,1e9.2)' ) rmn
        yphbot = ycor - 1.35d0*dy - 0.0075d0
        call tplot(xtext,yphbot,yy,9,1)
        write(yy, '(1p,1e9.2)' ) rmx
        yphbot = yphbot - 60.d0/(7.d0*xdv)
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

      call dplot( xleft , ycor - 1.35d0*dy                   , 3)
      call dplot( xright, ycor - 1.35d0*dy                   , 2)
      call dplot( xright, ycor - 1.35d0*dy - 60.d0/(7.d0*xdv), 2)
      call dplot( xleft , ycor - 1.35d0*dy - 60.d0/(7.d0*xdv), 2)
      call dplot( xleft , ycor - 1.35d0*dy                   , 2)

c     Add horizonal divider lines / ticks

      do i = 0,12
        call dplot(xleft        , ycor-1.35d0*dy-5.d0/7.d0*i/xdv,3)
        call dplot(xright+.005d0, ycor-1.35d0*dy-5.d0/7.d0*i/xdv,2)
      end do ! i

c     Place Minimum and Maximum values on plot

      if(ifrm.eq.0 .and. psmrk) then

c       Place min/max label on screen

        write(yy, '(7hMin  = ,1p1e9.2)' ) psmn
        yphbot = ycor - 1.35d0*dy - 65.d0/(7.d0*xdv) - 0.0075d0
        call tplot(1.12d0,yphbot,yy,16,1)

        write(yy, '(6hMax = ,1p1e9.2)' ) psmx
        yphbot = ycor - 1.35d0*dy - 10.d0/xdv - 0.0075d0
        call tplot(1.12d0,yphbot,yy,15,1)

        if(psmmx) then

c         Place marker for view minimum

          dx1 = 0.008d0
          x1  = 1.02d0
          x2  = ycor - 1.35d0*dy - 65.d0/(7.d0*xdv) - 0.0075d0

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

          x2  = ycor - 1.35d0*dy - 10.d0/xdv - 0.0075d0

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
