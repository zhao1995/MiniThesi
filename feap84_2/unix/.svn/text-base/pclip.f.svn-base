c$Id:$
      subroutine pclip(x1,y1,ipin)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add x11(1) and y11(1) for passing arguments      02/07/2014
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Clip definitions for plots

c      Inputs:
c         x1,y1     - Point to plot
c         ipin      - Pen up/down

c      Outputs:
c         none      - Outputs through common blocks
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdata2.h'
      include  'pdatap.h'
      include  'pdataq.h'
      include  'pdatps.h'
      include  'pdatxt.h'
      include  'plclip.h'
      include  'plflag.h'
      include  'x11f.h'

      integer   is(2), kmin(2),kmax(2)
      integer   i,j,k,l, ipen,ipin
      real*4    x11(1),y11(1), x(2,2) , xps(2,2)
      real*8    x1,y1, xd, scord, tol

      save

      data      kmin / 4,1/, kmax/ 2,3/
      data      tol  /1.d-06/

c     Set pen command

      ipen = abs(ipin)

c     Check pen command format

      if( ipen.eq.1 ) then

        npf    = 1
        ipan   = 1
        ienter = 0
        iexit  = 0

      elseif( ipen.eq.3 ) then

c       Check that point is in plot area

        if(((x1.lt.wmax(1)) .and. (x1.gt.wmin(1))) .and.
     +     ((y1.lt.wmax(2)) .and. (y1.gt.wmin(2)))) then
          clip = .true.
        else
          clip = .false.
        end if

c     Line drawing mode

      elseif(ipen.eq.2) then
        x(2,1) = real(x1)
        x(2,2) = real(y1)

c       Check that line can be in plot area

        if(((min(x(1,1),x(2,1)).lt.wmax(1))  .and.
     +      (max(x(1,1),x(2,1)).gt.wmin(1))) .and.
     +     ((min(x(1,2),x(2,2)).lt.wmax(2))  .and.
     +      (max(x(1,2),x(2,2)).gt.wmin(2)))) then

c         Loop over directions 1=x; 2=y

          clip  = .true.
          is(1) = 0
          is(2) = 0
          do i = 1,2
            if(clip) then
              if(abs(x(1,i)-x(2,i)) .gt. tol ) then
                j = mod(i,2) + 1

c               Loop over nodes to do a clip

                do k = 1,2
                  if(clip) then
                    l  = mod(k,2) + 1
                    xd = x(k,i)
                    if(xd.lt.wmin(i)) then
                      scord  = (wmin(i)-xd)/(x(l,i)-x(k,i))
                      if(scord.gt.0.0d0 .and. scord.le.1.0d0) then
                        is(k)  = kmin(i)
                        x(k,i) = real(wmin(i))
                        x(k,j) = x(k,j) + real(scord*(x(l,j) - x(k,j)))
                      else
                        clip = .false.
                      end if
                    elseif(xd.gt.wmax(i)) then
                      scord  = (wmax(i)-xd)/(x(l,i)-x(k,i))
                      if(scord.gt.0.0d0 .and. scord.le.1.0d0) then
                        is(k)  = kmax(i)
                        x(k,i) = real(wmax(i))
                        x(k,j) = x(k,j) + real(scord*(x(l,j) - x(k,j)))
                      else
                        clip = .false.
                      end if
                    endif
                  endif
                end do    ! k-loop
              end if
            end if
          end do    ! i-loop

c       Line is outside window

        else
          clip = .false.
        end if

      end if

c     Plot segment if clip = .true.

      if(clip) then

        xps(1,1) = x(1,1)*0.78125
        xps(1,2) = x(1,2)*0.78125
        xps(2,1) = x(2,1)*0.78125
        xps(2,2) = x(2,2)*0.78125

c       X11 clipping

        if(ipen .eq. 2) then

c         Draw lines

          if(ipan.eq.0) then
            x11(1) = xps(1,1)*xx(2)
            y11(1) =   x(1,2)*xx(3)
            if(screfl) then
              call gdx11(3,x11,y11)
            endif
            x11(1) = xps(2,1)*xx(2)
            y11(1) =   x(2,2)*xx(3)
            if(screfl) then
              call gdx11(4,x11,y11)
            endif

c         Check to add corners for fills if necessary

          elseif(ipan.gt.0) then
            if(is(1).eq.0) then
              if(is(2).ne.0) then
                iexit = is(2)
                if(ienter.eq.0) ienter = -1
              endif
            elseif(is(1).ne.0) then
              if(ienter.eq.0) then
                ienter = is(1)
                if(is(2).ne.0) then
                  iexit = is(2)
                else
                  iexit  = 0
                end if
              elseif(is(1).ne.iexit) then
                call pfclip(iexit, is(1))
                iexit  = 0
              else
                iexit  = 0
              endif
            end if

c           Add two points as necessary

            if(ipan.gt.1) then
              if(xp(ipan).ne.xps(1,1) .or. yp(ipan).ne.xps(1,2)) then
                ipan     = ipan + 1
                xp(ipan) = xps(1,1)
                yp(ipan) = xps(1,2)
              endif
            else
              xp(ipan) = xps(1,1)
              yp(ipan) = xps(1,2)
            endif
            if(xp(ipan).ne.xps(2,1) .or. yp(ipan).ne.xps(2,2)) then
              ipan = ipan + 1
              xp(ipan) = xps(2,1)
              yp(ipan) = xps(2,2)
            endif
          end if
        elseif(ipen .eq. 3) then
          xp(1) = xp(2)
          yp(1) = yp(2)
          xp(2) = 0.78125*real(x1)
          yp(2) = 0.78125*real(y1)
          ipan  = 0
        endif
      end if

c     Postcript

      if(hdcpy .and. ipan.eq.0 .and. ipen.eq.2) then

c       Plot line segment

        call fpplps(2, xps(1,1), xps(1,2))

      endif


c     Set x-node for next plot check

      x(1,1) = real(x1)
      x(1,2) = real(y1)

      end
