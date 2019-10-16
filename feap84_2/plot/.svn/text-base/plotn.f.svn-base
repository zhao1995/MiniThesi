c$Id:$
      subroutine plotn(v,x,ix,jc,k,kp,nd,nen,nen1,numnp,numel,
     &                 nxy,sve,ck)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Determine location of line for dplot/splot command

c      Inputs:
c         v(nd,*)   - Values at all nodes for plot quantity
c         x(ndm,*)  - Nodal coordinates
c         ix(nen1,*)- Element nodal connections
c         jc        - Plot number
c         k         - Component number for plot
c         kp        - Rescale indicator
c         nd        - Dimension of v array
c         nen       - Number of nodes/element
c         nen1      - Dimension of ix array
c         numnp     - Number of nodes
c         numel     - Number of elements
c         ck        - Plot type character

c      Outputs:
c         none      - Plot output to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pbody.h'
      include  'pdata1.h'
      include  'pdata2.h'
      include  'pdataq.h'
      include  'pdatps.h'
      include  'wdata.h'

      logical   flge,hdcsav,noerr,linein
      character ck*1, butn*1
      integer   iwinold,jc,k,kold,kp,nd,nen,nen1,numnp,numel
      integer   ma,n,npv,nsv,ne,npt
      real*8    x1,y1,xc,yc,scalo

c     integer   ix(nen1,numel),ixo(4),nxy(250)
      integer   ix(nen1,numel),ixo(4),nxy(*),sxy(2,3)
c     real*8    x(3,numnp),v(nd,1),sl(2),sc(2),sve(500)
      real*8    x(3,numnp),v(nd,1),sl(2),sc(2),sve(*),svs(4)

      save

c     Save window number

      iwinold = iwindow

c     Set center for plot

      xc     = s0(1)
      yc     = s0(2)

c     Put up a time plot

      if(ck.eq.'t') then
        if(k.eq.0) then
          k    = kold + 1
          kold = mod(k,12)
        endif

        call plotcl(v,x,v,1,numnp,nd,numnp,jc,k,kp,xc,yc,ck)

c     Display marker and input list of screen coordinates

      else

c       Save hardcopy status

        hdcsav = hdcpy
        hdcpy  = .false.

c       Set line color = 5 (yellow)

        call pppcol(5,1)

c       Request first end of line (A)

        linein   = kp.le.0
        if(linein) then

          ixy(1,1) = idx/2
          ixy(2,1) = idy/2
          ixy(1,2) = ixy(1,1)
          ixy(2,2) = ixy(2,1)

          write(*,*) 'Define line by placing markers at ends A and B'

          do while(linein)
            write(*,*) 'Mark end A with mouse + press left button.'
            call gin(x1,y1,noerr,butn)

            if(butn.eq.'l') then

              ixy(1,2) = nint(x1*dble(idx))
              ixy(2,2) = nint(y1*dble(idy))
              sve(1)   = x1
              sve(2)   = y1

              do n = 1,2
                ixy(n,3) = ixy(n,2)
                ixo(n)   = ixy(n,2)
              end do ! n

              call pltext(x1,y1,1,'A')

              linein = .false.
            else
              write(*,*) '-->Mark end A with LEFT button only'
            endif
          end do ! linein

C         Request other end of line (B)

          linein   = .true.
          do while(linein)
            write(*,*) 'Mark end B with mouse + press left button'
            call gin(x1,y1,noerr,butn)

            if(butn.eq.'l') then

              ixy(1,3) = nint(x1*dble(idx))
              ixy(2,3) = nint(y1*dble(idy))
              sve(3)   = x1
              sve(4)   = y1
              linein = .false.

              call pltext(x1,y1,1,'B')

c             Put line on mesh between A and B

              call dplot( sve(1), sve(2),3)
              call dplot( sve(3), sve(4),2)

            else
              write(*,*) '-->Mark end B with LEFT button only'
            endif
          end do ! while

c         Check if same line

          if((ixy(1,2).eq.ixy(1,3)).and.(ixy(2,2).eq.ixy(2,3))) then
            write(*,*) 'Zero length line: Reinput line'
            return
          endif

c         Save values for replots

          do n = 1,3
            sxy(1,n) = ixy(1,n)
            sxy(2,n) = ixy(2,n)
            svs(n)   = sve(n)
          end do ! n
          svs(4)   = sve(4)

c       Recover values for replots

        else

          do n = 1,3
            ixy(1,n) = sxy(1,n)
            ixy(2,n) = sxy(2,n)
            sve(n)   = svs(n)
          end do ! n
          sve(4)   = svs(4)
        endif

c       Save current values of plot scaling quantities

        scalo = scale
        do n = 1,2
          sc(n)    = s0(n)
          sl(n)    = sx(n)
          ixo(n+2) = ixy(n,3)
        end do ! n

        npt = 0

        do n = 1,numel
          ma = ix(nen1,n)
          if((ix(nen1-1,n).ge.nreg1  .and.
     &        ix(nen1-1,n).le.nreg2) .and.
     &    ((ma.gt.0 .and. maplt.eq.0) .or. ma.eq.maplt) ) then
            flge = .false.
            call findxy(ixo,sve,x,ix(1,n),nen,scalo,sc,sl,flge)
            if(flge) then
              npt = npt + 1
              nxy(npt) = n
            endif
          endif
        end do ! n

        if(npt.eq.0) then
          write(*,*) 'No points found: Reinput line to search'
          return
        endif

c       Plot stress curves

        if(jc.le.0) return

c       Restore hard copy status

        hdcpy = hdcsav

c       Determine plot values

        npv = 0
        nsv = 2*nen + 1

        do n = 1,npt
          ne = nxy(n)
          flge = .true.
          call findxy(ixo,sve,x,ix(1,ne),nen,scalo,sc,sl,flge)
          call findpt(v,sve,ix(1,ne),nd,nen, sve(nsv),npv)
        end do ! n

c       Sort values into increasing order

        call sortpt(sve(nsv),npv)

c       Plot values on line (limit plot number to 12)

        if(k.ge.0) then
          k    = kold + 1
          kold = mod(k,12)
        endif

c       Set window for plots

        if(ck.eq.'d' .or. ck.eq.'v' .or. ck.eq.'a') then
          iwindow = 2
          call plclos()
          call plopen()
        elseif(ck.eq.'s') then
          iwindow = 3
          call plclos()
          call plopen()
        endif
        call plotcl(v,sve(nsv),sve(nsv+1),2,npv,nd,numnp,
     &              jc,k,kp,xc,yc,ck)
      endif

c     Plot axes around line

      call plota(xc,yc)

      iwindow = iwinold

      end
