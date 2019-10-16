c$Id:$
      subroutine protin(prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Data input routine for rotational update types

c      Inputs:
c         prt       - Flag, print input data if true
c         prth      - Flag, print title/header if true

c      Outputs:
c         Rotational update information
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'crotas.h'
      include  'dstars.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'sdata.h'
      include  'comblk.h'

      logical   prt,prth,errck,pinput
      integer   i,imn,imx,inc, nn
      real*8    td(4)

      save

c     Output title and prompt information

      call prtitl(prth)
      if(prt) then
        write(iow,2000)
        if(ior.lt.0) then
          write(iow,2000)
        endif
      endif
      if(ior.lt.0) then
        write(*,3000)
        call pprint('   >')
      endif

c     Input first record

100   errck = pinput(td,4)
      if(errck) go to 100
      nn = nint(td(1))
      if (nn.eq.0) return
      imn = nint(td(2))
      imx = nint(td(3))

c     All directors use same update option

      if (imn.eq.0) then

        if(prt) then
          write(iow,2002) nn
          if(ior.lt.0) then
            write(*,2002) nn
          endif
        endif
        do i = 0 , numnp-1
          mr(np(81)+i) = nn
        end do ! i
        call updrot(hr,ndf,hr(np(82)),mr(np(81)),numnp,0)

c     Nodes use different update options

      else

        imn = nint(td(2)) + starnd
        imx = nint(td(3)) + starnd
        inc = nint(td(4))
        if (inc.eq.0) inc=1
        imx = min(imx,numnp)
        imn = max(imn,1)
        if (imx.eq.0) imx=imn

        if(prt) then
          write(iow,2001) nn,imn,imx,inc
          if(ior.lt.0) then
            write(*,2001) nn,imn,imx,inc
          endif
        endif

        do i = imn,imx,inc
          mr(np(81)+i-1) = nn
        end do ! i

c       Input other records: Stop on blank record.

200     errck = pinput(td,4)
        if(errck) go to 200
        nn = nint(td(1))
        if (nn.eq.0) return

        imn = nint(td(2)) + starnd
        imx = nint(td(3)) + starnd
        inc = nint(td(4))
        if (inc.eq.0) inc = 1
        imx = min(imx,numnp)
        imn = max(imn,1)

        if(prt) then
          write(iow,2001) nn,imn,imx,inc
          if(ior.lt.0) then
            write(*,2001) nn,imn,imx,inc
          endif
        endif

        do i = imn,imx,inc
          mr(np(81)+i-1) = nn
        end do ! i

        go to 200

      endif

c     Formats

2000  format('   R o t a t i o n a l   T y p e s'//
     &      '       Type  1-Node   2-Node  Inc')

2001  format(3i10,i4)

2002  format('       All directors updated by rotational type:',i3)

3000  format(' Input: option, 1st node, 2nd node, inc. ')

      end
