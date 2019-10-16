c$Id:$
      subroutine uftyplib(pstyp,nel,iel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    31/08/2008
c       1. Add set of ublknum                               19/11/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set plot date for user elements

c      Inputs:
c         pstyp   - Element topology
c         nel     - Number of element nodes
c         iel     - Element number

c      Output:
c         exord   - Number of plot
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iofile.h'
      include   'ublk1.h'

      integer    pstyp,nel,iel, uptyp

c     Reverse sign on pstyp

      uptyp   = -pstyp
      ublknum = max(1,uptyp/100)
      uptyp   = mod(uptyp,100)

c     User type 1

      if(uptyp.eq.1) then

        call uftyp01(uptyp,nel,iel)

c     User type 2

      elseif(uptyp.eq.2) then

        call uftyp02(uptyp,nel,iel)

c     User type 3

      elseif(uptyp.eq.3) then

        call uftyp03(uptyp,nel,iel)

c     User type 4

      elseif(uptyp.eq.4) then

        call uftyp04(uptyp,nel,iel)

c     User type 5

      elseif(uptyp.eq.5) then

        call uftyp05(uptyp,nel,iel)

      else

        write(  *,*) ' **ERROR** User plot number not between -1 and -5'
        write(ilg,*) ' **ERROR** User plot number not between -1 and -5'
        write(iow,*) ' **ERROR** User plot number not between -1 and -5'
        call plstop()

      endif

      end
