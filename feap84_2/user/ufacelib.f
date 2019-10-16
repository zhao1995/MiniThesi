c$Id:$
      subroutine ufacelib(pstyp,nel,iu,ufac)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    31/08/2008
c       1. Add check for ublknum                            19/11/2008
c       2. Add 'nel' on call to ufacelib                    08/12/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set plot date for user elements

c      Inputs:
c         pstyp    - Element topology
c         nel      - Number element nodes

c      Output:
c         iu(4,*)  - Quadrilateral face node lists
c         ufac     - Number of quad faces
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iofile.h'
      include   'ublk1.h'

      integer    pstyp,nel,iu(4,*),ufac, uptyp

c     Reverse sign on typ

      uptyp   = -pstyp
      ublknum = max(1,uptyp/100)
      uptyp   = mod(uptyp,100)

c     User type 1

      if(uptyp.eq.1) then

        call uface01(uptyp,nel,iu,ufac)

c     User type 2

      elseif(uptyp.eq.2) then

        call uface02(uptyp,nel,iu,ufac)

c     User type 3

      elseif(uptyp.eq.3) then

        call uface03(uptyp,nel,iu,ufac)

c     User type 4

      elseif(uptyp.eq.4) then

        call uface04(uptyp,nel,iu,ufac)

c     User type 5

      elseif(uptyp.eq.5) then

        call uface05(uptyp,nel,iu,ufac)

      else

        write(  *,*) ' **ERROR** User plot number not between -1 and -5'
        write(ilg,*) ' **ERROR** User plot number not between -1 and -5'
        write(iow,*) ' **ERROR** User plot number not between -1 and -5'
        call plstop()

      endif

      end
