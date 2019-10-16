c$Id:$
      subroutine coutmp(n,cp0)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Output pair data to file

c     Inputs:
c        n       - Pair number
c        cp0(*)  - Pair properties

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'c_0.h'
      include   'c_comnd.h'
      include   'c_dict.h'
      include   'c_geom.h'
      include   'c_mate.h'
      include   'c_pair.h'
      include   'c_tole.h'
      include   'cdata.h'
      include   'iodata.h'

      integer    n, typ,opp
      real*8     cp0(nr0,n0c3:*)

      write(ios,2000) n,cis(typ(3,ndrv)),nsurf1,nsurf2,
     &                  cis(opp(3,2,ifsolm)),cp0(3,2),cp0(4,2),
     &                  tlipen,tlopen,tlouts
      if(max(nmat1,nmat2).gt.0) then
        write(ios,2001) nmat1,nmat2
      endif

c     Formats

2000  format(/'PAIR',i5/'  ',a,2i8/
     &        '  SOLM ',a,1p,2e15.7/
     &        '  TOLE,,',1p,3e15.7)

2001  format( '  MATE,,',2i8)

      end
