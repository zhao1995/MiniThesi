c$Id:$
      subroutine coutm (cs0,cm0,cp0,ics,hic)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Extract type of surface into 'ntype'             16/01/2013
c       2. Add '0' to argument of call to setcomp           06/11/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Control for output of contact data to file

c     Inputs:

c     Outputs:
c       Data to output file for
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_comnd.h'
      include   'c_contac.h'

      integer    n, ofs, neps,dnope,nope, ntype
      integer    ics(*),hic((c_lp1+c_lp3),*)
      real*8     cs0(nr0,n0c1:nc01,*),cm0(nr0,n0c2:nc02,*)
      real*8     cp0(nr0,n0c3:nc03,*)

c     Surface outputs

      do n = 1,numcs
        ofs   = nint(abs(cs0(2,-1,n)))
        neps  = nint(abs(cs0(3,-1,n)))
        dnope = nint(abs(cs0(4,-1,n)))
        nope  = nint(abs(cs0(2, 0,n)))
        ntype = nint(abs(cs0(1, 0,n)))

        call coutms(n,ntype,ics(ofs),dnope,nope,neps)
      end do ! n

c     Pair outputs

      do n = 1,numcp
        call setcomp(n,cs0,cm0,cp0,hic,0)
        call coutmp(n,cp0)
      end do

      end
