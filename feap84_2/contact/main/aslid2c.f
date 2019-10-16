c$Id:$
      subroutine aslid2c(slid, cs0,ics, ns,lics0)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Auto surface descriptions in 2-d
c               Slideline numbering

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_contac.h'
      include   'c_comnd.h'
      include   'iofile.h'

      integer             slidn
      common     /aslids/ slidn

      integer    jj, nope,dnope,neps, ns,lics0,nd

      integer    slid(2,*), ics(*)
      real*8     cs0(nr0,n0c1:nc01,*)

      data       nope,dnope / 2, 4 /

      call cdebug0 ('    aslid2c',-1)

c     First slide line number

      nd    = ns - 1
      lics0 = slidn

      jj = 1
      ns = 0
100   ns = ns + 1

      neps    = slid(2,jj) - jj
      call cslid2d (slid(1,jj),cs0(1,n0c1,nd+ns),ics(ofssurf),
     &         ns,nope,dnope,neps,ofssurf)
      ofssurf = ofssurf + dnope*neps

      jj = slid(2,jj) + 1
      if(jj.lt.lics0) go to 100

      lics0 = ofssurf - 1

      end
