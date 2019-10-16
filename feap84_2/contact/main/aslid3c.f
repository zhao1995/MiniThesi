c$Id:$
      subroutine aslid3c(slid, cs0,ics, ns,lics0)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Auto surface descriptions in 2-d
c               Compute maximum number of surfaces

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_contac.h'
      include   'c_comnd.h'
      include   'cdata.h'
      include   'iofile.h'
      include   'pointer.h'

      integer             slidn
      common     /aslids/ slidn

      integer    nope,dnope,neps, ns,lics0,nd,nslid,nsld

      integer    slid(5,*), ics(*)
      real*8     cs0(nr0,n0c1:nc01,*)

      data       nope,dnope / 4, 4 /

      call cdebug0 ('    aslid3c',-1)

c     Compute maximum number of surfaces to generate

      nsld  = slidn
      nslid = 0
      do nd = 1,nsld
        nslid = max(slid(5,nd),nslid)
      end do ! nd

c     First slide line number

      nd = ns - 1

      ns = 0
100   ns = ns + 1

      call cslid3d (slid,cs0(1,n0c1,nd+ns),ics(ofssurf),
     &              ns,nope,dnope,neps,ofssurf)
      ofssurf = ofssurf + dnope*neps

      if(ns.lt.nslid) go to 100

      lics0 = ofssurf - 1

      end
