c$Id:$
      subroutine suchnach(seg,knoten,nachb,node,ix2,surpoin)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Anna Haraldsson             February 1998            1.0

c      Acronym:

c      Purpose: Search for neighbour to edge defined by node knoten

c      Inputs:
c         seg     - Segment number
c         knoten  - Node for search
c         nacb    - Neighbour node
c         node    -  ??
c         ix2(*)  - Facet node connections
c         surpoin - Surface points

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_geom.h'
      include  'c_keyh.h'
      include  'c_pair.h'
      include  'c_tole.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      integer   ix2(dnope2,*),surpoin(*)
      integer   anzkn,anzseg,i,j,knoten,nel2,nachb,node,seg

      save

c     Loop over all Segments i in Patch which are .ne. seg

      fp(1)  = np(191)+surpoin(nsurf2)
      anzkn  = mr(fp(1))
      fp(1)  = np(191)+mr(fp(1)+node+anzkn)
      anzseg = mr(fp(1))

      do i = 1,anzseg

c       if node part of Segment, nachb=i

        if (i.ne.seg) then
          nel2 = mr(fp(1)+i)
          do j = 1,nope2
            if (knoten.eq.ix2(j,nel2)) then
              nachb = i
              return
            endif
          end do ! j
        endif
      end do ! i

      end
