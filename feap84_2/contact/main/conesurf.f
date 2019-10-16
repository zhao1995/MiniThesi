c$Id:$
      subroutine conesurf(ix1,knotn,lockn,laenge,dnope,neps,nope)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Anna Haraldsson             February 1998            1.0

c      Acronym: Contact ONE SURFace

c      Purpose: Get patch of elements for each node of "ix1" array
c               and store in "knotn".
c               Length of array knotn = laenge

c      Inputs:
c         ix1(*)  - Facet node connections
c         lockn   - Local counting of surface
c         laenge  - Length of total knotn
c         dnope   - Dimension facet
c         neps    - Number of facets
c         nope    - Nodes per facet

c      Outputs:
c         knotn(*)- Node for search
c         knotn(lockn)                       = # surface nodes = anzkn
c         knotn(lockn+node)                  = global number of node
c         knotn(lockn+node+anzkn)            = pointer to patch
c         knotn(knotn(lockn+node+anzkn))     = # segments in patch
c         knotn(knotn(lockn+node+anzkn)+seg) = global number of segment
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_contac.h'
      include  'c_comnd.h'
      include  'c_dict.h'
      include  'iofile.h'

      logical   drin
      integer   anzkn,a,b,i,seg,dnope,neps,nope,node,lnod
      integer   lockn,laenge,einlaeng, knotn(*),ix1(dnope,neps)

      save

      anzkn    = 0
      einlaeng = 0

c     Loop over all nodes on "ix1" surface

      do i = 1,neps
        do a = 1,nope

c         Check for repeated nodes in list

          drin = .false.
          do b = 1,anzkn
            if (knotn(lockn+b).eq.ix1(a,i)) then
              drin = .true.
            endif
          end do ! b

c         Add new node to list

          if (.not.drin) then
            anzkn              = anzkn + 1
            knotn(lockn+anzkn) = ix1(a,i)
            einlaeng           = einlaeng + 1
          endif
        end do ! a
      end do ! i

      knotn(lockn)         = anzkn
      knotn(lockn+anzkn+1) = 2*anzkn + lockn

c     For each node set list of attached facets

      do node = 1,anzkn

        seg  = 0
        lnod = lockn + node

c       Look for a facet with number "node"

        do i = 1,neps
          do a = 1,nope
            if (knotn(lnod).eq.ix1(a,i)) then
              seg = seg + 1
              knotn(knotn(lnod+anzkn)+1+seg) = i
            endif
          end do ! a
        end do ! i

        knotn(knotn(lnod+anzkn)+1) = seg

        if (node.lt.anzkn) then
          knotn(lnod+anzkn+1) = knotn(lnod+anzkn) + seg + 1
        endif

      end do ! node

      lnod  = lockn + anzkn
      lockn = knotn(lockn+anzkn+anzkn) + seg + 1
      do node = 1,anzkn
        knotn(node+lnod) = knotn(node+lnod) + laenge
      end do ! node

      end
