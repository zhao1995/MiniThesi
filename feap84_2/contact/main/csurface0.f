c$Id:$
      subroutine csurface0(cs0,ics,surpoin)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Name changed from csurfaces                      04/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Anna Haraldsson             February 1998            1.0

c      Acronym:

c      Purpose: Get patch of elements for each node of "ics" array
c               and store in "knotn".

c      Inputs:
c         cs0(*)  - Surface parameters
c         ics(*)  - Surface facets
c         surpoin - Surface points

c      Outputs:
c         knotn(*)- Facets connected to nodes: Assigned pointer np(191)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_contac.h'
      include   'c_comnd.h'
      include   'c_dict.h'
      include   'pointer.h'
      include   'comblk.h'

      logical    setvar,palloc
      integer    surf,ofsurf,dnope,neps,nope, i, lockn,laenge
      integer    ics(*),surpoin(*)
      real*8     cs0(nr0,n0c1:nc01,*)

      save

c     Loop over individual contact surfaces

      laenge = 0
      do surf = 1,numcs

c       Extract data from surface table "CS0"

        ofsurf        = nint(cs0(2,-1,surf))
        nope          = nint(cs0(2, 0,surf))
        dnope         = nint(cs0(4,-1,surf))
        neps          = nint(cs0(3,-1,surf))
        surpoin(surf) = laenge

c       Allocate temporary array to store data for one surface

        setvar = palloc(136,'CTEM1',5*neps*dnope,1)

c       Check facets on each individual surface

        lockn  = 1
        call conesurf(ics(ofsurf),mr(np(136)),lockn,laenge,
     &                dnope,neps,nope)

c       Allocate space and store connection data

        setvar = palloc( 191,'KNOTN',laenge+lockn, 1)
        do i = 0,lockn-1
          mr(np(191)+laenge+i) = mr(np(136)+i)
        end do ! i

c       Save total length of array

        laenge = laenge + lockn

c       Destroy temporary array

        setvar = palloc( 136,'CTEM1',   0, 1)

      end do ! surf

      end
