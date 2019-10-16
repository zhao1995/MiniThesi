c$Id:$
      subroutine csurface1(cs0,ics,surpoin)

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
c         knotn(*)- Facets connected to nodes: Assigned pointer np(195)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_contac.h'
      include   'c_comnd.h'
      include   'c_dict.h'
      include   'cdata.h'
      include   'pointer.h'
      include   'comblk.h'

      logical    setvar,palloc
      integer    surf,ofsurf,dnope,neps,nope, i, nod, nfct
      integer    surnod,surptr,surfac
      integer    ics(*),surpoin(3,*)
      real*8     cs0(nr0,n0c1:nc01,*)

      save

c     Allocate temporary array to store data for surface

      setvar = palloc(136,'CTEM1',numnp,1)

c     Loop over individual contact surfaces

      surnod = 0
      surptr = 0
      surfac = 0
      do surf = 1,numcs

c       Extract data from surface table "CS0"

        ofsurf        = nint(cs0(2,-1,surf))
        nope          = nint(cs0(2, 0,surf))
        dnope         = nint(cs0(4,-1,surf))
        neps          = nint(cs0(3,-1,surf))
        surpoin(1,surf) = surnod
        surpoin(2,surf) = surptr
        surpoin(3,surf) = surfac

c       Check facets on each individual surface

c       nod  = 1
        call cont_con(ics(ofsurf),mr(np(136)),dnope,nope,neps,
     &                numnp, nod, nfct)

c       Allocate space and store connection data

        setvar = palloc( 195,'PNSEG',surnod+nod, 1)
        do i = 0,nod-1
          mr(np(195)+surnod+i) = mr(np(137)+i)
        end do ! i
        setvar = palloc( 193,'INSEG',surptr+nod+1, 1)
        do i = 0,nod
          mr(np(193)+surptr+i) = mr(np(138)+i)
        end do ! i
        setvar = palloc( 194,'CNSEG',surfac+nfct, 1)
        do i = 0,nfct-1
          mr(np(194)+surfac+i) = mr(np(139)+i)
        end do ! i

c       Destroy temporary allocation

        setvar = palloc( 137,'CTEM2',0,1)
        setvar = palloc( 138,'CTEM3',0,1)
        setvar = palloc( 139,'CTEM4',0,1)

c       Save total length of array

        surnod = surnod + nod
        surptr = surptr + nod + 1
        surfac = surfac + nfct

      end do ! surf

c     Put final points

      surpoin(1,numcs+1) = surnod
      surpoin(2,numcs+1) = surptr
      surpoin(3,numcs+1) = surfac

c     Destroy temporary array

      setvar = palloc( 136,'CTEM1',   0, 1)

      end
