c$Id:$
      subroutine pnumbl(ndm,nr,ns,nt,ntyp, nf,ng, flag)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Correct computation of 11, 14 & 15 node tets     12/09/2007
c       2. Add computation for 64-node brick                06/02/2009
c       3. Force line order based on 'ns' value             16/04/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generage blend blocks
c      Inputs:
c        ndm       - Spatial dimension of mesh
c        nr        - Number of 1-direction increments to generate
c        ns        - Number of 2-direction increments to generate
c        nt        - Number of 3-direction increments to generate
c        ntyp      - Element type to generate
c        flag      - 1-d generation

c      Outputs:
c        nf        - Number last element
c        ng        - Number last node
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   flag
      integer   ndm,nr,ns,nt,ntyp,nf,ng

      save

c     Check the 2-D types

      if(ntyp.lt.10 .or. abs(ntyp).eq.16) then
        if(flag) then
          ng = nr + 1
          if(ns.eq.1) then
            nf = nr
          else
            nf = nr/2
          endif
          nr = nr + 1
        else
          if(abs(ntyp).eq.16) then       ! 16-node quadrilateral
            nf = nr*ns/9
          elseif (ntyp.eq.0) then        !  4-node quadrilateral
            nf = nr*ns
          elseif (abs(ntyp).eq.7) then   ! 6/7-node triangles
            nf = (nr*ns)/2
          elseif (ntyp.ge.8) then        ! 8/9-node quadrilateral
            nf = (nr*ns)/4
          elseif (ntyp.lt.0) then        !  3-node crossed triangles
            nf = 4*nr*ns
          else                           !  3-node triangles
            nf = 2*nr*ns
          endif

c         Determine last node number to be generated

          nr = nr + 1
          ns = ns + 1
          if(ndm.eq.1) ns = 1
          ng = nr*ns
          if(ntyp.eq. -7) then            !  7-node triangles
            ng = ng + (nr-1)*(ns-1)/2
          elseif(ntyp .eq. -1) then       !  3-node crossed triangles
            ng = ng + (nr-1)*(ns-1)
          elseif(ntyp .eq.  8) then       !  8-node quadrilaterals
            ng = ng - ((nr-1)*(ns-1))/4
          endif
        endif

c     3-d generations

      elseif(ntyp.lt.20) then
        if(    ntyp.eq.11) then                 ! 4-node tetrahedron
          nf = (nr*ns*nt)*6
          ng = (nr+1)*(ns+1)*(nt+1)
        elseif(ntyp.eq.12 .or. ntyp.eq.14) then ! 20/27-node hexahedron
          nf = (nr*ns*nt)/8
          ng = (nr+1)*(ns+1)*(nt+1)
        elseif(ntyp.eq.13 .or. ntyp.eq.15) then ! 10/11-node tetrahedron
          nf = 6*(nr*ns*nt)/8
          ng = (nr+1)*(ns+1)*(nt+1)
          if(ntyp.eq.15) then                   ! 11-node tetrahedron
            ng = ng + nf
          endif
        elseif(ntyp.eq.17 .or. ntyp.eq.18) then ! 14/15-node tetrahedron
          nf = 6*(nr*ns*nt)/8
          ng = (nr+1)*(ns+1)*(nt+1)
          if(ntyp.eq.17) then
            ng = ng + 6*(nr*ns*nt)/4 + (nr*ns+ns*nt+nt*nr)/2
          else
            ng = ng + 9*(nr*ns*nt)/4 + (nr*ns+ns*nt+nt*nr)/2
          endif
        elseif(ntyp.eq.19) then                 ! 64-node hexahedron
          nf = (nr*ns*nt)/27
          ng = (nr+1)*(ns+1)*(nt+1)
        else                                    !  8-node hexahedron
          nf = nr*ns*nt
          ng = (nr+1)*(ns+1)*(nt+1)
        endif

c     Shell:

      elseif(ntyp.lt.30) then
        if (ntyp.eq.20) then
          nf = nr*ns
        elseif (ntyp.eq.27) then
          nf = (nr*ns)/2
        elseif (ntyp.ge.28) then
          nf = (nr*ns)/4
        else
          nf = 2*nr*ns
        endif

c       Determine last node number to be generated

        ng = (nr+1)*(ns+1)

c     Line:

      elseif(ntyp.lt.40) then
        if(ns.ge.1) then
          nf = nr/ns
        else
          nf = nr
        endif

        ng = nr + 1

      endif

      end
