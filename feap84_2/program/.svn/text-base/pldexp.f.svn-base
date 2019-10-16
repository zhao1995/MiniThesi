c$Id$
      subroutine pldexp(ldtyp, ldtab,ldnod,ldval, f)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/01/2009
c       1. Increase ldtab to store spin number of displ.    09/03/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Expand loads from tables

c      Inputs:
c        ldtyp          - Load type: 1 = force; 2 = displacement
c        ldtab(4,2,*)   - Pointer, length prop no. table
c        ldnod(*)       - List of nodes
c        ldval(ndf,*)   - List of values

c      Outputs:
c        f(ndf,*)       - Force/displacement values
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'pload1.h'
      include   'sdata.h'

      integer    ldtyp, j,nn, nod, nd,nl
      integer    ldtab(4,2,*),ldnod(*)
      real*8     ldval(ndf,*), f(ndf,*)

c     Set group

      nl = ldtab(2,ldtyp,ldnum)

      call pzero(f, nneq)

c     Do forces (1) or displacements (2)

      if(nl.gt.0) then

c       Loop over nodes in group and store in f

        nd = ldtab(1,ldtyp,ldnum)
        do nn = 1,nl
          nod = ldnod(nd+nn)
          do j = 1,ndf
            f(j,nod) = f(j,nod) + ldval(j,nd+nn)
          end do ! j
        end do ! nn
      endif

      end
