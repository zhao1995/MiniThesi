c$Id$
      subroutine pldout(ldtab,ldnod,ldval, iq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    05/04/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: output forces and displacements for load tables

c      Inputs:
c        ldtab(4,2,*)   - Pointer, length prop no. table
c        ldnod(*)       - List of nodes
c        ldval(ndf,*)   - List of values
c        iq(*)          - Node list

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iodata.h'
      include   'pload1.h'
      include   'prld1.h'
      include   'sdata.h'

      include   'pointer.h'
      include   'comblk.h'

      integer    i,j,n,nn, nod, nd,nl,npl
      integer    ldtab(4,2,*),ldnod(*), iq(*)
      real*8     ldval(ndf,*)

c     Loop over groups

      integer    iz
      data       iz / 0 /

      do n = 1,ldtot

c       Do forces (1) and displacements (2)

        do i = 1,2
          nl = ldtab(2,i,n)
          if(nl.gt.0) then

c           Set proportional load value for group

            npl = ldtab(3,i,n)
            if(npl.le.0) then
              write(ios,2001) 'LOAD'
            else
              write(ios,2001) 'LOAD PROP',npl
            endif

            if(i.eq.1) then
              write(ios,2001) 'FORCE conditions'
            else
              write(ios,2001) 'DISPlacement conditions'
            endif

c           Loop over nodes in group

            nd = ldtab(1,i,n)
            do nn = 1,nl
              nod = iq(ldnod(nd+nn)+1)
              write(ios,2002) nod,iz,(ldval(j,nd+nn),j=1,ndf)
            end do ! nn

            endif

        end do ! i
        write(ios,2001) 'LOAD END'
      end do ! n

c     formats

2001  format(/a:,i10)
2002  format(i9,i2,1p,14e14.6)

      end
