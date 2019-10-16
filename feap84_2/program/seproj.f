c$Id:$
      subroutine seproj(p,s,se,dt,st,ser,ix,nel,nen,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Assemble element projection quantities into global ones

c      Inputs:
c        p(*)    - Element nodal weights
c        s(*)    - Element projected quantities
c        se(*)   - Element error quantities
c        ix*)    - Element node numbers
c        nel     - Number element nodes
c        nen     - Number element nodes max
c        numnp   - Number global nodes

c      Outputs:
c        dt(*)   - Global nodal weights
c        st(*)   - Global projected quanties
c        ser(*)  - Global error quanties
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'strnum.h'
      include   'iofile.h'

      integer    nel,nen,numnp, nd,nn,i, ix(*)
      real*8     p(nen),s(nen,*),se(*),dt(numnp),st(numnp,*),ser(*)

      save

c     Lumped projection assembly over element nodes

      do nd = 1,nel
        nn = ix(nd)
        if(nn.gt.0) then
          if(istc.eq.1) then
            dt(nn)  = dt(nn)  + p(nd)
          endif
          ser(nn) = ser(nn) + se(nd)

c         Assemble integated stress

          do i = 1,iste
            st(nn,i) = st(nn,i) + s(nd,i)
          end do ! i
        endif
      end do ! nd

      end
