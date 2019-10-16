c$Id:$
      subroutine heproj(s,st,ix,nel,nen,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    04/01/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Assemble history projection quantities into global ones

c      Inputs:
c        s(*)    - User element projected quantities
c        ix*)    - Element node numbers
c        nel     - Number element nodes
c        nen     - Number element nodes max
c        numnp   - Number global nodes

c      Outputs:
c        st(*)   - Global projected quanties
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'eldatp.h'

      integer    nel,nen,numnp, nd,nn,i, ix(*)
      real*8     s(nen,*),st(numnp,*)

      save

c     Lumped projection assembly over element nodes

      do nd = 1,nel
        nn = ix(nd)
        if(nn.gt.0) then

c         Assemble integated stress

          do i = 1,plhmax
            st(nn,i) = st(nn,i) + s(nd,i)
          end do ! i
        endif
      end do ! nd

      end
