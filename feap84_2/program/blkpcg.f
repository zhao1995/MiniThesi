c$Id:$
      subroutine blkpcg(id,jp,ndf,numnp,neq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:       Generate profile for nodal blocks for use in
c                     PCG preconditioning.

c      Inputs:
c         id(ndf,*) - Equation numbers for each dof
c         ndf       - Number dof/node
c         numnp     - number of nodes in mesh
c         neq       - Number active equations

c      Outputs:
c         jp(neq)   - Pointers for nodal blocks
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'iofile.h'

      integer  i,i0,i1,i2, n, ndf,numnp,neq
      integer  id(ndf,*), jp(*)

      save

c     Profile for diagonal nodal blocks

      do n = 1,neq
        jp(n) = 0
      end do ! n

      do n = 1,numnp
        i0 = neq
        i1 = 0
        do i = 1,ndf
          if(id(i,n).gt.0) then
            i0 = min(i0,id(i,n))
            i1 = max(i1,id(i,n))
          end if
        end do ! i

c       Compute column heights

        i2 = 0
        do i = i0+1,i1
          i2 = i2 + 1
          jp(i) = max(jp(i),i2)
        end do ! i
      end do ! n

c     Convert to profile pointers

      do i = 2,neq
        jp(i) = jp(i) + jp(i-1)
      end do ! i

      if(ior.lt.0) then
        write(*,*) ' STORAGE FOR AUR =',jp(neq)
      endif

      end
