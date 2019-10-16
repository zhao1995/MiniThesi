c$Id:$
      subroutine pldir(ntyp,vec,dir,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    18/06/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Inputs:
c       ntyp(*)      - Active node flag
c       dir(
c       numnp        - Number of nodes

c     Outputs:
c       vec(numnp,9) - Vectors for director plots
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   numnp
      integer   ntyp(numnp)
      real*8    dir(3,3,6,numnp), vec(numnp,*)

      integer   n, i

c     Loop over nodes

      do n = 1,numnp

c       Check for active node

        if(ntyp(n).ge.0) then

c         Store undeformed vector

          do i = 1,3
            vec(n,i  ) = dir(i,1,6,n)
            vec(n,i+3) = dir(i,2,6,n)
            vec(n,i+6) = dir(i,3,6,n)
            vec(n,i+9) = 0.0d0
          end do ! i
        endif
      end do ! n

      end
