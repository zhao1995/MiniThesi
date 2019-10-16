c$Id:$
      subroutine mshckpt(cdamp,cmass,cstif,ip,ndf,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Shift entries in ip array by +1                  25/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Perform check to set active point values.

c      Inputs:
c         cdamp(*)       - Damping   terms
c         cmass(*)       - Mass      terms
c         cstif(*)       - Stiffness terms
c         ndf            - Number dof/node
c         numnp          - Number of nodes in mesh

c      Outputs:
c         ip(ndf+1,*)    - List of active nodes, used for graphics
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  ndf,numnp, i,n, ip(ndf+1,numnp)
      real*8   cdamp(ndf,numnp),cmass(ndf,numnp),cstif(ndf,numnp)

      save

      do n = 1,numnp
        do i = 1,ndf
          if(ip(i+1,n) .eq. 0 ) then
            if(cdamp(i,n).ne.0.0d0 .or.
     &         cmass(i,n).ne.0.0d0 .or.
     &         cstif(i,n).ne.0.0d0) then
              ip(i+1,n) = 1
            endif
          endif
        end do ! i
      end do ! n

      end
