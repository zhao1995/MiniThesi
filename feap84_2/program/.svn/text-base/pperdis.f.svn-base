c$Id:$
      subroutine pperdis(id,x,u)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    17/03/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set boundary displacement of periodic boundary problem

c      Inputs:
c         gradu(3,3)      - Displacement gradient (in common 'elpers')
c         id(ndf,numnp,2) - Equation/boundary codes
c         x(ndm,numnp)    - Reference coordinates

c     Outputs:
c         u(ndf,numnp)  - Force/Displacements at nodes
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'                 ! numnp,numel
      include   'elpers.h'                ! gradu(3,3)
      include   'prlod.h'                 ! prop
      include   'sdata.h'                 ! ndm,ndf

      integer    id(ndf,numnp,2)
      real*8     x(ndm,numnp), u(ndf,numnp)

      integer    i,j,n

      save

c     Loop over nodes and directions

      do n = 1,numnp
        do j = 1,ndm
          if(id(j,n,2).ne.0) then        ! Restrained dof
            u(j,n) = 0.0d0
            do i = 1,ndm
              u(j,n) = u(j,n) + gradu(j,i)*x(i,n)*prop
            end do ! i
          end if

        end do ! j
      end do ! n

      end
