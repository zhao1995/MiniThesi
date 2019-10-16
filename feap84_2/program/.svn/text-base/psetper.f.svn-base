c$Id:$
      subroutine psetper(elnk,x,u)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    17/03/2011
c       1. Add prop to set incremental displacements        12/04/2013
c          Remove unused dulnk statements
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Set boundary displacement values from deformation
c              gradient

c     Inputs:
c         elnk(ndm,numnp) - Linked edge indicators
c         x(ndm,numnp)    - Nodal coordinates

c     Outputs:
c         u(ndf,numnp,*)  - Displacements at nodes
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'                 ! numnp,numel
      include   'cdat2.h'                 ! dxlnk(3)
      include   'elpers.h'                ! gradu(3,3)
      include   'prlod.h'                 ! prop
      include   'sdata.h'                 ! ndm,ndf

      integer    elnk(ndm,numnp)
      real*8     x(ndm,numnp), u(ndf,numnp,*)
      real*8     ulnk(3)

      integer    j,k,n, e

      save

c     Loop over linked nodes to set displacements

      do n = 1,numnp
        do j = 1,ndm
          ulnk(j)  = 0.0d0
          u(j,n,4) = 0.0d0
        end do ! j
        do j = 1,ndm
          e = elnk(j,n)
          if(e.gt.0) then        ! Node n is linked to node elnk(j,n)
            do k = 1,ndm
              ulnk(j) = ulnk(j) + gradu(j,k)*(x(k,n) - x(k,e))*prop
            end do ! k
          end if
        end do ! j
        do j = 1,ndm
          e = elnk(j,n)
          if(e.gt.0) then
            u(j,n,4) = (u(j,e,1) - u(j,n,1)) + ulnk(j) ! Add to per nd
          endif
        end do ! j

      end do ! n

      end
