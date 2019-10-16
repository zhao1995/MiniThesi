      subroutine pidset(ip,ie,iedof,id,nty,ix,nie,nen,nen1,ndf,
     &                  numnp,numel,nummat)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    21/09/2007
c       1. Add call to 'uidset' to activate equations       20/12/2012
c       2. Separate ip array into parts for setting b.c.    25/01/2013
c          and parts to set node type for plots
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set final boundary restraints after ties/cbou, etc.

c      Inputs:
c         ie(nie,*)      - Material set assembly information
c         iedof(ndf,*,*) - Material set nodal assembly information
c         id(ndf,*)      - Boundary condition and equation number array
c         nty(*)         - Nodal type
c         ix(nen1,*)     - Element nodal connection lists
c         nie            - Dimension of ie array
c         nen            - Maximum number of nodes/element
c         nen1           - Dimension for ix array
c         ndf            - Number dof/node
c         numnp          - Number of nodes in mesh
c         numel          - Number of elemenst in mesh
c         nummat         - Number of material sets

c      Outputs:
c         ip(ndf+1,*)    - List of boundary restraints and unused dof's
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'prflag.h'
      include   'pointer.h'
      include   'comblk.h'

      integer    nie,nen,nen1,ndf,numnp,numel,nummat
      integer    n,i,j, ii,mg, ma
      integer    ip(ndf+1,*),  ie(nie,*), iedof(ndf,nen,*), id(ndf,*)
      integer    ix(nen1,*), nty(*)

c     Remove unused dof's using ie(nie,*) array

c     do n = 1,numnp
c       do j = 1,ndf
c         ip(j+1,n) = 0
c       end do ! j
c     end do ! n

c     Check nodes on each element for active dof's

      do n = 1,numel
        mg = ix(nen1,n)

c       Loop over material sets

        do ma = 1,nummat
          if(ie(nie-2,ma).eq.mg) then
            do i = 1,nen
              ii = ix(i,n)
              if(ii.gt.0) then
                do j = 1,ndf
                  if(iedof(j,i,ma).gt.0) then
                    ip(iedof(j,i,ma)+1,ii) = 1
                  endif
                end do ! j
              endif
            end do ! i
          endif
        end do ! ma
      end do ! n

c     Check for point mass, dampers, stiffness

      if(nmfl) then
        call mshckpt(hr(np(86)),hr(np(87)),hr(np(88)),ip,ndf,numnp)
      endif

c     Set b.c. restraints for unused dof's

      do n = 1,numnp
        do j = 1,ndf
          if(ip(j+1,n).eq.0) then
            id(j,n) = -1000
          else
            id(j,n) = abs(id(j,n))
          end if
        end do ! j
      end do ! n

c     Remove unused nodes - for graphics

      do n = 1,numel
        do i = 1,nen
          ii = ix(i,n)
          if(ii.gt.0) ip(1,ii) = 1
        end do ! i
      end do ! n

c     Set flat to indicate node is not used

      do n = 1,numnp
        if(ip(1,n) .eq. 0) then
          nty(n) = -1
        end if
      end do ! n

c     Fix all unspecified coordinate dof's

      do n = 1,numnp
        if(nty(n).lt.0) then
          do i = 1,ndf
            id(i,n) = -1001
          end do ! i
        endif
      end do ! n

c     User modification to active equations

      call uidset(id,ndf,numnp)

      end
