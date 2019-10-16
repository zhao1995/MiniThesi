c$Id:$
      subroutine nodel(ix,ixc,nnac,nnid,nd,ln,ne,numnp,numel,numcels,
     &                 nsum,nen,nen1,ncen,ncen1)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Formation of "nd" and "ln" arrays for profile
c               optimization.

c      Inputs:
c         ix(nen1,*)     - Element connection array
c         ixc(ncen1,*)   - Contact element connection array
c         nnac(*)        - Overlay markers for elements
c         nnid(*)        - Number active dof at nodes
c         numnp          - Number of nodes in mesh
c         numel          - Number of elements in mesh
c         numcels        - Number of contact elements in mesh
c         nen            - Maximum number nodes/element
c         nen1           - Dimension of ix  array
c         ncen           - Maximum number nodes/contact element
c         ncen1          - Dimension of ixc array

c      Outputs:
c         nd(*)          - Nodal-element array
c         ln(*)          - Location array
c         ne(*)          - Element order to minimize front
c         nsum           - Maximum value in location array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   l, m, n,nn,nen,nen1,ncen,ncen1,numnp,numel,numcels,nsum
      integer   nd(numnp),ix(nen1,*),ixc(ncen1,*),ln(numnp),ne(*)
      integer   nnac(numel),nnid(numnp)

      save

c     Formation of "nd" and "ln" arrays

      do m = 1,numnp
        nd(m) = 0
      end do ! m

c     Count elements or constraints attached to nodes

      do m = 1,numel

        if(ix(nen1-1,m).ge.0 .and. nnac(m).eq.0 ) then
          do l = 1,nen
            n = abs(ix(l,m))
            if(n.gt.0) then
              if(nnid(n).gt.0) then
                nd(n) = nd(n) + 1
              endif
            endif
          end do ! l
        end if

      end do ! m

      do m = 1,numcels

        do l = 1,ncen
          n = abs(ixc(l,m))
          if(n.gt.0) then
            if(nnid(n).gt.0) then
              nd(n) = nd(n) + 1
            endif
          endif
        end do ! l

      end do ! m

c     Form location array

      nsum = 0
      do n = 1,numnp
        nsum  = nsum + nd(n)
        ln(n) = nsum
      end do ! n

c     Form nodel-element array

      do m = 1,numel

        if(ix(nen1-1,m).ge.0 .and. nnac(m).eq.0 ) then
          do l = 1,nen
            n  = abs(ix(l,m))
            if(n.gt.0) then
              if(nnid(n).gt.0) then
                nd(n)  = nd(n) - 1
                nn     = ln(n) - nd(n)
                ne(nn) = m
              endif
            endif
          end do ! l
        endif

      end do ! m

      do m = 1,numcels

        do l = 1,ncen
          n  = abs(ixc(l,m))
          if(n.gt.0) then
            if(nnid(n).gt.0) then
              nd(n)  = nd(n) - 1
              nn     = ln(n) - nd(n)
              ne(nn) = m + numel
            endif
          endif
        end do ! l

      end do ! m

      end
