c$Id:$
      subroutine rignod(ixt,ie,ix,rben,nie,nen,nen1,numnp,numel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set up list of nodes attached to rigid bodies

c      Inputs:
c         ie(nie,*)  - Material set data values
c         ix(nen1,*) - Element nodal connections
c         rben(*)    - List of rigid body numbers for elements
c         nie        - Size of material set data array
c         nen        - Number of nodes/element
c         nen1       - Dimension of ix array
c         numnp      - Number of nodes in mesh
c         numel      - Number of elements in mesh

c      Outputs:
c         ixt(*)     - List of node/rigid body
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,n,nie,nen,nen1,numnp,numel
      integer   ixt(numnp),ie(nie,*),ix(nen1,numel),rben(numel)

      save

c     Zero ixt array

      do n = 1,numnp
        ixt(n) = 0
      end do ! n

c     Loop over elements

      do n = 1,numel

c       Check for material sets that are rigid

        if(ie(2,ix(nen1,n)).gt.0) then
          rben(n) = ie(2,ix(nen1,n))
        endif

c       Check for rigid body element

        if(rben(n).ge.1) then

c         Set nodes on this element to belong to rigid body

          do i = 1,nen
            if(ix(i,n).gt.0) then
              ixt(ix(i,n)) = rben(n)
            endif
          end do ! i

        end if
      end do ! n

      end
