c$Id:$
      subroutine p3econ(iface,elface,econ,nface,ifc1,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Determine the exterior faces connected to each node

c      Inputs:
c         iface(ifc1,*)  - Face node connections, materials, etc.
c         elface(*)      - Pointer to number of faces for each node
c         nface          - Number of element exterior faces
c         ifc1           - Dimension of iface array
c         numnp          - Number of nodes in mesh

c      Outputs:
c         econ(*)        - faces connected to each node
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nface,ifc1,numnp,ii,jj,i,n
      integer   iface(ifc1,nface),elface(*),econ(*)

      save

      do n = 1,elface(numnp+1)
        econ(n) = 0
      end do ! n

      do n = 1,nface
        if(iface(ifc1,n).gt.0) then
          if(iface(1,n).eq.iface(4,n) .and.
     &       iface(2,n).eq.iface(3,n)) then
          else
            do i = 1,4
              ii = iface(i,n)
              if(ii.gt.0) then

                do jj = elface(ii)+1,elface(ii+1)
                  if(econ(jj).eq.0) then
                    econ(jj) = n
                    go to 100
                  end if
                end do ! jj

              end if
100           continue
            end do ! i
          endif
        end if
      end do ! n

      end
