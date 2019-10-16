c$Id:$
      subroutine pout4e(nel,iface,elface,xl,norm,econ,nen1,ndm,n)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plot outline of a face in appropriate colors to
c               indicate surface and edge definitions for 3-d objects

c      Inputs:
c         nel           - Number nodes on face
c         iface(nen1,*) - Face node connections, materials, regions
c         elface(*)     - Pointer array for faces connected to nodes
c         xl(3,*)       - Nodal coordinates for face
c         econ(*)       - Face numbers connected to nodes
c         nen1          - Dimension of iface array
c         ndm           - Dimension of xl array
c         n             - Face number to compute

c      Outputs:
c         none          - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nel,nen1,ndm,ii,jj,i,j,k,n
      integer   iface(nen1,*),elface(*),econ(*)
      real*8    ang, dot, cosa, xe(3,2),xl(ndm,*),norm(3,*)

      save

c     Plot boundaries in appropriate color

      cosa = 0.95d0

      do i = 1,nel

        j = mod(i,nel) + 1

        do ii = 1,ndm
          xe(ii,1) = xl(ii,i)
          xe(ii,2) = xl(ii,j)
        end do ! ii

c       Find edges

        ii = iface(i,n)
        if(ii.gt.0) then

          do jj = elface(ii)+1,elface(ii+1)
            if(econ(jj).gt.0 .and. econ(jj).ne.n) then

              do k = 1,4
                if(iface(k,econ(jj)).eq.iface(j,n)) then

                  ang = abs(dot(norm(1,n),norm(1,econ(jj)),3))
                  if(ang.lt.cosa) then

c                   Draw Line

                    call plotl(xe(1,1),xe(2,1),xe(3,1),3)
                    call plotl(xe(1,2),xe(2,2),xe(3,2),2)
                    go to 200

                  end if

                end if
              end do ! k

            end if
          end do ! jj

        end if

200     continue

      end do ! i

      end
