c$Id:$
      subroutine aslid3da(ib,ip,ix,norm,ma,nen,nen1,numnp,numel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Auto surface descriptions in 3-d
c               Compute boundary connection array

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_contac.h'

      integer    ma,nen,nen1,numnp,numel
      integer    i, j,jj, n, nfaces
      real*8     norm13, norm24

      integer    ib(numnp),ip(numnp),ix(nen1,numel)
      integer    faces(4,6),ic(4)
      real*8     norm(3,numnp)

      save

      data       faces / 1,4,3,2, 1,2,6,5, 2,3,7,6,
     &                   3,4,8,7, 4,1,5,8, 5,6,7,8 /

      call cdebug0 ('    aslid3da',-1)

c     Loop over patches to compute boundary connection array

      do n = 1,numnp
        ib(n) = max(0,ib(n))
        ip(n) = 0
      end do ! n

      do n = 1,numel
        if(ma.eq.0 .or. ma.eq.ix(nen1,n)) then
          jj = 0
          do j = 1,min(8,nen)
            if(ix(j,n).gt.0) jj = jj + 1
          end do ! j

c         8-node brick elements

          if(jj.eq.8) then
            nfaces = 6
            do i = 1,nfaces
              jj = 0
              do j = 1,4
                ic(j) = ix(faces(j,i),n)
                jj    = jj + ib(ic(j))
              end do ! j
              if(jj.eq.4) then

c               Check normal to see if facet can be an interior surface!

                norm13 = norm(1,ic(1))*norm(1,ic(3))
     &                 + norm(2,ic(1))*norm(2,ic(3))
     &                 + norm(3,ic(1))*norm(3,ic(3))

                norm24 = norm(1,ic(2))*norm(1,ic(4))
     &                 + norm(2,ic(2))*norm(2,ic(4))
     &                 + norm(3,ic(2))*norm(3,ic(4))

                if(norm13 .gt. -0.25d0 .and. norm24 .gt. -0.25d0 ) then
                  do j = 1,4
                    ip(ic(j)) = ip(ic(j)) + 1
                  end do ! j
                endif
              endif
            end do ! i
          endif ! nen = 8
        endif
      end do ! n

      if(ifdb) then
        call iprint(ip,1,numnp,1,'Facets/Node')
      endif

c     Convert to pointer

      do n = 2,numnp
        ip(n) = ip(n) + ip(n-1)
      end do ! n

      if(ifdb) then
        call iprint(ip,1,numnp,1,'Facets/Node Pointer')
      endif

      end
