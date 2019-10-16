c$Id:$
      subroutine ctface3(dnope,ics,ix,ip,x,norm,nen,ndm, neps)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    30/03/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor           March 30, 2012            1.0

c      Acronym: Contact Read Facets type 3 (triangular facets)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      logical   facefl
      integer   dnope,neps,nen,ndm, i,j,n, ne,ni, nel
      integer   ics(dnope,*), ix(nen), ip(*), faces(3,5)
      real*8    x(ndm,*), norm(3,*), dnorm(3), dx(3), dy(3), len

      save

      data      faces / 1,3,2, 1,2,4, 2,3,4, 3,1,4, 1,2,3 /

      nel = 0
      do n = 1,nen
        if(ix(n).gt.0) nel = n
      end do ! n
      if(nel.eq.4) then
        ne = 4
        ni = 1
      elseif(nel.eq.3) then
        ne = 5
        ni = 4
      else
        ne = 0
        ni = 1
      end if
      do n = 1,ne,ni
        facefl = .true.
        do j = 1,3
          i = ix(faces(j,n))
          if(i.gt.0) then
            if(ip(i).eq.0) then
              facefl = .false.
            endif
          endif
        end do

c       Check if normal is directed in same way as nodal normals

        if(facefl) then
          do i = 1,3
            dx(i) = x(i,ix(faces(2,n)))-x(i,ix(faces(1,n)))
            dy(i) = x(i,ix(faces(3,n)))-x(i,ix(faces(1,n)))
          end do
          dnorm(1) = dx(2)*dy(3) - dx(3)*dy(2)
          dnorm(2) = dx(3)*dy(1) - dx(1)*dy(3)
          dnorm(3) = dx(1)*dy(2) - dx(2)*dy(1)

          len =  sqrt(dnorm(1)**2 + dnorm(2)**2 + dnorm(3)**2)

          i   = ix(faces(1,n))
          if(len.gt.0.0d0) then
            len = (norm(1,i)*dnorm(1) + norm(2,i)*dnorm(2)
     &           + norm(3,i)*dnorm(3)) / len

c           Store output list of surface nodes in ICS list

            if(len .gt. 0.0d0) then

c             Store output list of surface nodes in ICS list

              neps = neps + 1
              do i = 1,3
                ics(i,neps) = ix(faces(i,n))
              end do ! i
            endif
          endif
        endif

      end do ! n

      end
