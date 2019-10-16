c$Id:$
      subroutine rmas3d(xl,s,xmasc,jmasc,
     &                  xcg, tmasc, ndf,ndm,nel, nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 2-d quadrilateral and 3-d brick mass computations

c      Inputs:
c         xl(ndm,*)  - Element nodal coordinate
c         s(nst,*)   - Element mass array
c         xcg(3)     - Location of reference point for calculation
c         ndf        - Number dof/node
c         ndm        - Spatial dimension of mesh
c         nel        - Number of nodes on element
c         nst        - Dimension of s array

c      Outputs:
c         tmasc      - Mass of element
c         xmasc(3)   - Center of mass for element
c         jmasc(3,3) - Element inertia
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ndf,ndm,nel,nst, i,j,k,l, i1,j1,k1,k2
      real*8    a0, a1,xl(ndm,nel),s(nst,nst)
      real*8    tmasc, xmasc(3), jmasc(3,3), xcg(3)

      save

c     Consistent mass form

      tmasc = 0.0d0
      do i = 1,3
        do j = 1,3
          jmasc(j,i) = 0.0d0
        end do ! j
        xmasc(i) = 0.0d0
      end do ! i

      i1 = 1
      do i = 1,nel ! {
        a0 =  0.0d0
        j1 = 1
        do j = 1,nel ! {
          a0 = a0 + s(i1,j1)
          j1 = j1 + ndf
        end do ! j   }

c       Accumulate mass

        tmasc = tmasc + a0

c       Accumulate for center of mass

        do k = 1,ndm ! {
          xmasc(k) = xmasc(k) + (xl(k,i) - xcg(k))*a0
        end do ! k   }

c       Accumulate inertia tensor

        j1 = 1
        do j = 1,nel ! {
          if(ndm.eq.2) then
            a0 = ( xl(1,i)-xcg(1) )*( xl(1,j)-xcg(1) )*s(i1+1,j1+1)
     &         + ( xl(2,i)-xcg(2) )*( xl(2,j)-xcg(2) )*s(i1  ,j1  )
     &         - ( xl(2,i)-xcg(2) )*( xl(1,j)-xcg(1) )*s(i1  ,j1+1)
     &         - ( xl(1,i)-xcg(1) )*( xl(2,j)-xcg(2) )*s(i1+1,j1  )
            if(ndf.ge.3) then
              a0 = a0 + ( xl(1,i)-xcg(1) )*s(i1+1,j1+2)
     &                + ( xl(1,j)-xcg(1) )*s(i1+2,j1+1)
     &                - ( xl(2,i)-xcg(2) )*s(i1  ,j1+2)
     &                - ( xl(2,j)-xcg(2) )*s(i1+2,j1  ) + s(i1+2,j1+2)
            endif
            jmasc(3,3) = jmasc(3,3) + a0
          elseif(ndm.eq.3) then
            a0 = (( xl(1,i)-xcg(1) )*( xl(1,j)-xcg(1) )
     &         +  ( xl(2,i)-xcg(2) )*( xl(2,j)-xcg(2) )
     &         +  ( xl(3,i)-xcg(3) )*( xl(3,j)-xcg(3) ))*s(i1,j1)
            do k = 1,3 ! {
              jmasc(k,k) = jmasc(k,k) + a0
              a1 = (xl(k,i) - xcg(k))*s(i1,j1)
              do l = 1,3 ! {
                jmasc(k,l) = jmasc(k,l) - a1*(xl(l,j) - xcg(l))
              end do ! l   }
            end do ! k   }
          endif
          j1 = j1 + ndf
        end do ! j   }

        i1 = i1 + ndf
      end do ! i   }

c     Check for inertia part of mass matrix

c     if((ndm.eq.3 .and. ndf .ge. ndm+3)  .or.
c    &   (ndm.eq.2 .and. ndf .ge. ndm+1)) then

      if (ndm.eq.3 .and. ndf .ge. ndm+3)  then

c       Set maximum number of terms to add
        if(ndm.eq.2) then
          k1 = 3
          k2 = 0
        elseif(ndm.eq.3) then
          k1 = 1
          k2 = 3
        else
          k1 = 4
          k2 = 0
        endif

        i1 = k2
        do i = 1,nel
          j1 = k2
          do j = 1,nel
            do k = k1,3
              do l = k1,3
                jmasc(k,l) = jmasc(k,l) + s(i1+k,j1+l)
              end do ! l
            end do ! k
            j1 = j1 + ndf
          end do ! j
          i1 = i1 + ndf
        end do ! i

      endif

      end
