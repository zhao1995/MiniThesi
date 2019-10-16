      subroutine pttrans(ia,t,ul,p,s,nel,ndf,nst,jj)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    11/02/2013
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Set transformation data for element computations using
c               triad arrays

c      Inputs:
c         ia(3)         - Degrees of freedom to transform
c         t(3,3,*)      - 3x3 transformation array
c         ul(ndf,nen,*) - Input values in rotated form
c         p(ndf,*)      - Input of cartesian vector
c         s(nst,*)      - Input of cartesian matrix
c         nel           - Number nodes on element
c         ndf           - Number dofs/node
c         nst           - Matrix dimension
c         jj            - Switch: (1) Rotate ul; (2) s & p: (x) p

c      Outputs:
c         ul(ndf,nen,*) - Output values in cartesian form
c         p(ndf,*)      - Output of rotated vector
c         s(nst,*)      - Output of rotated matrix
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'cdata.h'

      integer    ia(3), nel,ndf,nst, jj
      integer    i,j,k, i1
      real*8     t(3,3,*),ul(ndf,nen,*),p(ndf,*),s(nst,nst), cs(3)

c     Loop over nodes

      i1 = 0
      do i = 1,nel

c       Check for rotation (t(1,1,*) > -100.0d0

        if(t(1,1,i).gt.-100.0d0) then

c         Rotate displacements to cartesian coordinates

          if(jj.eq.1) then
            do k = 1,3
              do j = 1,3
                cs(j) = t(j,1,i)*ul(ia(1),i,k)
     &                + t(j,2,i)*ul(ia(2),i,k)
     &                + t(j,3,i)*ul(ia(3),i,k)
              end do ! j
              do j = 1,3
                ul(ia(j),i,k) = cs(j)
              end do !j
            end do ! k

c         Rotate arrays to local form

          else

            do j = 1,3
              cs(j) = t(1,j,i)*p(ia(1),i)
     &              + t(2,j,i)*p(ia(2),i)
     &              + t(3,j,i)*p(ia(3),i)
            end do ! j
            do j = 1,3
              p(ia(j),i) = cs(j)
            end do ! j

c           Rotate matrix to local form

            if(jj.eq.2) then
              do k = 1,nst
                do j = 1,3
                  cs(j) = s(k,i1+ia(1))*t(1,j,i)
     &                  + s(k,i1+ia(2))*t(2,j,i)
     &                  + s(k,i1+ia(3))*t(3,j,i)
                end do ! j
                do j = 1,3
                  s(k,i1+ia(j)) = cs(j)
                end do ! j
              end do ! k

              do k = 1,nst
                do j = 1,3
                  cs(j) = t(1,j,i)*s(i1+ia(1),k)
     &                  + t(2,j,i)*s(i1+ia(2),k)
     &                  + t(3,j,i)*s(i1+ia(3),k)
                end do ! j
                do j = 1,3
                  s(i1+ia(j),k) = cs(j)
                end do ! j
              end do ! k
            endif ! jj.eq.2
          endif ! jj

        endif ! rotation check

        i1 = i1 + ndf
      end do ! i

      end
