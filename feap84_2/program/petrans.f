      subroutine petrans(ia,eang,ul,p,s,nel,ndf,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Set transformation data for element computations using
c               Euler angles

c      Inputs:
c         ia(3)     - Degrees of freedom to transform
c         eang(3,*) - Angle data for element nodes
c         ul(ndf,nen,*) - Input values in rotated form
c         p(ndf,*)      - Input of cartesian vector
c         s(nst,*)      - Input of cartesian matrix
c         nel           - Number nodes on element
c         ndf           - Number dofs/node
c         nst           - Matrix dimension
c         isw           - Switch: (1) Rotate ul; (2) Rotate p and s

c      Outputs:
c         ul(ndf,nen,*) - Output values in cartesian form
c         p(ndf,*)      - Output of rotated vector
c         s(nst,*)      - Output of rotated matrix
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'cdata.h'

      logical    etrans
      integer    ia(3), nel,ndf,nst, isw
      integer    i,j,k, i1
      real*8     eang(3,*),ul(ndf,nen,*),p(ndf,*),s(nst,nst)
      real*8     sn(3),cs(3), t(3,3)

c     Loop over nodes

      i1 = 0
      do i = 1,nel
        etrans = .false.
        do j = 1,3
          if(eang(j,i).ne.0.0d0) then
            call pdegree(eang(j,i), sn(j),cs(j))
            etrans = .true.
          else
            sn(j) = 0.0d0
            cs(j) = 1.0d0
          endif
        end do ! j

c       Perform transformations for non-zero rotations

        if(etrans) then
          t(1,1) =  cs(1)*cs(3) - sn(1)*sn(2)*sn(3)
          t(1,2) = -sn(1)*cs(2)
          t(1,3) =  cs(1)*sn(3) + sn(1)*sn(2)*cs(3)
          t(2,1) =  sn(1)*cs(3) + cs(1)*sn(2)*sn(3)
          t(2,2) =  cs(1)*cs(2)
          t(2,3) =  sn(1)*sn(3) - cs(1)*sn(2)*cs(3)
          t(3,1) = -cs(2)*sn(3)
          t(3,2) =  sn(2)
          t(3,3) =  cs(2)*cs(3)

c         Rotate displacements to cartesian coordinates

          if(isw.eq.1) then
            do k = 1,6
              do j = 1,3
                cs(j) = t(j,1)*ul(ia(1),i,k)
     &                + t(j,2)*ul(ia(2),i,k)
     &                + t(j,3)*ul(ia(3),i,k)
              end do ! j
              do j = 1,3
                ul(ia(j),i,k) = cs(j)
              end do !j
            end do ! k

c         Rotate arrays to local form

          else

            do j = 1,3
              cs(j) = t(1,j)*p(ia(1),i)
     &              + t(2,j)*p(ia(2),i)
     &              + t(3,j)*p(ia(3),i)
            end do ! j
            do j = 1,3
              p(ia(j),i) = cs(j)
            end do ! j

c           Rotate matrix to local form

            if(isw.eq.2) then
              do k = 1,nst
                do j = 1,3
                  cs(j) = s(k,i1+ia(1))*t(1,j)
     &                  + s(k,i1+ia(2))*t(2,j)
     &                  + s(k,i1+ia(3))*t(3,j)
                end do ! j
                do j = 1,3
                  s(k,i1+ia(j)) = cs(j)
                end do ! j
              end do ! k

              do k = 1,nst
                do j = 1,3
                  cs(j) = t(1,j)*s(i1+ia(1),k)
     &                  + t(2,j)*s(i1+ia(2),k)
     &                  + t(3,j)*s(i1+ia(3),k)
                end do ! j
                do j = 1,3
                  s(i1+ia(j),k) = cs(j)
                end do ! j
              end do ! k
            endif ! isw.eq.2
          endif ! isw
        endif ! etrans
        i1 = i1 + ndf
      end do ! i

      end
