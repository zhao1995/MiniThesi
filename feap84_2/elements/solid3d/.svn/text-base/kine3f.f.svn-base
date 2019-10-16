c$Id:$
      subroutine kine3f (shps,shpi,ul,ui,uin,f,f0,fi,df,detfi,
     &                   ndf,ndfi,nel,nen,fact,finc)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 3-D finite deformation kinematics

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      logical  finc
      integer  ndf,ndfi,nel,nen, i,j,k
      real*8   fact,factn,deti
      real*8   shps(4,*),shpi(3,*),ul(ndf,nen,*),ui(ndfi,*),uin(ndfi,*)
      real*8   f(3,3,*),fi(3,3),df(3,3),f0(3,3,*),detfi(*)

c     Compute compatible deformation gradient at t-n+1
c        F = I + GRAD u

      do i = 1,3
        do j = 1,3
          f(i,j,1)  = 0.0d0
          f(i,j,2)  = 0.0d0
          df(i,j) = 0.0d0
          do k = 1,nel
            f(i,j,1) = f(i,j,1) + ul(i,k,1)*shps(j,k)
            df(i,j)  = df(i,j)  + ul(i,k,2)*shps(j,k)
          end do ! k
          f(i,j,2) = f(i,j,1) - df(i,j)
        end do ! j
        f(i,i,1) = f(i,i,1) + 1.0d0
        f(i,i,2) = f(i,i,2) + 1.0d0
      end do ! i

c     Add enhanced deformation gradient if required.

      fact  = 0.0d0
      if(finc) then
        do j = 1,3
          do i = 1,3
            do k = 1,3
              f(i,j,1) = f(i,j,1) +  ui(i,k)*shpi(j,k)
              f(i,j,2) = f(i,j,2) + uin(i,k)*shpi(j,k)
            end do ! k
          end do ! i
        end do ! j

c       Compute fourth enhanced mode factor for def grad and shp fn.

        factn = 0.0d0
        do i = 1,3
          fact  = fact  +  ui(i,4)*shpi(i,4)
          factn = factn + uin(i,4)*shpi(i,4)
        end do ! i

        do i = 1,3
          do j = 1,3
            f(i,j,1) = f(i,j,1) + fact  * f0(i,j,1)
            f(i,j,2) = f(i,j,2) + factn * f0(i,j,2)
          end do ! j
        end do ! i

      endif

c     Invert F_n

      detfi(2) = f(1,1,2)*f(2,2,2)*f(3,3,2) + f(1,2,2)*f(2,3,2)*f(3,1,2)
     &         + f(1,3,2)*f(2,1,2)*f(3,2,2) - f(3,1,2)*f(2,2,2)*f(1,3,2)
     &         - f(3,2,2)*f(2,3,2)*f(1,1,2) - f(3,3,2)*f(2,1,2)*f(1,2,2)

c     Invert F_n+1

      detfi(1) = f(1,1,1)*f(2,2,1)*f(3,3,1) + f(1,2,1)*f(2,3,1)*f(3,1,1)
     &         + f(1,3,1)*f(2,1,1)*f(3,2,1) - f(3,1,1)*f(2,2,1)*f(1,3,1)
     &         - f(3,2,1)*f(2,3,1)*f(1,1,1) - f(3,3,1)*f(2,1,1)*f(1,2,1)

      deti    = 1.d0/detfi(1)
      fi(1,1) = (f(2,2,1)*f(3,3,1) - f(3,2,1)*f(2,3,1))*deti
      fi(1,2) =-(f(1,2,1)*f(3,3,1) - f(3,2,1)*f(1,3,1))*deti
      fi(1,3) = (f(1,2,1)*f(2,3,1) - f(2,2,1)*f(1,3,1))*deti
      fi(2,1) =-(f(2,1,1)*f(3,3,1) - f(3,1,1)*f(2,3,1))*deti
      fi(2,2) = (f(1,1,1)*f(3,3,1) - f(3,1,1)*f(1,3,1))*deti
      fi(2,3) =-(f(1,1,1)*f(2,3,1) - f(2,1,1)*f(1,3,1))*deti
      fi(3,1) = (f(2,1,1)*f(3,2,1) - f(3,1,1)*f(2,2,1))*deti
      fi(3,2) =-(f(1,1,1)*f(3,2,1) - f(3,1,1)*f(1,2,1))*deti
      fi(3,3) = (f(1,1,1)*f(2,2,1) - f(2,1,1)*f(1,2,1))*deti

      end
