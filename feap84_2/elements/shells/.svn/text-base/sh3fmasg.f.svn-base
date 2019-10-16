c$Id:$
      subroutine sh3fmasg( d, dthk, shp, xjw, s, lint, nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FMASG is the subroutine which computes the
c                        mass matrix for rigid body and eigen analyses.

c        Authors:        R. Taylor

c        Date:           March, 1995.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i , j , k , l , j1 , k1 , lint , nst
      real*8    rhoa  , rhoi  , dthk    , fac  , fac1
      real*8    d(*)  , s(nst,*) , shp(4,4), xjw(4)

c     Compute Mass Array:

      do i = 1,lint

c       Compute Weighted Density/Inertia

        rhoa = d(4)*dthk * xjw(i)
        rhoi = d(8)*rhoa*dthk**2/12.d0

c       Compute Mass

        j1 = 0
        do j = 1 , 4

c         Translational Mass

          fac = rhoa*shp(j,i)
          k1 = 0
          do k = 1 , 4
            fac1 = fac*shp(k,i)
            do l = 1 , 3
              s(j1+l,k1+l) = s(j1+l,k1+l) + fac1
            end do ! l
            k1   = k1 + 6
          end do ! k

c         Rotational Mass

          fac = rhoi*shp(j,i)
          do l = 4 , 6
            s(j1+l,j1+l) = s(j1+l,j1+l) + fac
          end do ! l

          j1 = j1 + 6
        end do ! j

      end do ! i

      end
