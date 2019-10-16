c$Id:$
      subroutine peigout(subtype,d,dp,dpp,dtl,tol,nv,itt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    03/09/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Eigevalue data ouputs

c      Inputs:
c         subtype   - Solver name
c         d(*)      - Eigenvalues
c         dp(*)     - Square root of values (rad/s & Hz.)
c         dpp(*)    - Period
c         dtl(*)    - Residuals
c         tol       - Solution tolerance
c         nv        - Number of values
c         itt       - Iterations to convergence

c      Outputs:
c         To files
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'evdata.h'
      include   'iofile.h'
      include   'pconstant.h'

      integer    itt, nc,nv, i,j,n
      character  subtype*(*), subname*10
      real*8     dm,dr,tol
      real*8     d(*),dp(*),dpp(*),dtl(*)

      save

c     Output solution values

      subname = subtype
      nc      = index(subname,' ') - 1

      write(iow,2000) subname(1:nc),itt,(d(n),n=1,nv)
      if(itt.gt.1) write(iow,2001) subname(1:nc),itt,(dtl(n),n=1,nv)
      if(imtyp.eq.1) then
        write(iow,2002) subname(1:nc),(dp(n),n=1,nv)
        dm = 0.5d0/pi
        write(iow,2003) subname(1:nc),(dp(n)*dm,n=1,nv)
        write(iow,2004) subname(1:nc)
        dr = 1.d0/dm
        do i = 0,nv-1,4
          do j = 1,min(4,nv-i)
            if(abs(dp(i+j)).gt.tol*10.d0) then
              dpp(j) = dr/dp(i+j)
            else
              dpp(j) = 0.0d0
            endif
          end do ! j
          write(iow,2005) (dpp(j),j=1,min(4,nv-i))
        end do ! i
      endif
      if(ior.lt.0) then
        write(*,2000) subname(1:nc),itt,(d(n),n=1,nv)
        if(itt.gt.1) write(*,2001) subname(1:nc),itt,(dtl(n),n=1,nv)
        if(imtyp.eq.1) then
          write(*,2002) subname(1:nc),(dp(n),n=1,nv)
          write(*,2003) subname(1:nc),(dp(n)*dm,n=1,nv)
          write(*,2004) subname(1:nc)
          do i = 0,nv-1,4
            do j = 1,min(4,nv-i)
              if(abs(dp(i+j)).gt.tol*10.d0) then
                dpp(j) = dr/dp(i+j)
              else
                dpp(j) = 0.0d0
              endif
            end do ! j
            write(*,2005) (dpp(j),j=1,min(4,nv-i))
          end do ! i
        endif
      endif

c     Formats

2000  format(/'  ',a,': Current eigenvalues, iteration',i4/
     +        (5x,1p,4d17.8))

2001  format(/'  ',a,': Current residuals,   iteration',i4/
     +        (5x,1p,4d17.8))

2002  format(/'  ',a,': Square root of eigenvalues (rad/t)'/
     +        (5x,1p,4d17.8))

2003  format(/'  ',a,': Square root of eigenvalues (Hz.)'/
     +        (5x,1p,4d17.8))

2004  format(/'  ',a,': Period in units of time      (T)')

2005  format(5x,1p,4d17.8)

      end
