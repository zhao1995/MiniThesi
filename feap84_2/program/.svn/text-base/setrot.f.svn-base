c$Id:$
      subroutine setrot ( x , mo, xlm , thk , numnp , option )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c      1. Recode computation of orthogonal matrix           01/06/2013
c      2. Add option for generation of lambda matrices      05/06/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Obtain directors from the original mesh x(1-3,*) and
c               second mesh xlm(1-3,3,*), normalize and set trans-
c               formation using the exponential map. The norm is
c               multiplied by 2 and stored in thk(*).

c      Inputs:
c         x(3,*)  - Nodal coordinates for mesh
c         mo(*)   - Rotational update type
c         numnp   - Number of nodes in mesh
c         option  - Rotation Matrix type

c      Outputs:
c         xlm(*)  - Normalized nodal directors
c         thk(*)  - Nodal thicknesses
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,n, numnp, mo(*)
      integer   imin, imax, option
      real*8    vmin, vmax, fac, tol
      real*8    x(3,*),xlm(9,6,*),thk(*),v(3)

      save

      data      tol / 1.0d-10 /

c     Loop over nodes

      do n = 1 , numnp

        if(mo(n).eq.-1 .or. mo(n).eq.-5) then

c         Compute Fiber:

          v(1) = xlm(7,1,n) - x(1,n)
          v(2) = xlm(8,1,n) - x(2,n)
          v(3) = xlm(9,1,n) - x(3,n)

c         Normalize:

          vmax = sqrt ( v(1)**2 + v(2)**2 + v(3)**2 )
          if (vmax.gt.1.d-6) then
            xlm(7,1,n) = v(1) / vmax
            xlm(8,1,n) = v(2) / vmax
            xlm(9,1,n) = v(3) / vmax
            thk(n)     = vmax * 2.d0
          else
            xlm(7,1,n) = 0.d0
            xlm(8,1,n) = 0.d0
            xlm(9,1,n) = 1.d0
            thk(n)     = 1.d0
          endif

c         Assemble [Lambda]:

          v(1) = xlm(7,1,n)
          v(2) = xlm(8,1,n)
          v(3) = xlm(9,1,n)

          if(option.eq.1) then

c           Compute location of maximum and minimum entries in v(:)

            vmin = abs(v(1))
            vmax = vmin
            imin = 1
            imax = 1
            do i = 2,3
              if(abs(v(i)).lt.vmin) then
                vmin = abs(v(i))
                imin = i
              endif
              if(abs(v(i)).gt.vmax) then
                vmax = abs(v(i))
                imax = i
              endif
            end do ! i

c           Location of intermediate component

            i = 6 - imin - imax

c           Compute first surface vector: xlm(:1)

            xlm(imax,1,n) = -v(imin)
            xlm(i   ,1,n) =  0.0d0
            xlm(imin,1,n) =  v(imax)

c         Use 1-2 components

          elseif(option.eq.2) then

            if(max(abs(v(1)),abs(v(2))).gt.1.d-03) then
              xlm(1,1,n) =  v(2)
              xlm(2,1,n) = -v(1)
              xlm(3,1,n) =  0.0d0
            else
              xlm(1,1,n) =  1.0d0
              xlm(2,1,n) =  0.0d0
              xlm(3,1,n) =  0.0d0
            endif

c         Original algorithm

          elseif(option.eq.3) then

            if (abs(v(3)).gt.tol) then
              fac        =  1.d0 / ( 1.d0 + v(3) )
              xlm(1,1,n) =  v(3) + fac*v(2)*v(2)
              xlm(2,1,n) =       - fac*v(1)*v(2)
              xlm(3,1,n) = -v(1)
c             xlm(4,1,n) =       - fac*v(2)*v(1)
c             xlm(5,1,n) =  v(3) + fac*v(1)*v(1)
c             xlm(6,1,n) = -v(2)
            else
              fac        =  1.d0 / ( 1.d0 - v(3) )
              xlm(1,1,n) = -v(3) + fac*v(2)*v(2)
              xlm(2,1,n) =       - fac*v(1)*v(2)
              xlm(3,1,n) =  v(1)
c             xlm(4,1,n) =         fac*v(2)*v(1)
c             xlm(5,1,n) =  v(3) - fac*v(1)*v(1)
c             xlm(6,1,n) = -v(2)
            endif

          endif

c         Scale first vector

          vmax = 1.d0/sqrt(xlm(1,1,n)**2
     &                    +xlm(2,1,n)**2
     &                    +xlm(3,1,n)**2)
          xlm(1,1,n) = xlm(1,1,n)*vmax
          xlm(2,1,n) = xlm(2,1,n)*vmax
          xlm(3,1,n) = xlm(3,1,n)*vmax

c         Compute second surface vector: xlm(:2)

          xlm(4,1,n) = v(2)*xlm(3,1,n) - v(3)*xlm(2,1,n)
          xlm(5,1,n) = v(3)*xlm(1,1,n) - v(1)*xlm(3,1,n)
          xlm(6,1,n) = v(1)*xlm(2,1,n) - v(2)*xlm(1,1,n)

          do i = 1,9
            xlm(i,2,n) = xlm(i,1,n)
            xlm(i,3,n) = xlm(i,1,n)
            xlm(i,6,n) = xlm(i,1,n)
          end do ! i
        endif
      end do ! n

      end
