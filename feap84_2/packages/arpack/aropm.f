c$Id:$
      subroutine aropm (n,u,v)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c     Author: D.S. Bindel 6/21/2005
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Mass interface to Lanczos/Arnoldi solvers

c      Inputs:
c         n    - Number of unknowns active
c         u(*) - Current vector

c      Outputs:
c         v(*) - Mass * u
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'compas.h'
      include   'evdata.h'
      include   'fdata.h'
      include   'ndata.h'
      include   'part0.h'
      include   'pointer.h'
      include   'comblk.h'

      integer    i, n
      real*8     u(n),v(n)

c     Consistent mass

      if(fl(1)) then

        do i = 1,n
          v(i) = 0.0d0
        end do ! i
        call caprod(hr(np(nx)),hr(np(nx)+n),u,v,
     &              mr(np(90)),mr(np(91)),  n)

c     Lumped mass is now active

      else

        do i = 1,n
          v(i) = hr(np(nx)+i-1)*u(i)
        end do ! i

      endif

      end
