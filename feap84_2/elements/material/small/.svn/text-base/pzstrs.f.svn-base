c$Id:$
      subroutine pzstrs(d,flux,eps,sig,dd)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add computation of psi from xref                 26/07/2009
c-----[--.----+----.----+----.-----------------------------------------]
c     Linear Piezo-electric Elastic Constitutive Model

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'elcoor.h'

      real*8          dc     ,cd     ,cc
      common /elmatl/ dc(6,3),cd(3,6),cc(3,3)

      integer   i,j
      real*8    psi, d(*), eps(6), sig(6), dd(6,6), beta(6), flux(3)

      save

c     Stress:

      if(nint(d(85)).gt.0) then
        psi = atan2(xref(2)-d(87),xref(1)-d(86))
      else
        psi = d(31)
      endif
      call dmat2d(d,psi,dd,beta)

      do i = 1,6
        sig(i) = 0.0d0
        do j = 1,6
          sig(i) = sig(i) + dd(i,j)*eps(j)
        end do ! j
      end do ! i

c     Add piezo-electric terms

      sig(1) = sig(1) - d(151)*flux(2)
      sig(2) = sig(2) - d(152)*flux(2)
      sig(3) = sig(3) - d(151)*flux(2)
      sig(4) = sig(4) - d(153)*flux(1)
      sig(5) = sig(5) - d(153)*flux(3)

c     Set up coupling and permeability array

      do i = 1,3
        do j = 1,6
          dc(j,i) = 0.0d0
          cd(i,j) = 0.0d0
        end do ! j
        do j = 1,3
          cc(j,i) = 0.0d0
        end do ! j
        cc(i,i) = d(156+i)
      end do ! i

      dc(1,2) = d(151)
      dc(2,2) = d(152)
      dc(3,2) = d(151)
      dc(4,1) = d(153)
      dc(5,3) = d(153)

      cd(2,1) = d(151)
      cd(2,2) = d(152)
      cd(2,3) = d(151)
      cd(1,4) = d(153)
      cd(3,5) = d(153)

      end
