c$Id:$
      subroutine bm2con(d,hn,h1,nh, cc,strs,def,def1,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set constitutive model

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'hdata.h'
      include  'comblk.h'

      integer   isw, i, nh,nlay
      real*8    d(*),hn(*),h1(*), cc(3,3,2),strs(3,2),def(3,2),def1(3)
      real*8    dl(3)

      save

c     Resultant model

      nlay = int(d(101))

      if(nlay.eq.0) then
        dl(1) = d(1)*d(32)
        dl(2) = d(1)*d(32)*d(37)*0.5d0/(1.d0+d(2))
        dl(3) = d(1)*d(33)

        do i = 1,9
          cc(i,1,1) = 0.0d0
          cc(i,1,2) = 0.0d0
        end do ! i

c       Elastic resultant model only

        do i = 1,3
          cc(i,i,1) = dl(i)
          cc(i,i,2) = dl(i)
          strs(i,1) = dl(i)*def(i,1)
          strs(i,2) = dl(i)*def(i,2)
        end do ! i

        if(nint(d(160)).eq.2) then ! Add augmented value
          strs(1,1) = strs(1,1) + hr(nh2)
          if(isw.eq.10) then ! Update for augmentation
            hr(nh2) = strs(1,1)
          endif
        endif

c     Thickness quadrature model

      else

        if(isw.eq.12) then
          call bm2ups(nlay,d,d(103),hn,h1,nh,def1, isw)
        else
          call bm2res(nlay,d,d(103),hn,h1,nh,def, strs,cc, isw)
        endif

      endif

      end
