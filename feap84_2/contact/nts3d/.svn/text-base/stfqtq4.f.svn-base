c$Id:$
      subroutine stfqtq4 (cp0,cm,ch2,ch1,ch3,tanm,resv)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0
c               Anna Haraldsson             February 1998            1.

c      Acronym: Contact MATerial Law 1

c      Purpose: compute parameters related to contact material model

c      Inputs :
c         cp0(*)  - Contact pair control data
c         cm(*)   - Contact materials data storage
c         ch1(*)  - Contact history variables (old)
c         ch2(*)  - Contact history variables (current)

c      Outputs:
c         ch2(*)  - Contact history variables (current)

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_contac.h'
      include  'c_keyh.h'
      include  'c_mate.h'
      include  'c_pair.h'
      include  'counts.h'
      include  'iofile.h'
      include  'ndata.h'

      integer   i,j,k,l,istgn,inew,jnew,istgt,istgnold
      real*8    cp0(nr0,n0c3:*),cm(*),ch1(*),ch2(*),ch3(*)
      real*8    tanm(18,18),resv(3,6)
      real*8    epsn,gn,fn,nvec(3),epst,ni(2),delgn(3,3),mue
      real*8    temp1

      save

      call cdebug0 ('     stfqtq4',-1)

      istgn    = nint(ch2(p1(4)))
      istgnold = nint(ch1(p1(4)))

c     Get geometrical variables

      gn      =  ch3(p3(9))
      nvec(1) =  ch3(p3(18))
      nvec(2) =  ch3(p3(18)+1)
      nvec(3) =  ch3(p3(18)+2)
      ni(1)   =  1.d0
      ni(2)   = -1.d0

      do j = 1,2
        do i = 1,3
          delgn(i,j) = nvec(i)*ni(j)
        end do ! i
      end do ! j

      epsn = cp0(3,2)
      fn   = gn*epsn
      if(ifsolm.eq.2) then
        fn        =  fn + ch2(p1(31))
        resv(1,6) = -gn*ch3(p3(16))
        do i = 1,3
          do j = 1,2
            inew          = 3*(j-1) + i
            tanm(inew,16) = delgn(i,j)*ch3(p3(16))
            tanm(16,inew) = tanm(inew,16)
          end do ! j
        end do ! i
      endif
      ch3(p3(14)) = fn

      temp1 = epsn*ch3(p3(16))
      do i = 1,3
        do j = 1,2
          inew = (3*(j-1))+i
          do k = 1,3
            do l = 1,2
              jnew            = (3*(l-1))+k
              tanm(inew,jnew) = temp1*delgn(i,j)*delgn(k,l)
            end do ! l
          end do ! k
        end do ! j
      end do ! i

      temp1 = temp1*gn
      do j = 1,2
        do i = 1,3
          resv(i,j) = -temp1*delgn(i,j)
        end do ! i
      end do ! j

      if (iffric.eq.1) then

        epst        = cp0(4,2)*ch3(p3(16))
        mue         = cm(1)
        istgt       = nint(ch2(p1(3)))
        ch3(p3( 5)) = 0.d0
        ch3(p3(15)) = 0.d0

        if (istgt.eq.1) then

          if(istgnold.eq.0) then

          else

            do j = 1,2
              do i = 1,3
                resv(i,j) = resv(i,j) - ch2(p1(26)+i-1)*epst*ni(j)
              end do ! i
            end do ! j

            ch3(p3( 5)) =-ch2(p1(26)+2)*epst*ni(1)
            ch3(p3(15)) = sqrt(ch2(p1(26  ))**2
     &                       + ch2(p1(26)+1)**2
     &                       + ch2(p1(26)+2)**2)*epst
          endif

        endif                  ! istgt.eq.1
      endif                    ! iffric.eq.1

      end
