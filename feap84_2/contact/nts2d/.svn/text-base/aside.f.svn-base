c$Id:$
      subroutine aside (ix1,nw,sorts,sidec)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Find history data sets around a given one

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_contac.h'
      include  'c_geom.h'
      include  'c_keyh.h'
      include  'c_mate.h'
      include  'c_pair.h'
      include  'c_tole.h'

      logical   unfill
      integer   ke,kc,kf,eleft,erigh,center
      integer   ix1(dnope1,*),nw,sorts(*),sidec(nset,*)

      save

c     Search left end element along contact surface

      ke = 1
      eleft = ix1(dnope1-1,ke)
      do while (eleft.ne.0)
        ke = ke+1
        eleft = ix1(dnope1-1,ke)
      end do

c     Form sorted vector

      erigh = ke
      do ke = 1,neps1
        sorts(ke) = erigh
        erigh = ix1(dnope1,erigh)
      end do

c     add last data set (it is always the last one)

      sorts(nset) = nset

c     Form table
c     Find center within the vector

      do ke = 1,nset
        kc = 1
        do while (sorts(kc).ne.ke)
          kc = kc+1
        end do
        center = kc

c       Store center and left-right

        sidec(ke,1) = sorts(center)
        unfill = .true.
        kc = 0
        kf = 1
        do while (unfill)
          kc = kc + 1
          if ((kc.gt.nset) .or. (kf.eq.nw)) then
            unfill = .false.
          else
            if ((center+kc).le.nset) then
              kf = kf+1
              sidec(ke,kf)   = sorts(center+kc)
            endif
            if ((center-kc).ge.1) then
              kf = kf+1
              sidec(ke,kf)   = sorts(center-kc)
            endif
          endif
        end do
      end do

c     Distance criterion has to be introduced

      if (ifdb) then
        write (69,*) 'SORTS NW NSET',nw,nset
        write (69,4000) (sorts(ke),ke=1,nset)

        write (69,*) 'SIDEC MATRIX'
        do ke=1,nset
          write (69,4000) (sidec(ke,kc),kc=1,nw)
        end do
      endif

4000  format (11i5)

      end
