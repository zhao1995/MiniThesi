c$Id:$
      subroutine pnumba(ix,ie, ib)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat1.h'
      include  'mdata.h'
      include  'sdata.h'

      integer   i,iel, ma, nd,nn
      integer   ix(nen1,*), ie(nie,*), ib(3,*)

      save

c     Zero array

      do nn = 1,numnp
        ib(1,nn) = 0
        ib(2,nn) = 0
      end do ! nn

c     Tag active nodes

      do nn = 1,numel
        do ma = 1,nummat
          if(ie(nie-2,ma) .eq. ix(nen1,nn)) then
            iel = ie(nie-1,ma)
            do i = 1,nen
              nd = ix(i,nn)
              if(nd.gt.0) then

c               User element values

                if(iel.gt.0) then
                  ib(1,nd) = ia(1,iel)
                  ib(2,nd) = ia(2,iel)

c               Program element values

                elseif(iel.lt.0) then
                  ib(1,nd) = ea(1,-iel)
                  ib(2,nd) = ea(2,-iel)
                endif ! iel
              endif ! nd
            end do ! i
          endif ! ie
        end do ! ma
      end do ! nn

c     Set third rotation node

      do nn = 1,numnp
        nd = 6 - ib(1,nn) - ib(2,nn)
        if(nd.gt.0 .and. nd.le.3) then
          ib(3,nn) = nd
        endif
      end do ! nn

      end
