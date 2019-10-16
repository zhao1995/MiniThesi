c$Id:$
      subroutine xrujnt(jnt,eqrb,rcg,rlam,rjtx,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Update joint parmaters for explicit integrations

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'rigid1.h'
      include  'rjoint.h'

      integer   isw, n, i,jtype,j1,j2, jnt(6,numjts),eqrb(nrbody)
      real*8    rcg(3,11,*),rlam(3,3,6,*),rjtx(3,nrbody)
      real*8    y1(3),y2(3), rca(3)

      save

      if(isw.eq.1) then

        do n = 1,numjts
          jtype = jnt(1,n)
          if(eqrb(jnt(2,n)).eq.-5 .or. eqrb(jnt(2,n)).eq.-6) then
            if(jtype.eq.1) then ! Spherical
              j1 = jnt(2,n)
              j2 = jnt(3,n)

c             Update position by constraint

              do i = 1,3

                y1(i)       = rlam(i,1,3,j1)*(rjtx(1,n) - rcg(1,1,j1))
     &                      + rlam(i,2,3,j1)*(rjtx(2,n) - rcg(2,1,j1))
     &                      + rlam(i,3,3,j1)*(rjtx(3,n) - rcg(3,1,j1))

                y2(i)       = rlam(i,1,3,j2)*(rjtx(1,n) - rcg(1,1,j2))
     &                      + rlam(i,2,3,j2)*(rjtx(2,n) - rcg(2,1,j2))
     &                      + rlam(i,3,3,j2)*(rjtx(3,n) - rcg(3,1,j2))

c               Average changes

                rca(i)      =  rcg(i,2,j1)
                rcg(i,2,j1) =  rcg(i,2,j2) + y2(i) - y1(i)
                rca(i)      = (rcg(i,2,j1) - rca(i))*0.5d0
                rcg(i,2,j1) =  rcg(i,2,j1) - rca(i)
                rcg(i,2,j2) =  rcg(i,2,j2) - rca(i)
              end do ! i
              do i = 1,3
                rca(i)      =  rcg(i,5,j1)
              end do ! i

c             Update velocity

              rcg(1,5,j1) = rcg(1,5,j2) + rlam(2,2,5,j2)*y2(3)
     &                                  - rlam(3,2,5,j2)*y2(2)
     &                                  - rlam(2,2,5,j1)*y1(3)
     &                                  + rlam(3,2,5,j1)*y1(2)

              rcg(2,5,j1) = rcg(2,5,j2) + rlam(3,2,5,j2)*y2(1)
     &                                  - rlam(1,2,5,j2)*y2(3)
     &                                  - rlam(3,2,5,j1)*y1(1)
     &                                  + rlam(1,2,5,j1)*y1(3)

              rcg(3,5,j1) = rcg(3,5,j2) + rlam(1,2,5,j2)*y2(2)
     &                                  - rlam(2,2,5,j2)*y2(1)
     &                                  - rlam(1,2,5,j1)*y1(2)
     &                                  + rlam(2,2,5,j1)*y1(1)
              do i = 1,3
                rca(i)      = (rcg(i,5,j1) - rca(i))*0.5d0
                rcg(i,5,j1) =  rcg(i,5,j1) - rca(i)
                rcg(i,5,j2) =  rcg(i,5,j2) - rca(i)
              end do ! i
            endif
          endif
        end do ! n
      endif

      end
