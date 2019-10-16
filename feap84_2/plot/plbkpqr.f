c$Id:$
      subroutine plbkpqr(np,iel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    06/02/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set 3-D Plot Sequence for 64-node brick elements

c      Inputs:
c         np        - Order of element (3 for 64 node brick)
c         iel       - Element number: > 0 for user    elements
c                                     < 0 for program elements

c      Outputs:
c         none      - Sequesnce returned in common /pdata6/
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdata5.h'
      include  'pdata6.h'

      integer   np, iel, i
      integer   ns(3),nt,nu, ns12

      save

c     Set control variables

      do i = 1,3
        ns(i) = np + 1
      end do ! i
      ns12 = ns(1)*ns(2)

c     Set number of points

      if(iel.gt.0) then

c       Trace around bottom

        do i = 1,ns(1)
          ipord(i,iel) = i
        end do ! i
        nt = ns(1)
        nu = ns(1)

        do i = 2,ns(2)
          nt            = nt + 1
          nu            = nu + ns(1)
          ipord(nt,iel) = nu
        end do ! i

        do i = ns(1)-1,1,-1
          nt            = nt + 1
          nu            = nu - 1
          ipord(nt,iel) = nu
        end do ! i

        do i = ns(2)-1,1,-1
          nt            = nt + 1
          nu            = nu - ns(1)
          ipord(nt,iel) = nu
        end do ! i

c       Up first 3-edge

        do i = 2,ns(3)
          nt            = nt + 1
          nu            = nu + ns12
          ipord(nt,iel) = nu
        end do ! i

c       Around top first edge

        do i = 2,ns(1)
          nt            = nt + 1
          nu            = nu + 1
          ipord(nt,iel) = nu
        end do ! i

c       Down-up second 3-edge

        do i = ns(3)-1,1,-1
          nt            = nt + 1
          nu            = nu - ns12
          ipord(nt,iel) = nu
        end do ! i

        do i = 2,ns(3)
          nt            = nt + 1
          nu            = nu + ns12
          ipord(nt,iel) = nu
        end do !

c       Second top edge

        do i = 2,ns(2)
          nt            = nt + 1
          nu            = nu + ns(1)
          ipord(nt,iel) = nu
        end do ! i

c       Down-up third 3-edge

        do i = ns(3)-1,1,-1
          nt            = nt + 1
          nu            = nu - ns12
          ipord(nt,iel) = nu
        end do ! i

        do i = 2,ns(3)
          nt            = nt + 1
          nu            = nu + ns12
          ipord(nt,iel) = nu
        end do !

c       Third top edge

        do i = ns(1)-1,1,-1
          nt            = nt + 1
          nu            = nu - 1
          ipord(nt,iel) = nu
        end do ! i

c       Down-up fourth 3-edge

        do i = ns(3)-1,1,-1
          nt            = nt + 1
          nu            = nu - ns12
          ipord(nt,iel) = nu
        end do ! i

        do i = 2,ns(3)
          nt            = nt + 1
          nu            = nu + ns12
          ipord(nt,iel) = nu
        end do !

c       Last top edge

        do i = ns(2)-1,1,-1
          nt            = nt + 1
          nu            = nu - ns(1)
          ipord(nt,iel) = nu
        end do ! i

        inord(iel) = nt

      elseif(iel.lt.0) then

c       Trace around bottom

        do i = 1,ns(1)
          epord(i,-iel) = i
        end do ! i
        nt = ns(1)
        nu = ns(1)

        do i = 2,ns(2)
          nt            = nt + 1
          nu            = nu + ns(1)
          epord(nt,-iel) = nu
        end do ! i

        do i = ns(1)-1,1,-1
          nt            = nt + 1
          nu            = nu - 1
          epord(nt,-iel) = nu
        end do ! i

        do i = ns(2)-1,1,-1
          nt            = nt + 1
          nu            = nu - ns(1)
          epord(nt,-iel) = nu
        end do ! i

c       Up first 3-edge

        do i = 2,ns(3)
          nt            = nt + 1
          nu            = nu + ns12
          epord(nt,-iel) = nu
        end do ! i

c       Around top first edge

        do i = 2,ns(1)
          nt            = nt + 1
          nu            = nu + 1
          epord(nt,-iel) = nu
        end do ! i

c       Down-up second 3-edge

        do i = ns(3)-1,1,-1
          nt            = nt + 1
          nu            = nu - ns12
          epord(nt,-iel) = nu
        end do ! i

        do i = 2,ns(3)
          nt            = nt + 1
          nu            = nu + ns12
          epord(nt,-iel) = nu
        end do !

c       Second top edge

        do i = 2,ns(2)
          nt            = nt + 1
          nu            = nu + ns(1)
          epord(nt,-iel) = nu
        end do ! i

c       Down-up third 3-edge

        do i = ns(3)-1,1,-1
          nt            = nt + 1
          nu            = nu - ns12
          epord(nt,-iel) = nu
        end do ! i

        do i = 2,ns(3)
          nt            = nt + 1
          nu            = nu + ns12
          epord(nt,-iel) = nu
        end do !

c       Third top edge

        do i = ns(1)-1,1,-1
          nt            = nt + 1
          nu            = nu - 1
          epord(nt,-iel) = nu
        end do ! i

c       Down-up fourth 3-edge

        do i = ns(3)-1,1,-1
          nt            = nt + 1
          nu            = nu - ns12
          epord(nt,-iel) = nu
        end do ! i

        do i = 2,ns(3)
          nt            = nt + 1
          nu            = nu + ns12
          epord(nt,-iel) = nu
        end do !

c       Last top edge

        do i = ns(2)-1,1,-1
          nt            = nt + 1
          nu            = nu - ns(1)
          epord(nt,-iel) = nu
        end do ! i

        exord(-iel) = nt

      endif

      end
