c$Id.$
      integer function bserchi(ivect,ivlen, node)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Binary search for node number 'node'

c     Inputs:
c       ivect(*) - Array to search
c       ivlen    - Length of array
c       node     - Number to locate

c     Output:
c       bserchi  - Location in ivect for 'node'
c                  0 = Not in list
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit    none

      integer     ivlen, node, is,ie,im

      integer     ivect(ivlen)

c     Check first and last node

      if(node.lt.ivect(1) .or. node.gt.ivect(ivlen)) then
        bserchi = 0

c     Serch list

      else
        is = 1
        ie = ivlen
        do while (ie.gt.is)

          if(node.eq.ivect(is) ) then
            bserchi = is
            return
          elseif(node.eq.ivect(ie) ) then
            bserchi = ie
            return
          else
            im = (ie+is)/2
            if(node.ge.ivect(is) .and. node.le.ivect(im)) then
              ie = im
            elseif(node.ge.ivect(im+1) .and. node.le.ivect(ie)) then
              is = im+1
            else
              bserchi = 0
              return
            endif
          endif
        end do ! while
        bserchi = 0
      endif

      end
