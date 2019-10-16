c$Id:$
      subroutine sisorti(nl,list,nk,key)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Order a list of integers in ascending sequence of their
c              keys using insertion sort

c     Input:
c     ------

c     nl   - length of list
c     list - a list of integers
c     nk   - length of key (nk must be ge nl)
c     key  - a list of integer keys

c     Output:
c     -------

c     nl   - unchanged
c     list - a list of integers sorted in ascending sequence of key
c     nk   - unchanged
c     key  - unchanged

c     Note:    efficient for short lists only (nl lt 20)
c     -----

c     Programmer:             Scott Sloan
c     -----------

c     Last modified:          1 March 1991     Scott Sloan
c     --------------

c     Copyright 1989:         Scott Sloan
c     ---------------         Department of Civil Engineering
c                             University of Newcastle
c                             NSW 2308
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      integer    i,j,t
      integer    nl,nk
      integer    value
      integer    key(nk)
      integer    list(nl)

      save

      do i = 2,nl
        t     = list(i)
        value = key(t)
        do j = i-1,1,-1
           if(value.ge.key(list(j)))then
             list(j+1) = t
             goto 20
           endif
           list(j+1) = list(j)
        end do ! j
        list(1) = t
   20   continue
      end do ! i

      end
