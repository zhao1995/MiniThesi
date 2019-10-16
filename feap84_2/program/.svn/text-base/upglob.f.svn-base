c$Id:$
      subroutine upglob(du,ug)

      implicit   none

      include   'pglob1.h'

      integer    n
      real*8     du(*), ug(*)

      do n = 1, geqnum
        ug(n) = ug(n) + du(n+gneq)
      end do ! n

      end
