c$Id:$
      logical function matchf(factyp,nec, mi,ix)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add match for tetrahedral elements               02/09/2007
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Set matching interfaces

c     Inputs:
c       factyp - Type of interface
c       nec    - Check unumber
c       mi     - Position in list to check
c       ix(*)  - List of element nodes

c     Outputs:
c       matchf - .true. if a match occurs
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'facset.h'
      include   'iofile.h'

      integer    factyp,nec,mi,mj,mk, ix(*), me(4,6), te(3,4), jchk(3)

      save

      data   me / 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 1,2,3,4, 5,6,7,8/
      data   te / 1,2,4, 2,3,4, 3,1,4, 3,2,1 /

      if(factyp.eq.1) then
        jchk(1) = ix(mi)
        matchf  = ichk(1).eq.jchk(1)
      elseif(factyp.eq.2) then
        mj      = mod(mi,nec) + 1
        jchk(1) = min(ix(mi),ix(mj))
        jchk(2) = max(ix(mi),ix(mj))
        matchf  = ichk(1).eq.jchk(1) .and. ichk(2).eq.jchk(2)

      elseif(factyp.eq.3) then
        mj      = mod(mi,nec) + 1
        mk      = mod(mj,nec) + 1
        jchk(1) = min(ix(mi),ix(mj),ix(mk))
        jchk(2) = max(ix(mi),ix(mj),ix(mk))
        jchk(3) = ix(mi) + ix(mj) + ix(mk) - ichk(1) - ichk(2)
        matchf  = ichk(1).eq.jchk(1) .and. ichk(2).eq.jchk(2)
     &                               .and. ichk(3).eq.jchk(3)
      elseif(factyp.eq.4) then
        jchk(1) = min(ix(me(1,mi)),ix(me(2,mi)),
     &                ix(me(3,mi)),ix(me(4,mi)))
        do mj = 1,4
          if(ix(me(1,mi)).eq.ichk(1)) then
            mk      = mod(mj+1,4) + 1
            jchk(2) = ix(me(mk,mi))
            exit
          endif
        end do ! mj
        matchf  = ichk(1).eq.jchk(1) .and. ichk(2).eq.jchk(2)
      elseif(factyp.eq.5) then
        jchk(1) = min(ix(te(1,mi)),ix(te(2,mi)),ix(te(3,mi)))
        jchk(3) = max(ix(te(1,mi)),ix(te(2,mi)),ix(te(3,mi)))
        do mj = 1,3
          if((ix(te(mj,mi)) .ne. jchk(1)) .and.
     &       (ix(te(mj,mi)) .ne. jchk(3))) then
            jchk(2) = ix(te(mj,mi))
            exit
          endif
        end do ! mj
        matchf  = ichk(1).eq.jchk(1) .and.
     &            ichk(2).eq.jchk(2) .and.
     &            ichk(3).eq.jchk(3)
      else
        matchf = .false.
        write(ilg,3000) factyp
      endif

c     Format

3000  format(' *ERROR* in MATCHF: No facetyp =',i5)
      end
