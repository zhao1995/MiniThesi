c$Id:$
      subroutine pmacr8 (lct,ct,j)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Command language instruction subprogram: Part 8 for NURBS

c     Inputs:
c        lct      - Command option for current command
c        ct(3)    - Command parameters for current command
c        j        - Number of command to execute

c     Outputs:
c        Depends on value of command j: None for serial version
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      integer    j
      character  text*15,lct*15
      real*8     ct(3)

      save

c     Solution command instruction subprogram - part 8.

      text = lct

      if(nint(ct(1)).ge.0) then

c       [extract block n_blk, k_dir inc_order
c       [extract patch n_blk, k_dir inc_order

        if(j.eq.1) then

          write(*,2000) ' EXTRact '

c       [insert block n_blk, k_dir u_knot n_times
c       [insert patch n_blk, k_dir u_knot n_times

        elseif(j.eq.2) then

          write(*,2000) ' INSErt '

        endif
      endif

c     Formats

2000  format('  *ERROR*',a,'command only available in NURBS version')

      end
