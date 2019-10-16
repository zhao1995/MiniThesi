c$Id:$
      subroutine comgeq(ir, kp, bycol,all, neq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    13/12/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Compute global equation numbers from elements

c      Inputs:
c         kp         -  Initial row entry
c         bycol      -  Storage by columns if true
c         all        -  All terms in row/col if true
c         neq        -  Total number of equations

c      Outputs:
c         ir(*)      -  Row number of each nonzero in stiffness matrix.
c         kp         -  Last entry in ir array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'part0.h'
      include  'pglob1.h'

      logical   addeq, bycol, all
      integer   neq, kp, i
      integer   ir(*)

      save

      if(gpart.eq.0 .or. npart.eq.gpart) then

c       Check if equation to be added

        if(all) then                           ! all terms
          addeq   = .true.
        elseif(bycol) then                     ! by columns
          addeq   = .false.
        else                                   ! by rows
          addeq   = .true.
        endif

c       Add equation to list

        if(addeq) then

c         New equation, add to list

          do i = gneq+1,neq
            kp     = kp + 1
            ir(kp) = i
          end do ! i
        endif

      endif

      end
