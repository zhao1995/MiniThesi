c$Id:$
      logical function naninfck(y,sizey,isw)

c     * * F E A P * * A Finite Element Analysis Program

c.... Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Check if double precision number is NAN or -INF or +INF
c              This is Fortran  replacement for C: isinf
c                                                  isnan
c       Inputs:
c          y(*)     : The array of values to be checked
c          sizey    : Length of list
c          isw      : Switch - 0: Initialize values; >0: Check values

c       Outputs:
c         naninfck  : A logical variable which is true or false
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      real*8   y(*)
      integer  sizey,isw,ic
      real*8   zero,Pinf,Minf

      save

c     Set values for checks

      if(isw.eq.0) then
        zero =  0.d0
        Pinf =  1.d0/zero   ! + Inf
        Minf = -1.d0/zero   ! - Inf
        naninfck = .true.

c     Check values

      else

        naninfck = .false.
        do ic = 1,sizey

c         if ( isnan(y(ic)) .ne. 0) then                  ! C  interface
          if(((y(ic) > 0.0) .eqv. (y(ic) <= 0.0))         ! NaN check
     &                      .and. (y(ic).ne.y(ic)) ) then
            naninfck = .true.
            return
c         elseif ( isinf(y(ic)) .ne. 0) then              ! C interface
          elseif((y(ic).eq.Pinf).or.(y(ic).eq.Minf)) then ! Inf check
            naninfck = .true.
            return
          end if
        end do ! ic

      endif ! isw

      end
