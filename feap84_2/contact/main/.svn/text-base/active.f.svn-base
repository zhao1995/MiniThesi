c$Id:$
      logical function active (var,nn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: ACTIVE history variables

c      Purpose: Allocate space in CH1, CH2, CH3 for requested variable

c      Inputs :
c         var     - VARiable name
c         nn      - # of memory location needed (.gt. 1 if vectors)

c      Outputs:
c         p1(*)   - memory location for variable inside ch1
c         p3(*)   - memory location for variable inside ch3
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_dict.h'
      include  'c_keyh.h'
      include  'c_pair.h'
      include  'iofile.h'

      logical   actfl
      character var*(*)
      integer   nn, i1,i3,ifwh,kv,kk,memgap,mingap,memvar

      save

c     Initialize counters of active variables
c     REMARK this data init work also with more than a problem in the
c            same file, because 'stop' reset always counters

      data      i1 /1/
      data      i3 /1/
      active = .false.

c     Output number of activated variables for ch1 and ch3
c     Re-initialize counters of active variables

      if (var.eq.'stop') then
        lh1  = i1 - 1
        lh3  = i3 - 1
        nset = nn                    ! Return nset in 'c_key.h'
        i1   = 1                     ! Note: Also returned as 'nn'
        i3   = 1

c       Search variables in ch1 list

      else
        ifwh = 0
        kv   = 0
        do while (ifwh.eq.0)
          kv = kv + 1
          if (kv.gt.c_lp1) then
            ifwh = 2
          elseif (var.eq.w1(kv,ndrv)) then

c           Unactivated variable

            if (p1(kv).eq.0) then
              p1(kv) = i1
              i1     = i1+nn

c           Resize previously activated variable

            else
              actfl = .true.
              do kk=1,c_lp1
                if (p1(kk).ne.0) then
                  memgap = p1(kk) - p1(kv)
                  if (memgap.gt.0) then
                    if(actfl) then
                      mingap =  memgap
                      actfl  = .false.
                    else
                      mingap = min(memgap,mingap)
                    endif
                  endif
                endif
              end do ! kk
              mingap = min((i1-p1(kv)),mingap)
              memvar = nn - mingap
              do kk=1,c_lp1
                if (p1(kk).gt.p1(kv)) then
                  p1(kk) = p1(kk) + memvar
                endif
              end do ! kk
              i1 = i1 + memvar
            endif
            ifwh = 1
          endif
        end do ! while

c       Search variable in ch3 list

        if (ifwh.eq.2) then
          ifwh = 0
          kv   = 0
          do while (ifwh.eq.0)
            kv = kv + 1
            if (kv.gt.c_lp3) then
              ifwh = 3
            elseif (var.eq.w3(kv,ndrv)) then

c             Unactivated variable

              if (p3(kv).eq.0) then
                p3(kv) = i3
                i3 = i3+nn

c             Resize previously activated variable

              else
                actfl = .true.
                do kk = 1,c_lp3
                  if (p3(kk).ne.0) then
                    memgap = p3(kk) - p3(kv)
                    if (memgap.gt.0) then
                      if(actfl) then
                        mingap =  memgap
                        actfl  = .false.
                      else
                        mingap = min(memgap,mingap)
                      endif
                    endif
                  endif
                end do ! kk
                mingap = min((i3-p3(kv)),mingap)
                memvar = nn - mingap
                do kk = 1,c_lp3
                  if (p3(kk).gt.p3(kv)) then
                    p3(kk) = p3(kk) + memvar
                  endif
                end do ! kk
                i3 = i3 + memvar
              endif
              ifwh = 1
            endif
          end do ! while
        endif

c       Variable not found

        if (ifwh.eq.3) then
          active = .true.
          write (  *,3000) var,ndrv
          write (  *,3001) (w1(kv,ndrv),kv=1,c_lp1)
          write (  *,3002)
          write (  *,3001) (w3(kv,ndrv),kv=1,c_lp3)
          write (ilg,3000) var,ndrv
          write (ilg,3001) (w1(kv,ndrv),kv=1,c_lp1)
          write (ilg,3002)
          write (ilg,3001) (w3(kv,ndrv),kv=1,c_lp3)
          call plstop()
        endif
      endif

c     Formats

3000  format (/' *ERROR* ACTIVE: Variable not defined: ',a8/
     &         ' in the variable list of contact driver #',i5//
     &         ' Variables currently defined in CH1 & CH2:')

3001  format (3x,8a9)

3002  format (/1x,'Variables currently defined in CH3:')

      end
