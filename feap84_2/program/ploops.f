c$Id:$
      subroutine ploops(lp_in,txt,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Control program for mesh LOOP-NEXT commands

c      Inputs:
c         lp_in      - Logical for starting loop
c         txt        - Text storing loop range value
c         isw        - Switch: isw = 1 for LOOP
c                              isw = 2 for NEXT

c      Outputs:
c         Depends on commands specified
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'iodata.h'
      include   'iofile.h'
      include   'setups.h'

      logical    lp_in, lopen
      character  txt*15, txl*15, lp_file*15, lp_fil0*15
      integer    isw, lp_lun,lp_num, lp_ior(0:9),lp_cur(9),lp_max(9)
      real*8     loopval

      save

c     [loop,#] - Loop start

      if(isw.eq.1) then

        if(lp_in) then

          lp_in     = .false.

c         Check rank number

          lp_fil0 = 'feaploop000.0'
          if(rank.gt.0) then
            if(rank.lt.10) then
              write(lp_fil0(11:11),'(i1)') rank
            elseif(rank.lt.100) then
              write(lp_fil0(10:11),'(i2)') rank
            else
              write(lp_fil0( 9:11),'(i3)') rank
            endif
          endif

c         Set loopvalue

          call setval(txt,15,loopval)
          lp_max(1) = max(1,nint(loopval))

c         Construct file(s) for loop commands

          call ploopin(lp_max(1))
          lp_ior(0) =  ior

          lopen     = .true.
          lp_lun    =  icl - 1
          do while(lopen)
            lp_lun  =  lp_lun + 1
            inquire(unit = lp_lun, opened = lopen)
          end do ! while

          ior       = lp_lun
          lp_file   = lp_fil0
          open(unit = ior, file = lp_file,status = 'old')
          rewind(ior)

          lp_num    = 0

        else

          call setval(txt,15,loopval)
          lp_num         = lp_num + 1
          lp_max(lp_num) = max(1,nint(loopval))
          lp_ior(lp_num) = ior
          lp_cur(lp_num) = 1

        endif

c     [next] - Loop end

      elseif(isw.eq.2) then

        if(.not.lp_in) then

          lp_cur(lp_num) = lp_cur(lp_num) + 1

          if(lp_cur(lp_num).gt.lp_max(lp_num)) then

            close(ior)
            lp_num = lp_num - 1
            ior    = lp_ior(lp_num)
            lp_lun = lp_lun - 1
            if(lp_num.eq.0) then
              lp_in = .true.
            endif

          else

            rewind(ior)
            read(ior,'(a)') txl ! prevent a re-read of the loop

          endif

c       Error

        else

          write(iow,2000)
          call plstop()

        endif
      endif

c     Output format

2000  format(/5x,'LOOP-NEXT error in PMESH: NEXT before LOOP.'/)

      end
