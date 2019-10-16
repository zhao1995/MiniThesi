c$Id:$
      subroutine ploopin(lp_max)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Adjust length of 'xxx' and 'yyy' to 256          21/12/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Data input routine for mesh description

c      Inputs:
c        lp_max  - Number of times to loop in root loop
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iodata.h'
      include   'iofile.h'
      include   'setups.h'

      logical    pcomp, lopen
      character  lp_file*15,lp_ext*1, xxx*256,yyy*256
      integer    lp_lun, lp_max, lp_num, lenx, loop_l(9)

      save

c     Set parameters

      lp_file = 'feaploop000.0'
      if(rank.gt.0) then
        if(rank.lt.10) then
          write(lp_file(11:11),'(i1)') rank
        elseif(rank.lt.100) then
          write(lp_file(10:11),'(i2)') rank
        else
          write(lp_file( 9:11),'(i3)') rank
        endif
      endif

      lp_lun  =  icl - 1
      lopen   = .true.
      do while(lopen)
        lp_lun = lp_lun + 1
        inquire(unit = lp_lun,opened = lopen)
      end do ! while

c     Start root file

      lp_num         = 1
      loop_l(lp_num) = lp_lun
      open(unit = lp_lun, file = lp_file, status = 'unknown')
      rewind(lp_lun)
      write(lp_lun,1000) 'loop',lp_max

c     Fill remaining file entries

      xxx   = 'start'
      do while(.not.pcomp(xxx,'next',4) .and. lp_num.gt.0)

        read(   ior,1000) yyy
        call pstrip(xxx,yyy,1)

c       New loop encountered: Reset file structures

        if(pcomp(xxx,'loop',4)) then

c         Establish new filename and save to current file

          write(lp_ext,'(i1)') lp_num
          lp_file(13:13) = lp_ext
          yyy(1:5)       = 'file='
          yyy(6:18)      = lp_file(1:13)
          write(lp_lun,1000) yyy(1:18)

c         Establish new file structure

          lopen = .true.
          do while(lopen)
            lp_lun       = lp_lun + 1
            inquire(unit=lp_lun,opened=lopen)
          end do ! while

          open(unit = lp_lun, file = lp_file, status = 'unknown')
          rewind(lp_lun)

          lp_num         = lp_num + 1
          if(lp_num.gt.9) then
            write(ilg,3000)
            write(iow,3000)
            call plstop()
          endif
          loop_l(lp_num) = lp_lun
        endif

c       Compress amount to be written

        lenx = 80
        do while(xxx(lenx:lenx).eq.' ' .and. lenx.gt.1)
          lenx = lenx - 1
        end do ! while
        write(lp_lun,1000) xxx(1:lenx)

c       End of loop structure encountered: Set back to previous

        if(pcomp(xxx,'next',4)) then

          if(lp_num.gt.1) then
            xxx = 'restart'
          endif

          close(lp_lun)
          lp_num = lp_num - 1
          lp_lun = loop_l(lp_num)
        endif
      end do ! while

c     Format structure

1000  format(a,i5)

3000  format(' *ERROR* PLOOPIN: Mesh loop levels nested deeper than 9')

      end
