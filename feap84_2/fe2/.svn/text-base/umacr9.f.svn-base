c$Id:$
      subroutine umacr9(lct,ctl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove 'prt' from argument list                  09/07/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Output FE^2 RVE driver files: Iname.01, etc.

c      Inputs:
c         lct       - Name of file "Iname" (limited to 15 characters)
c                     Use current job name if lct = '    '
c         ctl(1)    - Number of files to output

c      Outputs:
c         Files: Iname.01, etc. to Iname.nfile
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comfil.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'umac1.h'

      logical   pcomp
      character lct*15, ufinp*128, fext*2
      real*8    ctl(3)
      integer   n,nfile

      save

c     Set command word

      if(pcomp(uct,'mac9',4)) then      ! Usual    form
        uct = 'outf'                    ! Output 'FE2' files
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else                              ! Perform user operation

c       Specify number of files to output

        nfile = max(1,nint(ctl(1)))

        write(*,2000) nfile

        do n = 1,nfile
          if(pcomp(lct,'    ',4)) then
            ufinp = finp
          else
            ufinp = lct
          endif
          fext  = '00'
          if(n.lt.10) then
            write(fext(2:2),'(i1)') n
          else
            write(fext(1:2),'(i1)') n
          endif
          call addext(ufinp,fext,128,2)
          open(unit=ios,file = ufinp)
          write(ios,2001)
          close(unit=ios,status='keep')
        end do ! n

      endif

c     formats

2000  format(' --> Creating',i3,' RVE files')
2001  format('nocount'/'ufeap * * RVE for FE-squared'//
     &       'include solve_mpi'//'stop')

      end
