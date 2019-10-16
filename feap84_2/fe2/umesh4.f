c$Id:$
      subroutine umesh4(tx,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add tx(*) to argument list                       26/09/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Dummy user input routine

c      Inputs:
c         tx(*)  - Command line input data
c         prt    - Flag, output results if true

c      Outputs:
c         none   - Users are responsible for generating outputs
c                  through common blocks, etc.  See programmer
c                  manual for example.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'iofile.h'
      include  'sdata.h'
      include  'umac1.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   prt,pcomp
      character tx(*)*15

c     Set command

      if(pcomp(uct,'mes4',4)) then      ! Usual    form
        uct = 'tayl'                    ! Specify 'name'
      elseif(ucount) then               ! Count elements and nodes

      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else                              ! Perform user operation

c       Set Taylor boundary conditions on RVE's

        call umsh3bc(numnp,ndf,mr(np(31)))

        if(prt) write(iow,2000)

      endif

2000  format(5x,'Taylor condition: All boundary nodes fixed')

      end

      subroutine umsh3bc(numnp,ndf,id)

c     Set all boundary codes to fixed

      implicit   none

      integer    numnp,ndf
      integer    id(ndf,numnp,2)

      id(:,:,2) = 1

      end
