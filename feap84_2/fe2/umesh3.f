c$Id:$
      subroutine umesh3(tx,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Dimension td(1) in ucntnod and ucntelm           18/05/2012
c       2. Pass text(1) to pcomp checks                     15/08/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Input GiD data file

c      Inputs:
c         tx(*)  - Command line input data
c         prt    - Flag, output results if true

c      Outputs:
c         Mesh file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'chdata.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'sdata.h'
      include  'umac1.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   prt,pcomp, errck, tinput,vinput, filfl, exst
      character text(7)*15, tx(*)*15
      integer   iorsv, i,j
      real*8    td(1)

      save

c     Set command

      if(pcomp(uct,'mes3',4)) then      ! Usual    form
        uct = 'gid'                     ! Specify 'name'

      elseif(ucount) then               ! Count elements and nodes

        filfl   = .false.
        iorsv   = 0
        text(1) = 'start'
        do while(.not.pcomp(text(1),'end',3))
          errck = tinput(text,7,td,0)

c         Determine name of GiD data file

          if(pcomp(text(1),'file',4)) then
            do i = 256,1,-1
              if(xxx(i:i).ne.' ') then
                exit
              endif
            end do ! i
            do j = i,1,-1
              if(xxx(j:j).eq.'=' .or.
     &           xxx(j:j).eq.' ' .or.
     &           xxx(j:j).eq.',' ) then
                exit
              endif
            end do ! j
            j = j + 1
            inquire(file = xxx(j:i), exist = exst)
            if(exst) then
              iorsv  = ior
              ior    = ios
              filfl  = .true.
              open(unit = ior, file = xxx(j:i))
            else
              write(iow,*) ' No GiD input file'
              call plstop()
            endif

          elseif(pcomp(text(1),'mesh',4)) then
            errck = vinput(text(3),15,td(1),1)
            ndm = max(ndm,nint(td(1))-1)
            ndf = max(ndf,ndm)
            errck = vinput(text(7),15,td(1),1)
            nen = max(nen,nint(td(1)))
          elseif(pcomp(text(1),'coor',4)) then
            call ucntnod(unumnd)

          elseif(pcomp(text(1),'elem',4)) then
            call ucntelm(unumel)

c         Close GiD file and continue with inputs

c         elseif(pcomp(text(1),'end',3)) then
            ior   = iorsv
            filfl = .false.
            close(ios,status = 'keep')
            go to 100

          endif
        end do ! while

100     write(*,*) ' UNUMNP =',unumnd,' UNUMEL =',unumel

      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else                              ! Perform user operation

        filfl   = .false.
        iorsv   = 0
        text(1) = 'start'
        do while(.not.pcomp(text(1),'end',3))
          errck = tinput(text,7,td,0)

c         Determine name of GiD data file

          if(pcomp(text(1),'file',4)) then
            do i = 256,1,-1
              if(xxx(i:i).ne.' ') then
                exit
              endif
            end do ! i
            do j = i,1,-1
              if(xxx(j:j).eq.'=' .or.
     &           xxx(j:j).eq.' ' .or.
     &           xxx(j:j).eq.',' ) then
                exit
              endif
            end do ! j
            j = j + 1
            inquire(file = xxx(j:i), exist = exst)
            if(exst) then
              iorsv  = ior
              ior    = ios
              filfl  = .true.
              open(unit = ior, file = xxx(j:i))
            else
              write(iow,*) ' No GiD input file'
              call plstop()
            endif

          elseif(pcomp(text(1),'coor',4)) then
            call uinpnod(numnp,ndm,hr(np(43)),mr(np(190)))
            errck = tinput(text,1,td,0)  ! Read end coordinate record
            text = 'start'

          elseif(pcomp(text(1),'elem',4)) then
            call uinpelm(numel,nen,nen1,mr(np(33)))
            errck = tinput(text,1,td,0)  ! Read end element record
            text = 'start'

c         Close GiD file and continue with inputs

c         elseif(pcomp(text(1),'end',3)) then
            ior   = iorsv
            filfl = .false.
            close(ios,status = 'keep')
            go to 200

          endif
        end do ! while

200     write(*,*) ' NUMNP  =',unumnd,' NUMEL  =',unumel

      endif

      end

      subroutine  ucntnod(unumnd)

      implicit    none

      logical     errck, pcomp, tinput
      character   text*15
      integer     unumnd, nn
      real*8      td(1)

      save

      nn   = 0
      text = 'start'
      do while(.not.pcomp(text,'end',3))
        errck = tinput(text,1,td,0)
        nn = nn + 1
      end do ! while

      unumnd = nn - 1

      end

      subroutine  ucntelm(unumel)

      implicit    none

      logical     errck, pcomp, tinput
      character   text*15
      integer     unumel, nn
      real*8      td(1)

      save

      nn   = 0
      text = 'start'
      do while(.not.pcomp(text,'end',3))
        errck = tinput(text,1,td,0)
        nn = nn + 1
      end do ! while

      unumel = nn - 1

      end

      subroutine  uinpnod(numnp,ndm,x,ntyp)

      implicit    none

      include    'iofile.h'

      logical     errck, pinput
      integer     numnp, ndm, i,nn, ntyp(*)
      real*8      x(ndm,numnp), td(4)

      save

      do nn = 1,numnp
        errck = pinput(td,ndm+1)
        do i = 1,ndm
          x(i,nn) = td(i+1)
        end do ! i
        write(iow,2000) nn,(x(i,nn),i=1,ndm)
        ntyp(nn) = 0
      end do ! nn

2000  format(i8,1p,3e15.6)

      end

      subroutine  uinpelm(numel,nen,nen1, ix)

      implicit    none

      logical     errck, pinput
      integer     numel,nen,nen1, nn, ii
      integer     ix(nen1,numel)
      real*8      td(nen+2)

      save

      do nn = 1,numel
        errck = pinput(td,nen + 2)
        do ii = 1,nen
          ix(ii,nn) = nint(td(ii+1))     ! Nodal connection list
        end do ! ii
        ix(nen1,nn) = nint(td(nen+2))    ! Material number
      end do ! nn

      end
