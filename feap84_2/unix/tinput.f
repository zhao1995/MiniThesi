c$Id:$
      logical function tinput(tx,mt,d,nn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Use 'cinput' for keyboard inputs                 18/04/2005
c       2. Add set to '.true.' at line 59                   16/03/2008
c       3. Limit length of character array to 15            10/05/2008
c       4. Set lrec to 256                                  21/12/2008
c       5. Remove deletion of Ctrl-I (done in pstrip)       07/02/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Input routine for real data.: returns true on error
c               Data input device: Returns true on error

c      Inputs:
c         mt        - Number of text data items to extract from input
c         nn        - Number of real data items to extract from input
c                     N.B. Input performed by this function

c      Outputs:
c         tx(*)     - Values of text data input
c         d(*)      - Values of real data input
c         tinput    - Flag, returns true if error occurs during input
c-----[--.----+----.----+----.-----------------------------------------]
      implicit    none

      include    'chdata.h'
      include    'comfil.h'
      include    'iodata.h'
      include    'ioincl.h'
      include    'iofile.h'
      include    'iosave.h'

      logical     cinput
      logical     vinput, pcomp, cksep, first

      integer     i,j,k,mc,mm,mt,nn,lrec, iskip
c     21=19 Input variables      
      character   txl*15,tx(*)*15,tl(21)*15
      real*8      d(*)

c     Parameters
      integer     cw    ! Field width for parsing
      integer     fw    ! Field width for parsing
      integer     rec   ! Record length

      save

      data        lrec /256/
      data        cw   / 15/
      data        fw   / 15/
      data        rec  / 80/

c     Check on number of items, !! change, if number is changed !!


      if(mt+nn.gt.21) then
        tinput = .true.
        if(ior.lt.0) then
          write(*,2000) mt+nn
          return
        else
          write(iow,2000) mt+nn
          write(ilg,2000) mt+nn
          call plstop()
        endif
      else
        tinput = .false.
      endif

c     Initialize arrays

10    tx(1) = ' '
      do j = 1,nn
        d(j) = 0.0d0
      end do ! j
      record = ' '

11    if(ior.gt.0) then
        read (ior,1000,err=901,end=902) record
        irecrd(isf) = irecrd(isf) + 1
        iskip       = 1
      else
c       read (  *,1000,err=901,end=902) record
        if(.not.cinput()) then
          goto 902
        endif

        iskip       = 3
      endif

c     Strip control characters, leading blanks and comments

      yyy = record
      call pstrip(xxx,yyy,iskip)

12    do k = lrec,1,-1
        if(xxx(k:k).ne.' ') go to 13
      end do ! k
      k = 1
13    if(lsave) write(lfile,1000) xxx(1:k)
      if(lmate) write(  iwd,1001) xxx(1:k)

c     Load character variables

      if(pcomp(xxx,'incl',4)) then
        mm = max(2,mt)
      else
        mm = mt
      endif

      tl(1) = '    '
      if(mm.gt.0) then
        mc = 1
        txl = ' '
        j   = 0
        first = .false.

c       String text between double quotes (")

        do i = 1,lrec
          if(xxx(i:i).eq.'"' .or. first) then
            if(first) then
              if(xxx(i:i).eq.'"') then
                first    = .false.
              else
                j        = j + 1
                if(j.le.cw) then
                  txl(j:j) = xxx(i:i)
                endif
              endif
            else
              first      = .true.
            endif

c         Non-separator string data

          elseif(.not.cksep(xxx(i:i))) then
            j          = j + 1
            if(j.le.cw) then
              txl(j:j)   = xxx(i:i)
            endif

c         Separator encountered: save character string data

          else
            tl(mc)     = txl
            txl        = ' '
            j          = 0
            mc         = mc + 1
            if(mc.gt.mm) go to 14
          endif
        end do ! i
      else
        i = 0
      end if

14    do j = 1,mt
        tx(j) = tl(j)
      end do ! j

c     Change to an include file

      if(pcomp(tl(1),'incl',4)) then
        if(  pcomp(tl(2),'end',3) ) then
          if(ior.eq.icf) then
            call pincld('end')
          endif
          return
        else
          call pincld(tl(2))
          go to 10
        endif
      endif

c     Finish inputs for parameters

      call acheck(xxx,yyy,fw,rec,rec)
      zzz = xxx(i+1:rec)
      if(nn.gt.0 .and. .not.pcomp(tl(1),'proc',4)) then
        tinput = vinput(xxx(i+1:lrec),lrec-i,d,nn)
      else
        tinput = .false.
      end if

      return

c     Read error encoutered

901   call  errclr ('TINPUT')
      goto  11

c     EOF encountered

902   if(eofile) then
        tinput = .true.
      elseif(ior.eq.icf) then
        call pincld('end')
        tinput = .false.
        tx(1)  = ' '
      else
        call  endclr ('TINPUT',yyy)
        goto  12
      endif

c     Formats

1000  format(a)
1001  format(4x,a)

2000  format(' *ERROR* TINPUT: Too many items requested, limit = 21:',
     &       ' Requested',i8)

      end
