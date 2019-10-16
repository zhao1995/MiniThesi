c$Id:$
      subroutine plinka(filx,name,term)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Adjust length of 'yyy' to 256                    21/12/2008
c       2. Write ldnum,ldprf,ldflg                          02/01/2009
c       3. Write spnum,spflg                                09/03/2009
c       4. Add strip Ctrl-I from input records              07/02/2010
c       5. Set length of file extension to 8                25/01/2011
c       6. Remove set 'iclink = 0' at line = 68             28/10/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Read data and save on a file for future processing

c      Inputs:
c         filx      - Name of file extender for save file
c         name      - Type of data: 'set' or 'add'
c         term      - Terminator name: 'end' or '   '

c      Outputs:
c         none      - Data saved to disk
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdat2.h'
      include  'comfil.h'
      include  'conval.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'ioincl.h'
      include  'linka.h'
      include  'pload1.h'
      include  'print.h'
      include  'trdata.h'

      logical   pcomp,lsave,pinput,lopen
      character filx*(*),name*(*),term*(*)
      character fnamr*128,fext*8,yyy*256,type*4,fend*3
      integer   i,n
      real*8    td(1)

      save

      fnamr =  fsav
      fext  =  filx
      fend  =  term
      call addext(fnamr,fext,128,8)
      inquire(unit=iop,opened=lopen)
      if(lopen) then
        write(*,*) 'PLINKA: UNIT',iop,' is already open'
        call plstop()
      endif
      call opnfil(fext,fnamr,-1,iop,lsave)

c     Save values of current parameters

      type = name
      write(iop,1002) type,fincld(isf),irecrd(isf),prt,prth
      write(iop,1001) vvv
      write(iop,1001) tr,xr,trdet,x0
      write(iop,1003) ldnum,ldprp,spnum,ldflg,spflg
  10  if(ior.gt.0) then
        read(ior,1000,err=901,end=902) record
        irecrd(isf) = irecrd(isf) + 1
      else
        read(  *,1000) record
      endif
      yyy = record
      do n = 1,256
        if(ichar(yyy(n:n)).eq.13) then      ! Strip Ctrl-M
          yyy(n:n) = ' '
        elseif(ichar(yyy(n:n)).eq. 9) then  ! Strip Ctrl-I
          yyy(n:n) = ' '
        elseif(yyy(n:n).eq.'!') then        ! Strip comments
          yyy(n:256) = ' '
          go to 20
        endif
      end do ! n
      n = 256
  20  do i = n,1,-1
        if(yyy(i:i).ne.' ') go to 30
      end do ! i
      i = 1
  30  write(iop,1000) yyy(1:i)

      if(i.gt.1) then
        n = i
        do i = 1,n
          if(yyy(i:i).ne.' ') go to 40
        end do ! i
        i = 1
      endif
  40  if(.not.pcomp(yyy(i:i+2),fend,3)) then
        iclink = iclink + 1
        go to 10
      else
        close(iop)
      endif
      return

c     Read Error

 901  call errclr('PLINKA')
      return

c     Read EOF

 902  if(ior.eq.icf) then
        lsave = pinput(td,1)
        close(iop)
      else
        call endclr('PLINKA',yyy)
      endif

c     Format

1000  format (a)
1001  format (1p,4e20.12)
1002  format (a4,2x,a12,i8,2l5)
1003  format (3i8,2l3)

      end
