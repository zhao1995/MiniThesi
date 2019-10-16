c$Id:$
      subroutine pnumc ()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Assign length to cc of 4 characters              17/04/2007
c       2. Add skipcdat() for 'tran' option                 04/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Equivalence with pnums

c      Purpose: Determine maximum contact material set, maximum
c               contact surface set, and maximum contact pairs

c      Inputs:
c                 - From read of data file(s)

c      Outputs:
c         numcs   - # of contact surfaces
c         numcm   - # of contact material laws
c         numcp   - # of contact pairs
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_contac.h'
      include  'c_dict.h'
      include  'c_dicti.h'
      include  'chdata.h'
      include  'iofile.h'
      include  'iosave.h'

      logical   errck,pcomp,readfl,savefl,gettxtd
      character cc*4,c2*15,tx(2)*15
      integer   nn,k,ns
      real*8    td(15)

      save

      call cdebug0 ('  pnumc',-1)

c     clear command array

      do k = 1, c_ncc
        cck(k) = 0
      end do

c     Scan input data

      errck = gettxtd(tx,2,td,0,'skip')
      cc    = tx(1)(1:4)
      c2    = tx(2)

c     Loop till the end command is found

      do while (.not.pcomp(cc,'end',3))

c       [surf]ace data

        if (pcomp(cc,cwd(1),4)) then
          backspace (ior)
          errck = gettxtd(tx,1,td,1,'skip')
          cc = tx(1)(1:4)
          nn = nint(td(1))

          if(nn.le.0) then
            write(ilg,3001) nn
            write(iow,3001) nn
            call plstop()
          endif
          cck(1) = cck(1)+1

c         Skip surface data

          call skipcdat ()

c       [mate]rial data

        elseif (pcomp(cc,cwd(2),4)) then
          backspace (ior)
          errck = gettxtd(tx,1,td,1,'skip')
          cc = tx(1)(1:4)
          nn = nint(td(1))

          if(nn.le.0) then
            write(ilg,3002) nn
            write(iow,3002) nn
            call plstop()
          endif
          cck(2) = cck(2)+1

c         Skip material data

          call skipcdat ()

c       [pair] data

        elseif (pcomp(cc,cwd(3),4)) then
          backspace (ior)
          errck = gettxtd(tx,1,td,1,'skip')
          cc = tx(1)(1:4)
          nn = nint(td(1))

          if(nn.le.0) then
            write(ilg,3003) nn
            write(iow,3003) nn
            call plstop()
          endif
          cck(3) = cck(3)+1

c         Skip contactpair data

          call skipcdat ()

c       [auto] data

        elseif (pcomp(cc,cwd(4),4)) then
          backspace (ior)
          errck = gettxtd(tx,1,td,1,'skip')
          cc = tx(1)(1:4)
          nn = nint(td(1))

          if(nn.lt.0) then
            write(ilg,3003) nn
            write(iow,3003) nn
            call plstop()
          endif
          cck(4) = cck(4)+1

c         Count and set the data for slide lines

          call pautoc(nn,ns)
          cck(1) = cck(1) + ns

c       [read] data from a new file

        elseif (pcomp(cc,cwd(5),4)) then
          errck = readfl(c2)

c       [save] data in a file

        elseif (pcomp(cc,cwd(6),4)) then
          errck = savefl(c2)

c       [tran]sform data in a file

        elseif (pcomp(cc,cwd(7),4)) then

          call skipcdat ()

c       [help] not active

        elseif (pcomp(cc,cwd(8),4)) then

c         Skip data

          errck = gettxtd(tx,2,td,0,'skip')
          cc = tx(1)(1:4)
          do while (.not.pcomp(cc,'    ',4))
            errck = gettxtd(tx,2,td,0,'skip')
            cc = tx(1)(1:4)
          end do

c       [unu1] user commands

        elseif (pcomp(cc,cwd(9),4)) then
          backspace (ior)
          errck = gettxtd(tx,1,td,1,'skip')
          cc = tx(1)(1:4)
          nn = nint(td(1))

          if(nn.le.0) then
            write(  *,3000) 'unu1',nn
            write(ilg,3000) 'unu1',nn
            call plstop()
          endif
          cck(9) = cck(9)+1

c         Skip user commands

          call skipcdat ()

c       [unu2] user commands

        elseif (pcomp(cc,cwd(10),4)) then
          backspace (ior)
          errck = gettxtd(tx,1,td,1,'skip')
          cc = tx(1)(1:4)
          nn = nint(td(1))

          if(nn.le.0) then
            write(  *,3000) 'unu2',nn
            write(ilg,3000) 'unu2',nn
            call plstop()
          endif
          cck(10) = cck(10)+1

c         Skip user commands

          call skipcdat ()

c       [unu3] user commands

        elseif (pcomp(cc,cwd(11),4)) then
          backspace (ior)
          errck = gettxtd(tx,1,td,1,'skip')
          cc = tx(1)(1:4)
          nn = nint(td(1))

          if(nn.le.0) then
            write(*,*) 'ERROR - Wrong unu3 Number ',nn
            call plstop()
          endif
          cck(11) = cck(11)+1

c         Skip user commands

          call skipcdat ()

c       Unused blank line

        elseif (pcomp(cc,'    ',4) .or. pcomp(cc,'cont',4)) then
          continue

c       Wrong command (interactive mode must be set)

        else
          write (*,*) 'ERROR - Wrong Command ',cc
          call plstop()
        endif

        errck = gettxtd(tx,2,td,0,'skip')
        cc    = tx(1)(1:4)
        c2    = tx(2)
      end do ! while

c     [end] of contact data

      if(lsave) then
        write (ilg,3004)
        write (iow,3004)
        call plstop()
      endif

c     Reset position in data file

      rewind ior
      errck = gettxtd(tx,1,td,0,'skip')
      cc    = tx(1)(1:4)
      do while (.not.pcomp(cc,'cont',4))
        errck = gettxtd(tx,1,td,0,'skip')
        cc    = tx(1)(1:4)
      end do ! while

c     Copy information for three basic command types

      numcs = cck(1)
      numcm = cck(2)
      numcp = cck(3)

c     Formats

3000  format(' *ERROR* PNUMC: Wrong ',a,' Number',i5)
3001  format(' *ERROR* PNUMC: Negative or missing contact surface',
     &       ' number: # =',i10)
3002  format(' *ERROR* PNUMC: Negative or missing contact material',
     &       ' number: # =',i10)
3003  format(' *ERROR* PNUMC: Negative or missing contact pair',
     &       ' number: # =',i10)
3004  format(' *ERROR* PNUMC: No SAVE,END statement for input data.')

      end
