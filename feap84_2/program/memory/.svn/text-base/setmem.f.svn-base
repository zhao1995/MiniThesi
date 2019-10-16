c$Id:$
      logical function setmem(list,mlist,rlist,
     &                        num,name,length,precis)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2012: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set length of adr to 600                         19/11/2007
c       2. Move 'adr' to include w_int for type setting     25/01/2010
c       3. Redefine 'irp' array and use for various ipr     13/06/2010
c       4. Use 'point' on pointer loops                     21/03/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Define, delete, or resize a dictionary entry.
c               Pointer defined for integer (single) and real
c               (double precision arrays.

c      Inputs:
c         list       - Number of entries in variables
c         mlist(2,*) - Location entries for defined arrays
c         rlist(*)   - Admissible names for arrays
c         num        - Entry number for array (see below)
c         name       - Name of array          (see below)
c         length     - Length of array defined: =0 for delete
c         precis     - Precision of array: 1 = integers; ipr = reals
c                      N.B. if ipr = 1, all arrays padded by one word
c                           to pervent overlaps.

c      Outputs:
c         np(num)    - Pointer to first word of array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'allotd.h'
      include  'allotn.h'
      include  'cdata.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'psize.h'
      include  'errchk.h'
      include  'pointer.h'
      include  'comblk.h'
      include  'w_int.h'

      include  'p_point.h'

      logical   pcomp
      integer   list,num,length,precis,ip,ipa, n,i, iot, lensav
      integer   dicloc, mlist(2,list), irp(2,2)
      character name*(*),dname*5, rlist(list)*5

      save

      data      irp / 8,8,4,8 /, iot / 9 /

c     Find variable

      mmax = 0

      dname  = name
      n      = max(1,min(num,list))

c     Check match of number with name of array

      if(pcomp(dname,rlist(n),5)) then

c       Set pointer for arrays stored in Blank Common

        setmem = .true.

        if(mlist(1,n).eq.0) then
          if(length.le.0) then
            write(  *,*) '  *WARNING* Length allocation for:',dname,
     &                   ' Length =',length
            write(iow,*) '  *WARNING* Length allocation for:',dname,
     &                   ' Length =',length
          endif
          ip     = abs(precis)
          ipa    = irp(ip,ipr)

c         Use Malloc to allocate space for length*ipa bytes

          adr(n) = malloc(length*ipa)

c         Set pointer for array use

          if(ip.eq.1) then
            np(n) = 1 + (adr(n) - loc(mr(1))) / ipa
          else
            np(n) = 1 + (adr(n) - loc(hr(1))) / ipa
          endif

c         Add new array into dictionary

          ndict      = ndict+1
          mlist(1,n) = ndict
          mlist(2,n) = length

          if(ndict.le.200) then
            dict(ndict)   = dname
            ipoint(ndict) = length
            iprec(ndict)  = ip
            dlist(ndict)  = n
            if(num.le.llist) then
              ddict(ndict) = num
            else
              ddict(ndict) = num-llist
            endif
            pdict(ndict) = num
          else
            write(iow,2001) dname
            if(ior.lt.0) then
              write(*,2001) dname
            endif
          endif

c         Memory allocation available - Initialize to zero

          if(adr(n).ne.0) then
            if(ip.eq.1) then
              call pzeroi(mr(np(n)),length)
            else
              call pzero (hr(np(n)),length)
            endif

c         Sufficient memory does not exist - Write ERROR message & STOP

          else
            write(iow,2000) name,length
            if(ior.lt.0) then
              write(*,2000) name,length
            endif
            if(eralloc) then
              setmem = .false.
            else
              call plstop()
            endif

          endif

c       Pointer already exists: Delete if length = 0;

        elseif(length.eq.0) then

          dicloc = mlist(1,n)

          call free(adr(n))

          np(n)      = 0
          mlist(1,n) = 0
          mlist(2,n) = 0

          do i = dicloc,ndict-1
            dict(i)           = dict(i+1)
            ipoint(i)         = ipoint(i+1)
            iprec(i)          = iprec(i+1)
            dlist(i)          = dlist(i+1)
            ddict(i)          = ddict(i+1)
            pdict(i)          = pdict(i+1)
            mlist(1,dlist(i)) = mlist(1,dlist(i)) - 1
          end do ! i

c         Set last entry and reduce entries in dictionary

          dict(ndict)   = '     '
          ipoint(ndict) = 0
          dlist(ndict)  = 0
          ddict(ndict)  = 0
          pdict(ndict)  = 0
          iprec(ndict)  = 0

c         Reset dictionary lengths

          ndict = ndict - 1

c       Pointer already exists: resize in place

        elseif(length.ne.mlist(2,n)) then

          dicloc  = mlist(1,n)

c         Expand array if space available

          ip    =  abs(precis)
          ipa   =  irp(ip,ipr)

          if(length.gt.mlist(2,n)) then
            lensav = mlist(2,n)
          else
            lensav = length
          endif

c         Save current values

          open (unit = iot, file = 'scratch', form = 'unformatted')
          if(ip.eq.1) then
            write(iot) (mr(point),point = np(n),np(n)+lensav-1)
          else
            write(iot) (hr(point),point = np(n),np(n)+lensav-1)
          endif
          call free(adr(n))
          adr(n) = malloc(length*ipa)
          if(adr(n).ne.0) then
            rewind iot
            if(ip.eq.1) then
              np(n) = 1 + (adr(n) - loc(mr(1))) / ipa
              read(iot) (mr(point),point = np(n),np(n)+lensav-1)
              do i = mlist(2,n),length-1
                mr(np(n)+i) = 0
              end do ! i
            else
              np(n) = 1 + (adr(n) - loc(hr(1))) / ipa
              read(iot) (hr(point),point = np(n),np(n)+lensav-1)
              do i = mlist(2,n),length-1
                hr(np(n)+i) = 0.0d0
              end do ! i
            endif
          endif
          close(iot,status = 'delete')

c         Set new pointers

          mlist(2,n)      = length
          ipoint(dicloc)  = length

c         Cannot expand array, not enough space available

          if(adr(n).eq.0) then
            write(iow,2000) name,length
            if(ior.lt.0) then
              write(*,2000) name,length
            endif
            if(eralloc) then
              setmem = .false.
            else
              call plstop()
            endif
          endif

        endif

      else

c       Error indicator

        setmem = .false.
        if(num.le.llist) then
          write(  *,3000) num,dname
          write(iow,3000) num,dname
        else
          write(  *,3001) num-llist,dname
          write(iow,3001) num-llist,dname
        endif

      endif

c     Formats

2000  format(' **ERROR** Insufficient storage to allocate ',a/,
     &              11x,'Required size =',i12/,
     &       '           Check data or choose other option.')

2001  format(' *ERROR* No more room in dictionary for ',a5)

3000  format(' *ERROR* No allocation for array number',i4,' named: ',a)

3001  format(' *ERROR* No allocation for user array number',i4,
     &       ' named: ',a)

      end
