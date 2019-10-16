c$Id:$
      logical function setmem(list,mlist,rlist,
     &                        num,name,length,precis)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2013: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Introduce 'longv' to pass to the functions       10/06/2007
c       2. Compute allocation length in longest possible    27/04/2011
c          size
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
c         np(num)    - Pointer to first word of array in blank common
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'allotd.h'
      include  'allotn.h'
      include  'cdata.h'
      include  'iofile.h'
      include  'errchk.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'p_point.h'

      logical   pcomp
      integer   list,num,length,precis,ip,ipa, n,i, longv
      integer   dicloc, mlist(2,list), irp(2,2)
      character name*(*),dname*5, rlist(list)*5

      save

      data      irp / 3*1, 2 /

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
          ip    =  abs(precis)
          ipa   =  irp(ip,ipr)

c         Call allocate for requested "length" expressed in bytes
c         Compute length in integer*X recast to integer*4 for malloc
c         X is either integer4 or integer8

          point = length
          longv = int((point*ipa + mod(point*ipa,ipr))/ipr)
          call fallocfn(np(n),longv,ip,ipr)

c         Memory allocation available - Initialize to zero

          if(np(n).ne.0) then

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
            if(precis.eq.1) then
              call pzeroi(mr(np(n)),length)
            elseif(precis.eq.2) then
              call pzero (hr(np(n)),length)
            endif

c         Sufficient memory does not exist - Write ERROR message & STOP

          else
            write(iow,2000) dname,length*ipa
            if(ior.lt.0) then
              write(*,2000) dname,length*ipa
            endif
            if(eralloc) then
              setmem = .false.
            else
              call plstop()
            endif

          endif

c       Pointer already exists: Delete if length = 0;

        elseif(length.eq.0) then

          dicloc     =  mlist(1,n)

          call ffreefn(np(n),precis,ipr)

          np(n)      = 0
          mlist(1,n) = 0
          mlist(2,n) = 0

          do i = dicloc,ndict-1
            dict(i)     = dict(i+1)
            ipoint(i)   = ipoint(i+1)
            iprec(i)    = iprec(i+1)
            dlist(i)    = dlist(i+1)
            ddict(i)    = ddict(i+1)
            pdict(i)    = pdict(i+1)
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

c         Compute length in integer*X recast to integer*4 for malloc
c         X is either integer4 or integer8

          point = length
          longv = int((point*ipa + mod(point*ipa,ipr))/ipr)
          call freallocfn(np(n),longv,ip,ipr)

          if(np(n).ne.0) then

            if(ip.eq.1) then
              do i = mlist(2,n),length-1
                mr(np(n)+i) = 0
              end do ! i
            else
              do i = mlist(2,n),length-1
                hr(np(n)+i) = 0.0d0
              end do ! i
            endif
            mlist(2,n)     = length
            ipoint(dicloc) = length

c         Cannot expand array, not enough space available

          else
            write(iow,2000) dname,length*ipa
            if(ior.lt.0) then
              write(*,2000) dname,length*ipa
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

2000  format(' *ERROR* Insufficient memory for array ',a5,
     &       '. Need',i12,' integer words.')

2001  format(' *ERROR* No more room in dictionary for ',a5)

3000  format(' *ERROR* No allocation for array number',i4,' named: ',a)

3001  format(' *ERROR* No allocation for user array number',i4,
     &       ' named: ',a)

      end
