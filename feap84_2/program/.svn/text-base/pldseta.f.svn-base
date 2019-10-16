c$Id:$
      subroutine pldseta(ldtab,f, typ,ndf,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/01/2009
c       1. Increase ldtab to store spin number of displ.    09/03/2009
c          Save spnum for spflg = .true.
c       2. Correct set of 'npt' for load table pointer      24/08/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Move force/displacements to load tables

c      Inputs:
c        ldtab(4,2,*)  - Table with pointers, lengths, proportional no.
c        f(ndf,numnp)  - Expanded load/displacement values
c        typ           - Type: 1 = force; 2 = displacement
c        ndf           - Dof's/node (max)
c        numnp         - Number of nodes in mesh

c      Outputs:
c        ldnod(*)      - List of non-zero nodes
c        ldval(ndf,*)  - List of non-zero values
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'pload1.h'
      include   'pointer.h'
      include   'comblk.h'

      integer    typ,ndf,numnp, ldtab(4,2,*)
      real*8     f(ndf,numnp)

      logical    setvar, palloc
      integer    j
      integer    n, nld, npt, incnld, first, last

c     Determine number of nonzero components

      nld = 0
      do n = 1,numnp
        do j = 1,ndf
          if(f(j,n).ne.0.0d0) then
            nld = nld + 1
            exit
          endif
        end do ! j
      end do ! n

c     Allocate new space for node numbers: Compute pointer for start

      npt = 0
      do n = 1,ldtot
        if(n.eq.1) then
          npt = max(ldtab(1,1,n  ) + ldtab(2,1,n  ),
     &              ldtab(1,2,n  ) + ldtab(2,2,n  ), npt)
        else
          npt = max(ldtab(1,1,n-1) + ldtab(2,1,n-1),
     &              ldtab(1,2,n-1) + ldtab(2,2,n-1),
     &              ldtab(1,1,n  ) + ldtab(2,1,n  ),
     &              ldtab(1,2,n  ) + ldtab(2,2,n  ), npt)
        endif
      end do ! n

      if(nld.gt.ldtab(2,typ,ldnum)) then
        incnld = nld - ldtab(2,typ,ldnum)

c       Adjust table pointers values

        last = 0
        if(typ.eq.1 .and. ldtab(1,2,ldnum).gt.0) then
          last             = max(last,ldtab(1,2,ldnum)+ldtab(2,2,ldnum))
          ldtab(1,2,ldnum) = ldtab(1,2,ldnum) + incnld
        endif
        do n = ldnum+1,ldtot
          do j = 1,2
            if(ldtab(1,j,n).gt.0) then
              last         = max(last,ldtab(1,j,n)+ldtab(2,j,n))
              ldtab(1,j,n) = ldtab(1,j,n) + incnld
            endif
          end do ! j
        end do ! n

c       Increase storage for arrays

        first  = npt  + nld
        last   = max(first,last + incnld)
        setvar = palloc(266,'LDNOD',max(1,last)    ,1)
        setvar = palloc(267,'LDVAL',max(1,last*ndf),2)

c       Move data down in arrays

        do n = last-1,first,-1
          mr(np(266)+n)        = mr(np(266)+n-incnld)
          mr(np(266)+n-incnld) = 0
        end do ! n
        last   = ndf*last
        first  = ndf*first
        incnld = ndf*incnld
        do n = last-1,first,-1
          hr(np(267)+n)        = hr(np(267)+n-incnld)
          hr(np(267)+n-incnld) = 0.0d0
        end do ! n

c     New array

      else

        last = npt + nld
        setvar = palloc(266,'LDNOD',max(1,last)    ,1)
        setvar = palloc(267,'LDVAL',max(1,last*ndf),2)

      endif

c     Set new load table and array pointers

      call pldtabl(ldtab,ldnum,npt,typ, 1)
      call pldtabl(ldtab,ldnum,nld,typ, 2)

      if(spflg) then
        ldtab(4,2,ldnum) = spnum
      endif

c     Store values in arrays

      call pldsetb(f,ndf,numnp, mr(np(266)+npt), hr(np(267)+npt*ndf))

      end

      subroutine pldsetb(f,ndf,numnp, ldnod, ldval)

      implicit   none

      integer    ndf,numnp, ldnod(*)
      real*8     f(ndf,numnp), ldval(ndf,*)

      logical    flag
      integer    j,n,nn

c     Search for non-zero components

      nn = 0
      do n = 1,numnp
        flag = .false.
        do j = 1,ndf
          if(f(j,n).ne.0.0d0) then
            nn   = nn + 1
            flag = .true.
            exit
          endif
        end do ! j

c       Store values if non-zero component exists

        if(flag) then
          ldnod(nn) = n
          do j = 1,ndf
            ldval(j,nn) = f(j,n)
          end do ! j
        end if
      end do ! n

      end
