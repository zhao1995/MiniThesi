c$Id:$
      subroutine uiters(nnu,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  User utility for forming sparse arrays.

c     Inputs:
c       isw      Switch parameter (must be a negative number)

c     Outputs:
c       nnu      Length of allocated array: Integer parameter
c                Location arrays: Returned in pointers
c                   np(225) for # terms/column
c                   np(226) for location of terms
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'compac.h'
      include   'idptr.h'
      include   'pglob1.h'
      include   'sdata.h'

      include   'pointer.h'
      include   'comblk.h'

      logical    setvar,palloc
      integer    isw,kp,nnu, i

      save

c     Compute sparse storage for matrix

      if(isw.lt.0) then

        setvar = palloc(151,'USER1', max(numnp*ndf,neq), 1)
        call elcnt(numel,nen,nen1,mr(id31),mr(np(33)),mr(np(151)),1)

c       Check for contact elements

        if(numcels.gt.0 .and. np(168).ne.0) then
          call elcnt(numcels,ncen,ncen1,mr(id31),mr(np(168)),
     &               mr(np(151)),-1)
        endif
        call sumcnt(mr(np(151)),numnp*ndf,kp)

        setvar = palloc(152,'USER2', kp, 1)
        call pelcon(numel,nen,nen1,mr(np(33)),mr(id31),mr(np(151)),
     &              mr(np(152)),kp,1)
        if(numcels.gt.0 .and. np(168).ne.0) then
          call pelcon(numcels,ncen,ncen1,mr(np(168)),mr(id31),
     &                mr(np(151)),mr(np(152)),kp,-1)
        endif

c       Determine sparse matrix structure

        if(np(225).ne.0) then
          setvar = palloc(225,'USOL1',0, 1)
        endif
        setvar = palloc(225,'USOL1', 2*max(ndf*numnp,neq)+1, 1)
        setvar = palloc(153,'USER3',neq, 1)
        call comproa(numnp,nen,nen1,ndf,mr(np(33)),mr(id31),
     &               mr(np(151)),mr(np(152)),mr(np(153)),kp,
     &               ubycol,udiag,uall)
        setvar = palloc(153,'USER3', 0,  1)
        if(np(226).ne.0) then
          setvar = palloc(226,'USOL2',0, 1)
        endif
        setvar = palloc(226,'USOL2', kp, 1)
        call comprob(numnp,nen,nen1,ndf,mr(np(33)),mr(id31),
     &               mr(np(151)),mr(np(152)),mr(np(226)),mr(np(225)),
     &               ubycol,udiag,uall)
        nnu    = kp
        kcmplx = kp

c       Delete temporary arrays

        setvar = palloc(152,'USER2', 0, 1)
        setvar = palloc(151,'USER1', 0, 1)

c       Sort column entries

        call sortjc(mr(np(225)),mr(np(226)),neq)

c       Allocate permutation arrays

        if(isw.lt.-1) then
          setvar = palloc( 47, 'PNTER',neq,  1)
          setvar = palloc( 48, 'INVPT',neq,  1)
          do i = 1,neq
            mr(np(47)+i-1) = i
            mr(np(48)+i-1) = i
          end do ! j
        endif

      else
        write(*,*) '  *ERROR* Negative value required. ISW =',isw
      endif

      end
