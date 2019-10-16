c$Id:$
      subroutine iters(bkmax,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'idprtr.h' and change mp(31) to id31         01/04/2009
c       2. Add 'complx.h' and multiply nnm by ipc           24/07/2009
c       3. Force at least 1 word for OINMO                  28/07/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:    Controls solution by iterative schemes. Also sets
c                  storage for compressed arrays for direct solution
c                  by sparse solvers or in blocks with disk stores.

c      Inputs:
c         bkmax  - Maximum size of block for direct solution
c         isw    - Switch: isw =  1 for TANGn
c                          isw =  2 for MASSn
c                          isw =  3 for DAMPn
c                          isw = -1 for USER module

c      Outputs:
c         none   - Outputs stored in pointers.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'allotd.h'
      include  'cdata.h'
      include  'comblk.h'
      include  'compac.h'
      include  'compas.h'
      include  'complx.h'
      include  'debugs.h'
      include  'fdata.h'
      include  'idptr.h'
      include  'iofile.h'
      include  'ndata.h'
      include  'part0.h'
      include  'pglob1.h'
      include  'pointer.h'
      include  'sdata.h'
      include  'ssolve.h'

      logical   setvar,palloc
      character tname*5
      integer   isw,j,len,kp,bkmax, iptc
      real*4    tary(2), etime , tt

      save

c     Compute sparse storage for non-zero matrix

      if(isw.gt.0) then
        setvar = palloc(111,'TEMP1', max(numnp*ndf,neq), 1)
        call pzeroi(mr(np(111)),max(numnp*ndf,neq))
        call elcnt(numel,nen,nen1,mr(id31),mr(np(33)),mr(np(111)),1)
        call sumcnt(mr(np(111)),neq,kp)

        setvar = palloc(112,'TEMP2', kp, 1)
        call pelcon(numel,nen,nen1,mr(np(33)),mr(id31),
     &              mr(np(111)),mr(np(112)),kp,1)
      endif

c     1. TANGENT Formation

      if(isw.eq.1) then

        if(ittyp.eq.-2) then
          kbycol = .false.    ! Store by rows
          kdiag  = .true.
          kall   = .false.
          iptc   =  2*neq + 1
          nspo   =  0
        else
          kbycol = .true.     ! Store by columns
          kdiag  = .false.
          kall   = .false.
          iptc   =  neq + geqnum
        endif

c       Allow for assembly of global equations to tangent matrix

        gceflg = .true.

        setvar = palloc( 93,'OINC ', iptc , 1)
        if(np(94).ne.0) then
          setvar = palloc( 94,'OINO ',0, 1)
        endif
        setvar = palloc(113,'TEMP3',neq+geqnum, 1)
        call comproa(numnp,nen,nen1,ndf,mr(np(33)),mr(id31),
     &              mr(np(111)),mr(np(112)),mr(np(113)),kp,
     &              kbycol,kdiag,kall)
        setvar = palloc(113,'TEMP3',  0, 1)
        setvar = palloc( 94,'OINO ', kp, 1)
        call comprob(numnp,nen,nen1,ndf,mr(np(33)),mr(id31),
     &              mr(np(111)),mr(np(112)),mr(np(94)),mr(np(93)),
     &              kbycol,kdiag,kall)
        kcmplx = kp

        if(debug .and. ndebug.gt.1) then
          call ioprof(mr(np(93)),mr(np(94)),neq,'EQUATION STORAGE',1)
        endif

c       Delete temporary arrays

        setvar = palloc(112,'TEMP2', 0, 1)
        setvar = palloc(111,'TEMP1', 0, 1)

c       Set storage for sparse stiffness array

        write(tname,'(4hTANG,i1)') npart

        if(ittyp.eq.-1) then

          kp  = max(kp+neq,bkmax)
          len = kp

        elseif(ittyp.eq.-2) then

          len = 0
          call sortjc(mr(np(93)),mr(np(94)),neq)

c         Set 'p' and 'ip' to equation numbers

          setvar = palloc( 47, 'PNTER',neq,  1)
          setvar = palloc( 48, 'INVPT',neq,  1)
          do j = 1,neq
            mr(np(47)+j-1) = j
            mr(np(48)+j-1) = j
          end do ! j
          ddom = .true.

        else

          kp  = kp + neq
          len = 0

        endif

        setvar = palloc(npart,tname, kp+len, 2)

        na     = np(npart)
        nnr    = kp
        nau    = na  + neq
        nal    = nau + len
        numcels= 0
        compfl = .true.

c       Block diagonal preconditioner storage

        if(ittyp.eq.2) then
          setvar = palloc( 80,'JPR  ',neq, 1)
          call blkpcg(mr(id31),mr(np(80)),ndf,numnp,neq)
          setvar = palloc( 68,'AUR  ',mr(np(80)+neq-1),2)

c       Set profile preconditioner storage

        elseif(ittyp.gt.2) then
          setvar = palloc( 80,'JPR  ',neq, 1)
          call propcg(mr(np(93)),mr(np(94)),mr(np(80)),neq,ittyp)
          setvar = palloc( 68,'AUR  ',mr(np(80)+neq-1),2)
        endif

c     2. Consistent MASS Formation

      elseif(isw.eq.2) then
        compfl = .true.
        mbycol = .true.       ! Store by columns
        mdiag  = .false.
        mall   = .false.

c       No assembly of global equations to tangent matrix

        gceflg = .false.

c       Set 'p' and 'ip' to equation numbers

        if(np(48).eq.0) then
          setvar = palloc( 47, 'PNTER',neq,  1)
          setvar = palloc( 48, 'INVPT',neq,  1)
          do j = 1,neq
            mr(np(47)+j-1) = j
            mr(np(48)+j-1) = j
          end do ! j
        endif

        setvar =  palloc( 90,'OINMC', neq  , 1)
        if(np(91).ne.0) then
          setvar = palloc( 91,'OINMO',0, 1)
        endif
        setvar = palloc(113,'TEMP3',neq, 1)
        call comproa(numnp,nen,nen1,ndf,mr(np(33)),mr(id31),
     &               mr(np(111)),mr(np(112)),mr(np(113)),kp,
     &               mbycol,mdiag,mall)
        setvar = palloc(113,'TEMP3',  0, 1)
        setvar = palloc( 91,'OINMO', max(1,kp), 1)
        call comprob(numnp,nen,nen1,ndf,mr(np(33)),mr(id31),
     &               mr(np(111)),mr(np(112)),mr(np(91)),mr(np(90)),
     &               mbycol,mdiag,mall)
        mcmplx = kp

c       Delete temporary arrays

        setvar = palloc(112,'TEMP2', 0, 1)
        setvar = palloc(111,'TEMP1', 0, 1)

c       Allocate mass storage

        nnm = (kp  + neq)*ipc
        write(tname,'(4hCMAS,i1)') npart
        setvar = palloc(npart+8,tname, nnm, 2)

c     3. Consistent DAMP Formation

      elseif(isw.eq.3) then
        compfl = .true.
        dbycol = .true.       ! Store by columns
        ddiag  = .false.
        dall   = .false.

c       No assembly of global equations to tangent matrix

        gceflg = .false.

c       Set 'p' and 'ip' to equation numbers

        if(np(48).eq.0) then
          setvar = palloc( 47, 'PNTER',neq,  1)
          setvar = palloc( 48, 'INVPT',neq,  1)
          do j = 1,neq
            mr(np(47)+j-1) = j
            mr(np(48)+j-1) = j
          end do ! j
        endif

        setvar =  palloc(203,'OINDC', neq  , 1)
        if(np(204).ne.0) then
          setvar = palloc(204,'OINDO',0, 1)
        endif
        setvar = palloc(113,'TEMP3',neq, 1)
        call comproa(numnp,nen,nen1,ndf,mr(np(33)),mr(id31),
     &               mr(np(111)),mr(np(112)),mr(np(113)),kp,
     &               dbycol,ddiag,dall)
        setvar = palloc(113,'TEMP3',  0, 1)
        setvar = palloc(204,'OINDO', kp, 1)
        call comprob(numnp,nen,nen1,ndf,mr(np(33)),mr(id31),
     &               mr(np(111)),mr(np(112)),mr(np(204)),mr(np(203)),
     &               dbycol,ddiag,dall)
        dcmplx = kp

c       Delete temporary arrays

        setvar = palloc(112,'TEMP2', 0, 1)
        setvar = palloc(111,'TEMP1', 0, 1)

c       Allocate damping storage

        nnc = kp  + neq
        write(tname,'(4hDAMP,i1)') npart
        setvar = palloc(npart+16,tname, nnc, 2)

      elseif (isw.lt.0) then

        call uiters(kp,isw)
        bkmax  = kp

      endif

c     Output solution properties

      if(pfr) then
        tt = etime(tary)
        write(iow,2000) kp,tary
        if(ior.lt.0) then
          write(*,2000) kp,tary
        endif
      endif

c     Format

2000  format(10x,'Compressed Storage =',i9,20x,'t=',0p,2f9.2)

      end
