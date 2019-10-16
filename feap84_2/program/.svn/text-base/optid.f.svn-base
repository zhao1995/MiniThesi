c$Id:$
      subroutine optid()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Optimize profile numbering using algorithm of S. Sloan

c      Inputs:
c        Come from pointer arrays and common blocks.

c      Outputs:
c        Optimized list through pointer array 'mr(np(89))'.

c      Arrays used from FEAP
c         mr(np(21)) - JP: Column pointers for profile
c         mr(id31)   - ID: Equation number list
c         mr(np(33)) - IX: Element connection list
c         mr(np(34)) - LD: Local equation list
c         mr(np(89)) - NREN: Renumber list for nodes
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'compac.h'
      include  'complx.h'
      include  'fdata.h'
      include  'idptr.h'
      include  'region.h'
      include  'sdata.h'
      include  'part0.h'
      include  'part1.h'
      include  'part3.h'
      include  'rjoint.h'
      include  'iofile.h'
      include  'umac1.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   setvar, palloc
      integer   e2, oldpro, nesz, newpro, kp

      save

c     Wilson/Hoit Algorithm

      if(opthoit) then

c       Set initial numbering in renumber vector

        do kp = 0,numnp-1
          mr(np(89)+kp) = kp+1
        end do ! kp

c       Compute sparse structure for matrix

        setvar = palloc(111,'TEMP1', numel+numcels, 1)   ! NNAC
        call pzeroi(mr(np(111)),numel+numcels)

        setvar = palloc(112,'TEMP2', numnp+1, 1)         ! IC
        call optic(numnp,numel,numcels,nen,nen1,ncen,ncen1,
     &             mr(np(33)),mr(np(168)),mr(np(112)), kp)

        setvar = palloc(113,'TEMP3', kp, 1)              ! IP
        setvar = palloc(114,'TEMP4', numel+numcels, 1)   ! NNEL
        call pzeroi(mr(np(113)),kp)
        call opcon(numel,numcels,nen,nen1,ncen,ncen1,mr(np(33)),
     &              mr(np(168)),mr(np(112)),mr(np(113)),mr(np(114)) )

        setvar = palloc(115,'TEMP5', numel+numcels, 1) ! NNRM
        setvar = palloc(116,'TEMP6', numnp, 1)         ! NNID

        call optibc(mr(id31),mr(np(116)),ndf,numnp)

        setvar = palloc(117,'TEMP7', nummat, 1)        ! IMAT
        setvar = palloc(118,'TEMP8', nummat, 1)        ! IMRM
        call optrm(mr(np(33)),mr(np(168)),mr(np(112)),mr(np(113)),
     &             nen,nen1,ncen1,numel,numcels,numnp,nummat,
     &             mr(np(114)),mr(np(116)),mr(np(111)),mr(np(115)),
     &             mr(np(117)),mr(np(118)) )

        setvar = palloc(118,'TEMP8', 0,1)
        setvar = palloc(117,'TEMP7', 0,1)
        setvar = palloc(116,'TEMP6', 0,1)
        setvar = palloc(115,'TEMP5', 0,1)
        setvar = palloc(114,'TEMP4', 0,1)
        setvar = palloc(113,'TEMP3', 0,1)
        setvar = palloc(112,'TEMP2', 0,1)

c       Optimize profile using Wilson/Hoit method

        setvar = palloc(112,'TEMP2',numnp,1)                    ! ND
        setvar = palloc(113,'TEMP3',numnp,1)                    ! LN
        setvar = palloc(114,'TEMP4',numnp,1)                    ! NWD
        setvar = palloc(115,'TEMP5',max(numnp,numel+numcels),1) ! MSUM
        setvar = palloc(116,'TEMP6',numnp,1)                    ! NNID
        nesz   = numel*nen*nen
        setvar = palloc(117,'TEMP7',max(numnp,numel+numcels,nesz),1) !NE

        call optibc(mr(id31),mr(np(116)),ndf,numnp)

        call opnum(mr(np(33)),mr(np(168)),mr(np(112)),mr(np(113)),
     &             mr(np(117)),mr(np(114)),mr(np(115)),mr(np(89)),
     &             mr(np(111)),mr(np(116)),
     &             numnp,numel,numcels,nen,nen1,ncen,ncen1,.false.)

        setvar = palloc(117,'TEMP7',0,1)
        setvar = palloc(116,'TEMP6',0,1)
        setvar = palloc(115,'TEMP5',0,1)
        setvar = palloc(114,'TEMP4',0,1)
        setvar = palloc(113,'TEMP3',0,1)
        setvar = palloc(112,'TEMP2',0,1)
        setvar = palloc(111,'TEMP1',0,1)

c     Sloan Algorithm

      else

c       Allocate and initialize arrays

        nesz   = (numel+numcels)*nen*nen
        setvar = palloc( 232, 'USOL8', 4*numnp+1, 1 )
        setvar = palloc( 233, 'USOL9', numnp+1  , 1 )
        setvar = palloc( 234, 'USOL0', nesz     , 1 )
        call pzeroi(mr(np(232)),4*numnp+1)
        call pzeroi(mr(np(233)),  numnp+1)
        call pzeroi(mr(np(234)),  nesz   )

c       Perform optimization

        call sgraph(numnp,numel,nen,nen1,mr(np(33)), nesz, mr(np(234)),
     &              mr(np(233)), e2 )
        call slabel(numnp,e2,mr(np(234)),mr(np(233)),mr(np(232)),
     &              mr(np(232)+numnp),oldpro,newpro)
        call snren (numnp,mr(np(232)),mr(np(89)))

c       Delete scratch arrays

        setvar = palloc( 234, 'USOL0', 0, 1 )
        setvar = palloc( 233, 'USOL9', 0, 1 )
        setvar = palloc( 232, 'USOL8', 0, 1 )

      endif

      end
