c$Id:$
      subroutine pmeshmem()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    08/08/2011
c       1. Increase element storage for ix array            01/07/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Set memory allocation for mesh storage

c     Inputs:
c       Through common blocks

c     Outputs:
c       np(*)    - Pointer values to mesh data
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'complx.h'
      include   'contrl.h'
      include   'gltran.h'
      include   'idptr.h'
      include   'part0.h'
      include   'pload1.h'
      include   'pointer.h'
      include   'sdata.h'
      include   'setups.h'

      logical    setvar, palloc
      integer    i, l1,l2,l3,l4,l5,l6

c     Set pointers for allocation of mesh arrays

      nen1      = nen + 11
      nie       = 13
      nst       = max(nen*ndf + nad,1)
      nneq      = ndf*numnp
      do i = 1,ndf
        ndfl(i) = 0
        ndfo(i) = 0
        ndog(i) = 10
      end do ! i

c     Set pointers for allocation of mesh arrays

      if(cplxfl) then
        ipc = 2
      else
        ipc = 1
      end if

c     Allocate size for arrays for mesh and solution vecors

      l1   = ndm*numnp
      l2   = max(ndf*numnp,1)
      l3   = max(nen+1,7*nst,21)
      l4   = numnp*max(ndf,ndm)*ipc
      l5   = ndf*nen
      l6   = max(1,numel)

c     Allocate and zero arrays

      setvar = palloc( 26,'DR   ',l4          ,  2)
      setvar = palloc( 34,'LD   ',l3          ,  1)
      setvar = palloc( 35,'P    ',nst*3       ,  2)
      setvar = palloc( 36,'S    ',nst*nst*2   ,  2)
      setvar = palloc( 39,'TL   ',nen         ,  2)
      setvar = palloc( 41,'UL   ',nst*14      ,  2)
      setvar = palloc( 44,'XL   ',max(4,nen)*3,  2)
      setvar = palloc( 25,'D    ',nummat*ndd  ,  2)
      setvar = palloc( 32,'IE   ',nummat*nie  ,  1)
      setvar = palloc(240,'IEDOF',nummat*l5   ,  1)
      setvar = palloc( 31,'ID   ',l2*2        ,  1)
      setvar = palloc( 33,'IX   ',nen1*l6     ,  1)
      setvar = palloc(190,'NDTYP',numnp       ,  1)
      setvar = palloc(100,'RIXT ',numnp       ,  1)
      setvar = palloc(181,'RBEN ',l6          ,  1)
      setvar = palloc( 43,'X    ',l1          ,  2)
      setvar = palloc( 45,'ANG  ',numnp       ,  2)
      setvar = palloc( 46,'ANGL ',nen         ,  2)
      setvar = palloc( 27,'F    ',2*l2        ,  2)
      setvar = palloc( 28,'F0   ',4*l2        ,  2)
      setvar = palloc( 29,'FPRO ',2*l2        ,  1)
      setvar = palloc( 30,'FTN  ',4*l2        ,  2)
      setvar = palloc( 38,'T    ',numnp       ,  2)
      setvar = palloc( 40,'U    ',4*l2*ipc    ,  2)
      setvar = palloc( 89,'NREN ',numnp*2     ,  1)
      if(ldtot.gt.0) then
        setvar = palloc(265,'LDTAB',ldtot*12  ,  1)
      endif

c     Set ID address pointers

      id31    = np(31)
      idpt(1) = np(31)
      npid    = np(31)       ! ID
      npix    = np(33)       ! IX
      npuu    = np(40)       ! U
      npxx    = np(43)       ! X
      nprn    = np(89)       ! NREN
      npty    = np(190)      ! NDTYP

      end
