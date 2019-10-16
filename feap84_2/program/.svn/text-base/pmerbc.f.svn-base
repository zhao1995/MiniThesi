c$Id:$
      subroutine pmerbc(nopart,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Merge boundary condition, forces and displacements from
c               ties

c      Inputs:
c        nopart -  Partition flag
c        prt    -  Print results if true

c      Outputs:
c        none
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'cdata.h'
      include   'part0.h'
      include   'part3.h'
      include   'sdata.h'

      include   'pointer.h'
      include   'comblk.h'

      logical    nopart,prt, setvar,palloc
      integer    j

      save

c     Reset node numbers on elements to reflect any tied conditions

      call ptielm(mr(np(33)),mr(np(79)),nen,nen1,numel)

c     Output list of tied nodes

      setvar = palloc(111,'TEMP1',numnp, 1)
      call poutie(mr(np(111)),mr(np(33)),mr(np(190)),nen,nen1,
     &            numnp,numel,prt)
      setvar = palloc(111,'TEMP1',0, 1)

c     Set partion data

      if(nopart) then
        nopart = .false.
        npart  = 1
        call partpt(npart,tflp(npart),.false.)
      endif

c     Reset forces for tied nodes -- set to eliminate unused equations
     
      do j = 1,4
        if(.not.tflp(j)) then
          call partpt(j,tflp(j),.false.)
          call tiefor(mr(np(31)+nneq),mr(np(100)),hr(np(27)),mr(np(79)),
     &                ndf,numnp)
        endif
      end do ! j

c     Check contact surfaces for eliminated tied nodes

c     call contact (312)

      end
