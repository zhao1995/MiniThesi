c$Id:$
      subroutine pecmes()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add check on 'curv' inputs                       13/11/2008
c       2. Add 'espi'n option                               03/01/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Edge and coordinate values

c      Inputs:
c        none

c      Outputs:
c        none
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'cdata.h'
      include   'corset.h'
      include   'edgdat.h'
      include   'sdata.h'

      logical    setvar, palloc
      integer    nspin

      save

c     Set edge boundary codes, forces, displacements, and angles

      if(eanfl.or.ebcfl.or.edifl.or.efcfl.or.eprfl  .or.
     &   ebsfl.or.curfl) then
        call pedgin()
        eanfl = .false.
        ebcfl = .false.
        edifl = .false.
        efcfl = .false.
        eprfl = .false.
        ebsfl = .false.
        curfl = .false.
      endif

c     Set cordinate angles, boundary codes, forces, displacements,
c         proportional load types and surface loads

      if(boufl .or. surfl .or. angfl .or. disfl .or. cprfl .or.
     &   forfl .or. damfl .or. masfl .or. stifl .or. basfl .or.
     &   lfrfl                                            ) then
        call ploadc()
        boufl = .false.
        surfl = .false.
        angfl = .false.
        disfl = .false.
        cprfl = .false.
        forfl = .false.
        lfrfl = .false.
        damfl = .false.
        masfl = .false.
        stifl = .false.
        basfl = .false.
      endif

c     Set edge boundary codes, forces, displacements, and angles

      if(espfl) then
        setvar = palloc(266,'LDNOD', numnp, 1)
        call pedgex(nspin)
        setvar = palloc(266,'LDNOD', nspin    , 1)
        setvar = palloc(267,'LDVAL', nspin*ndf, 2)
        espfl = .false.
      endif

c     Set body forces and reactions

      if(reafl .or. intfl ) then
        call pintec()
        reafl = .false.
        intfl = .false.
      endif

      end
