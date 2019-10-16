c$Id:$
      subroutine pelnum(tx,iel,errck)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add St.Venant torsion                            19/04/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Get material set numbers for FEAP elements

c              Current FEAP Element Types
c                 Name     |     iel
c              ------------+-------------
c               Solid      |     -1
c               Truss      |     -2
c               Frame      |     -3
c               Plate      |     -4
c               Shell      |     -5
c               Membrane   |     -6
c               Gap        |     -7
c               Thermal    |     -8
c               Convection |     -9
c               Point      |     -10
c               Pressure   |     -11
c               Torsion    |     -12
c               Coupled:   |
c                  THMEch  !     -13

c      Inputs:
c         tx(2)  - Name of element type requested

c      Outputs:
c         iel    - Element type for request
c         errck  - Flag, true if request found
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      logical   pcomp, errck
      character tx(2)*15
      integer   iel

      save

      errck = .false.

c     Check based on coordinate dimension

      if    (pcomp(tx(1),'soli',4)) then
        iel = -1
      elseif(pcomp(tx(1),'trus',4)) then
        iel = -2
      elseif(pcomp(tx(1),'fram',4)) then
        iel = -3
      elseif(pcomp(tx(1),'plat',4)) then
        iel = -4
      elseif(pcomp(tx(1),'shel',4)) then
        iel = -5
      elseif(pcomp(tx(1),'memb',4)) then
        iel = -6
      elseif(pcomp(tx(1),'gap',3)) then
        iel = -7
      elseif(pcomp(tx(1),'ther',4)) then
        iel = -8
      elseif(pcomp(tx(1),'conv',4)) then
        iel = -9
      elseif(pcomp(tx(1),'poin',4)) then
        iel = -10
      elseif(pcomp(tx(1),'pres',4)) then
        iel = -11
      elseif(pcomp(tx(1),'tors',4)) then
        iel = -12
      elseif(pcomp(tx(1),'coup',4)) then
        if(pcomp(tx(2),'thme',4)) then
          iel = -13
        else
          write(  *,*) ' Unknown coupled material type:',tx(2)
          write(iow,*) ' Unknown coupled material type:',tx(2)
          call plstop()
        endif
      else
        errck = .true.
      endif

      end
