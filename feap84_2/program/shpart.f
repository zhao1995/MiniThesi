c$Id:$
      subroutine shpart()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Show partition data

c     Inputs: None

c     Outputs: None
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include  'cdata.h'
      include  'compas.h'
      include  'ddata.h'
      include  'fdata.h'
      include  'part0.h'

      save

      write(*,2000) npart,ittyp,neq,fl,noi,nrk,nrc,nrm

2000  format(/5x,'ACTIVE PARTITION = ',i3/
     &       /8x,'Solver type         = ',i8/
     &        8x,'Number of equations = ',i8/
     &        8x,'Solution flags      = ',12l2/
     &        8x,'Time integration    = ',i8/
     &       12x,'nrk, nrc, nrm   = ',3i5/)

      end
