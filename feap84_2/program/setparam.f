c$Id:$
      subroutine setparam(par, value, prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add fixed length character 'name'                04/01/2008
c       2. Set pconset to 'redo'                            15/06/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Command language instruction subprogram: Part 6
c               Set parameter 'par' to 'value'.

c      Inputs:
c        par     - Parameter name to set
c        value   - Value to assign
c        prt     - Echo value set if true

c      Outputs:
c        none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'comfil.h'

      character  par*(*),name*15
      logical    prt, redo, pconset
      real*8     value

c     Put data into 'record'

      name   = par
      record = ' '
      record( 1:14) = name(1:14)
      record(15:15) = '='
      write(record(16:30),'(1p,1e15.7)') value

      redo = pconset(prt)

      end
