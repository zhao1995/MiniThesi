c$Id:$
      subroutine prterr()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output error indicator information for accuracy
c               assessments on solution.

c      Inputs:
c         none

c      Outputs:
c         none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'errind.h'
      include  'adapt2.h'

      real*8    erind, enerr, errel

      save

c     Output error indicator values

      if(eproj.ne.0.0d0) then
        erind  = efem/eproj  - 1.0d0
      else
        erind = 0.0d0
      endif

      if(eener+eenere.ne.0.0d0) then
        enerr  = eenere/(eener + eenere)
      else
        enerr  = 0.0d0
      endif
      errel  = sqrt(enerr)

      write(iow,2000) efem,eproj,enerr, errel, erind,eerror,
     +                eenere, eener
      if(ior.lt.0) write(*,2000) efem,eproj,enerr, errel, erind,eerror,
     +                eenere, eener

c     Format

2000  format(/'   Finite Element Stress Measure =',1p,e24.16/
     +        '   Projected Stress Measure      =',1p,e24.16/
     +        '   Direct 1  Error  Measure      =',1p,e24.16/
     +        '   Relative  Error  Measure      =',1p,e24.16/
     +        '   Error Indicator               =',1p,e24.16/
     +        '   Direct Error Indicator        =',1p,e24.16/
     +        '   Projected Stress Measure   *  =',1p,e24.16/
     +        '   Energy    Error  Measure   *  =',1p,e24.16/)

      end
