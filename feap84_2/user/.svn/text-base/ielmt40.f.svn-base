c$Id:$
      subroutine ielmt40(d1,d2,u1,u2,x1,x2,t1,t2,ix1,ix2,intnod,s,r,isw)

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Interface element routine

c      Inputs:
c        d1(*),d2(*)   - Material parameters for elements 1 and 2
c        u1(*),u2(*)   - Solution parameters for elements 1 and 2
c        x1(*),x2(*)   - Nodal coordinates for elements 1 and 2
c        t1(*),t2(*)   - T values for elements 1 and 2
c        ix1(*),ix2(*) - Nodal connections for elements 1 and 2
c        intnod(*)     - Interface vertex identifier nodes
c        isw           - Function value to perform

c      Outputs:
c        s(nsts,nsts)  - Element matrix
c        r(ndf,nen,2)  - Element vector [can also be r(nsts)]
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'cdata.h'
      include   'sdata.h'
      include   'eldata.h'
      include   'ieldat.h'
      include   'iofile.h'

      integer    isw, ix1(*),ix2(*), intnod(*)
      real*8     d1(*),u1(ndf,nen,*),x1(ndm,*),t1(*)
      real*8     d2(*),u2(ndf,nen,*),x2(ndm,*),t2(*)
      real*8     s(nsts,nsts),r(ndf,nen,2)

c     Data input mode

      if(isw.eq.1) then

c       DO NOT MAKE CHANGES -- NO INPUTS ALLOWED AS YET

c     Interface check mode

      elseif(isw.eq.2) then

c     Residual and tangent computation

      elseif(isw.eq.3 .or. isw.eq.6) then

c     Output of results/ projections for plots

      elseif(isw.eq.4 .or. isw.eq.8) then

c     Mass/Geometric stiffness computation

      elseif(isw.eq.5) then

      endif

      end
