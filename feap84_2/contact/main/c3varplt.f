c$Id:$
      subroutine c3varplt (ix1,ch1,ch2,ch3,npair,csw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact 3d GEOmetry PLot

c      Purpose: Plot geometry of 3D contact pair

c      Inputs:
c         x(*)    - Nodal coordinates in deformed position
c         ix1(*)  - Element nodal connection list for surface 1
c         ix2(*)  - Element nodal connection list for surface 2

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_geom.h'
      include  'c_keyh.h'
      include  'c_pair.h'
      include  'cdat1.h'
      include  'pointer.h'
      include  'sdata.h'
      include  'comblk.h'

      integer   ix1(*),npair,csw
      real*8    ch1(*),ch2(*),ch3(*)

      logical   ifplt,setvar,palloc,initfl
      integer   chvec,chvar,nfl

      save

      call cdebug0 ('      c3varplt',-1)

c     Get plot flag and variable

      call setcplv (ifplt,chvec,chvar)

      if (ifplt) then

c       Allocate scratch vectors for variable to plot
c       WARNING check for free vector have to be performed

        setvar = palloc ( 136,'CTEM1',nset,2)

c       extract history variable

        call cextvar (ch1,ch2,ch3,chvec,chvar,hr(np(136)))

c       set flag for set range and profile

        if (csw.eq.308) then
          if (npair.eq.1) then
            initfl = .true.
          else
            initfl = .false.
          endif
          call crprint (hr(np(136)),initfl,nfl)
          setvar = palloc ( 136,'CTEM1',0,2)
        else
          call pltlns(hr(np(53)),ix1,hr(np(136)),nset,
     &                3,1,nope1,neps1,1,-chvar,7,.true.)
          setvar = palloc ( 136,'CTEM1',0,2)
        endif
      endif

      end
