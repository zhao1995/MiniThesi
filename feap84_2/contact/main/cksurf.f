c$Id:$
      subroutine cksurf(cs0,ics,ipos)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Coded/modified by:                           Date:          rel.:
c              Robert L. Taylor            June 15, 1997            1.0

c     Acronym: ChecK SURFace numbers after ties

c     Purpose: Place correct node number on contact surfaces to reflect
c              merges by the 'tie' command.

c     Inputs:
c       cs0(*,*)   - Offset and dimensioning information
c       ics(*)     - List of original nodes for each surface.
c       ipos(*)    - Renumber information from ties

c     Outputs:
c       ics(*)     - List of renumbered nodes for each surface.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_contac.h'
      include  'c_comnd.h'

      integer   n, ofs,dnope,ntype,nope,neps, nepsr
      integer   ics(*),ipos(*)
      real*8    cs0(nr0,n0c1:nc01,*)

      save

c     Check contact surface nodes for nodes eliminated by ties

      do n = 1,numcs

        ofs   = nint(cs0(2,-1,n))
        neps  = nint(cs0(3,-1,n))
        dnope = nint(cs0(4,-1,n))
        ntype = nint(cs0(1, 0,n))
        nope  = nint(cs0(2, 0,n))

        call ctiend(ics(ofs),neps,dnope,nope,ipos)

        if    (ntype.eq.1) then              ! 2-d facet surface
          call csurflr(ics(ofs),dnope,neps)
        elseif(ntype.eq.5) then              ! Point surface
          call psurflr(ics(ofs),neps,nepsr)
          if(nepsr.lt.neps) then
            cs0(3,-1,n) = nepsr
          endif
        endif

      end do ! n

      end
