c$Id:$
      subroutine prfrst(id,idl,ndf,mm,nad)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute column height for equation

c      Inputs:
c         id(*)    - Equation numbers for node
c         ndf      - Number dof/node

c      Outputs:
c         idl(*)   - List of terms connected to this equation
c         mm       - Minimum entry in idl
c         nad      - Number of entries in idl
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'part0.h'

      integer   ndf,mm,nad, j,jj
      integer   id(ndf),idl(*)

      save

c     Compute column height for equation

      do j = 1,ndf
        if(ndfp(j).eq.npart) then
          jj = id(j)
          if(jj.gt.0) then
            if(mm.eq.0) mm = jj
            mm = min(mm,jj)
            nad = nad + 1
            idl(nad) = jj
          endif
        endif
      end do ! j

      end
