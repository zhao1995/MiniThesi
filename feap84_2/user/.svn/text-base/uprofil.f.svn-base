c$Id:$
      subroutine uprofil(jp,idl,id,ix,iop,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute user additions to profile of global arrays

c      Inputs:
c        jp(*)  - Pointer array to row/column ends of profile (iop = 2)
c        idl(*) - Local temporary storage vector
c        id(*)  - Equation numbers for degree of freedoms     (iop = 1)
c        ix(*)  - Global node numbers on elements
c        iop    - Switch to control operation
c                  = 1 to set up equation numbers of dof's
c                  = 2 to compute the column/row lengths and profile.
c        prt    - Flag, print solution properties if true

c      Outputs:
c        jp(*)  - Modified pointer array to row/column ends of profile
c                 (iop = 2)
c        id(*)  - Modified equation numbers for degree of freedoms
c                 (iop = 1)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      logical    prt
      integer    jp(*),idl(*), id(*), ix(*), iop

c     Modify equation numbers

      if(iop.eq.1) then

c     Adjust column heights

      elseif(iop.eq.2) then

      endif

      end
