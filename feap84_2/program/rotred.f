c$Id:$
      subroutine rotred(rotyp,id,nty,ndf,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 5/6 Dof constraint check

c      rotyp:>0 = User rotation type: Program urotrd function

c            -1 = 5  dof - Delta Theta - Shells, quaterneon
c            -2 = 6  dof - Exponential map (Simo & Vu-Quoc)
c            -3 = 6  dof - Delta Theta - Cayley, quaterneon, E-M
c            -4 = 6  dof - Delta Theta - Using exponential map.
c            -5 = 6  dof - Delta Theta - Shells inter, Matrices

c      Inputs:
c         rotyp(*)  - Rotation type array for each node
c         nty(*)    - Nodal type
c         ndf       - Number dof/node
c         numnp     - Number of nodes in mesh

c      Outputs:
c         id(ndf,*) - Boundary condition array for rotational dof
c                     not used
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ndf,numnp, n
      integer   rotyp(numnp), id(ndf,numnp), nty(numnp)

      save

      if(ndf.ge.6) then

        do n = 1,numnp
          if(nty(n).ge.0) then
            if(rotyp(n).gt.0) then

c             User check functions

              call urotrd(rotyp(n),id(1,n),ndf)

c           Constrain 6-th dof for type 1

            elseif(rotyp(n).eq.-1) then

              id(6,n) = 1

            end if

          end if
        end do ! n

      end if

      end
