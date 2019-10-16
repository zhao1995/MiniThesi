c$Id:$
      subroutine prax(nty,vex,st,numnp,ndm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute principal vectors for plotting

c      Programmer: Sanjay Govindjee

c      Inputs:
c         nty(*)       - Nodal type
c         st(numnp,*)  - Projected nodal stress values
c         numnp        - Number of nodes in mesh
c         ndm          - Spatial dimension of mesh

c      Outputs:
c         vex(numnp,*) - Vectors for principal directions at nodes
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'strnum.h'

      integer   numnp,ndm, ii, rot, nty(*)
      real*8    vex(numnp,*),st(numnp,*),v(3,3),d(3),rdn

      save

      data      rdn/0.017453292519943d0/

      if(istv.gt.0) then

        do ii = 1,numnp
          if(nty(ii).ge.0) then

c           Three-dimensional

            if(ndm.eq.3 .or. istp.eq.8) then
              v(1,1)   = st(ii,1)
              v(2,2)   = st(ii,2)
              v(3,3)   = st(ii,3)
              v(1,2)   = st(ii,4)
              v(2,3)   = st(ii,5)
              v(1,3)   = st(ii,6)
              call eig3(v,d,rot)
              vex(ii,1) = v(1,1)
              vex(ii,2) = v(2,1)
              vex(ii,3) = v(3,1)
              vex(ii,4) = v(1,2)
              vex(ii,5) = v(2,2)
              vex(ii,6) = v(3,2)
              vex(ii,7) = v(1,3)
              vex(ii,8) = v(2,3)
              vex(ii,9) = v(3,3)
              vex(ii,10) = d(1)
              vex(ii,11) = d(2)
              vex(ii,12) = d(3)

c           Two-dimensional

            elseif(ndm.eq.2) then

              v(1,1)   = st(ii,1)
              v(2,1)   = st(ii,2)
              v(3,1)   = st(ii,3)
              v(1,2)   = st(ii,4)

              call pstr2d(v, d)

              d(3) = d(3)*rdn

              vex(ii,1) =  cos(d(3))
              vex(ii,2) =  sin(d(3))
              vex(ii,3) = -vex(ii,2)
              vex(ii,4) =  vex(ii,1)
              vex(ii,5) =  d(1)
              vex(ii,6) =  d(2)
            endif
          endif

        end do ! ii

      endif

      end
