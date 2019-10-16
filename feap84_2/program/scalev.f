c$Id:$
      subroutine scalev(v,pdf,ndm,ndf,numnp,vflag)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add vflag and tmax to check for eigenvalue scale 09/10/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Scale vector to have maximum element of +1.0

c      Inputs:
c         v(ndf,*) - Vector of values
c         pdf(*)   - DOF to scale on
c         ndm      - Space dimension of mesh
c         ndf      - DOF's/node (maximum)
c         numnp    - Number of nodes
c         vflag

c      Outputs:
c         v(ndf,*) - Unit vector of values
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   vflag
      integer   i,n,ndm,ndf,numnp, pdf(*)
      real*8    v(ndf,*),vmax,tmax

      save

c     Check for maximum in vector

      if(vflag) then
        tmax = 0.0d0
        do n = 1,numnp
          do i = 1,ndf
            tmax = max(tmax,abs(v(i,n)))
          end do ! i
        end do ! i
      endif

c     Locate maximum in dimension of plot

      vmax = 0.0d0
      do i = 1,ndm
        if(pdf(i).ge.1.and.pdf(i).le.ndf) then
          do n = 1,numnp
            vmax = max(vmax,abs(v(pdf(i),n)))
          end do ! n
        endif
      end do ! i

c     Check for zero scaling for eigenvector

      if(vflag .and. vmax.lt.1.d0*tmax) then
        vmax = 1.d0
      endif

c     Perform scaling

      if(vmax.gt.0.0d0) then
        vmax = 1.d0/vmax
        do n = 1,numnp
          do i = 1,ndf
            v(i,n) = v(i,n)*vmax
          end do ! i
        end do ! n
      else
        write(*,*) ' ** WARNING ** Zero length vector in SCALEV'
      endif

      end
