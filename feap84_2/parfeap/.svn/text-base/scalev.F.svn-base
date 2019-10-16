c$Id:$
      subroutine scalev(v,pdf,ndm,ndf,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Added petscksp.h                                 01/05/2007
c       2. Change 'include/finclude' to 'finclude'          23/01/2009
c       3. Remove common 'pfeapa' (values in 'setups')      05/02/2009
c       4. Update petsc includes to v3.1                    20/07/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Scale plot vector to have maximum element of +1.0

c      Inputs:
c         v(ndf,*) - Vector of values
c         pdf(*)   - DOF to scale on
c         ndm      - Space dimension of mesh
c         ndf      - DOF's/node (maximum)
c         numnp    - Number of nodes

c      Outputs:
c         v(ndf,*) - Unit vector of values
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

#     include  "finclude/petsc.h"
#     include  "pfeapb.h"
#     include  "pfeapc.h"
#     include  "setups.h"

      integer   i,n,ndm,ndf,numnp, pdf(*)
      real*8    v(ndf,*),vmax, vv, tdatabuf

      integer    uscale_msg, req
      parameter (uscale_msg=14)

c     Petsc values

      PetscErrorCode ierr

      integer        msg_stat(MPI_STATUS_SIZE)

      save

c     Locate maximum

      vv = 0.0d0
      do i = 1,ndm
        if(pdf(i).ge.1.and.pdf(i).le.ndf) then
          do n = 1,numnp
            vv = max(vv,abs(v(pdf(i),n)))
          end do ! n
        endif
      end do ! i
      vmax = vv

c     Get true scaling by comparing values from all processors

      if(pfeap_on) then
c       Sending data asynchronously

        do n = 0,ntasks-1
          if(n.ne.rank) then
            call MPI_Isend( vv, 1, MPI_DOUBLE_PRECISION, n,
     &                      uscale_msg, MPI_COMM_WORLD, req, ierr)
          endif
        end do ! n

c       Receiving data in blocking mode

        do n = 0,ntasks-1
          if(n.ne.rank) then
            call MPI_Recv( tdatabuf, 1, MPI_DOUBLE_PRECISION, n,
     &                     uscale_msg, MPI_COMM_WORLD, msg_stat, ierr)
            vmax = max(vmax,tdatabuf)
          endif
        end do ! n

c       Could fix above to improve efficiency

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

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
        write(*,*) 'Zero length vector in SCALEV'
      endif

      end
