c$Id:$
      subroutine psetb(b,getp,getv,senp,senv,id, ndf, rdatabuf,sdatabuf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Revise send/receive data add barrier             24/11/2006
c       2. Correct send/receive and add error messages      16/03/2007
c       3. Change 'include/finclude' to 'finclude'          23/01/2009
c       4. Remove common 'pfeapa' (values in 'setups')      05/02/2009
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Transfer PETSC vector to local arrays including ghost
c               node data via MPI messaging

c     Inputs:
c        getp(ntasks)    - Pointer array for getting ghost data
c        getv(sum(getp)) - Local node numbers for getting ghost data
c        senp(ntasks)    - Pointer array for sending ghost data
c        senv(sum(senp)) - Local node numbers to send out as ghost data
c        id(ndf,numnp)   - Local equation numbers
c        ndf             - dofs per node
c        rdatabuf(*)     - receive communication array
c        sdatabuf(*)     - send    communication array

c     Outputs:

c        b(neq) - Local solution vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "finclude/petsc.h"

#     include   "cdata.h"
#     include   "iofile.h"
#     include   "pfeapb.h"
#     include   "setups.h"

      integer    ndf
      integer    i, j, k, lnodei, eqi, soff, rbuf,sbuf,tbuf, idesp
      integer    getp(*),getv(*),senp(*),senv(*)
      integer    id(ndf,*)
      real*8     b(*), rdatabuf(*), sdatabuf(*)

      integer    usolve_msg, req
      parameter (usolve_msg=10)

c     Petsc values

      PetscErrorCode ierr

      integer        msg_stat(MPI_STATUS_SIZE)

c     Sending Data Asynchronously

      soff  = 0
      idesp = 0
      do i = 1, ntasks

        if (senp(i) .gt. 0) then
          sbuf   = soff
          do j = 1, senp(i)
            lnodei = senv(j + idesp)
            do k = 1, ndf
              eqi = id(k,lnodei)
              if (eqi .gt. 0) then
                sbuf           = sbuf + 1
                sdatabuf(sbuf) = b(eqi)
c               write(iow,*) 'O[',rank+1,'->',i,']:',
c    &                        sbuf,sdatabuf(sbuf)
c               write(  *,*) 'O[',rank+1,'->',i,']:',
c    &                        sbuf,sdatabuf(sbuf)
              endif
            end do ! k
          end do ! j
          idesp = idesp + senp(i)
          sbuf  = sbuf  - soff
c         write(iow,*) 'Bufflen[',rank+1,'->',i,']: ',sbuf
c         write(  *,*) 'Bufflen[',rank+1,'->',i,']: ',sbuf
          call MPI_Isend( sdatabuf(soff+1), sbuf, MPI_DOUBLE_PRECISION,
     &                    i-1, usolve_msg, MPI_COMM_WORLD, req, ierr )

c         Report send error

          if(ierr.ne.0) then
            write(iow,*) ' -->> Send Error[',rank+1,'->',i,']',ierr
            write(  *,*) ' -->> Send Error[',rank+1,'->',i,']',ierr
          endif

          soff = soff + sbuf
        endif

      end do ! i

c     Receiving Data in blocking mode

      idesp = 0
      do i = 1, ntasks
        if (getp(i) .gt. 0) then

c         Determine receive length

          tbuf = getp(i)*ndf
c         write(iow,*) 'Bufflen[',i,'->',rank+1,']: ',tbuf
c         write(  *,*) 'Bufflen[',i,'->',rank+1,']: ',tbuf
          call MPI_Recv( rdatabuf, tbuf, MPI_DOUBLE_PRECISION, i-1,
     &                   usolve_msg, MPI_COMM_WORLD,msg_stat, ierr)
          if(ierr.ne.0) then
             write(iow,*) 'Recv Error[',i,'->',rank+1,']',ierr
             write(  *,*) 'Recv Error[',i,'->',rank+1,']',ierr
          endif
          rbuf = 0
          do j = 1, getp(i)
            lnodei = getv(j + idesp)
            do k = 1, ndf
              eqi = id(k,lnodei)
c             write(*,*) 'eqi',eqi
c             write(iow,*) 'eqi',eqi
              if (eqi .gt. 0 ) then
                rbuf     = rbuf + 1
                b( eqi ) = rdatabuf( rbuf )
c               write(iow,*) 'I[',i,'->',rank+1,']: ',
c    &                        rdatabuf(rbuf)
c               write(  *,*) 'I[',i,'->',rank+1,']: ',
c    &                        rdatabuf(rbuf)
              endif
            end do ! k
          end do ! j
          idesp = idesp + getp(i)
        endif
      end do ! i

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      end
