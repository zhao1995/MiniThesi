c$Id:$
      subroutine rvesr(frvel,srvel, ums, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    10/12/2007
c       1. Move 'ik' inside 'i' loop                        17/12/2008
c       2. Add sbuf(7) as flat for first and last send to   04/03/2009
c          each processor. Increase size of sbuf to 17,
c          add variables nsbuf, nrbuf
c       3. Revise to send data to specified RVE types       13/04/2009
c       4. Add temperature to sbuf send from frvel.         13/06/2009
c          Increase dimension of frvel(13,*) & sbuf(19)
c       5. Increase dimension of sbuf(22) and rbuf(72)      10/05/2010
c       6. Increase dimension of sbuf(23); add sbuf(8) for
c          last send (case where only 1-pt/processor        13/12/2010
c       7. Limit lod of sbuf to 'dsend'                     08/05/2012
c       8. Remove ismap                                     18/05/2012
c       9. Reset size of sbuf to 24                         20/07/2012
c      10. Replace 'omacr1.h' by 'elpers.h'                 21/05/2013
c      11. Uset idum(2) to pass integer values              15/08/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Send deformation gradient (F) or strain (eps)  to
c              micro-scale problem.
c              Receive Cauchy stress and tangent moduli.

c     Input:
c        frvel(16,*)  - Array of deformation gradients:
c                     - 1    : Element number
c                     - 2    : Parameter to control micro-scale mesh
c                     - 3-11 : F(1:9)
c                     - 12   : detF
c                     - 13   : ta - temperature change
c                     - 14-16: Thermal gradient

c        ums(ncol,*)  - Send order: ncol = ntasks - 1 (No. RVE's)
c        isw          - Element switch parameter (3, 6, or 12)

c     Output:
c        srvel(45,*)  - Array to store Cauchy stress and material moduli
c                     -  1   : Element number
c                        2-7 : Cauchy stress - sigma(1:6)
c                        8-43: Tangent moduli ctan(1:6,1:6)
c                        44  : Scalar average (rho)
c                        45  : Convergence indicator
c        srvel(72,*)  -  1   : Element number
c                        2-7 : Cauchy stress - sigma(1:6)
c                        8-11: Thermal flux + - flux(1:4)
c                       12-60: Tangent moduli ctan(1:7,1:7)
c                       61-69: Thermal conductivity(1:3,1:3)
c                        70  : Density       - (rho)
c                        71  : Specific heat - (c)
c                        72  : Dissipation   - (r_d)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'counts.h'
      include   'debugs.h'
      include   'elpers.h'
      include   'hdatam.h'
      include   'iofile.h'
      include   'oumatl.h'
      include   'rdata.h'
      include   'sdata.h'
      include   'setups.h'
      include   'tdata.h'

      include   'mpif.h'

      logical    firstfl
      integer    isw, n,nn,mm,nsbuf,nrbuf
      integer    usr_msg, msg_stat(MPI_STATUS_SIZE), ierr
      integer    a, nproce
      integer    ums(ncol,nrow), idum(2)
      real*8     frvel(dsend,*), srvel(drecv,*), sbuf(24), rbuf(72)

      integer    n_err

      save

      data       n_err  / 0 /

c     Set values

      if(debug) then
        call udebug('   rvesr',isw)
      endif

c     Set size of buffers

      nsbuf = dsend + 7
      nrbuf = drecv

c     Set fixed buffer values

      sbuf(2) = nstep
      sbuf(3) = niter
      sbuf(4) = dt
      sbuf(5) = isw
      if(hflgu) then
        sbuf(6) = 1.0d0  ! True
      else
        sbuf(6) = 0.0d0  ! False
      endif

      firstfl = .true.
      do nn = 1,nrow

c       Send deformation gradient to each processor

        usr_msg = 12
        do mm = 1,ncol
          n = abs(ums(mm,nn))
          if(n.gt.0) then
            sbuf(1) = nint(frvel(1,n))
            sbuf(1) = abs(n)
            if(firstfl) then
              sbuf(7) = -1.0d0      ! First send to each processor
            else
              sbuf(7) =  0.0d0      ! Other sends to each processor
            endif
            if(ums(mm,nn).lt.0) then
              sbuf(8) = 1.0d0       ! Last send to a processor
            else
              sbuf(8) = 0.0d0       ! Other sends to a processor
            endif
            do a = 2,dsend
              sbuf(a+7) = frvel(a,n)
            end do ! a

c           Assign processor number

            nproce  = mm

            if(debug) then
              idum(1) = nproce
              idum(2) = nsbuf
              call iprint(idum(1),1,1,1,'PROCESSOR')
              call iprint(idum(2),1,1,1,'BUFFER SIZE')
              call mprint(sbuf   ,1,nsbuf,1,'SBUF')
            endif

c           Send message

            if(debug) then
              call udebug('     MPI_Send:Tgrad',usr_msg)
            endif
            call MPI_SSend( sbuf, nsbuf, MPI_DOUBLE_PRECISION, nproce,
     &                     usr_msg,  MPI_COMM_WORLD, ierr)
            if(ierr.ne.0) then
              write(*,*) ' IERR_send =',ierr
            endif
          endif
        end do ! mm
        firstfl = .false.

c       Receive stress and tangents from processors

        usr_msg = 13
        do mm = 1,ncol

          n = abs(ums(mm,nn))
          if(n.gt.0) then

c           Assign processor number

            nproce  = mm

c           Receive Kirchhoff stress and material moduli

            if(debug) then
              call udebug('     MPI_Recv:Flux',usr_msg)
            endif
            call MPI_Recv( rbuf, nrbuf, MPI_DOUBLE_PRECISION, nproce,
     &                     usr_msg,  MPI_COMM_WORLD, msg_stat, ierr)

            if(debug) then
              idum(1) = nproce
              idum(2) = nrbuf
              call iprint(idum(1),1,1,1,'PROCESSOR')
              call iprint(idum(2),1,1,1,'BUFFER SIZE')
              call mprint(rbuf(1),1,nrbuf,1,'RBUF')
              call mprint(rbuf(8),6,6,6,'CTAU_m')
            endif

c           Store stress and tangent moduli in srvel array

            do a = 1,nrbuf
              srvel(a,n) = rbuf(a)
            end do ! a

c           Check for error: 45 = stress; 17 = thermal problem

            if(nrbuf.eq.45 .or. nrbuf.eq.17) then
              if(rbuf(nrbuf).ne.0.0d0) then
                write(iow,4001) nn,mm,n,rbuf(nrbuf)
                write(  *,4001) nn,mm,n,rbuf(nrbuf)
                n_err = n_err + 1
                if(n_err.gt.100) then
                  write(*,*) ' --> ERROR: More than 100 no convergence'
                  call plstop()
                endif
              endif
            endif

          endif

        end do ! mm
      end do ! nn

c     Output result of point to files for external processing

      if(isw.eq.3 .and. n_pnt.gt.0) then
         write(90,5000) ttim,(frvel(a,n_pnt),a=3,11)
         write(91,5000) ttim,(srvel(a,n_pnt),a=2,7)
         write(92,5001) ttim,(srvel(a,n_pnt),a=8,43)
      endif

c     Format

4001  format(5x,'ERROR: RVE ='i5,' PROCESSOR =',i5,' POINT',i8/
     &       5x,'       NUMBER =',1p,1e12.4)

5000  format(1p,10e16.8)

5001  format(1p,7e16.8/(16x,1p,6e16.8))

      end
