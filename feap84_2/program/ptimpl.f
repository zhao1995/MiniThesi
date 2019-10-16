c$Id:$
      subroutine ptimpl

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Flush buffer for each write                      15/11/2008
c       2. Add function pflush to do flushes                25/02/2009
c       3. Set pltmfl to .true. before calls to formfe.     05/03/2009
c          Used to force stress retrieval from history data.
c       4. Separate 'id' & 'eq' on call to pload            27/04/2009
c       5. Add 'hist' option to tplot                       25/06/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Store time history plot information for converged step

c      Inputs:
c        none

c      Outputs:
c        none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'aceang.h'
      include   'arcler.h'
      include   'cdata.h'
      include   'complx.h'
      include   'elcount.h'
      include   'eltran.h'
      include   'endata.h'
      include   'fdata.h'
      include   'gltran.h'
      include   'hdatam.h'
      include   'idptr.h'
      include   'iofile.h'
      include   'part0.h'
      include   'part7.h'
      include   'print.h'
      include   'prlod.h'
      include   'ptdat1.h'
      include   'ptdat2.h'
      include   'ptdat3.h'
      include   'ptdat4.h'
      include   'ptdat5.h'
      include   'ptdat6.h'
      include   'ptdat7.h'
      include   'ptdat8.h'
      include   'ptdat9.h'
      include   'ptdata.h'
      include   'ptdatb.h'
      include   'ptdatc.h'
      include   'sdata.h'
      include   'tdata.h'
      include   'tdato.h'
      include   'pointer.h'
      include   'comblk.h'

      logical    err, errck, fl9sv(4), oflag
      integer    i, ii, iunit
      real*8     dtsav, ctl(25)

      save

      if(max(naplts,ncplts,ndplts,neplts,nlplts,nqplts,nrplts,
     &       nsplts,nhplts,ntplts,nmplts,nuplts,nvplts,nchplts)
     &                                                  .gt.0 ) then

        ntstep = ntstep + 1

c       Check for active output increment

        if(mod(ntstep-1,ntincr).eq.0) then

c         Set history update flag to false (no updates)

          hflgu  = .false.
          h3flgu = .false.

c         Set transient parameters for current

          if(fl(9)) call dsetci(.true.)
          do i = 1,3
            ctan(i) = gtan(i)
          end do ! i

          if(max(nmplts,nsplts,nhplts,nqplts,nrplts,ntplts,nuplts)
     &                                                      .gt.0) then
            if(.not.rfl) then
              if(ntplts.gt.0) then
                do i = 1,ntplts
                  tpl(i) = 0.0d0
                end do ! i
              endif
              dtsav = dt
              dt    = dtold
              do i = 1,4
                fl9sv(i) = flp(9,i)
                flp(9,i) = ofl9(i)
              end do ! i
              fl(9) = ofl9(npart)
              call pzero(hr(np(26)),nneq*ipc)
              pltmfl = .true.
              call formfe(np(40),np(26),np(26),np(26),
     &                   .false.,.true.,.false.,.true.,6,1,numel,1)
              pltmfl = .false.
              if(ntplts.gt.0) then
                call ptsumpl(hr(np(43)),hr(np(26)))
              endif
              rfl    = .true.
              dt     =  dtsav
              do i = 1,4
                flp(9,i) =  fl9sv(i)
              end do ! i
              fl(9) = flp(9,npart)
            end if
          end if

c         Set file name for displacements

          if(ndplts.gt.0) then
            iunit = 31
            oflag = ndplts.gt.100
            call pltmv(dpl,idpl,hr(np(40)),ndplts,1.d0)
            call ptmplt('dis', ttim, dpl,ndplts, ntstep, iunit, oflag)
          end if

c         Set file name for velocities

          if(nvplts.gt.0) then
            iunit = 36
            oflag = nvplts.gt.100
            call pltmv(dpl,ivpl,hr(np(42)),nvplts,1.d0)
            call ptmplt('vel', ttim, dpl,nvplts, ntstep, iunit, oflag)
          end if

c         Set file name for accelerations

          if(naplts.gt.0) then
            iunit = 41
            oflag = naplts.gt.100
            call pltmv(dpl,iapl,hr(np(42)+nneq),naplts,1.d0)
            call ptmplt('acc', ttim, dpl,naplts, ntstep, iunit, oflag)
          end if

c         Set file name for stresses

          if(nsplts.gt.0) then
            iunit = 46
            oflag = nsplts.gt.100
            call ptmplt('str', ttim, spl,nsplts,ntstep, iunit, oflag)
          end if

c         Set file name for history

          if(nhplts.gt.0) then
            iunit = 51
            oflag = nhplts.gt.100
            call ptmplt('his', ttim, hpl,nhplts,ntstep, iunit, oflag)
          end if

c         Set file name for reactions

          if(nrplts.gt.0) then
            iunit = 56
            oflag = nrplts.gt.100
            call pltmv(rpl,irpl,hr(np(26)),nrplts,-1.d0)
            call ptmplt('rea', ttim, rpl,nrplts, ntstep, iunit, oflag)
          end if

c         Set file name for energys/momenta

          if(neplts.gt.0) then

c           Compute energy from elements

c           epl(1),(2),(3): linear momentum components
c           epl(4),(5),(6): angular momentum components
c           epl(7):         kinetic energy
c           epl(8):         potential energy
c           epl(9):         work for external loads
c           epl(10):        total energy
c           epl(11):        angular momentum norm

            do i = 1,9
              epl(i) = 0.0d0
            end do ! i

c           Form nodal load vector

            call ploa1(ttim,dt)
            call pload(mr(np(31)+nneq),mr(id31),hr(np(40)),hr(np(30)),
     &                 hr(np(26)),prop*rlnew,.true.,.false.)
c           call ploade(mr(id31),hr(np(30)),hr(np(40)), epl(9))

c           Compute momentum and energy from elements and rigid bodies

            pltmfl = .true.
            call formfe(np(40),np(26),np(26),np(26),
     &                 .false.,.false.,.false.,.false.,13,1,numel,1)
            pltmfl = .false.

c           Compute total values and output results

c           epl(10) = epl(7) + epl(8) - epl(9)
            epl(10) = epl(7) + epl(8)
            epl(11) = sqrt(epl(4)**2 + epl(5)**2 + epl(6)**2)

            ang(1) = ang(2)

c           Determine if all momentum/energy are to be output to file

            errck = .false.
            err   = .false.
            do i = 1,neplts
              if(iepl(1,i).eq.0) then
                errck = .true.
                if(.not.err .and. iepl(2,i).gt.0
     &                      .and. ior.lt.0 .and. prnt) then
                  err = .true.
                endif
              endif
            end do ! i

c           Output total energy to screen

            if(err) write(*,*) '  TOTAL ENERGY =',epl(10)

c           Load active energy outputs

            if(errck) then
              do i = 1,11
                ctl(i) = epl(i)
              end do ! i
              ii = 11
            else
              do i = 1,neplts
                ctl(i) = epl(iepl(1,i))
              end do ! i
              ii = neplts
            endif
            iunit = 30
            oflag = .true.
            oflag = .false.
            call ptmplt('ene', ttim, ctl, ii   , ntstep, iunit, oflag)
            call pzero(epl,200)
          end if

c         Set file name for arclength

          if(nlplts.gt.0) then
            iunit = (nlplts -1)/20 + 61
            oflag =  nlplts.gt.10
            call pltmv(lpl,ilpl,hr(np(40)),nlplts,1.d0)
            call ptmplt('arc', rlnew*prop, lpl,nlplts, ntstep,
     &                  iunit, oflag)
          end if

c         Set file name for contacts

          if(ncplts.gt.0) then
            iunit = 62
            oflag = ncplts.gt.80
            do i = 1,ncplts
              cpl(i) = 0.0d0
            end do ! i
            call contact (206)
            call ptmplt('con', ttim, cpl,ncplts, ntstep, iunit, oflag)
          end if

c         Set file name for totals

          if(ntplts.gt.0) then
            iunit = 66
            oflag = ntplts.gt.100
            call ptmplt('sum', ttim, tpl,ntplts, ntstep, iunit, oflag)
          end if

c         Set file name for material states

          if(nmplts.gt.0) then
            ii = 0
            do i = 1,nmplts
              if(impl(1,i).ne.0) then
                mpl(ii+1) = nomats(1,impl(1,i))
                mpl(ii+2) = nomats(2,impl(1,i))
                ii        = ii + 2
              endif
            end do ! i
            iunit = 71
            oflag = ii.gt.100
            call ptmplt('mat', ttim, mpl, ii, ntstep, iunit, oflag)
          end if

c         Set file name for user stresses

          if(nuplts.gt.0) then
            iunit = 76
            oflag = nuplts.gt.100
            call ptmplt('use', ttim, upl,nuplts, ntstep, iunit, oflag)
          end if

c         Set file name for reaction sums

          if(nqplts.gt.0) then
            iunit = 81
            oflag = nqplts.gt.40
            do i = 1,nqplts
              call preacsm(hr(np(26)),ndf,numnp,iqpl(1,i),qpl(i) )
            end do ! i
            call ptmplt('rsm', ttim, qpl,nqplts, ntstep, iunit, oflag)
          end if

c         Set file name for contact history variables

          if(nchplts.gt.0) then
            iunit = 86
            oflag = nchplts.gt.100
            do i = 1,nchplts
              chpl(i) = 0.0d0
            end do ! i
            call contact (206)
            call ptmplt('chs', ttim, chpl,nchplts, ntstep, iunit, oflag)
          end if

          call pflush(iunit)

        end if  ! End of increment check

      end if

      end
