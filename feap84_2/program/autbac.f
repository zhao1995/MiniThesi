c$Id:$
      subroutine autbac ( dtnew )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Define npl(2) for safety                         06/09/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:   Back up solution and adjust to new time step size.

c      Inputs:
c         dtnew - New time step

c      Outputs:
c               - Loads and solution at new time step (dtnew)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'augdat.h'
      include  'cdata.h'
      include  'fdata.h'
      include  'idptr.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'print.h'
      include  'ndata.h'
      include  'part0.h'
      include  'part1.h'
      include  'prld1.h'
      include  'prlod.h'
      include  'rdata.h'
      include  'sdata.h'
      include  'tdata.h'
      include  'comblk.h'

      logical   err
      integer   i,k1,nparta, npl(2)
      real*8    dtnew,propld

      save

c     1.  Backup time step

      ttim   = ttim - dt
      npl(1) = 0
      npl(2) = 0
      if(npld.gt.0) prop = propld(ttim,npl)
      if(prnt) then
        if(npld.gt.1) then
          k1 = npld
        else
          k1 = 0
        endif
        write(iow,2002) ttim,prop,(i,prldv(i),i=1,k1)
        if(ior.lt.0) then
          write(*,2002) ttim,prop,(i,prldv(i),i=1,k1)
        endif
      endif
      rnmax = 0.0d0

c     Reinitialize solution and dynamic vectors for step

      nparta = npart
      do  i = 1,4
        if(nqp(i).gt.0) then
          npart = i
          err   = .false.
          call partpt(npart,err,.true.)

c         Back up dynamic vectors

          call update(mr(id31),hr(np(30)),hr(np(40)),hr(np(42)),
     &                hr(np(26)),fl(9),3)

          call contact (302)

        end if
      end do ! i

c     Restore active partition

      npart = nparta
      call partpt(npart,err,.true.)

c     Reinitialize history vectors for step

      call reshis(mr(np(33)+nen),nen1,numel,1, 2)

c     2.  Advance time to new start point

      if(dtnew .gt. 0.0d0) then

c       Increment time

        dt   = dtnew
        ttim = ttim + dt

        if(npld.gt.0) prop = propld(ttim,npl)
        if(prnt) then
          write(iow,2002) ttim,prop,(i,prldv(i),i=1,k1)
          if(ior.lt.0) then
            write(*,2002) ttim,prop,(i,prldv(i),i=1,k1)
          endif
        endif

c       Move interpolated force vector

        call pmove (hr(np(30)       ),hr(np(30)+  nneq), nneq)
        call pmove (hr(np(30)+2*nneq),hr(np(30)+3*nneq), nneq)

c       Scan partitions for dynamics and contact updates

        nparta = npart
        do  i = 1,4
          if(nqp(i).gt.0) then
            npart = i
            err   = .false.
            call partpt(npart,err,.true.)

c           Update dynamic vectors for transient step

            if(fl(9)) then
              call dsetci(.true.)
              call update(mr(id31),hr(np(30)),hr(np(40)),hr(np(42)),
     &                    hr(np(26)),fl(9),1)
            end if

c           Update slideline arrays

            call contact (301)

          end if
        end do ! i

c       Restore active partition

        npart = nparta
        call partpt(npart,err,.true.)

c       Zero displacement increments for time step

        call pzero(hr(np(40)+nneq),nneq+nneq)

c       Reset history variables

        call reshis(mr(np(33)+nen),nen1,numel,2, 1)

c       Reset parameters for solution step

        augf  = 1.0d0
        rnmax = 0.0d0
        fl( 8) = .false.
        fl(10) = .true.

      endif

c     Format

2002  format(/,'   Computing solution at time ',1p,1e11.4,
     &         ': Total proportional load ',1p,1e11.4:/
     &         '   Individual factors: '/(3x,4(i4,' =',1p,1e12.4)))

      end
