c$Id:$
      subroutine partpt(np1,tfl,pfl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Controls pointers to arrays and variables active
c               in each solution partition.

c      Inputs:
c         np1   - Partition number to activate
c         tfl   - Initial allocation flag
c         pfl   - Set proportional load table if true

c      Outputs:
c         Arrays through common blocks
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'arclel.h'
      include  'arclei.h'
      include  'arcler.h'
      include  'cdata.h'
      include  'compas.h'
      include  'complx.h'
      include  'ddata.h'
      include  'debugs.h'
      include  'endata.h'
      include  'eqsym.h'
      include  'fdata.h'
      include  'gltran.h'
      include  'idptr.h'
      include  'iofile.h'
      include  'mxsiz.h'
      include  'ndata.h'
      include  'part0.h'
      include  'part1.h'
      include  'part2.h'
      include  'part3.h'
      include  'part5.h'
      include  'part6.h'
      include  'part7.h'
      include  'part8.h'
      include  'part9.h'
      include  'pointer.h'
      include  'prlod.h'
      include  'rdata.h'
      include  'rdat0.h'
      include  'rjoint.h'
      include  'sdata.h'
      include  'tdatb.h'

      logical   pfl,tfl
      integer   i,j,np1,npo

      save

c     Set pointers for allocation of loads/partition data

      if(np1.ge.1 .and. np1.le.5) then
        if(tfl) then
          npo         = 0
          aenp(np1)   = 0.0d0
          aolp(np1)   = 0.0d0
          rnp(np1)    = 0.0d0
          arcp(np1)   = .false.
          flp( 1,np1) = .false.
          flp( 2,np1) = .false.
          flp( 3,np1) = .true.
          flp( 4,np1) = .true.
          flp( 5,np1) = .true.
          flp( 6,np1) = .true.
          flp( 7,np1) = .true.
          flp( 8,np1) = .false.
          flp( 9,np1) = .false.
          flp(10,np1) = .false.
          flp(11,np1) = .false.
          flp(12,np1) = .true.
          npld        = 0
          do i = 1,7
            ap(i,1)     = 0.d+0
            pa(i,1,np1) = 0.d+0
          end do ! i
          ap(2,1)      = 1.d+8
          ap(3,1)      = 1.d+0
          pa(2,1,np1)  = 1.d+8
          pa(3,1,np1)  = 1.d+0
          iexp(1)      = 0
          ik(1)        = 1
          nppld(np1)   = 0
          ipexp(1,np1) = 0
          ipk(1,np1)   = 1
          nrkp(np1)    = 0
          nrcp(np1)    = 0
          nrmp(np1)    = 0
          noip(np1)    = 0
          thetp(1,np1) = 0.0d0
          thetp(2,np1) = 0.0d0
          thetp(3,np1) = 1.0d0
          nmethp(np1)  = 1
          tfl         = .false.
        endif
        if(npo.ne.0) then
          aenp(npo)    = aengy
          aolp(npo)    = aold
          rnp (npo)    = rnmax
          arcp(npo)    = arcf
          do i = 1,12
            flp(i,npo) = fl(i)
          end do ! i
          do i = 1,4
            thetp(i,npo) = theta(i)
          end do ! i
          nrkp(npo)    = nrk
          nrcp(npo)    = nrc
          nrmp(npo)    = nrm
          noip(npo)    = noi
          nlp(npo)     = nl
          nmp(npo)     = nm
          nap(npo)     = na
          naup(npo)    = nau
          nalp(npo)    = nal
          ncp(npo)     = nc
          nxp(npo)     = nx
          nmethp(npo)  = nmeth
          nppld(npo)   = npld
          do i = 1,50
            do j = 1,7
              pa (j,i,npo) = ap(j,i)
            end do ! j
            ipexp(i,npo) = iexp(i)
            ipk(i,npo)   = ik(i)
          end do ! i
        endif
        do i = 1,12
          fl(i) = flp(i,np1)
        end do ! i
        do i = 1,4
          theta(i) = thetp(i,np1)
        end do ! i
        nrk      = nrkp(np1)
        nrc      = nrcp(np1)
        nrm      = nrmp(np1)
        noi      = noip(np1)
        nl       = nlp(np1)
        nm       = nmp(np1)
        na       = nap(np1)
        nau      = naup(np1)
        nal      = nalp(np1)
        nc       = ncp(np1)
        nx       = nxp(np1)
        nmeth    = nmethp(np1)
        if(pfl) then
          npld  = nppld(np1)
          do i = 1,50
            do j = 1,7
              ap(j,i) = pa (j,i,np1)
            end do ! j
            iexp(i) = ipexp(i,np1)
            ik(i)   = ipk(i,np1)
          end do ! i
        endif
        aengy = aenp(np1)
        aold  = aolp(np1)
        if(rnmax.eq.0.0d0) rnp(np1) = 0.0d0
        rnmax =  rnp(np1)
        arcf  =  arcp(np1)
        rfl   = .false.
        mxpro = mxprop(np1)
        mxneq = mxneqp(np1)
        neq   = nqp(np1)
        neqr  = nqr(np1)
        neqs  = nqs(np1)
        ittyp = nittyp(np1)
        itol  = nitolp(np1)
        atol  = natolp(np1)
c       nau   = na + neq*ipc
        if(npart.ne.np1) fl(11) = .false.
        gtan(1) = 1.0d0
        gtan(2) = 0.0d0
        gtan(3) = 0.0d0
        cc1     = 1.0d0
        cc2     = 1.0d0
        cc3     = 1.0d0

c       Set partition pointers

        do i = 1,ndf
          ndfp(i) = ndfst(i,np1)
          ndfg(i) = ndfp(i)
        end do ! i
        id31 = idpt(np1)

        npart   = np1
        npo     = np1
      else
        write(  *,3000) npart
        write(ilg,3000) npart
        call plstop()
      endif

3000  format(' *ERROR* PARTPT: Bad Partition Number',i3)

      end
