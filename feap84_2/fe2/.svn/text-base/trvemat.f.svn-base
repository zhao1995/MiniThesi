c$Id:$
      subroutine trvemat(d, ta,thgrad, hn,hn1,nh, tflux,dd,rhoc, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    19/11/2009
c       1. Set size of send/receive arrays                  07/02/2010
c       2. Set dsend to 6 and drecv to 16                   08/05/2012
c       3. Replace 'omacr1.h' by 'elpers.h'                 21/05/2013
c       4. Use dum(1) to pass rhoc                          15/08/2013
c       5. Set 'drecv' to 17 (for error message)            09/01/2014
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: User Constitutive Model for Thermal Exchanges

c     Input:
c          ta       -  Temperature
c          thgrad(3)-  Thermal gradient
c          hn(nh)   -  History terms at point: t_n
c          nh       -  Number of history terms
c          isw      -  Solution option from element

c     Output:
c          d(*)     -  Material averaged parameters
c          hn1(nh)  -  History terms at point: t_n+1
c          tflux(3) -  Flux at point.
c          dd(3,3)  -  Current material conductivity tangent moduli
c          rhoc     -  Density times specific heat
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'counts.h'
      include   'debugs.h'
      include   'eldata.h'
      include   'elpers.h'
      include   'hdatam.h'
      include   'iofile.h'
      include   'oelmt.h'
      include   'sdata.h'
      include   'setups.h'
      include   'tdata.h'

      include   'mpif.h'

      include   'pointer.h'
      include   'comblk.h'

      logical    setval,palloc
      integer    nh,isw
      integer    a,b
      real*8     ta, rhoc, dum(1)
      real*8     d(*), thgrad(3),hn(nh), hn1(nh), tflux(3),dd(3,3)

c     Compute and output stress (sig) and (moduli)

      save

      if(debug) then
        call udebug('   trvemat',nrecv+1)
        call udebug('   trvemat',isw)
      endif

c     Set values of isw to send information to micro-scale problem

      if(isw.eq.14) then

c       Store material number for each send

        setval = palloc(269,'RVEMA',nsend+1, 1)
        mr(np(269)+nsend) = ma

c       Count number of sends

        sendfl = .true.
        nsend  = nsend + 1

c       Set send receive sizes

        dsend = max(dsend,6)
        drecv = max(drecv,17)

      elseif(isw.eq.4 .or. isw.eq.8) then

        do a = 1,3
          tflux(a) = hn1(a)
        end do ! a

      elseif(isw.eq.3 .or. isw.eq.6 .or. isw.eq.12) then

c       Use current value of stress from array

        if(pltmfl) then

          do a = 1,3
            tflux(a) = hn1(a)
          end do ! a

c       Move deformation gradient to send buffer

        elseif(sendfl) then

          call utstore(hr(np(260)),hn, n, dsend,nrecv, thgrad, ta)

          if(debug) then
            call mprint(thgrad,1,3,1,'TGRAD_send_2')
          endif

c         This is a set to prevent adding non-zeros to tangent/residual

          do a = 1,3
            thgrad(a) = 0.0d0
            do b = 1,3
              dd(b,a) = 0.0d0
            end do ! b
          end do ! a

          rhoc = 0.0d0

c       Put thermal flux and moduli in arrays

        elseif(recvfl) then

          call utrecv(d, hr(np(261)), drecv,nrecv, tflux,dd, rhoc)

c         Save flux for outputs

          if(hflgu) then
            do a = 1,3
              hn1(a) = tflux(a)
            end do ! a
          endif

          if(debug) then
            call mprint(tflux,1,3,1,'FLUX_recv_2')
            call mprint(dd,3,3,3,'D_recv_2')
            dum(1) = rhoc
            call mprint(dum(1) ,1,1,1,'RHOC_recv_2')
          endif

        endif ! Receive

      endif

      end

      subroutine utstore(frvel, hn,n, dsend,nsend, thgrad, ta)

      implicit   none

      include   'iofile.h'

      integer    i, n,dsend,nsend
      real*8     ta, frvel(dsend,*), hn(*), thgrad(3)

      save

      nsend          = nsend + 1
      frvel(1,nsend) = n
      frvel(2,nsend) = hn(1)
      do i = 1,3
        frvel(i+2,nsend) = thgrad(i)
      end do ! i
      frvel(6,nsend) = ta

      end

      subroutine utrecv(d, srvel, drecv,nsend, tflux,dd, rhoc)

      implicit   none

      include   'eldata.h'

      integer    drecv,nsend,a,b,k
      real*8     d(*), srvel(drecv,*), tflux(*),dd(3,3), rhoc

      save

c     Set flux and tangent values from RVE

      nsend = nsend + 1
      do a = 1,3
        tflux(a) = -srvel(a+1,nsend)
      end do ! i
      k = 7
      do a = 1,3
        do b = 1,3
          k = k + 1
          dd(b,a) = srvel(k,nsend)
        end do ! b
      end do ! a

c     Set density and specific heat values from RVE

      d( 4) = srvel(5,nsend)
      d(64) = srvel(6,nsend)
      rhoc  = d(4)*d(64)

      end
