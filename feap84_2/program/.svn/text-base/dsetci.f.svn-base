c$Id:$
      subroutine dsetci(dtfl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add  Euler implicit method                       23/07/2008
c       2. Add pdt(1) for read                              01/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set integration constants for a dynamic analysis

c      Inputs:
c         dtfl    - Time increment flag

c      Outputs:
c         gtan(3) - Integration parameters for tangent
c         c1      - Integration parameters for updates
c         c2      - Integration parameters for updates
c         c3      - Integration parameters for updates
c         c4      - Integration parameters for updates
c         c5      - Integration parameters for updates
c         cc1     - Integration parameters for updates
c         cc2     - Integration parameters for updates
c         cc3     - Integration parameters for updates
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'counts.h'
      include  'ddata.h'
      include  'gltran.h'
      include  'iofile.h'
      include  'tdata.h'
      include  'tdatb.h'
      include  'tdato.h'

      character y
      logical   errck, pinput, dtfl
      integer   i
      real*8    beta,gamm,alph,pdt(1)

      save

c     Set default multiplication factors for residual and tangents

      do i = 1,3
        gtan(i) = 1.d0
      end do ! i
      if(flgexp) then
        expflg = .true.
      else
        expflg = .false.
      endif

c     Check that time increment is non-zero

      if( noi.ne.6 ) then
        if( dtfl. and. dt.eq.0.0d0 .and. ior.lt.0) then
          write(*,2001)
          call pprint('  (y or n) ->')
          read (*,1000) y
          if(y.eq.'y' .or. y.eq.'Y') then
            call pprint('    Specify time increment ->')
            errck = pinput(pdt,1)
            dt    = pdt(1)
            write(iow,2002) dt
            write(  *,2002) dt
          endif
        endif
        if( dt.eq.0.0d0) then
          if(dtfl) then
            write(iow,2003)
            if(ior.lt.0) then
              write(*,2003)
            endif
          endif

c         Set integration parameters to zero

          c1      = 0.0d0
          c2      = 0.0d0
          c3      = 0.0d0
          c4      = 0.0d0
          c5      = 0.0d0
          cc1     = 1.0d0
          cc2     = 1.0d0
          cc3     = 1.0d0

c         Set mass matrix multiplier to zero: permit initial conds.

          gtan(3) = 0.0d0
          return
        endif
      endif

c     Newmark-Beta parameters

      if(noi.eq.1) then

c       Retrieve Newmark parameters

        beta = theta(1)
        gamm = theta(2)

c       Compute integration constants 'c1' to 'c5' for current 'dt'

        c1  = 1.d0/(beta*dt*dt)
        c2  = gamm/(dt*beta)
        c3  = 1.d0 - 1.d0/(beta+beta)
        c4  = 1.d0 - gamm/beta
        c5  = (1.d0 - gamm/(beta+beta))*dt
        cc1 = 1.0d0
        cc2 = 1.0d0
        cc3 = 1.0d0
        if(idyn0.eq.3) then
c         expflg = .true.
        endif

c       Set update parameters for element

        gtan(2) = c2
        gtan(3) = c1

c     Backward Euler parameters

      elseif(noi.eq.2) then
        c1  = 1.d0/dt
        cc1 = 1.0d0
        cc2 = 1.0d0
        cc3 = 1.0d0

c       Set update parameters for element

        gtan(2) = c1
        gtan(3) = c1

c     Conserving HHT parameters

      elseif(noi.eq.3) then

c       Retrieve conserving algorithm parameters

        beta = theta(1)
        gamm = theta(2)
        alph = theta(3)

c       Compute integration constants 'c1' to 'c5' for current 'dt'

        c1  = 1.d0/(beta*dt*dt)
        c2  = gamm/(dt*beta)
        c3  = alph
        c4  = alph*c2
        c5  = 0.5d0/beta
        cc1 = 1.d0
        cc2 = alph
        cc3 = 1.0d0

c       Set update parameters for element

        gtan(1) = c3
        gtan(2) = c4
        gtan(3) = c1

c     Newmark explicit parameters

      elseif(noi.eq.4) then

        gamm = theta(2)

        c1  = 0.5d0*dt*dt
        c2  = gamm*dt
        cc1 = 0.0d0
        cc2 = 0.0d0
        cc3 = 1.0d0

c       Set update parameters for element

        gtan(1) = 0.0d0
        gtan(2) = 0.0d0
        gtan(3) = 1.0d0
        expflg  = .true.

c     Momentum conserving

      elseif(noi.eq.5) then

c       Retrieve conserving algorithm parameters

        beta = theta(1)
        gamm = theta(2)
        alph = theta(3)

c       Compute integration constants 'c1' to 'c5' for current 'dt'

        c1  = 1.d0/(beta*dt*dt)
        c2  = gamm/(dt*beta)
        c3  = alph
        c4  = alph*c2
        c5  = gamm*c1
        cc1 = 1.0d0
        cc2 = alph
        cc3 = 1.0d0

c       Set update parameters for element

        gtan(1) = c3
        gtan(2) = c4
        gtan(3) = c1

c     Static: Generalized Mid-point Configuration

      elseif(noi.eq.6) then

c       Retrieve algorithm parameters

        beta = 0.0d0
        gamm = 0.0d0
        alph = theta(3)

c       Compute integration constants 'c1' to 'c5' for current 'dt'

        c3  = alph
        cc1 = 1.0d0
        cc2 = alph
        cc3 = 1.0d0

c       Set update parameters for element

        gtan(1) = alph
        gtan(2) = 0.0d0
        gtan(3) = 0.0d0

c     First Order: Generalized Mid-point Configuration

      elseif(noi.eq.7) then

c       Retrieve algorithm parameters

        beta = 0.0d0
        gamm = 0.0d0
        alph = theta(3)

c       Compute integration constants 'c1' to 'c5' for current 'dt'

        c1  = 1.d0/dt
        c3  = alph
        cc1 = 1.0d0
        cc2 = alph
        cc3 = 1.0d0

c       Set update parameters for element

        gtan(1) = alph
        gtan(2) = c1
        gtan(3) = c1


c     Central Difference: Explicit

      elseif(noi.eq.8) then

c       Compute integration constants 'c1' to 'c5' for current 'dt'

        c1  = 0.5d0*(dt + dtold)
        cc1 = 0.0d0
        cc2 = 0.0d0
        cc3 = 1.0d0

c       Set update parameters for element

        gtan(1) = 0.0d0
        gtan(2) = 0.0d0
        gtan(3) = 0.0d0
        expflg  = .true.

c     Backward Difference Formula (BDF2) parameters

      elseif(noi.eq.9) then
        if(nstep-nstepa.gt.1) then
          c4  = (2.d0*dt + dtold)/(dt + dtold) ! alpha
        else
          c4  = 1.d0
        endif
        c1  = c4 / dt                          ! alpha / dt
        c2  = 1.d0/dt                          ! 1 / dt
        c3  = 1.d0 - c4                        ! 1 - alpha

        cc1 = 1.0d0
        cc2 = 1.0d0
        cc3 = 1.0d0

c       Set update parameters for element

        gtan(2) = c1
        gtan(3) = c1

c     Backward Euler for 2nd order problems

      elseif(noi.eq.10) then

c       Compute integration constants

        cc1 = 1.d0
        cc2 = 1.d0
        cc3 = 1.0d0

c       Set update parameters for element

        gtan(2) = 1.d0/dt
        gtan(3) = 1.d0/(dt*dt)

c     User integrator

      elseif(noi.eq.-1) then

        call usetci()

      endif

1000  format(a)
2001  format(' ** WARNING ** Current time increment is zero, change?')
2002  format('    Time increment reset to:',e15.5)
2003  format(' ** WARNING ** Current time increment and c(i) are zero')

      end
