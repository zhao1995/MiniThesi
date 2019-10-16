c$Id:$
      subroutine pmodal(dt,d,phi, fs, y, id,fn, u,ud, mf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Modal solution of linear transient problems
c                with possibility of active control

c      Inputs:
c         dt         - Time increment
c         d(*)       - Frequencies squared
c         phi(vneq,*)- Mass orthonormal eigenvectors
c         y(mf,3)    - Eigensolution at time t_n
c         id(*)      - Pointer array to active equations
c         fn(nneq,2) - Force/displacement vector at t_n+1 and t_n
c         mf         - Number of eigenpairs (converged)
c         neq        - Number of active components in eigenvectors

c      Outputs:
c         y(mf,3)    - Eigensolution at time t_n+1
c         u(*)       - Displacement at time t_n+1
c         ud(*,1)    - Velocity at time t_n+1
c         ud(*,2)    - Acceleration at time t_n+1
c         ud(*,nrk)  - Displacement at time t_n+a (some algorithms)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'ddata.h'
      include  'endata.h'
      include  'modcon.h'
      include  'pfeapb.h'
      include  'pointer.h'
      include  'sdata.h'
      include  'comblk.h'

      logical   setvar,palloc
      integer   i,j,jj, mp, mf, id(*)
      real*8    dt,omg,odt,cs,sn, xi,g0,g1, aa,bb, expt, m,c,k
      real*8    d(*),phi(vneq,*), fs(mf,*), y(mf,*), fn(nneq,*)
      real*8    u(*),ud(nneq,*)

      save

c     Set storge for modal quantities

      setvar = palloc(184,'CTEMP',mf,2)
      setvar = palloc(185,'KTEMP',mf,2)
      setvar = palloc(186,'MTEMP',mf,2)

c     Check to see if control values are available

      if(np(169).ne.0) then
        do i = 0,mf-1
          hr(np(186)+i) = hr(np(169)+i) + 1.0d0
          hr(np(184)+i) = hr(np(170)+i) + rayla0 + rayla1*d(i+1)
          hr(np(185)+i) = hr(np(171)+i) + d(i+1)
        end do ! i
      else
        do i = 0,mf-1
          hr(np(186)+i) = 1.0d0
          hr(np(184)+i) = rayla0 + rayla1*d(i+1)
          hr(np(185)+i) = d(i+1)
        end do ! i
      end if

c     Project eigen-forces

      do i = 1,mf
        fs(i,1) = fs(i,2)
        fs(i,2) = 0.0d0
      end do ! i

      do j = 1,nneq
        if(id(j).gt.0 .and. fn(j,1).ne.0.0d0) then
          jj = id(j)
          do i = 1,mf
            fs(i,2) = fs(i,2) + phi(jj,i)*fn(j,1)
          end do ! i
        endif
      end do ! j

c     Parallel version send/receive fs(*,2)

      if(pfeap_on) then

c       Sending fs(*,2) data asynchronously

        call pfeapsr(fs(1,2),fs(1,3),mf,.true.)

      endif ! pfeap_on

c     Set support displacements

      mp = 0
      if(noi.eq.5 .and. np(189).ne.0) then
        do i = 0,nneq-1
          mp = max(mp,mr(np(125)+i))
        end do ! i

c       Parallel communication

        if(pfeap_on) then
          call pfeapmi(mp)
        endif ! pfeap_on
        setvar = palloc(124,'PROBS',mp,1)
        do i = 0,nneq-1
          j = mr(np(125)+i)
          if(j.gt.0) then
            mr(np(124)+j-1) = mr(np(29)+i)
          endif
        end do ! i
        call basdis(mr(np(124)),hr(np(189)),dt,mp)
      endif

c     Advance projected solutions

      do i = 1,mf
        if(dt.gt.0.0d0) then

          k = hr(np(185)-1+i)
          m = hr(np(186)-1+i)
          c = hr(np(184)-1+i)

c         E-M Method: hard wired to beta=0.5, gamma=1.0, alpha=0.5

          if(noi.eq.5) then

            aa  = 2.d0/dt**2*m + c/dt + 0.5d0*k
            bb  = 0.5d0*(fs(i,1) + fs(i,2))
     &          + 2.d0/dt**2*m*(y(i,1) + dt*y(i,2)) + c/dt*y(i,2)
     &          - 0.5d0*k*y(i,1)

c           Check for multiple support excitations

            if(np(126).ne.0) then
              call baslod(hr(np(126)),hr(np(189)),mf,mp,i, bb)
            endif

            g0  = 2.d0*y(i,1)/dt + y(i,2)
            g1  = y(i,2)
            y(i,1) = bb/aa
            y(i,2) = 2.d0*y(i,1)/dt - g0
            y(i,3) = (y(i,2) - g1)/dt

          else

            omg = sqrt(abs(k/m))
            xi  = c/(2.0d0*omg*m)

            g1  = (fs(i,2) - fs(i,1))/(k*dt)
            g0  = (fs(i,1) - c*g1)/k

            aa  = y(i,1) - g0
            bb  = y(i,2) - g1 + xi*omg*aa

            if(xi.gt.1.d0) then

              odt = omg*sqrt(xi*xi - 1.d0)
              cs  = cosh(odt*dt)
              sn  = sinh(odt*dt)/odt

            elseif(xi.lt.1.d0) then

              odt = omg*sqrt(1.0d0 - xi*xi)
              cs  = cos(odt*dt)
              sn  = sin(odt*dt)/odt

            else

              odt = 0.0d0
              cs  = 1.d0
              sn  = dt

            endif

c           Solution for 't_n+1'

            expt   = exp(-xi*omg*dt)
            y(i,1) = g1*dt + g0 + (aa*cs + bb*sn)*expt
            y(i,2) = g1 - xi*omg*expt*(aa*cs + bb*sn)
     &             + expt*(sign(1.d0,xi-1.d0)*aa*odt*odt*sn + bb*cs)
            y(i,3) = (fs(i,2)-k*y(i,1)-c*y(i,2))/m

          endif

        else
          y(i,3) = (fs(i,2)-k*y(i,1)-c*y(i,2))/m
        endif

      end do ! i

c     Compute displacement, velocity, and acceleration

      call pmodis(id,phi,y,mr(np(125)),hr(np(127)),hr(np(189)),
     &            hr(np(27)+ndf*numnp),mp,neq,ndf,numnp,u,ud)

c     Set flag to force reaction/stress computation for tplot

      rfl = .false.

c     Exit

      end
