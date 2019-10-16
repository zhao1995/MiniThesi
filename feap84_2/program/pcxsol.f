c$Id:$
      subroutine pcxsol( ctl, prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Separate 'id' & 'eq' on call to pload            27/04/2009
c       2. Change 'nal' to 'point' (in p_point.h'           11/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Perform a forced periodic solutions for general
c              second order systems.  Solves:
c                M a + C v + K d = F exp ( i w t )
c              which is linear.

c              N.B. Problem must start with '*COMPLEX' command to set
c                   'ipc' to 2.

c     Inputs:
c              Command:  CXSOlve,,w,eta
c              ctl(1)  - Frequency: w (omega)
c              ctl(2)  - Damping factor (eta)

c     Outputs:
c              Solution vectors d, v, and a

c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'comblk.h'
      include   'compas.h'
      include   'complx.h'
      include   'eqsym.h'
      include   'evdata.h'
      include   'hdatam.h'
      include   'idptr.h'
      include   'iofile.h'
      include   'part0.h'
      include   'pointer.h'
      include   'p_point.h'
      include   'sdata.h'

      logical    setvar,palloc, prt, scplxfl
      character  tname*5
      integer    i, neqo,njp
      real*8     ctl(3), omega,eta, propom,engy

      save

      data       propom /1.d0/

      if(ipc.eq.2) then

c       Form stiffness, damping, and mass matrices

        if(np(npart+16).eq.0) then
          hflgu  = .false.
          h3flgu = .false.
          ittyp  = -3
          point  = 1
          neqo   = neq
          neq    = nneq
          neqs   = neq

c         Set equation pointers to form entire K, C, and M arrays.

          setvar = palloc(119,'TEMP9',ndf*numnp,1)
          do i = 0,ndf*numnp-1
            mr(np(119)+i) = mr(id31+i)
            mr(id31+i)  = i+1
          end do ! i

c         Compute sparse storage for arrays, allocate stiffness.

          call iters(0,2)

c         Allocate mass and damping storage

          write(tname,'(4hCMAS,i1)') npart
          setvar = palloc(npart+8, tname, nnm, 2)
          write(tname,'(4hDAMP,i1)') npart
          setvar = palloc(npart+16,tname, nnm, 2)

c         Zero arrays

          call pzero(hr(np(npart   )),nnm)
          call pzero(hr(np(npart+ 8)),nnm)
          call pzero(hr(np(npart+16)),nnm)

c         Save complex flag then set for assembly of real parts only

          scplxfl =  cplxfl
          cplxfl  = .false.
          imtyp   =  1

c         Form: Stiffness (isw=3); Damping (isw=9); Mass (isw=5)

          call formfe(np(40),np(26),np(npart   ),point,
     &               .true.,.false.,.false.,.false., 3,1,numel,1)
          call formfe(np(40),np(26),np(npart+ 8),point,
     &               .true.,.false.,.false.,.false., 5,1,numel,1)
          call formfe(np(40),np(26),np(npart+16),point,
     &               .true.,.false.,.false.,.false., 9,1,numel,1)

c         Restore active equation numbers to ID array

          neq  = neqo
          neqs = neq
          do i = 0,ndf*numnp-1
            mr(id31+i) = mr(np(119)+i)
          end do ! i
          setvar = palloc(119,'TEMP9',0,1)

c         Allocate storage for the profile complex tangent matrix

          njp = mr(np(npart+20) + neq - 1) + neq
          setvar = palloc(119,'TEMP9',njp*ipc,2)

c         Restore complex flag to do solutions

          cplxfl = scplxfl
        endif

c       Extract frequency

        omega  = ctl(1)
        eta    = ctl(2)

        if(prt) then
          write(iow,2000) omega,eta
          if(ior.lt.0) then
            write(*,2000) omega,eta
          endif
        endif

c       Assemble REAL force vector

        call pzero(hr(np(26)),ndf*numnp*ipc)
        call pload(mr(np(31)+nneq),mr(id31),hr(np(40)),hr(np(30)),
     &             hr(np(26)),propom,.true.,.false.)

c       Assemble complex tangent and modify rhs for current frequency

        call pzero(hr(np(119)),njp*ipc)
        call pcxasbl(omega,eta,mr(id31),hr(np(26)),hr(np(27)+nneq),
     &               hr(np(119)),hr(np(119)+neq),hr(np(119)+neq*2),
     &               hr(np(119)+neq+njp),mr(np(npart+20)),
     &               hr(np(npart)),hr(np(npart+8)),hr(np(npart+16)),
     &               mr(np(90)),mr(np(91)),nneq)

c       Factor and solve equations

        call datri(hr(np(119)+neq*ipc),hr(np(119)+neq*ipc),
     &             hr(np(119)),mr(np(npart+20)),neq,neq)
        call dasol(hr(np(119)+neq*ipc),hr(np(119)+neq*ipc),
     &             hr(np(119)),hr(np(26)),mr(np(npart+20)),
     &             neq,neq,engy,.false.)

c       Update solution

        call pzero(hr(np(40)),3*ndf*numnp*ipc)
        call update(mr(id31),hr(np(30)),hr(np(40)),hr(np(42)),
     &              hr(np(26)),.false.,2)

c     Error message

      else
        if(ipc.ne.2) then
          write(iow,3000)
          write(ilg,3000)
          if(ior.lt.0) write(*,3000)
          call plstop()
        endif
      endif

c     Formats

2000  format(/'     Frequency   (omega) =',1p,1e12.5/,
     &        '     Damping ratio (eta) =',1p,1e12.5/1x)

3000  format(' *ERROR* PCXSOL: Must declare *COMplex to use command')

      end
