c$Id:$
      subroutine presol(cfr, error)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Initialize data for different solvers

c      Inputs:
c         cfr    - Flag for unsymmetric solve

c      Outputs:
c         error  - Flag for errors detected
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'allotd.h'
      include   'cdata.h'
      include   'compas.h'
      include   'complx.h'
      include   'eqsym.h'
      include   'errchk.h'
      include   'iofile.h'
      include   'iodata.h'
      include   'ldata.h'
      include   'ndata.h'
      include   'part0.h'
      include   'part1.h'
      include   'pathn.h'
      include   'prflag.h'
      include   'setups.h'
      include   'pointer.h'
      include   'comblk.h'

      include   'p_point.h'

      character  tname*5
      logical    cfr, error, palloc, setvar, flags(5)
      integer    ip,kp,len

      save

c     Program initialization

      if(solver) then
        error    = .false.

c       Symmetric part allocation

        if(.not.compfl) then
          write(tname,'(4hTANG,i1)') npart
          if(np(npart).ne.0) then
            call pgetd(tname,point,kp,ip,setvar)
            len = kp*ipr
          else
            len = 0
          endif
          kp = (mr(np(20+npart)+neq-1)+neq)*ipc
          eralloc = .true.
          setvar = palloc(npart,tname,kp,2)
          eralloc = .false.
          if(.not.setvar) then
            write(iow,3007) 'Symmetric',kp*ipr
            if(ior.lt.0) then
              write(*,3007) 'Symmetric',kp*ipr
            endif
            if(ior.lt.0 .and. l.eq.2) then
              error = .true.
              return
            else
              call plstop()
            endif
          endif
          call pzero(hr(np(npart)),kp)

c         Unsymmetric part allocation

          if(cfr) then
            write(tname,'(4hUTAN,i1)') npart
            if(np(npart+4).ne.0) then
              call pgetd(tname,point,kp,ip,setvar)
              len = kp*ipr
            else
              len = 0
            endif
            kp = max(1,(mr(np(20+npart)+neq-1))*ipc)
            eralloc = .true.
            setvar = palloc(npart+4,tname,kp,2)
            eralloc = .false.
            if(.not.setvar) then
              write(iow,3007) 'Unsymmetric',kp*ipr
              if(ior.lt.0) then
                write(*,3007) 'Unsymmetric',kp*ipr
              endif
              if(ior.lt.0 .and. l.eq.2) then
                error = .true.
                return
              else
                call plstop()
              endif
            endif
            call pzero(hr(np(npart+4)),kp)
            na   = np(npart)
            nau  = na + neq*ipc
            nal  = np(npart+4)
            neqs = 1
          else
            na   = np(npart)
            nau  = na + neq*ipc
            nal  = nau
            neqs = neq
          endif
        else

c         Direct out-of-core case

          if(ittyp.eq.-1) then
            iual = ios
            iuau = ios
            fau  = 'Aupper'

            if(cfr) then
              fal  = 'Alower'
              neqs = 1
            else
              fal  = fau
              neqs = neq
            endif

            kp    = nnr - neq
            na    = np(npart)
            nau   = na + neq*ipc
            nal   = nau + kp

            call pzero(hr(nal),kp )

c         Sparse or Preconditioned CG case

          else
            neqs = neq
          endif

          call pzero(hr(np(npart)),nnr)

        endif

c       Save number symmetric equations for partition use

        nqs(npart) = neqs

c     User solver interface

      else
        flags(1) = .true.
        flags(2) = .false.
        flags(3) =  cfr
        flags(4) = .false.
        flags(5) = .false.
        call usolve(flags,hr(1)) ! N.B. hr(*) should not be modified
        error = flags(5)
      endif

c     Format

3007  format(' *WARNING* Insufficient Storage for ',a,' Profile:'/
     &       '           Array size =',i10,'.  Try DIREct,SPARse,',
     &       ' DIREct BLOCk or ITERation')

      end
