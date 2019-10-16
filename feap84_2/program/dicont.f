c$Id:$
      subroutine dicont(id,numnp,ndf,lflag)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add error output and format 3000,                23/11/2012
c          renumber formats.
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Provides information for displacement control
c               in arc length.

c      Inputs:
c         id(ndf,*) - Equation numbers for each dof
c         numnp     - number of nodal points in mesh
c         ndf       - Number dof/node
c         lflag     - If true, changes arc length

c      Outputs:
c         Equation number of assigned displacement
c         Factors to scale arc length control
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'arclel.h'
      include  'arclei.h'
      include  'arcler.h'
      include  'comfil.h'
      include  'ioincl.h'
      include  'iofile.h'

      character ch*1
      logical   errck, pinput
      integer   numnp,ndf,lflag
      integer   id(ndf,*)
      real*8    td(3)

      save

      if (lflag .ne. 0) go to 100

c     Read in if numerical damping desired or not

      if(ior.lt.0) then
        call pprint('Input: numerical damping (1 = no damping)->')
      endif
      errck = pinput(td,1)
      ndamp = nint(td(1))
      write (iow,2001) ndamp

c     Restart flag

      if (refl) go to 100

50    if (kflag.eq.4.or.kflag.eq.5) then
        if(ior.lt.0) then
          call pprint('Input: Node, DOF, Prescribed Displ.:')
        endif
        errck = pinput(td,3)
        nodis = nint(td(1))
        nddis = nint(td(2))
        alfa0 = td(3)
        if(ior.lt.0) then
          if(nodis.le.0 .or. nodis.gt.numnp .or.
     &       nddis.le.0 .or. nddis.gt.ndf  ) go to 50
        else
          if(nodis.le.0 .or. nodis.gt.numnp .or.
     &       nddis.le.0 .or. nddis.gt.ndf ) then
            write(iow,3000) nodis,nddis
            write(  *,3000) nodis,nddis
            call plstop()
          endif
        endif

c       Set equation number for displacement control

        ndis = id(nddis,nodis)
        write (iow,2000) nodis,nddis,alfa0

c       Check for error

        if(ndis.le.0) then
          if(ior.lt.0) then
            write(*,3001)
            go to 50
          else
            write(iow,3001)
            call plstop()
          endif
        endif
      endif

      return

c     For restart only
c     Any method (displacement control stiff.param. just for chance)

 100  write(iow,2002) rlnew,c0,cs01,cs02

c     Arc length method (any)

      if (kflag.lt.4.or.kflag.eq.6) then
        if(ior.lt.0) then
          write(*,2003) ds0,r
          call pprint('Keep arc-length and load-direction (y or n):')
          read (*,1000) ch
        else
          read (ior,1000,end=900) ch
          irecrd(isf) = irecrd(isf) + 1
          ch          = record(1:1)
        endif
        if(ch.eq.'n' .or. ch.eq.'N') then
          if(ior.lt.0) then
            call pprint('Input: new arc-length, new load direction->')
          endif
          errck = pinput(td,2)
          ds0 = td(1)
          r   = td(2)
          write(iow,2003) ds0,r
        endif
      endif

c     Displacement control

      if (kflag.eq.4.or.kflag.eq.5) then
        if(ior.lt.0) then
          write(*,2004) nodis,nddis,alfa0
          call pprint('Keep displacement control parameters (y or n): ')
          read (*,1000) ch
        else
          read (ior,1000,end=900) record
          irecrd(isf) = irecrd(isf) + 1
          ch          = record(1:1)
        endif
        if(ch.eq.'n' .or. ch.eq.'N') then
          go to 50
        endif
      endif
      return

c     EOF encountered

900   call  endclr ('DICONT',ch)

c     Formats

1000  format(a1)

2000  format(5x,'S i n g l e   D i s p l a c e m e n t   C o n t r o l'/
     &     /10x,'Node Number            =',i8,
     &     /10x,'Degree of Freedom      =',i8,
     &     /10x,'Displacement Increment =',1p,1e13.5)

2001  format(10x,'Numerical damping = ',i3,3x,'(1 = no damping)')

2002  format(/5x,'V a l u e s  for  R e s t a r t:',/,
     &       10x,'Current load level      = ',1p,1e12.4,/,
     &       10x,'S t i f f n e s s  parameter values ',/,
     &       10x,'Stiffness param first step  = ',1p,1e12.4,/,
     &       10x,'Stiffness param 1:prev.step = ',1p,1e12.4,/,
     &       10x,'Stiffness param 2:prev.step = ',1p,1e12.4,/)

2003  format(10x,'Specified arc length        = ',1p,1e12.4/
     *       10x,'Load direction              = ',1p,1e12.4)

2004  format(10x,'Node number                 = ',i3,/,
     &       10x,'Ndof number                 = ',i3,/,
     &       10x,'Prescribed Displacement     = ',1p,1e10.3,/)

3000  format('--> *ERROR* Displacement control for: Node =',i8,' DOF =',
     &       i3,' incorrect')

3001  format('--> *ERROR* Displacement control specified on restrained',
     &       ' node')

      end
