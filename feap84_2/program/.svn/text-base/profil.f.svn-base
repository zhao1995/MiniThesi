c$Id:$
      subroutine profil (jp,idl,id,ix,iop,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'nummat' to call of 'rstprf'                 20/07/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute profile of global arrays

c      Inputs:
c        ix(*)  - Element nodal connection list
c        iop    - Switch to control operation
c                  = 1 to set up equation numbers of dof's
c                  = 2 to compute the column/row lengths and profile.
c        prt    - Flag, print solution properties if true

c      Scratch:
c        idl(*) - Array to hold temporary information

c      Outputs:
c        id(*)  - Equation numbers for degree of freedoms     (iop = 1)
c        jp(*)  - Pointer array to row/column ends of profile (iop = 2)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'complx.h'
      include  'ddata.h'
      include  'iofile.h'
      include  'part0.h'
      include  'pointer.h'
      include  'rigid1.h'
      include  'rjoint.h'
      include  'sdata.h'
      include  'comblk.h'

      include  'p_point.h'

      logical   prt, setvar, palloc, getfl
      character dnam*5, itype(10)*10
      integer   i, iop, mi,mm,nad, ndft, jp(*),idl(*),id(*),ix(*)

      save

      data      itype /'Static   ','Newmark  ','BackODE_1',
     &                 'Cons. HHT','Explicit ','Cons Newm',
     &                 'Mid ODE_1','Mid ODE_1','Cent Diff',
     &                 'User     '/

c     Set up equation numbers

      if(iop.eq.1) then

        call seteq(id,mr(np(99)),mr(np(182)),mr(np(100)),mr(np(101)),
     &             prt)

c       User modifications

        call uprofil(jp,idl,id,ix,iop,prt)

c       Reset size of solution array if necessary

        setvar = palloc(26,'DR   ',max(neq,numnp*max(ndf*ipc,ndm)),2)

        write(dnam,'(a2,i1)') 'JP',npart
        call pgetd(dnam, point,nad,i, getfl)
        if(neq.gt.nad) then
          setvar = palloc( 20+npart, dnam, neq, 1)
        endif

c     Compute column heights

      elseif(iop.eq.2) then

        call rstprf(jp,idl,id,ix,mr(np(32)), mr(np(99)),mr(np(100)),
     &              mr(np(101)),ndf,nen1,nen,neq,numnp,numel,nummat)

c       Interface adjustments

        if(np(210).ne.0) then
          call isetprf(jp,idl,id,ix,mr(np(32)),mr(np(210)))
        endif

c       User modifications

        call uprofil(jp,idl,id,ix,iop,prt)

c       Compute diagonal pointers for profile

        call nwprof(jp,neq)

c       Output statistics

        if(neq.gt.0 .and. prt) then
          nad = jp(neq)
          mm = (nad+neq)/neq
          if(noi.ge.0) then
            mi = noi+1
          else
            mi = 10
          endif
          ndft = 0
          do i = 1,ndf
            if( ndfp(i) .eq. npart ) ndft = ndft + 1
          end do ! i
          write(iow,2001) ndm,numnp,ndft,numel,nadd,nummat,nrbody,
     &                    neq,numjts,nad,itype(mi),mm
          if(ior.lt.0)
     &      write(*,2001) ndm,numnp,ndft,numel,nadd,nummat,nrbody,
     &                    neq,numjts,nad,itype(mi),mm
        endif
      endif

c     Format

2001  format(/5x,'E q u a t i o n / P r o b l e m   S u m m a r y:'//
     & 7x,'Mesh dimension  (ndm) =',i10,' :  Number nodes     =',i16/
     & 7x,'Number dof/node (ndf) =',i10,' :  Number elements  =',i16/
     & 7x,'Element eqs.   (nadd) =',i10,' :  Number materials =',i16/
     & 7x,'Number rigid bodies   =',i10,' :  Number equations =',i16/
     & 7x,'Number joints         =',i10,' :  Number tang terms=',i16/
     & 7x,'Transient integrator  = ',a9,' :  Average col. ht. =',i16/)

      end
