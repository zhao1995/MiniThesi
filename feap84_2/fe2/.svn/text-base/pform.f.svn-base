c$Id:$
      subroutine pform(ul,xl,tl,ld,p,s,ie,d,id,x,ix,rben,f,t,jp,
     &                 u,ud,b,a,al,ndd,nie,ndf,ndm,nen1,nst,
     &                 aufl,bfl,alfl,dfl,isw,nn1,nn2,nn3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove interupt tests (intf and lintr)           17/11/2006
c       2. Remove unused format 1000                        13/12/2006
c       3. Add 'jsw.eq.3' to check for numerical tangent    10/02/2007
c       4. Add RVELM, recvfl, and nproc sets.               20/07/2008
c       5. Add check for RVE second pass on 'pltmfl'        06/03/2009
c       6. Increase dimension of FRVEL to 13                13/06/2009
c       7. Add set of nrvn for velocity at t_n              21/05/2011
c       8. Change 'a_avg' to 'v_avg'; zero v_rho and v_c    09/05/2012
c       9. Replace 'omacr1.h' by 'elpers.h'                 21/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute element arrays and assemble global arrays

c      Inputs:
c         ie(nie,*)   - Element information for material set
c         d(ndd,*)    - Material set parameters
c         id(ndf,*)   - Equation numbers for each active dof
c         x(ndm,*)    - Nodal coordinates of mesh
c         ix(nen1,*)  - Element nodal connections of mesh
c         rben(*)     - Rigid/modal element indicator
c         f(ndf,*,2)  - Nodal force and displacement values
c         t(*)        - Nodal temperature values
c         jp(*)       - Pointer array for row/columns of tangent
c         u(*)        - Nodal solution values
c         ud(*)       - Nodal rate values
c         ndd         - Dimension for d array
c         nie         - Dimension for ie array
c         ndf         - Number dof/node
c         ndm         - Spatial dimension of mesh
c         nen1        - Dimension for ix array
c         nst         - Dimension for element array
c         aufl        - Flag, assemble coefficient array if true
c         bfl         - Flag, assemble vector if true
c         alfl        - Flag, coefficient array unsymmetric if true
c         dfl         - Flag, assemble reactions if true
c         isw         - Switch to control quantity computed
c         nn1         - First element number to process
c         nn2         - Last element number to process
c         nn3         - Increment to nn1

c      Local element arrays:
c         ul(*)       - Element solution and rate values
c         xl(*)       - Element nodal coordinates
c         tl(*)       - Element nodal temperatures
c         ld(*)       - Element local/global equation numbers
c         p(nst,*)    - Element vector
c         s(nst,*)    - Element array

c      Outputs:
c         b(*)        - Global vector
c         a(*)        - Global matrix, diagonal and upper part
c         al(*)       - Global matrix, lower part
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'arcler.h'
      include   'cdata.h'
      include   'ddata.h'
      include   'debugs.h'
      include   'elcount.h'
      include   'eldata.h'
      include   'elpers.h'
      include   'iofile.h'
      include   'hdata.h'
      include   'hdatam.h'
      include   'oelmt.h'
      include   'part0.h'
      include   'pointer.h'
      include   'prld1.h'
      include   'prlod.h'
      include   'rigid1.h'
      include   'setups.h'
      include   'comblk.h'

      logical    aufl,bfl,alfl,dfl,efl,gfl,rfl,arotflg
      logical    rvfl, svfl, setval, palloc
      integer    isw, jsw, ksw
      integer    i, nn1, nn2, nn3, nst, nneq
      integer    ndf, ndm, ndd, nie, nen1
      integer    ld(*), ie(nie,*), id(ndf,*), ix(nen1,*), rben(*), jp(*)
      real*8     xl(ndm,*), p(nst,*), s(nst,*), d(ndd,*), ul(ndf,nen,*)
      real*8     x(ndm,*) ,f(ndf,numnp),u(ndf,*),ud(*),t(*),tl(*)
      real*8     prope, b(*), a(*), al(*)

      save

      if(debug) then
        call udebug(' pform',isw)
      endif

c     Initialize data

      v_avg  = 0.0d0
      v_rho  = 0.0d0
      v_c    = 0.0d0
      sig_33 = 0.0d0

      nproc  = 0
      svfl  = ntasks.gt.1 .and. rank.eq.0
      rvfl  = .false.
      if(isw.eq.14) then
        nsend = 0
      endif

c     Set element proportional loading value

      prope = (theta(3)*(prop - propo) + propo)
      if(isw.ne.23) then        ! Element arclength case
        prope = prope*rlnew
      endif

c     Set nh1, nh2, nh3 pointers for local history variables

      nh1 = np(50)
      nh2 = np(51)
      nh3 = np(52)

c     Set program and user material count parameters

      do i = 1,10
        nomats(1,i) = 0
        nomats(2,i) = 0
        unmats(1,i) = 0
        unmats(2,i) = 0
      end do ! i

c     Initialize RVELM array -- stores list of elements to send-receive

      if(isw.eq.14 .and. svfl) then
        if(np(259).eq.0) then
          setval = palloc(259,'RVELM',numel,1)
        endif
        do n = 1,numel
          mr(np(259)+n-1) = 0
        end do ! n
      endif

c     Set flags and parameters for solution

      iel = 0
      if(isw.eq.3 .and. bfl) then
        efl = .true.
      elseif(isw.eq.6 .and. .not.dfl) then
        efl = .true.
      else
        efl = .false.
      endif
      arotflg = aufl .or. bfl

c     Set flags for modal base excitations

      if(isw.eq.19) then
        if(bfl) efl = .true.
        jsw = 5
        ksw = 5
        gfl = .false.

c     Other cases

      else
        jsw = isw
        ksw = 3
        gfl = .true.
      endif

      nneq   = numnp*ndf
      nrkn   = nrk*nneq - nneq
      nrcn   = nrc*nneq - nneq
      nrmn   = nrm*nneq - nneq
      nrvn   = nrt*nneq - nneq - nneq

c     Check for a rigid body

      rfl = rbody .and. (nrbprt.eq.npart)

c     Loop over active elements: Compute any local elements and the
c     deformation gradient for send/receive microscale points.

      call pforma(ul,xl,tl,ld,p,s,ie,d,id,x,ix,rben,f,t,jp,
     &            u,ud,b,a,al,ndd,nie,ndf,ndm,nen1,nst,nneq,
     &            prope,arotflg,aufl,bfl,alfl,dfl,efl,gfl,rfl,
     &            rvfl,svfl,jsw,ksw,nn1,nn2,nn3)

c     Parallel send/receives for RVE

      if(svfl .and. nsend.gt.0 .and. .not.pltmfl) then
        if(isw.eq.14) then
          setval = palloc(260,'FRVEL', dsend*nsend, 2)
          setval = palloc(261,'SRVEL', drecv*nsend, 2)
        elseif(isw.eq.3 .or. isw.eq.6) then
          call rvesr(hr(np(260)),hr(np(261)),mr(np(270)), isw)

c         Form remaining micro-scale contributions

          call pforma(ul,xl,tl,ld,p,s,ie,d,id,x,ix,rben,f,t,jp,
     &                u,ud,b,a,al,ndd,nie,ndf,ndm,nen1,nst,nneq,
     &                prope,arotflg,aufl,bfl,alfl,dfl,efl,gfl,rfl,
     &                .true.,.false.,jsw,ksw,nn1,nn2,nn3)
        endif
      endif

      end
