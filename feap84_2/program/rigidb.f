c$Id:$
      subroutine rigidb (rsw,isw,flg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Rigid body description options

c      Inputs:
c        rsw    - Switch parameter for main level operation
c        isw    - Switch parameter for sub-level operation
c        flg    - Dynamic flag for updates

c      Outputs:
c        flg    - Flag for return of action/error
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_tanfl.h'
      include  'cdata.h'
      include  'crotas.h'
      include  'evdata.h'
      include  'modreg.h'
      include  'part0.h'
      include  'pointer.h'
      include  'rigid1.h'
      include  'rigid2.h'
      include  'rjoint.h'
      include  'sdata.h'
      include  'tdatb.h'
      include  'comblk.h'

      logical   flg, setvar,palloc
      integer   rsw,isw, i

      save

c     Initialize for start of problem

      if    (rsw.eq.0) then

        rbody  = .false.
        nmbody = 0
        nrbody = 0
        nrbdof = 0
        numjts = 0
        nrlds  = 0

c     ISW = 1: Input rigid body data and generate inertial properties

c     Create list of rigid material sets
c     [rigid,nrbdof,npart]  - rigid body specification,dof/rb,partition
c     [mat,density  - repeat for each rigid material]

c     ISW = 2: Create list of joints
c     [joint]  - joint specifications
c     [type,node1,node2,body1,body2] - type data
c       type = 'ball'&socket - needs node1,node2
c       type = 'revo'lute    - needs node1,node2,body1,body2

c     ISW = 3: Rigid body loads
c     [rloa]  - Rigid load specifications
c     [body,1-comp,2-comp,etc.]

c     ISW = 4: Rigid body boundary conditions
c     [rbou]  - Rigid b.c. specifications
c     [body,1-comp,2-comp,etc.]

      elseif(rsw.eq.1) then

        call rinput(isw,flg)

c     Form arrays for fem + rigid body

      elseif(rsw.eq.2) then

        if(rbody) then
          call pformr(mr(np(34)),hr(np(35)),hr(np(36)),  hr(np(43)),
     &                mr(np(99)),mr(np(20+npart)),hr(np(30)),hr(rpu),
     &                hr(rpb),hr(rpa),hr(rpl),uafl,dbfl,lafl,ddfl,
     &                ndm,ndf,nrbdof,isw)
        endif

c     Update A:

      elseif(rsw.eq.3) then

        if( rbody .and. npart.eq.nrbprt ) then

c       Get coordinates in reference configuration

          call rupdat(hr(np(26)),hr(np(95)),hr(np(108)),hr(np(104)),
     &                mr(np(99)),mr(np(96)),ndm,flg,isw)

          call uprigb(hr(np(40)),hr(np(42)),hr(np(43)),mr(np(100)),
     &                hr(np(95)),hr(np(104)),hr(np(108)),mr(np(96)),
     &                ndm,ndf,numnp,isw)
        endif

c     Update B:

      elseif(rsw.eq.4) then


c       Small deformation Master-Slave projections

        if(np(167).ne.0) then

          call ruplnk(hr(np(40)),hr(np(42)),hr(np(43)),mr(np(100)),
     &                mr(np(167)),3,ndm,ndf,numnp,flg)

c       Large deformation rigid bodies

        elseif( rbody .and. npart.eq.nrbprt ) then

          if(isw.ne.2) then
            call rupdat(hr(np(26)),hr(np(95)),hr(np(108)),hr(np(104)),
     &                  mr(np(99)),mr(np(96)),ndm,flg,isw)

            call uprigb(hr(np(40)),hr(np(42)),hr(np(43)),mr(np(100)),
     &                  hr(np(95)),hr(np(104)),hr(np(108)),mr(np(96)),
     &                  ndm,ndf,numnp,isw)
          endif

          call rupjnt(hr(np(26)),mr(np(101)),hr(np(103)),cc1,numjts,
     &                5,isw)

          if (frotas) then
            call rrupdt(hr(np(26)),hr(np(40)),mr(np(100)),
     &                  mr(np(99)),numnp,ndm,ndf)
          endif
        endif

c       Complete solution for modal unknowns

        if(nmbody.gt.0 .and. mf.gt.0) then

          do i = 1,nmbody

            if(isw.eq.1) then
              call dsolmod(hr(np(177)),hr(np(178)),hr(np(179)),
     &                     hr(np(180)),mf,3,.false.,.true.,.false.,i)
            elseif(isw.eq.2) then
              call dsolmod(hr(np(177)),hr(np(178)),hr(np(179)),
     &                     hr(np(180)),mf,3,.false.,.false.,.true.,i)
              setvar = palloc(183,'UMOD',neqmf*8,2)
              call upmodal(hr(np(40)),hr(np(77)),hr(np(180)),
     &                     hr(np(183)),mr(np(100)),mr(np(176)),
     &                     hr(np(104)),ndm,ndf,numnp,i,hr(np(42)))
            endif

          end do ! i
        endif

      endif

      end
