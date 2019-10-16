c$Id:$
      subroutine bplot(ct)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove arg 3 and arg 4 from call to perspz       09/01/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plot stresses on surface of beam element cross sections.

c      Inputs:
c         ct(*) - data indicators

c      Outputs:
c         To graphics screen
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'elpdat.h'
      include  'fdata.h'
      include  'hdatam.h'
      include  'iofile.h'
      include  'pdata2.h'
      include  'pdatay.h'
      include  'plflag.h'
      include  'ppers.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   labl, setvar, palloc, fl9sv
      real*8    ct(3),pmsh

      integer   nqy,nqz,jmax,iret,bndm,bndf,bnen,bnen1,bcmp
      integer   snumnp,bnumnp,bnumel,nfacold,stype

      integer   bln(2)

      save

c     Set size of surface mesh parameters

      elplt(1) = ct(1)
      stype    = nint(hr(np(25)+99))
      jmax     = 0

c     No shape specified

      if(stype .eq. 0) then

        write(*,3000)

c     Tubes

      elseif(stype .eq. 1) then

        write(*,*) ' *ERROR* Cannot plot for tubes'

c     Rectangles

      elseif(stype .eq. 2) then

        nqz  = nint(hr(np(25)+105))
        nqy  = nqz/10
        nqz  = mod(nqz,10)
        jmax = 2*(nqy + nqz) - 4  ! Hack for one block

        write(iow,2000) nqy,nqz
        if(ior.lt.0) then
          write(*,2000) nqy,nqz
        endif

c     Shaped Sections

      elseif(stype .eq. 3) then ! Wide flange
        jmax = 12
      elseif(stype .eq. 4) then ! Channel
        jmax = 10
      elseif(stype .eq. 5) then ! Angle
        jmax = 7
      elseif(stype .eq. 6) then ! Solid circle
        jmax = 8
      endif ! stype

      if(jmax.gt.0) then
        bndm   = 3
        bndf   = 1
        bnen   = 4
        bnen1  = 5
        bcmp   = 1
        bnumnp = jmax*numnp
        bnumel = jmax*numel

c       Resize the arrays for edge projections if necessary

        if(np(64).ne.0) then
          setvar = palloc( 64,'TEFAC',bnumnp+1,            1)
          setvar = palloc( 65,'TENRM',max(bnumel,bnumnp)*3,2)
          setvar = palloc( 63,'TECON',0, 1)  ! delete current storage
        endif
        if(np(199).ne.0) then
          setvar = palloc( 66,'VISN ',  bnumnp,1)
          setvar = palloc(199,'MXBEA',5*bnumel,1)
          setvar = palloc(200,'XBEAM',3*bnumnp,2)
          setvar = palloc(201,'SBEAM',  bnumnp,2)
          setvar = palloc(202,'WBEAM',   numnp,2)
        else
          call pzeroi(mr(np(199)),5*bnumel)
          call pzero (hr(np(200)),3*bnumnp)
          call pzero (hr(np(201)),  bnumnp)
          call pzero (hr(np(202)),   numnp)
        endif
        call pzeroi(mr(np(66)),  bnumnp)

        bln(1) = 0
        bln(2) = 1

c       Compute projected stresses

        hflgu  = .false.
        h3flgu = .false.
        fl9sv  = fl(9)
        fl(9)  = .false.   ! set to force plot on t_n+1 state
        call formfe(npuu,np(26),np(26),np(26),
     &             .false.,.false.,.false.,.false.,20,1,numel,1)
        fl(9)  = fl9sv

c       Divide by weights

        labl   = .true.
        call bdivwt(hr(np(200)),hr(np(201)),hr(np(202)),numnp,jmax)

        nfacold = nfac(1)
        setvar = palloc(113,'TEMP3',8*bnumel  ,1)
        call psetip(mr(np(113)),bnumel)

        if(ct(2).le.0.0d0) then
          pmsh = -1.d0
        else
          pmsh =  1.d0
        endif
        call plface(mr(np(199)),mr(np(113)),hr(np(200)),bndm,bnen1,
     &              bnumnp,bnumel,bln,pmsh)

c       Do hidden surface removal and edge definition

        if(kpers.gt.0) then
          ipb = nint(ct(3))
          call perspz(hr(np(200)),mr(np(199)),mr(np(113)),
     &                bnen1,bnen,bndm,bnumnp,bnumel)
          snumnp = numnp
          numnp  = bnumnp
          call p3edge(mr(np(199)),hr(np(200)),3,bnumel,bnen1)
          numnp = snumnp
          edgfl = .true.
        endif

c       Contours of dependent variables

        call bprint(hr(np(201)),bnumnp,bndf,iret)
        call bpltcon(hr(np(200)),mr(np(199)),mr(np(113)),hr(np(201)),
     &               bnumnp,bndm,bndf,bnen,bnen1,bcmp,labl)

        setvar = palloc(113,'TEMP3',0,1)
        nfac(1) = nfacold
      endif ! jmax > 0

c     Formats

2000  format('  Storage for surface mesh assigned:',
     &       ' nqy = ',i3,' nqz = ',i3)

3000  format(' *ERROR* BPLOT: Cannot plot for CROSS SECTION option')

      end
