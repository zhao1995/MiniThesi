c$Id:$
      subroutine pmacr1(lct,ct,j)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove interupt tests (intf and lintr)           17/11/2006
c       2. Add call to ckfixed for fixed dof checks         06/12/2006
c       3. Increase size of tdatabuf to 6 for buffering     16/03/2007
c       4. Add 'form,stre' option for plots                 29/03/2007
c       5. Exchange energy for reactions                    14/04/2007
c       6. Set length on assignment to cmtyp                17/04/2007
c       7. Add 'nummat' to 'rstprf' call                    21/07/2007
c       8. Add 'reac,nopr' for no output option             17/11/2007
c       9. Comment last two outputs for 'reac file'         19/02/2009
c      10. Add 'rffl' flag for file reaction outputs        22/02/2009
c      11. Avoid file output of reactions for ct(1,l) > 0   03/03/2009
c      12. Set 'pltmfl' to true for reaction outputs        22/03/2009
c      13. Separate 'id' & 'eq' on call to pload            27/04/2009
c      14. Change 'Maximum' to 'Initial' in format 2004     02/06/2009
c      15. Change pointer allocation for mass terms         24/07/2009
c      16. Add 'abs' to convergence test on residual        26/10/2009
c      17. Add output of convergence values to log file     16/02/2010
c      18. Add default to 'imaginary' prints (cplxpr)       09/08/2010
c      19. Add 'resnm' to 'rdata.h' for residual norm ck.   10/08/2010
c      20  Correct set of 'xtol' for 'stre coor' option     20/08/2010
c      21. Set length of file extension to 8                25/01/2011
c      22. Add 'hill'-mandel computation                    19/04/2011
c      23. Add set of storage for HSELM and HDNP hist plots 05/01/2012
c      24. Use dumar(1) instead of scalars                  09/01/2012
c      25. Set history plot/node array to zero              10/01/2012
c      26. Add numerical tangent for contact(223)           18/01/2013
c      27. Add stre on/off to allow quadrature stess output 14/03/2013
c      28. Add coords and ndm on call to prtrel             12/04/2013
c          Add displacement to prtrea and prtrel
c      29. Use 'finflg' to set gradu for hill-mandel        07/05/2013
c      30. Change sign on reactions to file                 29/11/2013
c      31. Save energy during Hill-Mandel solution          01/01/2014
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Command language instruction subprogram.  Mostly for
c               FEM arrays

c      Inputs:
c         lct(*)     - Command option
c         ct(3,*)    - Command parameters
c         j          - Command number in this routine

c      Outputs:
c         Depends on command number j
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'adapt1.h'
      include   'adapt2.h'
      include   'allotd.h'
      include   'arclel.h'
      include   'arclei.h'
      include   'arcler.h'
      include   'augdat.h'
      include   'auto1.h'
      include   'auto2.h'
      include   'cdata.h'
      include   'cdat1.h'
      include   'comfil.h'
      include   'compac.h'
      include   'compas.h'
      include   'complx.h'
      include   'counts.h'
      include   'c_values.h'
      include   'ddata.h'
      include   'elacts.h'
      include   'elcapt.h'
      include   'eldata.h'
      include   'eldatp.h'
      include   'elpers.h'
      include   'eltran.h'
      include   'endata.h'
      include   'eqslv.h'
      include   'eqsym.h'
      include   'errchk.h'
      include   'errind.h'
      include   'evdata.h'
      include   'fdata.h'
      include   'gltran.h'
      include   'hdatam.h'
      include   'idptr.h'
      include   'iodata.h'
      include   'iofile.h'
      include   'ldata.h'
      include   'modreg.h'
      include   'mxsiz.h'
      include   'ndata.h'
      include   'part0.h'
      include   'part1.h'
      include   'pbody.h'
      include   'pdata0.h'
      include   'pdata3.h'
      include   'pfeapb.h'
      include   'plist.h'
      include   'pmod2d.h'
      include   'pointer.h'
      include   'print.h'
      include   'prflag.h'
      include   'prlod.h'
      include   'prstrs.h'
      include   'pscal.h'
      include   'ptest.h'
      include   'rdata.h'
      include   'rdat0.h'
      include   'rdat1.h'
      include   'rigid1.h'
      include   'rjoint.h'
      include   'sdata.h'
      include   'setups.h'
      include   'ssolve.h'
      include   'strnum.h'
      include   'tdata.h'
      include   'xtout.h'
      include   'comblk.h'

      include   'p_formfe.h'
      include   'p_int.h'

      logical    fa,tr,cfr,pcomp,f8o,pcfl,exst,tfl,setvar,palloc,compsv
      logical    nomass,factor, cknon0, solvsv, unsfl, naninfck, hform
      logical    rffl, dtfl(19)
      character  tname*5, fext*8, fnamr*132, cmtyp*4, lct(*)*15
      integer    iops,mops,i,j, elist, neqms, kk,k1,k2,k3,k4,k5
      integer    gneq, toteq
      real*4     tary(2), etime , tt
      real*8     tops,ur,ee,etab,reln,step,propsv,aengysv, ct(3,*)
      real*8     dot,mnorm, rnaugm, tdatabuf(3,2),rnorms(3),dumar(1)
      real*8     td(10)

      save

      data       fa,tr/.false.,.true./ , pcfl/.true./

      data       dtfl / .true., .true., .true., .true., .true., .true.,
     &                 .false., .true., .true., .true., .true.,.false.,
     &                 .false., .true., .true., .true., .true., .true.,
     &                  .true./

c     Set history update flag to false (no updates)

      hflgu  = .false.
      h3flgu = .false.
      nomass = .false.

c     Set transient parameters for current

      if(fl(9)) then

c       Allow zero dt for initial acceleration computation

        if(.not.dtfl(j) .or. (j.eq.4 .and. pcomp(lct(l),'acce',4))) then
          call dsetci(.false.)
        else
          call dsetci(.true.)
        endif
      endif
      do i = 1,3
        ctan(i) = gtan(i)
      end do ! i

c     Transfer to correct process

      go to (1,2,2,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     &       4,24),j
c-----[--.----+----.----+----.-----------------------------------------]

c     Print stress values

c     [stre,,k1,k2,k3]        - output elmt values for k1 to k2 inc k3
c     [stre,all]              - output all element values
c     [stre,node,k1,k2,k3]    - output nodal stresses k1 to k2 inc k3
c                               (lumped projection)
c     [stre,gnod,k1,k2,k3]    - output global node stress k1-k2 inc k3
c                               (lumped projection)
c     [stre,cnod,k1,k2,k3]    - output nodal stresses k1 to k2 inc k3
c                               (consistent projection)
c     [stre,coor,nxt,xt,xtol] - output nodal stresses at x-nxt=xt+-xtol
c     [stre,cerr,k1,k2,k3]    - output nodal stresses k1 to k2 inc k3
c                               (consistent projection with errors)
c     [stre,erro,k1,k2,k3]    - output nodal stresses k1 to k2 inc k3
c                               (lumped projection with errors)
c     [stre,cont,k1,k2,k3]    - output contact values for pair k1
c                               from el k2 to el k3
c     [stre,<on,off>]         - turn on/off output at quadrature points

1     k1 = nint(ct(1,l))
      k2 = nint(ct(2,l))
      if(k2.eq.0) k2 = k1
      k3 = nint(ct(3,l))
      if(k3.eq.0) k3 = 1

      if(pcomp(lct(l),'off',3)) then
        qoutfl = .false.
        return
      elseif(pcomp(lct(l),'on',2)) then
        qoutfl = .true.
        return
      endif

c     Set flag for print using local numbers

      pfeap_gnod = .false.

      nxt = 0
      if (pcomp(lct(l),'node',4)  .or.
     &    pcomp(lct(l),'erro',4)  .or.
     &    pcomp(lct(l),'cerr',4)  .or.
     &    pcomp(lct(l),'cnod',4)  .or.
     &    pcomp(lct(l),'gnod',4)  .or.
     &    pcomp(lct(l),'coor',4)) then

c       Allocate array storage

        if (plfl) then
          setvar = palloc ( 58,'NDNP',numnp*npstr,2)
          setvar = palloc ( 57,'NDER',numnp*8    ,2)
          setvar = palloc ( 60,'NDNS',max(nen*npstr,nst*nst),2)
          setvar = palloc (207,'NSCR',numel      ,2)
          nper   = np(57)
          npnp   = np(58)
          plfl   = .false.
          if(histpltfl) then
            setvar = palloc(304,'HSELM',nen*hplmax  ,2)
            setvar = palloc(305,'HDNP ',numnp*hplmax,2)
          endif
        endif
        nph  = npnp
        ner  = nper

c       Set output limits

        if (pcomp(lct(l),'node',4)  .or.
     &      pcomp(lct(l),'erro',4)  .or.
     &      pcomp(lct(l),'cerr',4)  .or.
     &      pcomp(lct(l),'cnod',4)) then
          k1 = max(1,min(numnp,k1))
          k2 = max(1,min(numnp,k2))
          if(k2-k1.ne.0) k3 = sign(k3,k2-k1)
        elseif(pcomp(lct(l),'coor',4)) then
          nxt = max(1,min(k1,ndm))
          xt  = ct(2,l)
          if(ct(3,l).eq.0.0d0) then
            xtol = 1.0d-6*dxmsh(nxt)
          else
            xtol = abs(ct(3,l))
          endif
          k1 = 1
          k2 = numnp
          k3 = 1
        endif

c       Consistent nodal stress projection: Determine profile

        if ( (pcfl .and. pcomp(lct(l),'cnod',4)) .or.
     &       (pcfl .and. pcomp(lct(l),'cerr',4))      ) then
          setvar = palloc(197,'NS1  ',numnp,1)
          setvar = palloc(111,'TEMP1',numnp,1)
          do i = 1,numnp
            mr(np(111)+i-1) = i
          enddo ! i
          call rstprf(mr(np(197)),mr(np(34)),mr(np(111)),mr(np(33)),
     &                mr(np(32)),mr(np(99)),mr(np(100)),mr(np(101)),
     &                1,nen+4,nen,numnp,numnp,numel,nummat)
          setvar = palloc(111,'TEMP1', 0,1)
          call nwprof(mr(np(197)),numnp)
          setvar = palloc(198,'NS2  ',mr(np(197)+numnp-1)+numnp,2)
        endif

c       Error indicator projections to nodes

        fp(6) = npnp + numnp
        if(.not.fl(11) .or. pcomp(lct(l),'erro',4) .or.
     &                      pcomp(lct(l),'cerr',4) ) then
          if(pcomp(lct(l),'erro',4).or.pcomp(lct(l),'cerr',4)) then
            trifl = .true.
          else
            call pzero (hr(np(207)), numel)
          end if
          istv = npstr - 1
          call pzero (hr(npnp), npstr*numnp)
          call pzero (hr(nper),     8*numnp)
          if(histpltfl) then
            call pzero (hr(np(305)), numnp*hplmax)
          endif

c         Initialize caption array

          do i = 1,50
            ecapt(i) = '  '
          end do ! i

c         Consistent projection solution

          if(pcomp(lct(l),'cnod',4) .or.pcomp(lct(l),'cerr',4)) then
            fp(5)  = np(36)
            np(36) = np(60)
            call pzero (hr(np(198)),mr(np(197)+numnp-1)+numnp)
            call formfe(np(40),np(26),np(198),np(26),tr,fa,fa,fa,8,
     &                  1,numel,1)
            np(36) = fp(5)

            call datri(hr(np(198)+numnp),hr(np(198)+numnp),hr(np(198)),
     &                 mr(np(197)),numnp,numnp)
            do  i = 1,istv
              fp(5) = npnp + i*numnp
              if( cknon0(hr(fp(5)),numnp) ) then
                call dasol(hr(np(198)+numnp),hr(np(198)+numnp),
     &                     hr(np(198)),hr(fp(5)),mr(np(197)),
     &                     numnp,numnp,aengy,.false.)
              endif
            enddo ! i
            call pltstr(hr(npnp),hr(nper+numnp),hr(fp(6)),
     &                  numnp,ndm,fa)

c         Lumped projection solution

          else
            fp(5)  = np(36)
            np(36) = np(60)
            call formfe(np(40),np(26),np(26),np(26),fa,fa,fa,fa,8,
     &                  1,numel,1)
            np(36) = fp(5)
            call pltstr(hr(npnp),hr(nper+numnp),hr(fp(6)),
     &                  numnp,ndm,tr)
          endif
        endif

c       Set flag for global node prints

        if(pcomp(lct(l),'gnod',4)) then
          pfeap_gnod = .true.
        endif

c       Output sets

        call prtstr(hr(np(43)),hr(nper+numnp),hr(fp(6)),
     &              ndm,numnp,k1,k2,k3,prth)
        fl(11) = .true.

c     Output values as specified by contact elements

      elseif (pcomp(lct(l),'cont',4)) then
          cvaluei(1) = max(1,nint(ct(1,l)))
          cvalue (1) = ct(2,l)
          cvalue (2) = ct(3,l)
          call contact (204)

c     Output values as specified by elements

      else
        if (pcomp(lct(l),'all ',4)) then
          k1 = 1
          k2 = numel
          k3 = 1
        else
          k1 = max(1,min(numel,k1))
          k2 = max(1,min(numel,k2))
          if(k2-k1.ne.0) k3 = sign(k3,k2-k1)
        endif
        call formfe(np(40),np(26),np(26),np(26),fa,fa,fa,fa,4,k1,k2,k3)

      endif

      return

c     Form tangent stiffness: j= 2 (utan); j=3 (tang)
c     [utan]                     - form unsymmetric tangent
c     [tang]                     - form symmetric tangent

c     [utan,,1]                  -   " + form rhs and solve.
c     [tang,,1]                  -   " + form rhs and solve.

c     [utan,line,1,shift,value]  -   " with line search on value
c     [tang,line,1,shift,value]  -   " with line search on value

c     [tang,eigv,0,shift]        -   " with no mass/damping added
c     [tang,nume,0,shift]        - Form numerical differentiated tangent

2     castif = .true.
      camass = .false.
      cadamp = .false.
      compre = .false.
      msplt  =  0
      call contact(103)
      if( pcomp(lct(l),'nume',4) ) then
        write(iow,*)' Numerically differentiated tangent used'
        if(ior.lt.0) then
          write(*,*)' Numerically differentiated tangent used'
        endif
        ndflg = .true.
      else
        ndflg = .false.
      endif

c     Auto time step level set

      lvauto  = lv
      lvautoi = lv

c     Presolution of equations: Set sparse storage data

      if(neq.gt.0) then
        cfr = j.eq.2
        call presol(cfr, exst)
        if(exst) return
        jcmplx = kcmplx

c       Form residual for ct(l,1) > 0

        f8o = .false.
        if(ct(1,l).gt.0.0d0) then
          f8o    = .true.
          itract = itract + 1
          itrdea = itrdea + 1
          compre = .true.

c         Form interpolated load vector - include nodal/surface parts

          call ploa1(ttim,dt)
c         call pload(    id   ,    u    ,    ftn  ,     dr  ,
          call pload(mr(np(31)+nneq),mr(id31),hr(np(40)),hr(np(30)),
     &               hr(np(26)),prop*rlnew,tr,tr)
        endif

c       Construct shifts

        rfl   = .false.
        shflg = .false.
        shift = ct(2,l)
        if( pcomp(lct(l),'eigv',4)
     &      .or. (.not.fl(9).and.shift.ne.0.0d0)) then
          ctan(2) = 0.0d0
          if(prnt) write(iow,2000) shift
          if(ior.lt.0.and.prnt) write(*,2000) shift
          if(idenf) then
            call colred(hr(np(nx)),shift,neq, hr(na))
          else
            ctan(2) =  0.0d0
            ctan(3) = -shift*ctan(1)
            shflg   = .true.
          endif

c       Specified shifts or dynamics/static cases

        elseif(shift.ne.0.0d0) then
          ctan(2) =  0.0d0
          ctan(3) = -shift*ctan(1)
          shflg   = .true.
        else
          shift   = -ctan(3)
        endif

c       Add modal solution parts

        if(nmbody.gt.0 .and. mf.gt.0) then

          do i = 1,nmbody

            setvar = palloc(117,'TEMP7', ndf*numnp, 2)

            call pmodfc(hr(np(76)),hr(np(77)), mr(np(100)),mr(np(176)),
     &                  hr(np(104)),hr(np(117)),hr(np(30)),mf,ndf,
     &                  ndm,numnp,hr(np(177)),hr(np(178)),hr(np(179)),
     &                  hr(np(180)),i,mnorm)

            setvar = palloc(117,'TEMP7', 0, 2)

            call dsolmod(hr(np(177)),hr(np(178)),hr(np(179)),
     &             hr(np(180)),mf,3,tr,fa,fa,i)

          end do ! i

        else
          mnorm = 0.0d0
        endif

c       Compute flexible and rigid FE contributions to arrays

        hflgu  = f8o
        h3flgu = f8o
        call formfe(np(40),np(26),na,nal,tr,f8o,cfr,fa,3,1,numel,1)
        ndflg = .false.

c       Output residual norm

        tt = etime(tary)
        if(f8o) then
          if(pfeap_on) then
            if(pfeap_blk) then
              rnorm  = dot(hr(np(26)),hr(np(26)),ndf*numpn) + mnorm
            else
              rnorm  = dot(hr(np(26)),hr(np(26)),neq) + mnorm
            endif
            rnorms(1) = rnorm
            rnorms(2) = rnorm1
            rnorms(3) = rnormn
            call pfeapsr(rnorms,tdatabuf, 3,.true.)
            rnorm  = rnorms(1)
            rnorm1 = rnorms(2)
            rnormn = rnorms(3)
            toteq  = numteq
          else
            rnorm  = dot(hr(np(26)),hr(np(26)),neq) + mnorm
            toteq  = neq
          endif
          dumar(1) = rnorm
          if(naninfck(dumar,1,1)) then
            write(  *,*) ' *ERROR* Residual norm is NaN or Inf'
            write(ilg,*) ' *ERROR* Residual norm is NaN or Inf'
            write(iow,*) ' *ERROR* Residual norm is NaN or Inf'
            call plstop()
          endif
          rnorm  = sqrt(rnorm)
          gneq   = nneq
          compre = .false.
          if(rnmax.eq.0.0d0) then
            reln = 1.d0
            rel0 = rnorm
          else
            if(rel0.eq.0.d0) rel0 = 1.d0
            reln = rnorm/rel0
          endif
          resnm = reln  ! Norm for residual convergence
          if(prnt) write(iow,2001) rnorm,reln,tary
          if((ior.lt.0.and.prnt).or.echo) write(*,2001) rnorm,reln,tary
          fl(7) = .false.
          fl(8) = .false.
          if(nmeth.eq.2) then
            ee = aengy
          else
            ee = 0.0d0
          endif
          if(abs(aengy).lt.aold) aold = abs(aengy)
        else
          if(prnt) write(iow,2002) tary
          if((ior.lt.0.and.prnt) .or. echo) write(*,2002) tary
        endif

c       Set pointers then factor and solve equations

        fp(1)  = na
        fp(2)  = nau
        fp(3)  = nal
        fp(4)  = np(20+npart)
        factor = ct(1,l).ge.0.0d0

        call psolve(ittyp,hr(np(26)),fp,factor,f8o,cfr, prnt)

        if(factor) then

c         Timing for factorization

          if(tdiff .gt.0.0d0 .and. prnt .and. pfr .and.
     &            (ittyp.eq.-3 .or. ittyp.eq.-1)) then
            call datric(iops,mops,mr(fp(4)),neq)
            tops  = (dble(mops) + dble(iops)*1.d-6)/tdiff
            if(cfr   ) tops = 2.d0*tops
            if(cplxfl) tops = 4.d0*tops
            if(ior.lt.0) then
              write(*,2019) tops,tdiff
            endif
            write(iow,2019) tops,tdiff
          endif

        else
          write(iow,3002)
          if(ior.lt.0) write(*,3002)
        endif

c       Update solutions using line search or arc length method

        if(f8o) then

c         Update iteration counter

          niter = niter + 1
          iaugm = iaugm + 1
          lvsol = lv

c         Set parameter for augmented convergence test

          if(lvaug.gt.0) then
            if(naugm.eq.1) then
              if(rnmax.eq.0.0d0) then
                rnaugm = abs(aengy)
              else
                rnaugm = rnmax
              endif
            endif
            augmfl = .false.
            if(abs(aengy).lt.tol*rnaugm .and.
     &         nint(ct(1,lve(lv))).eq.1) then
              augmfl = .true.
            endif
          endif

c         Set initial values for a line search (conservative to
c         check initital iterate solution for possible line search)

          if(ct(3,l).le.0.0d0) ct(3,l) = 0.8d0
          if (rnmax.eq.0.0d0) then
            reln  = 1.d0
            rnmax = abs(aengy)
            ee    = 0.0d0
            aold  = rnmax/0.9999d0/ct(3,l)
            autcnv = .false.
            autr0 = rnorm
            do i = 1,10
              autr(i) = 0.0d0
            end do ! i
          else
            reln  = (aengy-ee)/rnmax
          endif
          do i = 9,1,-1
            autr(i+1) = autr(i)
          enddo ! i
          autr(1) = rnorm
          if(pfr) then
            write(iow,2004) rnmax,aengy-ee,reln,tol
            if(ior.lt.0) then
              write(*,2004) rnmax,aengy-ee,reln,tol
            endif
          endif

c         Convergence checks on energy and residual norms

          if(abs(aengy-ee).le.tol*rnmax .or. linear .or.   ! Energy norm
     &       abs(rnorm)/dble(toteq) .lt.                   ! Force meas.
     &       abs(rnorm1/rnormn)*sqrt(tol)*1.d-3     .or.
     &       resnm.lt.tolr                          .or.   ! Resid. norm
     &       abs(aengy-ee).lt.enzer                 ) then ! Zero energy

            if(lv.gt.1) then
              ct(1,lve(lv)) = ct(1,lvs(lv))
              l             = lve(lv) - 1
              floop(1)      = .false.
            endif
            autcnv = .true.

c         Scale back step and update the solution -- no line search

          elseif(testfl .and. niter.le.1) then

            do i = 0,nneq-1
              hr(np(26)+i) = hr(np(26)+i)*testva
            end do ! i

c         Perform line search and/or arclength updates

          elseif(.not.linear) then

c           Line search - linear search along solution direction

            if(pcomp(lct(l),'line',4)) then
              if(abs(aengy-ee).gt.ct(3,l)*aold) then
                setvar = palloc(111,'TEMP1',nneq        , 2)
                setvar = palloc(112,'TEMP2',nneq*(nrt+4), 2)

                step = 1.0d0
                call serchl(aold,mr(id31),np(111),np(40),
     &                      hr(np(26)),ct(3,l),hr(np(112)),gneq,step)

                setvar = palloc(112,'TEMP2',0, 2)
                setvar = palloc(111,'TEMP1',0, 2)
              endif
            endif
          endif

c         Perform arc length control

          if(arcf) then
            call arclen(hr(np(40)),hr(np(26)),hr(np(84)),hr(np(85)),
     &                  hr(np(27)),mr(id31),ndf,ttim)
          endif

c         Check solution method: subtract increment from new increment

          if(nmeth.eq.2) then
            do i = 0,nneq-1
              j = mr(id31+i) - 1
              if(j.ge.0) then
                hr(np(26)+j) = hr(np(26)+j) - hr(np(41)+nneq+i)
              endif
            end do ! i
          endif

c         Update solution

          call update(mr(id31),hr(np(30)),hr(np(40)),hr(np(42)),
     &                hr(np(26)),fl(9),2)

c         Output convergence information

          if(convfl) then
            write(ilg,2030) nstep,niter,rnorm,rnorm/rel0,
     &                                  aengy,aengy/rnmax
          endif

        endif

c     No active equation case

      else

        if(prnt .and. ior.lt.0) write(*,3006)
        call ploa1( ttim,dt)
        call pload (mr(np(31)+nneq),mr(id31),hr(np(40)),hr(np(30)),
     &              hr(np(26)),prop*rlnew,fa,fa)
        call update(mr(id31),hr(np(30)),hr(np(40)),hr(np(42)),
     &              hr(np(26)),fl(9),2)

c       Update stresses

        hflgu  = .true.
        h3flgu = .true.
        rnmax  = 1.d0
        niter  = niter  + 1
        call formfe(np(40),np(26),na,nal,fa,fa,fa,fa,6,1,numel,1)

c       Convergence set

        if(lv.gt.1) then
          ct(1,lve(lv)) = ct(1,lvs(lv))
          l             = lve(lv) - 1
          floop(1)      = .false.
        endif
        autcnv = .true.

      endif
      return

c     Compute residual for time step/iteration

c     [form]       - form rhs residual
c     [form,acce]  -    " + get initial acceleration if needed
c     [form,expl]  -    " + do explicit solution with lumped mass
c     [form,conv]  -    " + check residual for convergence.
c     [form,stre]  -    " + update history only
c     [resid]      - Form residual only

4     hform = .false.
      if(pcomp(lct(l),'stre',4)) then
        hform = .true.
      elseif(fl(8) .and. j.eq.4) then
        if(ior.lt.0) then
          write(*,3007)
        else
          write(iow,3007)
          write(ilg,3007)
          call plstop()
        endif
        return
      endif
      compre = .true.
      if(j.eq.4) then
        hflgu  = .true.
        h3flgu = .true.
        itract = itract + 1
        itrdea = itrdea + 1
      endif
      rfl    = .false.
      msplt  =  0

c     Compute current residual for loads

      call ploa1(ttim,dt)
      call pload(mr(np(31)+nneq),mr(id31),hr(np(40)),hr(np(30)),
     &           hr(np(26)),prop*rlnew,tr,fa)

c     Compute initial acceleration

      if(j.eq.4 .and. neq.gt.0 .and. pcomp(lct(l),'acce',4)) then

        if(fl(9)) then
          if(ior.lt.0) write(*,*) ' Forming initial acceleration'

c         Explicit integrator form uses lumped mass

          if(noi.eq.4 .or. noi.eq.8) then
            if(np(npart+12).ne.0) then
              fp(1) = np(npart+12)
              call formfe(np(40),np(26),na,nal,fa,tr,fa,fa,6,1,numel,1)
              call psolve(0,hr(np(26)),fp,.false.,.true.,.false., prnt)
            else
              if(ior.gt.0) then
                write(iow,3011)
                write(ilg,3011)
                call plstop()
              endif
              write(*,3011)
            endif

c         Implicit integrator form

          else
            castif = .true.
            camass = .false.
            cadamp = .false.
            compre = .false.
            msplt  =  0
            ndflg  = .false.
            cfr    = .false.
            call presol(cfr, exst)
            if(exst) return
            jcmplx = kcmplx

c           Form residual for initial acceleration

            ctan(1) = 0.0d0
            ctan(2) = 0.0d0
            ctan(3) = 1.0d0
            call formfe(np(40),np(26),na,nal,tr,tr,fa,fa,3,1,numel,1)

c           Set pointers then factor and solve equations

            fp(1)  = na
            fp(2)  = nau
            fp(3)  = nal
            fp(4)  = np(20+npart)
            factor = .true.

            call psolve(ittyp,hr(np(26)),fp,factor,.true.,.false., prnt)
          endif
          call pmove (hr(np(26)),hr(np(42)+nneq),neq)
          call pexpd(hr(np(42)+nneq),hr(np(26)),mr(id31),ndf,nneq)

c       Write error

        else
          write(ilg,3000)
          write(iow,3000)
          if(ior.lt.0) write(*,3000)
        endif

c     Form residual for FE arrays and contact

      else

        call formfe(np(40),np(26),np(26),np(26),fa,tr,fa,fa,6,1,numel,1)
        nform = nform + 1
        if(hform) return

      endif

      if(pfeap_on) then
        if(pfeap_blk) then
          rnorm  = dot(hr(np(26)),hr(np(26)),ndf*numpn)
        else
          rnorm  = dot(hr(np(26)),hr(np(26)),neq)
        endif
        rnorms(1) = rnorm
        rnorms(2) = rnorm1
        rnorms(3) = rnormn
        call pfeapsr(rnorms,tdatabuf, 3,.true.)
        rnorm  = rnorms(1)
        rnorm1 = rnorms(2)
        rnormn = rnorms(3)
        toteq  = numteq
      else
        rnorm = dot(hr(np(26)),hr(np(26)),neq)
        toteq  = neq
      endif
      dumar(1) = rnorm
      if(naninfck(dumar,1,1)) then
        write(  *,*) ' *ERROR* Residual norm is NaN or Inf'
        write(ilg,*) ' *ERROR* Residual norm is NaN or Inf'
        write(iow,*) ' *ERROR* Residual norm is NaN or Inf'
        call plstop()
      endif
      rnorm  = sqrt(rnorm)

c     Set residual array for auto dt use

      if(rnmax.eq.0.0d0) then
        reln  = 1.d0
        rel0  = rnorm
        if(j.eq.4) then
          autcnv = .false.
          autr0  = rnorm
          do i = 1,10
            autr(i) = 0.0d0
          end do ! i
        endif
      else
        if(rel0.eq.0.d0) rel0 = 1.d0
        reln  = rnorm/rel0
      endif
      if(j.eq.4) then
        do i = 9,1,-1
          autr(i+1) = autr(i)
        end do ! i
        autr(1) = rnorm
      endif

c     Output current residual norm

      if(prnt) then
        write(iow,2001) rnorm,reln
        if(ior.lt.0) then
          write(*,2001) rnorm,reln
        endif
      endif
      if(j.eq.23) return ! Residual form only
      fl(8)  = .true.
      compre = .false.

41    if(pcomp(lct(l),'expl',4) .or. nomass) then

c       Perform solution and update

        if(rnmax.eq.0.0d0) then
          rnmax  = abs(reln)
        endif
        nomass = .false.
        if(fl(2).and.fl(9)) then

c         Solve rigid body equations

          if(rbody) then
            call pformrx(mr(id31),mr(np(100)),ndf,numnp,neqms)
            call rasble (hr(np(43)),hr(np(40)),hr(np(109)),hr(np(26)),
     &                   mr(np(100)),mr(np(99)),hr(np(95)),hr(np(110)),
     &                   mr(np(101)),hr(np(102)),hr(np(97)))
          else
            neqms = neq
          endif

c         Solve flexible equations

          call piacel(hr(nl),hr(np(26)),hr(np(26)),neqms)

c         Update solution parameters

          call update(mr(id31),hr(np(30)),hr(np(40)),hr(np(42)),
     &                hr(np(26)),fl(9),2)
          fl(8) = .false.
          if(abs(reln).le.sqrt(tol)*rnmax .and. lv.gt.1) then
            ct(1,lve(lv)) = ct(1,lvs(lv))
            l             = lve(lv) - 1
            floop(1)      = .false.
            autcnv        = .true.
          endif

c       Write error

        else
          write(iow,3004)
          if(ior.lt.0) write(*,3004)
          nomass = .true.
          cmtyp  = 'lump'
          go to 51
        endif

c     Convergence check on residual

      elseif(pcomp(lct(l),'conv',4) .and. reln .lt. sqrt(tol)) then
        if(lv.gt.1) then
          ct(1,lve(lv)) = ct(1,lvs(lv))
          l             = lve(lv) - 1
          floop(1)      = .false.
        endif
        autcnv = .true.
      endif

      return

c     [mass],lump : Lumped mass matrix used
c     [mass],cons : Consistent  mass matrix used
c     [mass],unsy : Unsymmetric mass matrix used
c     [mass]      : Same as 'cons'istent

c     Form lumped mass approximation

5     imtyp  = 1
      cmtyp  = lct(l)(1:4)
      nomass = .false.
51    compsv =  compfl
      if(pcomp(cmtyp,'lump',4)) then
        imtyp = 1
        fl(1) = .false.
        fl(2) = .true.
        fl(5) = .false.
        call premas(fl,cmtyp,unsfl)
        if(prnt .and. ior.lt.0) write(*,2016)

c     Form consistent mass approximation

      else
        fl(1) = .true.
        fl(2) = .false.
        fl(6) = .false.
        call premas(fl,cmtyp,unsfl)
        if(prnt .and. ior.lt.0) then
          if(imtyp.eq.1) then
            if(unsfl) then
              write(*,2017) 'UNSYMMETRIC'
            else
              write(*,2017) 'CONSISTENT'
            endif
          elseif(imtyp.eq.2) then
            write(*,2018)
          endif
        endif
      endif

      jcmplx =  mcmplx ! Changed from jcmplx
      castif = .false.
      cadamp = .false.
      camass = .true.
      idenf  = .false.
      solvsv =  solver
      solver =  .true.
      pnl = nl
      pna = nm
      pnb = 0
      if(nl.eq.0) then
        pnl = pna
        pnu = nml ! Changed from pnl
      endif
      if(nm.eq.0) then
        pna = pnl
        pnu = pnl
      endif
      call formfe(np(40),pnl,pna,pnu,fl(1),fl(2),unsfl,fa,5,1,numel,1)
      if(rbody.and.fl(2)) then
        call pzero (hr(np(109)),ndf*numnp)
        call formfe(np(40),np(109),pna,pnu,fa,tr,unsfl,tr,5,1,numel,1)
      endif

      compfl =  compsv
      solver =  solvsv

c     Check that mass matrix has non-zero diagonal entries

      if(.not.pfeap_on) then
      if(imtyp.eq.1 .and. .not.cknon0(hr(np(nx)),neq)) then
        if(nfeqs.gt.0) then
          if(ior.lt.0) then
            write(  *,3008)
          else
            write(ilg,3008)
            write(iow,3008)
            call plstop()
          endif
        else
          if(ior.lt.0) then
            write(  *,3009)
          else
            write(iow,3009)
          endif
        endif
      endif
      endif
      if(nomass) then
        nomass = .false.
        go to 41
      endif
      return

c     Compute reactions and print

c     [reac,,k1,k2,k3]        - print reactions at nodes k1 to k2 inc k3
c     [reac,all]              - print all reactions
c     [reac,coor,nxt,xt,xtol] - print reactions at nodes x-nt=xt+-xtol
c     [reac,node,x1,x2,x3]    - print reactions at node (x1,x2,x3)
c     [reac,imag,k1,k2,k3]    - print complex imaginary reactions
c                               at nodes k1 to k2 inc k3
c     [reac,list,k1]          - print reactions in list k1
c     [reac,regi,k1,k2]       - compute reactions for region k1
c                               assign to proportional load  k2
c     [reac,file]             - compute reactions for active regions
c                               save to file: fsav.rea
c     [reac,nopr]             - prevents any prints

6     nxt = 0

c     Set output limits

      rffl = np(45).ne.0     ! Convert to Cartesian form for 'angl'es.
      tfl  = pfr
      if(.not.cplxpr) then
        if(cplxfl) then
          k4   = nneq
        else
          k4   = 0
        endif
      else
        k4   = 0
      endif
      if(.not.pcomp(lct(l),'list',4)) then
        if (pcomp(lct(l),'all ',4)) then
          k1 = 1
          k2 = numnp
          k3 = 1
        elseif(pcomp(lct(l),'coor',4)) then
          k1 = nint(abs(ct(1,l)))
          nxt = max(1,min(k1,ndm))
          xt  = ct(2,l)
          if(ct(3,l).eq.0.0d0) then
            xtol = 1.0d-6*dxmsh(nxt)
          else
            xtol = abs(ct(3,l))
          endif
          k1 = 1
          k2 = numnp
          k3 = 1
        elseif(pcomp(lct(l),'node',4)) then
          call pclnod(ct(1,l),hr(np(43)),mr(np(190)),ndm,numnp, k1)
          k2 = k1
          k3 = 1
        elseif(pcomp(lct(l),'file',4)) then
          k1      = 0        ! Used to output file instead of prints
          k5      = nint(abs(ct(1,l))) + 1  !  Avoid *.re* file if > 0
          rfl     = .false.  ! Force recompute of current reactions
          rffl    = .false.  ! Prevent transform to Cartesian
          ct(1,l) = -1.d0    ! Include applied loads on reactions
        else

c         Set range of print

          k1 = nint(abs(ct(1,l)))
          k1 = max(1,min(numnp,k1))
          k2 = nint(ct(2,l))
          if(k2.eq.0) k2 = k1
          k2 = max(1,min(numnp,k2))
          k3 = nint(ct(3,l))
          if(k3.eq.0) k3 = 1
          if(k2-k1.ne.0) k3 = sign(k3,k2-k1)
          if (pcomp(lct(l),'imag',4)) then
            if(cplxfl) then
              k4 = nneq
            else
              write(iow,3003)
              if(ior.lt.0) write(iow,3003)
            endif
          endif
          pfr = .true.
        endif
      endif

c     Compute new reactions

      if(.not.rfl) then
        call pzero(hr(np(26)),nneq*ipc)
        if(ct(1,l).lt.0.0d0) then
          call ploa1(ttim,dt)
          call pload(mr(np(31)+nneq),mr(id31),hr(np(40)),hr(np(30)),
     &               hr(np(26)),prop*rlnew,fa,fa)
!    &               hr(np(26)),prop*rlnew,tr,fa)
        endif
        pltmfl = .true.
        call formfe(np(40),np(26),np(26),np(26),fa,tr,fa,tr,6,1,numel,1)
        pltmfl = .false.
c       Compute current work: Energy = U x R
        ur = -dot(hr(np(26)),hr(np(40)),nneq)
        if(pfeap_on) then   ! Exchange with other processes
          dumar(1) = ur
          call pfeapsr(dumar,tdatabuf,1,.true.)
        endif
c       Convert reactions to cartesian form if necessary
        if(rffl) then
          call pcartre(hr(np(26)),hr(np(45)),ndf,numnp)
        endif
        if (pcomp(lct(l),'nopr',4)) then
          rfl = .true.
          pfr = tfl
          return
        endif

      endif

c     List outputs

      if(pcomp(lct(l),'list',4)) then
        k1 = nint(abs(ct(1,l)))
        k1 = max(1,min(3,k1))
        call prtrel(hr(np(26)),hr(np(43)),hr(np(40)),
     &              ttim,ndm,ndf,niols(k1),iolist(1,k1),prth)

c     Selected nodal outputs

      elseif(k1.gt.0) then
        call prtrea(hr(np(26)+k4),hr(np(43)),hr(np(40)),
     &              ndm,ndf,k1,k2,k3,prth)
      endif

c     Set flag to prevent recomputation of reactions for same state

      rfl = .true.
      pfr = tfl

c     Output work: Energy = U x R

      if(k1.gt.0) then

        write(iow,2005) ur
        if(ior.lt.0) write(*,2005) ur

c     Save reactions to file: fsav.reac0000 to fsav.reac9999

      else

c       Set reaction file number

        call setext('reac',k1,fext,.false.)
61      fnamr = fsav
        call addext(fnamr,fext,128,8)
        inquire(file=fnamr,exist=exst)
        if(exst) then
          k1 = k1 + 1
          call setext('reac',k1,fext,.false.)
          go to 61
        endif

c       Find maximum value of reactions and limit size of outputs

        ur = 0.0d0
        do i = 0,nneq-1
          ur = max(ur,abs(hr(np(26)+i)))
        end do ! i
        ur = sqrt(tol)*ur

c       Output reactions larger than sqrt(tol)*max to file
c       N.B.  Always save current state to file 'fsav.reacxxxx'
c             Unless ct(1,l) >= 1 or ct(2,l) .ge.2

        k5 = max(k5,min(2,nint(ct(2,l))))
        do kk = k5,2
          if(kk.eq.2) then
            fnamr = 'reaction'
            fext  = 'reac'
          endif
          call opnfil(fext,fnamr,-1,ios,cfr)
          call preaout(mr(np(31)),hr(np(26)), hr(np(41)), ur, ndf,numnp)
          close(ios,status = 'keep')
          write(iow,2007) fnamr
          if(ior.lt.0) write(*,2007) fnamr
        end do ! kk

c       Force recompute after file output

        rfl = .false.

      endif

      return

c     Check mesh for input errors

c     [chec]      - check mesh for errors
c     [chec,init] - check mesh for initialization of history data base
c                   changes

7     if(pcomp(lct(l),'init',4)) then
        k1     = 14
        hflgu  = .true.                ! Permit update on data base
        h3flgu = .true.
      else
        k1 = 2

c       Report restraint codes on dof's

        call ckfixed(mr(np(190)),mr(np(31)+ndf*numnp) )

      endif

c     Check individual elements for errors

      call formfe(np(40),np(26),np(26),np(26),fa,fa,fa,fa,k1,1,numel,1)

      return

c     Check stress errors

c     [erro]
c     [erro,stre]    --- measures stress error for adaptivity
c     [erro,ener]    --- measures energy error for adaptivity

c     Initialize and compute error indicators: isw = 11

8     if(fl(11)) then
        if(pcomp(lct(l),'stre',4)) ierr = 1
        if(pcomp(lct(l),'ener',4)) ierr = 2
        ertyp  = 1
        em1    = 0.0d0
        em2    = 0.0d0
        eener  = 0.0d0
        eenere = 0.0d0
        eerror = 0.0d0
        eproj  = 0.0d0
        efem   = 0.0d0
        arsq   = 0.0d0
        fp(5)  = np(36)
        np(36) = np(60)
        call formfe(np(40),np(26),np(26),np(26),fa,fa,fa,fa,11,
     &              1,numel,1)
        np(36) = fp(5)

c       Output error indicators

        call prterr

c       Compute remeshing parameters: isw = 11

        ertyp  = 2
        etab   = 1.00d0
        em1    = etab*sqrt((abs(efem+eerror))/dble(numel))
        em2    = etab*sqrt((abs(eener+eenere))/dble(numel))
        arsq   = 0.0d0
        eerror = 0.0d0
        eproj  = 0.0d0
        efem   = 0.0d0
        eener  = 0.0d0
        eenere = 0.0d0
        fp(5)  = np(36)
        np(36) = np(60)
        call formfe(np(40),np(26),np(26),np(26),fa,fa,fa,fa,11,
     &              1,numel,1)
        np(36) = fp(5)

c     Write error

      else
        write(  *,3001)
        write(ilg,3001)
        write(iow,3001)
      endif
      return

c --- [damp] form consistent damping matrix (isw=9)

9     compsv =  compfl
      if(compdp) then
        k1 = 0
        call iters(k1,3)
        compdp = .false.
      else
        compfl = .true.
      endif
      if(pcomp(lct(l),'unsy',4)) then
        if(prnt .and. ior.lt.0) write(*,2028) 'UNSYMMETRIC'
        unsfl  = .true.
        neqs   = 1
        write(tname,'(4hDAMP,i1)') npart
        setvar = palloc(npart+16,tname, nnc+nnc-neq, 2)
        nc     = np(npart+16)
        ncu    = nc  + neq
        ncl    = ncu + nnc - neq
        call pzero (hr(nc),nnc+nnc-neq)
      else
        if(prnt .and. ior.lt.0) write(*,2028) 'SYMMETRIC'
        unsfl  = .false.
        neqs   = neq
        nc     = np(npart+16)
        ncu    = nc  + neq
        ncl    = ncu
        call pzero (hr(nc),nnc)
      endif

      jcmplx =  kcmplx
      castif = .false.
      cadamp = .true.
      camass = .false.
      idenf  = .false.
      solvsv =  solver
      solver = .true.

      call formfe(np(40),nc,nc,ncl,tr,fa,unsfl,fa,9,1,numel,1)

      compfl =  compsv
      solver =  solvsv

      return

c --- [augm,,value,gap] perform nested update for augmented lagrangian
c                      'value' is used only for first iteration in
c                       each time step. (default value is 1.0 & normally
c                       is used.) 'gap' is maximum gap permitted for
c                       convergence.

c --- [augm,pena,value] reset augmented penalty parameter only

10    if(pcomp(lct(l),'pena',4)) then
        if(ct(1,l).gt.0.0d0) augf = ct(1,l)
        if(prnt) write(iow,2006) augf
        if(ior.lt.0.and.prnt) write(*,2006) augf
      elseif(pcomp(lct(l),'gap' ,3) .or. pcomp(lct(l),'    ',4)) then
        lvaug  = lv
        if(pcomp(lct(l),'gap' ,3)) then
          augg   = 0.0d0
          auggfl = .false.
        else
          auggfl = .true.
        endif
        if(rnmax.eq.0.0d0) then

c         New time step

          if(ct(1,l).gt.0.0d0) then
            augf = ct(1,l)
          else
            augf = 1.0d0
          endif
          if(prnt) write(iow,2006) augf
          if(ior.lt.0.and.prnt) write(*,2006) augf
        endif

c       Augment element values

        hflgu  = .true.
        h3flgu = .true.
        call formfe(np(40),np(26),np(26),np(26),fa,fa,fa,fa,10,
     &              1,numel,1)

c       Check gap convergence

        if(pcomp(lct(l),'gap',3)) then
          if(augg.lt.ct(2,l)) then
            auggfl = .true.
          endif
        endif

c       Continue with current time step

c       aold  = rnmax
c       aengy = rnmax
        naugm = naugm + 1
        iaugm = 0

      endif
      return

c --- [geom]        - Geometric stiffness formulation for eigenvalues
c     [geom,on/off] - Control to add element geometric stiffness terms

11    if    (pcomp(lct(l),'on', 2)) then
        gflag = .true.
        write(iow,2022)
        if(ior.lt.0) then
          write(*,2022)
        endif
      elseif(pcomp(lct(l),'off',3)) then
        gflag = .false.
        write(iow,2023)
        if(ior.lt.0) then
          write(*,2023)
        endif
      else
        imtyp = 2
        cmtyp = 'cons'
        go to 51
      endif
      return

c --- [dire]ct solution option
c     [direct]         - Profile solver (in-core) Symmetric/unsymmetric
c     [direct,blocked] - Profile solver (to disk) Symmetric/unsymmetric
c     [direct,sparse]  - Sparse (in-core) Symmetric

c     Blocked to disk

12    if(pcomp(lct(l),'bloc',4)) then
        ittyp  = -1
        k1     =  nint(ct(1,l))
        if(k1.gt.0) then
          if(ior.lt.0) write(*,2014) k1
          write(iow,2014) k1
        else
          if(ior.lt.0) write(*,2014)
          write(iow,2014)
        endif
        call iters(k1,1)

c     Direct sparse symmetric

      elseif(pcomp(lct(l),'spar',4)) then
        ittyp  = -2
        if(nint(ct(1,l)).eq.0) then
          domd = .true.
          if(ior.lt.0) write(*,2020)
          write(iow,2020)
        else
          domd = .false.
          if(ior.lt.0) write(*,2021)
          write(iow,2021)
        endif
        k1     =  0
        call iters(k1,1)

c     Direct profile

      else

c       Clean up memory if sparse or blocked used before

        if(np(67).ne.0) then
          setvar = palloc( 67,'SPTAN', 0, 2)
          setvar = palloc(307,'IPTAN', 0, 1)
        endif
        if(np(94).ne.0) then
          setvar = palloc( 93,'OINC ', 0, 1)
          setvar = palloc( 94,'OINO ', 0, 1)
        endif

c       Set flags and solution type

        k1 =  nint(ct(1,l))
        if(k1.eq.1) then
          soltyp = 1
        else
          soltyp = 2
        endif
        compfl = .false.
        ittyp  = -3
        if(ior.lt.0) write(*,2015)  soltyp
        write(iow,2015) soltyp
        k1     =  nint(ct(1,l))
      endif

c     Save solution type for this partition

      nittyp(npart) = ittyp
      return

c --- [iter]ation solution option

13    if(pcomp(lct(l),'ppcg',4)) then
        ittyp = max(3,2*ndf,int(ct(1,l)))
        if(ior.lt.0) write(*,2013)
        write(iow,2013)
      elseif(pcomp(lct(l),'bpcg',4)) then
        ittyp = 2
        if(ior.lt.0) write(*,2012)
        write(iow,2012)
      elseif(pcomp(lct(l),'tol',3)) then
        itol          = ct(1,l)
        nitolp(npart) = itol
        if(ct(2,l).ne.0.0d0) then
          atol          = ct(2,l)
          natolp(npart) = atol
        endif
        if(ct(3,l).ne.0.0d0) then
          dtol          = ct(3,l)
        endif
        if(ior.lt.0) write(*,2029) itol,atol,dtol
        write(iow,2029) itol,atol,dtol
        return
      else
        ittyp = 1
        if(ior.lt.0) write(*,2011)
        write(iow,2011)
      endif
      nittyp(npart) = ittyp
      if(nint(ct(3,l)).gt.0) then
        icgits = nint(ct(3,l))
      else
        icgits = neq
      endif
      k1 = 0
      call iters(k1,1)
      return

c     [expo]rt tangent and residual for coupling to external program

14    setvar = palloc(111,'TEMP1',numel    ,  1)

c     Establish list of elements with slaved nodes

      call pelbld(mr(np(33)),mr(id31),mr(np(111)),
     &            nen,nen1,ndf,numel,elist)

c     Assign variables to assemble the export tangent

      setvar = palloc(112,'TEMP2',nen       ,1)
      setvar = palloc(113,'TEMP3',nen       ,1)
      setvar = palloc(114,'TEMP4',neqg      ,2)
      setvar = palloc(115,'TEMP5',neqg*neqg ,2)
      setvar = palloc(116,'TEMP6',neqg*neq*2,2)

c     Compute R (TEMP4); H (TEMP5); G (TEMP6)

      call formgh(mr(np(111)),elist)

      call uexport(hr(np(114)),hr(np(115)),neqg)

      setvar = palloc(116,'TEMP6',0,2)
      setvar = palloc(115,'TEMP5',0,2)
      setvar = palloc(114,'TEMP4',0,2)
      setvar = palloc(113,'TEMP3',0,1)
      setvar = palloc(112,'TEMP2',0,1)
      setvar = palloc(111,'TEMP1',0,1)

      return

c     [impo]rt slaved node solution vector

15    call uimport(lct(l))
      return

c     [ntan]gent,mate,#mate - Numerical compute elmt tangents for mate#
c     [ntan]gent,elem,#elmt - Numerical compute elmt tangent for elmt
c     [ntan]gent,off,       - Turn off all numerical computes
c     [ntan]gent,cont,#pair - Numerical contact tangent for pair#

c     Compute tangent numerically for element 'k1'

16    if(pcomp(lct(l),'elem',4)) then
        k1 = max(1,min(numel,nint(ct(1,l))))
        write(iow,2008) k1
        if(ior.lt.0) then
          write(*,2008) k1
        endif
        call ptdiff(k1,ct(2,l),.true.)

c     Compute tangent numerically for material 'k1'

      elseif(pcomp(lct(l),'mate',4)) then
        k1 = max(1,min(nummat,nint(ct(1,l))))
        write(iow,2009) k1
        if(ior.lt.0) then
          write(*,2009) k1
        endif
        mr(np(32) + nie*k1 - 8) = 1  ! Set ie(nie-7,k1) material

c     Compute tangent numerically for contacts

      elseif(pcomp(lct(l),'cont',4)) then

        call contact (223)

c     Disable numerical tangent computation

      elseif(pcomp(lct(l),'off ',4)) then
        write(iow,2010)
        if(ior.lt.0) then
          write(*,2010)
        endif
        do i = 1,nummat
          mr(np(32) + nie*i - 8) = 0
        end do ! i
      endif
      return

c     [base] - Compute base modes for multiple support excitation

17    if(noi.eq.0 .or. noi.eq.6) then
        if(np(125).eq.0) then
          setvar = palloc(125,'NUBAS',ndf*numnp, 1)
        endif
        propsv = prop
        prop   = 0.0d0

c       Determine number of patterns

        kk = 0
        do i = 0,ndf*numnp-1
          kk = max(kk,mr(np(125)+i))
        end do ! i
        k1 = kk

c       Parallel

        if(pfeap_on) then

c         Communicate with other processors set vneq size

          call pfeapmi(kk)
          vneq = ndf*numnp
        else
          vneq = neq
        endif

c       Allocate storage for PHI_b, M_b, and temporary U

        setvar = palloc(127,'PHIBS',2*vneq*kk ,2)
        setvar = palloc(111,'TEMP1',ndf*numnp,2)

c       Loop over patterns

        do i = 1,kk

c         Set RHS non-zero pattern

          call baseld(ndf*numnp,i, mr(np(125)), hr(np(27)+ndf*numnp),
     &                hr(np(30)))

          fp(5) = np(127) + 2*vneq*(i-1)
          fp(6) = fp(5)      +   vneq

c         Compute PHI_b

          call pzero(hr(fp(5)),vneq)
          call formfe(np(111),fp(5),na,nal,fa,tr,fa,fa,6,1,numel,1)

          fp(1) = na
          fp(2) = nau
          fp(3) = nal
          fp(4) = np(20+npart)
          call psolve(ittyp,hr(fp(5)),fp,.false.,.true.,.true.,.false.)

c         Normalize base values by ctan(1)

          if(ctan(1).gt.0.0d0) then
            do k3 = 0,vneq-1
              hr(fp(5)+k3) = hr(fp(5)+k3)/ctan(1)
            end do ! k3
          endif

c         Set up M_b: N.B. Call with isw = 19 will use isw = 5 in
c                          elements to compute mass matrix.

          call pzero(hr(fp(6)),vneq)
          call formfe(np(111),fp(6),na,nal,fa,tr,fa,fa,19,1,numel,1)

        end do ! i

c       Destroy temporary storage, reset prop

        setvar = palloc(111,'TEMP1',0,2)

        prop   = propsv

c       Allocate storage for M_b, K_b, W_b and time derivatives.

        setvar = palloc(126,'MASBS',kk*mf,2)
        setvar = palloc(189,'WBASE',kk*3 ,2)
        setvar = palloc(111,'TEMP1',vneq ,2)

c       Loop over patterns

        do i = 1,kk

          fp(5) = np(127) + 2*vneq*(i-1)
          fp(6) = fp(5)      +   vneq

c         Form multiply of mass and static mode

          call basmat(hr(fp(5)),hr(fp(6)),hr(np(npart+8)),
     &                hr(np(npart+12)),hr(np(111)), vneq)

c         Finish projections

          call baspro(hr(np(111)),hr(np(77)),hr(np(126)),
     &                i,mf,neq,vneq,kk)
        end do ! i
        setvar = palloc(111,'TEMP1',0,2)

c       Interchange 'pmass' with all processors

        if(pfeap_on) then
          setvar = palloc(111,'TEMP1',mf*kk*2,2)
          call pfeapsr(hr(np(126)),hr(np(111)),mf*kk,.true.)
          setvar = palloc(111,'TEMP1',0,2)
        endif

c     Not a static mode

      else
        if(ior.lt.0) then
          write(  *,3010)
        else
          write(ilg,3010)
          write(iow,3010)
          call plstop()
        endif
      endif
      return

c     [jint],<grad>,node - J-integral calculation

18    hflgu  = .false.
      h3flgu = .false.
      do i = 1,3
        j_int(i) = 0.0d0
      end do ! i
      if(pcomp(lct(l),'grad',4)) then
        jshft = 1.d0
        write(iow,2025) 'Displacement)'
        if(ior.lt.0) then
          write(*,2025) 'Displacement'
        endif
      else
        jshft = 0.0d0
        write(iow,2025) 'Deformation'
        if(ior.lt.0) then
          write(*,2025) 'Deformation'
        endif
      endif
      call pzero(hr(np(26)),nneq)
      call formfe(np(40),np(26),na,nal,fa,tr,fa,tr,16,1,numel,1)

      write(iow,2024) (j_int(i),i=1,ndm)
      if(ior.lt.0 .and. prnt) then
        write(*,2024) (j_int(i),i=1,ndm)
      endif
      rfl = .true.
      go to 6

c     [zzhu],<off>,<ma> - Zienkiewicz-Zhu projections
c                       - ma .eq. 0: project over all materials
c                            .ne. 0: project over material 'ma'
c                       - 'off' - return to lumped projections

19    call zzprod(lct(l),ct(1,l))
      istv = npstr - 1
      return

c     Solve equations
c     [solv]
c     [solv,line,value] - use line search for energy ratios > value

20    if(neq.gt.0) then
        if(j.eq.20) then
          if(fl(4)) then
            if(ior.lt.0) then
              write(  *,3012) 'SOLV'
              return
            else
              write(ilg,3012) 'SOLV'
              write(iow,3012) 'SOLV'
              call plstop()
            endif
          endif

c         Test for RHS (form)

          if(.not.fl(8)) return

c         Set up solver call for a resolution

          fl(7)    = .false.
          fl(8)    = .false.
          if(abs(aengy).lt.aold) aold = abs(aengy)
          if(nmeth.eq.1) then
            ee = 0.0d0
          elseif(nmeth.eq.2) then
            ee = aengy
          endif

          fp(1)  = na
          fp(2)  = nau
          fp(3)  = nal
          fp(4)  = np(20+npart)
          call psolve(ittyp,hr(np(26)),fp,.false.,.true.,.true.,prnt)
        endif

c       Update iteration counter

        niter  = niter + 1
        iaugm  = iaugm + 1
        lvsol  = lv

c       Set parameter for augmented convergence test

        if(lvaug.gt.0) then
          if(naugm.eq.1) then
            if(rnmax.eq.0.0d0) then
              rnaugm = abs(aengy)
            else
              rnaugm = rnmax
            endif
          endif
          augmfl = .false.
          if(abs(aengy).lt.tol*rnaugm .and.
     &       nint(ct(1,lve(lv))).eq.1) then
            augmfl = .true.
          endif
        endif

c       Check convergence

        if(ct(1,l).le.0.0d0) ct(1,l) = 0.8d0
        if (rnmax.eq.0.0d0) then
          ee     = 0.0d0
          autcnv = .false.
          rnmax  = abs(aengy)
          reln   = 1.d0
          aold   = rnmax/ct(1,l)/0.9999d0
        else
          reln   = (aengy-ee)/rnmax
        endif
        if(pfr) then
          write(iow,2004) rnmax,aengy-ee,reln,tol
          if(ior.lt.0) then
            write(*,2004) rnmax,aengy-ee,reln,tol
          endif
        endif
        if(abs(aengy-ee).le.tol*rnmax .or. linear .or.
     &      abs(rnorm)/dble(toteq) .lt.
     &      abs(rnorm1/rnormn)*sqrt(tol)*1.d-3    .or.
     &      abs(aengy-ee).lt.enzer               ) then
          if(lv.gt.1) then
            ct(1,lve(lv)) = ct(1,lvs(lv))
            l             = lve(lv) - 1
            floop(1)      = .false.
          endif
          autcnv = .true.

c       Scale back step and update the solution -- no line search

        elseif(testfl .and. niter.le.1) then

          do i = 0,nneq-1
            hr(np(26)+i) = hr(np(26)+i)*testva
          end do ! i

c       Perform line search or arclength updates

        elseif(.not.linear) then

c         Line search

          if(pcomp(lct(l),'line',4).and.abs(aengy).gt.ct(1,l)*aold) then
            setvar = palloc(111,'TEMP1',nneq        , 2)
            setvar = palloc(112,'TEMP2',nneq*(nrt+4), 2)
            step = 1.0

            call serchl(aold,mr(id31),np(111),np(40),
     &                  hr(np(26)),ct(1,l),hr(np(112)),neq,step)

            setvar = palloc(112,'TEMP2', 0, 2)
            setvar = palloc(111,'TEMP1', 0, 2)
          endif
          if(arcf) then
            call arclen(hr(np(40)),hr(np(26)),hr(np(84)),hr(np(85)),
     &                  hr(np(27)),mr(id31),ndf,ttim)
          endif
        endif
      else
        call ploa1(ttim,dt )
        call pload(mr(np(31)+nneq),mr(id31),hr(np(40)),hr(np(30)),
     &             hr(np(26)),prop*rlnew,.false.,.false.)
      endif

c     Check solution method

      if(nmeth.eq.2) then
        do i = 0,nneq-1
          j = mr(id31+i) - 1
          if(j.gt.0) hr(np(26)+j) = hr(np(26)+j) - hr(np(40)+nneq+i)
        end do ! i
      end if
      call update(mr(id31),hr(np(30)),hr(np(40)),hr(np(42)),
     &            hr(np(26)),fl(9),2)
      fl(8) = .false.

c     Output convergence information

      if(convfl) then
        write(ilg,2030) nstep,niter, rel0,rnorm,rnmax,aengy
      endif
      return

c     [dsol]            - Diagonal solution routine (uses lumped mass)

21    if(np(npart+12).ne.0) then

c       Test for RHS (form)

        if(.not.fl(8)) return

c       Set pointers

        fp(1)  =  np(npart+12)
        call psolve(0,hr(np(26)),fp,.false.,.true.,.false., prnt)
        go to 20
      else
        if(ior.lt.0) then
          write(*,3011)
        endif
        write(iow,3011)
      endif
      return

c     [scal,<on,off>]            - Diagonal scaling of tangent matrix

22    if(np(npart+20).ne.0) then

        if(pcomp(lct(l),'off',3)) then
          if(np(npart+234).ne.0 .and.scale(npart)) then
            write(fnamr,'(a4,i1)') 'DSCA',npart
            setvar = palloc(npart+234,fnamr,  0,2)
            fl(4)  = .true.
          endif
          scale(npart) = .false.
          if(ior.lt.0) write(*,2026)
        else
          write(fnamr,'(a4,i1)') 'DSCA',npart
          setvar = palloc(npart+234,fnamr,neq,2)
          scale(npart) = .true.
          fl(4)        = .true.
          if(ior.lt.0) write(*,2027)
        endif
      endif
      return

c     [hill]-mandel computations
c     [hill tang]  - Compute tangent and stress
c     [hill stre]  - Compute stress only
c     [hill read]  - Read record from file
c     [hill clos]  - Close file

24    if(pcomp(lct(l),'read',4)) then

c       Check if inputs from file

        if(filflg) then
          inquire(file=hillfile, opened=exst)
          if(.not.exst) then
            open(unit=97, file=hillfile, access='sequential')
          endif
          read(97,*,end=8000) (td(i),i=1,10)
          write(iow,2031) td(1)
          ttim = td(1)
          if(finflg) then  ! Load from deformation gradient
            kk   = 1
            do j = 1,3
              do i = 1,3
                kk         = kk + 1
                gradu(i,j) = td(kk)
              end do ! i
            end do ! j
          else             ! Load from strains
            do i = 1,3
              gradu(i,i) = td(i+1)
            end do ! i
            gradu(1,2) = td(5)*0.5d0
            gradu(2,1) = gradu(1,2)
            gradu(2,3) = td(6)*0.5d0
            gradu(3,2) = gradu(2,3)
            gradu(3,1) = td(7)*0.5d0
            gradu(1,3) = gradu(3,2)
          endif
          call mprint(gradu,3,3,3,'GRAD U')
        endif
c     Close file if open
      elseif(pcomp(lct(l),'clos',4)) then
        inquire(file=hillfile, opened=exst)
        if(exst) then
          close(unit=97, status = 'keep')
          filflg = .false.
        endif
c     Compute Hill-Mandel projection of stress and moduli
      else
        aengysv = aengy   ! Save current energy of solution
        call phillmandel(lct(l),ct)
        aengy   = aengysv
      endif
      return

c     End of file read

8000  write(*,*) ' *WARNING* No more records on file'
      return

c     Formats

2000  format('   Shift to tangent matrix = ',1p,e12.5)
2001  format('   Residual norm = ',1p,2e17.7:,6x,'t=',0p,2f9.2)
2002  format(59x,'t=',0p,2f9.2)
2004  format('   Energy convergence test'/
     &       '    Initial   =',1p,e25.15,' Current   =',1p,e25.15/
     &       '    Relative  =',1p,e25.15,' Tolerance =',1p,e25.15)
2005  format('   Energy: Displacements * Reactions = ',1p,e25.15/1x)
2006  format('   Current Augmented Lagrangian Factor =',1p,e13.5)
2007  format('--> Reactions saved to file : ',a/)
2008  format('   Numerical tangent: Element =',i9/)
2009  format('   Numerical tangent: Material =',i9/)
2010  format('   Numerical tangent OFF for all materials.'/)
2011  format('   Iterative Solution: Diagonal Preconditioner'/)
2012  format('   Iterative Solution: Block Diagonal Preconditioner'/)
2013  format('   Iterative Solution: Profile Preconditioner'/)
2014  format('   Direct Solution: Blocked Profile':': Block =',i8/)
2015  format('   Direct Solution: Profile:',i2,' Column reduction'/)
2016  format('   Mass Type DIAGONAL (LUMP) Specified'/)
2017  format('   Mass Type ',a,' Specified'/)
2018  format('   Geometric Stiffness Specified'/)
2019  format('--> SOLVE AT',f9.2,' Mflops. Time=',f12.2)
2020  format('   Minimum degree ordering in solution')
2021  format('   No reordering during solution')
2022  format('   Geometric stiffness ON')
2023  format('   Geometric stiffness OFF')
2024  format('   J_integral(1) =',1p,1e13.6 /
     &       '   J_integral(2) =',1p,1e13.6:/
     &       '   J_integral(3) =',1p,1e13.6 /)
2025  format('   J_integral-Material Force Fracture Indicators'/
     &       '     (Using ',a,' Gradient'/)
2026  format( 10x,'No scaling of tangent matrix')
2027  format( 10x,'Scaling of tangent matrix for unit diagonals')
2028  format('   Damping Type ',a,' Specified'/)
2029  format('   Iteration: Residual tolerance = ',1p,1e13.6/
     &       '              Absolute tolerance = ',1p,1e13.6/
     &       '              Diverge  tolerance = ',1p,1e13.6)
2030  format(2i6,1p,4e16.8)
2031  format(5x,'Deformation Gradient at time =',1p,1e12.3)

3000  format(' *ERROR* PMACR1: Unable to compute starting acceleration,'
     &      ,' Static problem or density zero')

3001  format(' *ERROR* PMACR1: Compute STREss NODE before ERROr')

3002  format(' *WARNING* Unfactored tangent produced do not try'
     &      ,' normal solution.')

3003  format(' *WARNING* Complex mode not active: No imaginary outputs')

3004  format(' *WARNING* No mass defined: Default to diagonal mass.')

3006  format(' *WARNING* No active equations: Solution updated')

3007  format(' *ERROR* Use SOLVe before another FORM.'/)

3008  format(' *ERROR* PMACR1: No non-zero terms in mass matrix:',
     &       ' Check density value for materials')

3009  format(' *WARNING* No flexible equations, mass matrix not',
     &       ' required')

3010  format(' *ERROR* PMACR1: Cannot compute static modes.'/
     &       '         Modify algorithm by placing BASE command after'/
     &       '         TANGent command and before TRANsient,CONServing'/
     &       '         command.  If INTERactive solution mode can use'/
     &       '         TRANSient,OFF followed by TANGent and BASE.'/
     &       '         Then turn on the transient solution again using'/
     &       '         TRANSient,CONServing.'/)

3011  format(' *ERROR* DSOL: No lumped mass matrix available'/)

3012  format(' *ERROR* ',a,': No stiffness matrix, use TANG or UTAN'/)

      end
