c$Id:$
      subroutine pmacr (initf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove interupt tests (intf and lintr)           17/11/2006
c       2. Add 'forc' to pmacr3 options                     17/12/2006
c       3. Add 'glob'al and 'loca'l print options           13/03/2007
c       4. Dimension 'fint' to 128                          17/04/2007
c       5. Introduce 'macd' and 'maci' for special use      22/08/2007
c          Added to common 'pconstant.h'
c       6. Add restart from command line request.           12/10/2007
c          Needs 'prmptd.h' for 'reflg'.
c       7. Change 'Macro' to 'Command'                      07/01/2008
c       8. Set 'nocomp' to 1 for component outputs          24/07/2008
c       9. Set default for 'dyncon' to false                25/08/2008
c          (for dynamic contact solutions).
c      10. Remove 'pfr' from argument of pmacr7             13/12/2008
c      11. Initialize intlfl, idisfl, ivelfl, & tranfl      01/02/2009
c          to .false.
c      12. Add 'tie off' command to pmacr3.f                24/02/2009
c      13. Initialize nhplts                                25/06/2009
c      14. Add 'convfl' and outputs for monitory converges  16/02/2010
c      15. Add 'real' and 'imag' for print sets             09/08/2010
c      16. Add 'tolr' and 'resnm' initial set value         10/08/2010
c      17. Add 'hill'-mandel to pmacr1                      19/04/2011
c      18. Add umacr11 to umacr20, umacr0->umacr10          25/01/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Command language instruction subprogram:
c               Controls problem solution and output algorithms by
c               order of specifying commands in array wd.

c      Inputs:
c         initf     - Flag, Initialize solution data if true

c      Outputs:
c         none      - Routine may be called several times
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'arclel.h'
      include   'arclei.h'
      include   'arcler.h'
      include   'augdat.h'
      include   'auto1.h'
      include   'auto2.h'
      include   'cdata.h'
      include   'cdat1.h'
      include   'chdata.h'
      include   'comfil.h'
      include   'compac.h'
      include   'compas.h'
      include   'corfil.h'
      include   'counts.h'
      include   'ddata.h'
      include   'debugs.h'
      include   'dyndat.h'
      include   'elacts.h'
      include   'elauto.h'
      include   'eldata.h'
      include   'eltran.h'
      include   'endata.h'
      include   'fdata.h'
      include   'gltran.h'
      include   'hdata.h'
      include   'hlpdat.h'
      include   'iofile.h'
      include   'ldata.h'
      include   'modcon.h'
      include   'ndata.h'
      include   'part7.h'
      include   'pbody.h'
      include   'pconstant.h'
      include   'pdata2.h'
      include   'pdata3.h'
      include   'plflag.h'
      include   'plist.h'
      include   'pmatlab.h'
      include   'pmod2d.h'
      include   'prflag.h'
      include   'print.h'
      include   'prld1.h'
      include   'prlod.h'
      include   'prmptd.h'
      include   'ptdat1.h'
      include   'ptdat2.h'
      include   'ptdat3.h'
      include   'ptdat4.h'
      include   'ptdat5.h'
      include   'ptdat6.h'
      include   'ptdat7.h'
      include   'ptdat8.h'
      include   'ptdat9.h'
      include   'ptdata.h'
      include   'ptdatb.h'
      include   'ptest.h'
      include   'pview.h'
      include   'rdata.h'
      include   'rdat0.h'
      include   'region.h'
      include   'sdata.h'
      include   'setups.h'
      include   'tdata.h'
      include   'tdatb.h'
      include   'tdato.h'
      include   'umac1.h'
      include   'xtout.h'
      include   'pointer.h'
      include   'comblk.h'
      include   'comnds.h'

      integer    ncomd
      parameter (ncomd = 140)

      character  wd(ncomd)*4,lct(ncmds)*15
      integer    ed(ncomd),  jct(ncmds)

      logical    initf, errck,vinput,pcomp, tjump
      character  fint*128
      integer    nlp,nif,njp,ntm, i, j, k, ll, llast
      integer    nwd1,nwd2,nwd3,nwd4,nwd5,nwd6,nwd7,nwd8,nwdp
      integer    nw1,nw2,nw3,nw4,nw5,nw6,nw7,nw8,nwp
      real*4     tary(2), etime , tt
      real*8     td(4)

      save

c     List of entries in array wd for solution commands.
c     N.B.  The continue label 'n' indicates which pmacr'n'
c           subprogram contains command language statements.

      data wd/'stre','utan','tang','form','mass','reac','chec','erro',
     1        'damp','augm','geom','dire','iter','expo','impo','ntan',
     1        'base','jint','zzhu','solv','dsol','scal','resi','hill',

     2        'tol ','dt  ','loop','next','prop','data','time','prin',
     2        'nopr','beta','init','iden','newf','back','debu','line',
     2        'nonl','auto','meth','if  ','else','endi','tran','step',
     2        'jump','echo','conv','omeg',

     3        'disp','test','mesh','plot','subs','writ','read','cont',
     3        'rest','velo','acce','bfgs','arcl','save','paus','eige',
     3        'expl','memo','acti','deac','zero','epri','moda','opti',
     3        'eigv','rayl','cxso','broy','comp','rect','cyli','sphe',
     3        'forc','tie ','real','imag',

     4        'mac1','mac2','mac3','mac4','mac5','mac6','mac7','mac8',
     4        'mac9','ma10','ma11','ma12','ma13','ma14','ma15','ma16',
     4        'ma17','ma18','ma19','ma20',

     5        'outm','renu','show','scre','comm','smoo','set ','assi',
     5        'outp',

     6        'list','tplo','para','func','dync','part','mate','capt',
     6        'mono','flag','get ',

     7        'grap','outd','pets','psub','u_in','glis','gplo','glob',
     7        'loca',

     8        'elev','inse',

     p        'manu' /

      data ed/    0,     0,     0,     0,     0,     0,     0,     1,
     1            1,     1,     1,     3,     3,     1,     1,     1,
     1            1,     3,     1,     0,     3,     1,     3,     1,

     2            0,     0,     0,     0,     0,     1,     0,     0,
     2            0,     5,     0,     1,     1,     1,     0,     2,
     2            2,     1,     2,     1,     1,     1,     0,     1,
     2            2,     1,     1,     1,

     3            0,     0,     0,     0,     0,     1,     1,     1,
     3            1,     0,     0,     1,     1,     1,     2,     0,
     3            1,     1,     1,     1,     1,     0,     1,     1,
     3            1,     1,     1,     1,     2,     1,     1,     1,
     3            0,     1,     1,     1,

     4            5,     5,     5,     5,     5,     5,     5,     5,
     4            5,     5,     5,     5,     5,     5,     5,     5,
     4            5,     5,     5,     5,

     5            1,     1,     0,     3,     1,     3,     1,     1,
     5            1,

     6            0,     0,     0,     1,     2,     1 ,    2,     0,
     6            1,     1,     0,

     7            3,     3,     3,     3,     3,     3,     3,     3,
     7            3,

     8            3,     3,

     p            4 /

      data    nwd1,nwd2,nwd3,nwd4,nwd5,nwd6,nwd7,nwd8,nwdp
     &       /  24,  28,  36,  20,   9,  11,   9,   2,   1 /

      if(initf) then

c       Write Header to Solution file

        write(ilg,4000)

c       Set counter values

        nstep  = 0
        nsplt  = 0
        niter  = 0
        nform  = 0
        iform  = 0
        naugm  = 0
        titer  = 0
        taugm  = 0
        tform  = 0
        iaugm  = 0
        intvc  = 0
        iautl  = 0
        iview  = 0
        urest  = 0

        ntang  = 0
        nutan  = 0
        ncmas  = 0
        nlmas  = 0
        numas  = 0
        ncdam  = 0
        nudam  = 0
        nrfrm  = 0

c       Set initial values of parameters

        enzer  = 0.0d0
        aengy  = 0.0d0
        aold   = 0.0d0
        augf   = 1.0d0
        rnmax  = 0.0d0
        shift  = 0.0d0
        tol    = 1.d-16
        tolr   = 1.d-08
        itol   = 1.d-08
        atol   = 1.d-16
        dtol   = 1.d+16
        dt     = 0.0d0
        dtn    = 0.0d0
        dtold  = 0.0d0
        prop   = 1.0d0
        propo  = 1.0d0
        resnm  = 1.0d0
        ttim   = 0.0d0
        dtmax  = 0.0d0
        dtmin  = 0.0d0
        nocomp = 1
        pmaxit = 0

        call pzero(prldv,50)

c       Compute side lengths of mesh for tolerancing

        call pdxmsh(mr(np(190)),hr(np(43)))

c       Arc-length method

        rlnew  =  0.0d0
        timold = -1.0d0
        kflag  =  0
        maplt  =  0
        msplt  =  0
        noi    =  0
        nreg   = -1

c       Dynamic parameters

        call dparam(bpr,'init')
        cc1    = 1.d0
        cc2    = 1.d0
        cc3    = 1.d0
        do i = 1,3
          bpr(i)  = 0.0d0
          ctan(i) = 0.0d0
          gtan(i) = 0.0d0
        end do ! i
        ctan(1)  = 1.0d0
        gtan(1)  = 1.0d0
        convfl   = .false.
        dynflg   = .false.
        expflg   = .false.
        flgexp   = .false.
        floop(1) = .false.
        floop(2) = .false.
        ljump    = .false.
        gflag    = .true.
        modfl    = .true.
        intlfl   = .false.
        idisfl   = .false.
        ivelfl   = .false.
        qoutfl   = .true.
        tranfl   = .false.
        ndebug   = 0
        numint   = 5
        alpha    = 1.d0
        rayla0   = 0.d0
        rayla1   = 0.d0

        do i = 1,4
          steps(i) = 1.d0
          fstep(i) = .true.
        end do ! i

c       Save room for history storage (first check if interface)

        if(np(210).eq.0) then
          i = 1
        else
          i = 2
        endif

        theta(1) = 0.0d0
        theta(2) = 0.0d0
        theta(3) = 1.0d0

        zview0(1)= 0.0d0
        zview0(2)= 0.0d0
        zview0(3)= 0.0d0

        rmeas    = 0.0d0
        rvalu(1) = 0.0d0
        rvalu(2) = 0.0d0
        rvalu(3) = 0.0d0

c       Set default solution flags

        arcf   = .false.
        autofl = .false.
        autofli= .false.
        aratfl = .false.
        compfl = .false.
        compdp = .true.
        compms = .true.
        cadamp = .true.
        castif = .true.
        camass = .false.
        causer = .false.
        compre = .false.
        cview  = .false.
        debug  = .false.
        fl( 1) = .false.
        fl( 2) = .false.
        fl( 3) = .true.
        fl( 4) = .true.
        fl( 5) = .true.
        fl( 6) = .true.
        fl( 7) = .true.
        fl( 8) = .false.
        fl( 9) = .false.
        fl(10) = .true.
        fl(11) = .false.
        fl(12) = .true.
        do i = 1,4
          ofl9(i)  = fl(9)
          flp(9,i) = fl(9)
        end do ! i
        fops   = .true.
        dyncon = .false.
        echo   = .false.
        hadd   = .true.
        linear = .false.
        monofl = .false.
        solnfl = .false.
        pfl    = .false.
        pfr    = .true.
        plfl   = .true.
        prnt   = .true.
        refl   = .false.
        rfl    = .false.
        screfl = .true.
        testfl = .false.
        trifl  = .false.

c       Set integer parameters

        ittyp  = -3   ! default is incore profile solver
        icgits = neq
        itract = 0
        itrdea = 0
        kcmplx = 0
        mcmplx = 0
        dcmplx = 0
        ucmplx = 0
        li     = 0
        lvcn   = -1
        maxbl  = 0
        numact =-1
        numdea =-1
        npld   = 0
        naplts = 0
        ncplts = 0
        ndplts = 0
        neplts = 0
        nhplts = 0
        nlplts = 0
        nmplts = 0
        nqplts = 0
        nrplts = 0
        nsplts = 0
        ntplts = 0
        nuplts = 0
        nvplts = 0
        ntstep = 0
        nc     = 1
        nv     = 1
        nw     = 1
        niols(1) = 0
        niols(2) = 0
        niols(3) = 0
        npstr    = max(11,npstr+1)

        call pzero(epl,200)

c       Set initf to prevent reinitializing parameters

        initf  = .false.

c       Set user command names for default values

        nw4 = nwd1 + nwd2 + nwd3
        do j = 1,20
          i = nw4+j
          if(.not.pcomp(umacc(j),wd(i),4)) then
            wd(i) = umacc(j)
            ed(i) = 0
          endif
        end do ! i

      endif

c     Restart from command line execution

      if(reflg) then
        reflg = .false.
        write(iow,'(/a,a)') '  --> Restart from : ',fres
        if(ior.lt.0) then
          write(*,'(/a,a)') '  --> Restart from : ',fres
        endif
        call restrt(fres,ndm,ndf,nneq,1)
      endif

c     Set pointers to pmacr"n" subprograms

      nlp    = nwd1 + 3
      ntm    = nwd1 + 7
      nif    = nwd1 + 20
      njp    = nwd1 + 25
      nw1    = nwd1
      nw2    = nwd2 + nw1
      nw3    = nwd3 + nw2
      nw4    = nwd4 + nw3
      nw5    = nwd5 + nw4
      nw6    = nwd6 + nw5
      nw7    = nwd7 + nw6
      nw8    = nwd8 + nw7
      nwp    = nwdp + nw8

c     Input solution command statements


100   call pmacio (jct,lct,ct,wd,ed,nwp,nlp,nif,ll)
      if(ll.le.0) go to 300

c     Look for jump statment

      ljump = .false.
      tjump = .false.
      njump = 0
      vjump = 0
      lv    = 1
      do l = 2,ll-1
        if(jct(l).eq.nlp  ) lv = lv + 1
        if(jct(l).eq.nlp+1) lv = lv - 1
        if(jct(l).eq.ntm) then
          tjump = .true.
        endif
        if(tjump .and. jct(l).eq.njp) then
          njump = l
          vjump = lv
          cjump = lct(l)
        endif
      end do ! l

c     Execute solution command instruction program

      nh1 = np(50)
      nh2 = np(51)
      nh3 = np(52)
      lv = 0
      l = 1
200   j = jct(l)
      i = l - 1
      tt = etime(tary)
      if(j.ne.nlp .and. j.ne.nlp+1) then
        errck   = vinput(lzz(l),80,td,4)
        ct(1,l) = td(1)
        ct(2,l) = td(2)
        ct(3,l) = td(3)
        macd    = td(4)
        maci    = nint(td(4))
        write(yyy,2003) wd(j),lct(l),(ct(k,l),k=1,3)

c       Strip leading blanks and comments

        call pstrip(xxx,yyy,3)

c       Set yyy value for input use

        call acheck(xxx,yyy,15,75,75)

      endif
      if((l.ne.1.and.l.ne.ll).and.pfr) then
        if(prnt) then
          if(rank.eq.0) then
            write(iow,2001) i,wd(j),lct(l),(ct(k,l),k=1,3),tary
            if(ior.lt.0 .or. (ior.gt.0.and.debug) .or. echo) then
              write(*,2001) i,wd(j),lct(l),(ct(k,l),k=1,3),tary
            endif
          else
            write(iow,2004) rank,i,wd(j),lct(l),(ct(k,l),k=1,3),tary
            if(ior.lt.0 .or. (ior.gt.0.and.debug) .or. echo) then
              write(*,2004) rank,i,wd(j),lct(l),(ct(k,l),k=1,3),tary
            endif
          endif
          call pflush(iow)
        endif
      endif

c     Transfer to correct subprogram to execute command

      llast = ll
      if(j.le.nw1) then
        call pmacr1(lct,ct,j)
      elseif(j.le.nw2) then
        call pmacr2(lct,ct,j-nw1)
      elseif(j.le.nw3) then
        call pmacr3(lct,ct,j-nw2)
      elseif(j.le.nw4) then
        uct = wd(j)
        call pmacr4(ct(1,l),lct(l),j-nw3)
        inquire(file = finp, opened = errck)
      elseif(j.le.nw5) then
        call pmacr5(lct(l),ct(1,l),j-nw4)
      elseif(j.le.nw6) then
        call pmacr6(lct(l),ct(1,l),pfr,j-nw5)
      elseif(j.le.nw7) then
        call pmacr7(lct(l),ct(1,l),j-nw6)
      elseif(j.le.nw8) then
        call pmacr8(lct(l),ct(1,l),j-nw7)
      elseif(j.eq.nwp) then
        hlplev = max(-1,min(3,int(ct(1,l))))
      endif
      l = l + 1
      if(l.le.llast) go to 200
      if (ior.lt.0) go to 100
300   tt = etime(tary)
      write(iow,2000) tary
      if(ior.lt.0) write(*,2000) tary

c     Save restart information

      if(fl(7)) return
      if(ll.eq. -1) then
        fint = fsav
        call restrt(fint,ndm,ndf,nneq,2)
        if(ior.lt.0) then
          write(*,2002) fint
        endif
        write(iow,2002) fint
      endif

c     Formats

2000  format(' *End of Solution Execution*',31x,'t=',2f9.2)
2001  format(' *Command ',i3,' * ',a4,1x,a15,
     &   'v:',3g11.3/59x,'t=',2f9.2)
2002  format(/'           Saved  Restart  File: ',a)
2003  format (a4,',',a4,',',3(1p,1e14.7,','))
2004  format(' Rank[',i3,']*Command ',i3,' * ',a4,1x,a15,
     &   'v:',3g11.3/50x,'t=',2f9.2)

4000  format(/'  SOLUTION SUMMARY'/'  ----------------'//
     & '  Load    Total     Solution      Time     Residual Norm  ',
     & '      Energy Norm     CPU Time'/
     & '  Step  Tang+Forms      Time      Incr.  Initial     Final',
     & '   Initial     Final (Seconds)'/2x,86('-'))

      end
