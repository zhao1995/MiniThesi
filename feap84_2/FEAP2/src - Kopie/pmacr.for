       subroutine pmacr (ul,xl,tl,ld,p,s,ie,d,id,x,ix,f,t,jd,f0,u,
     1                dr,ct,ncmds,ndf,ndm,nen1,nst,prt,dir,plo)
c ----------------------------------------------------------------------
c.... Purpose: macro instruction subprogram
c              Controls problem solution and output algorithms by order
c              of specifying macro commands in array wd.
c
c....  Inputs:
c     nn    ul(nst,6)      - element displacements
c     n0    xl(ndm,nen)    - element coordinates
c     n1    tl(nen)        - element temperatures
c     n2    ld(nst)        - eq. numbers of dofs of actual element
c     n3    p(nst)         - element load vector
c     n4    s(nst,nst)     - element stiffness matrix
c     n5    ie(nie,numat)  - assembly information for material set nie=ndf+2
c     n6    d(ndd,numat)   - material set parameters
c     n7    id(ndf,numnp)  - equation numbers for each active dof
c     n8    x(ndm,numnp)   - nodal coordinates of mesh
c     n9    ix(nen1,numel) - element nodal connections of mesh
c     n10   f(ndf,numnp)   - load vector
c     n11   t(numnp)       - temperature vector
c     n12   jd(*)          - pointer array for row/columns of tangent
c     n13   f0(ndf,numnp)  - nodal initial force values
c     n14   u(ndf,3*numnp) - displacement vector
c     noff  dr(nneq)       - working array
c     nct   ct(3,ncmds)    - parameter array for 200 macros
c           ncmds          - number of commands
c           ndf            - number dof/node
c           ndm            - spatial dimension of mesh
c           nen            - max. number of nodes/element
c           nen1           - dimension for ix array: nen+4
c           nst            - dimension for element array: ndf*nen
c           prt            - print option
c           dir(10,knode,2)- array for director information            -
c           plo(10,nplo)   - array for plotting TPLO-data
c
c     nested loops in FEAP (up to nine loops possible)
c     lv                   - number of loop
c     lvs(lv)              - number of macro at loop lv starts ??
c     lve(lv)              - number of macro at loop lv ends   ??
c
c      Outputs:
c         none             - Routine may be called several times
c
c     COMMENTS
c     exit from PMACR via
c     ll=-1 exit,e                                    set in PMACIO
c     ll=-2 QUIT,q,RINP(irfeap=2), RINP,new(irfeap=1) set in PMACIO
c     ll=-2 REME,new(irfeap=3), REME,old(irfeap=3)    set in PMACR
c
c ----------------------------------------------------------------------
      USE arcext
      USE arcl
      USE augdat
      USE bdata
      USE cdata
      USE comfil
      USE conv
      USE damp1
      USE ddata
      USE debugs
      USE dirdat
      USE dyndat
      USE eig1
      USE endata
      USE ext2
      USE fdata
      USE fodata
      USE hdata
      USE hdatam
      USE idata
      USE iofile
      USE ldata
      USE macprt
      USE ndata
      USE pcent
      USE pdam
      USE pdata2
      USE pdata3
      USE plodf
      USE plodfa
      USE plodfb
      USE plslay
      USE prlod
      USE psize
      USE qload
      USE rdata
      USE soltyp
      USE subdt
      USE tdata
      USE timex
      USE uneig
      USE doalloc
      implicit real*8 (a-h,o-z)
c.....Declare variable types
      logical fa,tr
      logical prt
      logical pcomp

c.....Declare character types
      character fint*229
      character*4 titl

c.....Declare array types
      character*4 wd(113),wds(113),lct(ncmds)
      integer ld(*),   ie(*), id(*), ix(nen1,*),  jd(*), jct(ncmds)
      real*4  tary
      real*8  ct(3,*), ul(*), xl(*), tl(*), p(*), s(*), d(*)
      real*8  x(*),    f(*),  t(*),  u(*),  dr(*),f0(*),dir(*),plo(10,*)

c.... list of entries in array wd for macro commands
c       N.B.  The continue label 'n' indicates which pmacr'n'.f
c       file contains the macro command statements.
c
      data wd/'stre','utan','tang','form','lmas','cmas','reac','chec',
     1        'erro','damp','augm','geom','pbcg','jint','crit','show',
     1        'detk','summ','umas','auto','iimp','pgmr','prco','nsys',
     1        'updh','splo',
     *
     2        'tol ','dt  ','loop','next','prop','data','time','prin',
     2        'nopr','mate','tran','init','iden','newf','back','debu',
     2        'line','nonl','conv','rinp','if  ','else','endi','eas ',
     2        'dibc','beta',
     *
     3        'disp','solv','mesh','plot','subs','writ','read','cont',
     3        'rest','copy','velo','acce','bfgs','arcl','save','four',
     3        'fsol','fsum','paus','tplo','eigk','ueig','ext ','lamb',
     3        'curv','reme','eigi','eig1','edit','post','smoo','tec ',
     3        'parv','pola','lan ','sigq','epsq','dplo','feas',
     *
     4        'mac1','mac2','mac3','mac4','mac5','elem','man ','para',
     *
     5        'zcop','bsys',
     *
     6        'ytab','yang','ygra','ymsh','yevo','ytry',
     *
     7        'end','exit','help','hist','proc','quit'/

c.... macros which are in WD as dummys
c     END
c     EXIT
c     HELP
c     HIST
c     PROC
c     QUIT

c.... macros which are not documented
c     BSYS ok write COOR,ELEM,MATE to file Bname-binary  
c     NSYS Testversion Change solver 0-12, or calculate new profile
c     EAS  switch for EAS 
c     EDIT call user-defined editor
c     ZCOP  
c     ELEM menue for choosing elements 

      data nwd1,nwd2,nwd3,nwd4,nwd5,nwd6/26,26,39,8,2,6/
      data zero,one,tol_en/0.0d0,1.0d0,1.d-16/
      data fa,tr/.false.,.true./

c.... set initial values of parameters (which are saved in commons or data)
      wds    = wd   ! save due to to new names for mac1 etc
      aengy  = zero
      aold   = zero
      augf   = one
      rnmax  = zero
      shift  = zero
      tol    = tol_en
      dt     = zero
cw    omega  = zero
      prop   = one
      ttim   = zero
      propold = 0.d0       ! needed for convergence check
      intn   = 0
c.... arc-length method
      rlnew  = zero
      timold = zero
      kflag  = 0
      mu1    = 1
      mu2    = 1
      nop    = 0
      arcf   = fa
      refl   = fa
      foflg  = fa
      fostr  = fa
c.... extended system
      kflg   = fa
      kdig   = -6
      exeps  = one
      xmu    = zero
      kex    = 0
c.... eigenvalue computation using overrelaxation
      eigflg = fa
      nxeig1 = 0
      nxeig2 = 0
      nte    = 1
c.... complex eigenvalues
      !meigr  = 0
c.... determinant calculation
      det0  = 0.d0
c.... load defelection curve
      npldf  = 0
      mkflg  = 0
      flreac = fa
      nincp  = 1
      nincp0 = 1
c.... load QLOA
      propq = 1.d0
c.... if/else/endi
      li=0
c
      fl(12) = tr
      hadd   = tr
      linear = fa

c.... local history variable storage arrays for handling in pform nh1+nh2+nh3
c.... elemt history values are stored in
      if(nhmax.gt.0) then
        dbgtxt = 'PMACR: room for history: nh1,nhmax*ipr(hflgu)'
        call ralloc(eh1,nhmax,'History NH1',hflgu)
        dbgtxt = 'PMACR: room for history: nh2,nhmax*ipr(hflgu)'
        call ralloc(eh2,nhmax,'History NH2',hflgu)
      end if
      if(nh3max.gt.0) then
        dbgtxt = 'PMACR: room for history: nh3,nh3max*ipr(h3flgu)'
        call ralloc(eh3,nh3max,'History NH3',h3flgu)
      end if
c
c     calculate nodal director array (after h-array)
      if(ldir.eq.1) then
        call perform(0,1,1,18)
        call pdirec(x,dr,dir,numnp,numel,ndm,prt,u)
        call perform(1,1,3,18)
      end if
c
c
      hflgu  = fa
      h3flgu = fa
      flgc   = tr
      flgda  = tr
      fldyn  = tr
      pfl    = fa
      pfr    = tr
      plfl   = tr
      prnt   = tr
      rfl    = fa
      iclear = 0
      ipola  = 0
      iprd   = 0
      md     = 0  ! eigenvectors
      mf2    = 1  ! eigenvectors to store: mf2*mf mf=1=save mf=2 all
      !mvi    = 1  ! for inverse eigenvalue-iteration
      nneq   = ndf*numnp
      npld   = 0
      npstr  = 27         ! 1*dt + 25 stresses + 1*prin_stress
      nc     = 1  ! damping matrix diag
      ncu    = 1  ! damping matrix upper
      ncll   = 1  ! damping matrix lower
      nv     = 1  ! trans. terms velo,acce
      nw     = 1  !

      jtime  = 0

      nlp    = nwd1 + 3  ! Location of 'loop' command
      nif    = nwd1 + 21 ! Location of 'if'   command

      nw1    = nwd1
      nw2    = nwd2 + nw1
      nw3    = nwd3 + nw2
      nw4    = nwd4 + nw3
      nw5    = nwd5 + nw4
      nw6    = nwd6 + nw5
c.... check mesh or set more values via isw=2
      call formfe(u,dr,dr,dr,fa,fa,fa,fa,2,1,numel,1)
c.... input the macro commands
100   call pmacio (jct,lct,ct,wd,nw6,nlp,nif,ll,ncmds,asc1,asc2)

      if(ll.le.0) go to 300
c.... execute macro instruction program
      lv = 0
      l = 1
200   j = jct(l)
      i = l - 1
      call etimef(tary)
      if((l.ne.1.and.l.ne.ll).and.pfr) then
          write(iow,2001) i,wd(j),lct(l),(ct(k,l),k=1,3),tary
cww       if(ior.lt.0.and.prnt) then
cww         write(*,2001) i,wd(j),lct(l),(ct(k,l),k=1,3),tary
cww       end if
      end if
c.... transfer to correct subprogram to execute macro
      if(j.le.nw1)
     1     call pmacr1(id,ix,x,f,f0,jd,u,dr,lct,ct,ndf,nen1,ndm,nneq,j,
     2                 prt,plo,ld)
      if(j.ge.nw1+1.and.j.le.nw2)
     1     call pmacr2(id,ix,f,f0,u,dr,lct,ct,ndf,nen1,nneq,prt,j-nw1,
     2                 plo)
      if(j.ge.nw2+1.and.j.le.nw3)
     1     call pmacr3 (ul,xl,tl,ld,p,s,ie,d,id,x,ix,f0,f,t,jd,u,dr,
     2            lct,ct,ndf,ndm,nen1,nst,nneq,prt,j-nw2,plo,ll) ! ll for [reme,new+old] 
      if(j.ge.nw3+1.and.j.le.nw4)
     1     call pmacr4 (ul,xl,tl,ld,p,s,ie,d,id,x,ix,f0,f,t,jd,u,dr,
     2        lct(l),ct(1,l),ndf,ndm,nen1,nst,nneq,prt,wd,j-nw3,nw3,plo)
      if(j.ge.nw4+1.and.j.le.nw5)
     1     call pmacr5 (ul,xl,tl,ld,p,s,ie,d,id,x,ix,f0,f,t,jd,u,dr,
     2        lct(l),ct(1,l),ndf,ndm,nen1,nst,nneq,prt,wd,j-nw4,nw4)
      if(j.ge.nw5+1.and.j.le.nw6)
     1     call pmacr6 (ul,xl,tl,ld,p,s,ie,d,id,x,ix,f0,f,t,jd,u,dr,
     2        lct(l),ct(1,l),ndf,ndm,nen1,nst,nneq,prt,wd,j-nw5,nw5)

      if(ll.le.0) go to 300 ! only for [reme,new reme,old], other exits via pmacio 

      l = l + 1

      if (l.le.ll)  go to 200 ! end list of macros to do

      if (ior.lt.0) go to 100 ! new macro for inte

c.... continue with interactive mode (after batch mode)
      if(ior.gt.0) then
        do i=1,5  ! read up to 5 cards searching for 'inte'  or 'stop'
          read(ior,1000) titl
          if(pcomp(titl,'stop',4)) go to 300 ! end
          if(pcomp(titl,'inte',4)) then      ! switch to interactive mode
            close(ior)
            ior = -iabs(ior)
            go to 100
          end if
        end do
      end if
c.... end of calculation
300   call etimef(tary)
                   write(iow,2000) tary
      if(ior.lt.0) write(*  ,2000) tary
c.....return for quit
      if(ll.eq.-2) then
        wd = wds  ! reset wd in case of rinp and new names for mac1 etc in Batch
        return
      end if
c.... save restart information
cww for updh      if(fl(7)) return ! no solution
      fint = fsav
      call restrt(fint,u,ix,ndm,ndf,nen1,nneq,2,asc1,asc2)

      return
c.... formats
1000  format(a4)
2000  format(' *End of macro execution*',40x,'t=',f9.2)
2001  format(' *Macro ',i3,' *',2(a4,1x),
     1   'v1=',g10.3,' v2=',g10.3,' v3=',g10.3,' t=',f9.2)
      end
c
      subroutine pmacr1(id,ix,x,f,f0,jd,u,dr,lct,ct,ndf,nen1,ndm,nneq,j,
     +                  prt,plo,ld)
c ----------------------------------------------------------------------
c.... macro instruction subprogram 1                                   |
c     'stre','utan','tang','form','lmas','cmas','reac','chec',
c     'erro','damp','augm','geom','pbcg','jint','crit','show',
c     'detk','summ','umas','auto','iimp','pgmr','prco','nsys',
c     'updh','splo',
c ----------------------------------------------------------------------
      USE arcext
      USE arcl
      USE augdat
      USE cdata
      USE comfil
      USE conv
      USE damp1
      USE ddata
      USE debugs
      USE dirdat
      USE dspos
      USE dtauto
      USE dyndat
      USE eig1
      USE endata
      USE epsd1      ! for it1,it2,resife2 in FE2
      USE errin1
      USE errin2
      USE errin3
      USE evdata
      USE ext1
      USE ext2
      USE fdata
      USE hdata
      USE hdatam
      USE iimpd
      USE implstep
      USE iofile
      USE iscsr
      USE isgmr
      USE ispcg
      USE isprec
      USE jinteg
      USE ldata
      USE macprt
      USE ndata
      USE ndatx
      USE pcrit
      USE pdam
      USE pdata2
      USE pdata3
      USE plodf
      USE plodfa
      USE plodfb
      USE plong
      USE plslay
      USE pnodn
      USE prlod
      USE psize
      USE qload
      USE rdata
      USE rsum
      USE soltyp
      USE stepc
      USE strnam
      USE sumdt
      USE tdata
      USE doalloc
      implicit real*8 (a-h,o-z)
      logical fa,tr,cfr,pcomp,f8o,prt,cfrm
      logical dflg,flgd
      logical flgdu
      character*4 lct(*)
      character*229 fint1
      character yyy*80
      integer id(*),jd(*),ix(*)
      real*8  ct(3,*),x(*),f(*),f0(*),u(*),dr(*),plo(10,*),ld(*)
      real*8  td(3)
      real*4  tary,tary1
      real*8, allocatable, dimension(:) :: temp

      data fa,tr/.false.,.true./
c.... transfer to correct process
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     + 23,24,25,26),j
c
c.... print stress values
c     [stre,,n1,n2,v3]      - output element values for n1 to n2 inc n3
c     [stre,all]            - output all element values
c     [stre,lay,n1,n2,-v3]  - output stresses at GP for layer v3=k.m ,inc = 1
c     [stre,node,n1,n2,v3]  - output nodal stresses n1 to n2 inc v3
c     [stre,node,n1,n2,-v3] - output nodal stresses for layer v3=k.m ,inc = 1
c     [stre,dmag,n1,n2,n3]  - output stresses at damaged gauss points
c     [stre,ints,n1,n2]     - output interlaminar stresses at center of element n1 with n2 intervals=n2+1 values
c     [stre,pris,n1,n2,n3]  - output principal stress n3 at nodes for nodes n1 to n2,inc = 1
c     [stre,line]           - output nodal stresses along line defined with SPLO
1     n1 = ct(1,l)
      n2 = ct(2,l)
      if(n2.eq.0) n2 = n1
      n3 = ct(3,l)
      if(n3.eq.0) n3 = 1
      klay = 0
      mlay = 1
      propq = prop ! set actual load
      if(rlnew.gt.0.d0) propq = prop*rlnew

      if(pcomp(lct(l),'node',4).or.pcomp(lct(l),'line',4)) then
        n1 = max(1,min(numnp,n1))
        n2 = max(1,min(numnp,n2))
c       if(n2-n1.ne.0) n3 = isign(n3,n2-n1)
        u_om = 0.d0 ! energy for ERRO
cww        klay = 0
cww        mlay = 1
        if(n3.lt.0) then
          klay = abs(n3)     ! number   of layer
          play = dabs(ct(3,l))
          play = (play-klay)*10
          mlay = nint(play) ! position in layer
          if(mlay.eq.0) mlay=1
          n3 = 1
        end if
        if (plfl) then
          dbgtxt = 'PMACR: [stre] gen. array: np,numnp*npstr*ipr(plfl)'
          call ralloc(strea,numnp*npstr,'STRE',plfl)
        end if
        n4 = 1 + numnp
        if(.not.fl(11).or.klay.ne.0) then
          istv = npstr - 1
          strea = 0.d0
          hflgu  = fa
          h3flgu = fa
          call formfe(u,dr,dr,dr,fa,fa,fa,fa,8,1,numel,1)
          call pltstr(strea,strea(n4),numnp)
        end if
        if(pcomp(lct(l),'node',4)) then
          call prtstr(x,strea,strea(n4),numnp,ndm,n1,n2,n3,1,iabs(istv),
     +                gtie,0)
        else if(pcomp(lct(l),'line',4)) then
          call prtstr(x,strea,strea(n4),numnp,ndm,n1,numnp,1,1
     +                ,iabs(istv),gtie,1)
c         call prtstrl(x,strea(n4),numnp,ndm,1,iabs(istv-1))
        end if
        fl(11) = tr
      else if(pcomp(lct(l),'pris',4)) then
        n1 = max(1,min(numnp,n1))
        n2 = max(1,min(numnp,n2))
        if (plfl) then
          dbgtxt = 'PMACR: [stre] gen. array: np,numnp*npstr*ipr(plfl)'
          call ralloc(strea,numnp*npstr,'STRE',plfl)
        end if
        n4 = 1 + numnp
        if(.not.fl(11).or.klay.ne.0) then
          istv = npstr - 1
          strea = 0.d0
          hflgu  = fa
          h3flgu = fa
          call formfe(u,dr,dr,dr,fa,fa,fa,fa,8,1,numel,1)
          call pltstr(strea,strea(n4),numnp)
        end if
        call prinstr(strea(n4),n3,numnp)
        call prtstr(x,strea,strea(n4),numnp,ndm,n1,n2,1,26,26,gtie,0)
        fl(11) = tr
      else  ! not node
        if (pcomp(lct(l),'all ',4)) then
          n1 = 1
          n2 = numel
          n3 = 1
        else if (pcomp(lct(l),'ints',4)) then
          n1  = max(1,min(numel,n1))
          n2  = n1
          n3  = 1
          intn = n1
          intv = ct(2,l)
          if(intv.eq.0) intv = 10
          ipv  = 1
          ispv = 1
        else
          n1 = max(1,min(numel,n1))
          n2 = max(1,min(numel,n2))
          if(n2-n1.ne.0) n3 = isign(n3,n2-n1)
          if (pcomp(lct(l),'lay ',3)) then
            klay = abs(n3)     ! number   of layer
            play = dabs(ct(3,l))
            play = (play-klay)*10
            mlay = nint(play) ! position in layer
            if(mlay.eq.0) mlay=1
            n3 = 1
          end if
          if (pcomp(lct(l),'dmag',4)) iprd = 1
        end if
        hflgu  = fa
        h3flgu = fa
        call formfe(u,dr,dr,dr,fa,fa,fa,fa,4,n1,n2,n3)
        intn = 0
        iprd = 0
cww     call rdealloc(strea,plfl) ! ??
      end if
      return
c
c.... form tangent stiffness
c     [utan]                     - form unsymmetric tangent
c     [utan,,1]                  -   " + form rhs and solve.
c     [utan,line,1,shift,value]  -   " with line search on value
2     if(neq.eq.0) return
      if(istyp.eq.1.or.istyp.eq.2.or.istyp.eq.8) 
     +  stop 'actual solver not allowed with UTAN'
      na = 1
      if(istyp.eq.0) then
        nau = na+neq
        if(neq.eq.1) nau=1
        nal   = nau+jd(neq)
        gsize = 2*jd(neq)+neq
      else ! other solver
        nau   = 1
        nal   = 1
        gsize = jd(neq+1)-1 
      end if
      if(fl(3)) then
        dbgtxt = 'PMACR:[UTAN gen.array: na,gsize (fl(4))'
        call ralloc(gstiff,gsize,'UTAN',fl(3))
      end if
      goto 330
c
c     [tang]                     - form symmetric tangent
c     [tang,,n1]                 - n1>=0: factor the tangent matrix
c     [tang,,1]                  -   " + form rhs and solve.
c     [tang,line,1,shift,value]  -   " with line search on value
3     if(neq.eq.0) return
      na = 1
      if(istyp.eq.0) then
        nau = na+neq
        if(neq.eq.1) nau=1
        nal   = nau
        gsize = jd(neq)+neq
      else ! other solver
        nau   = 1
        nal   = 1
        gsize = jd(neq+1)-1 
      end if
      if(fl(4)) then
        dbgtxt = 'PMACR:[tang] gen.array: na,gsize (fl(4))'
        call ralloc(gstiff,gsize,'TANG',fl(4))
      end if
c.... from now for TANG/UTAN
330   gstiff=0.d0
      cfr = j.eq.2
      f8o   = fa
      if(ct(1,l).gt.0.0d0) then
        f8o = tr
        prop1 = prop
        if(nop.eq.5) then ! calculate load at n+1/2 for energy momentum conserving algorithm
          ttim5 = ttim - 0.5d0*dt
          prop1 = propld(ttim5,0)
        end if
        call ploa1(ttim,rlnew,prop1,propq)
        call ploads(u,dr,propq,fa,tr,cfr) ! not possible with nop.eq.5
        call pload(id,f,f0,dr,nneq,propq)
      end if
      rfl = fa
      shift = ct(2,l)
      if(.not.fl(9).and.shift.ne.0.0d0) then
                     write(iow,2000) shift
        if(ior.lt.0) write(*  ,2000) shift
        call daxpty(nxl,-shift,massm,gstiff)
        if(cfr.and.nxl.gt.neq)call daxpty(nxl-neq,-shift,massm(nxu)
     +                                    ,gstiff(nal))
      end if
      hflgu  = f8o
      h3flgu = f8o
      call etimef(tary)
      if(pfr)                write(iow,2014) tary
      if(pfr .and. ior.lt.0) write(*  ,2014) tary

      call formfe(u,dr,gstiff,gstiff(nal),tr,f8o
     +           ,cfr,fa,3,1,numel,1)

      call etimef(tary1)
      if(pfr)                write(iow,2015) tary1-tary  ! time for K
      if(pfr .and. ior.lt.0) write(*  ,2015) tary1-tary
c...  move residual in kc for extended system
      if(extflg) call pmove(dr,extkc,neq)
c...  eigenvalue computation during iteration
      if(eigflg .and. nte.eq.1) then

          call cor(gstiff(nal),gstiff(nau),gstiff,eigk1,eigk2
     1           ,massm,jd,neq,nzykel,omegax,epseig)

      end if
c...  modify stiffness matrix and residuum: K = c1M+c2C+K G = G - c3 M a - c4 C v
      if(fl(9)) then
          call darray(dr,gstiff(nal),gstiff(nau),gstiff,massm(nxll)
     1          ,massm(nxu),massm,dampm(ncll),dampm(ncu),dampm,trans,jd
     2          ,nneq,tr,fl(1),f8o,cfr)
      end if

      call etimef(tary)

      it1=ct(1,lve(lv)) ! actual  iterat
      it2=ct(1,lvs(lv)) ! maximum iterat

      if(f8o) then
        rnorm   = sqrt(ddot(neq,dr,1,dr,1))
        resife2 = rnorm ! for FE2
                     write(iow,2001) rnorm,it1,it2,tary
        if(ior.lt.0) write(*  ,2001) rnorm,it1,it2,tary
      else
                   write(iow,2002) tary
        if(ior.lt.0) write(*,2002) tary
      end if
c.... factor the tangent matrix
      if(ct(1,l).ge.0.0d0) then
          call datri(gstiff(nal),gstiff(nau),gstiff,jd,neq,cfr)

        if(istyp.le.4) then ! not for iterative solvers PCG/PGMRES
          call etimef(tary1)
          if(pfr)                write(iow,2006) tary1-tary ! time for factorization
          if(pfr .and. ior.lt.0) write(*  ,2006) tary1-tary
          tary = tary1
        end if
      else
        call drawmess(
     +  'Unfact. tangent produced, do not try normal solution',1,0)
      end if
c.... solve the equations
      if(f8o) then
        fl(7) = fa
        fl(8) = fa
        if(abs(aengy).lt.aold) aold = abs(aengy)

          call dasol (gstiff(nal),gstiff(nau),gstiff,dr,jd,neq,aengy)


        call etimef(tary1)
        if(pfr)                write(iow,2016) tary1-tary
        if(pfr .and. ior.lt.0) write(*  ,2016) tary1-tary ! time for solv

c....   set initial values for a line search (conservative to
c....   check initital iterate solution for possible line search)
        if(ct(3,l).le.0.0d0) ct(3,l) = 0.8d0
        if (rnmax.eq.0.0d0) then
          rnmax = abs(aengy)
          aold  = rnmax/0.9999d0/ct(3,l)
        end if

        if(pfr)              write(iow,2004) rnmax,aengy,tol
        if(pfr.and.ior.lt.0) write(*  ,2004) rnmax,aengy,tol
        rnormeps = 1.e-13  ! ww to stop if first iteration step is nearly 0
        if(abs(aengy).le.tol*rnmax .and. .not.linear
cww     1     .or. rnorm.eq.0.d0 ) then
     1     .or. rnorm.le.rnormeps ) then
          ct(1,lve(lv)) = ct(1,lvs(lv))
          l = lve(lv) - 1
          nte = 1
        else if(.not.linear) then
          nte = 0
          if(pcomp(lct(l),'line',4)) then
c....       line search - linear search along solution direction
            if(abs(aengy).gt.ct(3,l)*aold) then
              allocate(temp(4*nneq))
              temp = 0.d0
              ml2=nneq+1
c....         set the initial step size to full value and perform search
              step = 1.0
              call serchl(aold,f0,f,id,temp,u,dr,ct(3,l),temp(ml2),
     1                    neq,nneq,step,pfr)
              deallocate(temp)
            end if
          end if
        end if
        if (arcf) then
          if(arcfs) call stepcntl(dr,m,ct,1,ic)
          call arclen(u,dr,arclm1,arclm2,f,f0,gstiff(nal),gstiff(nau),
     +                gstiff,jd,id,nneq,neq,ttim,rlnew)
          if(arcfs) call stepcntl(dr,m,ct,2,ic)
          if(arcfs.and. ic.eq.1) then              ! do restart
            if(npld.gt.0) prop = propld(ttim,0)
            call update(id,f0,f,u,trans,dr,nneq,neq,fl(9),pfr,3)
            !call reshis(ix(nen+1),nen1,numel,1, 2)
            gh2 = gh1
            if(.not.clfl) call pmove(crit(nc3),crit(nc4),numnp*2)
            return
          end if
        end if
c....   extended system for computing stability points
        if(extflg) then
          if(kex.ne.0 .and. eflg) then
            call scaleh(extkh,neq,kex)
            eflg = fa
          end if

          call ext(dr,extkh,extkdh,extkc,extkd,extkz1,extkz2,extke,f,f0
     1             ,gstiff(nal),gstiff(nau),gstiff,u,jd,id,nneq,rlnew)
        end if
        call update(id,f0,f,u,trans,dr,nneq,neq,fl(9),pfr,2)
c....   control d.o.f. without mass:  not necessary WW/FG 14.10.03
cww        if(fl(9)) call vazero(trans,massm,neq,nneq,nrt)
      end if
      return
c
c.... form out of balance force for time step/iteration
c     [form]       - form rhs residual
c     [form,acel]  -    " + get initial acceleration if needed
c     [form,expl]  -    " + do explicit solution with lumped mass
c     [form,ener]  - calculate external energy  U^T*F e.g. = V^T int N^T*q dA for plates
4     if(fl(8))then          ! residuum calculated then solv
        if(nop.lt.7) return  ! residuum and solv for explicit
      end if
      rfl = fa
      prop1 = prop
      if(nop.eq.5) then ! calculate load at n+1/2 for energy momentum conserving algorithm
        ttim5 = ttim - 0.5d0*dt
        prop1 = propld(ttim5,0)
      end if
c.... load terms
      call pzero(dr,nneq)
      call ploa1(ttim,rlnew,prop1,propq)
      call ploads(u,dr,propq,fa,fa,fa)! not possible with nop.eq.5
      call pload(id,f,f0,dr,nneq,propq)
c
      if(pcomp(lct(l),'ener',4)) then
c....   calculate external energy U^T*F e.g. = V^T int N^T*q dA for plates,   load vector is compressed
        rnorm = 0.d0
        do  nr = 1,nneq
          jr = id(nr)
          if (jr.gt.0) rnorm = rnorm + u(nr) * dr(jr)
        end do
                     write(iow,2017) rnorm
        if(ior.lt.0) write(*  ,2017) rnorm
        return
      end if
      hflgu  = tr
      h3flgu = tr
c     Compute current residual
      call formfe(u,dr,dr,dr,fa,tr,fa,fa,6,1,numel,1)

      call etimef(tary)
      rnorm = sqrt(ddot(neq,dr,1,dr,1))
cww   if(pfr .and. ior.lt.0) write(*  ,2003) rnorm
      if(prnt.and.ior.lt.0) write(*  ,2001) rnorm,it1,it2,tary
      if(prnt.and.ior.lt.0) write(iow,2001) rnorm,it1,it2,tary
c...  extended system, move residual in kc
      if(extflg) call pmove(dr,extkc,neq)

c...  modify  residuum:  G = G - c3 M a - c4 C v
      if(fl(9)) then
          call darray(dr,gstiff(nal),gstiff(nau),gstiff,massm(nxll)
     1      ,massm(nxu),massm,dampm(ncll),dampm(ncu),dampm,trans,jd
     2      ,nneq,fa,fl(1),tr,fa)
      end if

c     current residual norm
      if(fl(9)) then ! added ww
        rnorm = sqrt(ddot(neq,dr,1,dr,1))
        if(nop.eq.7.or.nop.eq.8) then
          if(pfr .and. ior.lt.0) write(*  ,2013) rnorm,ttim,prop
        else
                       write(iow,2003) rnorm
          if(ior.lt.0) write(*  ,2003) rnorm
        end if
      end if ! added ww
      fl(8) = tr
      if(pcomp(lct(l),'acel',4)) then   ! a_0 = M^-1*(R_ext-R_int-D*v_0)
c       Compute initial acceleration
c....   initialize for consistent mass
        if(fl(1).and.fl(9)) then
          if(istyp.eq.0) then 
            gsize = jd(neq)+neq
            nau = 1 + neq
            if(neq.eq.1) nau=1
          else
            gsize = jd(neq+1)-1
            nau = 1
          end if
          if(fl(4)) then
            dbgtxt = 'PMACR: [form] gen. array: na,gsize (fl(4))'
            call ralloc(gstiff,gsize,'FORM',fl(4))
          end if
          gstiff=0.d0
          call pmove(massm,gstiff,gsize)
          call pmove(dr,trans(nw),neq)
          call datri(gstiff(nau),gstiff(nau),gstiff,jd,neq,fa)
          call dasol(gstiff(nau),gstiff(nau),gstiff,trans(nw)
     1              ,jd,neq,rnorm)
c....     initialize for lumped mass
        else if(fl(2).and.fl(9)) then
          call piacel(massm,dr,trans(nw),neq)
        else
c....     else write warning
          call drawmess(
     +    'Unable to comp. start. acceleration, No mass matrix defined',
     +     1,0)
        end if
      end if
      if(pcomp(lct(l),'expl',4)) then
c.....  Explicit solution, Perform solution and update
        if(.not.fl(2).and..not.fl(9)) then !
          call drawmess('Expl. integration: lumped mass matrix needed',
     +     1,0)
          return
        end if
c.....  solve
        call pisolv(massm,dampm,nc,dr,dr,neq)
c.....  update of arrays
        call update(id,f0,f,u,trans,dr,nneq,neq,fl(9),pfr,2)
      end if
      return
c
c.... form a lumped mass approximation
c     [lmas]
5     if(fl(5)) then
        dbgtxt = 'PMACR: [lmas] gen. array: nl,neq(fl(5))'
        call ralloc(massm,neq,'LMAS',fl(5))
      end if
      imtyp = 1
      nx  = 1 ! nl
      nxl = neq
      nxu = 1  + neq      
      go to 660
c
c.... form a consistent mass approximation
c     [umas]   form unsymmetric matrix
19    if(istyp.eq.1.or.istyp.eq.2.or.istyp.eq.8) 
     +  stop 'actual solver not allowed with UMAS'
      imtyp = 1
      if(istyp.eq.0) then
        gsize  = 2*jd(neq)+neq
        nxu  = neq + 1   
        if(neq.eq.1) nxu=1
        nxll = neq + jd(neq) + 1
      else
        gsize  = jd(neq+1)-1
        nxu  = 1   
        nxll = 1
      end if 
      if(fl(6)) then
        dbgtxt = 'PMACR: [umas] gen. array: nxll,gsize(fl(6))'
        call ralloc(massm,gsize,'UMAS',fl(6))
      end if
      goto 660 
c
c     [cmas]
6     if((istyp.eq.1.or.istyp.eq.2).and.fl(4)) then !necessary for profil!
        call drawmess('Compute one TANG before using CMAS.',1,0)
        return
      end if
      imtyp = 1
      if(istyp.eq.0) then
        gsize  = jd(neq)+neq
        nxu  = neq + 1
        if(neq.eq.1) nxu=1
        nxll = nxu  
      else
        gsize  = jd(neq+1)-1
        nxu  = 1
        nxll = 1  
      end if
      nxl = neq + jd(neq)

      if(fl(6)) then
        dbgtxt = 'PMACR: [cmas] gen. array: nm,gsize(fl(6))'
        call ralloc(massm,gsize,'CMAS',fl(6))
      end if

660   massm = 0.d0
      cfrm = j.eq.19 ! symmetry/unsymmetry in formfe
      fl(1) = j.ne.5     ! .not. lmas
      fl(2) = j.eq.5     !       lmas
      if(fl(2)) nxll = 1 ! in case of lmas
      if(fl(2)) nxu  = 1 ! in case of lmas
      hflgu  = fa
      h3flgu = fa
      call formfe(u,massm,massm,massm(nxll),fl(1),fl(2),cfrm,fa,5,
     +            1,numel,1)
c     check if mass terms occur
      if(j.eq.5) tmas = ddot(neq,massm,1,massm,1)
      if(istyp.eq.0) then
        if(j.eq.6) tmas=ddot(jd(neq)+neq,massm,1,massm,1)
        if(j.eq.19)tmas=ddot(jd(neq)+neq,massm,1,massm,1)
     +                 +ddot(jd(neq),massm(nxll),1,massm(nxll),1)
      else
        gsize  = jd(neq+1)-1
        if(j.eq.6) tmas=ddot(gsize,massm,1,massm,1)
        if(j.eq.19)tmas=ddot(gsize,massm,1,massm,1)
      end if
      tmas = sqrt(tmas)/neq
      epstmas = 1.e-14
      if(tmas.lt.epstmas) call drawmess(
     + 'No mass terms calculated or mass/dof<1.e-14, check CMAS/LMAS ',
     +  1,0)
      return
c
c.... compute reactions and print
c     in case of dynamic analysis: associated terms are added
c     -> sum R .ne. 0 due to bcs.!
c     [reac,,n1,n2,n3] - print reactions at nodes n1 to n2 inc n3
c     [reac,all]       - print all reactions
c     [reac,eps,n1,v2] - print reactions if abs(R(idf)).gt.tol
c                        n1=idf, v2=tol
c     [reac,sigq,v]    - calculate stresses from reactions
c     print always value for reaction force from RSUM
c
c-------------------------------------------------------
c      t=0:         loads from mate + single loads
c      after tang:  reactions
c-------------------------------------------------------
c
7      if ( pcomp(lct(l),'all ',4)) then
        n1   = 1
        n2   = numnp
        n3   = 1
        ieps = 0
      else if (pcomp(lct(l),'eps ',4)) then
        n1   = 1
        n2   = numnp
        n3   = 1
        irdf = ct(1,l)
        reps = ct(2,l)
        ieps = 1
      else if ( pcomp(lct(l),'sigq',4)) then
        n1   = 1
        n2   = numnp
        n3   = 1
        ieps = 0
      else
        n1   = abs(ct(1,l))
        n1   = max(1,min(numnp,n1))
        n2   = ct(2,l)
        if(n2.eq.0) n2 = n1
        n2   = max(1,min(numnp,n2))
        n3   = ct(3,l)
        if(n3.eq.0) n3 = 1
        if(n2-n1.ne.0) n3 = isign(n3,n2-n1)
        ieps = 0
      end if
cww   if(.not.rfl) then ! calculate not always
        call pzero(dr,nneq)
        if(fldyn) then
          dbgtxt = 'PMACR: [reac] gen. array: nneq*ipr,fldyn)'
          call ralloc(dynrea,nneq,'REAC',fldyn)
        end if
c....   calculate nodal loads
        dynrea = 0.d0
c....   load factor
        call ploa1(ttim,rlnew,prop,propq)
c....   dynamic loads F_dyn
        if(fl(9)) call ploadd(u,dynrea,massm(nxll),massm(nxu),massm
     +        ,dampm(ncll),dampm(ncu),dampm,trans,jd,neq,nneq,fl(1),1)
c....   external loads from LOAD
        call pload(id,f,f0,dynrea,nneq,propq)
c....   move to full vector to include b.c.
        call pmovec(id,dynrea,dr,nneq)
        hflgu  = fa
        h3flgu = fa
c....   external loads from SLOA/QLOA  (uncompressed,dfl=true)
        call ploads(u,dr,propq,tr,fa,fa)
c....   internal loads F_int (uncompressed,dfl=true)
        call formfe(u,dr,dr,dr,fa,tr,fa,tr,6,1,numel,1)
cww   end if
      if(npldf.eq.1) then     ! needed for [tplo,reac]
        flreac = tr
        jreac = ipl(1,1)
        react = dr(jreac)
        if(npldf.eq.1) then     ! needed for [tplo,def]
          call pzero (reacc,60)
          do 100 i =1,10        ! store reactions of defined nodes
            inode = noden(i)    ! write to disk in sr storpl1
            if(inode .eq.0) go to 100
            in1 = ndf*(inode-1) ! array-index in -dr- of node inode
            ndfmax = min(6,ndf) ! only up to 6 dofs
            do ii = 1,ndfmax
              reacc(i,ii) = dr(in1+ii)   ! reactions of node inode
            end do
100       continue
        end if
      end if
c
      if ( pcomp(lct(l),'sigq',4)) then 
c       calculate SIGQ
        call prtreas(x,dr,ndm,ndf,n1,n2,n3,.not.fl(4))
      else 
c       print reactions
        call prtrea(dr,ndf,n1,n2,n3,.not.fl(4),irdf,ieps,reps,gtie)
c       calculate/print sum of reactions
        if(nfs1.gt.0) call prtrsum(react,irpt,dr,ndf,ndfrs,nfs1,pfr,iow)
      end if
c
c-----------------------------------------------
cww      rfl = tr
c
cww      ee = ddot(nneq,dr,1,u,1)
cww   write(iow,2005) ee
cww   if(ior.lt.0) write(*,2005) ee

      return
c
c.... check mesh for input errors
c     [chec]
8     hflgu  = fa
      h3flgu = fa
      call formfe(u,dr,dr,dr,fa,fa,fa,fa,2,1,numel,1)
      return
c
c.... [erro,,n1] calculate and print distribution of errors
c..... with respect to 'n1' percent of energy (default n1=eval = 5 perc.)
c....  error--norms/values 
c....  1 energy-norm
c....  2 L_2-norm
c....  3 Y0 - only Elmt05    
c.... [erro,save] saves data on file fres_err
9     if(fl(11)) then
        call ralloc(e_ome,numerr*numel,'REMESH-erron/new',tflb1)! also remesh1
        iet(1) = 2          ! switch in element for 'macro>erro,,' mode
        eval = ct(1,l)
        if(eval.eq.0.0) eval = 5.
        ioerr = 0
        do i=1,3 ! numerr                         ! number of error indicators
          e_om(i)  = 0.0d0                        ! domain error
          e_bar(i) = eval/100.*sqrt(u_om(i)/numel)! u_om calc isw=8
          if(e_bar(i).eq.0.d0) e_bar(i)=1.d0
        end do
        if(pcomp(lct(l),'save',4)) then
          ioerr = 1
          fint1 = fres
          call addext(fint1,'err ')
          open(22,file=fint1,form='formatted',status='unknown')
          rewind(22)
        end if
        hflgu  = fa
        h3flgu = fa
        call formfe(u,dr,dr,dr,fa,fa,fa,fa,9,1,numel,1) ! error/element
        if(ioerr.eq.1) then
          close(22)
          ioerr = 0
        end if
        call prterr(numel) ! global values 
      else
        call drawmess('Compute nodal stresses before error',1,0)
      end if
      return
c
c.... [damp]      form consistent damping matrix (default)
c     [damp,cons] form consistent damping matrix
c     [damp,lump] form lumped damping matrix
c     [damp,ucon] form unsymmetric consistent damping matrix  (only for standard solver possible)
10    if(lct(l).eq.'lump') then
        if (flgc) then
          dbgtxt = 'PMACR: [damp,lump] gen. array: nc,neq (flgc)'
          call ralloc(dampm,neq,'DAMP,LUMP',flgc)
          nc=2
        end if
        flgd = tr
        dflg = fa
        flgda= fa
        flgdu= fa
      else if(lct(l).eq.'    '.or.lct(l).eq.'cons') then
        if(istyp.eq.0) then 
          gsize = jd(neq)+neq
          ncl   = jd(neq)+neq
          ncu   = 1+neq
          if(neq.eq.1) ncu=1
        else 
          gsize = jd(neq+1)-1
          ncl   = 1
          ncu   = 1
        end if         
        if (flgc) then
          dbgtxt='PMACR:[damp,cons] gen.array: nc,gsize (flgc)'
          call ralloc(dampm,gsize,'DAMP,CONS',flgc)
          nc=2
        end if
        ncll = ncu
        flgd = fa
        dflg = tr
        flgda= tr
        flgdu= fa  ! symmetry/unsymmetry in formfe
      else if(lct(l).eq.'ucon') then
        if(istyp.eq.1.or.istyp.eq.2.or.istyp.eq.8)
     +    stop 'actual solver not allowed with UDAMP'
        if(istyp.eq.0) then 
          gsize=2*jd(neq)+neq 
          ncl   = jd(neq)+neq
          ncu   = 1+neq      
          ncll  = ncu+jd(neq)
        else
          gsize=jd(neq+1)-1 
          ncl   = 1
          ncu   = 1
          ncll  = 1
        end if

        if (flgc) then
          dbgtxt='PMACR: [damp,ucon] gen. array: nc,gsize (flgc)'
          call ralloc(dampm,gsize,'DAMP,UCON',flgc)
          nc=2
        end if
        flgd = fa
        dflg = tr
        flgda= tr
        flgdu= tr ! symmetry/unsymmetry in formfe
      end if
c
      dampm = 0.d0
      hflgu  = fa
      h3flgu = fa
      call formfe(u,dampm,dampm,dampm(ncll),dflg,flgd,flgdu,fa,12,
     +            1,numel,1)
      return
c
c --- [augm,,v1,n1] perform nested update for augmented lagrangian
c              v1 - factor for augmentation
c              n1 - flag for constraint residual
11    if(rnmax.eq.0.0d0) then
c.....  new time step
        augf = 1.0d0
      else
c.....  continue with current time step - multiply by factor
        if(ct(1,l).gt.0.0d0) augf = augf*ct(1,l)
      end if
      n1 = ct(2,l)
      if(n1.ne.0) cplus = 0.d0
      hflgu  = tr
      h3flgu = tr
      call formfe(u,dr,dr,dr,fa,fa,fa,fa,10,1,numel,1)
      aold  = rnmax
      aengy = rnmax
      if(n1.ne.0) then
        cplus = sqrt(cplus)
                     write(*  ,2011) cplus
        if(ior.lt.0) write(iow,2011) cplus
      end if
      return
c
c     [geom] geometric stiffness formulation K_T - K_0 for eigenvalues
c     influence of solvers is not really checked!!!
12    imtyp = 2
      if(fl(4)) then
       call drawmess(
     + 'calculate Tangent matrix with [tang,,1] to have displacements]',
     +    1,0)
          return
      end if
      if(istyp.eq.0) then
        gsize = jd(neq)+neq
        nxl = neq + jd(neq)
        nxu = 1 + neq
        if(neq.eq.1) nxu=1
      else
        gsize = jd(neq+1)-1
        nxl = 1
        nxu = 1
      end if
      if(fl(6)) then
        dbgtxt = 'PMACR: [geom] gen. array: nm,gsize(fl(6))'
        call ralloc(massm,gsize,'GEOM',fl(6))
      end if
      fl(1)  = tr
      fl(2)  = fa
      hflgu  = fa
      h3flgu = fa
c.... calculate tangent stiffness matrix K_T in massm only sym!
      gstiff = 0.d0
      call pzero(dr,nneq)
      call ploa1(ttim,rlnew,prop1,propq)
c.... for terms with influence on KT: SLOA, QLOA with TEMP
      call ploads(u,dr,propq,fa,tr,fa) ! not possible with nop.eq.5  .. in na,nal!!
c     call pload(id,f,f0,dr,nneq,propq) ! not necessary 
      call formfe(u,dr,gstiff,gstiff(nal),tr,fa,fa,fa,3,1,numel,1)
c...  move na->nm
      call pmove(gstiff,massm,gsize)
c.... calculate linear stiffness matrix  K_0 = K_T with u=0 and lambda=0 in gstiff
      gstiff = 0.d0
c.... temporary array for u (5 if ext.sys)
      allocate(temp(5*nneq))
      temp=0.d0
      props = prop  ! set all load dependent terms to zero, e.g. temperature  -> K_0!
      prop = 0.d0
      call formfe(temp,dr,gstiff,gstiff(nal),tr,fa,fa,fa,3,1,numel,1)
      prop = props
c.... calculate geometrical matrix  -(K_U+K_G) = K_0-K_T 
      call matmin(gstiff,massm,gsize,1,massm)
c.... factor the linear stiffness matrix
      call datri (gstiff(nal),gstiff(nau),gstiff,jd,neq,fa)
      call drawmess(
     + 'Tangent matrix destroyed, recompute with TANG after SUBS/LAN',
     +    1,-2)
      deallocate(temp)
c
      return
c
c.... Set parameter for Pre-Cconditioned bi-Conjugated gradient method PBCG
c     [pbcg,    ,niter,tol]
c     [pbcg,iter,niter,tol]
c
13    itercg = ct(1,l)
      if(itercg.eq.0) itercg = 150
      tolcg=ct(2,l)
      if(tolcg.eq.0.d0) tolcg = 1.e-8
      itolcg=ct(3,l)
      itolcg=min(max(1,itolcg),4)
      return
c
c.... compute J-integral, material forces
c     [jint,,n1,n2,n3] - print material forces at nodes n1 to n2 inc n3
c     [jint,all]       - print all material forces
c     [jint,eps,n1,v2] - print material forces if abs(R(idf)).gt.tol
c                        n1=idf, v2=tol
14      if ( pcomp(lct(l),'all ',4)) then
        n1   = 1
        n2   = numnp
        n3   = 1
        ieps = 0
      else if (pcomp(lct(l),'eps ',4)) then
        n1 = 1
        n2 = numnp
        n3 = 1
        irdf = ct(1,l)
        reps = ct(2,l)
        ieps = 1
      else
        n1   = abs(ct(1,l))
        n1   = max(1,min(numnp,n1))
        n2   = ct(2,l)
        if(n2.eq.0) n2 = n1
        n2   = max(1,min(numnp,n2))
        n3   = ct(3,l)
        if(n3.eq.0) n3 = 1
        if(n2-n1.ne.0) n3 = isign(n3,n2-n1)
        ieps = 0
      end if
      call pzero(dr,nneq)
      hflgu  = fa
      h3flgu = fa
      call formfe(u,dr,dr,dr,fa,tr,fa,tr,16,1,numel,1)
      call prtjint(dr,ndf,n1,n2,n3,fa,irdf,ieps,reps)
c........ old version together with mesh>jint
cwwc     [jint,,n1]   :n1 = 1 for full problem, n1 = 2 for symmetry
cww14    n1 = ct(1,l)
cww      n1 = max(n1,1)
cww      fact = n1
cww      rint = 0.d0
cww      call cjint(u,dr,ajint,njint)
cww      rint = rint/fact
cww                   write(*  ,2007) rint
cww      if(ior.lt.0) write(iow,2007) rint
      return
c
c.... calculate fracture criterion at nodes for all layer boundaries
c     [crit]
c     [crit,prin,n1] - print fracture criterion
c                n1=1 -> ics/f   n1=2 ->f(i)i=1,ncs
c     crit(nc1) = dt(numnp)
c     crit(nc2) = st(numnp,ncs) = fail
c     crit(nc3) = crit1(numnp,2)  (io/fo)   converged values
c     crit(nc4) = crit2(numnp,2)  (io/fo)   actual    values
15    ipc=0
      if (pcomp(lct(l),'prin',4)) then
        ipc = 1
        if(ct(1,l).eq.2.) ipc = 2
      end if
      if (clfl) then
        dbgtxt = 'PMACR: [CRIT] gen. array: 5*numnp+ncs*numnp,clfl'
        call ralloc(crit,5*numnp+ncs*numnp,'Crit frac',clfl)
        nc1 = 1
        nc2 = nc1 + numnp
        nc3 = nc2 + numnp*ncs
        nc4 = nc3 + numnp*2
        clfl= .false.
c....   initial values for crack boundary  (ibc=0 )
        crit = 0.d0
      end if
c.... initial values for crack boundary  (fbc=0)
      call pzero (crit,numnp*(ncs+1))
c.... calculate stresses at nodes and layer boundaries
      hflgu  = fa
      h3flgu = fa
      call formfe(u,dr,dr,dr,fa,fa,fa,fa,17,1,numel,1)
cww      call pcktie(id,crit,crit(nc2),ndf,numnp)
c.... average of f at nodes
      call pltstr1(crit,crit(nc2),numnp,ncs)
c.... failure criterion at nodes and layer boundaries
      call pccrit(crit(nc2),crit(nc3),crit(nc4),numnp)
c.... print failure criterion at nodes for all layer boundaries
      if(ipc.gt.0) call prtcrit(crit(nc2),crit(nc4),numnp,ncs,ipc)
c.... initialize crack opening displacements
      call pinitv(crit(nc3),crit(nc4),u,ndf,numnp)
      return
c
c.... show current macro parameters and infos on nodes or elements
c     [show]
c     [show,coor,n1,n2,n3]
c     [show,boun,n1,n2,n3]
c     [show,forc,n1,n2,n3]
c     [show,elem,n1,n2,n3]
c     [show,node,n1,n2,n3]
c     [show,memo]
c     [show,pset]
16    n1=ct(1,l)
      n2=ct(2,l)
      n3=ct(3,l)
      if (pcomp(lct(l),'coor',4)) then
        call showcoor(x,ndm,numnp,n1,n2,n3)
      else if (pcomp(lct(l),'boun',4)) then
        call showboun(id,ndf,numnp,n1,n2,n3)
      else if (pcomp(lct(l),'forc',4)) then
        call showforc(f,ndf,numnp,n1,n2,n3)
      else if (pcomp(lct(l),'elem',4)) then
        call showelem(ix,nen,nen1,numel,n1,n2,n3)
      else if (pcomp(lct(l),'node',4)) then
        call shownode(x,id,f,ndm,ndf,numnp,n1)
      else if (pcomp(lct(l),'memo',4)) then
        call pseta_out(1)
        call pmemo_out(1)
      else
        call plshowm(arcf,rlnew)
      end if
      return
c
c.... calculate determinant of stiffness matrix
c     [detk,    ]  calculate determinant of K_T
c     [detk,init]  calculate determinant and set as initial value
c                  K_T must be factored!
17    continue
      if (pcomp(lct(l),'init',4)) det0=0.d0
      call detkt(gstiff,neq,0)

      write(iow,2008) detc,nneg
      write(*  ,2008) detc,nneg
cc    if(ior.lt.0.and.pfr) write(*  ,2008) detc,nneg
      return
c
c.... summarize results from previous calculations to use plot macros
c     [summ,file,fac]  - add disp,stre with 'fac'from file named 'file'
c     [summ,zero    ]  - zero summarized displ./stresses
c
18    fsum = ct(1,l)
      if(fsum.eq.0.d0) fsum=1.d0
c.... array for stresses to add
      if(flsum) then
        dbgtxt = 'PMACR: [summ] gen. array: nsum,numnp*npstr*ipr(flsum)'
        call ralloc(summ,numnp*npstr,'SUMM',flsum)
      end if
c.... original stress array, if not initialized
      npp = numnp*npstr
      if (plfl) then
        dbgtxt = 'PMACR: [stre] gen. array: np,numnp*npstr*ipr(plfl)'
        call ralloc(strea,numnp*npstr,'SUMM-STRE',plfl)
        fl(11) = tr  ! to prevent calc. of stresses
      end if
      call summariz(lct(l),u,dr,strea,summ,nneq,npp,fsum)
      return
c
c19   is before cmas(=6)
c
c.... AUTO dynamic time step control, ### only possible at END of time step  ###
c     [auto]
c     [auto,load]
c         calculate a typical value of load P_max
c         then htol should be chosen: h_tol=0.1*P_max
c     [auto,init,htol,dtmax]
c         htol:  Typical value of load or residual force in problem, def=0.10
c         dtmax: maximum timestep,   def=dt
c     [auto,para,dtdo,dtup1,dtup2]
c         dtdo (D_A):  decreasing factor, def=0.85
c         dtup1(D_G):  increasing factor, def=0.80
c         dtup2(D_M):  increasing factor, def=1.25
c         decreasing of time step:   dt_new  = dtdo *htol/hres*dt_old
c         increasing of time step:   dt_new1 = dtup1*htol/hres*dt_old
c                                    dt_new2 = dtup2*          dt_old
c                                    dt_new  = min(dt_new1,dt_new2)
c         itimax = 3   number of steps to start increasing of time step
c         nploc  = 0/1 plot/do not plot solution point
c
c         valid only for Newmark, HHT and Gen.Alpha
c
20    if(pcomp(lct(l),'impl',4)) then
        if(ttim.lt.1.002d0*dt) return
        tm = 0.d0
        call formfe(u,b,a,c,.false.,.false.,.false.,.false.,19,
     1              1,numel,1)
        if(ct(1,l).ne.0.d0)then
          eta = ct(1,l)
        else
          eta = 1.1d0
        end if
        if(ct(2,l).ne.0.d0)then
          xi = ct(2,l)
        else
          xi = 0.005d0
        end if
        if(ct(3,l).ne.0.d0) then
          red = ct(3,l)
        else
          red = 0.5d0
        end if
        tcrit = xi/tm*dt
        tcrit = dsqrt(tcrit)
        if(tcrit.gt.dt) then
          tcrit = min(tcrit,tstart)
          tcrit = min(tcrit,eta*dt)
        end if
        if(tcrit.le.dt*red)then
          dtnew = dt*0.5d0
          goto 102
        else
          dt = tcrit
        end if
        return
      end if


      if(.not.fl(9)) then
        call drawmess(
     +  'Problem not dynamic - no time step control possible',1,0)
        return
      end if
      if(nop.ne.1.and.nop.ne.3.and.nop.ne.4) then
        call drawmess(
     +  'time step control only possible for Newmark,HHT,Gen.Alpha',1,0)
        return
      end if
      if(fldyn) then
        dbgtxt = 'PMACR: [auto] gen. array: nneq*ipr,fldyn)'
        call ralloc(dynrea,nneq,'AUTO',fldyn)
      end if

c
      if(pcomp(lct(l),'load',4)) then
c....   calculate nodal loads for maximum entry
        call ploa1(ttim,rlnew,prop1,propq)
        call ploads(dynrea,dr,propq,fa,fa,fa)
        call pload(id,f,f0,dr,nneq,prop1*rlnew)
c       norm of load vector
        rdnorm = sqrt(ddot(nneq,dr,1,dr,1))
c....   calculate maximum entry
        call minmax(dr,nneq,hresmin,hresmax)
        hres = max(abs(hresmin),abs(hresmax))
        write(*,*) ' |P_n|     ',rdnorm
        write(*,*) ' |P_n|_max ',hres
        write(*,*) ' choose T=0.1*|P_n|_max '
        return
      end if
c
      if(pcomp(lct(l),'init',4)) then
         htol   = ct(1,l)
         dtmax  = ct(2,l)
         if (htol .eq.0.0d0) htol  = 0.1d0
         if (dtmax.eq.0.0d0) dtmax = dt
         dtdo  = 0.85d0
         dtup1 = 0.80d0
         dtup2 = 1.25d0
         write(*,*) 'htol,dtmax       ',htol,dtmax
         write(*,*) 'dtdo,dtup1,dtup2 ',dtdo,dtup1,dtup2
         nploc = 0
         iti   = 0
         itimax= 3
         iaback= 0
         dtnew=dtmax
         return
      end if

      if(pcomp(lct(l),'para',4)) then
         dtdo  = ct(1,l)
         dtup1 = ct(2,l)
         dtup2 = ct(3,l)
         if(dtdo .eq.0.d0) dtdo  = 0.85d0
         if(dtup1.eq.0.d0) dtup1 = 0.80d0
         if(dtup2.eq.0.d0) dtup2 = 1.25d0
         write(*,*) 'dtdo,dtup1,dtup2 ',dtdo,dtup1,dtup2
         return
      end if

c.... calculate residual at n+1/2 of actual time step
c.... u at n+1/2
      call du_at_n12(u,dynrea,trans,id,nneq,1)
c.... time at n+1/2
      ttim5 = ttim - 0.5d0*dt
c.... calculate nodal loads
      hflgu  = fa
      h3flgu = fa
c...  external loads F_ext->nneq
      call ploa1(ttim5,rlnew,prop1,propq)
      call ploads(dynrea,dr,propq,fa,fa,fa)
      call pload(id,f,f0,dr,nneq,prop1*rlnew)
c.... internal loads F_int at n+1/2 ->neq
      call formfe(dynrea,dr,dr,dr,fa,tr,fa,fa,6,1,numel,1)
c.... du at n+1/2
      call du_at_n12(u,dynrea,trans,id,nneq,2)
c.... dynamic loads F_dyn at n+1/2
      call ploadd(dynrea,dr,massm(nxll),massm(nxu),massm,dampm(ncll)
     +           ,dampm(ncu),dampm,trans,jd,neq,nneq,fl(1),2)
c     current residual norm
      rdnorm = sqrt(ddot(neq,dr,1,dr,1))
c.... calculate maximum element of current residual
      call minmax(dr,neq,hresmin,hresmax)
      hres = max(abs(hresmin),abs(hresmax))
cww   write(*,*) ' |G_n+1/2|       ',rdnorm
cww   write(*,*) ' S=|G_n+1/2|_max ',hres
c.... hres <0.1*htol:   high accuracy
c.... hres =    htol:   good accuracy
c.... hres = 10*htol: coarse accuracy
c
c.... calculate length of new time step
      dhtol = htol/hres
c.... limitation
cww      dhtol = min(dhtol,10.d0)  ! max factor down = 10
cww      dhtol = max(dhtol,0.1d0)  ! max factor up   = 10

      if (hres.gt.htol) then
c....   decrease time step
        dtnew = dtdo*dhtol*dt
cww??   dtnew = min(dtnew,dtmax)
        nploc = 1 ! do not plot values for TPLO
        iti   = 0
      else
        if (0.75d0*htol.gt.hres) then
          iti=iti+1
          if (iti.ge.itimax) then
c....       increase time step if itimax times low errors
            dt1 = dtup1*htol/hres*dt
            dt2 = dtup2          *dt
            dt3 = dtmax
            dtnew  = min(dt1,dt2,dt3)
          else
c....       hold time step with low errors
            dtnew=dt
          end if
        else
c....     hold time step
          dtnew=dt
          iti=0
        end if
        nploc = 0 ! plot values for TPLO
      end if
 102  continue
      if(dtnew.ge.dt) then
c....   next time step with increased or constant time increment
        iaback=0
      else
c....   go back to time t_n and repeat time step with reduced time increment
        iaback=iaback+1
        ttim = ttim - dt
c       write(iow,2002) ttim,prop
c       if(ior.lt.0.and.prnt) write(*,2002) ttim,prop
        if(npld.gt.0) prop = propld(ttim,0)
        rnmax = 0.0d0
c....   reinitialize solution and dynamic vectors at time t_n
        call update(id,f0,f,u,trans,dr,nneq,neq,fl(9),pfr,3)
c....   reinitialize history vectors for step
        !call reshis(ix(nen+1),nen1,numel,1, 2)
        gh2=gh1
c....   reinitialize crack values at nodes
        if(.not.clfl) call pmove(crit(nc3),crit(nc4),numnp*2)
c....   reset basis vectors
        if(ldir.eq.1) then
          mdir1 = 1 + knode*10
          call pmove(basea,basea(mdir1),knode*10)
        end if
      end if
c.... reset dt
cww   write(*,*) 'time',ttim
      write(*,2018) iaback,dt,dtnew
      dt = dtnew
      return
c
c.... initial imperfections
c     [iimp]
21    hflgu  = tr
      h3flgu = tr
      viimpx = ct(1,l)
      viimpy = ct(2,l)
      call formfe(u,dr,dr,dr,fa,fa,fa,fa,20,1,numel,1)
      return

c.... Set parameter for General Minimum Residual method PGMRES
c     [pgmr,    ,niter,tol]
c     [pgmr,iter,niter,tol]
c
22    itergmr = ct(1,l)
      if(itergmr.eq.0) itergmr = 150
      tolgmr=ct(2,l)
      if(tolgmr.eq.0.d0) tolgmr = 1.e-8
      return
c
c.... Set parameter for Preconditioner for PBCG/PGMRES
c     [prco,    ,ip,dtol,lfil]
c     [prco,iter,ip,dtol,lfil]
c
23    ippc = ct(1,l)
      if(ippc.gt.5.or.ippc.lt.1) ippc = 3
      tolpc = ct(2,l)
      if(tolpc.eq.0.d0) tolpc = 1.e-3
      lfil = ct(3,l)
      if(lfil.eq.0) lfil = 5
      lfil  = min(lfil,neq)
      return
c
c.... Change solver 0-12, or calculate new profile-TESTVERSION
c     [nsys,solver
24    istyp = ct(1,l)
      call profupd(prt)
      return

c
c.... Update History Variables on MICRO-Level
c     [updh,n1]
c      n1 = 1  use macro on MACRO-Level
c      n1 = 2  use macro on MICRO-Level  update:  call reshis once
25    iupd = ct(1,l)
      if(iupd.eq.1) then ! on macro level
        hflgu  = tr
        h3flgu = tr
        write(*,*) 'updh,1'
        call formfe(u,dr,dr,dr,fa,fa,fa,fa,15,1,numel,1)
      else if(iupd.eq.2) then ! on micro level
        !call reshis(ix(nen+1),nen1,numel,2, 1)
        gh1=gh2
cww     write(*,*) 'updh,2'
      end if
      return

c.... Plot nodal stresses along line
c     [Splo,set]  set line
c     [Splo] reset line to zero
26    if(pcomp(lct(l),'set',3)) then
c...    set point A and B of line
        call pzero(dscor,6)
        call pzero(dsdcor,3)
        if(ior.lt.0) write(*,3001)
        call dinput(td,3)
        do i=1,ndm
          dscor(i,1)=td(i)
        end do
        if(ior.lt.0) write(*,3002)
        call dinput(td,3)
        do i=1,ndm
          dscor(i,2)=td(i)
        end do
        do i=1,ndm
          dsdcor(i)=dscor(i,2)-dscor(i,1)
        end do
        dsdcor2=dot(dsdcor,dsdcor,ndm)
        if(ior.lt.0) then
          write(*,2019)
          do i=1,2
            write(*,2020) (dscor(ii,i),ii=1,3)
          end do
        end if
        write(iow,2019)
        do i=1,2
          write(iow,2020) (dscor(ii,i),ii=1,3)
        end do
      else
        call pzero(dscor,6)
        dsdcor2=0.d0
      end if
      return
c
c.... formats
2000  format('   Shift to stiffness of ',1pe12.5)
cww2001  format('   Residual norm = ',1pe15.7,31x,'t=',0pf9.2)
2001  format('   Residual norm = ',1pe15.7,' it=',i4,'/',i4,
     +       18x,'t=',0pf9.2)
2002  format(65x,'t=',0pf9.2)
2003  format('   Residual norm = ',1pe15.7)
2004  format('   Energy convergence test'/'    Maximum   =',1pe24.15/
     1       '    Current   =',1pe24.15  /'    Tolerance =',1pe24.15)
cww2005  format(/'  Energy: Displacements * Reactions = ',1pe15.7/1x)
2006  format('   Time for triangular decomposition',29x,'t=',0pf9.2)
2007  format('   J-integral value =',e15.7)
2008  format('   Determinant = ',g12.5,3x,'neg. diagonals =',i3)
2011  format('   Residual of Constraint:  cp =',e15.7)
2013  format('   Residuum ',1pe15.7,' time ',e12.5,' load ',e12.5)
2014  format('   Begin    stiffness matrix',37x,'t=',0pf9.2)
2015  format('   Time for stiffness matrix',37x,'t=',0pf9.4)
2016  format('   Time for solution of equations   ',29x,'t=',0pf9.2)
2017  format('   External energy for influence areas = ',1pe15.7)
2018  format('   AUTO: Rerun(0/i),dt,dtnew',i2,1x,e12.5,1x,e12.5)
2019  format('   Coordinates for SPLO/DPLO')
2020  format(3x,3e14.7)
3001  format('   Input x,y,z of point 1 for SPLO/DPLO')
3002  format('   Input x,y,z of point 2 for SPLO/DPLO')
      end
c
      subroutine pmacr2(id,ix,f,f0,u,dr,lct,ct,ndf,nen1,nneq,prt,j,plo)
c-----------------------------------------------------------------------
c.... macro instruction subprogram 2                                   |
c
c     'tol ','dt  ','loop','next','prop','data','time','prin',
c     'nopr','mate','tran','init','iden','newf','back','debu',
c     'line','nonl','conv','rinp','if  ','else','endi','eas ',
c     'dibc','beta',
c
c-----------------------------------------------------------------------
      USE arcl
      USE augdat
      USE aunr
      USE back
      USE back1
      USE cdata
      USE conv
      USE ddata
      USE debugs
      USE dirdat
      USE dtauto
      USE endata
      USE epsd1
      USE errin2
      USE errin3
      USE evdata
      USE ext2
      USE fdata
      USE foutp
      USE hdata
      USE implstep
      USE iofile
      USE ldata
      USE macprt
      USE mate
      USE ndata
      USE ndatx
      USE pcrit
      USE pdata2
      USE plodf
      USE plodfa
      USE plodfb
      USE plodfs
      USE plodfu
      USE plong
      USE prlod
      USE psize
      USE rdata
      USE slid1
      USE slid2
      USE slid3
      USE slid4
      USE soltyp
      USE stepc
      USE strnam
      USE subdt
      USE tdata
      USE yydata
      USE doalloc
      implicit real*8 (a-h,o-z)
      logical fa,tr,prt,pcomp,err,lexpr,ldummy
      character*4 lct(*),lctl(2),expr*80
      integer id(*),ix(*)
      real*8  f0(*),f(*),u(*),dr(*),ct(3,*),plo(10,*)
      real*8  ctl(25) ! wg IF org=ctl(3)
      real*8  adr1(nneq)
      real*8, allocatable, dimension(:) :: temp
cww   common /easdata/ ieas
c
      save ! IF

      data fa,tr/.false.,.true./

c.... transfer to correct process
      go to (1,2,3,4,5,6,7,8,8,10,11,12,13,14,15,16,17,18,19,20,21,22,
     +       23,24,25,26), j
c
c.... set solution tolerance
c     [tol,,value]
1     tol = ct(1,l)
      return
c
c.... set time increment
c     [dt,,dtnew,factor]
c     dt = dt*factor  or dt = dtnew
2     dto=dt
      if(ct(2,l) .ne.0.0d0) then
        dt = dt*ct(2,l)
      else
        dt = ct(1,l)
      end if
      tstart = dt
      if(dto.eq.0.d0) dto=dt
cww>> für Dynamik-Bathe
      dtold=dt
cww<<
      return
c
c.... set loop start indicators
c     [loop,,number]
3     lv = lv + 1
      lvs(lv) = l
      lve(lv) = ct(2,l)
      ct(1,lve(lv)) = 1.
c.... set problem type if not the outer default loop.
      if(l.ne.1) then
        if(linear) then
          if(ior.lt.0) write(  *,2006)
                       write(iow,2006)
        end if
        linear = fa
      end if
      return
c
c.... loop terminator control
c     [next]
4     n = ct(2,l)
      ct(1,l) = ct(1,l) + 1.0
      if(ct(1,l).gt.ct(1,n)) lv = lv - 1
      if(ct(1,l).le.ct(1,n)) l = n
      return
c
c.... input proportional load table
c     [prop,,number] - input "number" proportional loads
5     npld = ct(1,l)
      if(npld.gt.10) stop 'maximum no. of cards(=10) for PROP reached.'
      npld = max(1,min(npld,10))
      prop = propld (ttim,npld)
      return
c
c.... data command
c     [data,tol] ;  [data,dt] ; [data,four]
6     if(ior.lt.0) write(*,3000) lct(l)
      call pintio(yyy,10)
      read(yyy,1000,err=61) lctl,ctl
      if(.not.pcomp(lct(l),lctl(1),4)) go to 401
      if(pcomp(lctl(1),'tol ',4)) tol = ctl(1)
      if(pcomp(lctl(1),'dt  ',4)) dt  = ctl(1)
      if(pcomp(lctl(1),'four',4)) then
        call pcfour(ctl(2),ctl(2),ctl(1),ctl(2),ctl(1),
     1              ndf,numel,numnp,1)
      end if
      return
61    call  errclr ('PMACR2')
      go to 6
c
c.... [time],,<tmax>...increment time , quit after ttim > tmax
c.... [time],start.....set start increment
7     continue
csk   store u(t) for BACK in case of divergence, included only for static analysis! 
      if(.not.allocated(ustore)) then 
        call ralloc(ustore,nneq,'Store u(t)',ldummy)
      end if
      do i=1,nneq
        ustore(i)=u(i)
      end do
cww   Dynamik-Bathe-Algorithmus
      if(nop.eq.6) then
      dt=dtold
        if(int(theta(2)).eq.1) then !switching between newmark and 3 point central difference scheme
          dt       = theta(1)*dt
          theta(2) = 2.3d0        !value only for switching 2<theta(2)<2.5, because of type real
        else if(int(theta(2)).eq.2) then
          dt       = (1-theta(1))*dt
          theta(2) = 1.3d0
        end if
      end if

      ttim = ttim + dt
      if(ct(1,l).gt.0.0d0 .and. ttim.ge.ct(1,l)) then
        ct(1,lve(lv)) = ct(1,lvs(lv))
      end if
c...  set flag (ok, if a solution will be performed afer TIME)
      backstep=fa
      if(npld.gt.0) prop = propld(ttim,0)
c.... store values for load deflection curve
c.... plo( 1,nstedf) = loadfactor, plo( 2,nstedf) = displacement
c.... plo( 3,nstedf) = time      , plo( 4,nstedf) = velocity
c.... plo( 5,nstedf) = reac      , plo( 6,nstedf) = determinant
c.... plo( 7,nstedf) = acce      , plo( 8,nstedf) = stress
c.... plo( 9,nstedf) = use1      , plo(10,nstedf) = use2
c.... pf = multiplier for symmetry
      if(npldf.eq.1) then
       nincp = nincp - 1
       if(nincp.eq.0) then ! plot in increments
        nincp = nincp0
        if(nstedf.eq.nplo) then
          call drawmess(
     +    'max, no. of  plot points reached, change NPLO',1,0)
          return
        end if
c....   plot in case of convergence
cww     if(abs(aengy).gt.tol*rnmax*1.d9) then
cww       write(*,2007)
cww       if(nploc.eq.1) return
cww     end if

c....   plot only in case of new time step from auto
        if(htol.ne.0.d0.and. nploc.eq.1)  go to 72

        nstedf = nstedf + 1
        if(nstedf.eq.1) then  ! initial values, if any
          jj= ipl(1,1)
          plo(2,nstedf) = u(jj)
          if(ipl(2,2).lt.0) plo(2,nstedf) = -u(jj)
          if(fl(9)) then
                         nu = 1 +   nneq ! at n+1
            if(nop.eq.5) nu = 1 + 2*nneq ! at n+1/2
            call storpl(trans,trans(nu),plo(1,1),nstedf,nneq,id,jj)
          end if
        end if
        if(nstedf.gt.1) then ! store values
          if(arcf.or.extflg) then
            plo(1,nstedf) = propld(ttim-dt,0)*rlnew*pf
          else
            plo(1,nstedf) = propld(ttim-dt,0)*pf
          end if
          plo(6,nstedf) = detc
          jj= ipl(1,1)
          plo(2,nstedf) = u(jj)
          if(ipl(2,2).lt.0) plo(2,nstedf) = -u(jj)
          plo(3,nstedf) = ttim - dt
          if(fl(9)) then
                         nu = 1 +   nneq ! at n+1
            if(nop.eq.5) nu = 1 + 2*nneq ! at n+1/2
            call storpl(trans,trans(nu),plo(1,1),nstedf,nneq,id,jj)
          end if
          if(flreac) then
            plo(5,nstedf) = react*pf
            if(npldf1.eq.0) flreac = fa  ! false if no ldd
          end if

          if(nstri.ne.0)
     +    call storps(strea(1+numnp),plo,numnp,nstedf,nstri,nstrno)

          plo( 9,nstedf) = valuse1
          plo(10,nstedf) = valuse2

        end if
72      continue
c....   write selected node values to file
        nu = 1
        if(fl(9)) then
                       nu = 1 +   nneq ! at n+1
          if(nop.eq.5) nu = 1 + 2*nneq ! at n+1/2
c.....    temporary use of allocatable temp for a
          allocate(temp(nneq))
          call pmovec(id,trans,     dr,nneq) ! copy velo.
          call pmovec(id,trans(nu),temp,nneq) ! copy acce.
        end if
        if(npldf1.eq.1)then
         call storpl1 (u,dr,temp,plo(1,1),nneq,nstedf,id,ndf,pf)
          flreac = fa  ! false here not at react
        end if
        if(fl(9)) deallocate(temp)
       end if
      end if
      if(nop.eq.7.or.nop.eq.8) go to 71
      if(prnt)              write(iow,2002) ttim,prop
      if(ior.lt.0.and.prnt) write(*  ,2002) ttim,prop
71    augf  = 1.0d0
      rnmax = 0.0d0
c.... divergence check    !c.kne >
      iconv  = 0
      iconv1 = 0
      ineg = 0
      do i=1,5
        ener(i)=0.d0
      end do               !c.kne<

c.... set start increment
      if(pcomp(lct(l),'star',4))
     +   call update(id,f0,f,u,trans,dr,nneq,neq,fl(9),pfr,4)
c
c.... implicit/explicit integration
c
      if(nop.ge.1.and.nop.le.8) then
c....   update dynamic vectors for time step
        if(fl(9)) then
          call dsetci
          call update(id,f0,f,u,trans,dr,nneq,neq,fl(9),pfr,1)
c....     control d.o.f. without mass:  not necessary WW/FG 14.10.03
cww          call vazero(trans,massm,neq,nneq,nrt)
        end if
c....   update frictional slideline arrays
        if(contfl) then
          call conupd(cl28,cl00,cl01,cl02,cl20,cl21,
     1                cl22,cl23,cl24,cl25,cl26,cl27,nsl)
        end if
      end if
c.... zero displacement increment for time step
      if(pcomp(lct(l),'star',4)) then
        call pzero(u(nneq*2+1),nneq)      ! DDu
      else
        call pzero(u(nneq+1),nneq+nneq) ! Du+DDu
      end if
      fl(10) = tr
c.... reset history variables
      !call reshis(ix(nen+1),nen1,numel,2, 1)
      gh1=gh2
c.... reset crack values at nodes
      if(.not.clfl) call pmove(crit(nc4),crit(nc3),numnp*2)
c.... update basis vectors
      if(ldir.eq.1) then
        mdir1 = 1 + knode*10
        call pmove(basea(mdir1),basea,knode*10)
      end if
c
      return
c
c.... set print flag for macro print outputs
c     [prin] 8 ; [nopr] 9
c     [prin,lmas] ; [prin,iden] ; [prin,cmas] ; [prin,geom] ;
c     [prin,resi] ; [prin,tang] ; [prin,utan] ; [prin,damp] ;
c     [prin,eigv]
c     [prin,on  ] ; [prin,off ]
c     for matrices:
c     diagonal entries for SOLV,0
c     otherwise: first n entries
8     if(pcomp(lct(l),'    ',4)) then
        pfr = j.eq.8
      else
        if(pcomp(lct(l),'off',3)) prnt = fa
        if(pcomp(lct(l),'on',2))  prnt = tr
        if(pcomp(lct(l),'lmas',4).or.pcomp(lct(l),'iden',4)) then
          call mprint(massm,1,neq,1,' l-mass/iden ' )
        else if(pcomp(lct(l),'cmas',4).or.pcomp(lct(l),'geom',4)) then
          call mprint(massm,1,neq,1,' c-mass/geom ' )
        else if(pcomp(lct(l),'damp',4)) then
          call mprint(dampm,1,neq,1,' damp ' )
        else if(pcomp(lct(l),'tang',4).or.pcomp(lct(l),'utan',4)) then
          call mprint(gstiff,1,neq,1,' tang/utan ' )
        else if(pcomp(lct(l),'resi',4)) then
          call mprint(dr,1,neq,1,' residuals ' )
        else if(pcomp(lct(l),'eigv',4)) then
          if(mfmax.ne.0) then
            do ii = 1,mfmax
              call prtom(eigd(ii),ii)
            end do
          end if
        end if
      end if
      return
c
cww.... set device for plots   actuell: PC VGA
cww     [tekt,type]
cww10    if(ior.lt.0) then
cww       write(*,3001)
cww101    read(*,1000,err=102,end=103) lct(l)
cww       go to  104
cww102    call  errclr ('PMACR2')
cww       go to  101
cww103    call  endclr ('PMACR2',lct(l))
cww      end if
cww104   idev = 3
cww      if(pcomp(lct(l),'IBM',3)) idev = 1
cww      if(pcomp(lct(l),'HP',2)) idev = 2
cww      if(pcomp(lct(l),'PC',2)) idev = 3
cww      if(pcomp(lct(l),'WIN',3)) idev = 4
cww      if(idev.eq.1) lct(l) = 'IBM'
cww      if(idev.eq.2) lct(l) = 'HP'
cww      if(idev.eq.3) lct(l) = 'PC'
cww      if(idev.eq.4) lct(l) = 'WIN'
cww      iclear = 0
cww                            write(iow,2001) idev,lct(l)
cww      if(ior.lt.0.and.prnt) write(*  ,2001) idev,lct(l)
cww      return
c

c     [mate,def]  define how to reset material to elements
c     [mate,new,v1,n2,n3]  update of material after [ERRO
c     based on critical value v1 and error norm n2 
c     [mate,org]  reset to orginal material
c     [mate,save] save orginal material
10    if(flmat) then
        call ialloc(matenew,numel,'Mate Array new',flmat)
        call ialloc(mateorg,numel,'Mate Array org',flmat)
        call matco3(mateorg,ix,nen1,numel,nen1) ! automatic save  
      end if
      if(pcomp(lct(l),'def',3)) then
c...    define rules for resetting og material      

      else if(pcomp(lct(l),'new',3)) then
c....   calculate errors using MACR>ERRO,v1 (n1=percentage) before
c...    calculate  material information new      
        v1=ct(1,l)
        if(v1.eq.0.d0) v1=1.d0
        n2=ct(2,l)
        n2 = max(1,min(n2,numerr))
        n3=ct(3,l)
        if(n3.eq.0) n3=2
        n3 = max(1,min(n3,nummat))
        call mateerr(matenew,ix,e_ome,numel,nen1,v1,n2,n3)  
c       set material new: ma at nen1=nen+4 in ix(nen1,numel) 
        call matco2(matenew,ix,nen1,numel,nen1) 

      else if(pcomp(lct(l),'save',4)) then
c       save material information
        call matco3(mateorg,ix,nen1,numel,nen1) 

      else if(pcomp(lct(l),'org',3)) then
c       reset material information
        call matco2(mateorg,ix,nen1,numel,nen1) 
      end if

      return
c
c.... input parameters and initialize vectors for dynamic analysis
c     [trans,,v1,v2] = [trans,newm,beta,gamma] 1 (Newmark)
c                      [trans,back,    ,     ] 2 (Backward-Euler)
c                      [trans,hht,alpha,     ] 3 (HHT)
c                      [trans,alph,rho ,     ] 4 (Generalized-Alpha)
c                      [trans,expl,    ,     ] 5 (Explizit-RLT)
c                      [trans,expa,    ,     ] 5 (Explizit-Abaqus) open: different length of time steps! no damping!
c                      [trans,ener,    ,     ] 6 (Energy Momentum)
c                      [trans,off,     ,     ] -  stop dynamic calculation
c     possible: switch between static and dynamic analysis(only 1 algorithm!)
c
11    nrt = 1  ! number of used vectors in array urate
      nrk = 0  ! solution vector,rate vector ur-nrk
      nrc = 0  ! rate vector ur-nrc
      nrm = 1  ! rate vector ur-nrm
      if(pcomp(lct(l),'off',3).and.fl(9)) then
        fl(9)=fa ! stop dynamic calculation
        return
      end if
      if(.not.fl(1).and..not.fl(2)) then ! mass matrix needed
        call drawmess('Compute mass matrix before use TRANS!',1,0)
        return
      end if
      call dparam(ct(1,l),lct(l))  ! set parameter for dynamic analysis
      if(fl(9)) return
c.... set arrays for the first time
      dbgtxt = 'PMACR: [trans] gen. array: nv,nneq*nrt*ipr(fl(9))'
      call ralloc(trans,nneq*nrt,'TRANS',fl(9))
      nw = 1 + (nrm-1)*nneq
      fl(9) = tr  ! flag for dynamic analysis
      return
c
c.... input initial conditions for dynamic integration
c     [init,disp] for displacements
c     [init,rate] for velocities
c     [init,acce] for accelerations, see also FORM,ACEL!
12    if(fl(9).and.pcomp(lct(l),'rate',4)) then
        if(nop.eq.8) then
          nvv = 1 + 4*nneq ! for expl Abaqus
        else
          nvv=1
        end if
        call invec(ndf,id,trans(nvv),' rate   ',prt)
      end if
      if(fl(9).and.pcomp(lct(l),'acce',4)) then
        nu = 1 + nneq
        call invec(ndf,id,trans(nu),' acce   ',prt)
      end if
      if(pcomp(lct(l),'disp',4))
     1   call genvec(ndf,u,'displacement',prt,err,tr)
      return
c
c.... define an identity vector for stiffness eigen computation
c     [iden,,n1,n2]  - set dof n1 to n2 to unity
13    if(fl(5)) then
        dbgtxt = 'PMACR: [iden] gen. array: nl,neq(fl(5))'
        call ralloc(massm,neq,'IDEN',fl(5))
      end if
      fl(1) = fa
      fl(2) = tr
      imtyp = 1
      nx    = 1 ! nl
      nxl   = neq
      nxll  = 1 
      nxu   = 1 
      massm = 0.d0
      n = ct(1,l)
      n = max(1,(n-1)*ndf+1)
      i = ct(2,l)*ndf
      if(i.eq.0) i = neq
      call piden(massm,n,i)
      return
c
c.... update the current force vector f0
c     [newf]
14    call daxpty(nneq,prop,f,f0)
      return
c
c.... backup a time step
c     [back,,<dtnew>] - back-up to beginning of time step <and reset dt>.
15    if(backstep) then
                              write(iow,2014) ttim
        if(ior.lt.0.and.prnt) write(  *,2014) ttim

      else
        dtnew = ct(1,l)
cww>> für Dynamik-Bathe
        if(nop.eq.6) then
          if(dtnew.ne.0.d0) then
            dtold=dtnew
          end if
        end if
cww<<
        ttim  = ttim - dt
        if(npld.gt.0) prop = propld(ttim,0)
                              write(iow,2002) ttim,prop
        if(ior.lt.0.and.prnt) write(  *,2002) ttim,prop
        rnmax = 0.0d0
c....   reinitialize solution and dynamic vectors for step
        call update(id,f0,f,u,trans,dr,nneq,neq,fl(9),pfr,3)
c....   reinitialize history vectors for step
        !call reshis(ix(nen+1),nen1,numel,1, 2)
        gh2=gh1
c....   reinitialize crack values at nodes
        if(.not.clfl) call pmove(crit(nc3),crit(nc4),numnp*2)
c....   reset basis vectors
        if(ldir.eq.1) then
          mdir1 = 1 + knode*10
          call pmove(basea,basea(mdir1),knode*10)
        end if
c....   reset dt
        if(dtnew.ne.0.d0) dt = dtnew
c....   set flag
        backstep=tr
c....   for TPLO,  then reac,all necessary
cww     nstedf = nstedf - 1  ! if not point occurs twice
      end if
      return
c
c.... debug flag on/off
c     [debug,on] or [debug,off] or [debug]
c     [debug,par] open files for parallel debugging
c     names foutpar_0i i=1,8 units 51-58
c       direct print write(51+ip,*)  xxxx
c       with  ip = OMP_GET_THREAD_NUM()
c       print matrices with MPRINTP
c
16    if(pcomp(lct(l),'    ',4)) debug = 1
      if(pcomp(lct(l),'on',2))   debug = 1
      if(pcomp(lct(l),'off',3))  debug = 0
      if(pcomp(lct(l),'par',3)) then
        do idp=1,8
        write(foutpar(idp),'(a9,i1)') 'foutpar_0',idp
        open(unit=(50+idp),file=foutpar(idp),status='unknown')
      end do
      end if
      if(debug.eq.1) then
        if(ior.lt.0) write(  *,2003)
        if(ior.gt.0) write(iow,2003)
      else
        if(ior.lt.0) write(  *,2004)
        if(ior.gt.0) write(iow,2004)
      end if
      dbgtxt = ' '
      return
c
c.... linear problem - no test on convergence
c     [line]ar
17    if(.not.linear) then
        if(ior.lt.0) write(  *,2005)
                     write(iow,2005)
      end if
      linear = tr
      return
c
c.... non-linear problem - test on convergence
c     [nonl]inear
18    if(linear) then
        if(ior.lt.0) write(  *,2006)
                     write(iow,2006)
      end if
      linear = fa
      return
c
c.... check energy convergence within a iteration loop
c     [conv,loop,iconv  ] - stop current iteration loop if energy rises
c                           iconv times in a row
19    continue
      if(pcomp(lct(l),'loop',4)) then  ! check convergence inside loop
        if(iconv.eq.0) then
          iconv = ct(1,l)
          if(iconv.eq.0) iconv = 2    ! iconv default = 2
          if(iconv.gt.5) iconv = 5    ! iconv max = 5
            if(ior.lt.0)   write(*  ,5002) iconv
                           write(iow,5002) iconv
        end if
        if(iconv1 .eq.0 ) go to 1991 ! the energy in the first iteration
        if(iconv1 .lt.5 ) then      ! is not condidered
          ener(iconv1) = abs(aengy) ! (iconv1 is an iteration counter)
        else
          ener(5) = abs(aengy)
        end if
        irise = 0
        jj=min(5,iconv1)
        if(jj.gt.iconv) then       ! count how many times energy rises
          do ji=jj-iconv,jj-1
            if(ener(ji+1) .gt.ener(ji) ) irise = irise+1
          end do
        end if
        if(irise .ge. iconv) then           ! if energy rises
          ct(1,lve(lv)) = ct(1,lvs(lv))     ! stop current loop
                       write(iow,5000) ener
          if(ior.lt.0) write(*  ,5000) ener
          return
        else
          if(iconv1.ge.5) then       ! energy not rising keep on going
            do i=1,4
              ener(i) = ener(i+1)
            end do
            ener(5) = 0.d0
          end if
        end if
1991    iconv1=iconv1+1
        return
      end if
c
c.... time step control based on convergence behaviour
c.... using module dtauto variables
c     [conv,aini,itinc,itdec,dtmin]
c     [conv,auto,xsiup,xidown]
c
      if(pcomp(lct(l),'aini',4)) then
          iaback = int(ct(1,l)+0.3d0) ! max number of iterations for increasing dt
          itimax = int(ct(2,l)+0.3d0) ! max number of iterations for constant dt
          dtdo   = ct(3,l)      ! minimum dt
          dtmax  = 2.d0*dt      ! maximum dt equals given dt
          if(iaback.eq.0) iaback = 4
          if(itimax.eq.0) itimax = 10
          if(dtdo.eq.0)   dtdo   = dtmax/1.d3
          icstop = 0
          return
      end if
      if(pcomp(lct(l),'auto',4)) then
          if(abs(aengy).gt.tol*rnmax) then
            xsi = ct(1,l)
            if(xsi.eq.0.d0) xsi = 0.75d0
            if(xsi*dt.le.dtdo) xsi=dtdo/dt
            if(ior.lt.0) write(*,  5010) xsi*dt
                         write(iow,5010) xsi*dt
            go to 1992       ! do back-up
          end if
          icstop = 0
          xsidec = ct(1,l)
          xsiinc = ct(2,l)
          if(xsidec.eq.0.d0) xsidec = 0.75d0 ! default values for decrease
          if(xsiinc.eq.0.d0) xsiinc = 1.5d0  ! increase
          if(it1.le.iaback) then
              dt = xsiinc*dt
              if(dt.gt.dtmax) dt=dtmax
              write(iow,5006) dt
              if(ior.lt.0) write(*,5006) dt
          end if
          if(it1.gt.itimax) then
              dt = xsidec*dt
              if(dt.le.dtdo) dt=dtdo
              write(iow,5007) dt
              if(ior.lt.0) write(*,5007) dt
          end if
          return
      end if
c
c.... auto Newton-Raphson load step control
c     [conv,aunr]
c
      if(pcomp(lct(l),'aunr',4)) then
        if(.not.allocated(adr)) then ! store first displacement rate
          call ralloc(adr,nneq,'AUNR',anr)
          adr(1:nneq) = u(nneq+1:2*nneq)/dt
          return
        end if

        if(abs(aengy).gt.tol*rnmax) then
          xsi = 0.25d0
          write(iow,5014) xsi*dt
          if(ior.lt.0) write(*,5014) xsi*dt
          go to 1992       ! do back-up and reduce dt
        end if

c........calculate norm
        dtol = ct(1,l)
        if(dtol.eq.0.d0) dtol = 1.d-3

        do i=1,nneq
          adr1(i)=adr(i)-u(nneq+i)/dt
          adr1(i)=adr1(i)*dt/2.d0
        end do
        eth   = sqrt(ddot(nneq,adr1,1,adr1,1))
        unorm = sqrt(ddot(nneq,u,1,u,1))
        rth=eth/unorm
c
c  Berechnung des  Maschinenrundungsfehlers FMACHP., siehe subroutine pivot
        FMACHP=1.D0
 1990   FMACHP=0.5D0*FMACHP
        IF(MASCHD(1.D0+FMACHP).EQ.1)  GOTO  1990
        FMACHP=FMACHP*2.D0

        rth = max(rth,fmachp)
        if(rth.gt.dtol) then
          xsi = 0.8d0*dsqrt(dtol/rth)
          xsi = max(0.1d0,xsi)
          write(iow,5013) xsi*dt
          if(ior.lt.0) write(*,5013) xsi*dt
          go to 1992       ! do back-up and reduce dt
        end if
        adr(1:nneq) = u(nneq+1:2*nneq)/dt
        xsi = 0.8d0*dsqrt(dtol/rth)
        xsi = min(xsi,2.d0)
        dt  = xsi*dt
        if(xsi.gt.1.d0) then !increase dt with respect to dtmax
          if(ct(3,l).gt.0.d0.and.dt.gt.ct(3,l)) dt=ct(3,l)
          write(iow,5011) dt
        if(ior.lt.0) write(*,5011) dt
        else
          if(dt.lt.ct(2,l)) then !reduce dt
            go to 1992       ! do back-up
          end if
          write(iow,5012) dt
          if(ior.lt.0) write(*,5012) dt
        end if
        return
      end if
c
c.... check convergence at the end of the iteration step
c     [conv,,xsi,dp,n3] - back-up to beginning of time step reset dt.
c                       - stop if load increment is lower than dp
c     [conv,stop,  ,n3] - stop current load step loop if no convergence
c                         or if number of negative diagonals at the end
c                         of the iter. step is changed (if n3 .ne.0)
      if(abs(aengy).gt.tol*rnmax) then   ! check convergence
        if(ior.lt.0) write(*,  2009)
                     write(iow,2009)
        xsi = ct(1,l)
        if(pcomp(lct(l),'stop',4)) then    ! stop current load-step loop
          ct(1,lve(lv)) = ct(1,lvs(lv))
          xsi=1.d0
        end if
        go to 1992       ! do back-up
      else
c....   see if number of negative diagonals changed (if n3 .ne. 0)
        inn = ct(3,l)
        if(inn .ne.0) then
          if(nneg.ne.ineg) then
                         write(iow,5003) ineg,nneg
            if(ior.lt.0) write(  *,5003) ineg,nneg
            if(pcomp(lct(l),'stop',4)) then  ! stop  load-step loop
                           write(iow,5004)
              if(ior.lt.0) write(  *,5004)
              ct(1,lve(lv)) = ct(1,lvs(lv))
              xsi = 1.d0
            else
                           write(iow,5004)
              if(ior.lt.0) write(  *,5004)
              xsi = ct(1,l)
            end if
            go to 1992    ! do back-up
          end if
        end if
      end if
      if(arcf) then
        dprop   = rlnew*prop-propold
        propold = rlnew*prop
        rold    = rlnew
      else
        dprop   = prop-propold
        propold = prop
      end if
      go to 1993
c.... go back one step   (like macro back)
1992  if(xsi.le.0.d0 .or. xsi.gt.1.d0) then
        if(arcf)       write(*,2010)
        if(.not.arcf)  write(*,2011)
        read (*,1001) xsi
      end if
      if(arcf) then
        dtnew = dt        ! dt = 1,  set ds0 new
        timold=timold-dt
        ds0 = ds0*xsi
c        rlnew = rold - 1.d0
        rlnew = rold
      else
        dtnew = dt*xsi    ! set dt new
      end if
      ttim = ttim - dt                   ! reset time
      if(npld.gt.0) prop = propld(ttim,0)
                            write(iow,2002) ttim,prop
      if(ior.lt.0.and.prnt) write(  *,2002) ttim,prop
      rnmax = 0.d0
c.... reinitialize solution and dynamic vectors for step
      call update(id,f0,f,u,trans,dr,nneq,neq,fl(9),pfr,3)
c.... reinitialize history vectors for step
      !call reshis(ix(nen+1),nen1,numel,1, 2)
      gh2=gh1
c.... reinitialize crack values at nodes
      if(.not.clfl) call pmove(crit(nc3),crit(nc4),numnp*2)
      dt = max(0.0d0,dtnew)
      if(pcomp(lct(l),'auto',4).and.dt.le.dtdo) then
          dt=dtmax ! if min dt reached, try again with max dt
          if(icstop.gt.0) then
              ct(1,lve(lv)) = ct(1,lvs(lv))    ! stop load step
              write(iow,5009)
              if(ior.lt.0) write(*,5009)
              return
          end if
          write(iow,5008) dt
          if(ior.lt.0) write(*,5008) dt
          icstop = 1
      end if
      if(npld.gt.0) propn = propld(ttim+dt,0)
      dprop = 2.d0*(propn-propold)
      if(arcf) then
c        dprop = 2.d0*rlnew*propn-propold
        return
      end if
      if(pcomp(lct(l),'auto',4)) return                  ! lower bound limited by dtdo
      if(pcomp(lct(l),'aunr',4)) then
          if(dt.lt.ct(2,l)) then
          write(iow,*) 'Reached given dtmin, stopping loop'
          if(ior.lt.0) write(*,*) 'Reached given dtmin, stopping loop'
          ct(1,lve(lv)) = ct(1,lvs(lv))    ! stop load step
            end if
          return
      end if
c.... stop if change in load-increment is lower than an given value
1993  continue
      if(dprop.lt.ct(2,l)) then
        ct(1,lve(lv)) = ct(1,lvs(lv))    ! stop load step
        if(arcf) then
          write(iow,5005)            dprop,ct(2,l),rlnew*prop
          if(ior.lt.0) write(*,5005) dprop,ct(2,l),rlnew*prop
        else
          write(iow,5001)            dprop,ct(2,l),prop
          if(ior.lt.0) write(*,5001) dprop,ct(2,l),prop
        end if
        return
      end if
      return
c
c.... [rinp]     - read input file again
c.... [rinp,new] - read new input file
c     nothing more to do here, macro is executed in pmacio
20    continue
      return
c
c     [IF,expression]   - begin if/then/else
c     [ELSE,expression] - begin if/then/else
21    li      = li + 1
      lie(li) = int(ct(3,l)) ! pos of end if
      lexpr   = fa
22    if(lexpr) then
        l = lie(li)
      else
        expr = lct(l)
        call evalex(expr,ctl,v1,80,errck)
        if(v1.lt.0.0d0) then
          l     = int(ct(2,l)) ! pos of else
        else
          lexpr = tr
        end if
      end if
      return
c
c     [end if] - end if if/then/else
23    li = li - 1
      return

c
c.... EAS parameter  - switch on,off(=0,1) default=off, not active
c     [eas]
c     [eas,n1]
c      n1 = 1    activate use of reshis
c      n1 = 2    update:    call reshis once
c      n1 = 3 de-activate use of reshis

24    ieas=ct(1,l)
      return
cww24    if(pcomp(lct(l),'    ',4)) ieas  = 0
cww      if(pcomp(lct(l),'on',2))   ieas  = 1
cww      if(pcomp(lct(l),'off',3))  ieas  = 0
cww      if(ieas.eq.1) then
cww        if(ior.lt.0) write(  *,2013)
cww        if(ior.gt.0) write(iow,2013)
cww      else
cww        if(ior.lt.0) write(  *,2012)
cww        if(ior.gt.0) write(iow,2012)
cww      end if
cww      return
c
c.... set values of displacement array
c     [dibc,,n1,n2,n3]
c     n1: material for which displacements have to changed
25    n1 = ct(1,l)
      n2 = ct(2,l)
      call setdibc(ix,u,numnp,numel,ndf,nen,nen1,n1,n2)
      return
c.... dynamic analysis name is obsolete
c     [beta
26    goto 11

c.... error diagnostics
401   write(iow,4001)
      if(ior.lt.0) write(*,4001)
      return
c.... formats
1000  format(a4,6x,a4,6x,3f10.0)
1001  format(f10.0)
cww2001  format(' Plotting device number is',i2,' or tektronix ',a4)
2002  format('   Computing solution for time',e14.6/
     1       '   Proportional load value is ',e14.6)
2003  format('   Debug flag is set to .true. - printing is on')
2004  format('   Debug flag is set to .false. - printing is off')
2005  format('   Linear Problem -- no test on convergence with TOL.')
2006  format('   Non-linear Problem - test on convergence with TOL.')
cww2007  format('   No convergence for TPLO in last time step reached')
2009  format('   Solution not converged! Choose new factor xsi')
2010  format('   for new arclength: dsnew = xsi*ds:  ',$)
2011  format('   for new time step: dtnew = xsi*dt:  ',$)
cww2012  format('   EAS-paramters are not used')
cww2013  format('   EAS-paramters are used')
2014  format('   Only one BACK-macro allowed'
     1       '   Computing time is still    ',e14.6)
3000  format(' Input ',a4,' Macro >',$)
cww3001  format(' Input Plot Device Type: (1) 4010; (2) 4662;',
cww     1            ' (3) uvax; (4) 4100; (5) vt200;',
cww     2            ' (6) sun ; (7) EGA ; (8) VGA  .'/3x,'>',$)
4001  format(' **ERROR** Macro label mismatch on data command')
5000  format(/,'  current iteration loop has been terminated',
     +         ' due to energy-divergence',/,
     +         '  energy of the last 5 iterations :',
     +            5(/,35x,g15.7) )
5001  format (/,'  Iteration stoped : ',
     +          '  load-increment ',g12.5,' lower than',g10.3,/
     +        23x,'last converged prop was ',g12.5)
5002  format(/,'  Iteration-loop will stop if energy rises',i2,' times')
5003  format(/,'*number of negative diagonals changed from',i4,' to',i4)
5004  format('**current loop was stopped and back one step (no time)')
5005  format (/,'  Iteration stoped : ',
     +          '  load-increment ',g12.5,' lower than',g10.3,/
     +        23x,'last converged loadvalue was ',g12.5)
5006  format(/,'dt increased to: ',g12.5)
5007  format(/,'  dt reduced to: ',g12.5)
5008  format(/,'Reached dtmin, try again with dtmax: ',g12.5)
5009  format(/,'!!!! Reached dtmin the second time, stopping loop.
     + No convergence possible with these settings!!!')
5010  format('   Solution not converged! Reducing dt to: ',g12.5 )
5011  format('AUTONR: Increasing dt to: ',g12.5)
5012  format('AUTONR: Decreasing dt to: ',g12.5)
5013  format('AUTONR: Doing backstep and Decreasing dt to: ',g12.5)
5014  format('AUTONR: Not converged! Doing backstep and Decreasing dt to
     + : ',g12.5)
      end
c
      subroutine pmacr3(ul,xl,tl,ld,p,s,ie,d,id,x,ix,f0,f,t,jd,u,dr
     1,lct,ct,ndf,ndm,nen1,nst,nneq,prt,j,plo,llreme)
c-----------------------------------------------------------------------
c.... macro instruction subprogram 3                                   |
c
c     'disp','solv','mesh','plot','subs','writ','read','cont',
c     'rest','copy','velo','acce','bfgs','arcl','save','four',
c     'fsol','fsum','paus','tplo','eigk','ueig','ext ','lamb',
c     'curv','reme','eigi','eig1','edit','post','smoo','tec ',
c     'parv','pola','lan ','sigq','epsq','dplo','feas',
c
c-----------------------------------------------------------------------
      USE arcext
      USE arcl
      USE cdat1
      USE cdata
      USE comfil
      USE ddata
      USE debugs
      USE dspos
      USE eig1
      USE endata
      USE epsdh
      USE evdata
      USE ext1
      USE ext2
      USE fdata
      USE fe2dat
      USE feapprog
      USE fileno
      USE fodata
      USE hdata
      USE hdatam
      USE iofile
      USE isbfgs
      USE isbfgs1
      USE iscsr
      USE ldata
      USE mdata
      USE mdat2
      USE ndata
      USE ndatx
      USE nolink
      USE pcrit
      USE pdata2
      USE pindex
      USE plodf
      USE plodfa
      USE plodfb
      USE plodfs
      USE plong
      USE pnodn
      USE prlod
      USE psize
      USE rdata
      USE slid1
      USE slid2
      USE slid3
      USE slid4
      USE soltyp 
      USE stepc
      USE strnam
      USE subdt
      USE tdata
      USE tplomax
      USE uneig
      USE doalloc

      implicit real*8 (a-h,o-z)
      logical prt,pcomp,sfl,accrcy,tfl
      logical lexst,ldummy,flgdyn2
      logical lopen
      logical fa,tr
      character fint*229,lct(*)*4,y*1
      character*12 plotit(10),fint1*229,hash*1,dum*6
      character*80 yyy,extens*1
      character*4 cnum
      integer ld(*),ie(*),id(*),ix(*),jd(*)
      real*8  ct(3,*),ul(*),xl(*),tl(*),p(*),s(*),
     1      d(*),x(*),f0(*),f(*),t(*),u(*),dr(*),uu(6),plo(10,*)
      real*8  td(3)
      real*8 dnoden(10)
      real*8, allocatable, dimension(:) :: evect,evat !temporary eigenvector, eigenvalue arrays
      real*8, allocatable, dimension(:) :: eigg, eigdp, eigdt, eigp
      real*8, allocatable, dimension(:) :: glamb ! temporary array gstiff
      real*8, allocatable, dimension(:) :: temp ! temporary array
      real*8, allocatable, dimension(:) :: temp1, temp2, temp3 ! temporary array
     +        , eigt1,eigt2, eigt3, eigdh, eiggh, eigvh, eigm 
     +        , eigv1,eigd1
      integer, allocatable, dimension(:) :: itemp1, itemp2, itemp3
      common /pswit/  imod

      save tfl
      data plotit /'#  lambda   ','   displ    ','   time     ',
     1             '   veloc    ','   reac     ','   det      ',
     2             '   acce     ','   stress   ',
     3             '   valuse1  ','   valuse2  '/

      data fa,tr/.false.,.true./
c.... transfer to correct process
      go to (1,2,3,4,5,6,7,8,9,10,1,1,13,14,15,16,17,18,19,20,21,22,
     1     23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39), j


c
c.... print displacements
c     [disp,cart]      - set output to cart. directions
c     [disp,pola,n1]   - set output to polar directions acting on n1=i/k
c     [disp,,n1,n2,n3] - print displ. n1 to n2 step n3
c     [disp,all]       - print all displ.
c     [disp,line]      - print all displ.
c     [disp,eigv,n1,n2,n3] - print eigvectors n1 to n2 step n3
c     [disp,evex,n1,n2,n3] - print eigenvector ext. sys. n1 to n2 step n3
c     [disp,eigi,n1,n2,n3] - print eigenvector inv. ite. n1 to n2 step n3
1     n1 = ct(1,l)
      n2 = ct(2,l)
      if(n2.eq.0) n2 = n1
      n3 = ct(3,l)
      if(n3.eq.0) n3 = 1
      if(pcomp(lct(l),'cart',4)) then
        write(*,2011)
        ipola = 0
        return
      else if(pcomp(lct(l),'pola',4)) then
        ipola = ct(1,l)
        if(ipola.eq.12.or.ipola.eq.13.or.ipola.eq.23) then
          write(*,2012) ipola
        else
          write(*,2011)
          ipola = 0
        end if
        return
      end if
      if(pcomp(lct(l),'eigv',4)) go to 110
      if(pcomp(lct(l),'evex',4)) go to 120
      if(pcomp(lct(l),'eigi',4)) go to 130
      tfl = pfr
      if(pcomp(lct(l),'all ',4)) then
        n1 = 1
        n2 = numnp
        n3 = 1
      else
        n1 = max(1,min(n1,numnp))
        n2 = max(1,min(numnp,n2))
        if(n2-n1.ne.0) n3 = isign(n3,n2-n1)
        if(pcomp(lct(l),'rad',3)) irad = 1
cww     pfr = tr
      end if
      propl = prop
      if(arcf) propl = propl*rlnew
      if(j.eq.1) then
        if(pcomp(lct(l),'line',4)) then
          call prtdis(x,u,bang,ttim,propl,numnp,ndm,ndf,1,numnp,1,1,
     +    ipola,1)
        else
          call prtdis(x,u,bang,ttim,propl,numnp,ndm,ndf,n1,n2,n3,1,
     +    ipola,0)
        end if
      end if
c.... print velocities
c     [velo,,n1,n2,n3] - print veloc. n1 to n3 step n3
c     [velo,all]       - print all veloc.
      if(j.eq.11) then
        if(fl(9)) then
          call pmovec(id,trans,dr,nneq)
          call prtdis(x,dr,bang,ttim,propl,numnp,ndm,ndf,n1,n2,n3,2,
     +    ipola,0)
        end if
        if(.not.fl(9))
     +    call drawmess('Problem not dynamic - no output produced',1,0)
      end if
c.... print accelerations
c     [acce,,n1,n2,n3] - print accel. n1 to n3 step n3
c     [acce,all]       - print all accel.
      if(j.eq.12) then
        if(nop.eq.2) then
          call drawmess('No Accelerations for Euler Backward',1,0)
          return
        end if
                     nu = 1 +   nneq ! at n+1
        if(nop.eq.5) nu = 1 + 2*nneq ! at n+1/2
        if(fl(9)) then
          call pmovec(id,trans(nu),dr,nneq)
          call prtdis(x,dr,bang,ttim,propl,numnp,ndm,ndf,n1,n2,n3,3,
     +    ipola,0)
        end if
        if(.not.fl(9))
     +    call drawmess('Problem not dynamic - no output produced',1,0)
      end if
      pfr = tfl
      return
c.... print eigenvectors  eigv
110   n2 = max(1,min(n2,mf2*mf))
      n1 = max0(1,min0(mf2*mf,n1))
      if(n2-n1.ne.0) n3 = isign(n3,n2-n1)
      call pzero(dr,nneq)
      do nnn=n1,n2,n3
        kk = (nnn-1)*neq
        call pmovec(id,eigv(1+kk),dr,nneq)
        call prtdis(x,dr,bang,ttim,eigd(nnn),numnp,ndm,ndf,1,numnp,1,1,
     +  ipola,0)
      end do
      return
c.... print eigenvector  evex
120   n1 = ct(1,l)
      n2 = ct(2,l)
      n3 = ct(3,l)
      if(n1.eq.0) n1 = 1
      if(n2.eq.0) n2 = numnp
      if(n3.eq.0) n3 = 1
      call pzero(dr,nneq)
      call pmovec(id,extkh,dr,nneq)
      call prtdis(x,dr,bang,ttim,prop*rlnew,numnp,ndm,ndf,n1,n2,n3,1,
     +ipola,0)
      return
c.... print eigenvector  eigi
130   n1 = ct(1,l)
      n2 = ct(2,l)
      n3 = ct(3,l)
      if(n1.eq.0) n1 = 1
      if(n2.eq.0) n2 = numnp
      if(n3.eq.0) n3 = 1
      call pzero(dr,nneq)
      call pmovec(id,eigia,dr,nneq)
      call prtdis(x,dr,bang,ttim,evi,numnp,ndm,ndf,n1,n2,n3,1,ipola,
     +0)
      return
c
c.... solve the equations
c     [solv]
c     [solv,line,value] - use line search for energy ratios > value
2     if(neq.eq.0) go to 201
      if(fl(4)) then
        call drawmess('No stiffness matrix, use tang or utang',1,0)
cww     if(ior.lt.0) return
cww     if(ior.ge.0) stop
        return
      end if
      if(.not.fl(8)) return
      fl(7) = fa
      fl(8) = fa
      if(abs(aengy).lt.aold) aold = abs(aengy)
      call dasol (gstiff(nal),gstiff(nau),gstiff,dr,jd,neq,aengy)
      if (rnmax.eq.0.0d0) then
        rnmax = abs(aengy)
        if(ct(1,l).le.0.0d0) ct(1,l) = 0.8d0
        aold = rnmax/ct(1,l)/0.9999d0
      end if
      if(pfr)              write(iow,2000) rnmax,aengy,tol
      if(pfr.and.ior.lt.0) write(*  ,2000) rnmax,aengy,tol
      if(abs(aengy).le.tol*rnmax.and. .not.linear) then
        ct(1,lve(lv)) = ct(1,lvs(lv))
        l = lve(lv) - 1
      else if(.not.linear) then
c....   line search
        if(pcomp(lct(l),'line',4).and.abs(aengy).gt.ct(1,l)*aold) then
          allocate(temp1(nneq))
          allocate(temp2(3*nneq))
          temp1 = 0.d0
          temp2 = 0.d0
c....     set the initial step size to full value and perform the search
          step = 1.d0
          call serchl(aold,f0,f,id,temp1,u,dr,ct(1,l),temp2,
     1         neq,nneq,step,pfr)
          deallocate(temp1)
          deallocate(temp2)
        end if
        if (arcf) then
          if(arcfs) call stepcntl(dr,m,ct,1,ic)
            call arclen(u,dr,arclm1,arclm2,f,f0,gstiff(nal),gstiff(nau),
     +                  gstiff,jd,id,nneq,neq,ttim,rlnew)
          if(arcfs) call stepcntl(dr,m,ct,2,ic)
          if(arcfs.and. ic.eq.1) then              ! do restart
            if(npld.gt.0) prop = propld(ttim,0)
            call update(id,f0,f,u,trans,dr,nneq,neq,fl(9),pfr,3)
            !call reshis(ix(nen+1),nen1,numel,1, 2)
            gh2=gh1
            if(.not.clfl) call pmove(crit(nc3),crit(nc4),numnp*2)
            return
          end if
        end if
c....   extended system for computing stability points
        if(extflg) then
          if(kex.ne.0 .and. eflg) then
            call scaleh(extkh,neq,kex)
            eflg = fa
          end if

          call ext(dr,extkh,extkdh,extkc,extkd,extkz1,extkz2,extke,
     1        f,f0,gstiff(nal),gstiff(nau),gstiff,u,jd,id,nneq,rlnew)
        end if
      end if
201   call update(id,f0,f,u,trans,dr,nneq,neq,fl(9),pfr,2)
c.... control d.o.f. without mass:  not necessary WW/FG 14.10.03
cww      if(fl(9)) call vazero(trans,massm,neq,nneq,nrt)
      return
c
c.... modify mesh data (cannot change profile of stiffness/mass)
c     [mesh]
3     i = -1
      call pmesh(ld,ie,d,id,x,ix,f,t,f0,ndd,nie,ndf,ndm,nen1,i,prt)
cww      if (i.gt.0) go to 400
      return
c
c.... plot macro
c     [plot] - enter plot mode for interactive
c     [plot,optn] - see plot manual for optn's
4     continue
      rfl = fa
      call pplotf(ul,xl,tl,ld,p,s,ie,d,id,x,ix,f,f0,t,jd,u,dr,
     1            lct(l),ct(1,l),prop,ndf,ndm,nen1,nst,plo)
      return
c
c.... subspace eigencomputations
c     [subs,,n1,n2,n3]  - subspace for n1 eigenpairs , n2 guard vectors
c     n1 no of eigenpairs
c     n2 guard vectors   
c     n3 < 0            --> tols = tol    nits = 25
c     n3 = 0 (default)  --> tols = 1.d-12 nits = 25
c     n3 > 0            --> tols = tol    nits = n3
c     [subs,init]         - copy eigenvector to COR-method
c     [subs,prin,n1,n2]   - subspace for n1 eigenpairs - print matrices
c                         - n2 used to overwrite default no. guard vects.
c     [subs,set,n1] set no. of printed plotted eigenvalues/vectors 
c                         - n1=mf2=1 save n1 values(default),
c                         - n2=mf2=2 save mq=min(n1+n1,n1+n8) values
c                           NOT DOCUMENTED
c
5     continue
c...  set no of stored EV
      if(pcomp(lct(l),'set',3)) then
        mf2 = ct(1,l)
        mf2 = max(1,min(mf2,2))
        return
      end if
c...  check data for EV-Problem 
      call evcheck(mb,imas,ieverror)
      if(ieverror.eq.1) return
c...  initialize COR-method macr: [eigk]
      if(pcomp(lct(l),'init',4)) then
        if(.not.allocated(eigv))then
          call drawmess('Compute eigenvectors first',1,0)
          return
        end if
        call ralloc(eigk1,neq,'EIGK',eigflg)
        call pmove(eigv,eigk1,neq)
        return
      end if
c...  set data
      allocate(eigm(size(massm)))
      eigm=massm
      mf = ct(1,l)
      mad= ct(2,l)
      if(mad.le.0) mad = 8
      mf = min(neq,max(1,mf))
      mq = min(mf+mf,mf+mad,neq)
      call numass(eigm,neq,mq,jd,imas)
      if(mq.lt.mf.and.ior.gt.0) write(iow,2001) mq
      if(mq.lt.mf.and.ior.lt.0) write(*  ,2001) mq
      if(mq.eq.0) then
        call drawmess('0 nonzero lumped mass terms, check CMAS',1,0)
        return
      end if
      mf = min(mf,mq)
      sfl = pcomp(lct(l),'prin',4)
c.... the following defs depend on the value of mq and are only temporary used
      allocate(eigg (mq*(mq+1)/2))
      allocate(eigdp(mq         ))
      allocate(eigdt(mq         ))
      allocate(eigp (mq*mq      ))
      allocate(eigt1(neq        ))
      allocate(eigdh(mq*(mq+1)/2))
      allocate(eigvh(mq*neq     ))
      allocate(eigv1(mq*neq     ))
      allocate(eigd1(mq         ))

      eigg  = 0.d0
      eigdp = 0.d0
      eigdt = 0.d0
      eigp  = 0.d0
      eigt1 = 0.d0
      eigdh = 0.d0
      eigvh = 0.d0
      eigv1 = 0.d0
      eigd1 = 0.d0

      tols = 1.d-12
      nits = 25
      if(ct(3,l).ne.0) tols = tol
      if(ct(3,l).gt.0) nits = ct(3,l)

c.... Solve EV-problem
      write(iow,2026) mf,mq,mq
      write(*,  2026) mf,mq,mq
      call subsp(gstiff,eigm,eigv1,eigt1,eigg,eigdh,eigd1,eigdp,
     +     eigdt,eigp,eigvh,jd,mf,mq,neq,imas,shift,tols,sfl,nits)


c.... store 'save' evs      
      mfmax=mf 
      if(mf2.ne.1)mfmax=mq   
       
      call ralloc(eigv,neq*mfmax,'SUBS-EIGVEC',sfl)
      call ralloc(eigd,    mfmax,'SUBS-EIGVAL',sfl)

      call pmove(eigv1,eigv,mfmax*neq)
      call pmove(eigd1,eigd,mfmax)

c.... deallocating temp arrays
      deallocate(eigg )
      deallocate(eigdp)
      deallocate(eigdt)
      deallocate(eigp )
      deallocate(eigt1)
      deallocate(eigdh)
      deallocate(eigvh)
      deallocate(eigm)
      deallocate(eigv1)
      deallocate(eigd1)

      return
c
c.... write a file
c     [writ,file]  - open write file named 'file'
c     [writ,disp]  - write displacements to 'file'
c     [writ,stre]  - write nodal streses to 'file'
c     [writ,wind]  - rewind 'file'
c     [writ,clos]  - close 'file'
6     call writer(lct(l),u,nneq)
      return
c
c.... read a file
c     [read,file]  - open read file named 'file'
c     [read,disp]  - read displacements from 'file'
c     [read,stre]  - read nodal streses from 'file'
c     [read,wind]  - rewind 'file'
c     [read,clos]  - close 'file'
7     call reader(lct(l),u,nneq)
      return
c
c.... set contact flag
c     [cont,,n,sfacn,sfact]
c.....       new penalty factor normal      sfacn
c.....       new penalty factor tangential  sfact
8     contfl = tr
      if(pcomp(lct(l),'off',3)) contfl = fa
c.... reset penalty if requested
      n = ct(1,l)
      if(n.gt.0) then
        ls = 1
        lm = 1
        if(n.gt.1) then
          do 81 i = 1,n-1
            ls = ls + cl00(i)
            lm = lm + cl01(i)
81        continue
        end if
        call setcpn(cl00(n),cl01(n),ct(2,l),ct(3,l),cl08,cl34,
     +       cl09,cl35,ls,lm)
      end if
      return
c
c.... restart previously run problem
c     [rest,ext,n1,n2]
c     ext = name of restart file fint_ext, ext=4 character
c     n1 > 0 restart file is an ascii file
c     n2   shared memory FE^2
9     fint = fres
      if(.not.pcomp(lct(l),'    ',4)) then
        cnum     = lct(l)
        call addext(fint,cnum)
      end if
      if(prt)                write(iow,2020) fint
      if(prt .and. ior.lt.0) write(*  ,2020) fint
      call restrt(fint,u,ix,ndm,ndf,nen1,nneq,1,ct(1,l),ct(2,l))
      return
c
c.... copy solution into force values! ( for use to combine displacements)
c     [copy]
10    call pmove(u,f,nneq)
      fl(7) = fa
      return
c
c 11  velo -> 1
c
c 12  acce -> 1
c
c.... BFGS algorithm ** January 1986
c     [bfgs,,n1,v2] - BFGS soln n1 = no. steps; n2 = soln tol
c.... n1 - Max. of of iterations (nbfgs) <def = 15>
c.... v2 - stol for line search 0.5-0.9  <def = .8>
c
c.... allocate memory for arrays
c.... mbo : vector oldrsd       in bfgs algo. (iterat)
c     mbd : vector d            in bfgs algo. (iterat)
c     mbv : vector v            in bfgs algo. (iterat)
c     mbw : vector w            in bfgs algo. (iterat)
c     mbt : vector t            in bfgs algo. (iterat)
c     mbs : 15 store vectors s  in bfgs algo. (iterat)
13    if (fl(12)) then
cww     dbgtxt to be added!
        call ralloc(bfgsbo,  nneq,'BFGS-BO',fl(12))
        call ralloc(bfgsbd,  nneq,'BFGS-BD',fl(12))
        call ralloc(bfgsbv,   neq,'BFGS-BV',fl(12))
        call ralloc(bfgsbw,   neq,'BFGS-BW',fl(12))
        call ralloc(bfgsbt,3*nneq,'BFGS-BT',fl(12))
c       15 max. vectors for v, w in BfGS (store)  
        call ralloc(bfgsbs,15*neq,'BFGS-BS',fl(12))
      end if
c
      nbfgs = ct(1,l)
      stol  = ct(2,l)
      if (nbfgs.eq.0)    nbfgs = 15
      if (stol.eq.0.0d0) stol  = 0.8d0
      nbfgs = min(nbfgs,15)
      stol  = min(stol,0.9d0)
      call iterat(gstiff,jd,u,dr,bfgsbo,bfgsbd,bfgsbt,accrcy,
     1     bfgsbv,bfgsbw,prt,f0,f,id,nbfgs,stol,nneq)
      return
c
c.... arc-length method
c     [arcl,,kflag,lflag]  - set arc length parameters
c     [arcl,,,1]           - modify arc length parameters
c     [arcl,modi]          - modify arc length parameters
c     [arcl,add,n1,tau]    - add eigvenvector n1, amount = tau
c     [arcl,impf,n1,xsi]   - imperf=eigvenvector n1, amount = xsi
c     [arcl,chec,n1]       - check for bifurcation using eigv. n1
c     [arcl,step,n1,n2,n3] - use step contol
14    if(pcomp(lct(l),'add',3))  go to 143
      if(pcomp(lct(l),'chec',4)) go to 144
      if(pcomp(lct(l),'impf',4)) go to 145
      if(pcomp(lct(l),'off',3))  go to 141
      if(pcomp(lct(l),'on',2))   go to 141
      if(pcomp(lct(l),'step',4)) go to 142
c.... check dt-value
      if(dt.ne.1.0d0) then
        dt = 1.0d0
        if(ior.lt.0) write(*,2007)
      end if
      if(kflag.eq.0) kflag = ct(1,l)
      if(kflag.eq.0) kflag = 2
      lflag = ct(2,l)
      if(pcomp(lct(l),'modi',4)) lflag=1
      if(lflag.eq.0) then
        n1 = max(neq,nneq)
        if(mu1.eq.1) call ralloc(arclm1,n1,'ARCL-mu1',arcf)
        mu1=2 
cww     if(arcf) return
        if(mu2 .eq. 1) then
c         if(kflag .eq. 2 .or. kflag .eq. 3 .or. kflag .eq. 5) then
          call ralloc(arclm2,n1,'ARCL-mu2',arcf)
          mu2=2
cww       if(arcf) return
c         end if
        end if
        arcf = tr
      end if
c.... modify parameter
      call dicont(id,numnp,ndf,lflag)
      return
c.... turn arc length method on and off
141   if(pcomp(lct(l),'off',3)) arcf  = fa
      if(pcomp(lct(l),'off',3)) arcfs = fa
      if(pcomp(lct(l),'on',2))  arcf  = tr
      return
c.... set parameter for step control
142   continue             !  (default values are set in sr ini)
      if(ct(1,l) .ne. 0.d0)  itd  = ct(1,l)
      if(ct(2,l) .ne. 0.d0)  cmax = ct(2,l)
      if(ct(3,l) .ne. 0.d0)  rm   = ct(3,l)
      arcfs = tr
      write(iow,2016) itd,cmax,sp,rm
      write(*,  2016) itd,cmax,sp,rm
c     if(ior.lt.0) write(*,  2016) itd,cmax,sp,rm
      return
c.... add scaled eigenvector to displacement vector
143   tau = ct(2,l)
      n1  = ct(1,l)
      n1 = max(min(mf2*mf,n1),1)
      if(.not.allocated(eigv).and.ior.lt.0 .and. .not. extflg) then
        call drawmess('Compute eigenvectors first',1,0)
        return
      end if
      write(*,2006)
      read(*,1001) nmeth
      nmeth = max(1,nmeth)
      if(nmeth.ne.1) then
        if(tau.eq.0.0d0) then
          kk = (n1 - 1) * neq
          vphi = dotid(u,eigv(1+kk),id,nneq)
          phi2 = ddot(neq,eigv(1+kk),1,eigv(1+kk),1)
          v2   = ddot(nneq,u,1,u,1)
          tau  = 100.d0 * vphi / sqrt(v2*phi2) + 1.d0
          if(ior.lt.0) write(  *,2003) tau,n1
                       write(iow,2003) tau,n1
        end if
        kk = (n1 - 1) * neq
        call paddv(u(1),eigv(1+kk),nneq,tau,id)
      else
        if(tau.eq.0.0d0) then
          vphi = dotid(u,extkh,id,nneq)
          phi2 = ddot(neq,extkh,1,extkh,1)
          v2   = ddot(nneq,u,1,u,1)
          tau  = 100.d0 * vphi / sqrt(v2*phi2) + 1.d0
          if(ior.lt.0) write(  *,2003) tau,n1
                       write(iow,2003) tau,n1
        end if
        call paddv(u,extkh,nneq,tau,id)
      end if
      return
c.... check for bifurcation or limit point (only point loads)
144   n1 = ct(1,l)
      if(n1.eq.0) n1 = 1
      if(.not.allocated(eigv).and.ior.lt.0) then
        call drawmess('Compute eigenvectors first',1,0)
        return
      end if
      kk = (n1 - 1) *neq
      chec = dotid(f,eigv(1+kk),id,nneq)
      if(abs(chec).lt.1.e-02) then
        if(ior.lt.0) write(*,2004) n1,chec  ! bifurcation point
        write(iow,2004) n1,chec
      else
        if(ior.lt.0) write(*,2014) n1,chec  ! limit point
        write(iow,2014) n1,chec
      end if
      return
c.... add eigenvector as imperfection to displacement vector
145   xsi = ct(2,l)
      n1  = ct(1,l)
      n1 = max(min(mf2*mf,n1),1)
      if(.not.allocated(eigv).and.ior.lt.0) then
        call drawmess('Compute eigenvectors first',1,0)
        return
      end if
      kk = (n1 - 1) * neq
      call paddi(u(1),eigv(1+kk),nneq,xsi,id)
      if(ior.lt.0) write(  *,2009) n1,xsi
                   write(iow,2009) n1,xsi
      return
c
c.... save restart information for intermediate points
c     [save,ext,n1,n2]
c     ext = name of save file fint.ext, ext=4 character (not init or next!)
c     n1 > 0 save file is an ascii file
c     n2   shared memory FE^2
c     write numbered save files e.g.for each load step
c     [save,init,n1]  initialize isno for Fint.isno
c     [save,next,n1]  icount=icount+1
c                     save on next file     Fint.isno
C
15    fint = fres
      if(.not.pcomp(lct(l),'    ',4)) then
        if(pcomp(lct(l),'init',4)) then
          isno = ct(1,l)
          call dochar2(fint,ipos1)
          ncol = ipos(fint,229)
          write(*  ,2019) fint(1:ncol),isno+1
          write(iow,2019) fint(1:ncol),isno+1
          return
        else if(pcomp(lct(l),'next',4)) then
          isno = isno+1
          if(isno.ge.10000) then
            call drawmess('Rest-Files only up to 10000 possible',1,0)
            return
          end if
          write(cnum,'(i4)') isno
          if(cnum(1:1) .eq. ' ') cnum(1:1) = '0'
          if(cnum(2:2) .eq. ' ') cnum(2:2) = '0'
          if(cnum(3:3) .eq. ' ') cnum(3:3) = '0'
          call addext(fint,cnum)
        else
          cnum      = lct(l)
          call addext(fint,cnum)
        end if
      end if
      ncol = ipos(fint,229)
      write(*  ,2018) fint(1:ncol)
      write(iow,2018) fint(1:ncol)
      call restrt(fint,u,ix,ndm,ndf,nen1,nneq,2,ct(1,l),ct(2,l))
      return
c
c.... set up for a fourier series solution
c     [four,,nf,factor]
16    if(fouflg) then
        dbgtxt = 'PMACR: [four] gen. array: ifour1,16*numel*ipr(fouflg)'
        call ralloc(aifour1,8*numel,'FOUR1',fouflg)
        call ralloc(aifour2,8*numel,'FOUR2',fouflg)
      end if

      call pcfour(u,dr,aifour1,aifour2,ct(1,l),ndf,numel,numnp,1)
      return
c
c.... save a fourier series solution for the current harmonic
c     [fsol]
c     [fsol,,factor]
17    fostr = tr
      nfs   = 0
      call pcfour(u,dr,aifour1,aifour2,ct(1,l),ndf,numel,numnp,2)
cww   hflgu/h3flgu for history terms?
      call formfe(u,dr,dr,dr,fa,fa,fa,fa,4,1,numel,1)
      call pcfour(u,dr,aifour1,aifour2,ct(1,l),ndf,  nfs,numnp,3)
      fostr = fa
      return
c
c.... sum the displacements and stresses at the angle ct(1,l)
c     [fsum,,theta]
18    call pcfour(u,dr,aifour1,aifour2,ct(1,l),ndf,numel,numnp,4)
      call pcfour(u,dr,aifour1,aifour2,ct(1,l),ndf,  nfs,numnp,5)
      return
c
c.... pause on no convergence
c     [paus]
19    if(abs(aengy).ge.100.d0*rnmax) then
        if(ior.lt.0) then
          write(*,2005)
          read(*,1002) y
          if(y.eq.'y' .or. y.eq.'Y') return
          l = lve(1) - 1
        else
          l = lve(1) - 1
        end if
      end if
      return
c
c.... specify the load deflection curve plot values
c...  [tplo,init,n1,n2,n3] initialize l.d.c: n1 = node, n2 = dof, n3 = mult
c...                                     n2>0 -> +u n2<0 -u
c...  [tplo,mark,n1,n2,n3] mark curve in color 'n1' with marker 'n2'
c                          with increment 'n3'
c                          stop with tplo,mark
c...  [tplo,fact,n1,n2] set values new: n1=+-on 1.axis n2 = fact on 2.axis
c...  [tplo,xval,v1,v2] set xmin,xmax for plot on screen
c...  [tplo,yval,v1,v2] set ymin,ymax for plot on screen
c...  [tplo,save,v1]    save load deflection data on   disk on file fres(_v1)
c...  [tplo,read,v1]    read load deflection data from disk on file fres(_v1)
c
c...  [tplo,load]  1  plot load         vs displacement
c...  [tplo,disp]  2  plot displacement vs time
c...  [tplo,phas]  3  plot velocity     vs displacement (phase portrait)
c...  [tplo,reac]  4  plot reaction     vs displacement, reac can be from RSUM
c...  [tplo,detd]  5  plot determiant   vs displacement
c...  [tplo,velo]  6  plot velocity     vs time
c...  [tplo,acce]  7  plot acceleration vs time
c...  [tplo,loat]  8  plot load         vs time
c...  [tplo,reat]  9  plot reaction     vs time, reac can be from RSUM
c...  [tplo,dett] 10  plot determinant  vs time
c...  [tplo,tdis] 11  plot time         vs displacement
c...  [tplo,strd] 12  plot stress       vs displacement
c...  [tplo,strt] 13  plot stress       vs time
c...  [tplo,us1d] 14  plot value use1   vs displacement
c...  [tplo,us2d] 15  plot value use2   vs displacement
c...  [tplo,us1t] 16  plot value use1   vs time
c...  [tplo,us2t] 17  plot value use2   vs time
c
c...  [tplo,def,n1,n2,n3] define nodes for multiple plot-curves
c...   if no of nodes.lt.3: nodes= n1,n2,n3
c...   else                 up to 10 nodes interatively
c...  [tplo,set,n1,n2,n3] set active plot curve similar to [tplo,init]
c...  [tplo,conv,n1]  0(1) plot(not) in case of no convergence
c...                  1 = default in start
c...  [tplo,incr,n1]  set plot increment to n1 (def=1)
c...  [tplo,sset,n1,n2] set plot stress n1 at node n2
c
c     ipl(1,1)= absolute dof number to plot
c     ipl(2,1)= node to plot
c     ipl(2,2)= dof to plot
c
cc20  if(ior.lt.0) then
20    iswpl = 0
      imod  = 1
      if(pcomp(lct(l),'load',4)) then
        iswpl = 1
      else if(pcomp(lct(l),'disp',4)) then
        iswpl = 2
      else if(pcomp(lct(l),'phas',4)) then
        iswpl = 3
      else if(pcomp(lct(l),'reac',4)) then
        iswpl = 4
      else if(pcomp(lct(l),'detd',4)) then
        iswpl = 5
      else if(pcomp(lct(l),'velo',4)) then
        iswpl = 6
      else if(pcomp(lct(l),'acce',4)) then
        iswpl = 7
      else if(pcomp(lct(l),'loat',4)) then
        iswpl = 8
      else if(pcomp(lct(l),'reat',4)) then
        iswpl = 9
      else if(pcomp(lct(l),'dett',4)) then
        iswpl = 10
      else if(pcomp(lct(l),'tdis',4)) then
        iswpl = 11
      else if(pcomp(lct(l),'strd',4)) then
        iswpl = 12
      else if(pcomp(lct(l),'strt',4)) then
        iswpl = 13
      else if(pcomp(lct(l),'us1d',4)) then
        iswpl = 14
      else if(pcomp(lct(l),'us2d',4)) then
        iswpl = 15
      else if(pcomp(lct(l),'us1t',4)) then
        iswpl = 16
      else if(pcomp(lct(l),'us2t',4)) then
        iswpl = 17
      end if
c...  plot values
      if(iswpl.gt.0. and. iswpl.le.17) then
        if(nincp.eq.1) then ! plot in increments with respect to tplo,incr
          call modscal(1)
                       nu = 1 +   nneq ! at n+1
          if(nop.eq.5) nu = 1 + 2*nneq ! at n+1/2
                       nu = min(nu,size(trans))
                       istpos = 1+numnp
                       istpos = min(istpos,size(strea))
          call plotdf(plo,ipl(2,1),ipl(2,2),nstedf,mkflg,mmc,mmst,incmk,
     1         iswpl,u,trans,trans(nu),strea(istpos),id,nneq,numnp)
          call modscal(2)
          return
        end if
      end if

      if(pcomp(lct(l),'mark',4)) then
        mkflg = 1
        mmc   = ct(1,l)
        if(mmc.eq.0) mkflg = 0
        mmc   = max(1,min(mmc,7))
        mmst  = ct(2,l)
        mmst  = max(0,min(mmst,6))
        incmk = ct(3,l)
        incmk = max(1,incmk)
        return
      end if

      if(pcomp(lct(l),'save',4)) then
        fint1 = fres
        in1=int(abs(ct(1,l)))         ! [tplo,save,in1]
        in1=min(in1,9)
        if(in1.ne.0) then
           write(extens,'(i1)') in1
           call dochar2(fint1,ipos1)
           do ii = 1,7
             kk = ipos1+ii-1
             if(fint1(kk:kk).eq.' ') goto 202
           end do
202        fint1(kk  :kk  ) = '_'
           fint1(kk+1:kk+1) = extens
        end if
        call addext(fint1,'ldf ')
        if(idev.lt.4) open(12,file=fint1,form='formatted',
     +                status='unknown')
        if(idev.eq.4) call opencf(12,fint1)
        write(12,1011) '#',nstedf,ipl(2,1),ipl(2,2),mkflg,mmc,mmst,
     +                 incmk,pf
        write(12,1004) plotit
        write(12,1003) ((plo(i,jj),i=1,10),jj=1,nstedf)
        close(12)
        if(pfr)                write(iow,2023) ipl(2,1),ipl(2,2),pf
        if(pfr .and. ior.lt.0) write(*  ,2023) ipl(2,1),ipl(2,2),pf
        return
      end if

      if(pcomp(lct(l),'read',4)) then
        fint1 = fres
        in1=int(abs(ct(1,l)))         ! [tplo,read,in1]
        in1=min(in1,9)
        if(in1.ne.0) then
          write(extens,'(i1)') in1
          call dochar2(fint1,ipos1)
          do ii = 1,7
            kk = ipos1+ii-1
            if(fint1(kk:kk).eq.' ') goto 203
          end do
203       fint1(kk  :kk  ) = '_'
          fint1(kk+1:kk+1) = extens
        end if
        call addext(fint1,'ldf ')
        inquire(file=fint1,exist=lexst)
        if(.not.lexst) then
          call drawmess('No LDF-file available on disk',1,0)
        else
          if(idev.lt.4) open(12,file=fint1,form='formatted',
     +                  status='unknown')
          if(idev.eq.4) call opencf(12,fint1)
          read(12,1013) hash,dum,nstedf,dum,ipl(2,1),dum,ipl(2,2),
     +                  dum,mkflg,dum,mmc,dum,mmst,dum,incmk,dum,pf
          read(12,1004) plotit
          read(12,1003) ((plo(ii,jj),ii=1,10),jj=1,nstedf)
          close(12)
          ipl(1,1) = ndf*(ipl(2,1)-1) + abs(ipl(2,2))! absolute dof to plot
          write(*,'(a)') 'Values for [tplo] are set to'
          write(*,'(a8,i5,a4,i3,a6,f9.4)')
     +        'node nr=',ipl(2,1),' dof=',ipl(2,2),' fact=',pf
          npldf  = 1               ! set flag to on
        end if
c...    open file (r...ldd) for restart
        fint1       = fres
        call addext(fint1,'ldd ')
        inquire(file=fint1,exist=lexst,opened=lopen)
        if(.not.lexst) then
          call drawmess('No LDD-file available on disk',1,0)
        end if
        if(lopen) return
        if(lexst .and. .not. lopen) then
          if(idev.lt.4) open(25,file=fint1,form='formatted',
     +     status='unknown')
          if(idev.eq.4) call opencf(25,fint1)
          read(25,'(10i6)') noden
          npldf  = 1                ! set flag to on
          npldf1 = 1                ! set flag to on
          do iii=1,nplo*50          ! 5 values(t,d,r,v,a)*10 nodes
            read(25,*,end=124)      ! goto end of file
          end do
124       backspace(25)
          call drawmess(
     +    'Nodes for  [tplo] are defined through restart',1,0)
          call drawmess('Do not use [tplo,def] again !!!!!',1,0)
cww       write(*,*) 'The following nodes are available'
cww       write(*,'(10i8)') noden
C
c......   set default values if no ldf-file exists
          if(nstedf.eq.0) then
            lct(l)  = 'set'         ! set default values for
            ct(1,l) = noden(1)      ! macro [tplo,init,n1,n2]
            ct(2,l) = 1
            ct(3,l) = 1
            goto 20
          end if
          return
        end if
        return
      end if

      if(pcomp(lct(l),'fact',4)) then ! do not use with ldd-file -> values in ldd will be changed
        n2  = ipl(2,2)
        n2b = abs(n2)
        n3  = 1
        n3  = isign(n3,n2)
        do i = 1,nplo
          plo(1,i) = plo(1,i) / pf
          plo(5,i) = plo(5,i) / pf
          plo(2,i) = plo(2,i) * n3
        end do
        n1 = ct(1,l)
        n3 = 1
        n3 = isign(n3,n1)
        ipl(2,2) = n3 * n2b
        pf = ct(2,l)
        if(pf.eq.0.0d0) pf=1.0d0
        do i = 1,nplo
          plo(1,i) = plo(1,i) * pf
          plo(5,i) = plo(5,i) * pf
          plo(2,i) = plo(2,i) * n3
        end do
        return
      end if

      if(pcomp(lct(l),'sset',4)) then ! set stre nstri at node nstrno
        nstri  = ct(1,l)
        nstrno = ct(2,l)
        nstri =min(max(1,nstri),25)
        nstrno=min(max(1,nstrno),numnp)
      end if

      if(pcomp(lct(l),'xval',4)) then
        xmint = ct(1,l)
        xmaxt = ct(2,l)
        imaxx = 1
        if(dsqrt(xmint*xmint+xmaxt*xmaxt).lt.1.e-8) imaxx = 0
      end if

      if(pcomp(lct(l),'yval',4)) then
        ymint = ct(1,l)
        ymaxt = ct(2,l)
        imaxy = 1
        if(dsqrt(ymint*ymint+ymaxt*ymaxt).lt.1.e-8) imaxy = 0
      end if

      if(pcomp(lct(l),'init',4)) then
        n1 = ct(1,l)
        n2 = ct(2,l)
        n2b= abs(n2)
        n1 = max(min(n1,numnp),1)
        n2b= max(min(n2b,ndf),1)
        pf = ct(3,l)
        if(pf.eq.0.0d0) pf=1.0d0
        npldf  = 1
c       nstedf = 1
        nstedf = 0
        ipl(1,1) = ndf*(n1-1)+n2b
        ipl(2,1) = n1
        ipl(2,2) = n2
        call pzero ( plo,10*nplo)
        imaxx = 0
        imaxy = 0
      end if

      if(pcomp(lct(l),'def',3)) then
        noden(1) = abs(ct(1,l))    ! input nodes as [tplo,def,n1,n2,n3]
        noden(2) = abs(ct(2,l))
        noden(3) = abs(ct(3,l))
99      if(noden(1).eq.0) then     ! input interactive
          write(*,*) 'Define (max 10) node-numbers: n1,n2,....,n10'
          call dinput(dnoden,10)
          do ii=1,10
            noden(ii)=dnoden(ii)
          end do
        end if
        do i=1,10                          ! check for input errors
          if(noden(i).lt.0 .or. noden(i).gt. numnp) then
            call drawmess('wrong node number',1,0)
            goto 99
          end if
        end do
        call drawmess(
     + 'Each [tplo,def] destroys all previous informations !!!',1,-2)
        npldf1 = 1                          ! set flag to on
        fint1 = fres                        ! set filename and open file
        call addext(fint1,'ldd ')
        if(idev.lt.4) open(25,file=fint1,form='formatted',
     +                                 status='unknown')
        if(idev.eq.4) call opencf(25,fint1)
        write(25,'(10i6)') noden
        write(25,'(7x,a)') '0 time, lambda, detc,stress,valuse1,valuse2'
        write(25,'(7x,a)')
     +  '(d)isplacement, (r)eactions, (v)elocities, (a)ccelerations'
        write(25,'(7x,6(4x,a2,5x))') 'v1','v2','v3','v4','v5','v6'
        lct(l)  = 'init'         ! set default values for
        ct(1,l) = noden(1)       ! macro [tplo,init,n1,n2]
        ct(2,l) = 1
        ct(3,l) = 0
        write(*,'(a)') 'Values for [tplo] are reset to'
        write(*,'(a8,i5,a4,i3,a6,f9.4)')
     +        'node nr=',noden(1),' dof=',1,' fact=',0.0d0
        goto 20
      end if

      if(pcomp(lct(l),'set',3)) then
        if(npldf1.ne.1) then
          call drawmess
     +      ('No nodes are defined, use macr [tplo,def] first',1,0)
          return
        end if
        if(int(abs(ct(1,l))).eq.0) then    ! check for error in node-input
          write(*,'(a)')  'The following nodes are available'
          write(*,'(10i8)') noden
          return
        end if
        if(int(abs(ct(2,l))).gt.ndf) then   ! check for error in dof-input
          call drawmess('DOF not defined',1,0)
          return
        end if
        do in=1,10
         if( int(abs(ct(1,l))).eq.noden(in) ) goto 127
        end do
        write(*,'(a)') 'A plot for this node is not available'
        write(*,'(a)') 'The following nodes are available'
        write(*,'(10i8)') noden
        return
127     continue
        n1 = abs(ct(1,l))        ! node number
        n2 = ct(2,l)             ! dof
        n2b= abs(n2)
        n2b= min(n2b,6)          ! only up to 6 dof cww
        pf = ct(3,l)             ! mult.
        ipl(1,1) = ndf*(n1-1)+n2b
        ipl(2,1) = n1
        ipl(2,2) = n2
        if(pf.eq.0.0d0) pf=1.0d0
        write(*,'(a)') 'Values for [tplo] are reset to'
        write(*,'(a8,i5,a4,i3,a6,f9.4)')
     +        'node nr=',n1,' dof=',n2,' fact=',pf
        rewind(25)                ! rewind file  (R  .ldd) file
        read(25,*)                ! dummy: noden
        read(25,*)                ! dummy: time
        read(25,*)                ! dummy: displacement
        read(25,*)                ! dummy: val1 val2 ...
        iz = 0
        ndfmax = min(ndf,6)              ! only up to 6 dofs possible
        do ii=1,nplo*50                  ! do until end of file 5 values(t,d,r,v,a)*10 nodes = 50
          read(25,4000,end=123) hash,inode,(uu(iuu),iuu=1,ndfmax) ! read from file
          if(hash.eq.'t' .and. inode.eq.0) then
            iz = iz + 1
            plo(3,iz)  = uu(1)    ! time
            plo(1,iz)  = uu(2)*pf ! lambda*fact
            plo(6,iz)  = uu(3)    ! detc
            plo(8,iz)  = uu(4)    ! stre
            plo(9,iz)  = uu(5)    ! valuse1
            plo(10,iz) = uu(6)    ! valuse2
          end if
          if(inode.eq.n1) then
            if(hash.eq.'d') then
              plo(2,iz) = uu(n2b)                    ! displacements
              if(ipl(2,2).lt.0 ) plo(2,iz) = -uu(n2b)! disp including sign
            end if
            if(hash.eq.'r')  plo(5,iz) = uu(n2b)*pf  ! reactions
            if(hash.eq.'v')  plo(4,iz) = uu(n2b)     ! velocities
            if(hash.eq.'a')  plo(7,iz) = uu(n2b)     ! accelerations
          end if
        end do
123     continue
        if(nstedf.eq.0) nstedf = iz
        do in=1,10
          if(n1.eq.noden(in)) react = reacc(in,n2b) ! for plotdf1
        end do
        backspace(25)
      end if

      if(pcomp(lct(l),'conv',4)) then
        nploc = ct(1,l)
        if(nploc.ne.0) nploc = 1
      end if

      if(pcomp(lct(l),'incr',4)) then
        nincp0 = ct(1,l)
      end if
c
4000  format(a1,i6,6(1x,e10.4))
c
cc    else
        if(lct(l).eq.'disp'.or.lct(l).eq.'phas'.or.lct(l).eq.'load'.or.
     1     lct(l).eq.'reac'.or.lct(l).eq.'det '.or.lct(l).eq.'velo'.or.
     2     lct(l).eq.'acce'.and.ior.ge.0) write(iow,3007)
cc    end if
      return
c
c.... eigenvalue computation using overrelaxation
21    continue
c...  [eigk,chek,n1,n2,n3]  eigenvalue computation during dynamic solution
c...                        n1 - max. number of cycles
c...                        n2 - over relaxation parameter
c...                        n3 - Solution tolerance
c...  [eigk,init,n1,n2,n3]  initialize eigenvalue computation
c...  [eigk,off]            switch computation off
      nzykel = ct(1,l)
      if(nzykel.eq.0) nzykel = 20
      omegax = ct(2,l)
      if (omegax.eq.0.d0) omegax = 1.4d0
      epseig = ct(3,l)
      if(epseig.eq.0.d0)  epseig = 1.d-5
      if(pcomp(lct(l),'chec',4)) then
        call ralloc(eigk1,neq,'EIGK-nxeig1',eigflg)
        if(.not.allocated(eigk2)) allocate(eigk2(neq))   ! EIGK-nxeig2
c       call initei(eigk1,neq)
        eigflg =tr
        return
      end if
      if(pcomp(lct(l),'off',3)) then
        eigflg =fa
        return
      end if
      call ralloc(eigk1,neq,'EIGK-nxeig1',eigflg)
      if(.not.allocated(eigk2)) allocate(eigk2(neq))   ! EIGK-nxeig2
      if(.not.pcomp(lct(l),'init',4)) then
        call initei(eigk1,neq)
      end if

      call cor(gstiff(nal),gstiff(nau),gstiff,eigk1,eigk2,massm
     1         ,jd,neq,nzykel,omegax,epseig)
      return
c
c.... eigencomputations for general real matrices
c     [ueig]     - conj. complex eigenvalues and - vectors
c.... isw = 21 - uses routine assemb in pform to assemble
c                      full non-symmetric matrix eigmma
c.... common block 'ueign' initiated in pmacr
22    if(fl(4)) then
        call drawmess('No stiffness matrix, use tang or utang',1,0)
        return
      end if

      call ralloc(aeigr,neq,'UEIG-meigr',ldummy)
      call ralloc(aeigi,neq,'UEIG-meigi',ldummy)
      call ralloc(aeigv,neq*neq,'UEIG-meigv',ldummy)


      sfl = .false.
      allocate(temp(neq+neq*neq))
      temp=0.d0
      mscal = 1
      mcnt  = mscal + neq
      allocate(eigmma(neq*neq))
      eigmma=0.d0
      if(ms.gt.maxm) then
        write(yyy,3006) ms,maxm
        call drawmess(yyy,1,0)
cww     if(ior.lt.0) return
cww     if(ior.gt.0) stop
        return
      end if
cww   hflgu/h3flgu for history terms?

      call formfe(u,dr,gstiff,gstiff(nal),tr,fl(8),cfr,fa,21,1,numel,1)

      call ueigen(eigmma,aeigr,aeigi,aeigv,temp,temp(mcnt),
     1            neq)
      deallocate(temp)
      deallocate(eigmma)
      return
c
c.... extended system
c...  [ext,,mext,eps]     extended system initialization
c...  [ext,off]           switch extended system off
c...  [ext,on,mext,eps]   switch extended system on, choose appr. for phi
c...  [ext,eps,exeps,kex] change eps-
23    kflg = tr
      if(pcomp(lct(l),'eps',3)) then
        exeps = ct(1,l)
        kex   = ct(2,l)
        return
      end if
      mext = ct(1,l)
      mext = max(min(mext,5),1)
      eps  = ct(2,l)
      iev  = ct(3,l)
      if(abs(eps) .lt. tol) eps = 1.d-6
      eflg = tr
      if(pcomp(lct(l),'off',3)) then
        extflg = fa
        kflg   = fa
        xmu    = 0.d0
        kex    = 0
        return
      end if
      if(pcomp(lct(l),'on',2)) then
        if(arcf)  arcf  = fa
        if(arcfs) arcfs = fa
        extflg = tr
c       mext = ct(1,l)
c       mext = max(min(mext,5),1)
c       eps  = ct(2,l)
c       if(abs(eps) .lt. tol) eps = 1.d-6
        call pzero(extkh,neq)
        if(mext.le.4) then

          call pmove2(gstiff(nal),gstiff(nau),gstiff,dr,extkh,jd,neq
     +                ,mext)

          iev=0
        else if(mext.eq.5) then
c....     use eigenvector n3
          if(.not.allocated(eigv).and.ior.lt.0) then
            call drawmess('Compute eigenvectors first',1,0)
            extflg = fa
            kflg   = fa
            xmu    = 0.d0
            kex    = 0
            return
          end if
          iev = max(1,min(iev,mf2*mf))
          jev = (iev-1)*neq
          call pmove(eigv(1+jev),extkh,neq)
          call norm(extkh,extkh,neq)
        end if
        if(ior.lt.0) write(*,2008) mext,iev
      end if
      if(extflg) return
      call ralloc(  extkh,   neq,'EXT-kh ',extflg)
      call ralloc( extkdh,  nneq,'EXT-kdh',extflg)
      call ralloc(  extkc,   neq,'EXT-kc ',extflg)
      call ralloc(  extkd,   neq,'EXT-kd ',extflg)
      call ralloc( extkz1,3*nneq,'EXT-kz1',extflg)
      call ralloc( extkz2,   neq,'EXT-kz2',extflg)
      call ralloc(  extke,   neq,'EXT-ke ',extflg)
      extflg = tr
      if(mext.le.4) then
        call pmove2(gstiff(nal),gstiff(nau),gstiff,dr,extkh,jd,neq,mext)
        iev=0
      else if(mext.eq.5) then
c....   use eigenvector n3
        if(.not.allocated(eigv).and.ior.lt.0) then
          call drawmess('Compute eigenvectors first',1,0)
          extflg = fa
          kflg   = fa
          xmu    = 0.d0
          kex    = 0
          return
        end if
        iev = max(1,min(iev,2*mf))
        jev = (iev-1)*neq
        call pmove(eigv(1+jev),extkh,neq)
        call norm(extkh,extkh,neq)
      end if
      if(ior.lt.0) write(*,2008) mext,iev
      return
c
c.... calculate buckling factor Lambda_i from omega_i
c...  [lamb,,n1]
c...             n1  >0 = number of eigen vector
c...             n1 <=0 = all eigen vectors
24    continue
      if(mfmax.eq.0.and.ior.lt.0) then
        call drawmess('Compute eigenvalues/vectors first',1,0)
        return
      end if
      if(istyp.eq.0) then
        gsize = jd(neq)+neq
      else
        gsize = jd(neq+1)-1
      end if
c.... temporary array for K_T
      allocate(glamb(gsize))
      glamb = 0.d0
      glamb = gstiff ! save stiffness matrix
      hflgu  = fa
      h3flgu = fa
c.... linear stiffness matrix (K_T with u = 0)
      call drawmess('Computing linear stiffness matrix',1,-2)
c.... temporary array for u (5 if ext.sys)
      allocate(temp(5*nneq))
      temp = 0.d0
      call pzero (dr,nneq)
      gstiff = 0.d0
      props  = prop  ! set all load dependent terms to zero, e.g. temperature  -> K_0!
      prop   = 0.d0
      call formfe(temp,dr,gstiff,gstiff(nal),tr,fa,fa,fa,3,1,numel,1) ! 
      prop = props
      deallocate(temp)
c...  calculate load factors
      n1 = ct(1,l)
      if(n1.gt.0) then       ! n1
        n1e = min(mf2*mf,n1)
        n1a = n1e
      else                   ! all 
        n1e = mf2*mf
        n1a = 1   
      end if      
       if(ior.lt.0) write(  *,2024) 
                    write(iow,2024) 
      do n1 = n1a,n1e
        kk  = (n1 - 1) *neq
c....   K_L*phi
        dr(1:neq) = 0.d0
       call promul(gstiff(nal),gstiff(nal),gstiff,eigv(1+kk),dr,jd,neq,
     +             1,1)
c....   phiT*K_L*phi
        cl  = ddot(neq,eigv(1+kk),1,dr,1)
c....   omega*phiT*phi
        omg = eigd(n1) 
        cu  = omg*ddot(neq,eigv(1+kk),1,eigv(1+kk),1)
c....   test
        if(dabs(cl-cu).lt.1.e-6) then
          call drawmess('K_T = K_L, Calculate K_T again!',1,0)
          return
        end if
c....   buckling load factor
        flmb= cl/(cl-cu)
        aflmb = flmb*prop
        if(arcf.or.extflg) aflmb = flmb * prop * rlnew
        if(ior.lt.0) write(  *,2025) n1,omg,flmb,aflmb
                     write(iow,2025) n1,omg,flmb,aflmb
      end do
      gstiff = glamb ! restore stiffness matrix
      deallocate(glamb)
      return
c
c.... [curv] input parameters for curves
25    call curvin(m)
      return
c
c.... remesh mesh
c     [reme,(adap,unif),n1,n2]
c     [reme,adap]    adaptive Netzverfeinerung
c     [reme,univ]    gleichf. Netzverfeinerung
c     [reme,new]     lese input file n-file
c     [reme,old]     lese input file neu
c     n1= percantage, n2=type of indicator: 1=energy,2=L2
26    continue
      call mreme (d,id,ix,u,dr,lct,ct,ndm,ndf,nen1,llreme)
      return
c
c.... calculate lowest eigenvalue evi and eigenvector eigi
C     using inverse iteration for (K_T-w1)phi=0
c     [eigi,,nstep,tol]
C     nstep = max. number of iter. steps ( default=50 )
C     tol   = used iter. tolerance       ( default=1.d-08 )
27    continue
      tols   = 1.d-12
      nstep  = 50
      if(ct(1,l).ne.0.0) nstep = ct(1,l)
      if(ct(2,l).ne.0.0) tols  = ct(2,l)
      call ralloc(eigia,neq,'EIGI',ldummy)

      call eigi (gstiff(nal),gstiff(nau),gstiff,jd,neq,eigia,dr,nstep
     +            ,tols,evi)

      return
c
c.... call to a system routine, only for dos/win
c     [dos]
c28    if(idev.eq.4) then
c        yyy = ' '
c        yyy = lct(l)
c        if(yyy.eq.' ') then
c          write(*,2015)
c          read(*,1012) yyy
c        end if
c        call start_pprocess(yyy,'')
c      end if


c.... RSG eigencomputations for 1 element
c     [eig1,n1]
c     only with on element
28    mf = nneq
      n1 = ct(1,l)
      if(n1.eq.0) n1 = 1

c.... set arrays
      dbgtxt = 'PMACR: [eig1] eig.vector,gen. array: mv,mf*mf*ipr'
      dbgtxt = 'PMACR: [eig1] eig.value, gen. array: md,mmf*ipr'
      call ralloc(eigv,mf*mf,'EIG1-EIGVEC',sfl)
      call ralloc(eigd,   mf,'EIG1-EIGVAL',sfl)
      mfmax=mf
      
c.... temporary used
      allocate(temp1(mf))
      allocate(temp2(mf))
      allocate(temp3(mf*mf))
      temp1 = 0.d0
      temp2 = 0.d0
      temp3 = 0.d0
c
c.... Get element matrix for elmt 1
      call formfe(dr,dummy,dummy,dummy,fa,fa,fa,fa,3,1,1,1)

c.... Solve eigenvalue problem
      call elemevp(s,temp3,eigd,eigv,temp1,temp2,nen,ndf,nst,n1)
      deallocate(temp1)
      deallocate(temp2)
      deallocate(temp3)
      return
c
c.... call user-defined editor
c     [edit]
29    if(idev.eq.3.or.idev.eq.4) call start_pprocess(editor,'')
      return
c
c.... write  output-file for postprocessing
c     [post,init]     initialize outputfile, write system data
c     [post,init,1]   ... renumbering without tied nodes
c     [post,disp]     write displacements
c     [post,reac]     write nodal reactions
c     [post,eigv,n1]  write eigenvector n1
c     [post,stre]     write nodal stresses
c     [post,warp,n1]  write nodal warping values (wt:n1=1,wq:n1=2)
c     [post,close]    close files
30    continue
c.... initialize
      if(pcomp(lct(l),'init',4)) then
        if(lindex) then
          call ialloc(apost,numnp,'POST-init',lindex)
        end if
        nmesh = ct(1,l)
        if(nmesh.eq.1) then
          call newmesh(x,ndm,numnp,apost,numnpn)
        else
          nunmpn = numnp
        end if
        call fepost(1,x,ix,ndm,ndf,nen1,n1,apost,numnpn)
      end if
c.... displacements
      if(pcomp(lct(l),'disp',4)) then
        call pmove(u,dr,nneq)
        call panglb(dr,bang,numnp,ndf) ! global directions
        call pmovedt(gtie,dr,numnp,ndf) ! tied nodes
        call fepost(2,dr,ix,ndm,ndf,nen1,n1,apost,numnpn)
      end if
c.... reactions
      if(pcomp(lct(l),'reac',4)) then
c....   be sure to calculate reactions directly before post,reac
          if(.not.rfl) then
          call drawmess('Calculate first nodal reactions',1,0)
          return
        end if
        call pmovedt(gtie,dr,numnp,ndf)
        call fepost(3,dr,ix,ndm,ndf,nen1,n1,apost,numnpn)
      end if
c.... eigenvectors
      if(pcomp(lct(l),'eigv',4)) then
        if(.not.allocated(eigv)) then
          call drawmess('Compute Eigenvector first ',1,0)
          return
        end if
        n1 = ct(1,l)
        if(n1.gt.mf2*mf) then
          call drawmess('No. of Eigenvector not possible ',1,0)
          return
        end if
        call pzero(dr,nneq)
        call pmovec(id,eigv(1+(n1-1)*neq),dr,numnp*ndf)
        call panglb(dr,bang,numnp,ndf) ! global directions
        call pmovedt(gtie,dr,numnp,ndf) ! tied nodes
        call fepost(4,dr,ix,ndm,ndf,nen1,n1,apost,numnpn)
      end if
c.... stresses
      if(pcomp(lct(l),'stre',4)) then
        n4 = 1 + numnp
        if(.not.fl(11)) then
          call drawmess('Calculate first nodal stresses',1,0)
          return
        end if
        call fepost(5,strea(n4),ix,ndm,ndf,nen1,n1,apost,numnpn)
      end if
c.... warping EL12: stre,4=wq, stre,5=wt
      if(pcomp(lct(l),'warp',4)) then
        if(.not.fl(11)) then
          call drawmess('Calculate first nodal stresses',1,0)
          return
        end if
        n1 = ct(1,l)
        if(n1.eq.0.or.n1.eq.1) then
          n6 = 1 + numnp*(4+1) !=pos5, wt
        else if(n1.eq.2) then
          n6 = 1 + numnp*(3+1) !=pos4, wq
        else
          call drawmess('No warping values available',1,0)
          return
        end if
        call fepost(6,strea(n6),ix,ndm,ndf,nen1,n1,apost,numnpn)
      end if
c.... close files
      if(pcomp(lct(l),'clos',4))
     +   call fepost(7,dr,ix,ndm,ndf,nen1,n1,apost,numnpn)
      return
c
c.... smooth mesh
c     [smoo,xxxx,n1,n2,n3]
c     xxxx='oesp': smoothing without edge-or boundary-nodes (default)
c     xxxx='wesp': smoothing with    edge-or boundary-nodes
c     n1: smoothing starts from node n1 (default=1).
c     n2: number of iterations          (default=1).
c     n3: weight factor for equality of sidelengths (first two numbers:w1)
c                       and 90-degree-angles  (third and fourth number:w2)
c                       (default=0101).
31    continue
      if (pcomp(lct(l),'    ',4)) lct(l) = 'oesp'
      call smoo_control
     +     (lct,ct,ndf,ndm,nen1,nst,nneq,kmax,prt,rlnew,
     +        timold,kflag,l,m)
      return
c
c.... write  output-file for postprocessing with tecplot
c     [tec,init,n2]    initialize rfile.tec, n1=1 n2=ityp define element
c     [tec,write]      write all values (def)n1=3
c     [tec,eigv,n1,n2] write n2(=scal)*eigenvector n1 instead of displ.
c     [tec,close]      close file            n1=2
32    continue
c.... initialize
      n2 = 1 ! dummy
      n4 = 1 ! dummy for init,clos
      if(pcomp(lct(l),'init',4)) then
        n1 = 1
        n2 = ct(1,l)
      end if
c.... write data
      if(pcomp(lct(l),'writ',4)) n1 = 3
      if(pcomp(lct(l),'eigv',4)) n1 = 3
c.... close file
      if(pcomp(lct(l),'clos',4)) n1 = 2
c
      if(n1.eq.3) then
        if(.not.fl(11)) then
          call drawmess('Calculate first nodal stresses',1,0)
          return
        end if
        n4 = 1 + numnp ! be sure to have nodal stresses, otherwise tecpost could not write data!
c...    modify displacements/eigenvector
        if(pcomp(lct(l),'writ',4)) then
          call pmove(u,dr,nneq)
        else if(pcomp(lct(l),'eigv',4)) then
          nn1 = ct(1,l)
          cn2 = ct(2,l)
          if (cn2.eq.0.d0) cn2=1.d0
          if(nn1.gt.mf2*mf) then
            call drawmess('No. of Eigenvector not possible ',1,0)
            return
          end if
          call pmovec(id,eigv(1+(nn1-1)*neq),dr,nneq)
          call pmoves(dr,dr,nneq,cn2)      ! multiply by factor
        end if
        call panglb(dr,bang,numnp,ndf) ! global directions
        call pmovedt(gtie,dr,numnp,ndf) ! tied nodes
      end if
c
      call tecpost(x,ix,strea(n4),dr,nen1,ndm,ndf,n1,n2)
      return
c
c.... 'parv' write output-files for postprocessing with paraview 4.01
c     [parv,init,n1,n2] initalize + first step, n1=el.typ, n2=VTU-No.
c                       file names: Rfile.pvd   xx=1-99,yyyy=1-9999
c     [parv,next,<n1>]  next time step, n1=t    Rfile.mxx.tyyyy.vtu
c     [parv,eigv,n1,n2] n1=el.typ,n2=EV_no      Rfile.EVxx.vtu
c.... initialize and write first time step
33    if(.not.fl(11)) then
        if(pcomp(lct(l),'eigv',4)) goto 331
        call drawmess('Calculate first nodal stresses',1,0)
        return
331     continue
      end if

      n1 = ct(1,l)
      n4 = 1 + numnp ! be sure to have nodal stresses, otherwise paraview could not write data!
      if(fl(9).and.flparv) then ! dynamic arrays
        call ralloc(parvvel,nneq,'Parv-Velocities',flparv)
        call ralloc(parvacce,nneq,'Parv-Accelerations',flparv)
      end if
      if(pcomp(lct(l),'init',4)) then
        n3 = 1
        n2 = ct(2,l)
        call pmove(u,dr,nneq)
      else if(pcomp(lct(l),'next',4)) then
c....   write next time step
        n3 = 2
        call pmove(u,dr,nneq)
      else if(pcomp(lct(l),'eigv',4)) then
c....   write Eigenvector n1 instead of displacements
        n3 = 3
        n2 = ct(2,l)
        if(n2.eq.0) n2 = 1
        if(n2.gt.mf2*mf.or.n2.gt.99) then
          call drawmess('No. of Eigenvector not possible ',1,0)
          return
        end if
        call pmovec(id,eigv(1+(n2-1)*neq),dr,nneq)
      else
        write(*,'(a26)') ' Command must be specified!'
        write(*,'(a)')   ' [parv,init,n1,n2]'
        write(*,'(a)')   ' [parv,next,n1]'
        write(*,'(a)')   ' [parv,eigv,n1,n2]'
        return
      end if
c.... modify displacements
      call panglb(dr,bang,numnp,ndf)! global directions
      call pmovedt(gtie,dr,numnp,ndf)! tied nodes
c.... modify stresses
      call pmovest(gtie,strea(n4),numnp,abs(istv))! tied nodes
c.... velocities and accelerations
      if(fl(9)) then
        call pmovec(id,trans,parvvel,nneq)
        call panglb(parvvel,bang,numnp,ndf)! global directions
        call pmovedt(gtie,parvvel,numnp,ndf)! tied nodes
        flgdyn2=fa
        if(nop.ne.2) then ! no accelarations for Backward Euler
                       nu = 1 +   nneq ! at n+1
          if(nop.eq.5) nu = 1 + 2*nneq ! at n+1/2
          call pmovec(id,trans(nu),parvacce,nneq)
          call panglb(parvacce,bang,numnp,ndf)! global directions
          call pmovedt(gtie,parvacce,numnp,ndf)! tied nodes
          flgdyn2=tr
        end if
      end if
c.... modify in case of polar coordinates, to be set in PLOT
      if(ipola.ne.0) then
        call dispola(x,dr,ndm,ndf,numnp,ipola)
        if(fl(9)) then
          call dispola(x,parvvel,ndm,ndf,numnp,ipola)
          call dispola(x,parvacce,ndm,ndf,numnp,ipola)
        end if
      end if
c.... write values
      if(ixtie.eq.0) then
        call pvpost(x,ix,strea(n4),dr,parvvel,parvacce,fl(9),flgdyn2,
     +              nen1,ndm,ndf,n1,n2,n3)
      else
        call pvpost(x,tecon,strea(n4),dr,parvvel,parvacce,fl(9),flgdyn2,
     +              nen1,ndm,ndf,n1,n2,n3)
      end if
c
      return
c
c.... set polar directions e.g. for disp
c.... [pola,i] i=0 --> cartesian, i = 12,13,23 --> polar in 1/2 1/3 2/3
34    ipola = ct(1,l)
      if(ipola.eq.0.or.ipola.eq.12.or.ipola.eq.13.or.ipola.eq.23) then
        write(*,2012) ipola
      else
        write(*,2011)
        ipola = 0
      end if
      return
c
c.... Lanczos eigencomputations
c     [lan,,n1,n2,n3]  - lanczos for n1 eigenpairs , n2 guard vectors
c     n1 no of eigenpairs
c     n2 guard vectors   
c     n3 = 0 (default)  --> tols = 1.d-12
c     n3 > 0            --> tols = tol
35    continue
c...  check data for EV-Problem 
      call evcheck(mb,imas,ieverror)
      if(ieverror.eq.1) return 
c...  set data
      allocate(eigm(size(massm)))
      eigm=massm
      tols = 1.d-12
      nits = 1
      if(ct(3,l).gt.0) tols = tol
      mf = ct(1,l)
      mad= ct(2,l)
      if(mad.le.0) mad = 50
      mf = min(neq,max(1,mf))
      mq = min(mad,neq)
      if(mf.gt.mad) mq=min(mf,neq)
      call numass(eigm,neq,mq,jd,imas)
      if(mq.lt.mf.and.ior.gt.0) write(iow,2001) mq
      if(mq.lt.mf.and.ior.lt.0) write(*  ,2001) mq
      if(mq.eq.0) then
        call drawmess('0 nonzero lumped mass terms, check CMAS',1,0)
        return
      end if

      sfl = pcomp(lct(l),'prin',4)
c.... the following defs depend on the value of mq and are only temporary used
      allocate(eigg (mq*mq))
      allocate(eigdp(mq))
      allocate(eigdt(mq))
      allocate(eigp (mq*mq))
      allocate(eigt1(neq))
      allocate(eigt2(neq))
      allocate(eigt3(neq))
      allocate(eigdh(mq))
      allocate(eiggh(mq))
      allocate(eigvh(mq*neq))
      allocate(eigv1(mq*neq))
      allocate(eigd1(mq    ))
      eigg  = 0.d0
      eigdp = 0.d0
      eigdt = 0.d0
      eigp  = 0.d0
      eigt1 = 0.d0
      eigt2 = 0.d0
      eigt3 = 0.d0
      eigdh = 0.d0
      eiggh = 0.d0
      eigvh = 0.d0
      eigv1 = 0.d0
      eigd1 = 0.d0

c.... Solve EV-problem
      write(iow,2026) mf,2*mf,mq
      write(*,  2026) mf,2*mf,mq
      call lanczos(gstiff,eigm,eigv1,eigd1,eigg,eigdp,eigdt,eigp,
     +             eigt1,eigt2,eigt3,eigdh,eiggh,eigvh,jd,mf,mq,
     +             neq,imas,shift,tols,sfl,nits)

c.... store 'save' evs      
      mfmax=mf 
      if(mf2.ne.1)mfmax=2*mf   
      call ralloc(eigv,neq*mfmax,'LAN-EIGVEC',sfl)
      call ralloc(eigd,    mfmax,'LAN-EIGVAL',sfl)
       
      call pmove(eigv1,eigv,mfmax*neq)
      call pmove(eigd1,eigd,mfmax)

c...  deallocating temporary arrays
      deallocate(eigm)
      deallocate(eigg)
      deallocate(eigdp)
      deallocate(eigdt)
      deallocate(eigp)
      deallocate(eigt1)
      deallocate(eigt2)
      deallocate(eigt3)
      deallocate(eigdh)
      deallocate(eigv1)
      deallocate(eigd1)

      return
c
c.... Calculate C=A^TKA-G^TK_aa^-1G and Tau (on micro level)
c     [SIGQ]
c     nss is defined via module epsd
36    if(flgfe2) then 
        call ralloc( mfe2a,nen*ndf*nss,'FE2-Mat AL' ,flgfe2)
        call ralloc(mfe2g1,    neq*nss,'FE2-Mat G1' ,flgfe2)
        call ralloc(mfe2g2,    neq*nss,'FE2-Mat G2' ,flgfe2)
        call ralloc(mfe2g3,nen*ndf*nss,'FE2-Mat GL' ,flgfe2)
        call ralloc(mfe2ix,nen*ndf*nss,'FE2-Mat IX' ,flgfe2)
        call ralloc(mfe2ht,    nss*nss,'FE2-Mat HT' ,flgfe2)
        call ralloc(mfe2tau,       nss,'FE2-Mat TAU',flgfe2)
      end if
      call fe2micro(mfe2a,mfe2g1,mfe2g2,mfe2g3,mfe2ix,
     +        mfe2ht,mfe2tau,link2,link3,x,u,dr,s,p,xl,id,ix,
     +        numnp,numel,nen,neq,flnk)

      return
c
c.... Set prescribed values for displacements on predefined boundaries
c     [EPSQ,,iepsd
c            iepsd=0 read eps from file, ..=1 read eps from input-macro
37    iepsd  = ct(1,l)
      call epsq_ma(x,f,u,id,link2,ndm,ndf,numnp,flnk,prt)
      return
c
c.... Plot nodal displacements along line
c     [Dplo,set]  set line
c     [Dplo]      reset line to zero
38    if(pcomp(lct(l),'set',3)) then
c...    set point A and B of line
        call pzero(dscor,6)
        call pzero(dsdcor,3)
        if(ior.lt.0) write(*,3000)
        call dinput(td,3)
        do i=1,ndm
          dscor(i,1)=td(i)
        end do
        if(ior.lt.0) write(*,3001)
        call dinput(td,3)
        do i=1,ndm
          dscor(i,2)=td(i)
        end do
        do i=1,ndm
          dsdcor(i)=dscor(i,2)-dscor(i,1)
        end do
        dsdcor2=dot(dsdcor,dsdcor,ndm)
        if(ior.lt.0) then
          write(*,2021)
          do i=1,2
            write(*,2022) (dscor(ii,i),ii=1,3)
          end do
        end if
        write(iow,2021)
        do i=1,2
          write(iow,2022) (dscor(ii,i),ii=1,3)
        end do
      else
        call pzero(dscor,6)
        dsdcor2=0.d0
      end if
      return
c
c.... FEAST eigencomputations
c     [feast,,n1,n2,n3]  - subspace for n1 eigenvalues ,
c                        -n1>=1.5 save values
c                        - n2=EV_min  - n3=EV_max
c                        -n2=n3=0 chosse values from Gershgorin
c     only INTEL and CSR-Storage
c
39    call evcheck(mb,imas,ieverror)
      if(ieverror.eq.1) return
c...  set data
      if(idev.ne.3) then
c        'solver only possible for idev=3(INTEL)'
         return
      end if
c      if(istyp.ge.3.and.istyp.le.8) then !S-LU/Pardiso/PBCG/PGMRES: CSR

c...  set data
      allocate(eigm(size(massm)))
      eigm=massm
      mf  = ct(1,l)
      mq=neq
      call numass(eigm,neq,mq,jd,imas)
      if(mq.lt.mf.and.ior.gt.0) write(iow,2001) mq
      if(mq.lt.mf.and.ior.lt.0) write(*  ,2001) mq
      if(mq.eq.0) then
        call drawmess('0 nonzero lumped mass terms, check CMAS',1,0)
        return
      end if

      mf  = min(mf,mq)

      call ralloc(eigv,neq*mf,'FEAST-EIGVEC',sfl)
      call ralloc(eigd,    mf,'FEAST-EIGVAL',sfl)
      mfmax=mf 

c.... set temporay data
      allocate(temp1(mf))
      temp1=0.d0
      ms  = kmax
      allocate(itemp1(1))
      allocate(itemp2(1))
      if(imas.eq.2) then
        deallocate(itemp1)
        deallocate(itemp2)
        allocate(itemp1(neq+1))
        allocate(itemp2(neq))
        itemp1=0
        itemp2=0
        mib = ms
        mjb = mib + neq+1
        ms  = mjb + neq
        call pconsi(itemp1,neq+1)
        call pconsi(itemp2,neq)
      end if

c.... Solve EV-problem
      call feast(gstiff,eigm,jd,csrja,itemp1,itemp2,eigd,eigv,
     +     temp1,neq,imas,ct(2,l),ct(3,l),mf,info,ior,iow,prt)

c.... deallocating temp arrays
      deallocate(eigm)
      deallocate(temp1)
      deallocate(itemp1)
      deallocate(itemp2)

      return
c
c.... error diagnostics
400   call drawmess('Attempt to change profile during mesh',1,0)
      return
c.... formats
1001  format(i5)
1002  format(a1)
1003  format(10(1x,e11.5))
1004  format(10a12)
1011  format(a1,'Steps ',i5,' Node ',i5,' DOF  ',i5,' Mark ',i5,
     +          ' Color',i5,' Type ',i5,' Inkr ',i5,' Fact ',e12.5)
1012  format(a80)
1013  format(a1,a6,i5,a6,i5,a6,i5,a6,i5,a6,i5,a6,i5,a6,i5,a6,e12.5)
2000  format('   Energy convergence test'/'    Maximum   =',1pe24.15/
     1       '    Current   =',1pe24.15/'    Tolerance =',1pe24.15)
2001  format(' Number eigenvalues reduced to',i4,' by number of nonzero
     1lumped mass terms')
2002  format(' Requested number of eigenvalues exceeds previous request'
     1      /'    Request =',i4,' Maximum =',i4/1x)
2003  format(' Scaling factor tau = ',1pe15.5,' using phi',i2)
2004  format
     +(' Bifurcation check: f*phi',i2,' = ',1pe15.5,' (bifur. point)')
2014  format
     +(' Bifurcation check: f*phi',i2,' = ',1pe15.5,' (limit point)')
2005  format(' Solution is diverging. Continue (y or n)? >',$)
2006  format(' Use eigenvector from extended system   = "1" ',/,
     1       ' use eigenvector from subspace analysis = "2" ',/,
     2       ' Input ==>',$)
2007  format(' dt has been set to 1.0 for arclength method ')
2008  format(' extended system started with mext/iev = ',i3,i3)
2009  format(' Imperfection',/,
     1       ' Eigen vector No.',i3,' Ampliude ',1pe15.5)
2011  format(' Displacement output set to cartesian directions ')
2012  format(' Displacement output set to polar (',i3,' ) directions ')
2015  format(' Define System Command > ',$)
2016  format(' parameter for step control are set as ',/
     +        ' max.iter [n1]=',i3,';  mono [n2]=',f5.2,';  exp=',f5.2,
     +        '  safety [n3]=',f5.2)
2017  format('Memory for LAMB too large(Needed:',i8,'Available:',i8,')')
2018  format(' Save Data on file ',a)
2019  format(' Save Data on file ',a,' Start with number ', i4)
2020  format(' Restart with Data from file ',a)
2021  format(' Coordinates for SPLO/DPLO')
2022  format(1x,3e14.7)
2023  format('[tplo]-values are saved for: ',
     +        'node nr ',i5,' dof ',i3,' fact ',f9.4)
2024  format(' Buckling load factors',/,
     +       '   No.    Omega           Lambda(omega)   Lambda*lambda') 
2025  format(1x,i5,3(1x,1pe15.5))     
2026  format(5x,'desired no of EV:',i5,'  calculated no of EV:'
     +       ,i5,'/',i5)

3000  format(' Input x,y,z of point 1 for SPLO/DPLO')
3001  format(' Input x,y,z of point 2 for SPLO/DPLO')
3002  format('Memory for SUBS/LAN/FEAST too large (Needed:',
     +       i8,' Available:',i8,' )')
3006  format(
     +'Memory n-sym.eigencom.too large (Needed:',i8,'Available:',i8,')')
3007  format(' **ERROR** No load deflection plots in batch mode.')
      end
c
      subroutine pmacr4(ul,xl,tl,ld,p,s,ie,d,id,x,ix,f0,f,t,jp,u,dr
     1   ,lct,ct,ndf,ndm,nen1,nst,nneq,prt,wd,j,jsh,plo)
c-----------------------------------------------------------------------
c.... macro instruction subprogram 4                                   |
c
c     'mac1','mac2','mac3','mac4','mac5','elem','man ','para',
c
c-----------------------------------------------------------------------
      USE codat
      USE iofile
      USE pdata2
      implicit real*8 (a-h,o-z)
      logical prt,pcomp
      logical fa,tr
      character*4 lct,wd(1)
      integer ld(*),ie(*),id(*),ix(*),jp(*)
      real*8  ul(*),xl(*),tl(*),p(*),s(*),ct(3),
     1      d(*),x(*),f0(*),f(*),t(*),u(*),dr(*),plo(10,*)
      data fa,tr/.false.,.true./
c.... transfer to correct process
      ne = 0
      go to (1,2,3,4,5,6,7,8), j
c
c.... 'mac1' command
c     [mac1,,n1,n2,n3]
c     [mac1,name]  - rename 'mac1' to 'name'
1     if(lct.ne.'    ') then
        if(ior.gt.0) write(iow,2000) wd(j+jsh),lct
        if(ior.lt.0) write(  *,2000) wd(j+jsh),lct
        wd(j+jsh) = lct
      else
        call umacr1(ul,xl,tl,ld,p,s,ie,d,id,x,ix,f0,f,t,jp,u,dr,
     +              lct,ct,ndf,ndm,nen1,nst,nneq,ne,prt,plo)
      end if
      return
c
c.... 'mac2' command
c     [mac2,,n1,n2,n3]
c     [mac2,name]  - rename 'mac2' to 'name'
2     if(lct.ne.'    ') then
        if(ior.gt.0) write(iow,2000) wd(j+jsh),lct
        if(ior.lt.0) write(  *,2000) wd(j+jsh),lct
        wd(j+jsh) = lct
      else
        call umacr2(ul,xl,tl,ld,p,s,ie,d,id,x,ix,f0,f,t,jp,u,dr,
     +              lct,ct,ndf,ndm,nen1,nst,nneq,ne,prt,plo)
      end if
      return
c
c.... 'mac3' command
c     [mac3,,n1,n2,n3]
c     [mac3,name]  - rename 'mac3' to 'name'
3     if(lct.ne.'    ') then
        if(ior.gt.0) write(iow,2000) wd(j+jsh),lct
        if(ior.lt.0) write(  *,2000) wd(j+jsh),lct
        wd(j+jsh) = lct
      else
        call umacr3(ul,xl,tl,ld,p,s,ie,d,id,x,ix,f0,f,t,jp,u,dr,
     +              lct,ct,ndf,ndm,nen1,nst,nneq,ne,prt,plo)
      end if
      return
c
c.... 'mac4' command
c     [mac4,,n1,n2,n3]
c     [mac4,name]  - rename 'mac4' to 'name'
4     if(lct.ne.'    ') then
        if(ior.gt.0) write(iow,2000) wd(j+jsh),lct
        if(ior.lt.0) write(  *,2000) wd(j+jsh),lct
        wd(j+jsh) = lct
      else
        call umacr4(ul,xl,tl,ld,p,s,ie,d,id,x,ix,f0,f,t,jp,u,dr,
     +              lct,ct,ndf,ndm,nen1,nst,nneq,ne,prt,plo)
      end if
      return
c
c.... 'mac5' command
c     [mac5,,n1,n2,n3]
c     [mac5,name]  - rename 'mac5' to 'name'
5     if(lct.ne.'    ') then
        if(ior.gt.0) write(iow,2000) wd(j+jsh),lct
        if(ior.lt.0) write(  *,2000) wd(j+jsh),lct
        wd(j+jsh) = lct
      else
        call umacr5(ul,xl,tl,ld,p,s,ie,d,id,x,ix,f0,f,t,jp,u,dr,
     +              lct,ct,ndf,ndm,nen1,nst,nneq,ne,prt,plo)
      end if
      return
c
c     [elem]  list all available elements in feap
cww6    call helpel
6     return
c
c     [man ]  get the complete help manual
7     call pman(4,'titl')
      return
c
c.... [PARA] set parameter variables
8     if(pcomp(lct,'    ',4)) then
        coflg = tr
        call pconst(prt)
      else
        call pconstm(lct,prt)
        end if
      return
c
c.... format
2000  format(3x,'*WARNING* MACRO NAME ',a4,' CHANGED TO ',a4,'**')
      end
c
      subroutine pmacr5 (ul,xl,tl,ld,p,s,ie,d,id,x,ix,f0,f,t,jd,u,dr,
     +                  lct,ct,ndf,ndm,nen1,nst,nneq,prt,wd,j,jsh)
c-----------------------------------------------------------------------
c.... macro instruction subprogram 5                                   |
c
c     'zcop','bsys',
c
c-----------------------------------------------------------------------
      USE cdata
      implicit double precision (a-h,o-y)
c --- transfer to correct process
      go to (1,2), j
c
c.... [zcopy,i] i=1 for real, i=2 for complex
1     write(*,3000)
      return
c.... [BSYS] write COORdinates and ELEMents to file Bname-binary
2     call writerb(x,ix,ndm,numnp,nen,nen1,numel)
      return
c --- formats
3000  format(//5x,'Not Included in this FEAP Version')
      end
c
