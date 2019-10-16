      subroutine pplotf(ul,xl,tl,ld,p,s,ie,d,id,x,ix,f,f0,t,jp,b,dr,
     1                  lci,ct,prop,ndf,ndm,nen1,nst,plo)
c-----------------------------------------------------------------------
c.... Purpose: macro instruction subprogram
c              Controls plot by order of specifying macro commands
c              in array wd.
c
c.... Inputs:
c     nn    ul(nst,6)      - element displacements
c     n0    xl(ndm,nen)    - element coordinates
c     n1    tl(nen)        - element temperatures
c     n2    ld(nst)        -
c     n3    p(nst)         - element load vector
c     n4    s(nst,nst)     - element stiffness matrix
c     n5    ie(nie,numat)  - assembly information for material set nie=ndf+2
c     n6    d(ndd,numat)   - material set parameters
c     n7    id(ndf,numnp)  - equation numbers for each active dof
c     n8    x(ndm,numnp)   - nodal coordinates of mesh
c     n9    ix(nen1,numel) - element nodal connections of mesh
c     n10   f(ndf,numnp)   - load vector
c     n13   f0(ndf,numnp)  - nodal initial force values
c     n11   t(numnp)       - temperature vector
c     n12   jp(*)=jd       - pointer array for row/columns of tangent
c     n14   b(ndf,3*numnp) - displacement vector
c           dr(ndf,numnp)  - working vector
c           lci=lct        - actual plot macro
c           ct(3)          - data for actual plot macro
c           prop           - load factor
c           ndf            - number dof/node
c           ndm            - spatial dimension of mesh
c           nen            - max. number of nodes/element
c           nen1           - dimension for ix array: nen+4
c           nst            - dimension for element array: ndf*nen
c           plo(10,nplo)   - array for plotting TPLO-data
c
c      Outputs:            - desired plot arrays
c
c      Open:
c         me,kk  =10?
c 
c     update ww KIT 12/14
c-----------------------------------------------------------------------
      USE adap  
      USE arcl
      USE cdat1
      USE cdata
      USE comfil
      USE contval
      USE damp1
      USE ddata
      USE debugs
      USE dirdat
      USE eig1
      USE errin1
      USE errin2
      USE errin3
      USE ext1
      USE ext2
      USE fdata
      USE hdatam
      USE hidden
      USE hptext1
      USE iodata
      USE iofile
      USE iwinio
      USE lplot1
      USE mdata
      USE mdat2
      USE ndata
      USE ndatx
      USE nolink
      USE pback
      USE pcent
      USE pdam
      USE pdata1
      USE pdata2
      USE pdata3
      USE pdata4
      USE pdata7
      USE pdata8
      USE pdata10
      USE pdata11
      USE pdata12
      USE pdatah
      USE pdatap
      USE pdatas
      USE pftn77
      USE plodf
      USE plodfb
      USE plong
      USE plotdrp
      USE plotter
      USE plslay
      USE pltran
      USE pnodn
      USE ppers
      USE Psize
      USE qload
      USE rndata
      USE rsum
      USE strnam
      USE subdt
      USE tdata
      USE uneig
      USE ximp
      USE yydata
      USE doalloc
      implicit double precision (a-h,o-z)
      logical pcomp,outf,errv,hidm,hids
      logical tr,fa,clip,fdpl,defo
      logical isNurbs
      character*1  ctval4(15)
      character*2  extens
      character*4   lci,lct,lcc,wd,wdummy
      real*4 tsl,red,green,blue,tary, tary1
cww   integer*2 iyyy
      integer*2 icolo
      dimension ul(*),xl(*),tl(*),ld(*),p(*),s(*),ie(*),d(*),id(ndf,*),
     1     x(ndm,*),ix(*),f(ndf,*),t(*),jp(*),b(ndf,*),dr(*),ct(3),
     2     wd(102),ss(4),tt(4),f0(ndf,*),am(8),dra(3),plo(10,*)

c.... allocatable arrays
      integer, allocatable, dimension(:) :: mea,mza,mnix
      integer, allocatable, dimension(:) :: itm, indplo
      real*8,  allocatable, dimension(:) :: rtm, rndplo

c.... ***** add end quit in wd array as dummy ****
      common /gener/ kk,km,ketemp
cww   common /pdata9/ ifrm
cww   common /pswit/  imod
      common /rmsh1/  nxdrm,melem   !### delete when adap ok    

      common m(maxm)
         
      save   l,nsizt,iln,nxd,nxn,nne,nsym,clip,outf,cs,
     2     fdpl,msym,hidm,defo,nnope,cinv,hids,nnemax
c
c.... plot data command list
      data wd/'fram','wipe','fact','isom','cart','line','defm','disp',
     1        'outl','load','mesh','stre','node','boun','elem','zoom',
     2        'colo','mate','text','size','eigi','eigv','bord','scal',
     3        'axis','pers','adis','show','hide','ueig','movi','symm',
     4        'ceig','prin','tplo','flux','pris','forc','matn','logo',
     5        'reac','clip','rot0','rot1','rot2','rot3','evex','move',
     6        'titl','plof','hmsh','mono','dmag','pola','erro','back',
     7        'dplo','splo','rplo','eplo','xtic','hids','defo','tie ',
     8        'evan','base','init','xsca','rotm','velo','acce','ndii',
     9        'link','avel','aacc','aeig','aeve','draw','man ','rmsh',
     +        'slee','prof','cent','sect','copy','maxi','resi','angl',
     1        'rsum','pnod','pele','str1','jint','isec','traj','end',
     2        'quit','magn','pdis','ints','dimp','fill'/
      data list/102/
      data tr /.true./, fa /.false./
      data ss/.25d0,.75d0,.25d0,.75d0/,tt/.75d0,.75d0,.25d0,.25d0/
      data am /-1.d0,-.66d0,.0d0, .66d0,
     1          1.d0, .66d0,.0d0,-.66d0/


c.... macros which are in WD as dummys
c     TITL works
c     XTIC works

      if(pfl) go to 550
c.... definitions, same as [init]
c....   pseta
        pfl     = tr
        outf    = tr
        fdrp    = tr
        fdpl    = tr
        hidm    = tr
        hidi    = tr
        hidw1   = tr
        hidw2   = tr
        iln     = 1
c....   cart
        iso   = fa
c....   perspective
        kpers = 0
        call pzero (eold,3)
        call pzero (vold,3)
        vold(3) = 1.0d0
c....   rot
        call pzero(tra,9)
        do i = 1,3
          tra(i,i) = 1.0
          vr(i)    = 0.0
        end do
        call pzero(rotang,3)
c....   elements hide
        nxd   = nen1
        nxn   = nen
        nne   = numel
        nnemax  = 4*numel
        hide  = fa
        hids  = fa
c....   deformation
        defo  = fa
        cs    = 1.0
c....   zoom
        nzm1 =  0
        nzm2 =  0
        nzm3 =  0
        call pzero(xzm,6)
c....   clip
        clip = fa
        fact = 1.d0
c....   move
        deltax  = 0.d0
        deltay  = 0.d0
c....   xscal coor
        xfac(1) = 1.0d0
        xfac(2) = 1.0d0
        xfac(3) = 1.0d0
c....   sym
        nsym = 0
        msym = 0
c....   plot pola
        ipola = 0
c....   colors etc
        icgm    = 0
        imono   = 0
        iback   = 0
c...    Contour values
        icv     = 0
        call pzero (contv,3)
c....   maxi stress
        nmn     = 0
        nmx     = 0
cwwc....   frame
cww        ifrm    = 0
c....   plot/prin options
        ipgl    = 1
        nnope   = 0
        iclear  = 0
        iopl    = iop
        nexte   = 0
        lhpgl   = fa
        lps     = fa
c....   hptext
        isize   = 1
c....   set frame
        call frame(x,ndm,numnp,1,clip)
c...  end definitions
      if(ior.lt.0) then
c....   for interactive plot devices
        if(idev.eq.1) then
c.....    on PHIGS screen
c         call gpcrss(5,1,1,' ')
          call gpassw(1,5)
          if(icgm.eq.1) call gpassw(2,5)
        else if(idev.eq.2) then
c.....    on GKS screen
        else if(idev.eq.3) then
c.....    on   screen
        else if(idev.eq.4) then
c.....    on WINDOWS screen
        end if
      else
c....   for batch files
        ncplt = ipos(fplt,229)
        write(iow,2003) fplt(1:ncplt)
      end if

550   continue

c.... no plot for batch
      if(iplot.eq.0) return

      propq = prop ! set actual load
      if(rlnew.gt.0.d0) propq = prop*rlnew

c.... mixed array for el-no relation elements/faces m(nix)=ix or m(mi)
      mnixsize=max(5*nnemax,nen1*numel)
      if(.not. allocated(mnix)) allocate(mnix(mnixsize))
      if(hide) then
        mnix(1:5*nnemax)=mia(1:5*nnemax)
      else
        mnix(1:nen1*numel)=ix(1:nen1*numel)
      end if
c
c...  define array for hidden line algorithm
      if(hidw1) then 
        dbgtxt = 'PPLOTF: hidden line(IDIS)mh1,nnemax*1'
        call ialloc(idis,nnemax,'IDIS',hidw1)
      end if
      call plthid3(idis,nnemax) 
c
      cinv   = 1
      nsizt  = 1
      ijump  = 0
      ct4    = 0
      ipld   = 0
      icolo  = 7
      ibps   = 0
      ifor   = 13
      fopn   = fa
      call plopen
      call plclos
260   if(ior.lt.0.and.pcomp(lci,'    ',4)) then ! wait on plot macro
        ijump  = 1
        if(idev.lt.3) then
          if(icgm.eq.0) then
            if(nexte.eq.0) write(*,2004)   ! ijump plot
            if(nexte.ne.0) write(*,2014)   ! ijump prin
          else if(icgm.ne.0) then
            write(*,2013)                  ! ijump plof
          end if
        else if(idev.ge.3) then
          if(nexte.eq.0) write(*,2004)  ! ijump plot
          if(nexte.ne.0) write(*,2014)  ! ijump prin
        end if
        call pintio (yyy,15)
        read(yyy,1001,err=265) lct,lcc
        if(pcomp(lct,'end ',4).or.pcomp(lct,'quit',4).or.
     1     pcomp(lct,'q   ',1)) then
          call plclos
c...      close  hpgl or ps file, do not close use w.read/write
c         if(lhpgl) call hpglfile (5,dummy,dummy,dummy,dummy)
c         if(lps)   call psfile   (5,dummy,dummy,dummy,dummy,wdummy)
          return
        end if
      else   ! plot macro set under macro
        lct = lci
      end if
      if(ior.lt.0.and.pcomp(lct,'help',4)) then
        jflag = 0
        do 610 jj = 1,list
          if(pcomp(lcc,wd(jj),4)) then
          call pman(3,wd(jj))
          jflag = 1
          end if
610     continue
        if(jflag.eq.0) call phelpmp(wd,list,'PLOT',lcc)
        if(ijump.ne.0) go to 260
        return
      end if
      if(ijump.ne.0) then ! when direct from plot
c....   get macro instructions, ct4 only from plot possible
c        read(yyy,1000,err=265) lct,ct(1),ct(2),ct(3),ct4
        read(yyy,1000,err=265) lct,ct(1),ct(2),ct(3),ctval4
        call setval(ctval4(1:15),15,ct4)
      end if
      ipb = 0
      do 300 l = 1,list
        if(pcomp(lct,wd(l),4)) go to 301
300   continue
      write(yyy,'(a,a,a)') 'Macro  ',lct(1:4),'  does not exist'
      call drawmess(yyy,icolo,1)
      if(ijump.ne.0) go to 260
      return
c.... write macro on ps file
301   if(nexte.gt.0 .and. ihpgl.eq.1)
     +  call psfile(9,dummy,dummy,dummy,dummy,wd(l))
c
c----------------------------------------------------------------------
c            f  w  f  i  c  l  d  d  o  l  m  s  n  b  e  z  c  m  t  s |
c            r  i  a  s  a  i  e  i  u  o  e  t  o  o  l  o  o  a  e  i |
c            a  p  c  o  r  n  f  s  t  a  s  r  d  u  e  o  l  t  x  z |
c            m  e  t  m  t  e  m  p  l  d  h  e  e  n  m  m  o  e  t  e |
      go to( 1, 2, 3, 4, 4, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,
c----------------------------------------------------------------------
c            e  e  b  s  a  p  a  s  h  u  m  s  c  p  t  f  p  f  m  l |
c            i  i  o  c  x  e  d  h  i  e  o  y  e  r  p  l  r  o  a  o |
c            g  g  r  a  i  r  i  o  d  i  v  m  i  i  l  u  i  r  t  g |
c            i  v  d  l  s  s  s  w  e  g  i  m  g  n  o  x  s  c  n  o |
     1      21,22,23,24,25,26,10,28,29,30, 7,32,33,34,35,36,37,38,39,40,
c----------------------------------------------------------------------
c            r  c  r  r  r  r  e  m  t  p  h  m  d  p  e  b  d  s  r  e |
c            e  l  o  o  o  o  v  o  i  l  m  o  m  o  r  a  p  p  p  p |
c            a  i  t  t  t  t  e  v  t  o  s  n  a  l  r  c  l  l  l  l |
c            c  p  0  1  2  3  x  e  l  f  h  o  g  a  o  k  o  o  o  o |
     2      41,42,43,43,43,43,47,48,49,50,51,52,53,54,55,56,57,12,41,57,
c----------------------------------------------------------------------
c            x  h  d  t  e  d  i  x  r  v  a  n  l  a  a  a  a  d  m  r |
c            t  i  e  i  v  r  n  s  o  e  c  d  i  v  a  e  e  r  a  e |
c            i  d  f  e  a  e  i  c  t  l  c  i  n  e  c  i  v  a  n  m |
c            c  s  o     n  c  t  a  m  o  e  i  k  l  c  g  e  w     e |
     3      61,62,63,64,65,10,67,68,69,70,70,72,73,10,10,10,10,78,79,80,
c----------------------------------------------------------------------
c            s  p  c  s  c  m  r  a  r  p  p  s  j  i  t  e  q  m  p    |
c            l  r  e  e  o  a  e  n  s  n  e  t  i  s  r  n  u  a  d    |
c            e  o  n  c  p  x  s  g  u  o  l  r  n  e  a  d  i  g  i    |
c            e  f  t  t  y  i  i  l  m  d  e  1  t  c  j     t  n  s    |
     4      81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,10,
c----------------------------------------------------------------------
c            i   d   f |
c            n   i   i |
c            t   m   l |
c            s   p   l |
     5     100,101,102 ),l
c----------------------------------------------------------------------
c
c.... new frame not active, only reset scale
c.... [fram,v1] - defines plotting frame (v1- is the quadrant 1--4)
1     continue
cww1     ifrm = ct(1)
cww      ifrm = min(4,ifrm)
cww110   if(ifrm.ge.1) then
cww        scale = 0.5*scaleg*fact
cww        s0(1) = ss(ifrm)
cww        s0(2) = tt(ifrm)
cww      else
cww      end if

110   scale = scaleg*fact
      s0(1) = 0.5
      s0(2) = 0.5
      iclear = 0
      call plopen
      goto 200
c
c.... clean screen
c     [wipe]
c.... [wipe,v1]  v1: wipe between bounds def. by mouse cursor
c                v2: wipe outside bounds def. by mouse cursor
2     k1 = ct(1)
      if(k1.ne.0) then
        call pltwip(k1)
      else
        scale = scaleg*fact
        s0(1) = 0.5
        s0(2) = 0.5
        iclear = 0
        call plopen
      end if
      go to 200
c.... change plot scale
c.... [fact,v1] - changes plot scale (v1- factor)
3     if(ct(1).eq.0.0d0) then
        fact = 1.0
      else
        fact = ct(1)
      end if
      iclear = 0
      go to 110
c.... isometric view (l = 4) or cartesian view (l = 5, default value)
c.... [isom] - isometric view
c.... [cart] - cartesian view (default)
4     iso   = l.eq.4
      mnix(1:nen1*numel)=ix(1:nen1*numel)
      nxd   = nen1
      nxn   = nen
      nne   = numel
      hide  = fa
      kpers = 0
      call frame(x,ndm,numnp,1,clip)
      iclear = 0
      go to 200
c 5   [cart] under 4
c.... line type routine port
c.... [line,v1] - set line type (v1- defines type of line plotting)
6     iln = ct(1)
      if(iln.eq.0) iln = 1
      call plopen
      call plline(iln)
      go to 200
c.... plot deformed mesh
c.... [defm,v1,v2,v3] - plot deformed mesh (v1 - scaling of deformation)
c                                       (v2 - plot for material no. v2)
c                                       (v3 - color number,default = 5(green))
7     c = ct(1)
      k1 = ct(3)
      if(k1.le.0) k1 = 5
      if(c.eq.0.0d0) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      if(outf) dbgtxt = 'PPLOTF: [defm] gen. array: me,10*numnp*1(outf)'
      outf=fa
      if(.not.allocated(mea)) allocate(mea(10*numnp))
      call plopen
      call pppcol(k1)
      call pline(dr,ie,mnix,mea,numnp,nne,ndm,nxd,nxn,nie,ct(2),tr)
      if(lct.eq.'movi') then
        call pppcol(32)
        call pline(dr,ie,mnix,mea,numnp,nne,ndm,nxd,nxn,nie,ct(2),tr)
      end if
      go to 150
c.... displacement contours
c.... [disp,v1,v2,v3,v4]   plot filled /contour displacements
c                    v1 -    number of displ., neg. on deform. mesh
c                    v2 >  0 number of displ. values for plot
c                    v2 <= 0 filled plot
c                    v3.ne.0 contour plot without numbers
c                    v4 -  plot mesh in color v4
8     continue
      call etimef(tary)
      if(iadd.gt.0)goto 810
      k2 = ct(2)
      n = abs(k2)
      n = max(1,n)
      i = max(1,abs(int(ct(1))))
      i = min(i,ndf)
      ipb = ct(3)
      k5  = ct4
      if(fdrp)dbgtxt='PPLOTF:[disp] gen.array: numnp*ndf*ipr(fdrp)'
      call ralloc(drp,numnp*ndf,'pplotf-help-field',fdrp)
      call pmove(b,drp,numnp*ndf)
      call panglb(drp,bang,numnp,ndf)
      c = 0.d0
      if(ct(1).lt.0.0d0.or.defo) c = cs
810   continue
      if(icleas.eq.0) then
        if(ipola.ne.0) call dispola(x,drp,ndm,ndf,numnp,ipola)
        call rprint1(drp,ix,x,ndm,numnp,ndf,n,nen1,i,idev)
      end if
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      if(k2.le.0) then
        n = -i
      end if
      call pltcon(dr,id,ie,mnix,drp,idis,nie,ndm,ndf,nxd,nxn,nne,
     1            i,n,2,1,k5,cinv)
      icleas = 1
      call etimef(tary1)
      if(pfr)                write(iow,2018) tary1-tary  ! time for DISP
      if(pfr .and. ior.lt.0) write(*  ,2018) tary1-tary
      go to 150
c.... plot outline of parts
c.... [outl,v1,v2,v3] - outline of mesh(v1 < 0 on deformed mesh)
c                                      (v2 : outline for material v2)
c                                      (v3 : color for all)
9     if(outf) dbgtxt = 'PPLOTF: [outl] gen. array: me,10*numnp*1(outf)'
      outf=fa
      if(.not.allocated(mea)) allocate(mea(10*numnp))
      c = 0.0
      if(ct(1).lt.0.0d0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      call plopen
      k3 = ct(3)
      if(k3.eq.0) k3=1
      call pppcol(k3)
      call pline(dr,ie,mnix,mea,numnp,nne,ndm,nxd,nxn,nie,ct(2),fa)
      go to 150
c.... plot of loads,displ.,veloc.,acce.,director,eigenvectors, presc. d.
c                       (v1 < 0 on deformed mesh)
c                       (v2 .ne.0: vector tip to node)
c                       (v2  =  0: vector tip away from node (default)
c                       (v3  scaling factor(default=1)
c.... [load,v1,v2,v3] 10   plot loads (separately!)
c                       (v1 =  i:      force(i) on undeformed mesh)
c                       (v1 = -i:      force(i) on   deformed mesh)
c                       (v1 =  0:      all forces (1-ndm) on undeformed mesh)
c                       (v1 = (ndf+1): all forces (1-ndm) on undeformed mesh)
c                       (v1 =-(ndf+1): all forces (1-ndm) on   deformed mesh)
c                       (v1 = (ndf+2): all moments(4-5/6) on undeformed mesh)
c                       (v1 =-(ndf+2): all forces (4-5/6) on   deformed mesh)
c
c     [base,v1,v2,v3]   66   director
c     [adis,v1,v2,v3]   27   displacements
c     [avel,v1,v2,v3]   74   velocities
c     [aacc,v1,v2,v3]   75   accelerations
c.... [aeig,v1,v2,v3]   76   plot eigenvector v1
c.... [aevx,v1,v2,v3]   77   plot eigenvector of extended system
c.... [pdis,v1,v2,v3]   99   prescribed displacements
c
10    call plopen
      c   = 0.d0
      k1  = abs(int(ct(1)))
      k1  = min(k1,ndf+2)
      if(ct(1).lt.0.0d0.or.defo) c = cs
      k2  = ct(2) ! tip
      rk3 = dabs(ct(3))
      k4  = ct(3)
      if(ct(3).eq.0) rk3 = 1.0d0
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      if(fdrp)dbgtxt='PPLOTF:[load]gen. array: numnp*ndf*ipr(fdrp)'
      call ralloc(drp,numnp*ndf,'pplotf-help-field',fdrp)
      if(l.eq.10) then
c....   loads
        if(k1.le.3)     call pppcol(2) ! force
        if(k1.gt.3)     call pppcol(6) ! moment
        if(k1.eq.ndf+1) call pppcol(2) ! forces
        call ploads(b,drp,propq,tr,fa,fa)
        call setfor(f,f0,propq,numnp*ndf,drp )
        if(k4.lt.0) then
c.....    plot only loaded nodes
          call pltload(dr,drp,id,ndm,ndf,numnp,1)
        else
c.....    plot arrows(separate for each dof) + nodes
          call pltfor2(dr,drp,bang,id,ndm,ndf,numnp,k2,rk3,
     +                 tra,vr,1,k1,ix,nen1,1)
        end if
      else if(l.eq.99) then
c....   prescribed displacements
        call pppcol(5)
        call pmove(f,drp,numnp*ndf)
        if(k4.lt.0) then
c.....    plot only nodes with prescribed displacements
          call pltload(dr,drp,id,ndm,ndf,numnp,3)
        else
          call pltfor(dr,drp,bang,id,ndm,ndf,numnp,k2,rk3,tra,vr,
     +              2,ix,nen1,6,2)
        end if
      else if(l.eq.27) then
c....   displacements
        call pppcol(5)
        call pltfor(dr,b,bang,id,ndm,ndf,numnp,k2,rk3,tra,vr,2,
     +              ix,nen1,1,1)
      else if(l.eq.74) then
c....   velocities
        if(.not.fl(9)) then
          call drawmess('Problem not dynamic ',icolo,1)
          goto 200
        end if
        call pmovec(id,trans,drp,numnp*ndf)
        call panglb(drp,bang,numnp,ndf)
        call pppcol(4)
        call pltfor(dr,drp,bang,id,ndm,ndf,numnp,k2,rk3,tra,vr,
     +              2,ix,nen1,2,1)
      else if(l.eq.75) then
c....   accelerations
        if(.not.fl(9)) then
          call drawmess('Problem not dynamic ',icolo,1)
          goto 200
        end if
        nu = 1 + numnp*ndf
        call pmovec(id,trans(nu),drp,numnp*ndf)
        call panglb(drp,bang,numnp,ndf)
        call pppcol(2)
        call pltfor(dr,drp,bang,id,ndm,ndf,numnp,k2,rk3,tra,vr,
     +              2,ix,nen1,3,1)
      else if(l.eq.66) then
c....   director
        xm = dabs(ct(1))
        if(xm.eq.0.d0) xm = 1.d0
        k3 = ct(3)
        if(ldir.eq.1) call pltfor1(dr,ix,basea,xdir,knode,numnp,numel,
     +                             nen,nen1,ndm,xm,k2,k3,tra,vr,ipgl)
      else if(l.eq.76) then
c....   eigenvector
        k1 = min(mfmax,max(1,k1))
        if(mfmax.eq.0) then
          call drawmess('Compute Eigenvector first ',icolo,1)
          goto 200
        end if
        call pmovec(id,eigv(1+(k1-1)*neq),drp,numnp*ndf)
        call panglb(drp,bang,numnp,ndf)
        call pppcol(8)
        call pltfor(dr,drp,bang,id,ndm,ndf,numnp,k2,rk3,tra,vr,
     +              2,ix,nen1,3,1)
      else if(l.eq.77) then
c....   eigenvector from extended system
        if(.not.extflg) then
          call drawmess('Compute EV of ext. system first ',icolo,1)
          goto 200
        end if
        call pmovec(id,extkh,drp,numnp*ndf)
        call panglb(drp,bang,numnp,ndf)
        call pppcol(8)
        call pltfor(dr,drp,bang,id,ndm,ndf,numnp,k2,rk3,tra,vr,
     +              2,ix,nen1,4,1)
      end if
      go to 150
c.... plot mesh
c.... [mesh,v1,v2,v3] -(v1 < 0  on deformed mesh)
c                      (v2 - plot for material no. v2)
c                      (v3 - color number,default = 4(cyan))
11    dbgtxt = 'PPLOTF: [mesh] gen. array: me,10*numnp*1(outf)'
      outf=fa
      if(.not.allocated(mea)) allocate(mea(10*numnp))
      c = 0.0
      k1 = ct(3)
      if(k1.le.0) k1 = 4
      if(ct(1).lt.0.0d0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      call plopen
      call pppcol(k1)
      if(nne.eq.0) call drawmess('No ELEMENTS available',icolo,1)
      call pline(dr,ie,mnix,mea,numnp,nne,ndm,nxd,nxn,nie,ct(2),tr)
      go to 150
c.... stress contours
c.... [stre,v1,v2,v3,v4]    plot filled /contour stresses
c                    v1 -   number of stress, neg. on deform. mesh
c                    v2 > 0 contour plot
c                    v2 <=0 filled plot
c                    v3 < 0 filled plot for v3=k.m for layer k pos m
c                    v3 > 0 contour plot without numbers
c                    v4 -   plot mesh in color v4
c                    v1 = 27 --> main stresses, v3 = fact/length, v4 = fact/tip,
c                                               v2 = typ:2=2D/3=3D
c
c     #1 shell elements
c     stre,i,        plot stress resultant i
c     stre,i,,-k.m   plot stress i at layer k and gp m   (m=0 -> m=1)
c
c     #2 solid and solid shell elements (one element in thickness dir.)
c     without layer
c     stre,i,        plot stress i
c     with layer
c     stre,i,        plot sum of stress i through all layers,
c                    e.g. sum of damage through thickness
c     stre,i,,-k.m   plot stress i at layer k and gp m (m=0 -> m=1)
c
c     #3 solid and solid shell elements (more element in thickness dir.)
c     without layer
c     stre,i,        plot stress i
c     with layer
c     makes no sense
12    ipb = ct(3)
      if(plfl) then
        dbgtxt = 'PPLOTF: [stre] gen. array: np,numnp*npstr*ipr(plfl)'
        call ralloc(strea,numnp*npstr,'STRE',plfl)
      end if
      k1 = abs(ct(1))
      k1 = min(npstr,max(k1,1))
      k2 = ct(2)
      rk2= ct(2)
      if(rk2.eq.0) rk2 = 1.0d0
      k5 = ct4
      rk4= ct4
      if(rk4.eq.0) rk4 = 1.0d0
      klay = 0
      mlay = 1
      if(ipb.lt.0) then
        play = dabs(ct(3))
        klay = int(play)  ! number   of layer
        play = (play-klay)*10
        mlay = nint(play) ! position in layer
        if(mlay.eq.0) mlay=1
        ipb  = 0
      end if
      k3 = 1 + numnp
cww   if(.not.fl(11).or.(klay.ne.0.and.icleas.eq.0)) then
      if(.not.fl(11).or.(icleas.eq.0)) then
        ener1 = 0
        ener2 = 0
        strea = 0.d0
        hflgu  = fa
        h3flgu = fa
cww     call etimef(tary)
        call formfe(b,dr,dr,dr,fa,fa,fa,fa,8,1,numel,1)
        call pltstr(strea,strea(k3),numnp)
cww     call etimef(tary1)
cww     if(pfr .and. ior.lt.0) write(*  ,2016) tary1-tary
      end if
      if(k1.eq.npstr) then
        c = 0.0d0
        if(ct(1).lt.0.0d0.or.defo) c = cs
        rk3=ct(3)
        if(rk3.eq.0) rk3 = 1.0d0
        call plopen
        call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
        call pltmain(dr,strea(k3),ndm,numnp,tra,vr,k2,rk3,rk4,ipgl)
        go to 121
      end if
      k4 = k3 + numnp*(k1-1)
c.... for splo
      if(l.eq.58) goto 57
c     if(k1.eq.8) k4 = np
      if(icleas.eq.0) then
cwd   for isogeometric version call rprintIGA
        if (isNurbs()) then
            call rprintIGA(k1,ix,x,ndm,numnp,1,nen1,k2,idev,numel)
        else
            call rprint(strea(k4),ix,x,ndm,numnp,1,nen1,k2,idev)
        end if
      end if
      c = 0.0
      if(ct(1).lt.0.0d0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
cww   call frame(dr,ndm,numnp,1,clip)
      call plxtrn(dr,tra,vr,ndm,numnp)
      if(k2.le.0) then
        k2 = -k1
        k1 = 1
      end if
      call etimef(tary)
      call pltcon(dr,id,ie,mnix,strea(k4),idis,nie,ndm,1,nxd,nxn,nne,
     1            k1,k2,1,1,k5,cinv)
      call etimef(tary1)
      if(pfr .and. ior.lt.0) write(*  ,2017) tary1-tary ! time for STRE
      icleas = 1
121   fl(11) = tr
      goto 150
c.... plot/label nodes
c.... [node,v1,v2,v3] plot nodes
c                  (v1.lt.0 plot nodes on deform. mesh)
c                  (v2.gt.0 plot nodes numbers)
c                  (v2.lt.0 plot only node with number v2)
c                  (v3.ne.0 plot in color v3 (default 7 (yellow))
13    c = 0.0
      if(ct(1).lt.0.0d0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      k2 = ct(2)
      k3 = ct(3)
      if(k3.le.0) k3 = 7
      call plopen
      call pltnod(dr,ix,nen1,ndm,numnp,k2,k3)
      go to 150
c.... plot/label boundary restraints
c.... [boun,v1,v2] plot boundary constraints
c          abs(v1) = scaling factor (default = 1)
c              v1 < 0 plot bound. cond. on deform. mesh
c              v2 = dof to plot,   default = 0 (all)
14    c = 0.0
      if(ct(1).lt.0.0d0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      k2   = ct(2)
      k2   = min(max(0,k2),ndf)
      rk4  = 1.
      rrk4 = dabs(ct(1))
      if(rrk4.gt.0.0) rk4 = rrk4
      call plopen
      call pltbou(id,dr,bang,ix,nen1,ndm,ndf,numnp,k2,rk4,tra,vr)
      go to 150
c.... label elements with numbers
c.... [elem,v1,v2,v3] plot element numbers
c                  (v1 <  0 plot  on deform. mesh)
c                  (v2 >  0 plot  element v2)
c                  (v2 <  0 plot  element v2 on deform. mesh)
c                  (v3.ne.0 plot in color v3 (default 3 (blue))
15    c = 0.0
      if(ct(1).lt.0.0d0.or.ct(2).lt.0.0d0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      k2 = abs(ct(2))
      k3 = ct(3)
      if(k3.le.0) k3 = 3
      call plopen
      call pppcol(k3)
      if(nne.eq.0) call drawmess('No ELEMENTS available',icolo,1)
      call pltelm(dr,ie,mnix,idis,scale,nie,ndm,nxn,nxd,nne,k2)
      go to 150
c.... set window for zoom
c.... [zoom,v1,v2,v3] zoom system between nodes v1 and v2
c          v3.ne.0 zoom from (x,y,z)_a to (x,y,z)_b
16    nzm1 = ct(1)
      nzm2 = ct(2)
      nzm3 = ct(3)
      if(nzm3.gt.0) then
        write(*,2015)
        call dinput(xzm,6)
      end if
      iclear = 0
      clip = fa
160   c = 0.d0
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
c.... construct perspective projection if necessary
      if(kpers.ne.0) then
        call perspj(dr,dr,ndm,ndm,numnp,errv)
      end if
      isy = 1
      if(nsym.gt.0) isy = 3
      call frame(dr,ndm,numnp,isy,clip)
      go to 150
c.... set pencolor for plot (only works on some of the plots)
c.... [colo,v1,v2] set color to v1
c              v2 = +1 set color range red->blue(def.), v2=-1 blue->red
17    icol = ct(1)
      if(icol.eq.0) then
cww     iback=0
        imono=0
        iclear = 0
        call plopen
      end if
      call pppcol(icol)
      cinv = ct(2)
      if(cinv.eq.0) cinv = 1
      go to 200
c.... plot material of mesh filled
c.... [mate,v1,v2,v3] plot materials
c           v1: material number,neg. plot on def.mesh 
c           v2: color number for plot 
c           v3: >0 plot mesh on top 
c            
18    c = 0.0
      if(ct(1).lt.0.0d0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      k1 = abs(ct(1))
      k2 = ct(2)
      if(k2.lt.0) k2 = 1
      k3 = 0
cww   iclear = 0
      call plopen
      call pltmate ! Legend
      call pppcol(k2)
      call plot2d(ie,mnix,idis,b,dr,xl,nie,ndm,nxn,nxd,nne,k1,k2,k3)
      k4 = ct(3)
      if(k4.gt.0) then ! plot mesh
        ct(2)= k1    
        ct(3)= 32 
        call plclos
        goto 11
      end if
      go to 150
c.... put text on screen
c     [text,v1,v2,v3] put text on screen (v1: color)
c        not activ in plttxt     (v2: size  of(1-4) default = 1)
c                                        (v3: plot line v3 in r. window)
19    k1 = ct(1)
      k2 = ct(2)
      v3 = ct(3)
      if(k1.ne.0) call pppcol(k1)
      call plopen
ccc not activ      nsizt = min(4,max(1,k2))
      call plttxt(v3,nsizt)
      go to 200
c.... set text size
c     [size,v1] set size for all text including labels,legend
20    nsizt = ct(1)
      if(nsizt.eq.0) nsizt = 1
c     nsizt = min(4,max(1,nsizt))
      call pltsiz(nsizt)
      if(nexte.gt.0.and.ihpgl.eq.1) then ! scale text in postscript-file
        x1=nsizt
        call psfile(8,x1,dummy,dummy,dummy,wdummy)
      end if
      go to 200
c
c.... plot eigenvector from inverse iteration
c.... [eigi,v1,v2,v3,v4 (v1),v2,v3,v4 see eigv(22)
21    k1 = -1
      k11= ct(1)
      k5 = ct4
      if(.not.allocated(eigia)) then
        call drawmess('Compute EV from inv. it. first ',icolo,1)
        goto 200
      end if
      call pmovec(id,eigi,dr,numnp*ndf)
      goto 221
c
c.... plot eigenvector
c.... [eigv,v1,v2,v3,v4]
c                     v1 - no. of eigen vector
c                     v1 < 0  filled and hiddenline mesh for  v3 = 0
c                     v1 < 0  degree v3 on deformed(=ev) mesh v3 > 0
c           mesh: v3 = 0  eig. vect. as deformed mesh)
ccc                       v2      scaling factor (default from scale)
c                     v2 >=0  hidden line from 1 to numel
c                     v2 < 0  hidden line from  numel to 1
c           dof : v3 > 0  degree of freedom for filled plot
c                         v2 <=0  filled plot
c           dof : v3 < 0  plot contour lines without numbers
c                         v2 > 0 number of cont. lines for plot)
c                         v4 - plot mesh in color v4
22    k11 = ct(1)
      k1 = abs(ct(1))
      k1 = min(mfmax,max(1,k1))
      k5 = ct4
      if(mfmax.eq.0) then
        call drawmess('Compute Eigenvector first ',icolo,1)
        goto 200
      end if
      call pmovec(id,eigv(1+(k1-1)*neq),dr, numnp*ndf)
221   i = abs(ct(3))
      if(ct(3).lt.0) ipb = 1
      if(fdrp)dbgtxt='PPLOTF:[eigv]gen. array: numnp*ndf*ipr(fdrp)'
      call ralloc(drp,numnp*ndf,'pplotf-help-field',fdrp)
      call pmove(dr,drp,numnp*ndf)
      call psymm(drp,ndf,numnp,isym2,ndm,iadd+1,its)
c     c  = ct(2)
c     if(c.eq.0.0d0) c = cs
c.... scale
      if(i.eq.0) then ! full EV
        c = cs
      else
        if(k11.lt.0) then ! Single dof of EV
          c = cs
        else
          c = 0.d0
        end if
      end if
      call pdefm(x,drp,c,bang,ndm,ndf,numnp,drp,xfac)
      call plxtrn(drp,tra,vr,ndm,numnp)
      if(i.eq.0) then
        if(outf) dbgtxt='PPLOTF: [eigv] gen. array: me,10*numnp*1(outf)'
        outf=fa
        if(.not.allocated(mea)) allocate(mea(10*numnp))
        call plopen
        if(k11.ge.0) then
          call pppcol(8)
          call pline(drp,ie,mnix,mea,numnp,nne,ndm,nxd,nxn,
     1               nie,0.0d0,tr)
c    1               nie,ct(2),tr)
        else if(k11.lt.0) then
          k2 = ct(2)
          k22 = 1
          if(k2.lt.0) k22 = 2
          call plot2dh(ie,mnix,drp,xl,nie,ndm,nxn,nxd,nne,
     +                 32,27,k22,idis)
        end if
        if(l.eq.21) call pleigvt(evi,k1)
        if(l.eq.22) call pleigvt(eigd(k1),k1)
        if(l.eq.47) call pleigvt(propq,k1)
      else
        i  = min(i,ndf)
        k2 = ct(2)
        n  = abs(k2)
        call panglb(dr,bang,numnp,ndf)
        if(icleas.eq.0) then
          if(ipola.ne.0) call dispola(x,dr,ndm,ndf,numnp,ipola)
            call rprint1(dr,ix,x,ndm,numnp,ndf,n,nen1,i,idev)
        end if
        if(k2.le.0) then
          n = -i
        end if
        call pltcon(drp,id,ie,mnix,dr,idis,nie,ndm,ndf,nxd,nxn,nne, 
     1              i,n,5,k1,k5,cinv)
        icleas = 1
      end if
      go to 150
c.... toggle border
c.... [bord,v1] plot border in color v1 (default: 1 (white))
c....       v1 = 0 ->wipe + original state
c....         v1 < 0 do not plot border, logo, pers, rot
23    ibcol = ct(1)
      ibcol = max(1,min(ibcol,32))
      ibor  = ct(1)
      call plopen
      if(ibor.le.0) then
        iclear = 0
      else
        call plbord(ibcol)
      end if
      go to 200
c.... rescale for deformed window
c.... [scal,v1] scal for deformed mesh by factor cs (default: 1)
c               v2=..., v3=...
24    nzm1 = 0
      nzm2 = 0
      nzm3 = 0
      call pzero(xzm,6)
      jj   = 2
      cs   = ct(1)
      if(cs.eq.0.0d0) cs = 1.0
      call pdefm(x,b,cs,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
c     if(ct(2).eq.0.0d0.and.ndf.ge.ndm) call pwind(x,dr,ndm,ndf,numnp)
      if(ct(3).eq.0.0d0) jj = 1
c.... construct perspective projection if necessary
      if(kpers.ne.0) then
        call perspj(dr,dr,ndm,ndm,numnp,errv)
      end if
      if(nsym.gt.0) jj = 3
      call frame(dr,ndm,numnp,jj,clip)
      iclear = 0
      go to 150
c.... add axes to the plot
c.... [axis,v1,v2,v3] plot a set of axes at point v1,v2,v3 in color 2(red)
25    xm = 0.1*max(xmax(1)-xmin(1),xmax(2)-xmin(2),xmax(3)-xmin(3))
      call plopen
      call pppcol(2)
      call pltaxs(tra,vr,ct,ndm,xm)
      go to 200
c.... [pers,v1] input parameters for perspective projection
c           v1 = 1  include hide,1
c           v1 = 2  include hids,1
26    k1 = ct(1)
      if(ndm.lt.3) then
        kpers = 0
        call drawmess('PERSP only for 3-D poss. ',icolo,1)
        go to 200
      end if
      call maxcor(x,ndm,numnp)
      call perspe
      call perspj(x,dr,ndm,ndm,numnp,errv)
      iclear = 0
      call plopen
      if(errv) then
        kpers = 0
        go to 200
      else
        kpers = 1
      end if
      mnix(1:nen1*numel)=ix(1:nen1*numel)
      nxd   = nen1
      nxn   = nen
      nne   = numel
      call frame(dr,ndm,numnp,1,clip)
      if(k1.eq.1) goto 29
      if(k1.eq.2) then
        ct(1) = 1
        goto 62
      end if
      go to 200
c 27  [disp,v1,v2,v3] under 10
c.... [show] display current plot parameters
28    call plshowp(idev,ior,fplt,1,iln,nsizt,cs,fact,xfac,rotang,
     1     deltax,deltay,kpers,ipola,msym,aipma,nummat,iso,defo,iopl)
      go to 200
c
c.... [hide,v1] plot always visible faces of a 3d- mesh!
c....  hidden line have to be added to see results really correct
c.... v1 = 0 stop Hide
c.... v1 < 0 deformed mesh
c.... v2 .ne.0 show the hidden lines in color 8(magenta),not active
c.....v4 = plot color,not used
29    if(ct(1).eq.0) then
        hide = fa
        mnix(1:nen1*numel)=ix(1:nen1*numel)
        nxd   = nen1
        nxn   = nen
        nne   = numel
        goto 200
      end if
      if(kpers.eq.0) then
cww     call drawmess('Hide only for PERSP poss. ',icolo,1)
cww     goto 200
        call maxcor(x,ndm,numnp)
        call plthid4
      end if
      if(outf) dbgtxt = 'PPLOTF: [hide] gen. array: me,10*numnp*1(outf)'
      outf=fa
      if(.not.allocated(mea)) allocate(mea(10*numnp))
      if(hidm) dbgtxt = 'PPLOTF: [hidm]gen. array: mz,numel*1(hidm)'
      hidm = fa
      if(.not.allocated(mza)) allocate(mza(numel))
      hide = tr
      c    = 0.0
      if(ct(1).lt.0.0d0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
c.... only if mesh is plotted
cww      iclear = 0
cww      call plopen
cww      icol = 4
cww      if(ct(1).lt.0.0d0) icol=5
cww      call pppcol(icol)

c      izcol = ct4
c     if(izcol.gt.0) call pppcol(izcol)

c.... compute location of faces which face forward
      call plface(mza,mea,ix,dr,ndm,nen1,numnp,numel,
     1            iln,ct(2),nface,icol)

      if(nface.gt.nnemax) then !  nnemax=fixed maximum length
        call drawmess(
     1  'nface < nnemax not fulfilled, Plot not possible',icolo,1)
        go to 200
      end if

      if(hidi) then 
        dbgtxt = 'PPLOTF: [hidi]gen. array: mi,5*nnemax*1(hidi)'
        call ialloc(mia,5*nnemax,'HIDI',hidi)
      end if
      call plfacx(mza,mea,ix,mia,dr,ndm,nen1,numnp,numel)
c.... set plot mesh quantities (all values are reset by cart or hide,0!)
      if(.not.allocated(mnix))allocate(mnix(5*nnemax))
      mnix(1:5*nnemax)=mia(1:5*nnemax)
      nxd = 5     !nen1
      nxn = 4     !nen 
      nne = nface !numel 
      ct(1) = 0
      goto 62 ! hids
c.... plot unsymmetric eigenvector
c.... [ueig,v1,v2,v3,v4] needs to be installed similar to eigv!!
30    c = ct(2)
      if(c.eq.0.0d0.or.defo) c = cs
      k1 = ct(1)
      if(k1.eq.0) k1 = 1
      if(k1.gt.0) call pmovec(id,aeigv(1+(k1-1)*neq),dr,numnp*ndf)
      i = abs(ct(3))
      k5= ct4
      if(i.eq.0) then
        call psymm(dr,ndf,numnp,isym2,ndm,iadd+1,its)
        call pdefm(x,dr,c,bang,ndm,ndf,numnp,dr,xfac)
        call plxtrn(dr,tra,vr,ndm,numnp)
        if(outf)dbgtxt='PPLOTF: [ueig] gen. array: me,10*numnp*1(outf)'
        outf=fa
        if(.not.allocated(mea)) allocate(mea(10*numnp))
        call plopen
        call pppcol(8)
        call pline(dr,ie,mnix,mea,numnp,nne,ndm,nxd,nxn,nie,ct(2),tr)
      else
        i  = min(i,ndf)
        k2 = ct(2)
        n  = abs(k2)
        if(icleas.eq.0) then
cww       ... introduce dispola
          call rprint(dr(i),ix,x,ndm,numnp,ndf,nen1,n,idev)
        end if
        if(k2.le.0) then
          n = -i
        end if
        call pltcon(x,id,ie,mnix,dr,idis,nie,ndm,ndf,nxd,nxn,nne,
     1              i,n,5,k1,k5,cinv)
        icleas = 1
      end if
      go to 150
c 31  [movi] under 7
c.... set up parameters for adding quadrants in case of symmetry
c.... [symm,v1,v2] v1 = 0 - reset to no symmetry
c                  v1 = 1 - symmetry with respect to x-axis
c                  v1 = 2 - symmetry with respect to y-axis
c                  v1 = 3 - symmetry with respect to x/y-axis
c                  v1 = 4 - symmetry with respect to z in-xy-plane(Q1-4)
c                  v1 = 5 - symmetry with respect to z in-xy-plane(Q1)
c .................v2 ne.0 -> do not rescale and clear the window
32    call pzeroi(isym1,72)
      k1 = ct(1)
      k2 = ct(2)
      msym = k1
      if(k1.eq.0) then
c....   no symmetry
        nsym  = 0
        its   = 2
      else if(k1.eq.1) then
c....   symmetry with respect to x
        isym1(1,2,1) = 1
        isym1(2,2,1) = 1
        nsym  = 2
        its   = 1
      else if(k1.eq.2) then
c....   symmetry with respect to y
        isym1(1,1,1) = 1
        isym1(2,1,1) = 1
        nsym  = 2
        its   = 2
      else if(k1.eq.3) then
c....   symmetry with respect to x,y
        isym1(1,1,1) = 1
        isym1(2,2,1) = 1
        isym1(3,1,1) = 1
        isym1(4,2,1) = 1
        nsym  = 4
        its   = 2
      else if(k1.eq.4) then
c....   symmetry with respect to z in xy plane (with symm in Q1-4)
        if(ndm.gt.2) then
          isym1(1,1,1) = 1
          isym1(2,2,1) = 1
          isym1(3,1,1) = 1
          isym1(4,2,1) = 1
          isym1(4,3,1) = 1
          isym1(5,1,1) = 1
          isym1(6,2,1) = 1
          isym1(7,1,1) = 1
          isym1(8,2,1) = 1
          isym1(8,3,1) = 1
          nsym  = 8
          its   = 2
        end if
      else if(k1.eq.5) then
c....   symmetry with respect to z in xy plane (only symm to xy)
        if(ndm.gt.2) then
          isym1(1,3,1) = 1
          isym1(2,3,1) = 1
          nsym  = 2
          its   = 3
        end if
      end if
      if(k2.eq.0) then
c....   scale +x1/+x2 quadrant
        c = cs
        call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
        call plxtrn(dr,tra,vr,ndm,numnp)
        if(kpers.ne.0) call perspj(dr,dr,ndm,ndm,numnp,errv)
        call frame(dr,ndm,numnp,1,clip)
c....   scale +-x/+-y/+-z quadrant
        if(nsym.ne.0) then
          do 153 is = 1,nsym
            call psymm(dr,ndm,numnp,isym1,ndm,is,1)
            call plxtrn(dr,tra,vr,ndm,numnp)
            if(kpers.ne.0) call perspj(dr,dr,ndm,ndm,numnp,errv)
            call frame(dr,ndm,numnp,3,clip)
153       continue
        end if
        iclear = 0
      end if
      go to 200
c.... plot eigenvector
c.... [ceig,v1,v2,v3,v4] needs to be installed similar to eigv!!
33    c = ct(1)
      if(c.eq.0.0d0.or.defo) c = cs
      call pmovec(id,eigk1,dr, numnp*ndf)
      i  = abs(ct(3))
      k5 = ct4
      if(i.eq.0) then
        call psymm(dr,ndf,numnp,isym2,ndm,iadd+1,its)
        call pdefm(x,dr,c,bang,ndm,ndf,numnp,dr,xfac)
        call plxtrn(dr,tra,vr,ndm,numnp)
        if(outf) dbgtxt='PPLOTF: [ceig] gen. array: me,10*numnp*1(outf)'
        outf=fa
        if(.not.allocated(mea)) allocate(mea(10*numnp))
        call plopen
        call pppcol(8)
        call pline(dr,ie,mnix,mea,numnp,nne,ndm,nxd,nxn,nie,ct(1),tr)
      else
        i  = min(i,ndf)
        k2 = ct(1)
        n  = abs(k2)
        if(icleas.eq.0) then
cww       dispola !!
          call rprint(dr(i),ix,x,ndm,numnp,ndf,nen1,n,idev)
        end if
        if(k2.le.0) then
          n = -i
        end if
        call pltcon(x,id,ie,mnix,dr,idis,nie,ndm,ndf,nxd,nxn,nne,
     1              i,n,5,1,k5,cinv)
        icleas = 1
      end if
      go to 150
c.... [prin,v1,v2,v3]
c      plot  data on  HPGL/PS-file via own procedures
c      v1 =    abs of number of plotfile 1 - 99
c      v1 > 0  in color  for psfile
c      v1 < 0  in grey   for psfile
c      v1 = 0  close plotfile
c      v1 > 99 show existing plotfiles (only WIN)
c      v2 = 1  (default) plot on file fplt_v1.eps   for PS
c      v2 = 2            plot on file fplt_v1.pgl   for HPGL
c      v3 = 1  A4 portrait  (default) only for PS
c      v3 = 2  A4 landscape           only for PS
34    ipfile = ct(1)
      nexte  = abs(ipfile)
      imono  = 0
      if(ipfile.lt.0) imono = 1
      ihpgl  = ct(2)
      ihpgl  = max(1,min(ihpgl,2))
      if(ihpgl.eq.2) ipgl = 3      ! set no fill option for hpgl
      ilsc   = ct(3)
      ilsc   = max(1,min(ilsc,2))
c...  test for existing eps and hpgl files
      if(nexte.gt.99.and.idev.eq.4) then
        fpgl = fplt
        call prin_ps
        nexte = 0
      end if
c...  open file
      if(nexte.gt.0 .and.  nexte.le.99 .and. nnope.eq.0) then
        fpgl = fplt
        write(extens,'(i2)') nexte
        if(nexte.le.9) extens(1:1)='0'
        ncol = ipos(fpgl,229)
        if(ncol.gt.(229-7))
     +    write(*,*) 'Plot-filename too long, more than 122 characters!'
        fpgl(ncol+1:ncol+1) = '_'
        fpgl(ncol+2:ncol+2) = extens(1:1)
        fpgl(ncol+3:ncol+3) = extens(2:2)
        if(ihpgl.eq.1) then
          call addext(fpgl,'eps ')
          call psfile(1,31.d0,dummy,dummy,dummy,wdummy)
          if(.not.lps) nexte=0 ! in case of not overwrite existing file
          if(lps) call psfile(6, 0.d0,dummy,dummy,dummy,wdummy)
        else if(ihpgl.eq.2) then
          call addext(fpgl,'pgl ')
          call hpglfile(1,30.d0,0.d0,0.d0,0.d0)
        end if
        nnope  = 1
cww     iclear = 0    c for using splo etc!
      end if
c...  close file
      if(nexte.eq.0) then
        ipgl  = 1
        ihpgl = 1
        imono = 0
        nnope = 0
        if(lps)   call psfile   (5,dummy,dummy,dummy,dummy,wdummy)
        if(lhpgl) call hpglfile (5,dummy,dummy,dummy,dummy)
      end if
      go to 200
c.....[tplo,v1]
c.... load deflection curve, phase diagram, etc.
c.... iswpl =  1 --> plot load         vs displacement
c.... iswpl =  2 --> plot displacement vs time
c.... iswpl =  3 --> plot velocity     vs displacement (phase portrait)
c.... iswpl =  4 --> plot reaction     vs displacement
c.... iswpl =  5 --> plot determinant  vs displacement
c.... iswpl =  6 --> plot velocity     vs time
c.... iswpl =  7 --> plot acceleration vs time
c.... iswpl =  8 --> plot load         vs time
c.... iswpl =  9 --> plot reaction     vs time
c.... iswpl = 10 --> plot determinant  vs time
c.... iswpl = 11 --> plot time         vs displacement
c.... iswpl = 12 --> plot stre         vs displacement
c.... iswpl = 13 --> plot stre         vs time
c.... iswpl = 14 --> plot valuse1      vs displacement
c.... iswpl = 15 --> plot valuse2      vs displacement
c.... iswpl = 16 --> plot valuse2      vs time
c.... iswpl = 17 --> plot valuse2      vs time
35    continue
      ro = dot(rotang,rotang,3)
cww   if(iso) then
cww     call drawmess('Hit any key and repeat statement',icolo,1)
cww     goto 4
cww   else if(clip) then
cww     call drawmess('Hit any key and repeat statement',icolo,1)
cww     call pzero(ct,3)
cww     goto 16
cww   else if(ro.gt.1.e-8) then
cww     call drawmess('Hit any key and repeat statement',icolo,1)
cww     l  = 43
cww     goto 43
cww   end if
      iswpl  = ct(1)
cww   imod   = 0
      iclear = 0
      call plopen
      if(iswpl.gt.0. and. iswpl.le.17) then
        call modscal(1)
        call plotdf(plo,ipl(2,1),ipl(2,2),nstedf,mkflg,mmc,mmst,incmk,
     1             iswpl,b,trans,trans(1+numnp*ndf),strea(1+numnp),
     2             id,numnp*ndf,numnp)
        call modscal(2)
      end if
      go to 160
c.... 36:[flux,v1,v2,v3,v4]
c....       v1 = 1  plot heat-flux
c....       v1 = 2  plot principal stresses
c....       v1 < 0  plot on deformed mesh
c....       v2      2=2D[def],3=3D only for stresses, open for heat-flux
c....       v3      scaling factor for vector length
c....       v4      scaling factor for tip length (only 2D)
36    c   = ct(1)
      if(c.lt.0.0d0.or.defo) then
        c = cs
      else
        c = 0.d0
      end if
      k1  = abs(ct(1))
      k2  = ct(2)
      rk3 = ct(3)
      if(rk3.eq.0.d0) rk3 = 1.0d0
      rk4 = ct4
      if(rk4.eq.0.d0) rk4 = 1.0d0
      call plopen
      call pppcol(2)
      if(plfl) then
        dbgtxt = 'PPLOTF: [stre] gen. array: np,numnp*npstr*ipr(plfl)'
        call ralloc(strea,numnp*npstr,'STRE',plfl)
      end if
      !k5 = np + numnp*ipr   replaced by strea(1+numnp)
      if(.not.fl(11).or.(icleas.eq.0)) then
        strea = 0.d0
        hflgu  = fa
        h3flgu = fa
        call formfe(b,dr,dr,dr,fa,fa,fa,fa,8,1,numel,1)
        call pltstr(strea,strea(1+numnp),numnp)
      end if
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      if(k1 .eq. 1) call pltflu (dr,strea(1+numnp),ndm,numnp,tra,vr,rk3)
      if(k1 .eq. 2) call pltmain(dr,strea(1+numnp),ndm,numnp,tra,vr,k2
     +                   ,rk3,rk4,ipgl)
      fl(11) = tr
      icleas = 1
      go to 150
c
c.... principal stress contours
c.... [pris,v1,v2,v3,v4]    plot filled /contour stresses
c                    v1 -   number of principal stress, neg. on deform. mesh
c                    v2 > 0 contour plot
c                    v2 <=0 filled plot
c                    v3 < 0 filled plot for v3=k.m for layer k pos m
c                    v3 > 0 contour plot without numbers
c                    v4 -   plot mesh in color v4
37    ipb = ct(3)
      if(plfl) then
        dbgtxt = 'PPLOTF: [stre] gen. array: np,numnp*npstr*ipr(plfl)'
        call ralloc(strea,numnp*npstr,'STRE',plfl)
      end if
      k1 = abs(ct(1))
      k2 = ct(2)
      rk2= ct(2)
      if(rk2.eq.0) rk2 = 1.0d0
      k5 = ct4
      klay = 0
      mlay = 1
      if(ipb.lt.0) then
        play = dabs(ct(3))
        klay = int(play)  ! number   of layer
        play = (play-klay)*10
        mlay = nint(play) ! position in layer
        if(mlay.eq.0) mlay=1
        ipb  = 0
      end if
      k3 = 1 + numnp
      if(.not.fl(11).or.(icleas.eq.0)) then
        strea = 0.d0
        hflgu  = fa
        h3flgu = fa
c....   stresses at nodes
        call formfe(b,dr,dr,dr,fa,fa,fa,fa,8,1,numel,1)
c....   average stresses at nodes
        call pltstr(strea,strea(k3),numnp)
c....   principal stress at nodes
        call prinstr(strea(k3),k1,numnp)
      end if

      k4 = k3 + numnp*(npstr-2) ! plot st(26)
      if(icleas.eq.0) then
        call rprint(strea(k4),ix,x,ndm,numnp,1,nen1,k2,idev)
      end if
c.... set number to 26
      k1=npstr-1

      c = 0.0
      if(ct(1).lt.0.0d0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      if(k2.le.0) then
        k2 = -k1
        k1 = 1
      end if
      call pltcon(dr,id,ie,mnix,strea(k4),idis,nie,ndm,1,nxd,nxn,nne,
     1            k1,k2,1,1,k5,cinv)
      icleas = 1
      fl(11) = fa
      goto 150
c 38  [forc] .... [forc,v1,v2,v3] plot forces
c                 (v1: force number) v1=12,13 defines plane to plot
c                 (v2: color number for plot) weg wenn alle elemente umgestellt!!
c                 (v3: <0 layer number)
c                 (v3: >0 scaling factor)
38    c = 0.0
      if(ct(1).lt.0.0d0.or.defo) c = cs
      scal = c
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      k1 = abs(ct(1))
      k2 = ct(2)
      if(k2.lt.0) k2 = 1
      k3 = 0
      iclear = 0
      call plopen
      call pppcol(k2)
c.... define plot direction
      kfor = ct(1)
      if(kfor.eq.12.or.kfor.eq.13.or.kfor.eq.-12.or.kfor.eq.-13) then
        ifor = kfor
        goto 200
      end if
      if(kfor.eq.0) then
        ifor = 13
        goto 200
      end if
      if(ct(3).le.0.0) klayf = abs(ct(3))
      cfp  = 0.0
      nfp  = ct(1) ! no to plot including sign for plot on mesh/def.mesh
      flfp = tr
      k3  = 1
c.... calculate extreme values
      xmaxf = -1.e+38
      xminf =  1.e+38
      if(k2.eq.0) call pppcol(abs(nfp))
      if(k2.ne.0) call pppcol(k2)
      call plot2d(ie,mnix,idis,b,dr,xl,nie,ndm,nxn,nxd,nne,k1,k2,k3)
      flfp = fa
      if(cfp.gt.0.0) then
        ffo = 1.
        if(ct(3).gt.0.0) ffo = ct(3)
        cfp=ffo*.05*max(xmax(1)-xmin(1),xmax(2)-xmin(2))/cfp
      end if
c.... Set colors new and plot force-palette
      call setcol(2,0)
      call pltconf(idev,xmaxf,xminf,k1,1)
c
c.... plot forces
      call plot2d(ie,mnix,idis,b,dr,xl,nie,ndm,nxn,nxd,nne,k1,k2,k3)
c.... Set colors back
      call setcol(1,0)
      go to 200
c
c.... set material number for plotting parts of stresses, contours
c.....[matn,k1,k2,k3]  k1 < 0 : set all values to 0
c                      k1 = 0 : set all values to i (plot for all numbers)
c                      k1 > 0 : set values k1 to k2, inc = k3
39    k1=ct(1)
      k2=ct(2)
      k3=ct(3)
c
      if(k1.lt.0) then
        fl(11) = fa ! calculate stresses new, e.g. in case of coupling plate/disk
        iclear = 0  ! clear window
      end if
      fl(11) = fa   ! calculate stresses new, e.g. in case of coupling plate/disk
c
      call plmaset(aipma,nummat,k1,k2,k3)
      if(k1.ge.0) then
c....   for hide 3D
        if(hide) then
          call pzero(ct,3)
          goto 29
        end if
c....   for hids shell
        if(hids) then
          call pzero(ct,3)
          ct(1) = 1
          if(defo) ct(1) = -1
          goto 62
        end if
      end if
      goto 200
c.... [logo,v1] plot logo in color v1
40    icol = ct(1)
      if(icol.eq.0) icol=2
      call plopen
      call plogo(icol)
      go to 200
c.... plot reactions
c.... [reac,v1,v2,v3]   plot reactions (separately!)
c                       (v1 =  i:      force(i) on undeformed mesh)
c                       (v1 = -i:      force(i) on   deformed mesh)
c                       (v1 =  0:      all forces (1-ndm) on undeformed mesh)
c                       (v1 = (ndf+1): all forces (1-ndm) on undeformed mesh)
c                       (v1 =-(ndf+1): all forces (1-ndm) on   deformed mesh)
c                       (v1 = (ndf+2): all moments(4-5/6) on undeformed mesh)
c                       (v1 =-(ndf+2): all forces (4-5/6) on   deformed mesh)
c-------------------------------------------------------
c      t=0:         loads from mate + single loads
c      after tang:  reactions
c-------------------------------------------------------
41    continue
      if(iadd.gt.0)goto 410
      c  = 0.d0
      k1 = abs(int(ct(1)))
      k1 = min(k1,ndf+2)
      if(ct(1).lt.0.0d0.or.defo) c = cs
      k2 = ct(2) ! tip
      rk3 = dabs(ct(3))
      k4 = ct(3)
      if(ct(3).eq.0) rk3 = 1.0d0
      if(fdrp)dbgtxt='PPLOTF:[reac]gen. array: numnp*ndf*ipr(fdrp)'
      call ralloc(drp,numnp*ndf,'pplotf-help-field',fdrp)
      call pzero(dr,numnp*ndf)
c.... dynamic loads F_dyn
      if(fl(9)) call ploadd(b,dr,massm(nxll),massm(nxu),massm
     +    ,dampm(ncll),dampm(ncu),damp,trans,jp,neq,numnp*ndf,fl(1),1)
c.... external loads from LOAD
      call pload(id,f,f0,dr,numnp*ndf,propq)
c.... move to full vector to include b.c. (dfl=tr in formfe)
      call pmovec(id,dr,drp,numnp*ndf)
      hflgu  = fa
      h3flgu = fa
c.... external loads from SLOA/QLOA  (uncompressed,dfl=true)
      call ploads(b,drp,propq,tr,fa,fa)
c.... internal loads F_int (uncompressed,dfl=true)
      call formfe(b,drp,drp,drp,fa,tr,fa,tr,6,1,numel,1)
410   call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      if(l.eq.59) goto 57
      call plopen
      if(k4.lt.0) then
c.....plot only nodes with reactions
        call pltload(dr,drp,id,ndm,ndf,numnp,2)
      else
c.....  plot arrows(separate for each dof) + nodes
        if(.not.fl(4)) call pmoves(drp,drp,numnp*ndf,-1.d0) ! -R after tang
        call pltfor2(dr,drp,bang,id,ndm,ndf,numnp,k2,rk3,
     +                tra,vr,1,k1,ix,nen1,2)
      end if
c
cww   call pltrea(dr,drp,ndf,ndm,nr1,nr2,scfac,tra,vr,bang,ix,nen1)
      go to 150
c
c.... window with the cross hairs
c.... [clip] (possible with cart and rot) only one times if symm active
42    if(iadd.ne.0) goto 150
cww      if(iso.or.kpers.eq.1) then
cww        call drawmess('Clip not possible, Use CART+ROT ',icolo,1)
cww        goto 200
cww      end if
      nzm1  = -1
      clip = tr
      go to 160
c
c.... rotation of views, rescale to fit
c.... [rot..i,v1,v2] i=0 rescale, i=1-->rot_x i=2-->rot_y i=3-->rot_z
c....                    v1   rotation angle in degree
c....                    v2 > show axis on screen before plotting
43    i = l - 43
      if(ndm.lt.3) then
        call drawmess('Rot only for 3-D poss. ',icolo,1)
        goto 200
      end if
      if(i.gt.0) then
        rotang(i) = rotang(i)+ct(1)
        call pltrns(i,ct,tra,vr)
      else
c....   unity matrix
        call pzero(tra,9)
        do 430 i = 1,3
           tra(i,i) = 1.0
           vr(i)    = 0.0
430     continue
        call pzero(rotang,3)
        iclear = 0
        call plopen
      end if
      c = 0.0
      call pdefm(x,b, c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      call frame(dr,ndm,numnp,1,clip)
      iclear = 0
      k2=ct(2)
      if(k2.gt.0) then
c....   plot axis on screen
        iclear = 0
        call plopen
        if(ibor.gt.0) call pltrot(tra,ndm,rotang,1)
        iclear = 0
      end if
      go to 200
c
c 44  [rot1] under 43
c 45  [rot2] under 43
c 46  [rot3] under 43
c.... plot eigenvector  of extended system
c.... [evex,v1,v2,v3,v4 (v1),v2,v3,v4 see eigv(22)
47    k1 = -1
      k11= ct(1)
      k5 = ct4
      if(.not.extflg) then
        call drawmess('Compute EV of ext. system first ',icolo,1)
        goto 200
      end if
      call pmovec(id,extkh,dr,numnp*ndf)
      goto 221
c
c.... move plot on screen
c.... [move,v1,v2,v3,v4]
c     v1 =  no of axis [1,2]
c     v2 =  delta x(y) [0-1], def=0.02
c     v3 <=0 interactive with mouse, >0 only one step batch
c     v1 > 0 mesh, v1 < 0 deformed mesh
c     v3 = 1
c     v3 = 2 no mesh on screen
c     one mesh plot should be done before
c
48    iaxm   = abs(ct(1))
      if(iaxm.eq.0.d0) iaxm = 1
      ddxm  = ct(2)
      ixm   = ct(3)
      if(ddxm.eq.0.d0) ddxm = 0.02d0
      dxm     = 0.d0
      dym     = 0.d0
      if(iaxm.eq.1)dxm = ddxm
      if(iaxm.eq.2)dym = ddxm
      c=0.d0
      if(ct(1).lt.0.0d0.or.defo)  c = cs
      if(ixm.gt.0) go to 481
      if(idev.lt.4) goto 200
c.... interactive
      write(*,*)
     +'left/down = Left B., right/up = Right B., end = both B. of mouse'
      istat  = 0
480   iclear = 0
c.... move
      if(istat.eq.1)  then
c....   left/down
        deltax = deltax - dxm
        deltay = deltay - dym
      else if(istat.eq.2) then
c....   right/up
        deltax = deltax + dxm
        deltay = deltay + dym
      end if
c...  calculate mesh/defm
      if(ixm.ne.2) then
        call plopen
        call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
        call plxtrn(dr,tra,vr,ndm,numnp)
        call frame(dr,ndm,numnp,1,clip) ! not scale
c...    plot mesh/defm
        call pppcol(4)
        if(ct(1).lt.0.0d0.or.defo) call pppcol(5)
        if(.not.allocated(mea)) allocate(mea(10*numnp))
        call pline(dr,ie,mnix,mea,numnp,nne,ndm,nxd,nxn,nie,0.d0,tr)
      end if 
c...  decide
      call plmousec(istat)
      if(istat.lt.3) goto 480
      write(*,*) 'dx=',deltax,' dy=',deltay
      goto 200
c.... only one step
481   deltax = deltax + dxm
      deltay = deltay + dym
c...  calculate mesh/defm
      if(ixm.ne.2) then
        iclear = 0
        call plopen
        call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
        call plxtrn(dr,tra,vr,ndm,numnp)
c...    plot mesh/defm
        call pppcol(4)
        if(ct(1).lt.0.0d0.or.defo) call pppcol(5)
        if(.not.allocated(mea)) allocate(mea(10*numnp))
        call pline(dr,ie,mnix,mea,numnp,nne,ndm,nxd,nxn,nie,0.d0,tr)
      end if 
      goto 200
c
c     [titl] put title of problem on screen
49    call plopen
      call plttit
      go to 200
c***1 [plof,nplt,iplt,ldsc] write to cgm/gdf - file for IBM PHIGS
c                                 write to ps/cgm  - file for  HP GKS
c           nplt = 1...9 open  ps/cgm/gdf file fplt_nplt.ps(cgm,gdf)
c           nplt = (0)   close ps/cgm/gdf file
c           iplt = 1 --> GKS: ps(default)  PHIGS: cgm(bin)(default)
c           iplt = 2 --> GKS: cgm(bin)     PHIGS: gdf
c           ldsc = 1 --> A4 landscape only for PS (default)
c           ldsc = 2 --> A4 portrait  only for PS
c-------------------------------------------------------------------
c***2 [plof,nplt] write screen to File: fplt_nplt.pcx,  for FTN77
c           nplt = 1...9
c-------------------------------------------------------------------
50    nplt = ct(1)
      call drawmess('Use the macro PRIN for Plotting',icolo,1)
c     nplt = max(0,min(nplt,9))
c     if(idev.le.2) then
c       iplt = ct(2)
c       iplt = max(1,min(iplt,2))
c       ldsc = ct(3)
c       ldsc = max(1,min(ldsc,2))
c       call pltcgm(nplt,fplt,iplt,ldsc,istruc)
c       if(nplt.ne.0) iclear = 0      !begin new in plopen
c     else if(idev.eq.4) then
c       call pltpcx(nplt,fplt,idev)
c     end if
      goto 200
c
c.... plot mesh and deformed mesh with hidden line
c.... [hmsh,v1,v2,v3] -(v1 < 0  on deformed mesh)
c                      (v1 - +-1->forward(default) +-2->backward)
c                      (v2 - color number line ,default = 32(black) mesh
c                      (v3 - color number interior part,
c                            default = 4(cyan)-> mesh, 5(green)-> defm)
51    if(outf) dbgtxt = 'PPLOTF: [hmsh] gen. array: me,10*numnp*1(outf)'
      outf=fa
      if(.not.allocated(mea)) allocate(mea(10*numnp))
      c = 0.0
      if(ct(1).lt.0.0d0.or.defo) c = cs
      k4 = 1
      if(dabs(ct(1)).eq.2.0d0) k4 = 2
      k2 = ct(2)
      if(k2.eq.0) k2 = 32
      k3 = ct(3)
      if(k3.eq.0.and.ct(1).ge.0.0d0) k3 = 4
      if(k3.eq.0.and.ct(1).lt.0.0d0) k3 = 5
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      call plopen
      call plot2dh(ie,mnix,dr,xl,nie,ndm,nxn,nxd,nne,k2,k3,k4,idis)
      go to 150
c
c.... plot only in color black on white(for b&wplots)
c.... [mono,v1] -(v1 = 0  plot in standard colors)
c                (v1 > 0  plot only in black on white)
52    imono = 0
      icolo = 7
      if(ct(1).gt.0.0) then
        imono = 1
      end if
      call setcol(1,imono)
      iclear = 0
      call plopen
      goto 200
c
c.... plot damage values at gauss points
c.... [dmag] plot damage values at gauss points, symm not implemented
53    ipld = 1
      call plopen
      call pppcol(2)
      hflgu  = fa
      h3flgu = fa
      call formfe(b,dr,dr,dr,fa,fa,fa,fa,4,1,numel,1)
      ipld = 0
      goto 150
c
c.... plot 1.+ 2. displacement in cartesian(default) or polar directions
c.... [pola,i] i=0 --> cartesian, i = 12,13,23 --> polar in 1/2 1/3 2/3
54    ipola = ct(1)
      if(ipola.eq.0.or.ipola.eq.12.or.ipola.eq.13.or.ipola.eq.23) then
        write(*,2012) ipola
      else
        write(*,2011)
        ipola = 0
      end if
      goto 200
c
c.... [erro,v1,v2] plot distribution of errors
c.... v1 = eval 'eval' percent of energy (default = 5 percent)
c.... v2 = 1 energy-norm (default)
c.... v2 = 2 L_2   -norm
c.... v2 = 3 Y0   
55    if(.not.fl(11)) then
        call drawmess('Compute stresses before ERRO ',icolo,1)
        goto 200
      end if
      call ralloc(e_ome,numerr*numel,'REMESH-erron/new',tflb1)! also remesh1
      iet(1) = 1          ! switch in element for 'plot>erro' mode
      iet(2) = ct(2)      ! switch in element for error type
      iet(2) = min(max(1,iet(2)),numerr)
      eval = ct(1)
      if(eval.eq.0.0d0) eval = 5.d0
      iclear = 0
      call plopen
      call plterr(b,dr,ndf,numel)
      go to 200
c.... [back,v1] set background color (v1 = 1 black, v1 = 0 white)
56    iback = ct(1)
      red   = ct(2)
      green = ct(3)
      blue  = ct4
      call pltback(idev,iback,red,green,blue)
      iclear = 0
      call plopen
      goto 200
c
c.... plot values along a user defined line
c.... [dplo,k1,k2,k3] 57 plot displacement k1
c.... [splo,k1,k2,k3] 58 plot stress       k1
c.... [rplo,k1,k2,k3] 59 plot reaction   k1
c.... [eplo,k1,k2,k3] 60 plot displacement k1 of eigenvector k3
c                k1 < 0 plot on deformed mesh
c                k2 = 1 plot other value in same diagram  with    rescaling
c                k2 = 2 plot other value in same diagram  without rescaling
c                k2 < 0 print values on screen and in output file
c                k2 = 0 plot new   value in diagram
c                k3 < 0 stress for layer abs(k3)
57    continue
      numnp2=0.5*numnp
      if(fdpl) then
        dbgtxt = 'PPLOTF: [rdplo] gen. array: numnp*2(fdpl)'
        dbgtxt = 'PPLOTF: [idplo] gen. array: nnemax (fdpl)'
      end if
      fdpl = fa
      if(.not.allocated(rndplo)) allocate(rndplo(numnp))
      if(.not.allocated(indplo)) allocate(indplo(nnemax))

      k1 = abs(ct(1))
      k2 = ct(2)
      k3 = ct(3)
      c  = 0.0
      if(ct(1).lt.0.0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      if(l.eq.57) then
c....   plot curve of displacements
        if(fdrp)
     +  dbgtxt='PPLOTF:[dplo]gen. array: numnp*ndf*ipr(fdrp)'
        call ralloc(drp,numnp*ndf,'pplotf-help-field',fdrp)
        call pmove(b,drp,numnp*ndf)
        call panglb(drp,bang,numnp,ndf)
        if(ipola.ne.0) call dispola(x,drp,ndm,ndf,numnp,ipola)
        k1 = max(1,min(k1,ndf))
        call plotn(drp,dr,mnix,idis,rndplo,indplo,k1,k2,k3,
     1             ndf,nxn,nxd,ndm,numnp,numnp2,nne,s0(1),s0(2),'d')
      else if(l.eq.58) then
c....   plot curve of stresses
        call plotn(strea(k4),dr,mnix,idis,rndplo,indplo,k1,k2,k3,1,
     1             nxn,nxd,ndm,numnp,numnp2,nne,s0(1),s0(2),'s')
      else if(l.eq.59) then
c....   plot curve of reactions
        k1 = max(1,min(k1,ndf))
        call plotn(drp,dr,mnix,idis,rndplo,indplo,k1,k2,k3,
     1             ndf,nxn,nxd,ndm,numnp,numnp2,nne,s0(1),s0(2),'r')
      else if(l.eq.60) then
c....   plot curve of displacements of eigenvector k3
        k1 = max(1,min(k1,ndf))
        k3 = min(mfmax,max(1,k3))
        if(mfmax.eq.0) then
          call drawmess('Comp. Eigenvector first ',icolo,1)
          goto 200
        end if
        if(fdrp)
     +  dbgtxt='PPLOTF:[eplo]gen. array: numnp*ndf*ipr(fdrp)'
        call ralloc(drp,numnp*ndf,'pplotf-help-field',fdrp)
        call pmovec(id,eigv(1+(k3-1)*neq),drp,numnp*ndf)
        call panglb(drp,bang,numnp,ndf)
        if(ipola.ne.0) call dispola(x,drp,ndm,ndf,numnp,ipola)
        call plotn(drp,dr,mnix,idis,rndplo,indplo,k1,k2,k3,
     1             ndf,nxn,nxd,ndm,numnp,numnp2,nne,s0(1),s0(2),'e')
      end if
      goto 200
c 58  [splo] under 12, then 57
c 59  [rplo] under 41, then 57
c 60  [eplo] under 57
c
c.... [xtic,k1,k2,k3] plot axis with tics at vmin(i)
c            ki= number of tics on axis i
c            k1<0 on deformed mesh
c            k2<0 plot cubus of system region
61    k1=abs(ct(1))
      k2=ct(2)
      k3=ct(3)
      c = 0.0
      if(ct(1).lt.0.0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call pltcord(dr,ndm,numnp,k1,k2,k3)
      goto 200
c
c.... [hids,k1] sort elements due to distance to viewpoint
c               for hidden line plot of 2-d elements,not with symm!
c               for hidden line plot of 3-d elements
c                              (after calculating the faces with Hide)
c.... k1 = 0 reset to old values
c     k1 > 0 on          mesh
c     k1 < 0 on deformed mesh
62    k1 = ct(1)
      hids = fa
      if(k1.ne.0) hids = tr
      k1 = abs(ct(1))
      if(l.eq.29) go to 621 ! entry from 29 hide
      if(kpers.eq.0) then
        call maxcor(x,ndm,numnp)
        call plthid4
      end if
      c    = 0.0
      if(ct(1).lt.0.0d0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
621   if(hidw2) then 
        dbgtxt = 'PPLOTF: [hids] gen. array: DIST nnemax*2(hidw2)'
        call ralloc(dist,nnemax,'DIST',hidw2)
      end if
      call plthid1(mnix,x,xl,idis,dist,ndm,nxn,nxd,nne,k1)
      iclear = 0
      call plopen
      goto 200
c
c.... [defo,k1] k1=0 (default) plot all results on undeformed mesh
c               k1>0           plot all results on   deformed mesh
63    k1=ct(1)
      if(k1.eq.0) then
        defo=fa
      else
        defo=tr
      end if
      goto 200
c
c.... plot tied nodes
c.... [tie ,v1,v2,v3] plot tied nodes
c                  (v1.lt.0 plot nodes on deform. mesh)
c                  (v2.gt.0 plot nodes numbers)
c                  (v2.lt.0 plot only node with number v2)
c                  (v3.ne.0 plot in color v3 (default 2 (red))
64    c = 0.0
      if(ct(1).lt.0.0d0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      k2 = ct(2)
      k3 = ct(3)
      if(k3.le.0) k3 = 2
      call plopen
      call pppcol(k3)
      call plttie(dr,gtie,ndm,numnp,k2)
      go to 150
c.... plot eigenvector in animation mode (without symm)
c.... [evan,v1,v2,v3] v1  - no. of eigen vector
c                     v1 >= 0 --> mesh   v1 < 0 --> hmesh
c                     v2 >= 0 hidden line from 1 to numel
c                     v2 <  0 hidden line from  numel to 1
c....                 v3 = 'speed' default = 0
c....1: animate the mode shape 2 times  +1 to -1,inc = 2/8
c....2: animate the mode shape 6 times  +1 to -1,inc = 2
65    k3 = ct(1)
      k1 = abs(ct(1))
      k1 = min(mfmax,max(1,k1))
      if(mfmax.eq.0) then
        call drawmess('Comp. Eigenvector first ',icolo,1)
        goto 200
      end if
      iclear=0
      k4 = 1
      k2 = ct(2)
      if(k2.lt.0) k4 = 2
      if(outf) dbgtxt = 'PPLOTF: [evan] gen. array: me,10*numnp*1(outf)'
      outf=fa
      if(.not.allocated(mea)) allocate(mea(10*numnp))
      if(fdrp)dbgtxt='PPLOTF:[evan]gen. array: numnp*ndf*ipr(fdrp)'
      call ralloc(drp,numnp*ndf,'pplotf-help-field',fdrp)
      call pmovec(id,eigv(1+(k1-1)*neq),drp,numnp*ndf)
c.... animate the mode shape 2 times  +1 to -1,inc = 2/8
      do i  = 1,8
        c = cs*am(i)
        call pdefm(x,drp,c,bang,ndm,ndf,numnp,dr,xfac)
        call plxtrn(dr,tra,vr,ndm,numnp)
c....   plot
        call plopen
        if(k3.ge.0) then
          call pppcol(8)
          call pline(dr,ie,mnix,mea,numnp,nne,ndm,nxd,nxn,nie,0.0d0,tr)
        else if(k3.lt.0) then
          call plot2dh(ie,mnix,dr,xl,nie,ndm,nxn,nxd,nne,32,27,k4,idis)
        end if
        call plclos
c.....  sleep
        tsl = ct(3)
        tsl = max(0.0,tsl)
        call sleepw(tsl)
c....   wipe
        call plopen
        if(k3.ge.0) then
          call pppcol(31)
          call pline(dr,ie,mnix,mea,numnp,nne,ndm,nxd,nxn,nie,0.0d0,tr)
        else if(k3.lt.0) then
          call plot2dh(ie,mnix,dr,xl,nie,ndm,nxn,nxd,nne,32,32,k4,idis)
        end if
        call plclos
      enddo
c.... animate the mode shape 6 times  +1 to -1,inc = 2
      do ii = 1,6
      do i  = 1,5,4
       c = cs*am(i)
       call pdefm(x,drp,c,bang,ndm,ndf,numnp,dr,xfac)
       call plxtrn(dr,tra,vr,ndm,numnp)
c....  plot
       call plopen
       if(k3.ge.0) then
        call pppcol(8)
      call pline(dr,ie,mnix,mea,numnp,nne,ndm,nxd,nxn,nie,0.d0,tr)
       else if(k3.lt.0) then
      call plot2dh(ie,mnix,dr,xl,nie,ndm,nxn,nxd,nne,32,27,k4,idis)
       end if
       call plclos
c..... sleep
       call sleepw(tsl)
c....  wipe
       call plopen
       if(k3.ge.0) then
        call pppcol(31)
      call pline(dr,ie,mnix,mea,numnp,nne,ndm,nxd,nxn,nie,0.d0,tr)
       else if(k3.lt.0) then
      call plot2dh(ie,mnix,dr,xl,nie,ndm,nxn,nxd,nne,32,32,k4,idis)
       end if
       call plclos
      enddo
      enddo
c.... plot the final shape
      call plopen
      c = cs*am(5)
      call pdefm(x,drp,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      if(k3.ge.0) then
       call pppcol(8)
      call pline(dr,ie,mnix,mea,numnp,nne,ndm,nxd,nxn,nie,0.0d0,tr)
      else if(k3.lt.0) then
      call plot2dh(ie,mnix,dr,xl,nie,ndm,nxn,nxd,nne,32,27,k4,idis)
      end if
      call pleigvt(eigd(k1),k1)
      go to 200
c 66  [base,v1] under 10
c

c     set all plot quantities to initial values
c     [init,n1] set plot quantities to initial values
c     n1=0 reset all values
c     n1=1 interactive input
c     n1=2 set userdefined values for automatique input
67    k1 = ct(1)
      if(k1.eq.1) then
        icv = 1
      else if(k1.eq.2) then
        if(ior.lt.0) write(*,2005)
        icv = 0
        call dinput(contv,3)
      else  ! same definitions as at beginning!
c....   cart
        iso   = fa
c....   perspective
        kpers = 0
        call pzero (eold,3)
        call pzero (vold,3)
        vold(3) = 1.0d0
c....   rot
        call pzero(tra,9)
        do i = 1,3
          tra(i,i) = 1.0
          vr(i)    = 0.0
        end do
        call pzero(rotang,3)
c....   elements hide
        mnix(1:nen1*numel)=ix(1:nen1*numel)
        nxd   = nen1
        nxn   = nen
        nne   = numel
        nnemax  = 3*numel
        hide  = fa
        hids  = fa
c....   deformation
        defo  = fa
        cs    = 1.0
c....   zoom
        nzm1 =  0
        nzm2 =  0
        nzm3 =  0
        call pzero(xzm,6)
c....   clip
        clip = fa
        fact = 1.d0
c....   move
        deltax  = 0.d0
        deltay  = 0.d0
c....   xscal coor
        xfac(1) = 1.0d0
        xfac(2) = 1.0d0
        xfac(3) = 1.0d0
c....   sym
        nsym = 0
        msym = 0
c....   plot pola
        ipola = 0
c....   colors etc
        icgm    = 0
        imono   = 0
        iback   = 0
c...    Contour values
        icv     = 0
        call pzero (contv,3)
c....   maxi stress
        nmn     = 0
        nmx     = 0
c....   frame
        ifrm    = 0
c....   plot/prin options
        ipgl    = 1
        nnope   = 0
        iclear  = 0
        iopl    = iop
        nexte   = 0
        lhpgl   = fa
        lps     = fa
c....   hptext
        isize   = 1
c....   set frame
        c = 0.0
        call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
        call plxtrn(dr,tra,vr,ndm,numnp)
        call frame(dr,ndm,numnp,1,clip)
      end if
      goto 200
c
c     [xsca,v1,v2,v3]
c.... scale coordinates and displacements
68    do i = 1,3
        xfac(i) = ct(i)
        if(xfac(i).eq.0.0d0) xfac(i) = 1.0d0
      enddo
cww   if(c.eq.0.0d0) c = cs
      c=cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      call frame(dr, ndm,numnp,1,clip)
      iclear = 0
      goto 200
c
c.... rotation of views by mouse-click for mesh/defm (only idev=3/4)
c.... start from actual values of rot
c.... [rotm,v1,v2]
c                  v1 > 0 mesh, v1 < 0 deformed mesh
c                  v1=i --> rotation arround axis i
c                  v2   --> incremental  rotation angle(def.=5)
69    if(idev.eq.1.or.idev.eq.2) goto 200
      ira = dabs(ct(1))
      ira = max(1,min(3,ira))
      call pzero (dra,3)
      dalpha = ct(2)
      if(dalpha.eq.0.d0) dalpha = 5.d0
      istat  = 0
      c = 0.0d0
      if(ct(1).lt.0.0d0.or.defo)  c = cs
      if(outf) dbgtxt = 'PPLOTF: [rotm] gen. array: me,10*numnp*1(outf)'
      outf=fa
      if(.not.allocated(mea)) allocate(mea(10*numnp))
      if(ndm.lt.3) then
        call drawmess('Rot only for 3-D poss. ',icolo,1)
        goto 200
      end if
      if(idev.eq.3)
     +  write(*,*) '+ = left, - = right, end = double left button mouse'
      if(idev.eq.4)
     +  write(*,*) '+ = left, - = right, end = both buttons  mouse'

690   iclear = 0
c.... rotate
      if(istat.eq.1)  then
        dra(1) =  dalpha
        rotang(ira) = rotang(ira)+dalpha
      else if(istat.eq.2) then
        dra(1) = -dalpha
        rotang(ira) = rotang(ira)-dalpha
      end if
      call pltrns(ira,dra,tra,vr)
      call plopen
c...  calculate mesh/defm
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      call frame(dr,ndm,numnp,1,clip)
c...  plot mesh/defm
      call pppcol(4)
      if(ct(1).lt.0.0d0.or.defo) call pppcol(5)
      call pline(dr,ie,mnix,mea,numnp,nne,ndm,nxd,nxn,nie,0.d0,tr)
c...  decide
      call plmousec(istat)
      if(istat.lt.3) goto 690
      go to 200
c
c.... plot filled/contour of velocities and accelerations
c.70. [velo,v1,v2,v3,v4]
c.71. [acce,v1,v2,v3,v4]
c                    v1 -    number of velo/acce., neg. on deform. mesh
c                    v2 >  0 number of acce/velo. values for plot
c                    v2 <= 0 filled plot
c                    v3.ne.0 contour plot without numbers
c                    v4 -  plot mesh in color v4
70    continue
      if(.not.fl(9)) then
         call drawmess('Problem not dynamic ',icolo,1)
         goto 200
      end if
      if(iadd.gt.0)goto 710
      k2 = ct(2)
      n = abs(k2)
      n = max(1,n)
      i = max(1,abs(int(ct(1))))
      i = min(i,ndf)
      ipb = ct(3)
      k5  = ct4
      if(fdrp)dbgtxt='PPLOTF:[velo]gen. array: numnp*ndf*ipr(fdrp)'
      call ralloc(drp,numnp*ndf,'pplotf-help-field',fdrp)
                   nu = 1 +   numnp*ndf ! at n+1
      if(nop.eq.5) nu = 1 + 2*numnp*ndf ! at n+1/2
      if(l.eq.70) then
        call pmovec(id,trans,drp,numnp*ndf)
        immc = 3
      else if(l.eq.71) then
        if(nop.eq.2) then
          call drawmess('No Accelerations for Euler Backward',icolo,1)
          goto 150
        end if
        call pmovec(id,trans(nu),drp,numnp*ndf)
        immc = 4
      end if
      call panglb(drp,bang,numnp,ndf)
      c = 0.0
      if(ct(1).lt.0.0d0.or.defo) c = cs
710   continue
      if(icleas.eq.0) then
        if(ipola.ne.0) call dispola(x,drp,ndm,ndf,numnp,ipola)
        call rprint1(drp,ix,x,ndm,numnp,ndf,n,nen1,i,idev)
      end if
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      if(k2.le.0) then
        n = -i
      end if
      call pltcon(dr,id,ie,mnix,drp,idis,nie,ndm,ndf,nxd,nxn,nne,
     +            i,n,immc,1,k5,cinv)
      icleas = 1
      go to 150
c
c.... plot nodes+ numbers with neg.Dii etc
c.... [ndii,v1]    v1=1  nodes with Dii <  0 (default)
c                  v1=2  nodes with Dii =  0
c                  v1=3  nodes with Dii << 1 (loss of digits)
c                  v1 < 0 on deformed mesh
c72     continue
72    k1 = abs(ct(1))
      k1 = max(1,min(k1,3))
      c = 0.0
      if(ct(1).lt.0.0d0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      call plopen
      call pppcol(2)
      call pltnod1(dr,ndm,numnp,k1)
      go to 200
c
c.... plot linked nodes
c.... [link,v1,v2,v3]
c              plot linked nodes
c                  v1 < 0  all linked nodes
c                  v1 = 0  nodes+all link.cond.
c                  v1 = i  nodes+link.cond. i
c                  v2 > 0  node numbers
c                  v3  = scaling factor for link.cond. (default = 1)
73    if(.not.allocated(link1)) goto 150
      c = 0.0
      k1 = ct(1)
      k2 = ct(2)
      rk3 = ct(3)
      if(rk3.eq.0.d0) rk3=1.d0
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      call plopen
c.... plot nodes
      call pltlink1(dr,ndm,numnp,link1,k2)
      if(k1.lt.0) go to 150

c.... plot linking conditions
      k4=k1
      if(k1.eq. 0) k4 = 3
      call pppcol(k4)
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call pltlink2(link2,dr,bang,ndm,ndf,numnp,k1,rk3,tra,vr)
      go to 150
c
c.74 [avel] plot arrows velocity         --> 10
c.75 [aacc] plot arrows acceleration     --> 10
c.76 [aeig] plot arrows eigenvector      --> 10
c.77 [aeve] plot arrows eigenvector e.s. --> 10
c
c... [draw,v1,v2,v3] draw line  node1->node2 /interactive/coor1->coor2
c....      v2 = 1. node
c....      v3 = 2. node
c....      v1 = 0 node 1 to node 2 (default)
c....      v1 =-1 node 1 to node 2 on deformed mesh
c....      v1 = 1 interactiveley
c....      v1 = 2 with coordinates
78    k1 = ct(1)
      k2 = ct(2)
      k3 = ct(3)
      c = 0.0
      if(ct(1).lt.0.0d0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      call plopen
      call pltdraw(dr,ndm,numnp,k1,k2,k3)
      goto 200
c
c...  [man] get the complete help manual
79    call pman(4,'titl')
      goto 200
c
c.... plot refined mesh
c.... [rmsh,v1,v2,v3] - (v1 < 0 on deformed mesh)
c                             (v2 - plot for material no. v2)
c                             (v3 - color number,default = 5(green))
c     kk new number of nodes
80    continue
      allocate(rtm(ndm*kk))
      allocate(itm( 10*kk)) ! see mea, 10=estimate
      rtm = 0.d0
      itm = 0
      k1 = ct(3)
      if(k1.le.0) k1 = 5
      call matcop(adanx,kk,ndm,rtm)
      call plxtrn(rtm,tra,vr,ndm,kk)
      call plopen
      call pppcol(k1)
      call pline(rtm,ie,m(melem),itm,kk,ketemp,ndm,nxd,nxn,nie,ct(2),tr)
      deallocate(itm)
      deallocate(rtm)
      go to 150
c
c.... sleep(wait) for v1 seconds (for demos)
c.... [sleep,v1]
81    tsl = ct(1)
      tsl = max(0.0,tsl)
      if(idev.eq.1) call drawmess('Sleep not impl. on IBM',icolo,1)
      call sleepw(tsl)
      go to 200
c     Plot layout of upper profile
c.... [prof,v1]
c     [prof],  - Plot layout of upper profile
c     [prof],1 - Plot layout of total profile ony for solv,0
82    iclear  = 0
      call plopen
      call pltprf(jp,neq,ct(1).ne.0.0d0)
      go to 200
c
c     [cent]er,s0-1,s0-2  - center graphics in window, similar to move, v1=v2=0 then via mouse
83    if(abs(ct(1))+abs(ct(2)).ne.0.0d0) then
        s0(1) = ct(1)
        s0(2) = ct(2)
        if (ct(1).eq.-1.and.ct(2).eq.-1) then
          s0(1) = 0.5d0
          s0(2) = 0.5d0
        endif
      else
        call plopen
        call center_mouse(s0(1),s0(2),idev)
      end if
      iclear = 0
      go to 200
c
c     [sect]ion,v1,v2,v3 - plot values for cross sections
84    call secthp(ix,b,dr,ct,ct4,cs)
      call plclos
      call frame(x,ndm,numnp,1,clip)
      go to 200
c
c     [copy]   - copy plot window to clipboard for further use in windows programs
c     [copy,1] - copy total window to clipboard (default)
c     [copy,2] - copy part of window defined by mouse
85    k1 = ct(1)
      if(k1.eq.0) k1 = 1
      call copyclip(k1)
      go to 200
c
c     [maxi],n1 - plot the nodes with max and min value of the
c                 last filled plot
c     n1<0 on deformed mesh
86    c = 0.0
      if(ct(1).lt.0.0d0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      call plopen
      call pppcol(1)
      call pltmaxi(dr,ndm,numnp)
      go to 200
c.... plot residuum components
c.... [resi,v1,v2,v3,v4]   plot filled /contour residua
c                    v1 -    number of resid., neg. on deform. mesh
c                    v2 >  0 number of resid. values for plot
c                    v2 <= 0 filled plot
c                    v3.ne.0 contour plot without numbers
c                    v4 -  plot mesh in color v4
87    continue
      if(iadd.gt.0)goto 870
      c = 0.0
      if(ct(1).lt.0.0d0.or.defo) c = cs
      k2  = ct(2)
      n   = abs(k2)
      n   = max(1,n)
      i   = max(1,abs(int(ct(1))))
      i   = min(i,ndf)
      ipb = ct(3)
      k5  = ct4
      if(fdrp)dbgtxt='PPLOTF:[resi]gen. array: numnp*ndf*ipr(fdrp)'
      call ralloc(drp,numnp*ndf,'pplotf-help-field',fdrp)
c.... external loads F_ext
      call pzero(dr,neq)
      call ploads(b,dr,propq,fa,fa,fa)
      call pload(id,f,f0,dr,numnp*ndf,propq)
c.... dynamic loads F_dyn
      if(fl(9)) call ploadd(b,dr,massm(nxll),massm(nxu),massm
     +                    ,dampm(ncll),dampm(ncu),dampm,trans,jp,neq
     +                    ,numnp*ndf,fl(1),1)
c.... move to full vector to include b.c. (dfl=tr in formfe)
      call pmovec(id,dr,drp,numnp*ndf)
      hflgu  = fa
      h3flgu = fa
c.... internal loads F_int
      call formfe(b,drp,drp,drp,fa,tr,fa,tr,6,1,numel,1)
c.... test if tang has been used
      if(.not.fl(4)) then
c....   multiply all terms of R by -1
        call pmoves(drp,drp,numnp*ndf,-1.d0)
      end if
      call panglb(drp,bang,numnp,ndf)
870   continue
      if(icleas.eq.0) then
        if(ipola.ne.0)  call dispola(x,drp,ndm,ndf,numnp,ipola)
        call rprint1(drp,ix,x,ndm,numnp,ndf,n,nen1,i,idev)
      end if
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      if(k2.le.0) then
        n = -i
      end if
      call pltcon(dr,id,ie,mnix,drp,idis,nie,ndm,ndf,nxd,nxn,nne,
     1            i,n,6,1,k5,cinv)
      icleas = 1
      go to 150
c
c.... plot nodes with angles
c.... [angl,v1,v2,v3] plot nodes with angles
c                  (v1.lt.0 plot nodes on deform. mesh)
c                  (v2.gt.0 plot nodes numbers)
c                  (v2.lt.0 plot only node with number v2)
c                  (v3.ne.0 plot local basis,v3=length,def=1)
88    c = 0.0
      if(ct(1).lt.0.0d0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      k2 = ct(2)
      call plopen
      call pppcol(2)
      if(ct(3).eq.0) then
        call pltangl(dr,bang,ndm,numnp,k2)
      else
        call pltangl1(dr,bang,tra,vr,ndm,numnp,ipgl,ct(3))
      end if
      go to 150
c
c
c.... plot nodes used with rsum
c.... [rsum,v1,v2,v3] plot nodes used with rsum
c                  (v1.lt.0 plot nodes on deform. mesh)
c                  (v2.gt.0 plot nodes numbers)
c                  (v2.lt.0 plot only node with number v2)
c                  (v3.ne.0 plot in color v3 (default 2 (red))
89    c = 0.0
      if(ct(1).lt.0.0d0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      k2 = ct(2)
      k3 = ct(3)
      if(k3.le.0) k3 = 2
      call plopen
      call pppcol(k3)
      if(nfs1.eq.0) then
        call drawmess('No points available',6,1)
      else 
        call pltrsum(dr,irpt,ndm,numnp,k2)
      end if
      go to 150
c
c...  [pnod,v1,v2,v3] pick up node number graphically
c           v1 .lt.0. on deformed mesh
c           v2 .ne.0 write additional data
c           v3 .ne.0 write results on output file
c
90    continue
      if(ct(1).lt.0.0d0.or.defo) c = cs
      k1 = ct(1)
      k2 = ct(2)
      k3 = ct(3)
      if(plfl) then ! in case of use before PLOT,STRE
        dbgtxt = 'PPLOTF: [stre] gen. array: np,numnp*npstr*ipr(plfl)'
        call ralloc(strea,numnp*npstr,'STRE',plfl)
      end if
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      call plopen
      call pltpnod(x,dr,b,strea(1+numnp),bang,ix,nen1,ndm,ndf,numnp,
     +     fl(11),k2,k3,ipola)
      goto 200
c
c...  [pele] pick up element number graphically
91    continue
      if(ct(1).lt.0.0d0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      call plopen
      call pltpele(dr,mnix,ie,idis,nxn,nxd,nie,ndm,numnp,nne)
      goto 200
c
c     [str1] .... [str1,v1] plot stresses from 1 point integration without averaging
c                    v1 -   number of stress
c...  symm not possible, defm not possible
92    k1 = abs(ct(1))
      iclear = 0
      call plopen
      nfp = k1
c.... calculate extreme values
      flfp = tr
      xmaxf = -1.e+38
      xminf =  1.e+38
      hflgu  = fa
      h3flgu = fa
      call formfe(b,dr,dr,dr,fa,fa,fa,fa,14,1,numel,1)
c.... Plot str1-palette in new colors
      call setcol(2,0)
      call pltconf(idev,xmaxf,xminf,k1,2)
c.... plot stresses from 1 point integration
      flfp = fa
      call formfe(b,dr,dr,dr,fa,fa,fa,fa,14,1,numel,1)
c.... Set colors back
      call setcol(1,0)
      go to 200
c.... plot material forces
c.... [jint,v1,v2,v3]plot material forces (separately!)
c                    (v1 =  i:      force(i) on undeformed mesh)
c                    (v1 = -i:      force(i) on   deformed mesh)
c                    (v1 =  0:      all forces (1-ndm) on undeformed mesh)
c                    (v1 = (ndf+1): all forces (1-ndm) on undeformed mesh)
c                    (v1 =-(ndf+1): all forces (1-ndm) on   deformed mesh)
c                    (v1 = (ndf+2): all moments(4-5/6) on undeformed mesh)
c                    (v1 =-(ndf+2): all forces (4-5/6) on   deformed mesh)
c
93    continue
      if(iadd.gt.0)goto 931
      c = 0.0
      k1 = abs(int(ct(1)))
      k1 = min(k1,ndf+2)
      if(ct(1).lt.0.0d0.or.defo) c = cs
      k2 = ct(2) ! tip
      rk3 = dabs(ct(3))
      k4 = ct(3)
      if(ct(3).eq.0) rk3 = 1.0d0
      if(fdrp)dbgtxt='PPLOTF:[reac]gen. array: numnp*ndf*ipr(fdrp)'
      call ralloc(drp,numnp*ndf,'pplotf-help-field',fdrp)
      call pzero(dr,neq)
      hflgu  = fa
      h3flgu = fa
c...  calculate material forces
      call formfe(b,drp,drp,drp,fa,tr,fa,tr,16,1,numel,1)
931   call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plopen
      if(k1.le.3)     call pppcol(2) ! force
      if(k1.gt.3)     call pppcol(6) ! moment
      if(k1.eq.ndf+1) call pppcol(2) ! forces
      if(k4.lt.0) then
c.....plot only nodes with material forces
        call pltload(dr,drp,id,ndm,ndf,numnp,2)
      else
        call pltfor2(dr,drp,bang,id,ndm,ndf,numnp,k2,rk3,
     +                tra,vr,1,k1,ix,nen1,3)
      end if
      go to 150
c
c.... plot nodes on intersections
c.... [isec,v1]    (v1 < 0 on deformed mesh)
c                  (v2.gt.0 plot nodes numbers)
c
94    c = 0.0
      k2=ct(2)
      if(ct(1).lt.0.0d0.or.defo) c = cs
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      call plopen
      call pppcol(1)
      call pltnod2(dr,ndm,numnp,k2)
      go to 200
c
c.... 95:[traj,v1,v2,v3] plot stress trajectories
c           v1 =    no of trajectory
c....       v1 < 0  plot on deformed mesh
c....       v2      scaling factor for length
c....       v3      0,2=2D,3=3D
95    c   = ct(1)
      if(c.lt.0.0d0.or.defo) c = cs
      k1 = abs(c)
      rk2 = ct(2)
      if(rk2.eq.0.d0) rk2 = 1.0d0
      k3  = ct(3)
      call plopen
      if(plfl) then
        dbgtxt = 'PPLOTF: [traj] gen. array: np,numnp*npstr*ipr(plfl)'
        call ralloc(strea,numnp*npstr,'TRAJ',plfl)
      end if
      !k4 = np + numnp*ipr   replaced by strea(1+numnp)
      if(.not.fl(11)) then
        strea  = 0.d0
        hflgu  = fa
        h3flgu = fa
        call formfe(b,dr,dr,dr,fa,fa,fa,fa,8,1,numel,1)
        call pltstr(strea,strea(1+numnp),numnp)
      end if
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plttraj(dr,strea(1+numnp),ndm,numnp,tra,vr,rk2,k3,k1)
      fl(11) = tr
      go to 150
c
c
c96/97[END, QUIT] dummy only for wd array
96    go to 200
97    go to 200
c
c.... magnify by mouse-click for mesh/defm (only idev=4)
c.... [magn,v1]
c                  v1 > 0 mesh, v1 < 0 deformed mesh
c                 |v1|  incremental value (def.=0.05)
98    if(idev.ne.4) goto 200
      xmagn = dabs(ct(1))
      if(xmagn.eq.0.d0) xmagn = 0.05d0
      c = 0.0d0
      if(ct(1).lt.0.0d0.or.defo)  c = cs
      write(*,*) '+ = left B., - = right B., end = both B. of mouse'
      istat  = 0
980   iclear = 0
      call plopen
c.... magnify
      if(istat.eq.1)  then
c....   increase
        xfac(1) = xfac(1) + xmagn
        xfac(2) = xfac(2) + xmagn
        xfac(3) = xfac(3) + xmagn
      else if(istat.eq.2) then
c....   decrease
        xfac(1) = xfac(1) - xmagn
        xfac(2) = xfac(2) - xmagn
        xfac(3) = xfac(3) - xmagn
      end if
c...  calculate mesh/defm
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
c      call frame(dr,ndm,numnp,1,clip) ! not scale
c...  plot mesh/defm
      call pppcol(4)
      if(ct(1).lt.0.0d0.or.defo) call pppcol(5)
      call pline(dr,ie,mnix,mea,numnp,nne,ndm,nxd,nxn,nie,0.d0,tr)
c...  decide
      call plmousec(istat)
      if(istat.lt.3) goto 980
      goto 200
c
c.99  [pdis] prescribed displacements --> 10
c
c.... plot interlaminar stresses over thickness at center of element n1
c.... [ints,n1,n2,n3]
c                  n1 element number
c                  n2 number of stress
c                  n3 number of output intervals -> n3+1 output-points
100   intn = ct(1)
      ispv = ct(2)
      intv = ct(3)
      if(intv.eq.0) intv=10
      ipv  = 2
c.... calculate stresses in strea (special element formulation necessary!)
      if(plfl) then
        dbgtxt = 'PPLOTF: [ints] gen. array: np,numnp*npstr*ipr(plfl)'
        call ralloc(strea,numnp*npstr,'INTS',plfl)
      end if
      strea  = 0.d0
      hflgu  = fa
      h3flgu = fa
      call formfe(b,dr,dr,dr,fa,fa,fa,fa,4,intn,intn,1)
c.... plot stresses
      call plotxyl1(strea,npv,ispv,s0(1),s0(2))
      intn=0
      goto 200

c.... imperfect displacement contours, similar to DISP
c.... [dimp,v1,v2,v3,v4]   plot filled /contour displacements
c                    v1 -    number of displ., neg. on deform. mesh
c                    v2 >  0 number of displ. values for plot
c                    v2 <= 0 filled plot
c                    v3.ne.0 contour plot without numbers
c                    v4 -  plot mesh in color v4
101   continue
      if(iadd.gt.0)goto 1010
      k2 = ct(2)
      n = abs(k2)
      n = max(1,n)
      i = max(1,abs(int(ct(1))))
      i = min(i,ndf)
      ipb = ct(3)
      k5  = ct4
      if(fdrp)dbgtxt='PPLOTF:[dimp]gen. array: numnp*ndf*ipr(fdrp)'
      call ralloc(drp,numnp*ndf,'pplotf-help-field',fdrp)
      call pmove(b,drp,numnp*ndf)           ! v = u
      if(mimp.ne.1)
     + call daxpty1(numnp*ndf,1.d0,uimp,ndm,drp,ndf,numnp)! v = u+u_imp
      call panglb(drp,bang,numnp,ndf)
      c = 0.0
      if(ct(1).lt.0.0d0.or.defo) c = cs
1010  continue
      if(icleas.eq.0) then
        if(ipola.ne.0) call dispola(x,drp,ndm,ndf,numnp,ipola)
        call rprint1(drp,ix,x,ndm,numnp,ndf,n,nen1,i,idev)
      end if
      call pdefm(x,b,c,bang,ndm,ndf,numnp,dr,xfac)
      call plxtrn(dr,tra,vr,ndm,numnp)
      if(k2.le.0) then
        n = -i
      end if
      call pltcon(dr,id,ie,mnix,drp,idis,nie,ndm,ndf,nxd,nxn,nne,
     1            i,n,2,1,k5,cinv)
      icleas = 1
      go to 150

c.... [fill] --> new name MATE
102   goto 18

c.... add quadrants in case of symmetry
150   continue

      if(nsym.eq.0) goto 200
      iadd = iadd + 1
      call psymm(x,ndm,numnp,isym1,ndm,iadd,1)
      call psymm(b,ndf,numnp,isym1,ndm,iadd,1)
      if(iadd.eq.nsym) then
        iadd =  0
        go to 200
      end if
      call plclos
      go to 301  ! for symm repeat plot
c.... close actual plot macro
200   continue
      icleas = 0
      call plclos
      if(ijump.ne.0) go to 260
      call plclos
       
      return
c.... error message
265    call  errclr ('PPLOTF')
      go to 260
c.... formats
1000  format(a4,11x,3f15.0,15(a1))
1001  format(a4,11x,a4)
c2001  format(' ** WARNING ** no match on a plot,',a4,' request.')
cww2002  format(5x,'Output of plots will be to screen.'/)
2003  format(5x,'Output of plots will be to file: ',a,/)
cww2004  format(' Input PLOT instructions, use "help" for a list, ',
cww     1 'exit with "end".'/ '     Plot ',i3,'> ',$)
2004  format('Plot>',$)
2005  format(' Input No. of Lines/Colors, Min. Value, Max. Value for ',
     +       ' contur/filled Plots'/3x,'>',$)
2011  format(' Displacement output set to cartesian directions ')
2012  format(' Displacement output set to polar (',i3,' ) directions ')
cww2010  format(1x,'Plot >',$)
cww2013  format(' Input PLOT instructions, written on PS/CGM-file, ',
cww     1 'close with "PLOF".'/ '     Plof ',i3,'> ',$)
2013  format('Plof>',$)
cww2014  format(' Input PLOT instructions, written on PS/HPGL-file, ',
cww     1 'close with "PRIN".'/ '     Prin ',i3,'> ',$)
2014  format('Prin>',$)
2015  format(' Zoom between (x,y,z)_1 (x,y,z)_2:',$)
2016  format('   Time for caculating nodal stresses   ',25x,'t=',0pf9.4)
2017  format('   Time for plot of nodal stresses      ',25x,'t=',0pf9.4)
2018  format('   Time for plotting DISP               ',37x,'t=',0pf9.4)
      end
c
      block data syminpt

      end
