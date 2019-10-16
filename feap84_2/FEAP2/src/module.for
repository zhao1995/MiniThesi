c...  allocation data
      module allocdata
c-----------------------------------------------------------------------
c.... Purpose: allocation of integer array
c     allocname  - name of array
c     asizes     - size of array 
c     isalloc    - store true/false
c     at present - 500     
c
c-----------------------------------------------------------------------
      integer, parameter                :: maxnum = 500
      character(30), dimension(maxnum)  :: allocname
      integer*8, dimension(maxnum)      :: memloc
      real*8, dimension(maxnum)         :: asizes
      logical, dimension(maxnum)        :: isalloc
      end module


      module adap
c-----------------------------------------------------------------------
c.... Purpose: adapativity arrays for REME new by ww 
c-----------------------------------------------------------------------
cww      integer, allocatable, dimension(:) :: 
cww  +           adaike, adaikz, adaiael, adaianp,
cww  +   adaiek0,adaike0,adaikz0,adaiael0,adaianp0
cww  +        ,adamikno,adanewbc
cww   real*8, allocatable, dimension(:) :: adax0,adaerron0 

      real*8, allocatable, dimension(:) :: adanx,adanf,adanu,adaerrot
      logical tflb
      end module

      module arcext
c-----------------------------------------------------------------------
c.... Purpose: arcext.h
c-----------------------------------------------------------------------
      real*8  exeps,xmu
      integer kex,kdig
      logical kflg
      end module

      module arcl
c-----------------------------------------------------------------------
c.... Purpose: arcl.h
c-----------------------------------------------------------------------
      real*8  rlnew,timold,ds0,alfa0,c0,cs1,cs2,r,det0,detc,xn
      integer mu1,mu2,kflag,ite,ndis,ndamp,nodis,nddis
      logical arcf,refl
      real*8, allocatable, dimension(:) :: arclm1,arclm2
      end module

      module augdat
c-----------------------------------------------------------------------
c.... Purpose: augdat.h
c-----------------------------------------------------------------------
      real*8  augf,cplus
      end module

      module aunr
c-----------------------------------------------------------------------
c.... Purpose: aunr.h
c-----------------------------------------------------------------------
      logical anr
      integer nanr
      real*8, allocatable, dimension(:) :: adr
      end module

      module back
c-----------------------------------------------------------------------
c.... Purpose: back.h
c-----------------------------------------------------------------------
      logical backstep
      end module

      module back1
c-----------------------------------------------------------------------
c.... Purpose: back1.h
c     ustore(nneq) - save u(t) for BACK 
c-----------------------------------------------------------------------
      real*8, allocatable, dimension(:) :: ustore 
      end module

      module bdata
c-----------------------------------------------------------------------
c.... Purpose: bdata.h
c     head - title used in element
c-----------------------------------------------------------------------
      character*4 o,head(20)
      end module

      module bisw
c-----------------------------------------------------------------------
c.... Purpose: bisw.h
c-----------------------------------------------------------------------
      integer iswb
      end module

      module boun
c-----------------------------------------------------------------------
c.... Purpose: boun.h
c-----------------------------------------------------------------------
      real*8 x1,x2,x3,dx1,trb(3,3),vrb(3)
      end module

      module cdat1
c-----------------------------------------------------------------------
c.... Purpose: cdat1.h
c-----------------------------------------------------------------------
      integer ndd,nie
      end module

      module cdata
c-----------------------------------------------------------------------
c.... Purpose: cdata.h
c-----------------------------------------------------------------------
      integer numnp,numel,nummat,nen,neq,ipr
      end module

      module codat
c-----------------------------------------------------------------------
c.... Purpose: codat.h
c-----------------------------------------------------------------------
      logical coflg
      end module

      module colmap
c-----------------------------------------------------------------------
c.... Purpose: colmap.h
c-----------------------------------------------------------------------
      integer*2 icolmap (32)
      end module

      module comfil
c-----------------------------------------------------------------------
c.... Purpose: comfil.h
c     Filenames
c-----------------------------------------------------------------------
      character*229 finp,fout,fres,fsav,fplt
      end module

      module contval
c-----------------------------------------------------------------------
c.... Purpose: contval.h
c-----------------------------------------------------------------------
      real*8  contv(3)
      integer icv
      end module

      module conv
c-----------------------------------------------------------------------
c.... Purpose: conv.h
c-----------------------------------------------------------------------
      real*8  ener(5),propold
      integer iconv,iconv1,nneg,ineg
      end module

      module conval
c-----------------------------------------------------------------------
c.... Purpose: conval.h
c        vvv(i, 1-26) - Character a-z 
c        vvv(i,27-36) - Character 0-9 
c        vvv(i,    0) - Character ' ' 
c        vvv(1-26,k)  - Character a-z
c        www(26)      - Upper case letters with values assigned
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      real*8 vvv(26,0:36),www(26)
      end module

      module ddata
c-----------------------------------------------------------------------
c.... Purpose: ddata.h
c-----------------------------------------------------------------------
      real*8  theta(4)
      integer nrk,nrc,nrm,nrt,nop
      end module

      module damp1
c-----------------------------------------------------------------------
c.... Purpose: damp1.h
c-----------------------------------------------------------------------
      integer ncu,ncll
      logical flgc,flgda
      end module

      module debugs
c-----------------------------------------------------------------------
c.... Purpose: debugs.h
c-----------------------------------------------------------------------
      integer   debug
      character dbgtxt*80
      end module

      module dii
c-----------------------------------------------------------------------
c.... Purpose: dii.h
c     ndii(i,50) - i=1:negative,i02:zero,i=3:very smalldiagonal elements 
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      integer ndii(3,50),ii(3)
      end module

      module dirdat
c-----------------------------------------------------------------------
c.... Purpose: dirdat.h
c     basea(10,knode,2) - data for basevectors  replaces m(mdir)
c     knode             - nen*numel: at element nodes
c     knode             - numnp:     at global  nodes
c     xdir              - 1=typ, 2-4 add. data
c     ldir              - 0: not used, 1: used
c     mdir              - adress in m-array, later to delete
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      real*8  xdir(4)
      integer ldir,mdir,knode
      real*8, allocatable, dimension(:) :: basea
      end module

      module dspos
c-----------------------------------------------------------------------
c.... Purpose: dspos.h
c-----------------------------------------------------------------------
      real*8 dscor(3,2),dsdcor(3),dsdcor2,tolc
      end module

      module dtauto
c-----------------------------------------------------------------------
c.... Purpose: dtauto.h
c-----------------------------------------------------------------------
      real*8  htol,dtdo,dtup1,dtup2,dtmax
      integer itimax,iaback,icstop
      end module

      module dyndat
c-----------------------------------------------------------------------
c.... Purpose: dyndat.h
c     dynrea(neq) - dyn. reactions/help vector  replaces m(ndyn)
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      logical fldyn
      real*8, allocatable, dimension(:) :: dynrea 
      end module

      module edgdat
c-----------------------------------------------------------------------
c.... Purpose: edgdat.h
c-----------------------------------------------------------------------
      integer ne4,nde,ned,nedg(100),iedg(12,100)  ! ne1, ne2, ne3, ne5 removed
      integer, allocatable, dimension(:) :: edge1,edge2,edge3,edge4
      end module

      module eig1
c-----------------------------------------------------------------------
c.... Purpose: eig1.h
c-----------------------------------------------------------------------
      real*8  omegax,epseig
      integer nte,nzykel,nxeig1,nxeig2
      logical eigflg
      real*8, allocatable, dimension(:) :: eigk1,eigk2
      end module

      module eldata
c-----------------------------------------------------------------------
c.... Purpose: eldata.h
c-----------------------------------------------------------------------
      real*8  dm
      integer n,ma,mct,iel,nel
!$OMP THREADPRIVATE (dm,n,ma,mct,iel,nel)
      end module

      module endata
c-----------------------------------------------------------------------
c.... Purpose: endata.h
c-----------------------------------------------------------------------
      real*8  aengy,aold,rnorm
      logical rfl
      end module

      module epsdh
c-----------------------------------------------------------------------
c.... Purpose: epsdh.h
c-----------------------------------------------------------------------
      real*8  factx,facty,skfy,skfz,dvp,epsd(15),dxyz(3)
      integer iswm,nss,natyp,iepsd
      end module

      module epsd1
c-----------------------------------------------------------------------
c.... Purpose: epsdh1.h
c-----------------------------------------------------------------------
      real*8   resife2
      integer  it1,it2
      end module

      module epspu
c-----------------------------------------------------------------------
c.... Purpose: data for material changes in elmt35
c              calculate epsp(8) and store into epspg
c     epspg(numel,ngp,8) - global array
c
c-----------------------------------------------------------------------
      real*8  dmat1(8,8),dmat2(8,8),epsp(8),eps1(8),sig1(8)
c      real*8, allocatable, dimension(:,:,:) :: epspg
      real*8, allocatable, dimension(:) :: epspg
      logical fepspg
      end module

      module errblk
c-----------------------------------------------------------------------
c.... Purpose: errblk.h
c-----------------------------------------------------------------------
      integer iblk
      end module

      module errchk
c-----------------------------------------------------------------------
c.... Purpose: errchk.h
c-----------------------------------------------------------------------
      logical         errck
      end module

      module errin1
c-----------------------------------------------------------------------
c.... Purpose: errin1.h
c     u_om(3)  - energy element from stresses
c     e_om(3)  - energy element from stress differences=error
c     e_bar(3) - error element
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      real*8 u_om(3),e_om(3),e_bar(3),eval
      end module 

      module errin2
c-----------------------------------------------------------------------
c.... Purpose: errin2.h
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      integer iet(5),ioerr,nerr
      integer, parameter :: numerr = 3
c     numerr - 3 number of error indicators
c     nerr   - adress in m-array, later to delete 
      end module

      module errin3
c-----------------------------------------------------------------------
c.... Purpose: array e_ome
c     e_ome(numel,3) - error for elements
c     new ww KIT 11/14
c-----------------------------------------------------------------------
      real*8, allocatable, dimension(:) :: e_ome
      logical tflb1
      end module

      module errnam
c-----------------------------------------------------------------------
c.... Purpose: errnam.h
c     e_name(3) - error name
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      character*15 e_name(3)
      end module

      module evdata
c-----------------------------------------------------------------------
c.... Purpose: evdata.h
c     type of EV-problem
c-----------------------------------------------------------------------
      integer imtyp,ibuck
      end module

      module ext1
c-----------------------------------------------------------------------
c.... Purpose: ext1.h
c-----------------------------------------------------------------------
      integer mext
      real*8, allocatable, dimension(:) ::
     +        extkh,extkc,extkd,extke,extkdh,extkz1,extkz2
      end module

      module ext2
c-----------------------------------------------------------------------
c.... Purpose: ext2.h
c-----------------------------------------------------------------------
      real*8  eta
      integer nmode
      logical extflg,eflg
      end module

      module fdata
c-----------------------------------------------------------------------
c.... Purpose: fdata.h
c     fl - general flags
c-----------------------------------------------------------------------
      logical fl(12),pfr
      end module

      module fe2dat
c-----------------------------------------------------------------------
c.... Purpose: fe2dat.h
c-----------------------------------------------------------------------
      logical flgfe2
      real*8, allocatable, dimension(:) ::
     +        mfe2a,mfe2g1,mfe2g2,mfe2g3,mfe2ix,mfe2ht,mfe2tau
      end module

      module fe2mat
c-----------------------------------------------------------------------
c.... Purpose: fe2mat.h
c-----------------------------------------------------------------------
      integer matfe2
!$OMP THREADPRIVATE (matfe2)
      end module

      module fe2tran
c-----------------------------------------------------------------------
c.... Purpose: fe2tran.h
c     transfer arrays-10 entries at present
c     fcisi - input
c     fciso - output
c     fcisr - restart
c     fciss - save
c     fcisf - forward
c     fcisb - backward
c     fcish - save(hflg=false)
c-----------------------------------------------------------------------
      integer*4        irtyp,    irecl,    icgp
      character*229    fcisi(10),fciso(10),fcisr(10),fciss(10),
     +                 fcisf(10),fcisb(10),fcish(10),fcis,fresg
      end module

      module feapmpi
c-----------------------------------------------------------------------
c.... Purpose: feapmpi.h
c-----------------------------------------------------------------------
      logical parform
      end module

      module feapprog
c-----------------------------------------------------------------------
c.... Purpose: feapprog.h
c     path to different programs
c-----------------------------------------------------------------------
      character*229 feapsal,feapint,feaped,fnege,fproce,editor,
     +              psview,fpath,fintnet,fadobe,fgmesh,fcylt
      end module

      module fileno
c-----------------------------------------------------------------------
c.... Purpose: fileno.h
c-----------------------------------------------------------------------
      integer isno
      end module

      module fodata
c-----------------------------------------------------------------------
c.... Purpose: fodata.h
c-----------------------------------------------------------------------
      integer nf,nfs
      logical foflg,fostr,foout,fouflg
      real*8, allocatable, dimension(:) :: aifour1, aifour2
      end module

      module forccont
c-----------------------------------------------------------------------
c.... Purpose: forccont.h
c-----------------------------------------------------------------------
      real*8  vc(14)
      integer ipali(14),nc
      end module

      module fornam
c-----------------------------------------------------------------------
c.... Purpose: fornam.h
c     force names
c-----------------------------------------------------------------------
      character*15 forsus(11)
      end module

      module foutp
c-----------------------------------------------------------------------
c.... Purpose: foutp.h
c-----------------------------------------------------------------------
      character*229 foutpar(8)
      end module

      module hdata
c-----------------------------------------------------------------------
c.... Purpose: hdata.h
c     gh1( nhmax,numel) - global array history data save
c     gh2( nhmax,numel) - global array history data iterat 
c     gh3(nh3max,numel) - global array history data 
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      integer nh1,nh2,nh3
      real*8, allocatable, dimension(:) :: gh1,gh2,gh3
      end module

      module hdatam
c-----------------------------------------------------------------------
c.... Purpose: hdatam.h
c     eh1( nhmax) - local array history data save
c     eh2( nhmax) - local array history data iterat 
c     eh3(nh3max) - local array history data 
c      nhmax      - max length of h1/h2 in Problem
c     nh3max      - max length of    h3 in Problem
c     hflgu       - true: update h1/2-array false: no update h1/2-array
c     h3flgu      - true: update h3  -array false: no update h3  -array
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      integer nhmax,nh3max
      logical hflgu,h3flgu
      real*8, allocatable, dimension(:) :: eh1,eh2,eh3
      end module

      module hidden
c-----------------------------------------------------------------------
c.... Purpose: hidden.h
c-----------------------------------------------------------------------
      logical hidw1,hidw2,hidi  
      integer, allocatable, dimension(:) :: idis,mia
      real*8,  allocatable, dimension(:) :: dist
      end module

      module hpgl1
c-----------------------------------------------------------------------
c.... Purpose: hpgl1.h
c-----------------------------------------------------------------------
      real*8  clst(33),colps(3,32)
      integer lun,iprin
      data colps /
     +            0.00, 0.00, 0.00,   !  1 white->black
     +            1.00, 0.00, 0.00,   !  2 red
     +            0.00, 0.00, 1.00,   !  3 blue
     +            0.00, 1.00, 1.00,   !  4 cyan
     +            0.00, 1.00, 0.00,   !  5 green
     +            0.80, 0.20, 0.20,   !  6 brown
     +            1.00, 1.00, 0.00,   !  7 yellow
     +            1.00, 0.00, 1.00,   !  8 magenta
     +            0.75, 0.75, 0.75,   !  9 dark grey
     +            0.00, 0.00, 1.00,   ! 10 blue
     +            0.00, 1.00, 0.00,   ! 11 green
     +            0.00, 1.00, 1.00,   ! 12 cyan
     +            1.00, 0.00, 0.00,   ! 13 red
     +            1.00, 0.00, 1.00,   ! 14 magenta
     +            1.00, 1.00, 1.00,   ! 15 white
     +            0.00, 0.00, 0.00,   ! 16 black
     +            0.66, 0.00, 0.00,   ! 17 dark  red  17-30 for disp
     +            1.00, 0.00, 0.00,   ! 18 med.  red
     +            1.00, 0.25, 0.00,   ! 19 light red
     +            1.00, 0.50, 0.00,   ! 20 orange
     +            1.00, 0.75, 0.00,   ! 21 light orange
     +            1.00, 1.00, 0.00,   ! 22 yellow
     +            0.75, 1.00, 0.00,   ! 23 light yellow
     +            0.50, 1.00, 0.00,   ! 24 light green
     +            0.00, 1.00, 0.00,   ! 25 light cyan
     +            0.00, 1.00, 0.75,   ! 26 cyan
     +            0.00, 1.00, 1.00,   ! 27 dark  cyan
     +            0.00, 0.75, 1.00,   ! 28 light blue
     +            0.00, 0.50, 1.00,   ! 29 middle blue
     +            0.00, 0.00, 1.00,   ! 30 dark blue
     +            0.00, 0.00, 0.00,   ! 31 black
     +            0.00, 0.00, 0.00/   ! 32 black
      end module

      module hptext1
c-----------------------------------------------------------------------
c.... Purpose: hptext1.h
c-----------------------------------------------------------------------
      integer isize
      end module

      module idata
c-----------------------------------------------------------------------
c.... Purpose: idata.h
c     save macro data  limited to 200 macros!! 
c     vjs(3,201) -  ct Comnand options 
c      js(201)   - jct Command numbers to execute
c     ljs(201)   - lct Command options
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      integer     vjs(3,201),js(201),nn,nl1,nl2
      character*4 ljs(201)
      logical hadd
      end module idata

      module iimpd
c-----------------------------------------------------------------------
c.... Purpose: iimpd.h
c-----------------------------------------------------------------------
      real*8 viimpx,viimpy
      end module

      module implstep
c-----------------------------------------------------------------------
c.... Purpose: implstep.h
c-----------------------------------------------------------------------
      real*8 tm,tstart
      end module

      module inptc
c-----------------------------------------------------------------------
c.... Purpose: inptc.h
c-----------------------------------------------------------------------
      integer inptctrl,numnpic,numelic
      end module

      module iodata
c-----------------------------------------------------------------------
c.... Purpose: iodata.h
c-----------------------------------------------------------------------
      integer iop,ios,ird,iwd
      end module

      module iofile
c-----------------------------------------------------------------------
c.... Purpose: iofile.h
c-----------------------------------------------------------------------
      integer ior,iow
      end module

      module iosave
c-----------------------------------------------------------------------
c.... Purpose: iosave.h
c-----------------------------------------------------------------------
      integer lfile
      logical lread,lsave
      end module

      module isbfgs
c-----------------------------------------------------------------------
c.... Purpose: isbfgs.h
c     BFGS arrays
c     bfgsbo(nneq) 
c     bfgsbd(nneq) 
c     bfgsbv(neq) 
c     bfgsbw(neq) 
c     bfgsbt(nneq*3) 
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      real*8, allocatable, dimension(:) ::
     +        bfgsbo,bfgsbd,bfgsbv,bfgsbw,bfgsbt
      end module

      module isbfgs1
c-----------------------------------------------------------------------
c.... Purpose: isbfgs1.h
c     BFGS arrays
c     bfgsbs(neq*15) 
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      real*8, allocatable, dimension(:) :: bfgsbs
      end module

      module iscsr
c-----------------------------------------------------------------------
c.... Purpose: iscsr.h
c     csrja - pointer array csr storage
c     csrka - pointer array csr storage
c-----------------------------------------------------------------------
      integer  japt,kapt,isymcsr,ljacsr
      logical  ljapt,lkapt
      integer, allocatable, dimension(:) :: csrja, csrka
      end module

      module isecn
c-----------------------------------------------------------------------
c.... Purpose: isecn.h
c-----------------------------------------------------------------------
      integer isecno(1000),is
      end module

      module isgmr
c-----------------------------------------------------------------------
c.... Purpose: isgmr.h
c-----------------------------------------------------------------------
      real*8  tolgmr
      integer itergmr,im
      real*8,  allocatable, dimension(:) :: 
     +         rmgmrx,rmgmrss,rmgmrhh,rmgmrrs,rmgmrc,rmgmrs
      end module

      module isogeo
c-----------------------------------------------------------------------
c.... Purpose: isogeo.h
c-----------------------------------------------------------------------
      integer  NURnpatch,NURlenkv(2)
      logical  nurbs
      real*8   surface
      integer, allocatable, dimension(:) :: AInmpq, AIninc, AInien,
     +         AInipa,AInstre
      real*8,  allocatable, dimension(:) :: AInkv1, AInkv2
      end module

      module ispgmr
c-----------------------------------------------------------------------
c.... Purpose: ispgmr.h
c-----------------------------------------------------------------------
      real*8  tolpgmr
      integer iterpgmr,imp
      real*8,  allocatable, dimension(:) :: rmpgmrx,rmpgmrvv,rmpgmrrs,
     +                                      rmpgmrc,rmpgmrs
      end module

      module ispcg
c-----------------------------------------------------------------------
c.... Purpose: ispcg.h
c-----------------------------------------------------------------------
      real*8  tolcg
      integer itolcg,itercg

      real*8,  allocatable, dimension(:) :: 
     +                      amcgz,amcgzz,amcgr,amcgrr,amcgp,amcgpp,amcgx
      end module

      module isprec
c-----------------------------------------------------------------------
c.... Purpose: isprec.h
c-----------------------------------------------------------------------
      real*8   tolpc
      integer  ippc,lfil,ipcwk
      real*8,  allocatable, dimension(:) :: rmpcalu,rmpcwl,rmpcwu
      integer, allocatable, dimension(:) :: 
     +         impcjlu,impcju,impcjwl,impcjwu,impcjr,impclevs
      logical  lpreco
      end module

      module istat
c-----------------------------------------------------------------------
c.... Purpose: istat.h
c-----------------------------------------------------------------------
      character*6   wdi,wdo,wdr,wds
      character*229 newf
      end module

      module iwinio
c-----------------------------------------------------------------------
c.... Purpose: iwinio.h
c-----------------------------------------------------------------------
      integer*2 iwxs,iwys,iwxgs,iwygs
      integer   iwin(3)
      end module

      module jinteg
c-----------------------------------------------------------------------
c.... Purpose: jinteg.h
c-----------------------------------------------------------------------
      integer  njint
      real*8   rint
      logical  jflgu
      integer, allocatable, dimension(:) :: ajint
      end module

      module ldata
c-----------------------------------------------------------------------
c.... Purpose: ldata.h
c-----------------------------------------------------------------------
      integer l,lv,lvs(9),lve(9),li,lis(9),lie(9)
      end module

      module lplot1
c-----------------------------------------------------------------------
c.... Purpose: lplot1.h
c-----------------------------------------------------------------------
      logical lhpgl,lps
      integer ilsc
      end module

c...  add to make macl global accessible
      module maclg
c-----------------------------------------------------------------------
c.... Purpose: maclg
c-----------------------------------------------------------------------
      real*8, allocatable, dimension(:) :: macl
      end module

      module macprt
c-----------------------------------------------------------------------
c.... Purpose: macprt.h
c-----------------------------------------------------------------------
      logical prnt
      end module

      module mate
c-----------------------------------------------------------------------
c.... Purpose: material arrays
c     matenew(numel) - new mate array
c     mateold(numel) - old mate array
c     new ww KIT 04/15
c-----------------------------------------------------------------------
      integer, allocatable, dimension(:) :: matenew,mateorg
      logical flmat 
      end module

      module mdata
c-----------------------------------------------------------------------
c.... Purpose: mdata.h
c     edis(6*nst)       - Elmt-DISP-nn  
c     ecor(nen*ndm)     - Elmt-COOR-n0  
c     etem(nen)         - Elmt-TEMP-n1  
c     eeqn(nst)         - Elmt-EQNO-n2
c     epve(nst)         - Elmt-P   -n3
c     ekma(nst*nst)     - Elmt-K   -n4
c
c     nmat(nummat*nie)  - MAT-List- n5
c     edma(nummat*ndd)  - DMAT     -n6
c
c     psid(ndf*numnp)   - ID-Array -n7
c     coor(ndm*numnp)   - COOR     -n8
c     econ(nen1*numel)  - ELEM     -n9
c     gloa(numnp*ndf)   - LOAD     -n10
c     jdt12             - JD       -n12  
c     glo0(numnp*ndf)   - LOAD F0  -n13
c     gu  (3*numnp*ndf) - DISP     -n14
c     gtem(numnp)       - TEMP     -n11
cww   integer  nn,n0,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13! no longer needed
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      real*8, allocatable, dimension(:) :: edis,ecor,etem,epve,ekma,edma ! element vectors
      real*8, allocatable, dimension(:) :: coor,gloa,glo0,gu,gtem        ! global vectors
      integer,allocatable, dimension(:) :: eeqn,nmat,psid,econ
      integer,allocatable, dimension(:) :: jdt12                         ! JD-n12
      end module

      module mdat2
c-----------------------------------------------------------------------
c.... Purpose: mdat2.h
c     aang(nen)   - Elmt-ANGLE-n11a
c     bang(numnp) - Elmt-ANGLE-n11b
c     n11a,n11b   - not longer used
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      integer ia(2),itrot
      real*8, allocatable, dimension(:) :: aang, bang
      end module

      module mxasz
c-----------------------------------------------------------------------
c.... Purpose: mxasz.h
c-----------------------------------------------------------------------
      integer mxpro,mxneq,nren
      integer, allocatable, dimension(:) :: optin
      end module

      module ndata
c-----------------------------------------------------------------------
c.... Purpose: ndata.h
c     gstiff(gsize)    - stiffness matrix, size due to solver
c     dampm (gsize)    - damping matrix
c     massm ()         - lumped mass/iden matrix/geom. matrix
c     trans (nneq*nrt) - transient fields velo/acce, nrt due to algorithm
c     nal   - address lower part
c     nau   - address upper part
c     na    - address diagonal part 
c     nc    - address damping matrix
c     nl    - address 
c     nm    - address mass matrix 
c     nv    - address transient terms velo,acce 
c     nw    - address ?
c     gsize - size of K
c     most addresses are necessary
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      integer  na,nal,nau,nc,nl,nm,nv,nw,gsize
      real*8, allocatable, dimension(:) :: gstiff ! stiffness matrix
      real*8, allocatable, dimension(:) :: dampm  ! damping matrix
      real*8, allocatable, dimension(:) :: massm  ! lumped mass/iden matrix
      real*8, allocatable, dimension(:) :: trans  ! transient fields
      end module

      module ndatx
c-----------------------------------------------------------------------
c.... Purpose: ndatx.h
c-----------------------------------------------------------------------
      integer nx,nxl,nxu,nxll
      end module

      module nolink
c-----------------------------------------------------------------------
c.... Purpose: nolink.h
c-----------------------------------------------------------------------
      integer  nli1,nli2,nli3
      logical  flnk
      integer, allocatable, dimension(:) :: link1,link2,link3
      end module

      module pardi
c-----------------------------------------------------------------------
c.... Purpose: pardi.h
c     drpar(neq) - separate right hand side vector
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      real*8    dparm(64)
      integer*8 ipt(64)
      integer*4 iparm(64),nproc,maxfct,mnum,
     +          nrhs,mtype,iphase,msglvl,solver
      logical   lrhp
      real*8,   allocatable, dimension(:) :: drpar
      end module

      module pathn
c-----------------------------------------------------------------------
c.... Purpose: pathn.h
c-----------------------------------------------------------------------
      character*229 file(5)
      end module

      module pback
c-----------------------------------------------------------------------
c.... Purpose: pback.h
c-----------------------------------------------------------------------
      integer iback
      end module

      module pcent
c-----------------------------------------------------------------------
c.... Purpose: pcent.h
c-----------------------------------------------------------------------
      logical pfl
      end module

      module pclip
c-----------------------------------------------------------------------
c.... Purpose: pclip.h
c-----------------------------------------------------------------------
      real*8  xc(2,3)
      logical clip1
      end module

      module pcrit
c-----------------------------------------------------------------------
c.... Purpose: pcrit.h
c-----------------------------------------------------------------------
      integer nc1,nc2,nc3,nc4,ncs,icc
      logical clfl ! cww?? used also elsewhere??
      real*8, allocatable, dimension(:) :: crit
      end module

      module pdam
c-----------------------------------------------------------------------
c.... Purpose: pdam.h
c-----------------------------------------------------------------------
      integer iprd,ipld
      end module

      module pdata0
c-----------------------------------------------------------------------
c.... Purpose: pdata0.h
c-----------------------------------------------------------------------
      real*8 vmin(3),vmax(3),qq(3)
      end module

      module pdata1
c-----------------------------------------------------------------------
c.... Purpose: pdata1.h
c-----------------------------------------------------------------------
      logical iso
      real*8  scale,scaleg,s0(2),dx(2),sx(2),fact
      end module

      module pdata2
c-----------------------------------------------------------------------
c.... Purpose: pdata2.h
c-----------------------------------------------------------------------
      integer iclear,idev,iopl,imono,ipgl,icgm,ibor,ipola
      end module

      module pdata3
c-----------------------------------------------------------------------
c.... Purpose: pdata3.h
c-----------------------------------------------------------------------
      logical plfl
      integer npstr
      end module

      module pdata4
c-----------------------------------------------------------------------
c.... Purpose: pdata4.h
c-----------------------------------------------------------------------
      real*8  xmin(3),xmax(3),xzm(3,2)
      integer nzm1,nzm2,nzm3
      end module

      module pdata6
c-----------------------------------------------------------------------
c.... Purpose: pdata6.h
c     inord(100)    - element number to plot 
c     ipord(40,100) - nod connection how to plot mesh
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      integer inord(100),ipord(40,100)
      end module

      module pdata7
c-----------------------------------------------------------------------
c.... Purpose: pdata7.h
c-----------------------------------------------------------------------
      integer  ipb,iipma,ipla
      integer, allocatable, dimension(:) :: aipma
      end module

      module pdata8
c-----------------------------------------------------------------------
c.... Purpose: pdata8.h
c-----------------------------------------------------------------------
      integer iplot
      logical fopn
      end module

      module pdata10
c-----------------------------------------------------------------------
c.... Purpose: pdata10.h
c-----------------------------------------------------------------------
      real *8 cfp,xmaxf,xminf,scal
      integer nfp,klayf,ifor
      logical flfp
      end module

      module pdata11
c-----------------------------------------------------------------------
c.... Purpose: pdata11.h
c-----------------------------------------------------------------------
      real*8 deltax,deltay,xfac(3)
      end module

      module pdata12
c-----------------------------------------------------------------------
c.... Purpose: pdata12.h
c-----------------------------------------------------------------------
      integer ijump
      end module

      module pdatah
c-----------------------------------------------------------------------
c.... Purpose: pdatah.h
c-----------------------------------------------------------------------
      logical hide
      end module

      module pdatap
c-----------------------------------------------------------------------
c.... Purpose: pdatap.h
c-----------------------------------------------------------------------
      real*4  xp(400),xpp(200),ypp(200)
      integer ipan,istruc
      end module

      module pdatas
c-----------------------------------------------------------------------
c.... Purpose: pdatas.h
c-----------------------------------------------------------------------
      integer*4  isym1(8,3,3),isym2(8,3,3),icleas,iadd,its
      data isym1 /72*0/,icleas /0/, iadd /0/, its /2/                    !!!!! data originally in pplof line 2897
      data isym2 /9*0,1,14*0,
     2            0,1,1,0,0,1,1,0,0,0,1,1,0,0,1,1,0,0,0,0,1,1,1,1,
     3            17*0,1,6*0/
      end module


      module pftn77
c-----------------------------------------------------------------------
c.... Purpose: pftn77.h
c-----------------------------------------------------------------------
      integer*2 ixa,iya,icc,ixp(200),iyp(200)
      end module

      module pindex
c-----------------------------------------------------------------------
c.... Purpose: pindex.h
c-----------------------------------------------------------------------
      integer  nindex,numnpn
      logical  lindex
      integer, allocatable, dimension(:) :: apost
      end module

      module plinet
c-----------------------------------------------------------------------
c.... Purpose: plinet.h
c-----------------------------------------------------------------------
      integer ilinp
      end module

      module plodf
c-----------------------------------------------------------------------
c.... Purpose: plodf.h
c-----------------------------------------------------------------------
      integer*4 ipl(2,2),npldf,nstedf,mkflg,mmc,mmst,incmk,nploc
      real*8    pf
      end module

      module plodfa
c-----------------------------------------------------------------------
c.... Purpose: plodfa.h
c-----------------------------------------------------------------------
      real*8    reacc(10,6),react
      integer*4 noden(10),npldf1
      logical   flreac
      end module

      module plodfb
c-----------------------------------------------------------------------
c.... Purpose: plodf.h
c-----------------------------------------------------------------------
      integer*4 nplo
      integer*4 nincp,nincp0
      end module

      module plodfs
c-----------------------------------------------------------------------
c.... Purpose: plodfs.h
c-----------------------------------------------------------------------
      integer*4 nstri,nstrno
      end module

      module plodfu
c-----------------------------------------------------------------------
c.... Purpose: plodfu.h
c-----------------------------------------------------------------------
      real*8 valuse1,valuse2
      end module

      module plong
c-----------------------------------------------------------------------
c.... Purpose: plong.h
c-----------------------------------------------------------------------
      integer kmax
      end module

      module plotdrp
c-----------------------------------------------------------------------
c.... Purpose: hilfsfeld pplotf
c-----------------------------------------------------------------------
      logical fdrp
      real*8, allocatable, dimension(:) :: drp
      end module

      module plotter
c-----------------------------------------------------------------------
c.... Purpose: plotter.h
c-----------------------------------------------------------------------
      character*229 fpgl
      integer       nexte,ihpgl,ibps,icps
      real*8        xxp(200),yyp(200)
      end module

      module plslay
c-----------------------------------------------------------------------
c.... Purpose: plslay.h
c-----------------------------------------------------------------------
      integer klay,mlay,intn,intv,npv,ipv,ispv
      end module

      module pltran
c-----------------------------------------------------------------------
c.... Purpose: pltran.h
c-----------------------------------------------------------------------
      real*8 tra(3,3),vr(3),rotang(3)
      end module

      module pnodn
c-----------------------------------------------------------------------
c.... Purpose: pnodn.h
c-----------------------------------------------------------------------
      integer  nip,ntied,ixtie
      logical  flparv
      integer, allocatable, dimension(:) :: gtie,tecon
      real*8,  allocatable, dimension(:)  :: parvvel, parvacce
      end module

      module ppers
c-----------------------------------------------------------------------
c.... Purpose: ppers.h
c-----------------------------------------------------------------------
      real*8  eold(3),vold(3),e(3),enorm,fold,ff,tgold(3),tg(3),
     +        q(3,3),xlbda(3,3)
      integer kpers,kper1
      end module

      module printh
c-----------------------------------------------------------------------
c.... Purpose: print.h  new name
c-----------------------------------------------------------------------
      logical prt
      end module

      module prisdat
c-----------------------------------------------------------------------
c.... Purpose: prisdat.h
c-----------------------------------------------------------------------
      integer nptyp,nprip(8)
      end module

      module prlod
c-----------------------------------------------------------------------
c.... Purpose: prlod.h
c-----------------------------------------------------------------------
      real*8   prop,a(6,10)
      integer  iexp(10),ik(10),npld
      end module

      module proc
c-----------------------------------------------------------------------
c.... Purpose: proc.h
c-----------------------------------------------------------------------
      character*229 procpath
      end module

      module psethis
c-----------------------------------------------------------------------
c.... Purpose: psethis.h
c     arrays used by pseta, limited to 200,  finally 0!
c     update ww KIT 11/14
c
c-----------------------------------------------------------------------
      integer      ipset(200),lpset(200),kpset
      character*20 npset(200)
      end module

      module psize
c-----------------------------------------------------------------------
c.... Purpose: psize.h
c-----------------------------------------------------------------------
      integer  noff
      integer, parameter :: maxm = 3000000
      end module

      module ptext
c-----------------------------------------------------------------------
c.... Purpose: ptext.h
c-----------------------------------------------------------------------
      character*229 text1,text2
      end module

      module qload
c-----------------------------------------------------------------------
c.... Purpose: qload.h
c-----------------------------------------------------------------------
      real*8  propq
      integer mqloa
      real*8, allocatable, dimension(:) :: aqloa
      end module

      module rdata
c-----------------------------------------------------------------------
c.... Purpose: rdata.h
c-----------------------------------------------------------------------
      real*8  tol,rnmax,shift
      logical linear
      end module

      module rfeap
c-----------------------------------------------------------------------
c.... Purpose: rfeap.h
c-----------------------------------------------------------------------
      integer irfeap
      end module

      module rndata
c-----------------------------------------------------------------------
c.... Purpose: rndata.h
c-----------------------------------------------------------------------
      integer nmn,nmx
      end module

      module rpdata
c-----------------------------------------------------------------------
c.... Purpose: rpdata.h
c-----------------------------------------------------------------------
      real*8 rmn,rmx
      end module

      module rsum
c-----------------------------------------------------------------------
c.... Purpose: rsum.h
c-----------------------------------------------------------------------
      integer  ndfrs,nfs1
      integer, allocatable, dimension(:) :: irpt
      end module

      module sdata
c-----------------------------------------------------------------------
c.... Purpose: sdata.h
c-----------------------------------------------------------------------
      integer ndf,ndm,nen1,nst
      end module

      module sectio
c-----------------------------------------------------------------------
c.... Purpose: sectio.h
c-----------------------------------------------------------------------
      integer msec(13,10),numnps(2,10),numels(2,10),npstrs,igps
      end module

      module shmname
c-----------------------------------------------------------------------
c.... Purpose: shmname.h
c-----------------------------------------------------------------------
      character*25 shmbasename
      integer      lgp
      end module

      module sldata
c-----------------------------------------------------------------------
c.... Purpose: sldata.h
c-----------------------------------------------------------------------
      integer nums,iels(4,26),inods(26),ivals(26)
      end module

      module slid1
c-----------------------------------------------------------------------
c.... Purpose: slid1.h
c-----------------------------------------------------------------------
      integer, allocatable, dimension(:) :: cl00,cl01,cl02,cl04,cl05,
     +         cl06,cl07
      real*8,  allocatable,dimension(:)  :: cl03,cl08,cl09,cl10,cl11,
     +         cl12
      end module

      module slid2
c-----------------------------------------------------------------------
c.... Purpose: slid2.h
c-----------------------------------------------------------------------
      integer, allocatable, dimension(:) :: cl22,cl23
      real*8,  allocatable,dimension(:)  :: cl13,cl14,cl15,cl16,cl17,
     +                     cl18,cl19,cl20,cl21,cl24,cl25
      end module

      module slid3
c-----------------------------------------------------------------------
c.... Purpose: slid3.h
c-----------------------------------------------------------------------
      real*8   epsn
      integer  nsl,nsntl,nmntl
      logical  contfl
      integer, allocatable, dimension(:) :: cl26,cl27,cl28,cl29,cl31
      real*8,  allocatable,dimension(:)  :: cl30
      end module

      module slid4
c-----------------------------------------------------------------------
c.... Purpose: slid4.h
c-----------------------------------------------------------------------
      real*8, allocatable,dimension(:)   :: cl33,cl34,cl35
      end module

      module slid5
c-----------------------------------------------------------------------
c.... Purpose: slid5.h
c-----------------------------------------------------------------------
      integer naxi
      end module

      module slu
c-----------------------------------------------------------------------
c.... Purpose: slu.h
c       drslu(neq) - solution vector  
c     diagslu(neq) - diagonal entries
c-----------------------------------------------------------------------
      integer ifactors 
      real*8, allocatable,dimension(:) :: drslu,diagslu
      end module

      module smpak
c-----------------------------------------------------------------------
c.... Purpose: smpak.h
c     smsperm      -         permutation array SDRV 
c     smsperm1     - inverse permutation array SDRV
c     smsidr1(neq) - working array dasol2       
c     smsisp       - working array  m(isp)     SDRV  
c     smsrsp       - working array  m(irsp)    SDRV
c     update ww KIT 11/14
c-----------------------------------------------------------------------
cww   org   integer  iperm,iperm1,isp,irsp,nsp,idr1
      integer  nsp
      logical  lodrv
      integer, allocatable, dimension(:) :: smsperm,smsperm1,smsisp
      real*8,  allocatable, dimension(:) :: smsidr1,smsirsp
      end module


      module soltyp
c-----------------------------------------------------------------------
c.... Purpose: soltyp.h
c-----------------------------------------------------------------------
      real*8   ctis(6)
      integer  istyp
      end module

      module stepc
c-----------------------------------------------------------------------
c.... Purpose: stepc.h
c-----------------------------------------------------------------------
      real*8  sp,rm,cmax,unorm,cnorm,rold,cs1o,cs2o
      integer itd,iti
      logical arcfs
      end module

      module strnam
c-----------------------------------------------------------------------
c.... Purpose: strnam.h
c     strsus(26)          - stress names
c     strea((npstr*numnp) - stress values st and weights dt ersetzt m(np)
c     npstr=26+1 
c     strea contains
c     dt(numnp)           - weights
c     st(numnp,26)        - stresses
c     ipstv               - number of stresses
c     np                  - adress in m-array, later to delete     
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      integer      np,istv
      character*15 strsus(26)
      real*8, allocatable, dimension(:) :: strea  
      end module

      module subdt
c-----------------------------------------------------------------------
c.... Purpose: subdt.h
c     eigd(mfmax)     - Eigenvalues  from SUBS,LAN,FEAST,RSG     
c     eigv(mfmax*neq) - Eigenvectors from SUBS,LAN,FEAST,RSG
c     evi             - Eigenvalues  from EIGI     
c     eigi(neq)       - Eigenvector  from EIGI
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      real*8               evi
      real*8, allocatable, dimension(:) :: eigd,eigv,eigia
      integer              mf,mf2,mfmax   ! old: mv,md,mv,
      end module

      module sumdt
c-----------------------------------------------------------------------
c.... Purpose: sumdt.h
c-----------------------------------------------------------------------
      integer nsum
      logical flsum
      real*8, allocatable, dimension(:) :: summ ! ersetzt m(nsum)
      end module

      module tdata
c-----------------------------------------------------------------------
c.... Purpose: tdata.h
c-----------------------------------------------------------------------
      real*8 ttim,dt,dto,c1,c2,c3,c4,c5
      end module

      module timex
c-----------------------------------------------------------------------
c.... Purpose: timex.h
c-----------------------------------------------------------------------
      real*4  tma,tm
      integer jtime
      end module

      module tplomax
c-----------------------------------------------------------------------
c.... Purpose: tplomax.h
c-----------------------------------------------------------------------
      integer imaxx,imaxy
      real*8  xmint,xmaxt,ymint,ymaxt
      end module

      module uneig
c-----------------------------------------------------------------------
c.... Purpose: uneig.h
c-----------------------------------------------------------------------
      real*8,allocatable,dimension(:) :: eigmma,aeigr,aeigi,aeigv
      end module

      module vdata
c-----------------------------------------------------------------------
c.... Purpose: vdata.h
c-----------------------------------------------------------------------
      character*16 versn(3)
      end module

      module vpoint
c-----------------------------------------------------------------------
c.... Purpose: vpoint.h
c-----------------------------------------------------------------------
      real*8 vwpt(3),vwpto(3)
      end module

      module wincfg
c-----------------------------------------------------------------------
c.... Purpose: wincfg.h
c-----------------------------------------------------------------------
      real*8       fpx1,fpy1,fwx1,fwy1,fpx2,fpy2,fwx2,fwy2
      real*8       fpx3,fpy3,fwx3,fwy3
      character*40 title1,title2,title3
      end module

c...  working array
      module working
c-----------------------------------------------------------------------
c.... Purpose: 
c-----------------------------------------------------------------------
      real*8, allocatable, dimension(:) :: plo,dr
      end module

      module ximp
c-----------------------------------------------------------------------
c.... Purpose: ximp.h
c-----------------------------------------------------------------------
      integer  mimp
      logical  flimp
      real*8, allocatable, dimension(:) :: uimp
      end module

      module yltadr
c-----------------------------------------------------------------------
c.... Purpose: yltadr.h
c-----------------------------------------------------------------------
      integer         n00,n01,n02,n03,n04,n05,n06,n07,n08,n09,
     &                n10,n11,n12,n13,n14,n15,n16,n17,n18,n19,
     &                n20,n21,n22,n23,n24,n25,n26,n27,n28,nn
      end module

      module yltdata1
c-----------------------------------------------------------------------
c.... Purpose: yltdata1.h
c-----------------------------------------------------------------------
      real*8  alpha
      integer nkn,nth,nknr,nthr,ndx,idx(2),ny1,iter
      logical aff
      end module

      module yltdata2
c-----------------------------------------------------------------------
c.... Purpose: yltdata2.h
c-----------------------------------------------------------------------
      real*8 dlm(6,1000),angle(3)
      end module

      module yltdata3
c-----------------------------------------------------------------------
c.... Purpose: yltdata3.h
c-----------------------------------------------------------------------
      real*8   eps1,eps2,eps3,eps4,eps5,eps6
      integer  ihalt
      logical  halt
      end module

      module yltdata4
c-----------------------------------------------------------------------
c.... Purpose: yltdata4.h
c-----------------------------------------------------------------------
      real*8 alph
      end module

      module yltdata5
c-----------------------------------------------------------------------
c.... Purpose: yltdata5.h
c-----------------------------------------------------------------------
      integer npop,nprp,nrec,ncrt,imut,irec
      end module

      module yydata
c-----------------------------------------------------------------------
c.... Purpose: yydata.h
c-----------------------------------------------------------------------
      character yyy*80
      end module
