
c----------------------------------------------------------------------+
c                                                                      |
c     Special Subroutines for Yield-Line Element  jw071104             |
c                                                                      |
c     Part 1: Macro Instructions and Assemblation of Tableaux          |
c                                                                      |
c----------------------------------------------------------------------+


c=======================================================================

      subroutine pmacr6 (ul,xl,tl,ld,p,s,ie,d,id,x,ix,f0,f,t,jd,u,dr,
     &                  lct,ct,ndf,ndm,nen1,nst,nneq,prt,wd,j,jsh)

c----------------------------------------------------------------------+
c     Macro Instruction Subprogram 6  (jw-11/04)                       |
c     Macro Instructions for Yield-Line-Element using                  |
c     Simplex Optimization Algorithm                                   |
c----------------------------------------------------------------------+ 

c==== DECLARATION ======================================================
      USE cdat1
      USE cdata
      USE eldata
      USE fdata
      USE iofile
      USE mdat2
      USE psize
      USE tdata
c      USE sdata
      USE yltadr
      USE yltdata1
      USE yltdata3
      USE yltdata5

c---- mdat2      : rotation data
c     ia         :
c     itrot      :  
      implicit real*8 (a-h,o-z)
      implicit integer(i-n)

c---- Formal Parameters ------------------------------------------------

c---- cdata    : 
c     numnp    : number of mesh nodes
c     numel    : number of mesh elements
c     nummat   : number of material sets
c     nen      : maximum nodes/element
c     neq      : number active equations
c     ipr      : real variable precision  

c---- cdat1    :
c     ndd      :
c     nie      :

c---- dt      : search step length

c---- eldata     : element data
c     dm         : element proportional load
c     n          : current element number
c     ma         : current element material set
c     mct        : print counter
c     iel        : user element number
c     nel        : number nodes on current element  

c---- sdata      :
c     ndf        :
c     ndm        :
c     nen1       :
c     nst        :
             
c---- iofile     : input-output-file
c     ior        : input file unit number (<0 input keyboard)
c     iow        : output file unit number
      
c---- yltdata1 : data of ylt-element
c     nkn      : number of nodes
c     nth      : number of edges (theta) 
c     nknr     : reduced number of nodes (after condensation)
c     nthr     : reduced number of edges (after condensation)          
c     ndx      : number of derivations 
c     idx      : specify derivation: 1: node number, 2: 1=d/dx,2=d/dy 
c     ny1      : number of y_1 entries
c     alpha    : alpha history value
c     iter     : number of iteration    
c     aff      : affine distortion     
      
c---- yltdata3 : tolerances   
c     eps1     : 1st tolerance check
c     eps2     : 2nd tolerance check
c     eps3     : degereration criterion 
c     eps4     : alpha bracketing tolerance 
c     eps5     : step length for direct search 
c     eps6     : zero tolerance for sqp
c     halt     : iterations stopped
      
c---- yltdata5 : genetic data
c     npop     : number of individuals
c     nprp     : number of properties (=ndf)
c     nrec     : number of genes (=nkn)   
c     ncrt     : number of evaluation criterions
c     imut     : mutation switch
c     irec     : recombination switch  

      logical tfl 
c      character*4 wd(1),lct
      integer id(ndf,numnp), ix(6,numel)!,jd(*) , ld(*),ie(*)
      real*8  p(*),s(*), !ct(3),ul(*),xl(*),tl(*),
     &        d(ndd,*),u(ndf,numnp),x(ndm,numnp),f(*)!,dr(*),t(*),f0(*)
 
c      logical prt
c      character*4 lct(*)
c      integer id(*),jd(*),ix(*)
c      real*8  d(*),ct(3,*),x(*),f(*),f0(*),u(*),dr(*),plo(7,*)
  
c---- d(*)       : element data parameters
c---- id         : boundary conditions
c---- ix(nen)    : element global node numbers
c---- x          : node coordinates
c---- f          :
c---- f0         :
c---- jd         : 
c---- u          :
c---- dr         :
c---- lct        :
c---- ct         :
c---- ndf        :
c---- nen        :
c---- ndm        :
c---- nneq       :
c---- j          :
c---- prt        :
c---- plo        :

cc---- prt        : print option (nopr)
cc---- tfl        : error occurred in srt pseta
cc---- lct        : macro name
cc---- wd         : macro list
cc---- id         : boundary conditions
cc---- ix(nen)    : element global node numbers
cc---- jp         :
cc---- ul(ndf,nen): element nodal solution parameters
cc---- xl(ndm,nen): element nodal reference coordinates
cc---- tl(nen)    : element nodal temperature values
cc---- s(nst,nst) : element matrix
cc---- p(ndf,nen) : element vector
cc---- d(*)       : element data parameters
cc---- x          : node coordinates
cc---- f0         :
cc---- f          :
cc---- t          :
cc---- b          :
cc---- dr         :
cc---- ndf        : number unknowns per node (nodal dof)
cc---- ndm        : space dimension of mesh
cc---- nst        : size of element arrays (=ndf*nen)
cc---- isw        : task parameter


c---- fdata    :
c     fl(11)   : tecplot

      logical ltec
      
c---- halt       : iterations stopped     
c---- aff        : affine distortion done
c---- ltec       : tecplot-interface used

c---- Variables --------------------------------------------------------

c      integer lm(2,nth)
c      real*8 atb(1+nth,2*nkn+2*nth)
c      real*8 etb(nth,nkn),ftb(nkn)
c      real*8 ptb(1+nth),ctb(2*nkn+2*nth),mtb(2*nth)
c      real*8 vtb(2*nkn+2*nth)
c      real*8 btb(2*nth+2*nkn,2*nth+2*nkn)
c      real*8 th(nth)
      
cc      real*4 satb(1+nth,2*nkn+2*nth),svtb(38),sctb(1,2*nkn+2*nth)
cc      real*4 sptb(1+nth,1)
cc      real*4 sbtb(2*nth+2*nkn,2*nth+2*nkn)

c---- lm       : location matrix (nkn,nth)
c     atb      : table of constrains (th+th-w+w-,1+th) A
c     etb      : assembled element matrix E
c     ftb      : assembled element vector F 
c     ptb      : right hand side
c     ctb      : coefficients of objective function c
c     mtb      : assembled element moment vector m
c     vtb      : optimization variables (th+th-w+w-)
c---- th       : theta, rotations of yl's

c---- Variables --------------------------------------------------------

      real*8 xup,xlw,yup,ylw
      
c---- xup,ylw    : upper/lower coordinates to define step length       
 
      common m(maxm)
  
c==== PROGRAM ==========================================================

c---- Define SQP-Optimization Tolerances -------------------------------

      eps1=1.d-3             !1st tolerance check in yltupd
      eps2=1.d-5!1.d-3             !2nd tolerance check in yltupd
      eps3=1.d-5!1.d-2       !degeneration criterion, null-step, alpha
      eps4=6.d-3             !alpha-Bracketing-Interval
      eps5=1.d-2
      eps6=1.d-9 !27         !sqp-tolerance, small values of degeneration
      
      ltec=.true.

      
c---- Define Define Evolution Switches and Parameters ------------------

      maxpop=2000
      ncrt=5 !Evolution Crits: No,qu,qu_group,int_work,normalvector     

c---- Affine Distortion ------------------------------------------------

      if (.not.aff) then
        do i1=1,numnp
          x(1,i1)=d(7,ma)*x(1,i1)
          x(2,i1)=x(2,i1)
          x(3,i1)=0.d0
        enddo  
        aff=.true.
      endif
      
c---- After first call... ----------------------------------------------      

      if (ndx.gt.0) goto 10

c---- Count dof's for Derivation ---------------------------------------

      do i1=1,numnp
        do i2=1,2
          if (id(i2,i1).gt.0) ndx=ndx+1
        enddo
      enddo  

c---- Initialize Search Step Length ------------------------------------

      xup=x(1,1)
      xlw=x(1,1)
      yup=x(2,1)
      ylw=x(2,1)
      do i1=1,nkn
        if (x(1,i1).gt.xup) xup=x(1,i1)
        if (x(1,i1).lt.xlw) xlw=x(1,i1)
        if (x(2,i1).gt.yup) yup=x(2,i1)
        if (x(2,i1).lt.ylw) ylw=x(2,i1)
      enddo
      dt=min(dabs(xup-xlw),dabs(yup-ylw))
      dt=eps5*dt  

c---- Address Calculation (factor 2 for double precision) --------------

                                               ! start of ! previous array
c      n00=100000                               ! lm(*)    !
c      n01=n00+2*(2*nth)                        ! atb(*)   ! lm(2,nth)
c      n02=n01+2*((1+nth)*(2*nkn+2*nth))        ! etb(*)   ! atb(1+nth,2*nkn+2*nth) 
c      n03=n02+2*(nth*nkn)                      ! ftb      ! etb(nth,nkn)
c      n04=n03+2*(nkn)                          ! ptb      ! ftb(nkn)
c      n05=n04+2*(1+nth)                        ! ctb      ! ptb(1+nkn)
c      n06=n05+2*(2*nkn+2*nth)                  ! mtb      ! ctb(2*nkn+2*nth)
c      n07=n06+2*(2*nth)                        ! vtb      ! mtb(2*nth)
c      n08=n07+2*(2*nkn+2*nth)                  ! btb      ! vtb(2*nkn+2*nth)
c      n09=n08+2*((2*nth+2*nkn)*(2*nth+2*nkn))  ! th       ! btb(2*nth+2*nkn,2*nth+2*nkn)
c      n10=n09+2*(nth)                          ! a1       ! th(nth)
c      n11=n10+2*((nth+1)*(nth+1))              ! c1       ! a1(nth+1,nth+1)
c      n12=n11+2*(nth+1)                        ! y1       ! c1(nth+1)
c      n13=n12+2*(nth+1)                        ! ly1      ! y1(nth+1) 
c      n14=n13+1*(2*(nth+1))                    ! da1      ! ly1(2,nth+1) integer
c      n15=n14+2*((nth+1)*(nth+1)*ndx)          ! dc1      ! da1(nth+1,nth+1,ndx) 
c      n16=n15+2*((nth+1)*ndx)                  ! dqu(=dL) ! dc1(nth+1,ndx)
c      n17=n16+2*(ndx)                          ! lth      ! dqu(ndx) 
c      n18=n17+2*(3*nth)                        ! qu       ! lth(3,nth)
c      n19=n18+2*1                              ! alpha    ! qu(1)
c      n20=n19+2*1                              ! dtb      ! alpha(1)
c      n21=n20+2*(ndx*(2*nth+2*nkn))            ! vtbh     ! dtb(2*nth+2*nkn,ndx) orig.(ny1,ndx)
c      n22=n21+2*(2*nkn+2*nth)                  ! dquh     ! vtbh(2*nkn+2*nth)
c      n23=n22+2*(ndx)                          !          ! dquh(ndx)
c      
c      nn=n23
c      if (nn.gt.200000) then
c        write(*, '(a36)') '  WARNING: array overflow in pmacr6!'
c        write(*, '(i8,a36)') nn, ' entries required - stopping program'
c        stop  
c      endif
      
                                                              ! start of array
                                                              !
      call pseta(n00,(2*nth)                      ,ipr,tfl,'CYLT1-n00')   ! lm(2,nth)
      call pseta(n01,((1+nth)*(2*nkn+2*nth))      ,ipr,tfl,'CYLT1-n01')   ! atb(1+nth,2*nkn+2*nth) 
      call pseta(n02,(nth*nkn)                    ,ipr,tfl,'CYLT1-n02')   ! etb(nth,nkn)
      call pseta(n03,(nkn)                        ,ipr,tfl,'CYLT1-n03')   ! ftb(nkn)
      call pseta(n04,(1+nth)                      ,ipr,tfl,'CYLT1-n04')   ! ptb(1+nkn)
      call pseta(n05,(2*nkn+2*nth)                ,ipr,tfl,'CYLT1-n05')   ! ctb(2*nkn+2*nth)
      call pseta(n06,(2*nth)                      ,ipr,tfl,'CYLT1-n06')   ! mtb(2*nth)
      call pseta(n07,(2*nkn+2*nth)                ,ipr,tfl,'CYLT1-n07')   ! vtb(2*nkn+2*nth)
      call pseta(n08,((2*nth+2*nkn)*(2*nth+2*nkn)),ipr,tfl,'CYLT1-n08')   ! btb(2*nth+2*nkn,2*nth+2*nkn)
      call pseta(n09,(nth)                        ,ipr,tfl,'CYLT1-n09')   ! th(nth)
      call pseta(n10,((nth+1)*(nth+1))            ,ipr,tfl,'CYLT1-n10')   ! a1(nth+1,nth+1)
      call pseta(n11,(nth+1)                      ,ipr,tfl,'CYLT1-n11')   ! c1(nth+1)
      call pseta(n12,(nth+1)                      ,ipr,tfl,'CYLT1-n12')   ! y1(nth+1) 
      call pseta(n13,(2*(nth+1))                  ,1  ,tfl,'CYLT1-n13')   ! ly1(2,nth+1) integer
      call pseta(n14,((nth+1)*(nth+1)*ndx)        ,ipr,tfl,'CYLT1-n14')   ! da1(nth+1,nth+1,ndx) 
      call pseta(n15,((nth+1)*ndx)                ,ipr,tfl,'CYLT1-n15')   ! dc1(nth+1,ndx)
      call pseta(n16,(ndx)                        ,ipr,tfl,'CYLT1-n16')   ! dqu(ndx) 
      call pseta(n17,(3*nth)                      ,ipr,tfl,'CYLT1-n17')   ! lth(3,nth)
      call pseta(n18,1                            ,ipr,tfl,'CYLT1-n18')   ! qu(1)
      call pseta(n19,1                            ,ipr,tfl,'CYLT1-n19')   ! alpha(1)
      call pseta(n20,(ndx*(2*nth+2*nkn))          ,ipr,tfl,'CYLT1-n20')   ! dtb(2*nth+2*nkn,ndx) orig.(ny1,ndx)
      call pseta(n21,(2*nkn+2*nth)                ,ipr,tfl,'CYLT1-n21')   ! vtbh(2*nkn+2*nth)
      call pseta(n22,(ndx)                        ,ipr,tfl,'CYLT1-n22')   ! dquh(ndx)
      call pseta(n23,(ndm*nkn*maxpop)             ,ipr,tfl,'CYLT1-n23')   ! pop(ndm,nkn,npop)
c     call pseta(n24,1                            ,ipr,tfl,'CYLT1-n24')   ! war, reserved in srt yltsqp
      call pseta(nn,1                             ,ipr,tfl,'CYLT1-nn ')   

      nn=n23
      if (tfl) then
        write(*, '(a36)') '  WARNING: array overflow in pmacr6!'
        write(*, '(i24,a36)') nn, ' entries required - stopping program'
        stop  
      endif      

c---- Transfer to Process ----------------------------------------------

   10 continue
      goto (100,200,300,400,500,600), j

c==== [ytab] Assemble Optimization Tableau and execute Simplex =========

  100 continue

c---- Define preliminary parameters ------------------------------------

      iter=iter+1
      
c---- Assembly of Array + Element Check --------------------------------

c     call yltass(d(1,ma),u,x,ix,lm,id,f,s,p,atb,etb,
c     &           ftb,ptb,ctb,mtb,lth,30)
      call yltass(d(1,ma),u,x,ix,m(n00),id,f,s,p,m(n01),m(n02), 
     &            m(n03),m(n04),m(n05),m(n06),m(n17),30)                                                   
 
c---- Solve Simplex, calculate Ultimate Load --------------------------- 
     
c     call yltsqp(atb,btb,ctb,ptb,vtb,vtbh,isw)
      call yltsqp(m(n01),m(n08),m(n05),m(n04),m(n07),m(n21),1)
      
c      call yltevl(vtb,mtb,u,th,qu)
      call yltevl(m(n07),m(n06),u,m(n09),m(n18),1)             
                   
c---- Output Values ----------------------------------------------------

c      call yltprt(atb,ctb,vtb,mtb,lm,u,th,qu)
      call yltprt(m(n00),m(n09),m(n18),m(n16),1)  
      goto 7000
      
c==== [yang] Print Theta Values ========================================      

  200 continue    
c      call yltprt(lm,th,qu)
      call yltprt(m(n00),m(n09),m(n18),m(n16),2)

      goto 7000
      
c==== [ygra] Evaluate Gradient =========================================

  300 continue
    
c---- Criterion Already Met --------------------------------------------

      if (halt) goto 7000

c---- Error Control ----------------------------------------------------

      if (ny1.eq.0) write(*, '(a31)') 
     &                    '  * YLT Error: use [ytab] first'
      if (ny1.gt.(nth+1)) then
        write(*, '(a41)') '  * YLT Error: rank(A_1) < number of y_1!'
c        print*, ny1, (nth+1)   
        goto 7000
      endif  

c---- Assemble Arrays and perform Gradient Calculation -----------------

c       call yltgrd(d,u,x,ix,lm,id,f,s,p,ang,atb,etb,ftb,ptb,ctb,mtb,vtb,
c     &                  a1,c1,y1,da1,dc1,dqu,dquh,lth,dtb,isw)

       call yltgrd(d(1,ma),u,x,ix,m(n00),id,f,s,p,bang,m(n01),m(n02), 
     &             m(n03),m(n04),m(n05),m(n06),m(n07),
     &             m(n10),m(n11),m(n12),m(n14),m(n15),m(n16),m(n22),
     &             m(n17),m(n20),31)
     
c       if (lgrd) then
         call yltprt(m(n00),m(n09),m(n18),m(n16),3)
c       else
c         call ylttry(d(1,ma),u,x,ix,m(n00),id,f,s,p,m(n11b),m(n01),
c     &            m(n02),m(n03),m(n04),m(n05),m(n06),m(n07),m(n08),
c     &            m(n09),m(n16),m(n17),m(n18))
c       endif


c     BAUSTELLE: Gradientenberechnung. A: Spalten streichen,
c     Lösung. Einbau der Diskontinuitätsstelle irgendwo hier...???

      goto 7000           

c==== [ymsh] Renew Yield-Line Mesh with Gradient =======================
                                                                        
  400 continue                                                          

c---- Criterion Already Met --------------------------------------------

      if (halt) goto 7000

c     call yltupd(d(1,ma),u,x,ix,lm,id,f,s,p,ang,atb,etb,
c     &           ftb,ptb,ctb,mtb,vtb,btb,th,
c     &           a1,c1,y1,da1,dc1,dqu,
c     &           lth,qu,alpha,dtb,vtbh,dquh)
      
      call yltupd(d(1,ma),u,x,ix,m(n00),id,f,s,p,bang,m(n01),m(n02), 
     &            m(n03),m(n04),m(n05),m(n06),m(n07),m(n08),m(n09),
     &            m(n10),m(n11),m(n12),m(n14),m(n15),m(n16),
     &            m(n17),m(n18),m(n19),m(n20),m(n21),m(n22))

      goto 7000
      
c==== [yevo] Perform Genetic Algorithm =================================

  500 continue 
     
c---- Dimension Population ---------------------------------------------

      call gendim(id,ndf,nkn,npop,maxpop)
      
c---- Genetic Algorithm ------------------------------------------------

c      call yltevo(d,u,x,ix,lm,id,f,s,p,atb,etb,ftb,
c     &            ptb,ctb,mtb,vtb,lth,qu,pop,isw)

      call yltevo(d(1,ma),u,x,ix,m(n00),id,f,s,p,m(n01),m(n02),m(n03),
     &            m(n04),m(n05),m(n06),m(n07),m(n17),m(n18),m(n23),1)

      call yltprt(m(n00),m(n09),m(n18),m(n16),4)
      
      goto 7000
      
c==== [ytry] Direct Search =============================================

  600 continue

c      call ylttry(d,u,x,ix,lm,id,f,s,p,ang,atb,etb,ftb,ptb,ctb,
c     &            mtb,vtb,btb,th,dqu,lth,qu)
 
      call ylttry(d(1,ma),u,x,ix,m(n00),id,f,s,p,bang,m(n01),m(n02),
     &            m(n03),m(n04),m(n05),m(n06),m(n07),m(n08),m(n09),
     &            m(n16),m(n17),m(n18))
      
      goto 7000

c==== Restore Original Coordinates from Affine Distortion ==============

 7000 continue
      if (aff) then
        do i1=1,numnp
          x(1,i1)=(1.d0/d(7,ma))*x(1,i1)
          x(2,i1)=x(2,i1)
          x(3,i1)=0.d0
        enddo  
        aff=.false.
      endif

c==== Tecplot-Interface ================================================

      if (ltec.and.(j.eq.1)) then
c        fl(11)=.true.
        call ylttec(x,u,m(n18),ix)
      endif  
          
c==== Formats ==========================================================    

     
c==== End of Subroutine ================================================
      
c9999  continue
      return
      end    
      
 
c=======================================================================      
c=======================================================================      

      subroutine yltbou(d,x,ix,id,ang,ndf,ndm,prt,iang1,iang2)
      
c----------------------------------------------------------------------+
c     Perform Boundary Conditions of YLT-Element                       |
c----------------------------------------------------------------------+
 
c==== DECLARATION ======================================================

c---- Formal Parameters ------------------------------------------------

      USE cdat1
      USE cdata
      USE iofile
      USE yltdata1
      USE yltdata2
      implicit real*8 (a-h,o-z)
      implicit integer (i,n)
      
c---- yltdata1 : data of ylt-element
c     nkn      : number of nodes
c     nth      : number of edges (theta) 
c     nknr     : reduced number of nodes (after condensation)
c     nthr     : reduced number of edges (after condensation)          
c     ndx      : number of derivations 
c     idx      : specify derivation: 1: node number, 2: 1=d/dx,2=d/dy 
c     ny1      : number of y_1 entries
c     alpha    : alpha history value
c     iter     : number of iteration    
c     aff      : affine distortion     

c---- yltdata2 : data of ylt-element 
c     dlm      : edge nodes, m_pl (<0: free edge)
c     angle    : angle of 3 element nodes

c---- cdat1    :
c     ndd      :
c     nie      :

c---- cdata    : 
c     numnp    : number of mesh nodes
c     numel    : number of mesh elements
c     nummat   : number of material sets
c     nen      : maximum nodes/element
c     neq      : number active equations
c     ipr      : real variable precision  

c---- iofile     : input-output-file
c     ior        : input file unit number (<0 input keyboard)
c     iow        : output file unit number
    
      integer ix(6,numel),id(ndf,numnp)
      integer iang1,iang2

c---- ix(nen)    : element global node numbers 
c---- id         : boundary conditions  
c---- iang1    : angle of first dof  (not in use yet)
c---- iang2    : angle of second dof (not in use yet)
      
      real*8 d(*),x(ndm,numnp),ang(numnp)!ang(numnp),ang(ndm,*)

c---- d(*)     : element data parameters
c---- x        : coordinates
c---- alpha    : angle of line
c---- ang      : angle of node (in degree)


      logical prt
 
c---- prt      : print option

c---- Variables --------------------------------------------------------

      integer nove(numnp)

c---- nove      : node vector of supported nodes

      real*8 td(10),v1(3),v2(3),v3(3),a,len(3),tol
      
c---- td        : data parameters for input      
c---- v*        : vector      
c---- a         : auxiliary variable
c---- len       : length
c---- tol       : tolerance

c==== PROGRAM ========================================================== 
 
      call pzero(ang,numnp)
      tol=1.d-8 
      
c==== Read Boundary Conditions ========================================= 
 
      if(ior.lt.0) write(*,8110)
      
      i0=0
  100 continue
      call pzero(td,10)
      nove=0
      call dinput(td(1:9),9)  ! x1,y1,,x2,y2,,boun,mpu,mpo
      a=0.d0
      do i1=1,9
        a=a+dabs(td(i1))
      enddo  
      if (a.ge.0.d0.and.a.le.0.d0) goto 200       ! boun: 0=simp.supp. 1=free 

c---- Search For Nodes On Boundary Line And Calculate Angle ------------

      v1(1)=td(4)-td(1)
      v1(2)=td(5)-td(2)
      v1(3)=0.d0
      len(1)=dsqrt(v1(1)**2+v1(2)**2)
 
      if (v1(1).ge.0.d0.and.v1(1).le.0.d0) then
        alpha=-dacos(0.d0)
      else 
        alpha=datan(v1(2)/v1(1))
c        if (v1(1).lt.0.d0) alpha=-alpha
      endif
 
      do i1=1,numnp
        v1(1)=x(1,i1)-td(1)
        v1(2)=x(2,i1)-td(2)
        v1(3)=0.d0
        
        v2(1)=x(1,i1)-td(4)
        v2(2)=x(2,i1)-td(5)
        v2(3)=0.d0
        
        len(2)=dsqrt(v1(1)**2+v1(2)**2)
        len(3)=dsqrt(v2(1)**2+v2(2)**2)

        call vnorm(v1,a)
        if (a.gt.0.d0.or.a.lt.0.d0) v1=v1/a
        call vnorm(v2,a)
        if (a.gt.0.d0.or.a.lt.0.d0) v2=v2/a
        call vcross(v1,v2,v3)
        call vnorm(v3,a)
 
        if ((a.lt.tol.and.(len(2)+len(3)-len(1)).lt.tol).or.
     &      ((x(1,i1).ge.td(1).and.x(1,i1).le.td(1)).and.
     &       (x(2,i1).ge.td(2).and.x(2,i1).le.td(2))).or.
     &      ((x(1,i1).ge.td(4).and.x(1,i1).le.td(4)).and.
     &       (x(2,i1).ge.td(5).and.x(2,i1).le.td(5)))) then

          nove(i1)=1
c          do i2=1,ndm
           ang(i1)=alpha*90.d0/dacos(0.d0)
c          enddo
        endif  
      enddo 
 
c      if(prt)              write(iow,8120) td(1),td(2),td(3),
c     &                                     td(4),td(5),td(6),
c     &                                    dint(td(7)),td(8),td(9),-alpha
      if (prt.and.ior.lt.0) write(*  ,8120) td(1),td(2),td(3),
     &                                      td(4),td(5),td(6),  
     &                                dint(td(7)),td(8),td(9)   !BAUSTELLE: nopr-Option
                                           
c---- Search For Boundary Lines Between Nodes --------------------------

      do i2=1,numnp
        if (nove(i2).eq.1) then
          do i3=1,numnp
            if (i3.eq.i2) cycle
            if (nove(i3).eq.1) then
              do i1=1,numel
                do ia=1,3
                  ib=ia+1
                  ic=ia+2
                  if (ib.gt.3) ib=ib-3
                  if (ic.gt.3) ic=ic-3
            
                  if (ix(ia,i1).eq.i2.and.ix(ib,i1).eq.i3) then 
                    do i4=1,numel
                      if ((dlm(1,i4).ge.i2.and.
     &                     dlm(1,i4).le.i2.and.
     &                     dlm(2,i4).ge.i3.and.
     &                     dlm(2,i4).le.i3).or.
     &                    (dlm(1,i4).ge.i3.and.
     &                     dlm(1,i4).le.i3.and.
     &                     dlm(2,i4).ge.i2.and.
     &                     dlm(2,i4).le.i2)) then
                        cycle
                      else                          
                        i0=i0+1
                        dlm(1,i0)=i2
                        dlm(2,i0)=i3
                        dlm(3,i0)=td(7)
                        dlm(4,i0)=td(8)
                        dlm(5,i0)=td(9)
                        dlm(6,i0)=alpha !*90.d0/dacos(0.d0)
c                        if (dabs(alpha).eq.dacos(0.d0)) then
c                          id(2,i2)=1
c                          id(2,i3)=1
c                        else
c                          id(1,i2)=1
c                          id(1,i3)=1
c                        endif    
                        if (td(7).ge.1.and.td(7).le.1) then  !simp.supp.
                          id(3,i2)=1 
                          id(3,i3)=1
                        endif  
                        exit  
                      endif
                    enddo   
                  endif
                enddo
              enddo
            endif
          enddo
        endif
      enddo            
    
c      td(1:9)=0.d0
      goto 100

c==== Set Axes =========================================================
      
  200 continue
      iang1=2
      iang2=1
      
c---- Count Supported Boundaries ---------------------------------------

      d(1)=0.d0
      do i1=1,1000
        if (dlm(1,i1).gt.0.or.dlm(1,i1).lt.0) then
          d(1)=d(1)+1.d0
        else
          exit
        endif
      enddo           
c      nth=(3*numel+d(1))/2

c---- Set Boundary Conditions ------------------------------------------
 
      do i1=1,numnp
        do i2=1,numnp-1
          do i3=i2+1,numnp
            if ((int(dlm(1,i2)).eq.i1.and.int(dlm(2,i3)).eq.i1).or.
     &          (int(dlm(2,i2)).eq.i1.and.int(dlm(1,i3)).eq.i1)) then
              if (dlm(6,i2).ge.dlm(6,i3).and.
     &            dlm(6,i2).le.dlm(6,i3)) then
                if (dabs(dlm(6,i2)).ge.dacos(0.d0).and.
     &              dabs(dlm(6,i2)).le.dacos(0.d0)) then
                  id(1,i1)=1
                  ang(i1)=0.d0
                else  
                  id(2,i1)=1
                endif
              else
                id(1,i1)=1
                id(2,i1)=1
              endif
            endif
          enddo
        enddo
      enddo        
                    
 
c      do i1=1,numnp-1
c        if (dlm(1,i1).eq.0) exit
c        do i3=i1+1,numnp
c          i2=i3
c          if (dlm(1,i2).eq.0) i2=1 
c          if (int(dlm(2,i1)).eq.int(dlm(1,i2))) then !.or.
cc     &        (int(dlm(1,i1)).eq.int(dlm(2,i2)))) then
c            if (dlm(6,i1).ne.dlm(6,i2)) then !both DOFs fixed
c              id(1,int(dlm(2,i1)))=1
c              id(2,int(dlm(2,i1)))=1
c            else !only first DOF fixed, eventually turned  
c              id(2,int(dlm(2,i1)))=1
c            endif
c          endif
c        enddo
c      enddo    

cc---- Kontrollausgabe
c
c      print*, 'angle'
c      do i1=1,numnp
c        print*, ang(i1)
c      enddo      
c      print*, 'nove'
c      do i1=1,numnp
c        print*, nove(i1)
c      enddo      
     
c---- Correct Angle If Both DOFs Hold Anyway --------------------------

      do i1=1,numnp
        if (id(1,i1).eq.1.and.id(2,i1).eq.1) ang(i1)=0.d0        
      enddo  
        
c==== Formats ==========================================================

8110  format(6x,'Y i e l d   L i n e   B o u n d a r i e s')
8120  format(/,' BOUN generation with YBOU on ',a8,/,
     &         ' From         Point ',3(1x,e12.5),/,
     &         ' To           Point ',3(1x,e12.5),/,
     &         ' Free Edge (1=yes)  ',   i5,/,
     &         ' m_pl lower surface ',e12.5,/,
     &         ' m_pl upper surface ',e12.5)!,/,
!     &         ' DOF angle          ',e12.5
 
c==== End of Subroutine ================================================
      
c9999  continue
      return
      end   
 

c=======================================================================      
c=======================================================================      

      subroutine yltloa(d,x,ix,id,f,ndf,ndm,prt)
      
c----------------------------------------------------------------------+
c     Perform Area Load, Line and Point Load                           |
c----------------------------------------------------------------------+
 
c==== DECLARATION ======================================================


c---- Formal Parameters ------------------------------------------------

      USE cdat1
      USE cdata
      USE iofile
      USE yltdata1
      USE yltdata2
      implicit real*8 (a-h,o-z)
      implicit integer (i,n)

c---- yltdata1 : data of ylt-element
c     nkn      : number of nodes
c     nth      : number of edges (theta) 
c     nknr     : reduced number of nodes (after condensation)
c     nthr     : reduced number of edges (after condensation)          
c     ndx      : number of derivations 
c     idx      : specify derivation: 1: node number, 2: 1=d/dx,2=d/dy 
c     ny1      : number of y_1 entries
c     alpha    : alpha history value
c     iter     : number of iteration    
c     aff      : affine distortion     

c---- yltdata2 : data of ylt-element 
c     dlm      : edge nodes, m_pl (<0: free edge)
c     angle    : angle of 3 element nodes

c---- cdat1    :
c     ndd      :
c     nie      :

c---- cdata    : 
c     numnp    : number of mesh nodes
c     numel    : number of mesh elements
c     nummat   : number of material sets
c     nen      : maximum nodes/element
c     neq      : number active equations
c     ipr      : real variable precision  

c---- iofile     : input-output-file
c     ior        : input file unit number (<0 input keyboard)
c     iow        : output file unit number
    
      integer ix(6,numel),id(ndf,numnp)

c---- ix(nen)    : element global node numbers 
c---- id         : boundary conditions  
c---- icnt       : count entries
      
      real*8 d(*),x(ndm,numnp),f(ndf,*)

c---- d(*)     : element data parameters
c---- x        : coordinates
c---- f        : load

      logical prt
 
c---- prt      : print option

c---- Variables --------------------------------------------------------

      integer isl
      
c---- isl      : switch load (1=area,2=line,3=point)      

      real*8 td(10),fl(ndm+1,3,16),v(3)
      
c---- td       : data parameters for input      
c---- fl       : local loads (ndf,3points+1load,no.inputs) !no.inputs possible to change 
c---- v        : auxiliary vector 


c==== PROGRAM ==========================================================

      call pzero(fl,3*4*16)
      tol=1.d-8                 

c==== Read Loads ======================================================= 
 
      if(ior.lt.0) write(*,8100)
      
      i0=0
  100 continue
      i0=i0+1
c      if (i0.gt.16) then !no.inputs
c      
c      
c      endif
      call pzero(td,10)
      call dinput(td(1:9),9)  ! x1,y1,,x2,y2,,x3,y3,,
      a1=0.d0
      do i1=1,9
        a1=a1+dabs(td(i1))
      enddo  
c      if (a1.eq.0.d0) goto 200  
      
      fl(1,1,i0)=td(1)
      fl(2,1,i0)=td(2)
      fl(3,1,i0)=td(3)
      fl(1,2,i0)=td(4)
      fl(2,2,i0)=td(5)
      fl(3,2,i0)=td(6)
      fl(1,3,i0)=td(7)
      fl(2,3,i0)=td(8)
      fl(3,3,i0)=td(9)   
      call pzero(td,10)
      call dinput(td(1:3),3)  ! q1,q2,q3
      v(1)=td(1)
      v(2)=td(2)
      v(3)=td(3)
      a2=0.d0
      do i1=1,3
        a2=a2+dabs(td(i1))
      enddo
      if ((a1.ge.0.d0.and.a1.le.0.d0).and.(a2.ge.0.d0.and.a2.le.0.d0)) 
     &   goto 200
        
c      fl(1,4,i0)=td(1)
c      fl(2,4,i0)=td(2)
c      fl(3,4,i0)=td(3)

      isl=0
      do i1=1,3
        if (v(i1).gt.0.d0.or.v(i1).lt.0.d0) isl=isl+1
      enddo  
      
cjw  110 continue
cjw      isl=3
cjw      do i1=1,2
cjw        do i2=i1+1,3
cjw          if ((fl(1,i1,i0).eq.fl(1,i2,i0)).and.
cjw     &        (fl(2,i1,i0).eq.fl(2,i2,i0)).and.
cjw     &        (fl(3,i1,i0).eq.fl(3,i2,i0))) then
cjw            isl=isl-1
cjw            if (i1.eq.1.and.i2.eq.2) then  !correction of order, for print
cjw              do i3=1,3
cjw                fl(i3,2,i0)=fl(i3,3,i0)
cjw                fl(i3,3,i0)=fl(i3,1,i0)
cjw              enddo
cjw              if ((fl(1,2,i0).ne.fl(1,3,i0)).and.
cjw     &            (fl(2,2,i0).ne.fl(2,3,i0)).and.
cjw     &            (fl(3,2,i0).ne.fl(3,3,i0))) goto 110
cjw            endif    
cjw          endif
cjw        enddo
cjw      enddo
      
c---- Distribute Loads -------------------------------------------------

      if (isl.eq.2) then
        fl(4,1,i0)=(2.d0/3.d0)*v(1)+(1.d0/3.d0)*v(2)
        fl(4,2,i0)=(2.d0/3.d0)*v(2)+(1.d0/3.d0)*v(1)
        fl(4,3,i0)=0.d0
      elseif (isl.eq.3) then
         !BAUSTELLE Tetraeder
        fl(4,1,i0)=v(1)
        fl(4,2,i0)=v(2)
        fl(4,3,i0)=v(3) 
      else
        fl(4,1,i0)=v(1)
        fl(4,2,i0)=0.d0
        fl(4,3,i0)=0.d0
      endif    
            
c---- Print ------------------------------------------------------------      
      
c      if (prt.and.ior.lt.0) then         !BAUSTELLE: nopr-Option
        if (isl.eq.3) then
c          write(*  ,8110) fl(1,1,i0),fl(2,1,i0),fl(3,1,i0),
          write(iow,8110) fl(1,1,i0),fl(2,1,i0),fl(3,1,i0),
     &                    fl(1,2,i0),fl(2,2,i0),fl(3,2,i0),  
     &                    fl(1,3,i0),fl(2,3,i0),fl(3,3,i0),
     &                    fl(4,1,i0),fl(4,2,i0),fl(4,3,i0)
        elseif (isl.eq.2) then
c          write(*  ,8120) fl(1,1,i0),fl(2,1,i0),fl(3,1,i0),
          write(iow,8120) fl(1,1,i0),fl(2,1,i0),fl(3,1,i0),
     &                    fl(1,2,i0),fl(2,2,i0),fl(3,2,i0),  
     &                    fl(4,1,i0),fl(4,2,i0)
        else !elseif (isl.eq.1) then
c          write(*  ,8130) fl(1,1,i0),fl(2,1,i0),fl(3,1,i0),
          write(iow,8130) fl(1,1,i0),fl(2,1,i0),fl(3,1,i0),
     &                    fl(4,1,i0)
        endif              
c      endif

      goto 100
      
c==== Search for Corresponding Coordinates and Compute f ===============


  200 continue
      do i1=1,numnp
        f(3,i1)=0.d0
        do i2=1,i0/2 !16 !no.inputs
          do i3=1,3
            if (dabs(fl(1,i3,i2)-x(1,i1)).lt.tol.and.
     &          dabs(fl(2,i3,i2)-x(2,i1)).lt.tol.and.
     &          dabs(fl(3,i3,i2)-x(3,i1)).lt.tol) then
              f(3,i1)=f(3,i1)+fl(4,i3,i2)
              id(1,i1)=1
              id(2,i1)=1
            endif
          enddo
        enddo
      enddo      


c---- test

c      do i1=1,i0
c        do i2=1,3
c          print*,'fl',fl(1,i2,i1),fl(2,i2,i1),fl(3,i2,i1),fl(4,i2,i1)
c        enddo
c      enddo  
c
c      do i1=1,numnp
c      print*, 'f', f(3,i1)
c      enddo

c==== Formats ==========================================================

8100  format(6x,'Y i e l d   L i n e   L o a d s')
8110  format(/,' LOAD generation with YLOA on area',/,
     &         '            Point 1 ',3(1x,e12.5),/,
     &         '            Point 2 ',3(1x,e12.5),/,
     &         '            Point 3 ',3(1x,e12.5),/,
     &         '            Load  1 ',   e12.5,/,
     &         '            Load  2 ',   e12.5,/,
     &         '            Load  3 ',   e12.5)
8120  format(/,' LOAD generation with YLOA on line',/,
     &         '            Point 1 ',3(1x,e12.5),/,
     &         '            Point 2 ',3(1x,e12.5),/,
     &         '            Load  1 ',   e12.5,/,
     &         '            Load  2 ',   e12.5)
8130  format(/,' LOAD generation with YLOA on point',/,
     &         '              Point ',3(1x,e12.5),/,
     &         '              Load  ',   e12.5)
   
c==== End of Subroutine ================================================
      
c9999  continue
      return
      end    

      
c=======================================================================      
c=======================================================================      

      subroutine yltass(d,u,x,ix,lm,id,f,s,p,atb,etb,ftb,ptb,ctb,
     &                  mtb,lth,isw)
      
c----------------------------------------------------------------------+
c     Assembly of Optimization Tableau                                 |
c----------------------------------------------------------------------+

c==== DECLARATION ======================================================

c---- Formal Parameters ------------------------------------------------

      USE cdat1
      USE cdata
      USE eldata
      USE sdata
      USE yltdata1
      USE yltdata2
      implicit real*8 (a-h,o-z)
      implicit integer (i,n)

c---- yltdata1 : data of ylt-element
c     nkn      : number of nodes
c     nth      : number of edges (theta) 
c     nknr     : reduced number of nodes (after condensation)
c     nthr     : reduced number of edges (after condensation)          
c     ndx      : number of derivations 
c     idx      : specify derivation: 1: node number, 2: 1=d/dx,2=d/dy 
c     ny1      : number of y_1 entries
c     alpha    : alpha history value
c     iter     : number of iteration         

c---- yltdata2 : data of ylt-element 
c     dlm      : edge nodes, m_pl (<0: free edge)
c     angle    : angle of 3 element nodes

c---- cdat1    :
c     ndd      :
c     nie      :

c---- cdata    : 
c     numnp    : number of mesh nodes
c     numel    : number of mesh elements
c     nummat   : number of material sets
c     nen      : maximum nodes/element
c     neq      : number active equations
c     ipr      : real variable precision  

c---- eldata     : element data
c     dm         : element proportional load
c     n          : current element number
c     ma         : current element material set
c     mct        : print counter
c     iel        : user element number
c     nel        : number nodes on current element  

c---- sdata      :
c     ndf        :
c     ndm        :
c     nen1       :
c     nst        :

c      real*8  d(ndd,*),u(*),x(*) 
c      integer id(*), ix(*)
      integer id(ndf,nkn), ix(6,numel),lm(2,nth)
      real*8 d(ndd,*),u(ndf,*),x(ndm,numnp)
      real*8 f(ndf,numnp),s(nst,nst),p(nst)

c      integer ld(*), ie(nie,*), id(ndf,*), ix(nen1,*), jp(*)
c      real*8  xl(ndm,*), p(*), s(nst,*), d(ndd,*), x(ndm,*), f(ndf,*)
c      real*8  b(*), a(*), c(*), ul(nst,*), tl(*), t(*), u(ndf,*)
c      real*8  f0(ndf,*), ud(*) 

c---- ix(nen)    : element global node numbers
c---- id(*)      : boundary conditions
c---- d(*)       : element data parameters
c---- u(*)       : 
c---- x(*)       : node coordinates
c---- t(nen)     : 
c---- s(nst,nst) : element matrix
c---- p(ndf,nen) : element vector

cc      real*8 atb(1,1),etb(1,1),ftb(1,1),ptb(1,1)
cc      real*8 ctb(1,1),mtb(1,1)
cc      real*8 vtb(1,1)
c      real*8 atb(2*nkn+2*nth,1+nth),etb(nkn,nkn),ftb(1,nkn)
c      real*8 ptb(1,1+nth),ctb(1,1),mtb(1,nth)
c      real*8 vtb(1,2*nkn+2*nth)


c---- Variables --------------------------------------------------------

      logical nwe,lchk
      
c---- nwe      : new entry to lm
c---- aff      : affin distortion done
c---- lchk     : element check

      integer ikn(3), ith(3)
      
c---- ikn,ith  : node/rotation positions in tableau      

      real*8 xl(ndm,ndf),t(1)!,p(3),s(3,3)

c---- xl       : element coordinates (3 nodes, 2d, linear field)      
c---- t        : element nodal temperature values (dummy)

      real*8 atb(1+nth,2*nkn+2*nth),etb(nth,nkn),ftb(nkn)
      real*8 ptb(1+nth),ctb(2*nkn+2*nth),mtb(2*nth)
c      real*8 vtb(2*nkn+2*nth)
      real*8 lth(3,nth) !,tm(3,3)

c---- atb      : table of constrains (th+th-w+w-,1+th) A
c---- etb      : assembled element matrix E
c---- ftb      : assembled element vector F 
c---- ptb      : right hand side
c---- ctb      : coefficients of objective function c
c---- mtb      : assembled element moment vector m
c---- vtb      : optimization variables (th+th-w+w-)
c---- lth      : length of yield-lines (x,y,sqrt(x^2+y^2))
c---- tm       : transformation matrix

c-----dummy
      real*8 dummy(1)

c==== PROGRAM ==========================================================

c---- Initialization ---------------------------------------------------

c      call pzero(lm,2*nth)
      lm=0
c      call pzero(vtb,2*nkn+2*nth)
      call pzero(atb,(2*nkn+2*nth)*(1+nth))
      call pzero(etb,nth*nkn)
      call pzero(ftb,nkn)
      call pzero(ctb,2*nkn+2*nth)
      call pzero(mtb,2*nth)
      call pzero(ptb,1+nth)
      ptb(1)=-1.d0
      call pzero(lth,3*nth)
      if (isw.eq.2) then
        lchk=.true.
        isw=30
      else
        lchk=.false.
      endif    

c---- Establish Location Matrix nth to nkn -----------------------------

      i0=1
      do i1=1,numel
        do ia=1,3
          ib=ia+1
          ic=ia+2
          if(ib.gt.3) ib=ib-3
          if(ic.gt.3) ic=ic-3
          nwe=.true.
          
          do i2=1,i0
            if((lm(1,i2).eq.ix(ia,i1)).and.(lm(2,i2).eq.ix(ib,i1)).or.
     &         (lm(2,i2).eq.ix(ia,i1)).and.(lm(1,i2).eq.ix(ib,i1))) then
              nwe=.false.
            endif
          enddo
          if(nwe) then
            lm(1,i0)=ix(ia,i1)
            lm(2,i0)=ix(ib,i1)
            i0=i0+1
          endif   
        enddo
      enddo  

c---- Affine Distortion of Coordinates ---------------------------------

        
c      do i1=1,3
c        do i2=1,numnp
c          if (i1.eq.2) then
c            xa(i1,i2)=d(1,ma)*x(i1,i2)
c          else
c            xa(i1,i2)=x(i1,i2)
c          endif
c        enddo
c      enddo        
      
c---- Compute Length of Yield-Lines ------------------------------------

c      if (isw.eq.30) then
        do i1=1,nth
          lth(1,i1)=x(1,lm(1,i1))-x(1,lm(2,i1))
          lth(2,i1)=x(2,lm(1,i1))-x(2,lm(2,i1))
          lth(3,i1)=dsqrt(lth(1,i1)*lth(1,i1)+lth(2,i1)*lth(2,i1))
        enddo
      if (isw.eq.31) then  ! d/dx
c        if (BAUSTELLE: Schalter isotrop) then
          do i1=1,nth !BAUSTELLE: momentan isotrop; affine Verzerrung einbauen? -> Sawczuk...
c            lth(idx(2),i1)=(x(idx(2),lm(1,i1))-x(idx(2),lm(2,i1))) !*1.d0 ! ???BAUSTELLE: Ableitung dlx/dx ggf nochmal überdenken
c            if (lth(idx(2),i1).ne.0) lth(idx(2),i1)=lth(idx(2),i1)/
c     &                               dabs(lth(idx(2),i1))
c            if (idx(2).eq.1) lth(2,i1)=0.d0
c            if (idx(2).eq.2) lth(1,i1)=0.d0
            lth(2,i1)=0.d0
            lth(1,i1)=0.d0
            if (lm(1,i1).eq.idx(1)) then !.or.lm(2,i1).eq.idx(1))
              lth(3,i1)=(x(idx(2),lm(1,i1))-x(idx(2),lm(2,i1)))/
     &                   (lth(3,i1))     
            elseif (lm(2,i1).eq.idx(1)) then
              lth(3,i1)=(x(idx(2),lm(2,i1))-x(idx(2),lm(1,i1)))/
     &                   (lth(3,i1))
            else
              lth(3,i1)=0.d0
            endif
          enddo    

c        else !anisotrop

c          Idee: Affine Verzerrung

c        endif
      endif  

      
c---- Represent Boundary Conditions read by Element in dlm -------------
c     and form Objective Function
c     Calculate general mtb's and special boundary mtb's...

      do i1=1,nth   

c        if (d(3,ma).le.(0.d0).and.d(5,ma).le.(0.d0)) then !isotrop      

          mtb(i1)=lth(3,i1)*dabs(d(3,ma))
          mtb(nth+i1)=lth(3,i1)*dabs(d(5,ma))

c        else !orthotrop
        
c         Affine Verzerrung, s.Sawczuk.... (Protokoll 20.1.05)  

c        endif          
        
        do i2=1,nth !Kantenlagerung
          if ((lm(1,i1).eq.int(dlm(1,i2))).and.
     &        (lm(2,i1).eq.int(dlm(2,i2))).or.
     &        (lm(2,i1).eq.int(dlm(1,i2))).and.
     &        (lm(1,i1).eq.int(dlm(2,i2)))) then 
            lm(1,i1)=-lm(1,i1)  !Kantenlagerung allg. -> -lm(2,i) 
            if (int(dlm(3,i2)).eq.0) then  !Freier Rand -> -lm(1,i)
              lm(2,i1)=-lm(2,i1)
              if(int(dlm(4,i2)).ne.0.and.int(dlm(5,i2)).ne.0) then 
                lm(1,i1)=-lm(1,i1)   !freier Rand mit m_pl fuer Symmetriefall (wieder+)
              endif
            endif                     !m_pl Kanten
            mtb(i1)=lth(3,i1)*dlm(4,i2)
            mtb(nth+i1)=lth(3,i1)*dlm(5,i2)                      
          endif
        enddo 
      enddo
                 
c---- Assemble Optimization Tableau ------------------------------------

      do n=1,numel
       
        call pzero(xl,ndm*ndf)     
     
c------ Choose node coordinates out of x -------------------------------      
      
          do i1=1,nen !=3
            do i2=1,ndm !=3
              xl(i2,i1)=x(i2,ix(i1,n))
c              angle(i1)=ang(ix(i1,n))*dacos(0.d0)/90
            enddo
          enddo                       

c------ Check Element --------------------------------------------------

        if (lchk) then     
          call elmlib(d(1,ma),u,xl,ix,t,s,p,dummy,dummy,dummy,
     1                ndf,ndm,nst,iel,2)
        endif

c------ Call Element ---------------------------------------------------

        call elmlib(d(1,ma),u,xl,ix,t,s,p,dummy,dummy,dummy,
     1              ndf,ndm,nst,iel,isw)

c------ Find Positions in E-Matrix -------------------------------------

        do i2=1,nth
          do ia=1,3
            ib=ia+1
            ic=ia+2
            if(ib.gt.3) ib=ib-3
            if(ic.gt.3) ic=ic-3
          
            if((abs(lm(1,i2)).eq.ix(ic,n)).and.(abs(lm(2,i2)).eq.
     &                                                   ix(ib,n)).or.
     &         (abs(lm(2,i2)).eq.ix(ic,n)).and.(abs(lm(1,i2)).eq.
     &                                                   ix(ib,n))) then
              ith(ia)=i2                !!!BAUSTELLE: hier (?) muss präziser definiert werden: genauer Ort, th gegenüber! 040105
              ikn(ia)=ix(ia,n)
c              i0=i0+1
            endif
          enddo    
        enddo       
        
c------ Fill into E-Matrix ---------------------------------------------

        do i2=1,3
          ftb(ikn(i2))=ftb(ikn(i2))+p(i2)        !Element Load    !  BAUSTELLE: m-Vektor, etc, ->maple!
          do i3=1,numnp
            if (ikn(i2).eq.i3) then
              ftb(ikn(i2))=ftb(ikn(i2))+f(3,i3)  !Point Load from yloa      !BAUSTELLE 24.8.05, testen
            endif
          enddo  
          do i1=1,3
            etb(ith(i1),ikn(i2))=etb(ith(i1),ikn(i2))+s(i1,i2)          
          enddo  
        enddo    
      enddo

c      call plotmatrix(etb,nth,nkn)
c      call plotmatrix(ftb,1,nkn)
          
c---- Form reduced E-Matrix --------------------------------------------

      do i2=1,nkn
        if (id(3,i2).le.0) ftb(i2)=0.d0  
        do i1=1,nth
c          if ((id(3,i2).lt.0).or.((lm(2,i1).le.0.and.lm(1,i1).le.0).and.
c     &        (lm(1,i1)*lm(2,i1).gt.0)))      
          if ((id(3,i2).lt.0).or.(lm(2,i1).le.0.and.lm(1,i1).le.0))  !HIER ENTSCHEIDET SICH THETA!!!
     &      etb(i1,i2)=0.d0
        enddo
      enddo
      
c      call plotmatrix(etb,nth,nkn)  
c      call plotmatrix(ftb,1,nkn)

c---- Count reduced Entries --------------------------------------------

      nthr=0
      do i1=1,nth
        if ((lm(1,i1).gt.0).and.(lm(1,i1).gt.0)) nthr=nthr+1   
      enddo
      nknr=0
      do i1=1,nkn
        if (id(3,i1).gt.0) nknr=nknr+1
      enddo   

c---- Form Tableau of Constraints --------------------------------------

      if (isw.eq.30) then
        do i1=1,nth
            atb(i1+1,i1)=      1.d0
          if ((lm(2,i1).gt.0).and.(lm(2,i1).gt.0)) then 
            atb(i1+1,i1+nth)= -1.d0
          endif
        enddo
      endif
      do i1=1,nkn
        atb(1,2*nth+i1)=      ftb(i1)
        atb(1,2*nth+i1+nkn)= -ftb(i1)
        do i2=1,nth
          atb(1+i2,2*nth+i1)=     -etb(i2,i1)
          atb(1+i2,2*nth+i1+nkn)=  etb(i2,i1)
        enddo  
      enddo  
      
c      call plotmatrix(atb,1+nth,2*nkn+2*nth, '-atb-') 

cc---- Form Tableau of Constraints --------------------------------------
c
c      do i1=1,nth
c          atb(i1,i1)=      1.d0
c        if (lm(1,i1).gt.0) then 
c          atb(i1,i1+nth)= -1.d0
c        endif
c      enddo
c      do i1=1,nkn
c        atb(nth+1,2*nth+i1)=      ftb(1,i1)
c        atb(nth+1,2*nth+i1+nkn)= -ftb(1,i1)
c        do i2=1,nth
c          atb(i2,2*nth+i1)=     -etb(i2,i1)
c          atb(i2,2*nth+i1+nkn)=  etb(i2,i1)
c        enddo  
c      enddo  
c      
cc      call plotmatrix(atb,1+nth,2*nkn+2*nth) 
      
c---- Form Objective Function ------------------------------------------

      do i1=1,2*nth
        ctb(i1)=mtb(i1)
      enddo           
      
c      call plotmatrix(ctb,2*nth,1) 

c==== End of Subroutine ================================================
      
c999   continue
      return
      end


c=======================================================================      
c=======================================================================      

      subroutine yltsqp(atb,btb,ctb,ptb,vtb,vtbh,isw)
      
c----------------------------------------------------------------------+
c     Solution of Optimization Tableau                                 |
c----------------------------------------------------------------------+

c==== DECLARATION ======================================================

c---- Formal Parameters ------------------------------------------------

      USE iofile
      USE yltadr
      USE yltdata1
      USE yltdata2
      USE yltdata3
      implicit real*8 (a-h,o-z)
      implicit integer (i,n)

c---- yltdata1 : data of ylt-element
c     nkn      : number of nodes
c     nth      : number of edges (theta) 
c     nknr     : reduced number of nodes (after condensation)
c     nthr     : reduced number of edges (after condensation)          
c     ndx      : number of derivations 
c     idx      : specify derivation: 1: node number, 2: 1=d/dx,2=d/dy 
c     ny1      : number of y_1 entries
c     alpha    : alpha history value
c     iter     : number of iteration    
c     aff      : affine distortion     

c---- yltdata2 : data of ylt-element 
c     dlm      : edge nodes, m_pl (<0: free edge)
c     angle    : angle of 3 element nodes

c---- yltdata3 : tolerances   
c     eps1     : 1st tolerance check
c     eps2     : 2nd tolerance check
c     eps3     : degereration criterion 
c     eps4     : alpha bracketing tolerance 
c     eps5     : step length for direct search 
c     eps6     : zero tolerance for sqp
c     halt     : iterations stopped

      integer isw
      real*8 atb(1+nth,2*nkn+2*nth) !,etb(nth,nkn),ftb(nkn)
      real*8 ptb(1+nth),ctb(2*nkn+2*nth) !,mtb(2*nth)
      real*8 vtb(2*nkn+2*nth),vtbh(2*nkn+2*nth)
c      real*8 utb(nth+1+4*nth+4*nkn)   
c      real*8 war(3*(2*nth+2*nkn)*(2*nth+2*nkn)/2+10*(2*nth+2*nkn)
c     &           +2*(nth+1)+1)
      real*8 war(100000)
      real*8 xl(2*nth+2*nkn),xu(2*nth+2*nkn)
c      integer iwar(3*(2*nth+2*nkn)*(2*nth+2*nkn)/2+10*(2*nth+2*nkn)
c     &             +2*(nth+1)+1)
      integer iwar(100000)
      integer lwar

c---- atb      : table of constrains (th+th-w+w-,1+th) A
c---- etb      : assembled element matrix E
c---- ftb      : assembled element vector F 
c---- ptb      : right hand side
c---- ctb      : coefficients of objective function c
c---- mtb      : assembled element moment vector m
c---- vtb      : optimization variables (th+th-w+w-)
c---- vtbh   : previous vtb


c==== PROGRAM ==========================================================

c---- Initialization ---------------------------------------------------

c      call pzero(war,10000)
c      call pzero(iwar,10000)
c      call plotmatrix(ctb,1,2*nkn+2*nth) 
c      call plotmatrix(atb,1+nth,2*nkn+2*nth) 
      call pzero(btb,(2*nth+2*nkn)*(2*nth+2*nkn))


c---- Predefine Parameters ---------------------------------------------
          
      m=nth+1
      me=m
      mmax=nth+1    
      n=2*nth+2*nkn
      nmax=n
      mnn=m+n+n
      do i1=1,n
        xl(i1)=0.d0
        xu(i1)=1.d124
      enddo
c      ifail=0
c
c  100 continue    
c      if (ifail.eq.0) then
c        eps=1.d-8 !27
c      else
c        eps=eps*1.d-2
c      endif    

c---- Give previous vtb to vtbh --------------------------------------

      if (isw.ne.0) then !called by macro, else called by yltalp
        if (iter.gt.1) then
          do i1=1,2*nth+2*nkn
            vtbh(i1)=vtb(i1)
          enddo  
        else
          call pzero(vtbh,2*nth+2*nkn)  
        endif
      endif  

c---- Call SQP-Algorithm -----------------------------------------------

c   usage:
c
c      ql0001(m,me,mmax,n,nmax,mnn,c,d,a,b,xl,xu,x,u,iout,ifail,iprint,
c             war,lwar,iwar,liwar)

c   m :        total number of constraints.
c   me :       number of equality constraints.
c   mmax :     row dimension of a. mmax must be at least one and greater
c              than m.
c   n :        number of variables.
c   nmax :     row dimension of c. nmax must be greater or equal to n.
c   mnn :      must be equal to m + n + n.
c   c(nmax,nmax): objective function matrix which should be symmetric and
c              positive definite. if iwar(1) = 0, c is supposed to be the
c              choleskey-factor of another matrix, i.e. c is upper
c              triangular.
c   d(nmax) :  contains the constant vector of the objective function.
c   a(mmax,nmax): contains the data matrix of the linear constraints.
c   b(mmax) :  contains the constant data of the linear constraints.
c   xl(n),xu(n): contain the lower and upper bounds for the variables.
c   x(n) :     on return, x contains the optimal solution vector.
c   u(mnn) :   on return, u contains the lagrange multipliers. the first
c              m positions are reserved for the multipliers of the m
c              linear constraints and the subsequent ones for the 
c              multipliers of the lower and upper bounds. on successful
c              termination, all values of u with respect to inequalities 
c     

      iwar(1)=1
c      lwar=3*nmax*nmax/2 + 10*nmax + 2*mmax+1 
      lwar=3*(2*nth+2*nkn)*(2*nth+2*nkn)/2+10*(2*nth+2*nkn)+2*(nth+1)+1 
      lwar=min(lwar,100000)
      
      call ql0001(m,me,mmax,n,nmax,mnn,
     &           btb,ctb,atb,ptb,xl,xu,
     &           vtb,utb,iow,ifail,1,war,lwar,iwar,10000,eps6)

c      if (ifail.eq.2) goto 100

c---- Correction of small values ---------------------------------------

      do i1=1,nth
        if (dabs(vtb(i1)).ge.dabs(vtb(i1+nth))) vtb(i1+nth)=0.d0      
        if (dabs(vtb(i1)).lt.dabs(vtb(i1+nth))) vtb(i1)=0.d0      
      enddo
      do i1=2*nth+1,2*nth+nkn
        if (dabs(vtb(i1)).ge.dabs(vtb(i1+nkn))) vtb(i1+nkn)=0.d0      
        if (dabs(vtb(i1)).lt.dabs(vtb(i1+nkn))) vtb(i1)=0.d0      
      enddo
      
      do i1=1,2*nth+nkn     !raus wg Basiskontrolle !rein wg Singularität A_1, ny1...rein!!!
        if (dabs(vtb(i1)).lt.eps6) vtb(i1)=0.d0
      enddo  

cc---- test
c
c      print*,'ifail=',ifail,' eps=',eps,' nth=',nth,' nkn=',nkn  
c      do i1=1, 2*nth+2*nkn
c        print*, 'vtb(',i1,')=',vtb(i1),' ;mtb(',i1,')=',mtb(i1)
c      enddo

ctemp070505>>

c      ny1=0
c      do i1=1,2*nkn+2*nth
c        if (vtb(i1).ne.0.d0) ny1=ny1+1
c      enddo                             

c---- Count number of y_1 entries --------------------------------------

c      ny1=1
c      do i1=1,nth
c        if (lm(1,i1).gt.0) ny1=ny1+1
c      enddo 
c
c      print*, 'ny1 alt=',ny1

      if (isw.ne.0) then
        ny1=0
        do i1=1,2*nth+2*nkn
          if (vtb(i1).gt.0.d0.or.vtb(i1).lt.0.d0) ny1=ny1+1
        enddo
      endif 
       
c      print*, 'ny1=',ny1  
 
c      call plotmatrix(vtb,2*nth+2*nkn,1,'vtb')
 
c==== End of Subroutine ================================================
      
c999   continue
      return
      end


c=======================================================================      
c=======================================================================      

      subroutine yltevl(vtb,mtb,u,th,qu,isw)
      
c----------------------------------------------------------------------+
c     Evaluate YLT Results                                             |
c----------------------------------------------------------------------+

c==== DECLARATION ======================================================

c---- Formal Parameters ------------------------------------------------
      USE cdata
      USE sdata
      USE yltdata1
      implicit real*8 (a-h,o-z)
      implicit integer (i,n)

c---- cdata    : 
c     numnp    : number of mesh nodes
c     numel    : number of mesh elements
c     nummat   : number of material sets
c     nen      : maximum nodes/element
c     neq      : number active equations
c     ipr      : real variable precision  

c---- sdata      :
c     ndf        :
c     ndm        :
c     nen1       :
c     nst        :

c---- yltdata1 : data of ylt-element
c     nkn      : number of nodes
c     nth      : number of edges (theta) 
c     nknr     : reduced number of nodes (after condensation)
c     nthr     : reduced number of edges (after condensation)          
c     ndx      : number of derivations 
c     idx      : specify derivation: 1: node number, 2: 1=d/dx,2=d/dy 
c     ny1      : number of y_1 entries
c     alpha    : alpha history value
c     iter     : number of iteration    
c     aff      : affine distortion     

      integer isw
      real*8 vtb(2*nkn+2*nth),mtb(2*nth)
      real*8 u(ndf,numnp),th(nth),qu

c---- isw      : isw=0: dummy calculation for alpha; else: real calc.

c==== PROGRAM ==========================================================

      call pzero(th,nth)

c---- Transfer to Process ----------------------------------------------

c  100 continue
      goto (1000,2000,3000), isw   

c==== Give nodal Displacements from vtb to u(3,*) and th(*) ============
c     for alpha dummie (isw=2): skip

 1000 continue
c      if (isw.ne.0) then
        do i1=1,nkn
          if (vtb(2*nth+i1).gt.0.d0.or.vtb(2*nth+i1).lt.0.d0) 
     &      u(3,i1)=-vtb(2*nth+i1)
          if (vtb(2*nth+nkn+i1).gt.0.d0.or.vtb(2*nth+nkn+i1).lt.0.d0) 
     &      u(3,i1)=vtb(2*nth+nkn+i1)
        enddo     
        do i1=1,nth
          if (vtb(i1).gt.0.d0.or.vtb(i1).lt.0.d0) th(i1)=vtb(i1)
          if (vtb(nth+i1).gt.0.d0.or.vtb(nth+i1).lt.0.d0) 
     &      th(i1)=-vtb(nth+i1)
        enddo   
c      endif  
      goto 2000
        
c==== Calculate Ultimate Load ==========================================

 2000 continue
      qu=0.d0
      do i1=1,2*nth
        qu=qu+vtb(i1)*mtb(i1)
      enddo  
      return
      
c==== Evaluate YL-Pattern for Genetic Algorithm ========================

 3000 continue
      
c      return      
            
c==== End of Subroutine ================================================
      
c 9999 continue
      return
      end


c=======================================================================      
c=======================================================================      

      subroutine yltprt(lm,th,qu,dqu,isw)
      
c----------------------------------------------------------------------+
c     Print YLT Results                                                |
c----------------------------------------------------------------------+

c==== DECLARATION ======================================================

c---- Formal Parameters ------------------------------------------------
      USE cdata
      USE iofile
      USE sdata
      USE yltdata1
      implicit real*8 (a-h,o-z)
      implicit integer (i,n)

c---- cdata    : 
c     numnp    : number of mesh nodes
c     numel    : number of mesh elements
c     nummat   : number of material sets
c     nen      : maximum nodes/element
c     neq      : number active equations
c     ipr      : real variable precision  

c---- sdata      :
c     ndf        :
c     ndm        :
c     nen1       :
c     nst        :

c---- yltdata1 : data of ylt-element
c     nkn      : number of nodes
c     nth      : number of edges (theta) 
c     nknr     : reduced number of nodes (after condensation)
c     nthr     : reduced number of edges (after condensation)          
c     ndx      : number of derivations 
c     idx      : specify derivation: 1: node number, 2: 1=d/dx,2=d/dy 
c     ny1      : number of y_1 entries
c     alpha    : alpha history value
c     iter     : number of iteration    
c     aff      : affine distortion     

      integer lm(2,nth)
      
c---- lm       : location matrix yl -> nodes      

c      real*8 atb(1+nth,2*nkn+2*nth)
c      real*8 ctb(2*nkn+2*nth),mtb(2*nth)
c      real*8 vtb(2*nkn+2*nth)
      real*8 th(nth),qu
      real*8 dqu(ndx),dp(ndx),sca

c---- atb      : table of constrains (th+th-w+w-,1+th) A
c---- ctb      : coefficients of objective function c
c---- mtb      : assembled element moment vector m
c---- vtb      : optimization variables (th+th-w+w-)
c---- th       : rotation of yield lines
c---- qu       : ultimate load
c---- dqu      : gradient vector
c---- dp       : scaled gradient (not transmitted)
c---- sca      : scaling factor

c==== PROGRAM ==========================================================

c---- Transfer to Process ----------------------------------------------

      goto (100,200,300,400,500), isw

c==== [ytab] (Form Tableau) print Ultimate Load ========================

100   continue
                   write(iow,2010) iter  !BAUSTELLE: nopr-Option einbauen
      if(ior.lt.0) write(*  ,2010) iter
                   write(iow,2011) qu
      if(ior.lt.0) write(*  ,2011) qu
      return

c==== [yang] Print theta Values ========================================

200   continue
                   write(iow,2020) 
      if(ior.lt.0) write(*  ,2020) 
      do i1=1,nth        
                     write(iow,2021) abs(lm(1,i1)),abs(lm(2,i1)),th(i1)  !BAUSTELLE: nopr-Option einbauen
        if(ior.lt.0) write(*  ,2021) abs(lm(1,i1)),abs(lm(2,i1)),th(i1)
      enddo
      return
      
c==== [ygra] Print Gradient ============================================

300   continue    

c---- Scaled Gradient --------------------------------------------------

      sca=0.d0
      do i1=1,ndx
        if (dabs(dqu(i1)).gt.dabs(sca)) then
          sca=dqu(i1)
        endif
      enddo
      do i1=1,ndx
        dp(i1)=dqu(i1)/sca
      enddo      
      
c---- Print Gradient ---------------------------------------------------

                   write(iow,2030) 
      if(ior.lt.0) write(*  ,2030) 
      do i1=1,ndx        
                     write(iow,2031) i1,dqu(i1),dp(i1)  !BAUSTELLE: nopr-Option einbauen
        if(ior.lt.0) write(*  ,2031) i1,dqu(i1),dp(i1)
      enddo

      return  
      
c==== [yevo] Genetic Algorithm =========================================

  400 continue
                   write(iow,2010) iter  !BAUSTELLE: nopr-Option einbauen
      if(ior.lt.0) write(*  ,2010) iter
                   write(iow,2011) qu
      if(ior.lt.0) write(*  ,2011) qu
      return
 
c==== [ytry] Direct Search =============================================

  500 continue
      return      

c==== Formats ==========================================================    

 2010 format('  Simplex iteration step = ',1i3)
 2011 format('  Ultimate load q_u      = ',1pe15.7)
 2020 format(/'  y i e l d    l i n e    r o t a t i o n s'/
     &       '  line nodes',5x,'theta value')
 2021 format(3x,1i3,2x,1i3,4x,1pe15.7)
 2030 format(/'  Gradient evaluation'/
     &'  dof X     dq_u/dX           scaled')
 2031 format(3x,1i3,5x,1pe15.7,2x,1pe15.7)

c==== End of Subroutine ================================================
      
c 9999 continue
      return
      end
      

c=======================================================================      
c=======================================================================      

      subroutine ylttec(x,u,qu,ix)
      
c----------------------------------------------------------------------+
c     Tec Interface                                                    |
c----------------------------------------------------------------------+

c==== DECLARATION ======================================================

c---- Formal Parameters ------------------------------------------------
      USE cdata
      USE comfil
      USE iofile
      USE sdata
      USE yltdata1
      implicit real*8 (a-h,o-z)
      implicit integer (i,n)
      
      logical lylttec

c---- cdata    : 
c     numnp    : number of mesh nodes
c     numel    : number of mesh elements
c     nummat   : number of material sets
c     nen      : maximum nodes/element
c     neq      : number active equations
c     ipr      : real variable precision  

c---- sdata      :
c     ndf        :
c     ndm        :
c     nen1       :
c     nst        :

c---- yltdata1 : data of ylt-element
c     nkn      : number of nodes
c     nth      : number of edges (theta) 
c     nknr     : reduced number of nodes (after condensation)
c     nthr     : reduced number of edges (after condensation)          
c     ndx      : number of derivations 
c     idx      : specify derivation: 1: node number, 2: 1=d/dx,2=d/dy 
c     ny1      : number of y_1 entries
c     alpha    : alpha history value
c     iter     : number of iteration    
c     aff      : affine distortion     

c---- comfil   : common files
c     fsav     : save file -> tec

      common /ylttecdat/ lylttec
      
c---- ylttecdat: info about first call      

      real*8 x(ndm,numnp),u(ndm,numnp),qu
      integer ix(6,numel)
      
c---- x        : coordinates      
c---- u        : displacements
c---- qu       : limit load
c---- ix       : nodes of element      



c==== PROGRAM ==========================================================

      
      fsav='c:\temp\out.tec'
    
      if (.not.lylttec) then
        open(unit=23,file=fsav,status='replace',form='formatted')
        rewind(unit=23)
        write(unit=23,fmt='(a)') 'TITLE="CYLT"'       
        write(unit=23,fmt='(a)') 'VARIABLES="X" "Y" "U" "Q"'      
        write(unit=23,fmt='(a36,i3,a4,i3,a21)') 
     &   'ZONE T="CYLT", DATAPACKING=POINT, N=',numnp,', E=',numel,
     &   ', ZONETYPE=FETriangle' 
        lylttec=.true.
      else        
cww     next statement is ok for SAL, but not ok for INTEL        
cww     open(unit=23,file=fsav,  status='append',form='formatted')
        open(unit=23,file=fsav,position='append',form='formatted')
        write(unit=23,fmt='(a36,i3,a4,i3,a21)') 
     &   'ZONE T="CYLT", DATAPACKING=POINT, N=',numnp,', E=',numel,
     &   ', ZONETYPE=FETriangle' 
      endif

      do i1=1,numnp
        write(unit=23,fmt='(4f10.4)',err=990) x(1,i1),x(2,i1),u(3,i1),qu
      enddo
      write(unit=23,fmt='(a)')
      do i1=1,numel
        write(unit=23,fmt='(3i3)',err=990) ix(1,i1),ix(2,i1),ix(3,i1)
      enddo
      write(unit=23,fmt='(a)')

      close(unit=23)   
      goto 999   

c==== Error ============================================================

  990 continue
      print*, 'Error while writing Tecplot-file'

c==== End of Subroutine ================================================
      
  999 continue
      return
      end
