
c----------------------------------------------------------------------+
c                                                                      |
c     Special Subroutines for Yield-Line Element  jw020605             |
c                                                                      |
c     Part 3: Genetic Algorithm                                        |
c                                                                      |
c----------------------------------------------------------------------+


c=======================================================================      

      subroutine yltevo(d,u,x,ix,lm,id,f,s,p,atb,etb,ftb,ptb,ctb,mtb,
     &                  vtb,lth,qu,pop,isw)
      
c----------------------------------------------------------------------+
c     Perform Genetic Algorithm                                        |
c     Weicker: Evolutionaere Algorithmen, S.120                        |
c----------------------------------------------------------------------+

c==== DECLARATION ======================================================

c---- Formal Parameters ------------------------------------------------
      USE cdat1
      USE cdata
      USE eldata
      USE sdata
      USE yltdata1
      USE yltdata2
      USE yltdata5
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

c---- yltdata5 : genetic data
c     npop     : number of individuals
c     nprp     : number of properties (=ndf)
c     nrec     : number of genes (=nkn)   
c     ncrt     : number of evaluation criterions
c     imut     : mutation switch
c     irec     : recombination switch  

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

      integer id(ndf,nkn), ix(6,numel),lm(2,nth)
      real*8 d(ndd,*),u(ndf,*),x(ndm,numnp)
      real*8 f(ndf,numnp),s(nst,nst),p(nst) 
      real*8 atb(1+nth,2*nkn+2*nth),etb(nth,nkn),ftb(nkn)
      real*8 ptb(1+nth),ctb(2*nkn+2*nth),mtb(2*nth)
      real*8 vtb(2*nkn+2*nth)
      real*8 lth(3,nth),qu,pop(ndm,nkn,npop)

c---- atb      : table of constrains (th+th-w+w-,1+th) A
c---- etb      : assembled element matrix E
c---- ftb      : assembled element vector F 
c---- ptb      : right hand side
c---- ctb      : coefficients of objective function c
c---- mtb      : assembled element moment vector m
c---- vtb      : optimization variables (th+th-w+w-)
c---- lth      : length of yield-lines (x,y,sqrt(x^2+y^2))
c---- qu       : actual ultimate load
c---- pop      : population pop(nprp,nrec,npop)

c---- Variables --------------------------------------------------------

      integer max
      integer vmut(ndm,nkn),vmut1(ndm,nkn),dmut(nprp,nrec)
c      real*8 x0(ndm,numnp),que,qui,quh
      real*8 step(ndm) !,th(nth)
      real*8 vsor(ncrt+1,npop)
      real*8 pop1(ndm,nkn,npop),vsor1(ncrt+1,npop)
      real*8 pop2(ndm,nkn,npop),vsor2(ncrt+1,npop)
      real*8 par(ndm,nkn,2) !,psor(ncrt+1,2)
c      real*8 ind1(ndm,nkn),ind2(ndm,nkn)
      real*8 crs

c---- npop     : number of individuals in population
c---- imut     : mutational iteration
c---- irec     : recombination iteration
c---- max      : maximum of mutations/recombinations
c---- vmut     : allowed mutation steps (0=not allowed,1=allowed)
c---- dmut     : dependence of dof steps in mutation of gene properties
c               (0=no dependence,else=scaling factor)
c---- x0       : original mesh
c---- que      : external load
c---- qui      : internal load
c---- quh      : historic qu
c---- vevl     : evaluation vector
c---- step     : step length for mutation 
c---- th       : yield-line rotations
c---- vsor     : sorting vector for population (ind no, qu)
c---- cpop     : cloned population for recombination
c---- cvsor    : cloned evaluation vector
c---- opop     : old generation
c---- ovsor    : evaluation vector old generation
c---- par      : parents/children
c---- psor     : evaluation children
c---- ind*     : individual
c---- crs      : crossover probability

      logical loop
      
c---- loop     : new mutation required           
          
c==== PROGRAM ==========================================================

c---- Set Tolerance ----------------------------------------------------

      tol=1.d-8

c---- Set Switches -----------------------------------------------------

      iswmut=2 !mutation switch
      iswrec=2 !gene exchange (1=single, 2=sequential)
      iswcmb=1 !
      iswgen=1 !
      npai=4
      nseq=4  
      
c---- Set Dimensions And Initialization --------------------------------

      nrec=nkn
      nprp=ndm         
      
      iter=iter+1
      loop=.true.

c==== Genesis ==========================================================

c---- Save First Mesh Coordinates --------------------------------------

c      x0=x
      iloop=0
      ipop=0  !same individual value...
      iind=0

c---- New Recombination Vector -----------------------------------------
                    
      call pzero(vmut,ndm*nkn)
      do i1=1,ndm
        do i2=1,nkn
          if (id(i1,i2).lt.0) then 
            vmut(i1,i2)=0
          else
            vmut(i1,i2)=1     
          endif
        enddo   
      enddo  
      
      step(1)=0.d0
      step(2)=0.d0
      step(3)=1.d0

c---- Create Individuals -----------------------------------------------

 1000 continue
      do i1=1,npop
        do i2=1,nkn
          do i3=1,ndm
            pop(i3,i2,i1)=x(i3,i2)
          enddo
        enddo

        call randomr(0.d0,1.d0,a)

        do i2=1,nkn
          do i3=1,ndm
c            a=1.d0
            pop(i3,i2,i1)=pop(i3,i2,i1)+a*step(i3)*vmut(i3,i2)  
          enddo
        enddo    
      enddo
      
      iind=iind+1
c      print*, 'new individuals created', iind
      
c==== Algorithm ========================================================

 2000 continue

c---- Evaluate Population ----------------------------------------------
      
      call yltval(d,u,x,ix,lm,id,f,s,p,atb,etb,ftb,ptb,ctb,mtb,vtb,
     &                lth,pop,vsor,npop,ncrt)
                    

c---- Save Parent Generation -------------------------------------------

c      call pzero(pop1,ndm*nkn*npop)
      pop1=0
      pop1=pop
      vsor1=vsor
             
c---- Select Parent Individuals ----------------------------------------         
                
c      call pzero(pop2,ndm*nkn*npop)
      pop2=0
      
      do i1=1,int(npop/2)

c        call yltsel(pop,n1,n2,npop,ix,isw)

        call randomi(1,npop,n1)
 2010   continue
        call randomi(1,npop,n2)
        if (n2.eq.n1) goto 2010

        do i2=1,nkn
          do i3=1,ndm
            par(i3,i2,1)=pop(i3,i2,n1)
            par(i3,i2,2)=pop(i3,i2,n2)
          enddo
        enddo        
        
c---- Recombination or Mutation of Two Parents -------------------------

        call randomr(0.d0,1.d0,crs)
        
        if (vsor(3,n1).gt.crs.or.vsor(3,n2).gt.crs) then

          call genalg(par,vmut,step,ndm,nkn,2,npai,nseq,
     &                iswmut,iswrec,iswcmb,iswgen,2)
        endif

        call genalg(par,vmut,step,ndm,nkn,2,npai,nseq,
     &              iswmut,iswrec,iswcmb,iswgen,1) 
    
        
c---- Fill Into New Generation -----------------------------------------          

        do i2=1,nkn
          do i3=1,ndm
            pop2(i3,i2,2*i1-1)=par(i3,i2,1)
            pop2(i3,i2,2*i1)  =par(i3,i2,2)
          enddo
        enddo
      enddo 
     
      if (2*int(npop/2).lt.npop) then
        do i1=1,nkn
          do i2=1,ndm
            pop2(i2,i1,npop)=pop(i2,i1,1)
          enddo
        enddo
        do i1=1,ncrt+1
          vsor2(i1,npop)=vsor(i1,1)
        enddo  
      endif      
      
c---- Normalize --------------------------------------------------------

      do i1=1,npop
        a=0.d0
        do i2=1,nkn
          if (dabs(pop2(3,i2,i1)).gt.a) a=dabs(pop2(3,i2,i1))
        enddo
        if (a.gt.0.d0.or.a.lt.0.d0) then
          do i2=1,nkn
            if (vmut(3,i2).ne.0) pop2(3,i2,i1)=pop2(3,i2,i1)/a
          enddo
        endif
      enddo            
 
c---- Evaluate New Generation ------------------------------------------

      call yltval(d,u,x,ix,lm,id,f,s,p,atb,etb,ftb,ptb,ctb,mtb,vtb,
     &                lth,pop2,vsor2,npop,ncrt)

c---- If New Generation is Worse Take Old One --------------------------

      if (vsor2(2,1).lt.vsor1(2,1)) then ! neue Generation besser
         
        pop=pop2
        vsor=vsor2
        iloop=iloop+1
      
      else   ! alte Generation besser

        pop=pop1
        vsor=vsor1
        iloop=iloop+1
      
      endif
      
c---- If Individual Values Do Not Improve ------------------------------

      if (val.ge.vsor(2,1).and.val.le.vsor(2,1)) then
        ipop=ipop+1
      endif 
      val=vsor(2,1)     

c---- If Half Of Population Is Similar Then Fix Next Level And Mutate --

c      print*, 'ipop', ipop
c      print*, '  vsor-Kontrolle:        ',vsor(2,1),vsor(2,int(npop/2))
      if (ipop.gt.npop.or.
     &    dabs(vsor(2,1)-vsor(2,int(npop/2))).lt.tol) then
        
c        print*, 'bin in der schleife......'
        ipop=0
        vmut1=vmut
        do i1=1,numel
          ia=ix(1,i1)
          ib=ix(2,i1)
          ic=ix(3,i1)
          if (vmut1(3,ia).eq.0.or.vmut1(3,ib).eq.0) then
            vmut(3,ic)=0
          endif
          if (vmut1(3,ia).eq.0.or.vmut1(3,ic).eq.0) then
            vmut(3,ib)=0
          endif
          if (vmut1(3,ib).eq.0.or.vmut1(3,ic).eq.0) then
            vmut(3,ia)=0
          endif
        enddo
        n1=0
        do i1=1,nkn
          n1=n1+vmut(3,i1)
        enddo

        
c        print*, ' Next Level reached:', n1,n2
        
        if (n1.eq.0) then !next dimesion
        
c          print*, 'converged!'
          loop=.false.
        
        else
c          print*, 'not converged yet...'
          call randomr(-1.d0,1.d0,a)
          step(3)=step(3)+a
          loop=.true.
          
        endif  
      else
        goto 2000  
      endif      
 
c==== End of Evolution =================================================

c 4000 continue 
 
      a=0.d0
      do i1=1,npop
        if (vsor(2,i1).gt.a) a=vsor(2,i1)
      enddo  
      do i1=1,npop
        if (i1.eq.1.or.(vsor(2,i1).lt.a.and.vsor(2,i1).gt.0.d0)) then
          a=vsor(2,i1)
          qu=vsor(2,i1)
          n1=vsor(1,i1)
        endif  
      enddo
      do i1=1,nkn
        do i2=1,ndm
          x(i2,i1)=pop(i2,i1,n1)
        enddo
      enddo    

      a=0.d0
      do i1=1,nkn
        if (dabs(x(3,i1)).gt.a) a=dabs(x(3,i1))
      enddo  
      if (a.ge.0.d0.and.a.le.0.d0) a=1.d0
      
      if (loop) then
        goto 1000  
      else
        do i1=1,nkn
          u(ndf,i1)=-x(3,i1)/a
          x(ndm,i1)=0.d0
        enddo  
      endif  

c---- control

c        call plotmatrix(vsor,ncrt+1,npop,'vsor')
c      
c      print*, 'qu=', qu

c==== End of Subroutine ================================================
      
c 9999 continue
      return
      end   


c=======================================================================      

      subroutine yltsel(pop,n1,n2,npop,ix,isw)
      
c----------------------------------------------------------------------+
c     Search Parent Pair                                               |
c----------------------------------------------------------------------+

c==== DECLARATION ======================================================

c---- Formal Parameters ------------------------------------------------
      USE cdat1
      USE cdata
      USE eldata
      USE sdata
      USE yltdata1
      USE yltdata2
c      USE yltdata5
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

c---- yltdata5 : genetic data
c     npop     : number of individuals
c     nprp     : number of properties (=ndf)
c     nrec     : number of genes (=nkn)   
c     ncrt     : number of evaluation criterions
c     imut     : mutation switch
c     irec     : recombination switch  

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

      integer ix(6,numel)
      real*8 pop(ndm,nkn,npop) !,par(ndm,nkn,2)
      
c---- pop        : existing population
c---- par        : new selected pair


c---- Variables --------------------------------------------------------

      integer ipar(npop)
      real*8 sum !,vevl(2,npop)
      
c---- ipar       : array of selected partners
c---- sum        : sum of fitness 
c---- vevl       : scaled fitness     

c==== PROGRAM ==========================================================

c---- Find First Point -------------------------------------------------

  100 continue
      call randomi(1,npop,n1)
  
c---- Search Second Point ----------------------------------------------

      call pzero(ipar,npop)
      
      n0=0
      do i1=1,numel
        do ia=1,3
          ib=ia+1
          ic=ia+2
          if (ib.gt.3) ib=ib-3
          if (ic.gt.3) ic=ic-3
          if (ix(ia,i1).eq.n1) then
            ipar(n0+1)=ix(ib,i1)
            ipar(n0+2)=ix(ic,i1)
            n0=n0+2
          endif
        enddo
      enddo
      
      n0=0
      do i1=1,numel            
        if (ipar(i1).ne.0) n0=n0+1
      enddo  
      if (n0.eq.0) goto 100

c  200 continue      
      call randomi(1,n0,i2)
      n2=ipar(i2) !zweiter Senkpunkt

c---- Search For Potential Parents -------------------------------------

      call pzero(ipar,npop)
      n0=0
      do i1=1,npop
        if (pop(3,n2,i1).gt.0.d0.or.pop(3,n2,i1).lt.0.d0) then
          n0=n0+1
          ipar(n0)=i1
        endif  
      enddo
      if (n0.eq.0) goto 100
      
      call randomi(1,n0,n2)
c      do i1=1,ndm
c        do i2=1,nkn
c          par(i1,i2,1)=pop(i1,i2,n1)
c          par(i1,i2,2)=pop(i1,i2,n2)
c        enddo
c      enddo        

      print*, 'PARENTS  ',n1,n2

c==== End of Subroutine ================================================
      
c 9999 continue
      return
      end   


cc=======================================================================      
c
c      subroutine yltsmo(pop,npop,ix)
c      
cc----------------------------------------------------------------------+
cc     Smooth Elements's Slope                                          |
cc----------------------------------------------------------------------+
c
cc==== DECLARATION ======================================================
c
cc---- Formal Parameters ------------------------------------------------
c
c      USE cdat1
c      USE cdata
c      USE eldata
c      USE sdata
c      USE yltdata1
c      USE yltdata2
cc      USE yltdata5
c      implicit real*8 (a-h,o-z)
c      implicit integer (i,n)
c
c      
cc---- yltdata1 : data of ylt-element
cc     nkn      : number of nodes
cc     nth      : number of edges (theta) 
cc     nknr     : reduced number of nodes (after condensation)
cc     nthr     : reduced number of edges (after condensation)          
cc     ndx      : number of derivations 
cc     idx      : specify derivation: 1: node number, 2: 1=d/dx,2=d/dy 
cc     ny1      : number of y_1 entries
cc     alpha    : alpha history value
cc     iter     : number of iteration    
cc     aff      : affine distortion     
c      
c      
cc---- yltdata2 : data of ylt-element 
cc     dlm      : edge nodes, m_pl (<0: free edge)
cc     angle    : angle of 3 element nodes
c
c      
cc---- yltdata5 : genetic data
cc     npop     : number of individuals
cc     nprp     : number of properties (=ndf)
cc     nrec     : number of genes (=nkn)   
cc     ncrt     : number of evaluation criterions
cc     imut     : mutation switch
cc     irec     : recombination switch  
c
c
cc---- cdat1    :
cc     ndd      :
cc     nie      :
c
c
cc---- cdata    : 
cc     numnp    : number of mesh nodes
cc     numel    : number of mesh elements
cc     nummat   : number of material sets
cc     nen      : maximum nodes/element
cc     neq      : number active equations
cc     ipr      : real variable precision  
c
c      
cc---- eldata     : element data
cc     dm         : element proportional load
cc     n          : current element number
cc     ma         : current element material set
cc     mct        : print counter
cc     iel        : user element number
cc     nel        : number nodes on current element  
c
c
cc---- sdata      :
cc     ndf        :
cc     ndm        :
cc     nen1       :
cc     nst        :
c
c      integer ix(6,numel)
c      real*8 pop(ndm,nkn,npop),par(ndm,nkn,2)
c      
cc---- pop        : existing population
cc---- par        : new selected pair
c
c
cc---- Variables --------------------------------------------------------
c
cc==== PROGRAM ==========================================================
c
cc==== Collect Normal Vectors Of Elements ===============================
c
c      do i1=1,numel
c        
c
c      enddo
c
c
cc==== End of Subroutine ================================================
c      
c 9999 continue
c      return
c      end   


c=======================================================================      

      subroutine yltval(d,u,x,ix,lm,id,f,s,p,atb,etb,ftb,ptb,ctb,mtb,
     &                  vtb,lth,pop,vsor,npop,ncrt)
      
c----------------------------------------------------------------------+
c     Evaluate Yield-Line Mechanism                                    |
c----------------------------------------------------------------------+

c==== DECLARATION ======================================================

c---- Formal Parameters ------------------------------------------------
      USE cdat1
      USE cdata
      USE eldata
      USE mdat2
      USE psize
      USE sdata
      USE yltdata1
      USE yltdata2
c      USE yltdata5

c---- mdat2      : rotation data
c     ia         :
c     itrot      : 
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

c---- yltdata5 : genetic data
c     npop     : number of individuals
c     nprp     : number of properties (=ndf)
c     nrec     : number of genes (=nkn)   
c     ncrt     : number of evaluation criterions
c     imut     : mutation switch
c     irec     : recombination switch  

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

      integer id(ndf,nkn), ix(6,numel),lm(2,nth)
      real*8 d(*),u(ndf,*),x(ndm,numnp)
      real*8 f(ndf,numnp),s(nst,nst),p(nst) 
      real*8 atb(1+nth,2*nkn+2*nth),etb(nth,nkn),ftb(nkn)
      real*8 ptb(1+nth),ctb(2*nkn+2*nth),mtb(2*nth)
      real*8 vtb(2*nkn+2*nth)
      real*8 lth(3,nth),pop(ndm,nkn,npop),vsor(1+ncrt,npop)

c---- atb      : table of constrains (th+th-w+w-,1+th) A
c---- etb      : assembled element matrix E
c---- ftb      : assembled element vector F 
c---- ptb      : right hand side
c---- ctb      : coefficients of objective function c
c---- mtb      : assembled element moment vector m
c---- vtb      : optimization variables (th+th-w+w-)
c---- lth      : length of yield-lines (x,y,sqrt(x^2+y^2))
c---- pop      : population pop(nprp,nrec,npop)
c---- vsor     : sorting vector

c---- Variables --------------------------------------------------------

      real*8 th(nth),vevl(ncrt,npop),que,qui,tol
      real*8 pen1,pen2 !,pen3
      real*8 v1(3),v2(3),v3(3),vn1(3),vn2(3),vn3(3) !,lx,ly
      
c---- vevl     : evaluation vector      
c---- th       : rotations      
c---- que      : external work
c---- qui      : internal work
c---- tol      : tolerance
c---- pen      : penalties
c---- v        : vector
c---- vn       : normal vector
c---- lx       : length x
c---- ly       : length y

      common m(maxm)

c==== PROGRAM ==========================================================

c---- Set Tolerance ----------------------------------------------------

      tol=1.d-8

c==== Evaluate Individuals: 1st Crit ===================================    
 
c      call pzero(vevl,ncrt*npop)
c      call pzero(vsor,(ncrt+1)*npop)
      vevl=0
      vsor=0


c      call plotmatrix(vevl,npop,4,'vevl0 von yltval')  

      do i1=1,npop
        
c        call pzero(vtb,2*nth+2*nkn)
c        call pzero(th,nth)
        vtb=0
        th=0
        que=0.d0
        qui=0.d0
        pen1=0.d0
        pen2=0.d0
        
        do i2=1,nkn
          do i3=1,ndm
            x(i3,i2)=pop(i3,i2,i1) !coordinates
          enddo
        enddo      
        
c---- Normalize Displacement -------------------------------------------

        a=0.d0
        do i2=1,nkn
          if (dabs(x(3,i2)).gt.a) a=dabs(x(3,i2))
        enddo
        if (a.gt.0.d0.or.a.lt.0.d0) then
          do i2=1,nkn
            x(3,i2)=x(3,i2)/a
            pop(3,i2,i1)=x(3,i2)
          enddo
        endif                      
        
c---- Calculate Displacements w and Rotations theta --------------------        
        
        call yltass(d,u,x,ix,lm,id,f,s,p,bang,atb,etb,ftb,ptb,ctb,
     &              mtb,lth,30)

        do i2=1,nth
          do i3=1,nkn
            th(i2)=th(i2)+etb(i2,i3)*x(3,i3)
          enddo
        enddo
        do i2=1,nth 
          if (th(i2).gt.0.d0) then
            vtb(i2)    =th(i2)
          else                 
            vtb(nth+i2)=dabs(th(i2))
          endif
        enddo
        do i2=1,nkn
          if (x(3,i2).gt.0.d0) then
            vtb(2*nth+i2)    =x(3,i2)
          else                 
            vtb(2*nth+nkn+i2)=(x(3,i2))  
          endif
        enddo
        
c---- Calculate Internal Work ------------------------------------------        

        call yltevl(vtb,mtb,u,th,qui,1)

c---- Calculate External Work ------------------------------------------
                
        do i2=1,nkn
          que=que+x(3,i2)*ftb(i2)
        enddo
        
c---- Evaluate ---------------------------------------------------------

        if (que.gt.0.d0.or.que.lt.0.d0) vevl(1,i1)=qui/que        
        
c==== Calculate Internal Work Malus: 2nd Crit --------------------------

c  220   continue
        do i2=1,nth
          v1=0.d0
          v2=0.d0
          v3=0.d0
          vn1=0.d0
          vn2=0.d0
          vn3=0.d0
          n1=abs(lm(1,i2))
          n2=abs(lm(2,i2))
          if (n1.eq.0.or.n2.eq.0) exit 
c          le=dsqrt((x(1,n1)-x(1,n2))**2+(x(2,n1)-x(2,n2))**2)
c     &            +(x(3,n1)-x(3,n2))**2)
c          lx=dabs(x(1,n1)-x(1,n2))
c          ly=dabs(x(2,n1)-x(2,n2))
          
          do i3=1,numel
           
            do ia1=1,3
              ib1=ia1+1
              ic1=ia1+2
              if (ib1.gt.3) ib1=ib1-3
              if (ic1.gt.3) ic1=ic1-3
              
              if ((ix(ia1,i3).eq.n1.and.ix(ib1,i3).eq.n2).or.
     &            (ix(ia1,i3).eq.n2.and.ix(ib1,i3).eq.n1)) then
     
c                loc1=1
c                if (ix(ia1,i3).eq.n2.and.ix(ib1,i3).eq.n1) loc1=-1
                
                n3=ix(ic1,i3)
     
                do i4=i3+1,numel
                  do ia2=1,3
                    ib2=ia2+1                   
                    ic2=ia2+2              
                    if (ib2.gt.3) ib2=ib2-3
                    if (ic2.gt.3) ic2=ic2-3

                    if ((ix(ia2,i4).eq.n1.and.ix(ib2,i4).eq.n2).or.
     &                  (ix(ia2,i4).eq.n2.and.ix(ib2,i4).eq.n1)) then
                
c                      loc2=1
c                      if (ix(ia2,i4).eq.n2.and.ix(ib2,i4).eq.n1) loc2=-1

                      n4=ix(ic2,i4)
                      
                      do i5=1,3
                        v1(i5)=x(i5,n2)-x(i5,n1)
                        v2(i5)=x(i5,n3)-x(i5,n1)
                        v3(i5)=x(i5,n4)-x(i5,n1)
                      enddo
                      call vnorm(v1,a)
                      call vnorm(v2,a)
                      call vnorm(v3,a)
                      goto 230
                    endif
                  enddo
                enddo
              endif
            enddo
          enddo

  230     continue          
          call vcross(v2,v1,vn1)
          call vcross(v1,v3,vn2)
          call vnorm(vn1,a)
          call vnorm(vn2,a)
          if (vn1(3).lt.0.d0) vn1=-vn1
          if (vn2(3).lt.0.d0) vn2=-vn2
          call vcross(vn1,vn2,vn3)
          
c          if ((vn1(1).eq.0.d0.and.vn1(2).eq.0.d0.and.vn1(3).eq.0.d0).or.
c     &    (vn2(1).eq.0.d0.and.vn2(2).eq.0.d0.and.vn2(3).eq.0.d0)) cycle    
                                  
c          th(i2)=dacos((vn1(1)*vn2(1)+vn1(2)*vn2(2)+vn1(3)*vn2(3))/        
c     &           (dsqrt(vn1(1)**2+vn1(2)**2+vn1(3)**2)*
c     &            dsqrt(vn2(1)**2+vn2(2)**2+vn2(3)**2)))

          call vnorm(vn3,th(i2))
     
c             bux=th(i2)
          
c          call vnorm(vn3,a)
          
          if (((vn3(1)-v1(1)).lt.tol).and.
     &        ((vn3(2)-v1(2)).lt.tol).and.
     &        ((vn3(3)-v1(3)).lt.tol)) then !not negative yield-line or none

            pen1=pen1+dabs(mtb(i2)*th(i2)**4)
     
c            if (lm(2,i2).lt.0.d0) then
c              qui=qui+dsqrt(lx**2+ly**2)*dlm(4,i2)
c            else
c              qui=qui+dsqrt((lx*dabs(d(5)))**2+(ly*dabs(d(6)))**2)
c     &           *dabs(th(i2))      !BAUSTELLE: Zerlegung th?
c            endif
                                                    
          elseif (((vn3(1)+v1(1)).lt.tol).and.
     &            ((vn3(2)+v1(2)).lt.tol).and.
     &            ((vn3(3)+v1(3)).lt.tol)) then !not positive yield-line 

            pen1=pen1+dabs(mtb(nth+i2)*th(i2)**4)
                            
c            if (lm(2,i2).lt.0.d0) then
c              qui=qui+dsqrt(lx**2+ly**2)*dlm(4,i2)
c            else
c              qui=qui+dsqrt((lx*dabs(d(3)))**2+(ly*dabs(d(4)))**2)
c     &         *dabs(th(i2))
c            endif
                      
          endif
        enddo
        
c---- Evaluate ---------------------------------------------------------

        vevl(2,i1)=pen1        

c==== Punish Flat Elements: 3rd Crit ===================================

        do i2=1,numel
          do ia1=1,3
            ib1=ia1+1
            ic1=ia1+2
            if (ib1.gt.3) ib1=ib1-3
            if (ic1.gt.3) ic1=ic1-3
            
            do i3=1,3
              v1(i3)=x(i3,ix(ib1,i2))-x(i3,ix(ia1,i2))
              v2(i3)=x(i3,ix(ic1,i2))-x(i3,ix(ia1,i2))
            enddo
            
            call vcross(v1,v2,vn1)
            call vnorm(vn1,a)
            
            pen2=pen2+vn1(3)
          enddo
        enddo    

c---- Evaluate ---------------------------------------------------------

        vevl(3,i1)=pen2        
      enddo

cc==== Final Taxation: 4th Crit =========================================
c
c      do i1=1,npop
c        vevl(4,i1)=0.d0*vevl(1,i1)+1.d0*vevl(2,i1)+0.d0*vevl(3,i1)
c        
c      print*, vevl(2,i1),vevl(4,i1)
c      enddo      
      
c---- Normalize Criterions ---------------------------------------------

      a1=0.d0
      a2=0.d0
      a3=0.d0
      do i1=1,npop
        if (vevl(2,i1).gt.a1) a1=vevl(2,i1)
        if (vevl(3,i1).gt.a2) a2=vevl(3,i1)
        if (vevl(4,i1).gt.a3) a3=vevl(4,i1)
      enddo
c      if (a.ne.0.d0) then
        do i1=1,npop
          if (a1.gt.0.d0.or.a1.lt.0.d0) vevl(2,i1)=vevl(2,i1)/a1
          if (a2.gt.0.d0.or.a2.lt.0.d0) vevl(3,i1)=vevl(3,i1)/a2
          if (a3.gt.0.d0.or.a3.lt.0.d0) vevl(4,i1)=vevl(4,i1)/a3
        enddo

c---- Sort Individuals for Criterions ----------------------------------
  
      n1=1
      do i1=1,npop
        a=0.d0
        do i2=1,npop
          if (vevl(2,i2).gt.a) a=vevl(2,i2)
        enddo  
c        b=max(vevl(2,1:npop))
        do i2=1,npop
          if (vevl(2,i2).ge.a.and.vevl(2,i2).le.a) then
            vsor(1,n1)=i2
            vsor(2,n1)=vevl(1,i2)
            vsor(3,n1)=vevl(2,i2)
            vsor(4,n1)=vevl(3,i2)
            vsor(5,n1)=vevl(4,i2)
            vevl(1,i2)=0.d0
            vevl(2,i2)=0.d0
            n1=n1+1
            exit
          endif
        enddo
      enddo      

c      call plotmatrix(vsor,5,npop,'vsor von yltval')  
              

c==== End of Subroutine ================================================
      
c 999  continue
      return
      end   
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
      
