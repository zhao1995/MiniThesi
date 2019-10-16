 
c----------------------------------------------------------------------+
c                                                                      |
c     Special Subroutines for Yield-Line Element  jw071104             |
c                                                                      |
c     Part 2: Gradient Calculation, Mesh Adjustment and Parameters     |
c                                                                      |
c----------------------------------------------------------------------+


c=======================================================================      

      subroutine yltgrd(d,u,x,ix,lm,id,f,s,p,ang,atb,etb,ftb,ptb,ctb,
     &                  mtb,vtb,a1,c1y,y1,da1,dc1,dqu,dquh,lth,dtb,isw)
      
c----------------------------------------------------------------------+
c     Perform Gradient Calculation                                     |
c----------------------------------------------------------------------+

c==== DECLARATION ======================================================

c---- Formal Parameters ------------------------------------------------

      USE cdat1
      USE cdata
      USE eldata
      USE iofile
c      USE mdat2
      USE sdata
      USE tdata
      USE yltdata1
      USE yltdata2
      implicit real*8 (a-h,o-z)
      implicit integer (i,n)

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

c---- mdat2      : rotation data
c     ia         :
c     itrot      :  

c---- sdata      :
c     ndf        :
c     ndm        :
c     nen1       :
c     nst        :

c---- tdata    : time steps
c     dt       : step length

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

      integer id(ndf,nkn), ix(6,numel),lm(2,nth)
      real*8 d(ndd,*),u(ndf,*),x(ndm,numnp)
      real*8 f(ndf,numnp),s(nst,nst),p(nst),ang(numnp) 
      real*8 atb(1+nth,2*nkn+2*nth),etb(nth,nkn),ftb(nkn)
      real*8 ptb(1+nth),ctb(2*nkn+2*nth),mtb(2*nth)
      real*8 vtb(2*nkn+2*nth),dtb(ny1,ndx)
      real*8 lth(3,nth)
c      real*8 a1(nth+1,nth+1),c1y(nth+1),y1(nth+1)
      real*8 a1(ny1,ny1),c1y(ny1),y1(ny1)

c      integer iid(ndx)
      real*8 da1(ny1,ny1,ndx),dc1(ny1,ndx),dqu(ndx),dquh(ndx)

c---- ang      : angle of nodes
c---- lmy      : location matrix A_1 entries (rows,columns)
c---- atb      : table of constrains (th+th-w+w-,1+th) A
c---- ctb      : coefficients of objective function c
c---- mtb      : assembled element moment vector m
c---- vtb      : optimization variables (th+th-w+w-)
c---- dtb      : D-matrix of degeneration
c---- a1       : A_1 matrix
c---- c1y      : c_1 vector
c---- y1       : y_1 vector
c---- da1      : derivation array of A_1: dA_1/dX for 1 to ndx
c---- dc1      : derivation array of c_1
c---- dqu      : gradient array (dL(X)/dX)
c---- dquh     : previous gradient

      integer lmy(2,ny1),irg,isort(1+nth) !,isort2(1+nth),isgn
      real*8 atmp(1+nth,2*nth+2*nkn),vtmp(2*nth+2*nkn)!,ctmp(2*nth+2*nkn)
      real*8 a1i(ny1,ny1),dca(ny1,ny1),dcc(ny1)

c---- lmy      : location matrix for A_1 (rows,columns)
c---- irg      : rank
c---- isort    : sort vector for rows in ygauss
c---- isgn     : sign of Dg
c---- atmp     : auxiliary array to determine columns of A_1
c---- a1i      : inverse A_1
c---- dca      : auxiliary array for gradient calculation
c---- dcc      : auxiliary vector

c      logical lgrd
      
c---- lgrd     : exit if singular matrix A_1      

c==== PROGRAM ==========================================================

c---- Transmit former dqu to dquh --------------------------------------

      call pzero(dquh,ndx)
      if (iter.gt.1) then
        do i1=1,ndx
          dquh(i1)=dqu(i1)
        enddo  
      endif      

c---- Initialization ---------------------------------------------------

      call pzero(atmp,(2*nth+2*nkn)*(nth+1))
      call pzero(vtmp,2*nth+2*nkn)
      call pzero(a1,(nth+1)*(nth+1))
      call pzero(c1y,(nth+1))
c      call pzero(lmy,2*ny1)
      lmy=0
      call pzero(da1,ny1*ny1*ndx)
      call pzero(dc1,ny1*ndx)
      call pzero(dqu,ndx)
      call pzero(a1i,ny1*ny1)
      call pzero(dtb,ny1*ndx)
      if (iter.eq.1) dt=0.d0 !allows to define step length in yltupd
      
c---- Rearrange Columns Using vtb --------------------------------------

      i3=0
      do i1=1,2*nth+2*nkn
        if (vtb(i1).gt.0.d0.or.vtb(i1).lt.0.d0) then
          i3=i3+1
          do i2=1,nth+1
            atmp(i2,i3)=atb(i2,i1)
          enddo
          vtmp(i3)=vtb(i1)
c          ctmp(i3)=ctb(i1)
        endif
      enddo
      
c      call plotmatrix(vtb,2*nkn+2*nth,1,'vtb')
c      call plotmatrix(vtmp,2*nkn+2*nth,1,'vtmp')
c      call plotmatrix(atb ,1+nth,2*nkn+2*nth,'atb')
c      call plotmatrix(atmp,1+nth,2*nkn+2*nth,'atmp')            
c      call plotmatrix(ctb,2*nkn+2*nth,1,'ctb')
c      call plotmatrix(ctmp,2*nkn+2*nth,1,'ctmp')

c---- Bring atb Into Step Form -----------------------------------------

      call ygauss(atmp,1+nth,2*nkn+2*nth,isort,irg)
      
cc---- Exit if A_1 Singular ---------------------------------------------
c
c      if (irg.lt.ny1) then
c
cc        lgrd=.false.
c        goto 9999
c      
c      endif  

c---- Form Location Matrix lmy -----------------------------------------

      i2=0
      do i1=1,2*nth+2*nkn
        if (vtb(i1).gt.0.d0.or.vtb(i1).lt.0.d0) then
          i2=i2+1
          lmy(2,i2)=i1  !columns
        endif
      enddo  

      do i1=ny1,1,-1
        i3=0
        do i2=1,ny1
          if (isort(i2).gt.i3) then
            i3=isort(i2)
            i4=i2
          endif  
        enddo  
        lmy(1,i1)=i3
        isort(i4)=0
      enddo    
      
c---- Compose A_1, c_1, y_1 --------------------------------------------

      do i2=1,ny1
        do i1=1,ny1
          a1(i1,i2)=atb(lmy(1,i1),lmy(2,i2))
          a1i(i1,i2)=atb(lmy(1,i1),lmy(2,i2))
        enddo
        c1y(i2)=ctb(lmy(2,i2))
        y1(i2)=vtb(lmy(2,i2))        
      enddo
      
c---- Last Singularity Check -------------------------------------------

c      call plotmatrix(a1,ny1,ny1,'a1')
      a1i=a1
      call ygauss(a1i,ny1,ny1,isort,irg)
c      call plotmatrix(a1i,ny1,ny1,'a1nachgauss')
c      do i1=1,ny1
c        print*,'isort',isort(i1)
c      enddo  

c      a1i=a1
c      call gauss(ny1,ny1,a1i,dcc,det,ndet,eps3,ierr)
c
c      print*, 'det',det*10**ndet,det,ndet
c      if (dabs(det*10**ndet).lt.5.d-18) then
c                     write(iow,2010)  !BAUSTELLE: nopr-Option einbauen
c        if(ior.lt.0) write(*  ,2010) 
c        stop
c
c        dqu=dquh
c        goto 9999        
c      endif

c---- Calculate Inverse of A_1 -----------------------------------------

c      print*, 'irg', irg,ny1

c      a1i=a1
      call pivot(a1,ny1,ny1,a1i)
c      call invert(a1i,ny1,ny1)

c==== Preliminary Checks to Avoid Singularity of A_1 ===================

c---- Calculate Derivations dA1 ----------------------------------------

      i3=0
      do j1=1,nkn
        do j2=1,2
          if (id(j2,j1).gt.0) then
            i3=i3+1
            idx(1)=j1
            idx(2)=j2 

            call yltass(d,u,x,ix,lm,id,f,s,p,atb,etb,ftb,ptb,ctb,
     &                  mtb,lth,isw)
                                          

c---- Assemble dA_1, dc_1 and dy_1 -------------------------------------

            do i2=1,ny1
              do i1=1,ny1
                da1(i1,i2,i3)=atb(lmy(1,i1),lmy(2,i2)) 
              enddo
              dc1(i2,i3)=ctb(lmy(2,i2))
            enddo 
          endif

c---- Second Part of Rotated Angles ------------------------------------

          if (j2.eq.2.and.(ang(j1).gt.0.d0.or.ang(j1).lt.0.d0)) then

            idx(1)=j1
            idx(2)=j2 

            call yltass(d,u,x,ix,lm,id,f,s,p,atb,etb,ftb,ptb,ctb,
     &                  mtb,lth,isw)

            do i2=1,ny1
              do i1=1,ny1
                da1(i1,i2,i3)=da1(i1,i2,i3)*dcos(dacos(0.d0)
     &                                           *ang(j1)/9.d1)
     &          +atb(lmy(1,i1),lmy(2,i2))*dsin(dacos(0.d0)*ang(j1)/9.d1) 
              enddo
              dc1(i2,i3)=dc1(i2,i3)*dcos(dacos(0.d0)*ang(j1)/9.d1)
     &                  +ctb(lmy(2,i2))*dsin(dacos(0.d0)*ang(j1)/9.d1)
            enddo 

          endif
        enddo         
      enddo
      
c---- Calculate Gradient -----------------------------------------------
 
      do i1=1,ndx 
        call pzero(dca,ny1*ny1)
        call pzero(dcc,ny1)
        do i2=1,ny1
          do i3=1,ny1
            do i4=1,ny1
              dca(i2,i3)=dca(i2,i3)+a1i(i2,i4)*da1(i4,i3,i1)
            enddo
          enddo
        enddo
        do i2=1,ny1
          do i3=1,ny1
            dcc(i2)=dcc(i2)+c1y(i3)*dca(i3,i2)  
          enddo
        enddo
        do i2=1,ny1 ! D-Matrix dtb
          do i3=1,ny1
            dtb(i3,i1)=dtb(i3,i1)+dca(i3,i2)*y1(i2)
          enddo
        enddo
        do i2=1,ny1
          dqu(i1)=dqu(i1)+(dc1(i2,i1)-dcc(i2))*y1(i2)
        enddo
      enddo          
  
c      call plotmatrix(dqu,1,ndx)  
c      print*, 'lmy', lmy(1,3),!lmy(2,3)
c0601      call plotmatrix(dtb,ny1,ndx,'dtb') !lmy(2,3) für D ithav_96b, 27.3.05

c==== Formats ==========================================================

c 2010 format(/'  *** Warning: Determinant Zero in yltdeg, Change Step')  
                 
c==== End of Subroutine ================================================
      
c 9999 continue
      return
      end



c=======================================================================      
c=======================================================================      

      subroutine yltupd(d,u,x,ix,lm,id,f,s,p,ang,atb,etb,ftb,ptb,ctb,
     &                  mtb,vtb,btb,th,a1,c1y,y1,da1,dc1,dqu,lth,qu,
     &                  alp,dtb,vtbh,dquh) 

c----------------------------------------------------------------------+
c     Mesh Update Algorithm                                            |
c     McKeown's Algorithm                                              |
c----------------------------------------------------------------------+

c==== DECLARATION ======================================================

c---- Formal Parameters ------------------------------------------------

      USE cdat1
      USE cdata
      USE iofile
c      USE mdat2
      USE sdata
      USE tdata
      USE yltdata1
      USE yltdata3
      USE yltdata4
      implicit real*8 (a-h,o-z)
      implicit integer (i,n)

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

c---- sdata      :
c     ndf        :
c     ndm        :
c     nen1       :
c     nst        :

c---- mdat2      : rotation data
c     ia         :
c     itrot      :  

c---- tdata    : time steps
c     dt       : step length

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

      integer id(ndf,nkn),ix(6,numel),lm(2,nth) !,lmy(2,ny1)
      integer ideg(ny1),ndeg,ivtb(2*nkn+2*nth)

c---- ideg     : degeneration vector (entry=1=degenerated) 
c---- ndeg     : number of degenerated thetas  
c---- ivtb     : degeneration vector vtb
      
c      real*8 a1(ny1,ny1),c1y(ny1),y1(ny1),da1(ny1,ny1,ndx),dc1(ny1,ndx)    
c      real*8 x(ndm,numnp),xold(ndm,numnp),dqu(ndx),dp(ndx),sca,alp
c      
cc---- a1       : A_1 matrix
cc---- c1y      : c_1 vector
cc---- y1       : y_1 vector
cc---- da1      : derivation array of A_1: dA_1/dX for 1 to ndx
cc---- dc1      : derivation array of c_1
 
      real*8 d(ndd,*),u(ndf,*),x(ndm,numnp)
      real*8 f(ndf,numnp),s(nst,nst),p(nst),ang(numnp) 
      real*8 atb(1+nth,2*nkn+2*nth),etb(nth,nkn),ftb(nkn)
      real*8 ptb(1+nth),ctb(2*nkn+2*nth),mtb(2*nth)
      real*8 vtb(2*nkn+2*nth),dtb(ny1,ndx),vtbh(2*nkn+2*nth)
      real*8 lth(3,nth),th(nth)
c      real*8 a1(nth+1,nth+1),c1y(nth+1),y1(nth+1)
      real*8 a1(ny1,ny1),c1y(ny1),y1(ny1)

c      integer iid(ndx)
      real*8 da1(ny1,ny1,ndx),dc1(ny1,ndx),dqu(ndx),dquh(ndx)
c      real*8 dca(ny1,ny1),dcc(ny1)
c      real*8 a1i(ny1,ny1)
c      real*8 da(ny1,ny1),dc(ny1)
      real*8 xh(ndm,numnp),dp(ndx),sca,alp,bet
      real*8 dqu2(ndx)
c      real*8 a12(ny1,ny1),c12(ny1),y12(ny1)
c      real*8 da12(ny1,ny1,ndx),dc12(ny1,ndx),dtb2(ny1,ndx)
c      real*8 qaux(2),degth(nth),y1s(ny1)
      real*8 xs(ndm,numnp),qus

c---- lmy      : location matrix A_1 entries (rows,columns)
c---- ang      : angle of nodes
c---- atb      : table of constrains (th+th-w+w-,1+th) A
c---- ctb      : coefficients of objective function c
c---- mtb      : assembled element moment vector m
c---- vtb      : optimization variables (th+th-w+w-)
c---- vtbh     : previous vtb
c---- dtb      : complete D-matrix of degeneration
c---- a1       : A_1 matrix
c---- c1y      : c_1 vector
c---- y1       : y_1 vector
c---- da1      : derivation array of A_1: dA_1/dX for 1 to ndx
c---- dc1      : derivation array of c_1
c---- dca      : auxiliary array for gradient calculation
c---- dcc      : auxiliary vector
c---- dqu      : gradient array (dL(X)/dX)
c---- dquh     : previous gradient
c---- a1i      : inverse of A_1
c---- da       : auxiliary matrix for dA
c---- dc       : auxiliary vector for dc
c---- xh       : previous x
c---- dp       : normalized gradient
c---- sca      : scaling factor
c---- alp      : alpha value
c---- bet      : beta value
c---- dqu2     : gradient g' for degeneration check
c---- a12..dtb2: dummies for call yltgrd
c---- qaux     : auxiliary field for qu
c---- degth    : save theta-vector for degeneration
c---- xs       : save x-field
c---- qus      : save qu
c---- y1s      : save y1s

      integer bas(2*nth+2*nkn,2)  !,ideg(ny1),ndeg
      logical basis,degen !,lgrd

c---- bas      : vector for basis check
c---- ideg     : degeneration vector
c---- ndeg     : number of degeneration values       
c---- basis    : change in basis
c---- degen    : degeneration      
c---- lgrd     : successful gradient calculation


c==== PROGRAM ==========================================================

c==== Save Incoming Values =============================================

      xs=x
      xh=x
      qus=qu
      alph=qu
      dqu2=dquh
c      dquh=dqu
c      y1s=y1
      vtbh=vtb
      
c==== Initialization ===================================================

      degen=.false.      
      
c==== Check 1st Stop Criterion ========================================= c)

      if (iter.eq.1) then
        sca=0.d0
        do i1=1,ndx
          sca=sca+dqu(i1)**2
        enddo
        sca=dsqrt(sca)
        
        if (sca.lt.eps1) then 
          halt=.true. 
          if (ihalt.eq.0) ihalt=iter
                       write(iow,2010) ihalt  !BAUSTELLE: nopr-Option einbauen
          if(ior.lt.0) write(*  ,2010) ihalt 
          goto 9999
        endif

c==== Set Gradient ===================================================== d)

        dp=-dqu
      endif

c==== Compute Step Length alpha ======================================== e)

c   50 continue

c---- Calculate Step Length or Take User Defined One -------------------

      if (dt.ge.0.d0.and.dt.le.0.d0) then

        call yltalp(d,u,xh,ix,lm,id,f,s,p,ang,atb,etb,ftb,ptb,ctb,
     &              mtb,vtb,btb,th,dqu,lth,alp) 
c        if (iter.eq.1.and.alp.eq.0.d0) alp=1.d0           
      else
        alp=dt
      endif

c---- Update Mesh ------------------------------------------------------

      call yltmsh(x,ang,id,dqu,alp)
      
c==== Check 2nd Tolerance Criterion ==================================== f)

      sca=0.d0
      do i1=1,ndm
        do i2=1,numnp
          sca=sca+(x(i1,i2)-xh(i1,i2))**2
        enddo
      enddo
      sca=dsqrt(sca)
                          
      if (sca.lt.eps2.and.sca.gt.0.d0.and.iter.gt.1) then  
        halt=.true. 
        if (ihalt.eq.0) ihalt=iter
                     write(iow,2020) ihalt  !BAUSTELLE: nopr-Option einbauen
        if(ior.lt.0) write(*  ,2020) ihalt 
        goto 9999
      endif     

c==== Check Basis and Degeneration ===================================== g)

c---- Check Degeneration First -----------------------------------------

      degen=.false.    
      ideg=0
      ndeg=0
      isgn=0
      if (iter.ne.1) then   
        do i1=1,ny1
          if (dabs(y1(i1)).lt.eps3) then
            degen=.true.
            ideg(i1)=1
            ndeg=ndeg+1
          else
            ideg(i1)=0  
          endif
        enddo
        
        ivtb=0
        do i1=1,2*nth+2*nkn
          if (dabs(vtbh(i1)).lt.eps3.and.
     &       (vtbh(i1).gt.0.d0.or.vtbh(i1).lt.0.d0)) then
            ivtb(i1)=1
          else
            ivtb(i1)=0
          endif
        enddo      
        
        if (degen) then          
                       write(iow,2230)   !BAUSTELLE: nopr-Option einbauen
          if(ior.lt.0) write(*  ,2230)  

          do i2=1,2*nth+2*nkn
            if (ivtb(i2).eq.1) then
              vtb(i2)=eps3
            else
              vtb(i2)=vtbh(i2)
            endif
          enddo

          x=xh
          call yltass(d,u,x,ix,lm,id,f,s,p,atb,etb,ftb,
     &               ptb,ctb,mtb,lth,30)
          
          
          call yltsqp(atb,btb,ctb,ptb,vtb,vtb,0)

          call yltgrd(d,u,x,ix,lm,id,f,s,p,ang,atb,etb,ftb,ptb,ctb,
     &                mtb,vtb,a1,c1y,y1,da1,dc1,dqu,dquh,lth,dtb,31)            

          do i2=1,2*nth+2*nkn
            if (ivtb(i2).eq.1) then
              ivtb(i2)=0
              vtb(i2)=-eps3 !-vtbh(i2)
            else
              vtb(i2)=vtbh(i2)
            endif
          enddo
          
          x=xh
          call yltass(d,u,x,ix,lm,id,f,s,p,atb,etb,ftb,
     &               ptb,ctb,mtb,lth,30)
             
          call yltsqp(atb,btb,ctb,ptb,vtb,vtb,0)

          call yltgrd(d,u,x,ix,lm,id,f,s,p,ang,atb,etb,ftb,ptb,ctb,
     &                mtb,vtb,a1,c1y,y1,da1,dc1,dqu2,dquh,lth,dtb,31)            


          call yltdeg(dtb,dqu,dquh,dp,ideg,ndeg,isgn)
          
          if (isgn.eq.3) then
            do i1=1,ndx
              if (dabs(dp(i1)).lt.eps6) dp(i1)=0.d0
            enddo
          endif    

          do i2=1,2*nth+2*nkn
            if (ivtb(i2).eq.1) then
              ivtb(i2)=0
              vtb(i2)=0.d0 
            else
              vtb(i2)=vtbh(i2)
            endif
          enddo

          goto 300
        endif
      endif   

c---- Check Basis ------------------------------------------------------

c   70 continue
      basis=.false.

      if (iter.gt.1) then
        do i1=1,2*nkn+2*nth
c          vtbh(i1)=vtb(i1)                                                 !!! Luenberger-Check!
          if (vtbh(i1).ge.0.d0.and.vtbh(i1).le.0.d0) then
            bas(i1,1)=0
          else
            bas(i1,1)=1
          endif    
        enddo
              
        i2=0
        do i1=1,2*nth+2*nkn
          if (vtb(i1).ge.0.d0.and.vtb(i1).le.0.d0) then
            bas(i1,2)=0
          else
            bas(i1,2)=1
            i2=i2+1
c            lmy(2,i2)=i1
          endif
        enddo
        
        do i1=1,2*nth
          if ((bas(i1,1)-bas(i1,2)).ne.0) basis=.true.            
        enddo
      endif

      if (basis) then
        goto 200
      else
        goto 100
      endif    

c==== Decision Tree ====================================================

c---- [1] Non-Degenerate and Same Basis --------------------------------

  100 continue      
                     write(iow,2210)   !BAUSTELLE: nopr-Option einbauen
        if(ior.lt.0) write(*  ,2210)  
                
c---- Calculate beta ---------------------------------------------------

        if (iter.gt.1) then
          call yltbet(dqu,dquh,bet)

                       write(iow,2120) bet  !BAUSTELLE: nopr-Option einbauen
          if(ior.lt.0) write(*  ,2120) bet 
        else
          bet=1.d0
        endif    

c---- Calculate dp -----------------------------------------------------

        do i1=1,ndx
          dp(i1)=-dqu(i1)-bet*dquh(i1)
        enddo

c---- Calculate alpha --------------------------------------------------

        if (iter.gt.1.and.(dt.ge.0.d0.and.dt.le.0.d0)) then 
          call yltalp(d,u,xh,ix,lm,id,f,s,p,ang,atb,etb,ftb,ptb,ctb,
     &              mtb,vtb,btb,th,dp,lth,alp)        
        endif
        
                     write(iow,2110) alp  !BAUSTELLE: nopr-Option einbauen
        if(ior.lt.0) write(*  ,2110) alp 

c---- Normalized Gradient dp -------------------------------------------

        sca=0.d0
        do i1=1,ndx
          if (dabs(dp(i1)).gt.dabs(sca)) sca=dp(i1)         
        enddo
        do i1=1,ndx
          dp(i1)=dp(i1)/sca  
        enddo
 
c---- Update Mesh ------------------------------------------------------

        call yltmsh(x,ang,id,dp,alp)          
        xh=x        
        alph=qu
        dquh=dqu

        goto 500

cc---- [1] Non-Degenerate and Same Basis --------------------------------
c
c  101 continue
c                   write(iow,2210)   !BAUSTELLE: nopr-Option einbauen
c      if(ior.lt.0) write(*  ,2210)  
c                
cc---- Calculate beta ---------------------------------------------------
c
c      if (iter.gt.1) then
c        call yltbet(dqu,dquh,bet)
c
c                     write(iow,2120) bet  !BAUSTELLE: nopr-Option einbauen
c        if(ior.lt.0) write(*  ,2120) bet 
c      else
c        bet=1.d0
c      endif    
c
cc---- Calculate dp -----------------------------------------------------
c
c
c      do i1=1,ndx !hier noch dqu
c        dp(i1)=-dqu(i1)-bet*dquh(i1)!/sca
c      enddo
c
c      sca=0.d0  !...hier aber normalisiert im Ggs. zu [2]
c      do i1=1,ndx
c        if (dabs(dp(i1)).gt.dabs(sca)) sca=dp(i1)
c      enddo  
c      if (sca.eq.0.d0) sca=1.d0  
c
c      do i1=1,ndx
c        dp(i1)=dp(i1)/sca
c      enddo  
c
cc---- Calculate alpha --------------------------------------------------
c
c      if (iter.gt.1.and.dt.eq.0.d0) then 
c        call yltalp(d,u,xh,ix,lm,id,f,s,p,ang,atb,etb,ftb,ptb,ctb,         !Passt!
c     &              mtb,vtb,btb,th,dp,lth,alp)        
c      endif
c                   write(iow,2110) alp  !BAUSTELLE: nopr-Option einbauen
c      if(ior.lt.0) write(*  ,2110) alp 
c
c      if (iter.gt.1) then !geeicht nach ithav74, 20.2.06
c        x=xh
c        call yltmsh(x,ang,id,dp,alp)
c      endif 
c 
c      goto 500 

c---- [2] Non-Degenerate and Different Basis --------------------------- h)

  200 continue          
                   write(iow,2220)   !BAUSTELLE: nopr-Option einbauen
      if(ior.lt.0) write(*  ,2220)  

c---- Calculate beta ---------------------------------------------------

      if (iter.gt.1) then
        call yltbet(dqu,dquh,bet)

                     write(iow,2120) bet  !BAUSTELLE: nopr-Option einbauen
        if(ior.lt.0) write(*  ,2120) bet 
      else
        bet=1.d0
      endif    

c---- Calculate dp -----------------------------------------------------

      do i1=1,ndx !hier noch dqu
        dp(i1)=-dqu(i1)-bet*dquh(i1)
      enddo

ccjw>> einkommentiert: McKeown-Algorithmus. Offensichtlich aber anders...
cc---- Calculate dp -----------------------------------------------------
c
c      dp=-dqu
c 
ccjw<<  keine Normalisierung, DAS hat McKeown gemeint!!!!

 
c---- Calculate alpha --------------------------------------------------

      if (iter.gt.1.and.(dt.ge.0.d0.and.dt.le.0.d0)) then 
        call yltalp(d,u,xh,ix,lm,id,f,s,p,ang,atb,etb,ftb,ptb,ctb,        
     &              mtb,vtb,btb,th,dp,lth,alp)        
      endif

                   write(iow,2110) alp  !BAUSTELLE: nopr-Option einbauen
      if(ior.lt.0) write(*  ,2110) alp 

c---- Normalized Gradient dp -------------------------------------------

      sca=0.d0
      do i1=1,ndx
        if (dabs(dp(i1)).gt.dabs(sca)) sca=dp(i1)         
      enddo
      do i1=1,ndx
        dp(i1)=dp(i1)/sca  
      enddo
 
      if (iter.gt.1) then !geeicht nach ithav74, 20.2.06
        x=xh
        call yltmsh(x,ang,id,dp,alp)
      endif 

      goto 500
        
c---- [3] Degenerate --------------------------------------------------- i)

 300  continue
        
c---- Dg < 0 -----------------------------------------------------------

      if (isgn.eq.1) then !analog [2]

        if (iter.gt.1) then
          call yltbet(dqu,dquh,bet)
        
                       write(iow,2120) bet  !BAUSTELLE: nopr-Option einbauen
          if(ior.lt.0) write(*  ,2120) bet 
        else
          bet=1.d0
        endif    
                
        do i1=1,ndx !hier noch dqu
          dp(i1)=-dqu(i1)-bet*dquh(i1)
        enddo
                
        if (iter.gt.1.and.(dt.ge.0.d0.and.dt.le.0.d0)) then 
          call yltalp(d,u,xh,ix,lm,id,f,s,p,ang,atb,etb,ftb,ptb,ctb,        
     &                mtb,vtb,btb,th,dp,lth,alp)        
        endif
        
                     write(iow,2110) alp  !BAUSTELLE: nopr-Option einbauen
        if(ior.lt.0) write(*  ,2110) alp 
                
        sca=0.d0
        do i1=1,ndx
          if (dabs(dp(i1)).gt.dabs(sca)) sca=dp(i1)         
        enddo
        do i1=1,ndx
          dp(i1)=dp(i1)/sca  
        enddo
        
        if (iter.gt.1) then 
          x=xh
          call yltmsh(x,ang,id,dp,alp)
        endif 
        goto 500

c---- Dg' < 0 ----------------------------------------------------------

      elseif (isgn.eq.2) then

        if (iter.gt.1) then
          call yltbet(dqu,dquh,bet)
        
                       write(iow,2120) bet  !BAUSTELLE: nopr-Option einbauen
          if(ior.lt.0) write(*  ,2120) bet 
        else
          bet=1.d0
        endif    
         
        dqu=-dqu !nach McKeown-Buch, S.65 könnte dies g' sein        
        do i1=1,ndx !hier noch dqu
          dp(i1)=-dqu(i1)-bet*dquh(i1)
        enddo
                
        if (iter.gt.1.and.(dt.ge.0.d0.and.dt.le.0.d0)) then 
          call yltalp(d,u,xh,ix,lm,id,f,s,p,ang,atb,etb,ftb,ptb,ctb,        
     &                mtb,vtb,btb,th,dp,lth,alp)        
        endif
        
                     write(iow,2110) alp  !BAUSTELLE: nopr-Option einbauen
        if(ior.lt.0) write(*  ,2110) alp 
                
        sca=0.d0
        do i1=1,ndx
          if (dabs(dp(i1)).gt.dabs(sca)) sca=dp(i1)         
        enddo
        do i1=1,ndx
          dp(i1)=dp(i1)/sca  
        enddo
        
        if (iter.gt.1) then 
          x=xh
          call yltmsh(x,ang,id,dp,alp)
        endif 
        goto 500

c---- Dg > 0 and Dg' > 0 -----------------------------------------------

      else  

                     write(iow,2235) 
        if(ior.lt.0) write(*  ,2235) 
        do i1=1,ndx        
                       write(iow,2236) i1,dp(i1)  !BAUSTELLE: nopr-Option einbauen
          if(ior.lt.0) write(*  ,2236) i1,dp(i1)
        enddo

        if (iter.gt.1.and.(dt.ge.0.d0.and.dt.le.0.d0)) then 
          call yltalp(d,u,xh,ix,lm,id,f,s,p,ang,atb,etb,ftb,ptb,ctb,        
     &                mtb,vtb,btb,th,dp,lth,alp)        
        endif
        
                     write(iow,2110) alp  !BAUSTELLE: nopr-Option einbauen
        if(ior.lt.0) write(*  ,2110) alp 
 
c        alp=1.d0
        
        sca=0.d0
        do i1=1,ndx
          if (dabs(dp(i1)).gt.dabs(sca)) sca=dp(i1)         
        enddo
        do i1=1,ndx
          dp(i1)=dp(i1)/sca  
        enddo
        
        if (iter.gt.1) then !geeicht nach ithav74, 20.2.06
          x=xh
          call yltmsh(x,ang,id,dp,alp)
        endif 
        
        goto 500
        
      endif  
        
c>>>> nur noch als Archiv; kann weg...
c        if (iter.gt.1.and.dt.eq.0.d0) then 
c          call yltalp(d,u,x,ix,lm,id,f,s,p,ang,atb,etb,ftb,ptb,ctb,        
c     &                mtb,vtb,btb,th,dp,lth,alp)        
c        endif
c        
c                     write(iow,2110) alp  !BAUSTELLE: nopr-Option einbauen
c        if(ior.lt.0) write(*  ,2110) alp 
c        
c        sca=0.d0
c        do i1=1,ndx
c          if (dabs(dp(i1)).gt.dabs(sca)) sca=dp(i1)         
c        enddo
c        do i1=1,ndx
c          dp(i1)=dp(i1)/sca  
c        enddo
c        
c        if (iter.gt.1) then !geeicht nach ithav74, 20.2.06
c          x=xh
c          call yltmsh(x,ang,id,dp,alp)
c        endif 
c        
c        goto 500
c<<<<
        
c==== Definitive Actualization =========================================

  500 continue
      x=xs
c      if (alp.eq.0.d0) alp=eps3
      call yltmsh(x,ang,id,dp,alp)          
      xh=x        
      alph=qu
      dquh=dqu
 
      call yltass(d,u,x,ix,lm,id,f,s,p,atb,etb,ftb,
     &           ptb,ctb,mtb,lth,30)
      
      call yltsqp(atb,btb,ctb,ptb,vtb,vtb,0)
      
      call yltevl(vtb,mtb,u,th,qu,1)    

c---- Fine Justage -----------------------------------------------------

c      if (dabs(qu-qus).lt.eps3) then
c        if (dt.eq.0.d0) then
c          dt=eps3
c        else
c          dt=0.d0
c        endif
c      endif

      if (qus.lt.qu.and.iter.gt.2) then
                     write(iow,2250)  !BAUSTELLE: nopr-Option einbauen
        if(ior.lt.0) write(*  ,2250)
ccjw        if (dt.eq.0.d0) then
ccjw          dt=eps3
ccjw        else
ccjw          dt=0.d0
ccjw        endif
ccjw        print*, dt    
c        x=xs
c        qu=qus
      endif         
      
c==== Formats ==========================================================    

 2010 format(/'  1st tolerance criterion met in iteration no. ',1i3)
 2020 format(/'  2nd tolerance criterion met in iteration no. ',1i3)
c2030 format(/'   3rd tolerance criterion met after ',1i3, ' iterations')
 2110 format(/'  Parameter alpha        = ',1pe15.7)
 2120 format(/'  Parameter beta         = ',1pe15.7)
 2210 format(/'  [1] No Degeneration, no Change in Basis')  
 2220 format(/'  [2] No Degeneration, Change in Basis')
 2230 format(/'  [3] Degeneration')
 2235 format(/'  Scaled Gradient'/
     &'  dof X     dp/dX')
 2236 format(3x,1i3,5x,1pe15.7)
 2250 format(/'  Performing Null-Step')  


c==== End of Subroutine ================================================
      
9999  continue
      return
      end      


c=======================================================================      
c=======================================================================      

      subroutine yltmsh(x,ang,id,dqu,alp)

c----------------------------------------------------------------------+
c     Single Update of Mesh                                            |
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

      integer id(ndf,nkn) !, ix(6,numel)
      real*8 x(ndm,numnp),ang(numnp),dqu(ndx),dp(ndx),sca,alp!,tm(3,3)
 
c---- x        : node coordinates
c---- ang      : angles of nodes      
c---- dqu      : gradient array (dL(X)/dX)
c---- dp       : normalized gradient
c---- sca      : scaling factor
c---- alp      : alpha value
c---- tm       : transformation matrix
      

c==== PROGRAM ==========================================================

c---- Scale dqu --------------------------------------------------------

c      alp=0.4!1.d0  !BAUSTELLE -> alpha-Routine...
      sca=0.d0
      do i1=1,ndx
        if (dabs(dqu(i1)).gt.dabs(sca)) sca=dqu(i1)
      enddo
      if (sca.ge.0.d0.and.sca.le.0.d0) sca=1.d0
      do i1=1,ndx
        dp(i1)=dqu(i1)/sca
      enddo  

c---- Update Mesh ------------------------------------------------------

      i3=0
      do i1=1,nkn
        if (id(1,i1).gt.0) then
c          tm(1,1)= dcos((ang(i1))*(dacos(0.d0)/90.d0))
c          tm(1,2)= dsin((ang(i1))*(dacos(0.d0)/90.d0))
c          tm(1,3)= 0.d0
c          tm(2,1)=-dsin((ang(i1))*(dacos(0.d0)/90.d0))
c          tm(2,2)= dcos((ang(i1))*(dacos(0.d0)/90.d0))
c          tm(2,3)= 0.d0
c          tm(3,1)= 0.d0
c          tm(3,2)= 0.d0
c          tm(3,3)= 1.d0
          
          i3=i3+1
          x(1,i1)=x(1,i1)+(alp*dp(i3)*dcos(ang(i1)*(dacos(0.d0)/90.d0)))         
          x(2,i1)=x(2,i1)+(alp*dp(i3)*dsin(ang(i1)*(dacos(0.d0)/90.d0)))
        endif
        if (id(2,i1).gt.0) then
          i3=i3+1
          x(2,i1)=x(2,i1)+alp*dp(i3)
        endif
      enddo     
      


cjw      i3=0
cjw      do j1=1,nkn
cjw        do j2=1,2
cjw          if (id(j2,j1).gt.0) then
cjwc            idx(1)=j1
cjwc            idx(2)=j2 
cjw            do i1=1,ndm
cjw              do i2=1,nkn
cjw                if (i1.eq.j2.and.i2.eq.j1) then
cjw                  i3=i3+1 
cjw                  x(i1,i2)=x(i1,i2)+alp*dp(i3)   
cjw                endif
cjw              enddo
cjw            enddo
cjw          endif
cjw        enddo          
cjw      enddo      

c==== End of Subroutine ================================================
      
c999   continue
      return
      end      


c=======================================================================      
c=======================================================================      

      subroutine yltalp(d,u,x,ix,lm,id,f,s,p,ang,atb,etb,ftb,ptb,ctb,
     &                  mtb,vtb,btb,th,dqu,lth,alp)

c----------------------------------------------------------------------+
c     Line-Search Algorithm                                            |
c     Form new Mesh with Gradient Results                              |
c----------------------------------------------------------------------+

c==== DECLARATION ======================================================

c---- Formal Parameters ------------------------------------------------
      USE cdat1
      USE cdata
      USE eldata
      USE iofile
      USE mdat2
      USE sdata
      USE tdata
      USE yltdata1
c      USE yltdata2
      USE yltdata3
      USE yltdata4

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

c---- yltdata3 : tolerances   
c     eps1     : 1st tolerance check
c     eps2     : 2nd tolerance check
c     eps3     : degereration criterion 
c     eps4     : alpha bracketing tolerance 
c     eps5     : step length for direct search 
c     eps6     : zero tolerance for sqp
c     halt     : iterations stopped

c---- yltdata4 : parameter data and history
c     alph     : alpha history

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

c---- dt      : search step length

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
c---- f          : load vector
c---- s(nst,nst) : element matrix
c---- p(ndf,nen) : element vector  

      real*8 ang(numnp),atb(1+nth,2*nkn+2*nth),etb(nth,nkn),ftb(nkn)
      real*8 ptb(1+nth),ctb(2*nkn+2*nth),mtb(2*nth)
      real*8 vtb(2*nkn+2*nth),btb(2*nth+2*nkn,2*nth+2*nkn)
      real*8 lth(3,nth)
      real*8 dqu(ndx),alp

c---- ang      : angle of nodes
c---- atb      : table of constrains (th+th-w+w-,1+th) A
c---- etb      : assembled element matrix E
c---- ftb      : assembled element vector F 
c---- ptb      : right hand side
c---- ctb      : coefficients of objective function c
c---- mtb      : assembled element moment vector m
c---- vtb      : optimization variables (th+th-w+w-)
c---- btb      :
c---- lth      : length of yield-lines (x,y,sqrt(x^2+y^2))

      integer ital
      real*8 sca,a(3),z(3),xtemp(ndm,numnp),quot(ndx),maxalp
      real*8 al,an,au,am,zl,zn,zu,zm
c      real*8 ates(2,40)
      
c---- ital     : counter for alpha iterations     
c---- sca      : scale factor 
c---- quot     : quotient for evaluation max alpha
c---- maxalp   : maximal admissible alpha
c---- a        : alpha values
c---- z        : function values for alpha   
c---- al,an,au : lower,normal,upper alpha
c---- zl,zn,zu : lower,normal,upper function value         
c---- am,zm    : new middle values
c---- ates     : test alphas, values

      logical zero
      
c---- zero     : z becomes zero a runs out of area
c---- halt     : iterations stopped

c==== PROGRAM ==========================================================

      zero=.false.

c---- Compute Maximum alpha (see Bazaraa,S.404; Wolfe) -----------------

      quot=0.d0
      i0=0
      do i1=1,nkn  
        do i2=1,2 !ndf:(x,y)
          if (id(i2,i1).gt.0) then
            i0=i0+1
            if (dqu(i0).gt.0.d0.or.dqu(i0).lt.0.d0) 
     &        quot(i0)=dabs(x(i2,i1)/dqu(i0))
          endif   
        enddo
      enddo
        
      maxalp=0.d0
      do i1=1,ndx
        if (quot(i1).gt.maxalp) maxalp=quot(i1)
      enddo
      do i1=1,ndx
        if (quot(i1).lt.maxalp) maxalp=quot(i1)
      enddo
                 
c---- Scale Gradient ---------------------------------------------------

   10 continue
      sca=dt   !0.d0
      do i1=1,ndx
        if (dabs(dqu(i1)).gt.dabs(sca)) sca=dqu(i1)
      enddo
      
      if (zero) sca=-sca
      
c      do i1=1,ndx
c        dp(i1)=dqu(i1)/sca
c      enddo

cc---- Teststeps --------------------------------------------------------
c
c      
c      do i0=1,40  
c        
c        ates(1,i0)=(i0-20.d0)/10.d0
c        
c        
c        call pzero(xtemp,ndm*numnp)
c        
c        do i1=1,ndm
c          do i2=1,numnp
c            xtemp(i1,i2)=x(i1,i2)
c          enddo
c        enddo
c        
c        call yltmsh(xtemp,ang,id,dp,ates(1,i0))      
c        
c        call yltass(d,u,xtemp,ix,lm,id,f,s,p,atb,etb,ftb,
c     &             ptb,ctb,mtb,lth,30)
c        
c        call yltsqp(atb,btb,ctb,ptb,vtb,vtb,0)
c        
c        call yltevl(vtb,mtb,u,th,ates(2,i0),2)
c        print*, '-ates1-',ates(1,i0),ates(2,i0)
c      enddo
      
c---- Initialize alpha Factor ------------------------------------------

      ital=1
      if (iter.eq.1) then
        alpha=0.d0
c        sca=maxalp
c        sca=1.d0
      endif
c      sca=maxalp
      
      do i1=1,3
        a(i1)=0.d0
        z(i1)=0.d0
      enddo  
    
      if (sca.gt.alpha) then   !Momentan sca als oberer Wert; Alternative: 1.d0
        a(1)=sca
        a(2)=alpha
      elseif (sca.lt.alpha) then
        a(1)=alpha
        a(2)=sca
      else
        a(1)=0.5d0*alpha
        a(2)=1.5d0*alpha  
      endif                       
      
c---- test

c      print*, 'alpha bracketing interval:', a(1),a(2)

c---- Fill xtemp-Field -------------------------------------------------

c  100 continue
      do i3=1,2
        call pzero(xtemp,ndm*numnp)
        
        do i1=1,ndm
          do i2=1,numnp
            xtemp(i1,i2)=x(i1,i2)
          enddo
        enddo
        
        call yltmsh(xtemp,ang,id,dqu,a(i3))      
        
        call yltass(d,u,xtemp,ix,lm,id,f,s,p,atb,etb,ftb,
     &             ptb,ctb,mtb,lth,30)
        
        call yltsqp(atb,btb,ctb,ptb,vtb,vtb,0)
        
        call yltevl(vtb,mtb,u,th,z(i3),2)
        
        if (z(i3).ge.0.d0.and.z(i3).le.0.d0) then
          if (zero) then 
            goto 700
          else
            zero=.true.
            goto 10
          endif
        endif      
      enddo
      
c---- Calculate third Value --------------------------------------------

      if (z(1).lt.z(2)) then
        a(3)=2.d0*a(1)-a(2)

  110   continue
        call pzero(xtemp,ndm*numnp)
        
        do i1=1,ndm
          do i2=1,numnp
            xtemp(i1,i2)=x(i1,i2)
          enddo
        enddo

        call yltmsh(xtemp,ang,id,dqu,a(3))      
        
        call yltass(d,u,xtemp,ix,lm,id,f,s,p,atb,etb,ftb,
     &             ptb,ctb,mtb,lth,30)
        
        call yltsqp(atb,btb,ctb,ptb,vtb,vtb,0)
        
        call yltevl(vtb,mtb,u,th,z(3),2)
        
        if (z(i3).ge.0.d0.and.z(i3).le.0.d0) then
          if (zero) then 
            goto 700
          else
            zero=.true.
            goto 10
          endif
        endif      

        if (z(3).lt.z(1)) then 
          a(3)=2.d0*a(3)-a(1)          
          goto 110
        endif  
        
c        al=a(3)
c        an=a(1)
c        au=a(2)
c        zl=z(3)
c        zn=z(1)
c        zu=z(2)
        
      elseif (z(1).gt.z(2)) then
        a(3)=2.d0*a(2)-a(1)

  120   continue
        call pzero(xtemp,ndm*numnp)
        
        do i1=1,ndm
          do i2=1,numnp
            xtemp(i1,i2)=x(i1,i2)
          enddo
        enddo
  
        call yltmsh(xtemp,ang,id,dqu,a(3))      
        
        call yltass(d,u,xtemp,ix,lm,id,f,s,p,atb,etb,ftb,
     &             ptb,ctb,mtb,lth,30)
        
        call yltsqp(atb,btb,ctb,ptb,vtb,vtb,0)
        
        call yltevl(vtb,mtb,u,th,z(3),2)
        
        if (z(i3).ge.0.d0.and.z(i3).le.0.d0) then
          if (zero) then 
            goto 700
          else
            zero=.true.
            goto 10
          endif
        endif      

        if (z(3).lt.z(2)) then !s.o.
          a(3)=2.d0*a(3)-a(2)
          goto 120
        endif  
        
c        al=a(1)
c        an=a(2)
c        au=a(3)
c        zl=z(1)
c        zn=z(2)
c        zu=z(3)
                
      else
        call pzero(xtemp,ndm*numnp)
        
        do i1=1,ndm
          do i2=1,numnp
            xtemp(i1,i2)=x(i1,i2)
          enddo
        enddo
      
        a(3)=(a(1)+a(2))/2.d0
      
        call yltmsh(xtemp,ang,id,dqu,a(3))      
        
        call yltass(d,u,xtemp,ix,lm,id,f,s,p,atb,etb,ftb,
     &             ptb,ctb,mtb,lth,30)
        
        call yltsqp(atb,btb,ctb,ptb,vtb,vtb,0)
        
        call yltevl(vtb,mtb,u,th,z(3),2)
        
        if (z(i3).ge.0.d0.and.z(i3).le.0.d0) then
          if (zero) then 
            goto 700
          else
            zero=.true.
            goto 10
          endif
        endif      

c        al=a(1)
c        an=a(3)
c        au=a(2)
c        zl=z(1)
c        zn=z(3)
c        zu=z(2)

      endif     
 
c---- Define Bracketing Interval ---------------------------------------

  200 continue
      au=max(a(1),a(2),a(3)) 
      al=min(a(1),a(2),a(3))
      
      do i1=1,3
        if ((a(i1).gt.al.or.a(i1).lt.al).and.
     &      (a(i1).gt.au.or.a(i1).lt.au)) an=a(i1)
      enddo
      
      do i1=1,3
        if (a(i1).ge.al.and.a(i1).le.al) zl=z(i1)
        if (a(i1).ge.au.and.a(i1).le.au) zu=z(i1)
        if (a(i1).ge.an.and.a(i1).le.an) zn=z(i1)
      enddo    
     
c---- Calculate Mid Value am -------------------------------------------

c      alp=(al+au)/2.d0
c      print*, 'alpha...=',alp

      if ((zl*(an-au)+zn*(au-al)+zu*(al-an)).gt.0.d0.or.
     &    (zl*(an-au)+zn*(au-al)+zu*(al-an)).lt.0.d0) then
        am=(al+au)/2.d0-0.5*(((zl-zu)*(au-an)*(an-al))/
     &     (zl*(an-au)+zn*(au-al)+zu*(al-an))) 
      else
        am=(al+au)/2.d0
      endif  

c---- Calculate Function Value of am ----------------------------------- 

      call pzero(xtemp,ndm*numnp)
      
      do i1=1,ndm
        do i2=1,numnp
          xtemp(i1,i2)=x(i1,i2)
        enddo
      enddo
            
      call yltmsh(xtemp,ang,id,dqu,am)      
      
      call yltass(d,u,xtemp,ix,lm,id,f,s,p,atb,etb,ftb,
     &           ptb,ctb,mtb,lth,30)
      
      call yltsqp(atb,btb,ctb,ptb,vtb,vtb,0)
      
      call yltevl(vtb,mtb,u,th,zm,2)
        
        if (z(i3).ge.0.d0.and.z(i3).le.0.d0) then
          if (zero) then 
            goto 700
          else
            zero=.not.zero
            goto 10
          endif
        endif      

c---- Stop Criterion or Redefinition of Bracketing Interval ------------

      if (dabs(am-an).lt.eps4) then
c         alp=(am+an)/2.d0
        if (zm.lt.zn) then
          alp=am 
c          alpha=am         
        else              
          alp=an  
c          alpha=an        
        endif             
      else
        if (zn.lt.zm) then
          if (an.lt.am) then
            a(1)=al
            a(2)=an
            a(3)=am
            z(1)=zl
            z(2)=zn
            z(3)=zm
          else
            a(1)=am
            a(2)=an
            a(3)=au
            z(1)=zm
            z(2)=zn
            z(3)=zu
          endif
        else
          if (an.lt.am) then
            a(1)=an
            a(2)=am
            a(3)=au
            z(1)=zn
            z(2)=zm
            z(3)=zu
          else
            a(1)=al
            a(2)=am
            a(3)=an
            z(1)=zl
            z(2)=zm
            z(3)=zn
          endif
        endif         
        ital=ital+1
        goto 200
      endif    
          
cc---- Decision for alpha -----------------------------------------------
c
c      if (zm.lt.zn) then
c        alp=am
c      else
c        alp=an
c      endif    
c     

      if (dabs(alp).gt.maxalp.and.maxalp.gt.0.d0) alp=maxalp

      if (alp.gt.0.d0.or.alp.lt.0.d0) alpha=alp !max(alp,alpha)
 
c      print*, 'alpha value=', alp, alpha, ' iterations=',ital
      goto 999
      
c---- Error ------------------------------------------------------------

  700 continue
c      print*, '*  ERROR in yltalp - program stopped'
c      stop      
c      print*, '*  Error in yltalp - null step'

c      print*, 'maxalp=',maxalp      
c      sca=0.d0
c      do i1=1,3
c        if (dabs(a(i1)).gt.dabs(sca)) sca=
 
c      alp=0.d0
     
c==== End of Subroutine ================================================
      
999   continue
      return
      end      
      
      
c=======================================================================      
c=======================================================================      

      subroutine yltbet(dqu,dquh,bet)

c----------------------------------------------------------------------+
c     Calculate beta-Scalar                                            |
c----------------------------------------------------------------------+ 

c==== DECLARATION ======================================================

c---- Formal Parameters ------------------------------------------------

      USE yltdata1
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

      real*8 dqu(ndx),dquh(ndx),bet

c---- dqu      : gradient g
c---- dquh     : history gradient g

      real*8 norm1,norm2
      
c---- norm1    : normalized gradient g
c---- norm2    : normalized history gradient      


c==== PROGRAM ==========================================================

c---- Calculate bet (Geiger-Kanzow's Version) --------------------------

      norm1=0.d0
      norm2=0.d0
      
      do i1=1,ndx
        norm1=norm1+dqu(i1)**2
        norm2=norm2+dquh(i1)**2
      enddo
      
      bet=norm1/norm2  
      
c---- test

c      print*, 'beta value=', bet      

c==== End of Subroutine ================================================
      
c999   continue
      return
      end
      

c=======================================================================      
c=======================================================================      

      subroutine yltdeg(dtb,dqu,dquh,dp,ideg,ndeg,isgn)

c----------------------------------------------------------------------+
c     Calculate Degeneration-Case                                      |
c----------------------------------------------------------------------+

c==== DECLARATION ======================================================

c---- Formal Parameters ------------------------------------------------
      USE cdat1
      USE cdata
      USE eldata
      USE iofile
      USE sdata
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

c      integer id(ndf,nkn),ix(6,numel),lm(2,nth)!,lmy(2,ny1)
      integer ideg(ny1),ndeg,isgn
      
c---- ideg     : degeneration vector (entry=1=degenerated)
c---- ndeg     : number of degenerated thetas  
c---- isgn     : switch for sign of Dg      

c      real*8 x(ndm,numnp)
c      real*8 f(ndf,numnp),s(nst,nst),p(nst),ang(numnp) 
c      real*8 atb(1+nth,2*nkn+2*nth),etb(nth,nkn),ftb(nkn)
c      real*8 ptb(1+nth),ctb(2*nkn+2*nth),mtb(2*nth)
c      real*8 vtb(2*nkn+2*nth)!,dtb(ny1,ndx)!,vtbh(2*nkn+2*nth)
c      real*8 qu

      real*8 dtb(ny1,ndx),dqu(ndx),dquh(ndx),dp(ndx)
      real*8 d0(ndeg,ndx),d1(ndeg,ndeg),d2(ndeg,ndx-ndeg),d1i(ndeg,ndeg)      !,d1(ndx,ndx),d2(ndx,ny1-ndx),d1i(ndx,ndx)
      real*8 d0tmp(ndeg,ndx),t(ndx,ndx-ndeg),tt(ndx,ndx)
      real*8 ck(ndeg),ck2(ndeg) !,xh(ndm,numnp),quh
 
c---- dtb      : complete D-matrix of degeneration
c---- dqu      : gradient array (dL(X)/dX)
c---- dquh     : historic gradient
c---- dp       : normalized gradient
c---- d1       : first partition D1 of D
c---- d2       : second partition D2 of D
c---- dtb1i    : inverse of dtb1
c---- t        : projection matrix
c---- tt       : t*t^T
c---- ck       : vector D*dqu for check

      integer isort(ndx),irg
      
c---- isort    : new numbering of D columns in srt ygauss
c---- irg      : rank of matrix      

c==== PROGRAM ==========================================================

c---- Initialization ---------------------------------------------------

c      xh=x
c      quh=qu

c---- Initial Check ----------------------------------------------------

      if (ndeg.gt.ndx) print*, 'ERROR IN YLTDEG: ndeg > ndx!'

c---- Take D-Matrix from dtb -------------------------------------------

      i3=0
      do i1=1,ny1
        if (ideg(i1).eq.1) then
          i3=i3+1
          do i2=1,ndx
            d0(i3,i2)=-dtb(i1,i2)
            d0tmp(i3,i2)=-dtb(i1,i2)
          enddo
        endif  
      enddo

c      call plotmatrix(d0,ndeg,ndx,'d0')
c      call plotmatrix(dqu,1,ndx,'dqu')
c      call plotmatrix(dqu2,1,ndx,'dqu2')

c---- Correction of Small Values ---------------------------------------

      do i1=1,ndeg
        do i2=1,ndx
          if (dabs(d0(i1,i2)).lt.eps6) d0(i1,i2)=0.d0
          if (dabs(d0tmp(i1,i2)).lt.eps6) d0tmp(i1,i2)=0.d0
        enddo
      enddo    

c---- Check For Zero Vector --------------------------------------------

      d0tmp=d0
      do i1=1,ndeg
        do i2=1,ndx
          if (dabs(d0tmp(i1,i2)).lt.eps6) d0tmp(i1,i2)=0.d0
        enddo
        a0=0.d0
        do i2=1,ndx
          a0=a0+dabs(d0tmp(i1,i2))
        enddo
        if (a0.lt.eps6) then
          isgn=1
          goto 999
        endif     
      enddo

c---- Output -----------------------------------------------------------

      do i1=1,ndeg
                     write(iow,2010) i1
        if(ior.lt.0) write(*  ,2010) i1      
        do i2=1,ndx
                       write(iow,2011) d0(i1,i2)  !BAUSTELLE: nopr-Option einbauen
          if(ior.lt.0) write(*  ,2011) d0(i1,i2)
        enddo 
      enddo        
      
c---- Divide into D1 and D2 -------------------------------------------- 
      
c      call plotmatrix(d0tmp,ndeg,ndx,'d0tmp vorher')
      call ygauss(d0tmp,ndx,ndeg,isort,irg)
c      call plotmatrix(d0tmp,ndeg,ndx,'d0tmp nachher')
c      print*, 'irg=', irg
c      print*, 'isort', isort

c---- Check for Degeneration Area Criterion D*g < 0 --------------------

      do i1=1,ndeg
        ck(i1)=0.d0
        ck2(i1)=0.d0
        do i2=1,ndx
          ck(i1)=ck(i1)+d0(i1,i2)*dqu(i2)                         
          ck2(i1)=ck2(i1)+d0(i1,i2)*dquh(i2)
        enddo
      enddo  
      
cc---- Kontrollausgabe
c      
c      print*,'Degenerierung mit ndeg=',ndeg
c         do i1=1,ndeg
c          print*,'Dg =',ck(i1)
c          print*,'Dg2 =',ck2(i1)
c         enddo

c---- Select Different Cases of Degeneration --------------------------- 
      
c      isgn=0                 
      do i1=1,ndeg 
        if (ck(i1).lt.0.d0) then ! else goto adjacent region
c          isgn=isgn+1*10**(i1-1)
          isgn=1
          do i2=1,ndx
            dp(i2)=-dqu(i2)
          enddo  
        elseif (ck2(i1).lt.0.d0) then ! else follow border
c          isgn=isgn+2**10**(i1-1)
          isgn=3
          do i2=1,ndx
            dp(i2)=-dquh(i2)
          enddo  
        else  
          isgn=3
c          isgn=isgn+3*10**(i1-1)
        endif
      enddo
            
c---- Output and Return ------------------------------------------------      
      
      if (isgn.eq.1) then
                     write(iow,2110) ck
        if(ior.lt.0) write(*  ,2110) ck      
        goto 999
      elseif (isgn.eq.2) then
                     write(iow,2120) ck2
        if(ior.lt.0) write(*  ,2120) ck2     
        goto 999
      else
                     write(iow,2130) 
        if(ior.lt.0) write(*  ,2130)      
      endif


c---- If Steepest Descend Step Continue --------------------------------      
      
c  100 continue
      do i1=1,ndeg   
        sca=0.d0
        do i2=1,ndx
          if (dabs(d0(i1,isort(i2))).gt.sca) sca=d0(i1,isort(i2))
        enddo  
        
        if (sca.lt.eps6) then
          isgn=1
          goto 999
        endif  
        
        do i2=1,ndx
          if (i2.le.ndeg) then
            d1(i1,i2)=d0(i1,isort(i2))/sca           
          else
            d2(i1,i2-ndeg)=d0(i1,isort(i2))/sca
          endif
        enddo
      enddo        

c      call plotmatrix(d1,ndeg,ndeg,'d1')
c      call plotmatrix(d2,ndeg,ndx-ndeg,'d2')

c---- Check Singularity ------------------------------------------------

      d1i=d1
      call gauss(ndeg,ndeg,d1i,ck,det,ndet,eps3,ierr)
      
c      if (dabs(det).lt.5.d-18)
      
c      print*, 'deg',det*(10.**ndet)
      
c---- Compute Inverse of D1 --------------------------------------------

      d1i=d1      
      call invert(d1i,ndeg,ndeg)
c      call pivot(d1,ndeg,ndeg,d1i)
c      call plotmatrix(d1i,ndeg,ndeg,'d1i')
      
c---- Compute T-Matrix -------------------------------------------------

      call pzero(t,ndx*(ndx-ndeg))
      do i1=1,ndeg
        do i2=1,ndx-ndeg
          do i3=1,ndeg
            t(i1,i2)=t(i1,i2)-d1i(i1,i3)*d2(i3,i2) 
          enddo
        enddo
      enddo
      do i1=1,ndx-ndeg
        t(i1+ndeg,i1)=1.d0
      enddo
 
cc---- Normalize T ------------------------------------------------------ 
c      
c      do i1=1,ndeg
c        sca=0.d0
c        do i2=1,ndx-ndeg
c          if (dabs(t(i1,i2)).gt.dabs(sca)) sca=t(i1,i2)
c        enddo
c        do i2=1,ndx-ndeg
c          t(i1,i2)=t(i1,i2)/sca
cc          if (t(i1,i2).lt.eps6) t(i1,i2)=0.d0
c        enddo
c      enddo    
c      
c      call plotmatrix(t,ndx,ndx-ndeg,'t')
      
c---- Compute dp for Case Dg > 0 and Dg' > 0 ---------------------------

      call pzero(tt,ndx*ndx)
      call pzero(dp,ndx)
      do i1=1,ndx
        do i2=1,ndx
          do i3=1,ndx-ndeg
            tt(i1,i2)=tt(i1,i2)+t(i1,i3)*t(i2,i3)
          enddo
        enddo
      enddo

c      call plotmatrix(tt,ndx,ndx,'tt')

cc---- Normalize TT -----------------------------------------------------
c
c      sca=0.d0
c      do i1=1,ndeg
c        do i2=1,ndx
c          if (dabs(tt(i1,i2)).gt.dabs(sca)) sca=tt(i1,i2)
c        enddo
c      enddo    
c      do i1=1,ndeg
c        do i2=1,ndx
c          tt(i1,i2)=tt(i1,i2)/sca
c          if (i1.ne.i2) tt(i2,i1)=tt(i2,i1)/sca
c          if (dabs(tt(i1,i2)).lt.eps6) tt(i1,i2)=0.d0
c        enddo
c      enddo    
c        
c      
c      call plotmatrix(tt,ndx,ndx,'tt-norm')
      
      dp=0.d0
      do i1=1,ndx  
        do i2=1,ndx
          dp(isort(i1))=dp(isort(i1))+tt(i1,i2)*dqu(isort(i2))  !idee?       
c          dp(i1)=dp(i1)+tt(i1,i2)*dquh(i2)          
        enddo
      enddo      
            

cc---- Kontrollausgabe
c
c      do i1=1,ndx
c        print*, 'dp in yltdeg nach tt', dp(i1)
c      enddo  

c---- Scale dp ---------------------------------------------------------

      sca=0.d0
      do i1=1,ndx
        if (dabs(dp(i1)).gt.dabs(sca)) sca=dp(i1)
      enddo
      do i1=1,ndx
        dp(i1)=dp(i1)/sca
      enddo  

ccc---- Kontrollausgabe
c
c      do i1=1,ndx
c        print*, 'norm-dp in yltdeg nach tt', dp(i1)
c      enddo  


cctest>>>
c
c      dp=0.d0
c      do i1=1,ndeg
c        do i2=1,ndx
c          dp(i2)=dp(i2)+dqu(i2)*d0(i1,i2)
c        enddo  
c      enddo
c      sca=0.d0
c      do i1=1,ndx
c        if (dabs(dp(i1)).gt.dabs(sca)) sca=dp(i1)
c      enddo
c      do i1=1,ndx
c        dp(i1)=dp(i1)/sca
c      enddo  
c
cctest<<<




c==== Formats ==========================================================

2010  format(/'  Degeneration vector row ',1i3)
2011  format(3x,1pe15.7)
2110  format(/'  Degeneration check  Dg = ',1pe15.7)
2120  format(/'  Degeneration check  Dg`= ',1pe15.7)
2130  format(/'  Steepest Descend Step   ')
            
c==== End of Subroutine ================================================
      
999   continue
      return
      end
            
    
cc=======================================================================      
cc=======================================================================      
c
c      subroutine ylttry2(d,u,x,ix,lm,id,f,s,p,ang,atb,etb,ftb,ptb,ctb,
c     &                  mtb,vtb,btb,th,dqu,lth,qu)
c
cc----------------------------------------------------------------------+
cc     Direct Search                                                    |
cc----------------------------------------------------------------------+
c
cc==== DECLARATION ======================================================
c
cc---- Formal Parameters ------------------------------------------------
c
c      USE cdat1
c      USE cdata
c      USE eldata
c      USE iofile
c      USE mdat2   
c      USE sdata   
c      USE tdata
c      USE yltdata1
c      USE yltdata2
c      USE yltdata3
c      USE yltdata4
c      implicit real*8 (a-h,o-z)
c      implicit integer (i,n)
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
cc---- yltdata3 : tolerances   
cc     eps1     : 1st tolerance check
cc     eps2     : 2nd tolerance check
cc     eps3     : degereration criterion 
cc     eps4     : alpha bracketing tolerance 
cc     eps5     : step length for direct search 
cc     eps6     : zero tolerance for sqp
cc     halt     : iterations stopped
c
c
cc---- yltdata4 : parameter data and history
cc     alph     : alpha history
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
cc---- dt      : search step length
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
cc---- mdat2      : rotation data
cc     ia         :
cc     itrot      :  
c
c
cc---- sdata      :
cc     ndf        :
cc     ndm        :
cc     nen1       :
cc     nst        :
c
c
cc      real*8  d(ndd,*),u(*),x(*) 
cc      integer id(*), ix(*)
c      integer id(ndf,nkn), ix(6,numel),lm(2,nth)
c      real*8 d(ndd,*),u(ndf,*),x(ndm,numnp)
c      real*8 f(ndf,numnp),s(nst,nst),p(nst) 
c
cc      integer ld(*), ie(nie,*), id(ndf,*), ix(nen1,*), jp(*)
cc      real*8  xl(ndm,*), p(*), s(nst,*), d(ndd,*), x(ndm,*), f(ndf,*)
cc      real*8  b(*), a(*), c(*), ul(nst,*), tl(*), t(*), u(ndf,*)
cc      real*8  f0(ndf,*), ud(*) 
c
cc---- ix(nen)    : element global node numbers
cc---- id(*)      : boundary conditions
cc---- d(*)       : element data parameters
cc---- u(*)       : 
cc---- x(*)       : node coordinates
cc---- t(nen)     : 
cc---- f          : load vector
cc---- s(nst,nst) : element matrix
cc---- p(ndf,nen) : element vector  
c
c      real*8 ang(numnp),atb(1+nth,2*nkn+2*nth),etb(nth,nkn),ftb(nkn)
c      real*8 ptb(1+nth),ctb(2*nkn+2*nth),mtb(2*nth)
c      real*8 vtb(2*nkn+2*nth),btb(2*nth+2*nkn,2*nth+2*nkn)
c      real*8 lth(3,nth)
c      real*8 dqu(ndx),dp(ndx),alp
c
cc---- ang      : angle of nodes
cc---- atb      : table of constrains (th+th-w+w-,1+th) A
cc---- etb      : assembled element matrix E
cc---- ftb      : assembled element vector F 
cc---- ptb      : right hand side
cc---- ctb      : coefficients of objective function c
cc---- mtb      : assembled element moment vector m
cc---- vtb      : optimization variables (th+th-w+w-)
cc---- btb      :
cc---- lth      : length of yield-lines (x,y,sqrt(x^2+y^2))
c
c
cc---- Variables --------------------------------------------------------
c
c      integer ndir
c      
cc---- ndir     : search dof's (directions, 2-dimensional)      
c
c      real*8 xtemp(2,ndm,numnp),xh(ndm,numnp),xstep(2,numnp),quh
c      real*8 v1(3),v2(3),v3(3),av,z(2),tol, dec
c      real*8 tax(2,numnp)
c
cc---- xtemp*   : temporary coordinates (trial)
cc---- xh       : x history
cc---- quh      : qu history
cc---- xstep    : admitted dof steps
cc---- v*       : vector
cc---- av       : vector norm
cc---- z        : load factor
cc---- tol      : tolerance
cc---- dec      : decreasing factor of step length
cc---- tax      : taxation of ultimate load per test
c
c      logical admit(2),out
c      
cc---- admit    : admitted trial mesh  
cc---- out      : stop criterion
c
cc==== PROGRAM ==========================================================
c
cc---- Initialization ---------------------------------------------------
c
c      call pzero (xstep,2*numnp)
c      xh=x
c      quh=qu
c      tol=1.d-8
c      dec=0.1d0
c      ndir=2
c 
c                   write(iow,800) dt  !BAUSTELLE: nopr-Option einbauen
c      if(ior.lt.0) write(*  ,800) dt
c 
ccc---- Kontrollausgabe
cc
cc      do i1=1,numnp
cc        print*, 'ang', ang(i1)
cc      enddo   
c                   
cc==== Define Admitted Search Directions (dofs) =========================
c
cc      istep=0
c      do i1=1,numnp
c        do i2=1,ndir
c          if (id(i2,i1).lt.0) then
c            xstep(i2,i1)=0.d0
c          else
c            xstep(i2,i1)=dt !1.d0
cc            istep=istep+1
c          endif
c        enddo
c      enddo    
c
ccc---- Kontrollausgabe
cc
cc      print*
cc      do i1=1,numnp
cc        print*, xstep(1,i1), xstep(2,i1)
cc      enddo
c      
cc==== Trial Loop =======================================================
c
c      call pzero(xtemp,2*ndm*numnp)
c      
c  100 continue   
c      do i1=1,numnp
c  
c  110   continue
c  
c        do i2=1,ndir
c          if (xstep(i2,i1).eq.0.d0) cycle
c
c  120     continue
c          do i3=1,numnp
c            if (ang(i3).ne.0.d0) then   !rotated boundaries
c              xtemp(1,1,i3)=xh(1,i3)+xstep(1,i3)
c     &                      *dcos(dacos(0.d0)*ang(i3)/9.d1)
c              xtemp(2,1,i3)=xh(1,i3)-xstep(1,i3)
c     &                      *dcos(dacos(0.d0)*ang(i3)/9.d1)
c              xtemp(1,2,i3)=xh(2,i3)+xstep(1,i3)
c     &                      *dsin(dacos(0.d0)*ang(i3)/9.d1)
c              xtemp(2,2,i3)=xh(2,i3)-xstep(1,i3)
c     &                      *dsin(dacos(0.d0)*ang(i3)/9.d1)
c            else  
c              do i4=1,ndir                  
c                xtemp(1,i4,i3)=xh(i4,i3)+xstep(i4,i3)
c                xtemp(2,i4,i3)=xh(i4,i3)-xstep(i4,i3)
c              enddo
c            endif            
c          enddo  
c          
cc          do i3=1,numnp
cc            do i4=i,ndir
c              
c
cc-------- Test Direction and check Admittance --------------------------
c
c          do i0=1,2
c            v1=0.d0
c            v2=0.d0
c            v3=0.d0
c            av=0.d0
c            admit=.true.
c            do i3=1,numel
c              v1(i2)=xtemp(i0,i2,ix(2,i3))-xtemp(i0,i2,ix(1,i3))  
c              v2(i2)=xtemp(i0,i2,ix(2,i3))-xtemp(i0,i2,ix(3,i3))  
c              call vcross(v1,v2,v3)
cc              call vnorm(v3,av)
c              if (v3(3).gt.0.d0) admit(i0)=.false.
c            enddo  
c            z(i0)=qu            
c            
c            if (admit(i0)) then 
c              call yltmsh(xtemp(i0,1:ndm,1:numnp),ang,id,dqu,1.d0)    !BAUSTELLE 24.8. alpha checken...            
c              call yltass(d,u,xtemp(i0,1:ndm,1:numnp),ix,lm,id,f,s,p,
c     &                    atb,etb,ftb,ptb,ctb,mtb,lth,30)
c              call yltsqp(atb,btb,ctb,ptb,vtb,vtb,0)
c              call yltevl(vtb,mtb,u,th,z(i0),1)
c            endif
c          enddo  
c        
cc-------- If both Patterns not admitted, reduce Step Length ------------
c
c          if ((.not.(admit(1))).and.(.not.(admit(2)))) then
c            xstep(i2,i1)=xstep(i2,i1)*dec
c            call pzero(xtemp,2*ndm*numnp)
c            if (xstep(i2,i1).lt.tol) then
c              cycle
c            endif  
c            goto 120
c          endif      
c          
cc------ Decision for Basispoint: foreward/backward/stay ----------------
c
c          if (z(1).lt.z(2)) then
c            if (z(1).lt.quh) then
c              quh=z(1)
c              xh(1:ndm,1:numnp)=xtemp(1,1:ndm,1:numnp)
c            else
c              xstep(i2,i1)=xstep(i2,i1)*dec  
c            endif
c          elseif (z(2).lt.z(1)) then
c            if (z(2).lt.quh) then
c              quh=z(2)
c              xh(1:ndm,1:numnp)=xtemp(2,1:ndm,1:numnp)
c            else  
c              xstep(i2,i1)=xstep(i2,i1)*dec  
c            endif
c          else
c            if (z(1).lt.quh) then
c              quh=z(1)
c              xh(1:ndm,1:numnp)=xtemp(1,1:ndm,1:numnp)
c            else
c              xstep(i2,i1)=xstep(i2,i1)*dec  
c            endif              
c          endif 
c          
ccc---- Kontrollausgabe
cc
cc          print*,'xh', i1, i2
cc          do ii=1,numnp
cc            print*,xh(1,ii),xh(2,ii),xh(3,ii)
cc          enddo            
cc          print*,'xtemp1', i1, i2
cc          do ii=1,numnp
cc            print*,xtemp(1,1,ii),xtemp(1,2,ii),xtemp(1,3,ii)
cc          enddo            
cc          print*,'xtemp2', i1, i2
cc          do ii=1,numnp
cc            print*,xtemp(2,1,ii),xtemp(2,2,ii),xtemp(2,3,ii)
cc          enddo            
c
c         print*, 'numnp,xstep', i1,xstep(i2,i1)          
c          
c        enddo                
c      enddo
c      
cc---- Decision for Pattern Point ---------------------------------------
c
c      if (quh.lt.qu.and.(dabs(qu-quh).gt.tol)) then
c        qu=quh
c        x(1:ndm,1:numnp)=xh(1:ndm,1:numnp)
c                     write(iow,810)  !BAUSTELLE: nopr-Option einbauen
c        if(ior.lt.0) write(*  ,810)
c      else
c        out=.true.
c        do i1=1,numnp
c          do i2=1,ndir
c            if (xstep(i2,i1).gt.tol) then
c              out=.false.
c            endif 
c          enddo   
c        enddo
c        dt=dt*0.1d0
c        do i1=1,numnp
c          do i2=1,ndir
c            xstep(i2,i1)=xstep(i2,i1)*dec
cc            print*, 'step', i1,i2,step(i2,i1)
c          enddo
c        enddo
c        call pzero(xtemp,2*ndm*numnp)
c        if (.not.out) then
c          goto 100
c        else
c                       write(iow,820)  !BAUSTELLE: nopr-Option einbauen
c          if(ior.lt.0) write(*  ,820)
c        endif    
c      endif
c
ccc==== Replace Mesh =====================================================
cc
cc  200 continue 
cc      call plotmatrix(x,ndm,numnp,'x')        
c
cc==== Formats ==========================================================
c
c  800 format(/'  Direct search with step length =',e15.7)
c  810 format(/'  New pattern move achieved ')
c  820 format(/'  Tolerance of search step length  met in ylttry ')
c
cc==== End of Subroutine ================================================
c      
c999   continue
c      return
c      end
     
c=======================================================================      
c=======================================================================      

      subroutine ylttry(d,u,x,ix,lm,id,f,s,p,ang,atb,etb,ftb,ptb,ctb,
     &                  mtb,vtb,btb,th,dqu,lth,qu)

c----------------------------------------------------------------------+
c     Direct Search                                                    |
c----------------------------------------------------------------------+

c==== DECLARATION ======================================================

c---- Formal Parameters ------------------------------------------------
      USE cdat1
      USE cdata
      USE eldata
      USE iofile
      USE mdat2
      USE sdata
      USE tdata
      USE yltdata1
c      USE yltdata2
      USE yltdata3
      USE yltdata4

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

c---- yltdata3 : tolerances   
c     eps1     : 1st tolerance check
c     eps2     : 2nd tolerance check
c     eps3     : degereration criterion 
c     eps4     : alpha bracketing tolerance 
c     eps5     : step length for direct search 
c     eps6     : zero tolerance for sqp
c     halt     : iterations stopped

c---- yltdata4 : parameter data and history
c     alph     : alpha history

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
c---- f          : load vector
c---- s(nst,nst) : element matrix
c---- p(ndf,nen) : element vector  

      real*8 ang(numnp),atb(1+nth,2*nkn+2*nth),etb(nth,nkn),ftb(nkn)
      real*8 ptb(1+nth),ctb(2*nkn+2*nth),mtb(2*nth)
      real*8 vtb(2*nkn+2*nth),btb(2*nth+2*nkn,2*nth+2*nkn)
      real*8 lth(3,nth)
      real*8 dqu(ndx) !,dp(ndx),alp

c---- ang      : angle of nodes
c---- atb      : table of constrains (th+th-w+w-,1+th) A
c---- etb      : assembled element matrix E
c---- ftb      : assembled element vector F 
c---- ptb      : right hand side
c---- ctb      : coefficients of objective function c
c---- mtb      : assembled element moment vector m
c---- vtb      : optimization variables (th+th-w+w-)
c---- btb      :
c---- lth      : length of yield-lines (x,y,sqrt(x^2+y^2))


c---- Variables --------------------------------------------------------

      integer iswi(ndx) !,ndir
      
c---- ndir     : search dof's (directions, 2-dimensional) 
c---- iswi     : switch for search-direction ndx     

      real*8 xtemp(ndm,numnp),xh(ndm,numnp),xstep(2,numnp),quh
      real*8 v1(3),v2(3),v3(3) !,av,z(2),tol, dec
c      real*8 tax(2,numnp)

c---- xtemp*   : temporary coordinates (trial)
c---- xh       : x history
c---- quh      : qu history
c---- xstep    : admitted dof steps
c---- v*       : vector
c---- av       : vector norm
c---- z        : load factor
c---- tol      : tolerance
c---- dec      : decreasing factor of step length
c---- tax      : taxation of ultimate load per test

      logical admit !,out
      
c---- admit    : admitted trial mesh  
c---- out      : stop criterion

c==== PROGRAM ==========================================================

c---- Initialization ---------------------------------------------------

      call pzero (xstep,3*numnp)
      xh=x
      quh=qu
c      tol=1.d-8
c      dec=0.1d0
c      ndir=2
 
                   write(iow,800) dt  !BAUSTELLE: nopr-Option einbauen
      if(ior.lt.0) write(*  ,800) dt
 
cc---- Kontrollausgabe
c
c      do i1=1,numnp
c        print*, 'ang', ang(i1)
c      enddo   
                   
c==== Define Admitted Search Directions (dofs) =========================

      do i1=1,numnp
        do i2=1,2
          if (id(i2,i1).lt.0) then
            xstep(i2,i1)=0.d0
          else
            xstep(i2,i1)=dt !1.d0
          endif
        enddo
      enddo    

c==== Actual Function Value ============================================

      call yltass(d,u,xh(1:ndm,1:numnp),ix,lm,id,f,s,p,
     &            atb,etb,ftb,ptb,ctb,mtb,lth,30)
      call yltsqp(atb,btb,ctb,ptb,vtb,vtb,0)
      call yltevl(vtb,mtb,u,th,quh,1)
      qu=quh

c==== Trial Loops ======================================================

c  100 continue
      do i1=1,ndx
        iswi(i1)=-1
      enddo  
      n1=2 !1..ndx (Zehnerstelle)

  110 continue
      if (iswi(1).gt.1) then
        iswi(1)=0
  150   continue
        iswi(n1)=iswi(n1)+1
        if (iswi(ndx).gt.1) goto 200
        if (iswi(n1).gt.1) then
          iswi(n1)=0
          n1=n1+1
          goto 150
        endif  
        n1=2   
      endif

      xtemp=xh
      i0=0
      do i1=1,numnp
        do i2=1,2

          if (id(i2,i1).gt.0) then          
            if (ang(i1).gt.0.d0.or.ang(i1).lt.0.d0) then
              i0=i0+1
              xtemp(1,i1)=xh(1,i1)+iswi(i0)*xstep(1,i1)
     &                    *dcos(dacos(0.d0)*ang(i1)/9.d1)
              xtemp(2,i1)=xh(2,i1)+iswi(i0)*xstep(1,i1)
     &                    *dsin(dacos(0.d0)*ang(i1)/9.d1)
            else  
              i0=i0+1
              xtemp(i2,i1)=xh(i2,i1)+iswi(i0)*xstep(i2,i1)
            endif
          endif
        enddo
      enddo         
      
      do i1=1,2
        v1=0.d0
        v2=0.d0
        v3=0.d0
c        av=0.d0
        admit=.true.
        do i2=1,numel
          v1(i1)=xtemp(i1,ix(2,i2))-xtemp(i1,ix(1,i2))  
          v2(i1)=xtemp(i1,ix(2,i2))-xtemp(i1,ix(3,i2))  
          call vcross(v1,v2,v3)
c          call vnorm(v3,av)
          if (v3(3).gt.0.d0) admit=.false.
        enddo  
         
c         print*, 'admit', admit           
c                   
c             do ii=1,numnp
c             print*,'xtemp',xtemp(1,ii),xtemp(2,ii),xtemp(3,ii)      
c                   enddo
                   
        if (admit) then 
          call yltmsh(xtemp(1:ndm,1:numnp),ang,id,dqu,1.d0)    !BAUSTELLE 24.8. alpha checken...            
          call yltass(d,u,xtemp(1:ndm,1:numnp),ix,lm,id,f,s,p,
     &                atb,etb,ftb,ptb,ctb,mtb,lth,30)
          call yltsqp(atb,btb,ctb,ptb,vtb,vtb,0)
          call yltevl(vtb,mtb,u,th,qutemp,1)
        endif

        if (qutemp.lt.qu) then
          x=xtemp
          qu=qutemp            
        endif    
      enddo  

      iswi(1)=iswi(1)+1      
      goto 110            

c      print*, 'jetzt qu', qu

c---- Decision for Pattern Point ---------------------------------------

  200 continue

      if (qu.ge.quh.and.qu.le.quh) then
                       write(iow,830)  !BAUSTELLE: nopr-Option einbauen
          if(ior.lt.0) write(*  ,830)
      endif


c      if (quh.lt.qu.and.(dabs(qu-quh).gt.tol)) then
c        qu=quh
c        x(1:ndm,1:numnp)=xh(1:ndm,1:numnp)
c                     write(iow,810)  !BAUSTELLE: nopr-Option einbauen
c        if(ior.lt.0) write(*  ,810)
c      else
c        out=.true.
c        do i1=1,numnp
c          do i2=1,ndir
c            if (xstep(i2,i1).gt.tol) then
c              out=.false.
c            endif 
c          enddo   
c        enddo
c        dt=dt*0.1d0
c        do i1=1,numnp
c          do i2=1,ndir
c            xstep(i2,i1)=xstep(i2,i1)*dec
cc            print*, 'step', i1,i2,step(i2,i1)
c          enddo
c        enddo
cc        call pzero(xtemp,2*ndm*numnp)
c        if (.not.out) then
c          goto 100
c        else
c                       write(iow,820)  !BAUSTELLE: nopr-Option einbauen
c          if(ior.lt.0) write(*  ,820)
c        endif    
c      endif      

c==== Formats ==========================================================

  800 format(/'  Direct search with step length =',e15.7)
c  810 format(/'  New pattern move achieved ')
c  820 format(/'  Tolerance of search step length  met in ylttry ')
  830 format(/'  Minimum for step length reached in ylttry ')

c==== End of Subroutine ================================================
      
c999   continue
      return
      end
     
           