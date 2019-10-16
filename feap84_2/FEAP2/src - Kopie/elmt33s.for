      subroutine elmt33(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c**********************************************************************|
c.....linear Reissner Mindlin plate element                            |
c     with one-point integration and stabilization matrix              |
c     due to Belytschko, Tsay: IJNME 19(83)405-419                     |
c     A stabilization procedure for the quadrilateral plate element    |
c     with one-point quadrature                                        |
c.... Material parameters (input)                                      |
c      1: E      = Young's modulus                                     |
c      2: nu     = Poisson's ratio                                     |
c      3: h      = thickness of shell                                  |
c      4: q      = transverse load                                     |
c      5: rho = density                                                |
c      6: r_w    = stiffness value for w                               |
c      7: r_beta =  stiffness value for beta                           |
c                                                                      |   
c                                                                      |   
c      signs of stresses, see lecture CTWM                             |
c      signs of displacements FEM-Like w,phix,phiy                     |
c                                                                      |   
c      w.wagner uka 11/03           
c      w.wagner kit 11/13 update h1-3                                  |
c      w.wagner kit 11/13 QLOA                                         |
c**********************************************************************|
c      open/remarks                                                    |   
c      one-point part standard formulation                             |   
c      may be changed to closed solution of B/T                        |   
c      reac R=k*v-P                                                    |   
c      mass                                                            |   
c                                                                      |   
c**********************************************************************|
      USE bdata      
      USE cdata
      USE eldata
      USE iofile
      USE pdata7
      USE pdata10
      USE prlod
      USE pltran
      USE qload      
      USE strnam 
      implicit double precision (a-h,o-z)
      dimension xl(ndm,*),tl(*),ix(*),ul(ndf,*),s(nst,nst),p(nst),d(*),
     1    gp(2),pgb(4),wgb(4),h1(*),h2(*),h3(*),
     2    cartd(2,4),deriv(2,4),shp(4),gr(6),sig(5),db(3,3),ds(2,2),
     3    bbi(3,3),bbj(3,3),bsi(2,3),bsj(2,3),yl(3,4),
     4    h(4),b1(4),b2(4),x(4),y(4),g(4),gtg(4,4)
      dimension sh(12,12)
      data h /1,-1,1,-1/
      save maold,db,ds,c1q,c2q
c      
c.... go to correct array processor
      go to(1,2,3,4,5,3,2,4,2,2,3,2,2,2,2,2,2,2,2,2,2,22) isw
c.... input material properties
1     if(ior.lt.0) write(*,1001)
      call dinput(d,7)
                   write(iow,1002) (d(i),i=1,7)
      if(ior.lt.0) write(*  ,1002) (d(i),i=1,7)
      if(ndm.eq.3) ipla = 1
c.... description of stresses  
      strsus( 1) = '  MOMENT M_xx  '
      strsus( 2) = '  MOMENT M_xy  '
      strsus( 3) = '  MOMENT M_yy  '
      strsus( 4) = '               '
      strsus( 5) = '  MOMENT M_11  '
      strsus( 6) = '  MOMENT M_22  '
      strsus( 7) = '  ANGLE Phi_1  '
      strsus( 8) = '  SHEAR FORCE_x'
      strsus( 9) = '  SHEAR FORCE_y'
      do is =10,25
        strsus(is) = '               '
      end do
c...  material matrix
      call modplbt(db,ds,d)
      c1q = d(6)*5.d0/6.d0*(d(1)/(2.d0*(1.d0+d(2))))*d(3)**3/12.d0
      c2q = d(7)*d(1)*d(3)**3/(1.d0-d(2)**2)/192.d0
      maold = ma
2     return
c.... stiffness matrix
c.... elasticity matrix (if necessary)
3     if(ma.ne. maold) then
        call modplbt(db,ds,d)
        c1q = d(6)*5.d0/6.d0*(d(1)/(2.d0*(1.d0+d(2))))*d(3)**3/12.d0
        c2q = d(7)*d(1)*d(3)**3/(1.d0-d(2)**2)/192.d0
        maold = ma
      end if
c.... 1 point integrated standard matrix
      xsi = 0.d0
      eta = 0.d0
c.....shape functions,area etc.
      call sfr2plbt(deriv,eta,xsi,nel,shp)
      call jacoplbt(cartd,deriv,djacb,xl,n,nel,ndm)
      da = djacb*4.d0
c.....displacement gradients
      call gradplbt(cartd,gr,ul,ndf,nel)
c.....stresses and strains
      call strplbt(cartd,shp,ul,gr,db,ds,sig,nel,ndf)
      do 31  ino =1,4
        call bmatplbt(bbi,bsi,cartd,ino,shp,nel)
c.....  residual
        call mulplbt(bbi,sig(1),da,p,ndf,nst,ino,3,3,1)
        call mulplbt(bsi,sig(4),da,p,ndf,nst,ino,3,2,1)
c.....  external load vector due to transversal constant load 
        call qload33(qz,d,aqloa,numel,n,mqloa,propq,prop,isw)
        iq = (ino-1)*3+1
        p(iq) = p(iq) + shp(ino)*qz*da
        if(isw.eq.6) goto 31
        do 32 jno = ino,4 !  symmetry noted
          call bmatplbt(bbj,bsj,cartd,jno,shp,nel)
          call subplbt(bbi,bbj,da,db,s,ndf,nst,ino,jno,3,3,3,1,1)
   32     call subplbt(bsi,bsj,da,ds,s,ndf,nst,ino,jno,3,2,3,1,1)
   31 continue
c
c.... stabilization matrix
c.... vector b1
      b1(1)= 0.5d0*(xl(2,2)-xl(2,4)) 
      b1(2)= 0.5d0*(xl(2,3)-xl(2,1)) 
      b1(3)= 0.5d0*(xl(2,4)-xl(2,2)) 
      b1(4)= 0.5d0*(xl(2,1)-xl(2,3)) 
c.... vector b2
      b2(1)= 0.5d0*(xl(1,4)-xl(1,2)) 
      b2(2)= 0.5d0*(xl(1,1)-xl(1,3)) 
      b2(3)= 0.5d0*(xl(1,2)-xl(1,4)) 
      b2(4)= 0.5d0*(xl(1,3)-xl(1,1)) 
c.... vector x
      x(1)= xl(1,1)
      x(2)= xl(1,2)
      x(3)= xl(1,3)
      x(4)= xl(1,4)
c.... vector y
      y(1)= xl(2,1)
      y(2)= xl(2,2)
      y(3)= xl(2,3)
      y(4)= xl(2,4)
c.... area
      ax  = dot(x,b1,4)
c     ay  = dot(y,b2,4)
c.... h^T*x, h^T*y
      htx = dot(h,x,4)      
      hty = dot(h,y,4)      
c.... gamma
      do i=1,4
        g(i) = h(i) -( htx*b1(i) + hty*b2(i) )/ax
      end do
c.... parameter c1,c2
      b1tb2 = dot(b1,b1,4)+dot(b2,b2,4)
      c1 = c1q/ax/ax*b1tb2
      c2 = c2q/ax*b1tb2
c.... gamma^T*gamma
      do i=1,4
        do j=1,4
          gtg(i,j) = g(i)*g(j)
        end do
      end do
c.... store into sh 

      sh=0

      ii=0
      do i=1,4
        jj=0
        do j=1,4 
          sh(ii+1,jj+1) = sh(ii+1,jj+1) + c1*gtg(i,j)       
          sh(ii+2,jj+2) = sh(ii+2,jj+2) + c2*gtg(i,j)       
          sh(ii+3,jj+3) = sh(ii+3,jj+3) + c2*gtg(i,j)       
          jj = jj+3
        end do
        ii = ii+3
      end do

c...  add to S
      s=s+sh

c      
c.....lower part of K and Kh
      do 36 i = 1,nst
        do 36 j = 1,i
          sh(i,j) = sh(j,i)
   36     s (i,j) = s (j,i)

c.... calculate eigenvalues
c     call elemev(s,nel,ndf,nst,2)

c.... form additional residual part R=(-K*V+P)-KH*V   
      do 37 i = 1,nst
      do 37 j = 1,nst
         p(i) = p(i) - sh(i,j)*ul(j,1)
37    continue

      return
c
c.....stress resultants at center of element
4     istv = 9
c.... elasticity matrix (if necessary)
      if(ma.ne. maold) then
        call modplbt(db,ds,d)
        maold = ma
      end if
      xsi = 0.d0
      eta = 0.d0
c.....shape functions,area etc.
      call sfr2plbt(deriv,eta,xsi,nel,shp)
      call jacoplbt(cartd,deriv,djacb,xl,n,nel,ndm)
      da = djacb*4.d0
c.....displacement gradients
      call gradplbt(cartd,gr,ul,ndf,nel)
c.....stresses and strains
      call strplbt(cartd,shp,ul,gr,db,ds,sig,nel,ndf)
      if(isw.eq.4) then     
c.....  Output stresses  (M and Q)
c.....  coordinates
        do 41 idm = 1,2
          gp(idm) = 0.0
          do 41 ino = 1,nel
   41       gp(idm) = gp(idm) + xl(idm,ino) * shp(ino)
        mct = mct - 1
        if(mct.gt.0) go to 42
c.....  head
                     write(iow,2001) o,head
        if(ior.lt.0) write(*  ,2001) o,head
        mct = 50
c.....  print values
   42                write(iow,2002) n,ma,(gp(i),i=1,2),(sig(i),i=1,5)
        if(ior.lt.0) write(*  ,2002) n,ma,(gp(i),i=1,2),(sig(i),i=1,5)
      elseif(isw.eq.8) then     
c....   stresses for plot M,Q
        if(iplma(ma).eq.0)  return ! only if MATN
        call plotplbt(ix,strea,strea(1+numnp),shp,sig,nel,da,numnp)
      elseif(isw.eq.13) then     
c.....  plot stress resultants from center gauss-point without averaging 
        if(flfp) then
c....     calculate extreme values
          xmaxf = max(xmaxf,sig(nfp))
          xminf = min(xminf,sig(nfp))
        else
c....     plot stresses elementwise in one color      
c....     calculate color      
          call pppcolf(sig(nfp))
c....     plot stress      
c....     coordinates
          call pzero(yl,12)
          do i=1,ndm
            do j=1,4
              yl(i,j)=xl(i,j) 
            end do
          end do
c.....    transform for rot        
          call plxtrn(yl,tra,vr,3,4)
c......   plot element
          call plot9 (iel,ix,yl,ndm,nel,1)
        end if
      end if  
      return
c.... mass matrix  
5     ngb = 2
      dh = d(3)*d(3)/12.d0
c.....gauss points
      call gausplbt(ngb,pgb,wgb)
c.... loop gauss points 
      do 50 igb = 1,ngb
        do 50 jgb = 1,ngb
          xsi = pgb(igb)
          eta = pgb(jgb)
c.....    shape functions,area etc.
          call sfr2plbt(deriv,eta,xsi,nel,shp)
          call jacoplbt(cartd,deriv,djacb,xl,n,nel,ndm)
          da = djacb*wgb(igb)*wgb(jgb)*d(3)*d(5)
          do 51  ino =1,4
            ii=(ino-1)*ndf+1
            p(ii  ) = p(ii) + shp(ino)*da
            p(ii+1) = p(ii)*dh  
            p(ii+2) = p(ii)*dh  
            do 52  jno = ino,4
              jj=(jno-1)*ndf+1
              sm = shp(ino)*shp(jno)*da
              s(ii  ,jj  ) = s(ii  ,jj  ) + sm
              s(ii+1,jj+1) = s(ii+1,jj+1) + sm*dh
              s(ii+2,jj+2) = s(ii+2,jj+2) + sm*dh
   52       end do
   51   end do
   50 end do
c.....upper part of M
      do 56 i = 1,nst
        do 56 j = 1,i
   56     s(i,j) =s(j,i)
      return
c
c.... external load vector due to transversal constant load 
22    xsi = 0.d0
      eta = 0.d0
      call sfr2plbt(deriv,eta,xsi,nel,shp)
      call jacoplbt(cartd,deriv,djacb,xl,n,nel,ndm)
      da = djacb*4.d0
      do ino =1,4
        call qload33(qz,d,aqloa,numel,n,mqloa,propq,prop,isw)
        iq = (ino-1)*3+1
        p(iq) = p(iq) + shp(ino)*qz*da
      end do ! ino
      return
c
c
c.... formats
1001  format(' Input: E  nu  h  q rho ngb ngs',$)
1002  format(5x,'Materialdata Belyschtko/Tsay plate element:',/,
     + 5x,'elastic modulus............',   f15.6,/,
     + 5x,'poissons ratio................',f15.6,/,
     + 5x,'thickness.....................',f15.6,/,
     + 5x,'element load (transverse).....',f15.6,/,
     + 5x,'density ......................',f15.6,/,
     + 5x,'stability parameter r_w.......',f15.6,/,
     + 5x,'stability parameter r_beta....',f15.6)   
2001  format(a1,20a4,//,2x,'E L E M E N T   S T R E S S E S',//,
     1 2x,'EL',1x,'MAT',1x,'1-COORD',1x,'2-COORD',
     2 2X,'***MX***',3X,'***MY***',3X,'***MXY***',
     3 2X,'***QX***',3X,'***QY***',/)
2002  format(1x,i3,i3,2f8.3,1x,5e11.4)
      end
c----------------------------------------------------------------------
      subroutine bmatplbt(bb,bs,cartd,kk,shp,nel)
c***********************************************************************
c.....B-matrices for plate element
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension  bb(3,3),bs(2,3),cartd(2,nel),shp(nel)
      dnkdx = cartd(1,kk)
      dnkdy = cartd(2,kk)
c.....bb
        call  pzero(bb,9)
        bb(1,3) = dnkdx
        bb(2,2) =-dnkdy
        bb(3,3) = dnkdy
        bb(3,2) =-dnkdx
c.....bs
        call pzero(bs,6)
        bs(1,1) = dnkdx
        bs(1,3) = shp(kk)
        bs(2,1) = dnkdy
        bs(2,2) =-shp(kk)
      return
      end
c----------------------------------------------------------------------
      subroutine gradplbt(cartd,gr,ul,ndf,nel)
c***********************************************************************
c.....displacement gradients
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension cartd(2,nel),gr(6),ul(ndf,nel)
      call pzero(gr,6)
          do 20 ino = 1,nel
          dnidx = cartd(1,ino)
          dnidy = cartd(2,ino)
             do 20 idf = 1,ndf
             idfn = ndf + idf
             const = ul(idf,ino)
             gr(idf)  = gr(idf) + dnidx * const
   20        gr(idfn) = gr(idfn)+ dnidy * const
      return
      end
c----------------------------------------------------------------------
      subroutine strplbt(cartd,shp,ul,gr,db,ds,sig,nel,ndf)
c***********************************************************************
c.....Strains and stresses for plate
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension  cartd(2,nel),db(3,3),gr(6),ds(2,2),shp(nel),
     1           ul(ndf,nel),eb(3),es(2),sig(5)
      call pzero(sig,5)
c.....rotation at Gauss point
      xzrot=0.0
      yzrot=0.0
      nposn=ndf-1
      do 30 ino = 1,nel
         xzrot = xzrot+shp(ino)*ul(nposn,ino)
   30    yzrot = yzrot+shp(ino)*ul(ndf,ino)
c.....shear strains and shear forces
      es(1) = gr(1)+yzrot
      es(2) = gr(4)-xzrot
      sig(4)= ds(1,1)*es(1)
      sig(5)= ds(2,2)*es(2)
c.....bending strains and bending moments
      eb(1) = gr(3)
      eb(2) =-gr(5)
      eb(3) =-gr(2)+gr(6)
      sig(1)= db(1,1)*eb(1)+db(1,2)*eb(2)
      sig(2)= db(2,2)*eb(2)+db(2,1)*eb(1)
      sig(3)= db(3,3)*eb(3)
      return
      end
c----------------------------------------------------------------------
      subroutine jacoplbt(cartd,deriv,djacb,xl,n,nel,ndm)
c***********************************************************************
c.....Jacobi-Matrix and  cartesian derivatives of shape functions
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension xl(ndm,nel),cartd(2,nel),deriv(2,nel),
     +          xjaci(2,2),xjacm(2,2)
c.....Jacobi-matrix xjacm
        do 4 idm = 1,2
        do 4 jdm = 1,2
        xjacm(idm,jdm) = 0.0
          do 4 ino = 1,nel
    4     xjacm(idm,jdm)=xjacm(idm,jdm)+deriv(idm,ino)*xl(jdm,ino)
c.....Determinant  and Inverse of Jacobi-Matrix
        djacb=xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1)
        if(djacb) 6,6,8
    6     write(*,600) n
          stop
    8     xjaci(1,1)= xjacm(2,2)/djacb
          xjaci(2,2)= xjacm(1,1)/djacb
          xjaci(1,2)=-xjacm(1,2)/djacb
          xjaci(2,1)=-xjacm(2,1)/djacb
c.....cartesian derivatives
        do 10 idm = 1,2
          do 10 ino = 1,nel
          cartd(idm,ino) = 0.0
            do 10 jdm = 1,2
   10       cartd(idm,ino)=cartd(idm,ino)+xjaci(idm,jdm)*deriv(jdm,ino)
  600 format(//,' program halted in subroutine jacoplbt',/,11x,
     +          ' zero or negative area',/,10x,' element number',i5)
      return
      end
c----------------------------------------------------------------------
      subroutine mulplbt(bimat,bjvec,da,p,ndf,nst,ino,ncoli,nrowj,irow)
c***********************************************************************
c.....bi T * bj  * det j * weig(xsi) * weig(eta)
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension bimat(nrowj,ncoli),bjvec(nrowj),p(nst),sblod(3)
      ifrow = irow - 1
c.....bi T * bjv
        do 10 i = 1,ncoli
        sblod(i) = 0.0
          do 10 j = 1,nrowj
   10     sblod(i) = sblod(i) + bimat(j,i) * bjvec(j)
c.....store sblod in load vector
      ifrow = (ino-1)*ndf + ifrow
        do 30 i = 1,ncoli
        irsub = ifrow + i
   30   p(irsub) = p(irsub) - sblod(i)*da
      return
      end
c----------------------------------------------------------------------
      subroutine modplbt(db,ds,d)
c***********************************************************************
c.....Elasticity matrix for plate
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension db(3,3),ds(2,2),d(*)
      young = d(1)
      poiss = d(2)
      thick = d(3)
      cappa = 5.0d0/6.0d0
      const = (young*thick)/(1.0d0-poiss*poiss)
c.....Bending part
        const = const*thick*thick/12.0d0
        call pzero(db,9)
        db(1,1) = const
        db(2,2) = const
        db(1,2) = const*poiss
        db(2,1) = const*poiss
        db(3,3) = const*(1.0d0-poiss)*0.5d0
c.....Shear part
        const = 0.5d0*young*thick*cappa/(1.0d0+poiss)
        call pzero(ds,4)
        ds(1,1) = const
        ds(2,2) = const
      return
      end
c----------------------------------------------------------------------
      subroutine subplbt(bi,bj,da,dmat,s,ndf,nst,ino,jno,
     1 ncoli,nrowj,ncolj,irow,jcol)
c***********************************************************************
c.....bi T * d * bj * det j * weig(xsi) * weig(eta)
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension bi(nrowj,ncoli),bj(nrowj,ncolj),dmat(nrowj,nrowj),
     1 dbmat(3,3),s(nst,nst),sbstf(3,3)
      ifrow = irow - 1
      jfcol = jcol - 1
c.....d*bj
         do 10 j = 1,ncolj
         do 10 i = 1,nrowj
         dbmat(i,j) = 0.0
           do 10 k = 1,nrowj
   10      dbmat(i,j) = dbmat(i,j) + dmat(i,k) * bj(k,j)
c.....bi T *(d*bj)
         do 20 j = 1,ncolj
         do 20 i = 1,ncoli
         sbstf(i,j) = 0.0
           do 20 k = 1,nrowj
   20      sbstf(i,j) = sbstf(i,j) + bi(k,i) * dbmat(k,j)
c.....store into K_T
      ifrow = (ino-1)*ndf + ifrow
      jfcol = (jno-1)*ndf + jfcol
         do 30 i = 1,ncoli
         irsub = ifrow + i
           do 30 j = 1,ncolj
           jcsub = jfcol + j
   30      s(irsub,jcsub) = s(irsub,jcsub) + sbstf(i,j)*da
      return
       end
c----------------------------------------------------------------------
      subroutine sfr2plbt (deriv,eta,xsi,nel,shp)
c***********************************************************************
c.....Shape functions and derivatives for linear,quadratic
c.....lagrangian and serendipity isoparametric  2-d elements
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension deriv(2,nel),shp(nel)
      deins=1.0d0
      dzwei=2.0d0
      half=deins/dzwei
      viert=half/dzwei
      s = xsi
      t = eta
       st = s*t
c.....shp functions for 4-node element
        shp(1) = (deins-t-s+st)*viert
        shp(2) = (deins-t+s-st)*viert
        shp(3) = (deins+t+s+st)*viert
        shp(4) = (deins+t-s-st)*viert
c.....shp functions derivatives
        deriv(1,1) = (-deins+t)*viert
        deriv(1,2) = (+deins-t)*viert
        deriv(1,3) = (+deins+t)*viert
        deriv(1,4) = (-deins-t)*viert
        deriv(2,1) = (-deins+s)*viert
        deriv(2,2) = (-deins-s)*viert
        deriv(2,3) = (+deins+s)*viert
        deriv(2,4) = (+deins-s)*viert
      return
      end
c----------------------------------------------------------------------
      subroutine plotplbt(ix,dt,st,shp,sig,nel,da,numnp)
c***********************************************************************
c.....Plot   mx(1)   mxy(2)  my(3)      m1(5)  m2(6) phi_1(7)
c            qx(8)           qy(9)
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension dt(numnp),st(numnp,*),shp(9),sig(5),ix(1)
      do 10  j = 1,nel
        xsji = da*shp(j)
        ii = abs(ix(j))
        if(ii.eq.0) go to 10
        dt(ii) = dt(ii) + xsji
        st(ii,1) = st(ii,1) + sig(1)*xsji
        st(ii,2) = st(ii,2) + sig(3)*xsji
        st(ii,3) = st(ii,3) + sig(2)*xsji
        st(ii,8) = st(ii,8) + sig(4)*xsji
        st(ii,9) = st(ii,9) + sig(5)*xsji
10    continue
      return
      end
c----------------------------------------------------------------------
      subroutine gausplbt(ngs,pgp,wgp)
c***********************************************************************
c.....Gauss points  ngs = 1,2
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension pgp(2),wgp(2)
      do  i = 1,2
        pgp(i) = 0.0
        wgp(i) = 0.0
      end do
      if(ngs.gt.1) goto 2
      pgp(1) = 0.0
      wgp(1) = 2.0d0
      return
2     pgp(1) = -dsqrt(3.0d0)/3.0d0
      pgp(2) =  dsqrt(3.0d0)/3.0d0 
      wgp(1) =  1.0d0
      wgp(2) =  1.0d0
      return
      end
c
      subroutine qload33(qz,d,q,numel,n,mqloa,propq,prop,isw)
c----------------------------------------------------------
c.... set loads from macro qloa/mate
c----------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension q(numel,10),d(*)
      if(isw.eq.22) then
        qz=0.d0 
        if(mqloa.ne.1) qz = q(n,1)*propq 
      else
        qz = d(4)*prop
      end if
      return
      end
