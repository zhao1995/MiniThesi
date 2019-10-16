      subroutine elmt17(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c**********************************************************************|
c.....linear Dvorkin/Bathe plate element                               |
c     Version gamma_cart = J^-1*gamma_conv
c.... Material parameters (input)                                      |
c      1: E     = Young's modulus                                      |
c      2: nu    = Poisson's ratio                                      |
c      3: h     = thickness of shell                                   |
c      4: q     = transverse load                                      | 
c      5: rho   = density                                              | 
c      6: cappa = shear correction due to element size
c      signs of stresses, see lecture CTWM                             |
c      signs of displacements FEM-Like: w,phix,phiy                    |
c
c     card 2                                 
c     d(11) = f_cd                                     
c     d(12) = f_yd                                     
c     d(13) = d_x                                      
c     d(14) = d_y                                      


c      ng=1 for 4 prin stresses                                        | 
c      ng=2 for 8 plot stresses                                        | 
c                                                                      | 
c      w.wagner uka 10/03                   
c
c     provisorisch: 
c     Einbau Betonbemessung aus ELmt07, offen Bemessung für m_xy

c                          |
c**********************************************************************|
      USE bdata
      USE cdata
      USE eldata
      USE iofile
      USE pdata7
      USE pdata10
      USE pltran
      USE prisdat
      USE prlod
      USE qload      
      USE strnam 
      implicit double precision (a-h,o-z)
      dimension xl(ndm,4),tl(*),ix(*),ul(ndf,4),s(nst,nst),p(nst),d(*),
     1    pg(2),wg(2),gp(2),xjaci(2,2),
     2    cartd(2,4),deriv(2,4),shp(4),gr(6),sig(5),db(3,3),ds(2,2),
     3    bbi(3,3),bbj(3,3),bsi(2,3),bsj(2,3),yl(3,4),as(4)
      dimension h1(*),h2(*),h3(*)
c.... go to correct array processor
      go to(1,2,3,4,5,3,2,4,2,2,3,2,2,14,2,2,2,2,2,2,2,22) isw
c.... input material properties
1     if(ior.lt.0) write(*,1001)
      call dinput(d,6)
      if(d(6).eq.0.d0) d(6)=5.d0/6.d0  ! Def. cappa
      write(iow,1002) (d(i),i=1,6)
      if(ior.lt.0) write(*,1002) (d(i),i=1,6)
      if(ndm.eq.3) ipla = 1

c     d(11) = f_cd, d(12) = f_yd, d(13) = d_x, d(14) = d_y                                      
      call dinput(d(11),4)


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
      strsus(10) = 'STEEL as_x(Bot)'
      strsus(11) = 'STEEL as_x(Top)'                
      strsus(12) = 'STEEL as_y(Bot)'                
      strsus(13) = 'STEEL as_y(Top)'                
      do is =14,25
        strsus(is) = '               '
      end do
c...  names for principal moments
      nptyp = 3 
c
2     return
3     ng = 2
c.....dvorkin/bathe values
      call bdval(xl,ndm)
c.....gauss points
      call gausplbd(ng,pg,wg)
c.... elasticity matrix
      call shearfac17(d,xl,ndm,cappa)
      call modplbd(db,ds,d,cappa)
c.....stiffness matrix  
      do 30 ig = 1,ng
        do 30 jg = 1,ng
          xsi = pg(ig)
          eta = pg(jg)
c.....    shape functions,area etc.
          call sfr2plbd(deriv,xsi,eta,shp)
          call jacoplbd(cartd,deriv,djacb,xjaci,xl,n,ndm)
          da = djacb*wg(ig)*wg(jg)
c.....    displacement gradients
          call gradplbd(cartd,gr,ul,ndf)
c.....    stresses (M,Q)  and strains
          call strplbd(ul,shp,gr,db,ds,sig,ndf,xsi,eta,xjaci)
          do 31  ino =1,4
            call bmatplbd(bbi,bsi,cartd,deriv,xjaci,ino)
c.....      residual
            call mulplbd(bbi,sig(1),da,p,ndf,nst,ino,3,3,1)
            call mulplbd(bsi,sig(4),da,p,ndf,nst,ino,3,2,1)
c.....      external load vector due to transversal constant load 
            call qload17(qz,d,aqloa,numel,n,mqloa,propq,prop,isw)
            iq = (ino-1)*3+1
            p(iq) = p(iq) + shp(ino)*qz*da
            if(isw.eq.6) goto 31
            do 32 jno = ino,4
              call bmatplbd(bbj,bsj,cartd,deriv,xjaci,jno)
c.....        bbi T * db * bbj
              call subplbd(bbi,bbj,da,db,s,ndf,nst,ino,jno,3,3,3,1,1)
   32         call subplbd(bsi,bsj,da,ds,s,ndf,nst,ino,jno,3,2,3,1,1)
   31     continue
   30 continue
c.....lower part of K_T
      do 36 i = 1,nst
        do 36 j = 1,i
   36     s(i,j) =s(j,i)
c.....stress resultants 
4     istv = 9
      ng   = 1
      if(isw.eq.8) ng   = 2
      call bdval(xl,ndm)
      call gausplbd(ng,pg,wg)
c.... elasticity matrix
      call shearfac17(d,xl,ndm,cappa)
      call modplbd(db,ds,d,cappa)
      do 40 ig = 1,ng
      do 40 jg = 1,ng
        xsi = pg(ig)
        eta = pg(jg)
c.....  shape functions,area etc.
        call sfr2plbd(deriv,xsi,eta,shp)
        call jacoplbd(cartd,deriv,djacb,xjaci,xl,n,ndm)
        da = djacb*wg(ig)*wg(jg)
c.....  displacement gradients
        call gradplbd(cartd,gr,ul,ndf)
c.....  stresses (Q) and strains
        call strplbd(ul,shp,gr,db,ds,sig,ndf,xsi,eta,xjaci)
        if(isw.eq.4) then     
c.....    Output stresses  (M and Q)
c.....    coordinates
          do 41 idm = 1,2
            gp(idm) = 0.0
            do 41 ino = 1,4
   41         gp(idm) = gp(idm) + xl(idm,ino) * shp(ino)
            mct = mct - 1
            if(mct.gt.0) go to 42
c.....        head
                         write(iow,2001) o,head
            if(ior.lt.0) write(*  ,2001) o,head
            mct = 50
c.....      print values
   42                write(iow,2002) n,ma,(gp(i),i=1,2),(sig(i),i=1,5)
        if(ior.lt.0) write(*  ,2002) n,ma,(gp(i),i=1,2),(sig(i),i=1,5)
        elseif(isw.eq.8) then     
c....     stresses for plot M,Q
          if(iplma(ma).eq.0)  return ! only if MATN
c...      design for concrete
          if(d(11).gt.0.0d0) call design18(d,sig,as)
          call plotplbd(ix,strea,strea(1+numnp),shp,sig,as,da,numnp)
        endif  
   40 continue
      return
c
c.....mass matrix  
5     ng = 2
      dh = d(3)*d(3)/12.d0
c.....gauss points
      call gausplbd(ng,pg,wg)
c.... loop gauss points 
      do 50 ig = 1,ng
        do 50 jg = 1,ng
          xsi = pg(ig)
          eta = pg(jg)
c.....    shape functions,area etc.
          call sfr2plbd(deriv,xsi,eta,shp)
          call jacoplbd(cartd,deriv,djacb,xjaci,xl,n,ndm)
          da = djacb*wg(ig)*wg(jg)*d(3)*d(5)
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
   51     end do
   50 end do
c.....upper part of M
      do 56 i = 1,nst
        do 56 j = 1,i
   56     s(i,j) =s(j,i)
c
      return
c
c.....plot stress resultants from center gauss-point without averaging 
14    call gausplbd(1,pg,wg)
      call bdval(xl,ndm)
c.... elasticity matrix
      call shearfac17(d,xl,ndm,cappa)
      call modplbd(db,ds,d,cappa)
      xsi = pg(1)
      eta = pg(1)
c.... shape functions,area etc.
      call sfr2plbd(deriv,xsi,eta,shp)
      call jacoplbd(cartd,deriv,djacb,xjaci,xl,n,ndm)
      da = djacb*wg(1)*wg(1)
c.... displacement gradients
      call gradplbd(cartd,gr,ul,ndf)
c.... stresses 
      call strplbd(ul,shp,gr,db,ds,sig,ndf,xsi,eta,xjaci)
      if(nfp.eq.8) nfp=4 
      if(nfp.eq.9) nfp=5 
c
      if(flfp) then
c....   calculate extreme values
        xmaxf = max(xmaxf,sig(nfp))
        xminf = min(xminf,sig(nfp))
      else
c....   plot stresses elementwise in one color      
c....   calculate color      
        call pppcolf(sig(nfp))
c....   plot stress      
c....   coordinates
        call pzero(yl,12)
        do i=1,ndm
          do j=1,4
            yl(i,j)=xl(i,j) 
          end do
        end do
c..... transform for rot        
        call plxtrn(yl,tra,vr,3,4)
c...... plot element
        call plot9 (iel,ix,yl,ndm,4,1)
      endif
      return
c      

c.....external load vector due to transversal constant load 
22    ng = 2
      call gausplbd(ng,pg,wg)
      do ig = 1,ng
        do jg = 1,ng
          xsi = pg(ig)
          eta = pg(jg)
          call sfr2plbd(deriv,xsi,eta,shp)
          call jacoplbd(cartd,deriv,djacb,xjaci,xl,n,ndm)
          da = djacb*wg(ig)*wg(jg)
          do ino =1,4
            call qload17(qz,d,aqloa,numel,n,mqloa,propq,prop,isw)
            iq = (ino-1)*3+1
            p(iq) = p(iq) + shp(ino)*qz*da
          end do ! ino
        end do ! jg
      end do ! ig
      return
c 
c.... formats
1001  format(' Input: E  nu  h  q',$)
1002  format(5x,'Material data B/D plate element:',/,
     + 5x,'elastic modulus............',   f15.4,/,
     + 5x,'poissons ratio................',f12.4,/,
     + 5x,'thickness.....................',f12.4,/,
     + 5x,'element load (transverse).....',f12.4,/,
     + 5x,'density ......................',f12.4,/,
     + 5x,'SCF(<0:|VAL|*FE-SCF,>0:VAL)...',f12.4)
2001  format(a1,20a4,//,2x,'E L E M E N T   S T R E S S E S',//,
     1 2x,'EL',1x,'MAT',1x,'1-COORD',1x,'2-COORD',
     2 2X,'***MX***',3X,'***MY***',3X,'***MXY***',
     3 2X,'***QX***',3X,'***QY***',/)
2002  format(1x,i3,i3,2f8.3,1x,5e11.4)
      end
c----------------------------------------------------------------------
      subroutine bmatplbd(bb,bs,cartd,deriv,xjaci,i)
c***********************************************************************
c.....B-matrices for plate element
c***********************************************************************
      implicit double precision (a-h,o-z)
      common /bdval17/ xxsim(4),yxsim(4),xetal(4),yetal(4)
!$OMP THREADPRIVATE (/bdval17/)  
      dimension  bb(3,3),bs(2,3),cartd(2,4),xsii(4),etai(4),deriv(2,4),
     + xjaci(2,2),bsc(2,3) 
      data xsii /-1.d0, 1.d0,1.d0,-1.d0/
      data etai /-1.d0,-1.d0,1.d0, 1.d0/
      dnkdx = cartd(1,i)
      dnkdy = cartd(2,i)
      dnkds = deriv(1,i)
      dnkdt = deriv(2,i)
c
c.....bb
      call  pzero(bb,9)
      bb(1,3) =  dnkdx
      bb(2,2) = -dnkdy
      bb(3,3) =  dnkdy
      bb(3,2) = -dnkdx
c.....bs convective
      call pzero(bs,6)
      bsc(1,1) = dnkds                    
      bsc(1,3) = dnkds*xsii(i)*xxsim(i)
      bsc(1,2) =-dnkds*xsii(i)*yxsim(i)
      bsc(2,1) = dnkdt                 
      bsc(2,3) = dnkdt*etai(i)*xetal(i)
      bsc(2,2) =-dnkdt*etai(i)*yetal(i)
c.....bs cartesian gamma = J^-1 * gamma_conv
      call pzero(bs,6)
      do ii=1,2
        do k=1,3
          do j=1,2    
            bs(ii,k)=bs(ii,k)+xjaci(ii,j)*bsc(j,k)
          end do
        end do  
      end do 
c
      return
      end
c----------------------------------------------------------------------
      subroutine gradplbd(cartd,gr,ul,ndf)
c***********************************************************************
c.....displacement gradients
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension cartd(2,4),gr(6),ul(ndf,4)
      call pzero(gr,6)
          do 20 ino = 1,4
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
      subroutine strplbd(ul,shp,gr,db,ds,sig,ndf,xsi,eta,xjaci)
c***********************************************************************
c.....Strains and stresses for b/d-plate
c***********************************************************************
      implicit double precision (a-h,o-z)
      common /bdval17/ xxsim(4),yxsim(4),xetal(4),yetal(4)
!$OMP THREADPRIVATE (/bdval17/)  
      dimension  db(3,3),ds(2,2),gr(6),xjaci(2,2),ul(ndf,4),
     + eb(3),sig(5),shp(4)
      call pzero(sig,5)
c.....bending strains and bending moments
      eb(1) = gr(3)
      eb(2) =-gr(5)
      eb(3) =-gr(2)+gr(6)
      sig(1) = db(1,1)*eb(1)+db(1,2)*eb(2)
      sig(2) = db(2,2)*eb(2)+db(2,1)*eb(1)
      sig(3) = db(3,3)*eb(3)
c
c.....compatible shear strains in M,L
      esxsiB = (xxsim(1)*(ul(3,1) + ul(3,2))
     +       -  yxsim(1)*(ul(2,1) + ul(2,2)) + ul(1,2)-ul(1,1) )/2.d0  
      esxsiD = (xxsim(3)*(ul(3,3) + ul(3,4))
     +       -  yxsim(3)*(ul(2,3) + ul(2,4)) + ul(1,3)-ul(1,4) )/2.d0  
      esetaA = (xetal(1)*(ul(3,1) + ul(3,4))
     +       -  yetal(1)*(ul(2,1) + ul(2,4)) + ul(1,4)-ul(1,1) )/2.d0  
      esetaC = (xetal(2)*(ul(3,2) + ul(3,3))
     +       -  yetal(2)*(ul(2,2) + ul(2,3)) + ul(1,3)-ul(1,2) )/2.d0  
c.... convective shear strains Es
      esxsi = 0.5d0*( (1.d0-eta)*esxsiB + (1.d0+eta)*esxsiD )
      eseta = 0.5d0*( (1.d0-xsi)*esetaA + (1.d0+xsi)*esetaC )
c.....cartesian shear strains Ex=J-1*Es  
      esx = xjaci(1,1)*esxsi + xjaci(1,2)*eseta
      esy = xjaci(2,1)*esxsi + xjaci(2,2)*eseta
c.....cartesian shear forces Qx=Gx*Ex  
      sig(4)= ds(1,1)*esx 
      sig(5)= ds(2,2)*esy
      return
      end
c----------------------------------------------------------------------
      subroutine jacoplbd(cartd,deriv,djacb,xjaci,xl,n,ndm)
c***********************************************************************
c.....Jacobi-Matrix and  cartesian derivatives of shape functions
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension xl(ndm,4),cartd(2,4),deriv(2,4),xjaci(2,2),xjacm(2,2)
c.....Jacobi-matrix xjacm
        do 4 idm = 1,2
        do 4 jdm = 1,2
        xjacm(idm,jdm) = 0.0
          do 4 ino = 1,4
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
          do 10 ino = 1,4
          cartd(idm,ino) = 0.0
            do 10 jdm = 1,2
   10       cartd(idm,ino)=cartd(idm,ino)+xjaci(idm,jdm)*deriv(jdm,ino)
  600 format(' SR jacoplbd: zero or negative area for element ',i5)
      return
      end
c----------------------------------------------------------------------
      subroutine mulplbd(bimat,bjvec,da,p,ndf,nst,ino,ncoli,nrowj,irow)
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
      subroutine modplbd(db,ds,d,cappa)
c***********************************************************************
c.....Elasticity matrix for plate in cartesian coordinates
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension db(3,3),ds(2,2),d(*)
      young = d(1)
      poiss = d(2)
      thick = d(3)
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
      call pzero(ds,4)
      const = 0.5d0*young*thick*cappa/(1.0d0+poiss)
      ds(1,1) = const
      ds(2,2) = const
      return
      end
c----------------------------------------------------------------------
      subroutine subplbd(bi,bj,da,dmat,s,ndf,nst,ino,jno,
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
      subroutine sfr2plbd (deriv,xsi,eta,shp)
c***********************************************************************
c.....Shape functions and derivatives for linear isoparametric  2-d elements
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension deriv(2,4),shp(4)
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
      subroutine gausplbd(ng,pg,wg)
c***********************************************************************
c.....Gauss points  ng = 1,2
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension pg(2),wg(2)
      do 10 i = 1,2
        pg(i) = 0.0
   10   wg(i) = 0.0
      if(ng.gt.1) goto 2
        pg(1) = 0.0
        wg(1) = 2.0d0
        return
    2 continue
        pg(1) = -dsqrt(3.0d0)/3.0d0
        wg(1) =  1.0d0
        pg(2) = -pg(1)
        wg(2) =  wg(1)
      return
      end
c----------------------------------------------------------------------
      subroutine plotplbd(ix,dt,st,shp,sig,as,da,numnp)
c***********************************************************************
c.....Plot   mx(1)   mxy(2)  my(3)      m1(5)  m2(6) phi_1(7)
c            qx(8)           qy(9)
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension dt(numnp),st(numnp,*),shp(9),sig(5),ix(1),as(4)
        do 10  j = 1,4
           xsji = da*shp(j)
           ii = abs(ix(j))
           if(ii.eq.0) go to 10
             dt(ii) = dt(ii) + xsji
             st(ii, 1) = st(ii, 1) + sig(1)*xsji
             st(ii, 2) = st(ii, 2) + sig(3)*xsji
             st(ii, 3) = st(ii, 3) + sig(2)*xsji
             st(ii, 8) = st(ii, 8) + sig(4)*xsji
             st(ii, 9) = st(ii, 9) + sig(5)*xsji
             st(ii,10) = st(ii,10) + as (2)*xsji
             st(ii,11) = st(ii,11) + as (1)*xsji
             st(ii,12) = st(ii,12) + as (4)*xsji
             st(ii,13) = st(ii,13) + as (3)*xsji
10      continue
      return
      end
      subroutine bdval(xl,ndm)
c***********************************************************************
c.....Bathe/Dvorkin-values
c***********************************************************************
      implicit double precision (a-h,o-z)
      common /bdval17/ xxsim(4),yxsim(4),xetal(4),yetal(4)
!$OMP THREADPRIVATE (/bdval17/)  
      dimension xl(ndm,4)
c.... at M: B,D for all nodes BBDD
      xxsim(1)= 0.5d0*(xl(1,2)-xl(1,1))
      xxsim(2)= xxsim(1)
      xxsim(3)= 0.5d0*(xl(1,3)-xl(1,4))
      xxsim(4)= xxsim(3)
c
      yxsim(1)= 0.5d0*(xl(2,2)-xl(2,1))
      yxsim(2)= yxsim(1)
      yxsim(3)= 0.5d0*(xl(2,3)-xl(2,4))
      yxsim(4)= yxsim(3)
c
c.... at L: A,C for all nodes ACCA
      xetal(1)= 0.5d0*(xl(1,4)-xl(1,1))
      xetal(2)= 0.5d0*(xl(1,3)-xl(1,2))
      xetal(3)= xetal(2)
      xetal(4)= xetal(1)
c
      yetal(1)= 0.5d0*(xl(2,4)-xl(2,1))
      yetal(2)= 0.5d0*(xl(2,3)-xl(2,2))
      yetal(3)= yetal(2)
      yetal(4)= yetal(1)
c
      return
      end
c
      subroutine qload17(qz,d,q,numel,n,mqloa,propq,prop,isw)
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
c
      subroutine design17(d,sig,as)
c--------------------------------------------------------------------+
c     concrete design  (from general design diagram) DIN 1045-1      |
c     mue_Ed = m / (b*d**2*f_cd)                                     |
c     b      = 1                                                     |
c     f_cd   = alpha *f_ck/gamma_c                                   |
c     zeta   = 0.5 + sqrt(0.25 - 0.5*mue_Ed)       (Approximation)   |
c     z      = zeta*d                                                |  
c     ohne Druckbewehrung                                            |
c     as1    = m_Ed /(z*sig_s1d)                                     |
c                                                                    |
c     mit Druckbewehrung                                             |
c     m_Ed_lim   = mue_Ed_lim*(b*d**2*f_cd)                          |
c     Delta m_Ed = m_ed - m_Ed_lim                                   |
c     as1    = m_Ed_lim /(z*sig_s1d)+[Delta m_Ed /((d-d_2)*sig_s1d)] |
c     as2    =                       [Delta m_Ed /((d-d_2)*sig_s2d)] |
c                                                                    |
c     f_yd   = f_yk/gamma_s                                          |
c     Sigma_sd = f_yd                    für mue_Ed < mue_Ed_lim     |
c     Sigma_sd = f_yd(4.906-10.5 mue_Ed) für mue_Ed > mue_Ed_lim     |
c     mue_Ed_lim = 0.372                                             | 
c                                                                    |
c     load and material are given with respect to ultimate load
c     e.g. gamma_Q =1.35
c          gamma_s =1.15 -> f_yd
c          gamma_c =1.5  -> f_cd
c--------------------------------------------------------------------+
      implicit double precision (a-h,o-z)
      dimension d(*),sig(3),as(4)
      gmued  = 0.372d0
      f_cd   = d(11)
      f_yd   = d(12)
      dx     = d(13)
      dy     = d(14)
      h      = d(3)
      asxadd = 0.d0
      asyadd = 0.d0
      call pzero(as,4)
c.... design
c.....m_x, m_xy added not correct, see EC2  
      if(sig(1) .ge. 0.0d0) smx=      sig(1)  + dabs(sig(3))  
      if(sig(1) .lt. 0.0d0) smx= dabs(sig(1)) + dabs(sig(3)) 
      xmued = smx /(dx*dx*f_cd)
      val   = 0.25d0 - 0.5d0*xmued
      if(val.lt.0.d0) then ! provisorial 
        write(*,*) 'zeta  not available '  
        val=0.d0
      end if
      zeta  = 0.5d0  + dsqrt(val)
      sigma_sd = f_yd 
      if(xmued .gt. gmued) then
        sigma_sd = f_yd*(4.906-10.5*xmued)
        if(sigma_sd.lt.0.d0) then ! provisorial
          write(*,*)  'Sigma_sd < 0, load too high '  
          sigma_sd=1.e-3
        end if  
        dsmx  =  (xmued - gmued)*dx*dx*f_cd
        smx   =  smx  - dsmx
c....   approx. distance for steel d-d_2  
        dmd2x  = dx - (h-dx)
        asxadd = dsmx/(sigma_sd*dmd2x)
      end if
      asx = smx / (zeta*dx*sigma_sd)
      if(sig(1) .ge. 0.0d0) then
        as(2) = asx + asxadd
        as(1) =       asxadd
      elseif(sig(1) .lt. 0.0d0) then
        as(2) =       asxadd
        as(1) = asx + asxadd
      end if
c.....m_y m_xy added not correct, see EC2
      if(sig(2) .ge. 0.0d0) smy=      sig(2)  + dabs(sig(3))
      if(sig(2) .lt. 0.0d0) smy= dabs(sig(2)) + dabs(sig(3))
      ymued = smy /(dy*dy*f_cd)
      val   = 0.25d0 - 0.5d0*ymued
      if(val.lt.0.d0) then ! provisorial 
        write(*,*) 'zeta  not available '  
        val=0.d0
      end if
      zeta  = 0.5d0  + dsqrt(val)
      sigma_sd = f_yd 
      if(ymued .gt. gmued) then
        sigma_sd = f_yd*(4.906-10.5*ymued)
        if(sigma_sd.lt.0.d0) then ! provisorial
          write(*,*)  'Sigma_sd < 0, load too high '  
          sigma_sd=1.e-3
        end if  
        dsmy  =  (ymued - gmued)*dy*dy*f_cd
        smy   =  smy  - dsmy
c....   approx. distance for steel d-d_2  
        dmd2y  = dy - (h-dy)
        asyadd = dsmy/(sigma_sd*dmd2y)
      end if
      asy = smy / (zeta*dy*sigma_sd)
      if(sig(2) .ge. 0.0d0) then
        as(4) = asy + asyadd
        as(3) =       asyadd
      elseif(sig(2) .lt. 0.0d0) then
        as(4) =       asyadd
        as(3) = asy + asyadd
      end if
      return
      end
c
      subroutine shearfac17(d,xl,ndm,cappa)
c-----------------------------------------------------------------------
c.... shear correction factor due to size
c.... shear correction factor
c     d(6) <0: cappa = input value*SCF(Tessler) 
c     d(6) >0: cappa = input value 
c     d(6) =0: cappa = 1 (default)   
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  d(*),xl(ndm,*),yl(ndm)

      cappa = d(6)
      
      if(cappa.lt.0.d0)then
 
c....   Ae ... (longest element side)^2 
        xnu = d(2)   
        hs  = d(3)
        call pzero(yl,3) 
        ae = 0.d0
        do i = 1,4
          k = i+1
          if(i.eq.4) k=1
          do j = 1,ndm
            yl(j) = xl(j,k)-xl(j,i)
          end do 
          sl2 = dot(yl,yl,ndm)       
          ae = max(ae,sl2)
        end do  
        cappa  = dabs(d(6))/(1.d0+ae/(2.d0*hs*hs*(1.d0+xnu))) ! Tessler
      end if 
c      
      return
      end  
c
