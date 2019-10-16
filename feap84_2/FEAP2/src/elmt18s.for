      subroutine elmt18(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c-----------------------------------------------------------------------
c.....linear Reissner Mindlin plate element 'SRI' 4/8/9 nodes     
c.... Material parameters (input)
c     card 1                                 
c     d( 1) E     = Young's modulus                                 
c     d( 2) nu    = Poisson's ratio                                 
c     d( 3) h     = thickness of shell                              
c     d( 4) q     = transverse load                                 
c     d( 5) rho   = density                                         
c     d( 6) ngb   = gauss points bending def=2                      
c     d( 7) ngs   = gauss points shear   def=1                      
c     d( 8) cappa = shear correction due to element size
c
c     card 2                                 
c     d(11) = f_cd                                     
c     d(12) = f_yd                                     
c     d(13) = d_x                                      
c     d(14) = d_y                                      
c                                                                 
c      Be careful, this element  may lead to hourglass-modes for  
c      thin plates using a 2/1 or 3/2 integration.                
c                                                                 
c      thin plate: L/h>10 thick plate: L/h<= 10                   
c                                                                 
c      Integration: thin plate 2/1, thick plate 2/2               
c                              3/2              3/3               
c                                                                 
c      signs of stresses, see lecture CTWM                        
c      signs of displacements FEM-Like w,phix,phiy                
c                                                                 
c
c     provisorisch: 
c     Einbau Betonbemessung aus ELmt07, offen Bemessung für m_xy
c                                                                        
c      modified+updated w.wagner uka 10/03                            
c-----------------------------------------------------------------------
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
      dimension xl(ndm,*),tl(*),ix(*),ul(ndf,*),s(nst,nst),p(nst),d(*),
     1    pgb(4),wgb(4),pgs(4),wgs(4),gp(2),
     2    cartd(2,9),deriv(2,9),shp(9),gr(6),sig(5),db(3,3),ds(2,2),
     3    bbi(3,3),bbj(3,3),bsi(2,3),bsj(2,3),yl(3,4),as(4)
      dimension h1(*),h2(*),h3(*)
c.... go to correct array processor
      go to(1,2,3,4,5,3,2,4,2,2,3,2,2,14,2,2,2,2,2,2,2,22) isw
c.... input material properties
1     if(ior.lt.0) write(*,1001)
      call dinput(d,8)
      if(d(6).eq.0.d0) d(6)=2.d0
      if(d(7).eq.0.d0) d(7)=1.d0
      if(d(8).eq.0.d0) d(8)=5.d0/6.d0
                   write(iow,1002) (d(i),i=1,8)
      if(ior.lt.0) write(*  ,1002) (d(i),i=1,8)
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
3     ngb = d(6)
      call gauspl18(ngb,pgb,wgb)
c.... elasticity matrix
      call shearfac18(d,xl,ndm,cappa)
      call modpl18(db,ds,d,cappa)
c.....bending terms in K  
      do 30 igb = 1,ngb
        do 30 jgb = 1,ngb
          xsi = pgb(igb)
          eta = pgb(jgb)
c.....    shape functions,area etc.
          call sfr2pl18(deriv,eta,xsi,nel,shp)
          call jacopl18(cartd,deriv,djacb,xl,n,nel,ndm)
          da = djacb*wgb(igb)*wgb(jgb)
c.....    displacement gradients
          call gradpl18(cartd,gr,ul,ndf,nel)
c.....    stresses (M)  and strains
          call strpl18(cartd,shp,ul,gr,db,ds,sig,nel,ndf)
          do 31  ino =1,nel
            call bmatpl18(bbi,bsi,cartd,ino,shp,nel,1,0)
c.....      residual
            call mulpl18(bbi,sig(1),da,p,ndf,nst,ino,3,3,1)
c.....      external load vector due to transversal constant load 
            call qload18(qz,d,aqloa,numel,n,mqloa,propq,prop,isw)
            iq = (ino-1)*3+1
            p(iq) = p(iq) + shp(ino)*qz*da
            if(isw.eq.6) goto 31
            do 32 jno = ino,nel
              call bmatpl18(bbj,bsj,cartd,jno,shp,nel,1,0)
c.....        bbi T * db * bbj
   32         call subpl18(bbi,bbj,da,db,s,ndf,nst,ino,jno,3,3,3,1,1)
   31     continue
   30 continue
c.....Shear terms in K  
      ngs = d(7)
      call gauspl18(ngs,pgs,wgs)
      do 33 igs = 1,ngs
      do 33 jgs = 1,ngs
        xsi = pgs(igs)
        eta = pgs(jgs)
c.....  shape functions,area etc.
        call sfr2pl18(deriv,eta,xsi,nel,shp)
        call jacopl18(cartd,deriv,djacb,xl,n,nel,ndm)
        da = djacb*wgs(igs)*wgs(jgs)
c.....  displacement gradients
        call gradpl18(cartd,gr,ul,ndf,nel)
c.....  stresses (Q) and strains
        call strpl18(cartd,shp,ul,gr,db,ds,sig,nel,ndf)
c.....  bsi T * ds * bsj
        do 34 ino = 1,nel
          call bmatpl18(bbi,bsi,cartd,ino,shp,nel,0,1)
          call mulpl18(bsi,sig(4),da,p,ndf,nst,ino,3,2,1)
          if(isw.eq.6) goto 34
          do 35 jno = ino,nel
              call bmatpl18(bbj,bsj,cartd,jno,shp,nel,0,1)
   35         call subpl18(bsi,bsj,da,ds,s,ndf,nst,ino,jno,3,2,3,1,1)
   34   continue
   33 continue
c.....lower part of K_T
      do 36 i = 1,nst
        do 36 j = 1,i
   36     s(i,j) =s(j,i)
c.....stress resultants 
4     istv = 9
      ngs = d(7)
      call gauspl18(ngs,pgs,wgs)
c.... elasticity matrix
      call shearfac18(d,xl,ndm,cappa)
      call modpl18(db,ds,d,cappa)
      do 40 igs = 1,ngs
      do 40 jgs = 1,ngs
        xsi = pgs(igs)
        eta = pgs(jgs)
c.....  shape functions,area etc.
        call sfr2pl18(deriv,eta,xsi,nel,shp)
        call jacopl18(cartd,deriv,djacb,xl,n,nel,ndm)
        da = djacb*wgs(igs)*wgs(jgs)
c.....  displacement gradients
        call gradpl18(cartd,gr,ul,ndf,nel)
c.....  stresses (Q) and strains
        call strpl18(cartd,shp,ul,gr,db,ds,sig,nel,ndf)
        if(isw.eq.4) then     
c.....    Output stresses  (M and Q)
c.....    coordinates
          do 41 idm = 1,2
            gp(idm) = 0.0
            do 41 ino = 1,nel
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
          call plotpl18(ix,strea,strea(1+numnp),shp,sig,as,nel,da,
     +    numnp)
        endif  
   40 continue
      return
c
c.....mass matrix  
5     ngb = 2
      dh = d(3)*d(3)/12.d0
c.....gauss points
      call gauspl18(ngb,pgb,wgb)
c.... loop gauss points 
      do 50 igb = 1,ngb
        do 50 jgb = 1,ngb
          xsi = pgb(igb)
          eta = pgb(jgb)
c.....    shape functions,area etc.
          call sfr2pl18(deriv,eta,xsi,nel,shp)
          call jacopl18(cartd,deriv,djacb,xl,n,nel,ndm)
          da = djacb*wgb(igb)*wgb(jgb)*d(3)*d(5)
          do 51  ino =1,4
            j=(ino-1)*ndf+1
            p(j  ) = p(j) + shp(ino)*da
            p(j+1) = p(j)*dh  
            p(j+2) = p(j)*dh  
   51   end do
   50 end do
c.....lower part of K_T
c      do 56 i = 1,nst
c        do 56 j = 1,i
c   56     s(i,j) =s(j,i)
c
c
      return
c
c.....plot stress resultants from center gauss-point without averaging 
14    call gauspl18(1,pgs,wgs)
c.... elasticity matrix
      call shearfac18(d,xl,ndm,cappa)
      call modpl18(db,ds,d,cappa)
      xsi = pgs(1)
      eta = pgs(1)
c.... shape functions,area etc.
      call sfr2pl18(deriv,eta,xsi,nel,shp)
      call jacopl18(cartd,deriv,djacb,xl,n,nel,ndm)
      da = djacb*wgs(1)*wgs(1)
c.... displacement gradients
      call gradpl18(cartd,gr,ul,ndf,nel)
c.... stresses 
      call strpl18(cartd,shp,ul,gr,db,ds,sig,nel,ndf)
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
        call plot9 (iel,ix,yl,ndm,nel,1)
      endif
      return

c.....external load vector due to transversal constant load 
22    ngb = d(6)
      call gauspl18(ngb,pgb,wgb)

      do igb = 1,ngb
        do jgb = 1,ngb
          xsi = pgb(igb)
          eta = pgb(jgb)
c.....    shape functions,area etc.
          call sfr2pl18(deriv,eta,xsi,nel,shp)
          call jacopl18(cartd,deriv,djacb,xl,n,nel,ndm)
          da = djacb*wgb(igb)*wgb(jgb)
          do ino =1,nel
            call qload18(qz,d,aqloa,numel,n,mqloa,propq,prop,isw)
            iq = (ino-1)*3+1
            p(iq) = p(iq) + shp(ino)*qz*da
          end do ! ino
        end do ! jgb
      end do ! igb
      return
c
c.... formats
1001  format(' Input: E  nu  h  q rho ngb ngs cappa',$)
1002  format(5x,'Materialdata SRI plate element 18:',/,
     + 5x,'elastic modulus............',   f15.4,/,
     + 5x,'poissons ratio................',f12.4,/,
     + 5x,'thickness.....................',f12.4,/,
     + 5x,'element load (transverse).....',f12.4,/,
     + 5x,'density ......................',f12.4,/,
     + 5x,'number gauss points bending...',f12.4,/,
     + 5x,'number gauss points shear ....',f12.4,/,
     + 5x,'SCF(<0:|VAL|*FE-SCF,>0:VAL)...',f12.4)
2001  format(a1,20a4,//,2x,'E L E M E N T   S T R E S S E S',//,
     1 2x,'EL',1x,'MAT',1x,'1-COORD',1x,'2-COORD',
     2 2X,'***MX***',3X,'***MY***',3X,'***MXY***',
     3 2X,'***QX***',3X,'***QY***',/)
2002  format(1x,i3,i3,2f8.3,1x,5e11.4)
      end
c----------------------------------------------------------------------
      subroutine bmatpl18(bb,bs,cartd,kk,shp,nel,ifb,ifs)
c***********************************************************************
c.....B-matrices for plate element
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension  bb(3,3),bs(2,3),cartd(2,nel),shp(nel)
      dnkdx = cartd(1,kk)
      dnkdy = cartd(2,kk)
c.....bb
      if(ifb.eq.1) then
        call  pzero(bb,9)
        bb(1,3) = dnkdx
        bb(2,2) =-dnkdy
        bb(3,3) = dnkdy
        bb(3,2) =-dnkdx
      endif
c.....bs
      if(ifs.eq.1) then
        call pzero(bs,6)
        bs(1,1) = dnkdx
        bs(1,3) = shp(kk)
        bs(2,1) = dnkdy
        bs(2,2) =-shp(kk)
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine gradpl18(cartd,gr,ul,ndf,nel)
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
      subroutine strpl18(cartd,shp,ul,gr,db,ds,sig,nel,ndf)
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
      es(1) = (gr(1)+yzrot)
      es(2) = (gr(4)-xzrot)
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
      subroutine jacopl18(cartd,deriv,djacb,xl,n,nel,ndm)
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
  600 format(//,' program halted in subroutine jacopl18',/,11x,
     +          ' zero or negative area',/,10x,' element number',i5)
      return
      end
c----------------------------------------------------------------------
      subroutine mulpl18(bimat,bjvec,da,p,ndf,nst,ino,ncoli,nrowj,irow)
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
      subroutine modpl18(db,ds,d,cappa)
c***********************************************************************
c.....Elasticity matrix for plate
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
        const = 0.5d0*young*thick*cappa/(1.0d0+poiss)
        call pzero(ds,4)
        ds(1,1) = const
        ds(2,2) = const
      return
      end
c----------------------------------------------------------------------
      subroutine subpl18(bi,bj,da,dmat,s,ndf,nst,ino,jno,
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
      subroutine sfr2pl18 (deriv,eta,xsi,nel,shp)
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
      if(nel.gt.4) goto 10
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
   10 if(nel.gt.8) goto 20
        s2 = s*dzwei
        t2 = t*dzwei
        ss = s*s
        tt = t*t
        st = s*t
        sst = s*s*t
        stt = s*t*t
        st2 = s*t*dzwei
c.....shp functions for 8-node element
        shp(1) = (-deins+st+ss+tt-sst-stt)*viert
        shp(2) = (-deins-st+ss+tt-sst+stt)*viert
        shp(3) = (-deins+st+ss+tt+sst+stt)*viert
        shp(4) = (-deins-st+ss+tt+sst-stt)*viert
        shp(5) = (+deins-t-ss+sst)*half
        shp(6) = (+deins+s-tt-stt)*half
        shp(7) = (+deins+t-ss-sst)*half
        shp(8) = (+deins-s-tt+stt)*half
c.....shp functions derivatives
        deriv(1,1) = (+t+s2-st2-tt)*viert
        deriv(1,2) = (-t+s2-st2+tt)*viert
        deriv(1,3) = (+t+s2+st2+tt)*viert
        deriv(1,4) = (-t+s2+st2-tt)*viert
        deriv(2,1) = (+s+t2-ss-st2)*viert
        deriv(2,2) = (-s+t2-ss+st2)*viert
        deriv(2,3) = (+s+t2+ss+st2)*viert
        deriv(2,4) = (-s+t2+ss-st2)*viert
        deriv(1,5) =  -s+st
        deriv(1,7) =  -s-st
        deriv(2,6) =  -t-st
        deriv(2,8) =  -t+st
        deriv(1,6) = (+deins-tt)*half
        deriv(1,8) = (-deins+tt)*half
        deriv(2,5) = (-deins+ss)*half
        deriv(2,7) = (+deins-ss)*half
      return
  20  continue
        ss=s*s
        st=s*t
        tt=t*t
        s1=s+deins
        t1=t+deins
        s2=s*dzwei
        t2=t*dzwei
        s9=s-deins
        t9=t-deins
c.....shp functions and derivatives for 9 noded element
        shp(1)=viert*s9*st*t9
        shp(5)=half*(deins-ss)*t*t9
        shp(2)=viert*s1*st*t9
        shp(6)=half*s*s1*(deins-tt)
        shp(3)=viert*s1*st*t1
        shp(7)=half*(deins-ss)*t*t1
        shp(4)=viert*s9*st*t1
        shp(8)=half*s*s9*(deins-tt)
        shp(9)=(deins-ss)*(deins-tt)
        deriv(1,1)=viert*t*t9*(-deins+s2)
        deriv(1,5)=-st*t9
        deriv(1,2)=viert*(deins+s2)*t*t9
        deriv(1,6)=half*(deins+s2)*(deins-tt)
        deriv(1,3)=viert*(deins+s2)*t*t1
        deriv(1,7)=-st*t1
        deriv(1,4)=viert*(-deins+s2)*t*t1
        deriv(1,8)=half*(-deins+s2)*(deins-tt)
        deriv(1,9)=-s2*(deins-tt)
        deriv(2,1)=viert*(-deins+t2)*s*s9
        deriv(2,5)=half*(deins-ss)*(-deins+t2)
        deriv(2,2)=viert*s*s1*(-deins+t2)
        deriv(2,6)=-st*s1
        deriv(2,3)=viert*s*s1*(deins+t2)
        deriv(2,7)=half*(1-ss)*(deins+t2)
        deriv(2,4)=viert*s*s9*(deins+t2)
        deriv(2,8)=-st*s9
        deriv(2,9)=-t2*(deins-ss)
      return
      end
c----------------------------------------------------------------------
      subroutine gauspl18(ngs,pgp,wgp)
c***********************************************************************
c.....Gauss points  ngs = 1,3
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension pgp(4),wgp(4)
        do 10 i = 1,4
        pgp(i) = 0.0
   10   wgp(i) = 0.0
      if(ngs.gt.1) goto 2
        pgp(1) = 0.0
        wgp(1) = 2.0d0
        goto 16
    2 continue
      if(ngs.gt.2) goto 4
        pgp(1) = -dsqrt(3.0d0)/3.0d0
        wgp(1) =  1.0d0
      goto 6
    4   pgp(1) = -dsqrt(0.6d0)
        pgp(2) =  0.0d0
        wgp(1) = 5.0d0/9.0d0
        wgp(2) = 8.0d0/9.0d0
    6 kgs = ngs/2
        do 8 igash = 1,kgs
        jgash = ngs + 1 - igash
        pgp(jgash) = -pgp(igash)
    8   wgp(jgash) =  wgp(igash)
   16 return
      end
c----------------------------------------------------------------------
      subroutine plotpl18(ix,dt,st,shp,sig,as,nel,da,numnp)
c***********************************************************************
c.....Plot   mx(1)   mxy(2)  my(3)      m1(5)  m2(6) phi_1(7)
c            qx(8)           qy(9)
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension dt(numnp),st(numnp,*),shp(9),sig(5),ix(1),as(4)
        do 10  j = 1,nel
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
c
      subroutine qload18(qz,d,q,numel,n,mqloa,propq,prop,isw)
c----------------------------------------------------------
c.... set loads from macro qloa/mate
c----------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension q(numel,10),d(*)
      if(isw.eq.22) then
        qz = 0.d0 
        if(mqloa.ne.1) qz = q(n,1)*propq 
      else
        qz = d(4)*prop 
      end if 
      return
      end
c
      subroutine design18(d,sig,as)
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
      subroutine shearfac18(d,xl,ndm,cappa)
c-----------------------------------------------------------------------
c.... shear correction factor due to size
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  d(*),xl(ndm,4),yl(ndm)

      cappa = d(8)
       
      if(cappa.lt.0.d0) then
c....   Ae ... (longest element side)^2 
        xnu = d(2)   
        hs  = d(3)
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
        cappa  = dabs(d(8))/(1.d0+ae/(2.d0*hs*hs*(1.d0+xnu))) ! Tessler
      end if 
c      
      return
      end  
c
