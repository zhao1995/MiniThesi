      subroutine elmt13(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c-----------------------------------------------------------------------
c.....geometrical nonlinear shallow shell element 
c
c      1: E   = Young's modulus
c      2: nu  = Poisson's ratio
c      3: h   = thickness of shell
c      4: q   = transverse load
c      5: gam = spec. weight, M = gam*V, Gam = rho * g
c      6: lin = 0 linear,=1 nonlinear
c.... * lumped mass matrix [Lmas]
c..   * class.buckl. analysis: [ K_0 + K_U + l * K_G ] * phi = 0
c....   [Tang,cbuc], [Geom,cbuc], [Subs,,1]  (imtyp=ibuck=2)
c.... * lin.  buckl. analysis: [ K_0 + l * (K_U +K_G)] * phi = 0
c....   [Tang,lbuc], [Geom,lbuc], [Subs,,1]  (imtyp=ibuck=3)
c....   variable thickness if h = 0, input via temp, plot via stre,10
c.... * error analysis: active only bending terms
c-----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE eldata
      USE errin2
      USE errin3
      USE evdata
      USE iofile
      USE prlod
      USE prisdat
      USE qload
      USE strnam
      implicit double precision (a-h,o-z)
      dimension   xl(ndm,*),tl(*),ix(*),
     1 d(*),pgb(4),wgb(4),pgs(4),wgs(4),gp(2),cartd(2,9),
     2 deriv(2,9),s(nst,nst),shp(9),sig(8),
     3 dme(3,3),db(3,3),ds(2,2),bmi(3,2),bmj(3,2),bbi(3,3),bbj(3,3),
     4 bsi(2,3),bsj(2,3),bni(3,3),bnj(3,3),bvi(3,3),bvj(3,3),
     5 ul(ndf,*),gr(10),p(nst),e_ome13(numerr)
      dimension h1(*),h2(*),h3(*)
c
c.... go to correct array processor
      go to(1,2,3,3,5,3,2,3,3,2,3,2,2,2,2,2,2,2,2,2,2,22) isw
c.... input material properties
1     if(ior.lt.0) write(*,1001)
1001  format(' Input: E  nu  h  q gamma lin',$)
      call dinput(d,6)
      write(iow,1002) (d(i),i=1,6)
1002  format(5x,'materialdata:',/,
     + 5x,'elastic modulus............',g15.4,/,
     + 5x,'poissons ratio................',g12.4,/,
     + 5x,'thickness(h = 0 variable).....',g12.4,/,
     + 5x,'element load (transverse).....',g12.4,/,
     + 5x,'spec. weight (gamma)..........',g12.4,/,
     + 5x,'linear/nonlinear.....0/1......',g12.4)
      if(ior.lt.0) write(*,1002) (d(i),i=1,6)
c.... description of stresses  
      strsus( 1) = '  N-FORCE N_xx '
      strsus( 2) = '  N-FORCE N_yy '
      strsus( 3) = '  N-FORCE N_xy '
      strsus( 4) = '  MOMENT M_xx  '
      strsus( 5) = '  MOMENT M_yy  '
      strsus( 6) = '  MOMENT M_xy  '
      strsus( 7) = '  Q-FORCE Q_xz '
      strsus( 8) = '  Q-FORCE Q_yz '
      strsus( 9) = '  MOMENT M_1   '
      strsus(10) = '  MOMENT M_2   '
      strsus(11) = '  THICKNESS    '
      do is =12,25
        strsus(is) = '               '
      end do
c...  names for principal forces/moments
      nptyp = 5 
c...  position for principal forces/moments
      nprip(1)=1
      nprip(2)=3
      nprip(3)=2
      nprip(4)=4
      nprip(5)=6
      nprip(6)=5
      nprip(7)=7
      nprip(8)=8
c
2     return
c.... set number of printed/plotted stresses to 8
c.... a negative prevents computation of principle values.
3     istv = -11
      lin = d(6)
      ngb = 2
      if(nel.gt.4) ngb=3
      ngs = ngb-1
c      if(isw.eq.8.or.isw.eq.9) ngs=ngs+1
      call gausfs(ngb,pgb,wgb)
      call gausfs(ngs,pgs,wgs)
      if(isw.eq.4.or.isw.eq.8.or.isw.eq.9) goto 33
c.....Membrane and bending terms in K_T
      do 50 igb = 1,ngb
      do 50 jgb = 1,ngb
         xsi = pgb(igb)
         eta = pgb(jgb)
c.....   shape functions,area etc.
         call sfr2fs(deriv,eta,xsi,nel,shp)
         call jacofs(cartd,deriv,djacb,xl,n,nel,shp,ndm,iow)
         da = djacb*wgb(igb)*wgb(jgb)
c....    variable thickness
         thick = d(3)
         if(thick.eq.0.0d0) then
         do 30 i = 1,nel
            thick = thick + shp(i)*tl(i)
30       continue
         endif
c.....   Elasticity matrix
         call modfs(dme,db,ds,d,1,1,0,thick)
c.....derivation  zo,x , zo,y
         zkx = 0.0
         zky = 0.0
         do 27 ino = 1,nel
            zkx= zkx+ cartd(1,ino) * xl(3,ino)
   27       zky= zky+ cartd(2,ino) * xl(3,ino)
c.....displacement gradients
        call gradfs(cartd,gr,ul,ndf,nel)
c.....stresses (N and M)  and strains
      call strfs(cartd,shp,ul,gr,dme,db,ds,sig,nel,ndf,1,0,zkx,zky,lin)
c...... stiffness matrix for K_O for linear buckling analysis
      if(ibuck.eq.3) then
        do 320 ino =1,nel
        call bmatfs(bmi,bbi,bsi,cartd,ino,shp,nel,1,1,0)
        call bmafsnl(bvi,cartd,ino,zkx,zky,nel)
          do 420 jno = ino,nel
          call bmatfs(bmj,bbj,bsj,cartd,jno,shp,nel,1,1,0)
          call bmafsnl(bvj,cartd,jno,zkx,zky,nel)
c.....bbi T * db * bbj
          call subfs(bbi,bbj,da,db,s,ndf,nst,ino,jno,3,3,3,3,3)
c.....bmi T * dme * bm  (linear)
          call subfs(bmi,bmj,da,dme,s,ndf,nst,ino,jno,2,3,2,1,1)
c.....bmi T * dme * bvj
          call subfs(bmi,bvj,da,dme,s,ndf,nst,ino,jno,2,3,3,1,3)
c.....bvi T * dme * bmj
          call subfs(bvi,bmj,da,dme,s,ndf,nst,ino,jno,3,3,2,3,1)
c.....bvi T * dme * bvj
          call subfs(bvi,bvj,da,dme,s,ndf,nst,ino,jno,3,3,3,3,3)
  420     continue
  320    continue
c.... nonlinear stiffness matrix
      else
c.....add w + zo for bni
        wkx  =  zkx
        wky  =  zky
        if(lin.eq.1) then 
          wkx  =  zkx + gr(3)
          wky  =  zky + gr(8)
        end if
        do 32 ino =1,nel
        call bmatfs(bmi,bbi,bsi,cartd,ino,shp,nel,1,1,0)
        call bmafsnl(bni,cartd,ino,wkx,wky,nel)
c..... residual
        call mulfs(bni,sig(1),da,p,ndf,nst,ino,3,3,3)
        call mulfs(bmi,sig(1),da,p,ndf,nst,ino,2,3,1)
        call mulfs(bbi,sig(4),da,p,ndf,nst,ino,3,3,3)
c..... external load vector due to transversal constant load 
       call qload13(qz,d,aqloa,numel,n,mqloa,propq,prop,isw)
       iq = (ino-1)*ndf+3
       p(iq) = p(iq) + shp(ino)*qz*da
c
      if(isw.eq.6) goto 32
          do 42 jno = ino,nel
          call bmatfs(bmj,bbj,bsj,cartd,jno,shp,nel,1,1,0)
          call bmafsnl(bnj,cartd,jno,wkx,wky,nel)
c.....bbi T * db * bbj
          call subfs(bbi,bbj,da,db,s,ndf,nst,ino,jno,3,3,3,3,3)
c.....bmi T * dme * bm  (linear)
          call subfs(bmi,bmj,da,dme,s,ndf,nst,ino,jno,2,3,2,1,1)
c.....bmi T * dme * (bl+bv)
          call subfs(bmi,bnj,da,dme,s,ndf,nst,ino,jno,2,3,3,1,3)
c.....(bl+bv) trans * dme * bm
          call subfs(bni,bmj,da,dme,s,ndf,nst,ino,jno,3,3,2,3,1)
c.....(bl+bv) trans * dme * (bl+bv)
          call subfs(bni,bnj,da,dme,s,ndf,nst,ino,jno,3,3,3,3,3)
c.....initial stress matrix, only nonlinear case,ibuck=2-->K_0+K_U
          if(ibuck.eq.2) goto 42
          if(lin.eq.1) call ksifs(ino,jno,s,cartd,sig,nel,nst,ndf,da)
   42     continue
   32   continue
      endif
   50 continue
c.....Shear terms in K_T
   33 continue
c.... intial values for element errors
      e_ome13 = 0.0d0

       do 51 igs = 1,ngs
       do 51 jgs = 1,ngs
         xsi = pgs(igs)
         eta = pgs(jgs)
c.....   shape functions,area etc.
         call sfr2fs(deriv,eta,xsi,nel,shp)
         call jacofs(cartd,deriv,djacb,xl,n,nel,shp,ndm,iow)
         da = djacb*wgs(igs)*wgs(jgs)
c....    variable thickness
         thick = d(3)
         if(thick.eq.0.0d0) then
         do 54 i = 1,nel
            thick = thick + shp(i)*tl(i)
54       continue
         endif
c.....   Elasticity matrix
         call modfs(dme,db,ds,d,1,1,1,thick)
c.... derivation zo,x zo,y
       zkx = 0.0
       zky = 0.0
      do 53 ino = 1,nel
      zkx= zkx+ cartd(1,ino) * xl(3,ino)
   53 zky= zky+ cartd(2,ino) * xl(3,ino)
c.....displacement gradients
      call gradfs(cartd,gr,ul,ndf,nel)
c.....stresses (Q) and strains
      call strfs(cartd,shp,ul,gr,dme,db,ds,sig,nel,ndf,1,1,zkx,zky,lin)
      if(isw.eq.4) goto 46
      if(isw.eq.8) goto 82
      if(isw.eq.9) goto 90
c.....bsi T * ds * bsj
         do 31 ino = 1,nel
         call bmatfs(bmi,bbi,bsi,cartd,ino,shp,nel,0,0,1)
         call mulfs(bsi,sig(7),da,p,ndf,nst,ino,3,2,3)
         if(isw.eq.6) goto 31
           do 41 jno = ino,nel
           call bmatfs(bmj,bbj,bsj,cartd,jno,shp,nel,0,0,1)
   41      call subfs(bsi,bsj,da,ds,s,ndf,nst,ino,jno,3,2,3,3,3)
   31 continue
      goto 51
46    continue
c.....coordinates of Gauss points (for shear!!)
        do 44 idm = 1,2
        gp(idm) = 0.0
          do 44 ino = 1,nel
   44     gp(idm) = gp(idm) + xl(idm,ino) * shp(ino)
c.....Output stresses  (N M and Q)
      mct = mct - 1
      if(mct.gt.0) go to 430
        write(iow,2005)o,head
        if(ior.lt.0) write(*,2005)o,head
2005    format(a1,20a4,//,2x,'E L E M E N T   S T R E S S E S',//,
     1  2x,'EL',1x,'MAT',1x,'1-COORD',1x,'2-COORD',
     2  2X,'***NX***',3X,'***NY***',3X,'***NXY***',
     3  2X,'***MX***',3X,'***MY***',3X,'***MXY***',
     4  2X,'***QX***',3X,'***QY***',//)
        mct = 50
430     write(iow,2006)n,ma,(gp(i),i=1,2),(sig(i),i=1,8)
        if(ior.lt.0) write(*,2006)n,ma,(gp(i),i=1,2),(sig(i),i=1,8)
2006  format(1x,i3,i4,2f8.3,1x,8e11.4)
      goto 51
c.... stress plot N,M,Q
82    continue
c.... energy of stresses for error analysis
      call enerel13(d,sig,da)
c.... calculate nodal stresses
      if(iplma(ma).eq.0)  return ! only if MATN
      call plotfs(ix,strea,strea(1+numnp),shp,sig,nel,da,numnp,tl)
      goto 51
c.... energy of stress differences for error analysis
90    continue
      call sterr13(d,ix,strea(1+numnp),shp,sig,nel,da,numnp,e_ome13)
   51 continue
c.....plot/print/add errors
      if(isw.eq.9) then
        if (iet(1).ne.3) then
          nerr = 1  ! dummy for plot/print
        endif
        call elmterr(ix,xl,ndm,numel,e_ome13,e_ome)
      endif 
c.....lower part of K_T
       if(isw.eq.3.or.isw.eq.11) then
        do 52 i = 1,nst
           do 52 j = 1,i
   52      s(i,j) =s(j,i)
       endif
      return
c.... lumped mass or geometrical matrix
5     lin=d(6)
c       Number of Gauss points
        ngb = 2
        if (nel.gt.4) ngb = 3
c       Position and weight of Gauss points
        call gausfs(ngb,pgb,wgb)
      if(imtyp.eq.1) then
c....   lumped mass matrix
        gam  = d(5)
c....   8-node element produces negative mass terms
        if(nel.eq.8) stop 'element produces negative mass terms'
c....   loop over Gauss points
        do 501 igb = 1,ngb
           do 501 jgb = 1,ngb
c....   Coordinates Gauss points
           xsi = pgb(igb)
           eta = pgb(jgb)
c....   shape functions
           call sfr2fs(deriv,eta,xsi,nel,shp)
c....   Jacobian,Inverse and determinant
           call jacofs(cartd,deriv,djacb,xl,n,nel,shp,ndm,iow)
c....   area of element
           da = djacb*wgb(igb)*wgb(jgb)
c....   variable thickness
         thick = d(3)
         if(thick.eq.0.0d0) then
         do 55 i = 1,nel
            thick = thick + shp(i)*tl(i)
55       continue
         endif
c....   lumped mass matrix  
           do 520 ino = 1, nel
             iq   = (ino-1)*ndf
             dpma = shp(ino)*gam*da*thick
             dpmi = dpma*thick*thick/12.d0
             p(iq+1) = p(iq+1) + dpma
             p(iq+2) = p(iq+2) + dpma
             p(iq+3) = p(iq+3) + dpma
             p(iq+4) = p(iq+4) + dpmi
             p(iq+5) = p(iq+5) + dpmi
 520       continue
 501    continue
         if(n.eq.1) write(*,570)
 570    format(1x,'LMAS implemented/ CMAS not implemented for ELsme3_6')
        return
c....   Matrix K_U + K_G for  buckling analysis
      elseif(imtyp.eq.2.or.imtyp.eq.3) then
c....   loop over Gauss points
        do 550 igb = 1,ngb
          do 550 jgb = 1,ngb
c....      Coordinates Gauss points
           xsi = pgb(igb)
           eta = pgb(jgb)
c....      shape functions
           call sfr2fs(deriv,eta,xsi,nel,shp)
c....      Jacobian,Inverse and determinant
           call jacofs(cartd,deriv,djacb,xl,n,nel,shp,ndm,iow)
c....      area of element -K_sigma!
           da = -djacb*wgb(igb)*wgb(jgb)
c....      variable thickness
           thick = d(3)
           if(thick.eq.0.0d0) then
           do 551 i = 1,nel
             thick = thick + shp(i)*tl(i)
551        continue
           endif
c....      Elasticity matrix
           call modfs(dme,db,ds,d,1,1,0,thick)
c....      shallow shell parts  z,x , z,y
           zkx = 0.0
           zky = 0.0
           do 527 ino = 1,nel
              zkx= zkx+ cartd(1,ino) * xl(3,ino)
              zky= zky+ cartd(2,ino) * xl(3,ino)
  527      continue
c....      Displacement gradients
           call gradfs(cartd,gr,ul,ndf,nel)
c....   membrane and bending strains, normal forces and stress couples
        call strfs(cartd,shp,ul,gr,dme,db,ds,sig,nel,ndf,1,0,zkx,zky,
     +       lin)
c       Matrix K_sigma
           if (imtyp.eq.2) then
           do 532 ino =1,nel
            do 542 jno = ino,nel
c...........initial stress matrix
  542        if(lin.eq.1) call ksifs(ino,jno,s,cartd,sig,nel,nst,ndf,da)
  532      continue
           elseif(imtyp.eq.3) then
c....     K_U+K_G= int(B_m+b_v+b_nl)**T*D*(B_m+b_v+b_nl)dv - K_o + K_G
c
           wkx  =  gr(3)
           wky  =  gr(8)
c
           do 560 ino =1,nel
            call bmatfs(bmi,bbi,bsi,cartd,ino,shp,nel,1,0,0)
            call bmafsnl(bvi,cartd,ino,zkx,zky,nel)
            call bmafsnl(bni,cartd,ino,wkx,wky,nel)
            do 565 jno = ino,nel
             call bmatfs(bmj,bbj,bsj,cartd,jno,shp,nel,1,0,0)
             call bmafsnl(bvj,cartd,jno,zkx,zky,nel)
             call bmafsnl(bnj,cartd,jno,wkx,wky,nel)
c...........bmi T * dme * bnlj
             call subfs(bmi,bnj,da,dme,s,ndf,nst,ino,jno,2,3,3,1,3)
c...........bnliT * dme * bmj
             call subfs(bni,bmj,da,dme,s,ndf,nst,ino,jno,3,3,2,3,1)
c...........bvi T * dme * bnlj
             call subfs(bvi,bnj,da,dme,s,ndf,nst,ino,jno,3,3,3,3,3)
c...........bnliT * dme * bvj
             call subfs(bni,bvj,da,dme,s,ndf,nst,ino,jno,3,3,3,3,3)
c...........bnliT * dme * bnlj
             call subfs(bni,bnj,da,dme,s,ndf,nst,ino,jno,3,3,3,3,3)
c...........initial stress matrix
  565        call ksifs(ino,jno,s,cartd,sig,nel,nst,ndf,da)
  560      continue
           endif
  550   continue
c....   calculate lower part of K_Sigma + K_U
        do 522 i = 1,nst
          do 522 j = 1,i
  522        s(i,j) =s(j,i)
      endif
      return
c.... load vector
22    ngb = 2
      if(nel.gt.4) ngb=3
      ngs = ngb-1
      call gausfs(ngb,pgb,wgb)
      call gausfs(ngs,pgs,wgs)
      do 221 igb = 1,ngb
      do 221jgb = 1,ngb
         xsi = pgb(igb)
         eta = pgb(jgb)
c.....   shape functions,area etc.
         call sfr2fs(deriv,eta,xsi,nel,shp)
         call jacofs(cartd,deriv,djacb,xl,n,nel,shp,ndm,iow)
         da = djacb*wgb(igb)*wgb(jgb)
         do 222 ino =1,nel
c.....     external load vector due to transversal constant load 
           call qload13(qz,d,aqloa,numel,n,mqloa,propq,prop,isw)
           iq = (ino-1)*ndf+3
           p(iq) = p(iq) + shp(ino)*qz*da
222      continue 
221   continue
      return
      end
c----------------------------------------------------------------------
      subroutine bmatfs(bm,bb,bs,cartd,kk,shp,nel,ifm,ifb,ifs)
c***********************************************************************
c.....B-matrices for flat shell
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension  bm(3,2),bb(3,3),bs(2,3),cartd(2,nel),shp(nel)
      dnkdx = cartd(1,kk)
      dnkdy = cartd(2,kk)
c.....form bm
      if(ifm.eq.1) then
        do 1 irows = 1,3
        do 1 jcols = 1,2
    1   bm(irows,jcols) = 0.0
          bm(1,1) = dnkdx
          bm(2,2) = dnkdy
          bm(3,1) = dnkdy
          bm(3,2) = dnkdx
      endif
c.....form bb
      if(ifb.eq.1) then
        do 2 irows = 1,3
        do 2 jcols = 1,3
    2   bb(irows,jcols) = 0.0
          bb(1,3) =  dnkdx
          bb(2,2) = -dnkdy
          bb(3,2) = -dnkdx
          bb(3,3) =  dnkdy
      endif
c.....form bs
      if(ifs.eq.1) then
        do 3 irows = 1,2
        do 3 jcols = 1,3
    3   bs(irows,jcols) = 0.0
          bs(1,1) =  dnkdx
          bs(1,3) =  shp(kk)
          bs(2,1) =  dnkdy
          bs(2,2) = -shp(kk)
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine bmafsnl(bni,cartd,kk,wkx,wky,nel)
c***********************************************************************
c.....nonlinear B-matrix for flat shell
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension bni(3,3),cartd(2,nel)
      dnkdx = cartd(1,kk)
      dnkdy = cartd(2,kk)
        do 10 irows = 1,3
        do 10 jcols = 1,3
   10   bni(irows,jcols) = 0.0
           bni(1,1) = wkx * dnkdx
           bni(2,1) = wky * dnkdy
           bni(3,1) = wky * dnkdx + wkx * dnkdy
      return
      end
c----------------------------------------------------------------------
      subroutine ksifs(ino,jno,s,cartd,sig,nel,nst,ndf,da)
c***********************************************************************
c.... initial stress matrix K_s = g T * sig * g
c************************************************************************
      implicit double precision (a-h,o-z)
      dimension cartd(2,nel),s(nst,nst),sig(8)
      dndxi = cartd(1,ino)
      dndyi = cartd(2,ino)
      dndxj = cartd(1,jno)
      dndyj = cartd(2,jno)
      aksig = dndxi * sig(1) * dndxj + dndxi * sig(3) * dndyj +
     +        dndyi * sig(2) * dndyj + dndyi * sig(3) * dndxj
      aksig = aksig*da
c.....store in K_T at  3,3
      irow = (ino-1) * ndf + 3
      jcol = (jno-1) * ndf + 3
      s(irow,jcol) = s(irow,jcol) + aksig
      return
      end
c----------------------------------------------------------------------
      subroutine gradfs(cartd,gr,ul,ndf,nel)
c***********************************************************************
c.....displacement gradients
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension cartd(2,nel),gr(10),ul(ndf,nel)
        do 10 igr = 1,10
   10   gr(igr)   = 0.0
          do 20 ino = 1,nel
          dnidx = cartd(1,ino)
          dnidy = cartd(2,ino)
             do 20 idf = 1,5
             idfn = 5 + idf
             const = ul(idf,ino)
             gr(idf)  = gr(idf) + dnidx * const
   20        gr(idfn) = gr(idfn)+ dnidy * const
      return
      end
c----------------------------------------------------------------------
      subroutine strfs(cartd,shp,ul,gr,dme,db,ds,sig,nel,ndf,ifm,ifs,
     +  zkx,zky,lin)
c***********************************************************************
c.....Strains and stresses for flat shell
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension  cartd(2,nel),db(3,3),gr(10),ds(2,2),shp(nel),
     1           ul(ndf,nel),dme(3,3),em(3),eb(3),es(2),sig(8)
        do 10 isig = 1,8
   10   sig(isig)  = 0.0
      if(ifs.eq.1) then
c.....rotation at Gauss point
        xzrot = 0.0
        yzrot = 0.0
        nposn = 4
        nposn1= 5
          do 30 ino = 1,nel
          xzrot = xzrot+shp(ino)*ul(nposn,ino)
   30     yzrot = yzrot+shp(ino)*ul(nposn1,ino)
c.....shear strains and shear forces
        es(1)=gr(3)+yzrot
        es(2)=gr(8)-xzrot
        sig(7)=ds(1,1)*es(1)
        sig(8)=ds(2,2)*es(2)
      endif
      if(ifm.eq.1) then
c.....  normal strains and normal forces
        em(1) = gr(1)       +   gr(3)*zkx 
        em(2) = gr(7)       +   gr(8)*zky 
        em(3) = gr(2)+gr(6) +   gr(3)*zky + gr(8)*zkx 
        if(lin.eq.1) then
          em(1) = em(1) + gr(3)*gr(3) * 0.5d0
          em(2) = em(2) + gr(8)*gr(8) * 0.5d0
          em(3) = em(3) + gr(3)*gr(8)
        end if
        sig(1)= dme(1,1)*em(1)+dme(1,2)*em(2)
        sig(2)= dme(2,2)*em(2)+dme(2,1)*em(1)
        sig(3)= dme(3,3)*em(3)
c
c.....  bending strains and bending moments
        eb(1)= gr(5)
        eb(2)=-gr(9)
        eb(3)=-gr(4)+gr(10)
        sig(4) = db(1,1)*eb(1)+db(1,2)*eb(2)
        sig(5) = db(2,2)*eb(2)+db(2,1)*eb(1)
        sig(6) = db(3,3)*eb(3)
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine jacofs(cartd,deriv,djacb,xl,n,nel,shp,ndm,na)
c***********************************************************************
c.....Jacobi-Matrix and  cartesian derivatives of shape functions
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension xl(ndm,nel),
     +       cartd(2,nel),deriv(2,nel),shp(nel),xjaci(2,2),xjacm(2,2)
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
  600 format(/'Program Stop ELMT13- Jacofs',
     + 'zero or negative area in element',i5)
      return
      end
c----------------------------------------------------------------------
      subroutine mulfs(bimat,bjvec,da,p,ndf,nst,ino,ncoli,nrowj,irow)
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
      subroutine modfs(dme,db,ds,d,ifm,ifb,ifs,thick)
c***********************************************************************
c.....Elasticity matrix for flat shells
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension dme(3,3),db(3,3),ds(2,2),d(*)
      young = d(1)
      poiss = d(2)
c      thick = d(3)
      cappa = 5.0d0/6.0d0
      const = (young*thick)/(1.0d0-poiss*poiss)
c.....Membrane part
      if(ifm.eq.1) then
        do 1 irows = 1,3
        do 1 jcols = 1,3
    1   dme(irows,jcols) = 0.0
           dme(1,1) = const
           dme(2,2) = const
           dme(1,2) = const*poiss
           dme(2,1) = const*poiss
           dme(3,3) = const*(1.0d0-poiss)*0.5d0
      endif
c.....Bending part
      if(ifb.eq.1) then
      const = const*thick*thick/12.0d0
        do 2 irows = 1,3
        do 2 jcols = 1,3
    2   db(irows,jcols) = 0.0
          db(1,1) = const
          db(2,2) = const
          db(1,2) = const*poiss
          db(2,1) = const*poiss
          db(3,3) = const*(1.0d0-poiss)*0.5d0
      endif
c.....Shear part
      if(ifs.eq.1) then
      const = 0.5d0*young*thick*cappa/(1.0d0+poiss)
        do 3 irows = 1,2
        do 3 jcols = 1,2
    3   ds(irows,jcols) = 0.0
          ds(1,1) = const
          ds(2,2) = const
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine subfs(bi,bj,da,dmat,s,ndf,nst,ino,jno,
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
      subroutine sfr2fs (deriv,eta,xsi,nel,shp)
c***********************************************************************
c.....Shape functions and derivatives for linear,quadratic
c.....lagrangian and serendipity isoparametric  2-d elements
c.....version for  FEAP Numbering Owen/Hinton!!
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
      subroutine gausfs(ngs,pgp,wgp)
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
      subroutine plotfs(ix,dt,st,shp,sig,nel,da,numnp,tl)
c***********************************************************************
c.....Plot     nx(1 -> 1)   ny(2 -> 2)  nxy(3 -> 3)
c              mx(4 -> 4)   my(5 -> 5)  mxy(6 -> 6)
c              qx(7 -> 7)   qy(8 -> 8)
c              m1(* -> 9)   m2(* ->10)  c phi(* ->11)
c              h    -> 11
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension dt(numnp),st(numnp,*),stre(3),shp(9),sig(8),ix(*),tl(*)
        do 10  j = 1,nel
           xsji = da*shp(j)
           ii = abs(ix(j))
           if(ii.le.0) go to 10
             dt(ii) = dt(ii) + xsji
             st(ii,1) = st(ii,1) + sig(1)*xsji
             st(ii,2) = st(ii,2) + sig(2)*xsji
             st(ii,3) = st(ii,3) + sig(3)*xsji
             st(ii,4) = st(ii,4) + sig(4)*xsji
             st(ii,5) = st(ii,5) + sig(5)*xsji
             st(ii,6) = st(ii,6) + sig(6)*xsji
             st(ii,7) = st(ii,7) + sig(7)*xsji
             st(ii,8) = st(ii,8) + sig(8)*xsji
c.....  principal moments m_1 and m_2
           stre(1) = sig(4)
           stre(2) = sig(6)
           stre(3) = sig(5)
           call pstres(stre,hm1,hm2,phi)
           st(ii,9)  = st(ii,9)  + hm1*xsji
           st(ii,10) = st(ii,10) + hm2*xsji
c.....      variable thickness
            st(ii,11)  = st(ii,11)  + tl(j)*xsji
10      continue
      return
      end
c
      subroutine enerel13(d,sig,da)
c----------------------------------------------------------------------
c.....sum energy for flat shell element for error analysis
c----------------------------------------------------------------------
      USE errin1
      implicit double precision (a-h,o-z)
      dimension d(*),sig(8)
      e    = d(1)
      xnue = d(2)
      h    = d(3)
      chi  = 5./6.
      g    = e/(2.*(1.+xnue))

c.....Energy-Norm
c.... membrane energy
cww   u_om(1) = u_om(1) + 
cww  + da*1./(e*h)*(sig(1)*sig(1)-2.*xnue *sig(1)*sig(2)
cww  +     + da*sig(2)*sig(2)+2.*(1.+xnue)*sig(3)*sig(3))
c.... bending energy
      u_om(1) = u_om(1) +
     + da*12./(e*h*h*h)*(sig(4)*sig(4)-2.*xnue *sig(4)*sig(5)
     +           +da*sig(5)*sig(5)+2.*(1.+xnue)*sig(6)*sig(6))
c.... shear energy
cww   u_om(1) = u_om(1) + da/(chi*g*h)*(sig(7)*sig(7)+sig(8)*sig(8))
 
c.....L2-Norm
c.... membrane energy
cww   u_om(2) = u_om(2) + da*dot(sig(1),sig(1),3)
c.... bending energy
      u_om(2) = u_om(2) + da*12./(h*h)*dot(sig(4),sig(4),3)
c.... shear energy
cww   u_om(2) = u_om(2) + da*chi*(1.-xnue)*0.5*dot(sig(7),sig(7),2)

      return
      end
c
      subroutine sterr13(d,ix,st,shp,sig,nel,da,numnp,e_ome13)
c----------------------------------------------------------------------
c.....energy of stress differences for flat shell element
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),ix(*),st(numnp,*),shp(9),sig(8),sigp(8),dsig(8),
     +          e_ome13(2)
c.... calculate projected stresses at gauss point
      call pzero(sigp,8)
      do 10 i = 1,nel
        ll = iabs(ix(i))
        if(ll.ne.0) then
          do 15  j = 1,8
            sigp(j) = sigp(j) + shp(i)*st(ll,j)
15        continue
        endif
10    continue

      do i = 1,8
        dsig(i) = sigp(i)-sig(i)
      enddo
c.... update element energy
      e    = d(1)
      xnue = d(2)
      h    = d(3)
      chi  = 5./6.
      g    = e/(2.*(1.+xnue))

c.... Energy-Norm
c.... Membran
cww   e_ome13(1) = e_ome13(1) 
cww  +  + da/(e*h)*(dsig(1)*dsig(1)-2.*    xnue *dsig(1)*dsig(2)
cww  +             +dsig(2)*dsig(2)+2.*(1.+xnue)*dsig(3)*dsig(3))
c.... Bending
      e_ome13(1) = e_ome13(1) 
     + + da*12./(e*h*h*h)*(dsig(4)*dsig(4)-2.*    xnue *dsig(4)*dsig(5)
     +                    +dsig(5)*dsig(5)+2.*(1.+xnue)*dsig(6)*dsig(6))
c.... Shear
cww   e_ome13(1) = e_ome13(1) 
cww  + + da/(chi*g*h)*(dsig(7)*dsig(7)+dsig(8)*dsig(8))

c.... L2-Norm
c.... Membran
cww   e_ome13(2)=e_ome13(2)+da*dot(dsig(1),dsig(1),3)
c.... Bending
      e_ome13(2)=e_ome13(2)+da*12./(h*h)*dot(dsig(4),dsig(4),3)
c.... Shear
cww   e_ome13(2)=e_ome13(2)+da*chi*(1.-xnue)*0.5*dot(dsig(7),dsig(7),2)

      return
      end
c
      subroutine qload13(qz,d,q,numel,n,mqloa,propq,prop,isw)
c----------------------------------------------------------
c.... set loads from macro qloa/mate
c----------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension q(numel,10),d(*)
      qz = 0.d0 
      if(isw.eq.22) then
        if(mqloa.ne.1) qz = q(n,1)*propq 
      else
        qz = d(4)*prop 
      end if 
      return
      end
