      subroutine elmt22(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c----------------------------------------------------------------------
c.... plane stress element 4/8/9 nodes                          ELcfp2_8
c.....geometrically nonlinearc,elasto-plastic with lin. isotr. hardening
c     (c) w.wagner aug 93
c
c----------------------------------------------------------------------
c..... Remarks:
c      Dicke nicht notwendig, geht nur ein in Volumenintegration->dv
c      Aenderungen bzw. offen:
c      in ylds22: - >>> Abbruch fuer f = 0
c                       Abbruch fuer ddlam < tol
c                          --> keine quadr. Konvergenz
c                       Abbruch fuer ddlam < tol*dlam(RLT)
c                          --> keine Loesung fuer sig in jeweils 1.Iteration
c                 - >>> Zuwachs  ddlam von RLT, wichtig bei starken Aend.
c                 -     Aenderung 1.e-8 fuer dlam um Sig auf F zu setzen.
c                 - Anzahl der lokalen Iterationen  50
c     27/08/1993 WW
c     23/06/2009 WW Update
c     04/11/2013 WW Verwendung für FTW2,Einbau Cauchy/1.PK, QLOA
c----------------------------------------------------------------------
c.... Material parameters  and definition of array d
c
c   1     E   = Young's modulus
c   2     nu  = Poisson's ratio
c   3     h   = thickness
c   4     at  = alpha_t
c   5     theta= const. temperature (load)
c   6     lin = linear = 0  nonlinear = 1   def=0
c   7     initial yield stress Y_0          def=100*E
c   8     cp  = Hardening parameter         def=0
c   9     
c-------------------------------------
c  10    nh number of history terms at GP: 9=sig(3),ep(3),epq,yld,dlam0
c  11    e/(1-nu*nu)
c  12    e/(1-nu*nu) * nu
c  13    G
c----------------------------------------------------------------------
c.... definition of array h  for each gauss point
c   1     sig_x
c   2     sig_y
c   3     sig_xy
c   4     ep_x
c   5     ep_y
c   6     ep_xy
c   7     epq      = effective plastic strain
c   8     sigv/y_0
c   9     dlam0    = plastic multiplier of last iterate
c-----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE eldata
      USE hdata
      USE iofile
      USE qload      
      USE strnam
      implicit double precision(a-h,o-z)
      dimension xl(ndm,*),tl(*),ix(*),d(*),pg(4),wg(4),gp(2),shp(3,9),
     1          s(nst,nst),p(nst),dmat(3,3),bm(3,2),btd(2,3),
     2          ul(ndf,*),vl(2,9),epsi(3),epsn(3),gradi(4),gradn(4),
     3          deps(3),sig(3),epn(3),psig(4),ql(2)
      dimension h1(*),h2(*),h3(*)

      save ngaus,pg,wg

c.... go to correct array processor
      go to(1,2,3,4,2,3,2,4,2,2,3,2,2,2,2,2,2,2,2,2,2,22), isw

c.... input material properties
1     if(ior.lt.0) write(*,1002)
1002  format(' Input: E, nu, thick alpha_t theta linear(0/1) Y-0 cp'/)
      call dinput(d,8)
c...  default values      
      if(d(7).eq.0.d0) d(7)=d(1)*100.d0  
c
      write(iow,1001) (d(i),i=1,8)
      if(ior.lt.0) write(*,1001) (d(i),i=1,8)
1001  format(5x,'materialdata:',/,5x,'elastic modulus...',f10.3,
     1/,5x,'poissons ratio....',f10.3,/,5x,'thickness.........',f10.3,
     2/,5x,'alpha_t...........',f10.3,/,5x,'temperature.......',f10.3,
     3/,5x,'lin=0,nonlin=1....',f10.0,/,5x,'in.Yield stress Y0',f10.3,
     4/,5x,'hardening parameter',f10.3)
c.... define number of points for h-array, 9=sig(3),ep(3),epq,yld,dlam0
      ngaus=2
      if(nen.gt.4) ngaus=3
      call gauss22(ngaus,pg,wg)
      d(10) = 9
      nh1 = d(10)*ngaus*ngaus
c.... Elastizitaetsmatrix
      d(11)  = d(1)/(1.0d0-d(2)*d(2))
      d(12) = d(11)*d(2)
      d(13) = d(11)*(1.0d0-d(2))*0.5d0

c.... Spannungsnamen
      strsus( 1) = '  Stress S_xx  '  
      strsus( 2) = '  Stress S_yy  '
      strsus( 3) = '  Stress S_xy  '
      strsus( 4) = '  Strain Ep_xx '
      strsus( 5) = '  Strain Ep_yy '
      strsus( 6) = '  Strain Ep_xy '
      strsus( 7) = '  Strain Ep_q  '
      strsus( 8) = '  Yield  Y/Y_0 '
      strsus( 9) = '  DLambda_0_   '
      strsus(10) = '  Stress P_xx  '
      strsus(11) = '  Stress P_yy  '
      strsus(12) = '  Stress P_xy  '
      strsus(13) = '  Stress P_yx  '
      do is =14,25     
        strsus(is) = '               '
      end do          

2     return

3     lin = d(6)
      nh =  d(10)

c.... displacements at t_n  = u_i - du_i
      do 300 i = 1,nel
         do 300 j = 1,2
            vl(j,i) = ul(j,i) - ul(j,i+nen)
300   continue
c.....loop over gauss points
      nn   = 1
      do 30 igaus = 1,ngaus
       do 30 jgaus = 1,ngaus
        xsi = pg(igaus)
        eta = pg(jgaus)
c.....  shape functions and derivatives
        call sfr222(shp,eta,xsi,nel)
        call jaco22(shp,djacb,xl,n,nel,ndm,iow)
        dv = djacb*wg(igaus)*wg(jgaus)*d(3)
c.....  displacement gradients at  t_i
        call grad22(shp,gradi,ul,ndf,nel)
c.....  total strains at t_i
        call eps22(gradi,epsi,lin,d)
c.....  displacement gradients at  t_n
        call grad22(shp,gradn,vl,ndf,nel)
c.....  total strains at t_n
        call eps22(gradn,epsn,lin,d)
c.....  incremental strains at t_i
        deps(1) = epsi(1)-epsn(1)
        deps(2) = epsi(2)-epsn(2)
        deps(3) = epsi(3)-epsn(3)
c.....  stresses at gauss point
        call str22(h1(nn),h2(nn),nh,deps,dmat,d)
        nn2 = nn-1
        do 301 i = 1,3
          sig(i) = h2(nn2+i)
301     continue
c....   nonlinear stiffness matrix and residual
        i1=0
        do 31 ii=1,nel
c....    B- matrix
          call bmat22(bm,shp,ii,gradi,nel,lin)
c....     residual G = P - Bt*S and matrix Bt*D
          do 32 i = 1,ndf
             do  33 k = 1,3
             btd(i,k) = 0.0
             p(i1+i) = p(i1+i) - bm(k,i)*sig(k)*dv
                 do  34 j = 1,3
34               btd(i,k) = btd(i,k)+bm(j,i)*dmat(j,k)
33           continue
32        continue
        if(isw.eq.6) goto 305
c....   tangent stiffness matrix
        j1 = i1
        do 35 jj = ii,nel
          call bmat22(bm,shp,jj,gradi,nel,lin)
          do 36  i = 1,ndf
            do 37  j = 1,ndf
              do 38  k = 1,3
38              s(i1+i,j1+j) = s(i1+i,j1+j) + btd(i,k)*bm(k,j)*dv
37          continue
36        continue
c.....    initial stress matrix
          if(lin.ne.0) then
            f = shp(1,ii)*sig(1)*shp(1,jj) + shp(1,ii)*sig(3)*shp(2,jj)
     1        + shp(2,ii)*sig(2)*shp(2,jj) + shp(2,ii)*sig(3)*shp(1,jj)
            s(i1+1,j1+1) = s(i1+1,j1+1) + f*dv
            s(i1+2,j1+2) = s(i1+2,j1+2) + f*dv
          end if
          j1 = j1 + ndf
35      continue
305     i1 = i1 + ndf
31     continue
       nn = nn + nh
   30 continue
      if(isw.ne.3) return
c.....lower part of stiffness matrix
      do 39 i = 1,nst
       do 39 j = 1,i
   39   s(i,j) =s(j,i)

      return

c.... output stresses
4     lin = d(6)
      nh =  d(10)
c.... displacements at t_n  = u_i - du_i
      do 400 i = 1,nel
        do 400 j = 1,2
          vl(j,i) = ul(j,i) - ul(j,i+nen)
400   continue
c.....loop over gauss points
      nn   = 1
      do 40 igaus = 1,ngaus
       do 40 jgaus = 1,ngaus
        xsi = pg(igaus)
        eta = pg(jgaus)
c.....  shape functions and derivatives
        call sfr222(shp,eta,xsi,nel)
        call jaco22(shp,djacb,xl,n,nel,ndm,iow)
        dv = djacb*wg(igaus)*wg(jgaus)*d(3)
c.....  displacement gradients at  t_i
        call grad22(shp,gradi,ul,ndf,nel)
c.....  total strains at t_i
        call eps22(gradi,epsi,lin,d)
c.....  displacement gradients at  t_n
        call grad22(shp,gradn,vl,ndf,nel)
c.....  total strains at t_n
        call eps22(gradn,epsn,lin,d)
c.....  incremental strains at t_i
        deps(1) = epsi(1)-epsn(1)
        deps(2) = epsi(2)-epsn(2)
        deps(3) = epsi(3)-epsn(3)
c.....  stresses at gauss point
        call str22(h1(nn),h2(nn),nh,deps,dmat,d)
        nn2 = nn-1
        do 401 i = 1,3
          sig(i) = h2(nn2+i)
401     continue
        do 402 i = 1,3
          epn(i) = h2(nn2+i+3)
402     continue
        epq   = h2(nn2+7)
        yt    = h2(nn2+8)
        dlam0 = h2(nn2+9)
c....   1.PK stresse 
        call pstr22(gradi,sig,psig,lin)

        if(isw.eq.4) then
          do 403 idm = 1,2
            gp(idm) = 0.0d0
            do 403 inode = 1,nel
  403         gp(idm) = gp(idm) + xl(idm,inode)* shp(3,inode)
          mct = mct - 1
          if(mct.gt.0) go to 41
          write(iow,4000)o,head
          if(ior.lt.0) write(*,4000)o,head
4000   format(a1,20a4,//,2x,'E L E M E N T   S T R E S S E S ',/,
     1       1x,'EL',1x,'M',1x,'1-COR',1x,'2-COR',
     2       2X,'***SX***',4X,'***SY***',4X,'**SXY***',4X,'***PX***',
     3       4X,'***PY***',4X,'**PXY***',4X,'**PXY***',
     4       4X,'*Eps_pl*',4x,'**Y/Y0**')
          mct = 50
41               write(iow,4001)n,ma,(gp(i),i=1,2),(sig(i),i=1,3),
     +                          (psig(i),i=1,4),epq,yt
       if(ior.lt.0)write(*,4001)n,ma,(gp(i),i=1,2),(sig(i),i=1,3),
     +                          (psig(i),i=1,4),epq,yt
4001   format(1x,i2,i2,2f6.3,9e12.5)

        else if(isw.eq.8) then
c....     plot stresses
          istv = -9
          call spl22(ix,strea,strea(1+numnp),sig,psig,epn,epq,yt,dlam0,
     1               shp,nel,numnp,dv)
       end if
       nn = nn + nh
   40 continue
      return

c.... QLOA
22    call qload22(aqloa,ql,numel,n,mqloa,propq)
c.....loop over gauss points
      do igaus = 1,ngaus
        do jgaus = 1,ngaus
          xsi = pg(igaus)
          eta = pg(jgaus)
c.....    shape functions and derivatives
          call sfr222(shp,eta,xsi,nel)
          call jaco22(shp,djacb,xl,n,nel,ndm,iow)
          dv = djacb*wg(igaus)*wg(jgaus)*d(3)
c....     load vector
          i1=0
          do ii=1,nel
            p(i1+1) = p(i1+1) + ql(1)*shp(3,ii)*dv    
            p(i1+2) = p(i1+2) + ql(2)*shp(3,ii)*dv    
            i1 = i1 + ndf
          end do ! ii
        end do ! jgaus
      end do ! igaus
      return
      end
c
      subroutine bmat22(bm,shp,k,dgrad,nel,lin)
c-----------------------------------------------------------------------
c.....B-Matrix for geometrical nonlinear plain stress element
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension bm(3,2),shp(3,nel),dgrad(4)
      dndx = shp(1,k)
      dndy = shp(2,k)
      uxkx=dgrad(1)
      uykx=dgrad(2)
      uxky=dgrad(3)
      uyky=dgrad(4)
      do 10 i = 1,3
      do 10 j = 1,2
10    bm(i,j) = 0.0
c.... linear part
      bm(1,1) = dndx
      bm(2,2) = dndy
      bm(3,1) = dndy
      bm(3,2) = dndx
      if(lin.eq.0) return
c.... nonlinear part
      bm(1,1) = bm(1,1) + uxkx*dndx
      bm(1,2) = uykx*dndx
      bm(2,2) = bm(2,2) + uyky*dndy
      bm(2,1) = uxky*dndy
      bm(3,1) = bm(3,1) + uxkx*dndy + uxky*dndx
      bm(3,2) = bm(3,2) + uykx*dndy + uyky*dndx
      return
      end
c
      subroutine grad22(shp,dgrad,ul,ndf,nel)
c-----------------------------------------------------------------------
c.....displacement gradients
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension shp(3,nel),dgrad(4),ul(ndf,*)
      do 10 i = 1,4
   10 dgrad(i) = 0.0
      do 20 i = 1,nel
      dndx = shp(1,i)
      dndy = shp(2,i)
      do 20 idf = 1,ndf
      ip = ndf + idf
      dgrad(idf) = dgrad(idf)+dndx*ul(idf,i)
   20 dgrad(ip) =  dgrad(ip) +dndy*ul(idf,i)
      return
      end
c
      subroutine eps22(dgrad,eplan,lin,d)
c-----------------------------------------------------------------------
c.....strains for geometrically nonlinear plain stress element
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension dgrad(4),eplan(3),d(*)
      uxkx=dgrad(1)
      uykx=dgrad(2)
      uxky=dgrad(3)
      uyky=dgrad(4)
c.....strains (linear)
      eplan(1) = uxkx
      eplan(2) = uyky
      eplan(3) = uxky + uykx
c.....strains (thermal loading)
      etheta = d(4)*d(5)
      eplan(1) = eplan(1) - etheta
      eplan(2) = eplan(2) - etheta
c
      if(lin.eq.0) goto 20
c.....strains (nonlinear)
      eplan(1) = eplan(1) + 0.5d0 * ( uxkx*uxkx + uykx*uykx )
      eplan(2) = eplan(2) + 0.5d0 * ( uxky*uxky + uyky*uyky )
      eplan(3) = eplan(3) +         ( uxkx*uxky + uykx*uyky )
20    return
      end
c
      subroutine jaco22(shp,djacb,xl,n,nel,ndm,iow)
c-----------------------------------------------------------------------
c.....Jacobi-Matrix and  cartesian derivatives of shape functions
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension xl(ndm,nel),
     1       shp(3,nel),xjaci(2,2),xjacm(2,2),cartd(2,9)
c.....Jacobi-matrix xjacm
      do 4 idm = 1,2
      do 4 jdm = 1,2
      xjacm(idm,jdm) = 0.0
      do 4 inode = 1,nel
4     xjacm(idm,jdm)=xjacm(idm,jdm)+shp(idm,inode)*xl(jdm,inode)
c.....Determinant and Inverse of Jacobian
      djacb=xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1)
      if(djacb) 6,6,8
    6 write(iow,600) n
      write(  *,600) n
      stop
    8 continue
      xjaci(1,1)= xjacm(2,2)/djacb
      xjaci(2,2)= xjacm(1,1)/djacb
      xjaci(1,2)=-xjacm(1,2)/djacb
      xjaci(2,1)=-xjacm(2,1)/djacb
c.....cartesian derivatives
      do 10 idm = 1,2
      do 10 i = 1,nel
      cartd(idm,i) = 0.0
      do 10 jdm = 1,2
10    cartd(idm,i)=cartd(idm,i)+xjaci(idm,jdm)*shp(jdm,i)
      do 20 idm=1,2
      do 20 i=1,nel
20    shp(idm,i) = cartd(idm,i)
  600 format(1x,'program stop in subroutine jaco22',/,1x,
     +'zero or negative area for element number',i5)
      return
      end
c
      subroutine gauss22(ngaus,pg,wg)
c-----------------------------------------------------------------------
c.....Gauss points  ngaus = 1,3
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension pg(4),wg(4)
      do 10 i = 1,4
      pg(i) = 0.0
   10 wg(i) = 0.0
      if(ngaus.gt.1) goto 2
      pg(1) = 0.0
      wg(1) = 2.0d0
      return
    2 continue
      if(ngaus.gt.2) goto 4
      pg(1) = -dsqrt(3.0d0)/3.0d0
      wg(1) =  1.0d0
      goto 6
    4 pg(1) = -dsqrt(0.6d0)
      pg(2) =  0.0d0
      wg(1) = 5.0d0/9.0d0
      wg(2) = 8.0d0/9.0d0
    6 kgaus = ngaus/2
      do 8 ig = 1,kgaus
      jg = ngaus + 1 - ig
      pg(jg) = -pg(ig)
    8 wg(jg) =  wg(ig)
      return
      end
c
      subroutine sfr222 (shp,eta,xsi,nel)
c-----------------------------------------------------------------------
c.....shape functions and their derivatives for linear,quadratic
c.....lagrangian and serendipity isoparametric  2-d elements
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension shp(3,nel)
      deins=1.0d0
      dzwei=2.0d0
      halbe=deins/dzwei
      viert=halbe/dzwei
      s = xsi
      t = eta
      if(nel.gt.4) goto 10
      st = s*t
c.....shape functions for 4-node element
      shp(3,1) = (deins-t-s+st)*viert
      shp(3,2) = (deins-t+s-st)*viert
      shp(3,3) = (deins+t+s+st)*viert
      shp(3,4) = (deins+t-s-st)*viert
c.....derivatives
      shp(1,1) = (-deins+t)*viert
      shp(1,2) = (+deins-t)*viert
      shp(1,3) = (+deins+t)*viert
      shp(1,4) = (-deins-t)*viert
      shp(2,1) = (-deins+s)*viert
      shp(2,2) = (-deins-s)*viert
      shp(2,3) = (+deins+s)*viert
      shp(2,4) = (+deins-s)*viert
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
c.....shape functions for 8-node element
      shp(3,1) = (-deins+st+ss+tt-sst-stt)*viert
      shp(3,2) = (-deins-st+ss+tt-sst+stt)*viert
      shp(3,3) = (-deins+st+ss+tt+sst+stt)*viert
      shp(3,4) = (-deins-st+ss+tt+sst-stt)*viert
      shp(3,5) = (+deins-t-ss+sst)*halbe
      shp(3,6) = (+deins+s-tt-stt)*halbe
      shp(3,7) = (+deins+t-ss-sst)*halbe
      shp(3,8) = (+deins-s-tt+stt)*halbe
c.....derivatives
      shp(1,1) = (+t+s2-st2-tt)*viert
      shp(1,2) = (-t+s2-st2+tt)*viert
      shp(1,3) = (+t+s2+st2+tt)*viert
      shp(1,4) = (-t+s2+st2-tt)*viert
      shp(2,1) = (+s+t2-ss-st2)*viert
      shp(2,2) = (-s+t2-ss+st2)*viert
      shp(2,3) = (+s+t2+ss+st2)*viert
      shp(2,4) = (-s+t2+ss-st2)*viert
      shp(1,5) =  -s+st
      shp(1,7) =  -s-st
      shp(2,6) =  -t-st
      shp(2,8) =  -t+st
      shp(1,6) = (+deins-tt)*halbe
      shp(1,8) = (-deins+tt)*halbe
      shp(2,5) = (-deins+ss)*halbe
      shp(2,7) = (+deins-ss)*halbe
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
c.....shape functions for 9 noded element
      shp(3,1)=viert*s9*st*t9
      shp(3,5)=halbe*(deins-ss)*t*t9
      shp(3,2)=viert*s1*st*t9
      shp(3,6)=halbe*s*s1*(deins-tt)
      shp(3,3)=viert*s1*st*t1
      shp(3,7)=halbe*(deins-ss)*t*t1
      shp(3,4)=viert*s9*st*t1
      shp(3,8)=halbe*s*s9*(deins-tt)
      shp(3,9)=(deins-ss)*(deins-tt)
c.....derivatives
      shp(1,1)=viert*t*t9*(-deins+s2)
      shp(1,5)=-st*t9
      shp(1,2)=viert*(deins+s2)*t*t9
      shp(1,6)=halbe*(deins+s2)*(deins-tt)
      shp(1,3)=viert*(deins+s2)*t*t1
      shp(1,7)=-st*t1
      shp(1,4)=viert*(-deins+s2)*t*t1
      shp(1,8)=halbe*(-deins+s2)*(deins-tt)
      shp(1,9)=-s2*(deins-tt)
      shp(2,1)=viert*(-deins+t2)*s*s9
      shp(2,5)=halbe*(deins-ss)*(-deins+t2)
      shp(2,2)=viert*s*s1*(-deins+t2)
      shp(2,6)=-st*s1
      shp(2,3)=viert*s*s1*(deins+t2)
      shp(2,7)=halbe*(1-ss)*(deins+t2)
      shp(2,4)=viert*s*s9*(deins+t2)
      shp(2,8)=-st*s9
      shp(2,9)=-t2*(deins-ss)
      return
      end
c
      subroutine spl22(ix,dt,st,sig,psig,epn,epq,yt,dlam0,shp,nel,
     +                 numnp,dv)
c-----------------------------------------------------------------------
c.....plot stresses for nonlinear plain stress element
c      1 = S_x   2 = S_y   3 = S_xy  4 = Ep_x  5 = Ep_y  6 = Ep_xy
c      7 = E_pq  8 = y/y0  9 = dlam0
c     10 = P_x  11 = P_y  12 = P_xy 13 = P_yx
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension dt(numnp),st(numnp,*),ix(*),shp(3,9),sig(3),epn(3),
     +          psig(4)
      do 10 i = 1,nel
        xsji = dv*shp(3,i)
        ii = abs(ix(i))
        if(ii.eq.0) go to 10
        dt(ii) = dt(ii) + xsji
        st(ii, 1) = st(ii, 1) +  sig(1)*xsji
        st(ii, 2) = st(ii, 2) +  sig(2)*xsji
        st(ii, 3) = st(ii, 3) +  sig(3)*xsji
        st(ii, 4) = st(ii, 4) +  epn(1)*xsji
        st(ii, 5) = st(ii, 5) +  epn(2)*xsji
        st(ii, 6) = st(ii, 6) +  epn(3)*xsji
        st(ii, 7) = st(ii, 7) +     epq*xsji
        st(ii, 8) = st(ii, 8) +      yt*xsji
        st(ii, 9) = st(ii, 9) +   dlam0*xsji
        st(ii,10) = st(ii,10) + psig(1)*xsji
        st(ii,11) = st(ii,11) + psig(2)*xsji
        st(ii,12) = st(ii,12) + psig(3)*xsji
        st(ii,13) = st(ii,13) + psig(4)*xsji
10    continue
      return
      end
c
      subroutine str22(h1,h2,nh,deps,dmat,d)
c-----------------------------------------------------------------------
c.... stresses and tangent matrix for plane stress J-2 plasticity with
c     linear isotropic hardening
c
c.... 1. Input parameters
c
c        d      - array of material constants
c        deps   - incremental strains
c        epn    - plastic strains at t-n
c        epq    - effective plastic strain at t-n
c        sig    - stresses at t-n (in h1)
c        dlam0  - plastic multiplier at time t-n
c
c.... 2. Output parameters
c
c        epn    - plastic strains at t-n+1
c        epq    - effective platic strain at t-n+1
c        sig    - stresses at t-n+1 (in h2)
c        dmat   - "tangent" matrix at t-n+1
c        dlam0  - plastic multiplier at time t-n+1
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),h1(nh),h2(nh),deps(3),epn(3),sigtr(3),sig(3),
     1          dmat(3,3)
      common /eldata/ dm,n,ma,mct,iel,nel
      save /eldata/
      tol = 1.d-8
c.... move time n data to local array
c.... stresses at time n
      do 20 i = 1,3
  20  sig(i)  = h1(i)
c.... plastic strains at time t-n
      do 21 i = 1,3
  21  epn(i)  = h1(i+3)
c.... effective plastic strain
      epq = h1(7)
c.... plastic mutltiplier
      dlam0 = h1(9)
c.... trial stresses
      sigtr(1) = sig(1) + d(11)*deps(1) + d(12)*deps(2)
      sigtr(2) = sig(2) + d(12)*deps(1) + d(11)*deps(2)
      sigtr(3) = sig(3) + d(13)*deps(3)
c.... set up elastic tangent
      call pzero(dmat,9)
      dmat(1,1) = d(11)
      dmat(2,2) = d(11)
      dmat(1,2) = d(12)
      dmat(2,1) = d(12)
      dmat(3,3) = d(13)
c.... compute the yield state for the trial stresses
      sps=2.*((sigtr(1)**2-sigtr(1)*sigtr(2)+sigtr(2)**2)/3.
     1         +sigtr(3)**2)
      sigv = dsqrt(1.5*sps)
      yld  = d(7) + d(8)*epq
c.... for deps = 0 iteration is not necessary --> equilibrium!!
      enorm   = deps(1)*deps(1)+deps(2)*deps(2)+deps(3)*deps(3)
      if(enorm.lt.tol) dlam0 = 0.0
c
      if(sigv.gt.yld) then
c....   plastic correction necesarry
        call ylds22(d,sigtr,dmat,epn,epq,sigv,dlam0)
      end if
c.... set stresses etc. and  save time n+1 data from local array
      do 50 i = 1,3
50      h2(i) = sigtr(i)
      do 51 i = 1,3
51      h2(i+3) = epn(i)
      h2(7)   = epq
      h2(8)   = sigv/d(7)
      h2(9)   = dlam0
      return
      end
c
      subroutine ylds22(d,sigtr,dmat,epn,epq,sigv,dlam0)
c-----------------------------------------------------------------------
c.... plane stress plasticity routine for return map algorithm
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),sigtr(3),sig(3),dmat(3,3),ps(3),psh(3),epn(3)
      common /iofile/ ior,iow
c.... set parameters
      tol  = 1.d-8
      one3 = 1./3.
      two3 = 2./3.
c.... use plastic multiplier of last iterate as start value
      dlam = dlam0
      icnt = 0
100   continue
c.... calulate elements of matrix H**-1
      x4=1./d(1)+ two3*dlam
      x5=-d(2)/d(1)- one3*dlam
      x6=x4*x4-x5*x5
      x1=x4/x6
      x2=-x5/x6
      x3=d(13)/(1.+2.*d(13)*dlam)
c.... calulate stresses sig = h**-1*d**-1*sigtr
      s1 = (sigtr(1)-d(2)*sigtr(2))/d(1)
      s2 = (sigtr(2)-d(2)*sigtr(1))/d(1)
      s3 = sigtr(3)/d(13)
      sig(1) = x1*s1 + x2*s2
      sig(2) = x2*s1 + x1*s2
      sig(3) = x3*s3
c.... compute sig*P*sig
      sps = 2.*((sig(1)**2-sig(1)*sig(2)+sig(2)**2)/3.+sig(3)**2)
      sigv=dsqrt(1.5*sps)
c.... compute depq = sqrt(2/3*sig*P*sig)
      depq= dsqrt(two3*sps)
c.... compute effective plastic strain  dep = dlam*depq
      epg = epq + dlam*depq
c.... compute actual yield stress
      yld = d(7) + d(8)*epg
c.... compute yield criterion  f = 0.5*sps - 1/3 yld**2
      f = 0.5*sps - one3*yld**2
c.... basic values for derivative of yield criterion
c      c2 = depq = dsqrt(2/3*s*ps), c3 = ps * h**-1*ps
c.... calculate stresses ps,psh
      s1q = two3*sig(1)-one3*sig(2)
      s2q = two3*sig(2)-one3*sig(1)
      s3q = 2.*sig(3)
      ps(1)  = s1q
      ps(2)  = s2q
      ps(3)  = s3q
      psh(1) = x1*s1q + x2*s2q
      psh(2) = x2*s1q + x1*s2q
      psh(3) = x3*s3q
      c2 = depq
      c3 = dot(ps,psh,3)
c.... compute derivative of yield criterion
      fdl = -(1.-two3*two3*yld*d(8)*dlam/c2)*c3 -  two3*yld*d(8)*c2
c.... compute increment of dlam
      ddlam = -f/fdl
c.... update dlam
      if(ddlam.gt.0.0d0.or.abs(ddlam).lt.abs(dlam)) then
         dlam = dlam + ddlam
      else
         dlam = 0.5*abs(dlam)
      end if
c      dlam = dlam + ddlam
      icnt = icnt + 1
      if(icnt.gt.50) go to 110
c.... test for f
      if(abs(f).gt.tol) go to 100
c.... test for ddlam
c      if(abs(ddlam).gt.tol) go to 100
c.... test for ddlam  by RLT
c      if(abs(ddlam).gt.tol*abs(dlam)) go to 100
      go to 120
110   write(iow,2000)
2000  format(' * * Warning * * failure to converge in ylds22')
120   continue
c.... the term tol = 1.d-8 is to ensure that the resulting stress is
c.... outside the current yield surface.
c      dlam  = abs(dlam*(1.0 - tol))
c.... compute the elasto-plasticity matrix dmat
      c4 = two3*d(8)*sps/(1.-dlam*two3*d(8))+c3
      c4 = 1./c4
      call pzero(dmat,9)
      dmat(1,1) = x1 - c4*psh(1)**2
      dmat(1,2) = x2 - c4*psh(1)*psh(2)
      dmat(1,3) =    - c4*psh(1)*psh(3)
      dmat(2,2) = x1 - c4*psh(2)**2
      dmat(2,3) =    - c4*psh(2)*psh(3)
      dmat(3,3) = x3 - c4*psh(3)**2
      dmat(2,1) = dmat(1,2)
      dmat(3,1) = dmat(1,3)
      dmat(3,2) = dmat(2,3)
c.... set stresses (sigtr is moved into h2!)
      sigtr(1) = sig(1)
      sigtr(2) = sig(2)
      sigtr(3) = sig(3)
c.... set total plastic strains
      epn(1) = epn(1) + dlam*ps(1)
      epn(2) = epn(2) + dlam*ps(2)
      epn(3) = epn(3) + dlam*ps(3)
c...  set effective plastic strain
      epq  = epg
c...  set plastic mutilplier
      dlam0 = dlam
      return
      end
c
      subroutine pstr22(dgrad,sig,psig,lin)
c-----------------------------------------------------------------------
c.....1.PK stresses
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension dgrad(4),sig(3),psig(4)
      if(lin.eq.0) then
        psig(1)=sig(1)
        psig(2)=sig(2)
        psig(3)=sig(3)
        psig(4)=sig(3)
      else
        uxkx=dgrad(1)
        uykx=dgrad(2)
        uxky=dgrad(3)
        uyky=dgrad(4)
        psig(1)=(1.d0+uxkx)*sig(1)+uxky*sig(3) 
        psig(2)=(1.d0+uyky)*sig(2)+uykx*sig(3) 
        psig(3)=(1.d0+uxkx)*sig(3)+uxky*sig(2) 
        psig(4)=(1.d0+uyky)*sig(3)+uykx*sig(1) 
      end if
      return
      end 
c
      subroutine qload22(q,ql,numel,n,mqloa,propq)
c----------------------------------------------------------
c.... set loads from macro qloa
c----------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension q(numel,10),ql(2)
      call pzero(ql,2)
      if(mqloa.ne.1) then
        do i=1,2 
          ql(i) = q(n,i)*propq 
        end do
      end if  
      return
      end
c
