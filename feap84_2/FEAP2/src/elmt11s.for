      subroutine elmt11(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c-----------------------------------------------------------------------
c     (2+1) node geometrical nonlinear 3-D-beam element  
c     eccentrical Timoshenko beam theory                 
c     St. Venant torsion, mod. rotation                  
c     with coupling N-M and Q-MT                         
c                                                        
c--------------------------------------------------------
c     element number used for plot                       
c--------------------------------------------------------
c
c     # Verschiebungen und Schnittgrößen beziehen sich
c       auf gewähltes KOS
c     # Schnittgrößen können mit y_o,z_o auf andere Achse
c       transformiert werden. 
c     # Lasten beziehen sich auf gewähltes KOS
c
c     # 1.PK in sigtrats enthalten aber auskommentiert
c       d.h. 2.PK-groessen werden ausgegeben.
c
C     # nili Terme in  in phix noch enthalten
c     
c     # Schnittgrößen aus R = int B^T sigma
c
c     # Shear correction factor kappa fehlt. Gesetzt=1
c
c     # Shear correction factor FE - Bischoff/Bletzinger
c
c       kappaGA=kappaGA/(1+kappa*l^2/12*GA/EI), 
c
c--------------------------------------------------------|
c     Version 1: mit 3 Knoten: 3.Kn.= Hilfsknoten        |
c     definition of local y-axis by 3. node              |
c     e_yl = ((x)_3 - (x)_1) x e_xl                      |
c     e_zl = e_xl x e_yl                                 |
c--------------------------------------------------------|
c     Version 2: mit 2 Knoten: 3. Knoten entfaellt       |
c     Richtung 3 = (0,0,1) = e_Z                         |
c     e_yl =  e_xl x e_Z                                 |
c     e_zl = e_xl x e_yl                                 |
c     Sonderfall: e_xl = +-e_Z                           |
c     e_yl = (0,-1,0)  willkuerlich                       |
c--------------------------------------------------------|
c     D-Feldbelegung                                     |
c 1... d1( 1) = A         2... d( 1) = EA                |
c      d1( 2) = A              d( 2) = GA                |
c      d1( 3) = A              d( 3) = rhoA              |
c      d1( 4) = I_x            d( 4) = GI_x              |
c      d1( 5) = I_y            d( 5) = EI_y              |
c      d1( 6) = I_z            d( 6) = EI_z              |
c      d1( 7) = p_x            d( 7) = p_x               |
c      d1( 8) = p_y            d( 8) = p_y               |
c      d1( 9) = p_z            d( 9) = p_z               |
c      d1(10)= alpha           d(10)= alpha              |
c                                                        |
c 3... d2( 1) = lin            d(11) = lin               |
c      d2( 2) = rho            d(12) = EI_yz             |
c      d2( 3) = I_yz(=+intyzda)                          | 
c      d2( 4) = y_s            d(13) = y_s               |
c      d2( 5) = z_s            d(14) = z_s               |
c      d2( 6) = y_m            d(15) = y_m               |
c      d2( 7) = z_m            d(16) = z_m               |
c      d2( 8) = y_o            d(17) = y_o               |
c      d2( 9) = z_o            d(18) = z_o               |
c      d2(10) = shear corr     d(19) = scfs(0/1)         |
c               on/off                                   |
c--------------------------------------------------------+
      USE bdata
      USE cdata 
      USE eldata  
      USE evdata 
      USE fornam 
      USE iofile 
      USE pdata6 
      USE pdata10 
      USE pltran 
      USE prlod 
      USE qload       
      USE strnam 
      implicit double precision (a-h,o-z)
      logical zoom
      dimension ix(*),ixl(5),ixld(4),xl(ndm,*),tl(*),xll(3,4),sig(12),
     1          d(*),ul(ndf,*),s(nst,*),p(nst),shp(2,2),gr(8),dy(3),
     2          bm(6,6),btd(6,6),dmat(6,6),d1(10),d2(10),xdp(3),yl(3,2)
      dimension h1(*),h2(*),h3(*)
      common /tran3d/  dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3
!$OMP THREADPRIVATE (/tran3d/)  
      data ixl /1,2,3,4,1/,ixld /1,2,3,1/,eps/1.0e-5/
c-----------------
      ielno = 11
c-----------------
c.... transfer to correct processor
      go to (1,2,3,3,3,3,2,2,2,2,3,2,3,2,2,2,2,2,2,2,2,22), isw
      return
c.... input material properties
1     if(ior.lt.0) write(*,2004)
      call dinput(d1,10)
      if(ior.lt.0) write(*,2005)
      call dinput(d2,11)
      if(ior.lt.0) write(  *,2000)(d1(i),i=1,10),(d2(i),i=1,10)
                   write(iow,2000)(d1(i),i=1,10),(d2(i),i=1,10)
c.....store in d, multiply with E,G
      d( 1) = d1( 1)*d1(3)  ! EA  
      d( 2) = d1( 2)*d1(3)  ! GA  
      d( 3) = d2( 2)*d1(3)  ! rhoA
      d( 4) = d1( 2)*d1(4)  ! GI_x
      d( 5) = d1( 1)*d1(5)  ! EI_y
      d( 6) = d1( 1)*d1(6)  ! EI_z
      d( 7) = d1( 7)        ! p_x 
      d( 8) = d1( 8)        ! p_y 
      d( 9) = d1( 9)        ! p_z 
      d(10) = d1(10)        ! alpha
      d(11) = d2( 1)        ! lin  
      d(12) = d1( 1)*d2(3)  ! EI_yz
      d(13) = d2( 4)        ! y_s     
      d(14) = d2( 5)        ! z_s     
      d(15) = d2( 6)        ! y_m     
      d(16) = d2( 7)        ! z_m     
      d(17) = d2( 8)        ! y_o     
      d(18) = d2( 9)        ! z_o     
      d(19) = d2(10)        ! FE shear correction 0=on,1=off     
c.... define node numbering for plot mesh routine, see pltord
      inord(ielno)   = 3
      ipord(1,ielno) = 1
      ipord(2,ielno) = 2
      ipord(3,ielno) = 1
c.... description of stresses  
      forsus( 1) =  '  N-FORCE N_x  '
      forsus( 2) =  '  Q-FORCE Q_y  '
      forsus( 3) =  '  Q-FORCE Q_z  '
      forsus( 4) =  '  MOMENT  M_x  '
      forsus( 5) =  '  MOMENT  M_y  '
      forsus( 6) =  '  MOMENT  M_z  '
      forsus( 7) =  '  LOCAL   COS  '
      do i = 8,11         
         forsus(i) =  ' '
      enddo
2     return
c.... stiffness matrix and others
3     lin = d(11)
c.... compute dir. cosin terms for tr(3,3) and member length
      dx1 = xl(1,2) - xl(1,1)
      dx2 = xl(2,2) - xl(2,1)
      dx3 = xl(3,2) - xl(3,1)
      dl  = dsqrt(dx1*dx1+dx2*dx2+dx3*dx3)
      if(isw.eq.5) go to 5
c.... unit vector e_xL
      dx1 = dx1/dl
      dx2 = dx2/dl
      dx3 = dx3/dl
c.....Hilfsvektor x_3 - x_1 oder e_Z
      if(nel.eq.3) then ! 3 Knoten
        dz1 = xl(1,3) - xl(1,1)
        dz2 = xl(2,3) - xl(2,1)
        dz3 = xl(3,3) - xl(3,1)
      elseif(nel.eq.2) then ! 2 Knoten e_z = 0,0,1
        dz1 = 0.d0
        dz2 = 0.d0
        dz3 = 1.d0
      else
        stop '3D-beam-element with 2 or 3 nodes'
      endif
c.....e_yL =  (x_3-x_1) x e_xL  oder e_xL x e_Z = -e_Z x e_xL
      fac = -1.d0      
      if(nel.eq.3) fac = 1.d0
      dy1 = fac*(dz2*dx3 - dz3*dx2)
      dy2 = fac*(dz3*dx1 - dz1*dx3)
      dy3 = fac*(dz1*dx2 - dz2*dx1)
c.... length of e_yL
      dn  = dsqrt(dy1*dy1+dy2*dy2+dy3*dy3)
c.... special case e_xL = +-e_Z
      deps = 1.e-7*dl
      if(dn.lt.deps) then
        dn  =  1.d0
        dy1 =  0.d0
        dy2 = -1.d0
        dy3 =  0.d0
      endif
c.... unit vector e_yL
      dy1 = dy1/dn
      dy2 = dy2/dn
      dy3 = dy3/dn
c.... unit vector e_zL = e_xL x e_yL
      dz1 = dx2*dy3 - dx3*dy2
      dz2 = dx3*dy1 - dx1*dy3
      dz3 = dx1*dy2 - dx2*dy1
c.....rotate local coordinate system with angle alpha
      alpha = d(10)
      if(alpha.ne.0.d0) then
        if(nel.eq.3) stop '3D-beam-element: alpha only for 2 nodes'
        call rotate11(alpha)
      endif
c.... local displacement vector
      call trans11(s,p,ul,nst,ndf,2)
c.... shape functions  for one point integration
        shp(1,1) = -1.0d0/dl
        shp(1,2) =  1.0d0/dl
        shp(2,1) =  0.5d0
        shp(2,2) =  0.5d0
c.... material-matrix
      call dmat11(d,dmat,dl)
c.... stresses
      call stress11(sig,shp,ul,gr,dmat,lin)
c.... element load vector local
      call qload11(qxl,qyl,qzl,d,aqloa,numel,n,mqloa,propq,prop,isw)
      p(1) =  qxl*0.5d0*dl
      p(2) =  qyl*0.5d0*dl
      p(3) =  qzl*0.5d0*dl
      p(7) =  p(1)
      p(8) =  p(2)
      p(9) =  p(3)     
c.... loop for stiffness matrix
      i1=0
        do 31 ii=1,2
c....   B-matrix
          call bmat11(bm,shp,gr,ii,lin)
c....   residual G = P - Bt*S and matrix Bt*D
          do 32 i = 1,ndf
             do  33 k = 1,ndf
               btd(i,k) = 0.0
               p(i1+i) = p(i1+i) - bm(k,i)*sig(k)*dl
             do 34 j = 1,ndf
                btd(i,k) = btd(i,k)+bm(j,i) * dmat(j,k)
34           continue
33           continue
32        continue
        if(isw.eq.4.or.isw.eq.6.or.isw.eq.13) goto 39
          j1 = i1
          do 35 jj = ii,2
            call bmat11(bm,shp,gr,jj,lin)
            do 36  i = 1,ndf
              do 37  j = 1,ndf
                    do 38       k = 1,ndf
                      s(i1+i,j1+j) = s(i1+i,j1+j) + btd(i,k)*bm(k,j)*dl
38            continue
37            continue
36          continue
c....     K sigma
            if(lin.ne.0) call ksigma11(sig,shp,s,nst,dl,ii,jj,i1,j1)
            j1 = j1 + ndf
35        continue
39        i1 = i1 + ndf
31      continue
      if(isw.eq.6) goto 300
      if(isw.eq.4.or.isw.eq.13) goto 4
c.... lower part of stiffness matrix
      do i = 1,6
        do j = 1,6
          s(i+ndf,j) = s(j,i+ndf)
        enddo
      enddo
      call trans11(s,p,ul,nst,ndf,1)
300   call trans11(s,p,ul,nst,ndf,3)
      return
c.... output forces R = - G  = - (f - int bTS dv)  (nonlinear!)
4     continue
c.... stresses from R = -G
c.....copy 1-12, change 7-12  due to R = -G
      do i = 1,6 
        sig(i)   =  p(i)
        sig(i+6) = -p(i+6)
      enddo
c
c     Transformation
      call sigtra11(sig,d,gr,lin)
      if(isw.eq.13) goto 13
      mct = mct - 1
      if (mct.le.0) then
        if(ior.lt.0) write(*  ,2002) o,head,d(17),d(18)
                     write(iow,2002) o,head,d(17),d(18)
        mct = 50
      endif
                   write(iow,2003) n,ma,(sig(i), i=1, 6)
                   write(iow,2006)      (sig(i), i=7,12)
      if(ior.lt.0) write(*  ,2003) n,ma,(sig(i), i=1, 6)
      if(ior.lt.0) write(*  ,2006)      (sig(i), i=7,12)
      return
c
5     continue
c.... mass matrix (lumped only)
      pm   = d(3)*dl*0.5d0
      pr   = pm *0.5d0*(d(5)+d(6))/d(2)
      p(1) = pm
      p(2) = pm
      p(3) = pm
      p(4) = pr
      p(5) = pr 
      p(6) = pr 
      p(ndf+1) = pm
      p(ndf+2) = pm
      p(ndf+3) = pm
      p(ndf+4) = pr
      p(ndf+5) = pr
      p(ndf+6) = pr
      return
c.... plot stress resultants
13    klayf = 1
      if(iplma(ma).eq.0)       return ! only if MATN
      if(nfp.lt.1.or.nfp.gt.7) return
      nd2 = ndf + nfp
      s1  = sig(nfp)
      s2  = sig(nd2)
      call qload11(qxl,qyl,qzl,d,aqloa,numel,n,mqloa,propq,prop,isw)
      if(flfp) then
        if(nfp.eq.7) return
          ccfp  = max(abs(s1),abs(s2))
          ccfp1 = max(s1,s2)
          ccfp2 = min(s1,s2)
          xmaxf = max(xmaxf,ccfp1)
          xminf = min(xminf,ccfp2)
          cfp  = max(cfp,ccfp)
      else
c....   plot local axis
        if(nfp.eq.7) then
          call pppcol(2)
          call plloco11(xl,ndm,dl)
          return
        endif
c....   Randwerte
        if(abs(s1).lt.eps) s1 = 0.d0
        if(abs(s2).lt.eps) s2 = 0.d0
          if(abs(s1).lt.eps.and.abs(s2).lt.eps) return
c....   normal vector
        if(ifor.eq.12) then
          dp1 =  dy1
          dp2 =  dy2
          dp3 =  dy3
        elseif(ifor.eq.-12) then
          dp1 = -dy1
          dp2 = -dy2
          dp3 = -dy3
        elseif(ifor.eq.13) then
          dp1 =  dz1
          dp2 =  dz2
          dp3 =  dz3
        elseif(ifor.eq.-13) then
          dp1 = -dz1
          dp2 = -dz2
          dp3 = -dz3
        endif
c.....  schleife fuer alle Schnittgroessen
        nk = 6
c.....  Weginkrement
        do kk = 1,3
          dy(kk) = (xl(kk,2) - xl(kk,1))/float(nk)
        enddo
        dl1 =  dl/float(nk)
c
        do 135 ii = 1,nk
c.....  Randkoordinaten
        do kk = 1,3
          yl(kk,1) = xl(kk,1) + dy(kk) *(ii-1)
          yl(kk,2) = xl(kk,1) + dy(kk) * ii
        enddo
c.....  Randschnittgroessen
        dla = dl1 * (ii-1)
        dle = dl1 *  ii
        if(nfp.eq.1) then       ! N_x
          s1  = sig(1)  - qxl * dla
          s2  = sig(1)  - qxl * dle
        elseif(nfp.eq.2) then   ! Q_y
          s1  = sig(2)  - qyl * dla
          s2  = sig(2)  - qyl * dle
        elseif(nfp.eq.3) then   ! Q_z 
          s1  = sig(3)  - qzl * dla
          s2  = sig(3)  - qzl * dle
        elseif(nfp.eq.4) then   ! M_x
          s1 = sig(4)
          s2 = sig(4)
        elseif(nfp.eq.5) then   ! M_y
          s1  = sig(5) + sig(3) * dla - qzl * dla*dla*0.5d0
          s2  = sig(5) + sig(3) * dle - qzl * dle*dle*0.5d0
        else if(nfp.eq.6) then  ! M_z
          s1  = sig(6) - sig(2) * dla + qyl * dla*dla*0.5d0
          s2  = sig(6) - sig(2) * dle + qyl * dle*dle*0.5d0
        endif                                         
        if(abs(s1).lt.eps) s1 = 0.d0
        if(abs(s2).lt.eps) s2 = 0.d0
c
c.....  Berechne Farbe                                
        sm = (s1+s2)*0.5d0
        call pppcolf(sm)
c.....  plotte Schnittgroesse          
        if(s1.ge.0.d0.and.s2.ge.0.d0) then
          goto 130
        elseif(s1.le.0.d0.and.s2.le.0.d0) then
          goto 130
        else
c.....    Vorzeichenwechsel -> rechne zwischen (0,s1) und (dl,s2) lokal
          x1 = 0.d0
          x2 = dl1
          y1 = s1 
          y2 = s2 
c.....    gerade y = y1 + a1 * x   mit a1 = (y2-y1)/dl 
          a1  = (y2-y1)/dl1
c.....    Durchstosspunkt
          xd = -y1/a1
          if(xd.le.x1.or.xd.ge.x2) goto 135  ! ausserhalb
          xdp(1) =  yl(1,1) * (1.d0-xd/dl1) + yl(1,2) * xd/dl1
          xdp(2) =  yl(2,1) * (1.d0-xd/dl1) + yl(2,2) * xd/dl1
          xdp(3) =  yl(3,1) * (1.d0-xd/dl1) + yl(3,2) * xd/dl1
c....     plotte dreieck 1
          sm = s1*0.5d0
          call pppcolf(sm)
          call pzero(xll,12)
          xll(1,1) = yl(1,1)
          xll(2,1) = yl(2,1)
          xll(3,1) = yl(3,1)
          xll(1,2) = xdp(1)
          xll(2,2) = xdp(2)
          xll(3,2) = xdp(3)
          xll(1,3) = yl(1,1) + dp1*s1*cfp
          xll(2,3) = yl(2,1) + dp2*s1*cfp 
          xll(3,3) = yl(3,1) + dp3*s1*cfp
          call plxtrn(xll,tra,vr,3,4)
          if(zoom(xll,3,3)) call plot9s(ixld,xll,ndm,3)
c....     plotte dreieck 2
          sm = s2*0.5d0
          call pppcolf(sm)
          call pzero(xll,12)
          xll(1,1) = xdp(1)
          xll(2,1) = xdp(2)
          xll(3,1) = xdp(3)
          xll(1,2) = yl(1,2)
          xll(2,2) = yl(2,2)
          xll(3,2) = yl(3,2)
          xll(1,3) = yl(1,2) + dp1*s2*cfp
          xll(2,3) = yl(2,2) + dp2*s2*cfp 
          xll(3,3) = yl(3,2) + dp3*s2*cfp
          call plxtrn(xll,tra,vr,3,4)
          if(zoom(xll,3,3)) call plot9s(ixld,xll,ndm,3)
          goto 135
        endif
c....   plotte trapez
130     xll(1,1) = yl(1,1)
        xll(2,1) = yl(2,1)
        xll(3,1) = yl(3,1)
        xll(1,2) = yl(1,2)
        xll(2,2) = yl(2,2)
        xll(3,2) = yl(3,2)
        xll(1,3) = yl(1,2) + dp1*s2*cfp
        xll(2,3) = yl(2,2) + dp2*s2*cfp 
        xll(3,3) = yl(3,2) + dp3*s2*cfp
        xll(1,4) = yl(1,1) + dp1*s1*cfp
        xll(2,4) = yl(2,1) + dp2*s1*cfp 
        xll(3,4) = yl(3,1) + dp3*s1*cfp
        call plxtrn(xll,tra,vr,3,4)
        if(zoom(xll,3,4)) call plot9s(ixl,xll,ndm,4)
135     continue
      endif
      return

c.... loadvector for QLOA
22    continue
c.... compute dir. cosin terms for tr(3,3) and member length
      dx1 = xl(1,2) - xl(1,1)
      dx2 = xl(2,2) - xl(2,1)
      dx3 = xl(3,2) - xl(3,1)
      dl  = dsqrt(dx1*dx1+dx2*dx2+dx3*dx3)
c.... unit vector e_xL
      dx1 = dx1/dl
      dx2 = dx2/dl
      dx3 = dx3/dl
c.....Hilfsvektor x_3 - x_1 oder e_Z
      if(nel.eq.3) then ! 3 Knoten
        dz1 = xl(1,3) - xl(1,1)
        dz2 = xl(2,3) - xl(2,1)
        dz3 = xl(3,3) - xl(3,1)
      elseif(nel.eq.2) then ! 2 Knoten e_z = 0,0,1
        dz1 = 0.d0
        dz2 = 0.d0
        dz3 = 1.d0
      else
        stop '3D-beam-element with 2 or 3 nodes'
      endif
c.....e_yL =  (x_3-x_1) x e_xL  oder e_xL x e_Z = -e_Z x e_xL
      fac = -1.d0      
      if(nel.eq.3) fac = 1.d0
      dy1 = fac*(dz2*dx3 - dz3*dx2)
      dy2 = fac*(dz3*dx1 - dz1*dx3)
      dy3 = fac*(dz1*dx2 - dz2*dx1)
c.... length of e_yL
      dn  = dsqrt(dy1*dy1+dy2*dy2+dy3*dy3)
c.... special case e_xL = +-e_Z
      if(dn.eq.0.d0) then
        dn  =  1.d0
        dy1 =  0.d0
        dy2 = -1.d0
        dy3 =  0.d0
      endif
c.... unit vector e_yL
      dy1 = dy1/dn
      dy2 = dy2/dn
      dy3 = dy3/dn
c.... unit vector e_zL = e_xL x e_yL
      dz1 = dx2*dy3 - dx3*dy2
      dz2 = dx3*dy1 - dx1*dy3
      dz3 = dx1*dy2 - dx2*dy1
c.....rotate local coordinate system with angle alpha
      alpha = d(10)
      if(alpha.ne.0.d0) then
        if(nel.eq.3) stop '3D-beam-element: alpha only for 2 nodes'
        call rotate11(alpha)
      endif
c.... element load vector local
      call qload11(qxl,qyl,qzl,d,aqloa,numel,n,mqloa,propq,prop,isw)
      p(1) =  qxl*0.5d0*dl
      p(2) =  qyl*0.5d0*dl
      p(3) =  qzl*0.5d0*dl
      p(7) =  p(1)
      p(8) =  p(2)
      p(9) =  p(3)     

      call trans11(s,p,ul,nst,ndf,3)
      return

c.... format statements
2000  format(//5x,'3-D NONLINEAR  TIMOSHENKO BEAM  ELEMENT'/
     + 10x,'Elastic modulus      E ..........',e15.5/
     + 10x,'Shear   modulus      G ..........',e15.5/
     + 10x,'Area                 A...........',e15.5/
     + 10x,'Moment of inertia    I_x ........',e15.5/
     + 10x,'Moment of inertia    I_y ........',e15.5/
     + 10x,'Moment of inertia    I_z ........',e15.5/
     + 10x,'constant load        p_xl .......',e15.5/
     + 10x,'constant load        p_yl .......',e15.5/
     + 10x,'constant load        p_zl .......',e15.5/
     + 10x,'local rot.angle y->z (deg.)......',e15.5/ 
     + 10x,'geom.nonlinear (0=f,1=t) ........',e15.5/
     + 10x,'density = gamma/g    rho.. ......',e15.5/
     + 10x,'Moment of inertia    I_yz .......',e15.5/
     + 10x,'eccentricity         y_S ........',e15.5/
     + 10x,'eccentricity         z_S ........',e15.5/
     + 10x,'eccentricity         y_M ........',e15.5/
     + 10x,'eccentricity         z_M ........',e15.5/
     + 10x,'eccentricity         y_A (output)',e15.5/
     + 10x,'eccentricity         z_A (output)',e15.5/
     + 10x,'FE shear correction on/off=0/1...',e15.5)
2002  format(a1,20a4/ 
     1 3x,'3-D BEAM - STRESS RESULTANTS (2.PK) at y=',f8.3,'  z=',f8.3/,
     2  '   el mat','    N_x    ',1x,'    Q_y    ',1x,'    Q_z    ',1x,
     3              '    M_T    ',1x,'    M_y    ',1x,'    M_z    ')
2003  format(1x,i3,1x,i2,6(1x,g11.5))
2004  format('Input:E,G,A,I_x,I_y,I_z,p_x,p_y,p_z,alpha>',$)
2005  format('Input:rho,lin,I_yz,y_S,z_S,y_M,z_M,y_A,z_A,scf>',$)
2006  format(7x,6(1x,g11.5))
      end
c
      subroutine trans11(s,p,ul,nst,ndf,isw)
c------------------------------------------------------------------+
c.... isw: 1  transform matrix s(nst,nst)   Sg = Tt * Sl * T       |
c          2  transform vector ul(ndf,1)    ul = T  * ug           |
c          3  transform vector P(nst)       Pg = Tt * Pl           |
c          4  transform vector ul(ndf,1)    Ug = Tt * Ul           |
c------------------------------------------------------------------+
      implicit double precision (a-h,o-z)
      dimension s(nst,*),ss(6),st(6),ul(ndf,*),p(nst)
      common/tran3d/ tr(3,3)
!$OMP THREADPRIVATE (/tran3d/)  
c.... skip for an identity transformation
      if(tr(1,1)+tr(2,2)+tr(3,3).ge.2.999999d+00) return
      go to (1,2,3,4) isw
1     continue
c.... Sg = tr * Sl * T
      nsiz = ndf + ndf
c.... postmultiply local stiffness by transformation array
        j1 = 0
        do 120 k = 1,2
          do 110 i = 1,nsiz
            do 100 j = 1,6
              ss(j) = s(i,j+j1)
100         continue
            j2 = j1 + 3
            do 101 j = 1,3
              s(i,j+j1) = ss(1)*tr(1,j) + ss(2)*tr(2,j) + ss(3)*tr(3,j)
              s(i,j+j2) = ss(4)*tr(1,j) + ss(5)*tr(2,j) + ss(6)*tr(3,j)
101         continue
110       continue
          j1 = j1 + ndf
120     continue
c.... premultiply result by the transpose of the transformation array
        j1 = 0
        do 150 k = 1,2
          do 140 i = 1,nsiz
            do 130 j = 1,6
              ss(j) = s(j+j1,i)
130         continue
            j2 = j1 + 3
            do 131 j = 1,3
              s(j+j1,i) = tr(1,j)*ss(1) + tr(2,j)*ss(2) + tr(3,j)*ss(3)
              s(j+j2,i) = tr(1,j)*ss(4) + tr(2,j)*ss(5) + tr(3,j)*ss(6)
131         continue
140       continue
          j1 = j1 + ndf
150     continue
      return
2     continue
c.... Ul = T * Ug, (U_1l in ss, U_2l in st)
      do 200 i = 1,3
      ss(i)   = tr(i,1)*ul(1,1) + tr(i,2)*ul(2,1) + tr(i,3)*ul(3,1)
      ss(i+3) = tr(i,1)*ul(4,1) + tr(i,2)*ul(5,1) + tr(i,3)*ul(6,1)
      st(i)   = tr(i,1)*ul(1,2) + tr(i,2)*ul(2,2) + tr(i,3)*ul(3,2)
200   st(i+3) = tr(i,1)*ul(4,2) + tr(i,2)*ul(5,2) + tr(i,3)*ul(6,2)
      do 210 i = 1,6
      ul(i,1)   = ss(i)
210   ul(i,2)   = st(i)
      return
3     continue
c.... Pg = tr * Pl, (P_1g in ss, P_2g in st)
      do 320 i = 1,3
      ss(i)   = tr(1,i)*p(1)  + tr(2,i)*p(2)  + tr(3,i)*p(3)
      ss(i+3) = tr(1,i)*p(4)  + tr(2,i)*p(5)  + tr(3,i)*p(6)
      st(i)   = tr(1,i)*p(7)  + tr(2,i)*p(8)  + tr(3,i)*p(9)
320   st(i+3) = tr(1,i)*p(10) + tr(2,i)*p(11) + tr(3,i)*p(12)
      do 330 i = 1,6
      p(i)   = ss(i)
330   p(i+6) = st(i)
      return
4     continue
c.... Ug = Tt * Ul, (U_1g in ss, U_2g in st)
      do 400 i = 1,3
      ss(i)   = tr(1,i)*ul(1,1) + tr(2,i)*ul(2,1) + tr(3,i)*ul(3,1)
      ss(i+3) = tr(1,i)*ul(4,1) + tr(2,i)*ul(5,1) + tr(3,i)*ul(6,1)
      st(i)   = tr(1,i)*ul(1,2) + tr(2,i)*ul(2,2) + tr(3,i)*ul(3,2)
400   st(i+3) = tr(1,i)*ul(4,2) + tr(2,i)*ul(5,2) + tr(3,i)*ul(6,2)
      do 410 i = 1,6
      ul(i,1)   = ss(i)
410   ul(i,2)   = st(i)
      return
      end
c
      subroutine bmat11(bm,shp,gr,k,lin)
c----------------------------------------------------------------------+
c.... B-Matrix for eccentrical nonlinear 3d-Timoshenko beam
c----------------------------------------------------------------------+
      implicit double precision(a-h,o-z)
      dimension bm(6,6),shp(2,2),gr(8) 
      call pzero(bm,36)
      dn  = shp(2,k)
      dnx = shp(1,k)
      ukx  = gr(1) 
      vkx  = gr(2) 
      wkx  = gr(3) 
      pyx  = gr(4) 
      pzx  = gr(5) 
      pxx  = gr(6) 
      py   = gr(7) 
      pz   = gr(8) 
c.... b-mat
      bm(1,1) = dnx
      bm(2,2) = dnx
      bm(2,6) = -dn
      bm(3,3) = dnx
      bm(3,5) = dn
      bm(4,4) = dnx
      bm(5,5) = dnx
      bm(6,6) = -dnx
      if(lin.ne.0) then
        bm(1,2) = vkx * dnx
        bm(1,3) = wkx * dnx
        bm(4,5) = 0.5d0*( pz*dnx-pzx*dn)
        bm(4,6) = 0.5d0*(-py*dnx+pyx*dn)
      endif
      return
      end
c
      subroutine stress11(sig,shp,ul,gr,dmat,lin)
c----------------------------------------------------------------------+
c.... local stresses for eccentrical nonlinear 3d-Timoshenko beam
c----------------------------------------------------------------------+
      implicit double precision(a-h,o-z)
      dimension sig(12),eps(6),shp(2,2),gr(8),dmat(6,6),ul(6,*)
      call pzero(sig,12)
      ukx = shp(1,1)*ul(1,1) + shp(1,2)*ul(1,2)
      vkx = shp(1,1)*ul(2,1) + shp(1,2)*ul(2,2)
      wkx = shp(1,1)*ul(3,1) + shp(1,2)*ul(3,2)
      pxx = shp(1,1)*ul(4,1) + shp(1,2)*ul(4,2)
      pyx = shp(1,1)*ul(5,1) + shp(1,2)*ul(5,2)
      pzx = shp(1,1)*ul(6,1) + shp(1,2)*ul(6,2)
      py  = shp(2,1)*ul(5,1) + shp(2,2)*ul(5,2)
      pz  = shp(2,1)*ul(6,1) + shp(2,2)*ul(6,2)
c.... save data for bmat
      gr(1) = ukx
      gr(2) = vkx
      gr(3) = wkx
      gr(4) = pyx
      gr(5) = pzx
      gr(6) = pxx
      gr(7) = py
      gr(8) = pz
c.... strains     
      eps(1) = ukx 
      eps(2) = vkx - pz
      eps(3) = wkx + py
      eps(4) = pxx
      eps(5) = pyx
      eps(6) = -pzx
      if(lin.ne.0) then
        eps(1) = eps(1) + 0.5d0 * (vkx*vkx + wkx*wkx)
        eps(4) = eps(4) + 0.5d0 * (pyx*pz  - pzx*py)
      endif
c.... stresses  N_x, Q_y, Q_z, M_x, M_y, M_z
      do i=1,6
        do k=1,6
        sig(i) = sig(i) + dmat(i,k) * eps(k)
        enddo    
      enddo   
      return
      end
c
      subroutine ksigma11(sig,shp,s,nst,dl,ii,kk,i1,k1)
c----------------------------------------------------------------------+
c.... matrix K_sigma for eccentrical nonlinear 3d-Timoshenko beam
c----------------------------------------------------------------------+
      implicit double precision(a-h,o-z)
      dimension sig(12),shp(2,2),s(nst,nst)
c.... shape functions and derivatives
        dni   = shp(2,ii)
        dndxi = shp(1,ii)
        dnk   = shp(2,kk)
        dndxk = shp(1,kk)
c.... added terms
        s(i1+2,k1+2) = s(i1+2,k1+2)+dl* dndxi * dndxk * sig(1)
        s(i1+3,k1+3) = s(i1+3,k1+3)+dl* dndxi * dndxk * sig(1)
        s(i1+5,k1+6) = s(i1+5,k1+6)+dl*(dndxi*dnk-dni*dndxk)*sig(4)
        s(i1+6,k1+5) = s(i1+6,k1+5)+dl*(dni*dndxk-dndxi*dnk)*sig(4)
      return
      end
c
      subroutine dmat11(d,dmat,dl)
c----------------------------------------------------------------------+
c.... D-matrix for eccentrical nonlinear 3d-Timoshenko beam
c     with n-m + q-mt coupling 
c     dmat(2,4) zm!
c     dmat(3,4) ym!
c     dmat(4,4) ym,zm!
c     including shear correction on all GA terms! 
c----------------------------------------------------------------------+
      implicit double precision(a-h,o-z)
      dimension d(*),dmat(6,6)
      call pzero(dmat,36)
      scfy = 1.d0  
      scfz = 1.d0
      if(d(19).eq.0.d0) then 
        scfy = scfy/(1.d0+scfy*dl*dl/12.d0*d(2)/d(6)) 
        scfz = scfz/(1.d0+scfz*dl*dl/12.d0*d(2)/d(5)) 
      end if
      dmat(1,1) =  d(1)
      dmat(1,5) =  d(1) * d(14)
      dmat(1,6) =  d(1) * d(13)
      dmat(2,2) =  d(2) * scfy
      dmat(2,4) = -d(2) * d(16) * scfy
      dmat(3,3) =  d(2) * scfz 
      dmat(3,4) =  d(2) * d(15) * scfz
      dmat(4,4) =  d(4) + d(2)*(d(15)*d(15)*scfz + d(16)*d(16)* scfy) 
      dmat(5,5) =  d(5) + d(14) * d(14) * d(1)
      dmat(5,6) = d(12) + d(13) * d(14) * d(1)
      dmat(6,6) =  d(6) + d(13) * d(13) * d(1)
c...  symmetrie
      dmat(4,2) = dmat(2,4)    
      dmat(4,3) = dmat(3,4)    
      dmat(5,1) = dmat(1,5)    
      dmat(6,1) = dmat(1,6)
      dmat(6,5) = dmat(5,6)
c....
      return
      end
c
      subroutine sigtra11(sig,d,gr,lin)
c----------------------------------------------------------------------+
c.... transform local stresses for  eccentrical 3d-beam
c       d(17) = y_o     
c       d(18) = z_o     
c----------------------------------------------------------------------+
      implicit double precision(a-h,o-z)
      dimension sig(12),d(*),gr(8)
c     Spannungstransformation  ->1.PK
c      if(lin.ne.0) then
c        vkx    = gr(2)
c        wkx    = gr(3)
c        sig(2) = sig(2) + sig(1) * vkx
c        sig(3) = sig(3) + sig(1) * wkx
c        sig(8) = sig(8) + sig(7) * vkx
c        sig(9) = sig(9) + sig(7) * wkx
c      endif

c     Spannungstransformation auf Ausgabeachsen

      pos = dabs(d(17)) + dabs(d(18))
      if(pos.gt.0.d0) then
c       M_x = M_x0 + Q_y0*z_o - Q_z0*y_o
        sig( 4) =  sig(4) + sig(2)*d(18) - sig(3)*d(17)
        sig(10) = sig(10) + sig(8)*d(18) - sig(9)*d(17)
c       M_y = M_y0 - N_0*z_o
        sig( 5) =  sig(5) - sig(1)*d(18)                   
        sig(11) = sig(11) - sig(7)*d(18)                    
c       M_z = M_z0 + N_0*y_o
        sig( 6) =  sig(6) + sig(1)*d(17)                   
        sig(12) = sig(12) + sig(7)*d(17)                    
      endif
      return
      end        
c
      subroutine rotate11(alpha)
c--------------------------------------------------------+
c.... drehe KOS  e_yR =  c * e_yL  +  s * e_zL
c....            e_zR = -s * e_yL  +  c * e_zL
c--------------------------------------------------------+
      implicit double precision (a-h,o-z)
      common /tran3d/ dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3
!$OMP THREADPRIVATE (/tran3d/)  
c.... Winkel
      ra   = alpha*datan(1.d0)/45.d0
      cs   = dcos(ra)
      sn   = dsin(ra)
c.... e_yL rotated
      dhy1 = cs * dy1 + sn * dz1
      dhy2 = cs * dy2 + sn * dz2
      dhy3 = cs * dy3 + sn * dz3
c.... e_zL rotated
      dhz1 =-sn * dy1 + cs * dz1
      dhz2 =-sn * dy2 + cs * dz2
      dhz3 =-sn * dy3 + cs * dz3
c.... update
      dy1  = dhy1 
      dy2  = dhy2 
      dy3  = dhy3 
c
      dz1  = dhz1 
      dz2  = dhz2 
      dz3  = dhz3 
      return
      end
c
      subroutine plloco11(xl,ndm,dl)
c--------------------------------------------------------+
c.... Plot local basis
c--------------------------------------------------------+
      implicit double precision (a-h,o-z)
      dimension x0(3),xl(ndm,*)
c...  plot at midpoint of element (global coordinates!)
      x0(1) = 0.5d0*( xl(1,2) + xl(1,1) )
      x0(2) = 0.5d0*( xl(2,2) + xl(2,1) )
      x0(3) = 0.5d0*( xl(3,2) + xl(3,1) )
c.... length
      xm = 0.3d0 * dl
c...  plot axis 
      call pppcol(2)
      call pltaxs11(x0,ndm,xm)
      return
      end
c
      subroutine pltaxs11(x0,ndm,xm)
c----------------------------------------------------------------------
c.... draw vectors for axes, including rot-macro                      |
c----------------------------------------------------------------------
      USE pdata2
      USE pltran
      implicit double precision (a-h,o-z)
      logical zoom
      dimension xx(3,5),x0(3),x0z(3)
      common /tran3d/  tr(3,3)
!$OMP THREADPRIVATE (/tran3d/)  
c.... test if vectors are in plot-region
      x0z(1) = x0(1)
      x0z(2) = x0(2)
      x0z(3) = x0(3)
      call plxtrn(x0z,tra,vr,3,1)
      if(zoom(x0z,3,1)) then
c....   perspective projecton of axes
        do 120 m = 1,ndm
          call pzero(xx,15)
          do 100 n = 1,ndm
            fac1    = tr(m,n)*xm
            xx(n,1) = x0(n)
            xx(n,2) = xx(n,1) + fac1
100         xx(n,5) = xx(n,2)
          fac1 = tr(m,1)*xm
          fac2 = tr(m,2)*xm
          fac3 = tr(m,3)*xm
          xx(1,3) = xx(1,2) -.3*fac1 - .1*(fac2+fac3)
          xx(2,3) = xx(2,2) -.3*fac2 + .1*(fac1+fac3)
          xx(3,3) = xx(3,2) -.3*fac3 + .1*(fac1+fac2)
          xx(1,4) = xx(1,2) -.3*fac1 + .1*(fac2+fac3)
          xx(2,4) = xx(2,2) -.3*fac2 - .1*(fac1+fac3)
          xx(3,4) = xx(3,2) -.3*fac3 - .1*(fac1+fac2)
c.......  transform vector if necessary
                call plxtrn(xx,tra,vr,ndm,5)
c....     plot the vector
          call plotl(xx(1,1),xx(2,1),xx(3,1),3)
          call plotl(xx(1,2),xx(2,2),xx(3,2),2)
          call plotl(xx(1,2),xx(2,2),xx(3,2),ipgl)
          call plotl(xx(1,3),xx(2,3),xx(3,3),2)
          call plotl(xx(1,4),xx(2,4),xx(3,4),2)
          call plotl(xx(1,2),xx(2,2),xx(3,2),2)
          if(ipgl.eq.1) call clpan
          call plotl(xx(1,2),xx(2,2),xx(3,2),3)
120     continue
      endif
      return
      end
c
      subroutine qload11(qxl,qyl,qzl,d,q,numel,n,mqloa,propq,prop,isw)
c----------------------------------------------------------
c.... set loads from macro qloa or mate
c----------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension q(numel,10),d(*)
      if (isw.eq.22) then  ! load vector from QLOA
        qxl = 0.d0
        qyl = 0.d0 
        qzl = 0.d0 
        if(mqloa.ne.1) then 
          qxl = q(n,1)*propq 
          qyl = q(n,2)*propq 
          qzl = q(n,3)*propq 
        end if 
      else if (isw.eq.4.or.isw.eq.13) then ! load vector for STRE/FORC
        qxl = d(7)*prop
        qyl = d(8)*prop 
        qzl = d(9)*prop 
        if(mqloa.ne.1) then 
          qxl = qxl + q(n,1)*propq 
          qyl = qyl + q(n,2)*propq 
          qzl = qzl + q(n,3)*propq 
        end if 
      else ! load vector from MATE
        qxl = d(7)*prop
        qyl = d(8)*prop 
        qzl = d(9)*prop 
      end if  
      return
      end
