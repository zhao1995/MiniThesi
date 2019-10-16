      subroutine elmt04(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c--------------------------------------------------------+
c     (2+1) node 2nd order theory 3-D-beam element       |
c     Bernoulli theory, St. Venant torsion               |
c--------------------------------------------------------+
c                                                        |
c     Version 1: mit 3 Knoten: 3.Kn.= Hilfsknoten        |
c     definition of local y-axis by 3. node              |
c     e_yl = ((x)_3-(x)_1)  x  e_xl                      |
c     e_zl = e_xl x e_yl                                 |
c                                                        |
c     Version 2: mit 2 Knoten: 3. Knoten entfaellt       |
c     Richtung 3 = (0,0,1) = e_Z                         |
c     e_yl =  e_xl x e_Z                                 |
c     e_zl =  e_xl x e_yl                                |
c     Sonderfall: e_xl = +-e_Z                           |
c     e_yl = (0,-1,0)  willkürlich                       |
c--------------------------------------------------------|
c     D-Feldbelegung                                     |
c 1... d( 1) = E                                         |
c      d( 2) = G                                         |
c      d( 3) = A  ->EA                                   |
c      d( 4) = I_x->GI_x                                 |
c      d( 5) = I_y->EI_y                                 |
c      d( 6) = I_z->EI_z                                 |
c      d( 7) =      p_xL                                 |
c      d( 8) =      p_yL                                 |
c      d( 9) =      p_zL                                 |
c      d(10) =      p_xG  bezogen auf Stablänge im Raum  |
c      d(11) =      p_yG                                 |
c      d(12) =      p_zG                                 |
c      d(13) =      alpha Drehwinkel zu Hauptachse in Grad  
c                                                        | 
c 2... d(14) = Th.2O 0/1                                 |
c      d(15) = rho->rhoA                                 |
c                                                        | 
c--------------------------------------------------------+
c     Loads from QLOA
c
c      q01 =    p_xL                                 |
c      q02 =    p_yL                                 |
c      q03 =    p_zL                                 |
c      q04 =    p_xG  bezogen auf Stablänge im Raum  |
c      q05 =    p_yG                                 |
c      q06 =    p_zG                                 |
c
c--------------------------------------------------------+
c   * plot,forc from residual see elbfe2_4               |
c     forc,1 = N_x, ..,2 = Q_y, ..,3 = Q_z               |
c     forc,4 = M_x, ..,5 = M_y, ..,6 = M_z               |
c     forc,7 = local x-y-z COS                           |
c     fehl: plotte schnittgroessen blau = + rot = - und  |
c          Einbau fuer Trapezlasten                      |
c                                                        |
c   * Belastung: q_y=const. q_z = const.                 |
c     fehl: q_x=const. sowie alle als Trapezlast         |
c           q in globalen Richtungen                     |
c                                                        |
c--------------------------------------------------------|
c     residuum und schnittgroessen nur linear korrekt    |
c     linar  R = k*v-P  hier jedoch  R = k_T*v-P         |
c--------------------------------------------------------|
c                                                        |
c   * Stabilitaet (K+lambda*K_G)*Phi = 0                 |
c                                                        |
c   * Massenmatrix, nur lumped                           |
c                                                        |
c     (c) W. Wagner 8/92                                 |
c--------------------------------------------------------+
      USE bdata
      USE cdata
      USE eldata
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
     1    d(*),ul(ndf,*),s(nst,*),p(nst),xdp(3),dy(3),xlr(3,2)
      dimension h1(*),h2(*),h3(*)
      common /tran3d/  dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3
!$OMP THREADPRIVATE (/tran3d/)  
      data ixl /1,2,3,4,1/,ixld /1,2,3,1/,eps/1.0e-5/
      ielno = 4
c.... transfer to correct processor
      go to (1,2,3,3,3,3,2,2,2,2,3,2,3,2,2,2,2,2,2,2,2,3), isw
      return

c.... input material properties
1     if(ior.lt.0) write(*,3000)
3000  format('Input: E,G,A,Ix,Iy,Iz,qxL,qyL,qzL,qxG,qyG,qzG,a'/3x,'>',$)

      call dinput(d,13)

      if(ior.lt.0) write(*,3001)
3001  format('Input:lin, rho '/3x,'>',$)

      call dinput(d(14),2)

      if(ior.lt.0) write(*  ,2000) (d(i),i=1,15)
                         write(iow,2000) (d(i),i=1,15)
2000  format(5x,//' 3D - BERNOULLI-BEAM ELEMENT'//
     +  5x,'Elastic Modulus........... E .......',e13.5/
     +  5x,'Shear Modulus............. G .......',e13.5/
     +  5x,'Area...................... A .......',e13.5/
     +  5x,'Moment of inertia......... I_xx ....',e13.5/
     +  5x,'Moment of inertia......... I_yy ....',e13.5/
     +  5x,'Moment of inertia......... I_zz ....',e13.5/
     +  5x,'local  const. load........ q_xL.....',e13.5/
     +  5x,'local  const. load........ q_yL.....',e13.5/
     +  5x,'local  const. load........ q_zL.....',e13.5/
     +  5x,'global const. load........ q_xG.....',e13.5/
     +  5x,'global const. load........ q_yG.....',e13.5/
     +  5x,'global const. load........ q_zG.....',e13.5/
     +  5x,'local rot.angle y->z(deg.) alpha ...',e13.5/
     +  5x,'II.Order theory (0=false/1=true) ...',f6.1,/
     +  5x,'Density....................rho .....',e13.5)

      d(4)  = d(2)*d(4)  ! GI_x
      d(15) = d(3)*d(15) ! rhoA
      d(3)  = d(1)*d(3)  ! EA
      d(5)  = d(1)*d(5)  ! EI_y
      d(6)  = d(1)*d(6)  ! EI_z

c.... define node numbering for plot mesh routine, see pltord
      inord(ielno)   = 4
      ipord(1,ielno) = 1
      ipord(2,ielno) = 2
      ipord(3,ielno) = 2
      ipord(4,ielno) = 1

c.... description of stresses  
      forsus( 1) =  '  N-FORCE N_x  '
      forsus( 2) =  'Q-FORCE Q_y/R_y'
      forsus( 3) =  'Q-FORCE Q_z/R_z'
      forsus( 4) =  '  MOMENT  M_x  '
      forsus( 5) =  '  MOMENT  M_y  '
      forsus( 6) =  '  MOMENT  M_z  '
      forsus( 7) =  '  LOCAL   COS  '
      do i = 8,11         
         forsus(i) =  ' '
      enddo
2     return

c     tangent stiffness matrix
3     lin = d(14)

     
c.... compute direction cosine terms and member length
      dx1 = xl(1,2) - xl(1,1)
      dx2 = xl(2,2) - xl(2,1)
      dx3 = xl(3,2) - xl(3,1)
      dl  = sqrt(dx1*dx1+dx2*dx2+dx3*dx3)

      if(isw.eq.5) go to 5

      dl2 = dl*dl
      dl3 = dl*dl2

c.....e_xL
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
      dn  = sqrt(dy1*dy1+dy2*dy2+dy3*dy3)
c.... special case e_xL = +-e_Z
      deps = 1.e-7*dl
      if(dn.lt.deps) then
        dn  =  1.d0
        dy1 =  0.d0
        dy2 = -1.d0
        dy3 =  0.d0
      endif
c
      dy1 = dy1/dn
      dy2 = dy2/dn
      dy3 = dy3/dn

c.....e_zL
      dz1 = dx2*dy3 - dx3*dy2
      dz2 = dx3*dy1 - dx1*dy3
      dz3 = dx1*dy2 - dx2*dy1

c.....rotate local coordinate system with angle alpha
      alpha = d(13)
      if(alpha.ne.0.d0) then
        if(nel.eq.3) stop '3D-beam-element: alpha only for 2 nodes'
        call rotate04(alpha)
      endif
    
c.... local displacements
      call trans3b(s,p,ul,nst,ndf,2) 

c...  load vector   
      call qload04(qx,qy,qz,d,aqloa,propq,prop,numel,n,mqloa,isw)

c.... local load vector q_x
      if(qx.ne.0.d0) then
          p(1) = p(1) + qx * dl       * 0.5d0
          p(7) = p(7) + qx * dl       * 0.5d0
      endif

c.... local load vector q_y
      if(qy.ne.0.d0) then
          p(2) = p(2) + qy * dl       * 0.5d0
          p(8) = p(8) + qy * dl       * 0.5d0
          p(6) = p(6) + qy * dl * dl  /12.0d0
          p(12)= p(12)- qy * dl * dl  /12.0d0
      endif

c.... local load vector q_z
      if(qz.ne.0.d0) then
          p(3) = p(3) +  qz * dl      * 0.5d0
          p(9) = p(9) +  qz * dl      * 0.5d0
          p(5) = p(5) -  qz * dl * dl / 12.d0
          p(11)= p(11)+  qz * dl * dl / 12.d0
      endif

      if(isw.eq.22) go to 31


c.... compute axial force
       if(lin .eq. 1) then
        fn = d(3)*( ul(1,2)-ul(1,1))/dl
      end if
 
c.... compute axial stiffness terms
      i2 = ndf + 1
      s(1,1)   =  d(3)/dl
      s(i2,1)  = -s(1,1)
      s(1,i2)  = -s(1,1)
      s(i2,i2) =  s(1,1)

c.... compute torsional stiffness terms
      i2 = ndf + 4
      s(4,4)   = d(4)/dl
      s(i2,4)  = -s(4,4)
      s(4,i2)  = -s(4,4)
      s(i2,i2) =  s(4,4)

c.... compute the bending stiffness terms for z-displacements
      i1 = ndf + 3
      i2 = ndf + 5
      s(3,3)   = 12.d0*d(5)/dl3
      s(3,i1)  = -s(3,3)
      s(i1,3)  = -s(3,3)  ! fehlte ww 20.5.03
      s(i1,i1) =  s(3,3)
      s(5,i2)  = 2.d0*d(5)/dl
      s(5,5)   = 4.d0*d(5)/dl
      s(i2,5)  =  s(5,i2)
      s(i2,i2) =  s(5,5)
      s(3,5)   = -6.d0*d(5)/dl2
      s(5,3)   =  s(3,5)
      s(3,i2)  =  s(3,5)
      s(i2,3)  =  s(3,i2)
      s(5,i1)  = -s(3,5)
      s(i1,5)  =  s(5,i1)
      s(i1,i2) = -s(3,5)
      s(i2,i1) =  s(i1,i2)

c.... compute bending stiffness terms for y-displacement
      i1 = ndf + 2
      i2 = ndf + 6
      s(2,2)   = 12.*d(6)/dl3
      s(2,i1)  = -s(2,2)
      s(i1,2)  = -s(2,2)
      s(i1,i1) =  s(2,2)
      s(6,i2)  = 2.d0*d(6)/dl
      s(6,6)   = 4.d0*d(6)/dl
      s(i2,6)  =  s(6,i2)
      s(i2,i2) =  s(6,6)
      s(2,6)   = 6.d0*d(6)/dl2
      s(6,2)   =  s(2,6)
      s(2,i2)  =  s(2,6)
      s(i2,2)  =  s(2,i2)
      s(6,i1)  = -s(2,6)
      s(i1,6)  =  s(6,i1)
      s(i1,i2) = -s(2,6)
      s(i2,i1) =  s(i1,i2)

c.... compute nonlinear stiffness terms ( 2. order theory )
      if(lin.eq.1) then
        fn5  = 6.*fn/dl/5.
        fn10 = fn/10.
        fn15 = 2.*fn*dl/15.
        fn30 = fn*dl/30.
      else
        fn5  = 0.
        fn10 = 0.
        fn15 = 0.
        fn30 = 0.
      endif
      
      s(2,2)  = s(2,2)  + fn5
      s(2,6)  = s(2,6)  + fn10
      s(2,8)  = s(2,8)  - fn5
      s(2,12) = s(2,12) + fn10
      s(3,3)  = s(3,3)  + fn5
      s(3,5)  = s(3,5)  - fn10
      s(3,9)  = s(3,9)  - fn5
      s(3,11) = s(3,11) - fn10
      s(5,3)  = s(3,5)
      s(5,5)  = s(5,5)  + fn15
      s(5,9)  = s(5,9)  + fn10
      s(5,11) = s(5,11) - fn30
      s(6,2)  = s(2,6)
      s(6,6)  = s(6,6)  + fn15
      s(6,8)  = s(6,8)  - fn10
      s(6,12) = s(6,12) - fn30
      s(8,2)  = s(2,8)
      s(8,6)  = s(6,8)
      s(8,8)  = s(8,8)  + fn5
      s(8,12) = s(8,12) - fn10
      s(9,3)  = s(3,9)
      s(9,5)  = s(5,9)
      s(9,9)  = s(9,9)  + fn5
      s(9,11) = s(9,11) + fn10
      s(11,3) = s(3,11)
      s(11,5) = s(5,11)
      s(11,9) = s(9,11)
      s(11,11)= s(11,11)+ fn15
      s(12,2) = s(2,12)
      s(12,6) = s(6,12)
      s(12,8) = s(8,12)
      s(12,12)= s(12,12)+ fn15

      if(isw.eq.4 .or. isw.eq.13) go to 4
c
c.... calculate eigenvalues
cww   call elemev(s,nel,ndf,nst)
c


c.... compute residual force if necessary, p,s,ul are local 
      do 600 i = 1,nst
      do 600 j = 1,nst
600   p(i) = p(i) - s(i,j)*ul(j,1)


c.... transform to the global coordinate displacements
      call trans3b(s,p,ul,nst,ndf,1) ! s_G
31    call trans3b(s,p,ul,nst,ndf,3) ! P_G 

      return
      
c.... output forces  R = S = kv-f   
4     continue

c.... from loads S
      call pzero(sig,12)
      do i = 1,12 
        sig(i) = -p(i)
      enddo
c.... from local displacements k*v
      do 40 i = 1,2*ndf
        do 40 j = 1,2*ndf
40    sig(i) = sig(i) + s(i,j)*ul(j,1) 

c.... change to correct sign
      do i = 1,6
          sig(i)        = -sig(i)
      end do

      if(isw.eq.13) goto 13
      mct = mct - 3
      if (mct.le.0) then
         write(iow,2002) o,head
         if(ior.lt.0) write(*,2002) o,head
         mct = 50
      endif

c.... stresses with coordinates
c         write(iow,2003) n,ma,(xl(i,1),i=1,ndm),(sig(i),i=1,6),
c     +                        (xl(i,2),i=1,ndm),(sig(i),i=7,12)
c         if(ior.lt.0) write(*,2003) n,ma,(xl(i,1),i=1,ndm),
c     +                (sig(i),i=1,6),(xl(i,2),i=1,ndm),(sig(i),i=7,12)

c.... stresses without coordinates
         write(iow,2003) n,ma,(sig(i),i=1,6),(sig(i),i=7,12)
         if(ior.lt.0) write(*,2003) n,ma,
     +                (sig(i),i=1,6),(sig(i),i=7,12)
      return

c.... mass matrix 
5     continue
c.... mass matrix (lumped only)
      p(1) = d(15)*dl/2.d0
      p(2) = p(1)
      p(3) = p(1)
      p(ndf+1) = p(1)
      p(ndf+2) = p(1)
      p(ndf+3) = p(1)
cww   keine Transformation, da sonst neg. massen moeglich!!
cww   call trans3b(s,p,ul,nst,ndf,3) 
      return

c.... plot stress resultants
13    klayf = 1
      if(iplma(ma).eq.0)       return ! only if MATN
      if(nfp.lt.1.or.nfp.gt.7) return
      nd2 = ndf + nfp
      if(nfp.lt.7) then
        s1  = sig(nfp)
        s2  = sig(nd2)
      end if
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
          call plloco04(xl,ndm,dl)
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
          xlr(kk,1) = xl(kk,1) + dy(kk) *(ii-1)
          xlr(kk,2) = xl(kk,1) + dy(kk) * ii
        enddo
c.....  Randschnittgroessen
        dla = dl1 * (ii-1)
        dle = dl1 *  ii
        if(nfp.eq.1) then       ! N_x
          s1  = sig(1)  - qx * dla
          s2  = sig(1)  - qx * dle
        elseif(nfp.eq.2) then   ! Q_y
          s1  = sig(2)  - qy * dla
          s2  = sig(2)  - qy * dle
        elseif(nfp.eq.3) then   ! Q_z 
          s1  = sig(3)  - qz * dla
          s2  = sig(3)  - qz * dle
        elseif(nfp.eq.4) then   ! M_x
          s1 = sig(4)
          s2 = sig(4)
        elseif(nfp.eq.5) then   ! M_y
          s1  = sig(5) + sig(3) * dla - qz * dla*dla*0.5d0
          s2  = sig(5) + sig(3) * dle - qz * dle*dle*0.5d0
        else if(nfp.eq.6) then  ! M_z
          s1  = sig(6) - sig(2) * dla + qy * dla*dla*0.5d0
          s2  = sig(6) - sig(2) * dle + qy * dle*dle*0.5d0
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
          xdp(1) =  xlr(1,1) * (1.d0-xd/dl1) + xlr(1,2) * xd/dl1
          xdp(2) =  xlr(2,1) * (1.d0-xd/dl1) + xlr(2,2) * xd/dl1
          xdp(3) =  xlr(3,1) * (1.d0-xd/dl1) + xlr(3,2) * xd/dl1
c....     plotte dreieck 1
          sm = s1*0.5d0
          call pppcolf(sm)
          call pzero(xll,12)
          xll(1,1) = xlr(1,1)
          xll(2,1) = xlr(2,1)
          xll(3,1) = xlr(3,1)
          xll(1,2) = xdp(1)
          xll(2,2) = xdp(2)
          xll(3,2) = xdp(3)
          xll(1,3) = xlr(1,1) + dp1*s1*cfp
          xll(2,3) = xlr(2,1) + dp2*s1*cfp 
          xll(3,3) = xlr(3,1) + dp3*s1*cfp
          call plxtrn(xll,tra,vr,3,4)
          if(zoom(xll,3,3)) call plot9s(ixld,xll,ndm,3)
c....     plotte dreieck 2
          sm = s2*0.5d0
          call pppcolf(sm)
          call pzero(xll,12)
          xll(1,1) = xdp(1)
          xll(2,1) = xdp(2)
          xll(3,1) = xdp(3)
          xll(1,2) = xlr(1,2)
          xll(2,2) = xlr(2,2)
          xll(3,2) = xlr(3,2)
          xll(1,3) = xlr(1,2) + dp1*s2*cfp
          xll(2,3) = xlr(2,2) + dp2*s2*cfp 
          xll(3,3) = xlr(3,2) + dp3*s2*cfp
          call plxtrn(xll,tra,vr,3,4)
          if(zoom(xll,3,3)) call plot9s(ixld,xll,ndm,3)
          goto 135
        endif
c....   plotte trapez
130     xll(1,1) = xlr(1,1)
        xll(2,1) = xlr(2,1)
        xll(3,1) = xlr(3,1)
        xll(1,2) = xlr(1,2)
        xll(2,2) = xlr(2,2)
        xll(3,2) = xlr(3,2)
        xll(1,3) = xlr(1,2) + dp1*s2*cfp
        xll(2,3) = xlr(2,2) + dp2*s2*cfp 
        xll(3,3) = xlr(3,2) + dp3*s2*cfp
        xll(1,4) = xlr(1,1) + dp1*s1*cfp
        xll(2,4) = xlr(2,1) + dp2*s1*cfp 
        xll(3,4) = xlr(3,1) + dp3*s1*cfp
        call plxtrn(xll,tra,vr,3,4)
        if(zoom(xll,3,4)) call plot9s(ixl,xll,ndm,4)
135     continue
      endif
      return
c
c.... format statements
c.... with coordinates
c2002  format(a1,20a4,/,
c     *  13x,'3 - D beam: local(!)  s t r e s s  r e s u l a n t s',/,
c     *    'el mat',1x,'1-coor',1x,'2-coor',1x,'3-coor',2x,
c     *'N-x ',5x,'Q-y ',5x,'Q-z ',5x,'M-x ',5x,'M-y ',5x,'M-z ',/)
c2003  format(i2,i2,3f7.2,6e9.3/4x,3f7.2,6e9.3/)
c.... with coordinates
2002  format(a1,20a4,/,
     *  13x,'3 - D beam: local(!)  s t r e s s  r e s u l a n t s',/,
     *    ' el mat',4x,'N_x ',5x,'Q_y/R_y ',4x,'Q_z/R_z ',
     *              8x,'M_x ',8x,'M_y ',8x,'M_z ')
2003  format(i4,i3,6(1x,g11.5)/7x,6(1x,g11.5))
      end
c
      subroutine trans3b(s,p,ul,nst,ndf,isw)
c----------------------------------------------------------------------+
c.... isw: 1  transform matrix s(nst,nst)   Sg  = Tt * Sl * T
c          2  transform vector ul(ndf,1)    Ul = T  * Ug
c          3  transform vector P(nst)       Pg = Tt * Pl
c          4  transform vector ul(ndf,1)    Ug = Tt * Ul
c.... extended Version of fram3d/pcelm8b ww 8/91
c----------------------------------------------------------------------+
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
      subroutine rotate04(alpha)
c----------------------------------------------------------------------+
c.... drehe KOS  e_yR =  c * e_yL  +  s * e_zL
c....            e_zR = -s * e_yL  +  c * e_zL
c----------------------------------------------------------------------+
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
      subroutine plloco04(xl,ndm,dl)
c----------------------------------------------------------------------+
c.... Plot local basis
c----------------------------------------------------------------------+
      implicit double precision (a-h,o-z)
      dimension x0(3),xl(ndm,*),tr(3,3)
      common /tran3d/ dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3
!$OMP THREADPRIVATE (/tran3d/)  
c...  plot at midpoint of element (global coordinates!)
      x0(1) = 0.5d0*( xl(1,2) + xl(1,1) )
      x0(2) = 0.5d0*( xl(2,2) + xl(2,1) )
      x0(3) = 0.5d0*( xl(3,2) + xl(3,1) )
c.... length
      xm = 0.3d0 * dl
c.... transformation matrix 
      tr(1,1) = dx1
      tr(2,1) = dy1
      tr(3,1) = dz1

      tr(1,2) = dx2
      tr(2,2) = dy2
      tr(3,2) = dz2

      tr(1,3) = dx3
      tr(2,3) = dy3
      tr(3,3) = dz3

c...  plot axis 
      call pppcol(2)
      call pltaxs04(tr,x0,xm)
      return
      end
c
      subroutine pltaxs04(tr,x0,xm)
c----------------------------------------------------------------------
c.... draw vectors for axes, including rot-macro 
c.... t_1=black
c.... t_2=red
c.... t_3=blue
c----------------------------------------------------------------------
      USE pdata2
      USE pltran
      implicit double precision (a-h,o-z)
      logical zoom
      dimension xx(3,5),x0(3),x0z(3),tr(3,3)
c.... test if vectors are in plot-region
      x0z(1) = x0(1)
      x0z(2) = x0(2)
      x0z(3) = x0(3)
      call plxtrn(x0z,tra,vr,3,1)
      if(zoom(x0z,3,1)) then
c....   perspective projecton of axes
        do 120 m = 1,3  
	        call pppcol(m)
          call pzero(xx,15)
          do 100 n = 1,3  
            fac1    = tr(m,n)*xm
            xx(n,1) = x0(n)
            xx(n,2) = xx(n,1) + fac1
100       xx(n,5) = xx(n,2)
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
          call plxtrn(xx,tra,vr,3,5)
c....     plot the vector
          call plotl(xx(1,1),xx(2,1),xx(3,1),3)
          call plotl(xx(1,2),xx(2,2),xx(3,2),2)
          call plotl(xx(1,2),xx(2,2),xx(3,2),ipgl)
          call plotl(xx(1,3),xx(2,3),xx(3,3),2)
          call plotl(xx(1,4),xx(2,4),xx(3,4),2)
          call plotl(xx(1,2),xx(2,2),xx(3,2),2)
          if(ipgl.eq.1) call clpan
          call plotl(xx(1,2),xx(2,2),xx(3,2),3)
          call plabl(m)
120     continue
      end if
      return
      end
c      
      subroutine qload04(qx,qy,qz,d,q,propq,prop,numel,n,mqloa,isw)
c----------------------------------------------------------
c.... add loads local/global and from macro qloa
c.... P_L = P_L + T*P_G
c     P_L=d(7-9)+q(1-3) P_G=d(10-12)+q(4-6) 
c----------------------------------------------------------
      implicit double precision (a-h,o-z)
      common/tran3d/ tr(3,3)
!$OMP THREADPRIVATE (/tran3d/)  
      dimension ql(3),q(numel,10),d(*)
      call pzero(ql,3) 
      if(isw.eq.4.or.isw.eq.13) then 
c....   from mate
        do i = 1,3 
          ql(i) = d(6+i)*prop
          do j = 1,3  
            ql(i) = ql(i) + tr(i,j) * d(9+j)*prop
          end do
        end do  
      else if(isw.eq.22) then 
c....   from qloa
        if(mqloa.ne.1) then 
          do i = 1,3 
            ql(i) = ql(i) + q(n,i)*propq
            do j = 1,3  
              ql(i) = ql(i) + tr(i,j) * q(n,3+j)*propq
            end do
          end do  
        end if
      else
c....   from mate
        do i = 1,3 
          ql(i) = d(6+i)*prop
          do j = 1,3  
            ql(i) = ql(i) + tr(i,j) * d(9+j)*prop
          end do
        end do  
  
      end if 

      qx=ql(1)
      qy=ql(2)
      qz=ql(3)
      return
      end
