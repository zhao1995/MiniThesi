      subroutine elmt34(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c-----------------------------------------------------------------------
c..... flat shell element Elsmp3_7
c      * geometrically nonlinear (moderate rotations)
c      * elasto-plastic with isotropic linear hardeneing
c      * thickness integration: n layers with 2 point gauss quadrature
c      *   4-node shell element ->    TRAFO=1 
c      *   4-node plate element ->    TRAFO=0
c      * 8/9-node plate element ->    TRAFO has no influence 
c      * actual: Trafo=1 
c----------------------------------------------------------------------
c..... Remarks:
c      Aenderungen bzw. offen:
c      in ylds34: - >>> Abbruch fuer f = 0
c      Abbruch fuer ddlam < tol
c        --> keine quadr. Konvergenz
c      Abbruch fuer ddlam < tol*dlam(RLT)
c        --> keine Loesung fuer sig in jeweils 1.Iteration
c     - >>> Zuwachs  ddlam von RLT, wichtig bei starken Aend.
c     - Aenderung 1.e-8 fuer dlam um Sig auf F zu setzen.
c     - Anzahl der lokalen Iterationen  50
c      open problems:
c     -  plasticity only for sig, tau not included!
c
c     19.9.93 ELMP3_7
c
c     04/1994  WW IBNM UH same as Elsmp3_6 but rotations mod!
c                         output for resultants/layer stresses
c
c     01/2014 WW IBS KIT  Update auf aktuelle FEAP-Version
c           - Konvergenz ist quadratisch
c           - Rückrechnung korrigiert, Fehler bei Schichtplot in nn-Zählung korrigiert
c           - Trafo wahlweise
c
c----------------------------------------------------------------------
c
c....  Material parameters
c   1   E   = Young's modulus
c   2   nu  = Poisson's ratio
c   3   h   = thickness of shell
c   4   qn  = transverse load (local z-direction)
c   5   gam = spec. weight, M = gam*V, Gam = rho * g
c   6   nlay= no. of layers over thickness
c   7   initial yield stress Y-0
c   8   cp  = Hardening parameter
c  13   ilin= 0/1 lin/nonlinear
c  14   qz  = load in global z-direction
c   -   [stiffness parameter 6dof]
c-------------------------------------
c   9  nh
c  10  e/(1-nu*nu)
c  11  e/(1-nu*nu) * nu
c  12  G
c-----------------------------------------------------------------------
c.... definition of array h  for each gauss point and(!) each layer
c   1   sig_x
c   2   sig_y
c   3   sig_xy
c   4   ep_x
c   5   ep_y
c   6   ep_xy
c   7   epq    = effective plastic strain
c   8   sigv/y_0
c   9   dlam0    = plastic multiplier of last iterate
c-----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE eldata
      USE evdata
      USE hdata
      USE iofile
      USE plslay
      USE prisdat
      USE prlod
      USE strnam

      implicit double precision (a-h,o-z)
      dimension d(*),ul(ndf,*),xl(ndm,*),ix(*),tl(*),s(nst,*),p(*),
     1    pgb(4),wgb(4),pgs(4),wgs(4),pgz(4),wgz(4),gp(2),
     2    deriv(2,9),shape(3,9),sig(5),sigp(8),eps(3),
     3    dme(3,3),ds(2,2),bmi(3,5),bmj(3,5),bsi(2,3),bsj(2,3),
     4    gri(10),emi(3),ebi(3),es(2),grn(10),emn(3),ebn(3),
     5    yl(3,9),dul(6,9),vl(6,9),dvl(6,9),tra(3,3),ql(3),
     6    bmli(3,5),bmlj(3,5)
      dimension h1(*),h2(*),h3(*) 

      save ngb,pgb,wgb,ngs,pgs,wgs,ngz,pgz,wgz

      itrafo=1

c.... go to correct array processor
      go to(1,2,3,4,5,3,2,4,2,2,3,2,2) isw
      return 
c.... input material properties
1     if(ior.lt.0) write(*,1001)
1001  format(' Input: E nu h q gamma nlay y-0 Cp lin(0/1) qz)',$)
      call dinput(d,10)
      d(13) = d(9)
      d(14) = d(10)
                   write(iow,1002) (d(i),i=1,8),d(13),d(14)
      if(ior.lt.0) write(  *,1002) (d(i),i=1,8),d(13),d(14)
1002  format(5x,'ELMT34 Materialdata:',/,
     + 5x,'elastic modulus............',g15.4,/,
     + 5x,'poissons ratio................',g12.4,/,
     + 5x,'thickness(h = 0 variable).....',g12.4,/,
     + 5x,'element load q_n(transverse)..',g12.4,/,
     + 5x,'spec. weight (gamma)..........',g12.4,/,
     + 5x,'no of layers in z-dir(def.=3).',g12.4,/,
     + 5x,'in.Yield stress Y0............',g12.4,/,
     + 5x,'hardening parameter Cp = Hp...',g12.4,/,
     + 5x,'linear = 0,  nonlinear = 1....',g12.4,/,
     + 5x,'element load q_z..............',g12.4)
c.... set integration points
      ngb=2
c     nel ist noch nicht gegeben, nen richtig? bei mischung versch. elemente?
      if(nen.gt.4) ngb=3
      ngs = ngb-1
      ngz = 2
      call gauss34(ngb,pgb,wgb)
      call gauss34(ngs,pgs,wgs)
      call gauss34(ngz,pgz,wgz)
      nlay = d(6)
c.... define number of points for h-array, 9=sig(3),ep(3),epq,yld,dlam0
      d(9) = 9
c.... total number for h-array
      nh1 = d(9)*ngb*ngb*ngz*nlay
c.... elasticity matrix
      d(10)  = d(1)/(1.-d(2)*d(2))
      d(11)  = d(10)*d(2)
      d(12)  = d(10)*(1.-d(2))*0.5
c.... description of stresses not used  
      do i=10,26
       strsus(i) = '               '                
      end do  
c
2     return
3     nh   = d(9)
      nlay = d(6)
      lin  = d(13)
c.... displacements at t_n  = u_i - du_i
      do 350 i = 1,nel
        do 350 j = 1,6
          dul(j,i) = ul(j,i) - ul(j,i+nen)
350   continue
c.... transform coordinates and displacements
      call trans34(xl,yl,ul,vl,dul,dvl,tra,ndm,ndf,nel,itrafo)
c.....Membrane and bending terms --> ngb
      nn   = 1
      do 300 igb = 1,ngb
      do 300 jgb = 1,ngb
       xsi = pgb(igb)
       eta = pgb(jgb)
c..... shape functions,area etc.
       call sfr34(deriv,eta,xsi,nel,shape)
       call jaco34(shape,deriv,djacb,yl,n,nel,ndm,iow)
       da = djacb*wgb(igb)*wgb(jgb)
c..... displacement gradients at  t_i
       call grad34(shape,gri,vl,ndf,nel)
c..... total shell strains  e_m e_b  at t_i
       call stra34(shape,vl,gri,nel,ndf,emi,ebi,es,lin,1,0)
c..... displacement gradients at  t_n
       call grad34(shape,grn,dvl,ndf,nel)
c..... total shell strains  e_m e_b  at t_n
       call stra34(shape,dvl,grn,nel,ndf,emn,ebn,es,lin,1,0)
       wkx=0.0
       wky=0.0
       if(lin.eq.1) then
         wkx  =  gri(3)
         wky  =  gri(8)
       end if
c..... loop over all layers(integration through the thickness)
       dh  = d(3)/nlay
       dh5 = 0.5*dh
       h5  = 0.5*d(3)
       do 361 ilay = 1,nlay
        dz = (ilay-1)*dh-h5
        do 360 igz = 1,ngz
         zeta = pgz(igz)
         wz   = wgz(igz)*dh5
         zs   = dz + (1.+zeta)*dh5
         dv = da * wz
c.....   incremental normal strains at zs (E = e_m + zs*e_b)
         eps(1) = (emi(1)-emn(1)) + zs*(ebi(1)-ebn(1))
         eps(2) = (emi(2)-emn(2)) + zs*(ebi(2)-ebn(2))
         eps(3) = (emi(3)-emn(3)) + zs*(ebi(3)-ebn(3))
c.....   stresses and material matrix at integration point
         call str34(h1(nn),h2(nn),nh,eps,dme,d,sig,sigv)
c.....   calculate stresses at zs
c        call stre34(sig(1),dme,eps,sig(4),ds,es,1,0)
c....    nonlinear stiffness matrix
         do 303 ino =1,nel
          call bmat34(bmi,bsi,shape,zs,wkx,wky,lin,ino,nel,1,0)
c.....    residual
          call mul34(bmi,sig(1),dv,p,ndf,nst,ino,5,3,1)
          if(d(4).eq.0.0.and.d(14).eq.0.0) goto 353
c.....    loads in local directions
c.....    external constant load qz in z-direction
          qz = d(14)
          do 352 i = 1,3
352        ql(i) = tra(i,3)*qz
c.....    external constant load qn transversal to element
          qn = d(4)
          ql(3) = ql(3)+qn
c.....    add loads to load vector
          if(ilay.eq.1.and.igz.eq.1) then
           iq = (ino-1)*ndf
           do 354 i = 1,3
354        p(iq+i) = p(iq+i) + shape(3,ino)*ql(i)*da*prop
          end if
353       continue
          if(isw.eq.6) goto 303
          do 304 jno = ino,nel
           call bmat34(bmj,bsj,shape,zs,wkx,wky,lin,jno,nel,1,0)
c.....     bmi T * dme * bmj
           call sub34(bmi,bmj,dv,dme,s,ndf,nst,ino,jno,5,3,5,1,1)
c.....     initial stress matrix
           if(lin.eq.1) call ksig34(ino,jno,s,shape,sig,nel,nst,ndf,dv)
304       continue      ! JNO
303      continue       ! INO 
         nn = nn + nh
360     continue        ! IGZ
361    continue         ! ILAY
300   continue          ! IGB,JGB
c.....Shear terms   --> ngs
      do 320 igs = 1,ngs
      do 320 jgs = 1,ngs
       xsi = pgs(igs)
       eta = pgs(jgs)
c..... shape functions,area etc.
       call sfr34(deriv,eta,xsi,nel,shape)
       call jaco34(shape,deriv,djacb,yl,n,nel,ndm,iow)
       da = djacb*wgs(igs)*wgs(jgs)
c..... displacement gradients
       call grad34(shape,gri,vl,ndf,nel)
c..... shell strains e_s
       call stra34(shape,vl,gri,nel,ndf,emi,ebi,es,lin,0,1)
c..... loop over all layers(integration through the thickness)
       dh  = d(3)/nlay
       dh5 = 0.5*dh
       h5  = 0.5*d(3)
       do 371 ilay = 1,nlay
        dz = (ilay-1)*dh-h5
        do 370 igz = 1,ngz
         zeta = pgz(igz)
         wz   = wgz(igz)*dh5
         zs   = dz + (1.+zeta)*dh5
         dv   = da * wz
c.....   strains e_s are constant for z --> no calculation
c.....   material matrices
         call mod34(dme,ds,d,0,1)
c.....   calculate stresses at zs
         call stre34(sig(1),dme,eps,sig(4),ds,es,0,1)
c....    nonlinear stiffness matrix
         do 323 ino = 1,nel
          call bmat34(bmi,bsi,shape,zs,wkx,wky,lin,ino,nel,0,1)
c.....    residual
          call mul34(bsi,sig(4),dv,p,ndf,nst,ino,3,2,3)
          if(isw.eq.6) goto 323
          do 324 jno = ino,nel
           call bmat34(bmj,bsj,shape,zs,wkx,wky,lin,jno,nel,0,1)
324       call sub34(bsi,bsj,dv,ds,s,ndf,nst,ino,jno,3,2,3,3,3)
323      continue
370     continue
371    continue
320   continue
c.....lower part of K_T
      if(isw.eq.6) goto 328
      do 325 i = 1,nst
       do 325 j = 1,i
325   s(i,j) =s(j,i)
c
328    continue
c.... transform s and p to global coordinates
      call trans134(s,p,tra,nst,ndf,nel,itrafo)
      return
c.... print and plot stresses --> ngb
4     nh   = d(9)
      nlay = d(6)
      lin  = d(13)
c.... displacements at t_n  = u_i - du_i
      do 450 i = 1,nel
       do 450 j = 1,6
        dul(j,i) = ul(j,i) - ul(j,i+nen)
450   continue
c.... transform coordinates and displacements
      call trans34(xl,yl,ul,vl,dul,dvl,tra,ndm,ndf,nel,itrafo)
      nn = 1
      do 40 igb = 1,ngb
      do 40 jgb = 1,ngb
       xsi = pgb(igb)
       eta = pgb(jgb)
c..... shape functions,area etc.
       call sfr34(deriv,eta,xsi,nel,shape)
       call jaco34(shape,deriv,djacb,yl,n,nel,ndm,iow)
       da = djacb*wgb(igb)*wgb(jgb)
c..... coordinates of Gauss points
       call pzero(gp,2)
       do idm = 1,2
        do ino = 1,nel
         gp(idm) = gp(idm) + xl(idm,ino) * shape(3,ino)
        end do
       end do
c..... displacement gradients at t_i
       call grad34(shape,gri,vl,ndf,nel)
c..... shell strains e_m, e_b, e_s at t_i
       call stra34(shape,vl,gri,nel,ndf,emi,ebi,es,lin,1,1)
c..... displacement gradients at t_n
       call grad34(shape,grn,dvl,ndf,nel)
c..... shell strains e_m, e_b   at t_n
       call stra34(shape,dvl,grn,nel,ndf,emn,ebn,es,lin,1,0)
c..... loop over all layers(integration through the thickness)
       call pzero(sigp,8)
       dh  = d(3)/nlay
       dh5 = 0.5*dh
       h5  = 0.5*d(3)
       nzlay = 0
       do 48 ilay = 1,nlay
        dz = (ilay-1)*dh-h5
        do 43 igz = 1,ngz
         nzlay = nzlay+1
         zeta = pgz(igz)
         wz   = wgz(igz)*dh5
         zs   = dz + (1.+zeta)*dh5
c.....   incremental normal strains at zs (E = e_m + zs*e_b)
         eps(1) = (emi(1)-emn(1)) + zs*(ebi(1)-ebn(1))
         eps(2) = (emi(2)-emn(2)) + zs*(ebi(2)-ebn(2))
         eps(3) = (emi(3)-emn(3)) + zs*(ebi(3)-ebn(3))
c.....   stresses and material matrix at integration point
         call str34(h1(nn),h2(nn),nh,eps,dme,d,sig,sigv)
         if(nzlay.eq.klay) then
c...      plot/print stresses for layer  
c         call getstr34(h1(nn),sigp,8)
          call getstr34(h2(nn),sigp,8)
          if(isw.eq.4) then
           mct = mct - 1
           if(mct.gt.0) go to 45
                        write(iow,2007) o,head,klay
           if(ior.lt.0) write(*  ,2007) o,head,klay
2007       format(a1,20a4,/,2x,'E L E M E N T   S T R E S S E S',
     *         ' for layer',i5,/,
     1     2x,'El',1x,'Mat',1x,'1-Coord',1x,'2-Coord',
     2     2X,'Sigma_x ',3X,'Sigma_y ',3X,'Sigma_xy',3X,'Epspl_q ',/,
     3     26X,'Epspl_x ',3X,'Epspl_y ',3X,'Epspl_xy',3X,'Sigv/Y0 ',/)
           mct = 50
45                      write(iow,2006) n,ma,(gp(i),i=1,2),
     1                 (sigp(i),i=1,3),sigp(7),(sigp(i),i=4,6),sigp(8)
           if(ior.lt.0) write(*  ,2006) n,ma,(gp(i),i=1,2),
     1                 (sigp(i),i=1,3),sigp(7),(sigp(i),i=4,6),sigp(8)
          else if(isw.eq.8) then
c....      calculate nodal stresses
           istv = -9
           if(igb.eq.1.and.jgb.eq.1) then
            strsus( 1) = '  SIGMA_xx     '
            strsus( 2) = '  SIGMA_yy     '
            strsus( 3) = '  SIGMA_xy     '
            strsus( 4) = '  EPS_PL_xx    '
            strsus( 5) = '  EPS_PL_yy    '
            strsus( 6) = '  EPS_PL_xy    '
            strsus( 7) = '  EPS_EFF_PL   '
            strsus( 8) = '  SIGV / Y_0   '
            strsus( 9) = '  SIGV         '
            nptyp   =1 
            nprip(1)=1
            nprip(2)=3
            nprip(3)=2
            nprip(4)=4
            nprip(5)=5
            nprip(6)=6
            nprip(7)=7
            nprip(8)=8
           end if
           maxlay=ngz*nlay+1
           if(klay.le.maxlay)
     1     call plots34(ix,strea,strea(1+numnp),shape,sigp,sigv,
     2                  nel,da,numnp)
          end if
         end if
c.....   material matrix only for shear
         call mod34(dme,ds,d,0,1)
c.....   calculate only shear stresses
         call stre34(sig(1),dme,eps,sig(4),ds,es,0,1)
c.....   calculate stress resultants
         sigp(1) = sigp(1)+sig(1)*wz
         sigp(2) = sigp(2)+sig(2)*wz
         sigp(3) = sigp(3)+sig(3)*wz
         sigp(4) = sigp(4)+sig(1)*wz*zs
         sigp(5) = sigp(5)+sig(2)*wz*zs
         sigp(6) = sigp(6)+sig(3)*wz*zs
         sigp(7) = sigp(7)+sig(4)*wz
         sigp(8) = sigp(8)+sig(5)*wz
         nn = nn + nh
43      continue ! IGZ 
48     continue  ! ILAY
       if(klay.eq.0) then
c...    plot/print stress resultants  
        if(isw.eq.4) then
c....    print stress resultants  N M and Q
         mct = mct - 1
         if(mct.gt.0) go to 47
                      write(iow,2005) o,head
         if(ior.lt.0) write(*  ,2005) o,head
2005     format(a1,20a4,/,2x,'E L E M E N T   Stress Resultants',/,
     1        2x,'El',1x,'Mat',1x,'1-Coord',1x,'2-Coord',
     2    2X,'   N_x  ',3X,'   N_y  ',3X,'   N_xy ',3X,'   Q_y  ',/,
     3   26X,'   M_x  ',3X,'   M_y  ',3X,'   M_xy ',3X,'   Q_z  ',/)
         mct = 50
47                    write(iow,2006)n,ma,(gp(i),i=1,2),
     1               (sigp(i),i=1,3),sigp(7),(sigp(i),i=4,6),sigp(8)
         if(ior.lt.0) write(*  ,2006)n,ma,(gp(i),i=1,2),
     1               (sigp(i),i=1,3),sigp(7),(sigp(i),i=4,6),sigp(8)
2006     format(1x,i3,i4,2f8.3,1x,4e11.4,/,25x,4e11.4)
        else if(isw.eq.8) then
c....    calculate nodal stresses
         istv = -9
         sigv=0.d0
         if(igb.eq.1.and.jgb.eq.1) then
          strsus( 1) = '  N_xx         '
          strsus( 2) = '  N_yy         '
          strsus( 3) = '  N_xy         '
          strsus( 4) = '  M_xx         '
          strsus( 5) = '  M_yy         '
          strsus( 6) = '  M_xy         '
          strsus( 7) = '  Q_xz         '
          strsus( 8) = '  Q_yz         '
          strsus( 9) = '               '
          nptyp   =5 
          nprip(1)=1
          nprip(2)=3
          nprip(3)=2
          nprip(4)=4
          nprip(5)=6
          nprip(6)=5
          nprip(7)=7
          nprip(8)=8
         end if
         call plots34(ix,strea,strea(1+numnp),shape,sigp,sigv,
     1                nel,da,numnp)
        end if
       end if
40    continue

      return
5     continue
      if(imtyp.eq.1) then
c....  lumped mass matrix
       gam  = d(5)
c....  8-node element produces negative mass terms
       if(nel.eq.8) stop 'element produces negative mass terms'
c....  loop over Gauss points
       do 581 igb = 1,ngb
       do 581 jgb = 1,ngb
        xsi = pgb(igb)
        eta = pgb(jgb)
c....   shape functions
        call sfr34(deriv,eta,xsi,nel,shape)
c....   Jacobian,Inverse and determinant
        call jaco34(shape,deriv,djacb,yl,n,nel,ndm,iow)
c....   area of element
        da = djacb*wgb(igb)*wgb(jgb)
c....   lumped mass matrix  only for u,v,w!
        do 582 ino = 1, nel
         do 582 idf = 1, 3
          iq = (ino-1)*ndf+idf
 582      p(iq) = p(iq) + shape(3,ino)*gam*da*d(3)
 581   continue
c....  transform s and p to global coordinates
c      call trans134(s,p,tra,nst,ndf,nel,itrafo)  pruefen!!
      else if(imtyp.eq.2.or.imtyp.eq.3) then
c....  matrix k_u + k_sigma
       nh   = d(9)
       nlay = d(6)
       lin  = d(13)
c....  displacements at t_n  = u_i - du_i
       do 550 i = 1,nel
        do 550 j = 1,6
         dul(j,i) = ul(j,i) - ul(j,i+nen)
550      continue
c....  transform coordinates and displacements
       call trans34(xl,yl,ul,vl,dul,dvl,tra,ndm,ndf,nel,itrafo)
c..... Membrane       --> ngb
       nn   = 1
       do 500 igb = 1,ngb
       do 500 jgb = 1,ngb
        xsi = pgb(igb)
        eta = pgb(jgb)
c.....  shape functions,area etc.
        call sfr34(deriv,eta,xsi,nel,shape)
        call jaco34(shape,deriv,djacb,yl,n,nel,ndm,iow)
        da = djacb*wgb(igb)*wgb(jgb)
c.....  displacement gradients at  t_i
        call grad34(shape,gri,vl,ndf,nel)
c.....  total shell strains  e_m e_b  at t_i
        call stra34(shape,vl,gri,nel,ndf,emi,ebi,es,lin,1,0)
c.....  displacement gradients at  t_n
        call grad34(shape,grn,dvl,ndf,nel)
c.....  total shell strains  e_m e_b  at t_n
        call stra34(shape,dvl,grn,nel,ndf,emn,ebn,es,lin,1,0)
        wkx=0.0
        wky=0.0
        if(lin.eq.1) then
         wkx  =  gri(3)
         wky  =  gri(8)
        end if
c.....  loop over all layers(integration through the thickness)
        dh  = d(3)/nlay
        dh5 = 0.5*dh
        h5  = 0.5*d(3)
        do 561 ilay = 1,nlay
         dz = (ilay-1)*dh-h5
         do 560 igz = 1,ngz
          zeta = pgz(igz)
          wz   = wgz(igz)*dh5
          zs   = dz + (1.+zeta)*dh5
c......   -dv
          dv = -da * wz
c.....    incremental normal strains at zs (E = e_m + zs*e_b)
          eps(1) = (emi(1)-emn(1)) + zs*(ebi(1)-ebn(1))
          eps(2) = (emi(2)-emn(2)) + zs*(ebi(2)-ebn(2))
          eps(3) = (emi(3)-emn(3)) + zs*(ebi(3)-ebn(3))
c.....    stresses and material matrix at integration point
          call str34(h1(nn),h2(nn),nh,eps,dme,d,sig,sigv)
          do 503 ino =1,nel
           call bmatl34(bmli,shape,zs,ino)
           call bmatnl34(bmi,shape,wkx,wky,ino)
           do 504 jno = ino,nel
            call bmatl34(bmlj,shape,zs,jno)
            call bmatnl34(bmj,shape,wkx,wky,jno)
c.....      bmli T * dme * bmj + bmi T * dme * bmlj + bmi T * dme * bmj
            if(imtyp.eq.3) then
             call sub34(bmli,bmj,dv,dme,s,ndf,nst,ino,jno,5,3,5,1,1)
             call sub34(bmi,bmlj,dv,dme,s,ndf,nst,ino,jno,5,3,5,1,1)
             call sub34(bmi,bmj,dv,dme,s,ndf,nst,ino,jno,5,3,5,1,1)
            end if
c.....      initial stress matrix
            call ksig34(ino,jno,s,shape,sig,nel,nst,ndf,dv)
504        continue
503       continue
          nn = nn + nh
560      continue
561     continue
500    continue
c..... lower part of K_u + K_sigma
       do 525 i = 1,nst
        do 525 j = 1,i
525      s(i,j) =s(j,i)
c....  transform s and p to global coordinates
       call trans134(s,p,tra,nst,ndf,nel,itrafo)
      else if(imtyp.eq.4) then
c....  linear stiffness matrix (not implemented)
      end if
      return
      end
c
      subroutine bmat34(bm,bs,shape,zs,wkx,wky,nl,kk,nel,ifm,ifs)
c----------------------------------------------------------------------
c.....B-matrices for flat shell
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  bm(3,5),bs(2,3),shape(3,*)
      dndx = shape(1,kk)
      dndy = shape(2,kk)
      dn   = shape(3,kk)
c.....form bm  = bm0 + zs*bb
      if(ifm.eq.1) then
       call pzero(bm,15)
c....  membrane
       bm(1,1) = dndx
       bm(2,2) = dndy
       bm(3,1) = dndy
       bm(3,2) = dndx
c....  nonlinear
       if(nl.eq.1) then
        bm(1,3) = wkx*dndx
        bm(2,3) = wky*dndy
        bm(3,3) = wky*dndx + wkx*dndy
       end if
c....  bending
       bm(1,5) =  dndx*zs
       bm(2,4) = -dndy*zs
       bm(3,4) = -dndx*zs
       bm(3,5) =  dndy*zs
      end if
c.....form bs
      if(ifs.eq.1) then
       call pzero(bs,6)
       bs(1,1) =  dndx
       bs(1,3) =  dn
       bs(2,1) =  dndy
       bs(2,2) = -dn
      endif
      return
      end
c
      subroutine bmatl34(bm,shape,zs,kk)
c----------------------------------------------------------------------
c.....linear B-matrices for flat shell
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  bm(3,5),shape(3,*)
      dndx = shape(1,kk)
      dndy = shape(2,kk)
      dn   = shape(3,kk)
c.....form bm  = bm0 + zs*bb
      call pzero(bm,15)
c.... membrane
      bm(1,1) = dndx
      bm(2,2) = dndy
      bm(3,1) = dndy
      bm(3,2) = dndx
c.... bending
      bm(1,5) =  dndx*zs
      bm(2,4) = -dndy*zs
      bm(3,4) = -dndx*zs
      bm(3,5) =  dndy*zs
      return
      end
c
      subroutine bmatnl34(bm,shape,wkx,wky,kk)
c----------------------------------------------------------------------
c.....nonlinear part of B-matrices for flat shell
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  bm(3,5),shape(3,*)
      dndx = shape(1,kk)
      dndy = shape(2,kk)
      dn   = shape(3,kk)
      call pzero(bm,15)
      bm(1,3) = wkx*dndx
      bm(2,3) = wky*dndy
      bm(3,3) = wky*dndx + wkx*dndy
      return
      end
c
      subroutine ksig34(ino,jno,s,shape,sig,nel,nst,ndf,da)
c----------------------------------------------------------------------
c.... initial stress matrix K_s = g T * sig * g
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension shape(3,*),s(nst,*),sig(5)
      dndxi = shape(1,ino)
      dndyi = shape(2,ino)
      dndxj = shape(1,jno)
      dndyj = shape(2,jno)
      aksig = dndxi * sig(1) * dndxj + dndxi * sig(3) * dndyj +
     +        dndyi * sig(2) * dndyj + dndyi * sig(3) * dndxj
      aksig = aksig*da
c.....store in K_T at  3,3
      irow = (ino-1) * ndf + 3
      jcol = (jno-1) * ndf + 3
      s(irow,jcol) = s(irow,jcol) + aksig
      return
      end
c
      subroutine grad34(shape,gr,ul,ndf,nel)
c----------------------------------------------------------------------
c.....displacement gradients
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension shape(3,*),gr(10),ul(ndf,*)
      call pzero(gr,10)
      do 20 ino = 1,nel
       dnidx = shape(1,ino)
       dnidy = shape(2,ino)
       do 20 idf = 1,5
       idfn = 5 + idf
       const = ul(idf,ino)
       gr(idf)  = gr(idf) + dnidx * const
   20      gr(idfn) = gr(idfn)+ dnidy * const
      return
      end
c
      subroutine stra34(shape,ul,gr,nel,ndf,em,eb,es,nl,ifm,ifs)
c----------------------------------------------------------------------
c.....shell strains for flat shell
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  gr(10),shape(3,*),ul(ndf,*),em(3),eb(3),es(2)
      if(ifs.eq.1) then
c..... rotation at Gauss point
       xzrot = 0.0
       yzrot = 0.0
       nposn = 4
       nposn1= 5
       do 30 ino = 1,nel
        xzrot = xzrot+shape(3,ino)*ul(nposn, ino)
   30   yzrot = yzrot+shape(3,ino)*ul(nposn1,ino)
c..... shear strains
       es(1)=gr(3)+yzrot
        es(2)=gr(8)-xzrot
      end if
      if(ifm.eq.1) then
c..... membrane strains
       em(1) = gr(1)
       em(2) = gr(7)
       em(3) = gr(2)+ gr(6)
       if(nl.eq.1) then
        em(1) = em(1) + gr(3)*gr(3)/2.
        em(2) = em(2) + gr(8)*gr(8)/2.
        em(3) = em(3) + gr(3)*gr(8)
       end if
c..... bending strains
       eb(1) =  gr(5)
       eb(2) = -gr(9)
       eb(3) = -gr(4)+gr(10)
      end if
      return
      end
c
      subroutine jaco34(shape,deriv,djacb,xl,n,nel,ndm,na)
c----------------------------------------------------------------------
c.....Jacobi-Matrix and  cartesian derivatives of shape functions
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension xl(ndm,*),shape(3,*),deriv(2,*),xjaci(2,2),xjacm(2,2)
c.....Jacobi-matrix xjacm
      call pzero(xjacm,4)
      do   4 idm = 1,2
       do  4 jdm = 1,2
        do 4 ino = 1,nel
    4   xjacm(idm,jdm)=xjacm(idm,jdm)+deriv(idm,ino)*xl(jdm,ino)
c.....Determinant  and Inverse of Jacobi-Matrix
      djacb=xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1)
      if(djacb) 6,6,8
    6 write(*,600) n
      stop
    8 xjaci(1,1)= xjacm(2,2)/djacb
      xjaci(2,2)= xjacm(1,1)/djacb
      xjaci(1,2)=-xjacm(1,2)/djacb
      xjaci(2,1)=-xjacm(2,1)/djacb
c.....cartesian derivatives
      do 10 idm = 1,2
       do 10 ino = 1,nel
        shape(idm,ino) = 0.0
        do 10 jdm = 1,2
   10     shape(idm,ino)=shape(idm,ino)+xjaci(idm,jdm)*deriv(jdm,ino)
  600 format(' program halted in subroutine jaco34',/,
     1       ' zero or negative area for element number ',i5)
      return
      end
c
      subroutine mul34(bimat,bjvec,da,p,ndf,nst,ino,ncoli,nrowj,irow)
c----------------------------------------------------------------------
c.....bi T * bj  * det j * weig(xsi) * weig(eta)
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension bimat(nrowj,ncoli),bjvec(nrowj),p(nst),sblod(5)
      ifrow = irow - 1
c.....bi T * bjv
      do 10 i = 1,ncoli
       sblod(i) = 0.0
       do 10 j = 1,nrowj
   10   sblod(i) = sblod(i) + bimat(j,i) * bjvec(j)
c.....store sblod in load vector
      ifrow = (ino-1)*ndf + ifrow
      do 30 i = 1,ncoli
       irsub = ifrow + i
   30 p(irsub) = p(irsub) - sblod(i)*da
      return
      end
c
      subroutine mod34(dme,ds,d,ifm,ifs)
c----------------------------------------------------------------------
c.....Elasticity matrices for flat shells
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension dme(3,3),ds(2,2),d(*)
      e     = d(1)
      xnu   = d(2)
      cappa = 5.0d0/6.0d0
      c1    = e/(1.-xnu*xnu)
      c2    = e/(1.-xnu*xnu)*xnu
      c3    = 0.5*c1*(1.-xnu)
      c4    = cappa*0.5*e/(1.+xnu)
c.....normal   part
      if(ifm.eq.1) then
       call pzero(dme,9)
       dme(1,1) = c1
       dme(2,2) = c1
       dme(1,2) = c2
       dme(2,1) = c2
       dme(3,3) = c3
      end if
c.....Shear part
      if(ifs.eq.1) then
       call pzero(ds,4)
       ds(1,1) = c4
       ds(2,2) = c4
      end if
      return
      end
c
      subroutine sub34(bi,bj,da,dmat,s,ndf,nst,ino,jno,
     1 ncoli,nrowj,ncolj,irow,jcol)
c----------------------------------------------------------------------
c.....bi T * d * bj * det j * weig(xsi) * weig(eta)
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension bi(nrowj,ncoli),bj(nrowj,ncolj),dmat(nrowj,nrowj),
     1 dbmat(3,5),s(nst,nst),sbstf(5,5)
      ifrow = irow - 1
      jfcol = jcol - 1
c.....d*bj
      do 10 j = 1,ncolj
       do 10 i = 1,nrowj
        dbmat(i,j) = 0.0
        do 10 k = 1,nrowj
   10    dbmat(i,j) = dbmat(i,j) + dmat(i,k) * bj(k,j)
c.....bi T *(d*bj)
      do 20 j = 1,ncolj
       do 20 i = 1,ncoli
        sbstf(i,j) = 0.0
        do 20 k = 1,nrowj
   20    sbstf(i,j) = sbstf(i,j) + bi(k,i) * dbmat(k,j)
c.....store into K_T
      ifrow = (ino-1)*ndf + ifrow
      jfcol = (jno-1)*ndf + jfcol
      do 30 i = 1,ncoli
       irsub = ifrow + i
       do 30 j = 1,ncolj
        jcsub = jfcol + j
   30    s(irsub,jcsub) = s(irsub,jcsub) + sbstf(i,j)*da
      return
      end
c
      subroutine sfr34 (deriv,eta,xsi,nel,shape)
c----------------------------------------------------------------------
c.....Shape functions and derivatives for linear,quadratic
c.....lagrangian and serendipity isoparametric  2-d elements
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension deriv(2,*),shape(3,*)
      deins=1.0d0
      dzwei=2.0d0
      half=deins/dzwei
      viert=half/dzwei
      s = xsi
      t = eta
      if(nel.gt.4) goto 10
      st = s*t
c.....shape functions for 4-node element
      shape(3,1) = (deins-t-s+st)*viert
      shape(3,2) = (deins-t+s-st)*viert
      shape(3,3) = (deins+t+s+st)*viert
      shape(3,4) = (deins+t-s-st)*viert
c.....shape functions derivatives
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
c.....shape functions for 8-node element
      shape(3,1) = (-deins+st+ss+tt-sst-stt)*viert
      shape(3,2) = (-deins-st+ss+tt-sst+stt)*viert
      shape(3,3) = (-deins+st+ss+tt+sst+stt)*viert
      shape(3,4) = (-deins-st+ss+tt+sst-stt)*viert
      shape(3,5) = (+deins-t-ss+sst)*half
      shape(3,6) = (+deins+s-tt-stt)*half
      shape(3,7) = (+deins+t-ss-sst)*half
      shape(3,8) = (+deins-s-tt+stt)*half
c.....shape functions derivatives
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
c.....shape functions and derivatives for 9 noded element
      shape(3,1)=viert*s9*st*t9
      shape(3,5)=half*(deins-ss)*t*t9
      shape(3,2)=viert*s1*st*t9
      shape(3,6)=half*s*s1*(deins-tt)
      shape(3,3)=viert*s1*st*t1
      shape(3,7)=half*(deins-ss)*t*t1
      shape(3,4)=viert*s9*st*t1
      shape(3,8)=half*s*s9*(deins-tt)
      shape(3,9)=(deins-ss)*(deins-tt)
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
c
      subroutine gauss34(ngs,pgp,wgp)
c----------------------------------------------------------------------
c.....Gauss points  ngs = 1,4
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension pgp(4),wgp(4)
      call pzero(pgp,4)
      call pzero(wgp,4)
      if(ngs.gt.1) goto 2
      pgp(1) = 0.0
      wgp(1) = 2.0d0
      goto 16
    2 continue
      if(ngs.gt.2) goto 3
      pgp(1) = -dsqrt(3.0d0)/3.0d0
      wgp(1) =  1.0d0
      goto 6
    3 continue
      if(ngs.gt.3) goto 4
      pgp(1) = -dsqrt(0.6d0)
      pgp(2) =  0.0d0
      wgp(1) = 5.0d0/9.0d0
      wgp(2) = 8.0d0/9.0d0
      goto 6
    4 continue
      g = dsqrt(4.8d0)
      h = dsqrt(30.0d0)/36d0
      pgp(1) = -dsqrt((3.d0+g)/7.d0)
      pgp(2) = -dsqrt((3.d0-g)/7.d0)
      wgp(1) = 0.5d0 - h
      wgp(2) = 0.5d0 + h
    6 kgs = ngs/2
      do 8 igash = 1,kgs
       jgash = ngs + 1 - igash
       pgp(jgash) = -pgp(igash)
    8 wgp(jgash) =  wgp(igash)
   16 return
      end
c
      subroutine plots34(ix,dt,st,shape,sig,sigv,nel,da,numnp)
c----------------------------------------------------------------------
c.....extrpolate stresses to nodes  
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension dt(numnp),st(numnp,*),shape(3,*),sig(8),ix(*)
      do 10  j = 1,nel
       xsji = da*shape(3,j)
       ii = abs(ix(j))
       if(ii.eq.0) go to 10
       dt(ii) = dt(ii) + xsji
       st(ii,1) = st(ii,1) + sig(1)*xsji
       st(ii,2) = st(ii,2) + sig(2)*xsji
       st(ii,3) = st(ii,3) + sig(3)*xsji
       st(ii,4) = st(ii,4) + sig(4)*xsji
       st(ii,5) = st(ii,5) + sig(5)*xsji
       st(ii,6) = st(ii,6) + sig(6)*xsji
       st(ii,7) = st(ii,7) + sig(7)*xsji
       st(ii,8) = st(ii,8) + sig(8)*xsji
       st(ii,9) = st(ii,9) + sigv  *xsji
10    continue
      return
      end
c
      subroutine stre34(sign,dme,eps,sigt,ds,es,ifm,ifs)
c-----------------------------------------------------------------------
c.....Stresses for flat shell at position zs
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  sign(3),dme(3,3),eps(3),sigt(2),ds(2,2),es(2)
      if(ifm.eq.1) then
c..... normal stresses
       sign(1)= dme(1,1)*eps(1)+dme(1,2)*eps(2)
       sign(2)= dme(2,2)*eps(2)+dme(2,1)*eps(1)
       sign(3)= dme(3,3)*eps(3)
      end if
      if(ifs.eq.1) then
c..... shear stresses
       sigt(1)=ds(1,1)*es(1)
       sigt(2)=ds(2,2)*es(2)
      end if
      return
      end
c
      subroutine str34(h1,h2,nh,deps,dmat,d,sigout,sigv)
c-----------------------------------------------------------------------
c.... stresses and tangent matrix for plane stress J-2 plasticity with
c     linear isotropic hardening
c
c.... 1. Input parameters
c
c  d  - array of material constants
c  deps - incremental strains
c  epn  - plastic strains at t-n
c  epq  - effective plastic strain at t-n
c  sig  - stresses at t-n (in h1)
c  dlam0  - plastic multiplier at time t-n
c
c.... 2. Output parameters
c
c  epn  - plastic strains at t-n+1
c  epq  - effective platic strain at t-n+1
c  sig  - stresses at t-n+1 (in h2)
c  dmat - "tangent" matrix at t-n+1
c  dlam0  - plastic multiplier at time t-n+1
c-----------------------------------------------------------------------
      USE eldata
      implicit double precision (a-h,o-z)
      dimension d(*),h1(nh),h2(nh),deps(3),epn(3),sigtr(3),sig(3),
     1    dmat(3,3),sigout(3)
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

c.... plastic multiplier
      dlam0 = h1(9)

c.... trial stresses
      sigtr(1) = sig(1) + d(10)*deps(1) + d(11)*deps(2)
      sigtr(2) = sig(2) + d(11)*deps(1) + d(10)*deps(2)
      sigtr(3) = sig(3) + d(12)*deps(3)

c.... set up elastic tangent
      call pzero(dmat,9)
      dmat(1,1) = d(10)
      dmat(2,2) = d(10)
      dmat(1,2) = d(11)
      dmat(2,1) = d(11)
      dmat(3,3) = d(12)

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
c....  plastic correction necesarry
       call ylds34(d,sigtr,dmat,epn,epq,sigv,dlam0)
      end if
c.... set stresses etc. and  save time n+1 data from local array
      do 50 i = 1,3
50     h2(i) = sigtr(i)
      do 51 i = 1,3
51     h2(i+3) = epn(i)
      h2(7)   = epq
      h2(8)   = sigv/d(7)
      h2(9)   = dlam0
      do 52 i = 1,3
52     sigout(i) = sigtr(i)
      return
      end
c
      subroutine ylds34(d,sigtr,dmat,epn,epq,sigv,dlam0)
c-----------------------------------------------------------------------
c.... plane stress plasticity routine for return map algorithm
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension d(*),sigtr(3),sig(3),dmat(3,3),ps(3),psh(3),epn(3)
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
      x3=d(12)/(1.+2.*d(12)*dlam)

c.... calulate stresses sig = h**-1*d**-1*sigtr
      s1 = (sigtr(1)-d(2)*sigtr(2))/d(1)
      s2 = (sigtr(2)-d(2)*sigtr(1))/d(1)
      s3 = sigtr(3)/d(12)
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

c.... update dlam  from rlt
      if(ddlam.gt.0.0d0.or.abs(ddlam).lt.abs(dlam)) then
       dlam = dlam + ddlam
      else
       dlam = 0.5*abs(dlam)
      end if
c     dlam = dlam + ddlam
      icnt = icnt + 1
      if(icnt.gt.50) go to 110
c.... test for f
      if(abs(f).gt.tol) go to 100
c.... test for ddlam
c      if(abs(ddlam).gt.tol) go to 100
c.... test for ddlam   by RLT
c      if(abs(ddlam).gt.tol*abs(dlam)) go to 100
      go to 120
110   write(iow,2000)
2000  format(' * * Warning * * failure to converge in ylds34')
120   continue
c.... the term tol 1.d-8 is to ensure that the resulting stress is
c.... outside the current yield surface.  -->RLT!!!
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
      subroutine trans34(xl,yl,ul,vl,dul,dvl,t,ndm,ndf,nel,itrafo)
c---------------------------------------------------------------------------
c.... compute the transformation array and surface coords.
c     only for nel=4 and if desired(itrafo)
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension   x0(3),xl(ndm,4),yl(3,9),ul(ndf,4),vl(6,9),
     1      dul(6,9),dvl(6,9),t(3,3)

      if(nel.gt.4.or.itrafo.eq.0) then
        yl  = xl
        vl  = ul
        dvl =dul
        t=0.d0
        t(1,1)=1.d0
        t(2,2)=1.d0
        t(3,3)=1.d0

      else
c
c.... compute the inplane direction cosines (bisect diagonals)
      t=0
      do 100 i = 1,3
       t(1,i) = xl(i,3) - xl(i,1)
       t(2,i) = xl(i,2) - xl(i,4)
100   continue
      dl1 = sqrt(t(1,1)**2 + t(1,2)**2 + t(1,3)**2)
      dl2 = sqrt(t(2,1)**2 + t(2,2)**2 + t(2,3)**2)
      do 110 i = 1,3
       v1 = t(1,i)/dl1
       v2 = t(2,i)/dl2
       t(1,i) = v1 + v2
       t(2,i) = v1 - v2
110   continue
      dl1 = sqrt(t(1,1)**2 + t(1,2)**2 + t(1,3)**2)
      dl2 = sqrt(t(2,1)**2 + t(2,2)**2 + t(2,3)**2)
      do 120 i = 1,3
       t(1,i) = t(1,i)/dl1
       t(2,i) = t(2,i)/dl2
c....  compute the center (0,0) displacement
       x0(i) = 0.25*(xl(i,1) + xl(i,2) + xl(i,3) + xl(i,4))
120   continue
c.... compute the normal to the surface
      t(3,1) = t(1,2)*t(2,3) - t(2,2)*t(1,3)
      t(3,2) = t(1,3)*t(2,1) - t(2,3)*t(1,1)
      t(3,3) = t(1,1)*t(2,2) - t(2,1)*t(1,2)
c.... compute the projected middle surface coordinates
      do 140 i = 1,nel
      do 140 j = 1,3
       yl(j,i) = 0.0
       do 130 k = 1,3
        yl(j,i) = yl(j,i) + t(j,k)*(xl(k,i) - x0(k))
130     continue
140   continue
c.... set offset coordinates to zero if small compared to plan size
      htol =  0.0
      do 150 i = 1,nel
       htol = max(htol,abs(yl(1,i)),abs(yl(2,i)))
150   continue
      htol = htol*1.e-7
      do 160 i = 1,nel
       if(abs(yl(3,i)) .le. htol) yl(3,i) = 0.0
160   continue
c.... compute the transformation of displacements
      do 170 i = 1,nel
      do 170 j = 1,3
        vl(  j,i) = 0.0d0
        vl(3+j,i) = 0.0d0
       dvl(  j,i) = 0.0d0
       dvl(3+j,i) = 0.0d0
       do 180 k = 1,3
         vl(  j,i)  =  vl(  j,i) + t(j,k)* ul(  k,i)
         vl(3+j,i)  =  vl(3+j,i) + t(j,k)* ul(3+k,i)
        dvl(  j,i)  = dvl(  j,i) + t(j,k)*dul(  k,i)
        dvl(3+j,i)  = dvl(3+j,i) + t(j,k)*dul(3+k,i)
180    continue
170   continue

      end if
      return
      end
c
      subroutine trans134(s,p,t,nst,ndf,nel,itrafo)
c---------------------------------------------------------------------------
c.... transform the loads and stiffness to global coords.
c     only for nel=4 and if desired(itrafo)
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension s(nst,nst),p(nst),t(3,3),a(3,3),b(6)

      if(nel.gt.4.or.itrafo.eq.0) return

      i0 = 0
      do 170 ir = 1,nel
       do 110 ii = 1,3
        b(ii  ) = dot(t(1,ii),p(i0+1),3)
        b(ii+3) = dot(t(1,ii),p(i0+4),3)
110    continue
       do 111 ii = 1,6
        p(i0+ii) = b(ii)
111    continue
       j0 = i0
       do 160 jc = ir,nel
       i1 = i0
       do 150 i = 1,2
        j1 = j0
        do 140 j = 1,2
         do 120 ii = 1,3
          do 120 jj = 1,3
           a(jj,ii) = dot(t(1,ii),s(i1+1,jj+j1),3)
120      continue
         do 130 ii = 1,3
         do 130 jj = 1,3
          s(ii+i1,jj+j1) = dot(a(1,ii),t(1,jj),3)
130       continue
140     j1 = j1 + 3
150    i1 = i1 + 3
c....  compute the symmetric block
       if(ir.ne.jc) then
        do 155 i = 1,6
         do 155 j = 1,6
          s(j0+j,i0+i) = s(i0+i,j0+j)
155     continue
       end if
160    j0 = j0 + ndf
170   i0 = i0 + ndf
      return
      end
c
      subroutine getstr34(h,sig,n)
c---------------------------------------------------------------------------
c.... copy contents of h - field
c---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension h(*),sig(*)
      do i =1,n
        sig(i) = h(i)
      end do
      return
      end
