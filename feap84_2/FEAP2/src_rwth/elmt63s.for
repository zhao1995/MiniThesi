      subroutine elmt63(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c-----------------------------------------------------------------------
c
c.... geometrically nonlinear isogeometric shell element (based on El. 29 by FG/WW)
c     # enhanced strain formulation
c     # finite rotations, additive version 
c     # materials
c       1 small  strain isotropic J2 plasticity     
c       2 finite strain isotropic J2 plasticity     
c       3 Ogden´s material law     
c
c-----------------------------------------------------------------------
c.... input data
c.... 1. card: base system
c     1 d(01) = iRupd
c               1 = additive update of axial vector omega in every iteration
c               2 = multiplicative update of rotation matrix R in every iteration
c               3 = additive update within interation, multiplicative update after TIME       
c     2 d(02) = iRipl
c               0 = classic interpolation of director: interpolation after rotation, 
c                   global nodal basis systems for rotations and director
c               1 = classic interpolation of director: interpolation after rotation,                 
c                   global nodal basis systems for rotations and patch-wise basis systems for director      
c               2 = novel interpolation of director: rotation after interpolation,      
c                   global nodal basis systems for rotations and patch-wise basis systems for director  
c                   control points on kinks computed with T_3 and M with size 2,3 
c               3 = novel interpolation of director: rotation after interpolation
c                   omega is interpolated
c     3 d(03) = iCPupd
c               0   gauss point basis systems and control point basis systems are rotated  
c               1   only gauss point basis systems are rotated, 
c                   control point basis system are computed in every iteration from rotated gauss point systems
c                   This options works and perfectly computes nodal rotated control point systems. Results slightly
c                   change without showing significant change (untill now). Quadratic convergence in some cases lost
c     4 d(04) = 0   build only one part of stiffness matrix (only correct for symmetric matrices)
c               1   build complete stiffness matrix (required for unsymmetric stiffness matrix)      
c
c-----------------------------------------------------------------------
c.... 2. card: geometry etc.
c     1 d(11) = no. of layers
c     2 d(12) = <no. of Gauss Points/layer> def.=2
c     3 d(13) = shell thickness h_s
c     4 d(14) = z_h  = zeta (distance bottom surface to reference surface) 
c     5 d(15) = imat = type of material model
c     6 d(16) = ilin = 0=lin,1=nili (moderate),2=nili (finite) 
c     7 d(17) = ityp = 5=5dofs local rot,6=dofs global rot
c     8 d(18) = <nm/nb/ns enhanced parameter, each 0,4,5,7 as 3 digit number>
c     9 d(19) = shear correction factor 0=off, 1=on; for full locking set to 0 !!!!
c    10 d(20) = ANS switch 1=off, 0=on
c    11 d(21) = Number of Gauss points in one in-plane dimension 
c               (is overriden by value in base,21,x)        
c
c-----------------------------------------------------------------------
c.... 3. card: loads etc.
c     1 d(25) = element load q_n (transverse)
c     2 d(26) = element load q_x (global)
c     3 d(27) = element load q_y (global)
c     4 d(28) = element load q_z (global)
c     5 d(29) = alpha T
c     6 d(30) = temperature  bottom z-
c     7 d(31) = temperature  top    z+ 
c     8 d(32) = rho=gamma/g   per volume                              
c
c-----------------------------------------------------------------------
c.... 4. card: material data for type IMAT d(15)
c     1 d(40) = mat1   
c     2 d(41) = mat2   
c     3 d(42) = mat3   
c     4 d(43) = mat4   
c     .....
c-----------------------------------------------------------------------
c
c.... output of stresses see Subroutine plsn63 
c-----------------------------------------------------------------------
c
c.... allocation of h array  
c     h3)    rotational parameters:  3 at each node:                 3*4
c     h3)    mixed method:                                         npara
c     h2)    at each integration-point: nh=8  with 
c            1-6:     E_p = {E_p11, E_p22, E_p33, E_p12, E_p13, E_p23}
c            7:       e_v    equivalent plastic strain
c            8:       E_33
c                                             
c-----------------------------------------------------------------------
c
c     (c)  FG/WW  2009,  (c) WD  2011
c-----------------------------------------------------------------------
      USE arcl
      USE bdata
      USE cdata
      USE dirdat
      USE eldata
      USE hdata
      USE iofile
      USE qload
      USE plslay
      USE strnam
      USE isogeo
      USE errin1
      USE errin2
      USE errnam
      USE cdat1
      USE pdata10
      USE pltran
      USE prisdat
      USE prlod
      implicit double precision (a-h,o-z)
csk      common       m(1)
c
      parameter (neas_max=14)
c
      dimension xl(ndm,*),tl(*),d(*),ul(ndf,*),s(nst,*),p(*),ix(*),
     1          st(nen*ndf+neas_max,nen*ndf+neas_max),
     +          pt(nen*ndf+neas_max),h1(*),h2(*),h3(*),   
     2          shp(3,nen),shp1(3,4),sg(81),tg(81),wg(81),t0(3,3),
     3          dmat(8,8),sig(8),eps(8),epse(8),b(8,6,nen),btd(8,6),
     4          xs0(2,2),sx0(2,2),sx(2,2),xu(3,3),
     5          rd(3,3),twn(3,3,nen),rwn(3,3,nen),
     6       w(3,3,nen),x0(3,3),r0(3,3),te0(5,5),ts0(5,5),qs(2),siga(9),
     7          gam(4),ran(3,4),xan(3,4),bml(3,2,4),ddn(3,nen),
     8          pgz(5),wgz(5),sigp(6),epstr(5),cmat(5,5),
     9          gxy(8,neas_max),gtd(neas_max,8),eas(8),
     1          dd(3,3),w_1(3,3,nen),w_2(3,3,nen),
     2          twn_1(3,3,nen),twn_2(3,3,nen),iMsize(nen),t0exact(3,3,3)
      real*8 dwnold(3,nen),ddwn(3,nen),Mhq(3,3),Mh1(3,3),Mh2(3,3)

      dimension rinG(3,3,nen),rinP(3,3,nen),rwnPn(3,3,nen),
     1          rwnPnp1(3,3,nen),rwnGn(3,3,nen),rwnGnp1(3,3,nen),
     2          dwn(3,nen),dwg(9)
      save engyt
c
c.... go to correct array processor
      ngeas_max=nen*ndf+neas_max
      go to(1,2,3,3,5,3,7,3,3,2,2,2,2,2,2,2,2,18,2,2,2,22), isw

c.... input element properties
1      if(ior.lt.0) write(*,1000)
1000  format(' Input base system:',/
     +       ' igeo,add(1-3) >',$)
      call dinput(d,4)
c
      if(ior.lt.0) write(*,1001)
1001  format(' Input geometry etc.:',/ 
     +       ' nlay,ngp,h_s,z_h,imat,ilin,ityp,neas,stab,noans >',$)
      call dinput(d(11),11)
c
c     calculate EAS parameters 
      npeas= d(18)
      nm   = int(npeas/100) 
      nb1  = npeas-nm*100
      nb   = int(nb1/10)
      ns   = nb1-nb*10 
c
      if(ior.lt.0) write(*,1002)
1002  format(' Input loads etc.: q_n,q_x,q_y,q_z,a_t,t_-,t_+,rho >',$)
      call dinput(d(25),8)

c...  control input
      ngz = d(12)
      if(ngz.eq.0) d(12)=2.d0    ! no. of GP/layer 

c.... write element properties
      if(ior.lt.0) write(*  ,1004) 
     +         (d(i),i=1,4),(d(i),i=11,17),nm,nb,ns,d(19),d(20),
     +         d(21),(d(i),i=25,32) 
                   write(iow,1004) 
     +         (d(i),i=1,4),(d(i),i=11,17),nm,nb,ns,d(19),d(20),
     +         d(21),(d(i),i=25,32) 
1004  format(/5x,'Nonlinear isogeometric shell element ',
     +           ' element+material data:',//,
     + 5x,'iRupd 1=omega, 2=R, 3=omega+R@TIME... ',f8.0,/,
     + 5x,'iRipl 0=clas.,1=cl.+kink,2=novel,3=nv+',f8.0,/,
     + 5x,'iCPupd 0=comp.from.GP, 1=update.from.R',f8.0,/,
     + 5x,'unused........................... ',g12.4,//,
     + 5x,'no. of layers nlay ................. ',f8.0,/,
     + 5x,'no. of GP/layer (def.=2)............ ',f8.0,/,
     + 5x,'shell thickness h_s ............ ',g12.4,/,
     + 5x,'zeta_h ......................... ',g12.4,/,
     + 5x,'material type imat.................. ',f8.0,/,
     + 5x,'ilin 0=linear 1=moderate 2=finite... ',f8.0,/,
     + 5x,'ityp 5=5dofs local 6=6dofs global... ',f8.0,/,
     + 5x,'unused.............................  ',i4,'/',i1,'/',i1,/,
     + 5x,'shear correction factor 0=off,1=on.. ',g12.4,//,
     + 5x,'ANS  0=enabled,1=disabled........... ',f8.0,//,
     + 5x,'Number of 1D-surface-GP............. ',f8.0,//,
     + 5x,'element load q_n (local) ....... ',g12.4,/,
     + 5x,'element load q_x (global) ...... ',g12.4,/,
     + 5x,'element load q_y (global) ...... ',g12.4,/,
     + 5x,'element load q_z (global) ...... ',g12.4,/,
     + 5x,'alpha T ........................ ',g12.4,/,
     + 5x,'temperature  at z-.............. ',g12.4,/,
     + 5x,'temperature  at z+.............. ',g12.4,/,
     + 5x,'rho=gamma/g .................... ',g12.4)
c
c.... read material data for type imat 
      call matread63(d)
c
c.... initialize h-array    
c.... for eas:   kcr+kcc(neas*ngeas), rc(neas), alpha(neas)
      nlay = d(11)
      ngz  = d(12)                  ! integration points for each layer
      imat = d(15)
      neas = nm+nb+ns
      ngeas = nen*ndf + neas

c     check if number of integration points is given in base,21,x
      if (Iexactnormal.eq.1) then
        ngb = InumAreaGP
      else
        ngb = d(21)
        if (ngb.eq.0) then
c         if number of Gauss Points is not specified in input file  
          write(*,*) 'number of Gauss Points not specified->error in his
     +tory array WILL occur'
        end if
      end if

      npara = neas * (ngeas+2)      ! kcr+kcc(neas*ngeas), rc(neas), alpha(neas) 
      if (d(1).eq.1) then
c       + 3 rotational parameters at each node to store axial vector omega
        nhdr = npara + 3*nen        
        write(*,*) 
     +'   additive update of axial vector omega in every iteration'
        write(iow,*) 
     +'   additive update of axial vector omega in every iteration'
      else if (d(1).eq.2) then
        if (d(2).eq.2.or.d(2).eq.3.or.d(2).eq.1) then
c         + 9 entries at each node to store rotation matrix R
          nhdr = npara + 9*nen
          write(*,*) 
     +'   multiplicative update of rotation matrix R in every iteration'
          write(iow,*) 
     +'   multiplicative update of rotation matrix R in every iteration'
c        else if (d(2).eq.3) then
cc         no additional storage needed (3 rotational parameters each node  to store axial vector omega         
c          nhdr = npara + 3*nen 
c          write(*,*) 
c     +'   multiplicative update of rotation matrix R in every iteration 
c     +at the control points using Quaternions     '
c          write(iow,*) 
c     +'   multiplicative update of rotation matrix R in every iteration 
c     +at the control points using Quaternions     '
        end if
      else if (d(1).eq.3) then
c       + ? 
c        nhdr = npara + 9*nen  
      else !Standard case: set to additive update
c       + 3 rotational parameters at each node to store axial vector omega
        nhdr = npara + 3*nen             
      end if
      nh    = 8                     ! 1-6: E_p, 7: e_v, 8: E_33
      if(imat.eq.1.and.dabs(d(42)).le.1.d-18) nh = 0
      nh1   = nh*ngb*ngb*ngz*nlay  ! total number per element in h1/h2
c.... h3-field for rotation omega_n in all Gauss points
      if (d(2).eq.2) then
        write(*,*)   '   Rotation of INTERPOLATED basis systems - RIB'
        write(iow,*) '   Rotation of INTERPOLATED basis systems - RIB'
        if(d(1).eq.1) nhdr = nhdr + ngb*ngb*9 
        if(d(1).eq.2) nhdr = nhdr + ngb*ngb*9*3
      else if (d(2).eq.3) then
        write(*,*)   '   OMEGA SHELL - RIB'
        write(iow,*) '   OMEGA SHELL - RIB'
        if(d(1).eq.2) nhdr = nhdr + ngb*ngb*9*3
      else if (d(2).lt.2) then
        write(*,*)   '   Rotation of NODAL basis systems - RNB'
        write(iow,*) '   Rotation of NODAL basis systems - RNB'
      end if
      nh3   = nhdr                 ! total number per element in h3
csk   numerr =12
      e_name(1)  = 'H_1-norm w-w_ex'
      e_name(2)  = 'L_2-norm w-w_ex'
      e_name(3)  = 'sqrt(surface)  '
      e_name(4)  = 'H_1-norm w_ex  '
      e_name(5)  = 'L_2 ang.dev.dir'
      e_name(6)  = 'L_2 len.dev.dir'
      e_name(7)  = 'L_2 angl.dev.a1'
      e_name(8)  = 'L_2 leng.dev.a1'
      e_name(9)  = 'L_2 angl.dev.a2'
      e_name(10) = 'L_2 leng.dev.a2'
      e_name(11) = 'L_2 angl.dev.a3'
      e_name(12) = 'L_2 leng.dev.a3'
c     VERIFY IF e_name is long enough in errnam.h!      
c
2     return
3     iRupd  = d(1)
      iRipl  = d(2)
      iCPupd = d(3)
      isymm  = d(4)
      nlay  = d(11)
      imat  = d(15)
      ilin  = d(16) 
      ityp  = d(17) 
      npeas = d(18)       ! EAS-Input
      icappa = d(19)      ! shear correction factor on/off
cwd   set ANS switch 
      noans = d(20)
      
      nm    = int(npeas/100) 
      nb1   = npeas-nm*100
      nb    = int(nb1/10)
      ns    = nb1-nb*10 
      neas  = nm+nb+ns ! total no. EAS
c     check if number of integration points is given in base,21,x
      if (Iexactnormal.eq.1) then
        ngb = InumAreaGP
      else
        ngb = d(21)
        if (ngb.eq.0) then
c         if number of Gauss Points is not specified in input file  
          write(*,*) 'number of Gauss Points not specified->error in his
     +tory array WILL occur'
        end if
      end if

c.... get current patch and coressponding orders
      NURp = ngetNURnmpq(n,3,AInipa,AInmpq)
      NURq = ngetNURnmpq(n,4,AInipa,AInmpq)
c      if (ngb.eq.0) then
c       set number of Gauss Points to ceil(p/2)+1(see Hughes, Realli and Sangalli)
c       if number of Gauss Points is not specified in input file  
c        write(*,*) 'number of Gauss Points not specified->error in histo
c     +ry array WILL occur'
c        maxIGAorder=max(NURp,NURq)
c        ngb = ceiling(maxIGAorder/2.0d0)+1
c      end if
c      if (isw.eq.9) ngb = 7
      ngz   = d(12)
      nh    = 8
      ngeas = 4*ndf + neas
      dh  = d(13)/nlay
      dh5 = 0.5d0*dh
      ielas = 0 
      if(imat.eq.1.and.dabs(d(42)).lt.1.d-18)ielas = 1
      if(ielas.eq.1) nh = 0

      maxlay = ngz*nlay
      if(klay.gt.maxlay) klay = maxlay
c.... names for stresses
      if(isw.eq.8) then
        if(klay.eq.0) then
          call plsn63 (0)
        else
          call plsn63 (1)
        end if 
      end if
c      
c.... storage in h3-array: 
c     omega(3*4),  sc(neas*ngeas), pc(neas), alpha(neas)
c     
      nhom  = 1
      if (iRipl.ge.2) then
        if (iRupd.eq.1) then
          if (iRipl.eq.2)
     +     ndwn  = ngb*ngb*9+nhom
          if (iRipl.eq.3) ndwn = nhom
          nsc   = ndwn + 3*nen
        elseif (iRupd.eq.2)  then
c          if (iRipl.eq.3) then
c            ndwn = ngb*ngb*9+nhom
c            nsc = ndwn + 3*nen
c          else
            ndwn  = ngb*ngb*27+nhom
            nsc   = ndwn + 3*nen          
c          endif
        endif
      else
        if (iRupd.eq.1) then
          nsc = nhom + 3*nen
        elseif (iRupd.eq.2)  then
          nsc = nhom + 9*nen
        endif
      endif  
      npc   = nsc + neas*ngeas
      nalp  = npc + neas

c     For iRupd check if Rotation matrix has to be initialized
      if (iRupd.eq.2) then
        if (iRipl.lt.2)
     +  call InitRotationMatrix63(AInRupdinit,n,h3,numel,nen)
        if (iRipl.ge.2)
     +  call InitRotationMatrix63CPGP(AInRupdinit,n,h3(ndwn),
     +                            h3,numel,nen,ngb)
      end if  
c     HERE THE REAL ELEMENT STARTS
c.... compute NURBS coordinates of recent element      
      ni = ngetNURni(n,AIninc,AInien,AInipa,1)
      nj = ngetNURni(n,AIninc,AInien,AInipa,2)
c.... check if element has zero measure
      if (rNURknv1(ni).eq.rNURknv1(ni+1).or.
     +    rNURknv2(nj).eq.rNURknv2(nj+1)) return
c
      call pzero (pt,ngeas_max)
      call pzero (st,ngeas_max*ngeas_max)                    
      call pzero (eas,8)
      call IGAgauss (ngb,lint,sg,tg,wg)
c      call gaus1D (ngz,pgz,wgz)  !ONLY NEEDED FOR ANS
      call rjac63(d,xl,t0,xs0,sx0,ts0,te0,detj0,isw,ni,nj,NURp,NURq,ndm
     +             ,AInkv1,AInkv2)  !CHANGED
      call shearfac63(d,detj0,eps,sig,dmat,icappa,cappa,wcappa) !xl removed
      if(ielas.eq.1)call dmate63 (d,cappa,dmat)   !NO CHANGES 

      if(ldir.eq.0) then  !Macro BASE not called, only igeo in MATE entered
        call drawmess('CALL MACRO BASE',1,-2)
      else if(ldir.eq.1) then                        ! read from m(mdir)
        do k = 1,nen 
          if (knode.eq.numnp) then
csk         read from basea=m(mdir)              
            call pdirecread63(ix,basea,rinG(1,1,k),k,2,1,iRipl,
     +                        rinP(1,1,k),n)
          else
            write(*,*) 'error while loading nodal basis systems'
          end if
        end do
      end if
c     do Rodrigues update nodewise for standard formulation (iRipl=0)     
      if (iRipl.eq.0.or.iRipl.eq.1) then
        if (iRupd.eq.1)
     1    call updn63 (xl,ul,h3,twn,rwn,ddn,rinG,ran,        !CHANGED 
     2                 xan,gam,bml,w,ix,ndf,ndm,nen,ityp,ilin)
        if (iRupd.eq.2)
     1    call updn63Rupd2 (xl,ul,h3,twn,rwn,ddn,rinG,ran,        !CHANGED 
     2                 xan,gam,bml,w,ix,ndf,ndm,nen,ityp,ilin,dwn)
      end if
c      if (iRipl.eq.1)
c     1    call updn63 (xl,ul,h3,twn,rwn,ddn,rinG,ran,        !CHANGED 
c     2                 xan,gam,bml,w,ix,ndf,ndm,nen,ityp,ilin)
      if (iRipl.ge.2) then
c        if (iRipl.eq.2.and.iRupd.eq.1)
        if (iRupd.eq.1)
     1   call updCP63 (xl,ul,h3(ndwn),ix,ndf,ndm,nen,ityp,ilin,
     2                rinG,rinP,rwnGn,rwnGnp1,rwnPn,rwnPnp1)
c        if (iRipl.eq.3)
c     1   call updCP63Quart (xl,ul,h3(ndwn),ix,ndf,ndm,nen,ityp,ilin,
c     2                rinG,rinP,rwnGn,rwnGnp1,rwnPn,rwnPnp1,
c     3                dwnold,ddwn,iRupd)        
c        if (iRipl.eq.2.and.iRupd.eq.2)
        if (iRupd.eq.2.and.iRipl.eq.2)
     1   call updCP63Rupd2 (xl,ul,h3(ndwn),ix,ndf,ndm,nen,ityp,ilin,
     2                rinG,rinP,rwnGn,rwnGnp1,rwnPn,rwnPnp1)
        if (iRupd.eq.2.and.iRipl.eq.3)
     1   call updCP63OmegaMult (xl,ul,h3(ndwn),ix,ndf,ndm,nen,ityp,ilin,
     2                rinG,rinP,rwnGn,twn,rwnPn,rwnPnp1,ddwn)
        if (iCPupd.eq.1)
     1    call computeCPfromGP63 (xl,ul,h3,ix,ilin,rinG,rinP,lint,
     1                            rwnGn,rwnGnp1,rwnPn,rwnPnp1,
     2                            sg,tg,ndm,nen,ndf,
     3                            ni,nj,NURp,NURq,AInkv1,AInkv2,
     4                            Iexactnormal,AInbasissys,ngb,n)  
      end if
c
c.... loop over gauss points
      nhpl = 1
c      write(iow,*) n
      if (n.eq.10) then
c      write(iow,*) h3(ndwn:300)
      end if
      do 300 l = 1,lint
       call pzero (eps,8)
       call pzero (epse,8)
       xsi = sg(l)
       eta = tg(l)
       call shap63(xsi,eta,xl,t0,sx,shp,shp1,xsj,ndm,nen,ni,nj,NURp,
     +            NURq,AInkv1,AInkv2)    !CHANGED
       da = xsj*wg(l)
c      do Rodrigues update in Gauss points (iRipl.ne.0)
c      and determine deltadirector in Gauss point
       if (iRipl.ge.2) then
         if (Iexactnormal.eq.1) then
c          load exact basis system from IGAbasissystem to t0exact    
c           t0exact = IGAbasissystems(1:3,1:3,1:3,l,n)
           call IGAreadNormalm(AInbasissys,numel,ngb**2,t0exact,l,n)
c          test of implementation         
           
           call testtripod63vnorm(t0exact,
     +                       rinP,nen,shp,t0,da,isw,n,l)
c          call testdirector63vnorm
c     +      (t0exact(1:3,1:3,3),t0,da,isw)          
         else
           t0exact = 0.d0
         end if
        
c         if (iRipl.eq.4) call updGP63omInt(xl,ul,h3,twn,dd,
c     2                 w,ix,ndf,ndm,nen,ityp,ilin,
c     3                 l,shp,rd,r0,w_1,w_2,
c     4                 twn_1,twn_2,iMsize,rinG,rinP,
c     5                 rwnGn,rwnGnp1,rwnPn,rwnPnp1,t0exact,Iexactnormal,
c     6                 h3(ndwn))
         if (iRipl.eq.3) then
c           if (iRupd.eq.2) call updGP63omega(xl,ul,h3,twn,dd,  !QUATERNIONEN SACHE; GEHT NICHT GESCHEID
c     2                 w,ix,ndf,ndm,nen,ityp,ilin,
c     3                 l,shp,rd,r0,w_1,w_2,
c     4                 twn_1,twn_2,iMsize,rinG,rinP,
c     5                 rwnGn,rwnGnp1,rwnPn,rwnPnp1,t0exact,Iexactnormal,
c     6                 h3(ndwn),dwnold,ddwn)
           if (iRupd.eq.2) call updGP63omegaMult(xl,ul,dwg,twn,dd,
     2                 w,ix,ndf,ndm,nen,ityp,ilin,
     3                 l,shp,rd,r0,w_1,w_2,
     4                 iMsize,rinP,
     5                 t0exact,Iexactnormal,
     6                 ddwn,h3)
           if (iRupd.eq.1) call updGP63omegaAdd(xl,ul,dwg,twn,dd,
     2                 w,ix,ndf,ndm,nen,ityp,ilin,l,shp,rd,r0,w_1,
     3                 w_2,twn_1,twn_2,iMsize,rinG,rinP,
     5                 rwnGnp1,t0exact,Iexactnormal,h3(ndwn))              
         end if
         if (iRipl.eq.2) then
           if (iRupd.eq.1) call updGP63(xl,ul,h3,twn,dd,
     2                 w,ix,ndf,ndm,nen,ityp,ilin,
     3                 l,shp,rd,r0,w_1,w_2,
     4                 twn_1,twn_2,iMsize,rinG,rinP,
     5                 rwnGn,rwnGnp1,rwnPn,rwnPnp1,t0exact,Iexactnormal)
           if (iRupd.eq.2) call updGP63Rupd2(xl,ul,dwg,twn,dd,
     2                 w,ix,ndf,ndm,nen,ityp,ilin,
     3                 l,shp,rd,r0,w_1,w_2,
     4                 twn_1,twn_2,iMsize,rinG,rinP,
     5                 rwnGn,rwnGnp1,rwnPn,rwnPnp1,t0exact,Iexactnormal,
     6                 h3)
         end if
       end if
c      write(iow,*) rinG
       
       surface=surface+da    !   for test of surface area
       rj = detj0/xsj
c     
       if (iRipl.eq.1) then
         call stra63 (xl,ul,shp,xsi,eta,eas,eps,sx,gam,rinP,rwn,  !CHANGED
     1              ddn,xu,x0,rd,r0,ndm,ndf,ilin,noans,nen,Tgp,
     2              iRipl,dd)
       else
         call stra63 (xl,ul,shp,xsi,eta,eas,eps,sx,gam,rinG,rwn,  !CHANGED
     1              ddn,xu,x0,rd,r0,ndm,ndf,ilin,noans,nen,Tgp,
     2              iRipl,dd)
       end if
c       write(iow,*) x0(1:3,3)
c
       if (Iexactnormal.ne.1) then
         if (isw.eq.9) 
     +     call testtripod63vnormIn(rinP,nen,shp,t0,da,isw,n,l)
       end if
c      if (l.eq.1) write(iow,*) n
c      write(iow,*) r0(1,3)-t0(1,3),r0(2,3)-t0(2,3),r0(3,3)-t0(3,3)
c      write(iow,*) twn
c...   for isw=9 perform error analysis, discrete values for plate
!       if (isw.eq.9) then
!         x=x0(1,3)
!         y=x0(2,3)
!         ansol=1.d0/3.d0*x**3*(x-1.d0)**3*y**3*(y-1.d0)**3-
!     -         2*d(13)**2/(5.d0*(1.d0-d(41)))*(y**3*(y-1.d0)**3*
!     *         x*(x-1.d0)*(5.d0*x**2-5.d0*x+1.d0)+
!     +         x**3*(x-1.d0)**3*y*(y-1.d0)*(5.d0*y**2-5.d0*y+1.d0))
!         ansol_dx=(y**3*(y-1.d0)**3*x**2*(x-1.d0)**2*(2.d0*x-1)-
!     -            2*d(13)**2/(5.d0*(1.d0-d(41)))*(y**3*(y-1.d0)**3*
!     *            (20.d0*x**3-30.d0*x**2+12.d0*x-1)+(y**2-y)*
!     *            (5.d0*y*y-5.d0*y+1.d0)*3.d0*x*x*(x-1.d0)**2*
!     *            (2.d0*x-1.d0)))
!         y=x0(1,3)
!         x=x0(2,3) 
!         ansol_dy=(y**3*(y-1.d0)**3*x**2*(x-1.d0)**2*(2.d0*x-1)-
!     -            2*d(13)**2/(5.d0*(1.d0-d(41)))*(y**3*(y-1.d0)**3*
!     *            (20.d0*x**3-30.d0*x**2+12.d0*x-1)+(y**2-y)*
!     *            (5.d0*y*y-5.d0*y+1.d0)*3.d0*x*x*(x-1.d0)**2*
!     *            (2.d0*x-1.d0)))
!         do i = 1,3
!           do j = 1,3
!             xu(i,j) = 0.d0
!             do k  = 1,nen
!               xu(i,j) = xu(i,j) + ul(i,k)*shp(j,k)
!             end do
!           end do
!         end do 
!         cosol=xu(3,3)
!         cosol_dx=xu(3,1)*t0(1,1)+xu(3,2)*t0(1,2)
!         cosol_dy=xu(3,1)*t0(2,1)+xu(3,2)*t0(2,2)
!         
!c         cosol_dx=abs(xu(3,1))
!c         cosol_dy=abs(xu(3,2))
!         call err63(ansol-cosol,ansol_dx-cosol_dx,ansol_dy-cosol_dy,da,
!     +              ansol,ansol_dx,ansol_dy)
!         go to 300
!       end if       

       if (isw.eq.9) call testdirector63vnorm(r0,t0,da,isw)
       call strt63 (d,aqloa,numel,n,mqloa,propq,eps,epse,isw)   !NO CHANGES 
c        
      if(ielas.eq.1)then
        call mvmul(dmat,epse,8,8,sig)
        if(n.eq.1.and.l.eq.1)then
         engyt = 0.d0
        else
         engyt = engyt +  0.5d0*dot(sig(1),epse(1),8)*da
         detc  = engyt
        end if

        call mvmul(dmat,eps,8,8,sig)
cwd     compute the energetic norm for validation of nonlinear stuff   
        if (isw.eq.9) then
          energetic_norm = dot(sig(1),eps(1),8)*da
          call energeticnorm63(energetic_norm)
          go to 300
        end if
      else
       


c....   integration through the thickness

        call pzero (sig,8)
        call pzero (dmat,8*8)

        do 363 ilay = 1,nlay
           dz = d(14) + (ilay-1)*dh
           do 360 igz = 1,ngz
              wz   = wgz(igz)*dh5             
              zeta = pgz(igz)
              zs   = dz + (1.d0+zeta)*dh5
              call kine63(eps,x0,r0,zs,wcappa,epstr,dets)   !NO CHANGES 
              wz = wz*dets   
c
              call zerostr63(h1(nhpl),h2(nhpl),d,epstr,x0,   !NO CHANGES
     +                       sigp,cmat)

c.....        print/plot stresses in layer klay and GP mlay
              if (isw.eq.4 .or. isw.eq.8) then  
               if (ilay.eq.klay.and.igz.eq.mlay) then 
                 call trasig63 (sigp,t0,shp,rinP,nen)           !CHANGED
                 
                 if (isw.eq.4) call pout63 (xl,shp,sig,sigp, !NO CHANGES
     +                              h2(nhpl),ndm,maxlay)  
                 if (isw.eq.8) then   !PLOT STRESS ROUTINE
csk                   iAInoix = AInoix+(n-1)*(nen)
csk                   call plot63 (d,ix,m(np),m(np+numnp*ipr),     !TO BE CHANGED !!!!!!!!!!!!!!!!!!!!!!!!!
                   call plot63 (d,ix,strea,strea(1+numnp),     
     +               shp,sig,sigp,h2(nhpl),eps,numnp,maxlay,klay,da,nen,
     +                AInoix)
                 end if 
               end if
              end if

c.....        compute material matrix and stress resultants
              call dmat63 (dmat,cmat,sig,sigp,wz,zs,wcappa,cappa)  !NO CHANGES
              nhpl = nhpl + nh
360        continue
363     continue
       end if

       if (isw.eq.4 .or. isw.eq.8) then
        if (klay.eq.0) then         ! print/plot stress resultants
c....    transformation
         call transsig63(t0,xu,rd,rinP,shp,sig,siga,nen)    !CHANGED
         if (isw.eq.4) call pout63 (xl,shp,siga,sigp,       !NO CHANGES
     +                              h2,ndm,maxlay)
         if (isw.eq.8.and.iplma(ma).ne.0) then              !TO BE CHANGED !!!!!!!!!!!!!!!!!!!!!!!!!
csk            iAInoix = AInoix+(n-1)*(nen) 
           call plot63 (d,ix,strea,strea(1+numnp),shp,siga,
     +                   sigp,h2,eps,numnp,maxlay,klay,da,nen,
     +                    AInoix)
         end if
        end if
       end if
c         
       if (isw.eq.4.or.isw.eq.8) go to 300
c....  compute input values for Gww for Multiplicative update formulations (2,3)
       !if (iRupd.eq.2.and.iRipl.eq.3) then
       if (iRupd.eq.2.and.iRipl.ge.2) then
         call computeGww63_Mult_inputvalues(sig,xu,Mhq,Mh1,Mh2,da,rd)
       end if
c....  loop over rows, residual
       ir0 = 0
       if (isymm.eq.0) then
         do 310 inode = 1,nen
           jnode=1
           if (iRupd.eq.1) then 
             call bmat63 (inode,ir0,b,btd,dmat,st,pt,sig,qs,rwn,h3, 
     1       w,twn,bml,d,shp1,shp,sx,t0,ran,xan,xu,rd,da,ngeas_max,
     2       ndf,ilin,noans,x0,iRipl,l,w_1,w_2,isymm,jnode) 
           else if (iRupd.eq.2) then
             call bmat63 (inode,ir0,b,btd,dmat,st,pt,sig,qs,rwn,dwn, 
     1       w,twn,bml,d,shp1,shp,sx,t0,ran,xan,xu,rd,da,ngeas_max,
     2       ndf,ilin,noans,x0,iRipl,l,w_1,w_2,isymm,jnode)
           end if

 

c....      loop over columns (symmetry noted), stiffness matrix
           jc0 = 0       
           do jnode = 1,inode
             if (iRupd.eq.1) then
               if (iRipl.eq.3) then
                 ifake = 1
                 call stif63 (inode,jnode,ir0,jc0,b,btd,st,sig,qs,shp,
     1                shp1,w,ngeas_max,ndf,da,ilin,isw,noans,iRipl,xu,
     2                dwg,rd,twn,ifake,w_1,w_2,twn_1,twn_2,iMsize,iRupd,
     3               Mhq,Mh1,Mh2)
              else
                call stif63 (inode,jnode,ir0,jc0,b,btd,st,sig,qs,shp,
     1                shp1,w,ngeas_max,ndf,da,ilin,isw,noans,iRipl,
     2                xu,h3,rd,twn,l,w_1,w_2,twn_1,twn_2,iMsize,iRupd,
     3               Mhq,Mh1,Mh2)
              end if
             else
               if (iRipl.lt.2)
     1           call stif63 (inode,jnode,ir0,jc0,b,btd,st,sig,qs,shp, 
     1                shp1,w,ngeas_max,ndf,da,ilin,isw,noans,iRipl,
     2                xu,dwn,rd,twn,l,w_1,w_2,twn_1,twn_2,iMsize,iRupd,
     3               Mhq,Mh1,Mh2)
               if (iRipl.eq.2) then
                 ifake = 1
                 call stif63 (inode,jnode,ir0,jc0,b,btd,st,sig,qs,shp,
     1                shp1,w,ngeas_max,ndf,da,ilin,isw,noans,iRipl,xu,
     2                dwg,rd,twn,ifake,w_1,w_2,twn_1,twn_2,iMsize,iRupd,
     3               Mhq,Mh1,Mh2)
               end if
               if (iRipl.eq.3) then
!c....  compute input values for Gww for Multiplicative update formulations (2,3)
!       if (iRupd.eq.2.and.iRipl.eq.3) then
!         call computeGww63_Mult_inputvalues(sig,xu,Mhq,Mh1,Mh2,da,rd)
!       end if
                 ifake = 1
                 call stif63 (inode,jnode,ir0,jc0,b,btd,st,sig,qs,shp,
     1                shp1,w,ngeas_max,ndf,da,ilin,isw,noans,iRipl,xu,
     2                dwg,rd,twn,ifake,w_1,w_2,twn_1,twn_2,iMsize,iRupd,
     3               Mhq,Mh1,Mh2)
               end if
             end if
             jc0 = jc0 + ndf
          end do
310       ir0 = ir0 + ndf
        else if (isymm.eq.1) then
c....   loop over rows, residual        
          do 320 inode = 1,nen
c....       loop over columns, (complete) stiffness matrix
            jc0 = 0       
            do jnode = 1,nen
              if (iRupd.eq.1) then 
                call bmat63 (inode,ir0,b,btd,dmat,st,pt,sig,qs,rwn,h3, 
     1            w,twn,bml,d,shp1,shp,sx,t0,ran,xan,xu,rd,da,ngeas_max,
     2            ndf,ilin,noans,x0,iRipl,l,w_1,w_2,isymm,jnode) 
              else if (iRupd.eq.2) then
                call bmat63 (inode,ir0,b,btd,dmat,st,pt,sig,qs,rwn,dwn, 
     1            w,twn,bml,d,shp1,shp,sx,t0,ran,xan,xu,rd,da,ngeas_max,
     2            ndf,ilin,noans,x0,iRipl,l,w_1,w_2,isymm,jnode)
              end if

              if (iRupd.eq.1) then
                if (iRipl.eq.3) then
                 ifake = 1
                 call stif63 (inode,jnode,ir0,jc0,b,btd,st,sig,qs,shp,
     1                shp1,w,ngeas_max,ndf,da,ilin,isw,noans,iRipl,xu,
     2                dwg,rd,twn,ifake,w_1,w_2,twn_1,twn_2,iMsize,iRupd,
     3               Mhq,Mh1,Mh2)
                else
                  call stif63 (inode,jnode,ir0,jc0,b,btd,st,sig,qs,shp,
     1                shp1,w,ngeas_max,ndf,da,ilin,isw,noans,iRipl,
     2                xu,h3,rd,twn,l,w_1,w_2,twn_1,twn_2,iMsize,iRupd,
     3               Mhq,Mh1,Mh2)
                end if
              else
                if (iRipl.lt.2)
     1            call stif63 (inode,jnode,ir0,jc0,b,btd,st,sig,qs,shp, 
     1                shp1,w,ngeas_max,ndf,da,ilin,isw,noans,iRipl,
     2                xu,dwn,rd,twn,l,w_1,w_2,twn_1,twn_2,iMsize,iRupd,
     3               Mhq,Mh1,Mh2)
                if (iRipl.eq.2) then
                  ifake = 1
                  call stif63 (inode,jnode,ir0,jc0,b,btd,st,sig,qs,shp,
     1                shp1,w,ngeas_max,ndf,da,ilin,isw,noans,iRipl,xu,
     2                dwg,rd,twn,ifake,w_1,w_2,twn_1,twn_2,iMsize,iRupd,
     3               Mhq,Mh1,Mh2)
                end if
                if (iRipl.eq.3) then
                  ifake = 1
                  call stif63 (inode,jnode,ir0,jc0,b,btd,st,sig,qs,shp,
     1                shp1,w,ngeas_max,ndf,da,ilin,isw,noans,iRipl,xu,
     2                dwg,rd,twn,ifake,w_1,w_2,twn_1,twn_2,iMsize,iRupd,
     3               Mhq,Mh1,Mh2)
                end if
              end if
              jc0 = jc0 + ndf
            end do
320       ir0 = ir0 + ndf
        end if
300   continue
         if (isw.eq.4) go to 4
         if (isw.eq.8) go to 8

c
c.... copy stiffness matrix st to s and pt to p 
c     SYMMETRICAL
      if (isymm.eq.0) then
        do i = 1,nen*ndf
          p(i) = pt(i)
          do j = 1,i
            s(i,j) = st(i,j)
            s(j,i) = st(i,j)
          enddo
        enddo          
      else
c     NON-SYMMETRICAL     
        do i = 1,nen*ndf
          p(i) = pt(i)
          do j = 1,nen*ndf
            if (st(i,j)-st(j,i).gt.1d-9) then
              write(*,*) st(i,j)-st(j,i)
            end if
            s(i,j) = st(i,j)
          enddo
        enddo  
      end if
c     HERE THE REAL ELEMENTS STOPS; THE REST IS JUST POST-PROCESSING

           
     
      return
c
c.... output stresses   =>  isw 3!
4     return
c.... compute lumped mass matrix 
c     open: 5/6 dof?, value for dof 6? 
5     l = 2
      call pgauss(l,lint,sg,tg,wg)
      call rjac63(d,xl,t0,xs0,sx0,ts0,te0,detj0,isw)!to be changed
      hs = d(13)
      dh = hs*hs/12.d0
      do 50 l = 1,lint
        xsi = sg(l)
        eta = tg(l)
        call shap63(xsi,eta,xl,t0,sx,shp,shp1,xsj)
        da = 1.d0*xsj*hs*d(32)
        i1 = 0
        do 51 i = 1,nel
          dpf = shp(3,i)*da 
          p(i1+1) = p(i1+1) + dpf 
          p(i1+2) = p(i1+2) + dpf
          p(i1+3) = p(i1+3) + dpf
          p(i1+4) = p(i1+4) + dpf*dh
          p(i1+5) = p(i1+5) + dpf*dh
          if(ndf.gt.5) p(i1+6) = p(i1+6) + dpf*dh
          i1 = i1 + ndf
51      end do
50    end do
      return
c
c.... compute the surface tractions  (Input: SLOA)
7     return
c.... plot stresses   =>  isw 3!
8     return
c     error analysis  => isw 3
9     return  
c.... calculate initial nodal cartesian basis
18    continue
      igeo = d(1)
c.... get current patch and coressponding orders
      NURp = ngetNURnmpq(n,3,AInipa,AInmpq)
      NURq = ngetNURnmpq(n,4,AInipa,AInmpq)
c.... compute NURBS coordinates of recent element      
      ni = ngetNURni(n,AIninc,AInien,AInipa,1)
      nj = ngetNURni(n,AIninc,AInien,AInipa,2)
c.... check if element has zero measure
      if (rNURknv1(ni).eq.rNURknv1(ni+1).or.
     +    rNURknv2(nj).eq.rNURknv2(nj+1)) return
      call rjac63(d,xl,t0,xs0,sx0,ts0,te0,detj0,isw,ni,nj,NURp,NURq,ndm
     +             ,AInkv1,AInkv2)

c      call triad63 (igeo,xl,d,rin,ndm,nel)

      if (knode.eq.(nen*numel)) then
c....  BASE,1-2       
       do k = 1,nen 
csk     read from basea=m(mdir)              
        call pdirec2 (basea,t0,k,n,nen,1,1)
        call pdirec2 (basea,t0,k,n,nen,1,2)
       end do 
      else if (knode.eq.numnp) then
c....  BASE,0       
       do k = 1,4 
csk     read from basea=m(mdir)              
        call pdirec1(ix,basea,rinG(1,1,k),k,1,1)
        call pdirec1(ix,basea,rinG(1,1,k),k,1,2)
       end do 
      end if
      return
c.... compute the surface tractions  (Input: QLOA)  temperature loading
c.... compute the surface tractions  (Input: QLOA)  follower forces 
22    call testq63(d,aqloa,numel,n,iflg1,iflg2)
c
      if(iflg2.eq.1)then
      igeo  = d(1)
      nlay  = d(11)
      ilin  = d(16) 
      ityp  = d(17) 
      npeas = d(18)       ! EAS-Input
      nm    = int(npeas/100) 
      nb1   = npeas-nm*100
      nb    = int(nb1/10)
      ns    = nb1-nb*10 
      neas  = nm+nb+ns ! total no. EAS
      ngb   = 2
      ngz   = d(12)
      nh    = 8
      ngeas = 4*ndf + neas
      dh  = d(13)/nlay
      dh5 = 0.5d0*dh

      imat = d(15)
      ielas = 0 
      if(imat.eq.1.and.dabs(d(42)).lt.1.d-18)ielas = 1
      if(ielas.eq.1) nh = 0
      cappa = 1.d0
      if(ielas.eq.1)call dmate63 (d,cappa,dmat)
c      
c      
      nhom  = 1
c
      call pzero (pt,ngeas_max)
      call pzero (st,ngeas_max*ngeas_max)                    
      call pzero (eas,8)
      call pgauss (ngb,lint,sg,tg,wg)
      call gaus1D (ngz,pgz,wgz)
      call rjac63(d,xl,t0,xs0,sx0,ts0,te0,detj0,isw)
      if(ldir.eq.0) then 
        call triad63 (igeo,xl,d,rinG)            
      else if(ldir.eq.1) then 
        do k = 1,4 
csk       read from basea=m(mdir)              
          call pdirec1(ix,basea,rinG(1,1,k),k,2,1)
        end do 
      end if

      call updn63 (xl,ul,h3,twn,rwn,ddn,rinG,ran,
     1             xan,gam,bml,w,ix,ndf,ndm,nen,ityp,ilin)
c
c
c.... loop over gauss points
      nhpl = 1
      do 2200 l = 1,lint
       xsi = sg(l)
       eta = tg(l)
       call shap63(xsi,eta,xl,t0,sx,shp,shp1,xsj)
       da = xsj*wg(l)
       rj = detj0/xsj

       call stra63 (xl,ul,shp,xsi,eta,eas,eps,sx,gam,rinG,rwn,
     1              ddn,xu,x0,rd,r0,ndm,ndf,ilin)
       call pzero (eps,8)
       call strt63 (d,aqloa,numel,n,mqloa,propq,eps,epse,isw)

       if(ielas.eq.1)then
        call mvmul(dmat,eps,8,8,sig)
       else

c....   integration through the thickness
        call pzero (sig,8)
        call pzero (dmat,8*8)

        do 2261 ilay = 1,nlay
           dz = d(14) + (ilay-1)*dh
           do 2260 igz = 1,ngz
              wz   = wgz(igz)*dh5             
              zeta = pgz(igz)
              zs   = dz + (1.d0+zeta)*dh5
              call kine63(eps,x0,r0,zs,wcappa,epstr,dets)
              wz = wz*dets   
c
              call zerostr63(h1(nhpl),h2(nhpl),d,epstr,x0,
     +                       sigp,cmat)

c.....        compute material matrix and stress resultants
              call dmat63 (dmat,cmat,sig,sigp,wz,zs,wcappa,cappa)
              nhpl = nhpl + nh
2260        continue
2261     continue
       end if

c         
c....  loop over rows, residual
       ir0 = 0
       do 2210 inode = 1,4
        call bmat63 (inode,ir0,b,btd,dmat,st,pt,sig,qs,rwn,h3,
     1   w,twn,bml,d,shp1,shp,sx,t0,ran,xan,xu,rd,da,ngeas_max,ndf,ilin)
c....   loop over columns (symmetry noted), stiffness matrix
        jc0 = 0       
        do jnode = 1,inode
         call stif63 (inode,jnode,ir0,jc0,b,btd,st,sig,qs,shp,shp1,
     1                w,ngeas_max,ndf,da,ilin,isw)
         jc0 = jc0 + ndf
        end do
2210     ir0 = ir0 + ndf
2200   continue

c
      do i = 1,4*ndf
       p(i) = pt(i)
       do j = 1,i
        s(i,j) = st(i,j)
        s(j,i) = st(i,j)
       enddo
      enddo   

      endif
c
c     loads: conservative and follower forces
      if(iflg1.eq.1)then
       ngb   = 2
       call pgauss (ngb,lint,sg,tg,wg)
       call rjac63(d,xl,t0,xs0,sx0,ts0,te0,detj0,isw)
       do l = 1,lint
        xsi = sg(l)
        eta = tg(l)
        call shap63(xsi,eta,xl,t0,sx,shp,shp1,xsj)
        da = xsj*wg(l)
        call qload63
     1       (shp1,xl,ul,t0,aqloa,numel,n,mqloa,propq,ndf,nst,p,s,da)
       enddo
      endif
c       
      return
      end
c
      subroutine testq63(d,q,numel,n,iflg1,iflg2)
      implicit double precision (a-h,o-z)
      dimension d(*),q(numel,10)
c
c.... check load vector (follower forces ifg1=1, temperature loading iflg2=1)
c
      iflg1 = 0
      iflg2 = 0
      if(dabs(q(n,1)).gt.1.d-20)iflg1 = 1
      if(dabs(q(n,6))+dabs(q(n,7)).gt.1.d-20)iflg2 = 1
c
      return
      end
c
      subroutine qload63(shp1,xl,ul,t0,q,numel,n,mqloa,propq,
     1                  ndf,nst,p,s,da)
      USE prlod
      USE  iofile
      implicit double precision (a-h,o-z)
c      USE prlod
c      USE iofile
      dimension q(numel,10)
      dimension shp1(3,*),xl(3,*),ul(ndf,*),xu1(3,3),p(*),s(nst,*),
     1          t0(3,3),ql(3)
c
c     q(n,1) = element n, load q_n (transverse) 
c     q(n,2) = element n, load q_x (global) 
c     q(n,3) = element n, load q_y (global)
c     q(n,4) = element n, load q_z (global)
c     q(n,5) = element n, iltyp, conservative = 0,  follower load  = 1
c     
      if(mqloa.eq.1)return
      iltyp = q(n,5)
c     
c.... conservative loads
c
      if(iltyp.eq.0)then 
        ir0 = 0
        do inode = 1,4
          do i = 1,3
            ql(i) = t0(i,3)*q(n,1)
          end do
          qx = q(n,2)
          qy = q(n,3)
          qz = q(n,4)
          ql(1) = ql(1) + qx
          ql(2) = ql(2) + qy
          ql(3) = ql(3) + qz
          fact = propq * shp1(3,inode) * da
          do i = 1,3
           p(ir0+i) = p(ir0+i) + fact*ql(i)
          end do
        ir0 = ir0 + ndf
        end do
      endif
c
c.... follower loads
c
      if(iltyp.eq.1)then 
        
        do i = 1,3
         do j = 1,2
          xu1(i,j)  = 0.d0
           do k = 1,4
            xu1(i,j) = xu1(i,j) + (xl(i,k)+ul(i,k))*shp1(j,k)
           end do
         end do 
        end do
c       
        call vecp (xu1(1,1),xu1(1,2),xu1(1,3))
c       
        ir0 = 0
        do inode = 1,4
          fact = 0.d0
          fact=propq*q(n,1)*shp1(3,inode)
         do i = 1,3
          p(ir0+i) = p(ir0+i) + xu1(i,3)*fact
         end do 
        
         jc0 = 0
         do jnode = 1,4
          s1 = (xu1(1,1)*shp1(2,jnode) - xu1(1,2)*shp1(1,jnode))*fact
          s2 = (xu1(2,1)*shp1(2,jnode) - xu1(2,2)*shp1(1,jnode))*fact
          s3 = (xu1(3,1)*shp1(2,jnode) - xu1(3,2)*shp1(1,jnode))*fact
          s(ir0+2,jc0+3) = s(ir0+2,jc0+3) + s1
          s(ir0+3,jc0+1) = s(ir0+3,jc0+1) + s2
          s(ir0+1,jc0+2) = s(ir0+1,jc0+2) + s3
          s(ir0+3,jc0+2) = s(ir0+3,jc0+2) - s1
          s(ir0+1,jc0+3) = s(ir0+1,jc0+3) - s2
          s(ir0+2,jc0+1) = s(ir0+2,jc0+1) - s3
          jc0 = jc0 + ndf
         end do
        ir0 = ir0 + ndf
        end do
      endif 
c
      return
      end
c----------------------------------------------------------------------
      subroutine rjac63(d,xl,t0,xs0,sx0,ts0,te0,detj0,isw,ni,nj,
     +                  NURp,NURq,ndm,NURknv1,NURknv2)
c----------------------------------------------------------------------
      USE eldata
      USE isogeo
      implicit double precision (a-h,o-z)
      real*8 NURknv1(NURlenkv(1)),NURknv2(NURlenkv(2))
c
c.... compute jacobian matrix at element center
c
      dimension    d(*),
     1             xs0(2,2),sx0(2,2), xl(ndm,*),te0(5,*),ts0(5,*),
     2             g0(3,3),t0(3,3)            ,gs(3,3),x(3)
      dimension dR_dXi(2,(NURp+1)*(NURq+1)), shp(3,(NURp+1)*(NURq+1))
      real*8 N_Xi1(2,NURp+1),N_Xi2(2,NURq+1)
c      include 'eldata.h'
c      include 'isogeo.h'
      call pzero(g0,9)
c     calculate NURBS parametric coordinates at element center
      Xi1=0.5d0*(NURknv1(ni+1) + NURknv1(ni))
      Xi2=0.5d0*(NURknv2(nj+1) + NURknv2(nj))
c
c     compute univariate basis functions and their derivatives
      call DersBasisFuns(ni,Xi1,NURp,NURknv1,N_Xi1)
      call DersBasisFuns(nj,Xi2,NURq,NURknv2,N_Xi2)
c
      W=0;
      W_Xi1=0;
      W_Xi2=0;
      int=0
      do j=0,NURq
        do i=0,NURp
            int=int+1
            W=W+N_Xi1(1,NURp+1-i)*N_Xi2(1,NURq+1-j)*xl(4,int)
            W_Xi1=W_Xi1+N_Xi1(2,NURp+1-i)*N_Xi2(1,NURq+1-j)*xl(4,int)
            W_Xi2=W_Xi2+N_Xi1(1,NURp+1-i)*N_Xi2(2,NURq+1-j)*xl(4,int)
        end do
      end do
      int=0;
      do j=0,NURq
        do i=0,NURp
            int=int+1
            shp(3,int)=N_Xi1(1,NURp+1-i)*N_Xi2(1,NURq+1-j)
     +                 *xl(4,int)/W
            dR_dXi(1,int)=xl(4,int)*(N_Xi1(2,NURp+1-i)
     +                 *N_Xi2(1,NURq+1-j)*W-W_Xi1*N_Xi1(1,NURp+1-i)
     +                 *N_Xi2(1,NURq+1-j))/(W**2)
            dR_dXi(2,int)=xl(4,int)*(N_Xi1(1,NURp+1-i)
     +                 *N_Xi2(2,NURq+1-j)*W-W_Xi2*N_Xi1(1,NURp+1-i)
     +                 *N_Xi2(1,NURq+1-j))/(W**2)
        end do
      end do
c     gradient of mapping from parameter space dXi to physical space dx
c     dx_dXi(:,i) contains g_i (convective base vectors)
      x=0
      dx_dXi=0
      loc_num=0
      do j=0,NURq
        do i=0,NURp
            loc_num=loc_num+1
            do ia=1,3
                do ib=1,2
                    g0(ia,ib)=g0(ia,ib)+xl(ia,loc_num)
     +               *dR_dXi(ib,loc_num)         
                end do
                x(ia)=x(ia)+xl(ia,loc_num)*shp(3,loc_num)
            end do
        end do
      end do
c     construct g_3 from g_1 and g_2
      call vecp (g0(1,1),g0(1,2),g0(1,3))
c
c.... t0 ------ Taylor
!      do i = 1,3
!       t0(i,1) = xl(i,3) - xl(i,1)
!       t0(i,2) = xl(i,2) - xl(i,4)
!      end do
!      call norm (t0(1,1),t0(1,1),3)
!      call norm (t0(1,2),t0(1,2),3)
!      do i = 1,3
!       v1 = t0(i,1)
!       v2 = t0(i,2)
!       t0(i,1) = v1 + v2
!       t0(i,2) = v1 - v2
!      end do
!      call norm (t0(1,1),t0(1,1),3)
!      call norm (t0(1,2),t0(1,2),3)
!      call vecp (t0(1,1),t0(1,2),t0(1,3))
c.... t0-----  Hughes  lamina basis, Belytschko 541 
      call norm (t0(1,3),g0(1,3),3)
      call norm (gs(1,1),g0(1,1),3)
      call norm (gs(1,2),g0(1,2),3)
      do i = 1,3
       gs(i,1) = gs(i,1) + gs(i,2)
      end do
      call vecp (t0(1,3),gs(1,1),gs(1,2))
      do i = 1,3
       t0(i,1) = gs(i,1) - gs(i,2)
       t0(i,2) = gs(i,1) + gs(i,2)
      end do
      call norm (t0(1,1),t0(1,1),3)
      call norm (t0(1,2),t0(1,2),3)
c.... t0-----  FG
!      call norm (t0(1,3),g0(1,3),3)
!      call norm (t0(1,1),g0(1,1),3)
!      call vecp (t0(1,3),t0(1,1),t0(1,2))

      if(isw.eq.5) return

c.... construct jacobian and its inverse
      do i = 1,2
        do j = 1,2
         xs0(i,j) = dot(g0(1,i),t0(1,j),3)
        end do
      end do
      detj0 = xs0(1,1)*xs0(2,2)-xs0(1,2)*xs0(2,1)
      if(detj0.le.0.d0) write(*,*) 'detj negative in element  ', n
      sx0(1,1) = xs0(2,2)/detj0
      sx0(2,2) = xs0(1,1)/detj0
      sx0(1,2) =-xs0(1,2)/detj0
      sx0(2,1) =-xs0(2,1)/detj0
c
      call pzero(ts0,5*5)
      call pzero(te0,5*5)
      ts0(1,1) = xs0(1,1)*xs0(1,1)  
      ts0(2,1) = xs0(1,2)*xs0(1,2)  
      ts0(3,1) = xs0(1,1)*xs0(1,2)  
      ts0(1,2) = xs0(2,1)*xs0(2,1)  
      ts0(2,2) = xs0(2,2)*xs0(2,2)  
      ts0(3,2) = xs0(2,1)*xs0(2,2)  
      ts0(1,3) = 2.d0*xs0(1,1)*xs0(2,1)  
      ts0(2,3) = 2.d0*xs0(1,2)*xs0(2,2)  
      ts0(3,3) = xs0(1,1)*xs0(2,2) + xs0(1,2)*xs0(2,1)
      ts0(4,4) = xs0(1,1) 
      ts0(5,4) = xs0(1,2) 
      ts0(4,5) = xs0(2,1) 
      ts0(5,5) = xs0(2,2) 

      call invert (ts0,5,5)
      do i = 1,5
       do j = 1,5
        te0(i,j) = ts0(j,i)
       enddo 
      enddo 
c
      return
      end

      subroutine shap63onlyBasisFun(xsi,eta,xl,shp,ndm,nel,ni,nj,NURp
     +                 ,NURq,NURknv1,NURknv2)
      USE isogeo
      USE iofile
c----------------------------------------------------------------------
c
c      Purpose: Computes shape function and derivatives for
c               isogeometric surface elements
c
c      Inputs:
c         xsi        - Natural coordinates for point (unit element coordinates)
c         eta        - Natural coordinates for point (unit element coordinates)
c         x(ndm,*)  - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh
c         nel       - Number of nodes on element
c         ix(*)     - Nodes attached to element
c         flg       - Flag, compute global x/y derivatives if false,
c                           else derivatives are w/r natural coords.
c         ni        - number of highest non-zero basis function in 1-direction
c         nj        - number of highest non-zero basis function in 2-direction
c
c      Outputs:
c         shp(*)  - Shape functions and derivatives at point
c                     shp(1,i) =  N_i
c
c---------------------------------------------------------------------- 
      implicit double precision (a-h,o-z)
      dimension shp((NURp+1)*(NURq+1))
      dimension xl(ndm,nel)
      real*8 N_Xi1(2,NURp+1),N_Xi2(2,NURq+1)
      real*8 NURknv1(NURlenkv(1)),NURknv2(NURlenkv(2))
c      include 'isogeo.h'
c      include 'iofile.h'

c      
c     calculate NURBS parametric coordinates from unit element coordinates
      Xi1=0.5d0*((NURknv1(ni+1)-NURknv1(ni))*xsi+(NURknv1(ni+1)+
     + NURknv1(ni)))
      Xi2=0.5d0*((NURknv2(nj+1)-NURknv2(nj))*eta+(NURknv2(nj+1)+
     + NURknv2(nj)))
c
c     compute univariate basis functions and their derivatives
      call DersBasisFuns(ni,Xi1,NURp,NURknv1,N_Xi1)
      call DersBasisFuns(nj,Xi2,NURq,NURknv2,N_Xi2)
c
      W=0;
      int=0
      do j=0,NURq
        do i=0,NURp
            int=int+1
            W=W+N_Xi1(1,NURp+1-i)*N_Xi2(1,NURq+1-j)*xl(4,int)
        end do
      end do
      int=0;
      do j=0,NURq
        do i=0,NURp
            int=int+1
            shp(int)=N_Xi1(1,NURp+1-i)*N_Xi2(1,NURq+1-j)
     +                 *xl(4,int)/W
        end do
      end do

      return
      end
c
      subroutine shap63(xsi,eta,xl,t0,sx,shp,shp1,xsj,ndm,nel,ni,nj,NURp
     +                 ,NURq,NURknv1,NURknv2)
      USE isogeo
      USE iofile
c----------------------------------------------------------------------
c
c      Purpose: Computes shape function and derivatives for
c               isogeometric surface elements
c
c      Inputs:
c         xsi        - Natural coordinates for point (unit element coordinates)
c         eta        - Natural coordinates for point (unit element coordinates)
c         x(ndm,*)  - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh
c         nel       - Number of nodes on element
c         ix(*)     - Nodes attached to element
c         flg       - Flag, compute global x/y derivatives if false,
c                           else derivatives are w/r natural coords.
c         ni        - number of highest non-zero basis function in 1-direction
c         nj        - number of highest non-zero basis function in 2-direction
c
c      Outputs:
c         shp(3,*)  - Shape functions and derivatives at point
c                     shp(1,i) = dN_i/dx or dN_i/dxi_1
c                     shp(2,i) = dN_i/dy or dN_i/dxi_2
c                     shp(3,i) = N_i
c         xsj       - Jacobian determinant at point
c
c---------------------------------------------------------------------- 
      implicit double precision (a-h,o-z)
      dimension shp(3,(NURp+1)*(NURq+1)),dR_dXi(2,(NURp+1)*(NURq+1))
      dimension xl(ndm,nel)
      dimension dx_dXi(3,3)
      real*8 N_Xi1(2,NURp+1),N_Xi2(2,NURq+1)
      real*8 NURknv1(NURlenkv(1)),NURknv2(NURlenkv(2))
c      include 'isogeo.h'
c      include 'iofile.h'
      dimension    xs(2,2),sx(2,2),shp1(3,4),
     2		   gs(3,3),t0(3,3),t0gp(3,3)
      logical uset0local
      call pzero(dx_dXi,9)
c      
c     calculate NURBS parametric coordinates from unit element coordinates
      Xi1=0.5d0*((NURknv1(ni+1)-NURknv1(ni))*xsi+(NURknv1(ni+1)+
     + NURknv1(ni)))
      Xi2=0.5d0*((NURknv2(nj+1)-NURknv2(nj))*eta+(NURknv2(nj+1)+
     + NURknv2(nj)))
c
c     compute univariate basis functions and their derivatives
      call DersBasisFuns(ni,Xi1,NURp,NURknv1,N_Xi1)
      call DersBasisFuns(nj,Xi2,NURq,NURknv2,N_Xi2)
c
      W=0;
      W_Xi1=0;
      W_Xi2=0;
      int=0
      do j=0,NURq
        do i=0,NURp
            int=int+1
            W=W+N_Xi1(1,NURp+1-i)*N_Xi2(1,NURq+1-j)*xl(4,int)
            W_Xi1=W_Xi1+N_Xi1(2,NURp+1-i)*N_Xi2(1,NURq+1-j)*xl(4,int)
            W_Xi2=W_Xi2+N_Xi1(1,NURp+1-i)*N_Xi2(2,NURq+1-j)*xl(4,int)
        end do
      end do
      int=0;
      sum_shape=0.0d0
      do j=0,NURq
        do i=0,NURp
            int=int+1
            shp(3,int)=N_Xi1(1,NURp+1-i)*N_Xi2(1,NURq+1-j)
     +                 *xl(4,int)/W
            dR_dXi(1,int)=xl(4,int)*(N_Xi1(2,NURp+1-i)
     +                 *N_Xi2(1,NURq+1-j)*W-W_Xi1*N_Xi1(1,NURp+1-i)
     +                 *N_Xi2(1,NURq+1-j))/(W**2)
            dR_dXi(2,int)=xl(4,int)*(N_Xi1(1,NURp+1-i)
     +                 *N_Xi2(2,NURq+1-j)*W-W_Xi2*N_Xi1(1,NURp+1-i)
     +                 *N_Xi2(1,NURq+1-j))/(W**2)
            sum_shape=sum_shape+shp(3,int)
        end do
      end do

c      if (abs(sum_shape-1.0d0).ge.1.d-8) write(*,*) sum_shape
c     gradient of mapping from parameter space dXi to physical space dx
c     dx_dXi(:,i) contains g_i
c      dx_dXi=0
      loc_num=0
      do j=0,NURq
        do i=0,NURp
            loc_num=loc_num+1
            do ia=1,3
                do ib=1,2
                    dx_dXi(ia,ib)=dx_dXi(ia,ib)+xl(ia,loc_num)
     +               *dR_dXi(ib,loc_num)         
                end do
            end do
        end do
      end do
c     construct g_3 from g_1 and g_2
      call vecp (dx_dXi(1,1),dx_dXi(1,2),dx_dXi(1,3))
c.... area element
      xsj = dsqrt(dot(dx_dXi(1,3),dx_dXi(1,3),3))
     +      *(NURknv1(ni+1)-NURknv1(ni))/2
     +      *(NURknv2(nj+1)-NURknv2(nj))/2
c.... test: one additional transformation
c      dx_dXi(1:3,1) = dx_dXi(1:3,1)*(NURknv1(ni+1)-NURknv1(ni))/2
c      dx_dXi(1:3,2) = dx_dXi(1:3,2)*(NURknv2(nj+1)-NURknv2(nj))/2
c.... calculating t0 in every GP seems the best way for IGA
      uset0local = .true.    !absolutely needed for perpendicular element loads
      if (uset0local) then
c       compute lamina basis system (local cartesian), see Hughes, FEM, 2000, p. 385
        call norm (t0(1,3),dx_dXi(1,3),3)
        call norm (gs(1,1),dx_dXi(1,1),3)
        call norm (gs(1,2),dx_dXi(1,2),3)
        do i = 1,3
            gs(i,1) = gs(i,1) + gs(i,2)
        end do
        call vecp (t0(1,3),gs(1,1),gs(1,2))
        do i = 1,3
            t0(i,1) = gs(i,1) - gs(i,2)
            t0(i,2) = gs(i,1) + gs(i,2)
        end do
        call norm (t0(1,1),t0(1,1),3)
        call norm (t0(1,2),t0(1,2),3)
        !t0(1:3,1:3) = dx_dXi(1:3,1:3)
        !call plamina_base(t0)
      end if
c.... construct jacobian and its inverse
      do i = 1,2
        do j = 1,2
c         original laut Wagner/Gruttmann
          xs(i,j) = dot(dx_dXi(1,i),t0(1,j),3)
c         test          
          !xs(i,j) = dot(dx_dXi(1,j),t0(1,i),3)
          !sx(i,j) = dot(dx_dXi(1,i),t0(1,j),3)
        end do
      end do
      detj = xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
      sx(1,1) = xs(2,2)/detj
      sx(2,2) = xs(1,1)/detj
      sx(1,2) =-xs(1,2)/detj
      sx(2,1) =-xs(2,1)/detj
c.... form local derivatives
      do i = 1,(NURq+1)*(NURp+1)
        shp(1,i) = sx(1,1)*dR_dXi(1,i) + sx(1,2)*dR_dXi(2,i)
        shp(2,i) = sx(2,1)*dR_dXi(1,i) + sx(2,2)*dR_dXi(2,i)
      end do
      return
      end
c      
      subroutine easstr63 (xsi,eta,te0,rj,neas_max,neas,nm,nb,ns,
     1                     gxy,eas,alpha)
c-----------------------------------------------------------------------
c.... enhanced strain interpolation
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension gxy(8,*),te0(5,*),gxsi(3,7),g(3,7),eas(*),alpha(*)
c
      call pzero (gxsi,3*7)
      call pzero (gxy,8*neas_max)

c.... setup gxsi for max 7 parameters     
       gxsi(1,1) = xsi*rj
       gxsi(2,2) = eta*rj
       gxsi(3,3) = xsi*rj
       gxsi(3,4) = eta*rj
       gxsi(3,5) = xsi*eta*rj
       gxsi(1,6) = xsi*eta*rj
       gxsi(2,7) = xsi*eta*rj
c
      do i = 1,3
       do j = 1,7
         g(i,j) = 0.d0
        do k = 1,3
         g(i,j) =  g(i,j) + te0(i,k)*gxsi(k,j)
        enddo
       enddo
      enddo 
c.... membrane part               
      do i = 1,3
       do j = 1,nm
        gxy(i,j) = g(i,j)
       enddo
      enddo 
c.... bending part               
      do i = 1,3
       do j = 1,nb
        gxy(3+i,nm+j) = g(i,j)
       enddo
      enddo 
c.... shear part               
c
c.... enhanced strains
      do i = 1,8
       eas(i) = 0.d0
       do j = 1,neas
        eas(i) = eas(i) + gxy(i,j)*alpha(j)
       enddo
      enddo  
c
      return
      end
c
      subroutine easmat63 (gxy,gtd,pt,st,sig,dmat,b,da,neas,
     1                     neas_max,ngeas_max,ndf,inode,ir0,iflg)
c-----------------------------------------------------------------------
c.... matrix M~^T*D
c     residual -f~= int( M~^T *S )
c     stiffness H = int( M~^T *D* M~ ) (only left lower part!)
c     G = int( M~^T *D* B )
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension gxy(8,*),gtd(neas_max,*),pt(*),st(ngeas_max,*),
     1          sig(*),dmat(8,8),b(8,6,*)

      if (iflg.eq.2) goto 2
c.... residual -f~ and stiffness part M~^T *D
      do 10 i = 1,neas
       do 10 k = 1,8
         pt(4*ndf+i) = pt(4*ndf+i) - gxy(k,i)*sig(k)*da
10       continue

      do 15 i = 1,neas
       do 15 k = 1,8
        gtd(i,k) = 0.0d0
        do 15 j = 1,8
15       gtd(i,k) = gtd(i,k) + gxy(j,i)*dmat(j,k)*da

c.... stiffness H = int( M~^T *D* M~ )   lower part
      do 20 i = 1,neas
       do 20 j = 1,i
        do 20 k = 1,8
20       st(4*ndf+i,4*ndf+j) = st(4*ndf+i,4*ndf+j) + gtd(i,k)*gxy(k,j)
      return

c.... stiffness G = int( M~^T *D* B )
2     continue
      do 30 l = 1,neas
       do 30 i = 1,ndf
        do 30 k = 1,8
30       st(4*ndf+l,ir0+i) = st(4*ndf+l,ir0+i) + gtd(l,k)*b(k,i,inode)
      return
      end
c
      subroutine easend63(st,pt,s,p,skcr,prc,nst,ndf,ngeas_max,neas,nen)
c-----------------------------------------------------------------------
c.... static condensation and store matrices for back substitution
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension st(ngeas_max,*),pt(*),s(nst,*),p(*),skcr(neas,*),prc(*)

c.... static condensation
      if(neas.gt.0) call conden(st,pt,4*ndf,neas,ngeas_max,.false.)

c.... copy stiffness matrix st to s and pt to p 
      do i = 1,nen*ndf
       p(i) = pt(i)
       do j = 1,i
        s(i,j) = st(i,j)
        s(j,i) = st(i,j)
       enddo
      enddo   
       
      if(neas.eq.0)return

c.... store matrices for back substitution
      do 30 i = 1,neas
        ii = 4*ndf+i
        prc(i) = pt(ii)
        do 30 k = 1,4*ndf+neas
30        skcr(i,k) = st(ii,k)
      return
      end
c
      subroutine bmat63 (inode,ir0,b,btd,dmat,st,pt,sig,qs,rwn,dwn,w,
     1     twn,bml,d,shp1,shp,sx,t0,ran,xan,xu,rd,da,ngeas_max,ndf,ilin,
     2     noans,x0,useTgp,igauss,w_1,w_2,isymm,jnode)
      USE prlod
      USE eldata
      USE isogeo
c-----------------------------------------------------------------------
c.... b-matrix, residual vector, external load vector  
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      include 'prlod.h'
c      include 'eldata.h'
cwd   temporarily include for saving forcx and forcy
c      include 'isogeo.h'
      dimension d(*),shp1(3,*),shp(3,*),pt(*),st(ngeas_max,*),b(8,6,*),
     1          btd(8,*),dmat(8,8),xu(3,*),sx(2,2),rd(3,*),t0(3,*),
     2          rwn(3,3,*),fbb(3,2),ql(3),x0(3,3),
     3          sig(*),qs(*),wmh(3,3),xmh(3),dwn(3,*),ran(3,*),
     4          bml(3,2,*),twn(3,3,*),xan(3,*),w(3,3,*),wh(3,3),
     5          w_1(3,3,*),w_2(3,3,*) 
      dimension im(4),il(4),fbb11(3),fbb12(3),fbb21(3),
     1          fbb22(3)
      real*8 du(3,3),du1(3),du2(3)
      integer useTgp
      data im /-2, 2, 4,-4/, il /-1,-3, 3, 1/
c
      fl=0.d0;fm=0.d0;
c
      if (noans.eq.0) then
        mm = iabs(im(inode))
        ll = iabs(il(inode))
        fm = isign(1,im(inode))*shp1(1,inode)
        fl = isign(1,il(inode))*shp1(2,inode)
      end if
      call mttmul (w(1,1,inode),xu,3,3,2, fbb)
      if (useTgp.ge.2) then
        call mttmul (w_1(1,1,inode),xu(1,1),3,3,1,fbb11)
        call mttmul (w_1(1,1,inode),xu(1,2),3,3,1,fbb12)
        call mttmul (w_2(1,1,inode),xu(1,1),3,3,1,fbb21)
        call mttmul (w_2(1,1,inode),xu(1,2),3,3,1,fbb22)
      end if

c
c.... b-matrix: membrane (1,2,3), bending (4,5,6), shear (7,8) 
      do i = 1,3
c        b1 = ran(i,mm)*shp1(1,inode)
c        b2 = ran(i,ll)*shp1(2,inode)
        b(1,i,inode) = xu(i,1)*shp(1,inode)
        b(2,i,inode) = xu(i,2)*shp(2,inode)
        b(3,i,inode) = xu(i,1)*shp(2,inode) + xu(i,2)*shp(1,inode)
        b(4,i,inode) = rd(i,1)*shp(1,inode)
        b(5,i,inode) = rd(i,2)*shp(2,inode)
        b(6,i,inode) = rd(i,1)*shp(2,inode) + rd(i,2)*shp(1,inode)
        if (noans.eq.0) then
            b(7,i,inode) = sx(1,1)*b1 + sx(1,2)*b2
            b(8,i,inode) = sx(2,1)*b1 + sx(2,2)*b2
        else
            b(7,i,inode) = rd(i,3)*shp(1,inode)
            b(8,i,inode) = rd(i,3)*shp(2,inode)
        end if
      end do  

      do i = 1,ndf-3
        i3 = i+3
        b3 = fm*bml(i,1,inode)
        b4 = fl*bml(i,2,inode)
        b(1,i3,inode) = 0.0d0
        b(2,i3,inode) = 0.0d0
        b(3,i3,inode) = 0.0d0
        if (useTgp.lt.2) then
          b(4,i3,inode) = shp(1,inode)*fbb(i,1)
          b(5,i3,inode) = shp(2,inode)*fbb(i,2)
          b(6,i3,inode) = shp(1,inode)*fbb(i,2) + shp(2,inode)*fbb(i,1)
        else
          b(4,i3,inode) = fbb11(i)
          b(5,i3,inode) = fbb22(i)
          b(6,i3,inode) = fbb12(i) + fbb21(i)
        end if
!        temp = shp(1,inode)*fbb(i,1)-fbb11(i)
!        if(temp.gt.1d-10) then
!          write(*,*) temp
!        end if
!        temp = shp(2,inode)*fbb(i,2)-fbb22(i)
!        if(temp.gt.1d-10) then
!          write(*,*) temp
!        end if
        if (noans.eq.0) then
          b(7,i3,inode) = sx(1,1)*b3 + sx(1,2)*b4
          b(8,i3,inode) = sx(2,1)*b3 + sx(2,2)*b4
        else
          b(7,i3,inode) = shp(3,inode)*fbb(i,1)
          b(8,i3,inode) = shp(3,inode)*fbb(i,2)
        end if
      end do

c...  compute B(jnode) in case it is not yet computed
      if (isymm.eq.1.and.jnode.gt.inode) then
        if (noans.eq.0) then
          mm = iabs(im(jnode))
          ll = iabs(il(jnode))
          fm = isign(1,im(jnode))*shp1(1,jnode)
          fl = isign(1,il(jnode))*shp1(2,jnode)
        end if
        call mttmul (w(1,1,jnode),xu,3,3,2, fbb)
        if (useTgp.ge.2) then
          call mttmul (w_1(1,1,jnode),xu(1,1),3,3,1,fbb11)
          call mttmul (w_1(1,1,jnode),xu(1,2),3,3,1,fbb12)
          call mttmul (w_2(1,1,jnode),xu(1,1),3,3,1,fbb21)
          call mttmul (w_2(1,1,jnode),xu(1,2),3,3,1,fbb22)
        end if

c
c.... b-matrix: membrane (1,2,3), bending (4,5,6), shear (7,8) 
        do i = 1,3
c        b1 = ran(i,mm)*shp1(1,inode)
c        b2 = ran(i,ll)*shp1(2,inode)
          b(1,i,jnode) = xu(i,1)*shp(1,jnode)
          b(2,i,jnode) = xu(i,2)*shp(2,jnode)
          b(3,i,jnode) = xu(i,1)*shp(2,jnode) + xu(i,2)*shp(1,jnode)
          b(4,i,jnode) = rd(i,1)*shp(1,jnode)
          b(5,i,jnode) = rd(i,2)*shp(2,jnode)
          b(6,i,jnode) = rd(i,1)*shp(2,jnode) + rd(i,2)*shp(1,jnode)
          if (noans.eq.0) then
            b(7,i,jnode) = sx(1,1)*b1 + sx(1,2)*b2
            b(8,i,jnode) = sx(2,1)*b1 + sx(2,2)*b2
          else
            b(7,i,jnode) = rd(i,3)*shp(1,jnode)
            b(8,i,jnode) = rd(i,3)*shp(2,jnode)
          end if
        end do  

        do i = 1,ndf-3
          i3 = i+3
          b3 = fm*bml(i,1,jnode)
          b4 = fl*bml(i,2,jnode)
          b(1,i3,jnode) = 0.0d0
          b(2,i3,jnode) = 0.0d0
          b(3,i3,jnode) = 0.0d0
          if (useTgp.lt.2) then
            b(4,i3,jnode) = shp(1,jnode)*fbb(i,1)
            b(5,i3,jnode) = shp(2,jnode)*fbb(i,2)
            b(6,i3,jnode) = shp(1,jnode)*fbb(i,2)+ shp(2,jnode)*fbb(i,1)
          else
            b(4,i3,jnode) = fbb11(i)
            b(5,i3,jnode) = fbb22(i)
            b(6,i3,jnode) = fbb12(i) + fbb21(i)
          end if

          if (noans.eq.0) then
            b(7,i3,jnode) = sx(1,1)*b3 + sx(1,2)*b4
            b(8,i3,jnode) = sx(2,1)*b3 + sx(2,2)*b4
          else
            b(7,i3,jnode) = shp(3,jnode)*fbb(i,1)
            b(8,i3,jnode) = shp(3,jnode)*fbb(i,2)
          end if
        end do
      end if
      
      
c.... residual vector and  B_T * D
      do i = 1,ndf
        ir = ir0 + i
        if (jnode.eq.1) then
        pt(ir) = pt(ir) - dot(sig,b(1,i,inode),8)*da
        end if
c        if (isnan(pt(ir))) then
c          write(*,*) 'errrroorrr'
c        end if
       do j = 1,8
        btd(j,i) = dot(b(1,i,inode),dmat(1,j),8)*da
       enddo
      enddo 
c
c.... compute element load vector
      if(dot(d(25),d(25),4).gt.1.d-20) then
c..... external constant load perpendicular to element
        do i = 1,3
          ql(i) = t0(i,3)*d(25)
        end do
c..... external constant loads qx, qy, qz in x,y,z-direction
        qx = d(26)
        qy = d(27)
        qz = d(28)
!cwd     exact load: for tests
!        x = x0(1,3)
!        y = x0(2,3)
!        qz = d(40)*d(13)**3/(12.0d0*(1.0d0-d(41)*d(41)))*
!     *     (12.d0*y*(y-1.d0)*(5.d0*x*x-5.d0*x+1.d0)*
!     *     (2.d0*y*y*(y-1.d0)**2+x*(x-1.d0)*
!     *     (5.d0*y*y-5.d0*y+1.d0))+12.d0*x*(x-1.d0)*
!     *     (5.d0*y*y-5.d0*y+1.d0)*(2.d0*x*x*(x-1.d0)**2+
!     +     y*(y-1.d0)*(5.d0*x*x-5.d0*x+1.d0)))
!cwd   end TEST        
        ql(1) = ql(1) + qx
        ql(2) = ql(2) + qy
        ql(3) = ql(3) + qz
c..... add loads to load vector
        if (jnode.eq.1) then
        do i = 1,3
         pt(ir0+i) = pt(ir0+i) + prop*ql(i) * shp(3,inode)*da
        end do
        end if
       end if
c
c....  transform q to convective coordinate system
       qs(1) = sx(1,1)*sig(7) + sx(2,1)*sig(8)
       qs(2) = sx(1,2)*sig(7) + sx(2,2)*sig(8)

      if(ilin.eq.0.or.ilin.eq.1) return
      if(useTgp.ge.2) return
c
c....  geometric stiffness  (diagonal terms)
c      if(dot(dwn(1,inode),dwn(1,inode),3).gt.1.d-15)then 
       do i = 1,3
         !xmh is G_ßß from Gruttmann-script
         if (noans.eq.0) then
         xmh(i) = ( sig(4)*xu(i,1)*shp(1,inode) 
     1            + sig(5)*xu(i,2)*shp(2,inode)
     2            + sig(6)*xu(i,2)*shp(1,inode) 
     3            + sig(6)*xu(i,1)*shp(2,inode)
     4            + qs(1)*fm*xan(i,mm)+qs(2)*fl*xan(i,ll))*da
          else
          xmh(i) = ( sig(4)*xu(i,1)*shp(1,inode) 
     1            + sig(5)*xu(i,2)*shp(2,inode)
     2            + sig(6)*xu(i,2)*shp(1,inode) 
     3            + sig(6)*xu(i,1)*shp(2,inode)
     4            + sig(7)*xu(i,1)*shp(3,inode)
     5            + sig(8)*xu(i,2)*shp(3,inode))*da
          end if

       end do

       call faca63 (rwn(1,3,inode),xmh,dwn(1,inode), wmh)
       call mttmul (twn(1,1,inode),wmh,3,3,3, wh)
       call matmulf (wh,twn(1,1,inode),3,3,3, wmh)

       if (jnode.eq.1) then
       do i = 1,ndf-3
         ir2 = ir0+3+i
         do j = 1,i
           ic2 = ir0+3+j
           st(ir2,ic2) = st(ir2,ic2) + wmh(i,j)
         end do
       end do  
       end if

c
      return
      end

c
      subroutine boun63 (inode,ix,id,iflag,ndf)
c-----------------------------------------------------------------------
c.... bei glattem Elementuebergang: RB bei 6.Knotenfhg gesetzt,
c     d.h.  5dof         -> iflag = 1
c     sonst 6dof (Kante) -> iflag = 0
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension ix(*),id(ndf,*)

c.... Globale Knotennummer
      iglob = ix(inode)
      iflag = 0
      if(id(6,iglob).lt.0) iflag = 1
      return
      end
c
      subroutine shearfac63(d,detj0,eps,sig,dmat,icappa,cappa,wcappa)
c-----------------------------------------------------------------------
c.... shear correction factor
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  d(*),eps(*),sig(*),dmat(8,*),yl(3)

        if(icappa.eq.0) then
          cappa  =  5.0d0/6.d0 
cww       cappa  =  1000.d0   
        else  

          imat   = d(15) 
          xnu    = d(41)
          hs     = d(13)
          if(imat.eq.3) then
            gmod = 0.5d0*(d(40)*d(43)+d(41)*d(44)+d(42)*d(45))
            rlam = d(46) 
            xnu  = rlam/(2.d0*(rlam+gmod))
          end if       

c....     AE... element area
          ae = 4.d0*detj0 
       
          cappa  = (5.d0/6.d0)/(1.d0+5.d0*ae/(16.d0*hs*hs*(1.d0+xnu)))
        end if 
        wcappa = dsqrt(cappa) 
c      
      return
      end  
c

      subroutine dmate63 (d,cappa,dmat)
c-----------------------------------------------------------------------
c     d(13) = shell thickness h_s
c     d(14) = distance bottom surface to reference surface 
c     d(40) = elasticity modulus
c     d(41) = poisson's ratio
c     cmat  = elasticity matrix
c     dmat  = elasticity matrix for stress resultants
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension d(*),cmat(3,3),dmat(8,*)
c
      call pzero(cmat,3*3)
      call pzero(dmat,8*8)

      cmat(1,1) = d(40)/(1.0d0-d(41)*d(41))                  
      cmat(2,1) = cmat(1,1)*d(41)                               
      cmat(1,2) = cmat(2,1)                               
      cmat(2,2) = cmat(1,1)
      cmat(3,3) = cmat(1,1)*(1.0d0-d(41))*0.5d0                 
c      
      h  = d(13)
      cm = h
      zs = d(14) + 0.5d0*h
      cb = h**3/12.d0 + h*zs*zs  
      cmb= h*zs
      cs = h * cappa
c
c.... membrane/bending/coupling
      do i=1,3
        do j=1,3
          dmat(i  ,j  ) = cmat(i,j)*cm
          dmat(i+3,j  ) = cmat(i,j)*cmb
          dmat(i  ,j+3) = cmat(i,j)*cmb
          dmat(i+3,j+3) = cmat(i,j)*cb
        end do
      end do
c
c.... shear
      dmat(7,7) = 0.5d0*d(40)/(1.0d0+d(41))*cs  
      dmat(8,8) = dmat(7,7)

      return
      end
c
      subroutine kine63(eps,x0,r0,zs,wcappa,epstr,dets)
      implicit double precision (a-h,o-z)
      dimension  glu(3,3),x0(3,3),r0(3,3),eps(*),epstr(*)
c
c....  layer strains  E_n+1 = A * E  and  det shifter tensor
c
       epstr(1) = eps(1) + zs*eps(4) 
       epstr(2) = eps(2) + zs*eps(5) 
       epstr(3) = eps(3) + zs*eps(6) 
       epstr(4) = eps(7)*wcappa
       epstr(5) = eps(8)*wcappa
 
       do i = 1,3
        glu(i,1) = x0(i,1) + zs*r0(i,1)
        glu(i,2) = x0(i,2) + zs*r0(i,2)
       end do 

       call vecp(glu(1,1),glu(1,2),glu(1,3))  
       dets = dot(glu(1,3),r0(1,3),3)
       call vecp(x0(1,1),x0(1,2),glu(1,3))  
       dets = dabs(dets)/dsqrt(dot(glu(1,3),glu(1,3),3))
c
      return
      end
c
      subroutine dmat63 (dmat,cmat,sig,sigp,wz,zs,wcappa,cappa)
c-----------------------------------------------------------------------
c       sigp  stresses of layer point
c       sig   stress resultants   
c       cmat  consistent tangent matrix 
c       dmat  linearized stress resultants
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension dmat(8,*),cmat(5,*),sig(*),sigp(*)

      zs2 = zs*zs
      
      do 10 i=1,3
      do 10 j=1,i
        dmat(i  ,j  ) = dmat(i  ,j  ) + cmat(i,j)*wz
        dmat(i+3,j+3) = dmat(i+3,j+3) + cmat(i,j)*wz*zs2
10    continue

      do 20 i=1,3
        dmat(6,i) = dmat(6,i) + cmat(3,i)*wz*zs
      do 20 j=1,2
        dmat(j+3 ,i)  = dmat(j+3 ,i) + cmat(j,i)*wz*zs
        dmat(j+6 ,i)  = dmat(j+6 ,i) + cmat(j+3,i)*wz       *wcappa
        dmat(j+6,i+3) = dmat(j+6,i+3) + cmat(j+3,i)*wz*zs   *wcappa
20    continue

      dmat(7,7) = dmat(7,7) + cmat(4,4)*wz*cappa
      dmat(8,7) = dmat(8,7) + cmat(5,4)*wz*cappa
      dmat(8,8) = dmat(8,8) + cmat(5,5)*wz*cappa

c.... upper part of dmat 
      do 30 i = 1,8
      do 30 j = i+1,8
30      dmat(i,j) = dmat(j,i)

c...  compute stress resultants
      do 40 i=1,3
        sig(i  ) = sig(i  ) + sigp(i)*wz
        sig(i+3) = sig(i+3) + sigp(i)*wz*zs
40    continue
        sig(7)  = sig(7) + sigp(4)*wz*wcappa
        sig(8)  = sig(8) + sigp(5)*wz*wcappa
      return
      end
c
      subroutine faca63 (an,x,wn,aix)
c-----------------------------------------------------------------------
c     matrix M(x) for geometrical matrix 
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension an(3),x(3),wn(3),b(3),t(3),aix(3,3) 
c
      th2  = dot(wn,wn,3)
      th   = dsqrt(th2)
      call vecp (an,x,b)
      dotbw = dot(b,wn,3)

      if (th.lt.1.d-2)then
       c3   =  1.d0/6.d0  * (1.d0 + th2/60d0)
       c10  =  1.d0/6.d0  * (1.d0 + th2/30d0)
       c11  = -1.d0/360d0 * (1.d0 + th2/21d0)
      else
       sint = dsin(th)
       cost = dcos(th) - 1.d0
       c3   = (th*sint+2.d0*cost)/(th2*cost)
       c10  = (sint - th )/( 2.d0 *th *cost )
       c11  = (th2 + th*sint + 4.d0*cost)/( 2.d0*th2*th2*cost)
      end if

      c10  = c10 * dotbw  - dot(an,x,3)
      c11  = c11 * dotbw  

      do i = 1,3
       t(i) =  c11*wn(i) - c3*b(i) 
       do j = 1,i
        aix(i,j) = 0.5d0*(an(i)*x(j) + x(i)*an(j)) 
     1           + 0.5d0*(t(i)*wn(j) + wn(i)*t(j)) 
        aix(j,i) = aix(i,j)
       end do
       aix(i,i) = aix(i,i) + c10
      end do
      return
      end
c
      subroutine transsig63(t0,xu,rd,rin,shp,sig,siga,nen)
c-----------------------------------------------------------------------
c.... transform the stresses to basis from igeo
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension t0(3,3),glu(3,3),glo(3,3),gko(3,3),rl(3,3),xu(3,3),
     +          rd(3,3),ts(3,3),sig(*),siga(9),sigs(9),t(3,3),
     +          rin(3,3,nen),shp(3,nen) 
csk     +          rin(3,3,4),shp(3,nen) 
c
c.....covariant and contravariant base vectors
c
      do i=1,3
       glu(i,1) = xu(i,1)
       glu(i,2) = xu(i,2)
       glu(i,3) = rd(i,3)
      enddo
      do i=1,3
       do j=1,3
         gko(i,j) = dot(glu(1,i),glu(1,j),3)
       enddo 
      enddo
      call invert(gko,3,3)
      do i=1,3
       glo(i,1)=gko(1,1)*glu(i,1)+gko(1,2)*glu(i,2)+gko(1,3)*glu(i,3)
       glo(i,2)=gko(2,1)*glu(i,1)+gko(2,2)*glu(i,2)+gko(2,3)*glu(i,3)
       glo(i,3)=gko(3,1)*glu(i,1)+gko(3,2)*glu(i,2)+gko(3,3)*glu(i,3)
      end do
c.... Transformation coefficients  lambda_ij=rd(i)*gko(j)
      do i=1,3
       do j=1,3
        rl(i,j) = dot(rd(1,i),glo(1,j),3)
       enddo 
      enddo
c.... obtain stress resultants from effective stress resultants
c.... normal forces 1,2,3,4  bending moments 5,6,7 shear forces 8,9
      siga(1) = sig(1) + sig(4)*rl(1,1) + sig(6)*rl(2,1)
      siga(2) = sig(2) + sig(6)*rl(1,2) + sig(5)*rl(2,2)
      siga(3) = sig(3) + sig(4)*rl(1,2) + sig(6)*rl(2,2)
      siga(4) = sig(3) + sig(6)*rl(1,1) + sig(5)*rl(2,1)
      siga(5) = sig(4)
      siga(6) = sig(5)
      siga(7) = sig(6)
      siga(8) = sig(7) + sig(4)*rl(1,3) + sig(6)*rl(2,3)
      siga(9) = sig(8) + sig(6)*rl(1,3) + sig(5)*rl(2,3)
c
c.... local cartesian basis
c
      do i = 1,3
       do j = 1,3
        ts(i,j) = 0.d0 
        do k = 1,nen
         ts(i,j) = ts(i,j) + rin(i,j,k)*shp(3,k)
        enddo 
       enddo
      enddo
      call plamina_base(ts)
c
c.... Transformation matrix  T_ij=ts(i)*t0(j)
      do i = 1,3
       do j=1,3
        t(i,j) = dot(ts(1,i),t0(1,j),3) 
       end do
      end do
c        
c.... normal forces 
      sigs(1) =      t(1,1)*t(1,1)*siga(1) + t(1,2)*t(1,2)*siga(2)
     +        +      t(1,1)*t(1,2)*siga(3) + t(1,1)*t(1,2)*siga(4)
      sigs(2) =      t(2,1)*t(2,1)*siga(1) + t(2,2)*t(2,2)*siga(2)
     +        +      t(2,1)*t(2,2)*siga(3) + t(2,1)*t(2,2)*siga(4)
      sigs(3) =      t(1,1)*t(2,1)*siga(1) + t(1,2)*t(2,2)*siga(2)
     +        +      t(1,1)*t(2,2)*siga(3) + t(1,2)*t(2,1)*siga(4)
      sigs(4) =      t(1,1)*t(2,1)*siga(1) + t(1,2)*t(2,2)*siga(2)
     +        +      t(1,1)*t(2,2)*siga(4) + t(1,2)*t(2,1)*siga(3)
c.... bending moments 
      sigs(5) =      t(1,1)*t(1,1)*siga(5) + t(1,2)*t(1,2)*siga(6)
     +        + 2.d0*t(1,1)*t(1,2)*siga(7)
      sigs(6) =      t(2,1)*t(2,1)*siga(5) + t(2,2)*t(2,2)*siga(6)
     +        + 2.d0*t(2,1)*t(2,2)*siga(7)
      sigs(7) =      t(1,1)*t(2,1)*siga(5) + t(1,2)*t(2,2)*siga(6)
     +        +      t(1,1)*t(2,2)*siga(7) + t(1,2)*t(2,1)*siga(7)
c.... shear forces
      sigs(8) =      t(1,1)*siga(8) + t(1,2)*siga(9)
      sigs(9) =      t(2,1)*siga(8) + t(2,2)*siga(9)
c.... store
      call matcop (sigs,9,1, siga)          
      return
      end
c
      subroutine plsn63 (i)
c-----------------------------------------------------------------------
c.... define name for stresses in plot-output 
c     i = 0: stress resultants
c         1: stresses and history values at defined layer
c-----------------------------------------------------------------------
      USE strnam
      implicit  double precision (a-h,o-z)

      istv = -9
      do is =1,25
        strsus(is) = '               '
      end do

      if (i.eq.0) then
        strsus(1) = 'Resultant N_11 '
        strsus(2) = 'Resultant N_22 '
        strsus(3) = 'Resultant N_12 '
        strsus(4) = 'Resultant N_21 '
        strsus(5) = 'Resultant M_11 '
        strsus(6) = 'Resultant M_22 '
        strsus(7) = 'Resultant M_12 '
        strsus(8) = 'Resultant Q_13 '
        strsus(9) = 'Resultant Q_23 '
      else           
        strsus(1) = '  Sigma_11     '
        strsus(2) = '  Sigma_22     '
        strsus(3) = '  Sigma_12     '
        strsus(4) = '  Sigma_13     '
        strsus(5) = '  Sigma_23     '
        strsus(6) = 'v. Mises stress'
        strsus(7) = 'eq plast strain'
        strsus(8) = '  Eps_33     '
      end if
      return
      end
c
      subroutine zeroElementStressResults63 (NURstress,numel,n) 
c-----------------------------------------------------------------------   
c     set Stresses to Zero for current element
      implicit  double precision (a-h,o-z)
      real*8 NURstress(9*numel,11)
      do l = 1,9
        NURstress((n-1)*9+l,1:11)=0.0d0
      end do
      return
      end
c
      subroutine plot63 (d,ix,dt,st,shp,sig,sigp,h2,eps,
     +                   numnp,maxlay,klay,da,nen,ixold)
c-----------------------------------------------------------------------
c.... Plot stresses  at specified layer for general shell element
c     klay = 0:           stress resultants
c          = 1 to maxlay: stresses and history values at layer klay
c-----------------------------------------------------------------------
      implicit  double precision (a-h,o-z)
      dimension ix(*),dt(numnp),st(numnp,*),shp(3,*),sig(*),sigp(*),
     + h2(*),d(*),eps(*),ixold(*)
      if (klay.gt.0 .and. klay.le.maxlay) then
c....   print stresses and history values in layer klay
        do i = 1,5
          sig(i) = sigp(i)
        end do
        sig(6) = sigp(6)
        sig(7) = h2(7)
        sig(8) = h2(8)
      end if

c.... membrane energy   
c      hs     = d(13)
c      emod   = d(40)  
c      xnu    = d(41)
c      Dm     = emod*hs/(1.d0-xnu*xnu)
c      wm = 0.5d0*Dm*(eps(1)*eps(1)+ eps(2)*eps(2)+2.d0*xnu*eps(1)*eps(2)
c     1     +(1.d0-xnu)*eps(3)*eps(3)/2.d0) 
      
c...  store values
      do 20 j = 1,nen
        xsji = shp(3,j)*da
        ii = abs(ixold(j))
        !if (ix(j).ne.ixold(j)) then
        !   write(*,*) ixold(j)
        !end if
        if(ii.eq.0) go to 20
        dt(ii) = dt(ii) + xsji
        do 30 i = 1,9
30        st(ii,i)  = st(ii,i) + sig(i)*xsji
c          st(ii,10) = st(ii,10) + wm*xsji
20    continue
      return
      end
c
      subroutine pout63 (xl,shp,sig,sigp,epn,ndm,maxlay)
      USE iofile
      USE bdata
      USE plslay
      USE eldata
c-----------------------------------------------------------------------
c.... print stresses at Gauss points of each layer
c     klay = 0:           stress resultants
c          = 1 to maxlay: stresses and history values at layer klay
c
c     1. PK-Stresses/stress resultants with respect to deformed ref.surf.
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      include 'bdata.h'      
c      include 'plslay.h'
c      include 'iofile.h'
c      include 'eldata.h'
      dimension xl(ndm,*),shp(3,*),sig(*),sigp(*),gp(3),epn(*)

c.... coordinates of Gauss points
      call pzero(gp,3)
      do 10 idm = 1,3
       do 10 inod = 1,nel
10    gp(idm) = gp(idm) + xl(idm,inod) * shp(3,inod)

      if (klay .eq. 0) then
c....   print stress resultants
        mct = mct - 2
        if(mct.le.0) then
          mct = 50
          if(ior.lt.0) write (  *,1000) o,head
                       write (iow,1000) o,head
1000      format(a1,20a4,//,2x,'E L E M E N T  63 ',
     1       'Stress resultants',/, 
     2    2x,'Elmt',7x,'1-Coord',5x,'N_11',9x,'N_22',9x,'N_12',/,
     3    2x,' Mat',7x,'2-Coord',5x,'Q_11',9x,'Q_22',9x,'N_21',/,
     4    2x,'    ',7x,'3-Coord',5x,'M_11',9x,'M_22',9x,'M_12')
        end if 
                     write(iow,1010)  n,gp(1),(sig(i),i=1,3),
     1                               ma,gp(2), sig(8),sig(9),sig(4),
     2                                  gp(3),(sig(i),i=5,7)
        if(ior.lt.0) write(  *,1010)  n,gp(1),(sig(i),i=1,3),
     1                               ma,gp(2), sig(8),sig(9),sig(4),
     2                                  gp(3),(sig(i),i=5,7)

1010    format(1x,i5,1x,f13.4,3(1x,e12.5),/,
     1         1x,i5,1x,f13.4,3(1x,e12.5),/,
     2               7x,f13.4,3(1x,e12.5))

cww1000      format(a1,20a4,//,2x,'E L E M E N T  Stress Resultants',/,
cww     1    2x,'El',1x,'Mat',1x,'1-Coord',1x,'2-Coord',1x,'3-Coord',
cww     2    2X,'   N_11 ',3X,'   N_22 ',3X,'   N_12 ',3X,'   N_21 ',/,
cww     3   34X,'   M_11 ',3X,'   M_22 ',3X,'   M_12 ',3X,'        ',/,
cww     4   34X,'   Q_11 ',3X,'   Q_22  ',3X,'        ',3X,'        ',/)
cww        end if 
cww        write (iow,1010) n,ma,(gp(i),i=1,3),
cww     1                   (sig(i),i=1,4),(sig(i),i=5,7),sig(8),sig(9)
cww        if(ior.lt.0) write(*,1010) n,ma,(gp(i),i=1,3),
cww     1                   (sig(i),i=1,4),(sig(i),i=5,7),sig(8),sig(9)
cww1010    format(1x,i4,i2,3f8.3,1x,4e11.4,/,32x,3e11.4,/,32x,2e11.4)

      elseif (klay.gt.0 .and. klay.le.maxlay) then
c....   print stresses and history values in layer klay
        mct = mct - 2
        if(mct.le.0) then
          mct = 50
          if (ior.lt.0) write (  *,1001) o,head,klay
                        write (iow,1001) o,head,klay
1001      format(a1,20a4,//,2x,'E L E M E N T  63  Cauchy stresses',
     1    ' for layer',i5,/,
     2    2x,'Elmt',7x,'1-Coord',2x,'Sigma_11',5x,'Sigma_22',
     3                           5x,'Sigma_12',5x,'Sigma_13',/,
     4    2x,' Mat',7x,'2-Coord',2x,'Epspl_11',5x,'Epspl_22',
     5                           5x,'Epspl_33',5x,'Sigma_23',/,
     6    2x,'    ',7x,'3-Coord',2x,'Epspl_12',5x,'Epspl_13',
     7                           5x,'Epspl_23',5x,'ef.plast')
        end if    
                     write(iow,1011) n,gp(1),(sigp(i),i=1,4),
     1                              ma,gp(2),(epn(i), i=1,3),sigp(5),        
     2                                       (epn(i),i=4,7)
        if(ior.lt.0) write(  *,1011) n,gp(1),(sigp(i),i=1,4),
     1                              ma,gp(2),(epn(i), i=1,3),sigp(5),        
     2                                       (epn(i),i=4,7)

1011    format(1x,i5,1x,f13.4,4(1x,e12.5),/,
     1         1x,i5,1x,f13.4,4(1x,e12.5),/,
     2               7x,f13.4,4(1x,e12.5))


cww1001      format(a1,20a4,//,2x,'E L E M E N T  63   S T R E S S E S',
cww     *    ' for layer',i5,/,
cww     1    2x,'El',1x,'Mat',1x,'1-Coord',1x,'2-Coord',1x,'3-Coord',
cww     2    2X,'Sigma_1 ',3X,'Sigma_2 ',3X,'Sigma_12',3X,'Sigma_13',/,
cww     3   34X,'Epspl_1 ',3X,'Epspl_2 ',3X,'Epspl_3 ',3X,'Sigma_23',/,
cww     4   34X,'Epspl_12',3X,'Epspl_13',3X,'Epspl_23',3X,'ef.plast',/)
cww        end if    
cww        write (iow,1010) n,ma,(gp(i),i=1,3),
cww     1                  (sigp(i),i=1,4),(epn(i),i=1,3),sigp(5),
cww     2                  (epn(i),i=4,7)
cww        if(ior.lt.0) write(*,1010) n,ma,(gp(i),i=1,3),
cww     1                    (sigp(i),i=1,4),(epn(i),i=1,3),sigp(5),
cww     2                    (epn(i),i=4,7)
cww

      end if
      return
      end
c
      subroutine stif63 (inode,jnode,ir0,jc0,b,btd,
     1                 st,sig,qs,shp,shp1,w,ngeas_max,ndf,da,ilin,isw,
     2                 noans,iuseTgp,xu,dwg,rd,twn,igauss,HutT_1,HutT_2,
     3                 twn_1,twn_2,iMsize,iRupd,Mhq,Mh1,Mh2)
c-----------------------------------------------------------------------
c.... stiffness matrix
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  st(ngeas_max,*),sig(*),shp(3,*),shp1(3,*),w(3,3,*),
     1           id1(4,4),id2(4,4),qs(*),b(8,6,*),btd(8,*),xmh(3),
     2           xmh_1(3),xmh_2(3),wmh(3,3),wh(3,3),
     3           xu(3,*),dwg(9,*),rd(3,3),twn(3,3,*),
     4           Hutmuw(3,3),Hutmwu(3,3),HutT_1(3,3,*),HutT_2(3,3,*),
     5           twn_1(3,3,*),twn_2(3,3,*),xmh1(3),xmh2(3),Gww(3,3),
     6           iMsize(*)  
      double precision Mhq(3,3),Mh1(3,3),Mh2(3,3)
      data id1/ 1,1,0,0, 1,1,0,0, 0,0,1,1,  0,0,1,1/
      data id2/ 1,0,0,1, 0,1,1,0, 0,1,1,0,  1,0,0,1/
c
c.... material part     
      if(isw.eq.3)then
        do i  = 1,ndf
          ir   = ir0+i
          do j = 1,ndf
            jc  = jc0+j
            st(ir,jc) = st(ir,jc) + dot(btd(1,i),b(1,j,jnode),8)
          end do
        end do
      endif
      if(ilin.eq.0)return
c      
c.... geometrical part
      if (noans.eq.0) then
        if1  = id1(inode,jnode)
        if2  = id2(inode,jnode)
      end if
c.... sgd1: ^n_IK 
      sgd1 =  sig(1)*shp(1,inode)*shp(1,jnode) 
     1        + sig(2)*shp(2,inode)*shp(2,jnode)
     2        + sig(3)*shp(1,inode)*shp(2,jnode) 
     3        + sig(3)*shp(2,inode)*shp(1,jnode)
c.... compute G_uß and G_ßu part of G_IK
      if (iuseTgp.lt.2) then
c....   sgd2:  ^m_IK
        sgd2  = sig(4)*shp(1,inode)*shp(1,jnode) 
     1        + sig(5)*shp(2,inode)*shp(2,jnode)
     2        + sig(6)*shp(1,inode)*shp(2,jnode) 
     3        + sig(6)*shp(2,inode)*shp(1,jnode)
        if (noans.eq.0) then
          sgd3 = 0.5d0*(qs(1)*shp1(1,inode)*if1+qs(2)*shp1(2,inode)*if2)
          sgd4 = 0.5d0*(qs(1)*shp1(1,jnode)*if1+qs(2)*shp1(2,jnode)*if2)
        else
c...      sgd3: ^q^uß_IK, sgd4: ^q^ßu_IK
          sgd3 = sig(7)*shp(1,inode)*shp(3,jnode)
     1         + sig(8)*shp(2,inode)*shp(3,jnode)
          sgd4 = sig(7)*shp(1,jnode)*shp(3,inode)
     1         + sig(8)*shp(2,jnode)*shp(3,inode)
        end if
        guw = (sgd2+sgd3) * da
        gwu = (sgd2+sgd4) * da
      else
        Hutmuw=0.d0
        Hutmwu=0.d0
        do i = 1,3
          do j = 1,iMsize(jnode)
c           compute ^m_IK^uß + ^q_IK^1uß = G_IK^uß
            Hutmuw(i,j)=((shp(1,inode)*sig(4)+shp(2,inode)*sig(6))*
     *                   HutT_1(i,j,jnode) +
     +                   (shp(2,inode)*sig(5)+shp(1,inode)*sig(6))*
     *                   HutT_2(i,j,jnode) +
     +                   (shp(1,inode)*sig(7)+shp(2,inode)*sig(8))*
     *                   shp(3,jnode)*w(i,j,jnode))*da
          end do
       end do  
       do i = 1,3
          do j = 1,iMsize(inode)
c           compute ^m_IK^uß + ^q_IK^1uß = G_IK^uß
            Hutmwu(i,j)=((shp(1,jnode)*sig(4)+shp(2,jnode)*sig(6))*
     *                   HutT_1(i,j,inode) +
     +                   (shp(2,jnode)*sig(5)+shp(1,jnode)*sig(6))*
     *                   HutT_2(i,j,inode) +
     +                   (shp(1,jnode)*sig(7)+shp(2,jnode)*sig(8))*
     *                   shp(3,inode)*w(i,j,inode))*da        
          end do
        end do  
      end if
c
      guu =  sgd1       * da


c...  compute G_ßß part of G_IK for Gauss-Point Rodrigues
c     only works for disabled ANS
      if (iuseTgp.ge.2.and.ilin.ge.2) then
        if (iRupd.eq.2) then
          !do nothing
        else
          do i = 1,3
c...      compute vectors h(only for q), h1, h2, h_,1 and h_,2
!          xmh(i) = ( sig(4)*xu(i,1)*( shp(1,inode)*shp(3,jnode) 
!     +                               +shp(3,inode)*shp(1,jnode))
!     +              +sig(5)*xu(i,2)*( shp(2,inode)*shp(3,jnode)
!     +                               +shp(3,inode)*shp(2,jnode))
!     +              +sig(6)*xu(i,2)*( shp(1,inode)*shp(3,jnode)
!     +                               +shp(3,inode)*shp(1,jnode))
!     +              +sig(6)*xu(i,1)*( shp(2,inode)*shp(3,jnode)
!     +                               +shp(3,inode)*shp(2,jnode))
!     +              +sig(7)*xu(i,1)*shp(3,inode)*shp(3,jnode)
!     +              +sig(8)*xu(i,2)*shp(3,inode)*shp(3,jnode) )*da
          xmh1(i) = (sig(4)*xu(i,1) + sig(6)*xu(i,2))*da
          xmh2(i) = (sig(5)*xu(i,2) + sig(6)*xu(i,1))*da
          xmh(i)  = (sig(7)*xu(i,1)+sig(8)*xu(i,2))
     *              *shp(3,inode)*shp(3,jnode) *da
          xmh_1(i)= (sig(4)*xu(i,1)+sig(6)*xu(i,2))*
     *              shp(3,inode)*shp(3,jnode)*da
          xmh_2(i)= (sig(5)*xu(i,2)+sig(6)*xu(i,1))*
     *              shp(3,inode)*shp(3,jnode)*da
          end do
        end if
c...    compute G_IK^ww part of geometrical matrix
        Ilen = iMsize(inode)
        Jlen = iMsize(jnode)
        if (iRupd.eq.1)
     1    call computeGww63(Gww,dwg(1:3,igauss),dwg(4:6,igauss),
     1                    dwg(7:9,igauss),rd(1,3),rd(1,1),rd(1,2),
     2                    xmh,xmh_1,xmh_2,xmh1,xmh2,twn,twn_1,twn_2,
     3                    inode,jnode,shp,Ilen,Jlen)
c        Fast version. M matrices are computed in Gauss loop for speed
c        This is the fastest version. 2,3 works and is checked, 2,2 seems to work also        
         if (iRupd.eq.2)!.and.iuseTgp.eq.3)
     1    call computeGww63_Om_Mult_fast(Gww,twn,inode,jnode,shp,
     2                 Ilen,Jlen,Mhq,Mh1,Mh2)  
c        Quite fast version. M matrices are summed up before computing T_3I^T*M*T_3K
c        2,3 works and is checked, 2,2 seems to work also               
!         if (iRupd.eq.2.and.iuseTgp.eq.2)
!     1    call computeGww63_Om_Mult_test(Gww,rd(1,3),rd(1,1),
!     2                    rd(1,2),xmh,xmh_1,xmh_2,xmh1,xmh2,twn,twn_1,
!     3                    twn_2,inode,jnode,shp,Ilen,Jlen)
c        initial version with simplified M matrices for multiplicative update, works for
c        2,3 and 2,2 and is checked, but very slow         
!        call computeGww63_Om_Mult(Gww,rd(1,3),rd(1,1),rd(1,2),
!     2                    xmh,xmh_1,xmh_2,xmh1,xmh2,twn,twn_1,twn_2,
!     3                    inode,jnode,shp,Ilen,Jlen)
c...    store sum(M) to right place in stiffness matrix        
        do i = 1,ndf-3
          ir2 = ir0+3+i
          do j = 1,ndf-3
            ic2 = jc0+3+j
          st(ir2,ic2) = st(ir2,ic2) + Gww(i,j)
          end do
        end do  
      end if
      do i = 1,3
       ir1 =   ir0+i
       ic1 =   jc0+i
       st(ir1,ic1) = st(ir1,ic1) + guu
       do j = 1,ndf-3
         ir2 = ir0+3+j
         ic2 = jc0+3+j
         if (iuseTgp.lt.2) then
           st(ir1,ic2) = st(ir1,ic2) + w(i,j,jnode) * guw
           st(ir2,ic1) = st(ir2,ic1) + w(i,j,inode) * gwu
         else
           st(ir1,ic2) = st(ir1,ic2) + Hutmuw(i,j)
           st(ir2,ic1) = st(ir2,ic1) + Hutmwu(i,j)
         end if
       end do
      end do
c
      return
      end
c

      subroutine stra63 (xl,ul,shp,xsi,eta,eas,eps,sx,gam,rin,rwn,
     1                   ddn,xu,x0,rd,r0,ndm,ndf,ilin,noans,nen,T,
     2                   iuseTgp,dd)
c-----------------------------------------------------------------------
c     Green-Lagrangean strains of the reference surface
c     membrane(1,2,3),     curvatures(4,5,6),    shear(7,8) 
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension xl(ndm,*),ul(ndf,*),shp(3,*),ddn(3,nen),
     1          eps(*),rin(3,3,*),rwn(3,3,*),uu(3,3),dd(3,3),eas(*),
     2          x0(3,3),xu(3,3),r0(3,3),rd(3,3),gam(*),sx(2,2),T(3,2)
c
c.... current and initial geometry
c
      do i = 1,3
       do j = 1,3
        if (iuseTgp.lt.2) then
          rd(i,j)  = 0.d0
          r0(i,j)  = 0.d0
        end if
        xu(i,j)  = 0.d0
        x0(i,j)  = 0.d0
         do k = 1,nen
          if (iuseTgp.lt.2) then
            rd(i,j)  = rd(i,j)  + rwn(i,3,k) * shp(j,k)
            r0(i,j)  = r0(i,j)  + rin(i,3,k) * shp(j,k)
          end if
          xu(i,j)  = xu(i,j)  + (xl(i,k)+ul(i,k))*shp(j,k)
          x0(i,j)  = x0(i,j)  +  xl(i,k)         *shp(j,k)
         end do
        end do 
c       x0(i,3) = r0(i,3) !makes x0 the covariant base
      end do

      gam1 = 0.5d0*((1.d0-eta)*gam(2)+(1.d0+eta)*gam(4))
      gam2 = 0.5d0*((1.d0-xsi)*gam(1)+(1.d0+xsi)*gam(3))
   
c
      if (ilin.ne.0)then
       eps(1) = 0.5d0*(dot(xu(1,1),xu(1,1),3) - dot(x0(1,1),x0(1,1),3)) 
       eps(2) = 0.5d0*(dot(xu(1,2),xu(1,2),3) - dot(x0(1,2),x0(1,2),3)) 
       eps(3) =        dot(xu(1,1),xu(1,2),3) - dot(x0(1,1),x0(1,2),3)  
       eps(4) = dot(xu(1,1),rd(1,1),3) - dot(x0(1,1),r0(1,1),3) 
       eps(5) = dot(xu(1,2),rd(1,2),3) - dot(x0(1,2),r0(1,2),3) 
       eps(6) = dot(xu(1,1),rd(1,2),3) + dot(xu(1,2),rd(1,1),3)
     +         -dot(x0(1,1),r0(1,2),3) - dot(x0(1,2),r0(1,1),3) 
c
      else 
 
       do i = 1,3
        do j = 1,3
         uu(i,j) = 0.d0
         if (iuseTgp.lt.2) dd(i,j) = 0.d0
         xu(i,j) = x0(i,j)  
         rd(i,j) = r0(i,j)  
         do k  = 1,nen
           uu(i,j) = uu(i,j) + ul(i,k)*shp(j,k)
           if (iuseTgp.lt.2) then
             dd(i,j) = dd(i,j) + ddn(i,k)*shp(j,k)
           end if
          end do
        end do
       end do 
c
       eps(1) = dot(x0(1,1),uu(1,1),3) 
       eps(2) = dot(x0(1,2),uu(1,2),3) 
       eps(3) = dot(x0(1,1),uu(1,2),3) + dot(x0(1,2),uu(1,1),3) 
       eps(4) = dot(x0(1,1),dd(1,1),3) + dot(r0(1,1),uu(1,1),3)
       eps(5) = dot(x0(1,2),dd(1,2),3) + dot(r0(1,2),uu(1,2),3)
       eps(6) = dot(x0(1,1),dd(1,2),3) + dot(x0(1,2),dd(1,1),3)
     1        + dot(r0(1,1),uu(1,2),3) + dot(r0(1,2),uu(1,1),3)

      end if
c
c....  add enhanced strains
      eps(1) = eps(1) + eas(1) 
      eps(2) = eps(2) + eas(2) 
      eps(3) = eps(3) + eas(3) 
      eps(4) = eps(4) + eas(4) 
      eps(5) = eps(5) + eas(5) 
      eps(6) = eps(6) + eas(6) 
c....  transverse shear strains
cwd   for ANS enabled/disabled
      if (noans.eq.0) then ! ANS enabled
        eps(7) = sx(1,1)*gam1 + sx(1,2)*gam2 + eas(7)
        eps(8) = sx(2,1)*gam1 + sx(2,2)*gam2 + eas(8)
      elseif (noans.eq.1.and.ilin.ne.0) then ! ANS diabled
        eps(7) = dot(xu(1,1),rd(1,3),3) - dot(x0(1,1),r0(1,3),3) +eas(7)
        eps(8) = dot(xu(1,2),rd(1,3),3) - dot(x0(1,2),r0(1,3),3) +eas(8)
      elseif (noans.eq.1.and.ilin.eq.0) then ! ANS diabled
        eps(7) = dot(r0(1,3),uu(1,1),3)+dot(x0(1,1),dd(1,3),3) + eas(7) 
        eps(8) = dot(r0(1,3),uu(1,2),3)+dot(x0(1,2),dd(1,3),3) + eas(8)
      end if
c   
      return
      end
c
      subroutine strt63 (d,q,numel,n,mqloa,propq,eps,epse,isw)
      USE prlod
c-----------------------------------------------------------------------
c     - temperature strains
c     membrane(1,2,3),     curvatures(4,5,6),    shear(7,8) 
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      include 'prlod.h'
      dimension d(*),eps(*),epse(*),q(numel,10)
c
      alphat = d(29) 
      hs     = d(13)
      hm     = d(14)
      hp     = hm+hs   
      atn    = 0.d0
      atm    = 0.d0  

      if(isw.ne.22.and.mqloa.eq.1)then
        atn = prop * alphat * (d(30)*hp - d(31)*hm)/hs
        atm = prop * alphat * (d(31)    - d(30))   /hs
      endif
      if((isw.eq.22.or.isw.eq.4.or.isw.eq.8).and.mqloa.ne.1)then
        atn = propq * alphat * (q(n,6)*hp - q(n,7)*hm)/hs
        atm = propq * alphat * (q(n,7)    - q(n,6))   /hs
      endif

      if(mqloa.eq.1)then
        atne = prop * alphat * (d(30)*hp - d(31)*hm)/hs
        atme = prop * alphat * (d(31)    - d(30))   /hs
      else
        atne = propq * alphat * (q(n,6)*hp - q(n,7)*hm)/hs
        atme = propq * alphat * (q(n,7)    - q(n,6))   /hs
      endif 

c.... - temperature strains
      epse(1) = eps(1) - atne
      epse(2) = eps(2) - atne
      epse(3) = eps(3) 
      epse(4) = eps(4) - atme
      epse(5) = eps(5) - atme
      epse(6) = eps(6) 
      epse(7) = eps(7) 
      epse(8) = eps(8) 

      eps(1) = eps(1) - atn
      eps(2) = eps(2) - atn
      eps(4) = eps(4) - atm
      eps(5) = eps(5) - atm
c
      return
      end
c
      subroutine trasig63 (sigp,t0,shp,rin,nen)               
c-----------------------------------------------------------------------
c.... transform stresses         
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  sigp(*),t0(3,3),ts(3,3),t(3,3),
     1           tr(5,5),sigt(6),rin(3,3,3),shp(3,nen)
c
      do i = 1,3
       do j = 1,3
        ts(i,j) = 0.d0 
        do k = 1,nen
         ts(i,j) = ts(i,j) + rin(i,j,k)*shp(3,k)
        enddo 
       enddo
      enddo
      call plamina_base(ts)
c
      do i = 1,3
       do j = 1,3
        t(i,j) = dot(ts(1,i), t0(1,j), 3)
       enddo
      enddo 
c
      call pzero(tr,5*5)
      
      tr(1,1) = t(1,1)*t(1,1)
      tr(2,1) = t(2,1)*t(2,1)
      tr(3,1) = t(1,1)*t(2,1)

      tr(1,2) = t(1,2)*t(1,2)
      tr(2,2) = t(2,2)*t(2,2)
      tr(3,2) = t(1,2)*t(2,2)

      tr(1,3) = t(1,1)*t(1,2)*2.d0
      tr(2,3) = t(2,1)*t(2,2)*2.d0
      tr(3,3) = t(2,1)*t(1,2) + t(1,1)*t(2,2)

      tr(4,4) = t(1,1)
      tr(5,4) = t(2,1)

      tr(4,5) = t(1,2)
      tr(5,5) = t(2,2)

      call matmulf (tr,sigp,5,5,1, sigt)
      call matcop (sigt,5,1, sigp)          

      return
      end
c
      subroutine triad63 (igeo,xl,d,rin,ndm,nen,t0)
c-----------------------------------------------------------------------
c     initial cartesian system for specified geometry:
c     igeo =  1      ... arbitrary surface 
c             others ... see macro BASE
c
c     rin(3,3,4)     initial local triad at nodes
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  xl(ndm,*),g1(3,4),g2(3,4),rin(3,3,*),e3(3),d(*),add(3)
      dimension  t0(3,3)
      data e3 / 0.d0,0.d0,1.d0 /

          add(1)=d(2)
          add(2)=d(3)
          add(3)=d(4)

c
      call pzero(rin,9 * nen)
c
!      do 20 i = 1,3
!       g1(i,1) = xl(i,2)-xl(i,1)
!       g1(i,2) = g1(i,1)
!       g1(i,3) = xl(i,3)-xl(i,4)
!       g1(i,4) = g1(i,3)
!       g2(i,1) = xl(i,4)-xl(i,1)
!       g2(i,2) = xl(i,3)-xl(i,2)
!       g2(i,3) = g2(i,2)
!       g2(i,4) = g2(i,1)
!20    continue

      do k = 1,nen

       if(igeo.eq.1.or.igeo.eq.2) then
        rin(1:3,1,k)=t0(1:3,1)
        rin(1:3,2,k)=t0(1:3,2)
        rin(1:3,3,k)=t0(1:3,3)  
!        call vecp (g1(1,k),g2(1,k),rin(1,3,k))
!        call norm (rin(1,3,k),rin(1,3,k),3)
!        call norm (rin(1,1,k),g1(1,k),3)
!        call vecp (rin(1,3,k),rin(1,1,k),rin(1,2,k))
c        call norm (rin(1,2,k),g2(1,k),3)
c        call vecp (rin(1,2,k),rin(1,3,k),rin(1,1,k))
c
        if(igeo.eq.2) then  
          call norm(add,add,3)
c....     e_1
          rin(1,1,k) = add(1)  
          rin(2,1,k) = add(2) 
          rin(3,1,k) = add(3) 
c....     e_2
          call vecp (rin(1,3,k),rin(1,1,k),rin(1,2,k))
        end if
c
c....  others (analytical)
       else
        call pdirec3(rin(1,1,k),add,xl(1,k),igeo)
       end if
      
      end do
      
      return
      end
c
      subroutine updn63 (xl,ul,dwn,twn,rwn,ddn,rin,ran,xan,
     1                   gam,bml,w,ix,ndf,ndm,nen,ityp,ilin,iRupd)
      USE mdata
c-----------------------------------------------------------------------
c     update of nodal rotations
c     rin    initial local cartesian system   at element nodes
c     rwn    current nodal basis  R = R(om)       - " -
c     ran    current       director           at midside nodes
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      include 'mdata.h'      
csk      common m(1)
      dimension  xl(ndm,*),ul(ndf,*),ix(*),
     1           dwn(3,*),twn(3,3,*),rwn(3,3,*),ddn(3,*),rin(3,3,*),
     2           w(3,3,*),rh(3,3),rh1(3,3),rh2(3),gam(*),bml(3,2,*),
     3           ran(3,*),xan(3,*),ran0(3,4),xan0(3,4),
     4           uan(3,4),dan(3,4),im(4),il(4)
      data im /-2, 2, 4,-4/, il /-1,-3, 3, 1/
c
      do 100  k = 1,nen
c.... initial nodal cartesian system
       call matcop (rin(1,1,k),3,3, rh1)
       call matcop (ul(4,2*nen+k),3,1, rh2)          

c.... transformation of rotations 5dofs --> 6dofs
       if(ndf.eq.5) iflag = 1                          ! iflag = 1 : glatte Schale
csk       if(ndf.eq.6) call boun63(k,ix,m(n7),iflag,ndf)  ! iflag = 0 : Kante  
       if(ndf.eq.6) call boun63(k,ix,psid,iflag,ndf)  ! iflag = 0 : Kante  
       if(ityp.eq.6) iflag = 0     
c        iflag= 1 !test for iga_shell
       if (iflag.eq.1) then
         if(ilin.eq.2) then
           call updr63 (dwn(1,k),rh,ilin,1)
           call matmulf (rh,rin(1,1,k),3,3,3, rh1)
         end if
         call matmulf (rh1,ul(4,2*nen+k),3,2,1, rh2)
       end if
c.... update current axial vector Omega    
       call matadd (rh2,dwn(1,k),3,1, dwn(1,k))          

c.... update of current nodal basis
       if(ilin.eq.0)then                   ! linear

         call matcop (rin(1,1,k),3,3, rwn(1,1,k))
         call vecp (dwn(1,k),rin(1,3,k),ddn(1,k))

       else                                ! nonlinear

         call updr63 (dwn(1,k),rh,ilin,1)
         call matmulf (rh,rin(1,1,k),3,3,3, rwn(1,1,k))

       end if 

c....  T = W^T * H * T_3
       call updr63 (dwn(1,k),twn(1,1,k),ilin,3)
       if(ilin.eq.2)then                 !  finite rotations
         call skew (w(1,1,k),rwn(1,3,k))                      
         if (iflag.eq.1) then
           rh=0.d0
           call matmulf (twn(1,1,k),rwn(1,1,k),3,3,2, rh)          
c          call matmulf (twn(1,1,k),rwn(1,1,k),3,3,3, rh) !HIER NUR 2 SPALTEN???
           call matcop (rh,3,3, twn(1,1,k))
         end if

       else                              ! linear and moderate rotations
         call skew (w(1,1,k),rin(1,3,k))                      
         if (iflag.eq.1) call matcop (rin(1,1,k),3,3, twn(1,1,k))
       endif

       call mttmul (w(1,1,k),twn(1,1,k),3,3,3, rh)
       call matcop (rh,3,3, w(1,1,k))

         
100   continue
cwd   from here on only for ANS->not needed at the moment
c.... current director at midside nodes
      do 200 i = 1,3
        ran(i,1) = (rwn(i,3,4) + rwn(i,3,1))/2.d0
        ran(i,2) = (rwn(i,3,1) + rwn(i,3,2))/2.d0
        ran(i,3) = (rwn(i,3,2) + rwn(i,3,3))/2.d0
        ran(i,4) = (rwn(i,3,3) + rwn(i,3,4))/2.d0

        xan(i,1) = ( (xl(i,4)+ul(i,4))-(xl(i,1)+ul(i,1)) )/2.d0
        xan(i,2) = ( (xl(i,2)+ul(i,2))-(xl(i,1)+ul(i,1)) )/2.d0
        xan(i,3) = ( (xl(i,3)+ul(i,3))-(xl(i,2)+ul(i,2)) )/2.d0
        xan(i,4) = ( (xl(i,3)+ul(i,3))-(xl(i,4)+ul(i,4)) )/2.d0

        ran0(i,1) = (rin(i,3,4) + rin(i,3,1))/2.d0
        ran0(i,2) = (rin(i,3,1) + rin(i,3,2))/2.d0
        ran0(i,3) = (rin(i,3,2) + rin(i,3,3))/2.d0
        ran0(i,4) = (rin(i,3,3) + rin(i,3,4))/2.d0

        xan0(i,1) = ( xl(i,4) - xl(i,1) )/2.d0
        xan0(i,2) = ( xl(i,2) - xl(i,1) )/2.d0
        xan0(i,3) = ( xl(i,3) - xl(i,2) )/2.d0
        xan0(i,4) = ( xl(i,3) - xl(i,4) )/2.d0

       if(ilin.eq.0)then 

        uan(i,1) = (ul(i,4) - ul(i,1))/2.d0
        uan(i,2) = (ul(i,2) - ul(i,1))/2.d0
        uan(i,3) = (ul(i,3) - ul(i,2))/2.d0
        uan(i,4) = (ul(i,3) - ul(i,4))/2.d0

        dan(i,1) = (ddn(i,4) + ddn(i,1))/2.d0
        dan(i,2) = (ddn(i,1) + ddn(i,2))/2.d0
        dan(i,3) = (ddn(i,2) + ddn(i,3))/2.d0
        dan(i,4) = (ddn(i,3) + ddn(i,4))/2.d0
       end if

200   continue 

c
c.... shear strains at assumed strain points and factors b1m, b2l, 
      do 220 k = 1,4
        mm = iabs(im(k))
        ll = iabs(il(k))
        if(ilin.eq.0)then
         gam(k) = dot(uan(1,k),ran0(1,k),3) + dot(xan0(1,k),dan(1,k),3)
         call mttmul (w(1,1,k),xan0(1,mm),3,3,1, bml(1,1,k))
         call mttmul (w(1,1,k),xan0(1,ll),3,3,1, bml(1,2,k))
        else
         gam(k) = dot(xan(1,k),ran(1,k),3) - dot(xan0(1,k),ran0(1,k),3)
         call mttmul (w(1,1,k),xan(1,mm),3,3,1, bml(1,1,k))
         call mttmul (w(1,1,k),xan(1,ll),3,3,1, bml(1,2,k))
        end if
220   continue 
c
      return
      end
      subroutine updn63Rupd2 (xl,ul,Rnold,twn,rwn,ddn,rin,ran,xan,
     1                   gam,bml,w,ix,ndf,ndm,nen,ityp,ilin,dwn)
      USE mdata
c-----------------------------------------------------------------------
c     update of nodal rotations
c     rin    initial local cartesian system   at element nodes
c     rwn    current nodal basis  R = R(om)       - " -
c     ran    current       director           at midside nodes
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      include 'mdata.h'      
csk      common m(1)
      dimension  xl(ndm,*),ul(ndf,*),ix(*),
     1           dwn(3,*),twn(3,3,*),rwn(3,3,*),ddn(3,*),rin(3,3,*),
     2           w(3,3,*),rh(3,3),rh1(3,3),rh2(3),gam(*),bml(3,2,*),
     3           ran(3,*),xan(3,*),ran0(3,4),xan0(3,4),
     4           uan(3,4),dan(3,4),im(4),il(4),Rnold(1:3,1:3,nen)
      data im /-2, 2, 4,-4/, il /-1,-3, 3, 1/
c
      do 100  k = 1,nen
c.... initial nodal cartesian system
       call matcop (rin(1,1,k),3,3, rh1)
       call matcop (ul(4,2*nen+k),3,1, rh2)          

c.... transformation of rotations 5dofs --> 6dofs
       if(ndf.eq.5) iflag = 1                          ! iflag = 1 : glatte Schale
       if(ndf.eq.6) call boun63(k,ix,psid,iflag,ndf)  ! iflag = 0 : Kante  
csk       if(ndf.eq.6) call boun63(k,ix,m(n7),iflag,ndf)  ! iflag = 0 : Kante  
       if(ityp.eq.6) iflag = 0     
c        iflag= 1 !test for iga_shell
       if (iflag.eq.1) then
         if(ilin.eq.2) then
            call matmulf (Rnold(1:3,1:3,k),rin(1,1,k),3,3,3, rh1)
         end if
         call matmulf (rh1,ul(4,2*nen+k),3,2,1, rh2)
       end if
c.... compute increment of axial vector Omega since last iteration
      dwn(1:3,k) = rh2(1:3)         

c.... update of current nodal basis
       if(ilin.eq.0)then                   ! linear

         call matcop (rin(1,1,k),3,3, rwn(1,1,k))
         call vecp (dwn(1,k),rin(1,3,k),ddn(1,k))

       else                                ! nonlinear

         call updr63 (dwn(1,k),rh,ilin,1)
         call matmulf (rh,Rnold(1,1,k),3,3,3,rh1)
         call MATCOP(rh1,3,3,Rnold(1,1,k))
         call matmulf (rh1,rin(1,1,k),3,3,3, rwn(1,1,k))

       end if 

c....  T = W^T * H * T_3
c...   use H as derived         
       !call updr63 (dwn(1,k),twn(1,1,k),ilin,3)
c...   set H to unity
       call pzero(twn(1,1,k),9)
       twn(1,1,k) = 1.0d0
       twn(2,2,k) = 1.0d0
       twn(3,3,k) = 1.0d0
    
       if(ilin.eq.2)then                 !  finite rotations
         call skew (w(1,1,k),rwn(1,3,k))                      
         if (iflag.eq.1) then
           rh=0.d0
           call matmulf (twn(1,1,k),rwn(1,1,k),3,3,2, rh)          
c          call matmulf (twn(1,1,k),rwn(1,1,k),3,3,3, rh) !HIER NUR 2 SPALTEN???
           call matcop (rh,3,3, twn(1,1,k))
         end if

       else                              ! linear and moderate rotations
         call skew (w(1,1,k),rin(1,3,k))                      
         if (iflag.eq.1) call matcop (rin(1,1,k),3,3, twn(1,1,k))
       endif

       call mttmul (w(1,1,k),twn(1,1,k),3,3,3, rh)
       call matcop (rh,3,3, w(1,1,k))

         
100   continue
cwd   from here on only for ANS->not needed at the moment
c.... current director at midside nodes
      do 200 i = 1,3
        ran(i,1) = (rwn(i,3,4) + rwn(i,3,1))/2.d0
        ran(i,2) = (rwn(i,3,1) + rwn(i,3,2))/2.d0
        ran(i,3) = (rwn(i,3,2) + rwn(i,3,3))/2.d0
        ran(i,4) = (rwn(i,3,3) + rwn(i,3,4))/2.d0

        xan(i,1) = ( (xl(i,4)+ul(i,4))-(xl(i,1)+ul(i,1)) )/2.d0
        xan(i,2) = ( (xl(i,2)+ul(i,2))-(xl(i,1)+ul(i,1)) )/2.d0
        xan(i,3) = ( (xl(i,3)+ul(i,3))-(xl(i,2)+ul(i,2)) )/2.d0
        xan(i,4) = ( (xl(i,3)+ul(i,3))-(xl(i,4)+ul(i,4)) )/2.d0

        ran0(i,1) = (rin(i,3,4) + rin(i,3,1))/2.d0
        ran0(i,2) = (rin(i,3,1) + rin(i,3,2))/2.d0
        ran0(i,3) = (rin(i,3,2) + rin(i,3,3))/2.d0
        ran0(i,4) = (rin(i,3,3) + rin(i,3,4))/2.d0

        xan0(i,1) = ( xl(i,4) - xl(i,1) )/2.d0
        xan0(i,2) = ( xl(i,2) - xl(i,1) )/2.d0
        xan0(i,3) = ( xl(i,3) - xl(i,2) )/2.d0
        xan0(i,4) = ( xl(i,3) - xl(i,4) )/2.d0

       if(ilin.eq.0)then 

        uan(i,1) = (ul(i,4) - ul(i,1))/2.d0
        uan(i,2) = (ul(i,2) - ul(i,1))/2.d0
        uan(i,3) = (ul(i,3) - ul(i,2))/2.d0
        uan(i,4) = (ul(i,3) - ul(i,4))/2.d0

        dan(i,1) = (ddn(i,4) + ddn(i,1))/2.d0
        dan(i,2) = (ddn(i,1) + ddn(i,2))/2.d0
        dan(i,3) = (ddn(i,2) + ddn(i,3))/2.d0
        dan(i,4) = (ddn(i,3) + ddn(i,4))/2.d0
       end if

200   continue 

c
c.... shear strains at assumed strain points and factors b1m, b2l, 
      do 220 k = 1,4
        mm = iabs(im(k))
        ll = iabs(il(k))
        if(ilin.eq.0)then
         gam(k) = dot(uan(1,k),ran0(1,k),3) + dot(xan0(1,k),dan(1,k),3)
         call mttmul (w(1,1,k),xan0(1,mm),3,3,1, bml(1,1,k))
         call mttmul (w(1,1,k),xan0(1,ll),3,3,1, bml(1,2,k))
        else
         gam(k) = dot(xan(1,k),ran(1,k),3) - dot(xan0(1,k),ran0(1,k),3)
         call mttmul (w(1,1,k),xan(1,mm),3,3,1, bml(1,1,k))
         call mttmul (w(1,1,k),xan(1,ll),3,3,1, bml(1,2,k))
        end if
220   continue 
c
      return
      end      
c
      subroutine updCP63 (xl,ul,dwn,ix,ndf,ndm,nen,ityp,ilin,
     2                    rinG,rinP,rwnGn,rwnGnp1,rwnPn,rwnPnp1)
      USE mdata
c-----------------------------------------------------------------------
c     update of nodal rotations in Control points
c     rinG    initial global nodal basis system
c     rinP    initial patch nodal basis system
c     rwnGnp1 current global nodal basis system of current iteration step
c     rwnGn   current global nodal basis system of last iteration step
c     rwnPnp1 current patch nodal basis system of current iteration step
c     rwnPn   current patch nodal basis system of last iteration step
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      include 'mdata.h'      
csk      common m(1)
      dimension  xl(ndm,*),ul(ndf,*),ix(*),
     1           dwn(3,*),rwnGn(3,3,*),rinG(3,3,*),
     2           rh(3,3),rh1(3,3),rh2(3),rwnGnp1(3,3,*),
     3           rinP(3,3,*),rwnPn(3,3,*),rwnPnp1(3,3,*)
c
      do 100  k = 1,nen
c....   initial nodal cartesian system
        call matcop (rinG(1,1,k),3,3, rh1)
        call matcop (ul(4,2*nen+k),2,1, rh2)          

c....   determination if 5dofs or 6dofs
        if(ndf.eq.5) iflag = 1                          ! iflag = 1 : glatte Schale
        if(ndf.eq.6) call boun63(k,ix,psid,iflag,ndf)  ! iflag = 0 : Kante  
csk        if(ndf.eq.6) call boun63(k,ix,m(n7),iflag,ndf)  ! iflag = 0 : Kante  

        if(ilin.eq.2) then
          call updr63 (dwn(1,k),rh,ilin,1)
          call matmulf (rh,rinG(1,1,k),3,3,3, rh1)
          call matmulf (rh,rinG(1,1,k),3,3,3, rwnGn(1,1,k))
          call matmulf (rh,rinP(1,1,k),3,3,3, rwnPn(1,1,k))
        end if
        rh2=0.d0
        if (iflag.eq.1) call matmulf (rh1,ul(4,2*nen+k),3,2,1, rh2)
        if (iflag.eq.0) call matmulf (rh1,ul(4,2*nen+k),3,3,1, rh2)

c....   update current axial vector Omega    
        call matadd (rh2,dwn(1,k),3,1, dwn(1,k))          

c....   update of current nodal basis
        if(ilin.eq.0)then                   ! linear
          call matcop (rinG(1,1,k),3,3, rwnGn(1,1,k))
          call matcop (rinG(1,1,k),3,3, rwnGnp1(1,1,k))
          call matcop (rinP(1,1,k),3,3, rwnPn(1,1,k))
          call matcop (rinP(1,1,k),3,3, rwnPnp1(1,1,k))
        else                                ! nonlinear
          call updr63 (dwn(1,k),rh,ilin,1)
          call matmulf (rh,rinG(1,1,k),3,3,3, rwnGnp1(1,1,k))
          call matmulf (rh,rinP(1,1,k),3,3,3, rwnPnp1(1,1,k))
          if (ilin.eq.1) then
            call matmulf(rh,rinG(1,1,k),3,3,3, rwnGn(1,1,k))
            call matmulf(rh,rinP(1,1,k),3,3,3, rwnPn(1,1,k))
          end if
        end if 
100   continue
      return
      end
      
       subroutine updCP63OmegaMult (xl,ul,Rnold,ix,ndf,ndm,nen,ityp,
     2                    ilin,rinG,rinP,rwnGn,rwnGnp1,rwnPn,rwnPnp1,
     3                    ddwn)
      USE mdata
c-----------------------------------------------------------------------
c     update of nodal rotations in Control points
c     rinG    initial global nodal basis system
c     rinP    initial patch nodal basis system
c     rwnGnp1 current global nodal basis system of current iteration step
c     rwnGn   current global nodal basis system of last iteration step
c     rwnPnp1 current patch nodal basis system of current iteration step
c     rwnPn   current patch nodal basis system of last iteration step
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      include 'mdata.h'      
csk      common m(1)
      dimension  xl(ndm,*),ul(ndf,*),ix(*),
     1           Rnold(1:3,1:3,nen),rwnGn(3,3,*),rinG(3,3,*),
     2           rh(3,3),rh1(3,3),rh2(3),rwnGnp1(3,3,*),
     3           rinP(3,3,*),rwnPn(3,3,*),rwnPnp1(3,3,*),
     4           beta1(3),beta2(3),temp(3),beta_new(3),
     5           ddwn(3,*)
c
      do 100  k = 1,nen
c....   initial nodal cartesian system
        call matcop (rinG(1,1,k),3,3, rh1)
        call matcop (ul(4,2*nen+k),2,1, rh2)          

c....   determination if 5dofs or 6dofs
        if(ndf.eq.5) iflag = 1                          ! iflag = 1 : glatte Schale
csk     m(n7) => psid  
        if(ndf.eq.6) call boun63(k,ix,psid,iflag,ndf)  ! iflag = 0 : Kante  

        if(ilin.eq.2) then
          call matmulf (Rnold(1:3,1:3,k),rinG(1,1,k),3,3,3, rh1)
          call matcop (rh1,3,3,rwnGn(1,1,k))
          call matmulf (Rnold(1:3,1:3,k),rinP(1,1,k),3,3,3,rwnPn(1,1,k))
        end if
        rh2=0.d0
        if (iflag.eq.1) call matmulf (rh1,ul(4,2*nen+k),3,2,1, rh2)
        if (iflag.eq.0) call matmulf (rh1,ul(4,2*nen+k),3,3,1, rh2)
        
c...    only incremental values are passed to integration points        
        ddwn(1:3,k) = rh2(1:3)

c....   update of current nodal basis
        if(ilin.eq.0)then                   ! linear
          call matcop (rinG(1,1,k),3,3, rwnGn(1,1,k))
          call matcop (rinG(1,1,k),3,3, rwnGnp1(1,1,k))
          call matcop (rinP(1,1,k),3,3, rwnPn(1,1,k))
          call matcop (rinP(1,1,k),3,3, rwnPnp1(1,1,k))
        else                                ! nonlinear
          call updr63 (rh2,rh,ilin,1)
          call matmulf (rh,Rnold(1,1,k),3,3,3, rh1)
          call matmulf (rh1,rinG(1,1,k),3,3,3, rwnGnp1(1,1,k))
          call matmulf (rh1,rinP(1,1,k),3,3,3, rwnPnp1(1,1,k))
          call matcop (rh1,3,3,Rnold(1,1,k))
          if (ilin.eq.1) then
            call matcop(rwnGnp1(1,1,k),3,3, rwnGn(1,1,k))
            call matcop(rwnPnp1(1,1,k),3,3, rwnPn(1,1,k))
          end if
        end if 
100   continue
      return
      end
      
      subroutine updCP63Quart (xl,ul,dwn,ix,ndf,ndm,nen,ityp,ilin,
     2                    rinG,rinP,rwnGn,rwnGnp1,rwnPn,rwnPnp1,
     3                    dwnold,ddwn,iRupd)
c-----------------------------------------------------------------------
c     update of nodal rotations in Control points
c     rinG    initial global nodal basis system
c     rinP    initial patch nodal basis system
c     rwnGnp1 current global nodal basis system of current iteration step
c     rwnGn   current global nodal basis system of last iteration step
c     rwnPnp1 current patch nodal basis system of current iteration step
c     rwnPn   current patch nodal basis system of last iteration step
c-----------------------------------------------------------------------
      USE mdata      
      implicit double precision (a-h,o-z)
csk      common m(1)
      dimension  xl(ndm,*),ul(ndf,*),ix(*),
     1           dwn(3,*),rwnGn(3,3,*),rinG(3,3,*),
     2           rh(3,3),rh1(3,3),rh2(3),rwnGnp1(3,3,*),
     3           rinP(3,3,*),rwnPn(3,3,*),rwnPnp1(3,3,*),
     4           beta1(3),beta2(3),temp(3),beta_new(3),
     5           dwnold(3,*),ddwn(3,*)
      dimension test(3),trh1(3,3),trh2(3,3),trh3(3,3),trh4(3,3)
      real*8 tol
      tol=1.0d-13
c
      do 100  k = 1,nen
c....   store value of rotation vector from last step
        dwnold(1:3,k) = dwn(1:3,k)
c....   initial nodal cartesian system
        call matcop (rinG(1,1,k),3,3, rh1)
        call matcop (ul(4,2*nen+k),2,1, rh2)          

c....   determination if 5dofs or 6dofs
        if(ndf.eq.5) iflag = 1                          ! iflag = 1 : glatte Schale
csk     m(n7) => psid  
        if(ndf.eq.6) call boun63(k,ix,psid,iflag,ndf)  ! iflag = 0 : Kante  

        if(ilin.eq.2) then
          call updr63 (dwn(1,k),rh,ilin,1)
          call matmulf (rh,rinG(1,1,k),3,3,3, rh1)
          call matmulf (rh,rinG(1,1,k),3,3,3, rwnGn(1,1,k))
          call matmulf (rh,rinP(1,1,k),3,3,3, rwnPn(1,1,k))
        end if
        rh2=0.d0
        if (iflag.eq.1) call matmulf (rh1,ul(4,2*nen+k),3,2,1, rh2)
        if (iflag.eq.0) call matmulf (rh1,ul(4,2*nen+k),3,3,1, rh2)
        
        ddwn(1:3,k) = rh2(1:3)

        
        if (ilin.eq.2.and.iRupd.eq.3) then
        !if (ilin.eq.2) then
        
c....     tests
c          test(1:3) = dwn(1:3,k) + ddwn(1:3,k)      
c          call updr63 (dwn(1:3,k),trh1,ilin,1)
c          call updr63 (ddwn(1:3,k),trh2,ilin,1)
c          call matmulf (trh2,trh1,3,3,3,trh3)
        
        
c....     update current axial vector Omega using QUARTERNIONS
c....     taken from Sansour and Wagner 2003:
c....     Multiplicative updating of the rotation tensor in the finite element analysis of rods and shells – a path independent approach
c....     dwn(1:3,k) contains omega_I^(n-1)
c....     rh2(1:3) contains delta omega_I
          domega = dsqrt(dot(rh2,rh2,3))
          omega  = dsqrt(dot(dwn(1,k),dwn(1,k),3))
          alpha1 = dcos(domega/2.0d0)
          alpha2 = dcos(omega/2.0d0)
          if (dabs(domega).gt.tol) then
            beta1(1:3) = dsin(domega/2.0d0)/domega*rh2(1:3)
          else
            beta1(1:3) = 0.5d0*rh2(1:3)
          end if
          if (dabs(omega).gt.tol) then
            beta2(1:3) = dsin(omega/2.0d0)/omega*dwn(1:3,k)
          else
            beta2(1:3) = 0.5d0*dwn(1:3,k)
          end if
          call vcross(beta1,beta2,temp)
          beta_new(1:3) = alpha1*beta2(1:3)+alpha2*beta1(1:3)-temp(1:3)
c         TEST EIGENES UPDATE          
          alpha_new = alpha1*alpha2 - dot(beta1,beta2,3)
          if (dabs(alpha_new).gt.1.0d0) then
            write (*,*) 'not allowed: ERROR in updCP63Quart'
          end if
          omega_new = 2.0d0 * dacos(alpha_new)
          if (dabs(omega_new).gt.tol) then
            dwn(1:3,k)=beta_new(1:3) * omega_new / dsin(omega_new/2.0d0)
          else
            dwn(1:3,k) = 2.0d0*beta_new(1:3)
          end if
          
c         VERSION VON SANSOUR WAGNER          
          
c          beta_new_n = dsqrt(dot(beta_new,beta_new,3))
c          if (dabs(beta_new_n).gt.tol) then
c            if (dabs(beta_new_n).gt.1.0d0) then
c            write (*,*) 'not allowed'
c            end if
c            dwn(1:3,k) = beta_new(1:3)*2.0d0*dasin(beta_new_n)/
c     1           beta_new_n
c          else
c            dwn(1:3,k) = 2.0d0*beta_new(1:3)
c          end if
          
c          call updr63 (dwn(1:3,k),trh4,ilin,1)
c          call updr63 (test,trh2,ilin,1)
c         test
c          dwn(1:3,k) = test
          
        else
c....     standard additive update          
          call matadd (rh2,dwn(1,k),3,1, dwn(1,k))          
        end if
     
        
c....   update of current nodal basis
        if(ilin.eq.0)then                   ! linear
          call matcop (rinG(1,1,k),3,3, rwnGn(1,1,k))
          call matcop (rinG(1,1,k),3,3, rwnGnp1(1,1,k))
          call matcop (rinP(1,1,k),3,3, rwnPn(1,1,k))
          call matcop (rinP(1,1,k),3,3, rwnPnp1(1,1,k))
        else                                ! nonlinear
          call updr63 (dwn(1,k),rh,ilin,1)
          call matmulf (rh,rinG(1,1,k),3,3,3, rwnGnp1(1,1,k))
          call matmulf (rh,rinP(1,1,k),3,3,3, rwnPnp1(1,1,k))
          if (ilin.eq.1) then
            call matmulf(rh,rinG(1,1,k),3,3,3, rwnGn(1,1,k))
            call matmulf(rh,rinP(1,1,k),3,3,3, rwnPn(1,1,k))
          end if
        end if 
100   continue
      return
      end
c
      subroutine updCP63Rupd2 (xl,ul,Rnold,ix,ndf,ndm,nen,ityp,ilin,
     2                    rinG,rinP,rwnGn,rwnGnp1,rwnPn,rwnPnp1)
c-----------------------------------------------------------------------
c     update of nodal rotations in Control points
c     rinG    initial global nodal basis system
c     rinP    initial patch nodal basis system
c     rwnGnp1 current global nodal basis system of current iteration step
c     rwnGn   current global nodal basis system of last iteration step
c     rwnPnp1 current patch nodal basis system of current iteration step
c     rwnPn   current patch nodal basis system of last iteration step
c-----------------------------------------------------------------------
      USE mdata
      implicit double precision (a-h,o-z)
csk      common m(1)
      dimension  xl(ndm,*),ul(ndf,*),ix(*),
     1           Rnold(1:3,1:3,nen),rwnGn(3,3,*),rinG(3,3,*),
     2           rh(3,3),rh1(3,3),rh2(3),rwnGnp1(3,3,*),
     3           rinP(3,3,*),rwnPn(3,3,*),rwnPnp1(3,3,*)
c
      do 100  k = 1,nen
c....   initial nodal cartesian system
        call matcop (rinG(1,1,k),3,3, rh1)
        call matcop (ul(4,2*nen+k),2,1, rh2)          

c....   determination if 5dofs or 6dofs
        if(ndf.eq.5) iflag = 1                          ! iflag = 1 : glatte Schale
csk     m(n7) => psid  
        if(ndf.eq.6) call boun63(k,ix,psid,iflag,ndf)  ! iflag = 0 : Kante  

        if(ilin.eq.2) then
          call matmulf (Rnold(1:3,1:3,k),rinG(1,1,k),3,3,3, rh1)
c          call updr63 (dwn(1,k),rh,ilin,1)
c          call matmulf (rh,rinG(1,1,k),3,3,3, rh1)
          call matcop (rh1,3,3,rwnGn(1,1,k))
          call matmulf (Rnold(1:3,1:3,k),rinP(1,1,k),3,3,3,rwnPn(1,1,k))
        end if
        rh2=0.d0
        if (iflag.eq.1) call matmulf (rh1,ul(4,2*nen+k),3,2,1, rh2)
        if (iflag.eq.0) call matmulf (rh1,ul(4,2*nen+k),3,3,1, rh2)

c....   update current axial vector Omega    
c        call matadd (rh2,dwn(1,k),3,1, dwn(1,k))          

c....   update of current nodal basis
        if(ilin.eq.0)then                   ! linear
          call matcop (rinG(1,1,k),3,3, rwnGn(1,1,k))
          call matcop (rinG(1,1,k),3,3, rwnGnp1(1,1,k))
          call matcop (rinP(1,1,k),3,3, rwnPn(1,1,k))
          call matcop (rinP(1,1,k),3,3, rwnPnp1(1,1,k))
        else                                ! nonlinear
          call updr63 (rh2,rh,ilin,1)
          call matmulf (rh,Rnold(1,1,k),3,3,3, rh1)
          call matmulf (rh1,rinG(1,1,k),3,3,3, rwnGnp1(1,1,k))
          call matmulf (rh1,rinP(1,1,k),3,3,3, rwnPnp1(1,1,k))
          call matcop (rh1,3,3,Rnold(1,1,k))
          if (ilin.eq.1) then
            call matcop(rwnGnp1(1,1,k),3,3, rwnGn(1,1,k))
            call matcop(rwnPnp1(1,1,k),3,3, rwnPn(1,1,k))
          end if
        end if 
100   continue
      return
      end
c      
      subroutine computeCPfromGP63 (xl,ul,dwg,ix,ilin,rinG,rinP,lint,
     1                              rwnGn,rwnGnp1,rwnPn,rwnPnp1,
     2                            sg,tg,ndm,nen,ndf,
     3                            ni,nj,NURp,NURq,AInkv1,AInkv2,
     4                            Iexactnormal,AInbasissys,ngb,n ) 
      USE mdata
      USE iofile
c-----------------------------------------------------------------------
c     Compute nodal rotation vector for control points from Gauss point
c     rotation vector by solving a system of equations      
c     rinG    initial global nodal basis system
c     rinP    initial patch nodal basis system
c     rwnGnp1 current global nodal basis system of current iteration step
c     rwnGn   current global nodal basis system of last iteration step
c     rwnPnp1 current patch nodal basis system of current iteration step
c     rwnPn   current patch nodal basis system of last iteration step
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      include 'iofile.h'   
csk      common m(1)
      dimension rh0(3,3,nen),rh1(3,3,nen),rh2(3,nen)
      dimension rh(3,3),rh3(3,3),rhtest(3,3)
      dimension t0exact(3,3,3),tg(*),sg(*),dwg(9,*)
      dimension shp((NURp+1)*(NURq+1),lint)
      dimension xl(ndm,*),ul(ndf,*),ix(*),iMrow(nen),iMcol(nen)
      dimension rinG(3,3,nen),rinP(3,3,nen),rwnPn(3,3,nen),
     1          rwnPnp1(3,3,nen),rwnGn(3,3,nen),rwnGnp1(3,3,nen)
      real*8 MatM(3,3,nen), dbeta(3)
      real*8 matN(nen,nen),matN2(nen,nen),matB(nen,9)
      integer flag
      integer,pointer,dimension(:)   :: ipiv
      
    
      if (ilin.eq.0) then
        do k=1,nen
          call matcop (rinG(1,1,k),3,3, rwnGn(1,1,k))
          call matcop (rinG(1,1,k),3,3, rwnGnp1(1,1,k))
          call matcop (rinP(1,1,k),3,3, rwnPn(1,1,k))
          call matcop (rinP(1,1,k),3,3, rwnPnp1(1,1,k))
        end do
      end if
c     CAREFUL: #integration points MUST BE EQUAL (p+1)*(q+1)
      if (lint.ne.((NURp+1)*(NURq+1))) then
       write(*,*) 'ERROR: #integration points MUST BE EQUAL (p+1)*(q+1)'
      end if
c     initialize system of equations
      
c     set number of nodal rotations      
      do l= 1,nen
        if (ndf.eq.5) then
          iMrow(l) = 2
          iMcol(l) = 2
        else if (ndf.eq.6) then
csk       m(n7) => psid  
          call boun63(l,ix,psid,iflag,ndf)  ! iflag = 0 : Kante  
csk          call boun63(l,ix,m(n7),iflag,ndf) 
            if (iflag.eq.1) then
              iMrow(l) = 2
              iMcol(l) = 2
            else if (iflag.eq.0) then
              iMrow(l) = 2  !oder doch eher 2????
              iMcol(l) = 3
            end if

        else
          write(*,*) 'Wrong value for ndf: choose 5 or 6!!'
          write(iow,*) 'Wrong value for ndf: choose 5 or 6!!'
        end if
      end do
c     loop over all integration points to fill system      
      do l = 1,lint
        xsi = sg(l)
        eta = tg(l)
        call shap63onlyBasisFun(xsi,eta,xl,shp(1:nen,l),ndm,nen,ni,nj,
     +            NURp,NURq,AInkv1,AInkv2) 

        if (Iexactnormal.eq.1) then
          call IGAreadNormalm(AInbasissys,numel,ngb**2,t0exact,l,n)
          rh0(1:3,1:3,l)   = t0exact(1:3,3,1:3)
          rh1(1:3,1:3,l)   = t0exact(1:3,3,1:3)
        else
          do i = 1,3
            do j = 1,3
              rh1(i,j,l) = 0.0d0
              do l1 = 1,nen
                rh1(i,j,l) = rh1(i,j,l) + rinP(i,j,l1) * shp(l1,l) 
              end do
            end do
          end do
          call matcop(rh1(1:3,1:3,l),3,3,rh0(1:3,1:3,l))
c          write(*,*) 'USE EXACT NORMAL MODE!!!'
        end if
        if(ilin.eq.2) then
c....     calculate R_n(->rh) from value omega_gp(<-dwg).
c....     Omega_gp is stored in the h3 field and loaded        
          call updr63 (dwg(1:3,l),rh,ilin,1)
          call matmulf (rh,rh1(1:3,1:3,l),3,3,3, rh3)
          call matcop  (rh3,3,3,rh1(1:3,1:3,l))
        end if
c       write to matrix and vector
        do i = 1,nen
          matN(l,i) = shp(i,l)
        end do
        do i = 1,3
          matB(l,i) = rh1(i,1,l)
        end do
        do i = 1,3
          matB(l,i+3) = rh1(i,2,l)
        end do
        do i = 1,3
          matB(l,i+6) = rh1(i,3,l)
        end do
      end do
c     save matrix matN for later usage
      call matcop(matN,nen,nen,matN2)
c     solve matrix to get control point values
      allocate(ipiv(1:nen))
      call dgesv(nen,9,matN,nen,ipiv,matB,nen,flag)
      if (flag.ne.0) write(*,*) '!!!  error in MKL solver in computeCPfr
     1omGP63 on element level !!!'
      deallocate(ipiv)
      
c     write to rwnGn and rwnPn
      do i = 1,nen
        rwnPn(1:3,1,i) = matB(i,1:3)
        rwnPn(1:3,2,i) = matB(i,4:6)
        rwnPn(1:3,3,i) = matB(i,7:9)
        if (iMcol(i).eq.2) then
          rwnGn(1:3,1,i) = matB(i,1:3)
          rwnGn(1:3,2,i) = matB(i,4:6)
          rwnGn(1:3,3,i) = matB(i,7:9)
        end if
      end do
c     test
      !do i = 1,3
      !  do j = 1,3
      !    rhtest(i,j) = 0.0d0
      !    do l1 = 1,nen
      !      rhtest(i,j) = rhtest(i,j) + rwnPn(i,j,l1) * shp(l1,l) 
      !    end do
      !  end do
      !end do

      rh2=0.d0
      do j = 1,nen
        do i = 1,iMcol(j)
          rh2(i,j) = ul(3+i,2*nen+j)
        end do
      end do
c     loop over integration points      
      do l =1,lint
        
c....   determine M_I for all nodes based on rwnGn and rh1
        MatM=0.d0
        do l1 = 1,nen
          do i = 1,iMrow(l1)
            do j = 1,iMcol(l1)
              MatM(i,j,l1)   = dot(rh1(1,i,l),rwnGn(1,j,l1),3)
            end do
          end do
        end do
c...    calculate delta_beta in Gauss point
c...    transformation with A matrix necessary 
c...    delta_beta = sum(N_I*A_I*delta_beta_I)
        dbeta = 0.0d0
        do l1 = 1,nen
          do i = 1,iMrow(l1)
            do j = 1,iMcol(l1)
c...          rh2(:,l) is delta_beta for node l   
              dbeta(i) = dbeta(i) +
     +                        shp(l1,l)*MatM(i,j,l1)*rh2(j,l1)
            end do
          end do
        end do
c...    calculate delta_omega in Gauss point  
c...    the resulting additional rotation delta_omega is stored in rh3        
        call matmulf (rh1(1:3,1:3,l),dbeta(1),3,3,1, rh3(1,1))
c....   compute current axial vector Omega (NO UPDATE; THIS IS LATERON IN UPDGP63
        call matadd (rh3(1,1),dwg(1,l),3,1, rh3(1,1))  
c....   update of current nodal basis
        if(ilin.eq.0) then   
        
        else
          call updr63 (rh3(1:3,1),rh,ilin,1)
          call matmulf (rh,rh0(1:3,1:3,l),3,3,3, rh1(1:3,1:3,l))
        end if
c       write to vector
        do i = 1,3
          matB(l,i) = rh1(i,1,l)
        end do
        do i = 1,3
          matB(l,i+3) = rh1(i,2,l)
        end do
        do i = 1,3
          matB(l,i+6) = rh1(i,3,l)
        end do      
      end do
c     solve matrix to get control point values
      allocate(ipiv(1:nen))
      call dgesv(nen,9,matN2,nen,ipiv,matB,nen,flag)
      if (flag.ne.0) write(*,*) '!!!  error in MKL solver in computeCPfr
     1omGP63 on element level !!!'
      deallocate(ipiv)
      
c     write to rwnGn and rwnPn
      do i = 1,nen
        rwnPnp1(1:3,1,i) = matB(i,1:3)
        rwnPnp1(1:3,2,i) = matB(i,4:6)
        rwnPnp1(1:3,3,i) = matB(i,7:9)
        if (iMcol(i).eq.2) then
          rwnGnp1(1:3,1,i) = matB(i,1:3)
          rwnGnp1(1:3,2,i) = matB(i,4:6)
          rwnGnp1(1:3,3,i) = matB(i,7:9)
        end if
      end do      

      
      return
      end
c
c      
      subroutine updGP63 (xl,ul,dwg,twn,dd,w,ix,ndf,ndm,nen,ityp,ilin,
     2                   k,shp,rd,r0,HutT_1,HutT_2,twn_1,twn_2,iMcol,
     4                   rinG,rinP,rwnGn,rwnGnp1,rwnPn,rwnPnp1,
     5                   t0exact,Iexactnormal)
      USE mdata
      USE iofile
c-----------------------------------------------------------------------
c     update of rotations in Gauss points
c     rinG    initial global nodal basis system
c     rinP    initial patch nodal basis system
c     rwnGnp1 current global nodal basis system of current iteration step
c     rwnGn   current global nodal basis system of last iteration step
c     rwnPnp1 current patch nodal basis system of current iteration step
c     rwnPn   current patch nodal basis system of last iteration step
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      include 'iofile.h'   
csk      common m(1)
      dimension  xl(ndm,*),ul(ndf,*),ix(*),
     1           dwg(9,*),twn(3,3,*),dd(3,3),
     2           w(3,3,*),rh(3,3),rh1(3,3),rh2(3,nen),
     4           shp(3,*),dbeta(3,3),
     5           rd(3,3),r0(3,3),rh3(3,3),Aini(3,3),
     6           drhd1(3,3),drhd2(3,3),rh4(3,3),iMrow(nen),iMcol(nen),
     7           rh5(3,3),rh6(3,3),twn_1(3,3,nen),twn_2(3,3,nen)
      dimension rinG(3,3,nen),rinP(3,3,nen),rwnPn(3,3,nen),
     1          rwnPnp1(3,3,nen),rwnGn(3,3,nen),rwnGnp1(3,3,nen)
      double precision MatH(3,3),MatH_1(3,3),MatH_2(3,3),
     1                 MatW(3,3),MatW_1(3,3),MatW_2(3,3),
     2                 MatT3(3,3),MatT3_1(3,3),MatT3_2(3,3),
     3                 MatM(3,3,nen),MatM_1(3,3,nen),MatM_2(3,3,nen),
     4                 MatT(3,3),MatWHT_1(3,3),MatWHT_2(3,3),
     5                 MatW_1HT(3,3),MatW_2HT(3,3),MatWH_1T(3,3),
     6                 MatWH_2T(3,3),MatT_1(3,3),MatT_2(3,3),
     7                 HutT_1(3,3,nen),HutT_2(3,3,nen),aid1(3,3),
     8                 aid2(3,3),bound,t0exact(3,3,3)
      real*8 rhtest(3,3)
c
      ismoothrow = 2 !in CMAME2014 ist hier 2 eingesetzt!
      bound=1.d-5
c.... initial cartesian system in Gauss point: rh1
c.... nodal rotations delta_beta_I (nodewise rotations): rh2
c      write(iow,*) dwg(1:9,k)
c     if loaded systems ar used
      if (Iexactnormal.eq.1) then
        r0    = t0exact(1:3,1:3,3)
        rh1   = t0exact(1:3,3,1:3)
        drhd1 = t0exact(1:3,1,1:3)
        drhd2 = t0exact(1:3,2,1:3)
        rd = 0.0d0
      else
        do i = 1,3
          do j = 1,3
            rh1(i,j) = 0.0d0
            rd(i,j)  = 0.d0
            r0(i,j)  = 0.d0
            drhd1(i,j) = 0.d0
            drhd2(i,j) = 0.d0
            do l = 1,nen
              rh1(i,j) = rh1(i,j) + rinP(i,j,l) * shp(3,l) 
              r0(i,j)  = r0 (i,j) + rinP(i,3,l) * shp(j,l)
              drhd1(i,j) = drhd1(i,j) +rinP(i,j,l) *shp(1,l)
              drhd2(i,j) = drhd2(i,j) +rinP(i,j,l) *shp(2,l)
            end do
          end do
        end do
      end if
c      if(ndf.eq.5) iflag = 1  
      do l= 1,nen
        if (ndf.eq.5) then
          iMrow(l) = ismoothrow !2
          iMcol(l) = 2
        else if (ndf.eq.6) then
csk       m(n7) => psid  
          call boun63(l,ix,psid,iflag,ndf)  ! iflag = 0 : Kante  
csk       call boun63(l,ix,m(n7),iflag,ndf) 
            if (iflag.eq.1) then
              iMrow(l) = ismoothrow !2
              iMcol(l) = 2 
            else if (iflag.eq.0) then
              iMrow(l) = 2  ! Hier ist im CMAME2014 2 eingesetzt, 3 funktioniert bei free form surface besser, 2 beim channel section beam!?!?
              iMcol(l) = 3
            end if

        else
          write(*,*) 'Wrong value for ndf: choose 5 or 6!!'
          write(iow,*) 'Wrong value for ndf: choose 5 or 6!!'
        end if
      end do

      
      rh2=0.d0
      do l = 1,nen
        do i = 1,iMcol(l)
          rh2(i,l) = ul(3+i,2*nen+l)
        end do
      end do
      call matcop  (rh1,3,3,Aini)

c      if(ndf.eq.5) iflag = 1                 ! iflag = 1 : glatte Schale
c                                            ! iflag = 0 : Kante  
      if(ilin.eq.2) then
c....   calculate R_n(->rh) from value omega_gp(<-dwg).
c....   Omega_gp is stored in the h3 field and loaded        
        call updr63 (dwg(1:3,k),rh,ilin,1)
        call matmulf (rh,rh1,3,3,3, rh3)
        call matcop  (rh3,3,3,rh1)
c...    1-direction
        call updRalpha63(dwg(1,k),rh3,ilin,1)
c...    R,alp*Ai
        call matmulf (rh3,Aini,3,3,3,rh5)
c...    R*R_0,1
        call matmulf (rh ,drhd1,3,3,3,rh4)
c...    R,alp*Ai + (R*R_0,1)*(R_0^T*Ai)
        call matadd  (rh5(1,1),rh4(1,1),3,1,aid1(1,1))
        call matadd  (rh5(1,2),rh4(1,2),3,1,aid1(1,2))
        call matadd  (rh5(1,3),rh4(1,3),3,1,aid1(1,3))
c...    2-direction
        call updRalpha63(dwg(1,k),rh3,ilin,2)
c...    R,alp*Ai
        call matmulf (rh3,Aini,3,3,3,rh5)
c...    R*R_0,2
        call matmulf (rh ,drhd2,3,3,3,rh4)
c...    R,alp*Ai + (R*R_0,1)*(R_0^T*Ai)
        call matadd  (rh5(1,1),rh4(1,1),3,1,aid2(1,1))
        call matadd  (rh5(1,2),rh4(1,2),3,1,aid2(1,2))
        call matadd  (rh5(1,3),rh4(1,3),3,1,aid2(1,3))
      else
c...    determine derivation of R_0(stored in rh1) for M_I,alpha
        call matcop(drhd1,3,3,aid1)
        call matcop(drhd2,3,3,aid2)
      end if
c....   determine M_I for all nodes based on rwn and rh1
      MatM=0.d0
      MatM_1=0.d0
      MatM_2=0.d0
      do l = 1,nen
c        if (iMcol(l).eq.3) then
          do i = 1,iMrow(l)
            do j = 1,iMcol(l)
              MatM(i,j,l)   = dot(rh1(1,i),rwnGn(1,j,l),3)
              MatM_1(i,j,l) = dot(aid1(1,i),rwnGn(1,j,l),3)
              MatM_2(i,j,l) = dot(aid2(1,i),rwnGn(1,j,l),3)
            end do
          end do
c        else
c        do i = 1,iMrow(l)
c            do j = 1,iMcol(l)
c        MatM(i,j,l) = dot(rh1(1,i),rwnGn(1,j,l),3)
c     -                   -( dot(rh1(1,3),rwnGn(1,j,l),3)
c     *                    * dot(rh1(1,i),rwnGn(1,3,l),3)
c     /                    / dot(rh1(1,3),rwnGn(1,3,l),3))
c              MatM_1(i,j,l) = dot(aid1(1,i),rwnGn(1,j,l),3)
c     -                    -((( dot(aid1(1,3),rwnGn(1,j,l),3)
c     *                        *dot(rh1(1,i),rwnGn(1,3,l),3)
c     +                        +dot(rh1(1,3),rwnGn(1,j,l),3)
c     *                        *dot(aid1(1,i),rwnGn(1,3,l),3))
c     *                        *dot(rh1(1,3),rwnGn(1,3,l),3)
c     -                      -( dot(aid1(1,3),rwnGn(1,3,l),3)
c     *                        *dot(rh1(1,3),rwnGn(1,j,l),3)
c     *                        *dot(rh1(1,i),rwnGn(1,3,l),3)))
c     /                     /(dot(rh1(1,3),rwnGn(1,3,l),3)**2))    
c              MatM_2(i,j,l) = dot(aid2(1,i),rwnGn(1,j,l),3)
c     -                    -((( dot(aid2(1,3),rwnGn(1,j,l),3)
c     *                        *dot(rh1(1,i),rwnGn(1,3,l),3)
c     +                        +dot(rh1(1,3),rwnGn(1,j,l),3)
c     *                        *dot(aid2(1,i),rwnGn(1,3,l),3))
c     *                        *dot(rh1(1,3),rwnGn(1,3,l),3)
c     -                      -( dot(aid2(1,3),rwnGn(1,3,l),3)
c     *                        *dot(rh1(1,3),rwnGn(1,j,l),3)
c     *                        *dot(rh1(1,i),rwnGn(1,3,l),3)))
c     /                     /(dot(rh1(1,3),rwnGn(1,3,l),3)**2))
c           end do
c        end do
c        end if
      end do

c     test
      do it = 1,3
        do jt = 1,3
          rhtest(it,jt) = 0.0d0
          do l1 = 1,nen
            rhtest(it,jt) = rhtest(it,jt) + rwnPn(it,jt,l1) * shp(3,l1) 
          end do
        end do
      end do        
      
c...  calculate delta_beta in Gauss point
c...  transformation with A matrix necessary 
c...  delta_beta = sum(N_I*A_I*delta_beta_I)
      dbeta = 0.0d0
      do l = 1,nen
        do idev = 1,3
          do i = 1,iMrow(l)
            do j = 1,iMcol(l)
c...          rh2(:,l) is delta_beta for node l   
              dbeta(i,idev) = dbeta(i,idev) +
     +                        shp(idev,l)*MatM(i,j,l)*rh2(j,l)
              if (idev.eq.1) dbeta(i,idev) = dbeta(i,idev) +
     +                        shp(3,l)*MatM_1(i,j,l)*rh2(j,l)
              if (idev.eq.2) dbeta(i,idev) = dbeta(i,idev) +
     +                        shp(3,l)*MatM_2(i,j,l)*rh2(j,l)
            end do
          end do
        end do    
      end do
c...  calculate delta_omega in Gauss point
c...  with rh1 being the current coordinate system in current GP
c...  and the additional rotations delta_beta in current GP
c...  the resulting additional rotation delta_omega i stored in rh2
      do i = 1,3
        call matmulf (rh1,dbeta(1,i),3,3,1, rh2(1,i))
      end do
      do i = 1,3
        rh2(i,1) = rh2(i,1)+aid1(i,1)*dbeta(1,3)+aid1(i,2)*dbeta(2,3)
     +             +aid1(i,3)*dbeta(3,3)
        rh2(i,2) = rh2(i,2)+aid2(i,1)*dbeta(1,3)+aid2(i,2)*dbeta(2,3)
     +             +aid2(i,3)*dbeta(3,3)
      end do

c.... update current axial vector Omega    
      call matadd (rh2(1,3),dwg(1,k),3,1, dwg(1,k))  
      call matadd (rh2(1,1),dwg(4,k),3,1, dwg(4,k))     
      call matadd (rh2(1,2),dwg(7,k),3,1, dwg(7,k))             

c.... update of current nodal basis
      if(ilin.eq.0) then                   ! linear, correct!!!

         call matcop (r0,3,3, rd)
         call vecp (dwg(1,k),r0(1,3),dd(1,3))
         
         call vecp (dwg(4,k),r0(1,3),dd(1,1))
         call vecp (dwg(1,k),r0(1,1),rh3)
         call matadd(dd(1,1),rh3,3,1,dd(1,1))
         
         call vecp (dwg(7,k),r0(1,3),dd(1,2))
         call vecp (dwg(1,k),r0(1,2),rh3)
         call matadd(dd(1,2),rh3,3,1,dd(1,2))

         call matcop(drhd1(1,1),3,3,MatT3_1(1,1))
         call matcop(drhd2(1,1),3,3,MatT3_2(1,1))


      else                                ! nonlinear
c....   rd has to be computed
        call updr63 (dwg(1:3,k),rh,ilin,1)
c....   d = R * D; ->rh maps initial to current tripod 
        call matmulf (rh,Aini(1,3),3,3,1,rd(1,3))
c...    for small rotations:
c....   d,alph = skew(omega,alpha) * D + delta_R * D,alpha
c....   for finite rotations:ß
c....   ai,alph=R,alp*Ai+R*R_0,alp*R_0^T*Ai
c....   rh3: R,alph
c....   Aini(1,i): Ai (D=A3) (=R_0)
c....   rh: R
c....   drhdalpha:R_0,alpha
c....   rh4: R,0
c
c
c....   1-direction
        call updRalpha63(dwg(1,k),rh3,ilin,1)
c...    R,alp*Ai
        call matmulf(rh3,Aini,3,3,3,rh5)
c...    R*R_0,1
        call matmulf (rh ,drhd1,3,3,3,rh4)
c...    R,alp*Ai + (R*R_0,1)*(R_0^T*Ai)
        call matadd  (rh5(1,3),rh4(1,3),3,1,rd(1,1))
        call matadd  (rh5(1,1),rh4(1,1),3,1,MatT3_1(1,1))
        call matadd  (rh5(1,2),rh4(1,2),3,1,MatT3_1(1,2))
        call matadd  (rh5(1,3),rh4(1,3),3,1,MatT3_1(1,3))
        call matadd  (rh5(1,1),rh4(1,1),3,1,aid1(1,1))
        call matadd  (rh5(1,2),rh4(1,2),3,1,aid1(1,2))
        call matadd  (rh5(1,3),rh4(1,3),3,1,aid1(1,3))

c....   2-direction
        call updRalpha63(dwg(1,k),rh3,ilin,2)
c...    R,alp*Ai
        call matmulf (rh3,Aini,3,3,3,rh5)
c....   R*R_0,2
        call matmulf (rh ,drhd2,3,3,3,rh4)
c...    R,alp*D + (R*R_0,1)*(R_0^T*D)
        call matadd  (rh5(1,3),rh4(1,3),3,1,rd(1,2))
        call matadd  (rh5(1,1),rh4(1,1),3,1,MatT3_2(1,1))
        call matadd  (rh5(1,2),rh4(1,2),3,1,MatT3_2(1,2))
        call matadd  (rh5(1,3),rh4(1,3),3,1,MatT3_2(1,3))
        call matadd  (rh5(1,1),rh4(1,1),3,1,aid2(1,1))
        call matadd  (rh5(1,2),rh4(1,2),3,1,aid2(1,2))
        call matadd  (rh5(1,3),rh4(1,3),3,1,aid2(1,3))

c....   Gauss point tripod has to be updated 
        if (ilin.eq.2)   rh1= matmul (rh,Aini)

        
c....   determine M_I for all nodes based on rwnn and current a_i
        do l = 1,nen
c          if (iMcol(l).eq.3) then
            do i = 1,iMrow(l)
              do j = 1,iMcol(l)
                MatM(i,j,l)   = dot(rh1(1,i),rwnGnp1(1,j,l),3)
                MatM_1(i,j,l) = dot(aid1(1,i),rwnGnp1(1,j,l),3)
                MatM_2(i,j,l) = dot(aid2(1,i),rwnGnp1(1,j,l),3)
              end do
            end do
c          else
c          do i = 1,iMrow(l)
c              do j = 1,iMcol(l)
c          MatM(i,j,l) = dot(rh1(1,i),rwnGnp1(1,j,l),3)
c     -                     -( dot(rh1(1,3),rwnGnp1(1,j,l),3)
c     *                      * dot(rh1(1,i),rwnGnp1(1,3,l),3)
c     /                      / dot(rh1(1,3),rwnGnp1(1,3,l),3))
c                MatM_1(i,j,l) = dot(aid1(1,i),rwnGnp1(1,j,l),3)
c     -                      -((( dot(aid1(1,3),rwnGnp1(1,j,l),3)
c     *                           *dot(rh1(1,i),rwnGnp1(1,3,l),3)
c     +                          +dot(rh1(1,3),rwnGnp1(1,j,l),3)
c     *                          *dot(aid1(1,i),rwnGnp1(1,3,l),3))
c     *                          *dot(rh1(1,3),rwnGnp1(1,3,l),3)
c     -                        -( dot(aid1(1,3),rwnGnp1(1,3,l),3)
c     *                          *dot(rh1(1,3),rwnGnp1(1,j,l),3)
c     *                          *dot(rh1(1,i),rwnGnp1(1,3,l),3)))
c     /                       /(dot(rh1(1,3),rwnGnp1(1,3,l),3)**2))    
c                MatM_2(i,j,l) = dot(aid2(1,i),rwnGnp1(1,j,l),3)
c     -                      -((( dot(aid2(1,3),rwnGnp1(1,j,l),3)
c     *                          *dot(rh1(1,i),rwnGnp1(1,3,l),3)
c     +                          +dot(rh1(1,3),rwnGnp1(1,j,l),3)
c     *                          *dot(aid2(1,i),rwnGnp1(1,3,l),3))
c     *                          *dot(rh1(1,3),rwnGnp1(1,3,l),3)
c     -                        -( dot(aid2(1,3),rwnGnp1(1,3,l),3)
c     *                          *dot(rh1(1,3),rwnGnp1(1,j,l),3)
c     *                          *dot(rh1(1,i),rwnGnp1(1,3,l),3)))
c     /                       /(dot(rh1(1,3),rwnGnp1(1,3,l),3)**2)) 
c          end do
c            end do
c          end if
        end do
      end if 
c
c...  T = W^T * H * T_3 and
c...  T_,alpha = W_,a^T*H*T3+W^T*H_,a*T3+W^T*H*T3_,a
c
c...  create T3, T3_,alpha is done above
      MatT3(1:3,1) = rh1(1:3,1)
      MatT3(1:3,2) = rh1(1:3,2)
      MatT3(1:3,3) = rh1(1:3,3)
      if (ilin.eq.1) then
        call matcop(aid1(1,1),3,3,MatT3_1(1,1))
        call matcop(aid2(1,1),3,3,MatT3_2(1,1))
      end if
c
c...  compute H and H_,alpha (for ilin.ne.2 H=1)
      call updr63 (dwg(1:3,k),MatH,ilin,3)    !HIER?????? CHECKEN ob dwg(1,k) richtigß???
      call updHalpha63(dwg(1,k),MatH_1,ilin,1)
      call updHalpha63(dwg(1,k),MatH_2,ilin,2)

c...  compute W and W_,alpha
      if (ilin.eq.2) then
        call skew(MatW,rd(1,3))
        call skew(MatW_1,rd(1,1))
        call skew(MatW_2,rd(1,2))
      else
        call skew(MatW,r0(1,3))
        call skew(MatW_1,r0(1,1))
        call skew(MatW_2,r0(1,2))
      end if
c
c...  calculate T = W^T * H * T_3
      call mttmul (MatW,MatH,3,3,3, rh3)
      call matmulf(rh3,MatT3,3,3,3,MatT)
c...  calculate TM_I = W^T * H *T_3 *M_I and store it to w(1,1,I)
      call pzero(w,9*nen)
      do i = 1,nen
c        call matmulf (MatT,MatM(1,1,i),3,iMrow(i),iMrow(i),w(1,1,i))
        call matmulf (MatT,MatM(1,1,i),3,3,3,w(1,1,i))
      end do
c
c...  calculate W^T * H * T_3,a
      call matmulf(rh3,MatT3_1,3,3,3,MatWHT_1)
      call matmulf(rh3,MatT3_2,3,3,3,MatWHT_2)
c...  calculate W_,a^T*H*T3
      call mttmul (MatW_1,MatH,3,3,3, rh3)
      call matmulf(rh3,MatT3,3,3,3, MatW_1HT)
      call mttmul (MatW_2,MatH,3,3,3, rh3)
      call matmulf(rh3,MatT3,3,3,3, MatW_2HT)
c...  calculate W^T*H_,a*T3
      call mttmul (MatW,MatH_1,3,3,3,rh3)
      call matmulf(rh3,MatT3,3,3,3,MatWH_1T)
      call mttmul (MatW,MatH_2,3,3,3,rh3)
      call matmulf(rh3,MatT3,3,3,3,MatWH_2T)
c
c...  sum up above terms for T_,alpha
      call matadd(MatWHT_1,MatW_1HT,3,3,rh3)
      call matadd(rh3,MatWH_1T,3,3,MatT_1)
      call matadd(MatWHT_2,MatW_2HT,3,3,rh3)
      call matadd(rh3,MatWH_2T,3,3,MatT_2)
c
c...  calculate ^T_I,alpha and store it to HutT_,alpha
      HutT_1=0.d0
      HutT_2=0.d0
      do i = 1,nen
c...    1-direction
        call matmulf(MatT_1,MatM(1,1,i),3,3,iMcol(i),rh3)
        call msmul(rh3,3,iMcol(i),shp(3,i))
        call matcop(w(1,1,i),3,iMcol(i),rh4)
        call msmul(rh4,3,iMcol(i),shp(1,i))
        call matmulf(MatT,MatM_1(1,1,i),3,3,iMcol(i),rh5)
        call msmul(rh5,3,iMcol(i),shp(3,i))
        call matadd(rh3,rh4,3,iMcol(i),rh6)
        call matadd(rh6,rh5,3,iMcol(i),HutT_1(1,1,i))
!        call matcop(rh4,3,2,HutT_1(1,1,i))
c...    2-direction
        call matmulf(MatT_2,MatM(1,1,i),3,3,iMcol(i),rh3)
        call msmul(rh3,3,iMcol(i),shp(3,i))
        call matcop(w(1,1,i),3,iMcol(i),rh4)
        call msmul(rh4,3,iMcol(i),shp(2,i))
        call matmulf(MatT,MatM_2(1,1,i),3,3,iMcol(i),rh5)
        call msmul(rh5,3,iMcol(i),shp(3,i))
        call matadd(rh3,rh4,3,iMcol(i),rh6)
        call matadd(rh6,rh5,3,iMcol(i),HutT_2(1,1,i))
!        call matcop(rh4,3,2,HutT_2(1,1,i))
      end do
c
c...  for G_ßß part of geometrical matrix
c...  compute T_HTM_I = H * T_3 * M_I and store it to twn(1,1,I) 
      call matmulf(MatH,MatT3,3,3,3,rh)
      call pzero(twn,9*nen)
      call pzero(twn_1,9*nen)
      call pzero(twn_2,9*nen)
      do i = 1,nen
        call matmulf (rh,MatM(1,1,i),3,3,iMcol(i),twn(1,1,i))
      end do
c...  compute T_HTM_I,alpha and store it to twn_alpha(1,1,I)
c... 1-direction
      call matmulf(MatH_1,MatT3,3,3,3,rh3)
      call matmulf(MatH,MatT3_1,3,3,3,rh4)
      call matadd (rh3,rh4,3,3,rh5)
      do i= 1,nen
        rh3=0.d0
        do j = 1,iMrow(i)
          do l = 1,iMcol(i)
            rh3(j,l) = shp(1,i)*MatM(j,l,i)+shp(3,i)*MatM_1(j,l,i)
          end do
        end do
        call matmulf(rh5,MatM(1,1,i),3,3,iMcol(i),rh4)
        call msmul(rh4,3,iMcol(i),shp(3,i))
        call matmulf(rh,rh3,3,3,iMcol(i),rh6)
        call matadd(rh4,rh6,3,iMcol(i),twn_1(1,1,i))
      end do
c... 2-direction
      call matmulf(MatH_2,MatT3,3,3,3,rh3)
      call matmulf(MatH,MatT3_2,3,3,3,rh4)
      call matadd (rh3,rh4,3,3,rh5)
      do i= 1,nen
        rh3=0.d0
        do j = 1,iMrow(i)
          do l = 1,iMcol(i)
            rh3(j,l) = shp(2,i)*MatM(j,l,i)+shp(3,i)*MatM_2(j,l,i)
          end do
        end do
        call matmulf(rh5,MatM(1,1,i),3,3,iMcol(i),rh4)
        call msmul(rh4,3,iMcol(i),shp(3,i))
        call matmulf(rh,rh3,3,3,iMcol(i),rh6)
        call matadd(rh4,rh6,3,iMcol(i),twn_2(1,1,i))
      end do         
      return
      end
      
      subroutine updGP63omInt (xl,ul,dwg,twn,dd,w,ix,ndf,ndm,nen,ityp,
     2                   ilin,k,shp,rd,r0,HutT_1,HutT_2,twn_1,twn_2,
     4                   iMcol,rinG,rinP,rwnGn,rwnGnp1,rwnPn,rwnPnp1,
     5                   t0exact,Iexactnormal,dwn)
      USE mdata
      USE iofile
c-----------------------------------------------------------------------
c     update of rotations in Gauss points
c     rinG    initial global nodal basis system
c     rinP    initial patch nodal basis system
c     rwnGnp1 current global nodal basis system of current iteration step
c     rwnGn   current global nodal basis system of last iteration step
c     rwnPnp1 current patch nodal basis system of current iteration step
c     rwnPn   current patch nodal basis system of last iteration step
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      include 'iofile.h'   
csk      common m(1)
      dimension  xl(ndm,*),ul(ndf,*),ix(*),
     1           dwg(9,*),twn(3,3,*),dd(3,3),
     2           w(3,3,*),rh(3,3),rh1(3,3),rh2(3,nen),
     4           shp(3,*),dbeta(3,3),
     5           rd(3,3),r0(3,3),rh3(3,3),Aini(3,3),
     6           drhd1(3,3),drhd2(3,3),rh4(3,3),iMrow(nen),iMcol(nen),
     7           rh5(3,3),rh6(3,3),twn_1(3,3,nen),twn_2(3,3,nen)
      dimension rinG(3,3,nen),rinP(3,3,nen),rwnPn(3,3,nen),
     1          rwnPnp1(3,3,nen),rwnGn(3,3,nen),rwnGnp1(3,3,nen)
      double precision MatH(3,3),MatH_1(3,3),MatH_2(3,3),
     1                 MatW(3,3),MatW_1(3,3),MatW_2(3,3),
     2                 MatT3(3,3),MatT3_1(3,3),MatT3_2(3,3),
     3                 MatM(3,3,nen),MatM_1(3,3,nen),MatM_2(3,3,nen),
     4                 MatT(3,3),MatWHT_1(3,3),MatWHT_2(3,3),
     5                 MatW_1HT(3,3),MatW_2HT(3,3),MatWH_1T(3,3),
     6                 MatWH_2T(3,3),MatT_1(3,3),MatT_2(3,3),
     7                 HutT_1(3,3,nen),HutT_2(3,3,nen),aid1(3,3),
     8                 aid2(3,3),bound,t0exact(3,3,3),dwn(3,*)
      real*8 rhtest(3,3)
c
      ismoothrow = 2 !in CMAME2014 ist hier 2 eingesetzt!
      bound=1.d-5
c.... initial cartesian system in Gauss point: rh1
c.... nodal rotations delta_beta_I (nodewise rotations): rh2
c      write(iow,*) dwg(1:9,k)
c     if loaded systems ar used
      if (Iexactnormal.eq.1) then
        r0    = t0exact(1:3,1:3,3)
        rh1   = t0exact(1:3,3,1:3)
        drhd1 = t0exact(1:3,1,1:3)
        drhd2 = t0exact(1:3,2,1:3)
        rd = 0.0d0
      else
        do i = 1,3
          do j = 1,3
            rh1(i,j) = 0.0d0
            rd(i,j)  = 0.d0
            r0(i,j)  = 0.d0
            drhd1(i,j) = 0.d0
            drhd2(i,j) = 0.d0
            do l = 1,nen
              rh1(i,j) = rh1(i,j) + rinP(i,j,l) * shp(3,l) 
              r0(i,j)  = r0 (i,j) + rinP(i,3,l) * shp(j,l)
              drhd1(i,j) = drhd1(i,j) +rinP(i,j,l) *shp(1,l)
              drhd2(i,j) = drhd2(i,j) +rinP(i,j,l) *shp(2,l)
            end do
          end do
        end do
      end if
c      if(ndf.eq.5) iflag = 1  
      do l= 1,nen
        if (ndf.eq.5) then
          iMrow(l) = ismoothrow !2
          iMcol(l) = 2
        else if (ndf.eq.6) then
csk       m(n7) => psid  
          call boun63(l,ix,psid,iflag,ndf)  ! iflag = 0 : Kante  
csk       call boun63(l,ix,m(n7),iflag,ndf) 
            if (iflag.eq.1) then
              iMrow(l) = ismoothrow !2
              iMcol(l) = 2 
            else if (iflag.eq.0) then
              iMrow(l) = 2  ! Hier ist im CMAME2014 2 eingesetzt, 3 funktioniert bei free form surface besser, 2 beim channel section beam!?!?
              iMcol(l) = 3
            end if

        else
          write(*,*) 'Wrong value for ndf: choose 5 or 6!!'
          write(iow,*) 'Wrong value for ndf: choose 5 or 6!!'
        end if
      end do

      
      rh2=0.d0
      do l = 1,nen
        do i = 1,iMcol(l)
          rh2(i,l) = ul(3+i,2*nen+l)
        end do
      end do
      call matcop  (rh1,3,3,Aini)

c      if(ndf.eq.5) iflag = 1                 ! iflag = 1 : glatte Schale
c                                            ! iflag = 0 : Kante  
      if(ilin.eq.2) then
c....   calculate R_n(->rh) from value omega_gp(<-dwg).
c....   Omega_gp is stored in the h3 field and loaded        
        call updr63 (dwg(1:3,k),rh,ilin,1)
        call matmulf (rh,rh1,3,3,3, rh3)
        call matcop  (rh3,3,3,rh1)
c...    1-direction
        call updRalpha63(dwg(1,k),rh3,ilin,1)
c...    R,alp*Ai
        call matmulf (rh3,Aini,3,3,3,rh5)
c...    R*R_0,1
        call matmulf (rh ,drhd1,3,3,3,rh4)
c...    R,alp*Ai + (R*R_0,1)*(R_0^T*Ai)
        call matadd  (rh5(1,1),rh4(1,1),3,1,aid1(1,1))
        call matadd  (rh5(1,2),rh4(1,2),3,1,aid1(1,2))
        call matadd  (rh5(1,3),rh4(1,3),3,1,aid1(1,3))
c...    2-direction
        call updRalpha63(dwg(1,k),rh3,ilin,2)
c...    R,alp*Ai
        call matmulf (rh3,Aini,3,3,3,rh5)
c...    R*R_0,2
        call matmulf (rh ,drhd2,3,3,3,rh4)
c...    R,alp*Ai + (R*R_0,1)*(R_0^T*Ai)
        call matadd  (rh5(1,1),rh4(1,1),3,1,aid2(1,1))
        call matadd  (rh5(1,2),rh4(1,2),3,1,aid2(1,2))
        call matadd  (rh5(1,3),rh4(1,3),3,1,aid2(1,3))
      else
c...    determine derivation of R_0(stored in rh1) for M_I,alpha
        call matcop(drhd1,3,3,aid1)
        call matcop(drhd2,3,3,aid2)
      end if
c....   determine M_I for all nodes based on rwn and rh1
      MatM=0.d0
      MatM_1=0.d0
      MatM_2=0.d0
      do l = 1,nen
c        if (iMcol(l).eq.3) then
          do i = 1,iMrow(l)
            do j = 1,iMcol(l)
              MatM(i,j,l)   = dot(rh1(1,i),rwnGn(1,j,l),3)
              MatM_1(i,j,l) = dot(aid1(1,i),rwnGn(1,j,l),3)
              MatM_2(i,j,l) = dot(aid2(1,i),rwnGn(1,j,l),3)
            end do
          end do
c        else
c        do i = 1,iMrow(l)
c            do j = 1,iMcol(l)
c        MatM(i,j,l) = dot(rh1(1,i),rwnGn(1,j,l),3)
c     -                   -( dot(rh1(1,3),rwnGn(1,j,l),3)
c     *                    * dot(rh1(1,i),rwnGn(1,3,l),3)
c     /                    / dot(rh1(1,3),rwnGn(1,3,l),3))
c              MatM_1(i,j,l) = dot(aid1(1,i),rwnGn(1,j,l),3)
c     -                    -((( dot(aid1(1,3),rwnGn(1,j,l),3)
c     *                        *dot(rh1(1,i),rwnGn(1,3,l),3)
c     +                        +dot(rh1(1,3),rwnGn(1,j,l),3)
c     *                        *dot(aid1(1,i),rwnGn(1,3,l),3))
c     *                        *dot(rh1(1,3),rwnGn(1,3,l),3)
c     -                      -( dot(aid1(1,3),rwnGn(1,3,l),3)
c     *                        *dot(rh1(1,3),rwnGn(1,j,l),3)
c     *                        *dot(rh1(1,i),rwnGn(1,3,l),3)))
c     /                     /(dot(rh1(1,3),rwnGn(1,3,l),3)**2))    
c              MatM_2(i,j,l) = dot(aid2(1,i),rwnGn(1,j,l),3)
c     -                    -((( dot(aid2(1,3),rwnGn(1,j,l),3)
c     *                        *dot(rh1(1,i),rwnGn(1,3,l),3)
c     +                        +dot(rh1(1,3),rwnGn(1,j,l),3)
c     *                        *dot(aid2(1,i),rwnGn(1,3,l),3))
c     *                        *dot(rh1(1,3),rwnGn(1,3,l),3)
c     -                      -( dot(aid2(1,3),rwnGn(1,3,l),3)
c     *                        *dot(rh1(1,3),rwnGn(1,j,l),3)
c     *                        *dot(rh1(1,i),rwnGn(1,3,l),3)))
c     /                     /(dot(rh1(1,3),rwnGn(1,3,l),3)**2))
c           end do
c        end do
c        end if
      end do

!c     test
!      do it = 1,3
!        do jt = 1,3
!          rhtest(it,jt) = 0.0d0
!          do l1 = 1,nen
!            rhtest(it,jt) = rhtest(it,jt) + rwnPn(it,jt,l1) * shp(3,l1) 
!          end do
!        end do
!      end do        
      
!c...  calculate delta_beta in Gauss point
!c...  transformation with A matrix necessary 
!c...  delta_beta = sum(N_I*A_I*delta_beta_I)
!      dbeta = 0.0d0
!      do l = 1,nen
!        do idev = 1,3
!          do i = 1,iMrow(l)
!            do j = 1,iMcol(l)
!c...          rh2(:,l) is delta_beta for node l   
!              dbeta(i,idev) = dbeta(i,idev) +
!     +                        shp(idev,l)*MatM(i,j,l)*rh2(j,l)
!              if (idev.eq.1) dbeta(i,idev) = dbeta(i,idev) +
!     +                        shp(3,l)*MatM_1(i,j,l)*rh2(j,l)
!              if (idev.eq.2) dbeta(i,idev) = dbeta(i,idev) +
!     +                        shp(3,l)*MatM_2(i,j,l)*rh2(j,l)
!            end do
!          end do
!        end do    
!      end do
!c...  calculate delta_omega in Gauss point
!c...  with rh1 being the current coordinate system in current GP
!c...  and the additional rotations delta_beta in current GP
!c...  the resulting additional rotation delta_omega i stored in rh2
!      do i = 1,3
!        call matmulf (rh1,dbeta(1,i),3,3,1, rh2(1,i))
!      end do
!      do i = 1,3
!        rh2(i,1) = rh2(i,1)+aid1(i,1)*dbeta(1,3)+aid1(i,2)*dbeta(2,3)
!     +             +aid1(i,3)*dbeta(3,3)
!        rh2(i,2) = rh2(i,2)+aid2(i,1)*dbeta(1,3)+aid2(i,2)*dbeta(2,3)
!     +             +aid2(i,3)*dbeta(3,3)
!      end do
!
!c.... update current axial vector Omega    
!      call matadd (rh2(1,3),dwg(1,k),3,1, dwg(1,k))  
!      call matadd (rh2(1,1),dwg(4,k),3,1, dwg(4,k))     
!      call matadd (rh2(1,2),dwg(7,k),3,1, dwg(7,k))             

c...  compute Omgea from nodal values
      rh2 = 0.d0
      do l = 1,nen
        do idev=1,3
          do i = 1,3
            rh2(i,idev) = rh2(i,idev) + shp(idev,l) * dwn(i,l)
          end do
        end do
      end do
c...  store to value for this Gauss point
      call matcop(rh2(1,3),3,1,dwg(1,k))
      call matcop(rh2(1,1),3,1,dwg(4,k))
      call matcop(rh2(1,2),3,1,dwg(7,k))
      
c.... update of current nodal basis
      if(ilin.eq.0) then                   ! linear, correct!!!

         call matcop (r0,3,3, rd)
         call vecp (dwg(1,k),r0(1,3),dd(1,3))
         
         call vecp (dwg(4,k),r0(1,3),dd(1,1))
         call vecp (dwg(1,k),r0(1,1),rh3)
         call matadd(dd(1,1),rh3,3,1,dd(1,1))
         
         call vecp (dwg(7,k),r0(1,3),dd(1,2))
         call vecp (dwg(1,k),r0(1,2),rh3)
         call matadd(dd(1,2),rh3,3,1,dd(1,2))

         call matcop(drhd1(1,1),3,3,MatT3_1(1,1))
         call matcop(drhd2(1,1),3,3,MatT3_2(1,1))


      else                                ! nonlinear
c....   rd has to be computed
        call updr63 (dwg(1:3,k),rh,ilin,1)
c....   d = R * D; ->rh maps initial to current tripod 
        call matmulf (rh,Aini(1,3),3,3,1,rd(1,3))
c...    for small rotations:
c....   d,alph = skew(omega,alpha) * D + delta_R * D,alpha
c....   for finite rotations:ß
c....   ai,alph=R,alp*Ai+R*R_0,alp*R_0^T*Ai
c....   rh3: R,alph
c....   Aini(1,i): Ai (D=A3) (=R_0)
c....   rh: R
c....   drhdalpha:R_0,alpha
c....   rh4: R,0
c
c
c....   1-direction
        call updRalpha63(dwg(1,k),rh3,ilin,1)
c...    R,alp*Ai
        call matmulf(rh3,Aini,3,3,3,rh5)
c...    R*R_0,1
        call matmulf (rh ,drhd1,3,3,3,rh4)
c...    R,alp*Ai + (R*R_0,1)*(R_0^T*Ai)
        call matadd  (rh5(1,3),rh4(1,3),3,1,rd(1,1))
        call matadd  (rh5(1,1),rh4(1,1),3,1,MatT3_1(1,1))
        call matadd  (rh5(1,2),rh4(1,2),3,1,MatT3_1(1,2))
        call matadd  (rh5(1,3),rh4(1,3),3,1,MatT3_1(1,3))
        call matadd  (rh5(1,1),rh4(1,1),3,1,aid1(1,1))
        call matadd  (rh5(1,2),rh4(1,2),3,1,aid1(1,2))
        call matadd  (rh5(1,3),rh4(1,3),3,1,aid1(1,3))

c....   2-direction
        call updRalpha63(dwg(1,k),rh3,ilin,2)
c...    R,alp*Ai
        call matmulf (rh3,Aini,3,3,3,rh5)
c....   R*R_0,2
        call matmulf (rh ,drhd2,3,3,3,rh4)
c...    R,alp*D + (R*R_0,1)*(R_0^T*D)
        call matadd  (rh5(1,3),rh4(1,3),3,1,rd(1,2))
        call matadd  (rh5(1,1),rh4(1,1),3,1,MatT3_2(1,1))
        call matadd  (rh5(1,2),rh4(1,2),3,1,MatT3_2(1,2))
        call matadd  (rh5(1,3),rh4(1,3),3,1,MatT3_2(1,3))
        call matadd  (rh5(1,1),rh4(1,1),3,1,aid2(1,1))
        call matadd  (rh5(1,2),rh4(1,2),3,1,aid2(1,2))
        call matadd  (rh5(1,3),rh4(1,3),3,1,aid2(1,3))

c....   Gauss point tripod has to be updated 
        if (ilin.eq.2)   rh1= matmul (rh,Aini)

        
c....   determine M_I for all nodes based on rwnn and current a_i
        do l = 1,nen
c          if (iMcol(l).eq.3) then
            do i = 1,iMrow(l)
              do j = 1,iMcol(l)
                MatM(i,j,l)   = dot(rh1(1,i),rwnGnp1(1,j,l),3)
                MatM_1(i,j,l) = dot(aid1(1,i),rwnGnp1(1,j,l),3)
                MatM_2(i,j,l) = dot(aid2(1,i),rwnGnp1(1,j,l),3)
              end do
            end do
c          else
c          do i = 1,iMrow(l)
c              do j = 1,iMcol(l)
c          MatM(i,j,l) = dot(rh1(1,i),rwnGnp1(1,j,l),3)
c     -                     -( dot(rh1(1,3),rwnGnp1(1,j,l),3)
c     *                      * dot(rh1(1,i),rwnGnp1(1,3,l),3)
c     /                      / dot(rh1(1,3),rwnGnp1(1,3,l),3))
c                MatM_1(i,j,l) = dot(aid1(1,i),rwnGnp1(1,j,l),3)
c     -                      -((( dot(aid1(1,3),rwnGnp1(1,j,l),3)
c     *                           *dot(rh1(1,i),rwnGnp1(1,3,l),3)
c     +                          +dot(rh1(1,3),rwnGnp1(1,j,l),3)
c     *                          *dot(aid1(1,i),rwnGnp1(1,3,l),3))
c     *                          *dot(rh1(1,3),rwnGnp1(1,3,l),3)
c     -                        -( dot(aid1(1,3),rwnGnp1(1,3,l),3)
c     *                          *dot(rh1(1,3),rwnGnp1(1,j,l),3)
c     *                          *dot(rh1(1,i),rwnGnp1(1,3,l),3)))
c     /                       /(dot(rh1(1,3),rwnGnp1(1,3,l),3)**2))    
c                MatM_2(i,j,l) = dot(aid2(1,i),rwnGnp1(1,j,l),3)
c     -                      -((( dot(aid2(1,3),rwnGnp1(1,j,l),3)
c     *                          *dot(rh1(1,i),rwnGnp1(1,3,l),3)
c     +                          +dot(rh1(1,3),rwnGnp1(1,j,l),3)
c     *                          *dot(aid2(1,i),rwnGnp1(1,3,l),3))
c     *                          *dot(rh1(1,3),rwnGnp1(1,3,l),3)
c     -                        -( dot(aid2(1,3),rwnGnp1(1,3,l),3)
c     *                          *dot(rh1(1,3),rwnGnp1(1,j,l),3)
c     *                          *dot(rh1(1,i),rwnGnp1(1,3,l),3)))
c     /                       /(dot(rh1(1,3),rwnGnp1(1,3,l),3)**2)) 
c          end do
c            end do
c          end if
        end do
      end if 
c
c...  T = W^T * H * T_3 and
c...  T_,alpha = W_,a^T*H*T3+W^T*H_,a*T3+W^T*H*T3_,a
c
c...  create T3, T3_,alpha is done above
      MatT3(1:3,1) = rh1(1:3,1)
      MatT3(1:3,2) = rh1(1:3,2)
      MatT3(1:3,3) = rh1(1:3,3)
      if (ilin.eq.1) then
        call matcop(aid1(1,1),3,3,MatT3_1(1,1))
        call matcop(aid2(1,1),3,3,MatT3_2(1,1))
      end if
c
c...  compute H and H_,alpha (for ilin.ne.2 H=1)
      call updr63 (dwg(1:3,k),MatH,ilin,3)    !HIER?????? CHECKEN ob dwg(1,k) richtigß???
      call updHalpha63(dwg(1,k),MatH_1,ilin,1)
      call updHalpha63(dwg(1,k),MatH_2,ilin,2)

c...  compute W and W_,alpha
      if (ilin.eq.2) then
        call skew(MatW,rd(1,3))
        call skew(MatW_1,rd(1,1))
        call skew(MatW_2,rd(1,2))
      else
        call skew(MatW,r0(1,3))
        call skew(MatW_1,r0(1,1))
        call skew(MatW_2,r0(1,2))
      end if
c
c...  calculate T = W^T * H * T_3
      call mttmul (MatW,MatH,3,3,3, rh3)
      call matmulf(rh3,MatT3,3,3,3,MatT)
c...  calculate TM_I = W^T * H *T_3 *M_I and store it to w(1,1,I)
      call pzero(w,9*nen)
      do i = 1,nen
c        call matmulf (MatT,MatM(1,1,i),3,iMrow(i),iMrow(i),w(1,1,i))
        call matmulf (MatT,MatM(1,1,i),3,3,3,w(1,1,i))
      end do
c
c...  calculate W^T * H * T_3,a
      call matmulf(rh3,MatT3_1,3,3,3,MatWHT_1)
      call matmulf(rh3,MatT3_2,3,3,3,MatWHT_2)
c...  calculate W_,a^T*H*T3
      call mttmul (MatW_1,MatH,3,3,3, rh3)
      call matmulf(rh3,MatT3,3,3,3, MatW_1HT)
      call mttmul (MatW_2,MatH,3,3,3, rh3)
      call matmulf(rh3,MatT3,3,3,3, MatW_2HT)
c...  calculate W^T*H_,a*T3
      call mttmul (MatW,MatH_1,3,3,3,rh3)
      call matmulf(rh3,MatT3,3,3,3,MatWH_1T)
      call mttmul (MatW,MatH_2,3,3,3,rh3)
      call matmulf(rh3,MatT3,3,3,3,MatWH_2T)
c
c...  sum up above terms for T_,alpha
      call matadd(MatWHT_1,MatW_1HT,3,3,rh3)
      call matadd(rh3,MatWH_1T,3,3,MatT_1)
      call matadd(MatWHT_2,MatW_2HT,3,3,rh3)
      call matadd(rh3,MatWH_2T,3,3,MatT_2)
c
c...  calculate ^T_I,alpha and store it to HutT_,alpha
      HutT_1=0.d0
      HutT_2=0.d0
      do i = 1,nen
c...    1-direction
        call matmulf(MatT_1,MatM(1,1,i),3,3,iMcol(i),rh3)
        call msmul(rh3,3,iMcol(i),shp(3,i))
        call matcop(w(1,1,i),3,iMcol(i),rh4)
        call msmul(rh4,3,iMcol(i),shp(1,i))
        call matmulf(MatT,MatM_1(1,1,i),3,3,iMcol(i),rh5)
        call msmul(rh5,3,iMcol(i),shp(3,i))
        call matadd(rh3,rh4,3,iMcol(i),rh6)
        call matadd(rh6,rh5,3,iMcol(i),HutT_1(1,1,i))
!        call matcop(rh4,3,2,HutT_1(1,1,i))
c...    2-direction
        call matmulf(MatT_2,MatM(1,1,i),3,3,iMcol(i),rh3)
        call msmul(rh3,3,iMcol(i),shp(3,i))
        call matcop(w(1,1,i),3,iMcol(i),rh4)
        call msmul(rh4,3,iMcol(i),shp(2,i))
        call matmulf(MatT,MatM_2(1,1,i),3,3,iMcol(i),rh5)
        call msmul(rh5,3,iMcol(i),shp(3,i))
        call matadd(rh3,rh4,3,iMcol(i),rh6)
        call matadd(rh6,rh5,3,iMcol(i),HutT_2(1,1,i))
!        call matcop(rh4,3,2,HutT_2(1,1,i))
      end do
c
c...  for G_ßß part of geometrical matrix
c...  compute T_HTM_I = H * T_3 * M_I and store it to twn(1,1,I) 
      call matmulf(MatH,MatT3,3,3,3,rh)
      call pzero(twn,9*nen)
      call pzero(twn_1,9*nen)
      call pzero(twn_2,9*nen)
      do i = 1,nen
        call matmulf (rh,MatM(1,1,i),3,3,iMcol(i),twn(1,1,i))
      end do
c...  compute T_HTM_I,alpha and store it to twn_alpha(1,1,I)
c... 1-direction
      call matmulf(MatH_1,MatT3,3,3,3,rh3)
      call matmulf(MatH,MatT3_1,3,3,3,rh4)
      call matadd (rh3,rh4,3,3,rh5)
      do i= 1,nen
        rh3=0.d0
        do j = 1,iMrow(i)
          do l = 1,iMcol(i)
            rh3(j,l) = shp(1,i)*MatM(j,l,i)+shp(3,i)*MatM_1(j,l,i)
          end do
        end do
        call matmulf(rh5,MatM(1,1,i),3,3,iMcol(i),rh4)
        call msmul(rh4,3,iMcol(i),shp(3,i))
        call matmulf(rh,rh3,3,3,iMcol(i),rh6)
        call matadd(rh4,rh6,3,iMcol(i),twn_1(1,1,i))
      end do
c... 2-direction
      call matmulf(MatH_2,MatT3,3,3,3,rh3)
      call matmulf(MatH,MatT3_2,3,3,3,rh4)
      call matadd (rh3,rh4,3,3,rh5)
      do i= 1,nen
        rh3=0.d0
        do j = 1,iMrow(i)
          do l = 1,iMcol(i)
            rh3(j,l) = shp(2,i)*MatM(j,l,i)+shp(3,i)*MatM_2(j,l,i)
          end do
        end do
        call matmulf(rh5,MatM(1,1,i),3,3,iMcol(i),rh4)
        call msmul(rh4,3,iMcol(i),shp(3,i))
        call matmulf(rh,rh3,3,3,iMcol(i),rh6)
        call matadd(rh4,rh6,3,iMcol(i),twn_2(1,1,i))
      end do         
      return
      end
      

      subroutine updGP63Rupd2 (xl,ul,dwg,twn,dd,w,ix,ndf,ndm,nen,ityp,
     2                   ilin,k,shp,rd,r0,HutT_1,HutT_2,twn_1,twn_2,
     4                   iMcol,rinG,rinP,rwnGn,rwnGnp1,rwnPn,rwnPnp1,
     5                   t0exact,Iexactnormal,Riold)
      USE mdata
      USE iofile
c-----------------------------------------------------------------------
c     update of rotations in Gauss points
c     rinG    initial global nodal basis system
c     rinP    initial patch nodal basis system
c     rwnGnp1 current global nodal basis system of current iteration step
c     rwnGn   current global nodal basis system of last iteration step
c     rwnPnp1 current patch nodal basis system of current iteration step
c     rwnPn   current patch nodal basis system of last iteration step
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      include 'iofile.h'   
csk      common m(1)
      dimension  xl(ndm,*),ul(ndf,*),ix(*),
     1           dwg(9),twn(3,3,*),dd(3,3),
     2           w(3,3,*),rh(3,3),rh1(3,3),rh2(3,nen),
     4           shp(3,*),dbeta(3,3),
     5           rd(3,3),r0(3,3),rh3(3,3),Aini(3,3),
     6           drhd1(3,3),drhd2(3,3),rh4(3,3),iMrow(nen),iMcol(nen),
     7           rh5(3,3),rh6(3,3),twn_1(3,3,nen),twn_2(3,3,nen)
      dimension rinG(3,3,nen),rinP(3,3,nen),rwnPn(3,3,nen),
     1          rwnPnp1(3,3,nen),rwnGn(3,3,nen),rwnGnp1(3,3,nen)
      dimension Riold(3,3,3,*)
      double precision MatH(3,3),MatH_1(3,3),MatH_2(3,3),
     1                 MatW(3,3),MatW_1(3,3),MatW_2(3,3),
     2                 MatT3(3,3),MatT3_1(3,3),MatT3_2(3,3),
     3                 MatM(3,3,nen),MatM_1(3,3,nen),MatM_2(3,3,nen),
     4                 MatT(3,3),MatWHT_1(3,3),MatWHT_2(3,3),
     5                 MatW_1HT(3,3),MatW_2HT(3,3),MatWH_1T(3,3),
     6                 MatWH_2T(3,3),MatT_1(3,3),MatT_2(3,3),
     7                 HutT_1(3,3,nen),HutT_2(3,3,nen),aid1(3,3),
     8                 aid2(3,3),bound,t0exact(3,3,3)
      real*8 rhtest(3,3)
c
      ismoothrow = 3 !in CMAME2014 ist hier 2 eingesetzt!
      bound=1.d-5
c.... initial cartesian system in Gauss point: rh1
c.... nodal rotations delta_beta_I (nodewise rotations): rh2
c      write(iow,*) dwg(1:9,k)
c     if loaded systems ar used
      if (Iexactnormal.eq.1) then
        r0    = t0exact(1:3,1:3,3)
        rh1   = t0exact(1:3,3,1:3)
        drhd1 = t0exact(1:3,1,1:3)
        drhd2 = t0exact(1:3,2,1:3)
        rd = 0.0d0
      else
        do i = 1,3
          do j = 1,3
            rh1(i,j) = 0.0d0
            rd(i,j)  = 0.d0
            r0(i,j)  = 0.d0
            drhd1(i,j) = 0.d0
            drhd2(i,j) = 0.d0
            do l = 1,nen
              rh1(i,j) = rh1(i,j) + rinP(i,j,l) * shp(3,l) 
              r0(i,j)  = r0 (i,j) + rinP(i,3,l) * shp(j,l)
              drhd1(i,j) = drhd1(i,j) +rinP(i,j,l) *shp(1,l)
              drhd2(i,j) = drhd2(i,j) +rinP(i,j,l) *shp(2,l)
            end do
          end do
        end do
      end if
c      if(ndf.eq.5) iflag = 1  
      do l= 1,nen
        if (ndf.eq.5) then
          iMrow(l) = ismoothrow !2
          iMcol(l) = 2
        else if (ndf.eq.6) then
csk       m(n7) => psid  
          call boun63(l,ix,psid,iflag,ndf)  ! iflag = 0 : Kante  
csk       call boun63(l,ix,m(n7),iflag,ndf) 
            if (iflag.eq.1) then
              iMrow(l) = ismoothrow !2
              iMcol(l) = 2 
            else if (iflag.eq.0) then
              iMrow(l) = 3  ! nur bei 3 ist free form surface with kink lastschrittunabhängig!!! Hier ist im CMAME2014 2 eingesetzt
              iMcol(l) = 3
            end if

        else
          write(*,*) 'Wrong value for ndf: choose 5 or 6!!'
          write(iow,*) 'Wrong value for ndf: choose 5 or 6!!'
        end if
      end do

      
      rh2=0.d0
      do l = 1,nen
        do i = 1,iMcol(l)
          rh2(i,l) = ul(3+i,2*nen+l)
        end do
      end do
      call matcop  (rh1,3,3,Aini)

c      if(ndf.eq.5) iflag = 1                 ! iflag = 1 : glatte Schale
c                                            ! iflag = 0 : Kante  
      if(ilin.eq.2) then
c....   load R_n-1 from h3 field and multiply with initial basis system
        call matmulf (Riold(1:3,1:3,3,k),Aini,3,3,3,rh3)
        call matcop  (rh3,3,3,rh1)
c...    1-direction
c...    R,alp*Ai        
        call matmulf (Riold(1:3,1:3,1,k),Aini,3,3,3,rh3)
c...    R*Ai,alp
        call matmulf (Riold(1:3,1:3,3,k),drhd1,3,3,3,rh4)
        call matadd (rh3,rh4,3,3,aid1)
c...    2-direction
c...    R,alp*Ai        
        call matmulf (Riold(1:3,1:3,2,k),Aini,3,3,3,rh3)
c...    R*Ai,alp
        call matmulf (Riold(1:3,1:3,3,k),drhd2,3,3,3,rh4)
        call matadd (rh3,rh4,3,3,aid2)        
      else
c...    determine derivation of R_0(stored in rh1) for M_I,alpha
        call matcop(drhd1,3,3,aid1)
        call matcop(drhd2,3,3,aid2)
      end if
c....   determine M_I for all nodes based on rwn and rh1
      MatM=0.d0
      MatM_1=0.d0
      MatM_2=0.d0
      do l = 1,nen
c        if (iMcol(l).eq.3) then
          do i = 1,iMrow(l)
            do j = 1,iMcol(l)
              MatM(i,j,l)   = dot(rh1(1,i),rwnGn(1,j,l),3)
              MatM_1(i,j,l) = dot(aid1(1,i),rwnGn(1,j,l),3)
              MatM_2(i,j,l) = dot(aid2(1,i),rwnGn(1,j,l),3)
            end do
          end do
c        else
c        do i = 1,iMrow(l)
c            do j = 1,iMcol(l)
c        MatM(i,j,l) = dot(rh1(1,i),rwnGn(1,j,l),3)
c     -                   -( dot(rh1(1,3),rwnGn(1,j,l),3)
c     *                    * dot(rh1(1,i),rwnGn(1,3,l),3)
c     /                    / dot(rh1(1,3),rwnGn(1,3,l),3))
c              MatM_1(i,j,l) = dot(aid1(1,i),rwnGn(1,j,l),3)
c     -                    -((( dot(aid1(1,3),rwnGn(1,j,l),3)
c     *                        *dot(rh1(1,i),rwnGn(1,3,l),3)
c     +                        +dot(rh1(1,3),rwnGn(1,j,l),3)
c     *                        *dot(aid1(1,i),rwnGn(1,3,l),3))
c     *                        *dot(rh1(1,3),rwnGn(1,3,l),3)
c     -                      -( dot(aid1(1,3),rwnGn(1,3,l),3)
c     *                        *dot(rh1(1,3),rwnGn(1,j,l),3)
c     *                        *dot(rh1(1,i),rwnGn(1,3,l),3)))
c     /                     /(dot(rh1(1,3),rwnGn(1,3,l),3)**2))    
c              MatM_2(i,j,l) = dot(aid2(1,i),rwnGn(1,j,l),3)
c     -                    -((( dot(aid2(1,3),rwnGn(1,j,l),3)
c     *                        *dot(rh1(1,i),rwnGn(1,3,l),3)
c     +                        +dot(rh1(1,3),rwnGn(1,j,l),3)
c     *                        *dot(aid2(1,i),rwnGn(1,3,l),3))
c     *                        *dot(rh1(1,3),rwnGn(1,3,l),3)
c     -                      -( dot(aid2(1,3),rwnGn(1,3,l),3)
c     *                        *dot(rh1(1,3),rwnGn(1,j,l),3)
c     *                        *dot(rh1(1,i),rwnGn(1,3,l),3)))
c     /                     /(dot(rh1(1,3),rwnGn(1,3,l),3)**2))
c           end do
c        end do
c        end if
      end do

c     test
      do it = 1,3
        do jt = 1,3
          rhtest(it,jt) = 0.0d0
          do l1 = 1,nen
            rhtest(it,jt) = rhtest(it,jt) + rwnPn(it,jt,l1) * shp(3,l1) 
          end do
        end do
      end do        
      
c...  calculate delta_beta in Gauss point
c...  transformation with A matrix necessary 
c...  delta_beta = sum(N_I*A_I*delta_beta_I)
      dbeta = 0.0d0
      do l = 1,nen
        do idev = 1,3
          do i = 1,iMrow(l)
            do j = 1,iMcol(l)
c...          rh2(:,l) is delta_beta for node l   
              dbeta(i,idev) = dbeta(i,idev) +
     +                        shp(idev,l)*MatM(i,j,l)*rh2(j,l)
              if (idev.eq.1) dbeta(i,idev) = dbeta(i,idev) +
     +                        shp(3,l)*MatM_1(i,j,l)*rh2(j,l)
              if (idev.eq.2) dbeta(i,idev) = dbeta(i,idev) +
     +                        shp(3,l)*MatM_2(i,j,l)*rh2(j,l)
            end do
          end do
        end do    
      end do
c...  calculate delta_omega in Gauss point
c...  with rh1 being the current coordinate system in current GP
c...  and the additional rotations delta_beta in current GP
c...  the resulting additional rotation delta_omega i stored in rh2
      do i = 1,3
        call matmulf (rh1,dbeta(1,i),3,3,1, rh2(1,i))
      end do
      do i = 1,3
        rh2(i,1) = rh2(i,1)+aid1(i,1)*dbeta(1,3)+aid1(i,2)*dbeta(2,3)
     +             +aid1(i,3)*dbeta(3,3)
        rh2(i,2) = rh2(i,2)+aid2(i,1)*dbeta(1,3)+aid2(i,2)*dbeta(2,3)
     +             +aid2(i,3)*dbeta(3,3)
      end do

c.... The current axial vector Omega is now stored in rh2(1:3,3)
c.... Its derivatives in rh2(1:3,1) and rh2(1:3,2)
c.... Store this for further use as omega_n in dwg
      call matcop (rh2(1,3),3,1,dwg(1))
      call matcop (rh2(1,1),3,1,dwg(4))
      call matcop (rh2(1,2),3,1,dwg(7))
      

          

c.... update of current interpolated basis
      if(ilin.eq.0) then                   ! linear, HERE A FIELD TO STORE DWG OF LAST STEP IS MISSING; USE 1,2 INSTEAD FOR LINEAR!!!
        write(*,*) 'use additive version (command 1,2,0,0)'
         !call matcop (r0,3,3, rd)
         !call vecp (dwg(1),r0(1,3),dd(1,3))
         !
         !call vecp (dwg(4),r0(1,3),dd(1,1))
         !call vecp (dwg(1),r0(1,1),rh3)
         !call matadd(dd(1,1),rh3,3,1,dd(1,1))
         !
         !call vecp (dwg(7),r0(1,3),dd(1,2))
         !call vecp (dwg(1),r0(1,2),rh3)
         !call matadd(dd(1,2),rh3,3,1,dd(1,2))
         !
         !call matcop(drhd1(1,1),3,3,MatT3_1(1,1))
         !call matcop(drhd2(1,1),3,3,MatT3_2(1,1))


      else                                ! nonlinear
c....   now the current rotation matrix R_n has to be computed    
c....   R is computed to rh
        call updr63 (dwg(1:3),rh,ilin,1)
c....   update of Riold: R_n = R * R_n-1: stored in rh3
        call matmulf (rh,Riold(1:3,1:3,3,k),3,3,3,rh3)
c....   rd has to be computed: d = R_n * D
        call matmulf (rh3,Aini(1,3),3,3,1,rd(1,3))
        
c....   1-direction
        call updRalpha63(dwg(1),rh4,ilin,1)
c....   compute R_n,alp = R,alp * R_(n-1) + R * R_(n-1),alp
        call matmulf (rh4,Riold(1:3,1:3,3,k),3,3,3,rh5)
        call matmulf (rh,Riold(1:3,1:3,1,k),3,3,3,rh4)
        call matadd (rh4,rh5,3,3,Riold(1:3,1:3,1,k))
c....   2-direction
        call updRalpha63(dwg(1),rh4,ilin,2)
c....   compute R_n,alp = R,alp * R_(n-1) + R * R_(n-1),alp
        call matmulf (rh4,Riold(1:3,1:3,3,k),3,3,3,rh5)
        call matmulf (rh,Riold(1:3,1:3,2,k),3,3,3,rh4)
        call matadd (rh4,rh5,3,3,Riold(1:3,1:3,2,k))        
c....   now store the new value R_n for next iteration, R_n,alp is already stored    
        call matcop (rh3,3,3,Riold(1:3,1:3,3,k))
c
c....   compute derived director and basis vectors
c....   1-direction
c...    R_n,alp*Ai
        call matmulf(Riold(1:3,1:3,1,k),Aini,3,3,3,rh5)
c...    R_n*Ai,alp
        call matmulf (rh3 ,drhd1,3,3,3,rh4)
c...    R_n,alp*Ai + R_n*Ai,alp
        call matadd  (rh5(1,3),rh4(1,3),3,1,rd(1,1))
        call matadd  (rh5(1,1),rh4(1,1),3,1,MatT3_1(1,1))
        call matadd  (rh5(1,2),rh4(1,2),3,1,MatT3_1(1,2))
        call matadd  (rh5(1,3),rh4(1,3),3,1,MatT3_1(1,3))
        call matadd  (rh5(1,1),rh4(1,1),3,1,aid1(1,1))
        call matadd  (rh5(1,2),rh4(1,2),3,1,aid1(1,2))
        call matadd  (rh5(1,3),rh4(1,3),3,1,aid1(1,3))

c....   2-direction
c...    R_n,alp*Ai
        call matmulf (Riold(1:3,1:3,2,k),Aini,3,3,3,rh5)
c....   R_n*Ai,alp
        call matmulf (rh3 ,drhd2,3,3,3,rh4)
c...    R_n,alp*Ai + R_n*Ai,alp
        call matadd  (rh5(1,3),rh4(1,3),3,1,rd(1,2))
        call matadd  (rh5(1,1),rh4(1,1),3,1,MatT3_2(1,1))
        call matadd  (rh5(1,2),rh4(1,2),3,1,MatT3_2(1,2))
        call matadd  (rh5(1,3),rh4(1,3),3,1,MatT3_2(1,3))
        call matadd  (rh5(1,1),rh4(1,1),3,1,aid2(1,1))
        call matadd  (rh5(1,2),rh4(1,2),3,1,aid2(1,2))
        call matadd  (rh5(1,3),rh4(1,3),3,1,aid2(1,3))

c....   Gauss point tripod has to be updated 
        if (ilin.eq.2)   rh1= matmul (rh3,Aini)

        
c....   determine M_I for all nodes based on rwnn and current a_i
        do l = 1,nen
c          if (iMcol(l).eq.3) then
            do i = 1,iMrow(l)
              do j = 1,iMcol(l)
                MatM(i,j,l)   = dot(rh1(1,i),rwnGnp1(1,j,l),3)
                MatM_1(i,j,l) = dot(aid1(1,i),rwnGnp1(1,j,l),3)
                MatM_2(i,j,l) = dot(aid2(1,i),rwnGnp1(1,j,l),3)
              end do
            end do
c          else
c          do i = 1,iMrow(l)
c              do j = 1,iMcol(l)
c          MatM(i,j,l) = dot(rh1(1,i),rwnGnp1(1,j,l),3)
c     -                     -( dot(rh1(1,3),rwnGnp1(1,j,l),3)
c     *                      * dot(rh1(1,i),rwnGnp1(1,3,l),3)
c     /                      / dot(rh1(1,3),rwnGnp1(1,3,l),3))
c                MatM_1(i,j,l) = dot(aid1(1,i),rwnGnp1(1,j,l),3)
c     -                      -((( dot(aid1(1,3),rwnGnp1(1,j,l),3)
c     *                           *dot(rh1(1,i),rwnGnp1(1,3,l),3)
c     +                          +dot(rh1(1,3),rwnGnp1(1,j,l),3)
c     *                          *dot(aid1(1,i),rwnGnp1(1,3,l),3))
c     *                          *dot(rh1(1,3),rwnGnp1(1,3,l),3)
c     -                        -( dot(aid1(1,3),rwnGnp1(1,3,l),3)
c     *                          *dot(rh1(1,3),rwnGnp1(1,j,l),3)
c     *                          *dot(rh1(1,i),rwnGnp1(1,3,l),3)))
c     /                       /(dot(rh1(1,3),rwnGnp1(1,3,l),3)**2))    
c                MatM_2(i,j,l) = dot(aid2(1,i),rwnGnp1(1,j,l),3)
c     -                      -((( dot(aid2(1,3),rwnGnp1(1,j,l),3)
c     *                          *dot(rh1(1,i),rwnGnp1(1,3,l),3)
c     +                          +dot(rh1(1,3),rwnGnp1(1,j,l),3)
c     *                          *dot(aid2(1,i),rwnGnp1(1,3,l),3))
c     *                          *dot(rh1(1,3),rwnGnp1(1,3,l),3)
c     -                        -( dot(aid2(1,3),rwnGnp1(1,3,l),3)
c     *                          *dot(rh1(1,3),rwnGnp1(1,j,l),3)
c     *                          *dot(rh1(1,i),rwnGnp1(1,3,l),3)))
c     /                       /(dot(rh1(1,3),rwnGnp1(1,3,l),3)**2)) 
c          end do
c            end do
c          end if
        end do
      end if 
c
c...  T = W^T * H * T_3 and
c...  T_,alpha = W_,a^T*H*T3+W^T*H_,a*T3+W^T*H*T3_,a
c
c...  create T3, T3_,alpha is done above
      MatT3(1:3,1) = rh1(1:3,1)
      MatT3(1:3,2) = rh1(1:3,2)
      MatT3(1:3,3) = rh1(1:3,3)
      if (ilin.eq.1) then
        call matcop(aid1(1,1),3,3,MatT3_1(1,1))
        call matcop(aid2(1,1),3,3,MatT3_2(1,1))
      end if
c
c...  compute H and H_,alpha (for ilin.ne.2 H=1)
      !call updr63 (dwg(1:3),MatH,ilin,3)    !HIER?????? CHECKEN ob dwg(1,k) richtigß???
      !call updHalpha63(dwg(1),MatH_1,ilin,1)
      !call updHalpha63(dwg(1),MatH_2,ilin,2)
      ilintest = 0
c...  compute H and H_,alpha (for ilin.ne.2 H=1)
      call updr63 (dwg(1:3),MatH,ilintest,3)    !HIER?????? CHECKEN ob dwg(1,k) richtigß???
      call updHalpha63(dwg(1),MatH_1,ilintest,1)
      call updHalpha63(dwg(1),MatH_2,ilintest,2)
      MatH_1=0.d0
      MatH_2=0.d0

c...  compute W and W_,alpha
      if (ilin.eq.2) then
        call skew(MatW,rd(1,3))
        call skew(MatW_1,rd(1,1))
        call skew(MatW_2,rd(1,2))
      else
        call skew(MatW,r0(1,3))
        call skew(MatW_1,r0(1,1))
        call skew(MatW_2,r0(1,2))
      end if
c
c...  calculate T = W^T * H * T_3
      call mttmul (MatW,MatH,3,3,3, rh3)
      call matmulf(rh3,MatT3,3,3,3,MatT)
c...  calculate TM_I = W^T * H *T_3 *M_I and store it to w(1,1,I)
      call pzero(w,9*nen)
      do i = 1,nen
c        call matmulf (MatT,MatM(1,1,i),3,iMrow(i),iMrow(i),w(1,1,i))
        call matmulf (MatT,MatM(1,1,i),3,3,3,w(1,1,i))
      end do
c
c...  calculate W^T * H * T_3,a
      call matmulf(rh3,MatT3_1,3,3,3,MatWHT_1)
      call matmulf(rh3,MatT3_2,3,3,3,MatWHT_2)
c...  calculate W_,a^T*H*T3
      call mttmul (MatW_1,MatH,3,3,3, rh3)
      call matmulf(rh3,MatT3,3,3,3, MatW_1HT)
      call mttmul (MatW_2,MatH,3,3,3, rh3)
      call matmulf(rh3,MatT3,3,3,3, MatW_2HT)
c...  calculate W^T*H_,a*T3
      call mttmul (MatW,MatH_1,3,3,3,rh3)
      call matmulf(rh3,MatT3,3,3,3,MatWH_1T)
      call mttmul (MatW,MatH_2,3,3,3,rh3)
      call matmulf(rh3,MatT3,3,3,3,MatWH_2T)
c
c...  sum up above terms for T_,alpha
      call matadd(MatWHT_1,MatW_1HT,3,3,rh3)
      call matadd(rh3,MatWH_1T,3,3,MatT_1)
      call matadd(MatWHT_2,MatW_2HT,3,3,rh3)
      call matadd(rh3,MatWH_2T,3,3,MatT_2)
c
c...  calculate ^T_I,alpha and store it to HutT_,alpha
      HutT_1=0.d0
      HutT_2=0.d0
      do i = 1,nen
c...    1-direction
        call matmulf(MatT_1,MatM(1,1,i),3,3,iMcol(i),rh3)
        call msmul(rh3,3,iMcol(i),shp(3,i))
        call matcop(w(1,1,i),3,iMcol(i),rh4)
        call msmul(rh4,3,iMcol(i),shp(1,i))
        call matmulf(MatT,MatM_1(1,1,i),3,3,iMcol(i),rh5)
        call msmul(rh5,3,iMcol(i),shp(3,i))
        call matadd(rh3,rh4,3,iMcol(i),rh6)
        call matadd(rh6,rh5,3,iMcol(i),HutT_1(1,1,i))
!        call matcop(rh4,3,2,HutT_1(1,1,i))
c...    2-direction
        call matmulf(MatT_2,MatM(1,1,i),3,3,iMcol(i),rh3)
        call msmul(rh3,3,iMcol(i),shp(3,i))
        call matcop(w(1,1,i),3,iMcol(i),rh4)
        call msmul(rh4,3,iMcol(i),shp(2,i))
        call matmulf(MatT,MatM_2(1,1,i),3,3,iMcol(i),rh5)
        call msmul(rh5,3,iMcol(i),shp(3,i))
        call matadd(rh3,rh4,3,iMcol(i),rh6)
        call matadd(rh6,rh5,3,iMcol(i),HutT_2(1,1,i))
!        call matcop(rh4,3,2,HutT_2(1,1,i))
      end do
c
c...  for G_ßß part of geometrical matrix
c...  compute T_HTM_I = H * T_3 * M_I and store it to twn(1,1,I) 
      call matmulf(MatH,MatT3,3,3,3,rh)
      call pzero(twn,9*nen)
      call pzero(twn_1,9*nen)
      call pzero(twn_2,9*nen)
      do i = 1,nen
        call matmulf (rh,MatM(1,1,i),3,3,iMcol(i),twn(1,1,i))
      end do
c...  compute T_HTM_I,alpha and store it to twn_alpha(1,1,I)
c... 1-direction
      call matmulf(MatH_1,MatT3,3,3,3,rh3)
      call matmulf(MatH,MatT3_1,3,3,3,rh4)
      call matadd (rh3,rh4,3,3,rh5)
      do i= 1,nen
        rh3=0.d0
        do j = 1,iMrow(i)
          do l = 1,iMcol(i)
            rh3(j,l) = shp(1,i)*MatM(j,l,i)+shp(3,i)*MatM_1(j,l,i)
          end do
        end do
        call matmulf(rh5,MatM(1,1,i),3,3,iMcol(i),rh4)
        call msmul(rh4,3,iMcol(i),shp(3,i))
        call matmulf(rh,rh3,3,3,iMcol(i),rh6)
        call matadd(rh4,rh6,3,iMcol(i),twn_1(1,1,i))
      end do
c... 2-direction
      call matmulf(MatH_2,MatT3,3,3,3,rh3)
      call matmulf(MatH,MatT3_2,3,3,3,rh4)
      call matadd (rh3,rh4,3,3,rh5)
      do i= 1,nen
        rh3=0.d0
        do j = 1,iMrow(i)
          do l = 1,iMcol(i)
            rh3(j,l) = shp(2,i)*MatM(j,l,i)+shp(3,i)*MatM_2(j,l,i)
          end do
        end do
        call matmulf(rh5,MatM(1,1,i),3,3,iMcol(i),rh4)
        call msmul(rh4,3,iMcol(i),shp(3,i))
        call matmulf(rh,rh3,3,3,iMcol(i),rh6)
        call matadd(rh4,rh6,3,iMcol(i),twn_2(1,1,i))
      end do         
      return
      end
      
            subroutine updGP63omegaAdd (xl,ul,dwg,twn,dd,w,ix,ndf,ndm,
     1                   nen,ityp,
     2                   ilin,k,shp,rd,r0,HutT_1,HutT_2,twn_1,twn_2,
     4                   iMcol,rinG,rinP,rwnGnp1,
     5                   t0exact,Iexactnormal,dwn)
      USE mdata
      USE iofile
c-----------------------------------------------------------------------
c     update of rotations in Gauss points
c     rinG    initial global nodal basis system
c     rinP    initial patch nodal basis system
c     rwnGnp1 current global nodal basis system of current iteration step
c     rwnGn   current global nodal basis system of last iteration step
c     rwnPnp1 current patch nodal basis system of current iteration step
c     rwnPn   current patch nodal basis system of last iteration step
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      include 'iofile.h'   
csk      common m(1)
      dimension  xl(ndm,*),ul(ndf,*),ix(*),
     1           dwg(9),twn(3,3,*),dd(3,3),
     2           w(3,3,*),rh(3,3),rh1(3,3),rh2(3,nen),
     4           shp(3,*),dbeta(3,3),domega(3,3),
     5           rd(3,3),r0(3,3),rh3(3,3),Aini(3,3),
     6           drhd1(3,3),drhd2(3,3),rh4(3,3),iMrow(nen),iMcol(nen),
     7           rh5(3,3),rh6(3,3),twn_1(3,3,nen),twn_2(3,3,nen)
      dimension rinG(3,3,nen),rinP(3,3,nen),rwnGnp1(3,3,nen)
      double precision MatH(3,3),MatH_1(3,3),MatH_2(3,3),
     1                 MatW(3,3),MatW_1(3,3),MatW_2(3,3),
     4                 MatT_1(3,3),MatT_2(3,3),
     7                 HutT_1(3,3,nen),HutT_2(3,3,nen),
     8                 t0exact(3,3,3),dwn(3,*)
      double precision MatT3(3,3,nen)
 
c

c.... initial cartesian system in Gauss point: rh1
c.... nodal rotations delta_beta_I (nodewise rotations): rh2
c      write(iow,*) dwg(1:9,k)
c     if loaded systems are used
      if (Iexactnormal.eq.1) then
        r0    = t0exact(1:3,1:3,3)
        rh1   = t0exact(1:3,3,1:3)
        drhd1 = t0exact(1:3,1,1:3)
        drhd2 = t0exact(1:3,2,1:3)
        rd = 0.0d0
      else
        do i = 1,3
          do j = 1,3
            rh1(i,j) = 0.0d0
            rd(i,j)  = 0.d0
            r0(i,j)  = 0.d0
            drhd1(i,j) = 0.d0
            drhd2(i,j) = 0.d0
            do l = 1,nen
              rh1(i,j) = rh1(i,j) + rinP(i,j,l) * shp(3,l) 
              r0(i,j)  = r0 (i,j) + rinP(i,3,l) * shp(j,l)
              drhd1(i,j) = drhd1(i,j) +rinP(i,j,l) *shp(1,l)
              drhd2(i,j) = drhd2(i,j) +rinP(i,j,l) *shp(2,l)
            end do
          end do
        end do
      end if
      call matcop  (rh1,3,3,Aini)

      do l= 1,nen
        if (ndf.eq.5) then
          iMrow(l) = 3
          iMcol(l) = 2
        else if (ndf.eq.6) then
csk       m(n7) => psid  
          call boun63(l,ix,psid,iflag,ndf)  ! iflag = 0 : Kante  
csk       call boun63(l,ix,m(n7),iflag,ndf) 
            if (iflag.eq.1) then
              iMrow(l) = 3
              iMcol(l) = 2
            else if (iflag.eq.0) then
              iMrow(l) = 3  
              iMcol(l) = 3
            end if
        else
          write(*,*) 'Wrong value for ndf: choose 5 or 6!!'
          write(iow,*) 'Wrong value for ndf: choose 5 or 6!!'
        end if
      end do

c...  compute Omgea from total nodal values
      rh2 = 0.d0
      do l = 1,nen
        do idev=1,3
          do i = 1,3
            rh2(i,idev) = rh2(i,idev) + shp(idev,l) * dwn(i,l)
          end do
        end do
      end do
c...  write as value for this Gauss point
      call matcop(rh2(1,3),3,1,dwg(1))
      call matcop(rh2(1,1),3,1,dwg(4))
      call matcop(rh2(1,2),3,1,dwg(7))   
      
      
c.... update of current nodal basis
      if(ilin.eq.0) then                   ! linear, correct!!!

        call matcop (r0,3,3, rd)
        call vecp (dwg(1),r0(1,3),dd(1,3))
        
        call vecp (dwg(4),r0(1,3),dd(1,1))
        call vecp (dwg(1),r0(1,1),rh3)
        call matadd(dd(1,1),rh3,3,1,dd(1,1))
         
        call vecp (dwg(7),r0(1,3),dd(1,2))
        call vecp (dwg(1),r0(1,2),rh3)
        call matadd(dd(1,2),rh3,3,1,dd(1,2))
c...    write matrix T_3          
        do l = 1,nen
          if (iMcol(l).eq.2) then
             call matcop (rinG(1,1,l),3,2,MatT3(1,1,l))
             call pzero(MatT3(1,3,l),3)
          end if
          if (iMcol(l).eq.3) then
            call matcop (rinG(1,1,l),3,3,MatT3(1,1,l))
          end if             
        end do

      else                                ! nonlinear
c....   rd has to be computed
        call updr63 (dwg(1:3),rh,ilin,1)
c....   d = R * D; ->rh maps initial to current tripod 
        call matmulf (rh,Aini(1,3),3,3,1,rd(1,3))
c...    for small rotations:
c....   d,alph = skew(omega,alpha) * D + delta_R * D,alpha
c....   for finite rotations:ß
c....   ai,alph=R,alp*Ai+R*R_0,alp*R_0^T*Ai
c....   rh3: R,alph
c....   Aini(1,i): Ai (D=A3) (=R_0)
c....   rh: R
c....   drhdalpha:R_0,alpha
c....   rh4: R,0
c
c
c....   1-direction
        call updRalpha63(dwg(1),rh3,ilin,1)
c...    R,alp*Ai
        call matmulf(rh3,Aini,3,3,3,rh5)
c...    R*R_0,1
        call matmulf (rh ,drhd1,3,3,3,rh4)
c...    R,alp*Ai + (R*R_0,1)*(R_0^T*Ai)
        call matadd  (rh5(1,3),rh4(1,3),3,1,rd(1,1))

c....   2-direction
        call updRalpha63(dwg(1),rh3,ilin,2)
c...    R,alp*Ai
        call matmulf (rh3,Aini,3,3,3,rh5)
c....   R*R_0,2
        call matmulf (rh ,drhd2,3,3,3,rh4)
c...    R,alp*D + (R*R_0,1)*(R_0^T*D)
        call matadd  (rh5(1,3),rh4(1,3),3,1,rd(1,2))

c....   Gauss point tripod has to be updated 
        if (ilin.eq.2)   rh1= matmul (rh,Aini)

c...    write matrix T_3        
        do l = 1,nen
          if (iMcol(l).eq.2) then
             call matcop (rwnGnp1(1,1,l),3,2,MatT3(1,1,l))
             call pzero(MatT3(1,3,l),3)
          end if
          if (iMcol(l).eq.3) then
            call matcop (rwnGnp1(1,1,l),3,3,MatT3(1,1,l))
          end if             
        end do
      end if 
c
c...  T = W^T * H  and
c...  T_,alpha = W_,a^T*H + W^T*H_,a

c...  compute H and H_,alpha (for ilin.ne.2 H=1)
      call updr63 (dwg(1:3),MatH,ilin,3)    
      call updHalpha63(dwg(1),MatH_1,ilin,1)
      call updHalpha63(dwg(1),MatH_2,ilin,2)

c...  compute W and W_,alpha
      if (ilin.eq.2) then
        call skew(MatW,rd(1,3))
        call skew(MatW_1,rd(1,1))
        call skew(MatW_2,rd(1,2))
      else
        call skew(MatW,r0(1,3))
        call skew(MatW_1,r0(1,1))
        call skew(MatW_2,r0(1,2))
      end if
c
c...  calculate T = W^T * H * T_3
      call mttmul (MatW,MatH,3,3,3, rh3)
c...  calculate T_I = W^T * H *T_3I and store it to w(1,1,I)
      call pzero(w,9*nen)
      do i = 1,nen
        call matmulf (rh3,MatT3(1,1,i),3,3,iMcol(i),w(1,1,i))
      end do
c
c...  calculate W_,a^T*H
      call mttmul (MatW_1,MatH,3,3,3, rh3)
      call mttmul (MatW_2,MatH,3,3,3, rh5)
c...  calculate W^T*H_,a
      call mttmul (MatW,MatH_1,3,3,3,rh4)
      call mttmul (MatW,MatH_2,3,3,3,rh6)
c
c...  sum up above terms for T_,alpha
      call matadd(rh3,rh4,3,3,MatT_1)
      call matadd(rh5,rh6,3,3,MatT_2)
c
c...  calculate ^T_I,alpha and store it to HutT_,alpha
      HutT_1=0.d0
      HutT_2=0.d0
      do i = 1,nen
c...    1-direction
        call matmulf(MatT_1,MatT3(1,1,i),3,3,iMcol(i),rh3)
        call msmul(rh3,3,iMcol(i),shp(3,i))
        call matcop(w(1,1,i),3,iMcol(i),rh4)
        call msmul(rh4,3,iMcol(i),shp(1,i))
        call matadd (rh3,rh4,3,iMcol(i),HutT_1(1,1,i))
c...    2-direction
        call matmulf(MatT_2,MatT3(1,1,i),3,3,iMcol(i),rh3)
        call msmul(rh3,3,iMcol(i),shp(3,i))
        call matcop(w(1,1,i),3,iMcol(i),rh4)
        call msmul(rh4,3,iMcol(i),shp(2,i))
        call matadd(rh3,rh4,3,iMcol(i),HutT_2(1,1,i))
      end do
c
c...  for G_ßß part of geometrical matrix
c...  compute T_HT_I = H * T_3I and store it to twn(1,1,I) 

      call pzero(twn,9*nen)
      call pzero(twn_1,9*nen)
      call pzero(twn_2,9*nen)
      do i = 1,nen
          call matmulf (MatH,MatT3(1,1,i),3,3,iMcol(i),twn(1,1,i))
      end do
c...  compute T_HT_I,alpha and store it to twn_alpha(1,1,I)
c... 1-direction
      do i= 1,nen
        call matmulf (MatH_1,MatT3(1,1,i),3,3,iMcol(i),rh3)
        call msmul (rh3,3,iMcol(i),shp(3,i))
        call matmulf (MatH,MatT3(1,1,i),3,3,iMcol(i),rh4)
        call msmul (rh4,3,iMcol(i),shp(1,i))
        call matadd(rh3,rh4,3,iMcol(i),twn_1(1,1,i))
      end do
c... 2-direction
      do i= 1,nen
        call matmulf (MatH_2,MatT3(1,1,i),3,3,iMcol(i),rh3)
        call msmul (rh3,3,iMcol(i),shp(3,i))
        call matmulf (MatH,MatT3(1,1,i),3,3,iMcol(i),rh4)
        call msmul (rh4,3,iMcol(i),shp(2,i))
        call matadd(rh3,rh4,3,iMcol(i),twn_2(1,1,i))
      end do
        
      return
      end      
c
      subroutine updGP63omegaMult (xl,ul,dwg,twn,dd,w,ix,ndf,ndm,nen,
     2                   ityp,ilin,k,shp,rd,r0,HutT_1,HutT_2,
     4                   iMcol,rinP,t0exact,
     5                   Iexactnormal,ddwn,Riold)
      USE mdata
      USE iofile
c-----------------------------------------------------------------------
c     update of rotations in Gauss points
c     rinG    initial global nodal basis system
c     rinP    initial patch nodal basis system
c     rwnGnp1 current global nodal basis system of current iteration step
c     rwnGn   current global nodal basis system of last iteration step
c     rwnPnp1 current patch nodal basis system of current iteration step
c     rwnPn   current patch nodal basis system of last iteration step
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      include 'iofile.h'   
csk      common m(1)
      dimension  xl(ndm,*),ul(ndf,*),ix(*),
     1           dwg(9),twn(3,3,*),dd(3,3),
     2           w(3,3,*),rh(3,3),rh1(3,3),rh2(3,nen),
     4           shp(3,*),
     5           rd(3,3),r0(3,3),rh3(3,3),Aini(3,3),
     6           drhd1(3,3),drhd2(3,3),rh4(3,3),iMcol(nen),
     7           rh5(3,3),rh6(3,3)
      dimension rinP(3,3,nen)
      double precision MatW(3,3),MatW_1(3,3),MatW_2(3,3),
     7                 HutT_1(3,3,nen),HutT_2(3,3,nen),t0exact(3,3,3)
      dimension Riold(3,3,3,*)
      real*8 ddwn(3,*)
c

c.... initial cartesian system in Gauss point: rh1
c.... nodal rotations delta_beta_I (nodewise rotations): rh2
c      write(iow,*) dwg(1:9,k)
c     if loaded systems ar used
      if (Iexactnormal.eq.1) then
        r0    = t0exact(1:3,1:3,3)
        rh1   = t0exact(1:3,3,1:3)
        drhd1 = t0exact(1:3,1,1:3)
        drhd2 = t0exact(1:3,2,1:3)
        rd = 0.0d0
      else
        do i = 1,3
          do j = 1,3
            rh1(i,j) = 0.0d0
            rd(i,j)  = 0.d0
            r0(i,j)  = 0.d0
            drhd1(i,j) = 0.d0
            drhd2(i,j) = 0.d0
            do l = 1,nen
              rh1(i,j) = rh1(i,j) + rinP(i,j,l) * shp(3,l) 
              r0(i,j)  = r0 (i,j) + rinP(i,3,l) * shp(j,l)
              drhd1(i,j) = drhd1(i,j) +rinP(i,j,l) *shp(1,l)
              drhd2(i,j) = drhd2(i,j) +rinP(i,j,l) *shp(2,l)
            end do
          end do
        end do
      end if
      call matcop  (rh1,3,3,Aini)
c      if(ndf.eq.5) iflag = 1  
      do l= 1,nen
        if (ndf.eq.5) then
          iMcol(l) = 2
        else if (ndf.eq.6) then
csk       m(n7) => psid  
          call boun63(l,ix,psid,iflag,ndf)  ! iflag = 0 : Kante  
csk       call boun63(l,ix,m(n7),iflag,ndf) 
            if (iflag.eq.1) then
              iMcol(l) = 2
            else if (iflag.eq.0) then
              iMcol(l) = 3
            end if
        else
          write(*,*) 'Wrong value for ndf: choose 5 or 6!!'
          write(iow,*) 'Wrong value for ndf: choose 5 or 6!!'
        end if
      end do
      
c...  compute delta Omega from nodal values
      rh2 = 0.d0
      do l = 1,nen
        do idev=1,3
          do i = 1,3
            rh2(i,idev) = rh2(i,idev) + shp(idev,l) * ddwn(i,l)
          end do
        end do
      end do

c...  store to value for this Gauss point (resorting required for function updRalpha63)
c...  delta omega values
      call matcop(rh2(1,3),3,1,dwg(1))
      call matcop(rh2(1,1),3,1,dwg(4))
      call matcop(rh2(1,2),3,1,dwg(7))  

      
c.... update of current nodal basis
      if(ilin.eq.0) then                   ! linear
         write(*,*) 'use additive version (command 1,3,0,0)'
      else                                ! nonlinear
c...    perform multiplicative update in integration point
c...    compute R and deltaR from nodal values, then perform update
c...    R_n = R * R_(n-1)     
c....   R is computed to rh     
        call updr63 (dwg(1:3),rh,ilin,1)
c....   multiplicative update: R_n = R * R_n-1: stored in rh3
        call matmulf (rh,Riold(1:3,1:3,3,k),3,3,3,rh3)
    
c....   rd has to be computed
c....   d = R_n * D; ->rh maps initial to current tripod 
        call matmulf (rh3,Aini(1,3),3,3,1,rd(1,3))    

c...    Computation of derivatives of R_n:        
c...    R_n,alp = R,alp * R_(n-1) + R * R_(n-1),alp    

c....   1-direction

c....   R,alp is stored to rh4
        call updRalpha63(dwg(1),rh4,ilin,1)
c....   compute R_n,alp = R,alp * R_(n-1) + R * R_(n-1),alp
        call matmulf (rh4,Riold(1:3,1:3,3,k),3,3,3,rh5)
        call matmulf (rh,Riold(1:3,1:3,1,k),3,3,3,rh4)
        call matadd (rh4,rh5,3,3,Riold(1:3,1:3,1,k)) 
        
c....   2-direction
        call updRalpha63(dwg(1),rh4,ilin,2)
c....   compute R_n,alp = R,alp * R_(n-1) + R * R_(n-1),alp
        call matmulf (rh4,Riold(1:3,1:3,3,k),3,3,3,rh5)
        call matmulf (rh,Riold(1:3,1:3,2,k),3,3,3,rh4)
        call matadd (rh4,rh5,3,3,Riold(1:3,1:3,2,k))        
c....   now store the new value R_n for next iteration, R_n,alp is already stored    
        call matcop (rh3,3,3,Riold(1:3,1:3,3,k))

c....   compute derived director
c....   1-direction
c....   d,alp = R_n,alp * D + R_n * D,alp        
c...    R_n,alp*Ai
        call matmulf(Riold(1:3,1:3,1,k),Aini,3,3,3,rh5)
c...    R_n*Ai,alp
        call matmulf (rh3 ,drhd1,3,3,3,rh4)
c...    R_n,alp*Ai + R_n*Ai,alp
        call matadd  (rh4(1,3),rh5(1,3),3,1,rd(1,1))
  
        
c....   2-direction
c....   d,alp = R_n,alp * D + R_n * D,alp        
c...    R_n,alp*Ai
        call matmulf (Riold(1:3,1:3,2,k),Aini,3,3,3,rh5)
c...    R_n*Ai,alp
        call matmulf (rh3 ,drhd2,3,3,3,rh4)
c...    R_n,alp*Ai + R_n*Ai,alp
        call matadd  (rh4(1,3),rh5(1,3),3,1,rd(1,2))        

c...    for small rotations:
c....   d,alph = skew(omega,alpha) * D + delta_R * D,alpha

!c....   Gauss point tripod has to be updated 
!        if (ilin.eq.2)   rh1= matmul (rh3,Aini)

!c...    write matrix T_3        
!        do l = 1,nen
!          if (iMcol(l).eq.2) then
!             call matcop (rwnGnp1(1,1,l),3,2,MatT3(1,1,l))
!             call pzero(MatT3(1,3,l),3)
!          end if
!          if (iMcol(l).eq.3) then
!            call matcop (rwnGnp1(1,1,l),3,3,MatT3(1,1,l))
!          end if             
!        end do
      end if 
          
      
c...  compute W and W_,alpha
      if (ilin.eq.2) then
        call skew(MatW,rd(1,3))
        call skew(MatW_1,rd(1,1))
        call skew(MatW_2,rd(1,2))
      else
        call skew(MatW,r0(1,3))
        call skew(MatW_1,r0(1,1))
        call skew(MatW_2,r0(1,2))
      end if
c
c...  calculate T_I = W^T * T_3I and store it to w(1,1,I)
      call pzero(w,9*nen)
      do i = 1,nen
        call mttmul (MatW,twn(1,1,i),3,3,iMcol(i),w(1,1,i))
      end do
      
c...  calculate ^T_I,alpha and store it to HutT_,alpha
      HutT_1=0.d0
      HutT_2=0.d0
      do i = 1,nen
c...    1-direction
        call mttmul(MatW_1,twn(1,1,i),3,3,iMcol(i),rh3)
        call msmul(rh3,3,iMcol(i),shp(3,i))
        call matcop(w(1,1,i),3,iMcol(i),rh4)
        call msmul(rh4,3,iMcol(i),shp(1,i))
        call matadd (rh3,rh4,3,iMcol(i),HutT_1(1,1,i))
c...    2-direction
        call mttmul(MatW_2,twn(1,1,i),3,3,iMcol(i),rh3)
        call msmul(rh3,3,iMcol(i),shp(3,i))
        call matcop(w(1,1,i),3,iMcol(i),rh4)
        call msmul(rh4,3,iMcol(i),shp(2,i))
        call matadd(rh3,rh4,3,iMcol(i),HutT_2(1,1,i))
      end do


c
!c...  for G_ßß part of geometrical matrix
!c...  compute T_HT_I = T_3I and store it to twn(1,1,I) 
!c...  compute T_HT_I,alpha and store it to twn_alpha(1,1,I)
!      call pzero(twn,9*nen)
!      call pzero(twn_1,9*nen)
!      call pzero(twn_2,9*nen)
!      do i = 1,nen
!        call matcop (MatT3(1,1,i),3,iMcol(i),twn(1,1,i))
!c... 1-direction          
!        call matcop (MatT3(1,1,i),3,iMcol(i),twn_1(1,1,i))
!        call msmul (twn_1(1,1,i),3,iMcol(i),shp(1,i))
!c... 2-direction
!        call matcop (MatT3(1,1,i),3,iMcol(i),twn_2(1,1,i))
!        call msmul (twn_2(1,1,i),3,iMcol(i),shp(2,i))
!      end do
        
      return
      end      
c
      
      subroutine updGP63omega (xl,ul,dwg,twn,dd,w,ix,ndf,ndm,nen,ityp,
     2                   ilin,k,shp,rd,r0,HutT_1,HutT_2,twn_1,twn_2,
     4                   iMcol,rinG,rinP,rwnGn,rwnGnp1,rwnPn,rwnPnp1,
     5                   t0exact,Iexactnormal,dwn,dwnold,ddwn)
      USE mdata
      USE iofile
c-----------------------------------------------------------------------
c     update of rotations in Gauss points
c     rinG    initial global nodal basis system
c     rinP    initial patch nodal basis system
c     rwnGnp1 current global nodal basis system of current iteration step
c     rwnGn   current global nodal basis system of last iteration step
c     rwnPnp1 current patch nodal basis system of current iteration step
c     rwnPn   current patch nodal basis system of last iteration step
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      include 'iofile.h'   
csk      common m(1)
      dimension  xl(ndm,*),ul(ndf,*),ix(*),
     1           dwg(9,*),twn(3,3,*),dd(3,3),
     2           w(3,3,*),rh(3,3),rh1(3,3),rh2(3,nen),
     4           shp(3,*),dbeta(3,3),domega(3,3),
     5           rd(3,3),r0(3,3),rh3(3,3),Aini(3,3),
     6           drhd1(3,3),drhd2(3,3),rh4(3,3),iMrow(nen),iMcol(nen),
     7           rh5(3,3),rh6(3,3),twn_1(3,3,nen),twn_2(3,3,nen)
      dimension rinG(3,3,nen),rinP(3,3,nen),rwnPn(3,3,nen),
     1          rwnPnp1(3,3,nen),rwnGn(3,3,nen),rwnGnp1(3,3,nen)
      double precision MatH(3,3),MatH_1(3,3),MatH_2(3,3),
     1                 MatW(3,3),MatW_1(3,3),MatW_2(3,3),
     2                 MatT3(3,3),MatT3_1(3,3),MatT3_2(3,3),
     2                 MatT3_K(3,3),MatT3_K_12(3,3),MatT_K(3,3),
     3                 MatT_K_1(3,3),MatT_K_2(3,3), 
     3                 MatM(3,3,nen),MatM_1(3,3,nen),MatM_2(3,3,nen),
     4                 MatT(3,3),MatWHT_1(3,3),MatWHT_2(3,3),
     5                 MatW_1HT(3,3),MatW_2HT(3,3),MatWH_1T(3,3),
     6                 MatWH_2T(3,3),MatT_1(3,3),MatT_2(3,3),
     7                 HutT_1(3,3,nen),HutT_2(3,3,nen),aid1(3,3),
     8                 aid2(3,3),bound,t0exact(3,3,3),dwn(3,*)
      real*8 rhtest(3,3),dwnold(3,*),ddwn(3,*),dwgold(9),rh7(3,3),
     1       rh8(3,3)
c
      ismoothrow = 3
      iKinkpresent = 0
      bound=1.d-5
c.... initial cartesian system in Gauss point: rh1
c.... nodal rotations delta_beta_I (nodewise rotations): rh2
c      write(iow,*) dwg(1:9,k)
c     if loaded systems ar used
      if (Iexactnormal.eq.1) then
        r0    = t0exact(1:3,1:3,3)
        rh1   = t0exact(1:3,3,1:3)
        drhd1 = t0exact(1:3,1,1:3)
        drhd2 = t0exact(1:3,2,1:3)
        rd = 0.0d0
      else
        do i = 1,3
          do j = 1,3
            rh1(i,j) = 0.0d0
            rd(i,j)  = 0.d0
            r0(i,j)  = 0.d0
            drhd1(i,j) = 0.d0
            drhd2(i,j) = 0.d0
            do l = 1,nen
              rh1(i,j) = rh1(i,j) + rinP(i,j,l) * shp(3,l) 
              r0(i,j)  = r0 (i,j) + rinP(i,3,l) * shp(j,l)
              drhd1(i,j) = drhd1(i,j) +rinP(i,j,l) *shp(1,l)
              drhd2(i,j) = drhd2(i,j) +rinP(i,j,l) *shp(2,l)
            end do
          end do
        end do
      end if
      call matcop  (rh1,3,3,Aini)
c      if(ndf.eq.5) iflag = 1  
      do l= 1,nen
        if (ndf.eq.5) then
          iMrow(l) = ismoothrow !2
          iMcol(l) = 2
        else if (ndf.eq.6) then
csk       m(n7) => psid  
          call boun63(l,ix,psid,iflag,ndf)  ! iflag = 0 : Kante  
csk       call boun63(l,ix,m(n7),iflag,ndf) 
            if (iflag.eq.1) then
              iMrow(l) = ismoothrow !2
              iMcol(l) = 2
            else if (iflag.eq.0) then
              iMrow(l) = 3  
              iMcol(l) = 3
              iKinkpresent = 1
            end if

        else
          write(*,*) 'Wrong value for ndf: choose 5 or 6!!'
          write(iow,*) 'Wrong value for ndf: choose 5 or 6!!'
        end if
      end do
      
!c     If node on kink is present initialize matrices
!      if (iKinkpresent.eq.1) then
!        MatT3_K = 0.0d0
!        MatT3_K(1,1) = 1.0d0
!        MatT3_K(2,2) = 1.0d0
!        MatT3_K(3,3) = 1.0d0
!        MatT3_K_12 = 0.0d0
!      end if
!        
!      rh2=0.d0
!      do l = 1,nen
!        do i = 1,iMcol(l)
!          rh2(i,l) = ul(3+i,2*nen+l)
!        end do
!      end do
!      call matcop  (rh1,3,3,Aini)

!c      if(ndf.eq.5) iflag = 1                 ! iflag = 1 : glatte Schale
!c                                            ! iflag = 0 : Kante  
!      if(ilin.eq.2) then
!c....   calculate R_n(->rh) from value omega_gp(<-dwg).
!c....   Omega_gp is stored in the h3 field and loaded        
!        call updr63 (dwg(1:3,k),rh,ilin,1)
!        call matmulf (rh,rh1,3,3,3, rh3)
!        call matcop  (rh3,3,3,rh1)
!c...    1-direction
!        call updRalpha63(dwg(1,k),rh3,ilin,1)
!c...    R,alp*Ai
!        call matmulf (rh3,Aini,3,3,3,rh5)
!c...    R*R_0,1
!        call matmulf (rh ,drhd1,3,3,3,rh4)
!c...    R,alp*Ai + (R*R_0,1)*(R_0^T*Ai)
!        call matadd  (rh5(1,1),rh4(1,1),3,1,aid1(1,1))
!        call matadd  (rh5(1,2),rh4(1,2),3,1,aid1(1,2))
!        call matadd  (rh5(1,3),rh4(1,3),3,1,aid1(1,3))
!c...    2-direction
!        call updRalpha63(dwg(1,k),rh3,ilin,2)
!c...    R,alp*Ai
!        call matmulf (rh3,Aini,3,3,3,rh5)
!c...    R*R_0,2
!        call matmulf (rh ,drhd2,3,3,3,rh4)
!c...    R,alp*Ai + (R*R_0,1)*(R_0^T*Ai)
!        call matadd  (rh5(1,1),rh4(1,1),3,1,aid2(1,1))
!        call matadd  (rh5(1,2),rh4(1,2),3,1,aid2(1,2))
!        call matadd  (rh5(1,3),rh4(1,3),3,1,aid2(1,3))
!      else
!c...    determine derivation of R_0(stored in rh1) for M_I,alpha
!        call matcop(drhd1,3,3,aid1)
!        call matcop(drhd2,3,3,aid2)
!      end if
!c....   determine M_I for all nodes based on rwn and rh1
!      MatM=0.d0
!      MatM_1=0.d0
!      MatM_2=0.d0
!      do l = 1,nen
!        if (iMrow(l).eq.ismoothrow.and.iMcol(l).eq.2) then
!          do i = 1,iMrow(l)
!            do j = 1,iMcol(l)
!              MatM(i,j,l)   = dot(rh1(1,i),rwnGn(1,j,l),3)
!              MatM_1(i,j,l) = dot(aid1(1,i),rwnGn(1,j,l),3)
!              MatM_2(i,j,l) = dot(aid2(1,i),rwnGn(1,j,l),3)
!            end do
!          end do
!        else
!          do i = 1,iMrow(l)
!            do j = 1,iMcol(l)
!              if (i.eq.j) then
!                MatM(i,j,l)   = 1.0d0
!                MatM_1(i,j,l) = 0.0d0
!                MatM_2(i,j,l) = 0.0d0
!              else
!                MatM(i,j,l)   = 0.0d0
!                MatM_1(i,j,l) = 0.0d0
!                MatM_2(i,j,l) = 0.0d0
!              end if
!            end do
!          end do
!        end if
!      end do
!
!c     test
!      do it = 1,3
!        do jt = 1,3
!          rhtest(it,jt) = 0.0d0
!          do l1 = 1,nen
!            rhtest(it,jt) = rhtest(it,jt) + rwnPn(it,jt,l1) * shp(3,l1) 
!          end do
!        end do
!      end do        
      
!c...  calculate delta_beta in Gauss point
!c...  transformation with A matrix necessary 
!c...  delta_beta = sum(N_I*A_I*delta_beta_I)
!      dbeta = 0.0d0
!      domega = 0.0d0
!      do l = 1,nen
!        if (iMrow(l).eq.ismoothrow.and.iMcol(l).eq.2) then
!          do idev = 1,3
!            do i = 1,iMrow(l)
!              do j = 1,iMcol(l)
!c...            rh2(:,l) is delta_beta for node l   
!                dbeta(i,idev) = dbeta(i,idev) +
!     +                        shp(idev,l)*MatM(i,j,l)*rh2(j,l)
!                if (idev.eq.1) dbeta(i,idev) = dbeta(i,idev) +
!     +                        shp(3,l)*MatM_1(i,j,l)*rh2(j,l)
!                if (idev.eq.2) dbeta(i,idev) = dbeta(i,idev) +
!     +                        shp(3,l)*MatM_2(i,j,l)*rh2(j,l)
!              end do
!            end do
!          end do
!        elseif (iMrow(l).eq.3.and.iMcol(l).eq.3) then
!          do i=1,3  
!            domega(i,1) = domega(i,1) + shp(1,l)*rh2(i,l)
!            domega(i,2) = domega(i,2) + shp(2,l)*rh2(i,l)
!            domega(i,3) = domega(i,3) + shp(3,l)*rh2(i,l)
!          end do
!        end if
!      end do
!c...  calculate delta_omega in Gauss point
!c...  with rh1 being the current coordinate system in current GP
!c...  and the additional rotations delta_beta in current GP
!c...  the resulting additional rotation delta_omega i stored in rh2
!      do i = 1,3
!        call matmulf (rh1,dbeta(1,i),3,3,1, rh2(1,i))
!      end do
!      do i = 1,3
!        rh2(i,1) = rh2(i,1)+aid1(i,1)*dbeta(1,3)+aid1(i,2)*dbeta(2,3)
!     +             +aid1(i,3)*dbeta(3,3)
!        rh2(i,2) = rh2(i,2)+aid2(i,1)*dbeta(1,3)+aid2(i,2)*dbeta(2,3)
!     +             +aid2(i,3)*dbeta(3,3)
!      end do
!c.... add rotations due from omega to rotations due from beta
!      call matadd (domega,rh2,3,3,rh2)
!      
!c.... update current axial vector Omega    
!      call matadd (rh2(1,3),dwg(1,k),3,1, dwg(1,k))  
!      call matadd (rh2(1,1),dwg(4,k),3,1, dwg(4,k))     
!      call matadd (rh2(1,2),dwg(7,k),3,1, dwg(7,k))             

!c...  compute Omgea from nodal values
!      rh2 = 0.d0
!      do l = 1,nen
!        do idev=1,3
!          do i = 1,3
!            rh2(i,idev) = rh2(i,idev) + shp(idev,l) * dwn(i,l)
!          end do
!        end do
!      end do
c...  compute Omgea from nodal values, 
c...  separately for omega_n-1 and delta omega
c...  then perform multiplicative update (or additive?) in integration point
c...  compute R and deltaR from nodal values, then perform update!
c...  HIER WEITER MACHEN, DAS EINBAUEN
c...  BEI MULTI UPDATE WIRD R_,alpha so bestimmt!!!
c...  siehe line 4555 updgp63Rupd2
c...  R_n,alp = R,alp * R_(n-1) + R * R_(n-1),alp
c...  MÜSSTE BEI irupd=2 und iripl=2 implementiert sein!
c...  UND DANN H EHER MAL WEGLASSEN
      rh2 = 0.d0
      rh3 = 0.d0
      do l = 1,nen
        do idev=1,3
          do i = 1,3
            rh2(i,idev) = rh2(i,idev) + shp(idev,l) * ddwn(i,l)
            rh3(i,idev) = rh3(i,idev) + shp(idev,l) * dwnold(i,l)
          end do
        end do
      end do

c...  store to value for this Gauss point (resorting required for function updRalpha63)
c...  delta omega values
      call matcop(rh2(1,3),3,1,dwg(1,k))
      call matcop(rh2(1,1),3,1,dwg(4,k))
      call matcop(rh2(1,2),3,1,dwg(7,k))  
c...  omega_(n-1) values      
      call matcop(rh3(1,3),3,1,dwgold(1))
      call matcop(rh3(1,1),3,1,dwgold(4))
      call matcop(rh3(1,2),3,1,dwgold(7))  
      
c.... update of current nodal basis
      if(ilin.eq.0) then                   ! linear
c...    additive update of rotation vector
        call matadd (dwg(1:9,k),dwgold(1:9),9,9,dwgold(1:9))

        call matcop (r0,3,3, rd)
        call vecp (dwgold(1),r0(1,3),dd(1,3))
        
        call vecp (dwgold(4),r0(1,3),dd(1,1))
        call vecp (dwgold(1),r0(1,1),rh3)
        call matadd(dd(1,1),rh3,3,1,dd(1,1))
         
        call vecp (dwgold(7),r0(1,3),dd(1,2))
        call vecp (dwgold(1),r0(1,2),rh3)
        call matadd(dd(1,2),rh3,3,1,dd(1,2))
        MatM=0.d0
        MatM_1=0.d0
        MatM_2=0.d0
        do l = 1,nen
          do i = 1,iMrow(l)
            do j = 1,iMcol(l)
              MatM(i,j,l)   = rinG(i,j,l)
            end do
          end do
        end do

      else                                ! nonlinear
c...    perform multiplicative update in integration point
c...    compute R and deltaR from nodal values, then perform update
c...    R_n = R * R_(n-1)     
c....   R is computed to rh3      
        call updr63 (dwg(1:3,k),rh3,ilin,1)
c....   R_(n-1) is computed to rh4
        call updr63 (dwgold(1:3),rh4,ilin,1)
c....   multiplicative update: R_n = R * R_n-1: stored in rh
        call matmulf (rh3,rh4,3,3,3,rh)       
c....   rd has to be computed
c....   d = R_n * D; ->rh maps initial to current tripod 
        call matmulf (rh,Aini(1,3),3,3,1,rd(1,3))    

c...    Computation of derivatives of R_n:        
c...    R_n,alp = R,alp * R_(n-1) + R * R_(n-1),alp    

c....   1-direction

c....   R,alp is stored to rh5
        call updRalpha63(dwg(1,k),rh5,ilin,1)
c....   R_(n-1),alp is stored to rh6  
        call updRalpha63(dwgold(1),rh6,ilin,1)
c....   compute R_n,alp = R,alp * R_(n-1) + R * R_(n-1),alp and store to rh6   
        call matmulf (rh5,rh4,3,3,3,rh7)
        call matmulf (rh3,rh6,3,3,3,rh8)
        call matadd (rh7,rh8,3,3,rh6)
c....   compute derived director
c....   d,alp = R_n,alp * D + R_n * D,alp        
c...    R_n,alp*Ai
        call matmulf(rh6,Aini,3,3,3,rh7)
c...    R_n*Ai,alp
        call matmulf (rh ,drhd1,3,3,3,rh8)
c...    R_n,alp*Ai + R_n*Ai,alp
        call matadd  (rh7(1,3),rh8(1,3),3,1,rd(1,1))
  
        
c....   2-direction

c....   R,alp is stored to rh5
        call updRalpha63(dwg(1,k),rh5,ilin,2)
c....   R_(n-1),alp is stored to rh6  
        call updRalpha63(dwgold(1),rh6,ilin,2)
c....   compute R_n,alp = R,alp * R_(n-1) + R * R_(n-1),alp and store to rh6   
        call matmulf (rh5,rh4,3,3,3,rh7)
        call matmulf (rh3,rh6,3,3,3,rh8)
        call matadd (rh7,rh8,3,3,rh6)
c....   compute derived director
c....   d,alp = R_n,alp * D + R_n * D,alp        
c...    R_n,alp*Ai
        call matmulf(rh6,Aini,3,3,3,rh7)
c...    R_n*Ai,alp
        call matmulf (rh ,drhd2,3,3,3,rh8)
c...    R_n,alp*Ai + R_n*Ai,alp
        call matadd  (rh7(1,3),rh8(1,3),3,1,rd(1,2))        

c...    for small rotations:
c....   d,alph = skew(omega,alpha) * D + delta_R * D,alpha
!c
!c
!c....   1-direction
!        call updRalpha63(dwg(1,k),rh3,ilin,1)
!c...    R,alp*Ai
!        call matmulf(rh3,Aini,3,3,3,rh5)
!c...    R*R_0,1
!        call matmulf (rh ,drhd1,3,3,3,rh4)
!c...    R,alp*Ai + (R*R_0,1)*(R_0^T*Ai)
!        call matadd  (rh5(1,3),rh4(1,3),3,1,rd(1,1))
!        call matadd  (rh5(1,1),rh4(1,1),3,1,aid1(1,1))
!        call matadd  (rh5(1,2),rh4(1,2),3,1,aid1(1,2))
!        call matadd  (rh5(1,3),rh4(1,3),3,1,aid1(1,3))
!
!c....   2-direction
!        call updRalpha63(dwg(1,k),rh3,ilin,2)
!c...    R,alp*Ai
!        call matmulf (rh3,Aini,3,3,3,rh5)
!c....   R*R_0,2
!        call matmulf (rh ,drhd2,3,3,3,rh4)
!c...    R,alp*D + (R*R_0,1)*(R_0^T*D)
!        call matadd  (rh5(1,3),rh4(1,3),3,1,rd(1,2))
!        call matadd  (rh5(1,1),rh4(1,1),3,1,aid2(1,1))
!        call matadd  (rh5(1,2),rh4(1,2),3,1,aid2(1,2))
!        call matadd  (rh5(1,3),rh4(1,3),3,1,aid2(1,3))

c....   Gauss point tripod has to be updated 
        if (ilin.eq.2)   rh1= matmul (rh,Aini)

        
c....   determine M_I for all nodes based on rwnn and current a_i
        MatM=0.d0
        MatM_1=0.d0
        MatM_2=0.d0
        do l = 1,nen

            do i = 1,iMrow(l)
              do j = 1,iMcol(l)
                MatM(i,j,l)   = rwnGnp1(i,j,l)
              end do
            end do
       
        end do
          
      end if 
c
c...  T = W^T * H * T_3 and
c...  T_,alpha = W_,a^T*H*T3+W^T*H_,a*T3+W^T*H*T3_,a
c
c...  set T3 to unity and T3_,alpha to zero
      MatT3 = 0.0d0
      MatT3(1,1) = 1.0d0
      MatT3(2,2) = 1.0d0
      MatT3(3,3) = 1.0d0
      MatT3_1 = 0.0d0
      MatT3_2 = 0.0d0
      
c
c...  compute H and H_,alpha (for ilin.ne.2 H=1)
      !call updr63 (dwg(1:3,k),MatH,ilin,3)    !HIER?????? CHECKEN ob dwg(1,k) richtigß???
      !call updHalpha63(dwg(1,k),MatH_1,ilin,1)
      !call updHalpha63(dwg(1,k),MatH_2,ilin,2)
      ilinH= 0
c...  compute H and H_,alpha (for ilin.ne.2 H=1)
      call updr63 (dwg(1:3,k),MatH,ilinH,3)    !HIER?????? CHECKEN ob dwg(1,k) richtigß???
      call updHalpha63(dwg(1,k),MatH_1,ilinH,1)
      call updHalpha63(dwg(1,k),MatH_2,ilinH,2)
      MatH_1=0.d0
      MatH_2=0.d0

c...  compute W and W_,alpha
      if (ilin.eq.2) then
        call skew(MatW,rd(1,3))
        call skew(MatW_1,rd(1,1))
        call skew(MatW_2,rd(1,2))
      else
        call skew(MatW,r0(1,3))
        call skew(MatW_1,r0(1,1))
        call skew(MatW_2,r0(1,2))
      end if
c
c...  calculate T = W^T * H * T_3
      call mttmul (MatW,MatH,3,3,3, rh3)
      call matmulf(rh3,MatT3,3,3,3,MatT)
c...  calculate TM_I = W^T * H *T_3 *M_I and store it to w(1,1,I)
      call pzero(w,9*nen)
      do i = 1,nen
        call matmulf (MatT,MatM(1,1,i),3,3,3,w(1,1,i))
      end do
c
c...  calculate W^T * H * T_3,a
      call matmulf(rh3,MatT3_1,3,3,3,MatWHT_1)
      call matmulf(rh3,MatT3_2,3,3,3,MatWHT_2)
c...  calculate W_,a^T*H*T3
      call mttmul (MatW_1,MatH,3,3,3, rh3)
      call matmulf(rh3,MatT3,3,3,3, MatW_1HT)
      call mttmul (MatW_2,MatH,3,3,3, rh3)
      call matmulf(rh3,MatT3,3,3,3, MatW_2HT)
c...  calculate W^T*H_,a*T3
      call mttmul (MatW,MatH_1,3,3,3,rh3)
      call matmulf(rh3,MatT3,3,3,3,MatWH_1T)
      call mttmul (MatW,MatH_2,3,3,3,rh3)
      call matmulf(rh3,MatT3,3,3,3,MatWH_2T)
c
c...  sum up above terms for T_,alpha
      call matadd(MatWHT_1,MatW_1HT,3,3,rh3)
      call matadd(rh3,MatWH_1T,3,3,MatT_1)
      call matadd(MatWHT_2,MatW_2HT,3,3,rh3)
      call matadd(rh3,MatWH_2T,3,3,MatT_2)
      
c
c...  calculate ^T_I,alpha and store it to HutT_,alpha
      HutT_1=0.d0
      HutT_2=0.d0
      do i = 1,nen
c...    1-direction

          call matmulf(MatT_1,MatM(1,1,i),3,3,iMcol(i),rh3)
          call msmul(rh3,3,iMcol(i),shp(3,i))
          call matcop(w(1,1,i),3,iMcol(i),rh4)
          call msmul(rh4,3,iMcol(i),shp(1,i))
          call matmulf(MatT,MatM_1(1,1,i),3,3,iMcol(i),rh5)

        call msmul(rh5,3,iMcol(i),shp(3,i))
        call matadd(rh3,rh4,3,iMcol(i),rh6)
        call matadd(rh6,rh5,3,iMcol(i),HutT_1(1,1,i))
!        call matcop(rh4,3,2,HutT_1(1,1,i))
c...    2-direction

          call matmulf(MatT_2,MatM(1,1,i),3,3,iMcol(i),rh3)
          call msmul(rh3,3,iMcol(i),shp(3,i))
          call matcop(w(1,1,i),3,iMcol(i),rh4)
          call msmul(rh4,3,iMcol(i),shp(2,i))
          call matmulf(MatT,MatM_2(1,1,i),3,3,iMcol(i),rh5)

        call msmul(rh5,3,iMcol(i),shp(3,i))
        call matadd(rh3,rh4,3,iMcol(i),rh6)
        call matadd(rh6,rh5,3,iMcol(i),HutT_2(1,1,i))
!        call matcop(rh4,3,2,HutT_2(1,1,i))
      end do
c
c...  for G_ßß part of geometrical matrix
c...  compute T_HTM_I = H * T_3 * M_I and store it to twn(1,1,I) 

      call pzero(twn,9*nen)
      call pzero(twn_1,9*nen)
      call pzero(twn_2,9*nen)
      do i = 1,nen
          call matmulf (MatH,MatM(1,1,i),3,3,iMcol(i),twn(1,1,i))
      end do
c...  compute T_HTM_I,alpha and store it to twn_alpha(1,1,I)
c... 1-direction
      call matmulf(MatH_1,MatT3,3,3,3,rh3)
      call matmulf(MatH,MatT3_1,3,3,3,rh4)
      call matadd (rh3,rh4,3,3,rh5)
      do i= 1,nen

          rh3=0.d0
          do j = 1,iMrow(i)
            do l = 1,iMcol(i)
              rh3(j,l) = shp(1,i)*MatM(j,l,i)+shp(3,i)*MatM_1(j,l,i)
            end do
          end do
          call matmulf(rh5,MatM(1,1,i),3,3,iMcol(i),rh4)
          call msmul(rh4,3,iMcol(i),shp(3,i))
          call matmulf(MatH,rh3,3,3,iMcol(i),rh6)
          call matadd(rh4,rh6,3,iMcol(i),twn_1(1,1,i))

      end do
c... 2-direction
      call matmulf(MatH_2,MatT3,3,3,3,rh3)
      call matmulf(MatH,MatT3_2,3,3,3,rh4)
      call matadd (rh3,rh4,3,3,rh5)
      do i= 1,nen

          rh3=0.d0
          do j = 1,iMrow(i)
            do l = 1,iMcol(i)
              rh3(j,l) = shp(2,i)*MatM(j,l,i)+shp(3,i)*MatM_2(j,l,i)
            end do
          end do
          call matmulf(rh5,MatM(1,1,i),3,3,iMcol(i),rh4)
          call msmul(rh4,3,iMcol(i),shp(3,i))
          call matmulf(MatH,rh3,3,3,iMcol(i),rh6)
          call matadd(rh4,rh6,3,iMcol(i),twn_2(1,1,i))

      end do         
      return
      end      
c
      subroutine updRalpha63 (dwg,drw,ilin,idir)
c-----------------------------------------------------------------------
c     calculate current basis using rodriguez´ formula
c     this routine computes the derivation of R with respect to
c     local direction idir
c     works as well as Gruttmann Version
c-----------------------------------------------------------------------      
      implicit double precision (a-h,o-z)
      dimension dwg(9),drw(3,3),om(3,3),omom(3,3),omomalp(3,3),
     +          omalp(3,3),omalpom(3,3)

      tol = 1.d-16
      if (ilin.eq.1.or.ilin.eq.0) then
        call skew(drw,dwg(idir*3+1))
        return
      end if
  
c     compute Omega, Omega_,aplha and products of them
      call skew(om,dwg(1))
      call skew(omalp,dwg(idir*3+1))
      omom    = matmul (om,om   )
      omomalp = matmul (om,omalp)
      omalpom = matmul (omalp,om)
c     compute w and w_,alpha
      w     = dsqrt(dot(dwg(1),dwg(1),3))
      w2    = w*w
      w3    = w*w2
      sinw  = dsin(w)
      cosw  = dcos(w)
      walph = dot(dwg(idir*3+1),dwg(1),3) / w
c     compute factors c1-c4
      if (dabs(w).lt.tol) then
        c1 = 0.0d0
        c2 = 1.0d0
        c3 = 0.0d0
        c4 = 0.5d0
      else
        c1 = walph * (cosw*w-sinw) / w2
        c2 = sinw / w
        c3 = walph * (sinw*w-2.0d0+2.0d0*cosw) / w3
        c4 = (1.0d0-cosw) / w2
      end if
      do i= 1,3
        do j = 1,3
        drw(i,j) = c1*om(i,j)+c2*omalp(i,j)+c3*omom(i,j)+c4*omomalp(i,j)
     +           + c4*omalpom(i,j)
        end do
      end do
      return
      end

      subroutine updHalpha63 (dwg,drw,ilin,idir)
c-----------------------------------------------------------------------
c     calculate current basis using rodriguez´ formula
c     this routine computes the derivation of R with respect to
c     local direction idir
c     works as well as Gruttmann Version
c-----------------------------------------------------------------------      
      implicit double precision (a-h,o-z)
      dimension dwg(9),drw(3,3),om(3,3),omom(3,3),omomalp(3,3),
     +          omalp(3,3),omalpom(3,3)

      tol = 1.d-16
      if (ilin.eq.1.or.ilin.eq.0) then
        drw = 0.d0
        return
      end if
c     compute Omega, Omega_,aplha and products of them
      call skew(om,dwg(1))
      call skew(omalp,dwg(idir*3+1))
      omom    = matmul (om,om   )
      omomalp = matmul (om,omalp)
      omalpom = matmul (omalp,om)
c     compute w and w_,alpha
      w     = dsqrt(dot(dwg(1),dwg(1),3))
      w2    = w*w
      w3    = w*w2
      sinw  = dsin(w)
      cosw  = 1.0d0-dcos(w)
      walph = dot(dwg(idir*3+1),dwg(1),3) / w
c     compute factors c1-c4
      if (dabs(w).lt.tol) then
        c1 = 0.0d0
        c2 = 1.0d0/2.0d0
        c3 = 0.0d0
        c4 = 1.0d0/6.0d0
      else
        !HIERMIT SIND DIE RECHNUNGEN FÜR DIE CMAME PAPER GEMACHT; DIES IST ABER FALSCH
        !c1 = walph * (sinw*w-cosw*cosw) / w3 !HIER WAR DER FEHLER
        c1 = walph * (sinw*w-2.0d0*cosw) / w3
        c2 = cosw / w2
        c3 = walph * (cosw*w-3.0d0*w+3.0d0*sinw) / (w3*w)
        c4 = (w-sinw) / w3
      end if
      do i= 1,3
        do j = 1,3
        drw(i,j) = c1*om(i,j)+c2*omalp(i,j)+c3*omom(i,j)+c4*omomalp(i,j)
     +           + c4*omalpom(i,j)
        end do
      end do
      return
      end

      subroutine updRalpha63Gr (dwg,drw,ilin,idir)
c-----------------------------------------------------------------------
c     calculate current basis using rodriguez´ formula
c     this routine computes the derivation of R with respect to
c     local direction idir
c     see Gruttmann,Klinkel,Wagner 1995
c-----------------------------------------------------------------------      
      implicit double precision (a-h,o-z)
      dimension dwg(9),drw(3,3),om(3,3),omom(3,3),omomalp(3,3),
     +          omalp(3,3),omalpom(3,3),one(3,3)
      data one   / 1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0 /

      tol = 1.d-16
      if (ilin.eq.1.or.ilin.eq.0) then
        call skew(drw,dwg(idir*3+1))
        return
      end if

c     compute Omega, Omega_,aplha and products of them
      call skew(om,dwg(1))
      call skew(omalp,dwg(idir*3+1))
      call mtbmul(dwg(1),dwg(1),3,1,3, omom)
      call mtbmul(dwg(idir*3+1),dwg(1),3,1,3, omalpom)
      call mtbmul(dwg(1),dwg(idir*3+1),3,1,3, omomalp)
c     compute w and w_,alpha
      w     = dsqrt(dot(dwg(1),dwg(1),3))
      w2    = w*w
      w3    = w*w2
      sinw  = dsin(w)
      cosw  = dcos(w)
      walph = dot(dwg(idir*3+1),dwg(1),3) / w
c     compute factors c1-c4
      if (dabs(w).lt.tol) then
        c0 = 0.0d0
        c1 = 0.0d0
        c2 = 1.0d0
        c3 = 0.0d0
        c4 = 0.0d0
      else
        c0 = -sinw*walph
        c1 = walph * (cosw*w-sinw) / w2
        c2 = sinw / w
        c3 = walph * (sinw*w-2.0d0+2.0d0*cosw) / w3
        c4 = (1.0d0-cosw) / w2
      end if
      do i= 1,3
        do j = 1,3
        drw(i,j) = c1*om(i,j)+c2*omalp(i,j)+c3*omom(i,j)+c4*omomalp(i,j)
     +           + c4*omalpom(i,j)+c0*one(i,j)
        end do
      end do

      return
      end
c
      subroutine updr63 (dwn,rw,ilin,ifl)
c-----------------------------------------------------------------------
c     calculate current basis using rodriguez´ formula
c     ifl = 
c           1:  current basis  R(om)  
c           2:                 H^-1(om)         
c           3:                 H(om)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension   dwn(3),rw(3,3),om(3,3),om2(3,3), one(3,3)
      data one   / 1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0 /

      tol = 1.d-16
      theta = dsqrt(dot(dwn,dwn,3))
      !if (theta.ge.6.d0) then
      !write(*,*) theta
      !end if
      thet2 = theta*theta

      call skew (om,dwn)
      call matmulf (om,om, 3,3,3, om2)

      if (ifl .eq. 1) then                       ! R(om)
        if (dabs(theta) .lt. tol) then           ! theta -> 0
          c0 = 1.d0 
          c1 = 1.d0
          c2 = 1.d0/2.d0
        else
          c0 = 1.d0 
          c1 = dsin(theta)/theta
          c2 = (1.d0-dcos(theta))/thet2
        end if
        if(ilin.eq.0.or.ilin.eq.1)then      ! linear or moderate rotations
          c0 = 1.d0 
          c1 = 1.d0
          c2 = 0.d0
        endif  
      else if (ifl .eq. 2) then                  ! H^-1(om)  
        if (dabs( thet2*(dcos(theta)-1.d0) ) .lt. tol) then
          c0 = 1.d0 
          c1 =-1.d0/ 2.d0
          c2 = 1.d0/12.d0
        else
          c0 = 1.d0 
          c1 =-1.d0/2.d0
          c2 = 0.5d0* ( (2.d0*dcos(theta)+theta*dsin(theta)-2.d0)
     +                / (thet2*(dcos(theta)-1.d0)) )
        end if

      else if (ifl .eq. 3) then                  ! H(om)
        if (dabs(theta) .lt. tol) then           ! theta -> 0
          c0 = 1.d0 
          c1 = 1.d0/2.d0
          c2 = 1.d0/6.d0
        else
          c0 = 1.d0                    
          c1 = (1.d0 -dcos(theta))/thet2
          c2 = (theta-dsin(theta))/(theta*thet2)
        end if
        if(ilin.eq.0.or.ilin.eq.1)then     ! linear or moderate rotations               
          c0 = 1.d0 
          c1 = 0.d0
          c2 = 0.d0
        endif  
      end if
c
      do 10 i = 1,3
       do 10 j = 1,3
         rw(i,j) = c0*one(i,j) + c1*om(i,j) + c2*om2(i,j)
10    continue

      return
      end
c
      subroutine zerostr63(h1,h2,d,epstr,glu,sigp,cmat)
      USE iofile
c-----------------------------------------------------------------------
c     zero thickness normal stress iteration
c     im = 1:         Pegasus method       2: Newton´s method
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension h1(*),h2(*),d(*),ign(6),epstr(*),sigp(*),cmat(5,*)
     2         ,e3d(6),sig3d(7),Cmat3d(6,6),glu(3,3)
      
c      include 'iofile.h'
      data ign/1,2,4,5,6,3/
c
c.....setup starting values
c
      tol = 1.d-7
      iter = 0
      itermax = 25
      e3d(1) = epstr(1)
      e3d(2) = epstr(2)
      e3d(4) = epstr(3)
      e3d(5) = epstr(4)
      e3d(6) = epstr(5)
      e3d(3) = h1(8)
      teste  = dsqrt(dot(e3d,e3d,6))
c
      im = 2
      if(dabs(e3d(1)+e3d(2)).lt.1.d-15)im=2       
      if(im.eq.1)then
c
ccccccccccccccccccccc PEGASUS Method ccccccccccccccccccccccccccccccccccc
   
      xnu  = d(41)  
      x1   =  -xnu*(e3d(1)+e3d(2))        ! starting values
      x2   =  -    (e3d(1)+e3d(2))

      e3d(3) = x1
      call matelib63(h1,h2,d,e3d,glu,sig3d,cmat3d)
      f1   = sig3d(3)
      
      e3d(3) = x2
      call matelib63(h1,h2,d,e3d,glu,sig3d,cmat3d)
      f2   = sig3d(3)
cfg         if(f1*f2.gt.0)write(*,*)'Starting values Pegasus not ok'

cfgc.....   find admissible range
cfgc
cfg      de33 = 1.d-1
cfg
cfg      call matelib63(h1,h2,d,e3d,glu,sig3d,cmat3d)
cfg      do n = 1,20 
cfg       a  = e3d(3)
cfg       fa = sig3d(3)
cfg       e3d(3) = e3d(3) + de33
cfg       call matelib63(h1,h2,d,e3d,glu,sig3d,cmat3d)
cfg       fb = sig3d(3)
cfg       if(((fa*fb).gt.0.d0).and.(dabs(fb).gt.dabs(fa)))go to 1
cfg       if(fa*fb.lt.0.d0)go to 2
cfg      end do
cfg
cfg1     e3d(3) = h1(8)
cfg      call matelib63(h1,h2,d,e3d,glu,sig3d,cmat3d)
cfg      do n = 1,20 
cfg       b  = e3d(3)
cfg       fa = sig3d(3)
cfg       e3d(3) = e3d(3) - de33
cfg       call matelib63(h1,h2,d,e3d,glu,sig3d,cmat3d)
cfg       fb = sig3d(3)
cfg       if(fa*fb.lt.0.d0)go to 3
cfg      end do
cfg
cfg2     b = e3d(3)
cfg      go to 4    
cfg3     a = e3d(3)
cfg4     if(n.eq.10)write(*,*)'admissible range not found'
cfgc
cfgc      
cfg      x1   = a
cfg      e3d(3) = x1
cfg      call matelib63(h1,h2,d,e3d,glu,sig3d,cmat3d)
cfg      f1   = sig3d(3)
cfg      
cfg      x2   = b
cfg      e3d(3) = x2
cfg      call matelib63(h1,h2,d,e3d,glu,sig3d,cmat3d)
cfg      f2   = sig3d(3)
c....                           start of iteration
5     iter = iter + 1
      if (iter.gt.25) then
         write(  *,*) ' No convergence in zero normal stress iteration'
         write(iow,*) ' No convergence in zero normal stress iteration'
         if(f1*f2.gt.0)write(*,*)'Starting values Pegasus not ok'
         stop
      end if
        x3 = x2 - ((x2-x1)/(f2-f1))*f2
        e3d(3) = x3
        call matelib63(h1,h2,d,e3d,glu,sig3d,cmat3d)
        Czz = 1.d0/Cmat3d(3,3)       
        f3 = sig3d(3)
        sigv = dsqrt(sig3d(1)*sig3d(1)+sig3d(2)*sig3d(2))
        if(teste.lt.1.e-20)sigv=1.d0
        if (dabs(sig3d(3))/sigv.le.tol) go to 100

      if(f2*f3.lt.0.d0)then
         x1 = x2
         f1 = f2
      else
         f1 = f1*f2/(f2+f3) 
      end if
         x2 = x3
         f2 = f3
      go to 5  


cccccccccccccccccc NEWTON Iteration  ccccccccccccccccccccccccccccccccccc
c
      elseif(im.eq.2)then

101   iter = iter + 1
      if (iter.gt.itermax) then
         write(  *,*) ' No convergence in zero normal stress iteration'
         write(iow,*) ' No convergence in zero normal stress iteration'
         stop
      end if
      call matelib63(h1,h2,d,e3d,glu,sig3d,cmat3d)

      Czz = 1.d0/Cmat3d(3,3)
      e3d(3) = e3d(3) - sig3d(3)*Czz
      sigv = dsqrt(sig3d(1)*sig3d(1)+sig3d(2)*sig3d(2))
      if(teste.lt.1.e-20)sigv=1.d0
      if (dabs(sig3d(3))/sigv.le.tol) go to 100
      go to 101    
      endif 
c.....end of Newton iteration
c
c
c     store e33, sig, condensation of Cmat3d -> cmat
100   continue
      do i=1,5
       ii = ign(i)
       sigp(i) = sig3d(ii)
       do j=1,5
        jj = ign(j)
        cmat(i,j) = Cmat3d(ii,jj) - Cmat3d(ii,3)*Czz*Cmat3d(3,jj)  
       end do
      end do
c
      h2(8) = e3d(3)
c  
      sigp(6) = sig3d(7)
c
      return
      end
c
      subroutine matelib63(h1,h2,d,e,glu,sig,cmat)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension h1(*),h2(*),d(*),e(*),glu(3,3),sig(*),Cmat(6,*)
c
      imat = d(15)
c
      if (imat.eq.1) then
          call mate0163(h1,h2,d,e,glu,sig,cmat)
      elseif(imat.eq.2) then
          call mate0263(h1,h2,d,e,glu,sig,cmat)
      elseif(imat.eq.3) then
          call mate0363(h1,h2,d,e,glu,sig,cmat)
      end if
c
      return
      end
c
      subroutine mate0163(h1,h2,d,e,glu,sig,cmat)
      USE eldata
      USE iofile
c-----------------------------------------------------------------------
c     input :  h1,nh,d,glu
c     output:  h2,sig,cmat
c
c     h1(nh), h2(nh)                   history arrays
c     nh = 7 + 1  = 8                  1 hist. variable is necessary
c                                      for local iteration 
c                                      7 hist. variables are necessary
c                                      for plasticity 
c     d(*)                             contains material data
c                 40  E                Young's modulus
c                 41  v                Poisson's ratio
c                 42  y_o              yield stress
c                 43  y_i              yield stress infty
c                 44  xk               linear hardening modulus
c                 45  xd               exponential harderning
c     e(6), sig(6), cmat(6,6)          strains, stresses, tangent matrix
c     glu(3,3)                         base vectors of e in columns
c
c
c...  storage of epsilon and sigma
c     11  22  33  12  13  23
c
c...  storage of elasticity matrix
c     1111  1122  1133  1112  1113  1123  
c     2211  2222  2233  2212  2213  2223  
c     3311  3322  3333  3312  3313  3323  
c     1211  1222  1233  1212  1213  1223  
c     1311  1322  1333  1312  1313  1323 
c     2311  2322  2333  2312  2313  2323 
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),h1(*),h2(*),sig(*),cmat(6,*),e(*),ep(6),
     1          ee(6),sd(6),S(6),glu(3,3),ec(6),gko(3,3),glo(3,3),
     2          Te(6,6),to(3,3),t(3,3),Tet(6,6),xhelp(6,6),cort(6,6)
c    
c      include 'eldata.h'
c      include 'iofile.h'
      xE   = d(40)
      xn   = d(41)
      y_o  = d(42)
      y_i  = d(43)
      xk   = d(44)
      xd   = d(45)
c
c     get history variables at time t_n 
      if (dabs(y_o).lt.1.d-18)then
       call pzero(ep,6)
       a = 0.d0
      else
       do j = 1,6
        ep(j) = h1(j)
       end do
       a = h1(7)
      endif
c
c.....coefficient matrix
c      do 15 i=1,3
c         do 15 j=1,3
c            gko(i,j) = 0.d0
c            do 15 k=1,3
c            gko(i,j) = gko(i,j) + glu(k,i) * glu(k,j)
c15    continue
c      call invert(gko,3,3)
cc.....contravariant basis
c      do i=1,3
c         glo(i,1)=gko(1,1)*glu(i,1)+gko(1,2)*glu(i,2)+gko(1,3)*glu(i,3)
c         glo(i,2)=gko(2,1)*glu(i,1)+gko(2,2)*glu(i,2)+gko(2,3)*glu(i,3)
c         glo(i,3)=gko(3,1)*glu(i,1)+gko(3,2)*glu(i,2)+gko(3,3)*glu(i,3)
c      end do
cc.... transformation of E from convec into cart
c      call pzero(t,3*3)
c      do i=1,3
c         t(i,i)=1.d0
c      end do
c      do 20 l=1,3
c         do 20 m=1,3
c            to(l,m)=0.0d0
c            do 20 i=1,3
c               to(l,m)=to(l,m)+glo(i,l)*t(i,m)
c20    continue
cc     matrix Te=to(i,j)*to(k,l)
c      Te(1,1) =       to(1,1)*to(1,1)
c      Te(1,2) =       to(2,1)*to(2,1)
c      Te(1,3) =       to(3,1)*to(3,1)
c      Te(1,4) =       to(1,1)*to(2,1)
c      Te(1,5) =       to(1,1)*to(3,1)
c      Te(1,6) =       to(2,1)*to(3,1)
c      Te(2,1) =       to(1,2)*to(1,2)
c      Te(2,2) =       to(2,2)*to(2,2)
c      Te(2,3) =       to(3,2)*to(3,2)
c      Te(2,4) =       to(1,2)*to(2,2)
c      Te(2,5) =       to(1,2)*to(3,2)
c      Te(2,6) =       to(2,2)*to(3,2)
c      Te(3,1) =       to(1,3)*to(1,3)
c      Te(3,2) =       to(2,3)*to(2,3)
c      Te(3,3) =       to(3,3)*to(3,3) 
c      Te(3,4) =       to(1,3)*to(2,3)
c      Te(3,5) =       to(1,3)*to(3,3)
c      Te(3,6) =       to(2,3)*to(3,3)
c      Te(4,1) = 2.0d0*to(1,1)*to(1,2)
c      Te(4,2) = 2.0d0*to(2,1)*to(2,2)
c      Te(4,3) = 2.0d0*to(3,1)*to(3,2)
c      Te(4,4) =       to(1,1)*to(2,2) + to(2,1)*to(1,2)
c      Te(4,5) =       to(1,1)*to(3,2) + to(3,1)*to(1,2)
c      Te(4,6) =       to(2,1)*to(3,2) + to(3,1)*to(2,2)
c      Te(5,1) = 2.0d0*to(1,1)*to(1,3)           
c      Te(5,2) = 2.0d0*to(2,1)*to(2,3)
c      Te(5,3) = 2.0d0*to(3,1)*to(3,3)
c      Te(5,4) =       to(1,1)*to(2,3) + to(2,1)*to(1,3)
c      Te(5,5) =       to(1,1)*to(3,3) + to(3,1)*to(1,3)
c      Te(5,6) =       to(2,1)*to(3,3) + to(3,1)*to(2,3)
c      Te(6,1) = 2.0d0*to(1,2)*to(1,3)
c      Te(6,2) = 2.0d0*to(2,2)*to(2,3)
c      Te(6,3) = 2.0d0*to(3,2)*to(3,3)
c      Te(6,4) =       to(1,2)*to(2,3) + to(2,2)*to(1,3)
c      Te(6,5) =       to(1,2)*to(3,3) + to(3,2)*to(1,3)
c      Te(6,6) =       to(2,2)*to(3,3) + to(3,2)*to(2,3)
c.....transformation from the isopar. into cartes. space
c      call pzero(ec,6)
c      call mvmul(Te,e,6,6,ec)
c      
c.... trial strains  E^el = E - E^pl
      do i=1,6
         ee(i)=0.d0
c         ee(i)=ec(i)-ep(i)
         ee(i)=e(i)-ep(i)
      end do
c
c.....elasticity matrix
      call mat0163(cort,xE,xn)
c
c.....trial stresses   sig(i)=cmat(i,j)*ee(j)
      call mvmul(cort,ee,6,6,sig)
      if (abs(y_o).lt.1.d-18) go to 100
c
c.....yield stress 
      y = y_o + xk*a + (y_i-y_o)*(1-exp(-xd*a))
c
c.... compute g^trial
      onethird = 1.0d0/3.0d0
      twothird = 2.0d0/3.0d0
      sd(1) =  twothird*sig(1)-onethird*sig(2)-onethird*sig(3)
      sd(2) = -onethird*sig(1)+twothird*sig(2)-onethird*sig(3)
      sd(3) = -onethird*sig(1)-onethird*sig(2)+twothird*sig(3)
      sd(4) = sig(4)
      sd(5) = sig(5)
      sd(6) = sig(6)
      g_tr =dsqrt( 3.d0/2.d0*(sd(1)*sd(1)+sd(2)*sd(2)+sd(3)*sd(3))
     +            +     3.d0*(sd(4)*sd(4)+sd(5)*sd(5)+sd(6)*sd(6)) ) 
c      
c.....Yield condition   f = g^trial - y
      f = g_tr - y
      if(f.gt.1.d-10*y) call plas0163(d,sig,sd,f,g_tr,y,a,ep,cort)

      do i = 1,6
       h2(i) = ep(i)
      end do
      h2(7) = a

100   continue
c
c.....transpose Te(i,j) 
c      do 40 i=1,6
c        do 40 j=1,6
c        Tet(j,i)=Te(i,j)
c40    continue
c
c.....transformation of cort into cmat (cart. into convect.)
c      call pzero(xhelp,6*6)
c      call pzero(cmat,6*6)
c      call matmulf(Tet,cort,6,6,6,xhelp)
c      call matmulf(xhelp,Te,6,6,6,cmat)
c
c.....transformation of stresses
c      do i=1,6
c         S(i)=sig(i)
c         sig(i)=0.d0
c      end do
c      call mvmul(Tet,S,6,6,sig)
c
      call matcop(cort,6,6,cmat)
c
c
      sig(7) = dsqrt(0.5d0*((sig(1)-sig(2))*(sig(1)-sig(2))
     +                     +(sig(1)-sig(3))*(sig(1)-sig(3))
     +                     +(sig(2)-sig(3))*(sig(2)-sig(3)))
     +               + 3.d0*(sig(4)*sig(4)+sig(5)*sig(5)+sig(6)*sig(6)))
c
      return
      end
c
      subroutine plas0163(d,sig,sd,f,g_tr,y,a,ep,cmat)
      USE iofile
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),sig(*),sd(*),ep(*),cmat(6,6)
c      include 'iofile.h'
c
      xk = d(40)/(3.d0*(1.d0-2.d0*d(41)))
      xm = d(40)/(2.d0*(1.d0+d(41)))
      y_o  = d(42)
      y_i  = d(43)
      xh   = d(44)
      xd   = d(45)
c.....Startvalues
      gam = 0.d0
      newton = 0
c
c.....begin local Newton iteration...................................
100   continue
      newton = newton + 1 
      if (newton.gt.33) then
         write( * ,*) 'no convergence in plas0163 '
         write(iow,*) 'no convergence in plas0163'
         stop
      end if
      dgam = f /(3.d0*xm + xh + (y_i-y_o)*xd*exp(-xd*a))
      gam = gam + dgam
      a = a + dgam
      y = y_o + xh*a + (y_i-y_o)*(1.d0 - exp(-xd*a))
      f = g_tr - 3.d0*xm*gam - y

cfg      if(newton.gt.1)write(*,*)newton,f,1.d-10*y,a

c
      if (abs(f).gt.1.d-10*y) go to 100 
c.....end local Newton iteration.....................................
c
c.....elasto-plastic tangent operator
      call pzero(Cmat,6*6)
      dyda = xh + (y_i-y_o)*xd*exp(-xd*a)
      xnt = dsqrt(       sd(1)*sd(1)+sd(2)*sd(2)+sd(3)*sd(3)+
     +             2.d0*(sd(4)*sd(4)+sd(5)*sd(5)+sd(6)*sd(6)) )
      beta1 = 2.d0*xm*y/g_tr
      beta2 =((1.d0/(dyda/(3.d0*xm)+1.d0) + y/g_tr - 1.d0)*2.d0*xm)/
     +       (xnt*xnt)
      cmat(1,1) = xk + beta1*2.d0/3.d0
      cmat(1,2) = xk - beta1*1.d0/3.d0
      cmat(1,3) = xk - beta1*1.d0/3.d0
      cmat(2,1) = xk - beta1*1.d0/3.d0
      cmat(2,2) = xk + beta1*2.d0/3.d0
      cmat(2,3) = xk - beta1*1.d0/3.d0
      cmat(3,1) = xk - beta1*1.d0/3.d0
      cmat(3,2) = xk - beta1*1.d0/3.d0
      cmat(3,3) = xk + beta1*2.d0/3.d0
      cmat(4,4) = beta1*0.5d0
      cmat(5,5) = beta1*0.5d0   
      cmat(6,6) = beta1*0.5d0   
      do i=1,6
         do j=1,6
            cmat(i,j) = cmat(i,j) - beta2*sd(i)*sd(j)
         end do
      end do
c
c.....update stresses and plastic strains
      trsig = (sig(1)+sig(2)+sig(3))/3.d0
      fac = 1.d0
      do i=1,6
         sig(i) = y/g_tr * sd(i)
         if (i.le.3) sig(i) = sig(i) + trsig
         if (i.gt.3) fac = 2.d0
         ep(i)= ep(i) + 3.d0/2.d0*gam/g_tr*sd(i)*fac
      end do
      return
      end
c
      subroutine mat0163(cmat,e1,xnu)
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension cmat(6,6)
c
      call pzero(cmat,6*6)
c
c.... local isotropic elasticity matrix
      cc = e1/((1.0d0+xnu)*(1.0d0-2.0d0*xnu))
      tt = (1.0d0-2.0d0*xnu)/2.0d0
c
c.... material parameters 
      cmat(1,1) = cc*(1.0d0-xnu)
      cmat(1,2) = cc*xnu
      cmat(1,3) = cc*xnu
      cmat(2,2) = cc*(1.0d0-xnu)
      cmat(2,3) = cc*xnu
      cmat(3,3) = cc*(1.0d0-xnu)
      cmat(4,4) = cc*tt
      cmat(5,5) = cc*tt     
      cmat(6,6) = cc*tt     
c
c.... symmetry
      cmat(2,1) = cmat(1,2)
      cmat(3,1) = cmat(1,3)
      cmat(3,2) = cmat(2,3)
c
      return
      end
c
      subroutine mate0263(h1,h2,d,e,glu,S,cmat)
      USE iofile
      USE eldata
c-----------------------------------------------------------------------
c     finite strain isotropic J2 plasticity, nonl. isotropic hardening     
c
c     input data :  h1,d,e,glu
c     output data:  h2,S,cmat
c
c     h1(nh), h2(nh)                    history array
c     nh = 7 + 1  = 8                   1 hist. variable is necessary
c                                       for local iteration 
c                                       7 hist. variables are necessary
c                                       for plasticity 
c     d(*)                              contains material data
c                 40  E                 Youngs modulus
c                 41  v                 poissons ratio
c                 42  y_o               yield stress
c                 43  y_i               yield stress infty
c                 44  xk                linear hardening modulus
c                 45  xd                exponential hardening
c     e(6), S(6), cmat(6,6)             strains, stresses, tangent modulus
c     glu(3,3)                          base vectors of E in columns

c....  parameters     in h1
c      epn(6)    - plastic strains            at t_n
c      alpha     - equivalent plastic strain  at t_n
c
c
c...  storage of epsilon and sigma
c     11  22  33  12  13  23
c
c...  storage of elasticitymatrix
c     1111  1122  1133  1112  1131  1123
c     2211  2222  2233  2212  2231  2223
c     3311  3322  3333  3312  3331  3323
c     1211  1222  1233  1212  1231  1223
c     3111  3122  3133  3112  3131  3123
c     2311  2322  2333  2312  2331  2323
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      double precision L1(3,3),L2(3,3)
      dimension d(*),S(*),cmat(6,6)
      dimension G(3,3),glu(3,3),e(6),
     +   C(3,3),C_n(3,3),h1(*),h2(*),Cp(3,3),Cp_n(3,3),
     +   xlam(3),xN(3,3),fv1(3),fv2(3),
     +   ee(3),ep(3),e_tr(3),
     +   Ce(3,3),T_tr(3),Td_tr(3),sig(3),
     +   T1(6,3),T2(6,3),T1L1(6,3),T2L2(6,3),T1L1T1(6,6),T2L2T2(6,6)
c
c      include 'eldata.h'
c      include 'iofile.h'
c
c.... elasticity constants 
      xkappa  = d(40)/(3.d0*(1.d0-2.d0*d(41)))
      xmu  = d(40)/(2.d0*(1.d0+d(41)))
      y_o  = d(42)
      y_i  = d(43)
      xk   = d(44)
      xd   = d(45)
c
c.....covariant metric coefficients 
      do i=1,3
         do j=1,3
            G(i,j) = 0.d0
            do k=1,3
               G(i,j) = G(i,j) + glu(k,i) * glu(k,j)
            end do
         end do
      end do
      
      if((G(3,3)-1.d0).gt.1.d-15)write(*,*)'!!! G33 = ',G(3,3)

c
c.....Right Cauchy Green strain tensor
      do i=1,3
        C(i,i) = 2.d0*e(i) + G(i,i)
      end do
      C(1,2) = e(4) + G(1,2)
      C(1,3) = e(5) + G(1,3)
      C(2,3) = e(6) + G(2,3)
      C(2,1) = C(1,2) 
      C(3,1) = C(1,3)
      C(3,2) = C(2,3)
        
c
c.....covariant coefficients of the Right--Cauchy--Green tensor C
c     covariant coefficients of Plastic Metric tensor Cp
      do i = 1,3
         Cp(i,i) = h1(i)
      end do
      Cp(1,2) = h1(4)
      Cp(1,3) = h1(5)
      Cp(2,3) = h1(6)
      Cp(2,1) = h1(4)
      Cp(3,1) = h1(5)
      Cp(3,2) = h1(6)
      a       = h1(7)  
      detCp =  Cp(1,1)*(Cp(2,2)*Cp(3,3)-Cp(2,3)*Cp(3,2))
     +        -Cp(1,2)*(Cp(2,1)*Cp(3,3)-Cp(2,3)*Cp(3,1))
     +        +Cp(1,3)*(Cp(2,1)*Cp(3,2)-Cp(2,2)*Cp(3,1))
      detC  =  C(1,1)* (C(2,2)*C(3,3)-C(2,3)*C(3,2))
     +        -C(1,2)* (C(2,1)*C(3,3)-C(2,3)*C(3,1))
     +        +C(1,3)* (C(2,1)*C(3,2)-C(2,2)*C(3,1))
      if (detC.lt.1.e-14) then
         write(*,*) 'In mate02 det C < 0 ', detC
         write(iow,*) 'In mate02 det C <= 0 '
         stop
      end if
      do i=1,3
         do j=1,3
            if (abs(detCp).lt.1.e-14 ) then
               Cp(i,j) = G(i,j)
            end if
            if(dabs(C(i,j)).lt.1.e-13) C(i,j) = 0.d0
            if(dabs(Cp(i,j)).lt.1.e-13) Cp(i,j) = 0.d0
            Cp_n(i,j) = Cp(i,j)
            C_n(i,j) = C(i,j)
         end do
      end do
c
c.....Eigenvalue problem [C - lam_tr Cp] N = 0
      call pzero(xlam,3)
      call pzero(xN,3*3)
      call pzero(fv1,3)
      call pzero(fv2,3)
      ierr = 0
      call rsg(3,3,C,Cp,xlam,3,xN,fv1,fv2,ierr)
      do i=1,3
         if (xlam(i).lt.1.e-14) then
            write(*,*) 'xlam(',i,')=',xlam(i)
            stop
         end if
      end do
      if (abs(xlam(1)-xlam(2)).lt.1.e-15) xlam(2)=xlam(1)+1.e-15
      if (abs(xlam(1)-xlam(3)).lt.1.e-15) xlam(3)=xlam(1)+1.e-15
      if (abs(xlam(2)-xlam(3)).lt.1.e-15) xlam(3)=xlam(2)+1.e-15
c
c.....test of the eigenvectors
c     xN * (Cp xN) = 1
      do k=1,3
         xx = 0.d0
         do i=1,3
            xx = xx + xN(i,k)*
     +          (Cp_n(i,1)*xN(1,k)+Cp_n(i,2)*xN(2,k)+Cp_n(i,3)*xN(3,k))
         end do
         if (dabs(xx)-1.d0.gt.1.e-14) then
            write(*,*) 'N(',k,') Cp N(',k,') >1'
            stop
         end if
      end do
c
c.....logarithmic strains
      do i=1,3
         e_tr(i) = 0.5d0 * dlog(xlam(i))
      end do
c
c.....trial values
      do i=1,3
         do j=1,3
            fac=2.d0
            if(i.ne.j) fac=-1.d0
            Ce(i,j) = xkappa + 2.d0 * xmu * fac/3.d0
         end do
      end do
      do i=1,3
         T_tr(i) = 0.d0
         do j=1,3
            T_tr(i) = T_tr(i) + Ce(i,j) * e_tr(j)
         end do
      end do
      Td_tr(1)= 2.d0/3.d0*T_tr(1)-1.d0/3.d0*T_tr(2)-1.d0/3.d0*T_tr(3)
      Td_tr(2)=-1.d0/3.d0*T_tr(1)+2.d0/3.d0*T_tr(2)-1.d0/3.d0*T_tr(3)
      Td_tr(3)=-1.d0/3.d0*T_tr(1)-1.d0/3.d0*T_tr(2)+2.d0/3.d0*T_tr(3)
      g_tr = dsqrt(3.d0/2.d0*
     +       (Td_tr(1)*Td_tr(1)+Td_tr(2)*Td_tr(2)+Td_tr(3)*Td_tr(3)))
      y = y_o + xk*a + (y_i-y_o)*(1-exp(-xd*a))
      f = g_tr - y
c
c.....plastic corrector
      if (y_o.gt.1e-10.and.f.gt.1.d-10*y) then
         call plas0263(d,T_tr,Td_tr,g_tr,f,Ce,ep,a)
c.....   compute Cp_n+1
         do i=1,3
            ee(i) = e_tr(i) - ep(i)
         end do
         call pzero(Cp,3*3)
         do k=1,3
         fac = dexp(2.d0*ee(k)) / dexp(2.d0*e_tr(k))
            do i=1,3
               do j=1,3
               Cp(i,j) = Cp(i,j) + fac * xN(i,k)*xN(j,k)
               end do
            end do
         end do
         call invert(Cp,3,3)
      end if
      h2(1) = Cp(1,1)
      h2(2) = Cp(2,2)
      h2(3) = Cp(3,3)
      h2(4) = Cp(1,2)
      h2(5) = Cp(1,3)
      h2(6) = Cp(2,3)
      h2(7) = a
c
c.....von Mises stress
        S(7) = dsqrt(0.5d0*((T_tr(1)-T_tr(2))*(T_tr(1)-T_tr(2))
     1                     +(T_tr(1)-T_tr(3))*(T_tr(1)-T_tr(3))
     2                     +(T_tr(2)-T_tr(3))*(T_tr(2)-T_tr(3))))
c.....transformation of the stresses
      do i=1,3
         sig(i) = T_tr(i)/xlam(i)
         T1(1,i) = xN(1,i)*xN(1,i)
         T1(2,i) = xN(2,i)*xN(2,i)
         T1(3,i) = xN(3,i)*xN(3,i)
         T1(4,i) = xN(1,i)*xN(2,i)
         T1(5,i) = xN(1,i)*xN(3,i)
         T1(6,i) = xN(2,i)*xN(3,i)
      end do
      do i = 1,6
         S(i) = 0.d0
         do j=1,3
            S(i) = S(i) + T1(i,j) * sig(j)
         end do
      end do
c
cfg      if(itp.eq.0)return
      
c.....transformation of the elastoplastic tangent
      T2(1,1) = xN(1,1)*xN(1,2) + xN(1,2)*xN(1,1)
      T2(2,1) = xN(2,1)*xN(2,2) + xN(2,2)*xN(2,1)
      T2(3,1) = xN(3,1)*xN(3,2) + xN(3,2)*xN(3,1)
      T2(4,1) = xN(1,1)*xN(2,2) + xN(1,2)*xN(2,1)
      T2(5,1) = xN(1,1)*xN(3,2) + xN(1,2)*xN(3,1)
      T2(6,1) = xN(2,1)*xN(3,2) + xN(2,2)*xN(3,1)
      T2(1,2) = xN(1,1)*xN(1,3) + xN(1,3)*xN(1,1)
      T2(2,2) = xN(2,1)*xN(2,3) + xN(2,3)*xN(2,1)
      T2(3,2) = xN(3,1)*xN(3,3) + xN(3,3)*xN(3,1)
      T2(4,2) = xN(1,1)*xN(2,3) + xN(1,3)*xN(2,1)
      T2(5,2) = xN(1,1)*xN(3,3) + xN(1,3)*xN(3,1)
      T2(6,2) = xN(2,1)*xN(3,3) + xN(2,3)*xN(3,1)
      T2(1,3) = xN(1,2)*xN(1,3) + xN(1,3)*xN(1,2)
      T2(2,3) = xN(2,2)*xN(2,3) + xN(2,3)*xN(2,2)
      T2(3,3) = xN(3,2)*xN(3,3) + xN(3,3)*xN(3,2)
      T2(4,3) = xN(1,2)*xN(2,3) + xN(1,3)*xN(2,2)
      T2(5,3) = xN(1,2)*xN(3,3) + xN(1,3)*xN(3,2)
      T2(6,3) = xN(2,2)*xN(3,3) + xN(2,3)*xN(3,2)
c
      do i=1,3
         L1(i,i) = ( Ce(i,i) - 2.d0*T_tr(i) ) / (xlam(i)*xlam(i))
      end do
      L1(1,2) = Ce(1,2) / (xlam(1)*xlam(2))
      L1(1,3) = Ce(1,3) / (xlam(1)*xlam(3))
      L1(2,3) = Ce(2,3) / (xlam(2)*xlam(3))
      L1(2,1) = L1(1,2)
      L1(3,1) = L1(1,3)
      L1(3,2) = L1(2,3)
c
      call pzero(L2,3*3)
      L2(1,1) = (Sig(1)-Sig(2))/(xlam(1)-xlam(2))
      L2(2,2) = (Sig(1)-Sig(3))/(xlam(1)-xlam(3))
      L2(3,3) = (Sig(2)-Sig(3))/(xlam(2)-xlam(3))
c
      do i=1,6
         do j=1,3
            T1L1(i,j) = 0.d0
            T2L2(i,j) = 0.d0
            do k=1,3
            T1L1(i,j) = T1L1(i,j) + T1(i,k)*L1(k,j)
            T2L2(i,j) = T2L2(i,j) + T2(i,k)*L2(k,j)
            end do
         end do
      end do
      do i=1,6
         do j=1,6
            T1L1T1(i,j) = 0.d0
            T2L2T2(i,j) = 0.d0
            do k=1,3
               T1L1T1(i,j) = T1L1T1(i,j) + T1L1(i,k)*T1(j,k)
               T2L2T2(i,j) = T2L2T2(i,j) + T2L2(i,k)*T2(j,k)
            end do
            cmat(i,j) = T1L1T1(i,j) + T2L2T2(i,j)
         end do
      end do
c
      return
      end
c
      subroutine plas0263(d,T_tr,Td_tr,g_tr,f,Ce,ep,a)
      USE iofile
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),T_tr(3),Td_tr(3),Ce(3,3),ep(3)
c      include 'iofile.h'
c
      xkappa  = d(40)/(3.d0*(1.d0-2.d0*d(41)))
      xmu  = d(40)/(2.d0*(1.d0+d(41)))
      y_o  = d(42)
      y_i  = d(43)
      xh   = d(44)
      xd   = d(45)
c.....Startvalues
      gam = 0.d0
      newton = 0
c
c.....begin local Newton iteration.......................
100   continue
      newton = newton + 1
      if (newton.gt.33) then
         write( * ,*) 'no convergence in plas0263'
         write( * ,*) 'f=',f
         write(iow,*) 'no convergence in plas0263'
         stop
      end if
      dgam = f/(3.d0*xmu + xh + (y_i-y_o)*xd*exp(-xd*a))
      gam = gam + dgam
      a = a + dgam
      y = y_o + xh*a + (y_i-y_o)*(1.d0 - exp(-xd*a))
      f = g_tr - 3.d0*xmu*gam - y
      if (dabs(f).gt.1.d-10*y) go to 100
      gam = gam - 1.e-8
c.....end local Newton iteration.........................
c
c.....update stresses and plastic strains
      g = g_tr-3.d0*xmu*gam
      trT = (T_tr(1)+T_tr(2)+T_tr(3))/3.d0
      do i=1,3
         T_tr(i) = g/g_tr * Td_tr(i) + trT
         ep(i) = 3.d0/2.d0*gam/g_tr*Td_tr(i)
      end do
c
c.....elasto-plastic tangent operator
      dyda = xh + (y_i-y_o)*xd*exp(-xd*a)
      xnt = Td_tr(1)*Td_tr(1)+Td_tr(2)*Td_tr(2)+Td_tr(3)*Td_tr(3)
csk      beta1 = 2.d0*xmu*g/g_tr
      beta1 = 2.d0*xmu /(1.d0 + 3.d0 * gam * xmu / y)
      beta2 =((1.d0/(dyda/(3.d0*xmu)+1.d0) + g/g_tr - 1.d0)*2.d0*xmu)/
     +       xnt
      Ce(1,1) = xkappa + beta1*2.d0/3.d0
      Ce(1,2) = xkappa - beta1*1.d0/3.d0
      Ce(1,3) = xkappa - beta1*1.d0/3.d0
      Ce(2,1) = xkappa - beta1*1.d0/3.d0
      Ce(2,2) = xkappa + beta1*2.d0/3.d0
      Ce(2,3) = xkappa - beta1*1.d0/3.d0
      Ce(3,1) = xkappa - beta1*1.d0/3.d0
      Ce(3,2) = xkappa - beta1*1.d0/3.d0
      Ce(3,3) = xkappa + beta1*2.d0/3.d0
      do i=1,3
         do j=1,3
            Ce(i,j) = Ce(i,j) - beta2*Td_tr(i)*Td_tr(j)
         end do
      end do
c
      return
      end
c
      subroutine mate0363(h1,h2,d,e,glu,S,cmat)
      USE eldata
      USE iofile
c-----------------------------------------------------------------------
c     Ogden´s material law 
c
c     input data:  h1,nh,Epsc,glu
c     output data: h2,S,cmat
c     h1(8), h2(8)                     history array 
c                                      1-7 not used, 8 = e33
c
c     d(*)                             contains material data
c                 40  mue_1            shear modulus
c                 41  mue_2            shear modulus
c                 42  mue_3            shear modulus
c                 43  alpha_1          ! sum(mue_i*alpha_i) = 2mue 
c                 44  alpha_2          
c                 45  alpha_3          
c                 46  lame             Lamé constant
c     e(6), S(6), cmat(6,6)            strains, stresses, tangent matrix
c     glu(3,3)                         base vectors of E in columns
c
c
c...  storage of epsilon and sigma
c     11  22  33  12  13  23
c
c...  storage of elasticitymatrix
c     1111  1122  1133  1112  1131  1123  
c     2211  2222  2233  2212  2231  2223  
c     3311  3322  3333  3312  3331  3323  
c     1211  1222  1233  1212  1231  1223  
c     3111  3122  3133  3112  3131  3123 
c     2311  2322  2333  2312  2331  2323 
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      double precision L1(3,3),L2(3,3)
      dimension d(*),S(*),cmat(6,6),glu(3,3),e(6)
      dimension G(3,3),C(3,3),C_n(3,3),G_n(3,3),
     +   xmue(3),alpha(3),xlam(3),xN(3,3),fv1(3),fv2(3),sig(3),
     +   T1(6,3),T2(6,3),T1L1(6,3),T2L2(6,3),T1L1T1(6,6),T2L2T2(6,6)
c    
c      include 'eldata.h'
c      include 'iofile.h'
c
      xmue(1)  = d(40)
      xmue(2)  = d(41)
      xmue(3)  = d(42) 
      alpha(1) = d(43)  
      alpha(2) = d(44)  
      alpha(3) = d(45)
      xlame    = d(46)    
c
      do i=1,3
         do j=1,3
            G(i,j) = 0.d0
            do k=1,3
               G(i,j) = G(i,j) + glu(k,i) * glu(k,j)
            end do
         end do
      end do
c
c.....Right Cauchy Green strain tensor
      do i=1,3
         C(i,i) = 2.d0*e(i) + G(i,i)
      end do
      C(1,2) = e(4) + G(1,2)
      C(1,3) = e(5) + G(1,3)
      C(2,3) = e(6) + G(2,3)
      C(2,1) = e(4) + G(2,1)
      C(3,1) = e(5) + G(3,1)
      C(3,2) = e(6) + G(3,2)
c
c.....covariant coefficients of the Right--Cauchy--Green tensor C
c
      detC  =  C(1,1)* (C(2,2)*C(3,3)-C(2,3)*C(3,2))
     +        -C(1,2)* (C(2,1)*C(3,3)-C(2,3)*C(3,1))  
     +        +C(1,3)* (C(2,1)*C(3,2)-C(2,2)*C(3,1))  
      if (detC.lt.1.e-14) then 
         write(*,*) 'In Elas det C <= 0 '
         write(iow,*) 'In Elas det C <= 0 '
      end if
      do i=1,3
         do j=1,3
            if(dabs(C(i,j)).lt.1.e-13) C(i,j) = 0.d0
            if(dabs(G(i,j)).lt.1.e-13) G(i,j) = 0.d0
            C_n(i,j) = C(i,j)
            G_n(i,j) = G(i,j)
         end do
      end do
c
c.....Eigenvalue problem [C - lam_tr Cp] N = 0
      call pzero(xlam,3)
      call pzero(xN,3*3)
      call pzero(fv1,3)
      call pzero(fv2,3)
      ierr = 0
      call rsg(3,3,C,G,xlam,3,xN,fv1,fv2,ierr)
      do i=1,3
         if (xlam(i).lt.1.e-14) then
            write(*,*) 'xlam(',i,')=',xlam(i)
         end if
      end do
      if (abs(xlam(1)-xlam(2)).lt.1.e-15) xlam(2)=xlam(1)+1.e-15
      if (abs(xlam(1)-xlam(3)).lt.1.e-15) xlam(3)=xlam(1)+1.e-15
      if (abs(xlam(2)-xlam(3)).lt.1.e-15) xlam(3)=xlam(2)+1.e-15
c
c.....test the eigenvectors
c     xN * (Cp xN) = 1
      do k=1,3
         xx = 0.d0
         do i=1,3
            xx = xx + xN(i,k)*
     +         (G_n(i,1)*xN(1,k)+G_n(i,2)*xN(2,k)+G_n(i,3)*xN(3,k))
         end do
         if (dabs(xx)-1.d0.gt.1.e-14) then 
            write(*,*) 'N(',k,') G N(',k,') <>1'
            stop
         end if
      end do
c
c.....principal stresses
c     
      xJ       = xlam(1)*xlam(2)*xlam(3)  
      do i=1,3
         sig(i)=0.d0
         do k=1,3
          sig(i) = sig(i)+
     +     1.d0/xlam(i)*(xmue(k)*(xlam(i)**(alpha(k)/(2.d0))-1.d0))
         end do
       sig(i) = sig(i) + xlame/(2.d0*xlam(i))*(xJ-1.d0)
      end do
c
c.....transformation of the stresses      
      do i=1,3
         T1(1,i) = xN(1,i)*xN(1,i)
         T1(2,i) = xN(2,i)*xN(2,i)
         T1(3,i) = xN(3,i)*xN(3,i)
         T1(4,i) = xN(1,i)*xN(2,i)
         T1(5,i) = xN(1,i)*xN(3,i)
         T1(6,i) = xN(2,i)*xN(3,i)
      end do
      do i = 1,6
         S(i) = 0.d0
         do j=1,3  
            S(i) = S(i) + T1(i,j) * sig(j)
         end do
      end do

cfg      if(itp.eq.0)return

c
c.....transformation of the elastic tangent
      T2(1,1) = xN(1,1)*xN(1,2) + xN(1,2)*xN(1,1)
      T2(2,1) = xN(2,1)*xN(2,2) + xN(2,2)*xN(2,1)
      T2(3,1) = xN(3,1)*xN(3,2) + xN(3,2)*xN(3,1)
      T2(4,1) = xN(1,1)*xN(2,2) + xN(1,2)*xN(2,1)
      T2(5,1) = xN(1,1)*xN(3,2) + xN(1,2)*xN(3,1)
      T2(6,1) = xN(2,1)*xN(3,2) + xN(2,2)*xN(3,1)
      T2(1,2) = xN(1,1)*xN(1,3) + xN(1,3)*xN(1,1)
      T2(2,2) = xN(2,1)*xN(2,3) + xN(2,3)*xN(2,1)
      T2(3,2) = xN(3,1)*xN(3,3) + xN(3,3)*xN(3,1)
      T2(4,2) = xN(1,1)*xN(2,3) + xN(1,3)*xN(2,1)
      T2(5,2) = xN(1,1)*xN(3,3) + xN(1,3)*xN(3,1)
      T2(6,2) = xN(2,1)*xN(3,3) + xN(2,3)*xN(3,1)
      T2(1,3) = xN(1,2)*xN(1,3) + xN(1,3)*xN(1,2)
      T2(2,3) = xN(2,2)*xN(2,3) + xN(2,3)*xN(2,2)
      T2(3,3) = xN(3,2)*xN(3,3) + xN(3,3)*xN(3,2)
      T2(4,3) = xN(1,2)*xN(2,3) + xN(1,3)*xN(2,2)
      T2(5,3) = xN(1,2)*xN(3,3) + xN(1,3)*xN(3,2)
      T2(6,3) = xN(2,2)*xN(3,3) + xN(2,3)*xN(3,2)
c
      do i=1,3
         L1(i,i) = 0.d0
         do k=1,3
          L1(i,i) = L1(i,i)+
     +              xmue(k)/(xlam(i)**2.d0)*
     *              ((alpha(k)-2.d0)*xlam(i)**(alpha(k)/2.d0)+2.d0) 
         end do
       L1(i,i) = L1(i,i) + xlame/(xlam(i)**2.d0)  
      end do
      L1(1,2) = xlame*xlam(3)
      L1(1,3) = xlame*xlam(2)
      L1(2,3) = xlame*xlam(1)
      L1(2,1) = xlame*xlam(3)
      L1(3,1) = xlame*xlam(2) 
      L1(3,2) = xlame*xlam(1)
c
      call pzero(L2,3*3)
      L2(1,1) = (Sig(1)-Sig(2))/(xlam(1)-xlam(2))
      L2(2,2) = (Sig(1)-Sig(3))/(xlam(1)-xlam(3)) 
      L2(3,3) = (Sig(2)-Sig(3))/(xlam(2)-xlam(3)) 
c
      do i=1,6
         do j=1,3
            T1L1(i,j) = 0.d0
            T2L2(i,j) = 0.d0
            do k=1,3
            T1L1(i,j) = T1L1(i,j) + T1(i,k)*L1(k,j)
            T2L2(i,j) = T2L2(i,j) + T2(i,k)*L2(k,j)
            end do
         end do
      end do
      do i=1,6
         do j=1,6
            T1L1T1(i,j) = 0.d0
            T2L2T2(i,j) = 0.d0
            do k=1,3
               T1L1T1(i,j) = T1L1T1(i,j) + T1L1(i,k)*T1(j,k)
               T2L2T2(i,j) = T2L2T2(i,j) + T2L2(i,k)*T2(j,k)
            end do
            cmat(i,j) = T1L1T1(i,j) + T2L2T2(i,j)
         end do
      end do
c
c.....von Mises stresses 
c
        S(7) = dsqrt(0.5d0*((sig(1)-sig(2))*(sig(1)-sig(2))
     1                     +(sig(1)-sig(3))*(sig(1)-sig(3))
     2                     +(sig(2)-sig(3))*(sig(2)-sig(3))))
c
      return
      end
c      
      subroutine matread63(d)
      USE iofile
c-----------------------------------------------------------------------
c.... input of material data
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      include 'iofile.h'
      dimension d(*)
      imat = d(15)      
      
      goto(1,2,3) imat
c
c.... small  strain isotropic J2 plasticity     
1     if(ior.lt.0) write(*,1000)
1000  format('Input mat. parameter:E,nu,Y_0,Y_inf,H_lin,H_exp >',$)
      call dinput(d(40),6)
c
c.... write properties
      if(ior.lt.0) write(*  ,1001) (d(i),i=40,45)
                   write(iow,1001) (d(i),i=40,45) 
1001  format(/,5x,'Small  strain isotropic J2 plasticity',/,
     + 5x,'Youngs modulus E ......................... ',g12.4,/,
     + 5x,'Poissons ratio nu ........................ ',g12.4,/,
     + 5x,'initial Yield stress Y_0 (0 = elastisch).. ',g12.4,/,
     + 5x,'final   Yield stress Y_inf................ ',g12.4,/,
     + 5x,'lin. hardening parameter H_lin............ ',g12.4,/,
     + 5x,'exp. hardening parameter H_exp............ ',g12.4)
      return
c      
c.... finite strain isotropic J2 plasticity     
2     if(ior.lt.0) write(*,1002)
1002  format(' Input mat. parameter:E,nu,Y_0,Y_inf,H_lin,H_exp >',$)
      call dinput(d(40),6)
c
c.... write properties
      if(ior.lt.0) write(*  ,1003) (d(i),i=40,45)
                   write(iow,1003) (d(i),i=40,45) 
1003  format(/,5x,'Finite strain J2 isotropic plasticity',/,
     + 5x,'Youngs modulus E ............... ',g12.4,/,
     + 5x,'Poissons ratio nu .............. ',g12.4,/,
     + 5x,'initial Yield stress Y_0 ....... ',g12.4,/,
     + 5x,'final   Yield stress Y_inf...... ',g12.4,/,
     + 5x,'lin. hardening parameter H_lin.. ',g12.4,/,
     + 5x,'exp. hardening parameter H_exp.. ',g12.4)
      return
c      
c.... Ogden`s material law
3     if(ior.lt.0) write(*,1004)
1004  format(' Input mat. parameter: '
     + 'mue_1, mue_2, mue_3, alpha_1, alpha_2, alpha_3, lambda >',$)
      call dinput(d(40),7)
c
c.... write properties
      if(ior.lt.0) write(*  ,1005) (d(i),i=40,46)
                   write(iow,1005) (d(i),i=40,46) 
1005  format(/,5x,'Ogden`s material law',/,
     + 5x,'Shear modulus mue_1............. ',g12.4,/,
     + 5x,'Shear modulus mue_2............. ',g12.4,/,
     + 5x,'Shear modulus mue_3............. ',g12.4,/,
     + 5x,'Exponent alpha_1 ............... ',g12.4,/,
     + 5x,'Exponent alpha_2 ............... ',g12.4,/,
     + 5x,'Exponent alpha_3 ............... ',g12.4,/,
     + 5x,'Lame constant................... ',g12.4,/,
     + 5x,'!!! sum(mue_i*alpha_i) = 2mue !!!') 
c
      return      
      end

      subroutine err63(difference,d_dx,d_dy,da,val,val_dx,val_dy)
      USE errin1
c-----------------------------------------------------------------------
c.... input of material data
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      include 'errin1.h'
      e_om(1) = e_om(1) + (difference**2+d_dx**2+d_dy**2)*da
      e_om(2) = e_om(2) + difference**2*da
      e_om(3) = e_om(3) + da
      e_om(4) = e_om(4) + (val**2+val_dx**2+val_dy**2)*da
      return
      end
      subroutine energeticnorm63(error)
      USE errin1
c-----------------------------------------------------------------------
c.... input of material data
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      include 'errin1.h'
      e_om(1) = e_om(1) + error
      return
      end
c-----------------------------------------------------------------------      
      subroutine testtripod63vnorm(exact,rin,nen,shp,t0,da,isw,n,l)
      USE errin1
      USE iofile
c----------------------------------------------------------------------- 
c     tests loaded tripod versus tripod constructed in Gauss point
c     test loaded derivatives of tripod with interpolated tripods
c-----------------------------------------------------------------------       
      implicit double precision (a-h,o-z)
      real*8 r0(3,3),rin(3,3,nen),a(3,3,3),shp(3,nen),a3(3),t0(3,3)
      real*8 vec1(3),vec2(3),vec3(3),da,pi,exact(3,3,3)
      integer isw
c      include 'iofile.h'
c      include 'errin1.h'
      
      pi = dacos(0.0d0)*2
      
      do i = 1,3
       do j = 1,3
         do j2 = 1,3
           a(i,j,j2)  = 0.d0
           do k = 1,nen
             a(i,j,j2)  = a(i,j,j2)  + rin(i,j2,k) * shp(j,k)
           end do
         end do  
        end do 
      end do
c      if (isw.eq.3.and.l.eq.1) write(iow,*) n
      !write(iow,*) 'GP: ',l
      !write(iow,*) a(1:3,1,1)-exact(1:3,1,1)
      !write(iow,*) a(1:3,2,1)-exact(1:3,2,1)
      !write(iow,*) a(1:3,1,2)-exact(1:3,1,2)
      !write(iow,*) a(1:3,2,2)-exact(1:3,2,2)
      !write(iow,*) a(1:3,1,3)-exact(1:3,1,3)
      !write(iow,*) a(1:3,2,3)-exact(1:3,2,3)
      !write(iow,*) 'GP: ',l
      !write(iow,*) exact(1:3,1,1)
      !write(iow,*) exact(1:3,2,1)
      !write(iow,*) exact(1:3,1,2)
      !write(iow,*) exact(1:3,2,2)
      !write(iow,*) exact(1:3,1,3)
      !write(iow,*) exact(1:3,2,3)
      e_om(3) = e_om(3) + da

      do i=1,3
        r0len = sqrt(exact(1,3,i)**2+exact(2,3,i)**2+exact(3,3,i)**2)
        t0len = sqrt(t0(1,i)**2+t0(2,i)**2+t0(3,i)**2)
        dot = vskal(exact(1,3,i),t0(1,i))
        angle = dot/(r0len*t0len)
c        if (isw.eq.3.or.isw.eq.8) write(iow,*) r0len,t0len,angle
        if (isw.eq.9) then
          e_om(7+2*(i-1))=e_om(7+2*(i-1))+(angle-1.d0)**2*da
          e_om(8+2*(i-1))=e_om(8+2*(i-1))+(1.0d0-r0len)**2*da
        end if
      end do  
      return
      end
      
      subroutine testtripod63vnormIn(rin,nen,shp,t0,da,isw,n,l)
      USE errin1
      USE iofile
c----------------------------------------------------------------------- 
c     tests loaded tripod versus tripod constructed in Gauss point
c     test loaded derivatives of tripod with interpolated tripods
c-----------------------------------------------------------------------       
      implicit double precision (a-h,o-z)
      real*8 r0(3,3),rin(3,3,nen),a(3,3,3),shp(3,nen),a3(3),t0(3,3)
      real*8 vec1(3),vec2(3),vec3(3),da,pi
      integer isw
c      include 'iofile.h'
c      include 'errin1.h'
      
      pi = dacos(0.0d0)*2
      
      do i = 1,3
       do j = 1,3
         do j2 = 1,3
           a(i,j,j2)  = 0.d0
           do k = 1,nen
             a(i,j,j2)  = a(i,j,j2)  + rin(i,j2,k) * shp(j,k)
           end do
         end do  
        end do 
      end do
c      if (isw.eq.3.and.l.eq.1) write(iow,*) n
      !write(iow,*) 'GP: ',l
      !write(iow,*) a(1:3,1,1)-exact(1:3,1,1)
      !write(iow,*) a(1:3,2,1)-exact(1:3,2,1)
      !write(iow,*) a(1:3,1,2)-exact(1:3,1,2)
      !write(iow,*) a(1:3,2,2)-exact(1:3,2,2)
      !write(iow,*) a(1:3,1,3)-exact(1:3,1,3)
      !write(iow,*) a(1:3,2,3)-exact(1:3,2,3)
      !write(iow,*) 'GP: ',l
      !write(iow,*) exact(1:3,1,1)
      !write(iow,*) exact(1:3,2,1)
      !write(iow,*) exact(1:3,1,2)
      !write(iow,*) exact(1:3,2,2)
      !write(iow,*) exact(1:3,1,3)
      !write(iow,*) exact(1:3,2,3)
      e_om(3) = e_om(3) + da

      do i=1,3
        r0len = sqrt(a(1,3,i)**2+a(2,3,i)**2+a(3,3,i)**2)
        t0len = sqrt(t0(1,i)**2+t0(2,i)**2+t0(3,i)**2)
        dot = vskal(a(1,3,i),t0(1,i))
        angle = dot/(r0len*t0len)
c        if (isw.eq.3.or.isw.eq.8) write(iow,*) r0len,t0len,angle
        if (isw.eq.9) then
          e_om(7+2*(i-1))=e_om(7+2*(i-1))+(angle-1.d0)**2*da
          e_om(8+2*(i-1))=e_om(8+2*(i-1))+(1.0d0-r0len)**2*da
        end if
      end do  
      return
      end
      

      subroutine testdirector63(r0,t0,da,isw)
      USE errin1
      USE iofile
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      real*8 r0(3,3),t0(3,3), r0len,t0len,dotp
      real*8 angle,pi,da,const
      integer isw
c      include 'iofile.h'
c      include 'errin1.h'
      
      const = 1.d0-1d-16
      pi = dacos(0.0d0)*2
      r0len = sqrt(r0(1,3)**2+r0(2,3)**2+r0(3,3)**2)
      t0len = sqrt(t0(1,3)**2+t0(2,3)**2+t0(3,3)**2)
      dotp = vskal(r0(1,3),t0(1,3))
      angle = dotp/(r0len*t0len)
c     catches numerical errors!!!      
      if (angle.ge.1.d0) then
       angle=1.0d0
c      else if (abs(angle1)-const.ge.0) then
c       angle1=1.d0
      end if
      
      angle = dacos(angle)
      
      angle = angle*180.d0/pi

      if (isw.eq.9) then
        e_om(5) = e_om(5) + angle**2*da
        e_om(6) = e_om(6) + (1.0d0-r0len)**2*da
        e_om(3) = e_om(3) + da
      end if
      return
      end
      
      subroutine testdirector63vnorm(r0,t0,da,isw)
      USE errin1
      USE iofile
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      real*8 r0(3,3),t0(3,3), r0len,t0len,dotp
      real*8 angle,pi,da,const
      integer isw
c      include 'iofile.h'
c      include 'errin1.h'
      
      const = 1.d0-1d-16
      pi = dacos(0.0d0)*2
      r0len = sqrt(r0(1,3)**2+r0(2,3)**2+r0(3,3)**2)
      t0len = sqrt(t0(1,3)**2+t0(2,3)**2+t0(3,3)**2)
      dotp = vskal(r0(1,3),t0(1,3))
      angle = dotp/(r0len*t0len)

      if (isw.eq.9) then
        e_om(5) = e_om(5) + (angle-1.d0)**2*da
        e_om(6) = e_om(6) + (1.0d0-r0len)**2*da
        e_om(3) = e_om(3) + da
      end if
      return
      end

      subroutine transrot63(inode,shp,Tgp,AnGP,rwn,nen,iuseTgp,TA)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension shp(3,nen),Tgp(3,3),AnGP(2,2,nen),rwn(3,3,nen)
      dimension a1(3),a2(3),a3(3),g(3,3),TA(3,3),Agp(2,3)
      if (iuseTgp.eq.1) then
        do i= 1,3
          a1(i) = Tgp(i,2)
          a2(i) = -Tgp(i,1)
          Tgp(i,3) = 0
        end do
        do i = 1,3
          Agp(1,i) = dot(a1(1),rwn(1,i,inode),3)
          Agp(2,i) = dot(a2(1),rwn(1,i,inode),3)
        end do
        
      else if (iuseTgp.eq.2) then
!        do i= 1,3
!          a1(i) = Tgp(i,2)
!          a2(i) =  -Tgp(i,1)
!          a3(i) = Tgp(i,3)
!        end do
!        do i = 1,2
!          Agp(1,i) =  dot(a1(1),rwn(1,i,inode),3)
!     -             -( dot(a3(1),rwn(1,i,inode),3)
!     *               *dot(a1(1),rwn(1,3,inode),3)
!     /               /dot(a3(1),rwn(1,3,inode),3))
!          Agp(2,i) =  dot(a2(1),rwn(1,i,inode),3)
!     -             -( dot(a3(1),rwn(1,i,inode),3)
!     *             *  dot(a2(1),rwn(1,3,inode),3)
!     /               /dot(a3(1),rwn(1,3,inode),3))
!        end do
        do i =1,2
          do j = 1,2
            Agp(i,j) = AnGP(i,j,inode)
          end do
        end do
        
        Agp(2,3)=0.0d0
        Agp(1,3)=0.0d0
        
      end if
      call matmulf(Tgp,Agp,3,2,3,TA)
      return
      end

c      subroutine computeMandMalpha63(M,M_1,M_2,wn,wn_1,wn_2,d,d_1,d_2,h)
      subroutine computeMandMalpha63(M,wn,wn_1,wn_2,d,d_1,d_2,h,h_1,h_2)
c-----------------------------------------------------------------------
c     matrix M(h) and M_,alpha(h) for geometrical matrix 
c-----------------------------------------------------------------------
c...  output variables
      double precision M(3,3),M_1(3,3),M_2(3,3)
c...  input variables
      double precision wn(3),wn_1(3),wn_2(3),d(3),d_1(3),d_2(3),h(3)
      double precision h_1(3),h_2(3)
c...  internal variables
      double precision w2,w,dotbw,dotb1w,dotbw1,dotb2w,dotbw2,w_1,w_2
      double precision dotb_1w,dotb_2w
      double precision c3,c10,c11,c3_1,c3_2,c10_1,c10_2,c11_1,c11_2
      double precision c10M,c11M,c10M1,c10M2
      double precision sinw,cosw,coswm1,coswm2
      double precision b(3),b_1(3),b_2(3),t(3),t_1(3),t_2(3),b1(3),b2(3)
      integer i,j
c...  declaration of data type for used defined functions
      double precision dot


c     determine ||omega||, its square and b_,alpha and b and b.omega
      w2  = dot(wn,wn,3)
      w   = dsqrt(w2)
      if (w.lt.1.d-10) then
        w_1 = 0.d0
        w_2 = 0.d0
      else
        w_1 = dot(wn_1,wn,3)/w
        w_2 = dot(wn_2,wn,3)/w
      end if
      call vecp(d,h,b)
      call vecp(d,h_1,b1)
      call vecp(d,h_2,b2)
      call vecp(d_1,h_1,b_1)
      call vecp(d_2,h_2,b_2)
      dotbw  = dot(b,wn,3)
      dotb1w = dot(b1,wn,3)
      dotb_1w= dot(b_1,wn,3)
      dotbw1 = dot(b1,wn_1,3)
      dotb2w = dot(b2,wn,3)
      dotb_2w= dot(b_2,wn,3)
      dotbw2 = dot(b2,wn_2,3)

c     if omega is small use Taylor expansion series for factors
      if (w.lt.1.d-2) then
        c3    =  1.d0/6.d0  * (1.d0 + w2/60.d0)
        c10   =  1.d0/6.d0  * (1.d0 + w2/30.d0)
        c11   = -1.d0/360.d0 * (1.d0 + w2/21.d0)
        c3_1  =  1.d0/180.d0 * w * w_1
        c3_2  =  1.d0/180.d0 * w * w_2
        c10_1 =  1.d0/90.d0 * w * w_1
        c10_2 =  1.d0/90.d0 * w * w_2
        c11_1 =  1.d0/3780.d0 * w * w_1
        c11_2 =  1.d0/3780.d0 * w * w_2
      else
c       use exact formulas for factors
        sinw  = dsin(w)
        cosw  = dcos(w)
        coswm1 = dcos(w) - 1.d0
        coswm2 = coswm1*coswm1
        
c       factors for derived M
        c11 = (((-3.d0*sinw+2.d0*w+w*cosw)*w*coswm1)-
     -         ((4.d0*coswm1+w2+w*sinw)*(4.d0*coswm1-w*sinw)))/
     /        (2.d0*w2*w2*w*coswm2)
        c11_1 = c11 * w_1
        c11_2 = c11 * w_2
c        
        c10 = (w*coswm2-(sinw-w)*(coswm1-w*sinw))/
     /        (2.d0*w2*coswm2)
        c10_1 = c10 * w_1
        c10_2 = c10 * w_2
c
        c3 = ((w*cosw-sinw)*w*coswm1-(w*sinw+2.d0*coswm1)*
     *        (2.d0*coswm1-w*sinw))/(w2*w*coswm2)
        c3_1 = c3 * w_1
        c3_2 = c3 * w_2
c
c       factors for M
        c3   = (w*sinw+2.d0*coswm1)/(w2*coswm1)
        c10  = (sinw - w )/( 2.d0 *w *coswm1 )
        c11  = (w2 + w*sinw + 4.d0*coswm1)/( 2.d0*w2*w2*coswm1)
      end if
c     compute M
      c10M  = c10 * dotbw  - dot(d,h,3)
      c11M  = c11 * dotbw  

      do i = 1,3
        t(i) =  c11M*wn(i) - c3*b(i) 
        do j = 1,i
          M(i,j) = 0.5d0*(d(i)*h(j) + h(i)*d(j)) 
     1           + 0.5d0*(t(i)*wn(j) + wn(i)*t(j)) 
          M(j,i) = M(i,j)
        end do
        M(i,i) = M(i,i) + c10M
      end do
c
c     compute derived M
      c10M1 = c10_1*dotb1w + c10*dotb_1w + c10*dotbw1 - dot(d_1,h_1,3)
      c10M2 = c10_2*dotb2w + c10*dotb_2w + c10*dotbw2 - dot(d_2,h_2,3)

      do i = 1,3
        t_1(i)= -c3_1*b1(i) -c3*b_1(i) +c11_1*dotb1w*wn(i) 
     +          +c11*dotb_1w*wn(i) +c11*dotbw1*wn(i) +c11*dotb1w*wn_1(i)
        t_2(i)= -c3_2*b2(i) -c3*b_2(i) +c11_2*dotb2w*wn(i) 
     +          +c11*dotb_2w*wn(i) +c11*dotbw2*wn(i) +c11*dotb2w*wn_2(i)
        do j =1,i
          M_1(i,j) = 0.5d0*(d_1(i)*h_1(j) + h_1(i)*d_1(j))
     +             + 0.5d0*( t_1(i)*wn(j) + wn(i)*t_1(j)
     +                      +t(i)*wn_1(j) + wn_1(i)*t(j))
          M_2(i,j) = 0.5d0*(d_2(i)*h_2(j) + h_2(i)*d_2(j))
     +             + 0.5d0*( t_2(i)*wn(j) + wn(i)*t_2(j)
     +                      +t(i)*wn_2(j) + wn_2(i)*t(j))    
          M_1(j,i) = M_1(i,j)
          M_2(j,i) = M_2(i,j)
        end do
        M_1(i,i) = M_1(i,i) + c10M1
        M_2(i,i) = M_2(i,i) + c10M2
      end do
      do i = 1,3
        do j = 1,3
          M(i,j) = M(i,j) + M_1(i,j) + M_2(i,j)
        end do
      end do
      return
      end

      subroutine computeGww63(Gww,wn,wn_1,wn_2,d,d_1,d_2,h,h_1,h_2,
     1                        h1,h2,twn,twn_1,twn_2,inode,jnode,shp,
     2                        Ilen,Jlen)
c-----------------------------------------------------------------------
c     matrix M(h) and M_,alpha(h) for geometrical matrix 
c     and combination to G^ßß
c-----------------------------------------------------------------------
c...  output variables
      double precision Gww(3,3)
c...  input variables
      double precision wn(3),wn_1(3),wn_2(3),d(3),d_1(3),d_2(3),h(3)
      double precision h_1(3),h_2(3),h1(3),h2(3)
      double precision twn(3,3,*),twn_1(3,3,*),twn_2(3,3,*),shp(3,*)
      integer inode,jnode,Ilen,Jlen
c...  internal variables
      double precision M(3,3),M_1(3,3),M_2(3,3),Mh1(3,3),Mh2(3,3)
      double precision w2,w,dotbw,dotb1w,dotbw1,dotb2w,dotbw2,w_1,w_2
      double precision dotb_1w,dotb_2w
      double precision c3,c10,c11,c3_1,c3_2,c10_1,c10_2,c11_1,c11_2
      double precision c10M,c11M,c10M1,c10M2
      double precision sinw,cosw,coswm1,coswm2
      double precision b(3),b_1(3),b_2(3),t(3),t_1(3),t_2(3),b1(3),b2(3)
      double precision t1(3),t2(3)
      double precision temp22(Ilen,Jlen),temp23(Ilen,3),temp33(3,3)
      integer i,j
c...  declaration of data type for used defined functions
      double precision dot


c     determine ||omega||, its square and b_,alpha and b and b.omega
      w2  = dot(wn,wn,3)
      w   = dsqrt(w2)
      if (w.lt.1.d-10) then
        w_1 = 0.d0
        w_2 = 0.d0
      else
        w_1 = dot(wn_1,wn,3)/w
        w_2 = dot(wn_2,wn,3)/w
      end if
      call vecp(d,h,b)
      call vecp(d,h_1,b1)
      call vecp(d,h_2,b2)
      call vecp(d_1,h_1,b_1)
      call vecp(d_2,h_2,b_2)
      dotbw  = dot(b,wn,3)
      dotb1w = dot(b1,wn,3)
      dotb_1w= dot(b_1,wn,3)
      dotbw1 = dot(b1,wn_1,3)
      dotb2w = dot(b2,wn,3)
      dotb_2w= dot(b_2,wn,3)
      dotbw2 = dot(b2,wn_2,3)

c     if omega is small use Taylor expansion series for factors
      if (w.lt.1.d-2) then
        c3    =  1.d0/6.d0  * (1.d0 + w2/60.d0)
        c10   =  1.d0/6.d0  * (1.d0 + w2/30.d0)
        c11   = -1.d0/360.d0 * (1.d0 + w2/21.d0)
        c3_1  =  1.d0/180.d0 * w * w_1
        c3_2  =  1.d0/180.d0 * w * w_2
        c10_1 =  1.d0/90.d0 * w * w_1
        c10_2 =  1.d0/90.d0 * w * w_2
        c11_1 =  1.d0/3780.d0 * w * w_1
        c11_2 =  1.d0/3780.d0 * w * w_2
      else
c       use exact formulas for factors
        sinw  = dsin(w)
        cosw  = dcos(w)
        coswm1 = dcos(w) - 1.d0
        coswm2 = coswm1*coswm1
        
c       factors for derived M
        c11 = (((-3.d0*sinw+2.d0*w+w*cosw)*w*coswm1)-
     -         ((4.d0*coswm1+w2+w*sinw)*(4.d0*coswm1-w*sinw)))/
     /        (2.d0*w2*w2*w*coswm2)
        c11_1 = c11 * w_1
        c11_2 = c11 * w_2
c        
        c10 = (w*coswm2-(sinw-w)*(coswm1-w*sinw))/
     /        (2.d0*w2*coswm2)
        c10_1 = c10 * w_1
        c10_2 = c10 * w_2
c
        c3 = ((w*cosw-sinw)*w*coswm1-(w*sinw+2.d0*coswm1)*
     *        (2.d0*coswm1-w*sinw))/(w2*w*coswm2)
        c3_1 = c3 * w_1
        c3_2 = c3 * w_2
c
c       factors for M
        c3   = (w*sinw+2.d0*coswm1)/(w2*coswm1)
        c10  = (sinw - w )/( 2.d0 *w *coswm1 )
        c11  = (w2 + w*sinw + 4.d0*coswm1)/( 2.d0*w2*w2*coswm1)
      end if
c
c...  compute M(h)
      c10M  = c10 * dotbw  - dot(d,h,3)
      c11M  = c11 * dotbw  

      do i = 1,3
        t(i) =  c11M*wn(i) - c3*b(i) 
        do j = 1,i
          M(i,j) = 0.5d0*(d(i)*h(j) + h(i)*d(j)) 
     1           + 0.5d0*(t(i)*wn(j) + wn(i)*t(j)) 
          M(j,i) = M(i,j)
        end do
        M(i,i) = M(i,i) + c10M
      end do
c
c...  compute M(h1)
      call vecp(d,h1,b)
      dotbw  = dot(b,wn,3)
      c10M  = c10 * dotbw  - dot(d,h1,3)
      c11M  = c11 * dotbw  

      do i = 1,3
        t(i) =  c11M*wn(i) - c3*b(i) 
        do j = 1,i
          Mh1(i,j) = 0.5d0*(d(i)*h1(j) + h1(i)*d(j)) 
     1           + 0.5d0*(t(i)*wn(j) + wn(i)*t(j)) 
          Mh1(j,i) = Mh1(i,j)
        end do
        Mh1(i,i) = Mh1(i,i) + c10M
      end do
c
c...  compute M(h2)
      call vecp(d,h2,b)
      dotbw  = dot(b,wn,3)
      c10M  = c10 * dotbw  - dot(d,h2,3)
      c11M  = c11 * dotbw  

      do i = 1,3
        t(i) =  c11M*wn(i) - c3*b(i) 
        do j = 1,i
          Mh2(i,j) = 0.5d0*(d(i)*h2(j) + h2(i)*d(j)) 
     1           + 0.5d0*(t(i)*wn(j) + wn(i)*t(j)) 
          Mh2(j,i) = Mh2(i,j)
        end do
        Mh2(i,i) = Mh2(i,i) + c10M
      end do
c
c     compute derived M: M_,1(h_,1)
      c10M1 = c10_1*dotb1w + c10*dotb_1w + c10*dotbw1 - dot(d_1,h_1,3)
      c10M2 = c10_2*dotb2w + c10*dotb_2w + c10*dotbw2 - dot(d_2,h_2,3)

      do i = 1,3
        c11M  = c11 * dotb1w 
        t1(i) =  c11M*wn(i) - c3*b1(i)
        c11M  = c11 * dotb2w  
        t2(i) =  c11M*wn(i) - c3*b2(i) 
        t_1(i)= -c3_1*b1(i) -c3*b_1(i) +c11_1*dotb1w*wn(i) 
     +          +c11*dotb_1w*wn(i) +c11*dotbw1*wn(i) +c11*dotb1w*wn_1(i)
        t_2(i)= -c3_2*b2(i) -c3*b_2(i) +c11_2*dotb2w*wn(i) 
     +          +c11*dotb_2w*wn(i) +c11*dotbw2*wn(i) +c11*dotb2w*wn_2(i)
        do j =1,i
          M_1(i,j) = 0.5d0*(d_1(i)*h_1(j) + h_1(i)*d_1(j))
     +             + 0.5d0*( t_1(i)*wn(j) + wn(i)*t_1(j)
     +                      +t1(i)*wn_1(j) + wn_1(i)*t1(j))
          M_2(i,j) = 0.5d0*(d_2(i)*h_2(j) + h_2(i)*d_2(j))
     +             + 0.5d0*( t_2(i)*wn(j) + wn(i)*t_2(j)
     +                      +t2(i)*wn_2(j) + wn_2(i)*t2(j))    
          M_1(j,i) = M_1(i,j)
          M_2(j,i) = M_2(i,j)
        end do
        M_1(i,i) = M_1(i,i) + c10M1
        M_2(i,i) = M_2(i,i) + c10M2
      end do
      
c     
      Gww=0.d0
      call mttmul (twn_1(1,1,inode),Mh1,Ilen,3,3,temp23)
      call matmulf(temp23,twn(1,1,jnode),Ilen,3,Jlen,temp22)
      do i = 1,Ilen
        do j = 1, Jlen
          Gww(i,j) = Gww(i,j) + temp22(i,j)*shp(3,jnode)
        end do
      end do
c
      call mttmul (twn(1,1,inode),Mh1,Ilen,3,3,temp23)
      call matmulf(temp23,twn_1(1,1,jnode),Ilen,3,Jlen,temp22)
      do i = 1,Ilen
        do j = 1, Jlen
          Gww(i,j) = Gww(i,j) + temp22(i,j)*shp(3,inode)
        end do
      end do
c      
      call mttmul (twn_2(1,1,inode),Mh2,Ilen,3,3,temp23)
      call matmulf(temp23,twn(1,1,jnode),Ilen,3,Jlen,temp22)
      do i = 1,Ilen
        do j = 1, Jlen
          Gww(i,j) = Gww(i,j) + temp22(i,j)*shp(3,jnode)
        end do
      end do
c
      call mttmul (twn(1,1,inode),Mh2,Ilen,3,3,temp23)
      call matmulf(temp23,twn_2(1,1,jnode),Ilen,3,Jlen,temp22)
      do i = 1,Ilen
        do j = 1, Jlen
          Gww(i,j) = Gww(i,j) + temp22(i,j)*shp(3,inode)
        end do
      end do
c
      call matadd (M_1,M_2,3,3,temp33)
      call matadd (M,temp33,3,3,M)
      call mttmul (twn(1,1,inode),M,Ilen,3,3,temp23)
      call matmulf(temp23,twn(1,1,jnode),Ilen,3,Jlen,temp22)
      do i = 1,Ilen
        do j = 1, Jlen
          Gww(i,j) = Gww(i,j) + temp22(i,j)
        end do
      end do
      return
      end
      
      subroutine computeGww63_Mult_inputvalues(sig,xu,Mhq,Mh1,Mh2,da,d)
c-----------------------------------------------------------------------
c     matrix M(h) and M_,alpha(h) for geometrical matrix 
c     and combination to G^ßß
c-----------------------------------------------------------------------
c...  output variables
      double precision Mhq(3,3),Mh1(3,3),Mh2(3,3)
c...  input variables
      double precision sig(*),xu(3,*),da,d(3,3)
c...  internal variables
      double precision hq(3),h1(3),h2(3)
      integer i,j
c...  declaration of data type for used defined functions
      double precision dot
      
c     Compute vectors for M matrices
      hq(1:3) = (sig(7)*xu(1:3,1) + sig(8)*xu(1:3,2))*da
      h1(1:3) = (sig(4)*xu(1:3,1) + sig(6)*xu(1:3,2))*da
      h2(1:3) = (sig(5)*xu(1:3,2) + sig(6)*xu(1:3,1))*da
c...  compute Mhq = M(hq) + M,1(h1) + M,2(h2)
      do i = 1,3
        do j = 1,i
          Mhq(i,j) = 0.5d0*(d(i,3)*hq(j) + hq(i)*d(j,3)
     1                    + d(i,1)*h1(j) + h1(i)*d(j,1)
     2                    + d(i,2)*h2(j) + h2(i)*d(j,2)) 
          Mhq(j,i) = Mhq(i,j)
        end do
        Mhq(i,i) = Mhq(i,i) - dot(d(1:3,3),hq,3)
     1         - dot(d(1:3,1),h1,3) - dot(d(1:3,2),h2,3)              
      end do
c
c...  compute M(h1)
      do i = 1,3
        do j = 1,i
          Mh1(i,j) = 0.5d0*(d(i,3)*h1(j) + h1(i)*d(j,3)) 
          Mh1(j,i) = Mh1(i,j)
        end do
        Mh1(i,i) = Mh1(i,i) - dot(d(1:3,3),h1,3)
      end do
c
c...  compute M(h2)
      do i = 1,3
        do j = 1,i
          Mh2(i,j) = 0.5d0*(d(i,3)*h2(j) + h2(i)*d(j,3)) 
          Mh2(j,i) = Mh2(i,j)
        end do
        Mh2(i,i) = Mh2(i,i) - dot(d(1:3,3),h2,3)
      end do
c
!c     compute derived M: M_,1(h_,1)
!      do i = 1,3
!        do j =1,i
!          M_1(i,j) = 0.5d0*(d_1(i)*h_1(j) + h_1(i)*d_1(j))
!          M_2(i,j) = 0.5d0*(d_2(i)*h_2(j) + h_2(i)*d_2(j))
!          M_1(j,i) = M_1(i,j)
!          M_2(j,i) = M_2(i,j)
!        end do
!        M_1(i,i) = M_1(i,i) - dot(d_1,h_1,3)
!        M_2(i,i) = M_2(i,i) - dot(d_2,h_2,3)
!      end do

      return
      end
      
      subroutine computeGww63_Om_Mult_fast(Gww,twn,inode,jnode,shp,
     2                        Ilen,Jlen,Mhq,Mh1,Mh2)
c-----------------------------------------------------------------------
c     matrix M(h) and M_,alpha(h) for geometrical matrix 
c     and combination to G^ßß
c     Fast version. M matrices are computed in Gauss loop for speed
c     This is the fastest version. 2,3 works and is checked, 2,2 seems to work also         
c-----------------------------------------------------------------------
c...  output variables
      double precision Gww(3,3)
c...  input variables
      double precision Mhq(3,3),Mh1(3,3),Mh2(3,3)
      double precision twn(3,3,*),shp(3,*)
      integer inode,jnode,Ilen,Jlen
c...  internal variables
      double precision temp23(Ilen,3),M(3,3)
      double precision fMhq,fMh1,fMh2
      integer i,j,K
 
!      Gww=0.d0   


c     NEW; FAST
      fMhq=shp(3,inode)*shp(3,jnode)
      fMh1=shp(3,inode)*shp(1,jnode)+shp(1,inode)*shp(3,jnode)
      fMh2=shp(3,inode)*shp(2,jnode)+shp(2,inode)*shp(3,jnode)
      M(1:3,1:3)=Mhq(1:3,1:3)*fMhq+Mh1(1:3,1:3)*fMh1+Mh2(1:3,1:3)*fMh2
!      call mttmul (twn(1,1,inode),M,Ilen,3,3,temp23)
c     mttmul inlined for speed      
      DO 100 J=1,3
      DO 200 I=1,Ilen
  200 temp23(I,J) = twn(1,I,inode) * M(1,J)
      DO 100 K=2,3
      DO 100 I=1,Ilen
  100 temp23(I,J) = temp23(I,J) + twn(K,I,inode) * M(K,J)
!      call matmulf(temp23,twn(1,1,jnode),Ilen,3,Jlen,Gww(1:Ilen,1:Jlen))
c     matmul inlined for speed            
      DO 101 J=1,Jlen
      DO 201 I=1,Ilen
  201 Gww(I,J) = temp23(I,1) * twn(1,J,jnode)
      DO 101 K=2,3
      DO 101 I=1,Ilen
  101 Gww(I,J) = Gww(I,J) + temp23(I,K) * twn(K,J,jnode)
c     END NEW      
      return
      end      
      
      subroutine computeGww63_Om_Mult_test(Gww,d,d_1,d_2,h,h_1,
     1                        h_2,h1,h2,twn,twn_1,twn_2,inode,jnode,shp,
     2                        Ilen,Jlen)
c-----------------------------------------------------------------------
c     matrix M(h) and M_,alpha(h) for geometrical matrix 
c     and combination to G^ßß
c     Quite fast version. M matrices are summed up before computing T_3I^T*M*T_3K
c     2,3 works and is checked, 2,2 seems to work also           
c-----------------------------------------------------------------------
c...  output variables
      double precision Gww(3,3)
c...  input variables
      double precision d(3),d_1(3),d_2(3),h(3)
      double precision h_1(3),h_2(3),h1(3),h2(3)
      double precision twn(3,3,*),twn_1(3,3,*),twn_2(3,3,*),shp(3,*)
      integer inode,jnode,Ilen,Jlen
c...  internal variables
      double precision M(3,3),M_1(3,3),M_2(3,3),Mh1(3,3),Mh2(3,3)
      double precision temp22(Ilen,Jlen),temp23(Ilen,3),temp33(3,3)
      double precision temp23_2(Ilen,3),temp32(3,Jlen),temp32_2(3,Jlen)
      double precision temp22_2(Ilen,Jlen),temp22_3(Ilen,Jlen)
      double precision facMh1l,facMh2l,facMh1r,facMh2r
      integer i,j
c...  declaration of data type for used defined functions
      double precision dot


c
c...  compute M(h)
      do i = 1,3
        do j = 1,i
          M(i,j) = 0.5d0*(d(i)*h(j) + h(i)*d(j)) 
          M(j,i) = M(i,j)
        end do
        M(i,i) = M(i,i) - dot(d,h,3)
      end do
c
c...  compute M(h1)
      do i = 1,3
        do j = 1,i
          Mh1(i,j) = 0.5d0*(d(i)*h1(j) + h1(i)*d(j)) 
          Mh1(j,i) = Mh1(i,j)
        end do
        Mh1(i,i) = Mh1(i,i) - dot(d,h1,3)
      end do
c
c...  compute M(h2)
      do i = 1,3
        do j = 1,i
          Mh2(i,j) = 0.5d0*(d(i)*h2(j) + h2(i)*d(j)) 
          Mh2(j,i) = Mh2(i,j)
        end do
        Mh2(i,i) = Mh2(i,i) - dot(d,h2,3)
      end do
c
c     compute derived M: M_,1(h_,1)
      do i = 1,3
        do j =1,i
          M_1(i,j) = 0.5d0*(d_1(i)*h_1(j) + h_1(i)*d_1(j))
          M_2(i,j) = 0.5d0*(d_2(i)*h_2(j) + h_2(i)*d_2(j))
          M_1(j,i) = M_1(i,j)
          M_2(j,i) = M_2(i,j)
        end do
        M_1(i,i) = M_1(i,i) - dot(d_1,h_1,3)
        M_2(i,i) = M_2(i,i) - dot(d_2,h_2,3)
      end do
      
c     
      Gww=0.d0   

c      Sowohl MKL als auch matmul ist langsamer als das selbst programmierte!!!
!      temp32 = matmul (Mh1,twn(1:3,1:Jlen,jnode))
!      call dgemm('t','n',Ilen,Jlen,3,shp(3,jnode),
!     + twn_1(1:3,1:Ilen,inode),3,temp32,3,1.0d0,Gww,3)
!      call mttaml63 (twn_1(1:3,1:Ilen,inode),temp32,Ilen,3,Jlen,
!     + Gww,shp(3,jnode))



c     NEW; FAST
      facMh1l=shp(1,inode)*shp(3,jnode)
      facMh2l=shp(2,inode)*shp(3,jnode)
      facMh1r=shp(3,inode)*shp(1,jnode)
      facMh2r=shp(3,inode)*shp(2,jnode)
      M(1:3,1:3) = M(1:3,1:3) + M_1(1:3,1:3) + M_2(1:3,1:3)
     +           + Mh1(1:3,1:3)*facMh1l+Mh2(1:3,1:3)*facMh2l
     +           + Mh1(1:3,1:3)*facMh1r+Mh2(1:3,1:3)*facMh2r
      call mttmul (twn(1,1,inode),M,Ilen,3,3,temp23)
      call matmulf(temp23,twn(1,1,jnode),Ilen,3,Jlen,Gww(1:Ilen,1:Jlen))
c     END NEW      
   

!c     START ALTES VORGEHEN
!      call mttmul (twn_1(1,1,inode),Mh1,Ilen,3,3,temp23)
!      call mttmul (twn_2(1,1,inode),Mh2,Ilen,3,3,temp23_2)
!      temp23(1:Ilen,1:3) = temp23(1:Ilen,1:3) + temp23_2(1:Ilen,1:3)
!      M(1:3,1:3) = M(1:3,1:3) + M_1(1:3,1:3) + M_2(1:3,1:3)
!      call mttmul (twn(1,1,inode),M,Ilen,3,3,temp23_2)
!      temp23(1:Ilen,1:3) = temp23(1:Ilen,1:3)*shp(3,jnode) 
!     +                   + temp23_2(1:Ilen,1:3)
!      call matmulf(temp23,twn(1,1,jnode),Ilen,3,Jlen,temp22_2)
!
!      
!      call matmulf(Mh1,twn_1(1,1,jnode),3,3,Jlen,temp32)
!      call matmulf(Mh2,twn_2(1,1,jnode),3,3,Jlen,temp32_2)
!      temp32(1:3,1:Jlen) = temp32(1:3,1:Jlen) + temp32_2(1:3,1:Jlen)
!      call mttmul (twn(1,1,inode),temp32,Ilen,3,Jlen,temp22)
!      Gww(1:Ilen,1:Jlen) 
!     1      = Gww(1:Ilen,1:Jlen)+temp22(1:Ilen,1:Jlen)*shp(3,inode)
!     2                          +temp22_2(1:Ilen,1:Jlen)        
!c     ENDE ALTES VORGEHEN  
      
     
      
!            M(1:3,1:3) = M(1:3,1:3) + M_1(1:3,1:3) + M_2(1:3,1:3)
!      call mttmul (twn(1,1,inode),M,Ilen,3,3,temp23)
!      call matmulf(temp23,twn(1,1,jnode),Ilen,3,Jlen,temp22_3)
!
!      call mttmul (twn_1(1,1,inode),Mh1,Ilen,3,3,temp23)
!      call mttmul (twn_2(1,1,inode),Mh2,Ilen,3,3,temp23_2)
!      temp23(1:Ilen,1:3) = temp23(1:Ilen,1:3) + temp23_2(1:Ilen,1:3)
!      call matmulf(temp23,twn(1,1,jnode),Ilen,3,Jlen,temp22_2)
!!      Gww(1:Ilen,1:Jlen) 
!!     1      = Gww(1:Ilen,1:Jlen)+temp22(1:Ilen,1:Jlen)*shp(3,jnode)
!c
!      
!      call matmulf(Mh1,twn_1(1,1,jnode),3,3,Jlen,temp32)
!      call matmulf(Mh2,twn_2(1,1,jnode),3,3,Jlen,temp32_2)
!      temp32(1:3,1:Jlen) = temp32(1:3,1:Jlen) + temp32_2(1:3,1:Jlen)
!      call mttmul (twn(1,1,inode),temp32,Ilen,3,Jlen,temp22)
!      Gww(1:Ilen,1:Jlen) 
!     1      = Gww(1:Ilen,1:Jlen)+temp22(1:Ilen,1:Jlen)*shp(3,inode)
!     2                          +temp22_2(1:Ilen,1:Jlen)*shp(3,jnode) 
!     2                          +temp22_3(1:Ilen,1:Jlen)   

      return
      end
      
C----------------------------------------------------------------------
      SUBROUTINE MTTAML63 ( A,B,L,M,N, C, alpha )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.....MULTIPLIZIERT DIE MATRIZEN  C(L,N)=C(L,N)+AT(M,L)*B(M,N)
      DIMENSION A(M,L),B(M,N),C(L,N)
      DO 100 J=1,N
      DO 100 K=1,M
      DO 100 I=1,L
  100 C(I,J) = C(I,J) + A(K,I) * B(K,J) * alpha
      END
      
      subroutine computeGww63_Om_Mult(Gww,d,d_1,d_2,h,h_1,
     1                        h_2,h1,h2,twn,twn_1,twn_2,inode,jnode,shp,
     2                        Ilen,Jlen)
c-----------------------------------------------------------------------
c     matrix M(h) and M_,alpha(h) for geometrical matrix 
c     and combination to G^ßß
c     initial version with simplified M matrices for multiplicative update, works for
c     2,3 and 2,2 and is checked, but very slow     
c-----------------------------------------------------------------------
c...  output variables
      double precision Gww(3,3)
c...  input variables
      double precision d(3),d_1(3),d_2(3),h(3)
      double precision h_1(3),h_2(3),h1(3),h2(3)
      double precision twn(3,3,*),twn_1(3,3,*),twn_2(3,3,*),shp(3,*)
      integer inode,jnode,Ilen,Jlen
c...  internal variables
      double precision M(3,3),M_1(3,3),M_2(3,3),Mh1(3,3),Mh2(3,3)
      double precision temp22(Ilen,Jlen),temp23(Ilen,3),temp33(3,3)
      integer i,j
c...  declaration of data type for used defined functions
      double precision dot


c
c...  compute M(h)
      do i = 1,3
        do j = 1,i
          M(i,j) = 0.5d0*(d(i)*h(j) + h(i)*d(j)) 
          M(j,i) = M(i,j)
        end do
        M(i,i) = M(i,i) - dot(d,h,3)
      end do
c
c...  compute M(h1)
      do i = 1,3
        do j = 1,i
          Mh1(i,j) = 0.5d0*(d(i)*h1(j) + h1(i)*d(j)) 
          Mh1(j,i) = Mh1(i,j)
        end do
        Mh1(i,i) = Mh1(i,i) - dot(d,h1,3)
      end do
c
c...  compute M(h2)
      do i = 1,3
        do j = 1,i
          Mh2(i,j) = 0.5d0*(d(i)*h2(j) + h2(i)*d(j)) 
          Mh2(j,i) = Mh2(i,j)
        end do
        Mh2(i,i) = Mh2(i,i) - dot(d,h2,3)
      end do
c
c     compute derived M: M_,1(h_,1)
      do i = 1,3
        do j =1,i
          M_1(i,j) = 0.5d0*(d_1(i)*h_1(j) + h_1(i)*d_1(j))
          M_2(i,j) = 0.5d0*(d_2(i)*h_2(j) + h_2(i)*d_2(j))
          M_1(j,i) = M_1(i,j)
          M_2(j,i) = M_2(i,j)
        end do
        M_1(i,i) = M_1(i,i) - dot(d_1,h_1,3)
        M_2(i,i) = M_2(i,i) - dot(d_2,h_2,3)
      end do
      
c     
      Gww=0.d0   
      call mttmul (twn_1(1,1,inode),Mh1,Ilen,3,3,temp23)
      call matmulf(temp23,twn(1,1,jnode),Ilen,3,Jlen,temp22)
      do i = 1,Ilen
        do j = 1, Jlen
          Gww(i,j) = Gww(i,j) + temp22(i,j)*shp(3,jnode)
        end do
      end do
c
      call mttmul (twn(1,1,inode),Mh1,Ilen,3,3,temp23)
      call matmulf(temp23,twn_1(1,1,jnode),Ilen,3,Jlen,temp22)
      do i = 1,Ilen
        do j = 1, Jlen
          Gww(i,j) = Gww(i,j) + temp22(i,j)*shp(3,inode)
        end do
      end do
c      
      call mttmul (twn_2(1,1,inode),Mh2,Ilen,3,3,temp23)
      call matmulf(temp23,twn(1,1,jnode),Ilen,3,Jlen,temp22)
      do i = 1,Ilen
        do j = 1, Jlen
          Gww(i,j) = Gww(i,j) + temp22(i,j)*shp(3,jnode)
        end do
      end do
c
      call mttmul (twn(1,1,inode),Mh2,Ilen,3,3,temp23)
      call matmulf(temp23,twn_2(1,1,jnode),Ilen,3,Jlen,temp22)
      do i = 1,Ilen
        do j = 1, Jlen
          Gww(i,j) = Gww(i,j) + temp22(i,j)*shp(3,inode)
        end do
      end do
c
      call matadd (M_1,M_2,3,3,temp33)
      call matadd (M,temp33,3,3,M)
      call mttmul (twn(1,1,inode),M,Ilen,3,3,temp23)
      call matmulf(temp23,twn(1,1,jnode),Ilen,3,Jlen,temp22)
      do i = 1,Ilen
        do j = 1, Jlen
          Gww(i,j) = Gww(i,j) + temp22(i,j)
        end do
      end do
      return
      end
      

      subroutine pdirecread63(ix,dir,diri,i,isw,ifdir,iuseTgp,diriP,n)
      USE dirdat
      USE isogeo
      USE cdata
c-----------------------------------------------------------------------
c
c      Purpose: read/save nodal triad for global(!) node i on actual 
c               element for averaging of global nodal values
c
c      Inputs:
c
c         ix(nen1,*)      - Element nodal connections of mesh
c         dir(10,knode,2) - triad field 1-9 = dir, 10 = no.of el.per node
c         dir( , ,ifdir)  - save   
c         dir( , ,ifdir)  - test, update at TIME
c         diri(3,3)       - triad for node i, global values
c         diriP(3,3)      - triad for node i, values for current patch
c         i               - local  node number
c         isw             - 1=write 2=read
c         ifdir           - 1=save 2=test, update at TIME
c
c      Outputs:
c
c         dir(10,knode,2) - triad field 1-9 = dir, 10 = no.of el.per node
c
c.... ww bs uka 12/96+2/04+10/05
c-----------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      include 'dirdat.h'
c      include 'cdata.h'
c      include 'isogeo.h'
csk    common m(1)
      dimension dir(10,knode,2),diri(3,3),ix(*),diriP(3,3)
      if(isw.eq.1) then
c.....  write
        ii = abs(ix(i))
        if(ii.eq.0) return
        do k = 1,3
          dir(k  ,ii,ifdir) = dir(k  ,ii,ifdir) + diri(k,1)
          dir(k+3,ii,ifdir) = dir(k+3,ii,ifdir) + diri(k,2)
          dir(k+6,ii,ifdir) = dir(k+6,ii,ifdir) + diri(k,3)
        end do
        dir(10,ii,ifdir) = dir(10,ii,ifdir) + 1
      else if(isw.eq.2) then
c.....  read
        call pzero(diri,9)
        ii = abs(ix(i))
        if(ii.eq.0) return
        do k = 1,3
          diri(k,1) = dir(k  ,ii,ifdir)
          diri(k,2) = dir(k+3,ii,ifdir)
          diri(k,3) = dir(k+6,ii,ifdir)
        end do
!        diri(1,1)=0.577350269189626
!        diri(2,1)=0.577350269189626
!        diri(3,1)=0.577350269189626
!        diri(1,2)=-0.577350269189626
!        diri(2,2)=0.577350269189626
!        diri(3,2)=0.577350269189626
!        diri(1,3)=0.577350269189626
!        diri(2,3)=-0.577350269189626
!        diri(3,3)=0.577350269189626
        call plamina_base(diri(1,1))
        if (iuseTgp.ge.1) then
csk          call evaluateoldix(m(Inoix),numel,nen,i,n,iiold)
          call evaluateoldix(AInoix,numel,nen,i,n,iiold)
          call pzero(diriP,9)
          if(iiold.eq.0) then
            write(*,*) 'COULD NOT FIND ORIGINAL NODE NUMBER'
            return
          end if 
          if (IldirP.eq.1) then
            call IGAreadpatchdirectors(AImdirP,numnp,diriP,iiold)
          else
            do k = 1,3
              diriP(k,1) = dir(k  ,iiold,ifdir)
              diriP(k,2) = dir(k+3,iiold,ifdir)
              diriP(k,3) = dir(k+6,iiold,ifdir)
            end do
          end if
          if (iuseTgp.eq.1) then
!c           write global nodal basis system for A_\alpha          
!            diriP(1:3,1) = diri(1:3,1)
!            diriP(1:3,2) = diri(1:3,2)
!c           and keep patch nodal basis system for A_3=D    
c           write patch nodal basis system for A_3=D 
            diri(1:3,3) = diriP(1:3,3)
          end if
        end if
      end if
      return
      end
c
      subroutine InitRotationMatrix63(init,n,Rnold,numel,nen)
c-----------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer init(numel),i
      real*8 Rnold(3,3,nen)
      if (init(n).eq.0) then
        do i = 1,nen
          Rnold(1,1,i) = 1.d0
          Rnold(2,1,i) = 0.d0
          Rnold(3,1,i) = 0.d0
          Rnold(1,2,i) = 0.d0
          Rnold(2,2,i) = 1.d0
          Rnold(3,2,i) = 0.d0
          Rnold(1,3,i) = 0.d0
          Rnold(2,3,i) = 0.d0
          Rnold(3,3,i) = 1.d0          
        end do
        init(n) = 1
      end if  
      return
      end
      subroutine InitRotationMatrix63CPGP(init,n,Rnold,Riold,
     +                                    numel,nen,ngb)
c-----------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer init(numel),i
      real*8 Rnold(3,3,nen)
      real*8 Riold(3,3,3,ngb*ngb)
      if (init(n).eq.0) then
        do i = 1,nen
          Rnold(1,1,i) = 1.d0
          Rnold(2,1,i) = 0.d0
          Rnold(3,1,i) = 0.d0
          Rnold(1,2,i) = 0.d0
          Rnold(2,2,i) = 1.d0
          Rnold(3,2,i) = 0.d0
          Rnold(1,3,i) = 0.d0
          Rnold(2,3,i) = 0.d0
          Rnold(3,3,i) = 1.d0          
        end do
        do i = 1,ngb*ngb
          Riold(1,1,3,i) = 1.d0
          Riold(2,1,3,i) = 0.d0
          Riold(3,1,3,i) = 0.d0
          Riold(1,2,3,i) = 0.d0
          Riold(2,2,3,i) = 1.d0
          Riold(3,2,3,i) = 0.d0
          Riold(1,3,3,i) = 0.d0
          Riold(2,3,3,i) = 0.d0
          Riold(3,3,3,i) = 1.d0 
          do j1=1,2
            do j2=1,3
              do j3=1,3
                Riold(j3,j2,j1,i) = 0.0d0
              end do
            end do
          end do 
        end do
        init(n) = 1
      end if  
      return
      end      