      subroutine elmt62(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
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
c.... output of stresses see Subroutine plsn62 
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
     1          st(nst,nst),pt(nst),h1(*),h2(*),h3(*),   
     2          shp(3,nen),shp1(3,4),sg(81),tg(81),wg(81),t0(3,3),
     3          dmat(5,5),sig(5),eps(5),epse(5),b(5,6,nen),btd(5,6),
     4          xs0(2,2),sx0(2,2),sx(2,2),xu(3,3),
     5          rd(3,3),twn(3,3,nen),rwn(3,3,nen),
     6     w(3,3,nen),x0(3,3),r0(3,3),te0(5,5),ts0(5,5),qs(2),siga(9),
     7          gam(4),ran(3,4),xan(3,4),bml(3,2,4),ddn(3,nen),
     8          pgz(5),wgz(5),sigp(6),epstr(5),cmat(5,5),
     9          gxy(5,neas_max),gtd(neas_max,5),eas(5),
     1          dd(3,3),w_1(3,3,nen),w_2(3,3,nen),
     2          twn_1(3,3,nen),twn_2(3,3,nen),iMsize(nen),t0exact(3,3,3)
      real*8 dwnold(3,nen),ddwn(3,nen),Mhq(3,3),Mh1(3,3),Mh2(3,3)

      dimension rinG(3,3,nen),rinP(3,3,nen),rwnPn(3,3,nen),
     1          rwnPnp1(3,3,nen),rwnGn(3,3,nen),rwnGnp1(3,3,nen),
     2          dwn(3,nen),dwg(9)
      save engyt
c
c.... go to correct array processor
      ngeas_max=nen*ndf
      go to(1,2,3,3,2,3,2,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2), isw

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
c.... read material data  
c.... small  strain isotropic      
      if(ior.lt.0) write(*,1008)
1008  format('Input mat. parameter:E,nu >',$)
      call dinput(d(40),2)
c
c.... write properties
      if(ior.lt.0) write(*  ,1006) (d(i),i=40,41)
                   write(iow,1006) (d(i),i=40,41) 
1006  format(/,5x,'Small  strain isotropic',/,
     + 5x,'Youngs modulus E ......................... ',g12.4,/,
     + 5x,'Poissons ratio nu ........................ ',g12.4)
c
c.... initialize h-array    
      nh3   = 0                 ! total number per element in h3
c
      ngb = d(21)
      if (ngb.eq.0) then
c         if number of Gauss Points is not specified in input file  
          write(*,*) 'number of Gauss Points not specified->error in his
     +tory array WILL occur'
      end if

c      
c
2     return
3     continue
      isymm  = d(4)
      ngb = d(21)
c
c.... get current patch and coressponding orders
      NURp = ngetNURnmpq(n,3,AInipa,AInmpq)
      NURq = ngetNURnmpq(n,4,AInipa,AInmpq)
      ngz  = d(12)
c      
c.... names for stresses
      if(isw.eq.8) then
        strsus(1) = 'Resultant M_11 '
        strsus(2) = 'Resultant M_22 '
        strsus(3) = 'Resultant M_12 '
        strsus(4) = 'Resultant Q_1 '
        strsus(5) = 'Resultant Q_2 '
      end if
c      
c     HERE THE REAL ELEMENT STARTS
c.... compute NURBS coordinates of recent element      
      ni = ngetNURni(n,AIninc,AInien,AInipa,1)
      nj = ngetNURni(n,AIninc,AInien,AInipa,2)
c.... check if element has zero measure
      if (rNURknv1(ni).eq.rNURknv1(ni+1).or.
     +    rNURknv2(nj).eq.rNURknv2(nj+1)) return
c
      call pzero (pt,ngeas_max)            !sk needs to be changed
      call pzero (st,ngeas_max*ngeas_max)  !sk needs to be changed
      call IGAgauss (ngb,lint,sg,tg,wg)
      t0=0.d0
      t0(1,1)=1.d0
      t0(2,2)=1.d0
      t0(3,3)=1.d0
      cappa  =  5.0d0/6.d0 
      call dmate62 (d,cappa,dmat)
c
c
c.... loop over gauss points
      do 300 l = 1,lint
         xsi = sg(l)
         eta = tg(l)
         call shap62(xsi,eta,xl,t0,sx,shp,shp1,xsj,ndm,nen,ni,nj,NURp,
     +            NURq,AInkv1,AInkv2)    !CHANGED
         da = xsj*wg(l)
         surface=surface+da    !   for test of surface area
c        compute strains
         call stra62 (ul,shp,ndf,eps,nen)
         
c
         call mvmul(dmat,eps,5,5,sig)
c    
         if (isw.eq.4.or.isw.eq.8) go to 300
c....    loop over rows, residual
         ir0 = 0
         do 310 inode = 1,nen
            jnode=1
            call bmat62 (inode,ir0,b,btd,dmat,d,shp,pt,sig,da,ndf,jnode)
c....       loop over columns (symmetry noted), stiffness matrix
            jc0 = 0       
            do jnode = 1,inode
               call stif62 (inode,jnode,ir0,jc0,b,btd,st,ndf,nst,isw)
               jc0 = jc0 + ndf
            end do
            ir0 = ir0 + ndf
310      continue
c         
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
c.... plot stresses   =>  isw 3!
8     return
c     error analysis  => isw 3
9     return  
      end
c-------------------------------------------------------------------------------------

      subroutine shap62(xsi,eta,xl,t0,sx,shp,shp1,xsj,ndm,nel,ni,nj,NURp
     +                 ,NURq,NURknv1,NURknv2)
      USE isogeo
      USE iofile
c----------------------------------------------------------------------
c
c      Purpose: Computes shape function and derivatives for
c               isogeometric surface elements
c
c      Inputs:
c         xsi       - Natural coordinates for point (unit element coordinates)
c         eta       - Natural coordinates for point (unit element coordinates)
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh
c         nel       - Number of nodes on element
c         ix(*)     - Nodes attached to element
c         flg       - Flag, compute global x/y derivatives if false,
c                           else derivatives are w/r natural coords.
c         ni        - number of highest non-zero basis function in 1-direction
c         nj        - number of highest non-zero basis function in 2-direction
c
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
c      uset0local = .true.    !absolutely needed for perpendicular element loads
c      if (uset0local) then
cc       compute lamina basis system (local cartesian), see Hughes, FEM, 2000, p. 385
c        call norm (t0(1,3),dx_dXi(1,3),3)
c        call norm (gs(1,1),dx_dXi(1,1),3)
c        call norm (gs(1,2),dx_dXi(1,2),3)
c        do i = 1,3
c            gs(i,1) = gs(i,1) + gs(i,2)
c        end do
c        call vecp (t0(1,3),gs(1,1),gs(1,2))
c        do i = 1,3
c            t0(i,1) = gs(i,1) - gs(i,2)
c            t0(i,2) = gs(i,1) + gs(i,2)
c        end do
c        call norm (t0(1,1),t0(1,1),3)
c        call norm (t0(1,2),t0(1,2),3)
c        !t0(1:3,1:3) = dx_dXi(1:3,1:3)
c        !call plamina_base(t0)
c      end if
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
c
      subroutine bmat62 (inode,ir0,b,btd,dmat,d,shp,pt,sig,da,ndf,jnode)
      USE prlod
      USE eldata
      USE isogeo
c-----------------------------------------------------------------------
c.... b-matrix, residual vector, external load vector  
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),shp(3,*),b(5,3,*),btd(5,*),dmat(5,5),sig(*),
     1          pt(*),ql(3)
c
c.... b-matrix: membrane (1,2,3), bending (4,5,6), shear (7,8) 
      b(1,1,inode) = 0.d0
      b(1,2,inode) = shp(1,inode)
      b(1,3,inode) = 0.d0
      b(2,1,inode) = 0.d0
      b(2,2,inode) = 0.d0
      b(2,3,inode) = shp(2,inode)
      b(3,1,inode) = 0.d0
      b(3,2,inode) = shp(2,inode)
      b(3,3,inode) = shp(1,inode)
      b(4,1,inode) = shp(1,inode)
      b(4,2,inode) = shp(3,inode)
      b(4,3,inode) = 0.d0
      b(5,1,inode) = shp(1,inode)
      b(5,2,inode) = 0.d0
      b(5,3,inode) = shp(3,inode)
c      
c.... residual vector and  B_T * D
      do i = 1,ndf
        ir = ir0 + i
        if (jnode.eq.1) then
        pt(ir) = pt(ir) - dot(sig,b(1,i,inode),5)*da
        end if
        do j = 1,5
           btd(j,i) = dot(b(1,i,inode),dmat(1,j),5)*da
       enddo
      enddo 
c
c.... compute element load vector
c.... external constant loads qx, qy, qz in x,y,z-direction
      ql(1) = d(26)
      ql(2) = d(27)
      ql(3) = d(28)
c.... add loads to load vector
      if (jnode.eq.1) then
         do i = 1,3
            pt(ir0+i) = pt(ir0+i) + prop*ql(i) * shp(3,inode)*da
         end do
      end if
c
c
      return
      end

c
      subroutine boun62 (inode,ix,id,iflag,ndf)
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

      subroutine dmate62 (d,cappa,dmat)
c-----------------------------------------------------------------------
c     d(13) = shell thickness h_s
c     d(14) = distance bottom surface to reference surface 
c     d(40) = elasticity modulus
c     d(41) = poisson's ratio
c     cmat  = elasticity matrix
c     dmat  = elasticity matrix for stress resultants
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension d(*),cmat(3,3),dmat(5,*)
c
      call pzero(cmat,3*3)
      call pzero(dmat,5*5)
c
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
c          dmat(i  ,j  ) = cmat(i,j)*cm
c          dmat(i+3,j  ) = cmat(i,j)*cmb
c          dmat(i  ,j+3) = cmat(i,j)*cmb
          dmat(i,j) = cmat(i,j)*cb
        end do
      end do
c
c.... shear
      dmat(4,4) = 0.5d0*d(40)/(1.0d0+d(41))*cs  
      dmat(5,5) = dmat(5,5)
c
      return
      end
c
      subroutine kine62(eps,x0,r0,zs,wcappa,epstr,dets)
      implicit double precision (a-h,o-z)
      dimension  glu(3,3),x0(3,3),r0(3,3),eps(*),epstr(*)
c
c....  layer strains  E_n+1 = A * E  and  det shifter tensor
c
cjk     epstr(1) = eps(1) + zs*eps(4) 
cjk       epstr(2) = eps(2) + zs*eps(5) 
cjk       epstr(3) = eps(3) + zs*eps(6) 
cjk       epstr(4) = eps(7)*wcappa
cjk       epstr(5) = eps(8)*wcappa
 
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
cjk      subroutine dmat62 (dmat,cmat,sig,sigp,wz,zs,wcappa,cappa)
c-----------------------------------------------------------------------
c       sigp  stresses of layer point
c       sig   stress resultants   
c       cmat  consistent tangent matrix 
c       dmat  linearized stress resultants
c-----------------------------------------------------------------------
cjk      implicit double precision(a-h,o-z)
cjk      dimension dmat(5,*),cmat(5,*),sig(*),sigp(*)

cjk      zs2 = zs*zs
      
cjk       do 10 i=1,2
cjk       do 10 j=1,i
cjk        dmat(i  ,j  ) = dmat(i  ,j  ) + cmat(i,j)*wz
cjk        dmat(i+3,j+3) = dmat(i+3,j+3) + cmat(i,j)*wz*zs2
cjk 10    continue

cjk      do 20 i=1,3
cjk        dmat(6,i) = dmat(6,i) + cmat(3,i)*wz*zs
cjk      do 20 j=1,2
cjk        dmat(j+3 ,i)  = dmat(j+3 ,i) + cmat(j,i)*wz*zs
cjk        dmat(j+6 ,i)  = dmat(j+6 ,i) + cmat(j+3,i)*wz       *wcappa
cjk        dmat(j+6,i+3) = dmat(j+6,i+3) + cmat(j+3,i)*wz*zs   *wcappa
cjk 20    continue

cjk      dmat(7,7) = dmat(7,7) + cmat(4,4)*wz*cappa
cjk      dmat(8,7) = dmat(8,7) + cmat(5,4)*wz*cappa
cjk      dmat(8,8) = dmat(8,8) + cmat(5,5)*wz*cappa

c....  upper part of dmat 
cjk      do 30 i = 1,5
cjk      do 30 j = i+1,5
cjk 30    dmat(i,j) = dmat(j,i)

c...  compute stress resultants
cjk      do 40 i=1,3
cjk        sig(i  ) = sig(i  ) + sigp(i)*wz
cjk        sig(i+3) = sig(i+3) + sigp(i)*wz*zs
cjk 40    continue
cjk        sig(7)  = sig(7) + sigp(4)*wz*wcappa
cjk        sig(8)  = sig(8) + sigp(5)*wz*wcappa
cjk      return
cjk      end
c
 
c
      subroutine transsig62(t0,xu,rd,rin,shp,sig,siga,nen)
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
      subroutine plsn62 (i)
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
      subroutine zeroElementStressResults62 (NURstress,numel,n) 
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
      subroutine plot62 (d,ix,dt,st,shp,sig,sigp,h2,eps,
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

c
      subroutine stif62 (inode,jnode,ir0,jc0,b,btd,st,ndf,nst,isw)
c-----------------------------------------------------------------------
c.... stiffness matrix
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  st(nst,*),b(5,3,*),btd(5,*)
c
c.... material part     
      if(isw.eq.3)then
        do i  = 1,ndf
          ir   = ir0+i
          do j = 1,ndf
            jc  = jc0+j
            st(ir,jc) = st(ir,jc) + dot(btd(1,i),b(1,j,jnode),5)
          end do
        end do
      endif
c
      return
      end
c

      subroutine stra62 (ul,shp,ndf,eps,nen)
c-----------------------------------------------------------------------
c     strains of the reference surface
c     curvatures(1,2,3),    shear(4,5) 
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension ul(ndf,*),shp(3,*),eps(5)
c
      do k = 1,nen
          eps(1) = shp(1,k)*ul(2,k)
          eps(2) = shp(2,k)*ul(3,k)
          eps(3) = shp(2,k)*ul(2,k)+shp(1,k)*ul(3,k)
          eps(4) = shp(1,k)*ul(1,k)+shp(3,k)*ul(2,k)
          eps(5) = shp(2,k)*ul(1,k)+shp(3,k)*ul(3,k)
      end do
c   
      return
      end
c
      subroutine trasig62 (sigp,t0,shp,rin,nen)               
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