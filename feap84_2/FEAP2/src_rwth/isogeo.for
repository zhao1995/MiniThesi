c.... isogeo.for aus karlsruhe + fehlende funktionen aus feap_iga
c     rwth - 30.06.2015 sd
csk   common m auskommentiert
csk   pseta durch iallocate bzw. rallocate 
csk   Variablenfelder (siehe module isogeo)
csk   Integer Inmpq,Ininc, Inien, Inipa, Inoix, InRupdinit
csk   Real*8  Inpro, Inkv1, Inkv2, ImdirP, InRupdinit, Inbasissys, Instre
csk           durch AInmpq, .... ersetzen  
      
      subroutine initIGA(ix,nen1,numnp,numel,nen,ipr)
c--------------------------------------------------------------------
c.... initializes all important arrays and variables for IGA
c--------------------------------------------------------------------
c     Input:  from global variables NUR
c     Output: to global variables IGA
      USE isogeo
      USE doalloc
      implicit double precision (a-h,o-z)
      integer ix(nen1,*)
      logical flinc,flien,flipa,floix,flpro,flRupdinit
      call ialloc(AIninc,2*numnp,'INC',flinc)
c      allocate(IGA_INC(1:numnp,1:2))
c      allocate(IGA_IEN(1:nen,1:numel))
c      allocate(IGA_IPA(1:numel,1:3))
c      call pseta(Inoix,nen*numel,1,floix,'old ix')
c      call pseta(Inpro,3*numnp,ipr,flpro,'CP projection')
      call ialloc(AInien,nen*numel,'IEN',flien)
      call ialloc(AInipa,3*numel,'IPA',flipa)
      call ralloc(AInstre,9*numel,'NURStress',flinc)
      call ialloc(AInoix,nen*numel,'old ix',floix)
      call ralloc(AInpro,3*numel,'CP projection',flpro)
c      allocate(NURstress(numel*9,11))
      call buildINC_IEN(numnp,numel,nen,NURnpatch,AIninc,AInien,
     +                  AInipa,AInmpq,ix,nen1)
      call saveoldix(ix,nen1,nen,numel)
      call ialloc(AInRupdinit,numel,'Init R Multipl',flRupdinit)
      IldirP = 0
      nurbs = .true.
      surface=0.0d0
      Iexactnormal = 0
      end
c
      subroutine initIGA2(td,NURnmpq,NURnpatch,i)
c--------------------------------------------------------------------
c.... initializes all array for storage of n,m,p,q and material
c--------------------------------------------------------------------
c     Input:  from global variables NUR
c     Output: to global variables IGA
      implicit double precision (a-h,o-z)
      dimension td(6),NURnmpq(NURnpatch,6)
      NURnmpq(i,1) = td(1)  ! n of patch i
      NURnmpq(i,2) = td(2)  ! m of patch i
      NURnmpq(i,3) = td(3)  ! p of patch i
      NURnmpq(i,4) = td(4)  ! q of patch i
      NURnmpq(i,5) = td(5)  ! matnum of patch i
      NURnmpq(i,6) = td(6)  ! i_patch of patch i
      end
c
c--------------------------------------------------------------------
      integer function ngetnen()
c--------------------------------------------------------------------
      USE cdata
      implicit double precision (a-h,o-z)
      ngetnen = nen
      end
c--------------------------------------------------------------------
c
      integer function ngetNURnmpq(ielem,col,IGA_IPA,NURnmpq)
c--------------------------------------------------------------------
      USE cdata
      USE isogeo
      implicit double precision (a-h,o-z)
      integer ielem,col
      integer IGA_IPA(numel,3),NURnmpq(NURnpatch,6)
      ngetNURnmpq = NURnmpq(IGA_IPA(ielem,1),col)
      end
c
      integer function nNURnmpq(patch,col)
c--------------------------------------------------------------------
      USE isogeo
      USE psize
      implicit double precision (a-h,o-z)
      integer patch,col
      nNURnmpq =nNURnmpq1(patch,col,AInmpq)
      end
c
      integer function nNURnmpq1(patch,col,NURnmpq)
c--------------------------------------------------------------------
      USE isogeo
      implicit double precision (a-h,o-z)
      integer patch,col,NURnmpq(NURnpatch,6)
      nNURnmpq1 = NURnmpq(patch,col)
      end

c
      integer function ngetNURni(ielem,INC,IEN,IPA,dir)
c--------------------------------------------------------------------
      USE cdata
      USE isogeo
      implicit double precision (a-h,o-z)
      integer ielem,dir
      integer IPA(numel,3),INC(numnp,2),IEN(nen,numel)
      if (dir.eq.1) then
        ngetNURni = INC(IEN(1,ielem),1)+IPA(ielem,2)
      elseif (dir.eq.2) then
        ngetNURni = INC(IEN(1,ielem),2)+IPA(ielem,3)
      endif
      end
c
      double precision function rNURstress(node,val)
c--------------------------------------------------------------------
      USE isogeo
      USE psize
      implicit double precision (a-h,o-z)
      integer node,val
      rNURstress = rNURstress1(node,val,AInstre)
      end
c
      double precision function rNURstress1(node,val,NURstress)
c--------------------------------------------------------------------
      USE cdata
      implicit double precision (a-h,o-z)
      integer node,val
      real*8 NURstress(numel*9,11)
      rNURstress1 = NURstress(node,val)
      end
c
c--------------------------------------------------------------------
      double precision function rNURknv1(ni)
c--------------------------------------------------------------------
      USE isogeo
      USE psize
      implicit double precision (a-h,o-z)
      integer ni
      rNURknv1 = rNURknv11(ni,AInkv1)
      end
c
c--------------------------------------------------------------------
      double precision function rNURknv11(ni,knv)
c--------------------------------------------------------------------
      USE isogeo
      implicit double precision (a-h,o-z)
      integer ni
      real*8 knv(NURlenkv(1))
      rNURknv11 = knv(ni)
      end
c--------------------------------------------------------------------
      double precision function rNURknv2(ni)
c--------------------------------------------------------------------
      USE isogeo
      USE psize
      implicit double precision (a-h,o-z)
      integer ni
      rNURknv2 = rNURknv21(ni,AInkv2)
      end
c
c--------------------------------------------------------------------
      double precision function rNURknv21(ni,knv)
c--------------------------------------------------------------------
      USE isogeo
      implicit double precision (a-h,o-z)
      integer ni
      real*8 knv(NURlenkv(2))
      rNURknv21 = knv(ni)
      end
c
c--------------------------------------------------------------------
      logical function isNurbs()
c--------------------------------------------------------------------
      USE isogeo
      isNurbs = nurbs
      end function
c--------------------------------------------------------------------
      subroutine resetNurbs()
c--------------------------------------------------------------------
      USE isogeo
      nurbs = .false.
      return
      end
c
c--------------------------------------------------------------------
      subroutine setKv(pos,val,kv,lenkv)
c--------------------------------------------------------------------
      integer pos,lenkv
      real*8 val,kv(lenkv)
      kv(pos)=val
      return
      end
c
c--------------------------------------------------------------------
      subroutine setKnotVectLength(val,dir)
c--------------------------------------------------------------------
      USE isogeo
      integer val,dir
      NURlenkv(dir) = val
      return
      end
c--------------------------------------------------------------------
      subroutine buildINC_IEN(n_np,n_el,n_en,n_patch,INC,IEN,IPA,
     +                        nmpq,ix,nen1)
c--------------------------------------------------------------------
c.... initializes INC and IEN array
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer p,q, Ee,A,B,Bb
      integer ix(nen1,*),nmpq(n_patch,6)
      dimension INC(n_np,2),IEN(n_en,n_el),IPA(n_el,3)
      A=0;
      B=0;
      Bb=0;
      Ee=0;
      lknv1 = 0
      lknv2 = 0
      do 108 i_patch=1,n_patch
        n = nmpq(i_patch,1)
        m = nmpq(i_patch,2)
        p = nmpq(i_patch,3)
        q = nmpq(i_patch,4)
        do 107 j=1,m
          do 105 i=1,n
            A=A+1;
            INC(A,1)=i;
            INC(A,2)=j;
            if (i .ge. p+1 .and. j .ge. q+1) then
                Ee = Ee+1;
                IPA(Ee,1)=i_patch
                IPA(Ee,2)=lknv1
                IPA(Ee,3)=lknv2
                ix(nen1,Ee)=nmpq(i_patch,5)
                do 103 jloc=0,q
                    do 101 iloc=0,p
                        B = A-jloc*n-iloc;      !global function number
                        Bb = jloc*(p+1)+iloc+1;  !local function number
                        IEN(Bb,Ee)=B;
                        ix(Bb,Ee) = B
101                 end do
103             end do
            endif
105       end do
107     end do
        lknv1 = lknv1 + n+p+1
        lknv2 = lknv2 + m+q+1
108   end do
      end
c
      subroutine IGAshapeDir(Xi1,Xi2,shp,ni,nj,NURCPw,p,q,NURknv1,
     +                       NURknv2,nen,n,m)
c----------------------------------------------------------------------
c      Purpose: Computes shape function for DIRECTOR-CORRECTOR-ALGORITHM
c
c      Inputs:
c         Xi1       - NURBS coordinates for point
c         Xi2       - NURBS coordinates for point
c         NURCPw    - Nodal coordinates
c         nen       - Number of nodes on element
c         ni        - number of highest non-zero basis function in 1-direction
c         nj        - number of highest non-zero basis function in 2-direction
c
c      Outputs:
c         shp(*)  - Shape functions and derivatives at point
c                     shp(i) = N_i
c
c----------------------------------------------------------------------       
      implicit double precision (a-h,o-z)
      integer ni,nj,p,q,nen,n,m
      real*8 Xi1,Xi2,shp(nen),NURCPw(4,n,m),NURknv1(n+p+1)
      real*8 NURknv2(m+q+1),N_Xi1(2,p+1),N_Xi2(2,q+1)
c      
c     compute univariate basis functions and their derivatives
      call DersBasisFuns(ni,Xi1,p,NURknv1,N_Xi1)
      call DersBasisFuns(nj,Xi2,q,NURknv2,N_Xi2)
c
      W=0;
      do j=0,q
        do i=0,p
            W=W+N_Xi1(1,i+1)*N_Xi2(1,j+1)*NURCPw(4,ni-p+i,nj-q+j)
        end do
      end do
      int=0;
      do j=0,q
        do i=0,p
            int=int+1
            shp(int)=N_Xi1(1,p+1-i)*N_Xi2(1,q+1-j)
     +                 *NURCPw(4,ni-i,nj-j)/W
        end do
      end do
      return
      end
c
!      subroutine IGAshape2D(ss,tt,shp,xsj,ni,nj)
!c      subroutine IGAshape2D(ss,tt,x,shp,xsj,ndm,nel,ix,flg,ni,nj)
!c----------------------------------------------------------------------
!c     NOT IN USE: RELYING ON WEIGHTED COORDINATES IN NURCPw
!c      Purpose: Computes shape function and derivatives for
!c               isogeometric surface elements
!c
!c      Inputs:
!c         ss        - Natural coordinates for point (unit element coordinates)
!c         tt        - Natural coordinates for point (unit element coordinates)
!c         x(ndm,*)  - Nodal coordinates for element
!c         ndm       - Spatial dimension of mesh
!c         nel       - Number of nodes on element
!c         ix(*)     - Nodes attached to element
!c         flg       - Flag, compute global x/y derivatives if false,
!c                           else derivatives are w/r natural coords.
!c         ni        - number of highest non-zero basis function in 1-direction
!c         nj        - number of highest non-zero basis function in 2-direction
!c
!c      Outputs:
!c         shp(3,*)  - Shape functions and derivatives at point
!c                     shp(1,i) = dN_i/dx or dN_i/dxi_1
!c                     shp(2,i) = dN_i/dy or dN_i/dxi_2
!c                     shp(3,i) = N_i
!c         xsj       - Jacobian determinant at point
!c
!c----------------------------------------------------------------------
!      USE isogeo
!      implicit double precision (a-h,o-z)
!c      logical flg
!      dimension shp(3,IGAn_en),dR_dXi(IGAn_en,2)
!c      dimension xs(3,2),sx(2,2),ix(*),s(4),t(4),x(ndm,nel)
!      dimension xsjMat(2,2), dx_dXi(2,2)
!      real*8 N_Xi1(2,NURp+1),N_Xi2(2,NURq+1)
!c
!c     calculate NURBS parametric coordinates from unit element coordinates
!      Xi1=0.5d0*((NURknv1(ni+1)-NURknv1(ni))*ss+(NURknv1(ni+1)+
!     + NURknv1(ni)))
!      Xi2=0.5d0*((NURknv2(nj+1)-NURknv2(nj))*tt+(NURknv2(nj+1)+
!     + NURknv2(nj)))
!c
!c     compute univariate basis functions and their derivatives
!      call DersBasisFuns(ni,Xi1,NURp,NURknv1,N_Xi1)
!      call DersBasisFuns(nj,Xi2,NURq,NURknv2,N_Xi2)
!c
!      W=0;
!      W_Xi1=0;
!      W_Xi2=0;
!      do j=0,NURq
!        do i=0,NURp
!            W=W+N_Xi1(1,i+1)*N_Xi2(1,j+1)*NURCPw(4,ni-NURp+i,nj-NURq+j)
!            W_Xi1=W_Xi1+N_Xi1(2,i+1)*N_Xi2(1,j+1)
!     +            *NURCPw(4,ni-NURp+i,nj-NURq+j)
!            W_Xi2=W_Xi2+N_Xi1(1,i+1)*N_Xi2(2,j+1)
!     +            *NURCPw(4,ni-NURp+i,nj-NURq+j)
!        end do
!      end do
!      int=0;
!      do j=0,NURq
!        do i=0,NURp
!            int=int+1
!            shp(3,int)=N_Xi1(1,NURp+1-i)*N_Xi2(1,NURq+1-j)
!     +                 *NURCPw(4,ni-i,nj-j)/W
!            dR_dXi(int,1)=NURCPw(4,ni-i,nj-j)*(N_Xi1(2,NURp+1-i)
!     +                 *N_Xi2(1,NURq+1-j)*W-W_Xi1*N_Xi1(1,NURp+1-i)
!     +                 *N_Xi2(1,NURq+1-j))/(W**2)
!            dR_dXi(int,2)=NURCPw(4,ni-i,nj-j)*(N_Xi1(1,NURp+1-i)
!     +                 *N_Xi2(2,NURq+1-j)*W-W_Xi2*N_Xi1(1,NURp+1-i)
!     +                 *N_Xi2(1,NURq+1-j))/(W**2)
!        end do
!      end do
!c
!c     gradient of mapping from parameter space dXi to physical space dx
!      dx_dXi=0
!      loc_num=0
!      do j=0,NURq
!        do i=0,NURp
!            loc_num=loc_num+1
!            do ia=1,2
!                do ib=1,2
!                    dx_dXi(ia,ib)=dx_dXi(ia,ib)+NURCPw(ia,ni-i,nj-j)
!     +               *dR_dXi(loc_num,ib)/NURCPw(4,ni-i,nj-j)
!                end do
!            end do
!        end do
!      end do
!c
!c     determination of Jacobian matrix
!      xsjMat(:,1)=dx_dXi(:,1)*(NURknv1(ni+1)-NURknv1(ni))/2
!      xsjMat(:,2)=dx_dXi(:,2)*(NURknv2(nj+1)-NURknv2(nj))/2
!      xsj=xsjMat(1,1)*xsjMat(2,2)-xsjMat(1,2)*xsjMat(2,1)
!      if(xsj.le.0.0d0) then
!        call drawmess('negative Jacobian in shapef',1,-2)
!        return
!      end if
!c     compute inverse of dx_dXi
!      call invert(dx_dXi,2,2)
!      shp(1:2,:)=transpose(matmul(dR_dXi,dx_dXi))
!c      call Matmulf(dR_dXi,dXi_dx,IGAn_en,2,2,dR_dx)
!c      shp(1:2,:)=transpose(dR_dx)
!      return
!      end
c
      subroutine IGAshape2Dx(ss,tt,x,shp,xsj,ndm,nel,ni,nj,NURp,
     +                       NURq,NURknv1,NURknv2)
c----------------------------------------------------------------------
c
c      Purpose: Computes shape function and derivatives for
c               isogeometric surface elements
c
c      Inputs:
c         ss        - Natural coordinates for point (unit element coordinates)
c         tt        - Natural coordinates for point (unit element coordinates)
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
      USE isogeo
      implicit double precision (a-h,o-z)
c      logical flg
      dimension shp(3,(NURp+1)*(NURq+1)),dR_dXi((NURp+1)*(NURq+1),2)
      dimension x(ndm,nel)
      dimension xsjMat(2,2), dx_dXi(2,2)
      real*8 N_Xi1(2,NURp+1),N_Xi2(2,NURq+1)
      real*8 NURknv1(NURlenkv(1)),NURknv2(NURlenkv(2))
c
c     calculate NURBS parametric coordinates from unit element coordinates
      Xi1=0.5d0*((NURknv1(ni+1)-NURknv1(ni))*ss+(NURknv1(ni+1)+
     + NURknv1(ni)))
      Xi2=0.5d0*((NURknv2(nj+1)-NURknv2(nj))*tt+(NURknv2(nj+1)+
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
            W=W+N_Xi1(1,NURp+1-i)*N_Xi2(1,NURq+1-j)*x(4,int)
            W_Xi1=W_Xi1+N_Xi1(2,NURp+1-i)*N_Xi2(1,NURq+1-j)*x(4,int)
            W_Xi2=W_Xi2+N_Xi1(1,NURp+1-i)*N_Xi2(2,NURq+1-j)*x(4,int)
        end do
      end do
      int=0;
      do j=0,NURq
        do i=0,NURp
            int=int+1
            shp(3,int)=N_Xi1(1,NURp+1-i)*N_Xi2(1,NURq+1-j)
     +                 *x(4,int)/W
            dR_dXi(int,1)=x(4,int)*(N_Xi1(2,NURp+1-i)
     +                 *N_Xi2(1,NURq+1-j)*W-W_Xi1*N_Xi1(1,NURp+1-i)
     +                 *N_Xi2(1,NURq+1-j))/(W**2)
            dR_dXi(int,2)=x(4,int)*(N_Xi1(1,NURp+1-i)
     +                 *N_Xi2(2,NURq+1-j)*W-W_Xi2*N_Xi1(1,NURp+1-i)
     +                 *N_Xi2(1,NURq+1-j))/(W**2)
        end do
      end do
c
c     gradient of mapping from parameter space dXi to physical space dx
      dx_dXi=0
      loc_num=0
      do j=0,NURq
        do i=0,NURp
            loc_num=loc_num+1
            do ia=1,2
                do ib=1,2
                    dx_dXi(ia,ib)=dx_dXi(ia,ib)+x(ia,loc_num)
     +               *dR_dXi(loc_num,ib)
                end do
            end do
        end do
      end do
c
c     determination of Jacobian matrix
      xsjMat(:,1)=dx_dXi(:,1)*(NURknv1(ni+1)-NURknv1(ni))/2
      xsjMat(:,2)=dx_dXi(:,2)*(NURknv2(nj+1)-NURknv2(nj))/2
      xsj=xsjMat(1,1)*xsjMat(2,2)-xsjMat(1,2)*xsjMat(2,1)
      if(xsj.le.0.0d0) then
        call drawmess('negative Jacobian in shapef',1,-2)
        return
      end if
c     compute inverse of dx_dXi
      a=dx_dXi(1,1);b=dx_dXi(1,2);c=dx_dXi(2,1);d=dx_dXi(2,2)
      dx_dXi(1,1)=d;
      dx_dXi(1,2)=-b;
      dx_dXi(2,1)=-c;
      dx_dXi(2,2)=a;
      dx_dXi=dx_dXi/(a*d-b*c)
      shp(1:2,:)=transpose(matmul(dR_dXi,dx_dXi))
c      call Matmulf(dR_dXi,dXi_dx,IGAn_en,2,2,dR_dx)
c      shp(1:2,:)=transpose(dR_dx)
      return
      end
c
      subroutine DersBasisFuns(ni,Xi,p,KV,ders)
c----------------------------------------------------------------------
c
c      Purpose: Computes Basis Functions and their first derivative
c               for B-Splines, Algorithm A2.3 from Piegl/Tiller
c
c      Inputs:
c         ni        - number of highest non-zero basis function
c         Xi        - NURBS parametric coordinate
c         p         - order of NURBS
c         KV        - NURBS knot vector
c
c      Outputs:
c         ders      - Basis functions and their first derivative
c
c----------------------------------------------------------------------
c     input/output variable declaration
      integer p,ni
      real*8 ders(2,p+1),KV(*), Xi
c     functions needed within subroutine
      real*8 left(p+1),right(p+1),saved,ndu(p+1,p+1),a(2,p+1),temp,d
      integer j,r
c
      left=0d0; right=0d0; ndu=0d0; ndu(1,1)=1d0; ders=0d0; a=0d0
c
c     store basis function and knot differences to ndu
      do j=1,p
        left(j) = Xi-KV(ni+1-j);
        right(j) = KV(ni+j)-Xi;
        saved = 0d0;
        do r=0,j-1
            ndu(j+1,r+1)=right(r+1)+left(j-r);
            temp = ndu(r+1,j)/ndu(j+1,r+1);
            ndu(r+1,j+1)=saved+right(r+1)*temp;
            saved = left(j-r)*temp;
        end do
        ndu(j+1,j+1)=saved;
      end do
c
c     Load basis funtions
      do j=0,p
        ders(1,j+1)=ndu(j+1,p+1);
      end do
c     compute the first derivative
      do r=0,p
        a(1,1)=1d0;
c       here could be a loop for k to determine higher derivatives
        d=0d0;
        if (r>=1) then
            a(2,1)=a(1,1)/ndu(p+1,r)
            d=a(2,1)*ndu(r,p)
        end if
        if (r<=p-1) then
            a(2,2)=-a(1,1)/ndu(p+1,r+1)
            d=d+a(2,2)*ndu(r+1,p)
        end if
        ders(2,r+1)=d;
c       end of possible loop
      end do
c     Multiply with correct factors
c     for j=0:p
c       ders(2,j+1)=ders(2,j+1)*p;
c     end
      ders(2,:)=ders(2,:)*p;
c     Transpose results
c     ders=ders';
c     N_Xi=ders(1,:);
c     N_Xi=ders(2,:);
      return
      end
c
c
      subroutine IGAgaussResults(l,lint,sgN,tgN,w)
c----------------------------------------------------------------------
c     Gives back Gauss coordinates for edge and midside points for
c     calculating stress results for these points in IGA, Gauss weights are not needed
c
      implicit double precision (a-h,o-z)
      dimension sgN(*),tgN(*),w(*)
      lint=l*l
      sgN(1:3)=-1.0d0;sgN(4:6)=0.0d0;sgN(7:9)=1.0d0
      tgN(1)=-1.0d0;tgN(2)=0.0d0;tgN(3)=1.0d0
      tgN(4)=-1.0d0;tgN(5)=0.0d0;tgN(6)=1.0d0
      tgN(7)=-1.0d0;tgN(8)=0.0d0;tgN(9)=1.0d0
      w(1:9)=0.0d0
      return
      end
c
c
      subroutine IGAgauss(l,lint,r,z,w)
c----------------------------------------------------------------------
c
c      Purpose: Computes gauss integration values in right order for IGA
c               Values taken from Bathe, Tabelle 5.6 (p. 542) and from
c  http://projecteuclid.org/DPubS/Repository/1.0/Disseminate?view=body&id=pdf_1&handle=euclid.bams/1183504772
c      Inputs:
c         l       - Number of points/direction
c
c      Outputs:
c         lint    - Total number of points
c         r(*)    - 1-direction Gauss point
c         z(*)    - 2-direction Gauss point
c         w(*)    - Gauss weight
c
c----------------------------------------------------------------------
c      USE eldata
      implicit double precision (a-h,o-z)
      dimension r(*),z(*),w(*),gp(l),gw(l)
      lint = l*l
c
      select case (l)
      case (1)
        gp(1)=0d0;
        gw(1)=2d0;
      case (2)
        gp(1)=-0.577350269189626;
        gp(2)=-gp(1);
        gw(1)=1d0;
        gw(2)=1d0;
      case (3)
        gp(1)=-0.774596669241483;
        gp(2)=0d0;
        gp(3)=-gp(1);
        gw(1)=0.555555555555556;
        gw(2)=0.888888888888889;
        gw(3)=gw(1);
      case (4)
        gp(1)=-0.861136311594053;
        gp(2)=-0.339981043584856;
        gp(3)=-gp(2);
        gp(4)=-gp(1);
        gw(1)=0.347854845137454;
        gw(2)=0.652145154862546;
        gw(3)=gw(2);
        gw(4)=gw(1);
      case (5)
        gp(1)=-0.906179845938664;
        gp(2)=-0.538469310105683;
        gp(3)=0d0;
        gp(4)=-gp(2);
        gp(5)=-gp(1);
        gw(1)=0.236926885056189;
        gw(2)=0.478628670499366;
        gw(3)=0.568888888888889;
        gw(4)=gw(2);
        gw(5)=gw(1);
      case (6)
        gp(1)=-0.932469514203152
        gp(2)=-0.661209386466265
        gp(3)=-0.238619186083197
        gp(4)=0.238619186083197
        gp(5)=0.661209386466265
        gp(6)=0.932469514203152
        gw(1)=0.171324492379170
        gw(2)=0.360761573048139
        gw(3)=0.467913934572691
        gw(4)=0.467913934572691
        gw(5)=0.360761573048139
        gw(6)=0.171324492379170
      case (7)
        gp(1)=-0.949107912342759
        gp(2)=-0.741531185599394
        gp(3)=-0.405845151377397
        gp(4)=0.0d0
        gp(5)=-gp(3)
        gp(6)=-gp(2)
        gp(7)=-gp(1)
        gw(1)=0.129484966168870
        gw(2)=0.279705391489277
        gw(3)=0.381830050505119
        gw(4)=0.417959183673469
        gw(5)=0.381830050505119
        gw(6)=0.279705391489277
        gw(7)=0.129484966168870
      case (8)
        gp(1)=-0.960289856497536
        gp(2)=-0.796666477413627
        gp(3)=-0.525532409916329
        gp(4)=-0.183434642495650
        gp(5)=-gp(4)
        gp(6)=-gp(3)
        gp(7)=-gp(2)
        gp(8)=-gp(1)
        gw(1)=0.101228536290376
        gw(2)=0.222381034453374
        gw(3)=0.313706645877887
        gw(4)=0.362683783378362
        gw(5)=0.362683783378362
        gw(6)=0.313706645877887
        gw(7)=0.222381034453374
        gw(8)=0.101228536290376
      case (9)
        gp(1)=-0.968160239507626
        gp(2)=-0.836031107326636
        gp(3)=-0.613371432700590
        gp(4)=-0.324253423403809
        gp(5)=0.0d0
        gp(6)=-gp(4)
        gp(7)=-gp(3)
        gp(8)=-gp(2)
        gp(9)=-gp(1)
        gw(1)=0.081274388361574
        gw(2)=0.180648160694857
        gw(3)=0.260610696402935
        gw(4)=0.312347077040003
        gw(5)=0.330239355001260
        gw(6)=0.312347077040003
        gw(7)=0.260610696402935
        gw(8)=0.180648160694857
        gw(9)=0.081274388361574
      case default
        write (*,*) 'Integration order too high, reduce NURBS order!!'
      end select
      int=0
      do i=1,l
        do j=1,l
            int=int+1
            w(int)=gw(i)*gw(j)
            r(int)=gp(i)
            z(int)=gp(j)
         end do
      end do
c
      return
      end
c
!      subroutine physicalCoordinates(Xi1,Xi2,ni,nj,x1,x2,x3,ndm)
!c----------------------------------------------------------------------
!c     NOT IN USE AT THE MOMENT
!      USE isogeo
!      implicit double precision (a-h,o-z)
!      real*8 N_Xi1(2,NURp+1),N_Xi2(2,NURq+1),temp(ndm,NURq+1),Sw(ndm)
!c
!      call DersBasisFuns(ni,Xi1,NURp,NURknv1,N_Xi1)
!      call DersBasisFuns(nj,Xi2,NURq,NURknv2,N_Xi2)
!c
!      do indexXi1=0,NURq
!        temp(1:4,indexXi1+1) = 0
!        do indexXi2=0,NURp
!            temp(1:4,indexXi1+1) =temp(1:4,indexXi1+1)+
!     +       N_Xi1(1,indexXi2+1)
!     +      *NURCPw(1:4,ni-NURp+indexXi2,nj-NURq+indexXi1)
!        end do
!      end do
!      Sw = 0
!      do index=1,NURq+1
!        Sw = Sw + N_Xi2(1,index)*temp(1:4,index)
!      end do
!      x1 =Sw(1)/Sw(4)
!      x2 =Sw(2)/Sw(4)
!      x3 =Sw(3)/Sw(4)
!      return
!      end
c
      subroutine physCoor(Xi1,Xi2,ni,nj,x1,x2,x3,ndm,xl,NURp,NURq,
     +                    NURknv1,NURknv2,iflag,v,x4)
c----------------------------------------------------------------------
c     likely to be slower than physicalCoordinates, but works well
c     with local geometry array, deformed geometry is only available
c     in x and xl array, not in NURCPw!
      USE isogeo
      implicit double precision (a-h,o-z)
      real*8 N_Xi1(2,NURp+1),N_Xi2(2,NURq+1)
      dimension xl(ndm,(NURp+1)*(NURq+1)),v((NURp+1)*(NURq+1))
      dimension shp((NURp+1)*(NURq+1))
      real*8 NURknv1(NURlenkv(1)),NURknv2(NURlenkv(2))
c
      call DersBasisFuns(ni,Xi1,NURp,NURknv1,N_Xi1)
      call DersBasisFuns(nj,Xi2,NURq,NURknv2,N_Xi2)
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
      x1=0
      x2=0
      x3=0
      if (iflag.eq.4) x4 = 0
      int=0
      do j=0,NURq
        do i=0,NURp
            int = int+1
            x1 = x1 + shp(int)*xl(1,int)
            x2 = x2 + shp(int)*xl(2,int)
            x3 = x3 + shp(int)*xl(3,int)
            if (iflag.eq.4) x4 = x4 + shp(int)*v(int)
        end do
      end do
      return
      end
c-----------------------------------------------------------------------
      subroutine copyNURstressValues(targval,sourval,NURstress,numel)
      implicit double precision (a-h,o-z)
      integer targval,sourval
      real*8 NURstress(9*numel,11)
      NURstress(:,targval)=NURstress(:,sourval)
      return
      end
c-----------------------------------------------------------------------
      subroutine rprintIGA(dr,ix,x,ndm,numnp,ndf,nen1,idev,numel,
     +                     nen)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     Purpose: compute the profile of values for stre,cont etc.
c              values are in dr(numnp) 
c     Inputs:  
c       dr(ndf,numnp)- array to plot is dr(1,n) n=1..numnp  
c       ix(nen1,*)   - Element nodal connections of mesh
c       x(ndm,numnp) - coordinates
c         ndm        - Spatial dimension of mesh
c         numnp      - Number of nodes in mesh
c         ndf        - 1
c         nen1       - Dimension for ix array
csk         nfl        - no of different values to plot
c         idev       - Graphic type                                
c
c     Outputs:       - rmx, rmn, pr
c
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE contval
      USE iofile
      USE iwinio      
      USE pdata7
      USE rndata
      USE rpdata
      USE isogeo
      implicit double precision (a-h,o-z)
      dimension dr(ndf,*),x(ndm,numnp),pr(9)
      dimension xl_IGA(ndm,nen),v_IGA(nen),val(9*numel) 
      dimension Xi1fani(9),Xi1fanip(9)
      dimension Xi2fanj(9),Xi2fanjp(9) 
      integer ix(nen1,*)
csk      common m(1)
      data blank /-999.d0/
      save  im
      if(icv.eq.1.and.idev.eq.4.and.ior.lt.0)
cww     + call clwopen('Values of    PROFILE',1,iwys-190,630,230,1,2)
     + call clwopen('Values of    PROFILE',1,2)
      call pzero (pr,9)
      drv = 100./numnp


      if(imuse().eq.0) then
c....   plot in the range of all material numbers
        val = 0.0d0
        ival = 0
        do 100 n = 1,numel
c         creating local geometry array, needed for deformed mesh
          v_IGA = 0.0d0        
          do 1009 i = 1,nen
            ii = ix(i,n)
            if(ii.gt.0) then
              do j = 1,ndm
                xl_IGA(j,i) = x(j,ii)
              end do
              v_IGA(i) = dr(1,ii)
            end if
1009      continue
c....     get current patch and corresponding orders      
          NURp = ngetNURnmpq(n,3,AInipa,AInmpq)
          NURq = ngetNURnmpq(n,4,AInipa,AInmpq)
c....     compute NURBS coordinates of recent element      
          ni = ngetNURni(n,AIninc,AInien,AInipa,1)
          nj = ngetNURni(n,AIninc,AInien,AInipa,2)
c....     check if element has zero measure
          if (rNURknv1(ni+1).eq.rNURknv1(ni).or.rNURknv2(nj+1)
     +       .eq.rNURknv2(nj)) goto 100
c         calculate midside and corner points counter-clockwise
          Xi1fani =(/1.d0,0.5d0,0.d0,0.d0,0.d0,0.5d0,1.d0,1.d0,0.5d0/)
          Xi1fanip=(/0.d0,0.5d0,1.d0,1.d0,1.d0,0.5d0,0.d0,0.d0,0.5d0/)
          Xi2fanj =(/1.d0,1.d0,1.d0,0.5d0,0.d0,0.d0,0.d0,0.5d0,0.5d0/)
          Xi2fanjp=(/0.d0,0.d0,0.d0,0.5d0,1.d0,1.d0,1.d0,0.5d0,0.5d0/)
          do i = 1,9
            Xi1 = Xi1fani(i)*rNURknv1(ni)+Xi1fanip(i)*rNURknv1(ni+1)
            Xi2 = Xi2fanj(i)*rNURknv2(nj)+Xi2fanjp(i)*rNURknv2(nj+1)
            call physCoor(Xi1,Xi2,ni,nj,x1,x2,x3,ndm,xl_IGA,
     +                   NURp,NURq,AInkv1,AInkv2,4,v_IGA,x4)
            ival = ival+1
            val(ival) = x4
            if (ival.eq.1) then
              rmx = x4
              rmn = x4
              nmn = 1
              nmx = 1
            end if
            if(rmx.lt.x4) then
              rmx = x4
              nmx = n
            end if
            if(rmn.gt.x4) then
              rmn = x4
              nmn = n
            end if
          end do

100     continue

      else
c....   plot only in the range of specified material numbers(TO BE CHECKED!!!)
        rmx = 0.d0
        rmn = 0.d0
        nmn = 0
        nmx = 0
c....   min,max
        im = 0
        val = 0.0d0
        ival = 0
!        do n = 1,numnp
!          if(iplmano(ix,n,nen1).eq.0)  goto 110
!        end do
          
        do 110 n = 1,numel
c         creating local geometry array, needed for deformed mesh
          v_IGA = 0.0d0        
          do 1008 i = 1,nen
            ii = ix(i,new)
            if(ii.gt.0) then
              do j = 1,ndm
                xl_IGA(j,i) = x(j,ii)
              end do
              v_IGA(i) = dr(1,ii)
            end if
1008      continue
c....     get current patch and corresponding orders      
          NURp = ngetNURnmpq(n,3,AInipa,AInmpq)
          NURq = ngetNURnmpq(n,4,AInipa,AInmpq)
c....     compute NURBS coordinates of recent element      
          ni = ngetNURni(n,AIninc,AInien,AInipa,1)
          nj = ngetNURni(n,AIninc,AInien,AInipa,2)
c....     check if element has zero measure
          if (rNURknv1(ni+1).eq.rNURknv1(ni).or.rNURknv2(nj+1)
     +       .eq.rNURknv2(nj)) goto 110
c         calculate midside and corner points counter-clockwise
          Xi1fani =(/1.d0,0.5d0,0.d0,0.d0,0.d0,0.5d0,1.d0,1.d0,0.5d0/)
          Xi1fanip=(/0.d0,0.5d0,1.d0,1.d0,1.d0,0.5d0,0.d0,0.d0,0.5d0/)
          Xi2fanj =(/1.d0,1.d0,1.d0,0.5d0,0.d0,0.d0,0.d0,0.5d0,0.5d0/)
          Xi2fanjp=(/0.d0,0.d0,0.d0,0.5d0,1.d0,1.d0,1.d0,0.5d0,0.5d0/)

          do i = 1,9
            Xi1 = Xi1fani(i)*rNURknv1(ni)+Xi1fanip(i)*rNURknv1(ni+1)
            Xi2 = Xi2fanj(i)*rNURknv2(nj)+Xi2fanjp(i)*rNURknv2(nj+1)
            call physCoor(Xi1,Xi2,ni,nj,x1,x2,x3,ndm,xl_IGA,
     +                   NURp,NURq,AInkv1,AInkv2,4,v_IGA,x4)
            ival = ival+1
            val(ival) = x4
            if(rmx.lt.x4) then
              rmx = x4
              nmx = n
            end if
            if(rmn.gt.x4) then
              rmn = x4
              nmn = n
            end if
          end do
110     continue
      end if

      drm =  dabs(rmx-rmn)
      if(drm.lt.1.d-5*dabs(rmx).or.drm.lt.1.e-10) then
c....   nearly same values, no profile and plot, cont-->fill
          if(icv.eq.1) then
            if(ior.ge.0) write(iow,2001) rmn,rmx
            if(ior.lt.0) write(*  ,2001) rmn,rmx
          end if 
cww          nfl = 0
      else
c....   profile 
        do 200 n = 1,9*numel
          rs = (val(n) - rmn)/(rmx - rmn)
          do i = 1,9
            if(rs.ge.0.1d0*i) pr(i) = pr(i) + drv
          end do
200     continue
        if(icv.eq.1) then
          if(ior.ge.0) write(iow,2000) rmn,nmn,rmx,nmx,pr
          if(ior.lt.0) write(*  ,2000) rmn,nmn,rmx,nmx,pr
        end if
      end if
      return
2000  format(' Minimum is ',e15.5,' in element ',i6/,
     1       ' Maximum is ',e15.5,' in element ',i6,/,
     2  20x,'10%   20%   30%   40%   50%   60%   70%   80%   90%'/
     3       ' Profile above is:',10f6.1)
2001  format(' WARNING no profile, Minimum ',e15.5,' Maximum ',e15.5)
      end
c
c
      !subroutine pdireciga1(idtyp,x,dir,numnp,ndm,ifdir)
c-----------------------------------------------------------------------
c
c     Purpose: computes the director for all global nodes/control points
c
c     Inputs:  idtyp: 18: director averaged from linear control point mesh
c              x    : control points
c              dir  : director field storage
c              dir(10,knode,2) - triad field 1-9 = dir, 10 = no.of el.per node
c              dir( , ,1)      - save
c              dir( , ,2)      - test, update at TIME
c              ifdir           - 1=save 2=test, update at TIME
c
c     Outputs: -
c
c
c-----------------------------------------------------------------------
!      USE dirdat
!      USE iofile
!      USE isogeo
!      implicit double precision (a-h,o-z)
!      dimension dir(10,numnp,2),x(ndm,numnp)
!      dimension vec1(3,4),vec2(3,4),vec3(3,4)
!      integer nNURnmpq,patch
!      node = 0
!      do patch = 1,NURnpatch
!        n = nNURnmpq(patch,1)
!        m = nNURnmpq(patch,2)
!        do im=1,m
!          do in=1,n
!          node = node + 1
!            if ((im.eq.1.or.im.eq.m).and.(in.eq.1.or.in.eq.n)) then
!c           corner control point
!              nneighbor = 1
!              if (im.eq.1) nim = n
!              if (im.eq.m) nim =-n
!              if (im.eq.1) sim = 1
!              if (im.eq.m) sim =-1
!              if(in.eq.1) nin = 1
!              if(in.eq.n) nin =-1
!              if (idtyp.eq.18) then
!                vec1(1:3,1) = nin*(x(1:3,node+nin)-x(1:3,node))
!     +                        /x(4,node+nin)/x(4,node)
!                vec2(1:3,1) = sim*(x(1:3,node+nim)-x(1:3,node))
!     +                        /x(4,node+nim)/x(4,node)
!              else if (idtyp.eq.19) then
!                vec1(1:3,1) = nin*sim*(x(1:3,node+nim+nin)-x(1:3,node))
!                vec2(1:3,1) = sim*nin*(x(1:3,node+nin)-x(1:3,node+nim))
!              end if
!            else if (im.eq.1.or.im.eq.m) then
!c           edge control point direction 1
!              nneighbor = 2
!              if (im.eq.1) nim = n
!              if (im.eq.m) nim =-n
!              if (im.eq.1) sim = 1
!              if (im.eq.m) sim =-1
!              if (idtyp.eq.18) then
!                vec1(1:3,1) = (x(1:3,node+1)-x(1:3,node))
!     +                        /x(4,node+1)/x(4,node)
!                vec2(1:3,1) = sim*(x(1:3,node+nim)-x(1:3,node))
!     +                        /x(4,node+nim)/x(4,node)
!                vec1(1:3,2) = (x(1:3,node)-x(1:3,node-1))
!     +                        /x(4,node)/x(4,node-1)
!                vec2(1:3,2) = sim*(x(1:3,node+nim)-x(1:3,node))
!     +                        /x(4,node+nim)/x(4,node)
!              else if (idtyp.eq.19) then
!                vec1(1:3,1) = sim*(x(1:3,node+nim+1)-x(1:3,node))
!                vec2(1:3,1) = sim*(x(1:3,node+1)-x(1:3,node+nim))
!                vec1(1:3,2) = sim*(x(1:3,node+nim)-x(1:3,node-1))
!                vec2(1:3,2) = sim*(x(1:3,node)-x(1:3,node+nim-1))
!              end if
!            else if (in.eq.1.or.in.eq.n) then
!c           edge control point direction 2
!              nneighbor = 2
!              if (in.eq.1) nin = 1
!              if (in.eq.n) nin =-1
!              if (idtyp.eq.18) then
!                vec1(1:3,1) = nin*(x(1:3,node+nin)-x(1:3,node))
!     +                        /x(4,node+nin)/x(4,node)
!                vec2(1:3,1) = (x(1:3,node+n)-x(1:3,node))
!     +                        /x(4,node+n)/x(4,node)
!                vec1(1:3,2) = nin*(x(1:3,node+nin)-x(1:3,node))
!     +                        /x(4,node+nin)/x(4,node)
!                vec2(1:3,2) = (x(1:3,node)-x(1:3,node-n))
!     +                        /x(4,node)/x(4,node-n)
!              else if (idtyp.eq.19) then
!                vec1(1:3,1) = nin*(x(1:3,node+nin+n)-x(1:3,node))
!                vec2(1:3,1) = nin*(x(1:3,node+nin)-x(1:3,node+n))
!                vec1(1:3,2) = nin*(x(1:3,node+nin)-x(1:3,node-n))
!                vec2(1:3,2) = nin*(x(1:3,node+nin-n)-x(1:3,node))
!              end if
!            else
!c           regular control point
!              nneighbor = 4
!              if (idtyp.eq.18) then
!                vec1(1:3,1) = (x(1:3,node+1)-x(1:3,node))
!     +                        /x(4,node+1)/x(4,node)
!                vec2(1:3,1) = (x(1:3,node+n)-x(1:3,node))
!     +                        /x(4,node+n)/x(4,node)
!                vec1(1:3,2) = (x(1:3,node)-x(1:3,node-1))
!     +                        /x(4,node)/x(4,node-1)
!                vec2(1:3,2) = (x(1:3,node+n)-x(1:3,node))
!     +                        /x(4,node+n)/x(4,node)
!                vec1(1:3,3) = (x(1:3,node)-x(1:3,node-1))
!     +                        /x(4,node)/x(4,node-1)
!                vec2(1:3,3) = (x(1:3,node)-x(1:3,node-n))
!     +                        /x(4,node)/x(4,node-n)
!                vec1(1:3,4) = (x(1:3,node+1)-x(1:3,node))
!     +                        /x(4,node+1)/x(4,node)
!                vec2(1:3,4) = (x(1:3,node)-x(1:3,node-n))
!     +                        /x(4,node)/x(4,node-n)
!              else if (idtyp.eq.19) then
!                vec1(1:3,1) = x(1:3,node+1+n)-x(1:3,node)
!                vec2(1:3,1) = x(1:3,node+1)-x(1:3,node+n)
!                vec1(1:3,2) = x(1:3,node+n)-x(1:3,node-1)
!                vec2(1:3,2) = x(1:3,node)-x(1:3,node+n-1)
!                vec1(1:3,3) = x(1:3,node)-x(1:3,node-1-n)
!                vec2(1:3,3) = x(1:3,node-n)-x(1:3,node-1)
!                vec1(1:3,4) = x(1:3,node+1)-x(1:3,node-n)
!                vec2(1:3,4) = x(1:3,node-n+1)-x(1:3,node)
!              end if
!            end if
!            if (idtyp.eq.19) then
!            do i =1,nneighbor
!              call vecp (vec1(1:3,i),vec2(1:3,i),vec3(1:3,i))
!              if((vec3(1,i)+vec3(2,i)+vec3(3,i)).eq.0.0d0) then
!                write(iow,*) 'zero director, try n1=20 instead'
!                vec1(1,i)=1.0d0;vec1(2,i)=0.0d0;vec1(3,i)=0.0d0
!                vec2(1,i)=0.0d0;vec2(2,i)=1.0d0;vec2(3,i)=0.0d0
!                vec3(1,i)=0.0d0;vec3(2,i)=0.0d0;vec3(3,i)=1.0d0
!c                write(*,*) 'zero director due to control point coinciden
!c     +ce: results are likely to be corrupt, try n1=19 instead'
!              end if
!              call norm (vec3(1:3,i),vec3(1:3,i),3)
!              call norm (vec1(1:3,i),vec1(1:3,i),3)
!              call vecp (vec3(1:3,i),vec1(1:3,i),vec2(1:3,i))
!
!              do k = 1,3
!                dir(k  ,node,ifdir) = dir(k  ,node,ifdir) + vec1(k,i)
!c     +                                               /nneighbor
!                dir(k+3,node,ifdir) = dir(k+3,node,ifdir) + vec2(k,i)
!c     +                                               /nneighbor
!                dir(k+6,node,ifdir) = dir(k+6,node,ifdir) + vec3(k,i)
!c     +                                               /nneighbor
!              end do
!              dir(10,node,ifdir) = dir(10,node,ifdir) + 1
!            end do
!            end if
!            if(idtyp.eq.18) then
!c            call norm (vec2(1:3,1),vec2(1:3,1),3)
!c            call norm (vec1(1:3,1),vec1(1:3,1),3)
!c            vec1(1:3,1)=vec1(1:3,1)!*x(4,node)**nneighbor
!c            vec2(1:3,1)=vec2(1:3,1)!*x(4,node)**nneighbor
!            do i=2,nneighbor
!c              call norm (vec2(1:3,i),vec2(1:3,i),3)
!c              call norm (vec1(1:3,i),vec1(1:3,i),3)
!              vec1(1:3,1)=vec1(1:3,1)+vec1(1:3,i)!*x(4,node)**nneighbor
!              vec2(1:3,1)=vec2(1:3,1)+vec2(1:3,i)!*x(4,node)**nneighbor
!            end do
!            call norm (vec2(1:3,1),vec2(1:3,1),3)
!            call norm (vec1(1:3,1),vec1(1:3,1),3)
!            call vecp (vec1(1:3,1),vec2(1:3,1),vec3(1:3,1))
!            if((vec3(1,i)+vec3(2,i)+vec3(3,i)).eq.0.0d0) then
!                write(iow,*) 'zero director due to control point coincid
!     +ence: results are likely to be corrupt, try n1=19 instead'
!                vec1(1,1)=1.0d0;vec1(2,1)=0.0d0;vec1(3,1)=0.0d0
!                vec2(1,1)=0.0d0;vec2(2,1)=1.0d0;vec2(3,1)=0.0d0
!                vec3(1,1)=0.0d0;vec3(2,1)=0.0d0;vec3(3,1)=1.0d0
!            end if
!            call norm (vec3(1:3,1),vec3(1:3,1),3)
!            do k = 1,3
!                dir(k  ,node,ifdir) =  vec1(k,1)
!                dir(k+3,node,ifdir) =  vec2(k,1)
!                dir(k+6,node,ifdir) =  vec3(k,1)
!            end do
!            dir(10,node,ifdir) = dir(10,node,ifdir) + 1
!            end if
!          end do
!        end do
!      end do
!      return
!      end
!
!c-----------------------------------------------------------------------
!      subroutine pdireciga2(idtyp,x,dir,numnp,ndm,ifdir,tielist)
!c-----------------------------------------------------------------------
!c
!c     Purpose: computes the director for all global nodes/control points
!c              realized through the normed tangent to the NURBS lamina
!c              in closest surface point to the control point
!c              This requires a local Newton-Raphson iteration
!c              see Piegl/Tiller NURBS book p. 232
!c
!c     Inputs:  idtyp: 20: director from normal to shell lamina
!c              x    : global control points
!c              dir  : director field storage
!c              dir(10,knode,2) - triad field 1-9 = dir, 10 = no.of el.per node
!c              dir( , ,1)      - save
!c              dir( , ,2)      - test, update at TIME
!c              ifdir           - 1=save 2=test, update at TIME
!c
!c     Outputs: -
!c
!c
!c-----------------------------------------------------------------------
!      USE dirdat
!      USE iofile
!      USE isogeo
!
!c     input/output variable declaration
!      real*8 dir(10,numnp,2),x(ndm,*)
!      integer numnp,ndm,ifdir,idtyp
!c     functions needed within subroutine
!      real*8,pointer,dimension(:)    :: KV1,KV2
!      real*8,pointer,dimension(:,:,:):: CPw
!      integer p,q,n,m,sum_in,sum_im,patch,im,in,node,count,line,k
!      integer converged,nNURnmpq
!      integer tielist(numnp)
!      real*8 point(3),SKL(4,3,3),Aders(3,3,3),wders(1,3,3),r(3)
!      real*8 Xi1_min,Xi1_max,Xi2_min,Xi2_max,Xi1,Xi2,Xi1n,Xi2n
!      real*8 fuv,guv,fu,fv,gu,gv,det,inv11,inv22,inv12,inv21
!      real*8 eps1,eps2,pcoin,acoin1,acoin2,pchange(3)
!      real*8 rNURknv1,rNURknv2
!      double precision pcoin1
!      real*8 vec1(3),vec2(3),vec3(3),SmP(3)
!      real*8 dot
!      eps1 = 1.0d-13 !tolerance for point coincidence
!      eps2 = 1.0d-13 !zero cosine measure
!      sum_in = 0
!      sum_im = 0
!      node = 0
!      line = 0
!c     do for all NURBS patches
!      do patch = 1,NURnpatch
!        n = nNURnmpq(patch,1)
!        m = nNURnmpq(patch,2)
!        p = nNURnmpq(patch,3)
!        q = nNURnmpq(patch,4)
!        allocate(KV1(1:n+p+1))
!        allocate(KV2(1:m+q+1))
!        allocate(CPw(1:4,1:n,1:m))
!        do in=1,n+p+1
!          KV1(in) = rNURknv1(sum_in+in)
!        end do
!        do im=1,m+q+1
!          KV2(im) = rNURknv2(sum_im+im)
!        end do
!c       create CPw
!        do im=1,m
!          do in=1,n
!            line = line + 1
!            CPw(1:3,in,im) = x(1:3,tielist(line))*x(4,tielist(line))
!            CPw(4  ,in,im) = x(4,tielist(line))
!          end do
!        end do
!        do im=1,m
!          do in=1,n
!            node = node + 1
!            converged = 0
!            point(1:3) = x(1:3,node)
!c           utilize convex hull properties for start value
!            Xi1_min = rNURknv1(in+sum_in)
!            Xi1_max = rNURknv1(n+p+1+sum_in)
!            Xi2_min = rNURknv2(im+sum_im)
!            Xi2_max = rNURknv2(m+q+1+sum_im)
!c           choose start value
!            Xi1 = (Xi1_min+Xi1_max)/2.0d0
!            Xi2 = (Xi2_min+Xi2_max)/2.0d0
!            count=0
!            do while (converged.eq.0)
!              count = count + 1
!c             compute value, 1st and 2nd derivative
!              call SurfaceDerivsAlg1(n,p,KV1,m,q,KV2,CPw,Xi1,Xi2,d,SKL)
!              Aders(1:3,1:3,1:3) = SKL(1:3,1:3,1:3)
!              wders(1  ,1:3,1:3) = SKL(4  ,1:3,1:3)
!              call RatSurfaceDerivs(Aders,wders,SKL)
!c             check convergence criteria:
!              SmP = SKL(1:3,1,1)-point(1:3)
!              pcoin = sqrt(dot(SmP,SmP,3))
!              acoin1= sqrt(dot(SKL(1:3,2,1),SmP,3))/
!     +                sqrt(dot(SKL(1:3,2,1),SKL(1:3,2,1),3))/
!     +                sqrt(dot(SmP,SmP,3))
!              acoin2= sqrt(dot(SKL(1:3,1,2),SmP,3))/
!     +                sqrt(dot(SKL(1:3,1,2),SKL(1:3,1,2),3))/
!     +                sqrt(dot(SmP,SmP,3))
!c             1st: point coincidence
!              if (pcoin.le.eps1) then
!                converged = 1
!c             2nd: zero cosine
!              else if (acoin1.le.eps2.and.acoin2.le.eps2) then
!                converged = 1
!c             3rd: compute next values
!              else
!c               compute entries for system of equations
!                r = SKL(1:3,1,1) - point(1:3)
!                fuv = dot(r,SKL(1:3,2,1),3)
!                guv = dot(r,SKL(1:3,1,2),3)
!                fu  = dot(SKL(1:3,2,1),SKL(1:3,2,1),3)
!     +            + dot(r,SKL(1:3,3,1),3)
!                fv  = dot(SKL(1:3,2,1),SKL(1:3,1,2),3)
!     +            + dot(r,SKL(1:3,2,2),3)
!                gu  = dot(SKL(1:3,2,1),SKL(1:3,1,2),3)
!     +            + dot(r,SKL(1:3,2,2),3)
!                gv  = dot(SKL(1:3,1,2),SKL(1:3,1,2),3)
!     +            + dot(r,SKL(1:3,1,3),3)
!c               inversion
!                det = fu*gv-fv*gu
!                inv11 = gv/det
!                inv22 = fu/det
!                inv12 =-fv/det
!                inv21 =-gu/det
!c               compute next values
!                Xi1n = Xi1 - inv11*fuv - inv12*guv
!                Xi2n = Xi2 - inv21*fuv - inv22*guv
!c               check if next value is a possible value (perhaps use KV(1) to KV(end)???
!                if (Xi1n.le.Xi1_min) Xi1n=Xi1_min
!                if (Xi1n.ge.Xi1_max) Xi1n=Xi1_max
!                if (Xi2n.le.Xi2_min) Xi2n=Xi2_min
!                if (Xi2n.ge.Xi2_max) Xi2n=Xi2_max
!c                if (isnan(Xi1n)) Xi1n=Xi1
!c                if (isnan(Xi2n)) Xi2n=Xi2
!c              4th: check for change of parameters
!                pchange(1:3)=(Xi1n-Xi1)*SKL(1:3,2,1)
!     +                      +(Xi2n-Xi2)*SKL(1:3,1,2)
!                if (sqrt(dot(pchange,pchange,3)).le.(eps1*100)) then
!                  converged = 1
!                else if (count.ge.25) then
!                  converged = 1
!                  write (iow,*) 'NO LOCAL CONVERGENCE FOR DIRECTOR'
!                  write (*,*)   'NO LOCAL CONVERGENCE FOR DIRECTOR'
!                else
!                  Xi1 = Xi1n
!                  Xi2 = Xi2n
!                end if
!              end if
!            end do
!c           create director
!            vec1(1:3) = SKL(1:3,2,1)
!            vec2(1:3) = SKL(1:3,1,2)
!
!            call norm (vec2(1:3),vec2(1:3),3)
!            call norm (vec1(1:3),vec1(1:3),3)
!            call vecp(vec1(1:3),vec2(1:3),vec3(1:3))
!            if((vec3(1)+vec3(2)+vec3(3)).eq.0.0d0) then
!                write(iow,*) 'zero director, panic!'
!                vec1(1)=1.0d0;vec1(2)=0.0d0;vec1(3)=0.0d0
!                vec2(1)=0.0d0;vec2(2)=1.0d0;vec2(3)=0.0d0
!                vec3(1)=0.0d0;vec3(2)=0.0d0;vec3(3)=1.0d0
!            end if
!            call norm (vec3(1:3),vec3(1:3),3)
!            do k = 1,3
!                dir(k  ,node,ifdir) =  vec1(k)
!                dir(k+3,node,ifdir) =  vec2(k)
!                dir(k+6,node,ifdir) =  vec3(k)
!            end do
!            dir(10,node,ifdir) = dir(10,node,ifdir) + 1
!c           control
!c            inv11 = SKL(1,1,1)+5*vec3(1)
!c            inv22 = SKL(2,1,1)+5*vec3(2)
!c            write (iow,*) inv11,inv22
!          end do
!        end do
!        sum_in = sum_in + n+p+1
!        sum_im = sum_im + m+q+1
!        deallocate(KV1)
!        deallocate(KV2)
!        deallocate(CPw)
!      end do
!
!      return
!      end

c-----------------------------------------------------------------------
      subroutine pdireciga5(inumgp,x,dir,numnp,ndm,ifdir,numel,
     +           INC,IEN,IPA,nenglobal,prt,ipr)
c-----------------------------------------------------------------------
c
c     THIS IS THE FUNCTION TO COMPUTE THE NODAL BASIS SYSTEMS      
c     Purpose: computes the director for all global nodes/control points
c              realized through normal equation of an overdetermined equation
c              system
c     
c     Inputs:  idtyp: 21
c              x    : global control points
c              dir  : director field storage  
c              dir(10,knode,2) - triad field 1-9 = dir, 10 = no.of el.per node
c              dir( , ,1)      - save   
c              dir( , ,2)      - test, update at TIME
c              ifdir           - 1=save 2=test, update at TIME
c
c     Outputs: -
c
c    
c-----------------------------------------------------------------------
      USE dirdat
      USE isogeo
      USE iofile
      USE doalloc
c     input/output variable declaration
      real*8 dir(10,numnp,2),x(ndm,*)
      integer numnp,ndm,ifdir,inumgp,nenglobal,numel,ipr
      integer INC(numnp,2),IEN(nenglobal,numel),IPA(numel,3)
      logical prt,flbasissys
c      integer tielist(numnp)  
c     variables needed within subroutine
      integer patch,n,m,p,q,im,in,sum_in,sum_im,line,nnzel,nel,
     +        ngp,ngp2,numlines,lint,nze,ielem,ni,nj,l,nen,i,j,
     +        sum_cp,sum_el,eq_i,eq_j,flag,node,element
      real*8 Xi1,Xi2
      real*8 SKL(4,3,3),Aders(3,3,3),wders(1,3,3),blocal(9)
      real*8,pointer,dimension(:)    :: KV1,KV2,sg,tg,wg,shp
      real*8,pointer,dimension(:,:)  :: ntnlocal,ntblocal,ntn,ntb
      real*8,pointer,dimension(:,:,:):: CPw
      integer,pointer,dimension(:)   :: ipiv
c     functions called within
      real*8 rNURknv1,rNURknv2
      integer ngetNURni
c     end of variable declaration 
      sum_in = 0
      sum_im = 0
      line   = 0
      sum_cp = 0
      sum_el = 0
      node   = 0
c...  for inumgp.eq.0: Reference director D is interpolated    
      if (inumgp.eq.0) then
        
        Iexactnormal = 0
        if(prt) write (iow,*) 
     +   '    Ref. Director is interpolated with nodal directors'
        if(prt) write (*,*) 
     +   '    Ref. Director is interpolated with nodal directors'
      
      else if (inumgp.ge.1) then
        
        Iexactnormal = 1
        ngp2 = inumgp**2
        InumAreaGP = inumgp
        ngp = inumgp
        if(prt) write (iow,*) 
     +   '    Ref. Director is computed exactly from normal vector'
        if(prt) write (iow,*) 
     +   ngp2,'    Integration points per element'      
        if(prt) write (iow,*) 
     +   '    This overrides choice in elem!'
        if(prt) write (*,*) 
     +   '    Ref. Director is computed exactly from normal vector'
        if(prt) write (*,*) 
     +   inumgp**2,'    Integration points per element'      
        if(prt) write (*,*) 
     +   '    This overrides choice in elem!'
      
c        allocate(IGAbasissystems(1:3,1:3,1:3,1:ngp2,1:numel))
c        call pzero(IGAbasissystems,27*ngp2*numel)
        call ralloc(AInbasissys,27*numel*ngp2,'IGA basis',flbasissys)
c        call pseta(Inbasissys,27*numel*ngp2,ipr,flbasissys,'IGA basis')
c        call IGAsetZeroNormalm(AInbasissys,27*numel*ngp2)
        
        
      else
      
        write (iow,*) '    NO APPROPRIATE PARAMETER IN BASE,20,X'
        write (*  ,*) '    NO APPROPRIATE PARAMETER IN BASE,20,X'
      
      end if
c     do for all NURBS patches        
      do patch = 1,NURnpatch
        write (iow,*) 'patch: ',patch
        n = nNURnmpq(patch,1)
        m = nNURnmpq(patch,2)
        p = nNURnmpq(patch,3)
        q = nNURnmpq(patch,4)
        allocate(KV1(1:n+p+1))
        allocate(KV2(1:m+q+1))
        allocate(CPw(1:4,1:n,1:m))
        do in=1,n+p+1
          KV1(in) = rNURknv1(sum_in+in)
        end do
        do im=1,m+q+1
          KV2(im) = rNURknv2(sum_im+im)
        end do
c       create CPw
        do im=1,m
          do in=1,n
            line = line + 1
            CPw(1:3,in,im) = x(1:3,line)*x(4,line)
c            CPw(1:3,in,im) = x(1:3,tielist(line))*x(4,line)
            CPw(4  ,in,im) = x(4,line)
          end do
        end do
c       determine number of non-zero-sized elements nnzel and
c       number of elements nel in current patch
        call nnonzeroelements(KV1,KV2,n,m,p,q,nnzel,nel)
        if (Iexactnormal.eq.0) then
c         determine number of Gauss Points and Gauss weights
          ngp = max(p,q) + 1
c         ngp = min(ngp,9)
          ngp2= ngp**2
        end if  
        numlines = ngp2*nnzel
        nen = (p+1) * (q+1)
        allocate(sg(1:ngp2))
        allocate(tg(1:ngp2))
        allocate(wg(1:ngp2))
        allocate(shp(1:nen))
        allocate(ntnlocal(1:nen,1:nen))
        allocate(ntblocal(1:nen,9))
        allocate(ntn(1:n*m,1:n*m))
        allocate(ntb(1:n*m,1:9))
        call pzero(ntn,n*n*m*m)
        call pzero(ntb,n*m*9)
        call IGAgauss (ngp,lint,sg,tg,wg)
c       counter for non-zero-sized elements
        nze = 0
c       loop over all elements, including zero-sized
        do 1001 ielem = 1,nel
c         compute NURBS coordinates of recent element      
          ni = ngetNURni(ielem+sum_el,INC,IEN,IPA,1)-sum_in
          nj = ngetNURni(ielem+sum_el,INC,IEN,IPA,2)-sum_im
c         check if element has zero measure
          if (rNURknv1(ni).eq.rNURknv1(ni+1).or.
     +        rNURknv2(nj).eq.rNURknv2(nj+1)) goto 1001
          nze = nze + 1
          write (iow,*) ielem+sum_el
          do l = 1,lint
c           calculate NURBS parametric coordinates from unit element coordinates
            Xi1=0.5d0*((KV1(ni+1)-KV1(ni))*sg(l)+(KV1(ni+1)+KV1(ni)))
            Xi2=0.5d0*((KV2(nj+1)-KV2(nj))*tg(l)+(KV2(nj+1)+KV2(nj))) 
c           compute physical coordinates, 1st and 2nd derivative 
            call SurfaceDerivsAlg1(n,p,KV1,m,q,KV2,CPw,Xi1,Xi2,d,SKL)
            Aders(1:3,1:3,1:3) = SKL(1:3,1:3,1:3)
            wders(1  ,1:3,1:3) = SKL(4  ,1:3,1:3)
            call RatSurfaceDerivs(Aders,wders,SKL)
c            write (iow,*) SKL(1:3,1,1)
c           compute normal and store it in SKL(1:3,1,1), norm and save to blocal 
            call vecp(SKL(1:3,2,1),SKL(1:3,1,2),SKL(1:3,1,1))
            call norm(blocal(1:3),SKL(1:3,2,1),3)
            call norm(blocal(4:6),SKL(1:3,1,2),3)
            call norm(blocal(7:9),SKL(1:3,1,1),3)
c           construct lamina basis system
            call plamina_base(blocal)
            if (Iexactnormal.eq.1) then
c             store normal and basis system in IGAbasissystems              
c              call IGAcompStorNormal(SKL,blocal,l,ielem+sum_el)
              element = ielem+sum_el
              call IGAcompStorNormalm1(SKL,blocal,l,element,
     +                                AInbasissys,numel,ngp2)
            end if
c           compute basis functions of recent Gauss Point, into shp
            call IGAshapeDir(Xi1,Xi2,shp,ni,nj,CPw,p,q,KV1,KV2,nen,n,m)
            call mttmul(shp,shp,nen,1,nen,ntnlocal)
            call mttmul(shp,blocal,nen,1,9,ntblocal)
c           assembly to patch system of equation
            do i = 1,nen
              eq_i = IEN(i,ielem+sum_el) - sum_cp
              do j = 1,nen
                eq_j = IEN(j,ielem+sum_el) - sum_cp
                ntn(eq_i,eq_j) = ntn(eq_i,eq_j) + ntnlocal(i,j)
              end do
              do j = 1,9
                ntb(eq_i,j) = ntb(eq_i,j) + ntblocal(i,j)
              end do
            end do
          end do
1001    continue
c       solve system with INTEL MKL library
        allocate(ipiv(1:n*m))
        call dgesv(n*m,9,ntn,n*m,ipiv,ntb,n*m,flag)
        if (flag.ne.0) write(*,*) '!!!  error in MKL solver  !!!'
        deallocate(ipiv)
c       estimate condition number
c        call dgecon ('1',n*m
c       save values to dir array
        do i = 1,n*m
          node = node + 1
          do j = 1,9
            dir(j,node,ifdir) = ntb(i,j)
          end do
          dir(10,node,ifdir) = 1
c          call plamina_base(dir(1:9,node,ifdir))
c          write(iow,*) dir(7:9,node,ifdir)
        end do        
c       end current patch: deallocate all dynamic arrays
        deallocate(KV1,KV2,CPw,sg,tg,wg,shp,ntnlocal,ntblocal)
        deallocate(ntn,ntb)
        sum_cp = sum_cp + n*m
        sum_el = sum_el + nel
        sum_in = sum_in + n+p+1
        sum_im = sum_im + m+q+1
      end do
      if(prt) write(iow,2000) 
      if(prt) write(*  ,2000)
2000  format(/ '     Nodal basis system for NURBS surfaces determined by
     1 normal equation'/)     
      return
        end
c
      subroutine pdireciga6(idtyp,x,dir,numnp,ndm,ifdir,numel,prt)
c-----------------------------------------------------------------------
c
c     Purpose: computes the director for all global nodes/control points
c              realized through the normed tangent to the NURBS lamina
c              in closest surface point to the control point
c              This requires a local Newton-Raphson iteration
c              see Piegl/Tiller NURBS book p. 232
c     
c     Inputs:  idtyp: 20: director from normal to shell lamina
c              x    : global control points
c              dir  : director field storage  
c              dir(10,knode,2) - triad field 1-9 = dir, 10 = no.of el.per node
c              dir( , ,1)      - save   
c              dir( , ,2)      - test, update at TIME
c              ifdir           - 1=save 2=test, update at TIME
c
c     Outputs: -
c
c    
c-----------------------------------------------------------------------
      
      USE dirdat
      USE isogeo
      USE iofile
c     input/output variable declaration
      real*8 dir(10,numnp,2),x(ndm,*)
      integer numnp,ndm,ifdir,idtyp
      logical prt
c     functions needed within subroutine      
      real*8,pointer,dimension(:)    :: KV1,KV2
      real*8,pointer,dimension(:,:,:):: CPw
      integer p,q,n,m,sum_in,sum_im,patch,im,in,node,count,line,k
      integer converged,nNURnmpq
csk      integer tielist(numnp)
      real*8 point(3),SKL(4,3,3),Aders(3,3,3),wders(1,3,3),r(3)
      real*8 Xi1_min,Xi1_max,Xi2_min,Xi2_max,Xi1,Xi2,Xi1n,Xi2n
      real*8 fuv,guv,fu,fv,gu,gv,det,inv11,inv22,inv12,inv21
      real*8 eps1,eps2,pcoin,acoin1,acoin2,pchange(3)
      real*8 rNURknv1,rNURknv2,gs(3,2)
csk      double precision pcoin1
      real*8 vec1(3),vec2(3),vec3(3),SmP(3)
      real*8 dot
c     end of variable declaration
c
      eps1 = 1.0d-12 !tolerance for point coincidence
      eps2 = 1.0d-12 !zero cosine measure
      sum_in = 0
      sum_im = 0
      node = 0
      line = 0
c     do for all NURBS patches 
c
      do patch = 1,NURnpatch
        n = nNURnmpq(patch,1)
        m = nNURnmpq(patch,2)
        p = nNURnmpq(patch,3)
        q = nNURnmpq(patch,4)
        allocate(KV1(1:n+p+1))
        allocate(KV2(1:m+q+1))
        allocate(CPw(1:4,1:n,1:m))
        do in=1,n+p+1
          KV1(in) = rNURknv1(sum_in+in)
        end do
        do im=1,m+q+1
          KV2(im) = rNURknv2(sum_im+im)
        end do
c       create CPw
        do im=1,m
          do in=1,n
            line = line + 1
            CPw(1:3,in,im) = x(1:3,line)*x(4,line)
            CPw(4  ,in,im) = x(4,line)
          end do
        end do
        edgenode = 0
        do im=1,m
          do in=1,n
            node = node + 1
            converged = 0
            point(1:3) = x(1:3,node)
c           utilize convex hull properties for start value               
            Xi1_min = rNURknv1(in+sum_in)
            Xi1_max = rNURknv1(in+p+1+sum_in)
            Xi2_min = rNURknv2(im+sum_im)
            Xi2_max = rNURknv2(im+q+1+sum_im)
c           choose start value
            Xi1 = (Xi1_min+Xi1_max)/2.0d0
            Xi2 = (Xi2_min+Xi2_max)/2.0d0
            count=0
            do while (converged.eq.0)
              count = count + 1
c             compute value, 1st and 2nd derivative
              call SurfaceDerivsAlg1(n,p,KV1,m,q,KV2,CPw,Xi1,Xi2,d,SKL)
              Aders(1:3,1:3,1:3) = SKL(1:3,1:3,1:3)
              wders(1  ,1:3,1:3) = SKL(4  ,1:3,1:3)
              call RatSurfaceDerivs(Aders,wders,SKL)
c             check convergence criteria:
              SmP = SKL(1:3,1,1)-point(1:3)
              pcoin = sqrt(dot(SmP,SmP,3))
              acoin1= sqrt(dot(SKL(1:3,2,1),SmP,3))/
     +                sqrt(dot(SKL(1:3,2,1),SKL(1:3,2,1),3))/
     +                sqrt(dot(SmP,SmP,3))
              acoin2= sqrt(dot(SKL(1:3,1,2),SmP,3))/
     +                sqrt(dot(SKL(1:3,1,2),SKL(1:3,1,2),3))/
     +                sqrt(dot(SmP,SmP,3)) 
c             1st: point coincidence
              if (pcoin.le.eps1) then
                converged = 1
c             2nd: zero cosine
              else if (acoin1.le.eps2.and.acoin2.le.eps2) then
                converged = 1
c             3rd: compute next values
              else
c               compute entries for system of equations
                r = SKL(1:3,1,1) - point(1:3)
                fuv = dot(r,SKL(1:3,2,1),3)
                guv = dot(r,SKL(1:3,1,2),3)
                fu  = dot(SKL(1:3,2,1),SKL(1:3,2,1),3)
     +            + dot(r,SKL(1:3,3,1),3)   
                fv  = dot(SKL(1:3,2,1),SKL(1:3,1,2),3)
     +            + dot(r,SKL(1:3,2,2),3)
                gu  = dot(SKL(1:3,2,1),SKL(1:3,1,2),3)
     +            + dot(r,SKL(1:3,2,2),3)   
                gv  = dot(SKL(1:3,1,2),SKL(1:3,1,2),3)
     +            + dot(r,SKL(1:3,1,3),3)
c               inversion
                det = fu*gv-fv*gu
                inv11 = gv/det
                inv22 = fu/det
                inv12 =-fv/det
                inv21 =-gu/det
c               compute next values
                Xi1n = Xi1 - inv11*fuv - inv12*guv
                Xi2n = Xi2 - inv21*fuv - inv22*guv
c               check if next value is a possible value (perhaps use KV(1) to KV(end)???
                if (Xi1n.le.Xi1_min) Xi1n=Xi1_min
                if (Xi1n.ge.Xi1_max) Xi1n=Xi1_max
                if (Xi2n.le.Xi2_min) Xi2n=Xi2_min
                if (Xi2n.ge.Xi2_max) Xi2n=Xi2_max
                if (isnan(Xi1n)) Xi1n=Xi1
                if (isnan(Xi2n)) Xi2n=Xi2
c              4th: check for change of parameters
                pchange(1:3)=(Xi1n-Xi1)*SKL(1:3,2,1)
     +                      +(Xi2n-Xi2)*SKL(1:3,1,2)
                if (sqrt(dot(pchange,pchange,3)).le.(eps1*100)) then
                  converged = 1
                else if (count.ge.25) then
                  converged = 1
                  write (iow,*) 'NO LOCAL CONVERGENCE FOR DIRECTOR'
                  write (*,*)   'NO LOCAL CONVERGENCE FOR DIRECTOR'
                else
                  Xi1 = Xi1n
                  Xi2 = Xi2n
                end if
              end if
            end do
c           create director
            vec1(1:3) = SKL(1:3,2,1)
            vec2(1:3) = SKL(1:3,1,2)
c           compute lamina basis system (local cartesian), see Hughes, FEM, 2000, p. 385
            call norm (vec2(1:3),vec2(1:3),3)
            call norm (vec1(1:3),vec1(1:3),3)
            call vecp(vec1(1:3),vec2(1:3),vec3(1:3))
            if((vec3(1)+vec3(2)+vec3(3)).eq.0.0d0) then
                write(iow,*) 'zero director, panic!'
                vec1(1)=1.0d0;vec1(2)=0.0d0;vec1(3)=0.0d0
                vec2(1)=0.0d0;vec2(2)=1.0d0;vec2(3)=0.0d0
                vec3(1)=0.0d0;vec3(2)=0.0d0;vec3(3)=1.0d0
            end if
            vecIsNaN = 0
            do i=1,3
              if (isnan(vec3(i))) vecIsNan = 1
            end do   
            if (vecIsNaN.eq.1) then 
              write(iow,*) 'NaN-Director->set to e_i, ATTENTION'
              write(iow,*) node
              write(*  ,*) 'NaN-Director->set to e_i, ATTENTION'
              vec1(1)=1.0d0;vec1(2)=0.0d0;vec1(3)=0.0d0
              vec2(1)=0.0d0;vec2(2)=1.0d0;vec2(3)=0.0d0
              vec3(1)=0.0d0;vec3(2)=0.0d0;vec3(3)=1.0d0
            end if
            
            call norm (vec3(1:3),vec3(1:3),3)
            do i = 1,3
              gs(i,1) = vec1(i) + vec2(i)
            end do
            call vecp (vec3(1),gs(1,1),gs(1,2))
            do i = 1,3
              vec1(i) = gs(i,1) - gs(i,2)
              vec2(i) = gs(i,1) + gs(i,2)
            end do
            call norm (vec2(1:3),vec2(1:3),3)
            call norm (vec1(1:3),vec1(1:3),3)
c           end of orthonormalization
            do k = 1,3
              dir(k  ,node,ifdir) =  vec1(k)
              dir(k+3,node,ifdir) =  vec2(k)
              dir(k+6,node,ifdir) =  vec3(k)
            end do
            dir(10,node,ifdir) = dir(10,node,ifdir) + 1
          end do
        end do
        sum_in = sum_in + n+p+1
        sum_im = sum_im + m+q+1
        deallocate(KV1)
        deallocate(KV2)
        deallocate(CPw)
      end do
      if(prt) write(iow,2000) 
      if(prt) write(*  ,2000)
2000  format(/ '    Nodal director for NURBS surfaces determined by clos
     1est point projection'/)   
      return
      end      
c
c----------------------------------------------------------------------
      subroutine DersBasisFuns2(ni,Xi,p,KV,ders)
c----------------------------------------------------------------------
c
c      Purpose: Computes Basis Functions and their 1st and 2nd derivative
c               for B-Splines, Algorithm A2.3 from Piegl/Tiller
c
c      Inputs:
c         ni        - number of highest non-zero basis function
c         Xi        - NURBS parametric coordinate
c         p         - order of NURBS
c         KV        - NURBS knot vector
c
c      Outputs:
c         ders      - Basis functions and their first derivative
c
c----------------------------------------------------------------------
c     input/output variable declaration
      integer p,ni
      real*8 ders(3,p+1),KV(*), Xi
c     functions needed within subroutine
      real*8 left(p+1),right(p+1),saved,ndu(p+1,p+1),a(3,p+1),temp,d
      integer j,r,s1,s2,rk,pk,j1,j2
c
      left=0.d0; right=0.d0; ndu=0.d0; ndu(1,1)=1.d0; ders=0.d0; a=0.d0
c
c     store basis function and knot differences to ndu
      do j=1,p
        left(j) = Xi-KV(ni+1-j);
        right(j) = KV(ni+j)-Xi;
        saved = 0d0;
        do r=0,j-1
            ndu(j+1,r+1)=right(r+1)+left(j-r);
            temp = ndu(r+1,j)/ndu(j+1,r+1);
            ndu(r+1,j+1)=saved+right(r+1)*temp;
            saved = left(j-r)*temp;
        end do
        ndu(j+1,j+1)=saved;
      end do
c
c     Load basis funtions
      do j=0,p
        ders(1,j+1)=ndu(j+1,p+1);
      end do
c     compute the first derivative
      do r=0,p
        s1=0; s2=1
        a(1,1)=1d0;
c       loop for k to determine higher derivatives
        do k=1,2
          d=0d0; rk=r-k; pk=p-k
          if (r>=k) then
            a(s2+1,1)=a(s1+1,1)/ndu(pk+2,rk+1)
            d=a(s2+1,1)*ndu(rk+1,pk+1)
          end if
          if (rk >= -1) then
            j1=1
          else
            j1=-rk
          end if
          if (r-1 <= pk) then
            j2=k-1
          else
            j2=p-r
          end if
          do j=j1,j2
            a(s2+1,j+1)=(a(s1+1,j+1)-a(s1+1,j))/ndu(pk+2,rk+j+1)
            d = d+a(s2+1,j+1)*ndu(rk+j+1,pk+1)
          end do
          if (r<=pk) then
            a(s2+1,k+1)=-a(s1+1,k)/ndu(pk+2,r+1)
            d=d+a(s2+1,k+1)*ndu(r+1,pk+1)
          end if
          ders(k+1,r+1)=d;
          j=s1; s1=s2; s2=j
        end do
      end do
c     Multiply with correct factors
      r=p
      do k=1,2
        do j=0,p
          ders(k+1,j+1)=ders(k+1,j+1)*r;
        end do
        r = r*(p-k)
      end do
      return
      end
c
c----------------------------------------------------------------------
      subroutine SurfaceDerivsAlg2(p,KV1,q,KV2,xl,u,v,SKL,
     +                             uspan,vspan,ndm,d,Nu,Nv)
c----------------------------------------------------------------------
c
c      Purpose: Computes surface derivative
c               for B-Splines, Algorithm A3.6 from Piegl/Tiller
c
c      Inputs:
c         n,p,KV1,m,q,KV2,CPw,u,v
c
c      Outputs:
c         SKL
c
c----------------------------------------------------------------------
c     input/output variable declaration
      integer p,q,d,uspan,vspan
      real*8 KV1(*), KV2(*),xl(ndm,*),u,v,SKL(4,3,3)
c     variables needed within subroutine      
      real*8 Nu(3,p+1),Nv(3,q+1),temp(4,q+1)
      integer k,r,dd,l,s,du,dv
c     start: for p<2 set 2nd derivatives to zero
      d=min(2,d)
      du=min(d,p)
      do k=p+1,d
        do l=0,d-k
          SKL(1:4,k+1,l+1)=0.d0
        end do
      end do
      dv=min(d,q)
      do l=q+1,d
        do k=0,d-l
          SKL(1:4,k+1,l+1)=0.d0
        end do
      end do
c...  1D-Spline basis functions
      call DersBasisFuns2(uspan,u,p,KV1,Nu)

      call DersBasisFuns2(vspan,v,q,KV2,Nv)
c
      SKL=0.0d0

      do k=0,du
        do s=0,q     
          temp(1:4,s+1)=0.0d0
          do r=0,p
                 temp(1:3,s+1)=temp(1:3,s+1)+
     +            Nu(k+1,r+1)*xl(1:3,(p+1)*(q+1)-(s)*(p+1)-(r))
     *                        *xl(4,(p+1)*(q+1)-(s)*(p+1)-(r))
                 temp(4,s+1)=temp(4,s+1)+
     +            Nu(k+1,r+1)*xl(4,(p+1)*(q+1)-(s)*(p+1)-(r))
          end do
        end do
        dd = min(d-k,dv)
        do l = 0,dd
          SKL(1:4,k+1,l+1)=0.0d0
          do s = 0,q
            SKL(1:4,k+1,l+1)=SKL(1:4,k+1,l+1)+Nv(l+1,s+1)*temp(1:4,s+1)
          end do
        end do
      end do
      return
      end
c
c----------------------------------------------------------------------
      subroutine SurfaceDerivsAlg1(n,p,KV1,m,q,KV2,CPw,u,v,d,SKL)
c----------------------------------------------------------------------
c
c      Purpose: Computes surface derivative
c               for B-Splines, Algorithm A3.6 from Piegl/Tiller
c
c      Inputs:
c         n,p,KV1,m,q,KV2,CPw,u,v,d
c
c      Outputs:
c         SKL
c
c----------------------------------------------------------------------
c     input/output variable declaration
      integer n,p,m,q,d
      real*8 KV1(*), KV2(*),CPw(4,n,m),u,v,SKL(4,3,3)
c     variables needed within subroutine
      real*8 Nu(3,p+1),Nv(3,q+1),temp(4,q+1)
      integer k,r,dd,l,s,uspan,vspan
c     start: compute Spans and 1D-Spline basis functions
      call findspan(n,p,u,KV1,uspan)
      call DersBasisFuns2(uspan,u,p,KV1,Nu)
      call findspan(m,q,v,KV2,vspan)
      call DersBasisFuns2(vspan,v,q,KV2,Nv)
c
      SKL=0.0d0
      d=2
      do k=0,2
        do s=0,q
          temp(1:4,s+1)=0.0d0
          do r=0,p
            temp(1:4,s+1)=temp(1:4,s+1)+
     +            Nu(k+1,r+1)*CPw(1:4,uspan-p+r,vspan-q+s)
          end do
        end do
        dd = min(2-k,2)
        do l = 0,dd
          SKL(1:4,k+1,l+1)=0.0d0
          do s = 0,q
            SKL(1:4,k+1,l+1)=SKL(1:4,k+1,l+1)+Nv(l+1,s+1)*temp(1:4,s+1)
          end do
        end do
      end do
      return
      end
c----------------------------------------------------------------------
      subroutine RatSurfaceDerivs(Aders,wders,SKL)
c----------------------------------------------------------------------
c
c      Purpose: Computes surface derivative
c               for NURBS, Algorithm A4.4 from Piegl/Tiller
c
c      Inputs:
c         Aders, wders
c
c      Outputs:
c         SKL
c
c----------------------------------------------------------------------
c     input/output variable declaration
      real*8 Aders(3,3,3), wders(1,3,3), SKL(4,3,3)
c     variables needed within subroutine
      real*8 v(3),v2(3)
      integer i,j,k,l,Bin(3,3)
c     calculate binomial coefficients
      Bin(1,1)=1;Bin(2,1)=1;Bin(3,1)=1
      Bin(1,2)=-1;Bin(2,2)=1;Bin(3,2)=2
      Bin(1,3)=-1/4;Bin(2,3)=-1/2;Bin(3,3)=1
      SKL=0.0d0
      do k=0,2
        do l=0,2-k
          v = Aders(1:3,k+1,l+1)
          do j=1,l
            v = v - Bin(l+1,j+1)*wders(1,1,j+1)*SKL(1:3,k+1,l-j+1)
          end do
          do i =1,k
            v = v - Bin(k+1,i+1)*wders(1,i+1,1)*SKL(1:3,k-i+1,l+1)
            v2 = 0.0d0
            do j=1,l
              v2=v2+ Bin(l+1,j+1)*wders(1,i+1,j+1)*SKL(1:3,k-i+1,l-j+1)
            end do
            v = v - Bin(k+1,i+1)*v2
          end do
          SKL(1:3,k+1,l+1) = v/wders(1,1,1)
        end do
      end do
      return
      end
      subroutine findSpan(n,p,Xi,KV,span)
c----------------------------------------------------------------------
c
c     findSpan
c     KV: Knot Vector in given direction
c     n:  number of basis funtions in given direction
c     p:  order of basis functions in given direction
c     Xi: running variable
c     span: output
c
c     Algorithm FindSpan from Piegl/Tiller S. 68
c     changes because Piegl/Tiller start Knot Vector with Index 0
c----------------------------------------------------------------------
c     input/output variable declaration
      USE iofile
      real*8 KV(*), Xi
      integer n,p,span
c     variables needed within subroutine
      real*8 mid_real
      integer low,high,mid
c     start
      if (Xi>KV(n+p+1).or.Xi<KV(1)) then
        if (abs(Xi-KV(n+p+1))<1e-14) then
          span=KV(n+p+1);
        else if (abs(Xi-KV(1))<1e-14) then
          span=KV(1);
        else
          write(iow,*) 'error: value is outside Knot vector'
        end if
        return
      end if
      if (Xi.eq.KV(n+1)) then
        span = n
        return
      end if
      low = p+1;
      high = n+1;
      mid_real = (low+high)/2
      mid = floor(mid_real);
      do while (Xi < KV(mid) .or. Xi >= KV(mid+1))
        if (Xi< KV(mid)) then
          high = mid;
        else
          low = mid;
        end if
        mid_real = (low+high)/2.0d0
        mid = floor(mid_real);
      end do
      span = mid;
      return
      end


      subroutine getoldix(input,output,numel,nen)
      implicit double precision (a-h,o-z)
      integer input(nen,numel),output(nen,numel)
      output=input
      return
      end

      subroutine evaluateoldix(oldix,numel,nen,i,n,ii)
      integer oldix(nen,numel)
      integer nen,numel,i,n
      ii = abs(oldix(i,n))
      return
      end

      subroutine saveoldix(ix,nen1,nen,numel)
      USE isogeo
      integer nen1,numel,nen
      integer ix(nen1,numel)
csk      common m(1)
      call saveoldix1(ix,nen1,nen,numel,AInoix)
      return
      end

      subroutine saveoldix1(ix,nen1,nen,numel,oldix)
      implicit double precision (a-h,o-z)
      integer ix(nen1,numel),oldix(nen,numel)
      do i = 1,numel
        do j = 1,nen
          oldix(j,i)=ix(j,i)
        end do
      end do
      return
      end
      
      subroutine nnonzeroelements(KV1,KV2,n,m,p,q,nnzel,nel)
      integer n,m,p,q,nel,i,nel1,nel2
      real*8 KV1(n+p+1),KV2(m+q+1)
      
      nel1 = 0
      nel2 = 0
      do i = 1,n+p
        if (KV1(i).ne.KV1(i+1)) nel1 = nel1 + 1
      end do
      do i = 1,m+q
        if (KV2(i).ne.KV2(i+1)) nel2 = nel2 + 1
      end do
      nnzel = nel1 * nel2
      nel   = (n-p) * (m-q)
      return
      end

        
      subroutine setdrilldof(id,ip,ndf,numnp,numel,tol,
     1                       dir,prt,overwriteall,useglobalcart)
csk      subroutine setdrilldof(id,ix,ip,ndm,ndf,nen,nen1,numnp,numel,tol,
csk     1                       dir,prt,overwriteall,useglobalcart,ipr)
      USE iofile
      USE isogeo
      USE dirdat
      USE doalloc
c----------------------------------------------------------------------
c 
c     setdrilldof:
c     determines number of rotational degrees of freedom for every
c     control point.
c     2 for control points in smooth regions
c     3 for control points on kinks
c
c      Inputs:      
c         id(ndf,*)       - Equation numbers for each active dof
csk         ix(nen1,*)      - Element nodal connection list    
c         ip(*)           - Node number list for tied nodes     
csk         ndm             - Spatial dimension of mesh
c         ndf             - Number dof/node
csk         nen             - Number of nodes/element
csk         nen1            - Dimension of ix array
c         numnp           - Number of nodes in mesh
csk         numel           - Number of elements in mesh
c         tol             - tolerance angle in degrees    
c         dir(10,numnp,2) - Field of nodal directors   
c         prt             - write to screen about using ro56   
c         overwriteall    - for 0: free DOF 6, for 1: also fix DOF 6 
c     
c         useglobalcart   - for 1: use global cartesian system for 6 DOF nodes      
c      
c     determination by comparison of angle of directors of patch nodal
c     director and global nodal director        
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)       
      integer id(ndf,*),ip(numnp)
      integer list6dof(numnp)
      integer slavenode,masternode
      integer,pointer,dimension(:)  :: nodelist
      real*8  dir(10,numnp,2)
      real*8  a(3),b(3),dot
      integer counter,overwriteall,useglobalcart
      logical flag,prt,fldir2
      
csk      common m(1)
      PI=4.D0*DATAN(1.D0)
      tol_rad = abs(tol*PI/180.d0)
      list6dof = 0
      flag = .false.
      counter = 0

      if (useglobalcart.eq.1) then
c       allocate array for directors on patch level      
csk        call pseta(ImdirP,9*numnp,ipr,fldir2,'BASE-PW')
        call ralloc(AImdirP,9*numnp,'BASE-PW',fldir2)
        IldirP = 1
        call IGAcopydirectors(dir,AImdirP,numnp)
      end if
c     loop over all control points 
      do i = 1,numnp
        slavenode = i
        masternode = ip(i)
        if (slavenode.ne.masternode) then
          a(1:3)  = dir(7:9,slavenode ,1)
          b(1:3) = dir(7:9,masternode,1)
          c = dot(a,b,3)/dsqrt(dot(a,a,3)*dot(b,b,3))
          if (c.gt.1.d0.or.c.lt.-1.d0) then
            if (c.gt.1.d0.and.c.lt.1.0000000001d0) then
              c=1.0d0
            else if (c.lt.-1.d0.and.c.gt.-1.0000000001d0) then
              c=1.0d0
            else
              write (*,*) 'ERROR IN ro56'
            end if
          end if
          phi = abs(DACOS(c))
          if (phi.gt.tol_rad) then
            flag = .true.
            counter = counter + 1
c            list6dof(slavenode) = 1
            list6dof(masternode) = 1
c            write(*,*) phi*180/PI
          end if
        end if
      end do
c     now all master points on kinks are tagged for 6 DOFs
      allocate(nodelist(counter))
      counter = 0
      do i = 1,numnp
        if (list6dof(i).eq.1) then
          counter = counter + 1
          nodelist(counter) = i
c         set cartesian system for 6 DOF nodes
          if (useglobalcart.eq.1) then
            dir(1,i,1) = 1.d0
            dir(2,i,1) = 0.d0
            dir(3,i,1) = 0.d0
            dir(4,i,1) = 0.d0
            dir(5,i,1) = 1.d0
            dir(6,i,1) = 0.d0
            dir(7,i,1) = 0.d0
            dir(8,i,1) = 0.d0
            dir(9,i,1) = 1.d0
          end if
        else if (list6dof(i).eq.0.and.overwriteall.eq.1) then
          id(6,i) = 1
        end if
      end do
c     now only output
      write(iow,2000)
      if(prt) write(  *,2000)
      if(flag) then
        if(prt) write(  *,2001)
        write(iow,2001)
        if(prt) write(  *,2002) nodelist
        write(iow,2002) nodelist
        if (overwriteall.eq.1) then
          if(prt) write(  *,2005)
          write(iow,2005)  
        else if (overwriteall.eq.0) then
          if(prt) write(  *,2006)
          write(iow,2006)  
        end if
      else
        if (overwriteall.eq.1) then
          if(prt) write(*  ,2003)
          write(iow,2003)
        else if (overwriteall.eq.0) then
          if(prt) write(*  ,2007)
          write(iow,2007)
        end if
      end if
      write(iow,2004) tol
      if(prt) write(  *,2004) tol
      deallocate(nodelist)
      return
2000  format(/6x,' r o 5 6  - -  n o d a l   r o t a t i o n a l DOFs'/
     1    6x,'    Determination of number of rotational DOFs for every 
     2control point.'/)
2001  format(6x,'    The following nodes have been assigned 6 DOFs:'/)
2002  format(6x,10i5)
2003  format(10x,'No nodes have been assigned 6 DOFs, all 6th DOFs are
     1fixed!'/1x)
2004  format(/,6x,'ro56: Angle tolerance in degree=',1p,1e12.5/1x)      
2005  format(10x,'For all other nodes DOF 6 is fixed, regardless of valu
     1es specified before in boun etc.'/1x)   
2006  format(10x,'For all other nodes DOF 6 remains unchanged'/1x)
2007  format(10x,'No nodes have been assigned 6 DOFs, all 6th DOFs remai
     1n unchanged!'/1x)      
      end    
      
      subroutine IGAcopydirectors(dirG,dirP,numnp)
      implicit double precision (a-h,o-z)
      real*8 dirG(10,numnp,2),dirP(9,numnp)
      do i = 1,numnp
        do j = 1,9
          dirP(j,i)=dirG(j,i,1)
        end do
      end do
      return
      end
c----------------------------------------------------------------------     
      subroutine IGAreadpatchdirectors(dirP,numnp,diriP,node)
      implicit double precision (a-h,o-z)
      real*8 dirP(9,numnp),diriP(3,3)
      do k = 1,3
        diriP(k,1) = dirP(k  ,node)
        diriP(k,2) = dirP(k+3,node)
        diriP(k,3) = dirP(k+6,node)
      end do
      return
      end
c----------------------------------------------------------------------
!      subroutine IGAcompStorNormal(SKL,blocal,igp,ielem)
!c----------------------------------------------------------------------
!c     computes the normal vector and its derivatives from the derivatives
!c     of X. Transformation to local orthonormal basis systems and saved to
!c     global array   IGAbasissystems   
!c----------------------------------------------------------------------      
!      implicit double precision (a-h,o-z)
!c     input variables      
!      real*8 blocal(3,3), SKL(4,3,3)
!      integer igp,ielem
!c     internal variables      
!      real*8 N(3),N_1(3),N_2(3),Nhat(3),Nhat_1(3),Nhat_2(3)
!      real*8 normN_1,normN_2
!      real*8 length,temp1(3),temp2(3)
!      real*8 detj,xs(2,2),sx(2,2)
!      real*8 lXi1,eXi1(3),lXi2,eXi2(3),lXi1der(2),lXi2der(2)
!      real*8 eXi1d1(3),eXi1d2(3),eXi2d1(3),eXi2d2(3),ealpha(3)
!      real*8 ealphad1(3),ealphad2(3),ebetad1(3),ebetad2(3)
!      real*8 el1d1(3),el1d2(3),el2d1(3),el2d2(3),ebeta(3)
!      real*8 hel1d1(3),hel1d2(3),hel2d1(3),hel2d2(3)
!      real*8 el1(3),el2(3),hel1(3),hel2(3),lhel1,lhel2
!      real*8 lhel1d1,lhel1d2,lhel2d1,lhel2d2
!      integer i,j,k,l,m
!      include 'isogeo.h'
!
!c     normalized normal is stored in blocal(7:9)
!      Nhat(1:3) = SKL(1:3,1,1)
!      length = dsqrt(dot(Nhat,Nhat,3))
!      N(1:3) = blocal(1:3,3)
!
!c     derivative of Nhat_1
!      call vecp(SKL(1:3,3,1),SKL(1:3,1,2),temp1)
!      call vecp(SKL(1:3,2,1),SKL(1:3,2,2),temp2)
!      Nhat_1(1:3) = temp1(1:3)+temp2(1:3)
!c     derivative of Nhat_2
!      call vecp(SKL(1:3,2,2),SKL(1:3,1,2),temp1)
!      call vecp(SKL(1:3,2,1),SKL(1:3,1,3),temp2)
!      Nhat_2(1:3) = temp1(1:3)+temp2(1:3)
!c     derivative of norm       
!      normN_1 = dot(Nhat_1,Nhat,3)/length
!      normN_2 = dot(Nhat_2,Nhat,3)/length
!
!c     parametric derivatives      
!      do i = 1,3
!        N_1(i) = (length*Nhat_1(i)-normN_1*Nhat(i))/(length**2)
!      end do
!      do i = 1,3
!        N_2(i) = (length*Nhat_2(i)-normN_2*Nhat(i))/(length**2)
!      end do
!
!c.... construct jacobian and its inverse
!      do j = 1,2
!        xs(1,j) = dot(SKL(1:3,2,1),blocal(1,j),3)
!        xs(2,j) = dot(SKL(1:3,1,2),blocal(1,j),3)
!      end do
!      detj = xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
!      sx(1,1) = xs(2,2)/detj
!      sx(2,2) = xs(1,1)/detj
!      sx(1,2) =-xs(1,2)/detj
!      sx(2,1) =-xs(2,1)/detj 
!      
!c     storing values of director
!      do i = 1,3
!        IGAbasissystems(i,1,3,igp,ielem) = sx(1,1)*N_1(i)+sx(1,2)*N_2(i)
!        IGAbasissystems(i,2,3,igp,ielem) = sx(2,1)*N_1(i)+sx(2,2)*N_2(i)
!      end do
!      IGAbasissystems(1:3,3,3,igp,ielem) = N(1:3)    
!      
!c     basis vectors A_\alpha from lamina_base of derivatives
!      IGAbasissystems(1:3,3,1,igp,ielem) = blocal(1:3,1)
!      IGAbasissystems(1:3,3,2,igp,ielem) = blocal(1:3,2)
!      
!c     alternatively use Jacobian (small numerical error, thus not used)     
!c      do i = 1,3
!c        IGAbasissystems(i,3,1,igp,ielem) = 
!c     =   sx(1,1)*SKL(i,2,1)+sx(1,2)*SKL(i,1,2)
!c        IGAbasissystems(i,3,2,igp,ielem) = 
!c     =   sx(2,1)*SKL(i,2,1)+sx(2,2)*SKL(i,1,2)
!c      end do
!
!c     calculate derivatives of basis vectors: A_\alpha,\beta
!c     by derivating the lamina basis system
!      lXi1 = dsqrt(dot(SKL(1:3,2,1),SKL(1:3,2,1),3))
!      eXi1 = SKL(1:3,2,1)/lXi1
!      lXi2 = dsqrt(dot(SKL(1:3,1,2),SKL(1:3,1,2),3))
!      eXi2 = SKL(1:3,1,2)/lXi2
!      lXi1der(1) = (dot(SKL(1:3,3,1),SKL(1:3,2,1),3))/lXi1
!      lXi1der(2) = (dot(SKL(1:3,2,2),SKL(1:3,2,1),3))/lXi1
!      lXi2der(1) = (dot(SKL(1:3,2,2),SKL(1:3,1,2),3))/lXi2
!      lXi2der(2) = (dot(SKL(1:3,1,3),SKL(1:3,1,2),3))/lXi2
!      
!      eXi1d1 = (lXi1*SKL(1:3,3,1)-lXi1der(1)*SKL(1:3,2,1))/(lXi1**2)
!      eXi1d2 = (lXi1*SKL(1:3,2,2)-lXi1der(2)*SKL(1:3,2,1))/(lXi1**2)
!      eXi2d1 = (lXi2*SKL(1:3,2,2)-lXi2der(1)*SKL(1:3,1,2))/(lXi2**2)
!      eXi2d2 = (lXi2*SKL(1:3,1,3)-lXi2der(2)*SKL(1:3,1,2))/(lXi2**2)
!      
!      ealpha = eXi1 + eXi2
!      call vecp(N,ealpha,ebeta)
!      hel1=(ealpha-ebeta)
!      hel2=(ealpha+ebeta)
!      lhel1 = dsqrt(dot(hel1,hel1,3))
!      lhel2 = dsqrt(dot(hel2,hel2,3))
!      el1 = hel1 / lhel1
!      el2 = hel2 / lhel2
!      
!      ealphad1 = eXi1d1 + eXi2d1
!      ealphad2 = eXi1d2 + eXi2d2
!      
!      call vecp(N_1,ealpha,temp1)
!      call vecp(N,ealphad1,temp2)
!      ebetad1 = temp1 + temp2
!      call vecp(N_2,ealpha,temp1)
!      call vecp(N,ealphad2,temp2)
!      ebetad2 = temp1 + temp2
!      
!      hel1d1 = (ealphad1-ebetad1)
!      hel1d2 = (ealphad2-ebetad2)
!      hel2d1 = (ealphad1+ebetad1)
!      hel2d2 = (ealphad2+ebetad2)
!      
!      lhel1d1 = dot(hel1d1,hel1,3) / lhel1
!      lhel1d2 = dot(hel1d2,hel1,3) / lhel1
!      lhel2d1 = dot(hel2d1,hel2,3) / lhel2
!      lhel2d2 = dot(hel2d2,hel2,3) / lhel2
!      
!      el1d1 = (hel1d1*lhel1 - lhel1d1*hel1)/(lhel1**2)
!      el1d2 = (hel1d2*lhel1 - lhel1d2*hel1)/(lhel1**2)
!      el2d1 = (hel2d1*lhel2 - lhel2d1*hel2)/(lhel2**2)
!      el2d2 = (hel2d2*lhel2 - lhel2d2*hel2)/(lhel2**2)
!      
!c     transformation to local cartesian system
!      do i = 1,3
!      IGAbasissystems(i,1,1,igp,ielem)=sx(1,1)*el1d1(i)+sx(1,2)*el1d2(i)
!      IGAbasissystems(i,2,1,igp,ielem)=sx(2,1)*el1d1(i)+sx(2,2)*el1d2(i)
!      end do
!      do i = 1,3
!      IGAbasissystems(i,1,2,igp,ielem)=sx(1,1)*el2d1(i)+sx(1,2)*el2d2(i)
!      IGAbasissystems(i,2,2,igp,ielem)=sx(2,1)*el2d1(i)+sx(2,2)*el2d2(i)
!      end do
!      
!      
!     
!      
!!c     calculate derivatives of basis vectors: A_\alpha,\beta
!!      do i= 1,2
!!        do j = 1,2
!!          do k = 1,3
!!            IGAbasissystems(k,j,i,igp,ielem) =
!!     +      sx(i,1)*sx(j,1)*SKL(k,3,1)+
!!     +      sx(i,1)*sx(j,2)*SKL(k,2,2)+
!!     +      sx(i,2)*sx(j,1)*SKL(k,2,2)+            
!!     +      sx(i,2)*sx(j,2)*SKL(k,1,3)!+
!!c     +      sx(i,1)*sx(j,1)*SKL(k,2,1)+
!!c     +      sx(i,2)*sx(j,2)*SKL(k,1,2)            
!!c          do l = 1,2
!!c          do m = 1,2
!!c            IGAbasissystems(k,j,i,igp,ielem) =
!!c     =      IGAbasissystems(k,j,i,igp,ielem) + 
!!c            sx(i,l)*sx(j,m)*SKL(k,l,m
!!c          end do
!!c         end do
!!          end do
!!        end do
!!      end do
!      return
!      end
      
    
c----------------------------------------------------------------------
      subroutine IGAreadNormalm(IGAbasissystems,numel,ngp2,
     +                          output,igp,ielem)
      integer numel,ngp2,igp,ielem
      real*8 IGAbasissystems(3,3,3,ngp2,numel),output(3,3,3)
      output = IGAbasissystems(1:3,1:3,1:3,igp,ielem)
      return
      end
            
c----------------------------------------------------------------------
      subroutine IGAcompStorNormalm1(SKL,blocal,igp,ielem,AInbasissys,
     +                              numel,ngp2)
      implicit double precision (a-h,o-z)
      real*8 SKL(4,3,3),blocal(3,3),AInbasissys
      integer igp,ielem,numel,ngp2
csk      common m(1)
      call IGAcompStorNormalm(SKL,blocal,igp,ielem,
     +                                AInbasissys,numel,ngp2)
      return
      end
c----------------------------------------------------------------------
      subroutine IGAcompStorNormalm(SKL,blocal,igp,ielem,IGAbasissystems
     +                              ,numel,ngp2)
      USE isogeo
c----------------------------------------------------------------------
c     computes the normal vector and its derivatives from the derivatives
c     of X. Transformation to local orthonormal basis systems and saved to
c     global array   IGAbasissystems   
c----------------------------------------------------------------------      
      implicit double precision (a-h,o-z)
c     input variables      
      real*8 blocal(3,3), SKL(4,3,3)
      real*8 IGAbasissystems(3,3,3,ngp2,numel)
      integer igp,ielem,numel,ngp2
c     internal variables      
      real*8 N(3),N_1(3),N_2(3),Nhat(3),Nhat_1(3),Nhat_2(3)
      real*8 normN_1,normN_2
      real*8 length,temp1(3),temp2(3)
      real*8 detj,xs(2,2),sx(2,2)
      real*8 lXi1,eXi1(3),lXi2,eXi2(3),lXi1der(2),lXi2der(2)
      real*8 eXi1d1(3),eXi1d2(3),eXi2d1(3),eXi2d2(3),ealpha(3)
      real*8 ealphad1(3),ealphad2(3),ebetad1(3),ebetad2(3)
      real*8 el1d1(3),el1d2(3),el2d1(3),el2d2(3),ebeta(3)
      real*8 hel1d1(3),hel1d2(3),hel2d1(3),hel2d2(3)
csk      real*8 el1(3),el2(3),hel1(3),hel2(3),lhel1,lhel2
      real*8 hel1(3),hel2(3),lhel1,lhel2
      real*8 lhel1d1,lhel1d2,lhel2d1,lhel2d2
      integer i,j !csk ,k,l,m

c     normalized normal is stored in blocal(7:9)
      Nhat(1:3) = SKL(1:3,1,1)
      length = dsqrt(dot(Nhat,Nhat,3))
      N(1:3) = blocal(1:3,3)

c     derivative of Nhat_1
      call vecp(SKL(1:3,3,1),SKL(1:3,1,2),temp1)
      call vecp(SKL(1:3,2,1),SKL(1:3,2,2),temp2)
      Nhat_1(1:3) = temp1(1:3)+temp2(1:3)
c     derivative of Nhat_2
      call vecp(SKL(1:3,2,2),SKL(1:3,1,2),temp1)
      call vecp(SKL(1:3,2,1),SKL(1:3,1,3),temp2)
      Nhat_2(1:3) = temp1(1:3)+temp2(1:3)
c     derivative of norm       
      normN_1 = dot(Nhat_1,Nhat,3)/length
      normN_2 = dot(Nhat_2,Nhat,3)/length

c     parametric derivatives      
      do i = 1,3
        N_1(i) = (length*Nhat_1(i)-normN_1*Nhat(i))/(length**2)
      end do
      do i = 1,3
        N_2(i) = (length*Nhat_2(i)-normN_2*Nhat(i))/(length**2)
      end do

c.... construct jacobian and its inverse
      do j = 1,2
        xs(1,j) = dot(SKL(1:3,2,1),blocal(1,j),3)
        xs(2,j) = dot(SKL(1:3,1,2),blocal(1,j),3)
      end do
      detj = xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
      sx(1,1) = xs(2,2)/detj
      sx(2,2) = xs(1,1)/detj
      sx(1,2) =-xs(1,2)/detj
      sx(2,1) =-xs(2,1)/detj 
      
c     storing values of director
      do i = 1,3
        IGAbasissystems(i,1,3,igp,ielem) = sx(1,1)*N_1(i)+sx(1,2)*N_2(i)
        IGAbasissystems(i,2,3,igp,ielem) = sx(2,1)*N_1(i)+sx(2,2)*N_2(i)
      end do
      IGAbasissystems(1:3,3,3,igp,ielem) = N(1:3)    
      
c     basis vectors A_\alpha from lamina_base of derivatives
      IGAbasissystems(1:3,3,1,igp,ielem) = blocal(1:3,1)
      IGAbasissystems(1:3,3,2,igp,ielem) = blocal(1:3,2)
      
c     alternatively use Jacobian (small numerical error, thus not used)     
c      do i = 1,3
c        IGAbasissystems(i,3,1,igp,ielem) = 
c     =   sx(1,1)*SKL(i,2,1)+sx(1,2)*SKL(i,1,2)
c        IGAbasissystems(i,3,2,igp,ielem) = 
c     =   sx(2,1)*SKL(i,2,1)+sx(2,2)*SKL(i,1,2)
c      end do

c     calculate derivatives of basis vectors: A_\alpha,\beta
c     by derivating the lamina basis system
      lXi1 = dsqrt(dot(SKL(1:3,2,1),SKL(1:3,2,1),3))
      eXi1 = SKL(1:3,2,1)/lXi1
      lXi2 = dsqrt(dot(SKL(1:3,1,2),SKL(1:3,1,2),3))
      eXi2 = SKL(1:3,1,2)/lXi2
      lXi1der(1) = (dot(SKL(1:3,3,1),SKL(1:3,2,1),3))/lXi1
      lXi1der(2) = (dot(SKL(1:3,2,2),SKL(1:3,2,1),3))/lXi1
      lXi2der(1) = (dot(SKL(1:3,2,2),SKL(1:3,1,2),3))/lXi2
      lXi2der(2) = (dot(SKL(1:3,1,3),SKL(1:3,1,2),3))/lXi2
      
      eXi1d1 = (lXi1*SKL(1:3,3,1)-lXi1der(1)*SKL(1:3,2,1))/(lXi1**2)
      eXi1d2 = (lXi1*SKL(1:3,2,2)-lXi1der(2)*SKL(1:3,2,1))/(lXi1**2)
      eXi2d1 = (lXi2*SKL(1:3,2,2)-lXi2der(1)*SKL(1:3,1,2))/(lXi2**2)
      eXi2d2 = (lXi2*SKL(1:3,1,3)-lXi2der(2)*SKL(1:3,1,2))/(lXi2**2)
      
      ealpha = eXi1 + eXi2
      call vecp(N,ealpha,ebeta)
      hel1=(ealpha-ebeta)
      hel2=(ealpha+ebeta)
      lhel1 = dsqrt(dot(hel1,hel1,3))
      lhel2 = dsqrt(dot(hel2,hel2,3))
csk      el1 = hel1 / lhel1
csk      el2 = hel2 / lhel2
      
      ealphad1 = eXi1d1 + eXi2d1
      ealphad2 = eXi1d2 + eXi2d2
      
      call vecp(N_1,ealpha,temp1)
      call vecp(N,ealphad1,temp2)
      ebetad1 = temp1 + temp2
      call vecp(N_2,ealpha,temp1)
      call vecp(N,ealphad2,temp2)
      ebetad2 = temp1 + temp2
      
      hel1d1 = (ealphad1-ebetad1)
      hel1d2 = (ealphad2-ebetad2)
      hel2d1 = (ealphad1+ebetad1)
      hel2d2 = (ealphad2+ebetad2)
      
      lhel1d1 = dot(hel1d1,hel1,3) / lhel1
      lhel1d2 = dot(hel1d2,hel1,3) / lhel1
      lhel2d1 = dot(hel2d1,hel2,3) / lhel2
      lhel2d2 = dot(hel2d2,hel2,3) / lhel2
      
      el1d1 = (hel1d1*lhel1 - lhel1d1*hel1)/(lhel1**2)
      el1d2 = (hel1d2*lhel1 - lhel1d2*hel1)/(lhel1**2)
      el2d1 = (hel2d1*lhel2 - lhel2d1*hel2)/(lhel2**2)
      el2d2 = (hel2d2*lhel2 - lhel2d2*hel2)/(lhel2**2)
      
c     transformation to local cartesian system
      do i = 1,3
      IGAbasissystems(i,1,1,igp,ielem)=sx(1,1)*el1d1(i)+sx(1,2)*el1d2(i)
      IGAbasissystems(i,2,1,igp,ielem)=sx(2,1)*el1d1(i)+sx(2,2)*el1d2(i)
      end do
      do i = 1,3
      IGAbasissystems(i,1,2,igp,ielem)=sx(1,1)*el2d1(i)+sx(1,2)*el2d2(i)
      IGAbasissystems(i,2,2,igp,ielem)=sx(2,1)*el2d1(i)+sx(2,2)*el2d2(i)
      end do
      
      
     
      
!c     calculate derivatives of basis vectors: A_\alpha,\beta
!      do i= 1,2
!        do j = 1,2
!          do k = 1,3
!            IGAbasissystems(k,j,i,igp,ielem) =
!     +      sx(i,1)*sx(j,1)*SKL(k,3,1)+
!     +      sx(i,1)*sx(j,2)*SKL(k,2,2)+
!     +      sx(i,2)*sx(j,1)*SKL(k,2,2)+            
!     +      sx(i,2)*sx(j,2)*SKL(k,1,3)!+
!c     +      sx(i,1)*sx(j,1)*SKL(k,2,1)+
!c     +      sx(i,2)*sx(j,2)*SKL(k,1,2)            
!c          do l = 1,2
!c          do m = 1,2
!c            IGAbasissystems(k,j,i,igp,ielem) =
!c     =      IGAbasissystems(k,j,i,igp,ielem) + 
!c            sx(i,l)*sx(j,m)*SKL(k,l,m
!c          end do
!c         end do
!          end do
!        end do
!      end do
      return
      end