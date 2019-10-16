      subroutine elmt15(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c-----------------------------------------------------------------------
c      (2,3,4,5) node geometrically nonlinear 3D-beam element
c      Timoshenko beam theory, Reissner strains, finite rotations
c      Interface to mate3d08 (FE2)
c
c      (c)  F. Gruttmann 
c-----------------------------------------------------------------------
c.... Allocation of d-array with material and geometrical parameters
c-----------------------------------------------------------------------
c     1. card 
c     d(1) = E
c     d(2) = G
c     d(3) = A
c     d(4) = Iy
c     d(5) = Iz
c     d(6) = Iyz
c     d(7) = It
c     d(8) = ys
c     d(9) = zs
c     d(10)= ym
c     d(11)= zm
c     d(12)= alpha
c          local basis:  e_1  (= tangent vector) ,  e_Z=[0,0,1]
c                        e_2 = -e_1 x e_Z   (  if e_1 .ne. +-e_Z ) 
c                        e_2 = -e_Y         (  if e_1 .eq. +-e_Z )
c          rotate local coor. system with angle alpha (grad),
c          alpha is defined from 2-axis in direction to 3-axis
c
c-----------------------------------------------------------------------
c     2. card
c     d(13)= ilin (geom. nonl.: 0 = lin, 1 = mod, 2 = finite, 3 fin>360)
c     d(14)= iPlo (output: 0 = reference, 1 = current configuration)
c     d(15)= ish  (FE shear correction factor, 0=off,1=on)     
c     d(16)= rho  (density)
c     d(17)= imat (0= linear elasticity, 1 = FE2) 
c 
c          output of stress resultants at nodes               (isw  4)
c          compute stress resultants at nx points             (isw 13)
c          iPlo=0: related to reference configuration
c          iplo=1: related to current   configuration
c
c-----------------------------------------------------------------------
c     3. card: material data for type IMAT = d(17)= 8  (FE2)
c-----------------------------------------------------------------------
c
c.....  element loads p(i)(i=1..3) at point yp, zp  via qload  (isw 22)
c
c-----------------------------------------------------------------------
      USE bdata
      USE cdat1
      USE cdata
      USE dirdat
      USE eldata
      USE fe2mat    ! for matfe2=ma
      USE fornam
      USE iofile
      USE pdata6
      USE qload
      implicit double precision (a-h,o-z)
      dimension xl(ndm,*),tl(*),d(*),ul(ndf,*),s(nst,*),p(*),ix(*)
     +         ,shp(2,5),pg(5),wg(5),dmat(6,6),eps(6),sig(6),b(6,6,5)
     +         ,btd(6,6),ai(3,3),as(3,3),xs(3),rin(3,3,5),rwn(3,3,5)
     +         ,hwn(3,3,5),dwn(3,5),swn(3,3,5),pin(3,3,5),ddxyz(3)
     +         ,twh(3,3,3,5),ttt(3,3,6,15),h1(*),h2(*),h3(*),xgp(3)
     +         ,plout(10)
c-----------------------------------------------------------------------
      ielno=15
c-----------------------------------------------------------------------
c.... go to correct array processor
      go to(1,2,3,3,3,3,2,2,2,2,3,2,3,2,3,2,2,18,2,2,2,3), isw
 
c.... input material properties
1     if(ior.lt.0) write(*,1000)
1000  format('Input: E G A Iy Iz Iyz It ys zs ym zm alpha'/3x,'>',$)
      call dinput(d,12)
      if(ior.lt.0) write(*,1001)
1001  format('Input:ilin, iPlo, ish, rho '/3x,'>',$)
      call dinput(d(13),5)
      if(ior.lt.0) write(*  ,1002) (d(i),i=1,17)
                   write(iow,1002) (d(i),i=1,17)
1002  format(5x,//'  3D Finite Deformation Timoshenko Beam'//
     +  6x,'Elastic Modulus                             E     =',e13.5/
     +  6x,'Shear Modulus                               G     =',e13.5/
     +  6x,'Area                                        A     =',e13.5/
     +  6x,'Moment of inertia                           I_2   =',e13.5/
     +  6x,'Moment of inertia                           I_3   =',e13.5/
     +  6x,'Moment of inertia                           I_23  =',e13.5/
     +  6x,'Saint-Venant torsion stiffness              I_T   =',e13.5/
     +  6x,'y-coor. of the centroid                     y_s   =',e13.5/
     +  6x,'z-coor. of the centroid                     z_s   =',e13.5/
     +  6x,'y-coor. of the center of shear              y_m   =',e13.5/
     +  6x,'z-coor. of the center of shear              z_m   =',e13.5/
     +  6x,'local rotation angle y->z in degrees        alpha =',e13.5/
     +  6x,'0=lin.  1=mod.   2=finite   3=finite>360°   ilin  =',e13.5/
     +  6x,'Output rel. to ref./cur. config.(=0/1)      iPlo  =',e13.5/
     +  6x,'Shear correction factor (0=off,1=on)        ish   =',e13.5/
     +  6x,'density                                     rho   =',e13.5/
     +  6x,'material type:   0: elasticity,  1: FE2     imat  =',e13.5/)
      if (ndf.ne.6) stop '*** element 15 only with ndf = 6! ***'
c
      imat  = d(17)
      if(imat.eq.1)imat = 8
      mdx = 0
      nh = 0
      matfe2=ma  
      if(imat.ne.0)call matelib3d(h1,h2,nh,d(21),mdx,eps,sig,dmat,6,4,
     +                 plout,xgp,tgp,dv,1.d0,1.d0,1.d0,n,l,1,1,imat,isw)

      d(20) =  nh                      ! number of stored param. per GP
      nmax = 20 + mdx
      if(nmax.gt.ndd) then
         write(*,1006) ndd,nmax
1006     format(1x,'darray is set to  ',i4,' values',/,
     1          1x,'darray needs      ',i4,' values')
         stop
      end if


c.... define name for stresses in plot-output
      ilin=d(13)
      if(ilin.eq.0)d(14)=1.d0
      iPlo=d(14)
      if (iPlo.eq.0) then
        forsus( 1) = 'Force H_1  '
        forsus( 2) = 'Force V_2  '
        forsus( 3) = 'Force V_3  '
        forsus( 4) = 'Moment M_1 '
        forsus( 5) = 'Moment M_2 '
        forsus( 6) = 'Moment M_3 '
      else
        forsus( 1) = 'Force N_1  '
        forsus( 2) = 'Force Q_2  '
        forsus( 3) = 'Force Q_3  '
        forsus( 4) = 'Moment M_1 '
        forsus( 5) = 'Moment M_2 '
        forsus( 6) = 'Moment M_3 '
      endif
      forsus( 7) = ' '
      forsus( 8) = ' '
      forsus( 9) = ' '
      forsus(10) = ' '
      forsus(11) = ' '

c.... define node numbering for plot mesh routine
c.... for max. 5 nodes per element

      inord(ielno)  = 2*nen
      do i=1,nen
        ipord(i,   ielno) = i
        ipord(2*nen+1-i,ielno) = i
      end do
c.... cheque, if "base,1" is defined in ifile:
      if (ilin.eq.3.and.mdir.eq.1) then
        write(*,*) 'The following command is required in the ifile:'
        write(*,*) 
        write(*,*) 'base'
        write(*,*) '1'
        stop
      endif
2     return
c.... residuum vector and stiffness matrix
3     continue
      imat = d(17)
      if(imat.eq.1)imat = 8
      ngb=nel-1
      call gaus1d(ngb,pg,wg)
c.... update nodal values
      call updn15(d,xl,ul,shp,rin,rwn,hwn,dwn,swn,pin,twh,ttt,ndm,ndf)
      if (isw.eq.22) goto 22
      do i = 1,3
       ddxyz(i) = xl(i,1)-xl(i,nel)
      enddo
      dlx = dsqrt(ddxyz(1)*ddxyz(1)+ddxyz(2)*ddxyz(2)+ddxyz(3)*ddxyz(3)) 
      if(isw.eq.5) go to 5
      call dmat15 (d,dmat,dlx,scfy,scfz)
      xgp = 0
      tgp = 0
      detf = 1
      matfe2=ma  
c.... loop over gauss points
      do 300 l = 1,ngb
        xsi = pg(l)
        call shape1D(shp,xsi,xl,dl,ndm,nel)
        dxl = wg(l)*dl
c....   beam strains and stress resultants 
        call strain15(d,xl,ul,shp,rin,rwn,pin,ai,as,xs,eps,ndm,ndf)
        if(imat.eq.0)then 
          call mvmul(dmat,eps,6,6,sig )
        else
          call matelib3d(h1,h2,nh,d(21),mdx,eps,sig,dmat,6,4,plout,
     +                   xgp,tgp,dlx,detf,scfy,scfz,n,l,1,1,imat,isw)
        endif 
        if(isw.eq.15)go to 300 
c....   loop over rows
        ir0 = 0
        do 310 inode = 1,nel
          call bmat15(d,inode,ir0,p,s,shp,ai,as,xs,rwn,hwn,dwn,swn,
     +                twh,sig,b,btd,dmat,dxl,nst,isw)
          if (isw.eq. 4 .or. isw.eq. 6 .or. isw.eq.13) goto 310
c....     loop over columns (symmetry noted, lower part)
c....     tangent stiffness
          jc0 = 0
          do jnode = 1,inode
            call stif15(inode,jnode,ir0,jc0,d,sig,shp,b,btd,s,
     +                  twh,ttt,dxl,nst)
            jc0 = jc0 + ndf
          enddo
310     ir0 = ir0 + ndf
300   continue
      if (isw.eq.4 .or. isw.eq.13) goto  4
      if (isw.eq.6.or.isw.eq.15) return
c.... add upper part of stiffness matrix 
      call msym(s,nst,nst,1)
c        if(n.eq.numel)call mprint(s,12,12,nst,'s15 ')
      return

c.... output forces  R = P - BtS  at element nodes (isw=4)
c.... and plot stress resultants diagramm on frame (isw=13)
4     call pout15(d,aqloa,xl,ul,p,rin,rwn,hwn,pg,wg,ngb,ndf,ndm,isw,
     +            numel)
      return
c.... lumped mass matrix 
5     rhoA = d(16)*d(3)
      fact = 0.5d0*d(1)*(d(4)+d(5))/d(2)/d(3)     ! 0.5*E(Iy + Iz)/GA
      pm   = rhoA*dlx*0.5d0
      pr   = pm *fact
      if(nel.gt.2)stop 'mass matrix only for 2-node element'
      do i = 1,3 
       p(i) = pm
       p(3+i) = pr
       p(ndf+i) = pm
       p(ndf+3+i) = pr
      enddo
      return
c.... initial nodal cartesian basis
18    if (knode.eq.(nen*numel)) then
        alpha = d(12) 
        do k = 1,nel
          xsi = 2.d0/(nel-1) * (k - (nel+1)/2.d0)
          call triad15 (shp,xl(1,1),rin(1,1,k),xsi,alpha,ndm)
          call pdirec2 (basea,rin(1,1,k),k,n,nen,1,1)
          call pdirec2 (basea,rin(1,1,k),k,n,nen,1,2)
        enddo
      endif
      return

c.... element loads at loading point {p2,p3} (for qloa)
22    call qloa15(d,aqloa,xl,rwn,hwn,dwn,swn,pg,wg,ngb,p,s,numel,
     +            nst,ndm,ndf)
      return
      end
c
      subroutine strain15(d,xl,ul,shp,rin,rwn,pin,ai,as,xs,eps,ndm,ndf)
c-----------------------------------------------------------------------
c     calculate strains at reference axis
c
c     beam strains     :  eps| gamma_2| gamma_3| theta| kappa_2| kappa_3
c     stress resultants:  F_1| F_2    | F_3    | M_1  | M_2    | M_3
c
c-----------------------------------------------------------------------
c     Ai_k    = [rin(1,i,k),rin(2,i,k),rin(3,i,k)]^T         !..Kn-Werte
c     ai_k    = [rwn(1,i,k),rwn(2,i,k),rwn(3,i,k)]^T         !..Kn-Werte

c     x'      = [xs(1),xs(2),xs(3)]^T                        !..GP-Werte
c     X'      = [XXs(1),XXs(2),XXs(3)]^T                     !..GP-Werte

c     u'      = [us(1),us(2),us(3)]^T                        !..GP-Werte
c     phi     = [phi(1),phi(2),phi(3)]^T                     !..GP-Werte
c     phi'    = [phis(1),phis(2),phis(3)]^T                  !..GP-Werte

c     ai      = [ai(1,i),ai(2,i),ai(3,i)]^T                  !..GP-Werte
c     ai'     = [as(1,i),as(2,i),as(3,i)]^T                  !..GP-Werte

c     Ai      = [AAi(1,i),AAi(2,i),AAi(3,i)]^T               !..GP-Werte
c     Ai'     = [AAis(1,i),AAis(2,i),AAis(3,i)]^T            !..GP-Werte
c-----------------------------------------------------------------------
      USE eldata
      implicit double precision(a-h,o-z)
      dimension d(*),xl(ndm,*),ul(ndf,*),shp(2,*),eps(*)
     +         ,rin(3,3,*),rwn(3,3,*),pin(3,3,*),AAi(3,3),AAis(3,3)
     +         ,ai(3,3),as(3,3),xs(3),XXs(3),us(3),phi(3,3),phis(3,3)

      call pzero (XXs,3)
      call pzero (xs,3)
      call pzero (AAi,3*3)
      call pzero (AAis,3*3)
      call pzero (ai,3*3)
      call pzero (as,3*3)
      call pzero (us,3)
      call pzero (phi,3*3)
      call pzero (phis,3*3)

      do i = 1,3
        do j = 1,nel
          xs(i)  = xs(i)  + (xl(i,j)+ul(i,j))*shp(1,j)
          XXs(i) = XXs(i) + (xl(i,j)        )*shp(1,j)
        enddo
      enddo

      do i = 1,3
        do j = 1,3
          do k = 1,nel
            ai(i,j) =  ai(i,j) + rwn(i,j,k)*shp(2,k)
            as(i,j) =  as(i,j) + rwn(i,j,k)*shp(1,k)
            AAi(i,j) =  AAi(i,j) + rin(i,j,k)*shp(2,k)
            AAis(i,j) = AAis(i,j) + rin(i,j,k)*shp(1,k)
          enddo
        enddo
      enddo

      ilin=d(13)
      if (ilin.eq.0) then
        do i=1,3
          do k=1,nel
            us(i) = us(i)+ul(i,k)*shp(1,k)
            do j=1,3
              phi(i,j)  = phi(i,j)  + pin(i,j,k)*shp(2,k)
              phis(i,j) = phis(i,j) + pin(i,j,k)*shp(1,k)
            enddo
          enddo
        enddo 

        eps(1) = dot(us,AAi(1,1),3) + dot(XXs,phi(1,1),3)
        eps(2) = dot(us,AAi(1,2),3) + dot(XXs,phi(1,2),3)
        eps(3) = dot(us,AAi(1,3),3) + dot(XXs,phi(1,3),3)
        eps(4) = dot(phis(1,2),AAi(1,3),3) + dot(AAis(1,2),phi(1,3),3)
        eps(5) = dot(phis(1,3),AAi(1,1),3) + dot(AAis(1,3),phi(1,1),3)
        eps(6) = dot(phis(1,1),AAi(1,2),3) + dot(AAis(1,1),phi(1,2),3)

        call matcop(XXs,3,1, xs)
        call matcop(AAi,3,3, ai)
        call matcop(AAis,3,3, as)
      else
        eps(1) = dot(xs(1),ai(1,1),3)   - dot(XXs(1),AAi(1,1),3)
        eps(2) = dot(xs(1),ai(1,2),3)   - dot(XXs(1),AAi(1,2),3)
        eps(3) = dot(xs(1),ai(1,3),3)   - dot(XXs(1),AAi(1,3),3)
        eps(4) = dot(as(1,2),ai(1,3),3) - dot(AAis(1,2),AAi(1,3),3)
        eps(5) = dot(as(1,3),ai(1,1),3) - dot(AAis(1,3),AAi(1,1),3)
        eps(6) = dot(as(1,1),ai(1,2),3) - dot(AAis(1,1),AAi(1,2),3)
      endif

      return
      end
c
      subroutine bmat15(d,inode,ir0,p,s,shp,ai,as,xs,rwn,hwn,dwn,swn,
     +                  twh,sig,b,btd,dmat,dxl,nst,isw)
c-----------------------------------------------------------------------
c.... B-matrix, residuum vector p, Bt*D, Kg (diagonal part)
c-----------------------------------------------------------------------
      USE prlod
      implicit double precision (a-h,o-z)
      dimension d(*),p(*),s(nst,*),shp(2,*)
     +         ,ai(3,*),as(3,*),xs(*)
     +         ,rwn(3,3,*),hwn(3,3,*),dwn(3,*),swn(3,3,*)
     +         ,fb(3,9),sig(*),b(6,6,*),btd(6,*),dmat(6,*)
     +         ,twh(3,3,3,*),h1(3),h2(3),h3(3)
     +         ,xm1(3,3),xm2(3,3),xm3(3,3)

c.... parameter in matrix B
c       fb(i,j) = [ T1K^T x'    (*,1)--> epsilon
c                   T2K^T x'    (*,2)--> gamma_2
c                   T3K^T x'    (*,3)--> gamma_3
c                   T2K^T a3    (*,4)--> theta   (N_k')
c                   T3K^T a1    (*,5)--> kappa_2 (N_k')
c                   T1K^T a2    (*,6)--> kappa_3 (N_k') 
c                   T3K^T a2'   (*,7)--> theta   (N_k)
c                   T1K^T a3'   (*,8)--> kappa_2 (N_k)
c                   T2K^T a1'   (*,9)--> kappa_3 (N_k)  ]

      do i = 1,3 
        j = mod(i,3) + 1
        k = mod(j,3) + 1
        call mttmul (twh(1,1,i,inode),xs,3,3,1,fb(1,i))
        call mttmul (twh(1,1,j,inode),ai(1,k),3,3,1,fb(1,3+i))
        call mttmul (twh(1,1,k,inode),as(1,j),3,3,1,fb(1,6+i))
      enddo

c.... B-matrix 
c     strains: eps, gamma_2, gamma_3, theta, kappa_2, kappa_3
      do i  = 1,3
        i3 = i+3
        do j = 1,3
        j3 = j+3
         b(j,i,inode)   = ai(i,j)*shp(1,inode)
         b(j3,i,inode)  = 0.0d0
         b(j,i3,inode)  = fb(i,j)*shp(2,inode)
         b(j3,i3,inode) = fb(i,j3)*shp(1,inode) + fb(i,j+6)*shp(2,inode)
        enddo
      enddo

      do i = 1,3
c.... residual vector G = P - Bt*S
        ii = ir0 + i 
        p(ii) = p(ii) - dot(b(1,i,inode),sig(1),3) *dxl
        p(ii+3) = p(ii+3) - dot(b(1,i+3,inode),sig(1),6) *dxl
        do j = 1,6
c....     Bt * D
          btd(j,i) = dot(b(1,i,inode),dmat(1,j),3) *dxl
          btd(j,i+3) = dot(b(1,i+3,inode),dmat(1,j),6) *dxl
        enddo
      enddo

      ilin=d(13)
      if ((ilin.le.1).or.(isw.eq.4).or.(isw.eq.6).or.(isw.eq.13)) return

c.... Kg (diagonal part)
      do i=1,3
        h1(i)= (sig(1) * shp(2,inode) * xs(i)
     +       +  sig(5) * shp(2,inode) * as(i,3)    
     +       +  sig(6) * shp(1,inode) * ai(i,2))*dxl
        h2(i)= (sig(2) * shp(2,inode) * xs(i)
     +       +  sig(6) * shp(2,inode) * as(i,1)    
     +       +  sig(4) * shp(1,inode) * ai(i,3))*dxl
        h3(i)= (sig(3) * shp(2,inode) * xs(i)
     +       +  sig(4) * shp(2,inode) * as(i,2)    
     +       +  sig(5) * shp(1,inode) * ai(i,1))*dxl
      enddo
      
      k=inode
      call faca15(rwn(1,1,k),h1,hwn(1,1,k),dwn(1,k),swn(1,1,k),xm1)
      call faca15(rwn(1,2,k),h2,hwn(1,1,k),dwn(1,k),swn(1,1,k),xm2)
      call faca15(rwn(1,3,k),h3,hwn(1,1,k),dwn(1,k),swn(1,1,k),xm3)

      ir3 = ir0 + 3
      do i = 1,3
        ir2 = ir3 + i 
        do j= 1,i
          ic2 = ir3 + j 
          s(ir2,ic2) = s(ir2,ic2) + xm1(i,j) + xm2(i,j) + xm3(i,j)
        enddo
      enddo 

      return
      end
c
      subroutine stif15(inode,jnode,ir0,jc0,d,sig,shp,b,btd,s,
     +                  twh,ttt,dxl,nst)
c-----------------------------------------------------------------------
c.... compute tangential stiffness matrix
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  d(*),s(nst,*),sig(*),shp(2,*),btd(6,*),b(6,6,*)
     +          ,twh(3,3,3,*),ttt(3,3,6,*),guw(3),gwu(3),gww(3,2)

      ilin=d(13)
      if (ilin.eq.0) then
        call pzero(guw,3)
        call pzero(gwu,3)
        call pzero(gww,3*2)
      else
        do i = 1,3
         guw(i)   = sig(i) * shp(1,inode)*shp(2,jnode) *dxl
         gwu(i)   = sig(i) * shp(2,inode)*shp(1,jnode) *dxl
         gww(i,1) = sig(3+i) * shp(1,inode)*shp(2,jnode) *dxl
         gww(i,2) = sig(3+i) * shp(2,inode)*shp(1,jnode) *dxl
        enddo
      endif

      l = inode*(inode-1)/2 + jnode

      do i = 1,3
        ir1 = ir0+i
        ir2 = ir0+3+i
        do j= 1,3
          jc1 = jc0+j 
          jc2 = jc0+3+j

          s(ir1,jc1) = s(ir1,jc1) + dot(btd(1,i),b(1,j,  jnode),3)
          s(ir1,jc2) = s(ir1,jc2) + dot(btd(1,i),b(1,j+3,jnode),6)
     +               + guw(1)*twh(i,j,1,jnode) 
     +               + guw(2)*twh(i,j,2,jnode) 
     +               + guw(3)*twh(i,j,3,jnode)
          s(ir2,jc1) = s(ir2,jc1) + dot(btd(1,i+3),b(1,j,jnode),3)
     +               + gwu(1)*twh(j,i,1,inode) 
     +               + gwu(2)*twh(j,i,2,inode) 
     +               + gwu(3)*twh(j,i,3,inode)
          s(ir2,jc2) = s(ir2,jc2) + dot(btd(1,i+3),b(1,j+3,jnode),6)
     +               + gww(1,1)*ttt(i,j,1,l) + gww(1,2)*ttt(i,j,2,l)      
     +               + gww(2,1)*ttt(i,j,3,l) + gww(2,2)*ttt(i,j,4,l)      
     +               + gww(3,1)*ttt(i,j,5,l) + gww(3,2)*ttt(i,j,6,l)      
        enddo
      enddo

      return
      end
c
      subroutine dmat15 (d,dmat,dlx,scfy,scfz)
c-----------------------------------------------------------------------
c.... elasticity matrix for 3d beam 
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension d(*),dmat(6,6)
 
      E   = d(1)
      G   = d(2)
      A   = d(3)
      xI2 = d(4)
      xI3 = d(5)
      xI23= d(6)
      xIt = d(7)
      ys  = d(8)
      zs  = d(9)
      ym  = d(10)
      zm  = d(11)

      ish = d(15) 
      imat = d(17)

      EA = E*A
      GA = G*A
      EIy = E*xI2
      EIz = E*xI3
      scfy = 1.d0  
      scfz = 1.d0
      if(ish.ne.0)then
       scfy = scfy/(1.d0+scfy*dlx*dlx*GA/EIz/12.d0) 
       scfz = scfz/(1.d0+scfz*dlx*dlx*GA/EIy/12.d0) 
      endif

      if(imat.eq.0)then
      
       call pzero(dmat,36)
       xI = xIt + A*(ym*ym*scfz + zm*zm*scfy)
       
       dmat(1,1) =  EA
       dmat(1,5) =  EA*zs
       dmat(1,6) = -EA*ys
       dmat(2,2) = GA*scfy
       dmat(2,4) =-GA*zm*scfy
       dmat(3,3) = GA*scfz
       dmat(3,4) = GA*ym*scfz
       dmat(4,4) = G*xI
       dmat(5,5) = E*(xI2 + zs*zs*A)
       dmat(5,6) =-E*(xI23+ ys*zs*A)
       dmat(6,6) = E*(xI3 + ys*ys*A)
       
       do 10 i = 1,6
         do 10 j = i,6
10     dmat(j,i) = dmat(i,j)

      else
      
       scfy = dsqrt(scfy)
       scfz = dsqrt(scfz)
     
      endif
c
      return
      end
c
      subroutine faca15(an,x,hwn,wn,sn,aix)
c-----------------------------------------------------------------------
c     initial values M^=H^T*M(an,x)*H for geometrical matrix K_sigma
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension   an(3),x(3),hwn(3,3),wn(3),sn(3,3),aix(3,3)
     +           ,r(3),s(3),fh(3,3),one(3,3)
      data one   / 1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0 /

      call vecp (an,x, r)
      call mvmul(sn,r,3,3, s)
      theta = dsqrt(dot(wn,wn,3))
      if (dabs(theta) .lt. 1.0d-15)  then
         c = 1.d0/6.d0
      else
         c = (theta-dsin(theta))/(theta*theta*theta)
      endif
      
      c1 = dot(an,x,3)
      
      do 10 i = 1,3
      do 10 j = 1,3
       aix(i,j) = 0.5d0*(an(i)*x(j) + x(i)*an(j)) - c1*one(i,j)
10    continue
      
      c2 = c*dot(r,wn,3)

      call matmulf (aix,hwn, 3,3,3, fh)
      call mttmul (hwn, fh, 3,3,3, aix)

      do 20 i = 1,3
      do 20 j = 1,3
       aix(i,j) = aix(i,j)+0.5d0*(s(i)*wn(j) + wn(i)*s(j)) + c2*one(i,j)
20    continue

      return
      end
c
      subroutine triad15 (shp,xl,rin,xsi,phi,ndm)
c-----------------------------------------------------------------------
c     initial cartesian system 
c     rin(3,3)    initial local triad at node
c            e_1 in element direction
c            e_2 = -e_1 x  e_Z (=[0,0,1])
c            e_3 =  e_1 x  e_2
c     
c.... special case e_1 = +-e_Z   ==>   e_2 = -e_Y
c.... rotate local coordinate system with angle phi (grad)
c            (phi starts positiv from the 2-axis 
C                     in direction to the 3-axis
c            e_2R =  c * e_2  +  s * e_3
c            e_3R = -s * e_2  +  c * e_3
c-----------------------------------------------------------------------
      USE eldata
      implicit double precision (a-h,o-z)
      dimension  shp(2,*), xl(ndm,*), rin(3,3), dy(3), dz(3)

      CALL shape1D(shp,xsi,xl,dl,ndm,nel)

      do 100 i = 1,3
        rin(i,1) = 0.d0
         dz(i)   = 0.d0
        do 100 k = 1,nel
         rin(i,1) = rin(i,1) + xl(i,k)*shp(1,k)
100   continue

c.... e_Z = 0,0,1
      dz(3) = 1.d0

      call vnorm(rin(1,1),dxl)
      call vecp(dz,rin(1,1),rin(1,2))
      call vnorm(rin(1,2),dyl)

c.... special case e_1 = +-e_Z:
c                  e_2 = -e_Y
      if (dabs(dyl) .lt. 1.0d-10) then
         rin(2,2) = -1.d0
         call vnorm(rin(1,2),dyl)
      endif

      call vecp(rin(1,1),rin(1,2),rin(1,3))

c.... rotate local coordinate system with angle phi
      if (dabs(phi) .gt. 1.0d-10) then
        ra = phi*datan(1.d0)/45.d0
        cs = dcos(ra)
        sn = dsin(ra)
        do 200 i = 1,3
          dy(i)    =  cs*rin(i,2) + sn*rin(i,3)
          dz(i)    = -sn*rin(i,2) + cs*rin(i,3)
          rin(i,2) = dy(i)
          rin(i,3) = dz(i)
200     continue
      endif

      return
      end
c
      subroutine updn15(d,xl,ul,shp,rin,rwn,hwn,dwn,swn,pin,twh,ttt,
     +                  ndm,ndf)
c-----------------------------------------------------------------------
c     update of nodal rotations
c     rin    initial local cartesian system   at element nodes
c     rwn    current nodal basis R = R(omega)      - " -
c-----------------------------------------------------------------------
      USE cdata
      USE dirdat
      USE eldata
      implicit double precision (a-h,o-z)
      dimension d(*),xl(ndm,*),ul(ndf,*),shp(2,*),rh1(3,3),rh2(3,3)
     +         ,rin(3,3,*),rwn(3,3,*),hwn(3,3,*),dwn(3,*),swn(3,3,*)
     +         ,wnk(3,3,3),pin(3,3,*),twh(3,3,3,*),ttt(3,3,6,*)

      alpha = d(12)
      ilin  = d(13)
c.... initial nodal cartesian system
      do k = 1,nel
        xsi = 2.d0/(nel-1) * (k - (nel+1)/2.d0)
        call triad15 (shp,xl(1,1),rin(1,1,k),xsi,alpha,ndm)
c...    update of current nodal basis
        do  i = 1,3
         dwn(i,k) = ul(i+3,k)        
         if(ilin.eq.3)dwn(i,k) = ul(i+3,nen+k)
        enddo 

        if (ilin.eq.0) then
          do i=1,3
            call vecp(dwn(1,k),rin(1,i,k),pin(1,i,k))
          enddo
        endif

        call updr15(dwn(1,k),rh1,       ilin,1)
        call updr15(dwn(1,k),hwn(1,1,k),ilin,2)
        call updr15(dwn(1,k),swn(1,1,k),ilin,3)
        if(ilin.eq.3)then
         call pdirec2 (basea,rh2,k,n,nen,2,1)
         call matmulf (rh1,rh2,3,3,3,rwn(1,1,k))
         call pdirec2 (basea,rwn(1,1,k),k,n,nen,1,2)
        else
         call matmulf(rh1,rin(1,1,k),3,3,3,rwn(1,1,k))
        endif

c...    W_ni,    TmI = W_mI^T * H_I   m=1,2,3
        do i = 1,3 
         if (ilin.lt.2)call skew (wnk(1,1,i),rin(1,i,k))
         if (ilin.ge.2)call skew (wnk(1,1,i),rwn(1,i,k))
         call mttmul (wnk(1,1,i),hwn(1,1,k),3,3,3,twh(1,1,i,k))
        enddo

c....   TTT(*,*,KJ,l) = TK^T * TJ  for all node combinations l 
c              (KJ)=(23,32,31,13,12,21)
       do j = 1, k
        l = k*(k-1)/2 + j
         call mttmul (twh(1,1,2,k),twh(1,1,3,j), 3,3,3, ttt(1,1,1,l))
         call mttmul (twh(1,1,3,k),twh(1,1,2,j), 3,3,3, ttt(1,1,2,l))
         call mttmul (twh(1,1,3,k),twh(1,1,1,j), 3,3,3, ttt(1,1,3,l))
         call mttmul (twh(1,1,1,k),twh(1,1,3,j), 3,3,3, ttt(1,1,4,l))
         call mttmul (twh(1,1,1,k),twh(1,1,2,j), 3,3,3, ttt(1,1,5,l))
         call mttmul (twh(1,1,2,k),twh(1,1,1,j), 3,3,3, ttt(1,1,6,l))
       enddo 

      enddo

      return
      end
c
      subroutine updr15 (dwn,rw,ilin,ifl)
c-----------------------------------------------------------------------
c      
c     ifl = 1:  current basis  R(om)  
c           2:                 H(om)
c           3:                 S(om)
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension   dwn(3),rw(3,3),om(3,3),om2(3,3), one(3,3) ,c(3)
      data one   / 1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0 /

      c(1) = 1.d0 
      c(2) = 0.d0
      c(3) = 0.d0
      if(ilin.eq.1.and.ifl.eq.1)c(2) = 1.d0
      
      theta = dsqrt(dot(dwn,dwn,3))
      thet2 = theta*theta

      call skew(om,dwn)
      call matmulf (om,om, 3,3,3, om2)

      if(ilin.ge.2)then
      
       if (ifl .eq. 1) then                       ! R(om)
         if (dabs(theta) .lt. 1.0d-10) then       ! theta -> 0
           c(1) = 1.d0 
           c(2) = 1.d0
           c(3) = 1.d0/2.d0
         else
           c(1) = 1.d0 
           c(2) = dsin(theta)/theta
           c(3) = (1.d0-dcos(theta))/thet2
         endif
       
       else if (ifl .eq. 2) then                  ! H(om)  
         if (dabs(theta) .lt. 1.0d-10) then       ! theta -> 0
           c(1) = 1.d0 
           c(2) = 1.d0/2.d0
           c(3) = 1.d0/6.d0
         else
           c(1) = 1.d0 
           c(2) = (1.d0 -dcos(theta))/thet2
           c(3) = (theta-dsin(theta))/(theta*thet2)
         endif
       
       else if (ifl .eq. 3) then                  ! Si
         if (dabs(theta) .lt. 1.0d-10) then       ! theta -> 0
           c(1) = -1.d0/6.d0
           c(2) =  1.d0/12.d0
           c(3) = -1.d0/60.d0
         else
           c(1) = -(theta-dsin(theta))/(theta*thet2)
           c(2) = -( dsin(theta)/(theta*thet2)
     1            -2.d0*(1.d0 -dcos(theta))/(thet2*thet2) )
           c(3) = ( (1.d0 -dcos(theta))/(thet2*thet2)
     1            -3.d0*(theta-dsin(theta))/(theta*thet2*thet2) )
         endif
       endif
            
      endif      

c
      do i = 1,3
       do j = 1,3
         rw(i,j) = c(1)*one(i,j) + c(2)*om(i,j) + c(3)*om2(i,j)
       enddo
      enddo 

      return
      end
c
      subroutine qloa15(d,q,xl,rwn,hwn,dwn,swn,pg,wg,ngb,p,s,numel,
     1                   nst,ndm,ndf)
c-----------------------------------------------------------------------
c     calculate element load- and stiffness part
c       q(n,1): element n, y_p: load coordinates
c       q(n,2): element n, z_p: load coordinates
c       q(n,3): element n, q_X: load in global direction X
c       q(n,4): element n, q_Y: load in global direction Y
c       q(n,5): element n, q_Z: load in global direction Z
c-----------------------------------------------------------------------
      USE eldata
      USE prlod
      USE qload
      implicit double precision(a-h,o-z)
      dimension d(*),xl(ndm,*),shp(2,5),p(*),s(nst,*),pg(*),wg(*)
     +         ,rwn(3,3,*),hwn(3,3,*),dwn(3,*),swn(3,3,*)
     +         ,ap(3,3,5),dp(3,5),ql(3),rm(3,5),rm0(3,5),q(numel,10)

c.... load coordinates and load vector q (=const!)
      ilin=d(13)
      yp   = q(n,1)
      zp   = q(n,2)
      ql(1)= q(n,3)*propq
      ql(2)= q(n,4)*propq
      ql(3)= q(n,5)*propq
      yzp  = dsqrt((yp*yp+zp*zp)/d(3))
      if(yzp.gt.1.d-7.and.ilin.eq.1)stop 'stop in qloa15'

c.... distance dp and vector m = dp x q
      do k = 1,nel
        do i = 1,3
          dp(i,k) = yp*rwn(i,2,k) + zp*rwn(i,3,k)
        enddo
        call vecp(dp(1,k),ql,rm0(1,k))
        call mttmul (hwn(1,1,k),rm0(1,k),3,3,1,rm(1,k))
c....   initial values Ai(x)_n for geometrical matrix K_sigma
        call faca15(dp(1,k),ql,hwn(1,1,k),dwn(1,k),swn(1,1,k),ap(1,1,k))
      enddo

c.... loop over gauss points
      do l = 1,ngb
        xsi = pg(l)
        CALL shape1D(shp,xsi,xl,dl,ndm,nel)
        dxl = wg(l)*dl
c....   loop over rows
        do inode = 1,nel
c....     compute element load vector
          do i = 1,3
            ii = i + ndf*(inode-1)
            p(ii)   = p(ii)   +  ql(i)      * shp(2,inode)*dxl
            p(ii+3) = p(ii+3) + rm(i,inode) * shp(2,inode)*dxl
          enddo
          if (ilin.ge.2) then
c....       compute tangent stiffness, only for inode=jnode
            do i = 1,3
              ir = i+3+ndf*(inode-1)
              do j = 1,3
                ic = j+3+ndf*(inode-1)
                s(ir,ic) = s(ir,ic) - ap(i,j,inode)*shp(2,inode)*dxl
              enddo
            enddo
          endif
        enddo
      enddo
      return
      end
c
      subroutine pout15 (d,q,xl,ul,p,rin,rwn,hwn,pg,wg,ngb,ndf,ndm,isw,
     +                   numel)
c-----------------------------------------------------------------------
c     (A) get external loads (qloa)
c     (B) update residuum vector (add the load contributions)
c     (C) print local stress resultants in list and Ofile  (isw=4)
c     (D) plot local stress resultants on (un)deformed mesh (isw=13)
c         (with nx points)
c     iPlo=0   : stress resultants related to reference configuration
c     iPlo=else: stress resultants related to current   configuration
c-----------------------------------------------------------------------
      USE bdata
      USE eldata
      USE iofile
      USE pdata10
      USE prlod
      USE qload
      implicit double precision (a-h,o-z)
      parameter (nx=11)
      dimension  d(*),xl(ndm,*),ul(ndf,*),p(*),shp(2,5),rin(3,3,*)
     +          ,rwn(3,3,*),pl(6,5),ps(3),px0(6),px1(6),px2(6),px(6,nx)
     +          ,pxl(6,nx),rip(3,3,nx),xp0(3),xp(3,nx),pg(*),wg(*)
     +          ,q(numel,10),dp(3,5),rm(3,5),rm0(3,5),hwn(3,3,*)
      alpha=d(12)
      ilin=d(13)
      iPlo=d(14)
      iln =0
      if(ilin.ne.0)iln=1
      call pzero (rip,3*3*nx)
      call pzero (pl,6*nel)
      call pzero (pxl,6*nx)
      call pzero(ps,3)
      call pzero(rm,3*5) 
      yp = 0.d0
      zp = 0.d0

c.... (A) external loads for the current element (only constant loads!)
c     load coordinates yp and zp, and load vector q
      if(mqloa.ne.1)then
        yp    = q(n,1)
        zp    = q(n,2)
        ps(1) = q(n,3) * prop
        ps(2) = q(n,4) * prop
        ps(3) = q(n,5) * prop
      endif

c.... (B) compute the residual forces at end-nodes:
c     distance vector dp and vector , m = dp x q  , q (=const!)
      if(dot(ps,ps,3).gt.1.d-20) then
        do k=1,nel
          do i=1,3
            dp(i,k) = yp*rwn(i,2,k) + zp*rwn(i,3,k)
          enddo
          call vecp(dp(1,k),ps,rm0(1,k))
          call mttmul (hwn(1,1,k),rm0(1,k),3,3,1,rm(1,k))
        enddo
c      
        do l = 1,ngb
          xsi=pg(l)
          CALL shape1D(shp,xsi,xl,dl,ndm,nel)
          dxl=wg(l)*dl
          do k = 1,nel
            do i = 1,3
              j = i + ndf*(k-1)
              p(j)   = p(j)   + ps(i)   * shp(2,k)*dxl
              p(j+3) = p(j+3) + rm(i,k) * shp(2,k)*dxl
            enddo
          enddo
        enddo
      endif

c.... (C) print local stress resultants:
      if (isw.eq.4) then
c       transformation global to local direction at nodes
        do k = 1,nel
          do i = 1,3
            do j = 1,3
              jj = j + ndf*(k-1)
              if (iPlo.eq.0) then
                pl(i,  k) = pl(i,  k) - rin(j,i,k)*p(jj)
                pl(i+3,k) = pl(i+3,k) - rin(j,i,k)*p(jj+3)
              else
                pl(i,  k) = pl(i,  k) - rwn(j,i,k)*p(jj)
                pl(i+3,k) = pl(i+3,k) - rwn(j,i,k)*p(jj+3)
              endif
            enddo
          enddo
        enddo
c       print
        mct = mct - 1
        if (mct.le.0) then
          if (iPlo.eq.0) then
            if(ior.lt.0) write(*  ,4000) o,head
                         write(iow,4000) o,head
          else
            if(ior.lt.0) write(*  ,4010) o,head
                         write(iow,4010) o,head
          endif
          mct = 50
        endif
                     write(iow,4001) n,ma,(-pl(i,  1), i=1, 6)
                     write(iow,4002)      ( pl(i,nel), i=1, 6)
        if(ior.lt.0) write(*  ,4001) n,ma,(-pl(i,  1), i=1, 6)
        if(ior.lt.0) write(*  ,4002)      ( pl(i,nel), i=1, 6)
4000    format(a1,20a4/ 
     1      10x,'3-D beam element: STRESS RESULTANTS  ',/,
     2  '  el  mat','   H_1    ',1x,'   V_2    ',1x,'   V_3    ',1x,
     3              '   M_1    ',1x,'   M_2    ',1x,'   M_3    ')
4010    format(a1,20a4/ 
     1      10x,'3-D beam element: STRESS RESULTANTS  ',/,
     2  '  el  mat','   N_1    ',1x,'   Q_2    ',1x,'   Q_3    ',1x,
     3              '   M_1    ',1x,'   M_2    ',1x,'   M_3    ')
4001    format(1x,i3,2x,i3,6(1x,g10.4))
4002    format(9x,6(1x,g10.4))
        return
      endif

c.... (D) plot local stress resultants on (un)deformed mesh (isw=13)
c        compute the forces within the current element at nx-points
c        [distributed load possible! (ps(xsi))]
c     initialize start point at first node on element
      do i = 1,3
        do j = 1,3
           if (iplo.eq.0) rip(i,j,1) = rin(i,j,1)
           if (iplo.ne.0) rip(i,j,1) = rwn(i,j,1)
        enddo
      enddo
c
      xsi = -1.d0
      call shape1D(shp,xsi,xl,dl,ndm,nel)
      do i = 1,3
        px(i,  1) = p(i)
        px(i+3,1) = p(i+3)
        px1(i)    = ps(i)*dl
        px1(3+i)  = rm(i,1)*dl
        xp(i,1)   = xl(i,1)+ul(i,1)*iln
      enddo
c     next nx-1 points at element
      do nk = 1,nx-1
        xsi = -1.d0 + 2.d0/dble(nx-1)*nk
        CALL shape1D(shp,xsi,xl,dl,ndm,nel)
       if (iplo.eq.0) then
         call triad15(shp,xl,rip(1,1,nk+1),xsi,alpha,ndm)
       else
         do i = 1,3
           do j = 1,3
             do k = 1,nel
               rip(i,j,nk+1) = rip(i,j,nk+1)+shp(2,k)*rwn(i,j,k)
             end do ! k
           end do ! j
         end do ! i
       end if !iplo
c
        do i = 1,3
          px0(3+i) = 0.d0
           do k = 1,nel
             px0(3+i) = px0(3+i) + shp(2,k)*rm(i,k)
           end do
        end do
c
        do i = 1,3
          px2(i)   =  ps(i)*dl
          px2(3+i) =  px0(3+i)*dl
          px0(i)   = (px1(i) + px2(i))/2.d0
          px0(3+i) = (px1(3+i) + px2(3+i))/2.d0
          px(i,nk+1) = px(i,nk)-px0(i)*(2.d0/dble(nx-1))
          xp(i,nk+1) = 0.0d0
          do k = 1,nel
            xp(i,nk+1)= xp(i,nk+1) + (xl(i,k)+ul(i,k)*iln)*shp(2,k)   
          enddo
          xp0(i) = xp(i,nk) - xp(i,nk+1)
        enddo
        call vcross (xp0,px(1,nk), px(4,nk+1))
        call vcross (xp0,px0,px1)
        do i = 1,3
          px(3+i,nk+1)=px(3+i,nk+1)+px(3+i,nk)
     +                             -px1(i)*(1.d0/dble(nx-1))
     +                             -px0(3+i)*(2.d0/dble(nx-1))
          px1(i)   = px2(i)
          px1(3+i) = px2(3+i)
        enddo
      enddo
c     transformation to local direction 
      do k = 1,nx
        do i = 1,3
          do j = 1,3
            pxl(i,  k) = pxl(i,  k) + rip(j,i,k)*px(j,  k)
            pxl(i+3,k) = pxl(i+3,k) + rip(j,i,k)*px(j+3,k)
          enddo
        enddo
      enddo
c     plot
      call plot15(xp,pxl,rip,ndm,nx)
      return
      end
c
      subroutine plot15 (xp,p,rip,ndm,nx)
c-----------------------------------------------------------------------
c     plot forces at nx-points 
c-----------------------------------------------------------------------
      USE eldata
      USE pdata10
      USE pltran
      implicit double precision (a-h,o-z)
      logical zoom
      dimension  p(6,nx),xp(3,*),rip(3,3,*)
     1          ,ixl(5),ixld(4),xn(3),xd(3),xll(3,4)
      data ixl /1,2,3,4,1/, ixld /1,2,3,1/, eps/1.0e-5/

      klayf=1
      mfp = iabs(nfp)
      if(mfp.lt.1.or.mfp.gt.6) return

      do 100 inode = 1, nx-1
      inode1 = inode
      inode2 = inode+1
      s1 =  p(mfp,inode1)
      s2 =  p(mfp,inode2)
      if(dabs(s1).lt.eps) s1 = 0.d0
      if(dabs(s2).lt.eps) s2 = 0.d0
c.... check max/min values
      if(flfp) then
        ccfp  = max(dabs(s1),dabs(s2))
        ccfp1 = max(s1,s2)
        ccfp2 = min(s1,s2)
        xmaxf = max(xmaxf,ccfp1)
        xminf = min(xminf,ccfp2)
        cfp   = max(cfp,ccfp)
      else
       if(dabs(s1).lt.eps.and.dabs(s2).lt.eps) goto 100
c....  normal vector and length of dx
        ifact = sign(1,ifor)     ! ifor=+-12,+-13
        icolm = iabs(ifor)-10
        do 10 i = 1,3
          xn(i) = rip(i,icolm,inode1)*ifact  
          xd(i) = xp(i,inode2) - xp(i,inode1)
10      continue
        call vnorm (xd,dl)

c....  Berechne Vorzeichen
       sm = (s1+s2)*0.5d0
       call pppcolf(sm)
c....  Plotte Schnittgroesse
       if(s1.gt.0.d0.and.s2.gt.0.d0) then
         goto 210
       elseif(s1.lt.0.d0.and.s2.lt.0.d0) then
         goto 210
       else
c....  Vorzeichenwechsel -> rechne zwischen (0,s1) und (dl,s2) lokal
          x1 = 0.d0
          x2 = dl
          y1 = s1 
          y2 = s2 
c.....    gerade y = y1 + a1*x   mit a1 = (y2-y1)/dl 
          a1  = (y2-y1)/dl
c.....    Durchstosspunkt
          xdp = -y1/a1
          if(xdp.lt.x1 .or. xdp.gt.x2) goto 200  ! ausserhalb
c....     plotte dreieck 1
          sm = s1*0.5d0
          call pppcolf(sm)
          call pzero(xll,12)
          do 20 k = 1,3
            xd(k)    = xp(k,inode1)*(1.d0-xdp/dl) + xp(k,inode2)*xdp/dl
            xll(k,1) = xp(k,inode1)
            xll(k,2) = xd(k)
            xll(k,3) = xp(k,inode1) + xn(k)*s1*cfp
20        continue
          call plxtrn(xll,tra,vr,3,4)
          if(zoom(xll,3,3)) call plot9s(ixld,xll,ndm,3)
c....     plotte dreieck 2
          sm = s2*0.5d0
          call pppcolf(sm)
          call pzero(xll,12)
          do 30 k = 1,3
            xll(k,1) = xd(k)
            xll(k,2) = xp(k,inode2)
            xll(k,3) = xp(k,inode2) + xn(k)*s2*cfp
30        continue
          call plxtrn(xll,tra,vr,3,4)
          if(zoom(xll,3,3)) call plot9s(ixld,xll,ndm,3)
          goto 200
        endif
c....   plotte trapez
210     do 40 k = 1,3
          xll(k,1) = xp(k,inode1)
          xll(k,2) = xp(k,inode2)
          xll(k,3) = xp(k,inode2) + xn(k)*s2*cfp
          xll(k,4) = xp(k,inode1) + xn(k)*s1*cfp
40      continue
        call plxtrn(xll,tra,vr,3,4)
        if(zoom(xll,3,4)) call plot9s(ixl,xll,ndm,4)
200     continue
      endif
100   continue
      return
      end