      subroutine elmt29(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c-----------------------------------------------------------------------
c
c.... geometrically nonlinear four node shell element 
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
c     1 d(01) = igeo 1=arb., 2= arb.+given t_1, 3-17 see macro BASE
c     2 d(02) = <value 1 for igeo> 
c     3 d(03) = <value 2 for igeo> 
c     4 d(04) = <value 3 for igeo> 
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
c     9 d(19) = shear corr. factor: (<0 = FE-SCF, else > 0)
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
c.... output of stresses see Subroutine plsn29 
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
c     (c)  FG/WW  2015
c-----------------------------------------------------------------------
      USE bdata       
      USE cdata 
      USE dirdat 
      USE eldata 
      USE hdata 
      USE iofile 
      USE qload       
      USE plslay 
      USE strnam 
      USE plodfu 
      implicit double precision (a-h,o-z)

c
      parameter (neas_max=14, ngeas_max=24 + neas_max)
c
      dimension xl(ndm,*),tl(*),d(*),ul(ndf,*),s(nst,*),p(*),ix(*),
     1          st(ngeas_max,ngeas_max),pt(ngeas_max),h1(*),h2(*),h3(*),   
     2          shp(3,4),shp1(3,4),sg(16),tg(16),wg(16),t0(3,3),
     3          dmat(8,8),sig(8),eps(8),epse(8),b(8,6,4),btd(8,6),
     4          xs0(2,2),sx0(2,2),sx(2,2),xu(3,3),
     5          rd(3,3),rin(3,3,4),twn(3,3,4),rwn(3,3,4),w(3,3,4),
     6          x0(3,3),r0(3,3),te0(5,5),ts0(5,5),qs(2),siga(9),
     7          gam(4),ran(3,4),xan(3,4),bml(3,2,4),ddn(3,4),
     8          pgz(5),wgz(5),sigp(6),epstr(5),cmat(5,5),
     9          gxy(8,neas_max),gtd(neas_max,8),eas(8)        
c
c.... go to correct array processor
      go to(1,2,3,3,5,3,7,3,2,2,2,2,2,2,2,2,2,18,2,2,2,22), isw

c.... input element properties
1      if(ior.lt.0) write(*,1000)
1000  format(' Input base system:',/
     +       ' igeo,add(1-3) >',$)
      call dinput(d,4)
c
      if(ior.lt.0) write(*,1001)
1001  format(' Input geometry etc.:',/ 
     +       ' nlay,ngp,h_s,z_h,imat,ilin,ityp,neas,stab >',$)
      call dinput(d(11),9)
      if(dabs(d(19)).lt.1.d-12) d(19)=1.d0 ! default SCF
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
     +         (d(i),i=1,4),(d(i),i=11,17),nm,nb,ns,d(19),(d(i),i=25,32) 
                   write(iow,1004) 
     +         (d(i),i=1,4),(d(i),i=11,17),nm,nb,ns,d(19),(d(i),i=25,32) 
1004  format(/5x,'Nonlinear four-node shell element ',
     +           ' element+material data:',//,
     + 5x,'igeo 1=arb. others see macro BASE... ',f8.0,/,
     + 5x,'additional value 1 for igeo..... ',g12.4,/,
     + 5x,'additional value 2 for igeo..... ',g12.4,/,
     + 5x,'additional value 3 for igeo..... ',g12.4,//,
     + 5x,'no. of layers nlay ................. ',f8.0,/,
     + 5x,'no. of GP/layer (def.=2)............ ',f8.0,/,
     + 5x,'shell thickness h_s ............ ',g12.4,/,
     + 5x,'zeta_h ......................... ',g12.4,/,
     + 5x,'material type imat.................. ',f8.0,/,
     + 5x,'ilin 0=linear 1=moderate 2=finite... ',f8.0,/,
     + 5x,'ityp 5=5dofs local 6=6dofs global... ',f8.0,/,
     + 5x,'eas parameter nm/nb/ns each 0/4/5/7  ',i4,'/',i1,'/',i1,/,
     + 5x,'shear correction factor (<0 = FE-scf)',g12.4,//,
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
      call matread29(d)
c
c.... initialize h-array    
c.... for eas:   kcr+kcc(neas*ngeas), rc(neas), alpha(neas)
      nlay = d(11)
      ngz  = d(12)                  ! integration points for each layer
      imat = d(15)
      neas = nm+nb+ns
      ngeas = 4*ndf + neas

      ngb = 2                       ! integration points in plane
      npara = neas * (ngeas+2)      ! kcr+kcc(neas*ngeas), rc(neas), alpha(neas) 
      nhdr  = npara + 3*4           ! + 3 rotational parameter at each node
      nh    = 8                     ! 1-6: E_p, 7: e_v, 8: E_33
      if(imat.eq.1.and.dabs(d(42)).le.1.d-18) nh = 0
      nh1   = nh*ngb*ngb*ngz*nlay  ! total number per element in h1/h2
      nh3   = nhdr                 ! total number per element in h3
c
2     return
3     igeo  = d(1)
      nlay  = d(11)
      imat  = d(15)
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
      ielas = 0 
      if(imat.eq.1.and.dabs(d(42)).lt.1.d-18)ielas = 1
      if(ielas.eq.1) nh = 0

      maxlay = ngz*nlay
c      if(klay.gt.maxlay) klay = maxlay
c.... names for stresses
      if(isw.eq.8) then
       if(klay.eq.0) then
        call plsn29 (0)
       else
        call plsn29 (1)
       end if 
      end if
c      
c.... storage in h3-array: 
c     omega(3*4),  sc(neas*ngeas), pc(neas), alpha(neas)
c      
      nhom  = 1
      nsc   = nhom + 3*4
      npc   = nsc + neas*ngeas
      nalp  = npc + neas
c
      call pzero (pt,ngeas_max)
      call pzero (st,ngeas_max*ngeas_max)                    
      call pzero (eas,8)
      call pgauss (ngb,lint,sg,tg,wg)
      call gaus1D (ngz,pgz,wgz)
      call rjac29(d,xl,t0,xs0,sx0,ts0,te0,detj0,isw)
      call shearfac29(d,detj0,xl,cappa,wcappa)
      if(ielas.eq.1)call dmate29 (d,cappa,dmat)

      if(ldir.eq.0) then 
        call triad29 (igeo,xl,d,rin)            
      else if(ldir.eq.1) then          ! read from basea=m(mdir)
        do k = 1,4 
          if (knode.eq.numnp) call pdirec1(ix,basea,rin(1,1,k),k,2,1)  
          if (knode.eq.(nen*numel)) 
     +       call pdirec2 (basea,rin(1,1,k),k,n,nen,2,1)
        end do 
      end if

      if(neas.gt.0)
     1    call decond(h3(nsc),h3(npc),h3(nalp),ul(1,2*nen+1),neas,4*ndf)

      call updn29 (xl,ul,h3,twn,rwn,ddn,rin,ran,
     1             xan,gam,bml,w,ix,ndf,ndm,nen,ityp,ilin)
c
c
c.... loop over gauss points
      nhpl = 1
      do 300 l = 1,lint
       call pzero (eps,8)
       call pzero (epse,8)
       xsi = sg(l)
       eta = tg(l)
       call shap29(xsi,eta,xl,t0,sx,shp,shp1,xsj)
       da = xsj*wg(l)
       rj = detj0/xsj
       if(neas.gt.0)call easstr29 (xsi,eta,te0,rj,neas_max,neas,
     1                             nm,nb,ns,gxy,eas,h3(nalp))

       call stra29 (xl,ul,shp,xsi,eta,eas,eps,sx,gam,rin,rwn,
     1              ddn,xu,x0,rd,r0,ndm,ndf,ilin)

       call strt29 (d,aqloa,numel,n,mqloa,propq,eps,epse,isw)
        
       if(ielas.eq.1)then
        call mvmul(dmat,epse,8,8,sig)
        valuse1 = valuse1 + 0.5d0*dot(sig(1),epse(1),8)*da
        call mvmul(dmat,eps,8,8,sig)
       else

c....   integration through the thickness

        call pzero (sig,8)
        call pzero (dmat,8*8)

        do 361 ilay = 1,nlay
           dz = d(14) + (ilay-1)*dh
           do 360 igz = 1,ngz
              wz   = wgz(igz)*dh5             
              zeta = pgz(igz)
              zs   = dz + (1.d0+zeta)*dh5
              call kine29(eps,x0,r0,zs,wcappa,epstr,dets)
              wz = wz*dets   
c
              call zerostr29(h1(nhpl),h2(nhpl),d,epstr,x0,
     +                       sigp,cmat)

c.....        print/plot stresses in layer klay and GP mlay
              if (isw.eq.4 .or. isw.eq.8) then  
               if (ilay.eq.klay.and.igz.eq.mlay) then 
                 call trasig29 (sigp,t0,shp,rin)               
                 
                 if (isw.eq.4) call pout29 (xl,shp,sig,sigp,
     +                              h2(nhpl),ndm,maxlay)  
                 if (isw.eq.8) then 
                   call plot29 (d,ix,strea,strea(1+numnp),
     +               shp,sig,sigp,h2(nhpl),eps,numnp,maxlay,klay,da)
                 end if 
               end if
              end if

c.....        compute material matrix and stress resultants
              call dmat29 (dmat,cmat,sig,sigp,wz,zs,wcappa,cappa)
              nhpl = nhpl + nh
360        continue
361     continue
       end if

       if (isw.eq.4 .or. isw.eq.8) then
        if (klay.eq.0) then         ! print/plot stress resultants
c....    transformation
         call transsig29(t0,xu,rd,rin,shp,sig,siga)
         if (isw.eq.4) call pout29 (xl,shp,siga,sigp,
     +                              h2,ndm,maxlay)
         if (isw.eq.8.and.iplma(ma).ne.0) then 
            call plot29 (d,ix,strea,strea(1+numnp),shp,siga,
     +                   sigp,h2,eps,numnp,maxlay,klay,da)
         end if
        end if
       end if
c         
       if (isw.eq.4.or.isw.eq.8) go to 300
       if(neas.gt.0)call easmat29 (gxy,gtd,pt,st,sig,dmat,b,da,neas,
     1                             neas_max,ngeas_max,ndf,0,0,1)
c....  loop over rows, residual
       ir0 = 0
       do 310 inode = 1,4
        call bmat29 (inode,ir0,b,btd,dmat,st,pt,sig,qs,rwn,h3,
     1   w,twn,bml,d,shp1,shp,sx,t0,ran,xan,xu,rd,da,ngeas_max,ndf,ilin)
        if(neas.gt.0)call easmat29 (gxy,gtd,pt,st,sig,dmat,b,da,neas,
     1                           neas_max,ngeas_max,ndf,inode,ir0,2)

c....   loop over columns (symmetry noted), stiffness matrix
        jc0 = 0       
        do jnode = 1,inode
         call stif29 (inode,jnode,ir0,jc0,b,btd,st,sig,qs,shp,shp1,
     1                w,ngeas_max,ndf,da,ilin,isw)
         jc0 = jc0 + ndf
        end do
310     ir0 = ir0 + ndf
300   continue
         if (isw.eq.4) go to 4
         if (isw.eq.8) go to 8
c
      call easend29 (st,pt,s,p,h3(nsc),h3(npc),nst,ndf,ngeas_max,neas)

c      if(n.eq.1)call mprint(s,ndf,ndf,nst,'s29 ')            
     
      return
c
c.... output stresses   =>  isw 3!
4     return
c.... compute lumped mass matrix 
c     open: 5/6 dof?, value for dof 6? 
5     l = 2
      call pgauss(l,lint,sg,tg,wg)
      call rjac29(d,xl,t0,xs0,sx0,ts0,te0,detj0,isw)
      hs = d(13)
      dh = hs*hs/12.d0
      do 50 l = 1,lint
        xsi = sg(l)
        eta = tg(l)
        call shap29(xsi,eta,xl,t0,sx,shp,shp1,xsj)
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
c.... calculate initial nodal cartesian basis
18    continue
      igeo = d(1)
      call triad29 (igeo,xl,d,rin)

      if (knode.eq.(nen*numel)) then
c....  BASE,1-2       
       do k = 1,4 
        call pdirec2 (basea,rin(1,1,k),k,n,nen,1,1)
        call pdirec2 (basea,rin(1,1,k),k,n,nen,1,2)
       end do 
      else if (knode.eq.numnp) then
c....  BASE,0       
       do k = 1,4 
        call pdirec1(ix,basea,rin(1,1,k),k,1,1)
        call pdirec1(ix,basea,rin(1,1,k),k,1,2)
       end do 
      end if
      return
c.... compute the surface tractions  (Input: QLOA)  temperature loading
c.... compute the surface tractions  (Input: QLOA)  follower forces 
22    call testq29(d,amqloa,numel,n,iflg1,iflg2)
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
      call rjac29(d,xl,t0,xs0,sx0,ts0,te0,detj0,isw)
      call shearfac29(d,detj0,xl,cappa,wcappa)
      if(ielas.eq.1)call dmate29 (d,cappa,dmat)
c      
c      
      nhom  = 1
c
      call pzero (pt,ngeas_max)
      call pzero (st,ngeas_max*ngeas_max)                    
      call pzero (eas,8)
      call pgauss (ngb,lint,sg,tg,wg)
      call gaus1D (ngz,pgz,wgz)
      if(ldir.eq.0) then 
        call triad29 (igeo,xl,d,rin)            
      else if(ldir.eq.1) then 
        do k = 1,4 
          call pdirec1(ix,basea,rin(1,1,k),k,2,1)  ! read from basea=m(mdir)
        end do 
      end if

      call updn29 (xl,ul,h3,twn,rwn,ddn,rin,ran,
     1             xan,gam,bml,w,ix,ndf,ndm,nen,ityp,ilin)
c
c
c.... loop over gauss points
      nhpl = 1
      do 2200 l = 1,lint
       xsi = sg(l)
       eta = tg(l)
       call shap29(xsi,eta,xl,t0,sx,shp,shp1,xsj)
       da = xsj*wg(l)
       rj = detj0/xsj

       call stra29 (xl,ul,shp,xsi,eta,eas,eps,sx,gam,rin,rwn,
     1              ddn,xu,x0,rd,r0,ndm,ndf,ilin)
       call pzero (eps,8)
       call strt29 (d,aqloa,numel,n,mqloa,propq,eps,epse,isw)

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
              call kine29(eps,x0,r0,zs,wcappa,epstr,dets)
              wz = wz*dets   
c
              call zerostr29(h1(nhpl),h2(nhpl),d,epstr,x0,
     +                       sigp,cmat)

c.....        compute material matrix and stress resultants
              call dmat29 (dmat,cmat,sig,sigp,wz,zs,wcappa,cappa)
              nhpl = nhpl + nh
2260        continue
2261     continue
       end if

c         
c....  loop over rows, residual
       ir0 = 0
       do 2210 inode = 1,4
        call bmat29 (inode,ir0,b,btd,dmat,st,pt,sig,qs,rwn,h3,
     1   w,twn,bml,d,shp1,shp,sx,t0,ran,xan,xu,rd,da,ngeas_max,ndf,ilin)
c....   loop over columns (symmetry noted), stiffness matrix
        jc0 = 0       
        do jnode = 1,inode
         call stif29 (inode,jnode,ir0,jc0,b,btd,st,sig,qs,shp,shp1,
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
       call rjac29(d,xl,t0,xs0,sx0,ts0,te0,detj0,isw)
       do l = 1,lint
        xsi = sg(l)
        eta = tg(l)
        call shap29(xsi,eta,xl,t0,sx,shp,shp1,xsj)
        da = xsj*wg(l)
        call qload29
     1       (shp1,xl,ul,t0,aqloa,numel,n,mqloa,propq,ndf,nst,p,s,da)
       enddo
      endif
c       
      return
      end
c
      subroutine testq29(d,q,numel,n,iflg1,iflg2)
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
      subroutine qload29(shp1,xl,ul,t0,q,numel,n,mqloa,propq,
     1                  ndf,nst,p,s,da)
c
c     q(n,1) = element n, load q_n (transverse) 
c     q(n,2) = element n, load q_x (global) 
c     q(n,3) = element n, load q_y (global)
c     q(n,4) = element n, load q_z (global)
c     q(n,5) = element n, iltyp, conservative = 0,  follower load  = 1
c     
      USE prlod 
      USE iofile 
      implicit double precision (a-h,o-z)
      dimension q(numel,10)
      dimension shp1(3,*),xl(3,*),ul(ndf,*),xu1(3,3),p(*),s(nst,*),
     1          t0(3,3),ql(3)
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
c
      subroutine rjac29(d,xl,t0,xs0,sx0,ts0,te0,detj0,isw)
c
c.... compute jacobian matrix at element center
c
      USE eldata 
      implicit double precision (a-h,o-z)
      dimension    a1(4),a2(4),hh(4),d(*),
     1             xs0(2,2),sx0(2,2), xl(3,*),te0(5,*),ts0(5,*),
     2             g0(3,3),g1(3,3),t0(3,3)            !,gs(3,3)
      data a1    / -0.25d0,  0.25d0, 0.25d0, -0.25d0 /
      data a2    / -0.25d0, -0.25d0, 0.25d0,  0.25d0 /
      data hh    /  0.25d0, -0.25d0, 0.25d0, -0.25d0 /

c.... convective base vectors
      call pzero (g0,9)
      call pzero (g1,9)
      do i = 1,3
       do K = 1,4
        g0(i,1) = g0(i,1) + a1(K) * xl(i,K)
        g0(i,2) = g0(i,2) + a2(K) * xl(i,K)
        g1(i,1) = g1(i,1) + hh(K) * xl(i,K)
        g1(i,2) = g1(i,2) + hh(K) * xl(i,K)
       end do
      end do
      call vecp (g0(1,1),g0(1,2),g0(1,3))

c
c.... t0 ------ Taylor
      do i = 1,3
       t0(i,1) = xl(i,3) - xl(i,1)
       t0(i,2) = xl(i,2) - xl(i,4)
      end do
      call norm (t0(1,1),t0(1,1),3)
      call norm (t0(1,2),t0(1,2),3)
      do i = 1,3
       v1 = t0(i,1)
       v2 = t0(i,2)
       t0(i,1) = v1 + v2
       t0(i,2) = v1 - v2
      end do
      call norm (t0(1,1),t0(1,1),3)
      call norm (t0(1,2),t0(1,2),3)
      call vecp (t0(1,1),t0(1,2),t0(1,3))
c.... t0-----  Hughes  lamina basis, Belytschko 541 
c      call norm (t0(1,3),g0(1,3),3)
c      call norm (gs(1,1),g0(1,1),3)
c      call norm (gs(1,2),g0(1,2),3)
c      do i = 1,3
c       gs(i,1) = gs(i,1) + gs(i,2)
c      end do
c      call vecp (t0(1,3),gs(1,1),gs(1,2))
c      do i = 1,3
c       t0(i,1) = gs(i,1) - gs(i,2)
c       t0(i,2) = gs(i,1) + gs(i,2)
c      end do
c      call norm (t0(1,1),t0(1,1),3)
c      call norm (t0(1,2),t0(1,2),3)
c.... t0-----  FG
c      call norm (t0(1,3),g0(1,3),3)
c      call norm (t0(1,1),g0(1,1),3)
c      call vecp (t0(1,3),t0(1,1),t0(1,2))

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
c
      subroutine shap29(xsi,eta,xl,t0,sx,shp,shp1,xsj)
      implicit double precision (a-h,o-z)
c
c.... shape function routine for shells
c
      dimension    s(4),t(4),xs(2,2),sx(2,2),shp1(3,4),shp(3,*),
     2		   xl(3,*),g(3,3),t0(3,3)
      data s/-0.5,0.5,0.5,-0.5/,t/-0.5,-0.5,0.5,0.5/
c.... shape functions 
      do i = 1,4
        shp (3,i) = (0.5+s(i)*xsi)*(0.5+t(i)*eta)
        shp1(3,i) = shp(3,i)
        shp1(1,i) = s(i)*(0.5+t(i)*eta)
        shp1(2,i) = t(i)*(0.5+s(i)*xsi)
      end do
c.... convective base vectors 
      call pzero (g,9)
      do i = 1,3
        do k = 1,4
          g(i,1) = g(i,1) + shp1(1,k) * xl(i,k)
          g(i,2) = g(i,2) + shp1(2,k) * xl(i,k)
        end do
      end do
      call vecp (g(1,1),g(1,2),g(1,3))
c.... construct jacobian and its inverse
      do i = 1,2
        do j = 1,2
         xs(i,j) = dot(g(1,i),t0(1,j),3)
        end do
      end do
      detj = xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
      sx(1,1) = xs(2,2)/detj
      sx(2,2) = xs(1,1)/detj
      sx(1,2) =-xs(1,2)/detj
      sx(2,1) =-xs(2,1)/detj
c.... form local derivatives
      do i = 1,4
       shp(1,i) = sx(1,1)*shp1(1,i) + sx(1,2)*shp1(2,i)
       shp(2,i) = sx(2,1)*shp1(1,i) + sx(2,2)*shp1(2,i)
      end do
c.... area element
      xsj = dsqrt(dot(g(1,3),g(1,3),3))
      return
      end
c      
      subroutine easstr29 (xsi,eta,te0,rj,neas_max,neas,nm,nb,ns,
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
      subroutine easmat29 (gxy,gtd,pt,st,sig,dmat,b,da,neas,
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
      subroutine easend29(st,pt,s,p,skcr,prc,nst,ndf,ngeas_max,neas)
c-----------------------------------------------------------------------
c.... static condensation and store matrices for back substitution
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension st(ngeas_max,*),pt(*),s(nst,*),p(*),skcr(neas,*),prc(*)

c.... static condensation
      if(neas.gt.0) call conden(st,pt,4*ndf,neas,ngeas_max,.false.)

c.... copy stiffness matrix st to s and pt to p 
      do i = 1,4*ndf
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
      subroutine bmat29 (inode,ir0,b,btd,dmat,st,pt,sig,qs,rwn,dwn,w,
     1     twn,bml,d,shp1,shp,sx,t0,ran,xan,xu,rd,da,ngeas_max,ndf,ilin)
c-----------------------------------------------------------------------
c.... b-matrix, residual vector, external load vector  
c-----------------------------------------------------------------------
      USE prlod 
      USE eldata 
      implicit double precision (a-h,o-z)
      dimension d(*),shp1(3,*),shp(3,*),pt(*),st(ngeas_max,*),b(8,6,*),
     1          btd(8,*),dmat(8,8),xu(3,*),sx(2,2),rd(3,*),t0(3,*),
     2          rwn(3,3,*),fbb(3,2),ql(3),
     3          sig(*),qs(*),wmh(3,3),xmh(3),dwn(3,*),ran(3,*),
     4          bml(3,2,*),twn(3,3,*),xan(3,*),w(3,3,*),wh(3,3) 
      dimension im(4),il(4)
      data im /-2, 2, 4,-4/, il /-1,-3, 3, 1/
c
c
      mm = iabs(im(inode))
      ll = iabs(il(inode))
      fm = isign(1,im(inode))*shp1(1,inode)
      fl = isign(1,il(inode))*shp1(2,inode)
      call mttmul (w(1,1,inode),xu,3,3,2, fbb)
c
c.... b-matrix: membrane (1,2,3), bending (4,5,6), shear (7,8) 
      do i = 1,3
        b1 = ran(i,mm)*shp1(1,inode)
        b2 = ran(i,ll)*shp1(2,inode)
        b(1,i,inode) = xu(i,1)*shp(1,inode)
        b(2,i,inode) = xu(i,2)*shp(2,inode)
        b(3,i,inode) = xu(i,1)*shp(2,inode) + xu(i,2)*shp(1,inode)
        b(4,i,inode) = rd(i,1)*shp(1,inode)
        b(5,i,inode) = rd(i,2)*shp(2,inode)
        b(6,i,inode) = rd(i,1)*shp(2,inode) + rd(i,2)*shp(1,inode)
        b(7,i,inode) = sx(1,1)*b1 + sx(1,2)*b2
        b(8,i,inode) = sx(2,1)*b1 + sx(2,2)*b2
      end do  

      do i = 1,ndf-3
        i3 = i+3
        b3 = fm*bml(i,1,inode)
        b4 = fl*bml(i,2,inode)
        b(1,i3,inode) = 0.0d0
        b(2,i3,inode) = 0.0d0
        b(3,i3,inode) = 0.0d0
        b(4,i3,inode) = shp(1,inode)*fbb(i,1)
        b(5,i3,inode) = shp(2,inode)*fbb(i,2)
        b(6,i3,inode) = shp(1,inode)*fbb(i,2) + shp(2,inode)*fbb(i,1)
        b(7,i3,inode) = sx(1,1)*b3 + sx(1,2)*b4
        b(8,i3,inode) = sx(2,1)*b3 + sx(2,2)*b4
      end do

c.... residual vector and  B_T * D
      do i = 1,ndf
        ir = ir0 + i 
        pt(ir) = pt(ir) - dot(sig,b(1,i,inode),8)*da
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
        ql(1) = ql(1) + qx
        ql(2) = ql(2) + qy
        ql(3) = ql(3) + qz
c..... add loads to load vector
        do i = 1,3
         pt(ir0+i) = pt(ir0+i) + prop*ql(i) * shp(3,inode)*da
        end do
      end if
c
c....  transform q to convective coordinate system
       qs(1) = sx(1,1)*sig(7) + sx(2,1)*sig(8)
       qs(2) = sx(1,2)*sig(7) + sx(2,2)*sig(8)

      if(ilin.eq.0.or.ilin.eq.1) return
c
c....  geometric stiffness  (diagonal terms)
c      if(dot(dwn(1,inode),dwn(1,inode),3).gt.1.d-15)then 
       do i = 1,3
         xmh(i) = ( sig(4)*xu(i,1)*shp(1,inode) 
     1            + sig(5)*xu(i,2)*shp(2,inode)
     2            + sig(6)*xu(i,2)*shp(1,inode) 
     3            + sig(6)*xu(i,1)*shp(2,inode)
     4            + qs(1)*fm*xan(i,mm)+qs(2)*fl*xan(i,ll))*da
       end do
       call faca29 (rwn(1,3,inode),xmh,dwn(1,inode), wmh)
       call mttmul (twn(1,1,inode),wmh,3,3,3, wh)
       call matmulf (wh,twn(1,1,inode),3,3,3, wmh)
       do i = 1,ndf-3
         ir2 = ir0+3+i
         do j = 1,i
           ic2 = ir0+3+j
           st(ir2,ic2) = st(ir2,ic2) + wmh(i,j)
         end do
       end do    
c      end if
c
      return
      end
c
      subroutine boun29 (inode,ix,id,iflag,ndf)
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
      subroutine shearfac29(d,detj0,xl,cappa,wcappa)
c-----------------------------------------------------------------------
c.... shear correction factor
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  d(*),xl(3,*),yl(3)

      cappa = d(19)            

      if(cappa.lt.0.d0)then

c....     AE... element area
c          ae = 4.d0*detj0 
       
c....     AE... or square of longest element side
          call pzero(yl,3) 
          ae = 0.d0
          do i = 1,4
            k = i+1
            if(i.eq.4) k=1
            do j = 1,3
              yl(j) = xl(j,k)-xl(j,i)
            end do 
            sl2 = dot(yl,yl,3)       
            ae = max(ae,sl2)
          end do  

          hs  = d(13)
          xnu = 0  
          cappa  = dabs(d(19))/(1.d0+ae/(2.d0*hs*hs*(1.d0+xnu)))
      end if 
        wcappa = dsqrt(cappa) 
c      
      return
      end  
c
      subroutine dmate29 (d,cappa,dmat)
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
      subroutine kine29(eps,x0,r0,zs,wcappa,epstr,dets)
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
      subroutine dmat29 (dmat,cmat,sig,sigp,wz,zs,wcappa,cappa)
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
      subroutine faca29 (an,x,wn,aix)
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
      subroutine transsig29(t0,xu,rd,rin,shp,sig,siga)
c-----------------------------------------------------------------------
c.... transform the stresses to basis from igeo
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension t0(3,3),glu(3,3),glo(3,3),gko(3,3),rl(3,3),xu(3,3),
     +          rd(3,3),ts(3,3),sig(*),siga(9),sigs(9),t(3,3),
     +          rin(3,3,4),shp(3,4) 
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
        do k = 1,4
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
      subroutine plsn29 (i)
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
      subroutine plot29 (d,ix,dt,st,shp,sig,sigp,h2,eps,
     +                   numnp,maxlay,klay,da)
c-----------------------------------------------------------------------
c.... Plot stresses  at specified layer for general shell element
c     klay = 0:           stress resultants
c          = 1 to maxlay: stresses and history values at layer klay
c-----------------------------------------------------------------------
      implicit  double precision (a-h,o-z)
      dimension ix(*),dt(numnp),st(numnp,*),shp(3,*),sig(*),sigp(*),
     + h2(*),d(*),eps(*)
      if (klay.gt.0 .and. klay.le.maxlay) then
c....   print stresses and history values in layer klay
        do i = 1,5
          sig(i) = sigp(i)
        end do
        sig(6) = sigp(6)
        sig(7) = h2(7)
        sig(8) = h2(8)
      end if
c

c.... membrane energy   
c      hs     = d(13)
c      emod   = d(40)  
c      xnu    = d(41)
c      Dm     = emod*hs/(1.d0-xnu*xnu)
c      wm = 0.5d0*Dm*(eps(1)*eps(1)+ eps(2)*eps(2)+2.d0*xnu*eps(1)*eps(2)
c     1     +(1.d0-xnu)*eps(3)*eps(3)/2.d0) 

c...  store values
      do 20 j = 1,4
        xsji = shp(3,j)*da
        ii = abs(ix(j))
        if(ii.eq.0) go to 20
        dt(ii) = dt(ii) + xsji

        do 30 i = 1,9
30        st(ii,i)  = st(ii,i) + sig(i)*xsji
c          st(ii,10) = st(ii,10) + wm*xsji
20    continue
      return
      end
c
      subroutine pout29 (xl,shp,sig,sigp,epn,ndm,maxlay)
c-----------------------------------------------------------------------
c.... print stresses at Gauss points of each layer
c     klay = 0:           stress resultants
c          = 1 to maxlay: stresses and history values at layer klay
c
c     1. PK-Stresses/stress resultants with respect to deformed ref.surf.
c-----------------------------------------------------------------------
      USE bdata       
      USE plslay 
      USE iofile 
      USE eldata 
      implicit double precision (a-h,o-z)
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
1000      format(a1,20a4,//,2x,'E L E M E N T  29 ',
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
1001      format(a1,20a4,//,2x,'E L E M E N T  29  Cauchy stresses',
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


cww1001      format(a1,20a4,//,2x,'E L E M E N T  29   S T R E S S E S',
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
      subroutine stif29 (inode,jnode,ir0,jc0,b,btd,
     1                   st,sig,qs,shp,shp1,w,ngeas_max,ndf,da,ilin,isw)
c-----------------------------------------------------------------------
c.... stiffness matrix
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  st(ngeas_max,*),sig(*),shp(3,*),shp1(3,*),w(3,3,*),
     1           id1(4,4),id2(4,4),qs(*),b(8,6,*),btd(8,*)
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
      if1  = id1(inode,jnode)
      if2  = id2(inode,jnode)
        sgd1 =  sig(1)*shp(1,inode)*shp(1,jnode) 
     1        + sig(2)*shp(2,inode)*shp(2,jnode)
     2        + sig(3)*shp(1,inode)*shp(2,jnode) 
     3        + sig(3)*shp(2,inode)*shp(1,jnode)
        sgd2  = sig(4)*shp(1,inode)*shp(1,jnode) 
     1        + sig(5)*shp(2,inode)*shp(2,jnode)
     2        + sig(6)*shp(1,inode)*shp(2,jnode) 
     3        + sig(6)*shp(2,inode)*shp(1,jnode)
        sgd3 = 0.5d0*(qs(1)*shp1(1,inode)*if1+qs(2)*shp1(2,inode)*if2)
        sgd4 = 0.5d0*(qs(1)*shp1(1,jnode)*if1+qs(2)*shp1(2,jnode)*if2)
c
      guu =  sgd1       * da
      guw = (sgd2+sgd3) * da
      gwu = (sgd2+sgd4) * da

      do i = 1,3
       ir1 =   ir0+i
       ic1 =   jc0+i
       st(ir1,ic1) = st(ir1,ic1) + guu
       do j = 1,ndf-3
         ir2 = ir0+3+j
         ic2 = jc0+3+j
         st(ir1,ic2) = st(ir1,ic2) + w(i,j,jnode) * guw
         st(ir2,ic1) = st(ir2,ic1) + w(i,j,inode) * gwu
       end do
      end do
c
      return
      end
c
      subroutine stra29 (xl,ul,shp,xsi,eta,eas,eps,sx,gam,rin,rwn,
     1                   ddn,xu,x0,rd,r0,ndm,ndf,ilin)
c-----------------------------------------------------------------------
c     Green-Lagrangean strains of the reference surface
c     membrane(1,2,3),     curvatures(4,5,6),    shear(7,8) 
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension xl(ndm,*),ul(ndf,*),shp(3,*),ddn(3,*),
     1          eps(*),rin(3,3,*),rwn(3,3,*),uu(3,3),dd(3,3),eas(*),
     2          x0(3,3),xu(3,3),r0(3,3),rd(3,3),gam(*),sx(2,2)
c
c.... current and initial geometry
c
      do i = 1,3
       do j = 1,3
        rd(i,j)  = 0.d0
        r0(i,j)  = 0.d0
        xu(i,j)  = 0.d0
        x0(i,j)  = 0.d0
         do k = 1,4
          rd(i,j)  = rd(i,j)  + rwn(i,3,k) * shp(j,k)
          r0(i,j)  = r0(i,j)  + rin(i,3,k) * shp(j,k)
          xu(i,j)  = xu(i,j)  + (xl(i,k)+ul(i,k))*shp(j,k)
          x0(i,j)  = x0(i,j)  +  xl(i,k)         *shp(j,k)
         end do
       end do 
        x0(i,3) = r0(i,3)
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
         dd(i,j) = 0.d0
         xu(i,j) = x0(i,j)  
         rd(i,j) = r0(i,j)  
          do k  = 1,4
           uu(i,j) = uu(i,j) + ul(i,k)*shp(j,k)
           dd(i,j) = dd(i,j) + ddn(i,k)*shp(j,k)
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
      eps(7) = sx(1,1)*gam1 + sx(1,2)*gam2 + eas(7)
      eps(8) = sx(2,1)*gam1 + sx(2,2)*gam2 + eas(8)
c   
      return
      end
c
      subroutine strt29 (d,q,numel,n,mqloa,propq,eps,epse,isw)
c-----------------------------------------------------------------------
c     - temperature strains
c     membrane(1,2,3),     curvatures(4,5,6),    shear(7,8) 
c-----------------------------------------------------------------------
      USE prlod 
      implicit double precision (a-h,o-z)
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
      subroutine trasig29 (sigp,t0,shp,rin)               
c-----------------------------------------------------------------------
c.... transform stresses         
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  sigp(*),t0(3,3),ts(3,3),t(3,3),
     1           tr(5,5),sigt(6),rin(3,3,4),shp(3,4)
c
      do i = 1,3
       do j = 1,3
        ts(i,j) = 0.d0 
        do k = 1,4
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
      subroutine triad29 (igeo,xl,d,rin)
c-----------------------------------------------------------------------
c     initial cartesian system for specified geometry:
c     igeo =  1      ... arbitrary surface 
c             others ... see macro BASE
c
c     rin(3,3,4)     initial local triad at nodes
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  xl(3,*),g1(3,4),g2(3,4),rin(3,3,*),e3(3),d(*),add(3)
      data e3 / 0.d0,0.d0,1.d0 /

          add(1)=d(2)
          add(2)=d(3)
          add(3)=d(4)

c
      call pzero(rin,9 * 4)
c
      do 20 i = 1,3
       g1(i,1) = xl(i,2)-xl(i,1)
       g1(i,2) = g1(i,1)
       g1(i,3) = xl(i,3)-xl(i,4)
       g1(i,4) = g1(i,3)
       g2(i,1) = xl(i,4)-xl(i,1)
       g2(i,2) = xl(i,3)-xl(i,2)
       g2(i,3) = g2(i,2)
       g2(i,4) = g2(i,1)
20    continue

      do k = 1,4

       if(igeo.eq.1.or.igeo.eq.2) then  
        call vecp (g1(1,k),g2(1,k),rin(1,3,k))
        call norm (rin(1,3,k),rin(1,3,k),3)
        call norm (rin(1,1,k),g1(1,k),3)
        call vecp (rin(1,3,k),rin(1,1,k),rin(1,2,k))
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
      subroutine updn29 (xl,ul,dwn,twn,rwn,ddn,rin,ran,xan,
     1                   gam,bml,w,ix,ndf,ndm,nen,ityp,ilin)
c-----------------------------------------------------------------------
c     update of nodal rotations
c     rin    initial local cartesian system   at element nodes
c     rwn    current nodal basis  R = R(om)       - " -
c     ran    current       director           at midside nodes
c-----------------------------------------------------------------------
      USE mdata       
      implicit double precision (a-h,o-z)
      dimension  xl(ndm,*),ul(ndf,*),ix(*),
     1           dwn(3,*),twn(3,3,*),rwn(3,3,*),ddn(3,*),rin(3,3,*),
     2           w(3,3,*),rh(3,3),rh1(3,3),rh2(3),gam(*),bml(3,2,*),
     3           ran(3,*),xan(3,*),ran0(3,4),xan0(3,4),
     4           uan(3,4),dan(3,4),im(4),il(4)
      data im /-2, 2, 4,-4/, il /-1,-3, 3, 1/
c
      do 100  k = 1,4
c.... initial nodal cartesian system
       call matcop (rin(1,1,k),3,3, rh1)
       call matcop (ul(4,2*nen+k),3,1, rh2)          

c.... transformation of rotations 5dofs --> 6dofs
         if(ndf.eq.5) iflag = 1                          ! iflag = 1 : glatte Schale
         if(ndf.eq.6) call boun29(k,ix,psid,iflag,ndf)  ! iflag = 0 : Kante  
         if(ityp.eq.6) iflag = 0     
         if (iflag.eq.1) then
          if(ilin.eq.2) then
           call updr29 (dwn(1,k),rh,ilin,1)
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

        call updr29 (dwn(1,k),rh,ilin,1)
        call matmulf (rh,rin(1,1,k),3,3,3, rwn(1,1,k))

       end if 

c....  T = W^T * H * T_3
       call updr29 (dwn(1,k),twn(1,1,k),ilin,3)
       if(ilin.eq.2)then                 !  finite rotations
        call skew (w(1,1,k),rwn(1,3,k))                      
        if (iflag.eq.1) then
         call matmulf (twn(1,1,k),rwn(1,1,k),3,3,3, rh)
         call matcop (rh,3,3, twn(1,1,k))
        end if

       else                              ! linear and moderate rotations

        call skew (w(1,1,k),rin(1,3,k))                      
        if (iflag.eq.1) call matcop (rin(1,1,k),3,3, twn(1,1,k))
       endif

       call mttmul (w(1,1,k),twn(1,1,k),3,3,3, rh)
       call matcop (rh,3,3, w(1,1,k))
         
100   continue

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
      subroutine updr29 (dwn,rw,ilin,ifl)
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

      tol29 = 1.d-16
      theta = dsqrt(dot(dwn,dwn,3))
      thet2 = theta*theta

      call skew (om,dwn)
      call matmulf (om,om, 3,3,3, om2)

      if (ifl .eq. 1) then                       ! R(om)
        if (dabs(theta) .lt. tol29) then           ! theta -> 0
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
        if (dabs( thet2*(dcos(theta)-1.d0) ) .lt. tol29) then
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
        if (dabs(theta) .lt. tol29) then           ! theta -> 0
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
      subroutine zerostr29(h1,h2,d,epstr,glu,sigp,cmat)
c-----------------------------------------------------------------------
c     zero thickness normal stress iteration
c     im = 1:         Pegasus method       2: Newton´s method
c-----------------------------------------------------------------------
      USE iofile 
      implicit double precision (a-h,o-z)
      dimension h1(*),h2(*),d(*),ign(6),epstr(*),sigp(*),cmat(5,*)
     2         ,e3d(6),sig3d(7),Cmat3d(6,6),glu(3,3)
      data ign/1,2,4,5,6,3/
c
c.....setup starting values
c
      tol29 = 1.d-7
      iter  = 0
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
      call matelib29(h1,h2,d,e3d,glu,sig3d,cmat3d)
      f1   = sig3d(3)
      
      e3d(3) = x2
      call matelib29(h1,h2,d,e3d,glu,sig3d,cmat3d)
      f2   = sig3d(3)
cfg         if(f1*f2.gt.0)write(*,*)'Starting values Pegasus not ok'

cfgc.....   find admissible range
cfgc
cfg      de33 = 1.d-1
cfg
cfg      call matelib29(h1,h2,d,e3d,glu,sig3d,cmat3d)
cfg      do n = 1,20 
cfg       a  = e3d(3)
cfg       fa = sig3d(3)
cfg       e3d(3) = e3d(3) + de33
cfg       call matelib29(h1,h2,d,e3d,glu,sig3d,cmat3d)
cfg       fb = sig3d(3)
cfg       if(((fa*fb).gt.0.d0).and.(dabs(fb).gt.dabs(fa)))go to 1
cfg       if(fa*fb.lt.0.d0)go to 2
cfg      end do
cfg
cfg1     e3d(3) = h1(8)
cfg      call matelib29(h1,h2,d,e3d,glu,sig3d,cmat3d)
cfg      do n = 1,20 
cfg       b  = e3d(3)
cfg       fa = sig3d(3)
cfg       e3d(3) = e3d(3) - de33
cfg       call matelib29(h1,h2,d,e3d,glu,sig3d,cmat3d)
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
cfg      call matelib29(h1,h2,d,e3d,glu,sig3d,cmat3d)
cfg      f1   = sig3d(3)
cfg      
cfg      x2   = b
cfg      e3d(3) = x2
cfg      call matelib29(h1,h2,d,e3d,glu,sig3d,cmat3d)
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
        call matelib29(h1,h2,d,e3d,glu,sig3d,cmat3d)
        Czz = 1.d0/Cmat3d(3,3)       
        f3 = sig3d(3)
        sigv = dsqrt(sig3d(1)*sig3d(1)+sig3d(2)*sig3d(2))
        if(teste.lt.1.e-20)sigv=1.d0
        if (dabs(sig3d(3))/sigv.le.tol29) go to 100

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
      call matelib29(h1,h2,d,e3d,glu,sig3d,cmat3d)

      Czz = 1.d0/Cmat3d(3,3)
      e3d(3) = e3d(3) - sig3d(3)*Czz
      sigv = dsqrt(sig3d(1)*sig3d(1)+sig3d(2)*sig3d(2))
      if(teste.lt.1.e-20)sigv=1.d0
      if (dabs(sig3d(3))/sigv.le.tol29) go to 100
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
      subroutine matelib29(h1,h2,d,e,glu,sig,cmat)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension h1(*),h2(*),d(*),e(*),glu(3,3),sig(*),Cmat(6,*)
c
      imat = d(15)
c
      if (imat.eq.1) then
          call mate0129(h1,h2,d,e,glu,sig,cmat)
      elseif(imat.eq.2) then
          call mate0229(h1,h2,d,e,glu,sig,cmat)
      elseif(imat.eq.3) then
          call mate0329(h1,h2,d,e,glu,sig,cmat)
      end if
c
      return
      end
c
      subroutine mate0129(h1,h2,d,e,glu,sig,cmat)
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
      USE eldata 
      USE iofile 
      implicit double precision (a-h,o-z)
      dimension d(*),h1(*),h2(*),sig(*),cmat(6,*),e(*),ep(6),
     1          ee(6),sd(6),S(6),glu(3,3),ec(6),gko(3,3),glo(3,3),
     2          Te(6,6),to(3,3),t(3,3),Tet(6,6),xhelp(6,6),cort(6,6)
c    
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
      call mat0129(cort,xE,xn)
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
      if(f.gt.1.d-10*y) call plas0129(d,sig,sd,f,g_tr,y,a,ep,cort)

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
      subroutine plas0129(d,sig,sd,f,g_tr,y,a,ep,cmat)
c-----------------------------------------------------------------------
      USE iofile 
      implicit double precision (a-h,o-z)
      dimension d(*),sig(*),sd(*),ep(*),cmat(6,6)
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
         write( * ,*) 'no convergence in plas0129 '
         write(iow,*) 'no convergence in plas0129'
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
      subroutine mat0129(cmat,e1,xnu)
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
      subroutine mate0229(h1,h2,d,e,glu,S,cmat)
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
      USE eldata 
      USE iofile 
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
         call plas0229(d,T_tr,Td_tr,g_tr,f,Ce,ep,a)
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
      subroutine plas0229(d,T_tr,Td_tr,g_tr,f,Ce,ep,a)
c-----------------------------------------------------------------------
      USE iofile 
      implicit double precision (a-h,o-z)
      dimension d(*),T_tr(3),Td_tr(3),Ce(3,3),ep(3)
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
         write( * ,*) 'no convergence in plas0229'
         write( * ,*) 'f=',f
         write(iow,*) 'no convergence in plas0229'
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
      subroutine mate0329(h1,h2,d,e,glu,S,cmat)
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
      USE eldata 
      USE iofile 
      implicit double precision (a-h,o-z)
      double precision L1(3,3),L2(3,3)
      dimension d(*),S(*),cmat(6,6),glu(3,3),e(6)
      dimension G(3,3),C(3,3),C_n(3,3),G_n(3,3),
     +   xmue(3),alpha(3),xlam(3),xN(3,3),fv1(3),fv2(3),sig(3),
     +   T1(6,3),T2(6,3),T1L1(6,3),T2L2(6,3),T1L1T1(6,6),T2L2T2(6,6)
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
      subroutine matread29(d)
c-----------------------------------------------------------------------
c.... input of material data
c-----------------------------------------------------------------------
      USE iofile 
      implicit double precision (a-h,o-z)
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
