      subroutine elmt02(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c-----------------------------------------------------------------------
c.... finite rotation 2D-beam element, Theory: ww                    

c.... lumped and consistent/higher order (Hughes) mass matrices                            
c.... lumped damping matrix                                          

c.... plot,forc ueber residuum r = k*v-p                             
c.... plot,forc,1 = N, ..,2 = Q, ..,3 = M                            

c.... sloa: follower loads                                           
c.... qloa: follower loads                                           

c.... kappa = 5/6 
c.... Shear correction factor factor FE - Bischoff/Bletzinger
c     kappaGA=kappaGA/(1+kappa*l^2/12*GA/EI), 

c.... ityp = 1(0) fully nonlinear = default  Green                   
c....      = 2    fully nonlinear =          Reissner                
c....      = 3    moderate rotations         Green                   
c....      = 4    moderate rotations-simplified Green                
c....      = 5    linear                                             
c-----------------------------------------------------------------------
c     Input: E G A I rho eta alpha_T h ityp                            
c            pl pr nl nr ifol T_u T_o                                            
c                                                                    
c-----------------------------------------------------------------------
c.... d(1)  = E                                                      
c.... d(2)  = G                                                      
c.... d(3)  = A                                                      
c.... d(4)  = I                                                      
c.... d(5)  = rho                                                    
c.... d(6)  = eta                                                    
c.... d(7)  = alpha_T                                                
c.... d(8)  = h default=1                                            
c.... d(9)  = ityp                                                   

c.... d(10) = pl                                                     
c.... d(11) = pr                                                     
c.... d(12) = nl                                                     
c.... d(13) = nr                                                     

c.... d(15)= EA                                                      
c.... d(16)= EI                                                      
c.... d(17)= GA                                                      

c.... d(18)= rhoA                                                    
c.... d(19)= rhoI                                                    
c.... d(20)= etaA                                                    
c.... d(21)= etaI                                                    
c....                                                                
c-----------------------------------------------------------------------
c     Input for SLOA (isw=7)                                         
c     q01 = p                                                        
c-----------------------------------------------------------------------
c     Input for QLOA (isw=22)                                        
c     q01 = pl                                                       
c     q02 = pr                                                       
c     q03 = nl                                                       
c     q04 = nr                                                       
c     q05 = ifol: 0=no follower load 1= follower load  only for q    
c     q06 = T_u: Temperature at Bottom                               
c     q07 = T_o: Temperature at Top                                  
c     in case of follower load: q01=q02!!                            
c
c
c-----------------------------------------------------------------------
c.....(c) w.wagner  Version x-z                                      
c-----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE eldata
      USE fornam
      USE iofile
      USE pdata7
      USE pdata10
      USE plodfu
      USE prlod
      USE qload
      USE strnam
      USE tdata
      implicit double precision(a-h,o-z)
      dimension ix(*),ixl(5),xl(ndm,*),xu(ndm,2),tl(*),xll(3,4),
     +       d(*),gp(2),s(nst,nst),shp(2,2),dmat(3,3),btd(3,3),
     +       bm(3,3,2),sig(3),ul(ndf,*),vl(3,2),p(nst),gr(5),ps(6),
     +       ql(7),pt(6),st(6,6),sigt(3)
      dimension h1(*),h2(*),h3(*)

      data ixl/1,2,3,4,1/,eps/1.0e-5/

c.... go to correct array processor
      go to(1,2,3,4,5,3,7,2,2,2,3,12,4,2,2,2,2,2,2,2,2,22) isw
      return
c.... input material properties
1     if(ior.lt.0) write(*,1000)
1000  format(' Input: E G A I rho eta alpha_T h ityp'/3x,'>',$)
      call dinput(d,9)
    
c.... default G=E/2
      g = d(2)
      if(g.eq.0.d0) d(2) = 0.5d0*d(1)  
c.... default thickness, if no temperature loading
      thick = d(8)
      if(thick.eq.0.d0) d(8) = 1.d0  
c.... default ityp
      ityp = d(9)
      if(ityp.eq.0) ityp=1 
      d(9) = ityp

      if(ior.lt.0) write(*,1001)
1001  format(' Input: ql qr nl nr'/3x,'>',$)
      call dinput(d(10),4)

      if(ior.lt.0) write(*  ,1002) (d(i),i=1,13)
                   write(iow,1002) (d(i),i=1,13)
1002  format(5x,'Materialdata 2D-Finite Rotation Beam Element:',/,
     * 5x,'E.................',g12.4,/,
     * 5x,'G.................',g12.4,/,
     * 5x,'A.................',g12.4,/,
     * 5x,'I.................',g12.4,/,
     * 5x,'rho(gamma/g) .....',g12.4,/,
     * 5x,'eta(viscosity)....',g12.4,/,
     * 5x,'alpha_T ..........',g12.4,/,
     * 5x,'thickness h ......',g12.4,/,
     * 5x,'typ:1=FR_Gr,2=FR-Rei,3=MR1,4=MR2,5=lin..',f4.0,/,
     * 5x,'load q node 1 ....',g12.4,/,
     * 5x,'load q node 2 ....',g12.4,/,
     * 5x,'load n node 1 ....',g12.4,/,
     * 5x,'load n node 2 ....',g12.4)

      ea    = d(1)*d(3)                    
      ei    = d(1)*d(4)                  
      ga    = d(2)*d(3)       
      d(15) = ea           ! EA
      d(16) = ei           ! EI
      d(17) = ga           ! GA
      d(18) = d(5)*d(3)    ! rhoA
      d(19) = d(5)*d(4)    ! rhoI
      d(20) = d(6)*d(3)    ! etaA
      d(21) = d(6)*d(4)    ! etaI
      ipla  = 2
c.... description of stresses  
      forsus( 1) =  '  N-FORCE N_x  '
      forsus( 2) =  '  Q-FORCE Q_z  '
      forsus( 3) =  '  MOMENT  M_y  '
      do i = 4,11         
         forsus(i) =  ' '
      end do
2     return
c.... stiffness  matrix
3     ityp = d(9)

c.... Load terms
      call qload02(ql,d,aqloa,numel,n,mqloa,propq,prop,isw)

c.... length,angle,radius
      sn = xl(2,2)-xl(2,1)
      cs = xl(1,2)-xl(1,1)
      sl = dsqrt(cs*cs + sn*sn)
      sn = sn/sl
      cs = cs/sl
      dv = sl

c.... shape functions
      shp(1,1) = -1.0d0/sl
      shp(1,2) =  1.0d0/sl
      shp(2,1) =  0.5d0
      shp(2,2) =  0.5d0

c.... local  displacements
      do k = 1,2
        vl(1,k) = cs*ul(1,k) + sn*ul(2,k)
        vl(2,k) =-sn*ul(1,k) + cs*ul(2,k)
        vl(3,k) = ul(3,k)
      end do

c.... elasticity matrix
      ea = d(15)
      ei = d(16)
      ga = d(17)
      cappa = 5.d0/6.d0
      gac = cappa*ga/(1.d0+cappa*sl*sl/12.d0*ga/ei)

      call pzero(dmat,9)
      dmat(1,1) =  ea
      dmat(2,2) =  ei
      dmat(3,3) =  gac

c.... strains,stresses
      call stre02(sig,shp,vl,dmat,gr,epoten,ityp)

c.... Loadvector local
      p(1) = p(1) - ( 2.*ql(3) +    ql(4)) * sl  / 6.
      p(4) = p(4) - (    ql(3) + 2.*ql(4)) * sl  / 6.

      p(2) = p(2) - ( 2.*ql(1) +    ql(2)) * sl  / 6.
      p(5) = p(5) - (    ql(1) + 2.*ql(2)) * sl  / 6.

c.... stiffness matrix and residual
      i1=0
      do 31 ii=1,2
c....   B- matrix I
        call bmat02(bm,shp,vl,gr,ii,ityp)
c....   residual G = P - Bt*S and matrix Bt*D
        do 32 i = 1,3
          do  33 k = 1,3
            btd(i,k) = 0.0
            p(i1+i) = p(i1+i) - bm(k,i,ii)*sig(k)*dv
            do  34 j = 1,3
34          btd(i,k) = btd(i,k)+bm(j,i,ii)*dmat(j,k)
33        continue
32      continue
        if(isw.eq.6) go to 305
c....   tangent stiffness matrix
        j1 = 0
        do 35 jj = 1,ii
          do 36  i = 1,3
            do 37  j = 1,3
              do 38     k = 1,3
38            s(i1+i,j1+j) = s(i1+i,j1+j) + btd(i,k)*bm(k,j,jj)*dv
37          continue
36        continue
c....     K sigma
          call ksig02(sig,shp,gr,s,nst,dv,ii,jj,i1,j1,ityp)
          j1 = j1 + ndf
35      continue
305     i1 = i1 + ndf
31    continue
      if(isw.eq.6) go to 300

c.... upper part of stiffness matrix
      do 39 i = 1,3
        do 39 j = 1,3
39       s(j,i+ndf) = s(i+ndf,j)
      call trans02(s,cs,sn,nst,ndf,1)
300   continue
      call trans02(p,cs,sn,nst,ndf,2)
      return

c.....plot/print stresses at nodes
4     ityp = d(9)

c.... Load terms
      call qload02(ql,d,aqloa,numel,n,mqloa,propq,prop,isw)

c.... length,angle,radius
      sn = xl(2,2)-xl(2,1)
      cs = xl(1,2)-xl(1,1)
      sl = dsqrt(cs*cs + sn*sn)
      sn = sn/sl
      cs = cs/sl
      dv = sl

c.... shape functions
      shp(1,1) = -1.0d0/sl
      shp(1,2) =  1.0d0/sl
      shp(2,1) =  0.5d0
      shp(2,2) =  0.5d0

c.... local  displacements
      do k = 1,2
        vl(1,k) = cs*ul(1,k) + sn*ul(2,k)
        vl(2,k) =-sn*ul(1,k) + cs*ul(2,k)
        vl(3,k) = ul(3,k)
      end do

      ea = d(15)
      ei = d(16)
      ga = d(17)
      cappa = 5.d0/6.d0
      gac = cappa*ga/(1.d0+cappa*sl*sl/12.d0*ga/ei)

      call pzero(dmat,9)
      dmat(1,1) =  ea
      dmat(2,2) =  ei
      dmat(3,3) =  gac

      ifol=ql(5)

c.... load case: conservative loads 
      if(ifol.eq.0) then
        p(1) =  - ( 2.*ql(3) +    ql(4)) * sl  / 6.
        p(4) =  - (    ql(3) + 2.*ql(4)) * sl  / 6.

        p(2) = - ( 2.*ql(1) +    ql(2)) * sl  / 6.
        p(5) = - (    ql(1) + 2.*ql(2)) * sl  / 6.
      end if

      if(ifol.eq.1) then

c.... load case: non conservative loads 

c....   Normal vector 
        dn1 =      - (vl(2,2)-vl(2,1)) / sl
        dn2 = 1.d0 + (vl(1,2)-vl(1,1)) / sl

c....   residual constant load qf 
        qf = ql(1)
 
        p(1) =  qf * 0.5d0 * sl * dn1 
        p(2) =  qf * 0.5d0 * sl * dn2 
        p(4) =  qf * 0.5d0 * sl * dn1
        p(5) =  qf * 0.5d0 * sl * dn2

c.... add conservative loads n1,n2     
        p(1) =  p(1) - ( 2.*ql(3) +    ql(4)) * sl  / 6.
        p(4) =  p(4) - (    ql(3) + 2.*ql(4)) * sl  / 6.

      end if

c.... strains,stresses
      call stre02(sig,shp,vl,dmat,gr,epoten,ityp)

c.... load case: temperature
      at = d(7)
      he = d(8)
      ea = d(15)
      ei = d(16)
      tu = ql(6)
      to = ql(7)

      if((tu*tu+to*to).eq.0.d0) goto 42

      sig(1) = sig(1) - ea * (to+tu)*at*0.5d0 
      sig(2) = sig(2) - ei * (to-tu)*at/he 

c.... residual
42    i1=0
      do 41 ii=1,2
c....   B- matrix I
        call bmat02(bm,shp,vl,gr,ii,ityp)
c....   residual G = nodal forces -R = P - Bt*S 
        do i = 1,3
          do  k = 1,3
            p(i1+i) = p(i1+i) - bm(k,i,ii)*sig(k)*dv  
          end do  
        end do  
41     i1 = i1 + ndf
c     modify sign, here  P_xy=P-kv, S = -P_xy 
      ps(1) =  p(1)
      ps(2) = -p(2)
      ps(3) =  p(3) 
      ps(4) = -p(4)
      ps(5) =  p(5)
      ps(6) = -p(6) 

c.... modify with respect to reference/current configuration
      call rots02(gr,ps,ityp)

      if(isw.eq.13) go to 13
c.... Position
      gp(1) = shp(2,1)*xl(1,1) + shp(2,2)*xl(1,2)
      gp(2) = shp(2,1)*xl(2,1) + shp(2,2)*xl(2,2)

c.....Output Stresses (N, M, Q)
      mct = mct - 1
      if(mct.le.0) then
                     write(iow,4000)o,head
        if(ior.lt.0) write(*  ,4000)o,head
        mct = 50
      end if
4000  format(a1,20a4,//,
     +  2x,'2.PK-stress resultants 2D finite rotation beam',/,
     1  2x,'El',2x,'Mat',1x,'1-Coord',1x,'2-Coord',
     2  2X,'**Nx_1**',3X,'**Nx_2**',3X,'**Qz_1**',3X,'**Qz_2**',
     3  3X,'**My_1**',3X,'**My_2**',/)
                   write(iow,4001)n,ma,(gp(i),i=1,2),
     +                            ps(1),ps(4),ps(2),ps(5),ps(3),ps(6)
      if(ior.lt.0) write(  *,4001)n,ma,(gp(i),i=1,2),
     +                            ps(1),ps(4),ps(2),ps(5),ps(3),ps(6)
4001  format(1x,i4,i4,2f8.3,6(1x,e10.4))

      return

c.... mass matrices
5     sn = xl(2,2)-xl(2,1)
      cs = xl(1,2)-xl(1,1)
      sl = dsqrt(cs*cs + sn*sn)
      sn = sn/sl
      cs = cs/sl

c.... mass matrix lumped
      p(1) = 0.5d0 * d(18) * sl
      p(2) = p(1)
      p(3) = 0.5d0 * d(19) * sl
      p(4) = p(1)
      p(5) = p(2)
      p(6) = p(3)

c.... mass matrix consistent
      call pzero(s,6*6)

      ict = 1 
      if(ict.eq.1) then
c       consistent
        ra = d(18) * sl / 6.0d0
        ri = d(19) * sl / 6.0d0
        cm=2.d0 

      else if (ict.eq.2) then   
c       higher order mass matrix hughes p.446
        ra = d(18) * sl / 12.0d0
        ri = d(19) * sl / 12.0d0
        cm=5.d0 
      end if

c     riw=0 ! without rotational stiffness
      riw=1 ! with    rotational stiffness
 
      s(1,1) = ra*cm
      s(1,4) = ra
      s(2,2) = ra*cm
      s(2,5) = ra
      s(3,3) = ri*cm*riw
      s(3,6) = ri*riw
      s(4,1) = ra
      s(4,4) = ra*cm
      s(5,2) = ra
      s(5,5) = ra*cm
      s(6,3) = ri*riw
      s(6,6) = ri*cm*riw

      return

c...  Surface load:  constant follower loads 
c.... undeformed configuration
7     sn = xl(2,2)-xl(2,1)
      cs = xl(1,2)-xl(1,1)
      sl = dsqrt(cs*cs + sn*sn)
      sn = sn/sl
      cs = cs/sl

c.... local  displacements
      do k = 1,2
        vl(1,k) = cs*ul(1,k) + sn*ul(2,k)
        vl(2,k) =-sn*ul(1,k) + cs*ul(2,k)
        vl(3,k) = ul(3,k)
      end do

c.... Normal vector 
      dn1 =      - (vl(2,2)-vl(2,1)) / sl
      dn2 = 1.d0 + (vl(1,2)-vl(1,1)) / sl
c.... residual constant load qf 
      qf = d(1)
 
      p(1) =  qf * 0.5d0 * sl * dn1 
      p(2) =  qf * 0.5d0 * sl * dn2 
      p(4) =  qf * 0.5d0 * sl * dn1
      p(5) =  qf * 0.5d0 * sl * dn2
c.... make residual global
      call trans02(p,cs,sn,nst,ndf,2)
c.... tangent (symmetric terms)
      s(2,4) = - 0.5d0*qf
      s(4,2) = - 0.5d0*qf
      s(1,5) = + 0.5d0*qf
      s(5,1) = + 0.5d0*qf
c.... tangent (skewsymmetric terms) 
      s(1,2) =  - 0.5d0*qf
      s(2,1) =  + 0.5d0*qf
      s(4,5) =  + 0.5d0*qf
      s(5,4) =  - 0.5d0*qf
c.... make tangent global
      call trans02(s,cs,sn,nst,ndf,1)
      return

c.... damping matrix (only lumped)
12    sn = xl(2,2)-xl(2,1)
      cs = xl(1,2)-xl(1,1)
      sl = dsqrt(cs*cs + sn*sn)
      sn = sn/sl
      cs = cs/sl
      p(1) = 0.5d0 * d(20) * sl
      p(2) = p(1)
      p(3) = 0.5d0 * d(21) * sl
      p(4) = p(1)
      p(5) = p(2)
      p(6) = p(3)

c.... provisorial cdamp
      s(1,1) = p(1)
      s(2,2) = p(1)
      s(3,3) = p(3)
      s(4,4) = p(1)
      s(5,5) = p(1)
      s(6,6) = p(3)
      return

c.... plot N,Q,M diagrams on frame
13    mfp = abs(nfp)
      mfp = max(1,min(3,mfp))
      if(iplma(ma).eq.0)       return ! only if MATN
      nd2 = ndf+mfp
      p1 = ps(mfp)
      p2 = ps(nd2)
      klayf = 1
      if(flfp) then
        ccfp = max(abs(p1),abs(p2))
        ccfp1 = max(p1,p2)
        ccfp2 = min(p1,p2)
        xmaxf = max(xmaxf,ccfp1)
        xminf = min(xminf,ccfp2)
        cfp  = max(cfp,ccfp)
      else
        if(abs(p1).lt.eps.and.abs(p2).lt.eps) return
        smp=0.5d0*(p1+p2)  
        call pppcolf(smp)
c....   plot on mesh/deformed mesh
c....   global  displacements
        call trans02(ul,cs,sn,nst,ndf,2)
        sn = (xl(2,2)+scal*ul(2,2)) - (xl(2,1)+scal*ul(2,1))
        cs = (xl(1,2)+scal*ul(1,2)) - (xl(1,1)+scal*ul(1,1))
        sl = dsqrt(cs*cs + sn*sn)
        sn = sn/sl
        cs = cs/sl
        call pzero(xll,12)
        xll(1,1) = xl(1,1)+scal*ul(1,1)
        xll(2,1) = xl(2,1)+scal*ul(2,1)
        xll(1,2) = xl(1,nel)+scal*ul(1,nel)
        xll(2,2) = xl(2,nel)+scal*ul(2,nel)
        xll(1,3) = xl(1,nel)+scal*ul(1,nel) + sn*p2*cfp
        xll(2,3) = xl(2,nel)+scal*ul(2,nel) - cs*p2*cfp
        xll(1,4) = xl(1,1)  +scal*ul(1,1)   + sn*p1*cfp
        xll(2,4) = xl(2,1)  +scal*ul(2,1)   - cs*p1*cfp
        call plot9s(ixl,xll,3,4)
      end if
      return

c.... load vector and stiffness matrix from QLOA
22    ityp = d(9)
c.... Load terms
      call qload02(ql,d,aqloa,numel,n,mqloa,propq,prop,isw)

      ifol=ql(5)

      sn = xl(2,2)-xl(2,1)
      cs = xl(1,2)-xl(1,1)
      sl = dsqrt(cs*cs + sn*sn)
      sn = sn/sl
      cs = cs/sl
      dv = sl

c.... local  displacements
      do k = 1,2
        vl(1,k) = cs*ul(1,k) + sn*ul(2,k)
        vl(2,k) =-sn*ul(1,k) + cs*ul(2,k)
        vl(3,k) = ul(3,k)
      end do

c.... load case: conservative loads 
      if(ifol.eq.0) then
        p(1) =  - ( 2.*ql(3) +    ql(4)) * sl  / 6.
        p(4) =  - (    ql(3) + 2.*ql(4)) * sl  / 6.

        p(2) = - ( 2.*ql(1) +    ql(2)) * sl  / 6.
        p(5) = - (    ql(1) + 2.*ql(2)) * sl  / 6.

      end if

c.... load case: non-conservative loads 
      if(ifol.eq.1) then

c....   Normal vector 
        dn1 =      - (vl(2,2)-vl(2,1)) / sl
        dn2 = 1.d0 + (vl(1,2)-vl(1,1)) / sl

c....   residual constant load qf 
        qf = ql(1)
 
        p(1) =  qf * 0.5d0 * sl * dn1 
        p(2) =  qf * 0.5d0 * sl * dn2 
        p(4) =  qf * 0.5d0 * sl * dn1
        p(5) =  qf * 0.5d0 * sl * dn2

c....   add conservative loads
        p(1) =  p(1) - ( 2.*ql(3) +    ql(4)) * sl  / 6.
        p(4) =  p(4) - (    ql(3) + 2.*ql(4)) * sl  / 6.


c....   tangent (symmetric terms)
        s(2,4) = - 0.5d0*qf
        s(4,2) = - 0.5d0*qf
        s(1,5) = + 0.5d0*qf
        s(5,1) = + 0.5d0*qf

c....   tangent (skewsymmetric terms) 
        s(1,2) =  - 0.5d0*qf
        s(2,1) =  + 0.5d0*qf
        s(4,5) =  + 0.5d0*qf
        s(5,4) =  - 0.5d0*qf

      end if

c.... load case: temperature
      at = d(7)
      he = d(8)
      ea = d(15)
      ei = d(16)
      tu = ql(6)
      to = ql(7)

      if((tu*tu+to*to).eq.0.d0) goto 222

      call pzero(sigt,3)
      call pzero(pt,6)    
      call pzero(st,36)    
       
      sigt(1) =  ea * (to+tu)*at*0.5d0 
      sigt(2) =  ei * (to-tu)*at/he 

c.... shape functions
      shp(1,1) = -1.0d0/sl
      shp(1,2) =  1.0d0/sl
      shp(2,1) =  0.5d0
      shp(2,2) =  0.5d0

c.... elasticity matrix
      ea = d(15)
      ei = d(16)
      ga = d(17)
      cappa = 5.d0/6.d0
      gac = cappa*ga/(1.d0+cappa*sl*sl/12.d0*ga/ei)

      call pzero(dmat,9)
      dmat(1,1) =  ea
      dmat(2,2) =  ei
      dmat(3,3) =  gac

c.... strains,stresses only for gr
      call stre02(sig,shp,vl,dmat,gr,epoten,ityp)

c.... residual
      i1=0
      do 221 ii=1,2
c....   B- matrix I
        call bmat02(bm,shp,vl,gr,ii,ityp)
c....   residual G = P + BT*S_T 
        do i = 1,3
          do  k = 1,3
            pt(i1+i) = pt(i1+i) + bm(k,i,ii)*sigt(k)*dv
          end do   
        end do   
c....   and stiffness matrix K_sigma
        j1 = 0
        do jj = 1,ii
c....     K sigma
          call ksig02(sigt,shp,gr,st,nst,dv,ii,jj,i1,j1,ityp)
          j1 = j1 + ndf
        end do    
        i1 = i1 + ndf

221    continue

c.... upper part of stiffness matrix
      do i = 1,3
        do j = 1,3
          st(j,i+ndf) = st(i+ndf,j)
        end do
      end do

c.... add load vector and stiffness matrix for temperature
      do i = 1,6
        p(i) = p(i) + pt(i)
        do k = 1,6
          s(i,k) = s(i,k) - st(i,k)
        end do
      end do

c.... transform stiffness matrix
222   call trans02(s,cs,sn,nst,ndf,1)

c.... transform load vector
      call trans02(p,cs,sn,nst,ndf,2)
      
      return
      end
c
      subroutine bmat02(bm,shp,vl,gr,k,ityp)
c-----------------------------------------------------------------------
c.... nonlinear B-matrix for finite rotation beam
c     relation between 'beam theory-angle'and transformation angle
c     is defined by: x-y: sinf= cs    cosf= -sn
c     (c) w.wagner jan 83
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension bm(3,3,*),shp(2,2),vl(3,2),gr(5)
      call pzero(bm(1,1,k),9)
      dn  = shp(2,k)
      dns = shp(1,k)
      uks = gr(1)
      wks = gr(2)
      bks = gr(3)
      sb  = gr(4)
      cb  = gr(5)
      b   = shp(2,1)*vl(3,1) + shp(2,2)*vl(3,2)

      if(ityp.eq.1) then ! FR Green
        bm(1,1,k) = (1.0d0+uks)*dns
        bm(1,2,k) =  dns*wks
        bm(2,1,k) = -cb*bks*dns
        bm(2,2,k) = -sb*bks*dns
        bm(2,3,k) =(1.0d0+uks)*(sb*bks*dn-cb*dns)-wks*(cb*bks*dn+sb*dns)
        bm(3,1,k) = -sb*dns
        bm(3,2,k) =  cb*dns
        bm(3,3,k) = -((1.0d0+uks)*cb+wks*sb)*dn

      else if(ityp.eq.2) then ! FR Reissner
        bm(1,1,k) =  cb*dns
        bm(1,2,k) =  sb*dns
        bm(1,3,k) =  (wks*cb-(1.d0+uks)*sb)*dn
        bm(2,3,k) =  -dns
        bm(3,1,k) = -sb*dns
        bm(3,2,k) =  cb*dns
        bm(3,3,k) = -((1.0d0+uks)*cb+wks*sb)*dn

      else if(ityp.eq.3) then ! MR Green
        bm(1,1,k) = (1.0d0+uks)*dns
        bm(1,2,k) =  dns*wks
        bm(2,1,k) = -bks*dns
        bm(2,3,k) = -(1.0d0+uks)*dns
        bm(3,1,k) = -b*dns
        bm(3,2,k) =  dns
        bm(3,3,k) = -(1.0d0+uks)*dn

      else if(ityp.eq.4) then ! MR Green simplified
        bm(1,1,k) =  dns
        bm(1,2,k) =  dns*wks
        bm(2,3,k) = -dns
        bm(3,2,k) =  dns
        bm(3,3,k) = -dn

      else if(ityp.eq.5) then ! linear
        bm(1,1,k) =  dns
        bm(2,3,k) = -dns
        bm(3,2,k) =  dns
        bm(3,3,k) = -dn

      end if 
      return
      end
c
      subroutine stre02(sig,shp,vl,dmat,gr,epoten,ityp)
c-----------------------------------------------------------------------
c.... stresses for finite rotation beam
c     relation between 'shell theory-angle'and transformation angle
c     is defined by: x-y: sinf= cs    cosf= -sn
c     (c) w.wagner jan 83
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension sig(3),eps(3),shp(2,2),vl(3,2),dmat(3,3),gr(5)
c.... gradients
      uks = shp(1,1)*vl(1,1) + shp(1,2)*vl(1,2)
      wks = shp(1,1)*vl(2,1) + shp(1,2)*vl(2,2)
      bks = shp(1,1)*vl(3,1) + shp(1,2)*vl(3,2)
      b   = shp(2,1)*vl(3,1) + shp(2,2)*vl(3,2)
      sb  = dsin(b)
      cb  = dcos(b)
c.... save data for bmat  and ksigma
      gr(1) = uks
      gr(2) = wks
      gr(3) = bks
      gr(4) = sb
      gr(5) = cb
c...  strains
      if(ityp.eq.1) then ! FR Green
        eps(1) =  uks+0.5d0*(uks*uks+wks*wks)
        eps(2) = -((1.0d0+uks)*cb+wks*sb)*bks
        eps(3) = -(1.0d0+uks)*sb + wks*cb

      else if(ityp.eq.2) then ! FR Reissner
        eps(1) =  (1.d0+uks)*cb + wks*sb -1.d0
        eps(2) =  -bks
        eps(3) = -(1.0d0+uks)*sb + wks*cb

      else if(ityp.eq.3) then ! MR Green
        eps(1) =  uks+0.5d0*(uks*uks+wks*wks)
        eps(2) = -(1.0d0+uks)*bks
        eps(3) = -(1.0d0+uks)*b + wks

      else if(ityp.eq.4) then ! MR Green simplified
        eps(1) =  uks+0.5d0*(wks*wks)
        eps(2) = -bks
        eps(3) = -b + wks

      else if(ityp.eq.5) then ! linear
        eps(1) =  uks
        eps(2) = -bks
        eps(3) = -b + wks

      end if 

c.... stresses
        sig(1) = dmat(1,1)*eps(1)
        sig(2) = dmat(2,2)*eps(2)
        sig(3) = dmat(3,3)*eps(3)

c.... energy per lenght
      epoten=sig(1)*eps(1)+sig(2)*eps(2)+sig(3)*eps(3) 
cww ohne N      epoten=sig(2)*eps(2)+sig(3)*eps(3) 
      return
      end
c

      subroutine ksig02(sig,shp,gr,s,nst,dv,ii,kk,i1,k1,ityp)
c-----------------------------------------------------------------------
c.... matrix Ksigma for finite rotation beam
c     is defined by: x-y: sinf= cs    cosf= -sn
c     (c) w.wagner nov 88
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension sig(3),shp(2,2),gr(5),s(nst,nst)
c.... angle tranformation
c.... gradients
      uks = gr(1)
      wks = gr(2)
      bks = gr(3)
      sb  = gr(4)
      cb  = gr(5)
        dni   = shp(2,ii)
        dndsi = shp(1,ii)
        dnk   = shp(2,kk)
        dndsk = shp(1,kk)
c
      if(ityp.eq.1) then ! FR Green
c.... dof 1,1
        s(i1+1,k1+1) =  s(i1+1,k1+1) + sig(1)*dndsi*dndsk*dv
c.... dof 2,2
        s(i1+2,k1+2) =  s(i1+2,k1+2) + sig(1)*dndsi*dndsk*dv
c.... dof 3,3
        s(i1+3,k1+3) =  s(i1+3,k1+3)
     +               +  sig(2)*((1.0d0+uks)*cb+wks*sb)*bks*dni*dnk*dv
     +               +  sig(2)*((1.0d0+uks)*sb-wks*cb)*dni*dndsk*dv
     +               +  sig(2)*((1.0d0+uks)*sb-wks*cb)*dndsi*dnk*dv
     +               +  sig(3)*((1.0d0+uks)*sb-wks*cb)*dni*dnk*dv
c.... dof 3,1
        s(i1+3,k1+1) =  s(i1+3,k1+1)
     +               +  sig(2)*(sb*bks*dni - cb*dndsi)*dndsk*dv
     +               -  sig(3)*cb*dni*dndsk*dv
c.... dof 1,3
        s(i1+1,k1+3) =  s(i1+1,k1+3)
     +               +  sig(2)*dndsi*(sb*bks*dnk-cb*dndsk)*dv
     +               -  sig(3)*cb*dndsi*dnk*dv
c.... dof 3,2
        s(i1+3,k1+2) =  s(i1+3,k1+2)
     +               -  sig(2)*(cb*bks*dni+sb*dndsi)*dndsk*dv
     +               -  sig(3)*sb*dni*dndsk*dv
c.... dof 2,3
        s(i1+2,k1+3) =  s(i1+2,k1+3)
     +               -  sig(2)*dndsi*(cb*bks*dnk+sb*dndsk)*dv
     +               -  sig(3)*sb*dndsi*dnk*dv

      else if(ityp.eq.2) then ! FR Reissner

        s(i1+1,k1+3) =  s(i1+1,k1+3)
     +               -  (sig(1)*sb+sig(3)*cb)*dndsi*dnk*dv

        s(i1+3,k1+1) =  s(i1+3,k1+1)
     +               -  (sig(1)*sb+sig(3)*cb)*dni*dndsk*dv

        s(i1+2,k1+3) =  s(i1+2,k1+3)
     +               +  (sig(1)*cb-sig(3)*sb)*dndsi*dnk*dv

        s(i1+3,k1+2) =  s(i1+3,k1+2)
     +               +  (sig(1)*cb-sig(3)*sb)*dni*dndsk*dv

        s(i1+3,k1+3) =  s(i1+3,k1+3)
     +               +  (sig(1)*(-wks*sb-(1.d0+uks)*cb)
     +               +   sig(3)*(-wks*cb+(1.d0+uks)*sb) )*dni*dnk*dv   


      else if(ityp.eq.3) then ! MR Green
c.... dof 1,1
        s(i1+1,k1+1) =  s(i1+1,k1+1) + sig(1)*dndsi*dndsk*dv
c.... dof 2,2
        s(i1+2,k1+2) =  s(i1+2,k1+2) + sig(1)*dndsi*dndsk*dv
c.... dof 3,1
        s(i1+3,k1+1) =  s(i1+3,k1+1)
     +               -  sig(2)*dndsi*dndsk*dv
     +               -  sig(3)*dni*dndsk*dv
c.... dof 1,3
        s(i1+1,k1+3) =  s(i1+1,k1+3)
     +               -  sig(2)*dndsi*dndsk*dv
     +               -  sig(3)*dndsi*dnk*dv

      else if(ityp.eq.4) then ! MR Green simplified
c.... dof 2,2
        s(i1+2,k1+2) =  s(i1+2,k1+2) + sig(1)*dndsi*dndsk*dv

      else if(ityp.eq.5) then ! linear
c       nothing to do 

      end if
      return
      end
c
      subroutine rots02(gr,ps,ityp)
c-----------------------------------------------------------------------
c.... rotate Forces N,Q to current configuration
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension gr(5),ps(6)

      if(ityp.eq.5) return ! linear 
 
      uks=gr(1)
      wks=gr(2) 
 
c.....rotate to deformed vector d via beta 
      sb=gr(4)
      cb=gr(5)  
 
c.....rotate to deformed vector n via transformation from t_ref to t_cur 
c      ds=sqrt((1.d0+uks)*(1.d0+uks)+wks*wks)
c      sb=wks/ds   
c      cb=(1.d0+uks)/ds 
 
      ph1 = cb*ps(1)-sb*ps(2) 
      ph2 = cb*ps(2)+sb*ps(1) 
      ps(1)=ph1
      ps(2)=ph2
      ph4 = cb*ps(4)-sb*ps(5) 
      ph5 = cb*ps(5)+sb*ps(4) 
      ps(4)=ph4
      ps(5)=ph5

      return
      end  
c
      subroutine trans02(s,cs,snt,nst,ndf,itype)
c-----------------------------------------------------------------------
c.... itype: 1  Transform matrix s(nst,nst) from local to global
c            2  Transform vector s(nst,1)   from local to global
c            3  Transform vector s(nst,1)   from global to local
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension s(nst,*)
c
      if(cs.gt.0.999999d0) return
      sn = snt
      go to (1,2,3) itype
1     do 12 i = 1,nst,ndf
        j = i + 1
        do 11 n = 1,nst
          t      = s(n,i)*cs - s(n,j)*sn
          s(n,j) = s(n,i)*sn + s(n,j)*cs
          s(n,i) = t
11      continue
12    continue
      do 14 i = 1,nst,ndf
        j = i + 1
        do 13 n = 1,nst
          t      = s(i,n)*cs - s(j,n)*sn
          s(j,n) = s(i,n)*sn + s(j,n)*cs
          s(i,n) = t
13      continue
14    continue
      return
2     sn = -sn
3     do 31 i=1,nst,ndf
        j = i + 1
        t      =  s(i,1)*cs + s(j,1)*sn
        s(j,1) = -s(i,1)*sn + s(j,1)*cs
        s(i,1) =  t
31    continue
      return
      end
c
      subroutine qload02(ql,d,q,numel,n,mqloa,propq,prop,isw)
c----------------------------------------------------------
c
c      Purpose:  set loads from macros mate/qloa
c
c      Inputs:
c         d(*)        - d-array
c         q(numel,10) - loads for each element
c         numel       - Number of elements in mesh
c         n           - actual element
c         mqloa       - position where q is located in m
c         propq       - Current total proportional load level
c         isw         - macro to process
c
c      Outputs:
c
c         ql(1)      - load q1
c         ql(2)      - load q2
c         ql(3)      - load n1
c         ql(4)      - load n2
c         ql(5)      - ifol - 0 standard load,  1: follower load for q
c         ql(6)      - T_u
c         ql(7)      - T_o
c
c----------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension q(numel,10),d(*),ql(7)

      call pzero(ql,7)

      if(isw.eq.22) then
        if(mqloa.ne.1) then  
          ql(1) = q(n,1)*propq 
          ql(2) = q(n,2)*propq 
          ql(3) = q(n,3)*propq 
          ql(4) = q(n,4)*propq 
          ql(5) = q(n,5)
          ql(6) = q(n,6)*propq 
          ql(7) = q(n,7)*propq 
        end if
      else if(isw.eq.4.or.isw.eq.13) then
        ql(1) = d(10)*prop 
        ql(2) = d(11)*prop 
        ql(3) = d(12)*prop 
        ql(4) = d(13)*prop 
        if(mqloa.ne.1) then  
          ql(5)  = q(n,5)
          ifol   = ql(5)
          if(ifol.eq.0) then
            ql(1) = ql(1) + q(n,1)*propq 
            ql(2) = ql(2) + q(n,2)*propq 
          else
            ql(1) = q(n,1)*propq
            ql(2) = q(n,2)*propq 
          end if  
          ql(3) = ql(3) + q(n,3)*propq 
          ql(4) = ql(4) + q(n,4)*propq 
          ql(6) =         q(n,6)*propq 
          ql(7) =         q(n,7)*propq 
        end if
      else
        ql(1)  = d(10)*prop 
        ql(2)  = d(11)*prop 
        ql(3)  = d(12)*prop 
        ql(4)  = d(13)*prop 
      end if  

      return
      end
