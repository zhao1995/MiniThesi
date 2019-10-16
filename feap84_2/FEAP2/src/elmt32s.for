      subroutine elmt32(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c**********************************************************************|
c.....stabilized four node linear Mindlin-Reissner plate element       |
c.... Material parameters (input)                                      |
c      1: E     = Young's modulus                                      |
c      2: nu    = Poisson's ratio                                      |
c      3: h     = plate thickness                                      |
c      4: q     = transverse load                                      |
c      5: rho   = density                                              |
c      6: istab = 1/0   =  with / without stabilization                |
c      7: istyp = 1/0   =  stress resultants at nodes / center         |
c                                                                      |
c      signs of stresses, see lecture CTWM UKA or FEM TUD              |
c      signs of displacements:  w,phix,phiy                            |
c                                                                      |
c >>   cappa=5/6  isw=1+3!!                                            |
c      cappa=1000 -> Kirchhoff-theory                                  |
c                                                                      |
c      3=TANG,4=STRE,5=LMAS+CMAS,6=REAC,8=STRE,13=STR1                 | 
c                                                                      |   
c      F. Gruttmann TUD 2009                                           |
c**********************************************************************|
c      open  WW 10.12.03                                               |
c      mass matrix in closed form                                      |
c      winkler foundation                                              |
c      stress calculation at nodes or gp_center                        |
c**********************************************************************|
      USE bdata
      USE cdata
      USE eldata
      USE iofile
      USE pdata7
      USE pdata10
      USE prlod
      USE pltran
      USE strnam

      implicit double precision (a-h,o-z)
      dimension xl(ndm,*),tl(*),ix(*),ul(ndf,*),s(nst,nst),p(nst),d(*),
     1          shp(3,4),sig(5),ep(9),beta(9),gp(2),h1(*),h2(*),h3(*),
     2          b(5,3,4),im(4),il(4),xn(2,4),yl(3,4),
     3          xs(2,2),sx(2,2),g0(2,2),g1(2,2),a1(4),a2(4),hh(4),
     4          gm(2,2),hm(2,2),hs(2,2),rh(2,2),gs(2,3,4),gth(2,3),
     5          xsi(4),eta(4),sig1(5,4),sg(4),tg(4),wg(4),shp1(3,4)
      data a1    / -0.25d0,  0.25d0, 0.25d0, -0.25d0 /
      data a2    / -0.25d0, -0.25d0, 0.25d0,  0.25d0 /
      data hh    /  0.25d0, -0.25d0, 0.25d0, -0.25d0 /
      data im  /-2, 2, 4,-4/,  il /-1,-3, 3, 1/
      data xsi /-1, 1, 1,-1/, eta /-1,-1, 1, 1/
      save c00,c11,c12,c22,c33,c44,maold
c.... go to correct array processor
      go to(1,2,3,3,5,3,2,3,2,2,3,2,3) isw
c.... input material properties
1     if(ior.lt.0) write(*,1001)
      call dinput(d,6)
                   write(iow,1002) (d(i),i=1,6)
      if(ior.lt.0) write(  *,1002) (d(i),i=1,6)
      if(ndm.eq.3) ipla = 1
c.... elasticity matrix      
      c00 = d(1)*d(3)*d(3)*d(3)/12.d0
      c11 = c00/(1.0d0-d(2)*d(2))
      c22 = c11
      c12 = c11*d(2) 
      c33 = c11*(1.0d0-d(2))*0.5d0
      cappa = 5.d0/6.d0
c     cappa = 1000.d0
      c44 = 0.5d0*d(1)*d(3)/(1.0d0+d(2))*cappa
      maold = ma
c.... description of stresses  
      strsus( 1) = '  MOMENT M_xx  '
      strsus( 2) = '  MOMENT M_xy  '
      strsus( 3) = '  MOMENT M_yy  '
      strsus( 4) = '               '
      strsus( 5) = '  MOMENT M_11  '
      strsus( 6) = '  MOMENT M_22  '
      strsus( 7) = '  ANGLE Phi_1  '
      strsus( 8) = '  SHEAR FORCE_x'
      strsus( 9) = '  SHEAR FORCE_y'
      strsus(10) = '               '
      strsus(11) = '               '
c
2     return
3     istab = d(6)
      istyp = d(7)
c.... elasticity matrix      
      if(ma.ne. maold) then
        c00 = d(1)*d(3)*d(3)*d(3)/12.d0
        c11 = c00/(1.0d0-d(2)*d(2))
        c22 = c11
        c12 = c11*d(2) 
        c33 = c11*(1.0d0-d(2))*0.5d0
        cappa = 5.d0/6.d0
c       cappa = 1000.d0
        c44 = 0.5d0*d(1)*d(3)/(1.0d0+d(2))*cappa
        maold = ma
      end if 
      call pzero (g0,4)
      call pzero (g1,4)
      do i = 1,2
        xn(i,1) = ( xl(i,4) - xl(i,1) )/2.d0
        xn(i,2) = ( xl(i,2) - xl(i,1) )/2.d0
        xn(i,3) = ( xl(i,3) - xl(i,2) )/2.d0
        xn(i,4) = ( xl(i,3) - xl(i,4) )/2.d0
       do K = 1,4
        g0(i,1) = g0(i,1) + a1(K) * xl(i,K)
        g0(i,2) = g0(i,2) + a2(K) * xl(i,K)
        g1(i,1) = g1(i,1) + hh(K) * xl(i,K)
        g1(i,2) = g1(i,2) + hh(K) * xl(i,K)
       end do
      end do
c
      xs(1,1) = g0(1,1)
      xs(2,2) = g0(2,2)
      xs(1,2) = g0(2,1)
      xs(2,1) = g0(1,2)
      xsj = xs(1,1)*xs(2,2) - xs(1,2)*xs(2,1)
      if(xsj.le.0.d0)call drawmess('negative Jacobian',1,-2)
      da = 4.d0*xsj
      r1 = ( g0(1,1)*g1(2,2) - g0(2,1)*g1(1,1) ) / xsj
      r2 = ( g0(2,2)*g1(1,1) - g0(1,2)*g1(2,2) ) / xsj
      sx(1,1) = xs(2,2)/xsj
      sx(2,2) = xs(1,1)/xsj
      sx(1,2) =-xs(1,2)/xsj
      sx(2,1) =-xs(2,1)/xsj
c
c     compute stabilization array
      d11 = c11*da
      d12 = c12*da
      d22 = d11
      d33 = c33*da
      d44 = c44*da 
      factq = istab*c44*da*da/9.d0 
      factm = istab*c00*da*da/9.d0

      f11  = 1.d0-r2*r2/3.d0
      f22  = 1.d0-r1*r1/3.d0
      f12  =     -r1*r2/3.d0
      g11  = (xs(1,1)*xs(1,1)+xs(1,2)*xs(1,2)) 
      g12  = (xs(1,1)*xs(2,1)+xs(1,2)*xs(2,2)) 
      g22  = (xs(2,2)*xs(2,2)+xs(2,1)*xs(2,1)) 
      da3  = da/3.d0

      h44 = da3*f11*g11*g11
      h55 = da3*f22*g22*g22
      h45 = da3*f12*(g12*g12-d(2)*xsj*xsj)
      det = h44*h55 - h45*h45
      hm(1,1) =  factm*h55/det
      hm(2,2) =  factm*h44/det
      hm(1,2) = -factm*h45/det
      hm(2,1) =  hm(1,2)
      gm(1,1) = -xs(1,2)
      gm(2,2) =  xs(2,1)
      gm(1,2) =  xs(1,1)
      gm(2,1) = -xs(2,2)
      call matmulf (hm,gm,2,2,2,rh)
      call mttmul (gm,rh,2,2,2,hm)

      h11 = da3*f11*g11
      h22 = da3*f22*g22
      h12 = da3*f12*g12
      det = h11*h22 - h12*h12
      hs(1,1) = factq*h22/det
      hs(2,2) = factq*h11/det
      hs(1,2) =-factq*h12/det
      hs(2,1) = hs(1,2)
       
      do inode = 1,nel
       ir0 = 3*(inode-1)
       shp(1,inode) = sx(1,1)*a1(inode)+sx(1,2)*a2(inode)
       shp(2,inode) = sx(2,1)*a1(inode)+sx(2,2)*a2(inode)
      if(dabs(d(4)).gt.1.d-20) 
     1 p(ir0+1) = d(4)*prop*da*(0.75d0 + r1*a1(inode)+r2*a2(inode))/3.d0   
c
        mm  = iabs(im(inode))
        ll  = iabs(il(inode))
        b11 = isign(1,im(inode))*xn(1,mm)
        b12 = isign(1,im(inode))*xn(2,mm)
        b21 = isign(1,il(inode))*xn(1,ll)
        b22 = isign(1,il(inode))*xn(2,ll)

        gs(1,1,inode) =  hh(inode) - r2*a1(inode) - r1*a2(inode)
        gs(2,1,inode) =  gs(1,1,inode)
        gs(1,2,inode) = -b12*hh(inode)+b12*r2*a1(inode)+b22*r1*a2(inode)
        gs(1,3,inode) =  b11*hh(inode)-b11*r2*a1(inode)-b21*r1*a2(inode)
        gs(2,2,inode) = -b22*hh(inode)+b12*r2*a1(inode)+b22*r1*a2(inode)
        gs(2,3,inode) =  b21*hh(inode)-b11*r2*a1(inode)-b21*r1*a2(inode)
     
        b11 = b11*a1(inode)
        b12 = b12*a1(inode)
        b21 = b21*a2(inode)
        b22 = b22*a2(inode)

        b(4,1,inode) =  shp(1,inode)
        b(5,1,inode) =  shp(2,inode)
        b(4,2,inode) = -sx(1,1)*b12 - sx(1,2)*b22
        b(5,2,inode) = -sx(2,1)*b12 - sx(2,2)*b22
        b(4,3,inode) =  sx(1,1)*b11 + sx(1,2)*b21
        b(5,3,inode) =  sx(2,1)*b11 + sx(2,2)*b21

       if(isw.eq.3.or.isw.eq.6)then

        do i = 1,3
         gth(1,i) = gs(1,i,inode)*hs(1,1)+gs(2,i,inode)*hs(2,1)
         gth(2,i) = gs(1,i,inode)*hs(1,2)+gs(2,i,inode)*hs(2,2)
        end do

        do jnode = 1,inode
         jc0 = 3*(jnode-1)
         fm  = gs(1,1,inode)*gs(1,1,jnode)
           s(ir0+2,jc0+2) = d22*shp(2,inode)*shp(2,jnode)+hm(1,1)*fm
     1                    + d33*shp(1,inode)*shp(1,jnode)
           s(ir0+2,jc0+3) =-d12*shp(2,inode)*shp(1,jnode)+hm(1,2)*fm
     1                    - d33*shp(1,inode)*shp(2,jnode)
           s(ir0+3,jc0+2) =-d12*shp(1,inode)*shp(2,jnode)+hm(2,1)*fm
     1                    - d33*shp(2,inode)*shp(1,jnode)
           s(ir0+3,jc0+3) = d11*shp(1,inode)*shp(1,jnode)+hm(2,2)*fm
     1                    + d33*shp(2,inode)*shp(2,jnode)
         do i = 1,3
          do j = 1,3
           s(ir0+i,jc0+j) = s(ir0+i,jc0+j) 
     1     + (b(4,i,inode)*b(4,j,jnode)+b(5,i,inode)*b(5,j,jnode))*d44
     2     +  gth(1,i)*gs(1,j,jnode)+gth(2,i)*gs(2,j,jnode)
          end do
         end do 
        end do
      
       end if 
      
      end do   

      if(isw.eq.4.or.isw.eq.8.or.isw.eq.13)goto 4
c.....upper part of K
      do i = 2,nst
       do j = 1,i-1
         s(j,i) = s(i,j)
       end do
      end do    
c      if(isw.eq.6) then ! if line is comment: than R=0 after 2nd tang
c....   compute the element residual vector: R = P - K*v
        do i = 1,nst
          do j = 1,nst
            p(i) = p(i) - s(i,j)*ul(j,1)
          end do
        end do   
c      end if
      return 
c
c.....stress resultants (M and Q) at element center
4     istv = 9
      call pzero(ep,9)
      do k = 1,4
       ep(1)=ep(1)+shp(1,k)*ul(3,k) 
       ep(2)=ep(2)-shp(2,k)*ul(2,k) 
       ep(3)=ep(3)-shp(1,k)*ul(2,k)+shp(2,k)*ul(3,k)  
       ep(4)=ep(4)+b(4,1,k)*ul(1,k)+b(4,2,k)*ul(2,k)+b(4,3,k)*ul(3,k)  
       ep(5)=ep(5)+b(5,1,k)*ul(1,k)+b(5,2,k)*ul(2,k)+b(5,3,k)*ul(3,k)
       ep(6)=ep(6)+gs(1,1,k)*ul(2,k) 
       ep(7)=ep(7)+gs(1,1,k)*ul(3,k) 
       ep(8)=ep(8)+gs(1,1,k)*ul(1,k)+gs(1,2,k)*ul(2,k)+gs(1,3,k)*ul(3,k)
       ep(9)=ep(9)+gs(2,1,k)*ul(1,k)+gs(2,2,k)*ul(2,k)+gs(2,3,k)*ul(3,k) 
      end do
      beta(1) = c11*ep(1)+c12*ep(2)
      beta(2) = c12*ep(1)+c22*ep(2)
      beta(3) = c33*ep(3)
      beta(4) = c44*ep(4)
      beta(5) = c44*ep(5)
      beta(6) = (rh(1,1)*ep(6)+rh(1,2)*ep(7))/da3
      beta(7) = (rh(2,1)*ep(6)+rh(2,2)*ep(7))/da3
      beta(8) = (hs(1,1)*ep(8)+hs(1,2)*ep(9))/da3
      beta(9) = (hs(2,1)*ep(8)+hs(2,2)*ep(9))/da3

      xsiq = r1/3.d0
      etaq = r2/3.d0
      if(istyp.eq.0) then
c....   stresses at center
        sig(1) = beta(1) - xs(1,1)*xs(1,1)*etaq*beta(6) 
     1                   - xs(2,1)*xs(2,1)*xsiq*beta(7)
        sig(2) = beta(2) - xs(1,2)*xs(1,2)*etaq*beta(6) 
     1                   - xs(2,2)*xs(2,2)*xsiq*beta(7)
        sig(3) = beta(3) - xs(1,1)*xs(1,2)*etaq*beta(6) 
     1                   - xs(2,1)*xs(2,2)*xsiq*beta(7)
        sig(4) = beta(4) - xs(1,1)*etaq*beta(8) - xs(2,1)*xsiq*beta(9)
        sig(5) = beta(5) - xs(1,2)*etaq*beta(8) - xs(2,2)*xsiq*beta(9)
      elseif(istyp.eq.1) then
c....   stresses at nodes
        do i = 1,4
         sig1(1,i) = beta(1) + xs(1,1)*xs(1,1)*(eta(i)-etaq)*beta(6) 
     1                       + xs(2,1)*xs(2,1)*(xsi(i)-xsiq)*beta(7)
         sig1(2,i) = beta(2) + xs(1,2)*xs(1,2)*(eta(i)-etaq)*beta(6) 
     1                       + xs(2,2)*xs(2,2)*(xsi(i)-xsiq)*beta(7)
         sig1(3,i) = beta(3) + xs(1,1)*xs(1,2)*(eta(i)-etaq)*beta(6) 
     1                       + xs(2,1)*xs(2,2)*(xsi(i)-xsiq)*beta(7)
         sig1(4,i) = beta(4) + xs(1,1)*(eta(i)-etaq)*beta(8) 
     1                       + xs(2,1)*(xsi(i)-xsiq)*beta(9)
         sig1(5,i) = beta(5) + xs(1,2)*(eta(i)-etaq)*beta(8) 
     1                       + xs(2,2)*(xsi(i)-xsiq)*beta(9)
        end do
      end if 
      if(isw.eq.4) then     
c....   print stresses
        mct = mct - 1
        if(mct.gt.0) go to 42
                     write(iow,2001) o,head
        if(ior.lt.0) write(*  ,2001) o,head
        mct = 50
   42   continue
        if(istyp.eq.0) then
         gp(1) = 0.25d0*(xl(1,1)+xl(1,2)+xl(1,3)+xl(1,4)) 
         gp(2) = 0.25d0*(xl(2,1)+xl(2,2)+xl(2,3)+xl(2,4)) 
                      write(iow,2002) n,ma,(gp(i),i=1,2),(sig(i),i=1,5)
         if(ior.lt.0) write(*  ,2002) n,ma,(gp(i),i=1,2),(sig(i),i=1,5)
       elseif(istyp.eq.1) then
        mct = mct - 3
        do j = 1,4
         gp(1) = xl(1,j) 
         gp(2) = xl(2,j) 
                    write(iow,2002) n,ma,(gp(i),i=1,2),(sig1(i,j),i=1,5)
         if(ior.lt.0) write(*,2002) n,ma,(gp(i),i=1,2),(sig1(i,j),i=1,5)
        end do
       end if 
      else if(isw.eq.8) then     
c....  plot stresses
       if(istyp.eq.0) then
         call plot32  (ix,strea,strea(1+numnp),sig,da,numnp)
       else if(istyp.eq.1) then
         if(iplma(ma).eq.0)  return ! only if MATN
         call plot32_1(ix,strea,strea(1+numnp),sig1,da,r1,r2,numnp)
       end if
      else if(isw.eq.13) then     
c..... plot stress resultants from center gauss-point without averaging 
       if(flfp) then
c....   calculate extreme values
        xmaxf = max(xmaxf,sig(nfp))
        xminf = min(xminf,sig(nfp))
       else
c....   plot stresses elementwise in one color      
c....   color      
        call pppcolf(sig(nfp))
c....   coordinates
        call pzero(yl,12)
        do i=1,ndm
          do j=1,4
            yl(i,j) = xl(i,j) 
          end do
        end do
c.....  transform for rot        
        call plxtrn(yl,tra,vr,3,4)
c...... plot element
        call plot9 (iel,ix,yl,ndm,nel,1)
       end if
      end if 
      return
c
c.... compute lumped mass matrix (l=2=exakt)
5     l = 2
      call pgauss(l,lint,sg,tg,wg)
      dh = d(3)*d(3)/12.d0
      do 50 l = 1,lint
c....   compute shape functions
        call shape(sg(l),tg(l),xl,shp1,xsj,ndm,nel,ix,.false.)
        da = wg(l)*xsj*d(3)*d(5)
        i1 = 1
        do 51 i = 1,nel
          p(i1  ) = p(i1) + shp1(3,i)*da
          p(i1+1) = p(i1)*dh
          p(i1+2) = p(i1)*dh
          j1 = i1
          do 52 j = i,nel
            sm = shp1(3,i)*shp1(3,j)*da
            s(i1  ,j1  ) = s(i1  ,j1  ) + sm
            s(i1+1,j1+1) = s(i1+1,j1+1) + sm*dh
            s(i1+2,j1+2) = s(i1+2,j1+2) + sm*dh
            j1 = j1 + ndf
52        end do
          i1 = i1 + ndf
51       end do
50     end do
c.....upper part of M
      do 56 i = 1,nst
        do 56 j = 1,i
   56     s(i,j) =s(j,i)
c
      return
c
c.... formats
1001  format(' Input: E  nu  h  q rho istab istyp',$)
1002  format(5x,'Materialdata for plate element:',/,
     + 5x,'elastic modulus............................'   ,f15.4,/,
     + 5x,'poissons ratio................................',f12.4,/,
     + 5x,'thickness.....................................',f12.4,/,
     + 5x,'element load (transverse).....................',f12.4,/,
     + 5x,'density ......................................',f12.4,/,
     + 5x,'istab 0 = without, 1 = with stabilization.....',f4.0) 
cww  + 5x,'istyp 0/1 stress resultants at center/nodes..',f4.0)
2001  format(a1,20a4,//,2x,'E L E M E N T   S T R E S S E S',//,
     1 2x,'EL',1x,'MAT',1x,'1-COORD',1x,'2-COORD',
     2 2X,'***MX***',3X,'***MY***',3X,'***MXY***',
     3 2X,'***QX***',3X,'***QY***',/)
2002  format(1x,i3,i3,2f8.3,1x,5e11.4)
      end
c----------------------------------------------------------------------
      subroutine plot32(ix,dt,st,sig,da,numnp)
c***********************************************************************
c.....Plot   mx(1)   mxy(2)  my(3)      m1(5)  m2(6) phi_1(7)
c            qx(8)           qy(9)
c.... sig from element center
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension dt(numnp),st(numnp,*),sig(5),ix(1)
        xsji = da*0.25d0 
        do 10  j = 1,4
           ii = abs(ix(j))
           if(ii.eq.0) go to 10
             dt(ii)   = dt(ii)   + xsji
             st(ii,1) = st(ii,1) + sig(1)*xsji
             st(ii,2) = st(ii,2) + sig(3)*xsji
             st(ii,3) = st(ii,3) + sig(2)*xsji
             st(ii,8) = st(ii,8) + sig(4)*xsji
             st(ii,9) = st(ii,9) + sig(5)*xsji
10      continue
      return
      end
c----------------------------------------------------------------------
      subroutine plot32_1(ix,dt,st,sig1,da,r1,r2,numnp)
c***********************************************************************
c.....Plot   mx(1)   mxy(2)  my(3)      m1(5)  m2(6) phi_1(7)
c            qx(8)           qy(9)
c.... sig1 from element nodes
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension dt(numnp),st(numnp,*),sig1(5,4),ix(1),xsi(4),eta(4)
      data xsi /-1, 1, 1,-1/, eta /-1,-1, 1, 1/
        xsj0 = da*0.25d0 
        xsj1 = r1*xsj0
        xsj2 = r2*xsj0
        do 10  j = 1,4
           ii = abs(ix(j))
           if(ii.eq.0) go to 10
             xsji = xsj0 + xsi(j)*xsj1 + eta(j)*xsj2 
             dt(ii)   = dt(ii)   + xsji
             st(ii,1) = st(ii,1) + sig1(1,j)*xsji
             st(ii,2) = st(ii,2) + sig1(3,j)*xsji
             st(ii,3) = st(ii,3) + sig1(2,j)*xsji
             st(ii,8) = st(ii,8) + sig1(4,j)*xsji
             st(ii,9) = st(ii,9) + sig1(5,j)*xsji
10      continue
      return
      end
