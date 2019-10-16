      subroutine elmt06(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c--------------------------------------------------------------------
c     Plane stress/strain Pian Sumihara Element
c--------------------------------------------------------------------
c     Linear elastic material
c     Formed as a one-point element with rank-2 update
c
c
c     d( 1) = e/(1-nu^2) or e(1-nu)/((1+nu)(1-nu)] 
c     d( 2) = e/(1-nu^2) or e(1-nu)/((1+nu)(1-nu)] 
c     d( 3) = e/(1+nu) 
c     d( 4) = e/(2(1+nu)=G
c     d( 5) = 1/e or (1-nu)/G
c     d( 6) = -nu/e or nu/G
c     d( 7) =
c     d( 8) = e
c     d( 9) = nu 
c     d(10) = rho = gamma/g [F/L**3]/ [S/T**2]
c     d(11) = thick h (def=1)
c     d(12) = is = ityp 1,2
c     d(13) = 0 or nu
c     d(14) = -nu/e
c
c     Version for FEAP               ww 22.11.91
c--------------------------------------------------------------------
      USE bdata
      USE cdata
      USE eldata
      USE errin1
      USE errin2
      USE errin3
      USE iofile
      USE pdata10
      USE strnam
      implicit double precision (a-h,o-z)
      dimension   xl(ndm,*),sig(7),ix(*),d(*),ul(ndf,*),s(nst,*),p(*)
      dimension h1(*),h2(*),h3(*)
      character*4 wd(2)
      common /elcom061/ xs,xt,xh,ys,yt,yh,xj0,xj1,xj2,r1,r2,
     1          x1(4),x2(4),y1(4),y2(4),vol
!$OMP THREADPRIVATE (/elcom061/)  
      common /elcom062/ ax(4),bx(4),cx(4),ay(4),by(4),cy(4),a1(3),a2(3),
     1            beta(6)
!$OMP THREADPRIVATE (/elcom062/)  
      data wd/'ress','rain'/
c.... go to correct array processor
      go to(1,2,3,3,5,3,2,3,3,2,2,2,2,14), isw
      return
c.... input/output material properties
1     if(ior.lt.0) write(*,3000)
      call dinput(d(8),5)
      is = d(12)
      is = max(1,min(is,2))
      d(12)= is
      d(2) = d(9)*d(8)/(1.+d(9))/(1.-is * d(9))
      d(4) = d(8)/(1.+d(9))
      d(3) = d(4)/2.
      d(1) = d(2) + d(4)
      d(14) = -d(9)/d(8)
c.... set parameters for plane stress (is = 1)
      if(is.eq.1) then
        if(d(11).le.0.0d0) d(11) = 1.0
        d(5) =  1./d(8)
        d(6) = -d(9)/d(8)
        d(13)=  0.0
c....   set parameters for plane strain (is = 2)
      else
        d(11)=  1.0
        d(5) =  (1.-d(9))/d(4)
        d(6) = -d(9)/d(4)
        d(13)=  d(9)
      end if
      write(iow,2000) wd(is),(d(i),i=8,11)
      if(ior.lt.0) write(*,2000) wd(is),(d(i),i=8,11)
c.... description of stresses  
      strsus( 1) =  '  Stress  S_xx '
      strsus( 2) =  '  Stress  S_xy '
      strsus( 3) =  '  Stress  S_yy '
      strsus( 4) =  '  Stress  S_zz '
      strsus( 5) =  '  Stress  S_1  '
      strsus( 6) =  '  Stress  S_2  '
      strsus( 7) =  '  ANGLE phi_1  '
      do is =8,25
        strsus(is) = '               '
      end do
c
2     return
c.... Compute the Pian-Sumihara arrays for elastic
c
c.... compute jacobian
   3  xs = (-xl(1,1)+xl(1,2)+xl(1,3)-xl(1,4))/4.
      ys = (-xl(2,1)+xl(2,2)+xl(2,3)-xl(2,4))/4.
      xt = (-xl(1,1)-xl(1,2)+xl(1,3)+xl(1,4))/4.
      yt = (-xl(2,1)-xl(2,2)+xl(2,3)+xl(2,4))/4.
      xh = ( xl(1,1)-xl(1,2)+xl(1,3)-xl(1,4))/4.
      yh = ( xl(2,1)-xl(2,2)+xl(2,3)-xl(2,4))/4.
      xj0 = xs*yt - xt*ys
      xj1 = xs*yh - xh*ys
      xj2 = xh*yt - xt*yh
      ssa = xj1/xj0/3.
      tta = xj2/xj0/3.
c.... form stiffness for elastic part and compute the beta parameters
      call pian06(d,ul,s,p,nst,ndf,isw)
      if(isw.eq.4) go to 4
      if(isw.eq.8) go to 8
      if(isw.eq.9) go to 9
c
c.... Compute the residual (already above in pian06) 
c
c.... compute symetric part of s
      do 334 i = 1,nst-1
        do 332 j = i+1,nst
          s(j,i) = s(i,j)
 332    continue
 334  continue
      return
c
c.... compute the stresses
c
   4  is = d(12)
c.... compute the stresses at the center
      ssg = - ssa*beta(5)
      ttg = - tta*beta(4)
      sig(1) = beta(1) + a1(1)*ttg + a2(1)*ssg
      sig(3) = beta(2) + a1(2)*ttg + a2(2)*ssg
      sig(2) = beta(3) + a1(3)*ttg + a2(3)*ssg
      sig(4) = d(13)*(sig(1)+sig(2))
      call pstres(sig,p1,p2,p3)
      xx = (xl(1,1)+xl(1,2)+xl(1,3)+xl(1,4))/4.0
      yy = (xl(2,1)+xl(2,2)+xl(2,3)+xl(2,4))/4.0
      mct = mct - 1
      if(mct.gt.0) go to 450
      write(iow,2001) o,head,wd(is)
      if(ior.lt.0) write(*,2001) o,head,wd(is)
      mct = 25
450   write(iow,2002) n,ma,xx,yy,p1,p2,p3,(sig(i),i=1,4)
      if(ior.lt.0) write(*,2002) n,ma,xx,yy,p1,p2,p3,
     1 (sig(i),i=1,4)
      return
c
c.... compute a lumped mass matrix
c
   5  xs =(-xl(1,1)+xl(1,2)+xl(1,3)-xl(1,4))/4.
      ys =(-xl(2,1)+xl(2,2)+xl(2,3)-xl(2,4))/4.
      xt =(-xl(1,1)-xl(1,2)+xl(1,3)+xl(1,4))/4.
      yt =(-xl(2,1)-xl(2,2)+xl(2,3)+xl(2,4))/4.
      xst=(+xl(1,1)-xl(1,2)+xl(1,3)-xl(1,4))/4.
      yst=(+xl(2,1)-xl(2,2)+xl(2,3)-xl(2,4))/4.
      d11 = d(10)*d(11)
      xj0 = (xs*yt-xt*ys)*d11
      xj1 = (xs*yst-xst*ys)*d11/3.0
      xj2 = (xst*yt-xt*yst)*d11/3.0
      p(1) = xj0-xj1-xj2
      i1 = 1+ndf
      p(i1) = xj0+xj1-xj2
      i1 = i1 + ndf
      p(i1) = xj0+xj1+xj2
      i1 = i1 + ndf
      p(i1) = xj0-xj1+xj2
      do 510 i = 2,nst,ndf
510   p(i) = p(i-1)
      return
c
c.... compute the nodal stress values
   8  istv =  7
      ityp = d(12)
      if(iplma(ma).eq.0)  return ! only if MATN
      call stcn06(ix,d,ssa,tta,strea,strea(1+numnp),numnp,ityp)
      return
c.... compute the errors
   9  continue
      call ster06(ix,d,xl,ssa,tta,strea,strea(1+numnp),numnp,ityp,e_ome,
     +            numel,ndm)

      return
c.... calculate stresses sig(i) at center of element
c.... compute jacobian
  14  xs = (-xl(1,1)+xl(1,2)+xl(1,3)-xl(1,4))/4.
      ys = (-xl(2,1)+xl(2,2)+xl(2,3)-xl(2,4))/4.
      xt = (-xl(1,1)-xl(1,2)+xl(1,3)+xl(1,4))/4.
      yt = (-xl(2,1)-xl(2,2)+xl(2,3)+xl(2,4))/4.
      xh = ( xl(1,1)-xl(1,2)+xl(1,3)-xl(1,4))/4.
      yh = ( xl(2,1)-xl(2,2)+xl(2,3)-xl(2,4))/4.
      xj0 = xs*yt - xt*ys
      xj1 = xs*yh - xh*ys
      xj2 = xh*yt - xt*yh
      ssa = xj1/xj0/3.
      tta = xj2/xj0/3.
c.... form stiffness for elastic part and compute the beta parameters
      call pian06(d,ul,s,p,nst,ndf,isw)
      is = d(12)
c.... compute the stresses at the center
      ssg = - ssa*beta(5)
      ttg = - tta*beta(4)
      sig(1) = beta(1) + a1(1)*ttg + a2(1)*ssg
      sig(3) = beta(2) + a1(2)*ttg + a2(2)*ssg
      sig(2) = beta(3) + a1(3)*ttg + a2(3)*ssg
      sig(4) = d(13)*(sig(1)+sig(2))
      call pstres(sig,p1,p2,p3)
      sig(5) = p1
      sig(6) = p2
      sig(7) = p3
c
      if(nfp.gt.7 .or.nfp.lt.1) return
      if(flfp) then
c....   calculate extreme values
        xmaxf = max(xmaxf,sig(nfp))
        xminf = min(xminf,sig(nfp))
      else
c....   plot stresses 
c....   color      
        call pppcolf(sig(nfp)) 
c...... plot element
        call plot9(iel,ix,xl,ndm,nel,1)
      end if
      return
c
c.... format statements
2000  format(5x,'plane st',a4,' element'//10x,'modulus      =',e12.5/
     1  10x,'poisson ratio=', f8.5/10x,'mass density =',e12.5/
     2  10x,'thickness    =', f8.5)
2001  format(a1,20a4//5x,'plane st',a4,' stresses'//'  element',1x,
     1'material',5x,'1-coord',5x,'2-coord',8x,
     2 'sig1',8x,'sig2',3x,'angle'
     3 /38x,'s-11',8x,'s-12',8x,'s-22',8x,'s-33'/)
2002  format(2i9,2f12.4,2e12.4,f8.2/30x,4e12.4)
3000  format(' Input: e, nu, rho, thick, is (1=pl.stress,2=pl.strain)'
     1 /3x,'mate')
      end
c
      subroutine pian06(d,ul,s,p,nst,ndf,isw)
      implicit double precision (a-h,o-z)
c.... pian-sumihara stiffness matrix developed explicitly
c     last modified: 4/07/86
c
      dimension ul(ndf,*),s(nst,nst),p(*),d(*),rr(5),pl(8)
      common /elcom061/ xs,xt,xh,ys,yt,yh,xj0,xj1,xj2,r1,r2,
     1          x1(4),x2(4),y1(4),y2(4),vol
!$OMP THREADPRIVATE (/elcom061/)  
      common /elcom062/ ax(4),bx(4),cx(4),ay(4),by(4),cy(4),a1(3),a2(3),
     1          beta(6)
!$OMP THREADPRIVATE (/elcom062/)  
c..1.) set up stress interpolants for the 4-5 term
      r1 = xj1/xj0
      r2 = xj2/xj0
      a1(1) = xs*xs
      a1(2) = ys*ys
      a1(3) = xs*ys
      a2(1) = xt*xt
      a2(2) = yt*yt
      a2(3) = xt*yt
c..2.) set up shape function coefficients - jacobian weighted
      ax(1) = -yt + ys
      ax(2) =  yt + ys
      ax(3) = -ax(1)
      ax(4) = -ax(2)
      bx(1) = -yh - ys
      bx(2) = -bx(1)
      bx(3) =  yh - ys
      bx(4) = -bx(3)
      cx(1) =  yt + yh
      cx(2) = -yt + yh
      cx(3) = -cx(2)
      cx(4) = -cx(1)
      ay(1) =  xt - xs
      ay(2) = -xt - xs
      ay(3) = -ay(1)
      ay(4) = -ay(2)
      by(1) =  xh + xs
      by(2) = -by(1)
      by(3) = -xh + xs
      by(4) = -by(3)
      cy(1) = -xt - xh
      cy(2) =  xt - xh
      cy(3) = -cy(2)
      cy(4) = -cy(1)
c..3.) compute volume and stabilization h-array
      vol = 4.*xj0
      d11 = d(1)/vol
      d12 = d(2)/vol
      d33 = d(3)/vol
      hy  = vol*3.
      hx  = hy*d(5)
      h44 = hx*(1.-r2*r2/3.)*(a1(1)+a1(2))**2
      h55 = hx*(1.-r1*r1/3.)*(a2(1)+a2(2))**2
      h45 =-(r1*r2/3.)*(hx*(xs*xt+ys*yt)**2+d(6)*hy*(ys*xt-xs*yt)**2)
c..4.) Invert stabilization h-array
      hx  = h44*h55 - h45*h45
      hy  = h55/hx
      h55 = h44/hx
      h45 =-h45/hx
      h44 = hy
c..5.) Compute the current stress parameters
      call pzero(rr,5)
      do 50 j = 1,4
        hx = cx(j) - r2*ax(j)
        hy = cy(j) - r2*ay(j)
        x1(j) = a1(1)*hx + a1(3)*hy
        x2(j) = a1(2)*hy + a1(3)*hx
        hx = bx(j) - r1*ax(j)
        hy = by(j) - r1*ay(j)
        y1(j) = a2(1)*hx + a2(3)*hy
        y2(j) = a2(2)*hy + a2(3)*hx
        rr(1) = rr(1) + ax(j)*ul(1,j)
        rr(2) = rr(2) + ay(j)*ul(2,j)
        rr(3) = rr(3) + ay(j)*ul(1,j) + ax(j)*ul(2,j)
        rr(4) = rr(4) + x1(j)*ul(1,j) + x2(j)*ul(2,j)
        rr(5) = rr(5) + y1(j)*ul(1,j) + y2(j)*ul(2,j)
50    continue
      beta(1) = d11*rr(1) + d12*rr(2)
      beta(2) = d12*rr(1) + d11*rr(2)
      beta(3) = d33*rr(3)
      beta(4) = (h44*rr(4) + h45*rr(5))*3.
      beta(5) = (h45*rr(4) + h55*rr(5))*3.
c..6.) Form stiffness matrix for 1-pt (constant) terms
      if(isw.eq.3) then
        d11 = d11*d(11)
        d12 = d12*d(11)
        d33 = d33*d(11)
        i1 = 1
        do 110 i = 1,2
          bd11 = ax(i)*d11
          bd12 = ax(i)*d12
          bd13 = ay(i)*d33
          bd21 = ay(i)*d12
          bd22 = ay(i)*d11
          bd23 = ax(i)*d33
          j1   = i1
          do 100 j = i,2
            s(i1  ,j1  ) = bd11*ax(j) + bd13*ay(j)
            s(i1  ,j1+1) = bd12*ay(j) + bd13*ax(j)
            s(i1+1,j1  ) = bd21*ax(j) + bd23*ay(j)
            s(i1+1,j1+1) = bd22*ay(j) + bd23*ax(j)
            j1 = j1 + ndf
100       continue
          i1 = i1 + ndf
110     continue
c..7.)  Copy other parts from computed terms
        i1 = ndf + ndf
        do 120 i = 1,i1
          do 120 j = i,i1
            s(i ,j+i1) =-s(i,j)
            s(i+i1,j+i1) = s(i,j)
120     continue
        j1 = i1 + ndf
        s(2    ,i1+1) =-s(2,1)
        s(ndf+2,j1+1) =-s(ndf+2,ndf+1)
        s(ndf+1,i1+1) =-s(1,ndf+1)
        s(ndf+1,i1+2) =-s(2,ndf+1)
        s(ndf+2,i1+1) =-s(1,ndf+2)
        s(ndf+2,i1+2) =-s(2,ndf+2)
c..8.)  Add stabilization matrix
        h44 = h44*d(11)
        h45 = h45*d(11)
        h55 = h55*d(11)
        j1 = 1
        do 210 j = 1,4
          bd11 = h44*x1(j) + h45*y1(j)
          bd12 = h44*x2(j) + h45*y2(j)
          bd21 = h45*x1(j) + h55*y1(j)
          bd22 = h45*x2(j) + h55*y2(j)
          i1 = 1
          do 200 i = 1,j
            s(i1  ,j1  ) = s(i1  ,j1  ) + x1(i)*bd11 + y1(i)*bd21
            s(i1  ,j1+1) = s(i1  ,j1+1) + x1(i)*bd12 + y1(i)*bd22
            s(i1+1,j1  ) = s(i1+1,j1  ) + x2(i)*bd11 + y2(i)*bd21
            s(i1+1,j1+1) = s(i1+1,j1+1) + x2(i)*bd12 + y2(i)*bd22
            i1 = i1 + ndf
200      continue
         j1 = j1 + ndf
210     continue
      end if
c..9.) compute the residual force vector
      if(mod(isw,3).eq.0) then
c...a.) compute the constant part
        do 300 i = 1,5
          beta(i) = beta(i)*d(11)
  300   continue
        call pzero(pl,8)
        do 310 i = 1,2
          pl(2*i-1) = -(ax(i)*beta(1) + ay(i)*beta(3))
          pl(2*i+3) = -pl(2*i-1)
          pl(2*i  ) = -(ay(i)*beta(2) + ax(i)*beta(3))
          pl(2*i+4) = -pl(2*i  )
310     continue
c...b.) compute the stabilization part
        do 320 i = 1,4
          pl(2*i-1) = pl(2*i-1) - (x1(i)*beta(4) + y1(i)*beta(5))/3.0
          pl(2*i  ) = pl(2*i  ) - (x2(i)*beta(4) + y2(i)*beta(5))/3.0
320     continue
c...c.) copy into p vector (added ww 21.12.2006)      
        i1 = 1
        j1 = 1 
        do i = 1,4
          p(i1)   = pl(j1) 
          p(i1+1) = pl(j1+1) 
          i1 = i1+ndf
          j1 = j1+2
        end do
      end if
       
      return
      end
c
      subroutine stcn06(ix,d,ssa,tta,dt,st,nnp,ityp)
c ----------------------------------------------------------------------
c.... stress projection  and energy                                    |
c ----------------------------------------------------------------------
      USE errin1
      USE errin2
      implicit double precision (a-h,o-z)
      dimension   dt(nnp),st(nnp,*),ss(4),tt(4),ix(1),d(*),sig(4)
      common /elcom061/ xs,xt,xh,ys,yt,yh,xj0,xj1,xj2,r1,r2,
     1          x1(4),x2(4),y1(4),y2(4),vol
!$OMP THREADPRIVATE (/elcom061/)  
      common /elcom062/ ax(4),bx(4),cx(4),ay(4),by(4),cy(4),a1(3),a2(3),
     1          beta(6)
!$OMP THREADPRIVATE (/elcom062/)  
      data ss/-1.0,1.0,1.0,-1.0/,tt/-1.0,-1.0,1.0,1.0/
      gr = (1.+d(9))/d(8)
      sig=0
c.... compute stress projections
      do 200 jj = 1,4
        ll = abs(ix(jj))
        if(ll.gt.0) then
c....     compute weighted stresses at nodes
          xsj = xj0 + ss(jj)*xj1 + tt(jj)*xj2
          ssg = (ss(jj) - ssa)*beta(5)
          ttg = (tt(jj) - tta)*beta(4)
c....     accumulate jacobian weights in each element
          dt(ll)   = dt(ll)   + xsj
c....     11-stress contribution
          sig(1)   = (beta(1) + a1(1)*ttg + a2(1)*ssg) 
          st(ll,1) = st(ll,1) + sig(1)*xsj  
c....     22-stress contribution
          sig(3)   = (beta(2) + a1(2)*ttg + a2(2)*ssg) 
          st(ll,3) = st(ll,3) + sig(3)*xsj
c....     12-stress contribution
          sig(2)   = (beta(3) + a1(3)*ttg + a2(3)*ssg)
          st(ll,2) = st(ll,2) + sig(2)*xsj
c....     33-stress contribution (non-zero for plane strain only)
          sig(4)   = d(13)*(sig(1)+sig(3))
          st(ll,4) = st(ll,4) + sig(4)*xsj

          s4=sig(4)
          if(ityp.eq.1) s4=0.d0 ! plane stress
          u_om(1)=u_om(1)
     +         +(d(14)*(sig(1)+sig(3)+s4)**2+gr*sig(2)**2)*xsj
c....     element energy
          ntyp = 4
          if(ityp.eq.1) ntyp = 3
          do 110 i = 1,ntyp
            u_om(1) = u_om(1) + gr*sig(i)**2*xsj
            u_om(2) = u_om(2) +    sig(i)**2*xsj
110     continue

        end if
200   continue
      return
      end
c
      subroutine ster06(ix,d,xl,ssa,tta,dt,st,nnp,ityp,e_ome,numel,ndm)
c ----------------------------------------------------------------------
c.... error calculation                                                |
c ----------------------------------------------------------------------
      USE errin1
      USE errin2
      implicit double precision (a-h,o-z)
      dimension   dt(nnp),st(nnp,*),ss(4),tt(4),ix(*),d(*),sig(4)
     +         ,e_ome(numel,*),e_ome06(numerr),sigp(4),dsig(4),xl(ndm,*)
      common /elcom061/ xs,xt,xh,ys,yt,yh,xj0,xj1,xj2,r1,r2,
     1          x1(4),x2(4),y1(4),y2(4),vol
!$OMP THREADPRIVATE (/elcom061/)  
      common /elcom062/ ax(4),bx(4),cx(4),ay(4),by(4),cy(4),a1(3),a2(3),
     1          beta(6)
!$OMP THREADPRIVATE (/elcom062/)  
      data ss/-1.0,1.0,1.0,-1.0/,tt/-1.0,-1.0,1.0,1.0/
c.... intial values for element errors
      e_ome06 = 0.0d0

      gr = (1.+d(9))/d(8)
      sig=0
      sigp=0
      dsig=0
      do 200 jj = 1,4
c....   compute stresses at Gausspoints=nodes
        xsj = xj0 + ss(jj)*xj1 + tt(jj)*xj2
        ssg = (ss(jj) - ssa)*beta(5)
        ttg = (tt(jj) - tta)*beta(4)
        sig(1)  = (beta(1) + a1(1)*ttg + a2(1)*ssg) 
        sig(3)  = (beta(2) + a1(2)*ttg + a2(2)*ssg) 
        sig(2)  = (beta(3) + a1(3)*ttg + a2(3)*ssg)
        sig(4)  = d(13)*(sig(1)+sig(3))
c....   stresses at nodes from stress projection
        ll = abs(ix(jj))
        if(ll.gt.0) then
          sigp(1) = st(ll,1)    
          sigp(3) = st(ll,3)  
          sigp(2) = st(ll,2)  
          sigp(4) = st(ll,4)  
        end if 
c.....  stress differences
        do i = 1,4
          dsig(i) = sigp(i)- sig(i)
        end do 
c.....  element errors
        ntyp=4
        if(ityp.eq.1) ntyp=3
        do i = 1,ntyp
          e_ome06(1)= e_ome06(1) + gr*(dsig(i)**2)*xsj
          e_ome06(2)= e_ome06(2) +    (dsig(i)**2)*xsj
        end do 
c....   additional term
        s4=dsig(4)
        if(ityp.eq.1) s4=0.d0
        e_ome06(1)= e_ome06(1) 
     +         +(d(14)*(dsig(1)+dsig(3)+s4)**2+gr*dsig(2)**2)*xsj
200   continue
c.....plot/print/add errors
      call elmterr(ix,xl,ndm,numel,e_ome06,e_ome)
      return
      end
