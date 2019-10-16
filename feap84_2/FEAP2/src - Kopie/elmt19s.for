      subroutine elmt19(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c---------------------------------------------------------------------
c.... dkt-triangular plate bending element      Elple2_2
c.... El-Nummer für Plotten!
c     Eingabe  E,nu,h,q
c
c.... grober test am 4.12.91/1.3.93 ww
c.... Pruefe:geht stress projection besser, wenn S an Knoten berechnet?
c.... update sowie isw=13(str1) 1.10.2003 ww 
c---------------------------------------------------------------
c     Ref: Batoz, Bathe, Ho, "A Study of 3-Node Triangular Plate
c          Bending Elements," IJNME, v 15, no. 12, pp 1771-1812,
c          1980
c---------------------------------------------------------------------
      USE bdata
      USE cdata
      USE eldata
      USE iofile
      USE pdata6
      USE pdata7
      USE pdata10
      USE pltran
      USE prisdat
      USE prlod
      USE qload      
      USE strnam 
      implicit double precision (a-h,o-z)
      common /elm219/b(3),c(3),bd(3),cd(3),aa(3),bb(3),cc(3),dd(3),ee(3)
!$OMP THREADPRIVATE (/elm219/)  
      common /shp219/   xsj,bm(3,9)
!$OMP THREADPRIVATE (/shp219/)  
      dimension d(*),xl(ndm,*),ul(*),ix(*),tl(*),s(nst,*),p(nst),
     1 el(3,3),eps(3),sig(3),ipordl(4),yl(3,3),sigp(7)
      dimension h1(*),h2(*),h3(*)
      data ipordl /1,2,3,1/
c.... 4-times area coordinates stored in 'el' to reduce multiplications
      data el/2.0d0,0.0d0,2.0d0,2.0d0,2.0d0,0.0d0,0.0d0,2.0d0,2.0d0/
c------------------|
      ielno = 19    
c------------------|
c.... transfer to correct processor
      go to (1,2,3,4,5,3,2,4,2,2,2,2,3,2,2,2,2,2,2,2,2,22), isw
c.... input the material properties: E,nue,h,q
1     call dinput(d(4),4)
      thick = d(6)
      d(1) = d(4)/(1.-d(5)*d(5))*thick**3/12.0
      d(2) = d(5)*d(1)
      d(3) = 0.5*(d(1)-d(2))
      write(iow,2000) d(4),d(5),thick,d(7)
      if(ndm.eq.3) ipla = 1
c.... description of stresses  
      strsus( 1) = '  MOMENT M_xx  '
      strsus( 2) = '  MOMENT M_xy  '
      strsus( 3) = '  MOMENT M_yy  '
      strsus( 4) = '               '
      strsus( 5) = '  MOMENT M_11  '
      strsus( 6) = '  MOMENT M_22  '
      strsus( 7) = '  ANGLE Phi_1  '
      do is =8,25
        strsus(is) = '               '
      end do
c...  names for principal moments
      nptyp = 3 
c
c.... define node numbering for plot mesh routine, see pltord
      inord(ielno) = 4
      do 11 ii = 1,4
11    ipord(ii,ielno) = ipordl(ii)
c.... check element for errors
2     return
c.... compute the element tangent array
3     call jactri19(xl,ndm)
      xsj = xsj/6.0d0
c.... compute weighted jacobian and d-matrix constants
      d1 = d(1)*xsj
      d2 = d(2)*xsj
      d3 = d(3)*xsj
c.... compute the element load vector
      call qload19(qz,d,aqloa,numel,n,mqloa,propq,prop,isw)
      p(1) = qz*xsj
      p(4) = p(1)
      p(7) = p(1)
      do 310 l = 1,3
        call dktbm19(el(1,l))
c....   compute contribution to element stiffness for this point
        do 300 i = 1,9
          dn1 = d1*bm(1,i) + d2*bm(2,i)
          dn2 = d2*bm(1,i) + d1*bm(2,i)
          dn3 = d3*bm(3,i)
          do 300 j = i,9
300     s(i,j) = s(i,j) + dn1*bm(1,j) + dn2*bm(2,j) + dn3*bm(3,j)
310   continue
c.... make stiffness symmetric
      do 320 i = 2,9
        i1 = i - 1
        do 320 j = 1,i1
320   s(i,j) = s(j,i)
c      if(n.eq.1) call mprint(s,nst,nst,nst,' stiffness ')
c.... compute the element residual vector: R = P - K*v
      do 350 j = 1,nst
	      uu = ul(j)
	      if(uu.ne.0.0d0) then
	      do 340 i = 1,nst
	        p(i) = p(i) - s(i,j)*uu
340	    continue
	    end if
350   continue
      return
c.... compute and output the element variables
4     call jactri19(xl,ndm)
      xsj = xsj/6.0d0
      eps(1) = 4./3.d0
      eps(2) = eps(1)
      eps(3) = eps(1)
      call dktbm19(eps)
      do 410 i = 1,3
        eps(i) = 0.0
        do 410 j = 1,9
410   eps(i) = eps(i) + bm(i,j)*ul(j)
      sig(1) = -(d(1)*eps(1) + d(2)*eps(2)) ! -sign due to dkq
      sig(2) = -(d(2)*eps(1) + d(1)*eps(2))
      sig(3) = -(d(3)*eps(3)/2.0          )
      if(isw.eq.4) then
c....   print stresses
        xx = 1./3.*(xl(1,1)+xl(1,2)+xl(1,3))
        yy = 1./3.*(xl(2,1)+xl(2,2)+xl(2,3))
        mct = mct - 1
        if(mct.gt.0) go to 420
        mct = 50
                     write(iow,2002) o,head
        if(ior.lt.0) write(*  ,2002) o,head
420                  write(iow,2003) n,ma,xx,yy,sig
        if(ior.lt.0) write(*  ,2003) n,ma,xx,yy,sig
      else if(isw.eq.8) then
c....   plot  stresses
        istv= 7
        if(iplma(ma).eq.0)  return ! only if MATN
        call stcn19(ix,strea,strea(1+numnp),sig,xsj,numnp)
      else if(isw.eq.13) then
c.....  plot stress resultants without averaging 
        sigp(1)=sig(1)
        sigp(2)=sig(3)
        sigp(3)=sig(2)
        call pstres(sigp,sigp(5),sigp(6),sigp(7))
        if(nfp.gt.7.or.nfp.lt.1.or.nfp.eq.4) return
        strp = sigp(nfp) 
        if(flfp) then
c....     calculate extreme values
          xmaxf = max(xmaxf,strp)
          xminf = min(xminf,strp)
        else
c....     plot stresses 
c....     color      
          call pppcolf(strp) 
c....     coordinates 3d
          call pzero(yl,9)
          do i=1,ndm
            do j=1,3
              yl(i,j)=xl(i,j) 
            end do
          end do
c.....    transform ccordinates for rot        
          call plxtrn(yl,tra,vr,3,3)
c......   plot element
          call plot9(iel,ix,yl,3,nel,1)
        end if
      end if
      return
c.... compute element mass arrays
5     continue
      return
c.... compute the element residual vector
c6    computed under isw=3
c
c.... compute the load vector
22    call jactri19(xl,ndm)
      xsj = xsj/6.0d0
c.... compute the element load vector
      call qload19(qz,d,aqloa,numel,n,mqloa,propq,prop,isw)
      p(1) = qz*xsj
      p(4) = p(1)
      p(7) = p(1)
      return
c.... formats
2000  format(5x,'dkt triangular plate element'/9x,'material constants'
     1 //10x,12hmodulus    = ,e15.5/10x,12hpoisson-nu = ,e15.5/
     2 10x,12hthickness  = ,e15.5/10x,12hq-loading  = ,e15.5/)
2002  format(a1,20a4//5x,'dkt triangular plate moments'//
     1 1x,'elmt',2x,'mat',5x,'1-coord',5x,'2-coord',
     2 5x,'11-moment',5x,'22-moment',5x,'12-moment'/)
2003  format(2i5,2f12.3,3e14.5)
      end
c
      subroutine jactri19(xl,ndm)
c----------------------------------------------------------------------
c     calculate Jacobian for DKT-element
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension xl(ndm,*)
      common/elm219/b(3),c(3),bd(3),cd(3),aa(3),bb(3),cc(3),dd(3),ee(3)
!$OMP THREADPRIVATE (/elm219/)  
      common/shp219/xsj,bm(3,9)
!$OMP THREADPRIVATE (/shp219/)  
      do 100 i = 1,3
      j = mod(i,3) + 1
      k = mod(j,3) + 1
      b(i) = xl(2,k) - xl(2,j)
      c(i) = xl(1,j) - xl(1,k)
      sql = (b(i)+b(i))*(b(i)+b(i)) + (c(i)+c(i))*(c(i)+c(i))
      cs = c(i)/sql
      bs = b(i)/sql
      aa(i) = 6.*cs
      bb(i) = 3.*bs*c(i)
      cc(i) = c(i)*cs - b(i)*(bs+bs)
      dd(i) =-6.*bs
100   ee(i) = b(i)*bs - c(i)*(cs+cs)
      xsj = xl(1,1)*b(1) + xl(1,2)*b(2) + xl(1,3)*b(3)
      xsj = - xsj
      do 110 i = 1,3
      bd(i) = b(i)/xsj
110   cd(i) = c(i)/xsj
      return
      end
c
      subroutine dktbm19(el)
c----------------------------------------------------------------------
c     calculate B-Matrix for DKT-element
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension el(3),shn(2,3),
     1 a1(3),a2(3),b1(3),b2(3),c1(3),c2(3),d1(3),d2(3),e1(3),e2(3)
      common/elm219/b(3),c(3),bd(3),cd(3),aa(3),bb(3),cc(3),dd(3),ee(3)
!$OMP THREADPRIVATE (/elm219/)  
      common/shp219/xsj,bm(3,9)
!$OMP THREADPRIVATE (/shp219/)  
      do 100 i = 1,3
      j = mod(i,3) + 1
      k = mod(j,3) + 1
      shn(1,i) = bd(i)*(el(i) - 1.0)
      shn(2,i) = cd(i)*(el(i) - 1.0)
      shm1 = bd(j)*el(k) + bd(k)*el(j)
      shm2 = cd(j)*el(k) + cd(k)*el(j)
      a1(i) = aa(i)*shm1
      a2(i) = aa(i)*shm2
      b1(i) = bb(i)*shm1
      b2(i) = bb(i)*shm2
      c1(i) = cc(i)*shm1
      c2(i) = cc(i)*shm2
      d1(i) = dd(i)*shm1
      d2(i) = dd(i)*shm2
      e1(i) = ee(i)*shm1
100   e2(i) = ee(i)*shm2
      i1 = 1
      do 110 i = 1,3
      j = mod(i,3) + 1
      k = mod(j,3) + 1
      i2 = i1 + 1
      i3 = i2 + 1
      bm(1,i1) = a1(k) - a1(j)
      bm(1,i2) = b1(k) + b1(j)
      bm(1,i3) = c1(k) + c1(j) - shn(1,i)
      bm(2,i1) = d2(k) - d2(j)
      bm(2,i2) =-e2(k) - e2(j) + shn(2,i)
      bm(2,i3) =-b2(k) - b2(j)
      bm(3,i1) = a2(k) - a2(j) + d1(k) - d1(j)
      bm(3,i2) =-e1(k) - e1(j) + shn(1,i) - bm(2,i3)
      bm(3,i3) = c2(k) + c2(j) - shn(2,i) - bm(1,i2)
110   i1 = i1 + 3
      return
      end
c
      subroutine stcn19(ix,dt,st,sig,xsj,numnp)
c----------------------------------------------------------------------
c     plot stresses  Mx(1) Mxy(2) My(3) [M_1(5) M_2(6) phi_1(7)] 
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension ix(*),st(numnp,*),dt(numnp),sig(3)
        do 10 ll = 1,3
          ii = abs(ix(ll))
          if(ii.le.0) go to 10
          dt(ii)   = dt(ii)   +        xsj
          st(ii,1) = st(ii,1) + sig(1)*xsj
          st(ii,2) = st(ii,2) + sig(3)*xsj
          st(ii,3) = st(ii,3) + sig(2)*xsj
 10     continue
      return
      end
      subroutine qload19(qz,d,q,numel,n,mqloa,propq,prop,isw)
c----------------------------------------------------------
c.... set loads from macro qloa/mate
c----------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension q(numel,10),d(*)
      if(isw.eq.22) then
        qz = 0.d0 
        if(mqloa.ne.1) qz = q(n,1)*propq 
      else    
        qz = d(7)*prop 
      end if 
      return
      end
      
