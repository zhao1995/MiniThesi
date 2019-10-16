      subroutine elmt01(d,ul,xl,ix,tl,s,p,h1,h2,h3,ndf,ndm,nst,isw)
c---------------------------------------------------------------+
c     any dimensional truss element-geometrically linear 
c
c     Input E,A,rho=gamma/g
c
c     Material
c     d(1)=E
c     d(2)=A
c     d(3)=rho
c     d(4)=E*A
c     d(5)=rho*A
c
c     Options
c      3=TANG
c      4=STRE
c      5=LMAS
c     13=FORC forc,1 = N_x   forc,7 = eps_x  forc,8 = sigma_x 
c
c     Comments
c     WW KIT 11/2013: Plot,forc ist nur 2D ausgelegt, 3D geht,sofern  
c                     nichtr alle x-y Komponenten verschwinden 
c
c---------------------------------------------------------------+
      USE bdata
      USE cdata
      USE eldata
      USE fornam
      USE iofile
      USE pdata6
      USE pdata10
      USE pltran
      USE strnam
      implicit real*8 (a-h,o-z)
      integer*4 ix(*),ixl(5)
      real*8 xl(ndm,*),tl(*),xlp(3,4),sig(3),d(*),
     1       ul(ndf,*),vl(ndf*nel),s(nst,nst),p(nst),db(3),dx(3),xx(3)
      dimension h1(*),h2(*),h3(*)
c
      data ixl /1,2,3,4,1/,epsi/1.0e-5/
      ielno = 1
      nm = min(ndm,ndf)
      go to (1,2,3,4,5,4,2,2,2,2,2,2,4),isw
      return
c.... input material properties
1     if(ior.lt.0) write(*,3000)
      call dinput(d,3)
      d(4) = d(1)*d(2)
      d(5) = d(3)*d(2)
                   write(iow,2000) d(1),d(2),d(3)
      if(ior.lt.0) write(*  ,2000) d(1),d(2),d(3)
c.... define node numbering for plot mesh routine, see pltord
        inord(ielno)   = 4
        ipord(1,ielno) = 1
        ipord(2,ielno) = 2
        ipord(3,ielno) = 2
        ipord(4,ielno) = 1
c.... description of stresses  
      do i = 1,11         
         forsus(i) =  ' '
      end do
      forsus( 1) =  '  N-FORCE N_xx '
      forsus( 7) =  '  STRAIN  E_xx '
      forsus( 8) =  '  STRESS  S_xx '
2     return
c.... compute element stiffness matrix
3     xll = 0.0
      do 31 i = 1,nm
        dx(i) = xl(i,2) - xl(i,1)
        db(i) = d(4)*dx(i)
        xll = xll + dx(i)**2
31    continue
      xll = xll*sqrt(xll)
      do 32 i = 1,nm
        dx(i) = dx(i)/xll
32    continue
      i1 = 0
      do 36 ii = 1,2
        j1 = i1
        do 35 jj = ii,2
          do 33 i = 1,nm
          do 33 j = 1,nm
            s(i+i1,j+j1) = db(i)*dx(j)
33        continue
          j1 = j1 + ndf
          do 34 j = 1,nm
            dx(j) = -dx(j)
34        continue
35      continue
        i1 = i1 + ndf
36    continue
      do 37 i = 1,nm
      do 37 j = 1,nm
        s(i+ndf,j) = s(j,i+ndf)
37    continue
c.... form a residual
      call pmove(ul,vl,ndf*nel)
      do 38 i = 1,nst
      do 38 j = 1,nst
        p(i) = p(i) - s(i,j)*vl(j)
38    continue
c.... calculate eigenvalues
cww   call elemev(s,nel,ndf,nst)
      return
c.... output stress and strain in element
4     xll = 0.0
      eps = 0.0
      call pzero(xx,3)
      do 41 i = 1,nm
        dx(i) = xl(i,2) - xl(i,1)
        xll = xll + dx(i)**2
        eps = eps + dx(i)*(ul(i,2)-ul(i,1))
        xx(i) = (xl(i,2) + xl(i,1))/2.
41    continue
      eps    = eps/xll
      sig(1) = d(4)*eps
      sig(2) = eps
      sig(3) = d(1)*eps
c.... plot N,eps,sig on screen
      if(isw.eq.13) then
        klayf = 1
        ns = 0
        if(nfp.eq.1) ns=1
        if(nfp.eq.7) ns=2
        if(nfp.eq.8) ns=3
        if(ns.eq.0) return
        if(flfp) then
          ccfp  = sig(ns)
          xmaxf = max(xmaxf,ccfp)
          xminf = min(xminf,ccfp)
          cfp  = max(cfp,abs(ccfp))
        else
          if(abs(sig(ns)).lt.epsi) return
          call pzero(xlp,12)
c....     plot on mesh/deformed mesh
cww       scal01 = scal  ! plot on deformed mesh
          scal01 = 0.d0  ! plot on undeformed mesh      
          sn = (xl(2,2)+scal01*ul(2,2)) - (xl(2,1)+scal01*ul(2,1))
          cs = (xl(1,2)+scal01*ul(1,2)) - (xl(1,1)+scal01*ul(1,1))
          sl = dsqrt(cs*cs + sn*sn)
          if(sl.eq.0.d0) return ! provisorisch ww
	        sn = sn/sl
	        cs = cs/sl
          xlp(1,1) = xl(1,1) + scal01*ul(1,1)
          xlp(2,1) = xl(2,1) + scal01*ul(2,1)
          xlp(1,2) = xl(1,2) + scal01*ul(1,2)
          xlp(2,2) = xl(2,2) + scal01*ul(2,2)
          xlp(1,3) = xl(1,2) + scal01*ul(1,2) - sn*sig(ns)*cfp
          xlp(2,3) = xl(2,2) + scal01*ul(2,2) + cs*sig(ns)*cfp
          xlp(1,4) = xl(1,1) + scal01*ul(1,1) - sn*sig(ns)*cfp
          xlp(2,4) = xl(2,1) + scal01*ul(2,1) + cs*sig(ns)*cfp
          if(ndm.eq.3) then
            xlp(3,1) = xl(3,1) + scal01*ul(3,1)
            xlp(3,2) = xl(3,2) + scal01*ul(3,2)
            xlp(3,3) = xl(3,2) + scal01*ul(3,2)
            xlp(3,4) = xl(3,1) + scal01*ul(3,1)
          end if
c......   find color
          signs =sig(ns)
          call pppcolf(signs)
c          
          call plxtrn(xlp,tra,vr,ndm,4)
          call plot9s(ixl,xlp,3,4)
        end if
        return
      end if
      if(isw.eq.4) then
        mct = mct - 1
        if(mct.le.0) then
          write(iow,2001) o,head
          if(ior.lt.0) write(*,2001) o,head
          mct = 50
        end if
                     write(iow,2002) n,ma,xx,(sig(i),i=1,3)
        if(ior.lt.0) write(*  ,2002) n,ma,xx,(sig(i),i=1,3)
      else
c.... compute right hand side
        sign = sig(1)/sqrt(xll)
        do 42 i = 1,nm
          p(i) = dx(i)*sign
          p(i+ndf) = -p(i)
42      continue
      end if
      return
c.... compute element lumped mass matrix
5     xll = 0.0
      do 51 i = 1,nm
        xll = xll + (xl(i,2)-xl(i,1))**2
51    continue
      sm = d(5)*sqrt(xll)/2.0
      do 52 i = 1,nm
        p(i    ) = sm
        p(i+ndf) = sm
52    continue
      return
c.... formats
2000  format(5x,'T r u s s   E l e m e n t'//10x,'modulus =',e13.5/10x,
     1   'area    =',e13.5/10x,'density =',e13.5)
2001  format(a1,20a4,//5x,'T r u s s    E l e m e n t',//,' elem mate',
     1   3x,'1-coord',3x,'2-coord',3x,'3-coord',5x,'force',7x,'strain',
     2   7x,'stress')
2002  format(2i5,3f10.4,3e13.5)
3000  format(' Input: E, A, rho'/3x,'>',\)
      end

