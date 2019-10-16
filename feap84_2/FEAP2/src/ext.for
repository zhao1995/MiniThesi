      subroutine ext(du,h,dh,c,d,z1,z2,e,
     +         f,f0,al,au,ad,u,jdiag,id,nneq,rlnew)
c----------------------------------------------------------------------
c     algorithm for calculation of limit and symmetrical bifurcation
c            -    points via extended systems
c     almost Newton's method: using numerical derivatives of Kt
c
c     parameters
c     h  (neq)  - eigenvector associated with singular point
c     u  (nneq) - displacement and increments
c     dh (neq)  - displacement increment du3
c     du (neq)  - displacement increment du2
c     z2 (neq)  - displacement increment du1
c     z1 (nneq) - perturbed displacement field
c     d  (neq)  - load p   ==> h1 ==> p1
c     c  (neq)  - residual g   ==> h2 ==> p2
c     e  (neq)  - vector epsei ==> h3 ==> p3
c         commons
c     /ext1/   in pform, pmacr1 and pmacr3
c     /ext2/   in ext, pform, ploa1, pmacr and pmacr3
c     /arcext/ in datri,pform,pmacr,pmacr1,pmacr3
c
c         p. wriggers,  uc berkeley  June 4, 1988
c----------------------------------------------------------------------
      USE arcext
      USE cdata
      USE ext2
      USE iofile
      USE prlod
      implicit double precision(a-h,o-z)
      logical fa,tr
      dimension du(*),h(*),dh(*),c(*),d(*),z1(*),z2(*),f0(*),
     1    f(*),al(*),au(*),ad(*),u(*),jdiag(*),id(*),e(*)
      data fa,tr/.false.,.true./
c
c.... rescale eigenvector to norm=1
      fact = sqrt(ddot(neq,h,1,h,1))
      do 100 i = 1,neq
        h(i) = h(i) / fact
  100 continue
c.......calculate z2 for load level 1.0, z2 = z2! + ... fuer link
      call pzero(z2,neq) ! added Knebel
      do 10 i =1,nneq
        j = id(i)
        if (j.gt.0) z2(j) = z2(j) + f(i)*prop + f0(i)
10    continue
      call ploads(u,z2,prop,.false.,.false.,.false.)
      call pmove(z2,d,neq)
c.... solve for du1 = z2
      call dasol(al,au,ad,z2,jdiag,neq,energy)
c.... in case of singularity
      if(kex.ne.0) then
        call pzero(dh,neq)
        call pzero(e,neq)
        dh(kex) = exeps
        e(kex)  = exeps
        call dasol(al,au,ad,dh,jdiag,neq,energy)
      end if
c.... compute perturbed displacements (z1 = u + epsc*h)
      call peps(u,h,nneq,neq,epsc,z1,id)
c.... form tangent stiffness at z1 and comp. num. derivatives
      nmode = 1
cww   hflgu/h3flgu fuer history werte setzen?
      call formfe(z1,du,z2,dh,tr,tr,fa,fa,11,1,numel,1)
      nmode = 0
c.... compute c and d (c = h2, d = h1)
      do 15 i = 1,neq
        c(i) = + c(i)/epsc
        d(i) = + d(i)/epsc
15    continue
      if(kex.ne.0) then
        do 16 i = 1,neq
          e(i) = e(i)/epsc
16      continue
        c(kex) = c(kex) - du(kex)*exeps/epsc
        d(kex) = d(kex) - z2(kex)*exeps/epsc
        e(kex) = e(kex) - dh(kex)*exeps/epsc
      end if
c.... solve with right hand sides c,d (c = p2, d = p1)
      call dasol (al,au,ad,c,jdiag,neq,energy)
      call dasol (al,au,ad,d,jdiag,neq,energy)
      if(kex.ne.0) call dasol (al,au,ad,e,jdiag,neq,energy)
c.... Compute load increment, displacement increments and h
      if(kex.eq.0) then
c....   scalar product c * h and d * h
        cdh = ddot(neq,h,1,c,1)
        ddh = ddot(neq,h,1,d,1)
        hn  = sqrt(ddot(neq,h,1,h,1))
c....   calculate load increment
        dlamb = -(cdh - hn) / ddh
c....   calculate displacement and eigenvector increment
        do 20 i = 1,neq
          h(i)  = c(i)  + dlamb * d(i)
          du(i) = du(i) + dlamb * z2(i)
20      continue
      else
        do 30 i = 1,nneq
          if(id(i).eq.kex) then
            eiu = u(i)
            go to 31
          endif
30      continue
31      continue
        a11 = d(kex)
        a12 = e(kex)
        a21 = z2(kex)
        a22 = dh(kex) - 1.d0
        b1 = 1.d0 - c(kex) - dh(kex) - (xmu-eiu)*e(kex)
        b2 = xmu - eiu - du(kex) - (xmu-eiu)*dh(kex)
        deta = a11*a22 - a12*a21
        ahe  =   a22/deta
        a22  =   a11/deta
        a11  =   ahe
        a12  = - a12/deta
        a21  = - a21/deta
        dlamb = a11*b1 + a12*b2
        dmu   = a21*b1 + a22*b2
        do 32 i = 1,neq
          du(i) = du(i) + dlamb*z2(i) + (xmu+dmu-eiu)*dh(i)
          h(i)  = c(i)  + dlamb*d(i)  + (xmu+dmu-eiu)*e(i) + dh(i)
32      continue
        xmu = xmu + dmu
      end if
c.... update load level, displacements and eigenvector
      rlnew = rlnew + dlamb
c.... check whether limit or bifurcation point appears
      bitest = dotid(f,h,id,nneq)
      hn = sqrt(ddot(neq,h,1,h,1))
      bitest = bitest/hn
      if(abs(bitest).lt.1.d-2) then
        if(ior.lt.0) write(*  ,2004) bitest
                     write(iow,2004) bitest
      else
        if(ior.lt.0) write(*  ,2005) bitest
                     write(iow,2005) bitest
      end if
c.... print current state of iteration
      if(kex.eq.0) then
        hn  = hn - 1.0d0
      else
        hn  = h(kex) - 1.0d0
      end if
             write(iow,2000) rlnew*prop,dlamb*prop,xmu,hn
      if(ior.lt.0) write(*  ,2000) rlnew*prop,dlamb*prop,xmu,hn
c.... calculate determinant of K_T for tplo_det
      call detkt(ad,neq,0)
      return
c....  formats
2000  format(5x,'total load factor............ =',e15.6,/,
     1       5x,'load increment............... =',e15.6,/,
     2       5x,'xmu.......................... =',e15.6,/,
     3       5x,'norm of constraint:   |h|-1   =',e15.6)
2004  format(5x,'BIFURCATION point, ht*f       =',e15.6)
2005  format(5x,'LIMIT point,       ht*f       =',e15.6)
      end
c
      subroutine pmove2(al,au,ad,a,b,jdiag,neq,isw)
c----------------------------------------------------------------------
c.... initial vector eigenvector
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension al(*),au(*),ad(*),a(*),b(*),jdiag(*)
c.... move a-array into b-array
      go to (1,2,3,4),isw
c.... initial vector for h = (1,1,1....,1,1) /norm
1     continue
      do 11 n = 1,neq
 11     b(n) = 1.
      xn  = sqrt(ddot(neq,b,1,b,1))
      do 12 n = 1,neq
 12     b(n) = b(n) / xn
      return
c.... initial vector for h = du / norm(du)
2     continue
      xn  = sqrt(ddot(neq,a,1,a,1))
      do 21 n = 1,neq
 21      b(n) = a(n) / xn
      return
c.... initial vector for h = Kt-1 * (1,1,,,1)/norm
3     continue
      do 31 n = 1,neq
 31     b(n) = 1.
      xn  = sqrt(ddot(neq,b,1,b,1))
      do 32 n = 1,neq
 32     b(n) = b(n) / xn
      call dasol (al,au,ad,b,jdiag,neq,aengy)
      xn  = sqrt(ddot(neq,b,1,b,1))
      do 33 n = 1,neq
 33     b(n) = b(n) / xn
      return
4     continue
      do 41 n = 1,neq
 41     b(n) = 1.
      xn  = sqrt(ddot(neq,b,1,b,1))
      do 42 n = 1,neq
 42     b(n) = b(n) / xn
      do 44 ii = 1,10
        call dasol (al,au,ad,b,jdiag,neq,aengy)
        xn  = sqrt(ddot(neq,b,1,b,1))
      do 43 n = 1,neq
 43     b(n) = b(n) / xn
 44   continue
      return
      end
c
      subroutine writev1(a,nn,name)
      USE iofile
      implicit double precision(a-h,o-z)
      dimension a(*)
      character*4 name
c.... write array a
      write(iow,2000) name,nn
      write(iow,2001) (a(i),i=1,nn)
2000  format(5x,'vector ',a4,' lenght:',i4)
2001  format(5x,6(g10.4))
      return
      end
c
      subroutine peps(u,phi,nneq,neq,epsc,ueps,id)
      USE ext2
      implicit double precision(a-h,o-z)
      dimension u(*),phi(*),ueps(*),id(*)
      umax = 0.0d0
      phin = sqrt(ddot(neq,phi,1,phi,1))
      write(*,2000) phin
      do 10 i = 1,nneq
cww     if(abs(u(i)).gt.umax) umax = u(i) !falsch bei angl?
        if(abs(u(i)).gt.umax) umax = abs(u(i))
10    continue
      epsc = eps*umax/phin
cww   if (epsc.eq.0.d0) epsc = 1.0d-7 !falsch da epsc z.B. 1.d-20
      epsc = max(epsc,1.0d-7)
      do 20 i = 1,nneq
        ueps(i) = u(i)
        j = id(i)
        if (j .gt. 0) ueps(i)      = u(i) + epsc*phi(j)
        if (j .gt. 0) ueps(i+2*nneq) = epsc*phi(j)
20    continue
      return
2000  format(5x,'phi - norm :',g15.5)
      end
c
      subroutine pmult(s,p,nst,u,isw)
      implicit double precision (a-h,o-z)
      dimension s(nst,nst),p(nst),u(nst)
c
      call pzero(p,nst)
c
      if(isw.eq.1) then
        do 100 i = 1,nst
          do 100 j = 1,nst
            p(i) = p(i) + s(i,j) * u(j)
100     continue
      else if(isw.eq.2) then
        do 110 i = 1,nst
          do 110 j = 1,nst
            p(i) = p(i) - s(i,j) * u(j)
110     continue
      end if
      return
      end
c
      subroutine scaleh(phi,neq,kex)
c----------------------------------------------------------------------
c...  scale eigenvector to Phi(kex) = 1
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension phi(*)
      if(abs(phi(kex)) .lt. 1.d-8) then
        write(*,2000) kex,phi(kex)
       return
      end if
      fact = 1.d0/phi(kex)
      do 100 i = 1,neq
        phi(i) = phi(i)*fact
100   continue
      return
2000  format(' WARNING: use another initialization',/
     1       'kex = ',i4,'  phi(kex) =',g15.5)
      end
