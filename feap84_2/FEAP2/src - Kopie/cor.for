      subroutine cor(al,au,ad,x,c,xmas,jdiag,ne,nzykel,omega,epsa)
      USE iofile
      implicit double precision(a-h,o-z)
      dimension al(*),au(*),ad(*),c(*),x(*),jdiag(*),xmas(*)
      neq = abs(ne)
c.....parameter
      eps    = 1.0d-10
      npoint = 0
c.....initialisierung:
      xnorm = sqrt(ddot(neq,x,1,x,1))
      do 10 i = 1,neq
        x(i) = x(i)/xnorm
10    continue
      call pzero(c,neq)
      call promul(al,au,ad,x,c,jdiag,ne,1,1)
      a = ddot(neq,c,1,x,1)
      b = 0.d0
      do 100 i = 1,neq
        b = b + x(i)*x(i)*xmas(i)
100   continue
      if (b.eq.0.0d0)
     1  write(*,*)' Nenner ist 0 in subroutine COR'
      r=a/b
      write(*,*) ' Initial eigenvalue in COR:', r
c.....allgem. zyklus
      do 1100 k=1,nzykel
        jd = 0
        vpsi = 0.d0
        do 1000 j=1,neq
          if (k.eq.1.and.j.eq.1) then
            f = c(1)
            g = x(1)*xmas(1)
          else
c.........multiplikation von x mit der j-ten Zeile von A bzw. B
            call pzero(c,neq)
            call proml(al,au,ad,x,c,jdiag,jd,ne,.true.,j)
            f = c(j)
            g = x(j)*xmas(j)
          end if
          p = ad(j)
          q = xmas(j)
c.........Loesung des kleinen EWP
          x2 = b*q - g*g
          x1 =-b*p - a*q + 2*g*f
          x0 = a*p - f*f
c.........hier ist  det(t) = x2*t^2+x1*t+x0
          if (abs(x2).lt.1.d-20) then
            if(npoint.ne.1)write(*,*)' x2 ist 0 im subroutine cor '
            npoint = 1
          endif
          xp =-x1/x2/2.0d0
          xq = x0/x2
c.........hier ist t^2 - 2 * xp * t  + xq
          dis = xp*xp - xq
          if (dis .lt. 0.0d0)
     1       write(*,*) ' Dis ist kleiner 0 in subroutine cor '
          rs = xp - sqrt(dis)
c.........rs ist kleinster EW des kleinen EWP
c.........Fallunterscheidung:
          if( abs(p - rs*q).gt.eps )  then
            psi  =-omega*(f-rs*g)/(p-rs*q)
            x(j) = x(j)+psi
            if(abs(psi).gt.vpsi.and.vpsi.ge.0.d0) vpsi = abs(psi)
            a = a + 2*psi*f + psi*psi*p
            b = b + 2*psi*g + psi*psi*q
            r =a/b
          else
            if (abs(a-b*rs).gt.eps) then
              call pzero(x,neq)
              x(j) = 1.0d0
              a    = p
              b    = q
              r    = a/b
              vpsi = -1.d0
            endif
c...........letzter fall benoetigt keine Aenderung
          endif
1000    continue
c       write(*,*) '   Vpsi aus cor:', vpsi
        if(vpsi.ge.0.d0) then
        xnorm = sqrt(ddot(neq,x,1,x,1))
        vpsi = vpsi/xnorm
        if(vpsi.le.epsa) then
          if(ior.lt.0) write(*,3000) k,r,vpsi
          write(iow,3000) k,r,vpsi
          return
        endif
      endif
1100  continue
      if(ior.lt.0) write(*,*)'  Max. number of cycles UP COR exeeded '
      if(ior.lt.0) write(*,3000) k,r,vpsi
      write(iow,*)'  Max. number of cycles UP COR exeeded '
      write(iow,3000) k,r,vpsi
      return
3000  format(' Iteration:',i3,' Eigenvalue:',g12.5,' vpsi:',g12.5)
      end
c---------------------------------------------------------------------
      subroutine initei(x,neq)
      implicit double precision(a-h,o-z)
      dimension x(*)
      do 100 i=1,neq
        x(i) = 1.d0
100   continue
      return
      end
c------------------------------------------------------------------
      subroutine proml(al,au,ad,b,c,jdiag,jd,ne,add,j)
      implicit double precision (a-h,o-z)
      logical add
      dimension al(*),au(*),ad(*),b(*),c(*),jdiag(*)
c
c.... routine to form c_i = c_i +/- a_il*b_l where a is a square matrix
c.... stored in profile form, b is a vector, c_i a component of
c     a vector and jdiag locates the bottom of columns or rows in a.
c.... if add .true. : add to c, if add .false. : subtract from c
c
c.... multiply lower triangular part
      neq = abs(ne)
        js = jd
        jd = jdiag(j)
        if(jd.gt.js) then
          jh = jd - js
c         ab = dot(al(js+1),b(j-jh),jh)
          ab = ddot(jh,al(js+1),1,b(j-jh),1)
          if(add) then
            c(j) = c(j) + ab
          else
            c(j) = c(j) - ab
          endif
        endif
c.... do diagonal part
        if(add) then
          c(j) = c(j) + ad(j)*b(j)
        else
          c(j) = c(j) - ad(j)*b(j)
        endif
c.... multiply the upper triangular part
      jdd = jdiag(j)
      do 200 k = j+1,neq
        js  = jdd
        jdd = jdiag(k)
        jh  = jdd - js
        if(jh.ge.k-j) then
          bj  = b(k)
          if(.not.add) bj = -bj
          index = jdd-k+j+1
          if(au(index).ne.0.d0) c(j) = c(j) + bj*au(index)
        endif
200   continue
      return
      end
c
      subroutine pvmv(a,b,n,x,add)
      implicit double precision(a-h,o-z)
      dimension a(*), b(*), x(*)
      logical add
c
c     multiply:  x_i = x_i + a_i * b_i
c
      if(add) then
        do 100 i=1,n
          x(i) = x(i) + a(i)*b(i)
100     continue
      else
        do 200 i=1,n
          x(i) = x(i) - a(i)*b(i)
200     continue
      endif
      return
      end
