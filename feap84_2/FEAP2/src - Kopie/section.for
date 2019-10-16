      subroutine secthp(ix,u,dr,ctf,ct4,cs)
c-----------------------------------------------------------------------
c
c     Purpose:
c     main program for sectional values
c
c     Inputs:
c      ix(nen1,numel)  - nodes on elements
c      u(ndf,numnp)    - displacements
c      dr              - dummy
c      ctf(3),ct4      - ct(1-4)
c      cs              - scaling factor
c
c-----------------------------------------------------------------------
c.... adresses where values are  [length of field]
c     thinset
c     msec( 1,ma)   [  2*numels(1,ma)  ]      ix, 2-node mesh
c     msec( 2,ma)   [  2*numnps(1,ma)  ]      xy, 2-node mesh
c     msec( 3,ma)   [  1*numnps(1,ma)  ]       w, 2-node mesh
c     msec( 4,ma)   [  4*numels(1,ma)  ]     x+w, scale plot
c     msec( 5,ma)   [ 12*numels(1,ma)  ]     dat, input for each segment
c     thicset
c     msec( 6,ma)   [ 4*numels(2,ma)   ]     ixm, 4-node mesh
c     msec( 7,ma)   [ 2*numnps(2,ma)   ]     xym, 4-node mesh
c     msec( 8,ma)   [ 1*numnps(2,ma)   ]      wm, 4-node mesh
c     msec( 9,ma)   [ (npstrs+1)*numnps(2,ma)] npstrs+1 stresses
c     msec(10,ma)   [  2*3  ]                 centre g,s,r
c     msec(11,ma)   [  3*3  ]                 principal axis
c     msec(12,ma)   [    5  ]                 Iy,Iz,It,Cm,alpha
c     msec(13,ma)   [ nOrth ]                 see SR thicsec2
c-----------------------------------------------------------------------
c
c     Comments:
c
c     restrictions
c     ma=10 materials see common /sectio/
c     npstrs = 26
c
c     generation of arrays is general!
c
c-----------------------------------------------------------------------
      USE cdata
      USE psize
      USE sdata
      USE sectio
      implicit double precision (a-h,o-z)
      dimension ix(nen1,*),u(ndf,*),dr(*),ctf(3),ct(4)
      common m(maxm)

c.... GP for output(default)
      igps = 0

      do i=1,3
        ct(i) = ctf(i)
      end do
      ct(4) = ct4
      ma    = ct(3)
      if(ma.le.0) ma = 1
c.... for sect,12/13 = stre/flux
      if(ct(1).eq.12.or.ct(1).eq.13) then
        ma1 = ix(nen1,ma)       ! ix and mate-nr. of beam-element
        do i = 1,nen
	        ii = ix(i,ma)
	        if(ii.gt.0) nel = i
        end do
        ma = ma1
c....   GP for output
        if(ct(1).eq.12) igps = ct(4)
        if(ct(1).eq.13) igps = ct(3)
      end if

c.... plot values
      call section(u,dr,ct,m(msec(1,ma)),m(msec(2,ma)),m(msec(3,ma)),
     +      m(msec(4,ma)),m(msec(6,ma)),m(msec(7,ma)),m(msec(8,ma)),
     +      m(msec(9,ma)),m(msec(10,ma)),m(msec(11,ma)),m(msec(12,ma)),
     +      npstrs,nel,2,cs,
     +      numnps(1,ma),numels(1,ma),numnps(2,ma),numels(2,ma))

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine section(u,dr,ct,ix,x,w,xw,ixm,xm,wm,str,cgsr,vec,sec,
     +              npstrs,nel,nen,cs,numnp1,numel1,numnp2,numel2)
c-----------------------------------------------------------------------
c
c     Purpose: macro program
c-----------------------------------------------------------------------
c     u,dr,ct(4) as standard
c-----------------------------------------------------------------------
c     ix(2,*)       id-field          line    numel1
c     x (2,*)       coordinates       line    numnp1
c     w (*)         warping function  line
c     ixm(4,*)      id-field          mesh    numel2
c     xm (2,*)      coordinates       mesh    numnp2
c     wm (*)        warping function  mesh
c     str(npstrs,*) stresses          mesh
c     cgsr(2,3)  (1,1) center of gravity
c     cgsr(2,3)  (1,2) center of shear
c     cgsr(2,3)  (1,3) reference point
c     vec(3,3)    base vectors columnwise  values 0-1
c     sec(9)      section values
c-----------------------------------------------------------------------
c.... [sect,n1,n2,n3]
c                 |
c                 +> section (material) number
c     [sect,1]       mesh                            on line element
c     [sect,2,n2]    nodes (n2.ne.0-> with numbers)  on line element
c     [sect,3]       element numbers                 on line element
c     [sect,4]       center of gravity
c     [sect,5]       center of shear
c     [sect,6]       reference point
c     [sect,7,n2]    function w (n2 -> scale w)      on line element
c     [sect,8]       mesh                            on area element
c     [sect,9,n2]    nodes (n2.ne.0-> with numbers)  on area element
c     [sect,10]      element numbers                 on area element
c     [sect,11,n2]   function w (n2 -> mc lines <=0 filled)
c     [sect,12,n2,n3,n4] stresses n2 at element n3, gauss-point n4
c     [sect,13,n2,n3] flux at element n2, gauss-point n3, lenght n4
c----------------------------------------------------------------------
      USE pdata1
      USE pdata2
      implicit double precision (a-h,o-z)
      dimension u(*),dr(*),ct(*),ix(nen,*),x(2,*),w(*),xw(2,*),
     +          ixm(4,*),xm(2,*),wm(*),str(numnp2,*),
     +          vec(3,*),xvec(3,3),cgsr(2,*),sec(*)


      data xvec /1,0,0, 0,1,0, 0,0,1/
c
      isw1 = ct(1)
      isw2 = ct(2)
      isw3 = ct(3)
      isw4 = ct(4)
      rsw4 = ct(4)
      if(rsw4.eq.0.d0) rsw4=1.d0
c.... scale plot default
      if(numnp1.gt.0) then    ! for x+w on line mesh
         call plotse0(ix,x,w,xw,nen,numnp1,numel1,1.0d0)
      else                    ! for x   on full mesh
        call frame(xm,2,numnp2,1,.false.)
      end if

      scale = scale*cs
c.... open Graphic device
      if(isw1.ne.11.and.isw1.ne.12) call plopen

c.... action
      goto(1,2,3,4,5,6,7,8,9,10,11,12,13) isw1

c.... [sect,1] mesh on line element
1     call plotse1(ix,x,nen,numel1)
      return

c.... [sect,2,n2] nodes on line element (n2.ne.0-> with numbers)
2     call plotse2(x,numnp1,scale,isw2)
      return

c.... [sect,3] element numbers on line element
3     call plotse3(ix,x,nen,numel1,scale)
      return

c.... [sect,4] center of gravity
4     call plotse4(cgsr(1,1), vec,sec,scale,'GravityCent',isw1)
      return

c.... [sect,5] center of shear
5     call plotse4(cgsr(1,2),xvec,sec,scale,'Shear Cent.',isw1)
      return

c.... [sect,6] reference point
6     call plotse4(cgsr(1,3),xvec,sec,scale,'ReferencePt',isw1)
      return

c.... [sect,7] function w on line element
7     if(numnp1.gt.0)
     + call plotse7(ix,x,w,nen,numnp1,numel1,scale,ct(2))
      return

c.... [sect,8] mesh on area element
8     call plotse1(ixm,xm,4,numel2)
      return

c.... [sect,9,n2] nodes on area element (n2.ne.0-> with numbers)
9     call plotse2(xm,numnp2,scale,isw2)
      return

c.... [sect,10] element numbers on area element
10    call plotse3(ixm,xm,4,numel2,scale)
      return

c.... [sect,11,n2] funktion w on area element (n2 ->mc lines <=0 filled)
11    call plotse11(xm,ixm,wm,ndf2,numnp2,numel2,ct,idev)
      return

c.... [sect,12,n2,n3,n4] stresses n2 at element n3, gauss-point n4
12    call plotse12(xm,ixm,str,ct,u,dr,numnp2,numel2,nel,idev,
     +              isw2,isw3,isw4)
      return

c.... [sect,13,n2,n3] flux at element n2, gauss-point n3, length csw4
13    call plotse13(xm,str,ct,u,dr,numnp2,numel2,nel,idev,
     +              isw2,isw3,rsw4)

      return
      end
c
      subroutine plotse0(ix,x,w,dr,nen,numnp1,numel1,fact1)
c----------------------------------------------------------------------
c     Purpose: scale problem
c
c     Inputs:
c       ix(nen,*) id-field
c       x (2,*)   nodal coordinates
c       dr(2,*)   x + w  on nodes/element
c       w (*)     function
c
c----------------------------------------------------------------------
      USE pdata1
      implicit double precision (a-h,o-z)
      dimension ix(nen,*),x(2,*),w(*),dr(2,*)
c.... min / max  x+w
      if (dabs(fact1).lt.1e-10) fact1 = 1.d0
      wm = 0.d0
      do i = 1,numnp1
        wm = max(wm,dabs(w(i)))
      end do
      call frame(x,2,2*numel1,1,.false.)
      fact = 0.1d0 / scale / wm * fact1
c.... length of vectors
      call pzero(dr,2*2*numel1)
      i=0
    	do n = 1,numel1
c.....  nodes
	      i1      = ix(1,n)
	      i2      = ix(2,n)
c.....  direction
        dxq     = x(1,i2) - x(1,i1)
        dyq     = x(2,i2) - x(2,i1)
        ds = sqrt(dxq*dxq+dyq*dyq)
        sn = dyq/ds
        cs = dxq/ds
c....   coordinates of rectangle to plot
c....   point 1
        i=i+1
	      dr(1,i) = x(1,i1) - w(i1)*sn*fact
	      dr(2,i) = x(2,i1) + w(i1)*cs*fact
c....   point 2
        i=i+1
        dr(1,i) = x(1,i2) - w(i2)*sn*fact
        dr(2,i) = x(2,i2) + w(i2)*cs*fact
      end do
      call frame(dr,2,2*numel1,1,.false.)
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine plotse1(ix,x,nen,numel1)
c----------------------------------------------------------------------
c     Purpose: plot mesh of section   on line element
c
c     Inputs:
c       ix(nen,*) id-field
c       x (2,*)   nodal coordinates
c       xl(2,4)   plot  coordinates of element
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension ix(nen,*),x(2,*),xl(2,4)
c...  loop over elements to draw mesh elementwise
      call pzero(xl,8)
      call pppcol(4)
	    do n = 1,numel1
c....   coordinates of element
	      do i = 1,nen
	        ii = ix(i,n)
	        xl(1,i) = x(1,ii)
	        xl(2,i) = x(2,ii)
        end do
c.....  plot element
        call plotse1a(xl,nen)
     	end do
      return
      end
c
      subroutine plotse1a(xl,nen)
c-----------------------------------------------------------------------
c.... Purpose: plot elements of line element
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension xl(2,*)
      call plotl(xl(1,1),xl(2,1),0.d0,3)
      do i = 2,nen
	      call plotl(xl(1,i),xl(2,i),0.d0,2)
      end do
      call plotl(xl(1,1),xl(2,1),0.d0,2)
      return
      end
c
      subroutine plotse2(x,numnp1,scale,n1)
c-----------------------------------------------------------------------
c.... Purpose:  plot nodes of line element
c
c     Inputs
c       x (2,numnp1)   nodal coordinates
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension x(2,*)
      dx1 = 0.002/scale
      x3  = 0.d0
      call pppcol(7)
      do n = 1,numnp1
        x1 = x(1,n)
        x2 = x(2,n)
        call plotl(x1-dx1, x2+dx1, x3, 3)
        call plotl(x1-dx1, x2+dx1, x3, 1)
        call plotl(x1-dx1, x2-dx1, x3, 2)
        call plotl(x1+dx1, x2-dx1, x3, 2)
        call plotl(x1+dx1, x2+dx1, x3, 2)
        call plotl(x1-dx1, x2+dx1, x3, 2)
        call clpan
        call plotl(x1-12.d0*dx1, x2+2.d0*dx1, x3, 3)
        if(n1.ne.0) then
          call pppcol(1) ! plot number in black
          call plabl(n)
          call pppcol(7)
        end if
      end do
      return
      end
c
      subroutine plotse3(ix,x,nen,numel1,scale)
c-----------------------------------------------------------------------
c.... Purpose:  plot element numbers of line element
c
c     Inputs
c       ix(nen,numel1) id-field
c       x (2,*)        nodal coordinates
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension ix(nen,*),x(2,*)
      call pppcol(3)
c.... modify position
      dx1 = .005/scale
      do n = 1,numel1
c....   element center
	      x1 = 0.d0
	      x2 = 0.d0
	      do i = 1,nen
	        ii = ix(i,n)
	        x1 = x1 + x(1,ii)
	        x2 = x2 + x(2,ii)
        end do
	      x1 = x1/nen
	      x2 = x2/nen
c....   position of label
        x1 = x1 - dx1*5.d0
	      x2 = x2 - dx1
	      call plotl(x1,x2,0.d0,3)
c....   write element number
	      call plabl(n)
      end do
      return
      end
c
      subroutine plotse4(x,xvec,sec,scale,cc,isw)
c-----------------------------------------------------------------------
c.... purpose: plot one node of section
c     G=gravity,M=shear,R=reference
c
c-----------------------------------------------------------------------
      USE pdata4
      implicit double precision (a-h,o-z)
      dimension x(2,1),xpos(3),dum(3),xvec(3,3),sec(*)
      character*11 cc
      call pzero(dum,3)
      call pzero(xpos,3)
      dx1 = 0.002/scale
      sc1 = 0.1*max(xmax(1)-xmin(1),xmax(2)-xmin(2),xmax(3)-xmin(3))
      x3  = 0.d0
      call pppcol(2)
    	x1 = x(1,1)
	    x2 = x(2,1)
	    call plotl(x1-dx1, x2+dx1, x3, 3)
	    call plotl(x1-dx1, x2+dx1, x3, 1)
	    call plotl(x1-dx1, x2-dx1, x3, 2)
	    call plotl(x1+dx1, x2-dx1, x3, 2)
	    call plotl(x1+dx1, x2+dx1, x3, 2)
	    call plotl(x1-dx1, x2+dx1, x3, 2)
	    call clpan
      call plotl(x1-15.d0*dx1, x2-5.d0*dx1, x3, 3)
      call plablc(cc)
c.... plot coordinate system
      xpos(1) = x1
      xpos(2) = x2
      if (isw.ne.5) call pltaxs1(xvec,dum,xpos,2,sc1)
c.... plot values of center
      if (isw.ne.6) call pltwarp(x1,x2,sec,isw)
      return
      end
c
      subroutine plablc(yy)
c-----------------------------------------------------------------------
c     plot label on screen                                             |
c-----------------------------------------------------------------------
      USE hpgl1
      USE pdata2
      USE pdatap
      USE pftn77
      USE plotter
      implicit double precision (a-h,o-z)
      character yy*11
      if(idev.eq.1)  call gptx2(xp(3),11,yy)
      if(idev.eq.2)  call gtx(xpp(2),ypp(2),yy)
      if(idev.eq.3)  call draw_text(yy,ixa,iya,icc)
      if(idev.eq.4)  call draw_text(yy,ixa,iya,icc)
      if(nexte.gt.0.and.iprin.eq.1) then
        call hptext (xxp(2),yyp(2),yy,ihpgl)
      end if
      return
      end
c
      subroutine plotse7(ix,x,w,nen,numnp1,numel1,scale,fact1)
c----------------------------------------------------------------------
c     Purpose: plot warping function on line elements
c
c     Inputs:
c       ix(nen,*) id-field
c       x (2,*)   nodal coordinates
c       w (*)     function
c
c     Scratch
c       xl(2,4)
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension ix(nen,*),x(2,*),w(*),xl(2,4)
      if(nen.gt.2) return
c.... length of vectors
      if (dabs(fact1).lt.1e-10) fact1 = 1.d0
      wm = 0.d0
      do i = 1,numnp1
        wm = max(wm,dabs(w(i)))
      end do
      fact = 0.1d0 / scale / wm * fact1
c...  loop over elements to draw function elementwise
      call pzero(xl,8)
	    do 100 n = 1,numel1
c.....  nodes
	      i1 = ix(1,n)
	      i2 = ix(2,n)
c.....  direction
        dx = x(1,i2) - x(1,i1)
        dy = x(2,i2) - x(2,i1)
        ds = sqrt(dx*dx+dy*dy)
        sn = dy/ds
        cs = dx/ds
c....   coordinates of rectangle to plot
c....   point 1
        xl(1,1) = x(1,i1)
        xl(2,1) = x(2,i1)
c....   point 4
        xl(1,4) = xl(1,1) - w(i1)*sn*fact
        xl(2,4) = xl(2,1) + w(i1)*cs*fact
c....   point 2
        xl(1,2) = x(1,i2)
        xl(2,2) = x(2,i2)
c....   point 3
        xl(1,3) = xl(1,2) - w(i2)*sn*fact
        xl(2,3) = xl(2,2) + w(i2)*cs*fact
c...    Berechne Vorzeichen -> farbe + =>blau - =>rot
        w1 = w(i1)
        w2 = w(i2)
        if(w1.gt.0.d0.and.w2.gt.0.d0) then
	        call pppcol(3)
	        goto 110
        else if(w1.lt.0.d0.and.w2.lt.0.d0) then
          call pppcol(2)
	        goto 110
        else
c....     Vorzeichenwechsel
          if(dabs(w1).lt.1.d-10.and.dabs(w2).lt.1.d-10) goto 100
          x1 = 0.d0
          x2 = ds
          y1 = w1
          y2 = w2
c.....    Gerade y = y1 + a1*x   mit a1 = (y2-y1)/dl
          a1  = (y2-y1)/ds
c.....    Durchstosspunkt
          xdp = -y1/a1
c....     Plotte dreieck 1
          if(w1.gt.0.d0) call pppcol(3)
          if(w1.lt.0.d0) call pppcol(2)
c....     point 2
          xl(1,2) = x(1,i1)*(1.d0-xdp/ds) + x(1,i2)*xdp/ds
          xl(2,2) = x(2,i1)*(1.d0-xdp/ds) + x(2,i2)*xdp/ds
c....     point 3
          xl(1,3) = xl(1,2)
          xl(2,3) = xl(2,2)
c.....    plot function over element
          call plotse7a(xl)
c...      plotte dreieck 2
          if(w2.gt.0.d0) call pppcol(3)
          if(w2.lt.0.d0) call pppcol(2)
c....     point 1
          xl(1,1) = xl(1,2)
          xl(2,1) = xl(2,2)
c....     point 4
          xl(1,4) = xl(1,1)
          xl(2,4) = xl(2,1)
c....     point 2
	        xl(1,2) = x(1,i2)
	        xl(2,2) = x(2,i2)
c....     point 3
          xl(1,3) = xl(1,2) - w(i2)*sn*fact
          xl(2,3) = xl(2,2) + w(i2)*cs*fact
        end if
c....   plot function over element
110     call plotse7a(xl)
100   continue
      return
      end
c
      subroutine plotse7a(xl)
c-----------------------------------------------------------------------
c.... Purpose: plot function  over element
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension xl(2,*)
      call plotl(xl(1,1),xl(2,1),0.d0,3)
      call plotl(xl(1,1),xl(2,1),0.d0,1)
      do i = 2,4
	      call plotl(xl(1,i),xl(2,i),0.d0,2)
      end do
      call plotl(xl(1,1),xl(2,1),0.d0,2)
      call clpan
      return
      end
c
c
      subroutine pltaxs1(tr,vr,xi,ndm,ct)
c----------------------------------------------------------------------
c.... Purpose draw vectors for axes, including rot-macro
c----------------------------------------------------------------------
      USE pdata2
      implicit double precision (a-h,o-z)
      dimension tr(3,*),vr(*),xi(*),xx(3,5)
c.... compute plot location for axes and transform xi
c.... perspective projecton of axes
      do 120 m = 1,ndm
        call pzero(xx,15)
        do 100 n = 1,ndm
          fac1 = tr(m,n)*ct
          xx(n,1) = xi(n)
          xx(n,2) = xx(n,1) + fac1
          xx(n,5) = xx(n,2)
100     continue
        fac1 = tr(m,1)*ct
        fac2 = tr(m,2)*ct
        fac3 = tr(m,3)*ct
        xx(1,3) = xx(1,2) -.3*fac1 - .1*(fac2+fac3)
        xx(2,3) = xx(2,2) -.3*fac2 + .1*(fac1+fac3)
        xx(3,3) = xx(3,2) -.3*fac3 + .1*(fac1+fac2)
        xx(1,4) = xx(1,2) -.3*fac1 + .1*(fac2+fac3)
        xx(2,4) = xx(2,2) -.3*fac2 - .1*(fac1+fac3)
        xx(3,4) = xx(3,2) -.3*fac3 - .1*(fac1+fac2)
c....   plot the vector
        call plotl(xx(1,1),xx(2,1),xx(3,1),3)
        call plotl(xx(1,2),xx(2,2),xx(3,2),2)
        call plotl(xx(1,2),xx(2,2),xx(3,2),ipgl)
        call plotl(xx(1,3),xx(2,3),xx(3,3),2)
        call plotl(xx(1,4),xx(2,4),xx(3,4),2)
        call plotl(xx(1,2),xx(2,2),xx(3,2),2)
        if(ipgl.eq.1) call clpan
        call plotl(xx(1,2),xx(2,2),xx(3,2),3)
        call plabl(m+1)
120   continue
      return
      end
c
      subroutine pltwarp(ys,zs,sec,isw)
c-----------------------------------------------------------------------
c.... plot section properties on legend                                |
c-----------------------------------------------------------------------
      USE pdata2
      implicit double precision(a-h,o-z)
      character yy*15,c*2,ai*2
      dimension sec(*)
      x = 1.
      y = 0.76
      call drawtxt(1,x,y,1,1,15,'Distance from  ')
      y = 0.73
      call drawtxt(1,x,y,1,1,15,'Reference Point')
      if(isw.eq.4) c='S'
      if(isw.eq.5) c='M'
      write(yy,'(a1,a4,e10.4)') c,'_2: ',ys
      if(idev.eq.4) x = 0.97
      y = 0.64
      if(isw.eq.5) y = 0.55
      call drawtxt(1,x,y,1,1,15,yy)
      write(yy,'(a1,a4,e10.4)') c,'_3: ',zs
      y = 0.61
      if(isw.eq.5) y = 0.52
      call drawtxt(1,x,y,1,1,15,yy)
c
      if(isw.eq.4) then
      y=0.42
      call drawtxt(1,x,y,1,1,15,'Mome. of Inerta')
      y=0.38
      ai='I'
      write(yy,'(a1,a4,e10.4)') ai,'_2: ',sec(1)
      call drawtxt(1,x,y,1,1,15,yy)
      y=0.35
      ai='I'
      write(yy,'(a1,a4,e10.4)') ai,'_3: ',sec(2)
      call drawtxt(1,x,y,1,1,15,yy)
      y=0.30
      ai='a'
      write(yy,'(a1,a4,e10.4)') ai,'_0: ',sec(5)
      call drawtxt(1,x,y,1,1,15,yy)
      end if

      if(isw.eq.5) then
      y=0.22
      call drawtxt(1,x,y,1,1,15,'Torsional Prop.')
      y=0.18
      ai='I'
      write(yy,'(a1,a4,e10.4)') ai,'_t: ',sec(3)
      call drawtxt(1,x,y,1,1,15,yy)
      y=0.15
      ai='C'
      write(yy,'(a1,a4,e10.4)') ai,'_m: ',sec(4)
      call drawtxt(1,x,y,1,1,15,yy)
      end if
      return
      end
c
      subroutine plotse11(xm,ixm,wm,ndf2,numnp2,numel2,ct,idev)
c-----------------------------------------------------------------------
c.... Purpose: plot warping function on 4 node elements
c
c     Inputs:
c       ixm(4,*)   id-field
c       xm (2,*)   nodal coordinates
c       wm (*)     function
c
c     Comments:
c     open
cww   strsus(1) = '  WARPING WT   '
c     strsus(1) = '  WARPING WQ   '
c-----------------------------------------------------------------------
      USE strnam
      implicit double precision (a-h,o-z)
      character*15 strsus1
      dimension  xm(2,*),ixm(4,*),wm(*),ct(4),dummy(1)
c.....name of stress to plot
      strsus1   = strsus(1) ! save
cww   strsus(1) = '  WARPING WT   '
      strsus(1) = '  WARPING WQ   '
      ct(1)= 1    ! only one value possible
      ct(2)= -1   ! > 0 ->mc lines  <= 0 filled
      ct(3)= 0    ! contour plot
      ndf2 = 1    ! plot field is wm(ndf1,*)
      cinv = 1    ! set color range red -> blue
      ct4  = 0    !
      call rprint1(wm,ixm,xm,2,numnp2,ndf2,1,4,ndf2,idev)
      call pltcons(xm,ixm,wm,2,ndf2,numnp2,numel2,4,cinv,ct,ct4)
      strsus(1) = strsus1 ! reset
      return
      end
c
      subroutine plotse12(xm,ixm,str,ct,u,dr,numnp2,numel2,nel,idev,
     +                    n2,n3,n4)
c-----------------------------------------------------------------------
c.... Purpose: plot stresses on 4 node elements
c
c     Inputs:
c       ixm(4,*)   id-field
c       xm (2,*)   nodal coordinates
c       str (*)    stresses
c       n2         stress number
c       n3         element
c       n4         Gauss-Point
c
c-----------------------------------------------------------------------
      USE cdata
      USE hdatam
      USE sectio
      implicit double precision (a-h,o-z)
      logical fa
      dimension xm(2,*),ixm(4,*),str(numnp2,*),ct(*),u(*),dr(*),dummy(1)
      fa = .false.

c.... check values
      n2 = max(1,min(npstrs,n2))    ! stresses for plot
      n3 = max(1,min(numel,n3))     ! beam elements
      n4 = max(1,min(nel-1,n4))     ! gauss-points in beam element
      igps=n4                       ! gp to sectio.h
c
      ct(1)= n2      ! values for plot
      ct(2)= -n2     ! > 0 ->mc lines  <= 0 filled
      ct(3)= 0       ! contour plot
      ndf1 = 1       ! min/max value in plot
      cinv = 1       ! set color range red -> blue
      ct4  = 0       ! plot additionally mesh in color ct4
      nn  = 1
      istv = -npstrs
	    call pzero (str,(npstrs+1)*numnp2)
      hflgu  = .false.
      h3flgu = .false.
      call formfe(u,dr,dr,dr,fa,fa,fa,fa,20,n3,n3,nn)
      call pltstr1(str(1,1),str(1,2),numnp2,npstrs)
      call rprint1(str(1,n2+1),ixm,xm,2,numnp2,ndf1,1,4,ndf1,idev)
      call pltcons(xm,ixm,str(1,n2+1),2,ndf1,numnp2,numel2,
     +              4,cinv,ct,ct4)

      do i=1,numnp2
        write(16,*) i,str(i,5),str(i,6)
      end do

      return
      end
c
      subroutine plotse13(xm,str,ct,u,dr,numnp2,numel2,nel,idev,
     +                    n2,n3,r3)
c-----------------------------------------------------------------------
c.... Purpose: plot flux on 4 node elements
c
c     Inputs:
c       ixm(4,*)   id-field
c       xm (2,*)   nodal coordinates
c       str (*)    stresses
c       n2         element
c       n3         Gauss-Point
c       r3         length scaling for flux
c
c     Comments
c       Flux-vectors must be stored in str(i,5) and Str(i,6)!!
c       -> sig(4),sig(5)!  str(i,1)=dt!!
c
c-----------------------------------------------------------------------
      USE cdata
      USE hdatam
      USE pltran
      USE sectio
      implicit double precision (a-h,o-z)
      logical fa
      dimension xm(2,*),str(numnp2,*),ct(*),u(*),dr(*)
      fa = .false.

c.... check values
      n2 = max(1,min(numel,n2))       ! selected beam element
      n3 = max(1,min(nel,n3))         ! selected gauss-point
      igps=n3                         ! gp to sectio.h

c...  calculate
      istv = -npstrs
	    call pzero (str,(npstrs+1)*numnp2)
      hflgu  = .false.
      h3flgu = .false.
      call formfe(u,dr,dr,dr,fa,fa,fa,fa,20,n2,n2,1)
      call pltstr1(str(1,1),str(1,2),numnp2,npstrs)
      call pltflu(xm,str(1,5),2,numnp2,tra,vr,r3)

      return
      end

c
      subroutine pltcons(x,id,b,ndm,ndf,numnp1,numel1,nen1,cinv,ct,ct4)
c-----------------------------------------------------------------------
c
c.... Purpose: Plot contour values for displacements, stresses etc
c
c     Input:
c      x(ndm,numnp1)   - nodal coordinates
c      id(ndf,numnp1)  - equation numbers for each active dof
c      b(*)           - array to plot
c      ndm            - spatial dimension of mesh
c      ndf            - number dof/node
c      numnp1         - number of elements to plot
c      numel1         - number of elements to plot
c      nen1           - nen+4
c      cinv          - =  1 colors from -=red ->+=blue set by colo,,cinv
c                      = -1 colors from -=blue->+=red
c      ic [ct(1)]     - number of eg. stress to plot, e.g. stress_1
c      mc [ct(2)]     > 0 plot lines   (variable contour 1-14)
c      ipb[ct(3)]    .ne.0 contour plot without numbers
c                     < 0 filled plots (variable contour 1-13)
c      ct4            - plot additionally mesh in color ct4
c
c     Output:
c
c     Comments:
c       To DO: anpassen an SR PLTCON!!
c
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE contval
      USE errchk
      USE iofile
      USE pdata1
      USE pdata2
      USE plotter
      USE rpdata
      implicit double precision (a-h,o-z)
      character*1 y
      logical zoom,tvc(9,9),cont
      integer*4 ipal(14),ipali(14)
      real*8  xl(3,40),x(ndm,*),b(*),v(40),vc(14),vt(4),xt(3,4)
      real*8  ct(3)
      integer id(nen1,*),iplt(5),ilc(4)
      save cont,vl,vu,nc,nnc,vc,xmx,ymx,xmn,ymn,vmn,vmx
c.... plot nodes of element (fixed for 4)
      data iu /5/,iplt/1,2,3,4,1/
c----------------------------------------------------------------------
c
      ic    = ct(1)
      mc    = ct(2)
      ipb   = ct(3)
      k4    = ct4
      mmc   = 1
      iadd1 = 16
      call pzerol ( tvc , .true. , 81 )
      dx1   = .024/scale
c
1     continue
15    continue
      if(icv.eq.1) then ! interactive
        if(ior.lt.0) write(*,2003)
        call dinput(contv,3)
      end if
      nci   = contv(1)
      nc    = min(14,nci)
      if(nc.le.0) nc = 14
      nnc = nc
      if(mc.le.0) nc = nc-1
      vc(1) = contv(2)
      vc(2) = contv(3)
      if(errck) go to 15
      drm =  dabs(rmx-rmn)
      if(drm.lt.1.e-5*dabs(rmx).or.drm.lt.1.e-10) then
c....   no lines for nearly same values or zero values
        nnc = 1
        nc  = 0
        vc(1) = rmn
        vc(2) = rmx
      else if(vc(1).eq.vc(2)) then
c....	  no input(default)
        vc(1) = rmn
        vc(2) = rmx
      else
c....   input (value from screen)
        if(nc.gt.1) vc(nc) = vc(2)
      end if
c
c.... generate intermediate values
      if(nc.ge.1) then !
c....   limit values
        vmin = rmn
        vmax = rmx
        vc0=vc(1)
        vcn=vc(2)
        if(vc0.lt.vmin) vmin = vc0
        if(vcn.gt.vmax) vmax = vcn
c....   intermediate values: nc colors->nc-1 points, nc lines->nc points
        do i = 1,nc
          vc(i) = vc0 + (vcn-vc0)*i/(nc+1)
        end do
      end if
      if(mc.gt.0) cont = .true.
      if(mc.le.0) cont = .false.
      if (icv.eq.1) then ! interactive
        if(ior.ge.0) write(iow,2000) (vc(i),i=1,nc)
        if(ior.lt.0) write(*  ,2000) (vc(i),i=1,nc)
c....   if interactive offer chance to change inputs
        if(ior.lt.0) then
	        write(*,2002)
900	      read (*,1000,err=901,end=902) y
	        goto  903
901	      call  errclr ('PLTCONS')
	        goto  900
902	      call  endclr ('PLTCONS',y)
903	      if(y .eq. 'c' .or. y .eq. 'C') then
            if(idev.eq.4) call clwclose(1,2,0)
	          return
	        else if(y .ne. 'y' .and. y .ne. 'Y') then
	          go to 1
	        end if
          if(idev.eq.4) call clwclose(1,2,0)
        end if
      end if
c.... find max/min of plot variable
      j   = ic
      xmx = x(1,1)
      ymx = x(2,1)
      xmn = x(1,1)
      ymn = x(2,1)
      vmn = b(j)
      vmx = b(j)
      do 17 i = 1,numnp1
   	    xmx = max(x(1,i),xmx)
  	    ymx = max(x(2,i),ymx)
	      xmn = min(x(1,i),xmn)
	      ymn = min(x(2,i),ymn)
	      vmn = min(vmn,b(j))
	      vmx = max(vmx,b(j))
	      j   = j + ndf
17    continue
      delx  = (xmx-xmn)
      dely  = (ymx-ymn)
      deleps = 1.e-7    !cww  for values nearly zero?
      delx = max(delx,deleps)
      dely = max(dely,deleps)
      xmx = 8.2/delx
      ymx = 8.2/dely
      if(vmn.eq.vmx) then
cww     if(ior.ge.0) write(iow,2001)
	      if(ior.lt.0)
     +  call drawmess(' No plot - zero difference in values ',1,0)
	      return
      end if
c.... open plot and loop through elements
c.... set colortable for nnc.le.14
      call colpal(ipal,nnc)
c.... set direction for colors
      ncol = nc
      if(mc.gt.0) ncol = nc-1
      call pltcol(ipal,ipali,ncol,mmc,cinv)
      iclear = 0
      call plopen
      call pzero(xl,120)
c.... set colors for contour plots
      call setcol(2,0)
    	if(cont) then
	      call pltctx(vc,ic,nc,mmc,1,ipali)
	    else
	      call pltftx(vc,-mc,nc,mmc,1,ipali,vmin,vmax)
	    end if
      ic = max(1,min(ic,ndf))
c.... set for ps(prin) that no border is drawn
      ibps = 1
c.....loop over elements
      do 250 n = 1,numel1
c....	  check if element is in window and set values of vl and vu
	      vl = vmx
    	  vu = vmn
    	  ns = 0
	      do 110 i = 1,iu
	        ii = id(iplt(i),n)
	        if(ii.gt.0) then
	          ns = ns + 1
	          xl(1,ns) = x(1,ii)
	          xl(2,ns) = x(2,ii)
	          if(ndm.ge.3) xl(3,ns) = x(3,ii)
	          j = ndf*(ii-1) + ic
	          v(ns) = b(j)
	          vl = min(vl,v(ns))
	          vu = max(vu,v(ns))
	        end if
110	    continue
	      if ( zoom(xl,3,ns-1) ) then
c....	    get center values - assign to node 1 and 4 of triangle
          if(ns.gt.4) then
	          nsi = 1
	          xt(1,1) = 0.
	          xt(2,1) = 0.
	          xt(3,1) = 0.
	          vt(1) = 0.
	          do 120 i = 1,ns-1
	            xt(1,1) = xt(1,1) + xl(1,i)
	            xt(2,1) = xt(2,1) + xl(2,i)
	            if(ndm.ge.3) xt(3,1) = xt(3,1) + xl(3,i)
	            vt(1) = vt(1) + v(i)
120	        continue
	          xnn     = ns - 1
	          xt(1,1) = xt(1,1)/xnn
	          xt(2,1) = xt(2,1)/xnn
	          xt(3,1) = xt(3,1)/xnn
	          vt(1)   = vt(1)/xnn
  	      else
	          nsi     = 2
	          xt(1,1) = xl(1,1)
	          xt(2,1) = xl(2,1)
	          xt(3,1) = xl(3,1)
	          vt(1)   = v(1)
	        end if
	        xt(1,4) = xt(1,1)
	        xt(2,4) = xt(2,1)
	        xt(3,4) = xt(3,1)
	        vt(4) = vt(1)
c....	    loop over subtriangles
	        do 240 ii = nsi,ns-nsi
c....	      set other points on the triangle
	          xt(1,2) = xl(1,ii)
	          xt(2,2) = xl(2,ii)
	          xt(3,2) = xl(3,ii)
	          vt(2)   = v(ii)
	          xt(1,3) = xl(1,ii+1)
	          xt(2,3) = xl(2,ii+1)
	          xt(3,3) = xl(3,ii+1)
	          vt(3)   = v(ii+1)
	          if(cont) then
c....	        plot all contours which intersect this element
c....         define color table
	            do 230 nn = 1,nc
                vv = vc(nn)
                if(vv.ge.vl.and.vv.le.vu) then
                  call pppcol(iadd1+ipali(nn))
c....	            loop over sides of the triangle to find plot points
                  j = 3
		              do 220 i = 1,3
csa                 if(vv.eq.vt(i)) then
                    if(dabs(vv-vt(i)).lt.1.e-5 ) then
		                  x1 = xt(1,i)
		                  y1 = xt(2,i)
		                  z1 = xt(3,i)
		                  call plotl(x1,y1,z1,j)
		                  j = 2
                    else if((vt(i)-vv)*(vt(i+1)-vv).lt.0.0d0) then
		                   s = (vv - vt(i))/(vt(i+1)-vt(i))
		                  x1 = xt(1,i) + s*(xt(1,i+1) - xt(1,i))
		                  y1 = xt(2,i) + s*(xt(2,i+1) - xt(2,i))
		                  z1 = xt(3,i) + s*(xt(3,i+1) - xt(3,i))
		                  call plotl(x1,y1,z1,j)
		                  j = 2
		                end if
c....		            add lables on contour lines
                    if(ipb.eq.0 .and. j.eq.2) then
                      ivc = (x1-xmn)*xmx + 1
                      ivc = max(1,min(9,ivc))
                      jvc = (y1-ymn)*ymx + 1
                      jvc = max(1,min(9,jvc))
                      if(tvc(ivc,jvc)) then
                        tvc(ivc,jvc) = .false.
                        call plotl(x1-dx1,y1,z1,3)
                        call plabl(nn)
                        call plotl(x1,y1,z1,3)
                      end if
                    end if
220               continue
                end if
230           continue
            else
              call pltcor(3,ilc,vt,vc,nc)
              call pltefl(3,ilc,xt,vt,vc,nc,mmc,ipla,ipgl,idev,ipali)
            end if
240       continue
        end if
c....   plot element border in color ct4
        if(k4.gt.0) then
          call pppcol(k4)
          if ( zoom(xl,3,iu) ) then
            call plotl(xl(1,1),xl(2,1),xl(3,1),3)
            do ib = 2,iu
              call plotl(xl(1,ib),xl(2,ib),xl(3,ib),2)
            end do
          end if
        end if
250   continue
      ibps = 0
c.... set colors back
      call setcol(1,0)
      return
1000  format(a1)
2000  format('   ------ Contour Values for Plot ------'/
     +      (1x,7e11.4/1x,7e11.4))
cww2001  format(' ERROR: No plot - zero difference in values')
2002  format(' Input values correct? (y or n, c = cancel) > ',$)
2003  format(' Input No. of Lines/Colors, Min. Value, Max. Value for ',
     +       ' Plot'/3x,'>',$)
      end
c
      subroutine thinsec(d)
c-----------------------------------------------------------------------
c
c      Purpose: calulate data for thin cross sections
c
c      Inputs:

c      Outputs:
c
c-----------------------------------------------------------------------
      USE eldata
      USE psize
      USE sectio
      implicit double precision (a-h,o-z)
      common m(maxm)
      dimension d(*)

      call dinput(d(3),2)
      numnps(1,ma) = d(3)
      numels(1,ma) = d(4)

c.... set arrays for thin cross sections
      call thinset

c.... read input for thin cross sections
      call thininp(m(msec(5,ma)),numnps(1,ma),numels(1,ma),
     +                           numnps(2,ma),numels(2,ma))

c.... set arrays for thick cross sections
      call thicset(0)

c.... compute section properties
      call thinpro(d,m(msec(5,ma)),m(msec(11,ma)),m(msec(12,ma)),
     +             alpha,numels(1,ma))

c.... compute warping function
      call thinwar(d,m(msec(5,ma)),m(msec(3,ma)),m(msec(12,ma)),
     +             numnps(1,ma),numels(1,ma))

c.... save element data in m-array for plot
      call savesec(d,m(msec(5,ma)),m(msec(1,ma)),m(msec(2,ma)),
     +             m(msec(10,ma)),m(msec(12,ma)),numels(1,ma))

c.... create and save 4-node-mesh for num. integration
      call sectmesh(d,m(msec(5,ma)),m(msec(6,ma)),m(msec(7,ma)),
     1                m(msec(3,ma)),m(msec(8,ma)),numels(1,ma))


      return
      end
c
      subroutine thininp(dat,numnp1,numel1,numnp2,numel2)
c-----------------------------------------------------------------------
c.... read input for thin cross sections (dat[1..9])
c-----------------------------------------------------------------------
c     dat( 1,i): y-Koordinate des Schwerpunktes
c     dat( 2,i): z-Koordinate des Schwerpunktes
c     dat( 3,i): Laenge des Querschnittsegmentes
c     dat( 4,i): Dicke  des Querschnittsegmentes
c     dat( 5,i): Winkel, um den Segementachsen auf
c                Eingabeachsensystem gedreht werden muessen
c     dat( 6,i): Anfangsknoten des Querschnittsegmentes
c     dat( 7,i):     Endknoten des Querschnittsegmentes
c     dat( 8,i): Anzahl der Elemente ueber Segmentbreite
c     dat( 9,i): Anzahl der Elemente ueber Segmentlaenge
c     dat(10,i): Transformierte y-Koordinate (wird ueberschrieben)
c     dat(11,i): Transformierte z-Koordinate (wird ueberschrieben)
c     dat(12,i): Winkel, um den Segmentachsen auf
c                Hauptachsen des Gesamtquerschnitts gedreht werden
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension dat(12,*)
c
      numnp2 = 0
      numel2 = 0

      do i = 1,numel1
        call dinput(dat(1,i),9)
        ny = dat(8,i)
        nz = dat(9,i)
        numnp2 = numnp2 + (ny+1)*(nz+1)
        numel2 = numel2 + ny*nz
      end do

      return
      end
c
      subroutine thinpro(d,dat,vec,diq,alpha,numel1)
c-----------------------------------------------------------------------
c.... compute section properties
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),dat(12,*),vec(3,*),diq(*)

      tol = 1.d-10
      pi  = dacos(-1.0d0)

C...  Bestimmung der Koordinaten des Schwerpunktes

      ys = 0.0d0
      zs = 0.0d0
      A  = 0.0d0                ! Gesamtquerschnittsflche
      Ay = 0.0d0                ! Statisches Moment bzgl z-Achse
      Az = 0.0d0                ! Statisches Moment bzgl y-Achse
      Ait= 0.0d0                ! Torsionswiderstandsmoment -> It

      do 100 i = 1,numel1
         A  = A  + dat(3,i)*dat(4,i)
         Ay = Ay + dat(3,i)*dat(4,i)*dat(1,i)
         Az = Az + dat(3,i)*dat(4,i)*dat(2,i)
         c1 = 1.0d0-0.63d0*dat(4,i)/dat(3,i)+
     +        0.052d0*(dat(4,i)/dat(3,i))**5
         aIt= aIt+ dat(3,i)*dat(4,i)*dat(4,i)*dat(4,i)/3.0d0*c1
100   continue
      ys = Ay/A
      zs = Az/A

C...  Transformation der Koordinaten der Querschnittssegmente auf
C     System im Schwerpunkt parallel zum Eingabesystem

      Aiy = 0.0d0
      Aiz = 0.0d0
      Aiyz= 0.0d0

      do 110  i = 1,numel1
        dat(10,i) = dat(1,i) - ys
        dat(11,i) = dat(2,i) - zs
        dat(12,i) = dat(5,i)*pi/180.0d0
110   continue

C...  Berechnung der Flaechentraegheitsmomente bezogen auf Achsen in S

      call ftm(dat,numel1,aiy,aiz,aiyz)

      d(3) = A
      d(4) = aIy
      d(5) = aIz
      d(6) = aIyz
      d(7) = aIt
      d(9) = ys
      d(10)= zs

C     Bestimmung der Hauptachsen des Querschnitts
C     Bedingung : Iyz = 0

      if (dabs(aiyz).lt.tol) then
        alpha = 0.0d0
      else if (dabs(aiy-aiz).lt.tol) then
        alpha = pi/4.0d0
      else
       tgalpha = (2.0d0*aiyz)/(aiy-aiz)
       alpha = datan(tgalpha)/2.0d0
      end if

c...  Drehung der Achsen um Winkel Alpha auf Hauptachsen
C     und Ueberschreiben der verschobenen Koordinaten

      do 120 i = 1,numel1
        yy = dat(10,i)
        zz = dat(11,i)
        dat(10,i) = yy*dcos(alpha) + zz*dsin(alpha)
        dat(11,i) = zz*dcos(alpha) - yy*dsin(alpha)
        dat(12,i) = dat(12,i) + alpha
120   continue

c...  Berechnung der Flaechentraegheitsmomente bzgl Hauptachsen

      call transit(aiy,aiz,aiyz,alpha)

c...  groesseres Haupttraegheitsmoment --> Iy
      if (aiy.lt.aiz) then
        aaaiy = aiy
        aiy = aiz
        aiz = aaaiy
        do 130 i = 1,numel1
          dh = dat(10,i)
          dat(10,i) = -dat(11,i)
          dat(11,i) =  dh
          dat(12,i) = dat(12,i) - pi/2.0d0
130     continue
        alpha = alpha - pi/2.d0
      end if

c...  principal axis
      vec(1,1) = dcos(alpha)
      vec(2,1) = -dsin(alpha)
      vec(3,1) = 0.0d0
      vec(1,2) = dsin(alpha)
      vec(2,2) = dcos(alpha)
      vec(3,2) = 0.0d0
      vec(1,3) = 0.0d0
      vec(2,3) = 0.0d0
      vec(3,3) = 1.0d0
      call vnorm(vec(1,1),dummy)
      call vnorm(vec(1,2),dummy)

      diq(1) = aiy
      diq(2) = aiz
      diq(3) = ait
      diq(5) = alpha

      return
      end
c
      subroutine thinwar(d,dat,w,diq,numnp1,numel1)
c-----------------------------------------------------------------------
c.... compute section properties
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),dat(12,*),w(*),diq(*)

      tol = 1.d-10
      pi  = dacos(-1.0d0)

c...  Berechnung der Einheitsverwoelbung bei Drehung um Schwerpunkt
      call ewoelb(1,d,dat,w,ym,zm,numnp1,numel1)

c...  Bestimmung der Woelbmomente Rsy und Rsz
      Rsy  = 0.0d0
      Rsz  = 0.0d0
      do 100 i = 1,numel1
        ia = dat(6,i)
        ie = dat(7,i)
        call ecoor(i,dat,ya,za,ye,ze)
        call trapez(dat(3,i),w(ia),w(ie),za,ze,dRsy)
        Rsy = Rsy + dRsy * dat(4,i)
        call trapez(dat(3,i),w(ia),w(ie),ya,ye,dRsz)
        Rsz = Rsz + dRsz * dat(4,i)
100   continue

c...  Koordinaten des Schubmittelpunktes im Hauptachsensystem
c     und im Bezugssystem
      ym = -Rsy/diq(1)
      zm =  Rsz/diq(2)
      al = -diq(5)
      yme = ym*dcos(al) + zm*dsin(al) + d( 9)
      zme = zm*dcos(al) - ym*dsin(al) + d(10)

c...  Berechnung der Einheitsverwlbung bei Drehung um Schubmittelpunkt
c     Einheitsverwoelbung um Schwerpunkt wird berschrieben
      call ewoelb(2,d,dat,w,ym,zm,numnp1,numel1)

c...  Berechnung des Wlbwiderstandmomentes Iw
      wIw = 0.0d0
      do 110 i = 1,numel1
        ia = dat(6,i)
        ie = dat(7,i)
        call trapez(dat(3,i),w(ia),w(ie),w(ia),w(ie),dIw)
        wIw = wIw + dIw * dat(4,i)
110   continue

      d(8)  = wIw
      d(11) = yme
      d(12) = zme
      d(13) = ym
      d(14) = zm

      return
      end
c
      subroutine savesec(d,dat,ix,coor,cgs,diq,numel1)
c-----------------------------------------------------------------------
c.... save data in m-field
c-----------------------------------------------------------------------
      USE eldata
      implicit double precision (a-h,o-z)
      dimension d(*),dat(12,*),ix(2,*),coor(2,*),cgs(2,*),diq(*)

      pi = dacos(-1.0d0)

      do 10 i = 1,numel1
        ia = dat(6,i)
        ie = dat(7,i)
        wi = pi/180.0d0*dat(5,i)

        ix(1,i)    = ia
        ix(2,i)    = ie
        coor(1,ia) = dat(1,i) + 0.5d0*dat(3,i)*dsin(wi)
        coor(2,ia) = dat(2,i) + 0.5d0*dat(3,i)*dcos(wi)
        coor(1,ie) = dat(1,i) - 0.5d0*dat(3,i)*dsin(wi)
        coor(2,ie) = dat(2,i) - 0.5d0*dat(3,i)*dcos(wi)
10    continue

      cgs(1,1) = d( 9)
      cgs(2,1) = d(10)
      cgs(1,2) = d(11)
      cgs(2,2) = d(12)
      cgs(1,3) = 0.0d0
      cgs(2,3) = 0.0d0

      diq(4) = d(8)

      if (iel.eq.6) then
        d( 8) = d( 9)
        d( 9) = d(10)
        d(10) = d(11)
        d(11) = d(12)
      end if

      return
      end
c
      subroutine sectmesh(d,dat,ixm,xm,w,wm,numel1)
c-----------------------------------------------------------------------
c.... create 4-node mesh
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),dat(12,*),ixm(4,*),xm(2,*),w(*),wm(*)

      pi = dacos(-1.0d0)

      np = 0
      do 100 i = 1,numel1
        ia   = dat(6,i)
        ie   = dat(7,i)
        ny   = dat(8,i)
        nz   = dat(9,i)
        nnp  = (ny+1)*(nz+1)
        dylo = dat(4,i)/ny
        dzlo = dat(3,i)/nz
        dw   = (w(ia)-w(ie))/nz
        alpha= dat(5,i)*pi/180.0d0

c...    extended warping function: nonlinear fct. throught thickness
        call snvec(i,dat,rna,rne,d)
        difrn = rne - rna
        drn   = difrn/nz

c...    coordinates with respect to segments
        xm(1,1+np) = 0.5d0*dat(4,i)
        xm(2,1+np) = 0.5d0*dat(3,i)
        do 110 j = 1,nz+1
        rn = rna + (j-1)*drn
           do 120 k = 1,ny+1
              nn = np + (j-1)*(ny+1) + k
              xm(1,nn) = xm(1,1+np) - (k-1)*dylo
              xm(2,nn) = xm(2,1+np) - (j-1)*dzlo
              wm(nn)   = w(ia) - (j-1)*dw + rn*xm(1,nn)
120        continue
110     continue
        np = np + nnp
100   continue

c...  coordinates with respect to reference point
      np = 0
      do 200 i = 1,numel1
        ny   = dat(8,i)
        nz   = dat(9,i)
        nnp  = (ny+1)*(nz+1)
        alpha= dat(5,i)*pi/180.0d0
        do 220 j = 1,nnp
          y = xm(1,np+j)
          z = xm(2,np+j)
          xm(1,np+j) = y*dcos(alpha) + z*dsin(alpha) + dat(1,i)
          xm(2,np+j) = z*dcos(alpha) - y*dsin(alpha) + dat(2,i)
220     continue
        np = np + nnp
200   continue

c...  nodal values on each element
      np = 0
      ne = 0
      do 300 i = 1,numel1
        ny  = dat(8,i)
        nz  = dat(9,i)
        nnp = (ny+1)*(nz+1)
        nel = ny*nz
        do 310 k = 1,nz
          do 310 j = 1,ny
            ie = (k-1)*ny + j + ne
            ixm(1,ie) = np + j + (k-1)*(ny+1)
            ixm(2,ie) = ixm(1,ie) + 1
            ixm(3,ie) = ixm(1,ie) + ny+2
            ixm(4,ie) = ixm(1,ie) + ny+1
310     continue
        ne = ne + nel
        np = np + nnp
300   continue

      return
      end
c
      subroutine ewoelb(ib,d,dat,w,ym,zm,numnp1,numel1)
c-----------------------------------------------------------------------
C     Berechnung der Einheitsverwoelbung bzgl S bzw. M
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*),dat(12,*),w(*)

c     Berechnung der Abstaende S bzw.M - Querschnittssegment
C     fr ib=1: Drehung um S
C     fr ib=2: Drehung um M
C
c     Berechnung von ws=K0-Integral(rc*ds)
C
      skw  = 0.0d0
      w(1) = 0.0d0
      do 270 i = 1,numel1
        ia = dat(6,i)
        ie = dat(7,i)
        call ecoor(i,dat,ya,za,ye,ze)
        if (ib.eq.2) then
          ya = ya-ym
          ye = ye-ym
          za = za-zm
          ze = ze-zm
        end if
        call dabst(ya,za,ye,ze,rc)
        w(ie) = w(ia) - rc*dat(3,i)
        skw = skw + dat(4,i)*dat(3,i)*(w(ia)+w(ie))*0.5d0
270   continue
      sw = -skw/d(3)

c     Einheitsverwoelbung bezueglich S bzw. M

      do 280 i = 1,numnp1
        w(i) = w(i) + sw
280   continue

      return
      end
c
      subroutine trapez(abl,ai,aj,bi,bj,at)
c-----------------------------------------------------------------------
c     Ueberlagerung zweier Trapezflaechen
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      at = (ai*bi+(ai+aj)*(bi+bj)+aj*bj)*abl/6.0d0
      return
      end
c
      subroutine ftm(dat,numel1,aiy,aiz,aiyz)
c-----------------------------------------------------------------------
c     Berechnung der Flaechentraegheitsmomente 2.Ordnung in beliebigem
c     Koordinatensystem
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension dat(12,*)
      pi = dacos(-1.0d0)

      do 10  i = 1,numel1
        b = dat( 3,i)
        t = dat( 4,i)
        y = dat(10,i)
        z = dat(11,i)
        phi  = dat(5,i)*pi/180.0d0
        aiyh = t*b*b*b/12.0d0
        aizh = b*t*t*t/12.0d0
        aiiyz= 0.0d0
        call transit(aiyh,aizh,aiiyz,phi)
        sty = b*t*z*z
        stz = b*t*y*y
        Aiy = Aiy + aiyh + sty
        Aiz = Aiz + aizh + stz
        Aiyz = Aiyz + aiiyz - y*z*b*t
10    continue
      return
      end
c
      subroutine transit(aiy,aiz,aiyz,phi)
c-----------------------------------------------------------------------
c     Transformation der Traegheitsmomente auf gedrehtes Achsensystem
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)

      ad = aiy+aiz
      sub= aiy-aiz
      dev= aiyz
      c = dcos(2.0d0*phi)
      s = dsin(2.0d0*phi)
      aiy = 0.5d0*ad+0.5d0*sub*c+dev*s
      aiz = 0.5d0*ad-0.5d0*sub*c-dev*s
      aiyz= -0.5d0*sub*s+dev*c
      return
      end
c
      subroutine dabst(ya,za,ye,ze,rc)
c-----------------------------------------------------------------------
c     Berechnung des Abstandes 0-Punkt-Gerade (durch(ya,za) und (ye,ze))
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)

      di = dsqrt(ya*ya+za*za)
      sing = za/di
      cosg = ya/di
      ai  = ze-za
      bi  = ya-ye
      ci  = ya*(za-ze)+za*(ye-ya)
      dif = ze*cosg-ye*sing
      if (dabs(dif).lt.1.d-10) then
        qi = 0.0d0
      else
        qi = dif/dabs(dif)
      end if
      rc = qi*dabs(ci/dsqrt(ai*ai+bi*bi))
      return
      end
c
      subroutine ecoor(i,dat,ya,za,ye,ze)
c-----------------------------------------------------------------------
c     Berechnung der Segmentrandkoordinaten im Hauptachsensystem
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension dat(12,*)

      ya = dat(10,i) + 0.5d0*dat(3,i)*dsin(dat(12,i))
      za = dat(11,i) + 0.5d0*dat(3,i)*dcos(dat(12,i))
      ye = dat(10,i) - 0.5d0*dat(3,i)*dsin(dat(12,i))
      ze = dat(11,i) - 0.5d0*dat(3,i)*dcos(dat(12,i))
      return
      end
c
      subroutine snvec(i,dat,rna,rne,d)
c-----------------------------------------------------------------------
c     Berechnung des Normalenvektors an Segment (Koordinaten bezogen
c     auf Schubmittelpunkt) im Anfangs- und Endpunkt - nicht normiert
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension dat(12,*),d(*),xanvec(2,2),xenvec(2,2)

      call ecoor(i,dat,ya,za,ye,ze)
      ym = d(13)
      zm = d(14)
      sy = ye-ya
      sz = ze-za
      xanvec(1,1) = ya-ym
      xanvec(1,2) = za-zm
      xanvec(2,1) = ya+sz-ym
      xanvec(2,2) = za-sy-zm
      xenvec(1,1) = ye-ym
      xenvec(1,2) = ze-zm
      xenvec(2,1) = ye+sz-ym
      xenvec(2,2) = ze-sy-zm

c     Berechnung der Abstaende der Normalen im Anfangs-und Endpunkt vom
c     Schubmittelpunkt
      call dabst(xanvec(1,1),xanvec(1,2),xanvec(2,1),xanvec(2,2),rna)
      call dabst(xenvec(1,1),xenvec(1,2),xenvec(2,1),xenvec(2,2),rne)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine thicsec(d)
c-----------------------------------------------------------------------
c
c      Purpose: calulate data for thick cross sections
c
c      Inputs:

c      Outputs:
c
c-----------------------------------------------------------------------
      USE eldata
      USE psize
      USE sectio
      implicit double precision (a-h,o-z)
      character*229 fname
      common m(maxm)
      dimension d(*)
c
      numnps(1,ma) = 0
      numels(1,ma) = 0

c.... set arrays for thin cross sections
      call thinset

c.... read filename and numps2,numels from file
      call thicinp1(fname,numnps(2,ma),numels(2,ma))

c.... set arrays for thick cross sections
      call thicset(0)

c.... read input data for thick cross sections from file
      call thicinp2(fname,m(msec(6,ma)),m(msec(7,ma)),m(msec(8,ma)))

c.... compute section properties
      call thicpro_wq(d,m(msec(6,ma)),m(msec(7,ma)),m(msec(8,ma)),
     +             m(msec(10,ma)),m(msec(11,ma)),m(msec(12,ma)),
     +             numnps(2,ma),numels(2,ma),iel)
c

      return
      end

c-----------------------------------------------------------------------
c
      subroutine thicsec2(d,nOrth,nume,numn)
c-----------------------------------------------------------------------
c
c      Purpose: set arrays and read mesh data (thick cross section)
c
c      Inputs:
c        nOrth  -  total number of ortho-parameter
c        nume:  -  number of cross-section elements
c        numn:  -  number of cross-section nodes
c
c      Outputs:
c
c-----------------------------------------------------------------------
      USE eldata
      USE psize
      USE sectio
      implicit double precision (a-h,o-z)
      character*229 fname
      common m(maxm)
      dimension d(*)
c
      numnps(1,ma) = 0
      numels(1,ma) = 0

c.... set arrays for thin cross sections
      call thinset

c.... read numps2 and numels from ifile (mate)
      nen2=4
      ndf2=1
      ndm2=2
      numels(2,ma)=nume
      numnps(2,ma)=numn

c.... set arrays for thick cross sections
      call thicset(nOrth)

c.... read input data for thick cross sections from ifile (mate)
      call thicinp3(d,numnps(2,ma),numels(2,ma),nen2,
     +               m(msec(6,ma)),m(msec(7,ma)))

      return
      end
c
c----------------------------------------------------------------------
c
      subroutine thicinp1(fname,numnp,numel)
c----------------------------------------------------------------------
c
c     Purpose: read name of file and numnps/numel
c
c     Inputs:
c
c     Outputs:
c       fname
c       numnp  - no. of nodes
c       numel  - no. of elements
c
c     Comments:
c         iasc = 0   unformatted
c      >>>iasc = 1   formatted
c
c----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      logical lexist
      character*229 fname
      character*80 t
      data ip1/23/,lexist/.false./
c
c     iasc = 0  ! unformatted
      iasc = 1  ! formatted
c
c.... open files(at the moment including path!!) and read numnp,numel
      read(ior,'(a229)') fname
      inquire(file=fname,exist=lexist)
      if(lexist) then
        if(iasc.eq.0) then
          open(ip1,file=fname,status='unknown',form='unformatted')
        else if(iasc.eq.1) then
          open(ip1,file=fname,status='unknown',form='formatted')
        end if
      else
        stop 'Inputfile in SR thicinp does not exist!'
      end if
      if(iasc.eq.0) then
	      read(ip1) numnp,numel,nummat,ndm,ndf,nen
	    else if(iasc.eq.1) then
        read(ip1,1000) t,numnp,t,numel,t,nummat,t,ndm,t,ndf,t,nen
      end if
c.... close file
      close(ip1)
      return
c
c.... format statements
1000  format(a6,i5,a7,i5,a8,i5,a6,i5,a6,i5,a6,i5)
      end
c
      subroutine thicinp2(fname,ix,x,w)
c----------------------------------------------------------------------
c
c     Purpose: read data fields for thick cross sections
c
c     Inputs:
c       fname - name of data file
c
c     Outputs:
c       x(ndm,numnp)   -  coordinates
c       ix(numel,nen!) -  nodes/element
c       w(numnp)       -  nodal warping function
c
c     Comments:
c         iasc = 0   unformatted
c      >>>iasc = 1   formatted
c
c----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      logical lexist
      character*229 fname
      character*80 t,text
      dimension ix(4,*),x(*),w(*)
      data ip1/23/,lexist/.false./
c
c.... setze ausgabetyp
c     iasc = 0  ! unformatted
      iasc = 1  ! formatted
c

      inquire(file=fname,exist=lexist)
      if(lexist) then
        if(iasc.eq.0) then
          open(ip1,file=fname,status='unknown',form='unformatted')
        else if(iasc.eq.1) then
          open(ip1,file=fname,status='unknown',form='formatted')
        end if
      else
        stop 'Inputfile in SR thicinp does not exist!'
      end if

      if(iasc.eq.0) then
        read(ip1) numnp,numel,nummat,ndm,ndf,nen
c....   nodal coordinates:  field x (ndm,numnp)
   	    read(ip1) (x(i),i=1,ndm*numnp)
c....   element connections:  field ix (nen1,numel)
c       1-nen=nodes,nen+1=nt1,nen+2=nt2,nen+3=nt3,nen+4=mat.number)
        read(ip1)  ((ix(i,k),i=1,nen), du,du,du, k=1,numel)
c....   nodal warping values
        read(ip1) text
	      read(ip1) (w(i),i=1,ndf*numnp)
c
	    else if(iasc.eq.1) then
        read(ip1,1000) t,numnp,t,numel,t,nummat,t,ndm,t,ndf,t,nen
c....   nodal coordinates:  field x (ndm,numnp)
        do k = 1,numnp
          ia = (k-1)*ndm+1
          ie = ia+ndm-1
          read(ip1,1100) t,(x(i),i=ia,ie)
        end do
c....   element connections:  field ix (nen1,numel)
c       1-nen=nodes,nen+1=nt1,nen+2=nt2,nen+3=nt3,nen+4=mat.number)
        do k = 1,numel
          read(ip1,1200) t,(ix(i,k),i=1,nen)
        end do
c....   nodal warping values
        read(ip1,2000) text
        do k = 1,numnp
          read(ip1,1300)  t,w(k)
        end do
      end if
c
c.... close file
      close(ip1)
      return
c
c.... format statements
1000  format(a6,i5,a7,i5,a8,i5,a6,i5,a6,i5,a6,i5)
1100  format(a8,3f14.8)
1200  format(a8,7i6)
1300  format(a8,1e15.8)
2000  format(a80)
      end
c
      subroutine thicpro_wq(d,ixm,xm,wm,cgsr,xvec,sec,numnp,numel,iel)
c-----------------------------------------------------------------------
c     section properties
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension d(*),ixm(4,*),xm(2,*),wm(*),cgsr(2,*),xvec(3,*),sec(*)
     +         ,sg(4),tg(4),shp(3,4),ix(4),xl(2,4),wl(4),xa(10)
      data sg/-1.0, 1.0, 1.0,-1.0/
      data tg/-1.0,-1.0, 1.0, 1.0/
      data tol /1.d-5/
c
      g = 1.0d0/dsqrt(3.0d0)
      call pzero(xa,10)
      call pzero(cgsr,2*3)
      call pzero(xvec,3*3)
c
c...    loop over elements
        do 20 ielem = 1,numel
          do 100 i = 1,4
            ix(i)  = ixm(i,ielem)
            xl(1,i) = xm(1,ix(i))
            xl(2,i) = xm(2,ix(i))
            wl(i)   = wm(ix(i))
100       continue
c...      loop over gauss points
          do 30 l = 1,4
            call shape(sg(l)*g,tg(l)*g,xl,shp,xsj,2,4,ix,.false.)
            y  = 0.0d0
            z  = 0.0d0
            w  = 0.0d0
            wy = 0.0d0
            wz = 0.0d0

            do ii = 1,4
              y  = y  + shp(3,ii)*xl(1,ii)
              z  = z  + shp(3,ii)*xl(2,ii)
              w  = w  + shp(3,ii)*wl(ii)
              wy = wy + shp(1,ii)*wl(ii)
              wz = wz + shp(2,ii)*wl(ii)
            end do

            xa(1) = xa(1) +       xsj
            xa(2) = xa(2) + y   * xsj
            xa(3) = xa(3) + z   * xsj
            xa(4) = xa(4) + y*y * xsj
            xa(5) = xa(5) + z*z * xsj
            xa(6) = xa(6) + y*z * xsj
            xa(7) = xa(7) + w*y  * xsj
            xa(8) = xa(8) + w*z  * xsj
            xa(9) = xa(9) + w*w * xsj
            xa(10) = xa(10) + (wy*wy + wz*wz)*xsj

30        continue
20      continue
c
        ys = xa(2) / xa(1)
        zs = xa(3) / xa(1)
        xiys  = xa(5) - zs*zs*xa(1)
        xizs  = xa(4) - ys*ys*xa(1)
        xiyzs = xa(6) - ys*zs*xa(1)
        det   = xiys*xizs-xiyzs*xiyzs
        ym    = -(xa(8)*xizs-xa(7)*xiyzs)/det
        zm    =  (xa(7)*xiys-xa(8)*xiyzs)/det
        CM    =  xa(9)-xiys*ym*ym-xizs*zm*zm+2.d0*xiyzs*ym*zm
        xit   =  xa(4)+xa(5) - xa(10)
c
c...  store values in d()
      d( 3) = xa(1)
      d( 4) = xiys
      d( 5) = xizs
      d( 6) = xiyzs
      d( 7) = xit
      d( 8) = CM
      d( 9) = ys
      d(10) = zs
      d(11) = ym
      d(12) = zm
c...  store center in cgsr
      cgsr(1,1) = ys
      cgsr(2,1) = zs
      cgsr(1,2) = ym
      cgsr(2,2) = zm
c...  compute principal axes
      xad   = 0.5d0*(xiys-xizs)
      if (dabs(xiyzs).lt.tol .and. dabs(xad).lt.tol) then
        alpha = 0.d0
      else
        alpha = 0.5d0 * datan2(xiyzs,xad)
      end if
      xvec(1,1) = dcos(alpha)
      xvec(2,1) = dsin(alpha)
      xvec(1,2) =-dsin(alpha)
      xvec(2,2) = dcos(alpha)
      xvec(3,3) = 1.d0
      call vnorm(xvec(1,1),dummy)
      call vnorm(xvec(1,2),dummy)
c...  store sectional values
      pi = datan(1.0d0)*4.d0
      sec(1) = xiys
      sec(2) = xizs
      sec(3) = xit
      sec(4) = CM
      sec(5) = alpha*180.d0/pi
      return
      end
c
      subroutine thicpro_wt(d,ixm,xm,wm,cgsr,xvec,sec,numnp,numel,iel)
c-----------------------------------------------------------------------
c     section properties
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension d(*),ixm(4,*),xm(2,*),wm(*),cgsr(2,*),xvec(3,*),sec(*)
     +         ,sg(4),tg(4),shp(3,4),ix(4),xl(2,4),wl(4),xa(10)
      data sg/-1.0, 1.0, 1.0,-1.0/
      data tg/-1.0,-1.0, 1.0, 1.0/
      data tol /1.d-5/
c
      g = 1.0d0/dsqrt(3.0d0)
      call pzero(xa,10)
      call pzero(cgsr,2*3)
      call pzero(xvec,3*3)
c
      do 10 is = 1,2
c...    loop over elements
        do 20 ielem = 1,numel
          do 100 i = 1,4
            ix(i)  = ixm(i,ielem)
            xl(1,i) = xm(1,ix(i))
            xl(2,i) = xm(2,ix(i))
            wl(i)   = wm(ix(i))
100       continue
c...      loop over gauss points
          do 30 l = 1,4
            call shape(sg(l)*g,tg(l)*g,xl,shp,xsj,2,4,ix,.false.)
            y  = 0.0d0
            z  = 0.0d0
            w  = 0.0d0
            wy = 0.0d0
            wz = 0.0d0
            do ii = 1,4
              y  = y  + shp(3,ii)*xl(1,ii)
              z  = z  + shp(3,ii)*xl(2,ii)
              w  = w  + shp(3,ii)*wl(ii)
              wy = wy + shp(1,ii)*wl(ii)
              wz = wz + shp(2,ii)*wl(ii)
            end do
            if (is.eq.2) goto 2
            xa(1) = xa(1) +       xsj
            xa(2) = xa(2) + y   * xsj
            xa(3) = xa(3) + z   * xsj
            xa(4) = xa(4) + y*y * xsj
            xa(5) = xa(5) + z*z * xsj
            xa(6) = xa(6) + y*z * xsj
            xa(7) = xa(7) + wy  * xsj
            xa(8) = xa(8) + wz  * xsj
            xa(9) = xa(9) + w*w * xsj
            goto 30
2           yym = y - ym
            zzm = z - zm
            xa(10) = xa(10) + (yym*yym + zzm*zzm + wz*yym - wy*zzm)*xsj
30        continue
20      continue
        if (is.eq.2) goto 10
        ys = xa(2) / xa(1)
        zs = xa(3) / xa(1)
        xiys  = xa(5) - zs*zs*xa(1)
        xizs  = xa(4) - ys*ys*xa(1)
        xiyzs = xa(6) - ys*zs*xa(1)
        ym =  xa(8)/xa(1) + ys
        zm = -xa(7)/xa(1) + zs
10    continue
c
c...  store values in d()
      d( 3) = xa(1)
      d( 4) = xiys
      d( 5) = xizs
      d( 6) = xiyzs
      d( 7) = xa(10)
      d( 8) = xa(9)
      d( 9) = ys
      d(10) = zs
      d(11) = ym
      d(12) = zm
c...  store center in cgsr
      cgsr(1,1) = ys
      cgsr(2,1) = zs
      cgsr(1,2) = ym
      cgsr(2,2) = zm
c...  compute principal axes
      xad   = 0.5d0*(xiys-xizs)
      if (dabs(xiyzs).lt.tol .and. dabs(xad).lt.tol) then
        alpha = 0.d0
      else
        alpha = 0.5d0 * datan2(xiyzs,xad)
      end if
      xvec(1,1) = dcos(alpha)
      xvec(2,1) = dsin(alpha)
      xvec(1,2) =-dsin(alpha)
      xvec(2,2) = dcos(alpha)
      xvec(3,3) = 1.d0
      call vnorm(xvec(1,1),dummy)
      call vnorm(xvec(1,2),dummy)
c...  store sectional values
      pi = datan(1.0d0)*4.d0
      sec(1) = xiys
      sec(2) = xizs
      sec(3) = xa(10)
      sec(4) = xa(9)
      sec(5) = alpha*180.d0/pi
      return
      end
c
            subroutine thicset(nOrth)
c-----------------------------------------------------------------------
c
c      Purpose: set data fields for thick cross sections
c
c      Inputs:
c           numels2 - number of elements in section
c           numnps2 - number of nodes    in section
c
c      Outputs:
c      msec(12,10) - adresses for associated arrays
c
c      Limits:
c               ma  - 10 materials/section
c                   - 4 node element
c            npstrs - 25+1 stresses to plot
c
c-----------------------------------------------------------------------
      USE cdata
      USE eldata
      USE psize
      USE sectio
      implicit double precision (a-h,o-z)
      logical flg

      common m(maxm)

      npstrs = 26

c.... define m-array, input data for thick segments
      call pseta(msec( 6,ma),4*numels(2,ma),ipr,flg,'SECT-ix4') ! ixm, 4-node mesh
      call pseta(msec( 7,ma),2*numnps(2,ma),ipr,flg,'SECT-xy4') ! xym, 4-node mesh
      call pseta(msec( 8,ma),1*numnps(2,ma),ipr,flg,'SECT- w4') !  wm, 4-node mesh
      call pseta(msec( 9,ma),(npstrs+1)*numnps(2,ma),ipr,flg,'SECT-str')!  npstrs+1 stresses
      call pseta(msec(10,ma),2*3,ipr,flg,'SECT-cent')           !  centre g,s,r
      call pseta(msec(11,ma),3*3,ipr,flg,'SECT-paxi')           !  princip. axis
      call pseta(msec(12,ma),5,ipr,flg,'SECT-ftm ')             !  Iy,Iz,It,Cm,alph
      call pseta(msec(13,ma),nOrth,ipr,flg,'SECT-orth')         !  orth. parameter

      call pzero(m(msec( 9,ma)),(npstrs+1)*numnps(2,ma))

      return
      end

      subroutine thinset
c-----------------------------------------------------------------------
c
c      Purpose: set data fields for thin cross sections
c
c      Inputs:
c           numels1 - number of elements in section
c           numnps1 - number of nodes    in section
c
c      Outputs:
c      msec(12,10) - adresses for associated arrays
c
c      Limits:
c               ma - 10 materials/section
c                  - 2 node element
c
c-----------------------------------------------------------------------
      USE cdata
      USE eldata
      USE sectio
      implicit double precision (a-h,o-z)
      logical flg


c.... define m-array, input data for thin segments
      call pseta(msec(1,ma), 2*numels(1,ma),ipr,flg,'SECT-ix2') !  ix, 2-node mesh
      call pseta(msec(2,ma), 2*numnps(1,ma),ipr,flg,'SECT-xy2') !  xy, 2-node mesh
      call pseta(msec(3,ma), 1*numnps(1,ma),ipr,flg,'SECT- w2') !   w, 2-node mesh
      call pseta(msec(4,ma), 4*numels(1,ma),ipr,flg,'SECT-x+w') ! x+w, scale plot
      call pseta(msec(5,ma),12*numels(1,ma),ipr,flg,'SECT-dat') ! dat, read input

      return
      end

c----------------------------------------------------------------------
      subroutine thicinp3(d,numnp,numel,nen,ix,x)
c----------------------------------------------------------------------
c.... read data arrays for thick cross sections
c----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension ix(nen,*),x(*),d(*),dh(10)

c     nodal coordinates:
      write(iow,'(8x,A,8x,A,14x,A)') 'node','y','z'
      l = 1
      do k = 1,numnp
        call dinput(dh,2)
        write(iow,'(8x,i3,3x,2(e13.5,2x))') k,dh(1),dh(2)
        x(l)=dh(1)
        x(l+1)=dh(2)
        l=l+2
      end do

c     connectivity list:
      write(iow,'(8x,A,8x,4(A,2x))') 'Elem','node 1','node 2',
     +                                      'node 3','node 4'
      do k = 1,numel
        call dinput(dh,4)
        do i=1,nen
          ix(i,k)=dh(i)
        end do
        write(iow,'(8x,i3,9x,4(i4,4x))')
     +            k,ix(1,k),ix(2,k),ix(3,k),ix(4,k)
      end do

      return
      end
