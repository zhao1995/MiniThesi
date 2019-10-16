      subroutine adjchf(nsn,nmn,nsv,msr,id,jd,ndf,nsl)
      implicit double precision (a-h,o-z)
c
c     adjust column heights for contact elements
c
      dimension nsn(*),nmn(*),nsv(*),msr(*),id(ndf,*),jd(*)
c
      ls = 0
      lm = 0
      do 200 n = 1,nsl
	do 100 ni = 1, nsn(n)
	  ii = nsv(ls+ni)
	  do 90 ki = 1, 2
	    nii = id(ki,ii)
	    if (nii .gt. 0) then
	      do 85 nj = 1, nmn(n)
		jj = msr(lm+nj)
		do 80 kj = 1,2
		  njj = id(kj,jj)
		  if (njj .gt. 0) then
		    m = max(nii, njj)
		    jd(m) = max(jd(m), iabs(nii - njj))
		  endif
   80		continue
   85	      continue
	    endif
   90	  continue
  100	continue
	ls = ls + nsn(n)
	lm = lm + nmn(n)
  200 continue
c
      return
      end
c
      subroutine conaug(u,nsln,nmln,nsv,msr,stfsn,stfst,
     1	stfmn,stfmt,x0sn,x0mn,y0sn,y0mn,xdsn,xdmn,ydsn,ydmn,
     2	iwas,iwam,ctrs,ctrm,ifla,ndf,was,wam,fric,itcs,itcm,
     3	slens,slenm,ifpa,islt)
      USE slid3
      implicit double precision(a-h,o-z)
c
c ... routine to perform augmentation of contact Lagrangians
c ...  in a nested scheme
c

      dimension u(*),nsv(*),msr(*),x0sn(*),x0mn(*),y0sn(*),y0mn(*),
     1	xdsn(*),xdmn(*),ydsn(*),ydmn(*),iwas(*),iwam(*),ctrs(*),
     2	ctrm(*),ifla(*),nsln(*),nmln(*),stfsn(*),stfmn(*),was(*),wam(*),
     3	fric(*),itcs(*),stfst(*),stfmt(*),itcm(*),slens(*),slenm(*),
     4	ifpa(*),islt(*)
c
c ... update current positions of contact nodes
c
      call coorud(x0sn,y0sn,xdsn,ydsn,nsntl,u,nsv,ndf)
      call coorud(x0mn,y0mn,xdmn,ydmn,nmntl,u,msr,ndf)
c
c ... set pointer for arrays of different slidelines
c
      ls=1
      lm=1
      ls4=1
      lm4=1
      ls6=1
      lm6=1
      ls10=1
      lm10=1
c
c ... loop over slidelines
c
      do 20, n=1,nsl
	call subcaug(nsln(n),nmln(n),nsv(ls),msr(lm),stfmn(lm),
     1	 stfmt(lm),xdsn(ls),ydsn(ls),xdmn(lm),ydmn(lm),iwas(ls4),
     2	 ctrs(ls10),was(ls6),fric(n),itcs(ls4),slenm(lm),ifpa(n),
     3	 islt(n))
	if (ifla(n) .eq. 1) goto 10
	call subcaug(nmln(n),nsln(n),msr(lm),nsv(ls),stfsn(ls),
     1	 stfst(ls),xdmn(lm),ydmn(lm),xdsn(ls),ydsn(ls),iwam(lm4),
     2	 ctrm(lm10),wam(lm6),fric(n),itcm(lm4),slens(ls),ifpa(n),
     3	 islt(n))
   10 nm=nmln(n)
      ns=nsln(n)
      ls=ls+ns
      lm=lm+nm
      ls4=ls4+4*ns
      lm4=lm4+4*nm
      ls6=ls6+6*ns
      lm6=lm6+6*nm
      ls10=ls10+14*ns
      lm10=lm10+14*nm
   20 continue
      return
      end
c
      subroutine subcaug(nsn,nmn,nsv,msr,stfmn,stfmt,xds,yds,xdm,ydm,
     1	 iwa,ctr,wa,fric,itc,slen,ifpa,islt)
      implicit double precision (a-h,o-z)
c
c ... check gaps, update multipliers
c
      dimension nsv(*),msr(*),stfmn(*),xds(*),yds(*),xdm(*),ydm(*),
     1	    iwa(4,*),ctr(14,*),wa(6,*),stfmt(*),slen(*),itc(4,*)
c
c loop over all slave nodes
c
      do 200, i=1,nsn
c
c  get coordinates of slave,master nodes
c
      xs  = xds(i)
      ys  = yds(i)
      i1  = iwa(3,i)
      i2  = iwa(4,i)
      xkn = stfmn(i1)
      xkt = stfmt(i1)
      if(i2 .ge. 0) then
	x1  = xdm(i1)
	y1  = ydm(i1)
	x2  = xdm(i2)
	y2  = ydm(i2)
c
c  calculate geometrical quantities
c
	tx = x2-x1
	ty = y2-y1
	ax = xs-x1
	ay = ys-y1
	tl = sqrt(tx*tx+ty*ty)
	c  = tx/tl
	s  = ty/tl
	a0 = (ax*c+ay*s)/tl
	if(a0 .gt. 1.0 .or. a0 .lt. 0.0) then
	  ctr(9,i)  = 0.0
	  ctr(12,i) = 0.0
	  itc(4,i)  = 1
	  goto 200
	endif
      endif
      b0   = wa(3,i)
      gapt = 0.0
      nseg = i1-iwa(1,i)
      if(nseg .lt. 0)then
	if(i2 .lt. 0 .and. iwa(2,i) .lt. 0) then
	  do 150, k=1,-nseg
	    gapt = gapt - slen(i1+k-1)
  150	  continue
	elseif(i2 .gt. 0 .and. iwa(2,i) .lt. 0) then
	  gapt = -(1.-a0)*slen(i1)
	  do 155, k=1,-nseg-1
	    gapt = gapt - slen(i1+k)
  155	  continue
	elseif(i2 .lt. 0 .and. iwa(2,i) .gt. 0)then
	  do 160, k=1,-nseg
	    gapt = gapt - slen(i1+k-1)
  160	  continue
	  gapt = gapt - b0*slen(iwa(1,i))
	else
	  gapt = -(1.-a0)*slen(i1)
	  do 165, k=2,-nseg
	    gapt = gapt - slen(i1+k-1)
  165	  continue
	  gapt = gapt - b0*slen(iwa(1,i))
	endif
      else if (nseg .eq. 0) then
	if(iwa(2,i) .lt. 0) b0 = 0.0
	if(i2 .lt. 0)	    a0 = 0.0
	gapt = (a0-b0)*slen(i1)
      else
	if (i2 .lt. 0 .and. iwa(2,i) .lt. 0) then
	  do 170, k=1,nseg
	    gapt = gapt + slen(iwa(1,i)+k-1)
  170	  continue
	else if (i2 .gt. 0 .and. iwa(2,i) .lt. 0)then
	  do 175, k=1,nseg
	    gapt = gapt + slen(iwa(1,i)+k-1)
  175	  continue
	  gapt = gapt + a0*slen(i1)
	else if(i2 .lt. 0 .and. iwa(2,i) .gt. 0)then
	  gapt = (1.-b0)*slen(iwa(1,i))
	  do 180, k=1,nseg-1
	    gapt = gapt + slen(iwa(1,i)+k)
  180	  continue
	else
	  gapt = (1.-b0)*slen(iwa(1,i))
	  do 185, k=1,nseg-1
	    gapt = gapt + slen(iwa(1,i)+k)
  185	  continue
	  gapt = gapt + a0*slen(i1)
	endif
      endif
      if (i2 .ge. 0) then
c****************  segment geometry  ********************************
      gapn = -s*ax+c*ay
c
c  do augmented lagrangian update (normal)
c
      forn=xkn*gapn
      if(forn .le. ctr(9,i)) then
	ctr(9,i) = ctr(9,i) - forn
      else
	ctr(9,i) = 0.0
      endif
c
c  do tangential augmentation
c
      if(fric .eq. 0.0 .or. islt .eq. 3)then
	itc(4,i) = 0
	goto 200
      endif
      itc(4,i) = 1
      if(ctr(9,i) .eq. 0.0) then
	 ctr(12,i) = 0.0
	 goto 200
      endif
      ctr(12,i) = xkt*gapt+ctr(12,i)
      ftrial = abs(ctr(12,i)+wa(1,i))-fric*ctr(9,i)
      if(ftrial .gt. 0.0)then
	itc(4,i) = -2
	ctr(12,i) = ctr(12,i)/abs(ctr(12,i))*fric*ctr(9,i)-wa(1,i)
      endif
      else
c************************  corner geometry *************************
      i0=i1-1
      i2=i1+1
      x0=xdm(i0)
      y0=ydm(i0)
      x1=xdm(i1)
      y1=ydm(i1)
      x2=xdm(i2)
      y2=ydm(i2)
c  calculate normal gap
      tx1=x1-x0
      ty1=y1-y0
      tx2=x2-x1
      ty2=y2-y1
      tl1=sqrt(tx1*tx1+ty1*ty1)
      tl2=sqrt(tx2*tx2+ty2*ty2)
      alpha1=((xs-x0)*tx1+(ys-y0)*ty1)/tl1
      alpha2=((xs-x1)*tx2+(ys-y1)*ty2)/tl2
      gsign=(x2-x0)*ty1-(y2-y0)*tx1
      if(alpha1 .le. 1.0 .and. alpha2 .ge. 0.0)gsign=-gsign
      tx=ys-y1
      ty=x1-xs
      gapn=gsign/abs(gsign)*sqrt(tx*tx+ty*ty)
      forn=xkn*gapn
      if(forn .le. ctr(9,i)) then
	ctr(9,i)=ctr(9,i)-forn
      else
	ctr(9,i)=0.0
      endif
      if (fric .eq. 0.0 .or. islt .eq. 3) then
	itc(4,i)=0
	goto 200
      endif
      itc(4,i)=1
      if (ctr(9,i) .eq. 0.0)then
	ctr(12,i)=0.0
	goto 200
      endif
      ctr(12,i)=xkt*gapt+ctr(12,i)
      ftrial=abs(ctr(12,i)+wa(1,i))-fric*ctr(9,i)
      if (ftrial .gt. 0.0) then
	itc(4,i)=-2
	ctr(12,i)=ctr(12,i)/abs(ctr(12,i))*fric*ctr(9,i)-wa(1,i)
      endif
      endif
  200 continue
      return
      end
c
      subroutine contas(jd,id,u,nsln,nmln,islt,sfacn,sfact,ilocs,
     +	     ilocm,nsv,msr,stfsn,stfst,stfmn,stfmt,slens,slenm,
     +			x0sn,x0mn,y0sn,y0mn,xdsn,xdmn,ydsn,
     +			ydmn,was,wam,iwas,iwam,ctrs,ctrm,itcs,
     +			itcm,ifla,ifpa,fric,ndf)
      USE cdata
      USE slid3
      USE slid4
      implicit double precision(a-h,o-z)
c
c.... contact routine which looks for contact and release
c.... during the iteration process.
c
      dimension jd(*),id(*),nsln(*),nmln(*),islt(*),ilocs(*),
     +		ilocm(*),nsv(*),msr(*),stfsn(*),stfmn(*),slens(*),
     +		slenm(*),x0sn(*),x0mn(*),y0sn(*),y0mn(*),
     +		xdsn(*),xdmn(*),ydsn(*),ydmn(*),was(*),wam(*),
     +		iwas(*),iwam(*),itcs(*),itcm(*),ctrs(*),ctrm(*),
     +		sfacn(*),u(*),ifla(*),ifpa(*),fric(*),stfmt(*),
     +		stfst(*),sfact(*)
c
c.... insert updated displacements
c
      call coorud(x0sn,y0sn,xdsn,ydsn,nsntl,u,nsv,ndf)
      call coorud(x0mn,y0mn,xdmn,ydmn,nmntl,u,msr,ndf)
c
c.... set pointer for arrays of different slidelines
c
      ls  = 1
      lm  = 1
      ls4 = 1
      lm4 = 1
      ls6 = 1
      lm6 = 1
      ls10 = 1
      lm10 = 1
c
c.... loop over all slidelines
c
      do 20 n = 1, nsl
	  call subsl1(nsln(n),nmln(n),islt(n),fric(n),ilocs(ls),
     1	  nsv(ls),msr(lm),stfmn(lm),stfmt(lm),slenm(lm),xdsn(ls),
     2		  ydsn(ls),xdmn(lm),ydmn(lm),was(ls6),iwas(ls4),
     3		  ctrs(ls10),itcs(ls4),ifpa(n),jd,id,ndf)
	  if (islt(n) .lt. 3 .or. nsln(n) .eq. 1) go to 10
	  if (ifla(n) .eq. 1)			  go to 10
	  call subsl1(nmln(n),nsln(n),islt(n),fric(n),ilocm(lm),
     1	  msr(lm),nsv(ls),stfsn(ls),stfst(ls),slens(ls),xdmn(lm),
     2		  ydmn(lm),xdsn(ls),ydsn(ls),wam(lm6),iwam(lm4),
     3		  ctrm(lm10),itcm(lm4),ifpa(n),jd,id,ndf)
   10	  nm  = nmln(n)
	  ns  = nsln(n)
	  ls  = ls   + ns
	  lm  = lm   + nm
	  ls4 = ls4 + 4 * ns
	  lm4 = lm4 + 4 * nm
	  ls6 = ls6 + 6 * ns
	  lm6 = lm6 + 6 * nm
	  ls10 = ls10 + 14 * ns
	  lm10 = lm10 + 14 * nm
   20 continue
c
c.... set new pointer for the profil of global stiffness
c
      call nwprof(jd,neq)
c
      return
      end
c
      subroutine contat(ctrs,ctrm,stfsn,stfst,stfmn,stfmt,itcs,itcm,
     +jd,a,b,id,nsln,nmln,islt,ifla,fric,itan,nsl,ndf,neq,afl,rfl)
      USE rdata
      implicit double precision(a-h,o-z)
c
c.... contact tangent for augmented/penalty formulation
c
      logical afl,rfl
      dimension ctrs(*),ctrm(*),itcs(*),itcm(*),jd(*),a(*),b(*),
     + id(*),nsln(*),nmln(*),islt(*),ifla(*),stfsn(*),stfmn(*),
     + fric(*),itan(*),stfst(*),stfmt(*)
      ls   = 1
      lm   = 1
      ls4  = 1
      lm4  = 1
      ls10 = 1
      lm10 = 1
      do 20 n = 1, nsl
	call cotang(ctrs(ls10),stfsn(ls),stfst(ls),itcs(ls4)
     +	      ,jd,a,b,id,ndf,nsln(n),afl,rfl,fric(n),itan(n))
	if (islt(n) .lt. 3 .or. nsln(n) .eq. 1) go to 10
	if (ifla(n) .eq. 1)			go to 10
	call cotang(ctrm(lm10),stfmn(lm),stfmt(lm),itcm(lm4)
     +	      ,jd,a,b,id,ndf,nmln(n),afl,rfl,fric(n),itan(n))
   10	ns   = nsln(n)
	nm   = nmln(n)
	ls   = ls + ns
	lm   = lm + nm
	ls4  = ls4 + 4 * ns
	lm4  = lm4 + 4 * nm
	ls10 = ls10 + 14 * ns
	lm10 = lm10 + 14 * nm
   20 continue
c
      return
      end
c
      subroutine conupd(ifla,nsln,nmln,islt,was,wam,iwas,iwam,
     1	   ctrs,ctrm,itcs,itcm,nsl)
      implicit double precision(a-h,o-z)
c
c.... contact routine which initializes the t-n+1 data with t-n data
c
      dimension ifla(*),nsln(*),nmln(*),islt(*),was(*),wam(*),
     1		iwas(*),ctrs(*),ctrm(*),itcs(*),itcm(*),iwam(*)
c
c.... set pointer for arrays of different slidelines
c
      ls4 = 1
      lm4 = 1
      ls6 = 1
      lm6 = 1
      ls14=1
      lm14=1
c
c.... loop over all slidelines
c
      do 10 n = 1, nsl
	if(islt(n).eq.4) then
	  call subupd(nsln(n),was(ls6),iwas(ls4),ctrs(ls14),itcs(ls4))
	  if (ifla(n) .ne. 1) then
	    call subupd(nmln(n),wam(lm6),iwam(lm4),ctrm(lm14),itcm(lm4))
	  endif
	endif
	  lm4 = lm4 + 4 * nmln(n)
	  lm6 = lm6 + 6 * nmln(n)
	  lm14=lm14+14*nmln(n)
	  ls4 = ls4 + 4 * nsln(n)
	  ls6 = ls6 + 6 * nsln(n)
	  ls14= ls14+ 14* nsln(n)
   10 continue
c
      return
      end
c
      subroutine coorud(x0,y0,xd,yd,n,u,na,ndf)
      implicit double precision(a-h,o-z)
c
c.... insert updated displacements in the arrays
c.... dealing with contact nodes.
c
      dimension x0(*),y0(*),xd(*),yd(*),u(ndf,*),na(*)
      do 10 i = 1, n
	xd(i) = u(1,na(i)) + x0(i)
	yd(i) = u(2,na(i)) + y0(i)
   10 continue
      return
      end
c
      subroutine cotang(ctr,stfn,stft,itc,jd,a,b,id,ndf,nsn,ufl
     1			,rfl,fric,itan)
      USE cdata
      USE ndata
      implicit double precision(a-h,o-z)
c
c.... form contact tangent and residual
c
      logical lfl,ufl,rfl
      dimension ctr(14,*),itc(4,*),sc(6,6),jd(*),ld(6),a(*),b(*),
     +		id(ndf,*),cl(6),rl(6),t(6),stfn(*),stft(*)

c
      do 150 i = 1, nsn
	if (itc(1,i) .ne. 0) then
	 if(itc(3,i) .ge. 0)then
c
c.... set up connection array,do assy. for 3-node contact element
c
	  ld(1)  = id(1,itc(1,i))
	  ld(2)  = id(2,itc(1,i))
	  ld(3)  = id(1,itc(2,i))
	  ld(4)  = id(2,itc(2,i))
	  ld(5)  = id(1,itc(3,i))
	  ld(6)  = id(2,itc(3,i))
	  penn =  ctr(1,i)
	  pent =  ctr(13,i)
	  ss   = ctr(2,i)
	  cc   = ctr(3,i)
	  a0   = ctr(4,i)
	  a1   = 1.0d0 - a0
	  call pzero (sc,36)
	  call pzero (rl,6)
c  compute normal contribution to stiffness and residual
	    if (ufl) then
	      cl(1) = -ss
	      cl(2) =  cc
	      cl(3) =  ss*a1
	      cl(4) = -cc*a1
	      cl(5) =  ss*a0
	      cl(6) = -cc*a0
	      do 25 k = 1, 6
		facsc = cl(k)*penn
		do 20 kk = k, 6
		  sc(k,kk) = sc(k,kk) + facsc*cl(kk)
   20		continue
   25	      continue
	      t(1) =  cc
	      t(2) =  ss
	      t(3) = -cc*a1
	      t(4) = -ss*a1
	      t(5) = -cc*a0
	      t(6) = -ss*a0
	      cl(1) = 0.0d0
	      cl(2) = 0.0d0
	      cl(3) = ss
	      cl(4) =-cc
	      cl(5) =-ss
	      cl(6) = cc
	      gap   = ctr(11,i)
	      gp    = ctr(5,i)
	      do 35 k = 1, 6
		facsc = gp*t(k)
		facsc1= gp*cl(k)
		do 30 kk = k, 6
		  sc(k,kk) = sc(k,kk) - facsc*cl(kk) - facsc1*t(kk)
   30		continue
   35	      continue
	      gp =  gp*gap
	      do 45 k  = 3,6
		facsc = gp*cl(k)
		do 40 kk = k,6
		  sc(k,kk) = sc(k,kk) - facsc*cl(kk)
   40		continue
   45	      continue
	    endif
	    if(rfl) then
	      fx    = ctr(7,i)*ss
	      fy    =-ctr(7,i)*cc
	      rl(1) = fx
	      rl(2) = fy
	      rl(3) = -fx*a1
	      rl(4) = -fy*a1
	      rl(5) = -fx*a0
	      rl(6) = -fy*a0
	    endif
c
c  add in symmetric portion of frictional tangent
c
	    gpt=ctr(6,i)
	    ratlen=ctr(10,i)
	    if((itc(4,i) .ne. 0) .and. ufl) then
	    if(itc(4,i) .eq. 1) then
c  (frictional stick)
	      cl(1)=cc
	      cl(2)=ss
	      cl(3)=-cc*a1+gap*ss
	      cl(4)=-ss*a1-gap*cc
	      cl(5)=-cc*a0-gap*ss
	      cl(6)=-ss*a0+gap*cc
	      facm=ratlen*ratlen*pent
	      do 55,k=1,6
		facsc=facm*cl(k)
		do 50, kk=k,6
		  sc(k,kk)=sc(k,kk)+facsc*cl(kk)
   50		continue
   55	      continue
	      endif
c   add in symmetric part common to slip, stick conditions
	      facm=ratlen*gpt
	      facm1=facm*2.*gap
	      cl(1)=-ss
	      cl(2)=cc
	      cl(3)=a1*ss
	      cl(4)=-a1*cc
	      cl(5)=a0*ss
	      cl(6)=-a0*cc
	      t(1)=0.
	      t(2)=0.
	      t(3)=ss
	      t(4)=-cc
	      t(5)=-ss
	      t(6)=cc
	      do 65, k=1,6
		facsc=facm*cl(k)
		facsc1=facm*t(k)
		do 60, kk=k,6
		  sc(k,kk)=sc(k,kk)+facsc*t(kk)+facsc1*cl(kk)
   60		continue
   65	      continue
	      cl(1)=cc
	      cl(2)=ss
	      cl(3)=-a1*cc
	      cl(4)=-a1*ss
	      cl(5)=-a0*cc
	      cl(6)=-a0*ss
	      t(1)=0.0
	      t(2)=0.0
	      t(3)=-cc
	      t(4)=-ss
	      t(5)=cc
	      t(6)=ss
	      do 75, k=1,6
		facsc=facm*cl(k)
		facsc1=facm*t(k)
		do 70, kk=k,6
		  sc(k,kk)=sc(k,kk)-facsc*t(kk)-facsc1*cl(kk)
   70		continue
   75	      continue
	      cl(1)=0.0
	      cl(2)=0.0
	      cl(3)=ss
	      cl(4)=-cc
	      cl(5)=-ss
	      cl(6)=cc
	      do 85, k=1,6
		facsc=facm1*cl(k)
		facsc1=facm1*t(k)
		do 80, kk=k,6
		  sc(k,kk)=sc(k,kk)-facsc*t(kk)-facsc1*cl(kk)
   80		continue
   85	      continue
	      endif
c
c
c.... account for symmetry
c
	  do 95 k =1,6
	  do 95 kk=k,6
	    sc(kk, k) = sc( k,kk)
   95	  continue
c
c  add in non-symmetric portion of stiffness due to slip (for penalty form)
c
	  if(itan .eq. 1) goto 111
	  if((itc(4,i) .eq. -1) .and. ufl) then
	    temp=ctr(8,i)
	    if(temp .ne. 0)temp=temp/abs(temp)
	    facm=ratlen*penn*fric*temp
	    cl(1)=cc
	    cl(2)=ss
	    cl(3)=-a1*cc+gap*ss
	    cl(4)=-a1*ss-gap*cc
	    cl(5)=-a0*cc-gap*ss
	    cl(6)=-a0*ss+gap*cc
	    t(1)=-ss
	    t(2)=cc
	    t(3)=a1*ss
	    t(4)=-a1*cc
	    t(5)=a0*ss
	    t(6)=-a0*cc
	    do 110, k=1,6
	      facsc=facm*cl(k)
	      do 105, kk=1,6
		sc(k,kk)=sc(k,kk)-facsc*t(kk)
  105	      continue
  110	    continue
	  endif
  111	    continue
c
c  add in frictional contribution to residual
c
	  if((itc(4,i) .ne. 0) .and. rfl) then
	      forc=-ctr(8,i)*ratlen
	      rl(1)=rl(1)+forc*cc
	      rl(2)=rl(2)+forc*ss
	      rl(3)=rl(3)+forc*(-a1*cc+gap*ss)
	      rl(4)=rl(4)+forc*(-a1*ss-gap*cc)
	      rl(5)=rl(5)+forc*(-a0*cc-gap*ss)
	      rl(6)=rl(6)+forc*(-a0*ss+gap*cc)
	  endif
c
c.... add tangent/residual to global stiffness
c
	  lfl = .not.(nau .eq. nal)

        call dasbly(sc,rl,ld,jd,6,lfl,ufl,rfl,b,gstiff(nal),gstiff(nau)
     +             ,a)
	 else
c
c .... set up connection array, do assy for 2-node contact element
c
	  ld(1)=id(1,itc(1,i))
	  ld(2)=id(2,itc(1,i))
	  ld(3)=id(1,itc(2,i))
	  ld(4)=id(2,itc(2,i))
	  ld(5)=id(1,itc(2,i))
	  ld(6)=id(2,itc(2,i))
	  penn=ctr(1,i)
	  pent=ctr(13,i)
	  ss=ctr(2,i)
	  cc=ctr(3,i)
	  call pzero(sc,36)
	  call pzero(rl,6)
c  compute normal contribution to stiffness and residual
	  if(ufl)then
	    facsc = ctr(1,i)
	    do 115, k=1,4
	      sc(k,k) = sc(k,k)+facsc
  115	    continue
	    sc(1,3) = sc(1,3)-facsc
	    sc(2,4) = sc(2,4)-facsc
	    if (ctr(5,i) .ne. 0.0)then
	      cl(1) =-ss
	      cl(2) = cc
	      cl(3) = ss
	      cl(4) =-cc
	      do 120, k=1,4
		facsc = ctr(5,i)*cl(k)
		do 118, kk=k,4
		  sc(k,kk) = sc(k,kk)+facsc*cl(kk)
  118		continue
  120	      continue
	    endif
c
c... account for symmetry
c
	    do 125, k=1,4
	    do 125, kk=k,4
	      sc(kk,k)=sc(k,kk)
  125	    continue
	  endif
	  if(rfl) then
	    fx	  = ctr(7,i)*ss
	    fy	  =-ctr(7,i)*cc
	    rl(1) = fx
	    rl(2) = fy
	    rl(3) =-fx
	    rl(4) =-fy
	  endif
c
c ... add in frictional contribution to tangent
c
	  if (itc(4,i) .ne. 0 .and. ufl) then
	    cl(1) = cc
	    cl(2) = ss
	    cl(3) =-cc
	    cl(4) =-ss
	    t(1)  =-ss
	    t(2)  = cc
	    t(3)  = ss
	    t(4)  =-cc
	    do 130, k=1,4
	      facsc=ctr(6,i)*t(k)
	      do 128, kk=1,4
		sc(k,kk)=sc(k,kk)-facsc*cl(kk)
  128	      continue
  130	    continue
c	    if(itc(4,i) .ne. -1) goto 140
c------------- inserted by c.miehe
	    if(itc(4,i) .ne. -1) goto 141
c---------------------------------
	    temp = ctr(8,i)
	    if(temp .ne. 0) temp = temp/abs(temp)
	    facm = penn*fric*temp
	    do 140, k=1,4
	      facsc = facm*cl(k)
	      do 138, kk=1,4
		sc(k,kk) = sc(k,kk)-facsc*t(kk)
  138	      continue
  140	    continue
c------------- inserted by c.miehe
  141	    continue
c---------------------------------
	  endif
	  if(itc(4,i) .ne. 0 .and. rfl)then
	    forc  = -ctr(8,i)
	    rl(1) = rl(1)+forc*cc
	    rl(2) = rl(2)+forc*ss
	    rl(3) = rl(3)-forc*cc
	    rl(4) = rl(4)-forc*ss
	  endif
c ... add tangent/residual to global arrays
	  lfl=.not.(nau .eq. nal)

        call dasbly(sc,rl,ld,jd,6,lfl,ufl,rfl,b,gstiff(nal),gstiff(nau)
     +              ,a)
	 endif
	endif
  150 continue
      return
      end
c
      subroutine csegupd(u,nsln,nmln,nsv,msr,x0sn,x0mn,
     1	    y0sn,y0mn,xdsn,xdmn,ydsn,ydmn,slens,slenm,ndf)
      USE slid3
      implicit double precision (a-h,o-z)
c
c ... routine to update arrays containing segment lengths
c
      dimension u(*),nsln(*),nmln(*),nsv(*),msr(*),x0sn(*),
     1	   x0mn(*),y0sn(*),y0mn(*),xdsn(*),xdmn(*),ydsn(*),
     2	   ydmn(*),slens(*),slenm(*)
c
c ... update current positions of contact nodes
c
      call coorud(x0sn,y0sn,xdsn,ydsn,nsntl,u,nsv,ndf)
      call coorud(x0mn,y0mn,xdmn,ydmn,nmntl,u,msr,ndf)
c
c ... set pointer for arrays of different slidelines
c
      ls = 1
      lm = 1
c
c ... loop over slidelines
c
      do 20, n=1,nsl
	call subsegu(nsln(n),xdsn(ls),ydsn(ls),slens(ls))
	call subsegu(nmln(n),xdmn(lm),ydmn(lm),slenm(lm))
	nm = nmln(n)
	ns = nsln(n)
	ls = ls+ns
	lm = lm+nm
   20 continue
      return
      end
c
c
c
      subroutine subsegu(nmn,xdm,ydm,slen)
      implicit double precision(a-h,o-z)
      dimension xdm(*),ydm(*),slen(*)
      if(nmn .eq. 1) return
      nmax = nmn-1
      do 100, i=1,nmax
	x1 = xdm(i)
	y1 = ydm(i)
	x2 = xdm(i+1)
	y2 = ydm(i+1)
	slen(i) = sqrt((x2-x1)**2+(y2-y1)**2)
  100 continue
      return
      end
c
      subroutine incon(x,ndm,prt)
c-----------------------------------------------------------------------
c     input subroutine for contact data
c
c     parameters
c
c     nsl     - number of slidelines
c     nsntl   - total number of slave nodes
c     nmntl   - total number of master nodes
c     naxi    - flag (0 = plane, 1 = axisymmetric problem)
c     numel   - total number of elements
c     x(ndm,*)- coordinates of nodes
c
c     set pointer for allocation of contact data
c
c     l00 - number of slave nodes per slideline    (nsln ) - nsl
c     l01 - number of master nodes per slideline   (nmln ) - nsl
c     l02 - slide line type			   (ilst ) - nsl
c     l03 - penalty for slideline (normal)	   (sfacn) - nsl
c     l33 - penalty for slideline (tangential)	   (sfact) - nsl
c     l04 - location of slave nodes		   (ilocs) - nsntl
c     l05 - location of master nodes		   (ilocm) - nmntl
c     l06 - node number of slave node		   (nsv  ) - nsntl
c     l07 - node number of master node		   (msr  ) - nmntl
c     l08 - stiffness of slave surface (normal)    (stfsn) - nsntl
c     l34 - stiffness of slave surface (tang.)	   (stfst) - nsntl
c     l09 - stiffness of master surface (normal)   (stfmn) - nmntl
c     l35 - stiffness of master surface (tang.)    (stfmt) - nmntl
c     l10 - segment lengths, slave		   (slens) - nsntl
c     l11 - segment lengths, master		   (slenm) - nsntl
c     l12 - x-coord slave (reference config)	   (x0sn ) - nsntl
c     l13 - x-coord master(reference config)	   (x0mn ) - nmntl
c     l14 - y-coord slave (reference config)	   (y0sn ) - nsntl
c     l15 - y-coord master(reference config)	   (y0mn ) - nmntl
c     l16 - x-coord slave (current   config)	   (xdsn ) - nsntl
c     l17 - x-coord master(current   config)	   (xdmn ) - nmntl
c     l18 - y-coord slave (current   config)	   (ydsn ) - nsntl
c     l19 - y-coord master(current   config)	   (ydmn ) - nmntl
c     l20 - frict. force + surf. coord. slave	   (was  ) - 6*nsntl
c     l21 - frict. force + surf. coord. master	   (wam  ) - 6*nmntl
c     l22 - location of fric. surface	slave	   (iwas ) - 4*nsntl
c     l23 - location of fric. surface	master	   (iwam ) - 4*nmntl
c     l24 - real data for slave tangents	   (ctrs ) - 14*nsntl
c     l25 - real data for master tangents	   (ctrm)  - 14*nmntl
c     l26 - integer data for slave tangents	   (itcs)  - 4*nsntl
c     l27 - integer data for master master	   (itcm)  - 4*nmntl
c     l28 - flag to avoid change of master & slave (ifla)  - nsl
c		0 = two pass
c		1 = one pass
c     l29 - switch for determining solution method (ifpa)  - nsl
c		0 = penalty method
c		1 = augmented lagrangian method,simult. iteration
c		2 = augmented lagrangian method,nested iteration,
c			symmetrized treatment of friction
c		3 = augmented lagrangian method,nested iteration,
c			nonsymmetric treatment of friction
c     l30 - friction coefficients on each slideline(fric)  - nsl
c     l31 - input flag for contact tangent option	   (itan)  - nsl
c	    for use in "cotang" subroutine
c		0 = consistent tangent
c		1 = omit nonsymmetric frictional contribution
c
c-----------------------------------------------------------------------
      USE cdata
      USE iofile
      USE psize
      USE slid1
      USE slid2
      USE slid3
      USE slid4
      USE slid5
      USE doalloc
      implicit double precision(a-h,o-z)
      logical aafl,prt
      dimension x(ndm,*),td(4)
c
c.... input global parameters
c
      call dinput(td,4)
      nsl   = td(1)
      nsntl = td(2)
      nmntl = td(3)
      naxi  = td(4)
		   write(iow,20) nsl,nsntl,nmntl,naxi
      if(ior.lt.0) write(*  ,20) nsl,nsntl,nmntl,naxi
c
c.... set pointers (allocate space)
c
      call ialloc(cl00,nsl,'CONT-L100',aafl)
      call ialloc(cl01,nsl,'CONT-L101',aafl)
      call ialloc(cl02,nsl,'CONT-L102',aafl)
      call ralloc(cl03,nsl,'CONT-L103',aafl)
      call ralloc(cl33,nsl,'CONT-L133',aafl)
      call ialloc(cl28,nsl,'CONT-L128',aafl)
      call ialloc(cl29,nsl,'CONT-L129',aafl)
      call ralloc(cl30,nsl,'CONT-L130',aafl)
      call ialloc(cl31,nsl,'CONT-L131',aafl)

c.... input nsl - slave and master nodes (get sizes)
c
      call slcont(cl00,cl01,cl02,cl03,cl33,cl29,cl28,cl31,prt)
c
c.... set remaining pointers
c
      call ialloc(cl04,nsntl  ,'CONT-L104',aafl)
      call ialloc(cl05,nmntl  ,'CONT-L105',aafl)
      call ialloc(cl06,nsntl  ,'CONT-L106',aafl)
      call ialloc(cl07,nmntl  ,'CONT-L107',aafl)
      call ralloc(cl08,nsntl  ,'CONT-L108',aafl)
      call ralloc(cl34,nsntl  ,'CONT-L134',aafl)
      call ralloc(cl09,nmntl  ,'CONT-L109',aafl)
      call ralloc(cl35,nmntl  ,'CONT-L135',aafl)
      call ralloc(cl10,nsntl  ,'CONT-L110',aafl)
      call ralloc(cl11,nmntl  ,'CONT-L111',aafl)
      call ralloc(cl12,nsntl  ,'CONT-L112',aafl)
      call ralloc(cl13,nmntl  ,'CONT-L113',aafl)
      call ralloc(cl14,nsntl  ,'CONT-L114',aafl)
      call ralloc(cl15,nmntl  ,'CONT-L115',aafl)
      call ralloc(cl16,nsntl  ,'CONT-L116',aafl)
      call ralloc(cl17,nmntl  ,'CONT-L117',aafl)
      call ralloc(cl18,nsntl  ,'CONT-L118',aafl)
      call ralloc(cl19,nmntl  ,'CONT-L119',aafl)
      call ralloc(cl20,6*nsntl,'CONT-L120',aafl)
      call ralloc(cl21,6*nmntl,'CONT-L121',aafl)
      call ialloc(cl22,4*nsntl,'CONT-L122',aafl)
      call ialloc(cl23,4*nmntl,'CONT-L123',aafl)
      call ralloc(cl24,14*nsntl,'CONT-L124',aafl)
      call ralloc(cl25,14*nmntl,'CONT-L125',aafl)
      call ialloc(cl26,4*nsntl,'CONT-L126',aafl)
      call ialloc(cl27,4*nmntl,'CONT-L127',aafl)
c
c.... read and initialize slideline data (penalty)
c
      if (cl29(1) .le. 3) then
      	call slint0(cl00,cl01,cl02,cl03,cl33,cl04,cl05,
     +	  cl06,cl07,cl08,cl34,cl09,cl35,cl10,cl11,nsl,
     +	  x,ndm,numel,cl20,cl21,cl22,cl23,cl28,cl30,prt)
      else
		    write(iow,2002) cl29(1)
        if(ior.lt.0) write(*  ,2002) cl29(1)
	      stop
      end if
c
c.... initialize coordinates for contact nodes
c
      call slint1(x,ndm,cl06,cl12,cl14,nsntl)
      call slint1(x,ndm,cl07,cl13,cl15,nmntl)
c
c.... formats
c
   20 format(//5x,'c o n t a c t   i n f o r m a t i o n',
     +	     /,5x,' total # of slide lines  ',i4,
     +	     /,5x,' total # of slave nodes  ',i4,
     +	     /,5x,' total # of master nodes ',i4,
     +	     /,5x,' plane/axisym. (0/1)     ',i4)
 2002 format('*** INCON *** Unknown option -> ',i3)
c
      return
      end
c
      subroutine penck(x,ndm,nsn,nmn,islt,msr,nsv,iloc,slen,iwa,wa)
c
c.... check for penetration of slave nodes through master surface
c.... change coordinates if there is a penetration.
c
      USE iofile
      implicit double precision(a-h,o-z)
      dimension x(ndm,*),msr(*),nsv(*),iloc(*),slen(*),iwa(4,*),wa(6,*)
cww   data tol/.00000001d0/
c
c.... loop over slave nodes
c
      do 100 i = 1, nsn
	n  = nsv(i)
	xs = x(1,n)
	ys = x(2,n)
	l  = msr(1)
	al = sqrt((x(1,l)-xs)**2 + (x(2,l)-ys)**2)
	dm = al
	jm = 1
c
c.... find the closest master node "jm"
c
	if (al .ne. 0.0d0) then
	  do 20 j = 2,nmn
	    m = msr(j)
	    bl = sqrt((x(1,m)-xs)**2 + (x(2,m)-ys)**2)
	    if (bl .lt. dm) then
	      dm = bl
	      jm = j
	    endif
   20	  continue
	endif
cww50	i0 = max(  1,jm-1)
	i0 = max(  1,jm-1)
	i1 = jm
	i2 = min(nmn,jm+1)
	iloc(i) = jm
	k = msr(i0)
	l = msr(i1)
	m = msr(i2)
c
c.... decide on master segment, compute geometrical quantities, decide
c.... on type of geometry, test for penetration
c
	if(i1 .ne. i0 .and. i1 .ne. i2) then
	  tx1=x(1,l)-x(1,k)
	  ty1=x(2,l)-x(2,k)
	  tx2=x(1,m)-x(1,l)
	  ty2=x(2,m)-x(2,l)
	  tl1=sqrt(tx1*tx1+ty1*ty1)
	  tl2=sqrt(tx2*tx2+ty2*ty2)
	  tx1=tx1/tl1
	  ty1=ty1/tl1
	  tx2=tx2/tl2
	  ty2=ty2/tl2
	  alpha1=((xs-x(1,k))*tx1+(ys-x(2,k))*ty1)/tl1
	  alpha2=((xs-x(1,l))*tx2+(ys-x(2,l))*ty2)/tl2
	  if(alpha1 .gt. 1.0 .and. alpha2 .lt. 0.0)then
	    i1=i1
	    i2=-1
	    gsign=(x(1,m)-x(1,k))*ty1-(x(2,m)-x(2,k))*tx1
	    if(gsign .ne. 0.0)gsign=gsign/abs(gsign)
	    tx=ys-x(2,l)
	    ty=x(1,l)-xs
	    gapn=gsign*sqrt(tx*tx+ty*ty)
	  else if(alpha1 .lt. 1.0 .and. alpha2 .gt. 0.0)then
	    d1=(1.-alpha1)*tl1
	    d2=alpha2*tl2
	    if(d1 .gt. d2)then
	      i2=i1
	      i1=i0
	      gapn=-(xs-x(1,k))*ty1+(ys-x(2,k))*tx1
	      a0=alpha1
	      tl=tl1
	    else if (d1 .eq. d2) then
	      i1=i1
	      i2=-1
	      gsign=-(x(1,m)-x(1,k))*ty1+(x(2,m)-x(2,k))*tx1
	      if(gsign .ne. 0.0)gsign=gsign/abs(gsign)
	      tx=ys-x(2,l)
	      ty=x(1,l)-xs
	      gapn=gsign*sqrt(tx*tx+ty*ty)
	    else
	      i2=i2
	      i1=i1
	      gapn=-(xs-x(1,l))*ty2+(ys-x(2,l))*tx2
	      a0=alpha2
	      tl=tl2
	    endif
	  else if(alpha1 .le. 1.0 .and. alpha2 .le. 0.0)then
	    i2=i1
	    i1=i0
	    gapn=-(xs-x(1,k))*ty1+(ys-x(2,k))*tx1
	    a0=alpha1
	    tl=tl1
	  else
	    i1=i1
	    i2=i2
	    gapn=-(xs-x(1,l))*ty2+(ys-x(2,l))*tx2
	    a0=alpha2
	    tl=tl2
	  endif
	else if(i1 .eq. i0 .or. i1 .eq. i2) then
	  tx=x(1,m)-x(1,k)
	  ty=x(2,m)-x(2,k)
	  tl=sqrt(tx*tx+ty*ty)
	  tx=tx/tl
	  ty=ty/tl
	  alpha=((xs-x(1,k))*tx+(ys-x(2,k))*ty)/tl
	  if (alpha .lt. -1.d-4) then
	    iwa(3,i)=i0
	    iwa(4,i)=i2
	    wa(6,i)=0.0
	    wa(4,i)=0.0
	    goto 100
	  else if(alpha .gt. 1.0001) then
	    iwa(3,i)=i0
	    iwa(4,i)=i2
	    wa(6,i)=1.0
	    wa(4,i)=0.0
	    goto 100
	  else
	    i1=i0
	    i2=i2
	    gapn=-(xs-x(1,k))*ty+(ys-x(2,k))*tx
	    a0=alpha
	  endif
	endif
	if(gapn .gt. 0.0) then
	  iwa(3,i)=i1
	  iwa(4,i)=i2
	  if(i2 .ge. 0) wa(6,i)=a0
	  goto 100
	else
	  iwa(3,i)=i1
	  iwa(4,i)=i2
	  if(i2 .ge. 0) wa(6,i)=a0
	  if(gapn .eq. 0.0) goto 100
	endif
c
c... change coords of slave node to preclude penetration
c
	if (i2 .ge. 0) then
	  k=msr(i1)
	  l=msr(i2)
	  x(1,n)=(1.-a0)*x(1,k)+a0*x(1,l)
	  x(2,n)=(1.-a0)*x(2,k)+a0*x(2,l)
c	   if(abs(x(1,n)-xs)+abs(x(2,n)-ys) .gt. tol)then
	  if(abs(x(1,n)-xs)+abs(x(2,n)-ys) .gt. 0.0)then
	    write(iow,5000)n,xs,ys,x(1,n),x(2,n)
	    if(ior.lt.0) write(*,5000)n,xs,ys,x(1,n),x(2,n)
	  endif
	else
	  k=msr(i1)
	  x(1,n)=x(1,k)
	  x(2,n)=x(2,k)
c	   if(abs(x(1,n)-xs)+abs(x(2,n)-ys) .gt. tol)then
	  if(abs(x(1,n)-xs)+abs(x(2,n)-ys) .gt. 0.0)then
	    write(iow,5000)n,xs,ys,x(1,n),x(2,n)
	    if(ior.lt.0)write(*,5000)n,xs,ys,x(1,n),x(2,n)
	  endif
	endif
  100 continue
 5000 format(/3x,' Warning of coordinate change in slideline input'
     +	     /3x,' coordinates of node',i6,/5x,'xold=',e12.4,
     +	     5x,'yold=',e12.4,/5x,'xnew=',e12.4,5x,'ynew=',e12.4)
      return
      end
c
      subroutine reverm(nmn,msr)
      dimension msr(*)
c
c     changes the order of msr
c
      nmn1 = nmn + 1
      nmn2 = nmn / 2
      do 10 i = 1, nmn2
	  j = msr(i)
	  msr(i) = msr(nmn1-i)
	  msr(nmn1-i) = j
  10  continue
      return
      end
c
      subroutine slcont(nsln,nmln,islt,sfacn,sfact,ifpa,ifla,itan,prt)
c
c     read slideline control information
c
      USE iofile
      USE slid3
      implicit double precision(a-h,o-z)
      logical prt
      dimension nsln(*),nmln(*),islt(*),sfacn(*),sfact(*),ifla(*),
     +		ifpa(*),itan(*),td(8)
c
c.... input of control data
c
      nsntl = 0
      nmntl = 0
      do 10 i = 1, nsl
	call dinput(td,8)
	nsln(i)  = td(1)
	nmln(i)  = td(2)
	islt(i)  = td(3)
	ifla(i)  = td(4)
	sfacn(i) = td(5)
	sfact(i) = td(6)
	ifpa(i)  = td(7)
	itan(i)  = td(8)
c
c.... print heading
c
	if (prt) then
	  write(iow,50) i,nsln(i),nmln(i),islt(i),
     +			ifla(i),sfacn(i),sfact(i),ifpa(i),itan(i)
	  if(ior.lt.0) write(*,50) i,nsln(i),nmln(i),islt(i),
     +			ifla(i),sfacn(i),sfact(i),ifpa(i),itan(i)
	else
	  write(iow,60) nsl
	  if(ior.lt.0) write(*,60) nsl
	endif
c
	nsntl = nsntl + nsln(i)
	nmntl = nmntl + nmln(i)
   10 continue
c
c.... formats
c
   50 format(//,5x,' s l i d e - l i n e    d a t a',
     +	     /10x,' slide-line number . . . . . . . . . ',i4,
     +	     /15x,' # of slave nodes . . . . . . . ',i4,
     +	     /15x,' # of master nodes  . . . . . . ',i4,
     +	     /15x,' slide-line type  . . . . . . . ',i4,
     +	     /15x,'      1 - sliding only',
     +	     /15x,'      2 - slide-line tied',
     +	     /15x,'      3 - slide + void',
     +	     /15x,'      4 - friction + void',
     +	     /15x,' flag for double/single pass    ',i4,
     +	     /15x,'      0=double pass, 1=single pass',
     +	     /15x,' contact penalty, normal  . . . ',e12.3,
     +	     /15x,' contact penalty, tangential. . ',e12.3,
     +	     /15x,' Solution method  . . . . . . . ',i4,
     +	     /15x,'      0 - Penalty',
     +	     /15x,'      1 - Augmented Lagrangian, simult. iter.',
     +	     /15x,'      2 - Aug. Lagr., nested iter.,symm.',
     +	     /15x,'      3 - Aug. Lagr., nested iter.,nonsymm',
     +	     /15x,' Contact tangent  . . . . . . . ',i4,
     +	     /15x,'      0=Consistent,1=Omit nonsymm. part')
   60 format(//,' contact calculation with ',i3,' slidelines')
c
      return
      end
c
      subroutine slint0(nsln,nmln,islt,sfacn,sfact,ilocs,ilocm,nsv,msr,
     +	    stfsn,stfst,stfmn,stfmt,slens,slenm,nsl,x,ndm,numel,
     +	    was,wam,iwas,iwam,ifla,fric,prt)
      implicit double precision(a-h,o-z)
c
c     input of contact nodes and check of input data
c
      logical prt
      dimension nsln(*),nmln(*),islt(*),sfacn(*),ilocs(*),ilocm(*)
     +	       ,nsv(*),msr(*),stfsn(*),stfmn(*),slens(*),slenm(*)
     +	       ,x(ndm,*),iwas(*),iwam(*),was(*),wam(*),ifla(*),
     +		fric(*),sfact(*),stfst(*),stfmt(*)
      data ls,lm,ls4,lm4,ls6,lm6/1,1,1,1,1,1/
c
c.... loop over all slide-lines
c
      do 20 n = 1, nsl
c
c.... input of contact nodes and storage
c
	  call slnin(nsln(n),nmln(n),islt(n),msr(lm),nsv(ls),
     +		 fric(n),n,prt)
	  call reverm(nsln(n),nsv(ls))
c
c.... store penalty parameters
c
	  call setcpn(nsln(n),nmln(n),sfacn(n),sfact(n),stfsn,stfst,
     +		     stfmn,stfmt,ls,lm)
c
c.... scan over slave (and master) nodes to check input data
c
	  call penck(x,ndm,nsln(n),nmln(n),islt(n),msr(lm),nsv(ls),
     +		     ilocs(ls),slens(ls),iwas(ls4),was(ls6))
	  if((islt(n).ge.3).and.(nsln(n).ne.1).and.(ifla(n).ne.1))
     +	  call penck(x,ndm,nmln(n),nsln(n),islt(n),nsv(ls),msr(lm),
     +		     ilocm(lm),slenm(lm),iwam(lm4),wam(lm6))
c
c.... pointer for next slide-line
c
	  ls  = ls  + nsln(n)
	  lm  = lm  + nmln(n)
	  ls4 = ls4 + 4 * nsln(n)
	  lm4 = lm4 + 4 * nmln(n)
	  ls6 = ls6 + 6 * nsln(n)
	  lm6 = lm6 + 6 * nmln(n)
   20 continue
c
      return
      end
c
c
c
      subroutine setcpn(nsln,nmln,sfacn,sfact,stfsn,stfst,
     +	    stfmn,stfmt,ls,lm)
      implicit double precision (a-h,o-z)
      dimension stfsn(*),stfmn(*),stfst(*),stfmt(*)
	  do 10 i = 1, nsln
	    stfsn(ls+i-1) = sfacn
	    stfst(ls+i-1) = sfact
 10	  continue
	  do 15 i = 1, nmln
	    stfmn(lm+i-1) = sfacn
	    stfmt(lm+i-1) = sfact
 15	  continue
       return
       end
c
      subroutine slint1(x,ndm,ns,x0,y0,ntl)
      implicit double precision(a-h,o-z)
c
c     initialize coordinates of contact nodes
c
      dimension x(ndm,*),ns(*),x0(*),y0(*)
      do 10 i = 1, ntl
	x0(i) = x(1,ns(i))
	y0(i) = x(2,ns(i))
   10 continue
c
      return
      end
c
      subroutine slnin(nsn,nmn,mtt,msr,nsv,fric,nm,prt)
c
c     generate slide-line input
c
      USE iofile
      implicit double precision(a-h,o-z)
      logical prt
      dimension msr(*),nsv(*),td(2)
c
      ik = 0
      il = 0
c
c.... input of frictional constant
c
      call dinput(td,1)
      fric = td(1)
      if (prt) write(iow,180) nm,nsn,nmn,mtt,fric
      if (prt.and.ior.lt.0) write(*,180) nm,nsn,nmn,mtt,fric
      i = ik
c
c.... input and generation of slave nodes
c
10    call dinput(td,2)
      n    = td(1)
      nsvn = td(2)
      if (n .eq. 0) n = i - ik + 1
      n      = n + ik
      nsv(n) = nsvn
      nl     = n-i
      nk     = nl-1
      if(nk) 50,40,20
   20 ndif = nsv(n) - nsv(i)
      nd = ndif / nl
      if(nd*nl .ne. ndif) go to 140
      nk = nk + i - 1
      do 30 j = i, nk
	  j1 = j + 1
	  nsv(j1) = nsv(j) + nd
   30 continue
   40 if(nsn-n+ik) 130,60,50
   50 i=n
      go to 10
   60 i=il
c
c.... input and generation of master nodes
c
70    call dinput(td,2)
      n    = td(1)
      msrn = td(2)
      if(n .eq. 0) n=i-il+1
      n=n+il
      msr(n)=msrn
      nl=n-i
      nk=nl-1
      if(nk) 110,100,80
   80 ndif=msr(n)-msr(i)
      nd=ndif/nl
      if(nd*nl .ne. ndif) go to 140
      nk=nk+i-1
      do 90 j = i, nk
	  j1=j+1
	  msr(j1)=msr(j)+nd
   90 continue
  100 if(nmn-n+il) 130,120,110
  110 i=n
      go to 70
  120 continue
c
c.... write generated data
c
      if(prt) then
	write(iow,190) (i,nsv(i+ik),i=1,nsn)
	if(ior.lt.0) write(*,190) (i,nsv(i+ik),i=1,nsn)
	write(iow,200) (i,msr(i+il),i=1,nmn)
	if(ior.lt.0) write(*,200) (i,msr(i+il),i=1,nmn)
      endif
c
      return
c
c.... error in slide-line input
c
  130 write(iow,220) nm
      if(ior.lt.0) write(*,220) nm
      go to 150
  140 write(iow,210) nm
      if(ior.lt.0) write(*,210) nm
  150 continue
      stop
c
c.... formats
c
  180 format(// 5x,' s l i d e   l i n e ',i4,
     +	     /,10x,' # of slave nodes  ',i4,
     +	     /,10x,' # of master nodes ',i4,
     +	     /,10x,' slide line type   ',i4,
     +	     /,10x,' friction constant ',e14.7)
  190 format(/,5x,' slave  nodes '/(5x,i10,'.',i5))
  200 format(/,5x,' master nodes '/(5x,i10,'.',i5))
  210 format(/,' fatal input error on slide-line ',i5)
  220 format(/,' error on slide-line ',i5)
      end
c
      subroutine subsl0(xs,ys,xm,ym,iloc,nsn,nmn)
      implicit double precision(a-h,o-z)
c
c.... For each slave node, find the closest master node
c.... iloc contains the node number of this master node.
c
      dimension xs(*),ys(*),xm(*),ym(*),iloc(*)
      do 30 i = 1, nsn
	i1 = iloc(i)
   10	i0 = max(  1,i1-1)
	i2 = min(nmn,i1+1)
	al2 = (xm(i0)-xs(i))**2 + (ym(i0)-ys(i))**2
	cl2 = (xm(i2)-xs(i))**2 + (ym(i2)-ys(i))**2
	bl2 = (xm(i1)-xs(i))**2 + (ym(i1)-ys(i))**2
	if    (al2 .lt. bl2 .and. al2 .le. cl2) then
	  i1 = i0
	elseif(cl2 .lt. bl2) then
	  i1 = i2
	else
	  goto 20
	endif
	goto 10
   20	iloc(i) = i1
c     write(*,90) i,iloc(i)
   30 continue
c   90 format(2x,'iloc(' ,i4, ') = ',i4)
      return
      end
c
      subroutine subsl1(nsn,nmn,islt,fric,ilocs,nsv,msr,stfmn,stfmt,
     +	  slen,xds,yds,xdm,ydm,wa,iwa,ctr,itc,ifpa,jd,id,ndf)
c
c.... contact slidelines and voids. check for penetration,
c
      USE edgdat
      USE iofile
      USE psize
      USE slid5
      implicit double precision(a-h,o-z)
      dimension ilocs(*),nsv(*),msr(*),stfmn(*),slen(*),xds(*),
     +		yds(*),xdm(*),ydm(*),wa(6,*),iwa(4,*),jd(*),
     +		ctr(14,*),itc(4,*),id(ndf,*),idl(18),stfmt(*)
c
c.... die if there is only 1 master node
c
      if (nmn .eq. 1) then
	write(iow,2000)
	if(ior.lt.0) write(*,2000)
	stop
      endif
c
c.... find the closest master node for each slave node.
c
      if (islt .ne. 2) then
	call subsl0(xds,yds,xdm,ydm,ilocs,nsn,nmn)
      endif
c
c.... loop over all slave nodes
c
      do 280 i = 1,nsn
	do 3, k=1,3
	 itc(k,i)=0
    3	continue
c  itc(4,i) is set to 2 in subupd to ensure that the first iteration
c   in a time step for frictional problems is done in a stick mode
	if(itc(4,i) .eq. 2) then
	   itc(4,i)=1
	elseif(islt .eq. 4 .and. ifpa .ne. 2)then
	   itc(4,i)=0
	endif
	xs = xds(i)
	ys = yds(i)
	iopt = islt
	iloc = ilocs(i)
	if (iopt .ne. 2) then
c
c.... coordinates of master nodes
c
	  il = max(  1,iloc-1)
	  ir = min(nmn,iloc+1)
	  xc = xdm(iloc)
	  yc = ydm(iloc)
	  xl = xdm(il)
	  yl = ydm(il)
	  xr = xdm(ir)
	  yr = ydm(ir)
c
c  decide on master segment, compute geometrical quantities,
c  determine whether segment or corner geometry exists
c
	  if (iloc .ne. il .and. iloc .ne. ir) then
	    tx1=xc-xl
	    ty1=yc-yl
	    tx2=xr-xc
	    ty2=yr-yc
	    tl1=sqrt(tx1*tx1+ty1*ty1)
	    tl2=sqrt(tx2*tx2+ty2*ty2)
	    tx1=tx1/tl1
	    ty1=ty1/tl1
	    tx2=tx2/tl2
	    ty2=ty2/tl2
	    alpha1=((xs-xl)*tx1+(ys-yl)*ty1)/tl1
	    alpha2=((xs-xc)*tx2+(ys-yc)*ty2)/tl2
	    if(alpha1 .gt. 1.0 .and. alpha2 .lt.  0.0) then
	      i1=iloc
	      i2=-1
	      gsign=(xr-xl)*ty1-(yr-yl)*tx1
	      if(gsign .ne. 0.0)gsign=gsign/abs(gsign)
	      tx=ys-yc
	      ty=xc-xs
	      gapn=gsign*sqrt(tx*tx+ty*ty)
	      c=tx/gapn
	      s=ty/gapn
	    elseif (alpha1 .lt. 1.0 .and. alpha2 .gt. 0.0) then
	      d1=(1.-alpha1)*tl1
	      d2=alpha2*tl2
	      if(d1 .gt. d2) then
		i1=il
		i2=iloc
		gapn=-(xs-xl)*ty1+(ys-yl)*tx1
		a0=alpha1
		c=tx1
		s=ty1
		tl=tl1
	      elseif (d1 .eq. d2) then
		i1=iloc
		i2=-1
		gsign=-(xr-xl)*ty1+(yr-yl)*tx1
		if(gsign .ne. 0.0)gsign=gsign/abs(gsign)
		tx=ys-yc
		ty=xc-xs
		gapn=gsign*sqrt(tx*tx+ty*ty)
		c=tx/gapn
		s=ty/gapn
	      else
		i1=iloc
		i2=ir
		gapn=-(xs-xc)*ty2+(ys-yc)*tx2
		a0=alpha2
		c=tx2
		s=ty2
		tl=tl2
	      endif
	    elseif (alpha1 .le. 1.0 .and. alpha2 .le. 0.0) then
	      i1=il
	      i2=iloc
	      gapn=-(xs-xl)*ty1+(ys-yl)*tx1
	      a0=alpha1
	      c=tx1
	      s=ty1
	      tl=tl1
	    else
	      i1=iloc
	      i2=ir
	      gapn=-(xs-xc)*ty2+(ys-yc)*tx2
	      a0=alpha2
	      c=tx2
	      s=ty2
	      tl=tl2
	    endif
	  elseif(iloc .eq. il .or. iloc .eq. ir) then
	    tx=xr-xl
	    ty=yr-yl
	    tl=sqrt(tx*tx+ty*ty)
	    tx=tx/tl
	    ty=ty/tl
	    alpha=((xs-xl)*tx+(ys-yl)*ty)/tl
	    if (alpha .lt. -1.d-4 ) then
	      iwa(3,i)=il
	      iwa(4,i)=ir
	      wa(6,i)=0.0
	      wa(4,i)=0.0
	      goto 280
	    elseif (alpha .gt. 1.0001) then
	      iwa(3,i)=il
	      iwa(4,i)=ir
	      wa(6,i)=1.0
	      wa(4,i)=0.0
	      goto 280
	    else
	      i1=il
	      i2=ir
	      gapn=-(xs-xl)*ty+(ys-yl)*tx
	      a0=alpha
	      c=tx
	      s=ty
	    endif
	  endif
	  if(naxi.eq.0) then
	    xkn=stfmn(i1)
	    xkt=stfmt(i1)
	  elseif(naxi.eq.1) then
	    xkn = stfmn(i1)*abs(xl+xr)/2.d0
	    xkt = stfmt(i1)*abs(xl+xr)/2.d0
	  endif
	endif
c
c.... sliding, no release of contact node
c
	  if	 (iopt .eq. 1) then
	    ctr(5,i) = gapn/tl
	    forn     = xk*gapn
c.... penalty update
	    if (ifpa .le. 0) then
	      ctr(7,i) = forn
c.... augmented lagrangian update
	    else
	      ctr(9,i) = ctr(9,i) + forn
	      ctr(7,i) = ctr(9,i) + forn
	    endif
c
c.... tied sliding
c
	  elseif (iopt .eq. 2) then
c	   a0 = scoor(i)
	  a1 = 1.0d0-a0
	  i1 = iwa(1,i)
	  i2 = iwa(2,i)
	  x1 = xdm(i1)
	  x2 = xdm(i2)
	  y1 = ydm(i1)
	  y2 = ydm(i2)
	  tx = x2 - x1
	  ty = y2 - y1
	  tl = sqrt(tx*tx + ty*ty)
	  xc = a1*x1+a0*x2
	  yc = a1*y1+a0*y2
	  ax = xs-xc
	  ay = ys-yc
	  if(naxi.eq.0) then
	    xk = stfmn(i1)*tl
	  elseif(naxi.eq.1) then
	    xk = stfmn(i1)*tl*abs(x2+x1)/2.d0
	  endif
	  itc(4,i) =  1
	  c   = 1.0d0
	  s   = 0.0d0
	  ctr(5,i) =  ay/tl
	  ctr(6,i) = -ax/tl
	  fy = ay*xk
	  fx = ax*xk
c.... penalty update
	  if (ifpa .le. 0) then
	    ctr(7,i) =	fy
	    ctr(8,i) = -fx
c.... augmented lagrangian update
	  else
	    ctr( 9,i) = ctr( 9,i) + fy
	    ctr(10,i) = ctr(10,i) - fx
	    ctr( 7,i) = ctr( 9,i) + fy
	    ctr( 8,i) = ctr(10,i) - fx
	  endif
c
c.... sliding/release or sliding/friction
c
	  elseif(iopt .ge. 3) then
	   if(i2 .ge. 0) then
c.... segment geometry
	    ctr(11,i) = gapn/tl
	    forn     = gapn*xkn
c.... penalty update
	    if(ifpa.le.0) then
	      if (gapn .gt. 0.0) then
		iwa(3,i)=i1
		iwa(4,i)=i2
		wa(6,i)=a0
		wa(4,i)=0.0
		goto 280
	      else
		ctr(5,i)=forn/tl
		ctr(7,i) = forn
	      endif
c.... augmented lagrangian update
	    elseif(forn .le. ctr(9,i)) then
	      if(ifpa .eq. 1)ctr(9,i)=ctr(9,i)-forn
	      ctr(7,i)=forn-ctr(9,i)
	      ctr(5,i)=ctr(7,i)/tl
	    else
	      if(ifpa .eq. 1)ctr(9,i) = 0.0
	      wa(4,i)=0.0
	      wa(6,i)=a0
	      iwa(3,i)=i1
	      iwa(4,i)=i2
	      go to 280
	    endif
	   else
c.... corner geometry
	    forn=gapn*xkn
c.... penalty update
	    if(ifpa .le. 0)then
	      if(gapn .gt. 0.0)then
		iwa(3,i)=i1
		iwa(4,i)=i2
		wa(4,i)=0.0
		goto 280
	      else
		ctr(7,i)=forn
		ctr(5,i)=0.0
	      endif
c.... augmented Lagrangian update
	    else if (forn .le. ctr(9,i)) then
	      if(ifpa .eq. 1)ctr(9,i)=ctr(9,i)-forn
	      ctr(7,i)=forn-ctr(9,i)
	      ctr(5,i)=ctr(9,i)/gapn
	    else
	      if(ifpa .eq. 1)ctr(9,i)=0.0
	      wa(4,i)=0.0
	      iwa(3,i)=i1
	      iwa(4,i)=i2
	      goto 280
	    endif
	   endif
c
c     frictional forces-slideline type 4
c
	if(iopt .ne. 4)goto 260
	if(fric .eq. 0.0)goto 260
c.......  calculate tangential gap
	  b0=wa(3,i)
	  gapt=0.0
	  nseg=i1-iwa(1,i)
	  if(nseg .lt. 0)then
	    if(i2 .lt. 0 .and. iwa(2,i) .lt. 0)then
	      do 150, k=1,-nseg
	      gapt=gapt-slen(i1+k-1)
  150	      continue
	    elseif(i2 .gt. 0 .and. iwa(2,i) .lt. 0) then
	      gapt=-(1.-a0)*slen(i1)
	      do 155, k=1,-nseg-1
	      gapt=gapt-slen(i1+k)
  155	      continue
	    elseif(i2 .lt. 0 .and. iwa(2,i) .gt. 0)then
	      do 160, k=1,-nseg
	      gapt=gapt-slen(i1+k-1)
  160	      continue
	      gapt=gapt-b0*slen(iwa(1,i))
	    else
	      gapt=-(1.-a0)*slen(i1)
	      do 165, k=2,-nseg
	      gapt=gapt-slen(i1+k-1)
  165	      continue
	      gapt=gapt-b0*slen(iwa(1,i))
	    endif
	  else if (nseg .eq. 0) then
	    if(iwa(2,i) .lt. 0) b0=0.0
	    if(i2 .lt. 0) a0=0.0
	    gapt=(a0-b0)*slen(i1)
	  else
	    if (i2 .lt. 0 .and. iwa(2,i) .lt. 0) then
	      do 170, k=1,nseg
	      gapt=gapt+slen(iwa(1,i)+k-1)
  170	      continue
	    else if (i2 .gt. 0 .and. iwa(2,i) .lt. 0)then
	      do 175, k=1,nseg
	      gapt=gapt+slen(iwa(1,i)+k-1)
  175	      continue
	      gapt=gapt+a0*slen(i1)
	    else if(i2 .lt. 0 .and. iwa(2,i) .gt. 0)then
	      gapt=(1.-b0)*slen(iwa(1,i))
	      do 180, k=1,nseg-1
	      gapt=gapt+slen(iwa(1,i)+k)
  180	      continue
	    else
	      gapt=(1.-b0)*slen(iwa(1,i))
	      do 185, k=1,nseg-1
	      gapt=gapt+slen(iwa(1,i)+k)
  185	      continue
	      gapt=gapt+a0*slen(i1)
	    endif
	  endif
	endif
c
c.... penalty and nonsymmetric augmented lagrangians
c
	if(ifpa .eq. 0 .or. ifpa .eq. 3) then
	  fortr=xkt*gapt+wa(1,i)
	  if(ifpa .eq. 3) fortr=fortr+ctr(12,i)
c
c compute slip condition
c
	  formax=fric*ctr(7,i)
	  ftrial=abs(fortr)+formax
c check condition, do updates accordingly
c
	if(i2 .gt. 0)then
c segment geometry
	  if(ftrial .le. 0.0 .or. itc(4,i) .eq. 1) then
	    itc(4,i)=1
	    ctr(8,i)=fortr
	    ctr(6,i)=ctr(8,i)/tl
	  else
	    itc(4,i)=-1
	    ctr(8,i)=-formax*fortr/abs(fortr)
	    ctr(6,i)=ctr(8,i)/tl
	  endif
	  ctr(10,i)=slen(i1)/tl
	  wa(4,i)=ctr(8,i)
	  wa(6,i)=a0
	else
c corner geometry
	  if (ftrial .le. 0.0 .or. itc(4,i) .eq. 1)then
	     itc(4,i)=1
	     ctr(8,i)=fortr
	     ctr(6,i)=ctr(8,i)/gapn
	  else
	     itc(4,i)=-1
	     ctr(8,i)=-formax*fortr/abs(fortr)
	     ctr(6,i)=ctr(8,i)/gapn
	  endif
	     wa(4,i)=ctr(8,i)
	endif
	else
c symmetric  augmented lagrangians
	  if(i2 .gt. 0) then
c segment geometry
	    if(itc(4,i) .eq. 1) then
	       ctr(8,i)=xkt*gapt+ctr(12,i)+wa(1,i)
	       ctr(6,i)=ctr(8,i)/tl
	    else
	       ctr(8,i)=ctr(12,i)+wa(1,i)
	       ctr(6,i)=ctr(8,i)/tl
	    endif
	    ctr(10,i)=slen(i1)/tl
	    wa(4,i)=ctr(8,i)
	    wa(6,i)=a0
	  else
c corner geometry
	    if(itc(4,i) .eq. 1)then
	      ctr(8,i)=xkt*gapt+ctr(12,i)+wa(1,i)
	      ctr(6,i)=ctr(8,i)/gapn
	    else
	      ctr(8,i)=ctr(12,i)+wa(1,i)
	      ctr(6,i)=ctr(8,i)/gapn
	    endif
	    wa(4,i)=ctr(8,i)
	  endif
	endif
c
c.... set parameters for contact tangent/residual
c
  260	  continue
	  iwa(3,i)  =i1
	  iwa(4,i)  =i2
	  ctr(1,i)  = xkn
	  if(ifpa.gt.0 .and. i2.lt.0) ctr(1,i) = ctr(1,i)-ctr(9,i)/gapn
	  ctr(2,i)  = s
	  ctr(3,i)  = c
	  ctr(4,i)  = a0
	  ctr(13,i) = xkt
      nde1 = nde
      if(nde1.eq.0) nde1=1
    	if (i2 .ge. 0) then
c
c.... set up connection array for 3-node contact element
c
	    itc(1,i) = nsv(i)
  	  itc(2,i) = msr(i1)
	    itc(3,i) = msr(i2)
  	  call rstprf(jd,idl,id,itc(1,i),edge1,edge2,edge3,
     1		      nde1,ndf,3,3,0,1,1,.false.)
    	else
c
c.... set up connection array for 2-node contact element
c
	    itc(1,i)=nsv(i)
	    itc(2,i)=msr(i1)
	    itc(3,i)=-1
	    call rstprf(jd,idl,id,itc(1,i),edge1,edge2,edge3,
     1		      nde1,ndf,2,2,0,1,1,.false.)
	endif
  280 continue
 2000 format('*** subsl1 *** Only one master node - logic fails')
c 2001 format(2x,'zero master segment length for slave node',i4)
      return
      end
c
      subroutine subupd(nsn,was,iwas,ctrs,itcs)
      implicit double precision (a-h,o-z)
      dimension was(6,*),iwas(4,*),ctrs(14,*),itcs(4,*)
c.... move the time data for the two arrays
      do 10 n = 1,nsn
	was(1,n)  = was(4,n)
	was(2,n)  = was(5,n)
	was(3,n)  = was(6,n)
	iwas(1,n) = iwas(3,n)
	iwas(2,n) = iwas(4,n)
c.... do other updates
	itcs(4,n) = 2
c  (to ensure that the first iteration in a step is done with all nodes
c     sticking)
	ctrs(12,n) = 0.0
10    continue
      return
      end
