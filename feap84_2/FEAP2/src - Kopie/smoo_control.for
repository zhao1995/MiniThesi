      subroutine smoo_control
     +     (lct,ct,ndf,ndm,nen1,nst,nneq,ne,prt,rlnew,
     +        timold,kflag,l,m)
c......................................................................
c.... control type of smoothing  plane/shell
c     DA Peter Fellmoser
c......................................................................
      USE cdata
      USE iofile
      USE mdata
      USE plong
      implicit double precision (a-h,o-z)
      common /mopt1/ fer,ngr,iud,nmi,idim,ivnn,iern,izv,ib
      common /mopt2/ lgr,lmi,lmj,lwr
      common /mdat3/  n9a,n9b,n9c,n9d,n10a,n11c,n11d,n14a,numelo,numnpo
      common /msm1/ lop,lwe,lgl,lsw,lsi
      common /msm2/ iln,ila,ile,mel,ngl
      common /adap2/  lad
      character*4 lct(*)
      logical lgr,lmi,lmj,lwr,lad,lop,lwe,lgl,lsw,lsi
czr   dimension ct(2,*),m(*),dm(*)
      dimension ct(3,*),m(*)
c.... determine node-element connectivities for adaptivity and others
C.... felder m(ne) und m(nnum) werden nur temp. benoetigt
      call pseta(n9a,numel*(nen+1),1,tflb,'SMOOTH-ike')    ! --> ike
      call pseta(n9b,numel*(nen+1),1,tflb,'SMOOTH-ikz')    ! --> ikz
      ne = kmax
      nnum = Ne + numnp + MOD(numnp,IPR)
      call elnode (econ,m(n9a),m(n9b),nen1,numnp,numel,nmat,
     +             m(ne),m(nnum))
c
cpf.. transfer to correct process
cpf    go to (1,2), j
cpfc..  mesh optimization based on minimation of FE-error
cpf1   continue
cpf    if (lgr) call psetm(ngr,ne,numnp*(ipr*(2+ndm)),lgr)
cpf    call moptim(m(n8),lct,ct,ndf,ndm,nen1,nst,nneq,ne,prt,rlnew,
cpf   +      timold,kflag,l,m(ngr),m(ngr+numnp*ipr*ndm),
cpf   +        m(ngr+numnp*(ipr*(ndm+1))),m,dm)
cpf    return
cpfc..  geometrical mesh optimization
cpf2   continue

cpf   ivn =ct(1,l)
cpf   itnn=ct(2,l)
      ict31 = ct(3,l)/100
      ict32 = ct(3,l) - ict31*100
      if (abs(ict32).lt.0.00001) then
       ict32=1
      endif
      if (abs(ict31).lt.0.00001) then
       ict31=1
      endif
cpf   if (.not.lad) mel = 1
      if (ndm.eq.2) then
czr   call msmoot(m(n9a),m(n9b),m(n8),m(n9),ict31,ict32,lct(l),ivn,itnn,ndm,ndf,
czr  +      m(mel),m)
      call msmoot(m(n9a),m(n9b),coor,econ,ict31,ict32,lct,ct,
     +      ivn,itnn,ndm,ndf,m(mel),l,m)
      elseif (ndm.eq.3) then
cpf   write(*,1000) ndm
cpf   write(iow,1000) ndm
czr   if (lgl) call psetm(ngl,ne,numnp,lgl)
czr   call msmosh(m(n9a),m(n9b),m(n8),m(n9),lct(l),ivn,itnn,ndm,ndf,
c.the+    m(mel),m(ngl),m)
czr   +      m(n9d),m(ngl),m)
      call msmosh(m(n9a),m(n9b),coor,econ,lct,ivn,itnn,ndm,ndf,
     +            m(n9d),m(ngl),m,ict31,ict32,ct,l)
      else
        write(iow,1000) ndm
      endif
1000  format('mesh smoothing algorithm not implemented for ndm= ',i1)
      end
c -----------------------------------------------------------------------
czr   subroutine msmoot(ike,ikez,x,ix,lct,ivn,itnn,ndm,ndf,iael0,m)
      subroutine msmoot(ike,ikez,x,ix,ict31,ict32,lct,ct,
     +                    ivn,itnn,ndm,ndf,iael0,l,m)
c
c.... This macro command optimizes geometrically an given mesh
C
C     call: msmo,i1i2,c1,c2
C
C     i1 : weight factor for equality of sidelengths
C          (first two characters of character-chain) = integer i1
C     i2 : weight factor for 90-degree-angles
C          (third and fourth character of character-shain) = integer i2
C          (for instance i1i2 = 0210 means : i1=2 and i2=10
C     c1 = j1*1000+j2
C          j1 : control parameter : (loop,msmo,.. as adaptive loop necessary)
C            j1=1 : execute mesh smoothing only at the last adaptive step
C            j1=2 : execute mesh smoothing at first and last adaptive step
C            j1=3 : execute mesh smoothing at all adaptive refinement steps
C          j2 : number of vertices of the fe-model ,
C               (this vertices must be the first one at the nodal-point list)
C     c2 = k1*1000000 + k2*1000 + k3
C          k1 : number of newton steps in every local iteration
C          k2 : function choosing parameter (1-7 available)
C          k3 : number of iteration steps over all elements
C               (useful 1,2 or 3)
C
C
C                           /a**2 - b**2\ 2       /a**2 * b**2    \ 2
C             f1 = f3 = w1*| ----------- |  + w2*| ----------- - 1 |
C                           \(a** + b**2/         \ (a x b)**2    /
C
C             f2 = w1*(a**2-b**2)**2 + w2*(c**2-a**2-b**2)**2
C
C                      / a**2-b**2 \ 2              /   a*b   \
C             f4 = w1*|  ---------  |  + w2*arccos |  -------  |
C                      \ a**2+b**2 /                \ |a|*|b| /
C
C             f5 = (a**2-rk*c**2)**2 + (b**2-rk*c**2)**2
C                   rk = (sin(w/2))**2 , w = averanged angle
C
C             f6 = (a**2-rk*c**2)**2 + (b**2-rk*c**2)**2
C                   rk = (sin(w/2))**2 , w = desired angle (90 degrees)
C
C                      / a**2-b**2 \ 2              /   a*b     Pi \2
C             f7 = w1*|  ---------  |  + w2*arccos |  ------- - --  |
C                      \ a**2+b**2 /                \ |a|*|b|   2  /
C
C                      / a-b \ 2              /   a*b     Pi \ 2
C             f8 = w1*|  ---  |  + w2*arccos |  ------- - --  |
C                      \ a+b /                \ |a|*|b|   2  /
C
C                            4
C                     ( a-b )               /   a*b     Pi \ 2
C             f9 = w1* ----- 2 + w2*arccos |  ------- - --  |
C                     ( a*b )               \ |a|*|b|   2  /
C
C     N.b.: All boundaries must be defined by curves
C
      USE bdata
      USE cdata
      USE comfil
      USE iofile
      USE mdata
      USE mdat2
      implicit double precision (a-h,o-z)
czr   common /ycur1/ cpar(100,8),nrt(100,3),ic
      common /curvedat/  cpar(20,8),nrt(20,3),ic,nbe,nn3
      common /msm1/ lop,lwe,lgl,lsw,lsi
      common /msm2/ iln,ila,ile,mel,ngl
      common /fvc/ w1,w2,rhv(4,8),ieh
      common /adap2/  lad
      character*4 lct(*)
czr   character*4 o,head
czr   character*2 llct
      logical lop,lbn,lou,lf,lad,lwe,lne,lb,lsg,log,lfp,lgl
     +     ,lsw,lsi,lshel
     +     ,pcomp
      dimension ike(*),ikez(*),x(ndm,*),ix(nen+4,*),xm(3),m(*),der(3)
czr  +      xh(3),iael0(*),zfvh(2),fv1(2),fv2(2),grad(3)
     +      ,xh(3),iael0(*),zfvh(2),fv1(2),fv2(2)
     +      ,ct(3,*)
      data tol /1.0d-4/
c
c     write(*,113) 'ike........',(ike(k),k=1,10)
c     write(*,113) 'ikez.......',(ikez(k),k=1,10)
c     do i=1,nen+4
c      write(*,113) 'ix.........',(ix(i,k),k=1,10)
c     enddo
c     write(*,114) 'x..........',(x(1,k),k=1,10)
c     write(*,114) 'x..........',(x(2,k),k=1,10)
c     write(*,*) 'lct........',lct
c     write(*,*) 'ivn........',ivn
c     write(*,*) 'itnn.......',itnn
c     write(*,*) 'ndm........',ndm
c     write(*,*) 'ndf........',ndf
c     write(*,*) 'iae10......',iae10
c     write(*,*) 'm..........',m
c113   format(a,10(i3))
c114   format(a,10(f5.2,2x))
c.... log = only gradient iteration
      log=.false.
c.... lsg = switch to gradient iteration
      lsg=.false.
c.... w3 ... third weigth factor
cpf   w3=ivn/1000000
cpf   ivnn=ivn-w3*1000000
cpf   isq=ivnn/1000
cpf>
      isq=1
cpf<
cpf   iv=ivnn-isq*1000
cpf>
      iv=ct(1,l)-1
      if (abs(ct(1,l)).lt.tol)then
       iv=0
      endif
cpf<
      if(((isq.eq.1).and.(ila.eq.ile)).or.
     +   ((isq.eq.2).and.((ila.eq.1).or.(ila.eq.ile))).or.
     +    (isq.eq.3)) then
        lne = .not.((ila.eq.ile).or.(ila.eq.1))
        lou=.false.
czr     if ((fsav(1:3).ne.'nul').and.(.not.(lad))) then
czr       open (3,file=fsav,status='unknown',err=900)
czr       if (lop) rewind(3)
czr       lop=.false.
czr900    continue
czr       lou=.true.
czr       call rest(bang,m(n8),m(n9),m(n10),m(n6),m(n7),
czr   +              ndm,ndf,nen+4)
czr      endif
c
      rkap=0.7d0
      delta=0.1d0
      alp=0.5d0
czr     if (ndm.ne.2) then
czr       write(*,2000)
czr       write(iow,2000)
czr     endif
        lf=.true.
c..... w1 ... first weigth factor
czr     llct=lct(1:2)
czr     read(llct,4000,err=110) i
czr     if (i.ne.0) lf=.false.
czr     w1=dble(i)
        if (ict31.ne.0) lf=.false.
        w1=dble(ict31)
c.... w2 ... second weight factor
czr     llct=lct(3:4)
czr     read(llct,4000,err=110) i
czr     if (i.ne.0) lf=.false.
czr     w2=dble(i)
        if (ict32.ne.0) lf=.false.
      w2=dble(ict32)
cpf     mnv=itnn/1000000
cpf>
      mnv=1
cpf<
      if (mnv.eq.0) then
        write(*,9100)
        write(iow,9100)
        stop
      endif
c....
cpf   itnnn=itnn-mnv*1000000
cpf   iwf=itnnn/1000
cpf>
      iwf=10
cpf<
      if (iwf.eq.0) then
        write(*,9000)
        write(iow,9000)
        stop
      endif
cpf>
      if (pcomp(lct(l),'oesp',4) .or. pcomp(lct(l),'wesp',4)  ) then
        continue
      else
        write(*,*)'Enter: oesp wesp '
        goto 998
      endif
c.... number of iterations
cpf   itn=itnnn-iwf*1000
cpf>
      itn=ct(2,l)
      if (itn.eq.0) then
        itn=1
      endif
cpf<
      it=w1
      write(iow,8000) ict31,ict32,itn,mnv
      write(*,8000) ict31,ict32,itn,mnv
c.... loop - number of iterations
      do 100 it=1,itn
        lb = (mod(it,2).eq.0)
cl      do 165 i=iv+1,numnp
        do 165 iii=iv+1,numnp
          if (lb) then
            i=numnp-iii+iv+1
          else
            i=iii
          endif
cl end
          xm(1)=x(1,i)
          xm(2)=x(2,i)
          ieh=ikez(i+1)-ikez(i)
          if (ieh.gt.8) then
            write (iow,6000) i
            write (*,6000) i
            stop
          endif
c.... loop - neh: number of elements connected to node
          do 330 j=1,ieh
            ihv=ike(ikez(i)+j-1)
c.... ihv = node number
            if (lad) then
              if (lne.and.(iael0(ihv).lt.0)) goto 165
            endif
            do 660 k=1,nen
              if(ix(k,ihv).eq.i) then
                if(k.eq.1) then
                  il=ix(nen,ihv)
                  iu=ix(2,ihv)
                elseif(k.eq.nen) then
                  il=ix(nen-1,ihv)
                  iu=ix(1,ihv)
                else
                  il=ix(k-1,ihv)
                  iu=ix(k+1,ihv)
                endif
              endif
660         continue
            rhv(1,j) = x(1,il)
            rhv(2,j) = x(2,il)
            rhv(3,j) = x(1,iu)
            rhv(4,j) = x(2,iu)
330       continue
c
c.... distinguish boundary node and inner node
c
cpf>
      nc=0
cpf<
      if(ic.eq.0)
     +  write(*,*)'All boundaries must be defined by curves'
        do 825 jj=1,ic-1
            y=xm(2)
            diff=0.0d0
czr         call curve(xm(1),xm(2),0.0d0,m,jj,diff,nr)
            call curve(xm(1),xm(2),0.0d0,m,jj,diff)
czr         if(dabs(y-xm(2)).lt.tol) then
            if(dabs(diff).lt.1e-2) then
cpf>
              nc=nc+1
              njj=jj
czr>
              nr=jj
            endif
825     continue
          if (nc.eq.1)then
             lbn=.true.
          if (log) goto 700
            goto 495
          elseif (nc.ge.2) then
            goto 165
          endif
cpf<
czr          goto 500         !zr boundary nodes are not treated
cpf          endif
cpf825         continue
C.... case of inner node
      lbn=.false.
      xk=xm(1)
      yk=xm(2)
c.... simple middle value method
      if (iwf.ge.10) goto 400
c.... newton method
      if (log) goto 700
      do 200 nv=1,mnv
c.... check, if functional value decreases
      zfo  = f(xk,yk,iwf,1)
      f1  = f(xk,yk,iwf,2)
      f2  = f(xk,yk,iwf,3)
      a11  = f(xk,yk,iwf,4)
      a12  = f(xk,yk,iwf,5)
      a22  = f(xk,yk,iwf,6)
      det = a11*a22-a12**2
      rdet = abs(a11*a22)+a12**2
      if((abs(det).lt.1.0d-6*rdet).or.(abs(det).lt.1.0d-50))then
        write(*,3000) i
        if (lsg) goto  700
          goto 166
        endif
        rh1 = (f2*a12 - f1*a22)/det
        rh2 = (f1*a12 - f2*a11)/det
c       xkm2= xkm1
c       ykm2= ykm1
        xkm1= xk
        ykm1= yk
        xk  = xk + rh1
        yk  = yk + rh2
        zf  = f(xk,yk,iwf,1)
        if(zfo.lt.zf) then
          xk=xkm1
          yk=ykm1
c.... switch to gradient iteration
      write(*,5000) i,nv
      write(iow,5000) i,nv
      if (lsg) goto 700
        goto 166
      endif
      if (dsqrt(rh1**2+rh2**2).lt.1.0d-6*dsqrt(xk**2+yk**2))
     +  goto 500
200   continue
      goto 500
495   continue
c.... case of boundary node
C
cpf >
      if (pcomp(lct(l),'oesp',4) ) then
cpf     write(*,*)'oesp.....!'
        goto 165
      endif
cpf <
      xk=xm(1)
      yk=xm(2)
      lbn=.true.
cpf>
      log=.true.
cpf<
C
      if (log) goto 700
400         continue
      if ((((iwf.eq.10).or.(iwf.eq.12)).and.(.not.lbn))
     +      .or.(iwf.eq.11)) then
c.... simple middle value method / inner node
      xkm1=0.0d0
      ykm1=0.0d0
      do 300 ij=1,ieh
        x1=rhv(1,ij)
        y1=rhv(2,ij)
        x2=rhv(3,ij)
        y2=rhv(4,ij)
        rh1=(x1+x2-y1+y2)/2.0d0
        rh2=(y1+y2+x1-x2)/2.0d0
        f1=(x1+x2+y1-y2)/2.0d0
        f2=(y1+y2-x1+x2)/2.0d0
        if(dsqrt((rh1-xk)**2+(rh2-yk)**2).lt.dsqrt((f1-xk)**2+
     +       (f2-yk)**2)) then
          xkm1=xkm1+rh1
          ykm1=ykm1+rh2
        else
          xkm1=xkm1+f1
          ykm1=ykm1+f2
        endif
300   continue
      xk=xkm1/ieh
      yk=ykm1/ieh
cpf   if (lbn) call rcurve(xk,yk,jj,m)
      if (lbn) call rcurve(xk,yk,njj,m)
      x(1,i)=xk
      x(2,i)=yk
      goto 165
      elseif ((iwf.eq.12).and.(lbn)) then
        if (ieh.eq.2) then
          xk=rhv(1,1)
          yk=rhv(2,1)
czr         call curve(xm(1),xm(2),0.0d0,m,jj,diff,nr)
            call curve(xm(1),xm(2),0.0d0,m,jj,diff)
          if (abs(yk-rhv(2,1)).lt.1.0e-6) then
            xk=rhv(3,1)
            yk=rhv(4,1)
czr         call curve (xk,yk,dummy,m,jj,diff,nr)
            call curve (xk,yk,dummy,m,jj,diff)
            if (abs(yk-rhv(4,1)).lt.1.0e-6) then
              write(*,3100) i,jj
              write(iow,3100) i,jj
              stop
            endif
            yk=rhv(4,1)
          endif
          call rcurve(xk,yk,jj,m)
          goto 166
        else
          lfp=.true.
          do 310 ij=1,ieh
            do 320 ji=1,2
              if (lfp) then
                x1=rhv(2*(ji-1)+1,ij)
                y1=rhv(2*ji,ij)
czr             call curve (x1,y1,dummy,m,jj,diff,nr)
                call curve (x1,y1,dummy,m,jj,diff)
                if (abs(y1-rhv(2*ji,ij)).lt.1.0d-6) lfp=.false.
              else
                x2=rhv(2*(ji-1)+1,ij)
                y2=rhv(2*ji,ij)
                if ((dabs(x1-x2).lt.1.0d-50)
     +               .and.(dabs(y1-y2).lt.1.0d-50)) then
                   write(*,3300) jj,i
                   write(iow,3300) jj,i
                   stop
                endif
czr             call curve (x2,y2,dummy,m,jj,diff,nr)
                call curve (x2,y2,dummy,m,jj,diff)
                if (abs(y2-rhv(2*ji,ij)).lt.1.0e-6) then
                  xk=(x1+x2)/2.0d0
                  yk=(y1+y2)/2.0d0
                  call rcurve(xk,yk,jj,m)
                  goto 166
                endif
              endif
320         continue
310       continue
        endif
        write(*,3200) i,jj
        write(iow,3200) i,jj
        stop
      elseif (iwf.eq.13) then
c.... Gauss--Newton--method
      if (log) goto 700
      do 600 nv=1,mnv
c.... check, if function value decreases
        zfvo=0.0d0
        f1=0.0d0
        f2=0.0d0
        if (lbn) then
          xh(1)=xk
          xh(2)=yk
czr       call curved (xh,der,jj)
          write(*,*)'sr curved not implemented SR MSMOOT'
          do 140 ih=1,ieh
            x1=rhv(1,ih)
            y1=rhv(2,ih)
            x2=rhv(3,ih)
            y2=rhv(4,ih)
            call g(xk,yk,x1,y1,x2,y2,iwf,2,fv1)
            call g(xk,yk,x1,y1,x2,y2,iwf,3,fv2)
            call g(xk,yk,x1,y1,x2,y2,iwf,1,zfvh)
            zfvo=zfvo+zfvh(1)**2+zfvh(2)**2
            rh1=fv1(1)*der(1)+fv2(1)*der(2)
            rh2=fv1(2)*der(1)+fv2(2)*der(2)
            f1=f1+zfvh(1)*rh1+zfvh(2)*rh2
            f2=f2+rh1**2+rh2**2
140         continue
            a11=f1/f2
            rh1=a11*der(1)
            rh2=a11*der(2)
        else
            a11=0.0d0
            a12=0.0d0
            a22=0.0d0
            do 120 ih=1,ieh
              x1=rhv(1,ih)
              y1=rhv(2,ih)
              x2=rhv(3,ih)
              y2=rhv(4,ih)
              call g(xk,yk,x1,y1,x2,y2,iwf,1,zfvh)
              zfvo=zfvo+zfvh(1)**2+zfvh(2)**2
              call g(xk,yk,x1,y1,x2,y2,iwf,2,fv1)
              f1 = f1 + fv1(1)*zfvh(1)+fv1(2)*zfvh(2)
              call g(xk,yk,x1,y1,x2,y2,iwf,3,fv2)
              f2 = f2 + fv2(1)*zfvh(1)+fv2(2)*zfvh(2)
              a11  = a11 + fv1(1)**2+fv1(2)**2
              a12  = a12 + fv1(1)*fv2(1)+fv1(2)*fv2(2)
              a22  = a22 + fv2(1)**2+fv2(2)**2
120         continue
            det = a11*a22-a12**2
            rdet = abs(a11*a22)+a12**2
            if((abs(det).lt.1.d-6*rdet).or.(abs(det).lt.1.d-50))
     +          then
              write(*,3000) i
              if (lsg) goto  700
              goto 166
            endif
            rh1 = (f1*a22 - f2*a12)/det
            rh2 = (f2*a11 - f1*a12)/det
        endif
        xkm1= xk
        ykm1= yk
        xk  = xk - rh1
        yk  = yk - rh2
        if (lbn) call rcurve(xk,yk,jj,m)
          zfv=0.0d0
          do 130 ih=1,ieh
            x1=rhv(1,ih)
            y1=rhv(2,ih)
            x2=rhv(3,ih)
            y2=rhv(4,ih)
            call g(xk,yk,x1,y1,x2,y2,iwf,1,zfvh)
            zfv=zfv+zfvh(1)**2+zfvh(2)**2
130       continue
          if (zfvo.lt.zfv) then
            xk=xkm1
            yk=ykm1
c.... switch to gradient iteration
            write(*,5000) i,nv
            write(iow,5000) i,nv
            if (lsg) goto 700
              goto 166
            endif
            if
     +      (dsqrt(rh1**2+rh2**2).lt.1.0d-6*dsqrt(xk**2+yk**2))
     +       goto 500
600   continue
      goto 500
      endif
c.... case simple method and boundary node
      xk=xm(1)
      yk=xm(2)
      do 800 nv=1,mnv
        xh(1)=xk
        xh(2)=yk
c       call curveg ( xk,yk,z,jj,grad )
czr     call curved (xh,der,jj)
        write(*,*)'sr curved not implemented SR MSMOOT'
c....   der(1)=grad(2)
c....   der(2)=grad(1)
c....   der(3)=0
c....   call vnorm(der,dummy)
        zfo  = f(xk,yk,iwf,1)
        f1  = f(xk,yk,iwf,2)
        f2  = f(xk,yk,iwf,3)
        a11  = f(xk,yk,iwf,4)
        a12  = f(xk,yk,iwf,5)
        a22  = f(xk,yk,iwf,6)
c       a11 = a11*der(1)+a12*der(2)
c       a12 = a12*der(1)+a22*der(2)
cric c  a21 = der(2)
c       call curveg ( xk,yk,z,jj,grad )
c       a21 = grad(1)
c       a22 = grad(2)
cric c end(newb)            a22 = -der(1)
c       f1  = f1*der(1) + f2*der(2)
cric c  f2  = 0.0d0
c       call curve ( xk,yk,z,m,jj,di)
cric    f2 = -di
c       f2 = di
cric end(rnv)
c ric c end(newb)
c       det = a11*a22-a12*a21
c       rdet = abs(a11*a22)+abs(a12*a21)
        rdet = abs(a11*der(1)**2)+abs(2.0d0*a12*der(1)*der(2))+abs
     +            (a22*der(2)**2)
        fal = a11*der(1)**2+2.0d0*a12*der(1)*der(2)+a22*der(2)**2
c       if (abs(det).lt.10d-6*rdet) then
        if((abs(fal).lt.1.0d-6*rdet).or.(abs(fal).lt.1.0d-50))then
          write(*,3000) i
          write(iow,3000) i
          goto 800
        endif
c       rh1 = (f2*a12 - f1*a22)/det
c       rh2 = (f1*a21 - f2*a11)/det
c       xkm2= xkm1
c       ykm2= ykm1
        al=-(f1*der(1)+f2*der(2))/fal
        xkm1=xk
        ykm1=yk
        xk  = xk + al*der(1)
        yk  = yk + al*der(2)
        call rcurve(xk,yk,jj,m)
        zf  = f(xk,yk,iwf,1)
        if(zfo.lt.zf) then
          xk=xkm1
          yk=ykm1
          write(*,5005) i,nv
          write(iow,5005) i,nv
          if (lsg) goto 700
          goto 166
        endif
        if(dsqrt(rh1**2+rh2**2).lt.1.0d-6*dsqrt(xk**2+yk**2)) then
          goto 500
        endif
800   continue
500   continue
      goto 166
700   continue
c.... gradient method
c     xk=xm(1)
c     yk=xm(2)
cpf>
      xk=xm(1)
      yk=xm(2)
cpf<
C
      izv=0
      do 350 ita=1,mnv
C.... compute gradient
        zf1  = f(xk,yk,iwf,1)
        f1  = f(xk,yk,iwf,2)
        f2  = f(xk,yk,iwf,3)
c       rmas=dsqrt((xk-x1)**2+(yk-y1)**2)/dsqrt(f1**2+f2**2)/2.0d0
        rmas=10.0d10
C
        if (lbn) then
          xh(1)=xk
          xh(2)=yk
cpf       call curved (xh,der,jj)
cpf>
          xh(3)=0.0
          call curved (xh,der,njj,ityp,cpar)
cpf<
cpf       write(*,*)'sr curved not implemented SR MSMOOT!!!'
cric      f1 = f1*der(1)
          rab=f1*der(1)+f2*der(2)
          f1 = rab*der(1)
          f2 = rab*der(2)
cric end(ng)                f2 = f2*der(2)
        endif
C
c.... set local parameters for gradient algorithm and end-test
C
        gal=0.0d0
        dgxk2 = f1**2+f2**2
        if(ita.eq.1) ak=dgxk2/1.0d6
        if (dgxk2.lt.dmax1(ak,1.0d-30)) then
          goto 166
        elseif (ita.eq.1) then
          gam=abs(zf1/dgxk2/10.0d0)
        endif
c.... compute step-length by Goldstein-method
C
        rho=gofun(xk,yk,f1,f2,dgxk2,zf1,gam,rmas,izv,iwf)
111     continue
        if (rho.le.rkap) goto 555
        if (izv.gt.1) then
          write(iow,2009) gam
          write(*,2009) gam
          goto 701
        endif
        gam=gam/alp
        rho=gofun(xk,yk,f1,f2,dgxk2,zf1,gam,rmas,izv,iwf)
        goto 111
333     continue
        gam=(gal+gar)/2.0d0
        if (gam.lt.1.0d-40) then
cpf       write (iow,1000) i
cpf       write (*,1000) i
          goto 165
        endif
        rho=gofun(xk,yk,f1,f2,dgxk2,zf1,gam,rmas,izv,iwf)
555     continue
        if (rho.lt.delta) then
          gar=gam
          goto 333
        endif
c.... execute gradient step
C
701     continue
        xk=xk-gam*f1
        yk=yk-gam*f2
cpf.. Schnittpunkt Gerade und Kugel
        if ((ityp.eq.32).and.(lbn)) then
          b=(yk-cpar(njj,3))/(xk-cpar(njj,2))
          a=(cpar(njj,3)*xk-yk*cpar(njj,2))/(xk-cpar(njj,2))
          xd1=-1/(1+b**2)*(-cpar(njj,2)+b*a-b*cpar(njj,3))
     +        +sqrt((1/(1+b**2)*(-cpar(njj,2)+b*a-b*cpar(njj,3)
     +        ))**2-1/(1+b**2)*(cpar(njj,2)**2+a**2-2*a*
     +        cpar(njj,3)+cpar(njj,3)**2-cpar(njj,1)**2))
          xd2=-1/(1+b**2)*(-cpar(njj,2)+b*a-b*cpar(njj,3))
     +        -sqrt((1/(1+b**2)*(-cpar(njj,2)+b*a-b*cpar(njj,3)
     +        ))**2-1/(1+b**2)*(cpar(njj,2)**2+a**2-2*a*
     +        cpar(njj,3)+cpar(njj,3)**2-cpar(njj,1)**2))
          if(dabs(xd1-xk).lt.dabs(xd2-xk)) then
            xd=xd1
          else
            xd=xd2
          endif
          yd=b*xd+a
          xk=xd
          yk=yd
        endif
cpf.. Schnittpunkt Gerade und Kreis
        if ((ityp.eq.3).and.(lbn)) then
          b=(yk-cpar(njj,4))/(xk-cpar(njj,3))
          a=(cpar(njj,4)*xk-yk*cpar(njj,3))/(xk-cpar(njj,3))
          xd1=-1/(1+b**2)*(-cpar(njj,3)+b*a-b*cpar(njj,4))
     +        +sqrt((1/(1+b**2)*(-cpar(njj,3)+b*a-b*cpar(njj,4)
     +        ))**2-1/(1+b**2)*(cpar(njj,3)**2+a**2-2*a*
     +        cpar(njj,4)+cpar(njj,4)**2-cpar(njj,1)**2))
          xd2=-1/(1+b**2)*(-cpar(njj,3)+b*a-b*cpar(njj,4))
     +        -sqrt((1/(1+b**2)*(-cpar(njj,3)+b*a-b*cpar(njj,4)
     +        ))**2-1/(1+b**2)*(cpar(njj,3)**2+a**2-2*a*
     +        cpar(njj,4)+cpar(njj,4)**2-cpar(njj,1)**2))
          if(dabs(xd1-xk).lt.dabs(xd2-xk)) then
            xd=xd1
          else
            xd=xd2
          endif
            yd=b*xd+a
            xk=xd
            yk=yd
        endif
cpf<
cpf     if (lbn) call rcurve(xk,yk,jj,m)
        if (lbn) call rcurve(xk,yk,njj,m)
350     continue
166     continue
        if (lbn) then
cpf       call rcurve(xk,yk,jj,m)
          call rcurve(xk,yk,njj,m)
          do 370 ij=1,ieh
            x1=rhv(1,ij)
            y1=rhv(2,ij)
            x2=rhv(3,ij)
            y2=rhv(4,ij)
            if ((xk-x2)*(yk-y1)-(xk-x1)*(yk-y2).lt.0.0d0) then
cpf           write(*,3700) i
cpf           write(iow,3700) i
              goto 165
            endif
370       continue
        endif
        x(1,i)=xk
        x(2,i)=yk
165     continue
        if (lou) call rest
     +  (bang,coor,econ,gloa,edma,psid,ndm,ndf,nen+4)
100     continue
        if (lou) close(3)
      endif
      write(*,8800)
cryz
cpf>
      lsw=.false.
cpf<
      if (lsw) then
        lshel = .false.
czr     call quality (ike,ikez,x,ix,xm,ndm,lad,lne,lshel,isq,m)
        write(*,*)'sr quality not implemented SR MSMOOT'
        stop
      endif
cryz  end
998   continue
      return
czr110   call perror('msmo')
czr   stop
cpf1000  format('gradient iteration canceled at node',i5,' (gam to small)')
cpf2000  format
cpf     +('mesh optimization subroutine only for 2-d meshes available yet')
2009  format('gradient step executed with maximal allowed step length',
     +      g20.10)
3000  format('0 valued det. in newton-method (msmooth,',
     +        'node ',i5,' ',i3,'th iteration)')
3100  format('node',i5,' and further 2 nodes of 1 element on curve',i3)
3200  format('no complete set of boundary points find around node',i5,
     +        ' on curve',i3)
3300  format('2 same boundary points at curve',i3,' ; node=',i5)
cpf3700  format('neg. jacobi det. in element at node',i3,
cpf     +      ' ; reset to start value')
czr4000  format(i2)
5000  format('switched to gradient iteration at node',i5,
     +      ' ; newton step:',i4)
5005  format('switched to gradient iteration at boundary node',i5,
     +      ' ; newton step:',i4)
6000  format('error in macro command msmo : more then 8 elements at',
     +       ' node ',i5)
8000  format('macro msmo: function:  8    w1= ',i2,'  w2= ',i2,
     +       '  iterations: ',i3,' max. local steps: ',i3)
8800  format('end of macro msmo')
9000  format('new function-choosing-parameter: n2=fn*1000+itn',/,
     +       '(fn=number of function for minimum problem; itn=number ',
     +       'of global iterations)')
9100  format('new function-choosing-parameter: n2=mn*1000000+fn*1000',
     +        '+itn',/5x,'mn=maximal number of newton steps in local ',
     +        'iterations',/5x,'fn=number of function for minimum',
     +        ' problem ',/5x,'itn=number of global iterations',/)
      end
c......................................................................
      function gofun(xk,yk,f1,f2,dgxk2,zf1,gam,rmas,izv,iwf)
c
      implicit double precision (a-h,o-z)
c
      if (gam.gt.rmas) then
c       write(*,2008) gam,rmas
        gam=rmas
        izv=1
      endif
      gofun=(zf1-f(xk-gam*f1,yk-gam*f2,iwf,1))/(gam*dgxk2)
      if (gofun.ge.1.0d0) ib=1
ctestgrad      write(*,2000) gam,zf1,gofun
      return
c2000 format('gam = ',e13.6,'zf1 = ',e13.6,'   rho = ',e13.6)
c2008 format('steplength (',f13.6,') lowered to',f13.6)

      end
c......................................................................
      subroutine rcurve(xk,yk,jj,m)
c
      implicit double precision (a-h,o-z)
czr   common /ycur1/ cpar(100,8),nrt(100,3),ic
      common /curvedat/  cpar(20,8),nrt(20,3),ic,nbe,nn3
      dimension m(*)
      if (nrt(jj,2).eq.1) then
        x1 = cpar(jj,1)
        y1 = cpar(jj,2)
        x2 = cpar(jj,3)
        y2 = cpar(jj,4)
        d = x2-x1
        g = y2-y1
        izv=1
165     continue
        if (dabs(d).lt.1.0d-6*dabs(g)) then
          xk=x2-(y2-yk)*d/g
        else
          b = (y1*x2-y2*x1) / (x2-x1)
          a = (y2-y1) / (x2-x1)
          if (abs(a).lt.1.0d-10) then
            yk=b
          else
            yk= (a*xk+a**2*yk+b)/(a**2+1)
            xk= (yk-b)/a
          endif
        endif
        x1=xk
        y1=yk
        call curve (x1,y1,dummy,m,jj,diff)
        if(abs(y1-yk).gt.1.0d-3) then
          if (izv.lt.10) then
            goto 165
czr (never executed!)        izv=izv+1
          else
            write(*,*)'subroutine rcurve'
            stop
          endif
        endif
        return
      elseif (nrt(jj,2).eq.5) then
        a=cpar(jj,1)
        b=cpar(jj,2)
        d=cpar(jj,3)
        g=cpar(jj,4)
        h=cpar(jj,5)
        rj=cpar(jj,6)
330     continue
        p1=2.0d0*b*yk+d*xk+h
        p2=-2.0d0*a*xk-d*yk-g
        if(abs(p1).gt.abs(p2)) then
          p2=-2.0d0*xk*b*yk-d*xk**2-xk*h+d*yk**2+g*yk+2.0d0*a*xk*yk
          p3=-g-2.0d0*a*xk-d*yk
          r = a*p3**2/p1**2 + b - d*p3/p1
          p = -g*p3/p1 - d*p2/p1 + 2.0d0*a*p2*p3/(p1**2) + h
          q = a*p2**2/p1**2 + rj - g*p2/p1
          if (abs(r).gt.1.0d-10) then
            p=p/r
            q=q/r
            x1=-p/2.0d0
            x2=dsqrt((p/2)**2-q)
            y1=x1+x2
            y2=x1-x2
            x1=-(p2+y1*p3)/p1
            x2=-(p2+y2*p3)/p1
            p1=(xk-x1)**2+(yk-y1)**2
            p2=(xk-x2)**2+(yk-y2)**2
            if(p1.lt.p2) then
              xk=x1
              yk=y1
            else
              xk=x2
              yk=y2
            endif
          else
            yk=-q/p
            xk=-(p2+yk*p3)/p1
          endif
        else
          p1=p2
          p2=-d*xk**2 -xk*h +2.0d0*yk*a*xk -2.0d0*xk*b*yk +yk*g +d*yk**2
          p3=2.0d0*b*yk+d*xk+h
          r=b*p3**2/p1**2 + a - d*p3/p1
          p=-d*p2/p1 + 2.0d0*b*p2*p3/p1**2 - h*p3/p1 + g
          q=b*p2**2/p1**2 - h*p2/p1 + rj
          if (abs(r).gt.1.0d-10) then
            p=p/r
            q=q/r
            y1=-p/2.0d0
            y2=dsqrt((p/2)**2-q)
            x1=y1+y2
            x2=y1-y2
            y1=-(p2+x1*p3)/p1
            y2=-(p2+x2*p3)/p1
            p1=(xk-x1)**2+(yk-y1)**2
            p2=(xk-x2)**2+(yk-y2)**2
            if(p1.lt.p2) then
              xk=x1
              yk=y1
            else
              xk=x2
              yk=y2
            endif
          else
            xk=-q/p
            yk=-(p2+xk*p3)/p1
          endif
        endif
        x1=xk
        y1=yk
        call curve (x1,y1,dummy,m,jj,diff)
        if(abs(y1-yk).gt.1.0e-3) then
          goto 330
        endif
        return
      endif
      end
c......................................................................
      subroutine curved(xh,der,jj,ityp,nr)
c.... p. fellmoser
c     subroutine curve ( x,y,z,m,nr,diff )
      USE cdata
      USE iofile
      USE tdata
      implicit double precision (a-h,o-z)
C.... Auswertung von  definierten Funktionen
C     ityp 1  to 50  : standart curves     ( same in all feap version )
C     ityp 51 to 100 : user defined curves ( only in user version )
C     input :   x    = argument
C               y    = argument
C               z    = argument
C               nr   = curve number
C     output:   tol  = value  (only for 3-D)
C               y    = value  at x   (2-d)
C               z    = value  at x,y (3-d)
C               diff = 0 if x,y,z lies on curve (else diff<>0)
C     variables cpar = curve parametres
C.....................................
      common /curvedat/cpar(20,8),nrt1(20,3),ic,nbe,nn3
      dimension xx(3),yy(3),xy(3)
czr  +         ,m(*)
     +         ,der(3),derh(3),xh(3)
      x=0
czr   data pi /3.14159265359/
c
      do 100 i=1,20                   ! max 20 curves
        if( nrt1(i,1) .eq. jj ) then
          ityp = nrt1(i,2)             ! curve typ
          nn3  = nrt1(i,3)
         goto 200
        else
        endif
100   continue
                   write(iow,2000) jj     ! error curve not input
      if(ior.lt.0) write(*,  2000) jj
      stop
200   continue
c
C***********************************************************
c                  2 - D  Kurven
C***********************************************************
C
      if ( ityp .eq. 1 ) then
c
c.....Gerade durch zwei Punkte (x1,y1) (x2,y2)
cpf>
        x1 = cpar(jj,1)
        y1 = cpar(jj,2)
        x2 = cpar(jj,3)
        y2 = cpar(jj,4)
c
        derh(1)=cpar(jj,1)-cpar(jj,3)
        derh(2)=cpar(jj,2)-cpar(jj,4)
        derh(3)=0
c
        der(1)=derh(1)/sqrt(derh(1)**2+derh(2)**2)
        der(2)=derh(2)/sqrt(derh(1)**2+derh(2)**2)
        der(3)=0.0d0
cpf<
      elseif ( ityp .eq. 2 ) then
c
c.... Gerade durch einen Punkt (x1,y1) und Steigung b
c
        x1 = cpar(jj,1)
        y1 = cpar(jj,2)
        b  = cpar(jj,3)
        a  = y1 - b*x1
        y  = a + b*x
        z  = b
c
        derh(1)=cpar(jj,1)-xh(1)
        derh(2)=cpar(jj,2)-xh(2)
        derh(3)=0
c
        der(1)=derh(1)/sqrt(derh(1)**2+derh(2)**2)
        der(2)=derh(2)/sqrt(derh(1)**2+derh(2)**2)
        der(3)=0.0d0

c
      elseif ( ityp .eq. 3 ) then
c
C.... Halbkreis mit dem Radius r, Steuerparameter nx
c     und Mittelpunkt (x1,y1)
c
        r  = cpar(jj,1)
        xv = cpar(jj,2)
czr     if(abs(xv).ne.1.d0) xv=1.d0
        if( .not.(abs(xv).lt.1.d0.or.abs(xv).gt.1.d0)) xv=1.d0
        x0 = cpar(jj,3)
        y0 = cpar(jj,4)
cpf     if(x.gt.(x0+r) .or. x.lt.(x0-r)) then
cpf       y= 99999999.d0
cpf     else
cpf       y = y0 + xv* dsqrt(r*r-(x-x0)*(x-x0))
cpf     endif
cpf>
        derh(1)=-2.0d0*(xh(2)-cpar(jj,4))
        derh(2)=2.0d0*(xh(1)-cpar(jj,3))
c
        der(1)=derh(1)/sqrt(derh(1)**2+derh(2)**2)
        der(2)=derh(2)/sqrt(derh(1)**2+derh(2)**2)
        der(3)=0.0d0
cpf<
c
C****************************************************
c                   3 - D  Kurven
C****************************************************
C
      elseif ( ityp .eq. 32 ) then
c
C.....Kugel mit dem Radius r und Mittelpunkt (x0,y0,z0)
c
        r  = cpar(jj,1)
        x0 = cpar(jj,2)
        y0 = cpar(jj,3)
        z0 = cpar(jj,4)
cpf>
        derh(1)=-2.0d0*(xh(2)-cpar(jj,3))
        derh(2)=2.0d0*(xh(1)-cpar(jj,2))
        derh(3)=2.0d0*(xh(3)-cpar(jj,4))
c
        der(1)=derh(1)/sqrt(derh(1)**2+derh(2)**2)
        der(2)=derh(2)/sqrt(derh(1)**2+derh(2)**2)
        der(3)=0.0d0
cpf<
      elseif ( ityp .eq. 31 ) then
c
C.... Ebene mit der Normalen x0 y0 z0 und der Konstanten c0
c         (x0*x+y0*y+z0*z+c0=0)
cpf>
        der(1)=cpar(jj,2)
        der(2)=cpar(jj,1)
        der(3)=cpar(jj,3)
cpf<
      elseif ( ityp .eq. 30 ) then
c
C.....Gerade mit den 2 Punkten x1,y1,z1 und x2,y2,z2,
c
        x1 = cpar(jj,1)
        y1 = cpar(jj,2)
        z1 = cpar(jj,3)
        x2 = cpar(jj,4)
        y2 = cpar(jj,5)
        z2 = cpar(jj,6)
        xx(1)=x2-x1
        xx(2)=y2-y1
        xx(3)=z2-z1
c       call vnorm (xx,xl)
        yy(1)=x-x1
        yy(2)=y-y1
        yy(3)=z-z1
        call vcross(xx,yy,xy)
        if (dabs(xy(1))+dabs(xy(2))+dabs(xy(3)).lt.1e-5) then
          diff=0.d0
        else
          diff=10.d0
        endif
c
      elseif ( ityp .eq. 50 ) then
c
c.... curve defined via points
czr   call ctyp50 (x,y,m(nbe),z,nn3)
C
c.... user defined curves
      elseif ( ityp.gt.50 .and. ityp.le.100 ) then
czr     call cvuser ( x,y,z,ityp,m,nr,diff )
      else
                     write(iow,2001) ityp
        if(ior.lt.0) write(*,2001) ityp
        stop
      endif
      return
2000  format('***ERROR*** curve number ',i3,' in not been input')
2001  format('***ERROR*** curve type ',i3,' in not defined')
      end
c
c.....................................................................
      subroutine rest(x,iek,f,d,id,ndm,ndf,nen1)
c......................................................................
c.... copied out of fadaptz.for   (rene ziegler)
c......................................................................
      USE bdata
      USE cdata
      USE eldata
      implicit double precision (a-h,o-z)
      dimension x(ndm,*),iek(nen1,*),f(ndf,*),d(*),id(ndf,*)
      dimension idl(16)
      write(3,1000) (head(i),i=1,20)
      write(3,2000) numnp,numel,nummat,ndm,ndf,nen
c..
      write(3,'(a)')'    '
      write(3,'(a)')'coor'
      do 300 i=1,numnp
***     xll=dsqrt(x(1,i)**2+x(2,i)**2+x(3,i)**2)
***     write(*,*) 'laenge=',xll

300   write(3,3000)i,(x(ii,i),ii=1,ndm)
c..
      write(3,'(a)')'    '
      write(3,'(a)')'elem'
      do i=1,numel
        write(3,4000)i,iek(nen1,i),(iek(ii,i),ii=1,4)
      enddo
c..
      write(3,'(a)')'    '
      write(3,'(a)')'forc'
      do 500 i=1,numnp
        nsteuer=0
        do 510 ii=1,ndf
czr       if(f(ii,i).eq.0.d0)goto 510
          if(abs(f(ii,i)).lt.1e-8)goto 510
          nsteuer=1
510     continue
czr520  if(nsteuer.eq.1) write(3,5000)i,0,(f(ii,i),ii=1,ndf)
        if(nsteuer.eq.1) write(3,5000)i,0,(f(ii,i),ii=1,ndf)
500   continue
c....
      write(3,'(a)')'    '
      write(3,'(a)')'boun'
      do 600 i=1,numnp
        nsteuer=0
        do 610 ii=1,ndf
          if(id(ii,i).ge.0)goto 610
          nsteuer=1
610     continue
      if(nsteuer.eq.1)then
      do 620 ii=1,ndf
      idl(ii)=0
620   if(id(ii,i).lt.0) idl(ii)=1
      if(i.lt.10)then
      write(3,6001)i,0,(idl(ii),ii=1,ndf)
      elseif(i.lt.100)then
      write(3,6002)i,0,(idl(ii),ii=1,ndf)
      elseif(i.lt.1000)then
      write(3,6003)i,0,(idl(ii),ii=1,ndf)
      elseif(i.lt.10000)then
      write(3,6004)i,0,(idl(ii),ii=1,ndf)
      else
      write(3,6005)i,0,(idl(ii),ii=1,ndf)
      endif
      endif
600   continue
      if(iel.eq.11)then
CC    write(3,'(a)')'    '
CC    write(3,'(a)')'mate'
CC    if(iel.eq.11)then
CC    i1=idint(d(5))
CC    i2=idint(d(6))
**    write(3,7000)ma,iel
CC    write(3,7001)ma,1
CC    write(3,7100)d(16),d(17),d(4),i1,i2,0
      elseif(iel.gt.70.and.iel.lt.80)then
CC    write(3,'(a)')'    '
CC    write(3,'(a)')'mate'
CC    write(3,7000)ma,iel
CC    write(3,7101)d(16),d(17),0.0,0,0,0
      else
      write(3,*) '***error-no mate for adadaptivity'
      endif
c..
      write(3,'(a)')'    '
      write(3,'(a)')'end '
      write(3,'(a)')'    '
      write(3,'(a)')'macr'
      write(3,'(a)')'plon'
      write(3,'(a)')'plot,boun'
      write(3,'(a)')'plot,load'
      write(3,'(a)')'plot,mesh'
      write(3,'(a)')'end '
      write(3,'(a)')'    '

      return
1000  format(a4,1x,19(a4))
2000  format(6(i10,','))
3000  format(i10,',0',6(',',f12.4))
4000  format(5(i5,','),i5,',0')
5000  format(2(i10,','),6(f10.4,','))
6001  format(i1,',',15(i1,','))
6002  format(i2,',',15(i1,','))
6003  format(i3,',',15(i1,','))
6004  format(i4,',',15(i1,','))
6005  format(i5,',',15(i1,','))
czr7000  format(16(i5,','))
czr7001  format(16(i1,','))
czr7100  format(f10.0,',',2(f10.9,','),3(i10,','))
czr7101  format(f10.0,',',2(f10.9,','),3(i10,','))
      end
c
cpf   subroutine msmosh(ike,ikez,x,ix,lct,ivn,itnn,ndm,ndf,iael0,igl,m)
      subroutine msmosh(ike,ikez,x,ix,lct,ivn,itnn,ndm,ndf,iael0,igl,m,
     +                  ict31,ict32,ct,l)
c
c.... This macro command optimizes geometrically an given shell mesh
C
C     call: msmo,i1i2,c1,c2
C
C     i1 : weight factor for equality of sidelengths
C          (first two characters of charcter-chain) = integer i1
C     i2 : weight factor for 90-degree-angles
C          (third and fourth character of charcter-shain) = integer i2
C          (for instance i1i2 = 0210 means : i1=2 and i2=10
C     c1 = i3*1000000+j1*1000+j2
C          i3 : weight factor for out-of-plane-distorsion of element
C          j1 : control parameter : (loop,msmo,.. as adaptive loop necessary)
C             j1=1 : execute mesh smoothing only at the last adaptive step
C             j1=2 : execute mesh smoothing at first and last adaptive step
C             j1=3 : execute mesh smoothing at all adaptive refinement steps
C          j2 : number of vertices of the fe-model ,
C               (this vertices must be the first ones at the nodal-point list)
C     c2 : k1*1000000 + k2*1000 + k3
C          k1 : number of newton steps in every local iteration
C          k2 : function choosing parameter (1-7 available)
C
C                      / a-b \ 2              /   a*b     Pi \ 2
C             f1 = w1*|  ---  |  + w2*arccos |  ------- - --  |
C                      \ a+b /                \ |a|*|b|   2  /
C
C             f2 = simple averaging of single element optima
C
C          k3 : number of iteration steps over all elements
C               (useful 1,2 or 3)
C
C     N.b.: All boundaries must be defined by curves
C
      USE bdata
      USE cdata
      USE comfil
      USE iofile
      USE mdata
      USE mdat2
      implicit double precision (a-h,o-z)
czr   common /ycur1/ cpar(100,8),nrt(100,3),ic
      common /curvedat/  cpar(20,8),nrt(20,3),ic,nbe,nn3
      common /msm1/ lop,lwe,lgl,lsw,lsi
      common /msm2/ iln,ila,ile,mel,ngl
      common /fvd/ w1,w2,w3,ieh
      common /adap2 / lad
      common /mdat3/  n9a,n9b,n9c,n9d,n10a,n11c,n11d,n14a,numelo,numnpo
      common /sdata/ ndfx,ndmx,nen1,nstx
      character*4 lct(*)
cpf   character*2 llct
      logical lop,lbn,lou,lf,lad,lwe,lne,lb,lsg,log,lopn,lna,les
     +      ,lgl,lsw,lsi,lshel
     +      ,pcomp,lwrite
      dimension ike(*),ikez(*),x(ndm,*),ix(nen+4,*),xm(3),m(*),
     +      grad(3),ta(3),tb(3),iael0(*),hv1(3),hv2(3),hv3(3),hv4(3),
     +      rhv(9,16),igl(*),
     +      ct(3,*)
cTEMP
c.the      data reps /1.0d-3/
      data reps /5.0d-2/
c     reps =10* tol
cpf   logical for write-output
      lwrite =.true.
cryz  logical fr Abarbeitung der Schalenverschneidungsroutinen
      lsi =.true.
cryz end
cryz  logical fr Aufruf der subroutine quality
cpf   lsw =.true.
      lsw =.false.
cryz end
cryz
      if (lsi) then
      do 109 iz=1,numnp
           do 107 it=1,numnp
               deltax = x(1,iz) - x(1,it)
               deltay = x(2,iz) - x(2,it)
               deltaz = x(3,iz) - x(3,it)
c.the               if ((abs(deltax).lt.1.0E-4).and.(abs(deltay).lt.1.0E-4)
c.the     +              .and.(abs(deltaz).lt.1.0E-4).and.iz.ne.it) then
               if ((abs(deltax).lt.1.0E-5).and.(abs(deltay).lt.1.0E-5)
     +              .and.(abs(deltaz).lt.1.0E-5).and.iz.ne.it) then
                    igl(iz) = it
                    goto 109
               endif
 107     continue
            igl(iz) = 0
 109     continue
      endif
cryz end
      les=.false.
c.........log = only gradient iteration
      log=.false.
c.........lsg = switch to gradient iteration
      lsg=.false.
c ..............w3 ... third weigth factor
cpf   w3=ivn/1000000
cpf   ivnn=ivn-w3*1000000
cpf   isq=ivnn/1000
cpf>
      isq=1
cpf<
cpf   iv=ivnn-isq*1000
cpf>
      iv=ct(1,l)-1
      if (abs(ct(1,l)).lt.0.0001)then
       iv=0
      endif
cpf<
      if(((isq.eq.1).and.(ila.eq.ile)).or.
     +   ((isq.eq.2).and.((ila.eq.1).or.(ila.eq.ile))).or.
cric     +    (isq.eq.3)) then
     +    (isq.ge.3)) then
        lne = .not.((ila.eq.ile).or.(ila.eq.1))
        lna=.false.
cpf     if (fsav(1:3).ne.'nul') then
cpf       lou=.true.
cpf       inquire (file=fsav,opened = lopn)
cpf       if (.not.lopn) then
cpf         if (lad) then
cpf           open (3,file=fsav,status='unknown')
cpf         else
cpf           lna=.true.
cpf           tsav=fsav
cpf           fsav(1:1)='a'
cpf           open (3,file=fsav,status='unknown')
c              rewind(3)
cpf   call rest(m(n11b),m(n8),m(n9),m(n10),m(n6),m(n7),ndm,ndf,nen+4)
cpf         endif
cpf       endif
cpf     else
          lou=.false.
          lopn=.true.
cpf     endif

        rkap=0.7d0
        delta=0.1d0
        alp=0.5d0
        if (ndm.ne.3) then
          write(*,2000)
          write(iow,2000)
        endif
        lf=.true.
c .................w1 ... first weigth factor
cpf     llct=lct(1:2)
cpf     read(llct,4000,err=110) i
cpf     if (i.ne.0) lf=.false.
cpf     w1=dble(i)
c.... w2 second weight factor
cpf     llct=lct(3:4)
cpf     read(llct,4000,err=110) i
cpf     if (i.ne.0) lf=.false.
cpf     w2=dble(i)
cpf     mnv=itnn/1000000
cpf>
        mnv=1
cpf<
        if (mnv.eq.0) then
          write(*,9100)
          write(iow,9100)
          stop
        endif
cpf     itnnn=itnn-mnv*1000000
cpf     iwf=itnnn/1000
cpf>
        iwf=2
cpf<
        if (iwf.eq.0) then
          write(*,9000)
          write(iow,9000)
          stop
        endif
cpf>
        if(pcomp(lct(l),'wesp',4))then
      write(*,*)'3-d mesh optimization only for inner nodes available'
          return
        endif
        if(pcomp(lct(l),'oesp',4))then
          continue
        else
          write(*,*)'Enter new mesh optimization subroutine'
          return
        endif
cpf<
c.... number of iterations
cpf   itn=itnnn-iwf*1000
cpf>
      itn=ct(2,l)
      if (itn.eq.0) then
        itn=1
      endif
cpf<
      it=w1
cpf   write(iow,8000) iwf,it,i,itn,mnv
cpf   write(*,8000) iwf,it,i,itn,mnv
      write(iow,8000) ict31,ict32,itn,mnv
      write(*,8000) ict31,ict32,itn,mnv
c.... loop-number of iterations
      do 100 it=1,itn
        lb = (mod(it,2).eq.0)
cl        do 165 i=iv+1,numnp
          do 165 iii=iv+1,numnp
c.the>
c.Keine Glaettung von Punkten, die auf der Kante zweier Uebergangselemente
c liegen.
            nel = ikez(iii+1)-ikez(iii)
            if (lad.and.(nel.ne.4)) then
c
               do 166 ij=1,nel
                  it1=ikez(iii)
                  it2=ike(it1+ij-1)
                  it3=iael0 (it2)
166               if (it3.lt.0) goto 165
            endif
ctmp>
            if (iii.eq.229) then
              write(*,*)
            endif
c.the<
            if (lb) then
              i=numnp-iii+iv+1
            else
cric       i: number of node to optimize
              i=iii
            endif
cl end
            xm(1)=x(1,i)
            xm(2)=x(2,i)
            xm(3)=x(3,i)
            ieh=ikez(i+1)-ikez(i)
            if (ieh.gt.8) then
              write (iow,6000) i
              write (*,6000) i
              stop
            endif
c.... loop - neh: number of elements connected to node
            do 330 j=1,ieh
              ihv=ike(ikez(i)+j-1)
c.... ihv = number of element
              if (lad) then
                if (isq.ne.4) then
                  if (lne.and.(iael0(ihv).lt.0)) then
                    goto 165
                  endif
                endif
              endif
cric          if (lne.and.(iael0(ihv).lt.0)) goto 165
              do 660 k=1,nen
                if(ix(k,ihv).eq.i) then
                  ilhvh=mod(k+2,4)
                  iuhvh=mod(k,4)
                  iohvh=mod(k+1,4)
                  ilhv=ilhvh+1
                  iuhv=iuhvh+1
                  iohv=iohvh+1
                  il=ix(ilhv,ihv)
                  iu=ix(iuhv,ihv)
                  io=ix(iohv,ihv)
                endif
660           continue
cric          write coordinates of other nodes
              rhv(1,j) = x(1,il)
              rhv(2,j) = x(2,il)
              rhv(3,j) = x(3,il)
              rhv(4,j) = x(1,iu)
              rhv(5,j) = x(2,iu)
              rhv(6,j) = x(3,iu)
              rhv(7,j) = x(1,io)
              rhv(8,j) = x(2,io)
              rhv(9,j) = x(3,io)
330         continue
cryz
         if (igl(i).ne.0) then
            ii = igl(i)
            ieh2=ikez(ii+1)-ikez(ii)
            if (ieh2.gt.8) then
              write (iow,6000) ii
              write (*,6000) ii
              stop
            endif
            iehpu1 = ieh + 1
            iehpu2 = ieh + ieh2
            do 331 j=iehpu1,iehpu2
              jneu = j - iehpu1 + 1
              ihv2=ike(ikez(ii)+jneu-1)
c.... ihv2 = number of element
              if (lad) then
                if (isq.ne.4) then
                  if (lne.and.(iael0(ihv2).lt.0)) then
                    goto 165
                  endif
                endif
              endif
cric          if (lne.and.(iael0(ihv2).lt.0)) goto 165
              do 661 k=1,nen
                if(ix(k,ihv2).eq.ii) then
                  ilhvh=mod(k+2,4)
                  iuhvh=mod(k,4)
                  iohvh=mod(k+1,4)
                  ilhv=ilhvh+1
                  iuhv=iuhvh+1
                  iohv=iohvh+1
                  il2=ix(ilhv,ihv2)
                  iu2=ix(iuhv,ihv2)
                  io2=ix(iohv,ihv2)
                endif
661           continue
cric          write coordinates of other nodes
              rhv(1,j) = x(1,il2)
              rhv(2,j) = x(2,il2)
              rhv(3,j) = x(3,il2)
              rhv(4,j) = x(1,iu2)
              rhv(5,j) = x(2,iu2)
              rhv(6,j) = x(3,iu2)
              rhv(7,j) = x(1,io2)
              rhv(8,j) = x(2,io2)
              rhv(9,j) = x(3,io2)
331         continue
         ieh = iehpu2
         endif
cryz end
            xk=xm(1)
            yk=xm(2)
            zk=xm(3)
            icnr=0
            do 801 jj=1,ic-1
cpf           if (nrt(jj,3).eq.15) then
              if ((nrt(jj,3).eq.15).or.(nrt(jj,3).eq.17).or.
     +            (nrt(jj,3).eq.18)) then
                call curve (xm(1),xm(2),xm(3),m,jj,diff)
                if(dabs(diff).lt.reps) then
                  icnr=jj
                  goto 499
                endif
              endif
801         continue
            write(*,*) 'stop:subroutine msmosh'
            if (icnr.eq.0) stop
499         continue
c
c.... distinguish boundary node and inner node
c
            ircn=0
            lbn=.false.
            do 825 jj=1,ic-1
cpf           if ((nrt(jj,3).eq.14).or.
              if (((nrt(jj,3).eq.17).and.(jj.ne.icnr)).or.
c.the     +            ((nrt(jj,3).ge.16).and.(nrt(jj,3).le.25)).or.
     +            ((nrt(jj,3).eq.15).and.(jj.ne.icnr)).or.
     +            ((nrt(jj,3).eq.18).and.(jj.ne.icnr))) then
                call curve (xm(1),xm(2),xm(3),m,jj,diff)
c.the                if(dabs(diff).lt.reps) then
                if(dabs(diff).lt.1e-3) then
                  lbn=.true.
                  ircn=jj
                  goto 495
                endif
              endif
825         continue
495         continue
C
            if ((iwf.eq.1).or.((iwf.eq.3).and.(lbn)).or.(iwf.eq.4)) then
C
c.... compute enhanced start value
              if(iwf.eq.4)call simpav(xk,yk,zk,icnr,grad,ieh,rhv,ircn,m)
c.....Gauss-Newton method
              do 200 nv=1,mnv
                if (lbn) then
c.... case of boundary or shell intersection node
                  call curveg ( xk,yk,zk,icnr,grad )
                  call curveg ( xk,yk,zk,ircn,ta )
                  call vcross  ( grad,ta,tb )
                  f1=0.0d0
                  f2=0.0d0
                  zfo=0.0d0
                  do 567 ih=1,ieh
                    x1=rhv(1,ih)
                    y1=rhv(2,ih)
                    z1=rhv(3,ih)
                    x2=rhv(4,ih)
                    y2=rhv(5,ih)
                    z2=rhv(6,ih)
                    x3=rhv(7,ih)
                    y3=rhv(8,ih)
                    z3=rhv(9,ih)
             call sf(xk,yk,zk,x1,y1,z1,x2,y2,z2,x3,y3,z3,grad,iwf,2,hv2)
             call sf(xk,yk,zk,x1,y1,z1,x2,y2,z2,x3,y3,z3,grad,iwf,3,hv3)
             call sf(xk,yk,zk,x1,y1,z1,x2,y2,z2,x3,y3,z3,grad,iwf,4,hv4)
             call sf(xk,yk,zk,x1,y1,z1,x2,y2,z2,x3,y3,z3,grad,iwf,1,hv1)
                    do 234 jj=1,3
             f1=f1+hv1(jj)*(hv2(jj)*tb(1)+hv3(jj)*tb(2)+hv4(jj)*tb(3))
             f2=f2+(hv2(jj)*tb(1)+hv3(jj)*tb(2)+hv4(jj)*tb(3))**2
                      zfo=zfo+hv1(jj)**2
234                 continue
567               continue
                  a11=f1/f2
                  dxk=a11*tb(1)
                  dyk=a11*tb(2)
                  dzk=a11*tb(3)
                else
c.... case of inner node
c.... 2 orthogonal tangent vectors t1 and t2 of shell-surface
                  call curveg ( xk,yk,zk,icnr,grad )
                  tb(1)=xk-rhv(1,1)
                  tb(2)=yk-rhv(2,1)
                  tb(3)=zk-rhv(3,1)
                  call vcross(tb,grad,ta)
                  call vcross(ta,grad,tb)
                  call vnorm (grad,dummy)
c.... assemble tangential system
                  a11=0.0d0
                  a12=0.0d0
                  a22=0.0d0
                  f1=0.0d0
                  f2=0.0d0
                  zfo=0.0d0
                  do 122 ih=1,ieh
                    x1=rhv(1,ih)
                    y1=rhv(2,ih)
                    z1=rhv(3,ih)
                    x2=rhv(4,ih)
                    y2=rhv(5,ih)
                    z2=rhv(6,ih)
                    x3=rhv(7,ih)
                    y3=rhv(8,ih)
                    z3=rhv(9,ih)
             call sf(xk,yk,zk,x1,y1,z1,x2,y2,z2,x3,y3,z3,grad,iwf,2,hv1)
                    f1b1=hv1(1)*ta(1)
                    f2b1=hv1(2)*ta(1)
                    f3b1=hv1(3)*ta(1)
                    f1b2=hv1(1)*tb(1)
                    f2b2=hv1(2)*tb(1)
                    f3b2=hv1(3)*tb(1)
             call sf(xk,yk,zk,x1,y1,z1,x2,y2,z2,x3,y3,z3,grad,iwf,3,hv1)
                    f1b1=f1b1+hv1(1)*ta(2)
                    f2b1=f2b1+hv1(2)*ta(2)
                    f3b1=f2b1+hv1(3)*ta(2)
                    f1b2=f1b2+hv1(1)*tb(2)
                    f2b2=f2b2+hv1(2)*tb(2)
                    f3b2=f3b2+hv1(3)*tb(2)
             call sf(xk,yk,zk,x1,y1,z1,x2,y2,z2,x3,y3,z3,grad,iwf,4,hv1)
                    f1b1=f1b1+hv1(1)*ta(3)
                    f2b1=f2b1+hv1(2)*ta(3)
                    f3b1=f2b1+hv1(3)*ta(3)
                    f1b2=f1b2+hv1(1)*tb(3)
                    f2b2=f2b2+hv1(2)*tb(3)
                    f3b2=f3b2+hv1(3)*tb(3)
             call sf(xk,yk,zk,x1,y1,z1,x2,y2,z2,x3,y3,z3,grad,iwf,1,hv1)
                    zfo=zfo+hv1(1)**2+hv1(2)**2+hv1(3)**3
                    a11=a11+f1b1**2+f2b1**2+f3b1**2
                    a12=a12+f1b1*f1b2+f2b1*f2b2+f3b1*f3b2
                    a22=a22+f1b2**2+f2b2**2+f3b2**2
                    f1=f1+f1b1*hv1(1)+f2b1*hv1(2)+f3b1*hv1(3)
                    f2=f2+f1b2*hv1(1)+f2b2*hv1(2)+f3b2*hv1(3)
122               continue
c.... solve tangential system to rh1 and rh2
                  det = a11*a22-a12**2
                  rdet = abs(a11*a22)+a12**2
                  if (abs(det).lt.10d-6*rdet) then
                    write(*,3000) i
                    goto 165
                  endif
                  rh1 = (f1*a22 - f2*a12)/det
                  rh2 = (f2*a11 - f1*a12)/det
                  dxk = rh1*ta(1) +rh2*tb(1)
                  dyk = rh1*ta(2) +rh2*tb(2)
                  dzk = rh1*ta(3) +rh2*tb(3)
                endif
                xkm1= xk
                ykm1= yk
                zkm1= zk
                gam = 1.0d0
                izv=0
144             continue
                xk  = xk - gam*dxk
                yk  = yk - gam*dyk
                zk  = zk - gam*dzk
                call projec(xk,yk,zk,icnr,ircn,m)
                zf  = 0.0d0
                do 133 ih=1,ieh
                  x1=rhv(1,ih)
                  y1=rhv(2,ih)
                  z1=rhv(3,ih)
                  x2=rhv(4,ih)
                  y2=rhv(5,ih)
                  z2=rhv(6,ih)
                  x3=rhv(7,ih)
                  y3=rhv(8,ih)
                  z3=rhv(9,ih)
             call sf(xk,yk,zk,x1,y1,z1,x2,y2,z2,x3,y3,z3,grad,iwf,1,hv1)
                  zf=zf+hv1(1)**2+hv1(2)**2+hv1(3)**3
133             continue
                if (zfo.lt.zf) then
                  xk=xkm1
                  yk=ykm1
                  zk=zkm1
                  if (izv.le.5) then
                    izv=izv+1
                    gam=gam/2.0d0
                    goto 144
                  endif
c.... switch to gradient iteration
                  write(*,5000) i,nv
                  write(iow,5000) i,nv
                  goto 165
                endif
200           continue
            elseif (((iwf.eq.2).or.(iwf.eq.3)).and.(.not.lbn)) then
c.... simple averaging
          call simpav(xk,yk,zk,icnr,grad,ieh,rhv,ircn,m)
c              call curveg ( xk,yk,zk,icnr,grad )
c              call vnorm ( grad,dummy )
cc... projection to tangential plane
c              do 345 j=1,ieh
c                do 456 ii=1,2
c                  ij = 3*(ii-1)
c                  rshv = (rhv(1+ij,j)-xk)*grad(1) + (rhv(2+ij,
c     +                  j)-yk)*grad(2) + (rhv(3+ij,j)-zk)*grad(3)
c                  rhv(1+ij,j) = rhv(1+ij,j) - rshv*grad(1)
c                  rhv(2+ij,j) = rhv(2+ij,j) - rshv*grad(2)
c                  rhv(3+ij,j) = rhv(3+ij,j) - rshv*grad(3)
c456             continue
c345           continue
cC... 2 orthogonal tangent vectors
c              ta(1)=rhv(1,1)-xk
c              ta(2)=rhv(2,1)-yk
c              ta(3)=rhv(3,1)-zk
c              call vnorm (ta,dummy)
c              call vcross(ta,grad,tb)
c              xkm1=0.0d0
c              ykm1=0.0d0
c              do 300 j=1,ieh
cc... transformation to 2-dimensional problem
c          x1=(rhv(1,j)-xk)*ta(1)+(rhv(2,j)-yk)*ta(2)+(rhv(3,j)-zk)*ta(3)
c          y1=(rhv(1,j)-xk)*tb(1)+(rhv(2,j)-yk)*tb(2)+(rhv(3,j)-zk)*tb(3)
c          x2=(rhv(4,j)-xk)*ta(1)+(rhv(5,j)-yk)*ta(2)+(rhv(6,j)-zk)*ta(3)
c          y2=(rhv(4,j)-xk)*tb(1)+(rhv(5,j)-yk)*tb(2)+(rhv(6,j)-zk)*tb(3)
cC... 2-dimensional averaging
c                rh1=(x1+x2-y1+y2)/2.0d0
c                rh2=(y1+y2+x1-x2)/2.0d0
c                f1=(x1+x2+y1-y2)/2.0d0
c                f2=(y1+y2-x1+x2)/2.0d0
c                if(dsqrt(rh1**2+rh2**2).lt.dsqrt(f1**2+f2**2)) then
c                  xkm1=xkm1+rh1
c                  ykm1=ykm1+rh2
c                else
c                  xkm1=xkm1+f1
c                  ykm1=ykm1+f2
c                endif
c300           continue
c              rh1 = xkm1/ieh
c              rh2 = ykm1/ieh
c              xk = xk + rh1*ta(1) + rh2*tb(1)
c              yk = yk + rh1*ta(1) + rh2*tb(1)
c              zk = zk + rh1*ta(1) + rh2*tb(1)
c              call projec(xk,yk,zk,icnr,ircn,m)
              x(1,i)=xk
              x(2,i)=yk
              x(3,i)=zk
            endif
c            if (lbn) then
c              do 370 ij=1,ieh
c                x1=rhv(1,ij)
c                y1=rhv(2,ij)
c                x2=rhv(3,ij)
c                y2=rhv(4,ij)
c                if ((xk-x2)*(yk-y1)-(xk-x1)*(yk-y2).lt.0.0d0) then
c                  write(*,3700) i
c                  write(iow,3700) i
c                  goto 165
c                endif
c370           continue
c            endif
            x(1,i)=xk
            x(2,i)=yk
            x(3,i)=zk
cryz
         if (igl(i).ne.0) then
            ii = igl(i)
            x(1,ii)=xk
            x(2,ii)=yk
            x(3,ii)=zk
         endif
cryz end
c       if (lou) call rest(
c    #    m(n11b),m(n8),m(n9),m(n10),m(n6),m(n7),ndm,ndf,nen+4)
165       continue
c.the>
c..Neuberechnung der Belastungen, da Knoten sich verschoben haben;
c  wichtig bei Def. z.B. einer Streckenlast.
      if (lwrite) then
       write(*,*)'sr msmosh - check if linear loads are correct!'
       lwrite=.false.
      endif
      !call pzero (m(n10),numnp*ndf)
      gloa = 0.d0
      do 35 i=0,numnpo*ndf*ipr-1
35    m(n10+i)=m(n10a+i)
cpf   call pld (ndf,ndm,numel,nen,nen1,m(n8),m(n9),m(n10),m, numnp)
cpf   write(*,*)'sr pld not implemented SR MSMOSH!!!'
cpf   call sld (ndf,ndm,numel,nen,nen1,m(n8),m(n9),m(n10),m)
cpf   write(*,*)'sr sld not implemented SR MSMOSH!!!'
c.the<
        if (lou) call rest(
     #    bang,coor,econ,gloa,edma,psid,ndm,ndf,nen+4)
100     continue
        if (.not.lopn) close(3)
        if (lna) then
          close(3)
cpf       fsav=tsav
          fsav='test?'
          write(*,*)'subroutine msmosh - check fsav!'
        endif
      endif
      write(*,8800)
cryz
        if (lsw) then
          lshel =.true.
          call quality (ike,ikez,x,ix,xm,ndm,lad,lne,lshel,isq,m,iael0)
        endif
cryz end(testq)
      return
cpf110   call perror('msmo')
cpf   stop
cpf1009  format(a/)
2000  format('mesh optimization subroutine only for 2 or 3-dimensional',
     +      ' meshes available yet')
3000  format('0 valued det. in newton-method (msmooth,',
     +        'node ',i5,' ',i3,'th iteration)')
c3100  format('node',i5,' and further 2 nodes of 1 element on curve',i3)
c3200  format('no complete set of boundary points found around node',i5,
c     +        ' on curve',i3)
c3300  format('2 same boundary points at curve',i3,' ; node=',i5)
c3700  format('neg. jacobi det. in element at node',i3,
c     +      ' ; reset to start value')
cpf4000  format(i2)
5000  format('switched to gradient iteration at node',i5,
     +      ' ; newton step:',i4)
6000  format('error in macro command msmo : more then 8 elements at',
     +       ' node ',i5)
8000  format('macro msmo: function: *  w1= ',i2,'  w2= ',i2,
     +       '  iterations: ',i3,' max. newton steps: ',i3)
8800  format('end of macro msmo')
9000  format('new function-choosing-parameter: n2=fn*1000+itn',/,
     +       '(fn=number of function for minimum problem; itn=number ',
     +       'of global iterations)')
9100  format('new function-choosing-parameter: n2=mn*1000000+fn*1000',
     +        '+itn',/5x,'mn=maximal number of newton steps in local ',
     +        'iterations',/5x,'fn=number of function for minimum',
     +        ' problem ',/5x,'itn=number of global iterations',/)
      end
***********************************************************************
      subroutine sf(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,e,iwf,isp,sf1)
      USE iofile
      USE mdat2
      implicit double precision (a-h,o-z)
      common /fvd/ w1,w2,w3,ieh
      dimension e(3),sf1(3)
c
c.... function value f (isp=1), first derivatives df/dx (isp=2)
C     and df/dy (isp=3) and  df/dy (isp=4) of functions fi(x,y,z,w1,w2,w3):
C
C                       a**2-b**2
C     f1 = sf1(1) = w1* ---------
C                       a**2+b**2
C
C                                a*b     Pi
C     f2 = sf1(2) = w2*arccos  ------- - --
C                              |a|*|b|   2
C
C     f3 = sf1(3) = w3*a*e
C
      Pi=3.1415926d0
      go to(100,200,300,400), isp
100     continue
      t2 = x-x1
      t3 = t2**2
      t5 = y-y1
      t6 = t5**2
      t8 = z-z1
      t9 = t8**2
      t11 = sqrt(t3+t6+t9)
      t13 = x-x2
      t14 = t13**2
      t16 = y-y2
      t17 = t16**2
      t19 = z-z2
      t20 = t19**2
      t22 = sqrt(t14+t17+t20)
      sf1(1) = w1*(t11-t22)/(t11+t22)
      sf1(2) = w2*(acos((t2*t13+t5*t16+t8*t19)/t11/t22)-Pi/2)
      sf1(3) = w3*(t2*e(1)+t5*e(2)+t8*e(3))
      return
200     continue
c       c(1)=x1-x3
c       c(2)=y1-y3
c       c(3)=z1-z3
c       d(1)=x2-x3
c       d(2)=y2-y3
c       d(3)=z2-z3
c       call vcross(c,d,e)
c       call vnorm (e,dummy)
      t1 = -x1
      t2 = x+t1
      t3 = t2**2
      t5 = y-y1
      t6 = t5**2
      t8 = z-z1
      t9 = t8**2
      t10 = t3+t6+t9
      t11 = sqrt(t10)
      t12 = 1/t11
      t13 = 2*x
      t15 = t13-2*x1
      t17 = t12*t15/2
      t18 = -x2
      t19 = x+t18
      t20 = t19**2
      t22 = y-y2
      t23 = t22**2
      t25 = z-z2
      t26 = t25**2
      t27 = t20+t23+t26
      t28 = sqrt(t27)
      t29 = 1/t28
      t31 = t13-2*x2
      t32 = t29*t31
      t35 = t11+t28
      t41 = t35**2
      t53 = t2*t19+t5*t22+t8*t25
      t54 = t53**2
      t61 = sqrt(1-t54/t10/t27)
      t66 = t10**2
      t73 = t27**2
      sf1(1) = w1*(t17-t32/2)/t35-w1*(t11-t28)/t41*(t17+t32/2)
      sf1(2) =-w2/t61*((t13+t18+t1)*t12*t29-t53*t11/t66*t29*t15/2-t53*t1
     #2*t28/t73*t31/2)
      sf1(3) = w3*e(1)
      return
300     continue
c       c(1)=x1-x3
c       c(2)=y1-y3
c       c(3)=z1-z3
c       d(1)=x2-x3
c       d(2)=y2-y3
c       d(3)=z2-z3
c       call vcross(c,d,e)
c       call vnorm (e,dummy)
      t2 = x-x1
      t3 = t2**2
      t4 = -y1
      t5 = y+t4
      t6 = t5**2
      t8 = z-z1
      t9 = t8**2
      t10 = t3+t6+t9
      t11 = sqrt(t10)
      t12 = 1/t11
      t13 = 2*y
      t15 = t13-2*y1
      t17 = t12*t15/2
      t19 = x-x2
      t20 = t19**2
      t21 = -y2
      t22 = y+t21
      t23 = t22**2
      t25 = z-z2
      t26 = t25**2
      t27 = t20+t23+t26
      t28 = sqrt(t27)
      t29 = 1/t28
      t31 = t13-2*y2
      t32 = t29*t31
      t35 = t11+t28
      t41 = t35**2
      t53 = t2*t19+t5*t22+t8*t25
      t54 = t53**2
      t61 = sqrt(1-t54/t10/t27)
      t66 = t10**2
      t73 = t27**2
      sf1(1) = w1*(t17-t32/2)/t35-w1*(t11-t28)/t41*(t17+t32/2)
      sf1(2) =-w2/t61*((t13+t21+t4)*t12*t29-t53*t11/t66*t29*t15/2-t53*t1
     #2*t28/t73*t31/2)
      sf1(3) = w3*e(2)
      return
400     continue
c       c(1)=x1-x3
c       c(2)=y1-y3
c       c(3)=z1-z3
c       d(1)=x2-x3
c       d(2)=y2-y3
c       d(3)=z2-z3
c       call vcross(c,d,e)
c       call vnorm (e,dummy)
      t2 = x-x1
      t3 = t2**2
      t5 = y-y1
      t6 = t5**2
      t7 = -z1
      t8 = z+t7
      t9 = t8**2
      t10 = t3+t6+t9
      t11 = sqrt(t10)
      t12 = 1/t11
      t13 = 2*z
      t15 = t13-2*z1
      t17 = t12*t15/2
      t19 = x-x2
      t20 = t19**2
      t22 = y-y2
      t23 = t22**2
      t24 = -z2
      t25 = z+t24
      t26 = t25**2
      t27 = t20+t23+t26
      t28 = sqrt(t27)
      t29 = 1/t28
      t31 = t13-2*z2
      t32 = t29*t31
      t35 = t11+t28
      t41 = t35**2
      t53 = t2*t19+t5*t22+t8*t25
      t54 = t53**2
      t61 = sqrt(1-t54/t10/t27)
      t66 = t10**2
      t73 = t27**2
      sf1(1) = w1*(t17-t32/2)/t35-w1*(t11-t28)/t41*(t17+t32/2)
      sf1(2) =-w2/t61*((t13+t24+t7)*t12*t29-t53*t11/t66*t29*t15/2-t53*t1
     #2*t28/t73*t31/2)
      sf1(3) = w3*e(3)
      return
      end
**************************************************************************
      subroutine projec(xk1,yk1,zk1,xk,yk,zk,icnr,ircn,m)
      USE mdat2
      implicit double precision (a-h,o-z)
czr   common /ycur1/ cpar(100,8),nrt(100,3),ic
      common /curvedat/  cpar(20,8),nrt(20,3),ic,nbe,nn3
c.the>
czr   common /toleranz/tol
c.the<
cpf   dimension m(*),v(4),w(3),grad(3),gradr(3)
      dimension m(*)
c      if (nrt(icnr,2).eq.1) then
c        x1 = cpar(icnr,1)
c        y1 = cpar(icnr,2)
c        x2 = cpar(icnr,3)
c        y2 = cpar(icnr,4)
c        b = (y1*x2-y2*x1) / (x2-x1)
c        a = (y2-y1) / (x2-x1)
c165     continue
c        if (abs(a).lt.1.0d-10) then
c          yk=b
c        else
c          yk= (a*xk+a**2*yk+b)/(a**2+1)
c          xk= (yk-b)/a
c        endif
c        x1=xk
c        y1=yk
c        call curve (x1,y1,dummy,m,icnr,diff)
c        if(abs(y1-yk).gt.1.0e-3) then
c          goto 165
c        endif
c        return
c      elseif (nrt(icnr,2).eq.5) then
c      if (ircn.eq.0) then
c
c.... inner node
cpf>. Rückprojektion
      do 100 i=1,20                   ! max 20 curves
        if( nrt(i,1) .eq. icnr ) then
          ityp = nrt(i,2)             ! curve typ
          nn3  = nrt(i,3)
         goto 200
        else
        endif
100   continue
cpf                   write(iow,2000) icnr     ! error curve not input
cpf      if(ior.lt.0) write(*,  2000) icnr
      stop
200   continue
C
      if ( ityp .eq. 32 ) then
C.... Kugel mit dem Radius r und Mittelpunkt (x0,y0,z0)
c
        r  = cpar(icnr,1)
        x0 = cpar(icnr,2)
        y0 = cpar(icnr,3)
        z0 = cpar(icnr,4)
c
      if (((xk1-x0)**2+(yk1-y0)**2).gt.(r**2)) then
        xk1=xk
        yk1=yk
        zk1=zk
      else
        zk1=dsqrt(r**2-(xk1-x0)**2-(yk1-y0)**2)+z0
      endif
c
      elseif ( ityp .eq. 31 ) then
C.... Ebene mit der Normalen x0 y0 z0 und der Konstanten c0
c     (x0*x+y0*y+z0*z+c0=0)
c
       x0 = cpar(icnr,1)
       y0 = cpar(icnr,2)
       z0 = cpar(icnr,3)
       c0 = cpar(icnr,4)
c
       if ((dabs(y0).lt.0.00001).and.(dabs(z0).lt.0.00001)) then
         xk1=(-y0*yk1-z0*zk1-c0)/x0
       elseif  ((dabs(x0).lt.0.00001).and.(dabs(z0).lt.0.00001)) then
         yk1=(-x0*xk1-z0*zk1-c0)/y0
       else
         zk1=(-x0*xk1-y0*yk1-c0)/z0
       endif
c
      endif
      return
      end
cpf<
c
c
ccric        if (nrt(icnr,2).eq.34) then
ccpf     if (nrt(icnr,2).eq.165) then
c        if (nrt(icnr,2).eq.0) then
ccric (unvollst)
c          a=cpar(icnr,1)
c          b=cpar(icnr,2)
c          c=cpar(icnr,3)
c          d=cpar(icnr,4)
c          e=cpar(icnr,5)
c          f=cpar(icnr,6)
c          g=cpar(icnr,7)
c          h=cpar(icnr,8)
c          ri=cpar(icnr+1,1)
c          rj=cpar(icnr+1,2)
cc.the          do 660 i=1,10
c          do 660 i=1,100
c            t9x = -2*a*xk-d*yk-f*zk-g
c            t9y = -2*b*yk-d*xk-e*zk-h
c            t9z = -2*c*zk-e*yk-f*xk-ri
c            if ((abs(t9x).ge.abs(t9y)).and.(abs(t9x).ge.abs(t9z))) then
cc
cc... solve x-component-equation to lambda yk=v(1)+v(2)*xk, zk=v(3)+v(4)*xk
cc
c              ixyz = 1
c              t15 = (2*b*yk+d*xk+e*zk+h)/t9x
c              t24 = (2*c*zk+e*yk+f*xk+ri)/t9x
c              v(1) = yk+xk*t15
c              v(2) = -t15
c              v(3) = zk+xk*t24
c              v(4) = -t24
cc... enter yk and zk in geometry-function and
cc... factor out 1 (w(1)), xk (w(2)) and xk**2 (w(3))
cc
c              t1 = v(3)**2
c              t3 = v(1)**2
c              t28 = v(2)**2
c              t31 = v(4)**2
c              w(1) = c*t1+b*t3+h*v(1)+e*v(1)*v(3)+rj+ri*v(3)
c              w(2) = h*v(2)+f*v(3)+d*v(1)+e*v(1)*v(4)+e*v(2)*v(3)
c     +           +2*b*v(1)*v(2)+2*c*v(3)*v(4)+g+ri*v(4)
c              w(3) = a+d*v(2)+e*v(2)*v(4)+b*t28+f*v(4)+c*t31
cc
cC... solve reduced geometry-function
cC
c              if (abs(w(3)).gt.1.0d-10) then
c                t1 = 1/w(3)
c                t2 = -w(2)
c                t3 = w(2)**2
c                t7 = sqrt(t3-4*w(3)*w(1))
c                x1 = t1*(t2+t7)/2
c                x2 = t1*(t2-t7)/2
cc
cc... backsubstitution
cC
c                y1=v(1)+v(2)*x1
c                y2=v(1)+v(2)*x2
c                z1=v(3)+v(4)*x1
c                z2=v(3)+v(4)*x2
c                p1=(xk-x1)**2+(yk-y1)**2
c                p2=(xk-x2)**2+(yk-y2)**2
c                if(p1.lt.p2) then
c                  xk=x1
c                  yk=y1
c                  zk=z1
c                else
c                  xk=x2
c                  yk=y2
c                  zk=z2
c                endif
c              else
c                xk=-w(1)/w(2)
c                yk=v(1)+v(2)*xk
c                zk=v(3)+v(4)*xk
c              endif
c              goto 660
c           elseif((abs(t9y).ge.abs(t9x)).and.(abs(t9y).ge.abs(t9z)))then
cc
cc... solve y-component-equation to lambda xk=v(1)+v(2)*yk, zk=v(3)+v(4)*yk
cC
c              ixyz=2
c              t15 = (2*c*zk+e*yk+f*xk+ri)/t9y
c              t24 = (2*a*xk+d*yk+f*zk+g)/t9y
c              v(3) = zk+yk*t15
c              v(4) = -t15
c              v(1) = xk+yk*t24
c              v(2) = -t24
cc
cc... enter yk and zk in geometry-function and
cc... factor out 1 (w(1)), yk (w(2)) and yk**2 (w(3))
cC
c              t16 = v(1)**2
c              t21 = v(3)**2
c              t26 = v(2)**2
c              t31 = v(4)**2
c              w(1) = a*t16+g*v(1)+rj+f*v(3)*v(1)+c*t21+ri*v(3)
c              w(2) = g*v(2)+h+f*v(3)*v(2)+f*v(4)*v(1)+2*c*v(3)*v(4)
c     +           +2*a*v(1)*v(2)+d*v(1)+e*v(3)+ri*v(4)
c              w(3) = d*v(2)+a*t26+b+f*v(4)*v(2)+e*v(4)+c*t31
cc
cC... solve reduced geometry-function
cC
c              if (abs(w(3)).gt.1.0d-10) then
c                t1 = 1/w(3)
c                t2 = -w(2)
c                t3 = w(2)**2
c                t7 = sqrt(t3-4*w(3)*w(1))
c                y1 = t1*(t2+t7)/2
c                y2 = t1*(t2-t7)/2
cc
cc... backsubstitution
cC
c                x1=v(1)+v(2)*y1
c                x2=v(1)+v(2)*y2
c                z1=v(3)+v(4)*y1
c                z2=v(3)+v(4)*y2
c                p1=(xk-x1)**2+(yk-y1)**2+(zk-z1)**2
c                p2=(xk-x2)**2+(yk-y2)**2+(zk-z2)**2
c                if(p1.lt.p2) then
c                  xk=x1
c                  yk=y1
c                  zk=z1
c                else
c                  xk=x2
c                  yk=y2
c                  zk=z2
c                endif
c              else
c                yk=-w(1)/w(2)
c                xk=v(1)+v(2)*yk
c                zk=v(3)+v(4)*yk
c              endif
c              goto 660
c           elseif((abs(t9z).ge.abs(t9x)).and.(abs(t9z).ge.abs(t9y)))then
cc
cc... solve z-component-equation to lambda xk=v(1)+v(2)*zk, yk=v(3)+v(4)*zk
cC
c              ixyz=3
c              t15 = (2*a*xk+d*yk+f*zk+g)/t9z
c              t24 = (2*b*yk+d*xk+e*zk+h)/t9z
c              v(1) = xk+zk*t15
c              v(2) = -t15
c              v(3) = yk+zk*t24
c              v(4) = -t24
cc
cc... enter xk and yk in geometry-function and
cc... factor out 1 (w(1)), zk (w(2)) and zk**2 (w(3))
cC
c              t1 = v(1)**2
c              t4 = v(3)**2
c              t28 = v(2)**2
c              t31 = v(4)**2
c              w(1) = a*t1+g*v(1)+rj+b*t4+h*v(3)+d*v(1)*v(3)
c              w(2) = ri+g*v(2)+2*b*v(3)*v(4)+h*v(4)+e*v(3)
c     +           +d*v(1)*v(4)+d*v(2)*v(3)+2*a*v(1)*v(2)+f*v(1)
c              w(3) = e*v(4)+d*v(2)*v(4)+a*t28+f*v(2)+b*t31+c
cc
cC... solve reduced geometry-function
cC
c              if (abs(w(3)).gt.1.0d-10) then
c                t1 = 1/w(3)
c                t2 = -w(2)
c                t3 = w(2)**2
c                t7 = sqrt(t3-4*w(3)*w(1))
c                z1 = t1*(t2+t7)/2
c                z2 = t1*(t2-t7)/2
cc
cc... backsubstitution
cC
c                x1=v(1)+v(2)*z1
c                x2=v(1)+v(2)*z2
c                y1=v(3)+v(4)*z1
c                y2=v(3)+v(4)*z2
c                p1=(xk-x1)**2+(yk-y1)**2+(zk-z1)**2
c                p2=(xk-x2)**2+(yk-y2)**2+(zk-z2)**2
c                if(p1.lt.p2) then
c                  xk=x1
c                  yk=y1
c                  zk=z1
c                else
c                  xk=x2
c                  yk=y2
c                  zk=z2
c                endif
c              else
c                zk=-w(1)/w(2)
c                xk=v(1)+v(2)*zk
c                yk=v(3)+v(4)*zk
c              endif
c            else
c              stop
c            endif
c            call curve (xk,yk,zk,m,icnr,diff)
cc.the            if(abs(diff).lt.1.0e-12) goto 495
c            if(abs(diff).lt.0.000001) goto 495
c660       continue
c        else
cc... newton--method for arbitrary curves
c          call curve (xk,yk,zk,m,icnr,diff)
cc.the          if (abs(diff).lt.1.0e-12) goto 495
c          if (abs(diff).lt.0.000001) goto 495
cc.the          do 165 i=1,10
c          do 165 i=1,100
c            diffm1=diff
c            xkm1=xk
c            ykm1=yk
c            zkm1=zk
c            call curveg ( xk,yk,zk,icnr,grad )
c            gam=-diffm1/(grad(1)**2+grad(2)**2+grad(3)**2)
c            xk=xk+gam*grad(1)
c            yk=yk+gam*grad(2)
c            zk=zk+gam*grad(3)
c            call curve (xk,yk,zk,m,icnr,diff)
c            if (dabs(diffm1).lt.dabs(diff)) stop
cc.the          if (abs(diff).lt.1.0e-12) goto 495
c            if (abs(diff).lt.0.000001) goto 495
c165       continue
c        endif
c        stop
c      else
cc
cc... boundary node or shell intersection node
cC
c        do 200 nv=1,100
cc... check, whether functional value decreases
c          call curve (xk,yk,zk,m,icnr,f1)
c          call curve (xk,yk,zk,m,ircn,f2)
c          call curveg ( xk,yk,zk,icnr,grad )
c          call curveg ( xk,yk,zk,ircn,gradr )
c          zfo  = f1**2+f2**2
c          a11  = grad(1)**2+grad(2)**2+grad(3)**2
c          a12  = grad(1)*gradr(1)+grad(2)*gradr(2)+grad(3)*gradr(3)
c          a22  = gradr(1)**2+gradr(2)**2+gradr(3)**2
c          det = a11*a22-a12**2
c          rdet = abs(a11*a22)+a12**2
cc.the          if((abs(det).lt.1.0d-6*rdet).or.(abs(det).lt.1.0d-50))then
c          if((abs(det).lt.1.0d-7*rdet).or.(abs(det).lt.1.0d-50))then
c            stop 'stop in projec'
c          endif
c          rh1 = (f2*a12 - f1*a22)/det
c          rh2 = (f1*a12 - f2*a11)/det
c          xkm1= xk
c          ykm1= yk
c          zkm1= zk
c          xk  = xk + rh1*grad(1) + rh2*gradr(1)
c          yk  = yk + rh1*grad(2) + rh2*gradr(2)
c          zk  = zk + rh1*grad(3) + rh2*gradr(3)
c          call curve (xk,yk,zk,m,icnr,f1)
cc.the          call curve (xk,yk,zk,m,icnr,f2)
c          call curve (xk,yk,zk,m,ircn,f2)
c          zf  = f1**2 + f2**2
cc.the          if ((abs(f1).lt.1.0d-12).and.(abs(f2).lt.1.0d-12))
cc.the     +        goto 495
c          if ((abs(f1).lt.0.00000005).and.(abs(f2).lt.0.00000005))
c     +        goto 495
c          if(zfo.lt.zf) then
c            stop
c          endif
c200     continue
c        stop
c      endif
c495   continue
c      return
c      end
*************************************************************************
      subroutine simpav(xk,yk,zk,icnr,grad,ieh,rhv,ircn,m)
      implicit double precision (a-h,o-z)
c.... this subroutine locates a node quasioptimal by simple averaging
      dimension grad(3),m(*),rhv(9,8),ta(3),tb(3)
      common /curvedat/cpar(20,8),nrt1(20,3),ic,nbe,nn3
c
      x=0.0d0
       do 100 i=1,20                   ! max 20 curves
        if( nrt1(i,1) .eq. icnr ) then
          ityp = nrt1(i,2)             ! curve typ
          nn3  = nrt1(i,3)
         goto 200
        endif
100   continue
      stop
200   continue
c
      if (ityp.eq.132) then
c.... simple middle value method / inner node
              xkm1=0.0d0
              ykm1=0.0d0
              do 250 ij=1,ieh
                x1=rhv(1,ij)
                y1=rhv(2,ij)
                x2=rhv(4,ij)
                y2=rhv(5,ij)
                rh1=(x1+x2-y1+y2)/2.0d0
                rh2=(y1+y2+x1-x2)/2.0d0
                f1=(x1+x2+y1-y2)/2.0d0
                f2=(y1+y2-x1+x2)/2.0d0
                if(dsqrt((rh1-xk)**2+(rh2-yk)**2).lt.dsqrt((f1-xk)**2+
     +               (f2-yk)**2)) then
                  xkm1=xkm1+rh1
                  ykm1=ykm1+rh2
                else
                  xkm1=xkm1+f1
                  ykm1=ykm1+f2
                endif
250             continue
                xk1=xkm1/ieh
                yk1=ykm1/ieh
                zk1=0
                call projec(xk1,yk1,zk1,xk,yk,zk,icnr,ircn,m)
                xk=xk1
                yk=yk1
                zk=zk1
                goto 400
      endif
      call curveg ( xk,yk,zk,icnr,grad )
      call vnorm ( grad,dummy )
c.... projection to tangential plane
      do 345 j=1,ieh
        do 456 ii=1,2
          ij = 3*(ii-1)
          rshv = (rhv(1+ij,j)-xk)*grad(1) + (rhv(2+ij,
     +          j)-yk)*grad(2) + (rhv(3+ij,j)-zk)*grad(3)
          rhv(1+ij,j) = rhv(1+ij,j) - rshv*grad(1)
          rhv(2+ij,j) = rhv(2+ij,j) - rshv*grad(2)
          rhv(3+ij,j) = rhv(3+ij,j) - rshv*grad(3)
456     continue
345   continue
C.... 2 orthogonal tangent vectors
      ta(1)=rhv(1,1)-xk
      ta(2)=rhv(2,1)-yk
      ta(3)=rhv(3,1)-zk
      call vnorm (ta,dummy)
      call vcross(ta,grad,tb)
      xkm1=0.0d0
      ykm1=0.0d0
      do 300 j=1,ieh
c.... transformation to 2-dimensional problem
        x1=(rhv(1,j)-xk)*ta(1)+(rhv(2,j)-yk)*ta(2)+(rhv(3,j)-zk)*ta(3)
        y1=(rhv(1,j)-xk)*tb(1)+(rhv(2,j)-yk)*tb(2)+(rhv(3,j)-zk)*tb(3)
        x2=(rhv(4,j)-xk)*ta(1)+(rhv(5,j)-yk)*ta(2)+(rhv(6,j)-zk)*ta(3)
        y2=(rhv(4,j)-xk)*tb(1)+(rhv(5,j)-yk)*tb(2)+(rhv(6,j)-zk)*tb(3)
C.... 2-dimensional averaging
        rh1=(x1+x2-y1+y2)/2.0d0
        rh2=(y1+y2+x1-x2)/2.0d0
        f1=(x1+x2+y1-y2)/2.0d0
        f2=(y1+y2-x1+x2)/2.0d0
        if(dsqrt(rh1**2+rh2**2).lt.dsqrt(f1**2+f2**2)) then
          xkm1=xkm1+rh1
          ykm1=ykm1+rh2
        else
          xkm1=xkm1+f1
          ykm1=ykm1+f2
        endif
300   continue
      rh1 = xkm1/ieh
      rh2 = ykm1/ieh
      xk1 = xk + rh1*ta(1) + rh2*tb(1)
      yk1 = yk + rh1*ta(1) + rh2*tb(1)
      zk1 = zk + rh1*ta(1) + rh2*tb(1)
      call projec(xk1,yk1,zk1,xk,yk,zk,icnr,ircn,m)
      xk=xk1
      yk=yk1
      zk=zk1
400   return
      end
*************************************************************************
      subroutine curveg(xk,yk,zk,icnr,grad)
c.... p. fellmoser
      USE cdata
      USE iofile
      USE tdata
      implicit double precision (a-h,o-z)
      common /curvedat/cpar(20,8),nrt1(20,3),ic,nbe,nn3
      dimension grad(3)
c
      x=0.0d0
      do 100 i=1,20                   ! max 20 curves
        if( nrt1(i,1) .eq. icnr ) then
          ityp = nrt1(i,2)             ! curve typ
          nn3  = nrt1(i,3)
         goto 200
        else
        endif
100   continue
                   write(iow,2000) icnr     ! error curve not input
      if(ior.lt.0) write(*,  2000) icnr
      stop
200   continue

C------------------------------------------------------------------------
c                   3 - D  Kurven
C------------------------------------------------------------------------
C
      if ( ityp .eq. 32 ) then
c
C.... Kugel mit dem Radius r und Mittelpunkt (x0,y0,z0)
c
        r  = cpar(icnr,1)
        x0 = cpar(icnr,2)
        y0 = cpar(icnr,3)
        z0 = cpar(icnr,4)
c
cpf     grad(1)=-2.0d0*(yk-cpar(icnr,3))
cpf     grad(2)=2.0d0*(xk-cpar(icnr,2))
cpf     grad(3)=2.0d0*(zk-cpar(icnr,4))
        grad(1)=(xk-x0)
        grad(2)=(yk-y0)
        grad(3)=(zk-z0)
c
      elseif ( ityp .eq. 31 ) then
c
C.... Ebene mit der Normalen x0 y0 z0 und der Konstanten c0
c     (x0*x+y0*y+z0*z+c0=0)
c
        x0 = cpar(icnr,1)
        y0 = cpar(icnr,2)
        z0 = cpar(icnr,3)
        c0 = cpar(icnr,4)
c
       if ((dabs(y0).lt.0.00001).and.(dabs(z0).lt.0.00001)) then
         grad(1)=1
         grad(2)=-1
         grad(3)=0
c
       elseif ((dabs(x0).lt.0.00001).and.(dabs(z0).lt.0.00001)) then
         grad(1)=0
         grad(2)=1
         grad(3)=0
c
       else
c
cpf      grad(1)=cpar(icnr,2)
cpf      grad(2)=cpar(icnr,1)
cpf      grad(3)=-cpar(icnr,3)
         grad(1)=cpar(icnr,1)
         grad(2)=cpar(icnr,2)
         grad(3)=cpar(icnr,3)
       endif
c
      elseif ( ityp .eq. 50 ) then
c.... curve defined via points
czr     call ctyp50 (x,y,m(nbe),z,nn3)
C
c.... user defined curves
      elseif ( ityp.gt.50 .and. ityp.le.100 ) then
czr     call cvuser ( x,y,z,ityp,m,nr,diff )
      else
                     write(iow,2001) ityp
        if(ior.lt.0) write(*,2001) ityp
        stop
      endif
      return
2000  format('***ERROR*** curve number ',i3,' in not been input')
2001  format('***ERROR*** curve type ',i3,' in not defined')
      end
c
***********************************************************************
      subroutine quality(ike,ikez,x,ix,xm,ndm,lad,lne,lshel,isq,m,iael0)
      USE cdata
      USE iofile
      implicit double precision (a-h,o-z)
      dimension ike(*),ikez(*),x(ndm,*),ix(nen+4,*),xm(3),m(*),
     +      rhv(9,16),v1(3),v2(3),v3(3),v4(3),e(3),iael0(*)
      logical lad,lne,lshel
cryz
      pi = 3.14159265359
      nz = 0
      warm = 0
      iwarmax = 0
      iwarhvmax = 0
      sm = 0
      imax = 0
      ihvmax = 0
      smax = 0
      wm = 0
      iwmax = 0
      iwhvmax = 0
      wmax = 0
      if (lshel) then
      Do 200 i = 1,numnp
            xm(1)=x(1,i)
            xm(2)=x(2,i)
            xm(3)=x(3,i)
            ieh=ikez(i+1)-ikez(i)
            do 330 j=1,ieh
             ihv=ike(ikez(i)+j-1)
c.... ihv = number of element
              if (lad) then
                if (isq.ne.4) then
                  if (lne.and.(iael0(ihv).lt.0)) then
                    goto 165
                  endif
                endif
             endif
cric         if (lne.and.(iael0(ihv).lt.0)) goto 165
              do 660 k=1,nen
                if(ix(k,ihv).eq.i) then
                  ilhvh=mod(k+2,4)
                  iuhvh=mod(k,4)
                  iohvh=mod(k+1,4)
                   ilhv=ilhvh+1
                  iuhv=iuhvh+1
                  iohv=iohvh+1
                  il=ix(ilhv,ihv)
                  iu=ix(iuhv,ihv)
                  io=ix(iohv,ihv)
               endif
660           continue
cric          write coordinates of other nodes
             rhv(1,j) = x(1,il)
             rhv(2,j) = x(2,il)
             rhv(3,j) = x(3,il)
             rhv(4,j) = x(1,iu)
             rhv(5,j) = x(2,iu)
             rhv(6,j) = x(3,iu)
             rhv(7,j) = x(1,io)
             rhv(8,j) = x(2,io)
             rhv(9,j) = x(3,io)
330         continue
      do 390 j = 1,ieh
         ihv = ike(ikez(i)+(j-1))
              a = xm(1) - rhv(1,j)
              b = xm(2) - rhv(2,j)
              c = xm(3) - rhv(3,j)
              v1(1) = a
              v1(2) = b
              v1(3) = c
              sl = sqrt(a**2 + b**2 + c**2)
              a = xm(1) - rhv(4,j)
              b = xm(2) - rhv(5,j)
              c = xm(3) - rhv(6,j)
              v2(1) = a
              v2(2) = b
              v2(3) = c
              su = sqrt(a**2 + b**2 + c**2)
              v3(1) = rhv(1,j) - rhv(7,j)
              v3(2) = rhv(2,j) - rhv(8,j)
              v3(3) = rhv(3,j) - rhv(9,j)
              v4(1) = rhv(4,j) - rhv(7,j)
              v4(2) = rhv(5,j) - rhv(8,j)
              v4(3) = rhv(6,j) - rhv(9,j)
c.... Seitenl„ngenverh„ltnisse
              if ((sl/su).gt.(su/sl)) then
                 s = sl/su
              else
                 s = su/sl
              endif
c.... mittleres Seitenl„ngenverh„ltnis
              sm = sm + (s**2)
c.... maximales Seitenl„ngenverh„ltnis
              if (s.gt.smax) then
                 smax = s
                 imax = i
                 ihvmax = ihv
              endif
c.... Winkel berechnen
              alpharad = acos((vskal(v1,v2))/(sl*su))
              alpha = alpharad * (360/(2*pi))
              w = abs(90-alpha)
c.... mittleren Winkel
              wm = wm + (w**2)
c.... maximale Winkelabweichung
              if (w.gt.wmax) then
                 wmax = w
                 iwmax = i
                 iwhvmax = ihv
              endif
              nz = nz + 1
c.... warpage, nur bei Schalen
                 call vcross(v3,v4,e)
                 call vnorm(e,rnorm)
                 war = (vskal(v2,e))/(sl+su)
                 warm = warm + (war**2)
                 if (nz.eq.1) then
                    warmax = war
                 endif
                 if (war.ge.warmax) then
                    warmax = war
                    iwarmax = i
                    iwarhvmax = ihv
                 endif
390   continue
200   continue
              sm = sqrt(sm/nz)
              wm = sqrt(wm/nz)
              warm = sqrt(warm/nz)
              write(jfile,1000) sm, smax
              write(jfile,2000) wm, wmax
              if (lshel) then
                 write(jfile,3000) warm, warmax
              endif
cryz end
165   continue
      endif
cryz 2-dimensional
      If (.not.lshel) then
      Do 201 i = 1,numnp
            xm(1)=x(1,i)
            xm(2)=x(2,i)
            ieh=ikez(i+1)-ikez(i)
            do 331 j=1,ieh
             ihv=ike(ikez(i)+j-1)
c.... ihv = number of element
              if (lad) then
                if (isq.ne.4) then
                  if (lne.and.(iael0(ihv).lt.0)) then
                    goto 166
                  endif
                endif
             endif
cric              if (lne.and.(iael0(ihv).lt.0)) goto 166
              do 661 k=1,nen
                if(ix(k,ihv).eq.i) then
                  ilhvh=mod(k+2,4)
                  iuhvh=mod(k,4)
                  iohvh=mod(k+1,4)
                   ilhv=ilhvh+1
                  iuhv=iuhvh+1
                  iohv=iohvh+1
                  il=ix(ilhv,ihv)
                  iu=ix(iuhv,ihv)
                  io=ix(iohv,ihv)
               endif
661           continue
cric          write coordinates of other nodes
             rhv(1,j) = x(1,il)
             rhv(2,j) = x(2,il)
             rhv(4,j) = x(1,iu)
             rhv(5,j) = x(2,iu)
             rhv(7,j) = x(1,io)
             rhv(8,j) = x(2,io)
331         continue
      do 391 j = 1,ieh
         ihv = ike(ikez(i)+(j-1))
              a = xm(1) - rhv(1,j)
              b = xm(2) - rhv(2,j)
              v1(1) = a
              v1(2) = b
              sl = sqrt(a**2 + b**2)
              a = xm(1) - rhv(4,j)
              b = xm(2) - rhv(5,j)
              v2(1) = a
              v2(2) = b
              su = sqrt(a**2 + b**2)
              v3(1) = rhv(1,j) - rhv(7,j)
              v3(2) = rhv(2,j) - rhv(8,j)
              v4(1) = rhv(4,j) - rhv(7,j)
              v4(2) = rhv(5,j) - rhv(8,j)
c.... Seitenl„ngenverh„ltnisse
              if ((sl/su).gt.(su/sl)) then
                 s = sl/su
              else
                 s = su/sl
              endif
c.... mittleres Seitenl„ngenverh„ltnis
              sm = sm + (s**2)
c.... maximales Seitenl„ngenverh„ltnis
              if (s.gt.smax) then
                 smax = s
                 imax = i
                 ihvmax = ihv
              endif
c.... Winkel berechnen
              alpharad = acos((vskal(v1,v2))/(sl*su))
              alpha = alpharad * (360/(2*pi))
              w = abs(90-alpha)
c.... mittleren Winkel
              wm = wm + (w**2)
c.... maximale Winkelabweichung
              if (w.gt.wmax) then
                 wmax = w
                 iwmax = i
                 iwhvmax = ihv
              endif
              nz = nz + 1
391   continue
201   continue
              sm = sqrt(sm/nz)
              wm = sqrt(wm/nz)
              write(jfile,1000) sm, smax
              write(jfile,2000) wm, wmax
cryz end
166   continue
      endif
      return
1000  format ('Seitenl.-verh.: Mittel:',f12.5,' max:',f12.5)
2000  format ('Winkelabweichung: Mittel:',f12.5,' max:',f12.5)
3000  format ('warpage: Mittel:',f12.5,' max:',f12.5)
      end
c
      function f(x,y,iwf,isp)
      implicit double precision (a-h,o-z)
      common /fvc/ w1,w2,rhv(4,8),ieh
      go to(1,2,3,4,5,5,7,8,9,8,8), iwf
1     f = f1(x,y,iwf,isp)
      goto 10
2     f = f2(x,y,iwf,isp)
      goto 10
3     f = f3(x,y,iwf,isp)
      goto 10
4     f = f4(x,y,iwf,isp)
      goto 10
5     f = f5(x,y,iwf,isp)
      goto 10
7     f = f7(x,y,iwf,isp)
      goto 10
8     f = f8(x,y,iwf,isp)
      goto 10
9     f = f9(x,y,iwf,isp)
10    return
      end
C************************************************************************
      function f1(x,y,iwf,isp)
      USE iofile
      implicit double precision (a-h,o-z)
      common /fvc/ w1,w2,rhv(4,8),ieh
      e=1.0d0
      z=2.0d0
      v=4.0d0
      s=6.0d0
      a=8.0d0
c
c.....function value f (1), first derivatives df/dx (2) and df/dy (3)
C     and second derivatives d2f/dx2 (4) and d2f/dxdy (5) and d2f/dy2 (6)
C     of function f(x,y,w1,w2) =
C
C               /a**2 - b**2\ 2       /a**2 * b**2    \ 2
C           w1*| ----------- |  + w2*| ----------- - 1 |
C               \(a** + b**2/         \ (a x b)**2    /
C
      f1=0.0d0
      go to(11,12,13,14,15,16), isp
11    continue
        do 110 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
        f1 = f1 + w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2
     +         /((x-x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**2
     +        + w2*(((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)
     +         /((x-x1)*(y-y2)-(x-x2)*(y-y1))**2-e)**2
110   continue
      return
12    continue
        do 120 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
        f1 = f1 + 2.0d0*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)
     +       /((x-x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**2*(-2.0d0*x1
     +       +2.0d0*x2)-2.0d0*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)
     +       **2)**2/((x-x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3
     +       *(4.0d0*x-2.0d0*x1-2.0d0*x2)
        f1 = f1 +2.0d0*w2*(((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)
     +       /((x-x1)*(y-y2)-(x-x2)*(y-y1))**2-1.0d0)*((2.0d0*x
     +       -2.0d0*x1)*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-y2)-(x-x2)
     +       *(y-y1))**2+((x-x1)**2+(y-y1)**2)*(2.0d0*x-2*x2)/((x-x1)
     +       *(y-y2)-(x-x2)*(y-y1))**2-2*((x-x1)**2+(y-y1)**2)*((x-x2)
     +       **2+(y-y2)**2)/((x-x1)*(y-y2)-(x-x2)*(y-y1))**3*(-y2+y1))
120   continue
      return
13    continue
        do 130 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
        f1 = f1 + 2.0d0*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)
     +       /((x-x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**2*(-2.0d0*y1
     +       +z*y2)-2.0d0*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)
     +       **2/((x-x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3*(4.0d0*y
     +       -2.0d0*y1-z*y2)
        f1 = f1 + 2.0d0*w2*(((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)
     +       /((x-x1)*(y-y2)-(x-x2)*(y-y1))**2-1.0d0)*((2.0d0*y-2.0d0
     +       *y1)*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-y2)-(x-x2)*(y-y1))**2
     +       +((x-x1)**2+(y-y1)**2)*(2.0d0*y-2.0d0*y2)/((x-x1)*(y-y2)
     +       -(x-x2)*(y-y1))**2-2.0d0*((x-x1)**2+(y-y1)**2)*((x-x2)**2
     +       +(y-y2)**2)/((x-x1)*(y-y2)-(x-x2)*(y-y1))**3*(-x1+x2))
130   continue
      return
14     continue
        do 140 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      s1 = z*w1*(-z*x1+z*x2)**2/((x-x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2
     #)**2-a*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)/((x-x1)**2+(y-
     #y1)**2+(x-x2)**2+(y-y2)**2)**3*(-z*x1+z*x2)*(v*x-z*x1-z*x2)+s*w1*(
     #(x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)**2+(y-y1)**2+(
     #x-x2)**2+(y-y2)**2)**4*(v*x-z*x1-z*x2)**2
      s2 = s1-a*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)*
     #*2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3
      s3 = s2
      s5 = z*w2*((z*x-z*x1)*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-y2)-(x-x2)*
     #(y-y1))**2+((x-x1)**2+(y-y1)**2)*(z*x-z*x2)/((x-x1)*(y-y2)-(x-x2)*
     #(y-y1))**2-z*((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)/((x-x1)*(
     #y-y2)-(x-x2)*(y-y1))**3*(-y2+y1))**2
      s6 = z*w2*(((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-
     #y2)-(x-x2)*(y-y1))**2-e)*(z*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-y2)-(
     #x-x2)*(y-y1))**2+z*(z*x-z*x1)*(z*x-z*x2)/((x-x1)*(y-y2)-(x-x2)*(y-
     #y1))**2-v*(z*x-z*x1)*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-y2)-(x-x2)*(
     #y-y1))**3*(-y2+y1)+z*((x-x1)**2+(y-y1)**2)/((x-x1)*(y-y2)-(x-x2)*(
     #y-y1))**2-v*((x-x1)**2+(y-y1)**2)*(z*x-z*x2)/((x-x1)*(y-y2)-(x-x2)
     #*(y-y1))**3*(-y2+y1)+s*((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)
     #/((x-x1)*(y-y2)-(x-x2)*(y-y1))**4*(-y2+y1)**2)
      s4 = s5+s6
      f1 = f1 + s3+s4
140   continue
      return
15     continue
        do 150 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      s1 = z*w1*(-z*y1+z*y2)/((x-x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**
     #2*(-z*x1+z*x2)-v*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)/((x-
     #x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3*(-z*x1+z*x2)*(v*y-z*y1-z*
     #y2)-v*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)/((x-x1)**2+(y-y
     #1)**2+(x-x2)**2+(y-y2)**2)**3*(v*x-z*x1-z*x2)*(-z*y1+z*y2)
      s2 = s1+s*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)*
     #*2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**4*(v*x-z*x1-z*x2)*(v*y-z*y1-z*y
     #2)
      s3 = s2
      s5 = z*w2*((z*y-z*y1)*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-y2)-(x-x2)*
     #(y-y1))**2+((x-x1)**2+(y-y1)**2)*(z*y-z*y2)/((x-x1)*(y-y2)-(x-x2)*
     #(y-y1))**2-z*((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)/((x-x1)*(
     #y-y2)-(x-x2)*(y-y1))**3*(-x1+x2))*((z*x-z*x1)*((x-x2)**2+(y-y2)**2
     #)/((x-x1)*(y-y2)-(x-x2)*(y-y1))**2+((x-x1)**2+(y-y1)**2)*(z*x-z*x2
     #)/((x-x1)*(y-y2)-(x-x2)*(y-y1))**2-z*((x-x1)**2+(y-y1)**2)*((x-x2)
     #**2+(y-y2)**2)/((x-x1)*(y-y2)-(x-x2)*(y-y1))**3*(-y2+y1))
      s7 = z*w2
      s9 = ((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-y2)-(x
     #-x2)*(y-y1))**2-e
      s10 = (z*x-z*x1)*(z*y-z*y2)/((x-x1)*(y-y2)-(x-x2)*(y-y1))**2-z*(z*
     #x-z*x1)*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-y2)-(x-x2)*(y-y1))**3*(-x
     #1+x2)+(z*y-z*y1)*(z*x-z*x2)/((x-x1)*(y-y2)-(x-x2)*(y-y1))**2-z*((x
     #-x1)**2+(y-y1)**2)*(z*x-z*x2)/((x-x1)*(y-y2)-(x-x2)*(y-y1))**3*(-x
     #1+x2)-z*(z*y-z*y1)*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-y2)-(x-x2)*(y-
     #y1))**3*(-y2+y1)-z*((x-x1)**2+(y-y1)**2)*(z*y-z*y2)/((x-x1)*(y-y2)
     #-(x-x2)*(y-y1))**3*(-y2+y1)+s*((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-
     #y2)**2)/((x-x1)*(y-y2)-(x-x2)*(y-y1))**4*(-y2+y1)*(-x1+x2)
      s8 = s9*s10
      s6 = s7*s8
      s4 = s5+s6
      f1 = f1 + s3+s4
150   continue
      return
16    continue
        do 160 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      s1 = z*w1*(-z*y1+z*y2)**2/((x-x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2
     #)**2-a*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)/((x-x1)**2+(y-
     #y1)**2+(x-x2)**2+(y-y2)**2)**3*(-z*y1+z*y2)*(v*y-z*y1-z*y2)+s*w1*(
     #(x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)**2+(y-y1)**2+(
     #x-x2)**2+(y-y2)**2)**4*(v*y-z*y1-z*y2)**2
      s2 = s1-a*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)*
     #*2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3
      s3 = s2
      s5 = z*w2*((z*y-z*y1)*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-y2)-(x-x2)*
     #(y-y1))**2+((x-x1)**2+(y-y1)**2)*(2*y-z*y2)/((x-x1)*(y-y2)-(x-x2)*
     #(y-y1))**2-z*((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)/((x-x1)*(
     #y-y2)-(x-x2)*(y-y1))**3*(-x1+x2))**2
      s6 = z*w2*(((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-
     #y2)-(x-x2)*(y-y1))**2-e)*(z*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-y2)-(
     #x-x2)*(y-y1))**2+z*(z*y-z*y1)*(z*y-z*y2)/((x-x1)*(y-y2)-(x-x2)*(y-
     #y1))**2-v*(z*y-z*y1)*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-y2)-(x-x2)*(
     #y-y1))**3*(-x1+x2)+z*((x-x1)**2+(y-y1)**2)/((x-x1)*(y-y2)-(x-x2)*(
     #y-y1))**2-v*((x-x1)**2+(y-y1)**2)*(z*y-z*y2)/((x-x1)*(y-y2)-(x-x2)
     #*(y-y1))**3*(-x1+x2)+s*((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)
     #/((x-x1)*(y-y2)-(x-x2)*(y-y1))**4*(-x1+x2)**2)
      s4 = s5+s6
      f1 = f1 + s3+s4
160   continue
      fn=0.0d0
      return
      end
c***********************************************************************
      function f2(x,y,iwf,isp)
      USE iofile
      implicit double precision (a-h,o-z)
      common /fvc/ w1,w2,rhv(4,8),ieh
c.....function value f (1), first derivatives df/dx (2) and df/dy (3)
C     and second derivatives d2f/dx2 (4) and d2f/dxdy (5) and d2f/dy2 (6)
C     of function f(x,y,w1,w2) =
C
C     w1*(a**2-b**2)**2 + w2*(c**2-a**2-b**2)**2
C
      a=8.0d0
      pp=1.0d0
      f2=0.0d0
      zx4   = 0.0d0
      zx3   = 0.0d0
      zx2y2 = 0.0d0
      zx2y  = 0.0d0
      zx2   = 0.0d0
      zxy2  = 0.0d0
      zxy   = 0.0d0
      zx    = 0.0d0
      zy4   = 0.0d0
      zy3   = 0.0d0
      zy2   = 0.0d0
      zy    = 0.0d0
      z0    = 0.0d0
      do 270 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
        zx4  = zx4 + pp*4.0d0*w2
        zx3  = zx3 + pp*w2*(-8.0d0*x1-8.0d0*x2)
        zx2y2= zx2y2 + pp*8.0d0*w2
        zx2y = zx2y +pp*w2*(-8*y1-a*y2)
        zx2  = zx2 + pp*(w1*(-2.0d0*x1+2.0d0*x2)**2+w2*(-4.0d0*(x2
     +        -x1)**2-4.0d0*(y2-y1)**2+4.0d0*x1**2+4.0d0*y1**2
     +        +4.0d0*x2**2+4.0d0*y2**2+(2.0d0*x1+2.0d0*x2)**2))
        zxy2 = zxy2 - pp*4.0d0*w2*(2.0d0*x1+2.0d0*x2)
        zxy  = zxy + pp*(2.0d0*w1*(-2.0d0*y1+2.0d0*y2)*(-2.0d0*x1
     +        +2.0d0*x2)+2.0d0*w2*(2.0d0*y1+2.0d0*y2)*(2.0d0*x1
     +        +2.0d0*x2))
        zx   = zx +pp*(2.0d0*w1*(x1**2+y1**2-x2**2-y2**2)*
     +        (-2.0d0*x1+2.0d0*x2)+2.0d0*w2*((x2-x1)**2+(y2-y1)**2
     +        -x1**2-y1**2-x2**2-y2**2)*(2.0d0*x1+2.0d0*x2))
        zy4  = zy4 + pp*4*w2
        zy3  = zy3 + pp*w2*(-8.0d0*y1-8.0d0*y2)
        zy2  = zy2 + pp*(w1*(-2.0d0*y1+2.0d0*y2)**2+w2*(-4.0d0
     +        *(x2-x1)**2-4.0d0*(y2-y1)**2+4.0d0*x1**2+4.0d0*y1**2
     +        +4.0d0*x2**2+4.0d0*y2**2+(2.0d0*y1+2.0d0*y2)**2))
        zy   = zy + pp*(2.0d0*w1*(x1**2+y1**2-x2**2-y2**2)
     +        *(-2.0d0*y1+2.0d0*y2)+2.0d0*w2*((x2-x1)**2+(y2-y1)**2
     +        -x1**2-y1**2-x2**2-y2**2)*(2.0d0*y1+2.0d0*y2))
        z0   = z0 + pp*(w1*(x1**2+y1**2-x2**2-y2**2)**2+w2*((x2
     +        -x1)**2 + (y2-y1)**2-x1**2-y1**2-x2**2-y2**2)**2)
270   continue
      go to(21,22,23,24,25,26), isp
21     continue
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
        f2 = zx4*x**4 + zx3*x**3 + zx2y2*x**2*y**2 +
     +      zx2y*x**2*y + zx2*x**2 + zxy2*x*y**2 +
     +      zxy*x*y + zx*x + zy4*y**4 + zy3*y**3 +
     +      zy2*y**2 + zy*y + z0
      return
22     continue
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
        f2 = 4.0d0*zx4*x**3 + 3.0d0*zx3*x**2
     +            + 2.0d0*zx2y2*x*y**2 + 2.0d0*zx2y*x*y
     +            + 2.0d0*zx2*x + zxy2*y**2 + zxy*y + zx
      return
23     continue
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
        f2 = 2.0d0*zx2y2*x**2*y + zx2y*x**2 + 2.0d0*zxy2*x*y
     +            + zxy*x + 4.0d0*zy4*y**3 + 3.0d0*zy3*y**2
     +            + 2.0d0*zy2*y + zy
      return
24     continue
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
        f2 = 12.0d0*zx4*x**2 + 6.0d0*zx3*x + 2.0d0*zx2y2*y**2
     +            + 2.0d0*zx2y*y + 2.0d0*zx2
      return
25     continue
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
        f2 = 4.0d0*zx2y2*x*y + 2.0d0*zx2y*x + 2.0d0*zxy2*y +zxy
      return
26     continue
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
        f2 = 2.0d0*zx2y2*x**2 + 2.0d0*zxy2*x + 12.0d0*zy4*y**2
     +            + 6.0d0*zy3*y + 2.0d0*zy2
      return
      end
c*********************************************************************
      function f3(x,y,iwf,isp)
      USE iofile
      implicit double precision (a-h,o-z)
      common /fvc/ w1,w2,rhv(4,8),ieh
c
c.....function value f (1), first derivatives df/dx (2) and df/dy (3)
C     and second derivatives d2f/dx2 (4) and d2f/dxdy (5) and d2f/dy2 (6)
C     of function f(x,y,w1,w2) =
C
C               /a**2 - b**2\ 2       /a**2 * b**2    \ 2
C           w1*| ----------- |  + w2*| ----------- - 1 |
C               \(a** + b**2/         \ (a x b)**2    /
C
      f3=0.0d0
      go to(31,32,33,34,35,36), isp
31    continue
        do 310 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
        f3 = f3 + w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2
     +         /((x-x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**2
     +        + w2*(((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)
     +         /((x-x1)*(y-y2)-(x-x2)*(y-y1))**2-1)**2
310   continue
      return
32    continue
        do 320 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
        f3 = f3 + 2.0d0*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)
     +       /((x-x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**2*(-2.0d0*x1
     +       +2.0d0*x2)-2.0d0*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)
     +       **2)**2/((x-x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3
     +       *(4.0d0*x-2.0d0*x1-2.0d0*x2)
        f3 = f3 +2.0d0*w2*(((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)
     +       /((x-x1)*(y-y2)-(x-x2)*(y-y1))**2-1.0d0)*((2.0d0*x
     +       -2.0d0*x1)*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-y2)-(x-x2)
     +       *(y-y1))**2+((x-x1)**2+(y-y1)**2)*(2.0d0*x-2*x2)/((x-x1)
     +       *(y-y2)-(x-x2)*(y-y1))**2-2*((x-x1)**2+(y-y1)**2)*((x-x2)
     +       **2+(y-y2)**2)/((x-x1)*(y-y2)-(x-x2)*(y-y1))**3*(-y2+y1))
320   continue
      return
33    continue
        do 330 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
        f3 = f3 + 2.0d0*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)
     +       /((x-x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**2*(-2.0d0*y1
     +       +2*y2)-2.0d0*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)
     +       **2/((x-x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3*(4.0d0*y
     +       -2.0d0*y1-2*y2)
        f3 = f3 + 2.0d0*w2*(((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)
     +       /((x-x1)*(y-y2)-(x-x2)*(y-y1))**2-1.0d0)*((2.0d0*y-2.0d0
     +       *y1)*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-y2)-(x-x2)*(y-y1))**2
     +       +((x-x1)**2+(y-y1)**2)*(2.0d0*y-2.0d0*y2)/((x-x1)*(y-y2)
     +       -(x-x2)*(y-y1))**2-2.0d0*((x-x1)**2+(y-y1)**2)*((x-x2)**2
     +       +(y-y2)**2)/((x-x1)*(y-y2)-(x-x2)*(y-y1))**3*(-x1+x2))
330   continue
      return
34    continue
        do 340 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      s1 =
     +      2.0d0*w1*(-2.0d0*x1+2.0d0*x2)**2/((x-x1)**2+(y-y1)**2+(x-x2)
     +      **2+(y-y2)**2
     #)**2-8.0d0*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)/((x-x1)
     +      **2+(y-
     #y1)**2+(x-x2)**2+(y-y2)**2)**3*(-2.0d0*x1+2.0d0*x2)*(4.0d0*x-
     +      2.0d0*x1-2.0d0*x2)+6.0d0*w1*(
     #(x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)**2+(y-y1)**2+(
     #x-x2)**2+(y-y2)**2)**4*(4.0d0*x-2.0d0*x1-2.0d0*x2)**2
      s2 =
     +      s1-8.0d0*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/
     +      ((x-x1)*
     #*2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3
      s3 = s2
      s5 =
     +      2.0d0*w2*((2.0d0*x-2.0d0*x1)*((x-x2)**2+(y-y2)**2)/((x-x1)*
     +      (y-y2)-(x-x2)*
     #(y-y1))**2+((x-x1)**2+(y-y1)**2)*(2.0d0*x-2.0d0*x2)/((x-x1)*(y-y2)
     +      -(x-x2)*
     #(y-y1))**2-2.0d0*((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)/((x-
     +      x1)*(
     #y-y2)-(x-x2)*(y-y1))**3*(-y2+y1))**2
      s6 =
     +      2.0d0*w2*(((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)/((x-
     +      x1)*(y-
     #y2)-(x-x2)*(y-y1))**2-1.0d0)*(2.0d0*((x-x2)**2+(y-y2)**2)/((x-x1)*
     +      (y-y2)-(
     #x-x2)*(y-y1))**2+2.0d0*(2.0d0*x-2.0d0*x1)*(2.0d0*x-2.0d0*x2)/((x-
     +      x1)*(y-y2)-(x-x2)*(y-
     #y1))**2-4.0d0*(2.0d0*x-2.0d0*x1)*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-
     +      y2)-(x-x2)*(
     #y-y1))**3*(-y2+y1)+2.0d0*((x-x1)**2+(y-y1)**2)/((x-x1)*(y-y2)-(x-
     +      x2)*(
     #y-y1))**2-4.0d0*((x-x1)**2+(y-y1)**2)*(2.0d0*x-2.0d0*x2)/((x-x1)*
     +      (y-y2)-(x-x2)
     #*(y-y1))**3*(-y2+y1)+6.0d0*((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)
     +      **2)
     #/((x-x1)*(y-y2)-(x-x2)*(y-y1))**4*(-y2+y1)**2)+10.0d0
      s4 = s5+s6
      f3 = f3 + s3+s4
340   continue
      return
35    continue
        do 350 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      s1 =
     +      2.0d0*w1*(-2.0d0*y1+2.0d0*y2)/((x-x1)**2+(y-y1)**2+(x-x2)
     +      **2+(y-y2)**2)**
     #2*(-2.0d0*x1+2.0d0*x2)-4.0d0*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-
     +      y2)**2)/((x-
     #x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3*(-2.0d0*x1+2.0d0*x2)*
     +      (4.0d0*y-2.0d0*y1-2.0d0*
     #y2)-4.0d0*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)/((x-x1)**2+
     +      (y-y
     #1)**2+(x-x2)**2+(y-y2)**2)**3*(4.0d0*x-2.0d0*x1-2.0d0*x2)*(-
     +      2.0d0*y1+2.0d0*y2)
      s2 =
     +      s1+6.0d0*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/
     +      ((x-x1)*
     #*2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**4*(4.0d0*x-2.0d0*x1-2.0d0*x2)*
     +      (4.0d0*y-2.0d0*y1-2.0d0*y
     #2)
      s3 = s2
      s5 =
     +      2.0d0*w2*((2.0d0*y-2.0d0*y1)*((x-x2)**2+(y-y2)**2)/((x-x1)*
     +      (y-y2)-(x-x2)*
     #(y-y1))**2+((x-x1)**2+(y-y1)**2)*(2.0d0*y-2.0d0*y2)/((x-x1)*(y-y2)
     +      -(x-x2)*
     #(y-y1))**2-2.0d0*((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)/((x-
     +      x1)*(
     #y-y2)-(x-x2)*(y-y1))**3*(-x1+x2))*((2.0d0*x-2.0d0*x1)*((x-x2)**2+
     +      (y-y2)**2
     #)/((x-x1)*(y-y2)-(x-x2)*(y-y1))**2+((x-x1)**2+(y-y1)**2)*(2.0d0*x-
     +      2.0d0*x2
     #)/((x-x1)*(y-y2)-(x-x2)*(y-y1))**2-2.0d0*((x-x1)**2+(y-y1)**2)*
     +      ((x-x2)
     #**2+(y-y2)**2)/((x-x1)*(y-y2)-(x-x2)*(y-y1))**3*(-y2+y1))
      s7 = 2*w2
      s9 = ((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-y2)-(x
     #-x2)*(y-y1))**2-1.0d0
      s10 =
     +      (2.0d0*x-2.0d0*x1)*(2.0d0*y-2.0d0*y2)/((x-x1)*(y-y2)-(x-x2)*
     +      (y-y1))**2-2.0d0*(2.0d0*
     #x-2.0d0*x1)*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-y2)-(x-x2)*(y-y1))
     +      **3*(-x
     #1+x2)+(2.0d0*y-2.0d0*y1)*(2.0d0*x-2.0d0*x2)/((x-x1)*(y-y2)-(x-x2)*
     +      (y-y1))**2-2.0d0*((x
     #-x1)**2+(y-y1)**2)*(2.0d0*x-2.0d0*x2)/((x-x1)*(y-y2)-(x-x2)*(y-y1)
     +      )**3*(-x
     #1+x2)-2.0d0*(2.0d0*y-2.0d0*y1)*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-
     +      y2)-(x-x2)*(y-
     #y1))**3*(-y2+y1)-2.0d0*((x-x1)**2+(y-y1)**2)*(2.0d0*y-2.0d0*y2)/
     +      ((x-x1)*(y-y2)
     #-(x-x2)*(y-y1))**3*(-y2+y1)+6.0d0*((x-x1)**2+(y-y1)**2)*((x-x2)
     +      **2+(y-
     #y2)**2)/((x-x1)*(y-y2)-(x-x2)*(y-y1))**4*(-y2+y1)*(-x1+x2)
      s8 = s9*s10
      s6 = s7*s8
      s4 = s5+s6
      f3 = f3 + s3+s4
350   continue
      return
36    continue
        do 360 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      s1 =
     +      2.0d0*w1*(-2.0d0*y1+2.0d0*y2)**2/((x-x1)**2+(y-y1)**2+(x-x2)
     +      **2+(y-y2)**2
     #)**2-8.0d0*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)/((x-x1)
     +      **2+(y-
     #y1)**2+(x-x2)**2+(y-y2)**2)**3*(-2.0d0*y1+2.0d0*y2)*(4.0d0*y-
     +      2.0d0*y1-2.0d0*y2)+6.0d0*w1*(
     #(x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)**2+(y-y1)**2+(
     #x-x2)**2+(y-y2)**2)**4*(4.0d0*y-2.0d0*y1-2.0d0*y2)**2
      s2 =
     +      s1-8.0d0*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/
     +      ((x-x1)*
     #*2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3
      s3 = s2
      s5 =
     +      2.0d0*w2*((2.0d0*y-2.0d0*y1)*((x-x2)**2+(y-y2)**2)/((x-x1)*
     +      (y-y2)-(x-x2)*
     #(y-y1))**2+((x-x1)**2+(y-y1)**2)*(2.0d0*y-2.0d0*y2)/((x-x1)*(y-y2)
     +      -(x-x2)*
     #(y-y1))**2-2.0d0*((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)/((x-
     +      x1)*(
     #y-y2)-(x-x2)*(y-y1))**3*(-x1+x2))**2
      s6 =
     +      2.0d0*w2*(((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)**2)/((x-
     +      x1)*(y-
     #y2)-(x-x2)*(y-y1))**2-1.0d0)*(2.0d0*((x-x2)**2+(y-y2)**2)/((x-x1)*
     +      (y-y2)-(
     #x-x2)*(y-y1))**2+2.0d0*(2.0d0*y-2.0d0*y1)*(2.0d0*y-2.0d0*y2)/((x-
     +      x1)*(y-y2)-(x-x2)*(y-
     #y1))**2-4.0d0*(2.0d0*y-2.0d0*y1)*((x-x2)**2+(y-y2)**2)/((x-x1)*(y-
     +      y2)-(x-x2)*(
     #y-y1))**3*(-x1+x2)+2.0d0*((x-x1)**2+(y-y1)**2)/((x-x1)*(y-y2)-(x-
     +      x2)*(
     #y-y1))**2-4.0d0*((x-x1)**2+(y-y1)**2)*(2.0d0*y-2.0d0*y2)/((x-x1)*
     +      (y-y2)-(x-x2)
     #*(y-y1))**3*(-x1+x2)+6.0d0*((x-x1)**2+(y-y1)**2)*((x-x2)**2+(y-y2)
     +      **2)
     #/((x-x1)*(y-y2)-(x-x2)*(y-y1))**4*(-x1+x2)**2)
      s4 = s5+s6
      f3 = f3 + s3+s4
360   continue
      return
      end
c***********************************************************************
      function f4(x,y,iwf,isp)
      USE iofile
      implicit double precision (a-h,o-z)
      common /fvc/ w1,w2,rhv(4,8),ieh
c
c.....function value f (1), first derivatives df/dx (2) and df/dy (3)
C     and second derivatives d2f/dx2 (4) and d2f/dxdy (5) and d2f/dx2 (6)
C     of function f(x,y,w1,w2) =
C
C         / a**2-b**2 \ 2              /   a*b   \
C     w1*|  ---------  |  + w2*arccos |  -------  |
C         \ a**2+b**2 /                \ |a|*|b| /
C
      Pi=3.1415926d0
      f4=0.0d0
      go to(41,42,43,44,45,46), isp
41    continue
        do 410 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      t0 = w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)**2+(y
     #-y1)**2+(x-x2)**2+(y-y2)**2)**2+w2*acos(((x-x1)*(x-x2)+(y-y1)*(y-y
     #2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2+(y-y2)**2))**2
        f4 = f4 + t0
410   continue
      return
42    continue
        do 420 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      s1 = 2*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)/((x-x1)**2+(y-
     #y1)**2+(x-x2)**2+(y-y2)**2)**2*(-2*x1+2*x2)
      s2 = -2*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)**2
     #+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3*(4*x-2*x1-2*x2)-2*w2*acos(((x-x
     #1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2+
     #(y-y2)**2))/sqrt(1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-
     #y1)**2)/((x-x2)**2+(y-y2)**2))*((2*x-x2-x1)/sqrt((x-x1)**2+(y-y1)*
     #*2)/sqrt((x-x2)**2+(y-y2)**2)-((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((
     #(x-x1)**2+(y-y1)**2)**3)/sqrt((x-x2)**2+(y-y2)**2)*(2*x-2*x1)/2-((
     #x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)
     #**2+(y-y2)**2)**3)*(2*x-2*x2)/2)
c      f(1) = s1+s2
        f4 = f4 + s1+s2
420   continue
      return
43    continue
        do 430 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      s1 = 2*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)/((x-x1)**2+(y-
     #y1)**2+(x-x2)**2+(y-y2)**2)**2*(-2*y1+2*y2)
      s2 = -2*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)**2
     #+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3*(4*y-2*y1-2*y2)-2*w2*acos(((x-x
     #1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2+
     #(y-y2)**2))/sqrt(1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-
     #y1)**2)/((x-x2)**2+(y-y2)**2))*((2*y-y2-y1)/sqrt((x-x1)**2+(y-y1)*
     #*2)/sqrt((x-x2)**2+(y-y2)**2)-((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((
     #(x-x1)**2+(y-y1)**2)**3)/sqrt((x-x2)**2+(y-y2)**2)*(2*y-2*y1)/2-((
     #x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)
     #**2+(y-y2)**2)**3)*(2*y-2*y2)/2)
c      f(2) = s1+s2
        f4 = f4 + s1+s2
430   continue
      return
44    continue
        do 440 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      s1 = 2*w1*(-2*x1+2*x2)**2/((x-x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2
     #)**2-8*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)/((x-x1)**2+(y-
     #y1)**2+(x-x2)**2+(y-y2)**2)**3*(-2*x1+2*x2)*(4*x-2*x1-2*x2)+6*w1*(
     #(x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)**2+(y-y1)**2+(
     #x-x2)**2+(y-y2)**2)**4*(4*x-2*x1-2*x2)**2
      s2 = s1-8*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)*
     #*2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3
      s3 = s2+2*w2/(1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)
     #**2)/((x-x2)**2+(y-y2)**2))*((2*x-x2-x1)/sqrt((x-x1)**2+(y-y1)**2)
     #/sqrt((x-x2)**2+(y-y2)**2)-((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-
     #x1)**2+(y-y1)**2)**3)/sqrt((x-x2)**2+(y-y2)**2)*(2*x-2*x1)/2-((x-x
     #1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2
     #+(y-y2)**2)**3)*(2*x-2*x2)/2)**2
      s4 = s3
      s7 = w2*acos(((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**
     #2)/sqrt((x-x2)**2+(y-y2)**2))
      s9 = 1/(sqrt((1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)
     #**2)/((x-x2)**2+(y-y2)**2))**3))
      s10 = ((2*x-x2-x1)/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2+(y-y2)
     #**2)-((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)**3)/
     #sqrt((x-x2)**2+(y-y2)**2)*(2*x-2*x1)/2-((x-x1)*(x-x2)+(y-y1)*(y-y2
     #))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*x-2
     #*x2)/2)*(-2*((x-x1)*(x-x2)+(y-y1)*(y-y2))/((x-x1)**2+(y-y1)**2)/((
     #x-x2)**2+(y-y2)**2)*(2*x-x2-x1)+((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/(
     #(x-x1)**2+(y-y1)**2)**2/((x-x2)**2+(y-y2)**2)*(2*x-2*x1)+((x-x1)*(
     #x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)**2)/((x-x2)**2+(y-y2)**2
     #)**2*(2*x-2*x2))
      s8 = s9*s10
      s6 = s7*s8
      s8 = -2*w2
      s10 = acos(((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)
     #/sqrt((x-x2)**2+(y-y2)**2))
      s12 = 1/(sqrt(1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)
     #**2)/((x-x2)**2+(y-y2)**2)))
      s14 = 2/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2+(y-y2)**2)-(2*x-x
     #2-x1)/sqrt(((x-x1)**2+(y-y1)**2)**3)/sqrt((x-x2)**2+(y-y2)**2)*(2*
     #x-2*x1)-(2*x-x2-x1)/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y
     #2)**2)**3)*(2*x-2*x2)+3.0/4.0*((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((
     #(x-x1)**2+(y-y1)**2)**5)/sqrt((x-x2)**2+(y-y2)**2)*(2*x-2*x1)**2
      s13 = s14+((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)
     #**3)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*x-2*x1)*(2*x-2*x2)/2-((x-x1
     #)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)**3)/sqrt((x-x2)
     #**2+(y-y2)**2)+3.0/4.0*((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**
     #2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y2)**2)**5)*(2*x-2*x2)**2-((x-x1)*
     #(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y
     #-y2)**2)**3)
      s11 = s12*s13
      s9 = s10*s11
      s7 = s8*s9
      s5 = s6+s7
c      a(1,1) = s4+s5
      f4 = f4 + s4+s5
440   continue
      return
45    continue
        do 450 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      s1 = 2*w1*(-2*y1+2*y2)/((x-x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**
     #2*(-2*x1+2*x2)-4*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)/((x-
     #x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3*(4*x-2*x1-2*x2)*(-2*y1+2*
     #y2)-4*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)/((x-x1)**2+(y-y
     #1)**2+(x-x2)**2+(y-y2)**2)**3*(-2*x1+2*x2)*(4*y-2*y1-2*y2)
      s2 = s1+6*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)*
     #*2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**4*(4*x-2*x1-2*x2)*(4*y-2*y1-2*y
     #2)
      s4 = s2
      s6 = 2*w2
      s8 = 1/(1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)**2)/(
     #(x-x2)**2+(y-y2)**2))
      s9 = ((2*y-y2-y1)/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2+(y-y2)*
     #*2)-((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)**3)/s
     #qrt((x-x2)**2+(y-y2)**2)*(2*y-2*y1)/2-((x-x1)*(x-x2)+(y-y1)*(y-y2)
     #)/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*y-2*
     #y2)/2)*((2*x-x2-x1)/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2+(y-y2
     #)**2)-((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)**3)
     #/sqrt((x-x2)**2+(y-y2)**2)*(2*x-2*x1)/2-((x-x1)*(x-x2)+(y-y1)*(y-y
     #2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*x-
     #2*x2)/2)
      s7 = s8*s9
      s5 = s6*s7
      s3 = s4+s5
      s4 = s3
      s7 = w2*acos(((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**
     #2)/sqrt((x-x2)**2+(y-y2)**2))
      s9 = 1/(sqrt((1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)
     #**2)/((x-x2)**2+(y-y2)**2))**3))
      s10 = ((2*x-x2-x1)/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2+(y-y2)
     #**2)-((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)**3)/
     #sqrt((x-x2)**2+(y-y2)**2)*(2*x-2*x1)/2-((x-x1)*(x-x2)+(y-y1)*(y-y2
     #))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*x-2
     #*x2)/2)*(-2*((x-x1)*(x-x2)+(y-y1)*(y-y2))/((x-x1)**2+(y-y1)**2)/((
     #x-x2)**2+(y-y2)**2)*(2*y-y2-y1)+((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/(
     #(x-x1)**2+(y-y1)**2)**2/((x-x2)**2+(y-y2)**2)*(2*y-2*y1)+((x-x1)*(
     #x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)**2)/((x-x2)**2+(y-y2)**2
     #)**2*(2*y-2*y2))
      s8 = s9*s10
      s6 = s7*s8
      s8 = -2*w2
      s10 = acos(((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)
     #/sqrt((x-x2)**2+(y-y2)**2))
      s12 = 1/(sqrt(1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)
     #**2)/((x-x2)**2+(y-y2)**2)))
      s14 = -(2*x-x2-x1)/sqrt(((x-x1)**2+(y-y1)**2)**3)/sqrt((x-x2)**2+(
     #y-y2)**2)*(2*y-2*y1)/2-(2*x-x2-x1)/sqrt((x-x1)**2+(y-y1)**2)/sqrt(
     #((x-x2)**2+(y-y2)**2)**3)*(2*y-2*y2)/2-(2*y-y2-y1)/sqrt(((x-x1)**2
     #+(y-y1)**2)**3)/sqrt((x-x2)**2+(y-y2)**2)*(2*x-2*x1)/2+3.0/4.0*((x
     #-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)**5)/sqrt((x-
     #x2)**2+(y-y2)**2)*(2*x-2*x1)*(2*y-2*y1)
      s13 = s14+((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)
     #**3)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*x-2*x1)*(2*y-2*y2)/4-(2*y-y
     #2-y1)/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*
     #x-2*x2)/2+((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)
     #**3)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*x-2*x2)*(2*y-2*y1)/4+3.0/4.
     #0*((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x
     #-x2)**2+(y-y2)**2)**5)*(2*x-2*x2)*(2*y-2*y2)
      s11 = s12*s13
      s9 = s10*s11
      s7 = s8*s9
      s5 = s6+s7
c      a(1,2) = s4+s5
      f4 = f4 + s4+s5
450   continue
      return
46    continue
        do 460 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      s1 = 2*w1*(-2*y1+2*y2)**2/((x-x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2
     #)**2-8*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)/((x-x1)**2+(y-
     #y1)**2+(x-x2)**2+(y-y2)**2)**3*(-2*y1+2*y2)*(4*y-2*y1-2*y2)+6*w1*(
     #(x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)**2+(y-y1)**2+(
     #x-x2)**2+(y-y2)**2)**4*(4*y-2*y1-2*y2)**2
      s2 = s1-8*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)*
     #*2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3
      s3 = s2+2*w2/(1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)
     #**2)/((x-x2)**2+(y-y2)**2))*((2*y-y2-y1)/sqrt((x-x1)**2+(y-y1)**2)
     #/sqrt((x-x2)**2+(y-y2)**2)-((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-
     #x1)**2+(y-y1)**2)**3)/sqrt((x-x2)**2+(y-y2)**2)*(2*y-2*y1)/2-((x-x
     #1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2
     #+(y-y2)**2)**3)*(2*y-2*y2)/2)**2
      s4 = s3
      s7 = w2*acos(((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**
     #2)/sqrt((x-x2)**2+(y-y2)**2))
      s9 = 1/(sqrt((1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)
     #**2)/((x-x2)**2+(y-y2)**2))**3))
      s10 = ((2*y-y2-y1)/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2+(y-y2)
     #**2)-((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)**3)/
     #sqrt((x-x2)**2+(y-y2)**2)*(2*y-2*y1)/2-((x-x1)*(x-x2)+(y-y1)*(y-y2
     #))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*y-2
     #*y2)/2)*(-2*((x-x1)*(x-x2)+(y-y1)*(y-y2))/((x-x1)**2+(y-y1)**2)/((
     #x-x2)**2+(y-y2)**2)*(2*y-y2-y1)+((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/(
     #(x-x1)**2+(y-y1)**2)**2/((x-x2)**2+(y-y2)**2)*(2*y-2*y1)+((x-x1)*(
     #x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)**2)/((x-x2)**2+(y-y2)**2
     #)**2*(2*y-2*y2))
      s8 = s9*s10
      s6 = s7*s8
      s8 = -2*w2
      s10 = acos(((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)
     #/sqrt((x-x2)**2+(y-y2)**2))
      s12 = 1/(sqrt(1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)
     #**2)/((x-x2)**2+(y-y2)**2)))
      s14 = 2/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2+(y-y2)**2)-(2*y-y
     #2-y1)/sqrt(((x-x1)**2+(y-y1)**2)**3)/sqrt((x-x2)**2+(y-y2)**2)*(2*
     #y-2*y1)-(2*y-y2-y1)/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y
     #2)**2)**3)*(2*y-2*y2)+3.0/4.0*((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((
     #(x-x1)**2+(y-y1)**2)**5)/sqrt((x-x2)**2+(y-y2)**2)*(2*y-2*y1)**2
      s13 = s14+((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)
     #**3)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*y-2*y1)*(2*y-2*y2)/2-((x-x1
     #)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)**3)/sqrt((x-x2)
     #**2+(y-y2)**2)+3.0/4.0*((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**
     #2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y2)**2)**5)*(2*y-2*y2)**2-((x-x1)*
     #(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y
     #-y2)**2)**3)
      s11 = s12*s13
      s9 = s10*s11
      s7 = s8*s9
      s5 = s6+s7
c      a(2,2) = s4+s5
      f4 = f4 + s4+s5
460   continue
      return
      end
c***********************************************************************
      function f5(x,y,iwf,isp)
      USE iofile
      implicit double precision (a-h,o-z)
      common /fvc/ w1,w2,rhv(4,8),ieh
c
c.....function value f (1), first derivatives df/dx (2) and df/dy (3)
C     and second derivatives d2f/dx2 (4) and d2f/dxdy (5) and d2f/dx2 (6)
C     of function f(x,y,w1,w2) =
C
C     (a**2-rk*c**2)**2 + (b**2-rk*c**2)**2
C
C     rk = (sin(w/2))**2 , w = desired (averanged) angle
C
      if (iwf.eq.5) then
        w=0.0d0
        do 570 i=1,ieh
          x1=rhv(1,i)
          y1=rhv(2,i)
          x2=rhv(3,i)
          y2=rhv(4,i)
          a1=x-x1
          a2=y-y1
          b1=x-x2
          b2=y-y2
          w=w+dacos((a1*b1+a2*b2)/dsqrt((a1**2+a2**2)*(b1**2+b2**2)))
570     continue
        rk=(dsin(w/(2.0d0*ieh)))**2
      elseif (iwf.eq.6) then
        rk = 0.5d0
      else
        write(*,1000)
        write(jfile,1000)
        stop
      endif
      Pi=3.1415926d0
      f5=0.0d0
      go to(51,52,53,54,55,56), isp
51    continue
        do 510 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      t0 =((x-x1)**2+(y-y1)**2-rk*((x2-x1)**2+(y2-y1)**2))**2+((x-x2)**2
     #+(y-y2)**2-rk*((x2-x1)**2+(y2-y1)**2))**2
        f5 = f5 + t0
510   continue
      return
52    continue
        do 520 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      t0 = 2*((x-x1)**2+(y-y1)**2-rk*((x2-x1)**2+(y2-y1)**2))*(2*x-2*x1
     #)+2*((x-x2)**2+(y-y2)**2-rk*((x2-x1)**2+(y2-y1)**2))*(2*x-2*x2)
        f5 = f5 + t0
520   continue
      return
53    continue
        do 530 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      t0 = 2*((x-x1)**2+(y-y1)**2-rk*((x2-x1)**2+(y2-y1)**2))*(2*y-2*y1
     #)+2*((x-x2)**2+(y-y2)**2-rk*((x2-x1)**2+(y2-y1)**2))*(2*y-2*y2)
        f5 = f5 + t0
530   continue
      return
54    continue
        do 540 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      f5=f5 + 2*(2*x-2*x1)**2+4*(x-x1)**2+4*(y-y1)**2-8*rk*((x2-x1)**2+(
     #y2-y1)**2)+2*(2*x-2*x2)**2+4*(x-x2)**2+4*(y-y2)**2
540   continue
      return
55    continue
        do 550 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      f5 = f5 + 2*(2*y-2*y1)*(2*x-2*x1)+2*(2*y-2*y2)*(2*x-2*x2)
550   continue
      return
56    continue
        do 560 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      f5=f5 + 2*(2*y-2*y1)**2+4*(x-x1)**2+4*(y-y1)**2-8*rk*((x2-x1)**2+(
     #y2-y1)**2)+2*(2*y-2*y2)**2+4*(x-x2)**2+4*(y-y2)**2
560   continue
      return
1000  format('error in function f : parameter iwf out of range')
      end
c***********************************************************************
      function f7(x,y,iwf,isp)
      USE iofile
      implicit double precision (a-h,o-z)
      common /fvc/ w1,w2,rhv(4,8),ieh
c
c.....function value f (1), first derivatives df/dx (2) and df/dy (3)
C     and second derivatives d2f/dx2 (4) and d2f/dxdy (5) and d2f/dx2 (6)
C     of function f(x,y,w1,w2) =
C
C         / a**2-b**2 \ 2              /   a*b     Pi \2
C     w1*|  ---------  |  + w2*arccos |  ------- - --  |
C         \ a**2+b**2 /                \ |a|*|b|   2  /
C
      Pi=3.1415926d0
      f7=0.0d0
      go to(71,72,73,74,75,76), isp
71    continue
        do 710 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      t0 = w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)**2+(y
     #-y1)**2+(x-x2)**2+(y-y2)**2)**2+w2*(acos(((x-x1)*(x-x2)+(y-y1)*(y-
     #y2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2+(y-y2)**2))-Pi/2)**2
        f7 = f7 + t0
710   continue
      return
72    continue
        do 720 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      s1 = 2*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)/((x-x1)**2+(y-
     #y1)**2+(x-x2)**2+(y-y2)**2)**2*(-2*x1+2*x2)
      s2 = -2*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)**2
     #+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3*(4*x-2*x1-2*x2)-2*w2*(acos(((x-
     #x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2
     #+(y-y2)**2))-Pi/2)/sqrt(1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)
     #**2+(y-y1)**2)/((x-x2)**2+(y-y2)**2))*((2*x-x2-x1)/sqrt((x-x1)**2+
     #(y-y1)**2)/sqrt((x-x2)**2+(y-y2)**2)-((x-x1)*(x-x2)+(y-y1)*(y-y2))
     #/sqrt(((x-x1)**2+(y-y1)**2)**3)/sqrt((x-x2)**2+(y-y2)**2)*(2*x-2*x
     #1)/2-((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(
     #((x-x2)**2+(y-y2)**2)**3)*(2*x-2*x2)/2)
c      f(1) = s1+s2
        f7 = f7 + s1+s2
720   continue
      return
73    continue
        do 730 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      s1 = 2*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)/((x-x1)**2+(y-
     #y1)**2+(x-x2)**2+(y-y2)**2)**2*(-2*y1+2*y2)
      s2 = -2*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)**2
     #+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3*(4*y-2*y1-2*y2)-2*w2*(acos(((x-
     #x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2
     #+(y-y2)**2))-Pi/2)/sqrt(1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)
     #**2+(y-y1)**2)/((x-x2)**2+(y-y2)**2))*((2*y-y2-y1)/sqrt((x-x1)**2+
     #(y-y1)**2)/sqrt((x-x2)**2+(y-y2)**2)-((x-x1)*(x-x2)+(y-y1)*(y-y2))
     #/sqrt(((x-x1)**2+(y-y1)**2)**3)/sqrt((x-x2)**2+(y-y2)**2)*(2*y-2*y
     #1)/2-((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(
     #((x-x2)**2+(y-y2)**2)**3)*(2*y-2*y2)/2)
c      f(2) = s1+s2
        f7 = f7 + s1+s2
730   continue
      return
74    continue
        do 740 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      s1 = 2*w1*(-2*x1+2*x2)**2/((x-x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2
     #)**2-8*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)/((x-x1)**2+(y-
     #y1)**2+(x-x2)**2+(y-y2)**2)**3*(-2*x1+2*x2)*(4*x-2*x1-2*x2)+6*w1*(
     #(x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)**2+(y-y1)**2+(
     #x-x2)**2+(y-y2)**2)**4*(4*x-2*x1-2*x2)**2
      s2 = s1-8*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)*
     #*2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3
      s3 = s2+2*w2/(1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)
     #**2)/((x-x2)**2+(y-y2)**2))*((2*x-x2-x1)/sqrt((x-x1)**2+(y-y1)**2)
     #/sqrt((x-x2)**2+(y-y2)**2)-((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-
     #x1)**2+(y-y1)**2)**3)/sqrt((x-x2)**2+(y-y2)**2)*(2*x-2*x1)/2-((x-x
     #1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2
     #+(y-y2)**2)**3)*(2*x-2*x2)/2)**2
      s4 = s3
      s7 = w2*(acos(((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)*
     #*2)/sqrt((x-x2)**2+(y-y2)**2))-Pi/2)
      s9 = 1/(sqrt((1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)
     #**2)/((x-x2)**2+(y-y2)**2))**3))
      s10 = ((2*x-x2-x1)/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2+(y-y2)
     #**2)-((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)**3)/
     #sqrt((x-x2)**2+(y-y2)**2)*(2*x-2*x1)/2-((x-x1)*(x-x2)+(y-y1)*(y-y2
     #))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*x-2
     #*x2)/2)*(-2*((x-x1)*(x-x2)+(y-y1)*(y-y2))/((x-x1)**2+(y-y1)**2)/((
     #x-x2)**2+(y-y2)**2)*(2*x-x2-x1)+((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/(
     #(x-x1)**2+(y-y1)**2)**2/((x-x2)**2+(y-y2)**2)*(2*x-2*x1)+((x-x1)*(
     #x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)**2)/((x-x2)**2+(y-y2)**2
     #)**2*(2*x-2*x2))
      s8 = s9*s10
      s6 = s7*s8
      s8 = -2*w2
      s10 = acos(((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)
     #/sqrt((x-x2)**2+(y-y2)**2))-Pi/2
      s12 = 1/(sqrt(1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)
     #**2)/((x-x2)**2+(y-y2)**2)))
      s14 = 2/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2+(y-y2)**2)-(2*x-x
     #2-x1)/sqrt(((x-x1)**2+(y-y1)**2)**3)/sqrt((x-x2)**2+(y-y2)**2)*(2*
     #x-2*x1)-(2*x-x2-x1)/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y
     #2)**2)**3)*(2*x-2*x2)+3.0/4.0*((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((
     #(x-x1)**2+(y-y1)**2)**5)/sqrt((x-x2)**2+(y-y2)**2)*(2*x-2*x1)**2
      s13 = s14+((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)
     #**3)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*x-2*x1)*(2*x-2*x2)/2-((x-x1
     #)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)**3)/sqrt((x-x2)
     #**2+(y-y2)**2)+3.0/4.0*((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**
     #2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y2)**2)**5)*(2*x-2*x2)**2-((x-x1)*
     #(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y
     #-y2)**2)**3)
      s11 = s12*s13
      s9 = s10*s11
      s7 = s8*s9
      s5 = s6+s7
c      a(1,1) = s4+s5
      f7 = f7 + s4+s5
740   continue
      return
75    continue
        do 750 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      s1 = 2*w1*(-2*y1+2*y2)/((x-x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**
     #2*(-2*x1+2*x2)-4*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)/((x-
     #x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3*(4*x-2*x1-2*x2)*(-2*y1+2*
     #y2)-4*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)/((x-x1)**2+(y-y
     #1)**2+(x-x2)**2+(y-y2)**2)**3*(-2*x1+2*x2)*(4*y-2*y1-2*y2)
      s2 = s1+6*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)*
     #*2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**4*(4*x-2*x1-2*x2)*(4*y-2*y1-2*y
     #2)
      s4 = s2
      s6 = 2*w2
      s8 = 1/(1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)**2)/(
     #(x-x2)**2+(y-y2)**2))
      s9 = ((2*y-y2-y1)/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2+(y-y2)*
     #*2)-((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)**3)/s
     #qrt((x-x2)**2+(y-y2)**2)*(2*y-2*y1)/2-((x-x1)*(x-x2)+(y-y1)*(y-y2)
     #)/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*y-2*
     #y2)/2)*((2*x-x2-x1)/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2+(y-y2
     #)**2)-((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)**3)
     #/sqrt((x-x2)**2+(y-y2)**2)*(2*x-2*x1)/2-((x-x1)*(x-x2)+(y-y1)*(y-y
     #2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*x-
     #2*x2)/2)
      s7 = s8*s9
      s5 = s6*s7
      s3 = s4+s5
      s4 = s3
      s7 = w2*(acos(((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)*
     #*2)/sqrt((x-x2)**2+(y-y2)**2))-Pi/2)
      s9 = 1/(sqrt((1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)
     #**2)/((x-x2)**2+(y-y2)**2))**3))
      s10 = ((2*x-x2-x1)/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2+(y-y2)
     #**2)-((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)**3)/
     #sqrt((x-x2)**2+(y-y2)**2)*(2*x-2*x1)/2-((x-x1)*(x-x2)+(y-y1)*(y-y2
     #))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*x-2
     #*x2)/2)*(-2*((x-x1)*(x-x2)+(y-y1)*(y-y2))/((x-x1)**2+(y-y1)**2)/((
     #x-x2)**2+(y-y2)**2)*(2*y-y2-y1)+((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/(
     #(x-x1)**2+(y-y1)**2)**2/((x-x2)**2+(y-y2)**2)*(2*y-2*y1)+((x-x1)*(
     #x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)**2)/((x-x2)**2+(y-y2)**2
     #)**2*(2*y-2*y2))
      s8 = s9*s10
      s6 = s7*s8
      s8 = -2*w2
      s10 = acos(((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)
     #/sqrt((x-x2)**2+(y-y2)**2))-Pi/2
      s12 = 1/(sqrt(1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)
     #**2)/((x-x2)**2+(y-y2)**2)))
      s14 = -(2*x-x2-x1)/sqrt(((x-x1)**2+(y-y1)**2)**3)/sqrt((x-x2)**2+(
     #y-y2)**2)*(2*y-2*y1)/2-(2*x-x2-x1)/sqrt((x-x1)**2+(y-y1)**2)/sqrt(
     #((x-x2)**2+(y-y2)**2)**3)*(2*y-2*y2)/2-(2*y-y2-y1)/sqrt(((x-x1)**2
     #+(y-y1)**2)**3)/sqrt((x-x2)**2+(y-y2)**2)*(2*x-2*x1)/2+3.0/4.0*((x
     #-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)**5)/sqrt((x-
     #x2)**2+(y-y2)**2)*(2*x-2*x1)*(2*y-2*y1)
      s13 = s14+((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)
     #**3)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*x-2*x1)*(2*y-2*y2)/4-(2*y-y
     #2-y1)/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*
     #x-2*x2)/2+((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)
     #**3)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*x-2*x2)*(2*y-2*y1)/4+3.0/4.
     #0*((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x
     #-x2)**2+(y-y2)**2)**5)*(2*x-2*x2)*(2*y-2*y2)
      s11 = s12*s13
      s9 = s10*s11
      s7 = s8*s9
      s5 = s6+s7
c      a(1,2) = s4+s5
      f7 = f7 + s4+s5
750   continue
      return
76    continue
        do 760 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      s1 = 2*w1*(-2*y1+2*y2)**2/((x-x1)**2+(y-y1)**2+(x-x2)**2+(y-y2)**2
     #)**2-8*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)/((x-x1)**2+(y-
     #y1)**2+(x-x2)**2+(y-y2)**2)**3*(-2*y1+2*y2)*(4*y-2*y1-2*y2)+6*w1*(
     #(x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)**2+(y-y1)**2+(
     #x-x2)**2+(y-y2)**2)**4*(4*y-2*y1-2*y2)**2
      s2 = s1-8*w1*((x-x1)**2+(y-y1)**2-(x-x2)**2-(y-y2)**2)**2/((x-x1)*
     #*2+(y-y1)**2+(x-x2)**2+(y-y2)**2)**3
      s3 = s2+2*w2/(1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)
     #**2)/((x-x2)**2+(y-y2)**2))*((2*y-y2-y1)/sqrt((x-x1)**2+(y-y1)**2)
     #/sqrt((x-x2)**2+(y-y2)**2)-((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-
     #x1)**2+(y-y1)**2)**3)/sqrt((x-x2)**2+(y-y2)**2)*(2*y-2*y1)/2-((x-x
     #1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2
     #+(y-y2)**2)**3)*(2*y-2*y2)/2)**2
      s4 = s3
      s7 = w2*(acos(((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)*
     #*2)/sqrt((x-x2)**2+(y-y2)**2))-Pi/2)
      s9 = 1/(sqrt((1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)
     #**2)/((x-x2)**2+(y-y2)**2))**3))
      s10 = ((2*y-y2-y1)/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2+(y-y2)
     #**2)-((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)**3)/
     #sqrt((x-x2)**2+(y-y2)**2)*(2*y-2*y1)/2-((x-x1)*(x-x2)+(y-y1)*(y-y2
     #))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*y-2
     #*y2)/2)*(-2*((x-x1)*(x-x2)+(y-y1)*(y-y2))/((x-x1)**2+(y-y1)**2)/((
     #x-x2)**2+(y-y2)**2)*(2*y-y2-y1)+((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/(
     #(x-x1)**2+(y-y1)**2)**2/((x-x2)**2+(y-y2)**2)*(2*y-2*y1)+((x-x1)*(
     #x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)**2)/((x-x2)**2+(y-y2)**2
     #)**2*(2*y-2*y2))
      s8 = s9*s10
      s6 = s7*s8
      s8 = -2*w2
      s10 = acos(((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)
     #/sqrt((x-x2)**2+(y-y2)**2))-Pi/2
      s12 = 1/(sqrt(1-((x-x1)*(x-x2)+(y-y1)*(y-y2))**2/((x-x1)**2+(y-y1)
     #**2)/((x-x2)**2+(y-y2)**2)))
      s14 = 2/sqrt((x-x1)**2+(y-y1)**2)/sqrt((x-x2)**2+(y-y2)**2)-(2*y-y
     #2-y1)/sqrt(((x-x1)**2+(y-y1)**2)**3)/sqrt((x-x2)**2+(y-y2)**2)*(2*
     #y-2*y1)-(2*y-y2-y1)/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y
     #2)**2)**3)*(2*y-2*y2)+3.0/4.0*((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((
     #(x-x1)**2+(y-y1)**2)**5)/sqrt((x-x2)**2+(y-y2)**2)*(2*y-2*y1)**2
      s13 = s14+((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)
     #**3)/sqrt(((x-x2)**2+(y-y2)**2)**3)*(2*y-2*y1)*(2*y-2*y2)/2-((x-x1
     #)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)**2+(y-y1)**2)**3)/sqrt((x-x2)
     #**2+(y-y2)**2)+3.0/4.0*((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**
     #2+(y-y1)**2)/sqrt(((x-x2)**2+(y-y2)**2)**5)*(2*y-2*y2)**2-((x-x1)*
     #(x-x2)+(y-y1)*(y-y2))/sqrt((x-x1)**2+(y-y1)**2)/sqrt(((x-x2)**2+(y
     #-y2)**2)**3)
      s11 = s12*s13
      s9 = s10*s11
      s7 = s8*s9
      s5 = s6+s7
c      a(2,2) = s4+s5
      f7 = f7 + s4+s5
760   continue
      return
      end
**************************************************************************
      function f8(x,y,iwf,isp)
      USE iofile
      implicit double precision (a-h,o-z)
      common /fvc/ w1,w2,rhv(4,8),ieh
c
c.....function value f (1), first derivatives df/dx (2) and df/dy (3)
C     and second derivatives d2f/dx2 (4) and d2f/dxdy (5) and d2f/dx2 (6)
C     of function f(x,y,w1,w2) =
C
C         / a-b \ 2              /   a*b     Pi \ 2
C     w1*|  ---  |  + w2*arccos |  ------- - --  |
C         \ a+b /                \ |a|*|b|   2  /
C
      Pi=3.1415926d0
      f8=0.0d0
      go to(71,72,73,74,75,76), isp
71    continue
        do 710 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      t2 = x-x1
      t3 = t2**2
      t5 = y-y1
      t6 = t5**2
      t8 = sqrt(t3+t6)
      t10 = x-x2
      t11 = t10**2
      t13 = y-y2
      t14 = t13**2
      t16 = sqrt(t11+t14)
      t19 = (t8-t16)**2
      t21 = (t8+t16)**2
      t35 = (acos((t2*t10+t5*t13)/t8/t16)-Pi/2)**2
      t37 = w1*t19/t21+w2*t35
        f8 = f8 + t37
710   continue
      return
72    continue
        do 720 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      t1 = -x1
      t2 = x+t1
      t3 = t2**2
      t5 = y-y1
      t6 = t5**2
      t7 = t3+t6
      t8 = sqrt(t7)
      t9 = -x2
      t10 = x+t9
      t11 = t10**2
      t13 = y-y2
      t14 = t13**2
      t15 = t11+t14
      t16 = sqrt(t15)
      t18 = t8-t16
      t19 = t8+t16
      t20 = t19**2
      t22 = 1/t8
      t23 = 2*x
      t25 = t23-2*x1
      t27 = t22*t25/2
      t28 = 1/t16
      t30 = t23-2*x2
      t31 = t28*t30
      t38 = t18**2
      t49 = t2*t10+t5*t13
      t50 = t22*t28
      t55 = t49**2
      t62 = sqrt(1-t55/t7/t15)
      t66 = t7**2
      t73 = t15**2
      t85 = 2*w1*t18/t20*(t27-t31/2)-2*w1*t38/t20/t19*(t27+t31/2)-2*w2*(
     #acos(t49*t50)-Pi/2)/t62*((t23+t9+t1)*t50-t49*t8/t66*t28*t25/2-t49*
     #t22*t16/t73*t30/2)
        f8 = f8 + t85
720   continue
      return
73    continue
        do 730 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      t2 = x-x1
      t3 = t2**2
      t4 = -y1
      t5 = y+t4
      t6 = t5**2
      t7 = t3+t6
      t8 = sqrt(t7)
      t10 = x-x2
      t11 = t10**2
      t12 = -y2
      t13 = y+t12
      t14 = t13**2
      t15 = t11+t14
      t16 = sqrt(t15)
      t18 = t8-t16
      t19 = t8+t16
      t20 = t19**2
      t22 = 1/t8
      t23 = 2*y
      t25 = t23-2*y1
      t27 = t22*t25/2
      t28 = 1/t16
      t30 = t23-2*y2
      t31 = t28*t30
      t38 = t18**2
      t49 = t2*t10+t5*t13
      t50 = t22*t28
      t55 = t49**2
      t62 = sqrt(1-t55/t7/t15)
      t66 = t7**2
      t73 = t15**2
      t85 = 2*w1*t18/t20*(t27-t31/2)-2*w1*t38/t20/t19*(t27+t31/2)-2*w2*(
     #acos(t49*t50)-Pi/2)/t62*((t23+t12+t4)*t50-t49*t8/t66*t28*t25/2-t49
     #*t22*t16/t73*t30/2)
        f8 = f8 + t85
730   continue
      return
74    continue
        do 740 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      t1 = -x1
      t2 = x+t1
      t3 = t2**2
      t5 = y-y1
      t6 = t5**2
      t7 = t3+t6
      t8 = sqrt(t7)
      t9 = 1/t8
      t10 = 2*x
      t12 = t10-2*x1
      t14 = t9*t12/2
      t15 = -x2
      t16 = x+t15
      t17 = t16**2
      t19 = y-y2
      t20 = t19**2
      t21 = t17+t20
      t22 = sqrt(t21)
      t23 = 1/t22
      t25 = t10-2*x2
      t26 = t23*t25
      t28 = t14-t26/2
      t29 = t28**2
      t30 = t8+t22
      t31 = t30**2
      t32 = 1/t31
      t37 = t8-t22
      t39 = 1/t31/t30
      t41 = t14+t26/2
      t47 = t7**2
      t48 = 1/t47
      t49 = t8*t48
      t50 = t12**2
      t52 = -t49*t50/4
      t53 = t21**2
      t54 = 1/t53
      t55 = t22*t54
      t56 = t25**2
      t57 = t55*t56
      t65 = t37**2
      t66 = t31**2
      t68 = t41**2
      t81 = t2*t16+t5*t19
      t82 = t81**2
      t83 = 1/t7
      t84 = 1/t21
      t88 = 1-t82*t83*t84
      t90 = t10+t15+t1
      t91 = t9*t23
      t94 = t49*t23*t12
      t98 = t9*t55*t25
      t101 = t90*t91-t81*t94/2-t81*t98/2
      t102 = t101**2
      t109 = acos(t81*t91)-Pi/2
      t110 = sqrt(t88)
      t111 = t88**2
      t165 = 2*w1*t29*t32-8*w1*t37*t39*t28*t41+2*w1*t37*t32*(t52+t9+t57/
     #4-t23)+6*w1*t65/t66*t68-2*w1*t65*t39*(t52+t9-t57/4+t23)+2*w2/t88*t
     #102+w2*t109*t110/t111*t101*(-2*t81*t83*t84*t90+t82*t48*t84*t12+t82
     #*t83*t54*t25)-2*w2*t109/t110*(2*t91-t90*t94-t90*t98+3.0/4.0*t81*t8
     #/t47/t7*t23*t50+t81*t49*t55*t12*t25/2-t81*t49*t23+3.0/4.0*t81*t9*t
     #22/t53/t21*t56-t81*t9*t55)
      f8 = f8 + t165
740   continue
      return
75    continue
        do 750 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      t1 = -x1
      t2 = x+t1
      t3 = t2**2
      t4 = -y1
      t5 = y+t4
      t6 = t5**2
      t7 = t3+t6
      t8 = sqrt(t7)
      t9 = 1/t8
      t10 = 2*y
      t12 = t10-2*y1
      t14 = t9*t12/2
      t15 = -x2
      t16 = x+t15
      t17 = t16**2
      t18 = -y2
      t19 = y+t18
      t20 = t19**2
      t21 = t17+t20
      t22 = sqrt(t21)
      t23 = 1/t22
      t25 = t10-2*y2
      t26 = t23*t25
      t28 = t14-t26/2
      t29 = t8+t22
      t30 = t29**2
      t31 = 1/t30
      t32 = 2*x
      t34 = t32-2*x1
      t36 = t9*t34/2
      t38 = t32-2*x2
      t39 = t23*t38
      t41 = t36-t39/2
      t47 = t8-t22
      t49 = 1/t30/t29
      t51 = t14+t26/2
      t57 = t7**2
      t58 = 1/t57
      t59 = t8*t58
      t60 = t34*t12
      t62 = -t59*t60/4
      t63 = t21**2
      t64 = 1/t63
      t65 = t22*t64
      t66 = t38*t25
      t67 = t65*t66
      t75 = t36+t39/2
      t81 = t47**2
      t82 = t30**2
      t97 = t2*t16+t5*t19
      t98 = t97**2
      t99 = 1/t7
      t100 = 1/t21
      t104 = 1-t98*t99*t100
      t106 = t10+t18+t4
      t107 = t9*t23
      t110 = t59*t23*t12
      t114 = t9*t65*t25
      t118 = t32+t15+t1
      t121 = t59*t23*t34
      t125 = t9*t65*t38
      t128 = t118*t107-t97*t121/2-t97*t125/2
      t136 = acos(t97*t107)-Pi/2
      t137 = sqrt(t104)
      t138 = t104**2
      t194 = 2*w1*t28*t31*t41-4*w1*t47*t49*t41*t51+2*w1*t47*t31*(t62+t67
     #/4)-4*w1*t47*t49*t75*t28+6*w1*t81/t82*t75*t51-2*w1*t81*t49*(t62-t6
     #7/4)+2*w2/t104*(t106*t107-t97*t110/2-t97*t114/2)*t128+w2*t136*t137
     #/t138*t128*(-2*t97*t99*t100*t106+t98*t58*t100*t12+t98*t99*t64*t25)
     #-2*w2*t136/t137*(-t118*t110/2-t118*t114/2-t106*t121/2+3.0/4.0*t97*
     #t8/t57/t7*t23*t60+t97*t59*t65*t34*t25/4-t106*t125/2+t97*t59*t65*t3
     #8*t12/4+3.0/4.0*t97*t9*t22/t63/t21*t66)
      f8 = f8 + t194
750   continue
      return
76    continue
        do 760 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      t2 = x-x1
      t3 = t2**2
      t4 = -y1
      t5 = y+t4
      t6 = t5**2
      t7 = t3+t6
      t8 = sqrt(t7)
      t9 = 1/t8
      t10 = 2*y
      t12 = t10-2*y1
      t14 = t9*t12/2
      t16 = x-x2
      t17 = t16**2
      t18 = -y2
      t19 = y+t18
      t20 = t19**2
      t21 = t17+t20
      t22 = sqrt(t21)
      t23 = 1/t22
      t25 = t10-2*y2
      t26 = t23*t25
      t28 = t14-t26/2
      t29 = t28**2
      t30 = t8+t22
      t31 = t30**2
      t32 = 1/t31
      t37 = t8-t22
      t39 = 1/t31/t30
      t41 = t14+t26/2
      t47 = t7**2
      t48 = 1/t47
      t49 = t8*t48
      t50 = t12**2
      t52 = -t49*t50/4
      t53 = t21**2
      t54 = 1/t53
      t55 = t22*t54
      t56 = t25**2
      t57 = t55*t56
      t65 = t37**2
      t66 = t31**2
      t68 = t41**2
      t81 = t2*t16+t5*t19
      t82 = t81**2
      t83 = 1/t7
      t84 = 1/t21
      t88 = 1-t82*t83*t84
      t90 = t10+t18+t4
      t91 = t9*t23
      t94 = t49*t23*t12
      t98 = t9*t55*t25
      t101 = t90*t91-t81*t94/2-t81*t98/2
      t102 = t101**2
      t109 = acos(t81*t91)-Pi/2
      t110 = sqrt(t88)
      t111 = t88**2
      t165 = 2*w1*t29*t32-8*w1*t37*t39*t28*t41+2*w1*t37*t32*(t52+t9+t57/
     #4-t23)+6*w1*t65/t66*t68-2*w1*t65*t39*(t52+t9-t57/4+t23)+2*w2/t88*t
     #102+w2*t109*t110/t111*t101*(-2*t81*t83*t84*t90+t82*t48*t84*t12+t82
     #*t83*t54*t25)-2*w2*t109/t110*(2*t91-t90*t94-t90*t98+3.0/4.0*t81*t8
     #/t47/t7*t23*t50+t81*t49*t55*t12*t25/2-t81*t49*t23+3.0/4.0*t81*t9*t
     #22/t53/t21*t56-t81*t9*t55)
      f8 = f8 + t165
760   continue
      return
      end
**************************************************************************
      function f9(x,y,iwf,isp)
      USE iofile
      implicit double precision (a-h,o-z)
      common /fvc/ w1,w2,rhv(4,8),ieh
c
c.....function value f (1), first derivatives df/dx (2) and df/dy (3)
C     and second derivatives d2f/dx2 (4) and d2f/dxdy (5) and d2f/dx2 (6)
C     of function f(x,y,w1,w2) =
C
C                4
C         ( a-b )               /   a*b     Pi \ 2
C      w1* ----- 2 + w2*arccos |  ------- - --  |
C         ( a*b )               \ |a|*|b|   2  /
C
      Pi=3.1415926d0
      f9=0.0d0
      go to(71,72,73,74,75,76), isp
71    continue
        do 710 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      t2 = x-x1
      t3 = t2**2
      t5 = y-y1
      t6 = t5**2
      t7 = t3+t6
      t8 = sqrt(t7)
      t10 = x-x2
      t11 = t10**2
      t13 = y-y2
      t14 = t13**2
      t15 = t11+t14
      t16 = sqrt(t15)
      t19 = (t8-t16)**2
      t20 = t19**2
      t36 = (acos((t2*t10+t5*t13)/t8/t16)-Pi/2)**2
      t38 = w1*t20/t7/t15+w2*t36
        f9 = f9 + t38
710   continue
      return
72    continue
        do 720 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      t1 = -x1
      t2 = x+t1
      t3 = t2**2
      t5 = y-y1
      t6 = t5**2
      t7 = t3+t6
      t8 = sqrt(t7)
      t9 = -x2
      t10 = x+t9
      t11 = t10**2
      t13 = y-y2
      t14 = t13**2
      t15 = t11+t14
      t16 = sqrt(t15)
      t18 = t8-t16
      t19 = t18**2
      t21 = 1/t7
      t22 = 1/t15
      t23 = 1/t8
      t24 = 2*x
      t26 = t24-2*x1
      t29 = 1/t16
      t31 = t24-2*x2
      t40 = t19**2
      t41 = t7**2
      t42 = 1/t41
      t48 = t15**2
      t49 = 1/t48
      t57 = t2*t10+t5*t13
      t58 = t23*t29
      t63 = t57**2
      t68 = sqrt(1-t63*t21*t22)
      t87 = 4*w1*t19*t18*t21*t22*(t23*t26/2-t29*t31/2)-w1*t40*t42*t22*t2
     #6-w1*t40*t21*t49*t31-2*w2*(acos(t57*t58)-Pi/2)/t68*((t24+t9+t1)*t5
     #8-t57*t8*t42*t29*t26/2-t57*t23*t16*t49*t31/2)
        f9 = f9 + t87
720   continue
      return
73    continue
        do 730 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      t2 = x-x1
      t3 = t2**2
      t4 = -y1
      t5 = y+t4
      t6 = t5**2
      t7 = t3+t6
      t8 = sqrt(t7)
      t10 = x-x2
      t11 = t10**2
      t12 = -y2
      t13 = y+t12
      t14 = t13**2
      t15 = t11+t14
      t16 = sqrt(t15)
      t18 = t8-t16
      t19 = t18**2
      t21 = 1/t7
      t22 = 1/t15
      t23 = 1/t8
      t24 = 2*y
      t26 = t24-2*y1
      t29 = 1/t16
      t31 = t24-2*y2
      t40 = t19**2
      t41 = t7**2
      t42 = 1/t41
      t48 = t15**2
      t49 = 1/t48
      t57 = t2*t10+t5*t13
      t58 = t23*t29
      t63 = t57**2
      t68 = sqrt(1-t63*t21*t22)
      t87 = 4*w1*t19*t18*t21*t22*(t23*t26/2-t29*t31/2)-w1*t40*t42*t22*t2
     #6-w1*t40*t21*t49*t31-2*w2*(acos(t57*t58)-Pi/2)/t68*((t24+t12+t4)*t
     #58-t57*t8*t42*t29*t26/2-t57*t23*t16*t49*t31/2)
        f9 = f9 + t87
730   continue
      return
74    continue
        do 740 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      t1 = -x1
      t2 = x+t1
      t3 = t2**2
      t5 = y-y1
      t6 = t5**2
      t7 = t3+t6
      t8 = sqrt(t7)
      t9 = -x2
      t10 = x+t9
      t11 = t10**2
      t13 = y-y2
      t14 = t13**2
      t15 = t11+t14
      t16 = sqrt(t15)
      t18 = t8-t16
      t19 = t18**2
      t20 = 1/t7
      t21 = 1/t15
      t22 = 1/t8
      t23 = 2*x
      t25 = t23-2*x1
      t28 = 1/t16
      t30 = t23-2*x2
      t33 = t22*t25/2-t28*t30/2
      t34 = t33**2
      t40 = t19*t18
      t41 = t7**2
      t42 = 1/t41
      t49 = t15**2
      t50 = 1/t49
      t57 = t8*t42
      t58 = t25**2
      t61 = t16*t50
      t62 = t30**2
      t72 = t19**2
      t74 = 1/t41/t7
      t80 = t25*t30
      t92 = 1/t49/t15
      t104 = t2*t10+t5*t13
      t105 = t104**2
      t109 = 1-t105*t20*t21
      t111 = t23+t9+t1
      t112 = t22*t28
      t115 = t57*t28*t25
      t119 = t22*t61*t30
      t122 = t111*t112-t104*t115/2-t104*t119/2
      t123 = t122**2
      t130 = acos(t104*t112)-Pi/2
      t131 = sqrt(t109)
      t132 = t109**2
      t182 = 12*w1*t19*t20*t21*t34-8*w1*t40*t42*t21*t33*t25-8*w1*t40*t20
     #*t50*t33*t30+4*w1*t40*t20*t21*(-t57*t58/4+t22+t61*t62/4-t28)+2*w1*
     #t72*t74*t21*t58+2*w1*t72*t42*t50*t80-2*w1*t72*t42*t21+2*w1*t72*t20
     #*t92*t62-2*w1*t72*t20*t50+2*w2/t109*t123+w2*t130*t131/t132*t122*(-
     #2*t104*t20*t21*t111+t105*t42*t21*t25+t105*t20*t50*t30)-2*w2*t130/t
     #131*(2*t112-t111*t115-t111*t119+3.0/4.0*t104*t8*t74*t28*t58+t104*t
     #57*t61*t80/2-t104*t57*t28+3.0/4.0*t104*t22*t16*t92*t62-t104*t22*t6
     #1)
      f9 = f9 + t182
740   continue
      return
75    continue
        do 750 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      t1 = -x1
      t2 = x+t1
      t3 = t2**2
      t4 = -y1
      t5 = y+t4
      t6 = t5**2
      t7 = t3+t6
      t8 = sqrt(t7)
      t9 = -x2
      t10 = x+t9
      t11 = t10**2
      t12 = -y2
      t13 = y+t12
      t14 = t13**2
      t15 = t11+t14
      t16 = sqrt(t15)
      t18 = t8-t16
      t19 = t18**2
      t20 = 1/t7
      t21 = 1/t15
      t22 = 1/t8
      t23 = 2*x
      t25 = t23-2*x1
      t28 = 1/t16
      t30 = t23-2*x2
      t33 = t22*t25/2-t28*t30/2
      t34 = 2*y
      t36 = t34-2*y1
      t40 = t34-2*y2
      t43 = t22*t36/2-t28*t40/2
      t50 = t19*t18
      t51 = t7**2
      t52 = 1/t51
      t59 = t15**2
      t60 = 1/t59
      t67 = t8*t52
      t68 = t25*t36
      t71 = t16*t60
      t72 = t30*t40
      t87 = t19**2
      t89 = 1/t51/t7
      t96 = t25*t40
      t107 = t30*t36
      t113 = 1/t59/t15
      t121 = t2*t10+t5*t13
      t122 = t121**2
      t126 = 1-t122*t20*t21
      t128 = t34+t12+t4
      t129 = t22*t28
      t132 = t67*t28*t36
      t136 = t22*t71*t40
      t140 = t23+t9+t1
      t143 = t67*t28*t25
      t147 = t22*t71*t30
      t150 = t140*t129-t121*t143/2-t121*t147/2
      t158 = acos(t121*t129)-Pi/2
      t159 = sqrt(t126)
      t160 = t126**2
      t211 = 12*w1*t19*t20*t21*t33*t43-4*w1*t50*t52*t21*t33*t36-4*w1*t50
     #*t20*t60*t33*t40+4*w1*t50*t20*t21*(-t67*t68/4+t71*t72/4)-4*w1*t50*
     #t52*t21*t25*t43+2*w1*t87*t89*t21*t68+w1*t87*t52*t60*t96-4*w1*t50*t
     #20*t60*t30*t43+w1*t87*t52*t60*t107+2*w1*t87*t20*t113*t72+2*w2/t126
     #*(t128*t129-t121*t132/2-t121*t136/2)*t150+w2*t158*t159/t160*t150*(
     #-2*t121*t20*t21*t128+t122*t52*t21*t36+t122*t20*t60*t40)-2*w2*t158/
     #t159*(-t140*t132/2-t140*t136/2-t128*t143/2+3.0/4.0*t121*t8*t89*t28
     #*t68+t121*t67*t71*t96/4-t128*t147/2+t121*t67*t71*t107/4+3.0/4.0*t1
     #21*t22*t16*t113*t72)
      f9 = f9 + t211
750   continue
      return
76    continue
        do 760 i=1,ieh
        x1=rhv(1,i)
        y1=rhv(2,i)
        x2=rhv(3,i)
        y2=rhv(4,i)
      t2 = x-x1
      t3 = t2**2
      t4 = -y1
      t5 = y+t4
      t6 = t5**2
      t7 = t3+t6
      t8 = sqrt(t7)
      t10 = x-x2
      t11 = t10**2
      t12 = -y2
      t13 = y+t12
      t14 = t13**2
      t15 = t11+t14
      t16 = sqrt(t15)
      t18 = t8-t16
      t19 = t18**2
      t20 = 1/t7
      t21 = 1/t15
      t22 = 1/t8
      t23 = 2*y
      t25 = t23-2*y1
      t28 = 1/t16
      t30 = t23-2*y2
      t33 = t22*t25/2-t28*t30/2
      t34 = t33**2
      t40 = t19*t18
      t41 = t7**2
      t42 = 1/t41
      t49 = t15**2
      t50 = 1/t49
      t57 = t8*t42
      t58 = t25**2
      t61 = t16*t50
      t62 = t30**2
      t72 = t19**2
      t74 = 1/t41/t7
      t80 = t25*t30
      t92 = 1/t49/t15
      t104 = t2*t10+t5*t13
      t105 = t104**2
      t109 = 1-t105*t20*t21
      t111 = t23+t12+t4
      t112 = t22*t28
      t115 = t57*t28*t25
      t119 = t22*t61*t30
      t122 = t111*t112-t104*t115/2-t104*t119/2
      t123 = t122**2
      t130 = acos(t104*t112)-Pi/2
      t131 = sqrt(t109)
      t132 = t109**2
      t182 = 12*w1*t19*t20*t21*t34-8*w1*t40*t42*t21*t33*t25-8*w1*t40*t20
     #*t50*t33*t30+4*w1*t40*t20*t21*(-t57*t58/4+t22+t61*t62/4-t28)+2*w1*
     #t72*t74*t21*t58+2*w1*t72*t42*t50*t80-2*w1*t72*t42*t21+2*w1*t72*t20
     #*t92*t62-2*w1*t72*t20*t50+2*w2/t109*t123+w2*t130*t131/t132*t122*(-
     #2*t104*t20*t21*t111+t105*t42*t21*t25+t105*t20*t50*t30)-2*w2*t130/t
     #131*(2*t112-t111*t115-t111*t119+3.0/4.0*t104*t8*t74*t28*t58+t104*t
     #57*t61*t80/2-t104*t57*t28+3.0/4.0*t104*t22*t16*t92*t62-t104*t22*t6
     #1)
      f9 = f9 + t182
760   continue
      return
      end
**************************************************************************
      subroutine g(x,y,x1,y1,x2,y2,iwf,isp,fv)
      USE iofile
      implicit double precision (a-h,o-z)
      common /fvc/ w1,w2,rhv(4,8),ieh
      dimension fv(2)
c
c.....function value f (1), and first derivatives df/dx (2) and df/dy (3)
C     of functions
C
C                      / a-b \ 2
C     f1(x,y,w1) = w1*|  ---  |
C                      \ a+b /
C
C                             /   a*b     Pi \ 2
C     f2(x,y,w2) = w2*arccos |  ------- - --  |
C                             \ |a|*|b|   2  /
C
      Pi=3.1415926d0
      go to(11,12,13), isp
11    continue
      t2 = x-x1
      t3 = t2**2
      t5 = y-y1
      t6 = t5**2
      t8 = sqrt(t3+t6)
      t10 = x-x2
      t11 = t10**2
      t13 = y-y2
      t14 = t13**2
      t16 = sqrt(t11+t14)
      fv(1) = (w1*(t8-t16)/(t8+t16))
      fv(2) = (w2*(acos((t2*t10+t5*t13)/t8/t16)-Pi/2))
      return
12    continue
      t1 = -x1
      t2 = x+t1
      t3 = t2**2
      t5 = y-y1
      t6 = t5**2
      t7 = t3+t6
      t8 = sqrt(t7)
      t9 = 1/t8
      t10 = 2*x
      t12 = t10-2*x1
      t14 = t9*t12/2
      t15 = -x2
      t16 = x+t15
      t17 = t16**2
      t19 = y-y2
      t20 = t19**2
      t21 = t17+t20
      t22 = sqrt(t21)
      t23 = 1/t22
      t25 = t10-2*x2
      t26 = t23*t25
      t29 = t8+t22
      t35 = t29**2
      t46 = t2*t16+t5*t19
      t47 = t46**2
      t54 = sqrt(1-t47/t7/t21)
      t59 = t7**2
      t66 = t21**2
      fv(1) = w1*(t14-t26/2)/t29-w1*(t8-t22)/t35*(t14+t26/2)
      fv(2) = -w2/t54*((t10+t15+t1)*t9*t23-t46*t8/t59*t23*t12/2-t46*t9*t
     #22/t66*t25/2)
      return
13    continue
      t2 = x-x1
      t3 = t2**2
      t4 = -y1
      t5 = y+t4
      t6 = t5**2
      t7 = t3+t6
      t8 = sqrt(t7)
      t9 = 1/t8
      t10 = 2*y
      t12 = t10-2*y1
      t14 = t9*t12/2
      t16 = x-x2
      t17 = t16**2
      t18 = -y2
      t19 = y+t18
      t20 = t19**2
      t21 = t17+t20
      t22 = sqrt(t21)
      t23 = 1/t22
      t25 = t10-2*y2
      t26 = t23*t25
      t29 = t8+t22
      t35 = t29**2
      t46 = t2*t16+t5*t19
      t47 = t46**2
      t54 = sqrt(1-t47/t7/t21)
      t59 = t7**2
      t66 = t21**2
      fv(1) = w1*(t14-t26/2)/t29-w1*(t8-t22)/t35*(t14+t26/2)
      fv(2) = -w2/t54*((t10+t18+t4)*t9*t23-t46*t8/t59*t23*t12/2-t46*t9*t
     #22/t66*t25/2)
      return
      end
**************************************************************************
