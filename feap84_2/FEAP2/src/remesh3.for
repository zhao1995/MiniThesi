c ....................................................................SR.wboun
      subroutine wboun(nboue,ndf,inpf)
c ---------------------------------------------------------------------------
c.... write new boundary conditions
c ---------------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      common /gener/     kk,km,ke
      dimension nboue(ndf,*)
c....
      write(iow,999)'boun'
      do ni=1,kk
        nbsum = 0
        do nf=1,ndf
          nbsum = nbsum+nboue(nf,ni)
        enddo
        if (nbsum.ne.0) then
          write(iow,1001)  ni , (nboue(in,ni),in=1,ndf)
        endif
      enddo
999   format(/,A)
1001  format ( I4, ',0,', 6( I2 ,','))
      end
c................................................................end.SR.wboun
c....................................................................SR.wload
      subroutine wload(f,ndf,inpf)
c ---------------------------------------------------------------------------
c.... distribute surface loads
c ---------------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
czr      dimension x(ndm,*),f(ndf,*)
      dimension  f(ndf,*)
      common /gener/  kk,km,ke
czrc.... write forces for input file
czr      write(iow,999)'forc'
czr      do ikn = 1,kk
czr        if(ndf.eq.2) then
czr          if(abs(f(1,ikn)).gt.0.d0.or.abs(f(2,ikn)).gt.0.0d0)then
czr            write(iow,1001) ikn, (f(ifhg,ikn), ifhg=1,ndf)
czr          endif
czr        elseif(ndf.eq.3) then
czr          if(abs(f(1,ikn)).gt.0.d0.or.abs(f(2,ikn)).gt.0.d0
czr     +                          .or.abs(f(3,ikn)).gt.0.0d0)then
czr            write(iow,1001) ikn, (f(ifhg,ikn), ifhg=1,ndf)
czr          endif
czr        endif
czr      enddo
czr999   format(/,A)
czr1001  format ( I4, ',0,', 7(F9.4 ,','))
      end
c..............................................................end..SR.wload
c.................................................................SR.calload
      subroutine calload(x,iek,ike,ikz,f,ndm,ndf)
c ---------------------------------------------------------------------------
c.... distribute surface loads
c ---------------------------------------------------------------------------
      USE cdata
      USE iofile
      implicit double precision (a-h,o-z)
      common /curvedat/  cpar(20,8),nrt1(20,3),ic,nbe,nn3
c     common /gener/ numnp0,km,numel0
      common /gener/  kk,km,ke
      dimension iek(7,*),ike(*),ikz(*)
      dimension x(ndm,*),f(ndf,*)
      dimension m0(2),xl(10)
c.... calculate surface loads
      do 101 i=1,kk                                 !loop nodal points
        do il = 1,ndm
          xl(il) = x(il,i)
        enddo
       do 102 iif=1,ic-1                                   !loop functions
c.... load functions  nrt1(*,3)>30
        if(nrt1(iif,3).gt.30.and.nrt1(iif,3).lt.31+ndf)then
          if(ndm.lt.3)then
c.... 2-dimensional
            call curve(xl(1),y0,xl(2),m0,nrt1(iif,1),0.d0)
            if(abs(xl(2)-y0).lt.0.001)then
              xl0=dsqrt((cpar(iif,3)-cpar(iif,1))**2
     1                 +(cpar(iif,4)-cpar(iif,2))**2)
              xl1=dsqrt((xl(1)-cpar(iif,1))**2+(xl(2)-cpar(iif,2))**2)
              xl2=dsqrt((xl(1)-cpar(iif,3))**2+(xl(2)-cpar(iif,4))**2)
              pxx=(cpar(iif,7)*xl2+cpar(iif,8)*xl1)/(xl1+xl2)
              if((xl1.gt.xl0).or.(xl2.gt.xl0).or.
     1           (xl1.gt.xl0.and.xl2.gt.xl0)) goto 102
czr111196              call pload1 (iek,ike,ikz,f,x,x0,
czr111196     1                   ndm,ndf,nen,i+1,iif,pxx,xl0)
              call pload1 (iek,ike,ikz,f,x,
     1                   ndm,ndf,nen,i  ,iif,pxx,xl0)   ! warn****
            endif
          elseif(ndm.ge.3)then
c.... 3-dimensional
            call curve(xl(1),xl(2),xl(3),m0,nrt1(iif,1),diff)
czr         if(abs(diff).gt.0.0001)goto102
            if(abs(diff).lt.0.0001)then
              xl0=dsqrt((cpar(iif,4)-cpar(iif,1))**2
     1                 +(cpar(iif,5)-cpar(iif,2))**2
     2                 +(cpar(iif,6)-cpar(iif,3))**2)
              xl1=dsqrt((xl(1)-cpar(iif,1))**2+(xl(2)-cpar(iif,2))**2+
     1                  (xl(3)-cpar(iif,3))**2)
              xl2=dsqrt((xl(1)-cpar(iif,4))**2+(xl(2)-cpar(iif,5))**2+
     1                 (xl(3)-cpar(iif,6))**2)
              pxx= (cpar(iif,7)*xl2+cpar(iif,8)*xl1)/(xl1+xl2)
czr           if(nrt1(iif,2).eq.42) goto 103
              if((xl1.gt.xl0).or.(xl2.gt.xl0).or.
     1               (xl1.gt.xl0.and.xl2.gt.xl0)) goto 102
czr103        continue
czr  1        call pload2 (m(n9),m(n9a),m(n9b),m(n10),m(n8)
czr  2                        ,ndm,ndf,nen,i+1,iif,pxx,xl0)
              call pload2 (iek,  ike   ,ikz   ,f     ,x
     1                        ,ndm,ndf,nen,i  ,iif,pxx,xl0)
            endif
          endif
        endif
102    continue                                       !enddo functions
101   continue                                         !enddo nodal points
czrc.... write forces for input file
czr      write(iow,999)'forc'
czr      do ikn = 1,kk
czr        if(ndf.eq.2) then
czr          if(abs(f(1,ikn)).gt.0.d0.or.abs(f(2,ikn)).gt.0.0d0)then
czr            write(iow,1001) ikn, (f(ifhg,ikn), ifhg=1,ndf)
czr          endif
czr        elseif(ndf.eq.3) then
czr          if(abs(f(1,ikn)).gt.0.d0.or.abs(f(2,ikn)).gt.0.d0
czr     +                          .or.abs(f(3,ikn)).gt.0.0d0)then
czr            write(iow,1001) ikn, (f(ifhg,ikn), ifhg=1,ndf)
czr          endif
czr        endif
czr      enddo
czr999   format(/,A)
czr1001  format ( I4, ',0,', 7(F9.4 ,','))
      end
c..............................................................end.SR.calload
c..................................................................SR.calboun
      subroutine calboun(id,x0,nboue,ndm,ndf)
c ---------------------------------------------------------------------------
c.... calculate boundary conditions for newly generated nodes
c ---------------------------------------------------------------------------
      USE cdata
      USE iofile
      implicit double precision (a-h,o-z)
      integer id(ndf,*)
      common /curvedat/  cpar(20,8),nrt1(20,3),ic,nbe,nn3
      common /gener/     kk,km,ke
      dimension x0(ndm,*),nboue(ndf,*)
      dimension m0(2),xl(10),nbou(6)
c....
czr   write(iow,999)'boun'
      call pzeroi(nboue,ndf*kk)
c.... write old boundary into input file
      do ni=1,numnp
        do nf=1,ndf
          if (id(nf,ni).lt.0) then
            nbou(nf) = 1
            nboue(nf,ni) = 1
          else
            nbou(nf) = 0
c            nboue(nf,ni) = 0
          endif
        enddo
czr        nbsum = 0
czr        do isum=1,ndf
czr          nbsum = nbsum+nbou(isum)
czr        enddo
czr        if (nbsum.ne.0) then
czr          write(iow,1001)  ni , (nbou(in),in=1,ndf)
czr        endif
      enddo
c.... find and write new boundary
      do 101 i=1,(kk-numnp)                               ! loop np
        call pzeroi(nbou,6)
        call xcor1(xl,x0(1,i),ndm)
        do 102 iif=1,ic-1                                 ! loop curves
          if(nrt1(iif,3).gt.20.and.nrt1(iif,3).lt.21+ndf)then
            if(ndm.lt.3) then
c.... 2-dimensional
czr           call curve(xl(1),y0,0.d0,m0,nrt1(iif,1),0.d0)
              call curve(xl(1),y0,xl(2),m0,nrt1(iif,1),0.d0)
              if(abs(xl(2)-y0).lt.0.001)then
                xl0=dsqrt((cpar(iif,3)-cpar(iif,1))**2
     1                 +(cpar(iif,4)-cpar(iif,2))**2)
               xl1=dsqrt((xl(1)-cpar(iif,1))**2+(xl(2)-cpar(iif,2))**2)
               xl2=dsqrt((xl(1)-cpar(iif,3))**2+(xl(2)-cpar(iif,4))**2)
                if((xl1.gt.xl0).or.(xl2.gt.xl0).or.
     1            (xl1.gt.xl0.and.xl2.gt.xl0)) goto 102
c.... set boundary for dof
                nbou(nrt1(iif,3)-20) = 1
                nboue(nrt1(iif,3)-20,i+numnp) = 1
              endif
            elseif(ndm.ge.3) then
c.... 3-dimensional
              call curve(xl(1),xl(2),xl(3),m0,nrt1(iif,1),diff)
czr           if (abs(diff).gt.0.0001) goto102
              if(abs(diff).lt.0.0001)then
                if(nrt1(iif,2).eq.31) goto 103
                xl0=dsqrt((cpar(iif,4)-cpar(iif,1))**2
     1                 +(cpar(iif,5)-cpar(iif,2))**2
     2                 +(cpar(iif,6)-cpar(iif,3))**2)
                xl1=dsqrt((xl(1)-cpar(iif,1))**2+(xl(2)
     +                   -cpar(iif,2))**2+(xl(3)-cpar(iif,3))**2)
                xl2=dsqrt((xl(1)-cpar(iif,4))**2+(xl(2)
     +                   -cpar(iif,5))**2+(xl(3)-cpar(iif,6))**2)
czr             if(nrt1(iif,2).eq.42) goto 103
                if((xl1.gt.xl0).or.(xl2.gt.xl0).or.
     1               (xl1.gt.xl0.and.xl2.gt.xl0)) goto102
103             continue
                ntemp=i+numnp
c                write(*,*)  'b node 3D c/node',iif,ntemp
c.... set boundary for dof
                nbou(nrt1(iif,3)-20) = 1
                nboue(nrt1(iif,3)-20,i+numnp) = 1
              endif
            endif
          endif
102     continue
czrc.... write new boundary to input file
czr      nbsum = 0
czr      do isum=1,ndf
czr        nbsum = nbsum+nbou(isum)
czr      enddo
czr      if (nbsum.ne.0) then
czr        write(iow,1001)  i+numnp, (nbou(in),in=1,ndf)
czr      endif
101   continue
czr999   format(/,A)
czr1001  format ( I4, ',0,', 6( I2 ,','))
      end
c..............................................................end.SR.calboun
c...................................................................SR.pload1
czru  subroutine pload1(iek,ike,ikz,f,x,x0,
      subroutine pload1(iek,ike,ikz,f,x,
     1                  ndm,ndf,nen,i,iif,px1,xl0)
c ---------------------------------------------------------------------------
c.... Streckenlasten auf Knoten fuer 2D-Modelle
c ---------------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension iek(7,*),ike(*),ikz(*),f(ndf,*),x(ndm,*)
czru  dimension x0(ndm,*)
      dimension m0(2)
      common /curvedat/  cpar(20,8),nrt1(20,3),ic,nbe,nn3
c
c
      ist=0
c.... first element
      do 100 ii=ikz(i),ikz(i+1)-1
        iel=(ike(ii))
        do 200 iii=1,nen
          if(iek(iii,iel).eq.i)then
            ii1=iii-1
            if(ii1.eq.0) ii1=4
            i1=iek(ii1,iel)
czr         call curve(x(1,i1),y0,0.d0,   m0,nrt1(iif,1),0.d0)
            call curve(x(1,i1),y0,x(2,i1),m0,nrt1(iif,1),0.d0)
            if(abs(x(2,i1)-y0).lt.0.001)then
              ist=1
              goto 300
            else
              ii1=iii+1
              if(ii1.eq.5) ii1=1
              i1=iek(ii1,iel)
czr           call curve(x(1,i1),y0,0.d0,m0,nrt1(iif,1),0.d0)
              call curve(x(1,i1),y0,x(2,i1),m0,nrt1(iif,1),0.d0)
              if(abs(x(2,i1)-y0).lt.0.001)then
                ist=1
                goto 300
              endif
            endif
          endif
200     continue
100   continue
300   continue
      if(ist.ne.1) write(*,*)'stop in pload1'
      if(ist.ne.1) write(iow,*)'stop in pload1'
      if(ist.ne.1) stop
      xl1=dsqrt((x(1,i1)-x(1,i))**2+(x(2,i1)-x(2,i))**2)
      xl11=dsqrt((x(1,i1)-cpar(iif,1))**2+(x(2,i1)-cpar(iif,2))**2)
      xl12=dsqrt((x(1,i1)-cpar(iif,3))**2+(x(2,i1)-cpar(iif,4))**2)
c....
      nx=nrt1(iif,3)-30
      if((xl11.gt.xl0).or.(xl12.gt.xl0).or.
     1           (xl11.gt.xl0.and.xl12.gt.xl0))then
        xl=xl11
        if(xl11.gt.xl12)xl=xl12
        xl1=xl1-xl
        px2=cpar(iif,7)
        if(xl11.gt.xl12)px2=cpar(iif,8)
        f(nx,i)=f(nx,i)+xl1*(px1+px2)/2.d0
czr          write(iow,999)'forc  *first element  .gt.'
czr          write(iow,1000) i,(fx(i,inx),i=1,ndm)
czr          write(iow,*) i,(fx(i,inx),i=1,ndm)
czr          write(iow,*) nx,i,f(nx,i)
czr          write(*,*) nx,i,f(nx,i)
      else
        px2=(cpar(iif,7)*xl12+cpar(iif,8)*xl11)/(xl11+xl12)
        f(nx,i)=f(nx,i)+xl1*(px1/3.d0+px2/6.d0)
czr          write(iow,999)'forc  *first element  .le.'
czr          write(iow,1000) i,(fx(i,inx),i=1,ndm)
czr          write(iow,*) i,(fx(i,inx),i=1,ndm)
czr          write(iow,*) nx,i,f(nx,i)
czr          write(*,*) nx,i,f(nx,i)
      endif
c.... second element
      do 400 jj=ii+1,ikz(i+1)-1
        iel=(ike(jj))
        do 500 iii=1,nen
          if(iek(iii,iel).eq.i)then
            ii2=iii-1
            if(ii2.eq.0) ii2=4
            i2=iek(ii2,iel)
czr         call curve(x(1,i2),y0,0.d0,m0,nrt1(iif,1),0.d0)
            call curve(x(1,i2),y0, x(2,i2) ,m0,nrt1(iif,1),0.d0)
            if(abs(x(2,i2)-y0).lt.0.001)then
              if(i2.eq.i1)goto 400
              ist=2
              goto 600
            else
              ii2=iii+1
              if(ii2.eq.5) ii2=1
              i2=iek(ii2,iel)
czr           call curve(x(1,i2),y0,0.d0,m0,nrt1(iif,1),0.d0)
              call curve(x(1,i2),y0, x(2,i2) ,m0,nrt1(iif,1),0.d0)
              if(abs(x(2,i2)-y0).lt.0.001)then
                if(i2.eq.i1)goto 400
                ist=2
                goto 600
              endif
            endif
          endif
500     continue
400   continue
      goto 998    ! return
600   continue
      xl1=dsqrt((x(1,i2)-x(1,i))**2+(x(2,i2)-x(2,i))**2)
      xl11=dsqrt((x(1,i2)-cpar(iif,1))**2+(x(2,i2)-cpar(iif,2))**2)
      xl12=dsqrt((x(1,i2)-cpar(iif,3))**2+(x(2,i2)-cpar(iif,4))**2)
      if((xl11.gt.xl0).or.(xl12.gt.xl0).or.
     1           (xl11.gt.xl0.and.xl12.gt.xl0))then
        xl=xl11
        if(xl11.gt.xl12)xl=xl12
        xl1=xl1-xl
        px2=cpar(iif,7)
        if(xl11.gt.xl12)px2=cpar(iif,8)
        f(nx,i)=f(nx,i)+xl1*(px1+px2)/2.d0
      else
        px2=(cpar(iif,7)*xl12+cpar(iif,8)*xl11)/(xl11+xl12)
        f(nx,i)=f(nx,i)+xl1*(px1/3.d0+px2/6.d0)
      endif
998    return
czr999   format(/,A)
czr1001  format ('FHG   Knoten    Kraft ',/, 2(I5, 2x) , F10.2)
czr1000  format(I4, ',0,' , 3(F7.2,','))
      end
c...............................................................end.SR.pload1
c...................................................................SR.pload2
      subroutine pload2
     1      (iek,ike,ikz,f,x,
     2      ndm,ndf,nen,i,iif,px1,xl0)
c ---------------------------------------------------------------------------
c.... Streckenlasten auf Knoten fuer 3D-Modelle
c ---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension iek(7,*),ike(*),ikz(*),f(ndf,*),x(ndm,*)
      dimension m0(2)
czr   common /curvedat/  cpar(20,8),nrt1(20,3),ic,nn3
      common /curvedat/  cpar(20,8),nrt1(20,3),ic,nbe,nn3
      ist=0
c.... first element
      do 100 ii=ikz(i),ikz(i+1)-1
        iel=(ike(ii))
        do 200 iii=1,nen
          if(iek(iii,iel).eq.i)then
            ii1=iii-1
            if(ii1.eq.0) ii1=4
            i1=iek(ii1,iel)
CCCC        i1=abs(iek(ii1,iel))
cba         call curves(x(1,i1),x(2,i1),x(3,i1),m0,nrt1(iif,1),diff)
            call curve(x(1,i1),x(2,i1),x(3,i1),m0,nrt1(iif,1),diff)
            if(abs(diff).lt.0.001)then
              ist=1
              goto 300
            else
              ii1=iii+1
              if(ii1.eq.5) ii1=1
              i1=iek(ii1,iel)
CCCCC         i1=abs(iek(ii1,iel))
cba           call curves(x(1,i1),x(2,i1),x(3,i1),m0,nrt1(iif,1),diff)
              call curve(x(1,i1),x(2,i1),x(3,i1),m0,nrt1(iif,1),diff)
              if(abs(diff).lt.0.001)then
                ist=1
                goto 300
              endif
            endif
          endif
200     continue
100   continue
300   continue
crene  if(ist.ne.1)then stop
      if(ist.ne.1)then
        write(*,*)'Achtung Elemente ausserhalb'
        goto 998         !return
      endif
      xl1=dsqrt((x(1,i1)-x(1,i))**2+(x(2,i1)-x(2,i))**2+
     1           (x(3,i1)-x(3,i))**2)
c.... Gerade mit 2 Punkten
      if(nrt1(iif,2).eq.32)then
c...
        xl11=dsqrt((x(1,i1)-cpar(iif,1))**2+(x(2,i1)-cpar(iif,2))**2+
     1             (x(3,i1)-cpar(iif,3))**2)
        xl12=dsqrt((x(1,i1)-cpar(iif,4))**2+(x(2,i1)-cpar(iif,5))**2+
     1             (x(3,i1)-cpar(iif,6))**2)
        nx=nrt1(iif,3)-30
        if((xl11.gt.xl0).or.(xl12.gt.xl0).or.
     1           (xl11.gt.xl0.and.xl12.gt.xl0))then
          xl=xl11
          if(xl11.gt.xl12)xl=xl12
          xl1=xl1-xl
          px2=cpar(iif,7)
          if(xl11.gt.xl12)px2=cpar(iif,8)
          f(nx,i)=f(nx,i)+xl1*(px1+px2)/2.d0
        else
          px2=(cpar(iif,7)*xl12+cpar(iif,8)*xl11)/(xl11+xl12)
          f(nx,i)=f(nx,i)+xl1*(px1/3.d0+px2/6.d0)
        endif
c.... Kreis
      elseif(nrt1(iif,2).eq.42)then
        nx=nrt1(iif,3)-30
        px=cpar(iif,7)
        f(nx,i)=f(nx,i)+xl1*px/2.d0
      else
        write(*,*)'remesh3 line433'
        stop
      endif
c.... second element
      do 400 jj=ii+1,ikz(i+1)-1
        iel=(ike(jj))
        do 500 iii=1,nen
          if(iek(iii,iel).eq.i)then
            ii2=iii-1
            if(ii2.eq.0) ii2=4
            i2=iek(ii2,iel)
cba          call curves(x(1,i2),x(2,i2),x(3,i2),m0,nrt1(iif,1),diff)
            call curve(x(1,i2),x(2,i2),x(3,i2),m0,nrt1(iif,1),diff)
            if(abs(diff).lt.0.001)then
              if(i2.eq.i1)goto 400
              ist=2
              goto 600
            else
              ii2=iii+1
              if(ii2.eq.5) ii2=1
              i2=iek(ii2,iel)
cba           call curves(x(1,i2),x(2,i2),x(3,i2),m0,nrt1(iif,1),diff)
              call curve(x(1,i2),x(2,i2),x(3,i2),m0,nrt1(iif,1),diff)
              if(abs(diff).lt.0.001)then
                if(i2.eq.i1)goto 400
                ist=2
                goto 600
              endif
            endif
          endif
500     continue
400   continue
      goto 998   ! return
600   continue
      xl1=dsqrt((x(1,i2)-x(1,i))**2+(x(2,i2)-x(2,i))**2+
     1           (x(3,i2)-x(3,i))**2)
c.... Gerade mit 2 Punkten
      if(nrt1(iif,2).eq.32)then
          xl11=dsqrt((x(1,i2)-cpar(iif,1))**2+(x(2,i2)-cpar(iif,2))**2+
     1             (x(3,i2)-cpar(iif,3))**2)
          xl12=dsqrt((x(1,i2)-cpar(iif,4))**2+(x(2,i2)-cpar(iif,5))**2+
     1             (x(3,i2)-cpar(iif,6))**2)
        if((xl11.gt.xl0).or.(xl12.gt.xl0).or.
     1           (xl11.gt.xl0.and.xl12.gt.xl0))then
          xl=xl11
          if(xl11.gt.xl12)xl=xl12
          xl1=xl1-xl
          px2=cpar(iif,7)
          if(xl11.gt.xl12)px2=cpar(iif,8)
          f(nx,i)=f(nx,i)+xl1*(px1+px2)/2.d0
        else
          px2=(cpar(iif,7)*xl12+cpar(iif,8)*xl11)/(xl11+xl12)
          f(nx,i)=f(nx,i)+xl1*(px1/3.d0+px2/6.d0)
        endif
      elseif(nrt1(iif,2).eq.42)then
        nx=nrt1(iif,3)-30
        px=cpar(iif,7)
        f(nx,i)=f(nx,i)+xl1*px/2.d0
      else
        write(*,*)'remesh3 line492'
        stop
      endif
998   continue
      return
      end
c...............................................................end.SR.pload2
c.................................................................SR.conarray
      subroutine conarray(a,b,c,la,lb)
c ---------------------------------------------------------------------------
c.... attach b to the end of a on c
c.... la - length of a
c.... lb - length of b
c ---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension a(*),b(*),c(*)
      do i = 1,la+lb
        if (i.le.la) then
          c(i) = a(i)
        else
          c(i) = b(i-la)
        endif
      enddo
      end
c.............................................................end.SR.conarray
c....................................................................SR.wcurv
      subroutine wcurv(inpf)
c----------------------------------------------------------------------------
c.... writes curves to input file
c----------------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      common /curvedat/  cpar(20,8),nrt1(20,3),ic,nbe,nn3
c....
      write(inpf,999) 'curv'
c.... write old   curv  into input file
      do ni=1,ic-1
        write(inpf,1001)  (nrt1(ni,i),i=1,3)
        write(inpf,1002)  (cpar(ni,i),i=1,8)
      enddo
999   format(/, A )
1001  format(3(I4,','))
1002  format(8(F8.4,','))
      end
c................................................................end.SR.wcurv
c....................................................................SR.wmate
      subroutine wmate(inpf)
c----------------------------------------------------------------------------
c.... read old input file / write new input file
c----------------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      character*80 yyy
c....
      nior = abs(ior)
      rewind nior
c.... find  mate
      do i=1,10000000     ! do forever
        read(nior,1000,end=100) yyy
        if((yyy(1:4).eq.'mate').or.(yyy(1:4).eq.'MATE')) then
          write(inpf,999) yyy
          goto 200
        endif
      enddo
100   continue
        write(*,'(A)') '*** no mate found in input file'
200   continue
c.... read/write  mate
      do i=1,1000
        read(nior,1000,err=300,end=300) yyy
czr        if(yyy(1:10).eq.'          ') goto 300
        if (((ichar(yyy(1:1)).ge.65 ).and.(ichar(yyy(1:1)).lt.77 )) .or.
     1      ((ichar(yyy(1:1)).gt.77 ).and.(ichar(yyy(1:1)).le.90 )) .or.
     1      ((ichar(yyy(1:1)).ge.97 ).and.(ichar(yyy(1:1)).lt.109)) .or.
     1      ((ichar(yyy(1:1)).gt.109).and.(ichar(yyy(1:1)).le.122)))
     2      goto 300
        write(inpf,998) yyy
      enddo
300   continue
      write(inpf,998)'end '
30    continue    ! do forever
        read(nior,1000,end=130) yyy
        if((yyy(1:3).eq.'opt').or.(yyy(1:3).eq.'OPT')) then
          write(inpf,999) yyy
          goto 130
        endif
      goto 30
130   continue
      write(inpf,999)'inte'
      write(inpf,999)'stop'
      write(inpf,999)
998   format( A )
999   format(/, A )
1000  format(A80)
      end
c................................................................end.SR.wmate
c...................................................................SR.curve5
      subroutine curve5 (xx,nr,nxst)
c----------------------------------------------------------------------------
c----------------------------------------------------------------------------
      USE cdata
      implicit double precision (a-h,o-z)
      common /curvedat/  cpar(20,8),nrt1(20,3),ic,nbe,nn3
      dimension xx(6,5)
czr      dimension x0(3),a(3),b(3),c(3),d(4)
      dimension x0(3),a(3),b(3),     d(4)
czr      dimension xa(3),xaa(3),xm(3)
      dimension           xm(3)
      dimension f1(4),f2(4),f3(4),fx(3)
      tol=0.0001
c.... first middle plane function f2(1)*x+f2(2)*y+f2(3)*z+f2(4)=0
      f1(1)=2.d0*(xx(1,1)-xx(1,3))
      f1(2)=2.d0*(xx(2,1)-xx(2,3))
      f1(3)=2.d0*(xx(3,1)-xx(3,3))
      f1(4)=xx(1,3)**2+xx(2,3)**2+xx(3,3)**2
     1        -xx(1,1)**2-xx(2,1)**2-xx(3,1)**2
c.... second middle plane function f2(1)*x+f2(2)*y+f2(3)*z+f2(4)=0
      f2(1)=2.d0*(xx(1,2)-xx(1,4))
      f2(2)=2.d0*(xx(2,2)-xx(2,4))
      f2(3)=2.d0*(xx(3,2)-xx(3,4))
      f2(4)=xx(1,4)**2+xx(2,4)**2+xx(3,4)**2
     1        -xx(1,2)**2-xx(2,2)**2-xx(3,2)**2
c.... third middle plane function
      x0(1)=xx(1,2)
      x0(2)=xx(2,2)
      x0(3)=xx(3,2)
      a(1)=xx(1,3)-xx(1,2)
      a(2)=xx(2,3)-xx(2,2)
      a(3)=xx(3,3)-xx(3,2)
      b(1)=xx(1,1)-xx(1,2)
      b(2)=xx(2,1)-xx(2,2)
      b(3)=xx(3,1)-xx(3,2)
      call vcross(b,a,f3)
      call vnorm (f3,ylb)
      f3(4)=(f3(1)*x0(1)+f3(2)*x0(2)+f3(3)*x0(3))*(-1.d0)
*        write(*,*) f1(1),f1(2),f1(3),f1(4)
*        write(*,*) f2(1),f2(2),f2(3),f2(4)
*        write(*,*) f3(1),f3(2),f3(3),f3(4)
c.... calculate the line function of two planes x=d+lambda*c
      call vcross(f1,f2,b)
      call vnorm(b,ylb)
c.... calculate start point
      fx(1)=f1(4)*(-1.d0)
      fx(2)=f2(4)*(-1.d0)
      fx(3)=f3(4)*(-1.d0)
      call  i3plane(f1,f2,f3,fx,d(1),d(2),d(3))
*        write(*,*) b(1),b(2),b(3)
*        write(*,*) d(1),d(2),d(3)
      if(nxst.eq.1)then
c.... boundary geometry
        xx(1,5)=d(1)
        xx(2,5)=d(2)
        xx(3,5)=d(3)
      elseif(nxst.eq.2)then
c.... surface geometry
        call geopoi(b,d,tol,xm,nr)
        xx(1,5)=xm(1)
        xx(2,5)=xm(2)
        xx(3,5)=xm(3)
      endif
      return
      end
c...............................................................end.SR.curve5
c..................................................................SR.i3plane
      subroutine i3plane(f1,f2,f3,fx,x,y,z)
c----------------------------------------------------------------------------
c----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension f1(4),f2(4),f3(4),fx(3),x1(3,3),x1i(3,3)
      do i=1,3
        x1(1,i)=f1(i)
        x1(2,i)=f2(i)
        x1(3,i)=f3(i)
      enddo
      det=x1(1,1)*(x1(2,2)*x1(3,3)-x1(2,3)*x1(3,2))+
     1    x1(1,2)*(x1(2,3)*x1(3,1)-x1(2,1)*x1(3,3))+
     2    x1(1,3)*(x1(2,1)*x1(3,2)-x1(2,2)*x1(3,1))
czr      if(det.eq.0.d0) stop
      if(abs(det).lt.0.0001)
     +   write(*,*)'remesh3 line 671'
      if(abs(det).lt.0.0001) stop
      do 200 i=1,3
        i1=mod(i,3)+1
        i2=mod(i1,3)+1
        do 300 j=1,3
          j1=mod(j,3)+1
          j2=mod(j1,3)+1
          x1i(j,i)=(x1(i1,j1)*x1(i2,j2)-x1(i1,j2)*x1(i2,j1))/det
300     continue
200   continue
      x=x1i(1,1)*fx(1)+x1i(1,2)*fx(2)+x1i(1,3)*fx(3)
      y=x1i(2,1)*fx(1)+x1i(2,2)*fx(2)+x1i(2,3)*fx(3)
      z=x1i(3,1)*fx(1)+x1i(3,2)*fx(2)+x1i(3,3)*fx(3)
      return
      end
c..............................................................end.SR.i3plane
c----------------------------------------------------------------------------
c----------------------------------------------------------------------------

