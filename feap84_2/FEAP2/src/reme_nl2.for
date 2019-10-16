c .................................................................SR xgen_nl
      subroutine xgen_nl(nst,i,u0,u,hst0,x0,x,iek,ikz,ike, 
     +      mikno,nen1,nhist,ndm,ndf,
     +      ianp,ianp0,iael,iael0,naiter,erro, erron0)
c ---------------------------------------------------------------------------
c.... subdivide element with a cross
c.... interpolate displacements to new node
c.... nonlinear extension   interpolate displacements
c.... plastic extension   interpolate history data
c ---------------------------------------------------------------------------
      USE cdata
      implicit double precision (a-h,o-z)
czr   logical lav,lavx
czr   common /gene0/ kk0,ke0
czr   common /gefehl/ nfehler
czr   common /adap3/  n9e,asteu,fsteu,fmax,nerro
czr   dimension m0(2),xxh(4,9),ih(4),vn(3)
      logical lav
      common /shel1/  ndir,nlav,lav
      common /gener/ kk,km,ke
      common /curvedat/  cpar(20,8),nrt1(20,3),ic,nbe,nn3
      dimension x(ndm,*),x0(ndm,*),xx(6,5),
     +       u(ndf,*),u0(ndf,*),ux(6,5),           !.... nonlinear
     +       hst0(nhist,*),                        !.... plasti
c    +       hst0(nhist,*),stx(20,5),    
     +       iek(nen1,*),ikz(*),ike(*),
     +       nk(5),mikno(3,*),
     +       ianp(*),ianp0(*),iael0(*),iael(*),
     +       erro(*),erron0(*),
     +       xxh(4,9),ih(4)
czr  +       vn(3)     !   Direktor Ermittlung
      integer*2 l1(5),l2(5)
      data l1/-1,0,1,0,0/,l2/0,1,0,-1,0/
      necke=0
      nfunc1=0
      nfunc15=0
c.... curves on shell surface
      nshell=0
      call geos(iek,x,x0,i,ndm,nshell,nen1,xxh,ih,iif,1)
      if(nshell.eq.1) nfunc15=1
      if(nshell.eq.1) nfunc1=1
      do 200 j=0,4
      if(j.eq.4)goto 250
      ianp0(abs(iek(j+1,i)))=abs(ianp0(abs(iek(j+1,i))))
      ik1=abs(iek(j+1,i))
      if(j.eq.3)then
        ik2=abs(iek(1,i))
      else
        ik2=abs(iek(j+2,i))
      end if
      call kaelem(ik1,ik2,i2,i,ikz,ike)
      if(iek(j+1,i).gt.0)  goto 210
c.... bereits generierten kantenknoten ermitteln
      call miknof(ik1,ik2,km,kmi,mikno)
      nk(j+1)=kmi
      do jj=1,ndm
        xx(jj,j+1)= x0(jj,kmi-numnp)
      enddo
cnl>
      do jj=1,ndm
        ux(jj,j+1)= u0(jj,kmi-numnp)
      enddo
cnl<
      if(ndm.eq.6)then
        xx(1,j+1)=(xx(1,j+1)+xx(4,j+1))*0.5d0
        xx(2,j+1)=(xx(2,j+1)+xx(5,j+1))*0.5d0
        xx(3,j+1)=(xx(3,j+1)+xx(6,j+1))*0.5d0
      endif
c.... ike,ikz modifizieren
      call ikemod (iek,ike,ikz,kmi,i,0,ke,kk,j,1,nen1)
      goto 200
210   if(i2.eq.0)goto250
      do jj=1,4
        if(iek(jj,i2).eq.ik1.or.iek(jj,i2).eq.ik2) goto 230
      enddo
      nfehler=1
      write(*,*)'netzmo nicht kompatibel'
      goto 998                           ! return
c.... Nachbarelement negativ setzen
230   if( jj.eq.1) then
        if(abs(iek(2,i2)).eq.ik1.or.abs(iek(2,i2)).eq.ik2)then
          iek(jj,i2)=iek(jj,i2)*(-1)
        else
          iek(4,i2)=iek(4,i2)*(-1)
        endif
      else
        iek(jj,i2)=iek(jj,i2)*(-1)
      endif
c.... kantenmittelknoten abspeichern
      km=km+1
      mikno(1,km) =ik1
      mikno(2,km) =ik2
      mikno(3,km) =kk+1
c.... ike,ikz erweitern
250   call ikemod (iek,ike,ikz,kmi,i,i2,ke,kk,j,2,nen1)
c.... knotennummern 
      kk=kk+1
      if(j.lt.4)then
        ianp0(kk)=naiter*(-1)
      else
        ianp0(kk)=naiter
      endif
      nk(j+1)=kk
czr       if(nst.eq.0.and.fsteu.eq.1)then
czrc..     calculate nodal error for new nodes(adaf)
czr        call xerr(erron,erron0,abs(iek(1,i)),ecn1)
czr        call xerr(erron,erron0,abs(iek(2,i)),ecn2)
czr        call xerr(erron,erron0,abs(iek(3,i)),ecn3)
czr        call xerr(erron,erron0,abs(iek(4,i)),ecn4)
czr       erron0(kk-numnp)=1/4.*((1-l2(j+1))*(1-l1(j+1))*ecn1+
czr     +                        (1+l2(j+1))*(1-l1(j+1))*ecn2+
czr     +                        (1+l2(j+1))*(1+l1(j+1))*ecn3+
czr     +                        (1-l2(j+1))*(1+l1(j+1))*ecn4)
czr       endif
c.... curve interpolation
c.... nodal coordinates and displacements
      ncuv=0
      ncuvr=0
      if(j.eq.4) goto 800
      if(nshell.eq.1) ncuv=1
      jj=j+2
      if(jj.eq.5) jj=1
      ia =j+1
      ib =jj
      i1=abs(iek(j+1,i))
      i2=abs(iek(j+1,i))
      j1=abs(iek(jj,i))
      j2=abs(iek(jj,i))
      call geosr                   !... detect boundary geometry function 
     1     (ia,ib,iek,x,x0,i,ndm,ncuvr,nen1,xxh,ih,iifr,iif,nshell,1)
800   continue
      if(ncuv+ncuvr.ne.0)then
        if (ncuvr.ne.0) necke=necke+1
c....
        if(j+1.eq.1)then
          ia=1
          ib=2
        elseif(j+1.eq.2)then
          ia=2
          ib=3
        elseif(j+1.eq.3)then
          ia=3
          ib=4
        else
          ia=4
          ib=1
        endif
        nfunc1=1
        if(ncuv.eq.0.and.ncuvr.eq.1)then
c.... boundary geometry function
          if(ndm.eq.2)then                              !..... czr>
            xx(3,j+1)=0.0d0
cnl            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1),
cnl     1            nrt1(iifr,1),nrt1(iifr,1),      3   ,1)
            call curven_nl(ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1),
     +            u(1,ik1),u(1,ik2),ux(1,j+1),ndf,
     1            nrt1(iifr,1),nrt1(iifr,1),      3   ,1)
          else
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1           ,nrt1(iifr,1),nrt1(iifr,1),ndm,1)
          endif                                         !..... czr<
        elseif(ncuv.eq.1.and.ncuvr.eq.0)then
c.... only surface geometry function
          call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1        ,nrt1(iif,1),nrt1(iif,1),ndm,2)
        elseif(ncuv.eq.1.and.ncuvr.eq.2)then
c.... surface geometry function and plane boundary function
          call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1        ,nrt1(iif,1),nrt1(iifr,1),ndm,3)
        elseif(ncuv.eq.1.and.ncuvr.eq.3)then
c.... surface geometry function and boundary function
          call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1        ,nrt1(iif,1),nrt1(iifr,1),ndm,4)
        endif
c....
        if(ndm.eq.2)then              !.....  czr>
          do k=1,ndm
            x0(k,kk-numnp)=xx(k,j+1)
          enddo                       !.....  czr<
        elseif(ndm.eq.3) then
          do k=1,ndm
            x0(k,kk-numnp)=xx(k,j+1)
          enddo
        else if(ndm.eq.6)then
czr           call dirsea(xxh,ia,ib,vn,1,tn)
czr          do k=1,3
czr***        do 275 k=1,ndm
czr**         write(*,*)'keine Direktorermittrung warnung'
czr            x0(k+3,kk-numnp)=xx(k,j+1)-tn*vn(k)*0.5d0
czr            x0(k,kk-numnp)=xx(k,j+1)+tn*vn(k)*0.5d0
czrczr275      continue   
czr          enddo
        endif

cnl>
c.... evaluate nodal displacements        u_0 = Sum (N_i*–) 
        do k=1,ndf
c          u0(k,kk-numnp)=ux(k,j+1)
          u0(k,kk-numnp)=ux(k,j+1)
        enddo
cpl>
c.... evaluate nodal history data      st_0 = Sum (st_i*–) 
c        do k=1,nhist
cc....     node=abs(iek(1,i)
c          stcn1 = hst0(k , abs(iek(1,i)) )          
c          stcn2 = hst0(k , abs(iek(2,i)) )          
c          stcn3 = hst0(k , abs(iek(3,i)) )          
c          stcn4 = hst0(k , abs(iek(4,i)) )          
cc         stx(k,j+1)=1/4.*((1-l2(j+1))*(1-l1(j+1))*stcn1+
cc    +                     (1+l2(j+1))*(1-l1(j+1))*stcn2+
cc    +                     (1+l2(j+1))*(1+l1(j+1))*stcn3+
cc    +                     (1-l2(j+1))*(1+l1(j+1))*stcn4)
cc          hst0(k,kk)=stx(k,j+1)
c          hst0(k,kk)=1/4.*((1-l2(j+1))*(1-l1(j+1))*stcn1+
c     +                     (1+l2(j+1))*(1-l1(j+1))*stcn2+
c     +                     (1+l2(j+1))*(1+l1(j+1))*stcn3+
c     +                     (1-l2(j+1))*(1+l1(j+1))*stcn4)
c        enddo
cnl<
      else
c.... ncuv+ncuvr=0
        do k=1,ndm
          call xcor(x,x0,abs(iek(1,i)),k,xcn1,ndm)
          call xcor(x,x0,abs(iek(2,i)),k,xcn2,ndm)
          call xcor(x,x0,abs(iek(3,i)),k,xcn3,ndm)
          call xcor(x,x0,abs(iek(4,i)),k,xcn4,ndm)
          xx(k,j+1)=1/4.*((1-l2(j+1))*(1-l1(j+1))*xcn1+
     +                   (1+l2(j+1))*(1-l1(j+1))*xcn2+
     +                   (1+l2(j+1))*(1+l1(j+1))*xcn3+
     +                   (1-l2(j+1))*(1+l1(j+1))*xcn4)
          x0(k,kk-numnp)=xx(k,j+1)
        enddo
cnl>
c.... evaluate nodal displacements        u_0 = Sum (N_i*–) 
        do k=1,ndf
          call xcor(u,u0,abs(iek(1,i)),k,ucn1,ndf)
          call xcor(u,u0,abs(iek(2,i)),k,ucn2,ndf)
          call xcor(u,u0,abs(iek(3,i)),k,ucn3,ndf)
          call xcor(u,u0,abs(iek(4,i)),k,ucn4,ndf)
          ux(k,j+1)=1/4.*((1-l2(j+1))*(1-l1(j+1))*ucn1+
     +                    (1+l2(j+1))*(1-l1(j+1))*ucn2+
     +                    (1+l2(j+1))*(1+l1(j+1))*ucn3+
     +                    (1-l2(j+1))*(1+l1(j+1))*ucn4)
          u0(k,kk-numnp)=ux(k,j+1)
        enddo
cpl>
c.... evaluate nodal history data      st_0 = Sum (st_i*–) 
        do k=1,nhist
c....     node=abs(iek(1,i)
          stcn1 = hst0(k , abs(iek(1,i)) )          
          stcn2 = hst0(k , abs(iek(2,i)) )          
          stcn3 = hst0(k , abs(iek(3,i)) )          
          stcn4 = hst0(k , abs(iek(4,i)) )          
c          stx(k,j+1)=1/4.*((1-l2(j+1))*(1-l1(j+1))*stcn1+
c     +                     (1+l2(j+1))*(1-l1(j+1))*stcn2+
c     +                     (1+l2(j+1))*(1+l1(j+1))*stcn3+
c    +                     (1-l2(j+1))*(1+l1(j+1))*stcn4)
c          hst0(k,kk)=stx(k,j+1)
          hst0(k,kk)=1/4.*((1-l2(j+1))*(1-l1(j+1))*stcn1+
     +                     (1+l2(j+1))*(1-l1(j+1))*stcn2+
     +                     (1+l2(j+1))*(1+l1(j+1))*stcn3+
     +                     (1-l2(j+1))*(1+l1(j+1))*stcn4)
        enddo
cnl<
        if(nfunc15.eq.0.and.ndm.eq.6)then
          xx(1,j+1)=(xx(1,j+1)+xx(4,j+1))*0.5d0
          xx(2,j+1)=(xx(2,j+1)+xx(5,j+1))*0.5d0
          xx(3,j+1)=(xx(3,j+1)+xx(6,j+1))*0.5d0
        endif
        if(j+1.eq.5.and.nfunc1.eq.1)then
          if(nfunc15.eq.0)then
c.... only boundary geometry function
            if(ndm.eq.2)then                 !.....  czr>
              xx(3,j+1)=0.0d0
              do izx=1,5
                xx(3,izx)=0.0d0
              enddo
              call curve5 (xx,nrt1(iifr,1),1)
            else
              call curve5 (xx,nrt1(iifr,1),1)
            endif                            !.....  czr<
          elseif(nfunc15.eq.1)then
c.... surface geometry funtion
            call curve5 (xx,nrt1(iif,1),2)
          endif
          if(ndm.eq.2)then
            do k=1,ndm
              x0(k,kk-numnp)=xx(k,j+1)
            enddo
          elseif(ndm.eq.3) then
            do k=1,ndm
              x0(k,kk-numnp)=xx(k,j+1)
            enddo
          else if(ndm.eq.6)then
             write(*,*)'ndm.eq.6'
czr            call dirsea(xxh,ia,ib,vn,2,tn)
czr            do k=1,3
czr***           write(*,*)'keine Direktorermittlung Warnung'
czr              x0(k+3,kk-numnp)=xx(k,j+1)-tn*vn(k)*0.5d0
czr              x0(k,kk-numnp)=xx(k,j+1)+tn*vn(k)*0.5d0
czr            enddo
          endif
        endif
      endif
200   continue
c.... nachbarelemente markieren
      if(nst.eq.1) goto 700
      do 500j=1,4
        nn1=abs(ikz(abs(iek(j,i))))
        nn2=abs(ikz(abs(iek(j,i))+1))
        do 500 jj=nn1,nn2-1
          if(iek(nen1,ike(jj)).lt.0)goto 500
czr      if (fsteu.ne.1)then
c.... 'normal in BS' fsteu=0.0d0 in pcontr
          if(iek(4,ike(jj)).gt.numnp)goto 500
czr      elseif (fsteu.eq.1)then
czr       if(iek(4,ike(jj)).gt.kk0)goto 500
czr      endif
          iek(nen1,ike(jj))=iek(nen1,ike(jj))*(-1)
500   continue
C
c.... element knotenbeziehung neu schreiben
      iael0(i)=abs(iael(i))+1
700   k0=0
      do 400 k=2,4
        k0=k0+1
        ke=ke+1
        iael0(ke)=abs(iael0(i))
        iek(nen1,ke)=abs(iek(nen1,i))
        if(necke.eq.4) iek(nen1,ke)=iek(nen1,ke)*(-1.d0)
        iek(nen1-1,ke)=abs(iek(nen1-1,i))
        iek(nen1-2,ke)=abs(iek(nen1-2,i))
        iek(1,ke)=nk(k0)
        iek(2,ke)=abs(iek(k,i))
        iek(3,ke)=nk(k0+1)
        iek(4,ke)=nk(5)
czr       if(nst.eq.0.and.fsteu.eq.1)then
czrc.... calculate element error for new elements (adaf)
czr        call xerr(erron,erron0,abs(iek(1,ke)),ecn1)
czr        call xerr(erron,erron0,abs(iek(2,ke)),ecn2)
czr        call xerr(erron,erron0,abs(iek(3,ke)),ecn3)
czr        call xerr(erron,erron0,abs(iek(4,ke)),ecn4)
czr       erro(ke)=(ecn1+ecn2+ecn3+ecn4)/4.d0
czr       if(erro(ke).gt.fmax) fmax=erro(ke)
czr       elseif(fsteu.eq.1)then
czr       erro(ke)=-1.d0
czr       else
czr       endif
C
400   continue
c.... ike,ikz modifizieren (ecken)
      call ikemod (iek,ike,ikz,kmi,i,i2,ke,kk,j,3,nen1)
      iek(nen1,i)=abs(iek(nen1,i))
      if(necke.eq.4) iek(nen1,i)=iek(nen1,i)*(-1.d0)
      iek(2,i)=abs(iek(1,i))
      iek(3,i)=nk(1)
      iek(1,i)=nk(4)
      iek(4,i)=nk(5)
czr       if(nst.eq.0.and.fsteu.eq.1)then
czrc.... calculate element error for new elements (adaf)
czr        call xerr(erron,erron0,abs(iek(1,i)),ecn1)
czr        call xerr(erron,erron0,abs(iek(2,i)),ecn2)
czr        call xerr(erron,erron0,abs(iek(3,i)),ecn3)
czr        call xerr(erron,erron0,abs(iek(4,i)),ecn4)
czr       erro(i)=(ecn1+ecn2+ecn3+ecn4)/4.d0
czr       if(erro(i).gt.fmax) fmax=erro(i)
czr       else
czr       erro(i)=-1.d0
czr       endif
czrCCC   call exnode(nk,i)
998   return
      end
c ...............................................................end xgen_nl
      subroutine ygen_nl
     1      (nst,i,u0,u,hst0,x0,x,iek,ikz,ike,
     +      mikno,nen1,nhist,ndm,ndf,
     2      ianp,ianp0,iael,iael0,naiter,erro)
c ---------------------------------------------------------------------------
c.... ygen_nl  unterteilt das element mit einem y
c.... project displacements to new nodes
c ---------------------------------------------------------------------------
      USE cdata
      implicit double precision (a-h,o-z)
czr   logical lav,lavx
czr   dimension nk0(4),nk(5)
czr   dimension xxh(4,9),ih(4),xx(6,5),vn(3)
      dimension x(ndm,*),x0(ndm,*),xx(6,5),xxh(4,9),
     +      u(ndf,*),u0(ndf,*),ux(6,5),           !.... nonlinear
czr  +      hst0(nhist,*),stx(20,5),           !.... plasti
     +      hst0(nhist,*),           !.... plasti
     +      iek(nen1,*),ikz(*),ike(*),
     +      ianp(*),ianp0(*),iael0(*),iael(*),
     +      nk0(4),mikno(3,*),ih(4),
     +      erro(*)
czr     +      vn(3)
czr   common /gefehl/ nfehler
czr   common /adap3/  n9e,asteu,fsteu,fmax,nerro
      logical lav
      common /shel1/  ndir,nlav,lav
      common /gener/ kk,km,ke
      common /curvedat/ cpar(20,8),nrt(20,3),ic,nbe,nn3
c.... lokale Elementnumerierung
      if(naiter.eq.1) then
        k=2
      else
        do k=1,4
          if(ianp0(abs(iek(k,i))).lt.0 )goto 110
        enddo
      endif
110   nk0(1)=iek(k,i)
      do j=2,4
        k=k+1
        if(k.eq.5) k=1
        nk0(j)=iek(k,i)
      enddo
c.... curves on shell surface
      nfunc1=0
      nfunc15=0
      nshell=0
      do j=1,4
        ih(j)=abs(nk0(j))
      enddo
      call geos(iek,x,x0,i,ndm,nshell,nen1,xxh,ih,iif,2)
      if(nshell.eq.1) nfunc15=1
      if(nshell.eq.1) nfunc1=1
c.... curve 0.kante (keine Verfeinerung)
      ia=1
      ib=2
      j=0
90    ncuv=0
      ncuvr=0
      if(nshell.eq.1) ncuv=1
      call geosr                                            ! in SR ygenes
     1     (ia,ib,iek,x,x0,i,ndm,ncuvr,nen1,xxh,ih,iifr,iif,nshell,2)
      if(ncuv+ncuvr.ne.0)then
        nfunc1=1
        if(ncuv.eq.0.and.ncuvr.eq.1)then
c.... boundary geometry funtion
czr.. ndm.eq.2 added >
          if(ndm.eq.2)then
            xx(3,j+1)=0.0d0
            lav = .true.    ! undef
            nlav=2          ! undef
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1),
     1            nrt(iifr,1),nrt(iifr,1),      3   ,1)
          else
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1           ,nrt(iifr,1),nrt(iifr,1),ndm,1)
          endif
        elseif(ncuv.eq.1.and.ncuvr.eq.0)then
c.... only surface geometry function
          call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1        ,nrt(iif,1),nrt(iif,1),ndm,2)
        elseif(ncuv.eq.1.and.ncuvr.eq.2)then
c.... surface geometry funtion and plane boundary function
          call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1                 ,nrt(iif,1),nrt(iifr,1),ndm,3)
        elseif(ncuv.eq.1.and.ncuvr.eq.3)then
c.... surface geometry function and boundary function
          call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1                 ,nrt(iif,1),nrt(iifr,1),ndm,4)
        endif
      else
c.... ncuv+ncuvr=0
        do k=1,ndm
          call xcor(x,x0,ih(ia),k,xcn1,ndm)
          call xcor(x,x0,ih(ib),k,xcn2,ndm)
          xx(k,j+1)=0.5d0*(xcn1+xcn2)
        enddo 
cnl>
c.... evaluate nodal displacements        u_0 = Sum (N_i*–) 
        do k=1,ndf
          call xcor(u,u0,ih(ia),k,ucn1,ndf)
          call xcor(u,u0,ih(ib),k,ucn2,ndf)
          ux(k,j+1)=0.5d0*(ucn1+ucn2)
        enddo                              ! cnl<
        if(nfunc15.eq.0.and.ndm.eq.6)then
          xx(1,j+1)=(xx(1,j+1)+xx(4,j+1))*0.5d0
          xx(2,j+1)=(xx(2,j+1)+xx(5,j+1))*0.5d0
          xx(3,j+1)=(xx(3,j+1)+xx(6,j+1))*0.5d0
        endif
      endif
      if(j.eq.0)then
c.... 2. nicht zu verfeinernde Kante
        ia=4
        ib=1
        j=3
        goto 90
      endif
c.... 1.kante - verfeinerung
      ik1=abs(nk0(2))
      ik2=abs(nk0(3))
      call kaelem(ik1,ik2,i2,i,ikz,ike)
      if(nk0(2).lt.0)then
        call miknof(ik1,ik2,km,kmi,mikno)
        do jj=1,ndm
          xx(jj,2)= x0(jj,kmi-numnp)
        enddo
cnl..>
        do jj=1,ndf
          ux(jj,2)= u0(jj,kmi-numnp)
        enddo
        if(ndm.eq.6)then
          xx(1,2)=(xx(1,2)+xx(4,2))*0.5d0
          xx(2,2)=(xx(2,2)+xx(5,2))*0.5d0
          xx(3,2)=(xx(3,2)+xx(6,2))*0.5d0
        endif
c.... ike,ikz modifizieren
        call ikemod (iek,ike,ikz,kmi,i,0,ke,kk,0,7,nen1)
        km1=kmi
      else
c.... ike,ikz erweitern
        call ikemod (iek,ike,ikz,kmi,i,i2,ke,kk,1,4,nen1)
        if(i2.eq.0)goto 130
        do jj=1,4
          if(iek(jj,i2).eq.ik1.or.iek(jj,i2).eq.ik2)goto 230
        enddo
230     if( jj.eq.1) then
          if(abs(iek(2,i2)).eq.ik1.or.abs(iek(2,i2)).eq.ik2)then
            iek(jj,i2)=iek(jj,i2)*(-1)
          else
            iek(4,i2)=iek(4,i2)*(-1)
          endif
        else
          iek(jj,i2)=iek(jj,i2)*(-1)
        endif
c.... kantenmittelknoten abspeichern
        km=km+1
        mikno(1,km) =ik1
        mikno(2,km) =ik2
        mikno(3,km) =kk+1
c.... knotennummern und knotenkoordinaten
130     kk=kk+1
        ianp0(kk)=naiter*(-1)
        km1=kk
c.... curve interpolation
        ncuv=0
        i1=abs(ik1)
        i2=abs(ik1)
        j1=abs(ik2)
        j2=abs(ik2)
        ia=2
        ib=3
        j=1
        ncuv=0
        ncuvr=0
        if(nshell.eq.1) ncuv=1
        call geosr
     1     (ia,ib,iek,x,x0,i,ndm,ncuvr,nen1,xxh,ih,iifr,iif,nshell,2)
        if(ncuv+ncuvr.ne.0)then
          nfunc1=1
          if(ncuv.eq.0.and.ncuvr.eq.1)then
c.... boundary geometry funtion
            if(ndm.eq.2)then
              xx(3,j+1)=0.0d0
              lav = .true.    ! undef
              nlav=2          ! undef
              call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1),
     1            nrt(iifr,1),nrt(iifr,1),      3   ,1)
            else
              call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1           ,nrt(iifr,1),nrt(iifr,1),ndm,1)
            endif
          elseif(ncuv.eq.1.and.ncuvr.eq.0)then
c.... only surface geometry funtion
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1        ,nrt(iif,1),nrt(iif,1),ndm,2)
          elseif(ncuv.eq.1.and.ncuvr.eq.2)then
c.... surface geometry funtion and plane boundery function
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1        ,nrt(iif,1),nrt(iifr,1),ndm,3)
          elseif(ncuv.eq.1.and.ncuvr.eq.3)then
c.... surface geometry funtion and boundery function
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1        ,nrt(iif,1),nrt(iifr,1),ndm,4)
          endif
          if(ndm.eq.2)then
            do k=1,2
              x0(k,kk-numnp)=xx(k,j+1)
            enddo
          elseif(ndm.eq.3)then
            do k=1,3
               x0(k,kk-numnp)=xx(k,j+1)
            enddo
          elseif(ndm.eq.6)then
czr         call dirsea(xxh,ia,ib,vn,1,tn)
            write(*,*)'czr call dirsea(x '
            do k=1,3
cba           write(*,*)'keine Direktorermittrung Warnung'
czr              x0(k+3,kk-numnp)=xx(k,j+1)-tn*vn(k)*0.5d0
czr              x0(k,kk-numnp)=xx(k,j+1)+tn*vn(k)*0.5d0
            enddo
          endif
cnl>
czr       do k=1,ndf
czr         u0(k,kk-numnp)=ux(k,j+1)
czr       enddo
          do k=1,ndf
            call xcor(u,u0,ih(ia),k,ucn1,ndm)
            call xcor(u,u0,ih(ib),k,ucn2,ndm)
            ux(k,j+1)=0.5d0*(ucn1+ucn2)
            u0(k,kk-numnp)=(ucn1+ucn2)/2.d0
          enddo        ! cnl<
        else
c.... no geometry functions
          do k=1,ndm
            call xcor(x,x0,ih(ia),k,xcn1,ndm)
            call xcor(x,x0,ih(ib),k,xcn2,ndm)
            xx(k,j+1)=0.5d0*(xcn1+xcn2)
            x0(k,kk-numnp)=(xcn1+xcn2)/2.d0
          enddo
cnl>
          do k=1,ndf
            call xcor(u,u0,ih(ia),k,ucn1,ndm)
            call xcor(u,u0,ih(ib),k,ucn2,ndm)
            ux(k,j+1)=0.5d0*(ucn1+ucn2)
            u0(k,kk-numnp)=(ucn1+ucn2)/2.d0
          enddo
          if(nfunc15.eq.0.and.ndm.eq.6)then
              xx(1,j+1)=(xx(1,j+1)+xx(4,j+1))*0.5d0
              xx(2,j+1)=(xx(2,j+1)+xx(5,j+1))*0.5d0
              xx(3,j+1)=(xx(3,j+1)+xx(6,j+1))*0.5d0
              write(*,*)'u not evaluated for in ygenes
     +                   (nfunc15.eq.0.and.ndm.eq.6)'
          endif
        endif
      endif
c.... 2.kante
      ik1=abs(nk0(3))
      ik2=abs(nk0(4))
      call kaelem(ik1,ik2,i2,i,ikz,ike)
      if(nk0(3).lt.0)then
        call miknof(ik1,ik2,km,kmi,mikno)
        do jj=1,ndm
          xx(jj,3)= x0(jj,kmi-numnp)
        enddo
cnl>
        do jj=1,ndf
          ux(jj,3)= u0(jj,kmi-numnp)
        enddo
        if(ndm.eq.6)then
          xx(1,3)=(xx(1,3)+xx(4,3))*0.5d0
          xx(2,3)=(xx(2,3)+xx(5,3))*0.5d0
          xx(3,3)=(xx(3,3)+xx(6,3))*0.5d0
        endif
c.... ike,ikz modifizieren
        call ikemod (iek,ike,ikz,kmi,i,0,ke,kk,1,7,nen1)
        km2=kmi
      else
c.... ike,ikz erweitern
        call ikemod (iek,ike,ikz,kmi,i,i2,ke,kk,2,4,nen1)
        if(i2.eq.0)goto 131
        do jj=1,4
          if(iek(jj,i2).eq.ik1.or.iek(jj,i2).eq.ik2)goto 231
        enddo
231     if( jj.eq.1) then
          if(abs(iek(2,i2)).eq.ik1.or.abs(iek(2,i2)).eq.ik2)then
            iek(jj,i2)=iek(jj,i2)*(-1)
          else
            iek(4,i2)=iek(4,i2)*(-1)
          endif
        else
          iek(jj,i2)=iek(jj,i2)*(-1)
        endif
c.... kantenmittelknoten abspeichern
        km=km+1
        mikno(1,km) =ik1
        mikno(2,km) =ik2
        mikno(3,km) =kk+1
c.... knotennummern und knotenkoordinaten
131     kk=kk+1
        ianp0(kk)=naiter*(-1)
        km2=kk
c.... curve interpolation
        ncuv=0
        i1=abs(ik1)
        i2=abs(ik1)
        j1=abs(ik2)
        j2=abs(ik2)
        ia=3
        ib=4
        j=2
        ncuv=0
        ncuvr=0
        if(nshell.eq.1) ncuv=1
        call geosr
     1     (ia,ib,iek,x,x0,i,ndm,ncuvr,nen1,xxh,ih,iifr,iif,nshell,2)
        if(ncuv+ncuvr.ne.0)then
          nfunc1=1
        if(ncuv.eq.0.and.ncuvr.eq.1)then
c.... boundary geometry function
czr>
          if(ndm.eq.2)then
            xx(3,j+1)=0.0d0
            lav = .true.    ! undef
            nlav=2          ! undef
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1),
     1            nrt(iifr,1),nrt(iifr,1),  3 ,1)
          else
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1           ,nrt(iifr,1),nrt(iifr,1),ndm,1)
          endif
czr<
        elseif(ncuv.eq.1.and.ncuvr.eq.0)then
c.... only surface geometry function
          call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1        ,nrt(iif,1),nrt(iif,1),ndm,2)
        elseif(ncuv.eq.1.and.ncuvr.eq.2)then
c.... surface geometry function and plane boundery function
          call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1        ,nrt(iif,1),nrt(iifr,1),ndm,3)
        elseif(ncuv.eq.1.and.ncuvr.eq.3)then
c.... surface geometry function and boundary function
          call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1        ,nrt(iif,1),nrt(iifr,1),ndm,4)
        endif
         if(ndm.eq.2)then
            do k=1,2
              x0(k,kk-numnp)=xx(k,j+1)
            enddo
          elseif(ndm.eq.3)then
            do k=1,3
              x0(k,kk-numnp)=xx(k,j+1)
            enddo
          elseif(ndm.eq.6)then
czr         call dirsea(xxh,ia,ib,vn,1,tn)
            write(*,*)'czr call dirsea(x '
            do k=1,3
cba           write(*,*)'keine Direktorermittrung warnung'
czr              x0(k+3,kk-numnp)=xx(k,j+1)-tn*vn(k)*0.5d0
czr              x0(k,kk-numnp)=xx(k,j+1)+tn*vn(k)*0.5d0
            enddo
          endif
cnl.. 2406>
          do k=1,ndf
            call xcor(u,u0,ih(ia),k,ucn1,ndm)
            call xcor(u,u0,ih(ib),k,ucn2,ndm)
            ux(k,j+1)=0.5d0*(ucn1+ucn2)
            u0(k,kk-numnp)=(ucn1+ucn2)/2.d0
          enddo                                ! cnl<
        else
          do k=1,ndm
            call xcor(x,x0,ih(ia),k,xcn1,ndm)
            call xcor(x,x0,ih(ib),k,xcn2,ndm)
            xx(k,j+1)=0.5d0*(xcn1+xcn2)
            x0(k,kk-numnp)=(xcn1+xcn2)/2.d0
          enddo
C
cnl>   2401
            do k=1,ndf
              call xcor(u,u0,ih(ia),k,ucn1,ndm)
              call xcor(u,u0,ih(ib),k,ucn2,ndm)
              ux(k,j+1)=0.5d0*(ucn1+ucn2)
              u0(k,kk-numnp)=(ucn1+ucn2)/2.d0
            enddo
cnl<
C
          if(nfunc15.eq.0.and.ndm.eq.6)then
            xx(1,j+1)=(xx(1,j+1)+xx(4,j+1))*0.5d0
            xx(2,j+1)=(xx(2,j+1)+xx(5,j+1))*0.5d0
            xx(3,j+1)=(xx(3,j+1)+xx(6,j+1))*0.5d0
          endif
        endif
      endif
c.... koordinate elementmittelknoten ermitteln
      kk=kk+1
      kmm=kk
      ianp0(kk)=naiter
      if(nfunc1.eq.1)then
        if(nfunc15.eq.0)then
c.... only boundary geometry funtion
czr>
            if(ndm.eq.2)then
              do inxx=1,5
                xx(3,inxx)=0.0d0
              enddo
              call curve5 (xx,nrt(iifr,1),1)
            else
              call curve5 (xx,nrt(iifr,1),1)
            endif
        elseif(nfunc15.eq.1)then
c.... surface geometry function
          if(ndm.eq.2)then
            write(*,*)'*** surface geometry function on ndm=2 ' 
          else
            call curve5 (xx,nrt(iif,1),2)
          endif
        endif
        if(ndm.eq.2) then
          do k=1,ndm
            x0(k,kk-numnp)=xx(k,5)
          enddo
        elseif(ndm.eq.3) then
          do k=1,ndm
            x0(k,kk-numnp)=xx(k,5)
          enddo
        else if(ndm.eq.6)then
          ia=3
          ib=4
czr       call dirsea(xxh,ia,ib,vn,2,tn)
          do k=1,3
            write(*,*)'keine Direktorermittlung Warnung'
czr            x0(k+3,kk-numnp)=xx(k,5)-tn*vn(k)*0.5d0
czr            x0(k,kk-numnp)=xx(k,5)+tn*vn(k)*0.5d0
          enddo
        endif
cnl>
        do k=1,ndf
          u0(k,kk-numnp)=ux(k,j+1)
        enddo
cnl<
      else
        do k=1,ndm
          call xcor(x,x0,abs(iek(1,i)),k,xcn1,ndm)
          call xcor(x,x0,abs(iek(2,i)),k,xcn2,ndm)
          call xcor(x,x0,abs(iek(3,i)),k,xcn3,ndm)
          call xcor(x,x0,abs(iek(4,i)),k,xcn4,ndm)
          x0(k,kk-numnp)=1/4.*(xcn1+xcn2+xcn3+xcn4)
        enddo
cnl>
        do k=1,ndf
          call xcor(u,u0,abs(iek(1,i)),k,ucn1,ndm)
          call xcor(u,u0,abs(iek(2,i)),k,ucn2,ndm)
          call xcor(u,u0,abs(iek(3,i)),k,ucn3,ndm)
          call xcor(u,u0,abs(iek(4,i)),k,ucn4,ndm)
          u0(k,kk-numnp)=1/4.*(ucn1+ucn2+ucn3+ucn4)
        enddo
cnl<
      endif
c.... aufbau der element knotenbeziehung
      ke=ke+1
      iael0(ke)=(abs(iael(i))+1)*(-1)
c....
      if(naiter.eq.1) iael0(ke)=iael0(ke)*(-1)
      iek(1,ke)=abs(nk0(1))
      iek(2,ke)=abs(nk0(2))
      iek(3,ke)=km1
      iek(4,ke)=kmm
      iek(nen1,ke)=abs(iek(nen1,i))
      iek(nen1-1,ke)=iek(nen1-1,i)
      iek(nen1-2,ke)=iek(nen1-2,i)
c.... calculate element error for new elements(adaf)
czr   if(fsteu.eq.1)then
czr     erro(ke)=-1.d0
czr   endif
      ke=ke+1
      iael0(ke)=(abs(iael(i))+1)*(-1)
      if(naiter.eq.1) iael0(ke)=iael0(ke)*(-1)
      iek(1,ke)=abs(nk0(1))
      iek(2,ke)=kmm
      iek(3,ke)=km2
      iek(4,ke)=abs(nk0(4))
      iek(nen1,ke)=abs(iek(nen1,i))
      iek(nen1-1,ke)=iek(nen1-1,i)
      iek(nen1-2,ke)=iek(nen1-2,i)
c.... ike,ikz modifizieren (ecken)
      call ikemod
     1     (iek,ike,ikz,kmi,i,abs(nk0(2)),ke,kk,abs(nk0(4)),5,nen1)
      call ikemod(iek,ike,ikz,kmi,i,i2,ke,kk,abs(nk0(1)),6,nen1)
c.... calculate element error for new elements(adaf)
czr   if(fsteu.eq.1)then
czr     erro(ke)=-1.d0
czr   endif
      iael0(i)=abs(iael(i))+1
      iek(1,i)=kmm
      iek(2,i)=km1
      iek(3,i)=abs(nk0(3))
      iek(4,i)=km2
      iek(nen1,i)=abs(iek(nen1,i))
      return
      end
c...........................................................end.SR ygen_nl
      subroutine xymod_nl
     +      (nst,i,u0,u,hst0,x0,x,iek,ikz,ike,
     +      mikno,nen1,nhist,ndm,ndf,
     +      ianp,ianp0,iael,iael0,naiter,erro,kf)
c ---------------------------------------------------------------------------
c.... kymods modifiziert das y-Element in ein k-Element
c.... nonlinear 
c ---------------------------------------------------------------------------
      USE cdata
      implicit double precision (a-h,o-z)
czr   logical lav,lavx
czr   dimension nk0(4),nk(5),mikno(3,*)
czr   dimension erro(*),erron(*),erron0(*)
czr   dimension xxh(4,9),xx(6,5),ih(4),vn(3)
czr   common /gefehl/ nfehler
czr   common /adap3/  n9e,asteu,fsteu,fmax,nerro
      dimension iek(nen1,*),ikz(*),ike(*),mikno(3,*),ih(4),
     +      x(ndm,*),x0(ndm,*),xxh(4,9),xx(6,5),
     +      ianp(*),ianp0(*),iael0(*),iael(*),
     +      u(ndf,*),u0(ndf,*),ux(6,5),           !.... nonlinear
     +      hst0(nhist,*),                        !.... plasti
c    +      hst0(nhist,*),stx(20,5),              !.... plasti
     +      erro(*)
czr  +      vn(3)
czru>
      logical lav
      common /shel1/  ndir,nlav,lav
czru<
      common /gener/ kk,km,ke
      common /curvedat/ cpar(20,8),nrt(20,3),ic,nbe,nn3
      nn1=0
      nst1=0
      nst0=nst
      do k=2,4,2
        if(abs(iek(k,i)).gt.nn1 ) nn1=abs(iek(k,i))
      enddo
      nn=nn1
      nn2=abs(iek(1,i))
      call kaelem(nn1,nn2,ne2,i,ikz,ike)
      ne3=i
      if(ne2.gt.ne3)then
        nc=ne3
        ne3=ne2
        ne2=nc
      endif
c.... modify first element                           11111111
c.... curves on shell surface
      nfunc1=0
      nshell=0
c.... detect geometry functions on element surfaces      
      call geos(iek,x,x0,ne2,ndm,nshell,nen1,xxh,ih,iif,1)
      if(nshell.eq.1) nfunc1=1
c....
      ik1=abs(iek(1,ne2))
      ik2=abs(iek(2,ne2))
c.... ermittle das zweite Nachbarelement
      call kaelem(ik1,ik2,i2,ne2,ikz,ike)
      if(iek(1,ne2).lt.0)then
c.... find already existing center nodes      
        call miknof(ik1,ik2,km,kmi,mikno)
        call ikemod (iek,ike,ikz,kmi,ne2,0,ke,kk,0,7,nen1)
        km1=kmi
      else
c.... ike,ikz erweitern
        call ikemod (iek,ike,ikz,kmi,ne2,i2,ke,kk,1,4,nen1)
        if(i2.eq.0)goto 130
        do jj=1,4
          if(iek(jj,i2).eq.ik1.or.iek(jj,i2).eq.ik2)goto 230
        enddo
230     if( jj.eq.1) then
        if(abs(iek(2,i2)).eq.ik1.or.abs(iek(2,i2)).eq.ik2)then
          iek(jj,i2)=iek(jj,i2)*(-1)
        else
          iek(4,i2)=iek(4,i2)*(-1)
        endif
      else
        iek(jj,i2)=iek(jj,i2)*(-1)
      endif
c.... kantenmittelknoten abspeichern
      km=km+1
      mikno(1,km) =ik1
      mikno(2,km) =ik2
      mikno(3,km) =kk+1
      nst1=1
c.... knotennummern und knotenkoordinaten
130   kk=kk+1
      ianp0(kk)=naiter*(-1)
      km1=kk
czr   if(nst.eq.1.and.fsteu.eq.1)then
c.... calculate nodal error for new nodes(adaf)
czr     call xerr(erron,erron0,abs(ik1),ecn1)
czr     call xerr(erron,erron0,abs(ik2),ecn2)
czr     erron0(kk-numnp)=1/2.*(ecn1+ecn2)
czr   endif
c.... curve interpolation
      ncuv=0
      i1=abs(ik1)
      i2=abs(ik1)
      j1=abs(ik2)
      j2=abs(ik2)
      do jm=1,4
        if(abs(iek(jm,ne2)).eq.i1)ia=jm
        if(abs(iek(jm,ne2)).eq.j1)ib=jm
        if(abs(iek(jm,ne2)).eq.i1)j=jm-1
      enddo
      ncuv=0
      ncuvr=0
      if(nshell.eq.1) ncuv=1
      call geosr
     1     (ia,ib,iek,x,x0,ne2,ndm,ncuvr,nen1,xxh,ih,iifr,iif,nshell,1)
      if(ncuv+ncuvr.ne.0)then
        nfunc1=1
        if(ncuv.eq.0.and.ncuvr.eq.1)then
c.... boundary geometry function
czr>
          if(ndm.eq.2)then
            xx(3,j+1)=0.0d0
            lav = .true.    ! undef
            nlav=2          ! undef
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1),
     1            nrt(iifr,1),nrt(iifr,1),      3   ,1)
          else
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1           ,nrt(iifr,1),nrt(iifr,1),ndm,1)
          endif
czr<
          elseif(ncuv.eq.1.and.ncuvr.eq.0)then
c.... only surface geometry function
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1            ,nrt(iif,1),nrt(iif,1),ndm,2)
          elseif(ncuv.eq.1.and.ncuvr.eq.2)then
c.... surface geometry function and plane boundary function
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1           ,nrt(iif,1),nrt(iifr,1),ndm,3)
          elseif(ncuv.eq.1.and.ncuvr.eq.3)then
c.... surface geometry function and boundary function
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1           ,nrt(iif,1),nrt(iifr,1),ndm,4)
          endif
          if((ndm.eq.2).or.(ndm.eq.3))then
czr       if(ndm.eq.3)then
            do k=1,ndm
              x0(k,kk-numnp)=xx(k,j+1)
            enddo
          elseif(ndm.eq.6)then
czr            call dirsea(xxh,ia,ib,vn,1,tn)
                write(*,*)'czr call dirsea(x '
            do k=1,3
              write(*,*)'keine Direktorermittlung'
czr              x0(k+3,kk-numnp)=xx(k,j+1)-tn*vn(k)*0.5d0
czr              x0(k,kk-numnp)=xx(k,j+1)+tn*vn(k)*0.5d0
            enddo
          endif
cnl>
          do k=1,ndf
            u0(k,kk-numnp)=ux(k,j+1)
          enddo
cnl<
        else
czr>
          if(ndm.eq.2)then
            xx(3,j+1)=0.0d0
c            lav = .true.    ! undef
c            nlav=2          ! undef
          endif
czr<
cba test.... no geometry function ih(ia),ih(ib) not defined
          ih(ia)=abs(iek(ia,ne2))
          ih(ib)=abs(iek(ib,ne2))
c....
          do k=1,ndm
            call xcor(x,x0,ih(ia),k,xcn1,ndm)
            call xcor(x,x0,ih(ib),k,xcn2,ndm)
            xx(k,j+1)=0.5d0*(xcn1+xcn2)
            x0(k,kk-numnp)=(xcn1+xcn2)/2.d0
          enddo
cnl>
          do k=1,ndf
            call xcor(u,u0,ih(ia),k,ucn1,ndm)
            call xcor(u,u0,ih(ib),k,ucn2,ndm)
            ux(k,j+1)=0.5d0*(ucn1+ucn2)
            u0(k,kk-numnp)=(ucn1+ucn2)/2.d0
          enddo
cnl<
        endif
      endif
      neck=abs(iek(1,ne2))
      iael0(ne2)=abs(iael(ne2))
      iael(ne2)=abs(iael(ne2))
cba   iek0(1,ne2)=km1
      iek(1,ne2)=km1
      do k=2,nen1
        iek(k,ne2)=(iek(k,ne2))
      enddo
c....
czr     if(nst.eq.1.and.fsteu.eq.1)then
c.... calculate element error for new elements (adaf)
czr      call xerr(erron,erron0,abs(iek(1,ne2)),ecn1)
czr      call xerr(erron,erron0,abs(iek(2,ne2)),ecn2)
czr      call xerr(erron,erron0,abs(iek(3,ne2)),ecn3)
czr      call xerr(erron,erron0,abs(iek(4,ne2)),ecn4)
czr      erro(ne2)=(ecn1+ecn2+ecn3+ecn4)/4.d0
czr      erro(ne2) = 2.d0                ! test zr     
czr      if(erro(ne2).gt.fmax) fmax=erro(ne2)
czr      elseif(fsteu.eq.1)then
czr         erro(ne2)=-1.d0
czr      else
czr      endif
c.... modify second element                                 222222
c.... curves on shell surface
      nfunc1=0
      nshell=0
      call geos(iek,x,x0,ne3,ndm,nshell,nen1,xxh,ih,iif,1)
      if(nshell.eq.1) nfunc1=1
c..
      ik1=abs(iek(4,ne3))
      ik2=abs(iek(1,ne3))
c.... find neighboring element
      call kaelem(ik1,ik2,i2,ne3,ikz,ike)
      if(iek(4,ne3).lt.0)then
        call miknof(ik1,ik2,km,kmi,mikno)
        call ikemod (iek,ike,ikz,kmi,ne3,0,ke,kk,0,7,nen1)
        km2=kmi
      else
c.... ike,ikz erweitern
        call ikemod (iek,ike,ikz,kmi,ne3,i2,ke,kk,1,4,nen1)
        if(i2.eq.0)goto 131
        do jj=1,4
          if(iek(jj,i2).eq.ik1.or.iek(jj,i2).eq.ik2)goto 231
        enddo
231     if( jj.eq.1) then
          if(abs(iek(2,i2)).eq.ik1.or.abs(iek(2,i2)).eq.ik2)then
            iek(jj,i2)=iek(jj,i2)*(-1)
          else
            iek(4,i2)=iek(4,i2)*(-1)
          endif
        else
          iek(jj,i2)=iek(jj,i2)*(-1)
        endif
c.... kantenmittelknoten abspeichern
        km=km+1
        mikno(1,km) =ik1
        mikno(2,km) =ik2
        mikno(3,km) =kk+1
        nst1=1
c.... nodal numbering and coordinates
c.... increase number of nodes
131     kk=kk+1
        ianp0(kk)=naiter*(-1)
        km2=kk
c.... y-element will be modified for compatability  nst=1
        fsteu = 0            !          1:= adaf
        if(nst.eq.1.and.fsteu.eq.1)then
c.... calculate nodal error for new nodes(adaf)
czr       call xerr(erron,erron0,abs(ik1),ecn1)
czr       call xerr(erron,erron0,abs(ik2),ecn2)
czr       erron0(kk-numnp)=1/2.*(ecn1+ecn2)
       endif
c.... curve interpolation
        ncuv=0
        i1=abs(ik1)
        i2=abs(ik1)
        j1=abs(ik2)
        j2=abs(ik2)
        do jm=1,4
          if(abs(iek(jm,ne3)).eq.i1)ia=jm
          if(abs(iek(jm,ne3)).eq.j1)ib=jm
          if(abs(iek(jm,ne3)).eq.i1)j=jm-1
        enddo
        ncuv=0
        ncuvr=0
        if(nshell.eq.1) ncuv=1
c.... check for boundary geometry function
        call geosr
     1     (ia,ib,iek,x,x0,ne3,ndm,ncuvr,nen1,xxh,ih,iifr,iif,nshell,1)
        if(ncuv+ncuvr.ne.0)then
          nfunc1=1
          if(ncuv.eq.0.and.ncuvr.eq.1)then
c.... boundary geometry funtion exist
czr>
          if(ndm.eq.2)then
            xx(3,j+1)=0.0d0
            lav = .true.    ! undef
            nlav=2          ! undef
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1),
     1            nrt(iifr,1),nrt(iifr,1),      3   ,1)
          else
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1           ,nrt(iifr,1),nrt(iifr,1),ndm,1)
          endif
czr<
czr        call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
czr     1         ,nrt(iifr,1),nrt(iifr,1),ndm,1)
          elseif(ncuv.eq.1.and.ncuvr.eq.0)then
c.... only surface geometry funtion
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1        ,nrt(iif,1),nrt(iif,1),ndm,2)
          elseif(ncuv.eq.1.and.ncuvr.eq.2)then
c.... surface geometry funtion and plane boundery function
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1        ,nrt(iif,1),nrt(iifr,1),ndm,3)
          elseif(ncuv.eq.1.and.ncuvr.eq.3)then
c.... surface geometry funtion and boundery function
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1        ,nrt(iif,1),nrt(iifr,1),ndm,4)
          endif
czr       if(ndm.eq.3)then
          if((ndm.eq.2).or.(ndm.eq.3))then
            do k=1,ndm
              x0(k,kk-numnp)=xx(k,j+1)
            enddo
          else if(ndm.eq.6)then
czr         call dirsea(xxh,ia,ib,vn,1,tn)
            write(*,*)'czr call dirsea(x '
            do k=1,3
               write(*,*)'keine Direktorermittlung SR xymod_nl'
czr              x0(k+3,kk-numnp)=xx(k,j+1)-tn*vn(k)*0.5d0
czr              x0(k,kk-numnp)=xx(k,j+1)+tn*vn(k)*0.5d0
            enddo
          endif
        else
c....  ......................................................ncuv+ncuvr.eq.0
cba test.... no geometry function ih(ia),ih(ib) not defined
          ih(ia)=abs(iek(ia,ne3))
          ih(ib)=abs(iek(ib,ne3))
c....
          do k=1,ndm
            call xcor(x,x0,ih(ia),k,xcn1,ndm)
            call xcor(x,x0,ih(ib),k,xcn2,ndm)
            xx(k,j+1)=0.5d0*(xcn1+xcn2)
            x0(k,kk-numnp)=(xcn1+xcn2)/2.d0
          enddo
cnl>
          do k=1,ndf
            call xcor(u,u0,ih(ia),k,ucn1,ndm)
            call xcor(u,u0,ih(ib),k,ucn2,ndm)
            ux(k,j+1)=0.5d0*(ucn1+ucn2)
            u0(k,kk-numnp)=(ucn1+ucn2)/2.d0
          enddo
        endif
      endif
c.... koordinate elementmittelknoten ermitteln
      iael0(ne3)=abs(iael(ne3))
      iael(ne3)=abs(iael(ne3))
cba   iek0(1,ne3)=km2
      iek(1,ne3)=km2
      do k=2,nen1
        iek(k,ne3)=(iek(k,ne3))
      enddo
      iek(4,ne3)=abs(iek(4,ne3))
czr       if(nst.eq.1.and.fsteu.eq.1)then
c.... calculate element error for new elements (adaf)
czr         call xerr(erron,erron0,abs(iek(1,ne3)),ecn1)
czr         call xerr(erron,erron0,abs(iek(2,ne3)),ecn2)
czr         call xerr(erron,erron0,abs(iek(3,ne3)),ecn3)
czr         call xerr(erron,erron0,abs(iek(4,ne3)),ecn4)
czr         erro(ne3)=(ecn1+ecn2+ecn3+ecn4)/4.d0
czr         if(erro(ne3).gt.fmax) fmax=erro(ne3)
czr       elseif(fsteu.eq.1)then
czr         erro(ne3)=-1.d0
czr       else
czr       endif
c..
      if(iek(2,ne2).lt.0.or.iek(3,ne2).lt.0.or.
     1   iek(2,ne3).lt.0.or.iek(3,ne3).lt.0)then
        nst=2
        nst1=2
      else
        if(nst0.ne.1)then
          iek(nen1,ne2)=abs(iek(nen1,ne2))
          iek(nen1,ne3)=abs(iek(nen1,ne3))
        endif
      endif
      nis=ke
      k=0
cba310   if(nst.gt.0)then
310   if(nst1.gt.0)then
        if(ikz(neck)+k.eq.ikz(neck+1))goto 300
        if(ike(ikz(neck)+k).lt.nis)  nis=ike(ikz(neck)+k)
        if(ike(ikz(neck)+k).eq.ne2) goto 320
        if(ike(ikz(neck)+k).eq.ne3) goto 320
        iek(nen1,ike(ikz(neck)+k))=abs(iek(nen1,ike(ikz(neck)+k)))*(-1)
320     k=k+1
cba     ne2=ike(ikz(neck)+1)
cba     ne3=ike(ikz(neck)+2)
        goto310
300     continue
      endif
c.... modify third element        33333333
      call ikemod (iek,ike,ikz,neck,ne2,ne3,ke,kk,nn,8,nen1)
      ke=ke+1
      iael0(ke)= abs(iael0(ne3))
      iek(1,ke)=abs(neck)
      iek(2,ke)=km1
      iek(3,ke)=nn
      iek(4,ke)=km2
      iek(nen1,ke)=abs(iek(nen1,ne3))
      iek(nen1-1,ke)=iek(nen1-1,ne3)
      iek(nen1-2,ke)=iek(nen1-2,ne3)
      fsteu= 0                         !  test 
      if(nst0.eq.1.and.fsteu.eq.1)then
c.... calculate element error for new elements (adaf)
czr        call xerr(erron,erron0,abs(iek(1,ke)),ecn1)
czr        call xerr(erron,erron0,abs(iek(2,ke)),ecn2)
czr        call xerr(erron,erron0,abs(iek(3,ke)),ecn3)
czr        call xerr(erron,erron0,abs(iek(4,ke)),ecn4)
czr        erro(ke)=(ecn1+ecn2+ecn3+ecn4)/4.d0
czr        if(erro(ke).gt.fmax) fmax=erro(ke)
        write(*,*) '**** check reme_nl2.f'
      elseif(fsteu.eq.1.or.nst0.eq.1)then
        erro(ke)=-1.d0
      else
      endif
      ianp0(neck)=abs(ianp(neck))
      if(nst.eq.2.and.nst0.ne.1)then
        i=nis-1
      elseif(nst1.eq.1.and.nst0.ne.1)then
        i=nis-1
      else
        if(i.eq.ne2) i=i+1
cba     erro(ne2)=-1.d0
      endif
      
czr   if(kf+1.eq.ct2)then
      if(kf+1.eq.2)then        ! temp zr
czr     write(*,*)'* check xymod_nl ct2 = 2, ne2,ne3', ne2,ne3
        erro(ne2)=-1.d0
        erro(ne3)=-1.d0
cba     erro(ke)=-1.d0
      endif
      return
      end
c........................................................... end xymod_nl
      subroutine l_gener
     +     (nst,i,u0,u,hst0,x0,x,iek,ikz,ike, 
     +      mikno,nen1,nhist,ndm,ndf,
     +      ianp,ianp0,iael,iael0,naiter,erro, erron0)
c ---------------------------------------------------------------------------
c.... subdivide 2 point line element to remain compatible
c.... to neighboring 4 node element
c.... km : number of incompatible edge middle nodes
c ---------------------------------------------------------------------------
      USE cdata
      implicit double precision (a-h,o-z)
      common /gener/ kk,km,ke
      dimension x(ndm,*),x0(ndm,*)
     +       ,u(ndf,*),u0(ndf,*)           !.... nonlinear
     +       ,hst0(nhist,*)                        !.... plasti
     +       ,iek(nen1,*),ikz(*),ike(*)
     +       ,nk(5),mikno(3,*)
     +       ,ianp(*),ianp0(*),iael0(*),iael(*)
     +       ,erro(*),erron0(*)
c......................................................................
c.... find and treat new edge middle node
c     do j=0,2
      do j=0,0
c.... unmark node
czr     ianp0(abs(iek(j+1,i)))=abs(ianp0(abs(iek(j+1,i))))
c.... nodes of actual element
        ik1=abs(iek(j+1,i))
        ik2=abs(iek(j+2,i))
c.... find already generated edge middle node
        call miknof(ik1,ik2,km,kmi,mikno)
c.... save  edge middle node
        nk(1) = kmi
c.... modify node element connectivity of edge middle node
        call ikemod (iek,ike,ikz,kmi,i,0,ke,kk,j,9,nen1)
      enddo
c.... increase number of elements
      ke = ke + 1
c.... treat original node
c.... modify node element connectivity for one of the original nodes
      do inode=1,2
        if (iek(inode,i).gt. 0 ) then
          ik3 =  abs(iek(inode,i))
          call ikemod (iek,ike,ikz,ik3,i,0,ke,kk,j,10,nen1)
        endif
      enddo
c
c.... save new element nodal connectivities
      iael0(i)=abs(iael(i))+1
      iael0(ke)=abs(iael0(i))
      iek(nen1,ke)=abs(iek(nen1,i))
c
      iek(1,i)=ik1
      iek(2,i)=nk(1)
c....
      iek(1,ke)=nk(1)
      iek(2,ke)=ik2
c....
      goto 998
c......................................................................
998   return
      end
c........................................................... end l_gener
      subroutine histnp(ix,dt,hst,
     +      hfield,
     +      shp,numnp,dv)
c-----------------------------------------------------------------------
c.... project history variables to nodal points
c.... nodal identities  ix()
c.... shape functions   shp()
c.... number of history variables  nh
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      common /adap1/ naiter,mnph,mnphdt,nhist,nhsw
      dimension dt(numnp),hst(numnp,*),ix(*),shp(3,*),hfield(*)
      if (nhsw.eq.1) then
        do 10 i = 1,4
          xsji = dv*shp(3,i)
          ii = abs(ix(i))
          if(ii.le.0) go to 10
          dt(ii) = dt(ii) + xsji
          do j=1,nhist
            hst(ii,j) = hst(ii,j) + hfield(j)*xsji
          enddo
10      continue
      elseif (nhsw.eq.2 )then
        call pzero(hfield,nhist)
        do 20 i = 1,4
          ii = abs(ix(i))
          if(ii.le.0) go to 20
          do j=1,nhist
c           hfield(j) = hfield(j) + shp(3,i)*hst(ii,j)
            hfield(j) = hfield(j) + shp(3,i)*hst(ii,j)
          enddo
20      continue
c        write(iow,100) (hfield(j),j=1,nhist)
      endif
c100   format(10(f7.3,1x))
      end
c............................................................ end histnp
      subroutine haverage(hdt,h2,numnp)
c-----------------------------------------------------------------------
c.... nodal averaging h    (compare SR pltstr)                        |
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension hdt(numnp),h2(numnp,*)
      common /adap1/ naiter,mnph,mnphdt,nhist,nhsw
      do 100 ii = 1,numnp
        dh = hdt(ii)
c if(dh.ne.0.0d0) then
        if(abs(dh).gt.0.0d0) then
        do 200 ih = 1,nhist
          h2(ii,ih) = h2(ii,ih)/dh
200     continue
      endif
100   continue
      return
      end
c........................................................... end haverage
      subroutine rearrange(hst,hst0,nhist,iswitch)
c -----------------------------------------------------------------------
c.... rearrange nodal values on history field                          |
c -----------------------------------------------------------------------
      USE cdata
      implicit double precision (a-h,o-z)
      dimension hst(numnp,*),hst0(nhist,*)
      if (iswitch.eq.1) then
        do m=1,nhist
          do n=1,numnp
            hst0(m,n)=hst(n,m)
          enddo
        enddo
      elseif (iswitch.eq.2) then
        do i=1,numnp
          do j=1,nhist
czr            hst0(i,j)=hst(j,i)
            hst(i,j)=hst0(j,i)
          enddo
        enddo
      else
        write(*,*) 'error in rearrange'
      endif
      end
c.......................................................... end rearrange
c.......................................................... SR.rest_adap
      subroutine rest_adap(fres,b,hst0,ix,
     +      ndm,ndf,nen1,nhist,nneq,isw,asc)
c----------------------------------------------------------------------
c.... restart files - input/output  in ascii for iasc=1
c----------------------------------------------------------------------
      USE arcl
      USE cdata
      USE ddata
      USE dirdat
      USE edgdat
      USE fdata
      USE hdata
      USE iodata
      USE iofile
      USE mdata
      USE ndata
      USE pcrit
      USE prlod
      USE psize
      USE subdt
      USE tdata
      USE doalloc
      implicit double precision (a-h,o-z)
      logical exst,sfl
      character fres*229,yorn*1,y*80
      dimension b(*),ix(nen1,*),hst0(nhist,*)
      common m(maxm)  
      iasc = asc
C.... set number of displacement vectors to write or read
      iu = 1      ! short u vector ( only u )
c     iu = 3      ! long  u vector ( u, du, ddu)
c.... check file status
1     inquire(file=fres,exist=exst)
      if(.not.exst.and.isw.eq.1) then
        write(iow,3002) fres
        if(ior.lt.0) then
          write(*,3002) fres
          write(*,3003)
10        read (*,1000,err=11,end=12) yorn
          goto  13
11        call  errclr ('RESTRT')
          goto  10
12        call  endclr ('RESTRT',yorn)
13        if(yorn.eq.'y' .or. yorn.eq.'Y') then
            write(*,3004)
20          read (*,1001,err=21,end=22) fres
            goto  1
21          call  errclr ('RESTRT')
            goto  20
22          call  endclr ('RESTRT',fres)
            go to 1
          end if
        end if
        return
      end if
c.... open file
      if(iasc.eq.0) then
        if(exst) then
          open(ios,file=fres,form='unformatted',status='old',err=23)
        else
          open(ios,file=fres,form='unformatted',status='new')
        end if
      elseif(iasc.ne.0) then
        if(exst) then
          open(ios,file=fres,form='formatted'  ,status='old',err=23)
        else
          open(ios,file=fres,form='formatted'  ,status='new')
        end if
      end if
      goto 24
23    write(*,25)
      return
25    format(' Error, the restart/save file is not un-/formatted!')
24    rewind ios
      if(iasc.eq.0) then
c....   unformatted files
c....   read restart files
        write(*,*)'this version is not in use'
        if(isw.eq.1) then
          read(ios) nnpo,nnlo,nnmo,ndmo,ndfo,nrt,fl(9)
          if((nnpo.eq.numnp).and.(nnlo.eq.numel).and.(nnmo.eq.nummat)
     1   .and.(ndmo.eq.ndm).and.(ndfo.eq.ndf)) then
            read(ios) ttim,(b(i),i=1,nneq*iu) 
c....       eigenvalues and eigenvectors    
cww         read(ios) mf,mq,mfmax
cww         if(md.ne.0) then
cww           call pseta(md,mq,   ipr,sfl)
cww           call pseta(mv,mq*neq,ipr,sfl)
cww           read(ios) (m(i),i=md,md+mq*ipr)
cww           read(ios) (m(i),i=mv,mv+mq*neq*ipr)
cww         end if
c....       arc length values
            read(ios) prop,rlnew,c0,cs1,cs2,ds0,r,det0,xn
c....       transient fields v,a,...
            if(fl(9)) then
              call ralloc(trans,nrt*nneq,'REMENL-Trans',sfl) 
              read(ios) trans
            end if
c....       history fields
            read(ios) ttim,isgh1,isgh3
            if(size(gh1).ne.isgh1.or.size(gh3).ne.isgh3) then
              write(*,*) 'History data error'
              stop
            end if
            read(ios) gh1(1:size(gh1)),gh2(1:size(gh2)),gh3(1:size(gh3))
            if(nde*ned.gt.0) then
              ii = ne5 + 5*m(ne1+numnp-1)*nde*ipr -1
              read(ios) edge1,edge2,edge3,edge4
            end if
            refl = .true.
c....       crack values
            read(ios) clfl
            if(.not.clfl) then
              read(ios) ncs
              clfl = .true.
              nc1 = 1
              nc2 = nc1 + numnp
              nc3 = nc2 + numnp*ncs       ! REMENL nc1 .. nc4
              nc4 = nc3 + numnp*2
              ncmax = nc1+numnp*(ncs+5)
              call ralloc(crit,ncmax,'Crit frac',clfl)
              read(ios) (crit,i=nc1,ncmax)
              clfl=.false.
            end if
c.....      director values (only 1), all values defined in input-file!
            if(ldir.eq.1) then
              mdirmax = 10*knode
              read(ios) (basea(i),i=1,mdirmax)
              call pmove(basea,basea(1+knode),knode)
            end if
          else
            call drawmess(
     +      ' **ERROR** Incorrect information in a restart',1,0)
          end if
        elseif(isw.eq.2) then
c....     save information for restart during iteration
          write(ios) numnp,numel,nummat,ndm,ndf,nrt,fl(9)
          write(ios) ttim,(b(i),i=1,nneq*iu) 

cww       mq = min0(mf+mf,mf+8,neq)
cww       write(ios) mf,mq,mfmax
cww       if(md.ne.0) then
cww         write(ios) (m(i),i=md,md+mq*ipr)
cww         write(ios) (m(i),i=mv,mv+mq*neq*ipr)
cww       end if
c....     arc length values
          write(ios) prop,rlnew,c0,cs1,cs2,ds0,r,det0,xn
c....     transient fields v,a,...
          if(fl(9)) then
            write(ios) trans
          end if
c....     history fields:  numel*(2*nh1+nh3) 
          write(ios) ttim,size(gh1),size(gh3)
          write(ios) gh1(1:size(gh1)),gh2(1:size(gh2)),gh3(1:size(gh3))
c....     edge data
          if(nde*ned.gt.0) then
            ii = ne5 + 5*m(ne1+numnp-1)*nde*ipr -1  !cww??
            write(ios) edge1,edge2,edge3,edge4
          end if
c.....    crack values
          write(ios) clfl
          if(.not.clfl) then
            write(ios) ncs
            ncmax = nc1+numnp*(ncs+5)
            write(ios) crit
          end if
c.....    director values (only 1)
          if(ldir.eq.1) then
            mdirmax = 10*knode
            write(ios) (basea(i),i=1,mdirmax)
          end if
        end if
      elseif(iasc.ne.0) then
c....   formatted files
c....   read restart files
        if(isw.eq.1) then
c....     general values,displacements
          read(ios,4020) y
          read(ios,4002) nnpo,nnlo,nnmo,ndmo,ndfo,nrt,fl(9),ttim
          if((nnpo.eq.numnp).and.(nnlo.eq.numel).and.(nnmo.eq.nummat)
     1    .and.(ndmo.eq.ndm).and.(ndfo.eq.ndf)) then
            read(ios,4020) y
            call readf(b,ndf,iu*numnp,ios)  
c....       eigenvalues and eigenvectors
cww         read(ios,4020) y
cww         read(ios,4005) mf,mq,mfmax
cww         if(md.ne.0) then
cww           call pseta(md,mq,    ipr,sfl,'REMENL-EW')
cww           !call pseta(mv,mq*neq,ipr,sfl,'REMENL-EV')
cww           read(ios,4020) y
cww           call readfi(eigd,1,mq,ios)
cww                read(ios,4020) y
cww           call readf(eigv,mq,neq,ios)
cww         end if
cww         call ralloc(eigv,mq*neq,'REMENL-EV',sfl)
c....       arc length values
            read(ios,4020) y
            read(ios,4009) prop,rlnew,c0,cs1,cs2,ds0,r,det0,xn
            if(fl(9)) then
            if(fl(9)) then
c....         transient fields v,a,...
              call ralloc(trans,nrt*nneq,'REMENL-TRANS',sfl) !cww?? arcl
              read(ios,4020) y
              call readf(trans,ndf,nrt*numnp,ios)
            end if
              read(ios,4020) y
              call readf(trans,ndf,nrt*numnp,ios)
            end if
c....       history fields
            read(ios,4020) y
c            read(ios,4013) ttim  !cww?? passt write?
c            if(allocated(gh1)) then
c              call readf(gh1,1,size(gh1),ios)
c              call readf(gh2,1,size(gh2),ios)
c            end if
c            if(allocated(gh3)) then
c              call readf(gh3,1,size(gh3),ios)
c            end if
            read(ios,4013) nhist
            call readf(hst0,nhist,iu*numnp,ios) ! hst0 ??? ww  
c....       edge data
            if(nde*ned.gt.0) then
              ii = ne5 + 5*m(ne1+numnp-1)*nde*ipr -1
              read(ios,4020) y
              call readf(edge4(ne4),1,ii-ne5,ios)
            end if
            refl = .true.
c....       crack values
            read(ios,4020) y
            read(ios,4016) clfl
            if(.not.clfl) then
              read(ios,4005) ncs
              clfl = .true.
              nc1 = 1
              nc2 = nc1 + numnp
              nc3 = nc2 + numnp*ncs       ! REMENL nc1 .. nc4
              nc4 = nc3 + numnp*2
              call ralloc(crit,nc4+numnp*2,'Crit frac',clfl)
              call readf (crit,numnp,ncs+5,ios)
            end if
c.....      director values (only 1)
            if(ldir.eq.1) then
              read(ios,4020) y
              call readf(basea,10,knode,ios)
              call pmove(basea,basea(1+knode),knode)
            end if
          else
            call drawmess(
     +       ' **ERROR** Incorrect information in a restart',1,0)
          end if
        else if(isw.eq.2) then
c....     save information for restart during iteration
c....     general values,displacements
          write(ios,4001)
          write(ios,4002) numnp,numel,nummat,ndm,ndf,nrt,fl(9),ttim
          write(ios,4003)
          call writef(b,ndf,iu*numnp,ios) 
c....     eigenvalues and eigenvectors
cww       mq = min0(mf+mf,mf+8,neq)
cww       write(ios,4004)
cww       write(ios,4005) mf,mq,mfmax
cww       if(md.ne.0) then
cww         write(ios,4006)
cww         call writef(eigd,1,mq,ios)
cww         write(ios,4007)
cww         call writef(eigv,mq,neq,ios)
cww       end if
c....     arc length values
          write(ios,4008)
          write(ios,4009) prop,rlnew,c0,cs1,cs2,ds0,r,det0,xn
c....     transient fields v and a
          if(fl(9)) then
            write(ios,4010)
            call writef(trans,ndf,nrt*numnp,ios)
          end if
c....     history fields  numel*(2*nhmax+nh3max)   
          write(ios,4012)
c          write(ios,4013) ttim
c          call writef(gh1,1,size(gh1),ios)
c          call writef(gh2,1,size(gh2),ios)
c          call writef(gh3,1,size(gh3),ios)
          write(ios,4013) nhist
          call writef(hst0,nhist,iu*numnp,ios)
c....     edge data
          if(nde*ned.gt.0) then
            ii = ne5 + 5*m(ne1+numnp-1)*nde*ipr -1
            write(ios,4014)
            call writef(edge4(ne4),1,ii-ne5,ios)
          end if
c....     crack values
          write(ios,4015)
          write(ios,4016) clfl
          if(.not.clfl) then
            write(ios,4005) ncs
            call writef(crit,numnp,ncs+5,ios)
          end if
c.....    director values (only 1)
          if(ldir.eq.1) then
            write(ios,4017)
            call writef(basea,10,knode,ios)
          end if
        end if
      end if
c.... close the file
      close(ios)
      return
1000  format(a1)
1001  format(a229)
3002  format(' **ERROR** Restart file ',a229,/,' does not exist')
3003  format(11x,'Specify new name for restart file? (y or n) >',$)
3004  format(11x,'New Restart File Name >',$)
4001  format(' numnp,numel,nummat,ndm,ndf,nrt,fl(9),ttim')
4002  format(6i7,l5,e12.5)
4003  format(' displacements (nneq)')
4004  format(' mf,mq,mfmax')
4005  format(5i7)
4006  format(' eigenvalues  (neq)')
4007  format(' eigenvectors (neq)')
4008  format(' prop,rlnew,c0,cs1,cs2,ds0,r,det0,xn')
4009  format(9e12.5)
4010  format(' transient values v,a,... stored with neq, length=nneq')
4012  format(' history field - nodal projection (number of history',
     +      ' variables) ')
4013  format(i3)
4014  format(' edge data')
4015  format(' crack data')
4016  format(l10)
4017  format(' director data')
4020  format(a80)
      end
c............................................................end SR rest_adap
c.............................................................   SR curven_nl
      subroutine curven_nl(ia,ib,xx,i0,x,y,z,
     +      u1,u2,ux,ndf,
     +      nr,nrr,ndm,nxst)
c ---------------------------------------------------------------------------
c.... computes a new node on the edge of an element
c.... plugin for nonlinear problems
c ---------------------------------------------------------------------------
      USE cdata
      implicit double precision (a-h,o-z)
      logical lav,lavx
czr   logical lav
      common /shel1/  ndir,nlav,lav
      common /curvedat/ cpar(20,8),nrt(20,3),ic,nbe,nn3
      dimension xx(4,9),i0(4),x(3)
czr      dimension x0(3),a(3),b(3),c(3),d(4)
      dimension x0(3),a(3),b(3),d(4),
czr      dimension xa(3),xaa(3),
     +          xa(3),xl(6,4),xk(3,4),vn(3,4),tn(4),
     +          f1(4),f2(4), b1(3),b2(3),vn0(3),
     +          u1(ndf),u2(ndf),ux(6)          !,ux(6,5)     !.... nonlinear
      tol=0.0001
c.... interpolate displacements
      do iu=1,ndf
        ux(iu)=(u1(iu)+u2(iu))/2.d0
      enddo
c.... calculate director plane  
      x0(1)=(xx(ia,7)+xx(ib,7))/2.d0
      x0(2)=(xx(ia,8)+xx(ib,8))/2.d0
      x0(3)=(xx(ia,9)+xx(ib,9))/2.d0
      if(ndm.eq.6)then
czr        a(1)=((xx(ia,1)+xx(ib,1))/2.d0) -x0(1)
czr        a(2)=((xx(ia,2)+xx(ib,2))/2.d0) -x0(2)
czr        a(3)=((xx(ia,3)+xx(ib,3))/2.d0) -x0(3)
      elseif(ndm.eq.3)then
        do 40 i =1,4
          do 40 ii =1,3
            xk(ii,i)=xx(i,ii+6)
40      continue
        if(nxst.eq.2)then
czr       write(*,*)'nxst',nxst
          call dir(ia,xk,vn0,nr)
          vn(1,ia)=vn0(1)
          vn(2,ia)=vn0(2)
          vn(3,ia)=vn0(3)
***        write(*,*) vn(1,ia),vn(2,ia),vn(3,ia)
          call dir(ib,xk,vn0,nr)
          vn(1,ib)=vn0(1)
          vn(2,ib)=vn0(2)
          vn(3,ib)=vn0(3)
***       write(*,*) vn(1,ib),vn(2,ib),vn(3,ib)
        else
czr          lavx = lav
czr          nlavx=nlav
          lav =.true.       ! eingefuehrt wg. undef.
          nlav=2
          call midcom(4,xl,xk,vn,tn,1)
czr       lav =lavx
czr       nlav=nlavx
        endif
        a(1)=(vn(1,ia)+vn(1,ib))/2.d0
        a(2)=(vn(2,ia)+vn(2,ib))/2.d0
        a(3)=(vn(3,ia)+vn(3,ib))/2.d0
      else
        stop
      endif
      if(a(1).eq.0.d0.and.a(2).eq.0.d0.and.a(3).eq.0.d0) stop
      b(1)=xx(ib,7)-x0(1)
      b(2)=xx(ib,8)-x0(2)
      b(3)=xx(ib,9)-x0(3)
c.... director plane function f1(1)*x+f1(2)*y+f1(3)*z+f1(4)=0
      if(nxst.eq.1)then
c.... boundary geometry
        f1(1)=a(1)
        f1(2)=a(2)
        f1(3)=a(3)
        call vnorm (f1,ylb)
        f1(4)=(f1(1)*x0(1)+f1(2)*x0(2)+f1(3)*x0(3))*(-1.d0)
      else if(nxst.eq.2.or.nxst.eq.4)then
c.... surface geometry or first interpolation for surface and boundery
        if(nxst.eq.4)then
          ibxx=1
        endif
        call vcross(a,b,f1)
        call vnorm (f1,ylb)
        f1(4)=(f1(1)*x0(1)+f1(2)*x0(2)+f1(3)*x0(3))*(-1.d0)
      else if(nxst.eq.3)then
c.... surface geometry-plane boundary condition
        f1(1)=cpar(nrr,1)
        f1(2)=cpar(nrr,2)
        f1(3)=cpar(nrr,3)
        f1(4)=cpar(nrr,4)
      endif
c.... middle plane function f2(1)*x+f2(2)*y+f2(3)*z+f2(4)=0
      f2(1)=2.d0*(xx(ia,7)-xx(ib,7))
      f2(2)=2.d0*(xx(ia,8)-xx(ib,8))
      f2(3)=2.d0*(xx(ia,9)-xx(ib,9))
      f2(4)=xx(ib,7)**2+xx(ib,8)**2+xx(ib,9)**2
     1      -xx(ia,7)**2-xx(ia,8)**2-xx(ia,9)**2
cba      write(*,*) f1(1),f1(2),f1(3),f1(4)
cba      write(*,*) f2(1),f2(2),f2(3),f2(4)
c.... calculate the line function of two planes x=d+lambda*c
      call vcross(f1,f2,b)
      call vnorm(b,ylb)
      d(1)=(xx(ia,7)+xx(ib,7))/2.d0
      d(2)=(xx(ia,8)+xx(ib,8))/2.d0
      d(3)=(xx(ia,9)+xx(ib,9))/2.d0
cba      write(*,*) b(1),b(2),b(3)
cba      write(*,*) d(1),d(2),d(3)
100   call geopoi(b,d,tol,x,nr)
      if(nxst.eq.4)then
        if(ibxx.eq.1) goto 200
        diff0= abs(xa(1)-x(1))+abs(xa(2)-x(2))+abs(xa(3)-x(3))
        if(diff0.lt.tol) then
          nr=nr1
          goto 998                            ! return
        endif
200     xa(1)=x(1)
        xa(2)=x(2)
        xa(3)=x(3)
        ibxx=ibxx+1
        if(ibxx.eq.2) then
          nr1=nr
          nr2=nrr
          b1(1)=b(1)
          b1(2)=b(2)
          b1(3)=b(3)
          call vcross(f2,b,b2)
          call vnorm(b2,ylb)
        endif
        if(mod(ibxx,2).eq.0) then
          nr=nr2
          b(1)=b2(1)
          b(2)=b2(2)
          b(3)=b2(3)
        elseif(mod(ibxx,2).ne.0) then
          nr=nr1
          b(1)=b1(1)
          b(2)=b1(2)
          b(3)=b1(3)
        endif
        d(1)=x(1)
        d(2)=x(2)
        d(3)=x(3)
***     write(*,*) 'boundary iteration ',ibxx,diff0*1.d0
        if(ibxx.eq.150)goto 998
        goto 100
      endif
998   return
      end
c..............................................................end..SR curven
