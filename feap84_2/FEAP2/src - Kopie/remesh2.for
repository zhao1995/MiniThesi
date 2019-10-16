c ...................................................................SR kgenes
      subroutine kgenes(nst,i,x0,x,iek,ikz,ike,  mikno,nen1,ndm,
     1           ianp,ianp0,iael,iael0,naiter,erro, erron0)
c ---------------------------------------------------------------------------
c.... kgenes unterteilt das Element mit einem Kreuz
c ---------------------------------------------------------------------------
      USE cdata
      implicit double precision (a-h,o-z)
czru>
czr   logical lav,lavx
      logical lav
      common /shel1/  ndir,nlav,lav
czru<
      common /gener/ kk,km,ke
czr      common /gene0/ kk0,ke0
czr      common /gefehl/ nfehler
czr      common /adap3/  n9e,asteu,fsteu,fmax,nerro
      common /curvedat/  cpar(20,8),nrt(20,3),ic,nbe,nn3
      dimension x(ndm,*),iek(nen1,*),ikz(*),ike(*)
      dimension xx(6,5),nk(5),mikno(3,*)
      dimension x0(ndm,*),ianp(*),ianp0(*),iael0(*),iael(*)
czr   dimension erro(*),erron(*),erron0(*)
      dimension erro(*),         erron0(*)
czr   dimension m0(2),xxh(4,9),ih(4),vn(3)
      dimension       xxh(4,9),ih(4)
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
c.... knotennummern und knotenkoordinaten
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
      call geosr                                           ! in SR kgenes
     1     (ia,ib,iek,x,x0,i,ndm,ncuvr,nen1,xxh,ih,iifr,iif,nshell,1)
800   continue
      if(ncuv+ncuvr.ne.0)then
czr        write (*,*)'** (ncuv+ncuvr.ne.0) in kgenes'
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
c.... boundary geometry funtion
czr>
          if(ndm.eq.2)then
            xx(3,j+1)=0.0d0
czr111196            lav = .true.    ! undef
czr111196            nlav=2          ! undef
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1),
     1            nrt(iifr,1),nrt(iifr,1),      3   ,1)
          else
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1           ,nrt(iifr,1),nrt(iifr,1),ndm,1)
          endif
czr<
czr          call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
czr     1         ,nrt(iifr,1),nrt(iifr,1),ndm,1)
        elseif(ncuv.eq.1.and.ncuvr.eq.0)then
c.... only surface geometry function
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
c....
czr>
        if(ndm.eq.2)then
          do k=1,ndm
            x0(k,kk-numnp)=xx(k,j+1)
          enddo
czr<
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
        if(nfunc15.eq.0.and.ndm.eq.6)then
          xx(1,j+1)=(xx(1,j+1)+xx(4,j+1))*0.5d0
          xx(2,j+1)=(xx(2,j+1)+xx(5,j+1))*0.5d0
          xx(3,j+1)=(xx(3,j+1)+xx(6,j+1))*0.5d0
        endif
        if(j+1.eq.5.and.nfunc1.eq.1)then
          if(nfunc15.eq.0)then
c.... only boundary geometry funtion
czr>
czr         write(*,*)'call curve5 in SR kgener'
            if(ndm.eq.2)then
              xx(3,j+1)=0.0d0
              do izx=1,5
                xx(3,izx)=0.0d0
              enddo
              call curve5 (xx,nrt(iifr,1),1)
            else
              call curve5 (xx,nrt(iifr,1),1)
            endif
czr<
czr            call curve5 (xx,nrt(iifr,1),1)
          elseif(nfunc15.eq.1)then
c.... surface geometry funtion
              call curve5 (xx,nrt(iif,1),2)
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
C
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
C
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
c..................................................................end kgenes
c...................................................................SR ygenes
      subroutine ygenes
     1         (nst,i,x0,x,iek,ikz,ike   ,mikno,nen1,ndm,
     2           ianp,ianp0,iael,iael0,naiter,erro)
c ---------------------------------------------------------------------------
c.... ygenes unterteilt das element mit einem y
c ---------------------------------------------------------------------------
      USE cdata
      implicit double precision (a-h,o-z)
czru>
czr   logical lav,lavx
      logical lav
      common /shel1/  ndir,nlav,lav
czru<
      common /gener/ kk,km,ke
czr      common /gefehl/ nfehler
czr      common /adap3/  n9e,asteu,fsteu,fmax,nerro
      common /curvedat/ cpar(20,8),nrt(20,3),ic,nbe,nn3
      dimension x(ndm,*),iek(nen1,*),ikz(*),ike(*)
      dimension nk0(4),   mikno(3,*),erro(*)
czr      dimension nk0(4),nk(5),mikno(3,*),erro(*)
      dimension x0(ndm,*),ianp(*),ianp0(*),iael0(*),iael(*)
czr   dimension xxh(4,9),ih(4),xx(6,5),vn(3)
      dimension xxh(4,9),ih(4),xx(6,5)
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
c.... surface geometry funtion and plane boundary function
          call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1                 ,nrt(iif,1),nrt(iifr,1),ndm,3)
        elseif(ncuv.eq.1.and.ncuvr.eq.3)then
c.... surface geometry funtion and boundary function
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
            if(ndm.eq.2)then
              do k=1,2
                x0(k,kk-numnp)=xx(k,j+1)
              enddo
            elseif(ndm.eq.3)then
              do k=1,3
                 x0(k,kk-numnp)=xx(k,j+1)
              enddo
            else if(ndm.eq.6)then
czr              call dirsea(xxh,ia,ib,vn,1,tn)
czr              do k=1,3
czrcba             write(*,*)'keine Direktorermittrung Warnung'
czr                x0(k+3,kk-numnp)=xx(k,j+1)-tn*vn(k)*0.5d0
czr                x0(k,kk-numnp)=xx(k,j+1)+tn*vn(k)*0.5d0
czr              enddo
            endif
          else
            do k=1,ndm
              call xcor(x,x0,ih(ia),k,xcn1,ndm)
              call xcor(x,x0,ih(ib),k,xcn2,ndm)
              xx(k,j+1)=0.5d0*(xcn1+xcn2)
              x0(k,kk-numnp)=(xcn1+xcn2)/2.d0
            enddo
            if(nfunc15.eq.0.and.ndm.eq.6)then
              xx(1,j+1)=(xx(1,j+1)+xx(4,j+1))*0.5d0
              xx(2,j+1)=(xx(2,j+1)+xx(5,j+1))*0.5d0
              xx(3,j+1)=(xx(3,j+1)+xx(6,j+1))*0.5d0
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
     1            nrt(iifr,1),nrt(iifr,1),      3   ,1)
          else
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1           ,nrt(iifr,1),nrt(iifr,1),ndm,1)
          endif
czr<
czr        call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
czr     1         ,nrt(iifr,1),nrt(iifr,1),ndm,1)
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
        else if(ndm.eq.6)then
czr          call dirsea(xxh,ia,ib,vn,1,tn)
czr          do k=1,3
czrcba         write(*,*)'keine Direktorermittrung warnung'
czr            x0(k+3,kk-numnp)=xx(k,j+1)-tn*vn(k)*0.5d0
czr            x0(k,kk-numnp)=xx(k,j+1)+tn*vn(k)*0.5d0
czr          enddo
        endif
        else
          do k=1,ndm
            call xcor(x,x0,ih(ia),k,xcn1,ndm)
            call xcor(x,x0,ih(ib),k,xcn2,ndm)
            xx(k,j+1)=0.5d0*(xcn1+xcn2)
            x0(k,kk-numnp)=(xcn1+xcn2)/2.d0
          enddo
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
czr<
czr          call curve5 (xx,nrt(iifr,1),1)
        elseif(nfunc15.eq.1)then
c.... surface geometry funtion
czr>
          if(ndm.eq.2)then
czr            xx(3,j+1)=0.0d0
czr            do inxx=1,5
czr              xx(3,inxx)=0.0d0
czr            enddo
czr            call curve5 (xx,nrt(iif,1),2)
               write(*,*)'*** call curve5 (xx,nrt(iif,1),2)'
          else
            call curve5 (xx,nrt(iif,1),2)
          endif
czr<
czr          call curve5 (xx,nrt(iif,1),2)
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
czr          ia=3
czr          ib=4
czrczr          call dirsea(xxh,ia,ib,vn,2,tn)
czr          do k=1,3
czr***         write(*,*)'keine Direktorermittlung Warnung'
czr            x0(k+3,kk-numnp)=xx(k,5)-tn*vn(k)*0.5d0
czr            x0(k,kk-numnp)=xx(k,5)+tn*vn(k)*0.5d0
czr          enddo
        endif
      else
        do k=1,ndm
          call xcor(x,x0,abs(iek(1,i)),k,xcn1,ndm)
          call xcor(x,x0,abs(iek(2,i)),k,xcn2,ndm)
          call xcor(x,x0,abs(iek(3,i)),k,xcn3,ndm)
          call xcor(x,x0,abs(iek(4,i)),k,xcn4,ndm)
          x0(k,kk-numnp)=1/4.*(xcn1+xcn2+xcn3+xcn4)
        enddo
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
czr      if(fsteu.eq.1)then
czr        erro(ke)=-1.d0
czr      endif
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
czr      if(fsteu.eq.1)then
czr        erro(ke)=-1.d0
czr      endif
      iael0(i)=abs(iael(i))+1
      iek(1,i)=kmm
      iek(2,i)=km1
      iek(3,i)=abs(nk0(3))
      iek(4,i)=km2
      iek(nen1,i)=abs(iek(nen1,i))
      return
      end
c...............................................................end.SR ygenes
c...................................................................SR ogenes
      subroutine ogenes (nst,i,iek,nen1)
c ---------------------------------------------------------------------------
c.... ogenes demarkiert elemente ?
c ---------------------------------------------------------------------------
      USE cdata
      implicit double precision (a-h,o-z)
      dimension iek(nen1,*)
      if(nst.eq.1) iek (nen1,i)=abs(iek(nen1,i))
      return
      end
c...............................................................end.SR ogenes
c...................................................................SR kymods
      subroutine kymods
     1      (nst,i,x0,x,iek,ikz,ike,  mikno,nen1,ndm,
     2           ianp,ianp0,iael,iael0,naiter,
     3           erro  ,kf)
c ---------------------------------------------------------------------------
c.... kymods modifiziert das y-Element in ein k-Element
c ---------------------------------------------------------------------------
      USE cdata
      implicit double precision (a-h,o-z)
czru>
czr   logical lav,lavx
      logical lav
      common /shel1/  ndir,nlav,lav
czru<
      common /gener/ kk,km,ke
czr      common /gefehl/ nfehler
czr      common /adap3/  n9e,asteu,fsteu,fmax,nerro
      common /curvedat/ cpar(20,8),nrt(20,3),ic,nbe,nn3
      dimension x(ndm,*),iek(nen1,*),ikz(*),ike(*)
czr   dimension nk0(4),nk(5),mikno(3,*)
      dimension              mikno(3,*)
      dimension x0(ndm,*),ianp(*),ianp0(*),iael0(*),iael(*)
czr   dimension erro(*),erron(*),erron0(*)
      dimension erro(*)
czr   dimension xxh(4,9),xx(6,5),ih(4),vn(3)
      dimension xxh(4,9),xx(6,5),ih(4)
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
      call geos(iek,x,x0,ne2,ndm,nshell,nen1,xxh,ih,iif,1)
      if(nshell.eq.1) nfunc1=1
c....
      ik1=abs(iek(1,ne2))
      ik2=abs(iek(2,ne2))
      call kaelem(ik1,ik2,i2,ne2,ikz,ike)
      if(iek(1,ne2).lt.0)then
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
czr       if(nst.eq.1.and.fsteu.eq.1)then
c.... calculate nodal error for new nodes(adaf)
czr        call xerr(erron,erron0,abs(ik1),ecn1)
czr        call xerr(erron,erron0,abs(ik2),ecn2)
czr       erron0(kk-numnp)=1/2.*(ecn1+ecn2)
czr       endif
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
c.... boundary geometry funtion
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
c.... only surface geometry function
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1            ,nrt(iif,1),nrt(iif,1),ndm,2)
          elseif(ncuv.eq.1.and.ncuvr.eq.2)then
c.... surface geometry funtion and plane boundary function
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1           ,nrt(iif,1),nrt(iifr,1),ndm,3)
          elseif(ncuv.eq.1.and.ncuvr.eq.3)then
c.... surface geometry funtion and boundary function
            call curven (ia,ib,xxh,ih,xx(1,j+1),xx(2,j+1),xx(3,j+1)
     1           ,nrt(iif,1),nrt(iifr,1),ndm,4)
          endif
czr>
          if(ndm.eq.2)then
            do k=1,2
              x0(k,kk-numnp)=xx(k,j+1)
            enddo
czr<
          elseif(ndm.eq.3)then
            do k=1,3
              x0(k,kk-numnp)=xx(k,j+1)
            enddo
          else if(ndm.eq.6)then
czrczr            call dirsea(xxh,ia,ib,vn,1,tn)
czr            do k=1,3
czr**            write(*,*)'keine Direktorermittrung warnung'
czr              x0(k+3,kk-numnp)=xx(k,j+1)-tn*vn(k)*0.5d0
czr              x0(k,kk-numnp)=xx(k,j+1)+tn*vn(k)*0.5d0
czr            enddo
          endif
        else
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
        endif
      endif
      neck=abs(iek(1,ne2))
      iael0(ne2)=abs(iael(ne2))
      iael(ne2)=abs(iael(ne2))
*     iek0(1,ne2)=km1
      iek(1,ne2)=km1
      do k=2,nen1
        iek(k,ne2)=(iek(k,ne2))
      enddo
czr       if(nst.eq.1.and.fsteu.eq.1)then
c.... calculate element error for new elements (adaf)
czr        call xerr(erron,erron0,abs(iek(1,ne2)),ecn1)
czr        call xerr(erron,erron0,abs(iek(2,ne2)),ecn2)
czr        call xerr(erron,erron0,abs(iek(3,ne2)),ecn3)
czr        call xerr(erron,erron0,abs(iek(4,ne2)),ecn4)
czr       erro(ne2)=(ecn1+ecn2+ecn3+ecn4)/4.d0
czr       if(erro(ne2).gt.fmax) fmax=erro(ne2)
czr       elseif(fsteu.eq.1)then
czr       erro(ne2)=-1.d0
czr       else
czr       endif
c.... modify second element       222222
c.... curves on shell surface
      nfunc1=0
      nshell=0
      call geos(iek,x,x0,ne3,ndm,nshell,nen1,xxh,ih,iif,1)
      if(nshell.eq.1) nfunc1=1
c..
      ik1=abs(iek(4,ne3))
      ik2=abs(iek(1,ne3))
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
c.... knotennummern und knotenkoordinaten
131     kk=kk+1
        ianp0(kk)=naiter*(-1)
        km2=kk
czr       if(nst.eq.1.and.fsteu.eq.1)then
c.... calculate nodal error for new nodes(adaf)
czr        call xerr(erron,erron0,abs(ik1),ecn1)
czr        call xerr(erron,erron0,abs(ik2),ecn2)
czr       erron0(kk-numnp)=1/2.*(ecn1+ecn2)
czr       endif
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
        call geosr
     1     (ia,ib,iek,x,x0,ne3,ndm,ncuvr,nen1,xxh,ih,iifr,iif,nshell,1)
        if(ncuv+ncuvr.ne.0)then
          nfunc1=1
          if(ncuv.eq.0.and.ncuvr.eq.1)then
c.... boundary geometry funtion
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
czr>
          if(ndm.eq.2)then
            do k=1,2
              x0(k,kk-numnp)=xx(k,j+1)
            enddo
czr<
czr          if(ndm.eq.3)then
          elseif(ndm.eq.3)then
            do k=1,3
              x0(k,kk-numnp)=xx(k,j+1)
            enddo
          else if(ndm.eq.6)then
czrczr           call dirsea(xxh,ia,ib,vn,1,tn)
czr            do k=1,3
czr**         write(*,*)'keine Direktorermittrung warnung'
czr              x0(k+3,kk-numnp)=xx(k,j+1)-tn*vn(k)*0.5d0
czr              x0(k,kk-numnp)=xx(k,j+1)+tn*vn(k)*0.5d0
czr            enddo
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
        endif
      endif
c.... koordinate elementmittelknoten ermitteln
      iael0(ne3)=abs(iael(ne3))
      iael(ne3)=abs(iael(ne3))
*     iek0(1,ne3)=km2
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
**310   if(nst.gt.0)then
310   if(nst1.gt.0)then
        if(ikz(neck)+k.eq.ikz(neck+1))goto 300
        if(ike(ikz(neck)+k).lt.nis)  nis=ike(ikz(neck)+k)
        if(ike(ikz(neck)+k).eq.ne2) goto 320
        if(ike(ikz(neck)+k).eq.ne3) goto 320
        iek(nen1,ike(ikz(neck)+k))=abs(iek(nen1,ike(ikz(neck)+k)))*(-1)
320     k=k+1
*       ne2=ike(ikz(neck)+1)
*       ne3=ike(ikz(neck)+2)
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
czr       if(nst0.eq.1.and.fsteu.eq.1)then
c.... calculate element error for new elements (adaf)
czr         call xerr(erron,erron0,abs(iek(1,ke)),ecn1)
czr         call xerr(erron,erron0,abs(iek(2,ke)),ecn2)
czr         call xerr(erron,erron0,abs(iek(3,ke)),ecn3)
czr         call xerr(erron,erron0,abs(iek(4,ke)),ecn4)
czr         erro(ke)=(ecn1+ecn2+ecn3+ecn4)/4.d0
czr         if(erro(ke).gt.fmax) fmax=erro(ke)
czr       elseif(fsteu.eq.1.or.nst0.eq.1)then
czr         erro(ke)=-1.d0
czr       else
czr       endif
      ianp0(neck)=abs(ianp(neck))
      if(nst.eq.2.and.nst0.ne.1)then
        i=nis-1
      elseif(nst1.eq.1.and.nst0.ne.1)then
        i=nis-1
      else
        if(i.eq.ne2) i=i+1
cba     erro(ne2)=-1.d0
      endif
czr      if(kf+1.eq.ct2)then
czr     if(kf+1.eq.0.0d0)then        ! temp zr
        if(abs(kf+1).lt.0.00001)then        ! temp zr
          write(*,*)'**** check kymods ****'
          erro(ne2)=-1.d0
          erro(ne3)=-1.d0
***       erro(ke)=-1.d0
        endif
      return
      end
c............................................................end kymods
c......................................................................
      subroutine xcor(x,x0,node,k,xcn,ndm)
c ---------------------------------------------------------------------
c.... diese SR ermittelt die Knotenkoordinaten
c ---------------------------------------------------------------------
      USE cdata
      implicit double precision (a-h,o-z)
      dimension x(ndm,*),x0(ndm,*)
        if(node.le.numnp)then
          xcn=x(k,node)
        else
          xcn=x0(k,node-numnp)
        endif
        return
        end
c..............................................................SR xcor1
      subroutine xcor1(xl,x0,ndm)
c ---------------------------------------------------------------------
c.... diese SR transformiert die Knotenkoordinaten in
c.... lokalen Arbeitsspeicher
c ---------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension xl(10),x0(*)
      do 100 i=1,ndm
100   xl(i)=x0(i)
      return
      end
c......................................................................
      subroutine ikemod (iek,ike,ikz,ik,i,i2,ke,kk,j,nsw,nen1)
c ---------------------------------------------------------------------
c.... modify ike,ikz
c ---------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension ikz(*),ike(*),iek(nen1,*)
c.... kgener old
      if(nsw.eq.1)then
        n1=ikz(ik)
        n2=ikz(ik+1)
        do 100 ii=n1,n2-1
          if(ike(ii).eq.i     .and.j.ne.0) ike(ii)=ke+j
          if(ike(ii).eq.i*(-1).and.j.lt.3) ike(ii)=ke+j+1
          if(ike(ii).eq.i*(-1).and.j.eq.3) ike(ii)=i
100     continue
        goto 998       ! return
c.... kgener new
      elseif(nsw.eq.2)then
        if(j.eq.0)            ike(ikz(kk+1))=i
        if(j.gt.0.and.j.lt.4) ike(ikz(kk+1))=ke+j
        if(j.eq.3)            ike(ikz(kk+1)+1)=i
        if(j.lt.3)            ike(ikz(kk+1)+1)=ke+j+1
        if(i2.eq.0)then
          ikz(kk+2)=ikz(kk+1)+2
        else
          ike(ikz(kk+1)+2)=i2
          ike(ikz(kk+1)+3)=i2*(-1)
          ikz(kk+2)=ikz(kk+1)+4
        endif
        goto 998       ! return
c.... kgener oldmod
      elseif(nsw.eq.3)then
        ikz(kk+1)=ikz(kk)+4
        ike(ikz(kk))=i
        ike(ikz(kk)+1)=ke-2
        ike(ikz(kk)+2)=ke-1
        ike(ikz(kk)+3)=ke
        jj=2
        do 300 ii=2,4
          ix=abs(iek(ii,i))
          n1=ikz(ix)
          n2=ikz(ix+1)
          do 300 iii=n1,n2-1
            if(ike(iii).ne.i)goto 300
            ike(iii)=ke-jj
            jj=jj-1
300     continue
        goto 998       ! return
c.... ygener new
      elseif(nsw.eq.4)then
        ike(ikz(kk+1))=i
        ike(ikz(kk+1)+1)=ke+j
        if(i2.eq.0)then
          ikz(kk+2)=ikz(kk+1)+2
        else
          ike(ikz(kk+1)+2)=i2
          ike(ikz(kk+1)+3)=i2*(-1)
          ikz(kk+2)=ikz(kk+1)+4
        endif
        goto 998       ! return
c.... ygener oldmod
      elseif(nsw.eq.5)then
        ikz(kk+1)=ikz(kk)+3
        ike(ikz(kk))=i
        ike(ikz(kk)+1)=ke-1
        ike(ikz(kk)+2)=ke
        jj=1
        n1=ikz(i2)
        n2=ikz(i2+1)
        do 500 ii=1,2
        do 510 iii=n1,n2-1
        if(ike(iii).ne.i)goto 510
        ike(iii)=ke-jj
        jj=jj-1
510     continue
        n1=ikz(j)
        n2=ikz(j+1)
500     continue
        goto 998 ! return
c.... ygener oldmod
      elseif(nsw.eq.6)then
        do 600 ii=kk+1,j+1,-1
          ikz(ii)=ikz(ii)+1
600     continue
        do 610 ii=ikz(kk+1)-1,ikz(j+1),-1
          ike(ii)=ike(ii-1)
610     continue
        ike(ikz(j+1)-1)=ke
        n1=ikz(j)
        n2=ikz(j+1)
        do 620 iii=n1,n2-1
          if(ike(iii).ne.i)goto 620
          ike(iii)=ke-1
620     continue
        goto  998     ! return
c.... ygener old
      elseif(nsw.eq.7)then
        n1=ikz(ik)
        n2=ikz(ik+1)
        do 700 ii=n1,n2-1
          if(ike(ii).eq.i*(-1)) ike(ii)=ke+j+1
700     continue
        goto 998   ! return
c.... kymod new
      elseif(nsw.eq.8)then
        n1=ikz(ik)
        n2=ikz(ik+1)-1
        n3=ikz(j)
        n4=ikz(j+1)-1
        do 800 ii=n1,n2-1
          if(ike(ii).eq.i.or.ike(ii).eq.i2 ) goto 820
800     continue
820     ike(ii)=ke+1
        nst=0
        do 830 iii=ii+1,n2-1
          if(ike(iii).eq.i.or.ike(iii).eq.i2 ) nst=1
          if(nst.eq.1 ) ike(iii)=ike(iii+1)
830     continue
        do 840 ii =ik+1,j
          ikz(ii)=ikz(ii)-1
840     continue
        do 850 ii =n2,n3-2
          ike(ii)=ike(ii+1)
850     continue
        ike(n3-1)=ke+1
        goto 998                                               ! return
c.... modify   middle node   of line element
      elseif(nsw.eq.9)then
        n1=ikz(ik)                           ! new node of line element
        n2=ikz(ik+1)
        do ii=n1,n2-1                 ! j  local elem node (1 or 2)
          if(ike(ii).eq.i     .and.j.ne.0) ike(ii)=ke+j
          if(ike(ii).eq.i*(-1).and.j.lt.1) ike(ii)=ke+j+1
          if(ike(ii).eq.i*(-1).and.j.eq.1) ike(ii)=i
        enddo
        goto 998       ! return
c.... modify   original node  of line element
      elseif(nsw.eq.10)then
        n1=ikz(ik)                           ! new node of line element
        n2=ikz(ik+1)
        do ii=n1,n2-1
          if(ike(ii).eq.i) ike(ii) = ke      ! replace original element
        enddo
        goto 998       ! return
      endif
998   return
      end
c........................................................end..SR ikemod
      subroutine kaelem(k1,k2,ie2,ie1,ikz,ike)
c ---------------------------------------------------------------------------
c.... kaelem ermittelt das zweite nachbarelement zu einer kante
c.... ii1 : number of 'first' element connected to node k1
c.... ii2 : number of 'last'  element connected to node k1
c.... ie1 : actual element
c ---------------------------------------------------------------------------
      USE cdata
      implicit double precision (a-h,o-z)
      common /gener/ kk,km,ke
      dimension ike(*),ikz(*)
      ie2=0
      ii1=abs(ikz(k1))
      ii2=abs(ikz(k1+1))-1
      if(k1.eq.kk) ii2=nen*ke
      jj1=abs(ikz(k2))
      jj2=abs(ikz(k2+1))-1
      if(k2.eq.kk) jj2=nen*ke
      do 100 i=ii1,ii2       !loop elements connected k1
        do 100 j=jj1,jj2     !loop elements connected k2
          if(ike(i).eq.ie1)   goto 200
          if(ike(i).eq.ike(j))goto 300
200       continue
100   continue
      goto 998            ! return
300   ie2=ike(i)
998   return
      end
c..............................................................end..SR kaelem
c...................................................................SR miknof
      subroutine miknof (i1,i2,km,nn,mikno)
c ---------------------------------------------------------------------------
c.... miknof findet die bereits vorhandenen kantenmittelknoten
c.... km: number of incompatible element edges
c....     number of incompatible edge middle nodes
c ---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension mikno(3,*)
      do 100 i=1,km
        if(i1.eq.abs(mikno(1,i)).and. i2.eq.abs(mikno(2,i)))goto200
        if(i1.eq.abs(mikno(2,i)).and. i2.eq.abs(mikno(1,i)))goto200
100   continue
200   nn=mikno(3,i)
      do 300 j= i,km-1
        do 300 k=1,3
          mikno(k,j)=mikno(k,j+1)
300   continue
      km=km-1
      return
      end
c..............................................................end..SR miknof
c ---------------------------------------------------------------------------
      subroutine elnode(ix,ine,inez,nen1,nnp,nel,ie,num,ivec)
c ---------------------------------------------------------------------------
c.... 8.3.94 schle 07/98 zr
*---- Diese Routine ermittelt die Knoten-Element-Beziehungen INE aus
*---- den Element-Knoten-Beziehungen IX
*-----------------------------------------------------------------------------
*---- Felder       :
*---- IX (ien,iel) :  globaler Knoten, der dem lokalen Knoten IEN
*----                 des Elementes IEL entspricht
*---- NUM (node)   :  Anzahl der Elemente, mit denen der Knoten
*----                 NODE verbunden ist (maximal MAXNODE)
*---- IVEC         :  Zwischenablage der Elemente
*---- INEZ(i)      :  Anfangszeiger der Elementgruppe, die zum
*----                 i-ten Knoten geh”rt
*---- INE(k)       :  Endablage der Elemente in Vektorform
*---- maxnode      :  max. moegliche Anzahl der Knoten an einem Element
*-----------------------------------------------------------------------------

      USE cdata
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IX(NEN1,*),INE(*),INEZ(*),ie(8,*),num(*),ivec(*)
      maxnode = 10
      izeig   = 0
      do 100 i=1,nnp
100     num(i) = 0
      do 500    iel=1,nel
*---- nicht auszufhren fr Element 100 (gap-element)
        if( ie(7,ix(nen1,iel)) .eq.100 ) goto  500
        do 400 ien=1,nen
          node     = ix(ien,iel)
c.... for structures with different element types
Crol     neu eingefuegt:
           if (node.EQ.0) goto 500
          num(node)= num(node)+1
          izeig    = num(node)+(node-1)*maxnode
          ivec(izeig)=iel
400     continue
500   continue
      k=0
      inez(1)=1
      do 300 i=1,nnp
        ibeg=(i-1)*maxnode
         do 200 j=1,num(i)
           k=k+1
           ine(k)=ivec(ibeg+j)
200      continue
         if(i.gt.1) inez(i)=inez(i-1)+num(i-1)
300   continue
      inez(nnp+1)=k+1
      return
      end
c..............................................................end..SR elnode
c.....................................................................SR geos
      subroutine geos(iek,x,x0,n,ndm,ncuv,nen1,xx,i0,iif,nxst)
c ---------------------------------------------------------------------------
C.... detects geometry functions on element surfaces
c ---------------------------------------------------------------------------
      USE cdata
      USE iofile
      implicit double precision (a-h,o-z)
      common /curvedat/ cpar(20,8),nrt(20,3),ic,nbe,nn3
      dimension x(ndm,*),x0(ndm,*),iek(nen1,*)
      dimension xx(4,9),i0(4)
      dimension m0(2)
c....
      iif=0
      tol=0.0001
      do 50 i=1,ic-1
        if(nrt(i,3).eq.15)then
        do 70 ii=1,4
          if(nxst.eq.1) i0(ii)=abs(iek(ii,n))
          if(ndm.eq.6)then
c.... shell formulation with 6 dof
czr            do 80 iii=1,6
czr              call xcor(x,x0,i0(ii),iii,xa1,ndm)
czr              xx(ii,iii)=xa1
czr80          continue
czr            do 90 iii=1,3
czr              xx(ii,iii+6)=(xx(ii,iii)+xx(ii,iii+3))/2.d0
czr90          continue
          elseif(ndm.eq.3)then
c.... shell formulation with 3 dof
            do 91 iii=1,3
              call xcor(x,x0,i0(ii),iii,xa1,ndm)
              xx(ii,iii+6)=xa1
91          continue
          elseif(ndm.eq.2)then
c.... scheibe  formulation with 2 dof
            do iii=1,2
              call xcor(x,x0,i0(ii),iii,xa1,ndm)
              xx(ii,iii+6)=xa1
            enddo
             xx(ii,3+6)=0.0d0   ! dummy coordinate
          else
czr         write(*,*)  'ndm not specified - stop in SR geos'
            stop 'SR GEOS'
          endif
70      continue
        goto 60
        endif
50    continue
c.... no geometry function
      goto 120
czr      return
60    continue
c.... test geometry-functions
      do 100 i=1,ic-1
        if(nrt(i,3).eq.15)then
          nr=nrt(i,1)
          call curve(xx(1,7),xx(1,8),xx(1,9),m0,nrt(i,1),diff)
          if(abs(diff).gt.tol)goto 110
          call curve(xx(2,7),xx(2,8),xx(2,9),m0,nrt(i,1),diff)
          if(abs(diff).gt.tol)goto 110
          call curve(xx(3,7),xx(3,8),xx(3,9),m0,nrt(i,1),diff)
          if(abs(diff).gt.tol)goto 110
          call curve(xx(4,7),xx(4,8),xx(4,9),m0,nrt(i,1),diff)
          if(abs(diff).gt.tol)goto 110
          iif=i
          ncuv=1
        endif
110     continue
100   continue
120   continue
      return
      end
c................................................................end..SR geos
c....................................................................SR geosr
      subroutine geosr
     1   (ia,ib,iek,x,x0,n,ndm,ncuvr,nen1,xx,i0,iifr,iif,nshell,nxst)
c ---------------------------------------------------------------------------
c.... detects boundary functions for element edges
c ---------------------------------------------------------------------------
      USE cdata
      USE iofile
      implicit double precision (a-h,o-z)
      common /curvedat/ cpar(20,8),nrt(20,3),ic,nbe,nn3
      dimension x(ndm,*),x0(ndm,*),iek(nen1,*)
      dimension xx(4,9),i0(4)
      dimension m0(2)
c....
      tol=0.0001
      do 50 i=1,ic-1
        if(nrt(i,3).eq.16.and.nshell.eq.0)then
          do 70 ii=1,4
          if(nxst.eq.1) i0(ii)=abs(iek(ii,n))
          if(ndm.eq.6)then
c.... shell formulation with 6 dof
czr            do iii=1,6
czr              call xcor(x,x0,i0(ii),iii,xa1,ndm)
czr              xx(ii,iii)=xa1
czr            enddo
czr            do iii=1,3
czr              xx(ii,iii+6)=(xx(ii,iii)+xx(ii,iii+3))/2.d0
czr            enddo
          elseif(ndm.eq.3)then
c.... shell formulation with 3 dof
            do iii=1,3
              call xcor(x,x0,i0(ii),iii,xa1,ndm)
              xx(ii,iii+6)=xa1
            enddo
          elseif(ndm.eq.2)then
c.... scheibe  formulation with 2 dof
            do iii=1,2
              call xcor(x,x0,i0(ii),iii,xa1,ndm)
              xx(ii,iii+6)=xa1
            enddo
            xx(ii,3+6)=0.0d0     ! dummy coordinate
          else
czr         write(*,*)  'ndm not specified - stop in SR geos'
            stop 'SR GEOS'
          endif
70        continue
          goto 60
        endif
50    continue
c.... no geometry boundary function ...
60    continue
c.... test geometry-functions
      do 100 i=1,ic-1
        if(nrt(i,3).eq.16.and.nshell.eq.0.or.
     1     nrt(i,3).eq.17.and.nshell.eq.1.or.
     2     nrt(i,3).eq.18.and.nshell.eq.1.or.
     3     nrt(i,3).eq.15.and.nshell.eq.1.and.i.ne.iif)then
          nr=nrt(i,1)
          call curve(xx(ia,7),xx(ia,8),xx(ia,9),m0,nrt(i,1),diff)
          if(abs(diff).gt.tol)goto 110
          call curve(xx(ib,7),xx(ib,8),xx(ib,9),m0,nrt(i,1),diff)
          if(abs(diff).gt.tol)goto 110
          iifr=i
          ncuvr=1
          if(nrt(i,3).eq.18)ncuvr=2
          if(nrt(i,3).eq.17.or.nrt(i,3).eq.15) then
          ncuvr=3
          goto 120
        endif
      endif
110   continue
100   continue
120   continue
      return
      end
c...............................................................end..SR geosr
c...................................................................SR curven
      subroutine curven(ia,ib,xx,i0,x,y,z,nr,nrr,ndm,nxst)
c ---------------------------------------------------------------------------
c.... computes a new node on the edge of an element
c ---------------------------------------------------------------------------
      USE cdata
      implicit double precision (a-h,o-z)
      logical lav,lavx
czr   logical lav
      common /shel1/  ndir,nlav,lav
      common /curvedat/ cpar(20,8),nrt(20,3),ic,nbe,nn3
      dimension xx(4,9),i0(4),x(3)
czr      dimension x0(3),a(3),b(3),c(3),d(4)
      dimension x0(3),a(3),b(3),    d(4)
czr      dimension xa(3),xaa(3)
      dimension xa(3)
      dimension xl(6,4),xk(3,4),vn(3,4),tn(4)
      dimension f1(4),f2(4)
      dimension b1(3),b2(3),vn0(3)
      tol=0.0001
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
c          lavx = lav
c          nlavx=nlav
          lav =.true.       ! eingefuehrt wg. undef.
          nlav=2
          call midcom(4,xl,xk,vn,tn,1)
c          lav =lavx
c          nlav=nlavx
        endif
        a(1)=(vn(1,ia)+vn(1,ib))/2.d0
        a(2)=(vn(2,ia)+vn(2,ib))/2.d0
        a(3)=(vn(3,ia)+vn(3,ib))/2.d0
      else
        stop 'SR curven 1'
      endif
      if(a(1).eq.0.d0.and.a(2).eq.0.d0.and.a(3).eq.0.d0)
     +stop 'SR curven 2'
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
C...................................................................SR geopoi
      subroutine geopoi(b,x0,tol,xx,nr)
c ---------------------------------------------------------------------------
c.... iterative computation of a new node by functions
c ---------------------------------------------------------------------------
      USE cdata
      implicit double precision (a-h,o-z)
      common /curvedat/ cpar(20,8),nrt(20,3),ic,nbe,nn3
      dimension x0(3),b(3),xx(3),m0(2)
      diff=0.d0
      diffn=0.d0
      diffp=0.d0
      n1=-1
      i=0
      c=0.d0
100   continue
        if (abs(diff).gt.10000.d0) diff=100.d0
        if(i.eq.1) c=c+abs(diff)*n1/1000.d0
        if(i.gt.1) c=c+abs(diff)*n1*2
        x=x0(1)+c*b(1)
        y=x0(2)+c*b(2)
        z=x0(3)+c*b(3)
cba   call curves(x,y,z,m0,nr,diff)
        call curve(x,y,z,m0,nr,diff)
***   write(*,*) 'diff=',diff,'c=',c
        if(abs(diff).lt.tol/1000.d0) goto 200
        i=i+1
        if(i.gt.500) goto 700
        if(i.eq.1)then
          if(diff.lt.0.d0) cn=c
          if(diff.gt.0.d0) cp=c
          if(diff.lt.0.d0) diffn=diff
          if(diff.gt.0.d0) diffp=diff
          diff0=diff
      goto 100
        endif
        if(diff0*diff.lt.0.d0)goto 150
        if(i.eq.2.and.abs(diff0).lt.abs(diff))n1=1
        if(diff0.lt.0.d0.and.abs(diff).lt.abs(diffn))then
          diffn=diff
          cn=c
        endif
        if(diff0.gt.0.d0.and.abs(diff).lt.abs(diffp))then
          diffp=diff
          cp=c
        endif
      goto 100
150   continue
      if(diff.lt.0.d0)diffn=diff
      if(diff.gt.0.d0)diffp=diff
      if(diff.lt.0.d0)cn=c
      if(diff.gt.0.d0)cp=c
160   continue
      c=cn+(cp-cn)/(diffp-diffn)*(-1.d0)*diffn
***   write(*,*) c,cn,cp
      x=x0(1)+c*b(1)
      y=x0(2)+c*b(2)
      z=x0(3)+c*b(3)
cba   call curves(x,y,z,m0,nr,diff)
      call curve(x,y,z,m0,nr,diff)
      i=i+1
***   write(*,*) diff,c,diffn,diffp
      if(abs(diff).lt.tol/1000.d0) goto 200
c..
      if(diff.lt.0.d0.and.abs(diff).lt.abs(diffn))then
        diffn=diff
        cn=c
      endif
      if(diff.gt.0.d0.and.abs(diff).lt.abs(diffp))then
        diffp=diff
        cp=c
      endif
      i=i+1
      goto 160
200   continue
      xx(1)=x0(1)+c*b(1)
      xx(2)=x0(2)+c*b(2)
      xx(3)=x0(3)+c*b(3)
***   write(*,*)xx(1),xx(2),xx(3)
***   write(*,*)'nodal-iteration',i
      goto 998                     ! return
700   tol=-999.d0
998   continue
      return
      end
C..............................................................end..SR geopoi
c......................................................................SR dir
      subroutine dir (ia,xk,vn,nr)
c ---------------------------------------------------------------------------
c.... computes the director for a nodal point
c ---------------------------------------------------------------------------
      USE cdata
      implicit double precision (a-h,o-z)
      common /curvedat/ cpar(20,8),nrt(20,3),ic,nbe,nn3
      dimension xx(3,2),yy(3,2)
czr   dimension m0(2)
      dimension p1(3),p2(3),b(3)
      dimension p11(3),p22(3)
      dimension p33(3),p44(3)
      dimension y1(3,2),xk(3,4)
      dimension y2(3,2)
      dimension v1(3),v2(3),vn(3)
      dimension v3(3),v4(3)
      tol=0.001d0
      nnr=4
c
      xfak=1000.d0
      call clear(v1,1,3)
      call clear(v2,1,3)
      call clear(v3,1,3)
      call clear(v4,1,3)
      i1=ia+1
      if(i1.eq.5) i1=1
      i2=ia-1
      if(i2.eq.0) i2=4
      p1(1)=xk(1,i1)-xk(1,ia)
      p1(2)=xk(2,i1)-xk(2,ia)
      p1(3)=xk(3,i1)-xk(3,ia)
      p2(1)=xk(1,i2)-xk(1,ia)
      p2(2)=xk(2,i2)-xk(2,ia)
      p2(3)=xk(3,i2)-xk(3,ia)
      call vcross(p1,p2,b)
      call vnorm(b,xl)
c....
      y1(1,1)=xk(1,ia)+p1(1)/xfak
      y1(2,1)=xk(2,ia)+p1(2)/xfak
      y1(3,1)=xk(3,ia)+p1(3)/xfak
      y1(1,2)=xk(1,ia)-p1(1)/xfak
      y1(2,2)=xk(2,ia)-p1(2)/xfak
      y1(3,2)=xk(3,ia)-p1(3)/xfak
c....
      y2(1,1)=xk(1,ia)+p2(1)/xfak
      y2(2,1)=xk(2,ia)+p2(2)/xfak
      y2(3,1)=xk(3,ia)+p2(3)/xfak
      y2(1,2)=xk(1,ia)-p2(1)/xfak
      y2(2,2)=xk(2,ia)-p2(2)/xfak
      y2(3,2)=xk(3,ia)-p2(3)/xfak
      call geopoi(b,y1(1,1),tol,xx(1,1),nr)
      if(tol.lt.0.d0)then
        write(*,*)'no point'
        tol=0.001d0
        xx(1,1)=0.d0
        xx(2,1)=0.d0
        xx(3,1)=0.d0
      endif
      call geopoi(b,y1(1,2),tol,xx(1,2),nr)
      if(tol.lt.0.d0)then
        write(*,*)'no point'
        tol=0.001d0
        xx(1,2)=0.d0
        xx(2,2)=0.d0
        xx(3,2)=0.d0
      endif
      call geopoi(b,y2(1,1),tol,yy(1,1),nr)
      if(tol.lt.0.d0)then
        write(*,*)'no point'
        tol=0.001d0
        yy(1,1)=0.d0
        yy(2,1)=0.d0
        yy(3,1)=0.d0
      endif
      call geopoi(b,y2(1,2),tol,yy(1,2),nr)
      if(tol.lt.0.d0)then
        write(*,*)'no point'
        tol=0.001d0
        yy(1,2)=0.d0
        yy(2,2)=0.d0
        yy(3,2)=0.d0
      endif
      if(xx(1,1).eq.0.d0.and.xx(2,1).eq.0.d0.and.xx(3,1).eq.0.d0)then
        xl1=0.d0
      else
        p11(1)=xx(1,1)-xk(1,ia)
        p11(2)=xx(2,1)-xk(2,ia)
        p11(3)=xx(3,1)-xk(3,ia)
        xl1=dsqrt(p11(1)**2+p11(2)**2+p11(3)**2)
        call vnorm(p11,xl)
      endif
      if(xx(1,2).eq.0.d0.and.xx(2,2).eq.0.d0.and.xx(3,2).eq.0.d0)then
        xl2=0.d0
      else
        p22(1)=xx(1,2)-xk(1,ia)
        p22(2)=xx(2,2)-xk(2,ia)
        p22(3)=xx(3,2)-xk(3,ia)
        xl2=dsqrt(p22(1)**2+p22(2)**2+p22(3)**2)
        call vnorm(p22,xl)
      endif
      if(yy(1,1).eq.0.d0.and.yy(2,1).eq.0.d0.and.yy(3,1).eq.0.d0)then
        xl3=0.d0
      else
        p33(1)=yy(1,1)-xk(1,ia)
        p33(2)=yy(2,1)-xk(2,ia)
        p33(3)=yy(3,1)-xk(3,ia)
        xl3=dsqrt(p33(1)**2+p33(2)**2+p33(3)**2)
        call vnorm(p33,xl)
      endif
      if(yy(1,2).eq.0.d0.and.yy(2,2).eq.0.d0.and.yy(3,2).eq.0.d0)then
        xl4=0.d0
      else
        p44(1)=yy(1,2)-xk(1,ia)
        p44(2)=yy(2,2)-xk(2,ia)
        p44(3)=yy(3,2)-xk(3,ia)
        xl4=dsqrt(p44(1)**2+p44(2)**2+p44(3)**2)
        call vnorm(p44,xl)
      endif
      if(xl1.ne.0.d0.and.xl3.ne.0.d0)then
        anx1=1.d0
        call vcross(p11,p33,v1)
        call vnorm(v1,xl)
      endif
      if(xl2.ne.0.d0.and.xl3.ne.0.d0)then
        anx2=1.d0
        call vcross(p33,p22,v2)
        call vnorm(v2,xl)
      endif
      if(xl2.ne.0.d0.and.xl4.ne.0.d0)then
        anx3=1.d0
        call vcross(p22,p44,v3)
        call vnorm(v3,xl)
      endif
      if(xl1.ne.0.d0.and.xl4.ne.0.d0)then
        anx4=1.d0
        call vcross(p44,p11,v4)
        call vnorm(v4,xl)
      endif
      vn(1)=(v1(1)+v2(1)+v3(1)+v4(1))/(anx1+anx2+anx3+anx4)
      vn(2)=(v1(2)+v2(2)+v3(2)+v4(2))/(anx1+anx2+anx3+anx4)
      vn(3)=(v1(3)+v2(3)+v3(3)+v4(3))/(anx1+anx2+anx3+anx4)
      call vnorm(vn,xl)
      return
      end
c...................................................................SR midcom
      subroutine midcom( nk,xl,xk,vn,tn,ip )
c ---------------------------------------------------------------------------
c... nodal midplane coordinates,thickness and director
c ---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical lav
      common /el51/   vnn(3,4)
      common /shel1/  ndir,nlav,lav
      common /curvedat/ cpar(20,8),nrt(20,3),ic,nbe,nn3
      dimension xl(6,nk),xk(3,nk),vn(3,nk),tn(nk)
      dimension v1(3), v2(3) ,v3(3)
      dimension v0(3),m0(2)
      tol=0.0001d0
c
      if(lav) then
       if(nlav.eq.1)then
c... globale direktormittelung
        call matcop(xl, 3,4, xk)
        call matcop(vnn,3,4, vn)
c
       elseif(nlav.eq.2)then
c... direktor steht senkrecht auf dem knoten 1 keine mittelung
        if (ip.eq.0) call matcop(xl, 3,4, xk)

        v1(1) = xk(1,2) - xk(1,1)
        v1(2) = xk(2,2) - xk(2,1)
        v1(3) = xk(3,2) - xk(3,1)

        v2(1) = xk(1,4) - xk(1,1)
        v2(2) = xk(2,4) - xk(2,1)
        v2(3) = xk(3,4) - xk(3,1)
c
        call vcross (v1,v2,v3)
        call vnorm  (v3,dummy)
        do 90 i=1,4
          vn(1,i) = v3(1)
          vn(2,i) = v3(2)
          vn(3,i) = v3(3)
90      continue
c
       elseif(nlav.eq.3)then
c... direktor steht senkrecht auf dem jeweiligen knoten keine mittelung
C
        if (ip.eq.0) call matcop(xl, 3,4, xk)

c.. knoten 1
        v1(1) = xk(1,2) - xk(1,1)
        v1(2) = xk(2,2) - xk(2,1)
        v1(3) = xk(3,2) - xk(3,1)

        v2(1) = xk(1,4) - xk(1,1)
        v2(2) = xk(2,4) - xk(2,1)
        v2(3) = xk(3,4) - xk(3,1)
c
        call vcross (v1,v2,v3)
        call vnorm  (v3,dummy)
        vn(1,1) = v3(1)
        vn(2,1) = v3(2)
        vn(3,1) = v3(3)

c.. knoten2
        v1(1) = xk(1,3) - xk(1,2)
        v1(2) = xk(2,3) - xk(2,2)
        v1(3) = xk(3,3) - xk(3,2)

        v2(1) = xk(1,1) - xk(1,2)
        v2(2) = xk(2,1) - xk(2,2)
        v2(3) = xk(3,1) - xk(3,2)
c
        call vcross (v1,v2,v3)
        call vnorm  (v3,dummy)
        vn(1,2) = v3(1)
        vn(2,2) = v3(2)
        vn(3,2) = v3(3)

c.. knoten3
        v1(1) = xk(1,4) - xk(1,3)
        v1(2) = xk(2,4) - xk(2,3)
        v1(3) = xk(3,4) - xk(3,3)

        v2(1) = xk(1,2) - xk(1,3)
        v2(2) = xk(2,2) - xk(2,3)
        v2(3) = xk(3,2) - xk(3,3)
c
        call vcross (v1,v2,v3)
        call vnorm  (v3,dummy)
        vn(1,3) = v3(1)
        vn(2,3) = v3(2)
        vn(3,3) = v3(3)

c.. knoten4
        v1(1) = xk(1,1) - xk(1,4)
        v1(2) = xk(2,1) - xk(2,4)
        v1(3) = xk(3,1) - xk(3,4)

        v2(1) = xk(1,3) - xk(1,4)
        v2(2) = xk(2,3) - xk(2,4)
        v2(3) = xk(3,3) - xk(3,4)
c
        call vcross (v1,v2,v3)
        call vnorm  (v3,dummy)
        vn(1,4) = v3(1)
        vn(2,4) = v3(2)
        vn(3,4) = v3(3)
c
      elseif(nlav.eq.4)then
c... direktor steht senkrecht auf der geometriefunktion
c
       if (ip.eq.0) call matcop(xl, 3,4, xk)

        ncuv=0
        do 400 i=1,ic-1
          if(nrt(i,3).eq.15)then
            nr=nrt(i,1)
            call curve(xk(1,1),xk(2,1),xk(3,1),m0,nrt(i,1),diff)
            if(abs(diff).gt.tol)goto 410
            call curve(xk(1,2),xk(2,2),xk(3,2),m0,nrt(i,1),diff)
            if(abs(diff).gt.tol)goto 410
            call curve(xk(1,3),xk(2,3),xk(3,3),m0,nrt(i,1),diff)
            if(abs(diff).gt.tol)goto 410
            call curve(xk(1,4),xk(2,4),xk(3,4),m0,nrt(i,1),diff)
            if(abs(diff).gt.tol)goto 410
            iif=i
            ncuv=1
            goto 420
          endif
410   continue
400   continue
420   continue

       if (ncuv.eq.0)then
c.. no geometry function
        v1(1) = xk(1,2) - xk(1,1)
        v1(2) = xk(2,2) - xk(2,1)
        v1(3) = xk(3,2) - xk(3,1)
        v2(1) = xk(1,4) - xk(1,1)
        v2(2) = xk(2,4) - xk(2,1)
        v2(3) = xk(3,4) - xk(3,1)
        call vcross (v1,v2,v3)
        call vnorm  (v3,dummy)
        do 490 i=1,4
          vn(1,i) = v3(1)
          vn(2,i) = v3(2)
          vn(3,i) = v3(3)
490      continue
       else
c
c.. geometry function
        do 480 ia=1,4
          call dir(ia,xk,v0,iif)
         if(v0(1).eq.0.d0.and.v0(2).eq.0.d0.and.v0(3).eq.0.d0)then
          i1=ia+1
          if(i1.eq.5) i1=1
          i2=ia-1
          if(i2.eq.0) i2=4
          v1(1) = xk(1,i1) - xk(1,ia)
          v1(2) = xk(2,i1) - xk(2,ia)
          v1(3) = xk(3,i1) - xk(3,ia)
          v2(1) = xk(1,i2) - xk(1,ia)
          v2(2) = xk(2,i2) - xk(2,ia)
          v2(3) = xk(3,i2) - xk(3,ia)
          call vcross (v1,v2,v3)
          call vnorm  (v3,dummy)
          vn(1,ia) = v3(1)
          vn(2,ia) = v3(2)
          vn(3,ia) = v3(3)
          write(*,*)'warning director interpolation'
         else
           vn(1,ia)=v0(1)
           vn(2,ia)=v0(2)
           vn(3,ia)=v0(3)
         endif
480     continue
       endif
       endif
c...........................end...............
      else
        ii=0
        do 100 i=1,nk
        xk(1,i) = 0.5*(xl(1,i)+xl(4,i))
        xk(2,i) = 0.5*(xl(2,i)+xl(5,i))
        xk(3,i) = 0.5*(xl(3,i)+xl(6,i))
        vn(1,i) = xl(4,i)-xl(1,i)
        vn(2,i) = xl(5,i)-xl(2,i)
        vn(3,i) = xl(6,i)-xl(3,i)
        call vnorm  ( vn(1,i),tn(i) )
100     ii=ii+5
      endif
      end
c.............................................................end...SR midcom
c...................................................................SR curvem
      subroutine curvem(x1,x2,y1,y2,xx,yy,nr)
c ---------------------------------------------------------------------------
C.... this SR calculates the midpoint xx, yy of a function
c     between the points x1,y1 and x2,y2
C     for the curve number nr and the curve ityp
C     input :   x  = argument
C               nr = curve number
C     output:   y  = value at point x
C     variables cpar
C               nrt
c ---------------------------------------------------------------------------
      USE cdata
      USE iofile
      implicit double precision (a-h,o-z)
      common /curvedat/ cpar(20,8),nrt(20,3),ic,nbe,nn3
c
      do 100 i=1,40
        if( nrt(i,1) .eq. nr ) then
          ityp = nrt(i,2)
          nbeg = nrt(i,3)
          goto 200
        else
        endif
100   continue
      write(iow,2000) nr
      stop
200   continue
c
      if ( ityp .eq. 1 ) then
c
c.....Gerade durch zwei Punkte (x1,y1) (x2,y2)
c
       stop 'SR curvem 1'
      elseif ( ityp .eq. 2 ) then
c
C.....Gerade durch einen Punkt (x1,y1) und Steigung b
c
       stop 'SR curvem 2'
      elseif ( ityp .eq. 3 ) then
c
C.....Halbkreis mit dem Radius r, Steuerparameter nx und Mittelpunkt (x1,y1)
c
        r  = cpar(nr,1)
        xv = cpar(nr,2)
        if(abs(xv).ne.1.d0)xv=1.d0
        x0 = cpar(nr,3)
        y0 = cpar(nr,4)
        alf1=dacos((x1-x0)/r)
        alf2=dacos((x2-x0)/r)
        dalfn=(alf2-alf1)/2.d0+alf1
        xx=x0+dcos(dalfn)*r
        yy=y0+dsin(dalfn)*r*xv
        return
      elseif ( ityp .eq. 99 ) then
c
c....freie Kurve
C
       stop 'SR curven 3'
      else
      write(iow,2001) ityp
      stop 'SR curven 4'

      endif
      return
2000  format('***ERROR*** curve number ',i3,' in not been input')
2001  format('***ERROR*** curve type ',i3,' in not defined')
      end
c...............................................................end.SR curvem
c------------------------------------------------------------------------ende
