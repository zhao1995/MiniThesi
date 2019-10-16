      subroutine mate3d17(h1,h2,nh,d,md,Eps,Sig,aa,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c-----------------------------------------------------------------------
c     Purpose: calculate S and C for 
c     (c) f.k                                                   Okt,2013
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension h1(*),h2(*),d(*),Eps(6),Sig(6),eh(3),fv1(3),fv2(3),
     +          C3d(3,3),C3d2(3,3),ev(3,3),ewr(3),e(6),E3d(3,3),aa(6,6),
     +          aa3d(6,6),s3d(6),plout(10),xgp(3),delta(3,3)
c      
      if(isw.eq.1) then
c....   input data
      call mate3d16(h1,h2,nh,d,md,Eps,Sig,aa,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
      else
c     gestörte ew =0.d0 sonst 2.d0      
        xhlo=2.d0
        tol=1.d-15
        tol1=1.d-12
c...  speichern des rechten Cauchy-Green Tensors C 
c      do i=1,6
c        if (dabs(Eps(i)).lt.tol1) then
c          Eps(i)=0.0d0
c        endif
c      end do
      
        C3d(1,1) = Eps(1)*2.d0+1.d0
        C3d(2,2) = Eps(2)*2.d0+1.d0
        C3d(3,3) = Eps(3)*2.d0+1.d0
        C3d(1,2) = Eps(4)
        C3d(1,3) = Eps(5)
        C3d(2,3) = Eps(6)
        C3d(2,1) = C3d(1,2)
        C3d(3,1) = C3d(1,3)
        C3d(3,2) = C3d(2,3)
c...  speichern des rechten Cauchy-Green Tensors C
      do i=1,3
        do j=1,3
	   C3d2(i,j) = C3d(i,j)
	  end do
      end do
c      
c.... compute identity
      do i = 1,3
        do j = 1,3
          delta(i,j) = 0.d0
        end do
        delta(i,i) = 1.d0
      end do
c      
      call rsg(3,3,C3d,delta,ewr,1,ev,fv1,fv2,ierr)
      if (xhlo.ge.1.d0) then
              if (dabs(ewr(1)-ewr(2)).lt.(tol))then
              ewr(2)=(ewr(1)+ewr(2))/2.d0
              ewr(1)=ewr(2)
              end if
          if (dabs(ewr(1)-ewr(3)).lt.tol) then
              ewr(3)=(ewr(1)+ewr(3))/2.d0
              ewr(1)=ewr(3)
          end if
          if (dabs(ewr(2)-ewr(3)).lt.tol) then
              ewr(3)=(ewr(2)+ewr(3))/2.d0
              ewr(2)=ewr(3)
          end if
      else
      endif
c...  compute generalized strain tensor E = ln C
      do i = 1,3
          e(i) = dlog(ewr(i))
      end do
      if (xhlo.ge.1.d0) then
       if (abs(e(1)-e(2)).lt.(tol))then
              e(2)=(e(1)+e(2))/2.d0
              e(1)=e(2)
          end if
          if (abs(e(1)-e(3)).lt.tol) then
              e(3)=(e(1)+e(3))/2.d0
              e(1)=e(3)
          end if
          if (abs(e(2)-e(3)).lt.tol) then
              e(3)=(e(2)+e(3))/2.d0
              e(2)=e(3)
          end if
      else
      endif
c
      do i = 1,3
          eh(i)=e(i)
      end do     
c      
      do i = 1,3
        do j = 1,3
	  C3d(i,j)   = 0.d0
          do k = 1,3
            C3d(i,j)   = C3d(i,j)   + e(k)  *ev(i,k)*ev(j,k)
          enddo
        enddo
      enddo
c.... compute right Cauchy-Green tensor C = F^T F
      do i = 1,3
        do j = 1,3
            E3d(i,j) = C3d(i,j)/2.d0
        end do         
      end do
c
c.... restore generalized total green lagrangesche strain tensor        
       e(1) = E3d(1,1)
       e(2) = E3d(2,2)
       e(3) = E3d(3,3)
       e(4) = E3d(1,2)*2.d0 !Ingenieurgleitungen (Reihenfolge beachten)
       e(5) = E3d(1,3)*2.d0
       e(6) = E3d(2,3)*2.d0 
      call mate3d16(h1,h2,nh,d,md,e,Sig,aa,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c
      if (xhlo.ge.1.d0) then
          do i = 4,6
c            sig(i)=sig(i)*2.d0
          end do 
          call projsm17(ewr,ev,eh,sig,aa)
      else
          do i = 4,6
c            sig(i)=sig(i)*2.d0
          end do 
          call projsm17(ewr,ev,eh,sig,aa)   
      endif
c
      endif
c      
      return
      end
c
c
      subroutine projsm17(ewr,ev,eh,Sm3d,aa)
c---------------------------------------------------------------------72
c...    Projektion of stresses and moduli                 FG&JS 05.10.01
c       with special cases of equal eigenvalues: la_1 = la_2 = la_3 etc.
c       ewr(3)    eigen values of F^T * F           (I)
c       ev(3,3)   eigen vectors (normalized) of C   (I)
c       eh(3)     ln(ewr)                           (I)
c       Sm3d(6)   generalized  stresses             (I)
c       aa(6,6)   generalized moduli aa^m           (I)
c       SS3d(6)   2. PK stresses                    (O)
c       aa(6,6)   related moduli                    (O)
c-----------------------------------------------------------------------
      implicit none
      integer i,j,k
      real*8 ewr(3),ev(3,3),Sm3d(6),aa(6,6),ewrmi1(3),ewrmi2(3),eh(3)
      real*8 SS3d(6),xi(6,6),help4(6,6),xiaaxi(6,6),aa3djl1(6,6)
      real*8 TE(6,6),diagL(6),xih(6,6),trSNANB(6),rLmat(6,6),eei(3)
      real*8 fact,fact0,fact1,fact2,gama,fact7,fakt13
      real*8 gama112,gama221,gama113,gama331,gama223,gama332          
      real*8 xim,tol
c
      tol =1.d-15
c...  compute abbreviations 
      fact0 = 0.d0
      fact1 = -1.d0
      fact2 = -2.d0
      do i=1,3
        eei(i) = eh(i)/2.d0     
	  ewrmi1(i) =       ewr(i)**fact1
	  ewrmi2(i) = fact1*ewr(i)**fact2
      end do
c
c########ww: ftn95 salford requires integer exponents: 2 instead of 2.d0!!! 
      if(dabs(ewr(1)-ewr(2)) .gt. tol) then
        gama112 = ( ewrmi1(1) *(ewr(1)-ewr(2)) - 2.d0*(eei(1) - eei(2)))
     &           /(ewr(1) - ewr(2))**2
        gama221 = ( ewrmi1(2) *(ewr(2)-ewr(1)) - 2.d0*(eei(2) - eei(1)))
     &           /(ewr(2) - ewr(1))**2
      else
        gama112 = ewrmi2(1)/2.d0
        gama221 = ewrmi2(1)/2.d0
      endif
  
      if(dabs(ewr(1)-ewr(3)) .gt. tol) then
        gama113 = ( ewrmi1(1) *(ewr(1)-ewr(3)) - 2.d0*(eei(1) - eei(3)))
     &           /(ewr(1) - ewr(3))**2
        gama331 = ( ewrmi1(3) *(ewr(3)-ewr(1)) - 2.d0*(eei(3) - eei(1)))
     &           /(ewr(3) - ewr(1))**2
      else
        gama113 = ewrmi2(1)/2.d0
        gama331 = ewrmi2(1)/2.d0
      endif   

      if(dabs(ewr(2)-ewr(3)) .gt. tol) then     
        gama223 = ( ewrmi1(2) *(ewr(2)-ewr(3)) - 2.d0*(eei(2) - eei(3)))
     &           /(ewr(2) - ewr(3))**2
        gama332 = ( ewrmi1(3) *(ewr(3)-ewr(2)) - 2.d0*(eei(3) - eei(2)))
     &           /(ewr(3) - ewr(2))**2
      else
        gama223 = ewrmi2(2)/2.d0
        gama332 = ewrmi2(2)/2.d0
      endif   
c
      if(      dabs(ewr(1)-ewr(2)) .gt. tol 
     &   .and. dabs(ewr(2)-ewr(3)) .gt. tol
     &   .and. dabs(ewr(1)-ewr(3)) .gt. tol ) then
        gama = (ewr(1)*(eei(2)-eei(3))
     &         +ewr(2)*(eei(3)-eei(1))
     &         +ewr(3)*(eei(1)-eei(2)))
     &         /((ewr(1)-ewr(2)) *(ewr(2)-ewr(3)) *(ewr(3)-ewr(1)))
      elseif(dabs(ewr(1)-ewr(2)) .lt. tol) then
        gama = gama113/2.d0
      elseif(dabs(ewr(2)-ewr(3)) .lt. tol) then
        gama = gama221/2.d0
      elseif(dabs(ewr(3)-ewr(1)) .lt. tol) then
        gama = gama332/2.d0     
      endif
c
c.... Setup matrix TE
      TE(1,1) =       ev(1,1)*ev(1,1)
      TE(1,2) =       ev(2,1)*ev(2,1)
      TE(1,3) =       ev(3,1)*ev(3,1)
      TE(1,4) =       ev(1,1)*ev(2,1)
      TE(1,5) =       ev(1,1)*ev(3,1)
      TE(1,6) =       ev(2,1)*ev(3,1)
      TE(2,1) =       ev(1,2)*ev(1,2)
      TE(2,2) =       ev(2,2)*ev(2,2)
      TE(2,3) =       ev(3,2)*ev(3,2)
      TE(2,4) =       ev(1,2)*ev(2,2)
      TE(2,5) =       ev(1,2)*ev(3,2)
      TE(2,6) =       ev(2,2)*ev(3,2)
      TE(3,1) =       ev(1,3)*ev(1,3)
      TE(3,2) =       ev(2,3)*ev(2,3)
      TE(3,3) =       ev(3,3)*ev(3,3) 
      TE(3,4) =       ev(1,3)*ev(2,3)
      TE(3,5) =       ev(1,3)*ev(3,3)
      TE(3,6) =       ev(2,3)*ev(3,3)
      TE(4,1) = 2.0d0*ev(1,1)*ev(1,2)
      TE(4,2) = 2.0d0*ev(2,1)*ev(2,2)
      TE(4,3) = 2.0d0*ev(3,1)*ev(3,2)
      TE(4,4) =       ev(1,1)*ev(2,2) + ev(2,1)*ev(1,2)
      TE(4,5) =       ev(1,1)*ev(3,2) + ev(3,1)*ev(1,2)
      TE(4,6) =       ev(2,1)*ev(3,2) + ev(3,1)*ev(2,2)
      TE(5,1) = 2.0d0*ev(1,1)*ev(1,3)
      TE(5,2) = 2.0d0*ev(2,1)*ev(2,3)
      TE(5,3) = 2.0d0*ev(3,1)*ev(3,3)
      TE(5,4) =       ev(1,1)*ev(2,3) + ev(2,1)*ev(1,3)
      TE(5,5) =       ev(1,1)*ev(3,3) + ev(3,1)*ev(1,3)
      TE(5,6) =       ev(2,1)*ev(3,3) + ev(3,1)*ev(2,3)
      TE(6,1) = 2.0d0*ev(1,2)*ev(1,3)
      TE(6,2) = 2.0d0*ev(2,2)*ev(2,3)
      TE(6,3) = 2.0d0*ev(3,2)*ev(3,3)
      TE(6,4) =       ev(1,2)*ev(2,3) + ev(2,2)*ev(1,3)
      TE(6,5) =       ev(1,2)*ev(3,3) + ev(3,2)*ev(1,3)
      TE(6,6) =       ev(2,2)*ev(3,3) + ev(3,2)*ev(2,3)
C
c...  diag(L)
      diagL(1) = ewrmi1(1)
      diagL(2) = ewrmi1(2)
      diagL(3) = ewrmi1(3)
      if(dabs(ewr(1)-ewr(2)) .lt. tol) then
        diagL(4) = ewrmi1(1)/2.d0
      else
        diagL(4) = (eei(1)-eei(2))/(ewr(1)-ewr(2))
      endif
      if(dabs(ewr(1)-ewr(3)) .lt. tol) then
        diagL(5) = ewrmi1(1)/2.d0
      else
        diagL(5) = (eei(1)-eei(3))/(ewr(1)-ewr(3))
      endif
      if(dabs(ewr(2)-ewr(3)) .lt. tol) then
        diagL(6) = ewrmi1(2)/2.d0
      else  
        diagL(6) = (eei(2)-eei(3))/(ewr(2)-ewr(3))
      endif
c
c.... TE*^T diag(L) * TE
        do i = 1,6
          do j = 1,6
              xih(i,j) = diagL(i) * TE(i,j)
          end do
        end do
        do i = 1,6
          do j = 1,6
          xi(i,j)   = 0.d0
            do k = 1,6
              xi(i,j) = xi(i,j) + TE(k,i)* xih(k,j) 
            end do
          end do
        end do
c
c...  compute 2.PK-stresses: S = Xi : S^m
      do i=1,6
            SS3d(i)=0.d0
            fakt13=1.d0
            do j=1,6
            if (j.ge.4) fakt13=2.d0
            SS3d(i) =SS3d(i)+  xi(i,j)*Sm3d(j)*fakt13
            end do
      end do
c
c-----------------------------------------------------------------------
c...  compute Xi : IC^m : Xi = Xi * aa * Xi
c     compute help4 := aa * Xi
      do i=1,6
        do j=1,6
           help4(i,j) =  
     &         aa(i,1)*xi(1,j)+aa(i,2)*xi(2,j)+aa(i,3)*xi(3,j)+ 
     &        2.d0*(aa(i,4)*xi(4,j)+aa(i,5)*xi(5,j)+aa(i,6)*xi(6,j)) 
        end do
      end do
c...  compute Xi * help4
      do i=1,6
        do j=1,6
          xiaaxi(i,j) =  
     &     xi(i,1)*help4(1,j)+xi(i,2)*help4(2,j)+xi(i,3)*help4(3,j)+ 
     &  2.d0*(xi(i,4)*help4(4,j)+xi(i,5)*help4(5,j)+xi(i,6)*help4(6,j))
        end do
      end do
c
c...  compute S^m : 2*\partial_C (xi)
c     1. compute S^m : N_A \otmes N_B =: trSNANB(i)
      do i = 1,6
        trSNANB(i) = 0.d0
        do j = 1,6
          fact = 1.d0
          if(j .gt. 3) fact = 2.d0
          trSNANB(i) = trSNANB(i) + TE(i,j)*Sm3d(j)*fact
        end do
      end do
      do i = 4,6
        trSNANB(i) = trSNANB(i)/2.d0
      end do
c     2. compute L-matrix
        call pzero17jl(rLmat,6*6)
        rLmat(1,1) = 2.d0*trSNANB(1)*ewrmi2(1)
        rLmat(2,2) = 2.d0*trSNANB(2)*ewrmi2(2)
        rLmat(3,3) = 2.d0*trSNANB(3)*ewrmi2(3)
        rLmat(4,4) = trSNANB(1)*gama112 + trSNANB(2)*gama221
        rLmat(5,5) = trSNANB(1)*gama113 + trSNANB(3)*gama331
        rLmat(6,6) = trSNANB(2)*gama223 + trSNANB(3)*gama332
        rLmat(4,1) = 2.d0*trSNANB(4) * gama112
        rLmat(4,2) = 2.d0*trSNANB(4) * gama221
        rLmat(5,1) = 2.d0*trSNANB(5) * gama113
        rLmat(5,3) = 2.d0*trSNANB(5) * gama331
        rLmat(6,2) = 2.d0*trSNANB(6) * gama223
        rLmat(6,3) = 2.d0*trSNANB(6) * gama332
        rLmat(5,4) = 2.d0*trSNANB(6) * gama
        rLmat(6,4) = 2.d0*trSNANB(5) * gama
        rLmat(6,5) = 2.d0*trSNANB(4) * gama
c        
c      if(dabs(ewr(1)-ewr(2)) .lt. tol) then
c        rLmat(4,4) = (trSNANB(1)+ trSNANB(2))/2.d0*ewrmi2(1)
c        rLmat(4,1) = trSNANB(4) * ewrmi2(1)
c        rLmat(4,2) = trSNANB(4) * ewrmi2(1)
c        rLmat(5,4) =  trSNANB(6) * gama113
c        rLmat(6,4) =  trSNANB(5) * gama113
c        rLmat(6,5) =  trSNANB(4) * gama113
c      endif
c      if(dabs(ewr(1)-ewr(3)) .lt. tol) then
c        rLmat(5,5) = (trSNANB(1)+ trSNANB(3))/2.d0*ewrmi2(1)
c        rLmat(5,1) = trSNANB(5) * ewrmi2(1)
c        rLmat(5,3) = trSNANB(5) * ewrmi2(1)
c        rLmat(5,4) = trSNANB(6) * gama112
c        rLmat(6,4) = trSNANB(5) * gama112
c        rLmat(6,5) = trSNANB(4) * gama112
cc        
c      endif
c      if(dabs(ewr(2)-ewr(3)) .lt. tol) then
c        rLmat(6,6) = (trSNANB(2)+ trSNANB(3))/2.d0*ewrmi2(2)
c        rLmat(6,2) = trSNANB(6) * ewrmi2(2)
c        rLmat(6,3) = trSNANB(6) * ewrmi2(2)
c        rLmat(5,4) = trSNANB(6) * gama221
c        rLmat(6,4) = trSNANB(5) * gama221
c        rLmat(6,5) = trSNANB(4) * gama221
c      endif
c      if(      dabs(ewr(1)-ewr(2)) .lt. tol 
c     &   .and. dabs(ewr(2)-ewr(3)) .lt. tol) then
c        rLmat(5,4) = trSNANB(6)/2.d0*ewrmi2(1)
c        rLmat(6,4) = trSNANB(5)/2.d0*ewrmi2(1)
c        rLmat(6,5) = trSNANB(4)/2.d0*ewrmi2(1)
c      endif
c      
       do i = 1,6
         do j = i,6
             rLmat(i,j) = rLmat(j,i)
         end do 
       end do
c.... TE^T * Lmat * TE
        do i = 1,6
          do j = 1,6
            xih(i,j) = 0.d0 
            do k = 1,6
              xih(i,j) = xih(i,j) + rLmat(i,k) * TE(k,j)
            end do
          end do
        end do
        do i = 1,6
          do j = 1,6
          aa3djl1(i,j)   = 0.d0
            do k = 1,6
              aa3djl1(i,j) = aa3djl1(i,j) + TE(k,i)* xih(k,j) 
            end do
          end do
        end do

c...  compute aa = aa3djl1 + xiaaxi   
      do i=1,6
        do j=1,6        
	  aa(i,j) = aa3djl1(i,j) +xiaaxi(i,j)
c	  aa(i,j) = xiaaxi(i,j)
	end do
      end do	
c
c...  restore 2nd PK-stresses SS3d in Sm3d-array
      do i=1,6
        Sm3d(i) = SS3d(i)
      end do
c      
      return
      end
c
c
      subroutine projsmd17(ewr,ev,eh,Sm3d,aa)
c---------------------------------------------------------------------72
c...    Projektion of stresses and moduli                 FG&JS 05.10.01
c       with special cases of equal eigenvalues: la_1 = la_2 = la_3 etc.
c       ewr(3)    eigen values of F^T * F           (I)
c       ev(3,3)   eigen vectors (normalized) of C   (I)
c       eh(3)     ln(ewr)                           (I)
c       Sm3d(6)   generalized  stresses             (I)
c       aa(6,6)   generalized moduli aa^m           (I)
c       SS3d(6)   2. PK stresses                    (O)
c       aa(6,6)   related moduli                    (O)
c-----------------------------------------------------------------------
      implicit none
      integer i,j,k
      real*8 ewr(3),ev(3,3),Sm3d(6),aa(6,6),ewrmi1(3),ewrmi2(3),eh(3)
      real*8 SS3d(6),xi(6,6),help4(6,6),xiaaxi(6,6),aa3djl1(6,6)
      real*8 TE(6,6),diagL(6),xih(6,6),trSNANB(6),rLmat(6,6),eei(3)
      real*8 fact,fact0,fact1,fact2,gama
      real*8 gama112,gama221,gama113,gama331,gama223,gama332          
      real*8 xim,tol
c
      tol =1.d-15
c...  compute abbreviations 
      fact0 = 0.d0
      fact1 = -1.d0
      fact2 = -2.d0
      do i=1,3
        eei(i) = eh(i)/2.d0     
	  ewrmi1(i) =       ewr(i)**fact1
	  ewrmi2(i) = fact1*ewr(i)**fact2
      end do
c
c########ww: ftn95 salford requires integer exponents: 2 instead of 2.d0!!! 

        gama112 = ( ewrmi1(1) *(ewr(1)-ewr(2)) - 2.d0*(eei(1) - eei(2)))
     &           /(ewr(1) - ewr(2))**2
        gama221 = ( ewrmi1(2) *(ewr(2)-ewr(1)) - 2.d0*(eei(2) - eei(1)))
     &           /(ewr(2) - ewr(1))**2
        gama113 = ( ewrmi1(1) *(ewr(1)-ewr(3)) - 2.d0*(eei(1) - eei(3)))
     &           /(ewr(1) - ewr(3))**2
        gama331 = ( ewrmi1(3) *(ewr(3)-ewr(1)) - 2.d0*(eei(3) - eei(1)))
     &           /(ewr(3) - ewr(1))**2   
        gama223 = ( ewrmi1(2) *(ewr(2)-ewr(3)) - 2.d0*(eei(2) - eei(3)))
     &           /(ewr(2) - ewr(3))**2
        gama332 = ( ewrmi1(3) *(ewr(3)-ewr(2)) - 2.d0*(eei(3) - eei(2)))
     &           /(ewr(3) - ewr(2))**2
        gama = (ewr(1)*(eei(2)-eei(3))
     &         +ewr(2)*(eei(3)-eei(1))
     &         +ewr(3)*(eei(1)-eei(2)))
     &         /((ewr(1)-ewr(2)) *(ewr(2)-ewr(3)) *(ewr(3)-ewr(1)))
c
c.... Setup matrix TE
      TE(1,1) =       ev(1,1)*ev(1,1)
      TE(1,2) =       ev(2,1)*ev(2,1)
      TE(1,3) =       ev(3,1)*ev(3,1)
      TE(1,4) =       ev(1,1)*ev(2,1)
      TE(1,5) =       ev(1,1)*ev(3,1)
      TE(1,6) =       ev(2,1)*ev(3,1)
      TE(2,1) =       ev(1,2)*ev(1,2)
      TE(2,2) =       ev(2,2)*ev(2,2)
      TE(2,3) =       ev(3,2)*ev(3,2)
      TE(2,4) =       ev(1,2)*ev(2,2)
      TE(2,5) =       ev(1,2)*ev(3,2)
      TE(2,6) =       ev(2,2)*ev(3,2)
      TE(3,1) =       ev(1,3)*ev(1,3)
      TE(3,2) =       ev(2,3)*ev(2,3)
      TE(3,3) =       ev(3,3)*ev(3,3) 
      TE(3,4) =       ev(1,3)*ev(2,3)
      TE(3,5) =       ev(1,3)*ev(3,3)
      TE(3,6) =       ev(2,3)*ev(3,3)
      TE(4,1) = 2.0d0*ev(1,1)*ev(1,2)
      TE(4,2) = 2.0d0*ev(2,1)*ev(2,2)
      TE(4,3) = 2.0d0*ev(3,1)*ev(3,2)
      TE(4,4) =       ev(1,1)*ev(2,2) + ev(2,1)*ev(1,2)
      TE(4,5) =       ev(1,1)*ev(3,2) + ev(3,1)*ev(1,2)
      TE(4,6) =       ev(2,1)*ev(3,2) + ev(3,1)*ev(2,2)
      TE(5,1) = 2.0d0*ev(1,1)*ev(1,3)
      TE(5,2) = 2.0d0*ev(2,1)*ev(2,3)
      TE(5,3) = 2.0d0*ev(3,1)*ev(3,3)
      TE(5,4) =       ev(1,1)*ev(2,3) + ev(2,1)*ev(1,3)
      TE(5,5) =       ev(1,1)*ev(3,3) + ev(3,1)*ev(1,3)
      TE(5,6) =       ev(2,1)*ev(3,3) + ev(3,1)*ev(2,3)
      TE(6,1) = 2.0d0*ev(1,2)*ev(1,3)
      TE(6,2) = 2.0d0*ev(2,2)*ev(2,3)
      TE(6,3) = 2.0d0*ev(3,2)*ev(3,3)
      TE(6,4) =       ev(1,2)*ev(2,3) + ev(2,2)*ev(1,3)
      TE(6,5) =       ev(1,2)*ev(3,3) + ev(3,2)*ev(1,3)
      TE(6,6) =       ev(2,2)*ev(3,3) + ev(3,2)*ev(2,3)
C
c...  diag(L)
      diagL(1) = ewrmi1(1)
      diagL(2) = ewrmi1(2)
      diagL(3) = ewrmi1(3)
      diagL(4) = (eei(1)-eei(2))/(ewr(1)-ewr(2))/2.d0
      diagL(5) = (eei(1)-eei(3))/(ewr(1)-ewr(3))/2.d0
      diagL(6) = (eei(2)-eei(3))/(ewr(2)-ewr(3))/2.d0
c
c.... TE*^T diag(L) * TE
        do i = 1,6
          do j = 1,6
              xih(i,j) = diagL(i) * TE(i,j)
          end do
        end do
        do i = 1,6
          do j = 1,6
          xi(i,j)   = 0.d0
            do k = 1,6
              xi(i,j) = xi(i,j) + TE(k,i)* xih(k,j) 
            end do
          end do
        end do
c
c...  compute 2.PK-stresses: S = Xi : S^m
      do i=1,6
            SS3d(i)=0.d0
            do j=1,6
            SS3d(i) =SS3d(i)+  xi(i,j)*Sm3d(j)
            end do
      end do
c
c-----------------------------------------------------------------------
c...  compute Xi : IC^m : Xi = Xi * aa * Xi
c     compute help4 := aa * Xi
      do i=1,6
        do j=1,6
           help4(i,j) =  
     &         aa(i,1)*xi(1,j)+aa(i,2)*xi(2,j)+aa(i,3)*xi(3,j)+ 
     &         2.d0*(aa(i,4)*xi(4,j)+aa(i,5)*xi(5,j)+aa(i,6)*xi(6,j)) 
        end do
      end do
c...  compute Xi * help4
      do i=1,6
        do j=1,6
          xiaaxi(i,j) =  
     &     xi(i,1)*help4(1,j)+xi(i,2)*help4(2,j)+xi(i,3)*help4(3,j)+ 
     &  2.d0*(xi(i,4)*help4(4,j)+xi(i,5)*help4(5,j)+xi(i,6)*help4(6,j)) 
        end do
      end do
c
c...  compute S^m : 2*\partial_C (xi)
c     1. compute S^m : N_A \otmes N_B =: trSNANB(i)
      do i = 1,6
        trSNANB(i) = 0.d0
        do j = 1,6
          fact = 1.d0
          if(j .gt. 3) fact = 2.d0
          trSNANB(i) = trSNANB(i) + TE(i,j)*Sm3d(j)*fact
        end do
      end do
      do i = 4,6
        trSNANB(i) = trSNANB(i)/2.d0
      end do
c     2. compute L-matrix
        call pzero17jl(rLmat,6*6)
        rLmat(1,1) = -2.d0*trSNANB(1)*ewrmi2(1)
        rLmat(2,2) = -2.d0*trSNANB(2)*ewrmi2(2)
        rLmat(3,3) = -2.d0*trSNANB(3)*ewrmi2(3)
        rLmat(4,4) = trSNANB(1)*gama112 + trSNANB(2)*gama221
        rLmat(5,5) = trSNANB(1)*gama113 + trSNANB(3)*gama331
        rLmat(6,6) = trSNANB(2)*gama223 + trSNANB(3)*gama332
        rLmat(4,1) = 2.d0*trSNANB(4) * gama112
        rLmat(4,2) = 2.d0*trSNANB(4) * gama221
        rLmat(5,1) = 2.d0*trSNANB(5) * gama113
        rLmat(5,3) = 2.d0*trSNANB(5) * gama331
        rLmat(6,2) = 2.d0*trSNANB(6) * gama223
        rLmat(6,3) = 2.d0*trSNANB(6) * gama332
        rLmat(5,4) = 2.d0*trSNANB(6) * gama
        rLmat(6,4) = 2.d0*trSNANB(5) * gama
        rLmat(6,5) = 2.d0*trSNANB(4) * gama
       do i = 1,6
         do j = i,6
             rLmat(i,j) = rLmat(j,i)
         end do 
       end do
c.... TE^T * Lmat * TE
        do i = 1,6
          do j = 1,6
            xih(i,j) = 0.d0 
            do k = 1,6
              xih(i,j) = xih(i,j) + rLmat(i,k) * TE(k,j)
            end do
          end do
        end do
        do i = 1,6
          do j = 1,6
          aa3djl1(i,j)   = 0.d0
            do k = 1,6
              aa3djl1(i,j) = aa3djl1(i,j) + TE(k,i)* xih(k,j) 
            end do
          end do
        end do

c...  compute aa = aa3djl1 + xiaaxi
      do i=1,6
        do j=1,6
	  aa(i,j) = aa3djl1(i,j) + xiaaxi(i,j)
	end do
      end do	
c
c...  restore 2nd PK-stresses SS3d in Sm3d-array
      do i=1,6
        Sm3d(i) = SS3d(i)
      end do
c      
      return
      end
c

c--------------------------------------------------------------------72
      Subroutine dyad1_17 (v,vov)
      implicit double precision (a-h,o-z)
      dimension v(3),vov(6)
c
      vov(1) = v(1)*v(1)
      vov(2) = v(2)*v(2)
      vov(3) = v(3)*v(3)
      vov(4) = v(2)*v(3)
      vov(5) = v(3)*v(1)
      vov(6) = v(1)*v(2)
c
      return
      end
c--------------------------------------------------------------------72
c
      SUBROUTINE tqli17(d,e,n,np,z)
      INTEGER n,np
      DOUBLE PRECISION d(np),e(np),z(np,np)
CU    USES pythag17
      INTEGER i,iter,k,l,m
      DOUBLE PRECISION b,c,dd,f,g,p,r,s,pythag17
      do 11 i=2,n
        e(i-1)=e(i)
11    continue
      e(n)=0.d0
      do 15 l=1,n
        iter=0
1       do 12 m=l,n-1
          dd=abs(d(m))+abs(d(m+1))
          if (abs(e(m))+dd.eq.dd) goto 2
12      continue
        m=n
2       if(m.ne.l)then
          if(iter.eq.30)pause 'too many iterations in tqli17'
          iter=iter+1
          g=(d(l+1)-d(l))/(2.d0*e(l))
          r=pythag17(g,1.d0)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.d0
          c=1.d0
          p=0.d0
          do 14 i=m-1,l,-1
            f=s*e(i)
            b=c*e(i)
            r=pythag17(f,g)
            e(i+1)=r
            if(r.eq.0.d0)then
              d(i+1)=d(i+1)-p
              e(m)=0.d0
              goto 1
            endif
            s=f/r
            c=g/r
            g=d(i+1)-p
            r=(d(i)-g)*s+2.d0*c*b
            p=s*r
            d(i+1)=g+p
            g=c*r-b
C     Omit lines from here ...
            do 13 k=1,n
              f=z(k,i+1)
              z(k,i+1)=s*z(k,i)+c*f
              z(k,i)=c*z(k,i)-s*f
13          continue
C     ... to here when finding only eigenvalues.
14        continue
          d(l)=d(l)-p
          e(l)=g
          e(m)=0.d0
          goto 1
        endif
15    continue
      return
      END
c
c
c
      FUNCTION pythag17(a,b)
      DOUBLE PRECISION a,b,pythag17
      DOUBLE PRECISION absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag17=absa*sqrt(1.d0+(absb/absa)**2)
      else
        if(absb.eq.0.d0)then
          pythag17=0.d0
        else
          pythag17=absb*sqrt(1.d0+(absa/absb)**2)
        endif
      endif
      return
      END
c
c
c      
      SUBROUTINE tred17(a,n,np,d,e)
      INTEGER n,np
      DOUBLE PRECISION a(np,np),d(np),e(np)
      INTEGER i,j,k,l
      DOUBLE PRECISION f,g,h,hh,scale
      do 18 i=n,2,-1
        l=i-1
        h=0.d0
        scale=0.d0
        if(l.gt.1)then
          do 11 k=1,l
            scale=scale+abs(a(i,k))
11        continue
          if(scale.eq.0.d0)then
            e(i)=a(i,l)
          else
            do 12 k=1,l
              a(i,k)=a(i,k)/scale
              h=h+a(i,k)**2
12          continue
            f=a(i,l)
            g=-sign(sqrt(h),f)
            e(i)=scale*g
            h=h-f*g
            a(i,l)=f-g
            f=0.d0
            do 15 j=1,l
C     Omit following line if finding only eigenvalues
              a(j,i)=a(i,j)/h
              g=0.d0
              do 13 k=1,j
                g=g+a(j,k)*a(i,k)
13            continue
              do 14 k=j+1,l
                g=g+a(k,j)*a(i,k)
14            continue
              e(j)=g/h
              f=f+e(j)*a(i,j)
15          continue
            hh=f/(h+h)
            do 17 j=1,l
              f=a(i,j)
              g=e(j)-hh*f
              e(j)=g
              do 16 k=1,j
                a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
16            continue
17          continue
          endif
        else
          e(i)=a(i,l)
        endif
        d(i)=h
18    continue
C     Omit following line if finding only eigenvalues.
      d(1)=0.d0
      e(1)=0.d0
      do 24 i=1,n
C     Delete lines from here ...
        l=i-1
        if(d(i).ne.0.d0)then
          do 22 j=1,l
            g=0.d0
            do 19 k=1,l
              g=g+a(i,k)*a(k,j)
19          continue
            do 21 k=1,l
              a(k,j)=a(k,j)-g*a(k,i)
21          continue
22        continue
        endif
C     ... to here when finding only eigenvalues.
        d(i)=a(i,i)
C     Also delete lines from here ...
        a(i,i)=1.d0
        do 23 j=1,l
          a(i,j)=0.d0
          a(j,i)=0.d0
23      continue
C     ... to here when finding only eigenvalues.
24    continue
      return
      END
c
c
c
      subroutine pzero17jl(vvv,ib)
      implicit  none
      integer   ia,ib
      real*8    vvv(ib)
      save
      do ia = 1,ib
        vvv(ia) = 0.0d0
      end do
      end
c
c
c
      subroutine invert217(xa,ne,ndm)
c.... invert matrix
      implicit double precision (a-h,o-z)
      integer ne,p(20)
      dimension xa(ndm,ndm)
      if(ndm.gt.20) then
        write(*,*) 'invert2: ndm > 20. stop'
        stop
      endif
      n   = ne
      xtol = 1.d-8
      do 100 k = 1,n
        xmax = 0.d0
        p(k) = 0
        do 110 i = k,n
          xs = 0.d0
          do 120 j = k,n
            xs = xs + dabs( xa(i,j) )
  120     continue
          xq = dabs( xa(i,k) ) / xs
          if (xq.gt.xmax) then
            xmax = xq
            p(k) = i
          endif
  110   continue
        if (xmax.lt.xtol) stop
        if (p(k).ne.k) then
          do 130 j = 1,n
            xh         = xa(k,j)
            xa(k,j)    = xa(p(k),j)
            xa(p(k),j) = xh
  130     continue
        endif
        xpivot = xa(k,k)
        do 140 j = 1,n
          if (j.ne.k) then
            xa(k,j) = - xa(k,j) / xpivot
            do 145 i = 1,n
              if (i.ne.k) xa(i,j) = xa(i,j) + xa(i,k)*xa(k,j)
  145       continue
          endif
  140   continue
        do 150 i = 1,n
          xa(i,k) = xa(i,k) / xpivot
  150   continue
        xa(k,k) = 1.d0 / xpivot
  100 continue
      do 160 k = n-1,1,-1
        if (p(k).ne.k) then
          do 170 i = 1,n
            xh         = xa(i,k)
            xa(i,k)    = xa(i,p(k))
            xa(i,p(k)) = xh
  170     continue
        endif
  160   continue
      return
      end
c
c
c
c
      subroutine deta3x317(F,detf)
c...  compute determinant of 3 x 3 matrix      
      implicit none
      real*8 F(3,3),detf
      detf =        F(1,1)*(F(2,2)*F(3,3)-F(2,3)*F(3,2))
      detf = detf - F(1,2)*(F(2,1)*F(3,3)-F(2,3)*F(3,1)) 
      detf = detf + F(1,3)*(F(2,1)*F(3,2)-F(2,2)*F(3,1))
      return
      end
c