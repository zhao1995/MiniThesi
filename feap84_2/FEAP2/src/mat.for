C -----------------------------------------------------------------------
c.... Package of basic routines  
c
c     Vskal         Scalar product of to vectors v1(3), v2(3) 
c     Vcross        Cross  product of to vectors v1(3), v2(3) 
c     Vnorm         Norm of vectors v(3) 
c     Matadd        C = A+B
c     Mvadd         C = A+f*B
c     Matmin        C = A-B
c     Matmulf       C = A*B
c     Mvmul         C = A*b        
c     Mttmul        C = A^T*B
c     Mtbmul        C = A*B^T
c     Mttaml        C = C+A^T*B   
c     Mttmml        C = C-A^T*B
c     Clear         A = 0
c     cleari        I = 0 
c     Matcop        B = A  
c     Msmul         A = f*A
c     Mstiff        COMPUTE CONTRIBUTION TO ELEMENT STIFFNESS FOR ONE POINT
c     Msym          Symmetrize A  
c     Gauss         Gauss-scher Algorithmus   
c     Minmax        Find min and max of A
c     Matcoi        M = N
c     Pivot         Invert A with method of pivoting
c     daxpty        a = a + f * b 
c     skew          setup a skew-symmetric matrix w from vector a
c     conden        Condensation   of internal DOF's  using elimination
c     decond        Decondensation of internal DOF's using elimination 
c     conden1       De-/Condensation of internal dofs using invert
c     atri          Dreieckzerlegung  A = LDU 
c     asol          Gleichungslöser   A x = b  mit  A = LDU 
c     pmulst        TeilMultipl. mit Trafo K1=T^T*K(i,j), K2=K(i,j)*T
c
c
c.... history
C.....MAT   .FOR....GENERAL VERSION   : 12. 1.1989
C.....              LAST MODIFICATION : 17. 7.1990  KLA
C.....              MODIFIED 21.AUGUST 1989 SCHW/KLA
C.....              ADAPTED TO KNEBELS VERSION 19.10.92 SF
C-----------------------------------------------------------------------
      FUNCTION VSKAL ( V1,V2 )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION   V1(3),V2(3)
      VSKAL =  V1(1)*V2(1) + V1(2)*V2(2) + V1(3)*V2(3)
      END
C----------------------------------------------------------------------
      SUBROUTINE VCROSS  ( V1,V2, VN )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION   V1(3),V2(3),VN(3)
      VN(1) = V1(2)*V2(3) - V1(3)*V2(2)
      VN(2) = V1(3)*V2(1) - V1(1)*V2(3)
      VN(3) = V1(1)*V2(2) - V1(2)*V2(1)
      END
C----------------------------------------------------------------------
      SUBROUTINE VNORM ( V,RNORM )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION   V(3)
      RNORM = SQRT ( V(1)*V(1) + V(2)*V(2) + V(3)*V(3) )
      IF ( RNORM .LT. 1.E-15 )  RETURN
      V(1) = V(1) / RNORM
      V(2) = V(2) / RNORM
      V(3) = V(3) / RNORM
      END
C----------------------------------------------------------------------
      SUBROUTINE MATADD ( A,B,M,N, C )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*),B(*),C(*)
      DO 100 I=1,M*N
  100 C(I) = A(I) + B(I)
      END
C----------------------------------------------------------------------
      SUBROUTINE MAVADD ( A,B,M,N,FAK, C )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*),B(*),C(*)
      DO 100 I=1, M*N
  100 C(I) = A(I) + FAK*B(I)
      END
C----------------------------------------------------------------------
      SUBROUTINE MATMIN ( A,B,M,N, C )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*),B(*),C(*)
      DO 100 I=1,M*N
  100 C(I) = A(I) - B(I)
      END
C----------------------------------------------------------------------
      SUBROUTINE MATMULF ( A,B,L,M,N, C )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C.....MULTIPLIZIERT DIE MATRIZEN  C(L,N)=A(L,M)*B(M,N)

      DIMENSION A(L,M),B(M,N),C(L,N)
      DO 100 J=1,N
      DO 200 I=1,L
  200 C(I,J) = A(I,1) * B(1,J)
      DO 100 K=2,M
      DO 100 I=1,L
  100 C(I,J) = C(I,J) + A(I,K) * B(K,J)
      END
C----------------------------------------------------------------------
      SUBROUTINE MVMUL ( A,B,L,M, C )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.....MULTIPLIZIERT DEN VEKTOR MIT DER MATRIX C(L)=A(L,M)*B(M)
      DIMENSION A(L,M),B(M),C(L)
      DO 200 I=1,L
  200 C(I) = A(I,1) * B(1)
      DO 100 K=2,M
      DO 100 I=1,L
  100 C(I) = C(I) + A(I,K) * B(K)
      END
C----------------------------------------------------------------------
      SUBROUTINE MTTMUL ( A,B,L,M,N, C )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.....MULTIPLIZIERT DIE MATRIZEN  C(L,N)=AT(M,L)*B(M,N)
      DIMENSION A(M,L),B(M,N),C(L,N)
      DO 100 J=1,N
      DO 200 I=1,L
  200 C(I,J) = A(1,I) * B(1,J)
      DO 100 K=2,M
      DO 100 I=1,L
  100 C(I,J) = C(I,J) + A(K,I) * B(K,J)
      END
C----------------------------------------------------------------------
      SUBROUTINE MTBMUL ( A,B,L,M,N, C )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(L,M),B(N,M),C(L,N)
C.....MULTIPLIZIERT DIE MATRIZEN  C(L,N)=A(L,M)*BT(M,N)
      DO 100 J=1,N
      DO 200 I=1,L
  200 C(I,J) = A(I,1) * B(J,1)
      DO 100 K=2,M
      DO 100 I=1,L
  100 C(I,J) = C(I,J) + A(I,K) * B(J,K)
      END
C----------------------------------------------------------------------
      SUBROUTINE MTTAML ( A,B,L,M,N, C )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.....MULTIPLIZIERT DIE MATRIZEN  C(L,N)=C(L,N)+AT(M,L)*B(M,N)
      DIMENSION A(M,L),B(M,N),C(L,N)
      DO 100 J=1,N
      DO 100 K=1,M
      DO 100 I=1,L
  100 C(I,J) = C(I,J) + A(K,I) * B(K,J)
      END
C----------------------------------------------------------------------
      SUBROUTINE MTTMML ( A,B,L,M,N, C )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.....MULTIPLIZIERT DIE MATRIZEN  C(L,N)=C(L,N)-AT(M,L)*B(M,N)
      DIMENSION A(M,L),B(M,N),C(L,N)
      DO 100 J=1,N
      DO 100 K=1,M
      DO 100 I=1,L
  100 C(I,J) = C(I,J) - A(K,I) * B(K,J)
      END
C----------------------------------------------------------------------
      SUBROUTINE CLEAR ( A,M,N )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*)
      DO 100 I=1,M*N
  100 A(I) = 0.D0
      END
C----------------------------------------------------------------------
      SUBROUTINE CLEARI ( K,M,N )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION K(*)
      DO 100 I=1,M*N
  100 K(I) = 0
      END
C----------------------------------------------------------------------
      SUBROUTINE MATCOP ( A,M,N,B )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*),B(*)
      DO 100 I=1,M*N
  100 B(I) = A(I)
      END
C----------------------------------------------------------------------
      SUBROUTINE MSMUL ( A,M,N,S )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*)
      DO 100 I=1,M*N
  100 A(I) = A(I)*S
      END
C--------------------------------------------------------------------
      SUBROUTINE MSTIFF ( NST,B,C, S )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.... COMPUTE CONTRIBUTION TO ELEMENT STIFFNESS FOR ONE POINT, REPLACES :
C.... CALL MATMULF ( C,B,5,5,NST, B2 )
C.... CALL MTTAML ( B2,B,NST,5,NST, S )
C
      DIMENSION  B(5,NST),C(5,5),S(NST,NST)
C
      DO 400 J=1,5
      DO 400 I=1,NST
      H = B(1,I)*C(1,J)
      DO 300 K=2,5
300   H = H + B(K,I)*C(K,J)
      DO 400 L=I,NST
400   S(I,L) = S(I,L) + H*B(J,L)
      END
C----------------------------------------------------------------------
      SUBROUTINE MSYM ( A,K,N,ityp )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.... MAKES MATRIX A SYMMETRIC
      DIMENSION A(K,K)
      if(ityp.eq.1)then        ! ADD UPPER PART
      DO 500 I=1,N
      DO 500 L=1,I
500   A(L,I) = A(I,L)
      elseif(ityp.eq.2)then    ! ADD LOWER PART
      DO 600 I=1,N
      DO 600 L=1,I
600   A(I,L) = A(L,I)
      endif
      END
C----------------------------------------------------------------------
      SUBROUTINE GAUSS (M,N,A,B,DET,NDET,EPS,IERR)
************************************************************************
*     GLEICHUNGSLOESER NACH DEM GAUSS'SCHEN ALGORITHMUS   (A*X=L)      *
*     M=VEREINBARTE FELDGROESSEN A(M,M),B(M)                           *
*     N=ANZAHL DER UNBEKANNTEN (= TATSAECHLICHE FELDGROESSEN)          *
*     A=MATRIX                                                         *
*     B=LAST- UND LOESUNGSVEKTOR                                       *
*     DET=DETERMINANTE                                                 *
*     NDET=EXPONENT VON 10;DET MUSS MIT 10**NDET MULTIPLIZIERT WERDEN! *
*                                                                      *
*   AUFSTELLER       : DIPL.ING.D.H.MAIER , UNIVERSITAET KARLSRUHE     *
*   AUFGESTELLT      :  5.02.1985                                      *
*   LETZTE AENDERUNG : 15.04.1988  GE                                  *
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,M), B(M)

      IERR = 0
      DET=1.D0
      NDET=0

C     FALLS GLEICHUNGSSYSTEM DER GROESSE 0 :
      IF (N.LT.1)    RETURN

      DO 60 K=1,N-1
*
C     AUFSUCHEN DES PIVOTELEMENTES
*
        PIVOT=DABS(A(K,K))
        IZEILE=0
        DO 10 I=K+1,N
          IF (DABS(A(I,K)).LE.PIVOT) GO TO 10
          PIVOT=DABS(A(I,K))
          IZEILE=I
10        CONTINUE
        IF (IZEILE.EQ.0) GO TO 30
*
C   BEI JEDEM ZEILENTAUSCH AENDERT SICH DAS VORZEICHEN VON DET
*
        DET=-1.0*DET
        DO 20 I=K,N
          ZWSPEI=A(IZEILE,I)
          A(IZEILE,I)=A(K,I)
20        A(K,I)=ZWSPEI
        ZWSPEI=B(IZEILE)
        B(IZEILE)=B(K)
        B(K)=ZWSPEI
30      VER=DABS(A(K,K)/A(1,1))
C  KONTROLLAUSGABE  +++++++++++++++++++++++++++++++++++++++++++++++++++
C     WRITE(*,*) '('' K= '',I2,'' VER = '',G12.5,''  A(K,K)= '',G12.5,
C    &        '' A(1,1)= '',G12.5,'' IZEILE = '',I3//)',K,VER,A(K,K),
C    &        A(1,1),IZEILE
C  KONTROLLAUSGABE  +++++++++++++++++++++++++++++++++++++++++++++++++++

        IF (VER.LT.EPS) GOTO 90
*
C     ELIMINATION
*
        DO 50 I=K+1,N
          FAKTOR=A(I,K)/A(K,K)
          A(I,K)=0.D0
          DO 40 L=K+1,N
40          A(I,L)=A(I,L)-FAKTOR*A(K,L)
50        B(I)=B(I)-FAKTOR*B(K)
60      CONTINUE
*
C     BERECHNUNG DER UNBEKANNTEN UND DER DETERMINANTE
*
      VER=DABS(A(N,N)/A(1,1))
C  KONTROLLAUSGABE  +++++++++++++++++++++++++++++++++++++++++++++++++++
C     WRITE(*,*)'('' K= '',I2,'' VER = '',G12.5,''  A(K,K)= '',G12.5,
C    &        '' A(1,1)= '',G12.5,'' IZEILE = '',I3//)',K,VER,A(K,K),
C    &        A(1,1),IZEILE
C  KONTROLLAUSGABE  +++++++++++++++++++++++++++++++++++++++++++++++++++

      IF(VER.LT.EPS) GOTO 90
      DET=DET*A(N,N)
      B(N)=B(N)/A(N,N)
      DO 80 K=N-1,1,-1
      IF(DABS(DET).LT.1.D50) GOTO 100
           DET = DET/1.D50
           NDET = NDET+50
100     DET=DET*A(K,K)
        S=0.D0
        DO 70 I=K+1,N
70        S=S+B(I)*A(K,I)
80      B(K)=(B(K)-S)/A(K,K)
      RETURN

90    IERR = K
      RETURN

C     FORMATE:
      END
c---------------------------------------------------------------------
      SUBROUTINE MINMAX(ARRAY,NPTS,BMIN,BMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C       FIND MINIMUM AND MAXIMUM OF THE ARRAY 
C
        DIMENSION ARRAY(NPTS)
        BMIN = ARRAY(1)
        BMAX = BMIN
        DO 100 I=2,NPTS
        BMIN = dMIN1(BMIN,ARRAY(I))
100     BMAX = dMAX1(BMAX,ARRAY(I))
        RETURN
        END
c---------------------------------------------------------------------
      subroutine matcoi ( ia,m,n,ib )
      implicit double precision (a-h,o-z)
      dimension ia(*),ib(*)
      do 100 i=1,m*n
  100 ib(i) = ia(i)
      end
c---------------------------------------------------------------------
      subroutine matco2 ( ia,ib,m,n,ipm)
c...  copy ia(i) into ib(ipm,i)  
      implicit double precision (a-h,o-z)
      dimension ia(n),ib(m,n)
        do i=1,n
          ib(ipm,i) = ia(i)
        end do
      end
c---------------------------------------------------------------------
      subroutine matco3 ( ia,ib,m,n,ipm)
c...  copy ib(ipm,i) into ia(i)  
      implicit double precision (a-h,o-z)
      dimension ia(n),ib(m,n)
        do i=1,n
          ia(i)=ib(ipm,i)
        end do
      end
c----------------------------------------------------------------------
      SUBROUTINE PIVOT (A,LDA,N,B)
c-----------------------------------------------------------------------
c                                                                      *
c     Kommentar SK: genauer als INVERT bei schlechter Kondition!       * 
c                                                                      *
c     Dieses Unterprogramm berechnet die Inverse einer reellen         *
c     quadratischen NxN-Matrix  A  mit dem Austauschverfahren          *
c     (auch Methode des Pivotisierens genannt).                        *
c                                                                      *
c     EINGABEPARAMETER:                                                *
c                                                                      *
c  A       2-dim.  Feld (1:LDA,1:N), das die zu invertierende          *
c          Matrix  A enthält.                                          *
c  LDA     führende Dimension von A, wie im rufenden                   *
c          Programm vereinbart.                                        *
c  N       Ordnung der Matrix A.                                       *
c                                                                      *
c  AUSGABEPARAMETER:                                                   *
c                                                                      *
c  B       2-dim.  Feld (1:LDA,1:N), das die Inverse von A             *
c          enthält.                                                    *
c  S1,S2   Kontrollgrenzen; S1 ist die Summe der Beträge               *
c          der   Diagonalelemente  der  Matrix   A*B-E   (E=Einheits-  *
c          matrix),  S2 ist die Summe der Beträge aller                *
c          übrigen   Elemente;                                         *
c          A*B-E   ist  theoretisch  die  Nullmatrix.                  *
c  IERR    = 1, inverse Matrix von A gefunden.                         *
c          = 2, Matrix A ist singulär, die Inverse existiert           *
c               nicht.                                                 *
c  VAL     Pivotelement, falls A als singulär erkannt wurde.           *
c                                                                      *
c  HILFSFELDER:                                                        *
c                                                                      *
c  MX,MY   1-dim.  INTEGER Felder (1:N).                               *
c                                                                      *
c-----------------------------------------------------------------------
c  Autor       Gisela Engeln-Müllges                                   *
c  Datum       18.05.1987                                              *
c  Quellcode   FORTRAN 77                                              *
c-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)                                  
      DOUBLE PRECISION A(LDA,N),B(LDA,N)                                    
      INTEGER MX(N),MY(N)                                                   
c
c  Berechnung des  Maschinenrundungsfehlers FMACHP.                        
      FMACHP=1.D0
    5 FMACHP=0.5D0*FMACHP
      IF(MASCHD(1.D0+FMACHP).EQ.1)  GOTO  5
      FMACHP=FMACHP*2.D0                                               
c
c  Umspeichern der Matrix  A  auf  B.
      DO 10 I=1,N                                                         
        DO 20 J=1,N
          B(I,J) = A(I,J)                                                 
   20   CONTINUE                                                          
   10 CONTINUE                                                            
c
c  Vorbesetzen der Pivotvektoren MX  und MY mit  Null.                     
      DO 30 I=1,N                                                        
         MX(I)=0                                                        
         MY(I)=0
   30 CONTINUE                                                            
c
c  Bestimmung des  Pivotelementes.                                
      DO 40 I=1,N
         PIVO= 0.0D0
         DO 50 IX=1,N
             IF (MX(IX).EQ.0)  THEN                                     
                DO 60 IY=1,N
                      IF (MY(IY).EQ.0)  THEN
                         IF (DABS(B(IX,IY)).GT.DABS(PIVO))  THEN
                            PIVO=B(IX,IY)                          
                            NX=IX
                            NY=IY                                  
                         END IF                               
                      END IF                               
   60           CONTINUE
             END IF
   50    CONTINUE
c  falls das Pivotelement Null  ist,  ist die Matrix singulär.
csa      IF (DABS(PIVO).LT.1.0d-36)  THEN
         IF (DABS(PIVO).LT. 4.*FMACHP)  THEN
csa            VAL=PIVO                                                   
csa            IERR=2                                                     
csa            RETURN                                                     
             write(*,1000)
             stop
         END  IF              
c   merken der  Indizes  des  Pivotelementes.               
         MX(NX)=NY
         MY(NY)=NX
c   Berechnung  der Matrixelemente gemäß der Rechenregeln für
c   einen Austauschschritt.
         DUMMY=1.D0/PIVO
         DO 70 J=1,N
            IF(J.NE.NX) THEN
               FACTOR=B(J,NY)*DUMMY
               DO 80 K=1,N
                  B(J,K)=B(J,K)-B(NX,K)*FACTOR
   80          CONTINUE
               B(J,NY)=-FACTOR
            END  IF
   70    CONTINUE
         DO 90 K=1,N
            B(NX,K)=B(NX,K)*DUMMY
   90    CONTINUE
         B(NX,NY)=DUMMY
   40 CONTINUE
c
c   Zeilen- und Spaltenvertauschungen rückgängig machen.
      DO 100 I=1,N-1
        DO 110 M=I,N
           IF(MX(M).EQ.I) GOTO 120
  110   CONTINUE
  120   J=M
        IF(J.NE.I) THEN
           DO 130 K=1,N
              H=B(I,K)
              B(I,K)=B(J,K)
              B(J,K)=H
  130      CONTINUE
           MX(J)=MX(I)
           MX(I)=I
        END  IF
        DO 140 M=I,N
           IF(MY(M).EQ.I) GOTO 150
  140   CONTINUE
  150   J=M
        IF(J.NE.I) THEN
           DO 160 K=1,N
              H=B(K,I)
              B(K,I)=B(K,J)
              B(K,J)=H
  160      CONTINUE
           MY(J)=MY(I)
           MY(I)=I
        END  IF
  100 CONTINUE
c
c  Bildung der Differenzmatrix S=(A*B-E), E=Einheitsmatrix.
c  Bildung der Summe S1 der Beträge der Diagonalelemente von S
c  und der Summe S2 der Beträge aller übrigen Glieder.
c  Theoretisch müßten S1 und S2 Null sein.
      S1=0.D0                                                              
      S2=0.D0                                                              
      DO 170 I=1,N                                                         
         DO 180 J=1,N                                                      
              H=0.D0                                                       
              DO 190 K=1,N                                                 
                H=H+A(I,K)*B(K,J)                                          
 190          CONTINUE                                                     
              IF(I.EQ.J) THEN                                              
                S1=S1+DABS(H-1.)                                           
              ELSE                                                         
                S2=S2+DABS(H)                                              
              END IF                                                       
 180     CONTINUE                                                          
 170  CONTINUE                                                              
csk      IERR=1                                                               
1000  format(1x,'stop in SR pivot, diagonal of matrix is zero ' 
     +          'within inversion')
      RETURN                                                               
      END                                                                  
c
      integer function maschd(x)
      double precision x
      maschd=0
      if (1.d0.lt.x) maschd=1
      return
      end
c
      subroutine daxpty(n,t,x,y)
c-------------------------------------------------------------------
c     Purpose:  y = y + t * x 
c-------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension x(*), y(*)
      do k=1,n
        y(k) = y(k) + x(k)*t
      end do    
      return
      end
c
      subroutine daxpty1(n,t,x,ndm,y,ndf,numnp)
c-------------------------------------------------------------------
c     Purpose:  y = y + t * x 
c-------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension x(ndm,numnp), y(ndf,numnp) !ndm<ndf
      if(ndm.gt.ndf) stop 'SR daxpty1' 
      do k=1,numnp
        do i=1,ndm
          y(i,k) = y(i,k) + x(i,k)*t
        end do
      end do    
      return
      end
c
      subroutine skew (w,a)
c-----------------------------------------------------------------------
c           setup a skew-symmetric matrix w from vector a
c                            --               --
c                            |  0   -a3    a2  |
c                            |                 |
c             w = skew a  =  |  a3   0    -a1  |
c                            |                 |
c                            | -a2   a1    0   |
c                            --               --
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension w(3,3), a(3)
c
      w(1,1) =  0.d0
      w(2,1) =  a(3)
      w(3,1) = -a(2)
      w(1,2) = -a(3)
      w(2,2) =  0.d0
      w(3,2) =  a(1)
      w(1,3) =  a(2)
      w(2,3) = -a(1)
      w(3,3) =  0.d0
c
      return
      end
c
      subroutine conden(sc,pc,nseq,neas,nc,flg)
c----------------------------------------------------------------------
c     Condensation of internal DOF's      (see R.D. Cook, page 187)
c
c        (Krr Krc) Dr = Rr
c        (Kcr Kcc) Dc = Rc
c
c        (Krr-Krc Kcc^-1 Kcr) Dr = Rr - Krc Kcc^-1 Rc
c
c     Dim:  Krr(nseq,nseq)  Krc(nseq,neas)  Kcc(neas,neas)
c           Rr (nseq)                       Rc (neas)
c        sc : full matrix
c        pc : full vector
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension sc(nc,nc),pc(nc)
      logical flg
c.... condense matrix and right-hand side vector
c.... lower triangle of sc
      do 10 ieq = 1,neas
        keq = nseq + neas + 1 - ieq
        if(sc(keq,keq).le.0.d0) write(*,*) 'Dii<0 in CONDEN eq. :',keq
        do 20 irow = 1,keq-1
          dum = sc(keq,irow) / sc(keq,keq)
          do 30 jcol = 1,irow
            sc(irow,jcol) = sc(irow,jcol) - sc(keq,jcol) * dum
30        continue
          pc(irow) = pc(irow) - pc(keq) * dum
20      continue
10    continue
      if(flg) return
c.... upper triangle of sc by symmetry
      do 40 k=1,keq-1
        do 40 l=1,k
40        sc(l,k) = sc(k,l)
      return
      end
c
      subroutine decond(skcr,prc,alpha,dv,neas,nseq)
c----------------------------------------------------------------------
c     Decondensation of internal DOF's      (see R.D. Cook, page 187)
c
c        (Krr Krc) Dr = Rr
c        (Kcr Kcc) Dc = Rc
c
c        Dc=-Kcc^-1(Kcr Dr - Rc)
c
c     Kcr(neas,nseq)   Kcc(neas,neas)   --> scr(neas,nseq+neas)
c                      prc(neas)        
c     dv(nseq)         alpha(neas) 
c     dvl(nseq+neas) (only local)   
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension skcr(neas,*), prc(*), dv(*), alpha(*), dvl(nseq+neas)
c     

c.... test for K_cc=0
      tol = 1.d-10
      sum = 0.0d0
      do i = 1,neas
         j = nseq+i
         dsum = skcr(i,j)
          sum = sum + dsum*dsum
      end do
      if(sum.lt.tol)  return

c
      do i = 1,nseq
       dvl(i) = dv(i)
      end do
c
      do ieq = 1,neas
        keq = nseq + ieq
        dum = 0.d0
        do irow = 1,keq-1
          dum = dum + skcr(ieq,irow)*dvl(irow)
        end do ! irow
        if(skcr(ieq,keq).le.0.d0) write(*,*)'D<0 in DECOND',keq
        dvl(keq) = (prc(ieq) - dum) / skcr(ieq,keq)
        alpha(ieq) = alpha(ieq) + dvl(keq)
      end do !ieq
c
      return
      end
c
      subroutine conden1(st,pt,s,p,sc,pc,alpha,vl,nst,ndf,ngeas_max,
     +                   neas,ityp)
c-----------------------------------------------------------------------
c.... static condensation and back substitution for symmetric matrices
c     using INVERT
c     1=decondensation 
c     2=  condensation 
c
c-----------------------------------------------------------------------
      USE eldata
      implicit double precision(a-h,o-z)
      dimension st(ngeas_max,*),pt(*),s(nst,*),p(*),sc(neas,*),pc(*),
     +          alpha(*),vl(*),ph(neas),sh(neas)
c
      nst1 = nel*ndf
c
      if(ityp.eq.1)then
       
c....   back substitution
        do i = 1,neas
          do j = 1,nst1
            pc(i) = pc(i) - sc(i,j)*vl(j)
          end do
        end do
        do i = 1,neas
          do j = 1,neas
            alpha(i) = alpha(i) + sc(i,nst1+j)*pc(j)
          end do
        end do
 
      else if(ityp.eq.2)then
       
c....   store in sc and pc
        do i = 1,neas
          ii = nst1+i
          pc(i) = pt(ii)
          do j = 1,nst1
            sc(i,j) = st(ii,j)
          end do
          do j = 1,i
            jj = nst1+j
            sc(i,jj) = st(ii,jj)
            sc(j,ii) = sc(i,jj)
          end do
        end do
        call invert(sc(1,nst1+1),neas,neas)
       
c.... condensation
        call mvmul (sc(1,nst1+1),pc,neas,neas,ph)
        do i = 1,nst1
          p(i) = pt(i) - ddot(neas,sc(1,i),1,ph,1)
          call mvmul (sc(1,nst1+1),sc(1,i),neas,neas,sh)
          do j = 1,i
            s(i,j) = st(i,j) - ddot(neas,sh,1,sc(1,j),1)
            s(j,i) = s(i,j)
          end do
        end do
 
      end if
c     
      return
      end
c
      subroutine atri(a,neq,meq,det)
c-----------------------------------------------------------------------
c.... dreieckzerlegung  A = LDU' 
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension a(meq,meq)
c.... schleife über zeilen von L und spalten von U -- start mit 2.ter
c.... wegen A(1,1) = d(1,1)
      if(neq.le.1) return
      do 130 j = 2,neq
c....   schleife über spalten von L und zeilen von U
        do 110 i = 1,j-1
          if(i.gt.1) then
            do 100 k = 1,i-1
c....         berechne U(i,j)
              a(i,j) = a(i,j) - a(i,k)*a(k,j)
c....         berechne L(j,i)
              a(j,i) = a(j,i) - a(j,k)*a(k,i)
100           continue
          endif
110     continue
c....   reduziere diagonalterm zu D und berechne U' und L'
        do 120 i = 1,j-1
c....     berechne U'(i,j)
          a(i,j) = a(i,j)/a(i,i)
c....     reduziere diagonale D(j,j)
          a(j,j) = a(j,j) - a(j,i)*a(i,j)
c....     berechne L'(j,i)
          a(j,i) = a(j,i)/a(i,i)
120     continue
130   continue
      det = 1.d0
      do 140 i = 1,neq
140     det = det*a(i,i)
c      print 1000,det
      return
c 1000 format(1x,'Determinante ...............        det = ',e13.4)
      end
c
      subroutine asol(a,b,neq,meq)
c.... gleichungslöser    A x = b  mit  A = LDU' ------------ 
      implicit double precision (a-h,o-z)
      dimension a(meq,meq),b(meq)
      data tol/1.0d-30/
c.... vorwärtseinsetzen für   L'z = b
      if(neq.gt.1) then
        do 110 i = 2,neq
          do 100 k = 1,i-1
            b(i) = b(i) - a(i,k)*b(k)
100       continue
110     continue
      end if
c.... division durch diagonalen zur lsg.  Dy = z
      do 200 i = 1,neq
        if(dabs(a(i,i)).lt.tol) then
          write(*,1000)  i,i,a(i,i),tol
          return
        end if
        b(i) = b(i)/a(i,i)
200   continue
c.... rückwärtseinsetzen durch spalten von U' zur lsg.  U'x = y
      if(neq.gt.1) then
        do 310 i = neq,2,-1
          do 300 k = 1,i-1
            b(k) = b(k) - a(k,i)*b(i)
300       continue
310     continue
      end if
      return
 1000 format(
     +'Gl.system singulär, A(',i3,',',i3,')=',e12.4,'<tol=',e12.4)
      end
c
      subroutine pmultst (s,t,ndim,nr,nc,nt,ir,jc,itp)
c-----------------------------------------------------------------------
c.... multiply column/row --bloc with transformation matrix T(nt,nt)
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension sh(ndim,ndim),t(nt,nt),s(ndim,*)
c
      do i = 1,nr
       do j = 1,nc
        sh(i,j) = s(ir+i,jc+j)
       enddo
      enddo 
c
      if (itp.eq.1) then       ! T^T S 
        do i = 1,nt
         do j = 1,nc
          s(ir+i,jc+j) = 0.d0 
          do k = 1,nt
           s(ir+i,jc+j) = s(ir+i,jc+j) + t(k,i)*sh(k,j)
          enddo 
         enddo
        enddo
      elseif (itp.eq.2) then        ! S T 
        do i = 1,nr
         do j = 1,nt
          s(ir+i,jc+j) = 0.d0 
          do k = 1,nt
           s(ir+i,jc+j) = s(ir+i,jc+j) + sh(i,k)*t(k,j)
          enddo 
         enddo
        enddo
      endif  
c
      return
      end

