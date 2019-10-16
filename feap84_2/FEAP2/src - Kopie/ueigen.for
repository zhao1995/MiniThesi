C
      SUBROUTINE UEIGEN(MAT,WERTR,WERTI,EIVEC,SKAL,CNT,ND)
      USE iofile
      IMPLICIT	DOUBLE PRECISION (A-H,O-Z)
      INTEGER CNT(*),EIGEN,BASIS,RES,LOW,HIGH,IMA
      DOUBLE PRECISION MAT
      DIMENSION MAT(ND,*),WERTR(*),WERTI(1),EIVEC(ND,*),SKAL(*)

C..................................................................
      DATA BASIS/2/
      M = ND
C..................................................................
C     Berechnung der Eigenwerte und Eigenvektoren
C...................................................................
      RES = EIGEN(BASIS,ND,M,MAT,SKAL,EIVEC,WERTR,WERTI,CNT,LOW,HIGH)
C...................................................................
C
      IF (RES .GT. 0) THEN
	RES = RES-400
C...................................................................
C     AUSGABE EINER ABBRUCHMELDUNG GEMAESS RES (0/401/402/403)
C...................................................................
	GOTO (30,40,50),RES
   30	STOP 'FEHLER: DIE ORDNUNG DER MATRIX IST KLEINER ALS 1 !'
   40	STOP 'ABBRUCH: DIE MATRIX IST DIE NULLMATRIX !'
   50	continue
	write(*,2001)
	return
      ELSE
C..................................................................
C     AUSGABE DER LOESUNG
C..................................................................
	IMA	  = 0
	DO 60	I = 1,M
cww im1 no init.value  IF ( DABS(WERTI(I)).GT.1.D-12 ) IMA = IM1 + 2
   60	CONTINUE
c	 IF ( IMA.LT.1 )  THEN
C..................................................................
C     NORMIERUNG DER REELLEN EIGENVEKTOREN
C..................................................................
c	   DO 70   I = 1,M
c	     BETRAG   = 0.D0
c	     DO 75   J = 1,M
c   75	     BETRAG    = BETRAG    + EIVEC(I,J) * EIVEC(I,J)
c	     BETRAG    = SQRT(BETRAG)
c	     DO 76   J = 1,M
c   76	     EIVEC(I,J)=EIVEC(I,J)/BETRAG
c   70	   CONTINUE
c	 ENDIF
C..................................................................
	WRITE(iow,2015)
	if(ior.le.0) WRITE(*,2015)
	WRITE(iow,2020)(WERTR(I),WERTI(I),I=1,M)
	if(ior.le.0) WRITE(*,2020)(WERTR(I),WERTI(I),I=1,M)
c	 WRITE(iow,2030)
c	 if(ior.le.0) WRITE(*,2030)
c	 DO 20 I=1,M
c	   WRITE(iow,2040)(EIVEC(I,J),J=1,M)
c	   if(ior.le.0) WRITE(*,2040)(EIVEC(I,J),J=1,M)
c20	 CONTINUE
      ENDIF
      RETURN
C
c 1900 FORMAT(/,1X,'ELEMENT:',I3,'  GAUSSPUNKT:',I1)
c 2000 FORMAT(14H			 ,/,
c     +       1X,'ZU BERECHNENDE MATRIX:',/)
 2001 format(1x,'ABBRUCH: MAX.SCHRITTZAHL FUER QR-VERFAHREN',
     1	    ' UEBERSCHRITTEN!')
c 2005 FORMAT('EL:',I3,'  GP:',I2,'-NICHTPOSITIVER EIGENWERT: ',D15.5)
c 2010 FORMAT(6(1X,D12.4))
 2015 FORMAT(/,1X,'-BERECHNETE EIGENWERTE:    REALTEIL    I  ',
     +	     'IMAGINAERTEIL',/,25X,31(1H-))
 2020 FORMAT(24X,D15.5,' I',D15.5)
c 2025 FORMAT(/,1X,'-ELEMENT:', I3,'   GP:', I1,'   DET(AA) :',D15.5)
c 2030 FORMAT(/,1X,'-MATRIX DER NORMALISIERTEN EIGENVEKTOREN:',/)
c 2040 FORMAT(6(1X,D12.4))
c 2051 FORMAT('VOLLST. STOFFMATRIX ')
c 2052 FORMAT('ANTIM.  STOFFMATRIX ')
c 2053 FORMAT('SYM.    STOFFMATRIX ')
      END
c
      INTEGER FUNCTION EIGEN
     *	 (BASIS,ND,N,MAT,SKAL,EIVEC,WERTR,WERTI,CNT,LOW,HIGH)
C
C*****************************************************************
C								 *
C     ZWECK DES PROGRAMMS:					 *
C     ====================					 *
C     DIESES FUNKTIONSUNTERPROGRAMM VOM TYP INTEGER BERECHNET	 *
C     ALLE EIGENWERTE UND EIGENVEKTOREN EINER BELIEBIGEN	 *
C     REELLEN MATRIX MAT.					 *
C     DIE EIGENWERTE WERDEN IN DEN FELDERN WERTR(1:N)		 *
C     (REALTEIL) UND WERTI(1:N) (IMAGINAERTEIL) UND DIE 	 *
C     EIGENVEKTOREN IM FELD EIVEC(1:N,1:N) ABGESPEICHERT.	 *
C								 *
C     EINGABEPARAMETER: 					 *
C     ================= 					 *
C     N:	   ORDNUNG DER MATRIX MAT			 *
C     ND:	   FUEHRENDE DIMENSION DER FELDER MAT UND	 *
C		   EIVEC, WIE IM RUFENDEN PROGRAMM VEREINBART	 *
C     MAT:	   EIN (N,N)-FELD VOM TYP DOUBLE PRECISION,	 *
C		   DAS DIE MATRIX ENTHAELT, DEREN EIGENWERTE	 *
C		   UND EIGENVEKTOREN BERECHNET WERDEN SOLLEN	 *
C     BASIS:	   DIE BASIS DER GLEITKOMMADARSTELLUNG DES	 *
C		   VERWENDETEN RECHNERS (MEISTENS 2 ODER 16)	 *
C								 *
C     AUSGABEPARAMETER: 					 *
C     ================= 					 *
C     MAT:	   DER OBERE TEIL DIESES (N,N)-FELDES ENTHAELT	 *
C		   DIE EIGENVEKTOREN DER QUASI-DREIECKSMATRIX,	 *
C		   DIE VOM QR-VERFAHREN ERZEUGT WIRD.		 *
C     WERTR,WERTI: ZWEI (N,1)-FELDER VOM TYP DOUBLE		 *
C		   PRECISION, DIE REALTEIL UND IMAGINAER-	 *
C		   TEIL DER EIGENWERTE AUFNEHMEN		 *
C     EIVEC:	   EIN (N,N)-FELD VOM TYP DOUBLE PRECISION,	 *
C		   DAS DIE NORMALISIERTEN EIGENVEKTOREN DER	 *
C		   URSPRUENGLICHEN VOLLEN MATRIX AUFNIMMT.	 *
C		   FALLS DER I-TE EIGENWERT REELL IST, DANN	 *
C		   IST DIE I-TE SPALTE VON EIVEC DER DAZUGE-	 *
C		   HOERIGE REELLE EIGENVEKTOR. FALLS DIE	 *
C		   EIGENWERTE I UND I+1 EIN KOMPLEXES		 *
C		   PAAR BILDEN, GEBEN I-TE UND (I+1)-TE 	 *
C		   SPALTE REAL- UND IMAGINAERTEIL DESJENIGEN	 *
C		   EIGENVEKTORS AN, DER ZU DEM EIGENWERT MIT	 *
C		   POSITIVEM IMAGINAERTEIL GEHOERT.		 *
C     CNT:	   EIN (N,1)-FELD VOM TYP INTEGER, DAS DIE	 *
C		   ZAHL DER ITERATIONSSCHRITTE FUER JEDEN	 *
C		   EIGENWERT AUFNIMMT. FALLS ZWEI EIGENWERTE	 *
C		   ALS PAAR GLEICHZEITIG GEFUNDEN WERDEN, WIRD	 *
C		   DIE ZAHL DER ITERATIONSSCHRITTE MIT EINEM	 *
C		   POSITIVEN VORZEICHEN FUER DEN ERSTEN UND	 *
C		   EINEM NEGATIVEN VORZEICHEN FUER DEN		 *
C		   ZWEITEN EIGENWERT EINGETRAGEN.		 *
C     SKAL:	   EIN (N,1)-FELD VOM TYP DOUBLE PRECISION, DAS  *
C		   DIE INFORMATION UEBER DIE DURCHGEFUEHRTEN	 *
C		   VERTAUSCHUNGEN UND DIE SKALIERUNGSFAKTOREN	 *
C		   ENTHAELT.					 *
C								 *
C     RUECKGABEWERTE:						 *
C     ===============						 *
C     0:      KEIN FEHLER					 *
C     401:    DIE ORDNUNG N DER MATRIX MAT IST KLEINER ALS 1.	 *
C     402:    MAT IST DIE NULLMATRIX.				 *
C     403:    DIE MAXIMALE SCHRITTZAHL VON FUER DAS		 *
C	      QR-VERFAHREN IST UEBER- SCHRITTEN, OHNE DASS	 *
C	      ALLE EIGENWERTE BERECHNET WERDEN KONNTEN. 	 *
C								 *
C     LOKALE GROESSEN:						 *
C     ================						 *
C     ONE,TWO,HALF: GLEITKOMMAKONSTANTEN			 *
C     EPS:	    MASCHINENGENAUIGKEIT			 *
C     TEMP:	    HILFSVARIABLE				 *
C								 *
C----------------------------------------------------------------*
C								 *
C  BENOETIGTE UNTERPROGRAMME: BALAN, ELMHES, ELMTRA, HQR2,	 *
C			      BALBAK, NORMAL, SWAP, COMDIV,	 *
C			      COMABS				 *
C								 *
C								 *
C  QUELLEN : 1. MARTIN, R. S. UND WILKINSON, J. H.,		 *
C		SIEHE [MART70]. 				 *
C	     2. PARLETT, B. N. UND REINSCH, C., SIEHE [PARL69].  *
C	     3. PETERS, G. UND WILKINSON, J. H., SIEHE [PETE70]. *
C								 *
C*****************************************************************
C
      INTEGER BASIS,ND,N,CNT(N),LOW,HIGH
      DOUBLE PRECISION MAT(ND,N),SKAL(N),EIVEC(ND,N),WERTR(N),
     *		       WERTI(N)
      DOUBLE PRECISION ONE,TWO,HALF
      PARAMETER (ONE = 1.0D0,TWO = 2.0D0,HALF = 0.5D0)
      INTEGER RES,BALAN,ELMHES,ELMTRA,HQR2,BALBAK,NORMAL
      DOUBLE PRECISION EPS,TEMP
C
C     BERECHNUNG DER MASCHINENGENAUIGKEIT EPS (D. H. DER KLEINSTEN
C     POSITIVEN MASCHINENZAHL, FUER DIE AUF DEM RECHNER GILT:
C     1 + EPS > 1):
C
      TEMP = TWO
      EPS = ONE
   10 IF (ONE .LT. TEMP) THEN
	  EPS = EPS * HALF
	  TEMP = ONE + EPS
	  GOTO 10
	  ENDIF
      EPS = TWO * EPS
      RES = BALAN(ND,N,MAT,SKAL,LOW,HIGH,BASIS)
	 IF (RES .NE. 0) THEN
	    EIGEN = RES + 100
	    RETURN
	 ENDIF
      RES = ELMHES(ND,N,LOW,HIGH,MAT,CNT)
	 IF (RES .NE. 0) THEN
	    EIGEN = RES + 200
	    RETURN
	 ENDIF
      RES = ELMTRA(ND,N,LOW,HIGH,MAT,CNT,EIVEC)
	 IF (RES .NE. 0) THEN
	    EIGEN = RES + 300
	    RETURN
	 ENDIF
      RES = HQR2(ND,N,LOW,HIGH,MAT,WERTR,WERTI,EIVEC,CNT,EPS)
	 IF (RES .NE. 0) THEN
	    EIGEN = RES + 400
	    RETURN
	 ENDIF
      RES = BALBAK(ND,N,LOW,HIGH,SKAL,EIVEC)
	 IF (RES .NE. 0) THEN
	    EIGEN = RES + 500
	    RETURN
	 ENDIF
      RES = NORMAL(ND,N,EIVEC,WERTI)
	 IF (RES .NE. 0) THEN
	    EIGEN = RES + 600
	    RETURN
	 ENDIF
      EIGEN = 0
      RETURN
      END
C
C
      INTEGER FUNCTION BALAN (ND,N,MAT,SKAL,LOW,HIGH,BASIS)
C
C*****************************************************************
C								 *
C     ZWECK DES PROGRAMMS:					 *
C     ====================					 *
C     DIE PROZEDUR BALAN BALANCIERT EINE GEGEBENE REELLE MATRIX  *
C     IN DER 1-NORM AUS.					 *
C								 *
C     EINGABEPARAMETER: 					 *
C     ================= 					 *
C     ND:	FUEHRENDE DIMENSION DER MATRIZEN, WIE SIE IM	 *
C		HAUPTPROGRAMM VEREINBART WURDEN 		 *
C     N:	DIE ORDNUNG DER GEGEBENEN MATRIX		 *
C     MAT:	EIN (1:N,1:N)-FELD, DAS DIE KOMPONENTEN DER	 *
C		GEGEBENEN MATRIX ENTHAELT			 *
C     BASIS:	DIE BASIS DER GLEITKOMMADARSTELLUNG AUF 	 *
C		DER MASCHINE					 *
C								 *
C     AUSGABEPARAMETER: 					 *
C     ================= 					 *
C     MAT:	DIE AUSBALANCIERTE MATRIX			 *
C     LOW,HIGH: ZWEI INTEGERZAHLEN, FUER DIE MAT(I,J)		 *
C		GLEICH NULL IST, FALLS GILT:			 *
C		1. I>J UND					 *
C		2. J=1,...LOW-1 ODER I=HIGH+1,...N		 *
C     SKAL:	EIN (1:N)-FELD, DAS DIE INFORMATIONEN UEBER	 *
C		DIE DURCHGEFUEHRTEN VERTAUSCHUNGEN UND DIE	 *
C		SKALIERUNGSFAKTOREN ENTHAELT.			 *
C								 *
C     RUECKGABEWERT:						 *
C     ==============						 *
C     0:	     KEIN FEHLER				 *
C								 *
C     LOKALE GROESSEN:						 *
C     ================						 *
C     ZERO,ONE,PT95: GLEITKOMMAKONSTANTEN			 *
C     I,J,K,L:	     ZAEHLVARIABLEN				 *
C     B2:	     QUADRAT VON BASIS				 *
C     R,C,F,G,S:     HILFSVARIABLEN ZUR BERECHNUNG VON		 *
C		     ZEILENNORMEN, IHREN KEHRWERTEN UND 	 *
C		     AEHNLICHEM 				 *
C								 *
C----------------------------------------------------------------*
C								 *
C  BENOETIGTE UNTERPROGRAMME: SWAP				 *
C								 *
C								 *
C  QUELLEN : PARLETT, B. N. UND REINSCH, C., SIEHE [PARL69].	 *
C								 *
C*****************************************************************
C								 *
C  AUTOR     : JUERGEN DIETEL					 *
C  DATUM     : 10.04.1987					 *
C  QUELLCODE : FORTRAN 77					 *
C								 *
C*****************************************************************
C
      INTEGER ND,N,LOW,HIGH,BASIS
      DOUBLE PRECISION SKAL(N),MAT(ND,N)
      DOUBLE PRECISION ZERO,ONE,PT95
      PARAMETER (ZERO = 0.0D0,ONE = 1.0D0,PT95 = 0.95D0)
      INTEGER I,J,K,L,B2
      DOUBLE PRECISION R,C,F,G,S
C
C     DIE NORM VON MAT(1:N,1:N) REDUZIEREN DURCH EXAKTE AEHNLICH-
C     KEITSTRANSFORMATIONEN, DIE IN SKAL(1:N) ABGESPEICHERT WERDEN
C
      B2 = BASIS*BASIS
      L = 1
      K = N
C
C     NACH ZEILEN MIT EINEM ISOLIERTEN EIGENWERT SUCHEN UND SIE
C     NACH UNTEN SCHIEBEN
C
   10 DO 50 J=K,1,-1
	 R = ZERO
	 DO 20 I=1,K
   20	    IF (I .NE. J) R = R+ABS(MAT(J,I))
	 IF (R .EQ. ZERO) THEN
	    SKAL(K) = J
	    IF (J .NE. K) THEN
	       DO 30 I=1,K
   30		  CALL SWAP(MAT(I,J),MAT(I,K))
	       DO 40 I=L,N
   40		  CALL SWAP(MAT(J,I),MAT(K,I))
	    ENDIF
	    K = K-1
	    GOTO 10
	 ENDIF
   50	 CONTINUE
C
C     NACH SPALTEN MIT EINEM ISOLIERTEN EIGENWERT SUCHEN UND SIE
C     NACH LINKS SCHIEBEN
C
   60 DO 100 J=L,K
	 C = ZERO
	 DO 70 I=L,K
   70	    IF (I .NE. J) C = C+ABS(MAT(I,J))
	 IF (C .EQ. ZERO) THEN
	    SKAL(L) = J
	    IF (J .NE. L) THEN
	       DO 80 I=1,K
   80		  CALL SWAP(MAT(I,J),MAT(I,L))
	       DO 90 I=L,N
   90		  CALL SWAP(MAT(J,I),MAT(L,I))
	    ENDIF
	    L = L+1
	    GOTO 60
	 ENDIF
  100	 CONTINUE
C
C     NUN DIE TEILMATRIX IN DEN ZEILEN L BIS K AUSBALANCIEREN
C
      LOW = L
      HIGH = K
      DO 110 I=L,K
  110	 SKAL(I) = ONE
  120 DO 180 I=L,K
	 C = ZERO
	 R = ZERO
	 DO 130 J=L,K
	    IF (J .NE. I) THEN
	       C = C+ABS(MAT(J,I))
	       R = R+ABS(MAT(I,J))
	    ENDIF
  130	    CONTINUE
	 G = R/BASIS
	 F = ONE
	 S = C+R
  140	 IF (C .LT. G) THEN
	    F = F*BASIS
	    C = C*B2
	    GOTO 140
	 ENDIF
	 G = R*BASIS
  150	 IF (C .GE. G) THEN
	    F = F/BASIS
	    C = C/B2
	    GOTO 150
	 ENDIF
	 IF ((C+R)/F .LT. PT95*S) THEN
	    G = ONE/F
	    SKAL(I) = SKAL(I)*F
	    DO 160 J=L,N
  160	       MAT(I,J) = MAT(I,J)*G
	    DO 170 J=1,K
  170	       MAT(J,I) = MAT(J,I)*F
	    GOTO 120
	 ENDIF
  180	 CONTINUE
      BALAN = 0
      END
C
C
      INTEGER FUNCTION BALBAK(ND,N,LOW,HIGH,SKAL,EIVEC)
C
C*****************************************************************
C								 *
C     ZWECK DES PROGRAMMS:					 *
C     ====================					 *
C     BALBAK FUEHRT EINE RUECKTRANSFORMATION ALLER RECHTSEIGEN-  *
C     VEKTOREN EINER AUSBALANCIERTEN MATRIX IN DIE EIGENVEKTOREN *
C     DER ORIGINALMATRIX DURCH, VON DER DIE BALANCIERTE MATRIX	 *
C     DURCH AUFRUF DER PROZEDUR BALAN ABGELEITET WURDE. 	 *
C								 *
C     EINGABEPARAMETER: 					 *
C     ================= 					 *
C     ND:	FUEHRENDE DIMENSION DER MATRIZEN, WIE SIE IM	 *
C		HAUPTPROGRAMM VEREINBART WURDEN 		 *
C     N:	DIE ORDNUNG DER EIGENVEKTOREN (ZAHL DER 	 *
C		KOMPONENTEN)					 *
C     LOW,HIGH: ZWEI INTEGERZAHLEN, DIE VON DER PROZEDUR	 *
C		BALAN STAMMEN					 *
C     SKAL:	AUSGABEVEKTOR DER PROZEDUR BALAN		 *
C     EIVEC:	EIN (1:N,1:N)-FELD, VON DEM JEDE SPALTE EINEN	 *
C		EIGENVEKTOR (ODER SEINEN REALTEIL ODER SEINEN	 *
C		IMAGINAERTEIL) DER AUSBALANCIERTEN MATRIX	 *
C		DARSTELLT					 *
C								 *
C     AUSGABEPARAMETER: 					 *
C     ================= 					 *
C     EIVEC:	DIE ENTSPRECHENDEN EIGENVEKTOREN (ODER REAL-	 *
C		TEILE ODER IMAGINAERTEILE) DER			 *
C		URSPRUENGLICHEN MATRIX				 *
C								 *
C     RUECKGABEWERT:						 *
C     ==============						 *
C     0:	KEIN FEHLER					 *
C								 *
C     LOKALE VARIABLEN: 					 *
C     ================= 					 *
C     I,J,K:	HILFSVARIABLEN ZUR INDEXBILDUNG 		 *
C     S:   :	SKALIERUNGSWERT 				 *
C								 *
C----------------------------------------------------------------*
C								 *
C  BENOETIGTE UNTERPROGRAMME: SWAP				 *
C								 *
C								 *
C  QUELLEN : PARLETT, B. N. UND REINSCH, C., SIEHE [PARL69].	 *
C								 *
C*****************************************************************
C								 *
C  AUTOR     : JUERGEN DIETEL					 *
C  DATUM     : 10.04.1987					 *
C  QUELLCODE : FORTRAN 77					 *
C								 *
C*****************************************************************
C
      INTEGER ND,N,LOW,HIGH
      DOUBLE PRECISION SKAL(N),EIVEC(ND,N)
      INTEGER I,J,K
      DOUBLE PRECISION S
C
      DO 20 I=LOW,HIGH
	 S = SKAL(I)
C
C	 LINKSEIGENVEKTOREN WERDEN ZURUECKTRANSFORMIERT, INDEM MAN
C	 DIE VORIGE ANWEISUNG ERSETZT DURCH: 'S = 1.0D0/SKAL(I)'
C
	 DO 10 J=1,N
   10	    EIVEC(I,J) = EIVEC(I,J)*S
   20	 CONTINUE
      DO 40 I=LOW-1,1,-1
	 K=SKAL(I)
	 DO 30 J=1,N
   30	    CALL SWAP(EIVEC(I,J),EIVEC(K,J))
   40	 CONTINUE
      DO 60 I=HIGH+1,N
	 K=SKAL(I)
	 DO 50 J=1,N
   50	    CALL SWAP(EIVEC(I,J),EIVEC(K,J))
   60	 CONTINUE
      BALBAK = 0
      END
C
C
      INTEGER FUNCTION ELMHES(ND,N,LOW,HIGH,MAT,PERM)
C
C*****************************************************************
C								 *
C     ZWECK DES PROGRAMMS:					 *
C     ====================					 *
C     GEGEBEN IST EINE UNSYMMETRISCHE MATRIX A(1:N,1:N). DANN	 *
C     REDUZIERT DIESE PROZEDUR DIE TEILMATRIX DER ORDNUNG	 *
C     HIGH-LOW+1, DIE BEIM ELEMENT A(LOW,LOW) BEGINNT UND BEIM	 *
C     ELEMENT A(HIGH,HIGH) ENDET, AUF HESSENBERGFORM H DURCH	 *
C     NICHTORTHOGONALE ELEMENTARTRANSFORMATIONEN. DIE TEILMATRIX *
C     WIRD MIT H UEBERSCHRIEBEN, WOBEI DIE EINZELHEITEN DER	 *
C     TRANSFORMATIONEN IN DEM UEBRIGBLEIBENDEN DREIECK UNTERHALB *
C     VON H UND IN DEM FELD PERM ABGESPEICHERT WERDEN.		 *
C								 *
C     EINGABEPARAMETER: 					 *
C     ================= 					 *
C     ND:	FUEHRENDE DIMENSION DER MATRIZEN, WIE SIE IM	 *
C		HAUPTPROGRAMM VEREINBART WURDEN 		 *
C     N:	ORDNUNG DER VOLLEN MATRIX A			 *
C     LOW,HIGH: AUSGABEPARAMETER EINER PROZEDUR, DIE		 *
C		A AUFBEREITET (SIEHE 'PARLETT, B. N., AND C.     *
C		REINSCH: BALANCING A MATRIX FOR CALCULATION	 *
C		OF EIGENVALUES AND EIGENVECTORS. NUMERISCHE	 *
C		MATHEMATIK 13 (1969), SEITEN 293 - 304, 	 *
C		PROZEDUR BALANCE'). FALLS A NICHT DERART         *
C		AUFBEREITET IST, SETZE LOW:=1, HIGH:=N. 	 *
C     MAT:	DIE (N,N)-MATRIX A, NORMALERWEISE IN		 *
C		AUFBEREITETER FORM (SIEHE OBEN) 		 *
C								 *
C     AUSGABEPARAMETER: 					 *
C     ================= 					 *
C     MAT:	EIN (N,N)-FELD, DAS ZU EINEM TEIL AUS DER	 *
C		ABGELEITETEN OBEREN HESSENBERGMATRIX BESTEHT;	 *
C		DIE GROESSE N(I,R+1), DIE BEI DER REDUKTION	 *
C		EINE ROLLE SPIELT, WIRD IM (I,R)-ELEMENT	 *
C		ABGESPEICHERT.					 *
C     PERM:	EIN INTEGERFELD, DAS DIE BEI DER REDUKTION	 *
C		AUSGEFUEHRTEN ZEILEN- UND SPALTENVERTAU-	 *
C		SCHUNGEN BESCHREIBT				 *
C								 *
C     RUECKGABEWERT:						 *
C     ==============						 *
C     0:	KEIN FEHLER					 *
C								 *
C     LOKALE GROESSEN:						 *
C     ================						 *
C     ZERO,ONE: GLEITKOMMAKONSTANTEN				 *
C     I,J,M:	ZAEHLVARIABLEN					 *
C     X,Y:	HILFSVARIABLEN ZUR AUFNAHME VON MATRIXELEMENTEN  *
C		UND ZWISCHENERGEBNISSEN 			 *
C								 *
C----------------------------------------------------------------*
C								 *
C  BENOETIGTE UNTERPROGRAMME: SWAP				 *
C								 *
C								 *
C  QUELLEN : MARTIN, R. S. UND WILKINSON, J. H., SIEHE [MART70]. *
C								 *
C*****************************************************************
C								 *
C  AUTOR     : JUERGEN DIETEL					 *
C  DATUM     : 10.04.1987					 *
C  QUELLCODE : FORTRAN 77					 *
C								 *
C*****************************************************************
C
      INTEGER N,LOW,HIGH,PERM(N)
      DOUBLE PRECISION MAT(ND,N)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO = 0.0D0)
      INTEGER I,J,M
      DOUBLE PRECISION X,Y
C
      DO 70 M=LOW+1,HIGH-1
	 I = M
	 X = ZERO
	 DO 10 J=M,HIGH
	    IF (ABS(MAT(J,M-1)) .GT. ABS(X)) THEN
	       X = MAT(J,M-1)
	       I=J
	    ENDIF
   10	    CONTINUE
	 PERM(M) = I
	 IF (I .NE. M) THEN
C
C	    ZEILEN UND SPALTEN VON MAT VERTAUSCHEN
C
	    DO 20 J=M-1,N
   20	       CALL SWAP (MAT(I,J),MAT(M,J))
	    DO 30 J=1,HIGH
   30	       CALL SWAP (MAT(J,I),MAT(J,M))
	 ENDIF
	 IF (X .NE. ZERO) THEN
	    DO 60 I=M+1,HIGH
	       Y = MAT(I,M-1)
	       IF (Y .NE. ZERO) THEN
		  Y = Y/X
		  MAT(I,M-1) = Y
		  DO 40 J=M,N
   40		     MAT(I,J) = MAT(I,J)-Y*MAT(M,J)
		  DO 50 J=1,HIGH
   50		     MAT(J,M) = MAT(J,M)+Y*MAT(J,I)
	       ENDIF
   60	       CONTINUE
	 ENDIF
   70	 CONTINUE
      ELMHES = 0
      END
C
C
      INTEGER FUNCTION ELMTRA(ND,N,LOW,HIGH,MAT,PERM,H)
C
C*****************************************************************
C								 *
C     ZWECK DES PROGRAMMS:					 *
C     ====================					 *
C     DIEJENIGE MATRIX IN DEM FELD H(1:N,1:N) ABSPEICHERN, DIE	 *
C     SICH AUS DEN INFORMATIONEN ERGIBT, DIE DIE PROZEDUR ELMHES *
C     HINTERLASSEN HAT IM UNTEREN DREIECK DER HESSENBERGMATRIX	 *
C     H, UND ZWAR IM FELD MAT(1:N,1:N) UND IM INTEGERFELD	 *
C     PERM(1:N) 						 *
C								 *
C     EINGABEPARAMETER: 					 *
C     ================= 					 *
C     ND:	FUEHRENDE DIMENSION DER MATRIZEN, WIE SIE IM	 *
C		HAUPTPROGRAMM VEREINBART WURDEN 		 *
C     N:	ORDNUNG DER HESSENBERGMATRIX H			 *
C     LOW,HIGH: INTEGERZAHLEN, DIE VON DER PROZEDUR BALAN	 *
C		ERZEUGT WURDEN (FALLS SIE VERWANDT WURDE;	 *
C		ANDERNFALLS SETZE LOW:=1, HIGH:=N.)		 *
C     PERM:	EIN VON ELMHES ERZEUGTES (N,1)-INTEGERFELD	 *
C     MAT:	EIN (N,N)-FELD, DAS VON ELMHES ERZEUGT WURDE	 *
C		UND DIE HESSENBERGMATRIX H UND DIE		 *
C		MULTIPLIKATOREN ENTHAELT, DIE BENUTZT WURDEN,	 *
C		UM ES AUS DER ALLGEMEINEN MATRIX		 *
C		A ZU ERZEUGEN					 *
C								 *
C     AUSGABEPARAMETER: 					 *
C     ================= 					 *
C     H:	DASJENIGE (N,N)-FELD, DAS DIE AEHNLICHKEITS-	 *
C		TRANSFORMATION VON A IN H DEFINIERT		 *
C								 *
C     RUECKGABEWERT:						 *
C     ==============						 *
C     0:	     KEIN FEHLER				 *
C								 *
C     LOKALE GROESSEN:						 *
C     ================						 *
C     ZERO,ONE: GLEITKOMMAKONSTANTEN				 *
C     I,J,K:	INDEXVARIABLEN					 *
C								 *
C----------------------------------------------------------------*
C								 *
C  BENOETIGTE UNTERPROGRAMME: KEINE				 *
C								 *
C								 *
C  QUELLEN : PETERS, G. UND WILKINSON, J. H., SIEHE [PETE70].	 *
C								 *
C*****************************************************************
C								 *
C  AUTOR     : JUERGEN DIETEL					 *
C  DATUM     : 10.04.1987					 *
C  QUELLCODE : FORTRAN 77					 *
C								 *
C*****************************************************************
C

      INTEGER ND,N,LOW,HIGH,PERM(N)
      DOUBLE PRECISION MAT(ND,N),H(ND,N)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO = 0.0D0,ONE = 1.0D0)
      INTEGER I,J,K
      DO 20 I=1,N
	 DO 10 J=1,N
   10	    H(I,J) = ZERO
	 H(I,I) = ONE
   20	 CONTINUE
C
      DO 50 I=HIGH-1,LOW+1,-1
	 J=PERM(I)
	 DO 30 K=I+1,HIGH
   30	    H(K,I) = MAT(K,I-1)
	 IF (I .NE. J) THEN
	    DO 40 K=I,HIGH
	       H(I,K)=H(J,K)
	       H(J,K)=ZERO
   40	       CONTINUE
	    H(J,I) = ONE
	 ENDIF
   50	 CONTINUE
      ELMTRA = 0
      END
C
C
      INTEGER FUNCTION HQR2 (ND,N,LOW,HIGH,H,WERTR,WERTI,
     *			     EIVEC,CNT,EPS)
C
C*****************************************************************
C								 *
C     ZWECK DES PROGRAMMS:					 *
C     ====================					 *
C     FINDET DIE EIGENWERTE UND EIGENVEKTOREN EINER REELLEN	 *
C     MATRIX, DIE, AUF OBERE HESSENBERGFORM REDUZIERT, IM FELD	 *
C     H(1:N,1:N) STEHT, WOBEI DAS PRODUKT DER BISHER DURCHGE-	 *
C     FUEHRTEN TRANSFORMATIONEN IM FELD EIVEC(1:N,1:N) STEHT.	 *
C     DIE REAL- UND DIE IMAGINAERTEILE DER EIGENWERTE WERDEN IN  *
C     DEN FELDERN WERTR(1:N),WERTI(1:N) UND DIE EIGENVEKTOREN IM *
C     FELD EIVEC(1:N,1:N) GEBILDET, WO NUR EIN KOMPLEXER VEKTOR, *
C     NAEMLICH DER ZU DER WURZEL MIT POSITIVEM IMAGINAERTEIL GE- *
C     HOERIGE, FUER JEDES KOMPLEXE PAAR VON EIGENWERTEN ERZEUGT  *
C     WIRD. LOW UND HIGH SIND ZWEI INTEGERZAHLEN, DIE BEIM AUS-  *
C     BALANCIEREN ENTSTEHEN, WO EIGENWERTE IN DEN POSITIONEN 1	 *
C     BIS LOW-1 UND HIGH+1 BIS N ISOLIERT WERDEN.  FALLS KEINE	 *
C     AUSBALANCIERUNG DURCHGEFUEHRT WURDE, SETZE LOW:=1,HIGH:=N. *
C     DAS UNTERPROGRAMM BRICHT MIT EINER FEHLERMELDUNG AB, FALLS *
C     IRGENDEIN EIGENWERT MEHR ALS MAXSTP ITERATIONSSCHRITTE	 *
C     BENOETIGT.						 *
C								 *
C     EINGABEPARAMETER: 					 *
C     ================= 					 *
C     N:	ORDNUNG DER HESSENBERGMATRIX H			 *
C     LOW,HIGH: VON BALAN ERZEUGTE INTEGERZAHLEN, FALLS 	 *
C		BALAN BENUTZT WURDE. ANSONSTEN SETZE		 *
C		LOW:=1, HIGH:=N.				 *
C     EPS:	DIE KLEINSTE ZAHL AUF DEM COMPUTER, FUER	 *
C		DIE GILT: 1 + EPS > 1.				 *
C     H:	EIN (N,N)-FELD, DAS DIE MATRIX H IN IHREN	 *
C		RELEVANTEN TEILEN ENTHAELT			 *
C     EIVEC:	EIN (N,N)-FELD, DAS DIE MATRIX ENTHAELT, DIE	 *
C		DIE AEHNLICHKEITSTRANSFORMATION VON A IN H	 *
C		DEFINIERT. (ES WIRD VON ELMTRA ERZEUGT.)	 *
C		FALLS H DIE URSPRUENGLICHE MATRIX IST,		 *
C		SETZE EIVEC := EINHEITSMATRIX.			 *
C								 *
C     AUSGABEPARAMETER: 					 *
C     ================= 					 *
C     H:	   DER OBERE TEIL DIESES (N,N)-FELDES ENT-	 *
C		   HAELT DIE EIGENVEKTOREN DER QUASI-		 *
C		   DREIECKSMATRIX, DIE VOM QR-VERFAHREN 	 *
C		   ERZEUGT WIRD.				 *
C     WERTR,WERTI: ZWEI (N,1)-FELDER, DIE REAL- UND IMAGI-	 *
C		   NAERTEIL DER EIGENWERTE AUFNEHMEN		 *
C     CNT:	   EIN (N,1)-INTEGERFELD, DAS DIE ZAHL DER	 *
C		   ITERATIONSSCHRITTE FUER JEDEN EIGENWERT	 *
C		   AUFNIMMT. FALLS ZWEI EIGENWERTE ALS PAAR	 *
C		   GLEICHZEITIG GEFUNDEN WERDEN, DANN WIRD	 *
C		   DIE ZAHL DER ITERATIONSSCHRITTE MIT EINEM	 *
C		   POSITIVEN VORZEICHEN FUER DEN ERSTEN UND	 *
C		   EINEM NEGATIVEN VORZEICHEN FUER DEN		 *
C		   ZWEITEN EIGENWERT EINGETRAGEN.		 *
C     EIVEC:	   EIN (N,N)-FELD, DAS DIE NICHTNORMALISIER-	 *
C		   TEN EIGENVEKTOREN DER URSPRUENGLICHEN	 *
C		   VOLLEN MATRIX AUFNIMMT (FALLS H NICHT DIE	 *
C		   AUSGANGSMATRIX WAR). FALLS DER I-TE EIGEN-	 *
C		   WERT REELL IST, DANN IST DIE I-TE SPALTE	 *
C		   VON EIVEC DER DAZUGEHOERIGE REELLE EIGEN-	 *
C		   VEKTOR. FALLS DIE EIGENWERTE I UND I+1 EIN	 *
C		   KOMPLEXES PAAR BILDEN, GEBEN I-TE UND	 *
C		   (I+1)-TE SPALTE REAL- UND IMAGINAERTEIL	 *
C		   DESJENIGEN EIGENVEKTORS AN, DER ZU DEM	 *
C		   EIGENWERT MIT POSITIVEM IMAGINAERTEIL	 *
C		   GEHOERT.					 *
C								 *
C     RUECKGABEWERTE:						 *
C     ===============						 *
C     0:	   KEIN FEHLER					 *
C     1:	   DIE PARAMETER N, LOW ODER HIGH HABEN 	 *
C		   UNERLAUBTE WERTE.				 *
C     2:	   ALLE EIGENVEKTOREN SIND DER NULLVEKTOR.	 *
C     3:	   DIE MAXIMALE SCHRITTZAHL IST UEBERSCHRITTEN.  *
C								 *
C     LOKALE GROESSEN:						 *
C     ================						 *
C     ZERO,ONE,TWO,PT75,PT4375: WICHTIGE GLEITKOMMAKONSTANTEN	 *
C     MAXSTP:			KONSTANTE FUER DIE MAXIMALE	 *
C				SCHRITTZAHL			 *
C     I,J,K,L,M,N,NA,EN:	INDEXVARIABLEN			 *
C     ITER:			SCHRITTZAEHLER			 *
C     P,Q,R,S,T,W,X,Y,Z,NORM,					 *
C	RA,SA,VR,VI:		HILFSVARIABLEN FUER GLEITKOMMA-  *
C				BERECHNUNGEN			 *
C								 *
C----------------------------------------------------------------*
C								 *
C  BENOETIGTE UNTERPROGRAMME: COMDIV				 *
C								 *
C								 *
C  QUELLEN : PETERS, G. UND WILKINSON, J. H., SIEHE [PETE70].	 *
C								 *
C*****************************************************************
C								 *
C  AUTOR     : JUERGEN DIETEL					 *
C  DATUM     : 10.04.1987					 *
C  QUELLCODE : FORTRAN 77					 *
C								 *
C*****************************************************************
C
      INTEGER ND,N,LOW,HIGH,CNT(N)
      DOUBLE PRECISION H(ND,N),EIVEC(ND,N),WERTR(N),WERTI(N),EPS
      DOUBLE PRECISION ZERO,ONE,TWO,PT75,PT4375
      INTEGER MAXSTP
      PARAMETER (ZERO = 0.0D0,ONE = 1.0D0,TWO = 2.0D0,
     *		 PT75 = 0.75D0,PT4375 = 0.4375D0,MAXSTP = 100)
      INTEGER I,J,K,L,M,NA,EN,ITER
      DOUBLE PRECISION P,Q,R,S,T,W,X,Y,Z,NORM,RA,SA,VR,VI
C
C     FEHLER 1: DIE PARAMETER N, LOW ODER HIGH HABEN UNERLAUBTE
C     WERTE:
C
      IF (N .LT. 1 .OR. LOW .LT. 1 .OR. HIGH .GT. N) THEN
	 HQR2 = 1
	 RETURN
      ENDIF
C
C     VORBESETZUNG FUER DIE BEI DER AUSBALANCIERUNG GEFUNDENEN
C     ISOLIERTEN EIGENWERTE:
C
      DO 10 I=1,N
	 IF (I .LT. LOW .OR. I .GT. HIGH) THEN
	    WERTR(I) = H(I,I)
	    WERTI(I) = ZERO
	    CNT(I) = 0
	 ELSE
	    CNT(I) = -1
	 ENDIF
   10	 CONTINUE
C
      EN = HIGH
      T = ZERO
   15	 IF (EN .LT. LOW) GOTO 333
	 ITER = 0
	 NA = EN-1
C
C	    NACH EINEM EINZELNEN KLEINEN SUBDIAGONALELEMENT
C	    SUCHEN:
C
   20	    DO 30 L=EN,LOW+1,-1
   30	       IF (ABS(H(L,L-1)) .LE. EPS*
     *		   (ABS(H(L-1,L-1))+ABS(H(L,L)))) GOTO 40
   40	    X = H(EN,EN)
	    IF (L .EQ. EN) THEN
C
C	       EINE WURZEL GEFUNDEN:
C
	       WERTR(EN) = X + T
	       H(EN,EN) = WERTR(EN)
	       WERTI(EN) = ZERO
	       CNT(EN) = ITER
	       EN = NA
	       GOTO 15
	    ENDIF
C
	    Y = H(NA,NA)
	    W = H(EN,NA)*H(NA,EN)
	    IF (L .EQ. NA) THEN
C
C	       ZWEI WURZELN GEFUNDEN:
C
	       P = (Y-X)/TWO
	       Q = P*P+W
	       Z = SQRT(ABS(Q))
	       H(EN,EN) = X+T
	       X = H(EN,EN)
	       H(NA,NA) = Y+T
	       CNT(EN) = -ITER
	       CNT(NA) = ITER
	       IF (Q .GE. ZERO) THEN
C
C		  EIN REELLES PAAR GEFUNDEN:
C
		  IF (P .LT. ZERO) Z = -Z
		  Z = P+Z
		  WERTR(NA) = X+Z
		  WERTR(EN) = X-W/Z
		  WERTI(NA) = ZERO
		  WERTI(EN) = ZERO
		  X = H(EN,NA)
		  R = SQRT(X*X+Z*Z)
		  P = X/R
		  Q = Z/R
C
C		  ZEILENMODIFIKATION:
C
		  DO 50 J=NA,N
		     Z = H(NA,J)
		     H(NA,J) = Q*Z+P*H(EN,J)
		     H(EN,J) = Q*H(EN,J)-P*Z
   50		     CONTINUE
C
C		  SPALTENMODIFIKATION:
C
		  DO 60 I=1,EN
		     Z = H(I,NA)
		     H(I,NA) = Q*Z+P*H(I,EN)
		     H(I,EN) = Q*H(I,EN)-P*Z
   60		     CONTINUE
C
C		  AKKUMULATION:
C
		  DO 70 I=LOW,HIGH
		     Z = EIVEC(I,NA)
		     EIVEC(I,NA) = Q*Z+P*EIVEC(I,EN)
		     EIVEC(I,EN) = Q*EIVEC(I,EN)-P*Z
   70		     CONTINUE
	       ELSE
C
C		  KOMPLEXES PAAR:
C
		  WERTR(NA) = X+P
		  WERTR(EN) = WERTR(NA)
		  WERTI(NA) = Z
		  WERTI(EN) = -Z
	       ENDIF
	       EN = EN-2
	       GOTO 15
	    ENDIF
C
	    IF (ITER .EQ. MAXSTP) THEN
C
C	       FEHLER 3: MAXIMALE SCHRITTZAHL UEBERSCHRITTEN:
C
	       CNT(EN) = MAXSTP+1
	       HQR2 = 3
	       RETURN
	    ENDIF
	    IF (MOD(ITER,10) .EQ. 0 .AND. ITER .NE. 0) THEN
C
C	       EINEN UNGEWOEHNLICHEN SHIFT DURCHFUEHREN:
C
	       T = T+X
	       DO 80 I=LOW,EN
   80		  H(I,I) = H(I,I)-X
	       S = ABS(H(EN,NA))+ABS(H(NA,EN-2))
	       X = PT75*S
	       Y = X
	       W = -PT4375*S*S
	    ENDIF
	    ITER = ITER+1
C
C	    NACH ZWEI AUFEINANDERFOLGENDEN KLEINEN
C	    SUBDIAGONALELEMENTEN SUCHEN:
C
	    DO 90 M=EN-2,L,-1
	       Z = H(M,M)
	       R = X-Z
	       S = Y-Z
	       P = (R*S-W)/H(M+1,M)+H(M,M+1)
	       Q = H(M+1,M+1)-Z-R-S
	       R = H(M+2,M+1)
	       S = ABS(P)+ABS(Q)+ABS(R)
	       P = P/S
	       Q = Q/S
	       R = R/S
	       IF (M .EQ. L) GOTO 100
	       IF (ABS(H(M,M-1))*(ABS(Q)+ABS(R)) .LE. EPS*ABS(P)*
     *		   (ABS(H(M-1,M-1))+ABS(Z)+ABS(H(M+1,M+1))))
     *		   GOTO 100
   90	       CONTINUE
  100	    DO 110 I=M+2,EN
  110	       H(I,I-2) = ZERO
	    DO 120 I=M+3,EN
  120	       H(I,I-3) = ZERO
C
C	    EIN DOPPELTER QR-SCHRITT, DER DIE ZEILEN L BIS EN UND
C	    DIE SPALTEN M BIS EN DES GANZEN FELDES BETRIFFT:
C
	    DO 200 K=M,NA
	       IF (K .NE. M) THEN
		  P = H(K,K-1)
		  Q = H(K+1,K-1)
		  IF (K .NE. NA) THEN
		     R = H(K+2,K-1)
		  ELSE
		     R = ZERO
		  ENDIF
		  X = ABS(P)+ABS(Q)+ABS(R)
		  IF (X .EQ. ZERO) GOTO 200
		  P = P/X
		  Q = Q/X
		  R = R/X
	       ENDIF
	       S = SQRT(P*P+Q*Q+R*R)
	       IF (P .LT. ZERO) S = -S
	       IF (K .NE. M) THEN
		  H(K,K-1) = -S*X
	       ELSEIF (L .NE. M) THEN
		  H(K,K-1) = -H(K,K-1)
	       ENDIF
	       P = P+S
	       X = P/S
	       Y = Q/S
	       Z = R/S
	       Q = Q/P
	       R = R/P
C
C	       ZEILENMODIFIKATION:
C
	       DO 130 J=K,N
		  P = H(K,J)+Q*H(K+1,J)
		  IF (K .NE. NA) THEN
		     P = P+R*H(K+2,J)
		     H(K+2,J) = H(K+2,J)-P*Z
		  ENDIF
		  H(K+1,J) = H(K+1,J)-P*Y
		  H(K,J) = H(K,J)-P*X
  130		  CONTINUE
	       J = MIN(K+3,EN)
C
C	       SPALTENMODIFIKATION:
C
	       DO 140 I=1,J
		  P = X*H(I,K)+Y*H(I,K+1)
		  IF (K .NE. NA) THEN
		     P = P+Z*H(I,K+2)
		     H(I,K+2) = H(I,K+2)-P*R
		  ENDIF
		  H(I,K+1) = H(I,K+1)-P*Q
		  H(I,K) = H(I,K)-P
  140		  CONTINUE
C
C	       TRANSFORMATIONEN AKKUMULIEREN:
C
	       DO 150 I=LOW,HIGH
		  P = X*EIVEC(I,K)+Y*EIVEC(I,K+1)
		  IF (K .NE. NA) THEN
		     P = P+Z*EIVEC(I,K+2)
		     EIVEC(I,K+2) = EIVEC(I,K+2)-P*R
		  ENDIF
		  EIVEC(I,K+1) = EIVEC(I,K+1)-P*Q
		  EIVEC(I,K) = EIVEC(I,K)-P
  150		  CONTINUE
  200	       CONTINUE
	    GOTO 20
C
C
C     ALLE WURZELN GEFUNDEN, NUN WIRD RUECKTRANSFORMIERT:
C
C     1-NORM VON H BESTIMMEN:
C
  333 NORM = ZERO
      K=1
      DO 201 I=1,N
	 DO 101 J=K,N
  101	    NORM = NORM+ABS(H(I,J))
  201	 K = I
      IF (NORM .EQ. ZERO) THEN
C	 FEHLER 2: 1-NORM VON H IST GLEICH 0:
	 HQR2 = 2
	 RETURN
      ENDIF
C
C     RUECKTRANSFORMATION:
C
      DO 207 EN=N,1,-1
	 P = WERTR(EN)
	 Q = WERTI(EN)
	 NA = EN - 1
	 IF (Q .EQ. ZERO) THEN
C
C	    REELLER VEKTOR:
C
	    M = EN
	    H(EN,EN) = ONE
	    DO 63 I=NA,1,-1
	       W = H(I,I)-P
	       R = H(I,EN)
	       DO 38 J=M,NA
   38		  R = R+H(I,J)*H(J,EN)
	       IF (WERTI(I) .LT. ZERO) THEN
		  Z = W
		  S = R
	       ELSE
		  M = I
		  IF (WERTI(I) .EQ. ZERO) THEN
		     IF (W .NE. ZERO) THEN
			H(I,EN) = -R/W
		     ELSE
			H(I,EN) = -R/(EPS*NORM)
		     ENDIF
		  ELSE
C
C		     LOESE DAS GLEICHUNGSSYSTEM:
C		     [ W   X ] [ H(I,EN)   ]   [ -R ]
C		     [	     ] [	   ] = [    ]
C		     [ Y   Z ] [ H(I+1,EN) ]   [ -S ]
C
		     X = H(I,I+1)
		     Y = H(I+1,I)
		     Q = (WERTR(I)-P)*(WERTR(I) - P)+
     *			 WERTI(I)*WERTI(I)
		     T = (X*S-Z*R)/Q
		     H(I,EN) = T
		     IF (ABS(X) .GT. ABS(Z)) THEN
			H(I+1,EN) = (-R-W*T)/X
		     ELSE
			H(I+1,EN) = (-S-Y*T)/Z
		     ENDIF
		  ENDIF
	       ENDIF
   63	       CONTINUE
	 ELSEIF (Q .LT. ZERO) THEN
C
C	    KOMPLEXER VEKTOR, DER ZU LAMBDA = P - I * Q GEHOERT:
C
	    M = NA
	    IF (ABS(H(EN,NA)) .GT. ABS(H(NA,EN))) THEN
	       H(NA,NA) = -(H(EN,EN)-P)/H(EN,NA)
	       H(NA,EN) = -Q/H(EN,NA)
	    ELSE
	       CALL COMDIV(-H(NA,EN),ZERO,H(NA,NA)-P,Q,
     *			 H(NA,NA),H(NA,EN))
	    ENDIF
	    H(EN,NA) = ONE
	    H(EN,EN) = ZERO
	    DO 190 I=NA-1,1,-1
	       W = H(I,I)-P
	       RA = H(I,EN)
	       SA = ZERO
	       DO 75 J=M,NA
		  RA = RA+H(I,J)*H(J,NA)
		  SA = SA+H(I,J)*H(J,EN)
   75		  CONTINUE
	       IF (WERTI(I) .LT. ZERO) THEN
		  Z = W
		  R = RA
		  S = SA
	       ELSE
		  M = I
		  IF (WERTI(I) .EQ. ZERO) THEN
		      CALL COMDIV(-RA,-SA,W,Q,H(I,NA),H(I,EN))
		  ELSE
C
C		     LOESE DIE KOMPLEXEN GLEICHUNGEN:
C	    [ W+Q*I   X   ] [H(I,NA)+H(I,EN)*I	  ]   [-RA-SA*I]
C	    [		  ] [			  ] = [        ]
C	    [	Y   Z+Q*I ] [H(I+1,NA)+H(I+1,EN)*I]   [-R-S*I  ]
C
		     X = H(I,I+1)
		     Y = H(I+1,I)
		     VR = (WERTR(I)-P)*(WERTR(I)-P)+
     *			  WERTI(I)*WERTI(I)-Q*Q
		     VI = TWO*Q*(WERTR(I)-P)
		     IF (VR .EQ. ZERO .AND. VI .EQ. ZERO) VR =
     *			EPS*NORM*
     *			(ABS(W)+ABS(Q)+ABS(X)+ABS(Y)+ABS(Z))
		     CALL COMDIV(X*R-Z*RA+Q*SA,X*S-Z*SA-Q*RA,
     *				 VR,VI,H(I,NA),H(I,EN))
		     IF (ABS(X) .GT. ABS(Z)+ABS(Q)) THEN
			H(I+1,NA) = (-RA-W*H(I,NA)+Q*H(I,EN))/X
			H(I+1,EN) = (-SA-W*H(I,EN)-Q*H(I,NA))/X
		     ELSE
		       CALL COMDIV(-R-Y*H(I,NA),-S-Y*H(I,EN),Z,Q,
     *				 H(I+1,NA),H(I+1,EN))
		     ENDIF
		  ENDIF
	       ENDIF
  190	       CONTINUE
	 ENDIF
  207	 CONTINUE
C
C     ZU DEN ISOLIERTEN WURZELN GEHOERIGE VEKTOREN:
C
      DO 230 I=1,N
	 IF (I .LT. LOW .OR. I .GT. HIGH) THEN
	    DO 220 J=I+1,N
  220	       EIVEC(I,J) = H(I,J)
	 ENDIF
  230	 CONTINUE
C
C     MIT DER TRANSFORMATIONSMATRIX MULTIPLIZIEREN, UM DIE
C     VEKTOREN DER URSPRUENGLICHEN VOLLEN MATRIX ZU ERHALTEN:
C
      DO 300 J=N,LOW,-1
	 IF (J .LE. HIGH) THEN
	    M = J
	 ELSE
	    M = HIGH
	 ENDIF
	 L = J-1
	 IF (WERTI(J) .LT. ZERO) THEN
	    DO 330 I=LOW,HIGH
	       Y = ZERO
	       Z = ZERO
	       DO 320 K=LOW,M
		  Y = Y+EIVEC(I,K)*H(K,L)
		  Z = Z+EIVEC(I,K)*H(K,J)
  320		  CONTINUE
	       EIVEC(I,L) = Y
	       EIVEC(I,J) = Z
  330	       CONTINUE
	 ELSE
	    IF (WERTI(J) .EQ. ZERO) THEN
	       DO 350 I=LOW,HIGH
		  Z = ZERO
		  DO 340 K=LOW,M
  340		     Z =Z+EIVEC(I,K)*H(K,J)
  350		  EIVEC(I,J) = Z
	    ENDIF
	 ENDIF
  300	 CONTINUE
C
C     RUECKGABE VON 0: KEIN FEHLER:
C
      HQR2 = 0
      END
C
C
      SUBROUTINE COMDIV (AR,AI,BR,BI,RESR,RESI)
C
C*****************************************************************
C								 *
C     ZWECK DES PROGRAMMS:					 *
C     ====================					 *
C     KOMPLEXE DIVISION: RESR+I*RESI := (AR+I*AI)/(BR+I*BI).	 *
C     (DIESE PROZEDUR SOLLTE NICHT MIT				 *
C      BR=BI=0 AUFGERUFEN WERDEN.)				 *
C								 *
C     EINGABEPARAMETER: 					 *
C     ================= 					 *
C     AR,AI:	 REAL- UND IMAGINAERTEIL DES DIVIDENDEN 	 *
C     BR,BI:	 REAL- UND IMAGINAERTEIL DES DIVISORS		 *
C								 *
C     AUSGABEPARAMETER: 					 *
C     ================= 					 *
C     RESR,RESI: REAL- UND IMAGINAERTEIL DES QUOTIENTEN 	 *
C								 *
C     LOKALE GROESSEN:						 *
C     ================						 *
C     ZERO:		 GLEITKOMMAKONSTANTE 0			 *
C     TEMP1,TEMP2,TEMP3: HILFSVARIABLEN ZUR SPEICHERUNG VON	 *
C			 ZWISCHENERGEBNISSEN			 *
C								 *
C----------------------------------------------------------------*
C								 *
C  BENOETIGTE UNTERPROGRAMME: KEINE				 *
C								 *
C								 *
C  QUELLEN : MARTIN, R. S. UND WILKINSON, J. H., SIEHE [MART70]. *
C								 *
C*****************************************************************
C								 *
C  AUTOR     : JUERGEN DIETEL					 *
C  DATUM     : 10.04.1987					 *
C  QUELLCODE : FORTRAN 77					 *
C								 *
C*****************************************************************
C
      DOUBLE PRECISION AR,AI,BR,BI,RESR,RESI
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.0D0)
      DOUBLE PRECISION TEMP1,TEMP2,TEMP3
C
      IF (BR .EQ. ZERO .AND. BI .EQ. ZERO) THEN
	 RESR = ZERO
	 RESI = ZERO
	 RETURN
      ENDIF
      IF (ABS(BR) .GT. ABS(BI)) THEN
	 TEMP1 = BI/BR
	 TEMP2 = TEMP1*BI+BR
	 TEMP3 = (AR+TEMP1*AI)/TEMP2
	 RESI = (AI-TEMP1*AR)/TEMP2
	 RESR = TEMP3
      ELSE
	 TEMP1 = BR/BI
	 TEMP2 = TEMP1*BR+BI
	 TEMP3 = (TEMP1*AR+AI)/TEMP2
	 RESI = (TEMP1*AI-AR)/TEMP2
	 RESR = TEMP3
      ENDIF
      END
C
C
      DOUBLE PRECISION FUNCTION COMABS(AR,AI)
C
C*****************************************************************
C								 *
C     ZWECK DES PROGRAMMS:					 *
C     ====================					 *
C     BERECHNUNG DES BETRAGES DER KOMPLEXEN ZAHL AR+I*AI:	 *
C     COMABS:=SQRT(AR*AR+AI*AI) 				 *
C								 *
C     EINGABEPARAMETER: 					 *
C     ================= 					 *
C     AR,AI: REAL- UND IMAGINAERTEIL DER KOMPLEXEN ZAHL, DEREN	 *
C	     BETRAG ZU BERECHNEN IST				 *
C								 *
C     AUSGABEPARAMETER: 					 *
C     ================= 					 *
C     KEINE							 *
C								 *
C     RUECKGABEWERT:						 *
C     ==============						 *
C     BETRAG DES KOMPLEXEN PARAMETERS				 *
C								 *
C     LOKALE GROESSEN:						 *
C     ================						 *
C     ZERO,ONE:    KONSTANTEN					 *
C     TEMP1,TEMP2: HILFSVARIABLEN ZUR SPEICHERUNG VON		 *
C		   ZWISCHENERGEBNISSEN				 *
C								 *
C----------------------------------------------------------------*
C								 *
C  BENOETIGTE UNTERPROGRAMME: SWAP				 *
C								 *
C								 *
C  QUELLEN : MARTIN, R. S. UND WILKINSON, J. H., SIEHE [MART70]. *
C								 *
C*****************************************************************
C								 *
C  AUTOR     : JUERGEN DIETEL					 *
C  DATUM     : 10.04.1987					 *
C  QUELLCODE : FORTRAN 77					 *
C								 *
C*****************************************************************
C
      DOUBLE PRECISION AR,AI
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO = 0.0D0,ONE = 1.0D0)
      DOUBLE PRECISION TEMP1,TEMP2
      TEMP1 = ABS(AR)
      TEMP2 = ABS(AI)
      IF (AR .EQ. ZERO .OR. AI .EQ. ZERO) THEN
	 COMABS = ZERO
	 RETURN
      ENDIF
      IF (TEMP2 .GT. TEMP1) CALL SWAP (TEMP1,TEMP2)
      IF (TEMP2 .EQ. ZERO) THEN
	 COMABS = TEMP1
      ELSE
	 COMABS = TEMP1*SQRT(ONE+(TEMP2/TEMP1)**2)
      ENDIF
      END
C
C
      INTEGER FUNCTION NORMAL (ND,N,V,WI)
C
C*****************************************************************
C								 *
C     ZWECK DES PROGRAMMS:					 *
C     ====================					 *
C     NORMAL NORMALISIERT DIE EIGENVEKTOREN IN DER MAXIMUMNORM.  *
C								 *
C     EINGABEPARAMETER: 					 *
C     ================= 					 *
C     ND:     FUEHRENDE DIMENSION DER MATRIX V, WIE SIE IM	 *
C	      HAUPTPROGRAMM VEREINBART WURDE			 *
C     N:      DIE ORDNUNG DER MATRIX V				 *
C     V:      EIN (N,N)-FELD VOM TYP DOUBLE PRECISION, DAS	 *
C	      SPALTENWEISE DIE EIGENVEKTOREN ENTHAELT		 *
C	      (SIEHE EIGEN)					 *
C     WI:     EIN FELD MIT N KOMPONENTEN VOM TYP		 *
C	      DOUBLE PRECISION, DESSEN KOMPONENTEN DIE		 *
C	      IMAGINAERTEILE DER EIGENWERTE SIND		 *
C								 *
C     AUSGABEPARAMETER: 					 *
C     ================= 					 *
C     V:      MATRIX DER NORMALISIERTEN EIGENVEKTOREN		 *
C								 *
C     LOKALE GROESSEN:						 *
C     ================						 *
C     ZERO,ONE: GLEITKOMMAKONSTANTEN 0 UND 1			 *
C     I,J:	INDEXVARIABLEN					 *
C     MAXI:	HILFSVARIABLE ZUR BERECHNUNG DER REELLEN	 *
C		VEKTORNORM					 *
C     TR,TI:	HILFSVARIABLEN ZUR BERECHNUNG DER KOMPLEXEN	 *
C		VEKTORNORM					 *
C								 *
C----------------------------------------------------------------*
C								 *
C  BENOETIGTE UNTERPROGRAMME: COMABS, COMDIV			 *
C								 *
C*****************************************************************
C								 *
C  AUTOR     : JUERGEN DIETEL					 *
C  DATUM     : 10.04.1987					 *
C  QUELLCODE : FORTRAN 77					 *
C								 *
C*****************************************************************
C
      INTEGER ND,N
      DOUBLE PRECISION V(ND,N),WI(N)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO = 0.0D0,ONE = 1.0D0)
      INTEGER I,J
      DOUBLE PRECISION MAXI,TR,TI,COMABS
C
      J = 1
   10 IF (J .GT. N) GOTO 80
	 IF (WI(J) .EQ. ZERO) THEN
	     MAXI = V(1,J)
	     DO 15 I=2,N
   15		IF (ABS(V(I,J)) .GT. ABS(MAXI)) MAXI = V(I,J)
	     IF (MAXI .NE. ZERO) THEN
		MAXI = ONE/MAXI
		DO 20 I=1,N
   20		   V(I,J) = V(I,J)*MAXI
	     ENDIF
	     J = J+1
	 ELSE
	    TR = V(1,J)
	    TI = V(1,J+1)
	    DO 30 I=2,N
	       IF (COMABS(V(I,J),V(I,J+1)) .GT. COMABS(TR,TI))
     *	       THEN
		  TR = V(I,J)
		  TI = V(I,J+1)
	       ENDIF
   30	       CONTINUE
	    IF (TR .NE. ZERO .OR. TI .NE. ZERO) THEN
	       DO 40 I=1,N
   40		  CALL COMDIV (V(I,J),V(I,J+1),TR,TI,
     *			       V(I,J),V(I,J+1))
	    ENDIF
	    J = J+2
	 ENDIF
	 GOTO 10
   80 NORMAL = 0
      END
C
C
      SUBROUTINE SWAP(X,Y)
C
C*****************************************************************
C								 *
C     ZWECK DES PROGRAMMS:					 *
C     ====================					 *
C     SWAP VERTAUSCHT DIE WERTE DER BEIDEN DOUBLE PRECISION-	 *
C     VARIABLEN X UND Y.					 *
C								 *
C*****************************************************************
C
      DOUBLE PRECISION X,Y,TEMP
      TEMP = X
      X = Y
      Y = TEMP
      END
c
      subroutine assemb(s,ix,id,ie,ma,neq,nen1,ndf,nie,nen,n,ns,a)
      implicit double precision (a-h,o-z)
c
c.... assemble unsymmetric full matrix 'uneigv'
c
      dimension a(neq,*),s(ns,ns),ix(nen1,*),id(ndf,*),ie(nie,*)
c.... loop through the rows to perform the assembly
      do 200 i = 1,nen
	ii = ix(i,n)
	do 201 j = 1,ndf
	  jj = ie(j,ma)
	  if(jj .le. 0) go to 201
	  jr = (i-1)*ndf + j
	  irg = id(jj,ii)
	  if(irg.le.0) go to 201
	  do 202 k = 1,nen
	    kk = ix(k,n)
	    do 203 l = 1,ndf
	      ll = ie(l,ma)
	      if(ll .le. 0) go to 203
	      lc = (k-1)*ndf + l
	      icg = id(ll,kk)
	      if (icg.le.0) go to 203
!$OMP ATOMIC
	      a(irg,icg) = a(irg,icg) + s(jr,lc)
203	    continue
202	  continue
201	continue
200   continue
      return
      end
