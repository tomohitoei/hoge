C
c      implicit none
      implicit real*8(a-h,o-z)
      
      OPEN (UNIT=07,STATUS='OLD'    ,IOSTAT=IOSS,FILE='d-sig')
      OPEN (UNIT=10,STATUS='UNKNOWN',IOSTAT=IOSS,FILE='ans-sig')
C
C *** 1991  7/23 20:05  M GA N YORI MO SAKI NI NARU YOO NI KAETA
C *-* 1990  7/16  9:05  RG KEISUU WO PRINT SURU YOO NI SHITA
C *+* 1989  5/24 16:30  NATTOKU NO IKU YOO NI KANSEI SHITA  T.HAYASHIDA
C *** KYOYOO OORYOKU-DO HOO NI SHITAGATTE OORYOKU-DO WO KEISAN SURU ***
C aaaa
C ======================================================================
CC    WRITE(7,710)
CC710 FORMAT(1H1,/,1H0,25X,'M',7X,'N',7X,'S     B     H     D    DD     '
CC   1,'AS     ASD     U      UD  SIG-S  SIGSD  SIG-C  SIG-C2 '
CC   2,'SC2/SC  X   KTRY',/,1H )
C ======================================================================
      WRITE(10,610)
  610 FORMAT(1H1,/,1H0,5X,'M',7X,'N',7X,'S     B     H     D    DD     '
     1,'AS     ASD     U   SIG-S  SIGSD  SIG-C   TAU    '
     2,'TAU0    X      C      S     Z',/,1H )
C
      WRITE(6,*) ' ENTER N '
      READ(5,*) SN
C
C
    1 DUMMY=0.0
      READ(7,503,END=90) INUM,INNUM,IPOS,ISTEP,
     1                   BOM,BN,BOQ,B ,H ,D ,DD ,AS,ASD,U,UD
  503 FORMAT(4I5,11F12.0)
      BM=DABS(BOM)
      BQ=DABS(BOQ)
C ======================================================================
CC    WRITE(7,750) BM,BN,BQ,B,H,D,DD,AS,ASD,U,UD
C ======================================================================
C
CC    SUBROUTINE CALSIG(BN,BM,BQ,B0,H0,D0,DD0,AS,ASD,U,UD
CC   1                 ,SIGS,SIGSD,SIGC,TAU,TAU02,TAU01,X)
      AN=BN*1.000
      AM=BM*100.*1.000
      AQ=BQ*1.000
CC    B=B0*100.
CC    H=H0*100.
CC    D=D0*100.
CC    DD=DD0*100.
C>>>  SN=15.0
      CALL ALLONA(0.,SIGS,SIGSD,SIGC,TAU,TAU02,TAU01,X,DUM,DUM,DUM)
      IF(AS.LT..001) AS=0.001
C
      AAVAL=B*H
      AIVAL=B*H*H*H/12.
      ASIGNC=AN/AAVAL
      ASIGMC=(AM/AIVAL)*H/2.
      ASIGNC=DABS(ASIGNC)
C
      HW6=H/6.
      HW2=H/2.
      C=D-HW2
      CD=HW2-DD
      HW2PC=HW2+C
      HW2MC=HW2-C
      HW2PCD=HW2+CD
      HW2MCD=HW2-CD
      XLOW=0.
      XUP=D
C                        AMD1 : RG KEISUU C,S,Z WO MOTOMERU TAME NO M
      AMD1=AM+AN*C
C
      IF(ASIGNC.LT..0001) GO TO 5
      IF(AN.GT.0.) GO TO 30
C                                       N < 0  NO BAAI
      E=-AM/AN
      EPC=E+C
      EMC=E-C
      EPCD=E+CD
      EMCD=E-CD
      IF(E.LE.C) GO TO 22
C
      AK1=1.0
      AK2=-3.*(E+HW2)
      AK3=-6.*SN*(ASD*EPCD+AS*EMC)/B
      AK4=6.*SN*(ASD*EPCD*DD+AS*EMC*D)/B
C ======================================================================
CC    WRITE(7,704) AK1,AK2,AK3,AK4
CC704 FORMAT(1H ,10X,'AK1,AK2,AK3,AK4 ... ',4(1PE10.3))
C ======================================================================
      CALL CARDA1(AK1,AK2,AK3,AK4,X,XLOW,XUP,9,KKOS)
      GO TO 8
C                                       N = 0  NO BAAI  NO X NO KEISAN
    5 DUMMY=0.0
      AN=0.0
      KTRY=99
      IF(ASD.GT..001) GO TO 6
      A1=1.+2.*B*D/(SN*AS)
      A2=DSQRT(A1)
      X=-SN*AS*(1.-A2)/B
      GO TO 7
    6 DUMMY=0.0
      A1=SN*(AS+ASD)/B
      A2=A1*A1+2.*SN*(D*AS+DD*ASD)/B
      X=-A1+DSQRT(A2)
    7 DUMMY=0.0
      A2=.5*B*X*(HW2-X/3.)+SN*ASD*(HW2-DD)*(X-DD)/X
     1                    +SN*AS*(D-HW2)*(D-X)/X
      SIGC=AM/A2
      SIGS=SN*SIGC*(C+HW2-X)/X
c      SIGSD=SN*SIGC*(X-HW2+C)/X
      SIGSD=SN*SIGC*(X-DD)/X
      GO TO 12
C                                       N < 0  NO BAAI
    8 DUMMY=0.0
      KTRY=0
   10 DUMMY=0.0
      A1=.5*B*X*X+SN*ASD*(X-DD)-SN*AS*(D-X)
      SIGC=X*AN/A1
      SIGS=SN*SIGC*(C+HW2-X)/X
c      SIGSD=SN*SIGC*(X-HW2+C)/X
      SIGSD=SN*SIGC*(X-DD)/X
      A2=.5*B*X*(HW2-X/3.)+SN*ASD*(HW2-DD)*(X-DD)/X
     1                    +SN*AS*(D-HW2)*(D-X)/X
      SIGC2=AM/A2
      A1=0.0
      IF(SIGC.LT.-1.E-6.OR.SIGC.GT.1.E-6) A1=SIGC2/SIGC
C ======================================================================
CC    WRITE(7,750) BM,BN,BQ,B,H,D,DD,AS,ASD,U,UD,SIGS,SIGSD,SIGC,SIGC2
CC   1            ,A1,X,KTRY
CC750 FORMAT(1H ,3F8.1,3F6.0,F5.0,2F7.2,2F7.1,2F7.1,2F7.2,F7.4,F7.2,I4)
C ======================================================================
      IF(KTRY.GE.1.AND.SIGC.LT.0.) GO TO 11
      IF(KTRY.GE.1) GO TO 12
      IF(A1.LT.0.99 .OR.A1.GT.1.01 ) GO TO 14
      GO TO 12
C
   11 DUMMY=0.0
      SIGC=SIGC2
      SIGS=SN*SIGC*(C+HW2-X)/X
c      SIGSD=SN*SIGC*(X-HW2+C)/X
      SIGSD=SN*SIGC*(X-DD)/X
C                                N < 0  &  N = 0  NO BAAI NO  TAU & TAU0
   12 DUMMY=0.0
c      AJ=7./8.
      AJ=1./1.15
      X1=X/D
c      IF(X1.GT.1.E-6.AND.X1.LT.(1.-1.E-6)) AJ=1.-X1/3.
      TAU=AQ/(B*AJ*D)
      A1=AQ/AAVAL
c      IF(DABS(A1).LT..003) TAU=0.0
      TAU02=AQ/(U*AJ*D)
CC    TAU01=0.0
      GO TO 50
C                                       TRAIAL KEISAN
   14 DUMMY=0.0
      KTRY=0
      Y=0.45*H
   16 DUMMY=0.0
      FY=AK1*Y*Y*Y+AK2*Y*Y+AK3*Y+AK4
      DY=3.*AK1*Y*Y+2.*AK2*Y+AK3
C ======================================================================
CC    YY2=0.0
CC    IF(DABS(DY).GE.1.E-8) YY2=Y-FY/DY
CC    YY3=YY2-Y
CC    WRITE(7,716) KTRY,Y,FY,DY,YY2,YY3
CC716 FORMAT(1H ,12X,'KTRY,Y,FY,DY,Y2,(Y2-Y)=-FY/DY ... ',I3,5(1PE10.3))
C ======================================================================
      IF(DABS(DY).LT.1.E-8) GO TO 18
      Y2=Y-FY/DY
      IF(DABS((Y2-Y)).LE.1.E-7) GO TO 18
      Y=Y2
      KTRY=KTRY+1
      IF(KTRY.GT.50) GO TO 18
      GO TO 16
   18 DUMMY=0.0
      X=Y
      GO TO 10
C                                       ZEN-DANMEN HIPPARI
   22 DUMMY=0.0
      IF(ASD.LT..001) ASD=0.001
      AL0=CD+C
      AL2=-EMC
      AL1=EPCD
      AN1=-AN*AL2/AL0
      AN2=-AN*AL1/AL0
      SIGS=AN2/AS
      SIGSD=-AN1/ASD
CC    SIGC=0.0
CC    SIGC2=0.0
CC    X=0.0
      GO TO 12
C                                       N > 0  NO BAAI
   30 DUMMY=0.0
      E=AM/AN
      IF(DABS(E).LT..005.AND.DABS(AM).LT..01) GO TO 44
      IF(DABS(E).LE.HW6) GO TO 46
      EPC=E+C
      EMC=E-C
      EPCD=E+CD
      EMCD=E-CD
C
      AK1=1.0
      AK2=3.*(E-HW2)
      AK3=6.*SN*(ASD*EMCD+AS*EPC)/B
      AK4=-6.*SN*(ASD*EMCD*DD+AS*EPC*D)/B
C ======================================================================
CC    WRITE(7,704) AK1,AK2,AK3,AK4
C ======================================================================
      CALL CARDA1(AK1,AK2,AK3,AK4,X,XLOW,XUP,9,KKOS)
C
      KTRY=0
   32 DUMMY=0.0
      A1=.5*B*X*X+SN*ASD*(X-DD)-SN*AS*(D-X)
      SIGC=X*AN/A1
      SIGS=SN*SIGC*(C+HW2-X)/X
c      SIGSD=SN*SIGC*(X-HW2+C)/X
      SIGSD=SN*SIGC*(X-DD)/X
      A2=.5*B*X*(HW2-X/3.)+SN*ASD*(HW2-DD)*(X-DD)/X
     1                    +SN*AS*(D-HW2)*(D-X)/X
      SIGC2=AM/A2
      A1=0.0
      IF(SIGC.LT.-1.E-6.OR.SIGC.GT.1.E-6) A1=SIGC2/SIGC
C ======================================================================
CC    WRITE(7,750) BM,BN,BQ,B,H,D,DD,AS,ASD,U,UD,SIGS,SIGSD,SIGC,SIGC2
CC   1            ,A1,X,KTRY
C ======================================================================
      IF(KTRY.GE.1.AND.SIGC.LT.0.) GO TO 33
      IF(KTRY.GE.1) GO TO 34
      IF(A1.LT.0.99 .OR.A1.GT.1.01 ) GO TO 36
      GO TO 34
C
   33 DUMMY=0.0
      SIGC=SIGC2
      SIGS=SN*SIGC*(C+HW2-X)/X
c      SIGSD=SN*SIGC*(X-HW2+C)/X
      SIGSD=SN*SIGC*(X-DD)/X
C                                N > 0  &  HIPPARI-AS > 0  NO TAU & TAU0
   34 DUMMY=0.0
CC    GI=B*X*X/2.+(SN-1.)*ASD*(X-DD)-(SN-1.)*AS*(D-X)      '89 4/3 DUMMY
CC    AI=B*X+(SN-1.)*(AS+ASD)                              '89 4/3 DUMMY
CC    V0=GI/AI                                             '89 4/3 DUMMY
CC    TAU=((AQ/AI)*(B/2.)*(X-V0)*(X-V0)+SN*ASD*(X-DD-V0))  '89 4/3 DUMMY
CC   1                               /(B*V0*(E-H/2.+X-V0)) '89 4/3 DUMMY
c      AJ=7./8.
      AJ=1./1.15
      X1=X/D
c      IF(X1.GT.1.E-6.AND.X1.LT.(1.-1.E-6)) AJ=1.-X1/3.
CC    TAU2=AQ/(B*AJ*D)                                     '89 4/3 DUMMY
      TAU=AQ/(B*AJ*D)
CC    IF(TAU2.GT.TAU) TAU=TAU2                             '89 4/3 DUMMY
      A1=AQ/AAVAL
c      IF(DABS(A1).LT..003) TAU=0.0
CC    TAU01=(SN*ASD*AQ/(UD*GI))*(GI-(X-DD)*AI)/(GI-(X+E-H/2.)*AI) '89 4/
CC    TAU02=-(SN*AS*AQ/(U*GI))*(GI+(D-X)*AI)/(GI-(X+E-H/2.)*AI)   '89 4/
      TAU02=AQ/(U*AJ*D)
      TAU01=0.0
      GO TO 50
C                                       TRIAL KEISAN
   36 DUMMY=0.0
      KTRY=0
      Y=0.6*H
   38 DUMMY=0.0
      FY=AK1*Y*Y*Y+AK2*Y*Y+AK3*Y+AK4
      DY=3.*AK1*Y*Y+2.*AK2*Y+AK3
C ======================================================================
CC    YY2=0.0
CC    IF(DABS(DY).GE.1.E-8) YY2=Y-FY/DY
CC    YY3=YY2-Y
CC    WRITE(7,716) KTRY,Y,FY,DY,YY2,YY3
C ======================================================================
      IF(DABS(DY).LT.1.E-8) GO TO 40
      Y2=Y-FY/DY
      IF(DABS((Y2-Y)).LE.1.E-7) GO TO 40
      Y=Y2
      KTRY=KTRY+1
      IF(KTRY.GT.50) GO TO 40
      GO TO 38
   40 DUMMY=0.0
      X=Y
      GO TO 32
C                                       M = 0  NO BAAI
   44 DUMMY=0.0
      X=0.0
      SIGC=AN/B/H
      GO TO 34
C                             N > 0  & ZEN DANMEN ASSHUKU NO BAAI
   46 DUMMY=0.0
      X=0.0
      AI=B*H+(SN-1.)*(AS+ASD)
      Y0=(.5*B*H*H+(SN-1.)*(AS*D+ASD*DD))/AI
      HMY0=H-Y0
      DMY0=D-Y0
      Y0MDD=Y0-DD
      AII=(Y0*Y0*Y0+HMY0*HMY0*HMY0)*B/3.
     1    +SN*(AS*DMY0*DMY0+ASD*Y0MDD*Y0MDD)
      SIGC=AN/AI+AM*Y0/AII
      SIGCD=AN/AI-AM*HMY0/AII
      SIGS=-SN*(AN/AI-AM*DMY0/AII)
      SIGSD=SN*(AN/AI+AM*Y0MDD/AII)
      GO TO 34
C                                       KEKKA NO PRINT
   50 DUMMY=0.0
C                                       RG KEISUU NO KEISAN
      ACVAL=0.0
      ASVAL=0.0
      AZVAL=0.0
      A1=AMD1/(B*D*D)
      IF(DABS(A1).LT.1.E-3) GO TO 52
      ACVAL=SIGC/A1
      ASVAL=SIGS/(SN*A1)
      AZVAL=3./(3.-X/D)
   52 DUMMY=0.0
C
      SIGSD=-SIGSD
      if (0.lt.BOM) go to 12345
      AAATEMP=SIGS
      SIGS=SIGSD
      SIGSD=AAATEMP
12345 DUMMY=0.0
c
      WRITE(10,650) INUM,INNUM,IPOS,ISTEP,
     1              BM,BN,BQ,B,H,D,DD,AS,ASD,U,
     2              SIGS*10, SIGSD*10, SIGC*10, TAU*10, TAU02*10,
     3              X, ACVAL, ASVAL, AZVAL,BOM,BOQ
c  650 FORMAT(1H ,4I5, 3F8.1,3F6.0,F5.0,2F7.2, F7.1,2F7.1,4F7.2,2F7.2,F7.3)
  650 FORMAT(4I5, 21E15.8)
C
CC    WRITE(10,650) BM,BN,BQ,B,H,D,DD,AS,ASD,U,UD,SIGS,SIGSD,SIGC  * 1990
CC   1            ,TAU,TAU02,TAU01,X,KTRY                         * 8/16
CC650 FORMAT(1H ,3E8.2,3F6.0,F5.0,2F7.2,2F7.1,2F7.1,5F7.2,I4)     *DUMMY
CC    RETURN
      GO TO 1
C
   90 DUMMY=0.0
      WRITE(6,190)
  190 FORMAT('  *** END OF JOB ***')
      STOP
      END
C --------------------------- 3-JI HOOTEI SHIKI  ( CARDANO METHOD )
      SUBROUTINE CARDA1(A3,A2,A1,A0,XANS,XANMIN,XANMAX,KANCOD,KKOS)
      implicit real*8(a-h,o-z)
      real*8 A3,A2,A1,A0,XANS,XANMIN,XANMAX
      real*8 U3,V3,AL,BE,PAI,TH
      DIMENSION XX(3),XX2(3)
      PAI=3.14159265
      DO 28 I=1,3
      XX(I)=0.0
      XX2(I)=0.0
   28 CONTINUE
C                                  SHIKI NO HENKEI
      XANS=0.0
      A=A2/A3
      B=A1/A3
      C=A0/A3
      A4=A/3.
C                                  HANBETSU SHIKI
      P=B/3.-A*A/9.
      Q=C-A*B/3.+2.*A*A*A/27.
      D=Q*Q+4.*P*P*P
      DSQR=0.0
      IF(D) 1,3,2
C                                  3-JIKKON
    1 TH=ATAN(-DSQRT(-D)/Q)
      Z=2.*DSQRT(-P)
      X1=Z*DCOS(TH/3.)-A4
      X2=-Z*DCOS((PAI-TH)/3.)-A4
      X3=-Z*DCOS((PAI+TH)/3.)-A4
      XX(1)=X1
      XX(2)=X2
      XX(3)=X3
      DO 5 I=1,2
      IP1=I+1
      DO 5 J=IP1,3
      IF(XX(I).LE.XX(J)) GO TO 5
      W=XX(I)
      XX(I)=XX(J)
      XX(J)=W
    5 CONTINUE
      KOSUU=3
      GO TO 8
C                                  1-JIKKON & 2-KYOKON
    2 DSQR=DSQRT(D)
    3 U3=.5*(-Q+DSQR)
      V3=.5*(-Q-DSQR)
      CUB=1./3.
      AL=0.0
      BE=0.0
      IF(U3.GT.0.) AL=U3**CUB
      IF(U3.LT.0.) AL=-(DABS(U3))**CUB
      IF(V3.GT.0.) BE=V3**CUB
      IF(V3.LT.0.) BE=-(DABS(V3))**CUB
      X1=AL+BE-A4
      X2=-.5*(AL+BE)-A4
      Y2=.5*DSQRT(3.0D0)*DABS((AL-BE))
      XX(1)=X1
      KOSUU=1
      GO TO 8
C                                       KOTAE NO SENBETSU
    8 DUMMY=0.0
      KKOS=0
      DO 6 I=1,KOSUU
      XX2(I)=0.0
      IF(XX(I).LT.XANMIN.OR.XX(I).GT.XANMAX) GO TO 6
      KKOS=KKOS+1
      XX2(KKOS)=XX(I)
    6 CONTINUE
C ======================================================================
CC    WRITE(7,710) XANS,KOSUU,(XX(I),I=1,3),KKOS,(XX2(I),I=1,3)
C ======================================================================
CC    IF(KKOS.GE.1) XANS=XX2(KKOS) ............ '89  4/ 3  DUMMY
      IF(KKOS.EQ.1) XANS=XX2(1)
      IF(KKOS.LE.0.AND.KANCOD.LE.-1) GO TO 12
      IF(KKOS.LE.0.AND.KANCOD.GE.1)  GO TO 18
      IF(KKOS.EQ.1) GO TO 10
      KKOSM1=KKOS-1
      DO 30 I=1,KKOSM1
      IP1=I+1
      DO 30 J=2,KKOS
      IF(XX2(I).LT.XX2(J)) GO TO 30
      A1=XX2(I)
      XX2(I)=XX2(J)
      XX2(J)=A1
   30 CONTINUE
      XANS=XX2(1)
      GO TO 10
C                                                 KANCOD =< -1
   12 DUMMY=0.0
      DO 14 I=1,KOSUU
      J=KOSUU+1-I
      IF(XX(J).GT.XANMIN) GO TO 14
      GO TO 16
   14 CONTINUE
   16 IF(XX(J).LE.XANMIN) XANS=XX(J)
      GO TO 10
C                                                 KANCOD >= 1
   18 DUMMY=0.0
      DO 20 I=1,KOSUU
      J=I
      IF(XX(J).LT.XANMAX) GO TO 20
      GO TO 22
   20 CONTINUE
   22 IF(XX(J).GE.XANMAX) XANS=XX(J)
      GO TO 10
C
   10 DUMMY=0.0
C ======================================================================
CC    WRITE(7,710) XANS,KOSUU,(XX(I),I=1,3),KKOS,(XX2(I),I=1,3)
CC710 FORMAT(1H ,'  XANS,KOSUU,XX,KKOS,XX2   :',F10.3,2(I5,3F10.3))
C ======================================================================
      RETURN
      END
C *** SUBETE ONAJI ATAI TO SURU ***
      SUBROUTINE ALLONA(A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10)
      real*8 A1,A2,A3,A4,A5,A6,A7,A8,A9,A10
      A1 =A0
      A2 =A0
      A3 =A0
      A4 =A0
      A5 =A0
      A6 =A0
      A7 =A0
      A8 =A0
      A9 =A0
      A10=A0
      RETURN
      END
