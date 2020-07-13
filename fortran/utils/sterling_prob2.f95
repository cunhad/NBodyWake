      PROGRAM BinomialFactorial_FMPack

      USE FMZM

      IMPLICIT NONE
!             (IM) for multiple precision integer

      TYPE (FM),SAVE :: P_W,P_NW,PROB1,PROB2,STERLING,PI,FN_VOL,FN_OUTL,TEST
        TYPE (FM),SAVE :: SUM_W,SUM_NW 
      TYPE (IM),SAVE :: N_VOL,N_OUTL,BINOM
!      LOGICAL :: C1
      INTEGER :: ITER,K,NVOL,NOUTL,KOUT, NERROR
        LOGICAL, EXTERNAL :: FMCOMPARE
!      CHARACTER(80)  :: ST1

!      KOUT = 18
!      OPEN (KOUT,FILE='Prob_out_NW_ZA_NN.out')

!      NERROR = 0
!      CALL FM_SETVAR(' KW = 18 ')



!                 Set precision to give at least 60 significant digits.
      
!      CALL FM_SET(60)
      CALL FM_SET(4000)

!      I1 = TO_IM('6000000000')
!      I2 = TO_IM('6000000')
!      NVOL      = 400
!      N_VOL      = TO_IM(NVOL)
       N_VOL      = TO_IM('450000000') !from 60000000 on starts to be slow
!      N_VOL       = TO_IM('90000000') !from 60000000 on starts to be slow
!      N_VOL      = TO_IM('50000') !from 60000000 on starts to be slow
!      N_OUTL     = TO_IM('600000')
!      N_OUTL      = TO_IM('450000')
!      N_OUTL     = TO_IM('600000')
      FN_VOL     = TO_FM('450000000')           
!      FN_VOL     = TO_FM('90000000')
!      FN_VOL     = TO_FM('50000')           
!      FN_OUTL    = TO_FM('600000')
!      FN_OUTL    = TO_FM('450000')
!      FN_OUTL    = TO_FM('600000')
!      BINOM      = TO_IM('20')
       P_W     = TO_FM('0.00299479')
!       P_NW    = TO_FM('0.001041666')  
!       P_W     = TO_FM('0.004904514')    !z=2.9
       P_NW    = TO_FM('0.001041666')    !z=2.9


        SUM_W   = TO_FM('0')
        SUM_NW   = TO_FM('0')



!      IF (I1 > 300) THEN
!         CALL IM_PRINT(I1)
!      ENDIF


      CALL FM_PI(PI)


      K = 1400000
!      K = 490000
!      K = 96000
!      K = 100
!      K = 1 
      DO ITER = 0,K
!      DO ITER = 400000,480000
        N_OUTL      = TO_IM(ITER)
        FN_OUTL     = TO_FM(ITER) 
!        BINOM = BINOMIAL(N_VOL,N_OUTL)        
        BINOM = BINOMIAL(FN_VOL,FN_OUTL)        
!      BINOM = SQRT(FN_VOL/(2*PI*(FN_VOL-FN_OUTL)*FN_OUTL))&
!        *((FN_VOL/(FN_VOL-FN_OUTL))**FN_VOL)&
!               *(((FN_VOL-FN_OUTL)/FN_OUTL)**FN_OUTL)
        PROB2   = (P_W**FN_OUTL)*((1-P_W)**(FN_VOL-FN_OUTL))*BINOM
        PROB1   = (P_NW**FN_OUTL)*((1-P_NW)**(FN_VOL-FN_OUTL))*BINOM
!        PROB2   = (P_W**FN_OUTL)*((1-P_W)**(FN_VOL-FN_OUTL))*STERLING
!        TEST    = PROB1/PROB2
!        SUM_W=SUM_W+PROB2
!        SUM_NW=SUM_NW+PROB1

                IF (FMCOMPARE(TO_FM(1)-SUM_NW,'<',SUM_W)) THEN
!                        IF (ITER > 0) THEN
                                CALL IM_PRINT(N_OUTL)
                                CALL FM_PRINT(SUM_NW)
                                CALL FM_PRINT(SUM_W)                   
!                                EXIT
!                        END IF
                END IF

        SUM_W=SUM_W+PROB2
        SUM_NW=SUM_NW+PROB1
                IF (FMCOMPARE(TO_FM(1)-SUM_NW,'<',SUM_W)) THEN
!                        IF (ITER > 0) THEN
                                CALL IM_PRINT(N_OUTL)
                                CALL FM_PRINT(SUM_NW)
                                CALL FM_PRINT(SUM_W)
                                EXIT
!                        END IF
                END IF

!        CALL IM_PRINT(BINOM)
!        CALL FM_PRINT(STERLING)
!        CALL FM_PRINT(1-TEST)
!        CALL IM_PRINT(N_OUTL)
!        CALL FM_PRINT(1-SUM_W)
!     CALL FM_PRINT(PROB1)

      ENDDO  


!      I3 = BINOMIAL(I1,I2)

!      X3=X1*X2*I3
!      X3=I3
!      CALL IM_PRINT(I3)
!      CALL FM_PRINT(X3)






      END PROGRAM BinomialFactorial_FMPack
