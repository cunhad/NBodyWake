      PROGRAM BinomialFactorial_FMPack

      USE FMZM

      IMPLICIT NONE
!             (IM) for multiple precision integer

      TYPE (FM), SAVE :: X1, X2, X3, X4
      TYPE (IM), SAVE :: I1, I2, I3
      CHARACTER(80)  :: ST1

!                 Set precision to give at least 60 significant digits.

!      CALL FM_SET(60)

!      X1 = TO_FM('3.12')
!      CALL FM_FORM('F65.60',X1,ST1)
!      WRITE (*   ,"(/I10,4X,A)") 0,TRIM(ST1)
 
!     I1 = TO_IM('600000000')
!      I2 = TO_IM('600000')
!      I3 = BINOMIAL(I1,I2)
!      CALL IM_FORM('I65',I3,ST1)
!      WRITE (*   ,"(' p =',A)") 0,TRIM(ST1)
!      CALL IM_PRINT(I3)        
!      CALL IM_FORM('I65',I3,ST1)
!      WRITE (*   ,"(' p =',A)") TRIM(ST1)

      
      CALL FM_SET(600)

      I1 = TO_IM('6000000000')
!      I1 = TO_IM('6000')
      I2 = TO_IM('6000000')

      X1 = (TO_FM('0.0026'))**(I2)
      X2 = (1-TO_FM('0.0026'))**(I1-I2)


      I3 = BINOMIAL(I1,I2)

!      X3=X1*X2*I3
      X3=I3
!      CALL IM_PRINT(I3)
      CALL FM_PRINT(X3)
      


      END PROGRAM BinomialFactorial_FMPack
