        PROGRAM PRINT_FUNC

        USE FMZM

        IMPLICIT NONE


        TYPE (FM):: F,X

        CALL FM_SET(800)

        X=TO_FM('50')
        F=ERF(X/sqrt(TO_FM('2')))

        CALL FM_PRINT(F)
        CALL FM_PRINT(TO_FM(1)-F)
        

        END PROGRAM PRINT_FUNC
