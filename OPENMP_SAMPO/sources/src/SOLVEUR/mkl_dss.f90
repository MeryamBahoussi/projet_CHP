!*******************************************************************************
!                              INTEL CONFIDENTIAL
!   Copyright(C) 2004-2005 Intel Corporation. All Rights Reserved.
!   The source code contained  or  described herein and all documents related to
!   the source code ("Material") are owned by Intel Corporation or its suppliers
!   or licensors.  Title to the  Material remains with  Intel Corporation or its
!   suppliers and licensors. The Material contains trade secrets and proprietary
!   and  confidential  information of  Intel or its suppliers and licensors. The
!   Material  is  protected  by  worldwide  copyright  and trade secret laws and
!   treaty  provisions. No part of the Material may be used, copied, reproduced,
!   modified, published, uploaded, posted, transm1itted, distributed or disclosed
!   in any way without Intel's prior express written permission.
!   No license  under any  patent, copyright, trade secret or other intellectual
!   property right is granted to or conferred upon you by disclosure or delivery
!   of the Materials,  either expressly, by implication, inducement, estoppel or
!   otherwise.  Any  license  under  such  intellectual property  rights must be
!   express and approved by Intel in writing.
!
!*******************************************************************************
!   Content : MKL DSS Fortran-90 header file
!
!           Contains main datatypes, routines and constants definition.
!           For CDECL use only.
!
!*******************************************************************************
!DEC$ IF .NOT. DEFINED( __MKL_DSS_F90 )

!DEC$ DEFINE __MKL_DSS_F90 


MODULE MKL_DSS_PRIVATE
  TYPE MKL_DSS_HANDLE
     INTEGER(KIND=8) DUMMY
  END TYPE MKL_DSS_HANDLE
END MODULE MKL_DSS_PRIVATE

MODULE MKL_DSS

  USE MKL_DSS_PRIVATE

  INTEGER, PARAMETER :: MKL_DSS_DEFAULTS = 0

  !
  ! Message level option definitions
  !

  INTEGER(KIND=4), PARAMETER :: MKL_DSS_MSG_LVL_SUCCESS = -2147483647
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_MSG_LVL_DEBUG   = -2147483646
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_MSG_LVL_INFO    = -2147483645
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_MSG_LVL_WARNING = -2147483644
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_MSG_LVL_ERROR   = -2147483643
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_MSG_LVL_FATAL   = -2147483642

  !
  ! Termination level option definitions
  !

  INTEGER(KIND=4), PARAMETER :: MKL_DSS_TERM_LVL_SUCCESS = 1073741832
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_TERM_LVL_DEBUG   = 1073741840
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_TERM_LVL_INFO    = 1073741848
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_TERM_LVL_WARNING = 1073741856
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_TERM_LVL_ERROR   = 1073741864
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_TERM_LVL_FATAL   = 1073741872

  !
  ! Structure option definitions
  !

  INTEGER(KIND=4), PARAMETER :: MKL_DSS_SYMMETRIC           = 536870976
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_SYMMETRIC_STRUCTURE = 536871040
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_NON_SYMMETRIC       = 536871104

  !
  ! Reordering option definitions
  !

  INTEGER(KIND=4), PARAMETER :: MKL_DSS_AUTO_ORDER    = 268435520
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_MY_ORDER      = 268435584
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_OPTION1_ORDER = 268435648

  !
  ! Factorization option definitions
  !

  INTEGER(KIND=4), PARAMETER :: MKL_DSS_POSITIVE_DEFINITE           = 134217792
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_INDEFINITE                  = 134217856
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_HERMITIAN_POSITIVE_DEFINITE = 134217920
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_HERMITIAN_INDEFINITE        = 134217984

  !
  ! Return status values
  !

  INTEGER(KIND=4), PARAMETER :: MKL_DSS_SUCCESS         = 0
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_ZERO_PIVOT      = -1
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_OUT_OF_MEMORY   = -2
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_FAILURE         = -3
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_ROW_ERR         = -4
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_COL_ERR         = -5
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_TOO_FEW_VALUES  = -6
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_TOO_MANY_VALUES = -7
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_NOT_SQUARE      = -8
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_STATE_ERR       = -9
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_INVALID_OPTION  = -10
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_OPTION_CONFLICT = -11
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_MSG_LVL_ERR     = -12
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_TERM_LVL_ERR    = -13
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_STRUCTURE_ERR   = -14
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_REORDER_ERR     = -15
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_VALUES_ERR      = -16
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_STATISTICS_INVALID_MATRIX = -17
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_STATISTICS_INVALID_STATE  = -18
  INTEGER(KIND=4), PARAMETER :: MKL_DSS_STATISTICS_INVALID_STRING = -19

  !
  ! Function prototypes for DSS routines
  !

  INTERFACE
     FUNCTION DSS_CREATE( HANDLE, OPT )
       USE MKL_DSS_PRIVATE
       TYPE(MKL_DSS_HANDLE), INTENT(OUT)   :: HANDLE
       INTEGER(KIND=4),       INTENT(IN)    :: OPT
       INTEGER(KIND=4)                      :: DSS_CREATE
     END FUNCTION DSS_CREATE
  END INTERFACE

  INTERFACE
     FUNCTION DSS_DEFINE_STRUCTURE( HANDLE, OPT, ROWINDEX, NROWS, RCOLS, COLUMNS, NNONZEROS )
       USE MKL_DSS_PRIVATE
       TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
       INTEGER(KIND=4),       INTENT(IN)    :: OPT
       INTEGER(KIND=4),       INTENT(IN)    :: NROWS
       INTEGER(KIND=4),       INTENT(IN)    :: RCOLS
       INTEGER(KIND=4),       INTENT(IN)    :: NNONZEROS
       INTEGER(KIND=4),       INTENT(IN)    :: ROWINDEX( * ) ! * = MIN(NROWS, NCOLS)+1
       INTEGER(KIND=4),       INTENT(IN)    :: COLUMNS( * ) ! * = NNONZEROS
       INTEGER(KIND=4)                      :: DSS_DEFINE_STRUCTURE
     END FUNCTION DSS_DEFINE_STRUCTURE
  END INTERFACE

  INTERFACE
     FUNCTION DSS_REORDER( HANDLE, OPT, PERM )
       USE MKL_DSS_PRIVATE
       TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
       INTEGER(KIND=4),       INTENT(IN)    :: OPT
       INTEGER(KIND=4),       INTENT(IN)    :: PERM( * )
       INTEGER(KIND=4)                      :: DSS_REORDER
     END FUNCTION DSS_REORDER
  END INTERFACE

  INTERFACE
     FUNCTION DSS_FACTOR_REAL( HANDLE, OPT, RVALUES )
       USE MKL_DSS_PRIVATE
       TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
       INTEGER(KIND=4),       INTENT(IN)    :: OPT
       REAL(KIND=8),          INTENT(IN)    :: RVALUES( * )
       INTEGER(KIND=4)                      :: DSS_FACTOR_REAL
     END FUNCTION DSS_FACTOR_REAL
  END INTERFACE

  INTERFACE
     FUNCTION DSS_FACTOR_COMPLEX( HANDLE, OPT, RVALUES )
       USE MKL_DSS_PRIVATE
       TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
       INTEGER(KIND=4),       INTENT(IN)    :: OPT
       COMPLEX(KIND=8),       INTENT(IN)    :: RVALUES( * )
       INTEGER(KIND=4)                      :: DSS_FACTOR_COMPLEX
     END FUNCTION DSS_FACTOR_COMPLEX
  END INTERFACE

  INTERFACE
     FUNCTION DSS_SOLVE_REAL( HANDLE, OPT, RRHSVALUES, NRHS, RSOLVALUES )
       USE MKL_DSS_PRIVATE
       TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
       INTEGER(KIND=4),       INTENT(IN)    :: OPT
       INTEGER(KIND=4),       INTENT(IN)    :: NRHS
       REAL(KIND=8),          INTENT(IN)    :: RRHSVALUES( * )
       REAL(KIND=8),          INTENT(OUT)   :: RSOLVALUES( * )
       INTEGER(KIND=4)                      :: DSS_SOLVE_REAL
     END FUNCTION DSS_SOLVE_REAL
  END INTERFACE

  INTERFACE
     FUNCTION DSS_SOLVE_COMPLEX( HANDLE, OPT, RRHSVALUES, NRHS, RSOLVALUES )
       USE MKL_DSS_PRIVATE
       TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
       INTEGER(KIND=4),       INTENT(IN)    :: OPT
       INTEGER(KIND=4),       INTENT(IN)    :: NRHS
       COMPLEX(KIND=8),       INTENT(IN)    :: RRHSVALUES( * )
       COMPLEX(KIND=8),       INTENT(OUT)   :: RSOLVALUES( * )
       INTEGER(KIND=4)                      :: DSS_SOLVE_COMPLEX
     END FUNCTION DSS_SOLVE_COMPLEX
  END INTERFACE

  INTERFACE
     FUNCTION DSS_DELETE( HANDLE, OPT )
       USE MKL_DSS_PRIVATE
       TYPE(MKL_DSS_HANDLE), INTENT(IN)    :: HANDLE
       INTEGER(KIND=4),       INTENT(IN)    :: OPT
       INTEGER(KIND=4)                      :: DSS_DELETE
     END FUNCTION DSS_DELETE
  END INTERFACE

  INTERFACE
     FUNCTION DSS_STATISTICS( HANDLE, OPT, STAT, RET )
       USE MKL_DSS_PRIVATE
       TYPE(MKL_DSS_HANDLE), INTENT(IN)    :: HANDLE
       INTEGER(KIND=4),       INTENT(IN)    :: OPT
       INTEGER(KIND=4),       INTENT(IN)    :: STAT( * )
       REAL(KIND=8),          INTENT(OUT)   :: RET( * )
       INTEGER(KIND=4)                      :: DSS_STATISTICS
     END FUNCTION DSS_STATISTICS
  END INTERFACE

  INTERFACE
     SUBROUTINE MKL_CVT_TO_NULL_TERMINATED_STR( DESTSTR, DESTLEN, SRCSTR )
       INTEGER(KIND=4),       INTENT(OUT)   :: DESTSTR( * )
       INTEGER(KIND=4),       INTENT(IN)    :: DESTLEN
       CHARACTER(*),          INTENT(IN)    :: SRCSTR
     END SUBROUTINE MKL_CVT_TO_NULL_TERMINATED_STR
  END INTERFACE

END MODULE MKL_DSS


!DEC$ ENDIF

