# This module finds an installed GSL package.
#
# It sets the following variables:
#  GSL_FOUND         - Set to false, or undefined, if vigra isn't found.
#  GSL_INCLUDE_DIR   - The CppUnit include directory.
#  GSL_LIBRARIES     - the GSL library

FIND_PATH(GSL_INCLUDE_DIR
          NAMES gsl_blas.h
          PATHS /opt/local/include/gsl
                /usr/include/gsl)
FIND_LIBRARY(GSL_LIBRARY
          NAMES gsl
          PATHS /opt/local/lib
                /usr/local/lib
                /export/home/tobinder/local/lib)
FIND_LIBRARY(GSL_CBLAS_LIBRARY
          NAMES gslcblas
          PATHS /opt/local/lib
                /usr/local/lib
                /export/home/tobinder/local/lib)

IF (GSL_INCLUDE_DIR AND GSL_LIBRARY AND GSL_CBLAS_LIBRARY)
    SET(GSL_FOUND TRUE)
    SET(GSL_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
ENDIF (GSL_INCLUDE_DIR AND GSL_LIBRARY AND GSL_CBLAS_LIBRARY)

IF (NOT GSL_FOUND)
    MESSAGE(FATAL_ERROR "Could not find GSL")
ENDIF (NOT GSL_FOUND)
