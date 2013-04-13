# This module finds an installed PLPLOT package.
#
# It sets the following variables:
#  PLPLOT_FOUND         - Set to false, or undefined, if vigra isn't found.
#  PLPLOT_INCLUDE_DIR   - The CppUnit include directory.
#  PLPLOT_LIBRARIES     - the PLPLOT library

FIND_PATH(PLPLOT_INCLUDE_DIR
          NAMES plplot.h
          PATHS /opt/local/include/plplot
                /usr/include/plplot
                /export/home/tobinder/local/include/plplot)
FIND_LIBRARY(PLPLOT_LIBRARY1
          NAMES plplotcxxd
          PATHS /opt/local/lib
                /usr/local/lib
                /export/home/tobinder/local/lib)
FIND_LIBRARY(PLPLOT_LIBRARY2
          NAMES plplotd
          PATHS /opt/local/lib
                /usr/local/lib
                /export/home/tobinder/local/lib)
FIND_LIBRARY(PLPLOT_LIBRARY3
          NAMES ltdl
          PATHS /opt/local/lib
                /usr/local/lib
                /export/home/tobinder/local/lib)
FIND_LIBRARY(PLPLOT_LIBRARY4
          NAMES dl
          PATHS /opt/local/lib
                /usr/local/lib
                /export/home/tobinder/local/lib)
FIND_LIBRARY(PLPLOT_LIBRARY5
          NAMES m
          PATHS /opt/local/lib
                /usr/local/lib
                /export/home/tobinder/local/lib)
FIND_LIBRARY(PLPLOT_LIBRARY6
          NAMES csirocsa
          PATHS /opt/local/lib
                /usr/local/lib
                /export/home/tobinder/local/lib)
FIND_LIBRARY(PLPLOT_LIBRARY7
          NAMES freetype
          PATHS /opt/local/lib
                /usr/local/lib
                /export/home/tobinder/local/lib)

IF (PLPLOT_INCLUDE_DIR AND PLPLOT_LIBRARY1 AND PLPLOT_LIBRARY2 AND PLPLOT_LIBRARY3 AND PLPLOT_LIBRARY4 AND
    PLPLOT_LIBRARY5 AND PLPLOT_LIBRARY6 AND PLPLOT_LIBRARY7)
    SET(PLPLOT_FOUND TRUE)
    SET(PLPLOT_LIBRARIES ${PLPLOT_LIBRARY1} ${PLPLOT_LIBRARY2} ${PLPLOT_LIBRARY3} ${PLPLOT_LIBRARY4}
                         ${PLPLOT_LIBRARY5} ${PLPLOT_LIBRARY6} ${PLPLOT_LIBRARY7})
ENDIF (PLPLOT_INCLUDE_DIR AND PLPLOT_LIBRARY1 AND PLPLOT_LIBRARY2 AND PLPLOT_LIBRARY3 AND PLPLOT_LIBRARY4 AND
       PLPLOT_LIBRARY5 AND PLPLOT_LIBRARY6 AND PLPLOT_LIBRARY7)

IF (NOT PLPLOT_FOUND)
    MESSAGE(FATAL_ERROR "Could not find PLPLOT")
ENDIF (NOT PLPLOT_FOUND)
