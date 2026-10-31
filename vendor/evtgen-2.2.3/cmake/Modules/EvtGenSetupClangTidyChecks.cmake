
########################################################################
# Copyright 1998-2022 CERN for the benefit of the EvtGen authors       #
#                                                                      #
# This file is part of EvtGen.                                         #
#                                                                      #
# EvtGen is free software: you can redistribute it and/or modify       #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# EvtGen is distributed in the hope that it will be useful,            #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with EvtGen.  If not, see <https://www.gnu.org/licenses/>.     #
########################################################################


if (EVTGEN_RUN_CLANG_TIDY)
    find_program( CLANG_TIDY_PATH NAMES clang-tidy )
    if (${CLANG_TIDY_PATH} STREQUAL CLANG_TIDY_PATH-NOTFOUND)
        message(WARNING "EvtGen: clang-tidy not found, disabling checks during build")
        unset(CMAKE_CXX_CLANG_TIDY)
    else()
        # NB have to have these two separate set commands (rather than doing one before the if and then doing a list append) because of the FORCE
        if (EVTGEN_CLANG_TIDY_FIX)
            set(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY_PATH} --checks=-*,${EVTGEN_CLANG_TIDY_CHECKS} --fix CACHE STRING "Used to set the CXX_CLANG_TIDY property on each target" FORCE)
        else()
            set(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY_PATH} --checks=-*,${EVTGEN_CLANG_TIDY_CHECKS} CACHE STRING "Used to set the CXX_CLANG_TIDY property on each target" FORCE)
        endif()
    endif()
else()
    unset(CMAKE_CXX_CLANG_TIDY)
endif()
