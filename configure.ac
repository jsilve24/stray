# autoconf file inspiration from 
# https://github.com/USCbiostats/software-dev/wiki/Setting-up-optional-OpenMP-support

AC_PREREQ([2.69])
AC_INIT(stray, m4_esyscmd_s([awk '/^Version:/ {print $2}' DESCRIPTION]))

## Set R_HOME, respecting an environment variable if one is set 
: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])   
fi
## Use R to set CXX and CXXFLAGS
CXX=$(${R_HOME}/bin/R CMD config CXX)
CXXFLAGS=$("${R_HOME}/bin/R" CMD config CXXFLAGS)
CXX11=$(${R_HOME}/bin/R CMD config CXX11)
CXX11FLAGS=$("${R_HOME}/bin/R" CMD config CXX11FLAGS)

## We are using C++
AC_LANG(C++)
AC_REQUIRE_CPP

## below pilfered from configure.ac in RcppArmadillo
## Check for Apple LLVM

AC_MSG_CHECKING([for macOS])
RSysinfoName=$("${R_HOME}/bin/Rscript" --vanilla -e 'cat(Sys.info()[["sysname"]])')

if test x"${RSysinfoName}" == x"Darwin"; then
    AC_MSG_RESULT([found])
    AC_MSG_CHECKING([for macOS Apple compiler])

    apple_compiler=$($CXX --version 2>&1 | grep -i -c -e 'apple llvm')

    if test x"${apple_compiler}" == x"1"; then
        AC_MSG_RESULT([found])
        AC_MSG_WARN([OpenMP unavailable and turned off.])
        can_use_openmp="no"
    else
        AC_MSG_RESULT([not found])
        AC_MSG_CHECKING([for clang compiler])
        clang_compiler=$($CXX --version 2>&1 | grep -i -c -e 'clang ')

        if test x"${clang_compiler}" == x"1"; then
            AC_MSG_RESULT([found])
            AC_MSG_CHECKING([for OpenMP compatible version of clang])
            clang_version=$(${CXX} -v 2>&1 | awk '/^.*clang version/ {print $3}')

            case ${clang_version} in
                4.*|5.*|6.*|7.*|8.*)
                    AC_MSG_RESULT([found and suitable])
                    can_use_openmp="yes"
                ;;
                *)
                    AC_MSG_RESULT([not found])
                    AC_MSG_WARN([OpenMP unavailable and turned off.])
                    can_use_openmp="no"
                ;;
            esac
        else
            AC_MSG_RESULT([not found])
            AC_MSG_WARN([unsupported macOS build detected; if anything breaks, you keep the pieces.])
        fi
    fi
fi


## Default the OpenMP flag to the empty string.
## If and only if OpenMP is found, expand to $(SHLIB_OPENMP_CXXFLAGS)
openmp_flag=""
## Set the fallback, by default it is nope
have_openmp="FALSE"

if test x"${can_use_openmp}" == x"yes"; then
    AC_MSG_CHECKING([for OpenMP])
    ## if R has -fopenmp we should be good
    allldflags=$(${R_HOME}/bin/R CMD config --ldflags)
    hasOpenMP=$(echo ${allldflags} | grep -- -fopenmp)
    if test x"${hasOpenMP}" == x""; then
        AC_MSG_RESULT([missing])
        have_openmp="FALSE"
    else   
        AC_MSG_RESULT([found and suitable])
        have_openmp="TRUE"
        openmp_flag='$(SHLIB_OPENMP_CXXFLAGS)'
    fi
fi


## now use all these
AC_SUBST([OPENMP_FLAG], ["${openmp_flag}"])
AC_SUBST([OPENMP_TF], ["${have_openmp}"])
AC_CONFIG_FILES([src/Makevars R/open_mp_status.R])
AC_OUTPUT
