library(testthat)
library(stray)

#Sys.setenv(KMP_DUPLICATE_LIB_OK="TRUE")
test_check("stray")
#Sys.unsetenv("KMP_DUPLICATE_LIB_OK")
