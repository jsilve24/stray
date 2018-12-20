library(testthat)
library(mongrel)

#Sys.setenv(KMP_DUPLICATE_LIB_OK="TRUE")
test_check("mongrel")
#Sys.unsetenv("KMP_DUPLICATE_LIB_OK")
