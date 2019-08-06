matrix_maintain_dim <- function(x){
  if (is.vector(x)) return(matrix(x))
  return(x)
}


#' Log-Ratio transforms for orthus objects
#' @param x orthus data array (e.g., first s rows are multinomial parameters or log-ratios)
#' @param y orthus data array (e.g., first s rows are multinomial parameters or log-ratios)
#' @param s first s rows of x are transformed
#' @param V transformation matrix (defines transform) 
#' @param d for ALR, which component (integer position) to take as reference
#' (default is ncol(x)) for alrInv corresponds to column position in untransformed
#' matrix.
# @param inv for ALR and CLR, transformation matrix is different forward and inverse
# @param D the number of parts (e.g., number of columns in untransformed data)
#' @name orthus_lr_transforms
NULL

#' @rdname orthus_lr_transforms
#' @export
oglr <- function(x,s, V){
  x.star <- glr_array(x[1:s,,,drop=F],V, parts=1)
  d.star <- dim(x.star)[1]
  n <-  dim(x)[1] + d.star - s
  y <- array(0, dim=c(n, dim(x)[2], dim(x)[3]))
  y[1:d.star,,] <- x.star
  y[(d.star+1):n,,] <- x[(s+1):dim(x)[1],,,drop=F]
  return(y)
}

#' @rdname orthus_lr_transforms
#' @export
oglrInv <- function(x, s, V){
  x.star <- glrInv_array(x[1:s,,,drop=F],V, coords=1)
  d.star <- dim(x.star)[1]
  n <-  dim(x)[1] + d.star - s
  y <- array(0, dim=c(n, dim(x)[2], dim(x)[3]))
  y[1:d.star,,] <- x.star
  y[(d.star+1):n,,] <- x[(s+1):dim(x)[1],,,drop=F]
  return(y)
}

#' @rdname orthus_lr_transforms
#' @export
oalr <- function(x, s, d=NULL){
  added_dim <- FALSE
  if (length(dim(x))==2) {x <- add_array_dim(x,3); added_dim=TRUE}
  if (is.null(d)) d <- s
  B <- create_alr_base(s, d, inv=FALSE)
  y <- oglr(x, s, B)
  if (added_dim) return(matrix_maintain_dim(y[,,1]))
  y
}

#' @rdname orthus_lr_transforms
#' @export
oalrInv <- function(y, s, d=NULL){
  added_dim <- FALSE
  if (length(dim(y))==2) {y <- add_array_dim(y,3); added_dim=TRUE}
  if (is.null(d)) d <- s+1
  B <- create_alr_base(s+1, d, inv=TRUE)
  x <- oglrInv(y, s, B)
  if (added_dim) return(matrix_maintain_dim(x[,,1]))
  x
}

#' @rdname orthus_lr_transforms
#' @export
oilr <- function(x, s, V=NULL){
  added_dim <- FALSE
  if (length(dim(x))==2) { x <- add_array_dim(x,3); added_dim=TRUE}
  if (is.null(V)) V <- create_default_ilr_base(s)
  y <- oglr(x, s, V)
  if (added_dim) return(matrix_maintain_dim(y[,,1]))
  y
}

#' @rdname orthus_lr_transforms
#' @export
oilrInv <- function(y, s, V=NULL){
  added_dim <- FALSE
  if (length(dim(y))==2) {y <- add_array_dim(y,3); added_dim=TRUE}
  if (is.null(V)) V <- create_default_ilr_base(s+1)
  x <- oglrInv(y, s, V)
  if (added_dim) return(matrix_maintain_dim(x[,,1]))
  x
}

#' @rdname orthus_lr_transforms
#' @export
oclr <- function(x, s){
  added_dim <- FALSE
  if (length(dim(x))==2) {x <- add_array_dim(x,3); added_dim=TRUE}
  y <- oglr(x, s, create_clr_base(s))
  if (added_dim) return(matrix_maintain_dim(y[,,1]))
  y
}

#' @rdname orthus_lr_transforms
#' @export
oclrInv <- function(x, s){
  added_dim <- FALSE
  if (length(dim(x))==2) {x <- add_array_dim(x,3); added_dim=TRUE}
  x.star <- clrInv_array(x[1:s,,], coords=1)
  d.star <- dim(x.star)[1]
  n <-  dim(x)[1] + d.star - s
  y <- array(0, dim=c(n, dim(x)[2], dim(x)[3]))
  y[1:d.star,,] <- x.star
  y[(d.star+1):n,,] <- x[(s+1):dim(x)[1],,]
  if (added_dim) return(matrix_maintain_dim(y[,,1]))
  y
}


# Covariance Matricies ----------------------------------------------------


#' Convert orthus covariance matricies between representations
#'
#' @param Sigma covariance matrix arrat in specified transformed space 
#'   (dim(Sigma)[3]=iter)
#' @param s first s rows and colums of Sigma are transformed
#' @param V ILR contrast matrix (i.e., transformation matrix of ILR)
#' @param V1 ILR contrast matrix of basis Sigma is already in
#' @param V2 ILR contrast matrix of basis Sigma is desired in
#' @param d1 alr reference element Sigma is already expressed with respec to
#' @param d2 alr reference element Sigma is to be expressed with respect to
#'
#' @return matrix
#' @name convert_orthus_covariance
NULL

#' @rdname convert_orthus_covariance
#' @export
oilrvar2ilrvar <- function(Sigma, s, V1, V2){
  added_dim <- FALSE
  if (length(dim(Sigma))==2) {Sigma <- add_array_dim(Sigma,3); added_dim=TRUE}
  for (i in 1:dim(Sigma)[3]){
    one <- 1:s; two <- (s+1):dim(Sigma)[1]
    Sigma[1:s, 1:s, i] <-  t(V2) %*% V1 %*% Sigma[one,one,i] %*% t(V1) %*% V2
    Sigma[one,two,i] <- t(V2) %*% V1 %*% Sigma[one,two,i]
    Sigma[two,one,i] <- t(Sigma[one,two,i])
  }
  if (added_dim) return(matrix_maintain_dim(Sigma[,,i]))
  return(Sigma)
}

#' @rdname convert_orthus_covariance
#' @export
oilrvar2clrvar <- function(Sigma, s, V){
  added_dim <- FALSE
  if (length(dim(Sigma))==2) {Sigma <- add_array_dim(Sigma,3); added_dim=TRUE}
  d <- dim(Sigma)
  d[1] <-  d[2]<- d[1]+1
  O <- array(0, dim=d)
  for (i in 1:dim(Sigma)[3]){
    one <- 1:s; two <- (s+1):dim(Sigma)[1]
    O[1:(s+1), 1:(s+1), i] <- V %*% Sigma[one,one,i] %*% t(V)
    O[1:(s+1), two+1,i] <-   V %*% Sigma[one,two,i]
    O[two+1,1:(s+1),i] <- t(O[1:(s+1), two+1,i])
    O[two+1,two+1,i] <- Sigma[two,two,i]
  }
  if (added_dim) return(matrix_maintain_dim(O[,,1]))
  return(O)
}

#' @rdname convert_orthus_covariance
#' @export
oclrvar2ilrvar <- function(Sigma, s, V){
  added_dim <- FALSE
  if (length(dim(Sigma))==2) {Sigma <- add_array_dim(Sigma,3); added_dim=TRUE}
  d <- dim(Sigma)
  d[1] <- d[2] <- d[1]-1
  O <- array(0, dim=d)
  for (i in 1:dim(Sigma)[3]){
    one <- 1:s; two <- (s+1):dim(Sigma)[1]
    O[1:(s-1), 1:(s-1), i] <- t(V) %*% Sigma[one,one,i] %*% V
    O[1:(s-1),two-1,i] <-   t(V) %*% Sigma[one,two,i]
    O[two-1,1:(s-1),i] <- t(O[1:(s-1),two-1,i])
    O[two-1,two-1,i] <- Sigma[two,two,i]
  }
  if (added_dim) return(matrix_maintain_dim(O[,,1]))
  return(O)
}

#' @rdname convert_orthus_covariance
#' @export
oalrvar2clrvar <- function(Sigma, s, d1){
  added_dim <- FALSE
  if (length(dim(Sigma))==2) {Sigma <- add_array_dim(Sigma,3); added_dim=TRUE}
  d <- dim(Sigma)
  d[1] <- d[2] <- d[1]+1
  O <- array(0, dim=d)
  G1 <- driver::create_alr_base(s+1, d1, inv=TRUE) - 1/(s+1)
  for (i in 1:dim(Sigma)[3]){
    one <- 1:s; two <- (s+1):dim(Sigma)[1]
    O[1:(s+1), 1:(s+1), i] <- G1 %*% Sigma[one,one,i] %*% t(G1)
    O[1:(s+1), two+1,i] <-   G1 %*% Sigma[one,two,i]
    O[two+1,1:(s+1),i] <- t(O[1:(s+1), two+1,i])
    O[two+1,two+1,i] <- Sigma[two,two,i]
  }
  if (added_dim) return(matrix_maintain_dim(O[,,1]))
  return(O)
}

#' @rdname convert_orthus_covariance
#' @export
oclrvar2alrvar <- function(Sigma, s, d2){
  added_dim <- FALSE
  if (length(dim(Sigma))==2) {Sigma <- add_array_dim(Sigma,3); added_dim=TRUE}
  d <- dim(Sigma)
  d[1] <- d[2] <- d[1]-1
  O <- array(0, dim=d)
  G1 <- driver::create_alr_base(s, d2, inv=FALSE)
  for (i in 1:dim(Sigma)[3]){
    one <- 1:s; two <- (s+1):dim(Sigma)[1]
    O[1:(s-1), 1:(s-1), i] <- t(G1) %*% Sigma[one,one,i] %*% G1
    O[1:(s-1),two-1,i] <-   t(G1) %*% Sigma[one,two,i]
    O[two-1,1:(s-1),i] <- t(O[1:(s-1),two-1,i])
    O[two-1,two-1,i] <- Sigma[two,two,i]
  }
  if (added_dim) return(matrix_maintain_dim(O[,,1]))
  return(O)
}


#' @rdname convert_orthus_covariance
#' @export
oalrvar2alrvar <- function(Sigma, s, d1, d2){
  O <- oalrvar2clrvar(Sigma, s, d1)
  oclrvar2alrvar(O, s+1, d2)
}

#' @rdname convert_orthus_covariance
#' @export
oalrvar2ilrvar <- function(Sigma, s, d1, V2){
  O <- oalrvar2clrvar(Sigma, s, d1)
  oclrvar2ilrvar(O, s+1, V2)
}

#' @rdname convert_orthus_covariance
#' @export
oilrvar2alrvar <- function(Sigma, s, V1, d2){
  O <- oilrvar2clrvar(Sigma, s, V1)
  oclrvar2alrvar(O, s+1, d2)
}








