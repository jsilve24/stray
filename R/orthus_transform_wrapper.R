#' Log-Ratio transforms for orthus objects
#' @param x orthus data array (e.g., first s rows are multinomial parameters or log-ratios)
#' @param s first s rows of x are transformed
#' @param V transformation matrix (defines transform) 
#' @param d for ALR, which component (integer position) to take as reference
#' (default is ncol(x)) for alrInv corresponds to column position in untransformed
#' matrix.
#' @param inv for ALR and CLR, transformation matrix is different forward and inverse
#' @param D the number of parts (e.g., number of columns in untransformed data)
#' @name orthus_lr_transforms
NULL

#' @rdname orthus_lr_transforms
#' @export
oglr <- function(x,s, V){
  x.star <- glr_array(x[1:s,,],V, parts=1)
  d.star <- dim(x.star)[1]
  n <-  dim(x)[1] + d.star - s
  y <- array(0, dim=c(n, dim(x)[2], dim(x)[3]))
  y[1:d.star,,] <- x.star
  y[(d.star+1):n,,] <- x[(s+1):dim(x)[1],,]
  return(y)
}

#' @rdname orthus_lr_transforms
#' @export
oglrInv <- function(x, s, V){
  x.star <- glrInv_array(x[1:s,,],V, coords=1)
  d.star <- dim(x.star)[1]
  n <-  dim(x)[1] + d.star - s
  y <- array(0, dim=c(n, dim(x)[2], dim(x)[3]))
  y[1:d.star,,] <- x.star
  y[(d.star+1):n,,] <- x[(s+1):dim(x)[1],,]
  return(y)
}

#' @rdname orthus_lr_transforms
#' @export
oalr <- function(x, s, d=NULL){
  if (is.null(d)) d <- s
  B <- create_alr_base(s, d, inv=FALSE)
  oglr(x, s, B)
}

#' @rdname orthus_lr_transforms
#' @export
oalrInv <- function(y, s, d=NULL){
  if (is.null(d)) d <- s+1
  B <- create_alr_base(s+1, d, inv=TRUE)
  oglrInv(y, s, B)
}

#' @rdname orthus_lr_transforms
#' @export
oilr <- function(x, s, V=NULL){
  if (is.null(V)) V <- create_default_ilr_base(s)
  oglr(x, s, V)
}

#' @rdname orthus_lr_transforms
#' @export
oilrInv <- function(y, s, V=NULL){
  if (is.null(V)) V <- create_default_ilr_base(s+1)
  oglrInv(y, s, V)
}

#' @rdname orthus_lr_transforms
#' @export
oclr <- function(x, s){
  oglr(x, s, create_clr_base(s))
}

#' @rdname orthus_lr_transforms
#' @export
oclrInv <- function(x, s){
  x.star <- clrInv_array(x[1:s,,], coords=1)
  d.star <- dim(x.star)[1]
  n <-  dim(x)[1] + d.star - s
  y <- array(0, dim=c(n, dim(x)[2], dim(x)[3]))
  y[1:d.star,,] <- x.star
  y[(d.star+1):n,,] <- x[(s+1):dim(x)[1],,]
  return(y)
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
  for (i in 1:dim(Sigma)[3]){
    one <- 1:s; two <- (s+1):dim(Sigma)[1]
    Sigma[1:s, 1:s, i] <-  t(V2) %*% V1 %*% Sigma[one,one,i] %*% t(V1) %*% V2
    Sigma[one,two,i] <- t(V2) %*% V1 %*% Sigma[one,two,i]
    Sigma[two,one,i] <- t(Sigma[one,two,i])
  }
  return(Sigma)
}

#' @rdname convert_orthus_covariance
#' @export
oilrvar2clrvar <- function(Sigma, s, V){
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
  return(O)
}

#' @rdname convert_orthus_covariance
#' @export
oclrvar2ilrvar <- function(Sigma, s, V){
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
  return(O)
}

#' @rdname convert_orthus_covariance
#' @export
oalrvar2clrvar <- function(Sigma, s, d1){
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
  return(O)
}

#' @rdname convert_orthus_covariance
#' @export
oclrvar2alrvar <- function(Sigma, s, d2){
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
  return(O)
}


#' @rdname convert_orthus_covariance
#' @export
oalrvar2alrvar <- function(Sigma, s, d1, d2){
  O <- oalrvar2clrvar(Sigma, s, d1)
  oclrvar2alrvar(O, d2)
}

#' @rdname convert_orthus_covariance
#' @export
oalrvar2ilrvar <- function(Sigma, s, d1, V2){
  O <- oalrvar2clrvar(Sigma, s, d1)
  oclrvar2ilrvar(O, V2)
}

#' @rdname convert_orthus_covariance
#' @export
oilrvar2alrvar <- function(Sigma, s, V1, d2){
  O <- oilrvar2clrvar(Sigma, s, V1)
  oclrvar2alrvar(O, d2)
}








