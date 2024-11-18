hmd.mx2 =  function (country, username, password, label = country, Ontario = FALSE) 
{
  if(Ontario ==F)
  {
  path <- paste("https://www.mortality.org/hmd/", country, "/STATS/", 
                "Mx_1x1.txt", sep = "")
  userpwd <- paste(username, ":", password, sep = "")
  txt <- RCurl::getURL(path, userpwd = userpwd)
  con <- textConnection(txt)
  mx <- try(utils::read.table(con, skip = 2, header = TRUE, 
                              na.strings = "."), TRUE)
  close(con)
  if (class(mx) == "try-error") 
    stop("Connection error at www.mortality.org. Please check username, password and country label.")
  path <- paste("https://www.mortality.org/hmd/", country, "/STATS/", 
                "Exposures_1x1.txt", sep = "")
  userpwd <- paste(username, ":", password, sep = "")
  txt <- RCurl::getURL(path, userpwd = userpwd)
  con <- textConnection(txt)
  pop <- try(utils::read.table(con, skip = 2, header = TRUE, 
                               na.strings = "."), TRUE)
  close(con)
  if (class(pop) == "try-error") 
    stop("Exposures file not found at www.mortality.org")
    }else if( Ontario ==T)
    {
      path <- "https://www.prdh.umontreal.ca/BDLC/data/ont/Mx_1x1.txt"
      txt <- RCurl::getURL(path)
      #tcon <- textConnection(txt) #original line - if ontario == T 'con' is not defined, so it needs to be here as well as line 9
      con <- textConnection(txt)
      mx <- try(utils::read.table(con, skip = 2, header = TRUE, 
                                  na.strings = "."), TRUE)
      close(con)
      if (class(mx) == "try-error") 
        stop("Connection error at www.mortality.org. Please check username, password and country label.") 
      path <- "https://www.prdh.umontreal.ca/BDLC/data/ont/Population.txt"
      txt <- RCurl::getURL(path)
      con <- textConnection(txt)
      pop <- try(utils::read.table(con, skip = 2, header = TRUE, 
                                   na.strings = "."), TRUE)
      close(con)
      if (class(pop) == "try-error") 
        stop("Exposures file not found at www.mortality.org")   
    }
  obj <- list(type = "mortality", label = label, lambda = 0)
  obj$year <- sort(unique(mx[, 1]))
  n <- length(obj$year)
  m <- length(unique(mx[, 2]))
  obj$age <- mx[1:m, 2]
  mnames <- names(mx)[-c(1, 2)]
  n.mort <- length(mnames)
  obj$rate <- obj$pop <- list()
  for (i in 1:n.mort) {
    obj$rate[[i]] <- matrix(mx[, i + 2], nrow = m, ncol = n)
    obj$rate[[i]][obj$rate[[i]] < 0] <- NA
    obj$pop[[i]] <- matrix(pop[, i + 2], nrow = m, ncol = n)
    obj$pop[[i]][obj$pop[[i]] < 0] <- NA
    dimnames(obj$rate[[i]]) <- dimnames(obj$pop[[i]]) <- list(obj$age, 
                                                              obj$year)
  }
  names(obj$pop) = names(obj$rate) <- tolower(mnames)
  obj$age <- as.numeric(as.character(obj$age))
  if (is.na(obj$age[m])) 
    obj$age[m] <- 2 * obj$age[m - 1] - obj$age[m - 2]
  return(structure(obj, class = "demogdata"))
}
