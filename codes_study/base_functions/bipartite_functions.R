### H2 #####
H2fun2=function (web, H2_integer = TRUE) 
{
  if (H2_integer & any((web%%1) != 0)) 
    stop("web does not contain integers! maybe you should set H2_integer to FALSE")
  tot <- sum(web)
  rs <- rowSums(web)
  cs <- colSums(web)
  H2uncorr = -sum(web/tot * log(web/tot), na.rm = TRUE)
  exexpec <- outer(rs, cs/tot)
  if (!H2_integer) {
    newweb <- exexpec
    H2_max <- -sum(newweb/tot * log(newweb/tot), na.rm = TRUE)
  }
  else {
    expec <- matrix(0, nrow(web), ncol(web))
    difexp <- exexpec - expec
    newweb <- floor(exexpec)
    webfull <- matrix("no", nrow(web), ncol(web))
    while (sum(newweb) < tot) {
      webfull[which(rowSums(newweb) == rs), ] <- "yo"
      webfull[, which(colSums(newweb) == cs)] <- "yo"
      OK <- webfull == "no"
      smallestpossible <- newweb == min(newweb[OK])
      greatestdif <- max(difexp[smallestpossible & OK])
      bestone <- which(OK & smallestpossible & difexp == 
                         greatestdif)
      if (length(bestone) > 1) 
        bestone <- sample(bestone, 1)
      newweb[bestone] <- newweb[bestone] + 1
      difexp <- exexpec - newweb
    }
    H2_max <- -sum(newweb/tot * log(newweb/tot), na.rm = TRUE)
    if (max(exexpec) > 0.3679 * tot) {
      for (tries in 1:500) {
        newmx <- newweb
        difexp <- exexpec - newmx
        greatestdif <- difexp == min(difexp)
        if (length(which(greatestdif)) > 1) {
          largestvalue = newmx == max(newmx[greatestdif])
          first <- greatestdif & largestvalue
        }
        else {
          first = greatestdif
        }
        newmx[first][1] <- newmx[first][1] - 1
        throw = which(rowSums(first) > 0)[1]
        thcol = which(colSums(first) > 0)[1]
        mr = max(difexp[throw, ])
        mc = max(difexp[, thcol])
        if (mr >= mc) {
          scnd = which(difexp[throw, ] == mr)[1]
          newmx[throw, scnd] = newmx[throw, scnd] + 1
          thrd = which(difexp[, scnd] == min(difexp[, 
                                                    scnd]))[1]
          newmx[thrd, scnd] = newmx[thrd, scnd] - 1
          newmx[thrd, thcol] = newmx[thrd, thcol] + 1
        }
        else {
          scnd = which(difexp[, thcol] == mc)[1]
          newmx[scnd, thcol] = newmx[scnd, thcol] + 1
          thrd = which(difexp[scnd, ] == min(difexp[scnd, 
          ]))[1]
          newmx[scnd, thrd] = newmx[scnd, thrd] - 1
          newmx[throw, thrd] = newmx[throw, thrd] + 1
        }
      }
      newweb <- newmx
    }
  }
  H2_max.improved <- -sum(newweb/tot * log(newweb/tot), na.rm = TRUE)
  H2_max <- ifelse(H2_integer && H2_max >= H2_max.improved, 
                   H2_max, H2_max.improved)
  newweb <- matrix(0, length(rs), length(cs))
  rsrest = rs
  csrest = cs
  while (round(sum(rsrest), 10) != 0) {
    newweb[which(rsrest == max(rsrest))[1], which(csrest == 
                                                    max(csrest))[1]] = min(c(max(rsrest), max(csrest)))
    rsrest = rs - rowSums(newweb)
    csrest = cs - colSums(newweb)
  }
  Pnew <- newweb/sum(newweb)
  H2_min <- -sum(Pnew * log(Pnew), na.rm = TRUE)
  if (H2uncorr < H2_min) 
    H2_min <- H2uncorr
  if (H2_max < H2uncorr) 
    H2_max <- H2uncorr
  H_2prime <- (H2_max - H2uncorr)/(H2_max - H2_min)
  c(H2 = H_2prime, H2min = round(H2_min, 3), H2max = round(H2_max, 
                                                           3), H2uncorr = round(H2uncorr, 3))
}
### nest.smdm ####
nest.smdm2=function (x, constraints = NULL, weighted = FALSE, decreasing = "fill", 
                     sort = TRUE) 
{
  if (!is.null(constraints) & length(unique(constraints)) == 
      1) {
    warning("Only one module. Nestedness calculated only for the entire matrix")
    constraints = NULL
  }
  if (is.element(NA, constraints) | is.element(NaN, constraints)) {
    warning("NA or NaN in constraints. Nestedness calculated only for the entire matrix")
    constraints = NULL
  }
  if (!is.null(constraints) & length(constraints) != nrow(x) + 
      ncol(x)) {
    stop("constraints vector is not of the same length that network vertices")
  }
  if (weighted == FALSE & any(x != 0 & x != 1)) {
    x[x > 0] = 1
    warning("binary metric applied")
  }
  if (decreasing != "fill" & decreasing != "abund") {
    stop("decreasing should be fill or abund")
  }
  if (!is.null(constraints)) {
    constraints = as.character(constraints)
  }
  if (is.null(rownames(x))) {
    xrnames <- paste("R", 1:nrow(x), "")
    rownames(x) <- xrnames
  }
  if (is.null(colnames(x))) {
    xcnames <- paste("C", 1:ncol(x), "")
    colnames(x) <- xcnames
  }
  unweightednodf = function(x, constraints) {
    if (sort == TRUE) {
      tab0 = x[sort(rowSums(x), index = TRUE, decreasing = TRUE)$ix, 
               sort(colSums(x), index = TRUE, decreasing = TRUE)$ix]
    }
    else {
      tab0 = x
    }
    MTrow = rowSums(tab0)
    Nrow = matrix(rep(NA, times = nrow(tab0)^2), nrow(tab0), 
                  nrow(tab0))
    dimnames(Nrow) = list(rownames(tab0), rownames(tab0))
    for (jrow in 2:nrow(tab0)) {
      for (irow in 1:(jrow - 1)) {
        if (MTrow[jrow] >= MTrow[irow]) {
          Nrow[jrow, irow] = 0
        }
        else {
          S = 0
          for (i in 1:ncol(tab0)) {
            if (tab0[jrow, i] == 1 & tab0[jrow, i] == 
                tab0[irow, i]) {
              S = S + 1
            }
          }
          Nrow[jrow, irow] = S * 100/MTrow[jrow]
        }
      }
    }
    Nrow = Nrow[rownames(x), rownames(x)]
    NODFrow = mean(Nrow, na.rm = TRUE)
    MTcol = colSums(tab0)
    Ncol = matrix(rep(NA, times = ncol(tab0)^2), ncol(tab0), 
                  ncol(tab0))
    dimnames(Ncol) = list(colnames(tab0), colnames(tab0))
    for (jcol in 2:ncol(tab0)) {
      for (icol in 1:(jcol - 1)) {
        if (MTcol[jcol] >= MTcol[icol]) {
          Ncol[jcol, icol] = 0
        }
        else {
          S = 0
          for (i in 1:nrow(tab0)) {
            if (tab0[i, jcol] == 1 & tab0[i, jcol] == 
                tab0[i, icol]) {
              S = S + 1
            }
          }
          Ncol[jcol, icol] = S * 100/MTcol[jcol]
        }
      }
    }
    Ncol = Ncol[colnames(x), colnames(x)]
    NODFcol = mean(Ncol, na.rm = TRUE)
    NODFmatrix = mean(c(Ncol, Nrow), na.rm = TRUE)
    if (!is.null(constraints)) {
      rowcons = cbind(rownames(x), constraints[1:nrow(x)])
      tabrcons = table(rowcons[, 1], rowcons[, 2])
      distrcons = dist(tabrcons, method = "binary")
      distrcons = as.matrix(distrcons)
      distrcons = distrcons[rownames(x), rownames(x)]
      rm(rowcons, tabrcons)
      SM_Nrow = 0
      SM_nrow = 0
      DM_Nrow = 0
      DM_nrow = 0
      for (i in 1:nrow(x)) {
        for (j in 1:nrow(x)) {
          if (!is.na(Nrow[i, j])) {
            if (distrcons[i, j] == 0) {
              SM_Nrow = SM_Nrow + Nrow[i, j]
              SM_nrow = SM_nrow + 1
            }
            else {
              DM_Nrow = DM_Nrow + Nrow[i, j]
              DM_nrow = DM_nrow + 1
            }
          }
        }
      }
      NODF_SM_row = SM_Nrow/SM_nrow
      NODF_DM_row = DM_Nrow/DM_nrow
      colcons = cbind(colnames(x), constraints[(nrow(x) + 
                                                  1):length(constraints)])
      tabccons = table(colcons[, 1], colcons[, 2])
      distccons = dist(tabccons, method = "binary")
      distccons = as.matrix(distccons)
      distccons = distccons[colnames(x), colnames(x)]
      rm(colcons, tabccons)
      SM_Ncol = 0
      SM_ncol = 0
      DM_Ncol = 0
      DM_ncol = 0
      for (i in 1:ncol(x)) {
        for (j in 1:ncol(x)) {
          if (!is.na(Ncol[i, j])) {
            if (distccons[i, j] == 0) {
              SM_Ncol = SM_Ncol + Ncol[i, j]
              SM_ncol = SM_ncol + 1
            }
            else {
              DM_Ncol = DM_Ncol + Ncol[i, j]
              DM_ncol = DM_ncol + 1
            }
          }
        }
      }
      NODF_SM_col = SM_Ncol/SM_ncol
      NODF_DM_col = DM_Ncol/DM_ncol
      NODF_SM_matrix = (SM_Nrow + SM_Ncol)/(SM_nrow + SM_ncol)
      NODF_DM_matrix = (DM_Nrow + DM_Ncol)/(DM_nrow + DM_ncol)
      return(list(NODFrow = NODFrow, NODFcol = NODFcol, 
                  NODFmatrix = NODFmatrix, NODF_SM_row = NODF_SM_row, 
                  NODF_DM_row = NODF_DM_row, NODF_SM_col = NODF_SM_col, 
                  NODF_DM_col = NODF_DM_col, NODF_SM_matrix = NODF_SM_matrix, 
                  NODF_DM_matrix = NODF_DM_matrix))
    }
    else {
      return(list(NODFrow = NODFrow, NODFcol = NODFcol, 
                  NODFmatrix = NODFmatrix))
    }
  }
  weightednodf = function(x, constraints) {
    if (sort == TRUE) {
      tab0 = x[sort(rowSums(x != 0), index = TRUE, decreasing = TRUE)$ix, 
               sort(colSums(x != 0), index = TRUE, decreasing = TRUE)$ix]
    }
    else {
      tab0 = x
    }
    MTrow = rowSums(tab0)
    Frow = rowSums(tab0 != 0)
    Nrow = matrix(rep(NA, times = nrow(tab0)^2), nrow(tab0), 
                  nrow(tab0))
    dimnames(Nrow) = list(rownames(tab0), rownames(tab0))
    for (jrow in 2:nrow(tab0)) {
      for (irow in 1:(jrow - 1)) {
        if (Frow[jrow] >= Frow[irow]) {
          Nrow[jrow, irow] = 0
        }
        else {
          S = 0
          for (i in 1:ncol(tab0)) {
            if (tab0[jrow, i] != 0 & tab0[jrow, i] < 
                tab0[irow, i]) {
              S = S + 1
            }
          }
          Nrow[jrow, irow] = S * 100/Frow[jrow]
        }
      }
    }
    Nrow = Nrow[rownames(x), rownames(x)]
    NODFrow = mean(Nrow, na.rm = TRUE)
    MTcol = colSums(tab0)
    Fcol = colSums(tab0 != 0)
    Ncol = matrix(rep(NA, times = ncol(tab0)^2), ncol(tab0), 
                  ncol(tab0))
    dimnames(Ncol) = list(colnames(tab0), colnames(tab0))
    for (jcol in 2:ncol(tab0)) {
      for (icol in 1:(jcol - 1)) {
        if (Fcol[jcol] >= Fcol[icol]) {
          Ncol[jcol, icol] = 0
        }
        else {
          S = 0
          for (i in 1:nrow(tab0)) {
            if (tab0[i, jcol] != 0 & tab0[i, jcol] < 
                tab0[i, icol]) {
              S = S + 1
            }
          }
          Ncol[jcol, icol] = S * 100/Fcol[jcol]
        }
      }
    }
    Ncol = Ncol[colnames(x), colnames(x)]
    NODFcol = mean(Ncol, na.rm = TRUE)
    NODFmatrix = mean(c(Ncol, Nrow), na.rm = TRUE)
    if (!is.null(constraints)) {
      rowcons = cbind(rownames(x), constraints[1:nrow(x)])
      tabrcons = table(rowcons[, 1], rowcons[, 2])
      distrcons = dist(tabrcons, method = "binary")
      distrcons = as.matrix(distrcons)
      distrcons = distrcons[rownames(x), rownames(x)]
      rm(rowcons, tabrcons)
      SM_Nrow = 0
      SM_nrow = 0
      DM_Nrow = 0
      DM_nrow = 0
      for (i in 1:nrow(x)) {
        for (j in 1:nrow(x)) {
          if (!is.na(Nrow[i, j])) {
            if (distrcons[i, j] == 0) {
              SM_Nrow = SM_Nrow + Nrow[i, j]
              SM_nrow = SM_nrow + 1
            }
            else {
              DM_Nrow = DM_Nrow + Nrow[i, j]
              DM_nrow = DM_nrow + 1
            }
          }
        }
      }
      NODF_SM_row = SM_Nrow/SM_nrow
      NODF_DM_row = DM_Nrow/DM_nrow
      colcons = cbind(colnames(x), constraints[(nrow(x) + 
                                                  1):length(constraints)])
      tabccons = table(colcons[, 1], colcons[, 2])
      distccons = dist(tabccons, method = "binary")
      distccons = as.matrix(distccons)
      distccons = distccons[colnames(x), colnames(x)]
      rm(colcons, tabccons)
      SM_Ncol = 0
      SM_ncol = 0
      DM_Ncol = 0
      DM_ncol = 0
      for (i in 1:ncol(x)) {
        for (j in 1:ncol(x)) {
          if (!is.na(Ncol[i, j])) {
            if (distccons[i, j] == 0) {
              SM_Ncol = SM_Ncol + Ncol[i, j]
              SM_ncol = SM_ncol + 1
            }
            else {
              DM_Ncol = DM_Ncol + Ncol[i, j]
              DM_ncol = DM_ncol + 1
            }
          }
        }
      }
      NODF_SM_col = SM_Ncol/SM_ncol
      NODF_DM_col = DM_Ncol/DM_ncol
      NODF_SM_matrix = (SM_Nrow + SM_Ncol)/(SM_nrow + SM_ncol)
      NODF_DM_matrix = (DM_Nrow + DM_Ncol)/(DM_nrow + DM_ncol)
      return(list(WNODFrow = NODFrow, WNODFcol = NODFcol, 
                  WNODFmatrix = NODFmatrix, WNODF_SM_row = NODF_SM_row, 
                  WNODF_DM_row = NODF_DM_row, WNODF_SM_col = NODF_SM_col, 
                  WNODF_DM_col = NODF_DM_col, WNODF_SM_matrix = NODF_SM_matrix, 
                  WNODF_DM_matrix = NODF_DM_matrix))
    }
    else {
      return(list(WNODFrow = NODFrow, WNODFcol = NODFcol, 
                  WNODFmatrix = NODFmatrix))
    }
  }
  weightednoda = function(x, constraints) {
    if (sort == TRUE) {
      tab0 = x[sort(rowSums(x), index = TRUE, decreasing = TRUE)$ix, 
               sort(colSums(x), index = TRUE, decreasing = TRUE)$ix]
    }
    else {
      tab0 <- x
    }
    MTrow = rowSums(tab0)
    Frow = rowSums(tab0 != 0)
    Nrow = matrix(rep(NA, times = nrow(tab0)^2), nrow(tab0), 
                  nrow(tab0))
    dimnames(Nrow) = list(rownames(tab0), rownames(tab0))
    for (jrow in 2:nrow(tab0)) {
      for (irow in 1:(jrow - 1)) {
        if (MTrow[jrow] >= MTrow[irow]) {
          Nrow[jrow, irow] = 0
        }
        else {
          S = 0
          for (i in 1:ncol(tab0)) {
            if (tab0[jrow, i] != 0 & tab0[jrow, i] < 
                tab0[irow, i]) {
              S = S + 1
            }
          }
          Nrow[jrow, irow] = S * 100/Frow[jrow]
        }
      }
    }
    Nrow = Nrow[rownames(x), rownames(x)]
    NODArow = mean(Nrow, na.rm = TRUE)
    MTcol = colSums(tab0)
    Fcol = colSums(tab0 != 0)
    Ncol = matrix(rep(NA, times = ncol(tab0)^2), ncol(tab0), 
                  ncol(tab0))
    dimnames(Ncol) = list(colnames(tab0), colnames(tab0))
    for (jcol in 2:ncol(tab0)) {
      for (icol in 1:(jcol - 1)) {
        if (MTcol[jcol] >= MTcol[icol]) {
          Ncol[jcol, icol] = 0
        }
        else {
          S = 0
          for (i in 1:nrow(tab0)) {
            if (tab0[i, jcol] != 0 & tab0[i, jcol] < 
                tab0[i, icol]) {
              S = S + 1
            }
          }
          Ncol[jcol, icol] = S * 100/Fcol[jcol]
        }
      }
    }
    Ncol = Ncol[colnames(x), colnames(x)]
    NODAcol = mean(Ncol, na.rm = TRUE)
    NODAmatrix = mean(c(Ncol, Nrow), na.rm = TRUE)
    if (!is.null(constraints)) {
      rowcons = cbind(rownames(x), constraints[1:nrow(x)])
      tabrcons = table(rowcons[, 1], rowcons[, 2])
      distrcons = dist(tabrcons, method = "binary")
      distrcons = as.matrix(distrcons)
      distrcons = distrcons[rownames(x), rownames(x)]
      rm(rowcons, tabrcons)
      SM_Nrow = 0
      SM_nrow = 0
      DM_Nrow = 0
      DM_nrow = 0
      for (i in 1:nrow(x)) {
        for (j in 1:nrow(x)) {
          if (!is.na(Nrow[i, j])) {
            if (distrcons[i, j] == 0) {
              SM_Nrow = SM_Nrow + Nrow[i, j]
              SM_nrow = SM_nrow + 1
            }
            else {
              DM_Nrow = DM_Nrow + Nrow[i, j]
              DM_nrow = DM_nrow + 1
            }
          }
        }
      }
      NODA_SM_row = SM_Nrow/SM_nrow
      NODA_DM_row = DM_Nrow/DM_nrow
      colcons = cbind(colnames(x), constraints[(nrow(x) + 
                                                  1):length(constraints)])
      tabccons = table(colcons[, 1], colcons[, 2])
      distccons = dist(tabccons, method = "binary")
      distccons = as.matrix(distccons)
      distccons = distccons[colnames(x), colnames(x)]
      rm(colcons, tabccons)
      SM_Ncol = 0
      SM_ncol = 0
      DM_Ncol = 0
      DM_ncol = 0
      for (i in 1:ncol(x)) {
        for (j in 1:ncol(x)) {
          if (!is.na(Ncol[i, j])) {
            if (distccons[i, j] == 0) {
              SM_Ncol = SM_Ncol + Ncol[i, j]
              SM_ncol = SM_ncol + 1
            }
            else {
              DM_Ncol = DM_Ncol + Ncol[i, j]
              DM_ncol = DM_ncol + 1
            }
          }
        }
      }
      NODA_SM_col = SM_Ncol/SM_ncol
      NODA_DM_col = DM_Ncol/DM_ncol
      NODA_SM_matrix = (SM_Nrow + SM_Ncol)/(SM_nrow + SM_ncol)
      NODA_DM_matrix = (DM_Nrow + DM_Ncol)/(DM_nrow + DM_ncol)
      return(list(WNODArow = NODArow, WNODAcol = NODAcol, 
                  WNODAmatrix = NODAmatrix, WNODA_SM_row = NODA_SM_row, 
                  WNODA_DM_row = NODA_DM_row, WNODA_SM_col = NODA_SM_col, 
                  WNODA_DM_col = NODA_DM_col, WNODA_SM_matrix = NODA_SM_matrix, 
                  WNODA_DM_matrix = NODA_DM_matrix))
    }
    else {
      return(list(WNODArow = NODArow, WNODAcol = NODAcol, 
                  WNODAmatrix = NODAmatrix))
    }
  }
  if (decreasing == "abund") {
    return(weightednoda(x, constraints))
  }
  if (decreasing == "fill") {
    if (weighted == F) {
      return(unweightednodf(x, constraints))
    }
    if (weighted == TRUE) {
      return(weightednodf(x, constraints))
    }
  }
}
### restricted null ####
PosteriorProb2=function (web, R.partitions = NULL, C.partitions = NULL, Prior.Pij = "degreeprob", 
                         conditional.level = "modules") 
{
  M <- web
  if (conditional.level == "matrix") {
    R.partitions <- rep(1, nrow(M))
    C.partitions <- rep(1, ncol(M))
  }
  if (conditional.level == "modules" & (is.null(R.partitions) | 
                                        is.null(C.partitions))) 
    stop("When using conditional.level='modules', R- and C-partitions must be given, based on modularity. See example.")
  if (!is.matrix(M)) {
    stop("M is not a matrix.")
  }
  if (0 %in% rowSums(M) | 0 %in% colSums(M)) {
    stop("M is degenerated. There are rows and/or columns without interactions in the matrix. Remove them before proceeding.")
  }
  if (!is.numeric(R.partitions) | !is.numeric(C.partitions)) {
    stop("Partitions are not numeric.")
  }
  if (length(R.partitions) != nrow(M) | length(C.partitions) != 
      ncol(M)) {
    stop("Partitions and matrix dimensions have different sizes.")
  }
  if (!(conditional.level %in% c("matrix", "modules", "areas"))) {
    stop("conditional.level should be 'matrix','modules' or 'areas'.")
  }
  if (!(Prior.Pij %in% c("degreeprob", "equiprobable", "degreeprob.byarea"))) {
    stop("Pij.probs should be 'equiprobable' or 'degreeprob' or 'degreeprob.byarea.")
  }
  nr <- dim(M)[1]
  nc <- dim(M)[2]
  Matrix.mod <- array(0, dim = c(nr, nc, 3))
  for (rr in 1:nr) {
    for (cc in 1:nc) {
      Matrix.mod[rr, cc, 1] <- ifelse(R.partitions[rr] == 
                                        C.partitions[cc], 1, 0)
      Matrix.mod[rr, cc, 2] <- R.partitions[rr]
      Matrix.mod[rr, cc, 3] <- C.partitions[cc]
    }
  }
  if (Prior.Pij == "equiprobable") {
    Pi <- rep(1/nr, times = nr)
    Pj <- rep(1/nc, times = nc)
    Prior.Pij.species <- tcrossprod(Pi, Pj)
  }
  if (Prior.Pij == "degreeprob") {
    Pi <- rowSums(M)/sum(rowSums(M))
    Pj <- colSums(M)/sum(colSums(M))
    Prior.Pij.species <- tcrossprod(Pi, Pj)
  }
  if (Prior.Pij == "degreeprob.byarea") {
    Prior.Pij.species <- M
    RMod <- sort(unique(R.partitions))
    CMod <- sort(unique(C.partitions))
    for (rr in RMod) {
      for (cc in CMod) {
        M.rr.cc <- matrix(M[R.partitions == rr, C.partitions == 
                              cc], sum(1 * (R.partitions == rr)), sum(1 * 
                                                                        (C.partitions == cc)))
        Pi.rr.cc <- rowSums(M.rr.cc)/sum(rowSums(M.rr.cc))
        Pj.rr.cc <- colSums(M.rr.cc)/sum(colSums(M.rr.cc))
        Prior.Pij.species[R.partitions == rr, C.partitions == 
                            cc] <- tcrossprod(Pi.rr.cc, Pj.rr.cc)
      }
    }
  }
  if (conditional.level == "matrix") {
    Post.Pij <- Prior.Pij.species
  }
  else {
    Prior.Pij.area <- matrix(NA, nr, nc)
    Cond.Pij.area <- matrix(NA, nr, nc)
    if (conditional.level == "modules") {
      WMod.prior <- sum(Prior.Pij.species[Matrix.mod[, 
                                                     , 1] == 1])
      OMod.prior <- sum(Prior.Pij.species[Matrix.mod[, 
                                                     , 1] == 0])
      Prior.Pij.area[Matrix.mod[, , 1] == 1] <- WMod.prior
      Prior.Pij.area[Matrix.mod[, , 1] == 0] <- OMod.prior
      WMod.cond <- sum(M[Matrix.mod[, , 1] == 1])/sum(M)
      OMod.cond <- sum(M[Matrix.mod[, , 1] == 0])/sum(M)
      Cond.Pij.area[Matrix.mod[, , 1] == 1] <- WMod.cond
      Cond.Pij.area[Matrix.mod[, , 1] == 0] <- OMod.cond
    }
    if (conditional.level == "areas") {
      RMod <- sort(unique(R.partitions))
      CMod <- sort(unique(C.partitions))
      for (rr in RMod) {
        for (cc in CMod) {
          WArea.prior <- sum(Prior.Pij.species[Matrix.mod[, 
                                                          , 2] == rr & Matrix.mod[, , 3] == cc])
          Prior.Pij.area[Matrix.mod[, , 2] == rr & Matrix.mod[, 
                                                              , 3] == cc] <- WArea.prior
          WArea.cond <- sum(M[Matrix.mod[, , 2] == rr & 
                                Matrix.mod[, , 3] == cc])/sum(M)
          Cond.Pij.area[Matrix.mod[, , 2] == rr & Matrix.mod[, 
                                                             , 3] == cc] <- WArea.cond
        }
      }
    }
    Post.Pij <- Prior.Pij.species * (Cond.Pij.area/Prior.Pij.area)
  }
  return(Post.Pij = Post.Pij)
}
restrictednull2= function (web, Prior.Pij = "degreeprob", conditional.level = "modules", 
                           N = 10, print.null = FALSE, allow.degeneration = FALSE, return.nonrm.species = FALSE, 
                           connectance = TRUE, byarea = FALSE, R.partitions = NULL, 
                           C.partitions = NULL, NCores=1) 
{
  M <- web
  Pij.Prob <- PosteriorProb2(web = M, R.partitions = R.partitions, 
                             C.partitions = C.partitions, Prior.Pij = Prior.Pij, conditional.level = conditional.level)
  if (!is.matrix(M)) {
    stop("M is not a matrix.")
  }
  if (0 %in% rowSums(M) | 0 %in% colSums(M)) {
    stop("M is degenerated.")
  }
  if (!is.matrix(Pij.Prob)) {
    stop("Pij is not a matrix.")
  }
  if (T %in% c(Pij.Prob < 0)) {
    stop("Pij must contain only values >= 0.")
  }
  if (nrow(M) != nrow(Pij.Prob) | ncol(M) != ncol(Pij.Prob)) {
    stop("Dimensions of M and Pij.Prob must be identical.")
  }
  if (byarea) {
    if (is.null(C.partitions) | is.null(R.partitions)) {
      stop("Partitions missing")
    }
    if (length(unique(c(length(R.partitions), nrow(M), nrow(Pij.Prob)))) != 
        1) {
      stop("The number of elements of R.partition should be the same as the number of rows of M and Pij.prob.")
    }
    if (length(unique(c(length(C.partitions), ncol(M), ncol(Pij.Prob)))) != 
        1) {
      stop("The number of elements of C.partition should be the same as the number of column of M and Pij.prob.")
    }
    if (!identical(sort(unique(R.partitions)), sort(unique(C.partitions)))) {
      stop("The number and labels of modules in R.partition and C.partition must be the same")
    }
  }
  if (N <= 0 | !is.numeric(N)) {
    stop("Number of nulls, N, must be > 0")
  }
  if (!is.logical(connectance)) {
    stop("Connectance must be logical (T or F)")
  }
  if (!is.logical(allow.degeneration)) {
    stop("allow.degeneration must be logical (T or F)")
  }
  if (!is.logical(return.nonrm.species)) {
    stop("return.nonrm.species must be logical (T or F)")
  }
  if (!is.logical(byarea)) {
    stop("byarea must be logical (T or F)")
  }
  nr <- dim(M)[1]
  nc <- dim(M)[2]
  if (byarea) {
    Matrix.area <- array(0, dim = c(nr, nc, 2))
    for (rr in 1:nr) {
      for (cc in 1:nc) {
        Matrix.area[rr, cc, 1] <- R.partitions[rr]
        Matrix.area[rr, cc, 2] <- C.partitions[cc]
      }
    }
  }
  else {
    Matrix.area <- array(1, dim = c(nr, nc, 2))
    R.partitions <- rep(1, nrow(M))
    C.partitions <- rep(1, ncol(M))
  }
  cl <- makeCluster(NCores)
  registerDoSNOW(cl)
  NullMatrices=foreach(j=1:N,.verbose = T)%dopar%{
    R.part <- sort(unique(as.vector(Matrix.area[, , 1])))
    C.part <- sort(unique(as.vector(Matrix.area[, , 2])))
    finalmat <- matrix(NA, nr, nc)
    for (R.p in R.part) {
      for (C.p in C.part) {
        M.a <- as.matrix(M[R.partitions == R.p, C.partitions == 
                             C.p])
        Pij.a <- Pij.Prob[R.partitions == R.p, C.partitions == 
                            C.p]
        r.a <- dim(M.a)[1]
        c.a <- dim(M.a)[2]
        P.a <- P1.a <- Pij.a
        finalmat.a <- matrix(0, r.a, c.a)
        if (allow.degeneration == FALSE & R.p == C.p) {
          D.int.finalmat.a <- 0
          while (D.int.finalmat.a < sum(dim(M.a))) {
            sel <- sample(1:length(M.a), 1, prob = P.a)
            finalmat.a[sel] <- 1
            P.a[outer(as.numeric(rowSums(finalmat.a) > 
                                   0), as.numeric(colSums(finalmat.a) > 0)) == 
                  1] <- 0
            D.int.finalmat.a <- sum(rowSums(finalmat.a) > 
                                      0) + sum(colSums(finalmat.a) > 0)
          }
        }
        conn.remain <- sum(M.a > 0) - sum(finalmat.a > 
                                            0)
        if (conn.remain > 0) {
          if (connectance == T) {
            if (length(which(finalmat.a == 0)) == 1) {
              add <- which(finalmat.a == 0)
            }
            else {
              add <- sample(which(finalmat.a == 0), conn.remain, 
                            prob = P1.a[finalmat.a == 0], replace = FALSE)
            }
          }
          else {
            add <- sample(1:length(finalmat.a), conn.remain, 
                          prob = P1.a, replace = TRUE)
          }
          for (add1 in add) {
            finalmat.a[add1] <- finalmat.a[add1] + 1
          }
        }
        int.remain <- (sum(M.a) - sum(finalmat.a))
        if (int.remain > 0) {
          if (length(which(finalmat.a > 0)) == 1) {
            add <- rep(which(finalmat.a > 0), int.remain)
          }
          else {
            add <- sample(which(finalmat.a > 0), int.remain, 
                          prob = P1.a[which(finalmat.a > 0)], replace = TRUE)
          }
          finalmat.a[as.numeric(names(table(add)))] <- finalmat.a[as.numeric(names(table(add)))] + 
            (table(add))
        }
        finalmat[R.partitions == R.p, C.partitions == 
                   C.p] <- finalmat.a
      }
    }
    R2keep <- which(rowSums(finalmat) != 0)
    C2keep <- which(colSums(finalmat) != 0)
    finalmat2 <- finalmat[R2keep, C2keep]
    if (return.nonrm.species) {
      output = list(NullMatrix = finalmat2, 
                                R.Kept = R2keep, C.Kept = C2keep)
    }
    else {
      output = finalmat2
    }
    if (print.null) {
      print(nn)
    }
    output
  }
  stopCluster(cl)
  return(NullMatrices)
}
#### propnull ####
propnull2= function (M, N=10, NCores=1){
  S=sum(M)
  PROB=outer(X = rowSums(M),Y = colSums(M),FUN = function(X,Y){(X/S)*(Y/S)})
  # PROB = probabilities based on relative marginal sums
  cl <- makeCluster(NCores)
  registerDoSNOW(cl)
  output=foreach(j=1:N,.verbose = T)%dopar%{
    nullM=matrix(0,nrow=nrow(M),ncol=ncol(M))
    INT=sample(x = 1:length(M), size = S, replace=T, prob = PROB) #sorting interactions
    P=as.numeric(names(table(INT)))
    INTs=as.numeric(table(INT))
    nullM[P]=INTs #including in the null matrix
    nullM=nullM[rowSums(nullM)>0,colSums(nullM)>0] #removing possible empty rows/columns
    nullM
  }
  stopCluster(cl)
  output
}
### DIRT LPA wb + ####
# LPA_wb_plus.R
# Label propagation algorithm for weighted bipartite networks that finds modularity.
# Contains the LPAwb+ and the DIRTLPAwb+ algorithms
# Author :  Stephen Beckett ( https://github.com/sjbeckett/weighted-modularity-LPAwbPLUS )
# MIT License


LPA_wb_plus2 <- function(MATRIX,initialmoduleguess=NA) {
  
  #Make sure the smallest matrix dimension represent the red labels by making them the rows (If matrix is transposed here, will be transposed back at the end)
  flipped = 0
  if(dim(MATRIX)[1] > dim(MATRIX)[2]) {
    MATRIX = t(MATRIX)
    flipped = 1
  }
  
  Matsum = sum(MATRIX)
  col_marginals = colSums(MATRIX)
  row_marginals = rowSums(MATRIX)
  BMatrix = BarbersMatrix(MATRIX)
  
  #initiliase labels
  bluelabels = rep(NA,dim(MATRIX)[2])
  
  if (is.na(initialmoduleguess))
    redlabels = 1:dim(MATRIX)[1]
  else
    redlabels = sample(1:(initialmoduleguess+1),dim(MATRIX)[1],replace=TRUE)
  
  #Run Phase 1: Locally update lables to maximise Qb
  outlist = StageOne_LPAwbdash(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels)
  redlabels = outlist[[1]]
  bluelabels = outlist[[2]]
  Qb_now = outlist[[3]]
  
  #Run Phase 2: Connect divisions from top-down if it improves Qb, then run
  #phase 1 again. Repeat till Qb cannot be improved.
  outlist = StageTwo_LPAwbdash(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels,Qb_now)
  redlabels = outlist[[1]]
  bluelabels = outlist[[2]]
  Qb_now = outlist[[3]]
  
  if(flipped==1) { #If matrix was flipped, swap row and column labels
    holder = redlabels
    redlabels = bluelabels
    bluelabels = holder
  }
  
  return(list(Row_labels=redlabels, Col_labels=bluelabels, modularity=Qb_now))
}

DIRT_LPA_wb_plus2 <- function(MATRIX,mini=4,reps=10) {
  A=LPA_wb_plus2(MATRIX)
  
  mods=length(unique(A[[1]]))
  
  if((mods-mini) > 0) {
    for(aa in mini:mods) {
      for(bb in 1:reps) {
        B=LPA_wb_plus2(MATRIX,aa)
        if(B[[3]]>A[[3]])
          A=B
      }
    }
  }
  
  return(list(Row_labels=A[[1]], Col_labels=A[[2]], modularity=A[[3]]))
}





BarbersMatrix <- function(MATRIX) {
  return(MATRIX - (cbind(rowSums(MATRIX))%*%rbind(colSums(MATRIX)))/sum(MATRIX))
}


WEIGHTEDMODULARITY <- function(BMatrix,Matsum,redlabels,bluelabels) {
  #see equation 8
  holdsum = 0
  
  for (rr in 1:length(redlabels)) {
    for (cc in 1:length(bluelabels)) {
      kroneckerdelta = redlabels[rr] == bluelabels[cc]
      holdsum = holdsum + BMatrix[rr,cc] * kroneckerdelta
    }
  }
  return(holdsum/Matsum)
}


TRACE <- function(MATRIX) { return(sum(diag(MATRIX))) }



WEIGHTEDMODULARITY2 <- function(BMatrix,Matsum,redlabels,bluelabels) {
  #see equation 9
  UNIred = unique(redlabels)
  Lred = length(UNIred)
  UNIblu = unique(bluelabels)
  Lblu = length(UNIblu)
  LABELMAT1 = matrix(0,Lred,length(redlabels))
  LABELMAT2 = matrix(0,length(bluelabels),Lblu)
  
  for(aa in 1:length(redlabels))
    LABELMAT1[which(UNIred == redlabels[aa]),aa] = 1
  
  for(aa in 1:length(bluelabels))
    LABELMAT2[aa,which(UNIblu == bluelabels[aa])] = 1
  
  return(TRACE(LABELMAT1 %*% BMatrix %*% LABELMAT2)/Matsum)
  
}



DIVISION <- function(redlabels,bluelabels) {
  
  divisionsFound <- intersect(redlabels,bluelabels)
  
  return(divisionsFound)
}



StageTwo_LPAwbdash <- function(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels, Qb_now) {
  
  divisionsFound = DIVISION(redlabels,bluelabels)
  NUMdiv = length(divisionsFound)
  IterateFlag = 1
  while(IterateFlag == 1) {
    CombinedDivisionsThisTime = 0
    if(NUMdiv > 1) {
      for(div1check in 1:(NUMdiv-1)) {
        Mod1 = divisionsFound[div1check]
        for(div2check in (div1check+1):NUMdiv) {
          CHECK_RED = redlabels
          CHECK_RED[which(redlabels==Mod1)] = divisionsFound[div2check]
          CHECK_BLUE = bluelabels
          CHECK_BLUE[which(bluelabels==Mod1)] = divisionsFound[div2check]
          
          
          QQ = WEIGHTEDMODULARITY2(BMatrix,Matsum,CHECK_RED,CHECK_BLUE)			
          if(QQ > Qb_now) { #If a division will improve modularity - find the best way to do this
            FoundBetter = 0
            for(aa in 1:NUMdiv) {
              CHECK_RED2 = redlabels
              CHECK_RED2[which(redlabels==divisionsFound[aa])] = Mod1
              CHECK_BLUE2 = bluelabels
              CHECK_BLUE2[which(bluelabels==divisionsFound[aa])] = Mod1
              if(WEIGHTEDMODULARITY2(BMatrix,Matsum,CHECK_RED2,CHECK_BLUE2) > QQ) {
                FoundBetter = 1
              }
              CHECK_RED2 = redlabels
              CHECK_RED2[which(redlabels==divisionsFound[aa])] = divisionsFound[div2check]
              CHECK_BLUE2 = bluelabels
              CHECK_BLUE2[which(bluelabels==divisionsFound[aa])] = divisionsFound[div2check]
              if(WEIGHTEDMODULARITY2(BMatrix,Matsum,CHECK_RED2,CHECK_BLUE2) > QQ) {
                FoundBetter = 1
              }
            }
            if(FoundBetter == 0) { #If no better configuration found - JOIN.
              redlabels = CHECK_RED
              bluelabels = CHECK_BLUE
              CombinedDivisionsThisTime = CombinedDivisionsThisTime + 1
            }
          }
        }
      }
      if(CombinedDivisionsThisTime == 0) {#If no divisions were joined move on
        IterateFlag = 0
      }
    }
    else {
      IterateFlag = 0		
    }
    
    outlist=StageOne_LPAwbdash(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels)  ##
    redlabels = outlist[[1]]
    bluelabels = outlist[[2]]
    Qb_now = outlist[[3]]
    divisionsFound = DIVISION(redlabels,bluelabels)
    NUMdiv = length(divisionsFound)
  }
  
  return(list(redlabels, bluelabels, Qb_now))
}




StageOne_LPAwbdash <- function(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels) {
  #Create storage containers for total marginals attached to each red(row)
  #label and blue(column) label
  
  BLUELABELLENGTH=length(bluelabels)
  REDLABELLENGTH=length(redlabels)
  TotalRedDegrees = rep(NA,max(redlabels))
  TotalBlueDegrees = rep(NA,max(BLUELABELLENGTH,REDLABELLENGTH))
  
  #Fill up these containers according to current labels
  #Red
  for(aa in 1:REDLABELLENGTH) {
    if(is.na(TotalRedDegrees[redlabels[aa]])) {
      TotalRedDegrees[redlabels[aa]] = row_marginals[aa]
    }
    else {
      TotalRedDegrees[redlabels[aa]] = TotalRedDegrees[redlabels[aa]] + row_marginals[aa]
    }
  }
  
  #Blue
  if(sum(is.na(bluelabels)) != BLUELABELLENGTH) { #occurs first time through as blue nodes unlabelled
    for(bb in 1:BLUELABELLENGTH) {
      if(is.na(TotalBlueDegrees[bluelabels[bb]])) {
        TotalBlueDegrees[bluelabels[bb]] = col_marginals[bb]
      }
      else {
        TotalBlueDegrees[bluelabels[bb]] = TotalBlueDegrees[bluelabels[bb]] + col_marginals[bb]
      }
    }
  }
  else {
    TotalBlueDegrees = rep(0,max(BLUELABELLENGTH,REDLABELLENGTH))
  }
  
  
  #locally maximise modularity!
  outlist = LOCALMAXIMISATION(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels,TotalRedDegrees,TotalBlueDegrees)
  redlabels = outlist[[1]]
  bluelabels = outlist[[2]]
  Qb_now = outlist[[3]]
  
  return(list(redlabels, bluelabels, Qb_now))
  
}





LOCALMAXIMISATION <-  function(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels,TotalRedDegrees,TotalBlueDegrees) {
  
  #Find score for current partition
  QbAfter = WEIGHTEDMODULARITY2(BMatrix,Matsum,redlabels,bluelabels)
  
  if(is.na(QbAfter)) { QbAfter = -999 }
  
  IterateFlag = 1
  while(IterateFlag == 1) {
    #Save old information
    QbBefore = QbAfter
    old_redlabels = redlabels
    old_bluelabels = bluelabels
    old_TRD = TotalRedDegrees
    old_TBD = TotalBlueDegrees
    
    #Update Blue Nodes using red node information (see equation 10)
    bluelabelchoices = unique(redlabels)
    
    for(bb in 1:length(bluelabels)) {
      if(is.na(bluelabels[bb]) == FALSE) {
        TotalBlueDegrees[bluelabels[bb]] = TotalBlueDegrees[bluelabels[bb]] - col_marginals[bb] 
      }
      changebluelabeltest = rep(NA,length(bluelabelchoices))
      
      for(ww in 1:length(bluelabelchoices)) {
        changebluelabeltest[ww] = sum( (redlabels == bluelabelchoices[ww]) * MATRIX[,bb])  -  col_marginals[bb]*TotalRedDegrees[bluelabelchoices[ww]]/Matsum 
      }
      
      #assign new label based on maximisation of above condition  
      
      labels = which(changebluelabeltest == max(changebluelabeltest,na.rm =TRUE))
      newlabelindex = labels[sample(1:length(labels),1)]
      bluelabels[bb] = bluelabelchoices[newlabelindex[1]]
      if(bluelabels[bb] > length(TotalBlueDegrees)) {
        TotalBlueDegrees[bluelabels[bb]] = 0
      }
      
      #Update total marginals on new labelling
      TotalBlueDegrees[bluelabels[bb]] = TotalBlueDegrees[bluelabels[bb]] + col_marginals[bb]
    }
    
    #Now update red node labels based on blue node information (see equation 10)
    redlabelchoices = unique(bluelabels)
    
    
    for(aa in 1:length(redlabels)) {
      TotalRedDegrees[redlabels[aa]] = TotalRedDegrees[redlabels[aa]] - row_marginals[aa]
      changeredlabeltest = rep(NA,length(redlabelchoices))
      
      for(ww in 1:length(redlabelchoices)) {
        changeredlabeltest[ww] = sum( (bluelabels == redlabelchoices[ww]) * MATRIX[aa,])  -  row_marginals[aa]*TotalBlueDegrees[redlabelchoices[ww]]/Matsum  
      }
      
      #assign new label based on maximisation of above condition
      labels = which(changeredlabeltest == max(changeredlabeltest,na.rm = TRUE))
      newlabelindex = labels[sample(1:length(labels),1)]
      redlabels[aa] = redlabelchoices[newlabelindex[1]]
      
      if(redlabels[aa] > length(TotalRedDegrees)) {
        TotalRedDegrees[redlabels[aa]] = 0
      }
      TotalRedDegrees[redlabels[aa]] = TotalRedDegrees[redlabels[aa]] + row_marginals[aa]
    }
    
    
    #Find the new modularity score based on node label updates.
    QbAfter = WEIGHTEDMODULARITY(BMatrix,Matsum,redlabels,bluelabels)
    
    #If this modularity is not as good as previous stop iterating and
    #use that previous best information
    
    if(QbAfter <= QbBefore) {
      redlabels = old_redlabels
      bluelabels = old_bluelabels
      TotalRedDegrees = old_TRD
      TotalBlueDegrees = old_TBD
      IterateFlag = 0
    }
    
  }
  
  Qb_now = QbAfter
  
  
  return(list(redlabels, bluelabels, Qb_now))
}






### connectance ####
connectance2=function (M=NULL){
  Mbin=M
  Mbin[M>0]=1
  return(sum(Mbin)/(nrow(Mbin)*ncol(Mbin)))
  }