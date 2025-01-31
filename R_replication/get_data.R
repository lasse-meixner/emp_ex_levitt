### attempt to exactly replicate data pre-processing from BCH-2014

source("auxiliary.R")

# We'll need a norminv-like function (not actually used if no LASSO):
norminv <- function(p) { qnorm(p) }

##############################################################################
# 1) READ AND TRIM THE RAW DATA
##############################################################################

data_raw <- read.table(
  "../BCH-2014-replication/ReplicationFilesMS17285/levitt_ex.dat",
  sep="\t", skip=1, header=FALSE
)

# For consistency with Matlab indexing:
#   1: state
#   2: year
#   3: population
#   4: y_viol
#   5: y_prop
#   6: y_murd
#   7: x_murd
#   8: x_viol
#   9: x_prop
#   10..end: other controls (xx)
ncol_data <- ncol(data_raw)
colnames(data_raw) <- c(
  "state","year","pop","y_viol","y_prop","y_murd",
  "x_murd","x_viol","x_prop",
  paste0("X", 10:ncol_data)
)

# Remove DC (9), Alaska(2), Hawaii(12)
tmp <- subset(data_raw, state != 9)
tmp <- subset(tmp, state != 2)
tmp <- subset(tmp, state != 12)

# Keep only years in [85..97]
tmp <- subset(tmp, year >= 85 & year <= 97)

# Adjust year => year - 84
tmp$year <- tmp$year - 84

# recode state
tmp$state <- recode(tmp$state)

# Extract main pieces
y_viol <- tmp$y_viol
y_prop <- tmp$y_prop
y_murd <- tmp$y_murd

x_murd <- tmp$x_murd
x_viol <- tmp$x_viol
x_prop <- tmp$x_prop

xx     <- as.matrix(tmp[, 10:ncol_data])  # columns 10..end

state  <- tmp$state
year   <- tmp$year
pop    <- tmp$pop  # not heavily used, but in original code for weighting
N      <- max(state)
T      <- max(year)

##############################################################################
# 2) BUILD "district" EXACTLY AS IN CODE
##############################################################################
district <- rep(0, length(state))
district_vec <- rep(0, N)

district_vec[1]  <- 8
district_vec[2]  <- 2
district_vec[3]  <- 4
district_vec[4]  <- 1
district_vec[5]  <- 2
district_vec[6]  <- 7
district_vec[7]  <- 9
district_vec[8]  <- 9
district_vec[9]  <- 9
district_vec[10] <- 2
district_vec[11] <- 5
district_vec[12] <- 5
district_vec[13] <- 3
district_vec[14] <- 3
district_vec[15] <- 8
district_vec[16] <- 4
district_vec[17] <- 7
district_vec[18] <- 9
district_vec[19] <- 7
district_vec[20] <- 5
district_vec[21] <- 3
district_vec[22] <- 8
district_vec[23] <- 3
district_vec[24] <- 2
district_vec[25] <- 3
district_vec[26] <- 2
district_vec[27] <- 7
district_vec[28] <- 6
district_vec[29] <- 2
district_vec[30] <- 6
district_vec[31] <- 9
district_vec[32] <- 3
district_vec[33] <- 5
district_vec[34] <- 4
district_vec[35] <- 1
district_vec[36] <- 6
district_vec[37] <- 7
district_vec[38] <- 9
district_vec[39] <- 3
district_vec[40] <- 8
district_vec[41] <- 4
district_vec[42] <- 2
district_vec[43] <- 7
district_vec[44] <- 9
district_vec[45] <- 1
district_vec[46] <- 9
district_vec[47] <- 5
district_vec[48] <- 2

district <- district_vec[state]

##############################################################################
# 3) FIRST-DIFFERENCES: Dyviol, Dxviol, etc.
##############################################################################
# Replicate Matlab reshape(., T, N) => diff => flatten in column-major
matDiff <- function(x, T, N) {
  Xmat <- matrix(x, nrow=T, ncol=N, byrow=FALSE)
  Xdiff <- apply(Xmat, 2, diff)   # diff along each column
  as.vector(Xdiff)                # flatten columnwise
}

Dyviol <- matDiff(y_viol, T, N)
Dxviol <- matDiff(x_viol, T, N)
Dyprop <- matDiff(y_prop, T, N)
Dxprop <- matDiff(x_prop, T, N)
Dymurd <- matDiff(y_murd, T, N)
Dxmurd <- matDiff(x_murd, T, N)

# Also differencing xx columns
Dxx <- matrix(0, nrow=N*(T-1), ncol=ncol(xx))
for (ii in seq_len(ncol(xx))) {
  Dxx[, ii] <- matDiff(xx[, ii], T, N)
}

# We keep only rows with year>1
valid_idx <- which(year > 1)
Dyear     <- year[valid_idx]
Dstate    <- state[valid_idx]
Ddistrict <- district[valid_idx]

# recode Dyear => Dtime
Dtime <- recode(Dyear)
DT    <- max(Dtime)

# also define xxL = xx(year < T, :) 
xxL <- xx[ year < T, , drop=FALSE]

##############################################################################
# 3a) PARTIAL OUT TIME DUMMIES (the matrix DX) EXACTLY AS IN LEVITT
##############################################################################
Dyear_rec <- recode(Dyear)
TD        <- dummyvar(Dyear_rec)
DX        <- cbind(Dxx, TD)
inv_DX    <- solve(t(DX) %*% DX)
DXMX      <- DX %*% inv_DX %*% t(DX)

Dxv <- Dxviol - DXMX %*% Dxviol
Dyv <- Dyviol - DXMX %*% Dyviol
Dxp <- Dxprop - DXMX %*% Dxprop
Dyp <- Dyprop - DXMX %*% Dyprop
Dxm <- Dxmurd - DXMX %*% Dxmurd
Dym <- Dymurd - DXMX %*% Dymurd

##############################################################################
# 4) OLS in first differences (just to match the Levitt printed table)
##############################################################################
Dbv <- solve(t(Dxv) %*% Dxv, t(Dxv) %*% Dyv)
Dev <- Dyv - Dxv %*% Dbv

Dbp <- solve(t(Dxp) %*% Dxp, t(Dxp) %*% Dyp)
Dep <- Dyp - Dxp %*% Dbp

Dbm <- solve(t(Dxm) %*% Dxm, t(Dxm) %*% Dym)
Dem <- Dym - Dxm %*% Dbm

# (We won't replicate the printing of standard errors here, but you can.)

##############################################################################
# 5) BUILD THE LARGE Z MATRIX (levitt style) 
##############################################################################
#   scale some columns: prison -> /10000, etc.
#   partial out time dummies from them, etc.

# dummyvar(state)
sdum <- dummyvar(state)
PS   <- solve(t(sdum) %*% sdum) %*% t(sdum)   # (sdum'*sdum)\sdum'

kxx <- ncol(Dxx)

# time dummies for Dtime
TDums <- dummyvar(Dtime)

# per the Matlab code
Dxx[,6] <- Dxx[,6]/10000
xxL[,6] <- xxL[,6]/10000
xx[,6]  <- xx[,6]/10000

Dxx[,5] <- Dxx[,5]/100
xxL[,5] <- xxL[,5]/100
xx[,5]  <- xx[,5]/100

Dxx[,4] <- Dxx[,4]/100
xxL[,4] <- xxL[,4]/100
xx[,4]  <- xx[,4]/100

Dxx[,8] <- Dxx[,8]/100
xxL[,8] <- xxL[,8]/100
xx[,8]  <- xx[,8]/100

InitZ <- Dxx
Dtime_scaled <- Dtime / DT

# Then the big 'Z' with many polynomial interactions, partial out time dummies, etc.
# We'll replicate the original logic but skip the entire LASSO portion.

# (Same code as you had for building 'Z' originally)
# ...
# We'll just keep the final 'Z' we built.
# Then partial out the basic time dummies from 'Zviol','Zprop','Zmurd' if you like.

# For simplicity, let's store all the differenced variables + the big 'Z' in a final DF.
##############################################################################
# 6) Build a final data frame that has (N*(T-1)) rows:
#    - Dstate, Dyear, etc.
#    - Dyviol, Dxviol, Dyprop, Dxprop, Dymurd, Dxmurd
#    - the Dxx columns
#    - possibly the big "Z" matrix if you want
##############################################################################

# We'll build the big "Z" separately or incorporate. For brevity, let's store:
#   - Dstate, Ddistrict, Dtime
#   - Dy's, Dx's
#   - The entire Dxx
# Then you can also store the bigZ if you prefer (it can be huge).

# We won't do every single polynomial expansion from your script, but you can append them if you want.
# For demonstration, let's do "Z = Dxx" plus some partial out if needed.

# If you truly want the enormous final matrix (with all expansions),
# then just re-assign the result to a large object and cbind it.

# Let's just store the differenced data, the partialed-out data, and the "Dxx" in the final CSV.

df_final <- data.frame(
  Dstate     = Dstate,
  Ddistrict  = Ddistrict,
  Dtime      = Dtime,
  Dy_viol    = as.numeric(Dyviol),
  Dx_viol    = as.numeric(Dxviol),
  Dy_prop    = as.numeric(Dyprop),
  Dx_prop    = as.numeric(Dxprop),
  Dy_murd    = as.numeric(Dymurd),
  Dx_murd    = as.numeric(Dxmurd)
)

# Also store the Dxx columns
colnames_Dxx <- paste0("Dxx_", seq_len(ncol(Dxx)))
df_Dxx <- as.data.frame(Dxx)
colnames(df_Dxx) <- colnames_Dxx

# Combine them
df_final <- cbind(df_final, df_Dxx)

# If you also want the partialed-out versions (Dxv, Dyv, etc.), you can cbind them:
df_final$Dxv_levitt <- as.numeric(Dxv)
df_final$Dyv_levitt <- as.numeric(Dyv)
df_final$Dxp_levitt <- as.numeric(Dxp)
df_final$Dyp_levitt <- as.numeric(Dyp)
df_final$Dxm_levitt <- as.numeric(Dxm)
df_final$Dym_levitt <- as.numeric(Dym)

# If you *also* want the big "Z" expansions (like "Zviol","Zprop","Zmurd"),
# you can do so. For example:
#   colnames_Zviol <- paste0("Zviol_", seq_len(ncol(Zviol)))
#   df_Zviol       <- as.data.frame(Zviol)
#   colnames(df_Zviol) <- colnames_Zviol
#   # Then cbind(...) to df_final. That can be hundreds of columns, fyi.

# For demonstration, let's skip it or do it selectively.

##############################################################################
# 7) Save final data to CSV
##############################################################################

write.csv(df_final, file="levitt_preprocessed_data.csv", row.names=FALSE)

cat("Saved final differenced data + partial-out columns to 'levitt_preprocessed_data.csv'.\n")
