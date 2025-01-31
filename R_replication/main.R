source("auxiliary.R")

# We'll need a norminv-like function:
norminv <- function(p) { qnorm(p) }

##############################################################################
# 2) MAIN SCRIPT replicating the Matlab code
##############################################################################
#  NOTE: Adjust path/file name if needed
data_raw <- read.table("../BCH-2014-replication/ReplicationFilesMS17285/levitt_ex.dat", sep="\t", skip=1, header=FALSE)

#  For consistency with the Matlab indexing:
#    Column 1: state
#    Column 2: year
#    Column 3: population
#    Column 4: y_viol
#    Column 5: y_prop
#    Column 6: y_murd
#    Column 7: x_murd
#    Column 8: x_viol
#    Column 9: x_prop
#    Columns 10..end: xx (other controls)
ncol_data <- ncol(data_raw)
colnames(data_raw) <- c(
  "state","year","pop","y_viol","y_prop","y_murd",
  "x_murd","x_viol","x_prop",
  paste0("X", 10:ncol_data)
)

# Remove DC (state=9), Alaska(2), Hawaii(12)
tmp <- subset(data_raw, state != 9)
tmp <- subset(tmp, state != 2)
tmp <- subset(tmp, state != 12)

# Keep only years in [85..97]
tmp <- subset(tmp, year >= 85 & year <= 97)

# Adjust year => year - 84
tmp$year <- tmp$year - 84

# recode state
tmp$state <- recode(tmp$state)

# Store variables
y_viol <- tmp$y_viol
y_prop <- tmp$y_prop
y_murd <- tmp$y_murd
x_murd <- tmp$x_murd
x_viol <- tmp$x_viol
x_prop <- tmp$x_prop
xx     <- as.matrix(tmp[, 10:ncol_data])  # col 10..end

state <- tmp$state
year  <- tmp$year
pop   <- tmp$pop  # not heavily used but in original code for weighting
N     <- max(state)
T     <- max(year)

# Build "district" exactly as in code
district <- rep(0, length(state))
# map from the Matlab code: district(state == 1) = 8, etc.
# but note we have re-labeled states, so be careful that recode() changes IDs.
# We'll do the same approach:
district_vec <- rep(0, N)
# the lines in Matlab: 
#   district(state == 1) = 8; etc.
# translate them literally but referencing the *re-labeled* state:
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
# 2a) Compute First Differences: reshape(...,T,N), diff along rows, flatten
##############################################################################
# helper to replicate: reshape(x, T, N) in column-major => diff => flatten
matDiff <- function(x, T, N) {
  # In Matlab: reshape(x, T, N) => T-by-N (columns are states)
  # Then diff along dim=1 => (T-1)-by-N
  # Then reshape => N*(T-1)-by-1 column vector
  # In R, 'matrix(x, nrow=T, ncol=N, byrow=FALSE)' also fills in column-major.
  Xmat <- matrix(x, nrow=T, ncol=N, byrow=FALSE)
  Xdiff <- apply(Xmat, 2, function(col) diff(col))
  as.vector(Xdiff)  # flatten column-wise
}

Dyviol <- matDiff(y_viol, T, N)
Dxviol <- matDiff(x_viol, T, N)
Dyprop <- matDiff(y_prop, T, N)
Dxprop <- matDiff(x_prop, T, N)
Dymurd <- matDiff(y_murd, T, N)
Dxmurd <- matDiff(x_murd, T, N)

Dxx <- matrix(0, nrow=N*(T-1), ncol=ncol(xx))
for (ii in seq_len(ncol(xx))) {
  Dxx[, ii] <- matDiff(xx[, ii], T, N)
}

# We also define t1 = find(year > 1). In Matlab, year is the vector (for each obs).
# But here we still have year for each obs. We want to keep only obs with year>1.
# Because the difference from year=1 to year=2 is the first difference, etc.
valid_idx <- which(year > 1)
Dyear <- year[valid_idx]
Dstate <- state[valid_idx]
Ddistrict <- district[valid_idx]

# recode Dyear => Dtime
Dtime <- recode(Dyear)
DT <- max(Dtime)

# We also define xxL = xx(year < T, :) => those obs where original year < T
xxL <- xx[ year < T, , drop=FALSE]

##############################################################################
# 2b) Build DX = [Dxx dummyvar(recode(Dyear))]
##############################################################################
Dyear_rec <- recode(Dyear)
TD <- dummyvar(Dyear_rec) # time dummies
DX <- cbind(Dxx, TD)

Dk <- ncol(DX) + 1

# Partial out control variables
#  DXMX = DX/(DX'*DX); in matrix form => H = X*(X'X)^(-1)*X'
inv_DX <- solve(t(DX) %*% DX)
DXMX <- DX %*% inv_DX %*% t(DX)

Dxv <- Dxviol - DXMX %*% Dxviol
Dyv <- Dyviol - DXMX %*% Dyviol
Dxp <- Dxprop - DXMX %*% Dxprop
Dyp <- Dyprop - DXMX %*% Dyprop
Dxm <- Dxmurd - DXMX %*% Dxmurd
Dym <- Dymurd - DXMX %*% Dymurd

# OLS in first differences
Dbv <- solve(t(Dxv) %*% Dxv, t(Dxv) %*% Dyv)
Dev <- Dyv - Dxv %*% Dbv

Dbp <- solve(t(Dxp) %*% Dxp, t(Dxp) %*% Dyp)
Dep <- Dyp - Dxp %*% Dbp

Dbm <- solve(t(Dxm) %*% Dxm, t(Dxm) %*% Dym)
Dem <- Dym - Dxm %*% Dbm

# Clustered SE by state
XpXinv_v <- solve(t(Dxv) %*% Dxv)
Dsv <- cluster_se(cbind(Dxv), Dev, XpXinv_v, Dstate, Dk)

XpXinv_p <- solve(t(Dxp) %*% Dxp)
Dsp <- cluster_se(cbind(Dxp), Dep, XpXinv_p, Dstate, Dk)

XpXinv_m <- solve(t(Dxm) %*% Dxm)
Dsm <- cluster_se(cbind(Dxm), Dem, XpXinv_m, Dstate, Dk)

cat("Difference Violence - Level\n")
print(c(Dbv[1], Dsv[1]))
cat("\nDifference Property - Level\n")
print(c(Dbp[1], Dsp[1]))
cat("\nDifference Murder - Level\n")
print(c(Dbm[1], Dsm[1]))
cat("\n")

# ACF
acfM <- acfcomp2(Dem, 20)
acfP <- acfcomp2(Dep, 20)
acfV <- acfcomp2(Dev, 20)
cat("Autocorrelation 1-3\n")
print(rbind(acfV[2:4], acfP[2:4], acfM[2:4]))
cat("\n")

##############################################################################
# 3) Make a big set of potential control variables (Z, etc.)
##############################################################################
# replicate:
#   sdum   = dummyvar(state)
sdum <- dummyvar(state)
PS <- solve(t(sdum) %*% sdum) %*% t(sdum)  # (sdum'*sdum)\sdum'

kxx <- ncol(Dxx)

# time dummies for Dtime
TDums <- dummyvar(Dtime)

# The code then does scaling:
#   Dxx(:,6) = Dxx(:,6)/10000; xxL(:,6) = xxL(:,6)/10000; xx(:,6) = xx(:,6)/10000
#   Dxx(:,5) = Dxx(:,5)/100; ...
#   Dxx(:,4) = Dxx(:,4)/100; ...
#   Dxx(:,8) = Dxx(:,8)/100; ...
# etc.
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

# The giant Z matrix (Z = [InitZ xxL xxL^2 ...] exactly as in the code)
Z <- cbind(
  InitZ,
  xxL,
  xxL^2,
  (Dtime_scaled*1)[,drop=FALSE] %*% t(rep(1,kxx)) * Dxx,
  (Dtime_scaled^2)[,drop=FALSE] %*% t(rep(1,kxx)) * Dxx,
  Dxx^2,
  (Dtime_scaled*1)[,drop=FALSE] %*% t(rep(1,kxx)) * (Dxx^2),
  (Dtime_scaled^2)[,drop=FALSE] %*% t(rep(1,kxx)) * (Dxx^2),
  # Next set: (Dxx(:,1)*ones(1,kxx-1)).*Dxx(:,2:8) etc. 
  # We'll replicate each block carefully:
  #
  #  (Dxx[,1]*ones(1,kxx-1)) * Dxx[,2:8]
  #  means for each row, multiply Dxx[row,1] * Dxx[row, col] for col in [2..8].
  #
  # We'll systematically build these to match the code's pattern.
  #
  # -- block1: Dprison*Dpolice, Dprison*Dur, ...
  #    i.e. pairwise products of columns 1 with columns 2..8
  (Dxx[,1,drop=FALSE] %*% matrix(1, nrow=1, ncol=kxx-1)) * Dxx[,2:kxx,drop=FALSE],
  #  (Dxx[,2]*ones(1,kxx-2)).*Dxx(:,3:8)
  (Dxx[,2,drop=FALSE] %*% matrix(1, nrow=1, ncol=kxx-2)) * Dxx[,3:kxx,drop=FALSE],
  (Dxx[,3,drop=FALSE] %*% matrix(1, nrow=1, ncol=kxx-3)) * Dxx[,4:kxx,drop=FALSE],
  (Dxx[,4,drop=FALSE] %*% matrix(1, nrow=1, ncol=kxx-4)) * Dxx[,5:kxx,drop=FALSE],
  (Dxx[,5,drop=FALSE] %*% matrix(1, nrow=1, ncol=kxx-5)) * Dxx[,6:kxx,drop=FALSE],
  (Dxx[,6,drop=FALSE] %*% matrix(1, nrow=1, ncol=kxx-6)) * Dxx[,7:kxx,drop=FALSE],
  (Dxx[,7,drop=FALSE] %*% matrix(1, nrow=1, ncol=kxx-7)) * Dxx[,8:kxx,drop=FALSE],
  #
  # plus the same but multiplied by Dtime, Dtime^2, etc.
  ((Dxx[,1]*Dtime_scaled) %*% matrix(1, nrow=1, ncol=kxx-1)) * Dxx[,2:kxx],
  ((Dxx[,2]*Dtime_scaled) %*% matrix(1, nrow=1, ncol=kxx-2)) * Dxx[,3:kxx],
  ((Dxx[,3]*Dtime_scaled) %*% matrix(1, nrow=1, ncol=kxx-3)) * Dxx[,4:kxx],
  ((Dxx[,4]*Dtime_scaled) %*% matrix(1, nrow=1, ncol=kxx-4)) * Dxx[,5:kxx],
  ((Dxx[,5]*Dtime_scaled) %*% matrix(1, nrow=1, ncol=kxx-5)) * Dxx[,6:kxx],
  ((Dxx[,6]*Dtime_scaled) %*% matrix(1, nrow=1, ncol=kxx-6)) * Dxx[,7:kxx],
  ((Dxx[,7]*Dtime_scaled) %*% matrix(1, nrow=1, ncol=kxx-7)) * Dxx[,8:kxx],
  ((Dxx[,1]*(Dtime_scaled^2)) %*% matrix(1, nrow=1, ncol=kxx-1)) * Dxx[,2:kxx],
  ((Dxx[,2]*(Dtime_scaled^2)) %*% matrix(1, nrow=1, ncol=kxx-2)) * Dxx[,3:kxx],
  ((Dxx[,3]*(Dtime_scaled^2)) %*% matrix(1, nrow=1, ncol=kxx-3)) * Dxx[,4:kxx],
  ((Dxx[,4]*(Dtime_scaled^2)) %*% matrix(1, nrow=1, ncol=kxx-4)) * Dxx[,5:kxx],
  ((Dxx[,5]*(Dtime_scaled^2)) %*% matrix(1, nrow=1, ncol=kxx-5)) * Dxx[,6:kxx],
  ((Dxx[,6]*(Dtime_scaled^2)) %*% matrix(1, nrow=1, ncol=kxx-6)) * Dxx[,7:kxx],
  ((Dxx[,7]*(Dtime_scaled^2)) %*% matrix(1, nrow=1, ncol=kxx-7)) * Dxx[,8:kxx],
  #
  # kron(Dxx(Dyear-1 == 1,[1:6,8]), ones(DT,1)) ...
  # We mimic exactly. First we find rows where (Dyear-1==1) => Dyear==2
  # Then replicate (kron) each row across DT times
  # We'll define a small helper:
  {
    rowidx <- which(Dyear == 2)
    blockA <- Dxx[rowidx, c(1:6,8), drop=FALSE]
    # kron(blockA, ones(DT,1)) => replicate each row of blockA vertically DT times
    # A standard approach in R: do.call(rbind, replicate(DT, blockA, simplify=FALSE))
    do.call(rbind, replicate(DT, blockA, simplify=FALSE))
  },
  # same ^2
  {
    rowidx <- which(Dyear == 2)
    blockA <- Dxx[rowidx, c(1:6,8), drop=FALSE]^2
    do.call(rbind, replicate(DT, blockA, simplify=FALSE))
  },
  # kron(xx(year ==1,:), ones(DT,1)), etc.
  {
    rowidx <- which(year == 1)
    blockB <- xx[rowidx,,drop=FALSE]
    do.call(rbind, replicate(DT, blockB, simplify=FALSE))
  },
  {
    rowidx <- which(year == 1)
    blockB <- xx[rowidx,,drop=FALSE]^2
    do.call(rbind, replicate(DT, blockB, simplify=FALSE))
  },
  # then .*(Dtime) or .*(Dtime.^2), etc.
  {
    rowidx <- which(Dyear == 2)
    blockA <- Dxx[rowidx, c(1:6,8), drop=FALSE]
    bigA   <- do.call(rbind, replicate(DT, blockA, simplify=FALSE))
    bigA * (Dtime_scaled)
  },
  {
    rowidx <- which(Dyear == 2)
    blockA <- (Dxx[rowidx, c(1:6,8), drop=FALSE])^2
    bigA   <- do.call(rbind, replicate(DT, blockA, simplify=FALSE))
    bigA * (Dtime_scaled)
  },
  {
    rowidx <- which(year == 1)
    blockB <- xx[rowidx,,drop=FALSE]
    bigB   <- do.call(rbind, replicate(DT, blockB, simplify=FALSE))
    bigB * (Dtime_scaled)
  },
  {
    rowidx <- which(year == 1)
    blockB <- xx[rowidx,,drop=FALSE]^2
    bigB   <- do.call(rbind, replicate(DT, blockB, simplify=FALSE))
    bigB * (Dtime_scaled)
  },
  {
    rowidx <- which(Dyear == 2)
    blockA <- Dxx[rowidx, c(1:6,8), drop=FALSE]
    bigA   <- do.call(rbind, replicate(DT, blockA, simplify=FALSE))
    bigA * (Dtime_scaled^2)
  },
  {
    rowidx <- which(Dyear == 2)
    blockA <- (Dxx[rowidx, c(1:6,8), drop=FALSE])^2
    bigA   <- do.call(rbind, replicate(DT, blockA, simplify=FALSE))
    bigA * (Dtime_scaled^2)
  },
  {
    rowidx <- which(year == 1)
    blockB <- xx[rowidx,,drop=FALSE]
    bigB   <- do.call(rbind, replicate(DT, blockB, simplify=FALSE))
    bigB * (Dtime_scaled^2)
  },
  {
    rowidx <- which(year == 1)
    blockB <- xx[rowidx,,drop=FALSE]^2
    bigB   <- do.call(rbind, replicate(DT, blockB, simplify=FALSE))
    bigB * (Dtime_scaled^2)
  },
  #   kron(PS*xx, ones(DT,1)), (PS*xx)^2, etc.
  {
    px <- PS %*% xx
    bigC <- do.call(rbind, replicate(DT, px, simplify=FALSE))
    bigC
  },
  {
    px <- (PS %*% xx)^2
    bigC <- do.call(rbind, replicate(DT, px, simplify=FALSE))
    bigC
  },
  {
    px <- PS %*% xx
    bigC <- do.call(rbind, replicate(DT, px, simplify=FALSE))
    bigC * (Dtime_scaled)
  },
  {
    px <- PS %*% xx
    bigC <- do.call(rbind, replicate(DT, px, simplify=FALSE))
    bigC * (Dtime_scaled^2)
  },
  {
    px <- (PS %*% xx)^2
    bigC <- do.call(rbind, replicate(DT, px, simplify=FALSE))
    bigC * (Dtime_scaled)
  },
  {
    px <- (PS %*% xx)^2
    bigC <- do.call(rbind, replicate(DT, px, simplify=FALSE))
    bigC * (Dtime_scaled^2)
  }
)

# Now remove time dummy means from the "InitZ":
#  InitZ = InitZ - TDums*((TDums'*TDums)\(TDums'*InitZ))
#  i.e. partial out time dummies from each column
Htd <- TDums %*% solve(t(TDums) %*% TDums) %*% t(TDums)
InitZ_adj <- InitZ - Htd %*% InitZ  # but the code sets only for NameZ, we do for the big Z?

# The Matlab actually does: "InitZ = InitZ - TDums*((TDums'*TDums)\(TDums'*InitZ));"
# Then it sets Z = [ ... ] but it never uses the updated "InitZ"?  
# Actually, in the code, "InitZ" is assigned but isn't used further except in the final line:
#   "InitZ = InitZ - TDums*((TDums'*TDums)\(TDums'*InitZ));"
# Then "Z" is formed from the original Dxx etc. 
# So let's stay consistent: We'll replicate that single line for "InitZ" but keep "Z" as built above:
InitZ <- InitZ_adj

NameZ <- c(
  "Dprison","Dpolice","Dur","Dinc","Dpov","Dafdc","Dgun","Dbeer",
  "Lprison","Lpolice","Lur","Linc","Lpov","Lafdc","Lgun","Lbeer",
  "Lprison^2","Lpolice^2","Lur^2","Linc^2","Lpov^2","Lafdc^2","Lgun^2","Lbeer^2",
  "Dprison*t","Dpolice*t","Dur*t","Dinc*t","Dpov*t","Dafdc*t","Dgun*t","Dbeer*t",
  "Dprison*t^2","Dpolice*t^2","Dur*t^2","Dinc*t^2","Dpov*t^2","Dafdc*t^2","Dgun*t^2","Dbeer*t^2",
  "Dprison^2","Dpolice^2","Dur^2","Dinc^2","Dpov^2","Dafdc^2","Dgun^2","Dbeer^2",
  "Dprison^2*t","Dpolice^2*t","Dur^2*t","Dinc^2*t","Dpov^2*t","Dafdc^2*t","Dgun^2*t","Dbeer^2*t",
  "Dprison^2*t^2","Dpolice^2*t^2","Dur^2*t^2","Dinc^2*t^2","Dpov^2*t^2","Dafdc^2*t^2","Dgun^2*t^2","Dbeer^2*t^2",
  "Dprison*Dpolice","Dprison*Dur","Dprison*Dinc","Dprison*Dpov","Dprison*Dafdc","Dprison*Dgun","Dprison*Dbeer",
  "Dpolice*Dur","Dpolice*Dinc","Dpolice*Dpov","Dpolice*Dafdc","Dpolice*Dgun","Dpolice*Dbeer",
  "Dur*Dinc","Dur*Dpov","Dur*Dafdc","Dur*Dgun","Dur*Dbeer",
  "Dinc*Dpov","Dinc*Dafdc","Dinc*Dgun","Dinc*Dbeer",
  "Dpov*Dafdc","Dpov*Dgun","Dpov*Dbeer",
  "Dafdc*Dgun","Dafdc*Dbeer",
  "Dgun*Dbeer",
  "Dprison*Dpolice*t","Dprison*Dur*t","Dprison*Dinc*t","Dprison*Dpov*t","Dprison*Dafdc*t","Dprison*Dgun*t","Dprison*Dbeer*t",
  "Dpolice*Dur*t","Dpolice*Dinc*t","Dpolice*Dpov*t","Dpolice*Dafdc*t","Dpolice*Dgun*t","Dpolice*Dbeer*t",
  "Dur*Dinc*t","Dur*Dpov*t","Dur*Dafdc*t","Dur*Dgun*t","Dur*Dbeer*t",
  "Dinc*Dpov*t","Dinc*Dafdc*t","Dinc*Dgun*t","Dinc*Dbeer*t",
  "Dpov*Dafdc*t","Dpov*Dgun*t","Dpov*Dbeer*t",
  "Dafdc*Dgun*t","Dafdc*Dbeer*t",
  "Dgun*Dbeer*t",
  "Dprison*Dpolice*t^2","Dprison*Dur*t^2","Dprison*Dinc*t^2","Dprison*Dpov*t^2","Dprison*Dafdc*t^2","Dprison*Dgun*t^2","Dprison*Dbeer*t^2",
  "Dpolice*Dur*t^2","Dpolice*Dinc*t^2","Dpolice*Dpov*t^2","Dpolice*Dafdc*t^2","Dpolice*Dgun*t^2","Dpolice*Dbeer*t^2",
  "Dur*Dinc*t^2","Dur*Dpov*t^2","Dur*Dafdc*t^2","Dur*Dgun*t^2","Dur*Dbeer*t^2",
  "Dinc*Dpov*t^2","Dinc*Dafdc*t^2","Dinc*Dgun*t^2","Dinc*Dbeer*t^2",
  "Dpov*Dafdc*t^2","Dpov*Dgun*t^2","Dpov*Dbeer*t^2",
  "Dafdc*Dgun*t^2","Dafdc*Dbeer*t^2",
  "Dgun*Dbeer*t^2",
  "Dprison0","Dpolice0","Dur0","Dinc0","Dpov0","Dafdc0","Dbeer0",
  "Dprison0^2","Dpolice0^2","Dur0^2","Dinc0^2","Dpov0^2","Dafdc0^2","Dbeer0^2",
  "Lprison0","Lpolice0","Lur0","Linc0","Lpov0","Lafdc0","Lgun0","Lbeer0",
  "Lprison0^2","Lpolice0^2","Lur0^2","Linc0^2","Lpov0^2","Lafdc0^2","Lgun0^2","Lbeer0^2",
  "Dprison0*t","Dpolice0*t","Dur0*t","Dinc0*t","Dpov0*t","Dafdc0*t","Dbeer0*t",
  "Dprison0^2*t","Dpolice0^2*t","Dur0^2*t","Dinc0^2*t","Dpov0^2*t","Dafdc0^2*t","Dbeer0^2*t",
  "Lprison0*t","Lpolice0*t","Lur0*t","Linc0*t","Lpov0*t","Lafdc0*t","Lgun0*t","Lbeer0*t",
  "Lprison0^2*t","Lpolice0^2*t","Lur0^2*t","Linc0^2*t","Lpov0^2*t","Lafdc0^2*t","Lgun0^2*t","Lbeer0^2*t",
  "Dprison0*t^2","Dpolice0*t^2","Dur0*t^2","Dinc0*t^2","Dpov0*t^2","Dafdc0*t^2","Dbeer0*t^2",
  "Dprison0^2*t^2","Dpolice0^2*t^2","Dur0^2*t^2","Dinc0^2*t^2","Dpov0^2*t^2","Dafdc0^2*t^2","Dbeer0^2*t^2",
  "Lprison0*t^2","Lpolice0*t^2","Lur0*t^2","Linc0*t^2","Lpov0*t^2","Lafdc0*t^2","Lgun0*t^2","Lbeer0*t^2",
  "Lprison0^2*t^2","Lpolice0^2*t^2","Lur0^2*t^2","Linc0^2*t^2","Lpov0^2*t^2","Lafdc0^2*t^2","Lgun0^2*t^2","Lbeer0^2*t^2",
  "prisonBar","policeBar","urBar","incBar","povBar","afdcBar","gunBar","beerBar",
  "prisonBar^2","policeBar^2","urBar^2","incBar^2","povBar^2","afdcBar^2","gunBar^2","beerBar^2",
  "prisonBar*t","policeBar*t","urBar*t","incBar*t","povBar*t","afdcBar*t","gunBar*t","beerBar*t",
  "prisonBar*t^2","policeBar*t^2","urBar*t^2","incBar*t^2","povBar*t^2","afdcBar*t^2","gunBar*t^2","beerBar*t^2",
  "prisonBar^2*t","policeBar^2*t","urBar^2*t","incBar^2*t","povBar^2*t","afdcBar^2*t","gunBar^2*t","beerBar^2*t",
  "prisonBar^2*t^2","policeBar^2*t^2","urBar^2*t^2","incBar^2*t^2","povBar^2*t^2","afdcBar^2*t^2","gunBar^2*t^2","beerBar^2*t^2"
)

##############################################################################
# Now build Zviol, Zprop, Zmurd, partial out time dummies, etc.
##############################################################################
buildZcrime <- function(Z, DxCrime, Dyear, year, x_crime, Dtime, DT, NameZcrime_label) {
  # replicate the last lines:
  #   [Z kron(DxCrime(Dyear-1==1),ones(DT,1)) ... ]
  
  cidx <- which(Dyear == 2)  # Dyear-1 == 1 => Dyear==2
  blockDX <- DxCrime[cidx]
  blockDX <- matrix(blockDX, ncol = 1)  # format as matrix 
  blockDX_mat <- do.call(rbind, replicate(DT, blockDX, simplify=FALSE))
  
  cidx2 <- which(year == 1)
  x_crime_block <- x_crime[cidx2]
  x_crime_block <- matrix(x_crime_block, ncol = 1)
  x_crime_mat <- do.call(rbind, replicate(DT, x_crime_block, simplify=FALSE))
  
  # so Z + 6 extra columns in the same structure
  #   kron(DxCrime(...)) => blockDX_mat, also squared, etc.
  #   plus multiply by Dtime, etc.
  # We'll do exactly the same pattern of 12 appended columns (the code has 6 pairs).
  
  # columns: 
  #   'DxV0','DxV0^2','DxV0*t','DxV0^2*t','DxV0*t^2','DxV0^2*t^2',
  #   'xV0','xV0^2','xV0*t','xV0^2*t','xV0*t^2','xV0^2*t^2'
  # for violence. Similarly for P or M.
  
  blockDX2 <- blockDX_mat^2
  dt <- Dtime
  dt2 <- (Dtime^2)
  
  tailZ <- cbind(
    blockDX_mat,
    blockDX2,
    blockDX_mat * dt,
    blockDX2 * dt,
    blockDX_mat * dt2,
    blockDX2 * dt2,
    #
    x_crime_mat,
    x_crime_mat^2,
    x_crime_mat * dt,
    (x_crime_mat^2) * dt,
    x_crime_mat * dt2,
    (x_crime_mat^2) * dt2
  )
  
  cat("Z: ", dim(Z), "tailZ: ", dim(tailZ), "\n")
  Z_out <- cbind(Z, tailZ)
  # build name vector
  nm_tail <- c("DxC0","DxC0^2","DxC0*t","DxC0^2*t","DxC0*t^2","DxC0^2*t^2",
               "xC0","xC0^2","xC0*t","xC0^2*t","xC0*t^2","xC0^2*t^2")
  # rename them (DxV0, etc.) by substituting the letter in "Crime"
  # but we can do a direct replacement, or we can let the user pass them
  nm_out <- c(NameZ, NameZcrime_label)
  return(list(Z=Z_out, Name=nm_out))
}

NV <- c("DxV0","DxV0^2","DxV0*t","DxV0^2*t","DxV0*t^2","DxV0^2*t^2",
        "xV0","xV0^2","xV0*t","xV0^2*t","xV0*t^2","xV0^2*t^2")
NP <- c("DxP0","DxP0^2","DxP0*t","DxP0^2*t","DxP0*t^2","DxP0^2*t^2",
        "xP0","xP0^2","xP0*t","xP0^2*t","xP0*t^2","xP0^2*t^2")
NM <- c("DxM0","DxM0^2","DxM0*t","DxM0^2*t","DxM0*t^2","DxM0^2*t^2",
        "xM0","xM0^2","xM0*t","xM0^2*t","xM0*t^2","xM0^2*t^2")

outV <- buildZcrime(Z, Dxviol, Dyear, year, x_viol, Dtime_scaled, DT, NV)
Zviol <- outV$Z
NameZviol <- outV$Name

outP <- buildZcrime(Z, Dxprop, Dyear, year, x_prop, Dtime_scaled, DT, NP)
Zprop <- outP$Z
NameZprop <- outP$Name

outM <- buildZcrime(Z, Dxmurd, Dyear, year, x_murd, Dtime_scaled, DT, NM)
Zmurd <- outM$Z
NameZmurd <- outM$Name

# Partial out time dummies:
Htd2 <- TDums %*% solve(t(TDums) %*% TDums) %*% t(TDums)
Zviol <- Zviol - Htd2 %*% Zviol
Zprop <- Zprop - Htd2 %*% Zprop
Zmurd <- Zmurd - Htd2 %*% Zmurd

DxV <- Dxviol - Htd2 %*% Dxviol
DxP <- Dxprop - Htd2 %*% Dxprop
DxM <- Dxmurd - Htd2 %*% Dxmurd
DyV <- Dyviol - Htd2 %*% Dyviol
DyP <- Dyprop - Htd2 %*% Dyprop
DyM <- Dymurd - Htd2 %*% Dymurd

maxIter <- 100

# findNonCollinear
keepZviol <- findNonCollinear(Zviol, tol = 1e-16)
Zviol <- Zviol[, keepZviol, drop=FALSE]
NameZviol <- NameZviol[keepZviol]

keepZprop <- findNonCollinear(Zprop, tol = 1e-16)
Zprop <- Zprop[, keepZprop, drop=FALSE]
NameZprop <- NameZprop[keepZprop]

keepZmurd <- findNonCollinear(Zmurd, tol = 1e-16)
Zmurd <- Zmurd[, keepZmurd, drop=FALSE]
NameZmurd <- NameZmurd[keepZmurd]

##############################################################################
# 4) LASSO - Violence
##############################################################################
DxV <- as.numeric(DxV)
StV <- Zviol * (DxV)
UpsV <- sqrt(colSums(StV^2)/(DT*N))
# lambda = 1.1*2*sqrt(DT*N)*norminv(1-.05/(2*size(Zviol,2)))
lambdaV <- 1.1*2*sqrt(DT*N)*norminv(1 - 0.05/(2*ncol(Zviol)))

PInitV <- LassoShooting2(Zviol, DxV, lambdaV, UpsV, verbose=1)
IndInitV <- abs(PInitV) > 0
ZV <- Zviol[, IndInitV, drop=FALSE]
RefResidV <- DxV - Zviol %*% PInitV
efsV <- DxV - ZV %*% solve(crossprod(ZV), crossprod(ZV, DxV))
IndRefV <- IndInitV

kk <- 1
RefResidV <- as.numeric(RefResidV)
StRefV <- Zviol * RefResidV
UpsRefV <- sqrt(colSums(StRefV^2)/(DT*N))
while (sqrt(sum((UpsRefV - UpsV)^2)) > 1e-4 && kk < maxIter) {
  cat(kk, "\n")
  PRefV <- LassoShooting2(Zviol, DxV, lambdaV, UpsRefV, verbose=1)
  IndRefV <- abs(PRefV) > 0
  ZV <- Zviol[, IndRefV, drop=FALSE]
  bfsV <- solve(crossprod(ZV), crossprod(ZV, DxV))
  efsV <- DxV - ZV %*% bfsV
  UpsV <- UpsRefV
  StRefV <- Zviol * as.numeric((DxV - Zviol %*% PRefV))
  UpsRefV <- sqrt(colSums(StRefV^2)/(DT*N))
  kk <- kk+1
}
sfsV <- cluster_se(ZV, efsV, solve(crossprod(ZV)), Dstate, ncol(ZV))

# Structural eq
DyV <- as.numeric(DyV)
StyV <- Zviol * DyV
UpsyV <- sqrt(colSums(StyV^2)/(DT*N))
PInityV <- LassoShooting2(Zviol, DyV, lambdaV, UpsyV, verbose=1)
IndInityV <- abs(PInityV) > 0
ZyV <- Zviol[, IndInityV, drop=FALSE]
RefResidyV <- DyV - Zviol %*% PInityV
efsyV <- DyV - ZyV %*% solve(crossprod(ZyV), crossprod(ZyV, DyV)) # Note: This throws an error, since ZyV is zero dimensional: lasso throws them all out.
IndRefyV <- IndInityV

kk <- 1
RefResidyV <- as.numeric(RefResidyV)
StRefyV <- Zviol * RefResidyV
UpsRefyV <- sqrt(colSums(StRefyV^2)/(DT*N))
while (sqrt(sum((UpsRefyV - UpsyV)^2)) > 1e-4 && kk < maxIter) {
  cat(kk, "\n")
  PRefyV <- LassoShooting2(Zviol, DyV, lambdaV, UpsRefyV, verbose=1)
  IndRefyV <- abs(PRefyV) > 0
  ZyV <- Zviol[, IndRefyV, drop=FALSE]
  bfsyV <- solve(crossprod(ZyV), crossprod(ZyV, DyV))
  efsyV <- DyV - ZyV %*% bfsyV
  UpsyV <- UpsRefyV
  StRefyV <- Zviol * as.numeric((DyV - Zviol %*% PRefyV))
  UpsRefyV <- sqrt(colSums(StRefyV^2)/(DT*N))
  kk <- kk + 1
}
sfsyV <- cluster_se(ZyV, efsyV, solve(crossprod(ZyV)), Dstate, ncol(ZyV))

IndUnionV <- pmax(IndRefyV, IndRefV)
ZUnionV <- Zviol[, IndUnionV>0, drop=FALSE]

R2_1V <- sum((ZUnionV %*% solve(crossprod(ZUnionV), crossprod(ZUnionV, DxV)))^2)/sum(DxV^2)
R2_2V <- sum((ZUnionV %*% solve(crossprod(ZUnionV), crossprod(ZUnionV, DyV)))^2)/sum(DyV^2)

ZviolLASSO <- cbind(DxV, ZUnionV)
bviolLASSO <- solve(crossprod(ZviolLASSO), crossprod(ZviolLASSO, DyV))
eviolLASSO <- DyV - ZviolLASSO %*% bviolLASSO
ResviolLASSO <- cluster_se(ZviolLASSO, eviolLASSO, solve(crossprod(ZviolLASSO)), Dstate, ncol(ZviolLASSO))

##############################################################################
# 4b) LASSO - Property
##############################################################################
DxP <- as.numeric(DxP)
StP <- Zprop * DxP
UpsP <- sqrt(colSums(StP^2)/(DT*N))
lambdaP <- 1.1*2*sqrt(DT*N)*norminv(1 - 0.05/(2*ncol(Zprop)))
PInitP <- LassoShooting2(Zprop, DxP, lambdaP, UpsP, verbose=1)
IndInitP <- abs(PInitP) > 0
ZP <- Zprop[, IndInitP, drop=FALSE]
RefResidP <- DxP - Zprop %*% PInitP
efsP <- DxP - ZP %*% solve(crossprod(ZP), crossprod(ZP, DxP))
IndRefP <- IndInitP

kk <- 1
RefResidP <- as.numeric(RefResidP)
StRefP <- Zprop * RefResidP
UpsRefP <- sqrt(colSums(StRefP^2)/(DT*N))
while (sqrt(sum((UpsRefP - UpsP)^2)) > 1e-4 && kk < maxIter) {
  cat(kk, "\n")
  PRefP <- LassoShooting2(Zprop, DxP, lambdaP, UpsRefP, verbose=1)
  IndRefP <- abs(PRefP) > 0
  ZP <- Zprop[, IndRefP, drop=FALSE]
  bfsP <- solve(crossprod(ZP), crossprod(ZP, DxP))
  efsP <- DxP - ZP %*% bfsP
  UpsP <- UpsRefP
  StRefP <- Zprop * as.numeric((DxP - Zprop %*% PRefP))
  UpsRefP <- sqrt(colSums(StRefP^2)/(DT*N))
  kk <- kk+1
}
sfsP <- cluster_se(ZP, efsP, solve(crossprod(ZP)), Dstate, ncol(ZP))

# Structural eq
DyP <- as.numeric(DyP)
StyP <- Zprop * DyP
UpsyP <- sqrt(colSums(StyP^2)/(DT*N))
PInityP <- LassoShooting2(Zprop, DyP, lambdaP, UpsyP, verbose=1)
IndInityP <- abs(PInityP) > 0
ZyP <- Zprop[, IndInityP, drop=FALSE]
RefResidyP <- DyP - Zprop %*% PInityP
efsyP <- DyP - ZyP %*% solve(crossprod(ZyP), crossprod(ZyP, DyP))
IndRefyP <- IndInityP

kk <- 1
RefResidyP <- as.numeric(RefResidyP)
StRefyP <- Zprop * RefResidyP
UpsRefyP <- sqrt(colSums(StRefyP^2)/(DT*N))
while (sqrt(sum((UpsRefyP - UpsyP)^2)) > 5e-4 && kk < maxIter) {
  cat(kk, "\n")
  PRefyP <- LassoShooting2(Zprop, DyP, lambdaP, UpsRefyP, verbose=1)
  IndRefyP <- abs(PRefyP) > 0
  ZyP <- Zprop[, IndRefyP, drop=FALSE]
  bfsyP <- solve(crossprod(ZyP), crossprod(ZyP, DyP))
  efsyP <- DyP - ZyP %*% bfsyP
  UpsyP <- UpsRefyP
  StRefyP <- Zprop * as.numeric((DyP - Zprop %*% PRefyP))
  UpsRefyP <- sqrt(colSums(StRefyP^2)/(DT*N))
  kk <- kk+1
}
sfsyP <- cluster_se(ZyP, efsyP, solve(crossprod(ZyP)), Dstate, ncol(ZyP))

IndUnionP <- pmax(IndRefyP, IndRefP)
ZUnionP <- Zprop[, IndUnionP>0, drop=FALSE]
R2_1P <- sum((ZUnionP %*% solve(crossprod(ZUnionP), crossprod(ZUnionP, DxP)))^2)/sum(DxP^2)
R2_2P <- sum((ZUnionP %*% solve(crossprod(ZUnionP), crossprod(ZUnionP, DyP)))^2)/sum(DyP^2)

ZpropLASSO <- cbind(DxP, ZUnionP)
bpropLASSO <- solve(crossprod(ZpropLASSO), crossprod(ZpropLASSO, DyP))
epropLASSO <- DyP - ZpropLASSO %*% bpropLASSO
RespropLASSO <- cluster_se(ZpropLASSO, epropLASSO, solve(crossprod(ZpropLASSO)), Dstate, ncol(ZpropLASSO))

##############################################################################
# 4c) LASSO - Murder
##############################################################################
DxM <- as.numeric(DxM)
StM <- Zmurd * DxM
UpsM <- sqrt(colSums(StM^2)/(DT*N))
lambdaM <- 1.1*2*sqrt(DT*N)*norminv(1 - 0.05/(2*ncol(Zmurd)))
PInitM <- LassoShooting2(Zmurd, DxM, lambdaM, UpsM, verbose=1)
IndInitM <- abs(PInitM) > 0
ZM <- Zmurd[, IndInitM, drop=FALSE]
RefResidM <- DxM - Zmurd %*% PInitM
efsM <- DxM - ZM %*% solve(crossprod(ZM), crossprod(ZM, DxM))
IndRefM <- IndInitM

kk <- 1
RefResidM <- as.numeric(RefResidM)
StRefM <- Zmurd * RefResidM
UpsRefM <- sqrt(colSums(StRefM^2)/(DT*N))
while (sqrt(sum((UpsRefM - UpsM)^2)) > 1e-4 && kk < maxIter) {
  cat(kk, "\n")
  PRefM <- LassoShooting2(Zmurd, DxM, lambdaM, UpsRefM, verbose=1)
  IndRefM <- abs(PRefM) > 0
  ZM <- Zmurd[, IndRefM, drop=FALSE]
  bfsM <- solve(crossprod(ZM), crossprod(ZM, DxM))
  efsM <- DxM - ZM %*% bfsM
  UpsM <- UpsRefM
  StRefM <- Zmurd * as.numeric((DxM - Zmurd %*% PRefM))
  UpsRefM <- sqrt(colSums(StRefM^2)/(DT*N))
  kk <- kk+1
}
sfsM <- cluster_se(ZM, efsM, solve(crossprod(ZM)), Dstate, ncol(ZM))

# Structural eq
DyM <- as.numeric(DyM)
StyM <- Zmurd * DyM
UpsyM <- sqrt(colSums(StyM^2)/(DT*N))
PInityM <- LassoShooting2(Zmurd, DyM, lambdaM, UpsyM, verbose=1)
IndInityM <- abs(PInityM) > 0
ZyM <- Zmurd[, IndInityM, drop=FALSE]
RefResidyM <- DyM - Zmurd %*% PInityM
efsyM <- DyM - ZyM %*% solve(crossprod(ZyM), crossprod(ZyM, DyM))
IndRefyM <- IndInityM

kk <- 1
RefResidyM <- as.numeric(RefResidyM)
StRefyM <- Zmurd * RefResidyM
UpsRefyM <- sqrt(colSums(StRefyM^2)/(DT*N))
while (sqrt(sum((UpsRefyM - UpsyM)^2)) > 1e-4 && kk < maxIter) {
  cat(kk, "\n")
  PRefyM <- LassoShooting2(Zmurd, DyM, lambdaM, UpsRefyM, verbose=1)
  IndRefyM <- abs(PRefyM) > 0
  ZyM <- Zmurd[, IndRefyM, drop=FALSE]
  bfsyM <- solve(crossprod(ZyM), crossprod(ZyM, DyM))
  efsyM <- DyM - ZyM %*% bfsyM
  UpsyM <- UpsRefyM
  StRefyM <- Zmurd * as.numeric((DyM - Zmurd %*% PRefyM))
  UpsRefyM <- sqrt(colSums(StRefyM^2)/(DT*N))
  kk <- kk+1
}
sfsyM <- cluster_se(ZyM, efsyM, solve(crossprod(ZyM)), Dstate, ncol(ZyM))

IndUnionM <- pmax(IndRefyM, IndRefM)
ZUnionM <- Zmurd[, IndUnionM>0, drop=FALSE]
R2_1M <- sum((ZUnionM %*% solve(crossprod(ZUnionM), crossprod(ZUnionM, DxM)))^2)/sum(DxM^2)
R2_2M <- sum((ZUnionM %*% solve(crossprod(ZUnionM), crossprod(ZUnionM, DyM)))^2)/sum(DyM^2)

ZmurdLASSO <- cbind(DxM, ZUnionM)
bmurdLASSO <- solve(crossprod(ZmurdLASSO), crossprod(ZmurdLASSO, DyM))
emurdLASSO <- DyM - ZmurdLASSO %*% bmurdLASSO
ResmurdLASSO <- cluster_se(ZmurdLASSO, emurdLASSO, solve(crossprod(ZmurdLASSO)), Dstate, ncol(ZmurdLASSO))

acfML <- acfcomp2(emurdLASSO, 20)
acfPL <- acfcomp2(epropLASSO, 20)
acfVL <- acfcomp2(eviolLASSO, 20)

cat("Difference Violence - Level - LASSO\n")
cat("Selected Variables - Abortion\n")
cat(paste(NameZviol[IndRefV], collapse=", "),"\n\n")
cat("Selected Variables - Crime\n")
cat(paste(NameZviol[IndRefyV], collapse=", "),"\n\n")
cat("Crime equation\n")
cat(c(bviolLASSO[1], ResviolLASSO[1]), "\n\n")

cat("Difference Property - Level - LASSO\n")
cat("Selected Variables - Abortion\n")
cat(paste(NameZprop[IndRefP], collapse=", "),"\n\n")
cat("Selected Variables - Crime\n")
cat(paste(NameZprop[IndRefyP], collapse=", "),"\n\n")
cat("Crime equation\n")
cat(c(bpropLASSO[1], RespropLASSO[1]), "\n\n")

cat("Difference Murder - Level - LASSO\n")
cat("Selected Variables - Abortion\n")
cat(paste(NameZmurd[IndRefM], collapse=", "),"\n\n")
cat("Selected Variables - Crime\n")
cat(paste(NameZmurd[IndRefyM], collapse=", "),"\n\n")
cat("Crime equation\n")
cat(c(bmurdLASSO[1], ResmurdLASSO[1]), "\n\n")

cat("First Order Autocorrelation\n")
cat(c(acfVL[2], acfPL[2], acfML[2]),"\n\n")

##############################################################################
# 5) Estimates with Donohue & Levitt variables in amelioration set
##############################################################################
IndUnionV <- IndUnionV
IndUnionP <- IndUnionP
IndUnionM <- IndUnionM

# The code forcibly sets first 8 indexes to 1
IndUnionV[1:8] <- 1
IndUnionP[1:8] <- 1
IndUnionM[1:8] <- 1

ZUnionV <- Zviol[, IndUnionV>0, drop=FALSE]
ZUnionP <- Zprop[, IndUnionP>0, drop=FALSE]
# The code has "ZUnionM = Zprop(:,IndUnionM);"? Possibly a bug in the original code.
# We replicate exactly the line in the Matlab script:
#  "ZUnionM = Zprop(:,IndUnionM);" 
#  but that looks suspicious. Probably it should be "Zmurd"?
# The given code says:
#   ZUnionM = Zprop(:,IndUnionM)
# That is presumably a small slip in the original. We'll do exactly what it says:
ZUnionM <- Zprop[, IndUnionM>0, drop=FALSE]

ZviolLASSO2 <- cbind(DxV, ZUnionV)
bviolLASSO2 <- solve(crossprod(ZviolLASSO2), crossprod(ZviolLASSO2, DyV))
eviolLASSO2 <- DyV - ZviolLASSO2 %*% bviolLASSO2
ResviolLASSO2 <- cluster_se(ZviolLASSO2, eviolLASSO2, solve(crossprod(ZviolLASSO2)), Dstate, ncol(ZviolLASSO2))

ZpropLASSO2 <- cbind(DxP, ZUnionP)
bpropLASSO2 <- solve(crossprod(ZpropLASSO2), crossprod(ZpropLASSO2, DyP))
epropLASSO2 <- DyP - ZpropLASSO2 %*% bpropLASSO2
RespropLASSO2 <- cluster_se(ZpropLASSO2, epropLASSO2, solve(crossprod(ZpropLASSO2)), Dstate, ncol(ZpropLASSO2))

ZmurdLASSO2 <- cbind(DxM, ZUnionM)
bmurdLASSO2 <- solve(crossprod(ZmurdLASSO2), crossprod(ZmurdLASSO2, DyM))
emurdLASSO2 <- DyM - ZmurdLASSO2 %*% bmurdLASSO2
ResmurdLASSO2 <- cluster_se(ZmurdLASSO2, emurdLASSO2, solve(crossprod(ZmurdLASSO2)), Dstate, ncol(ZmurdLASSO2))

acfVLp <- acfcomp2(eviolLASSO2, 20)
acfPLp <- acfcomp2(epropLASSO2, 20)
acfMLp <- acfcomp2(emurdLASSO2, 20)

cat("Difference Violence - LASSO + Original\n")
cat("Crime equation\n")
cat(c(bviolLASSO2[1], ResviolLASSO2[1]), "\n\n")

cat("Difference Property - LASSO + Original\n")
cat("Crime equation\n")
cat(c(bpropLASSO2[1], RespropLASSO2[1]), "\n\n")

cat("Difference Murder - LASSO + Original\n")
cat("Crime equation\n")
cat(c(bmurdLASSO2[1], ResmurdLASSO2[1]), "\n\n")

cat("Autocorrelation 1-3\n")
cat(rbind(acfVLp[2:4], acfPLp[2:4], acfMLp[2:4]),"\n\n")

##############################################################################
# 6) Estimates without any selection
##############################################################################
xVall <- cbind(DxV, Zviol)
bVall <- solve(crossprod(xVall), crossprod(xVall, DyV))
sVall <- cluster_se(xVall, DyV - xVall %*% bVall, solve(crossprod(xVall)), Dstate, ncol(xVall))

xPall <- cbind(DxP, Zprop)
bPall <- solve(crossprod(xPall), crossprod(xPall, DyP))
sPall <- cluster_se(xPall, DyP - xPall %*% bPall, solve(crossprod(xPall)), Dstate, ncol(xPall))

xMall <- cbind(DxM, Zmurd)
bMall <- solve(crossprod(xMall), crossprod(xMall, DyM))
sMall <- cluster_se(xMall, DyM - xMall %*% bMall, solve(crossprod(xMall)), Dstate, ncol(xMall))

cat("No Selection\n")
cat("Violence\n")
cat(c(bVall[1], sVall[1]), "\n")
cat("Property\n")
cat(c(bPall[1], sPall[1]), "\n")
cat("Murder\n")
cat(c(bMall[1], sMall[1]), "\n")