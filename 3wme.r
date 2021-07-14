#Berehnung des Durchmesserzuwachses für Fichte in einem Jahr:

x <- read.table(header=TRUE, text="
   d    h   nrep spGrp
45.3 34.0  24.82     1
35.8 31.5  39.74     1
50.3 28.8  20.13     1
30.2 24.3  55.84     1
16.8 11.0 180.45     1
13.7  7.5 271.35     1
 7.7  4.7 470.87     1
47.0 32.2  23.06     2
50.3 31.7  20.13     3
39.4 33.1  32.81     4
50.4 30.7  20.05     5
57.5 30.1  15.40     6")

#spGrp: 1=PiAb, 2=AbAl, 3=FaSy, 4=LaDe, 5=OTHER, 6=PiSy

c0 <- 10 #Bonitaet

x$d0 <- ifelse(x$h > 1.3, x$d + 1.3, x$h) #Stammdurchmesser am Boden

#Standflaeche
funA <- function(y, r=c(1, 1, 0.92, 0.72, 1.6, 0.34), e=1.6) {
  10000 * (1.3 + y$d0)^e*r[y$spGrp] / sum((1.3 + y$d0)^e*r[y$spGrp] * y$nrep)
}
x$stfl <- funA(x)
x$stfl <- pmin(x$h^2 * pi, x$stfl) #Maximale Standflaeche

funNrep <- function(nrep) { #Nrep der Nachbarn bzw. Konkurrenten
  z <- matrix(nrep, length(nrep), length(nrep))
  diag(z) <- pmax(0, diag(z) - 1)
  z}

#Konkurrenz durch Hoenhenunterschied
z <- funNrep(x$nrep) * outer(x$h, x$h*1.07, "-")
z[z<0] <- 0
x$hhSp <- do.call(cbind, lapply(1:6, function(k) {
  i <- which(x$spGrp==k)
  if(length(i)>0) colSums(z[i,, drop = FALSE]) else 0}))

#Konkurrenz durch Kreisflaechenunterschied
z <- funNrep(x$nrep) * outer(x$d0^2, (0.63*x$d0)^2, "-") * pi / 40000
z[z<0] <- 0
x$cbaSp <- do.call(cbind, lapply(1:6, function(k) {
  i <- which(x$spGrp==k)
  if(length(i)>0) colSums(z[i,, drop = FALSE]) else 0}))

#Schaetzwert fuer Zuwachspottential
zM <- exp(2.77635 - 2.02020 * log(1.3 + x$d0) + 0.72712 * log(c0) + 2.76468 * log(x$h) - 2.59996 * x$h/(1.3 + x$d0) - 0.19262 * x$h/c0)

#Schaetzwert fuer Zuwachs je ha
x$ghaMax <- 160.88389592 * (1 - 1/(1+x$h*0.04631407))
x$gha <- (1.3+x$d0)^2*pi/4 / x$stfl
f <- function(b, c0, h0, d0, c1=-8.638e+00, c2=-6.676e+00, c3=-1.309e+01, c4=-1.259e+01, c5=-1.076e+00, c6=-3.474e+00, c7=8.790e+00, c8=5.506e+00, c9=1.921e+01, c10=1.631e+01) {
  .. <- 0.5 + 1/(1 + exp(c1 + c5*c0 + c6*h0/c0 + c7*log(d0) + c8*log(h0/d0) ))
  . <- (0.5 + 0.5/(1+exp(c2 + c3*log(c0) + c4*h0/c0 + c9*log(d0) + c10*log(h0/d0) )))
  1 - exp(-.. * log(.^2/(. - b)^2))
}
z <- zM / 0.95 * f(pmin(1, x$gha / x$ghaMax), c0, x$h, x$d0) *
  1/(1 + (x$hhSp %*% c(2.405e-04, 1.971e-04, 5.396e-04, 1.307e-04, 1.526e-04, 2.251e-04))^9.053e-01) *
  1/(1 + (x$cbaSp %*% c(8.977e-03, 8.815e-03, 3.440e-02, 1.332e-02, 4.094e-03, 8.550e-03))^1.156)
z[x$spGrp != 1] <- NA
z[] <- pmin(c0 * 6.255492, z) #Maximaler Zuwachs
z <- sqrt(x$d0^2 + (z * x$stfl / x$h * 4 / pi)) - x$d0 #Schaetzwert umrechnen in Durchmesserzuwachs
x <- cbind(x, id=z)

x[,c("d","h","spGrp","id")]
#     d    h spGrp        id
#1  45.3 34.0     1 0.3923514
#2  35.8 31.5     1 0.3488206
#3  50.3 28.8     1 0.4045024
#4  30.2 24.3     1 0.2891614
#5  16.8 11.0     1 0.2269061
#6  13.7  7.5     1 0.1849388
#7   7.7  4.7     1 0.1441662
#8  47.0 32.2     2        NA
#9  50.3 31.7     3        NA
#10 39.4 33.1     4        NA
#11 50.4 30.7     5        NA
#12 57.5 30.1     6        NA
