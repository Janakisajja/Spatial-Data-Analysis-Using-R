f <- system.file("WILDPOT.txt", package="rspatial")
f
d <- read.table(f, header=TRUE)

d <- readLines(f, encoding='UTF-8')

dd <- strsplit(d, '\t')

table(sapply(dd, length))


fun <- function(x) {
  r <- rep("", 22)
  r[1:length(x)] <- x
  r
}


ddd <- lapply(dd, fun)

v <- do.call(rbind, ddd)
head(v)


colnames(v) <- v[1,]
v <- v[-1,]

v <- data.frame(v, stringsAsFactors=FALSE)


for (i in c('LongD', 'LongM', 'LongS', 'LatD', 'LatM', 'LatS')) {
  v[, i] <- as.numeric(v[,i])
}
v$lon <- -1 * (v$LongD + v$LongM / 60 + v$LongS / 3600)
v$lat <- v$LatD + v$LatM / 60 + v$LatS / 3600

v$lat[v$LatH == 'S'] <- -1 * v$lat[v$LatH == 'S']
head(v)

library(raster)
library(rspatial)
cn <- sp_data('pt_countries')
proj4string(cn) <- CRS("+proj=longlat +datum=WGS84")
class(cn)

plot(cn, xlim=c(-120, -40), ylim=c(-40,40), axes=TRUE)
points(v$lon, v$lat, cex=.5, col='red')

sp <- SpatialPoints( v[, c('lon', 'lat')],
                     proj4string=CRS("+proj=longlat +datum=WGS84") )
sp <- SpatialPointsDataFrame(sp, v)

table(v$COUNTRY)

v$COUNTRY <- toupper(v$COUNTRY)
table(v$COUNTRY)

sp$COUNTRY <- toupper(sp$COUNTRY)

ov <- over(sp, cn)
colnames(ov) <- 'name'
head(ov)

v <- cbind(v, ov)
table(v$COUNTRY)


v$name[is.na(v$name)] <- ''

v$name[v$name=="UNITED STATES, THE"] <- "UNITED STATES"
v$name[v$name=="BRASIL"] <- "BRAZIL"
i <- which(toupper(v$name) != v$COUNTRY)
i

plot(cn, xlim=c(-120, -40), ylim=c(-40,40), axes=TRUE)
points(sp, cex=.25, pch='+', col='blue')
points(sp[i,], col='red', pch='x', cex=1.5)

spc <- tapply(v$SPECIES, sp$COUNTRY, function(x)length(unique(x)) )
spc <- data.frame(COUNTRY=names(spc), nspp = spc)

cn <- merge(cn, spc, by='COUNTRY')
print(spplot(cn, 'nspp', col.regions=rev(terrain.colors(25))))

tb <- table(v[ c('COUNTRY', 'SPECIES')])
dim(tb)

tb[,2:3]

install.packages('rgdal')

library(rgdal)

projection(cn) <- "+proj=longlat +datum=WGS84"

laea <- CRS("+proj=laea  +lat_0=0 +lon_0=-80")
clb <- spTransform(cn, laea)
pts <- spTransform(sp, laea)
plot(clb, axes=TRUE)
points(pts, col='red', cex=.5)

sp
pts
Uni <- unique(pts@data[["SPECIES"]])
Uni

r <- raster(clb)

res(r) <- 50000

rich <- rasterize(pts, r, 'SPECIES', function(x, ...) length(unique(na.omit(x))))
plot(rich, main = "Species Richness at 50km resolution")
plot(clb, add=TRUE)

obs <- rasterize(pts, r, field='SPECIES', fun=function(x, ...)length((na.omit(x))) )
plot(obs,  main = "Number of observations at 50km resolution")
plot(clb, add=TRUE)

r1 <- raster(clb)

res(r1) <- 100000

rich <- rasterize(pts, r1, 'SPECIES', function(x, ...) length(unique(na.omit(x))))
plot(rich, main = "Species Richness at 100km resolution")
plot(clb, add=TRUE)

obs <- rasterize(pts, r1, field='SPECIES', fun=function(x, ...)length((na.omit(x))) )
plot(obs,  main = "Number of observations at 100km resolution")
plot(clb, add=TRUE)

r2 <- raster(clb)

res(r2) <- 250000

rich <- rasterize(pts, r2, 'SPECIES', function(x, ...) length(unique(na.omit(x))))
plot(rich, main = "Species Richness at 250km resolution")
plot(clb, add=TRUE)

obs <- rasterize(pts, r2, field='SPECIES', fun=function(x, ...)length((na.omit(x))) )
plot(obs,  main = "Number of Observations at 250km resolution")
plot(clb, add=TRUE)

r3 <- raster(clb)

res(r3) <- 500000

rich <- rasterize(pts, r3, 'SPECIES', function(x, ...) length(unique(na.omit(x))))
plot(rich, main = "Species Richness at 500km resolution")
plot(clb, add=TRUE)

obs <- rasterize(pts, r3, field='SPECIES', fun=function(x, ...)length((na.omit(x))) )
plot(obs,  main = "Number of Observations at 500km resolution")
plot(clb, add=TRUE)

#Shannon Diversity
vv <- as.vector(table(v$SPECIES))
p <- vv / sum(vv)
p
install.packages("SciViews")
library(SciViews)
H = -sum(p * ln(p))
H
res(r) <- 200000

plot(H)
plot(clb, add=TRUE)




