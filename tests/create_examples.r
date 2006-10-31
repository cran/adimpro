library(adimpro)

img <- read.image(system.file("img/wias.tif",package="adimpro"))
img.noise <- img
img.noise$img <- img$img + array(rnorm(30000,0,5000),c(100,100,3))
img.noise$img[img.noise$img>65535]<-65535
imghat<-awsimage(img.noise,hmax=15,graph=TRUE)
show.image(img.noise)
write.image(img.noise,file="wias-noise.tif")

a <- array(0,c(100,100,3))
a[30:50,30:50,2] <- 65535
attr(a,"type") <- "rgb"
attr(a,"gamma") <- TRUE
attr(a,"wb") <- "UNKNOWN"
show.image(a)
write.image(a,file="green_box.tif")

a<-array(0,c(100,100))
a[11:100,] <- runif(9000,32768,65535)
attr(a,"type") <- "greyscale"
attr(a,"gamma") <- TRUE
attr(a,"wb") <- "UNKNOWN"
show.image(a)
write.image(a,file="unbalanced_histo.tif")


a <- array(0,c(100,100,3))
a[50,,1] <- 65535
attr(a,"type") <- "rgb"
attr(a,"gamma") <- TRUE
attr(a,"wb") <- "UNKNOWN"
show.image(a)
write.image(a,file="red_line.tif")
