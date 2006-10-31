library(adimpro)

cat("Testing R-package adimpro\n\n")



cat("1. Reading files\n")

cat("1.1. Reading images in RAW format ... ")
#file.raw <- system.file("img/CRW_2844.CRW",package="adimpro")
#img.raw <- read.raw(file.raw)
cat("no file ... OK\n")

cat("1.2. Reading images in PPM format ... ")
file.ppm <- system.file("img/wias.ppm",package="adimpro")
img.ppm <- read.image(file.ppm)
cat("no file OK\n")




cat("2. Transforming color spaces ...")
img.yuv <- rgb2yuv(img.ppm)
img.yiq <- rgb2yiq(img.ppm)
img.hsi <- rgb2hsi(img.ppm)
img.xyz <- rgb2xyz(img.ppm)
img.grey <- rgb2grey(img.ppm)
X11(height=3,width=10)
par(mfrow=c(1,5),mar=c(2,2,1,0.1))
show.image(img.yuv,main="YUV")
show.image(img.yiq,main="YIQ")
show.image(img.hsi,main="HSI")
show.image(img.xyz,main="XYZ")
show.image(img.grey,main="GREYSCALE")
cat("OK\n")
dev.off()




cat("3. Showing images ...")
X11(height=3,width=6)
#show.image(img.raw)
show.image(img.ppm)
cat("OK\n")
dev.off()




cat("4. Writing images ...")
write.image(img.ppm,file="test_8bit.tif")
#write.image(img.raw,file="test_16bit.tif")
write.image(img.grey,file="test_grey.tif")
cat("OK\n")




cat("5. Adaptive smoothing ...")
#imghat <- awsimage(img.raw,graph=TRUE)
#imghat <- awspimage(img.raw,graph=TRUE)
cat("OK\n")





cat("6. Info on images ...")
#plot(imghat)
#plot(img.raw)
#summary(imghat)
#summary(img.raw)
cat("OK\n")



cat("7. Egde detection ...")
img.edge <- edges(img.ppm)
show.image(img.edge)
cat("OK\n")




cat("8. Histogram equalization ...")
img.eq <- histogram.equalization(img.ppm)
show.image(img.eq)
cat("OK\n")




cat("9. Rotate image ...")
img.rot <- rotate.image(img.ppm)
show.image(img.rot)
cat("OK\n")




cat("10. Shrink image ...")
img.shrink <- shrink.image(img.raw)
show.image(img.shrink)
cat("OK\n")

