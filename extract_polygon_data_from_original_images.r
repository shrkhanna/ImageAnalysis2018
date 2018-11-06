#########################################
##   EXTRACT VECTOR DATA FROM IMAGES   ##
##      example with AVIRIS data       ##
##     Developed by Shruti Khanna      ##
##  Date last updated: Apr 10, 2018    ##
#########################################


#########  CHANGE  #########
# working directory with the R program file
dir_progr = "X:/delta_sav/code/R/extract_data/"
######## END CHANGE ########

# set working directory
setwd(dir_progr)

# install the necessary packages and load them
# install.packages(c("raster","sp","rgdal","randomForest","asbio"))
# install.packages("rgeos")
require("raster")
require("sp")
require("rgdal")
require("rgeos")


# name all necessary folders and files

#########  CHANGE  #########

# image directory
dir_image = "X:/scratch/ang/Delta2017/masked_coreg_orig/"
# training/test data directory without the last slash front slash (very important)
dir_shape = "X:/delta_sav/vector/field_data/2017/data3"
# output directory for both csv files
dir_out   = "X:/delta_sav/raster/analysis/indx_waterline/201711/"
# suffix of folders and images to be processed
fldpre = "Delta"
imgpat = "rot_mskd.img"
# name of training shapefile without the .shp extension
name_train  = "2017_water_sav_final_forSVI"
# name of training csv files
name_outcsv = paste(dir_out, "R_WatSAV_spectra_201711.csv", sep="")


######## END CHANGE ########


# get list of folders
fld_list <- list.files(dir_image, pattern = imgpat,full.names = TRUE)
no_folders <- length(fld_list)


# img_list <- subset(img_list, !grepl(".aux.xml", img_list))
#img_list = fld_list


#########  CHANGE  #########
# start at this file
stfile = 1
# end at the last file
enfile = no_folders
# total number of bands in the image
bnumbr = 425
# mask value
mskval = -9999
# extract starting from this band
stband = 20
# extract ending at this band
enband = 150
######## END CHANGE ########

# IMPORTANT information about your training and test vectors
# make sure your training and test data have exactly the same fields
# create a new field called ORIG_FID in your training and test shapefiles that is the same as your FID
# the FID number will change through this script but the ORIG_FID should stay the same

# read the field data shapefile: readOGR(dsn = destination folder, layer = filename)
vecname <- readOGR(dsn = dir_shape, layer = name_train)

# get names of field data columns
fields <- colnames(vecname@data)

# create basic band names: "band #"
bnames = paste('band', as.character(seq(1, bnumbr, by=1)), sep="")
ttl_cols = 1 + bnumbr + length(fields)

# initialize the flag for capturing first intersecting polygons
first = 0

# RESOLVE ERROR FROM I = 26
#############################

for (i in stfile:enfile) {
  
  # create the correct image name with path
  #img_list = list.files(paste(dir_image, fld_list[i], sep=''), pattern = imgpat)
  #img_list = subset(img_list, !grepl(".hdr", img_list))
  #img_list = subset(img_list, !grepl("coreg", img_list))
  
  # x = length(img_list)
  # if (x < 1) { 
  #   next
  # }
  #img_name = paste(dir_image, fld_list[i], '/', img_list, sep="")
  img_name = fld_list[i]
  # read the input image (give entire path)
  input_image = brick(img_name)
  #input_image@crs
  # get subset of points that intersect with the image
  vint = intersect(vecname, as(extent(input_image), 'SpatialPolygons'))
  # check projection between vector and raster: vint@proj4string & input_image@crs
  
  if (!is.null(vint)) {
    numpoly = nrow(vint)
    if (first == 0) {
      # first file found where vector data intersected with image
      # turn flag "ON" so that master array can be initialized
      first = i
    }
  } else {
    numpoly = 0
  }

  if (!is.null(vint)) {
    
    # Overlay vector data on image and extract bands as a data frame
    exvec = extract(input_image, vint, cellnumbers=TRUE)
    # Make a vector of no. of pixels extracted per polygon
    freq = t(as.data.frame(lapply(exvec, length)))
    # divide by number of columns to get number of rows
    freq = freq/(bnumbr+1)
    # Total number of pixels extracted from the image
    total = sum(freq)
    
    # convert large list into a single list
    subvec = NULL
    for (j in 1:numpoly) {
      subvec = rbind(subvec, as.data.frame(exvec[[j]]))
    }    
    names(subvec) = c("FID", bnames)
    
    # add columns to vector data containing number of pixels in each polygon
    vint_freq = as.data.frame(cbind(vint@data, freq))
    # add name of column containing the pixel count
    names(vint_freq) = c(fields, "freq")
    # make this new vector data frame have the same number of rows as the extracted data (duplicate rows for no. of pixels)
    vector = vint_freq[rep(row.names(vint_freq), vint_freq$freq),]
    
    # bind extracted data to the vector data which are now equal length
    subvecb = cbind(subvec, vector)
    # give correct column headings to the composite file
    names(subvecb) = c("FID", bnames, fields, "freq")
    # remove all masked rows from result data frame
    subvecb = base::subset(subvecb, !((subvecb$band20 == mskval) & (subvecb$band30 == mskval) & (subvecb$band40 == mskval) & (subvecb$band50 == mskval)))
    # select only the bands we want: 1, 33 to 100, 434 to ttl_cols
    subvcol = cbind(subvecb[,1], subvecb[,stband:enband], subvecb[,bnumbr:ttl_cols])

    if (i == first) {

      # initialize the master data frames for test and training
      mastrvec <- subvcol

    } else {
  
      if (!is.null(vint)) {
    
        # remove all fids in following file which have been extracted before
        subvcol = subvcol[!(subvcol$ORIG_FID %in% mastrvec$ORIG_FID),]
        # append new extracted data to existing data (previously extracted)
        mastrvec <- rbind(mastrvec, subvcol)
      }
    }
    
    print(i)

  } # end if no polygons condition
}  # end i for loop

# write table to csv
write.csv(mastrvec, file=name_outcsv)


