#' Combine snap objects
#'
#' Takes two snap objects and combines them. The original function found in 
#' github was throwing an error because Matrix::rBind is defunct, therefore I 
#' eddited the code to use base::rbind instead of Matrix::rBind.
#' https://github.com/cole-trapnell-lab/monocle3/issues/509
#'
#' @param obj1 a snap object
#' @param obj2 a snap object
#' @return a combined snap object
#' @export
snapRbind <- function(obj1, obj2){
  if(!is.snap(obj1)){stop(paste("Error @snapRbind: obj1 is not a snap object!", sep=""))};
  if(!is.snap(obj2)){stop(paste("Error @snapRbind: obj2 is not a snap object!", sep=""))};
  
  # barcode from obj1 and obj2
  barcode1 = paste(obj1@file, obj1@barcode, sep=".");
  barcode2 = paste(obj2@file, obj2@barcode, sep=".");	
  
  # check barcode name, if there exists duplicate barcode raise error and exist
  if(length(unique(c(barcode1, barcode2))) < length(barcode1) + length(barcode2)){
    stop("Error: @snapRbind: identifcal barcodes found in obj1 and obj2!")
  }
  rm(barcode1, barcode2)
  gc()
  
  # check meta data
  if(nrow(obj1@metaData) > 0 && nrow(obj2@metaData) > 0){
    metaData = rbind(obj1@metaData, obj2@metaData);		
  }else{
    metaData = data.frame();
  }
  
  # check feature
  feature1 = obj1@feature;
  feature2 = obj2@feature;
  if((length(feature1) == 0) != (length(feature2) == 0)){
    stop("different feature found in obj1 and obj2!")
  }else{
    if(length(feature1) > 0){
      if(FALSE %in% (feature1$name == feature2$name)){
        stop("Error: @snapRbind: different feature found in obj1 and obj2!")
      }
      feature = feature1;					
    }else{
      feature = feature1;								
    }
  }
  gc()
  
  # check peak
  peak1 = obj1@peak;
  peak2 = obj2@peak;
  if((length(peak1) == 0) != (length(peak2) == 0)){
    stop("different peak found in obj1 and obj2!")
  }else{
    if(length(peak1) > 0){
      if(FALSE %in% (peak1$name == peak2$name)){
        stop("Error: @snapRbind: different feature found in obj1 and obj2!")
      }
      peak = peak1;					
    }else{
      peak = peak1;								
    }
  }
  rm(peak1, peak2)
  gc()
  
  # check bmat	
  bmat1 = obj1@bmat;
  bmat2 = obj2@bmat;
  if((length(bmat1) == 0) != (length(bmat2) == 0)){
    stop("bmat has different dimentions in obj1 and obj2!")
  }else{
    bmat = rbind(bmat1, bmat2);
  }
  rm(bmat1, bmat2)
  gc()
  
  # check gmat	
  gmat1 = obj1@gmat;
  gmat2 = obj2@gmat;
  if((length(gmat1) == 0) != (length(gmat2) == 0)){
    stop("gmat has different dimentions in obj1 and obj2!")
  }else{
    gmat = rbind(gmat1, gmat2);
  }
  rm(gmat1, gmat2)
  gc()
  
  # check pmat	
  pmat1 = obj1@pmat;
  pmat2 = obj2@pmat;
  if((length(pmat1) == 0) != (length(pmat2) == 0)){
    stop("pmat has different dimentions in obj1 and obj2!")
  }else{
    pmat = rbind(pmat1, pmat2);
  }
  rm(pmat1, pmat2)
  gc()
  
  
  # check gmat	
  dmat1 = obj1@smat@dmat;
  dmat2 = obj2@smat@dmat;
  
  if((length(dmat1) == 0) != (length(dmat2) == 0)){
    stop("dmat has different dimentions in obj1 and obj2!")
  }else{
    dmat = rbind(dmat1, dmat2);
  }
  rm(dmat1, dmat2)
  gc()
  
  res = newSnap();
  res@feature = feature;
  res@barcode = c(obj1@barcode, obj2@barcode);
  res@file = c(obj1@file, obj2@file);
  res@sample = c(obj1@sample, obj2@sample);
  res@metaData = metaData;
  res@bmat = bmat;
  res@pmat = pmat;
  res@peak = peak;
  res@gmat = gmat;
  res@smat@dmat = dmat;
  res@smat@sdev = obj1@smat@sdev;
  return(res)
}

# Does not solve the Matrix::rBind error.
assignInNamespace(x='snapRbind', snapRbind, ns='SnapATAC')
# assignInNamespace(x='readBins', rbind, ns='Matrix')
