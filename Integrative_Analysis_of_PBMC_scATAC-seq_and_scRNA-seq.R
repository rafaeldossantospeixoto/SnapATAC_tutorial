# Install packages
# install.packages('doSNOW')
# library(devtools)
# install_github("r3fang/SnapATAC")

library(ggplot2)

# Step 1. Barcode selection -----------------------------------------------

library(SnapATAC);

# Load supplemental code to fix errors with the package
# load(supplemental_code.R)

snap.files = c(
  "PBMC_ATAC_RNA/atac_pbmc_5k_nextgem.snap", 
  "PBMC_ATAC_RNA/atac_pbmc_10k_nextgem.snap"
);
sample.names = c(
  "PBMC 5K",
  "PBMC 10K"
);
barcode.files = c(
  "PBMC_ATAC_RNA/atac_pbmc_5k_nextgem_singlecell.csv",
  "PBMC_ATAC_RNA/atac_pbmc_10k_nextgem_singlecell.csv"
);
x.sp.ls = lapply(seq(snap.files), function(i){
  createSnap(
    file=snap.files[i],
    sample=sample.names[i]
  );
})
names(x.sp.ls) = sample.names;
barcode.ls = lapply(seq(snap.files), function(i){
  barcodes = read.csv(
    barcode.files[i], 
    head=TRUE
  );
  # remove NO BAROCDE line
  barcodes = barcodes[2:nrow(barcodes),];
  barcodes$logUMI = log10(barcodes$passed_filters + 1);
  barcodes$promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1);
  barcodes
})
plots = lapply(seq(snap.files), function(i){
  p1 = ggplot(
    barcode.ls[[i]], 
    aes(x=logUMI, y=promoter_ratio)) + 
    geom_point(size=0.3, col="grey") +
    theme_classic()	+
    ggtitle(sample.names[[i]]) +
    ylim(0, 1) + xlim(0, 6) + 
    labs(x = "log10(UMI)", y="promoter ratio")
  p1
})
plots

x.sp.ls

# for both datasets, we identify usable barcodes using [3.5-5] for log10(UMI) and [0.4-0.8] for promoter ratio as cutoff.
cutoff.logUMI.low = c(3.5, 3.5);
cutoff.logUMI.high = c(5, 5);
cutoff.FRIP.low = c(0.4, 0.4);
cutoff.FRIP.high = c(0.8, 0.8);
barcode.ls = lapply(seq(snap.files), function(i){
  barcodes = barcode.ls[[i]];
  idx = which(
    barcodes$logUMI >= cutoff.logUMI.low[i] & 
      barcodes$logUMI <= cutoff.logUMI.high[i] & 
      barcodes$promoter_ratio >= cutoff.FRIP.low[i] &
      barcodes$promoter_ratio <= cutoff.FRIP.high[i]
  );
  barcodes[idx,]
});
x.sp.ls = lapply(seq(snap.files), function(i){
  barcodes = barcode.ls[[i]];
  x.sp = x.sp.ls[[i]];
  barcode.shared = intersect(x.sp@barcode, barcodes$barcode);
  x.sp = x.sp[match(barcode.shared, x.sp@barcode),];
  barcodes = barcodes[match(barcode.shared, barcodes$barcode),];
  x.sp@metaData = barcodes;
  x.sp
})
names(x.sp.ls) = sample.names;
x.sp.ls

# combine two snap object
x.sp = Reduce(snapRbind, x.sp.ls);
x.sp@metaData["sample"] = x.sp@sample;
x.sp

table(x.sp@sample);


# Step 2. Add cell-by-bin matrix ------------------------------------------

# Error with Matrix::rBind function.
x.sp = addBmatToSnap(x.sp, bin.size=5000);



