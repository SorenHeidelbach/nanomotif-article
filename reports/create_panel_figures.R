
######## 
pacman::p_load(
  "patchwork",
  "grid",
  "gridExtra",
  "png",
  "imager"
)


# Load the images
img1 <- readPNG("figures/benchmark_heatmap_bin_consensus.png")
img2 <- readPNG("figures/benchmark_heatmap_both.png")
img4 <- readPNG("figures/zymo_and_mono_heatmap.png")

# Resize images if necessary, for example, img1 is the largest, keep as is,
# resize others proportionally or to specific dimensions
# img2 <- image_scale(img2, "x400") # Example resizing command

# Create a plot window
par(mfrow=c(2,2), mar=c(1,1,1,1))

plot.new()
# Plot each image
{
grid.raster(img1)
grid.raster(img2)
grid.raster(img4)
}
