## Prepare_ESCMID_figures.R  -  make submission-sized and poster-sized copies of the figures
## Run in RStudio on the laptop. Needs the magick package: install.packages("magick")
##
## Puts everything in:  F:/Ub/Paper_one_revision/Draft_manis/new_abstract_sub_sweden
##   abstract/  web-sized PNGs for the ESCMID upload (longest side 1000 px, well under 3 MB)
##   poster/    full-resolution copies for printing
## ---------------------------------------------------------------------------------------

SRC <- c(getOption("cholera.proj", "."),                       # where Fig4 was written
         file.path(PROJ))           # add any other folder holding the figures
OUT <- file.path(FIG_DIR, "conference")

## figure file name (without extension) -> the name to save it under
FIGS <- c("Fig4_occurrence_size_maps_v2"        = "ESCMID_Fig1_occurrence_size",
          "Fig5_profiles_communities_classes"   = "ESCMID_Fig2_profiles_communities",
          "Fig6_central_zones"                  = "ESCMID_Fig3_central_zones")

MAX_PX <- 1000        # ESCMID resizes anything larger, so 1000 is the useful maximum
MAX_MB <- 3           # ESCMID file size limit

if (!requireNamespace("magick", quietly = TRUE)) stop("Please run: install.packages(\"magick\")")
library(magick)
dir.create(file.path(OUT, "abstract"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "poster"),   recursive = TRUE, showWarnings = FALSE)

find_fig <- function(stem) {
  hits <- unlist(lapply(SRC, function(d) list.files(d, pattern = paste0("^", stem, "\\.(png|tiff|tif)$"),
                                                    full.names = TRUE, ignore.case = TRUE)))
  if (!length(hits)) return(NA_character_)
  hits[order(!grepl("\\.png$", hits, ignore.case = TRUE))][1]   # prefer the PNG
}

for (i in seq_along(FIGS)) {
  stem <- names(FIGS)[i]; new <- FIGS[[i]]
  f <- find_fig(stem)
  if (is.na(f)) { message("NOT FOUND: ", stem, " - copy it into one of the SRC folders"); next }
  im <- image_read(f)
  inf <- image_info(im)

  ## poster copy: the original, unchanged
  file.copy(f, file.path(OUT, "poster", basename(f)), overwrite = TRUE)

  ## abstract copy: longest side MAX_PX, shrunk further if the file is still too large
  small <- image_resize(im, paste0(MAX_PX, "x", MAX_PX, ">"))
  out_png <- file.path(OUT, "abstract", paste0(new, ".png"))
  image_write(small, out_png, format = "png")
  if (file.size(out_png) / 1e6 > MAX_MB)
    image_write(image_resize(im, "800x800>"), out_png, format = "png")

  s <- image_info(image_read(out_png))
  message(sprintf("%-34s source %5.0fx%-5.0f ->  abstract %4.0fx%-4.0f  %.2f MB",
                  new, inf$width, inf$height, s$width, s$height, file.size(out_png) / 1e6))
}

cat("\nFigure titles for the ESCMID form (each well under 500 characters):\n\n")
cat("Figure 1: Zone-level probability of an outbreak per month and expected cases per outbreak month, Ethiopia, 2015-2024, estimated by a Bayesian spatiotemporal hurdle model.\n\n")
cat("Figure 2: Latent profiles, spatial risk communities and classes of cholera activity in Ethiopia, 2015-2024, with stable and uncertain zones distinguished by shading.\n\n")
cat("Figure 3: Eigenvector centrality of zones within the six spatial risk communities, Ethiopia, 2015-2024.\n\n")
cat("Saved in:", OUT, "\n")
