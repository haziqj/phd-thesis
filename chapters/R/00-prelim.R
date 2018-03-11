## ---- prelim ----

# Packages
library(iprior)
library(iprobit)
library(ipriorBVS)
library(knitr)
library(devtools)
library(R2MLwiN)
library(lme4)
library(jmcm)
library(caret)
library(tidyverse)
library(directlabels)
library(reshape2)
library(kableExtra)
library(cowplot)
library(ElemStatLearn)
library(spatstat)
library(stringr)
library(caret)
library(kernlab)

# ggplot2 theme
ggplot2::theme_set(theme_bw())

# knitr settings
chapter.no <- substr(current_input(), start = 0, stop = 2)
knitr::opts_chunk$set(fig.align = "center", prompt = TRUE, myspac = TRUE,
                      fig.path = paste0("figure/", chapter.no, "-"))
knitr::knit_theme$set("bclear")
options(prompt = "R> ", continue = "+  ", width = 70,
        knitr.table.format = "latex")
knit_hooks$set(myspac = function(before, options, envir) {
  if (before) return("\\singlespacing")
})

# BibLaTeX with Biber backend
system(paste("biber", sub("\\.Rnw$", "", current_input())))

# Cut str
options(str = strOptions(strict.width = "cut"), width = 70)

# Move .tex file to correct chapter
move_tex_to_chapter <- function() {
  this.file <- sub("\\.Rnw$", "\\.tex", current_input())
  chapter.no <- substr(this.file, start = 0, stop = 2)
  file.copy(file.path(this.file), file.path(paste0("../", chapter.no), this.file),
            overwrite = TRUE)
}

# Move figures
move_fig_to_chapter <- function() {
  files <- list.files("figure/")
  files <- files[grep(chapter.no, files)]
  figure.path <- paste0("../", chapter.no, "/figure")
  if (!file.exists(figure.path)) {
    dir.create(file.path(figure.path))
  }
  file.copy(file.path("figure", files), file.path(figure.path, files),
            overwrite = TRUE)
  file.copy(file.path("figure", files), file.path("../../figure", files),
            overwrite = TRUE)
}

# Combine both
move_files_to_chapter <- function() {
  move_fig_to_chapter()
  move_tex_to_chapter()
}

# Delete auxiliary files
delete_files <- function() {
  file.types <- c("*.bbl", "*.bcf", "*.glo", "*.glo-abr", "*.idx", "*.ilg",
                  "*.ind", "*.ist", "*.log", "*.run.xml", "*.slo", "*.pdf",
                  "*.synctex.gz", "*concordance.tex", "*.mw")
  for (i in seq_along(file.types)) {
    file.remove(dir(pattern = file.types[i], full.names = TRUE))
  }
}
# delete_files()

# # Function to determine even numbers
# isEven <- function(x) x %% 2 == 0

