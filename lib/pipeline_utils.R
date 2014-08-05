library(limma)
library(maDB)
library(RColorBrewer)
library('yaml')
library('outliers')

getType<-function(config) {
	if(config$algorithms$median) {
		return("agilent.median");
	} else {
		return("agilent")
	}
}

getUniqueGeneNames<-function(config) {
	if(config$algorithms$median) {
		return(matrix(unique(Slides.norm$genes$ProbeName),ncol=1))
	} else {
		return(matrix(unique(Slides.norm$genes$GeneName),ncol=1))
	}
}

rm_outliers  <-function(dat, name) {
	punk <- outlier(dat)
	outlier_mask = dat==punk
	test=dixon.test(dat)
	if(test$p.value < config$dixon$cutoff) {
		cat("Removed outlier ", punk, " from dataset ", 
				name, "(p=", sprintf('%0.5e',test$p.value), ')\n')
		return(dat[!outlier_mask])
	} else {
		return(dat)
	}
}

validateSlides <- function(slideList) {
  for(slide in slideList) {
    stopifnot(!is.null(slide$name))
    stopifnot(!is.null(slide$file))
    stopifnot(!is.null(slide$IP))
  }
}

# validateSlides <- function(slideList) {
#   for(slide in slideList) {
#     stopifnot(!is.null(slide$name))
#     stopifnot(!is.null(slide$raw_file))
#     stopifnot(!is.null(slide$swap_file))
#   }
# }

writeMetadata <- function(slideList, outfile) {
  validateSlides(slideList)
  f=file(outfile, 'w')
  write("Name\tFileName\tCy5\tCy3",f)
  for(slide in slideList) {
    if(slide$IP == "cy5")
      s <- sprintf("%s\t%s\t%s\t%s", slide$name, slide$file, slide$name, "control")
    else
      s <- sprintf("%s\t%s\t%s\t%s", slide$name, slide$file, "control", slide$name)
    write(s, f)
  }
  unlink(f)
}