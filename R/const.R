const <- function(x, print = TRUE){
  
  const.labels <- c("uncons", "nonneg", "period", "pernon",
                    "smooth", "smonon", "smoper", "smpeno",
                    "orthog", "ortnon", "ortsmo", "orsmpe",
                    "moninc", "monnon", "monsmo", "mosmno",
                    "unimod", "uninon", "uniper", "unpeno", 
                    "unismo", "unsmno", "unsmpe", "unsmpn")
  const.details <- c("unconstrained", 
                     "non-negative", 
                     "periodic", 
                     "periodic and non-negative",
                     "smooth", 
                     "smooth and non-negative",
                     "smooth and periodic",
                     "smooth, periodic, and non-negative",
                     "orthogonal", 
                     "orthogonal and non-negative",
                     "orthogonal and smooth", 
                     "orthogonal, smooth, and periodic",
                     "monotonic increasing", 
                     "monotonic increasing and non-negative",
                     "monotonic increasing and smooth", 
                     "monotonic increasing, smooth, and non-negative",
                     "unimodal",
                     "unimodal and non-negative",
                     "unimodal and periodic",
                     "unimodal, periodic, and non-negative",
                     "unimodal and smooth",
                     "unimodal, smooth, and non-negative",
                     "unimodal, smooth, and periodic",
                     "unimodal, smooth, periodic, and non-negative")
  contab <- data.frame(paste(const.labels, ": "), const.details)
  
  if(missing(x)){
    
    if(print) {
      cat(t(cbind(contab, "\n")), sep = "")
    } else {
      return(data.frame(label = const.labels, details = const.details))
    }
    
  } else {
    
    id <- pmatch(x, const.labels, duplicates.ok = TRUE)
    if(print){
      cat(t(cbind(contab, "\n")[id,]), sep = "")
    } else {
      return(data.frame(label = const.labels[id], details = const.details[id]))
    }
    
  } # end if(missing(x))
  
}