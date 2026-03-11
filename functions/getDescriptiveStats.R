getDescriptiveStats <- function(dat, metric){
  data.frame(num_O = sum(dat$obs == 1), 
             num_U = sum(dat$obs == 0),
             sum_O = sum(dat[[metric]][dat$obs == 1]), 
             sum_U = sum(dat[[metric]][dat$obs == 0]),
             mean_O = mean(dat[[metric]][dat$obs == 1]), 
             mean_U = mean(dat[[metric]][dat$obs == 0]),
             sd_O = sd(dat[[metric]][dat$obs == 1]),
             sd_U = sd(dat[[metric]][dat$obs == 0]),
             freq_O = mean(dat[[metric]][dat$obs == 1]>0),
             freq_U = mean(dat[[metric]][dat$obs == 0]>0) 
  )
}
