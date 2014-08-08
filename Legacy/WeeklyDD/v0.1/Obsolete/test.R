dyn.load(file.path(paste("/home/mkienzle/mystuff/Programs/C++/DelayDifference/WeeklyDD/", "test", .Platform$dynlib.ext, sep="")))

test <- function(x)
  .C("test",
     as.double(x), # a vector of doubles representing the variable on which selectivity depend on (length or weight)
     result = double(length(x)))$result
