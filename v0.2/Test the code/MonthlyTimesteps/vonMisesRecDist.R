vonMisesRecDist <- function(a,b){


library(circular)
a <- circular(a)
b <- circular(b)

# Interval of time you are interested in
intervals <- rep(NA, 12)
boundaries <- circular(seq(-pi, pi, length =13))
boundaries[13] <- boundaries[13] - 1e-3

  # Calculate von Mises probability associated with each interval

for(i in 1:12){
      # Calculate probability in each interval using the cumulative distribution function

      if(boundaries[i+1] <= (a - pi) ){
	intervals[i] = pvonmises(boundaries[i+1] + 2 * pi, a, b) - pvonmises(boundaries[i] + 2 * pi, a, b);
      }

      if(boundaries[i] <= (a - pi) && (a - pi) < boundaries[i+1]){
	intervals[i] = pvonmises(a + pi, a, b) - pvonmises(boundaries[i] + 2 * pi, a, b)
	  + pvonmises(boundaries[i+1], a, b) - pvonmises(a - pi, a, b);
}

      if(boundaries[i] > (a - pi) && (a + pi) > boundaries[i+1]){
	intervals[i] = pvonmises(boundaries[i+1], a, b) - pvonmises(boundaries[i], a, b);
}

      if(boundaries[i] < (a + pi) && (a + pi) <= boundaries[i+1]){
	intervals[i] = pvonmises(a + pi, a, b) - pvonmises(boundaries[i], a, b)
	  + pvonmises( boundaries[i+1] - 2 * pi, a,b ) - pvonmises( a - pi, a, b);
}

      if(boundaries[i] >= (a + pi)){
	intervals[i] = pvonmises(boundaries[i+1] - 2 * pi, a, b) - pvonmises(boundaries[i] - 2 * pi, a, b);
}

      # sometimes small negative probability appears 
      if(intervals[i] < 0.0)
	# if they are small, replace them by zero
	if(abs(intervals[i]) < 1e-5){ intervals[i] = 0.0;}

      	}
 
    return(intervals)

}
