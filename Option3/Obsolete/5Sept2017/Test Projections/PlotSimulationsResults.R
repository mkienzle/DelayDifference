# weekly.data <- read.csv("Results/Simulation/Projections_fortest.txt")
# par(mfrow=c(4,1), mai = c(0,0,0,0))
# with(weekly.data, plot(timestep, biomass)); abline(v = 52 * seq(1,500));
# with(weekly.data, plot(timestep, recruitment)); abline(v = 52 * seq(1,500));
# with(weekly.data, plot(timestep, catch)); abline(v = 52 * seq(1,500));
# with(weekly.data, plot(timestep, effort)); abline(v = 52 * seq(1,500));

yearly.data <- read.csv("Results/Simulation/Projections_fortest_byYear.txt")

yearly.data <- within(yearly.data,
	    Effort <- as.factor(Effort))
	    
library(ggplot2)
p1 <- ggplot(yearly.data, aes(x=Effort, y=Catch)) + geom_boxplot() + facet_wrap( ~ Recruitment) +
   scale_x_discrete( breaks = seq(0,5e4, 1e4), labels = seq(0,5e4, 1e4))
print(p1)
ggsave("Results/Graphics/Test-CatchVsEffortForVariousLevelOfConstantRecruitment.pdf")
