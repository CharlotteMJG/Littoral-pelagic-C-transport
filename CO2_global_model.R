#### Global CO2 transport for pH<7.2 

# Ensure working directory is set to the location of this script

library(ggplot2)
library(gridExtra)
library(grid)
library(deSolve)
source("Support_functions.R")


d <- read.csv("global.vegetation.cover.csv", h = T)
d<-d[c("size.bin","depth_max", "veg.cover", "global_surface", "mean_area", "r_max","kCO2", "kCH4", "floating", "emergent")]

# global_surface and size_bin in km2, r_max and depth_max in m
# global_surface from Holgerson et al 2016 and Verpoorter et al 2014
# veg.cover in % from Zhang et al. 2017 and Vaderboncoeur et al. 2011
# submerged, emergent, floating: % presence of each vegetation type
# kCO2 and kCH4 mean gas exchange coefficient for CO2 and CH4 respectively (m d-1)
# kCO2 given in Raymond et al. 2013 and kCH4 calculated from kCO2 according to Jaehne et al 1987

d$size.bin <- factor(d$size.bin, levels = c("inf 0.1", "0.1-1", ".1-10", "10-100", "sup 100"))
d <- d[order(d$size.bin, decreasing = F),]
d$global_surface <- as.numeric(d$global_surface)

#For a conservative estimate of vegetation cover:
#for the size bins >1 and <100, we divide the vegetation cover by two
#for the size bin >100, we divide the cover by four
d$veg.cover[d$size.bin==".1-10"]<-d$veg.cover[d$size.bin==".1-10"]/2
d$veg.cover[d$size.bin=="10-100"]<-d$veg.cover[d$size.bin=="10-100"]/2
d$veg.cover[d$size.bin=="sup 100"]<-d$veg.cover[d$size.bin=="sup 100"]/4

### Converting the littoral CO2 concentration from ppm to uM at 20 °C:
Temp <- 20
KHCO2 <- exp(-58.0931 + (90.5069 * (100 / (Temp + 273.15))) +
               (22.294 * log((Temp + 273.15) / 100))) # (mol/L/atm), Weiss 1974
CO2_littoral_ppm<-3552# value used for the maximum estimate. 1681 is used for the mimimum estimate
CO2_littoral <- CO2_littoral_ppm * KHCO2 #uM

##### Transport model: assume dC/dt to be 0 (i.e. steady state) -----
# Model to describe CO2 concentration in a circular lake, based on
# Peeters et al. (2019, doi:10.1038/s41598-018-36530-w)

# We here consider the CO2 and CH4 that is derived from the littoral zone. We assume
# that most of it is lost by evasion to the atmosphere. Equation describes 2nd
# order derivative, see example from Soetaert et al., 2010.

# Function valid for lakes with pH < 7.2, CO2 evasion will lead to CO2 decrease.
C_flux <- function(r, state, parameters){
  with(as.list(c(state, parameters)), {
    Hd = calc_hd(Vd)
   	depth = H_r_johansson(r, r_max = r_max, max_depth = max_depth, max_surf_layer = max_surf_layer, Hd = Hd)
    dHdr = dHdr_johansson(r, r_max = r_max, max_depth = max_depth, max_surf_layer = max_surf_layer, Hd = Hd)
    Fatm = v_gas * (state[1] - C_eq)
	r_s = sqrt((1 - veg_emergent / 100) * r_max^2) 
  if(r < r_s){
    d2C = (r * Fatm / Kh - dHdr * r * state[2] - depth * state[2]) / (depth * r + small_nr)    
  }else{
    # Monod equation to force dC towards 0 (factor 10 to ensure 0 is approached faster)
    d2C = -10 * state[2] / (k_monod + state[2])
  }
   list(c(state[2], d2C), c(d2C = d2C, depth = depth, dHdr = dHdr))
  })
}

C_littoral <- matrix(data = CO2_littoral, nrow = nrow(d)) # Littoral CO2 conc, uM 
Ccenter <- matrix(data = NA, nrow = nrow(d)) # CO2 conc in the center of the lake, uM
Cintegrated <- matrix(data = NA, nrow = nrow(d)) # CO2 average conc, uM
Fatm_integrated <- matrix(data = NA, nrow = nrow(d)) # Lake average CO2 atmospheric flux, that is derived from the littoral zone, mol/m2/s
r_step <- c(0.5, 2, 5, 12, 50)
d <- cbind(d, Ccenter, Cintegrated, Fatm_integrated, C_littoral, r_step)
p <- list()

# Loop over the size bins, running the transport model for each bin
for(i in 1:5){
	max_depth <- d$depth_max[i]
	r_max <- d$r_max[i] # m
	veg <- d$veg.cover[i] # % vegetation cover
	veg_emergent <- d$veg.cover[i] * (d$emergent[i] + d$floating[i]) / 100 # % emergent and floating veg cover
	Kh <- 3.2 * 10^-4 * r_max^(1.10) # m2/s, based on equation S4 from Peeters, with L = r_max
	v_gas <- d$kCO2[i] / (3600 * 24) # Temperature adjusted k for CO2 in m/d in Raymond 2013 (Table S4), v_gas is in m/s
	C_eq <- 400 * KHCO2 # uM, using CO2=400ppm for atmospheric concentration
	
	parms <- c(r_max = r_max,             # Radius of the lake (m)
	           max_depth = max_depth,     # Maximum depth of the lake (m)     
	           max_surf_layer = 6,          # Maximum depth of the littoral zone (m) (or surface mixed layer)¨
	           v_gas = v_gas,   			    # m/s
	           C_eq = C_eq,               # for CO2 uM
	           Kh = Kh,               	  # m2/s, based on equation S4
	           Vd = 0.8,                  # Volume development parameter, Johansson et al. 2007
	           small_nr = 1E-10,          # To avoid division by 0
	           veg = veg,
	           veg_emergent = veg_emergent,
	           k_monod = 0.0001)
	
	# R-values to generate output for - ensure that r_max is included
	r_values <- seq(from = 0, to = parms[["r_max"]], by = d$r_step[i]) # Values of the radius
	if(!(parms[["r_max"]] %in% r_values)){
	  r_values <- c(r_values, parms[["r_max"]])
	}
	
	# Initial condition in centre of lake - define automatically using optim
	get_littoral_c_abs_error <- function(Cc, obj = CO2_littoral){
	  init = c(C = Cc,
	           dCdr = 0)
	  out = ode(y = init, times = r_values, func = C_flux, parms = parms)
	  abs(out[nrow(out), "C"][[1]] - obj)
	}
	
	C_init <- optim(c(Cc = 10), get_littoral_c_abs_error, method = "Brent",
	                lower = 0.1, upper = CO2_littoral)$par
	
	init <- c(C = C_init, # state variable initial conditions, CO2 concentration in the middle of the lake
	          dCdr = 0) # Initial condition of derivative
	
	# Run the model and plot output
	out <- ode(y = init, times = r_values, func = C_flux, parms = parms)
	df_out <- data.frame(out)
	
	# Wrap in local() is required to create multiple plots in a loop, as r_max is changed.
	p[[i]] <- local({
	  r_max = r_max
	  ggplot(df_out, aes(x=time, y=C/C[time == r_max]))+
	    geom_point()+
	    theme_bw() +
	    theme(panel.border = element_blank(),
	          panel.grid.major = element_blank(),
	          panel.grid.minor = element_blank(),
	          axis.line = element_line(colour = "black"),
	          legend.text=element_text(size=10),
	          legend.position="none",
	          axis.text=element_text(size=9),
	          axis.title=element_text(size=10),
	          plot.title = element_text(size=10)) +
	    scale_x_reverse() +
	    annotate(geom = "text", x = r_max, y = 1.1,
	             label = paste("r[max] == ",round(r_max, 0)),
	             hjust = 0, vjust = 1, parse = TRUE) +
	    ylab(NULL) +
	    xlab(NULL) 
	})
	d$C_littoral[i] <- df_out$C[df_out$time == r_max]  
	d$Ccenter[i] <- df_out$C[1]
	d$Cintegrated[i] <- mean(df_out$C)
	r_s_tot <- sqrt((1 - veg / 100) * r_max^2) # Calculating r_s_tot to identify the boundary of the littoral habitat (including all vegetation types)

	df_out_pelagic <- subset(df_out, time < r_s_tot) # Excluding vegetated littoral zone
	d$Fatm_pelagic[i] <- v_gas * (mean(df_out_pelagic$C) - C_eq) * 10^3 * 10^-6 # uM*m/s *10^3*10^-6=> mol/m3*m/s =>mol/m2/s
	d$Fatm_integrated[i] <- v_gas * (mean(df_out$C) - C_eq) * 10^3 * 10^-6 # =>mol/m2/s 
}

pdf("Proportion CO2 decrease.pdf") 
grid.arrange(arrangeGrob(p[[1]],p[[2]], p[[3]], p[[4]], p[[5]], ncol = 2,
                         heights = c(2, 2, 2),
                         bottom = textGrob("r (m)",
                                         gp = gpar(col = "black", fontsize = 14)),
                         left = textGrob(expression(paste("C/C"["littoral"])),
                                         gp = gpar(col = "black", fontsize = 14), rot = 90)))
dev.off() 

ratio<- 2#gas concentration in vegetated areas/gas concentration in unvegetated areas
time <- 60*60*24*365 # number of s in 1 yr
M <- 12 # molar mass C in g/mol
surface <- 10^6 # number of m2 in 1 km2
d$Ftot_pelagic <- d$global_surface * (1 - d$veg.cover/100) *
	d$Fatm_pelagic * M * surface * time /ratio/ 10^15 # PgC/yr, CO2 atmospheric flux in the pelagic zone originating from the littoral zone,  *10^6 because global surface is in km2
d$Ftot_integrated <- d$global_surface *
	d$Fatm_integrated * M * surface * time /ratio/ 10^15 # PgC/yr, CO2 atmospheric flux of the lake originating from the littoral zone, taking whole lake integrated CO2 concentrations
d$pCO2Center <- d$Ccenter/KHCO2 # Concentration in uatm 
d$pCO2integrated <- d$Cintegrated / KHCO2 # Concentration in uatm
d$propCcenter_littoral <- d$Ccenter / d$C_littoral# Proportion of center CO2 concentration derived from CO2 in the littoral zone
d

summary_table <- d[, c("size.bin", "pCO2Center", "propCcenter_littoral","Ftot_pelagic", "Ftot_integrated", "C_littoral")]
SUM<-c("size.bin"=NA,"pCO2Center"=NA,"propCcenter_littoral"=NA, "Ftot_pelagic"=sum(d$Ftot_pelagic),"Ftot_integrated"=sum(d$Ftot_integrated), "C_littoral"=NA)
summary_table2<-rbind(summary_table, SUM)

round_df <- function(df, digits){
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))

  df[,nums] <- round(df[,nums], digits = digits)

  (df)
}
summary_table2 <- round_df(summary_table2, 3)
write.csv(summary_table2, file = "CO2_global_model.csv")
