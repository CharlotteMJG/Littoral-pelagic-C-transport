# Functions supporting the transport model scripts

# Equations based on:
# Johansson, H., Brolin, A.A., Håkanson, L., 2007. New Approaches to the
#    Modelling of Lake Basin Morphometry. Environmental Modeling & Assessment
#    12(3) 213-228 doi:10.1007/s10666-006-9069-z.

# Function to calculate depth based on radius, with a maximum surface layer depth
H_r_johansson <- function(r, r_max, max_depth, max_surf_layer, Hd){
  depth = max_depth * log((r / r_max) * (1 - 1/Hd) + 1/Hd, base = 1/Hd)
  
  if(depth >= max_surf_layer){
    return(max_surf_layer)
  }else{
    return(depth)
  }
}

# Analytical solution of H_r_johansson
dHdr_johansson <- function(r, r_max, max_depth, max_surf_layer, Hd){
  depth = max_depth * log((r / r_max) * (1 - 1/Hd) + 1/Hd, base = 1/Hd)
  dHdr = max_depth * ((1 - 1/Hd) / r_max) / (((r / r_max) * (1 - 1/Hd) + 1/Hd) * log(1/Hd))
  
  if(depth >= max_surf_layer){
    return(0.0)
  }else{
    return(dHdr)
  }
}

calc_hd <- function(Vd){
  # See equations 15-18 in Johansson et al.
  b = log10(Vd)
  
  if(Vd >= 0.2 & Vd <= 0.55){
    hd = 10^(1.6 * b^4 - 2.1 * b^3 - 1.2 * b^2 - 3.92 * b)
  }else if(Vd >= 0.55 & Vd <= 1.85){
    hd = 10^(-34 * b^5 - 18 * b^4 - 6.3 * b^3 - 1.9 * b^2 - 4 * b)
  }else if(Vd >= 1.85 & Vd <= 2.40){
    hd = 10^(-19261 * b^5 + 27179 * b^4 - 15506 * b^3 + 4432.6 * b^2 - 639.51 * b + 36.442)
  }else if(Vd >= 2.40 & Vd <= 2.70){
    hd = 10^(-3077645 * b^5 + 6028981 * b^4 - 4728714.5 * b^3 + 1855729.6 * b^2 - 364329.44 * b + 28621.94)
  }else{
    stop("No equation provided to calculate Hd")
  }
}
