##############################################
##############################################
#
# R Code supplementing the paper
# The mechanics of predator-prey interactions: first principles of physics predict predator-prey size ratios
# by Portalier, Fussmann, Loreau, Cherif

# April 2018
#
##############################################
##############################################
#

#### READ ME #################################
#
# The following code defines the function fnSpecies
# that calculates species traits
# The function takes body mass and a vector of several parameters

#### Physical parameters taken as arguments:
# 1. Body density (kg.m-3)
# 2. Medium density (kg.m-3)
# 3. Medium dynamic viscosity (N.s.m-2)
# 4. Acceleration due to gravity (m.s-2)
#
#### Parameters defined within the function
# see main text or supplementary methods for references
# 1. Ash-free dry mass to wet mass ratio
# 2. Energy (J.kg-1) to ash-free dry mass ratio
# 3. Bite diameter (m) at reference size
# 4. Reference mass (kg) for bite size
# 5. Detection distance (m) at reference size
# 6. Reference mass (kg) for detection distance
#
#### Function returns a vector with
# 1. Body volume (m3)
# 2. Body mass (kg)
# 3. Body radius (m)
# 4. Section surface of the body (m2)
# 5. Stroke period (s)
# 6. Energetic content (J)
# 7. Bite size (kg)
# 8. Bite time (s)
# 9. Detection distance (m)
# 10. Species-specific speed (m.s-1)
# 11. Motion cost per time at species-specific speed (J.s-1)
#
##############################################

fnSpecies = function(BodyMass,param){
  
  ## unwrap physical parameters
  BodyDensity=param[1];
  MediumDensity=param[2];
  DynamicViscosity=param[3];
  Gravity=param[4];
  
  ## Define other parameters
  # Ash-free dry mass to wet mass ratio
  DryMassRatio=0.16;
  # Energy (J.kg-1) to ash-free dry mass ratio
  EnergyDryMassRatio=23000000; 
  # Bite diameter (m) at reference size
  ReferenceBiteSize=0.00026; 
  # Reference mass (kg) for bite size
  BiteSizeReferenceMass=2.9; 
  # Detection distance (m) at reference size
  ReferenceDetectionDistance=0.255; 
  # Reference mass (kg) for detection distance
  ReferenceSizeDetectiondistance=0.0376;
  
  # response vector
  Species=rep(0,11);
  
  # Body volume
  Species[1]=BodyMass/BodyDensity;
  # Body mass
  Species[2]=BodyMass;
  # Body radius
  Species[3]=(Species[1]*3.0/4.0/pi)^(1.0/3.0);
  # Section surface
  Species[4]= (Species[3]^2.0)*pi;
  # Stroke period
  Species[5]=0.35*(BodyMass^0.25)#0.008*(BodyMass^0.25);
  # Energy content
  Species[6]=BodyMass*DryMassRatio*EnergyDryMassRatio;
  # Bite size
  BiteRadius=(ReferenceBiteSize*((BodyMass/BiteSizeReferenceMass)^0.32))/2.0;
  Species[7]=4.0/3.0*(BiteRadius^3)*pi*BodyDensity;
  # Bite time
  Species[8]= 0.1*(Species[7])^2;
  # Detection distance
  Species[9]=ReferenceDetectionDistance*((Species[1]/ReferenceSizeDetectiondistance)^0.49);
  ## Species-specific speed calculation
  # Maximal muscular output
  Maxforce=55*BodyMass;
  # Time step for computation
  TimeStep=Species[5]/100;
  # Switch variable for fnVelocity function
  # i_Switch = 1: the output of fnVelocity returns a vector with horizontal and vertical optimized forces
  i_Switch=1;
  # Boundaries of parameter space for optimization: 0 < force <= Maxforce
  lower=0;
  upper=Maxforce;
  # Vector for fnVelocity function (see corresponding file for details)
  p=c(Species[1],BodyMass,Species[3],Species[4],Species[5],Gravity,MediumDensity,DynamicViscosity,TimeStep,i_Switch,lower,upper);
  # Starting points for vertical and horizontal forces
  x0=c(Maxforce/10, Maxforce/10);
  # Optimization: returns optimized vertical and horizontal forces (x$par)
  x=optim(par=x0,fn=fnVelocity,p=p,control=c(maxit=1000),method="Nelder-Mead");
  # i_Switch = 3: fnVelocity will retun a vector with speed and motion cost, for optimized force values
  i_Switch=3;
  p=c(Species[1],BodyMass,Species[3],Species[4],Species[5],Gravity,MediumDensity,DynamicViscosity,TimeStep,i_Switch,lower,upper);
  Res=fnVelocity(x$par,p);
  # Species-specific speed
  Species[10]=Res[1];
  # Motion cost per time
  Species[11]=Res[2];
  return(Species)
}
