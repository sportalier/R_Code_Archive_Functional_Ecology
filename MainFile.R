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
# The following code is used to run the model defined in the paper
# It calls several functions defined in other files
# It defines body masses for n species
# Then, the code computes for each species: 
# 1. A vector of traits (see fnSpecies function for details)
# 2. An array of values for capture sequence (see fnMotion file for details)
# Several n*n matrices are computed to define all aspects of predator-prey relationship
# 1. Search cost
# 2. Capture cost
# 3. Handling cost
# 4. Energy content
# 5. Net gains
# 6. Net gains per mass
#
# Details on calculation and full references for values are provided in the main text or in the supplementary methods
#
##############################################

source('fnSpecies.R')
source('fnMotion.R')

# Physical parameters (sea water at 20 C) (change values for air, see main text)
BodyDensity=1080 # (kg.m-3)
MediumDensity=1026.95 # (kg.m-3)
DynamicViscosity=0.00139 # (N.s.m-2) 
Gravity=9.8 # (m.s-2)

param=c(BodyDensity,MediumDensity,DynamicViscosity,Gravity)

######## Species body mass
# The example is done over 4 species only
# But any number of species should work
BodyMass=c(1e-8,1e-5,1e-3,1e1)
NumberofSpecies=4

######## Species traits
SpeciesTraits=vector('list',NumberofSpecies)

for(i in 1:NumberofSpecies){
  SpeciesTraits[[i]]=fnSpecies(BodyMass[i],param)
}

######## Array for capture sequence
SpeciesArray=vector('list',NumberofSpecies)

for(i in 1:NumberofSpecies){
  # Retrieve species traits
  Sp=SpeciesTraits[[i]];
  # Maximal muscular output
  Maxforce=55*Sp[2][[1]];
  # Time step for computation
  TimeStep=Sp[5][[1]]/100;
  # Switch variable for fnVelocity function
  # 2: the output of fnVelocity returns a vector with horizontal and vertical optimized forces for capture sequence
  i_Switch=2;
  # Boundaries of parameter space for optimization: 0 < force <= Maxforce
  lower=0;
  upper=Maxforce;
  # Vector for fnVelocity function (see fnMotion file for details)
  p=c(Sp[1][[1]],Sp[2][[1]],Sp[3][[1]],Sp[4][[1]],Sp[5][[1]],Gravity,MediumDensity,DynamicViscosity,TimeStep,i_Switch,lower,upper);
  # Starting points for vertical and horizontal forces
  x0=c(Maxforce/5, Maxforce/5);
  # Optimization: returns optimized vertical and horizontal forces (x$par)
  x=optim(par=x0,fn=fnVelocity,p=p,control=c(maxit=10000),method="Nelder-Mead");
  # i_Switch = 4: fnVelocity will retun an array with speed, distance, work , and time, for optimized force values
  i_Switch=4;
  p=c(Sp[1][[1]],Sp[2][[1]],Sp[3][[1]],Sp[4][[1]],Sp[5][[1]],Gravity,MediumDensity,DynamicViscosity,TimeStep,i_Switch,lower,upper);
  SpeciesArray[[i]]=fnVelocity(x$par,p);
}

######## field metabolic rate
fieldmetab=12.51*BodyMass^0.75

######## Searching cost
# Defines a matrix for searching cost (output)
# and a intermediate matrix for encounter rate (used for searching time)  
searchingcost=matrix(NA,nrow=NumberofSpecies,ncol=NumberofSpecies)
EncounterRate=matrix(NA,nrow=NumberofSpecies,ncol=NumberofSpecies)

for (i in 1:NumberofSpecies){
  Predator=SpeciesTraits[[i]]
  for (j in 1:NumberofSpecies){
    Prey=SpeciesTraits[[j]]
    ### prey abundance = 1% of space volume (see main text)
    Abundance=round(0.01/Prey[1][[1]])
    # If prey is very large, the minimal abundance should be 1 individual (not a fraction of individual)
    if (Abundance<1.0){
      Abundance=1
    }
    ### encounter rate
    # If the predator moves faster than the prey
    if (Predator[10][[1]]>Prey[10][[1]]){
      EncounterRate[i,j]=(Abundance*(pi*Predator[9][[1]]^2.0)*((Prey[10][[1]]^2.0)+3.0*(Predator[10][[1]]^2.0)))/(3.0*Predator[10][[1]])
    }else{ # if the prey moves faster than the predator
      EncounterRate[i,j]=(Abundance*(pi*Predator[9][[1]]^2.0)*((Predator[10][[1]]^2.0)+3.0*(Prey[10][[1]]^2.0)))/(3.0*Prey[10][[1]])
    }
    ### searching cost
    # cost = (motion cost at species-specific speed per time + metabolic rate per time) * searching time
    # searching time = 1 / encounter rate
    searchingcost[i,j]=(Predator[11][[1]]+fieldmetab[i])*(1/EncounterRate[i,j])
  }
}

######## Capture cost
# Defines a matrix for capture cost (output)
# and a intermediate matrix for capture probability 
# total capture cost = (capture cost per attempt + metabolic cost) * number of attempt
# number of attempt = 1 / capture probability
capturecost=matrix(NA,nrow=NumberofSpecies,ncol=NumberofSpecies)
captureprobability=matrix(NA,nrow=NumberofSpecies,ncol=NumberofSpecies)

for (i in 1:NumberofSpecies){
  PredatorArray=SpeciesArray[[i]]
  Predator=SpeciesTraits[[i]]
  for (j in 1:NumberofSpecies){
    PreyArray=SpeciesArray[[j]]
    Prey=SpeciesTraits[[j]]
    # prey detection distance (initial distance when chase begins)
    DetectiondistancePrey=Prey[9][[1]]
    # time vectors
    predatortime=PredatorArray[,4]
    preytime=PreyArray[,4]
    # time steps
    predatortimestep=predatortime[2]-predatortime[1]
    preytimestep=preytime[2]-preytime[1]
    # distance vectors
    predatordistance=PredatorArray[,2]
    preydistance=PreyArray[,2]
    # determine maximal time
    maxtime=max(predatortime[length(predatortime)],preytime[length(preytime)])
    
    ### adjust for different duration of the capture sequence between the predator and its prey
    # the predator can run over a longer preriod of time
    if (predatortime[length(predatortime)] > preytime[length(preytime)]){
      # missing time values
      missingtime=predatortime[length(predatortime)]-preytime[length(preytime)]
      missingtimevector=rep(0,missingtime,by=preytimestep)
      missingtimevector=missingtimevector+preytime[length(preytime)]
      preytime=c(preytime,missingtimevector)
      # missing distance values (prey does not move anymore)
      preydistance=c(preydistance,rep(preydistance[length(preydistance)],length(missingtimevector)))
    }else{ # the prey can run over a longer preriod of time
      # missing time values
      missingtime=preytime[length(preytime)]-predatortime[length(predatortime)]
      missingtimevector=rep(0,missingtime,by=predatortimestep)
      missingtimevector=missingtimevector+predatortime[length(predatortime)]
      predatortime=c(predatortime,missingtimevector)
      # missing distance values (predator does not move anymore)
      predatordistance=c(predatordistance,rep(predatordistance[length(predatordistance)],length(missingtimevector)))
    }
    
    ### find a crossing point
    preydistance=preydistance+DetectiondistancePrey
    # fit polynomial functions
    P1=lm(predatordistance~poly(predatortime,5))
    P2=lm(preydistance~poly(preytime,5))
    P1f <- function(x) coef(P1)[1] +x*coef(P1)[2]+(x^2)*coef(P1)[3]+(x^3)*coef(P1)[4] +(x^4)*coef(P1)[5]+(x^5)*coef(P1)[6];
    P2f <- function(x) coef(P2)[1] +x*coef(P2)[2]+(x^2)*coef(P2)[3]+(x^3)*coef(P2)[4] +(x^4)*coef(P2)[5]+(x^5)*coef(P2)[6];
    # check if coefficients are similar (i.e., parallel curves)
    i_check=abs(coef(P1)-coef(P2))
    # parallel curves (no crossing point)
    if (max(i_check[2:6])<1e-15){
      distanceroot=NA
    }else{ # if not parallel
      # determine the first root
      root=tryCatch({
        mod=optimize(f=function(x) abs(P2f(x)-P1f(x)), interval=c(0,1000),tol=1e-6)
        root=mod$minimum[[1]]
      }, error=function(err){
        # totally unrealistic value in case or error
        root=1e6
      })
      # check whether or not curves cross (a root exists)
      if (root<1e6){
        R1=P1f(root)
        R2=P2f(root)
        distanceroot=R2-R1
        distanceroot=distanceroot[[1]]
      }else{
        distanceroot=NA
      }
    }
    # check if cross at the beginning (predator size > prey detection distance)
    if (P1f(0)>P2f(0)){ 
      distanceroot=0;
    }
    # cross or close enough
    threshold=Predator[3][[1]]/100
    if (is.na(distanceroot)==F && distanceroot < threshold){
      # the predator can run over this period of time
      if (root<=predatortime[length(predatortime)]){
        indextimepredator=which(predatortime>=root)[1]
        indextimeprey=which(preytime>=root)[1]
        capturecostperattempt=PredatorArray[indextimepredator,3]+fieldmetab[i]*root
        # predator speed > 0
        if (PredatorArray[indextimepredator,1]>0){
          captureprobability[i,j]=1.0/(1.0+(PreyArray[indextimeprey,1]/PredatorArray[indextimepredator,1]))
          capturecost[i,j]=capturecostperattempt/captureprobability[i,j]
        }else{ # predator speed = 0
          captureprobability[i,j]=0.0
          capturecost[i,j]=NA
        }
        
      }else{ # the predator stops before
        captureprobability[i,j]=0.0
        capturecost[i,j]=NA
      }
    }else{ # no crossing point
      captureprobability[i,j]=0.0
      capturecost[i,j]=NA
    }
  }
}

######## Handling cost
# Defines a matrix for handling cost (output)
handlingcost=matrix(NA,nrow=NumberofSpecies,ncol=NumberofSpecies)

for (i in 1:NumberofSpecies){
  Predator=SpeciesTraits[[i]]
  for (j in 1:NumberofSpecies){
    Prey=SpeciesTraits[[j]]
    #### handling time calculation
    ### consumption time = number of bites * bite time
    # if predator bite size < prey size
    if (Predator[7][[1]]<Prey[2][[1]]){
      consumptiontime=(Prey[2][[1]]/Predator[7][[1]])*Predator[8][[1]]
    }else{ # prey consumed in one bite
      consumptiontime=Predator[8][[1]]
    }
    ### digestion time
    digestiontime=2.3e4*(Prey[2][[1]]/Predator[2][[1]])*(Predator[2][[1]]^0.25)#2.3e4*(Prey[2][[1]]/Predator[7][[1]])*(Predator[2][[1]]^0.25)
    ### handling time
    handlingtime=consumptiontime+digestiontime
    
    #### handling cost calculation
    # Cumulated body volume
    Totalvolume=Predator[1][[1]]+Prey[1][[1]]
    # Cumulated body mass
    Totalmass=Predator[2][[1]]+Prey[2][[1]]
    # Max body radius
    Maxradius=max(Predator[3][[1]],Prey[3][[1]])
    # Max body section surface
    Maxbodysurface=max(Predator[4][[1]],Prey[4][[1]])
    # Time step for computation
    TimeStep=Predator[5][[1]]/100
    # Predator body length
    BodyLength=Predator[3][[1]]*2
    # Vector for fnHandlingMotion function (see corresponding file for details)
    p=c(Totalvolume,Totalmass,Maxradius,Maxbodysurface,Predator[5][[1]],Gravity,MediumDensity,DynamicViscosity,TimeStep,BodyLength)
    # Vertical force
    VerticalForce=55*Predator[2][[1]]
    # Work per time
    Workpt=fnHandlingMotion(VerticalForce,p)
    ### handling cost
    # function returns a cost (feasible interaction)
    if (is.na(Workpt)==F){
      handlingcost[i,j]=(Workpt+fieldmetab[i])*handlingtime
    }else{ # function returns NA (non feasible interaction)
      handlingcost[i,j]=NA
    }
  }
}

######## Energy given by the prey (in a matrix form)
energy=matrix(NA,nrow=NumberofSpecies,ncol=NumberofSpecies)

for (j in 1:NumberofSpecies){
  Prey=SpeciesTraits[[j]]
  energy[,j]=Prey[6][[1]]
}

######## Calculation of net energetic gain for each predator-prey interaction
# Defines predation net gain
# and predation net gain per kg of predator
netgain=energy-(searchingcost+capturecost+handlingcost)

netgainperkg=matrix(NA,nrow=NumberofSpecies,ncol=NumberofSpecies)
for (i in 1:NumberofSpecies){
  Predator=SpeciesTraits[[i]]
  netgainperkg[i,]=netgain[i,]/Predator[2][[1]]
}


