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
# The following code defines 10 functions used in the calculation of motion speed and cost 
# fnDragCoef
# fnActiveAscending
# fnPassiveAscending
# fnRootPassiveMotion
# fnPassiveDescending
# fnRootPassiveDescending
# fnActiveForwardMotion
# fnPassiveForwardMotion
# fnVelocity 
# fnHandlingMotion
#
#
#### fnDragCoef is used to compute the drag coefficient
#
# The function takes body radius, speed, and a vector of physical parameters (Medium density and Medium dynamic viscosity)
# 
# The function returns a single value: drag coefficient
#
#
#### fnActiveAscending is used to compute the derivative of speed with respect to time during active ascending phase of motion (vertical plan)
# This ODE is called by the solver lsoda within fnVelocity function
#
# The function takes time, speed, and a vector of parameters
# Parameters: Body mass, Body volume, Body radius, Body section surface, Vertical muscular force, Medium density, Medium dynamic viscosity, and Acceleration due to gravity
#
# The function returns the value of the derivative of speed with respect to time (acceleration) as a list
#
#
#### fnPassiveAscending is used to compute the derivative of speed with respect to time during the passive (inertial) ascending phase (vertical plan)
# This ODE is called by the solver lsodar within fnVelocity function
#
# The function takes time, speed, and a vector of parameters
# Parameters: Body mass, Body volume, Body radius, Body section surface, Medium density, Medium dynamic viscosity, and Acceleration due to gravity
#
# The function returns the value of the derivative of speed with respect to time (acceleration) as a list
#
#### fnRootPassiveMotion is used to trigger a stop event to the computation (for passive ascending and passive forward motion phases)
# This root function is called by the solver lsodar within fnVelocity function
# 
# The function takes time, speed, and a vector of parameters (same as fnPassiveAscending, fnPassiveDescending and fnPassiveForwardMotion functions)
#
# It returns 0 if speed = 0: the individual body stops, which triggers the end of computation
# It retuns 200 otherwise
# 
#### fnPassiveDescending is used to compute the derivative of speed with respect to time during the passive (inertial) descending phase (vertical plan)
# This ODE is called by the solver lsodar within fnVelocity function
# The function takes time, speed, and a vector of parameters (same as fnPassiveAscending and fnPassiveForwardMotion functions)
#
# The function returns the value of the derivative of speed with respect to time (acceleration) as a list
#
#### fnRootPassiveDescending is used to trigger a stop event to the computation (for passive descending motion phase)
# This root function is called by the solver lsodar within fnVelocity function
# 
# The function takes time, speed, and a vector of parameters (same as fnPassiveAscending, fnPassiveDescending and fnPassiveForwardMotion functions)
#
# It returns 0 if derivative of speed (acceleration) = 0: the system is assumed to be at steady state
# It returns 200 otherwise
#
#### fnActiveForwardMotion is used to compute the derivative of speed with respect to time during the active forward motion phase (horizontal plan)
#
# The function takes time, speed, and a vector of parameters
# Parameters: Body mass, Body volume, Body radius, Body section surface, Horizontal muscular force, Medium density, Medium dynamic viscosity, and Acceleration due to gravity
#
# The function returns the value of the derivative of speed with respect to time (acceleration) as a list
#
#### fnPassiveForwardMotion is used to compute the derivative of speed with respect to time during the passive (inertial) forward motion phase (horizontal plan)
#
# The function takes time, speed, and a vector of parameters (same as fnPassiveAscending and fnPassiveDescending functions)
#
# The function returns the value of the derivative of speed with respect to time (acceleration) as a list
#
#### fnVelocity is used to compute speed, forces, or motion cost
# It is the only function called from outside
#
# The function takes two vectors as arguments (x and p)
#
## x is a vector containing vertical and horizontal forces 
# 1. Vertical muscular force (N)
# 2. Horizontal muscular force (N)
#
## p is vector containing biological and physical parameters
# 1. Body volume (m3)
# 2. Body mass (kg)
# 3. Body radius (m)
# 4. Body section surface (m2)
# 5. Stroke period (s)
# 6. Acceleration due to gravity (m.s-2)
# 7. Medium density (kg.m-3)
# 8. Medium dynamic viscosity (N.s.m-2)
# 9. Time step for calculation (s)
# 10. Switch parameter (see below)
# 11. Lower bound for optimization
# 12. Upper bound for optimization
#
## The function returns different outputs according to the value of the switch parameter
#
## If i_swith = 1
# fnVelocity returns the ratio between work and horizontal speed: this value has to be minimized for species-specific speed (search sequence)
# This option is used during optimization procedure for species-specific speed
#
## If i_swith = 2
# fnVelocity returns a vector with horizontal speed and motion cost 
# It is used by fnSpecies function to define species-specific speed and motion cost at species-specific speed, after optimization
#
## If i_swith = 3
# fnVelocity returns 1 / total horizontal distance travelled: this value has to be minimized for predator jump (capture sequence)
# This option is used during optimization procedure of capture sequence
#
## If i_swith = 4
# fnVelocity returns an array of 4 columns used during capture sequence (used after optimization procedure)
# This array is used for the calculation of capture cost
# The 4 colums give for each time step:
# 1. Instantaneous horizontal speed (m.s-1)
# 2. Cumulated horizontal distance travelled since the beginning (m)
# 3. Cumulated work (J)
# 4. Cumulated time (s)
#
# 
#### fnHandlingMotion is used to compute handling cost per time specificaly
# 
# It uses the same code than fnVelocity function with 2 important differences
#
# 1. Only vertical motion is considered. The predator is assumed not to move forwards
# 2. The vertical distance travelled corresponds to body length (upwards and downwards)
# overall, the predator remains at the same vertical position
#
# The function returns the motion cost per time unit (J.s-1)
#
#
##############################################
library(deSolve)

#### fnDragCoef function
fnDragCoef=function(Radius, Speed,param){
  
  # unwrap physical parameters
  MediumDensity=param[1]
  DynamicViscosity=param[2]
  
  # Avoid negative values
  if (Speed<0.0){
    Speed=-Speed
  }
  # Reynolds' number
  Reynolds = MediumDensity*Speed*Radius*2/DynamicViscosity
  # Drag coefficient
  if (Reynolds>0.0){
    DragCoef=(0.352+(0.124+24.0/Reynolds)^0.5)^2.0
  }else{ # in case of speed = 0
    DragCoef=0.0
  }
  return(DragCoef)
}

#### fnActiveAscending function
fnActiveAscending=function(t,x,parms){
  
  # unwrap parameters
  BodyMass=parms[1]
  BodyVolume=parms[2]
  BodyRadius=parms[3]
  BodySectionSurface=parms[4]
  VerticalForce=parms[5]
  MediumDensity=parms[6]
  DynamicViscosity=parms[7]
  Gravity=parms[8]
  
  # compute drag coefficient
  paramDragCoef=c(MediumDensity,DynamicViscosity)
  DragCoef=fnDragCoef(BodyRadius,x,paramDragCoef)
  
  # compute derivative
  y=VerticalForce/BodyMass-Gravity+MediumDensity*BodyVolume*Gravity/BodyMass-0.5*MediumDensity*BodySectionSurface*DragCoef*(x^2.0)/BodyMass;
  return(list(y))
}

#### fnPassiveAscending function
fnPassiveAscending=function(t,x,parms){
  
  # unwrap parameters
  BodyMass=parms[1]
  BodyVolume=parms[2]
  BodyRadius=parms[3]
  BodySectionSurface=parms[4]
  MediumDensity=parms[5]
  DynamicViscosity=parms[6]
  Gravity=parms[7]
  
  # compute drag coefficient
  paramDragCoef=c(MediumDensity,DynamicViscosity)
  DragCoef=fnDragCoef(BodyRadius,x,paramDragCoef)
  
  # compute derivative
  y=-Gravity+MediumDensity*BodyVolume*Gravity/BodyMass-0.5*MediumDensity*BodySectionSurface*DragCoef*(x^2.0)/BodyMass;
  return(list(y))
}

#### fnRootPassiveMotion function
fnRootPassiveMotion=function(t,x,parms){
  # check if x<=0: individual body stops
  if (x<=0){
    # triggers end of computation
    y=0
  }else{
    y=200
  }
  return(y)
}

#### fnPassiveDescending function
fnPassiveDescending=function(t,x,parms){
  
  # unwrap parameters
  BodyMass=parms[1]
  BodyVolume=parms[2]
  BodyRadius=parms[3]
  BodySectionSurface=parms[4]
  MediumDensity=parms[5]
  DynamicViscosity=parms[6]
  Gravity=parms[7]
  
  # compute drag coefficient
  paramDragCoef=c(MediumDensity,DynamicViscosity)
  DragCoef=fnDragCoef(BodyRadius,x,paramDragCoef)
  
  # compute derivative
  y=Gravity-MediumDensity*BodyVolume*Gravity/BodyMass-0.5*MediumDensity*BodySectionSurface*DragCoef*(x^2.0)/BodyMass;
  return(list(y))
}

#### fnRootPassiveDescending function
fnRootPassiveDescending=function(t,x,parms){
  # check value of the derivative
  y1=fnPassiveDescending(t,x,parms)
  # if derivative < 1e-9: triggers end of computation (the system is assumed to be at steady state)
  if (y1[[1]]<1e-9){
    y=0
  }else{
    y=200
  }
  return(y)
}

#### fnActiveForwardMotion function
fnActiveForwardMotion=function(t,x,parms){
  
  # unwrap parameters
  BodyMass=parms[1]
  BodyVolume=parms[2]
  BodyRadius=parms[3]
  BodySectionSurface=parms[4]
  HorizontalForce=parms[5]
  MediumDensity=parms[6]
  DynamicViscosity=parms[7]
  Gravity=parms[8]
  
  # compute drag coefficient
  paramDragCoef=c(MediumDensity,DynamicViscosity)
  DragCoef=fnDragCoef(BodyRadius,x,paramDragCoef)
  
  # compute derivative
  y=HorizontalForce/BodyMass-0.5*MediumDensity*BodySectionSurface*DragCoef*(x^2.0)/BodyMass;
  return(list(y))
}

#### fnPassiveForwardMotion function
fnPassiveForwardMotion=function(t,x,parms){
  
  # unwrap parameters
  BodyMass=parms[1]
  BodyVolume=parms[2]
  BodyRadius=parms[3]
  BodySectionSurface=parms[4]
  MediumDensity=parms[5]
  DynamicViscosity=parms[6]
  Gravity=parms[7]
  
  # compute drag coefficient
  paramDragCoef=c(MediumDensity,DynamicViscosity)
  DragCoef=fnDragCoef(BodyRadius,x,paramDragCoef)
  
  # compute derivative
  y=-0.5*MediumDensity*BodySectionSurface*DragCoef*(x^2.0)/BodyMass;
  return(list(y))
}

#### fnVelocity function
fnVelocity=function(x,p){
  library(deSolve)
  # unwrap x
  VerticalForce=x[1]
  HorizontalForce=x[2]
  # unwrap p
  BodyVolume=p[1]
  BodyMass=p[2]
  BodyRadius=p[3]
  BodySectionSurface=p[4]
  StrokePeriod=p[5]
  Gravity=p[6]
  MediumDensity=p[7]
  DynamicViscosity=p[8]
  TimeStep=p[9]
  i_Switch=p[10]
  lowerbound=p[11]
  upperbound=p[12]
  #TimeStep=TimeStep/2
  
  # Parameter vector for fnDragCoef function
  paramdrag=c(MediumDensity,DynamicViscosity)
  
  # Define a flag to avoid going out of parameter bounds
  # 0: values within the defined range
  # 1: values out of bounds
  flag=0
  
  # i_Switch = 0 or 2: species-specific speed
  # Check for validity of values 
  if (i_Switch==0 || i_Switch==2){
    total=x[1]+x[2]
    if (total>upperbound){
      flag=1 
    }
  }
  
  if (x[1]>upperbound || x[2]>upperbound){
    flag=1
  }
  if (x[1]<lowerbound || x[2]<lowerbound){
    flag=1
  }
  
  if (flag==0){
    ##### Vertical component
    ## Active phase
    # initial speed = 0
    x0=0 
    # time vector for computation
    timevector=seq(0,StrokePeriod,by=TimeStep)
    # vector of parameters for computation
    paramactiveascending=c(BodyMass,BodyVolume,BodyRadius,BodySectionSurface,VerticalForce,MediumDensity,DynamicViscosity,Gravity)
    # solve speed with respect to time during ascending active phase (stroke period)
    res1=lsoda(x0,timevector,fnActiveAscending,parms=paramactiveascending,rtol=1e-8,atol=1e-8,maxsteps=10000)
    # vertical speed at the end of the active phase
    verticalspeed=res1[nrow(res1),2]
    verticalspeed=verticalspeed[[1]]
    # Instaneous distance travelled
    instantdistancevertical=res1[,2]*TimeStep
    # Cumulated vertical distance
    totaldistance=cumsum(instantdistancevertical)
    # Total vertical distance travelled during active phase
    activedistancetravelled=totaldistance[length(totaldistance)]
    # Check if vertical force allows effective motion (or if it is too weak)
    if (activedistancetravelled<0){
      flag=1
    }
  }
  
  if (flag==0){
    ## Inertial ascending phase
    # Speed = vertical speed at the end of the active phase
    x0=verticalspeed
    timevector=seq(0,TimeStep*10000,by=TimeStep)
    # Vector of parameters for active phases (without forces)
    parampassivephase=c(BodyMass,BodyVolume,BodyRadius,BodySectionSurface,MediumDensity,DynamicViscosity,Gravity)
    # solve speed with respect to time during ascending passive (inertial) phase
    # A stop event is triggered when speed = 0 (individual body stops)
    if (x0 >0){
      res2=lsodar(x0,timevector,fnPassiveAscending,parms=parampassivephase,rootfunc=fnRootPassiveMotion,rtol=1e-8,atol=1e-8,maxsteps=10000)
      # Total vertical distance travelled during intertial ascending phase 
      instantdistancevertical=res2[,2]*TimeStep
      totaldistance=cumsum(instantdistancevertical)
      passivedistanceupwards=totaldistance[length(totaldistance)-1]
      # Duration of inertial ascending phase
      passivetimeupwards=res2[nrow(res2)-1,1]
      passivetimeupwards=passivetimeupwards[[1]]
    }else{
      passivedistanceupwards=0;
      passivetimeupwards=0;
    }
    # Total distance travelled during ascending (active + passive) phases 
    totalascendingdistance=activedistancetravelled+passivedistanceupwards
    
    ## Inertial descending phase
    x0=0
    timevector=seq(0,TimeStep*10000,by=TimeStep)
    # solve speed with respect to time during descending passive (inertial) phase
    # A stop event is triggered when the derivative of speed = 0 (steady state: terminal sinking velocity)
    res3=lsodar(x0,timevector,fnPassiveDescending,parms=parampassivephase,rootfunc=fnRootPassiveDescending,rtol=1e-8,atol=1e-8,maxsteps=10000)
    # Total vertical distance travelled during intertial descending phase 
    instantdistancevertical=res3[,2]*TimeStep
    totaldistance=cumsum(instantdistancevertical)
    # Check if distance travelled > total ascending distance travelled
    if (totaldistance[length(totaldistance)] > totalascendingdistance){
      # Search when distance travelled = total ascending distance travelled
      index=which(totaldistance >= totalascendingdistance)
      passivetimedownwards=res3[index[1],1]
      passivetimedownwards=passivetimedownwards[[1]]
    }else{
      # Extrapolate when distance travelled = total ascending distance travelled (knowing that the individual body has reached its terminal velocity)
      remainingdistance=totalascendingdistance-totaldistance[length(totaldistance)]
      remainingtime=remainingdistance/res3[nrow(res3),2][[1]]
      passivetimedownwards=res3[nrow(res3),1][[1]]+remainingtime
    }
    # Total duration of vertical motion
    totaltime=StrokePeriod+passivetimeupwards+passivetimedownwards
    
    ##### Horizontal component
    ## Active phase
    x0=0
    timevector=seq(0,StrokePeriod,by=TimeStep)
    # vector of parameters (same as active ascending phase, except Horizontal muscular force)
    paramactiveforward=c(BodyMass,BodyVolume,BodyRadius,BodySectionSurface,HorizontalForce,MediumDensity,DynamicViscosity,Gravity)
    # solve speed with respect to time during active phase (stroke period)
    res4=lsoda(x0,timevector,fnActiveForwardMotion,parms=paramactiveforward,rtol=1e-8,atol=1e-8,maxsteps=10000)
    # horizontal speed at the end of the active phase
    horizontalspeed=res4[nrow(res4),2]
    horizontalspeed=horizontalspeed[[1]]
    # Total horizontal distance travelled during active motion phase forwards
    instantdistancehorizontal=res4[,2]*TimeStep
    totaldistance=cumsum(instantdistancehorizontal)
    activedistancetravelledforwards=totaldistance[length(totaldistance)]
    # Check if horizontal force allows effective motion (or if it is too weak)
    if (activedistancetravelledforwards<0){
      flag=1
    }
  }
  
  if (flag==0){
    ## Inertial phase
    x0=horizontalspeed
    timevector=seq(0,totaltime,by=TimeStep)
    if (x0 >0){
      # solve speed with respect to time during passive (inertial) motion phase forwards
      # A stop event is triggered when speed = 0 (individual body stops)
      res5=lsodar(x0,timevector,fnPassiveForwardMotion,parms=parampassivephase,rootfunc=fnRootPassiveMotion,rtol=1e-8,atol=1e-8,maxsteps=10000)
      # Total horizontal distance travelled during active motion phase forwards
      instantdistancehorizontal=res5[,2]*TimeStep
      totaldistance=cumsum(instantdistancehorizontal)
      passivedistancetravelledforwards=totaldistance[length(totaldistance)]
    }else{
      passivedistancetravelledforwards=0
    }
    # Total distance travelled forwards
    totaldistanceforwards=activedistancetravelledforwards+passivedistancetravelledforwards
  }
  
  #### Outputs
  
  ### Case 1: i_Switch = 1
  # Optimization procedure for species-specific speed
  # The function returns the ratio between work and speed
  if (i_Switch==1){
    if (flag==0){
      #Work=(activedistancetravelled*VerticalForce+activedistancetravelledforwards*HorizontalForce)
      Workpt=(activedistancetravelled*VerticalForce+activedistancetravelledforwards*HorizontalForce)/(totaltime)
      HorizontalSpeed=totaldistanceforwards/totaltime
      #Output=Work/HorizontalSpeed
      Output=Workpt/totaldistanceforwards
    }else{
      # Values out of bounds lead to huge penalty
      Output=1e9
    }
  }
  
  ### Case 2: i_Switch = 2
  # Optimization procedure for capture sequence
  # The function returns - total horizontal distance (i.e., the maximal total horizontal distance)
  if (i_Switch==2){
    if (flag==0){
      Output=-totaldistanceforwards
    }else{
      # Values out of bounds lead to huge penalty
      Output=1e9
    }
  }
  
  ### Case 3: i_Switch = 3
  # After optimization, function returns species-specific speed and motion cost per time
  if (i_Switch==3){
    Workpt=(activedistancetravelled*VerticalForce+activedistancetravelledforwards*HorizontalForce)/(totaltime)
    HorizontalSpeed=totaldistanceforwards/totaltime
    Output=c(HorizontalSpeed, Workpt)
  }
  
  ### Case 4: i_Switch = 4
  # After optimization, function returns an array for capture sequence
  if (i_Switch==4){
    if (flag==0){
      ## Instantneous speed
      instantspeed=c(res4[-1,2],res5[-nrow(res5),2],0)
      ## cumulated horizontal distance
      # active phase
      instantdistancehorizontal=res4[-1,2]*TimeStep
      totaldistanceactive=cumsum(instantdistancehorizontal)
      finaldistanceactive=totaldistanceactive[length(totaldistanceactive)]
      # passive phase
      instantdistancehorizontal=res5[,2]*TimeStep
      totaldistancepassive=cumsum(instantdistancehorizontal)
      totaldistancepassive=totaldistancepassive+finaldistanceactive
      # cumulated distance
      cumulatedhorizontaldistance=c(totaldistanceactive,totaldistancepassive)
      
      ## cumulated work
      cumulatedworkactivephase=totaldistanceactive*HorizontalForce
      finalactivework=cumulatedworkactivephase[length(cumulatedworkactivephase)]
      # passive phase: instaneous work = 0, thus cumulated work = final work at the end of active phase
      cumulatedworkpassivephase=rep(finalactivework,length(totaldistancepassive))
      totalcumulatedwork=c(cumulatedworkactivephase,cumulatedworkpassivephase)
      
      ## cumulated time
      finaltimeactive=res4[nrow(res4),1][[1]]
      cumulatedtimepassive=res5[,1]
      cumulatedtimepassive=cumulatedtimepassive+finaltimeactive
      cumulatedtime=c(res4[-1,1],cumulatedtimepassive)
      
      ## adjust size of speed vector
      maxsize=max(length(instantspeed),length(cumulatedhorizontaldistance),length(totalcumulatedwork),length(cumulatedtime))
      if (length(instantspeed)<maxsize){
        missingdata=maxsize-length(instantspeed)
        instantspeed=c(instantspeed,rep(0,missingdata))
      }
      
      ## data frame
      Output=data.frame(instantspeed,cumulatedhorizontaldistance,totalcumulatedwork,cumulatedtime)
      names(Output)=c('InstantaneousSpeed','CumulatedHorizontalDistance','CumulatedWork','CumulatedTime')
    }else{
      Output=NA
    }
  }
  
  ### return Output
  return(Output)
}

#### fnHandlingMotion function
fnHandlingMotion=function(VerticalForce,p){
  # unwrap p
  BodyVolume=p[1]
  BodyMass=p[2]
  BodyRadius=p[3]
  BodySectionSurface=p[4]
  StrokePeriod=p[5]
  Gravity=p[6]
  MediumDensity=p[7]
  DynamicViscosity=p[8]
  TimeStep=p[9]
  BodyLength=p[10]
  
  # Parameter vector for fnDragCoef function
  paramdrag=c(MediumDensity,DynamicViscosity)
  
  ## Active phase
  flag=0
  # initial speed = 0
  x0=0 
  # time vector for computation
  timevector=seq(0,StrokePeriod,by=TimeStep)
  # vector of parameters for computation
  paramactiveascending=c(BodyMass,BodyVolume,BodyRadius,BodySectionSurface,VerticalForce,MediumDensity,DynamicViscosity,Gravity)
  # preliminary check: is the predator able to produce enough force to lift itself and the prey? Without drag (i.e., speed = 0)
  test=fnActiveAscending(0,0,paramactiveascending)
  # if test < 0: the predator cannot lift its body and that of prey
  if (test<0){
    flag=1
  }
  
  if (flag==0){
    # solve speed with respect to time during ascending active phase (stroke period)
    res1=lsoda(x0,timevector,fnActiveAscending,parms=paramactiveascending,rtol=1e-8,atol=1e-8,maxsteps=10000)
    # vertical speed at the end of the active phase
    verticalspeed=res1[nrow(res1),2]
    verticalspeed=verticalspeed[[1]]
    # Instaneous distance travelled
    instantdistancevertical=res1[,2]*TimeStep
    # Cumulated vertical distance
    totaldistance=cumsum(instantdistancevertical)
    # Total vertical distance travelled during active phase
    activedistancetravelled=totaldistance[length(totaldistance)]
    # Check if vertical force allows effective motion (or if it is too weak)
    if (activedistancetravelled<0){
      flag=1
    }
  }
  
  # Check if distance travelled > body length
  if (flag==0){
    if (activedistancetravelled > BodyLength){
      index=which(totaldistance >= BodyLength)
      activetimeupwards=res1[index[1],1]
      activetimeupwards=activetimeupwards[[1]]
      activedistancetravelled=totaldistance[index[1]]
    }else{
      activetimeupwards=StrokePeriod
    }
  }
  
  ## Inertial ascending phase
  if (flag==0){
    # Speed = vertical speed at the end of the active phase
    x0=verticalspeed
    timevector=seq(0,TimeStep*1000,by=TimeStep)
    # Vector of parameters for active phases (without forces)
    parampassivephase=c(BodyMass,BodyVolume,BodyRadius,BodySectionSurface,MediumDensity,DynamicViscosity,Gravity)
    # solve speed with respect to time during ascending passive (inertial) phase
    # A stop event is triggered when speed = 0 (individual body stops)
    if (x0 >0){
      res2=lsodar(x0,timevector,fnPassiveAscending,parms=parampassivephase,rootfunc=fnRootPassiveMotion,rtol=1e-8,atol=1e-8,maxsteps=10000)
      # Total vertical distance travelled during intertial ascending phase 
      instantdistancevertical=res2[,2]*TimeStep
      totaldistance=cumsum(instantdistancevertical)
      passivedistanceupwards=totaldistance[length(totaldistance)-1]
      # Duration of inertial ascending phase
      passivetimeupwards=res2[nrow(res2)-1,1]
      passivetimeupwards=passivetimeupwards[[1]]
    }else{
      passivedistanceupwards=0;
      passivetimeupwards=0;
    }
    
    # Total distance travelled during ascending (active + passive) phases 
    totalascendingdistance=activedistancetravelled+passivedistanceupwards
    
    ## Inertial descending phase
    x0=0
    timevector=seq(0,TimeStep*1000,by=TimeStep)
    # solve speed with respect to time during descending passive (inertial) phase
    # A stop event is triggered when the derivative of speed = 0 (steady state: terminal sinking velocity)
    res3=lsodar(x0,timevector,fnPassiveDescending,parms=parampassivephase,rootfunc=fnRootPassiveDescending,rtol=1e-8,atol=1e-8,maxsteps=10000)
    # Total vertical distance travelled during intertial descending phase 
    instantdistancevertical=res3[,2]*TimeStep
    totaldistance=cumsum(instantdistancevertical)
    # Check if distance travelled > total ascending distance travelled
    if (totaldistance[length(totaldistance)] > totalascendingdistance){
      # Search when distance travelled = total ascending distance travelled
      index=which(totaldistance >= totalascendingdistance)
      passivetimedownwards=res3[index[1],1]
      passivetimedownwards=passivetimedownwards[[1]]
    }else{
      # Extrapolate when distance travelled = total ascending distance travelled (knowing that the individual body has reached its terminal velocity)
      remainingdistance=BodyLength-totaldistance[length(totaldistance)]
      remainingtime=remainingdistance/res3[nrow(res3),2][[1]]
      passivetimedownwards=res3[nrow(res3),1][[1]]+remainingtime
    }
    
    # Total duration of vertical motion
    totaltime=activetimeupwards+passivetimeupwards+passivetimedownwards
  }
  
  #### Output
  # The predator is able to lift itself and the prey
  if (flag==0){
    # Work per time (J.s-1)
    Output=(activedistancetravelled*VerticalForce)/(totaltime)
  }else{ # The predator cannot lift itself and the prey
    Output=NA
  }
  
  return(Output)
}

