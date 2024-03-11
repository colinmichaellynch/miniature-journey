#Clear the workspace
rm(list = ls())

#load all necessary packages 
library(ggplot2)
library(nlme)
library(ggpubr)
library(fitdistrplus)
library(poweRlaw)

#set the working directory - set this to where you plan on running the code
setwd("~/Documents/Temporal Sampling")

#import data. Columns represent ants, the names give the ant's color code combination
data = read.csv("rawDataTask.csv")

### --- Estimating minimum sample size for accurately reconstructing rare states: the 1% rule --- ###

#set constants
maxSampleSize = 3500 #highest sample size we explore
replicates = 100 #number of simulations we run per sample size
CThreshold = .99 #selection criterion for certainty

#initialize vectors
MinSampleSize = c() #vector which will contain certainty estimate for sample size
MinSampleSizeSimple = c() #vector for the binomial sample size estimate

#first, we iterate through each ant, represented as a column in the dataframe
for(i in 1:ncol(data)){
  
  NumericData = match(data[,i], unique(data[,i])) #convert strings to numbers for faster processing
  TaskNumber = length(unique(data[,i]))
  TheoreticalVariance = 2*(TaskNumber-1) #get theoretical variance for ant i
  TrueCounts = table(NumericData)
  TrueProbs = as.numeric(TrueCounts)/sum(TrueCounts) #find smallest probability that any state will occur
  StatesPresent = 1:TaskNumber
  MinSampleSizeSimple[i] = 1/min(TrueProbs) #binomial estimate, equation 1
  
  certainty = c() #initialize vector which will contain certainty estimates for each sample size
  
  #iterate through all possible sample sizes
  for(j in 1:maxSampleSize){
    
    ChiSquared = c()
    
    for(k in 1:replicates){
      
      #this entire function randomly samples from the behavioral data of ant i and creates a counts table which has as many states as what the ant actually has. This corrects for the fact that there should be a set number of states but any sample may
      RandSample = sample(NumericData, j, replace = TRUE)
      ObservedCounts = table(RandSample)
      SampleStates = as.numeric(names(ObservedCounts))
      if(length(SampleStates) != length(StatesPresent)){
        for(t in 1:length(StatesPresent)){
          SampleStates = as.numeric(names(ObservedCounts))
          SampleStates = SampleStates[is.na(SampleStates)==FALSE]
          BooleanLog = c()
          if((StatesPresent[t] %in% SampleStates)==FALSE){
            BooleanLog = as.numeric(StatesPresent[t]) < SampleStates
            if(all(BooleanLog == TRUE)){
              ObservedCounts = c(0, ObservedCounts)
            } else if(all(BooleanLog == FALSE)) {
              ObservedCounts = c(ObservedCounts, 0)
            } else {
              index = min(which(BooleanLog))-1
              ObservedCounts = c(ObservedCounts[1:index], 0, ObservedCounts[(index+1):length(ObservedCounts)])
            }
          }
        }
      }
      
      #calculate the chi squared value
      ObservedCounts = as.numeric(ObservedCounts)
      ExpectedCounts = j*TrueProbs
      ChiSquared[k] = sum(((ObservedCounts-ExpectedCounts)^2)/ExpectedCounts)
      
    }
    
    #find the certainty at this sample size
    certainty[j] = TheoreticalVariance/var(ChiSquared)
    #these print what stage the simulations are at so you will know how far along we are. We include a counter in every for loop so users know how far along they are in a simulation
    print(i/9)
    print(j/maxSampleSize)
    
  }
  
  fit = loess(certainty ~ c(1:maxSampleSize))
  prediction = predict(fit, c(1:maxSampleSize))
  MinSampleSize[i] = min(which(prediction > CThreshold))
  
  #We choose this iteratation to display figure 2, as this is a fairly prototypical example of what this process looks like
  if(i == 4){
    dataGraph = data.frame(SampleSize = 1:maxSampleSize, Certainty = certainty, MinSampleSize = MinSampleSize[i], cross = CThreshold)
    p1 = ggplot(dataGraph, aes(SampleSize, Certainty)) +
      geom_line() +
      geom_smooth(method = 'loess', formula = 'y ~ x') +
      geom_vline(data = dataGraph, aes(xintercept = MinSampleSize), linetype = "dashed", size = 1.1, col = "darkgoldenrod") +
      geom_hline(data = dataGraph, aes(yintercept = cross), linetype = "dashed", size = 1.1, col = "darkgoldenrod") +
      theme_bw() + xlab("Sample Size") + theme(text = element_text(size=12.5))
    print(p1)
  }
  
}

mean(MinSampleSizeSimple)
mean(MinSampleSize)
confint(lm(MinSampleSize~MinSampleSizeSimple),'MinSampleSizeSimple',level=0.95) #here we compare the slope and the intercept of the linear model which compares the two estimates of sample size 

### --- Validating optimization results with discrete events simulations ---###

#Here, we show how we run our discrete events simulations, how we measure different behavioral features, errors associated with those features, and how we execute piecewise continuous sampling

#clear everything but data
rm(list=setdiff(ls(), "data"))

#define functions
getmode = function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#set constants
SampleSizeVector = 575:825 #this captures the mean of the binomial sample size estimate
K = 9 #the number of tasks
T = 11041 #number of seconds the virtual and real ants are watched for
IValues = c(1, 2, 4, 8, 16, 32, 64, 128) #set of all interval numbers we iterate through
simulations = 250 #total number of discrete event simulations 
penalty = 0.2418923*0.06589 #This is the penatly for increasing I. The final sum of errors for each value of I would be error + penatly*log(I). See line 601 for how we got 0.2418923
IValues2 = log(IValues) #used for penalty 

#create bout length vectors for each task
dataAnts = as.matrix(data)
dataAnts = rbind(dataAnts, NA) #as the RLE function creates a vector, we want to ensure that the bout lengths across ants aren't concatenated with one another, so we put in some NA's
dataAnts = as.numeric(as.factor(dataAnts))
runLength = rle(dataAnts)
stateVector = as.numeric(na.omit(unique(runLength$values)))
x = 1:K

#initialize vectors
proportionOfStatesError = c()
boutLengthError = c()
sampleSize = c()
IMat = c()
degreeOfSpecializationSim = c()
meanBoutLengthSim = c()
numberOfTaskPerformedSim = c()
counterVector = 0
id = c()

for(i in 1:simulations){
  
  #discrete event simulation
  stateVector = sample(stateVector) #by shuffling the states, we ensure that different states will have different probabilities of occuring every simulation 
  d = runif(1)
  pr_x = (((1-d)^(x-1))*d)/(1-((1-d)^(K-1))) #equation 2
  pr_x = pr_x/sum(pr_x) #ensure probabilities sum to 1
  pr_x = cumsum(pr_x) 
  
  behaviorVec = c()
  state = stateVector[1]
  while(length(behaviorVec)<=T){ #here we choose states, and then draw from the corresponding bout length vector until the behavior vector is longer than T. We then truncate this vector so it is of length T 
    stateBouts = runLength$lengths[runLength$values==state]
    stateBouts = na.omit(stateBouts)
    behaviorVec = c(behaviorVec, rep(state, sample(stateBouts, 1)))
    state2 = sum(runif(1) > pr_x)+1
    while(state2==state){
      state2 = sum(runif(1) > pr_x)+1
    }
    state = state2
  }
  
  behaviorVec = behaviorVec[1:T]
  
  #summarize data for later derivations
  counts = table(behaviorVec)
  countDF = data.frame(names = names(counts), value = as.numeric(counts))
  StatesPresent = as.numeric(names(counts))
  if(length(StatesPresent)<9){ #we check for states tat are missing from the samples and raw datasets, ensuring that they have a value of '0'
    missingStates = x[!(x %in% StatesPresent)]
    missingRows = data.frame(names = missingStates, value = 0)
    countDF = rbind(countDF, missingRows)
    countDF = countDF[order(countDF$names), ]
  }
  #true proportion of states/bout lengths for the simulated ants
  proportionOfStates = countDF$value/T
  
  runLengthAnt = rle(behaviorVec)
  runLengthAnt = data.frame(state = runLengthAnt$values, bouts = runLengthAnt$lengths)
  
  boutLengths = aggregate(runLengthAnt$bouts~runLengthAnt$state, FUN = mean)
  boutLengths = data.frame(names = boutLengths$`runLengthAnt$state`, value = boutLengths$`runLengthAnt$bouts`)
  if(length(StatesPresent)<9){ #we also make sure that missing states have corresponding bout lengths of 0 
    boutLengths = rbind(boutLengths, missingRows)
    boutLengths = boutLengths[order(boutLengths$names), ]
  }
  boutLengths = boutLengths$value
  
  #piecewise continuous sampling
  
  for(j in 1:length(IValues)){
    
    I = IValues[j]
    
    for(n in SampleSizeVector){
      
      counterVector = counterVector + 1
      w = ceiling(n/I) #w is the length of each interval or window. it has to be a whole number, so we round up
      if(w == 0){
        w = 1
      }
      
      secondsVec = 1:(11041-w) #these are the seconds we can sample from, this prevents us from sampling over what is possible (we can't sample past T = 11,041)
      starts = sort(sample(x = secondsVec, size = I, replace = FALSE)) #the start of each interval
      #the next function checks to see if the intervals will overlap each other. If they do, then we select new starting points until there is no overlap
      differencesBool = diff(starts)<w
      counter = 0
      
      #this while loop finds non overlaping indexes that we can sample from later 
      while(any(differencesBool)){
        diffIndex = min(which(differencesBool == TRUE))
        invervalRand = sample(c(diffIndex, diffIndex+1), size=1)
        starts[invervalRand] = sample(x = secondsVec, size = 1)
        starts = sort(starts)
        differencesBool = diff(starts)<w
        counter = counter+1
        if(counter > 1000){
          differencesBool = FALSE
        }
      }
      
      #finally, the sample index allows us to sample from the behavioral dataset
      sampleIndex = c()
      for(m in 1:I){
        sampleIndex = c(sampleIndex, starts[m]:(starts[m]+w-1), NA)
      }
      
      sampleIndex = sampleIndex[1:(n+I)]
      sample = behaviorVec[sampleIndex]
      
      runLengthSample = rle(sample)
      runLengthSample = data.frame(state = runLengthSample$values, bouts = runLengthSample$lengths)
      
      #sample proportion of states/bout lengths for the simulated ants
      countsSample = table(sample)
      countsSampleDF = data.frame(names = names(countsSample), value = as.numeric(countsSample))
      StatesPresentSample = as.numeric(names(countsSample))
      if(length(StatesPresentSample)<9){
        missingStates = x[!(x %in% StatesPresentSample)]
        missingRows = data.frame(names = missingStates, value = 0)
        countsSampleDF = rbind(countsSampleDF, missingRows)
        countsSampleDF = countsSampleDF[order(countsSampleDF$names), ]
      }
      proportionOfStatesSample = countsSampleDF$value/sum(countsSample)
      
      boutLengthsSample = aggregate(runLengthSample$bouts~runLengthSample$state, FUN = mean)
      boutLengthsSample = data.frame(names = boutLengthsSample$`runLengthSample$state`, value = boutLengthsSample$`runLengthSample$bouts`)
      if(length(StatesPresentSample)<9){
        boutLengthsSample = rbind(boutLengthsSample, missingRows)
        boutLengthsSample = boutLengthsSample[order(boutLengthsSample$names), ]
      }
      boutLengthsSample = boutLengthsSample$value
      
      proportionOfStatesError[counterVector] = sum((proportionOfStatesSample-proportionOfStates)^2)
      boutLengthError[counterVector] = sum((boutLengthsSample-boutLengths)^2)
      sampleSize[counterVector] = n
      IMat[counterVector] = I
      
      #finally we measure characteristics of simulated ant such as degree of specialization 
      degreeOfSpecializationSim[counterVector] = sum(abs(proportionOfStates-1/K))/(2-(2/K))
      meanBoutLengthSim[counterVector] = mean(runLengthAnt$bouts)
      numberOfTaskPerformedSim[counterVector] = length(StatesPresent)
      id[counterVector] = i
      
      print(counterVector/(length(SampleSizeVector)*length(IValues)*simulations))
      
    }
    
  }
  
}

dataSimulated = data.frame(id = id, proportionOfStatesError = proportionOfStatesError, boutLengthError = boutLengthError, sampleSize = sampleSize, I = IMat, degreeOfSpecializationSim = degreeOfSpecializationSim, meanBoutLengthSim = meanBoutLengthSim, numberOfTaskPerformedSim = numberOfTaskPerformedSim)

### next, we repeat the above measurements for real ants

#convert behavior to a numeric code
dataAnts = as.numeric(as.factor(as.matrix(data)))
dataAnts = matrix(dataAnts, nrow = 11041, ncol = 9)

proportionOfStatesError = c()
boutLengthError = c()
sampleSize = c()
IMat = c()
degreeOfSpecialization = c()
meanBoutLength = c()
numberOfTaskPerformed = c()
counterVector = 0
id = c()

for(i in 1:9){
  
  behaviorVec = dataAnts[,i]
  
  counts = table(behaviorVec)
  countDF = data.frame(names = names(counts), value = as.numeric(counts))
  StatesPresent = as.numeric(names(counts))
  if(length(StatesPresent)<9){
    missingStates = x[!(x %in% StatesPresent)]
    missingRows = data.frame(names = missingStates, value = 0)
    countDF = rbind(countDF, missingRows)
    countDF = countDF[order(countDF$names), ]
  }
  proportionOfStates = countDF$value/T
  
  runLengthAnt = rle(behaviorVec)
  runLengthAnt = data.frame(state = runLengthAnt$values, bouts = runLengthAnt$lengths)
  
  boutLengths = aggregate(runLengthAnt$bouts~runLengthAnt$state, FUN = mean)
  boutLengths = data.frame(names = boutLengths$`runLengthAnt$state`, value = boutLengths$`runLengthAnt$bouts`)
  if(length(StatesPresent)<9){
    boutLengths = rbind(boutLengths, missingRows)
    boutLengths = boutLengths[order(boutLengths$names), ]
  }
  boutLengths = boutLengths$value
  
  for(j in 1:length(IValues)){
    
    I = IValues[j]
    
    for(n in SampleSizeVector){
      
      counterVector = counterVector + 1
      w = ceiling(n/I)
      if(w == 0){
        w = 1
      }
      
      secondsVec = 1:(11041-w)
      starts = sort(sample(x = secondsVec, size = I, replace = FALSE))
      differencesBool = diff(starts)<w
      counter = 0
      
      while(any(differencesBool)){
        diffIndex = min(which(differencesBool == TRUE))
        invervalRand = sample(c(diffIndex, diffIndex+1), size=1)
        starts[invervalRand] = sample(x = secondsVec, size = 1)
        starts = sort(starts)
        differencesBool = diff(starts)<w
        counter = counter+1
        if(counter > 1000){
          differencesBool = FALSE
        }
      }
      
      sampleIndex = c()
      for(m in 1:I){
        sampleIndex = c(sampleIndex, starts[m]:(starts[m]+w-1), NA)
      }
      
      sampleIndex = sampleIndex[1:(n+I)]
      sample = behaviorVec[sampleIndex]
      
      runLengthSample = rle(sample)
      runLengthSample = data.frame(state = runLengthSample$values, bouts = runLengthSample$lengths)
      
      countsSample = table(sample)
      countsSampleDF = data.frame(names = names(countsSample), value = as.numeric(countsSample))
      StatesPresentSample = as.numeric(names(countsSample))
      if(length(StatesPresentSample)<9){
        missingStates = x[!(x %in% StatesPresentSample)]
        missingRows = data.frame(names = missingStates, value = 0)
        countsSampleDF = rbind(countsSampleDF, missingRows)
        countsSampleDF = countsSampleDF[order(countsSampleDF$names), ]
      }
      proportionOfStatesSample = countsSampleDF$value/sum(countsSample)
      
      boutLengthsSample = aggregate(runLengthSample$bouts~runLengthSample$state, FUN = mean)
      boutLengthsSample = data.frame(names = boutLengthsSample$`runLengthSample$state`, value = boutLengthsSample$`runLengthSample$bouts`)
      if(length(StatesPresentSample)<9){
        boutLengthsSample = rbind(boutLengthsSample, missingRows)
        boutLengthsSample = boutLengthsSample[order(boutLengthsSample$names), ]
      }
      boutLengthsSample = boutLengthsSample$value
      
      proportionOfStatesError[counterVector] = sum((proportionOfStates-proportionOfStatesSample)^2)
      boutLengthError[counterVector] = sum((boutLengthsSample-boutLengths)^2)
      sampleSize[counterVector] = n
      IMat[counterVector] = I
      
      degreeOfSpecialization[counterVector] = sum(abs(proportionOfStates-1/K))/(2-(2/K))
      meanBoutLength[counterVector] = mean(runLengthAnt$bouts)
      numberOfTaskPerformed[counterVector] = length(StatesPresent)
      
      id[counterVector] = i
      
      print(counterVector/(length(SampleSizeVector)*length(IValues)*9))
      
    }
    
  }
  
}

dataReal = data.frame(id = id, proportionOfStatesError = proportionOfStatesError, boutLengthError = boutLengthError, sampleSize = sampleSize, I = IMat, degreeOfSpecialization = degreeOfSpecialization, meanBoutLength = meanBoutLength, numberOfTaskPerformed = numberOfTaskPerformed)

### finding the optimal I value for each simulated ant ###

IMat = c()
errorSum = c()
DOSVec = c()
DurationVec = c()
TaskNumberVec = c()

#find realistic values of each behavioral metric from real ant data 
DOSRange = range(dataReal$degreeOfSpecialization)
ADRange = range(dataReal$meanBoutLength)
StatesRange = range(dataReal$numberOfTaskPerformed)

for(i in 1:simulations){
  
  #check if simulated ant looks like real ant 
  dataTemp = subset(dataSimulated, id == i & degreeOfSpecializationSim > DOSRange[1] & degreeOfSpecializationSim < DOSRange[2] & meanBoutLengthSim > ADRange[1] & meanBoutLengthSim < ADRange[2] & numberOfTaskPerformedSim >= StatesRange[1] & numberOfTaskPerformedSim <= StatesRange[2])
  
  #if the simulated ant doesn't look like the real ants, delete it 
  if(nrow(dataTemp) == 0){
    IMat[i] = NA
    errorSum[i] = NA
    DOSVec[i] = NA
    DurationVec[i] = NA
    TaskNumberVec[i] = NA
  } else {
    
    #standardize the errors by making sure that duration error (bout length error) is on the same scale as proportion of states error 
    coeff = mean(dataTemp$boutLengthError)/mean(dataTemp$proportionOfStatesError)
    dataTemp$normDurationError = dataTemp$boutLengthError/coeff
    
    #find average error for each value of I across sample sizes 
    aggregatePropError = aggregate(proportionOfStatesError~I, data = dataTemp, mean)
    aggregateDurationError = aggregate(normDurationError~I, data = dataTemp, mean)
    aggregateErrorSum = aggregatePropError$proportionOfStatesError + aggregateDurationError$normDurationError+IValues2*penalty #add penalty term to error sum 
    errorSum[i] = min(aggregateErrorSum)
    IIndex = which.min(aggregateErrorSum) #find the value of I that minimizes the sum 
    IMat[i]= IValues[IIndex]
    
    aggregateDOS = aggregate(degreeOfSpecializationSim~I, data = dataTemp, mean)
    aggregateDuration = aggregate(meanBoutLengthSim~I, data = dataTemp, mean)
    aggregateTaskNumber = aggregate(numberOfTaskPerformedSim~I, data = dataTemp, mean)
    DOSVec[i] = aggregateDOS$degreeOfSpecializationSim[IIndex] #track characteristics of sample that could have influenced the optimal value of I
    DurationVec[i] = aggregateDuration$meanBoutLengthSim[IIndex]
    TaskNumberVec[i] = aggregateTaskNumber$numberOfTaskPerformedSim[IIndex]
    
    print(i/simulations)
  }
  
}

dataOptimalSim = data.frame(DOS = DOSVec, AverageTaskDuration = DurationVec, TaskNumber = TaskNumberVec, I = IMat)
dataOptimalSim = na.omit(dataOptimalSim)
dataOptimalSim$logI = log(dataOptimalSim$I)

#count how often each value of I is optimal
table(dataOptimalSim$I)

### finding the optimal I value for each real ant ###

IMat = c()
errorSum = c()
DOSVec = c()
DurationVec = c()
TaskNumberVec = c()

for(i in 1:9){
  
  dataTemp = subset(dataReal, id == i)
  
  coeff = mean(dataTemp$boutLengthError)/mean(dataTemp$proportionOfStatesError)
  dataTemp$normDurationError = dataTemp$boutLengthError/coeff
  
  aggregatePropError = aggregate(proportionOfStatesError~I, data = dataTemp, mean)
  aggregateDurationError = aggregate(normDurationError~I, data = dataTemp, mean)
  aggregateErrorSum = aggregatePropError$proportionOfStatesError + aggregateDurationError$normDurationError+IValues2*penalty
  errorSum[i] = min(aggregateErrorSum)
  IIndex = which.min(aggregateErrorSum)
  IMat[i]= IValues[IIndex]
  
  aggregateDOS = aggregate(degreeOfSpecialization~I, data = dataTemp, mean)
  aggregateDuration = aggregate(meanBoutLength~I, data = dataTemp, mean)
  aggregateTaskNumber = aggregate(numberOfTaskPerformed~I, data = dataTemp, mean)
  DOSVec[i] = aggregateDOS$degreeOfSpecialization[IIndex]
  DurationVec[i] = aggregateDuration$meanBoutLength[IIndex]
  TaskNumberVec[i] = aggregateTaskNumber$numberOfTaskPerformed[IIndex]
  
  print(i/9)
  
}

dataOptimalReal = data.frame(DOS = DOSVec, AverageTaskDuration = DurationVec, TaskNumber = TaskNumberVec, I = IMat)
dataOptimalReal$logI = log(dataOptimalReal$I)

cor1 = cor.test(dataOptimalReal$I, dataOptimalReal$DOS, method= 'spearman')
cor2 = cor.test(dataOptimalReal$I, dataOptimalReal$AverageTaskDuration, method= 'spearman')
cor3 = cor.test(dataOptimalReal$I, dataOptimalReal$TaskNumber, method= 'spearman')

p.adjust(c(cor1$p.value, cor2$p.value,cor3$p.value), method = "holm")

cor1 = cor.test(dataOptimalSim$I, dataOptimalSim$DOS, method= 'spearman')
cor2 = cor.test(dataOptimalSim$I, dataOptimalSim$AverageTaskDuration, method= 'spearman')
cor3 = cor.test(dataOptimalSim$I, dataOptimalSim$TaskNumber, method= 'spearman')

p.adjust(c(cor1$p.value, cor2$p.value,cor3$p.value), method = "holm")

#fig 5
p1 = ggplot(dataOptimalSim, aes(x = DOS, y = I)) + geom_point() + theme_bw() + theme(text = element_text(size=12.5)) + xlab("Degree of Specialization") + ylab("# Intervals (I)") + ggtitle("A") + ylim(0, 16) + geom_point(data = dataOptimalReal, aes(x = DOS, y = I), col = "gold", size = 2)

p2 = ggplot(dataOptimalSim, aes(x = AverageTaskDuration, y = I)) + geom_point() + theme_bw() + theme(text = element_text(size=12.5)) + xlab("Mean Bout Length (seconds)") + ylab("# Intervals (I)") + ggtitle("B")  + ylim(0, 16) + geom_point(data = dataOptimalReal, aes(x = AverageTaskDuration, y = I), col = "gold", size = 2) + geom_smooth(data = dataOptimalReal, aes(x = AverageTaskDuration, y = I), method='lm', formula= y~x, color = "gold", se = F)

p3 = ggplot(dataOptimalSim, aes(x = TaskNumber, y = I)) + geom_point() + theme_bw() + theme(text = element_text(size=12.5)) + xlab("Number of Tasks Performed") + ylab("# Intervals (I)")+ geom_jitter(width = 0.25) + ggtitle("C")  + ylim(0, 16)+ geom_point(data = dataOptimalReal, aes(x = TaskNumber, y = I), col = "gold", size = 2)

ggarrange(p1, p2, p3, nrow = 3, ncol = 1)

getmode(dataOptimalSim$I)
getmode(dataOptimalReal$I)

### Visualize Tradeoff Over I - fig 4 ###

#first we will visualize this tradeoff for simulated ants

IMat = array(numeric(), c(simulations, length(IValues)))
PropStatesError = array(numeric(), c(simulations, length(IValues)))
DurationError = array(numeric(), c(simulations, length(IValues)))

for(i in 1:simulations){
  for(j in 1:length(IValues)){
    
    IMat[i,j] = IValues[j]
    
    dataTemp = subset(dataSimulated, I == IValues[j] & id == i & degreeOfSpecializationSim > DOSRange[1] & degreeOfSpecializationSim < DOSRange[2] & meanBoutLengthSim > ADRange[1] & meanBoutLengthSim < ADRange[2] & numberOfTaskPerformedSim >= StatesRange[1] & numberOfTaskPerformedSim <= StatesRange[2])
    
    #average error for each I value within each simulation but across sample sizes 
    PropStatesError[i,j] = mean(dataTemp$proportionOfStatesError)
    DurationError[i,j] = mean(dataTemp$boutLengthError)
    
  }
}

visualDataProp = data.frame(I = c(IMat), PropError = c(PropStatesError))
visualDataProp1 = aggregate(PropError~I, data = visualDataProp, mean)
visualDataProp2 = aggregate(PropError~I, data = visualDataProp, FUN= function(x) sqrt(var(x)/length(x)))
visualDataProp1$StandardErrorProp = visualDataProp2$PropError

visualDataDur = data.frame(I = c(IMat), DurationError = c(DurationError))
visualDataDur1 = aggregate(DurationError~I, data = visualDataDur, mean)
visualDataDur2 = aggregate(DurationError~I, data = visualDataDur, FUN= function(x) sqrt(var(x)/length(x)))
visualDataDur1$StandardErrorDur = visualDataDur2$DurationError

drops = c("I")
visualDataDur1 = visualDataDur1[ , !(names(visualDataDur1) %in% drops)]
visualDataFull = cbind(visualDataProp1, visualDataDur1)

#we use the following coefficient to ensure that each type of error is on the same scale
coeff = mean(visualDataFull$DurationError)/mean(visualDataFull$PropError)
visualDataFull$Sum = visualDataFull$PropError+(visualDataFull$DurationError/coeff) + IValues2*penalty #The penalty is added to the log of I 

p2 = ggplot(visualDataFull, aes(x=as.factor(I), y=PropError)) +
  geom_point(size = 2, col = "5ab4ac") + geom_errorbar(aes(ymin=PropError-StandardErrorProp, ymax=PropError+StandardErrorProp), width=.2, col = "5ab4ac") + geom_path(aes(y = PropError), group = 1, size = 1, col = "5ab4ac") +
  geom_point(aes(x=as.factor(I), y=DurationError/coeff), size = 2, col = "#d8b365") + geom_errorbar(aes(ymin=(DurationError-StandardErrorDur)/coeff, ymax=(DurationError+StandardErrorDur)/coeff), width=.2, col = "#d8b365") + geom_path(aes(y=DurationError/coeff), group = 1, size = 1, col = "#d8b365") +
  scale_y_continuous(name = "Proporion of States Error", sec.axis = sec_axis(~.*coeff, name="Mean Bout Length Error")) + xlab("Number of Intervals (I)")+ theme_bw() + theme(text = element_text(size=12.5)) + ggtitle("B") + geom_point(aes(x=as.factor(I), y=Sum), size = 2, col = "#999999") + geom_path(aes(y=Sum), group = 1, size = 1, col = "#999999")

### now we visualize error across I for real ants ### 

IMat = array(numeric(), c(9, length(IValues)))
PropStatesError = array(numeric(), c(9, length(IValues)))
DurationError = array(numeric(), c(9, length(IValues)))

for(i in 1:9){
  for(j in 1:length(IValues)){
    
    IMat[i,j] = IValues[j]
    
    dataTemp = subset(dataReal, I == IValues[j] & id == i)
    
    PropStatesError[i,j] = mean(dataTemp$proportionOfStatesError)
    DurationError[i,j] = mean(dataTemp$boutLengthError)
    
  }
}

visualDataProp = data.frame(I = c(IMat), PropError = c(PropStatesError))
visualDataProp1 = aggregate(PropError~I, data = visualDataProp, mean)
visualDataProp2 = aggregate(PropError~I, data = visualDataProp, FUN= function(x) sqrt(var(x)/length(x)))
visualDataProp1$StandardErrorProp = visualDataProp2$PropError

visualDataDur = data.frame(I = c(IMat), DurationError = c(DurationError))
visualDataDur1 = aggregate(DurationError~I, data = visualDataDur, mean)
visualDataDur2 = aggregate(DurationError~I, data = visualDataDur, FUN= function(x) sqrt(var(x)/length(x)))
visualDataDur1$StandardErrorDur = visualDataDur2$DurationError

drops = c("I")
visualDataDur1 = visualDataDur1[ , !(names(visualDataDur1) %in% drops)]
visualDataFull = cbind(visualDataProp1, visualDataDur1)

coeff = mean(visualDataFull$DurationError)/mean(visualDataFull$PropError)

#range in error for real ants, used to compute penalty: 
visualDataFull$Sum = visualDataFull$PropError+(visualDataFull$DurationError/coeff) 
range(visualDataFull$Sum)[2]-range(visualDataFull$Sum)[1]
visualDataFull$Sum = visualDataFull$Sum + IValues2*penalty

p1 = ggplot(visualDataFull, aes(x=as.factor(I), y=PropError)) +
  geom_point(size = 2, col = "5ab4ac") + geom_errorbar(aes(ymin=PropError-StandardErrorProp, ymax=PropError+StandardErrorProp), width=.2, col = "5ab4ac") + geom_path(aes(y = PropError), group = 1, size = 1, col = "5ab4ac") +
  geom_point(aes(x=as.factor(I), y=DurationError/coeff), size = 2, col = "#d8b365") + geom_errorbar(aes(ymin=(DurationError-StandardErrorDur)/coeff, ymax=(DurationError+StandardErrorDur)/coeff), width=.2, col = "#d8b365") + geom_path(aes(y=DurationError/coeff), group = 1, size = 1, col = "#d8b365") +
  scale_y_continuous(name = "Proporion of States Error", sec.axis = sec_axis(~.*coeff, name="Mean Bout Length Error")) + xlab("Number of Intervals (I)")+ theme_bw() + theme(text = element_text(size=12.5)) + ggtitle("A") + geom_point(aes(x=as.factor(I), y=Sum), size = 2, col = "#999999") + geom_path(aes(y=Sum), group = 1, size = 1, col = "#999999")

ggarrange(p1, p2, nrow = 1, ncol = 2)

### now evaluate error sums over I for activity to validate results from the task-based dataset. We go through the same steps as before, where we use piecewise continuous sampling to estimate to estimate error 

dataActivity = read.csv("rawDataActivity.csv")

#show that there are strong differences in both datasets. We first collect the behavioral attributes of the ants from the task-based dataset 

taskDOS = dataOptimalReal$DOS 
taskDuration = dataOptimalReal$AverageTaskDuration
taskNumber = dataOptimalReal$TaskNumber

#and next we will estiamte these attributes for these ants in the activity dataset, but we will do this in the context of the piecewise continuous sampling function as we did before 
dataAnts = as.numeric(as.factor(as.matrix(dataActivity)))
dataAnts = matrix(dataAnts, nrow = 11041, ncol = 9)

proportionOfStatesError = c()
boutLengthError = c()
sampleSize = c()
IMat = c()
degreeOfSpecialization = c()
meanBoutLength = c()
numberOfTaskPerformed = c()
counterVector = 0
id = c()

activityDOS = c()
activityDuration = c()
activityNumber = c() 

for(i in 1:9){
  
  behaviorVec = dataAnts[,i]
  
  counts = table(behaviorVec)
  countDF = data.frame(names = names(counts), value = as.numeric(counts))
  StatesPresent = as.numeric(names(counts))
  proportionOfStates = countDF$value/T
  
  runLengthAnt = rle(behaviorVec)
  runLengthAnt = data.frame(state = runLengthAnt$values, bouts = runLengthAnt$lengths)
  
  boutLengths = aggregate(runLengthAnt$bouts~runLengthAnt$state, FUN = mean)
  boutLengths = data.frame(names = boutLengths$`runLengthAnt$state`, value = boutLengths$`runLengthAnt$bouts`)
  boutLengths = boutLengths$value
  
  for(j in 1:length(IValues)){
    
    I = IValues[j]
    
    for(n in SampleSizeVector){
      
      counterVector = counterVector + 1
      w = ceiling(n/I)
      if(w == 0){
        w = 1
      }
      
      secondsVec = 1:(11041-w)
      starts = sort(sample(x = secondsVec, size = I, replace = FALSE))
      differencesBool = diff(starts)<w
      counter = 0
      
      while(any(differencesBool)){
        diffIndex = min(which(differencesBool == TRUE))
        invervalRand = sample(c(diffIndex, diffIndex+1), size=1)
        starts[invervalRand] = sample(x = secondsVec, size = 1)
        starts = sort(starts)
        differencesBool = diff(starts)<w
        counter = counter+1
        if(counter > 1000){
          differencesBool = FALSE
        }
      }
      
      sampleIndex = c()
      for(m in 1:I){
        sampleIndex = c(sampleIndex, starts[m]:(starts[m]+w-1), NA)
      }
      
      sampleIndex = sampleIndex[1:(n+I)]
      sample = behaviorVec[sampleIndex]
      
      runLengthSample = rle(sample)
      runLengthSample = data.frame(state = runLengthSample$values, bouts = runLengthSample$lengths)
      
      countsSample = table(sample)
      countsSampleDF = data.frame(names = names(countsSample), value = as.numeric(countsSample))
      StatesPresentSample = as.numeric(names(countsSample))
      proportionOfStatesSample = countsSampleDF$value/sum(countsSample)
      
      boutLengthsSample = aggregate(runLengthSample$bouts~runLengthSample$state, FUN = mean)
      boutLengthsSample = data.frame(names = boutLengthsSample$`runLengthSample$state`, value = boutLengthsSample$`runLengthSample$bouts`)
      boutLengthsSample = boutLengthsSample$value
      
      proportionOfStatesError[counterVector] = sum((proportionOfStates-proportionOfStatesSample)^2)
      boutLengthError[counterVector] = sum((boutLengthsSample-boutLengths)^2)
      sampleSize[counterVector] = n
      IMat[counterVector] = I
      
      degreeOfSpecialization[counterVector] = sum(abs(proportionOfStates-1/3))/(2-(2/3))
      meanBoutLength[counterVector] = mean(runLengthAnt$bouts)
      numberOfTaskPerformed[counterVector] = length(StatesPresent)
      
      id[counterVector] = i
      
      print(counterVector/(length(SampleSizeVector)*length(IValues)*9))
      
    }
    
  }
  
  activityDOS[i] =  sum(abs(proportionOfStates-1/3))/(2-(2/3))
  activityDuration[i] = mean(runLengthAnt$bouts)
  activityNumber[i] = length(StatesPresent)
  
}

dataReal = data.frame(id = id, proportionOfStatesError = proportionOfStatesError, boutLengthError = boutLengthError, sampleSize = sampleSize, I = IMat, degreeOfSpecialization = degreeOfSpecialization, meanBoutLength = meanBoutLength, numberOfTaskPerformed = numberOfTaskPerformed)

#show differences in task and activity data by performing paired t-tests between activity measurements of behavioral attributes as well as the task-based measuremements 
t1 = t.test(activityDOS, taskDOS, paired = TRUE)
t2 = t.test(activityDuration, taskDuration, paired = TRUE)
t3 = t.test(activityNumber, taskNumber, paired = TRUE)

#adjust p-values to control for multiple comparisons 
p.adjust(c(t1$p.value, t2$p.value, t3$p.value), method = 'holm')

dataGraph = data.frame(DOS = c(activityDOS, taskDOS), Duration = c(activityDuration, taskDuration), stateNumber = c(activityNumber, taskNumber), Type = c(rep(1, 9), rep(2, 9)), ID = 1:9)

#Fig S1
p1 = ggplot(dataGraph) +
  geom_boxplot(aes(x = as.numeric(Type), y = DOS, group = Type))+
  geom_point(aes(x = as.numeric(Type), y = DOS)) +
  geom_line(aes(x  = as.numeric(Type), y = DOS, group = ID)) +
  scale_x_continuous(breaks = c(1,2), labels = c("Activity", "Task"))+
  xlab("") + ylab("Degree of Specialization") + theme_bw() + theme(text = element_text(size=12.5)) + ggtitle("A")

p2 = ggplot(dataGraph) +
  geom_boxplot(aes(x = as.numeric(Type), y = Duration, group = Type))+
  geom_point(aes(x = as.numeric(Type), y = Duration)) +
  geom_line(aes(x  = as.numeric(Type), y = Duration, group = ID)) +
  scale_x_continuous(breaks = c(1,2), labels = c("Activity", "Task"))+
  xlab("") + ylab("Mean Bout Length") + theme_bw() + theme(text = element_text(size=12.5)) + ggtitle("B")

p3 = ggplot(dataGraph) +
  geom_boxplot(aes(x = as.numeric(Type), y = stateNumber, group = Type))+
  geom_point(aes(x = as.numeric(Type), y = stateNumber)) +
  geom_line(aes(x  = as.numeric(Type), y = stateNumber, group = ID)) +
  scale_x_continuous(breaks = c(1,2), labels = c("Activity", "Task"))+
  xlab("") + ylab("Number of States") + theme_bw() + theme(text = element_text(size=12.5)) + ggtitle("C")

ggarrange(p1, p2, p3, nrow = 1)

#now visualize error over I for activity in real ants 

IMat = array(numeric(), c(9, length(IValues)))
PropStatesError = array(numeric(), c(9, length(IValues)))
DurationError = array(numeric(), c(9, length(IValues)))

for(i in 1:9){
  for(j in 1:length(IValues)){
    
    IMat[i,j] = IValues[j]
    
    dataTemp = subset(dataReal, I == IValues[j] & id == i)
    
    PropStatesError[i,j] = mean(dataTemp$proportionOfStatesError)
    DurationError[i,j] = mean(dataTemp$boutLengthError)
    
  }
}

visualDataProp = data.frame(I = c(IMat), PropError = c(PropStatesError))
visualDataProp1 = aggregate(PropError~I, data = visualDataProp, mean)
visualDataProp2 = aggregate(PropError~I, data = visualDataProp, FUN= function(x) sqrt(var(x)/length(x)))
visualDataProp1$StandardErrorProp = visualDataProp2$PropError

visualDataDur = data.frame(I = c(IMat), DurationError = c(DurationError))
visualDataDur1 = aggregate(DurationError~I, data = visualDataDur, mean)
visualDataDur2 = aggregate(DurationError~I, data = visualDataDur, FUN= function(x) sqrt(var(x)/length(x)))
visualDataDur1$StandardErrorDur = visualDataDur2$DurationError

drops = c("I")
visualDataDur1 = visualDataDur1[ , !(names(visualDataDur1) %in% drops)]
visualDataFull = cbind(visualDataProp1, visualDataDur1)

coeff = mean(visualDataFull$DurationError)/mean(visualDataFull$PropError)

visualDataFull$Sum = visualDataFull$PropError+(visualDataFull$DurationError/coeff) 
visualDataFull$Sum = visualDataFull$Sum + IValues2*penalty

#Fig S2
ggplot(visualDataFull, aes(x=as.factor(I), y=PropError)) +
  geom_point(size = 2, col = "5ab4ac") + geom_errorbar(aes(ymin=PropError-StandardErrorProp, ymax=PropError+StandardErrorProp), width=.2, col = "5ab4ac") + geom_path(aes(y = PropError), group = 1, size = 1, col = "5ab4ac") +
  geom_point(aes(x=as.factor(I), y=DurationError/coeff), size = 2, col = "#d8b365") + geom_errorbar(aes(ymin=(DurationError-StandardErrorDur)/coeff, ymax=(DurationError+StandardErrorDur)/coeff), width=.2, col = "#d8b365") + geom_path(aes(y=DurationError/coeff), group = 1, size = 1, col = "#d8b365") +
  scale_y_continuous(name = "Proporion of States Error", sec.axis = sec_axis(~.*coeff, name="Mean Bout Length Error")) + xlab("Number of Intervals (I)")+ theme_bw() + theme(text = element_text(size=12.5)) + geom_point(aes(x=as.factor(I), y=Sum), size = 2, col = "#999999") + geom_path(aes(y=Sum), group = 1, size = 1, col = "#999999")

### Measuring penalty for increasing the interval number I ###

rm(list = ls())

#Here we load data 
dataAgreement = read.csv("agreementData.csv")
mdl = lme(percentAgreement ~ log(I),random=~1|Colony,data=dataAgreement)
summary(mdl)
anova(mdl)

#figure s3
ggplot(dataAgreement, aes(x = log(I), y = percentAgreement)) + geom_point(size = 2) + theme_bw() + xlab("log(I)") + ylab("% Agreement") + theme(text = element_text(size=12.5)) + geom_smooth(method = 'lm', se = F)

### Next we fit different types of distributions for the bout lengths of different tasks ### 

rm(list = ls())

data = read.csv("rawDataTask.csv")
data = c(data$BBB, data$GPW, data$GWB, data$WWR, data$YBR, data$YGR, data$YWR, data$YYR, data$YYW)
boutLengths = rle(data)
states = unique(data)

#we iterate through each state, randomly sample (with replacement) 50 bouts from that tasks's list of bout lengths, and then we fit different distributions and measure which distribution fits the best and whether this fit is significant 
for(i in 1:length(states)){
  dataTemp = boutLengths$lengths[boutLengths$values==states[i]]
  if(length(dataTemp>50)){
    dataTemp = sample(dataTemp, 50)
  }
  fitPois = fitdist(dataTemp, "pois")
  fitGeom = fitdist(dataTemp, "geom")
  fitNBin = fitdist(dataTemp, "nbinom")
  gofstat(list(fitPois, fitGeom, fitNBin), fitnames = c("pois","geom","nbinom"))
  m_pl = displ$new(dataTemp)
  bs = bootstrap_p(m_pl, threads = 4)
  bs$p
}

boutLengths = data.frame(BoutLength = boutLengths$lengths, Task = boutLengths$values)
distType = c()
for(i in 1:nrow(boutLengths)){
  if(boutLengths$Task[i]=="Trash Maintenance"){
    distType[i] = "Negative Binomial"
  } else if(boutLengths$Task[i]=="Allogrooming" | boutLengths$Task[i] == "Exploring" | boutLengths$Task[i] == "Nest Maintenance"){
    distType[i] = "Geometric"
  } else {
    distType[i] = "Power Law"
  }
}

boutLengths$distType = distType

#figure S4
ggplot(boutLengths, aes(x = BoutLength, fill = distType)) + geom_histogram() + facet_wrap(~Task, scales = "free") + theme_bw() + xlab("Bout Length") + ylab("Count") + theme(text = element_text(size=12.5)) + labs(fill = "Distribution Type") 

### Expected Shape of Bout Length Error Curve ### 

#first, we create a vector from which we later sample. First, we sequentially sample from a geometric distribution, first sampling from a distribution whose mean is the same as that of brood care, and then we sample from another distribution whose mean bout length is the average of all other tasks.
rm(list = ls())
n = 660
bout1 = 82
bout2 = 36
iterations = 1000
sims = 100
behaviorVecGeom = c()
for(i in 1:iterations){
  behaviorVecGeom = c(behaviorVecGeom, rep(1, rgeom(1, p = 1/bout1)))
  behaviorVecGeom = c(behaviorVecGeom, rep(2, rgeom(1, p = 1/bout2)))
}

sampleMetric = rle(behaviorVecGeom)
sampleMetric = data.frame(lengths = sampleMetric$lengths, values = sampleMetric$values)
sampleMetric1 = subset(sampleMetric, values==1)
sampleMetric2 = subset(sampleMetric, values==2)

#we repeat this process, but instead we sample from a power law distribution instead using the poweRlaw package 
m_pl1 = displ$new(sampleMetric1$lengths) #fit separately to bout1 and 2
est = estimate_xmin(m_pl1)
m_pl1$setXmin(est)

m_pl2 = displ$new(sampleMetric2$lengths) #fit separately to bout1 and 2
est = estimate_xmin(m_pl2)
m_pl2$setXmin(est)

behaviorVecPower = c()

for(i in 1:iterations){
  behaviorVecPower = c(behaviorVecPower, rep(1, dist_rand(m_pl1, n = 1)))
  behaviorVecPower = c(behaviorVecPower, rep(2, dist_rand(m_pl2, n = 1)))
}

#now we will sample from these behavioral vectors using piecewise continuous sampling, although here we allow for potential overlap so that simulation-based estimates of error could align with closed form solution. First we sample from the geometric behavioral vector
IVec = c(1, 2, 4, 8, 16, 32, 64, 128)
wVec = round(n/IVec)
meanOverall = c()
meanProp1 = c()
simulations = 250
for(k in 1:length(wVec)){
  windowSize = wVec[k]
  I = IVec[k]
  meanEstimate = matrix(, simulations, I)
  for(j in 1:simulations){
    for(i in 1:I){
      startingPoint = sample(1:length(behaviorVecGeom), 1, replace=T)
      sample = behaviorVecGeom[startingPoint:(windowSize+startingPoint-1)]
      sampleMetric = rle(sample)
      sampleMetric = data.frame(lengths = sampleMetric$lengths, values = sampleMetric$values)
      sampleMetric = subset(sampleMetric, values==1)
      meanEstimate[j, i] = mean(sampleMetric$lengths)
      if(is.nan(meanEstimate[j, i])){
        meanEstimate[j, i] = 0
      }
    }
  }
  
  meanOverall[k] = mean(meanEstimate)
  
}

boutLengthErrorGeom = (meanOverall-mean(sampleMetric1$lengths))^2
boutLengthDiffGeom = meanOverall-mean(sampleMetric1$lengths)

#now we sample from the power law behavioral vector 

meanOverall = c()
meanProp1 = c()
for(k in 1:length(wVec)){
  windowSize = wVec[k]
  I = IVec[k]
  meanEstimate = matrix(, simulations, I)
  numberOfEvents = matrix(, simulations, I)
  for(j in 1:simulations){
    for(i in 1:I){
      startingPoint = sample(1:length(behaviorVecPower), 1, replace=T)
      sample = behaviorVecPower[startingPoint:(windowSize+startingPoint-1)]
      sampleMetric = rle(sample)
      sampleMetric = data.frame(lengths = sampleMetric$lengths, values = sampleMetric$values)
      sampleMetric = subset(sampleMetric, values==1)
      meanEstimate[j, i] = mean(sampleMetric$lengths)
      if(is.nan(meanEstimate[j, i])){
        meanEstimate[j, i] = 0
      }
    }
  }
  
  meanOverall[k] = mean(meanEstimate)
  
}

boutLengthErrorPower = (meanOverall-mean(sampleMetric1$lengths))^2
boutLengthDiffPower = meanOverall-mean(sampleMetric1$lengths)

dataGraph1 = data.frame(I = IVec, Error = boutLengthErrorGeom)
dataGraph2 = data.frame(I = IVec, Error = boutLengthErrorPower)
dataGraph3 = data.frame(I = IVec, Diff = boutLengthDiffGeom)
dataGraph4 = data.frame(I = IVec, Diff = boutLengthDiffPower)

#figure 5

p1 = ggplot(dataGraph1, aes(x = I, y = Error)) + geom_point(size = 2.5) + geom_line(size = 1) + theme_bw() + xlab("# Intervals (I)") + ylab("Bout Length Error \nfor Task 1 (E)") + theme(text = element_text(size=12.5)) + ggtitle("A")
p2 = ggplot(dataGraph2, aes(x = I, y = Error)) + geom_point(size = 2.5) + geom_line(size = 1) + theme_bw() + xlab("# Intervals (I)") + ylab("Bout Length Error \nfor Task 1 (E)") + theme(text = element_text(size=12.5)) + ggtitle("B")
p3 = ggplot(dataGraph3, aes(x = I, y = Diff)) + geom_point(size = 2.5) + geom_line(size = 1) + theme_bw() + xlab("# Intervals (I)") + ylab("Estimated Bout Mean - \n True Bout Mean") + theme(text = element_text(size=12.5)) + ggtitle("C") + geom_hline(aes(yintercept = 0), color = "gray", linetype=2, size = 1)
p4 = ggplot(dataGraph4, aes(x = I, y = Diff)) + geom_point(size = 2.5) + geom_line(size = 1) + theme_bw() + xlab("# Intervals (I)") + ylab("Estimated Bout Mean - \n True Bout Mean") + theme(text = element_text(size=12.5)) + ggtitle("D") + geom_hline(aes(yintercept = 0), color = "gray", linetype=2, size = 1)

ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)

### fig s6: we first plot the closed form solution for error given these average bout lengths
lone = bout1
ltwo = bout2
IVec = seq(0, 128, by=.01)
wVec = n/IVec
E = c()
for(i in 1:length(wVec)){
  if(ltwo > (wVec[i]-1)){
    E[i] = ((-lone^2-lone*ltwo+lone*wVec[i])^2)/((lone+ltwo)^2)
  } else {
    E[i] = ((-lone^2+lone)^2)/((wVec[i]-1+lone)^2)
  }
}

dataGraph = data.frame(I = IVec, Error = E)
dataGraph2 = data.frame(I = c(1, 2, 4, 8, 16, 32, 64, 128), E = boutLengthErrorGeom)

ggplot(dataGraph, aes(x = I, y = E)) + geom_line(size = 1) + geom_point(data = dataGraph2, aes(x = I, y = E), color = "gold", size = 2.5) + theme_bw() + geom_vline(aes(xintercept = 660/(36 + 1)), color = "gray", linetype=2, size = 1) + xlab("Number of Intervals (I)") + ylab("Bout Length Error \nfor Task 1 (E)") + theme(text = element_text(size=12.5))

IVec = c(1, 2, 4, 8, 16, 32, 64, 128)
wVec = n/IVec
EPred = c()
for(i in 1:length(wVec)){
  if(ltwo > (wVec[i]-1)){
    EPred[i] = ((-lone^2-lone*ltwo+lone*wVec[i])^2)/((lone+ltwo)^2)
  } else {
    EPred[i] = ((-lone^2+lone)^2)/((wVec[i]-1+lone)^2)
  }
}

Rsquared = 1-sum((boutLengthErrorGeom-EPred)^2)/sum((boutLengthErrorGeom-mean(boutLengthErrorGeom))^2)
