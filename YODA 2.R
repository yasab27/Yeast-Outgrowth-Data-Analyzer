# YEAST OUTGROWTH DATA ANALYZER

# This code is designed to generate viability curves for strains based
# on yeast outgrowth data from a bioscreen CLS machine. It inputs
# ODs from a specific day, normalizes them with respect to the the blank
# value, and then uses the information from multiple days to generate
# a viability curve.

# Yasa Baig -- 6/7/18 -- Harvard University

#--------------------------------------------------------------------------------------#


# - - - - - - - - - - - - - - - - - - - RUN - - - - - - - - - - - - - - - - - - - - - -#
# =====================================================================================#
# - - - - - - - - - - - - - - - - - - - RUN - - - - - - - - - - - - - - - - - - - - - -#


# DATA INPUT:

# Names of each strain of yeast being analyzed. Generally the BY4743 strain will refer to
# a control name of yeast.
#
# TO-DO: Have the program be able to yield the name of the strains from a generalized
# file
strainNames <- c("BY4743", "DOA1", "ERG5", "EDE1", "MIR1")

# Vector of the day time points (Day 2, Day 3, etc.)
timePointTimes <- c(2, 4, 6, 9, 11, 13)

# Generates a vector containing matrices of all the data from each individual day point
allDataList <- generateAllDataList(timePointTimes)
allDataList

#--------------------------------------------------------------------------------------#

# ANALYSIS:
survivalMatrix <- generateAllSurvivalCurves(allDataList, strainNames, timePointTimes)
survivalMatrix
plotSurvivalCurves(survivalMatrix, timePointTimes, strainNames)

survivalIntegrals <- computeSurvivalIntegrals(allDataList, strainNames, timePointTimes)

strainSurvivalIntegrals <- computeStrainSurvivalIntegrals(survivalIntegrals, strainNames)
strainSurvivalIntegrals



# - - - - - - - - - - - - - - - - - - -/RUN - - - - - - - - - - - - - - - - - - - - - -#
# =====================================================================================#
# - - - - - - - - - - - - - - - - - - -/RUN - - - - - - - - - - - - - - - - - - - - - -#


#--------------------------------------------------------------------------------------#

# HELPER FUNCTIONS:

# DATA INPUT HELPERS:

# INPUT:
# -timePointTimes: The individual days (2, 3, etc. ) which will be used for naming each component
# matrix in the final return list to once again order increased help in reading the data.
#
# OPERATIION:
# This function reads all of the data from .txt files using a helper function in order
# to create one large list of all the data with appropriate names for each matrix
# and column in each matrix.
#
# OUTPUT:
# Returns a list containing matrices for each individual day.
#
generateAllDataList <- function(timePointTimes)
{
  # Calls on a helper function for each individual day to generate data for that day.
  # For the time being, this is hard coded in (each individual day matrix is created "by hand")
  #
  # TO-DO: Use a large data structure and iteration to be able to create the final matrix
  # without having to use hard coding. This will likely require some form of GUI for file input. (R-Shiny? )
  day2Data <- generateDayDataMatrix("Day2.txt")
  day4Data <- generateDayDataMatrix("Day4.txt")
  day6Data <- generateDayDataMatrix("Day6.txt")
  day9Data <- generateDayDataMatrix("Day9.txt")
  day11Data <- generateDayDataMatrix("Day11.txt")
  day13Data <- generateDayDataMatrix("Day13.txt")
  
  # Create a list containing each and every data for each well from each and every day and store it in one
  # large list. This large data is all cleaned (normalized with respect to the blank and has the time converted to
  # seconds).
  allDaysDataList <-
    list(day2Data, day4Data, day6Data, day9Data, day11Data, day13Data)
  
  # Assigning a name to each element of the list which corresponds to the Day Time Point for each set of Data.
  # This allows for much more easily human interpretable data.
  names(allDaysDataList) <-
    paste("Day", timePointTimes, "Cleaned Data")
  
  # Return the cleaned "SUPER LIST" containing matrices for each day point.
  allDaysDataList
}

# INPUT:
# -fileName: (Char)
#
# OPERATIION:
# This function will read the data from one specific day, convert the initial column from String literals to seconds (00:00:05 -> 5)
# to allow for future calculations. Each column will also be named with the appropriate name based on the strain name. It is assumed
# each strain is in triplicate. Additionally, normalizes each OD with respect to the blank value.
#
# OUTPUT:
# Returns a Data Frame of normalized OD times for one day with appropriate column names and times converted to seconds.
#
generateDayDataMatrix <- function(fileName)
{
  # Importing Data from DAY 2 Documen
  
  # Inputting a Data Frame by reading the table on the .txt document.
  testWellDF = read.table(fileName,  fill = TRUE, header = TRUE)
  
  # Acquiring the last the column of time measurements for the day as char
  # values in order to convert them into numeric second values.
  timeVector <- testWellDF[, "Time"]
  timeVectorCha <- as.character(timeVector)
  
  # Converting characters to numeric times in seconds by splitting the strings and
  # then multiplying by teh appropriate constants.
  numericTimeVector <- sapply(strsplit(timeVectorCha, ":"),
                              function(x) {
                                x <- as.numeric(x)
                                3600 * x[1] + 60 * x[2] + x[3]
                              })
  
  # Assigning the initial column to the new numeric times (eliminating the old 00:00:00 format)
  # This is done to allow for later calculations.
  testWellDF$Time <- numericTimeVector
  
  # Normalizing each well by subtractin the blank value for each row (lineage). This compensates for the amount of "background"
  # OD for each well, allowing for only the yeast concentration changes be responsible for any changes.
  
  # Iterate through each and every row and subtract from each column (except for the time column) the BLANK value for that row
  # Store this new matrix of normalized well ODS in "testWellDFN" (Test Well Data Frame Normalized)
  testWellDFN <-
    testWellDF[1:nrow(testWellDF), 2:ncol(testWellDF)] - testWellDF[1:nrow(testWellDF), ncol(testWellDF)]
  
  # Appending the "Time" Values from the initial time point into the initial column of the Data Frame
  # Initializing it as a new vector first allows for the column name to be "Time" as opposed to a longer, less helpful
  # Data Frame Subname.
  Time <- testWellDF$Time
  testWellDFN <- cbind(Time, testWellDFN)
  
  # Eliminating the empty "BLANK" column which only contains the number "0" multiple times.
  testWellDFN <- testWellDFN[, 1:ncol(testWellDFN) - 1]
  
  # Return the entire normalized and cleaned Data Frame  as a numeric matrix to the client function
  testWellMatrixN <- data.matrix(testWellDFN, rownames.force = NA)
  testWellMatrixN
}

#--------------------------------------------------------------------------------------#

# ANALYSIS HELPERS:

# INPUT:
# -allDataList: A list of matrices containing the cleaned data for each and every day time point (Day 2, etc. )
#
# OPERATION:
# This functions computes the survival integral for a particular survival curve utilizing a trapezoidal
# sum approximation. 
#
# OUTPUT:
# Outputs a matrix of both the survival integrals for each strain and the standard deviation of each 
# strain.
computeSurvivalIntegrals <- function(allDataList, strainNames, timePointTimes)
{
  # Generate a matrix of all threshold times (times for specific outgrowth curves to reach 0.3 OD)
  # These will be used to for comparison in time shifts in order to generate a survival curve.
  totalThresholdMatrix <- generateThresholdTimeMatrixTotal(allDataList, timePointTimes)
  totalThresholdMatrix
  
  # Generate a matrix of all the doubling times for each well, calculated from the data in the
  # day 2 data matrix. 
  doublingTimeMatrix <- calculateDoublingTimesMatrix(allDataList[[1]])
  doublingTimeMatrix
  
  # Generate the first column of the matrix for the survival curve. 
  survivalMatrix <- generateSurvivalCurve( doublingTimeMatrix[,1], totalThresholdMatrix[,1] )
  
  # Complete generation of the survival Curve for the remaining well
  for(i in 2: ncol(doublingTimeMatrix))
  {
    survivalMatrix <- cbind(survivalMatrix, generateSurvivalCurve( doublingTimeMatrix[,i], totalThresholdMatrix[,i] ))
  }
  
  # Assign the names of the new survival matrix to have the well names.
  wellNameVector <- colnames(doublingTimeMatrix)
  colnames(survivalMatrix) <- wellNameVector
  survivalMatrix
  
  # Generate a list of survival integrals for each well
  
  survivalIntegralMatrix <- computeSurvivalIntegral(survivalMatrix[,1],timePointTimes)
  
  for(i in 2:ncol(survivalMatrix))
  {
    survivalIntegralMatrix <- cbind(survivalIntegralMatrix, computeSurvivalIntegral(survivalMatrix[,i],timePointTimes))
  }
  
  colnames(survivalIntegralMatrix) <- wellNameVector
  survivalIntegralMatrix
}

# Input the survival values for each time point for one well
# and compute the survival integral for that well. 
computeSurvivalIntegral <- function(wellSurvivalValues, timePointTimes)
{
  # First, the x-axis (time Point Length) differences must be found for use in the trapezoidal
  # sum approximation. This is done by subtracting consecutive points. 
  
  # First the "left points" must be acquired (everyPoint but the last one)
  leftPoints <- timePointTimes[1:(length(timePointTimes) - 1) ]

  # Next the "right points" must be acquired (every point by first one)
  rightPoints <- timePointTimes[-1]

  # The differences must now be found to get the "heights" for each trapezoid
  trapezoidHeights <- rightPoints - leftPoints

  # The same process must now be repeated for the "heights" of the trapzeoids. 
  
  # The left heights (EVERY POINT BUT THE LAST)
  leftHeights <- wellSurvivalValues[-length(wellSurvivalValues)]

  # The right heights (Every Point but the first)
  rightHeights <- wellSurvivalValues[-1]

  # To get the height for the trapezoid, we now must average the value of the consecutive poitns
  averagedTrapezoidHeights <- (leftHeights+rightHeights)/2

  # Now the areas of each subtrapezoid is calcuclated
  componentTrapezoidHeights <- trapezoidHeights * averagedTrapezoidHeights

  # To get the entire survival integral, now we must just add up all the component trapezoids
  survivalIntegral <- sum(componentTrapezoidHeights)
  survivalIntegral
  
}

# Average the data in sets of three in order to yield the average survival integral and 
# standard deviation for each strain
computeStrainSurvivalIntegrals <- function(survivalIntegrals, strainNames)
{
  # Start by averaging the first three triplicate points
  # Store the standard deviations as well
  survivalIntegrals
  
  averageStrainSurvivalIntegrals <- mean(survivalIntegrals[,1:3])
  
  sdStrainSurivalIntegrals <- sd(survivalIntegrals[,1:3])

  # Calculate the averages for the remaining three means. 
  for( x in seq( 4, (ncol(survivalIntegrals)-2) ,3 ) )
  {
    lowBound <- x
    upBound <- x+2
    holderMatrix <- mean(survivalIntegrals[, lowBound:upBound ])
    averageStrainSurvivalIntegrals <- cbind(averageStrainSurvivalIntegrals,holderMatrix  )
    holderSDMatrix <- sd(survivalIntegrals[, lowBound:upBound ])
    sdStrainSurivalIntegrals <- cbind(sdStrainSurivalIntegrals,holderSDMatrix  )
    
  }
  
  colnames(averageStrainSurvivalIntegrals) <- strainNames
  strainSurvivalStats <- rbind(averageStrainSurvivalIntegrals, sdStrainSurivalIntegrals)
  
  rownames(strainSurvivalStats) <- c("Mean SI", "SD")
  strainSurvivalStats
}

# INPUT:
# -surivalMatrix: This is the final super calculated survival matrix containing the averaged triplicate
# value of survivability at each time point. 
#
# OPERATION:
# This functions generates survival plots for each and every strain utilizing the precalculated survival curves
#
# OUTPUT:
# Outputs plots on the bottom right.
plotSurvivalCurves <- function(survivalMatrix, timePointTimes, nameVector)
{
  for(i in 1:length(nameVector))
  {
   plotSurvivalCurve(survivalMatrix[,i], timePointTimes, nameVector[i])
  }
}

# Computes the same but only for one well. This function is used iteratively
# in the prior function in order to allow for easier procedural abstraction.
plotSurvivalCurve <- function(wellSurvivalVector, timePointTimes, strainName)
{
  
  plot(
    timePointTimes,
    wellSurvivalVector,
    main = paste(strainName, "Survival Curve"),
    xlab = "Time Point (Day) ",
    ylab = "Survivability Percentage",
    ylim = c(0,100),
    pch = 6,
    col = 3
  )
  polygon(c(2,2,timePointTimes, timePointTimes[length(timePointTimes)]),c(0,100,wellSurvivalVector,0),col="brown")
  
  lines(timePointTimes, wellSurvivalVector, type = "o")
}

# INPUT:
# -allDataList: A list of matrices containing the cleaned data for each and every day time point (Day 2, etc. )
# -strainNames: A char vector of all the strain names.
#
# OPERATION:
# This functions generates survival plot data for each and every strain utilizing the the outgrowth data from each well.
#
# OUTPUT:
# Outputs plots on the bottom right.
generateAllSurvivalCurves <- function(allDataList, strainNames, timePointTimes)
{
  # Generate a matrix of all threshold times (times for specific outgrowth curves to reach 0.3 OD)
  # These will be used to for comparison in time shifts in order to generate a survival curve.
  totalThresholdMatrix <- generateThresholdTimeMatrixTotal(allDataList, timePointTimes)
  totalThresholdMatrix
  
  # Generate a matrix of all the doubling times for each well, calculated from the data in the
  # day 2 data matrix. 
  doublingTimeMatrix <- calculateDoublingTimesMatrix(allDataList[[1]])
  doublingTimeMatrix
  
  # Generate the first column of the matrix for the survival curve. 
  survivalMatrix <- generateSurvivalCurve( doublingTimeMatrix[,1], totalThresholdMatrix[,1] )
  
  # Complete generation of the survival Curve for the remaining well
  for(i in 2: ncol(doublingTimeMatrix))
  {
    survivalMatrix <- cbind(survivalMatrix, generateSurvivalCurve( doublingTimeMatrix[,i], totalThresholdMatrix[,i] ))
  }
  
  # Assign the names of the new survival matrix to have the well names.
  wellNameVector <- colnames(doublingTimeMatrix)
  colnames(survivalMatrix) <- wellNameVector
  survivalMatrix
  
  # Average every three consecutive columns in order to generate an averaged matrix for a specific strain
  # since consecutive wells are in triplicate.
  
  # Start by first finding the first strain
  strainMatrix <- rowMeans(survivalMatrix[,1:3])
  
  # Calculate the averages for the remaining three means. 
  for( x in seq( 4, (ncol(survivalMatrix)-2) ,3 ) )
  {
    lowBound <- x
    upBound <- x+2
    holderMatrix <- rowMeans(survivalMatrix[, lowBound:upBound ])
    strainMatrix <- cbind(strainMatrix,holderMatrix  )
  }
  
  # Assign the names of the column matrix to the strains of yeast.
  colnames(strainMatrix) <- strainNames
  
  # Return the strain matrix which contains the survival percentages at each time point. 
  strainMatrix
}

# INPUT:
# -doublingTime: The doubling time for the well who's survival values are being calculated.
# -threshHoldWellVector: A vector containing the threshold times (times to reach 0.3 OD for each time point) for one speciic well.
#
# OPERATION:
# This functions generates a matrix of survivability percentages for each time point for the specific well. 
#
# OUTPUT:
# A matrix of survivability values for the specific well in question, calculated from doubling times and 
# the threshold well vectors.
generateSurvivalCurve <- function(doublingTime, threshHoldWellVector)
{
  # Calculate the time shifts for each time point by subtracting the distance from 
  # the day 2 (first data point) threshold time. 
  timeShiftVector <- threshHoldWellVector - threshHoldWellVector[1]
  
  # Generate the survivability percentage at each time point by using the formula
  # Survivability = (1/2^(TS/DT)) * 100
  survivabilityVector <- (1/2^(timeShiftVector/doublingTime))*100
  survivabilityVector
  
}

# INPUT:
# -dayTwoData: The matrix of initial values from the day 2 time point which will be used to calculate the
# doubling times for the entire matrix. 
#
# OPERATION:
# Uses the doubling time formula to calculcate the doubling time for all wells of yeast and then averages each 
# doubling time to yield an accurate doubling time for the yeast. 
#
# OUTPUT:
# Returns a matrix representing the doubling time for each well. 
calculateDoublingTimesMatrix <- function(dayTwoData)
{
  doublingTimeMatrix <- calculateDoublingTimesWell(dayTwoData[,1], dayTwoData[,2])
  
  for(i in 3:ncol(dayTwoData))
  {
    newDoublingTime <- calculateDoublingTimesWell( dayTwoData[,1] , dayTwoData[,i])
    doublingTimeMatrix <- cbind(doublingTimeMatrix,newDoublingTime)
  }
  
  wellNameVector <- colnames(dayTwoData)
  colnames(doublingTimeMatrix) <- wellNameVector[-1]
  rownames(doublingTimeMatrix) <- "DTs"
  doublingTimeMatrix
}

# INPUT:
# -wellVector: The vector of a particular well OD times taken from a larger day 2 dataset. 
# -timeVector: The vector of time in seconds for the day 2 data.
#
# OPERATION:
# Uses the doubling time formula to calculcate the doubling time for a specific well of yeast and then averages each 
# doubling time to yield an accurate doubling time for the yeast. 
#
# OUTPUT:
# Returns a numeric representing the doubling time for a specific well. 
calculateDoublingTimesWell <- function(timeVector, wellVector)
{
  # In order to calculate doubling time, we must calculate the individual doubling times for consecutive times between 
  # 0.2 and 0.5 ODs (data before and after these points are considered to be innacurate).
  
  # First we must determine the indices between 0.2 OD and 0.5 OD 
  
  # Finding the lower bound. This corresponds to the smallest value larger or greater than0 0. 2 OD.
  lowerBoundIndex <- min(which(wellVector >= 0.2))
  
  # Finding the upper bound. This corresponds to the largest value below or equal to 0.5 OD.
  upperBoundIndex <- max(which(wellVector <= 0.5))
  
  # Acquire the data values from the first the starting index to the value right before the final index
  lowerODRange <- wellVector[lowerBoundIndex:(upperBoundIndex-1)]
  lowerTimeRange <- timeVector[lowerBoundIndex:(upperBoundIndex-1)]
  
  # Acquire the data values from one after the starting index to the final index
  upperODRange <- wellVector[(lowerBoundIndex+1):upperBoundIndex]
  upperTimeRange <- timeVector[(lowerBoundIndex+1):upperBoundIndex]
  
  # Create a vector of the doubling times calculated between consecutive OD measurements in the predefined range.
  doublingTimeVector <- ((log(2))*(upperTimeRange - lowerTimeRange))/(log(upperODRange) - log(lowerODRange))
  
  # Average the doubling time values to yield one final average doubling time
  wellDoublingTime <- mean(doublingTimeVector)
  
  # Return the well doubling time in seconds. In general, times will range from 85 - 90 minutes for most
  # wild type yeast. 
  wellDoublingTime
}
  

# INPUT:
# -allDataList: A list of matrices containing the cleaned data for each and every day time point (Day 2, etc. )
#
# OPERATION:
# Generates a matrix of all threshold times for each and every day utilizing the below helper functions. 
#
# OUTPUT:
# Returns a matrix for each threshold tiem for each well for each time point, each time point appropriately labeled. 
generateThresholdTimeMatrixTotal <- function(allDataList , timePointTimes)
{
  # Generate a matrix of every threshold for each well for each time point by starting at the first day.
  allTimePointThresholdMatrix <- generateThresholdTimeMatrix(allDataList[[1]])
  
  # Iterating through every other member of the list, adding a new row to the new matrix, until finally a matrix containing 
  # each and every threshold time is generated and added to the list.
  for( i in 2:length(allDataList))
  {
    allTimePointThresholdMatrix <- rbind(allTimePointThresholdMatrix, generateThresholdTimeMatrix(allDataList[[i]]))
  }
  
  # Name the rows with appropriate day time point (Day 2, Day 4, etc.)
  rownames(allTimePointThresholdMatrix) <- paste("Day",timePointTimes)
  allTimePointThresholdMatrix
}

# INPUT:
# -timePointMatrix: A matrix of outgrowth data which is cleaned (generated from an earlier function). The times should be in seconds
# and all data should be normalized.
#
# OPERATION:
# This functions finds when each well in a specific time point reaches it reaches the threshold comparison OD (0.3 OD).
# This information will be later used in order to generate survival curves for specific straings of yeast.
#
# OUTPUT:
# Outputs a matrix for well for the time point matrix reaches 0.3 OD in seconds.
generateThresholdTimeMatrix <- function(timePointMatrix)
{
  # Acquire the left column time vector. This is the same for each well in the matrix, so it can be considered
  # an external variable. Note this is always the first column of the matrix, so it can be retrieved using [,1]
  timeVector <- timePointMatrix[, 1]
  timePointMatrix
  # Calls another helper function which inputs the time vector along with one well. Finds the threshold time
  # for that specific well.
  
  # Generate the threshold time for the first well in that specific time point.
  allWellThresholdMatrix <- findThresholdTime(timeVector, timePointMatrix[,2])
  
  # Iteratue through the remaining wells, generating a vector containing the threshold time for each and every well. 
  for(i in 3:ncol(timePointMatrix))
  {
    allWellThresholdMatrix <- c(allWellThresholdMatrix, findThresholdTime(timeVector, timePointMatrix[,i]))
  }
  # Assign the names of each well to the new vector of all Threshold times. 
  nameVector <- colnames(timePointMatrix)
  nameVector <- nameVector[-1]
  
  names(allWellThresholdMatrix) <- nameVector
  allWellThresholdMatrix
}

# INPUT:
# -timeVector: A vector containing time points in seconds corresponding to Well OD Readings.
# -wellVector: A vector containing OD readings for a specific well for a specific time point.
#
# OPERATION:
# This function will calculate when a specific well reaches 0.3 OD (in seconds) by computing a logarithmic
# regression around the two points braketing 0.3 OD in the well. As the outgrowth of the yeast is approximately
# exponential during the first phase of growth, a logarithmic regression approximates the time very well.
#
# OUTPUT:
# Returns the time for which a specific well reaches 0.3 OD.
findThresholdTime <- function(timeVector, wellVector)
{
  
  #print(cbind(timeVector, wellVector))
  
  # The first step is to find the values braketing 0.3 OD in order to compute the logarithmic-linear regression.
  # Finding the index of the highest value below 0.3 OD.
  lowerBoundIndex <- max(which(wellVector < 0.3))
  
  #print(lowerBoundIndex)
  
  # Finding the index of the smallest value above 0.3 OD.
  upperBoundIndex <- min(which(wellVector > 0.3))
  
  #print(upperBoundIndex)
  
  # Defining the independent and dependent variables which will be used to generate
  # the linear model
  
  # First, the time values must be acquired utilizing the indexes calculcated above
  timeLowerBound <- timeVector[lowerBoundIndex]
  timeUpperBound <- timeVector[upperBoundIndex]
  timeDomain <- c(timeLowerBound, timeUpperBound)
  # Next the corresponding OD values should be acquired at the same vectors
  ODLowerBound <- wellVector[lowerBoundIndex]
  ODUpperBound <- wellVector[upperBoundIndex]
  
  # For the purposes of linear modeling, the natural log of both of the OD Values
  # will be stored into the final vector, as this will result in a more natural
  # regression model.
  ODRange <- c(log(ODLowerBound), log(ODUpperBound))
  
  # Generating a linear regression model of ln(OD) vs time (s). The coefficients and
  # intercepts of this model will be used to calculcate the specific time for when
  # the well reaches 0.3 OD.
  
  # Note that this regression predicts time from OD values, and not the reverse. This
  # is because in order to minimize error, once must use a least-squares regression model
  # which only minimizes error in terms of the response variable, which in this case
  # must be time as that is the value we are trying to predict.
  OcularDensityTimeRegression <- lm(timeDomain ~ ODRange)
  
  # Predicting the time at which the OD reading will reach ln(0.3). This corresponds to the threshold time
  # for this particular well.
  thresholdTime <-
    predict(OcularDensityTimeRegression, data.frame(ODRange = log(0.3)))
  thresholdTime
}

displayDayGrowCurvesCombined <- function(dayMatrix)
{
  plot(
    dayMatrix[, 1],
    dayMatrix[, 2],
    main = "Combined Strains",
    xlab = "Time (s) ",
    ylab = "Concentration (OD)",
    pch = 2,
    col = 2
  )
  
  for (i in 3:(ncol(dayMatrix)))
  {
    points(dayMatrix[, 1], dayMatrix[, i], pch = i, col = i)
  }
  
}
