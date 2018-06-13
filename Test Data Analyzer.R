nameVector <- c("BY4743", "DOA1", "ERG5", "EDE1", "MIR1")
timePointVector <- c(2, 4, 6, 9, 11, 13)

day2Data <- generateDayGrowthCurves("Day2.txt", nameVector)

day4Data <- generateDayGrowthCurves("Day4.txt", nameVector)

day6Data <- generateDayGrowthCurves("Day6.txt", nameVector)

day9Data <- generateDayGrowthCurves("Day9.txt", nameVector)

day11Data <- generateDayGrowthCurves("Day11.txt", nameVector)

day13Data <- generateDayGrowthCurves("Day13.txt", nameVector)

totalDataVector <-
  list(day2Data, day4Data, day6Data, day9Data, day11Data, day13Data)

generateAllTimeProgressions(totalDataVector, timePointVector, nameVector)

displayDayGrowCurvesCombined(day2Data)

displayDayGrowCurves(day2Data)
day2Data

BY4743Matrix <-
  generateStrainMatrix(totalDataVector, timePointVector, 2, nameVector[1])

calculateDoublingTimeStrain(day2Data, 2)
generateSurvivalCurve(BY4743Matrix, timePointVector, nameVector[1])

generateAllSurvivalCurves(totalDataVector, timePointVector, nameVector)


generateAllSurvivalCurves(myTotalDataVector, myTimePointVector, myNameVector)
#INPUTS A A TXT DOCUMENT OF THE DATA TAB-DELIMITED WITH SPACING REMOVED BETWEEN THE NAMES
#OF WELLS
generateDayGrowthCurves <- function(file, nameVector)
{
  #Importing Data from DAT 2 Document
  setwd("D:/R Stuff/Test Data")
  
  testWellDF = read.table(file,  fill = TRUE, header = TRUE)
  testWellDF
  #Getting initial times as character
  timeVector <- testWellDF[, "Time"]
  timeVectorCha <- as.character(timeVector)
  
  #Converting characters to numeric times in seconds
  numericTimeVector <- sapply(strsplit(timeVectorCha, ":"),
                              function(x) {
                                x <- as.numeric(x)
                                3600 * x[1] + 60 * x[2] + x[3]
                              })
  
  #Convert the time column into pure seconds format
  testWellDF$Time <- numericTimeVector
  
  print(nameVector)
  
  #Getting blank vector
  blankVector <- testWellDF[, "BLANK"]
  
  #Getting Matrix from Data Frame
  testWellMatrix <- data.matrix(testWellDF)
  
  
  #Creating normalized matrix to Contain post blank subtraction data
  normalizedTestWellMatrix <- testWellMatrix
  testWellMatrix
  #Iterating through all rows and substracting out all of the blank values
  for (row in 1:nrow(testWellMatrix))
  {
    for (col in 2:ncol(testWellMatrix))
    {
      normalizedTestWellMatrix[row, col] <-
        testWellMatrix[row, col] - blankVector[row]
    }
  }
  normalizedTestWellMatrix
  #Name Vector containing the strains of yeast/deletino
  
  
  #for(i in 2:(ncol(testWellMatrix)-1))
  #{
  #plot(normalizedTestWellMatrix[,1], normalizedTestWellMatrix[,i], main=paste("Well 1",(i-1)),
  #    xlab="Time (s) ", ylab="Concentration (OD)", pch=1)
  #}
  
  #Normalized averageNTWS where the triplicate samples are joined together
  averageNTWS <-
    matrix(
      nrow = nrow(normalizedTestWellMatrix),
      ncol = ncol(normalizedTestWellMatrix) / 3 + 1,
      byrow = TRUE
    )
  
  dimnames(averageNTWS) = list(1:nrow(averageNTWS), c("Time", nameVector))
  
  #First column is time vector
  for (i in 1:nrow(normalizedTestWellMatrix))
  {
    averageNTWS[i, 1] <- normalizedTestWellMatrix[i, 1]
  }
  
  #Inputing averaged triplicate data:
  for (row in 1:nrow(averageNTWS))
  {
    for (col in 2:ncol(averageNTWS))
    {
      averageNTWS[row, col] = (
        normalizedTestWellMatrix[row, 3 * col - 4] + normalizedTestWellMatrix[row, 3 *
                                                                                col - 3] + normalizedTestWellMatrix[row, 3 * col - 2]
      ) / 3.0
    }
  }
  
  #Return the matrix containing the averaged, normalized data
  averageNTWS
  
}

displayDayGrowCurves <- function(dayMatrix)
{
  nameVector <- colnames(dayMatrix)
  
  for (i in 2:ncol(dayMatrix))
  {
    plot(
      dayMatrix[, 1],
      dayMatrix[, i],
      main = paste("Strain", nameVector[i]),
      xlab = "Time (s) ",
      ylab = "Concentration (OD)",
      pch = 1,
      col = i
    )
  }
  
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

generateStrainCurveProgression <-
  function(TimePointDataVector,
           timeVector,
           strainIndex,
           nameVector)
  {
    strainPointMatrix <-
      matrix(
        nrow = nrow(TimePointDataVector[[1]]),
        ncol = 2 * ncol(TimePointDataVector[[1]]),
        byrow = TRUE
      )
    for (i in seq(1, length(TimePointDataVector), 1))
    {
      dayMatrix <- TimePointDataVector[[i]]
      strainPointMatrix[, 2 * i - 1] <- dayMatrix[, 1]
      strainPointMatrix[, 2 * i] <- dayMatrix[, strainIndex]
    }
    
    plot(
      strainPointMatrix[, 1],
      strainPointMatrix[, 2],
      main = paste(nameVector[(strainIndex - 1)], "Progression"),
      xlab = "Time (s) ",
      ylab = "Concentration (OD)",
      pch = 2,
      col = 2
    )
    
    for (i in 3:(ncol(strainPointMatrix) - 1))
    {
      points(strainPointMatrix[, i],
             strainPointMatrix[, i + 1],
             pch = i,
             col = i)
    }
    
    
  }

generateAllTimeProgressions <-
  function(totalDataVector,
           timePointVector,
           nameVector)
  {
    for (i in 1:length(nameVector) + 1)
    {
      generateStrainCurveProgression(totalDataVector, timePointVector, i, nameVector)
    }
  }

calculateDoublingTimeStrain <- function(day2Data, strainIndex)
{
  doublingTimeSums <- 0
  counter <- 1
  print(day2Data)
  startRow <- 0
  finalRow <- 0
  initialPointFound <- FALSE
  
  #Find first value over/equal to 0.3 OD
  while (!initialPointFound)
  {
    if (day2Data[counter, strainIndex] >= 0.2) {
      initialPointFound <- TRUE
      
      startRow <- counter
      
    }
    else
    {
      counter <- counter + 1
    }
  }
  
  counter <- 1
  
  finalPointFound <- FALSE
  
  #Find last value below or equal to 0.5
  while (!finalPointFound)
  {
    if (day2Data[counter, strainIndex] >= 0.5) {
      finalPointFound <- TRUE
      
    }
    else
    {
      finalRow <- counter
      
      counter <- counter + 1
    }
  }
  
  print(startRow)
  print(finalRow)
  
  totalDoublingTimeSum <- 0
  counter <- 0
  for (i in startRow:(finalRow - 1))
  {
    firstPoint <- day2Data[i, strainIndex]
    secondPoint <- day2Data[i + 1, strainIndex]
    
    firstTime <- day2Data[i, 1]
    secondTime <- day2Data[i + 1, 1]
    
    componentDoublingTime <-
      (log(2)) * (secondTime - firstTime) / (log(secondPoint) - log(firstPoint))
    totalDoublingTimeSum <-
      totalDoublingTimeSum + componentDoublingTime
    counter = counter + 1
    
  }
  
  print((1 / 60) * totalDoublingTimeSum / counter)
  #Return Average Doubling Time in Seconds
  totalDoublingTimeSum / counter
}

generateLinearizedGrowthCurve <-
  function(dayMatrix, strainIndex, name, day)
  {
    plot(
      dayMatrix[, 1],
      log(dayMatrix[, strainIndex]),
      main = paste("Linearized", name, "Day", day),
      xlab = "Time (s) ",
      ylab = "ln(Concentration) ln((OD))",
      pch = 2,
      col = 2
    )
  }

calculateComparisonODTime <- function(dayMatrix, strainIndex)
{
  #Finds the time at which the a particualr growth curve reaches 0.3 OD by
  # utilizing a linearized regression between two points (since the relationship is exponential, the log is taken)
  
  print(dayMatrix)
  #Find Point Right below 0.3 OD
  counter <- 1
  initialPointFound <- FALSE
  startRow <- 1
  #Find last value under to 0.3 OD
  while (!initialPointFound)
  {
    if (dayMatrix[counter, strainIndex] > 0.3) {
      initialPointFound <- TRUE
      
    }
    else
    {
      startRow <- counter
      
      counter <- counter + 1
    }
  }
  
  print(startRow)
  
  initialPointFound <- FALSE
  
  counter <- 1
  endRow <- 1
  
  #Find first value over/equal to 0.3 OD
  while (!initialPointFound)
  {
    if (day2Data[counter, strainIndex] > 0.3) {
      initialPointFound <- TRUE
      
      endRow <- counter
      
    }
    else
    {
      counter <- counter + 1
    }
  }
  
  print(endRow)
  
  OD1 <- dayMatrix[startRow, strainIndex]
  OD2 <- dayMatrix[endRow, strainIndex]
  
  time1 <- dayMatrix[startRow, 1]
  time2 <- dayMatrix[endRow, 1]
  
  regressionSlope <- (log(OD2) - log(OD1)) / (time2 - time1)
  timePoint0.3 <- (time1 + (log(0.3) - log(OD1)) / (regressionSlope))
  timePoint0.3
}

generateStrainMatrix <-
  function(TimePointDataVector,
           timeVector,
           strainIndex,
           name)
  {
    strainPointMatrix <-
      matrix(
        nrow = nrow(TimePointDataVector[[1]]),
        ncol = 2 * ncol(TimePointDataVector[[1]]),
        byrow = TRUE
      )
    for (i in seq(1, length(TimePointDataVector), 1))
    {
      dayMatrix <- TimePointDataVector[[i]]
      strainPointMatrix[, 2 * i - 1] <- dayMatrix[, 1]
      strainPointMatrix[, 2 * i] <- dayMatrix[, strainIndex]
    }
    
    # Generate a column name vector
    columnNames <-
      vector(mode = "character", length = ncol(strainPointMatrix))
    for (i in seq(1, length(columnNames), 2))
    {
      columnNames[i] <- name
    }
    
    for (i in seq(2, length(columnNames), 2))
    {
      columnNames[i] <- paste("Day", timeVector[i / 2])
    }
    
    colnames(strainPointMatrix) <- columnNames
    
    print(strainPointMatrix)
    
    strainPointMatrix
  }


generateSurvivalCurve <-  function(strainMatrix, timeVector, name)
{
  # Calculate the Doubling Time fromt the Strain's day 2 Data
  day2Data <-
    matrix(nrow = nrow(strainMatrix),
           ncol = 2,
           byrow = TRUE)
  
  day2Data[, 1:2] <- strainMatrix[, 1:2]
  
  doublingTime <- calculateDoublingTimeStrain(day2Data, 2)
  
  #Calculcate the time when each day reaches 0.3 OD
  day2TimeToReach0.3OD <- calculateComparisonODTime(day2Data, 2)
  print(day2TimeToReach0.3OD)
  thresholdTimes <- c(day2TimeToReach0.3OD)
  for (i in seq(3, ncol(strainMatrix), 2))
  {
    dayXData <-
      matrix(nrow = nrow(strainMatrix),
             ncol = 2,
             byrow = TRUE)
    dayXData[, 1:2] <- strainMatrix[, (i):(i + 1)]
    timeTo0.3OD <- calculateComparisonODTime(dayXData, 2)
    thresholdTimes <- c(thresholdTimes, timeTo0.3OD)
  }
  
  print(thresholdTimes)
  
  #Generate Time Shift Vector
  timeShiftVector <- thresholdTimes - thresholdTimes[1]
  
  survivabilityVector <- 1 / (2 ^ (timeShiftVector / doublingTime)) * 100
  print(survivabilityVector)
  
  surivivalMatrix <- cbind(timeVector, survivabilityVector)
  
  
  plot(
    surivivalMatrix[, 1],
    surivivalMatrix[, 2],
    main = paste(name, "Survival Curve"),
    xlab = "Time Point (Day) ",
    ylab = "Survivability Percentage",
    pch = 2,
    col = 2
  )
  lines(surivivalMatrix[, 1], surivivalMatrix[, 2], type = "o")
  
  surivivalMatrix
}