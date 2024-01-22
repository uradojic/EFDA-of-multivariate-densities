############################################################################
#setwd("~/Library/CloudStorage/OneDrive-Personal/Notes/school/bc/statnice/bakalarska_prace/lr_ddc_algorithm")
############################################################################
#
#
#     Read me
#
#
#  D: The data matrix (required input argument).
#
#  alpha: Specifies the quantile to be used for flagging cellwise outliers.
# Defaults to 0.95.
#
#  numDiscrete: A variable with values >= numDiscrete will be 
# considered discrete and ignored for flagging cellwise outliers.
# Defaults to 3.
#
#  pOutLR: Compositional parts are flagged as cellwise outliers if the 
# percentage of outlying logratios involving that part >= pOutLR.
# Defaults to 0.3.
#
#  pOutRow: Observations are flagged as rowwise outliers if the 
# percentage of outlying cells >= pOutRow. (Or if they are flagged by DDC 
# as rowwise outliers in the data with logratios). Defaults to 0.75.
#
#  showVals: Takes the values "D", "R", "DR", "RD" or NULL and determines whether
# or not to show the entries of the data matrix (D) or the residuals (R) 
# in the cellmap or no values (NULL). If choosing "DR" or "RD" then both
# cellMaps will be ploted side by side. Defaults to "D".
#
#  roundD: If "D" is selected for showVals, then use rounded values of data 
# matrix (D), rounded to integer by default.
#
#  showrows: Indices of the rows to be shown, defaults to NULL which means all
# rows are shown. Defaults to NULL.
#
#  showcolumns: Indices of the columns to be shown, defaults to NULL which
# means all columns are shown. Defaults to NULL.
#
#  nrowsinblock: How many rows are combined in a block. Defaults to 1.
#
#  ncolumnsinblock: How many columns are combined in a block. Defaults to 1.
#
#  rowtitle: Title for the rows. Defaults to NULL.
#
#  columntitle: Title for the columns. Defaults to NULL.
#
#  columnangle: Angle of the column labels. Defaults to 90.
#
#  sizetitles: Size of row title and column title. Defaults to 3.
#
#  adjustrowlabels: Adjust row labels: 0=left, 0.5=centered, 1=right.
# Defaults to 1.
#
#  adjustcolumnlabels: Adjust column labels: 0=left, 0.5=centered, 1=right.
# Defaults to 1.
#
#  datasetName: If value != NULL, then the selected text with chosen parameters
# will be shown in cellmap and used as title of pdf file. Defaults to NULL.
#
# drawCircles: Whether or not to draw circles indicating rowwise outliers.
# Defaults to TRUE.
#
#  mTitle: Main title of the cellMap. Defaults contains datasetName with used
# parameters.
#
#  showCellMap: If TRUE then cellmap will be plotted. Defaults to TRUE.
#
#  pdfCellMap: If TRUE then cellmap will be saved as pdf file with title
# by default mTitle. Defaults to FALSE.
#
#  pdf_width: Width of pdf file. Defaults to 20.
#
#  pdf_height: Height of pdf file. Defaults to 30.
#
############################################################################


############################################################################
#
#
#    LR-DDC (Logratio Deviating Data Cells)
#
#
LR_DDC = function(
  D,
  alpha = 0.95,
  numDiscrete = 3,
  pOutLR = 0.3,
  pOutRow = 0.75,
  showVals = "D",
  roundVal = 0,
  showrows = NULL,
  showcolumns = NULL,
  nrowsinblock = 1,
  ncolumnsinblock = 1,
  rowtitle = NULL,
  columntitle = NULL,
  columnangle = 90,
  sizetitles = 3,
  adjustrowlabels = 1,
  adjustcolumnlabels = 1,
  datasetName = NULL,
  drawCircles = TRUE,
  mTitle = paste(
    "LR-DDC ", datasetName, ": ", used_parameters,
    sep = ""
    ),
  showCellMap = TRUE,
  pdfCellMap = FALSE,
  pdf_width = 20,
  pdf_height = 30
)
  {
#
# Import libraries
#
  library(robCompositions)
  library(cellWise)
  library(robustHD)
  require(gridExtra)
#
#
#  
# Load input data as matrix and save parameters
#
  X = as.matrix(D)
  Xrow = nrow(X)
  Xcol = ncol(X)
#
# Prepare column names of matrix X
#
  namesX = colnames(X)
#
# If matrix has no colnames then name columns are X1, X2, ... XN
#
  if (is.null(namesX)) {
    namesX = paste0("X",seq_len(Xcol))
    colnames(X) = namesX
  } 
#
# Add space at the end of each column name
# so that no name is a subset of another name
#
  list_names = c()
  for (name in namesX){
    name = paste(name, "") 
    list_names = append(list_names, name)
  }
  list_names = list_names
#  
# All combinations of column names with and without space
#
  combination_namesX = combn(
    namesX,
    2,
    simplify = FALSE
    )
  combination_uniq_namesX = combn(
    list_names,
    2,
    simplify = FALSE
    )
#
# Compute pairwise logratios
#
  Xlr = lapply(
    combination_namesX,
    function(j) log(X[, j[1]] / X[, j[2]])
    )
  Xlr = do.call(
    cbind,
    Xlr
    )
  colnames(Xlr) = sapply(
    combination_uniq_namesX,
    function(j) paste(j, collapse = "/")
    )
#  
# Detect outlying cells in pairwise logratios using DDC on dataframe
# made of pairwise logratios and the origin matrix
#
  Xddc = data.frame(Xlr, check.names = FALSE)
  ddc = DDC(
    Xddc, 
    DDCpars = list(
      tolProb = alpha,
      numDiscrete = numDiscrete,
      silent = TRUE
      )
    )
#
# Convert index of outlying cells to array indices
# computed on logratio and origin values
#
  indcells = arrayInd(
    ddc$indcells,
    dim(ddc$remX)
    )
  indcells = data.frame(
    row = indcells[, 1],
    col = colnames(ddc$remX)[indcells[, 2]],
    stringsAsFactors = FALSE
    )
#  
# Indices of outlying cells in the pairwise logratios
#
  indicesXlr = split(
    indcells$row,
    indcells$col
    )
#
# Find indices of outlying cells in compositional parts 
# (the percentage of the outlying logratios involving that part >= pOutLR)
#
  indicesX = sapply(list_names,
                    namesXlr = names(indicesXlr),
                    simplify = FALSE,
                    function(part, namesXlr){
                      
    # 1. Find name from X in logratios names, e.g. acit -> acit/furfu
    # 2. Select all indcells found in all logratios 
    # 3. Make indcells as uniqeu numbers
    # 4. Flag cells if at least pOutLR*NR_Of_Variables are flagged in logratios
                      
    whichLR = grep(
      part,
      namesXlr,
      fixed = TRUE
      ) 
    indicesLR = as.factor(
      unlist(
        indicesXlr[whichLR],
        use.names = FALSE
        )
      )
    indices = as.integer(levels(indicesLR))
    thresholdLR = pOutLR * length(whichLR)
    indices[tabulate(indicesLR) >= thresholdLR]
  }
  )
#
# Find outlying rows 
# 1. Flag whole row if row has more than thresholdRow of outliying values
# 2. Union flagged rows with rows flagged by DDC
# 3. Return indices of cellwise and rowwise outliers
#
  indicesXR = c(
    indicesX,
    recursive = TRUE,
    use.names = FALSE
    )
  indicesXR = as.factor(indicesXR)
  thresholdRow = pOutRow * (Xcol)
  indicesRow = as.integer(levels(indicesXR))
  indicesRow = indicesRow[tabulate(indicesXR) >= thresholdRow]
  indicesRow = union(as.integer(ddc$indrows), indicesRow)
  
  outliers = list(
    indicesX = indicesX,
    indicesRow = indicesRow
    )
#
# Reindex indicesX to know which values of matrix X are outlying (by rows)
#
  for (index in 1:Xcol){
    outliers$indicesX[[index]] = outliers$indicesX[[index]] + Xrow * (index - 1)
  }
#
# Create matrix of standard residuals from pairwise logratios computed by DDC
#
  stdResidPairs = ddc$stdResid[,1:ncol(Xlr)]
  colnames(stdResidPairs) = colnames(Xlr)

  stdResid = matrix(nrow = Xrow, ncol = Xcol)
  for (index in 1:Xcol){
    
    # 1. Union pairwise logratios using mean for each combination of list_names
    # 2. Union logratios where list_name is in numerator with opposite
    # value of logratios where list_name is in denominator
    # 3. Calculate mean of rows and add as column to standard residuals matrix
    
    numerator_log_name = paste(list_names[index], "/", sep = "")
    denominator_log_name = paste("/" , list_names[index], sep = "")
    idx_numenator_in_logpairs = which(
      grepl(
        numerator_log_name,
        colnames(stdResidPairs)
        )
      )
    idx_denominator_in_logpairs = which(
      grepl(
        denominator_log_name,
        colnames(stdResidPairs)
        )
      )
    union_listname_logpairs = cbind(
      stdResidPairs[,idx_numenator_in_logpairs],
      -stdResidPairs[,idx_denominator_in_logpairs]
      )
    stdResid[,index] = apply(
      union_listname_logpairs,
      1, 
      mean #also median would be possible
      )
  }

  colnames(stdResid) = namesX
#
# Create cellMap with flagged outlying rows and cells
# Prepare data
#
  indrows = outliers$indicesRow
  indcells = unlist(outliers$indicesX)
  columnlabels = namesX
  rowlabels = rownames(X)
  used_parameters = paste(
    "alpha = ", alpha,
    ", pOutLR = ", pOutLR,
    ", pOutRow = ", pOutRow,
    sep = ""
    )
#
# CellMap
# Higher values than predicted are colored red
# Lower values than predicted are colored blue
# Colors of values are determined by values in stdResid matrix
#
  if (is.null(showVals)==FALSE){
    if (showVals %in% c("RD", "DR")){
      ggpDDC_D = cellMap(
        D = round(D, roundVal),
        R = round(stdResid, roundVal),
        indcells = indcells,
        indrows = indrows,
        columnlabels = columnlabels,
        rowlabels = rowlabels,
        
        showVals = "D",
        showrows = showrows,
        showcolumns = showcolumns,
        nrowsinblock = nrowsinblock,
        ncolumnsinblock = ncolumnsinblock,
        
        rowtitle = rowtitle,
        columntitle = columntitle,
        columnangle = columnangle,
        sizetitles = sizetitles,
        adjustrowlabels = adjustrowlabels,
        adjustcolumnlabels = adjustcolumnlabels,
        drawCircles = drawCircles,
        
        mTitle = paste(mTitle,": absolute", sep=""),
      )
      ggpDDC_R = cellMap(
        D = round(D, roundVal),
        R = round(stdResid, roundVal),
        indcells = indcells,
        indrows = indrows,
        columnlabels = columnlabels,
        rowlabels = rowlabels,
        
        showVals = "R",
        showrows = showrows,
        showcolumns = showcolumns,
        nrowsinblock = nrowsinblock,
        ncolumnsinblock = ncolumnsinblock,
        
        rowtitle = rowtitle,
        columntitle = columntitle,
        columnangle = columnangle,
        sizetitles = sizetitles,
        adjustrowlabels = adjustrowlabels,
        adjustcolumnlabels = adjustcolumnlabels,
        drawCircles = drawCircles,
        
        mTitle = paste(mTitle,": residuals", sep=""),
      )
    }
  } else {
    ggpDDC = cellMap(
      D = round(D, roundVal),
      R = round(stdResid, roundVal),
      indcells = indcells,
      indrows = indrows,
      columnlabels = columnlabels,
      rowlabels = rowlabels,
      
      showVals = showVals,
      showrows = showrows,
      showcolumns = showcolumns,
      nrowsinblock = nrowsinblock,
      ncolumnsinblock = ncolumnsinblock,
      
      rowtitle = rowtitle,
      columntitle = columntitle,
      columnangle = columnangle,
      sizetitles = sizetitles,
      adjustrowlabels = adjustrowlabels,
      adjustcolumnlabels = adjustcolumnlabels,
      drawCircles = drawCircles,
      
      mTitle = mTitle,
    )
  }
#
# Plot CellMap only if show_cellmap is true
#
  if (showCellMap == TRUE){
    if (is.null(showVals)==FALSE){
      if (showVals %in% c("RD", "DR")){
        grid.arrange(ggpDDC_D, ggpDDC_R, ncol=2)
      }
    } else {
      grid.arrange(ggpDDC, ncol=1)
    }
  }
#
# Creating pdf file if pdf_cellmap is true
#
  if (pdfCellMap == TRUE){
    pdf(
      paste(
        datasetName, " (",used_parameters,").pdf",
        sep = ""
        ),
      width = pdf_width,
      height = pdf_height
      )
    if (is.null(showVals)==FALSE){
      if (showVals %in% c("RD", "DR")){
        gridExtra::grid.arrange(ggpDDC_D, ggpDDC_R, ncol=2)
        }
    } else {
      gridExtra::grid.arrange(ggpDDC, nrow = 1) 
    }
    dev.off()
  }
#
# Return indices of outlying cells and rows from DDC logratio algorithm
#
  return = list(
    indcells,
    indrows
    )
#
  return(return)
}
############################################################################