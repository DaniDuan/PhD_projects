library(flowViz)
data(GvHD)
head(pData(GvHD))
### Only using patient 6 as an example
xyplot(`SSC-H` ~ `FSC-H` | Visit, data=GvHD, subset=Patient=="6")
GvHD = GvHD[pData(GvHD)$Patient == 6]
# transform some of the fluorescence channels to an adequate log-like scale
# asinh function which can deal with negative values
# look in FlowCore manual 7.5 Transformation Filters
tf = transformList(from=colnames(GvHD)[3:7], tfun = asinh) 
GvHD = tf %on% GvHD 
