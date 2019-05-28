data=read.delim("e141_Plant aboveground biomass data.txt",header=T)
PlotSpecies=read.csv2("PlantedSpecies.csv",header=T)

dim(data)
names(data)

sum(data$CountOfSpecies==4)
table(data$Date)
table(data$Date, data$Species)
Year=c()
for(i in 1:(dim(data)[1]))
Year=c(Year, strsplit(as.character(data[i,"Date"]), "/")[[1]][3])
Month=c()
for(i in 1:(dim(data)[1]))
Month=c(Month, strsplit(as.character(data[i,"Date"]), "/")[[1]][1])

table(Year[Month==8], data$Species[Month==8], data$CountOfSpecies[Month==8])

###corrections
table(data$CO2.Treatment)
data[data$CO2.Treatment=="Cenriched ","CO2.Treatment"]="Cenrich"
unique(data$CO2.Treatment)

data$Species=gsub("^\\s+|\\s+$", "", data$Species)
sort(unique(data$Species))
data[data$Species=="bromus inermis","Species"]="Bromus inermis"
data[data$Species=="poa pratensis","Species"]="Poa pratensis"
sort(unique(data$Species))

unique(data[data$CountOfSpecies==0,"Plot"])

#### reformat data
ReFormattedData=c()
ListAll=sort(unique(data$Species))
Bzzz=rep(0, length(ListAll))
names(Bzzz)=ListAll
length(ListAll)

for(samp in unique(data$Sampling..))
for(pl in unique(data$Plot))
{
gz=data[(data$Sampling..==samp)&(data$Plot==pl),]
if(dim(gz)[1]>0)
	{
	for(i in 1:13) if(length(unique(gz[,i]))>1) break()
	Newline=c(gz[1,1:13],Bzzz)
	#names(Newline)
	for (k in 1:(dim(gz)[1]))
	Newline[gz[k,"Species"]]=as.numeric(Newline[gz[k,"Species"]])+gz[k,"Aboveground.Biomass..g.m.2."]
	ReFormattedData=rbind(ReFormattedData,data.frame(Newline),make.row.names = F)
	}
}

dim(ReFormattedData)

RefData=ReFormattedData[,c(1:7,12:13)]

Month=c()
for(i in 1:(dim(RefData)[1]))
Month=c(Month, strsplit(as.character(RefData[i,"Date"]), "/")[[1]][1])

PlanTabl=c()
for(i in 1:(dim(RefData)[1]))
PlanTabl=rbind(PlanTabl, PlotSpecies[PlotSpecies$Plot==RefData[i,"Plot"],3:18])

sum(rowSums(PlanTabl)==RefData[,"CountOfSpecies"])
dim(PlanTabl)

RefData=cbind(RefData,Month,PlanTabl,ReFormattedData[,c(16,      
	18,21:23,25,30,31,51,55:56,72,76,83,89,90)])

Tot=rowSums(matrix(as.numeric(as.matrix(ReFormattedData[,14:98])), nrow=dim(ReFormattedData)[1],
	ncol=dim(ReFormattedData[,14:98])[2]))
Litter=as.numeric(ReFormattedData[,60])+as.numeric(ReFormattedData[,61])
Other=rowSums(matrix(as.numeric(as.matrix(ReFormattedData[,c(14:15,17,19:20,24,26:29,
	32:50,52:54,57:59,62:71,73:75,77:82,84:88,91:98)])), nrow=dim(ReFormattedData)[1],
	ncol=(dim(ReFormattedData[,14:98])[2]-18)))
RefData=cbind(RefData,Other,Litter)

sum(Tot==rowSums(matrix(as.numeric(as.matrix(RefData[,27:44])), nrow=dim(RefData)[1],
	ncol=dim(RefData[,27:44])[2])))
dim(RefData)[1]

which(round(Tot,5)!=round(rowSums(matrix(as.numeric(as.matrix(RefData[,27:44])), nrow=dim(RefData)[1],
	ncol=dim(RefData[,27:44])[2])),5))

#remove plots with no species planted
RefData=RefData[RefData$CountOfSpecies>0,]

#remove water and warming treatment plots
RefData=RefData[RefData$Water.Treatment.!="H2Oneg",]
RefData=RefData[RefData$Temp.Treatment.!="HTelv",]

dim(RefData)
write.csv(RefData, file="BioCONFormatted.csv",row.names = F)