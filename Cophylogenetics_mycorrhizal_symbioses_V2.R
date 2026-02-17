###### Script for analyzing cophylogenetic signal and phylogenetic congruence ######

rm(list=ls())

# def dataset, name_fungi_tree, fungi_outgroup, name_plant_tree, plant_outgroup

# define variables after running bash code using Mafft, Trimal, IqTree
name_fungi_tree=paste0("alignment_",dataset,"_fungi_trimal.fasta.treefile")
name_plant_tree=paste0("alignment_",dataset,"_plant_trimal.fasta.treefile")

setwd(paste0("C:/fantine/Nextcloud/Immex_2023/data/data_",dataset))

library(ape)
library(paco)
library(phytools)


###### Step 1 : Data ######

fungi_tree <- ape::read.tree(name_fungi_tree)
fungi_tree <- root(fungi_tree, fungi_outgroup) #root the tree with outgroup
fungi_tree <- drop.tip(fungi_tree, tip=fungi_outgroup) #delete outgroup
fungi_tree <- ladderize(fungi_tree)
plot(fungi_tree)

plant_tree <- read.tree(name_plant_tree)
plant_tree <- root(plant_tree, plant_outgroup)
plant_tree <- drop.tip(plant_tree, tip=plant_outgroup)
plant_tree <- ladderize(plant_tree)
plot(plant_tree)

network<-read.table(name_network, sep=";", header=T) # !! network must have plants as rows and fungi as columns, network<-t(network) if it is not the case

!plant_tree$tip.label %in% rownames(network) #check if plants are in network : true if they aren't
!fungi_tree$tip.label %in% colnames(network)


plant_tree$tip.label[!plant_tree$tip.label %in% rownames(network)] #list of plants which are not in network
plant_tree <- drop.tip(plant_tree, tip=plant_tree$tip.label[!plant_tree$tip.label %in% rownames(network)]) #tree with only plants of network
plant_tree <- ladderize(plant_tree) 
plot(plant_tree) 

# idem for fungi
!fungi_tree$tip.label %in% colnames(network)
fungi_tree$tip.label[!fungi_tree$tip.label %in% colnames(network)]
fungi_tree <- drop.tip(fungi_tree, tip=fungi_tree$tip.label[!fungi_tree$tip.label %in% colnames(network)])
fungi_tree <- ladderize(fungi_tree)
plot(fungi_tree)

network <- network[plant_tree$tip.label,fungi_tree$tip.label]

##### Step 2 : Parafit and Paco #######
res_parafit <- parafit(cophenetic(plant_tree),cophenetic(fungi_tree),network,nperm=10000, correction="lingoes")$p.global
#cophenetic(plant_tree) pour convertir arbre en matrice des distances

D<-prepare_paco_data(cophenetic(plant_tree),cophenetic(fungi_tree),network)
D <- add_pcoord(D)
data_paco <- PACo(D,nperm=10000, symmetric=TRUE)
res_paco <- c(data_paco$gof$p, 1- data_paco$gof$ss)
# gof$p is pvalue, 1-(gof$ss) is Rsquare

res<- data.frame()
res <- rbind(res,c(res_parafit,res_paco))
colnames(res) <- c("p-value Parafit","p-value Paco","Rcarré Paco")
write.table(res,paste0("C:/fantine/Nextcloud/Immex_2023/results/GlobalFit_",dataset,".csv"),sep=";",row.names=FALSE,quote=FALSE)

##### Step 3 : before empress #####

# Problem: polytomies are not supported by eMPRess
## Approach n°1 : pick one interaction at random or select the most abundant interaction
network <- t(network)
list_links <- reshape2::melt(as.matrix(network)) # Var1 corresponds to symbionts, Var2 to hosts
list_links <- list_links[list_links$value>0,] # remove no interaction

# pick one interaction at random, and set 0 to others
for (symbiont in unique(list_links$Var1)){
  list_links$value[which(list_links$Var1==symbiont)] <- 0
  list_links$value[as.numeric(sample(size=1, as.character(which(list_links$Var1==symbiont))))] <- 1
}

list_links <- list_links[list_links$value>0,] # remove no interaction


#creation of table of links
write.table(paste0(list_links$Var1,":", list_links$Var2), paste0("links_empress_", dataset,".txt"), col.names=F, row.names=F, quote=F)

# remove the tips with no interactions and creation of new trees
tree_hosts_interact <- drop.tip(plant_tree, tip=plant_tree$tip.label[which(!plant_tree$tip.label %in% list_links$Var2)])
tree_hosts_interact$node.label <- NULL
tree_hosts_interact <- multi2di(tree_hosts_interact)
write.tree(tree_hosts_interact, paste0("tree_hosts_empress_", dataset,".tre"))
tree_parasites_interact <- drop.tip(fungi_tree, tip=fungi_tree$tip.label[which(!fungi_tree$tip.label %in% list_links$Var1)])
tree_parasites_interact$node.label <- NULL
tree_parasites_interact <- multi2di(tree_parasites_interact)
write.tree(tree_parasites_interact, paste0("tree_parasites_empress_", dataset,".tre"))


## Approach n°2 : randomly simulate bifurcating sub-trees

# network <- t(network) # we need fungi in columns and plants in rows
position <- min(c(0.001, fungi_tree$edge.length))

num_duplicate <- 0
while (!all(colSums(network)==1)){
  ind_species <- which(colSums(network)>1)[sample(length(which(colSums(network)>1)), size=1)]
  species <- colnames(network)[ind_species]
  num_duplicate <- num_duplicate + 1
  
  network <- cbind(network, rep(0, nrow(network)))
  colnames(network)[ncol(network)] <- paste0("duplicate_",num_duplicate)
  list_index <- which(network[,species]>0)
  index <- list_index[sample(length(list_index),size=1)]
  network[index, ncol(network)] <- 1
  network[index, species] <- 0
  
  fungi_tree <- bind.tip(fungi_tree, tip.label=paste0("duplicate_",num_duplicate), edge.length=0.001, where=which(fungi_tree$tip.label==species), position=position)
  
}

if (!is.ultrametric(fungi_tree)) fungi_tree <- force.ultrametric(fungi_tree, method="extend")

list_links <- reshape2::melt(as.matrix(network))
list_links <- list_links[list_links$value>0,]

write.table(paste0(list_links$Var2,":", list_links$Var1), paste0("links_bifurcations_empress_", dataset, "_", seed, ".txt"), col.names=F, row.names=F, quote=F)
fungi_tree$node.label <- NULL
fungi_tree <- multi2di(fungi_tree)
write.tree(fungi_tree, paste0("tree_parasites_bifurcations_empress_", dataset, "_", seed, ".tre"))
plant_tree$node.label <- NULL
plant_tree <- multi2di(plant_tree)
write.tree(plant_tree, paste0("tree_hosts_bifurcations_empress_", dataset, "_", seed, ".tre"))

##### Step 4 : empress #####
#cf code bash empress

##### Step 5 : using empress data #####
bifurcation <- TRUE
suffix <- if (bifurcation) "_bifurcations" else ""
empress=data.frame()

parameters<-data.frame()
parameters<-rbind(c(1,1,1),c(4,1,1),c(2,1,2),c(4,2,1),c(2,3,1))
colnames(parameters)<-c("d","t","l")

for (i in 1:5) {
  d=parameters[i,1]
  t=parameters[i,2]
  l=parameters[i,3]
  costs=paste0("d",d,"_t",t,"_l",l)
  
  reconciliation <- read.table(paste0("C:/fantine/Nextcloud/Immex_2023/data/data_",dataset,"/recon_",dataset,"_",costs,suffix,"_output.csv"),sep=",")
  cospeciation <- length(reconciliation$V3[reconciliation$V3=="Cospeciation"])/length(reconciliation$V3[reconciliation$V3=="Transfer"]) #indicate if we have more co-speciation events than transfer events
  
  graph <- read.table(paste0("C:/fantine/Nextcloud/Immex_2023/data/data_",dataset,"/pvalue_",dataset,"_",costs,suffix,"_output.svg"), comment.char = "", fill=TRUE, sep=";")
  pvalue_empress <- graph$V1[grep("p-value", graph$V1)] #isolate the row of the document with the p-value
  pvalue_empress <- gsub("    <!-- p-value = ", "", pvalue_empress)
  pvalue_empress <- as.numeric(gsub(" -->", "", pvalue_empress)) #isolate the p-value
  
  
  empress <- rbind(empress,c(costs,cospeciation ,pvalue_empress))
  colnames(empress) <- c("Costs","Ratio cospeciation / transfer events","P-value")
}

write.table(empress,paste0("C:/fantine/Nextcloud/Immex_2023/results/Empress_",dataset,suffix,".csv"),sep=";",row.names=FALSE,quote=FALSE)


##### Step 6 : Representation of network ######

network <- t(network)
melt_network <- reshape2::melt(as.matrix(network))
melt_network <- melt_network[which(melt_network$value>0),]

library(igraph)

g <- graph_from_edgelist(as.matrix(melt_network[,1:2]), directed = FALSE)

edge.attributes(g)$weight <- melt_network[,3]
E(g)$weight # same as edge.attributes(g)$weight

V(g)$type <- rep("Fungi", length(V(g))) # same as vertex.attributes(g)$type

## Merckx
V(g)$type[names(V(g)) %in% c("EU420988.1")] <- "A. foertheriana"
V(g)$type[names(V(g)) %in% c("EU420989.1")] <- "A. gesnerioides"
V(g)$type[names(V(g)) %in% c( "EU420990.1")] <- "A. hydra"
V(g)$type[names(V(g)) %in% c("EU420991.1")] <- "A. korupensis"
V(g)$type[names(V(g)) %in% c("EU420992.1")] <- "A. winkleri"


# define color and shape mappings.
shape <- c(rep("circle", 5), "square")
size <- c(rep(10, 5), 5)
col <- c( "#f4d03f", "#5499c7", "#16a085","#ec7063","#bb8fce","#616a6b")
names(size) <- names(shape) <- names(col) <- c("A. foertheriana", "A. gesnerioides", "A. hydra", "A. korupensis", "A. winkleri", "Fungi")

## Zhao

V(g)$type[names(V(g))=="Burmannia_championii"] <- "B. championii"
V(g)$type[names(V(g))=="Burmannia_longifolia"] <- "B. longifolia"
V(g)$type[names(V(g))=="Burmannia_wallichii"] <- "B. wallichii"
V(g)$type[names(V(g))=="Burmannia_itoana"] <- "B. itoana"
V(g)$type[names(V(g))=="Burmannia_disticha"] <- "B. disticha"
V(g)$type[names(V(g))=="Burmannia_filamentosa"] <- "B. filamentosa"
V(g)$type[names(V(g))=="Burmannia_coelestis"] <- "B. coelestis"
V(g)$type[names(V(g))=="Burmannia_cryptopetala"] <- "B. cryptopetala"
V(g)$type[names(V(g))=="Burmannia_oblonga"] <- "B. oblonga"
V(g)$type[names(V(g))=="Burmannia_nepalensis"] <- "B. nepalensis"


# define color and shape mappings.
shape <- c(rep("circle", 10), "square")
size <- c(rep(10, 10), 5)
col <- c( "#f4d03f", "#85c1e9", "#148f77","#ec7063","#9b59b6","#7e5109","#58d68d","#e6b0aa","#1f618d","#f39c12","#616a6b")
names(size) <- names(shape) <- names(col) <- c("B. championii","B. longifolia","B. wallichii","B. itoana","B. disticha","B. filamentosa","B. coelestis","B. cryptopetala","B. oblonga","B. nepalensis","Fungi")

## Hayward

V(g)$type[names(V(g))=="Guapira_fragrans"] <- "Guapira_fragrans"
V(g)$type[names(V(g))== "Pisonia_grandis" ] <-  "Pisonia_grandis" 
V(g)$type[names(V(g))=="Pisonia_albida"] <- "Pisonia_albida"
V(g)$type[names(V(g))=="Pisonia_taina"] <- "Pisonia_taina"
V(g)$type[names(V(g))=="Pisonia_subcordata"] <- "Pisonia_subcordata"
V(g)$type[names(V(g))=="Pisonia_rotundata"] <- "Pisonia_rotundata"
V(g)$type[names(V(g))=="Pisonia_aculeata"] <- "Pisonia_aculeata"
V(g)$type[names(V(g))=="Pisonia_sechellarum"] <- "Pisonia_sechellarum"
V(g)$type[names(V(g))=="Pisonia_umbellifera"] <- "Pisonia_umbellifera"
V(g)$type[names(V(g))=="Pisonia_brunoniana"] <- "Pisonia_brunoniana"
V(g)$type[names(V(g))=="Pisonia_sandwicensis"] <- "Pisonia_sandwicensis"
V(g)$type[names(V(g))=="Neea_buxifolia"] <- "Neea_buxifolia"
V(g)$type[names(V(g))=="Guapira_discolor"] <- "Guapira_discolor"


# define color and shape mappings.
shape <- c(rep("circle", 13), "square")
size <- c(rep(10, 13), 5)
col <- c( "#f4d03f", "#85c1e9", "#148f77","#ec7063","#2e86c1","#f5cba7","#9b59b6","#7e5109","#58d68d","#e6b0aa","#1f618d","#f39c12","#DAF7A6","#616a6b")
names(size) <- names(shape) <- names(col) <- c("Guapira_fragrans", "Pisonia_grandis" ,"Pisonia_albida","Pisonia_taina","Pisonia_subcordata","Pisonia_rotundata","Pisonia_aculeata","Pisonia_sechellarum","Pisonia_umbellifera","Pisonia_brunoniana","Pisonia_sandwicensis","Neea_buxifolia","Guapira_discolor","Fungi")

#Toju
V(g)$type[names(V(g))=="Arcterica_nana"] <- "Arcterica_nana"
V(g)$type[names(V(g))== "Arctous_alpina" ] <-  "Arctous_alpina"
V(g)$type[names(V(g))=="Empetrum_nigrum"] <- "Empetrum_nigrum"
V(g)$type[names(V(g))=="Harrimanella_stelleriana"] <- "Harrimanella_stelleriana"
V(g)$type[names(V(g))=="Loiseleuria_procumbens"] <- "Loiseleuria_procumbens"
V(g)$type[names(V(g))=="Gaultheria_pyroloides"] <- "Gaultheria_pyroloides"
V(g)$type[names(V(g))=="Phyllodoce_nipponica"] <- "Phyllodoce_nipponica"
V(g)$type[names(V(g))=="Cladothamnus_bracteatus"] <- "Cladothamnus_bracteatus"
V(g)$type[names(V(g))=="Rhododendron_aureum"] <- "Rhododendron_aureum"
V(g)$type[names(V(g))=="Vaccinium_ovalifolium"] <- "Vaccinium_ovalifolium"
V(g)$type[names(V(g))=="Vaccinium_vitis-idaea"] <- "Vaccinium_vitis-idaea"
V(g)$type[names(V(g))=="Pinus_pumila"] <- "Pinus_pumila"
V(g)$type[names(V(g))=="Vaccinium_uliginosum"] <- "Vaccinium_uliginosum"
V(g)$type[names(V(g))=="Vaccinium_hirtum"] <- "Vaccinium_hirtum"
V(g)$type[names(V(g))=="Vaccinium_smallii"] <- "Vaccinium_smallii"
V(g)$type[names(V(g))=="Phyllodoce_aleutica"] <- "Phyllodoce_aleutica"

# define color and shape mappings.
shape <- c(rep("circle", 16), "square")
size <- c(rep(10, 16), 5)
col <- c( "#f4d03f", "#85c1e9", "#148f77","#ec7063","#2e86c1","#f5cba7","#9b59b6","#7e5109","#58d68d","#e6b0aa","#1f618d","#f39c12", "#daf7a6","#a93226","#1b2631","#c39bd3","#616a6b")
names(size) <- names(shape) <- names(col) <- c("Arcterica_nana", "Arctous_alpina" ,"Empetrum_nigrum","Gaultheria_pyroloides","Harrimanella_stelleriana","Loiseleuria_procumbens","Phyllodoce_nipponica","Cladothamnus_bracteatus","Rhododendron_aureum","Vaccinium_ovalifolium","Vaccinium_vitis-idaea","Pinus_pumila","Vaccinium_uliginosum","Vaccinium_hirtum","Vaccinium_smallii","Phyllodoce_aleutica","Fungi")

# Opik
V(g)$type[grep(" ",names(V(g)))] <- "Plant"
shape <- c("circle","square")
size <- c(10,5)
col <- c("#4EA72E","#83625B")
names(size) <- names(shape) <- names(col) <- c("Plant","Fungi")

## create pdf
pdf(paste0("network_graph_",dataset,".pdf"))
set.seed(1)
plot(g, vertex.color = col[V(g)$type], vertex.shape = shape[V(g)$type], vertex.label=NA, vertex.size=size[V(g)$type],
     edge.width=sqrt(E(g)$weight),  layout= layout_with_fr(g, weights = E(g)$weight,niter = 5000 ) )
plot(1,1,col="white")
legend("bottomleft", pch=19,names(col),col=col,cex=0.8)
dev.off()


# 2 phylogenetic trees with links
library(phytools)
network <- t(network)
list_links <- reshape2::melt(as.matrix(network)) # Var1 corresponds to symbionts, Var2 to hosts
list_links <- list_links[list_links$value>0,] # remove no interaction
list_links <- list_links[,-3]
list_links <- as.matrix(list_links)
list_links[,"Var1"] <- as.character(list_links[,"Var1"])
list_links[,"Var1"] <- gsub(" ","",list_links[,"Var1"])
list_links[,"Var1"] <- paste0("VT",list_links[,"Var1"])
fungi_tree$tip.label <- paste0("VT",fungi_tree$tip.label)

#cophyloplot(plant_tree,fungi_tree,list_links[which(list_links[,1] %in% fungi_tree$tip.label),])
trees_network<-phytools::cophylo(fungi_tree,plant_tree,assoc=list_links,link.type="curved")
png("trees_Opik.png",width=14000,height=6000)
plot(trees_network,link.type="straight",link.lty="solid",
     fsize=c(13,17),link.lwd=3)
dev.off()
png("plant_tree_Opik.png",width=3000,height=4000)
plot(trees_network$trees[[2]],show.tip.label = F,edge.color = "#4EA72E",edge.width=20)
dev.off()
png("fungi_tree_Opik.png",width=3000,height=4000)
plot(trees_network$trees[[1]],show.tip.label = F,edge.color = "#83625B",edge.width=20)
dev.off()

###### Script Merckx 2008 ######
dataset="Merckx"
network <- data.frame(cbind(c(1,1,1,0,0,0,0,0,0,0,0,0,0,0,0),c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0),c(0,0,0,0,1,1,1,1,0,0,0,0,0,0,0),c(0,0,0,0,0,0,0,0,1,1,0,0,0,0,0),c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,0)))
rownames(network) <- c("EU417584.1","EU417588.1","EU417592.1","EU417580.1","EU417620.1","EU417624.1","EU417628.1","EU417632.1","EU417596.1","EU417600.1","EU417604.1","EU417608.1","EU417612.1","EU417616.1","AY635831.1")
colnames(network) <- c("EU420988.1","EU420989.1","EU420990.1","EU420991.1","EU420992.1")
fungi_outgroup="Z14007.1"
plant_outgroup="DQ786081.1"
name_fungi_tree=paste0("alignment_",dataset,"_Glomeromycotina_trimal.fasta.treefile")
name_plant_tree=paste0("alignment_",dataset,"_Afrothismia_trimal.fasta.treefile")

##### Script Perez 2020 #####
dataset="global_scale"
name_fungi_tree="tree_VT_18S_barcode_PB_LN_GTR_constrained.tre"
name_plant_tree="plant_tree_phylomatic_zanne_grafted_familly.tre"
fungi_tree <- read.nexus(name_fungi_tree)
name_network="network_all.csv"
# trees have roots already
# too much time for a simple computer

##### Script Martos #####
library("RPANDA")
data(mycorrhizal_network)
dataset="Martos"
network <- mycorrhizal_network[[1]] # interaction matrix with orchids in columns and fungi in rows
network<-t(network)
plant_tree <- mycorrhizal_network[[2]] # phylogenetic tree (phylo object)
fungi_tree <- mycorrhizal_network[[3]] # phylogenetic tree (phylo object)
#all plants and fungi are in network, and trees have roots already

### Script Barret ####
dataset="Barret"

# put together data of the 2 genes for plants
gene1<-read.table("AppendixS4_1.interleaved.txt")
gene2<-read.table("AppendixS4_2.interleaved.txt")
gene1and2<-gene1
gene1and2$V2<-paste0(gene1and2$V2,gene2$V2)
write.table(gene1and2,"AppendixS4.interleaved.txt",quote=FALSE,row.names=FALSE,sep="\t")

# transformation nexus format -> fasta format for fungi
alignment_fungi<- read.dna("AppendixS3.interleaved.txt")
write.dna(alignment_fungi,paste0("alignment_",dataset,"_fungi.fasta"),format="fasta")

#defining variables
name_fungi_tree="alignment_Barret_fungi.fasta.treefile"
name_plant_tree="alignment_Barret_plant.fasta.treefile"
fungi_outgroup="UDB000216_Thel._terrestris"
plant_outgroup="Aplectrum"

#building network
plant<-plant_tree$tip.label
fungi<-fungi_tree$tip.label
plant<-sort(plant)
fungi<-sort(fungi)
network<-data.frame(diag(107))
rownames(network) <- plant
colnames(network) <- fungi

##### Script Dowie #######
dataset="Dowie"
library(rentrez)

# import data
data<-read.table(paste0("data_",dataset,".csv"),sep=";")

#creation of fasta file for plants and fungi, and creation of network
for (i in 1:70){
  genbank_id <- data$V3[i]
  
  genbank_record <- entrez_fetch(db="nucleotide", id=genbank_id, rettype="fasta", retmode="text")
  
  write.table(genbank_record,file=paste0("sequences_",dataset,"_plant.fasta"),append=TRUE,quote=FALSE,row.names = FALSE,col.names=FALSE,eol="")
}

for (i in 1:70){
  genbank_id <- data$V6[i]
  
  genbank_record <- entrez_fetch(db="nucleotide", id=genbank_id, rettype="fasta", retmode="text")
  
  write.table(genbank_record,file=paste0("sequences_",dataset,"_fungi.fasta"),append=TRUE,quote=FALSE,row.names = FALSE,col.names=FALSE,eol="")
}

#define variables (after running bash code using Mafft, Trimal, IqTree)
name_fungi_tree=paste0("alignment_",dataset,"_fungi_trimal.fasta.treefile")
name_plant_tree=paste0("alignment_",dataset,"_plant_trimal.fasta.treefile")


for (i in 1:70){
  data$V3[i]<-paste0(data$V3[i],".1")
  data$V6[i]<-paste0(data$V6[i],".1")
}

network<-data.frame(diag(70))
rownames(network) <- data$V3
colnames(network) <- data$V6

#midpoint rooting the trees
library(phangorn)
fungi_tree<-midpoint(fungi_tree)
fungi_tree <- ladderize(fungi_tree)
plot(fungi_tree)
plant_tree<-midpoint(plant_tree)
plant_tree <- ladderize(plant_tree)
plot(plant_tree)

##### Script Otero #####
dataset="Otero"

library(rentrez)

# import data
data<-read.table(paste0("data_",dataset,".csv"),sep=";")

#creation of fasta file for plants and fungi (think about delete the 2 duplicate sequences in plant fasta file)
for (i in 1:34){
  genbank_id <- data$V6[i]
  genbank_record <- entrez_fetch(db="nucleotide", id=genbank_id, rettype="gb", retmode="text")
  sequence <- strsplit(split="\n        1 ", genbank_record)[[1]][2]
  sequence <- paste0(sequence, collapse = "")
  sequence <- gsub("1|2|3|4|5|6|7|8|9|0", "", sequence)
  sequence <- gsub("\n", "", sequence, fixed=TRUE)
  sequence <- gsub("//", "", sequence, fixed=TRUE)
  sequence <- gsub(" ", "", sequence, fixed=TRUE)
  write.table(paste0(">",genbank_id,"\n",sequence,"-"),file=paste0("sequences_",dataset,"_plant.fasta"),append=TRUE,quote=FALSE,row.names = FALSE,col.names=FALSE,eol="")
  
  genbank_id <- data$V7[i]
  genbank_record <- entrez_fetch(db="nucleotide", id=genbank_id, rettype="gb", retmode="text")
  sequence <- strsplit(split="\n        1 ", genbank_record)[[1]][2]
  sequence <- paste0(sequence, collapse = "")
  sequence <- gsub("1|2|3|4|5|6|7|8|9|0", "", sequence)
  sequence <- gsub("\n", "", sequence, fixed=TRUE)
  sequence <- gsub("//", "", sequence, fixed=TRUE)
  sequence <- gsub(" ", "", sequence, fixed=TRUE)
  write.table(paste0(sequence,"\n"),file=paste0("sequences_",dataset,"_plant.fasta"),append=TRUE,quote=FALSE,row.names = FALSE,col.names=FALSE,eol="")
}

for (i in 1:34){
  genbank_id <- data$V8[i]
  genbank_record <- entrez_fetch(db="nucleotide", id=genbank_id, rettype="gb", retmode="text")
  sequence <- strsplit(split="\n        1 ", genbank_record)[[1]][2]
  sequence <- paste0(sequence, collapse = "")
  sequence <- gsub("1|2|3|4|5|6|7|8|9|0", "", sequence)
  sequence <- gsub("\n", "", sequence, fixed=TRUE)
  sequence <- gsub("//", "", sequence, fixed=TRUE)
  sequence <- gsub(" ", "", sequence, fixed=TRUE)
  write.table(paste0(">",genbank_id,"\n",sequence,"-"),file=paste0("sequences_",dataset,"_fungi.fasta"),append=TRUE,quote=FALSE,row.names = FALSE,col.names=FALSE,eol="")
  
  genbank_id <- data$V9[i]
  genbank_record <- entrez_fetch(db="nucleotide", id=genbank_id, rettype="gb", retmode="text")
  sequence <- strsplit(split="\n        1 ", genbank_record)[[1]][2]
  sequence <- paste0(sequence, collapse = "")
  sequence <- gsub("1|2|3|4|5|6|7|8|9|0", "", sequence)
  sequence <- gsub("\n", "", sequence, fixed=TRUE)
  sequence <- gsub("//", "", sequence, fixed=TRUE)
  sequence <- gsub(" ", "", sequence, fixed=TRUE)
  write.table(paste0(sequence,"\n"),file=paste0("sequences_",dataset,"_fungi.fasta"),append=TRUE,quote=FALSE,row.names = FALSE,col.names=FALSE,eol="")
}

#creation of network
network<-data.frame(diag(31)) # there are 2 plants (‘GQ405650’, ‘GQ405652’) which have 2 fungi, and the 33 others have 1 fungi. The last plant is extra-group.
rownames(network) <- data$V6[c(1:31)]
colnames(network) <- data$V8[c(1:31)]
network$GQ405558 <- rep(0,times=31)
network$DQ028814 <- rep(0,times=31)
network["GQ405650","GQ405558"] <- 1
network["GQ405652","DQ028814"] <- 1

#define variables (after running bash code using Mafft, Trimal, IqTree)
name_fungi_tree=paste0("alignment_",dataset,"_fungi_trimal.fasta.treefile")
name_plant_tree=paste0("alignment_",dataset,"_plant_trimal.fasta.treefile")
fungi_outgroup="GQ405535"
plant_outgroup="GQ405626"


##### Script Zhao #####
dataset="Zhao"
name_fungi_tree="OTU-tree-bootstrap1000.tre"
name_plant_tree="tree-plant Export.nwk"
plant_outgroup="Aletris_lutea"

# !!two outgroups for fungi !!
fungi_tree <- ape::read.tree(name_fungi_tree)
fungi_tree <- root(fungi_tree, node = getMRCA(tip=c("Otu000085","Otu000093"), fungi_tree))
fungi_tree <- drop.tip(fungi_tree, tip=c("Otu000085","Otu000093"))
fungi_tree <- ladderize(fungi_tree)
plot(fungi_tree)

# Building of network
abundance <- read.table("S3.csv",sep=";",row.names=1)

# 1. attribute a species to each sample
species=abundance[["Species","V2"]]
for (i in 2:length(abundance["Species",])+1) {
  if (abundance[["Species",paste0("V",i)]]=="") {
    abundance[["Species",paste0("V",i)]]=species
  }  
  else {
    species=abundance[["Species",paste0("V",i)]]
  }
}

# 2. building of a new dataframe of mean abundance for each species
abundance<-data.frame(t(abundance))
abundance_species<-data.frame(matrix(0,length(unique(abundance$Species)),ncol(abundance)-3))
row.names(abundance_species)<-unique(abundance$Species)
colnames(abundance_species)<-colnames(abundance)[grep("Otu",colnames(abundance))]

for (i in rownames(abundance_species)) {
  abundance_one_species<-abundance[abundance$Species==i,]
  for (j in colnames(abundance_species)) {
    abundance_species[[i,j]]<-mean(as.numeric(abundance_one_species[,j]))
  }
  
}

#rename columns and row as in tree files
library(dplyr)
abundance_species <- abundance_species %>% rename_with(~sub("Otu", "Otu00", .x), starts_with("Otu"))  
rownames(abundance_species)<-sub(pattern = "B. ", replacement="Burmannia_",x=rownames(abundance_species))

# 3.building of the network
prevalence<-abundance_species
library(ggplot2)
#ggplot(abundance_species,aes(x=Otu0002))+geom_density()
threshold=0.05
prevalence[prevalence>threshold]<-1
prevalence[prevalence<threshold]<-0
sum(rowSums(prevalence)!=0)
sum(colSums(prevalence)!=0)
network<-prevalence[colSums(prevalence)!=0]

colnames(network) %in% fungi_tree$tip.label

# Just before preparing empress : for each fungi, select the plant where the fungi is the most abundant
most_abundant_species <- abundance_species
for (j in colnames(most_abundant_species)) {
  for (i in rownames(most_abundant_species)) {
    if (most_abundant_species[[i,j]]<max(most_abundant_species[,j])){
      most_abundant_species[[i,j]]=0
    }
  }
}

prevalence_max<-most_abundant_species
prevalence_max[prevalence_max>threshold]<-1
prevalence_max[prevalence_max<threshold]<-0
network<-prevalence_max[colSums(prevalence)!=0]
network <- network[plant_tree$tip.label,fungi_tree$tip.label]

      


##### Script Arifin #####
dataset="Arifin"
name_network="Network6.csv"
network<-read.table(name_network, sep=",", header=T,row.names=1)
name_plant_tree="plant_tree.treefile" # already rooted
name_fungi_tree="alignment_Arifin_fungi_trimal.fasta.treefile"
fungi_outgroup="OTU_4" ## !!! do not delete it, it is not really the outgroup !!!

## A special Step 3 before empress : Beware, here there are more plants than fungi so we chose to consider plants as symbionts and fungi as host !
## Approach n°1 : pick one interaction at random
list_links <- reshape2::melt(as.matrix(network)) # Var1 corresponds to symbionts, Var2 to hosts
list_links <- list_links[list_links$value>0,] # remove no interaction

# pick one interaction at random, and set 0 to others.
for (symbiont in unique(list_links$Var1)){
  list_links$value[which(list_links$Var1==symbiont)] <- 0
  list_links$value[as.numeric(sample(size=1, as.character(which(list_links$Var1==symbiont))))] <- 1
}
list_links <- list_links[list_links$value>0,] # remove no interaction

# remove the tips with no interactions and creation of new trees
tree_hosts_interact <- drop.tip(fungi_tree, tip=fungi_tree$tip.label[which(!fungi_tree$tip.label %in% list_links$Var2)])
tree_hosts_interact$node.label <- NULL
tree_hosts_interact <- multi2di(tree_hosts_interact)
write.tree(tree_hosts_interact, paste0("tree_hosts_empress_", dataset,".tre"))
tree_parasites_interact <- drop.tip(plant_tree, tip=plant_tree$tip.label[which(!plant_tree$tip.label %in% list_links$Var1)])
tree_parasites_interact$node.label <- NULL
tree_parasites_interact <- multi2di(tree_parasites_interact)
write.tree(tree_parasites_interact, paste0("tree_parasites_empress_", dataset,".tre"))


## Approach n°2 : randomly simulate bifurcating sub-trees
seed<-1
network <- t(network) # we need plant in columns and fungi in rows
position <- min(c(0.001, plant_tree$edge.length))

num_duplicate <- 0
while (!all(colSums(network)==1)){
  ind_species <- which(colSums(network)>1)[sample(length(which(colSums(network)>1)), size=1)]
  species <- colnames(network)[ind_species]
  num_duplicate <- num_duplicate + 1
  
  network <- cbind(network, rep(0, nrow(network)))
  colnames(network)[ncol(network)] <- paste0("duplicate_",num_duplicate)
  list_index <- which(network[,species]>0)
  index <- list_index[sample(length(list_index),size=1)]
  network[index, ncol(network)] <- 1
  network[index, species] <- 0
  
  plant_tree <- bind.tip(plant_tree, tip.label=paste0("duplicate_",num_duplicate), edge.length=0.001, where=which(plant_tree$tip.label==species), position=position)
  
}

if (!is.ultrametric(plant_tree)) plant_tree <- force.ultrametric(plant_tree, method="extend")

list_links <- reshape2::melt(as.matrix(network))
list_links <- list_links[list_links$value>0,]

write.table(paste0(list_links$Var2,":", list_links$Var1), paste0("links_bifurcations_empress_", dataset, "_", seed, ".txt"), col.names=F, row.names=F, quote=F)
fungi_tree$node.label <- NULL
fungi_tree <- multi2di(fungi_tree)
write.tree(fungi_tree, paste0("tree_hosts_bifurcations_empress_", dataset, "_", seed, ".tre"))
plant_tree$node.label <- NULL
plant_tree <- multi2di(plant_tree)
write.tree(plant_tree, paste0("tree_parasites_bifurcations_empress_", dataset, "_", seed, ".tre"))


##### Script Hayward #####
dataset="Hayward"

fungi_accessions<-read.table("genbank_accessions_fungi.csv",sep=";",head=T)
for (i in 1:length(fungi_accessions$Accession.number)){
  fungi_accessions$Accession.number[i]<-paste0(fungi_accessions$Accession.number[i],".1")
}

# building of fasta files for fungi
library(rentrez)
for (i in 1:96){
  genbank_id <- fungi_accessions$Accession.number[i]
  
  genbank_record <- entrez_fetch(db="nucleotide", id=genbank_id, rettype="fasta", retmode="text")
  
  write.table(genbank_record,file=paste0("sequences_",dataset,"_fungi.fasta"),append=TRUE,quote=FALSE,row.names = FALSE,col.names=FALSE,eol="")
}

# building of fasta files for plants with ITS and trnL genes
plant_accessions<-read.table("genbank_accessions_plant.csv",sep=";",head=T)
plant_accessions$Plant<-lapply(plant_accessions$Plant,function(x) gsub(" ","_",x))
for (i in 1:14){
  genbank_id <- plant_accessions$ITS[i]
  genbank_record <- entrez_fetch(db="nucleotide", id=genbank_id, rettype="gb", retmode="text")
  sequence <- strsplit(split="\n        1 ", genbank_record)[[1]][2]
  sequence <- paste0(sequence, collapse = "")
  sequence <- gsub("1|2|3|4|5|6|7|8|9|0", "", sequence)
  sequence <- gsub("\n", "", sequence, fixed=TRUE)
  sequence <- gsub("//", "", sequence, fixed=TRUE)
  sequence <- gsub(" ", "", sequence, fixed=TRUE)
  write.table(paste0(">",plant_accessions$Plant[i],"\n",sequence,"-"),file=paste0("sequences_",dataset,"_plant.fasta"),append=TRUE,quote=FALSE,row.names = FALSE,col.names=FALSE,eol="")
  
  genbank_id <- plant_accessions$trnL[i]
  genbank_record <- entrez_fetch(db="nucleotide", id=genbank_id, rettype="gb", retmode="text")
  sequence <- strsplit(split="\n        1 ", genbank_record)[[1]][2]
  sequence <- paste0(sequence, collapse = "")
  sequence <- gsub("1|2|3|4|5|6|7|8|9|0", "", sequence)
  sequence <- gsub("\n", "", sequence, fixed=TRUE)
  sequence <- gsub("//", "", sequence, fixed=TRUE)
  sequence <- gsub(" ", "", sequence, fixed=TRUE)
  write.table(paste0(sequence,"\n"),file=paste0("sequences_",dataset,"_plant.fasta"),append=TRUE,quote=FALSE,row.names = FALSE,col.names=FALSE,eol="")
}

# Mafft, Trimal,IqTree
plant_outgroup="Bougainvillea_sp."

# height outgroups for fungi
fungi_tree <- ape::read.tree(name_fungi_tree)
fungi_tree <- root(fungi_tree, node = getMRCA(tip=c("KF836022.1","KF836023.1","KF836024.1","KF836025.1","KF836026.1","KF836027.1","KF836028.1","KF836029.1"), fungi_tree))
fungi_tree <- drop.tip(fungi_tree, tip=c("KF836022.1","KF836023.1","KF836024.1","KF836025.1","KF836026.1","KF836027.1","KF836028.1","KF836029.1"))
fungi_tree <- ladderize(fungi_tree)
plot(fungi_tree)

# building of network
fungi_accessions$Host.plant<-lapply(fungi_accessions$Host.plant,function(x) gsub(" ","_",x))
fungi_accessions$Host.plant<-lapply(fungi_accessions$Host.plant,function(x) gsub("Neea_sp.","Neea_buxifolia",x))
fungi_accessions$Host.plant<-lapply(fungi_accessions$Host.plant,function(x) gsub("Guapira_sp.","Guapira_discolor",x))
network<-data.frame(matrix(0,length(plant_tree$tip.label),length(fungi_tree$tip.label)),row.names = plant_tree$tip.label)
colnames(network)=fungi_tree$tip.label
for (i in 1:(length(fungi_accessions$Accession.number)-8)) {
  network[[fungi_accessions[[i,"Host.plant"]],fungi_accessions[[i,"Accession.number"]]]]<-1
}

##### Script Shefferson #####
dataset="Shefferson"

library(rentrez)
# building of a file for plants

for (i in 988:999){
  genbank_id <- paste0("HM140",i)
  
  genbank_record <- entrez_fetch(db="nucleotide", id=genbank_id, rettype="fasta", retmode="text")
  
  write.table(paste0(strsplit(genbank_record,",")[[1]][1],"\n"),file=paste0("GenbankAccessions_",dataset,"_plant.csv"),append=TRUE,quote=FALSE,row.names = FALSE,col.names=FALSE,eol="")
}
for (i in 1000:1006){
  genbank_id <- paste0("HM14",i)
  
  genbank_record <- entrez_fetch(db="nucleotide", id=genbank_id, rettype="fasta", retmode="text")
  
  write.table(paste0(strsplit(genbank_record,",")[[1]][1],"\n"),file=paste0("GenbankAccessions_",dataset,"_plant.csv"),append=TRUE,quote=FALSE,row.names = FALSE,col.names=FALSE,eol="")
}
for (i in c("HM141077","HM151402","AJ539519","AF366896")){
  genbank_record <- entrez_fetch(db="nucleotide", id=i, rettype="fasta", retmode="text")
  
  write.table(paste0(strsplit(genbank_record,",")[[1]][1],"\n"),file=paste0("GenbankAccessions_",dataset,"_plant.csv"),append=TRUE,quote=FALSE,row.names = FALSE,col.names=FALSE,eol="")
}

# build plant tree with V.Phylomaker2 : build table of plants and build tree
library("V.PhyloMaker2")
plant_list <- read.csv(paste0("GenbankAccessions_",dataset,"_plant.csv"),sep=" ", header=F)
plant_list$V1 <- paste0(plant_list$V2," ",plant_list$V3)
plant_list$V3 <- "Orchidaceae"
plant_list <- plant_list[, 1:3]
colnames(plant_list) <- c("species","genus","family")
plant_tree <- phylo.maker(sp.list = plant_list, tree = GBOTB.extended.TPL, nodes = nodes.info.1.TPL, scenarios = "S3")$scenario.3
write.tree(plant_tree,"plant_tree_Shefferson.tre")

# building of the fasta file for fungi
for (i in 1007:1067){
  genbank_id <- paste0("HM14",i)
  
  genbank_record <- entrez_fetch(db="nucleotide", id=genbank_id, rettype="fasta", retmode="text")
  
  write.table(genbank_record,file=paste0("sequences_",dataset,"_fungi.fasta"),append=TRUE,quote=FALSE,row.names = FALSE,col.names=FALSE,eol="")
}

# plant tree
name_plant_tree="plant_tree_Shefferson.tre"
plant_tree <- ape::read.tree(name_plant_tree)
plant_tree <- ladderize(plant_tree)
plot(plant_tree)

# do not use Trimal for this alignment
name_fungi_tree=paste0("alignment_",dataset,"_fungi.fasta.treefile")
fungi_tree <- ape::read.tree(name_fungi_tree)
fungi_tree <- root(fungi_tree, node = getMRCA(tip=c("HM141067.1","HM141050.1"), fungi_tree))
plot(fungi_tree)

# building of the network
network <- data.frame(matrix(0,length(plant_tree$tip.label),length(fungi_tree$tip.label)),row.names = plant_tree$tip.label)
colnames(network)=fungi_tree$tip.label

# complete the network
library(rentrez)
for (i in 1007:1067){
  genbank_id <- paste0("HM14",i,".1")
  host <- entrez_fetch(db = "nuccore", id = genbank_id, rettype = "gb", retmode = "text")
  host <- strsplit(host,"host=")[[1]][2]
  host <- strsplit(host,"\n")[[1]][1]
  host <- gsub("\"","",host)
  host <- strsplit(host," var")[[1]][1]
  host <- gsub(" ", "_", host)
  network[[host,genbank_id]]<-1
}

###### Script Shefferson_Ceratobasidiaceae ######
dataset="Shefferson_Ceratobasidiaceae"

# same fasta file for plants
# same tree for plants

# building of the fasta file for fungi
for (i in 1007:1011){
  genbank_id <- paste0("HM14",i)
  
  genbank_record <- entrez_fetch(db="nucleotide", id=genbank_id, rettype="fasta", retmode="text")
  
  write.table(genbank_record,file=paste0("sequences_",dataset,"_fungi.fasta"),append=TRUE,quote=FALSE,row.names = FALSE,col.names=FALSE,eol="")
}
for (i in 1013:1046){
  genbank_id <- paste0("HM14",i)
  
  genbank_record <- entrez_fetch(db="nucleotide", id=genbank_id, rettype="fasta", retmode="text")
  
  write.table(genbank_record,file=paste0("sequences_",dataset,"_fungi.fasta"),append=TRUE,quote=FALSE,row.names = FALSE,col.names=FALSE,eol="")
}
genbank_record <- entrez_fetch(db="nucleotide", id="AJ427405", rettype="fasta", retmode="text")
write.table(genbank_record,file=paste0("sequences_",dataset,"_fungi.fasta"),append=TRUE,quote=FALSE,row.names = FALSE,col.names=FALSE,eol="")

# fungi tree
fungi_outgroup="AJ427405.1"

# same for building of the network
# complete the network
library(rentrez)
for (i in 1007:1011){
  genbank_id <- paste0("HM14",i,".1")
  host <- entrez_fetch(db = "nuccore", id = genbank_id, rettype = "gb", retmode = "text")
  host <- strsplit(host,"host=")[[1]][2]
  host <- strsplit(host,"\n")[[1]][1]
  host <- gsub("\"","",host)
  host <- strsplit(host," var")[[1]][1]
  host <- gsub(" ", "_", host)
  network[[host,genbank_id]]<-1
}
for (i in 1013:1046){
  genbank_id <- paste0("HM14",i,".1")
  host <- entrez_fetch(db = "nuccore", id = genbank_id, rettype = "gb", retmode = "text")
  host <- strsplit(host,"host=")[[1]][2]
  host <- strsplit(host,"\n")[[1]][1]
  host <- gsub("\"","",host)
  host <- strsplit(host," var")[[1]][1]
  host <- gsub(" ", "_", host)
  network[[host,genbank_id]]<-1
}

##### Script VanGalen #####
dataset="VanGalen"

# fungi tree
name_fungi_tree="NPH_18802_Supporting Information Notes S3.nwk"
fungi_tree <- read.tree(name_fungi_tree)
for (i in fungi_tree$tip.label){
  name <- strsplit(i,"_")[[1]][3]
  fungi_tree$tip.label[fungi_tree$tip.label==i] <- paste0("asv_",name)
}
# delete fungi which are not in the network

# plant tree
set.seed(1)
plant_tree <- pbtree(n=4,tip.label=c("N.fusca","N.cliffortioides","N.solandri","N.menziesii"))
plot(plant_tree)

# build network with the supplementary data
network <- read.table("TableS4.csv", sep = ";")
network <- cbind(network,rep(0,length(network$V1)),rep(0,length(network$V1)),rep(0,length(network$V1)),rep(0,length(network$V1)))
colnames(network) <- c("fungi","N. fusca","N. cliffortioides","N. solandri","N. menziesii")
network <- network[grep("asv_|N. ",network$fungi),]
row.names(network) <- network$fungi
network <- network[,-1]
species=network[[1,"fungi"]]
for (i in rownames(network)){
  if (startsWith(i,"asv_")){
    network[[i,species]] <- 1
  }
  else{
    species=i
  }
}
network <- network[grep("asv_",rownames(network)),]
colnames(network) <- c("N.fusca","N.cliffortioides","N.solandri","N.menziesii")
network <- t(network)
write.table(network,"Network.csv",sep=";",quote=FALSE)

#
network <- read.table("Network.csv",sep=";")

###### Script Toju #####
dataset="Toju"

# Assess the most abundant fungi orders
fungi_clades <- read.table("fungi_clades.csv",sep=";",header=T)
nrow(fungi_clades[fungi_clades$Order=="",])
for (i in (1:nrow(fungi_clades))) {
  if (fungi_clades[[i,"Order"]]==""){
    fungi_clades[[i,"Order"]]<-fungi_clades[[i,"Order_Unite"]] #using of Unite data when Claident data is not available
  }
}
nrow(fungi_clades[fungi_clades$Order=="",])
clades_abundance <- data.frame(unique(fungi_clades$Order),rep(0,length(unique(fungi_clades$Order))),rep(0,length(unique(fungi_clades$Order)))) # dataframe of percentage of reads per order
colnames(clades_abundance) <- c("Order","Reads_percentage","Nb_fungi")
nb_reads <- sum(fungi_clades$Number.of.reads.in.the.sample.level.matrix)
for (i in (1:nrow(fungi_clades))) {
  clades_abundance[clades_abundance$Order==fungi_clades[[i,"Order"]],][,"Reads_percentage"] <- clades_abundance[clades_abundance$Order==fungi_clades[[i,"Order"]],][,"Reads_percentage"]+fungi_clades[[i,"Number.of.reads.in.the.sample.level.matrix"]]*100/nb_reads
  clades_abundance[clades_abundance$Order==fungi_clades[[i,"Order"]],][,"Nb_fungi"] <- clades_abundance[clades_abundance$Order==fungi_clades[[i,"Order"]],][,"Nb_fungi"]+1
}
clades_abundance <- clades_abundance[order(clades_abundance$Reads_percentage,decreasing=T),]

# fasta files for fungi
fungi=c("Helotiales","Sebacinales","Agaricales")
all_fungi_sequences <-  seqinr::read.fasta(paste0("sequences_",dataset,"_fungi.fasta"))
for (i in fungi){
  order <- fungi_clades[fungi_clades$Order==i,"Fungus"]
  order <- c(order,"F_7") # a Chaetothyriales as extragroup
  order_sequences <- all_fungi_sequences[which(names(all_fungi_sequences) %in% order)]
  write.dna(order_sequences, paste0("sequences_",dataset,"_",i,".fasta"), format = "fasta", nbcol = -1,colsep="",colw=10000)
}

#alignment and building of tree for fungi

# build plant tree with V.Phylomaker2 : build table of plants and build tree
library("V.PhyloMaker2")
genus_list <- rep(0,length(row.names(network)))
i=1
for (species in row.names(network)) {
  name <- strsplit(species,"_")[[1]][1]
  genus_list[i] <- name
  i=i+1
}
plant_list <- data.frame(row.names(network),genus_list,rep("Ericaceae",length(genus_list)))
colnames(plant_list) <- c("species","genus","family")
plant_list[[13,"family"]] <- "Pinaceae"
plant_tree <- phylo.maker(sp.list = plant_list, tree = GBOTB.extended.TPL, nodes = nodes.info.1.TPL, scenarios = "S3")$scenario.3
write.tree(plant_tree,"plant_tree_Toju.tre")

#load network
fungi="Helotiales"
backbone="-backbone"
fungi="Sebacinales"
backbone="-backbone"
fungi="Agaricales"
backbone=""
network <- read.table("Network_Toju.csv",sep=";")
network[network>0] <- 1
#subsample network
order <- fungi_clades[fungi_clades$Order==fungi,"Fungus"]
length(intersect(order,names(network)))
network <- network[,intersect(order,names(network))]

# load trees
name_plant_tree <- "plant_tree_Toju.tre"
name_fungi_tree <- paste0("alignment_",dataset,"_",fungi,backbone,"_trimal.fasta.treefile")
fungi_outgroup <- "F_7"

#for global fit methods
write.table(res,paste0("C:/fantine/Nextcloud/Immex_2023/results/GlobalFit_",dataset,"_",fungi,backbone,".csv"),sep=";",row.names=FALSE,quote=FALSE)

# Just before preparing empress : for each fungi, select the plant where the fungi is the most abundant
network <- read.table("Network_Toju.csv",sep=";")
order <- fungi_clades[fungi_clades$Order==fungi,"Fungus"]
length(intersect(order,names(network)))
network <- network[,intersect(order,names(network))]
for (j in colnames(network)) {
  for (i in rownames(network)) {
    if (network[[i,j]]<max(network[,j])){
      network[[i,j]]=0
    }
  }
}
network[network>0] <- 1


dataset <- paste0("Toju_",fungi)

# Empress
bifurcation <- TRUE
suffix <- if (bifurcation) "_bifurcations" else ""
fungi=c("Helotiales","Sebacinales","Agaricales")
dataset="Toju"
for (order in fungi){
  empress=data.frame()
  
  parameters<-data.frame()
  parameters<-rbind(c(1,1,1),c(4,1,1),c(2,1,2),c(4,2,1),c(2,3,1))
  colnames(parameters)<-c("d","t","l")
  
  for (i in 1:5) {
    d=parameters[i,1]
    t=parameters[i,2]
    l=parameters[i,3]
    costs=paste0("d",d,"_t",t,"_l",l)
    
    reconciliation <- read.table(paste0("C:/fantine/Nextcloud/Immex_2023/data/data_",dataset,"/recon_",dataset,"_",order,"_",costs,suffix,"_output.csv"),sep=",")
    cospeciation <- length(reconciliation$V3[reconciliation$V3=="Cospeciation"])/length(reconciliation$V3[reconciliation$V3=="Transfer"]) #indique si l'on a plus d'évenements de co-speciation que de tranferts
    
    graph <- read.table(paste0("C:/fantine/Nextcloud/Immex_2023/data/data_",dataset,"/pvalue_",dataset,"_",order,"_",costs,suffix,"_output.svg"), comment.char = "", fill=TRUE, sep=";")
    pvalue_empress <- graph$V1[grep("p-value", graph$V1)] #isole la ligne du doc avec la pvalue
    pvalue_empress <- gsub("    <!-- p-value = ", "", pvalue_empress)
    pvalue_empress <- as.numeric(gsub(" -->", "", pvalue_empress)) #isole la pvalue
    
    
    empress <- rbind(empress,c(costs,cospeciation ,pvalue_empress))
    colnames(empress) <- c("Costs","Ratio cospeciation / transfer events","P-value")
  }
  
  write.table(empress,paste0("C:/fantine/Nextcloud/Immex_2023/results/Empress_",dataset,"_",order,suffix,".csv"),sep=";",row.names=FALSE,quote=FALSE)
}

##### Script Sepp #####
dataset="Sepp"

# build plant tree with V.Phylomaker2 : build table of plants and build tree
library("V.PhyloMaker2")
plant_list <- read.table("Plant_clades.csv", sep=";", header=T)
plant_list <- head(plant_list, -3) #if you want to delete ferns and shrubs
i=1
for (species in plant_list$Species) {
  plant_list[[i,"Species"]] <- paste0(plant_list[[i,"Genus"]]," ",species)
  i=i+1
}
plant_tree <- phylo.maker(sp.list = plant_list, tree = GBOTB.extended.TPL, nodes = nodes.info.1.TPL, scenarios = "S3")$scenario.3
write.tree(plant_tree,paste0("plant_tree_",dataset,".tre"))

# list of fungi and their host (/!\ very long)
library(rentrez)
links <- data.frame()
for (i in 1301:4406){
  genbank_id <- paste0("LT98",i)
  information <- entrez_fetch(db = "nuccore", id = genbank_id, rettype = "gb", retmode = "text")
  fungi <- strsplit(information,"DEFINITION ")[[1]][2]
  fungi <- strsplit(fungi,",")[[1]][2]
  fungi <- strsplit(fungi," ")[[1]]
  fungi <- fungi[length(fungi)]
  host <- strsplit(information,"host=")[[1]][2]
  host <- strsplit(host,"\n")[[1]][1]
  host <- gsub("\"","",host)
  host <- strsplit(host," var")[[1]][1]
  links <- rbind(links,c(host,fungi))
}
colnames(links) <- c("Plant","Fungi")
write.table(links,"network.csv", row.names=F,quote=F, sep=";")

# building of network
links <- read.table("network.csv",sep=";",header=T)
length(unique(links$Fungi))
length(unique(links$Fungi[grep("VT",links$Fungi)]))
length(unique(links$Plant))
#gsub(" ","_",unique(links$Plant))[!gsub(" ","_",unique(links$Plant)) %in% plant_tree$tip.label]

network <- data.frame(matrix(0,length(unique(links$Plant)),length(unique(links$Fungi[grep("VT",links$Fungi)]))))
rownames(network) <- unique(links$Plant)
colnames(network)=unique(links$Fungi[grep("VT",links$Fungi)])
for (i in (1:length(links$Fungi))) {
  if (grepl("VT",links[[i,"Fungi"]])){
    network[[links[[i,"Plant"]],links[[i,"Fungi"]]]] <- network[[links[[i,"Plant"]],links[[i,"Fungi"]]]] + 1
  }
}
row.names(network) <- gsub(" ", "_", row.names(network))
network <- network[rownames(network) != "Frangula_alnus",]
network <- network[rownames(network) != "Juniperus_communis",]
network <- network[rownames(network) != "Asplenium_ruta-muraria",]

#for global fit methods
network[network>0] <- 1
network <- network[, colSums(network) > 0]

# for Empress
for (j in colnames(network)) {
  for (i in rownames(network)) {
    if (network[[i,j]]<max(network[,j])){
      network[[i,j]]=0
    }
  }
}
network[network>0] <- 1

# building of tree file
fungi_tree <- read.nexus("C:/fantine/Nextcloud/Immex_2023/data/data_global_scale/tree_VT_18S_barcode_PB_LN_GTR_constrained.tre")
fungi_tree$tip.label <- gsub("X0000","",fungi_tree$tip.label)
fungi_tree$tip.label <- gsub("X000","",fungi_tree$tip.label)
fungi_tree$tip.label <- gsub("X00","",fungi_tree$tip.label)
#tree with only plants in network


name_plant_tree <- paste0("plant_tree_",dataset,".tre")


##### Script Opik ######
dataset="Opik"

#same script than Sepp for plant tree
name_plant_tree <- paste0("plant_tree_",dataset,".tre")

#network
network <- read.table("network.csv", sep=";", header=T,row.names=1)
network <- t(network)
network[is.na(network)] <- 0
network <- as.data.frame(network)
#some fungi have not the same name in the tree
network$"30" <- network$"30" + network$"33" + network$"34" + network$"37"
network$"143" <- network$"143"+network$"141"
network$"417" <- network$"162"
network$"441" <- network$"145"
network$"104" <- network$"198"
network <- network[,!names(network) %in% c("33","34","37","141","162","145","198")]

# building of fungi tree
fungi_tree <- read.nexus("C:/fantine/Nextcloud/Immex_2023/data/data_global_scale/tree_VT_18S_barcode_PB_LN_GTR_constrained.tre")
fungi_tree$tip.label <- gsub("VTX0000","",fungi_tree$tip.label)
fungi_tree$tip.label <- gsub("VTX000","",fungi_tree$tip.label)
fungi_tree$tip.label <- gsub("VTX00","",fungi_tree$tip.label)
#tree with only fungi in network

#for global fit methods
network[network>0] <- 1

# for Empress
for (j in colnames(network)) {
  for (i in rownames(network)) {
    if (network[[i,j]]<max(network[,j])){
      network[[i,j]]=0
    }
  }
}
network[network>0] <- 1

         