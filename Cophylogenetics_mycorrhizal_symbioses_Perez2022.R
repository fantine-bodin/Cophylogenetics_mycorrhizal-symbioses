###### Script for the network of Perez 2022 ######

###### Step 1: Observation of data ######


rm(list=ls())

dataset="La_Reunion"
setwd(paste0("C:/fantine/Nextcloud/Immex_2023/data/data_",dataset))

library(ape)
library(paco)

#### Load the data

fungi="Sebacinales"
marker="ITS"

fungi="Glomeromycotina"
marker="18S"

fungi="Cantharellales"
marker="ITS"

fungi="Helotiales"
marker="ITS"

site="Dimetile"
site="Plaine des palmistes"
site="Grand brûlé"

fungi_names=c("Cantharellales","Glomeromycotina","Helotiales","Mucoromycotina","Sebacinales")
fungi_markers=c("ITS-backbone","18S","ITS-backbone","18S","ITS-backbone")
fungi_table<-data.frame(fungi_names,fungi_markers)
sites=c("Grand brûlé","Plaine des palmistes","Dimetile")


OTU="OTU97"

arbre_reseau <- function (fungi,marker,site){
    tree <- ape::read.tree(paste0("alignment_",fungi,"_",marker,"_",OTU,"_trimal.fasta.treefile"))
    
    #look at the structure of the tree (ape)
    tree$tip.label #names of tips
    tree$edge #structure of the tree: for each line, number of the node in column 1 and to which node it is linked in column 2
    tree$edge.length
    plot(tree)
    nodelabels()
    
    if (length(grep("outgroup", tree$tip.label)) == 0) {
    } else if (length(grep("outgroup", tree$tip.label)) > 1) {
      tree <- root(tree, node = getMRCA(tree, tip = tree$tip.label[grep("outgroup", tree$tip.label)]))
    } else {
      tree <- root(tree, outgroup = tree$tip.label[grep("outgroup", tree$tip.label)])
    } #root the tree with the outgroup
    tree <- drop.tip(tree, tip=tree$tip.label[grep("outgroup", tree$tip.label)]) #delete outgroup
    tree <- ladderize(tree)
    
    plot(tree)
    
    plant_tree <- read.tree("tree_plants_phylomatic_La_Reunion.tre")
    plot(plant_tree)
    plant_tree$tip.label <- gsub(pattern = "_", replacement = ".", plant_tree$tip.label) #in plant names, replace spaces by plants
    plant_tree <- ladderize(plant_tree)
    plot(plant_tree)
    
    
    network <- read.table(paste0("network_relative_incidences_",fungi,"_",OTU,"_",site,".csv"), sep=";", header=T)
    
    network[network>0] <- 1
    
    
    colnames(network) %in% plant_tree$tip.label
    rownames(network) %in% tree$tip.label
    
    
    
    !plant_tree$tip.label %in% colnames(network) #check if plants are in network : true if they are not
    plant_tree$tip.label[!plant_tree$tip.label %in% colnames(network)] #list of plants which are not in network
    
    plant_tree <- drop.tip(plant_tree, tip=plant_tree$tip.label[!plant_tree$tip.label %in% colnames(network)]) #tree with only plants in network
    plant_tree <- ladderize(plant_tree) # reorganize branches in a unique way
    plot(plant_tree) 
    
    # same for fungi
    !tree$tip.label %in% rownames(network)
    tree$tip.label[!tree$tip.label %in% rownames(network)]
    tree <- drop.tip(tree, tip=tree$tip.label[!tree$tip.label %in% rownames(network)])
    tree <- ladderize(tree)
    plot(tree)
    
    # order network and tree 
    network <- network[tree$tip.label,plant_tree$tip.label]
    network <- t(network)
    return(list(parasite=tree,host=plant_tree,table=network))
    
    
}
    
###### Step 2: Building of a table for GlobalFit ######

globalfit=data.frame()

for (i in 1:5) {
  fungi=fungi_table[i,1]
  marker=fungi_table[i,2]
  
  for (j in 1:3){
    
    site=sites[j]
    
    L<-arbre_reseau(fungi,marker,site)
    tree <- L$parasite
    plant_tree <- L$host
    network <- L$table
    
    
    res_parafit <- parafit(cophenetic(plant_tree),cophenetic(tree),network,nperm=10000, correction="lingoes")$p.global
    #cophenetic(plant_tree) to convert a tree in a distance matrix
    
    D<-prepare_paco_data(cophenetic(plant_tree),cophenetic(tree),network)
    D <- add_pcoord(D)
    res_paco <- c(PACo(D,nperm=1000, symmetric=TRUE)$gof$p, 1- PACo(D,nperm=1000, symmetric=TRUE)$gof$ss)
    # gof$p is pvalue, 1-(gof$ss) is Rsquare
    
    globalfit <- rbind(globalfit,c(site, fungi, res_parafit,res_paco))
    
  }}

colnames(globalfit) <- c("Site","Fungi","Parafit_p-value","PACo_p-value","PACo_Rcarré")

write.table(globalfit,"GlobalFit_La_Reunion-backbone.csv",sep=";",row.names=FALSE,quote=FALSE)


##### Step 3: Preparation of empress #####

# https://sites.google.com/g.hmc.edu/empress/home
# https://academic.oup.com/bioinformatics/article/37/16/2481/5995312

# sample the most abundant host per symbiont species, at random if several host have the same abundance

set.seed(3)


for (i in 1:5) {
  fungi=fungi_table[i,1]
  marker=fungi_table[i,2]
  
  for (j in 1:3){
    
    site=sites[j]
    
    L<-arbre_reseau(fungi,marker,site)
    tree <- L$parasite
    plant_tree <- L$host
    network <- read.table(paste0("network_relative_incidences_",fungi,"_",OTU,"_",site,".csv"), sep=";", header=T)
    
    
    for (c in colnames(network)) {
      for (r in rownames(network)) {
        if (network[[r,c]]<max(network[r,])){
          network[[r,c]]=0
        }
      }
    }
    
    network[network>0] <- 1
    network<-network[colSums(network)!=0]
    
    plant_tree <- drop.tip(plant_tree, tip=plant_tree$tip.label[!plant_tree$tip.label %in% colnames(network)]) #arbre avec seulement plantes de network
    plant_tree <- ladderize(plant_tree) # reordonne branches de manière unique
    plot(plant_tree) 
    
    # ordonner réseau et arbre 
    network <- network[tree$tip.label,plant_tree$tip.label]
    
    
list_links <- reshape2::melt(as.matrix(network))
# Var1 corresponds to symbionts, Var2 to hosts

list_links <- list_links[list_links$value>0,] # remove no interaction

# pick one interaction at random, and set 0 to others
for (symbiont in unique(list_links$Var1)){
  list_links$value[which(list_links$Var1==symbiont)] <- 0
  list_links$value[as.numeric(sample(size=1, as.character(which(list_links$Var1==symbiont))))] <- 1
}

list_links <- list_links[list_links$value>0,] # remove no interaction

#creation of table of links
write.table(paste0(list_links$Var1,":", list_links$Var2), paste0("links_empress_", dataset,"_", fungi, "_", site, ".txt"), col.names=F, row.names=F, quote=F)

# remove the tips with no interactions and creation of new trees
tree_hosts_interact <- drop.tip(plant_tree, tip=plant_tree$tip.label[which(!plant_tree$tip.label %in% list_links$Var2)])
write.tree(tree_hosts_interact, paste0("tree_hosts_empress_", dataset,"_", fungi, "_", site, ".tre"))
tree_parasites_interact <- drop.tip(tree, tip=tree$tip.label[which(!tree$tip.label %in% list_links$Var1)])

tree_parasites_interact$node.label <- NULL
write.tree(tree_parasites_interact, paste0("tree_parasites_empress_", dataset,"_", fungi, "_", site, ".tre"))

}}

##### Step 4: empress #####
#cf code bash empressLaReunion

#### Step 5: Treating of empress data #####

empress=data.frame()

for (i in 1:5) {
  fungi=fungi_table[i,1]
  marker=fungi_table[i,2]
  
  for (j in 1:3){
    
    site=sites[j]
    
    reconciliation <- read.table(paste0("C:/fantine/Nextcloud/Immex_2023/data/data_",dataset,"/recon_",dataset,"_",fungi,"_",site,"_output.csv"),sep=",")
    cospeciation <- length(reconciliation$V3[reconciliation$V3=="Cospeciation"])/length(reconciliation$V3[reconciliation$V3=="Transfer"]) #indique si l'on a plus d'évenements de co-speciation que de tranferts

    graph <- read.table(paste0("C:/fantine/Nextcloud/Immex_2023/data/data_",dataset,"/pvalue_",dataset,"_",fungi,"_",site,"_output.svg"), comment.char = "", fill=TRUE, sep=";")
    pvalue_empress <- graph$V1[grep("p-value", graph$V1)] #isole la ligne du doc avec la pvalue
    pvalue_empress <- gsub("    <!-- p-value = ", "", pvalue_empress)
    pvalue_empress <- as.numeric(gsub(" -->", "", pvalue_empress)) #isole la pvalue
    
    empress <- rbind(empress,c(site, fungi, cospeciation ,pvalue_empress))
  }}

colnames(empress) <- c("Site","Fungi","Ratio cospeciation / transfer events","P-value")

write.table(empress,paste0("C:/fantine/Nextcloud/Immex_2023/results/Empress_",dataset,"backbone.csv"),sep=";",row.names=FALSE,quote=FALSE)


##### Size of the network ######
size=data.frame()

for (i in 1:5) {
  fungi=fungi_table[i,1]
  marker=fungi_table[i,2]
  
  for (j in 1:3){
    
    site=sites[j]
    
    L<-arbre_reseau(fungi,marker,site)
    network <- L$table
    
    size <- rbind(size,c(fungi,site,length(rownames(network)),length(colnames(network))))
    colnames(size) <- c("Fungi","Site","Nb of plants","Nb of fungi")
  }}

write.table(size,"Networks_sizes_La_Reunion.csv",sep=";",row.names=FALSE,quote=FALSE)
