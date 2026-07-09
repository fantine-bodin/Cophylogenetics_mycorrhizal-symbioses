## Questioning the Evidence for Host–Symbiont Codiversification in Mycorrhizal Symbioses



This repository contains the code and data used in the study "Questioning the Evidence for Host-Symbiont Codiversification in Mycorrhizal Symbioses". It provides all the scripts required to reproduce the cophylogenetic analyses, phylogenetic reconstructions, and associated datasets used throughout the study.





# Available scripts: 



The repository includes the R scripts used to perform cophylogenetic analyses on the mycorrhizal associations described in 13 published studies. 



For each mycorrhizal network, the analysis pipeline consists of three main steps. First, the data are formatted to generate the host and fungal phylogenies together with the corresponding association matrices. Second, cophylogenetic signal is assessed using the global-fit methods ParaFit and PACo. Third, phylogenetic congruence is evaluated using the event-based method eMPRess, which requires formatting the input files, running the analyses through Bash scripts, and importing and processing the reconciliation results. The complete workflow was applied independently to every mycorrhizal network included in the study.



The repository also contains the Bash scripts used to run eMPRess. This software reconstructs the evolutionary history of host and symbiont lineages using a maximum-parsimony framework that models four types of evolutionary events: cospeciation, duplication, host transfer, and loss. Throughout this study, eMPRess was run using an event-cost scheme assigning a cost of 0 to cospeciation, 4 to duplication, and 1 to both host transfer and loss events.



For mycorrhizal networks for which phylogenetic trees were not already available, representative DNA sequences were retrieved for both plant hosts and fungal symbionts. The repository includes the scripts used to align these sequences with MAFFT, trim the alignments with trimAl, and reconstruct phylogenetic trees with IQ-TREE.





# Available data: 



The repository also contains the data associated with each mycorrhizal network, including the plant and fungal phylogenetic trees and the corresponding host-symbiont association matrices. Most reconstructed phylogenies are provided as files named `alignment\_<Name>\_plant\_trimal.fasta.treefile` or `alignment\_<Name>\_fungi\_trimal.fasta.treefile`. These files are in \*\*Newick format and include branch support values\*\*. 

The sequence alignments used for phylogenetic reconstruction are also provided. Plant phylogenies named `plant\_tree\_<Name>` were generated using V.PhyloMaker2, whereas some fungal phylogenies, named `tree\_backbone\_<Name>`, were refined using a backbone-based approach. The alignments used for these backbone reconstructions are provided as `alignment\_all\_fungi\_trimal.fasta`. Finally, the host-symbiont association data are available either as CSV files or as eMPRess-formatted files named `links\_empress\_<Name>.txt`.





