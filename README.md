# Mols2020
Code for the project Mols2020 - published in the following paper:

Thomassen, E. E., Sigsgaard, E. E., Jensen, M. R., Olsen, K., Hansen, M. D. D., Svenning, J.-C., & Thomsen, P. F. (2023). Contrasting seasonal patterns in diet and dung-associated invertebrates of feral cattle and horses in a rewilding area. Molecular Ecology, 00, 1– 21. doi:10.1111/mec.16847

# Abstract
Trophic rewilding is increasingly applied in restoration efforts, with the aim of reintroducing the ecological functions provided by large-bodied mammals and thereby promote self-regulating, biodiverse ecosystems. However, empirical evidence for the effects of megafauna introductions on the abundance and richness of other organisms such as plants and invertebrates, and the mechanisms involved still need strengthening. In this study, we use environmental DNA (eDNA) metabarcoding of dung from co-existing feral cattle and horses to assess the seasonal variation in plant diet and dung-associated arthropods and nematodes. We found consistently high diet richness of horses, with low seasonal variability, while the generally lower dietary diversity of cattle increased substantially during summer. Intriguingly, season-specific diets differed, with a greater proportion of trees in the horses' diet during winter, where cattle relied more on shrubs. Graminoids were predominantly found in the diet of horses, but were generally underrepresented compared to previous studies, possibly due to the high prevalence of forbs in the study area. Dung- associated arthropod richness was higher for cattle, largely due to a high richness of flies during summer. Several species of dung- associated arthropods were found primarily in dung from one of the two herbivores, and our data confirmed known patterns of seasonal activity. Nematode richness was constantly higher for horses, and nematode communities were markedly different between the two species. Our results demonstrate complementary effects of cattle and horses through diet differences and dung-associated invertebrate communities, enhancing our understanding of large herbivore effects on vegetation and associated biodiversity. These results are directly applicable for decision- making in rewilding projects, suggesting biodiversity-benefits by inclusion of functionally different herbivores.

## Content
The main directory consist of three subdirectories - one for each primer set used for amplifying DNA, and one for creating maps and plotting sample locations for figure 1 in the associated paper.

ITS: Directory containing all scripts and data files for the amplification of the ITS gene targeting plants (diet)
COI: Directory containing all scripts and data files for the amplification of the COI gene targeting invertebrates
Map: Directory containing a script for plotting sample locations and creating maps - and the generated map files


Each directory includes 2 directories:

"MetaBarFlow", which contain the scripts used to perform the initial bioinformatic treatment of the sequences obtained from the NovaSeq platform. The raw data files are available on DRYAD: https://doi.org/10.5061/dryad.9p8cz8wmb

The scripts in the MetaBarFlow directory origins from the MetaBarFlow pipeline: 

Sigsgaard, E. E., Soraggi, S., Jensen, M. R., Repollés, A. G., Thomassen, E. E., & Thomsen, P. F. (2022). MetaBarFlow (Version 0.1.0) [Computer software]. https://doi.org/10.5281/zenodo.6006700

"Local_analysis", which contains R-scripts for initial filtering and statistical analyses of the data
