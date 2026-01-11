# Gene-network meta-analysis of prioritized genes from two cohorts of ME/CFS patients

## Methods

### Data sources

The following two independently published ME/CFS-associated gene modules were used as primary inputs:

| Number of cases | Sequencing Method | Gene-Mapping Method    | Criteria | Reference |
|-----------------|:------------------|:-----------------------|:---------|:----------|
|464              |WGS                | Deep Learning          | ICC-IOM  |([Zhang S 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12047926/))|
|14767            |Axiom UKB array    | Combinatorial analysis | CCC-IOM  |([Sardell JM 2025](https://www.medrxiv.org/content/10.64898/2025.12.01.25341362v1))|

### Protein-protein interactions and gene identifiers

Protein–protein interaction (PPI) data were obtained from STRING v12.0 (Homo sapiens). If not already present locally, these datasets were automatically downloaded and cached. Interactions were filtered using a combined score threshold of 0.4, corresponding to medium confidence interactions. All PPI analyses were performed locally.
NCBI Entrez Gene IDs were assigned using a local NCBI gene_info database, which is downloaded automatically. 

### Gene-network analysis

PPIs from STRING v12.0 are stored in a squared, symmetric matrix. Component analysis, degree assignment, and plotting is performed using the R package igraph. STRING is a PPI database that integrates evidence from gene expression, proteomics, curated
databases, and text mining ([Szklarczyk D 2023](https://pubmed.ncbi.nlm.nih.gov/36370105/)).

### Over-representation analysis

Pathway, tissue, disease, and cellular component enrichment analysis were performed by over-representation analysis (ORA) employing the following R tools:

| R package      | R function | Database | Reference |
|:---------------|:-----------|:---------|:----------|
|ClusterProfile  |enrichKEGG()|KEGG      |([Yu G 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3339379/)) |
|ReactomePA      |enrichPA()  |Reactome  |([Yu G 2016](https://europepmc.org/article/med/26661513)) |
|DOSE            |enrichDO()  |DiseaseOntology | ([Guangchuang Y 2015](https://pubmed.ncbi.nlm.nih.gov/25677125/)) |
|TissueEnrich    |teEnrichment() | Human Protein Atlas | ([Ashish J 2019](https://pubmed.ncbi.nlm.nih.gov/30346488/)) |
|ClusterProfile  |enrichGO()  |GeneOntoogy | ([Yu G 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3339379/)) |

The background universe consisted of all STRING genes successfully mapped to Entrez IDs (size = 19,338). Multiple testing correction was applied using the Benjamini–Hochberg procedure, and only terms with adjusted p ≤ 0.05 were retained. 
Visualization of genes mapped to terms, was performed by function GOchord of package GOPlot.  
