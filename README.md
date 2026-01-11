# Gene-network meta-analysis of prioritized genes from two cohorts of ME/CFS patients

## Methods

### Data sources

The following two independently published ME/CFS-associated gene modules were used as primary inputs:

| Number of cases | Sequencing Method | Gene-Mapping Method    | Genes | Criteria | Reference |
|----------------:|:------------------|:-----------------------|---------|:---------|:----------|
|464              |WGS                | Deep Learning          | 115 | ICC-IOM  |([Zhang S 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12047926/))|
|14767            |Axiom UKB array    | Combinatorial analysis | 259 | CCC-IOM  |([Sardell JM 2025](https://www.medrxiv.org/content/10.64898/2025.12.01.25341362v1))|

### Protein-protein interactions and gene identifiers

Protein–protein interaction (PPI) data were obtained from STRING v12.0 (Homo sapiens). If not already present locally, these datasets were automatically downloaded and cached. Interactions were filtered using a combined score threshold of 0.4, corresponding to medium confidence interactions. All PPI analyses were performed locally.
NCBI Entrez Gene IDs were assigned using a local NCBI gene_info database, which is downloaded automatically. 

### Gene-network analysis

PPIs from STRING v12.0 are stored in a squared, symmetric matrix. Component analysis, degree assignment, and plotting is performed using the R package igraph. STRING is a PPI database that integrates evidence from gene expression, proteomics, curated databases, and text mining ([Szklarczyk D 2023](https://pubmed.ncbi.nlm.nih.gov/36370105/)).

### Over-representation analysis

Pathway, tissue, disease, and cellular component enrichment analysis were performed by over-representation analysis (ORA) employing the following R tools:

| R package      | R function | Database | Reference |
|:---------------|:-----------|:---------|:----------|
|ClusterProfile  |enrichKEGG()|KEGG      |([Yu G 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3339379/)) |
|ReactomePA      |enrichPA()  |Reactome  |([Yu G 2016](https://europepmc.org/article/med/26661513)) |
|DOSE            |enrichDO()  |DiseaseOntology | ([Guangchuang Y 2015](https://pubmed.ncbi.nlm.nih.gov/25677125/)) |
|TissueEnrich    |teEnrichment() | Human Protein Atlas | ([Ashish J 2019](https://pubmed.ncbi.nlm.nih.gov/30346488/)) |
|ClusterProfile  |enrichGO()  |GeneOntoogy | ([Yu G 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3339379/)) |

The background universe consisted of all STRING genes successfully mapped to Entrez IDs (size = 19,338). Multiple testing correction was applied using the Benjamini–Hochberg procedure, and only terms with adjusted p ≤ 0.05 were retained. Visualization of genes mapped to terms, was performed by function GOchord of package GOPlot.  

## Results

### Gene overlap and gene-network representation

The union of the two gene sets generates a list of 369 genes (Figure 1). The overlap between the two gene sets is of four genes, namely MAX, NEDD9, DLGAP2, and GABBR. While the overlap is not statistically significant (p-value of 0.064, by hypergeometric test on a universe of 19,338 genes), the two gene sets converge to a connected graph of 276 genes that includes 112 of the 115 genes prioritized by Zhang and colleagues and 168 of the 259 genes (64%) prioritized by Precision Life. 

![All_genes_graph](https://github.com/user-attachments/assets/2b5d4f66-04e7-4a2d-9157-ccf6edbbe124)

<p align="left">
  <em>Figure 1. Protein-protein interactions for the genes prioritized by Precision Life (PL) by combinatorial analysis using DecodeME GWAS data (green) and by Zhang and colleges by deep learning on WGS from 464 ME/CFS cases (red). Overlapping genes are in white. The edges represents STRING interactions with a probability above 0.4 and their length is proportional to the probability of interaction. This lists includes a total of 369 elements. </em>
</p>

The list of 369 genes, with NCBI Entrez ID and STRING preferred name, is included in this repository (see `all_genes.tsv`). Also, the columns indicate the component they belong to and a score derived from their degree (average PPI score) normalized between zero and one. The same list is included in a format that can be used on Cytoscape (`All_genes_cytoscape.tsv`).

### Tissue enrichment analysis

The results of over-representation analysis of the 369 genes of Figure 1 over the HUman Protein Atlas database, gives Figure 2. 

![Tissue_ORA](https://github.com/user-attachments/assets/fb20ec34-3cf1-4d1e-8100-f1df5bb116fc)

<p align="left">
  <em>Figure 2. Tissue-enrichment analysis by hypergeometric test of the 369 genes of Figure 1 against the Human Protein Atlas (HPA) database, with a background universe of 19699 genes. </em>
</p>

### Cellular Component analysis

We used Gene Ontology Cellular Component (GO CC) to further specify the result of tissue-enrichment. Below and in Figure 3, we report the top 10 results from GO CC enrichment. For the complete results, including also KEGG, Reactome, and Dsease Ontology, see file `ORA_Disease_Module.tsv`.

| ID          | Description                                     | GeneRatio | pvalue           | p.adjust         
|------------|-------------------------------------------------|-----------|-----------------|-----------------|
| GO:0098984 | neuron to neuron synapse                        | 26/350    | 1.23e-07        | 5.29e-05        | 
| GO:0099572 | postsynaptic specialization                     | 24/350    | 7.86e-07        | 1.56e-04        | 
| GO:0019774 | proteasome core complex, beta-subunit complex  | 5/350     | 1.08e-06        | 1.56e-04        | 
| GO:0032279 | asymmetric synapse                              | 22/350    | 3.86e-06        | 4.15e-04        | 
| GO:0014069 | postsynaptic density                            | 21/350    | 6.70e-06        | 5.76e-04        | 
| GO:0000502 | proteasome complex                              | 8/350     | 3.03e-05        | 0.00194         | 
| GO:0005839 | proteasome core complex                          | 5/350     | 3.16e-05        | 0.00194         |
| GO:0043198 | dendritic shaft                                 | 6/350     | 1.25e-04        | 0.00608         | 
| GO:0031594 | neuromuscular junction                           | 8/350     | 1.27e-04        | 0.00608         | 
| GO:0071339 | MLL1 complex                                    | 5/350     | 2.91e-04        | 0.0125          | 

![GO enrichment](https://github.com/user-attachments/assets/a0162a7a-c52f-4dd6-9738-83e21e5addf3)

<p align="left">
  <em>Figure 3. Gene Ontology Cellular Component over-representation analysis of the 369 genes of Figure 1, with a background universe of 19699 genes. Top ten results, after Benjamini-Hochberg correction. See `ORA_Disease_Module.tsv` for the full list significant results. </em>
</p>

### Disease Ontology and KEGG

Over-representation analysis (ORA) with the Disease Ontology database identifies the following two phenotypes as significantly overlapping with the disease module.

| ID          | Description                                     | GeneRatio | pvalue           | p.adjust         |
|------------|-------------------------------------------------|-----------|-----------------|-----------------|
| DOID:0110734 | neurodegeneration with brain iron accumulation | 5/206     | 7.64e-05        | 0.0274          | 
| DOID:1561   | cognitive disorder                              | 28/206    | 8.28e-05        | 0.0274          | 

ORA on KEGG identifies as Spincerebellar ataxia as top enrichment (see Figure 4).

<img width="1467" height="1736" alt="hsa05017 Disease_ORA" src="https://github.com/user-attachments/assets/34455da4-8886-4197-b516-a5c1dbaab432" />

<p align="left">
  <em>Figure 4. Top over-representation analysis on KEGG. Nodes in red indicate the overlap with the disease module of 369 genes. See `ORA_Disease_Module.tsv` for full list of significant results. </em>
</p>

### Comparison with proteomic results

