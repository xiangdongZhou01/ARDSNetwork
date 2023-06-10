# ARDSNetwork
“Network Analysis Identifies A Gene Biomarker Panel for Sepsis-Induced Acute Respiratory Distress Syndrome”

R scripts

Step.0 Data pre-processing
output: human.filter.expr.RData

step.1 Differential expression analysis
output: human.DEGs.RData

step.2 Weighted gene co-expression network analysis (WGCNA)
output: CytoscapeInput-edges-darkgrey.txt, CytoscapeInput-nodes-darkgrey.txt

Notes: The input data files (including expression file and metadata) are too big to upload. Everyone can download the matrix from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32707.
