# Functional Classifier

In this project was trying to find whether some samples were enriched in specific functions, i.e., a method that allows for the classification of organisms according to their functional profile (KO based).

To do this I applied a machine learning approach:

1. Annotate reference proteomes with Mantis (N=5885)
2. Classify them based on enriched function (above q95) - now we have a baseline functional profile for each function enrichment 
3. Train model based on KEGG functional profiles and cluster them into groups (based on enriched functions)
4. Classify isolates based on their functional profile



This project is not cleaned at all. I provided some results to show an example but everything is at a prototype stage.
This was discontinued after my 2nd CET (PhD progress meeting) as I had to focus on other projects. 