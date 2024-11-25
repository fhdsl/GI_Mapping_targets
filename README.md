# gimap - Genetic Interaction MAPping for dual target CRISPR screens

## Background on paired guide CRISPR
Some genes have "backup copies" - these are called paralog genes. Think of it like having a spare tire in your car - if one fails, the other can often do a similar job. This redundancy makes it tricky to study what these genes do, because if you knock out just one gene, its partner might pick up the slack.

CRISPR allows genes to be knocked out so we can see what their function is. But because of paralogs it can be hard to parse out what genes are actually involved. So instead of just targeting one gene at a time, paired guide CRISPR screening, allows us to knockout two genes at a time.
It's particularly useful for understanding:
- What happens when you disable both backup copies of a gene
- How different genes might work together
- Which gene combinations are essential for cell survival

## What is gimap?
gimap - is software tool that helps make sense of paired CRISPR screening data. Here's what it does:
1. Takes data from paired CRISPR screens that has been pre-processed by the pgmap software. 
2. The input data will have cell counts for how well cells grow (or don't grow) when different genes or pairs of genes are disabled
3. gimap can take this data and helps identify interesting patterns, like:
   - When disabling two genes together is more devastating than you'd expect from disabling them individually (called synthetic lethality)
   - When genes work together cooperatively
   - Which gene combinations might be important in diseases like various canceres

gimap can help find meaningful patterns in complex genetic experiments. It's particularly focused on analyzing data from the Berger Lab paired CRISPR screening library called pgPEN (paired guide RNAs for genetic interaction mapping).

The gimap package is based off of the original code and research from the Berger Lab stored in this repository: https://github.com/FredHutch/GI_mapping

## Prerequisites

In order to run this pipeline you will need R and to install the `gimap` package and its dependencies. In R you can run this to install the package:
```
install.packages("remotes")
remotes::install_github("FredHutch/gimap")
```

## Getting Started Tutorial

Now you can [go to our quick start tutorial to get started!](https://fredhutch.github.io/gimap/articles/quick-start.html)

We also have tutorial examples that show how to run timepoint or treatment experimental set ups with gimap:

- [Timepoint example](https://fredhutch.github.io/gimap/articles/timepoint-example.html)
- [Treatment example](https://fredhutch.github.io/gimap/articles/treatment-example.html)

Follow the steps there that will walk you through the example data. Then you can tailor that tutorial to use your own data.

## Citations: 
- https://pubmed.ncbi.nlm.nih.gov/34469736/
- https://www.nature.com/articles/s43586-021-00093-4

![See metrics about this repository here](https://cauldron.io/project/8779/stats.svg)
