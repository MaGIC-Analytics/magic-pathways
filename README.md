# MaGIC Gene Set Enrichment Tool

![GitHub last commit](https://img.shields.io/github/last-commit/MaGIC-Analytics/magic-pathways)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![made with Shiny](https://img.shields.io/badge/R-Shiny-blue)](https://shiny.rstudio.com/)

Run gene set enrichment and over-representation analysis against MSigDB, Reactome, and other curated collections using clusterProfiler. Supports human and mouse annotations and produces standard enrichment plots (dotplot, ridgeplot, cnetplot).

## Running the App
This Shiny App has been built in to a docker container for easy deployment. You can build the image yourself (and thereby customize any ports you need) after downloading it:
```
docker build -t enrichment .
docker run -d --rm -p 8080:8080 enrichment
#Or for testing docker run -t -i --rm -p 8080:8080 enrichment
```
And it should be hosted at localhost:8080
