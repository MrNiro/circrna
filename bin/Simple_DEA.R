#!/usr/bin/Rscript

# Automated differential expression analysis script for nf-core/circrna
# Relies on Stirngtie prepde.py outputs
# Would be fun to adapt to STAR + Kallisto for own use in future.

get_args <- function(){

    argp <- arg_parser(
            description="this script is designed to automate DESeq2 analysis for circRNA nextflow script",
            hide.opts=TRUE)

    argp <- add_argument(
            parser=argp,
            arg="gene_counts",
            short="g",
            help="gene_count_matrix.csv produced by prepDE.py downstream of Stringtie quant with -e flag",
            default="gene_count_matrix.csv")

    argp <- add_argument(
            parser=argp,
            arg="phenotype",
            short="p",
            help="file containing sample metadata information to guide the design",
            default="phenotype.csv")

    argp <- add_argument(
            parser=argp,
            arg="circRNA",
            short="c",
            help="circRNA counts matrix",
            default="circRNA_matrix.txt")

    argp <- add_argument(
            parser=argp,
            arg="species",
            short="s",
            help="species ID",
            default="hsa")

    argp <- add_argument(
            parser=argp,
            arg="map",
            short="m",
            help="ensDB",
            "default"="ensemblDatabase.txt")

    argv <- parse_args(
            parser=argp,
            argv = commandArgs(trailingOnly = TRUE))

    return(argv)

}

giveError <- function(message){
    cat(paste("\n", message, sep=""))
    quit()
}

usage <- function(){giveError("USAGE: DEA.R <gene_counts.csv> <phenotype.txt> <circRNA_matrix.txt> <species id> <ensembl_map>")}


stage_data <- function(gene_counts, phenotype, circRNA, species, map){

    inputdata <- list()

    dbmap <- read.table(map, sep="\t", header=T, quote="", stringsAsFactors=FALSE)
    gene_mat <- read.csv(gene_counts, row.names="gene_id", check.names=F)
    circ <- read.table(circRNA, sep ="\t", header = T, stringsAsFactors=FALSE)

    # Merge circRNA genomic loci to ID
    circ$circ <- with(circ, paste0(Chr, sep=":", Start, sep="-", Stop, sep=":", Strand))
    rownames(circ) <- circ$circ
    circ <- subset(circ, select=-c(Chr, Start, Stop, Strand, circ))

    ## add pseudocount of 1
    gene_mat <- gene_mat + 1
    circ <- circ + 1

    inputdata$pheno <- checkinputdata(phenotype)
    cols <- rownames(inputdata$pheno)

    if(identical(rownames(inputdata$pheno), colnames(gene_mat))){
        circ <- circ[, cols,]
    }else{
        giveError(c("Samples in phenotype file do not match sequencing sample names.\n",
                    "Please check that phenotype samples match gene_count_matrix.csv headers.\n",
                    "*Make sure they are sorted in alphabetical order:",
                    "tail -n +2 phenotype.txt | sort -k1,1\n\n"))
    }

    inputdata$gene <- gene_mat
    inputdata$circ <- circ
    inputdata$design <- makedesign(inputdata$pheno)
    inputdata$species <- species
    inputdata$map <- dbmap

    inputdata$gene <- ens2symbol(inputdata$gene, inputdata)

    return(inputdata)
}


checkinputdata <- function(phenotype){

    # Stage phenotype file
    pheno <- read.csv(phenotype, row.names=1, header = T, stringsAsFactors = T)

    # Check if there are at least 3 replicates (DESeq2 fails if < 3)
    if(min(table(pheno$condition)) >= 3){
        print("Suitable sample size for DE analysis")
    }else{
        giveError("Not enough samples per condition to perform DE analysis!")
    }

    # Rename sex to gender, how progressive!
    if("sex" %in% names(pheno)){
        print("Renaming sex to gender in phenotype file")
        rename <- gsub("sex", "gender", names(pheno))
        names(pheno) <- rename
    }

    # Check gender is only male, female, unknown
    if ("gender" %in% names(pheno)) {
        if (! all(unique(pheno$gender) %in% c("m", "f", "u"))) {
            giveError("SAMPLEINFO ERROR:\nOnly the values m [male], f [female] and u [unknown] are supported in field <gender>.\n")
        }
    }

    ## check if all columns are factors. If numeric, convert to factor.
    factor_cols <- sapply(pheno, is.factor)
    if(all(factor_cols) == TRUE){
        print("All columns in phenotype are factors and suitable for analysis.")
    }else{
        numeric_cols <- sapply(pheno, is.numeric)
        names <- colnames(pheno)[numeric_cols]
        print(paste0("Column(s) ", names, " is numeric. Converting to factor."))
        pheno[numeric_cols] <- as.data.frame(lapply(pheno[numeric_cols], factor))
        final_check <- sapply(pheno, is.factor)
        if(all(final_check) == TRUE){
            print("Finished coverting to factor")
        }else{
            giveError("Error in converting to factors. See checkinputdata function.")
        }
    }

    return(pheno)

}



makedesign <- function(phenotype){

    # Covariates i.e explanatory variables.
    covariates <- names(phenotype)[which(!names(phenotype) %in% c("condition"))]
    design <- formula(paste("~", paste(c(covariates, "condition"), sep="", collapse=" + ")))
    return(design)

}




ens2symbol <- function(mat, inputdata){

    ## designed to work on input gene_count_matrix.csv file
    ## everything else downstream no longer needs to be converted

    ## figure out if working with ENS, or ENS IDs

    mat <- as.data.frame(mat)
    map <- inputdata$map
    species <- inputdata$species

    if(all(grepl(pattern="^ENSG", rownames(mat)))){
        filter = "ensembl_gene_id"
            if(all(grepl(pattern=".", rownames(mat)))){
                filter = "ensembl_gene_id_version"
            }
    }else{
        filter = "external_gene_name"
    }

    ## set up Mart
    mart_call <- as.character(subset(map$command, map$species == species))
    print("ENS2SYMBOL")
    mart <- eval(str2expression(mart_call))

    ## now go about converting ENS2SYMBOL
    if(filter == "ensembl_gene_id"){

        mat$ensembl_gene_id <- rownames(mat)
        info <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),
                        filters = c("ensembl_gene_id"),
                        values = mat$ensembl_gene_id,
                        mart = mart,
                        useCache=FALSE)

        tmp <- merge(mat, info, by="ensembl_gene_id")
        tmp$external_gene_name <- make.names(tmp$external_gene_name, unique = T)
        rownames(tmp) <- tmp$external_gene_name
        tmp <- subset(tmp, select=-c(ensembl_gene_id, external_gene_name))

        mat <- tmp
        print("input mat ensembl gene id detected and converted")
        return(mat)
    }else if(filter == "ensembl_gene_id_version"){

        mat$ensembl_gene_id_version <- rownames(mat)
        info <- getBM(attributes=c("ensembl_gene_id_version","external_gene_name"),
                    filters = c("ensembl_gene_id_version"),
                    values = mat$ensembl_gene_id_version,
                    mart = mart,
                    useCache=FALSE)

        tmp <- merge(mat, info, by="ensembl_gene_id_version")
        tmp$external_gene_name <- make.names(tmp$external_gene_name, unique = T)
        rownames(tmp) <- tmp$external_gene_name
        tmp <- subset(tmp, select=-c(ensembl_gene_id_version, external_gene_name))

        mat <- tmp
        print("input mat ensembl gene id version detected and converted")
        return(mat)
    }else{
        print("NO change made to input mat ")
        return(mat)
    }

}

# Data type provided at end of script to activate RNA-Seq / circRNA analysis.
DESeq2 <- function(inputdata, data_type){

    if(data_type == "RNA-Seq"){
        outdir <- "RNA-Seq/"

        dds <- DESeqDataSetFromMatrix(
        countData=inputdata$gene,
        colData=inputdata$pheno,
        design = inputdata$design)

        levels <- as.character(unique(inputdata$pheno$condition))
        # for(level in levels){
        reference <- "CTRL"
        contrasts <- levels[levels != paste0(reference)]
        dds$condition <- relevel(dds$condition, ref = paste0(reference))
        dds <- DESeq(dds, quiet=TRUE)

        DESeq2_plots(dds, outdir)

        for(var in contrasts){
            contrast <- paste(var, "vs", reference, sep="_")
            DEG <- getDESeqDEAbyContrast(dds, contrast, reference, var, outdir, inputdata)
        }
        # }
    }else if(data_type == "circRNA"){
        outdir <- "circRNA/"

        ## use gene sizeFactors
        tmp <- DESeqDataSetFromMatrix(
        countData=inputdata$gene,
        colData=inputdata$pheno,
        design = inputdata$design)
        tmp <- DESeq(tmp, quiet=TRUE)

        sizefactors <- sizeFactors(tmp)
        rm(tmp)

        dds <- DESeqDataSetFromMatrix(
        countData=inputdata$circ,
        colData=inputdata$pheno,
        design = inputdata$design)

        levels <- as.character(unique(inputdata$pheno$condition))
        reference <- "CTRL"
        contrasts <- levels[levels != paste0(reference)]
        dds$condition <- relevel(dds$condition, ref = paste0(reference))
        dds <- DESeq(dds, quiet=TRUE)
        sizeFactors(dds) <- sizefactors

        for(var in contrasts){
            contrast <- paste(var, "vs", reference, sep="_")
            DEG <- getDESeqDEAbyContrast(dds, contrast, reference, var, outdir)
        }
    }else{
        giveError("Data type not provided correctly, check end of script")
    }
    return(DEG)
}


getDESeqDEAbyContrast <- function(dds, contrast, reference, var, outdir, inputdata) {

    res <- results(dds, filterFun=ihw, alpha=0.05,  contrast=c("condition", var, reference))
    cat('\n\nSummary data from DESeq2 for ', contrast, ':', sep="")
    summary(res)

    dir <- paste(outdir, contrast, sep="")
    dir.create(dir)

    res_df <- as.data.frame(res)
    out_res_df <- tibble::rownames_to_column(res_df, "ID")
    write.table(out_res_df, file.path(dir, paste("DESeq2", contrast, "whole_differential_expression.txt", sep="_")), sep="\t", row.names=F, quote=F)
    write.table(res, file.path(dir, paste("DESeq2", contrast, "RES_differential_expression.txt", sep="_")), sep="\t", row.names=F, quote=F)
}


options(error=function()traceback(2))
suppressPackageStartupMessages(library("argparser"))
#suppressPackageStartupMessages(library("BiocParallel"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("DESeq2"))
# suppressPackageStartupMessages(library("dplyr"))
#suppressPackageStartupMessages(library("edgeR"))
# suppressPackageStartupMessages(library("EnhancedVolcano"))
#suppressPackageStartupMessages(library("EnsDb.Hsapiens.v86"))
#suppressPackageStartupMessages(library("genefilter")) #for rowVars
# suppressPackageStartupMessages(library("ggplot2"))
# suppressPackageStartupMessages(library("ggpubr"))
#suppressPackageStartupMessages(library("ggrepel"))
#suppressPackageStartupMessages(library("ggfortify"))
# suppressPackageStartupMessages(library("gplots"))
# suppressPackageStartupMessages(library("IHW"))
#suppressPackageStartupMessages(library("limma"))
#suppressPackageStartupMessages(library("parallel"))
# suppressPackageStartupMessages(library("PCAtools"))
# suppressPackageStartupMessages(library("pheatmap"))
# suppressPackageStartupMessages(library("RColorBrewer"))
#suppressPackageStartupMessages(library("readr"))
#suppressPackageStartupMessages(library("Rsubread"))
#suppressPackageStartupMessages(library("tximport"))
#suppressPackageStartupMessages(library("VennDiagram"))

arg <- get_args()

inputdata <- stage_data(arg$gene_counts, arg$phenotype, arg$circRNA, arg$species, arg$map)
dir.create("RNA-Seq")
dir.create("circRNA")
y <- DESeq2(inputdata, "circRNA")
