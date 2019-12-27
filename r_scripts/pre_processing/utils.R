# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples


# some useful functions


#import this package if using auc() to calculate areas under curves
#suppressPackageStartupMessages(library("MESS"))






rmv.duplicatedSamp <- function(mat, cor.lim){
	dup_col <- colnames(mat)[duplicated(colnames(mat))]

	to_keep <- sapply( dup_col, function(x) cor(na.exclude(mat[ , which(x == colnames(mat)) ]))[1,2] >= cor.lim )

	to_remove <- dup_col[!to_keep]

	to_remove <- as.numeric( sapply( to_remove, function(x) which(x == colnames(mat)) ) )

	return(mat[,-to_remove])
}



avg.duplicatedSamp <- function(mat){
	dup_col <- colnames(mat)[duplicated(colnames(mat))]

	avg_mat <- sapply( dup_col, function(x) rowMeans(mat[ , which(x == colnames(mat)) ], na.rm = TRUE) )
	colnames(avg_mat) <- dup_col

	mat <- mat[,-as.numeric( sapply( dup_col, function(x) which(x == colnames(mat)) ) )]

	return(cbind(mat, avg_mat))
}



gkd_norm <- function(x, min = -1000, max = 1000){
	pdf <- density(x, na.rm=T, from=min, to=max)
	#f <- approxfun(pdf)
	
	v <- rep(NA, length(x))

	for( i in 1:length(x) ){
		if ( !is.na(x[i]) ) { v[i] <- log( integrate.xy( pdf$x, pdf$y, min, x[i] ) / integrate.xy( pdf$x, pdf$y, x[i], max ) ) }
	}
	
	#{ v[i] <- log( integrate( f, -10, x[i] )$value / integrate( f, x[i], 10 )$value ) }
	
	return(v)
}




uniprot2genename <- function(fasta_file){
	fasta_file <- readAAStringSet(filepath = fasta_file, format = "fasta")
	fasta_file <- as.vector(fasta_file)

	uniprot2genename <- list()

	for(rec in names(fasta_file)){
		if( grepl(pattern = "OS=Homo sapiens", x = rec, fixed = T) & grepl(pattern = "GN=", x = rec, fixed = T) ){
			uniprot = strsplit(rec, split="|", fixed=T)[[1]][2]
          	genename = strsplit(strsplit(rec, split=" GN=", fixed=T)[[1]][2], split=" ", fixed=T)[[1]][1]
            accname = strsplit(strsplit(rec, split="|", fixed=T)[[1]][3], split=" ", fixed=T)[[1]][1]
        	uniprot2genename[[uniprot]] <- c(genename, accname)
		}
	}

	return(uniprot2genename)
}





