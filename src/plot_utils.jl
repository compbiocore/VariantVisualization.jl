"""
    genotype_heatmap2(input::Array{Any,2},title::AbstractString,filename,sample_names,gt_chromosome_labels,y_axis_label_option,x_axis_label_option,save_ext,chrom_label_info)
generate heatmap of genotype data.
"""
function genotype_heatmap2(input,title,filename,sample_names,gt_chromosome_labels,y_axis_label_option,x_axis_label_option,save_ext,chrom_label_info)

    chromosome_label_array=make_chromosome_labels(chrom_label_info)

    chrom=gt_chromosome_labels[:,1]
    pos=gt_chromosome_labels[:,2]

    @rput input
    @rput title
    @rput filename
    @rput sample_names
    @rput chrom
    @rput pos
    @rput save_ext
    @rput chromosome_label_array
    @rput x_axis_label_option

    if y_axis_label_option == "positions"

            reval("

            suppressPackageStartupMessages(library(heatmaply))

            genotypes <- c('no call', 'homozygous reference',
                'heterozygous variant', 'homozygous variant')

            genotype_index <- as.matrix(input) + 1
            genotype_text <- matrix(
                paste('Genotype:', genotypes[genotype_index]),
                ncol = ncol(input)
            )

            d=paste(chrom,pos,sep=',')

            colnames(input)<-sample_names
            row.names(input)<-d

            h=heatmaply(
            #plot_method='plotly',
            input,
            dend=FALSE,
            showticklabels = if (x_axis_label_option=='true') {c(TRUE)} else {c(FALSE,TRUE)},
            limits=c(0,3),
            custom_hovertext = genotype_text,
            Rowv=NULL,
            Colv=NULL,
            label_names = c('Position', 'Sample ID', 'Genotype Value'),
            main = title,
            xlab = 'Sample IDs',
            ylab = 'Chromosomal Positions',
            key.title = 'Genotype'
            )

            h[['x']][['data']][[length(h[['x']][['data']])]][['marker']][['colorbar']][['ticktext']] <- genotypes

            f<-filename
            if (save_ext == 'html') {htmlwidgets::saveWidget(h,file.path(normalizePath(dirname(f)),basename(f)))} else {withr::with_dir(dirname(f), orca(h, basename(f)))

            }

            ")

        elseif y_axis_label_option == "hover_positions"

            reval("

            suppressPackageStartupMessages(library(heatmaply));

            d=paste(chrom,pos,sep=',')

            genotypes <- c('no call', 'homozygous reference',
                'heterozygous variant', 'homozygous variant')

            genotype_index <- as.matrix(input) + 1
            genotype_text <- matrix(
                paste('Genotype:', genotypes[genotype_index]),
                ncol = ncol(input)
            )

            colnames(input)<-sample_names
            row.names(input)<-d

            h=heatmaply(
            input,
            dend=FALSE,
            showticklabels= if (x_axis_label_option=='true') {c(TRUE,FALSE)} else {FALSE},
            limits=c(0,3),
            custom_hovertext = genotype_text,
            Rowv=NULL,
            Colv=NULL,
            label_names = c('Position', 'Sample ID', 'Genotype Value'),
            main = title,
            xlab = 'Sample IDs',
            ylab = 'Chromosomal Positions',
            key.title = 'Genotype'

            )

            h[['x']][['data']][[length(h[['x']][['data']])]][['marker']][['colorbar']][['ticktext']] <- genotypes

            f<-filename
            if (save_ext == 'html') {htmlwidgets::saveWidget(h,file.path(normalizePath(dirname(f)),basename(f)))} else {withr::with_dir(dirname(f), orca(h, basename(f)))
            }

            ")

        elseif y_axis_label_option == "chromosomes"

            reval("

            #Convert array for chromosome row_side_colors colorbar from list to matrix and name with colname
            chromosome_label_array2 <- matrix(unlist(chromosome_label_array))
            colnames(chromosome_label_array2)<-'Chromosome'

            genotypes <- c('no call', 'homozygous reference',
                'heterozygous variant', 'homozygous variant')

            genotype_index <- as.matrix(input) + 1
            genotype_text <- matrix(
                paste('Genotype:', genotypes[genotype_index]),
                ncol = ncol(input)
            )

            #prepare chromosome positions hover labels
            d=paste(chrom,pos,sep=',')

            #m <- data.matrix(input)
            m<-input

            colnames(m)<-sample_names
            row.names(m)<-d

            suppressPackageStartupMessages(library(heatmaply));

            h=heatmaply(
            m,
            dend=FALSE,
            showticklabels= if (x_axis_label_option=='true') {c(TRUE,FALSE)} else {FALSE},
            #plot_method='plotly',
            colorbar_xanchor = 'left',
            colorbar_yanchor = 'top',
            colorbar_xpos = 1,
            colorbar_ypos = 1,
            row_side_colors=chromosome_label_array2,
            row_side_palette = c(
            '1'= '#d3d7cf',
            '2'='#babdb6',
            '3'='#fce94f',
            '4'='#edd400',
            '5'='#c4a000',
            '6'='#8ae234',
            '7'='#4e9a06',
            '8'='#fcaf3e',
            '9'='#f57900',
            '10'='#ce5c00',
            '11'='#e9b96e',
            '12'='#c17d11',
            '13'='#8f5902',
            '14'='#729fcf',
            '15'='#3465a4',
            '16'='#204a87',
            '17'='#ad7fa8',
            '18'='#75507b',
            '19'='#5c3566',
            '20'='#888a85',
            '21'='#555753',
            '22'='#2e3436',
            'X'='#ef2929',
            'Y'='#cc0000',
            'M'='#a40000'),
            limits=c(0,3),
            custom_hovertext = genotype_text,
            Rowv=NULL,
            Colv=NULL,
            label_names = c('Position', 'Sample ID', 'Genotype Value'),
            main = title,
            xlab = 'Sample IDs',
            ylab = 'Chromosomal Positions',
            key.title = 'Genotype'
            )

            h[['x']][['data']][[length(h[['x']][['data']])]][['marker']][['colorbar']][['ticktext']] <- genotypes

            f<-filename
            if (save_ext == 'html') {htmlwidgets::saveWidget(h,file.path(normalizePath(dirname(f)),basename(f)))} else {withr::with_dir(dirname(f), orca(h, basename(f)))
            }
            ")

        else
            println("--y_axis_labels is not a valid option. Choose positions or hover_positions")
        end
end

"""
   genotype_heatmap_with_groups(input::Array{Int64,2},title::String,filename,sample_names,gt_chromosome_labels,pheno_matrix,trait_labels,y_axis_label_option,x_axis_label_option,save_ext,chrom_label_info)
generate heatmap of genotype data.
"""
function genotype_heatmap_with_groups(input,title,filename,sample_names,gt_chromosome_labels,pheno_matrix,trait_labels,y_axis_label_option,x_axis_label_option,save_ext,chrom_label_info)

    chromosome_label_array=make_chromosome_labels(chrom_label_info)

    pheno_rows=size(pheno_matrix,1)-1
    pheno_cols=size(pheno_matrix,2)

    chrom=gt_chromosome_labels[:,1]
    pos=gt_chromosome_labels[:,2]

    @rput input
    @rput title
    @rput filename
    @rput sample_names
    @rput chrom
    @rput pos
    @rput pheno_matrix
    @rput trait_labels
    @rput save_ext
    @rput chromosome_label_array
    @rput x_axis_label_option
    @rput pheno_rows
    @rput pheno_cols

    if y_axis_label_option == "positions"

        reval("

        #set object genotypes to hold titles for legend. Set genotype_index and genotype_index to build labeled legend
        genotypes <- c('no call', 'homozygous reference',
            'heterozygous variant', 'homozygous variant')

        genotype_index <- as.matrix(input) + 1
        genotype_text <- matrix(
            paste('Genotype:', genotypes[genotype_index]),
            ncol = ncol(input)
        )

        #prepare chromosome positions hover labels
        d=paste(chrom,pos,sep=',')

        #transpose phenotype matrix to correct automatic transposition in heatmaply function.
        pheno_matrix=t(pheno_matrix)
        #Remove sample names from pheno matrix
        pheno_matrix=pheno_matrix[,-1]
        #transpose trait labels to match dims of pheno_matrix
        trait_labels=t(trait_labels)
        #convert pheno_matrix from type list to matrix
        pheno_matrix2 <- matrix(unlist(pheno_matrix, use.names=FALSE), ncol = pheno_rows, byrow = FALSE)
        #assign sample names as rownames and trait_labels as col_names on pheno_matrix because pheno_matrix is automatically transposed in heatmaply function
        row.names(pheno_matrix2)<-sample_names
        colnames(pheno_matrix2)<-trait_labels

        #assign sample ids to input matrix column names and chromosome positions to row names for automatically generated tick labels
        colnames(input)<-sample_names
        row.names(input)<-d

        suppressPackageStartupMessages(library(heatmaply))

        h=heatmaply(
        plot_method='plotly',
        input,
        dend=FALSE,
        showticklabels= if (x_axis_label_option=='true') {c(TRUE)} else {c(FALSE,TRUE)},
        limits=c(0,3),
        custom_hovertext = genotype_text,
        Rowv=NULL,
        Colv=NULL,
        col_side_colors=pheno_matrix2,
        col_side_palette = c('1'= 'mediumpurple2',
        '2'='steelblue2'),
        label_names = c('Position', 'Sample ID', 'Genotype Value'),
        main = title,
        xlab = 'Sample IDs',
        ylab = 'Chromosomal Positions',
        key.title = 'Genotype'
        )

        h[['x']][['data']][[length(h[['x']][['data']])]][['marker']][['colorbar']][['ticktext']] <- genotypes

        f<-filename
        if (save_ext == 'html') {htmlwidgets::saveWidget(h,file.path(normalizePath(dirname(f)),basename(f)))} else {withr::with_dir(dirname(f), orca(h, basename(f)))

        }
        ")

    elseif y_axis_label_option == "hover_positions"

        reval("

        #set object genotypes to hold titles for legend. Set genotype_index and genotype_index to build labeled legend
        genotypes <- c('no call', 'homozygous reference',
            'heterozygous variant', 'homozygous variant')

        genotype_index <- as.matrix(input) + 1
        genotype_text <- matrix(
            paste('Genotype:', genotypes[genotype_index]),
            ncol = ncol(input)
        )

        #prepare chromosome positions hover labels
        d=paste(chrom,pos,sep=',')

        #transpose phenotype matrix to correct automatic transposition in heatmaply function.
        pheno_matrix=t(pheno_matrix)
        #Remove sample names from pheno matrix
        pheno_matrix=pheno_matrix[,-1]
        #transpose trait labels to match dims of pheno_matrix
        trait_labels=t(trait_labels)
        #convert pheno_matrix from type list to matrix
        pheno_matrix2 <- matrix(unlist(pheno_matrix, use.names=FALSE), ncol = pheno_rows, byrow = FALSE)
        #assign sample names as rownames and trait_labels as col_names on pheno_matrix because pheno_matrix is automatically transposed in heatmaply function
        row.names(pheno_matrix2)<-sample_names
        colnames(pheno_matrix2)<-trait_labels

        #assign sample ids to input matrix column names and chromosome positions to row names for automatically generated tick labels
        colnames(input)<-sample_names
        row.names(input)<-d

        suppressPackageStartupMessages(library(heatmaply))

        h=heatmaply(
        plot_method='plotly',
        input,
        dend=FALSE,
        showticklabels = if (x_axis_label_option=='true') {c(TRUE,FALSE)} else {FALSE},
        limits=c(0,3),
        custom_hovertext = genotype_text,
        Rowv=NULL,
        Colv=NULL,
        showticklabels=FALSE,
        col_side_colors=pheno_matrix2,
        col_side_palette = c('1'= 'mediumpurple2',
        '2'='steelblue2'),
        label_names = c('Position', 'Sample ID', 'Genotype Value'),
        main = title,
        xlab = 'Sample IDs',
        ylab = 'Chromosomal Positions',
        key.title = 'Genotype'
        )

        h[['x']][['data']][[length(h[['x']][['data']])]][['marker']][['colorbar']][['ticktext']] <- genotypes

        f<-filename
        if (save_ext == 'html') {htmlwidgets::saveWidget(h,file.path(normalizePath(dirname(f)),basename(f)))} else {withr::with_dir(dirname(f), orca(h, basename(f)))

        }

        ")

    elseif y_axis_label_option == "chromosomes"

        reval("

        #Convert array for chromosome row_side_colors colorbar from list to matrix and name with colname
        chromosome_label_array2 <- matrix(unlist(chromosome_label_array))
        colnames(chromosome_label_array2)<-'Chromosome'

        #set object genotypes to hold titles for legend. Set genotype_index and genotype_index to build labeled legend
        genotypes <- c('no call', 'homozygous reference',
            'heterozygous variant', 'homozygous variant')

        genotype_index <- as.matrix(input) + 1
        genotype_text <- matrix(
            paste('Genotype:', genotypes[genotype_index]),
            ncol = ncol(input)
        )

        #prepare chromosome positions hover labels
        d=paste(chrom,pos,sep=',')

        #transpose phenotype matrix to correct automatic transposition in heatmaply function.
        pheno_matrix=t(pheno_matrix)
        #Remove sample names from pheno matrix
        pheno_matrix=pheno_matrix[,-1]
        #transpose trait labels to match dims of pheno_matrix
        trait_labels=t(trait_labels)
        #convert pheno_matrix from type list to matrix
        pheno_matrix2 <- matrix(unlist(pheno_matrix, use.names=FALSE), ncol = pheno_rows, byrow = FALSE)
        #assign sample names as rownames and trait_labels as col_names on pheno_matrix because pheno_matrix is automatically transposed in heatmaply function
        row.names(pheno_matrix2)<-sample_names
        colnames(pheno_matrix2)<-trait_labels

        #assign sample ids to input matrix column names and chromosome positions to row names for automatically generated tick labels
        colnames(input)<-sample_names
        row.names(input)<-d

        #load heatmaply package and suppress printing start up message
        suppressPackageStartupMessages(library(heatmaply))

        #define heatmap function for grouped genotype matrix with chromosome colorbar
        h=heatmaply(
        input,
        dend=FALSE,
        plot_method='plotly',
        showticklabels= if (x_axis_label_option=='true') {c(TRUE,FALSE)} else {FALSE},
        colorbar_xanchor = 'left',
        colorbar_yanchor = 'top',
        colorbar_xpos = -0.5,
        colorbar_ypos = 0.5,
        row_side_colors=(chromosome_label_array2),
        subplot_widths=c(0.95,0.05),
        subplot_heights=c(0.05,0.95),
        row_side_palette = c(
        '1'= '#d3d7cf',
        '2'='#babdb6',
        '3'='#fce94f',
        '4'='#edd400',
        '5'='#c4a000',
        '6'='#8ae234',
        '7'='#4e9a06',
        '8'='#fcaf3e',
        '9'='#f57900',
        '10'='#ce5c00',
        '11'='#e9b96e',
        '12'='#c17d11',
        '13'='#8f5902',
        '14'='#729fcf',
        '15'='#3465a4',
        '16'='#204a87',
        '17'='#ad7fa8',
        '18'='#75507b',
        '19'='#5c3566',
        '20'='#888a85',
        '21'='#555753',
        '22'='#2e3436',
        'X'='#ef2929',
        'Y'='#cc0000',
        'M'='#a40000'),
        limits=c(0,3),
        custom_hovertext = genotype_text,
        Rowv=NULL,
        Colv=NULL,
        col_side_colors=pheno_matrix2,
        col_side_palette = c('1'='mediumpurple2','2'='steelblue2'),
        label_names = c('Position', 'Sample ID', 'Genotype Value'),
        main = title,
        xlab = 'Sample IDs',
        ylab = 'Chromosomal Positions',
        key.title = 'Genotype'
        )

        #set custom ticktext to colorbar genotype legend
        #h[['x']][['data']][[length(h[['x']][['data']])]][['marker']][['colorbar']][['ticktext']] <- genotypes

        #save file to output directory in correct format
        f<-filename
        if (save_ext == 'html') {htmlwidgets::saveWidget(h,file.path(normalizePath(dirname(f)),basename(f)))} else {withr::with_dir(dirname(f), orca(h, basename(f)))

        }

        ")

    else
        println("--y_axis_labels is not a valid option. Choose positions or hover_positions")
    end

end

"""
    dp_heatmap2(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String}, sample_names,chr_pos_tuple_list_rev,y_axis_label_option,x_axis_label_option,save_ext,chrom_label_info,dp_limit)
generate heatmap of read depth data.
"""
function dp_heatmap2(input,title,filename,sample_names,gt_chromosome_labels,y_axis_label_option,x_axis_label_option,save_ext,chrom_label_info,dp_limit)

    chromosome_label_array=make_chromosome_labels(chrom_label_info)

    chrom=gt_chromosome_labels[:,1]
    pos=gt_chromosome_labels[:,2]

    @rput input
    @rput title
    @rput filename
    @rput sample_names
    @rput chrom
    @rput pos
    @rput save_ext
    @rput chromosome_label_array
    @rput x_axis_label_option
    @rput dp_limit

    if y_axis_label_option == "positions"

        reval("

        d=paste(chrom,pos,sep=',')

        colnames(input)<-sample_names
        row.names(input)<-d

        suppressPackageStartupMessages(library(heatmaply));

        read_depth_colors <- c(
              '#ffffff',
              '#0037ff'
        )

        h=heatmaply(
        #plot_method='plotly',
        input,
        dend=FALSE,
        showticklabels= if (x_axis_label_option=='true') {c(TRUE)} else {c(FALSE,TRUE)},
        colorbar_len=0.7,
        colors = read_depth_colors,
        limits = c(-1, dp_limit),
        Rowv=NULL,
        Colv=NULL,
        label_names = c('Position', 'Sample ID', 'Read Depth'),
        main = title,
        xlab = 'Sample IDs',
        ylab = 'Chromosomal Positions',
        key.title = 'Read Depth'
        )

        f<-filename
        if (save_ext == 'html') {htmlwidgets::saveWidget(h,file.path(normalizePath(dirname(f)),basename(f)))} else {withr::with_dir(dirname(f), orca(h, basename(f)))

        }

        ")

    elseif y_axis_label_option == "hover_positions"

        reval("

        d=paste(chrom,pos,sep=',')

        colnames(input)<-sample_names
        row.names(input)<-d

        suppressPackageStartupMessages(library(heatmaply));

        read_depth_colors <- c(
              '#ffffff',
              '#0037ff'
        )

        #read_depth_colors <- c(
        #      '#ffffff',
        #      '#47aee1',
        #      '#0037ff'
        #)

        h=heatmaply(
        #plot_method='plotly',
        input,
        dend=FALSE,
        colorbar_len=0.7,
        colors = read_depth_colors,
        limits = c(-1, dp_limit),
        showticklabels= if (x_axis_label_option=='true') {c(TRUE,FALSE)} else {FALSE},
        Rowv=NULL,
        Colv=NULL,
        label_names = c('Position', 'Sample ID', 'Read Depth'),
        main = title,
        xlab = 'Sample IDs',
        ylab = 'Chromosomal Positions',
        key.title = 'Read Depth'
        )

        f<-filename
        if (save_ext == 'html') {htmlwidgets::saveWidget(h,file.path(normalizePath(dirname(f)),basename(f)))} else {withr::with_dir(dirname(f), orca(h, basename(f)))

        }

        ")

    elseif y_axis_label_option == "chromosomes"

        reval("

        chromosome_label_array2 <- matrix(unlist(chromosome_label_array))
        colnames(chromosome_label_array2)<-'Chromosome'

        colnames(input)<-sample_names

        suppressPackageStartupMessages(library(heatmaply));

        read_depth_colors <- c(
              '#ffffff',
              '#0037ff'
        )

        h=heatmaply(
        input,
        dend=FALSE,
        plot_method ='plotly',
        colorbar_len=0.7,
        colors = read_depth_colors,
        limits = c(-1, dp_limit),
        colorbar_xanchor = 'left',
        colorbar_yanchor = 'top',
        colorbar_xpos = 1,
        colorbar_ypos = 1,
        showticklabels= if (x_axis_label_option=='true') {c(TRUE,FALSE)} else {FALSE},
        row_side_colors=chromosome_label_array2,
        row_side_palette = c(
        '1'= '#d3d7cf',
        '2'='#babdb6',
        '3'='#fce94f',
        '4'='#edd400',
        '5'='#c4a000',
        '6'='#8ae234',
        '7'='#4e9a06',
        '8'='#fcaf3e',
        '9'='#f57900',
        '10'='#ce5c00',
        '11'='#e9b96e',
        '12'='#c17d11',
        '13'='#8f5902',
        '14'='#729fcf',
        '15'='#3465a4',
        '16'='#204a87',
        '17'='#ad7fa8',
        '18'='#75507b',
        '19'='#5c3566',
        '20'='#888a85',
        '21'='#555753',
        '22'='#2e3436',
        'X'='#ef2929',
        'Y'='#cc0000',
        'M'='#a40000'),
        Rowv=NULL,
        Colv=NULL,
        label_names = c('Position', 'Sample ID', 'Read Depth'),
        main = title,
        xlab = 'Sample IDs',
        ylab = 'Chromosomal Positions',
        key.title = 'Read Depth'
        )

        f<-filename
        if (save_ext == 'html') {htmlwidgets::saveWidget(h,file.path(normalizePath(dirname(f)),basename(f)))} else {withr::with_dir(dirname(f), orca(h, basename(f)))

        }

        ")

    else
        println("--y_axis_labels is not a valid option. Choose positions or hover_positions")
    end

end

#=
read_depth_colors <- c(
      'dark_blue' = '#0037ff',
      'light_blue' = '#f37735',
      'white' = '#ffffff'
)

=#

"""
    dp_heatmap2_with_groups(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev,y_axis_label_option,x_axis_label_option,save_ext,chrom_label_info,dp_limit)
generate heatmap of read depth data with grouped samples.
"""
function dp_heatmap2_with_groups(input,title,filename,sample_names,gt_chromosome_labels,pheno_matrix,trait_labels,y_axis_label_option,x_axis_label_option,save_ext,chrom_label_info,dp_limit)

    chromosome_label_array=make_chromosome_labels(chrom_label_info)

    pheno_rows=size(pheno_matrix,1)-1
    pheno_cols=size(pheno_matrix,2)

    chrom=gt_chromosome_labels[:,1]
    pos=gt_chromosome_labels[:,2]

    @rput input
    @rput title
    @rput filename
    @rput sample_names
    @rput chrom
    @rput pos
    @rput pheno_matrix
    @rput trait_labels
    @rput save_ext
    @rput chromosome_label_array
    @rput x_axis_label_option
    @rput pheno_rows
    @rput pheno_cols
    @rput dp_limit

    if y_axis_label_option == "positions"

        reval("

        #prepare chromosome positions hover labels
        d=paste(chrom,pos,sep=',')

        #transpose phenotype matrix to correct automatic transposition in heatmaply function.
        pheno_matrix=t(pheno_matrix)
        #Remove sample names from pheno matrix
        pheno_matrix=pheno_matrix[,-1]
        #transpose trait labels to match dims of pheno_matrix
        trait_labels=t(trait_labels)
        #convert pheno_matrix from type list to matrix
        pheno_matrix2 <- matrix(unlist(pheno_matrix, use.names=FALSE), ncol = pheno_rows, byrow = FALSE)
        #assign sample names as rownames and trait_labels as col_names on pheno_matrix because pheno_matrix is automatically transposed in heatmaply function
        row.names(pheno_matrix2)<-sample_names
        colnames(pheno_matrix2)<-trait_labels

        #assign sample ids to input matrix column names and chromosome positions to row names for automatically generated tick labels
        colnames(input)<-sample_names
        row.names(input)<-d

        #define read_depth_colors
        read_depth_colors <- c(
              '#ffffff',
              '#0037ff'
        )

        suppressPackageStartupMessages(library(heatmaply))

        h=heatmaply(
        plot_method='plotly',
        input,
        dend=FALSE,
        showticklabels= if (x_axis_label_option=='true') {c(TRUE)} else {c(FALSE,TRUE)},
        colorbar_len=0.7,
        colors = read_depth_colors,
        limits = c(-1, dp_limit),        Rowv=NULL,
        Colv=NULL,
        col_side_colors=pheno_matrix2,
        col_side_palette = c('1'= 'mediumpurple2',
        '2'='steelblue2'),
        label_names = c('Position', 'Sample ID', 'Read Depth'),
        main = title,
        xlab = 'Sample IDs',
        ylab = 'Chromosomal Positions',
        key.title = 'Read Depth'
        )

        f<-filename
        if (save_ext == 'html') {htmlwidgets::saveWidget(h,file.path(normalizePath(dirname(f)),basename(f)))} else {withr::with_dir(dirname(f), orca(h, basename(f)))

        }


        ")

    elseif y_axis_label_option == "hover_positions"

        reval("

        #prepare chromosome positions hover labels
        d=paste(chrom,pos,sep=',')

        #transpose phenotype matrix to correct automatic transposition in heatmaply function.
        pheno_matrix=t(pheno_matrix)
        #Remove sample names from pheno matrix
        pheno_matrix=pheno_matrix[,-1]
        #transpose trait labels to match dims of pheno_matrix
        trait_labels=t(trait_labels)
        #convert pheno_matrix from type list to matrix
        pheno_matrix2 <- matrix(unlist(pheno_matrix, use.names=FALSE), ncol = pheno_rows, byrow = FALSE)
        #assign sample names as rownames and trait_labels as col_names on pheno_matrix because pheno_matrix is automatically transposed in heatmaply function
        row.names(pheno_matrix2)<-sample_names
        colnames(pheno_matrix2)<-trait_labels

        #assign sample ids to input matrix column names and chromosome positions to row names for automatically generated tick labels
        colnames(input)<-sample_names
        row.names(input)<-d

        #define read_depth_colors
        read_depth_colors <- c(
              '#ffffff',
              '#0037ff'
        )

        suppressPackageStartupMessages(library(heatmaply))

        h=heatmaply(
        plot_method='plotly',
        input,
        dend=FALSE,
        showticklabels= if (x_axis_label_option=='true') {c(TRUE,FALSE)} else {FALSE},
        colorbar_len=0.7,
        colors = read_depth_colors,
        limits = c(-1, dp_limit),
        Rowv=NULL,
        Colv=NULL,
        col_side_colors=pheno_matrix2,
        col_side_palette = c('1'= 'mediumpurple2',
        '2'='steelblue2'),
        label_names = c('Position', 'Sample ID', 'Read Depth'),
        main = title,
        xlab = 'Sample IDs',
        ylab = 'Chromosomal Positions',
        key.title = 'Read Depth'
        )

        f<-filename
        if (save_ext == 'html') {htmlwidgets::saveWidget(h,file.path(normalizePath(dirname(f)),basename(f)))} else {withr::with_dir(dirname(f), orca(h, basename(f)))

        }

        ")

    elseif y_axis_label_option == "chromosomes"

        reval("

        #Convert array for chromosome row_side_colors colorbar from list to matrix and name with colname
        chromosome_label_array2 <- matrix(unlist(chromosome_label_array))
        colnames(chromosome_label_array2)<-'Chromosome'

        #prepare chromosome positions hover labels
        d=paste(chrom,pos,sep=',')

        #transpose phenotype matrix to correct automatic transposition in heatmaply function.
        pheno_matrix=t(pheno_matrix)
        #Remove sample names from pheno matrix
        pheno_matrix=pheno_matrix[,-1]
        #transpose trait labels to match dims of pheno_matrix
        trait_labels=t(trait_labels)
        #convert pheno_matrix from type list to matrix
        pheno_matrix2 <- matrix(unlist(pheno_matrix, use.names=FALSE), ncol = pheno_rows, byrow = FALSE)
        #assign sample names as rownames and trait_labels as col_names on pheno_matrix because pheno_matrix is automatically transposed in heatmaply function
        row.names(pheno_matrix2)<-sample_names
        colnames(pheno_matrix2)<-trait_labels

        #assign sample ids to input matrix column names and chromosome positions to row names for automatically generated tick labels
        colnames(input)<-sample_names
        row.names(input)<-d

        #load heatmaply package and suppress printing start up message
        suppressPackageStartupMessages(library(heatmaply))

        #define read_depth_colors
        read_depth_colors <- c(
              '#ffffff',
              '#0037ff'
        )

        h=heatmaply(
        plot_method='plotly',
        input,
        dend=FALSE,
        showticklabels= if (x_axis_label_option=='true') {c(TRUE,FALSE)} else {FALSE},
        colorbar_len=0.7,
        colors = read_depth_colors,
        limits = c(-1, dp_limit),
        colorbar_xanchor = 'left',
        colorbar_yanchor = 'top',
        colorbar_xpos = 1,
        colorbar_ypos = 1,
        row_side_colors=chromosome_label_array2,
        row_side_palette = c(
        '1'= '#d3d7cf',
        '2'='#babdb6',
        '3'='#fce94f',
        '4'='#edd400',
        '5'='#c4a000',
        '6'='#8ae234',
        '7'='#4e9a06',
        '8'='#fcaf3e',
        '9'='#f57900',
        '10'='#ce5c00',
        '11'='#e9b96e',
        '12'='#c17d11',
        '13'='#8f5902',
        '14'='#729fcf',
        '15'='#3465a4',
        '16'='#204a87',
        '17'='#ad7fa8',
        '18'='#75507b',
        '19'='#5c3566',
        '20'='#888a85',
        '21'='#555753',
        '22'='#2e3436',
        'X'='#ef2929',
        'Y'='#cc0000',
        'M'='#a40000'),
        Rowv=NULL,
        Colv=NULL,
        col_side_colors=pheno_matrix2,
        col_side_palette = c('1'= 'mediumpurple2',
        '2'='steelblue2'),
        label_names = c('Position', 'Sample ID', 'Read Depth'),
        main = title,
        xlab = 'Sample IDs',
        ylab = 'Chromosomal Positions',
        key.title = 'Read Depth'
        )

        f<-filename
        if (save_ext == 'html') {htmlwidgets::saveWidget(h,file.path(normalizePath(dirname(f)),basename(f)))} else {withr::with_dir(dirname(f), orca(h, basename(f)))

        }

        ")

    else
        println("--y_axis_labels is not a valid option. Choose positions or hover_positions")
    end

end

#PlotlyJS scatter plots

"""
    avg_sample_dp_scatter(sample_avg_list::Array{Float64,1},sample_names)
generate line chart of average read depths of each sample.
"""
function avg_sample_dp_scatter(sample_avg_list::Array{Float64,1},sample_names)

    sample_name_indices = collect(1:1:size(sample_names,2))

    trace = scatter(;x=1:size(sample_avg_list,1), y=sample_avg_list,mode="markers")
    layout = Layout(title="Average Sample Read Depth",xaxis=attr(title="Samples",ticktext=sample_names,tickvals=sample_name_indices,tickangle=45,tickfont_size=5),yaxis=attr(title="Average Read Depth"),showticklabels=False)
    plot(trace,layout)
end

"""
           avg_variant_dp_line_chart(variant_avg_list::Array{Float64,1},chr_pos_tuple_list)
generate line chart of average read depths of each variant.
"""
function avg_variant_dp_line_chart(variant_avg_list::Array{Float64,1},chr_pos_tuple_list)

    chr_pos_tuple_indices = collect(1:1:size(chr_pos_tuple_list,1))

    trace = scatter(;x=1:size(variant_avg_list,1), y=variant_avg_list,mode="lines")
    layout = Layout(title="Average Variant Read Depth",xaxis=attr(title="Variant Positions",ticktext=chr_pos_tuple_list,tickvals=chr_pos_tuple_indices,tickangle=45,tickfont_size=5,showticklabels=False),yaxis=attr(title="Average Read Depth"))
    plot(trace,layout)
end
