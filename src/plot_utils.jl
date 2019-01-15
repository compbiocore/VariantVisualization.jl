"""
    genotype_heatmap2(input::Array{Any,2},title::AbstractString,chrom_label_info,sample_names,chr_pos_tuple_list_rev,y_axis_label_option,save_ext)
generate heatmap of genotype data.
"""
function genotype_heatmap2(input,title,filename,sample_names,gt_chromosome_labels,y_axis_label_option,save_ext)

    @rput input
    @rput title
    @rput filename
    @rput sample_names

    chrom=gt_chromosome_labels[:,1]
    pos=gt_chromosome_labels[:,2]

    @rput chrom
    @rput pos
    @rput save_ext

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
            input,
            dend=FALSE,
            #limits=c(0,3),
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
            if (save_ext == 'html') {htmlwidgets::saveWidget(h,file.path(normalizePath(dirname(f)),basename(f)))} else {orca(h, file.path(normalizePath(dirname(f)),basename(f)))
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
            limits=c(0,3),
            custom_hovertext = genotype_text,
            Rowv=NULL,
            Colv=NULL,
            showticklabels=FALSE,
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
   genotype_heatmap_with_groups(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev,y_axis_label_option,save_ext)
generate heatmap of genotype data.
"""
function genotype_heatmap_with_groups(input,title,filename,sample_names,gt_chromosome_labels,pheno_matrix,trait_labels,y_axis_label_option,save_ext)

    chrom=gt_chromosome_labels[:,1]
    pos=gt_chromosome_labels[:,2]

    pheno_matrix=pheno_matrix[:,:]
    pheno_rows=size(pheno_matrix,1)

    @rput input
    @rput title
    @rput filename
    @rput sample_names
    @rput chrom
    @rput pos
    @rput pheno_matrix
    @rput pheno_rows
    @rput trait_labels
    @rput save_ext


    if y_axis_label_option == "positions"

        reval("

        genotypes <- c('no call', 'homozygous reference',
            'heterozygous variant', 'homozygous variant')

        genotype_index <- as.matrix(input) + 1
        genotype_text <- matrix(
            paste('Genotype:', genotypes[genotype_index]),
            ncol = ncol(input)
        )

        pheno_matrix=t(pheno_matrix)
        pheno_matrix=pheno_matrix[,-1]
        trait_labels=t(trait_labels)

        if (pheno_rows == 2) {pheno_matrix2 <- matrix(unlist(pheno_matrix))} else {pheno_matrix2 <-pheno_matrix}

        d=paste(chrom,pos,sep=',')

        colnames(input)<-sample_names
        row.names(input)<-d

        row.names(pheno_matrix)<-sample_names
        colnames(pheno_matrix)<-trait_labels

        suppressPackageStartupMessages(library(heatmaply))

        h=heatmaply(
        input,
        dend=FALSE,
        limits=c(0,3),
        custom_hovertext = genotype_text,
        Rowv=NULL,
        Colv=NULL,
        col_side_colors=pheno_matrix,
        col_side_palette = c('1'= 'mediumpurple2',
        '2'='steelblue2'),
        label_names = c('Position', 'Sample ID', 'Genotype Value'),
        main = title,
        xlab = 'Sample IDs',
        ylab = 'Chromosomal Positions',
        key.title = 'Genotype',
        file = filename
        )

        h[['x']][['data']][[length(h[['x']][['data']])]][['marker']][['colorbar']][['ticktext']] <- genotypes

        f<-filename
        if (save_ext == 'html') {htmlwidgets::saveWidget(h,file.path(normalizePath(dirname(f)),basename(f)))} else {orca(h, file.path(normalizePath(dirname(f)),basename(f)))
        }
        ")

    elseif y_axis_label_option == "hover_positions"

        reval("

        genotypes <- c('no call', 'homozygous reference',
            'heterozygous variant', 'homozygous variant')

        genotype_index <- as.matrix(input) + 1
        genotype_text <- matrix(
            paste('Genotype:', genotypes[genotype_index]),
            ncol = ncol(input)
        )

        pheno_matrix=t(pheno_matrix)
        pheno_matrix=pheno_matrix[,-1]
        trait_labels=t(trait_labels)

        if (pheno_rows == 2) {pheno_matrix2 <- matrix(unlist(pheno_matrix))} else {pheno_matrix2 <-pheno_matrix}

        d=paste(chrom,pos,sep=',')

        colnames(input)<-sample_names
        row.names(input)<-d

        row.names(pheno_matrix2)<-sample_names
        colnames(pheno_matrix2)<-trait_labels

        suppressPackageStartupMessages(library(heatmaply))

        h=heatmaply(
        input,
        dend=FALSE,
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
        key.title = 'Genotype',
        file = filename
        )

        h[['x']][['data']][[length(h[['x']][['data']])]][['marker']][['colorbar']][['ticktext']] <- genotypes

        f<-filename
        if (save_ext == 'html') {htmlwidgets::saveWidget(h,file.path(normalizePath(dirname(f)),basename(f)))} else {orca(h, file.path(normalizePath(dirname(f)),basename(f)))
        }

        ")

    else
        println("--y_axis_labels is not a valid option. Choose positions or hover_positions")
    end


end

"""
    dp_heatmap2(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String}, sample_names,chr_pos_tuple_list_rev,y_axis_label_option,save_ext)
generate heatmap of read depth data.
"""
function dp_heatmap2(input,title,filename,sample_names,gt_chromosome_labels,y_axis_label_option,save_ext)

    @rput input
    @rput title
    @rput filename
    @rput sample_names
    @rput save_ext


    chrom=gt_chromosome_labels[:,1]
    pos=gt_chromosome_labels[:,2]

    @rput chrom
    @rput pos

    if y_axis_label_option == "positions"

        reval("

        d=paste(chrom,pos,sep=',')

        colnames(input)<-sample_names
        row.names(input)<-d

        suppressPackageStartupMessages(library(heatmaply));

        h=heatmaply(
        input,
        dend=FALSE,
        scale_fill_gradient_fun = scale_fill_gradient2(low = 'rgb(153,231,255)', high = 'rgb(0, 97, 255)', midpoint = 50, limits = c(0, 100)),
        #limits=c(0,100),
        Rowv=NULL,
        Colv=NULL,
        label_names = c('Position', 'Sample ID', 'Read Depth'),
        main = title,
        xlab = 'Sample IDs',
        ylab = 'Chromosomal Positions',
        key.title = 'Read Depth'
        )

        f<-filename
        if (save_ext == 'html') {htmlwidgets::saveWidget(h,file.path(normalizePath(dirname(f)),basename(f)))} else {orca(h, file.path(normalizePath(dirname(f)),basename(f)))
        }

        ")

    elseif y_axis_label_option == "hover_positions"

        reval("

        d=paste(chrom,pos,sep=',')

        colnames(input)<-sample_names
        row.names(input)<-d

        suppressPackageStartupMessages(library(heatmaply));

        read_depth_colors <- c(
              '#0037ff',
              '#f37735',
              '#ffffff'
        )

        h=heatmaply(
        input,
        dend=FALSE,
        scale_fill_gradient_fun = scale_fill_gradientn(colours = read_depth_colors,values = c(0,0.01,1), space = 'Lab', na.value = 'grey50'),
        limits = c(0, 100),
        showticklabels=FALSE,
        Rowv=NULL,
        Colv=NULL,
        label_names = c('Position', 'Sample ID', 'Read Depth'),
        main = title,
        xlab = 'Sample IDs',
        ylab = 'Chromosomal Positions',
        key.title = 'Read Depth'
        )

        f<-filename
        if (save_ext == 'html') {htmlwidgets::saveWidget(h,file.path(normalizePath(dirname(f)),basename(f)))} else {orca(h, file.path(normalizePath(dirname(f)),basename(f)))
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
    dp_heatmap2_with_groups(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev,y_axis_label_option,save_ext)
generate heatmap of read depth data with grouped samples.
"""
function dp_heatmap2_with_groups(input,title,filename,sample_names,gt_chromosome_labels,pheno_matrix,trait_labels,y_axis_label_option,save_ext)

    chrom=gt_chromosome_labels[:,1]
    pos=gt_chromosome_labels[:,2]

    pheno_rows=size(pheno_matrix,1)

    @rput input
    @rput title
    @rput filename
    @rput sample_names
    @rput chrom
    @rput pos
    @rput pheno_matrix
    @rput pheno_rows
    @rput trait_labels
    @rput save_ext

    if y_axis_label_option == "positions"

        reval("

        pheno_matrix=t(pheno_matrix)
        pheno_matrix=pheno_matrix[,-1]
        trait_labels=t(trait_labels)

        if (pheno_rows == 2) {pheno_matrix2 <- matrix(unlist(pheno_matrix))} else {pheno_matrix2 <-pheno_matrix}

        d=paste(chrom,pos,sep=',')

        colnames(input)<-sample_names
        row.names(input)<-d

        row.names(pheno_matrix2)<-sample_names
        colnames(pheno_matrix2)<-trait_labels

        suppressPackageStartupMessages(library(heatmaply))

        h=heatmaply(
        input,
        dend=FALSE,
        scale_fill_gradient_fun = scale_fill_gradient2(low = 'white', high = 'blue', midpoint = 20,
        limits=c(0,100),
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
        if (save_ext == 'html') {htmlwidgets::saveWidget(h,file.path(normalizePath(dirname(f)),basename(f)))} else {orca(h, file.path(normalizePath(dirname(f)),basename(f)))
        }


        ")

    elseif y_axis_label_option == "hover_positions"

        reval("

        pheno_matrix=t(pheno_matrix)
        pheno_matrix=pheno_matrix[,-1]
        trait_labels=t(trait_labels)

        if (pheno_rows == 2) {pheno_matrix2 <- matrix(unlist(pheno_matrix))} else {pheno_matrix2 <-pheno_matrix}

        d=paste(chrom,pos,sep=',')

        colnames(input)<-sample_names
        row.names(input)<-d

        row.names(pheno_matrix2)<-sample_names
        colnames(pheno_matrix2)<-trait_labels

        (library(heatmaply))

        h=heatmaply(
        input,
        dend=FALSE,
        limits=c(0,100),
        Rowv=NULL,
        Colv=NULL,
        showticklabels=FALSE,
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
        if (save_ext == 'html') {htmlwidgets::saveWidget(h,file.path(normalizePath(dirname(f)),basename(f)))} else {orca(h, file.path(normalizePath(dirname(f)),basename(f)))
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
