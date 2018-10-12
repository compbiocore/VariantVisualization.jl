
<a id='Heatmap-functions-documentation-1'></a>

# Heatmap functions documentation


Specify output directory for saving all graphics. If directory doesn't exist, it creates the directory within the working directory. Defaults to "output."     "–output_directory", "-o"     ... –line_chart [option] ...


```
Specify file format you wish to save all graphics as (eg. pdf, html, png). [REQUIRED]
"--save_format", "-s"    
... -s [option] ...
Possible options:
html
pdf
svg
png   

Generate heatmap of distribution of values in genotype field across samples. Defaults to genotype.
"--heatmap", "-m"
... -m [option] ...
Possible options:
genotype
read_depth

Specify filename for saving heatmap. #test this, positional argument?
"heatmap_title"

Generate scatter plots of average read depths across either samples or variants.
... --line_chart [option] ...
Possible options:
sample
variant
```

