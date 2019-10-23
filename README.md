# pipeComp

`pipeComp` is a simple framework to facilitate the comparison of pipelines involving various steps and parameters. Given a `PipelineDefinition`, a set of alternative parameters (which might include different subroutines) and benchmark datasets, the `runPipeline` function proceeds through all combinations arguments, avoiding recomputing the same step twice and compiling evaluations on the fly to avoid storing potentially large intermediate data.

## PipelineDefinition

The `PipelineDefinition` S4 class represents pipelines as, minimally, a set of functions consecutively executed on the output of the previous one, and optionally accompanied by evaluation and aggregation functions. As simple pipeline can be constructed as follows:

```{r, eval=FALSE}
my_pip <- PipelineDefinition( list( step1=function(x, param1){
                                      # do something with x and param1
                                      x
                                    },
                                    step2=function(x, method1, param2){
                                      get(method1)(x, param2)
                                    },
                                    step3=function(x, param3){
                                      x <- some_fancy_function(x, param3)
                                      # the functions can also output evaluation
                                      # through the `intermediate_return` slot:
                                      e <- my_evaluation_function(x)
                                      list( x=x, intermediate_return=e)
                                    }
                                  ))
```

The PipelineDefinition can also include descriptions of each step or evaluation and aggregation functions. See the `?PipelineDefinition` for more information, or `scrna_seurat_pipeline` for a more complex example:

```{r}
pipDef <- scrna_seurat_pipeline()
pipDef
```
```
A PipelineDefinition object with the following steps:
  - doublet(x, doubletmethod) *
Takes a SCE object with the `phenoid` colData column, passes it through the 
function `doubletmethod`, and outputs a filtered SCE.
  - filtering(x, filt) *
Takes a SCE object, passes it through the function `filt`, and outputs a 
filtered Seurat object.
  - normalization(x, norm)
Passes the object through function `norm` to return the object with the 
normalized and scale data slots filled.
  - selection(x, sel, selnb)
Returns a seurat object with the VariableFeatures filled with `selnb` features 
using the function `sel`.
  - dimreduction(x, dr, maxdim) *
Returns a seurat object with the PCA reduction with up to `maxdim` components
using the `dr` function.
  - clustering(x, clustmethod, dims, k, steps, resolution, min.size) *
Uses function `clustmethod` to return a named vector of cell clusters.
```

A number of generic methods are implemented on the object, including `show`, `names`, `length`, `[`, `as.list`.

## Running pipelines

### Preparing the other arguments

`runPipeline` requires 3 main arguments: i) the pipelineDefinition, ii) the list of alternative parameters values to try, and iii) the list of benchmark datasets.

Functions can be passed as arguments through their name (if they are loaded in the environment).

```{r}
alternatives <- list(
  doubletmethod=c("none", "scDblFinder_wrapper"),
  filt=c("lenientfilter"),
  norm=c("norm.seurat", "norm.seuratvst", "norm.scran"),
  sel=c("sel.vst"),
  selnb=2000,
  dr=c("seurat.pca"),
  clustmethod=c("clust.seurat"),
  maxdim=30,
  dims=c(10, 15, 20, 30),
  k=20,
  steps=4,
  resolution=c(0.01, 0.1, 0.2, 0.3, 0.5, 0.8, 1, 1.2, 2),
  min.size=50
)
```

### Running the analyses

```{r}
res <- runPipeline( datasets, alternatives, pipDef, nthreads=3,
                    output.prefix="myfolder/" )
```

### Exploring the metrics

Data can be explored manually or plotted using wrapper functions designed for each step at which benchmark data is gathered. For example:

```{r}
scrna_evalPlot_DR(res, scale=FALSE)
```

<img src="inst/docs/dr_stats_example.png"/>

```{r}
scrna_evalPlot_clust(res)
```

<img src="inst/docs/clust_stats_example.png"/>