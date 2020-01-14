#install.packages('BiocManager')
#BiocManager::install(version = "3.10")
#BiocManager::install(c("rhdf5"))
#devtools::install_github("pachterlab/sleuth")

packages <- c('sleuth','ggplot2','magrittr','dplyr','amap')
load.status <- lapply(packages, require, character.only = TRUE)
mapply(list, packages, load.status)

setwd(paste0('/Users/rensholmer/Documents/',
             'Data/transcriptome/parasponia',
             '/noot_samples/'))

filter.fun <- function (row, min.reads = 5, min.prop = .23)
{
  mean(row >= min.reads) >= min.prop
}

format.genename <- function(name)
{
  sub('.1', '',sub('_asm01_ann01','',name), fixed=TRUE)
}

find.patterns <- function(genes, name, dat, palette='Set3')
{
  clust.data <- dplyr::filter(dat, target_id %in% genes) %>%
    tidyr::unite(sample_, genotype, tissue) %>%
    dplyr::group_by(target_id) %>%
    dplyr::mutate(sample = make.unique(sample_)) %>%
    dplyr::ungroup() %>%
    reshape2::dcast(target_id ~ sample, value.var = 'tpm')

  rownames(clust.data) = clust.data[,1]
  clust.data = clust.data[,-1]

  transformed.data <- t(apply(clust.data,1,function(gene){
    return(gene / sum(gene))
  }))

  set.seed(3)
  k <- 2
  c <- amap::Kmeans(transformed.data, k, nstart = 10, iter.max = 100)

  centers <- data.frame(t(c$centers))
  colnames(centers) <- paste0(
    rep('pattern ', k),
    seq(k),
    ' ',
    rep('(n=',k),
    c$size,
    rep(')',k)
  )
  centers$genotype <- sapply(strsplit(rownames(centers), split = '_'), function(x){ return(x[1]) })
  centers$tissue <- sapply(strsplit(rownames(centers), split = '_'), function(x){
    return(strsplit(x[2], split = '.', fixed = TRUE)[[1]][1])
  })

  plot.dat <- reshape2::melt(centers)

  p <- ggplot(plot.dat, aes(x = tissue, y = value, group = tissue, fill = tissue)) +
    geom_bar( stat = 'summary', fun.y = 'mean',
              color = 'black', width = .7, position = 'dodge') +
    scale_fill_brewer(palette = palette) +
    geom_errorbar(stat = 'summary', fun.data = 'mean_se',
                  position = position_dodge(width = .5), width = .4) +
    geom_point(shape = 21, fill = 'white', size = 3,
               position = position_dodge(width = .5)) +
    facet_grid(variable ~ genotype) +
    theme_linedraw() +
    labs(
      title = paste0(name,' responsive genes'),
      y = 'TPM',
      x = 'Tissue'
    )
  ggsave(paste0('internode.', name, '.patterns.pdf'), device='pdf')

  write.table(file=paste0('internode.',name,'_pattern1.genes.txt'),
              x=format.genename(names(c$cluster[c$cluster == 1])),
              quote=FALSE, row.names=FALSE, col.names=FALSE)

  write.table(file=paste0('internode.',name,'_pattern2.genes.txt'),
              x=format.genename(names(c$cluster[c$cluster == 2])),
              quote=FALSE, row.names=FALSE, col.names=FALSE)
  p
}

s2c <- data.frame(sample = seq(1, 39),
                  genotype = c(
                    rep('A5', 9), rep('A10', 9), rep('ctr44', 9),
                    rep('noot1', 6), rep('ctr44', 6)
                  ),
                  noot = c(
                    rep('noot1', 18), rep('ctr44', 9),
                    rep('noot1', 6), rep('ctr44', 6)
                  ),
                  tissue = c(
                    rep(c(rep('in4', 3), rep('in2', 3), rep('in3', 3)), 3),
                    rep(c(rep('nodule', 3), rep('root', 3)), 2)
                  ),
                  path = as.character(seq(1, 39)),
                  stringsAsFactors = FALSE)

nodule.samples <- subset(s2c, sample >= 28)
stem.samples <- subset(s2c, sample <= 27)


so <- sleuth_prep(stem.samples, filter_fun = filter.fun)

so <- sleuth_fit(so, ~1, 'intercept')
so <- sleuth_fit(so, ~genotype, 'genotype')
so <- sleuth_fit(so, ~tissue, 'tissue')
so <- sleuth_fit(so, ~genotype + tissue, 'additive')
so <- sleuth_fit(so, ~genotype * tissue, 'interaction')

so <- sleuth_lrt(so, 'additive', 'interaction')
so <- sleuth_lrt(so, 'intercept', 'tissue')
so <- sleuth_lrt(so, 'intercept', 'genotype')
so <- sleuth_lrt(so, 'genotype','interaction')

interaction.lrt <- sleuth_results(so, 'additive:interaction', 'lrt', show_all = TRUE)
interaction.sig <- dplyr::filter(interaction.lrt, qval <= 0.01)

tissue.lrt <- sleuth_results(so, 'intercept:tissue', 'lrt', show_all = TRUE)
tissue.sig <- dplyr::filter(tissue.lrt, qval <= 0.01)

genotype.lrt <- sleuth_results(so, 'intercept:genotype', 'lrt',
                               show_all = TRUE)
genotype.sig <- dplyr::filter(genotype.lrt, qval <= 0.01)

dat <- kallisto_table(so)

find.patterns(genotype.sig$target_id, 'genotype', dat, 'Blues')
find.patterns(tissue.sig$target_id, 'tissue', dat, 'Greens')

gene_id <- tissue.sig$target_id[1]

d <- dplyr::filter(dat, target_id == gene_id)

ggplot(d, aes(x = target_id, y = tpm, group = tissue, fill = tissue)) +
  geom_bar( stat = 'summary', fun.y = 'mean',
            color = 'black', width = .5, position = 'dodge') +
  geom_errorbar(stat = 'summary', fun.data = 'mean_se',
                position = position_dodge(width = .5), width = .4) +
  geom_point(shape = 21, fill = 'white', size = 3,
             position = position_dodge(width = .5)) +
  facet_grid(~genotype) +
  theme_linedraw() +
  theme(
    legend.text = element_text(),
    axis.text.x = element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank()) +
  ylab("TPM") +
  xlab(gene_id)





