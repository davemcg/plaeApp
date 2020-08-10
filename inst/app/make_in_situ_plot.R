## Function to read in original images
get_image <- function(file) {
  image_read(file.path(paste0('www/insitu_layers/ret_',file, ".png")))
}

## Function to recolor each individual cell layer based on expression
recolor <- function(ret_layer, color){
  if (length(color) == 0) {
    recolored_layer <- ret_layer
  } else {
    recolored_layer <- ret_layer %>% image_colorize(75,color)
  }
  return(recolored_layer)
}

## Function to gather gene expression data in table
get_insitu_table <- function(input, db, meta_filter) {

  ### Pull the data for the gene of interest
  gene <- input$insitu_Gene
  grouping_features <- "CellType_predict"

  if (input$insitu_filter_cat !=''){
    validate(
      need(input$insitu_filter_on != '', "Please select at least one feature to filter on")
    )}

  ### Filter expression table, if filters active
  if (input$insitu_filter_cat != ''){
    filt_cat <- input$insitu_filter_cat
    filt_on <- input$insitu_filter_on
    full_table <- db %>% tbl('grouped_stats') %>%
      filter(Gene == gene) %>%
      filter(!!as.symbol(filt_cat) %in% filt_on) %>%
      group_by_at(vars(one_of(c('Gene', grouping_features)))) %>%
      summarise(cpm = sum(cpm * cell_exp_ct) / sum(cell_exp_ct),
                cell_exp_ct = sum(cell_exp_ct, na.rm = TRUE)) %>%
      as_tibble() %>%
      tidyr::drop_na() %>%
      full_join(., meta_filter %>%
                  group_by_at(vars(one_of(grouping_features))) %>%
                  summarise(Count = n())) %>%
      mutate(cell_exp_ct = ifelse(is.na(cell_exp_ct), 0, cell_exp_ct)) %>%
      mutate(`%` = round((cell_exp_ct / Count) * 100, 2),
             Expression = round(cpm * (`%` / 100), 2)) %>%
      select_at(vars(one_of(c('Gene', grouping_features, 'cell_exp_ct', 'Count', '%', 'Expression')))) %>%
      arrange(-Expression)
  }

  ### Or make expression table without filtering if none selected
  else {
    full_table <- db %>% tbl('grouped_stats') %>%
      filter(Gene == gene) %>%
      group_by_at(vars(one_of(c('Gene', grouping_features)))) %>%
      summarise(cpm = sum(cpm * cell_exp_ct) / sum(cell_exp_ct),
                cell_exp_ct = sum(cell_exp_ct, na.rm = TRUE)) %>%
      as_tibble() %>%
      tidyr::drop_na() %>%
      full_join(., meta_filter %>%
                  group_by_at(vars(one_of(grouping_features))) %>%
                  summarise(Count = n())) %>%
      mutate(cell_exp_ct = ifelse(is.na(cell_exp_ct), 0, cell_exp_ct)) %>%
      mutate(`%` = round((cell_exp_ct / Count) * 100, 2),
             Expression = round(cpm * (`%` / 100), 2)) %>%
      select_at(vars(one_of(c('Gene', grouping_features, 'cell_exp_ct', 'Count', '%', 'Expression')))) %>%
      arrange(-Expression)
  }
}


make_insitu_plot <- function(input, scEiaD_2020_v01, meta_filter){
  ### Load unedited images
  amacrine <- get_image('amacrine')
  artery <- get_image('artery')
  astrocyte <- get_image('astrocyte')
  axons <- get_image('axons')
  layer_labels <- get_image('background')
  bipolar <- get_image('bipolar')
  bruch <- get_image('bruch')
  choriocap <- get_image('choriocap')
  cones <- get_image('cones')
  horizontal <- get_image('horizontal')
  cell_labels <- get_image('labels')
  melanocytes <- get_image('melanocytes')
  microglia <- get_image('microglia')
  muller <- get_image('muller')
  rgc <- get_image('rgc')
  rods <- get_image('rods')
  rpe <- get_image('rpe')
  sclera <- get_image('sclera')
  vein <- get_image('vein')

  full_table <- get_insitu_table(input, scEiaD_2020_v01, meta_filter)

  ### Create mini table with color codes
  p <- full_table %>%
    select(CellType_predict,Expression) %>%
    tidyr::drop_na() %>%
    arrange(Expression)

  ### Convert expression to color scale
  p$col <- viridis(length(p$Expression))

  ### Generate legend plot
  leg_lab <- seq(min(p$Expression),max(p$Expression),l=5)
  leg_lab[] <- lapply(leg_lab, round,2)
  legend <- image_graph(width = 300, height = 600, res = 96)
  legend_image <- as.raster(matrix(rev(p$col), ncol=1))
  plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = expression(bold("Log"[2] * "(CPM + 1)")), cex.main=1.5)
  text(x=2, y = seq(0,1,l=5), labels = leg_lab[], cex=1.5)
  rasterImage(legend_image, 0, 0, 1,1)
  dev.off()

  ### Recolor each layer based on expression
  amacrine <- recolor(amacrine, p$col[which(p$CellType_predict == "Amacrine Cells")])
  artery <- recolor(artery,  p$col[which(p$CellType_predict == "Artery")])
  astrocyte <- recolor(astrocyte,  p$col[which(p$CellType_predict == "Astrocytes")])
  axons <- recolor(axons, p$col[which(p$CellType_predict == "Axons")])
  bipolar <- recolor(bipolar, p$col[which(p$CellType_predict == "Bipolar Cells")])
  bruch <- recolor(bruch, p$col[which(p$CellType_predict == "Bruch Membrane")])
  cones <- recolor(cones, p$col[which(p$CellType_predict == "Cones")])
  choriocap <- recolor(choriocap, p$col[which(p$CellType_predict == "Choriocapillaris")])
  horizontal <- recolor(horizontal, p$col[which(p$CellType_predict == "Horizontal Cells")])
  melanocytes <- recolor(melanocytes, p$col[which(p$CellType_predict == "Melanocytes")])
  microglia <- recolor(microglia, p$col[which(p$CellType_predict == "Microglia")])
  muller <- recolor(muller, p$col[which(p$CellType_predict == "Muller Glia")])
  rpe <- recolor(rpe, p$col[which(p$CellType_predict == "RPE")])
  rgc <- recolor(rgc, p$col[which(p$CellType_predict == "Retinal Ganglion Cells")])
  rods <- recolor(rods, p$col[which(p$CellType_predict == "Rods")])
  sclera <- recolor(sclera, p$col[which(p$CellType_predict == "Sclera")])
  vein <- recolor(vein,  p$col[which(p$CellType_predict == "Vein")])

  ### Merge the recolored layers into single image
  retina_insitu <- c(layer_labels, amacrine, artery, astrocyte,choriocap, bipolar, bruch, cones, horizontal, melanocytes, microglia, muller, axons, rpe, rgc, rods, sclera, vein, cell_labels)
  ret_img <- retina_insitu %>%
    image_mosaic() %>%
    image_flatten()

  ### Append the legend to the side and write a temporary file with the complete image
  tmpfile <- image_append(c(ret_img,legend), stack=FALSE) %>%
    image_write(tempfile(fileext='png'), format = 'png')

  ### Return the location of the temporary file
  return(list(src = tmpfile,
              height = input$insitu_height,
              contentType = "image/png"))
}
