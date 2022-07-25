library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(Seurat)
library(miniUI)
library(layer)
library(grid)
library(Cairo)
library(shiny)

# object <- readRDS("u1.rds")
object <- U1.Seurat
feature = "ACTA2"
slot = 'data'
alpha = c(0.1,1)

SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

# Feature plot palettes
#
FeaturePalettes <- list(
  'Spatial' = SpatialColors(n = 100),
  'Seurat' = c('lightgrey', 'blue')
)

# Set inital data values
message("Finding assay keys")
assay.keys <- Key(object = object)[Seurat::Assays(object = object)]
keyed <- sapply(X = assay.keys, FUN = grepl, x = feature)
assay <- if (any(keyed)) {
  names(x = which(x = keyed))[1]
} else {
  DefaultAssay(object = object)
}
message("Finding features")
features <- sort(x = rownames(x = GetAssayData(
  object = object,
  slot = slot,
  assay = assay
)))
feature.label <- 'Feature to visualize'
message("Finding Assays to use")
assays.use <- vapply(
  X = Seurat::Assays(object = object),
  FUN = function(x) {
    return(!IsMatrixEmpty(x = GetAssayData(
      object = object,
      slot = slot,
      assay = x
    )))
  },
  FUN.VALUE = logical(length = 1L)
)
assays.use <- sort(x = Seurat::Assays(object = object)[assays.use])
# Setup gadget UI
ui <- miniPage(
  miniButtonBlock(miniTitleBarButton(
    inputId = 'done',
    label = 'Done',
    primary = TRUE
  )),
  miniContentPanel(
    fillRow(
      sidebarPanel(
        sliderInput(
          inputId = 'alpha',
          label = 'Alpha intensity',
          min = 0,
          max = max(alpha),
          value = min(alpha),
          step = 0.01,
          width = '100%'
        ),
        sliderInput(
          inputId = 'pt.size',
          label = 'Point size',
          min = 0,
          max = 5,
          value = 1.6,
          step = 0.1,
          width = '100%'
        ),
        selectInput(
          inputId = 'assay',
          label = 'Assay',
          choices = assays.use,
          selected = assay,
          selectize = FALSE,
          width = '100%'
        ),
        selectInput(
          inputId = 'feature',
          label = feature.label,
          choices = features,
          selected = feature,
          selectize = FALSE,
          width = '100%'
        ),
        selectInput(
          inputId = 'palette',
          label = 'Color scheme',
          choices = names(x = FeaturePalettes),
          selected = 'Spatial',
          selectize = FALSE,
          width = '100%'
        ),
        width = '100%'
      ),
      plotOutput(outputId = 'plot', height = '100%'),
      flex = c(1, 4)
      
    )
  )
)
# Prepare plotting data
message("Prepare plotting data...")
image <- object@images$slice1
cells.use <- Cells(object)
coords <- GetTissueCoordinates(image)
feature.data <- FetchData(
  object = object,
  vars = feature,
  cells = cells.use,
  slot = slot
)
plot.data <- cbind(coords, feature.data)
server <- function(input, output, session) {
  message("Initialize reactive values")
  plot.env <- reactiveValues(
    data = plot.data,
    feature = feature,
    palette = 'Spatial'
  )
  # Observe events
  observeEvent(
    eventExpr = input$done,
    handlerExpr = stopApp(returnValue = plot.env$plot)
  )
  observe(x = {
    assay <- input$assay
    feature.use <- input$feature
    features.assay <- sort(x = rownames(x = GetAssayData(
      object = object,
      slot = slot,
      assay = assay
    )))
    feature.use <- ifelse(
      test = feature.use %in% features.assay,
      yes = feature.use,
      no = features.assay[1]
    )
    updateSelectInput(
      session = session,
      inputId = 'assay',
      label = 'Assay',
      choices = assays.use,
      selected = assay
    )
    updateSelectInput(
      session = session,
      inputId = 'feature',
      label = feature.label,
      choices = features.assay,
      selected = feature.use
    )
  })
  observe(x = {
    feature.use <- input$feature
    try(
      expr = {
        feature.data <- FetchData(
          object = object,
          vars = paste0(Key(object = object[[input$assay]]), feature.use),
          cells = cells.use,
          slot = slot
        )
        colnames(x = feature.data) <- feature.use
        plot.env$data <- cbind(coords, feature.data)
        plot.env$feature <- feature.use
      },
      silent = TRUE
    )
  })
  observe(x = {
    plot.env$palette <- input$palette
  })
  # Create plot
  message("Creating Plot")
  output$plot <- renderPlot(expr = {
    plot.env$plot <- SingleSpatialPlot(
      data = plot.env$data,
      image = image,
      col.by = plot.env$feature,
      pt.size.factor = input$pt.size,
      crop = TRUE,
      alpha.by = plot.env$feature
    ) +
      # scale_fill_gradientn(name = plot.env$feature, colours = cols) +
      scale_fill_gradientn(name = plot.env$feature, colours = FeaturePalettes[[plot.env$palette]]) +
      theme(legend.position = 'top') +
      scale_alpha(range = c(input$alpha, 1)) +
      guides(alpha = FALSE)
    plot.env$plot
  })
}

SingleSpatialPlot <- function(
  data,
  image,
  cols = NULL,
  image.alpha = 1,
  pt.alpha = NULL,
  crop = TRUE,
  pt.size.factor = NULL,
  stroke = 0.25,
  col.by = NULL,
  alpha.by = NULL,
  cells.highlight = NULL,
  cols.highlight = c('#DE2D26', 'grey50'),
  geom = c('spatial', 'interactive', 'poly'),
  na.value = 'grey50'
) {
  geom <- match.arg(arg = geom)
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find '", col.by, "' in data, not coloring", call. = FALSE, immediate. = TRUE)
    col.by <- NULL
  }
  col.by <- col.by %iff% paste0("`", col.by, "`")
  alpha.by <- alpha.by %iff% paste0("`", alpha.by, "`")
  if (!is.null(x = cells.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = pt.size.factor,
      cols.highlight = cols.highlight[1],
      col.base = cols.highlight[2]
    )
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- 'highlight'
    levels(x = data$ident) <- c(order, setdiff(x = levels(x = data$ident), y = order))
    data <- data[order(data$ident), ]
  }
  plot <- ggplot(data = data, aes_string(
    x = colnames(x = data)[2],
    y = colnames(x = data)[1],
    fill = col.by,
    alpha = alpha.by
  ))
  plot <- switch(
    EXPR = geom,
    'spatial' = {
      if (is.null(x = pt.alpha)) {
        plot <- plot + geom_spatial(
          point.size.factor = pt.size.factor,
          data = data,
          image = image,
          image.alpha = image.alpha,
          crop = crop,
          stroke = stroke,
        )
      } else {
        plot <- plot + geom_spatial(
          point.size.factor = pt.size.factor,
          data = data,
          image = image,
          image.alpha = image.alpha,
          crop = crop,
          stroke = stroke,
          alpha = pt.alpha
        )
      }
      plot + coord_fixed() + theme(aspect.ratio = 1)
    },
    'interactive' = {
      plot + geom_spatial_interactive(
        data = tibble(grob = list(GetImage(object = image, mode = 'grob'))),
        mapping = aes_string(grob = 'grob'),
        x = 0.5,
        y = 0.5
      ) +
        geom_point(mapping = aes_string(color = col.by)) +
        xlim(0, ncol(x = image)) +
        ylim(nrow(x = image), 0) +
        coord_cartesian(expand = FALSE)
    },
    'poly' = {
      data$cell <- rownames(x = data)
      data[, c('x', 'y')] <- NULL
      data <- merge(
        x = data,
        y = GetTissueCoordinates(object = image, qhulls = TRUE),
        by = "cell"
      )
      plot + geom_polygon(
        data = data,
        mapping = aes_string(fill = col.by, group = 'cell')
      ) + coord_fixed() + theme_cowplot()
      
    },
    stop("Unknown geom, choose from 'spatial' or 'interactive'", call. = FALSE)
  )
  if (!is.null(x = cells.highlight)) {
    plot <- plot + scale_fill_manual(values = cols.highlight)
  }
  if (!is.null(x = cols) && is.null(x = cells.highlight)) {
    if (length(x = cols) == 1 && (is.numeric(x = cols) || cols %in% rownames(x = brewer.pal.info))) {
      scale <- scale_fill_brewer(palette = cols, na.value = na.value)
    } else if (length(x = cols) == 1 && (cols %in% c('alphabet', 'alphabet2', 'glasbey', 'polychrome', 'stepped'))) {
      colors <- DiscretePalette(length(unique(data[[col.by]])), palette = cols)
      scale <- scale_fill_manual(values = colors, na.value = na.value)
    } else {
      scale <- scale_fill_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot + NoAxes() + theme(panel.background = element_blank())
  return(plot)
}

GeomSpatialInteractive <- ggproto(
  "GeomSpatialInteractive",
  Geom,
  setup_data = function(self, data, params) {
    data <- ggproto_parent(parent = Geom, self = self)$setup_data(data, params)
    data
  },
  draw_group = function(data, panel_scales, coord) {
    vp <- viewport(x = data$x, y = data$y)
    g <- editGrob(grob = data$grob[[1]], vp = vp)
    # Replacement for ggname
    g$name <- grobName(grob = g, prefix = 'geom_spatial_interactive')
    return(g)
    # return(ggname(prefix = "geom_spatial", grob = g))
  },
  required_aes = c("grob","x","y")
)

geom_spatial_interactive <-  function(
  mapping = NULL,
  data = NULL,
  stat = "identity",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = FALSE,
  ...
) {
  layer(
    geom = GeomSpatialInteractive,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

GeomSpatial <- ggproto(
  "GeomSpatial",
  Geom,
  required_aes = c("x", "y"),
  extra_params = c("na.rm", "image", "image.alpha", "crop"),
  default_aes = aes(
    shape = 21,
    colour = "black",
    point.size.factor = 1.0,
    fill = NA,
    alpha = NA,
    stroke = 0.25
  ),
  setup_data = function(self, data, params) {
    data <- ggproto_parent(Geom, self)$setup_data(data, params)
    # We need to flip the image as the Y coordinates are reversed
    data$y = max(data$y) - data$y + min(data$y)
    data
  },
  draw_key = draw_key_point,
  draw_panel = function(data, panel_scales, coord, image, image.alpha, crop) {
    # This should be in native units, where
    # Locations and sizes are relative to the x- and yscales for the current viewport.
    if (!crop) {
      y.transform <- c(0, nrow(x = image)) - panel_scales$y.range
      data$y <- data$y + sum(y.transform)
      panel_scales$x$continuous_range <- c(0, ncol(x = image))
      panel_scales$y$continuous_range <- c(0, nrow(x = image))
      panel_scales$y.range <- c(0, nrow(x = image))
      panel_scales$x.range <- c(0, ncol(x = image))
    }
    z <- coord$transform(
      data.frame(x = c(0, ncol(x = image)), y = c(0, nrow(x = image))),
      panel_scales
    )
    # Flip Y axis for image
    z$y <- -rev(z$y) + 1
    wdth <- z$x[2] - z$x[1]
    hgth <- z$y[2] - z$y[1]
    vp <- viewport(
      x = unit(x = z$x[1], units = "npc"),
      y = unit(x = z$y[1], units = "npc"),
      width = unit(x = wdth, units = "npc"),
      height = unit(x = hgth, units = "npc"),
      just = c("left", "bottom")
    )
    img.grob <- GetImage(object = image)
    
    img <- editGrob(grob = img.grob, vp = vp)
    # spot.size <- slot(object = image, name = "spot.radius")
    spot.size <- Radius(object = image)
    coords <- coord$transform(data, panel_scales)
    pts <- pointsGrob(
      x = coords$x,
      y = coords$y,
      pch = data$shape,
      size = unit(spot.size, "npc") * data$point.size.factor,
      gp = gpar(
        col = alpha(colour = coords$colour, alpha = coords$alpha),
        fill = alpha(colour = coords$fill, alpha = coords$alpha),
        lwd = coords$stroke)
    )
    vp <- viewport()
    gt <- gTree(vp = vp)
    if (image.alpha > 0) {
      if (image.alpha != 1) {
        img$raster = as.raster(
          x = matrix(
            data = alpha(colour = img$raster, alpha = image.alpha),
            nrow = nrow(x = img$raster),
            ncol = ncol(x = img$raster),
            byrow = TRUE)
        )
      }
      gt <- addGrob(gTree = gt, child = img)
    }
    gt <- addGrob(gTree = gt, child = pts)
    # Replacement for ggname
    gt$name <- grobName(grob = gt, prefix = 'geom_spatial')
    return(gt)
    # ggplot2:::ggname("geom_spatial", gt)
  }
)

# influenced by: https://stackoverflow.com/questions/49475201/adding-tables-to-ggplot2-with-facet-wrap-in-r
# https://ggplot2.tidyverse.org/articles/extending-ggplot2.html
#' @importFrom ggplot2 layer
#'
#'
geom_spatial <-  function(
  mapping = NULL,
  data = NULL,
  image = image,
  image.alpha = image.alpha,
  crop = crop,
  stat = "identity",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  ...
) {
  layer(
    geom = GeomSpatial,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, image = image, image.alpha = image.alpha, crop = crop, ...)
  )
}

shinyApp(ui = ui, server = server)