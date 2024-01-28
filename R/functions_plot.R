#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_plot.R  
#' @description R script containing all functions relative to data
#               visualisation
#' @author Julien Barrere
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Plot spatial distribution of all species tested
#' @param data_gbif gbif occurence data for all species studied
#' @param file.in name (including path) of the file to save
plot_species_distribution <- function(data_gbif, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # List of all species present
  species.in <- unique(data_gbif$species)
  species.in <- species.in[order(species.in)]
  
  # Initialize final list
  list.plots <- list()
  
  # Loop on all species
  for(i in 1:length(species.in)){
    # Print species
    print(paste0("Making plot for ", species.in[i]))
    
    # Data for plotting
    data.i <- data_gbif %>%
      filter(species == species.in[i]) %>%
      dplyr::select(obsid = occurrenceid, species, longitude = decimallongitude, 
                    latitude = decimallatitude) %>%
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
    # Plot for species i
    plot.i <- ne_countries(scale = "medium", returnclass = "sf") %>%
      ggplot(aes(geometry = geometry)) +
      geom_sf(show.legend = F, fill = "#242423", size = 0.01) +
      geom_sf(data = data.i, 
              shape = 16, size = 0.1, color = "red") + 
      coord_sf(ylim = c(-56, 80), expand = TRUE) + 
      theme(panel.background = element_rect(color = 'black', fill = 'white'), 
            panel.grid = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank(), 
            plot.title = element_text(hjust = 0.5)) + 
      ggtitle(species.in[i])
    # Add to the list
    eval(parse(text = paste0("list.plots$", gsub("\\ ", "\\.", species.in[i]), " <- plot.i")))
  }
  
  # Make the final plot
  plot.out <- cowplot::plot_grid(plotlist = list.plots)
  
  # Save the plot
  ggsave(file.in, plot.out, width = 40, height = 20, units = "cm", dpi = 600)
  
  # Return the name of the file saved
  return(file.in)
  
}


#' Plot spatial distribution of all species tested through the density of gbif obs
#' @param data_gbif gbif occurence data for all species studied
#' @param file.in name (including path) of the file to save
plot_species_distribution_density  <- function(data_gbif, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # List of all species present
  species.in <- unique(data_gbif$species)
  species.in <- species.in[order(species.in)]
  
  # Country shapefile
  countries <- ne_countries(scale = "medium", returnclass = "sf") 
  
  # Grid of the world large enough to include all occurence points
  grid <- st_make_grid(
    (data_gbif %>%
       dplyr::select(obsid = occurrenceid, species, longitude = decimallongitude, 
                     latitude = decimallatitude) %>%
       st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant")), 
    n = c(200, 200)
  ) %>%
    st_intersection(countries) %>%
    st_sf() %>%
    mutate(ID = as.character(c(1:dim(.)[1])))
  
  # Initialize final list
  list.plots <- list()
  
  # Loop on all species
  for(i in 1:length(species.in)){
    # Print species
    print(paste0("Making plot for ", species.in[i]))
    
    # Data for plotting
    data.i <- data_gbif %>%
      filter(species == species.in[i]) %>%
      dplyr::select(obsid = occurrenceid, species, longitude = decimallongitude, 
                    latitude = decimallatitude) %>%
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
    
    # Assign each tree to a grid cell
    tree_in_grid.i <- st_join(data.i, grid, join = st_within)
    
    # Calculate the number of trees per grid cell
    tree_per_cell.i <- count(as_tibble(tree_in_grid.i), ID) %>%
      filter(!is.na(ID))
    
    # Add tree number to the grid, and calculate density
    grid.i <- grid %>%
      left_join(tree_per_cell.i, by = "ID") %>%
      filter(!is.na(n)) %>%
      mutate(area = as.numeric(st_area(.))) %>%
      mutate(density = n/area)
    
    # Scaling for the color gradient
    qn = quantile(grid.i$density, c(0:19)/19, na.rm = TRUE)
    qn01 <- as.numeric(scales::rescale(c(qn, range(grid.i$density))))[c(1:(length(qn)))]
    
    # Plot for species i
    plot.i <- ne_countries(scale = "medium", returnclass = "sf") %>%
      ggplot(aes(geometry = geometry)) +
      geom_sf(show.legend = F, fill = "#242423", size = 0.01) +
      geom_sf(data = grid.i,
              color = "NA", aes(fill = density)) + 
      scale_fill_gradientn(colours = colorRampPalette(c("#F1A7A9", "#DD2C2F", "#9C191B"))(20),
                           values = qn01) + 
      coord_sf(ylim = c(-56, 80), expand = TRUE) + 
      theme(panel.background = element_rect(color = 'black', fill = 'white'), 
            panel.grid = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank(), 
            plot.title = element_text(hjust = 0.5), 
            legend.position = "none") + 
      ggtitle(species.in[i])
    
    # Add to the list
    eval(parse(text = paste0("list.plots$", gsub("\\ ", "\\.", species.in[i]), " <- plot.i")))
  }
  
  # Make the final plot
  plot.out <- cowplot::plot_grid(plotlist = list.plots, ncol = 5)
  
  # Save the plot
  ggsave(file.in, plot.out, width = 28, height = 26, units = "cm", dpi = 600)
  
  # Return the name of the file saved
  return(file.in)
  
}






#' Plot the native distribution blocks
#' @param file.in name (including path) of the file to save
plot_native_distribution_blocks <- function(file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Make the final plot
  plot.out <- ne_countries(scale = "medium", returnclass = "sf") %>%
    ggplot(aes(geometry = geometry)) +
    geom_sf(show.legend = F, color = "#6C757D", fill = "#CED4DA", size = 0.01) +
    geom_rect(aes(xmin = -10, xmax = 160, ymin = 25, ymax = 75), color = "#22577A", fill = "NA") +
    geom_text(aes(x = 75, y =  80, label = "Eurasian native block"), color = "#22577A") +
    geom_rect(aes(xmin = -160, xmax = -30, ymin = 0, ymax = 80), color = "#386641", fill = "NA") +
    geom_text(aes(x = -95, y =  -5, label = "American native block"), color = "#386641") +
    geom_rect(aes(xmin = 110, xmax = 154, ymin = -45, ymax = -8), color = "#540B0E", fill = "NA") +
    geom_text(aes(x = 131, y =  -3, label = "Australian native block"), color = "#540B0E")  + 
    coord_sf(ylim = c(-56, 80), expand = TRUE) + 
    xlab("") + ylab("") +
    theme(panel.background = element_rect(color = 'black', fill = 'white'), 
          panel.grid = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          plot.title = element_text(hjust = 0.5))
  
  # Save the plot
  ggsave(file.in, plot.out, width = 20, height = 8, units = "cm", dpi = 600)
  
  # Return the name of the file saved
  return(file.in)
  
}
