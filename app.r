#library(shiny)
library(rgbif)
library(sf)
library(dbscan)
library(leaflet)
library(ape)
library(dplyr)
library(geojsonio)
library(htmltools)
library(shinyWidgets)
library(maps)
library(RColorBrewer)
library(picante)
library(ggtree)
library(ggplot2)
library(maptools)
library(rmapshaper)
library(raster)
library(geodata)



# Function to get phylogeny
get_phy <- function(phy_choice) {
  if (phy_choice == "ALLMB") {
    phy <- 'https://github.com/FePhyFoFum/big_seed_plant_trees/releases/download/v0.1/v0.1.zip'
    download.file(phy, file.path(tempdir(), 'v01.zip'))
    unzip(file.path(tempdir(), 'v01.zip'), exdir = tempdir())
    tree <- ape::read.tree(file.path(tempdir(), 'v0.1', 'ALLMB.tre'))
  } else if (phy_choice == "ALLOTB") {
    phy <- 'https://github.com/FePhyFoFum/big_seed_plant_trees/releases/download/v0.1/v0.1.zip'
    download.file(phy, file.path(tempdir(), 'v01.zip'))
    unzip(file.path(tempdir(), 'v01.zip'), exdir = tempdir())
    tree <- ape::read.tree(file.path(tempdir(), 'v0.1', 'ALLOTB.tre'))
  } else if (phy_choice == "GBMB") {
    phy <- 'https://github.com/FePhyFoFum/big_seed_plant_trees/releases/download/v0.1/v0.1.zip'
    download.file(phy, file.path(tempdir(), 'v01.zip'))
    unzip(file.path(tempdir(), 'v01.zip'), exdir = tempdir())
    tree <- ape::read.tree(file.path(tempdir(), 'v0.1', 'GBMB.tre'))
  } else if (phy_choice == "GBOTB") {
    phy <- 'https://github.com/FePhyFoFum/big_seed_plant_trees/releases/download/v0.1/v0.1.zip'
    download.file(phy, file.path(tempdir(), 'v01.zip'))
    unzip(file.path(tempdir(), 'v01.zip'), exdir = tempdir())
    tree <- ape::read.tree(file.path(tempdir(), 'v0.1', 'GBOTB.tre'))
  } else {
    stop("Invalid choice. Please specify 'ALLMB', 'ALLOTB', 'GBMB', or 'GBOTB'.")
  }
  
  return(tree)
}

# Function to subset the tree

subset_tree <- function(phy, index_of_interest, num_species = 200) {
  # Calculate the start and end indices for the subset
  start_index <- max(1, index_of_interest - floor(num_species/2))
  end_index <- min(length(phy$tip.label), start_index + num_species - 1)
  
  # Get the tip labels to keep
  tips_to_keep <- phy$tip.label[start_index:end_index]
  
  # Drop tips not in the desired range
  subtree <- drop.tip(phy, setdiff(phy$tip.label, tips_to_keep))
  
  # Remove tip labels except for the species of interest
  subtree$tip.label <- ifelse(subtree$tip.label %in% tips_to_keep, subtree$tip.label, "")
  
  return(subtree)
}

find_closest_relatives <- function(phy, species_name, num_relatives = 200) {
  # Assuming species names are in the tip labels of the phylogeny
  tip_labels <- phy$tip.label
  
  # Find the index of the species name
  index <- grep(species_name, tip_labels)
  
  if (length(index) == 0) {
    stop("Species not found in the phylogeny.")
  }
  
  # Get the distances from the species to all other tips
  distances <- cophenetic(phy)[index, ]
  
  # Find the closest relatives
  closest_relatives <- order(distances)[1:num_relatives]
  
  return(tip_labels[closest_relatives])
}


get_evol_dist_statement <- function(evol_dist_df, species_name) {
  evol_dist <- evol_dist_df$w[species_name == evol_dist_df$Species]
  percentile <- quantile(evol_dist_df$w, 0.90)
  distinctiveness <- ifelse(evol_dist >= percentile, "Evolutionary Distinct", "Evolutionary Indistinct")
  return(distinctiveness)
}

get_evol_dist_dataframe <- function(evol_dist_df, percentile = 0.90) {
  evol_distinct <- evol_dist_df$w >= quantile(evol_dist_df$w, percentile)
  species_names <- evol_dist_df$Species
  evol_dist_df <- data.frame(Species = species_names, EvolutionaryDistinct = evol_distinct)
  return(evol_dist_df)
}

plot_phylogeny_with_color <- function(phy, evol_dist_df, species_name) {
  # Convert phy to ggtree object
  ggtree_obj <- ggtree(phy, layout = "circular")
  
  # Match tip labels and evolutionary distinctiveness
  tip_labels <- phy$tip.label
  matched_indices <- match(tip_labels, evol_dist_df$Species)
  
  # Check if any tip labels are not matched
  unmatched_labels <- tip_labels[is.na(matched_indices)]
  if (length(unmatched_labels) > 0) {
    cat("Unmatched tip labels:", paste(unmatched_labels, collapse = ", "), "\n")
  }
  
  # Add tip colors based on evolutionary distinctiveness
  ggtree_obj <- ggtree_obj %<+% evol_dist_df +
    geom_tippoint(aes(color = EvolutionaryDistinct), size = 3)

  
  # Plot the tree
 circular_phylogeny <- ggtree_obj + 
    geom_tiplab(aes(label = ifelse(label == species_name, gsub("_", " ", label), "")), size = 3, hjust = -0.5, vjust = 0.5) +
    theme(legend.position = c(1.6, 0.7),      # Adjust legend position to the right
          legend.justification = c(1, 0.5))  +
    labs(color = "Evolutionary Distinctiveness")

  
  # Print the circular phylogeny
  print(circular_phylogeny)
}



# Function to create range polygons
create_range_polygons <- function(data_frame, num_cores = 1, min_points = 5, min_distance = 1, gbif_limit = 2000) {
  rownames(data_frame) <- gsub("_", " ", rownames(data_frame))
  
  species_names <- data_frame$species_name
  
  range_results <- lapply(species_names, function(species_name) {
    message("Processing species:", species_name)
    
    # Remove leading numbers from species name using regular expression
    cleaned_species_name <- gsub("^\\d+\\.", "", species_name)
    
    gbif_data <- rgbif::occ_search(
      scientificName = cleaned_species_name,
      hasCoordinate = TRUE,
      limit = gbif_limit
    )
    
    if (length(gbif_data$data) == 0) {
      message(paste("No data found for", cleaned_species_name))
      return(NULL) # Return NULL if no data found
    }
    
    # Check if "occurrenceID" column is present
    if ("occurrenceID" %in% colnames(gbif_data$data)) {
      # Filter based on specified criteria
      gbif_data$data <- gbif_data$data %>%
        filter(
          !is.na(decimalLatitude),
          !is.na(decimalLongitude),
          decimalLatitude >= -90 & decimalLatitude <= 90,
          decimalLongitude >= -180 & decimalLongitude <= 180,
          decimalLatitude != decimalLongitude,
          !duplicated(gbif_data$data$occurrenceID),
          !is.na(countryCode)
        )
    } else {
      # Handle the case where "occurrenceID" column is not present
      message("The 'occurrenceID' column is not present in the data frame for ", cleaned_species_name)
      return(NULL) # Return NULL if no data found
    }
    
    if (nrow(gbif_data$data) == 0) {
      message(paste("No data left after filtering for", cleaned_species_name))
      return(NULL) # Return NULL if no data found
    }
    
    gbif_sf <- st_as_sf(gbif_data$data, coords = c("decimalLongitude", "decimalLatitude"))
    gbif_coords <- st_coordinates(gbif_sf)
    gbif_coords <- as.data.frame(gbif_coords)
    cluster_result <- dbscan(gbif_coords, eps = min_distance, minPts = min_points)
    gbif_coords$Cluster <- cluster_result$cluster
    sdf <- st_as_sf(gbif_coords, coords = c("X", "Y"), crs = 4326)
    
    sdf_filtered <- sdf[sdf$Cluster != 0, ]
    
    # Merge overlapping polygons within each cluster
    sdf_list <- lapply(split(sdf_filtered, sdf_filtered$Cluster), function(cluster_data) {
      merged_data <- st_union(cluster_data)
      if (!st_is_valid(merged_data)) {
        merged_data <- st_make_valid(merged_data)
      }
      merged_data
    })
    
    # Calculate convex hulls
    convex_hulls <- lapply(sdf_list, function(merged_data) {
      st_convex_hull(merged_data)
    })
    
    return(convex_hulls)
  })
  
  # Check if all species have no data, if so, return NULL
  if(all(sapply(range_results, is.null))) {
    return(NULL)
  }
  
  return(range_results)
}



# Get continent polygon boundaries
get_continent_sf <- function(url) {
  # Read the geojson file
  countries <- geojsonio::geojson_read(url, what = "sp")
  
  # Convert to sf
  countries_sf <- st_as_sf(countries)
  
  # Filter out countries with valid continent information
  continent_countries <- countries_sf[!is.na(countries_sf$continent), ]
  
  # Group by continent and combine geometries
  continent_polygons <- aggregate(continent_countries["geometry"], by=list(continent_countries$continent), FUN = function(x) st_union(x))
  
  # Create a new sf object with continent and geometry
  continent_sf <- st_sf(continent = continent_polygons$Group.1, geometry = continent_polygons$geometry)
  
  return(continent_sf)
}

# Function to clip range polygons to continent boundaries 
clip_polygons_to_land_shiny <- function(polygons_list, continent_sf) {
  # Check if the input is NULL
  if (is.null(polygons_list) || length(polygons_list) == 0) {
    return(NULL)
  }
  
  clipped_polygons_list <- lapply(polygons_list, function(convex_hulls) {
    clipped_polygons <- lapply(convex_hulls, function(ch) {
      tryCatch({
        # Set the CRS of the input polygon to match the CRS of the continent polygons
        st_crs(ch) <- st_crs(continent_sf)
        
        # Perform the clipping operation
        land_polygons <- ms_clip(ch, continent_sf)
        
        # Ensure validity and clean the geometry
        if (!is.null(land_polygons) && length(land_polygons) > 0) {
          valid_polygons <- st_make_valid(land_polygons)
          if (st_is_empty(valid_polygons)) {
            return(NULL)
          } else {
            return(valid_polygons)
          }
        } else {
          return(NULL)
        }
      }, error = function(e) {
        message(paste("Error clipping polygon:", e$message))
        return(NULL)
      })
    })
    # Filter out NULL results
    clipped_polygons <- Filter(Negate(is.null), clipped_polygons)
    return(clipped_polygons)
  })
  
  return(clipped_polygons_list)
}

# Function to calculate range size
calculate_range_size <- function(clipped_polygons_list) {
  if (is.null(clipped_polygons_list) || length(clipped_polygons_list) == 0) {
    return(0)
  }
  
  # Initialize a vector to store area sizes
  areas <- numeric()
  
  # Iterate over each clipped polygon
  for (clipped_polygons in clipped_polygons_list) {
    for (clipped_polygon in clipped_polygons) {
      tryCatch({
        # Check if the clipped polygon is valid
        if (st_is_valid(clipped_polygon)) {
          # Calculate the area and append it to the areas vector
          area <- st_area(clipped_polygon) / 1e6  # Convert to square kilometers
          areas <- c(areas, area)
        }
      }, error = function(e) {
        # Handle invalid geometries
        warning("Error calculating area for clipped polygon: ", conditionMessage(e))
      })
    }
  }
  
  # Sum the areas to get the total range size
  range_size <- sum(areas)
  return(range_size)
}



# New functions 

search_species <- function(dataframe, species_column) {
  # Assume dataframe has a column named 'species_name' containing species names
  
  # Get list of .csv files in the folder
  folder_path <- "Invasive_dataframes"
  csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  # Initialize a list to store results
  search_results <- vector("list", length = nrow(dataframe))
  
  # Loop through each species in the specified column
  for (i in seq_along(dataframe[[species_column]])) {
    species_name <- dataframe[[species_column]][i]
    
    # Initialize a vector to store country names where the species is found
    found_in_files <- character(0)
    
    # Loop through each .csv file
    for (file in csv_files) {
      # Read the file
      data <- read.csv(file, fill = TRUE)
      
      # Check if the "accepted_name.species" column exists
      if ("accepted_name.species" %in% colnames(data)) {
        # Check if the species_name exists in the column "accepted_name.species"
        if (any(data$accepted_name.species == species_name)) {
          # Append the country name (filename without extension) to found_in_files
          found_in_files <- c(found_in_files, tools::file_path_sans_ext(basename(file)))
        }
      }
    }
    
    # Store the search result for this species
    search_results[[i]] <- found_in_files
  }
  
  # Add a new column to the existing dataframe with the search results
  dataframe$found_in_files <- search_results
  
  return(dataframe)
}

# Function to get invasive countries
get_invasive_countries <- function(invasive_species_list) {
  # Extract unique invasive countries from the search results
  invasive_countries <- unique(unlist(invasive_species_list$found_in_files))
  
  # If "United States of America" is in the list, change it to "United States"
  if ("United States of America" %in% invasive_countries) {
    invasive_countries[invasive_countries == "United States of America"] <- "United States"
  }
  
  return(invasive_countries)
}

# Function to get country polygons for invasive countries
get_invasive_country_polygons <- function(countries, invasive_country_names) {
  # Filter countries data frame to get polygons of invasive countries
  invasive_country_polygons <- countries[countries$name %in% invasive_country_names, ]
  
  return(invasive_country_polygons)
}

# Function to plot convex hulls on the Leaflet map
plot_convex_hulls_map <- function(map, convex_hulls, invasive_country_sf) {
  if (!is.null(convex_hulls) && length(convex_hulls) > 0) {
    # Remove existing legend if present
    map <- map %>% clearControls()
    
    for (species_convex_hulls in convex_hulls) {
      if (!is.null(species_convex_hulls) && length(species_convex_hulls) > 0) {
        for (i in 1:length(species_convex_hulls)) {
          hull <- species_convex_hulls[[i]]
          
          if (!is.null(hull) && st_is_valid(hull) &&
              (st_geometry_type(hull) == "POLYGON" || st_geometry_type(hull) == "MULTIPOLYGON")) {
            # Ensure the hull is valid
            hull <- st_make_valid(hull)
            
            # Check if any point of the polygon falls within the invasive polygons
            intersects <- st_intersects(hull, invasive_country_sf, sparse = FALSE)
            
            if (any(intersects)) {  # If any point intersects, color red
              map <- map %>% addPolygons(data = hull, color = "red", fillColor = "red", group = "Invasive Range")
            } else {
              map <- map %>% addPolygons(data = hull, color = "blue", fillColor = "blue", group = "Native Range")
            }
          } else {
            message(paste("Skipping invalid or NULL polygon at index", i))
          }
        }
      } else {
        message("Skipping empty or NULL species_convex_hulls")
      }
    }
    
    # Add legend
    map <- map %>%
      addLegend(position = "bottomleft",
                colors = c("red", "blue"),
                labels = c("Invasive Range", "Native Range"))
  } else {
    # If no convex hulls to plot, just return the base map
    map <- map %>% addTiles() %>% setView(lng = 0, lat = 0, zoom = 2)
  }
  
  return(map)
}


source_and_clean_occurrences_decade <- function(species_name, gbif_limit_per_decade = 500) {
  # Clean species name
  cleaned_species_name <- gsub(" ", "%20", species_name)
  
  # Initialize a list to store decade-wise occurrence data
  decade_occurrences <- list()
  
  # Fetch occurrence data for each decade
  for (decade in 196:202) {
    # Generate year range for the decade
    year_start <- decade * 10
    year_end <- year_start + 9
    
    # Fetch occurrence data for the decade
    gbif_data <- rgbif::occ_search(
      scientificName = cleaned_species_name,
      limit = gbif_limit_per_decade,
      year = paste(year_start, ",", year_end)
    )
    
    # Check if data is not null
    if (!is.null(gbif_data$data) && nrow(gbif_data$data) > 0) {
      # Filter based on specified criteria
      cleaned_data <- gbif_data$data %>%
        filter(
          !is.na(decimalLatitude),
          !is.na(decimalLongitude),
          decimalLatitude >= -90 & decimalLatitude <= 90,
          decimalLongitude >= -180 & decimalLongitude <= 180,
          decimalLatitude != decimalLongitude,
          !duplicated(occurrenceID),
          !is.na(countryCode)
        )
      
      # Add the decade information
      cleaned_data$decade <- paste(decade, "s")
      
      # Store data for this decade as a dataframe
      decade_occurrences[[paste(decade)]] <- as.data.frame(cleaned_data)
    } else {
      message("No occurrence data found for decade", paste(decade, "s"))
    }
  }
  
  return(decade_occurrences)
}


plot_occurrences_with_leaflet <- function(data) {
  # Extract latitude and longitude
  coordinates <- data[, c("decimalLatitude", "decimalLongitude")]
  
  # Apply DBSCAN clustering
  dbscan_result <- dbscan(coordinates, eps = 0.1, minPts = 5)
  
  # Add cluster information to data
  data$cluster <- dbscan_result$cluster
  
  # Calculate the size of bubbles based on the number of points in each cluster
  cluster_sizes <- table(dbscan_result$cluster)
  max_cluster_size <- max(cluster_sizes)
  
  # Define the size of the largest bubble
  max_bubble_size <- 20  # Adjust as needed
  
  # Scale the bubble sizes to maintain the proportion relative to the largest cluster
  bubble_sizes <- sqrt(cluster_sizes) / sqrt(max_cluster_size) * max_bubble_size
  
  # Create Leaflet map
  map <- leaflet() %>%
    addTiles() %>%
    addCircleMarkers(
      data = data,
      lng = ~decimalLongitude,
      lat = ~decimalLatitude,
      color = "red",  # Color can be customized
      radius = bubble_sizes[dbscan_result$cluster],
      opacity = 0.8,
      popup = paste("Cluster:", data$cluster),
      group = "Occurrences"
    ) %>%
    addLayersControl(
      overlayGroups = "Occurrences",
      options = layersControlOptions(collapsed = FALSE)
    )
  
  return(map)
}


plot_occurrences_by_decade <- function(occurrences_by_decade) {
  leaflet_maps <- list()
  for (i in seq_along(occurrences_by_decade)) {
    data <- occurrences_by_decade[[i]]
    leaflet_map <- plot_occurrences_with_leaflet(data)
    leaflet_map <- div(leaflet_map, class = paste0("map-div map-", i))
    leaflet_maps[[i]] <- leaflet_map
  }
  return(leaflet_maps)
}


# Define the plot_occurrences_pop function
plot_occurrences_pop <- function(decade_occurrences, epsilon = 0.075, minPts = 2) {
  # Generate colors based on decades
  unique_decades <- unique(names(decade_occurrences))
  num_decades <- length(unique_decades)
  colors <- rainbow(num_decades)
  
  # Get world boundaries
  world <- map("world", fill = TRUE, plot = FALSE)
  
  # Create Leaflet map
  pop_map <- leaflet() %>%
    addTiles(urlTemplate = "https://cartodb-basemaps-{s}.global.ssl.fastly.net/dark_all/{z}/{x}/{y}.png") %>%
    addPolygons(data = world, 
                fillColor = "transparent", 
                color = "white", 
                weight = 1,
                opacity = 1,
                fillOpacity = 0) %>%
    setView(lng = 0, lat = 0, zoom = 2) # Set initial view
  
  # Combine all decades' data into one dataframe
  all_data <- bind_rows(decade_occurrences, .id = "decade")
  
  # Calculate clusters using DBSCAN on combined data
  dbscan_result <- dbscan(all_data[, c("decimalLongitude", "decimalLatitude")], eps = epsilon, minPts = minPts)
  
  if (max(dbscan_result$cluster) > 1) {
    # Connect points within each cluster
    for (cluster_id in unique(dbscan_result$cluster[dbscan_result$cluster > 0])) {
      cluster_points <- all_data[dbscan_result$cluster == cluster_id, ]
      
      cluster_points <- as.data.frame(cluster_points) # Convert to data frame to ensure columns are numeric
      cluster_points$decimalLatitude <- as.numeric(cluster_points$decimalLatitude)
      cluster_points$decimalLongitude <- as.numeric(cluster_points$decimalLongitude)
      
      # Sort points by decade
      cluster_points <- cluster_points[order(as.numeric(cluster_points$decade)), ]
      
      # Find the color of the first point connected by each line
      first_point_decade <- unique(cluster_points$decade)[1]
      color <- colors[match(first_point_decade, unique_decades)]
      
      for (i in 1:(nrow(cluster_points) - 1)) {
        point1 <- c(cluster_points[i, "decimalLatitude"], cluster_points[i, "decimalLongitude"])
        point2 <- c(cluster_points[i + 1, "decimalLatitude"], cluster_points[i + 1, "decimalLongitude"])
        
        # Find the earliest decade of the connected points
        connected_points_decades <- unique(c(cluster_points$decade[i], cluster_points$decade[i + 1]))
        earliest_decade <- min(connected_points_decades)
        
        color <- colors[match(earliest_decade, unique_decades)]
        
        pop_map <- addPolylines(pop_map, lng = c(point1[2], point2[2]), lat = c(point1[1], point2[1]), color = color, weight = 2)
      }
    }
  } else {
    message("No clusters found.")
  }
  
  # Loop through each decade's data
  for (i in 1:num_decades) {
    decade_data <- decade_occurrences[[unique_decades[i]]]
    
    if (nrow(decade_data) > 0) {
      # Add circle markers for each occurrence point
      pop_map <- addCircleMarkers(
        pop_map,
        data = decade_data,
        lng = ~decimalLongitude,
        lat = ~decimalLatitude,
        radius = 4, # Constant point size
        color = colors[i],
        stroke = FALSE,
        fillOpacity = 0.7,
        popup = ~paste("Decade:", unique_decades[i])
      )
    } else {
      message("No occurrence data found for decade", unique_decades[i])
    }
  }
  
  # Create color legend
  legend_colors <- colors
  legend_labels <- unique_decades
  legend_title <- "Decade"
  
  pop_map <- addLegend(
    pop_map,
    position = "bottomleft",
    colors = legend_colors,
    labels = paste(legend_labels, "0", sep = ""),
    #labels = legend_labels,
    title = legend_title
  )
  
  return(pop_map)
}



# Define custom sidebarPanel function
sidebarPanel2 <- function (..., out = NULL, width = 4) {
  div(class = paste0("col-sm-", width), 
      tags$form(class = "well", ...),
      out
  )
}


calculate_average_metric <- function(cluster_points) {
  total_distance_change <- 0
  total_year_change <- 0
  valid_pairs <- 0
  
  # Loop through each pair of consecutive points in the cluster
  for (i in 1:(nrow(cluster_points) - 1)) {
    # Extract latitude, longitude, and year for two consecutive points
    lat1 <- cluster_points[i, "decimalLatitude"]
    lon1 <- cluster_points[i, "decimalLongitude"]
    lat2 <- cluster_points[i + 1, "decimalLatitude"]
    lon2 <- cluster_points[i + 1, "decimalLongitude"]
    year1 <- as.numeric(cluster_points[i, "year"])
    year2 <- as.numeric(cluster_points[i + 1, "year"])
    
    # Calculate geographic distance between the two points
    distance <- geosphere::distHaversine(c(lon1, lat1), c(lon2, lat2))
    
    # Calculate the difference in years
    year_diff <- year2 - year1
    
    # Check if the years are different (to avoid division by zero)
    if (year_diff != 0) {
      # Accumulate total distance change and year change
      total_distance_change <- total_distance_change + distance
      total_year_change <- total_year_change + year_diff
      valid_pairs <- valid_pairs + 1
    }
  }
  
  # Check if there are valid pairs
  if (valid_pairs > 0) {
    # Calculate the average metric by dividing total distance change by total year change
    return(total_distance_change / total_year_change)
  } else {
    # If no valid pairs found, return NA
    return(NA)
  }
}



calculate_cluster_metrics <- function(decade_occurrences, epsilon = 0.075, minPts = 2) {
  # Combine all decades' data into one dataframe
  all_data <- bind_rows(decade_occurrences, .id = "decade")
  
  # Check for latitude values exceeding 90
  if (any(all_data$decimalLatitude > 90)) {
    cat("Problematic latitude values:\n")
    print(all_data$decimalLatitude[all_data$decimalLatitude > 90])
    stop("Latitude values exceed 90.")
  }
  
  # Calculate clusters using DBSCAN on combined data
  dbscan_result <- dbscan(all_data[, c("decimalLongitude", "decimalLatitude")], eps = epsilon, minPts = minPts)
  
  cluster_metrics <- list()
  
  if (max(dbscan_result$cluster) > 1) {
    # Calculate metrics for each cluster
    cluster_distances <- c()
    for (cluster_id in unique(dbscan_result$cluster[dbscan_result$cluster > 0])) {
      cluster_points <- all_data[dbscan_result$cluster == cluster_id, ]
      
      # Sort points by year
      cluster_points <- cluster_points[order(as.numeric(cluster_points$year)), ]
      
      # Calculate average distance for cluster
      average_distance <- calculate_average_metric(cluster_points)
      cluster_distances <- c(cluster_distances, average_distance)
    }
    
    # Calculate average distance across clusters
    avg_distance_across_clusters <- mean(cluster_distances, na.rm = TRUE)
    cluster_metrics <- avg_distance_across_clusters
  } else {
    message("No clusters found.")
    cluster_metrics <- NA
  }
  
  return(cluster_metrics)
}

# Function for calculating native versus invasive range sizes 
calculate_range_size_separated <- function(clipped_polygons_list, invasive_country_sf) {
  if (is.null(clipped_polygons_list) || length(clipped_polygons_list) == 0) {
    return(list(invasive_range_size = 0, native_range_size = 0))
  }
  
  # Convert clipped_polygons_list to sf objects
  clipped_poly_sf <- lapply(clipped_polygons_list[[1]], function(poly) {
    if (inherits(poly, "sfc_MULTIPOLYGON")) {
      st_sf(st_sfc(poly))
    } else if (inherits(poly, "sfc_POLYGON")) {
      st_sf(st_sfc(poly))
    } else {
      stop("Unsupported geometry type")
    }
  })
  
  # Convert invasive_country_sf to sf object if not already
  inv_country_sf <- st_as_sf(invasive_country_sf)
  
  # Check for invalid polygons in invasive_country_sf
  invalid_indices <- which(!st_is_valid(inv_country_sf))
  if (length(invalid_indices) > 0) {
    print("Warning: Invalid polygons detected in invasive_country_sf.")
    print(paste("Invalid indices:", invalid_indices))
    
    # Attempt to make polygons valid
    inv_country_sf <- st_make_valid(inv_country_sf)
    
    # Check again for invalid polygons
    invalid_indices <- which(!st_is_valid(inv_country_sf))
    if (length(invalid_indices) > 0) {
      stop("Could not make all polygons in invasive_country_sf valid.")
    } else {
      print("All polygons in invasive_country_sf are now valid.")
    }
  } else {
    print("All polygons in invasive_country_sf are valid.")
  }
  
  # Initialize vectors to store area sizes
  invasive_areas <- numeric()
  native_areas <- numeric()
  
  # Iterate over each clipped polygon
  for (clipped_polygons in clipped_polygons_list) {
    for (clipped_polygon in clipped_polygons) {
      tryCatch({
        # Check if the clipped polygon is valid
        if (st_is_valid(clipped_polygon)) {
          # Check intersection with invasive country polygons
          intersects <- st_intersects(clipped_polygon, inv_country_sf, sparse = FALSE)
          
          if (any(intersects)) {
            # Calculate area and add to invasive areas
            area <- st_area(clipped_polygon) / 1e6  # Convert to square kilometers
            invasive_areas <- c(invasive_areas, area)
          } else {
            # Calculate area and add to native areas
            area <- st_area(clipped_polygon) / 1e6  # Convert to square kilometers
            native_areas <- c(native_areas, area)
          }
        }
      }, error = function(e) {
        # Handle invalid geometries
        warning("Error calculating area for clipped polygon: ", conditionMessage(e))
      })
    }
  }
  
  # Sum the areas to get the total range sizes
  invasive_range_size <- sum(invasive_areas)
  native_range_size <- sum(native_areas)
  
  # Return as a list
  return(list(invasive_range_size = invasive_range_size, native_range_size = native_range_size))
}

calculate_species_bbox <- function(occurrences) {
  if (is.null(occurrences)) return(NULL)
  
  # If occurrences is a list of species occurrences
  if (inherits(occurrences, "list")) {
    # Extract the first species' bounding box
    first_species <- occurrences[[1]]
  } else {
    # Assume occurrences is a single species' occurrence
    first_species <- occurrences
  }
  
  # Calculate bounding box
  bbox <- sf::st_bbox(first_species)
  
  # Create a list with the extracted values
  species_bbox <- list(
    xmin = bbox$xmin,
    ymin = bbox$ymin,
    xmax = bbox$xmax,
    ymax = bbox$ymax
  )
  
  # Return the list
  return(species_bbox)
}

download_and_crop_worldclim <- function(species_bbox, vars, output_dir) {
  # Extract numeric values from species_bbox list
  xmin <- species_bbox[["xmin"]]
  ymin <- species_bbox[["ymin"]]
  xmax <- species_bbox[["xmax"]]
  ymax <- species_bbox[["ymax"]]
  
  # Download WorldClim data for the specified variables and crop to the extent
  rasters <- worldclim_global(var = "bio", lon = c(xmin, xmax), lat = c(ymin, ymax), path = output_dir, version = "2.1", res = 10)
  
  # Crop the downloaded rasters to the extent
  extent_bbox <- extent(xmin, xmax, ymin, ymax)
  rasters_cropped <- lapply(rasters, function(raster) {
    crop(raster, extent_bbox)
  })
  
  return(rasters_cropped)
}

output_directory <- tempdir()


convert_to_RasterLayerList <- function(rasters_list) {
  raster_layers <- list()
  
  for (i in seq_along(rasters_list)) {
    raster_obj <- raster(rasters_list[[i]])
    raster_layers[[i]] <- raster_obj
  }
  
  return(raster_layers)
}

get_occurrences <- function(data_frame, num_cores = 1, min_points = 5, min_distance = 1, gbif_limit = 2000) {
  rownames(data_frame) <- gsub("_", " ", rownames(data_frame))
  
  species_names <- data_frame$species_name
  
  range_results <- lapply(species_names, function(species_name) {
    message("Processing species:", species_name)
    
    # Remove leading numbers from species name using regular expression
    cleaned_species_name <- gsub("^\\d+\\.", "", species_name)
    
    gbif_data <- rgbif::occ_search(
      scientificName = cleaned_species_name,
      hasCoordinate = TRUE,
      limit = gbif_limit
    )
    
    if (length(gbif_data$data) == 0) {
      message(paste("No data found for", cleaned_species_name))
      return(NULL) # Return NULL if no data found
    }
    
    # Check if "occurrenceID" column is present
    if ("occurrenceID" %in% colnames(gbif_data$data)) {
      # Filter based on specified criteria
      gbif_data$data <- gbif_data$data %>%
        dplyr::filter(
          !is.na(decimalLatitude),
          !is.na(decimalLongitude),
          decimalLatitude >= -90 & decimalLatitude <= 90,
          decimalLongitude >= -180 & decimalLongitude <= 180,
          decimalLatitude != decimalLongitude,
          !duplicated(occurrenceID),
          !is.na(countryCode)
        )
    } else {
      # Handle the case where "occurrenceID" column is not present
      message("The 'occurrenceID' column is not present in the data frame for ", cleaned_species_name)
      return(NULL) # Return NULL if no data found
    }
    
    if (nrow(gbif_data$data) == 0) {
      message(paste("No data left after filtering for", cleaned_species_name))
      return(NULL) # Return NULL if no data found
    }
    
    gbif_sf <- sf::st_as_sf(gbif_data$data, coords = c("decimalLongitude", "decimalLatitude"))
    gbif_coords <- sf::st_coordinates(gbif_sf)
    gbif_coords <- as.data.frame(gbif_coords)
    cluster_result <- dbscan::dbscan(gbif_coords, eps = min_distance, minPts = min_points)
    gbif_coords$Cluster <- cluster_result$cluster
    sdf <- sf::st_as_sf(gbif_coords, coords = c("X", "Y"), crs = 4326)
    
    sdf_filtered <- sdf[sdf$Cluster != 0, ]
    
    # Return cleaned coordinates dataframe
    return(sdf_filtered)
  })
  
  # Check if all species have no data, if so, return NULL
  if(all(sapply(range_results, is.null))) {
    return(NULL)
  }
  
  # Set names for the list elements based on species names
  range_results <- setNames(range_results, species_names)
  
  return(range_results)
}

calculate_avg_var <- function(occurrences_list, vars) {
  
  # Initialize a list to store results for each species
  results_list <- list()
  
  # Iterate over each species occurrence dataframe
  for (species_name in names(occurrences_list)) {
    message("Processing species ", species_name)
    
    # Get species occurrence coordinates
    species_occ <- occurrences_list[[species_name]]
    
    # Initialize a list to store average max temperatures for each raster
    avg_var <- numeric(length = length(vars))
    
    # Iterate over each raster
    for (j in seq_along(vars)) {
      raster <- vars[[j]]
      
      # Extract raster values at occurrence points
      raster_values <- extract(raster, species_occ)
      
      # Compute mean of non-NA values
      avg_var[j] <- mean(raster_values, na.rm = TRUE)
    }
    
    # Store results for the current species
    results_list[[species_name]] <- avg_var
  }
  
  return(results_list)
}

avg_var_to_dataframe_bio_19 <- function(avg_var_results) {
  # Initialize an empty dataframe
  df <- data.frame(
    species_name = character(),
    stringsAsFactors = FALSE
  )
  
  # Define the specific indices or variables in the order you want
  variable_names <- c(
    "Annual mean temperature", "Max temperature of warmest month",
    "Min temperature of coldest month", "Temperature Annual Range",
    "Annual precipitation", "Precipitation of wettest month",
    "Precipitation of driest month"
  )
  
  # Iterate over each species in avg_var_results
  for (species_name in names(avg_var_results)) {
    # Extract average variable values for the current species
    avg_vars <- avg_var_results[[species_name]]
    
    # Create a new row for the species in the dataframe
    species_row <- data.frame(
      species_name = species_name,
      t(avg_vars)  # Transpose to convert to row format
    )
    
    # Append the row to the dataframe
    df <- rbind(df, species_row)
  }
  
  # Rename the columns of the dataframe to represent variable names
  colnames(df)[-1] <- variable_names
  
  return(df)
}


ui <- fluidPage(
  titlePanel("EcoSMART"),
  tags$div(
    style = "font-size: 16px; text-align: left; color: #666; margin-bottom: 14px;",
    "Statistics, Modeling, Analysis of Ranges and Traits"
  ),
  sidebarLayout(
    sidebarPanel(
      textInput("species_name", "Enter Species Name"),
      actionButton("submit", "Submit"),
      div(id = "progress", style = "position: fixed; bottom: 0; width: 100%;")
    ),
    mainPanel(
      fluidRow(
        column(12,
          h4(textOutput("species_label"), style = "color: #333"),
          br(),
          tabsetPanel(
            tabPanel("Range Clusters", 
                     leafletOutput("map"),
                     htmlOutput("range_size_html")),
            tabPanel("Occurrences by Decade",
                     sliderTextInput("decade", "Decade", choices = seq(1960, 2020, by = 10), width = "100%", animate = TRUE),
                     uiOutput("map_output")),
            tabPanel("Occurrences with Clusters",
                     leafletOutput("cluster_map"),
                     textOutput("avg_dist")),
            tabPanel("Summary Statistics",
                     htmlOutput("summary_range_size_html"),
                     textOutput("summary_avg_dist"),
                     textOutput("avg_temp"),
                     textOutput("temp_range"),
                     textOutput("avg_precipitation")
            ),
            tabPanel("Climate Data",
                     tableOutput("bio_var_table")
            )
          )
        )
      )
    )
  )
)


server <- function(input, output, session) {
  
  species_name <- reactiveVal()
  convex_hulls_list <- reactiveVal()
  range_size <- reactiveVal(0)
  invasive_range_size <- reactiveVal(0)
  native_range_size <- reactiveVal(0)
  range_calc_in_progress <- reactiveVal(FALSE)
  invasive_country_sf <- reactiveVal()
  occurrences_by_decade <- reactiveVal(NULL)
  maps_by_decade <- reactiveVal(NULL)
  avg_dist <- reactiveVal("")
  bio_vars <- reactiveVal(NULL)
  avg_temp <- reactiveVal(0)
  temp_range <- reactiveVal(0)
  avg_precipitation <- reactiveVal(0)
  
  # Create initial leaflet map
  output$map <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%
      setView(lng = 0, lat = 0, zoom = 2)
  })
  
  observeEvent(input$submit, {
    withProgress(
      message = "Processing data...",
      detail = "This may take a while...",
      value = 0,
      {
        species_name(input$species_name)
        
        # Update progress bar value
        incProgress(0.1, detail = "Fetching data...")
        
        # Preprocess species name to replace underscores with spaces
        species_name_processed <- gsub("_", " ", input$species_name)
        
        data <- data.frame(species_name = c(species_name_processed))

        
        folder_path <- "~/Desktop/Invasive_dataframes"
        
        # Read the GeoJSON file containing country polygons if not already read
        countries <- geojsonio::geojson_read("https://d2ad6b4ur7yvpq.cloudfront.net/naturalearth-3.3.0/ne_50m_admin_0_countries.geojson", what = "sp")
        
        # Call the function to search for species in invasive countries
        invasive_species_list <- search_species(data, "species_name")

        # Update progress bar value
        incProgress(0.2, detail = "Mapping ranges...")

        continent_bounds <- get_continent_sf("https://d2ad6b4ur7yvpq.cloudfront.net/naturalearth-3.3.0/ne_50m_admin_0_countries.geojson")

        # Create range polygons for species
        convex_hulls_old <- create_range_polygons(invasive_species_list)

        # Clip polygons to land
        convex_hulls <- clip_polygons_to_land_shiny(convex_hulls_old, continent_bounds)
        
       
        # Calculate range size
        if (!is.null(convex_hulls)) {
          range_size(round(calculate_range_size(convex_hulls), 3))
        } else {
          range_size(0)
        }
        
        # Get invasive countries
        invasive_countries <- get_invasive_countries(invasive_species_list)
        
        # Get country polygons for invasive countries
        invasive_country_polygons <- get_invasive_country_polygons(countries, invasive_countries)
        
        # Convert to sf object
        invasive_country_sf_data <- st_as_sf(invasive_country_polygons)
        invasive_country_sf_data <- st_make_valid(invasive_country_sf_data)
        invasive_country_sf(invasive_country_sf_data)
        
        # Calculate separated range sizes
        if (!is.null(convex_hulls) && !is.null(invasive_country_sf_data)) {
          range_sizes <- calculate_range_size_separated(convex_hulls, invasive_country_sf_data)
          invasive_range_size(round(range_sizes$invasive_range_size, 3))
          native_range_size(round(range_sizes$native_range_size, 3))
        } else {
          invasive_range_size(0)
          native_range_size(0)
        }
        
        # Plot convex hulls on the map
        if (!is.null(convex_hulls)) {
          leafletProxy("map") %>%
            clearShapes() %>%
            plot_convex_hulls_map(convex_hulls, invasive_country_sf_data)
        } else {
          leafletProxy("map") %>%
            clearShapes() %>%
            addTiles() %>%
            setView(lng = 0, lat = 0, zoom = 2)
        }
        


        bio_occurrences <- get_occurrences(data, num_cores = 1, min_points = 5, min_distance = 1, gbif_limit = 2000)

        bbox <- calculate_species_bbox(bio_occurrences)

        output_directory <- tempdir()

        rasters <- download_and_crop_worldclim(bbox, vars, output_directory)

        raster_list <- convert_to_RasterLayerList(rasters)

        selected_rasters <- list(raster_list[[1]], raster_list[[5]], raster_list[[6]], raster_list[[7]], raster_list[[12]], raster_list[[13]], raster_list[[14]])
 
        avg_bio_var_list <- calculate_avg_var(bio_occurrences, selected_rasters)
        avg_bio_var_df <- avg_var_to_dataframe_bio_19(avg_bio_var_list)
        

        bio_vars(avg_bio_var_df)

        avg_temp(avg_bio_var_df$`Annual mean temperature`)
        temp_range(avg_bio_var_df$`Temperature Annual Range`)
        avg_precipitation(avg_bio_var_df$`Annual precipitation`)

        # Update progress bar value
        incProgress(0.3, detail = "Fetching occurrences by decade...")
        
        # Fetch occurrences by decade
        occurrences_by_decade_val <- source_and_clean_occurrences_decade(input$species_name)

        # Update progress bar value
        incProgress(0.6, detail = "Rendering maps...")
        
        # Generate maps for each decade if convex_hulls is not NULL
        if (!is.null(convex_hulls)) {
          maps_by_decade_val <- plot_occurrences_by_decade(occurrences_by_decade_val)
          
          # Update reactive values
          occurrences_by_decade(occurrences_by_decade_val)
          maps_by_decade(maps_by_decade_val)
        } else {
          occurrences_by_decade(NULL)
          maps_by_decade(NULL)
        }
        
        # Reset slider to 1960
        updateSliderTextInput(session, "decade", selected = 1960)
        
        # Update progress bar value
        incProgress(1, detail = "Processing complete...")
        
      }
    )
  })
  
  # Calculate average distance if occurrences are available
  observe({
    if (!is.null(occurrences_by_decade())) {
      avg_dist_statement <- calculate_cluster_metrics(occurrences_by_decade())
      avg_dist(avg_dist_statement)
    } else {
      avg_dist(NULL)
    }
  })
  

  output$avg_dist <- renderText({
    if (!is.null(avg_dist())) {
     avg_distance <- avg_dist()
      if (!is.null(avg_distance) && !is.na(avg_distance)) {
       avg_distance_m_per_year <- round(avg_distance, 2)  # Round to two decimal places
       paste("Average dispersal is", avg_distance_m_per_year, "m per year")
     } else {
        ""
     }
   } else {
    "Average dispersal is 0m per year"
   }
  })

  # Display range size

    output$range_size_html <- renderUI({
    req(range_size(), native_range_size(), invasive_range_size())
    HTML(paste("Total Range Size (km^2): ", range_size(), "<br>",
               "&emsp;- Native Range Size (km^2): ", native_range_size(), "<br>",
               "&emsp;- Invasive Range Size (km^2): ", invasive_range_size()))
  })

  
  # Reactive value for selected decade
  selected_decade <- reactiveVal(1960)
  
  # Update decade slider input
  observe({
    updateSliderTextInput(
      session = session,
      inputId = "decade",
      choices = seq(1960, 2020, by = 10),
      selected = selected_decade()
    )
  })
  
  # Display map for the selected decade
  output$map_output <- renderUI({
    req(maps_by_decade())
    selected_map <- maps_by_decade()[[(selected_decade() - 1960) / 10 + 1]]  # Adjust indexing to match slider values
    div(id = "map-container", selected_map)
  })
  
  # Update selected decade when slider changes
  observeEvent(input$decade, {
    selected_decade(input$decade)
  })
  
  # Plot occurrences with clusters on the map
  output$cluster_map <- renderLeaflet({
    req(occurrences_by_decade(), invasive_country_sf())
    
    decade_occurrences <- occurrences_by_decade()
    if (!is.null(decade_occurrences)) {
      plot_occurrences_pop(decade_occurrences)
    } else {
      leaflet() %>%
        addTiles() %>%
        setView(lng = 0, lat = 0, zoom = 2)
    }
  })
  
  # Plot occurrences by decade
  output$occurrences_by_decade_map <- renderLeaflet({
    req(occurrences_by_decade(), invasive_country_sf())
    
    decade_occurrences <- occurrences_by_decade()
    if (!is.null(decade_occurrences)) {
      plot_occurrences_by_decade_map(decade_occurrences)
    } else {
      leaflet() %>%
        addTiles() %>%
        setView(lng = 0, lat = 0, zoom = 2)
    }
  })
  
  output$species_label <- renderText({
    req(input$submit)
    paste("Species Name: ", input$species_name)
  })
  
  # Summary Statistics Outputs
  output$summary_range_size_html <- renderUI({
    req(range_size(), native_range_size(), invasive_range_size())
    HTML(paste("Total Range Size (km^2): ", range_size(), "<br>",
               "&emsp;- Native Range Size (km^2): ", native_range_size(), "<br>",
               "&emsp;- Invasive Range Size (km^2): ", invasive_range_size()))
  })

  output$bio_var_table <- renderTable({
    req(bio_vars())
    bio_vars()
  })

  output$avg_temp <- renderText({
  req(avg_temp())
  paste("Average Temperature (Cº): ", round(avg_temp(), 2))
})

output$temp_range <- renderText({
  req(temp_range())
  paste("Temperature Range (Cº): ", round(temp_range(), 2))
})

output$avg_precipitation <- renderText({
  req(avg_precipitation())
  paste("Average Precipitation (mm): ", round(avg_precipitation(), 2))
})
  
  output$summary_avg_dist <- renderText({
    req(avg_dist())
    avg_distance <- avg_dist()
    if (!is.null(avg_distance) && !is.na(avg_distance)) {
      avg_distance_m_per_year <- round(avg_distance, 2)
      paste("Average Dispersal Rate (m/year):", avg_distance_m_per_year)
    } else {
      "Average Dispersal Rate (m/year): 0"
    }
  })
}



# Run the application
shinyApp(ui = ui, server = server)


