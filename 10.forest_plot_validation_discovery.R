# Forest plot comparing discovery and validation cohort effects
# For 15 proteins significant across obs+MR+hyprcoloc in all three groups
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Set paths
base_path <- "D:/OneDrive/工作/1.工作/blood pressure"

# Define models to process
models <- c("model_3")

# Loop through each model
for (current_model in models) {
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("Processing:", current_model, "\n")
  cat(rep("=", 80), "\n\n", sep = "")
  
  discovery_path <- file.path(base_path, current_model)
  validation_path <- file.path(base_path, paste0("sg_", current_model))
  
  # Check if paths exist
  if (!dir.exists(discovery_path)) {
    cat("Warning: Discovery path does not exist:", discovery_path, "\n")
    next
  }
  if (!dir.exists(validation_path)) {
    cat("Warning: Validation path does not exist:", validation_path, "\n")
    next
  }
  
  # Read the significant proteins from the heatmap CSV file
  protein_file <- file.path(discovery_path, "heatmap","Level3_Regression_MR_Coloc", "circular_heatmap_data.csv")
  
  if (!file.exists(protein_file)) {
    cat("Warning: Protein file does not exist:", protein_file, "\n")
    next
  }
  
  sig_proteins <- read.csv(protein_file, stringsAsFactors = FALSE)

  
  # Get unique protein names
  unique_proteins <- unique(sig_proteins$Entrez.Gene.Name)
  cat("Number of unique proteins:", length(unique_proteins), "\n")
  cat("Proteins:", paste(unique_proteins, collapse = ", "), "\n\n")
  
  # Read discovery cohort results
  discovery_obs_file <- file.path(discovery_path, "obs_all_results.txt")
  validation_obs_file <- file.path(validation_path, "obs_all_results.txt")
  
  if (!file.exists(discovery_obs_file)) {
    cat("Warning: Discovery obs file does not exist:", discovery_obs_file, "\n")
    next
  }
  if (!file.exists(validation_obs_file)) {
    cat("Warning: Validation obs file does not exist:", validation_obs_file, "\n")
    next
  }
  
  discovery_obs <- read.table(
    discovery_obs_file,
    header = TRUE, sep = "\t", stringsAsFactors = FALSE
  )
  
  # Read validation cohort results
  validation_obs <- read.table(
    validation_obs_file,
    header = TRUE, sep = "\t", stringsAsFactors = FALSE
  )

  
  # Filter for the significant proteins
  discovery_filtered <- discovery_obs %>%
    filter(Entrez.Gene.Name %in% unique_proteins)
  
  validation_filtered <- validation_obs %>%
    filter(Entrez.Gene.Name %in% unique_proteins)
  
  # Calculate confidence intervals (95% CI)
  discovery_filtered <- discovery_filtered %>%
    mutate(
      lower_ci = Estimate - 1.96 * se,
      upper_ci = Estimate + 1.96 * se,
      cohort = "Discovery"
    )
  
  validation_filtered <- validation_filtered %>%
    mutate(
      lower_ci = Estimate - 1.96 * se,
      upper_ci = Estimate + 1.96 * se,
      cohort = "Validation"
    )
  
  # Combine datasets
  combined_data <- bind_rows(discovery_filtered, validation_filtered)
  
  # Remove duplicates: keep only one result per protein-trait-cohort combination
  combined_data <- combined_data %>%
    distinct(Entrez.Gene.Name, traits, cohort, .keep_all = TRUE)
  
  # Truncate estimates and CIs to -1 to 1 range
  # Keep track of which values were truncated
  combined_data <- combined_data %>%
    mutate(
      # Mark if truncation occurred
      truncated_lower = lower_ci < -1,
      truncated_upper = upper_ci > 1,
      truncated_estimate = Estimate < -1 | Estimate > 1,
      
      # Apply truncation
      Estimate_original = Estimate,
      lower_ci_original = lower_ci,
      upper_ci_original = upper_ci,
      
      Estimate = pmin(pmax(Estimate, -1), 1),
      lower_ci = pmin(pmax(lower_ci, -1), 1),
      upper_ci = pmin(pmax(upper_ci, -1), 1)
    )
  
  # Create a combined label for protein and trait
  combined_data <- combined_data %>%
    mutate(
      protein_trait = paste(Entrez.Gene.Name, traits, sep = " - "),
      cohort = factor(cohort, levels = c("Discovery", "Validation"))
    )
  
  # Print summary
  cat("Discovery cohort entries:", nrow(discovery_filtered), "\n")
  cat("Validation cohort entries:", nrow(validation_filtered), "\n")
  cat("Total entries for plotting:", nrow(combined_data), "\n\n")
  
  # Define trait order based on user specification
  trait_order <- c("sbp", "dbp","pp","map","hypertension")

  # Define trait labels for display
  trait_labels <- c(
    "sbp" = "SBP",
    "dbp" = "DBP",
    "pp" = "PP",
    "map" = "MAP",
    "hypertension" = "Hypertension"
  )
  
  # Set traits as factor with specified order
  combined_data$traits <- factor(combined_data$traits, levels = trait_order)
  
  # Create output directory for validation plots
  validation_dir <- file.path(discovery_path, "validation")
  if (!dir.exists(validation_dir)) {
    dir.create(validation_dir, recursive = TRUE)
    cat("Created validation directory:", validation_dir, "\n")
  }

  
  # # Loop through each protein and create individual forest plot
  # for (protein in unique_proteins) {
  #   cat("\nCreating forest plot for protein:", protein, "\n")
  #   
  #   # Filter data for current protein
  #   protein_data <- combined_data %>%
  #     filter(Entrez.Gene.Name == protein) %>%
  #     arrange(traits)
  #   
  #   if (nrow(protein_data) == 0) {
  #     cat("  No data for protein:", protein, "\n")
  #     next
  #   }
  #   
  #   cat("  Traits found:", paste(unique(protein_data$traits), collapse = ", "), "\n")
  #   
  #   # Create forest plot for this protein
  #   p <- ggplot(protein_data, aes(x = Estimate, y = traits, color = cohort)) +
  #     geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
  #     geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), 
  #                    height = 0.3, 
  #                    position = position_dodge(width = 0.5),
  #                    size = 0.8) +
  #     geom_point(position = position_dodge(width = 0.5), 
  #                size = 3, 
  #                shape = 18) +
  #     # Add triangle markers for truncated values
  #     geom_point(data = protein_data %>% filter(truncated_lower),
  #                aes(x = -1, y = traits),
  #                shape = 17, size = 2, alpha = 0.7,
  #                position = position_dodge(width = 0.5)) +
  #     geom_point(data = protein_data %>% filter(truncated_upper),
  #                aes(x = 1, y = traits),
  #                shape = 17, size = 2, alpha = 0.7,
  #                position = position_dodge(width = 0.5)) +
  #     scale_color_manual(
  #       values = c("Discovery" = "#337CB2", "Validation" = "#D64643"),
  #       name = "Cohort"
  #     ) +
  #     scale_y_discrete(limits = rev(trait_order), labels = trait_labels[rev(trait_order)], drop = FALSE) +
  #     scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  #     labs(
  #       x = "Effect Size (95% CI) [Truncated at ±1]",
  #       y = "Trait",
  #       title = paste0("Forest Plot: ", protein, " (", current_model, ")"),
  #       subtitle = "Discovery vs Validation Cohort Effects"
  #     ) +
  #     theme_bw() +
  #     theme(
  #       axis.text.y = element_text(size = 10, face = "bold"),
  #       axis.text.x = element_text(size = 10),
  #       axis.title = element_text(size = 11, face = "bold"),
  #       legend.position = "top",
  #       legend.title = element_text(face = "bold"),
  #       plot.title = element_text(size = 13, face = "bold"),
  #       plot.subtitle = element_text(size = 10),
  #       panel.grid.major.y = element_line(color = "gray90"),
  #       panel.grid.minor = element_blank()
  #     )
  #   
  #   # Save plot for this protein
  #   output_file_pdf <- file.path(validation_dir, paste0("forest_plot_", protein, ".pdf"))
  #   output_file_png <- file.path(validation_dir, paste0("forest_plot_", protein, ".png"))
  #   
  #   ggsave(output_file_pdf, p, width = 8, height = 6, dpi = 300, limitsize = FALSE)
  #   ggsave(output_file_png, p, width = 8, height = 6, dpi = 300, limitsize = FALSE)
  #   
  #   cat("  Saved:", basename(output_file_pdf), "\n")
  # }
  # 
  # cat("\nAll individual forest plots saved to:", validation_dir, "\n")

  
  # Optional: Create a combined plot with all proteins (faceted by protein)
  cat("\nCreating combined overview plot...\n")
  p_combined <- ggplot(combined_data, aes(x = Estimate, y = traits, color = cohort)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
    geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), 
                   height = 0.3, 
                   position = position_dodge(width = 0.5),
                   size = 0.6) +
    geom_point(position = position_dodge(width = 0.5), 
               size = 2, 
               shape = 18) +
    # Add triangle markers for truncated values
    geom_point(data = combined_data %>% filter(truncated_lower),
               aes(x = -1, y = traits),
               shape = 17, size = 1.5, alpha = 0.6,
               position = position_dodge(width = 0.5)) +
    geom_point(data = combined_data %>% filter(truncated_upper),
               aes(x = 1, y = traits),
               shape = 17, size = 1.5, alpha = 0.6,
               position = position_dodge(width = 0.5)) +
    scale_color_manual(
      values = c("Discovery" = "#337CB2", "Validation" = "#D64643"),
      name = "Cohort"
    ) +
    scale_y_discrete(
      limits = rev(trait_order),
      labels = trait_labels[rev(trait_order)],
      drop = FALSE,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    facet_wrap(~Entrez.Gene.Name, ncol = 4) +
    labs(
      x = "Effect Size (95% CI) [Truncated at ±1]",
      y = "Trait",
      title = paste0("Forest Plot Overview: All Proteins (", current_model, ")"),
      subtitle = "Discovery vs Validation Cohort Effects (Triangles indicate truncated CIs)"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 7),
      axis.title = element_text(size = 10, face = "bold"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9),
      strip.text = element_text(face = "bold", size = 8),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.spacing.x = grid::unit(0.4, "lines"),
      panel.spacing.y = grid::unit(0.2, "lines")
    )
  
  # Save combined overview plot
  output_file_combined <- file.path(validation_dir, "forest_plot_all_proteins_overview.pdf")
  output_file_combined_png <- file.path(validation_dir, "forest_plot_all_proteins_overview.png")
  
  n_proteins <- length(unique_proteins)
  plot_height <- max(5.5, ceiling(n_proteins / 3) * 1.8)
  
  ggsave(output_file_combined, p_combined, width = 12, height = plot_height, dpi = 300, limitsize = FALSE)
  ggsave(output_file_combined_png, p_combined, width = 12, height = plot_height, dpi = 300, limitsize = FALSE)
  
  cat("Combined overview plot saved to:\n")
  cat("  PDF:", output_file_combined, "\n")
  cat("  PNG:", output_file_combined_png, "\n")

  
  # Create validation cohort overview plot (proteins with shared significant traits)
  cat("\nIdentifying proteins with shared significant traits (P < 0.05 in both cohorts)...\n")
  
  validation_data <- combined_data %>%
    filter(cohort == "Validation")
  discovery_data <- combined_data %>%
    filter(cohort == "Discovery")
  
  validation_significant <- validation_data %>%
    filter(!is.na(P_value) & P_value < 0.05) %>%
    select(Entrez.Gene.Name, traits)
  discovery_significant <- discovery_data[discovery_data$adjust_P < 0.05,] %>%
    select(Entrez.Gene.Name, traits)
  
  shared_significant <- inner_join(
    discovery_significant,
    validation_significant,
    by = c("Entrez.Gene.Name", "traits")
  ) %>%
    distinct()
  
  validated_proteins <- shared_significant %>%
    group_by(Entrez.Gene.Name) %>%
    summarise(
      n_validated_traits = n_distinct(traits),
      validated_traits = paste(sort(unique(traits)), collapse = ", "),
      .groups = "drop"
    ) %>%
    arrange(desc(n_validated_traits), Entrez.Gene.Name)
  
  cat(sprintf("Found %d proteins with ≥1 trait significant in both cohorts:\n", nrow(validated_proteins)))
  for (i in 1:nrow(validated_proteins)) {
    cat(sprintf("  %s: %d validated trait(s) (%s)\n",
                validated_proteins$Entrez.Gene.Name[i], 
                validated_proteins$n_validated_traits[i],
                validated_proteins$validated_traits[i]))
  }
  
  if (nrow(validated_proteins) > 0) {
    # Filter combined data for validated proteins only
    validated_combined_data <- combined_data %>%
      filter(Entrez.Gene.Name %in% validated_proteins$Entrez.Gene.Name)
    
    # Create validated proteins overview plot
    cat("\nCreating validated proteins overview plot...\n")
    p_validated <- ggplot(validated_combined_data, aes(x = Estimate, y = traits, color = cohort)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
      geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), 
                     height = 0.3, 
                     position = position_dodge(width = 0.5),
                     size = 0.6) +
      geom_point(position = position_dodge(width = 0.5), 
                 size = 2, 
                 shape = 18) +
      # Add triangle markers for truncated values
      geom_point(data = validated_combined_data %>% filter(truncated_lower),
                 aes(x = -1, y = traits),
                 shape = 17, size = 1.5, alpha = 0.6,
                 position = position_dodge(width = 0.5)) +
      geom_point(data = validated_combined_data %>% filter(truncated_upper),
                 aes(x = 1, y = traits),
                 shape = 17, size = 1.5, alpha = 0.6,
                 position = position_dodge(width = 0.5)) +
      scale_color_manual(
        values = c("Discovery" = "#337CB2", "Validation" = "#D64643"),
        name = "Cohort"
      ) +
      scale_y_discrete(
        limits = rev(trait_order),
        labels = trait_labels[rev(trait_order)],
        drop = FALSE,
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
      facet_wrap(~Entrez.Gene.Name, ncol = 4) +
      labs(
        x = "Effect Size (95% CI) [Truncated at ±1]",
        y = "Trait",
        title = paste0("Forest Plot Overview: Validated Proteins (", current_model, ")"),
        subtitle = sprintf("Proteins with ≥1 trait significant in both cohorts (n=%d) - Triangles indicate truncated CIs", nrow(validated_proteins))
      ) +
      theme_bw() +
      theme(
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 10, face = "bold"),
        legend.position = "top",
        legend.title = element_text(face = "bold"),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9),
        strip.text = element_text(face = "bold", size = 8),
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        panel.spacing.x = grid::unit(0.4, "lines"),
        panel.spacing.y = grid::unit(0.2, "lines")
      )
    
    # Save validated proteins overview plot
    output_file_validated <- file.path(validation_dir, "forest_plot_all_proteins_overview_validated.pdf")
    output_file_validated_png <- file.path(validation_dir, "forest_plot_all_proteins_overview_validated.png")
    
    n_validated <- nrow(validated_proteins)
    plot_height_validated <- max(5.5, ceiling(n_validated / 3) * 1.8)
    
    ggsave(output_file_validated, p_validated, width = 12, height = plot_height_validated, dpi = 300, limitsize = FALSE)
    ggsave(output_file_validated_png, p_validated, width = 12, height = plot_height_validated, dpi = 300, limitsize = FALSE)
    
    cat("\nValidated proteins overview plot saved to:\n")
    cat("  PDF:", output_file_validated, "\n")
    cat("  PNG:", output_file_validated_png, "\n")
    
    # Export validated proteins list
    validated_list_file <- file.path(validation_dir, "validated_proteins_list.txt")
    write.table(validated_proteins, validated_list_file, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("  Validated proteins list:", validated_list_file, "\n")

    # Prepare plots splitting proteins by effect direction
    cat("\nPreparing effect-direction plots for validated proteins...\n")
    trait_selection <- validated_combined_data %>%
      group_by(Entrez.Gene.Name, traits) %>%
      summarise(
        validation_p = {
          vals <- P_value[cohort == "Validation"]
          if (length(vals) == 0 || all(is.na(vals))) NA_real_ else min(vals, na.rm = TRUE)
        },
        discovery_p = {
          vals <- adjust_P[cohort == "Discovery"]
          if (length(vals) == 0 || all(is.na(vals))) NA_real_ else min(vals, na.rm = TRUE)
        },
        .groups = "drop"
      ) %>%
      mutate(primary_p = coalesce(validation_p, discovery_p, Inf)) %>%
      group_by(Entrez.Gene.Name) %>%
      slice_min(order_by = primary_p, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(Entrez.Gene.Name, traits)

    selected_trait_data <- validated_combined_data %>%
      semi_join(trait_selection, by = c("Entrez.Gene.Name", "traits"))

    protein_effect_sign <- selected_trait_data %>%
      group_by(Entrez.Gene.Name) %>%
      summarise(
        effect_sign = {
          est_discovery <- Estimate_original[cohort == "Discovery"]
          est_discovery <- est_discovery[!is.na(est_discovery)]
          reference_est <- if (length(est_discovery) > 0) {
            est_discovery[1]
          } else {
            est_validation <- Estimate_original[cohort == "Validation"]
            est_validation <- est_validation[!is.na(est_validation)]
            if (length(est_validation) > 0) est_validation[1] else NA_real_
          }
          if (is.na(reference_est) || reference_est == 0) {
            NA_character_
          } else if (reference_est > 0) {
            "Positive"
          } else {
            "Negative"
          }
        },
        .groups = "drop"
      ) %>%
      filter(!is.na(effect_sign))

    positive_proteins <- protein_effect_sign %>%
      filter(effect_sign == "Positive") %>%
      pull(Entrez.Gene.Name)

    negative_proteins <- protein_effect_sign %>%
      filter(effect_sign == "Negative") %>%
      pull(Entrez.Gene.Name)

    positive_effect_data <- validated_combined_data %>%
      filter(Entrez.Gene.Name %in% positive_proteins)

    negative_effect_data <- validated_combined_data %>%
      filter(Entrez.Gene.Name %in% negative_proteins)

    build_effect_plot <- function(data, sign_label) {
      if (nrow(data) == 0) {
        return(NULL)
      }

      n_panels <- length(unique(data$Entrez.Gene.Name))
      data <- data %>%
        mutate(traits = factor(traits, levels = trait_order))

      plot_title <- paste0("Validated Proteins (", sign_label, " Effects) – ", current_model)
      plot_subtitle <- "Proteins meeting shared significance; triangles mark truncated CIs"

      p_effect <- ggplot(data, aes(x = Estimate, y = traits, color = cohort)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
        geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci),
                       height = 0.3,
                       position = position_dodge(width = 0.5),
                       size = 0.6) +
        geom_point(position = position_dodge(width = 0.5),
                   size = 2,
                   shape = 18) +
        geom_point(data = data %>% filter(truncated_lower),
                   aes(x = -1, y = traits),
                   shape = 17, size = 1.5, alpha = 0.6,
                   position = position_dodge(width = 0.5)) +
        geom_point(data = data %>% filter(truncated_upper),
                   aes(x = 1, y = traits),
                   shape = 17, size = 1.5, alpha = 0.6,
                   position = position_dodge(width = 0.5)) +
        scale_color_manual(
          values = c("Discovery" = "#337CB2", "Validation" = "#D64643"),
          name = "Cohort"
        ) +
        scale_y_discrete(
          limits = rev(trait_order),
          labels = trait_labels[rev(trait_order)],
          drop = FALSE,
          expand = expansion(mult = c(0.02, 0.02))
        ) +
        scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
        facet_wrap(~Entrez.Gene.Name, ncol = 2) +
        labs(
          x = "Effect Size (95% CI) [Truncated at ±1]",
          y = "Trait",
          title = plot_title,
          subtitle = plot_subtitle
        ) +
        theme_bw() +
        theme(
          axis.text.y = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          axis.title = element_text(size = 10, face = "bold"),
          legend.position = "top",
          legend.title = element_text(face = "bold"),
          plot.title = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 9),
          strip.text = element_text(face = "bold", size = 8),
          panel.grid.major.y = element_line(color = "gray90"),
          panel.grid.minor = element_blank(),
          panel.spacing.x = grid::unit(0.4, "lines"),
          panel.spacing.y = grid::unit(0.2, "lines")
        )

      list(
        plot = p_effect,
        height = max(5.5, ceiling(n_panels / 2) * 1.8),
        sign_label = sign_label
      )
    }

    positive_plot <- build_effect_plot(positive_effect_data, "Positive")
    negative_plot <- build_effect_plot(negative_effect_data, "Negative")

    if (is.null(positive_plot)) {
      cat("  No positive effect proteins with shared significance found.\n")
    } else {
      output_pdf_pos <- file.path(validation_dir, "forest_plot_validated_positive_effects.pdf")
      output_png_pos <- file.path(validation_dir, "forest_plot_validated_positive_effects.png")
      ggsave(output_pdf_pos, positive_plot$plot, width = 12, height = positive_plot$height, dpi = 300, limitsize = FALSE)
      ggsave(output_png_pos, positive_plot$plot, width = 12, height = positive_plot$height, dpi = 300, limitsize = FALSE)
      cat("  Positive effect plot saved to:\n")
      cat("    PDF:", output_pdf_pos, "\n")
      cat("    PNG:", output_png_pos, "\n")
    }

    if (is.null(negative_plot)) {
      cat("  No negative effect proteins with shared significance found.\n")
    } else {
      output_pdf_neg <- file.path(validation_dir, "forest_plot_validated_negative_effects.pdf")
      output_png_neg <- file.path(validation_dir, "forest_plot_validated_negative_effects.png")
      ggsave(output_pdf_neg, negative_plot$plot, width = 12, height = negative_plot$height, dpi = 300, limitsize = FALSE)
      ggsave(output_png_neg, negative_plot$plot, width = 12, height = negative_plot$height, dpi = 300, limitsize = FALSE)
      cat("  Negative effect plot saved to:\n")
      cat("    PDF:", output_pdf_neg, "\n")
      cat("    PNG:", output_png_neg, "\n")
    }

    if (!is.null(positive_plot) || !is.null(negative_plot)) {
      combined_plot <- (
        if (!is.null(positive_plot)) positive_plot$plot else patchwork::plot_spacer()
      ) + (
        if (!is.null(negative_plot)) negative_plot$plot else patchwork::plot_spacer()
      )

      combined_plot <- combined_plot +
        plot_layout(ncol = 2, guides = "collect") &
        theme(legend.position = "top")

      combined_height <- max(
        ifelse(is.null(positive_plot), 0, positive_plot$height),
        ifelse(is.null(negative_plot), 0, negative_plot$height),
        5.5
      )

      output_pdf_combined <- file.path(validation_dir, "forest_plot_validated_effects_side_by_side.pdf")
      output_png_combined <- file.path(validation_dir, "forest_plot_validated_effects_side_by_side.png")

      ggsave(output_pdf_combined, combined_plot, width = 12, height = combined_height, dpi = 300, limitsize = FALSE)
      ggsave(output_png_combined, combined_plot, width = 12, height = combined_height, dpi = 300, limitsize = FALSE)

      cat("  Combined positive/negative effect plot saved to:\n")
      cat("    PDF:", output_pdf_combined, "\n")
      cat("    PNG:", output_png_combined, "\n")
    }
  } else {
    cat("No proteins met the shared significance criterion\n")
  }
  
  # Export combined data for reference
  output_data <- file.path(validation_dir, "forest_plot_data.txt")
  write.table(combined_data, output_data, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("\nData exported to:", output_data, "\n")
  
  cat("\n=== Summary for", current_model, "===\n")
  cat("Total proteins plotted:", length(unique_proteins), "\n")
  cat("Validation cohort proteins retained for overview (shared significance):", ifelse(exists("validated_proteins"), nrow(validated_proteins), 0), "\n")
  cat("Total trait categories:", length(trait_order), "\n")
  cat("Plots saved to:", validation_dir, "\n")
}

cat("\n", rep("=", 80), "\n", sep = "")
cat("All models processed successfully!\n")
cat(rep("=", 80), "\n", sep = "")
