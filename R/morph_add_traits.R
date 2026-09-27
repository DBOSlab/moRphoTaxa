#' Add morphological trait columns to a specimen data frame
#'
#' @author
#' João Dornelas & Domingos Cardoso
#'
#' @description
#' Creates a standardized morphological trait matrix for morphometric analyses.
#' The function can either append morphological trait columns to an existing
#' specimen-occurrence data frame or create an empty trait matrix when
#' `base_df = NULL`.
#'
#' Trait columns can include a general set of qualitative and quantitative
#' characters applicable across flowering plants, as well as additional
#' taxon-specific characters selected from predefined morphological presets.
#' Trait names may be returned as full descriptive names, abbreviated trait
#' codes, or a combination of both.
#'
#' @details
#' `morph_add_traits()` is intended to prepare standardized data matrices for
#' downstream morphometric analyses in **moRphoTaxa**. The function
#' distinguishes between two major classes of morphological characters:
#'
#' * **Qualitative traits**, such as leaf shape, indumentum, floral morphology,
#'   fruit type, and seed morphology.
#' * **Quantitative traits**, such as lengths, widths, areas, numbers, and
#'   other measurements of vegetative and reproductive structures.
#'
#' Default trait sets contain characters broadly applicable to flowering
#' plants. Additional taxon-specific traits can be requested independently for
#' qualitative and quantitative characters using `quali_specific` and
#' `quanti_specific`. Taxon-specific presets contain only traits that are not
#' already part of the default sets.
#'
#' Available taxon-specific presets are `"araceae"`, `"asteraceae"`,
#' `"bryophyta"`, `"cactaceae"`, `"caesalpinioideae"`, `"lycopodiopsida"`,
#' `"orchidaceae"`, `"papilionoideae"`, `"passifloraceae"`, `"pinophyta"`,
#' `"poaceae"` and `"polypodiopsida"`.
#'
#' When an existing data frame is supplied through `base_df`, all original
#' columns are retained and the selected morphological trait columns are added.
#' Existing trait columns are preserved, whereas missing trait columns are
#' initialized with `NA`.
#'
#' When `base_df = NULL`, the function returns an empty data frame containing
#' only the selected morphological trait columns. This can be used as a
#' template for manual morphological data entry.
#'
#' Trait columns are arranged according to a standardized sequence of plant
#' structures, from vegetative to reproductive characters.
#'
#' If `save = TRUE`, the resulting data frame is additionally written to an
#' Excel (`.xlsx`) file.
#'
#' @param base_df A data frame containing specimen or occurrence records to
#' which morphological trait columns will be added. If `NULL`, an empty data
#' frame containing only the selected morphological trait columns is created.
#' Default is `NULL`.
#'
#' @param trait_name Character string specifying how morphological trait column
#' names should be represented. One of:
#' \itemize{
#'   \item `"trait full name"`: use descriptive trait names, such as
#'   `"leaf_length"`;
#'   \item `"trait code"`: use abbreviated trait codes, such as `"LEAFlng"`;
#'   \item `"both"`: combine the descriptive name and abbreviation, for example
#'   `"leaf_length/LEAFlng"`.
#' }
#' If omitted, `"trait full name"` is used.
#'
#' @param quali_default Logical. If `TRUE` (default), includes the default set
#' of qualitative morphological traits broadly applicable to flowering plants.
#' Set to `FALSE` to exclude these characters. Automatically treated as
#' `FALSE` when `quali_specific` refers to a non-angiosperm group (currently
#' `"bryophyta"`, `"lycopodiopsida"`, `"pinophyta"` or `"polypodiopsida"`),
#' since these default traits (petal, sepal, fruit, seed, etc.) do not apply
#' to those lineages.
#'
#' @param quanti_default Logical. If `TRUE` (default), includes the default set
#' of quantitative morphological traits broadly applicable to flowering plants.
#' Set to `FALSE` to exclude these characters. Automatically treated as
#' `FALSE` when `quanti_specific` refers to a non-angiosperm group (currently
#' `"bryophyta"`, `"lycopodiopsida"`, `"pinophyta"` or `"polypodiopsida"`),
#' since these default traits do not apply to those lineages.
#'
#' @param quali_specific Character string specifying a taxon-specific preset
#' of qualitative morphological traits to add (case-insensitive). See
#' *Details* for the available presets. If `NULL` (default), no
#' taxon-specific qualitative traits are added.
#'
#' @param quanti_specific Character string specifying a taxon-specific preset
#' of quantitative morphological traits to add (case-insensitive). See
#' *Details* for the available presets. If `NULL` (default), no
#' taxon-specific quantitative traits are added.
#'
#' @param save Logical. If `TRUE` (default), writes the resulting morphological
#' data matrix to an Excel (`.xlsx`) file. If `FALSE`, the data frame is
#' returned without writing a file.
#'
#' @param file_name Character string specifying the file path (including file
#' name and `.xlsx` extension) used to save the resulting data frame when
#' `save = TRUE`. If `NULL` (default), a file is automatically saved inside an
#' `output_data/` folder (created if it does not already exist), with a file
#' name generated from the current date and time. Ignored if `save = FALSE`.
#'
#' @param verbose Logical. If `TRUE` (default), displays informative messages
#' during execution, including the location of exported files.
#'
#' @return
#' A data frame containing the original columns of `base_df`, when supplied,
#' together with the selected morphological trait columns. Missing trait
#' values are initialized as `NA`. If `base_df = NULL`, an empty data frame
#' containing only the selected trait columns is returned.
#'
#' If `save = TRUE`, the same data frame is also written to an Excel file.
#'
#' @examples
#' \dontrun{
#' # Add default qualitative and quantitative traits
#' traits <- morph_add_traits(
#'   base_df = specimen_data,
#'   trait_name = "trait full name",
#'   save = FALSE
#' )
#'
#' # Add Passifloraceae-specific qualitative and quantitative traits
#' traits_passifloraceae <- morph_add_traits(
#'   base_df = specimen_data,
#'   trait_name = "both",
#'   quali_specific = "passifloraceae",
#'   quanti_specific = "passifloraceae",
#'   save = FALSE
#' )
#'
#' # Add Poaceae-specific traits on top of the default traits
#' traits_poaceae <- morph_add_traits(
#'   base_df = specimen_data,
#'   trait_name = "both",
#'   quali_specific = "poaceae",
#'   quanti_specific = "poaceae",
#'   save = FALSE
#' )
#'
#' # Add Araceae-specific traits on top of the default traits
#' traits_araceae <- morph_add_traits(
#'   base_df = specimen_data,
#'   trait_name = "both",
#'   quali_specific = "araceae",
#'   quanti_specific = "araceae",
#'   save = FALSE
#' )
#'
#' # Add conifer-specific traits (angiosperm defaults are excluded automatically)
#' traits_pinophyta <- morph_add_traits(
#'   base_df = specimen_data,
#'   trait_name = "both",
#'   quali_specific = "pinophyta",
#'   quanti_specific = "pinophyta",
#'   save = FALSE
#' )
#'
#' # Add lycophyte-specific traits (angiosperm defaults are excluded automatically)
#' traits_lycopodiopsida <- morph_add_traits(
#'   base_df = specimen_data,
#'   trait_name = "both",
#'   quali_specific = "lycopodiopsida",
#'   quanti_specific = "lycopodiopsida",
#'   save = FALSE
#' )
#'
#' # Create an empty template containing only morphological traits
#' trait_template <- morph_add_traits(
#'   base_df = NULL,
#'   trait_name = "trait code",
#'   quali_default = TRUE,
#'   quanti_default = TRUE,
#'   save = FALSE
#' )
#' }
#'
#' @importFrom openxlsx write.xlsx
#'
#' @export
#'

morph_add_traits <- function(base_df = NULL,
                           trait_name = c("trait full name",
                                          "trait code",
                                          "both"),
                           quali_default = TRUE,
                           quanti_default = TRUE,
                           quali_specific = NULL,
                           quanti_specific = NULL,
                           save = TRUE,
                           file_name = NULL,
                           verbose = TRUE) {

#_______________________________________________________________________________

# Qualitative traits ####

# Default traits ####
quali_default_traits <- c(
  # vegetative traits
  "stem_indumentum" = "STEMind",
  "stem_type" = "STEMtyp",
  "petiole_indumentum" = "PETIind",
  "petiole_form" = "PETIfrm",
  "petiole_surface" = "PETIsrf",
  "stipule_form" = "STIPfrm",
  "stipule_base" = "STIPbas",
  "stipule_margin" = "STIPmrg",
  "stipule_apex" = "STIPapx",
  "stipule_indumentum" = "STIPind",
  "stipule_persistence" = "STIPper",
  "leaf_arrangement" = "LEAFarr",
  "leaf_division" = "LEAFdiv",
  "leaf_texture" = "LEAFtxt",
  "leaf_form" = "LEAFfrm",
  "leaf_base" = "LEAFbas",
  "leaf_margin" = "LEAFmrg",
  "leaf_apex" = "LEAFapx",
  "leaf_indumentum" = "LEAFind",
  "venation_type" = "VENAtyp",

  # reproductive traits
  "inflorescence_type" = "INFLtyp",
  "inflorescence_indumentum" = "INFLind",
  "bract_disposition" = "BRCTdsp",
  "bract_form" = "BRCTfrm",
  "bract_base" = "BRCTbas",
  "bract_margin" = "BRCTmrg",
  "bract_apex" = "BRCTapx",
  "bract_color" = "BRCTclr",
  "bract_texture" = "BRCTtxt",
  "bract_indumentum" = "BRCTind",
  "bract_venation" = "BRCTven",
  "bracteole_form" = "BRCLfrm",
  "sepal_form" = "SEPLfrm",
  "sepal_fusion" = "SEPLfus",
  "sepal_color" = "SEPLclr",
  "sepal_indumentum" = "SEPLind",
  "sepal_symmetry" = "SEPLsym",
  "sepal_apex" = "SEPLapx",
  "sepal_margin" = "SEPLmrg",
  "sepal_texture" = "SEPLtxt",
  "sepal_lobe_form" = "SEPLlobfrm",
  "petal_apex" = "PETLapx",
  "petal_margin" = "PETLmrg",
  "petal_orientation" = "PETLori",
  "petal_form" = "PETLfrm",
  "petal_fusion" = "PETLfus",
  "petal_color" = "PETLclr",
  "petal_indumentum" = "PETLind",
  "petal_symmetry" = "PETLsym",
  "petal_aestivation" = "PETLaes",
  "petal_lobe_form" = "PETLlobfrm",
  "ovary_position" = "OVARpos",
  "ovary_indumentum" = "OVARind",
  "placentation_type" = "PLCNtyp",
  "style_form" = "STYLfrm",
  "stigma_form" = "STIGfrm",
  "stamen_fusion" = "ANDRfus",
  "filament_form" = "FLMTfrm",
  "anther_form" = "ANTHfrm",
  "anther_attachment" = "ANTHatt",
  "anther_dehiscence" = "ANTHdeh",
  "anther_color" = "ANTHclr",
  "fruit_type" = "FRUTtyp",
  "fruit_dehiscence" = "FRUTdeh",
  "fruit_form" = "FRUTfrm",
  "fruit_color" = "FRUTclr",
  "fruit_texture" = "FRUTtxt",
  "fruit_indumentum" = "FRUTind",
  "seed_form" = "SEEDfrm",
  "seed_color" = "SEEDclr",
  "seed_surface" = "SEEDsrf"
)

# Specific traits ####

# Araceae ####
quali_araceae <- c(
  # vegetative traits
  "cataphyll_persistence" = "CATAper",
  "petiole_sheath_form" = "PSHTfrm",
  "geniculum_form" = "GENIfrm",
  "primary_lateral_vein_arrangement" = "VENAlatarr",
  "collective_vein_type" = "VENAcol",

  # reproductive traits
  "peduncle_form" = "PDNCfrm",
  "spathe_form" = "SPTHfrm",
  "spathe_apex" = "SPTHapx",
  "spathe_outer_color" = "SPTHclrout",
  "spathe_inner_color" = "SPTHclrinn",
  "spathe_orientation" = "SPTHori",
  "spadix_form" = "SPDXfrm",
  "spadix_attachment" = "SPDXatt",
  "spadix_appendix_form" = "SPDXappfrm",
  "spadix_appendix_color" = "SPDXappclr",
  "flower_arrangement" = "FLWRarr",
  "perigone_form" = "PRGNfrm",
  "perigone_fusion" = "PRGNfus",
  "staminode_form" = "STMDfrm",
  "ovary_form" = "OVARfrm"
)

# Asteraceae ####
quali_asteraceae <- c(
  # reproductive traits
  "capitulescence_type" = "CPEStyp",
  "capitule_type" = "CAPTtyp",
  "phyllary_seriation" = "PHYLser",
  "phyllary_form" = "PHYLfrm",
  "phyllary_margin" = "PHYLmrg",
  "phyllary_apex" = "PHYLapx",
  "phyllary_texture" = "PHYLtxt",
  "phyllary_indumentum" = "PHYLind",
  "phyllary_color" = "PHYLclr",
  "receptacle_form" = "RECPfrm",
  "receptacle_paleae_form" = "RECPplfrm",
  "ray_floret_symmetry" = "RAYFsym",
  "ray_floret_color" = "RAYFclr",
  "ray_floret_apex" = "RAYFapx",
  "disc_floret_symmetry" = "DISFsym",
  "disc_floret_color" = "DISFclr",
  "disc_petal_lobe_form" = "DISFlobfrm",
  "anther_base_form" = "ANTHbasfrm",
  "anther_appendage_form" = "ANTHappfrm",
  "style_branch_form" = "STYLbrfrm",
  "style_branch_apex" = "STYLbrapx",
  "cypsela_form" = "CYPSfrm",
  "cypsela_ribbing" = "CYPSrib",
  "cypsela_indumentum" = "CYPSind",
  "cypsela_carpopodium_form" = "CYPScarfrm",
  "pappus_type" = "PAPPtyp",
  "pappus_color" = "PAPPclr",
  "pappus_persistence" = "PAPPper"
)

# Bromeliaceae ####
quali_bromeliaceae <- c(
  # vegetative traits
  "leaf_sheath_form" = "SHTHfrm",
  "leaf_sheath_color" = "SHTHclr",
  "leaf_spine_orientation" = "LSPNori",
  "leaf_spine_color" = "LSPNclr",
  "trichome_scale_form" = "TRCHfrm",
  "trichome_scale_margin" = "TRCHmrg",
  "trichome_scale_color" = "TRCHclr",
  "trichome_scale_distribution" = "TRCHdis",

  # reproductive traits
  "scape_form" = "SCAPfrm",
  "scape_bract_form" = "SCBRfrm",
  "scape_bract_disposition" = "SCBRdsp",
  "scape_bract_apex" = "SCBRapx",
  "scape_bract_margin" = "SCBRmrg",
  "scape_bract_color" = "SCBRclr",
  "scape_bract_indumentum" = "SCBRind",
  "primary_bract_form" = "PBRCfrm",
  "primary_bract_apex" = "PBRCapx",
  "primary_bract_margin" = "PBRCmrg",
  "primary_bract_color" = "PBRCclr",
  "primary_bract_indumentum" = "PBRCind",
  "floral_bract_form" = "FBRCfrm",
  "floral_bract_apex" = "FBRCapx",
  "floral_bract_margin" = "FBRCmrg",
  "floral_bract_color" = "FBRCclr",
  "floral_bract_indumentum" = "FBRCind",
  "petal_appendage_form" = "PETLapdfrm",
  "nectary_form" = "NCTRfrm",
  "seed_appendage_form" = "SEEDapdfrm"
)

# Bryophyta ####
quali_bryophyta <- c(
  # vegetative traits
  "rhizoid_color" = "RHIZclr",
  "phyllid_insertion" = "PHYDins",
  "phyllid_form" = "PHYDfrm",
  "phyllid_apex" = "PHYDapx",
  "phyllid_base" = "PHYDbas",
  "phyllid_margin" = "PHYDmrg",
  "costa_number" = "COSTnmb",
  "costa_ending" = "COSTend",
  "phyllid_cell_shape" = "PCELshp",
  "phyllid_cell_wall" = "PCELwal",
  "alar_cell_differentiation" = "ALARdif",

  # reproductive traits
  "seta_color" = "SETAclr",
  "seta_surface" = "SETAsrf",
  "capsule_orientation" = "CAPSori",
  "capsule_form" = "CAPSfrm",
  "operculum_form" = "OPCLfrm",
  "peristome_type" = "PERItyp",
  "peristome_ornamentation" = "PERIorn",
  "calyptra_form" = "CLPTfrm",
  "calyptra_surface" = "CLPTsrf",
  "spore_ornamentation" = "SPORorn",
  "spore_color" = "SPORclr"
)

# Cactaceae ####
quali_cactaceae <- c(
  # vegetative traits
  "stem_segmentation" = "STEMseg",
  "rib_arrangement" = "RIBSarr",
  "areole_wool_color" = "ARLWclr",
  "spine_arrangement" = "SPINarr",
  "central_spine_form" = "CSPNfrm",
  "radial_spine_form" = "RSPNfrm",
  "spine_color" = "SPINclr",
  "glochid_color" = "GLOCclr",
  "tubercle_form" = "TUBRfrm",

  # reproductive traits
  "perigone_form" = "PRGNfrm",
  "perigone_fusion" = "PRGNfus",
  "perigone_color" = "PRGNclr",
  "perigone_indumentum" = "PRGNind",
  "pericarpel_form" = "PRCPfrm",
  "fruit_scale_form" = "FSCLfrm"
)

# Caesalpinioideae ####
quali_caesalpinioideae <- c(
  # vegetative traits
  "leaflet_disposition" = "LFLTdsp",
  "leaflet_form" = "LFLTfrm",
  "pinnule_disposition" = "PINLdsp",
  "pinnule_form" = "PINLfrm"
)

# Eriocaulaceae ####
quali_eriocaulaceae <- c(
  # vegetative traits
  "scape_ribbing" = "SCAPrib",
  "scape_sheath_form" = "SCAPshtfrm",
  "scape_indumentum" = "SCAPind",
  "scape_color" = "SCAPclr",

  # reproductive traits
  "capitule_shape" = "CAPTshp",
  "capitule_color" = "CAPTclr",
  "involucral_bract_seriation" = "INVLser",
  "involucral_bract_form" = "INVLfrm",
  "involucral_bract_texture" = "INVLtxt",
  "involucral_bract_color" = "INVLclr",
  "involucral_bract_indumentum" = "INVLind",
  "receptacle_indumentum" = "RECPind",
  "palea_form" = "PALEfrm",
  "palea_margin" = "PALEmrg",
  "palea_color" = "PALEclr",
  "staminate_sepal_fusion" = "STFLseplfus",
  "staminate_petal_fusion" = "STFLpetlfus",
  "staminate_petal_gland_form" = "STFLglafrm",
  "staminate_petal_gland_color" = "STFLglaclr",
  "pistillate_sepal_fusion" = "PSFLseplfus",
  "pistillate_petal_fusion" = "PSFLpetlfus",
  "pistillate_petal_gland_form" = "PSFLglafrm",
  "pistillate_petal_gland_color" = "PSFLglaclr",
  "stigma_branch_form" = "STIGbrfrm",
  "seed_appendage_arrangement" = "SEEDapdarr",
  "seed_appendage_form" = "SEEDapdfrm"
)

# Euphorbiaceae ####
quali_euphorbiaceae <- c(
  # vegetative traits
  "latex_color" = "LATXclr",
  "stem_cross_section_form" = "STEMxsc",
  "leaf_gland_position" = "LGLDpos",
  "leaf_gland_form" = "LGLDfrm",
  "spine_origin" = "SPINori",
  "spine_form" = "SPINfrm",
  "spine_color" = "SPINclr",

  # reproductive traits
  "cyathium_arrangement" = "CYATarr",
  "cyathium_symmetry" = "CYATsym",
  "cyathium_involucre_form" = "CYATinvfrm",
  "cyathium_involucre_indumentum" = "CYATinvind",
  "cyathium_involucre_lobe_form" = "CYATlobfrm",
  "cyathium_involucre_lobe_margin" = "CYATlobmrg",
  "cyathium_gland_form" = "CYATglafrm",
  "cyathium_gland_color" = "CYATglaclr",
  "cyathium_gland_appendage_form" = "CYATglaappfrm",
  "cyathophyll_form" = "CYPHfrm",
  "bract_gland_form" = "BRCTglafrm",
  "staminate_perianth_form" = "STFLperfrm",
  "pistillate_perianth_form" = "PSFLperfrm",
  "pistillate_flower_pedicel_form" = "PSFLpedfrm",
  "filament_articulation" = "FLMTart",
  "style_branch_form" = "STYLbrfrm",
  "style_branch_apex" = "STYLbrapx",
  "fruit_lobing" = "FRUTlob",
  "coccus_form" = "COCCfrm",
  "columella_form" = "COLUfrm",
  "caruncle_form" = "CARUfrm",
  "caruncle_color" = "CARUclr"
)

# Lecythidaceae ####
quali_lecythidaceae <- c(
  # reproductive traits
  "corolla_androecium_calyptra_form" = "CALYfrm",
  "androecial_ring_symmetry" = "ANDRsym",
  "androecial_hood_form" = "HOODfrm",
  "androecial_hood_orientation" = "HOODori",
  "androecial_hood_color" = "HOODclr",
  "pyxidium_operculum_form" = "PYXOfrm",
  "pyxidium_operculum_position" = "PYXOpos",
  "pericarp_texture" = "PERCtxt",
  "seed_aril_form" = "ARILfrm",
  "seed_wing_form" = "SEEDwngfrm",
  "testa_texture" = "TESTtxt",
  "testa_ornamentation" = "TESTorn",
  "funicle_form" = "FUNIfrm"
)

# Lycopodiopsida ####
quali_lycopodiopsida <- c(
  # vegetative traits
  "stem_habit" = "STEMhab",
  "stem_branching" = "STEMbrn",
  "rhizophore_form" = "RHPHfrm",
  "corm_form" = "CORMfrm",
  "microphyll_arrangement" = "MPHLarr",
  "microphyll_form" = "MPHLfrm",
  "microphyll_apex" = "MPHLapx",
  "microphyll_base" = "MPHLbas",
  "microphyll_margin" = "MPHLmrg",
  "microphyll_texture" = "MPHLtxt",
  "microphyll_color" = "MPHLclr",
  "microphyll_terminal_seta_form" = "MPHLsetfrm",
  "ligule_form" = "LIGLfrm",

  # reproductive traits
  "strobilus_form" = "STRBfrm",
  "strobilus_position" = "STRBpos",
  "strobilus_peduncle_form" = "STRBpedfrm",
  "sporophyll_form" = "SPHLfrm",
  "sporophyll_apex" = "SPHLapx",
  "sporophyll_margin" = "SPHLmrg",
  "sporophyll_color" = "SPHLclr",
  "sporangium_form" = "SPRGfrm",
  "sporangium_position" = "SPRGpos",
  "sporangium_color" = "SPRGclr",
  "velum_extent" = "VELUext",
  "spore_ornamentation" = "SPORorn",
  "spore_color" = "SPORclr",
  "megaspore_ornamentation" = "MSPOorn",
  "megaspore_color" = "MSPOclr",
  "microspore_ornamentation" = "MISPorn",
  "microspore_color" = "MISPclr"
)

# Moraceae ####
quali_moraceae <- c(
  # reproductive traits
  "receptacle_form" = "RECPfrm",
  "receptacle_indumentum" = "RECPind",
  "ostiole_form" = "OSTLfrm",
  "ostiole_scale_arrangement" = "OSTLscaarr",
  "staminate_flower_arrangement" = "STFLarr",
  "staminate_perianth_form" = "STFLperfrm",
  "staminate_perianth_fusion" = "STFLperfus",
  "pistillate_flower_arrangement" = "PSFLarr",
  "pistillate_perianth_form" = "PSFLperfrm",
  "pistillate_perianth_fusion" = "PSFLperfus",
  "gall_flower_type" = "GLFLtyp",
  "stigma_branch_form" = "STIGbrfrm",
  "syncarp_form" = "SYNCfrm",
  "syncarp_color" = "SYNCclr",
  "syncarp_texture" = "SYNCtxt",
  "drupelet_arrangement" = "DRUParr",
  "achenium_enclosure" = "ACHEenc"
)

# Orchidaceae ####
quali_orchidaceae <- c(
  # vegetative traits
  "pseudobulb_form" = "PSBLfrm",
  "pseudobulb_color" = "PSBLclr",
  "root_velamen_color" = "VELAclr",

  # reproductive traits
  "labellum_form" = "LABLfrm",
  "labellum_margin" = "LABLmrg",
  "labellum_color" = "LABLclr",
  "column_form" = "COLMfrm",
  "anther_cap_form" = "ANTCfrm",
  "pollinia_form" = "POLNfrm"
)

# Papilionoideae ####
quali_papilionoideae <- c(
  # vegetative traits
  "leaflet_disposition" = "LFLTdsp",
  "leaflet_form" = "LFLTfrm",

  # reproductive traits
  "standard_petal_form" = "STNDfrm",
  "wing_petal_form" = "WINGfrm",
  "keel_petal_form" = "KEELfrm",
  "valve_consistency" = "VALVcon",
  "hilum_form" = "HILUfrm"
)

# Passifloraceae ####
quali_passifloraceae <- c(
  # vegetative traits
  "petiole_nectary_form" = "PTNCfrm",
  "foliar_nectary_form" = "FLNCfrm",
  "tendril_consistency" = "TNDRcon",

  # reproductive traits
  "hypanthium_form" = "HIPTfrm",
  "coronal_filament_form" = "CRFLfrm",
  "coronal_filament_color" = "CRFLclr",
  "operculum_form" = "OPERfrm",
  "operculum_texture" = "OPERtxt",
  "operculum_orientation" = "OPERori",
  "operculum_color" = "OPERclr",
  "limen_form" = "LIMNfrm",
  "limen_texture" = "LIMNtxt"
)

# Pinophyta ####
quali_pinophyta <- c(
  # vegetative traits
  "bark_texture" = "BARKtxt",
  "bark_color" = "BARKclr",
  "crown_form" = "CRWNfrm",
  "branching_pattern" = "BRNCpat",
  "twig_form" = "TWIGfrm",
  "twig_indumentum" = "TWIGind",
  "needle_arrangement" = "NEEDarr",
  "leaf_type" = "LEAFtyp",
  "needle_form" = "NEEDfrm",
  "needle_apex" = "NEEDapx",
  "needle_margin" = "NEEDmrg",
  "needle_stomatal_band_position" = "NEEDstmpos",
  "needle_color" = "NEEDclr",

  # reproductive traits
  "pollen_cone_arrangement" = "POLCarr",
  "pollen_cone_color" = "POLCclr",
  "seed_cone_position" = "SEECpos",
  "seed_cone_form" = "SEECfrm",
  "seed_cone_color" = "SEECclr",
  "seed_cone_dehiscence" = "SEECdeh",
  "cone_scale_form" = "SCALfrm",
  "cone_scale_apex" = "SCALapx",
  "cone_scale_texture" = "SCALtxt",
  "cone_scale_indumentum" = "SCALind",
  "bract_scale_form" = "BRSCfrm",
  "apophysis_form" = "APOPfrm",
  "umbo_form" = "UMBOfrm",
  "umbo_position" = "UMBOpos",
  "aril_color" = "ARILclr",
  "seed_wing_form" = "SEEDwngfrm",
  "seed_wing_articulation" = "SEEDwngart"
)

# Piperaceae ####
quali_piperaceae <- c(
  # vegetative traits
  "leaf_pellucid_gland_distribution" = "LGLDpel",
  "leaf_asymmetry" = "LEAFasy",

  # reproductive traits
  "spike_orientation" = "SPIKori",
  "spike_disposition" = "SPIKdsp",
  "floral_bract_form" = "FBRCfrm",
  "floral_bract_margin" = "FBRCmrg",
  "stamen_arrangement" = "STMNarr",
  "anther_dehiscence_orientation" = "ANTHdehori",
  "stigma_arrangement" = "STIGarr",
  "stigma_surface" = "STIGsrf",
  "fruit_ornamentation" = "FRUTorn"
)

# Poaceae ####
quali_poaceae <- c(
  # vegetative traits
  "culm_habit" = "CULMhab",
  "culm_internode_solidity" = "CULMsol",
  "leaf_sheath_margin" = "SHTHmrg",
  "leaf_sheath_indumentum" = "SHTHind",
  "auricle_form" = "AURCfrm",
  "ligule_type" = "LIGLtyp",
  "ligule_margin" = "LIGLmrg",

  # reproductive traits
  "spikelet_arrangement" = "SPKTarr",
  "spikelet_form" = "SPKTfrm",
  "spikelet_indumentum" = "SPKTind",
  "lower_glume_form" = "GLM1frm",
  "lower_glume_apex" = "GLM1apx",
  "lower_glume_indumentum" = "GLM1ind",
  "upper_glume_form" = "GLM2frm",
  "upper_glume_apex" = "GLM2apx",
  "upper_glume_indumentum" = "GLM2ind",
  "callus_form" = "CALSfrm",
  "callus_indumentum" = "CALSind",
  "lemma_form" = "LEMMfrm",
  "lemma_apex" = "LEMMapx",
  "lemma_texture" = "LEMMtxt",
  "lemma_indumentum" = "LEMMind",
  "lemma_color" = "LEMMclr",
  "awn_position" = "AWNNpos",
  "awn_form" = "AWNNfrm",
  "palea_form" = "PALEfrm",
  "palea_keel_form" = "PALEkelfrm",
  "lodicule_form" = "LODIfrm",
  "caryopsis_form" = "CARYfrm",
  "caryopsis_color" = "CARYclr",
  "caryopsis_pericarp_adherence" = "CARYadh",
  "hilum_form" = "HILUfrm"
)

# Polypodiopsida ####
quali_polypodiopsida <- c(
  # vegetative traits
  "rhizome_scale_form" = "RHSCfrm",
  "rhizome_scale_margin" = "RHSCmrg",
  "rhizome_scale_color" = "RHSCclr",
  "frond_division" = "FRONdiv",
  "frond_stalk_indumentum" = "FSTKind",
  "frond_stalk_color" = "FSTKclr",
  "pinna_form" = "PINAfrm",
  "pinna_apex" = "PINAapx",
  "pinna_base" = "PINAbas",
  "pinna_margin" = "PINAmrg",
  "pinnule_form" = "PINLfrm",
  "pinnule_apex" = "PINLapx",
  "pinnule_base" = "PINLbas",
  "pinnule_margin" = "PINLmrg",
  "venation_pattern" = "VENApat",

  # reproductive traits
  "sorus_shape" = "SORUshp",
  "sorus_position" = "SORUpos",
  "indusium_form" = "INDUfrm",
  "indusium_attachment" = "INDUatt",
  "indusium_color" = "INDUclr",
  "sporangium_arrangement" = "SPRGarr",
  "annulus_position" = "ANNLpos",
  "spore_shape" = "SPORshp",
  "spore_ornamentation" = "SPORorn",
  "spore_color" = "SPORclr"
)

# Rubiaceae ####
quali_rubiaceae <- c(
  # vegetative traits
  "stipule_position" = "STIPpos",
  "stipule_fusion" = "STIPfus",
  "colleter_form" = "COLLfrm",
  "colleter_position" = "COLLpos",
  "domatia_form" = "DOMAfrm",
  "domatia_position" = "DOMApos",

  # reproductive traits
  "inflorescence_position" = "INFLpos",
  "floral_morph" = "FMRPtyp",
  "petal_tube_form" = "PETLtubfrm",
  "petal_tube_inner_indumentum" = "PETLtubind",
  "petal_throat_indumentum" = "PETLthrind",
  "stamen_insertion" = "ANDRins",
  "style_indumentum" = "STYLind",
  "stigma_position" = "STIGpos",
  "disc_form" = "DISCfrm",
  "disc_indumentum" = "DISCind",
  "placenta_form" = "PLCNfrm",
  "fruit_crown_form" = "FRUTcrnfrm",
  "pyrene_form" = "PYRNfrm",
  "pyrene_surface" = "PYRNsrf",
  "seed_wing_form" = "SEEDwngfrm"
)

#_______________________________________________________________________________
# Quantitative traits ####

# Default traits ####

quanti_default_traits <- c(
  # vegetative traits
  "internode_length" = "STEMintlng",
  "petiole_length" = "PETIlng",
  "petiole_width" = "PETIwid",
  "stipule_length" = "STIPlng",
  "stipule_width" = "STIPwid",
  "leaf_length" = "LEAFlng",
  "leaf_apex_width" = "LEAFapxwid",
  "leaf_middle_width" = "LEAFmidwid",
  "leaf_base_width" = "LEAFbaswid",
  "lateral_vein_number" = "VENAlatnmb",
  "lateral_vein_angle" = "VENAlatang",

  # reproductive traits
  "peduncle_length" = "PDNClng",
  "inflorescence_length" = "INFLlng",
  "bract_length" = "BRCTlng",
  "bract_width" = "BRCTwid",
  "bract_area" = "BRCTare",
  "bracteole_length" = "BRCLlng",
  "bracteole_width" = "BRCLwid",
  "pedicel_length" = "PDCLlng",
  "pedicel_width" = "PDCLwid",
  "sepal_number" = "SEPLnmb",
  "sepal_length" = "SEPLlng",
  "sepal_width" = "SEPLwid",
  "sepal_lobe_length" = "SEPLloblng",
  "sepal_lobe_width" = "SEPLlobwid",
  "petal_number" = "PETLnmb",
  "petal_length" = "PETLlng",
  "petal_width" = "PETLwid",
  "petal_tube_length" = "PETLtublng",
  "petal_lobe_length" = "PETLloblng",
  "petal_lobe_width" = "PETLlobwid",
  "stamen_number" = "STMNnmb",
  "filament_length" = "FLMTlng",
  "filament_width" = "FLMTwid",
  "anther_length" = "ANTHlng",
  "anther_width" = "ANTHwid",
  "ovary_length" = "OVARlng",
  "ovary_maximum_width" = "OVARwid",
  "locule_number" = "LOCLnmb",
  "ovule_number_per_locule" = "OVULnmb",
  "style_length" = "STYLlng",
  "style_width" = "STYLwid",
  "fruit_length" = "FRUTlng",
  "fruit_width" = "FRUTwid",
  "seed_number" = "SEEDnmb",
  "seed_length" = "SEEDlng",
  "seed_width" = "SEEDwid"
)

# Specific traits ####

# Araceae ####
quanti_araceae <- c(
  # vegetative traits
  "petiole_sheath_length" = "PSHTlng",
  "geniculum_length" = "GENIlng",
  "leaf_lobe_number" = "LEAFlobnmb",

  # reproductive traits
  "peduncle_width" = "PDNCwid",
  "spathe_length" = "SPTHlng",
  "spathe_width" = "SPTHwid",
  "spathe_tube_length" = "SPTHtublng",
  "spathe_blade_length" = "SPTHbldlng",
  "spadix_length" = "SPDXlng",
  "spadix_width" = "SPDXwid",
  "spadix_stipe_length" = "SPDXstplng",
  "spadix_female_zone_length" = "SPDXfemlng",
  "spadix_sterile_zone_length" = "SPDXsterlng",
  "spadix_male_zone_length" = "SPDXmallng",
  "spadix_appendix_length" = "SPDXapplng",
  "flower_number" = "FLWRnmb",
  "perigone_length" = "PRGNlng",
  "stigma_width" = "STIGwid"
)

# Asteraceae ####
quanti_asteraceae <- c(
  # reproductive traits
  "capitulescence_length" = "CPESlng",
  "capitulescence_width" = "CPESwid",
  "capitule_diameter" = "CAPTdia",
  "capitule_number" = "CAPTnmb",
  "phyllary_series_number" = "PHYLsernmb",
  "phyllary_length" = "PHYLlng",
  "phyllary_width" = "PHYLwid",
  "receptacle_diameter" = "RECPdia",
  "receptacle_paleae_length" = "RECPpllng",
  "ray_floret_number" = "RAYFnmb",
  "ray_floret_petal_length" = "RAYFlng",
  "ray_floret_petal_width" = "RAYFwid",
  "disc_floret_number" = "DISFnmb",
  "disc_petal_length" = "DISFlng",
  "disc_petal_tube_length" = "DISFtublng",
  "disc_petal_lobe_length" = "DISFloblng",
  "anther_appendage_length" = "ANTHapplng",
  "style_branch_length" = "STYLbrlng",
  "cypsela_length" = "CYPSlng",
  "cypsela_width" = "CYPSwid",
  "cypsela_carpopodium_length" = "CYPScarlng",
  "pappus_length" = "PAPPlng",
  "pappus_bristle_number" = "PAPPbrsnmb",
  "pappus_series_number" = "PAPPsernmb"
)

# Bromeliaceae ####
quanti_bromeliaceae <- c(
  # vegetative traits
  "leaf_sheath_length" = "SHTHlng",
  "leaf_sheath_width" = "SHTHwid",
  "leaf_blade_length" = "BLADlng",
  "leaf_blade_width" = "BLADwid",
  "leaf_spine_length" = "LSPNlng",
  "leaf_spine_spacing" = "LSPNspc",
  "trichome_scale_diameter" = "TRCHdia",
  "trichome_scale_density" = "TRCHden",

  # reproductive traits
  "scape_length" = "SCAPlng",
  "scape_diameter" = "SCAPdia",
  "scape_bract_number" = "SCBRnmb",
  "scape_bract_length" = "SCBRlng",
  "scape_bract_width" = "SCBRwid",
  "primary_bract_length" = "PBRClng",
  "primary_bract_width" = "PBRCwid",
  "floral_bract_length" = "FBRClng",
  "floral_bract_width" = "FBRCwid",
  "nectary_length" = "NCTRlng",
  "seed_appendage_length" = "SEEDapdlng"
)

# Bryophyta ####
quanti_bryophyta <- c(
  # vegetative traits
  "phyllid_length" = "PHYDlng",
  "phyllid_width" = "PHYDwid",
  "costa_length" = "COSTlng",
  "costa_width_at_base" = "COSTwidbas",
  "phyllid_cell_length" = "PCELlng",
  "phyllid_cell_width" = "PCELwid",
  "alar_cell_number" = "ALARnmb",

  # reproductive traits
  "seta_length" = "SETAlng",
  "seta_width" = "SETAwid",
  "capsule_length" = "CAPSlng",
  "capsule_width" = "CAPSwid",
  "operculum_length" = "OPCLlng",
  "peristome_teeth_number" = "PERInmb",
  "peristome_teeth_length" = "PERIlng",
  "calyptra_length" = "CLPTlng",
  "spore_diameter" = "SPORdia",
  "elater_length" = "ELATlng"
)

# Cactaceae ####
quanti_cactaceae <- c(
  # vegetative traits
  "stem_diameter" = "STEMdia",
  "rib_number" = "RIBSnmb",
  "areole_distance" = "ARLDdst",
  "central_spine_number" = "CSPNnmb",
  "central_spine_length" = "CSPNlng",
  "radial_spine_number" = "RSPNnmb",
  "radial_spine_length" = "RSPNlng",
  "tubercle_length" = "TUBRlng",

  # reproductive traits
  "perigone_length" = "PRGNlng",
  "perigone_diameter" = "PRGNdia",
  "pericarpel_length" = "PRCPlng",
  "pericarpel_diameter" = "PRCPdia"
)

# Leguminosae Caesalpinioideae ####
quanti_caesalpinioideae <- c(
  # vegetative traits
  "distance_between_leaflets" = "LFLTdst",
  "pulvinus_length" = "PLVNlng",
  "pulvinus_width" = "PLVNwid",
  "pulvinulus_length" = "PVNLlng",
  "pulvinulus_width" = "PVNLwid",
  "rachis_length" = "RACHlng",
  "leaflet_number" = "LFLTnmb",
  "leaflet_length" = "LFLTlng",
  "leaflet_width" = "LFLTwid",
  "rachilla_length" = "RCHLlng",
  "distance_between_pinnules" = "PINLdst",
  "pinnule_number" = "PINLnmb",
  "pinnule_length" = "PINLlng",
  "pinnule_width" = "PINLwid"
)

# Eriocaulaceae ####
quanti_eriocaulaceae <- c(
  # vegetative traits
  "scape_length" = "SCAPlng",
  "scape_diameter" = "SCAPdia",
  "scape_rib_number" = "SCAPribnmb",
  "scape_sheath_length" = "SCAPshtlng",
  "leaf_number_per_rosette" = "LEAFnmb",

  # reproductive traits
  "capitule_diameter" = "CAPTdia",
  "capitule_height" = "CAPThgt",
  "involucral_bract_series_number" = "INVLsernmb",
  "involucral_bract_length" = "INVLlng",
  "involucral_bract_width" = "INVLwid",
  "palea_length" = "PALElng",
  "palea_width" = "PALEwid",
  "staminate_flower_number_per_capitulum" = "STFLnmb",
  "staminate_sepal_length" = "STFLsepllng",
  "staminate_petal_length" = "STFLpetllng",
  "pistillate_flower_number_per_capitulum" = "PSFLnmb",
  "pistillate_sepal_length" = "PSFLsepllng",
  "pistillate_petal_length" = "PSFLpetllng",
  "style_branch_number" = "STYLbrnmb",
  "style_branch_length" = "STYLbrlng",
  "seed_appendage_length" = "SDAPlng"
)

# Euphorbiaceae ####
quanti_euphorbiaceae <- c(
  # vegetative traits
  "trichome_ray_number" = "TRCHraynmb",
  "trichome_ray_length" = "TRCHraylng",
  "leaf_gland_number" = "LGLDnmb",
  "leaf_gland_length" = "LGLDlng",
  "leaf_gland_width" = "LGLDwid",
  "spine_number_per_node" = "SPINnmb",
  "spine_length" = "SPINlng",

  # reproductive traits
  "cyathium_length" = "CYATlng",
  "cyathium_width" = "CYATwid",
  "cyathium_involucre_length" = "CYATinvlng",
  "cyathium_involucre_lobe_number" = "CYATlobnmb",
  "cyathium_involucre_lobe_length" = "CYATloblng",
  "cyathium_gland_number" = "CYATglanmb",
  "cyathium_gland_length" = "CYATglalng",
  "cyathium_gland_width" = "CYATglawid",
  "cyathium_gland_appendage_length" = "CYATglaapplng",
  "cyathophyll_length" = "CYPHlng",
  "cyathophyll_width" = "CYPHwid",
  "staminate_flower_number_per_cyathium" = "STFLnmb",
  "staminate_flower_pedicel_length" = "STFLpedlng",
  "pistillate_flower_pedicel_length" = "PSFLpedlng",
  "style_branch_number" = "STYLbrnmb",
  "style_branch_length" = "STYLbrlng",
  "coccus_length" = "COCClng",
  "coccus_width" = "COCCwid",
  "columella_length" = "COLUlng",
  "caruncle_length" = "CARUlng",
  "caruncle_width" = "CARUwid"
)

# Lecythidaceae ####
quanti_lecythidaceae <- c(
  # reproductive traits
  "calyptra_length" = "CALYPlng",
  "calyptra_width" = "CALYPwid",
  "androecial_ring_width" = "ANDRwid",
  "androecial_hood_length" = "HOODlng",
  "androecial_hood_width" = "HOODwid",
  "staminode_number" = "STMDnmb",
  "fertile_stamen_number" = "STMNfernmb",
  "pyxidium_operculum_diameter" = "PYXOdia",
  "pyxidium_operculum_height" = "PYXOhgt",
  "pericarp_thickness" = "PERCthk",
  "aril_length" = "ARILlng",
  "aril_width" = "ARILwid",
  "testa_thickness" = "TESTthk",
  "funicle_length" = "FUNIlng"
)

# Lycopodiopsida ####
quanti_lycopodiopsida <- c(
  # vegetative traits
  "stem_width" = "STEMwid",
  "stem_branching_angle" = "STEMbrnang",
  "stem_branching_number" = "STEMbrnnmb",
  "rhizophore_length" = "RHPHlng",
  "rhizophore_diameter" = "RHPHdia",
  "corm_diameter" = "CORMdia",
  "microphyll_length" = "MPHLlng",
  "microphyll_width" = "MPHLwid",
  "microphyll_number_per_cm" = "MPHLnmb",
  "lateral_microphyll_length" = "MPHLlatlng",
  "lateral_microphyll_width" = "MPHLlatwid",
  "median_microphyll_length" = "MPHLmedlng",
  "median_microphyll_width" = "MPHLmedwid",
  "microphyll_terminal_seta_length" = "MPHLsetlng",
  "ligule_length" = "LIGLlng",

  # reproductive traits
  "strobilus_number" = "STRBnmb",
  "strobilus_length" = "STRBlng",
  "strobilus_width" = "STRBwid",
  "strobilus_peduncle_length" = "STRBpedlng",
  "sporophyll_number_per_strobilus" = "SPHLnmb",
  "sporophyll_length" = "SPHLlng",
  "sporophyll_width" = "SPHLwid",
  "sporangium_length" = "SPRGlng",
  "sporangium_width" = "SPRGwid",
  "velum_coverage" = "VELUcov",
  "spore_diameter" = "SPORdia",
  "megaspore_diameter" = "MSPOdia",
  "megaspore_number_per_sporangium" = "MSPOnmb",
  "microspore_diameter" = "MISPdia"
)

# Moraceae ####
quanti_moraceae <- c(
  # vegetative traits
  "stipule_scar_length" = "STIPscrlng",
  "cystolith_length" = "CYSTlng",

  # reproductive traits
  "receptacle_diameter" = "RECPdia",
  "receptacle_wall_thickness" = "RECPwalthk",
  "ostiole_diameter" = "OSTLdia",
  "ostiole_scale_number" = "OSTLscanmb",
  "staminate_flower_number" = "STFLnmb",
  "staminate_perianth_lobe_number" = "STFLperlobnmb",
  "pistillate_flower_number" = "PSFLnmb",
  "pistillate_perianth_lobe_number" = "PSFLperlobnmb",
  "syncarp_length" = "SYNClng",
  "syncarp_diameter" = "SYNCdia",
  "drupelet_number" = "DRUPnmb",
  "drupelet_length" = "DRUPlng",
  "achenium_length" = "ACHElng",
  "achenium_width" = "ACHEwid"
)

# Orchidaceae ####
quanti_orchidaceae <- c(
  # vegetative traits
  "pseudobulb_length" = "PSBLlng",
  "pseudobulb_width" = "PSBLwid",

  # reproductive traits
  "labellum_length" = "LABLlng",
  "labellum_width" = "LABLwid",
  "column_length" = "COLMlng",
  "column_width" = "COLMwid",
  "anther_cap_length" = "ANTClng",
  "pollinia_number" = "POLNnmb",
  "spur_length" = "SPURlng",
  "spur_width" = "SPURwid"
)

# Leguminosae Papilionoideae ####
quanti_papilionoideae <- c(
  # vegetative traits
  "distance_between_leaflets" = "LFLTdst",
  "pulvinus_length" = "PLVNlng",
  "pulvinus_width" = "PLVNwid",
  "pulvinulus_length" = "PVNLlng",
  "pulvinulus_width" = "PVNLwid",
  "rachis_length" = "RACHlng",
  "leaflet_number" = "LFLTnmb",
  "basal_leaflet_length" = "LFLTbaslng",
  "basal_leaflet_maximum_width" = "LFLTbaswid",
  "basal_leaflet_area" = "LFLTbasare",
  "middle_leaflet_length" = "LFLTmidlng",
  "middle_leaflet_maximum_width" = "LFLTmidwid",
  "middle_leaflet_area" = "LFLTmidare",
  "apex_leaflet_length" = "LFLTapxlng",
  "apex_leaflet_maximum_width" = "LFLTapxwid",
  "apex_leaflet_area" = "LFLTapxare",

  # reproductive traits
  "standard_petal_length" = "STNDlng",
  "standard_petal_width" = "STNDwid",
  "wing_petal_length" = "WINGlng",
  "wing_petal_width" = "WINGwid",
  "keel_petal_length" = "KEELlng",
  "keel_petal_width" = "KEELwid",
  "stipe_length" = "STPElng",
  "stipe_width" = "STPEwid",
  "valve_length" = "VALVlng",
  "valve_width" = "VALVwid",
  "seed_chamber_number" = "SDCHnmb",
  "seed_chamber_length" = "SDCHlng",
  "seed_chamber_width" = "SDCHwid",
  "hilum_length" = "HILUlng",
  "hilum_width" = "HILUwid"
)

# Passifloraceae ####
quanti_passifloraceae <- c(
  # vegetative traits
  "stipule_area" = "STIPare",
  "petiole_nectary_number" = "PTNCnmb",
  "petiole_nectary_length" = "PTNClng",
  "petiole_nectary_width" = "PTNCwid",
  "foliar_nectary_number" = "FLNCnmb",
  "foliar_nectary_length" = "FLNClng",
  "foliar_nectary_width" = "FLNCwid",

  # reproductive traits
  "hypanthium_length" = "HIPTlng",
  "hypanthium_width" = "HIPTwid",
  "coronal_filament_series_number" = "CRFLnmb",
  "coronal_filament_length" = "CRFLlng",
  "coronal_filament_width" = "CRFLwid",
  "operculum_length" = "OPERlng",
  "limen_length" = "LIMNlng",
  "limen_width" = "LIMNwid",
  "androgynophore_length" = "ANGYlng",
  "androgynophore_width" = "ANGYwid"
)

# Pinophyta ####
quanti_pinophyta <- c(
  # vegetative traits
  "twig_diameter" = "TWIGdia",
  "bud_length" = "BUDSlng",
  "needle_length" = "NEEDlng",
  "needle_width" = "NEEDwid",
  "needle_number_per_fascicle" = "FASCnmb",
  "fascicle_sheath_length" = "FASClng",
  "needle_stomatal_line_number" = "NEEDstmnmb",
  "needle_resin_canal_number" = "NEEDrcnmb",

  # reproductive traits
  "pollen_cone_length" = "POLClng",
  "pollen_cone_width" = "POLCwid",
  "seed_cone_length" = "SEEClng",
  "seed_cone_width" = "SEECwid",
  "seed_cone_peduncle_length" = "SEECpedlng",
  "cone_scale_number" = "SCALnmb",
  "cone_scale_length" = "SCALlng",
  "cone_scale_width" = "SCALwid",
  "bract_scale_length" = "BRSClng",
  "bract_scale_width" = "BRSCwid",
  "apophysis_length" = "APOPlng",
  "apophysis_width" = "APOPwid",
  "umbo_length" = "UMBOlng",
  "aril_length" = "ARILlng",
  "seed_body_length" = "SEEDbodlng",
  "seed_body_width" = "SEEDbodwid",
  "seed_wing_length" = "SEEDwnglng",
  "seed_wing_width" = "SEEDwngwid",
  "seed_number_per_scale" = "SEEDscanmb"
)

# Piperaceae ####
quanti_piperaceae <- c(
  # vegetative traits
  "leaf_pellucid_gland_density" = "LGLDden",

  # reproductive traits
  "spike_length" = "SPIKlng",
  "spike_diameter" = "SPIKdia",
  "spike_peduncle_length" = "SPIKpedlng",
  "flower_density_per_cm" = "FLWRden",
  "floral_bract_diameter" = "FBRCdia",
  "floral_bract_stipe_length" = "FBRCstplng"
)

# Poaceae ####
quanti_poaceae <- c(
  # vegetative traits
  "culm_height" = "CULMhgt",
  "culm_diameter" = "CULMdia",
  "culm_internode_length" = "CULMintlng",
  "leaf_sheath_length" = "SHTHlng",
  "auricle_length" = "AURClng",
  "ligule_length" = "LIGLlng",
  "leaf_blade_length" = "BLADlng",
  "leaf_blade_width" = "BLADwid",

  # reproductive traits
  "inflorescence_branch_number" = "INFLbrnmb",
  "spikelet_length" = "SPKTlng",
  "spikelet_width" = "SPKTwid",
  "floret_number" = "FLRTnmb",
  "rachilla_length" = "RCHLlng",
  "lower_glume_length" = "GLM1lng",
  "lower_glume_width" = "GLM1wid",
  "lower_glume_vein_number" = "GLM1nrvnmb",
  "upper_glume_length" = "GLM2lng",
  "upper_glume_width" = "GLM2wid",
  "upper_glume_vein_number" = "GLM2nrvnmb",
  "callus_length" = "CALSlng",
  "lemma_length" = "LEMMlng",
  "lemma_width" = "LEMMwid",
  "lemma_vein_number" = "LEMMnrvnmb",
  "awn_length" = "AWNNlng",
  "palea_length" = "PALElng",
  "palea_width" = "PALEwid",
  "lodicule_length" = "LODIlng",
  "caryopsis_length" = "CARYlng",
  "caryopsis_width" = "CARYwid",
  "hilum_length" = "HILUlng"
)

# polypodiopsida ####
quanti_polypodiopsida <- c(
  # vegetative traits
  "rhizome_diameter" = "RHIZdia",
  "rhizome_scale_length" = "RHSClng",
  "rhizome_scale_width" = "RHSCwid",
  "frond_stalk_length" = "FSTKlng",
  "frond_stalk_width" = "FSTKwid",
  "rachis_length" = "RACHlng",
  "pinna_number" = "PINAnmb",
  "pinna_length" = "PINAlng",
  "pinna_width" = "PINAwid",
  "pinnule_number" = "PINLnmb",
  "pinnule_length" = "PINLlng",
  "pinnule_width" = "PINLwid",

  # reproductive traits
  "sorus_diameter" = "SORUdia",
  "sorus_number" = "SORUnmb",
  "spore_length" = "SPORlng",
  "spore_width" = "SPORwid"
)

# Rubiaceae ####
quanti_rubiaceae <- c(
  # vegetative traits
  "colleter_length" = "COLLlng",
  "domatia_number_per_leaf" = "DOMAnmb",
  "domatia_diameter" = "DOMAdia",

  # reproductive traits
  "sepal_tube_length" = "SEPLtublng",
  "petal_tube_mouth_width" = "PETLtubmouwid",
  "stamen_insertion_height" = "ANDRinshgt",
  "stamen_exsertion_length" = "STMNexslng",
  "style_exsertion_length" = "STYLexslng",
  "stigma_lobe_number" = "STIGlobnmb",
  "stigma_lobe_length" = "STIGloblng",
  "disc_height" = "DISChgt",
  "fruit_crown_length" = "FRUTcrnlng",
  "pyrene_number" = "PYRNnmb",
  "pyrene_length" = "PYRNlng",
  "seed_wing_length" = "SEEDwnglng"
)

#_______________________________________________________________________________

# Lists for specific trait groups
quali_specific_list <- list(
  araceae = quali_araceae,
  asteraceae = quali_asteraceae,
  bromeliaceae = quali_bromeliaceae,
  bryophyta = quali_bryophyta,
  cactaceae = quali_cactaceae,
  caesalpinioideae = quali_caesalpinioideae,
  eriocaulaceae = quali_eriocaulaceae,
  euphorbiaceae = quali_euphorbiaceae,
  lecythidaceae = quali_lecythidaceae,
  lycopodiopsida = quali_lycopodiopsida,
  moraceae = quali_moraceae,
  orchidaceae = quali_orchidaceae,
  papilionoideae = quali_papilionoideae,
  passifloraceae = quali_passifloraceae,
  pinophyta = quali_pinophyta,
  piperaceae = quali_piperaceae,
  poaceae = quali_poaceae,
  polypodiopsida = quali_polypodiopsida,
  rubiaceae = quali_rubiaceae
)

quanti_specific_list <- list(
  araceae = quanti_araceae,
  asteraceae = quanti_asteraceae,
  bromeliaceae = quanti_bromeliaceae,
  bryophyta = quanti_bryophyta,
  cactaceae = quanti_cactaceae,
  caesalpinioideae = quanti_caesalpinioideae,
  eriocaulaceae = quanti_eriocaulaceae,
  euphorbiaceae = quanti_euphorbiaceae,
  lecythidaceae = quanti_lecythidaceae,
  lycopodiopsida = quanti_lycopodiopsida,
  moraceae = quanti_moraceae,
  orchidaceae = quanti_orchidaceae,
  papilionoideae = quanti_papilionoideae,
  passifloraceae = quanti_passifloraceae,
  pinophyta = quanti_pinophyta,
  piperaceae = quanti_piperaceae,
  poaceae = quanti_poaceae,
  polypodiopsida = quanti_polypodiopsida,
  rubiaceae = quanti_rubiaceae
)

#_______________________________________________________________________________

# Argument checks ####

# Resolve `trait_name` (falls back to the first option when omitted)
trait_name <- match.arg(trait_name)

# Validate presets (case-insensitive); fail loudly instead of silently
# returning no traits
if (!is.null(quali_specific)) {
  quali_specific <- match.arg(tolower(quali_specific),
                              choices = names(quali_specific_list))
}

if (!is.null(quanti_specific)) {
  quanti_specific <- match.arg(tolower(quanti_specific),
                               choices = names(quanti_specific_list))
}

# Presets that belong to non-angiosperm lineages.
non_angiosperm_presets <- c("bryophyta", "lycopodiopsida", "pinophyta",
                            "polypodiopsida")

# Auto-disable angiosperm default traits when a non-angiosperm preset is requested
if (!is.null(quali_specific) && quali_specific %in% non_angiosperm_presets) {
  if (quali_default && verbose) {
    message(
      "`quali_specific = \"", quali_specific, "\"` is a non-angiosperm group; ",
      "the angiosperm `quali_default` traits do not apply and have been excluded. ",
      "Set `quali_default = FALSE` explicitly to silence this message."
    )
  }
  quali_default <- FALSE
}

if (!is.null(quanti_specific) && quanti_specific %in% non_angiosperm_presets) {
  if (quanti_default && verbose) {
    message(
      "`quanti_specific = \"", quanti_specific, "\"` is a non-angiosperm group; ",
      "the angiosperm `quanti_default` traits do not apply and have been excluded. ",
      "Set `quanti_default = FALSE` explicitly to silence this message."
    )
  }
  quanti_default <- FALSE
}

#_______________________________________________________________________________

# Assemble trait names and codes ####

quali_traits <- c()
quali_codes <- c()
quanti_traits <- c()
quanti_codes <- c()

# Add default qualitative traits
if (quali_default) {
  quali_traits <- switch(trait_name,
                         "trait full name" = names(quali_default_traits),
                         "trait code" = unname(quali_default_traits),
                         "both" = paste(names(quali_default_traits),
                                        quali_default_traits, sep = "/")
  )
  quali_codes <- unname(quali_default_traits)
}

# Add specific qualitative traits
if (!is.null(quali_specific)) {
  traits_fam <- quali_specific_list[[quali_specific]]
  quali_traits <- c(quali_traits, switch(trait_name,
                                         "trait full name" = names(traits_fam),
                                         "trait code" = unname(traits_fam),
                                         "both" = paste(names(traits_fam),
                                                        traits_fam, sep = "/")
  ))
  quali_codes <- c(quali_codes, unname(traits_fam))
}

# Add default quantitative traits
if (quanti_default) {
  quanti_traits <- switch(trait_name,
                          "trait full name" = names(quanti_default_traits),
                          "trait code" = unname(quanti_default_traits),
                          "both" = paste(names(quanti_default_traits),
                                         quanti_default_traits, sep = "/")
  )
  quanti_codes <- unname(quanti_default_traits)
}

# Add specific quantitative traits
if (!is.null(quanti_specific)) {
  traits_fam <- quanti_specific_list[[quanti_specific]]
  quanti_traits <- c(quanti_traits, switch(trait_name,
                                           "trait full name" = names(traits_fam),
                                           "trait code" = unname(traits_fam),
                                           "both" = paste(names(traits_fam),
                                                          traits_fam, sep = "/")
  ))
  quanti_codes <- c(quanti_codes, unname(traits_fam))
}

# Combine all traits, keeping each label paired with its code
all_traits_raw <- c(quali_traits, quanti_traits)
all_codes_raw <- c(quali_codes, quanti_codes)

keep <- !duplicated(all_traits_raw)
all_traits <- all_traits_raw[keep]
all_codes <- all_codes_raw[keep]

#_______________________________________________________________________________

# Load or create data frame

# Start with an empty df_traits
df_traits <- data.frame(matrix(ncol = length(all_traits), nrow = 0))
colnames(df_traits) <- all_traits

if (is.null(base_df)) {
  base_df <- df_traits
}

# Add any missing trait columns
missing_traits <- setdiff(all_traits, names(base_df))
base_df[missing_traits] <- NA

# Define character prefix order

character_order <- c(
  # vegetative traits — general (angiosperm)
  "STEM","LATX", "CULM", "VELA", "PSBL", "CATA", "TNDR", "STIP", "PLVN",
  "SHTH", "AURC", "LIGL", "RIBS", "ARLW", "ARLD", "SPIN", "CSPN", "RSPN",
  "GLOC", "TUBR", "PVNL", "PETI", "PSHT", "GENI", "PTNC", "BLAD", "LEAF",
  "LSPN", "CYST", "TRCH", "LGLD", "COLL", "DOMA", "FLNC", "RACH", "LFLT",
  "RCHL", "PINL", "VENA",

  # vegetative traits — moss-specific
  "RHIZ", "PHYD", "COST", "PCEL", "ALAR",

  # vegetative traits — conifer-specific
  "CRWN", "BARK", "BRNC", "TWIG", "BUDS", "FASC", "NEED",

  # vegetative traits — fern-specific
  "RHSC", "FRON", "FSTK", "PINA", "PINL",

  # vegetative traits — lycophyte-specific
  "RHPH", "CORM", "MPHL",

  # reproductive traits — general (angiosperm)
  "PDNC", "INFL", "SPIK", "CYAT", "CYPH", "OSTL", "FMRP", "SPTH", "SPDX",
  "FLWR", "STFL", "PSFL", "GLFL", "PRGN", "SPKT", "GLM1", "GLM2", "FLRT",
  "CALS", "LEMM", "AWNN", "PALE", "LODI", "CPES", "CAPT", "INVL", "PHYL",
  "RECP", "BRCT", "SCAP", "SCBR", "PBRC", "FBRC", "BRCL", "PDCL", "HIPT",
  "SEPL", "RAYF", "DISF", "PAPP", "PETL", "LABL", "SPUR", "STND", "WING",
  "KEEL", "CRFL", "OPER", "LIMN", "ANGY", "STPE", "COLM", "PRCP", "OVAR",
  "NCTR", "DISC", "PLCN", "LOCL", "OVUL", "STIG", "STYL", "ANDR", "CALY",
  "HOOD", "STMN", "STMD", "FLMT", "ANTH", "ANTC", "POLN", "FRUT", "CYPS",
  "ACHE", "PYXO", "PERC", "SYNC", "DRUP", "COCC", "COLU", "PYRN", "FSCL",
  "CARY", "VALV", "SDCH", "CARU", "SEED", "TEST", "FUNI", "HILU",

  # reproductive traits — moss-specific
  "SETA", "CAPS", "OPCL", "PERI", "CLPT", "SPOR", "ELAT",

  # reproductive traits — conifer-specific
  "POLC", "SEEC", "SCAL", "BRSC", "APOP", "UMBO", "ARIL",

  # reproductive traits — fern-specific
  "SORU", "INDU", "SPRG", "ANNL",

  # reproductive traits — lycophyte-specific
  "STRB", "SPHL", "MSPO", "MISP", "VELU"
)

# Warn if any trait code uses a prefix that is not in `character_order`
unknown_prefix <- setdiff(unique(substr(all_codes, 1, 4)), character_order)
if (length(unknown_prefix) > 0 && verbose) {
  message(
    "Trait prefixes not listed in `character_order` (columns placed last): ",
    paste(unknown_prefix, collapse = ", ")
  )
}

# Determine correct column order: non-trait columns first, then traits by prefix
trait_order <- order(match(substr(all_codes, 1, 4), character_order))

ordered_cols <- c(
  setdiff(names(base_df), all_traits),
  all_traits[trait_order]
)

df_traits <- base_df[, ordered_cols, drop = FALSE]

# Define output file name/path if not supplied by the user
if (save && is.null(file_name)) {
  if (!dir.exists("output_data/")) {
    dir.create("output_data/")
  }
  file_name <- file.path("output_data/",
                         paste0("trait_data_",
                                format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"))
}

# Export to Excel

if (save) {
  openxlsx::write.xlsx(df_traits, file_name)
  if (verbose) {
    message("Trait DataFrame exported to: ", file_name)
  }
}

return(df_traits)
}
