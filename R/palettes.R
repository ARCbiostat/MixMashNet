#' Default layer palettes for MixMashNet
#'
#' Internal palette bank used to color layers and communities.
#' Each palette contains 20 colors.
#'
#' @keywords internal
#' @noRd
.mixmashnet_palettes <- list(

  palette1 = c(
    "#38E5F4", "#A3AF0E", "#E9192D", "#7D4DB6", "#D95D7D",
    "#DFA357", "#43B883", "#2792ED", "#D448E8", "#ED6145",
    "#C52059", "#830860", "#631EBE", "#3194B9", "#E91898",
    "#D09537", "#F770A8", "#D67BF6", "#949FFD", "#52DB9E"
  ),

  palette2 = c(
    "#C068F0", "#BA378A", "#D86A70", "#4658B7", "#9ACC57",
    "#D740AB", "#FA8E4B", "#656E27", "#288DA1", "#459957",
    "#FBB0B2", "#AD8324", "#225C5E", "#CB5339", "#1F6F7A",
    "#34BDF3", "#5586EA", "#744ABF", "#E4539B", "#924513"
  ),

  palette3 = c(
    "#87ACC1", "#8770A0", "#754447", "#9F6A82", "#26331D",
    "#575635", "#8E7955", "#C29D87", "#44222E", "#5C462D",
    "#41271B", "#57704F", "#19323B", "#B99ABB", "#355C5D",
    "#5B8777", "#92B08F", "#392543", "#514D7C", "#647DA5"
  ),

  palette4 = c(
    "#5F00F1", "#0073A1", "#008A7F", "#4F9A00", "#B68900",
    "#FF6640", "#E63946", "#C0A8FF", "#95CDFF", "#00F2F4",
    "#98FFA0", "#FFF1BA", "#8B007B", "#870015", "#503300",
    "#263400", "#002A21", "#001A23", "#112E45", "#3A5157"
  ),

  palette5 = c(
    "#798100", "#A58E00", "#D59600", "#FF9F4C", "#FFBDAA",
    "#008A94", "#00A595", "#00C182", "#54DB00", "#C6DF00",
    "#8941FF", "#6D7FFF", "#61A5FF", "#6EC2FF", "#91DAFF",
    "#E3002F", "#FF1F69", "#FF659B", "#FF91BC", "#FFB6D6"
  ),

  palette6 = c(
    "#9976AE", "#7E566B", "#987878", "#4B526D", "#9EB08B",
    "#996585", "#CB8E6C", "#585C43", "#4E7781", "#5D7F62",
    "#D9A8A8", "#8C7852", "#354748", "#916257", "#3E5F68",
    "#6AA7C5", "#647BA7", "#5B5373", "#A37085", "#694634"
  ),

  palette7 = c(
    "#5A0339", "#B21554", "#FA8283", "#A73916", "#5AEDF6",
    "#6DE251", "#E3A426", "#66560C", "#828516", "#8BBC29",
    "#093C2E", "#11645C", "#189096", "#EC692F", "#9F1199",
    "#082D65", "#3E3AC9", "#855AF4", "#C489FA", "#F2C1F7"
  ),

  palette8 = c(
    "#660031", "#006A67", "#559900", "#AEB900", "#E0E500",
    "#BD8500", "#A93E00", "#FF2A4F", "#FF94AD", "#FFDAE6",
    "#005543", "#00854D", "#00B913", "#8CE000", "#404C00",
    "#0046B7", "#007ECD", "#00B2ED", "#59E4FF", "#00C5BE"
  ),

  palette9 = c(
    "#3A006C", "#7B00A0", "#C800BE", "#FF5BBE", "#A20040",
    "#EA0048", "#FFCDC2", "#933F00", "#BC6F00", "#E1A100",
    "#FED800", "#8D9000", "#B5F300", "#00D7B2", "#79FFF7",
    "#006176", "#0078AA", "#389BFF", "#ABC2FF", "#EFEFFF"
  ),

  palette10 = c(
    "#D0816A", "#4D3F1E", "#527245", "#5E9F98", "#9AC2DC",
    "#91C485", "#265350", "#46789A", "#9A93CD", "#DEBAD2",
    "#163144", "#564789", "#A8618E", "#D29489", "#E2CFAF",
    "#542441", "#884A3C", "#9C8351", "#9DBC89", "#C1E6E0"
  )
)


#' @keywords internal
#' @noRd
.get_layer_palette <- function(i, palette_bank = .mixmashnet_palettes) {
  k <- length(palette_bank)
  if (k < 1L) stop(".get_layer_palette(): empty palette bank.", call. = FALSE)
  palette_bank[[ ((i - 1L) %% k) + 1L ]]
}
