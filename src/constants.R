

MOD_TYPE_PRETTY <- c(
  a = "m6",
  m = "m5",
  `21839` = "m4"
)
MOD_TYPE_SINGLE <- c(
  "6mA" = "a",
  "5mC" = "m",
  "4mC" = "21839"
)


PLOT_COLORS <- c(
  light_blue="#B9CFF0",
  dark_blue="#404080",
  light_green="#9FE2BF",
  dark_green="#69B3A2",
  medium_blue="#5E5EAF"
)


MOD_CODE_TO_PRETTY <- c(
  "m" = "5mC",
  "a" = "6mA",
  "21839" = "4mC"
)
MOD_PRETTY_TO_CODE <- setNames(names(MOD_CODE_TO_PRETTY), MOD_CODE_TO_PRETTY)