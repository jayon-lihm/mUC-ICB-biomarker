
############ graphic setting #############
col = c(frameshift = "red", missense = "dodgerblue", stop_lost="brown", stop_gain="black", other="moccasin")

alter_fun = list(
  background = alter_graphic("rect", fill = "#CCCCCC"),
  frameshift = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["frameshift"], col = NA)),
  missense = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["missense"], col = NA)),
  stop_lost = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["stop_lost"], col = NA)),
  stop_gain = function(x, y, w, h) 
    grid.segments(x - w*0.9*0.5, y - h*0.9*0.5, x + w*0.9*0.5, y + h*0.9*0.5, gp = gpar(lwd = 2)),
  other = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["other"], col = NA)) )


col_response <- c("CR"="blue", "PR" = "green",
                  "SD"="brown", "PD" = "red",
                  "NE"= "grey")

################ graphic setting END #############

