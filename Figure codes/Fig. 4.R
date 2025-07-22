#---------------------------------------
# Fig. 4
#---------------------------------------
pacman::p_load(tidyverse, readxl, networkD3, htmlwidgets)

# data
data <- read_excel("/Source Data/Source Data Fig. 4.xlsx")

# node and color
all_nodes <- unique(c(data$Exposure, data$Mediator, data$Outcome))
nodes <- data.frame(name = all_nodes)

color_map <- c(
  "LHR"="#ec7686", "Height"="#61d6f7", "HDL"="#8052d8", "APOA"="#c6c5e7",
  "HbA1c"="#d192fb", "AST/ALT"="#fdb9ba", "TG"="#0a086d", "UA"="#bedf50",
  "CRE"="#47bca9", "GLU"="#adb9eb", "NEUT"="#5640a7", "APOB"="#a3e2b3",
  "CRP"="#406bde", "CYS"="#c3f887", "ALT"="#ffdc4a", "WBC"="#c34779",
  "TBIL"="#831962", "TC"="#d543ae", "CAD"="#7595bb", "MI"="#f0b042",
  "T2D"="#1ce8bf", "HT"="#0eb1ee", "HbA1c."="#d192fb", "AST/ALT."="#fdb9ba",
  "TG."="#0a086d", "APOB."="#a3e2b3", "ALT."="#ffdc4a"
)
nodes$color <- color_map[nodes$name]

# link
links <- bind_rows(
  data %>% select(source = Exposure, target = Mediator, value = MP),
  data %>% select(source = Mediator, target = Outcome, value = MP)
) %>%
  group_by(source, target) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  mutate(
    IDsource = match(source, nodes$name) - 1,
    IDtarget = match(target, nodes$name) - 1
  )

colourScale <- sprintf(
  'd3.scaleOrdinal().domain(["%s"]).range(["%s"])',
  paste(nodes$name, collapse = '","'),
  paste(nodes$color, collapse = '","')
)

# sankey plot
sankey <- sankeyNetwork(
  Links = links, Nodes = nodes,
  Source = "IDsource", Target = "IDtarget", Value = "value",
  NodeID = "name", units = "T", fontSize = 12,
  fontFamily = "Arial", nodeWidth = 15, nodePadding = 20,
  width = 1200, height = 800, iterations = 0,
  colourScale = colourScale
)

p <- htmlwidgets::onRender(sankey, '
function(el, x) {
  var svg = d3.select(el).select("svg");

  function validID(name) {
    return name ? name.replace(/[^a-zA-Z0-9]/g, "_") : "unknown";
  }

  svg.selectAll(".link").each(function(d) {
    var gid = "grad-" + validID(d.source.name) + "-" + validID(d.target.name);
    var grad = svg.append("defs").append("linearGradient")
      .attr("id", gid)
      .attr("gradientUnits", "userSpaceOnUse")
      .attr("x1", d.source.x + d.source.dx / 2)
      .attr("x2", d.target.x + d.target.dx / 2)
      .attr("y1", d.source.y + d.sy + d.dy / 2)
      .attr("y2", d.target.y + d.ty + d.dy / 2);

    var c1 = d3.select(el).selectAll(".node").filter(n => n.name === d.source.name).select("rect").style("fill");
    var c2 = d3.select(el).selectAll(".node").filter(n => n.name === d.target.name).select("rect").style("fill");

    grad.append("stop").attr("offset", "0%").attr("stop-color", c1);
    grad.append("stop").attr("offset", "100%").attr("stop-color", c2);

    d3.select(this).style("stroke", "url(#" + gid + ")").style("stroke-opacity", 0.3);
  });

  svg.selectAll(".node rect").style("stroke", "none");

  svg.selectAll(".node text")
    .filter(d => ["HDL","APOA","AST/ALT","UA","TG","HbA1c","APOB","CYS","ALT","GLU","NEUT",
                  "AST/ALT.","TG.","WBC","CRP","ALT.","APOB.","TBIL","HbA1c."].includes(d.name))
    .attr("x", -6)
    .style("text-anchor", "end");
}
')

# save
saveWidget(p, "Fig. 4.html", selfcontained = TRUE)
