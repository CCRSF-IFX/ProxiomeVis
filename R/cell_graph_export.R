cell_graph_appearance <- function(background_color = "#D3D3D3", highlight_color = "#C00000",
    background_opacity = 64 / 255, highlight_opacity = 1, background_radius = 0.003, highlight_radius = 0.003) {
  colors <- c(background_color, highlight_color)
  if (!all(vapply(list(background_color, highlight_color), function(x) is.character(x) && length(x) == 1L, logical(1))) ||
      anyNA(colors) || !all(grepl("^#[0-9a-fA-F]{6}$", colors)))
    cell_graph_error("Use a six-digit hex color, such as #D3D3D3.")
  values <- list(background_opacity, highlight_opacity, background_radius, highlight_radius)
  if (!all(vapply(values, function(x) is.numeric(x) && length(x) == 1L && is.finite(x), logical(1))) ||
      any(c(background_opacity, highlight_opacity) < 0 | c(background_opacity, highlight_opacity) > 1) ||
      any(c(background_radius, highlight_radius) < 0.0001 | c(background_radius, highlight_radius) > 0.05))
    cell_graph_error("Opacity must be 0–1; dot radius must be 0.0001–0.05.")
  list(background_color = toupper(background_color), highlight_color = toupper(highlight_color),
       background_opacity = background_opacity, highlight_opacity = highlight_opacity,
       background_radius = background_radius, highlight_radius = highlight_radius)
}

write_cell_graph_glb <- function(asset, marker, glb, appearance = cell_graph_appearance()) {
  validate_cell_graph_asset(asset)
  appearance <- do.call(cell_graph_appearance, appearance)
  nodes <- asset$nodes
  if (!is.character(marker) || length(marker) != 1L || is.na(marker) ||
      !marker %in% unlist(nodes$markers, use.names = FALSE))
    cell_graph_error("Choose one highlight marker present in the displayed cell.")
  if (nrow(nodes) > 250000L) cell_graph_error("GLB export supports up to 250,000 nodes.")
  if (any(abs(as.matrix(nodes[c("x", "y", "z")])) > 3e38))
    cell_graph_error("Coordinates exceed the GLB float32 range.")
  id <- asset$cell_id
  red <- vapply(nodes$markers, function(m) marker %in% m, logical(1))
  span <- max(vapply(nodes[c("x", "y", "z")], function(x) diff(range(x)), numeric(1)))
  if (span == 0) span <- 1
  v <- rbind(c(1,1,1),c(1,-1,-1),c(-1,1,-1),c(-1,-1,1))/sqrt(3)
  f <- rbind(c(1L,2L,3L),c(1L,4L,2L),c(1L,3L,4L),c(2L,4L,3L))
  chunks <- views <- accessors <- primitives <- list()
  offset <- 0L
  for (h in c(FALSE,TRUE)) {
    ids <- which(red==h)
    if (!length(ids)) next
    n <- length(ids)
    xyz <- as.matrix(nodes[ids,c('x','y','z')])
    radius <- span * if (h) appearance$highlight_radius else appearance$background_radius
    positions <- xyz[rep(seq_len(n),each=4),]+v[rep(1:4,n),]*radius
    position_bytes <- writeBin(as.double(t(positions)),raw(),size=4L,endian='little')
    # Bounds describe the float32 coordinates actually stored in the file.
    positions <- matrix(readBin(position_bytes,'double',n=n*12L,size=4L,endian='little'),ncol=3,byrow=TRUE)
    indices <- rep(as.integer(t(f))-1L,n)+rep((seq_len(n)-1L)*4L,each=12L)
    stopifnot(min(indices)==0L,max(indices)==n*4L-1L,length(indices)==n*12L)
    index_bytes <- writeBin(as.integer(indices),raw(),size=4L,endian='little')
    k <- length(views)
    views[[k+1L]] <- list(buffer=0L,byteOffset=offset,byteLength=length(position_bytes),target=34962L)
    views[[k+2L]] <- list(buffer=0L,byteOffset=offset+length(position_bytes),byteLength=length(index_bytes),target=34963L)
    accessors[[k+1L]] <- list(bufferView=k,componentType=5126L,count=n*4L,type='VEC3',
      min=unname(apply(positions,2,min)),max=unname(apply(positions,2,max)))
    accessors[[k+2L]] <- list(bufferView=k+1L,componentType=5125L,count=n*12L,type='SCALAR')
    primitives[[length(primitives)+1L]] <- list(attributes=list(POSITION=k),indices=k+1L,material=as.integer(h),mode=4L)
    chunks[[k+1L]] <- position_bytes
    chunks[[k+2L]] <- index_bytes
    offset <- offset+length(position_bytes)+length(index_bytes)
  }
  # Material factors are linear RGB; convert the threejs sRGB hex colors.
  linear <- function(x) ifelse(x<=0.04045,x/12.92,((x+0.055)/1.055)^2.4)
  material <- function(name,rgb,alpha) list(name=name,
    pbrMetallicRoughness=list(baseColorFactor=c(linear(rgb),alpha),metallicFactor=0,roughnessFactor=1),
    alphaMode=if(alpha<1) 'BLEND' else 'OPAQUE')
  document <- list(asset=list(version='2.0',generator='ProxiomeVis'),scene=0L,
    scenes=list(list(nodes=list(0L))),nodes=list(list(mesh=0L,name=id)),
    meshes=list(list(primitives=primitives)),
    materials=list(material('Other nodes',as.numeric(grDevices::col2rgb(appearance$background_color))/255,appearance$background_opacity),
      material(marker,as.numeric(grDevices::col2rgb(appearance$highlight_color))/255,appearance$highlight_opacity)),
    buffers=list(list(byteLength=offset)),bufferViews=views,accessors=accessors,
    extras=list(cell_id=id,highlight_marker=marker, appearance=appearance,node_count=nrow(nodes),
      interpretation='Abstract graph layout, not measured physical positions.'))
  json <- charToRaw(enc2utf8(as.character(jsonlite::toJSON(document,auto_unbox=TRUE,digits=NA))))
  json <- c(json,rep(charToRaw(' '),(4L-length(json)%%4L)%%4L))
  con <- file(glb,'wb')
  tryCatch({
    writeBin(as.integer(c(0x46546C67,2L,28L+length(json)+offset,length(json),0x4E4F534A)),con,size=4L,endian='little')
    writeBin(json,con)
    writeBin(as.integer(c(offset,0x004E4942)),con,size=4L,endian='little')
    for (chunk in chunks) writeBin(chunk,con)
  },finally=close(con))
  invisible(glb)
}
