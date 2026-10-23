Shiny.addCustomMessageHandler('proxiome-cell-graph-reset', ({id}) => {
  const el = document.getElementById(id);
  if (!el) return;
  el.proxiomeCell = el.proxiomeCamera = null;
  el.proxiomeGifCancel?.();
  const graph = HTMLWidgets.getInstance(el)?.widget;
  if (graph) graph.idle = true;
  el.proxiomeDispose?.();
  el.proxiomeDispose = () => {};
  el.querySelectorAll('canvas').forEach(canvas => canvas.remove());
  for (const suffix of ['_reset', '_png', '_gif']) {
    const button = document.getElementById(id + suffix);
    button.disabled = true;
    button.onclick = null;
  }
  for (const suffix of ['model_glb', 'nodes_csv', 'provenance_json']) {
    const link = document.getElementById(id.replace(/plot$/, suffix));
    link.inert = true;
    link.classList.add('disabled');
    link.setAttribute('aria-disabled', 'true');
  }
  document.getElementById(id + '_gif_status').textContent = 'Display a cell before exporting.';
});

window.proxiomeCellGraph = function(el, data) {
  const graph = HTMLWidgets.getInstance(el).widget;
  if (!el.proxiomeInitialized) {
    el.proxiomeInitialized = true;
    const initialize = graph.init;
    // CRAN threejs 0.3.4 recreates WebGL resources on Display without disposing them.
    graph.init = function(...args) {
      el.proxiomeGifCancel?.();
      el.proxiomeCamera = {cell: el.proxiomeCell, position: graph.camera.position.clone(),
        up: graph.camera.up.clone(), target: graph.controls.target.clone()};
      el.proxiomeDispose();
      graph.idle = true;
      initialize.apply(graph, args);
    };
    const observer = new ResizeObserver(() => {
      if (el.clientWidth && el.clientHeight) graph.controls.handleResize();
    });
    observer.observe(el);
    $(el).one('shiny:unbound.proxiomeGraph', () => {
      el.proxiomeGifCancel?.();
      graph.idle = true;
      observer.disconnect();
      el.proxiomeDispose();
    });
  }
  const renderer = graph.renderer;
  const controls = graph.controls;
  // Residual damping otherwise moves the camera again immediately after Reset.
  controls.staticMoving = true;
  // Trackball wheel events do not wake the package's idle render loop.
  controls.addEventListener('start', () => {
    if (graph.idle) requestAnimationFrame(() => graph.animate());
  });
  const groups = [graph.pointgroup, graph.linegroup];
  el.proxiomeDispose = () => {
    controls.dispose();
    groups.forEach(group => group.traverse(object => {
      object.geometry?.dispose();
      if (object.material) {
        object.material.uniforms?.texture1?.value?.dispose();
        object.material.dispose();
      }
    }));
    renderer.dispose();
    renderer.forceContextLoss();
  };
  const previous = el.proxiomeCamera;
  if (previous?.cell === data.main) {
    graph.camera.position.copy(previous.position);
    graph.camera.up.copy(previous.up);
    controls.target.copy(previous.target);
    controls.update();
    graph.animate();
  }
  el.proxiomeCamera = null;
  el.proxiomeCell = data.main;
  for (const suffix of ['_reset', '_png']) document.getElementById(el.id + suffix).disabled = false;
  for (const suffix of ['model_glb', 'nodes_csv', 'provenance_json']) {
    const link = document.getElementById(el.id.replace(/plot$/, suffix));
    link.inert = false;
    link.classList.remove('disabled');
    link.setAttribute('aria-disabled', 'false');
  }
  window.proxiomeCellGraphGif(el, graph, data);
  document.getElementById(el.id + '_reset').onclick = () => {
    controls.reset();
    graph.animate();
  };
  document.getElementById(el.id + '_png').onclick = () => {
    graph.animate();
    const link = document.createElement('a');
    link.download = 'proxiomevis-cell-3d.png';
    link.href = renderer.domElement.toDataURL('image/png');
    link.click();
  };
};
