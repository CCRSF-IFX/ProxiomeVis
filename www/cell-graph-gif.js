// Resolve locally, including when Shiny is hosted under a URL prefix.
const proxiomeGifWorkerURL = new URL('cell-graph-gif-worker.js', document.currentScript.src);

window.proxiomeCellGraphGif = function(el, graph, data) {
  const control = suffix => document.getElementById(el.id + suffix);
  const button = control('_gif'), cancel = control('_gif_cancel');
  const progress = control('_gif_progress'), status = control('_gif_status');
  const inputs = ['_gif_seconds', '_gif_fps', '_gif_pixels'].map(control);
  button.disabled = false;
  status.textContent = 'Ready to export the displayed cell.';
  button.onclick = async () => {
    const [seconds, fps, size] = inputs.map(input => Number(input.value));
    const frames = seconds * fps;
    if (![seconds, fps, size].every(Number.isInteger) || seconds < 1 || seconds > 20 ||
        fps < 1 || fps > 30 || size < 128 || size > 1024 ||
        frames > 300 || frames * size * size > 120000000) {
      status.textContent = 'Use whole numbers: 1–20 seconds, 1–30 FPS, 128–1024 pixels. Limit: 300 frames and 120 million total pixels; reduce duration, FPS, or size.';
      return;
    }
    const display = document.getElementById(el.id.replace(/plot$/, 'display'));
    const previousDisabled = display.disabled;
    const appearance = document.getElementById(el.id.replace(/plot$/, 'apply_appearance'));
    const appearanceDisabled = appearance?.disabled;
    if (appearance) appearance.disabled = true;
    button.disabled = display.disabled = true;
    inputs.forEach(input => { input.disabled = true; });
    cancel.hidden = progress.hidden = false;
    progress.max = frames;
    progress.value = 0;
    status.textContent = 'Preparing GIF…';
    let worker, renderer, scene, rejectFrame, cancelled = false;
    const stop = () => {
      cancelled = true;
      worker?.terminate();
      rejectFrame?.(new Error('Export cancelled.'));
    };
    el.proxiomeGifCancel = cancel.onclick = stop;
    try {
      worker = new Worker(proxiomeGifWorkerURL, {type: 'module'});
      const send = (message, transfers = []) => new Promise((resolve, reject) => {
        const timer = setTimeout(() => reject(new Error('GIF encoder timed out. Try a smaller export.')), 60000);
        rejectFrame = error => { clearTimeout(timer); reject(error); };
        worker.onerror = () => rejectFrame(new Error('GIF encoder could not load or failed. Reload the app and try again.'));
        worker.onmessage = ({data}) => {
          clearTimeout(timer);
          if (data.error) reject(new Error(data.error)); else resolve(data);
        };
        worker.postMessage(message, transfers);
      });
      // Independent resources let cleanup release the export without touching the live graph.
      scene = new THREE.Scene();
      scene.add(graph.pointgroup.clone());
      scene.traverse(object => {
        if (object.geometry) object.geometry = object.geometry.clone();
        if (object.material) object.material = object.material.clone();
      });
      const camera = graph.camera.clone();
      const target = graph.controls.target.clone(), up = camera.up.clone().normalize();
      const offset = camera.position.clone().sub(target);
      camera.aspect = 1;
      camera.updateProjectionMatrix();
      renderer = new THREE.WebGLRenderer({antialias: true});
      renderer.setPixelRatio(1);
      renderer.setSize(size, size, false);
      renderer.setClearColor(0xffffff, 1);
      const canvas = document.createElement('canvas');
      canvas.width = canvas.height = size;
      const context = canvas.getContext('2d', {willReadFrequently: true});
      for (let frame = 0; frame < frames; frame++) {
        if (cancelled) throw new Error('Export cancelled.');
        camera.position.copy(offset).applyAxisAngle(up, 2 * Math.PI * frame / frames).add(target);
        camera.lookAt(target);
        renderer.render(scene, camera);
        context.drawImage(renderer.domElement, 0, 0);
        const rgba = context.getImageData(0, 0, size, size).data;
        // GIF delays are centiseconds; distribute rounding to preserve total duration.
        const delay = 10 * (Math.round((frame + 1) * 100 / fps) - Math.round(frame * 100 / fps));
        await send({rgba: rgba.buffer, size, delay}, [rgba.buffer]);
        progress.value = frame + 1;
        status.textContent = 'Encoding frame ' + (frame + 1) + ' / ' + frames + '…';
      }
      const {bytes} = await send({finish: true});
      if (cancelled) throw new Error('Export cancelled.');
      const url = URL.createObjectURL(new Blob([bytes], {type: 'image/gif'}));
      const link = document.createElement('a');
      link.download = 'proxiomevis-' + String(data.main).replace(/[^a-zA-Z0-9_-]/g, '-') + '-rotating.gif';
      link.href = url;
      link.click();
      setTimeout(() => URL.revokeObjectURL(url), 60000);
      status.textContent = 'GIF downloaded: ' + seconds + ' s, ' + fps + ' FPS, ' + size + ' × ' + size + ' pixels.';
    } catch (error) {
      status.textContent = cancelled ? 'Export cancelled.' : 'GIF export failed: ' + error.message;
    } finally {
      worker?.terminate();
      scene?.traverse(object => {
        object.geometry?.dispose();
        object.material?.uniforms?.texture1?.value?.dispose();
        object.material?.dispose();
      });
      renderer?.dispose();
      renderer?.forceContextLoss();
      if (el.proxiomeGifCancel === stop) el.proxiomeGifCancel = null;
      button.disabled = !el.isConnected;
      if (appearance) appearance.disabled = appearanceDisabled;
      display.disabled = previousDisabled;
      inputs.forEach(input => { input.disabled = false; });
      cancel.hidden = progress.hidden = true;
    }
  };
};
