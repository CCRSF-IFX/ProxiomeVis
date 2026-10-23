// Prepared hybrid demo on :3838 and Chromium with CDP on :9223 (WebGL enabled).
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';

const debugUrl = process.env.PROXIOME_CDP_URL || 'http://127.0.0.1:9223';
const appUrl = process.env.PROXIOME_SMOKE_URL || 'http://127.0.0.1:3838';
const tabs = await (await fetch(`${debugUrl}/json`)).json();
const socket = new WebSocket(tabs.find(t => t.type === 'page').webSocketDebuggerUrl);
await new Promise(resolve => socket.addEventListener('open', resolve, { once: true }));
let serial = 0;
const pending = new Map();
socket.onmessage = ({ data }) => {
  const message = JSON.parse(data);
  pending.get(message.id)?.(message);
  pending.delete(message.id);
};
const call = (method, params = {}) => new Promise((resolve, reject) => {
  const id = ++serial;
  const timeout = setTimeout(() => reject(new Error(`Timed out: ${method}`)), 30000);
  pending.set(id, message => {
    clearTimeout(timeout);
    if (message.error) reject(new Error(JSON.stringify(message.error)));
    else resolve(message.result);
  });
  socket.send(JSON.stringify({ id, method, params }));
});
const evaluate = async expression => {
  const value = await call('Runtime.evaluate', { expression, returnByValue: true, awaitPromise: true });
  if (value.exceptionDetails) throw new Error(JSON.stringify(value.exceptionDetails));
  return value.result.value;
};
const delay = ms => new Promise(resolve => setTimeout(resolve, ms));
async function waitFor(expression) {
  for (let n = 0; n < 180; n++) {
    if (await evaluate(expression)) return;
    await delay(250);
  }
  throw new Error(`Condition not reached: ${expression}\n${await evaluate("document.body.innerText.slice(-5000)")}`);
}
async function settled() {
  await delay(400);
  await waitFor("!document.documentElement.classList.contains('shiny-busy')");
  assert.deepEqual(await evaluate("[...document.querySelectorAll('.shiny-output-error:not(.shiny-output-error-validation)')].filter(x => x.offsetParent).map(x => x.textContent)"), []);
}
async function tab(label) {
  await evaluate(`[...document.querySelectorAll('a[data-bs-toggle=tab]')].find(x => x.innerText === ${JSON.stringify(label)}).click()`);
  await settled();
}
async function select(id, value) {
  await evaluate(`document.getElementById(${JSON.stringify(id)}).selectize.setValue(${JSON.stringify(value)})`);
  await settled();
  await waitFor(`JSON.stringify(Shiny.shinyapp.$inputValues[${JSON.stringify(id)}]) === JSON.stringify(${JSON.stringify(value)})`);
}
async function download(id) {
  await waitFor(`document.getElementById(${JSON.stringify(id)})?.getAttribute('href')?.includes('download')`);
  return evaluate(`fetch(document.getElementById(${JSON.stringify(id)}).href).then(async r => { if (!r.ok) throw new Error(r.status); return await r.text(); })`);
}
async function run() {
  await evaluate("document.getElementById('report-analysis-run').click()");
  await waitFor("document.getElementById('report-analysis-run').disabled");
  // A live worker must not prevent navigation in the same session.
  await tab('Overview / QC');
  await tab('New comparison');
  await waitFor("document.getElementById('report-analysis-status').textContent.startsWith('Completed.')");
  await settled();
}
async function screenshot(name) {
  await mkdir('output', { recursive: true });
  const shot = await call('Page.captureScreenshot', { format: 'png' });
  await writeFile(`output/${name}.png`, Buffer.from(shot.data, 'base64'));
}
async function checkControlContrast() {
  const failures = await evaluate(`(() => {
    const rgb = value => value.match(/[\\d.]+/g).map(Number);
    const luminance = c => c.slice(0,3).map(v => v/255).map(v => v <= .04045 ? v/12.92 : ((v+.055)/1.055)**2.4)
      .reduce((sum,v,i) => sum + v*[.2126,.7152,.0722][i],0);
    return [...document.querySelectorAll('.btn,.alert')].filter(e => e.offsetParent && !e.disabled).flatMap(e => {
      let bg = [255,255,255];
      for(let parent=e;parent;parent=parent.parentElement) {
        const candidate = rgb(getComputedStyle(parent).backgroundColor);
        if(candidate.length===3 || candidate[3]===1) { bg=candidate; break; }
      }
      const fg=rgb(getComputedStyle(e).color), a=luminance(fg), b=luminance(bg);
      const ratio=(Math.max(a,b)+.05)/(Math.min(a,b)+.05);
      return ratio >= 4.5 ? [] : [{text:e.textContent.trim().slice(0,40),ratio}];
    });
  })()`);
  assert.deepEqual(failures, []);
}

try {
  await call('Page.bringToFront');
  await call('Emulation.setDeviceMetricsOverride', { width: 1600, height: 1100, deviceScaleFactor: 1, mobile: false });
  await call('Page.navigate', { url: appUrl });
  await waitFor("document.getElementById('report-metric_cells')?.textContent.trim() === '4,024'");
  await tab('New comparison');
  await select('report-analysis-donor', 'd1');
  await select('report-analysis-celltype', 'CD4 T');
  await select('report-analysis-target', 'PHA');
  await select('report-analysis-reference', 'unstim');
  await evaluate("document.getElementById('report-analysis-run').click()");
  await waitFor("document.getElementById('report-analysis-status').textContent.includes('Confirm the exploratory')");
  await evaluate("document.getElementById('report-analysis-acknowledge').click()");
  await settled();

  for (const readout of ['abundance', 'self_proximity', 'colocalization']) {
    await select('report-analysis-readout', readout);
    await run();
    const provenance = JSON.parse(await download('report-analysis-provenance_json'));
    assert.equal(provenance.request.readout, readout);
    assert.equal(provenance.request.donor, 'd1');
    assert.equal(provenance.request.celltype, 'CD4 T');
    assert.match(provenance.provenance.bundle_md5, /^[a-f0-9]{32}$/);
    assert.match(await download('report-analysis-results_csv'), /exploratory_component_level/);
    await evaluate("document.getElementById('report-analysis-open_result').click()");
    await tab('Volcano');
    await waitFor("document.getElementById('report-report_volcano_plot')?.data?.some(t => t.mode === 'markers' && t.x.length)");
    await tab('New comparison');
  }
  console.log('PASS: three real background analyses, confirmation, same-session navigation, CSV/provenance');
  await checkControlContrast();

  await screenshot('hybrid-analysis');

  await evaluate("document.getElementById('report-analysis-run').click()");
  await waitFor("document.getElementById('report-analysis-run').disabled");
  await evaluate("document.getElementById('report-analysis-cancel').click()");
  await waitFor("document.getElementById('report-analysis-status').textContent.includes('cancelled')");
  await run();
  console.log('PASS: cancellation releases the worker and another run succeeds');

  await tab('Cell 3D');
  await waitFor("!!document.getElementById('report-graph-markers').selectize.options.CD45");
  assert.equal(await evaluate("!!document.querySelector('#report-graph-plot canvas')"), false);
  assert.equal(await evaluate("document.getElementById('report-graph-plot_gif').disabled"), true);
  await select('report-graph-markers', 'CD45');
  await select('report-graph-cell', 'd1_PHA_004e4ae004879a14');
  assert.equal(await evaluate("document.getElementById('report-graph-markers').selectize.settings.maxItems"), 1);
  assert.equal(await evaluate("!!document.querySelector('#report-graph-plot canvas')"), false);
  await evaluate("document.getElementById('report-graph-display').click()");
  await waitFor("HTMLWidgets.getInstance(document.getElementById('report-graph-plot'))?.widget?.pointgroup?.children.length > 0");
  await waitFor("document.querySelector('#report-graph-plot canvas') !== null");
  await waitFor("HTMLWidgets.getInstance(document.getElementById('report-graph-plot')).widget.pointgroup.children.some(t => t.geometry.labels.some(x => x.includes('CD45')))");
  const graph = JSON.parse(await download('report-graph-provenance_json'));
  assert.equal(graph.schema_version, 3);
  assert.equal(graph.highlight_marker, 'CD45');
  await select('report-graph-markers', 'CD4');
  assert.deepEqual(JSON.parse(await download('report-graph-provenance_json')), graph);
  await select('report-graph-pair_1', 'CD4');
  await select('report-graph-pair_2', 'CD45');
  await waitFor("document.getElementById('report-graph-ranking_status').textContent.includes('Ranking complete:')");
  const ranking = await download('report-graph-ranking_csv');
  assert.equal(ranking.trim().split('\n').length, 11);
  await select('report-graph-pair_1', 'CD45');
  await select('report-graph-pair_2', 'CD4');
  assert.equal(await download('report-graph-ranking_csv'), ranking);
  await select('report-graph-rank_view', 'both');
  assert.equal((await download('report-graph-ranking_csv')).trim().split('\n').length, 21);
  assert.deepEqual(JSON.parse(await download('report-graph-provenance_json')), graph);
  assert.match(graph.source_pxl_md5, /^[a-f0-9]{32}$/);
  assert.equal(await evaluate("document.getElementById('report-graph-edges') === null"), true);
  const nodeCsv = (await download('report-graph-nodes_csv')).trim().split('\n');
  assert.equal(nodeCsv[0], '"cell_id","x","y","z","markers"');
  const nodes = nodeCsv.length - 1;
  assert.equal(await evaluate("HTMLWidgets.getInstance(document.getElementById('report-graph-plot')).widget.pointgroup.children.reduce((n,t) => n + t.geometry.attributes.position.count,0)"), nodes);
  assert.equal(await evaluate("HTMLWidgets.getInstance(document.getElementById('report-graph-plot')).widget.pointgroup.children.every(t => t.type === 'Points')"), true);
  assert.equal(await evaluate(`(() => {
    const points=HTMLWidgets.getInstance(document.getElementById('report-graph-plot')).widget.pointgroup.children[0];
    const colors=points.geometry.attributes.color.array, texture=points.material.uniforms.texture1.value.image.data;
    let background=0, highlighted=0;
    for (let i=0;i<colors.length;i+=4) {
      if (Math.abs(colors[i+3]-64/255)<1e-6 && Math.abs(colors[i]-211/255)<1e-6) background++;
      if (colors[i+3]===1 && Math.abs(colors[i]-192/255)<1e-6 && colors[i+1]===0 && colors[i+2]===0) highlighted++;
    }
    return background>0 && highlighted>0 && background+highlighted===colors.length/4 &&
      !texture.some((v,i)=>i%4===3 && v>0 && texture[i-1]===0 && texture[i-2]===0 && texture[i-3]===0);
  })()`), true, 'Faint grey background, opaque red highlights, and no black dot outlines');
  await checkControlContrast();
  await evaluate("document.getElementById('report-graph-plot').closest('.card').querySelector('.bslib-full-screen-enter').click()");
  await waitFor("document.getElementById('report-graph-plot').closest('.card').dataset.fullScreen === 'true'");
  await settled();
  // Wait for the widget and its trackball controls to adopt the fullscreen size.
  await waitFor(`(() => {
    const plot = document.getElementById('report-graph-plot');
    const canvas = plot.querySelector('canvas');
    const graph = HTMLWidgets.getInstance(plot).widget;
    if (!canvas || Math.abs(graph.controls.screen.width - plot.clientWidth) > 1 ||
        Math.abs(graph.controls.screen.height - plot.clientHeight) > 1) return false;
    const r = canvas.getBoundingClientRect();
    return r.width > 0 && r.height > 0 &&
      document.elementFromPoint(r.x + r.width * .45, r.y + r.height * .5) === canvas;
  })()`);
  await evaluate("new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)))");
  const cameraBefore = await evaluate("JSON.stringify(HTMLWidgets.getInstance(document.getElementById('report-graph-plot')).widget.camera.position)");
  const graphBox = await evaluate(`(() => {
    const r=document.querySelector('#report-graph-plot canvas').getBoundingClientRect();
    return {x:r.x+r.width*.45,y:r.y+r.height*.5};
  })()`);
  await call('Input.dispatchMouseEvent', {type:'mouseMoved',...graphBox});
  await call('Input.dispatchMouseEvent', {type:'mousePressed',...graphBox,button:'left',buttons:1,clickCount:1});
  await delay(100);
  for (let step = 1; step <= 8; step++) {
    await call('Input.dispatchMouseEvent', {type:'mouseMoved',x:graphBox.x+90*step/8,y:graphBox.y+35*step/8,button:'left',buttons:1});
    await delay(75);
  }
  await call('Input.dispatchMouseEvent', {type:'mouseReleased',x:graphBox.x+90,y:graphBox.y+35,button:'left',clickCount:1});
  await waitFor(`JSON.stringify(HTMLWidgets.getInstance(document.getElementById('report-graph-plot')).widget.camera.position) !== ${JSON.stringify(cameraBefore)}`);
  await delay(1000);
  const rotatedCamera = await evaluate("JSON.stringify(HTMLWidgets.getInstance(document.getElementById('report-graph-plot')).widget.camera.position)");
  await call('Input.dispatchMouseEvent', {type:'mouseWheel',...graphBox,deltaX:0,deltaY:-120});
  await waitFor(`JSON.stringify(HTMLWidgets.getInstance(document.getElementById('report-graph-plot')).widget.camera.position) !== ${JSON.stringify(rotatedCamera)}`);
  await evaluate("document.getElementById('report-graph-plot_reset').click()");
  await waitFor(`JSON.stringify(HTMLWidgets.getInstance(document.getElementById('report-graph-plot')).widget.camera.position) === ${JSON.stringify(cameraBefore)}`);
  await screenshot('hybrid-cell-3d');
  await evaluate("document.querySelector('.bslib-full-screen-exit').click()");
  await select('report-graph-pair_1', '');
  await select('report-graph-pair_2', '');
  await select('report-graph-cell', 'd1_PHA_004e4ae004879a14');
  await select('report-graph-markers', 'CD45');
  for (let repeat = 0; repeat < 3; repeat++) {
    const savedCamera = await evaluate(`(()=>{
      const g=HTMLWidgets.getInstance(document.getElementById('report-graph-plot')).widget;
      g.camera.position.multiplyScalar(1.1);g.controls.update();g.animate();
      return JSON.stringify(g.camera.position);
    })()`);
    await evaluate("window.previousGraphRenderer=HTMLWidgets.getInstance(document.getElementById('report-graph-plot')).widget.renderer");
    await evaluate("document.getElementById('report-graph-display').click()");
    await waitFor("HTMLWidgets.getInstance(document.getElementById('report-graph-plot')).widget.renderer !== window.previousGraphRenderer");
    await waitFor("window.previousGraphRenderer.getContext().isContextLost()");
    assert.equal(await evaluate("JSON.stringify(HTMLWidgets.getInstance(document.getElementById('report-graph-plot')).widget.camera.position)"), savedCamera);
    assert.equal(await evaluate("document.querySelectorAll('#report-graph-plot canvas').length"), 1);
  }
  await evaluate("document.getElementById('report-graph-plot').scrollIntoView({block:'center'})");
  const hoverPoint = await evaluate(`(()=>{
    const el=document.getElementById('report-graph-plot'),g=HTMLWidgets.getInstance(el).widget;
    const p=g.pointgroup.children[0].geometry.attributes.position;
    const v=new THREE.Vector3().fromBufferAttribute(p,Math.floor(p.count/2)).project(g.camera);
    const r=el.getBoundingClientRect();return {x:r.left+(v.x+1)*r.width/2,y:r.top+(1-v.y)*r.height/2};
  })()`);
  await call('Input.dispatchMouseEvent', {type:'mouseMoved',...hoverPoint});
  await waitFor("document.querySelector('#report-graph-plot .infobox').textContent.includes('Markers:')");
  const pngHasPoints = await evaluate(`(async()=>{
    let url;const click=HTMLAnchorElement.prototype.click;
    try {HTMLAnchorElement.prototype.click=function(){url=this.href;};document.getElementById('report-graph-plot_png').click();}
    finally {HTMLAnchorElement.prototype.click=click;}
    const bitmap=await createImageBitmap(await (await fetch(url)).blob());
    const canvas=document.createElement('canvas');canvas.width=bitmap.width;canvas.height=bitmap.height;
    const context=canvas.getContext('2d');context.drawImage(bitmap,0,0);
    const pixels=context.getImageData(0,0,canvas.width,canvas.height).data;
    return pixels.some((v,i)=>i%4===3&&v>0&&pixels[i-1]<240);
  })()`);
  assert.equal(pngHasPoints, true);
  console.log('PASS: repeated Display releases old WebGL contexts; marker hover and nonblank PNG work');
  console.log('PASS: prepared 3D graph, all nodes, coloring, exports, full-screen, mouse rotation and zoom');

  const gifState = await evaluate(`(() => {
    const el=document.getElementById('report-graph-plot'), g=HTMLWidgets.getInstance(el).widget;
    return JSON.stringify({position:g.camera.position,up:g.camera.up,target:g.controls.target,
      width:g.renderer.domElement.width,height:g.renderer.domElement.height});
  })()`);
  assert.deepEqual(await evaluate("['seconds','fps','pixels'].map(x=>Number(document.getElementById('report-graph-plot_gif_'+x).value))"), [6,15,600]);
  await evaluate("document.getElementById('report-graph-plot_gif_seconds').value=0;document.getElementById('report-graph-plot_gif').click()");
  await waitFor("document.getElementById('report-graph-plot_gif_status').textContent.includes('Use whole numbers')");
  await evaluate("document.getElementById('report-graph-plot_gif_seconds').value=20;document.getElementById('report-graph-plot_gif_fps').value=30;document.getElementById('report-graph-plot_gif').click()");
  await waitFor("document.getElementById('report-graph-plot_gif_status').textContent.includes('Limit: 300 frames')");
  await evaluate(`(() => {
    document.getElementById('report-graph-plot_gif_seconds').value=2;
    document.getElementById('report-graph-plot_gif_fps').value=7;
    document.getElementById('report-graph-plot_gif_pixels').value=256;
    const click=HTMLAnchorElement.prototype.click;
    window.restoreGifClick=()=>{HTMLAnchorElement.prototype.click=click;};
    const Renderer=THREE.WebGLRenderer;
    window.restoreGifRenderer=()=>{THREE.WebGLRenderer=Renderer;};
    THREE.WebGLRenderer=function(...args){const renderer=new Renderer(...args);window.gifRenderer=renderer;return renderer;};
    HTMLAnchorElement.prototype.click=function(){
      if(this.download.endsWith('.gif')){
        window.gifFilename=this.download;
        fetch(this.href).then(r=>r.arrayBuffer()).then(b=>{window.downloadedGif=b;});
      } else click.call(this);
    };
    document.getElementById('report-graph-plot_gif').click();
  })()`);
  assert.equal(await evaluate("document.getElementById('report-graph-display').disabled"), true);
  await waitFor("!!window.downloadedGif");
  assert.equal(await evaluate("window.gifRenderer.getContext().isContextLost()"), true);
  await evaluate("window.restoreGifClick()");
  const gifCheck = await evaluate(`(async () => {
    const bytes=new Uint8Array(window.downloadedGif);
    const decoder=new ImageDecoder({data:bytes,type:'image/gif'});
    await decoder.tracks.ready;
    let duration=0, first, changed=false, red=false, dimensions=true;
    for(let i=0;i<decoder.tracks.selectedTrack.frameCount;i++){
      const {image}=await decoder.decode({frameIndex:i});
      dimensions &&= image.displayWidth===256 && image.displayHeight===256;
      duration+=image.duration;
      const canvas=document.createElement('canvas');canvas.width=canvas.height=256;
      const ctx=canvas.getContext('2d');ctx.drawImage(image,0,0);image.close();
      const pixels=ctx.getImageData(0,0,256,256).data;
      if(!first)first=pixels;else changed ||= pixels.some((v,j)=>v!==first[j]);
      red ||= pixels.some((v,j)=>j%4===0 && v>80 && v>pixels[j+1]*1.5 && v>pixels[j+2]*1.5);
    }
    const result={header:String.fromCharCode(...bytes.slice(0,6)),frames:decoder.tracks.selectedTrack.frameCount,
      duration,dimensions,changed,red,filename:window.gifFilename};
    decoder.close();return result;
  })()`);
  assert.equal(gifCheck.header, 'GIF89a');
  assert.equal(gifCheck.frames, 14);
  assert.equal(gifCheck.duration, 2000000);
  assert.equal(gifCheck.dimensions && gifCheck.changed && gifCheck.red, true);
  assert.match(gifCheck.filename, /d1_PHA_004e4ae004879a14-rotating.gif$/);
  await writeFile('output/hybrid-cell-rotating.gif', Buffer.from(await evaluate("Array.from(new Uint8Array(window.downloadedGif))")));
  await evaluate("document.getElementById('report-graph-plot_gif').click();document.getElementById('report-graph-plot_gif_cancel').click()");
  await waitFor("document.getElementById('report-graph-plot_gif_status').textContent === 'Export cancelled.'");
  assert.equal(await evaluate("document.getElementById('report-graph-display').disabled || document.getElementById('report-graph-plot_gif').disabled"), false);
  assert.equal(await evaluate("window.gifRenderer.getContext().isContextLost()"), true);
  await evaluate(`(() => {
    const OriginalWorker=window.Worker;
    try {window.Worker=class{constructor(){throw new Error('Simulated worker failure');}};
      document.getElementById('report-graph-plot_gif').click();
    } finally {window.Worker=OriginalWorker;}
  })()`);
  await waitFor("document.getElementById('report-graph-plot_gif_status').textContent.includes('Simulated worker failure')");
  assert.equal(await evaluate("document.getElementById('report-graph-display').disabled || document.getElementById('report-graph-plot_gif').disabled"), false);
  assert.equal(await evaluate(`(() => {
    const el=document.getElementById('report-graph-plot'),g=HTMLWidgets.getInstance(el).widget;
    return JSON.stringify({position:g.camera.position,up:g.camera.up,target:g.controls.target,
      width:g.renderer.domElement.width,height:g.renderer.domElement.height});
  })()`), gifState);
  await evaluate("window.restoreGifRenderer()");
  console.log('PASS: GIF settings, decoded frames/duration/dimensions/red highlights, rotation, cancellation, failure recovery, unchanged view');

  await call('Emulation.setDeviceMetricsOverride', { width: 390, height: 844, deviceScaleFactor: 1, mobile: true });
  await tab('New comparison');
  assert.ok(await evaluate('document.documentElement.scrollWidth <= window.innerWidth + 2'));
  await screenshot('hybrid-mobile');
  console.log('PASS: mobile navigation/layout without horizontal page overflow');
} finally {
  socket.close();
}
