// Serve output/renderer-benchmark on :3840. PROXIOME_CHROME is the Chromium executable.
import assert from 'node:assert/strict';
import {writeFile, mkdtemp, rm} from 'node:fs/promises';
import {spawn} from 'node:child_process';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
const base = process.env.PROXIOME_BENCH_URL || 'http://127.0.0.1:3840';
const {cases, versions} = await (await fetch(`${base}/cases.json`)).json();
let socket;
let serial = 0;
const pending = new Map();
const onmessage = ({data}) => {
  const value = JSON.parse(data);
  pending.get(value.id)?.(value);
  pending.delete(value.id);
};
const call = (method, params = {}) => new Promise((resolve, reject) => {
  const id = ++serial;
  const timer = setTimeout(() => reject(new Error(`Timeout: ${method}`)), 60000);
  pending.set(id, value => {
    clearTimeout(timer);
    if (value.error) reject(new Error(JSON.stringify(value.error)));
    else resolve(value.result);
  });
  socket.send(JSON.stringify({id, method, params}));
});
const evaluate = async expression => {
  const result = await call('Runtime.evaluate', {expression, returnByValue: true, awaitPromise: true});
  if (result.exceptionDetails) throw new Error(JSON.stringify(result.exceptionDetails));
  return result.result.value;
};
const delay = ms => new Promise(r => setTimeout(r, ms));
async function ready() {
  for (let i = 0; i < 200; i++) {
    if (await evaluate('!!window.benchReady')) return;
    await delay(100);
  }
  throw new Error('Widget did not render');
}
async function metrics() {
  return Object.fromEntries((await call('Performance.getMetrics')).metrics.map(x => [x.name, x.value]));
}
const results = [];
  for (const item of cases) {
    // Fresh browser per case avoids retaining old WebGL contexts and back/forward pages.
    const profile = await mkdtemp(join(tmpdir(), 'proxiome-renderer-'));
    await writeFile(join(profile, 'SwiftShader.ini'), '[Processor]\nThreadCount=4\n');
    const chrome = process.env.PROXIOME_CHROME || 'chromium';
    const flags = [
      '--headless', '--no-sandbox', '--single-process', '--no-zygote',
      '--disable-extensions', '--disable-background-networking', '--disable-dev-shm-usage',
      '--disable-features=BackForwardCache', '--num-raster-threads=1',
      '--js-flags=--single-threaded',
      '--use-angle=swiftshader', '--enable-unsafe-swiftshader',
      '--remote-debugging-port=0', `--user-data-dir=${profile}`, 'about:blank'
    ];
    const cpus = process.env.PROXIOME_BENCH_CPUS;
    const browser = spawn(cpus ? 'taskset' : chrome,
      cpus ? ['-c', cpus, chrome, ...flags] : flags,
      {cwd: profile, stdio: ['ignore', 'ignore', 'pipe']});
    const exited = new Promise(resolve => browser.once('close', resolve));
    try {
      const endpoint = await new Promise((resolve, reject) => {
        let output = '';
        const timer = setTimeout(() => reject(new Error('Chromium startup timeout')), 30000);
        browser.once('error', error => {clearTimeout(timer); reject(error);});
        browser.once('exit', code => {
          clearTimeout(timer);
          if (code) console.error(output.slice(-2000));
          reject(new Error(`Chromium exited (${code}): ${output.slice(-2000)}`));
        });
        browser.stderr.on('data', data => {
          output += data.toString();
          const match = output.match(/DevTools listening on (ws:\/\/[^\s]+)/);
          if (match) {clearTimeout(timer); resolve(new URL(match[1]));}
        });
      });
      let page;
      for (let attempt = 0; attempt < 100 && !page; attempt++) {
        const tabs = await (await fetch(`http://${endpoint.host}/json`)).json();
        page = tabs.find(t => t.type === 'page');
        if (!page) await delay(100);
      }
      assert.ok(page, 'Chromium must create its initial page');
      socket = new WebSocket(page.webSocketDebuggerUrl);
      await new Promise(r => socket.addEventListener('open', r, {once: true}));
      socket.onmessage = onmessage;
      socket.onclose = () => {
        for (const finish of pending.values()) finish({error: {message: 'Browser disconnected'}});
        pending.clear();
      };
      await call('Performance.enable');
      await call('Emulation.setDeviceMetricsOverride', {width: 1100, height: 800, deviceScaleFactor: 1, mobile: false});
      // Warm the renderer's dependencies using the smallest cell.
      const renderer = item.renderer === 'plotly' ? 'plotly' : 'threejs';
      await evaluate('window.benchReady=null');
      await call('Page.navigate', {url: `${base}/${cases.find(x => x.renderer === renderer).file}`});
      await ready();
      await evaluate('window.benchReady=null');
      await call('Page.navigate', {url: `${base}/${item.file}`});
      await ready();
      const initial = await evaluate(`(()=>{
        const el=document.getElementById('bench');
        const three=HTMLWidgets.getInstance(el)?.widget;
        const nodes=three ? three.pointgroup.children.reduce((n,t)=>n+(t.geometry.attributes?.position?.count||0),0) : el.data.reduce((n,t)=>n+t.x.length,0);
        const gl=el.querySelector('canvas').getContext('webgl') || el.querySelector('canvas').getContext('webgl2');
        const ext=gl?.getExtension('WEBGL_debug_renderer_info');
        return {ready_ms:window.benchReady,nodes,gpu:ext?gl.getParameter(ext.UNMASKED_RENDERER_WEBGL):null};
      })()`);
      assert.equal(initial.nodes, item.nodes);
      await delay(700);
      await call('HeapProfiler.collectGarbage');
      const before = await metrics();
      await delay(1000);
      const idle = await metrics();
      const cameraBefore = await evaluate("(()=>{const e=document.getElementById('bench');const t=HTMLWidgets.getInstance(e)?.widget;return JSON.stringify(t?t.camera.position:e._fullLayout.scene._scene.getCamera());})()");
      await evaluate(`(()=>{window.framesSeen=[];window.captureFrames=true;function frame(t){if(!window.captureFrames)return;window.framesSeen.push(t);requestAnimationFrame(frame);}requestAnimationFrame(frame);})()`);
      const start = performance.now();
      await call('Input.dispatchMouseEvent', {type: 'mousePressed', x: 500, y: 350, button: 'left', buttons: 1, clickCount: 1});
      let step = 0;
      while (performance.now() - start < 2000) {
        await call('Input.dispatchMouseEvent', {type: 'mouseMoved', x: 500 + 120 * Math.sin(step / 5), y: 350 + 60 * Math.cos(step / 5), button: 'left', buttons: 1});
        step++;
        await delay(30);
      }
      await call('Input.dispatchMouseEvent', {type: 'mouseReleased', x: 500, y: 350, button: 'left', buttons: 0, clickCount: 1});
      const frameStats = await evaluate(`(()=>{window.captureFrames=false;const a=window.framesSeen;const dt=a.slice(1).map((t,i)=>t-a[i]).sort((a,b)=>a-b);return {raf_fps:(a.length-1)*1000/(a.at(-1)-a[0]),p95_frame_ms:dt[Math.floor(dt.length*.95)]};})()`);
      const cameraAfter = await evaluate("(()=>{const e=document.getElementById('bench');const t=HTMLWidgets.getInstance(e)?.widget;return JSON.stringify(t?t.camera.position:e._fullLayout.scene._scene.getCamera());})()");
      assert.notEqual(cameraBefore, cameraAfter, `${item.renderer}: mouse rotation must change the camera`);
      await delay(1000);
      const restBefore = await metrics();
      await delay(1000);
      const restAfter = await metrics();
      await call('HeapProfiler.collectGarbage');
      const after = await metrics();
      const row = {...item, ...initial, ...frameStats,
        idle_task_ms: (idle.TaskDuration - before.TaskDuration) * 1000,
        post_drag_idle_task_ms: (restAfter.TaskDuration - restBefore.TaskDuration) * 1000,
        js_heap_mb: after.JSHeapUsedSize / 1024 ** 2};
      results.push(row);
      console.log(JSON.stringify(row));
      await writeFile('output/renderer-benchmark/browser-results.json', JSON.stringify({versions, results}, null, 2));
    } finally {
      socket?.close();
      browser.kill('SIGTERM');
      await exited;
      await rm(profile, {recursive: true, force: true});
    }
  }
