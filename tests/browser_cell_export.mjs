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

try {
  await call('Page.navigate', {url: appUrl});
  await waitFor("!!window.Shiny?.shinyapp?.$socket");
  await waitFor("[...document.querySelectorAll('a[data-bs-toggle=tab]')].some(x=>x.innerText==='Cell 3D')");
  await tab('Cell 3D');
  await waitFor("!!document.getElementById('report-graph-markers').selectize.options.CD82");
  await select('report-graph-markers', 'CD82');
  await waitFor("!!Shiny.shinyapp.$inputValues['report-graph-cell']");
  await evaluate("document.getElementById('report-graph-display').click()");
  await waitFor("HTMLWidgets.getInstance(document.getElementById('report-graph-plot'))?.widget?.pointgroup?.children.length>0");
  await settled();
  const readModel = () => evaluate(`fetch(document.getElementById('report-graph-model_glb').href).then(r=>r.arrayBuffer()).then(b=>{
    const v=new DataView(b); if(v.getUint32(0,true)!==0x46546c67) throw Error('Not GLB');
    return JSON.parse(new TextDecoder().decode(new Uint8Array(b,20,v.getUint32(12,true))));
  })`);
  const before = await readModel();
  assert.equal(before.extras.highlight_marker,'CD82');
  assert.equal(before.materials[0].alphaMode,'OPAQUE');
  assert.equal(before.extras.appearance.background_radius,0.005);
  await evaluate(`(()=>{
    for(const [key,value] of Object.entries({background_color:'#123456',highlight_color:'#00FF00',
      background_opacity:0.6,highlight_opacity:1,background_radius:0.005,highlight_radius:0.007,
      glb_background_opacity:0.8,glb_highlight_opacity:1,glb_background_radius:0.008,glb_highlight_radius:0.009})) {
      const e=document.getElementById('report-graph-'+key);
      e.value=value; e.dispatchEvent(new Event('change',{bubbles:true}));
    }
  })()`);
  await settled();
  assert.deepEqual((await readModel()).extras.appearance,before.extras.appearance);
  await evaluate("document.getElementById('report-graph-apply_appearance').click()");
  await waitFor(`(()=>{
    const g=HTMLWidgets.getInstance(document.getElementById('report-graph-plot')).widget;
    return g.pointgroup.children.some(p=>p.geometry.attributes.color.array.some((v,i,a)=>
      i%4===0 && v===0 && a[i+1]===1 && a[i+2]===0));
  })()`);
  const after = await readModel();
  assert.equal(after.extras.appearance.highlight_radius,0.009);
  assert.equal(after.extras.appearance.background_opacity,0.8);
  assert.deepEqual(after.materials[1].pbrMetallicRoughness.baseColorFactor,[0,1,0,1]);
  await evaluate(`(()=>{
    for(const [key,value] of Object.entries({seconds:1,fps:1,pixels:128}))
      document.getElementById('report-graph-plot_gif_'+key).value=value;
    document.getElementById('report-graph-plot_gif').click();
  })()`);
  await waitFor("document.getElementById('report-graph-plot_gif_status').textContent.startsWith('GIF downloaded:')");
  assert.equal(await evaluate("document.getElementById('report-graph-apply_appearance').disabled"),false);
  console.log('PASS: actual GLB download, pending/applied appearance, updated WebGL colors, and GIF export');
} finally {socket.close();}
