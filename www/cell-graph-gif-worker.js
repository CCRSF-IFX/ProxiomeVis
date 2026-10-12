import { GIFEncoder, quantize, applyPalette } from './vendor/gifenc-1.0.3.js';

const gif = GIFEncoder();
self.onmessage = ({data}) => {
  try {
    if (data.finish) {
      gif.finish();
      const bytes = gif.bytes();
      self.postMessage({bytes}, [bytes.buffer]);
    } else {
      const rgba = new Uint8Array(data.rgba);
      const palette = quantize(rgba, 256);
      gif.writeFrame(applyPalette(rgba, palette), data.size, data.size,
        {palette, delay: data.delay, repeat: 0, dispose: 1});
      if (gif.bytesView().byteLength > 64 * 1024 * 1024) {
        throw new Error('GIF exceeds 64 MB. Reduce duration, FPS, or image size.');
      }
      self.postMessage({frame: true});
    }
  } catch (error) {
    self.postMessage({error: error.message});
  }
};
