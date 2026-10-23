$(function () {
  Shiny.addCustomMessageHandler("precomputed-spatial-size", function (message) {
    const region = document.getElementById(message.id);
    if (!region) return;
    // Leave at least 16 pixels per protein, plus space for labels and legends.
    region.style.setProperty("--spatial-panel-min", `${Math.max(480, message.proteins * 16 + 260)}px`);
    window.dispatchEvent(new Event("resize"));
  });
  Shiny.addCustomMessageHandler("precomputed-spatial-availability", function (message) {
    const tab = document.getElementById(message.id)?.querySelector('a[data-value="Spatial patterns"]');
    if (!tab) return;
    tab.classList.toggle("disabled", message.disabled);
    tab.setAttribute("aria-disabled", String(message.disabled));
    if (message.disabled) {
      tab.setAttribute("tabindex", "-1");
      tab.title = "Select Self-proximity or Colocalization to view the spatial matrix.";
    } else {
      tab.removeAttribute("tabindex");
      tab.removeAttribute("title");
    }
  });
});
