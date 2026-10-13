$(function () {
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
