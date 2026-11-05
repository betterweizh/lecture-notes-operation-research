/* Open folded blocks (Proof / Solution dropdowns, <details>) that contain
   the target of a link, so that following a link to an equation, table or
   figure inside a folded solution shows it instead of a closed bar. Runs on
   page load and whenever the URL fragment changes. */
(function () {
  function openTarget() {
    var id = decodeURIComponent(window.location.hash.slice(1));
    if (!id) return;
    var target = document.getElementById(id);
    if (!target) return;
    var opened = false;
    for (var el = target; el; el = el.parentElement) {
      if (el.tagName === "DETAILS" && !el.open) {
        el.open = true;
        opened = true;
      }
    }
    if (opened) target.scrollIntoView();
  }
  document.addEventListener("DOMContentLoaded", openTarget);
  window.addEventListener("hashchange", openTarget);
})();
