/*
 * Language switch between the editions of the book.
 *
 * - Remembers the edition of the current page in localStorage, so the root
 *   page (/index.html) opens the reader's edition next time.
 * - On a click on a switch link (.lang-switch__link), stores the target
 *   edition and keeps the URL hash when it is a label that both editions
 *   share (equations, theorems, figures, sections, ...).
 *
 * The storage key and the current edition come from the attributes of this
 * <script> tag (set by apply_language() in docs/_shared/conf_common.py).
 * Every storage access is wrapped in try/catch (private mode, blocked
 * storage).
 */
(function () {
  "use strict";

  var script = document.currentScript;
  var data = (script && script.dataset) || {};
  var KEY = data.storageKey || "orln-preferred-language";
  var current =
    data.currentLang ||
    (document.documentElement.lang || "en").toLowerCase().split("-")[0];

  // Anchors that are identical in every edition.
  var SHARED_HASH =
    /^#(equation-eq|eq|thm|lem|def|prop|cor|asm|alg|ex|rem|fig|tab|sec)-/;

  function store(lang) {
    try {
      window.localStorage.setItem(KEY, lang);
    } catch (e) {
      /* storage unavailable: ignore */
    }
  }

  // Remember the edition of the page being read.
  store(current);

  document.addEventListener("click", function (event) {
    var target = event.target;
    var link = target && target.closest ? target.closest(".lang-switch__link") : null;
    if (!link) {
      return;
    }
    var lang = link.getAttribute("data-lang");
    if (lang) {
      store(lang);
    }
    var hash = window.location.hash;
    if (hash && SHARED_HASH.test(hash)) {
      event.preventDefault();
      window.location.href = link.href.split("#")[0] + hash;
    }
  });
})();
