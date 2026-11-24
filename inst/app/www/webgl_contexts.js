// Browsers keep only ~16 WebGL contexts per page and drop the oldest when more are created. A WebGL
// plotly plot (toWebGL(): the projections and the gene scatters) holds two, and the plots of every
// visited tab stayed alive, so after 4-5 tabs the plots of the first ones were blank when the user
// came back to them. Free the WebGL plots that are hidden (in another tab) and draw them again, from
// the value shiny last sent, when they are shown.
(function() {
    function widgets() {
        return document.querySelectorAll(".plotly.html-widget");
    }

    // the last value of each plotly output: what its binding needs to draw it again
    $(document).on("shiny:value", function(e) {
        var el = e.target;
        if (el && el.classList && el.classList.contains("plotly")) {
            el.mcviewValue = e.value;
            el.mcviewReleased = false;
        }
    });

    function release(el) {
        if (el.mcviewReleased || !el.mcviewValue || !window.Plotly) return;
        if (!el.querySelector("canvas.gl-canvas")) return;
        try {
            // the data and the layout (with the user's zoom) stay in el.mcviewValue.x
            window.Plotly.purge(el);
            el.mcviewReleased = true;
        } catch (err) {
            console.warn("mcview: could not release the plot " + el.id, err);
        }
    }

    function restore(el) {
        if (!el.mcviewReleased) return;
        el.mcviewReleased = false;
        var adapter = $(el).data("shiny-output-binding");
        if (!adapter || !adapter.binding || !el.mcviewValue) return;
        try {
            // as when shiny sends the value: draws the plot and registers its events
            adapter.binding.renderValue(el, el.mcviewValue);
        } catch (err) {
            console.warn("mcview: could not draw the plot " + el.id, err);
        }
    }

    function sync() {
        var els = widgets();
        // free first, so that the contexts are available to the plots that are drawn
        for (var i = 0; i < els.length; i++) {
            if (els[i].offsetParent === null) release(els[i]);
        }
        for (var j = 0; j < els.length; j++) {
            if (els[j].offsetParent !== null) restore(els[j]);
        }
    }

    $(document).on("shown.bs.tab", function() {
        setTimeout(sync, 0);
    });

    // A plot that is redrawn several times (e.g. a new comparison) makes new contexts before the browser
    // collects those of its previous draws, so the browser can still drop the context of a plot that is
    // shown. Draw such a plot again. (The canvases of freed plots are not in the page: only the context
    // of a live plot gets here. The event doesn't bubble, hence the capturing listener.)
    document.addEventListener("webglcontextlost", function(e) {
        var el = e.target && e.target.closest && e.target.closest(".plotly.html-widget");
        if (!el || el.mcviewReleased || el.mcviewRedraw || !el.mcviewValue || el.offsetParent === null) return;
        el.mcviewRedraw = setTimeout(function() {
            el.mcviewRedraw = null;
            if (el.offsetParent === null || !window.Plotly) return;
            console.info("mcview: drawing the plot " + el.id + " again, its WebGL context was lost");
            try {
                window.Plotly.purge(el);
                el.mcviewReleased = true;
                restore(el);
            } catch (err) {
                console.warn("mcview: could not draw the plot " + el.id, err);
            }
        }, 200);
    }, true);
})();
