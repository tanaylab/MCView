// Gene choices of virtualSelect inputs (see update_gene_select in R/utils_selectors.R): the
// gene list is sent once per session and reused by every gene selector.
(function() {
    var choices = {};
    Shiny.addCustomMessageHandler("mcview-gene-choices", function(msg) {
        choices[msg.key] = msg.options;
    });
    Shiny.addCustomMessageHandler("mcview-update-gene-select", function(msg) {
        var el = document.getElementById(msg.id);
        var binding = el && $(el).data("shiny-input-binding");
        if (!binding) {
            return;
        }
        var data = {};
        if (choices[msg.key]) data.options = choices[msg.key];
        if (msg.hasOwnProperty("value")) data.value = msg.value;
        if (msg.hasOwnProperty("label")) data.label = msg.label;
        binding.receiveMessage(el, data);
    });
})();
