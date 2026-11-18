// Gene choices of virtualSelect inputs (see update_gene_select in R/utils_selectors.R): the
// gene list is sent once per session and reused by every gene selector.
(function() {
    var choices = {};
    Shiny.addCustomMessageHandler("mcview-gene-choices", function(msg) {
        // the labels that differ from the gene name come by index
        var label = msg.genes.slice();
        for (var i in msg.labels) label[+i] = msg.labels[i];
        // the format shinyWidgets sends for named choices
        choices[msg.key] = {type: "transpose", choices: {label: label, value: msg.genes}};
    });
    Shiny.addCustomMessageHandler("mcview-update-gene-select", function(msg) {
        var el = document.getElementById(msg.id);
        var binding = el && $(el).data("shiny-input-binding");
        if (!binding) {
            return;
        }
        var data = {};
        var options = choices[msg.key];
        if (options && msg.plain_names) {
            options = {type: "transpose", choices: {label: options.choices.value, value: options.choices.value}};
        }
        if (options) data.options = options;
        if (msg.hasOwnProperty("value")) data.value = msg.value;
        if (msg.hasOwnProperty("label")) data.label = msg.label;
        binding.receiveMessage(el, data);
    });
})();
