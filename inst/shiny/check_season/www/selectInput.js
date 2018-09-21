$(document).on("keydown", function(e) {
    var obj = $(".selectized").selectize().last()[0].selectize;

    var val = obj.getValue();
    var ind = obj.options[val]["$order"];

    var code = e.keyCode;
    var pos;
    if (code == 37) { // || code == 38
        pos = -1;
    } else if (code == 39) { // || code == 40
        pos = 1;
    }

    if (pos) {
        var newval, item;
        var newind = ind + pos;

        // look for adjacent item
        if (newind >= 0 && newind < obj.order) {
            for (var i in obj.options) {
                item = obj.options[i];
                if (item.$order === ind + pos) {
                    newval = item.value;
                    break;
                }
            }
        }

        if (newval) obj.setValue(newval);
    }
    Shiny.onInputChange("keycode", [e.which, e.timeStamp]);
});
