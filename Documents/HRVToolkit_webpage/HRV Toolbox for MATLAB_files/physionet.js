$(document).ready(function(){
// debug ("--start--");


// initialise plugin
var $pn = $('.sf-menu').superfish({
	//add options here if required
//popUpSelector: "table",
//disableHI: true,
autoArrows: false,

});

/*$("#example-wrapper").find("table,img").attr ("role", "presentation")
.find("img").removeAttr("alt");
*/

makeAccessible ($pn, "a11y-", {
open: function ($node) {
$node
//.children("a:first")
.superfish("show");
}, // open

afterOpen: function ($node) {
}, // afterOpen

close: function ($node) {
$node
//.children("a:first")
.superfish("hide");
} // close

}); // makeAccessible



}); // ready

//alert ("complex.js loaded");
