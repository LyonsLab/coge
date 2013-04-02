function set_input(id, text) {
	var obj = $('#'+id);
	var val = obj.val();
	if (!text) {
		text = 'Search';
	}	
	
	if ( val && val != text ) {
		obj.css({"font-style": "normal", "color": "black"});
	}
	else {
		if (obj.is(":focus")) {
			obj.val('');
			obj.css({"font-style": "normal", "color": "black"});
		}
		else {
			obj.val(text);
			obj.css({"font-style": "italic", "color": "gray"});
		}
	}

	// Init event handlers
	$('#'+id).unbind('blur').bind('blur', function() { set_input(id, text); })
			 .unbind('focus').bind('focus', function() { set_input(id, text); })
			 .unbind('click').bind('click', function() { set_input(id, text); });
}