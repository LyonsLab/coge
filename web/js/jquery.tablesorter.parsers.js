$.tablesorter.addParser({
	id: 'genericLink',
	is: function(s) {
		var exp = /<a[^>]*">(.*)<\/a>/.test(s);
		if(exp) this.type = (/<a[^>]*">(\d+)<\/a>/.test(s)) ? 'numeric' : 'text';
		return exp;
	},
	format: function(s) {
		var val = s.toLowerCase().match(/<a[^>]*">(.*)<\/a>/)[1];
		// format your data for normalization
		return (this.type == 'text') ?  $.trim(val) : $.tablesorter.formatFloat(val);   
	},
	type: 'text'
});

$.tablesorter.addParser({
	id: 'currency',
	is: function(s) {
		return /^[£$€?.]/.test(s);
	},
	format: function(s) {
		return $.tablesorter.formatFloat(s.replace(new RegExp(/[^0-9.]/g),''));
	},
	type: 'numeric'
});

$.tablesorter.addParser({
	id: 'integer',
	is: function(s) {
		return /^\d+$/.test(s);
	},
	format: function(s) {
		return $.tablesorter.formatFloat(s);
	},
	type: 'numeric'
});

$.tablesorter.addParser({
	id: 'floating',
	is: function(s) {
		return /^(\+|-)?[0-9]+\.[0-9]+((E|e)(\+|-)?[0-9]+)?$/.test(s);
	},
	format: function(s) {
		return $.tablesorter.formatFloat(s.replace(new RegExp(/,|\./),''));
	},
	type: 'numeric'
});

$.tablesorter.addParser({
	id: 'ipAddress',
	is: function(s) {
		return /^\d{2,3}[\.]\d{2,3}[\.]\d{2,3}[\.]\d{2,3}$/.test(s);
	},
	format: function(s) {
		var a = s.split('.');
		var r = '';
		for (var i = 0, item; item = a[i]; i++) {
		   if(item.length == 2) {
				r += '0' + item;
		   } else {
				r += item;
		   }
		}
		return $.tablesorter.formatFloat(s);
	},
	type: 'numeric'
});

$.tablesorter.addParser({
	id: 'url',
	is: function(s) {
		return /(https?|ftp|file):\/\//.test(s);
	},
	format: function(s) {
		return jQuery.trim(s.replace(new RegExp(/(https?|ftp|file):\/\//),''));
	},
	type: 'text'
});

$.tablesorter.addParser({
	id: 'isoDate',
	is: function(s) {
		return /^\d{4}[\/-]\d{1,2}[\/-]\d{1,2}$/.test(s);
	},
	format: function(s) {
		return parseFloat((s != "") ? new Date(s.replace(new RegExp(/-/g),'/')).getTime() : "0");
	},
	type: 'numeric'
});

$.tablesorter.addParser({
	id: 'usLongDate',
	is: function(s) {
		return /^[A-Za-z]{3,10}\.? [0-9]{1,2}, ([0-9]{4}|'?[0-9]{2}) (([0-2]?[0-9]:[0-5][0-9])|([0-1]?[0-9]:[0-5][0-9]\s(AM|PM)))$/.test(s);
	},
	format: function(s) {
		return $.tablesorter.formatFloat(new Date(s).getTime());
	},
	type: 'numeric'
});
/*
$.tablesorter.addParser({
	id: 'shortDate',
	is: function(s) {
		return /\d{1,2}[\/-]\d{1,2}[\/-]\d{2,4}/.test(s);
	},
	format: function(s,defaults) {
		s = s.replace(new RegExp(/-/g),'/');
		if(defaults.dateFormat == "mm/dd/yyyy" || defaults.dateFormat == "mm-dd-yyyy") {
			// reformat the string in ISO format
			s = s.replace(new RegExp(/(\d{1,2})[\/-](\d{1,2})[\/-](\d{4})/), '$3/$1/$2');
		} else if(defaults.dateFormat == "dd/mm/yyyy" || defaults.dateFormat == "dd-mm-yyyy") {
			// reformat the string in ISO format
			s = s.replace(new RegExp(/(\d{1,2})[\/-](\d{1,2})[\/-](\d{4})/), '$3/$2/$1');
		} else if(defaults.dateFormat == "dd/mm/yy" || defaults.dateFormat == "dd-mm-yy") {
			s = s.replace(new RegExp(/(\d{1,2})[\/-](\d{1,2})[\/-](\d{2})/), '$1/$2/$3');	
		}
		return $.tableSorter.utils.formatFloat((new Date(s)).getTime());
	},
	type: 'numeric'
});
*/
$.tablesorter.addParser({
    id: 'time',
    is: function(s) {
        return /^(([0-2]?[0-9]:[0-5][0-9])|([0-1]?[0-9]:[0-5][0-9]\s(am|pm)))$/.test(s);
    },
    format: function(s) {
        return $.tableSorter.utils.formatFloat((new Date("2000/01/01 " + s)).getTime());
    },
  type: 'numeric'
});