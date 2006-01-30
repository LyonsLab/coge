	// Flooble.com's Animated Text script. Will animate a specified 
	// bit of text (determined by the ID of containing tag) by 
	// highlighting it with specified color one character at a time 
	// in a moving pattern.
	//
	// Summary of use: 
	//     call animate(tagID, color); where "tagID" is the ID 
	//     of the tag that contains text to be animated,
	//     and "color" is the color to use to highlight the text with.
	//
	// For more information, and detailed instructions, see 
	//     http://www.flooble.com/scripts/animate.php
	//
	// Copyright (c) 2002 by Animus Pactum Consulting Inc.
	// This script comes with no warranties whatsoever. 
	// Animus Pactum Consulting will not be responsible
	// for any damages resulting from its use.

        var ie4 = false;
        if(document.all) {
                ie4 = true; 
        }       
        function setContent(name, value) {
                var d;  
                if (ie4) { 
                        d = document.all[name];
                } else {
                        d = document.getElementById(name);
                }       
                d.innerHTML = value;    
        }       

	function getContent(name) {
		var d;
                if (ie4) {
                        d = document.all[name];
                } else {
                        d = document.getElementById(name);
                }
                return d.innerHTML;
	}

        function setColor(name, value) {
                var d;  
                if (ie4) { 
                        d = document.all[name];
                } else {
                        d = document.getElementById(name);
                }
                d.style.color = value;  
        }

	function getColor(name) {
                var d;
                if (ie4) {
                        d = document.all[name];
                } else {
                        d = document.getElementById(name);
                }
                return d.style.color;
        }

        function animate(name, col) {
		var value = getContent(name);
		if (value.indexOf('<span') >= 0) { return; }
		var length = 0;
                var str = '';
		var ch;
		var token = '';
		var htmltag = false;	
                for (i = 0; i < value.length; i++) {
			ch = value.substring(i, i+1);
			if (i < value.length - 1) { nextch = value.substring(i+1, i+2); } else { nextch = ' '; }
			token += ch;
			if (ch == '<' && '/aAbBpPhHiIoOuUlLtT'.indexOf(nextch) >= 0) { htmltag = true; }
			if (ch == '>' && htmltag) { htmltag = false; }
			if (!htmltag && ch.charCodeAt(0) > 30 && ch != ' ' && ch != '\n') {		
                        	str += '<span id="' + name + '_' + length + '">' + token + '</span>';
				token = '';
				length++;
			}
                }
                setContent(name, str);
                command = 'animateloop(\'' + name + '\', ' + length + ', 0, 1, \'' + col + '\')';
                setTimeout(command , 100);
        }

        function animateloop(name, length, ind, delta, col) {
		var next = ind + delta;
		if (next >= length) { delta = delta * -1; next = ind + delta; }
		if (next < 0) { delta = delta * -1; next = ind + delta; }
                setColor(name + '_' + ind, getColor(name + '_' + next));
                setColor(name + '_' + next, col);
                command = 'animateloop(\'' + name + '\', ' + length + ', ' + next + ', ' + delta + ', \'' + col + '\')';
                setTimeout(command , 100);
        }
