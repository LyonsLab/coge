//Global Variables
var user_is_admin = false;

$(function () {
	// Initialize CoGe web services
    coge.services.init({
    	baseUrl: API_BASE_URL,
    	userName: USER_NAME
    });
    
    // See if the current user is an admin
    $.ajax({
		data: {
			fname: 'user_is_admin',
		},
		success: function(data) {
			user_is_admin = data == 1;
		}
	}).done(function() {
		search_stuff(SEARCH_TEXT);
	});
});

function search_stuff(search_term) {
	if (!search_term || search_term.length <= 2) {
		$("#noresult").html('Please specify a search term longer than 2 characters').show();
		$("#loading").hide();
		return;
	}
	
	$("#noresult").hide();
	$('masterTable').css('display', 'none');
	$("#loading").show();
	
	coge.services.search_global(search_term)
		.done(function(response) {
			var obj = response;
			
			if (!obj || !obj.results)
				return;

			function add_item(url, obj, num) {
				let html = '<tr class="';
				html += num % 2 ? 'odd' : 'even';
				html += '"><td>';
				if (url) {
					html += '<a target="_blank" href="';
					html += url;
					html += obj.id;
					html += '"';
				} 
				else
					html += '<span';
				if (obj.description) {
					html += ' title=';
					html += JSON.stringify(obj.description);
				}
				if (obj.deleted == 1)
					html += ' style="color:red;"';
				html += '>';
				if (obj.favorite)
					html += '&#11088;&nbsp;&nbsp;';
				if (obj.certified)
					html += '&#x2705;&nbsp;&nbsp;';
				if (obj.restricted)
					html += '&#x1f512;&nbsp;&nbsp;';
				html += obj.name;
				html += url ? '</a>' : '</span>';
				if (obj.type == 'genome')
					html += '</td><td><span class="coge-button" onclick="window.open(&quot;GenomeView.pl?embed=0&amp;gid=' + obj.id + '&quot;)">Browse</span>';
				else if (obj.type == 'feature')
					html += '</td><td>' + obj.feature_type + '</td><td>' + obj.genome;
				html += '</td><td>';
				html += obj.id;
				html += '</td></tr>';
				return html;
			}
			
			var org_count = 0, gen_count = 0, exp_count = 0, feature_count = 0, note_count = 0, user_group_count = 0;
			var org_list = "", gen_list = "", exp_list = "", feature_list = "", note_list = "", user_group_list = "";

			for (var i = 0; i < obj.results.length; i++) {
				var o = obj.results[i];

				if (o.type == "organism")
					org_list += add_item('OrganismView.pl?oid=', o, org_count++);
				else if (o.type == "genome")
					gen_list += add_item('GenomeInfo.pl?gid=', o, gen_count++);
				else if (o.type == "experiment")
					exp_list += add_item('ExperimentView.pl?eid=', o, exp_count++);
				else if (o.type == "feature")
					feature_list += add_item('FeatView.pl?fid=', o, feature_count++);
				else if (o.type == "notebook")
					note_list += add_item('NotebookView.pl?lid=', o, note_count++);	
				else if (o.type == "user_group")
					user_group_list += add_item(null, o, user_group_count++);
			}

			//Populate the html with the results
			$("#loading").show();
			$('masterTable').css('display', 'block');
			$(".result").fadeIn( 'fast');
			
			if (org_count + gen_count + exp_count + feature_count + note_count + user_group_count == 0)
				$('#noresult').html('No matching results found').show();

			function setup_results(type, icon, count, cols, rows) {
				var div = $('#' + type.replace(' ', '_'));
				if (count > 0) {
					div.show();
					div.children().first().children().first().html('<img src="picts/' + icon + '" style="width:15px" /> ' + type + 's: ' + count);
					var table = div.find('table');
					var tr = $('<tr></tr>').appendTo($('<thead></thead>').appendTo(table));
					cols.forEach(function(col) {
						tr.append($('<th>' + col + '</th>'))
					});
					table.append($('<tbody>' + rows + '</tbody>'));
				}
			}
			setup_results('Experiment', 'testtube-icon.png', exp_count, ['name', 'id'], exp_list);
			setup_results('Feature', 'feature-icon.png', feature_count, ['name', 'type', 'genome', 'id'], feature_list);
			setup_results('Genome', 'dna-icon.png', gen_count, ['name', 'EPIC-CoGe', 'id'], gen_list);
			setup_results('Notebook', 'notebook-icon.png', note_count, ['name', 'id'], note_list);
			setup_results('Organism', 'Organism.svg', org_count, ['name', 'id'], org_list);
			setup_results('User Group', 'group-icon.png', user_group_count, ['name', 'id'], user_group_list);
			
			$("#loading").hide();
			$("#masterTable").show();
			if (!user_is_admin)
				$(".access").hide();
			else
				$(".access").show();
		})
		.fail(function() {
			$("#loading").hide();
			$("#masterTable").html("An error occured. Please reload the page and try again.");
		});
}

function toggle_results(div, table) {
	div = $(div);
	var img = div.children().last();
	if (img.attr('src') == 'picts/arrow-right-icon.png') {
        img.attr('src', 'picts/arrow-down-icon.png');
		if (!$.fn.dataTable.isDataTable(table))
			table.DataTable({info: false, lengthChange: false, oLanguage: { sSearch: 'Filter:' }, paging: false});
    } else
		img.attr('src', 'picts/arrow-right-icon.png');
	div.next().fadeToggle('fast');
}
