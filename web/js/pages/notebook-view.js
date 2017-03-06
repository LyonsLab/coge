function show_dialog(id, title, html, width, height) {
	var d = $('#'+id);
	if (title) { d.dialog("option", "title", title); }
	if (html) { d.html(html); }
	if (width) {
		d.dialog('option', 'width', width);
		d.dialog('option', 'minWidth', width);
	} else
		width = d.dialog('option', 'width');
	if (height) { d.dialog('option', 'height', height); }
	else { height = d.dialog('option', 'height') };
	var xpos = $(window).width()/2 - width/2;
	var ypos = 100;//$(window).height()/2 - height/2; // hardcode height because jquery not correctly reporting $(window).height()
	d.dialog('option', 'position', [xpos, 100]);
	d.dialog('open');
}

function edit_list_info () {
	$.ajax({
		data: {
			fname: 'edit_list_info',
			lid: NOTEBOOK_ID,
		},
		success : function(data) {
			var obj = JSON.parse(data);
			show_dialog('list_info_edit_box', '', obj.output, '31em');
		},
	});
}

function update_list_info (){
	var name = $('#edit_name').val();
	if (!name) {
		alert('Please specify a name.');
		return;
	}

	var desc = $('#edit_desc').val();
	//var type = $('#edit_type').val(); // mdb removed 12/14/16 COGE-800

	$.ajax({
		data: {
			fname: 'update_list_info',
			lid: NOTEBOOK_ID,
			name: name,
			desc: desc,
			//type: type // mdb removed 12/14/16 COGE-800
		},
		success : function(val) {
			get_list_info();
			$("#list_info_edit_box").dialog('close');
		},
	});
}

function get_list_info() {
	$.ajax({
		data: {
			fname: 'get_list_info',
			lid: NOTEBOOK_ID
		},
		success : function (data) {
			$('#list_info').html(data);
		}
	});
}

function make_list_public () {
	$.ajax({
		data: {
			fname: 'make_list_public',
			lid: NOTEBOOK_ID,
		},
		success : function(val) {
			get_list_info();
		}
	});
}

function make_list_private () {
	$.ajax({
		data: {
			fname: 'make_list_private',
			lid: NOTEBOOK_ID,
		},
		success : function(val) {
			get_list_info();
		},
	});
}

function add_list_items (opts) {
	$.ajax({
		data: {
			fname: 'add_list_items',
			lid: NOTEBOOK_ID,
		},
		success : function(data) {
			var obj = JSON.parse(data);
			show_dialog('list_contents_edit_box', '', obj.output, 650, '600');
		},
	});
}

function add_selected_items (select_id){
	var num_items = $('#' + select_id).find('option:selected').length;
	$('#' + select_id).find('option:selected').each(
		function() {
			var item_spec = $(this).attr("value");
			$.ajax({
				data: {
					fname: 'add_item_to_list',
					lid: NOTEBOOK_ID,
					item_spec : item_spec,
				},
				success :
					function(data) {
						if (data != 1) { alert(data); }
						else {
							if (--num_items == 0) { // only do update on last item
								get_list_contents();
							}
						}
					},
			});
		}
	);
	$("#list_contents_edit_box").dialog('close');
}

function remove_list_item (obj, opts) {
	var item_id = opts.item_id;
	var item_type = opts.item_type;

	$(obj).closest('tr').children().animate({opacity:0});

	$.ajax({
		data: {
			fname: 'remove_list_item',
			lid: NOTEBOOK_ID,
			item_id: item_id,
			item_type: item_type,
		},
		success : function(val) {
			get_list_contents();
		},
	});
}

function get_list_contents() {
	$.ajax({
		data: {
			fname: 'get_list_contents',
			lid: NOTEBOOK_ID,
		},
		success : set_contents
	});
}

function wait_to_search (search_func, search_term) {  //TODO use version in utils
	if (!search_term || search_term.length > 2) {
		pageObj.search_term = search_term;
		if (pageObj.time) {
			clearTimeout(pageObj.time);
		}

		// FIXME: could generalize by passing select id instead of separate search_* functions
		pageObj.time = setTimeout(
			function() {
				search_func(pageObj.search_term);
			},
			500
		);
	}
}

// FIXME: the search functions below are all the same, consolidate?

function search_mystuff () { //TODO migrate to web services
	var search_term = $('#edit_mystuff_search').attr('value');

	$("#wait_mystuff").animate({opacity:1});
	$("#select_mystuff_items").html("<option disabled='disabled'>Searching...</option>");
	pageObj.timestamp['mystuff'] = new Date().getTime();

	$.ajax({
		data: {
			fname: 'search_mystuff',
			lid: NOTEBOOK_ID,
			search_term: search_term,
			timestamp: pageObj.timestamp['mystuff']
		},
		success : function(val) {
			var items = JSON.parse(val);
			if (items.timestamp == pageObj.timestamp['mystuff']) {
				$("#select_mystuff_items").html(items.html);
				$("#wait_mystuff").animate({opacity:0});
			}
		},
	});
}

function search_genomes () { //TODO migrate to web services
	var search_term = $('#edit_genome_search').val();

	$("#wait_genome").animate({opacity:1});
	$("#select_genome_items").html("<option disabled='disabled'>Searching...</option>");
	pageObj.timestamp['genomes'] = new Date().getTime();

	$.ajax({
		data: {
			fname: 'search_genomes',
			lid: NOTEBOOK_ID,
			search_term: search_term,
			timestamp: pageObj.timestamp['genomes']
		},
		success : function(val) {
			var items = JSON.parse(val);
			if (items.timestamp == pageObj.timestamp['genomes']) {
				$("#select_genome_items").html(items.html);
				$("#wait_genome").animate({opacity:0});
			}
		},
	});
}

function search_experiments (search_term) { //TODO migrate to web services
	var search_term = $('#edit_experiment_search').val();

	$("#wait_experiment").animate({opacity:1});
	$("#select_experiments_items").html("<option disabled='disabled'>Searching...</option>");
	pageObj.timestamp['experiments'] = new Date().getTime();

	$.ajax({
		data: {
			fname: 'search_experiments',
			lid: NOTEBOOK_ID,
			search_term: search_term,
			timestamp: pageObj.timestamp['experiments']
		},
		success : function(val) {
			var items = JSON.parse(val);
			if (items.timestamp == pageObj.timestamp['experiments']) {
				$("#select_experiment_items").html(items.html);
				$("#wait_experiment").animate({opacity:0});
			}
		},
	});
}

function search_features () { //TODO migrate to web services
	var search_term = $('#edit_feature_search').attr('value');

	$("#wait_feature").animate({opacity:1});
	$("#select_feature_items").html("<option disabled='disabled'>Searching...</option>");
	pageObj.timestamp['features'] = new Date().getTime();

	$.ajax({
		data: {
			fname: 'search_features',
			lid: NOTEBOOK_ID,
			search_term: search_term,
			timestamp: pageObj.timestamp['features']
		},
		success : function(val) {
			var items = JSON.parse(val);
			if (items.timestamp == pageObj.timestamp['features']) {
				$("#select_feature_items").html(items.html);
				$("#wait_feature").animate({opacity:0});
			}
		},
	});
}

function search_users (search_term) {
    $.ajax({
        data: {
            jquery_ajax: 1,
            fname: 'search_users',
            search_term: search_term,
            timestamp: new Date().getTime()
        },
        success : function(data) {
            var obj = JSON.parse(data);
            if (obj && obj.items) {
                $("#edit_user").autocomplete({source: obj.items});
                $("#edit_user").autocomplete("search");
            }
        },
    });
}

function search_lists () { //TODO migrate to web services
	var search_term = $('#edit_list_search').attr('value');

	$("#wait_list").animate({opacity:1});
	$("#select_list_items").html("<option disabled='disabled'>Searching...</option>");
	pageObj.timestamp['lists'] = new Date().getTime();

	$.ajax({
		data: {
			fname: 'search_lists',
			lid: NOTEBOOK_ID,
			search_term: search_term,
			timestamp: pageObj.timestamp['lists']
		},
		success : function(val) {
			var items = JSON.parse(val);
			if (items.timestamp == pageObj.timestamp['lists']) {
				$("#select_list_items").html(items.html);
				$("#wait_list").animate({opacity:0});
			}
		},
	});
}

function delete_list () {
	$.ajax({
		data: {
			fname: 'delete_list',
			lid: NOTEBOOK_ID
		},
		success : function(val) {
			location.reload();
		},
	});
}

function send_list_to() {
	var action = $('#checked_action').val();

	$.ajax({
		data: {
			fname: action,
			lid: NOTEBOOK_ID
		},
		success : function(val) {
			var items = JSON.parse(val);
			if (items.alert) {
				alert(items.alert);
			}
			if (items.url) {
				window.open(items.url, '_blank');
			}
		}
	});
}

function toggle_favorite(img) {
	$.ajax({
		data: {
			fname: 'toggle_favorite',
			nid: NOTEBOOK_ID
		},
		success :  function(val) {
			$(img).attr({ src: (val == '0' ? "picts/star-hollow.png" : "picts/star-full.png") });
		}
	});
}

class Contents {
	constructor() {
		this.div = $('#list_contents');
	}
	add_row(row, name, page, id_name, node_type, tbody) {
		let tr = $('<tr><td>' + name + '</td></tr>').appendTo(tbody);
		tr.append($('<td><span class="link" onclick="window.open(\'' + page + '.pl?' + id_name + '=' + row[0] + '\')">' + row[1] + '</span></td>'));
		let td = $('<td></td>').appendTo(tr);
		if (row.length > 2 && row[2])
			td.text(row[2]);
		if (this.user_can_edit)
			tr.append($('<td style="padding-left:20px;"><span onClick="remove_list_item(this, {lid: ' + NOTEBOOK_ID + ', item_type: \'' + node_type + '\', item_id: ' + row[0] + '});" class="link ui-icon ui-icon-closethick"></span></td>'));
	}
	build() {
		this.div.empty();
		let table = $('<table id="contents" class="dataTable compact hover stripe border-top border-bottom" style="font-size:14px;margin:0;"></table>').appendTo(this.div);
		let tr = $('<tr><th>Type</th><th>Name</th><th>Date</th></tr>').appendTo($('<thead></thead>').appendTo(table));
		if (this.user_can_edit)
			$('<th>Remove</th>').appendTo(tr);
		let tbody = $('<tbody></tbody>').appendTo(table);
		this.genomes.forEach(function(row) {
			this.add_row(row, 'Genome', 'GenomeInfo', 'gid', 2, tbody);
		}.bind(this));
		this.experiments.forEach(function(row) {
			this.add_row(row, 'Experiment', 'ExperimentView', 'eid', 3, tbody);
		}.bind(this));
		this.features.forEach(function(row) {
			this.add_row(row, 'Feature', 'FeatView', 'fid', 4, tbody);
		}.bind(this));
		let options = {
			columnDefs: [{
				targets: 0,
				render: function(data, type, full, meta) {
					if (type != 'display')
						return data;
//					if (meta.row == 0 || this.type != data) {
//						this.type = data;
						return '<span class="title5" style="font-size:inherit">' + data + '</span>'; //s (' + (data == 'Genome' ? this.genomes.length : data == 'Experiment' ? this.experiments.length : this.features.length) + ')</span>';
//					}
//					return '';
				}.bind(this)
			},{
				targets: 1,
				render: function(data, type, full, meta) {
					if (type == 'sort')
						return data.replace('\uD83D\uDD12 ', ''); // remove lock char
					return data;
				}
			}],
			pageLength: 20
		};
		let num_entries = this.genomes.length + this.experiments.length + this.features.length;
		if (num_entries == 0) {
			options.language.infoEmpty = '';
		}
		if (num_entries <= 20) {
			options.paging = false;
			options.searching = false;
		}
		if (this.user_can_edit)
			options.columnDefs.push({
				targets: 3,
				orderable: false
			})
		table.DataTable(options);
		if (this.user_can_edit)
			this.div.append($('<div class="padded"><span class="coge-button" onClick="add_list_items();"><span class="ui-icon ui-icon-plus"></span>Add</span></div>'));
	}
	set(json) {
		this.experiments = json.experiments;
		this.features = json.features;
		this.genomes = json.genomes;
		this.user_can_edit = json.user_can_edit;
		this.build();
	}
}

var contents = new Contents();
function set_contents(json) {
	if (typeof json == 'string')
		json = JSON.parse(json);
	contents.set(json);
}

function case_insensitive_sort(a, b) {
	let _a = a.toLowerCase();
	if (_a.startsWith('&#x1f512; '))
		_a = _a.substr(6);
	let _b = b.toLowerCase();
	if (_b.startsWith('&#x1f512; '))
		_b = _b.substr(6);
	return _a < _b ? -1 : _a > _b ? 1 : 0;
}
