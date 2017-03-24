function show_dialog(id, title, html, width, height) {
	var d = $('#'+id);
	if (title) { d.dialog("option", "title", title); }
	if (html) { d.html(html); }
	if (width) {
		d.dialog('option', 'width', width);
		d.dialog('option', 'minWidth', width);
	}
	else
		width = d.dialog('option', 'width');
	if (height) { d.dialog('option', 'height', height); }
	else { height = d.dialog('option', 'height') };
	var xpos = $(window).width()/2 - width/2;
	var ypos = 100;//$(window).height()/2 - height/2; // hardcode height because jquery not correctly reporting $(window).height()
	d.dialog('option', 'position', [xpos, 100]);
	d.dialog('open');
}

function edit_notebook_info() {
	$.ajax({
		data: {
			fname: 'edit_notebook_info',
			lid: NOTEBOOK_ID,
		},
		success : function(data) {
			var obj = JSON.parse(data);
			show_dialog('notebook_info_edit_box', '', obj.output, '31em');
		},
	});
}

function update_notebook_info() {
	var name = $('#edit_name').val();
	if (!name) {
		alert('Please specify a name.');
		return;
	}

	var desc = $('#edit_desc').val();
	//var type = $('#edit_type').val(); // mdb removed 12/14/16 COGE-800

	$.ajax({
		data: {
			fname: 'update_notebook_info',
			lid: NOTEBOOK_ID,
			name: name,
			desc: desc,
			//type: type // mdb removed 12/14/16 COGE-800
		},
		success : function(val) {
			get_notebook_info();
			$("#notebook_info_edit_box").dialog('close');
		},
	});
}

function get_notebook_info() {
	$.ajax({
		data: {
			fname: 'get_notebook_info',
			nid: NOTEBOOK_ID
		},
		success : function (data) {
			$('#notebook_info').html(data);
		}
	});
}

function make_notebook_public(public) {
    if (typeof public === 'undefined')
        public = true;
	coge.services.update('notebook', NOTEBOOK_ID, {metadata: {restricted: !public}}).done(get_notebook_info);
}

function add_list_items(opts) {
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

function add_selected_items(select_id){
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

function remove_list_items() {
	var selected_rows = contents.grid.getSelectedRows();
	var item_list = contents.grid.getSelectedItemList();
	if (item_list) {
	    contents.busy();
        $.ajax({
            data: {
                fname: 'remove_list_items',
                lid: NOTEBOOK_ID,
                item_list: item_list
            },
            success : function(val) {
                get_list_contents().done( function() { contents.busy(false) } );
            },
        });
	}
}

function get_list_contents() {
	return coge.services.fetch_notebook(NOTEBOOK_ID)
		.done(function(response) { // success
		    if (response && response.items) {
		        // Update contents table
		        contents.set(response.items);

		        // Filter unsupported types in SendTo menu
		        if (sendTo) {
		            // Create a list of unique types of items
		            var types = response.items
                        .map(function(item) { return item.type})
                        .filter(function(value, index, self) {
                            return self.indexOf(value) === index;
                        });
		            sendTo.filter(types);
		        }
            }
        })
		.fail(function() { // error
			//TODO
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

function search_mystuff() { //TODO migrate to web services
	var search_term = $('#edit_mystuff_search').val();

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

function search_genomes() { //TODO migrate to web services
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

function search_experiments(search_term) { //TODO migrate to web services
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

function search_features() { //TODO migrate to web services
	var search_term = $('#edit_feature_search').val();

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

function search_users(search_term) {
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

function search_lists() { //TODO migrate to web services
	var search_term = $('#edit_list_search').val();

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

function delete_list() {
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

function snp_merge(method) {
        var experiments = contents.grid.dataTable.api().rows().data() //TODO add getRows() to DataGrid
            .filter(function(item) { return (item.data_type === 'polymorphism') })
            .map(function(item) { return item.id })
            .toArray();

        if (typeof experiments === 'undefined' || experiments.length == 0) {
            alert("No SNP experiments exist in this notebook");
            return;
        }

		coge.progress.begin({
			title: "Merging (" + method + " method) SNP experiments ...",
			width: '60%',
	    	height: '50%'
		});

		// Build request
		var request = {
			type: 'merge_snps',
			requester: {
				page: PAGE_NAME,
				url: PAGE_NAME + "?nid=" + NOTEBOOK_ID
			},
			parameters: {
			    method: method,
				experiments: experiments,
				metadata: {
				    name: "test" // FIXME
				}
			}
		};

		// Submit request
		coge.services.submit_job(request)
			.done(function(response) {
		  		if (!response) {
		  			coge.progress.failed("Error: empty response from server");
		  			return;
		  		}
		  		else if (!response.success || !response.id) {
		  			coge.progress.failed("Error: failed to start workflow");
		  			return;
		  		}

		        // Start status update
		  		window.history.pushState({}, "Title", PAGE_NAME + "?nid=" + NOTEBOOK_ID + "&wid=" + response.id); // Add workflow id to browser URL
		  		coge.progress.update(response.id, response.site_url);
		    })
		    .fail(function(jqXHR, textStatus, errorThrown) {
		    	coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
		    });
}

class SendToMenu {
    constructor(params) {
        if (!params.elementId) {
	        console.error('SendToMenu: please specify target element');
	        return;
        }
        this.elementId  = params.elementId;
        this.element = $("#" + this.elementId);
        this.width = params.width || '10em';

        this.render();
    }

    render() {
        var self = this;

        var menu = this.element.menu().css({ position: 'absolute', width: self.width });
        menu
            .position({
                my: "left top",
                at: "left bottom",
                of: "#send_button"
            });
        menu.on("mouseleave", function() { menu.hide(); } );

        // Setup click events on menu options
        menu.children('li')
            .addClass("link")
            .unbind()
            .click(function() {
                var action = $(this).data("action");
                self.submit(action);
                self.hide();
            });
    }

    toggle() {
        if (this.element.is(":visible"))
            this.hide();
        else
            this.show();
    }

    show() {
        this.element.show();
    }

    hide() {
        this.element.hide();
    }

    filter(types) { // array of string types, e.g. ['genome', 'experiment']
        var self = this;

        if (!(types instanceof Array)) {
            console.error('SendTo.filter requires an array argument');
            return;
        }

        if (types && types.length) {
            self.element.children('li').hide();
            self.element.children('li[data-type=all]').show();
            types.forEach(function(type) {
                self.element.children('li[data-type=' + type + ']').show();
            });
        }
    }

    submit(fname) {
        $.ajax({
            dataType: 'json',
            data: {
                fname: fname,
                lid: NOTEBOOK_ID
            },
            success : function(response) {
                if (response && response.url)
                    window.open(response.url);
            }
        });
    }
}

class NotebookContents { // based on ContentPanel, perhaps we should just use that class instead
	constructor(params) {
	    var self = this;

	    if (!params.elementId) {
	        console.error('NotebookContents: please specify target element');
	        return;
        }
        this.elementId = params.elementId;
        this.title     = params.title || "Contents";
        this.height    = params.height || 400;
        this.titleElementId = params.titleElementId;
        this.buttonPanelElementId = params.buttonPanelElementId;

        this.grid = new DataGrid({
			elementId: self.elementId,
			height:    self.height,
			options: { // dataTable options
                language: {
                    emptyTable: '<i>This notebook is empty.  Click the "+" button to add genomes, experiments, or features.</i>'
                }
			},
			columns: [
	            { 	title: "Type",
	            	targets: 0,
	            	//orderable: false,
	            	type: "string",
	            	data: null, // use full data object
	            	width: "5em",
	            	render: function(data, type, row, meta) {
	            		return data.type;
	            	}
	            },
	            { 	title: "Name",
	            	targets: 1,
	            	type: "html",
	            	data: null, // use full data object
	            	render: function(data, type, row, meta) {
	            	    var url;
	            	    if (data.type == 'genome') url = "GenomeInfo.pl?gid=" + data.id;
	            	    else if (data.type == 'experiment') url = "ExperimentView.pl?eid=" + data.id;
	            	    else if (data.type == 'feature') url = "FeatView.pl?fid=" + data.id;
	            		return "<a href='" + url + "' target='_blank'>" + data.info + "</a>";
	            	}
	            },
	            { 	title: "Date added",
	            	targets: 2,
	            	type: "date",
	            	data: null, // use full data object
					orderSequence: [ 'desc', 'asc' ],
	            	width: "10em",
	            	render: function(data, type, row, meta) {
	            		return ("date" in data ? data.date : '');
	            	}
	            },
			],
			selectionCallback: function(selectedItems) {
                self.renderButtons(selectedItems);
			}
		});
    }

	renderButtons(selectedItems) {
	    if (this.buttonPanelElementId) {
	        var buttonPanel = $('#'+this.buttonPanelElementId);

            var deleteBtn = buttonPanel.find('#delete_button');
            if (selectedItems && selectedItems.length > 0)
                deleteBtn.removeClass('coge-disabled');
            else
                deleteBtn.addClass('coge-disabled');

            var sendBtn = buttonPanel.find('#send_button');
            if (this.grid.getNumRows() > 0)
                sendBtn.removeClass('coge-disabled');
            else
                sendBtn.addClass('coge-disabled');
        }
	}

    renderTitle() {
        if (this.titleElementId) {
            var title = this.title + '&nbsp;&nbsp;<span class="small info">[' + this.grid.getNumRowsDisplayed() + ']</span>';
            $('#'+this.titleElementId).html(title);
        }
    }

	busy(on) {
	    if (typeof on === 'undefined' || on)
            $('#'+this.elementId).addClass('coge-disabled');
        else
            $('#'+this.elementId).removeClass('coge-disabled');
	}

	show(visible) {
	    if (typeof visible === "undefined" || visible)
	        $('#'+this.elementId).removeClass("invisible").show();
        else
            $('#'+this.elementId).hide();
	}

	set(data) {
		this.grid.update(data).redraw();
		this.renderButtons();
		this.renderTitle();
	}
}
