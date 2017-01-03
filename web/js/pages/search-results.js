$(function () {
	// Initialize CoGe web services
    coge.services.init({
    	baseUrl: API_BASE_URL,
    	userName: USER_NAME
    });

    // Initiate search
    search_stuff(SEARCH_TEXT);

    // Define views in the Content Panel
	var views = {
		organism: {
			title: 'Organisms',
			displayType: 'grid',
			dataTypes: ['organism'],
			operations: ['share', 'organize', 'favorite', 'delete', 'sendto']
		},
		genome: {
			title: 'Genomes',
			displayType: 'grid',
			dataTypes: ['genome'],
			operations: ['share', 'organize', 'favorite', 'delete', 'sendto']
		},
		feature: {
			title: 'Features',
			displayType: 'grid',
			dataTypes: ['feature'],
			operations: ['share', 'organize', 'favorite', 'delete', 'sendto']
		},
		experiment: {
			title: 'Experiments',
			displayType: 'grid',
			dataTypes: ['experiment'],
			operations: ['share', 'organize', 'favorite', 'delete', 'sendto']
		},
		notebook: {
			title: 'Notebooks',
			displayType: 'grid',
			dataTypes: ['notebook'],
			operations: ['share', 'favorite', 'delete', 'sendto', 'add']
		},
		group: {
			title: 'User Groups',
			displayType: 'grid',
			dataTypes: ['group'],
			operations: ['share', 'favorite', 'delete', 'sendto', 'add']
		}
	};

	// Initialize the main panels
	infoPanel = new InfoPanel({
		elementId: 'info_panel',
		defaultInfo: default_info
	});

	contentPanel = new ContentPanel({
		elementId: 'contents_panel',
		views: views,
		grid: new DataGrid({
			element: $('#contents_panel').find('.grid'),
			height: $(window).height() - 210, // this depends on the height of the header/footer and should be passed in as an argument
			columns: [
                { 	title: "",
	            	targets: 0,
	            	orderable: false,
	            	type: "string",
	            	width: 15,
	            	data: null, // use full data object
	            	render: function(data, type, row, meta) {
	            		return data.getFlags();
	            	}
	            },
				{ 	title: "Name",
					targets: 1,
					type: "html",
					data: null, // use full data object
					render: function(data, type, row, meta) {
						return data.getDescription();
					}
				}
			],
			selectionCallback: function(items) {
			    infoPanel.busy().update(items);
				//update_icons(items);
			}
		})
	});

	tocPanel = new TocPanel({
		elementId: 'toc_panel',
		selection: function(typeId) {
			contentPanel
			    .setView(typeId)
			    .render();
			contentPanel.grid.search(''); // clear search filter
			infoPanel.update(null);
			//update_icons(null);
			//$('#search_input').val(''); //FIXME move into ContentPanel
		}
	});
});

function search_stuff(search_term) {
	if (!search_term || search_term.length <= 2) {
		$("#msg").html('Please specify a search term longer than 2 characters').show();
		$("#loading").hide();
		return;
	}
	
	$("#msg,#bottom-panel").hide();
	$("#loading").show();

	coge.services.search_global(search_term)
		.done(function(response) {
			if (!response || !response.results || !response.results.length) {
				$("#loading").hide();
				$("#msg").html('No matching results').show();
				return;
			}

            // Index results by type
            var resultsByType = [];
			for (var i = 0; i < response.results.length; i++) {
                var o = response.results[i];
                if (!resultsByType[o.type])
                    resultsByType[o.type] = [];
                resultsByType[o.type].push(o);
			}
			//console.log(resultsByType);

			for (var type in resultsByType) {
			    contentPanel.setData(type, resultsByType[type]);
			    tocPanel.setCount(type, resultsByType[type].length);
			}

            var firstType = Object.keys(resultsByType)[0];
            if (Object.keys(resultsByType).length == 1) {
                tocPanel.selectItemType(firstType);
                $('#toc_panel').hide();
            }
            else {
                // Remove TOC items with no results
                $('#toc_panel').find('li > span').each(function (index, el) {
                    var type = $(el).attr('data-type');
                    if (!resultsByType.hasOwnProperty(type))
                        $(el).hide();
                });

                // Get starting page from URL and initialize TOC panel
                var toc_id = coge.utils.getURLParameter('p');
                if (toc_id && toc_id != 'null' && resultsByType.hasOwnProperty(toc_id))
                    tocPanel.selectItemType(toc_id);
			    else
			        tocPanel.selectItemType(firstType);
			    $('#toc_panel').show();
            }

			$("#loading").hide();
			$("#bottom-panel").fadeIn();
		})
		.fail(function() {
			$("#loading").hide();
			$("#bottom-panel").html("An error occurred. Please reload the page and try again.");
		});
}

function default_info() {
	switch(contentPanel.getView()) {
	    case 'organism':
        case 'genome':
        case 'feature':
		case 'experiment':
		case 'notebook':
		case 'group':
			return "<p>These are data items that match your search specification.</p>"
				+ "<p><b>Hover over</b> an item to view additional info.</p>"
                + "<p><b>Single-click</b> to select one or more items to " + (USER_ID ? "share, organize, delete, or" : "") + " send to one of CoGe's tools. Use <b>Ctrl-click</b> to select multiple items.</p>"
                + "<p><b>Double-click</b> an item for a detailed view of the item.</p>";
	}
}

function open_item(url) {
	var selected_rows = contentPanel.grid.getSelectedRows();
	if (selected_rows && selected_rows.length == 1) {
		selected_rows.every(function() {
		    console.log('open_item: ' + url);
			this.data().open(url);
		});
    }
}

/*
 * Data Grid Row
 */

class DataGridRow { //FIXME duplicated in user.js
    constructor(data, type) {
        $.extend(this, data);
    }

    getID() {
        return this.id;
    }

    getFlags(opts) {
    	if (this.type == 'genome' ||
    		this.type == 'experiment' ||
    		this.type == 'notebook' ||
    		this.type == 'favorite')
    	{
    		var noSpaces = (opts && opts.noSpaces);

    		var flags = '';
    		if (!noSpaces || this.favorite == '1')
    			flags = '<span style="color:goldenrod;visibility:' +
	    			(this.favorite == '1' ? 'visible' : 'hidden') +
	    			'">&#9733;</span>&nbsp;';
    		if (this.restricted == '1')
	    		flags += '&#x1f512;' + '&nbsp;';

    		return flags;
    	}
    	return '';
    }

    getDescription() {
        if (this.type == 'organism')
            return this.name;
    	if (this.type == 'genome')
    		return this._formatGenome();
        if (this.type == 'feature')
    		return this._formatFeature();
    	if (this.type == 'experiment')
    		return this._formatExperiment();
    	if (this.type == 'notebook')
    		return this._formatNotebook();
    	if (this.type == 'group')
    		return this._formatGroup();
        return this.name;
    }

    _formatGenome() {
    	var certified = '<span class="glyphicon glyphicon-ok coge-certified-icon"></span> <span class="coge-small-text">Certified Genome<span>';
    	var descStr =
    	   	this.name + // this is actually the genome->info
    	   	(this.certified == '1' ? '&nbsp;&nbsp;' + certified : '');
    	return descStr;
    }

    _formatFeature() {
    	var descStr =
    	   	this.name +
    	   	' (' + this.feature_type + ')' +
    	   	' ' + this.genome;
    	return descStr;
    }

    _formatExperiment() {
    	var descStr =
    	   	this.name; // this is actually the experiment->info
    	return descStr;
    }

    _formatNotebook() {
    	var descStr =
    		this.name; // this is actually the list->info
    	return descStr;
    }

    _formatGroup() {
    	var descStr =
    		this.name; // this is actually the group->info
    	return descStr;
    }

    getInfo() {
    	console.log('DataGridRow.getInfo');
    	var self = this;

    	return coge.utils.ajaxWithTimestamp({
    		dataType: 'json',
    		url: 'User.pl',
    		data: {
    			fname: 'get_item_info',
    			item_id: self.id,
    			item_type: self.type,
    		}
    	}).pipe(function(data) {
    	    //console.log('getInfo response');
    	    //console.log(data);
    		if (data)
				return data.html;
    		return;
    	});
    }

    getLink() {
    	var type = this.type;

        if (type == 'organism')
            return 'OrganismView.pl?oid=' + this.id;
    	else if (type == 'genome')
    		return 'GenomeInfo.pl?gid=' + this.id;
    	else if (type == 'experiment')
    		return 'ExperimentView.pl?eid=' + this.id;
    	else if (type == 'notebook')
    		return 'NotebookView.pl?nid=' + this.id;
    	else
    		return this.link;
    }

    open(url) {
        if (this.type == 'group' && !url)
            group_dialog();
        else if (this.type == 'analyses' || this.type == 'loads')
            window.open(this.link, '_blank');
        else {
            var title = this.getDescription();
            var link  = url || this.getLink();
            var flags = this.getFlags({noSpaces: 1});
            title = flags + ' ' + title + "<br><a class='xsmall' style='color:#eeeeee;' href='" + link + "' target='_blank'>[Open in new tab]</a> ";
            link = link + "&embed=1";
            console.log('DataGrid.openItem: ' + link);
            var height = $(window).height() * 0.8;
            var d = $('<div class="dialog_box"><iframe src="'+link+'" height="100%" width="100%" style="border:none;"/></div>')
                .dialog({
                    //title: title,
                    width: '80%',
                    height: height,
                    open: function() { // mdb added 10/16/16 -- fix html in dialog title bar for jQuery 3.1.1 update
                        $(this).prev().find("span.ui-dialog-title").append('<span>'+title+'</span>');
                    }
                })
                .dialog('open');
        }
    }
}
