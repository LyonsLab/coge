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
		elementId: 'info_panel'
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
			    .update(typeId)
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
	
	$("#msg").hide();
	$('masterTable').css('display', 'none');
	$("#loading").show();

	coge.services.search_global(search_term)
		.done(function(response) {
			if (!response || !response.results || !response.results.length) {
				$("#loading").hide();
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
			}

            var firstType = Object.keys(resultsByType)[0];
            if (Object.keys(resultsByType).length == 1) {
                tocPanel.selectItemType(firstType);
                $('#toc_panel').hide();
            }
            else {
			    tocPanel.selectItemType(firstType);
			    $('#toc_panel').show();
            }

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

/*
 * Data Grid Row
 */

class DataGridRow {
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
//    	if (this.type == 'genome' || this.type == 'favorite')
//    		return this._formatGenome();
//    	if (this.type == 'experiment')
//    		return this._formatExperiment();
//    	if (this.type == 'notebook')
//    		return this._formatNotebook();
//    	if (this.type == 'group')
//    		return this._formatGroup();
        return this.name;
    }

    _formatGenome() {
    	var icon = '<img src="picts/dna-icon.png" width="15" height="15" style="vertical-align:middle;"/> ';
    	var certified = '<span class="glyphicon glyphicon-ok coge-certified-icon"></span> <span class="coge-small-text">Certified Genome<span>';
    	var descStr =
    		icon +
    	   	(this.organism ? this.organism : '') +
    	   	(this.name ? ' (' + this.name + ')' : '') +
    	   	(this.description ? ': ' + this.description : '') +
    	   	' (v' + this.version + ', id' + this.id + ')' +
    	   	(this.certified == '1' ? '&nbsp;&nbsp;' + certified : '');
    	return descStr;
    }

    _formatExperiment() {
    	var descStr =
    		'<img src="picts/testtube-icon.png" width="15" height="15" style="vertical-align:middle;"/> ' +
    	   	this.name +
    	   	(this.description ? ': ' + this.description : '') +
    	   	' (v' + this.version + ', id' + this.id + ')';
    	return descStr;
    }

    _formatNotebook() {
    	var descStr =
    		'<img src="picts/notebook-icon.png" width="15" height="15" style="vertical-align:middle;"/> ' +
    		this.name +
    		(this.description ? ': ' + this.description : '') +
    		(this.type_name ? ' (' + this.type_name + ')' : '');
    	return descStr;
    }

    _formatGroup() {
    	var descStr =
    		'<img src="picts/group-icon.png" width="15" height="15" style="vertical-align:middle;"/> ' +
    		this.name +
    		(this.description ? ': ' + this.description : '');;
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
    	    console.log('getInfo response');
    	    console.log(data);
    		if (data)
				return data.html;
    		return;
    	});
    }

    getLink() {
    	var type = this.type;

    	if (type == 'genome')
    		return 'GenomeInfo.pl?gid=' + this.id;
    	else if (type == 'experiment')
    		return 'ExperimentView.pl?eid=' + this.id;
    	else if (type == 'notebook')
    		return 'NotebookView.pl?nid=' + this.id;
    	else
    		return this.link;
    }
}
