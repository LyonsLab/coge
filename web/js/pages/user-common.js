/*
 * user-common.js
 *
 * Functions in common between user.js and search-results.js
 *
 * TODO: These functions will eventually be transformed into classes/widgets.
 *
 */

var timers = [];
function wait_to_search (search_func, search_term) { //TODO move into module
	if (!search_term || search_term.length > 2) {
		pageObj.search_term = search_term;
		if (timers['search']) {
			clearTimeout(timers['search']);
		}

		timers['search'] = window.setTimeout(
			function() {
				search_func(pageObj.search_term);
			},
			500
		);
	}
}

function update_icons(items) { //TODO move into ContentPanel
	if ( items && items.length > 0) 
		$('.item-button:not(#add_button)').removeClass('coge-disabled');
	else
		$('.item-button:not(#add_button)').addClass('coge-disabled');
}

function open_item(url, title) {
	var selected = contentPanel.grid.getSelectedItems();
	if (selected && selected.length == 1)
		selected[0].open(url, title);
    else
        new DataGridRow({}, 'null').open(url, title); // kludge
}

function favorite_items() {
	var selected_rows = contentPanel.grid.getSelectedRows();
	var item_list = contentPanel.grid.getSelectedItemList();
	if (item_list) {
		// Optimistically toggle favorite in UI
		selected_rows.every(function() {
			var d = this.data();
			d.favorite = (d.favorite == '0' ? '1' : '0');
			this.data(d);
		});
		contentPanel.grid.redraw();
		infoPanel.update(null);
		
		// Toggle favorite on server
		$.ajax({
		    url: 'User.pl',
			data: {
				fname: 'favorite_items',
				item_list: item_list
			},
			success : function(data) {
				// TODO handle error
			}
		});
	}
}

function delete_items() {
	var selected_rows = contentPanel.grid.getSelectedRows();
	var item_list = contentPanel.grid.getSelectedItemList();
	if (item_list) {
		$.ajax({
		    url: 'User.pl',
			data: {
				fname: 'delete_items',
				item_list: item_list
			},
			success : function(data) {
				selected_rows.every(function() {
					var d = this.data();
					d.deleted = '1';
					this.data(d);
				});
				contentPanel.grid.redraw();
				infoPanel.update(null);
			}
		});
	}
}

function undelete_items() {
	var selected_rows = contentPanel.grid.getSelectedRows();
	var item_list = contentPanel.grid.getSelectedItemList();
	if (item_list) {
		$.ajax({
		    url: 'User.pl',
			data: {
				fname: 'undelete_items',
				item_list: item_list
			},
			success : function(data) {
				selected_rows.every(function() {
					var d = this.data();
					d.deleted = '0';
					this.data(d);
				});
				contentPanel.grid.redraw();
				infoPanel.update(null);
			}
		});
	}
}

function add_to_notebook_dialog() {
	var selected = contentPanel.grid.getSelectedItems();
	if (selected.length) {
		// Initialize add-to-notebook dialog
		window.setTimeout(search_notebooks, 1000);
		
		// Open dialog window
		$('#add_to_notebook_dialog').dialog({width:500}).dialog('open');
	}
}

function search_notebooks () {
	var search_term = $('#notebook_search_input').val();

	$("#wait_notebook").animate({opacity:1});
	$("#notebook_select").html("<option disabled='disabled'>Searching...</option>");

	coge.utils.ajaxWithTimestamp({
	    url: 'User.pl',
		data: {
			fname: 'search_notebooks',
			search_term: search_term,
		},
	}).pipe(function(data) {
	    console.log('search_notebooks');
	    console.log(data);
        if (data) {
            var obj = jQuery.parseJSON(data);
            $("#notebook_select").html(obj.html);
            $("#wait_notebook").animate({opacity:0});
        }
    });
}

function add_items_to_notebook() {
	var item_list = contentPanel.grid.getSelectedItemList();
	var nid = $('#notebook_select').find('option:selected').val();
	if (nid && item_list) {
		$.ajax({
		    url: 'User.pl',
			data: {
				fname: 'add_items_to_notebook',
				nid: nid,
				item_list: item_list,
			},
			success : function(data) {
				$('#add_to_notebook_dialog').dialog('close');
			}
		});
	}
	else {
		$('#add_to_notebook_dialog').dialog('close');
	}
}

function share_dialog() {
	var item_list = contentPanel.grid.getSelectedItemList();
	if (item_list) {
		$.ajax({
		    url: 'User.pl',
			data: {
				fname: 'get_share_dialog',
				item_list: item_list,
			},
			success : function(data) {
				if (data)
					$('#share_dialog').html(data).dialog({width:'31em'}).dialog('open');
			}
		});
	}
}

function make_items_public(make_public) {
	var item_list = contentPanel.grid.getSelectedItemList();
	if (item_list) {
		$.ajax({
		    url: 'User.pl',
			data: {
				fname: 'make_items_public',
				item_list: item_list,
				make_public: make_public
			},
			success : function(data) {
				if (data)
					$('#share_dialog').html(data);
			}
		});
	}
}

function remove_items_from_user_or_group(target_item) {
	var item_list = contentPanel.grid.getSelectedItemList();
	if (target_item && item_list) {
		$.ajax({
		    url: 'User.pl',
			data: {
				fname: 'remove_items_from_user_or_group',
				target_item: target_item,
				item_list: item_list,
			},
			success : function(data) {
				if (data)
					$('#share_dialog').html(data);
			}
		});
	}
}

function add_items_to_user_or_group() {
	var item_list = contentPanel.grid.getSelectedItemList();
	var target_item = $('#share_input').data('select_id');
	var role_id = $('#share_role_select').val();
	if (target_item && item_list) {
		$.ajax({
		    url: 'User.pl',
			data: {
				fname: 'add_items_to_user_or_group',
				target_item: target_item,
				role_id: role_id,
				item_list: item_list,
			},
			success : function(data) {
				if (data) 
					$('#share_dialog').html(data);
			}
		});
	}
}

function search_share () {
	var search_term = $('#share_input').val();
	//$("#wait_notebook").animate({opacity:1});

	coge.utils.ajaxWithTimestamp({
	    url: 'User.pl',
		data: {
			fname: 'search_share',
			search_term: search_term,
		}
	}).pipe(function(data) {
        var obj = jQuery.parseJSON(data);
        if (obj && obj.items)
            $("#share_input").autocomplete({source: obj.items}).autocomplete("search");
    });
}

function search_group () { // FIXME dup of above routine but for group dialog
	var search_term = $('#group_input').val();

	coge.utils.ajaxWithTimestamp({
	    url: 'User.pl',
		data: {
			fname: 'search_share',
			search_term: search_term,
		}
	}).pipe(function(data) {
        var obj = jQuery.parseJSON(data);
        if (obj && obj.items)
            $("#group_input").autocomplete({source: obj.items}).autocomplete("search");
    });
}

function group_dialog() {
	var item_list = contentPanel.grid.getSelectedItemList();
	if (item_list) {
		$.ajax({
		    url: 'User.pl',
			data: {
				fname: 'get_group_dialog',
				item_list: item_list,
			},
			success : function(data) {
				if (data)
					$('#group_dialog').html(data).dialog({width:500}).dialog('open');
			}
		});
	}
}

function change_group_role() {
	var item_list = contentPanel.grid.getSelectedItemList();
	var role_id = $('#group_role_select').val();
	if (role_id && item_list) {
		$.ajax({
		    url: 'User.pl',
			data: {
				fname: 'change_group_role',
				target_items: item_list,
				role_id: role_id,
			},
			success : function(data) {
				if (data)
					$('#group_dialog').html(data);
			}
		});
	}
}

function add_users_to_group() {
	var item_list = contentPanel.grid.getSelectedItemList();
	var new_item = $('#group_input').data('select_id');
	if (new_item && item_list) {
		$.ajax({
		    url: 'User.pl',
			data: {
				fname: 'add_users_to_group',
				target_items: item_list,
				new_item: new_item,
			},
			success : function(data) {
				if (data) 
					$('#group_dialog').html(data);
			}
		});
	}
}

function remove_user_from_group(user_id) {
	var item_list = contentPanel.grid.getSelectedItemList();
	if (user_id && item_list) {
		$.ajax({
		    url: 'User.pl',
			data: {
				fname: 'remove_user_from_group',
				target_items: item_list,
				user_id: user_id,
			},
			success : function(data) {
				if (data) 
					$('#group_dialog').html(data);
			}
		});
	}
}

function edit_dialog() {
	if (contentPanel.getView() == 'group') {
		group_dialog();
	}
//	else if (contentPanel.getView() == 'notebook') {
//		add_to_notebook_dialog();
//	}
}

function send_menu() {
	var menu = $("#send_menu");

	// Positioning is done here instead of onload to work-around misplacement problem
	// due to contents title missing on page load.
	if (!pageObj.positionMenu) {
		menu.position({
			my: "left top",
			at: "left bottom",
			of: "#send_button"
		});
		pageObj.positionMenu = 1;
	}

	if (menu.is(":visible")) {
		menu.hide();
	}
	else {
		menu.show();
		menu.one("mouseleave", function() { menu.hide(); } );
	}
}

function send_items_to(page_name, format) {
	var item_list = contentPanel.grid.getSelectedItemList();
	if (item_list) {
		$.ajax({
		    url: 'User.pl',
			data: {
				fname: 'send_items_to',
				page_name: page_name,
				format: format,
				item_list: item_list,
			},
			success : function(url) {
				if (url)
					window.open(url);
			}
		});
	}
}

function toggle_star(img, id) {
	$.ajax({
	    url: 'User.pl',
		data: {
			fname: 'toggle_star',
			log_id: id,
		},
		success :  function(val) {
			$(img).attr({ src: (val == 0 ? "picts/star-hollow.png" : "picts/star-full.png") });
		}
	});
}