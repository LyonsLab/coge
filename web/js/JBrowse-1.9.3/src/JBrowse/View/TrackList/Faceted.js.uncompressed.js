require({cache:{
'dijit/layout/AccordionContainer':function(){
require({cache:{
'url:dijit/layout/templates/AccordionButton.html':"<div data-dojo-attach-event='onclick:_onTitleClick' class='dijitAccordionTitle' role=\"presentation\">\n\t<div data-dojo-attach-point='titleNode,focusNode' data-dojo-attach-event='onkeypress:_onTitleKeyPress'\n\t\t\tclass='dijitAccordionTitleFocus' role=\"tab\" aria-expanded=\"false\"\n\t\t><span class='dijitInline dijitAccordionArrow' role=\"presentation\"></span\n\t\t><span class='arrowTextUp' role=\"presentation\">+</span\n\t\t><span class='arrowTextDown' role=\"presentation\">-</span\n\t\t><img src=\"${_blankGif}\" alt=\"\" class=\"dijitIcon\" data-dojo-attach-point='iconNode' style=\"vertical-align: middle\" role=\"presentation\"/>\n\t\t<span role=\"presentation\" data-dojo-attach-point='titleTextNode' class='dijitAccordionText'></span>\n\t</div>\n</div>\n"}});
define("dijit/layout/AccordionContainer", [
	"require",
	"dojo/_base/array", // array.forEach array.map
	"dojo/_base/declare", // declare
	"dojo/_base/event", // event.stop
	"dojo/_base/fx", // fx.Animation
	"dojo/dom", // dom.setSelectable
	"dojo/dom-attr", // domAttr.attr
	"dojo/dom-class", // domClass.remove
	"dojo/dom-construct", // domConstruct.place
	"dojo/dom-geometry",
	"dojo/keys", // keys
	"dojo/_base/lang", // lang.getObject lang.hitch
	"dojo/sniff", // has("ie") has("dijit-legacy-requires")
	"dojo/topic", // publish
	"../focus",			// focus.focus()
	"../_base/manager",	// manager.defaultDuration
	"dojo/ready",
	"../_Widget",
	"../_Container",
	"../_TemplatedMixin",
	"../_CssStateMixin",
	"./StackContainer",
	"./ContentPane",
	"dojo/text!./templates/AccordionButton.html"
], function(require, array, declare, event, fx, dom, domAttr, domClass, domConstruct, domGeometry,
			keys, lang, has, topic, focus, manager, ready,
			_Widget, _Container, _TemplatedMixin, _CssStateMixin, StackContainer, ContentPane, template){

	// module:
	//		dijit/layout/AccordionContainer


	// Design notes:
	//
	// An AccordionContainer is a StackContainer, but each child (typically ContentPane)
	// is wrapped in a _AccordionInnerContainer.   This is hidden from the caller.
	//
	// The resulting markup will look like:
	//
	//	<div class=dijitAccordionContainer>
	//		<div class=dijitAccordionInnerContainer>	(one pane)
	//				<div class=dijitAccordionTitle>		(title bar) ... </div>
	//				<div class=dijtAccordionChildWrapper>   (content pane) </div>
	//		</div>
	//	</div>
	//
	// Normally the dijtAccordionChildWrapper is hidden for all but one child (the shown
	// child), so the space for the content pane is all the title bars + the one dijtAccordionChildWrapper,
	// which on claro has a 1px border plus a 2px bottom margin.
	//
	// During animation there are two dijtAccordionChildWrapper's shown, so we need
	// to compensate for that.


	var AccordionButton = declare("dijit.layout._AccordionButton", [_Widget, _TemplatedMixin, _CssStateMixin], {
		// summary:
		//		The title bar to click to open up an accordion pane.
		//		Internal widget used by AccordionContainer.
		// tags:
		//		private

		templateString: template,

		// label: String
		//		Title of the pane
		label: "",
		_setLabelAttr: {node: "titleTextNode", type: "innerHTML" },

		// title: String
		//		Tooltip that appears on hover
		title: "",
		_setTitleAttr: {node: "titleTextNode", type: "attribute", attribute: "title"},

		// iconClassAttr: String
		//		CSS class for icon to left of label
		iconClassAttr: "",
		_setIconClassAttr: { node: "iconNode", type: "class" },

		baseClass: "dijitAccordionTitle",

		getParent: function(){
			// summary:
			//		Returns the AccordionContainer parent.
			// tags:
			//		private
			return this.parent;
		},

		buildRendering: function(){
			this.inherited(arguments);
			var titleTextNodeId = this.id.replace(' ','_');
			domAttr.set(this.titleTextNode, "id", titleTextNodeId+"_title");
			this.focusNode.setAttribute("aria-labelledby", domAttr.get(this.titleTextNode, "id"));
			dom.setSelectable(this.domNode, false);
		},

		getTitleHeight: function(){
			// summary:
			//		Returns the height of the title dom node.
			return domGeometry.getMarginSize(this.domNode).h;	// Integer
		},

		// TODO: maybe the parent should set these methods directly rather than forcing the code
		// into the button widget?
		_onTitleClick: function(){
			// summary:
			//		Callback when someone clicks my title.
			var parent = this.getParent();
				parent.selectChild(this.contentWidget, true);
				focus.focus(this.focusNode);
		},

		_onTitleKeyPress: function(/*Event*/ evt){
			return this.getParent()._onKeyPress(evt, this.contentWidget);
		},

		_setSelectedAttr: function(/*Boolean*/ isSelected){
			this._set("selected", isSelected);
			this.focusNode.setAttribute("aria-expanded", isSelected ? "true" : "false");
			this.focusNode.setAttribute("aria-selected", isSelected ? "true" : "false");
			this.focusNode.setAttribute("tabIndex", isSelected ? "0" : "-1");
		}
	});

	var AccordionInnerContainer = declare("dijit.layout._AccordionInnerContainer", [_Widget, _CssStateMixin], {
		// summary:
		//		Internal widget placed as direct child of AccordionContainer.containerNode.
		//		When other widgets are added as children to an AccordionContainer they are wrapped in
		//		this widget.

/*=====
		// buttonWidget: Function|String
		//		Class to use to instantiate title
		//		(Wish we didn't have a separate widget for just the title but maintaining it
		//		for backwards compatibility, is it worth it?)
		 buttonWidget: null,
=====*/

/*=====
		// contentWidget: dijit/_WidgetBase
		//		Pointer to the real child widget
		contentWidget: null,
=====*/

		baseClass: "dijitAccordionInnerContainer",

		// tell nested layout widget that we will take care of sizing
		isLayoutContainer: true,

		buildRendering: function(){
			// Builds a template like:
			//	<div class=dijitAccordionInnerContainer>
			//		Button
			//		<div class=dijitAccordionChildWrapper>
			//			ContentPane
			//		</div>
			//	</div>

			// Create wrapper div, placed where the child is now
			this.domNode = domConstruct.place("<div class='" + this.baseClass +
				"' role='presentation'>", this.contentWidget.domNode, "after");

			// wrapper div's first child is the button widget (ie, the title bar)
			var child = this.contentWidget,
				cls = lang.isString(this.buttonWidget) ? lang.getObject(this.buttonWidget) : this.buttonWidget;
			this.button = child._buttonWidget = (new cls({
				contentWidget: child,
				label: child.title,
				title: child.tooltip,
				dir: child.dir,
				lang: child.lang,
				textDir: child.textDir,
				iconClass: child.iconClass,
				id: child.id + "_button",
				parent: this.parent
			})).placeAt(this.domNode);

			// and then the actual content widget (changing it from prior-sibling to last-child),
			// wrapped by a <div class=dijitAccordionChildWrapper>
			this.containerNode = domConstruct.place("<div class='dijitAccordionChildWrapper' style='display:none'>", this.domNode);
			domConstruct.place(this.contentWidget.domNode, this.containerNode);
		},

		postCreate: function(){
			this.inherited(arguments);

			// Map changes in content widget's title etc. to changes in the button
			var button = this.button;
			this._contentWidgetWatches = [
				this.contentWidget.watch('title', lang.hitch(this, function(name, oldValue, newValue){
					button.set("label", newValue);
				})),
				this.contentWidget.watch('tooltip', lang.hitch(this, function(name, oldValue, newValue){
					button.set("title", newValue);
				})),
				this.contentWidget.watch('iconClass', lang.hitch(this, function(name, oldValue, newValue){
					button.set("iconClass", newValue);
				}))
			];
		},

		_setSelectedAttr: function(/*Boolean*/ isSelected){
			this._set("selected", isSelected);
			this.button.set("selected", isSelected);
			if(isSelected){
				var cw = this.contentWidget;
				if(cw.onSelected){ cw.onSelected(); }
			}
		},

		startup: function(){
			// Called by _Container.addChild()
			this.contentWidget.startup();
		},

		destroy: function(){
			this.button.destroyRecursive();

			array.forEach(this._contentWidgetWatches || [], function(w){ w.unwatch(); });

			delete this.contentWidget._buttonWidget;
			delete this.contentWidget._wrapperWidget;

			this.inherited(arguments);
		},

		destroyDescendants: function(/*Boolean*/ preserveDom){
			// since getChildren isn't working for me, have to code this manually
			this.contentWidget.destroyRecursive(preserveDom);
		}
	});

	var AccordionContainer = declare("dijit.layout.AccordionContainer", StackContainer, {
		// summary:
		//		Holds a set of panes where every pane's title is visible, but only one pane's content is visible at a time,
		//		and switching between panes is visualized by sliding the other panes up/down.
		// example:
		//	|	<div data-dojo-type="dijit/layout/AccordionContainer">
		//	|		<div data-dojo-type="dijit/layout/ContentPane" title="pane 1">
		//	|		</div>
		//	|		<div data-dojo-type="dijit/layout/ContentPane" title="pane 2">
		//	|			<p>This is some text</p>
		//	|		</div>
		//	|	</div>

		// duration: Integer
		//		Amount of time (in ms) it takes to slide panes
		duration: manager.defaultDuration,

		// buttonWidget: [const] String
		//		The name of the widget used to display the title of each pane
		buttonWidget: AccordionButton,

/*=====
		// _verticalSpace: Number
		//		Pixels of space available for the open pane
		//		(my content box size minus the cumulative size of all the title bars)
		_verticalSpace: 0,
=====*/
		baseClass: "dijitAccordionContainer",

		buildRendering: function(){
			this.inherited(arguments);
			this.domNode.style.overflow = "hidden";		// TODO: put this in dijit.css
			this.domNode.setAttribute("role", "tablist");	// TODO: put this in template
		},

		startup: function(){
			if(this._started){ return; }
			this.inherited(arguments);
			if(this.selectedChildWidget){
				this.selectedChildWidget._wrapperWidget.set("selected", true);
			}
		},

		layout: function(){
			// Implement _LayoutWidget.layout() virtual method.
			// Set the height of the open pane based on what room remains.

			var openPane = this.selectedChildWidget;

			if(!openPane){ return;}

			// space taken up by title, plus wrapper div (with border/margin) for open pane
			var wrapperDomNode = openPane._wrapperWidget.domNode,
				wrapperDomNodeMargin = domGeometry.getMarginExtents(wrapperDomNode),
				wrapperDomNodePadBorder = domGeometry.getPadBorderExtents(wrapperDomNode),
				wrapperContainerNode = openPane._wrapperWidget.containerNode,
				wrapperContainerNodeMargin = domGeometry.getMarginExtents(wrapperContainerNode),
				wrapperContainerNodePadBorder = domGeometry.getPadBorderExtents(wrapperContainerNode),
				mySize = this._contentBox;

			// get cumulative height of all the unselected title bars
			var totalCollapsedHeight = 0;
			array.forEach(this.getChildren(), function(child){
	            if(child != openPane){
					// Using domGeometry.getMarginSize() rather than domGeometry.position() since claro has 1px bottom margin
					// to separate accordion panes.  Not sure that works perfectly, it's probably putting a 1px
					// margin below the bottom pane (even though we don't want one).
					totalCollapsedHeight += domGeometry.getMarginSize(child._wrapperWidget.domNode).h;
				}
			});
			this._verticalSpace = mySize.h - totalCollapsedHeight - wrapperDomNodeMargin.h
				- wrapperDomNodePadBorder.h - wrapperContainerNodeMargin.h - wrapperContainerNodePadBorder.h
				- openPane._buttonWidget.getTitleHeight();

			// Memo size to make displayed child
			this._containerContentBox = {
				h: this._verticalSpace,
				w: this._contentBox.w - wrapperDomNodeMargin.w - wrapperDomNodePadBorder.w
					- wrapperContainerNodeMargin.w - wrapperContainerNodePadBorder.w
			};

			if(openPane){
				openPane.resize(this._containerContentBox);
			}
		},

		_setupChild: function(child){
			// Overrides _LayoutWidget._setupChild().
			// Put wrapper widget around the child widget, showing title

			child._wrapperWidget = AccordionInnerContainer({
				contentWidget: child,
				buttonWidget: this.buttonWidget,
				id: child.id + "_wrapper",
				dir: child.dir,
				lang: child.lang,
				textDir: child.textDir,
				parent: this
			});

			this.inherited(arguments);
		},

		addChild: function(/*dijit/_WidgetBase*/ child, /*Integer?*/ insertIndex){
			// Overrides _LayoutWidget.addChild().
			if(this._started){
				// Adding a child to a started Accordion is complicated because children have
				// wrapper widgets.  Default code path (calling this.inherited()) would add
				// the new child inside another child's wrapper.

				// First add in child as a direct child of this AccordionContainer
				var refNode = this.containerNode;
				if(insertIndex && typeof insertIndex == "number"){
					var children = _Widget.prototype.getChildren.call(this);	// get wrapper panes
					if(children && children.length >= insertIndex){
						refNode = children[insertIndex-1].domNode;
						insertIndex = "after";
					}
				}
				domConstruct.place(child.domNode, refNode, insertIndex);

				if(!child._started){
					child.startup();
				}

				// Then stick the wrapper widget around the child widget
				this._setupChild(child);

				// Code below copied from StackContainer
				topic.publish(this.id+"-addChild", child, insertIndex);	// publish
				this.layout();
				if(!this.selectedChildWidget){
					this.selectChild(child);
				}
			}else{
				// We haven't been started yet so just add in the child widget directly,
				// and the wrapper will be created on startup()
				this.inherited(arguments);
			}
		},

		removeChild: function(child){
			// Overrides _LayoutWidget.removeChild().

			// Destroy wrapper widget first, before StackContainer.getChildren() call.
			// Replace wrapper widget with true child widget (ContentPane etc.).
			// This step only happens if the AccordionContainer has been started; otherwise there's no wrapper.
			if(child._wrapperWidget){
				domConstruct.place(child.domNode, child._wrapperWidget.domNode, "after");
				child._wrapperWidget.destroy();
				delete child._wrapperWidget;
			}

			domClass.remove(child.domNode, "dijitHidden");

			this.inherited(arguments);
		},

		getChildren: function(){
			// Overrides _Container.getChildren() to return content panes rather than internal AccordionInnerContainer panes
			return array.map(this.inherited(arguments), function(child){
				return child.declaredClass == "dijit.layout._AccordionInnerContainer" ? child.contentWidget : child;
			}, this);
		},

		destroy: function(){
			if(this._animation){
				this._animation.stop();
			}
			array.forEach(this.getChildren(), function(child){
				// If AccordionContainer has been started, then each child has a wrapper widget which
				// also needs to be destroyed.
				if(child._wrapperWidget){
					child._wrapperWidget.destroy();
				}else{
					child.destroyRecursive();
				}
			});
			this.inherited(arguments);
		},

		_showChild: function(child){
			// Override StackContainer._showChild() to set visibility of _wrapperWidget.containerNode
			child._wrapperWidget.containerNode.style.display="block";
			return this.inherited(arguments);
		},

		_hideChild: function(child){
			// Override StackContainer._showChild() to set visibility of _wrapperWidget.containerNode
			child._wrapperWidget.containerNode.style.display="none";
			this.inherited(arguments);
		},

		_transition: function(/*dijit/_WidgetBase?*/ newWidget, /*dijit/_WidgetBase?*/ oldWidget, /*Boolean*/ animate){
			// Overrides StackContainer._transition() to provide sliding of title bars etc.

			if(has("ie") < 8){
				// workaround animation bugs by not animating; not worth supporting animation for IE6 & 7
				animate = false;
			}

			if(this._animation){
				// there's an in-progress animation.  speedily end it so we can do the newly requested one
				this._animation.stop(true);
				delete this._animation;
			}

			var self = this;

			if(newWidget){
				newWidget._wrapperWidget.set("selected", true);

				var d = this._showChild(newWidget);	// prepare widget to be slid in

				// Size the new widget, in case this is the first time it's being shown,
				// or I have been resized since the last time it was shown.
				// Note that page must be visible for resizing to work.
				if(this.doLayout && newWidget.resize){
					newWidget.resize(this._containerContentBox);
				}
			}

			if(oldWidget){
				oldWidget._wrapperWidget.set("selected", false);
				if(!animate){
					this._hideChild(oldWidget);
				}
			}

			if(animate){
				var newContents = newWidget._wrapperWidget.containerNode,
					oldContents = oldWidget._wrapperWidget.containerNode;

				// During the animation we will be showing two dijitAccordionChildWrapper nodes at once,
				// which on claro takes up 4px extra space (compared to stable AccordionContainer).
				// Have to compensate for that by immediately shrinking the pane being closed.
				var wrapperContainerNode = newWidget._wrapperWidget.containerNode,
					wrapperContainerNodeMargin = domGeometry.getMarginExtents(wrapperContainerNode),
					wrapperContainerNodePadBorder = domGeometry.getPadBorderExtents(wrapperContainerNode),
					animationHeightOverhead = wrapperContainerNodeMargin.h + wrapperContainerNodePadBorder.h;

				oldContents.style.height = (self._verticalSpace - animationHeightOverhead) + "px";

				this._animation = new fx.Animation({
					node: newContents,
					duration: this.duration,
					curve: [1, this._verticalSpace - animationHeightOverhead - 1],
					onAnimate: function(value){
						value = Math.floor(value);	// avoid fractional values
						newContents.style.height = value + "px";
						oldContents.style.height = (self._verticalSpace - animationHeightOverhead - value) + "px";
					},
					onEnd: function(){
						delete self._animation;
						newContents.style.height = "auto";
						oldWidget._wrapperWidget.containerNode.style.display = "none";
						oldContents.style.height = "auto";
						self._hideChild(oldWidget);
					}
				});
				this._animation.onStop = this._animation.onEnd;
				this._animation.play();
			}

			return d;	// If child has an href, promise that fires when the widget has finished loading
		},

		// note: we are treating the container as controller here
		_onKeyPress: function(/*Event*/ e, /*dijit/_WidgetBase*/ fromTitle){
			// summary:
			//		Handle keypress events
			// description:
			//		This is called from a handler on AccordionContainer.domNode
			//		(setup in StackContainer), and is also called directly from
			//		the click handler for accordion labels
			if(this.disabled || e.altKey || !(fromTitle || e.ctrlKey)){
				return;
			}
			var c = e.charOrCode;
			if((fromTitle && (c == keys.LEFT_ARROW || c == keys.UP_ARROW)) ||
					(e.ctrlKey && c == keys.PAGE_UP)){
				this._adjacent(false)._buttonWidget._onTitleClick();
				event.stop(e);
			}else if((fromTitle && (c == keys.RIGHT_ARROW || c == keys.DOWN_ARROW)) ||
					(e.ctrlKey && (c == keys.PAGE_DOWN || c == keys.TAB))){
				this._adjacent(true)._buttonWidget._onTitleClick();
				event.stop(e);
			}
		}
	});

	// Back compat w/1.6, remove for 2.0
	if(has("dijit-legacy-requires")){
		ready(0, function(){
			var requires = ["dijit/layout/AccordionPane"];
			require(requires);	// use indirection so modules not rolled into a build
		});
	}

	// For monkey patching
	AccordionContainer._InnerContainer = AccordionInnerContainer;
	AccordionContainer._Button = AccordionButton;

	return AccordionContainer;
});

},
'dijit/layout/StackContainer':function(){
define("dijit/layout/StackContainer", [
	"dojo/_base/array", // array.forEach array.indexOf array.some
	"dojo/cookie", // cookie
	"dojo/_base/declare", // declare
	"dojo/dom-class", // domClass.add domClass.replace
	"dojo/has",	// has("dijit-legacy-requires")
	"dojo/_base/lang",	// lang.extend
	"dojo/ready",
	"dojo/topic", // publish
	"../registry",	// registry.byId
	"../_WidgetBase",
	"./_LayoutWidget",
	"dojo/i18n!../nls/common"
], function(array, cookie, declare, domClass, has, lang, ready, topic,
			registry, _WidgetBase, _LayoutWidget){

// module:
//		dijit/layout/StackContainer

// Back compat w/1.6, remove for 2.0
if(has("dijit-legacy-requires")){
	ready(0, function(){
		var requires = ["dijit/layout/StackController"];
		require(requires);	// use indirection so modules not rolled into a build
	});
}

var StackContainer = declare("dijit.layout.StackContainer", _LayoutWidget, {
	// summary:
	//		A container that has multiple children, but shows only
	//		one child at a time
	//
	// description:
	//		A container for widgets (ContentPanes, for example) That displays
	//		only one Widget at a time.
	//
	//		Publishes topics [widgetId]-addChild, [widgetId]-removeChild, and [widgetId]-selectChild
	//
	//		Can be base class for container, Wizard, Show, etc.
	//
	//		See `StackContainer.ChildWidgetProperties` for details on the properties that can be set on
	//		children of a `StackContainer`.

	// doLayout: Boolean
	//		If true, change the size of my currently displayed child to match my size
	doLayout: true,

	// persist: Boolean
	//		Remembers the selected child across sessions
	persist: false,

	baseClass: "dijitStackContainer",

/*=====
	// selectedChildWidget: [readonly] dijit._Widget
	//		References the currently selected child widget, if any.
	//		Adjust selected child with selectChild() method.
	selectedChildWidget: null,
=====*/

	buildRendering: function(){
		this.inherited(arguments);
		domClass.add(this.domNode, "dijitLayoutContainer");
		this.containerNode.setAttribute("role", "tabpanel");
	},

	postCreate: function(){
		this.inherited(arguments);
		this.connect(this.domNode, "onkeypress", this._onKeyPress);
	},

	startup: function(){
		if(this._started){ return; }

		var children = this.getChildren();

		// Setup each page panel to be initially hidden
		array.forEach(children, this._setupChild, this);

		// Figure out which child to initially display, defaulting to first one
		if(this.persist){
			this.selectedChildWidget = registry.byId(cookie(this.id + "_selectedChild"));
		}else{
			array.some(children, function(child){
				if(child.selected){
					this.selectedChildWidget = child;
				}
				return child.selected;
			}, this);
		}
		var selected = this.selectedChildWidget;
		if(!selected && children[0]){
			selected = this.selectedChildWidget = children[0];
			selected.selected = true;
		}

		// Publish information about myself so any StackControllers can initialize.
		// This needs to happen before this.inherited(arguments) so that for
		// TabContainer, this._contentBox doesn't include the space for the tab labels.
		topic.publish(this.id+"-startup", {children: children, selected: selected});

		// Startup each child widget, and do initial layout like setting this._contentBox,
		// then calls this.resize() which does the initial sizing on the selected child.
		this.inherited(arguments);
	},

	resize: function(){
		// Overrides _LayoutWidget.resize()
		// Resize is called when we are first made visible (it's called from startup()
		// if we are initially visible). If this is the first time we've been made
		// visible then show our first child.
		if(!this._hasBeenShown){
			this._hasBeenShown = true;
			var selected = this.selectedChildWidget;
			if(selected){
				this._showChild(selected);
			}
		}
		this.inherited(arguments);
	},

	_setupChild: function(/*dijit/_WidgetBase*/ child){
		// Overrides _LayoutWidget._setupChild()

		this.inherited(arguments);

		domClass.replace(child.domNode, "dijitHidden", "dijitVisible");

		// remove the title attribute so it doesn't show up when i hover
		// over a node
		child.domNode.title = "";
	},

	addChild: function(/*dijit/_WidgetBase*/ child, /*Integer?*/ insertIndex){
		// Overrides _Container.addChild() to do layout and publish events

		this.inherited(arguments);

		if(this._started){
			topic.publish(this.id+"-addChild", child, insertIndex);	// publish

			// in case the tab titles have overflowed from one line to two lines
			// (or, if this if first child, from zero lines to one line)
			// TODO: w/ScrollingTabController this is no longer necessary, although
			// ScrollTabController.resize() does need to get called to show/hide
			// the navigation buttons as appropriate, but that's handled in ScrollingTabController.onAddChild().
			// If this is updated to not layout [except for initial child added / last child removed], update
			// "childless startup" test in StackContainer.html to check for no resize event after second addChild()
			this.layout();

			// if this is the first child, then select it
			if(!this.selectedChildWidget){
				this.selectChild(child);
			}
		}
	},

	removeChild: function(/*dijit/_WidgetBase*/ page){
		// Overrides _Container.removeChild() to do layout and publish events

		this.inherited(arguments);

		if(this._started){
			// this will notify any tablists to remove a button; do this first because it may affect sizing
			topic.publish(this.id + "-removeChild", page);	// publish
		}

		// If all our children are being destroyed than don't run the code below (to select another page),
		// because we are deleting every page one by one
		if(this._descendantsBeingDestroyed){ return; }

		// Select new page to display, also updating TabController to show the respective tab.
		// Do this before layout call because it can affect the height of the TabController.
		if(this.selectedChildWidget === page){
			this.selectedChildWidget = undefined;
			if(this._started){
				var children = this.getChildren();
				if(children.length){
					this.selectChild(children[0]);
				}
			}
		}

		if(this._started){
			// In case the tab titles now take up one line instead of two lines
			// (note though that ScrollingTabController never overflows to multiple lines),
			// or the height has changed slightly because of addition/removal of tab which close icon
			this.layout();
		}
	},

	selectChild: function(/*dijit/_WidgetBase|String*/ page, /*Boolean*/ animate){
		// summary:
		//		Show the given widget (which must be one of my children)
		// page:
		//		Reference to child widget or id of child widget

		page = registry.byId(page);

		if(this.selectedChildWidget != page){
			// Deselect old page and select new one
			var d = this._transition(page, this.selectedChildWidget, animate);
			this._set("selectedChildWidget", page);
			topic.publish(this.id+"-selectChild", page);	// publish

			if(this.persist){
				cookie(this.id + "_selectedChild", this.selectedChildWidget.id);
			}
		}

		return d;		// If child has an href, promise that fires when the child's href finishes loading
	},

	_transition: function(newWidget, oldWidget /*===== ,  animate =====*/){
		// summary:
		//		Hide the old widget and display the new widget.
		//		Subclasses should override this.
		// newWidget: dijit/_WidgetBase
		//		The newly selected widget.
		// oldWidget: dijit/_WidgetBase
		//		The previously selected widget.
		// animate: Boolean
		//		Used by AccordionContainer to turn on/off slide effect.
		// tags:
		//		protected extension
		if(oldWidget){
			this._hideChild(oldWidget);
		}
		var d = this._showChild(newWidget);

		// Size the new widget, in case this is the first time it's being shown,
		// or I have been resized since the last time it was shown.
		// Note that page must be visible for resizing to work.
		if(newWidget.resize){
			if(this.doLayout){
				newWidget.resize(this._containerContentBox || this._contentBox);
			}else{
				// the child should pick it's own size but we still need to call resize()
				// (with no arguments) to let the widget lay itself out
				newWidget.resize();
			}
		}

		return d;	// If child has an href, promise that fires when the child's href finishes loading
	},

	_adjacent: function(/*Boolean*/ forward){
		// summary:
		//		Gets the next/previous child widget in this container from the current selection.

		// TODO: remove for 2.0 if this isn't being used.   Otherwise, fix to skip disabled tabs.

		var children = this.getChildren();
		var index = array.indexOf(children, this.selectedChildWidget);
		index += forward ? 1 : children.length - 1;
		return children[ index % children.length ]; // dijit/_WidgetBase
	},

	forward: function(){
		// summary:
		//		Advance to next page.
		return this.selectChild(this._adjacent(true), true);
	},

	back: function(){
		// summary:
		//		Go back to previous page.
		return this.selectChild(this._adjacent(false), true);
	},

	_onKeyPress: function(e){
		topic.publish(this.id+"-containerKeyPress", { e: e, page: this});	// publish
	},

	layout: function(){
		// Implement _LayoutWidget.layout() virtual method.
		var child = this.selectedChildWidget;
		if(child && child.resize){
			if(this.doLayout){
				child.resize(this._containerContentBox || this._contentBox);
			}else{
				child.resize();
			}
		}
	},

	_showChild: function(/*dijit/_WidgetBase*/ page){
		// summary:
		//		Show the specified child by changing it's CSS, and call _onShow()/onShow() so
		//		it can do any updates it needs regarding loading href's etc.
		// returns:
		//		Promise that fires when page has finished showing, or true if there's no href
		var children = this.getChildren();
		page.isFirstChild = (page == children[0]);
		page.isLastChild = (page == children[children.length-1]);
		page._set("selected", true);

		domClass.replace(page.domNode, "dijitVisible", "dijitHidden");

		return (page._onShow && page._onShow()) || true;
	},

	_hideChild: function(/*dijit/_WidgetBase*/ page){
		// summary:
		//		Hide the specified child by changing it's CSS, and call _onHide() so
		//		it's notified.
		page._set("selected", false);
		domClass.replace(page.domNode, "dijitHidden", "dijitVisible");

		page.onHide && page.onHide();
	},

	closeChild: function(/*dijit/_WidgetBase*/ page){
		// summary:
		//		Callback when user clicks the [X] to remove a page.
		//		If onClose() returns true then remove and destroy the child.
		// tags:
		//		private
		var remove = page.onClose(this, page);
		if(remove){
			this.removeChild(page);
			// makes sure we can clean up executeScripts in ContentPane onUnLoad
			page.destroyRecursive();
		}
	},

	destroyDescendants: function(/*Boolean*/ preserveDom){
		this._descendantsBeingDestroyed = true;
		this.selectedChildWidget = undefined;
		array.forEach(this.getChildren(), function(child){
			if(!preserveDom){
				this.removeChild(child);
			}
			child.destroyRecursive(preserveDom);
		}, this);
		this._descendantsBeingDestroyed = false;
	}
});

StackContainer.ChildWidgetProperties = {
	// summary:
	//		These properties can be specified for the children of a StackContainer.

	// selected: Boolean
	//		Specifies that this widget should be the initially displayed pane.
	//		Note: to change the selected child use `dijit/layout/StackContainer.selectChild`
	selected: false,

	// disabled: Boolean
	//		Specifies that the button to select this pane should be disabled.
	//		Doesn't affect programmatic selection of the pane, nor does it deselect the pane if it is currently selected.
	disabled: false,

	// closable: Boolean
	//		True if user can close (destroy) this child, such as (for example) clicking the X on the tab.
	closable: false,

	// iconClass: String
	//		CSS Class specifying icon to use in label associated with this pane.
	iconClass: "dijitNoIcon",

	// showTitle: Boolean
	//		When true, display title of this widget as tab label etc., rather than just using
	//		icon specified in iconClass
	showTitle: true
};

// Since any widget can be specified as a StackContainer child, mix them
// into the base widget class.  (This is a hack, but it's effective.)
// This is for the benefit of the parser.   Remove for 2.0.  Also, hide from doc viewer.
lang.extend(_WidgetBase, /*===== {} || =====*/ StackContainer.ChildWidgetProperties);

return StackContainer;
});

},
'url:dijit/layout/templates/AccordionButton.html':"<div data-dojo-attach-event='onclick:_onTitleClick' class='dijitAccordionTitle' role=\"presentation\">\n\t<div data-dojo-attach-point='titleNode,focusNode' data-dojo-attach-event='onkeypress:_onTitleKeyPress'\n\t\t\tclass='dijitAccordionTitleFocus' role=\"tab\" aria-expanded=\"false\"\n\t\t><span class='dijitInline dijitAccordionArrow' role=\"presentation\"></span\n\t\t><span class='arrowTextUp' role=\"presentation\">+</span\n\t\t><span class='arrowTextDown' role=\"presentation\">-</span\n\t\t><img src=\"${_blankGif}\" alt=\"\" class=\"dijitIcon\" data-dojo-attach-point='iconNode' style=\"vertical-align: middle\" role=\"presentation\"/>\n\t\t<span role=\"presentation\" data-dojo-attach-point='titleTextNode' class='dijitAccordionText'></span>\n\t</div>\n</div>\n",
'dojox/grid/EnhancedGrid':function(){
define("dojox/grid/EnhancedGrid", [
	"dojo/_base/kernel",
	"../main",
	"dojo/_base/declare",
	"dojo/_base/lang",
	"dojo/_base/array",
	"dojo/_base/sniff",
	"dojo/dom",
	"dojo/dom-geometry",
	"./DataGrid",
	"./DataSelection",
	"./enhanced/_PluginManager",
	"./enhanced/plugins/_SelectionPreserver",//default loaded plugin
	"dojo/i18n!./enhanced/nls/EnhancedGrid"
], function(dojo, dojox, declare, lang, array, has, dom, domGeometry,
	DataGrid, DataSelection, _PluginManager, _SelectionPreserver, nls){

dojo.experimental("dojox.grid.EnhancedGrid");

var EnhancedGrid = declare("dojox.grid.EnhancedGrid", DataGrid, {
	// summary:
	//		Provides enhanced features based on DataGrid
	//
	// description:
	//		EnhancedGrid features are implemented as plugins that could be loaded on demand.
	//		Explicit dojo.require() is needed to use these feature plugins.
	//
	// example:
	//		A quick sample to use EnhancedGrid features:
	//
	//	   Step 1. Load EnhancedGrid and required features
	// |   <script type="text/javascript">
	// |		dojo.require("dojox.grid.EnhancedGrid");
	// |		dojo.require("dojox.grid.enhanced.plugins.DnD");
	// |		dojo.require("dojox.grid.enhanced.plugins.Menu");
	// |		dojo.require("dojox.grid.enhanced.plugins.NestedSorting");
	// |		dojo.require("dojox.grid.enhanced.plugins.IndirectSelection");
	// |	</script>
	//
	//		Step 2. Use EnhancedGrid
	//
	//		- Via HTML markup
	//
	// |	<div dojoType="dojox.grid.EnhancedGrid" ...
	// |		plugins="{nestedSorting: true, dnd: true, indirectSelection: true,
	// |		menus:{headerMenu:"headerMenuId", rowMenu:"rowMenuId", cellMenu:"cellMenuId",
	// |		selectedRegionMenu:"selectedRegionMenuId"}}">
	// |			...
	// |	</div>
	//
	//		- Or via JavaScript
	//
	// |	<script type="text/javascript">
	// |		var grid = new dojox.grid.EnhancedGrid({plugins : {nestedSorting: true, dnd: true, indirectSelection: true,
	// |	               menus:{headerMenu:"headerMenuId", rowMenu:"rowMenuId", cellMenu:"cellMenuId",selectedRegionMenu:"selectedRegionMenuId"}},
	// |			       ... }, dojo.byId('gridDiv'));
	// |		grid.startup();
	// |	</script>
	//
	//
	//		Plugin Support
	//		[Note: Plugin support is still experimental]
	//
	//		You can either customize the default plugins or add new ones, more details please see
	//
	//		- dojox.grid.enhanced._PluginManager
	//		- dojox.grid.enhanced._Plugin
	//		- dojox.grid.enhanced.plugins.*

	// plugins: Object
	//		Plugin properties, e.g. {nestedSorting: true, dnd: true, ...}
	plugins: null,

	// pluginMgr: Object
	//		Singleton plugin manager
	pluginMgr: null,

	// _pluginMgrClass: Object
	//		Default plugin manager class
	_pluginMgrClass: _PluginManager,

	postMixInProperties: function(){
		this._nls = nls;
		this.inherited(arguments);
	},
	postCreate: function(){
		//create plugin manager
		this.pluginMgr = new this._pluginMgrClass(this);
		this.pluginMgr.preInit();
		this.inherited(arguments);
		this.pluginMgr.postInit();
	},
	plugin: function(/*String*/name){
		// summary:
		//		An easier way for getting a plugin, e.g. grid.plugin('dnd')
		return this.pluginMgr.getPlugin(name);
	},
	startup: function(){
		this.inherited(arguments);
		this.pluginMgr.startup();
	},
	createSelection: function(){
		this.selection = new dojox.grid.enhanced.DataSelection(this);
	},
	canSort: function(colIndex, field){
		// summary:
		//		Overwritten
		return true;
	},
	doKeyEvent: function(e){
		// summary:
		//		Overwritten, see _Grid.doKeyEvent()
		try{
			var view = this.focus.focusView;
			view.content.decorateEvent(e);
			if(!e.cell){ view.header.decorateEvent(e); }
		}catch(e){}
		this.inherited(arguments);
	},
	doApplyCellEdit: function(inValue, inRowIndex, inAttrName){
		// summary:
		//		Overwritten, see DataGrid.doApplyCellEdit()
		if(!inAttrName){
			this.invalidated[inRowIndex] = true;
			return;
		}
		this.inherited(arguments);
	},
	mixin: function(target, source){
		var props = {};
		for(var p in source){
			if(p == '_inherited' || p == 'declaredClass' || p == 'constructor' ||
				source['privates'] && source['privates'][p]){
				continue;
			}
			props[p] = source[p];
		}
		lang.mixin(target, props);
	},
	_copyAttr: function(idx, attr){
		// summary:
		//		Overwritten, see DataGrid._copyAttr()
		//		Fix cell TAB navigation for single click editing
		if(!attr){ return; }
		return this.inherited(arguments);
	},
	_getHeaderHeight: function(){
		// summary:
		//		Overwritten, see _Grid._getHeaderHeight()
		//		Should include borders/margins of this.viewsHeaderNode
		this.inherited(arguments);
		return domGeometry.getMarginBox(this.viewsHeaderNode).h;
	},
	_fetch: function(start, isRender){
		// summary:
		//		Overwritten, see DataGrid._fetch()
		if(this.items){
			return this.inherited(arguments);
		}
		start = start || 0;
		if(this.store && !this._pending_requests[start]){
			if(!this._isLoaded && !this._isLoading){
				this._isLoading = true;
				this.showMessage(this.loadingMessage);
			}
			this._pending_requests[start] = true;
			try{
				var req = {
					start: start,
					count: this.rowsPerPage,
					query: this.query,
					sort: this.getSortProps(),
					queryOptions: this.queryOptions,
					isRender: isRender,
					onBegin: lang.hitch(this, "_onFetchBegin"),
					onComplete: lang.hitch(this, "_onFetchComplete"),
					onError: lang.hitch(this, "_onFetchError")
				};
				this._storeLayerFetch(req);
			}catch(e){
				this._onFetchError(e, {start: start, count: this.rowsPerPage});
			}
		}
		return 0;
	},
	_storeLayerFetch: function(req){
		// summary:
		//		Extracted fetch specifically for store layer use
		this.store.fetch(req);
	},
	getCellByField: function(field){
		return array.filter(this.layout.cells, function(cell){
			return cell.field == field;
		})[0];
	},
	onMouseUp: function(e){	},
	createView: function(){
		// summary:
		//		Overwrite: rewrite getCellX of view.header
		var view = this.inherited(arguments);
		if(has('mozilla')){
			var ascendDom = function(inNode, inWhile){
				for(var n = inNode; n && inWhile(n); n = n.parentNode){}
				return n;
			};//copied from dojox.grid._Builder
			var makeNotTagName = function(inTagName){
				var name = inTagName.toUpperCase();
				return function(node){ return node.tagName != name; };
			};//copied from dojox.grid._Builder

			var func = view.header.getCellX;
			view.header.getCellX = function(e){
				var x = func.call(view.header, e);
				var n = ascendDom(e.target, makeNotTagName("th"));
				if(n && n !== e.target && dom.isDescendant(e.target, n)){ x += n.firstChild.offsetLeft; }
				return x;
			};
		}
		return view;
	},
	destroy: function(){
		// summary:
		//		Destroy all resources
		delete this._nls;
		this.pluginMgr.destroy();
		this.inherited(arguments);
	}
});

declare("dojox.grid.enhanced.DataSelection", DataSelection, {
	constructor: function(grid){
		if(grid.keepSelection){
			if(this.preserver){
				this.preserver.destroy();
			}
			this.preserver = new _SelectionPreserver(this);
		}
	},
	_range: function(inFrom, inTo){
		this.grid._selectingRange = true;
		this.inherited(arguments);
		this.grid._selectingRange = false;
		this.onChanged();
	},
	deselectAll: function(inItemOrIndex){
		this.grid._selectingRange = true;
		this.inherited(arguments);
		this.grid._selectingRange = false;
		this.onChanged();
	}
});

EnhancedGrid.markupFactory = function(props, node, ctor, cellFunc){
	return dojox.grid._Grid.markupFactory(props, node, ctor,
					lang.partial(DataGrid.cell_markupFactory, cellFunc));
};

EnhancedGrid.registerPlugin = function(clazz, props){
	_PluginManager.registerPlugin(clazz, props);
};

return EnhancedGrid;

});

},
'dojox/grid/DataGrid':function(){
define("dojox/grid/DataGrid", [
	"../main",
	"dojo/_base/array",
	"dojo/_base/lang",
	"dojo/_base/json",
	"dojo/_base/sniff",
	"dojo/_base/declare",
	"./_Grid",
	"./DataSelection",
	"dojo/_base/html"
], function(dojox, array, lang, json, has, declare, _Grid, DataSelection, html){

/*=====
declare("dojox.grid.__DataCellDef", dojox.grid.__CellDef, {
	constructor: function(){
		// field: String?
		//		The attribute to read from the dojo.data item for the row.
		// fields: String[]?
		//		An array of fields to grab the values of and pass as an array to the grid
		// get: Function?
		//		function(rowIndex, item?){} rowIndex is of type Integer, item is of type
		//		Object.  This function will be called when a cell requests data.  Returns
		//		the unformatted data for the cell.
	}
});
=====*/

/*=====
declare("dojox.grid.__DataViewDef", dojox.grid.__ViewDef, {
	constructor: function(){
		// cells: dojox.grid.__DataCellDef[]|Array[dojox.grid.__DataCellDef[]]?
		//		The structure of the cells within this grid.
		// defaultCell: dojox.grid.__DataCellDef?
		//		A cell definition with default values for all cells in this view.  If
		//		a property is defined in a cell definition in the "cells" array and
		//		this property, the cell definition's property will override this
		//		property's property.
	}
});
=====*/

var DataGrid = declare("dojox.grid.DataGrid", _Grid, {
	store: null,
	query: null,
	queryOptions: null,
	fetchText: '...',
	sortFields: null,
	
	// updateDelay: int
	//		Time, in milliseconds, to delay updates automatically so that multiple
	//		calls to onSet/onNew/onDelete don't keep rerendering the grid.  Set
	//		to 0 to immediately cause updates.  A higher value will result in
	//		better performance at the expense of responsiveness of the grid.
	updateDelay: 1,

/*=====
	// structure: dojox.grid.__DataViewDef|dojox.grid.__DataViewDef[]|dojox.grid.__DataCellDef[]|Array[dojox.grid.__DataCellDef[]]
	//		View layout definition.
	structure: '',
=====*/

	// You can specify items instead of a query, if you like.  They do not need
	// to be loaded - but the must be items in the store
	items: null,
	
	_store_connects: null,
	_by_idty: null,
	_by_idx: null,
	_cache: null,
	_pages: null,
	_pending_requests: null,
	_bop: -1,
	_eop: -1,
	_requests: 0,
	rowCount: 0,

	_isLoaded: false,
	_isLoading: false,
	
	// keepSelection: Boolean
	//		Whether keep selection after sort, filter etc.
	keepSelection: false,	
	
	postCreate: function(){
		this._pages = [];
		this._store_connects = [];
		this._by_idty = {};
		this._by_idx = [];
		this._cache = [];
		this._pending_requests = {};

		this._setStore(this.store);
		this.inherited(arguments);
	},
	
	destroy: function(){
		this.selection.destroy();
		this.inherited(arguments);
	},

	createSelection: function(){
		this.selection = new DataSelection(this);
	},

	get: function(inRowIndex, inItem){
		// summary:
		//		Default data getter.
		// description:
		//		Provides data to display in a grid cell. Called in grid cell context.
		//		So this.cell.index is the column index.
		// inRowIndex: Integer
		//		Row for which to provide data
		// returns:
		//		Data to display for a given grid cell.
		
		if(inItem && this.field == "_item" && !this.fields){
			return inItem;
		}else if(inItem && this.fields){
			var ret = [];
			var s = this.grid.store;
			array.forEach(this.fields, function(f){
				ret = ret.concat(s.getValues(inItem, f));
			});
			return ret;
		}else if(!inItem && typeof inRowIndex === "string"){
			return this.inherited(arguments);
		}
		return (!inItem ? this.defaultValue : (!this.field ? this.value : (this.field == "_item" ? inItem : this.grid.store.getValue(inItem, this.field))));
	},

	_checkUpdateStatus: function(){
		if(this.updateDelay > 0){
			var iStarted = false;
			if(this._endUpdateDelay){
				clearTimeout(this._endUpdateDelay);
				delete this._endUpdateDelay;
				iStarted = true;
			}
			if(!this.updating){
				this.beginUpdate();
				iStarted = true;
			}
			if(iStarted){
				var _this = this;
				this._endUpdateDelay = setTimeout(function(){
					delete _this._endUpdateDelay;
					_this.endUpdate();
				}, this.updateDelay);
			}
		}
	},
	
	_onSet: function(item, attribute, oldValue, newValue){
		this._checkUpdateStatus();
		var idx = this.getItemIndex(item);
		if(idx>-1){
			this.updateRow(idx);
		}
	},
	
	_createItem: function(item, index){
		var idty = this._hasIdentity ? this.store.getIdentity(item) : json.toJson(this.query) + ":idx:" + index + ":sort:" + json.toJson(this.getSortProps());
		var o = this._by_idty[idty] = { idty: idty, item: item };
		return o;
	},

	_addItem: function(item, index, noUpdate){
		this._by_idx[index] = this._createItem(item, index);
		if(!noUpdate){
			this.updateRow(index);
		}
	},

	_onNew: function(item, parentInfo){
		this._checkUpdateStatus();
		var rowCount = this.get('rowCount');
		this._addingItem = true;
		this.updateRowCount(rowCount+1);
		this._addingItem = false;
		this._addItem(item, rowCount);
		this.showMessage();
	},

	_onDelete: function(item){
		this._checkUpdateStatus();
		var idx = this._getItemIndex(item, true);

		if(idx >= 0){
			// When a row is deleted, all rest rows are shifted down,
			// and migrate from page to page. If some page is not
			// loaded yet empty rows can migrate to initialized pages
			// without refreshing. It causes empty rows in some pages, see:
			// http://bugs.dojotoolkit.org/ticket/6818
			// this code fix this problem by reseting loaded page info
			this._pages = [];
			this._bop = -1;
			this._eop = -1;

			var o = this._by_idx[idx];
			this._by_idx.splice(idx, 1);
			delete this._by_idty[o.idty];
			this.updateRowCount(this.get('rowCount')-1);
			if(this.get('rowCount') === 0){
				this.showMessage(this.noDataMessage);
			}
		}
		if(this.selection.isSelected(idx)){
			this.selection.deselect(idx);
			this.selection.selected.splice(idx, 1);
		}
	},

	_onRevert: function(){
		this._refresh();
	},

	setStore: function(store, query, queryOptions){
		if(this._requestsPending(0)){
			return;
		}
		this._setQuery(query, queryOptions);
		this._setStore(store);
		this._refresh(true);
	},
	
	setQuery: function(query, queryOptions){
		if(this._requestsPending(0)){
			return;
		}
		this._setQuery(query, queryOptions);
		this._refresh(true);
	},
	
	setItems: function(items){
		this.items = items;
		this._setStore(this.store);
		this._refresh(true);
	},
	
	_setQuery: function(query, queryOptions){
		this.query = query;
		this.queryOptions = queryOptions || this.queryOptions;
	},

	_setStore: function(store){
		if(this.store && this._store_connects){
			array.forEach(this._store_connects, this.disconnect, this);
		}
		this.store = store;

		if(this.store){
			var f = this.store.getFeatures();
			var h = [];

			this._canEdit = !!f["dojo.data.api.Write"] && !!f["dojo.data.api.Identity"];
			this._hasIdentity = !!f["dojo.data.api.Identity"];

			if(!!f["dojo.data.api.Notification"] && !this.items){
				h.push(this.connect(this.store, "onSet", "_onSet"));
				h.push(this.connect(this.store, "onNew", "_onNew"));
				h.push(this.connect(this.store, "onDelete", "_onDelete"));
			}
			if(this._canEdit){
				h.push(this.connect(this.store, "revert", "_onRevert"));
			}

			this._store_connects = h;
		}
	},

	_onFetchBegin: function(size, req){
		if(!this.scroller){ return; }
		if(this.rowCount != size){
			if(req.isRender){
				this.scroller.init(size, this.keepRows, this.rowsPerPage);
				this.rowCount = size;
				this._setAutoHeightAttr(this.autoHeight, true);
				this._skipRowRenormalize = true;
				this.prerender();
				this._skipRowRenormalize = false;
			}else{
				this.updateRowCount(size);
			}
		}
		if(!size){
			this.views.render();
			this._resize();
			this.showMessage(this.noDataMessage);
			this.focus.initFocusView();
		}else{
			this.showMessage();
		}
	},

	_onFetchComplete: function(items, req){
		if(!this.scroller){ return; }
		if(items && items.length > 0){
			//console.log(items);
			array.forEach(items, function(item, idx){
				this._addItem(item, req.start+idx, true);
			}, this);
			this.updateRows(req.start, items.length);
			if(req.isRender){
				this.setScrollTop(0);
				this.postrender();
			}else if(this._lastScrollTop){
				this.setScrollTop(this._lastScrollTop);
			}
			if(has("ie")){
				html.setSelectable(this.domNode, this.selectable);
			}	
		}
		delete this._lastScrollTop;
		if(!this._isLoaded){
			this._isLoading = false;
			this._isLoaded = true;
		}
		this._pending_requests[req.start] = false;
	},

	_onFetchError: function(err, req){
		console.log(err);
		delete this._lastScrollTop;
		if(!this._isLoaded){
			this._isLoading = false;
			this._isLoaded = true;
			this.showMessage(this.errorMessage);
		}
		this._pending_requests[req.start] = false;
		this.onFetchError(err, req);
	},

	onFetchError: function(err, req){
	},

	_fetch: function(start, isRender){
		start = start || 0;
		if(this.store && !this._pending_requests[start]){
			if(!this._isLoaded && !this._isLoading){
				this._isLoading = true;
				this.showMessage(this.loadingMessage);
			}
			this._pending_requests[start] = true;
			//console.log("fetch: ", start);
			try{
				if(this.items){
					var items = this.items;
					var store = this.store;
					this.rowsPerPage = items.length;
					var req = {
						start: start,
						count: this.rowsPerPage,
						isRender: isRender
					};
					this._onFetchBegin(items.length, req);
					
					// Load them if we need to
					var waitCount = 0;
					array.forEach(items, function(i){
						if(!store.isItemLoaded(i)){ waitCount++; }
					});
					if(waitCount === 0){
						this._onFetchComplete(items, req);
					}else{
						var onItem = function(item){
							waitCount--;
							if(waitCount === 0){
								this._onFetchComplete(items, req);
							}
						};
						array.forEach(items, function(i){
							if(!store.isItemLoaded(i)){
								store.loadItem({item: i, onItem: onItem, scope: this});
							}
						}, this);
					}
				}else{
					this.store.fetch({
						start: start,
						count: this.rowsPerPage,
						query: this.query,
						sort: this.getSortProps(),
						queryOptions: this.queryOptions,
						isRender: isRender,
						onBegin: lang.hitch(this, "_onFetchBegin"),
						onComplete: lang.hitch(this, "_onFetchComplete"),
						onError: lang.hitch(this, "_onFetchError")
					});
				}
			}catch(e){
				this._onFetchError(e, {start: start, count: this.rowsPerPage});
			}
		}
	},

	_clearData: function(){
		this.updateRowCount(0);
		this._by_idty = {};
		this._by_idx = [];
		this._pages = [];
		this._bop = this._eop = -1;
		this._isLoaded = false;
		this._isLoading = false;
	},

	getItem: function(idx){
		var data = this._by_idx[idx];
		if(!data||(data&&!data.item)){
			this._preparePage(idx);
			return null;
		}
		return data.item;
	},

	getItemIndex: function(item){
		return this._getItemIndex(item, false);
	},
	
	_getItemIndex: function(item, isDeleted){
		if(!isDeleted && !this.store.isItem(item)){
			return -1;
		}

		var idty = this._hasIdentity ? this.store.getIdentity(item) : null;

		for(var i=0, l=this._by_idx.length; i<l; i++){
			var d = this._by_idx[i];
			if(d && ((idty && d.idty == idty) || (d.item === item))){
				return i;
			}
		}
		return -1;
	},

	filter: function(query, reRender){
		this.query = query;
		if(reRender){
			this._clearData();
		}
		this._fetch();
	},

	_getItemAttr: function(idx, attr){
		var item = this.getItem(idx);
		return (!item ? this.fetchText : this.store.getValue(item, attr));
	},

	// rendering
	_render: function(){
		if(this.domNode.parentNode){
			this.scroller.init(this.get('rowCount'), this.keepRows, this.rowsPerPage);
			this.prerender();
			this._fetch(0, true);
		}
	},

	// paging
	_requestsPending: function(inRowIndex){
		return this._pending_requests[inRowIndex];
	},

	_rowToPage: function(inRowIndex){
		return (this.rowsPerPage ? Math.floor(inRowIndex / this.rowsPerPage) : inRowIndex);
	},

	_pageToRow: function(inPageIndex){
		return (this.rowsPerPage ? this.rowsPerPage * inPageIndex : inPageIndex);
	},

	_preparePage: function(inRowIndex){
		if((inRowIndex < this._bop || inRowIndex >= this._eop) && !this._addingItem){
			var pageIndex = this._rowToPage(inRowIndex);
			this._needPage(pageIndex);
			this._bop = pageIndex * this.rowsPerPage;
			this._eop = this._bop + (this.rowsPerPage || this.get('rowCount'));
		}
	},

	_needPage: function(inPageIndex){
		if(!this._pages[inPageIndex]){
			this._pages[inPageIndex] = true;
			this._requestPage(inPageIndex);
		}
	},

	_requestPage: function(inPageIndex){
		var row = this._pageToRow(inPageIndex);
		var count = Math.min(this.rowsPerPage, this.get('rowCount') - row);
		if(count > 0){
			this._requests++;
			if(!this._requestsPending(row)){
				setTimeout(lang.hitch(this, "_fetch", row, false), 1);
				//this.requestRows(row, count);
			}
		}
	},

	getCellName: function(inCell){
		return inCell.field;
		//console.log(inCell);
	},

	_refresh: function(isRender){
		this._clearData();
		this._fetch(0, isRender);
	},

	sort: function(){
		this.edit.apply();
		this._lastScrollTop = this.scrollTop;
		this._refresh();
	},

	canSort: function(){
		return (!this._isLoading);
	},

	getSortProps: function(){
		var c = this.getCell(this.getSortIndex());
		if(!c){
			if(this.sortFields){
				return this.sortFields;
			}
			return null;
		}else{
			var desc = c["sortDesc"];
			var si = !(this.sortInfo>0);
			if(typeof desc == "undefined"){
				desc = si;
			}else{
				desc = si ? !desc : desc;
			}
			return [{ attribute: c.field, descending: desc }];
		}
	},

	styleRowState: function(inRow){
		// summary:
		//		Perform row styling
		if(this.store && this.store.getState){
			var states=this.store.getState(inRow.index), c='';
			for(var i=0, ss=["inflight", "error", "inserting"], s; s=ss[i]; i++){
				if(states[s]){
					c = ' dojoxGridRow-' + s;
					break;
				}
			}
			inRow.customClasses += c;
		}
	},

	onStyleRow: function(inRow){
		this.styleRowState(inRow);
		this.inherited(arguments);
	},

	// editing
	canEdit: function(inCell, inRowIndex){
		return this._canEdit;
	},

	_copyAttr: function(idx, attr){
		var row = {};
		var backstop = {};
		var src = this.getItem(idx);
		return this.store.getValue(src, attr);
	},

	doStartEdit: function(inCell, inRowIndex){
		if(!this._cache[inRowIndex]){
			this._cache[inRowIndex] = this._copyAttr(inRowIndex, inCell.field);
		}
		this.onStartEdit(inCell, inRowIndex);
	},

	doApplyCellEdit: function(inValue, inRowIndex, inAttrName){
		this.store.fetchItemByIdentity({
			identity: this._by_idx[inRowIndex].idty,
			onItem: lang.hitch(this, function(item){
				var oldValue = this.store.getValue(item, inAttrName);
				if(typeof oldValue == 'number'){
					inValue = isNaN(inValue) ? inValue : parseFloat(inValue);
				}else if(typeof oldValue == 'boolean'){
					inValue = inValue == 'true' ? true : inValue == 'false' ? false : inValue;
				}else if(oldValue instanceof Date){
					var asDate = new Date(inValue);
					inValue = isNaN(asDate.getTime()) ? inValue : asDate;
				}
				this.store.setValue(item, inAttrName, inValue);
				this.onApplyCellEdit(inValue, inRowIndex, inAttrName);
			})
		});
	},

	doCancelEdit: function(inRowIndex){
		var cache = this._cache[inRowIndex];
		if(cache){
			this.updateRow(inRowIndex);
			delete this._cache[inRowIndex];
		}
		this.onCancelEdit.apply(this, arguments);
	},

	doApplyEdit: function(inRowIndex, inDataAttr){
		var cache = this._cache[inRowIndex];
		/*if(cache){
			var data = this.getItem(inRowIndex);
			if(this.store.getValue(data, inDataAttr) != cache){
				this.update(cache, data, inRowIndex);
			}
			delete this._cache[inRowIndex];
		}*/
		this.onApplyEdit(inRowIndex);
	},

	removeSelectedRows: function(){
		// summary:
		//		Remove the selected rows from the grid.
		if(this._canEdit){
			this.edit.apply();
			var fx = lang.hitch(this, function(items){
				if(items.length){
					array.forEach(items, this.store.deleteItem, this.store);
					this.selection.clear();
				}
			});
			if(this.allItemsSelected){
				this.store.fetch({
							query: this.query,
							queryOptions: this.queryOptions,
							onComplete: fx});
			}else{
				fx(this.selection.getSelected());
			}
		}
	}
});

DataGrid.cell_markupFactory = function(cellFunc, node, cellDef){
	var field = lang.trim(html.attr(node, "field")||"");
	if(field){
		cellDef.field = field;
	}
	cellDef.field = cellDef.field||cellDef.name;
	var fields = lang.trim(html.attr(node, "fields")||"");
	if(fields){
		cellDef.fields = fields.split(",");
	}
	if(cellFunc){
		cellFunc(node, cellDef);
	}
};

DataGrid.markupFactory = function(props, node, ctor, cellFunc){
	return _Grid.markupFactory(props, node, ctor,
					lang.partial(DataGrid.cell_markupFactory, cellFunc));
};

return DataGrid;

});

},
'dojox/grid/_Grid':function(){
require({cache:{
'url:dojox/grid/resources/_Grid.html':"<div hidefocus=\"hidefocus\" role=\"grid\" dojoAttachEvent=\"onmouseout:_mouseOut\">\n\t<div class=\"dojoxGridMasterHeader\" dojoAttachPoint=\"viewsHeaderNode\" role=\"presentation\"></div>\n\t<div class=\"dojoxGridMasterView\" dojoAttachPoint=\"viewsNode\" role=\"presentation\"></div>\n\t<div class=\"dojoxGridMasterMessages\" style=\"display: none;\" dojoAttachPoint=\"messagesNode\"></div>\n\t<span dojoAttachPoint=\"lastFocusNode\" tabindex=\"0\"></span>\n</div>\n"}});
define("dojox/grid/_Grid", [
	"dojo/_base/kernel",
	"../main",
	"dojo/_base/declare",
	"./_Events",
	"./_Scroller",
	"./_Layout",
	"./_View",
	"./_ViewManager",
	"./_RowManager",
	"./_FocusManager",
	"./_EditManager",
	"./Selection",
	"./_RowSelector",
	"./util",
	"dijit/_Widget",
	"dijit/_TemplatedMixin",
	"dijit/CheckedMenuItem",
	"dojo/text!./resources/_Grid.html",
	"dojo/string",
	"dojo/_base/array",
	"dojo/_base/lang",
	"dojo/_base/sniff",
	"dojox/html/metrics",
	"dojo/_base/html",
	"dojo/query",
	"dojo/dnd/common",
	"dojo/i18n!dijit/nls/loading"
], function(dojo, dojox, declare, _Events, _Scroller, _Layout, _View, _ViewManager,
	_RowManager, _FocusManager, _EditManager, Selection, _RowSelector, util, _Widget,
	 _TemplatedMixin, CheckedMenuItem, template, string, array, lang, has, metrics, html, query){

	// NOTE: this is for backwards compatibility with Dojo 1.3
	if(!dojo.isCopyKey){
		dojo.isCopyKey = dojo.dnd.getCopyKeyState;
	}
	/*=====
	dojox.grid.__CellDef = {
		// name: String?
		//		The text to use in the header of the grid for this cell.
		// get: Function?
		//		function(rowIndex){} rowIndex is of type Integer.  This
		//		function will be called when a cell	requests data.  Returns the
		//		unformatted data for the cell.
		// value: String?
		//		If "get" is not specified, this is used as the data for the cell.
		// defaultValue: String?
		//		If "get" and "value" aren't specified or if "get" returns an undefined
		//		value, this is used as the data for the cell.  "formatter" is not run
		//		on this if "get" returns an undefined value.
		// formatter: Function?
		//		function(data, rowIndex){} data is of type anything, rowIndex
		//		is of type Integer.  This function will be called after the cell
		//		has its data but before it passes it back to the grid to render.
		//		Returns the formatted version of the cell's data.
		// type: dojox.grid.cells._Base|Function?
		//		TODO
		// editable: Boolean?
		//		Whether this cell should be editable or not.
		// hidden: Boolean?
		//		If true, the cell will not be displayed.
		// noresize: Boolean?
		//		If true, the cell will not be able to be resized.
		// width: Integer|String?
		//		A CSS size.  If it's an Integer, the width will be in em's.
		// colSpan: Integer?
		//		How many columns to span this cell.  Will not work in the first
		//		sub-row of cells.
		// rowSpan: Integer?
		//		How many sub-rows to span this cell.
		// styles: String?
		//		A string of styles to apply to both the header cell and main
		//		grid cells.  Must end in a ';'.
		// headerStyles: String?
		//		A string of styles to apply to just the header cell.  Must end
		//		in a ';'
		// cellStyles: String?
		//		A string of styles to apply to just the main grid cells.  Must
		//		end in a ';'
		// classes: String?
		//		A space separated list of classes to apply to both the header
		//		cell and the main grid cells.
		// headerClasses: String?
		//		A space separated list of classes to apply to just the header
		//		cell.
		// cellClasses: String?
		//		A space separated list of classes to apply to just the main
		//		grid cells.
		// attrs: String?
		//		A space separated string of attribute='value' pairs to add to
		//		the header cell element and main grid cell elements.
	};
	=====*/

	/*=====
	dojox.grid.__ViewDef = {
		// noscroll: Boolean?
		//		If true, no scrollbars will be rendered without scrollbars.
		// width: Integer|String?
		//		A CSS size.  If it's an Integer, the width will be in em's. If
		//		"noscroll" is true, this value is ignored.
		// cells: dojox.grid.__CellDef[]|Array[dojox.grid.__CellDef[]]?
		//		The structure of the cells within this grid.
		// type: String?
		//		A string containing the constructor of a subclass of
		//		dojox.grid._View.  If this is not specified, dojox.grid._View
		//		is used.
		// defaultCell: dojox.grid.__CellDef?
		//		A cell definition with default values for all cells in this view.  If
		//		a property is defined in a cell definition in the "cells" array and
		//		this property, the cell definition's property will override this
		//		property's property.
		// onBeforeRow: Function?
		//		function(rowIndex, cells){} rowIndex is of type Integer, cells
		//		is of type Array[dojox.grid.__CellDef[]].  This function is called
		//		before each row of data is rendered.  Before the header is
		//		rendered, rowIndex will be -1.  "cells" is a reference to the
		//		internal structure of this view's cells so any changes you make to
		//		it will persist between calls.
		// onAfterRow: Function?
		//		function(rowIndex, cells, rowNode){} rowIndex is of type Integer, cells
		//		is of type Array[dojox.grid.__CellDef[]], rowNode is of type DOMNode.
		//		This function is called	after each row of data is rendered.  After the
		//		header is rendered, rowIndex will be -1.  "cells" is a reference to the
		//		internal structure of this view's cells so any changes you make to
		//		it will persist between calls.
	};
	=====*/

	var _Grid = declare('dojox.grid._Grid',
		[ _Widget, _TemplatedMixin, _Events ],
		{
		// summary:
		//		A grid widget with virtual scrolling, cell editing, complex rows,
		//		sorting, fixed columns, sizeable columns, etc.
		//
		// description:
		//		_Grid provides the full set of grid features without any
		//		direct connection to a data store.
		//
		//		The grid exposes a get function for the grid, or optionally
		//		individual columns, to populate cell contents.
		//
		//		The grid is rendered based on its structure, an object describing
		//		column and cell layout.
		//
		// example:
		//		A quick sample:
		//
		//		define a get function
		//	|	function get(inRowIndex){ // called in cell context
		//	|		return [this.index, inRowIndex].join(', ');
		//	|	}
		//
		//		define the grid structure:
		//	|	var structure = [ // array of view objects
		//	|		{ cells: [// array of rows, a row is an array of cells
		//	|			[
		//	|				{ name: "Alpha", width: 6 },
		//	|				{ name: "Beta" },
		//	|				{ name: "Gamma", get: get }]
		//	|		]}
		//	|	];
		//
		//	|	<div id="grid"
		//	|		rowCount="100" get="get"
		//	|		structure="structure"
		//	|		dojoType="dojox.grid._Grid"></div>

		templateString: template,

		// classTag: String
		//		CSS class applied to the grid's domNode
		classTag: 'dojoxGrid',

		// settings
		// rowCount: Integer
		//		Number of rows to display.
		rowCount: 5,

		// keepRows: Integer
		//		Number of rows to keep in the rendering cache.
		keepRows: 75,

		// rowsPerPage: Integer
		//		Number of rows to render at a time.
		rowsPerPage: 25,

		// autoWidth: Boolean
		//		If autoWidth is true, grid width is automatically set to fit the data.
		autoWidth: false,
		
		// initialWidth: String
		//		A css string to use to set our initial width (only used if autoWidth
		//		is true).  The first rendering of the grid will be this width, any
		//		resizing of columns, etc will result in the grid switching to
		//		autoWidth mode.  Note, this width will override any styling in a
		//		stylesheet or directly on the node.
		initialWidth: "",

		// autoHeight: Boolean|Integer
		//		If autoHeight is true, grid height is automatically set to fit the data.
		//		If it is an integer, the height will be automatically set to fit the data
		//		if there are fewer than that many rows - and the height will be set to show
		//		that many rows if there are more
		autoHeight: '',

		// rowHeight: Integer
		//		If rowHeight is set to a positive number, it will define the height of the rows
		//		in pixels. This can provide a significant performance advantage, since it
		//		eliminates the need to measure row sizes during rendering, which is one
		//		the primary bottlenecks in the DataGrid's performance.
		rowHeight: 0,
		
		// autoRender: Boolean
		//		If autoRender is true, grid will render itself after initialization.
		autoRender: true,

		// defaultHeight: String
		//		default height of the grid, measured in any valid css unit.
		defaultHeight: '15em',
		
		// height: String
		//		explicit height of the grid, measured in any valid css unit.  This will be populated (and overridden)
		//		if the height: css attribute exists on the source node.
		height: '',

		// structure: dojox.grid.__ViewDef|dojox.grid.__ViewDef[]|dojox.grid.__CellDef[]|Array[dojox.grid.__CellDef[]]
		//		View layout defintion.
		structure: null,

		// elasticView: Integer
		//	Override defaults and make the indexed grid view elastic, thus filling available horizontal space.
		elasticView: -1,

		// singleClickEdit: boolean
		//		Single-click starts editing. Default is double-click
		singleClickEdit: false,

		// selectionMode: String
		//		Set the selection mode of grid's Selection.  Value must be 'single', 'multiple',
		//		or 'extended'.  Default is 'extended'.
		selectionMode: 'extended',

		// rowSelector: Boolean|String
		//		If set to true, will add a row selector view to this grid.  If set to a CSS width, will add
		//		a row selector of that width to this grid.
		rowSelector: '',

		// columnReordering: Boolean
		//		If set to true, will add drag and drop reordering to views with one row of columns.
		columnReordering: false,

		// headerMenu: dijit.Menu
		//		If set to a dijit.Menu, will use this as a context menu for the grid headers.
		headerMenu: null,

		// placeholderLabel: String
		//		Label of placeholders to search for in the header menu to replace with column toggling
		//		menu items.
		placeholderLabel: "GridColumns",
		
		// selectable: Boolean
		//		Set to true if you want to be able to select the text within the grid.
		selectable: false,
		
		// Used to store the last two clicks, to ensure double-clicking occurs based on the intended row
		_click: null,
		
		// loadingMessage: String
		//		Message that shows while the grid is loading
		loadingMessage: "<span class='dojoxGridLoading'>${loadingState}</span>",

		// errorMessage: String
		//		Message that shows when the grid encounters an error loading
		errorMessage: "<span class='dojoxGridError'>${errorState}</span>",

		// noDataMessage: String
		//		Message that shows if the grid has no data - wrap it in a
		//		span with class 'dojoxGridNoData' if you want it to be
		//		styled similar to the loading and error messages
		noDataMessage: "",
		
		// escapeHTMLInData: Boolean
		//		This will escape HTML brackets from the data to prevent HTML from
		//		user-inputted data being rendered with may contain JavaScript and result in
		//		XSS attacks. This is true by default, and it is recommended that it remain
		//		true. Setting this to false will allow data to be displayed in the grid without
		//		filtering, and should be only used if it is known that the data won't contain
		//		malicious scripts. If HTML is needed in grid cells, it is recommended that
		//		you use the formatter function to generate the HTML (the output of
		//		formatter functions is not filtered, even with escapeHTMLInData set to true).
		escapeHTMLInData: true,
		
		// formatterScope: Object
		//		An object to execute format functions within.  If not set, the
		//		format functions will execute within the scope of the cell that
		//		has a format function.
		formatterScope: null,
		
		// editable: boolean
		//		indicates if the grid contains editable cells, default is false
		//		set to true if editable cell encountered during rendering
		editable: false,

		// summary: String
		//		Customizable summary descriptions which will be added to grid.domNode
		summary: '',
		_setSummaryAttr: 'domNode',
		
		// sortInfo: [private] Number
		sortInfo: 0,

		// _placeholders: [private] Array
		_placeholders: null,

		// _layoutClass: Object
		//	The class to use for our layout - can be overridden by grid subclasses
		_layoutClass: _Layout,

		// initialization
		buildRendering: function(){
			this.inherited(arguments);
			if(!this.domNode.getAttribute('tabIndex')){
				this.domNode.tabIndex = "0";
			}
			this.createScroller();
			this.createLayout();
			this.createViews();
			this.createManagers();

			this.createSelection();

			this.connect(this.selection, "onSelected", "onSelected");
			this.connect(this.selection, "onDeselected", "onDeselected");
			this.connect(this.selection, "onChanged", "onSelectionChanged");

			metrics.initOnFontResize();
			this.connect(metrics, "onFontResize", "textSizeChanged");
			util.funnelEvents(this.domNode, this, 'doKeyEvent', util.keyEvents);
			if (this.selectionMode != "none") {
				this.domNode.setAttribute("aria-multiselectable", this.selectionMode == "single" ? "false" : "true");
			}

			html.addClass(this.domNode, this.classTag);
			if(!this.isLeftToRight()){
				html.addClass(this.domNode, this.classTag+"Rtl");
			}
			if(this.rowHeight > 0){
				html.addClass(this.viewsNode, this.classTag + "FixedRowHeight");
			}
		},
		
		postMixInProperties: function(){
			this.inherited(arguments);
			var messages = dojo.i18n.getLocalization("dijit", "loading", this.lang);
			this.loadingMessage = string.substitute(this.loadingMessage, messages);
			this.errorMessage = string.substitute(this.errorMessage, messages);
			if(this.srcNodeRef && this.srcNodeRef.style.height){
				this.height = this.srcNodeRef.style.height;
			}
			// Call this to update our autoheight to start out
			this._setAutoHeightAttr(this.autoHeight, true);
			this.lastScrollTop = this.scrollTop = 0;
		},
		
		postCreate: function(){
			this._placeholders = [];
			this._setHeaderMenuAttr(this.headerMenu);
			this._setStructureAttr(this.structure);
			this._click = [];
			this.inherited(arguments);
			if(this.domNode && this.autoWidth && this.initialWidth){
				this.domNode.style.width = this.initialWidth;
			}
			if (this.domNode && !this.editable){
				// default value for aria-readonly is false, set to true if grid is not editable
				html.attr(this.domNode,"aria-readonly", "true");
			}
		},

		destroy: function(){
			this.domNode.onReveal = null;
			this.domNode.onSizeChange = null;

			// Fixes IE domNode leak
			delete this._click;

			if(this.scroller){
				this.scroller.destroy();
				delete this.scroller;
			}
			this.edit.destroy();
			delete this.edit;
			this.views.destroyViews();
			if(this.focus){
				this.focus.destroy();
				delete this.focus;
			}
			if(this.headerMenu&&this._placeholders.length){
				array.forEach(this._placeholders, function(p){ p.unReplace(true); });
				this.headerMenu.unBindDomNode(this.viewsHeaderNode);
			}
			this.inherited(arguments);
		},

		_setAutoHeightAttr: function(ah, skipRender){
			// Calculate our autoheight - turn it into a boolean or an integer
			if(typeof ah == "string"){
				if(!ah || ah == "false"){
					ah = false;
				}else if (ah == "true"){
					ah = true;
				}else{
					ah = window.parseInt(ah, 10);
				}
			}
			if(typeof ah == "number"){
				if(isNaN(ah)){
					ah = false;
				}
				// Autoheight must be at least 1, if it's a number.  If it's
				// less than 0, we'll take that to mean "all" rows (same as
				// autoHeight=true - if it is equal to zero, we'll take that
				// to mean autoHeight=false
				if(ah < 0){
					ah = true;
				}else if (ah === 0){
					ah = false;
				}
			}
			this.autoHeight = ah;
			if(typeof ah == "boolean"){
				this._autoHeight = ah;
			}else if(typeof ah == "number"){
				this._autoHeight = (ah >= this.get('rowCount'));
			}else{
				this._autoHeight = false;
			}
			if(this._started && !skipRender){
				this.render();
			}
		},

		_getRowCountAttr: function(){
			return this.updating && this.invalidated && this.invalidated.rowCount != undefined ?
				this.invalidated.rowCount : this.rowCount;
		},
		
		textSizeChanged: function(){
			this.render();
		},

		sizeChange: function(){
			this.update();
		},

		createManagers: function(){
			// summary:
			//		create grid managers for various tasks including rows, focus, selection, editing

			// row manager
			this.rows = new _RowManager(this);
			// focus manager
			this.focus = new _FocusManager(this);
			// edit manager
			this.edit = new _EditManager(this);
		},

		createSelection: function(){
			// summary:	Creates a new Grid selection manager.

			// selection manager
			this.selection = new Selection(this);
		},

		createScroller: function(){
			// summary:
			//		Creates a new virtual scroller
			this.scroller = new _Scroller();
			this.scroller.grid = this;
			this.scroller.renderRow = lang.hitch(this, "renderRow");
			this.scroller.removeRow = lang.hitch(this, "rowRemoved");
		},

		createLayout: function(){
			// summary:
			//		Creates a new Grid layout
			this.layout = new this._layoutClass(this);
			this.connect(this.layout, "moveColumn", "onMoveColumn");
		},

		onMoveColumn: function(){
			this.update();
		},
		
		onResizeColumn: function(/*int*/ cellIdx){
			// Called when a column is resized.
		},

		// views
		createViews: function(){
			this.views = new _ViewManager(this);
			this.views.createView = lang.hitch(this, "createView");
		},

		createView: function(inClass, idx){
			var c = lang.getObject(inClass);
			var view = new c({ grid: this, index: idx });
			this.viewsNode.appendChild(view.domNode);
			this.viewsHeaderNode.appendChild(view.headerNode);
			this.views.addView(view);
			html.attr(this.domNode, "align", this.isLeftToRight() ? 'left' : 'right');
			return view;
		},

		buildViews: function(){
			for(var i=0, vs; (vs=this.layout.structure[i]); i++){
				this.createView(vs.type || dojox._scopeName + ".grid._View", i).setStructure(vs);
			}
			this.scroller.setContentNodes(this.views.getContentNodes());
		},

		_setStructureAttr: function(structure){
			var s = structure;
			if(s && lang.isString(s)){
				dojo.deprecated("dojox.grid._Grid.set('structure', 'objVar')", "use dojox.grid._Grid.set('structure', objVar) instead", "2.0");
				s=lang.getObject(s);
			}
			this.structure = s;
			if(!s){
				if(this.layout.structure){
					s = this.layout.structure;
				}else{
					return;
				}
			}
			this.views.destroyViews();
			this.focus.focusView = null;
			if(s !== this.layout.structure){
				this.layout.setStructure(s);
			}
			this._structureChanged();
		},

		setStructure: function(/* dojox.grid.__ViewDef|dojox.grid.__ViewDef[]|dojox.grid.__CellDef[]|Array[dojox.grid.__CellDef[]] */ inStructure){
			// summary:
			//		Install a new structure and rebuild the grid.
			dojo.deprecated("dojox.grid._Grid.setStructure(obj)", "use dojox.grid._Grid.set('structure', obj) instead.", "2.0");
			this._setStructureAttr(inStructure);
		},
		
		getColumnTogglingItems: function(){
			// summary:
			//		returns an array of dijit.CheckedMenuItem widgets that can be
			//		added to a menu for toggling columns on and off.
			var items, checkedItems = [];
			items = array.map(this.layout.cells, function(cell){
				if(!cell.menuItems){ cell.menuItems = []; }

				var self = this;
				var item = new CheckedMenuItem({
					label: cell.name,
					checked: !cell.hidden,
					_gridCell: cell,
					onChange: function(checked){
						if(self.layout.setColumnVisibility(this._gridCell.index, checked)){
							var items = this._gridCell.menuItems;
							if(items.length > 1){
								array.forEach(items, function(item){
									if(item !== this){
										item.setAttribute("checked", checked);
									}
								}, this);
							}
							checked = array.filter(self.layout.cells, function(c){
								if(c.menuItems.length > 1){
									array.forEach(c.menuItems, "item.set('disabled', false);");
								}else{
									c.menuItems[0].set('disabled', false);
								}
								return !c.hidden;
							});
							if(checked.length == 1){
								array.forEach(checked[0].menuItems, "item.set('disabled', true);");
							}
						}
					},
					destroy: function(){
						var index = array.indexOf(this._gridCell.menuItems, this);
						this._gridCell.menuItems.splice(index, 1);
						delete this._gridCell;
						CheckedMenuItem.prototype.destroy.apply(this, arguments);
					}
				});
				cell.menuItems.push(item);
				if(!cell.hidden) {
					checkedItems.push(item);
				}
				return item;
			}, this); // dijit.CheckedMenuItem[]
			if(checkedItems.length == 1) {
				checkedItems[0].set('disabled', true);
			}
			return items;
		},

		_setHeaderMenuAttr: function(menu){
			if(this._placeholders && this._placeholders.length){
				array.forEach(this._placeholders, function(p){
					p.unReplace(true);
				});
				this._placeholders = [];
			}
			if(this.headerMenu){
				this.headerMenu.unBindDomNode(this.viewsHeaderNode);
			}
			this.headerMenu = menu;
			if(!menu){ return; }

			this.headerMenu.bindDomNode(this.viewsHeaderNode);
			if(this.headerMenu.getPlaceholders){
				this._placeholders = this.headerMenu.getPlaceholders(this.placeholderLabel);
			}
		},

		setHeaderMenu: function(/* dijit.Menu */ menu){
			dojo.deprecated("dojox.grid._Grid.setHeaderMenu(obj)", "use dojox.grid._Grid.set('headerMenu', obj) instead.", "2.0");
			this._setHeaderMenuAttr(menu);
		},
		
		setupHeaderMenu: function(){
			if(this._placeholders && this._placeholders.length){
				array.forEach(this._placeholders, function(p){
					if(p._replaced){
						p.unReplace(true);
					}
					p.replace(this.getColumnTogglingItems());
				}, this);
			}
		},

		_fetch: function(start){
			this.setScrollTop(0);
		},

		getItem: function(inRowIndex){
			return null;
		},
		
		showMessage: function(message){
			if(message){
				this.messagesNode.innerHTML = message;
				this.messagesNode.style.display = "";
			}else{
				this.messagesNode.innerHTML = "";
				this.messagesNode.style.display = "none";
			}
		},

		_structureChanged: function() {
			this.buildViews();
			if(this.autoRender && this._started){
				this.render();
			}
		},

		hasLayout: function() {
			return this.layout.cells.length;
		},

		// sizing
		resize: function(changeSize, resultSize){
			// summary:
			//		Update the grid's rendering dimensions and resize it
			
			// Calling sizeChange calls update() which calls _resize...so let's
			// save our input values, if any, and use them there when it gets
			// called.  This saves us an extra call to _resize(), which can
			// get kind of heavy.
			
			this._pendingChangeSize = changeSize;
			this._pendingResultSize = resultSize;
			this.sizeChange();
		},

		_getPadBorder: function() {
			this._padBorder = this._padBorder || html._getPadBorderExtents(this.domNode);
			return this._padBorder;
		},

		_getHeaderHeight: function(){
			var vns = this.viewsHeaderNode.style, t = vns.display == "none" ? 0 : this.views.measureHeader();
			vns.height = t + 'px';
			// header heights are reset during measuring so must be normalized after measuring.
			this.views.normalizeHeaderNodeHeight();
			return t;
		},
		
		_resize: function(changeSize, resultSize){
			// Restore our pending values, if any
			changeSize = changeSize || this._pendingChangeSize;
			resultSize = resultSize || this._pendingResultSize;
			delete this._pendingChangeSize;
			delete this._pendingResultSize;
			// if we have set up everything except the DOM, we cannot resize
			if(!this.domNode){ return; }
			var pn = this.domNode.parentNode;
			if(!pn || pn.nodeType != 1 || !this.hasLayout() || pn.style.visibility == "hidden" || pn.style.display == "none"){
				return;
			}
			// useful measurement
			var padBorder = this._getPadBorder();
			var hh = undefined;
			var h;
			// grid height
			if(this._autoHeight){
				this.domNode.style.height = 'auto';
			}else if(typeof this.autoHeight == "number"){
				h = hh = this._getHeaderHeight();
				h += (this.scroller.averageRowHeight * this.autoHeight);
				this.domNode.style.height = h + "px";
			}else if(this.domNode.clientHeight <= padBorder.h){
				if(pn == document.body){
					this.domNode.style.height = this.defaultHeight;
				}else if(this.height){
					this.domNode.style.height = this.height;
				}else{
					this.fitTo = "parent";
				}
			}
			// if we are given dimensions, size the grid's domNode to those dimensions
			if(resultSize){
				changeSize = resultSize;
			}
			if(!this._autoHeight && changeSize){
				html.marginBox(this.domNode, changeSize);
				this.height = this.domNode.style.height;
				delete this.fitTo;
			}else if(this.fitTo == "parent"){
				h = this._parentContentBoxHeight = (this._parentContentBoxHeight > 0 ? this._parentContentBoxHeight : html._getContentBox(pn).h);
				this.domNode.style.height = Math.max(0, h) + "px";
			}
			
			var hasFlex = array.some(this.views.views, function(v){ return v.flexCells; });

			if(!this._autoHeight && (h || html._getContentBox(this.domNode).h) === 0){
				// We need to hide the header, since the Grid is essentially hidden.
				this.viewsHeaderNode.style.display = "none";
			}else{
				// Otherwise, show the header and give it an appropriate height.
				this.viewsHeaderNode.style.display = "block";
				if(!hasFlex && hh === undefined){
					hh = this._getHeaderHeight();
				}
			}
			if(hasFlex){
				hh = undefined;
			}

			// NOTE: it is essential that width be applied before height
			// Header height can only be calculated properly after view widths have been set.
			// This is because flex column width is naturally 0 in Firefox.
			// Therefore prior to width sizing flex columns with spaces are maximally wrapped
			// and calculated to be too tall.
			this.adaptWidth();
			this.adaptHeight(hh);

			this.postresize();
		},

		adaptWidth: function() {
			// summary:
			//		sets width and position for views and update grid width if necessary
			// tags:
			//		private
			var doAutoWidth = (!this.initialWidth && this.autoWidth);
			var w = doAutoWidth ? 0 : this.domNode.clientWidth || (this.domNode.offsetWidth - this._getPadBorder().w),
				vw = this.views.arrange(1, w);
			this.views.onEach("adaptWidth");
			if(doAutoWidth){
				this.domNode.style.width = vw + "px";
			}
		},

		adaptHeight: function(inHeaderHeight){
			// summary:
			//		measures and normalizes header height, then sets view heights, and then updates scroller
			//		content extent
			// tags:
			//		private
			var t = inHeaderHeight === undefined ? this._getHeaderHeight() : inHeaderHeight;
			var h = (this._autoHeight ? -1 : Math.max(this.domNode.clientHeight - t, 0) || 0);
			this.views.onEach('setSize', [0, h]);
			this.views.onEach('adaptHeight');
			if(!this._autoHeight){
				var numScroll = 0, numNoScroll = 0;
				var noScrolls = array.filter(this.views.views, function(v){
					var has = v.hasHScrollbar();
					if(has){ numScroll++; }else{ numNoScroll++; }
					return (!has);
				});
				if(numScroll > 0 && numNoScroll > 0){
					array.forEach(noScrolls, function(v){
						v.adaptHeight(true);
					});
				}
			}
			if(this.autoHeight === true || h != -1 || (typeof this.autoHeight == "number" && this.autoHeight >= this.get('rowCount'))){
				this.scroller.windowHeight = h;
			}else{
				this.scroller.windowHeight = Math.max(this.domNode.clientHeight - t, 0);
			}
		},

		// startup
		startup: function(){
			if(this._started){return;}
			this.inherited(arguments);
			if(this.autoRender){
				this.render();
			}
		},

		// render
		render: function(){
			// summary:
			//	Render the grid, headers, and views. Edit and scrolling states are reset. To retain edit and
			//	scrolling states, see Update.

			if(!this.domNode){return;}
			if(!this._started){return;}

			if(!this.hasLayout()) {
				this.scroller.init(0, this.keepRows, this.rowsPerPage);
				return;
			}
			//
			this.update = this.defaultUpdate;
			this._render();
		},

		_render: function(){
			this.scroller.init(this.get('rowCount'), this.keepRows, this.rowsPerPage);
			this.prerender();
			this.setScrollTop(0);
			this.postrender();
		},

		prerender: function(){
			// if autoHeight, make sure scroller knows not to virtualize; everything must be rendered.
			this.keepRows = this._autoHeight ? 0 : this.keepRows;
			this.scroller.setKeepInfo(this.keepRows);
			this.views.render();
			this._resize();
		},

		postrender: function(){
			this.postresize();
			this.focus.initFocusView();
			// make rows unselectable
			html.setSelectable(this.domNode, this.selectable);
		},

		postresize: function(){
			// views are position absolute, so they do not inflate the parent
			if(this._autoHeight){
				var size = Math.max(this.views.measureContent()) + 'px';
				
				this.viewsNode.style.height = size;
			}
		},

		renderRow: function(inRowIndex, inNodes){
			// summary:
			//		used internally to render rows
			// tags:
			//		private
			this.views.renderRow(inRowIndex, inNodes, this._skipRowRenormalize);
		},

		rowRemoved: function(inRowIndex){
			// summary:
			//		used internally to remove rows
			// tags:
			//		private
			this.views.rowRemoved(inRowIndex);
		},

		invalidated: null,

		updating: false,

		beginUpdate: function(){
			// summary:
			//		Use to make multiple changes to rows while queueing row updating.
			// NOTE: not currently supporting nested begin/endUpdate calls
			this.invalidated = [];
			this.updating = true;
		},

		endUpdate: function(){
			// summary:
			//		Use after calling beginUpdate to render any changes made to rows.
			this.updating = false;
			var i = this.invalidated, r;
			if(i.all){
				this.update();
			}else if(i.rowCount != undefined){
				this.updateRowCount(i.rowCount);
			}else{
				for(r in i){
					this.updateRow(Number(r));
				}
			}
			this.invalidated = [];
		},

		// update
		defaultUpdate: function(){
			// note: initial update calls render and subsequently this function.
			if(!this.domNode){return;}
			if(this.updating){
				this.invalidated.all = true;
				return;
			}
			//this.edit.saveState(inRowIndex);
			this.lastScrollTop = this.scrollTop;
			this.prerender();
			this.scroller.invalidateNodes();
			this.setScrollTop(this.lastScrollTop);
			this.postrender();
			//this.edit.restoreState(inRowIndex);
		},

		update: function(){
			// summary:
			//		Update the grid, retaining edit and scrolling states.
			this.render();
		},

		updateRow: function(inRowIndex){
			// summary:
			//		Render a single row.
			// inRowIndex: Integer
			//		Index of the row to render
			inRowIndex = Number(inRowIndex);
			if(this.updating){
				this.invalidated[inRowIndex]=true;
			}else{
				this.views.updateRow(inRowIndex);
				this.scroller.rowHeightChanged(inRowIndex);
			}
		},

		updateRows: function(startIndex, howMany){
			// summary:
			//		Render consecutive rows at once.
			// startIndex: Integer
			//		Index of the starting row to render
			// howMany: Integer
			//		How many rows to update.
			startIndex = Number(startIndex);
			howMany = Number(howMany);
			var i;
			if(this.updating){
				for(i=0; i<howMany; i++){
					this.invalidated[i+startIndex]=true;
				}
			}else{
				for(i=0; i<howMany; i++){
					this.views.updateRow(i+startIndex, this._skipRowRenormalize);
				}
				this.scroller.rowHeightChanged(startIndex);
			}
		},

		updateRowCount: function(inRowCount){
			// summary:
			//		Change the number of rows.
			// inRowCount: int
			//		Number of rows in the grid.
			if(this.updating){
				this.invalidated.rowCount = inRowCount;
			}else{
				this.rowCount = inRowCount;
				this._setAutoHeightAttr(this.autoHeight, true);
				if(this.layout.cells.length){
					this.scroller.updateRowCount(inRowCount);
				}
				this._resize();
				if(this.layout.cells.length){
					this.setScrollTop(this.scrollTop);
				}
			}
		},

		updateRowStyles: function(inRowIndex){
			// summary:
			//		Update the styles for a row after it's state has changed.
			this.views.updateRowStyles(inRowIndex);
		},
		getRowNode: function(inRowIndex){
			// summary:
			//		find the rowNode that is not a rowSelector
			if (this.focus.focusView && !(this.focus.focusView instanceof _RowSelector)){
					return this.focus.focusView.rowNodes[inRowIndex];
			}else{ // search through views
				for (var i = 0, cView; (cView = this.views.views[i]); i++) {
					if (!(cView instanceof _RowSelector)) {
						return cView.rowNodes[inRowIndex];
					}
				}
			}
			return null;
		},
		rowHeightChanged: function(inRowIndex){
			// summary:
			//		Update grid when the height of a row has changed. Row height is handled automatically as rows
			//		are rendered. Use this function only to update a row's height outside the normal rendering process.
			// inRowIndex: Integer
			//		index of the row that has changed height

			this.views.renormalizeRow(inRowIndex);
			this.scroller.rowHeightChanged(inRowIndex);
		},

		// fastScroll: Boolean
		//		flag modifies vertical scrolling behavior. Defaults to true but set to false for slower
		//		scroll performance but more immediate scrolling feedback
		fastScroll: true,

		delayScroll: false,

		// scrollRedrawThreshold: int
		//	pixel distance a user must scroll vertically to trigger grid scrolling.
		scrollRedrawThreshold: (has('ie') ? 100 : 50),

		// scroll methods
		scrollTo: function(inTop){
			// summary:
			//		Vertically scroll the grid to a given pixel position
			// inTop: Integer
			//		vertical position of the grid in pixels
			if(!this.fastScroll){
				this.setScrollTop(inTop);
				return;
			}
			var delta = Math.abs(this.lastScrollTop - inTop);
			this.lastScrollTop = inTop;
			if(delta > this.scrollRedrawThreshold || this.delayScroll){
				this.delayScroll = true;
				this.scrollTop = inTop;
				this.views.setScrollTop(inTop);
				if(this._pendingScroll){
					window.clearTimeout(this._pendingScroll);
				}
				var _this = this;
				this._pendingScroll = window.setTimeout(function(){
					delete _this._pendingScroll;
					_this.finishScrollJob();
				}, 200);
			}else{
				this.setScrollTop(inTop);
			}
		},

		finishScrollJob: function(){
			this.delayScroll = false;
			this.setScrollTop(this.scrollTop);
		},

		setScrollTop: function(inTop){
			this.scroller.scroll(this.views.setScrollTop(inTop));
		},

		scrollToRow: function(inRowIndex){
			// summary:
			//		Scroll the grid to a specific row.
			// inRowIndex: Integer
			//		grid row index
			this.setScrollTop(this.scroller.findScrollTop(inRowIndex) + 1);
		},

		styleRowNode: function(inRowIndex, inRowNode){
			// summary:
			//		styling (used internally to style individual parts of a row)
			// tags:
			//		private
			if(inRowNode){
				this.rows.styleRowNode(inRowIndex, inRowNode);
			}
		},
		
		// called when the mouse leaves the grid so we can deselect all hover rows
		_mouseOut: function(e){
			this.rows.setOverRow(-2);
		},
	
		// cells
		getCell: function(inIndex){
			// summary:
			//		Retrieves the cell object for a given grid column.
			// inIndex: Integer
			//		Grid column index of cell to retrieve
			// returns:
			//		a grid cell
			return this.layout.cells[inIndex];
		},

		setCellWidth: function(inIndex, inUnitWidth){
			this.getCell(inIndex).unitWidth = inUnitWidth;
		},

		getCellName: function(inCell){
			// summary:
			//		Returns the cell name of a passed cell
			return "Cell " + inCell.index; // String
		},

		// sorting
		canSort: function(inSortInfo){
			// summary:
			//		Determines if the grid can be sorted
			// inSortInfo: Integer
			//		Sort information, 1-based index of column on which to sort, positive for an ascending sort
			//		and negative for a descending sort
			// returns: Boolean
			//		True if grid can be sorted on the given column in the given direction
		},

		sort: function(){
		},

		getSortAsc: function(inSortInfo){
			// summary:
			//		Returns true if grid is sorted in an ascending direction.
			inSortInfo = inSortInfo == undefined ? this.sortInfo : inSortInfo;
			return Boolean(inSortInfo > 0); // Boolean
		},

		getSortIndex: function(inSortInfo){
			// summary:
			//		Returns the index of the column on which the grid is sorted
			inSortInfo = inSortInfo == undefined ? this.sortInfo : inSortInfo;
			return Math.abs(inSortInfo) - 1; // Integer
		},

		setSortIndex: function(inIndex, inAsc){
			// summary:
			//		Sort the grid on a column in a specified direction
			// inIndex: Integer
			//		Column index on which to sort.
			// inAsc: Boolean
			//		If true, sort the grid in ascending order, otherwise in descending order
			var si = inIndex +1;
			if(inAsc != undefined){
				si *= (inAsc ? 1 : -1);
			} else if(this.getSortIndex() == inIndex){
				si = -this.sortInfo;
			}
			this.setSortInfo(si);
		},

		setSortInfo: function(inSortInfo){
			if(this.canSort(inSortInfo)){
				this.sortInfo = inSortInfo;
				this.sort();
				this.update();
			}
		},

		// DOM event handler
		doKeyEvent: function(e){
			e.dispatch = 'do' + e.type;
			this.onKeyEvent(e);
		},

		// event dispatch
		//: protected
		_dispatch: function(m, e){
			if(m in this){
				return this[m](e);
			}
			return false;
		},

		dispatchKeyEvent: function(e){
			this._dispatch(e.dispatch, e);
		},

		dispatchContentEvent: function(e){
			this.edit.dispatchEvent(e) || e.sourceView.dispatchContentEvent(e) || this._dispatch(e.dispatch, e);
		},

		dispatchHeaderEvent: function(e){
			e.sourceView.dispatchHeaderEvent(e) || this._dispatch('doheader' + e.type, e);
		},

		dokeydown: function(e){
			this.onKeyDown(e);
		},

		doclick: function(e){
			if(e.cellNode){
				this.onCellClick(e);
			}else{
				this.onRowClick(e);
			}
		},

		dodblclick: function(e){
			if(e.cellNode){
				this.onCellDblClick(e);
			}else{
				this.onRowDblClick(e);
			}
		},

		docontextmenu: function(e){
			if(e.cellNode){
				this.onCellContextMenu(e);
			}else{
				this.onRowContextMenu(e);
			}
		},

		doheaderclick: function(e){
			if(e.cellNode){
				this.onHeaderCellClick(e);
			}else{
				this.onHeaderClick(e);
			}
		},

		doheaderdblclick: function(e){
			if(e.cellNode){
				this.onHeaderCellDblClick(e);
			}else{
				this.onHeaderDblClick(e);
			}
		},

		doheadercontextmenu: function(e){
			if(e.cellNode){
				this.onHeaderCellContextMenu(e);
			}else{
				this.onHeaderContextMenu(e);
			}
		},

		// override to modify editing process
		doStartEdit: function(inCell, inRowIndex){
			this.onStartEdit(inCell, inRowIndex);
		},

		doApplyCellEdit: function(inValue, inRowIndex, inFieldIndex){
			this.onApplyCellEdit(inValue, inRowIndex, inFieldIndex);
		},

		doCancelEdit: function(inRowIndex){
			this.onCancelEdit(inRowIndex);
		},

		doApplyEdit: function(inRowIndex){
			this.onApplyEdit(inRowIndex);
		},

		// row editing
		addRow: function(){
			// summary:
			//		Add a row to the grid.
			this.updateRowCount(this.get('rowCount')+1);
		},

		removeSelectedRows: function(){
			// summary:
			//		Remove the selected rows from the grid.
			if(this.allItemsSelected){
				this.updateRowCount(0);
			}else{
				this.updateRowCount(Math.max(0, this.get('rowCount') - this.selection.getSelected().length));
			}
			this.selection.clear();
		}

	});

	_Grid.markupFactory = function(props, node, ctor, cellFunc){
		var widthFromAttr = function(n){
			var w = html.attr(n, "width")||"auto";
			if((w != "auto")&&(w.slice(-2) != "em")&&(w.slice(-1) != "%")){
				w = parseInt(w, 10)+"px";
			}
			return w;
		};
		// if(!props.store){ console.debug("no store!"); }
		// if a structure isn't referenced, do we have enough
		// data to try to build one automatically?
		if(	!props.structure &&
			node.nodeName.toLowerCase() == "table"){

			// try to discover a structure
			props.structure = query("> colgroup", node).map(function(cg){
				var sv = html.attr(cg, "span");
				var v = {
					noscroll: (html.attr(cg, "noscroll") == "true") ? true : false,
					__span: (!!sv ? parseInt(sv, 10) : 1),
					cells: []
				};
				if(html.hasAttr(cg, "width")){
					v.width = widthFromAttr(cg);
				}
				return v; // for vendetta
			});
			if(!props.structure.length){
				props.structure.push({
					__span: Infinity,
					cells: [] // catch-all view
				});
			}
			// check to see if we're gonna have more than one view

			// for each tr in our th, create a row of cells
			query("thead > tr", node).forEach(function(tr, tr_idx){
				var cellCount = 0;
				var viewIdx = 0;
				var lastViewIdx;
				var cView = null;
				query("> th", tr).map(function(th){
					// what view will this cell go into?

					// NOTE:
					//		to prevent extraneous iteration, we start counters over
					//		for each row, incrementing over the surface area of the
					//		structure that colgroup processing generates and
					//		creating cell objects for each <th> to place into those
					//		cell groups.  There's a lot of state-keepking logic
					//		here, but it is what it has to be.
					if(!cView){ // current view book keeping
						lastViewIdx = 0;
						cView = props.structure[0];
					}else if(cellCount >= (lastViewIdx+cView.__span)){
						viewIdx++;
						// move to allocating things into the next view
						lastViewIdx += cView.__span;
						var lastView = cView;
						cView = props.structure[viewIdx];
					}

					// actually define the cell from what markup hands us
					var cell = {
						name: lang.trim(html.attr(th, "name")||th.innerHTML),
						colSpan: parseInt(html.attr(th, "colspan")||1, 10),
						type: lang.trim(html.attr(th, "cellType")||""),
						id: lang.trim(html.attr(th,"id")||"")
					};
					cellCount += cell.colSpan;
					var rowSpan = html.attr(th, "rowspan");
					if(rowSpan){
						cell.rowSpan = rowSpan;
					}
					if(html.hasAttr(th, "width")){
						cell.width = widthFromAttr(th);
					}
					if(html.hasAttr(th, "relWidth")){
						cell.relWidth = window.parseInt(html.attr(th, "relWidth"), 10);
					}
					if(html.hasAttr(th, "hidden")){
						cell.hidden = (html.attr(th, "hidden") == "true" || html.attr(th, "hidden") === true/*always boolean true in Chrome*/);
					}

					if(cellFunc){
						cellFunc(th, cell);
					}

					cell.type = cell.type ? lang.getObject(cell.type) : dojox.grid.cells.Cell;

					if(cell.type && cell.type.markupFactory){
						cell.type.markupFactory(th, cell);
					}

					if(!cView.cells[tr_idx]){
						cView.cells[tr_idx] = [];
					}
					cView.cells[tr_idx].push(cell);
				});
			});
		}

		return new ctor(props, node);
	};

	return _Grid;

});

},
'dojox/grid/_Events':function(){
define("dojox/grid/_Events", [
	"dojo/keys",
	"dojo/dom-class",
	"dojo/_base/declare",
	"dojo/_base/event",
	"dojo/_base/sniff"
], function(keys, domClass, declare, event, has){

return declare("dojox.grid._Events", null, {
	// summary:
	//		_Grid mixin that provides default implementations for grid events.
	// description:
	//		Default synthetic events dispatched for _Grid. dojo.connect to events to
	//		retain default implementation or override them for custom handling.
	
	// cellOverClass: String
	//		css class to apply to grid cells over which the cursor is placed.
	cellOverClass: "dojoxGridCellOver",
	
	onKeyEvent: function(e){
		// summary:
		//		top level handler for Key Events
		this.dispatchKeyEvent(e);
	},

	onContentEvent: function(e){
		// summary:
		//		Top level handler for Content events
		this.dispatchContentEvent(e);
	},

	onHeaderEvent: function(e){
		// summary:
		//		Top level handler for header events
		this.dispatchHeaderEvent(e);
	},

	onStyleRow: function(inRow){
		// summary:
		//		Perform row styling on a given row. Called whenever row styling is updated.
		// inRow: Object
		//		Object containing row state information: selected, true if the row is selcted; over:
		//		true of the mouse is over the row; odd: true if the row is odd. Use customClasses and
		//		customStyles to control row css classes and styles; both properties are strings.
		// example:
		// |	onStyleRow({ selected: true, over:true, odd:false })
		var i = inRow;
		i.customClasses += (i.odd?" dojoxGridRowOdd":"") + (i.selected?" dojoxGridRowSelected":"") + (i.over?" dojoxGridRowOver":"");
		this.focus.styleRow(inRow);
		this.edit.styleRow(inRow);
	},
	
	onKeyDown: function(e){
		// summary:
		//		Grid key event handler. By default enter begins editing and applies edits, escape cancels an edit,
		//		tab, shift-tab, and arrow keys move grid cell focus.
		if(e.altKey || e.metaKey){
			return;
		}
		var colIdx;
		switch(e.keyCode){
			case keys.ESCAPE:
				this.edit.cancel();
				break;
			case keys.ENTER:
				if(!this.edit.isEditing()){
					colIdx = this.focus.getHeaderIndex();
					if(colIdx >= 0) {
						this.setSortIndex(colIdx);
						break;
					}else {
						this.selection.clickSelect(this.focus.rowIndex, dojo.isCopyKey(e), e.shiftKey);
					}
					event.stop(e);
				}
				if(!e.shiftKey){
					var isEditing = this.edit.isEditing();
					this.edit.apply();
					if(!isEditing){
						this.edit.setEditCell(this.focus.cell, this.focus.rowIndex);
					}
				}
				if (!this.edit.isEditing()){
					var curView = this.focus.focusView || this.views.views[0];  //if no focusView than only one view
					curView.content.decorateEvent(e);
					this.onRowClick(e);
					event.stop(e);
				}
				break;
			case keys.SPACE:
				if(!this.edit.isEditing()){
					colIdx = this.focus.getHeaderIndex();
					if(colIdx >= 0) {
						this.setSortIndex(colIdx);
						break;
					}else {
						this.selection.clickSelect(this.focus.rowIndex, dojo.isCopyKey(e), e.shiftKey);
					}
					event.stop(e);
				}
				break;
			case keys.TAB:
				this.focus[e.shiftKey ? 'previousKey' : 'nextKey'](e);
				break;
			case keys.LEFT_ARROW:
			case keys.RIGHT_ARROW:
				if(!this.edit.isEditing()){
					var keyCode = e.keyCode;  // IE seems to lose after stopEvent when modifier keys
					event.stop(e);
					colIdx = this.focus.getHeaderIndex();
					if (colIdx >= 0 && (e.shiftKey && e.ctrlKey)){
						this.focus.colSizeAdjust(e, colIdx, (keyCode == keys.LEFT_ARROW ? -1 : 1)*5);
					}
					else{
						var offset = (keyCode == keys.LEFT_ARROW) ? 1 : -1;
						if(this.isLeftToRight()){ offset *= -1; }
						this.focus.move(0, offset);
					}
				}
				break;
			case keys.UP_ARROW:
				if(!this.edit.isEditing() && this.focus.rowIndex !== 0){
					event.stop(e);
					this.focus.move(-1, 0);
				}
				break;
			case keys.DOWN_ARROW:
				if(!this.edit.isEditing() && this.focus.rowIndex+1 != this.rowCount){
					event.stop(e);
					this.focus.move(1, 0);
				}
				break;
			case keys.PAGE_UP:
				if(!this.edit.isEditing() && this.focus.rowIndex !== 0){
					event.stop(e);
					if(this.focus.rowIndex != this.scroller.firstVisibleRow+1){
						this.focus.move(this.scroller.firstVisibleRow-this.focus.rowIndex, 0);
					}else{
						this.setScrollTop(this.scroller.findScrollTop(this.focus.rowIndex-1));
						this.focus.move(this.scroller.firstVisibleRow-this.scroller.lastVisibleRow+1, 0);
					}
				}
				break;
			case keys.PAGE_DOWN:
				if(!this.edit.isEditing() && this.focus.rowIndex+1 != this.rowCount){
					event.stop(e);
					if(this.focus.rowIndex != this.scroller.lastVisibleRow-1){
						this.focus.move(this.scroller.lastVisibleRow-this.focus.rowIndex-1, 0);
					}else{
						this.setScrollTop(this.scroller.findScrollTop(this.focus.rowIndex+1));
						this.focus.move(this.scroller.lastVisibleRow-this.scroller.firstVisibleRow-1, 0);
					}
				}
				break;
			default:
				break;
		}
	},
	
	onMouseOver: function(e){
		// summary:
		//		Event fired when mouse is over the grid.
		// e: Event
		//		Decorated event object contains reference to grid, cell, and rowIndex
		e.rowIndex == -1 ? this.onHeaderCellMouseOver(e) : this.onCellMouseOver(e);
	},
	
	onMouseOut: function(e){
		// summary:
		//		Event fired when mouse moves out of the grid.
		// e: Event
		//		Decorated event object that contains reference to grid, cell, and rowIndex
		e.rowIndex == -1 ? this.onHeaderCellMouseOut(e) : this.onCellMouseOut(e);
	},
	
	onMouseDown: function(e){
		// summary:
		//		Event fired when mouse is down inside grid.
		// e: Event
		//		Decorated event object that contains reference to grid, cell, and rowIndex
		e.rowIndex == -1 ? this.onHeaderCellMouseDown(e) : this.onCellMouseDown(e);
	},
	
	onMouseOverRow: function(e){
		// summary:
		//		Event fired when mouse is over any row (data or header).
		// e: Event
		//		Decorated event object contains reference to grid, cell, and rowIndex
		if(!this.rows.isOver(e.rowIndex)){
			this.rows.setOverRow(e.rowIndex);
			e.rowIndex == -1 ? this.onHeaderMouseOver(e) : this.onRowMouseOver(e);
		}
	},
	onMouseOutRow: function(e){
		// summary:
		//		Event fired when mouse moves out of any row (data or header).
		// e: Event
		//		Decorated event object contains reference to grid, cell, and rowIndex
		if(this.rows.isOver(-1)){
			this.onHeaderMouseOut(e);
		}else if(!this.rows.isOver(-2)){
			this.rows.setOverRow(-2);
			this.onRowMouseOut(e);
		}
	},
	
	onMouseDownRow: function(e){
		// summary:
		//		Event fired when mouse is down inside grid row
		// e: Event
		//		Decorated event object that contains reference to grid, cell, and rowIndex
		if(e.rowIndex != -1)
			this.onRowMouseDown(e);
	},

	// cell events
	onCellMouseOver: function(e){
		// summary:
		//		Event fired when mouse is over a cell.
		// e: Event
		//		Decorated event object contains reference to grid, cell, and rowIndex
		if(e.cellNode){
			domClass.add(e.cellNode, this.cellOverClass);
		}
	},
	
	onCellMouseOut: function(e){
		// summary:
		//		Event fired when mouse moves out of a cell.
		// e: Event
		//		Decorated event object which contains reference to grid, cell, and rowIndex
		if(e.cellNode){
			domClass.remove(e.cellNode, this.cellOverClass);
		}
	},
	
	onCellMouseDown: function(e){
		// summary:
		//		Event fired when mouse is down in a header cell.
		// e: Event
		//		Decorated event object which contains reference to grid, cell, and rowIndex
	},

	onCellClick: function(e){
		// summary:
		//		Event fired when a cell is clicked.
		// e: Event
		//		Decorated event object which contains reference to grid, cell, and rowIndex
		this._click[0] = this._click[1];
		this._click[1] = e;
		if(!this.edit.isEditCell(e.rowIndex, e.cellIndex)){
			this.focus.setFocusCell(e.cell, e.rowIndex);
		}
		// in some cases click[0] is null which causes false doubeClicks. Fixes #100703
		if(this._click.length > 1 && this._click[0] == null){
			this._click.shift();
		}
		this.onRowClick(e);
	},

	onCellDblClick: function(e){
		// summary:
		//		Event fired when a cell is double-clicked.
		// e: Event
		//		Decorated event object contains reference to grid, cell, and rowIndex
		var event;
		if(this._click.length > 1 && has('ie')){
			event = this._click[1];
		}else if(this._click.length > 1 && this._click[0].rowIndex != this._click[1].rowIndex){
			event = this._click[0];
		}else{
			event = e;
		}
		this.focus.setFocusCell(event.cell, event.rowIndex);
		this.onRowClick(event);
		this.edit.setEditCell(event.cell, event.rowIndex);
		this.onRowDblClick(e);
	},

	onCellContextMenu: function(e){
		// summary:
		//		Event fired when a cell context menu is accessed via mouse right click.
		// e: Event
		//		Decorated event object which contains reference to grid, cell, and rowIndex
		this.onRowContextMenu(e);
	},

	onCellFocus: function(inCell, inRowIndex){
		// summary:
		//		Event fired when a cell receives focus.
		// inCell: Object
		//		Cell object containing properties of the grid column.
		// inRowIndex: Integer
		//		Index of the grid row
		this.edit.cellFocus(inCell, inRowIndex);
	},

	// row events
	onRowClick: function(e){
		// summary:
		//		Event fired when a row is clicked.
		// e: Event
		//		Decorated event object which contains reference to grid, cell, and rowIndex
		this.edit.rowClick(e);
		this.selection.clickSelectEvent(e);
	},

	onRowDblClick: function(e){
		// summary:
		//		Event fired when a row is double clicked.
		// e: Event
		//		decorated event object which contains reference to grid, cell, and rowIndex
	},

	onRowMouseOver: function(e){
		// summary:
		//		Event fired when mouse moves over a data row.
		// e: Event
		//		Decorated event object which contains reference to grid, cell, and rowIndex
	},

	onRowMouseOut: function(e){
		// summary:
		//		Event fired when mouse moves out of a data row.
		// e: Event
		//		Decorated event object contains reference to grid, cell, and rowIndex
	},
	
	onRowMouseDown: function(e){
		// summary:
		//		Event fired when mouse is down in a row.
		// e: Event
		//		Decorated event object which contains reference to grid, cell, and rowIndex
	},

	onRowContextMenu: function(e){
		// summary:
		//		Event fired when a row context menu is accessed via mouse right click.
		// e: Event
		//		Decorated event object which contains reference to grid, cell, and rowIndex
		event.stop(e);
	},

	// header events
	onHeaderMouseOver: function(e){
		// summary:
		//		Event fired when mouse moves over the grid header.
		// e: Event
		//		Decorated event object contains reference to grid, cell, and rowIndex
	},

	onHeaderMouseOut: function(e){
		// summary:
		//		Event fired when mouse moves out of the grid header.
		// e: Event
		//		Decorated event object which contains reference to grid, cell, and rowIndex
	},

	onHeaderCellMouseOver: function(e){
		// summary:
		//		Event fired when mouse moves over a header cell.
		// e: Event
		//		Decorated event object which contains reference to grid, cell, and rowIndex
		if(e.cellNode){
			domClass.add(e.cellNode, this.cellOverClass);
		}
	},

	onHeaderCellMouseOut: function(e){
		// summary:
		//		Event fired when mouse moves out of a header cell.
		// e: Event
		//		Decorated event object which contains reference to grid, cell, and rowIndex
		if(e.cellNode){
			domClass.remove(e.cellNode, this.cellOverClass);
		}
	},
	
	onHeaderCellMouseDown: function(e) {
		// summary:
		//		Event fired when mouse is down in a header cell.
		// e: Event
		//		Decorated event object which contains reference to grid, cell, and rowIndex
	},

	onHeaderClick: function(e){
		// summary:
		//		Event fired when the grid header is clicked.
		// e: Event
		//		Decorated event object which contains reference to grid, cell, and rowIndex
	},

	onHeaderCellClick: function(e){
		// summary:
		//		Event fired when a header cell is clicked.
		// e: Event
		//		Decorated event object which contains reference to grid, cell, and rowIndex
		this.setSortIndex(e.cell.index);
		this.onHeaderClick(e);
	},

	onHeaderDblClick: function(e){
		// summary:
		//		Event fired when the grid header is double clicked.
		// e: Event
		//		Decorated event object which contains reference to grid, cell, and rowIndex
	},

	onHeaderCellDblClick: function(e){
		// summary:
		//		Event fired when a header cell is double clicked.
		// e: Event
		//		Decorated event object which contains reference to grid, cell, and rowIndex
		this.onHeaderDblClick(e);
	},

	onHeaderCellContextMenu: function(e){
		// summary:
		//		Event fired when a header cell context menu is accessed via mouse right click.
		// e: Event
		//		Decorated event object which contains reference to grid, cell, and rowIndex
		this.onHeaderContextMenu(e);
	},

	onHeaderContextMenu: function(e){
		// summary:
		//		Event fired when the grid header context menu is accessed via mouse right click.
		// e: Event
		//		Decorated event object which contains reference to grid, cell, and rowIndex
		if(!this.headerMenu){
			event.stop(e);
		}
	},

	// editing
	onStartEdit: function(inCell, inRowIndex){
		// summary:
		//		Event fired when editing is started for a given grid cell
		// inCell: Object
		//		Cell object containing properties of the grid column.
		// inRowIndex: Integer
		//		Index of the grid row
	},

	onApplyCellEdit: function(inValue, inRowIndex, inFieldIndex){
		// summary:
		//		Event fired when editing is applied for a given grid cell
		// inValue: String
		//		Value from cell editor
		// inRowIndex: Integer
		//		Index of the grid row
		// inFieldIndex: Integer
		//		Index in the grid's data store
	},

	onCancelEdit: function(inRowIndex){
		// summary:
		//		Event fired when editing is cancelled for a given grid cell
		// inRowIndex: Integer
		//		Index of the grid row
	},

	onApplyEdit: function(inRowIndex){
		// summary:
		//		Event fired when editing is applied for a given grid row
		// inRowIndex: Integer
		//		Index of the grid row
	},

	onCanSelect: function(inRowIndex){
		// summary:
		//		Event to determine if a grid row may be selected
		// inRowIndex: Integer
		//		Index of the grid row
		// returns: Boolean
		//		true if the row can be selected
		return true;
	},

	onCanDeselect: function(inRowIndex){
		// summary:
		//		Event to determine if a grid row may be deselected
		// inRowIndex: Integer
		//		Index of the grid row
		// returns: Boolean
		//		true if the row can be deselected
		return true;
	},

	onSelected: function(inRowIndex){
		// summary:
		//		Event fired when a grid row is selected
		// inRowIndex: Integer
		//		Index of the grid row
		this.updateRowStyles(inRowIndex);
	},

	onDeselected: function(inRowIndex){
		// summary:
		//		Event fired when a grid row is deselected
		// inRowIndex: Integer
		//		Index of the grid row
		this.updateRowStyles(inRowIndex);
	},

	onSelectionChanged: function(){
	}
});
});
},
'dojox/grid/_Scroller':function(){
define("dojox/grid/_Scroller", [
	"dijit/registry",
	"dojo/_base/declare",
	"dojo/_base/lang",
	"./util",
	"dojo/_base/html"
], function(dijitRegistry, declare, lang, util, html){

	var indexInParent = function(inNode){
		var i=0, n, p=inNode.parentNode;
		while((n = p.childNodes[i++])){
			if(n == inNode){
				return i - 1;
			}
		}
		return -1;
	};
	
	var cleanNode = function(inNode){
		if(!inNode){
			return;
		}
		dojo.forEach(dijitRegistry.toArray(), function(w){
			if(w.domNode && html.isDescendant(w.domNode, inNode, true)){
				w.destroy();
			}
		});
	};

	var getTagName = function(inNodeOrId){
		var node = html.byId(inNodeOrId);
		return (node && node.tagName ? node.tagName.toLowerCase() : '');
	};
	
	var nodeKids = function(inNode, inTag){
		var result = [];
		var i=0, n;
		while((n = inNode.childNodes[i])){
			i++;
			if(getTagName(n) == inTag){
				result.push(n);
			}
		}
		return result;
	};
	
	var divkids = function(inNode){
		return nodeKids(inNode, 'div');
	};

	return declare("dojox.grid._Scroller", null, {
		constructor: function(inContentNodes){
			this.setContentNodes(inContentNodes);
			this.pageHeights = [];
			this.pageNodes = [];
			this.stack = [];
		},
		// specified
		rowCount: 0, // total number of rows to manage
		defaultRowHeight: 32, // default height of a row
		keepRows: 100, // maximum number of rows that should exist at one time
		contentNode: null, // node to contain pages
		scrollboxNode: null, // node that controls scrolling
		// calculated
		defaultPageHeight: 0, // default height of a page
		keepPages: 10, // maximum number of pages that should exists at one time
		pageCount: 0,
		windowHeight: 0,
		firstVisibleRow: 0,
		lastVisibleRow: 0,
		averageRowHeight: 0, // the average height of a row
		// private
		page: 0,
		pageTop: 0,
		// init
		init: function(inRowCount, inKeepRows, inRowsPerPage){
			switch(arguments.length){
				case 3: this.rowsPerPage = inRowsPerPage;
				case 2: this.keepRows = inKeepRows;
				case 1: this.rowCount = inRowCount;
				default: break;
			}
			this.defaultPageHeight = (this.grid.rowHeight > 0 ? this.grid.rowHeight : this.defaultRowHeight) * this.rowsPerPage;
			this.pageCount = this._getPageCount(this.rowCount, this.rowsPerPage);
			this.setKeepInfo(this.keepRows);
			this.invalidate();
			if(this.scrollboxNode){
				this.scrollboxNode.scrollTop = 0;
				this.scroll(0);
				this.scrollboxNode.onscroll = lang.hitch(this, 'onscroll');
			}
		},
		_getPageCount: function(rowCount, rowsPerPage){
			return rowCount ? (Math.ceil(rowCount / rowsPerPage) || 1) : 0;
		},
		destroy: function(){
			this.invalidateNodes();
			delete this.contentNodes;
			delete this.contentNode;
			delete this.scrollboxNode;
		},
		setKeepInfo: function(inKeepRows){
			this.keepRows = inKeepRows;
			this.keepPages = !this.keepRows ? this.keepPages : Math.max(Math.ceil(this.keepRows / this.rowsPerPage), 2);
		},
		// nodes
		setContentNodes: function(inNodes){
			this.contentNodes = inNodes;
			this.colCount = (this.contentNodes ? this.contentNodes.length : 0);
			this.pageNodes = [];
			for(var i=0; i<this.colCount; i++){
				this.pageNodes[i] = [];
			}
		},
		getDefaultNodes: function(){
			return this.pageNodes[0] || [];
		},
		// updating
		invalidate: function(){
			this._invalidating = true;
			this.invalidateNodes();
			this.pageHeights = [];
			this.height = (this.pageCount ? (this.pageCount - 1)* this.defaultPageHeight + this.calcLastPageHeight() : 0);
			this.resize();
			this._invalidating = false;
		},
		updateRowCount: function(inRowCount){
			this.invalidateNodes();
			this.rowCount = inRowCount;
			// update page count, adjust document height
			var oldPageCount = this.pageCount;
			if(oldPageCount === 0){
				//We want to have at least 1px in height to keep scroller.  Otherwise with an
				//empty grid you can't scroll to see the header.
				this.height = 1;
			}
			this.pageCount = this._getPageCount(this.rowCount, this.rowsPerPage);
			if(this.pageCount < oldPageCount){
				for(var i=oldPageCount-1; i>=this.pageCount; i--){
					this.height -= this.getPageHeight(i);
					delete this.pageHeights[i];
				}
			}else if(this.pageCount > oldPageCount){
				this.height += this.defaultPageHeight * (this.pageCount - oldPageCount - 1) + this.calcLastPageHeight();
			}
			this.resize();
		},
		// implementation for page manager
		pageExists: function(inPageIndex){
			return Boolean(this.getDefaultPageNode(inPageIndex));
		},
		measurePage: function(inPageIndex){
			if(this.grid.rowHeight){
				return ((inPageIndex + 1) * this.rowsPerPage > this.rowCount ?
					this.rowCount - inPageIndex * this.rowsPerPage :
					this.rowsPerPage) * this.grid.rowHeight;
					 
			}
			var n = this.getDefaultPageNode(inPageIndex);
			return (n && n.innerHTML) ? n.offsetHeight : undefined;
		},
		positionPage: function(inPageIndex, inPos){
			for(var i=0; i<this.colCount; i++){
				this.pageNodes[i][inPageIndex].style.top = inPos + 'px';
			}
		},
		repositionPages: function(inPageIndex){
			var nodes = this.getDefaultNodes();
			var last = 0;

			for(var i=0; i<this.stack.length; i++){
				last = Math.max(this.stack[i], last);
			}
			//
			var n = nodes[inPageIndex];
			var y = (n ? this.getPageNodePosition(n) + this.getPageHeight(inPageIndex) : 0);
			for(var p=inPageIndex+1; p<=last; p++){
				n = nodes[p];
				if(n){
					if(this.getPageNodePosition(n) == y){
						return;
					}
					this.positionPage(p, y);
				}
				y += this.getPageHeight(p);
			}
		},
		installPage: function(inPageIndex){
			for(var i=0; i<this.colCount; i++){
				this.contentNodes[i].appendChild(this.pageNodes[i][inPageIndex]);
			}
		},
		preparePage: function(inPageIndex, inReuseNode){
			var p = (inReuseNode ? this.popPage() : null);
			for(var i=0; i<this.colCount; i++){
				var nodes = this.pageNodes[i];
				var new_p = (p === null ? this.createPageNode() : this.invalidatePageNode(p, nodes));
				new_p.pageIndex = inPageIndex;
				nodes[inPageIndex] = new_p;
			}
		},
		// rendering implementation
		renderPage: function(inPageIndex){
			var nodes = [];
			var i, j;
			for(i=0; i<this.colCount; i++){
				nodes[i] = this.pageNodes[i][inPageIndex];
			}
			for(i=0, j=inPageIndex*this.rowsPerPage; (i<this.rowsPerPage)&&(j<this.rowCount); i++, j++){
				this.renderRow(j, nodes);
			}
		},
		removePage: function(inPageIndex){
			for(var i=0, j=inPageIndex*this.rowsPerPage; i<this.rowsPerPage; i++, j++){
				this.removeRow(j);
			}
		},
		destroyPage: function(inPageIndex){
			for(var i=0; i<this.colCount; i++){
				var n = this.invalidatePageNode(inPageIndex, this.pageNodes[i]);
				if(n){
					html.destroy(n);
				}
			}
		},
		pacify: function(inShouldPacify){
		},
		// pacification
		pacifying: false,
		pacifyTicks: 200,
		setPacifying: function(inPacifying){
			if(this.pacifying != inPacifying){
				this.pacifying = inPacifying;
				this.pacify(this.pacifying);
			}
		},
		startPacify: function(){
			this.startPacifyTicks = new Date().getTime();
		},
		doPacify: function(){
			var result = (new Date().getTime() - this.startPacifyTicks) > this.pacifyTicks;
			this.setPacifying(true);
			this.startPacify();
			return result;
		},
		endPacify: function(){
			this.setPacifying(false);
		},
		// default sizing implementation
		resize: function(){
			if(this.scrollboxNode){
				this.windowHeight = this.scrollboxNode.clientHeight;
			}
			for(var i=0; i<this.colCount; i++){
				//We want to have 1px in height min to keep scroller.  Otherwise can't scroll
				//and see header in empty grid.
				util.setStyleHeightPx(this.contentNodes[i], Math.max(1,this.height));
			}
			
			// Calculate the average row height and update the defaults (row and page).
			var needPage = (!this._invalidating);
			if(!needPage){
				var ah = this.grid.get("autoHeight");
				if(typeof ah == "number" && ah <= Math.min(this.rowsPerPage, this.rowCount)){
					needPage = true;
				}
			}
			if(needPage){
				this.needPage(this.page, this.pageTop);
			}
			var rowsOnPage = (this.page < this.pageCount - 1) ? this.rowsPerPage : ((this.rowCount % this.rowsPerPage) || this.rowsPerPage);
			var pageHeight = this.getPageHeight(this.page);
			this.averageRowHeight = (pageHeight > 0 && rowsOnPage > 0) ? (pageHeight / rowsOnPage) : 0;
		},
		calcLastPageHeight: function(){
			if(!this.pageCount){
				return 0;
			}
			var lastPage = this.pageCount - 1;
			var lastPageHeight = ((this.rowCount % this.rowsPerPage)||(this.rowsPerPage)) * this.defaultRowHeight;
			this.pageHeights[lastPage] = lastPageHeight;
			return lastPageHeight;
		},
		updateContentHeight: function(inDh){
			this.height += inDh;
			this.resize();
		},
		updatePageHeight: function(inPageIndex, fromBuild, fromAsynRendering){
			if(this.pageExists(inPageIndex)){
				var oh = this.getPageHeight(inPageIndex);
				var h = (this.measurePage(inPageIndex));
				if(h === undefined){
					h = oh;
				}
				this.pageHeights[inPageIndex] = h;
				if(oh != h){
					this.updateContentHeight(h - oh);
					var ah = this.grid.get("autoHeight");
					if((typeof ah == "number" && ah > this.rowCount)||(ah === true && !fromBuild)){
						if(!fromAsynRendering){
							this.grid.sizeChange();
						}else{//fix #11101 by using fromAsynRendering to avoid deadlock
							var ns = this.grid.viewsNode.style;
							ns.height = parseInt(ns.height) + h - oh + 'px';
							this.repositionPages(inPageIndex);
						}
					}else{
						this.repositionPages(inPageIndex);
					}
				}
				return h;
			}
			return 0;
		},
		rowHeightChanged: function(inRowIndex, fromAsynRendering){
			this.updatePageHeight(Math.floor(inRowIndex / this.rowsPerPage), false, fromAsynRendering);
		},
		// scroller core
		invalidateNodes: function(){
			while(this.stack.length){
				this.destroyPage(this.popPage());
			}
		},
		createPageNode: function(){
			var p = document.createElement('div');
			html.attr(p,"role","presentation");
			p.style.position = 'absolute';
			//p.style.width = '100%';
			p.style[this.grid.isLeftToRight() ? "left" : "right"] = '0';
			return p;
		},
		getPageHeight: function(inPageIndex){
			var ph = this.pageHeights[inPageIndex];
			return (ph !== undefined ? ph : this.defaultPageHeight);
		},
		// FIXME: this is not a stack, it's a FIFO list
		pushPage: function(inPageIndex){
			return this.stack.push(inPageIndex);
		},
		popPage: function(){
			return this.stack.shift();
		},
		findPage: function(inTop){
			var i = 0, h = 0;
			for(var ph = 0; i<this.pageCount; i++, h += ph){
				ph = this.getPageHeight(i);
				if(h + ph >= inTop){
					break;
				}
			}
			this.page = i;
			this.pageTop = h;
		},
		buildPage: function(inPageIndex, inReuseNode, inPos){
			this.preparePage(inPageIndex, inReuseNode);
			this.positionPage(inPageIndex, inPos);
			// order of operations is key below
			this.installPage(inPageIndex);
			this.renderPage(inPageIndex);
			// order of operations is key above
			this.pushPage(inPageIndex);
		},
		needPage: function(inPageIndex, inPos){
			var h = this.getPageHeight(inPageIndex), oh = h;
			if(!this.pageExists(inPageIndex)){
				this.buildPage(inPageIndex, (!this.grid._autoHeight/*fix #10543*/ && this.keepPages&&(this.stack.length >= this.keepPages)), inPos);
				h = this.updatePageHeight(inPageIndex, true);
			}else{
				this.positionPage(inPageIndex, inPos);
			}
			return h;
		},
		onscroll: function(){
			this.scroll(this.scrollboxNode.scrollTop);
		},
		scroll: function(inTop){
			this.grid.scrollTop = inTop;
			if(this.colCount){
				this.startPacify();
				this.findPage(inTop);
				var h = this.height;
				var b = this.getScrollBottom(inTop);
				for(var p=this.page, y=this.pageTop; (p<this.pageCount)&&((b<0)||(y<b)); p++){
					y += this.needPage(p, y);
				}
				this.firstVisibleRow = this.getFirstVisibleRow(this.page, this.pageTop, inTop);
				this.lastVisibleRow = this.getLastVisibleRow(p - 1, y, b);
				// indicates some page size has been updated
				if(h != this.height){
					this.repositionPages(p-1);
				}
				this.endPacify();
			}
		},
		getScrollBottom: function(inTop){
			return (this.windowHeight >= 0 ? inTop + this.windowHeight : -1);
		},
		// events
		processNodeEvent: function(e, inNode){
			var t = e.target;
			while(t && (t != inNode) && t.parentNode && (t.parentNode.parentNode != inNode)){
				t = t.parentNode;
			}
			if(!t || !t.parentNode || (t.parentNode.parentNode != inNode)){
				return false;
			}
			var page = t.parentNode;
			e.topRowIndex = page.pageIndex * this.rowsPerPage;
			e.rowIndex = e.topRowIndex + indexInParent(t);
			e.rowTarget = t;
			return true;
		},
		processEvent: function(e){
			return this.processNodeEvent(e, this.contentNode);
		},
		// virtual rendering interface
		renderRow: function(inRowIndex, inPageNode){
		},
		removeRow: function(inRowIndex){
		},
		// page node operations
		getDefaultPageNode: function(inPageIndex){
			return this.getDefaultNodes()[inPageIndex];
		},
		positionPageNode: function(inNode, inPos){
		},
		getPageNodePosition: function(inNode){
			return inNode.offsetTop;
		},
		invalidatePageNode: function(inPageIndex, inNodes){
			var p = inNodes[inPageIndex];
			if(p){
				delete inNodes[inPageIndex];
				this.removePage(inPageIndex, p);
				cleanNode(p);
				p.innerHTML = '';
			}
			return p;
		},
		// scroll control
		getPageRow: function(inPage){
			return inPage * this.rowsPerPage;
		},
		getLastPageRow: function(inPage){
			return Math.min(this.rowCount, this.getPageRow(inPage + 1)) - 1;
		},
		getFirstVisibleRow: function(inPage, inPageTop, inScrollTop){
			if(!this.pageExists(inPage)){
				return 0;
			}
			var row = this.getPageRow(inPage);
			var nodes = this.getDefaultNodes();
			var rows = divkids(nodes[inPage]);
			for(var i=0,l=rows.length; i<l && inPageTop<inScrollTop; i++, row++){
				inPageTop += rows[i].offsetHeight;
			}
			return (row ? row - 1 : row);
		},
		getLastVisibleRow: function(inPage, inBottom, inScrollBottom){
			if(!this.pageExists(inPage)){
				return 0;
			}
			var nodes = this.getDefaultNodes();
			var row = this.getLastPageRow(inPage);
			var rows = divkids(nodes[inPage]);
			for(var i=rows.length-1; i>=0 && inBottom>inScrollBottom; i--, row--){
				inBottom -= rows[i].offsetHeight;
			}
			return row + 1;
		},
		findTopRow: function(inScrollTop){
			var nodes = this.getDefaultNodes();
			var rows = divkids(nodes[this.page]);
			for(var i=0,l=rows.length,t=this.pageTop,h; i<l; i++){
				h = rows[i].offsetHeight;
				t += h;
				if(t >= inScrollTop){
					this.offset = h - (t - inScrollTop);
					return i + this.page * this.rowsPerPage;
				}
			}
			return -1;
		},
		findScrollTop: function(inRow){
			var rowPage = Math.floor(inRow / this.rowsPerPage);
			var t = 0;
			var i, l;
			for(i=0; i<rowPage; i++){
				t += this.getPageHeight(i);
			}
			this.pageTop = t;
			this.page = rowPage;//fix #10543
			this.needPage(rowPage, this.pageTop);

			var nodes = this.getDefaultNodes();
			var rows = divkids(nodes[rowPage]);
			var r = inRow - this.rowsPerPage * rowPage;
			for(i=0,l=rows.length; i<l && i<r; i++){
				t += rows[i].offsetHeight;
			}
			return t;
		},
		dummy: 0
	});
});

},
'dojox/grid/util':function(){
define("dojox/grid/util", [
	"../main",
	"dojo/_base/lang",
	"dojo/dom"
], function(dojox, lang, dom){

	var dgu = lang.getObject("grid.util", true, dojox);

/*=====
dgu = {
	// summary:
	//		grid utility library
};
=====*/

	dgu.na = '...';
	dgu.rowIndexTag = "gridRowIndex";
	dgu.gridViewTag = "gridView";


	dgu.fire = function(ob, ev, args){
		var fn = ob && ev && ob[ev];
		return fn && (args ? fn.apply(ob, args) : ob[ev]());
	};
	
	dgu.setStyleHeightPx = function(inElement, inHeight){
		if(inHeight >= 0){
			var s = inElement.style;
			var v = inHeight + 'px';
			if(inElement && s['height'] != v){
				s['height'] = v;
			}
		}
	};
	
	dgu.mouseEvents = [ 'mouseover', 'mouseout', /*'mousemove',*/ 'mousedown', 'mouseup', 'click', 'dblclick', 'contextmenu' ];

	dgu.keyEvents = [ 'keyup', 'keydown', 'keypress' ];

	dgu.funnelEvents = function(inNode, inObject, inMethod, inEvents){
		var evts = (inEvents ? inEvents : dgu.mouseEvents.concat(dgu.keyEvents));
		for (var i=0, l=evts.length; i<l; i++){
			inObject.connect(inNode, 'on' + evts[i], inMethod);
		}
	};

	dgu.removeNode = function(inNode){
		inNode = dom.byId(inNode);
		inNode && inNode.parentNode && inNode.parentNode.removeChild(inNode);
		return inNode;
	};
	
	dgu.arrayCompare = function(inA, inB){
		for(var i=0,l=inA.length; i<l; i++){
			if(inA[i] != inB[i]){return false;}
		}
		return (inA.length == inB.length);
	};
	
	dgu.arrayInsert = function(inArray, inIndex, inValue){
		if(inArray.length <= inIndex){
			inArray[inIndex] = inValue;
		}else{
			inArray.splice(inIndex, 0, inValue);
		}
	};
	
	dgu.arrayRemove = function(inArray, inIndex){
		inArray.splice(inIndex, 1);
	};
	
	dgu.arraySwap = function(inArray, inI, inJ){
		var cache = inArray[inI];
		inArray[inI] = inArray[inJ];
		inArray[inJ] = cache;
	};

	return dgu;

});
},
'dojox/grid/_Layout':function(){
define("dojox/grid/_Layout", [
	"dojo/_base/kernel",
	"../main",
	"dojo/_base/declare",
	"dojo/_base/array",
	"dojo/_base/lang",
	"dojo/dom-geometry",
	"./cells",
	"./_RowSelector"
], function(dojo, dojox, declare, array, lang, domGeometry){

return declare("dojox.grid._Layout", null, {
	// summary:
	//	Controls grid cell layout. Owned by grid and used internally.
	constructor: function(inGrid){
		this.grid = inGrid;
	},
	// flat array of grid cells
	cells: [],
	// structured array of grid cells
	structure: null,
	// default cell width
	defaultWidth: '6em',

	// methods
	moveColumn: function(sourceViewIndex, destViewIndex, cellIndex, targetIndex, before){
		var source_cells = this.structure[sourceViewIndex].cells[0];
		var dest_cells = this.structure[destViewIndex].cells[0];

		var cell = null;
		var cell_ri = 0;
		var target_ri = 0;

		for(var i=0, c; c=source_cells[i]; i++){
			if(c.index == cellIndex){
				cell_ri = i;
				break;
			}
		}
		cell = source_cells.splice(cell_ri, 1)[0];
		cell.view = this.grid.views.views[destViewIndex];

		for(i=0, c=null; c=dest_cells[i]; i++){
			if(c.index == targetIndex){
				target_ri = i;
				break;
			}
		}
		if(!before){
			target_ri += 1;
		}
		dest_cells.splice(target_ri, 0, cell);

		var sortedCell = this.grid.getCell(this.grid.getSortIndex());
		if(sortedCell){
			sortedCell._currentlySorted = this.grid.getSortAsc();
		}

		this.cells = [];
		cellIndex = 0;
		var v;
		for(i=0; v=this.structure[i]; i++){
			for(var j=0, cs; cs=v.cells[j]; j++){
				for(var k=0; c=cs[k]; k++){
					c.index = cellIndex;
					this.cells.push(c);
					if("_currentlySorted" in c){
						var si = cellIndex + 1;
						si *= c._currentlySorted ? 1 : -1;
						this.grid.sortInfo = si;
						delete c._currentlySorted;
					}
					cellIndex++;
				}
			}
		}
		
		//Fix #9481 - reset idx in cell markup
		array.forEach(this.cells, function(c){
			var marks = c.markup[2].split(" ");
			var oldIdx = parseInt(marks[1].substring(5));//get old "idx"
			if(oldIdx != c.index){
				marks[1] = "idx=\"" + c.index + "\"";
				c.markup[2] = marks.join(" ");
			}
		});
		
		this.grid.setupHeaderMenu();
		//this.grid.renderOnIdle();
	},

	setColumnVisibility: function(columnIndex, visible){
		var cell = this.cells[columnIndex];
		if(cell.hidden == visible){
			cell.hidden = !visible;
			var v = cell.view, w = v.viewWidth;
			if(w && w != "auto"){
				v._togglingColumn = domGeometry.getMarginBox(cell.getHeaderNode()).w || 0;
			}
			v.update();
			return true;
		}else{
			return false;
		}
	},
	
	addCellDef: function(inRowIndex, inCellIndex, inDef){
		var self = this;
		var getCellWidth = function(inDef){
			var w = 0;
			if(inDef.colSpan > 1){
				w = 0;
			}else{
				w = inDef.width || self._defaultCellProps.width || self.defaultWidth;

				if(!isNaN(w)){
					w = w + "em";
				}
			}
			return w;
		};

		var props = {
			grid: this.grid,
			subrow: inRowIndex,
			layoutIndex: inCellIndex,
			index: this.cells.length
		};

		if(inDef && inDef instanceof dojox.grid.cells._Base){
			var new_cell = lang.clone(inDef);
			props.unitWidth = getCellWidth(new_cell._props);
			new_cell = lang.mixin(new_cell, this._defaultCellProps, inDef._props, props);
			return new_cell;
		}

		var cell_type = inDef.type || inDef.cellType || this._defaultCellProps.type || this._defaultCellProps.cellType || dojox.grid.cells.Cell;
		if(lang.isString(cell_type)){
			cell_type = lang.getObject(cell_type);
		}

		props.unitWidth = getCellWidth(inDef);
		return new cell_type(lang.mixin({}, this._defaultCellProps, inDef, props));
	},
	
	addRowDef: function(inRowIndex, inDef){
		var result = [];
		var relSum = 0, pctSum = 0, doRel = true;
		for(var i=0, def, cell; (def=inDef[i]); i++){
			cell = this.addCellDef(inRowIndex, i, def);
			result.push(cell);
			this.cells.push(cell);
			// Check and calculate the sum of all relative widths
			if(doRel && cell.relWidth){
				relSum += cell.relWidth;
			}else if(cell.width){
				var w = cell.width;
				if(typeof w == "string" && w.slice(-1) == "%"){
					pctSum += window.parseInt(w, 10);
				}else if(w == "auto"){
					// relative widths doesn't play nice with auto - since we
					// don't have a way of knowing how much space the auto is
					// supposed to take up.
					doRel = false;
				}
			}
		}
		if(relSum && doRel){
			// We have some kind of relWidths specified - so change them to %
			array.forEach(result, function(cell){
				if(cell.relWidth){
					cell.width = cell.unitWidth = ((cell.relWidth / relSum) * (100 - pctSum)) + "%";
				}
			});
		}
		return result;
	
	},

	addRowsDef: function(inDef){
		var result = [];
		if(lang.isArray(inDef)){
			if(lang.isArray(inDef[0])){
				for(var i=0, row; inDef && (row=inDef[i]); i++){
					result.push(this.addRowDef(i, row));
				}
			}else{
				result.push(this.addRowDef(0, inDef));
			}
		}
		return result;
	},
	
	addViewDef: function(inDef){
		this._defaultCellProps = inDef.defaultCell || {};
		if(inDef.width && inDef.width == "auto"){
			delete inDef.width;
		}
		return lang.mixin({}, inDef, {cells: this.addRowsDef(inDef.rows || inDef.cells)});
	},
	
	setStructure: function(inStructure){
		this.fieldIndex = 0;
		this.cells = [];
		var s = this.structure = [];

		if(this.grid.rowSelector){
			var sel = { type: dojox._scopeName + ".grid._RowSelector" };

			if(lang.isString(this.grid.rowSelector)){
				var width = this.grid.rowSelector;

				if(width == "false"){
					sel = null;
				}else if(width != "true"){
					sel['width'] = width;
				}
			}else{
				if(!this.grid.rowSelector){
					sel = null;
				}
			}

			if(sel){
				s.push(this.addViewDef(sel));
			}
		}

		var isCell = function(def){
			return ("name" in def || "field" in def || "get" in def);
		};

		var isRowDef = function(def){
			if(lang.isArray(def)){
				if(lang.isArray(def[0]) || isCell(def[0])){
					return true;
				}
			}
			return false;
		};

		var isView = function(def){
			return (def !== null && lang.isObject(def) &&
					("cells" in def || "rows" in def || ("type" in def && !isCell(def))));
		};

		if(lang.isArray(inStructure)){
			var hasViews = false;
			for(var i=0, st; (st=inStructure[i]); i++){
				if(isView(st)){
					hasViews = true;
					break;
				}
			}
			if(!hasViews){
				s.push(this.addViewDef({ cells: inStructure }));
			}else{
				for(i=0; (st=inStructure[i]); i++){
					if(isRowDef(st)){
						s.push(this.addViewDef({ cells: st }));
					}else if(isView(st)){
						s.push(this.addViewDef(st));
					}
				}
			}
		}else if(isView(inStructure)){
			// it's a view object
			s.push(this.addViewDef(inStructure));
		}

		this.cellCount = this.cells.length;
		this.grid.setupHeaderMenu();
	}
});
});
},
'dojox/grid/cells':function(){
define("dojox/grid/cells", ["../main", "./cells/_base"], function(dojox){
	return dojox.grid.cells;
});
},
'dojox/grid/cells/_base':function(){
define("dojox/grid/cells/_base", [
	"dojo/_base/kernel",
	"dojo/_base/declare",
	"dojo/_base/lang",
	"dojo/_base/event",
	"dojo/_base/connect",
	"dojo/_base/array",
	"dojo/_base/sniff",
	"dojo/dom",
	"dojo/dom-attr",
	"dojo/dom-construct",
	"dijit/_Widget",
	"../util"
], function(dojo, declare, lang, event, connect, array, has, dom, domAttr, domConstruct, _Widget, util){

	var _DeferredTextWidget = declare("dojox.grid._DeferredTextWidget", _Widget, {
		deferred: null,
		_destroyOnRemove: true,
		postCreate: function(){
			if(this.deferred){
				this.deferred.addBoth(lang.hitch(this, function(text){
					if(this.domNode){
						this.domNode.innerHTML = text;
					}
				}));
			}
		}
	});

	var focusSelectNode = function(inNode){
		try{
			util.fire(inNode, "focus");
			util.fire(inNode, "select");
		}catch(e){// IE sux bad
		}
	};
	
	var whenIdle = function(/*inContext, inMethod, args ...*/){
		setTimeout(lang.hitch.apply(dojo, arguments), 0);
	};

	var BaseCell = declare("dojox.grid.cells._Base", null, {
		// summary:
		//		Represents a grid cell and contains information about column options and methods
		//		for retrieving cell related information.
		//		Each column in a grid layout has a cell object and most events and many methods
		//		provide access to these objects.
		styles: '',
		classes: '',
		editable: false,
		alwaysEditing: false,
		formatter: null,
		defaultValue: '...',
		value: null,
		hidden: false,
		noresize: false,
		draggable: true,
		//private
		_valueProp: "value",
		_formatPending: false,

		constructor: function(inProps){
			this._props = inProps || {};
			lang.mixin(this, inProps);
			if(this.draggable === undefined){
				this.draggable = true;
			}
		},

		_defaultFormat: function(inValue, callArgs){
			var s = this.grid.formatterScope || this;
			var f = this.formatter;
			if(f && s && typeof f == "string"){
				f = this.formatter = s[f];
			}
			var v = (inValue != this.defaultValue && f) ? f.apply(s, callArgs) : inValue;
			if(typeof v == "undefined"){
				return this.defaultValue;
			}
			if(v && v.addBoth){
				// Check if it's a deferred
				v = new _DeferredTextWidget({deferred: v},
									domConstruct.create("span", {innerHTML: this.defaultValue}));
			}
			if(v && v.declaredClass && v.startup){
				return "<div class='dojoxGridStubNode' linkWidget='" +
						v.id +
						"' cellIdx='" +
						this.index +
						"'>" +
						this.defaultValue +
						"</div>";
			}
			return v;
		},
		
		// data source
		format: function(inRowIndex, inItem){
			// summary:
			//		provides the html for a given grid cell.
			// inRowIndex: int
			//		grid row index
			// returns:
			//		html for a given grid cell
			var f, i=this.grid.edit.info, d=this.get ? this.get(inRowIndex, inItem) : (this.value || this.defaultValue);
			d = (d && d.replace && this.grid.escapeHTMLInData) ? d.replace(/&/g, '&amp;').replace(/</g, '&lt;') : d;
			if(this.editable && (this.alwaysEditing || (i.rowIndex==inRowIndex && i.cell==this))){
				return this.formatEditing(d, inRowIndex);
			}else{
				return this._defaultFormat(d, [d, inRowIndex, this]);
			}
		},
		formatEditing: function(inDatum, inRowIndex){
			// summary:
			//		formats the cell for editing
			// inDatum: anything
			//		cell data to edit
			// inRowIndex: int
			//		grid row index
			// returns:
			//		string of html to place in grid cell
		},
		// utility
		getNode: function(inRowIndex){
			// summary:
			//		gets the dom node for a given grid cell.
			// inRowIndex: int
			//		grid row index
			// returns:
			//		dom node for a given grid cell
			return this.view.getCellNode(inRowIndex, this.index);
		},
		getHeaderNode: function(){
			return this.view.getHeaderCellNode(this.index);
		},
		getEditNode: function(inRowIndex){
			return (this.getNode(inRowIndex) || 0).firstChild || 0;
		},
		canResize: function(){
			var uw = this.unitWidth;
			return uw && (uw!=='auto');
		},
		isFlex: function(){
			var uw = this.unitWidth;
			return uw && lang.isString(uw) && (uw=='auto' || uw.slice(-1)=='%');
		},
		// edit support
		applyEdit: function(inValue, inRowIndex){
			if(this.getNode(inRowIndex)){
				this.grid.edit.applyCellEdit(inValue, this, inRowIndex);
			}
		},
		cancelEdit: function(inRowIndex){
			this.grid.doCancelEdit(inRowIndex);
		},
		_onEditBlur: function(inRowIndex){
			if(this.grid.edit.isEditCell(inRowIndex, this.index)){
				//console.log('editor onblur', e);
				this.grid.edit.apply();
			}
		},
		registerOnBlur: function(inNode, inRowIndex){
			if(this.commitOnBlur){
				connect.connect(inNode, "onblur", function(e){
					// hack: if editor still thinks this editor is current some ms after it blurs, assume we've focused away from grid
					setTimeout(lang.hitch(this, "_onEditBlur", inRowIndex), 250);
				});
			}
		},
		//protected
		needFormatNode: function(inDatum, inRowIndex){
			this._formatPending = true;
			whenIdle(this, "_formatNode", inDatum, inRowIndex);
		},
		cancelFormatNode: function(){
			this._formatPending = false;
		},
		//private
		_formatNode: function(inDatum, inRowIndex){
			if(this._formatPending){
				this._formatPending = false;
				// make cell selectable
				if(!has('ie')){
					dom.setSelectable(this.grid.domNode, true);
				}
				this.formatNode(this.getEditNode(inRowIndex), inDatum, inRowIndex);
			}
		},
		//protected
		formatNode: function(inNode, inDatum, inRowIndex){
			// summary:
			//		format the editing dom node. Use when editor is a widget.
			// inNode: dom node
			//		dom node for the editor
			// inDatum: anything
			//		cell data to edit
			// inRowIndex: int
			//		grid row index
			if(has('ie')){
				// IE sux bad
				whenIdle(this, "focus", inRowIndex, inNode);
			}else{
				this.focus(inRowIndex, inNode);
			}
		},
		dispatchEvent: function(m, e){
			if(m in this){
				return this[m](e);
			}
		},
		//public
		getValue: function(inRowIndex){
			// summary:
			//		returns value entered into editor
			// inRowIndex: int
			//		grid row index
			// returns:
			//		value of editor
			return this.getEditNode(inRowIndex)[this._valueProp];
		},
		setValue: function(inRowIndex, inValue){
			// summary:
			//		set the value of the grid editor
			// inRowIndex: int
			//		grid row index
			// inValue: anything
			//		value of editor
			var n = this.getEditNode(inRowIndex);
			if(n){
				n[this._valueProp] = inValue;
			}
		},
		focus: function(inRowIndex, inNode){
			// summary:
			//		focus the grid editor
			// inRowIndex: int
			//		grid row index
			// inNode: dom node
			//		editor node
			focusSelectNode(inNode || this.getEditNode(inRowIndex));
		},
		save: function(inRowIndex){
			// summary:
			//		save editor state
			// inRowIndex: int
			//		grid row index
			this.value = this.value || this.getValue(inRowIndex);
			//console.log("save", this.value, inCell.index, inRowIndex);
		},
		restore: function(inRowIndex){
			// summary:
			//		restore editor state
			// inRowIndex: int
			//		grid row index
			this.setValue(inRowIndex, this.value);
			//console.log("restore", this.value, inCell.index, inRowIndex);
		},
		//protected
		_finish: function(inRowIndex){
			// summary:
			//		called when editing is completed to clean up editor
			// inRowIndex: int
			//		grid row index
			dom.setSelectable(this.grid.domNode, false);
			this.cancelFormatNode();
		},
		//public
		apply: function(inRowIndex){
			// summary:
			//		apply edit from cell editor
			// inRowIndex: int
			//		grid row index
			this.applyEdit(this.getValue(inRowIndex), inRowIndex);
			this._finish(inRowIndex);
		},
		cancel: function(inRowIndex){
			// summary:
			//		cancel cell edit
			// inRowIndex: int
			//		grid row index
			this.cancelEdit(inRowIndex);
			this._finish(inRowIndex);
		}
	});
	BaseCell.markupFactory = function(node, cellDef){
		var formatter = lang.trim(domAttr.get(node, "formatter")||"");
		if(formatter){
			cellDef.formatter = lang.getObject(formatter)||formatter;
		}
		var get = lang.trim(domAttr.get(node, "get")||"");
		if(get){
			cellDef.get = lang.getObject(get);
		}
		var getBoolAttr = function(attr, cell, cellAttr){
			var value = lang.trim(domAttr.get(node, attr)||"");
			if(value){ cell[cellAttr||attr] = !(value.toLowerCase()=="false"); }
		};
		getBoolAttr("sortDesc", cellDef);
		getBoolAttr("editable", cellDef);
		getBoolAttr("alwaysEditing", cellDef);
		getBoolAttr("noresize", cellDef);
		getBoolAttr("draggable", cellDef);

		var value = lang.trim(domAttr.get(node, "loadingText")||domAttr.get(node, "defaultValue")||"");
		if(value){
			cellDef.defaultValue = value;
		}

		var getStrAttr = function(attr, cell, cellAttr){
			var value = lang.trim(domAttr.get(node, attr)||"")||undefined;
			if(value){ cell[cellAttr||attr] = value; }
		};
		getStrAttr("styles", cellDef);
		getStrAttr("headerStyles", cellDef);
		getStrAttr("cellStyles", cellDef);
		getStrAttr("classes", cellDef);
		getStrAttr("headerClasses", cellDef);
		getStrAttr("cellClasses", cellDef);
	};

	var Cell = declare("dojox.grid.cells.Cell", BaseCell, {
		// summary:
		//		grid cell that provides a standard text input box upon editing
		constructor: function(){
			this.keyFilter = this.keyFilter;
		},
		// keyFilter: RegExp
		//		optional regex for disallowing keypresses
		keyFilter: null,
		formatEditing: function(inDatum, inRowIndex){
			this.needFormatNode(inDatum, inRowIndex);
			return '<input class="dojoxGridInput" type="text" value="' + inDatum + '">';
		},
		formatNode: function(inNode, inDatum, inRowIndex){
			this.inherited(arguments);
			// FIXME: feels too specific for this interface
			this.registerOnBlur(inNode, inRowIndex);
		},
		doKey: function(e){
			if(this.keyFilter){
				var key = String.fromCharCode(e.charCode);
				if(key.search(this.keyFilter) == -1){
					event.stop(e);
				}
			}
		},
		_finish: function(inRowIndex){
			this.inherited(arguments);
			var n = this.getEditNode(inRowIndex);
			try{
				util.fire(n, "blur");
			}catch(e){}
		}
	});
	Cell.markupFactory = function(node, cellDef){
		BaseCell.markupFactory(node, cellDef);
		var keyFilter = lang.trim(domAttr.get(node, "keyFilter")||"");
		if(keyFilter){
			cellDef.keyFilter = new RegExp(keyFilter);
		}
	};

	var RowIndex = declare("dojox.grid.cells.RowIndex", Cell, {
		name: 'Row',

		postscript: function(){
			this.editable = false;
		},
		get: function(inRowIndex){
			return inRowIndex + 1;
		}
	});
	RowIndex.markupFactory = function(node, cellDef){
		Cell.markupFactory(node, cellDef);
	};

	var Select = declare("dojox.grid.cells.Select", Cell, {
		// summary:
		//		grid cell that provides a standard select for editing

		// options: Array
		//		text of each item
		options: null,

		// values: Array
		//		value for each item
		values: null,

		// returnIndex: Integer
		//		editor returns only the index of the selected option and not the value
		returnIndex: -1,

		constructor: function(inCell){
			this.values = this.values || this.options;
		},
		formatEditing: function(inDatum, inRowIndex){
			this.needFormatNode(inDatum, inRowIndex);
			var h = [ '<select class="dojoxGridSelect">' ];
			for (var i=0, o, v; ((o=this.options[i]) !== undefined)&&((v=this.values[i]) !== undefined); i++){
				v = v.replace ? v.replace(/&/g, '&amp;').replace(/</g, '&lt;') : v;
				o = o.replace ? o.replace(/&/g, '&amp;').replace(/</g, '&lt;') : o;
				h.push("<option", (inDatum==v ? ' selected' : ''), ' value="' + v + '"', ">", o, "</option>");
			}
			h.push('</select>');
			return h.join('');
		},
		_defaultFormat: function(inValue, callArgs){
			var v = this.inherited(arguments);
			// when 'values' and 'options' both provided and there is no cutomized formatter,
			// then we use 'options' as label in order to be consistent
			if(!this.formatter && this.values && this.options){
				var i = array.indexOf(this.values, v);
				if(i >= 0){
					v = this.options[i];
				}
			}
			return v;
		},
		getValue: function(inRowIndex){
			var n = this.getEditNode(inRowIndex);
			if(n){
				var i = n.selectedIndex, o = n.options[i];
				return this.returnIndex > -1 ? i : o.value || o.innerHTML;
			}
		}
	});
	Select.markupFactory = function(node, cell){
		Cell.markupFactory(node, cell);
		var options = lang.trim(domAttr.get(node, "options")||"");
		if(options){
			var o = options.split(',');
			if(o[0] != options){
				cell.options = o;
			}
		}
		var values = lang.trim(domAttr.get(node, "values")||"");
		if(values){
			var v = values.split(',');
			if(v[0] != values){
				cell.values = v;
			}
		}
	};

	var AlwaysEdit = declare("dojox.grid.cells.AlwaysEdit", Cell, {
		// summary:
		//		grid cell that is always in an editable state, regardless of grid editing state
		alwaysEditing: true,
		_formatNode: function(inDatum, inRowIndex){
			this.formatNode(this.getEditNode(inRowIndex), inDatum, inRowIndex);
		},
		applyStaticValue: function(inRowIndex){
			var e = this.grid.edit;
			e.applyCellEdit(this.getValue(inRowIndex), this, inRowIndex);
			e.start(this, inRowIndex, true);
		}
	});
	AlwaysEdit.markupFactory = function(node, cell){
		Cell.markupFactory(node, cell);
	};

	var Bool = declare("dojox.grid.cells.Bool", AlwaysEdit, {
		// summary:
		//		grid cell that provides a standard checkbox that is always on for editing
		_valueProp: "checked",
		formatEditing: function(inDatum, inRowIndex){
			return '<input class="dojoxGridInput" type="checkbox"' + (inDatum ? ' checked="checked"' : '') + ' style="width: auto" />';
		},
		doclick: function(e){
			if(e.target.tagName == 'INPUT'){
				this.applyStaticValue(e.rowIndex);
			}
		}
	});
	Bool.markupFactory = function(node, cell){
		AlwaysEdit.markupFactory(node, cell);
	};

	return BaseCell;

});
},
'dojox/grid/_RowSelector':function(){
define("dojox/grid/_RowSelector", [
	"dojo/_base/declare",
	"./_View"
], function(declare, _View){

return declare('dojox.grid._RowSelector', _View, {
	// summary:
	//	Custom grid view. If used in a grid structure, provides a small selectable region for grid rows.
	defaultWidth: "2em",
	noscroll: true,
	padBorderWidth: 2,
	buildRendering: function(){
		this.inherited('buildRendering', arguments);
		this.scrollboxNode.style.overflow = "hidden";
		this.headerNode.style.visibility = "hidden";
	},
	getWidth: function(){
		return this.viewWidth || this.defaultWidth;
	},
	buildRowContent: function(inRowIndex, inRowNode){
		var w = this.contentWidth || 0;
		inRowNode.innerHTML = '<table class="dojoxGridRowbarTable" style="width:' + w + 'px;height:1px;" border="0" cellspacing="0" cellpadding="0" role="presentation"><tr><td class="dojoxGridRowbarInner">&nbsp;</td></tr></table>';
	},
	renderHeader: function(){
	},
	updateRow: function(){
	},
	resize: function(){
		this.adaptHeight();
	},
	adaptWidth: function(){
		// Only calculate this here - rather than every call to buildRowContent
		if(!("contentWidth" in this) && this.contentNode && this.contentNode.offsetWidth > 0){
			this.contentWidth = this.contentNode.offsetWidth - this.padBorderWidth;
		}
	},
	// styling
	doStyleRowNode: function(inRowIndex, inRowNode){
		var n = [ "dojoxGridRowbar dojoxGridNonNormalizedCell" ];
		if(this.grid.rows.isOver(inRowIndex)){
			n.push("dojoxGridRowbarOver");
		}
		if(this.grid.selection.isSelected(inRowIndex)){
			n.push("dojoxGridRowbarSelected");
		}
		inRowNode.className = n.join(" ");
	},
	// event handlers
	domouseover: function(e){
		this.grid.onMouseOverRow(e);
	},
	domouseout: function(e){
		if(!this.isIntraRowEvent(e)){
			this.grid.onMouseOutRow(e);
		}
	}
});
});
},
'dojox/grid/_View':function(){
require({cache:{
'url:dojox/grid/resources/View.html':"<div class=\"dojoxGridView\" role=\"presentation\">\n\t<div class=\"dojoxGridHeader\" dojoAttachPoint=\"headerNode\" role=\"presentation\">\n\t\t<div dojoAttachPoint=\"headerNodeContainer\" style=\"width:9000em\" role=\"presentation\">\n\t\t\t<div dojoAttachPoint=\"headerContentNode\" role=\"row\"></div>\n\t\t</div>\n\t</div>\n\t<input type=\"checkbox\" class=\"dojoxGridHiddenFocus\" dojoAttachPoint=\"hiddenFocusNode\" role=\"presentation\" />\n\t<input type=\"checkbox\" class=\"dojoxGridHiddenFocus\" role=\"presentation\" />\n\t<div class=\"dojoxGridScrollbox\" dojoAttachPoint=\"scrollboxNode\" role=\"presentation\">\n\t\t<div class=\"dojoxGridContent\" dojoAttachPoint=\"contentNode\" hidefocus=\"hidefocus\" role=\"presentation\"></div>\n\t</div>\n</div>\n"}});
define("dojox/grid/_View", [
	"dojo",
	"dijit/registry",
	"../main",
	"dojo/_base/declare",
	"dojo/_base/array",
	"dojo/_base/lang",
	"dojo/_base/connect",
	"dojo/_base/sniff",
	"dojo/query",
	"dojo/_base/window",
	"dojo/text!./resources/View.html",
	"dojo/dnd/Source",
	"dijit/_Widget",
	"dijit/_TemplatedMixin",
	"dojox/html/metrics",
	"./util",
	"dojo/_base/html",
	"./_Builder",
	"dojo/dnd/Avatar",
	"dojo/dnd/Manager"
], function(dojo, dijit, dojox, declare, array, lang, connect, has, query,
	win, template, Source, _Widget, _TemplatedMixin, metrics, util, html, _Builder, Avatar, Manager){

	// a private function
	var getStyleText = function(inNode, inStyleText){
		return inNode.style.cssText == undefined ? inNode.getAttribute("style") : inNode.style.cssText;
	};

	// some public functions
	var _View = declare('dojox.grid._View', [_Widget, _TemplatedMixin], {
		// summary:
		//		A collection of grid columns. A grid is comprised of a set of views that stack horizontally.
		//		Grid creates views automatically based on grid's layout structure.
		//		Users should typically not need to access individual views directly.
		//
		// defaultWidth: String
		//		Default width of the view
		defaultWidth: "18em",

		// viewWidth: String
		//		Width for the view, in valid css unit
		viewWidth: "",

		templateString: template,

		classTag: 'dojoxGrid',
		marginBottom: 0,
		rowPad: 2,

		// _togglingColumn: int
		//		Width of the column being toggled (-1 for none)
		_togglingColumn: -1,
		
		// _headerBuilderClass: Object
		//		The class to use for our header builder
		_headerBuilderClass: _Builder._HeaderBuilder,
		
		// _contentBuilderClass: Object
		//		The class to use for our content builder
		_contentBuilderClass: _Builder._ContentBuilder,
		
		postMixInProperties: function(){
			this.rowNodes = {};
		},

		postCreate: function(){
			this.connect(this.scrollboxNode,"onscroll","doscroll");
			util.funnelEvents(this.contentNode, this, "doContentEvent", [ 'mouseover', 'mouseout', 'click', 'dblclick', 'contextmenu', 'mousedown' ]);
			util.funnelEvents(this.headerNode, this, "doHeaderEvent", [ 'dblclick', 'mouseover', 'mouseout', 'mousemove', 'mousedown', 'click', 'contextmenu' ]);
			this.content = new this._contentBuilderClass(this);
			this.header = new this._headerBuilderClass(this);
			//BiDi: in RTL case, style width='9000em' causes scrolling problem in head node
			if(!this.grid.isLeftToRight()){
				this.headerNodeContainer.style.width = "";
			}
		},

		destroy: function(){
			html.destroy(this.headerNode);
			delete this.headerNode;
			for(var i in this.rowNodes){
				this._cleanupRowWidgets(this.rowNodes[i]);
				html.destroy(this.rowNodes[i]);
			}
			this.rowNodes = {};
			if(this.source){
				this.source.destroy();
			}
			this.inherited(arguments);
		},

		// focus
		focus: function(){
			if(has('ie') || has('webkit') || has('opera')){
				this.hiddenFocusNode.focus();
			}else{
				this.scrollboxNode.focus();
			}
		},

		setStructure: function(inStructure){
			var vs = (this.structure = inStructure);
			// FIXME: similar logic is duplicated in layout
			if(vs.width && !isNaN(vs.width)){
				this.viewWidth = vs.width + 'em';
			}else{
				this.viewWidth = vs.width || (vs.noscroll ? 'auto' : this.viewWidth); //|| this.defaultWidth;
			}
			this._onBeforeRow = vs.onBeforeRow||function(){};
			this._onAfterRow = vs.onAfterRow||function(){};
			this.noscroll = vs.noscroll;
			if(this.noscroll){
				this.scrollboxNode.style.overflow = "hidden";
			}
			this.simpleStructure = Boolean(vs.cells.length == 1);
			// bookkeeping
			this.testFlexCells();
			// accomodate new structure
			this.updateStructure();
		},
		
		_cleanupRowWidgets: function(inRowNode){
			// Summary:
			//		Cleans up the widgets for the given row node so that
			//		we can reattach them if needed
			if(inRowNode){
				array.forEach(query("[widgetId]", inRowNode).map(dijit.byNode), function(w){
					if(w._destroyOnRemove){
						w.destroy();
						delete w;
					}else if(w.domNode && w.domNode.parentNode){
						w.domNode.parentNode.removeChild(w.domNode);
					}
				});
			}
		},
		
		onBeforeRow: function(inRowIndex, cells){
			this._onBeforeRow(inRowIndex, cells);
			if(inRowIndex >= 0){
				this._cleanupRowWidgets(this.getRowNode(inRowIndex));
			}
		},
		
		onAfterRow: function(inRowIndex, cells, inRowNode){
			this._onAfterRow(inRowIndex, cells, inRowNode);
			var g = this.grid;
			array.forEach(query(".dojoxGridStubNode", inRowNode), function(n){
				if(n && n.parentNode){
					var lw = n.getAttribute("linkWidget");
					var cellIdx = window.parseInt(html.attr(n, "cellIdx"), 10);
					var cellDef = g.getCell(cellIdx);
					var w = dijit.byId(lw);
					if(w){
						n.parentNode.replaceChild(w.domNode, n);
						if(!w._started){
							w.startup();
						}
						dojo.destroy(n);
					}else{
						n.innerHTML = "";
					}
				}
			}, this);
		},

		testFlexCells: function(){
			// FIXME: cheater, this function does double duty as initializer and tester
			this.flexCells = false;
			for(var j=0, row; (row=this.structure.cells[j]); j++){
				for(var i=0, cell; (cell=row[i]); i++){
					cell.view = this;
					this.flexCells = this.flexCells || cell.isFlex();
				}
			}
			return this.flexCells;
		},

		updateStructure: function(){
			// header builder needs to update table map
			this.header.update();
			// content builder needs to update markup cache
			this.content.update();
		},

		getScrollbarWidth: function(){
			var hasScrollSpace = this.hasVScrollbar();
			var overflow = html.style(this.scrollboxNode, "overflow");
			if(this.noscroll || !overflow || overflow == "hidden"){
				hasScrollSpace = false;
			}else if(overflow == "scroll"){
				hasScrollSpace = true;
			}
			return (hasScrollSpace ? metrics.getScrollbar().w : 0); // Integer
		},

		getColumnsWidth: function(){
			var h = this.headerContentNode;
			return h && h.firstChild ? (h.firstChild.offsetWidth || html.style(h.firstChild, 'width')) : 0; // Integer
		},

		setColumnsWidth: function(width){
			this.headerContentNode.firstChild.style.width = width + 'px';
			if(this.viewWidth){
				this.viewWidth = width + 'px';
			}
		},

		getWidth: function(){
			return this.viewWidth || (this.getColumnsWidth()+this.getScrollbarWidth()) +'px'; // String
		},

		getContentWidth: function(){
			return Math.max(0, html._getContentBox(this.domNode).w - this.getScrollbarWidth()) + 'px'; // String
		},

		render: function(){
			this.scrollboxNode.style.height = '';
			this.renderHeader();
			if(this._togglingColumn >= 0){
				this.setColumnsWidth(this.getColumnsWidth() - this._togglingColumn);
				this._togglingColumn = -1;
			}
			var cells = this.grid.layout.cells;
			var getSibling = lang.hitch(this, function(node, before){
				!this.grid.isLeftToRight() && (before = !before);
				var inc = before?-1:1;
				var idx = this.header.getCellNodeIndex(node) + inc;
				var cell = cells[idx];
				while(cell && cell.getHeaderNode() && cell.getHeaderNode().style.display == "none"){
					idx += inc;
					cell = cells[idx];
				}
				if(cell){
					return cell.getHeaderNode();
				}
				return null;
			});
			if(this.grid.columnReordering && this.simpleStructure){
				if(this.source){
					this.source.destroy();
				}
				
				// Create the top and bottom markers
				var bottomMarkerId = "dojoxGrid_bottomMarker";
				var topMarkerId = "dojoxGrid_topMarker";
				if(this.bottomMarker){
					html.destroy(this.bottomMarker);
				}
				this.bottomMarker = html.byId(bottomMarkerId);
				if(this.topMarker){
					html.destroy(this.topMarker);
				}
				this.topMarker = html.byId(topMarkerId);
				if (!this.bottomMarker) {
					this.bottomMarker = html.create("div", {
						"id": bottomMarkerId,
						"class": "dojoxGridColPlaceBottom"
					}, win.body());
					this._hide(this.bottomMarker);

					
					this.topMarker = html.create("div", {
						"id": topMarkerId,
						"class": "dojoxGridColPlaceTop"
					}, win.body());
					this._hide(this.topMarker);
				}
				this.arrowDim = html.contentBox(this.bottomMarker);

				var headerHeight = html.contentBox(this.headerContentNode.firstChild.rows[0]).h;
				
				this.source = new Source(this.headerContentNode.firstChild.rows[0], {
					horizontal: true,
					accept: [ "gridColumn_" + this.grid.id ],
					viewIndex: this.index,
					generateText: false,
					onMouseDown: lang.hitch(this, function(e){
						this.header.decorateEvent(e);
						if((this.header.overRightResizeArea(e) || this.header.overLeftResizeArea(e)) &&
							this.header.canResize(e) && !this.header.moveable){
							this.header.beginColumnResize(e);
						}else{
							if(this.grid.headerMenu){
								this.grid.headerMenu.onCancel(true);
							}
							// IE reports a left click as 1, where everything else reports 0
							if(e.button === (has('ie') < 9 ? 1 : 0)){
								Source.prototype.onMouseDown.call(this.source, e);
							}
						}
					}),
					onMouseOver: lang.hitch(this, function(e){
						var src = this.source;
						if(src._getChildByEvent(e)){
							Source.prototype.onMouseOver.apply(src, arguments);
						}
					}),
					_markTargetAnchor: lang.hitch(this, function(before){
						var src = this.source;
						if(src.current == src.targetAnchor && src.before == before){ return; }
						if(src.targetAnchor && getSibling(src.targetAnchor, src.before)){
							src._removeItemClass(getSibling(src.targetAnchor, src.before), src.before ? "After" : "Before");
						}
						Source.prototype._markTargetAnchor.call(src, before);
						
						var target = before ? src.targetAnchor : getSibling(src.targetAnchor, src.before);
						var endAdd = 0;

						if (!target) {
							target = src.targetAnchor;
							endAdd = html.contentBox(target).w + this.arrowDim.w/2 + 2;
						}

						var pos = html.position(target, true);
						var left = Math.floor(pos.x - this.arrowDim.w/2 + endAdd);

						html.style(this.bottomMarker, "visibility", "visible");
						html.style(this.topMarker, "visibility", "visible");
						html.style(this.bottomMarker, {
							"left": left + "px",
							"top" : (headerHeight + pos.y) + "px"
						});

						html.style(this.topMarker, {
							"left": left + "px",
							"top" : (pos.y - this.arrowDim.h) + "px"
						});

						if(src.targetAnchor && getSibling(src.targetAnchor, src.before)){
							src._addItemClass(getSibling(src.targetAnchor, src.before), src.before ? "After" : "Before");
						}
					}),
					_unmarkTargetAnchor: lang.hitch(this, function(){
						var src = this.source;
						if(!src.targetAnchor){ return; }
						if(src.targetAnchor && getSibling(src.targetAnchor, src.before)){
							src._removeItemClass(getSibling(src.targetAnchor, src.before), src.before ? "After" : "Before");
						}
						this._hide(this.bottomMarker);
						this._hide(this.topMarker);
						Source.prototype._unmarkTargetAnchor.call(src);
					}),
					destroy: lang.hitch(this, function(){
						connect.disconnect(this._source_conn);
						connect.unsubscribe(this._source_sub);
						Source.prototype.destroy.call(this.source);
						if(this.bottomMarker){
							html.destroy(this.bottomMarker);
							delete this.bottomMarker;
						}
						if(this.topMarker){
							html.destroy(this.topMarker);
							delete this.topMarker;
						}
					}),
					onDndCancel: lang.hitch(this, function(){
						Source.prototype.onDndCancel.call(this.source);
						this._hide(this.bottomMarker);
						this._hide(this.topMarker);
					})
				});

				this._source_conn = connect.connect(this.source, "onDndDrop", this, "_onDndDrop");
				this._source_sub = connect.subscribe("/dnd/drop/before", this, "_onDndDropBefore");
				this.source.startup();
			}
		},
		
		_hide: function(node){
			html.style(node, {
				top: "-10000px",
				"visibility": "hidden"
			});
		},

		_onDndDropBefore: function(source, nodes, copy){
			if(Manager.manager().target !== this.source){
				return;
			}
			this.source._targetNode = this.source.targetAnchor;
			this.source._beforeTarget = this.source.before;
			var views = this.grid.views.views;
			var srcView = views[source.viewIndex];
			var tgtView = views[this.index];
			if(tgtView != srcView){
				srcView.convertColPctToFixed();
				tgtView.convertColPctToFixed();
			}
		},

		_onDndDrop: function(source, nodes, copy){
			if(Manager.manager().target !== this.source){
				if(Manager.manager().source === this.source){
					this._removingColumn = true;
				}
				return;
			}
			this._hide(this.bottomMarker);
			this._hide(this.topMarker);

			var getIdx = function(n){
				return n ? html.attr(n, "idx") : null;
			};
			var w = html.marginBox(nodes[0]).w;
			if(source.viewIndex !== this.index){
				var views = this.grid.views.views;
				var srcView = views[source.viewIndex];
				var tgtView = views[this.index];
				if(srcView.viewWidth && srcView.viewWidth != "auto"){
					srcView.setColumnsWidth(srcView.getColumnsWidth() - w);
				}
				if(tgtView.viewWidth && tgtView.viewWidth != "auto"){
					tgtView.setColumnsWidth(tgtView.getColumnsWidth());
				}
			}
			var stn = this.source._targetNode;
			var stb = this.source._beforeTarget;
			!this.grid.isLeftToRight() && (stb = !stb);
			var layout = this.grid.layout;
			var idx = this.index;
			delete this.source._targetNode;
			delete this.source._beforeTarget;
			
			layout.moveColumn(
				source.viewIndex,
				idx,
				getIdx(nodes[0]),
				getIdx(stn),
				stb);
		},

		renderHeader: function(){
			this.headerContentNode.innerHTML = this.header.generateHtml(this._getHeaderContent);
			if(this.flexCells){
				this.contentWidth = this.getContentWidth();
				this.headerContentNode.firstChild.style.width = this.contentWidth;
			}
			util.fire(this, "onAfterRow", [-1, this.structure.cells, this.headerContentNode]);
		},

		// note: not called in 'view' context
		_getHeaderContent: function(inCell){
			var n = inCell.name || inCell.grid.getCellName(inCell);
			if(/^\s+$/.test(n)){
				n = '&nbsp;'//otherwise arrow styles will be messed up
			}
			var ret = [ '<div class="dojoxGridSortNode' ];
			
			if(inCell.index != inCell.grid.getSortIndex()){
				ret.push('">');
			}else{
				ret = ret.concat([ ' ',
							inCell.grid.sortInfo > 0 ? 'dojoxGridSortUp' : 'dojoxGridSortDown',
							'"><div class="dojoxGridArrowButtonChar">',
							inCell.grid.sortInfo > 0 ? '&#9650;' : '&#9660;',
							'</div><div class="dojoxGridArrowButtonNode" role="presentation"></div>',
							'<div class="dojoxGridColCaption">']);
			}
			ret = ret.concat([n, '</div></div>']);
			return ret.join('');
		},

		resize: function(){
			this.adaptHeight();
			this.adaptWidth();
		},

		hasHScrollbar: function(reset){
			var hadScroll = this._hasHScroll||false;
			if(this._hasHScroll == undefined || reset){
				if(this.noscroll){
					this._hasHScroll = false;
				}else{
					var style = html.style(this.scrollboxNode, "overflow");
					if(style == "hidden"){
						this._hasHScroll = false;
					}else if(style == "scroll"){
						this._hasHScroll = true;
					}else{
						this._hasHScroll = (this.scrollboxNode.offsetWidth - this.getScrollbarWidth() < this.contentNode.offsetWidth );
					}
				}
			}
			if(hadScroll !== this._hasHScroll){
				this.grid.update();
			}
			return this._hasHScroll; // Boolean
		},

		hasVScrollbar: function(reset){
			var hadScroll = this._hasVScroll||false;
			if(this._hasVScroll == undefined || reset){
				if(this.noscroll){
					this._hasVScroll = false;
				}else{
					var style = html.style(this.scrollboxNode, "overflow");
					if(style == "hidden"){
						this._hasVScroll = false;
					}else if(style == "scroll"){
						this._hasVScroll = true;
					}else{
						this._hasVScroll = (this.scrollboxNode.scrollHeight > this.scrollboxNode.clientHeight);
					}
				}
			}
			if(hadScroll !== this._hasVScroll){
				this.grid.update();
			}
			return this._hasVScroll; // Boolean
		},
		
		convertColPctToFixed: function(){
			// Fix any percentage widths to be pixel values
			var hasPct = false;
			this.grid.initialWidth = "";
			var cellNodes = query("th", this.headerContentNode);
			var fixedWidths = array.map(cellNodes, function(c, vIdx){
				var w = c.style.width;
				html.attr(c, "vIdx", vIdx);
				if(w && w.slice(-1) == "%"){
					hasPct = true;
				}else if(w && w.slice(-2) == "px"){
					return window.parseInt(w, 10);
				}
				return html.contentBox(c).w;
			});
			if(hasPct){
				array.forEach(this.grid.layout.cells, function(cell, idx){
					if(cell.view == this){
						var cellNode = cell.view.getHeaderCellNode(cell.index);
						if(cellNode && html.hasAttr(cellNode, "vIdx")){
							var vIdx = window.parseInt(html.attr(cellNode, "vIdx"));
							this.setColWidth(idx, fixedWidths[vIdx]);
							html.removeAttr(cellNode, "vIdx");
						}
					}
				}, this);
				return true;
			}
			return false;
		},

		adaptHeight: function(minusScroll){
			if(!this.grid._autoHeight){
				var h = (this.domNode.style.height && parseInt(this.domNode.style.height.replace(/px/,''), 10)) || this.domNode.clientHeight;
				var self = this;
				var checkOtherViewScrollers = function(){
					var v;
					for(var i in self.grid.views.views){
						v = self.grid.views.views[i];
						if(v !== self && v.hasHScrollbar()){
							return true;
						}
					}
					return false;
				};
				if(minusScroll || (this.noscroll && checkOtherViewScrollers())){
					h -= metrics.getScrollbar().h;
				}
				util.setStyleHeightPx(this.scrollboxNode, h);
			}
			this.hasVScrollbar(true);
		},

		adaptWidth: function(){
			if(this.flexCells){
				// the view content width
				this.contentWidth = this.getContentWidth();
				this.headerContentNode.firstChild.style.width = this.contentWidth;
			}
			// FIXME: it should be easier to get w from this.scrollboxNode.clientWidth,
			// but clientWidth seemingly does not include scrollbar width in some cases
			var w = this.scrollboxNode.offsetWidth - this.getScrollbarWidth();
			if(!this._removingColumn){
				w = Math.max(w, this.getColumnsWidth()) + 'px';
			}else{
				w = Math.min(w, this.getColumnsWidth()) + 'px';
				this._removingColumn = false;
			}
			var cn = this.contentNode;
			cn.style.width = w;
			this.hasHScrollbar(true);
		},

		setSize: function(w, h){
			var ds = this.domNode.style;
			var hs = this.headerNode.style;

			if(w){
				ds.width = w;
				hs.width = w;
			}
			ds.height = (h >= 0 ? h + 'px' : '');
		},

		renderRow: function(inRowIndex){
			var rowNode = this.createRowNode(inRowIndex);
			this.buildRow(inRowIndex, rowNode);
			//this.grid.edit.restore(this, inRowIndex);
			return rowNode;
		},

		createRowNode: function(inRowIndex){
			var node = document.createElement("div");
			node.className = this.classTag + 'Row';
			if (this instanceof dojox.grid._RowSelector){
				html.attr(node,"role","presentation");
			}else{
				html.attr(node,"role","row");
				if (this.grid.selectionMode != "none") {
					node.setAttribute("aria-selected", "false"); //rows can be selected so add aria-selected prop
				}
			}
			node[util.gridViewTag] = this.id;
			node[util.rowIndexTag] = inRowIndex;
			this.rowNodes[inRowIndex] = node;
			return node;
		},

		buildRow: function(inRowIndex, inRowNode){
			
			this.buildRowContent(inRowIndex, inRowNode);
		  	
			this.styleRow(inRowIndex, inRowNode);
		  
		 
		},

		buildRowContent: function(inRowIndex, inRowNode){
			inRowNode.innerHTML = this.content.generateHtml(inRowIndex, inRowIndex);
			if(this.flexCells && this.contentWidth){
				// FIXME: accessing firstChild here breaks encapsulation
				inRowNode.firstChild.style.width = this.contentWidth;
			}
			util.fire(this, "onAfterRow", [inRowIndex, this.structure.cells, inRowNode]);
		},

		rowRemoved:function(inRowIndex){
			if(inRowIndex >= 0){
				this._cleanupRowWidgets(this.getRowNode(inRowIndex));
			}
			this.grid.edit.save(this, inRowIndex);
			delete this.rowNodes[inRowIndex];
		},

		getRowNode: function(inRowIndex){
			return this.rowNodes[inRowIndex];
		},

		getCellNode: function(inRowIndex, inCellIndex){
			var row = this.getRowNode(inRowIndex);
			if(row){
				return this.content.getCellNode(row, inCellIndex);
			}
		},

		getHeaderCellNode: function(inCellIndex){
			if(this.headerContentNode){
				return this.header.getCellNode(this.headerContentNode, inCellIndex);
			}
		},

		// styling
		styleRow: function(inRowIndex, inRowNode){
			inRowNode._style = getStyleText(inRowNode);
			this.styleRowNode(inRowIndex, inRowNode);
		},

		styleRowNode: function(inRowIndex, inRowNode){
			if(inRowNode){
				this.doStyleRowNode(inRowIndex, inRowNode);
			}
		},

		doStyleRowNode: function(inRowIndex, inRowNode){
			this.grid.styleRowNode(inRowIndex, inRowNode);
		},

		// updating
		updateRow: function(inRowIndex){
			var rowNode = this.getRowNode(inRowIndex);
			if(rowNode){
				rowNode.style.height = '';
				this.buildRow(inRowIndex, rowNode);
			}
			return rowNode;
		},

		updateRowStyles: function(inRowIndex){
			this.styleRowNode(inRowIndex, this.getRowNode(inRowIndex));
		},

		// scrolling
		lastTop: 0,
		firstScroll:0,
		_nativeScroll: false,

		doscroll: function(inEvent){
			if(has('ff') >= 13){
				this._nativeScroll = true;
			}
			//var s = dojo.marginBox(this.headerContentNode.firstChild);
			var isLtr = this.grid.isLeftToRight();
			if(this.firstScroll < 2){
				if((!isLtr && this.firstScroll == 1) || (isLtr && this.firstScroll === 0)){
					var s = html.marginBox(this.headerNodeContainer);
					if(has('ie')){
						this.headerNodeContainer.style.width = s.w + this.getScrollbarWidth() + 'px';
					}else if(has('mozilla')){
						//TODO currently only for FF, not sure for safari and opera
						this.headerNodeContainer.style.width = s.w - this.getScrollbarWidth() + 'px';
						//this.headerNodeContainer.style.width = s.w + 'px';
						//set scroll to right in FF
						this.scrollboxNode.scrollLeft = isLtr ?
							this.scrollboxNode.clientWidth - this.scrollboxNode.scrollWidth :
							this.scrollboxNode.scrollWidth - this.scrollboxNode.clientWidth;
					}
				}
				this.firstScroll++;
			}
			this.headerNode.scrollLeft = this.scrollboxNode.scrollLeft;
			// 'lastTop' is a semaphore to prevent feedback-loop with setScrollTop below
			var top = this.scrollboxNode.scrollTop;
			if(top !== this.lastTop){
				this.grid.scrollTo(top);
			}
			this._nativeScroll = false;
		},

		setScrollTop: function(inTop){
			// 'lastTop' is a semaphore to prevent feedback-loop with doScroll above
			this.lastTop = inTop;
			if(!this._nativeScroll){
				//fix #15487
				this.scrollboxNode.scrollTop = inTop;
			}
			return this.scrollboxNode.scrollTop;
		},

		// event handlers (direct from DOM)
		doContentEvent: function(e){
			if(this.content.decorateEvent(e)){
				this.grid.onContentEvent(e);
			}
		},

		doHeaderEvent: function(e){
			if(this.header.decorateEvent(e)){
				this.grid.onHeaderEvent(e);
			}
		},

		// event dispatch(from Grid)
		dispatchContentEvent: function(e){
			return this.content.dispatchEvent(e);
		},

		dispatchHeaderEvent: function(e){
			return this.header.dispatchEvent(e);
		},

		// column resizing
		setColWidth: function(inIndex, inWidth){
			this.grid.setCellWidth(inIndex, inWidth + 'px');
		},

		update: function(){
			if(!this.domNode){
				return;
			}
			this.content.update();
			this.grid.update();
			//get scroll after update or scroll left setting goes wrong on IE.
			//See trac: #8040
			var left = this.scrollboxNode.scrollLeft;
			this.scrollboxNode.scrollLeft = left;
			this.headerNode.scrollLeft = left;
		}
	});

	var _GridAvatar = declare("dojox.grid._GridAvatar", Avatar, {
		construct: function(){
			var dd = win.doc;

			var a = dd.createElement("table");
			a.cellPadding = a.cellSpacing = "0";
			a.className = "dojoxGridDndAvatar";
			a.style.position = "absolute";
			a.style.zIndex = 1999;
			a.style.margin = "0px"; // to avoid dojo.marginBox() problems with table's margins
			var b = dd.createElement("tbody");
			var tr = dd.createElement("tr");
			var td = dd.createElement("td");
			var img = dd.createElement("td");
			tr.className = "dojoxGridDndAvatarItem";
			img.className = "dojoxGridDndAvatarItemImage";
			img.style.width = "16px";
			var source = this.manager.source, node;
			if(source.creator){
				// create an avatar representation of the node
				node = source._normalizedCreator(source.getItem(this.manager.nodes[0].id).data, "avatar").node;
			}else{
				// or just clone the node and hope it works
				node = this.manager.nodes[0].cloneNode(true);
				var table, tbody;
				if(node.tagName.toLowerCase() == "tr"){
					// insert extra table nodes
					table = dd.createElement("table");
					tbody = dd.createElement("tbody");
					tbody.appendChild(node);
					table.appendChild(tbody);
					node = table;
				}else if(node.tagName.toLowerCase() == "th"){
					// insert extra table nodes
					table = dd.createElement("table");
					tbody = dd.createElement("tbody");
					var r = dd.createElement("tr");
					table.cellPadding = table.cellSpacing = "0";
					r.appendChild(node);
					tbody.appendChild(r);
					table.appendChild(tbody);
					node = table;
				}
			}
			node.id = "";
			td.appendChild(node);
			tr.appendChild(img);
			tr.appendChild(td);
			html.style(tr, "opacity", 0.9);
			b.appendChild(tr);

			a.appendChild(b);
			this.node = a;

			var m = Manager.manager();
			this.oldOffsetY = m.OFFSET_Y;
			m.OFFSET_Y = 1;
		},
		destroy: function(){
			Manager.manager().OFFSET_Y = this.oldOffsetY;
			this.inherited(arguments);
		}
	});

	var oldMakeAvatar = Manager.manager().makeAvatar;
	Manager.manager().makeAvatar = function(){
		var src = this.source;
		if(src.viewIndex !== undefined && !html.hasClass(win.body(),"dijit_a11y")){
			return new _GridAvatar(this);
		}
		return oldMakeAvatar.call(Manager.manager());
	};

	return _View;

});
},
'url:dojox/grid/resources/View.html':"<div class=\"dojoxGridView\" role=\"presentation\">\n\t<div class=\"dojoxGridHeader\" dojoAttachPoint=\"headerNode\" role=\"presentation\">\n\t\t<div dojoAttachPoint=\"headerNodeContainer\" style=\"width:9000em\" role=\"presentation\">\n\t\t\t<div dojoAttachPoint=\"headerContentNode\" role=\"row\"></div>\n\t\t</div>\n\t</div>\n\t<input type=\"checkbox\" class=\"dojoxGridHiddenFocus\" dojoAttachPoint=\"hiddenFocusNode\" role=\"presentation\" />\n\t<input type=\"checkbox\" class=\"dojoxGridHiddenFocus\" role=\"presentation\" />\n\t<div class=\"dojoxGridScrollbox\" dojoAttachPoint=\"scrollboxNode\" role=\"presentation\">\n\t\t<div class=\"dojoxGridContent\" dojoAttachPoint=\"contentNode\" hidefocus=\"hidefocus\" role=\"presentation\"></div>\n\t</div>\n</div>\n",
'dojox/html/metrics':function(){
define("dojox/html/metrics", ["dojo/_base/kernel","dojo/_base/lang", "dojo/_base/sniff", "dojo/ready", "dojo/_base/unload",
		"dojo/_base/window", "dojo/dom-geometry"],
  function(kernel,lang,has,ready,UnloadUtil,Window,DOMGeom){
	var dhm = lang.getObject("dojox.html.metrics",true);
	var dojox = lang.getObject("dojox");

	//	derived from Morris John's emResized measurer
	dhm.getFontMeasurements = function(){
		// summary:
		//		Returns an object that has pixel equivilents of standard font size values.
		var heights = {
			'1em':0, '1ex':0, '100%':0, '12pt':0, '16px':0, 'xx-small':0, 'x-small':0,
			'small':0, 'medium':0, 'large':0, 'x-large':0, 'xx-large':0
		};
	
		if(has("ie")){
			//	we do a font-size fix if and only if one isn't applied already.
			//	NOTE: If someone set the fontSize on the HTML Element, this will kill it.
			Window.doc.documentElement.style.fontSize="100%";
		}
	
		//	set up the measuring node.
		var div=Window.doc.createElement("div");
		var ds = div.style;
		ds.position="absolute";
		ds.left="-100px";
		ds.top="0";
		ds.width="30px";
		ds.height="1000em";
		ds.borderWidth="0";
		ds.margin="0";
		ds.padding="0";
		ds.outline="0";
		ds.lineHeight="1";
		ds.overflow="hidden";
		Window.body().appendChild(div);
	
		//	do the measurements.
		for(var p in heights){
			ds.fontSize = p;
			heights[p] = Math.round(div.offsetHeight * 12/16) * 16/12 / 1000;
		}
		
		Window.body().removeChild(div);
		div = null;
		return heights; 	//	object
	};

	var fontMeasurements = null;
	
	dhm.getCachedFontMeasurements = function(recalculate){
		if(recalculate || !fontMeasurements){
			fontMeasurements = dhm.getFontMeasurements();
		}
		return fontMeasurements;
	};

	var measuringNode = null, empty = {};
	dhm.getTextBox = function(/* String */ text, /* Object */ style, /* String? */ className){
		var m, s;
		if(!measuringNode){
			m = measuringNode = Window.doc.createElement("div");
			// Container that we can set contraints on so that it doesn't
			// trigger a scrollbar.
			var c = Window.doc.createElement("div");
			c.appendChild(m);
			s = c.style;
			s.overflow='scroll';
			s.position = "absolute";
			s.left = "0px";
			s.top = "-10000px";
			s.width = "1px";
			s.height = "1px";
			s.visibility = "hidden";
			s.borderWidth = "0";
			s.margin = "0";
			s.padding = "0";
			s.outline = "0";
			Window.body().appendChild(c);
		}else{
			m = measuringNode;
		}
		// reset styles
		m.className = "";
		s = m.style;
		s.borderWidth = "0";
		s.margin = "0";
		s.padding = "0";
		s.outline = "0";
		// set new style
		if(arguments.length > 1 && style){
			for(var i in style){
				if(i in empty){ continue; }
				s[i] = style[i];
			}
		}
		// set classes
		if(arguments.length > 2 && className){
			m.className = className;
		}
		// take a measure
		m.innerHTML = text;
		var box = DOMGeom.position(m);
		// position doesn't report right (reports 1, since parent is 1)
		// So we have to look at the scrollWidth to get the real width
		// Height is right.
		box.w = m.parentNode.scrollWidth;
		return box;
	};

	//	determine the scrollbar sizes on load.
	var scroll={ w:16, h:16 };
	dhm.getScrollbar=function(){ return { w:scroll.w, h:scroll.h }; };

	dhm._fontResizeNode = null;

	dhm.initOnFontResize = function(interval){
		var f = dhm._fontResizeNode = Window.doc.createElement("iframe");
		var fs = f.style;
		fs.position = "absolute";
		fs.width = "5em";
		fs.height = "10em";
		fs.top = "-10000px";
		fs.display = "none";
		if(has("ie")){
			f.onreadystatechange = function(){
				if(f.contentWindow.document.readyState == "complete"){
					f.onresize = f.contentWindow.parent[dojox._scopeName].html.metrics._fontresize;
				}
			};
		}else{
			f.onload = function(){
				f.contentWindow.onresize = f.contentWindow.parent[dojox._scopeName].html.metrics._fontresize;
			};
		}
		//The script tag is to work around a known firebug race condition.  See comments in bug #9046
		f.setAttribute("src", "javascript:'<html><head><script>if(\"loadFirebugConsole\" in window){window.loadFirebugConsole();}</script></head><body></body></html>'");
		Window.body().appendChild(f);
		dhm.initOnFontResize = function(){};
	};

	dhm.onFontResize = function(){};
	dhm._fontresize = function(){
		dhm.onFontResize();
	};

	UnloadUtil.addOnUnload(function(){
		// destroy our font resize iframe if we have one
		var f = dhm._fontResizeNode;
		if(f){
			if(has("ie") && f.onresize){
				f.onresize = null;
			}else if(f.contentWindow && f.contentWindow.onresize){
				f.contentWindow.onresize = null;
			}
			dhm._fontResizeNode = null;
		}
	});

	ready(function(){
		// getScrollbar metrics node
		try{
			var n=Window.doc.createElement("div");
			n.style.cssText = "top:0;left:0;width:100px;height:100px;overflow:scroll;position:absolute;visibility:hidden;";
			Window.body().appendChild(n);
			scroll.w = n.offsetWidth - n.clientWidth;
			scroll.h = n.offsetHeight - n.clientHeight;
			Window.body().removeChild(n);
			//console.log("Scroll bar dimensions: ", scroll);
			delete n;
		}catch(e){}

		// text size poll setup
		if("fontSizeWatch" in kernel.config && !!kernel.config.fontSizeWatch){
			dhm.initOnFontResize();
		}
	});
	return dhm;
});
},
'dojox/grid/_Builder':function(){
define("dojox/grid/_Builder", [
	"../main",
	"dojo/_base/array",
	"dojo/_base/lang",
	"dojo/_base/window",
	"dojo/_base/event",
	"dojo/_base/sniff",
	"dojo/_base/connect",
	"dojo/dnd/Moveable",
	"dojox/html/metrics",
	"./util",
	"dojo/_base/html"
], function(dojox, array, lang, win, event, has, connect, Moveable, metrics, util, html){

	var dg = dojox.grid;

	var getTdIndex = function(td){
		return td.cellIndex >=0 ? td.cellIndex : array.indexOf(td.parentNode.cells, td);
	};
	
	var getTrIndex = function(tr){
		return tr.rowIndex >=0 ? tr.rowIndex : array.indexOf(tr.parentNode.childNodes, tr);
	};
	
	var getTr = function(rowOwner, index){
		return rowOwner && ((rowOwner.rows||0)[index] || rowOwner.childNodes[index]);
	};

	var findTable = function(node){
		for(var n=node; n && n.tagName!='TABLE'; n=n.parentNode){}
		return n;
	};
	
	var ascendDom = function(inNode, inWhile){
		for(var n=inNode; n && inWhile(n); n=n.parentNode){}
		return n;
	};
	
	var makeNotTagName = function(inTagName){
		var name = inTagName.toUpperCase();
		return function(node){ return node.tagName != name; };
	};

	var rowIndexTag = util.rowIndexTag;
	var gridViewTag = util.gridViewTag;

	// base class for generating markup for the views
	var _Builder = dg._Builder = lang.extend(function(view){
		if(view){
			this.view = view;
			this.grid = view.grid;
		}
	},{
		view: null,
		// boilerplate HTML
		_table: '<table class="dojoxGridRowTable" border="0" cellspacing="0" cellpadding="0" role="presentation"',

		// Returns the table variable as an array - and with the view width, if specified
		getTableArray: function(){
			var html = [this._table];
			if(this.view.viewWidth){
				html.push([' style="width:', this.view.viewWidth, ';"'].join(''));
			}
			html.push('>');
			return html;
		},
		
		// generate starting tags for a cell
		generateCellMarkup: function(inCell, inMoreStyles, inMoreClasses, isHeader){
			var result = [], html;
			if(isHeader){
				var sortInfo = inCell.index != inCell.grid.getSortIndex() ? "" : inCell.grid.sortInfo > 0 ? 'aria-sort="ascending"' : 'aria-sort="descending"';
				if (!inCell.id){
					inCell.id = this.grid.id + "Hdr" + inCell.index;
				}
				// column headers are not editable, mark as aria-readonly=true
				html = ['<th tabIndex="-1" aria-readonly="true" role="columnheader"', sortInfo, 'id="', inCell.id, '"'];
			}else{
				// cells inherit grid aria-readonly property; default value for aria-readonly is false(grid is editable)
				// if grid is editable (had any editable cells), mark non editable cells as aria-readonly=true
				// if no editable cells, grid's aria-readonly value will have been set to true and cells will inherit
				var editInfo = this.grid.editable && !inCell.editable ? 'aria-readonly="true"' : "";
				html = ['<td tabIndex="-1" role="gridcell"', editInfo];
			}
			if(inCell.colSpan){
				html.push(' colspan="', inCell.colSpan, '"');
			}
			if(inCell.rowSpan){
				html.push(' rowspan="', inCell.rowSpan, '"');
			}
			html.push(' class="dojoxGridCell ');
			if(inCell.classes){
				html.push(inCell.classes, ' ');
			}
			if(inMoreClasses){
				html.push(inMoreClasses, ' ');
			}
			// result[0] => td opener, style
			result.push(html.join(''));
			// SLOT: result[1] => td classes
			result.push('');
			html = ['" idx="', inCell.index, '" style="'];
			if(inMoreStyles && inMoreStyles[inMoreStyles.length-1] != ';'){
				inMoreStyles += ';';
			}
			html.push(inCell.styles, inMoreStyles||'', inCell.hidden?'display:none;':'');
			if(inCell.unitWidth){
				html.push('width:', inCell.unitWidth, ';');
			}
			// result[2] => markup
			result.push(html.join(''));
			// SLOT: result[3] => td style
			result.push('');
			html = [ '"' ];
			if(inCell.attrs){
				html.push(" ", inCell.attrs);
			}
			html.push('>');
			// result[4] => td postfix
			result.push(html.join(''));
			// SLOT: result[5] => content
			result.push('');
			// result[6] => td closes
			result.push(isHeader?'</th>':'</td>');
			return result; // Array
		},

		// cell finding
		isCellNode: function(inNode){
			return Boolean(inNode && inNode!=win.doc && html.attr(inNode, "idx"));
		},
		
		getCellNodeIndex: function(inCellNode){
			return inCellNode ? Number(html.attr(inCellNode, "idx")) : -1;
		},
		
		getCellNode: function(inRowNode, inCellIndex){
			for(var i=0, row; ((row = getTr(inRowNode.firstChild, i)) && row.cells); i++){
				for(var j=0, cell; (cell = row.cells[j]); j++){
					if(this.getCellNodeIndex(cell) == inCellIndex){
						return cell;
					}
				}
			}
			return null;
		},
		
		findCellTarget: function(inSourceNode, inTopNode){
			var n = inSourceNode;
			while(n && (!this.isCellNode(n) || (n.offsetParent && gridViewTag in n.offsetParent.parentNode && n.offsetParent.parentNode[gridViewTag] != this.view.id)) && (n!=inTopNode)){
				n = n.parentNode;
			}
			return n!=inTopNode ? n : null;
		},
		
		// event decoration
		baseDecorateEvent: function(e){
			e.dispatch = 'do' + e.type;
			e.grid = this.grid;
			e.sourceView = this.view;
			e.cellNode = this.findCellTarget(e.target, e.rowNode);
			e.cellIndex = this.getCellNodeIndex(e.cellNode);
			e.cell = (e.cellIndex >= 0 ? this.grid.getCell(e.cellIndex) : null);
		},
		
		// event dispatch
		findTarget: function(inSource, inTag){
			var n = inSource;
			while(n && (n!=this.domNode) && (!(inTag in n) || (gridViewTag in n && n[gridViewTag] != this.view.id))){
				n = n.parentNode;
			}
			return (n != this.domNode) ? n : null;
		},

		findRowTarget: function(inSource){
			return this.findTarget(inSource, rowIndexTag);
		},

		isIntraNodeEvent: function(e){
			try{
				return (e.cellNode && e.relatedTarget && html.isDescendant(e.relatedTarget, e.cellNode));
			}catch(x){
				// e.relatedTarget has permission problem in FF if it's an input: https://bugzilla.mozilla.org/show_bug.cgi?id=208427
				return false;
			}
		},

		isIntraRowEvent: function(e){
			try{
				var row = e.relatedTarget && this.findRowTarget(e.relatedTarget);
				return !row && (e.rowIndex==-1) || row && (e.rowIndex==row.gridRowIndex);
			}catch(x){
				// e.relatedTarget on INPUT has permission problem in FF: https://bugzilla.mozilla.org/show_bug.cgi?id=208427
				return false;
			}
		},

		dispatchEvent: function(e){
			if(e.dispatch in this){
				return this[e.dispatch](e);
			}
			return false;
		},

		// dispatched event handlers
		domouseover: function(e){
			if(e.cellNode && (e.cellNode!=this.lastOverCellNode)){
				this.lastOverCellNode = e.cellNode;
				this.grid.onMouseOver(e);
			}
			this.grid.onMouseOverRow(e);
		},

		domouseout: function(e){
			if(e.cellNode && (e.cellNode==this.lastOverCellNode) && !this.isIntraNodeEvent(e, this.lastOverCellNode)){
				this.lastOverCellNode = null;
				this.grid.onMouseOut(e);
				if(!this.isIntraRowEvent(e)){
					this.grid.onMouseOutRow(e);
				}
			}
		},
		
		domousedown: function(e){
			if (e.cellNode)
				this.grid.onMouseDown(e);
			this.grid.onMouseDownRow(e);
		}
	});

	// Produces html for grid data content. Owned by grid and used internally
	// for rendering data. Override to implement custom rendering.
	var _ContentBuilder = dg._ContentBuilder = lang.extend(function(view){
		_Builder.call(this, view);
	},_Builder.prototype,{
		update: function(){
			this.prepareHtml();
		},

		// cache html for rendering data rows
		prepareHtml: function(){
			var defaultGet=this.grid.get, cells=this.view.structure.cells;
			for(var j=0, row; (row=cells[j]); j++){
				for(var i=0, cell; (cell=row[i]); i++){
					cell.get = cell.get || (cell.value == undefined) && defaultGet;
					cell.markup = this.generateCellMarkup(cell, cell.cellStyles, cell.cellClasses, false);
					if (!this.grid.editable && cell.editable){
						this.grid.editable = true;
					}
				}
			}
		},

		// time critical: generate html using cache and data source
		generateHtml: function(inDataIndex, inRowIndex){
			var
				html = this.getTableArray(),
				v = this.view,
				cells = v.structure.cells,
				item = this.grid.getItem(inRowIndex);

			util.fire(this.view, "onBeforeRow", [inRowIndex, cells]);
			for(var j=0, row; (row=cells[j]); j++){
				if(row.hidden || row.header){
					continue;
				}
				html.push(!row.invisible ? '<tr>' : '<tr class="dojoxGridInvisible">');
				for(var i=0, cell, m, cc, cs; (cell=row[i]); i++){
					m = cell.markup; cc = cell.customClasses = []; cs = cell.customStyles = [];
					// content (format can fill in cc and cs as side-effects)
					m[5] = cell.format(inRowIndex, item);
					// classes
					m[1] = cc.join(' ');
					// styles
					m[3] = cs.join(';');
					// in-place concat
					html.push.apply(html, m);
				}
				html.push('</tr>');
			}
			html.push('</table>');
			return html.join(''); // String
		},

		decorateEvent: function(e){
			e.rowNode = this.findRowTarget(e.target);
			if(!e.rowNode){return false;}
			e.rowIndex = e.rowNode[rowIndexTag];
			this.baseDecorateEvent(e);
			e.cell = this.grid.getCell(e.cellIndex);
			return true; // Boolean
		}
	});

	// Produces html for grid header content. Owned by grid and used internally
	// for rendering data. Override to implement custom rendering.
	var _HeaderBuilder = dg._HeaderBuilder = lang.extend(function(view){
		this.moveable = null;
		_Builder.call(this, view);
	},_Builder.prototype,{
		_skipBogusClicks: false,
		overResizeWidth: 4,
		minColWidth: 1,
		
		update: function(){
			if(this.tableMap){
				this.tableMap.mapRows(this.view.structure.cells);
			}else{
				this.tableMap = new dg._TableMap(this.view.structure.cells);
			}
		},

		generateHtml: function(inGetValue, inValue){
			var html = this.getTableArray(), cells = this.view.structure.cells;
			
			util.fire(this.view, "onBeforeRow", [-1, cells]);
			for(var j=0, row; (row=cells[j]); j++){
				if(row.hidden){
					continue;
				}
				html.push(!row.invisible ? '<tr>' : '<tr class="dojoxGridInvisible">');
				for(var i=0, cell, markup; (cell=row[i]); i++){
					cell.customClasses = [];
					cell.customStyles = [];
					if(this.view.simpleStructure){
						if(cell.draggable){
							if(cell.headerClasses){
								if(cell.headerClasses.indexOf('dojoDndItem') == -1){
									cell.headerClasses += ' dojoDndItem';
								}
							}else{
								cell.headerClasses = 'dojoDndItem';
							}
						}
						if(cell.attrs){
							if(cell.attrs.indexOf("dndType='gridColumn_") == -1){
								cell.attrs += " dndType='gridColumn_" + this.grid.id + "'";
							}
						}else{
							cell.attrs = "dndType='gridColumn_" + this.grid.id + "'";
						}
					}
					markup = this.generateCellMarkup(cell, cell.headerStyles, cell.headerClasses, true);
					// content
					markup[5] = (inValue != undefined ? inValue : inGetValue(cell));
					// styles
					markup[3] = cell.customStyles.join(';');
					// classes
					markup[1] = cell.customClasses.join(' '); //(cell.customClasses ? ' ' + cell.customClasses : '');
					html.push(markup.join(''));
				}
				html.push('</tr>');
			}
			html.push('</table>');
			return html.join('');
		},

		// event helpers
		getCellX: function(e){
			var n, x = e.layerX;
			if(has('mozilla') || has('ie') >= 9){
				n = ascendDom(e.target, makeNotTagName("th"));
				x -= (n && n.offsetLeft) || 0;
				var t = e.sourceView.getScrollbarWidth();
				if(!this.grid.isLeftToRight()/*&& e.sourceView.headerNode.scrollLeft < t*/){
					//fix #11253
					table = ascendDom(n,makeNotTagName("table"));
					x -= (table && table.offsetLeft) || 0;
				}
				//x -= getProp(ascendDom(e.target, mkNotTagName("td")), "offsetLeft") || 0;
			}
			n = ascendDom(e.target, function(){
				if(!n || n == e.cellNode){
					return false;
				}
				// Mozilla 1.8 (FF 1.5) has a bug that makes offsetLeft = -parent border width
				// when parent has border, overflow: hidden, and is positioned
				// handle this problem here ... not a general solution!
				x += (n.offsetLeft < 0 ? 0 : n.offsetLeft);
				return true;
			});
			return x;
		},

		// event decoration
		decorateEvent: function(e){
			this.baseDecorateEvent(e);
			e.rowIndex = -1;
			e.cellX = this.getCellX(e);
			return true;
		},

		// event handlers
		// resizing
		prepareResize: function(e, mod){
			do{
				var i = e.cellIndex;
				e.cellNode = (i ? e.cellNode.parentNode.cells[i+mod] : null);
				e.cellIndex = (e.cellNode ? this.getCellNodeIndex(e.cellNode) : -1);
			}while(e.cellNode && e.cellNode.style.display == "none");
			return Boolean(e.cellNode);
		},

		canResize: function(e){
			if(!e.cellNode || e.cellNode.colSpan > 1){
				return false;
			}
			var cell = this.grid.getCell(e.cellIndex);
			return !cell.noresize && cell.canResize();
		},

		overLeftResizeArea: function(e){
			// We are never over a resize area if we are in the process of moving
			if(html.hasClass(win.body(), "dojoDndMove")){
				return false;
			}
			//Bugfix for crazy IE problem (#8807).  IE returns position information for the icon and text arrow divs
			//as if they were still on the left instead of returning the position they were 'float: right' to.
			//So, the resize check ends up checking the wrong adjacent cell.  This checks to see if the hover was over
			//the image or text nodes, then just ignored them/treat them not in scale range.
			if(has('ie')){
				var tN = e.target;
				if(html.hasClass(tN, "dojoxGridArrowButtonNode") ||
					html.hasClass(tN, "dojoxGridArrowButtonChar") ||
					html.hasClass(tN, "dojoxGridColCaption")){
					return false;
				}
			}

			if(this.grid.isLeftToRight()){
				return (e.cellIndex>0) && (e.cellX > 0 && e.cellX < this.overResizeWidth) && this.prepareResize(e, -1);
			}
			var t = e.cellNode && (e.cellX > 0 && e.cellX < this.overResizeWidth);
			return t;
		},

		overRightResizeArea: function(e){
			// We are never over a resize area if we are in the process of moving
			if(html.hasClass(win.body(), "dojoDndMove")){
				return false;
			}
			//Bugfix for crazy IE problem (#8807).  IE returns position information for the icon and text arrow divs
			//as if they were still on the left instead of returning the position they were 'float: right' to.
			//So, the resize check ends up checking the wrong adjacent cell.  This checks to see if the hover was over
			//the image or text nodes, then just ignored them/treat them not in scale range.
			if(has('ie')){
				var tN = e.target;
				if(html.hasClass(tN, "dojoxGridArrowButtonNode") ||
					html.hasClass(tN, "dojoxGridArrowButtonChar") ||
					html.hasClass(tN, "dojoxGridColCaption")){
					return false;
				}
			}

			if(this.grid.isLeftToRight()){
				return e.cellNode && (e.cellX >= e.cellNode.offsetWidth - this.overResizeWidth);
			}
			return (e.cellIndex>0) && (e.cellX >= e.cellNode.offsetWidth - this.overResizeWidth) && this.prepareResize(e, -1);
		},

		domousemove: function(e){
			//console.log(e.cellIndex, e.cellX, e.cellNode.offsetWidth);
			if(!this.moveable){
				var c = (this.overRightResizeArea(e) ? 'dojoxGridColResize' : (this.overLeftResizeArea(e) ? 'dojoxGridColResize' : ''));
				if(c && !this.canResize(e)){
					c = 'dojoxGridColNoResize';
				}
				html.toggleClass(e.sourceView.headerNode, "dojoxGridColNoResize", (c == "dojoxGridColNoResize"));
				html.toggleClass(e.sourceView.headerNode, "dojoxGridColResize", (c == "dojoxGridColResize"));
				if(c){
					event.stop(e);
				}
			}
		},

		domousedown: function(e){
			if(!this.moveable){
				if((this.overRightResizeArea(e) || this.overLeftResizeArea(e)) && this.canResize(e)){
					this.beginColumnResize(e);
				}else{
					this.grid.onMouseDown(e);
					this.grid.onMouseOverRow(e);
				}
				//else{
				//	this.beginMoveColumn(e);
				//}
			}
		},

		doclick: function(e) {
			if(this._skipBogusClicks){
				event.stop(e);
				return true;
			}
			return false;
		},

		// column resizing
		colResizeSetup: function(/*Event Object*/e, /*boolean*/ isMouse ){
			//Set up the drag object for column resizing
			// Called with mouse event in case of drag and drop,
			// Also called from keyboard shift-arrow event when focus is on a header
			var headContentBox = html.contentBox(e.sourceView.headerNode);
			
			if(isMouse){  //IE draws line even with no mouse down so separate from keyboard
				this.lineDiv = document.createElement('div');

				var vw = html.position(e.sourceView.headerNode, true);
				var bodyContentBox = html.contentBox(e.sourceView.domNode);
				//fix #11340
				var l = e.pageX;
				if(!this.grid.isLeftToRight() && has('ie') < 8){
					l -= metrics.getScrollbar().w;
				}
				html.style(this.lineDiv, {
					top: vw.y + "px",
					left: l + "px",
					height: (bodyContentBox.h + headContentBox.h) + "px"
				});
				html.addClass(this.lineDiv, "dojoxGridResizeColLine");
				this.lineDiv._origLeft = l;
				win.body().appendChild(this.lineDiv);
			}
			var spanners = [], nodes = this.tableMap.findOverlappingNodes(e.cellNode);
			for(var i=0, cell; (cell=nodes[i]); i++){
				spanners.push({ node: cell, index: this.getCellNodeIndex(cell), width: cell.offsetWidth });
				//console.log("spanner: " + this.getCellNodeIndex(cell));
			}

			var view = e.sourceView;
			var adj = this.grid.isLeftToRight() ? 1 : -1;
			var views = e.grid.views.views;
			var followers = [];
			for(var j=view.idx+adj, cView; (cView=views[j]); j=j+adj){
				followers.push({ node: cView.headerNode, left: window.parseInt(cView.headerNode.style.left) });
			}
			var table = view.headerContentNode.firstChild;
			var drag = {
				scrollLeft: e.sourceView.headerNode.scrollLeft,
				view: view,
				node: e.cellNode,
				index: e.cellIndex,
				w: html.contentBox(e.cellNode).w,
				vw: headContentBox.w,
				table: table,
				tw: html.contentBox(table).w,
				spanners: spanners,
				followers: followers
			};
			return drag;
		},
		beginColumnResize: function(e){
			this.moverDiv = document.createElement("div");
			html.style(this.moverDiv,{position: "absolute", left:0}); // to make DnD work with dir=rtl
			win.body().appendChild(this.moverDiv);
			html.addClass(this.grid.domNode, "dojoxGridColumnResizing");
			var m = (this.moveable = new Moveable(this.moverDiv));

			var drag = this.colResizeSetup(e,true);

			m.onMove = lang.hitch(this, "doResizeColumn", drag);

			connect.connect(m, "onMoveStop", lang.hitch(this, function(){
				this.endResizeColumn(drag);
				if(drag.node.releaseCapture){
					drag.node.releaseCapture();
				}
				this.moveable.destroy();
				delete this.moveable;
				this.moveable = null;
				html.removeClass(this.grid.domNode, "dojoxGridColumnResizing");
			}));

			if(e.cellNode.setCapture){
				e.cellNode.setCapture();
			}
			m.onMouseDown(e);
		},

		doResizeColumn: function(inDrag, mover, leftTop){
			var changeX = leftTop.l;
			var data = {
				deltaX: changeX,
				w: inDrag.w + (this.grid.isLeftToRight() ? changeX : -changeX),//fix #11341
				vw: inDrag.vw + changeX,
				tw: inDrag.tw + changeX
			};
			
			this.dragRecord = {inDrag: inDrag, mover: mover, leftTop:leftTop};
			
			if(data.w >= this.minColWidth){
				if (!mover) { // we are using keyboard do immediate resize
					this.doResizeNow(inDrag, data);
				}
				else{
					html.style(this.lineDiv, "left", (this.lineDiv._origLeft + data.deltaX) + "px");
				}
			}
		},

		endResizeColumn: function(inDrag){
			if(this.dragRecord){
				var leftTop = this.dragRecord.leftTop;
				var changeX = this.grid.isLeftToRight() ? leftTop.l : -leftTop.l;
				// Make sure we are not under our minimum
				// http://bugs.dojotoolkit.org/ticket/9390
				changeX += Math.max(inDrag.w + changeX, this.minColWidth) - (inDrag.w + changeX);
				if(has('webkit') && inDrag.spanners.length){
					// Webkit needs the pad border extents back in
					changeX += html._getPadBorderExtents(inDrag.spanners[0].node).w;
				}
				var data = {
					deltaX: changeX,
					w: inDrag.w + changeX,
					vw: inDrag.vw + changeX,
					tw: inDrag.tw + changeX
				};
				// Only resize the columns when the drag has finished
				this.doResizeNow(inDrag, data);
				delete this.dragRecord;
			}
			
			html.destroy(this.lineDiv);
 			html.destroy(this.moverDiv);
			html.destroy(this.moverDiv);
			delete this.moverDiv;
			this._skipBogusClicks = true;
			inDrag.view.update();
			this._skipBogusClicks = false;
			this.grid.onResizeColumn(inDrag.index);
		},
		doResizeNow: function(inDrag, data){
			inDrag.view.convertColPctToFixed();
			if(inDrag.view.flexCells && !inDrag.view.testFlexCells()){
				var t = findTable(inDrag.node);
				if(t){
					(t.style.width = '');
				}
			}
			var i, s, sw, f, fl;
			for(i=0; (s=inDrag.spanners[i]); i++){
				sw = s.width + data.deltaX;
				if(sw > 0){
					s.node.style.width = sw + 'px';
					inDrag.view.setColWidth(s.index, sw);
				}
			}
			if(this.grid.isLeftToRight() || !has('ie')){//fix #11339
				for(i=0; (f=inDrag.followers[i]); i++){
					fl = f.left + data.deltaX;
					f.node.style.left = fl + 'px';
				}
			}
			inDrag.node.style.width = data.w + 'px';
			inDrag.view.setColWidth(inDrag.index, data.w);
			inDrag.view.headerNode.style.width = data.vw + 'px';
			inDrag.view.setColumnsWidth(data.tw);
			if(!this.grid.isLeftToRight()){
				inDrag.view.headerNode.scrollLeft = inDrag.scrollLeft + data.deltaX;
			}
		}
	});

	// Maps an html table into a structure parsable for information about cell row and col spanning.
	// Used by HeaderBuilder.
	dg._TableMap = lang.extend(function(rows){
		this.mapRows(rows);
	},{
		map: null,

		mapRows: function(inRows){
			// summary:
			//		Map table topography

			//console.log('mapRows');
			// # of rows
			var rowCount = inRows.length;
			if(!rowCount){
				return;
			}
			// map which columns and rows fill which cells
			this.map = [];
			var row;
			for(var k=0; (row=inRows[k]); k++){
				this.map[k] = [];
			}
			for(var j=0; (row=inRows[j]); j++){
				for(var i=0, x=0, cell, colSpan, rowSpan; (cell=row[i]); i++){
					while(this.map[j][x]){x++;}
					this.map[j][x] = { c: i, r: j };
					rowSpan = cell.rowSpan || 1;
					colSpan = cell.colSpan || 1;
					for(var y=0; y<rowSpan; y++){
						for(var s=0; s<colSpan; s++){
							this.map[j+y][x+s] = this.map[j][x];
						}
					}
					x += colSpan;
				}
			}
			//this.dumMap();
		},

		dumpMap: function(){
			for(var j=0, row, h=''; (row=this.map[j]); j++,h=''){
				for(var i=0, cell; (cell=row[i]); i++){
					h += cell.r + ',' + cell.c + '   ';
				}
			}
		},

		getMapCoords: function(inRow, inCol){
			// summary:
			//		Find node's map coords by it's structure coords
			for(var j=0, row; (row=this.map[j]); j++){
				for(var i=0, cell; (cell=row[i]); i++){
					if(cell.c==inCol && cell.r == inRow){
						return { j: j, i: i };
					}
					//else{console.log(inRow, inCol, ' : ', i, j, " : ", cell.r, cell.c); };
				}
			}
			return { j: -1, i: -1 };
		},
		
		getNode: function(inTable, inRow, inCol){
			// summary:
			//		Find a node in inNode's table with the given structure coords
			var row = inTable && inTable.rows[inRow];
			return row && row.cells[inCol];
		},
		
		_findOverlappingNodes: function(inTable, inRow, inCol){
			var nodes = [];
			var m = this.getMapCoords(inRow, inCol);
			//console.log("node j: %d, i: %d", m.j, m.i);
			for(var j=0, row; (row=this.map[j]); j++){
				if(j == m.j){ continue; }
				var rw = row[m.i];
				//console.log("overlaps: r: %d, c: %d", rw.r, rw.c);
				var n = (rw?this.getNode(inTable, rw.r, rw.c):null);
				if(n){ nodes.push(n); }
			}
			//console.log(nodes);
			return nodes;
		},
		
		findOverlappingNodes: function(inNode){
			return this._findOverlappingNodes(findTable(inNode), getTrIndex(inNode.parentNode), getTdIndex(inNode));
		}
	});

	return {
		_Builder: _Builder,
		_HeaderBuilder: _HeaderBuilder,
		_ContentBuilder: _ContentBuilder
	};
});
},
'dojox/grid/_ViewManager':function(){
define("dojox/grid/_ViewManager", [
	"dojo/_base/declare",
	"dojo/_base/sniff",
	"dojo/dom-class"
], function(declare, has, domClass){

return declare('dojox.grid._ViewManager', null, {
	// summary:
	//		A collection of grid views. Owned by grid and used internally for managing grid views.
	// description:
	//		Grid creates views automatically based on grid's layout structure.
	//		Users should typically not need to access individual views or the views collection directly.
	constructor: function(inGrid){
		this.grid = inGrid;
	},

	defaultWidth: 200,

	views: [],

	// operations
	resize: function(){
		this.onEach("resize");
	},

	render: function(){
		this.onEach("render");
	},

	// views
	addView: function(inView){
		inView.idx = this.views.length;
		this.views.push(inView);
	},

	destroyViews: function(){
		for(var i=0, v; v=this.views[i]; i++){
			v.destroy();
		}
		this.views = [];
	},

	getContentNodes: function(){
		var nodes = [];
		for(var i=0, v; v=this.views[i]; i++){
			nodes.push(v.contentNode);
		}
		return nodes;
	},

	forEach: function(inCallback){
		for(var i=0, v; v=this.views[i]; i++){
			inCallback(v, i);
		}
	},

	onEach: function(inMethod, inArgs){
		inArgs = inArgs || [];
		for(var i=0, v; v=this.views[i]; i++){
			if(inMethod in v){
				v[inMethod].apply(v, inArgs);
			}
		}
	},

	// layout
	normalizeHeaderNodeHeight: function(){
		var rowNodes = [];
		for(var i=0, v; (v=this.views[i]); i++){
			if(v.headerContentNode.firstChild){
				rowNodes.push(v.headerContentNode);
			}
		}
		this.normalizeRowNodeHeights(rowNodes);
	},

	normalizeRowNodeHeights: function(inRowNodes){
		var h = 0;
		var currHeights = [];
		if(this.grid.rowHeight){
			h = this.grid.rowHeight;
		}else{
			if(inRowNodes.length <= 1){
				// no need to normalize if we are the only one...
				return;
			}
			for(var i=0, n; (n=inRowNodes[i]); i++){
				// We only care about the height - so don't use marginBox.  This
				// depends on the container not having any margin (which it shouldn't)
				// Also - we only look up the height if the cell doesn't have the
				// dojoxGridNonNormalizedCell class (like for row selectors)
				if(!domClass.contains(n, "dojoxGridNonNormalizedCell")){
					currHeights[i] = n.firstChild.offsetHeight;
					h =  Math.max(h, currHeights[i]);
				}
			}
			h = (h >= 0 ? h : 0);
	
			//Work around odd FF3 rendering bug: #8864.
			//A one px increase fixes FireFox 3's rounding bug for fractional font sizes.
			if((has('mozilla') || has('ie') > 8 ) && h){h++;}
		}
		for(i=0; (n=inRowNodes[i]); i++){
			if(currHeights[i] != h){
				n.firstChild.style.height = h + "px";
			}
		}
	},
	
	resetHeaderNodeHeight: function(){
		for(var i=0, v, n; (v=this.views[i]); i++){
			n = v.headerContentNode.firstChild;
			if(n){
				n.style.height = "";
			}
		}
	},

	renormalizeRow: function(inRowIndex){
		var rowNodes = [];
		for(var i=0, v, n; (v=this.views[i])&&(n=v.getRowNode(inRowIndex)); i++){
			n.firstChild.style.height = '';
			rowNodes.push(n);
		}
		this.normalizeRowNodeHeights(rowNodes);
	},

	getViewWidth: function(inIndex){
		return this.views[inIndex].getWidth() || this.defaultWidth;
	},

	// must be called after view widths are properly set or height can be miscalculated
	// if there are flex columns
	measureHeader: function(){
		// need to reset view header heights so they are properly measured.
		this.resetHeaderNodeHeight();
		this.forEach(function(inView){
			inView.headerContentNode.style.height = '';
		});
		var h = 0;
		// calculate maximum view header height
		this.forEach(function(inView){
			h = Math.max(inView.headerNode.offsetHeight, h);
		});
		return h;
	},

	measureContent: function(){
		var h = 0;
		this.forEach(function(inView){
			h = Math.max(inView.domNode.offsetHeight, h);
		});
		return h;
	},

	findClient: function(inAutoWidth){
		// try to use user defined client
		var c = this.grid.elasticView || -1;
		// attempt to find implicit client
		if(c < 0){
			for(var i=1, v; (v=this.views[i]); i++){
				if(v.viewWidth){
					for(i=1; (v=this.views[i]); i++){
						if(!v.viewWidth){
							c = i;
							break;
						}
					}
					break;
				}
			}
		}
		// client is in the middle by default
		if(c < 0){
			c = Math.floor(this.views.length / 2);
		}
		return c;
	},

	arrange: function(l, w){
		var i, v, vw, len = this.views.length, self = this;
		// find the client
		var c = (w <= 0 ? len : this.findClient());
		// layout views
		var setPosition = function(v, l){
			var ds = v.domNode.style;
			var hs = v.headerNode.style;

			if(!self.grid.isLeftToRight()){
				ds.right = l + 'px';
				// fixed rtl, the scrollbar is on the right side in FF < 4
				if(has('ff') < 4){
					hs.right = l + v.getScrollbarWidth() + 'px';
				}else{
					hs.right = l + 'px';
				}
				if(!has('webkit')){
					hs.width = parseInt(hs.width, 10) - v.getScrollbarWidth() + 'px';					
				}
			}else{
				ds.left = l + 'px';
				hs.left = l + 'px';
			}
			ds.top = 0 + 'px';
			hs.top = 0;
		};
		// for views left of the client
		//BiDi TODO: The left and right should not appear in BIDI environment. Should be replaced with
		//leading and tailing concept.
		for(i=0; (v=this.views[i])&&(i<c); i++){
			// get width
			vw = this.getViewWidth(i);
			// process boxes
			v.setSize(vw, 0);
			setPosition(v, l);
			if(v.headerContentNode && v.headerContentNode.firstChild){
				vw = v.getColumnsWidth()+v.getScrollbarWidth();
			}else{
				vw = v.domNode.offsetWidth;
			}
			// update position
			l += vw;
		}
		// next view (is the client, i++ == c)
		i++;
		// start from the right edge
		var r = w;
		// for views right of the client (iterated from the right)
		for(var j=len-1; (v=this.views[j])&&(i<=j); j--){
			// get width
			vw = this.getViewWidth(j);
			// set size
			v.setSize(vw, 0);
			// measure in pixels
			vw = v.domNode.offsetWidth;
			// update position
			r -= vw;
			// set position
			setPosition(v, r);
		}
		if(c<len){
			v = this.views[c];
			// position the client box between left and right boxes
			vw = Math.max(1, r-l);
			// set size
			v.setSize(vw + 'px', 0);
			setPosition(v, l);
		}
		return l;
	},

	// rendering
	renderRow: function(inRowIndex, inNodes, skipRenorm){
		var rowNodes = [];
		for(var i=0, v, n, rowNode; (v=this.views[i])&&(n=inNodes[i]); i++){
			rowNode = v.renderRow(inRowIndex);
			n.appendChild(rowNode);
			rowNodes.push(rowNode);
		}
		if(!skipRenorm){
			this.normalizeRowNodeHeights(rowNodes);
		}
	},
	
	rowRemoved: function(inRowIndex){
		this.onEach("rowRemoved", [ inRowIndex ]);
	},
	
	// updating
	updateRow: function(inRowIndex, skipRenorm){
		for(var i=0, v; v=this.views[i]; i++){
			v.updateRow(inRowIndex);
		}
		if(!skipRenorm){
			this.renormalizeRow(inRowIndex);
		}
	},
	
	updateRowStyles: function(inRowIndex){
		this.onEach("updateRowStyles", [ inRowIndex ]);
	},
	
	// scrolling
	setScrollTop: function(inTop){
		var top = inTop;
		for(var i=0, v; v=this.views[i]; i++){
			top = v.setScrollTop(inTop);
			// Work around IE not firing scroll events that cause header offset
			// issues to occur.
			if(has('ie') && v.headerNode && v.scrollboxNode){
				v.headerNode.scrollLeft = v.scrollboxNode.scrollLeft;
			}
		}
		return top;
		//this.onEach("setScrollTop", [ inTop ]);
	},
	
	getFirstScrollingView: function(){
		// summary:
		//		Returns the first grid view with a scroll bar
		for(var i=0, v; (v=this.views[i]); i++){
			if(v.hasHScrollbar() || v.hasVScrollbar()){
				return v;
			}
		}
		return null;
	}
});
});
},
'dojox/grid/_RowManager':function(){
define("dojox/grid/_RowManager", [
	"dojo/_base/declare",
	"dojo/_base/lang",
	"dojo/dom-class"
], function(declare, lang, domClass){

	var setStyleText = function(inNode, inStyleText){
		if(inNode.style.cssText == undefined){
			inNode.setAttribute("style", inStyleText);
		}else{
			inNode.style.cssText = inStyleText;
		}
	};

	return declare("dojox.grid._RowManager", null, {
		//	Stores information about grid rows. Owned by grid and used internally.
		constructor: function(inGrid){
			this.grid = inGrid;
		},
		linesToEms: 2,
		overRow: -2,
		// styles
		prepareStylingRow: function(inRowIndex, inRowNode){
			return {
				index: inRowIndex,
				node: inRowNode,
				odd: Boolean(inRowIndex&1),
				selected: !!this.grid.selection.isSelected(inRowIndex),
				over: this.isOver(inRowIndex),
				customStyles: "",
				customClasses: "dojoxGridRow"
			};
		},
		styleRowNode: function(inRowIndex, inRowNode){
			var row = this.prepareStylingRow(inRowIndex, inRowNode);
			this.grid.onStyleRow(row);
			this.applyStyles(row);
		},
		applyStyles: function(inRow){
			var i = inRow;

			i.node.className = i.customClasses;
			var h = i.node.style.height;
			setStyleText(i.node, i.customStyles + ';' + (i.node._style||''));
			i.node.style.height = h;
		},
		updateStyles: function(inRowIndex){
			this.grid.updateRowStyles(inRowIndex);
		},
		// states and events
		setOverRow: function(inRowIndex){
			var last = this.overRow;
			this.overRow = inRowIndex;
			if((last!=this.overRow)&&(lang.isString(last) || last >= 0)){
				this.updateStyles(last);
			}
			this.updateStyles(this.overRow);
		},
		isOver: function(inRowIndex){
			return (this.overRow == inRowIndex && !domClass.contains(this.grid.domNode, "dojoxGridColumnResizing"));
		}
	});
});
},
'dojox/grid/_FocusManager':function(){
define("dojox/grid/_FocusManager", [
	"dojo/_base/array",
	"dojo/_base/lang",
	"dojo/_base/declare",
	"dojo/_base/connect",
	"dojo/_base/event",
	"dojo/_base/sniff",
	"dojo/query",
	"./util",
	"dojo/_base/html"
], function(array, lang, declare, connect, event, has, query, util, html){

// focus management
return declare("dojox.grid._FocusManager", null, {
	// summary:
	//		Controls grid cell focus. Owned by grid and used internally for focusing.
	//		Note: grid cell actually receives keyboard input only when cell is being edited.
	constructor: function(inGrid){
		this.grid = inGrid;
		this.cell = null;
		this.rowIndex = -1;
		this._connects = [];
		this._headerConnects = [];
		this.headerMenu = this.grid.headerMenu;
		this._connects.push(connect.connect(this.grid.domNode, "onfocus", this, "doFocus"));
		this._connects.push(connect.connect(this.grid.domNode, "onblur", this, "doBlur"));
		this._connects.push(connect.connect(this.grid.domNode, "mousedown", this, "_mouseDown"));
		this._connects.push(connect.connect(this.grid.domNode, "mouseup", this, "_mouseUp"));
		this._connects.push(connect.connect(this.grid.domNode, "oncontextmenu", this, "doContextMenu"));
		this._connects.push(connect.connect(this.grid.lastFocusNode, "onfocus", this, "doLastNodeFocus"));
		this._connects.push(connect.connect(this.grid.lastFocusNode, "onblur", this, "doLastNodeBlur"));
		this._connects.push(connect.connect(this.grid,"_onFetchComplete", this, "_delayedCellFocus"));
		this._connects.push(connect.connect(this.grid,"postrender", this, "_delayedHeaderFocus"));
	},
	destroy: function(){
		array.forEach(this._connects, connect.disconnect);
		array.forEach(this._headerConnects, connect.disconnect);
		delete this.grid;
		delete this.cell;
	},
	_colHeadNode: null,
	_colHeadFocusIdx: null,
	_contextMenuBindNode: null,
	tabbingOut: false,
	focusClass: "dojoxGridCellFocus",
	focusView: null,
	initFocusView: function(){
		this.focusView = this.grid.views.getFirstScrollingView() || this.focusView || this.grid.views.views[0];
		this._initColumnHeaders();
	},
	isFocusCell: function(inCell, inRowIndex){
		// summary:
		//		states if the given cell is focused
		// inCell: object
		//		grid cell object
		// inRowIndex: int
		//		grid row index
		// returns:
		//		true of the given grid cell is focused
		return (this.cell == inCell) && (this.rowIndex == inRowIndex);
	},
	isLastFocusCell: function(){
		if(this.cell){
			return (this.rowIndex == this.grid.rowCount-1) && (this.cell.index == this.grid.layout.cellCount-1);
		}
		return false;
	},
	isFirstFocusCell: function(){
		if(this.cell){
			return (this.rowIndex === 0) && (this.cell.index === 0);
		}
		return false;
	},
	isNoFocusCell: function(){
		return (this.rowIndex < 0) || !this.cell;
	},
	isNavHeader: function(){
		// summary:
		//		states whether currently navigating among column headers.
		// returns:
		//		true if focus is on a column header; false otherwise.
		return (!!this._colHeadNode);
	},
	getHeaderIndex: function(){
		// summary:
		//		if one of the column headers currently has focus, return its index.
		// returns:
		//		index of the focused column header, or -1 if none have focus.
		if(this._colHeadNode){
			return array.indexOf(this._findHeaderCells(), this._colHeadNode);
		}else{
			return -1;
		}
	},
	_focusifyCellNode: function(inBork){
		var n = this.cell && this.cell.getNode(this.rowIndex);
		if(n){
			html.toggleClass(n, this.focusClass, inBork);
			if(inBork){
				var sl = this.scrollIntoView();
				try{
					if(has("webkit") || !this.grid.edit.isEditing()){
						util.fire(n, "focus");
						if(sl){ this.cell.view.scrollboxNode.scrollLeft = sl; }
					}
				}catch(e){}
			}
		}
	},
	_delayedCellFocus: function(){
		if(this.isNavHeader()||!this.grid.focused){
				return;
		}
		var n = this.cell && this.cell.getNode(this.rowIndex);
		if(n){
			try{
				if(!this.grid.edit.isEditing()){
					html.toggleClass(n, this.focusClass, true);
					if(this._colHeadNode){
						this.blurHeader();
					}
					util.fire(n, "focus");
				}
			}
			catch(e){}
		}
	},
	_delayedHeaderFocus: function(){
		if(this.isNavHeader()){
			this.focusHeader();
			//this.grid.domNode.focus();
		}
	},
	_initColumnHeaders: function(){
		array.forEach(this._headerConnects, connect.disconnect);
		this._headerConnects = [];
		var headers = this._findHeaderCells();
		for(var i = 0; i < headers.length; i++){
			this._headerConnects.push(connect.connect(headers[i], "onfocus", this, "doColHeaderFocus"));
			this._headerConnects.push(connect.connect(headers[i], "onblur", this, "doColHeaderBlur"));
		}
	},
	_findHeaderCells: function(){
		// This should be a one liner:
		//	query("th[tabindex=-1]", this.grid.viewsHeaderNode);
		// But there is a bug in query() for IE -- see trac #7037.
		var allHeads = query("th", this.grid.viewsHeaderNode);
		var headers = [];
		for (var i = 0; i < allHeads.length; i++){
			var aHead = allHeads[i];
			var hasTabIdx = html.hasAttr(aHead, "tabIndex");
			var tabindex = html.attr(aHead, "tabIndex");
			if (hasTabIdx && tabindex < 0) {
				headers.push(aHead);
			}
		}
		return headers;
	},
	_setActiveColHeader: function(/*Node*/colHeaderNode, /*Integer*/colFocusIdx, /*Integer*/ prevColFocusIdx){
		//console.log("setActiveColHeader() - colHeaderNode:colFocusIdx:prevColFocusIdx = " + colHeaderNode + ":" + colFocusIdx + ":" + prevColFocusIdx);
		this.grid.domNode.setAttribute("aria-activedescendant",colHeaderNode.id);
		if (prevColFocusIdx != null && prevColFocusIdx >= 0 && prevColFocusIdx != colFocusIdx){
			html.toggleClass(this._findHeaderCells()[prevColFocusIdx],this.focusClass,false);
		}
		html.toggleClass(colHeaderNode,this.focusClass, true);
		this._colHeadNode = colHeaderNode;
		this._colHeadFocusIdx = colFocusIdx;
		this._scrollHeader(this._colHeadFocusIdx);
	},
	scrollIntoView: function(){
		var info = (this.cell ? this._scrollInfo(this.cell) : null);
		if(!info || !info.s){
			return null;
		}
		var rt = this.grid.scroller.findScrollTop(this.rowIndex);
		// place cell within horizontal view
		if(info.n && info.sr){
			if(info.n.offsetLeft + info.n.offsetWidth > info.sr.l + info.sr.w){
				info.s.scrollLeft = info.n.offsetLeft + info.n.offsetWidth - info.sr.w;
			}else if(info.n.offsetLeft < info.sr.l){
				info.s.scrollLeft = info.n.offsetLeft;
			}
		}
		// place cell within vertical view
		if(info.r && info.sr){
			if(rt + info.r.offsetHeight > info.sr.t + info.sr.h){
				this.grid.setScrollTop(rt + info.r.offsetHeight - info.sr.h);
			}else if(rt < info.sr.t){
				this.grid.setScrollTop(rt);
			}
		}

		return info.s.scrollLeft;
	},
	_scrollInfo: function(cell, domNode){
		if(cell){
			var cl = cell,
				sbn = cl.view.scrollboxNode,
				sbnr = {
					w: sbn.clientWidth,
					l: sbn.scrollLeft,
					t: sbn.scrollTop,
					h: sbn.clientHeight
				},
				rn = cl.view.getRowNode(this.rowIndex);
			return {
				c: cl,
				s: sbn,
				sr: sbnr,
				n: (domNode ? domNode : cell.getNode(this.rowIndex)),
				r: rn
			};
		}
		return null;
	},
	_scrollHeader: function(currentIdx){
		var info = null;
		if(this._colHeadNode){
			var cell = this.grid.getCell(currentIdx);
			if(!cell){ return; }
			info = this._scrollInfo(cell, cell.getNode(0));
		}
		if(info && info.s && info.sr && info.n){
			// scroll horizontally as needed.
			var scroll = info.sr.l + info.sr.w;
			if(info.n.offsetLeft + info.n.offsetWidth > scroll){
				info.s.scrollLeft = info.n.offsetLeft + info.n.offsetWidth - info.sr.w;
			}else if(info.n.offsetLeft < info.sr.l){
				info.s.scrollLeft = info.n.offsetLeft;
			}else if(has('ie') <= 7 && cell && cell.view.headerNode){
				// Trac 7158: scroll dojoxGridHeader for IE7 and lower
				cell.view.headerNode.scrollLeft = info.s.scrollLeft;
			}
		}
	},
	_isHeaderHidden: function(){
		// summary:
		//		determine if the grid headers are hidden
		//		relies on documented technique of setting .dojoxGridHeader { display:none; }
		// returns: Boolean
		//		true if headers are hidden
		//		false if headers are not hidden
		
		var curView = this.focusView;
		if (!curView){
			// find one so we can determine if headers are hidden
			// there is no focusView after adding items to empty grid (test_data_grid_empty.html)
			for (var i = 0, cView; (cView = this.grid.views.views[i]); i++) {
				if(cView.headerNode ){
					curView=cView;
					break;
				}
			}
		}
		return (curView && html.getComputedStyle(curView.headerNode).display == "none");
	},
	colSizeAdjust: function (e, colIdx, delta){ // adjust the column specified by colIdx by the specified delta px
		var headers = this._findHeaderCells();
		var view = this.focusView;
		if(!view || !view.header.tableMap.map){
			for(var i = 0, cView; (cView = this.grid.views.views[i]); i++){
				// find first view with a tableMap in order to work with empty grid
				if(cView.header.tableMap.map){
					view=cView;
					break;
				}
			}
		}
		var curHeader = headers[colIdx];
		if (!view || (colIdx == headers.length-1 && colIdx === 0)){
			return; // can't adjust single col. grid
		}
		view.content.baseDecorateEvent(e);
		// need to adjust event with header cell info since focus is no longer on header cell
		e.cellNode = curHeader; //this.findCellTarget(e.target, e.rowNode);
		e.cellIndex = view.content.getCellNodeIndex(e.cellNode);
		e.cell = (e.cellIndex >= 0 ? this.grid.getCell(e.cellIndex) : null);
		if (view.header.canResize(e)){
			var deltaObj = {
				l: delta
			};
			var drag = view.header.colResizeSetup(e,false);
			view.header.doResizeColumn(drag, null, deltaObj);
			view.update();
		}
	},
	styleRow: function(inRow){
		return;
	},
	setFocusIndex: function(inRowIndex, inCellIndex){
		// summary:
		//		focuses the given grid cell
		// inRowIndex: int
		//		grid row index
		// inCellIndex: int
		//		grid cell index
		this.setFocusCell(this.grid.getCell(inCellIndex), inRowIndex);
	},
	setFocusCell: function(inCell, inRowIndex){
		// summary:
		//		focuses the given grid cell
		// inCell: object
		//		grid cell object
		// inRowIndex: int
		//		grid row index
		if(inCell && !this.isFocusCell(inCell, inRowIndex)){
			this.tabbingOut = false;
			if (this._colHeadNode){
				this.blurHeader();
			}
			this._colHeadNode = this._colHeadFocusIdx = null;
			this.focusGridView();
			this._focusifyCellNode(false);
			this.cell = inCell;
			this.rowIndex = inRowIndex;
			this._focusifyCellNode(true);
		}
		// even if this cell isFocusCell, the document focus may need to be rejiggered
		// call opera on delay to prevent keypress from altering focus
		if(has('opera')){
			setTimeout(lang.hitch(this.grid, 'onCellFocus', this.cell, this.rowIndex), 1);
		}else{
			this.grid.onCellFocus(this.cell, this.rowIndex);
		}
	},
	next: function(){
		// summary:
		//	focus next grid cell
		if(this.cell){
			var row=this.rowIndex, col=this.cell.index+1, cc=this.grid.layout.cellCount-1, rc=this.grid.rowCount-1;
			if(col > cc){
				col = 0;
				row++;
			}
			if(row > rc){
				col = cc;
				row = rc;
			}
			if(this.grid.edit.isEditing()){ //when editing, only navigate to editable cells
				var nextCell = this.grid.getCell(col);
				if (!this.isLastFocusCell() && (!nextCell.editable ||
					this.grid.canEdit && !this.grid.canEdit(nextCell, row))){
					this.cell=nextCell;
					this.rowIndex=row;
					this.next();
					return;
				}
			}
			this.setFocusIndex(row, col);
		}
	},
	previous: function(){
		// summary:
		//	focus previous grid cell
		if(this.cell){
			var row=(this.rowIndex || 0), col=(this.cell.index || 0) - 1;
			if(col < 0){
				col = this.grid.layout.cellCount-1;
				row--;
			}
			if(row < 0){
				row = 0;
				col = 0;
			}
			if(this.grid.edit.isEditing()){ //when editing, only navigate to editable cells
				var prevCell = this.grid.getCell(col);
				if (!this.isFirstFocusCell() && !prevCell.editable){
					this.cell=prevCell;
					this.rowIndex=row;
					this.previous();
					return;
				}
			}
			this.setFocusIndex(row, col);
		}
	},
	move: function(inRowDelta, inColDelta) {
		// summary:
		//		focus grid cell or  simulate focus to column header based on position relative to current focus
		// inRowDelta: int
		//		vertical distance from current focus
		// inColDelta: int
		//		horizontal distance from current focus

		var colDir = inColDelta < 0 ? -1 : 1;
		// Handle column headers.
		if(this.isNavHeader()){
			var headers = this._findHeaderCells();
			var savedIdx = currentIdx = array.indexOf(headers, this._colHeadNode);
			currentIdx += inColDelta;
			while(currentIdx >=0 && currentIdx < headers.length && headers[currentIdx].style.display == "none"){
				// skip over hidden column headers
				currentIdx += colDir;
			}
			if((currentIdx >= 0) && (currentIdx < headers.length)){
				this._setActiveColHeader(headers[currentIdx],currentIdx, savedIdx);
			}
		}else{
			if(this.cell){
				// Handle grid proper.
				var sc = this.grid.scroller,
					r = this.rowIndex,
					rc = this.grid.rowCount-1,
					row = Math.min(rc, Math.max(0, r+inRowDelta));
				if(inRowDelta){
					if(inRowDelta>0){
						if(row > sc.getLastPageRow(sc.page)){
							//need to load additional data, let scroller do that
							this.grid.setScrollTop(this.grid.scrollTop+sc.findScrollTop(row)-sc.findScrollTop(r));
						}
					}else if(inRowDelta<0){
						if(row <= sc.getPageRow(sc.page)){
							//need to load additional data, let scroller do that
							this.grid.setScrollTop(this.grid.scrollTop-sc.findScrollTop(r)-sc.findScrollTop(row));
						}
					}
				}
				var cc = this.grid.layout.cellCount-1,
				i = this.cell.index,
				col = Math.min(cc, Math.max(0, i+inColDelta));
				var cell = this.grid.getCell(col);
				while(col>=0 && col < cc && cell && cell.hidden === true){
					// skip hidden cells
					col += colDir;
					cell = this.grid.getCell(col);
				}
				if (!cell || cell.hidden === true){
					// don't change col if would move to hidden
					col = i;
				}
				//skip hidden row|cell
				var n = cell.getNode(row);
				if(!n && inRowDelta){
					if((row + inRowDelta) >= 0 && (row + inRowDelta) <= rc){
						this.move(inRowDelta > 0 ? ++inRowDelta : --inRowDelta, inColDelta);
					}
					return;
				}else if((!n || html.style(n, "display") === "none") && inColDelta){
					if((col + inColDelta) >= 0 && (col + inColDelta) <= cc){
						this.move(inRowDelta, inColDelta > 0 ? ++inColDelta : --inColDelta);
					}
					return;
				}
				this.setFocusIndex(row, col);
				if(inRowDelta){
					this.grid.updateRow(r);
				}
			}
		}
	},
	previousKey: function(e){
		if(this.grid.edit.isEditing()){
			event.stop(e);
			this.previous();
		}else if(!this.isNavHeader() && !this._isHeaderHidden()) {
			this.grid.domNode.focus(); // will call doFocus and set focus into header.
			event.stop(e);
		}else{
			this.tabOut(this.grid.domNode);
			if (this._colHeadFocusIdx != null) { // clear grid header focus
				html.toggleClass(this._findHeaderCells()[this._colHeadFocusIdx], this.focusClass, false);
				this._colHeadFocusIdx = null;
			}
			this._focusifyCellNode(false);
		}
	},
	nextKey: function(e) {
		var isEmpty = (this.grid.rowCount === 0);
		if(e.target === this.grid.domNode && this._colHeadFocusIdx == null){
			this.focusHeader();
			event.stop(e);
		}else if(this.isNavHeader()){
			// if tabbing from col header, then go to grid proper.
			this.blurHeader();
			if(!this.findAndFocusGridCell()){
				this.tabOut(this.grid.lastFocusNode);
			}
			this._colHeadNode = this._colHeadFocusIdx= null;
		}else if(this.grid.edit.isEditing()){
			event.stop(e);
			this.next();
		}else{
			this.tabOut(this.grid.lastFocusNode);
		}
	},
	tabOut: function(inFocusNode){
		this.tabbingOut = true;
		inFocusNode.focus();
	},
	focusGridView: function(){
		util.fire(this.focusView, "focus");
	},
	focusGrid: function(inSkipFocusCell){
		this.focusGridView();
		this._focusifyCellNode(true);
	},
	findAndFocusGridCell: function(){
		// summary:
		//		find the first focusable grid cell
		// returns: Boolean
		//		true if focus was set to a cell
		//		false if no cell found to set focus onto
		
		var didFocus = true;
		var isEmpty = (this.grid.rowCount === 0); // If grid is empty this.grid.rowCount == 0
		if (this.isNoFocusCell() && !isEmpty){
			var cellIdx = 0;
			var cell = this.grid.getCell(cellIdx);
			if (cell.hidden) {
				// if first cell isn't visible, use _colHeadFocusIdx
				// could also use a while loop to find first visible cell - not sure that is worth it
				cellIdx = this.isNavHeader() ? this._colHeadFocusIdx : 0;
			}
			this.setFocusIndex(0, cellIdx);
		}
		else if (this.cell && !isEmpty){
			if (this.focusView && !this.focusView.rowNodes[this.rowIndex]){
				// if rowNode for current index is undefined (likely as a result of a sort and because of #7304)
				// scroll to that row
				this.grid.scrollToRow(this.rowIndex);
			}
			this.focusGrid();
		}else {
			didFocus = false;
		}
		this._colHeadNode = this._colHeadFocusIdx= null;
		return didFocus;
	},
	focusHeader: function(){
		var headerNodes = this._findHeaderCells();
		var saveColHeadFocusIdx = this._colHeadFocusIdx;
		if (this._isHeaderHidden()){
			// grid header is hidden, focus a cell
			this.findAndFocusGridCell();
		}
		else if (!this._colHeadFocusIdx) {
			if (this.isNoFocusCell()) {
				this._colHeadFocusIdx = 0;
			}
			else {
				this._colHeadFocusIdx = this.cell.index;
			}
		}
		this._colHeadNode = headerNodes[this._colHeadFocusIdx];
		while(this._colHeadNode && this._colHeadFocusIdx >=0 && this._colHeadFocusIdx < headerNodes.length &&
				this._colHeadNode.style.display == "none"){
			// skip over hidden column headers
			this._colHeadFocusIdx++;
			this._colHeadNode = headerNodes[this._colHeadFocusIdx];
		}
		if(this._colHeadNode && this._colHeadNode.style.display != "none"){
			// Column header cells know longer receive actual focus.  So, for keyboard invocation of
			// contextMenu to work, the contextMenu must be bound to the grid.domNode rather than the viewsHeaderNode.
			// unbind the contextmenu from the viewsHeaderNode and to the grid when header cells are active.  Reset
			// the binding back to the viewsHeaderNode when header cells are no longer acive (in blurHeader) #10483
			if (this.headerMenu && this._contextMenuBindNode != this.grid.domNode){
				this.headerMenu.unBindDomNode(this.grid.viewsHeaderNode);
				this.headerMenu.bindDomNode(this.grid.domNode);
				this._contextMenuBindNode = this.grid.domNode;
			}
			this._setActiveColHeader(this._colHeadNode, this._colHeadFocusIdx, saveColHeadFocusIdx);
			this._scrollHeader(this._colHeadFocusIdx);
			this._focusifyCellNode(false);
		}else {
			// all col head nodes are hidden - focus the grid
			this.findAndFocusGridCell();
		}
	},
	blurHeader: function(){
		html.removeClass(this._colHeadNode, this.focusClass);
		html.removeAttr(this.grid.domNode,"aria-activedescendant");
		// reset contextMenu onto viewsHeaderNode so right mouse on header will invoke (see focusHeader)
		if (this.headerMenu && this._contextMenuBindNode == this.grid.domNode) {
			var viewsHeader = this.grid.viewsHeaderNode;
			this.headerMenu.unBindDomNode(this.grid.domNode);
			this.headerMenu.bindDomNode(viewsHeader);
			this._contextMenuBindNode = viewsHeader;
		}
	},
	doFocus: function(e){
		// trap focus only for grid dom node
		if(e && e.target != e.currentTarget){
			event.stop(e);
			return;
		}
		// don't change focus if clicking on scroller bar
		if(this._clickFocus){
			return;
		}
		// do not focus for scrolling if grid is about to blur
		if(!this.tabbingOut){
			this.focusHeader();
		}
		this.tabbingOut = false;
		event.stop(e);
	},
	doBlur: function(e){
		event.stop(e);	// FF2
	},
	doContextMenu: function(e){
	//stop contextMenu event if no header Menu to prevent default/browser contextMenu
		if (!this.headerMenu){
			event.stop(e);
		}
	},
	doLastNodeFocus: function(e){
		if (this.tabbingOut){
			this._focusifyCellNode(false);
		}else if(this.grid.rowCount >0){
			if (this.isNoFocusCell()){
				this.setFocusIndex(0,0);
			}
			this._focusifyCellNode(true);
		}else {
			this.focusHeader();
		}
		this.tabbingOut = false;
		event.stop(e);	 // FF2
	},
	doLastNodeBlur: function(e){
		event.stop(e);	 // FF2
	},
	doColHeaderFocus: function(e){
		this._setActiveColHeader(e.target,html.attr(e.target, "idx"),this._colHeadFocusIdx);
		this._scrollHeader(this.getHeaderIndex());
		event.stop(e);
	},
	doColHeaderBlur: function(e){
		html.toggleClass(e.target, this.focusClass, false);
	},
	_mouseDown: function(e){
		// a flag indicating grid is being focused by clicking
		this._clickFocus = dojo.some(this.grid.views.views, function(v){
			return v.scrollboxNode === e.target;
		});
	},
	_mouseUp: function(e){
		this._clickFocus = false;
	}
});
});
},
'dojox/grid/_EditManager':function(){
define("dojox/grid/_EditManager", [
	"dojo/_base/lang",
	"dojo/_base/array",
	"dojo/_base/declare",
	"dojo/_base/connect",
	"dojo/_base/sniff",
	"./util"
], function(lang, array, declare, connect, has, util){

return declare("dojox.grid._EditManager", null, {
	// summary:
	//		Controls grid cell editing process. Owned by grid and used internally for editing.
	constructor: function(inGrid){
		// inGrid: dojox.Grid
		//		The dojox.Grid this editor should be attached to
		this.grid = inGrid;
		if(has('ie')){
			this.connections = [connect.connect(document.body, "onfocus", lang.hitch(this, "_boomerangFocus"))];
		}else{
			this.connections = [connect.connect(this.grid, 'onBlur', this, 'apply')];
		}
	},
	
	info: {},

	destroy: function(){
		array.forEach(this.connections, connect.disconnect);
	},

	cellFocus: function(inCell, inRowIndex){
		// summary:
		//		Invoke editing when cell is focused
		// inCell: cell object
		//		Grid cell object
		// inRowIndex: Integer
		//		Grid row index
		if(this.grid.singleClickEdit || this.isEditRow(inRowIndex)){
			// if same row or quick editing, edit
			this.setEditCell(inCell, inRowIndex);
		}else{
			// otherwise, apply any pending row edits
			this.apply();
		}
		// if dynamic or static editing...
		if(this.isEditing() || (inCell && inCell.editable && inCell.alwaysEditing)){
			// let the editor focus itself as needed
			this._focusEditor(inCell, inRowIndex);
		}
	},

	rowClick: function(e){
		if(this.isEditing() && !this.isEditRow(e.rowIndex)){
			this.apply();
		}
	},

	styleRow: function(inRow){
		if(inRow.index == this.info.rowIndex){
			inRow.customClasses += ' dojoxGridRowEditing';
		}
	},

	dispatchEvent: function(e){
		var c = e.cell, ed = (c && c["editable"]) ? c : 0;
		return ed && ed.dispatchEvent(e.dispatch, e);
	},

	// Editing
	isEditing: function(){
		// summary:
		//		Indicates editing state of the grid.
		// returns: Boolean
		//	 	True if grid is actively editing
		return this.info.rowIndex !== undefined;
	},

	isEditCell: function(inRowIndex, inCellIndex){
		// summary:
		//		Indicates if the given cell is being edited.
		// inRowIndex: Integer
		//		Grid row index
		// inCellIndex: Integer
		//		Grid cell index
		// returns: Boolean
		//	 	True if given cell is being edited
		return (this.info.rowIndex === inRowIndex) && (this.info.cell.index == inCellIndex);
	},

	isEditRow: function(inRowIndex){
		// summary:
		//		Indicates if the given row is being edited.
		// inRowIndex: Integer
		//		Grid row index
		// returns: Boolean
		//	 	True if given row is being edited
		return this.info.rowIndex === inRowIndex;
	},

	setEditCell: function(inCell, inRowIndex){
		// summary:
		//		Set the given cell to be edited
		// inRowIndex: Integer
		//		Grid row index
		// inCell: Object
		//		Grid cell object
		if(!this.isEditCell(inRowIndex, inCell.index) && this.grid.canEdit && this.grid.canEdit(inCell, inRowIndex)){
			this.start(inCell, inRowIndex, this.isEditRow(inRowIndex) || inCell.editable);
		}
	},

	_focusEditor: function(inCell, inRowIndex){
		util.fire(inCell, "focus", [inRowIndex]);
	},

	focusEditor: function(){
		if(this.isEditing()){
			this._focusEditor(this.info.cell, this.info.rowIndex);
		}
	},

	// implement fix for focus boomerang effect on IE
	_boomerangWindow: 500,
	_shouldCatchBoomerang: function(){
		return this._catchBoomerang > new Date().getTime();
	},
	_boomerangFocus: function(){
		//console.log("_boomerangFocus");
		if(this._shouldCatchBoomerang()){
			// make sure we don't utterly lose focus
			this.grid.focus.focusGrid();
			// let the editor focus itself as needed
			this.focusEditor();
			// only catch once
			this._catchBoomerang = 0;
		}
	},
	_doCatchBoomerang: function(){
		// give ourselves a few ms to boomerang IE focus effects
		if(has('ie')){this._catchBoomerang = new Date().getTime() + this._boomerangWindow;}
	},
	// end boomerang fix API

	start: function(inCell, inRowIndex, inEditing){
		if(!this._isValidInput()){
			return;
		}
		this.grid.beginUpdate();
		this.editorApply();
		if(this.isEditing() && !this.isEditRow(inRowIndex)){
			this.applyRowEdit();
			this.grid.updateRow(inRowIndex);
		}
		if(inEditing){
			this.info = { cell: inCell, rowIndex: inRowIndex };
			this.grid.doStartEdit(inCell, inRowIndex);
			this.grid.updateRow(inRowIndex);
		}else{
			this.info = {};
		}
		this.grid.endUpdate();
		// make sure we don't utterly lose focus
		this.grid.focus.focusGrid();
		// let the editor focus itself as needed
		this._focusEditor(inCell, inRowIndex);
		// give ourselves a few ms to boomerang IE focus effects
		this._doCatchBoomerang();
	},

	_editorDo: function(inMethod){
		var c = this.info.cell;
		//c && c.editor && c.editor[inMethod](c, this.info.rowIndex);
		if(c && c.editable){
			c[inMethod](this.info.rowIndex);
		}
	},

	editorApply: function(){
		this._editorDo("apply");
	},

	editorCancel: function(){
		this._editorDo("cancel");
	},

	applyCellEdit: function(inValue, inCell, inRowIndex){
		if(this.grid.canEdit(inCell, inRowIndex)){
			this.grid.doApplyCellEdit(inValue, inRowIndex, inCell.field);
		}
	},

	applyRowEdit: function(){
		this.grid.doApplyEdit(this.info.rowIndex, this.info.cell.field);
	},

	apply: function(){
		// summary:
		//		Apply a grid edit
		if(this.isEditing() && this._isValidInput()){
			this.grid.beginUpdate();
			this.editorApply();
			this.applyRowEdit();
			this.info = {};
			this.grid.endUpdate();
			this.grid.focus.focusGrid();
			this._doCatchBoomerang();
		}
	},

	cancel: function(){
		// summary:
		//		Cancel a grid edit
		if(this.isEditing()){
			this.grid.beginUpdate();
			this.editorCancel();
			this.info = {};
			this.grid.endUpdate();
			this.grid.focus.focusGrid();
			this._doCatchBoomerang();
		}
	},

	save: function(inRowIndex, inView){
		// summary:
		//		Save the grid editing state
		// inRowIndex: Integer
		//		Grid row index
		// inView: Object
		//		Grid view
		var c = this.info.cell;
		if(this.isEditRow(inRowIndex) && (!inView || c.view==inView) && c.editable){
			c.save(c, this.info.rowIndex);
		}
	},

	restore: function(inView, inRowIndex){
		// summary:
		//		Restores the grid editing state
		// inRowIndex: Integer
		//		Grid row index
		// inView: Object
		//		Grid view
		var c = this.info.cell;
		if(this.isEditRow(inRowIndex) && c.view == inView && c.editable){
			c.restore(this.info.rowIndex);
		}
	},
	
	_isValidInput: function(){
		var w = (this.info.cell || {}).widget;		
		if(!w || !w.isValid){
			//no validation needed
			return true;
		}		
		w.focused = true;
		return w.isValid(true);
	}
});
});
},
'dojox/grid/Selection':function(){
define("dojox/grid/Selection", [
	"dojo/_base/declare",
	"dojo/_base/array",
	"dojo/_base/lang",
	"dojo/dom-attr"
], function(declare, array, lang, domAttr){

return declare("dojox.grid.Selection", null, {
	// summary:
	//		Manages row selection for grid. Owned by grid and used internally
	//		for selection. Override to implement custom selection.

	constructor: function(inGrid){
		this.grid = inGrid;
		this.selected = [];

		this.setMode(inGrid.selectionMode);
	},

	mode: 'extended',

	selected: null,
	updating: 0,
	selectedIndex: -1,
	rangeStartIndex: -1,

	setMode: function(mode){
		if(this.selected.length){
			this.deselectAll();
		}
		if(mode != 'extended' && mode != 'multiple' && mode != 'single' && mode != 'none'){
			this.mode = 'extended';
		}else{
			this.mode = mode;
		}
	},

	onCanSelect: function(inIndex){
		return this.grid.onCanSelect(inIndex);
	},

	onCanDeselect: function(inIndex){
		return this.grid.onCanDeselect(inIndex);
	},

	onSelected: function(inIndex){
	},

	onDeselected: function(inIndex){
	},

	//onSetSelected: function(inIndex, inSelect) { };
	onChanging: function(){
	},

	onChanged: function(){
	},

	isSelected: function(inIndex){
		if(this.mode == 'none'){
			return false;
		}
		return this.selected[inIndex];
	},

	getFirstSelected: function(){
		if(!this.selected.length||this.mode == 'none'){ return -1; }
		for(var i=0, l=this.selected.length; i<l; i++){
			if(this.selected[i]){
				return i;
			}
		}
		return -1;
	},

	getNextSelected: function(inPrev){
		if(this.mode == 'none'){ return -1; }
		for(var i=inPrev+1, l=this.selected.length; i<l; i++){
			if(this.selected[i]){
				return i;
			}
		}
		return -1;
	},

	getSelected: function(){
		var result = [];
		for(var i=0, l=this.selected.length; i<l; i++){
			if(this.selected[i]){
				result.push(i);
			}
		}
		return result;
	},

	getSelectedCount: function(){
		var c = 0;
		for(var i=0; i<this.selected.length; i++){
			if(this.selected[i]){
				c++;
			}
		}
		return c;
	},

	_beginUpdate: function(){
		if(this.updating === 0){
			this.onChanging();
		}
		this.updating++;
	},

	_endUpdate: function(){
		this.updating--;
		if(this.updating === 0){
			this.onChanged();
		}
	},

	select: function(inIndex){
		if(this.mode == 'none'){ return; }
		if(this.mode != 'multiple'){
			this.deselectAll(inIndex);
			this.addToSelection(inIndex);
		}else{
			this.toggleSelect(inIndex);
		}
	},

	addToSelection: function(inIndex){
		if(this.mode == 'none'){ return; }
		if(lang.isArray(inIndex)){
			array.forEach(inIndex, this.addToSelection, this);
			return;
		}
		inIndex = Number(inIndex);
		if(this.selected[inIndex]){
			this.selectedIndex = inIndex;
		}else{
			if(this.onCanSelect(inIndex) !== false){
				this.selectedIndex = inIndex;
				var rowNode = this.grid.getRowNode(inIndex);
				if(rowNode){
					domAttr.set(rowNode, "aria-selected", "true");
				}
				this._beginUpdate();
				this.selected[inIndex] = true;
				//this.grid.onSelected(inIndex);
				this.onSelected(inIndex);
				//this.onSetSelected(inIndex, true);
				this._endUpdate();
			}
		}
	},

	deselect: function(inIndex){
		if(this.mode == 'none'){ return; }
		if(lang.isArray(inIndex)){
			array.forEach(inIndex, this.deselect, this);
			return;
		}
		inIndex = Number(inIndex);
		if(this.selectedIndex == inIndex){
			this.selectedIndex = -1;
		}
		if(this.selected[inIndex]){
			if(this.onCanDeselect(inIndex) === false){
				return;
			}
			var rowNode = this.grid.getRowNode(inIndex);
			if(rowNode){
				domAttr.set(rowNode, "aria-selected", "false");
			}
			this._beginUpdate();
			delete this.selected[inIndex];
			//this.grid.onDeselected(inIndex);
			this.onDeselected(inIndex);
			//this.onSetSelected(inIndex, false);
			this._endUpdate();
		}
	},

	setSelected: function(inIndex, inSelect){
		this[(inSelect ? 'addToSelection' : 'deselect')](inIndex);
	},

	toggleSelect: function(inIndex){
		if(lang.isArray(inIndex)){
			array.forEach(inIndex, this.toggleSelect, this);
			return;
		}
		this.setSelected(inIndex, !this.selected[inIndex]);
	},

	_range: function(inFrom, inTo, func){
		var s = (inFrom >= 0 ? inFrom : inTo), e = inTo;
		if(s > e){
			e = s;
			s = inTo;
		}
		for(var i=s; i<=e; i++){
			func(i);
		}
	},

	selectRange: function(inFrom, inTo){
		this._range(inFrom, inTo, lang.hitch(this, "addToSelection"));
	},

	deselectRange: function(inFrom, inTo){
		this._range(inFrom, inTo, lang.hitch(this, "deselect"));
	},

	insert: function(inIndex){
		this.selected.splice(inIndex, 0, false);
		if(this.selectedIndex >= inIndex){
			this.selectedIndex++;
		}
	},

	remove: function(inIndex){
		this.selected.splice(inIndex, 1);
		if(this.selectedIndex >= inIndex){
			this.selectedIndex--;
		}
	},

	deselectAll: function(inExcept){
		for(var i in this.selected){
			if((i!=inExcept)&&(this.selected[i]===true)){
				this.deselect(i);
			}
		}
	},

	clickSelect: function(inIndex, inCtrlKey, inShiftKey){
		if(this.mode == 'none'){ return; }
		this._beginUpdate();
		if(this.mode != 'extended'){
			this.select(inIndex);
		}else{
			if(!inShiftKey || this.rangeStartIndex < 0){
				this.rangeStartIndex = inIndex;
			}
			if(!inCtrlKey){
				this.deselectAll(inIndex);
			}
			if(inShiftKey){
				this.selectRange(this.rangeStartIndex, inIndex);
			}else if(inCtrlKey){
				this.toggleSelect(inIndex);
			}else{
				this.addToSelection(inIndex);
			}
		}
		this._endUpdate();
	},

	clickSelectEvent: function(e){
		this.clickSelect(e.rowIndex, dojo.isCopyKey(e), e.shiftKey);
	},

	clear: function(){
		this._beginUpdate();
		this.deselectAll();
		this._endUpdate();
	}
});
});
},
'url:dojox/grid/resources/_Grid.html':"<div hidefocus=\"hidefocus\" role=\"grid\" dojoAttachEvent=\"onmouseout:_mouseOut\">\n\t<div class=\"dojoxGridMasterHeader\" dojoAttachPoint=\"viewsHeaderNode\" role=\"presentation\"></div>\n\t<div class=\"dojoxGridMasterView\" dojoAttachPoint=\"viewsNode\" role=\"presentation\"></div>\n\t<div class=\"dojoxGridMasterMessages\" style=\"display: none;\" dojoAttachPoint=\"messagesNode\"></div>\n\t<span dojoAttachPoint=\"lastFocusNode\" tabindex=\"0\"></span>\n</div>\n",
'dojox/grid/DataSelection':function(){
define("dojox/grid/DataSelection", [
	"dojo/_base/declare",
	"./_SelectionPreserver",
	"./Selection"
], function(declare, _SelectionPreserver, Selection){
	
return declare("dojox.grid.DataSelection", Selection, {
	constructor: function(grid){
		if(grid.keepSelection){
			this.preserver = new _SelectionPreserver(this);
		}
	},
	
	destroy: function(){
		if(this.preserver){
			this.preserver.destroy();
		}
	},
	
	getFirstSelected: function(){
		var idx = Selection.prototype.getFirstSelected.call(this);

		if(idx == -1){ return null; }
		return this.grid.getItem(idx);
	},

	getNextSelected: function(inPrev){
		var old_idx = this.grid.getItemIndex(inPrev);
		var idx = Selection.prototype.getNextSelected.call(this, old_idx);

		if(idx == -1){ return null; }
		return this.grid.getItem(idx);
	},

	getSelected: function(){
		var result = [];
		for(var i=0, l=this.selected.length; i<l; i++){
			if(this.selected[i]){
				result.push(this.grid.getItem(i));
			}
		}
		return result;
	},

	addToSelection: function(inItemOrIndex){
		if(this.mode == 'none'){ return; }
		var idx = null;
		if(typeof inItemOrIndex == "number" || typeof inItemOrIndex == "string"){
			idx = inItemOrIndex;
		}else{
			idx = this.grid.getItemIndex(inItemOrIndex);
		}
		Selection.prototype.addToSelection.call(this, idx);
	},

	deselect: function(inItemOrIndex){
		if(this.mode == 'none'){ return; }
		var idx = null;
		if(typeof inItemOrIndex == "number" || typeof inItemOrIndex == "string"){
			idx = inItemOrIndex;
		}else{
			idx = this.grid.getItemIndex(inItemOrIndex);
		}
		Selection.prototype.deselect.call(this, idx);
	},

	deselectAll: function(inItemOrIndex){
		var idx = null;
		if(inItemOrIndex || typeof inItemOrIndex == "number"){
			if(typeof inItemOrIndex == "number" || typeof inItemOrIndex == "string"){
				idx = inItemOrIndex;
			}else{
				idx = this.grid.getItemIndex(inItemOrIndex);
			}
			Selection.prototype.deselectAll.call(this, idx);
		}else{
			this.inherited(arguments);
		}
	}
});
});
},
'dojox/grid/_SelectionPreserver':function(){
define("dojox/grid/_SelectionPreserver", [
	"dojo/_base/declare",
	"dojo/_base/connect",
	"dojo/_base/lang",
	"dojo/_base/array"
], function(declare, connect, lang, array){

return declare("dojox.grid._SelectionPreserver", null, {
	// summary:
	//		Preserve selections across various user actions.
	//
	// description:
	//		When this feature is turned on, Grid will try to preserve selections across actions, e.g. sorting, filtering etc.
	//
	//		Precondition - Identifier(id) is required for store since id is the only way for differentiating row items.
	//		Known issue - The preserved selections might be inaccurate if some unloaded rows are previously selected by range(e.g.SHIFT + click)
	//
	// example:
	// |	//To turn on this - please set 'keepSelection' attribute to true
	// |	<div dojoType="dojox.grid.DataGrid" keepSelection = true .../>
	// |	<div dojoType="dojox.grid.TreeGrid" keepSelection = true .../>
	// |	<div dojoType="dojox.grid.LazyTreeGrid" keepSelection = true .../>
	
	constructor: function(selection){
		this.selection = selection;
		var grid = this.grid = selection.grid;
		this.reset();
		this._connects = [
			connect.connect(grid, '_setStore', this, 'reset'),
			connect.connect(grid, '_addItem', this, '_reSelectById'),
			connect.connect(selection, 'onSelected', lang.hitch(this, '_selectById', true)),
			connect.connect(selection, 'onDeselected', lang.hitch(this, '_selectById', false)),
			connect.connect(selection, 'deselectAll', this, 'reset')
		];
	},
	destroy: function(){
		this.reset();
		array.forEach(this._connects, connect.disconnect);
		delete this._connects;
	},
	reset: function(){
		this._selectedById = {};
	},
	_reSelectById: function(item, index){
		// summary:
		//		When some rows is fetched, determine whether it should be selected.
		if(item && this.grid._hasIdentity){
			this.selection.selected[index] = this._selectedById[this.grid.store.getIdentity(item)];
		}
	},
	_selectById: function(toSelect, inItemOrIndex){
		// summary:
		//		Record selected rows by ID.
		if(this.selection.mode == 'none' || !this.grid._hasIdentity){ return; }
		var item = inItemOrIndex, g = this.grid;
		if(typeof inItemOrIndex == "number" || typeof inItemOrIndex == "string"){
			var entry = g._by_idx[inItemOrIndex];
			item = entry && entry.item;
		}
		if(item){
			this._selectedById[g.store.getIdentity(item)] = !!toSelect;
		}
		return item;
	}
});
});
},
'dojox/grid/enhanced/_PluginManager':function(){
define("dojox/grid/enhanced/_PluginManager", [
	"dojo/_base/kernel",
	"dojo/_base/lang",
	"dojo/_base/declare",
	"dojo/_base/array",
	"dojo/_base/connect",
	"./_Events",
	"./_FocusManager",
	"../util"
], function(dojo, lang, declare, array, connect, _Events, _FocusManager, util){

var _PluginManager = declare("dojox.grid.enhanced._PluginManager", null, {
	// summary:
	//		Singleton plugin manager
	// description:
	//		Plugin manager is responsible for:
	//
	//		1. Loading required plugins
	//		2. Handling collaboration and dependencies among plugins
	//
	//		Some plugin dependencies:
	//
	//		- "columnReordering" attribute won't work when either DnD or Indirect Selections plugin is on.
		
	// _options: Object
	//		Normalized plugin options
	_options: null,

	// _plugins: Array
	//		Plugin list
	_plugins: null,

	// _connects: Array
	//		Connection list
	_connects: null,

	constructor: function(inGrid){
		this.grid = inGrid;
		this._store = inGrid.store;
		this._options = {};
		this._plugins = [];
		this._connects = [];
		this._parseProps(this.grid.plugins);
		
		inGrid.connect(inGrid, "_setStore", lang.hitch(this, function(store){
			if(this._store !== store){
				this.forEach('onSetStore', [store, this._store]);
				this._store = store;
			}
		}));
	},
	startup: function(){
		this.forEach('onStartUp');
	},
	preInit: function(){
		// summary:
		//		Load appropriate plugins before DataGrid.postCreate().
		//		See EnhancedGrid.postCreate()
		this.grid.focus.destroy();
		this.grid.focus = new _FocusManager(this.grid);
		new _Events(this.grid);//overwrite some default events of DataGrid
		this._init(true);
		this.forEach('onPreInit');
	},
	postInit: function(){
		// summary:
		//		Load plugins after DataGrid.postCreate() - the default phase when plugins are created
		//		See EnhancedGrid.postCreate()
		this._init(false);
		
		array.forEach(this.grid.views.views, this._initView, this);
		this._connects.push(connect.connect(this.grid.views, 'addView', lang.hitch(this, this._initView)));
			
		if(this._plugins.length > 0){
			var edit = this.grid.edit;
			if(edit){ edit.styleRow = function(inRow){}; }
		}
		this.forEach('onPostInit');
	},
	forEach: function(func, args){
		array.forEach(this._plugins, function(p){
			if(!p || !p[func]){ return; }
			p[func].apply(p, args ? args : []);
		});
	},
	_parseProps: function(plugins){
		// summary:
		//		Parse plugins properties
		// plugins: Object
		//		Plugin properties defined by user
		if(!plugins){ return; }
		
		var p, loading = {}, options = this._options, grid = this.grid;
		var registry = _PluginManager.registry;//global plugin registry
		for(p in plugins){
			if(plugins[p]){//filter out boolean false e.g. {p:false}
				this._normalize(p, plugins, registry, loading);
			}
		}
		//"columnReordering" attribute won't work when DnD plugin is turned on.
		if(options.dnd){
			options.columnReordering = false;
		}
		
		//mixin all plugin properties into Grid
		lang.mixin(grid, options);
	},
	_normalize: function(p, plugins, registry, loading){
		// summary:
		//		Normalize plugin properties especially the dependency chain
		// p: String
		//		Plugin name
		// plugins: Object
		//		Plugin properties set by user
		// registry: Object
		//		The global plugin registry
		// loading: Object
		//		Map for checking process state
		if(!registry[p]){ throw new Error('Plugin ' + p + ' is required.');}
		
		if(loading[p]){ throw new Error('Recursive cycle dependency is not supported.'); }
		
		var options = this._options;
		if(options[p]){ return options[p]; }
		
		loading[p] = true;
		//TBD - more strict conditions?
		options[p] = lang.mixin({}, registry[p], lang.isObject(plugins[p]) ? plugins[p] : {});
		
		var dependencies = options[p]['dependency'];
		if(dependencies){
			if(!lang.isArray(dependencies)){
				dependencies = options[p]['dependency'] = [dependencies];
			}
			array.forEach(dependencies, function(dependency){
				if(!this._normalize(dependency, plugins, registry, loading)){
					throw new Error('Plugin ' + dependency + ' is required.');
				}
			}, this);
		}
		delete loading[p];
		return options[p];
	},
	_init: function(pre){
		// summary:
		//		Find appropriate plugins and load them
		// pre: Boolean
		//		True - preInit | False - postInit(by default)
		var p, preInit, options = this._options;
		for(p in options){
			preInit = options[p]['preInit'];
			if((pre ? preInit : !preInit) && options[p]['class'] && !this.pluginExisted(p)){
				this.loadPlugin(p);
			}
		}
	},
	loadPlugin: function(name){
		// summary:
		//		Load required plugin("name")
		// name: String
		//		Plugin name
		// returns: Object
		//		The newly loaded plugin
		var option = this._options[name];
		if(!option){ return null; } //return if no plugin option
		
		var plugin = this.getPlugin(name);
		if(plugin){ return plugin; } //return if plugin("name") already existed
		
		var dependencies = option['dependency'];
		array.forEach(dependencies, function(dependency){
			if(!this.loadPlugin(dependency)){
				throw new Error('Plugin ' + dependency + ' is required.');
			}
		}, this);
		var cls = option['class'];
		delete option['class'];//remove it for safety
		plugin = new this.getPluginClazz(cls)(this.grid, option);
		this._plugins.push(plugin);
		return plugin;
	},
	_initView: function(view){
		// summary:
		//		Overwrite several default behavior for each views(including _RowSelector view)
		if(!view){ return; }
		//add more events handler - _View
		util.funnelEvents(view.contentNode, view, "doContentEvent", ['mouseup', 'mousemove']);
		util.funnelEvents(view.headerNode, view, "doHeaderEvent", ['mouseup']);
	},
	pluginExisted: function(name){
		// summary:
		//		Check if plugin("name") existed
		// name: String
		//		Plugin name
		// returns: Boolean
		//		True - existed | False - not existed
		return !!this.getPlugin(name);
	},
	getPlugin: function(name){
		// summary:
		//		Get plugin("name")
		// name: String
		//		Plugin name
		// returns: Object
		//		Plugin instance
		var plugins = this._plugins;
		name = name.toLowerCase();
		for(var i = 0, len = plugins.length; i < len; i++){
			if(name == plugins[i]['name'].toLowerCase()){
				return plugins[i];
			}
		}
		return null;
	},
	getPluginClazz: function(clazz){
		// summary:
		//		Load target plugin which must be already required (require(..))
		// clazz: Class|String
		//		Plugin class
		if(lang.isFunction(clazz)){
			return clazz;//return if it's already a clazz
		}
		var errorMsg = 'Please make sure Plugin "' + clazz + '" is existed.';
		try{
			var cls = lang.getObject(clazz);
			if(!cls){ throw new Error(errorMsg); }
			return cls;
		}catch(e){
			throw new Error(errorMsg);
		}
	},
	isFixedCell: function(cell){
		// summary:
		//		See if target cell(column) is fixed or not.
		// cell: Object
		//		Target cell(column)
		// returns: Boolean
		//		True - fixed| False - not fixed

		//target cell can use Boolean attributes named "isRowSelector" or "fixedPos" to mark it's a fixed cell(column)
		return cell && (cell.isRowSelector || cell.fixedPos);
	},
	destroy: function(){
		// summary:
		//		Destroy all resources
		array.forEach(this._connects, connect.disconnect);
		this.forEach('destroy');
		if(this.grid.unwrap){
			this.grid.unwrap();
		}
		delete this._connects;
		delete this._plugins;
		delete this._options;
	}
});

_PluginManager.registerPlugin = function(clazz, props){
		// summary:
		//		Register plugins - TODO, a better way rather than global registry?
		// clazz: String
		//		Full class name, e.g. "dojox.grid.enhanced.plugins.DnD"
		// props: Object?
		//		Plugin properties e.g. {"dependency": ["nestedSorting"], ...}
	if(!clazz){
		console.warn("Failed to register plugin, class missed!");
		return;
	}
	var cls = _PluginManager;
	cls.registry = cls.registry || {};
	cls.registry[clazz.prototype.name]/*plugin name*/ = lang.mixin({"class": clazz}, (props ? props : {}));
};

return _PluginManager;

});

},
'dojox/grid/enhanced/_Events':function(){
define("dojox/grid/enhanced/_Events", [
	"dojo/_base/kernel",
	"dojo/_base/declare",
	"dojo/keys",
	"dojo/_base/html",
	"dojo/_base/event",
	"dojox/grid/_Events"
], function(dojo, declare, keys, html, event, _Events){

return declare("dojox.grid.enhanced._Events", null, {
	// summary:
	//		Overwrite some default events of DataGrid
	//
	// description:
	//		Methods are copied or replaced for overwriting, this might be refined once
	//		an overall plugin architecture is set up for DataGrid.

	// _events: Object
	//		Method map cached from dojox.grid._Events().
	_events: null,

	// headerCellActiveClass: String
	//		css class to apply to grid header cells when activated(mouse down)
	headerCellActiveClass: 'dojoxGridHeaderActive',
	
	// cellActiveClass: String
	//		css class to apply to grid content cells when activated(mouse down)
	cellActiveClass: 'dojoxGridCellActive',
	
	// rowActiveClass: String
	//		css class to apply to grid rows when activated(mouse down)
	rowActiveClass: 'dojoxGridRowActive',

	constructor: function(inGrid){
		//TODO - extend dojox.grid._Events rather than mixin for 1.8
		this._events = new _Events();
		//mixin "this" to Grid
		inGrid.mixin(inGrid, this);
	},
	dokeyup: function(e){
		// summary:
		//		Grid key up event handler.
		// e: Event
		//		Un-decorated event object
		this.focus.currentArea().keyup(e);
	},
	onKeyDown: function(e){
		// summary:
		//		Overwritten, see dojox.grid._Events.onKeyDown();
		if(e.altKey || e.metaKey){ return; }
		var focus = this.focus;
		var editing = this.edit.isEditing();
		switch(e.keyCode){
			case keys.TAB:
				if(e.ctrlKey){ return; }
				focus.tab(e.shiftKey ? -1:1,e);
				break;
			case keys.UP_ARROW:
			case keys.DOWN_ARROW:
				if(editing){ return; }
				focus.currentArea().move(e.keyCode == keys.UP_ARROW ? -1 : 1, 0, e);
				break;
			case keys.LEFT_ARROW:
			case keys.RIGHT_ARROW:
				if(editing){ return; }
				var offset = (e.keyCode == keys.LEFT_ARROW) ? 1 : -1;
				if(html._isBodyLtr()){ offset *= -1; }
				focus.currentArea().move(0, offset, e);
				break;
			case keys.F10:
				if(this.menus && e.shiftKey){
					this.onRowContextMenu(e);
				}
				break;
			default:
				focus.currentArea().keydown(e);
				break;
		}
	},
	//TODO - make the following events more reasonalble - e.g. more accurate conditions
	//events for row selectors
	domouseup: function(e){
		if(e.cellNode){
			this.onMouseUp(e);
		}else{
			this.onRowSelectorMouseUp(e);
		}
	},
	domousedown: function(e){
		if(!e.cellNode){
			this.onRowSelectorMouseDown(e);
		}
	},
	onMouseUp: function(e){
		// summary:
		//		New - Event fired when mouse is up inside grid.
		// e: Event
		//		Decorated event object that contains reference to grid, cell, and rowIndex
		this[e.rowIndex == -1 ? "onHeaderCellMouseUp" : "onCellMouseUp"](e);
	},
	onCellMouseDown: function(e){
		// summary:
		//		Overwritten, see dojox.grid._Events.onCellMouseDown()
		html.addClass(e.cellNode, this.cellActiveClass);
		html.addClass(e.rowNode, this.rowActiveClass);
	},
	onCellMouseUp: function(e){
		// summary:
		//		New - Event fired when mouse is up inside content cell.
		// e: Event
		//		Decorated event object that contains reference to grid, cell, and rowIndex
		html.removeClass(e.cellNode, this.cellActiveClass);
		html.removeClass(e.rowNode, this.rowActiveClass);
	},
	onCellClick: function(e){
		// summary:
		//		Overwritten, see dojox.grid._Events.onCellClick()

		//invoke dojox.grid._Events.onCellClick()
		this._events.onCellClick.call(this, e);
		//move mouse events to the focus manager.
		this.focus.contentMouseEvent(e);//TODO
	},
	onCellDblClick: function(e){
		// summary:
		//		Overwritten, see dojox.grid._Events.onCellDblClick()
		if(this.pluginMgr.isFixedCell(e.cell)){ return; }
		if(this._click.length > 1 && (!this._click[0] || !this._click[1])){
			this._click[0] = this._click[1] = e;
		}
		//invoke dojox.grid._Events.onCellDblClick()
		this._events.onCellDblClick.call(this, e);
	},
	onRowClick: function(e){
		// summary:
		//		Overwritten, see dojox.grid._Events.onRowClick()
		this.edit.rowClick(e);
		if(!e.cell || !this.plugin('indirectSelection')){
			this.selection.clickSelectEvent(e);
		}
	},
	onRowContextMenu: function(e){
		// summary:
		//		Overwritten, see dojox.grid._Events.onRowContextMenu()
		if(!this.edit.isEditing() && this.menus){
			this.showMenu(e);
		}
	},
	onSelectedRegionContextMenu: function(e){
		// summary:
		//		New - Event fired when a selected region context menu is accessed via mouse right click.
		// e: Event
		//		Decorated event object which contains reference to grid and info of selected
		//		regions(selection type - row|column, selected index - [...])
		if(this.selectedRegionMenu){
			this.selectedRegionMenu._openMyself({
				target: e.target,
				coords: e.keyCode !== keys.F10 && "pageX" in e ? {
					x: e.pageX,
					y: e.pageY
				} : null
			});
			event.stop(e);
		}
	},
	onHeaderCellMouseOut: function(e){
		// summary:
		//		Overwritten, see dojox.grid._Events.onHeaderCellMouseOut()
		if(e.cellNode){
			html.removeClass(e.cellNode, this.cellOverClass);
			html.removeClass(e.cellNode, this.headerCellActiveClass);
		}
	},
	onHeaderCellMouseDown: function(e){
		// summary:
		//		Overwritten, see dojox.grid._Events.onHeaderCellMouseDown()
		if(e.cellNode){//TBD - apply to selection region for nested sorting?
			html.addClass(e.cellNode, this.headerCellActiveClass);
		}
	},
	onHeaderCellMouseUp: function(e){
		// summary:
		//		New event
		if(e.cellNode){
			html.removeClass(e.cellNode, this.headerCellActiveClass);
		}
	},
	onHeaderCellClick: function(e){
		// summary:
		//		Overwritten, see dojox.grid._Events.onHeaderCellClick()
		
		//move focus to header.
		this.focus.currentArea("header");
		//invoke dojox.grid._Events.onHeaderCellClick()
		if(!e.cell.isRowSelector){
			this._events.onHeaderCellClick.call(this, e);
		}
		//move mouse events to the focus manager.
		this.focus.headerMouseEvent(e);
	},
	onRowSelectorMouseDown: function(e){
		this.focus.focusArea("rowHeader", e);
	},
	
	onRowSelectorMouseUp: function(e){},
	
	//triggered in _View, see Selector plugin
	onMouseUpRow: function(e){
		if(e.rowIndex != -1){
			this.onRowMouseUp(e);
		}
	},
	onRowMouseUp: function(e){}
});
});
},
'dojox/grid/enhanced/_FocusManager':function(){
define("dojox/grid/enhanced/_FocusManager", [
	"dojo/_base/kernel",
	"dojo/_base/lang",
	"dojo/_base/declare",
	"dojo/_base/array",
	"dojo/_base/connect",
	"dojo/_base/event",
	"dojo/_base/sniff",
	"dojo/_base/html",
	"dojo/keys",
	"dijit/a11y",
	"dijit/focus",
	"../_FocusManager"
], function(dojo, lang, declare, array, connect, event, has, html, keys, dijitA11y, dijitFocus, _FocusManager){

var _FocusArea = declare("dojox.grid.enhanced._FocusArea", null, {
	// summary:
	//		This is a friend class of _FocusManager
/*=====
		// name: string
		//		Name of this area.
		name: "",

		onFocus: function(event, step){
			// summary:
			//		Called when this area logically gets focus.
			// event: Event object
			//		May be unavailable, should check before use.
			// step: Integer
			//		The distance in the tab sequence from last focused area to this area.
			// returns:
			//		whether this area is successfully focused. If not, the next area will get focus.

 			return true;
		},
		
		onBlur: function(event, step){
			// summary:
			//		Called when this area logically loses focus.
			// event: Event object
			//		May be unavailable, should check before use.
			// step: Integer
			//		The distance in the tab sequence from this area to the area to focus.
			// returns:
			//		If Boolean, means whether this area has successfully blurred. If not, the next area to focus is still this one.
			//		If String, means the next area to focus is given by this returned name.

			return true;
		},
		
		onMove: function(rowStep, colStep, event){
			// summary:
			//		Called when focus is moving around within this area.
			// rowStep: Integer
			// colStep: Integer
			// event: Event object
			//		May be unavailable, should check before use.
		},
		
		onKey: function(event, isBubble){
			// summary:
			//		Called when some key is pressed when focus is logically in this area.
			// event: Event object
			// isBubble: Boolean
			//		Whether is in bubble stage (true) or catch stage (false).
			// returns:
			//		If you do NOT want the event to propagate any further along the area stack, return exactly false.
			//		So if you return nothing (undefined), this event is still propagating.
			return true;
		},
		
		getRegions: function(){
			// summary:
			//		Define the small regions (dom nodes) in this area.
			// returns:
			//		Array of dom nodes.
		},
		
		onRegionFocus: function(event){
			// summary:
			//		Connected to the onfocus event of the defined regions (if any)
		},
		
		onRegionBlur: function(event){
			// summary:
			//		Connected to the onblur event of the defined regions (if any)
		},
=====*/
	constructor: function(area, focusManager){
		this._fm = focusManager;
		this._evtStack = [area.name];
		var dummy = function(){return true;};
		area.onFocus = area.onFocus || dummy;
		area.onBlur = area.onBlur || dummy;
		area.onMove = area.onMove || dummy;
		area.onKeyUp = area.onKeyUp || dummy;
		area.onKeyDown = area.onKeyDown || dummy;
		lang.mixin(this, area);
	},
	move: function(rowStep, colStep, evt){
		if(this.name){
			var i, len = this._evtStack.length;
			for(i = len - 1; i >= 0; --i){
				if(this._fm._areas[this._evtStack[i]].onMove(rowStep, colStep, evt) === false){
					return false;
				}
			}
		}
		return true;
	},
	_onKeyEvent: function(evt, funcName){
		if(this.name){
			var i, len = this._evtStack.length;
			for(i = len - 1; i >= 0; --i){
				if(this._fm._areas[this._evtStack[i]][funcName](evt, false) === false){
					return false;
				}
			}
			for(i = 0; i < len; ++i){
				if(this._fm._areas[this._evtStack[i]][funcName](evt, true) === false){
					return false;
				}
			}
		}
		return true;
	},
	keydown: function(evt){
		return this._onKeyEvent(evt, "onKeyDown");
	},
	keyup: function(evt){
		return this._onKeyEvent(evt, "onKeyUp");
	},
	contentMouseEventPlanner: function(){
		return 0;
	},
	headerMouseEventPlanner: function(){
		return 0;
	}
});

return declare("dojox.grid.enhanced._FocusManager", _FocusManager, {
	_stopEvent: function(evt){
		try{
			if(evt && evt.preventDefault){
				event.stop(evt);
			}
		}catch(e){}
	},
	
	constructor: function(grid){
		this.grid = grid;
		this._areas = {};
		this._areaQueue = [];
		this._contentMouseEventHandlers = [];
		this._headerMouseEventHandlers = [];
		this._currentAreaIdx = -1;
		this._gridBlured = true;
		this._connects.push(connect.connect(grid, "onBlur", this, "_doBlur"));
		this._connects.push(connect.connect(grid.scroller, "renderPage", this, "_delayedCellFocus"));
		
		this.addArea({
			name: "header",
			onFocus: lang.hitch(this, this.focusHeader),
			onBlur: lang.hitch(this, this._blurHeader),
			onMove: lang.hitch(this, this._navHeader),
			getRegions: lang.hitch(this, this._findHeaderCells),
			onRegionFocus: lang.hitch(this, this.doColHeaderFocus),
			onRegionBlur: lang.hitch(this, this.doColHeaderBlur),
			onKeyDown: lang.hitch(this, this._onHeaderKeyDown)
		});
		this.addArea({
			name: "content",
			onFocus: lang.hitch(this, this._focusContent),
			onBlur: lang.hitch(this, this._blurContent),
			onMove: lang.hitch(this, this._navContent),
			onKeyDown: lang.hitch(this, this._onContentKeyDown)
		});
		this.addArea({
			name: "editableCell",
			onFocus: lang.hitch(this, this._focusEditableCell),
			onBlur: lang.hitch(this, this._blurEditableCell),
			onKeyDown: lang.hitch(this, this._onEditableCellKeyDown),
			onContentMouseEvent: lang.hitch(this, this._onEditableCellMouseEvent),
			contentMouseEventPlanner: function(evt, areas){ return -1; }
		});
		this.placeArea("header");
		this.placeArea("content");
		this.placeArea("editableCell");
		this.placeArea("editableCell","above","content");
	},
	destroy: function(){
		for(var name in this._areas){
			var area = this._areas[name];
			array.forEach(area._connects, connect.disconnect);
			area._connects = null;
			if(area.uninitialize){
				area.uninitialize();
			}
		}
		this.inherited(arguments);
	},
	addArea: function(area){
		if(area.name && lang.isString(area.name)){
			if(this._areas[area.name]){
				//Just replace the original area, instead of remove it, so the position does not change.
				array.forEach(area._connects, connect.disconnect);
			}
			this._areas[area.name] = new _FocusArea(area, this);
			if(area.onHeaderMouseEvent){
				this._headerMouseEventHandlers.push(area.name);
			}
			if(area.onContentMouseEvent){
				this._contentMouseEventHandlers.push(area.name);
			}
		}
	},
	getArea: function(areaName){
		return this._areas[areaName];
	},
	_bindAreaEvents: function(){
		var area, hdl, areas = this._areas;
		array.forEach(this._areaQueue, function(name){
			area = areas[name];
			if(!area._initialized && lang.isFunction(area.initialize)){
				area.initialize();
				area._initialized = true;
			}
			if(area.getRegions){
				area._regions = area.getRegions() || [];
				array.forEach(area._connects || [], connect.disconnect);
				area._connects = [];
				array.forEach(area._regions, function(r){
					if(area.onRegionFocus){
						hdl = connect.connect(r, "onfocus", area.onRegionFocus);
						area._connects.push(hdl);
					}
					if(area.onRegionBlur){
						hdl = connect.connect(r, "onblur", area.onRegionBlur);
						area._connects.push(hdl);
					}
				});
			}
		});
	},
	removeArea: function(areaName){
		var area = this._areas[areaName];
		if(area){
			this.ignoreArea(areaName);
			var i = array.indexOf(this._contentMouseEventHandlers, areaName);
			if(i >= 0){
				this._contentMouseEventHandlers.splice(i, 1);
			}
			i = array.indexOf(this._headerMouseEventHandlers, areaName);
			if(i >= 0){
				this._headerMouseEventHandlers.splice(i, 1);
			}
			array.forEach(area._connects, connect.disconnect);
			if(area.uninitialize){
				area.uninitialize();
			}
			delete this._areas[areaName];
		}
	},
	currentArea: function(areaName, toBlurOld){
		// summary:
		//		Set current area to the one areaName refers.
		// areaName: String
		var idx, cai = this._currentAreaIdx;
		if(lang.isString(areaName) && (idx = array.indexOf(this._areaQueue, areaName)) >= 0){
			if(cai != idx){
				this.tabbingOut = false;
				if(toBlurOld && cai >= 0 && cai < this._areaQueue.length){
					this._areas[this._areaQueue[cai]].onBlur();
				}
				this._currentAreaIdx = idx;
			}
		}else{
			return (cai < 0 || cai >= this._areaQueue.length) ?
				new _FocusArea({}, this) :
				this._areas[this._areaQueue[this._currentAreaIdx]];
		}
		return null;
	},
	placeArea: function(name, pos, otherAreaName){
		// summary:
		//		Place the area refered by *name* at some logical position relative to an existing area.
		// example:
		//		placeArea("myarea","before"|"after",...)
		//		placeArea("myarea","below"|"above",...)
		if(!this._areas[name]){ return; }
		var idx = array.indexOf(this._areaQueue,otherAreaName);
		switch(pos){
			case "after":
				if(idx >= 0){ ++idx; }
				//intentional drop through
			case "before":
				if(idx >= 0){
					this._areaQueue.splice(idx,0,name);
					break;
				}
				//intentional drop through
			default:
				this._areaQueue.push(name);
				break;
			case "above":
				var isAbove = true;
				//intentional drop through
			case "below":
				var otherArea = this._areas[otherAreaName];
				if(otherArea){
					if(isAbove){
						otherArea._evtStack.push(name);
					}else{
						otherArea._evtStack.splice(0,0,name);
					}
				}
		}
	},
	ignoreArea: function(name){
		this._areaQueue = array.filter(this._areaQueue,function(areaName){
			return areaName != name;
		});
	},
	focusArea: function(/* int|string|areaObj */areaId,evt){
		var idx;
		if(typeof areaId == "number"){
			idx = areaId < 0 ? this._areaQueue.length + areaId : areaId;
		}else{
			idx = array.indexOf(this._areaQueue,
				lang.isString(areaId) ? areaId : (areaId && areaId.name));
		}
		if(idx < 0){ idx = 0; }
		var step = idx - this._currentAreaIdx;
		this._gridBlured = false;
		if(step){
			this.tab(step, evt);
		}else{
			this.currentArea().onFocus(evt, step);
		}
	},
	tab: function(step,evt){
		//console.log("===========tab",step,"curArea",this._currentAreaIdx,"areaCnt",this._areaQueue.length);
		this._gridBlured = false;
		this.tabbingOut = false;
		if(step === 0){
			return;
		}
		var cai = this._currentAreaIdx;
		var dir = step > 0 ? 1:-1;
		if(cai < 0 || cai >= this._areaQueue.length){
			cai = (this._currentAreaIdx += step);
		}else{
			var nextArea = this._areas[this._areaQueue[cai]].onBlur(evt,step);
			if(nextArea === true){
				cai = (this._currentAreaIdx += step);
			}else if(lang.isString(nextArea) && this._areas[nextArea]){
				cai = this._currentAreaIdx = array.indexOf(this._areaQueue,nextArea);
			}
		}
		//console.log("target area:",cai);
		for(; cai >= 0 && cai < this._areaQueue.length; cai += dir){
			this._currentAreaIdx = cai;
			if(this._areaQueue[cai] && this._areas[this._areaQueue[cai]].onFocus(evt,step)){
				//console.log("final target area:",this._currentAreaIdx);
				return;
			}
		}
		//console.log("tab out");
		this.tabbingOut = true;
		if(step < 0){
			this._currentAreaIdx = -1;
			dijitFocus.focus(this.grid.domNode);
		}else{
			this._currentAreaIdx = this._areaQueue.length;
			dijitFocus.focus(this.grid.lastFocusNode);
		}
	},
	_onMouseEvent: function(type, evt){
		var lowercase = type.toLowerCase(),
			handlers = this["_" + lowercase + "MouseEventHandlers"],
			res = array.map(handlers, function(areaName){
				return {
					"area": areaName,
					"idx": this._areas[areaName][lowercase + "MouseEventPlanner"](evt, handlers)
				};
			}, this).sort(function(a, b){
				return b.idx - a.idx;
			}),
			resHandlers = array.map(res, function(handler){
				return res.area;
			}),
			i = res.length;
		while(--i >= 0){
			if(this._areas[res[i].area]["on" + type + "MouseEvent"](evt, resHandlers) === false){
				return;
			}
		}
	},
	contentMouseEvent: function(evt){
		this._onMouseEvent("Content", evt);
	},
	headerMouseEvent: function(evt){
		this._onMouseEvent("Header", evt);
	},
	initFocusView: function(){
		// summary:
		//		Overwritten
		this.focusView = this.grid.views.getFirstScrollingView() || this.focusView || this.grid.views.views[0];
		this._bindAreaEvents();
	},
	isNavHeader: function(){
		// summary:
		//		Overwritten
		//		Check whether currently navigating among column headers.
		// returns:
		//		true - focus is on a certain column header | false otherwise
		return this._areaQueue[this._currentAreaIdx] == "header";
	},
	previousKey: function(e){
		// summary:
		//		Overwritten
		this.tab(-1,e);
	},
	nextKey: function(e){
		// summary:
		//		Overwritten
		this.tab(1,e);
	},
	setFocusCell: function(/* Object */inCell, /* Integer */inRowIndex){
		// summary:
		//		Overwritten - focuses the given grid cell
		if(inCell){
			this.currentArea(this.grid.edit.isEditing() ? "editableCell" : "content", true);
			//This is very slow when selecting cells!
			//this.focusGridView();
			this._focusifyCellNode(false);
			this.cell = inCell;
			this.rowIndex = inRowIndex;
			this._focusifyCellNode(true);
		}
		this.grid.onCellFocus(this.cell, this.rowIndex);
	},
	doFocus: function(e){
		// summary:
		//		Overwritten
		//		trap focus only for grid dom node
		//		do not focus for scrolling if grid is about to blur
		if(e && e.target == e.currentTarget && !this.tabbingOut){
			if(this._gridBlured){
				this._gridBlured = false;
				if(this._currentAreaIdx < 0 || this._currentAreaIdx >= this._areaQueue.length){
					this.focusArea(0, e);
				}else{
					this.focusArea(this._currentAreaIdx, e);
				}
			}
		}else{
			this.tabbingOut = false;
		}
		event.stop(e);
	},
	_doBlur: function(){
		this._gridBlured = true;
	},
	doLastNodeFocus: function(e){
		// summary:
		//		Overwritten
		if(this.tabbingOut){
			this.tabbingOut = false;
		}else{
			this.focusArea(-1, e);
		}
	},
	_delayedHeaderFocus: function(){
		// summary:
		//		Overwritten
		if(this.isNavHeader() && !has('ie')){
			this.focusHeader();
		}
	},
	_delayedCellFocus: function(){
		// summary:
		//		Overwritten
		
		//If focus header here, the page will scroll to grid when the grid is created.
		//this.focusArea("header");
	},
	_changeMenuBindNode: function(oldBindNode, newBindNode){
		var hm = this.grid.headerMenu;
		if(hm && this._contextMenuBindNode == oldBindNode){
			hm.unBindDomNode(oldBindNode);
			hm.bindDomNode(newBindNode);
			this._contextMenuBindNode = newBindNode;
		}
	},
	//---------------Header Area------------------------------------------
	focusHeader: function(evt, step){ //need a further look why these changes to parent's
		// summary:
		//		Overwritten
		var didFocus = false;
		this.inherited(arguments);
		if(this._colHeadNode && html.style(this._colHeadNode, 'display') != "none"){
			dijitFocus.focus(this._colHeadNode);
			this._stopEvent(evt);
			didFocus = true;
		}
		return didFocus;
	},
	_blurHeader: function(evt,step){
		// summary:
		//		Overwritten
		if(this._colHeadNode){
			html.removeClass(this._colHeadNode, this.focusClass);
		}
		html.removeAttr(this.grid.domNode,"aria-activedescendant");
		// reset contextMenu onto viewsHeaderNode so right mouse on header will invoke (see focusHeader)
		this._changeMenuBindNode(this.grid.domNode,this.grid.viewsHeaderNode);
		//moved here from nextKey
		this._colHeadNode = this._colHeadFocusIdx = null;
		return true;
	},
	_navHeader: function(rowStep, colStep, evt){
		var colDir = colStep < 0 ? -1 : 1,
			savedIdx = array.indexOf(this._findHeaderCells(), this._colHeadNode);
		if(savedIdx >= 0 && (evt.shiftKey && evt.ctrlKey)){
			this.colSizeAdjust(evt, savedIdx, colDir * 5);
			return;
		}
		this.move(rowStep, colStep);
	},
	_onHeaderKeyDown: function(e, isBubble){
		if(isBubble){
			var dk = keys;
			switch(e.keyCode){
				case dk.ENTER:
				case dk.SPACE:
					var colIdx = this.getHeaderIndex();
					if(colIdx >= 0 && !this.grid.pluginMgr.isFixedCell(e.cell)/*TODO*/){
						this.grid.setSortIndex(colIdx, null, e);
						event.stop(e);
					}
					break;
			}
		}
		return true;
	},
	_setActiveColHeader: function(){
		// summary:
		//		Overwritten
		this.inherited(arguments);
		//EDG now will decorate event on header key events, if no focus, the cell will be wrong
		dijitFocus.focus(this._colHeadNode);
	},
	//---------------Content Area------------------------------------------
	findAndFocusGridCell: function(){
		// summary:
		//		Overwritten
		this._focusContent();
	},
	_focusContent: function(evt,step){
		var didFocus = true;
		var isEmpty = (this.grid.rowCount === 0); // If grid is empty this.grid.rowCount == 0
		if(this.isNoFocusCell() && !isEmpty){
			//skip all the hidden cells
			for(var i = 0, cell = this.grid.getCell(0); cell && cell.hidden; cell = this.grid.getCell(++i)){}
			this.setFocusIndex(0, cell ? i : 0);
		}else if(this.cell && !isEmpty){
			if(this.focusView && !this.focusView.rowNodes[this.rowIndex]){
				// if rowNode for current index is undefined (likely as a result of a sort and because of #7304)
				// scroll to that row
				this.grid.scrollToRow(this.rowIndex);
				this.focusGrid();
			}else{
				this.setFocusIndex(this.rowIndex, this.cell.index);
			}
		}else{
			didFocus = false;
		}
		if(didFocus){ this._stopEvent(evt); }
		return didFocus;
	},
	_blurContent: function(evt,step){
		this._focusifyCellNode(false);
		return true;
	},
	_navContent: function(rowStep, colStep, evt){
		if((this.rowIndex === 0 && rowStep < 0) || (this.rowIndex === this.grid.rowCount - 1 && rowStep > 0)){
			return;
		}
		this._colHeadNode = null;
		this.move(rowStep, colStep, evt);
		if(evt){
			event.stop(evt);
		}
	},
	_onContentKeyDown: function(e, isBubble){
		if(isBubble){
			var dk = keys, s = this.grid.scroller;
			switch(e.keyCode){
				case dk.ENTER:
				case dk.SPACE:
					var g = this.grid;
					if(g.indirectSelection){ break; }
					g.selection.clickSelect(this.rowIndex, connect.isCopyKey(e), e.shiftKey);
					g.onRowClick(e);
					event.stop(e);
					break;
				case dk.PAGE_UP:
					if(this.rowIndex !== 0){
						if(this.rowIndex != s.firstVisibleRow + 1){
							this._navContent(s.firstVisibleRow - this.rowIndex, 0);
						}else{
							this.grid.setScrollTop(s.findScrollTop(this.rowIndex - 1));
							this._navContent(s.firstVisibleRow - s.lastVisibleRow + 1, 0);
						}
						event.stop(e);
					}
					break;
				case dk.PAGE_DOWN:
					if(this.rowIndex + 1 != this.grid.rowCount){
						event.stop(e);
						if(this.rowIndex != s.lastVisibleRow - 1){
							this._navContent(s.lastVisibleRow - this.rowIndex - 1, 0);
						}else{
							this.grid.setScrollTop(s.findScrollTop(this.rowIndex + 1));
							this._navContent(s.lastVisibleRow - s.firstVisibleRow - 1, 0);
						}
						event.stop(e);
					}
					break;
			}
		}
		return true;
	},
	//------------------editable content area-------------------------
	_blurFromEditableCell: false,
	_isNavigating: false,
	_navElems: null,
	_focusEditableCell: function(evt,step){
		var didFocus = false;
		if(this._isNavigating){
			didFocus = true;
		}else if(this.grid.edit.isEditing() && this.cell){
			if(this._blurFromEditableCell || !this._blurEditableCell(evt, step)){
				this.setFocusIndex(this.rowIndex,this.cell.index);
				didFocus = true;
			}
			this._stopEvent(evt);
		}
		return didFocus;
	},
	_applyEditableCell: function(){
		try{
			this.grid.edit.apply();
		}catch(e){
			console.warn("_FocusManager._applyEditableCell() error:", e);
		}
	},
	_blurEditableCell: function(evt,step){
		this._blurFromEditableCell = false;
		if(this._isNavigating){
			var toBlur = true;
			if(evt){
				var elems = this._navElems;
				var firstElem = elems.lowest || elems.first;
				var lastElem = elems.last || elems.highest || firstElem;
				var target = has('ie') ? evt.srcElement : evt.target;
				toBlur = target == (step > 0 ? lastElem : firstElem);
			}
			if(toBlur){
				this._isNavigating = false;
				html.setSelectable(this.cell.getNode(this.rowIndex), false);
				return "content";
			}
			return false;
		}else if(this.grid.edit.isEditing() && this.cell){
			if(!step || typeof step != "number"){ return false; }
			var dir = step > 0 ? 1 : -1;
			var cc = this.grid.layout.cellCount;
			for(var cell, col = this.cell.index + dir; col >= 0 && col < cc; col += dir){
				cell = this.grid.getCell(col);
				if(cell.editable){
					this.cell = cell;
					this._blurFromEditableCell = true;
					return false;
				}
			}
			if((this.rowIndex > 0 || dir == 1) && (this.rowIndex < this.grid.rowCount || dir == -1)){
				this.rowIndex += dir;
				//this.cell = this.grid.getCell(0); //There must be an editable cell, so this is not necessary.
				for(col = dir > 0 ? 0 : cc - 1; col >= 0 && col < cc; col += dir){
					cell = this.grid.getCell(col);
					if(cell.editable){
						this.cell = cell;
						break;
					}
				}
				this._applyEditableCell();
				return "content";
			}
		}
		return true;
	},
	_initNavigatableElems: function(){
		this._navElems = dijitA11y._getTabNavigable(this.cell.getNode(this.rowIndex));
	},
	_onEditableCellKeyDown: function(e, isBubble){
		var dk = keys,
			g = this.grid,
			edit = g.edit,
			editApplied = false,
			toPropagate = true;
		switch(e.keyCode){
			case dk.ENTER:
				if(isBubble && edit.isEditing()){
					this._applyEditableCell();
					editApplied = true;
					event.stop(e);
				}
				//intentional drop through
			case dk.SPACE:
				if(!isBubble && this._isNavigating){
					toPropagate = false;
					break;
				}
				if(isBubble){
					if(!this.cell.editable && this.cell.navigatable){
						this._initNavigatableElems();
						var toFocus = this._navElems.lowest || this._navElems.first;
						if(toFocus){
							this._isNavigating = true;
							html.setSelectable(this.cell.getNode(this.rowIndex), true);
							dijitFocus.focus(toFocus);
							event.stop(e);
							this.currentArea("editableCell", true);
							break;
						}
					}
					if(!editApplied && !edit.isEditing() && !g.pluginMgr.isFixedCell(this.cell)){
						edit.setEditCell(this.cell, this.rowIndex);
					}
					if(editApplied){
						this.currentArea("content", true);
					}else if(this.cell.editable && g.canEdit()){
						this.currentArea("editableCell", true);
					}
				}
				break;
			case dk.PAGE_UP:
			case dk.PAGE_DOWN:
				if(!isBubble && edit.isEditing()){
					//prevent propagating to content area
					toPropagate = false;
				}
				break;
			case dk.ESCAPE:
				if(!isBubble){
					edit.cancel();
					this.currentArea("content", true);
				}
		}
		return toPropagate;
	},
	_onEditableCellMouseEvent: function(evt){
		if(evt.type == "click"){
			var cell = this.cell || evt.cell;
			if(cell && !cell.editable && cell.navigatable){
				this._initNavigatableElems();
				if(this._navElems.lowest || this._navElems.first){
					var target = has('ie') ? evt.srcElement : evt.target;
					if(target != cell.getNode(evt.rowIndex)){
						this._isNavigating = true;
						this.focusArea("editableCell", evt);
						html.setSelectable(cell.getNode(evt.rowIndex), true);
						dijitFocus.focus(target);
						return false;
					}
				}
			}else if(this.grid.singleClickEdit){
				this.currentArea("editableCell");
				return false;
			}
		}
		return true;
	}
});
});
},
'dojox/grid/enhanced/plugins/_SelectionPreserver':function(){
define("dojox/grid/enhanced/plugins/_SelectionPreserver", [
	"dojo/_base/declare",
	"dojo/_base/lang",
	"dojo/_base/connect",
	'../../_SelectionPreserver'
], function(declare, lang, connect, _SelectionPreserver){

return declare("dojox.grid.enhanced.plugins._SelectionPreserver", _SelectionPreserver, {
	// summary:
	//		Preserve selections across various user actions.
	//
	// description:
	//		Extends dojox.grid._SelectionPreserver adding a bit more support to make selection persistence working well
	//		with various EnhancedGrid features, e.g. filtering, nested sorting, pagination, select all etc.
	//
	//		Precondition - Identifier(id) is required for store, as id is used for differentiating row items.
	//		Known issue - The preserved selections might be inaccurate if some unloaded rows are previously selected by range(e.g.SHIFT + click)
	//
	// example:
	// |	//To turn on this - set 'keepSelection' attribute to true
	// |	<div dojoType="dojox.grid.EnhancedGrid" keepSelection = true .../>
	
	constructor: function(selection){
		var grid = this.grid;
		grid.onSelectedById = this.onSelectedById;
		this._oldClearData = grid._clearData;
		var self = this;
		grid._clearData = function(){
			self._updateMapping(!grid._noInternalMapping);
			self._trustSelection = [];
			self._oldClearData.apply(grid, arguments);
		};
		this._connects.push(
			connect.connect(selection, 'selectRange', lang.hitch(this, '_updateMapping', true, true, false)),
			connect.connect(selection, 'deselectRange', lang.hitch(this, '_updateMapping', true, false, false)),
			connect.connect(selection, 'deselectAll', lang.hitch(this, '_updateMapping', true, false, true))
		);
	},
	destroy: function(){
		this.inherited(arguments);
		this.grid._clearData = this._oldClearData;
	},
	reset: function(){
		this.inherited(arguments);
		this._idMap = [];
		this._trustSelection = [];
		this._defaultSelected = false;
	},
	_reSelectById: function(item, index){
		// summary:
		//		Overwritten
		var s = this.selection, g = this.grid;
		if(item && g._hasIdentity){
			var id = g.store.getIdentity(item);
			if(this._selectedById[id] === undefined){
				if(!this._trustSelection[index]){
					s.selected[index] = this._defaultSelected;
				}
			}else{
				s.selected[index] = this._selectedById[id];
			}
			this._idMap.push(id);
			g.onSelectedById(id, index, s.selected[index]);
		}
	},
	_selectById: function(toSelect, inItemOrIndex){
		// summary:
		//		Overwritten
		if(!this.inherited(arguments)){
			this._trustSelection[inItemOrIndex] = true;
		}
	},
	onSelectedById: function(id, rowIndex, value){},
	
	_updateMapping: function(trustSelection, isSelect, isForAll, from, to){
		// summary:
		//		This function try to keep the selection info updated when range selection is performed.
		//
		//		1. Calculate how many unloaded rows are there;
		//		2. update _selectedById data if grid.selection._selected can be trusted, so loaded but unselected rows can
		//			be properly recorded.
		var s = this.selection, g = this.grid, flag = 0, unloaded = 0, i, id;
		for(i = g.rowCount - 1; i >= 0; --i){
			if(!g._by_idx[i]){
				++unloaded;
				flag += s.selected[i] ? 1 : -1;
			}else{
				id = g._by_idx[i].idty;
				if(id && (trustSelection || this._selectedById[id] === undefined)){
					this._selectedById[id] = !!s.selected[i];
				}
			}
		}
		if(unloaded){
			this._defaultSelected = flag > 0;
		}
		if(!isForAll && from !== undefined && to !== undefined){
			isForAll = !g.usingPagination && Math.abs(to - from + 1) === g.rowCount;
		}
		// When deselectAll, make sure every thing is deselected, even if it was selected but not loaded now.
		// This occurs only when pagination's "All" is used.
		if(isForAll && (!g.usingPagination || g.selectionMode === 'single')){
			for(i = this._idMap.length - 1; i >= 0; --i){
				this._selectedById[this._idMap[i]] = isSelect;
			}
		}
	}
});
});

},
'dojox/grid/enhanced/plugins/IndirectSelection':function(){
define("dojox/grid/enhanced/plugins/IndirectSelection", [
	"dojo/_base/declare",
	"dojo/_base/array",
	"dojo/_base/event",
	"dojo/_base/lang",
	"dojo/_base/html",
	"dojo/_base/window",
	"dojo/_base/connect",
	"dojo/_base/sniff",
	"dojo/query",
	"dojo/keys",
	"dojo/string",
	"../_Plugin",
	"../../EnhancedGrid",
	"../../cells/dijit"
], function(declare, array, evt, lang, html, win, connect, has, query, keys, string, _Plugin, EnhancedGrid){

var gridCells = lang.getObject("dojox.grid.cells");

var RowSelector = declare("dojox.grid.cells.RowSelector", gridCells._Widget, {
	// summary:
	//		 Common attributes & functions for row selectors(Radio|CheckBox)

	// inputType: String
	//		Input type - Radio|CheckBox
	inputType: "",
	
	// map: Object
	//		Cache div refs of radio|checkbox to avoid querying each time
	map: null,
	
	// disabledMap: Object
	//		Cache index of disabled rows
	disabledMap: null,
	
	// isRowSelector: Boolean
	//		Marker of indirectSelection cell(column)
	isRowSelector: true,

	// _connects: Array
	//		List of all connections.
	_connects: null,
	
	// _subscribes: Array
	//		List of all subscribes.
	_subscribes: null,

	// checkedText: String
	//		Checked character for high contrast mode
	checkedText: '&#10003;',

	// unCheckedText: String
	//		Unchecked character for high contrast mode
	unCheckedText: 'O',

	constructor: function(){
		this.map = {}; this.disabledMap = {}; this.disabledCount= 0;
		this._connects = []; this._subscribes = [];
		this.inA11YMode = html.hasClass(win.body(), "dijit_a11y");
		
		this.baseClass = "dojoxGridRowSelector dijitReset dijitInline dijit" + this.inputType;
		this.checkedClass = " dijit" + this.inputType + "Checked";
		this.disabledClass = " dijit" + this.inputType + "Disabled";
		this.checkedDisabledClass = " dijit" + this.inputType + "CheckedDisabled";
		this.statusTextClass = " dojoxGridRowSelectorStatusText";//a11y use

		this._connects.push(connect.connect(this.grid, 'dokeyup', this, '_dokeyup'));
		this._connects.push(connect.connect(this.grid.selection, 'onSelected', this, '_onSelected'));
		this._connects.push(connect.connect(this.grid.selection, 'onDeselected', this, '_onDeselected'));
		this._connects.push(connect.connect(this.grid.scroller, 'invalidatePageNode', this, '_pageDestroyed'));
		this._connects.push(connect.connect(this.grid, 'onCellClick', this, '_onClick'));
		this._connects.push(connect.connect(this.grid, 'updateRow', this, '_onUpdateRow'));
	},
	formatter: function(data, rowIndex, scope){
		// summary:
		//		Overwritten, see dojox.grid.cells._Widget
		var _this = scope;
		var clazz = _this.baseClass;
		var checked = !!_this.getValue(rowIndex);
		var disabled = !!_this.disabledMap[rowIndex];//normalize 'undefined'
		
		if(checked){
			clazz += _this.checkedClass;
			if(disabled){ clazz += _this.checkedDisabledClass; }
		}else if(disabled){
			clazz += _this.disabledClass;
		}
		return ["<div tabindex = -1 ",
				"id = '" + _this.grid.id + "_rowSelector_" + rowIndex + "' ",
				"name = '" + _this.grid.id + "_rowSelector' class = '" + clazz + "' ",
				"role = " + _this.inputType + " aria-checked = '" + checked + "' aria-disabled = '" + disabled +
				"' aria-label = '" + string.substitute(_this.grid._nls["indirectSelection" + _this.inputType], [rowIndex + 1]) + "'>",
				"<span class = '" + _this.statusTextClass + "'>" + (checked ? _this.checkedText : _this.unCheckedText) + "</span>",
				"</div>"].join("");
	},
	setValue: function(rowIndex, inValue){
		// summary:
		//		Overwritten, see dojox.grid.cells._Widget
		//		Simply return, no action
	},
	getValue: function(rowIndex){
		// summary:
		//		Overwritten, see dojox.grid.cells._Widget
		return this.grid.selection.isSelected(rowIndex);
	},
	toggleRow: function(index, value){
		// summary:
		//		toggle checked | unchecked state for given row
		// index: Integer
		//		Row index
		// value: Boolean
		//		True - checked | False - unchecked
		this._nativeSelect(index, value);
	},
	setDisabled: function(index, disabled){
		// summary:
		//		toggle disabled | enabled state for given row
		// idx: Integer
		//		Row index
		// disabled: Boolean
		//		True - disabled | False - enabled
		if(index < 0){ return; }
		this._toggleDisabledStyle(index, disabled);
	},
	disabled: function(index){
		// summary:
		//		Check if one row is disabled
		return !!this.disabledMap[index];
	},
	_onClick: function(e){
		// summary:
		//		When mouse click on the selector cell, select/deselect the row.
		if(e.cell === this){
			this._selectRow(e);
		}
	},
	_dokeyup: function(e){
		// summary:
		//		Event handler for key up event
		//		- from dojox.grid.enhanced._Events.dokeyup()
		// e: Event
		//		Key up event
		if(e.cellIndex == this.index && e.rowIndex >= 0 && e.keyCode == keys.SPACE){
			this._selectRow(e);
		}
	},
	focus: function(rowIndex){
		// summary:
		//		Set focus to given row
		// rowIndex: Integer
		//		Target row
		var selector = this.map[rowIndex];
		if(selector){ selector.focus(); }
	},
	_focusEndingCell: function(rowIndex, cellIndex){
		// summary:
		//		Set focus to the ending grid cell(rowIndex,cellIndex) when swipe selection finished
		// rowIndex: Integer
		//		Row index
		// cellIndex: Integer
		//		Column index
		var cell = this.grid.getCell(cellIndex);
		this.grid.focus.setFocusCell(cell, rowIndex);
	},
	_nativeSelect: function(index, value){
		// summary:
		//		Use grid's native selection
		this.grid.selection[value ? 'select' : 'deselect'](index);
	},
	_onSelected: function(index){
		// summary:
		//		Triggered when a row is selected
		this._toggleCheckedStyle(index, true);
	},
	_onDeselected: function(index){
		// summary:
		//		Triggered when a row is deselected
		this._toggleCheckedStyle(index, false);
	},
	_onUpdateRow: function(index){
		// summary:
		//		Clear cache when row is re-built.
		delete this.map[index];
	},
	_toggleCheckedStyle: function(index, value){
		// summary:
		//		Change css styles for checked | unchecked
		var selector = this._getSelector(index);
		if(selector){
			html.toggleClass(selector, this.checkedClass, value);
			if(this.disabledMap[index]){
				html.toggleClass(selector, this.checkedDisabledClass, value);
			}
			selector.setAttribute("aria-checked", value);
			if(this.inA11YMode){
				selector.firstChild.innerHTML = (value ? this.checkedText : this.unCheckedText);
			}
		}
	},
	_toggleDisabledStyle: function(index, disabled){
		// summary:
		//		Change css styles for disabled | enabled
		var selector = this._getSelector(index);
		if(selector){
			html.toggleClass(selector, this.disabledClass, disabled);
			if(this.getValue(index)){
				html.toggleClass(selector, this.checkedDisabledClass, disabled);
			}
			selector.setAttribute("aria-disabled", disabled);
		}
		this.disabledMap[index] = disabled;
		if(index >= 0){
			this.disabledCount += disabled ? 1 : -1;
		}
	},
	_getSelector: function(index){
		// summary:
		//		Find selector for given row caching it if 1st time found
		var selector = this.map[index];
		if(!selector){//use accurate query for better performance
			var rowNode = this.view.rowNodes[index];
			if(rowNode){
				selector = query('.dojoxGridRowSelector', rowNode)[0];
				if(selector){ this.map[index] = selector; }
			}
		}
		return selector;
	},
	_pageDestroyed: function(pageIndex){
		// summary:
		//		Explicitly empty map cache when a page destroyed
		//		See dojox.grid._Scroller.invalidatePageNode()
		// pageIndex: Integer
		//		Index of destroyed page
		var rowsPerPage = this.grid.scroller.rowsPerPage;
		var start = pageIndex * rowsPerPage, end = start + rowsPerPage - 1;
		for(var i = start; i <= end; i++){
			if(!this.map[i]){continue;}
			html.destroy(this.map[i]);
			delete this.map[i];
		}
		//console.log("Page ",pageIndex, " destroyed, Map=",this.map);
	},
	destroy: function(){
		for(var i in this.map){
			html.destroy(this.map[i]);
			delete this.map[i];
		}
		for(i in this.disabledMap){ delete this.disabledMap[i]; }
		array.forEach(this._connects, connect.disconnect);
		array.forEach(this._subscribes, connect.unsubscribe);
		delete this._connects;
		delete this._subscribes;
		//console.log('Single(Multiple)RowSelector.destroy() executed!');
	}
});

var SingleRowSelector = declare("dojox.grid.cells.SingleRowSelector", RowSelector, {
	// summary:
	//		IndirectSelection cell(column) for single selection mode, using styles of dijit.form.RadioButton
	inputType: "Radio",

	_selectRow: function(e){
		// summary:
		//		Select the target row
		// e: Event
		//		Event fired on the target row
		var index = e.rowIndex;
		if(this.disabledMap[index]){ return; }
		this._focusEndingCell(index, e.cellIndex);
		this._nativeSelect(index, !this.grid.selection.selected[index]);
	}
});

var MultipleRowSelector = declare("dojox.grid.cells.MultipleRowSelector", RowSelector, {
	// summary:
	//		Indirect selection cell for multiple or extended mode, using dijit.form.CheckBox

	// inputType: String
	inputType: "CheckBox",
	
	// swipeStartRowIndex: Integer
	//		Start row index for swipe selection
	swipeStartRowIndex: -1,

	// swipeMinRowIndex: Integer
	//		Min row index for swipe selection
	swipeMinRowIndex: -1,
	
	// swipeMinRowIndex: Integer
	//		Max row index for swipe selection
	swipeMaxRowIndex: -1,
	
	// toSelect: Boolean
	//		new state for selection
	toSelect: false,
	
	// lastClickRowIdx: Integer
	//		Row index for last click, used for range selection via Shift + click
	lastClickRowIdx: -1,
		
	unCheckedText: '&#9633;',

	constructor: function(){
		this._connects.push(connect.connect(win.doc, 'onmouseup', this, '_domouseup'));
		this._connects.push(connect.connect(this.grid, 'onRowMouseOver', this, '_onRowMouseOver'));
		this._connects.push(connect.connect(this.grid.focus, 'move', this, '_swipeByKey'));
		this._connects.push(connect.connect(this.grid, 'onCellMouseDown', this, '_onMouseDown'));
		if(this.headerSelector){//option set by user to add a select-all checkbox in column header
			this._connects.push(connect.connect(this.grid.views, 'render', this, '_addHeaderSelector'));
			this._connects.push(connect.connect(this.grid, '_onFetchComplete', this, '_addHeaderSelector'));
			this._connects.push(connect.connect(this.grid, 'onSelectionChanged', this, '_onSelectionChanged'));
			this._connects.push(connect.connect(this.grid, 'onKeyDown', this, function(e){
				if(e.rowIndex == -1 && e.cellIndex == this.index && e.keyCode == keys.SPACE){
					this._toggletHeader();//TBD - a better way
				}
			}));
		}
	},
	toggleAllSelection:function(checked){
		// summary:
		//		Toggle select all|deselect all
		// checked: Boolean
		//		True - select all, False - deselect all
		var grid = this.grid, selection = grid.selection;
		if(checked){
			selection.selectRange(0, grid.rowCount-1);
		}else{
			selection.deselectAll();
		}
	},
	_onMouseDown: function(e){
		if(e.cell == this){
			this._startSelection(e.rowIndex);
			evt.stop(e);
		}
	},
	_onRowMouseOver: function(e){
		// summary:
		//		Event fired when mouse moves over a data row(outside of this column).
		//		- from dojox.grid.enhanced._Events.onRowMouseOver()
		// e: Event
		//		Decorated event object which contains reference to grid, cell, and rowIndex
		this._updateSelection(e, 0);
	},
	_domouseup: function(e){
		// summary:
		//		Event handler for mouse up event - from dojo.doc.domouseup()
		// e: Event
		//		Mouse up event
		if(has('ie')){
			this.view.content.decorateEvent(e);//TODO - why only e in IE hasn't been decorated?
		}
		var inSwipeSelection = e.cellIndex >= 0 && this.inSwipeSelection() && !this.grid.edit.isEditRow(e.rowIndex);
		if(inSwipeSelection){
			this._focusEndingCell(e.rowIndex, e.cellIndex);
		}
		this._finishSelect();
	},
	_dokeyup: function(e){
		// summary:
		//		Event handler for key up event
		//		- from dojox.grid.enhanced._Events.dokeyup()
		// e: Event
		//		Key up event
		this.inherited(arguments);
		if(!e.shiftKey){
			this._finishSelect();
		}
	},
	_startSelection: function(rowIndex){
		// summary:
		//		Initialize parameters to start a new swipe selection
		// rowIndex: Integer
		//		Index of the start row
		this.swipeStartRowIndex = this.swipeMinRowIndex = this.swipeMaxRowIndex = rowIndex;
		this.toSelect = !this.getValue(rowIndex);
	},
	_updateSelection: function(e, delta){
		// summary:
		//		Update row selections, fired during a swipe selection
		// e: Event
		//		Event of the current row,
		// delta: Integer
		//		Row index delta, used for swipe selection via Shift + Arrow key
		//		0: not via key, -1 : Shift +  Up, 1 : Shift + Down
		if(!this.inSwipeSelection()){ return; }
		
		var byKey = delta !== 0;//whether via Shift + Arrow Key
		var currRow = e.rowIndex, deltaRow = currRow - this.swipeStartRowIndex + delta;
		if(deltaRow > 0 && this.swipeMaxRowIndex < currRow + delta){
			this.swipeMaxRowIndex = currRow + delta;
		}
		if(deltaRow < 0 && this.swipeMinRowIndex > currRow + delta){
			this.swipeMinRowIndex = currRow + delta;
		}

		var min = deltaRow > 0 ? this.swipeStartRowIndex : currRow + delta;
		var max = deltaRow > 0 ? currRow + delta : this.swipeStartRowIndex;
		for(var i = this.swipeMinRowIndex; i <= this.swipeMaxRowIndex; i++){
			if(this.disabledMap[i] || i < 0){ continue; }
			if(i >= min && i <= max){//deltaRow != 0 || this.toSelect
				this._nativeSelect(i, this.toSelect);
			}else if(!byKey){
				this._nativeSelect(i, !this.toSelect);
			}
		}
	},
	_swipeByKey: function(rowOffset, colOffset, e){
		// summary:
		//		Update row selections, fired when Shift + Cursor is used for swipe selection
		//		See dojox.grid.enhanced._Events.onKeyDown
		// e: Event
		//		Event of the current row,
		// rowOffset: Integer
		//		Row offset, used for swipe selection via Shift + Cursor
		//		-1 : Shift +  Up, 1 : Shift + Down
		if(!e || rowOffset === 0 || !e.shiftKey || e.cellIndex != this.index ||
			this.grid.focus.rowIndex < 0){ //TBD - e.rowIndex == 0 && delta == -1
			return;
		}
		var rowIndex = e.rowIndex;
		if(this.swipeStartRowIndex < 0){
			//A new swipe selection starts via Shift + Arrow key
			this.swipeStartRowIndex = rowIndex;
			if(rowOffset > 0){//Shift + Down
				this.swipeMaxRowIndex = rowIndex + rowOffset;
				this.swipeMinRowIndex = rowIndex;
			}else{//Shift + UP
				this.swipeMinRowIndex = rowIndex + rowOffset;
				this.swipeMaxRowIndex = rowIndex;
			}
			this.toSelect = this.getValue(rowIndex);
		}
		this._updateSelection(e, rowOffset);
	},
	_finishSelect: function(){
		// summary:
		//		Reset parameters to end a swipe selection
		this.swipeStartRowIndex = -1;
		this.swipeMinRowIndex = -1;
		this.swipeMaxRowIndex = -1;
		this.toSelect = false;
	},
	inSwipeSelection: function(){
		// summary:
		//		Check if during a swipe selection
		// returns: Boolean
		//		Whether in swipe selection
		return this.swipeStartRowIndex >= 0;
	},
	_nativeSelect: function(index, value){
		// summary:
		//		Overwritten
		this.grid.selection[value ? 'addToSelection' : 'deselect'](index);
	},
	_selectRow: function(e){
		// summary:
		//		Select the target row or range or rows
		// e: Event
		//		Event fired on the target row
		var rowIndex = e.rowIndex;
		if(this.disabledMap[rowIndex]){ return; }
		evt.stop(e);
		this._focusEndingCell(rowIndex, e.cellIndex);
		
		var delta = rowIndex - this.lastClickRowIdx;
		var newValue = !this.grid.selection.selected[rowIndex];
		if(this.lastClickRowIdx >= 0 && !e.ctrlKey && !e.altKey && e.shiftKey){
			var min = delta > 0 ? this.lastClickRowIdx : rowIndex;
			var max = delta > 0 ? rowIndex : this.lastClickRowIdx;
			for(var i = min; i >= 0 && i <= max; i++){
				this._nativeSelect(i, newValue);
			}
		}else{
			this._nativeSelect(rowIndex, newValue);
		}
		this.lastClickRowIdx = rowIndex;
	},
	getValue: function(rowIndex){
		// summary:
		//		Overwritten
		if(rowIndex == -1){//header selector
			var g = this.grid;
			return g.rowCount > 0 && g.rowCount <= g.selection.getSelectedCount();
		}
		return this.inherited(arguments);
	},
	_addHeaderSelector: function(){
		// summary:
		//		Add selector in column header for selecting|deselecting all
		var headerCellNode = this.view.getHeaderCellNode(this.index);
		if(!headerCellNode){ return; }
		html.empty(headerCellNode);
		var g = this.grid;
		var selector = headerCellNode.appendChild(html.create("div", {
			'aria-label': g._nls["selectAll"],
			"tabindex": -1, "id": g.id + "_rowSelector_-1", "class": this.baseClass, "role": "Checkbox",
			"innerHTML": "<span class = '" + this.statusTextClass +
				"'></span><span style='height: 0; width: 0; overflow: hidden; display: block;'>" +
				g._nls["selectAll"] + "</span>"
		}));
		this.map[-1] = selector;
		var idx = this._headerSelectorConnectIdx;
		if(idx !== undefined){
			connect.disconnect(this._connects[idx]);
			this._connects.splice(idx, 1);
		}
		this._headerSelectorConnectIdx = this._connects.length;
		this._connects.push(connect.connect(selector, 'onclick', this, '_toggletHeader'));
		this._onSelectionChanged();
	},
	_toggletHeader: function(){
		// summary:
		//		Toggle state for head selector
		if(!!this.disabledMap[-1]){ return; }
		this.grid._selectingRange = true;
		this.toggleAllSelection(!this.getValue(-1));
		this._onSelectionChanged();
		this.grid._selectingRange = false;
	},
	_onSelectionChanged: function(){
		// summary:
		//		Update header selector anytime selection changed
		var g = this.grid;
		if(!this.map[-1] || g._selectingRange){ return; }
		g.allItemsSelected = this.getValue(-1);
		this._toggleCheckedStyle(-1, g.allItemsSelected);
	},
	_toggleDisabledStyle: function(index, disabled){
		// summary:
		//		Overwritten
		this.inherited(arguments);
		if(this.headerSelector){
			var allDisabled = (this.grid.rowCount == this.disabledCount);
			if(allDisabled != !!this.disabledMap[-1]){//only if needed
				arguments[0] = -1;
				arguments[1] = allDisabled;
				this.inherited(arguments);
			}
		}
	}
});

var IndirectSelection = declare("dojox.grid.enhanced.plugins.IndirectSelection", _Plugin, {
	// summary:
	//		A handy way for adding check boxe/radio button for rows, and selecting rows by swiping(or keyboard)

	// description:
	//		For better rendering performance, div(images) are used to simulate radio button|check boxes
	//
	// example:
	// |	<div dojoType="dojox.grid.EnhancedGrid" plugins="{indirectSelection: true}" ...></div>
	//		or
	// |	<div dojoType="dojox.grid.EnhancedGrid" plugins="{indirectSelection: {name: 'xxx', width:'30px', styles:'text-align: center;'}}" ...></div>

	// name: String
	//		Plugin name
	name: "indirectSelection",
	
	constructor: function(){
		// Hook layout.setStructure(), so that indirectSelection is always included
		var layout = this.grid.layout;
		this.connect(layout, 'setStructure', lang.hitch(layout, this.addRowSelectCell, this.option));
	},
	addRowSelectCell: function(option){
		// summary:
		//		Add indirectSelection cell(mapped to a column of radio button|check boxes)
		if(!this.grid.indirectSelection || this.grid.selectionMode == 'none'){
			return;
		}
		var rowSelectCellAdded = false, inValidFields = ['get', 'formatter', 'field', 'fields'],
		defaultCellDef = {type: MultipleRowSelector, name: '', width:'30px', styles:'text-align: center;'};
		if(option.headerSelector){ option.name = ''; }//mutual conflicting attrs

		if(this.grid.rowSelectCell){//remove the existed one
			this.grid.rowSelectCell.destroy();
		}
		
		array.forEach(this.structure, function(view){
			var cells = view.cells;
			if(cells && cells.length > 0 && !rowSelectCellAdded){
				var firstRow = cells[0];
				if(firstRow[0] && firstRow[0].isRowSelector){
					console.debug('addRowSelectCell() - row selector cells already added, return.');
					rowSelectCellAdded = true;
					return;
				}
				var selectDef, cellType = this.grid.selectionMode == 'single' ? SingleRowSelector : MultipleRowSelector;
				selectDef = lang.mixin(defaultCellDef, option, {type: cellType, editable: false, notselectable: true, filterable: false, navigatable: true, nosort: true});
				array.forEach(inValidFields, function(field){//remove invalid fields
					if(field in selectDef){ delete selectDef[field]; }
				});
				if(cells.length > 1){ selectDef.rowSpan = cells.length; }//for complicate layout
				array.forEach(this.cells, function(cell, i){
					if(cell.index >= 0){
						cell.index += 1;
						//console.debug('cell '+ (cell.index - 1) +  ' is updated to index ' + cell.index);
					}else{
						console.warn('Error:IndirectSelection.addRowSelectCell()-  cell ' + i + ' has no index!');
					}
				});
				var rowSelectCell = this.addCellDef(0, 0, selectDef);
				rowSelectCell.index = 0;
				firstRow.unshift(rowSelectCell);
				this.cells.unshift(rowSelectCell);
				this.grid.rowSelectCell = rowSelectCell;
				rowSelectCellAdded = true;
			}
		}, this);
		this.cellCount = this.cells.length;
	},
	destroy: function(){
		this.grid.rowSelectCell.destroy();
		delete this.grid.rowSelectCell;
		this.inherited(arguments);
	}
});

EnhancedGrid.registerPlugin(IndirectSelection/*name:'indirectSelection'*/, {"preInit": true});

return IndirectSelection;
});

},
'dojox/grid/enhanced/_Plugin':function(){
define("dojox/grid/enhanced/_Plugin", [
	"dojo/_base/kernel",
	"dojo/_base/lang",
	"dojo/_base/declare",
	"dojo/_base/array",
	"dojo/_base/connect",
    "../EnhancedGrid"
], function(dojo, lang, declare, array, connect){
	
return declare("dojox.grid.enhanced._Plugin", null, {
	// summary:
	//		Base class for all plugins.
	// description:
	//		Provides common plugin functionality and basic life cycle management.
	//
	//		Each concrete plugin must have a name field and is responsible for registering itself to the global plugin registry
	//		e.g. for dnd plugin:
	// |		dojox.grid.EnhancedGrid.registerPlugin("dnd" /*plugin name*/,
	// |												dojox.grid.enhanced.plugins.DnD /*full class name of a plugin*/
	// |												{"preInit": false, "dependency": ["nestedSorting"]} /*properties*/);
	//
	//		[Keywords] of plugin properties (case sensitive):
	//
	//		- "preInit": boolean, whether a plugin should be created before EnhancedGrid.postCreate(),
	//		   false by default(plugins are created after EnhancedGrid.postCreate()).
	//		- "dependency": array or string, plugin(s) indicated by "dependency" will be created before the current one.
	//		   Note: recursive cycle dependencies are not supported e.g. following dependency is invalid:
	//		   pluginA -> pluginB -> pluginA
	//
	// example:
	//		1. Customize default DnD plugin
	//
	// |	declare("mygrid.MyDnD", dojox.grid.enhanced.plugins.DnD, {
	// |		name:"dnd" //still reuse the plugin name
	// |		constructor: function(inGrid, option){ ... }
	// |	});
	// |	dojox.grid.EnhancedGrid.registerPlugin("dnd", mygrid.MyDnD);
	//
	//		2. Add new plugin - PluginA
	//
	// |	declare("mygrid.PluginA", dojox.grid.enhanced._Plugin, {
	// |		name: "pA",
	// |		constructor: function(inGrid, option){ ... }
	// |	});
	// |	dojox.grid.EnhancedGrid.registerPlugin("pA",mygrid.PluginA);
	//
	//		3. Use plugins
	//
	// |	dojo.require("mygrid.MyDnD");
	// |	dojo.require("mygrid.PluginA");
	// |
	// |	<script type="text/javascript">
	// |		var grid = new dojox.grid.EnhancedGrid(
	// |		{plugins: {dnd:true, pA:true}, ... }, dojo.byId("gridDiv"));
	// |		grid.startup();
	// |	</script>

	// name: String
	//		Plugin name, e.g. 'nestedSorting', 'dnd'...
	name: 'plugin',
	
	// grid: Object
	//		Grid that the plugin belongs to
	grid: null,

	// option: Object
	//		Plugin properties - leveraged with default and user specified properties.
	//		e.g. for dnd plugin, it may look like {"class": dojox.grid.enhanced.plugins.DnD, "dependency": ["nestedSorting"], ...}
	option: {},

	// _connects: Array
	//		List of all connections.
	_connects: [],
	
	// _subscribes: Array
	//		List of all subscribes.
	_subscribes: [],

	// privates: Object
	//		Private properties/methods shouldn't be mixin-ed anytime.
	privates: {},
	
	constructor: function(inGrid, option){
		this.grid = inGrid;
		this.option = option;
		this._connects = [];
		this._subscribes = [];
		this.privates = lang.mixin({},dojox.grid.enhanced._Plugin.prototype);
		this.init();
	},
	
	init: function(){},
	
	onPreInit: function(){},
	
	onPostInit: function(){},
	
	onStartUp: function(){},
	
	connect: function(obj, event, method){
		// summary:
		//		Connects specified obj/event to specified method of this object.
		// example:
		//	|	var plugin = new dojox.grid.enhanced._Plugin(grid,"myPlugin",{...});
		//	|	// when foo.bar() is called, call the listener in the scope of plugin
		//	|	plugin.connect(foo, "bar", function(){
		//	|		console.debug(this.xxx());//"this" - plugin scope
		//	|	});
		var conn = connect.connect(obj, event, this, method);
		this._connects.push(conn);
		return conn;
	},
	disconnect: function(handle){
		// summary:
		//		Disconnects handle and removes it from connection list.
		array.some(this._connects, function(conn, i, conns){
			if(conn == handle){
				connect.disconnect(handle);
				conns.splice(i, 1);
				return true;
			}
			return false;
		});
	},
	subscribe: function(topic, method){
		// summary:
		//		Subscribes to the specified topic and calls the specified method
		//		of this object.
		// example:
		//	|	var plugin = new dojox.grid.enhanced._Plugin(grid,"myPlugin",{...});
		//	|	// when /my/topic is published, call the subscriber in the scope of plugin
		//	|	// with passed parameter - "v"
		//	|	plugin.subscribe("/my/topic", function(v){
		//	|		console.debug(this.xxx(v));//"this" - plugin scope
		//	|	});
		var subscribe = connect.subscribe(topic, this, method);
		this._subscribes.push(subscribe);
		return subscribe;
	},
	unsubscribe: function(handle){
		// summary:
		//		Un-subscribes handle and removes it from subscriptions list.
		array.some(this._subscribes, function(subscribe, i, subscribes){
			if(subscribe == handle){
				connect.unsubscribe(handle);
				subscribes.splice(i, 1);
				return true;
			}
			return false;
		});
	},
	onSetStore: function(store){
		// summary:
		//		Called when store is changed.
	},
	destroy: function(){
		// summary:
		//		Destroy all resources.
		array.forEach(this._connects, connect.disconnect);
		array.forEach(this._subscribes, connect.unsubscribe);
		delete this._connects;
		delete this._subscribes;
		delete this.option;
		delete this.privates;
		//console.log('Plugin [', this.name, '].destroy() executed!');
	}
});

//Each plugin is responsible for registering itself
// e.g. for DnD plugin(name:'dnd'):
// |	dojox.grid.EnhancedGrid.registerPlugin(dojox.grid.enhanced.plugins.DnD/*class*/,
// |		{"dependency": ["nestedSorting"]}/*Optional - properties*/);

});

},
'dojox/grid/cells/dijit':function(){
define("dojox/grid/cells/dijit", [
	"dojo/_base/kernel",
	"../../main",
	"dojo/_base/declare",
	"dojo/_base/array",
	"dojo/_base/lang",
	"dojo/_base/json",
	"dojo/_base/connect",
	"dojo/_base/sniff",
	"dojo/dom",
	"dojo/dom-attr",
	"dojo/dom-construct",
	"dojo/dom-style",
	"dojo/dom-geometry",
	"dojo/data/ItemFileReadStore",
	"dijit/form/DateTextBox",
	"dijit/form/TimeTextBox",
	"dijit/form/ComboBox",
	"dijit/form/CheckBox",
	"dijit/form/TextBox",
	"dijit/form/NumberSpinner",
	"dijit/form/NumberTextBox",
	"dijit/form/CurrencyTextBox",
	"dijit/form/HorizontalSlider",
	"dijit/form/_TextBoxMixin",
	"dijit/Editor",
	"../util",
	"./_base"
], function(dojo, dojox, declare, array, lang, json, connect, has, dom, domAttr, domConstruct, domStyle,
	domGeometry, ItemFileReadStore, DateTextBox, TimeTextBox, ComboBox, CheckBox, TextBox,
	NumberSpinner, NumberTextBox, CurrencyTextBox, HorizontalSlider, _TextBoxMixin, Editor, util, BaseCell){
		
// TODO: shouldn't it be the test file's job to require these modules,
// if it is using them?  Most of these modules aren't referenced by this file.
	
	var _Widget = declare("dojox.grid.cells._Widget", BaseCell, {
		widgetClass: TextBox,
		constructor: function(inCell){
			this.widget = null;
			if(typeof this.widgetClass == "string"){
				dojo.deprecated("Passing a string to widgetClass is deprecated", "pass the widget class object instead", "2.0");
				this.widgetClass = lang.getObject(this.widgetClass);
			}
		},
		formatEditing: function(inDatum, inRowIndex){
			this.needFormatNode(inDatum, inRowIndex);
			return "<div></div>";
		},
		getValue: function(inRowIndex){
			return this.widget.get('value');
		},
		_unescapeHTML: function(value){
			return (value && value.replace && this.grid.escapeHTMLInData) ? 
					value.replace(/&lt;/g, '<').replace(/&amp;/g, '&') : value;
		},
		setValue: function(inRowIndex, inValue){
			if(this.widget&&this.widget.set){
				inValue = this._unescapeHTML(inValue);
				//Look for lazy-loading editor and handle it via its deferred.
				if(this.widget.onLoadDeferred){
					var self = this;
					this.widget.onLoadDeferred.addCallback(function(){
						self.widget.set("value",inValue===null?"":inValue);
					});
				}else{
					this.widget.set("value", inValue);
				}
			}else{
				this.inherited(arguments);
			}
		},
		getWidgetProps: function(inDatum){
			return lang.mixin(
				{
					dir: this.dir,
					lang: this.lang
				},
				this.widgetProps||{},
				{
					constraints: lang.mixin({}, this.constraint) || {}, //TODO: really just for ValidationTextBoxes
					required: (this.constraint || {}).required,
					value: this._unescapeHTML(inDatum)
				}
			);
		},
		createWidget: function(inNode, inDatum, inRowIndex){
			return new this.widgetClass(this.getWidgetProps(inDatum), inNode);
		},
		attachWidget: function(inNode, inDatum, inRowIndex){
			inNode.appendChild(this.widget.domNode);
			this.setValue(inRowIndex, inDatum);
		},
		formatNode: function(inNode, inDatum, inRowIndex){
			if(!this.widgetClass){
				return inDatum;
			}
			if(!this.widget){
				this.widget = this.createWidget.apply(this, arguments);
			}else{
				this.attachWidget.apply(this, arguments);
			}
			this.sizeWidget.apply(this, arguments);
			this.grid.views.renormalizeRow(inRowIndex);
			this.grid.scroller.rowHeightChanged(inRowIndex, true/*fix #11101*/);
			this.focus();
			return undefined;
		},
		sizeWidget: function(inNode, inDatum, inRowIndex){
			var p = this.getNode(inRowIndex);
			dojo.marginBox(this.widget.domNode, {w: domStyle.get(p, 'width')});
		},
		focus: function(inRowIndex, inNode){
			if(this.widget){
				setTimeout(lang.hitch(this.widget, function(){
					util.fire(this, "focus");
					if(this.focusNode && this.focusNode.tagName === "INPUT"){
						_TextBoxMixin.selectInputText(this.focusNode);
					}
				}), 0);
			}
		},
		_finish: function(inRowIndex){
			this.inherited(arguments);
			util.removeNode(this.widget.domNode);
			if(has('ie')){
				dom.setSelectable(this.widget.domNode, true);
			}
		}
	});
	_Widget.markupFactory = function(node, cell){
		BaseCell.markupFactory(node, cell);
		var widgetProps = lang.trim(domAttr.get(node, "widgetProps")||"");
		var constraint = lang.trim(domAttr.get(node, "constraint")||"");
		var widgetClass = lang.trim(domAttr.get(node, "widgetClass")||"");
		if(widgetProps){
			cell.widgetProps = json.fromJson(widgetProps);
		}
		if(constraint){
			cell.constraint = json.fromJson(constraint);
		}
		if(widgetClass){
			cell.widgetClass = lang.getObject(widgetClass);
		}
	};

	var ComboBox = declare("dojox.grid.cells.ComboBox", _Widget, {
		widgetClass: ComboBox,
		getWidgetProps: function(inDatum){
			var items=[];
			array.forEach(this.options, function(o){
				items.push({name: o, value: o});
			});
			var store = new ItemFileReadStore({data: {identifier:"name", items: items}});
			return lang.mixin({}, this.widgetProps||{}, {
				value: inDatum,
				store: store
			});
		},
		getValue: function(){
			var e = this.widget;
			// make sure to apply the displayed value
			e.set('displayedValue', e.get('displayedValue'));
			return e.get('value');
		}
	});
	ComboBox.markupFactory = function(node, cell){
		_Widget.markupFactory(node, cell);
		var options = lang.trim(domAttr.get(node, "options")||"");
		if(options){
			var o = options.split(',');
			if(o[0] != options){
				cell.options = o;
			}
		}
	};

	var DateTextBox = declare("dojox.grid.cells.DateTextBox", _Widget, {
		widgetClass: DateTextBox,
		setValue: function(inRowIndex, inValue){
			if(this.widget){
				this.widget.set('value', new Date(inValue));
			}else{
				this.inherited(arguments);
			}
		},
		getWidgetProps: function(inDatum){
			return lang.mixin(this.inherited(arguments), {
				value: new Date(inDatum)
			});
		}
	});
	DateTextBox.markupFactory = function(node, cell){
		_Widget.markupFactory(node, cell);
	};

	var CheckBox = declare("dojox.grid.cells.CheckBox", _Widget, {
		widgetClass: CheckBox,
		getValue: function(){
			return this.widget.checked;
		},
		setValue: function(inRowIndex, inValue){
			if(this.widget&&this.widget.attributeMap.checked){
				this.widget.set("checked", inValue);
			}else{
				this.inherited(arguments);
			}
		},
		sizeWidget: function(inNode, inDatum, inRowIndex){
			return;
		}
	});
	CheckBox.markupFactory = function(node, cell){
		_Widget.markupFactory(node, cell);
	};

	var Editor = declare("dojox.grid.cells.Editor", _Widget, {
		widgetClass: Editor,
		getWidgetProps: function(inDatum){
			return lang.mixin({}, this.widgetProps||{}, {
				height: this.widgetHeight || "100px"
			});
		},
		createWidget: function(inNode, inDatum, inRowIndex){
			// widget needs its value set after creation
			var widget = new this.widgetClass(this.getWidgetProps(inDatum), inNode);
			// use onLoadDeferred because onLoad may have already fired
			widget.onLoadDeferred.then(lang.hitch(this, 'populateEditor'));
			return widget;
		},
		formatNode: function(inNode, inDatum, inRowIndex){
			this.content = inDatum;
			this.inherited(arguments);
			if(has('mozilla')){
				// FIXME: seem to need to reopen the editor and display the toolbar
				var e = this.widget;
				e.open();
				if(this.widgetToolbar){
					domConstruct.place(e.toolbar.domNode, e.editingArea, "before");
				}
			}
		},
		populateEditor: function(){
			this.widget.set('value', this.content);
			this.widget.placeCursorAtEnd();
		}
	});
	Editor.markupFactory = function(node, cell){
		_Widget.markupFactory(node, cell);
		var h = lang.trim(domAttr.get(node, "widgetHeight")||"");
		if(h){
			if((h != "auto")&&(h.substr(-2) != "em")){
				h = parseInt(h, 10)+"px";
			}
			cell.widgetHeight = h;
		}
	};

	return dojox.grid.cells.dijit;

});

},
'dojo/data/ItemFileReadStore':function(){
define("dojo/data/ItemFileReadStore", ["../_base/kernel", "../_base/lang", "../_base/declare", "../_base/array", "../_base/xhr",
	"../Evented", "./util/filter", "./util/simpleFetch", "../date/stamp"
], function(kernel, lang, declare, array, xhr, Evented, filterUtil, simpleFetch, dateStamp){

// module:
//		dojo/data/ItemFileReadStore

var ItemFileReadStore = declare("dojo.data.ItemFileReadStore", [Evented],{
	// summary:
	//		The ItemFileReadStore implements the dojo/data/api/Read API and reads
	//		data from JSON files that have contents in this format --
	// |	{ items: [
	// |		{ name:'Kermit', color:'green', age:12, friends:['Gonzo', {_reference:{name:'Fozzie Bear'}}]},
	// |		{ name:'Fozzie Bear', wears:['hat', 'tie']},
	// |		{ name:'Miss Piggy', pets:'Foo-Foo'}
	// |	]}
	//		Note that it can also contain an 'identifier' property that specified which attribute on the items
	//		in the array of items that acts as the unique identifier for that item.

	constructor: function(/* Object */ keywordParameters){
		// summary:
		//		constructor
		// keywordParameters:
		//		{url: String} {data: jsonObject} {typeMap: object}
		//		The structure of the typeMap object is as follows:
		// |	{
		// |		type0: function || object,
		// |		type1: function || object,
		// |		...
		// |		typeN: function || object
		// |	}
		//		Where if it is a function, it is assumed to be an object constructor that takes the
		//		value of _value as the initialization parameters.  If it is an object, then it is assumed
		//		to be an object of general form:
		// |	{
		// |		type: function, //constructor.
		// |		deserialize:	function(value) //The function that parses the value and constructs the object defined by type appropriately.
		// |	}

		this._arrayOfAllItems = [];
		this._arrayOfTopLevelItems = [];
		this._loadFinished = false;
		this._jsonFileUrl = keywordParameters.url;
		this._ccUrl = keywordParameters.url;
		this.url = keywordParameters.url;
		this._jsonData = keywordParameters.data;
		this.data = null;
		this._datatypeMap = keywordParameters.typeMap || {};
		if(!this._datatypeMap['Date']){
			//If no default mapping for dates, then set this as default.
			//We use the dojo/date/stamp here because the ISO format is the 'dojo way'
			//of generically representing dates.
			this._datatypeMap['Date'] = {
				type: Date,
				deserialize: function(value){
					return dateStamp.fromISOString(value);
				}
			};
		}
		this._features = {'dojo.data.api.Read':true, 'dojo.data.api.Identity':true};
		this._itemsByIdentity = null;
		this._storeRefPropName = "_S"; // Default name for the store reference to attach to every item.
		this._itemNumPropName = "_0"; // Default Item Id for isItem to attach to every item.
		this._rootItemPropName = "_RI"; // Default Item Id for isItem to attach to every item.
		this._reverseRefMap = "_RRM"; // Default attribute for constructing a reverse reference map for use with reference integrity
		this._loadInProgress = false; //Got to track the initial load to prevent duelling loads of the dataset.
		this._queuedFetches = [];
		if(keywordParameters.urlPreventCache !== undefined){
			this.urlPreventCache = keywordParameters.urlPreventCache?true:false;
		}
		if(keywordParameters.hierarchical !== undefined){
			this.hierarchical = keywordParameters.hierarchical?true:false;
		}
		if(keywordParameters.clearOnClose){
			this.clearOnClose = true;
		}
		if("failOk" in keywordParameters){
			this.failOk = keywordParameters.failOk?true:false;
		}
	},

	url: "",	// use "" rather than undefined for the benefit of the parser (#3539)

	//Internal var, crossCheckUrl.  Used so that setting either url or _jsonFileUrl, can still trigger a reload
	//when clearOnClose and close is used.
	_ccUrl: "",

	data: null,	// define this so that the parser can populate it

	typeMap: null, //Define so parser can populate.

	// clearOnClose: Boolean
	//		Parameter to allow users to specify if a close call should force a reload or not.
	//		By default, it retains the old behavior of not clearing if close is called.  But
	//		if set true, the store will be reset to default state.  Note that by doing this,
	//		all item handles will become invalid and a new fetch must be issued.
	clearOnClose: false,

	// urlPreventCache: Boolean
	//		Parameter to allow specifying if preventCache should be passed to the xhrGet call or not when loading data from a url.
	//		Note this does not mean the store calls the server on each fetch, only that the data load has preventCache set as an option.
	//		Added for tracker: #6072
	urlPreventCache: false,

	// failOk: Boolean
	//		Parameter for specifying that it is OK for the xhrGet call to fail silently.
	failOk: false,

	// hierarchical: Boolean
	//		Parameter to indicate to process data from the url as hierarchical
	//		(data items can contain other data items in js form).  Default is true
	//		for backwards compatibility.  False means only root items are processed
	//		as items, all child objects outside of type-mapped objects and those in
	//		specific reference format, are left straight JS data objects.
	hierarchical: true,

	_assertIsItem: function(/* dojo/data/api/Item */ item){
		// summary:
		//		This function tests whether the item passed in is indeed an item in the store.
		// item:
		//		The item to test for being contained by the store.
		if(!this.isItem(item)){
			throw new Error(this.declaredClass + ": Invalid item argument.");
		}
	},

	_assertIsAttribute: function(/* attribute-name-string */ attribute){
		// summary:
		//		This function tests whether the item passed in is indeed a valid 'attribute' like type for the store.
		// attribute:
		//		The attribute to test for being contained by the store.
		if(typeof attribute !== "string"){
			throw new Error(this.declaredClass + ": Invalid attribute argument.");
		}
	},

	getValue: function(	/* dojo/data/api/Item */ item,
						   /* attribute-name-string */ attribute,
						   /* value? */ defaultValue){
		// summary:
		//		See dojo/data/api/Read.getValue()
		var values = this.getValues(item, attribute);
		return (values.length > 0)?values[0]:defaultValue; // mixed
	},

	getValues: function(/* dojo/data/api/Item */ item,
						/* attribute-name-string */ attribute){
		// summary:
		//		See dojo/data/api/Read.getValues()

		this._assertIsItem(item);
		this._assertIsAttribute(attribute);
		// Clone it before returning.  refs: #10474
		return (item[attribute] || []).slice(0); // Array
	},

	getAttributes: function(/* dojo/data/api/Item */ item){
		// summary:
		//		See dojo/data/api/Read.getAttributes()
		this._assertIsItem(item);
		var attributes = [];
		for(var key in item){
			// Save off only the real item attributes, not the special id marks for O(1) isItem.
			if((key !== this._storeRefPropName) && (key !== this._itemNumPropName) && (key !== this._rootItemPropName) && (key !== this._reverseRefMap)){
				attributes.push(key);
			}
		}
		return attributes; // Array
	},

	hasAttribute: function(	/* dojo/data/api/Item */ item,
							   /* attribute-name-string */ attribute){
		// summary:
		//		See dojo/data/api/Read.hasAttribute()
		this._assertIsItem(item);
		this._assertIsAttribute(attribute);
		return (attribute in item);
	},

	containsValue: function(/* dojo/data/api/Item */ item,
							/* attribute-name-string */ attribute,
							/* anything */ value){
		// summary:
		//		See dojo/data/api/Read.containsValue()
		var regexp = undefined;
		if(typeof value === "string"){
			regexp = filterUtil.patternToRegExp(value, false);
		}
		return this._containsValue(item, attribute, value, regexp); //boolean.
	},

	_containsValue: function(	/* dojo/data/api/Item */ item,
								 /* attribute-name-string */ attribute,
								 /* anything */ value,
								 /* RegExp?*/ regexp){
		// summary:
		//		Internal function for looking at the values contained by the item.
		// description:
		//		Internal function for looking at the values contained by the item.  This
		//		function allows for denoting if the comparison should be case sensitive for
		//		strings or not (for handling filtering cases where string case should not matter)
		// item:
		//		The data item to examine for attribute values.
		// attribute:
		//		The attribute to inspect.
		// value:
		//		The value to match.
		// regexp:
		//		Optional regular expression generated off value if value was of string type to handle wildcarding.
		//		If present and attribute values are string, then it can be used for comparison instead of 'value'
		return array.some(this.getValues(item, attribute), function(possibleValue){
			if(possibleValue !== null && !lang.isObject(possibleValue) && regexp){
				if(possibleValue.toString().match(regexp)){
					return true; // Boolean
				}
			}else if(value === possibleValue){
				return true; // Boolean
			}
		});
	},

	isItem: function(/* anything */ something){
		// summary:
		//		See dojo/data/api/Read.isItem()
		if(something && something[this._storeRefPropName] === this){
			if(this._arrayOfAllItems[something[this._itemNumPropName]] === something){
				return true;
			}
		}
		return false; // Boolean
	},

	isItemLoaded: function(/* anything */ something){
		// summary:
		//		See dojo/data/api/Read.isItemLoaded()
		return this.isItem(something); //boolean
	},

	loadItem: function(/* object */ keywordArgs){
		// summary:
		//		See dojo/data/api/Read.loadItem()
		this._assertIsItem(keywordArgs.item);
	},

	getFeatures: function(){
		// summary:
		//		See dojo/data/api/Read.getFeatures()
		return this._features; //Object
	},

	getLabel: function(/* dojo/data/api/Item */ item){
		// summary:
		//		See dojo/data/api/Read.getLabel()
		if(this._labelAttr && this.isItem(item)){
			return this.getValue(item,this._labelAttr); //String
		}
		return undefined; //undefined
	},

	getLabelAttributes: function(/* dojo/data/api/Item */ item){
		// summary:
		//		See dojo/data/api/Read.getLabelAttributes()
		if(this._labelAttr){
			return [this._labelAttr]; //array
		}
		return null; //null
	},

	filter: function(/* Object */ requestArgs, /* item[] */ arrayOfItems, /* Function */ findCallback){
		// summary:
		//		This method handles the basic filtering needs for ItemFile* based stores.
		var items = [],
			i, key;

		if(requestArgs.query){
			var value,
				ignoreCase = requestArgs.queryOptions ? requestArgs.queryOptions.ignoreCase : false;

			//See if there are any string values that can be regexp parsed first to avoid multiple regexp gens on the
			//same value for each item examined.  Much more efficient.
			var regexpList = {};
			for(key in requestArgs.query){
				value = requestArgs.query[key];
				if(typeof value === "string"){
					regexpList[key] = filterUtil.patternToRegExp(value, ignoreCase);
				}else if(value instanceof RegExp){
					regexpList[key] = value;
				}
			}
			for(i = 0; i < arrayOfItems.length; ++i){
				var match = true;
				var candidateItem = arrayOfItems[i];
				if(candidateItem === null){
					match = false;
				}else{
					for(key in requestArgs.query){
						value = requestArgs.query[key];
						if(!this._containsValue(candidateItem, key, value, regexpList[key])){
							match = false;
						}
					}
				}
				if(match){
					items.push(candidateItem);
				}
			}
			findCallback(items, requestArgs);
		}else{
			// We want a copy to pass back in case the parent wishes to sort the array.
			// We shouldn't allow resort of the internal list, so that multiple callers
			// can get lists and sort without affecting each other.  We also need to
			// filter out any null values that have been left as a result of deleteItem()
			// calls in ItemFileWriteStore.
			for(i = 0; i < arrayOfItems.length; ++i){
				var item = arrayOfItems[i];
				if(item !== null){
					items.push(item);
				}
			}
			findCallback(items, requestArgs);
		}
	},

	_fetchItems: function(	/* Object */ keywordArgs,
							  /* Function */ findCallback,
							  /* Function */ errorCallback){
		// summary:
		//		See dojo/data/util.simpleFetch.fetch()
		var self = this;

		if(this._loadFinished){
			this.filter(keywordArgs, this._getItemsArray(keywordArgs.queryOptions), findCallback);
		}else{
			//Do a check on the JsonFileUrl and crosscheck it.
			//If it doesn't match the cross-check, it needs to be updated
			//This allows for either url or _jsonFileUrl to he changed to
			//reset the store load location.  Done this way for backwards
			//compatibility.  People use _jsonFileUrl (even though officially
			//private.
			if(this._jsonFileUrl !== this._ccUrl){
				kernel.deprecated(this.declaredClass + ": ",
					"To change the url, set the url property of the store," +
						" not _jsonFileUrl.  _jsonFileUrl support will be removed in 2.0");
				this._ccUrl = this._jsonFileUrl;
				this.url = this._jsonFileUrl;
			}else if(this.url !== this._ccUrl){
				this._jsonFileUrl = this.url;
				this._ccUrl = this.url;
			}

			//See if there was any forced reset of data.
			if(this.data != null){
				this._jsonData = this.data;
				this.data = null;
			}

			if(this._jsonFileUrl){
				//If fetches come in before the loading has finished, but while
				//a load is in progress, we have to defer the fetching to be
				//invoked in the callback.
				if(this._loadInProgress){
					this._queuedFetches.push({args: keywordArgs, filter: lang.hitch(self, "filter"), findCallback: lang.hitch(self, findCallback)});
				}else{
					this._loadInProgress = true;
					var getArgs = {
						url: self._jsonFileUrl,
						handleAs: "json-comment-optional",
						preventCache: this.urlPreventCache,
						failOk: this.failOk
					};
					var getHandler = xhr.get(getArgs);
					getHandler.addCallback(function(data){
						try{
							self._getItemsFromLoadedData(data);
							self._loadFinished = true;
							self._loadInProgress = false;

							self.filter(keywordArgs, self._getItemsArray(keywordArgs.queryOptions), findCallback);
							self._handleQueuedFetches();
						}catch(e){
							self._loadFinished = true;
							self._loadInProgress = false;
							errorCallback(e, keywordArgs);
						}
					});
					getHandler.addErrback(function(error){
						self._loadInProgress = false;
						errorCallback(error, keywordArgs);
					});

					//Wire up the cancel to abort of the request
					//This call cancel on the deferred if it hasn't been called
					//yet and then will chain to the simple abort of the
					//simpleFetch keywordArgs
					var oldAbort = null;
					if(keywordArgs.abort){
						oldAbort = keywordArgs.abort;
					}
					keywordArgs.abort = function(){
						var df = getHandler;
						if(df && df.fired === -1){
							df.cancel();
							df = null;
						}
						if(oldAbort){
							oldAbort.call(keywordArgs);
						}
					};
				}
			}else if(this._jsonData){
				try{
					this._loadFinished = true;
					this._getItemsFromLoadedData(this._jsonData);
					this._jsonData = null;
					self.filter(keywordArgs, this._getItemsArray(keywordArgs.queryOptions), findCallback);
				}catch(e){
					errorCallback(e, keywordArgs);
				}
			}else{
				errorCallback(new Error(this.declaredClass + ": No JSON source data was provided as either URL or a nested Javascript object."), keywordArgs);
			}
		}
	},

	_handleQueuedFetches: function(){
		// summary:
		//		Internal function to execute delayed request in the store.
		
		//Execute any deferred fetches now.
		if(this._queuedFetches.length > 0){
			for(var i = 0; i < this._queuedFetches.length; i++){
				var fData = this._queuedFetches[i],
					delayedQuery = fData.args,
					delayedFilter = fData.filter,
					delayedFindCallback = fData.findCallback;
				if(delayedFilter){
					delayedFilter(delayedQuery, this._getItemsArray(delayedQuery.queryOptions), delayedFindCallback);
				}else{
					this.fetchItemByIdentity(delayedQuery);
				}
			}
			this._queuedFetches = [];
		}
	},

	_getItemsArray: function(/*object?*/queryOptions){
		// summary:
		//		Internal function to determine which list of items to search over.
		// queryOptions: The query options parameter, if any.
		if(queryOptions && queryOptions.deep){
			return this._arrayOfAllItems;
		}
		return this._arrayOfTopLevelItems;
	},

	close: function(/*dojo/data/api/Request|Object?*/ request){
		// summary:
		//		See dojo/data/api/Read.close()
		if(this.clearOnClose &&
			this._loadFinished &&
			!this._loadInProgress){
			//Reset all internalsback to default state.  This will force a reload
			//on next fetch.  This also checks that the data or url param was set
			//so that the store knows it can get data.  Without one of those being set,
			//the next fetch will trigger an error.

			if(((this._jsonFileUrl == "" || this._jsonFileUrl == null) &&
				(this.url == "" || this.url == null)
				) && this.data == null){
				console.debug(this.declaredClass + ": WARNING!  Data reload " +
					" information has not been provided." +
					"  Please set 'url' or 'data' to the appropriate value before" +
					" the next fetch");
			}
			this._arrayOfAllItems = [];
			this._arrayOfTopLevelItems = [];
			this._loadFinished = false;
			this._itemsByIdentity = null;
			this._loadInProgress = false;
			this._queuedFetches = [];
		}
	},

	_getItemsFromLoadedData: function(/* Object */ dataObject){
		// summary:
		//		Function to parse the loaded data into item format and build the internal items array.
		// description:
		//		Function to parse the loaded data into item format and build the internal items array.
		// dataObject:
		//		The JS data object containing the raw data to convery into item format.
		// returns: Array
		//		Array of items in store item format.

		// First, we define a couple little utility functions...
		var addingArrays = false,
			self = this;

		function valueIsAnItem(/* anything */ aValue){
			// summary:
			//		Given any sort of value that could be in the raw json data,
			//		return true if we should interpret the value as being an
			//		item itself, rather than a literal value or a reference.
			// example:
			// 	|	false == valueIsAnItem("Kermit");
			// 	|	false == valueIsAnItem(42);
			// 	|	false == valueIsAnItem(new Date());
			// 	|	false == valueIsAnItem({_type:'Date', _value:'1802-05-14'});
			// 	|	false == valueIsAnItem({_reference:'Kermit'});
			// 	|	true == valueIsAnItem({name:'Kermit', color:'green'});
			// 	|	true == valueIsAnItem({iggy:'pop'});
			// 	|	true == valueIsAnItem({foo:42});
			return (aValue !== null) &&
				(typeof aValue === "object") &&
				(!lang.isArray(aValue) || addingArrays) &&
				(!lang.isFunction(aValue)) &&
				(aValue.constructor == Object || lang.isArray(aValue)) &&
				(typeof aValue._reference === "undefined") &&
				(typeof aValue._type === "undefined") &&
				(typeof aValue._value === "undefined") &&
				self.hierarchical;
		}

		function addItemAndSubItemsToArrayOfAllItems(/* dojo/data/api/Item */ anItem){
			self._arrayOfAllItems.push(anItem);
			for(var attribute in anItem){
				var valueForAttribute = anItem[attribute];
				if(valueForAttribute){
					if(lang.isArray(valueForAttribute)){
						var valueArray = valueForAttribute;
						for(var k = 0; k < valueArray.length; ++k){
							var singleValue = valueArray[k];
							if(valueIsAnItem(singleValue)){
								addItemAndSubItemsToArrayOfAllItems(singleValue);
							}
						}
					}else{
						if(valueIsAnItem(valueForAttribute)){
							addItemAndSubItemsToArrayOfAllItems(valueForAttribute);
						}
					}
				}
			}
		}

		this._labelAttr = dataObject.label;

		// We need to do some transformations to convert the data structure
		// that we read from the file into a format that will be convenient
		// to work with in memory.

		// Step 1: Walk through the object hierarchy and build a list of all items
		var i,
			item;
		this._arrayOfAllItems = [];
		this._arrayOfTopLevelItems = dataObject.items;

		for(i = 0; i < this._arrayOfTopLevelItems.length; ++i){
			item = this._arrayOfTopLevelItems[i];
			if(lang.isArray(item)){
				addingArrays = true;
			}
			addItemAndSubItemsToArrayOfAllItems(item);
			item[this._rootItemPropName]=true;
		}

		// Step 2: Walk through all the attribute values of all the items,
		// and replace single values with arrays.  For example, we change this:
		//		{ name:'Miss Piggy', pets:'Foo-Foo'}
		// into this:
		//		{ name:['Miss Piggy'], pets:['Foo-Foo']}
		//
		// We also store the attribute names so we can validate our store
		// reference and item id special properties for the O(1) isItem
		var allAttributeNames = {},
			key;

		for(i = 0; i < this._arrayOfAllItems.length; ++i){
			item = this._arrayOfAllItems[i];
			for(key in item){
				if(key !== this._rootItemPropName){
					var value = item[key];
					if(value !== null){
						if(!lang.isArray(value)){
							item[key] = [value];
						}
					}else{
						item[key] = [null];
					}
				}
				allAttributeNames[key]=key;
			}
		}

		// Step 3: Build unique property names to use for the _storeRefPropName and _itemNumPropName
		// This should go really fast, it will generally never even run the loop.
		while(allAttributeNames[this._storeRefPropName]){
			this._storeRefPropName += "_";
		}
		while(allAttributeNames[this._itemNumPropName]){
			this._itemNumPropName += "_";
		}
		while(allAttributeNames[this._reverseRefMap]){
			this._reverseRefMap += "_";
		}

		// Step 4: Some data files specify an optional 'identifier', which is
		// the name of an attribute that holds the identity of each item.
		// If this data file specified an identifier attribute, then build a
		// hash table of items keyed by the identity of the items.
		var arrayOfValues;

		var identifier = dataObject.identifier;
		if(identifier){
			this._itemsByIdentity = {};
			this._features['dojo.data.api.Identity'] = identifier;
			for(i = 0; i < this._arrayOfAllItems.length; ++i){
				item = this._arrayOfAllItems[i];
				arrayOfValues = item[identifier];
				var identity = arrayOfValues[0];
				if(!Object.hasOwnProperty.call(this._itemsByIdentity, identity)){
					this._itemsByIdentity[identity] = item;
				}else{
					if(this._jsonFileUrl){
						throw new Error(this.declaredClass + ":  The json data as specified by: [" + this._jsonFileUrl + "] is malformed.  Items within the list have identifier: [" + identifier + "].  Value collided: [" + identity + "]");
					}else if(this._jsonData){
						throw new Error(this.declaredClass + ":  The json data provided by the creation arguments is malformed.  Items within the list have identifier: [" + identifier + "].  Value collided: [" + identity + "]");
					}
				}
			}
		}else{
			this._features['dojo.data.api.Identity'] = Number;
		}

		// Step 5: Walk through all the items, and set each item's properties
		// for _storeRefPropName and _itemNumPropName, so that store.isItem() will return true.
		for(i = 0; i < this._arrayOfAllItems.length; ++i){
			item = this._arrayOfAllItems[i];
			item[this._storeRefPropName] = this;
			item[this._itemNumPropName] = i;
		}

		// Step 6: We walk through all the attribute values of all the items,
		// looking for type/value literals and item-references.
		//
		// We replace item-references with pointers to items.  For example, we change:
		//		{ name:['Kermit'], friends:[{_reference:{name:'Miss Piggy'}}] }
		// into this:
		//		{ name:['Kermit'], friends:[miss_piggy] }
		// (where miss_piggy is the object representing the 'Miss Piggy' item).
		//
		// We replace type/value pairs with typed-literals.  For example, we change:
		//		{ name:['Nelson Mandela'], born:[{_type:'Date', _value:'1918-07-18'}] }
		// into this:
		//		{ name:['Kermit'], born:(new Date(1918, 6, 18)) }
		//
		// We also generate the associate map for all items for the O(1) isItem function.
		for(i = 0; i < this._arrayOfAllItems.length; ++i){
			item = this._arrayOfAllItems[i]; // example: { name:['Kermit'], friends:[{_reference:{name:'Miss Piggy'}}] }
			for(key in item){
				arrayOfValues = item[key]; // example: [{_reference:{name:'Miss Piggy'}}]
				for(var j = 0; j < arrayOfValues.length; ++j){
					value = arrayOfValues[j]; // example: {_reference:{name:'Miss Piggy'}}
					if(value !== null && typeof value == "object"){
						if(("_type" in value) && ("_value" in value)){
							var type = value._type; // examples: 'Date', 'Color', or 'ComplexNumber'
							var mappingObj = this._datatypeMap[type]; // examples: Date, dojo.Color, foo.math.ComplexNumber, {type: dojo.Color, deserialize(value){ return new dojo.Color(value)}}
							if(!mappingObj){
								throw new Error("dojo.data.ItemFileReadStore: in the typeMap constructor arg, no object class was specified for the datatype '" + type + "'");
							}else if(lang.isFunction(mappingObj)){
								arrayOfValues[j] = new mappingObj(value._value);
							}else if(lang.isFunction(mappingObj.deserialize)){
								arrayOfValues[j] = mappingObj.deserialize(value._value);
							}else{
								throw new Error("dojo.data.ItemFileReadStore: Value provided in typeMap was neither a constructor, nor a an object with a deserialize function");
							}
						}
						if(value._reference){
							var referenceDescription = value._reference; // example: {name:'Miss Piggy'}
							if(!lang.isObject(referenceDescription)){
								// example: 'Miss Piggy'
								// from an item like: { name:['Kermit'], friends:[{_reference:'Miss Piggy'}]}
								arrayOfValues[j] = this._getItemByIdentity(referenceDescription);
							}else{
								// example: {name:'Miss Piggy'}
								// from an item like: { name:['Kermit'], friends:[{_reference:{name:'Miss Piggy'}}] }
								for(var k = 0; k < this._arrayOfAllItems.length; ++k){
									var candidateItem = this._arrayOfAllItems[k],
										found = true;
									for(var refKey in referenceDescription){
										if(candidateItem[refKey] != referenceDescription[refKey]){
											found = false;
										}
									}
									if(found){
										arrayOfValues[j] = candidateItem;
									}
								}
							}
							if(this.referenceIntegrity){
								var refItem = arrayOfValues[j];
								if(this.isItem(refItem)){
									this._addReferenceToMap(refItem, item, key);
								}
							}
						}else if(this.isItem(value)){
							//It's a child item (not one referenced through _reference).
							//We need to treat this as a referenced item, so it can be cleaned up
							//in a write store easily.
							if(this.referenceIntegrity){
								this._addReferenceToMap(value, item, key);
							}
						}
					}
				}
			}
		}
	},

	_addReferenceToMap: function(/*item*/ refItem, /*item*/ parentItem, /*string*/ attribute){
		// summary:
		//		Method to add an reference map entry for an item and attribute.
		// description:
		//		Method to add an reference map entry for an item and attribute.
		// refItem:
		//		The item that is referenced.
		// parentItem:
		//		The item that holds the new reference to refItem.
		// attribute:
		//		The attribute on parentItem that contains the new reference.

		//Stub function, does nothing.  Real processing is in ItemFileWriteStore.
	},

	getIdentity: function(/* dojo/data/api/Item */ item){
		// summary:
		//		See dojo/data/api/Identity.getIdentity()
		var identifier = this._features['dojo.data.api.Identity'];
		if(identifier === Number){
			return item[this._itemNumPropName]; // Number
		}else{
			var arrayOfValues = item[identifier];
			if(arrayOfValues){
				return arrayOfValues[0]; // Object|String
			}
		}
		return null; // null
	},

	fetchItemByIdentity: function(/* Object */ keywordArgs){
		// summary:
		//		See dojo/data/api/Identity.fetchItemByIdentity()

		// Hasn't loaded yet, we have to trigger the load.
		var item,
			scope;
		if(!this._loadFinished){
			var self = this;
			//Do a check on the JsonFileUrl and crosscheck it.
			//If it doesn't match the cross-check, it needs to be updated
			//This allows for either url or _jsonFileUrl to he changed to
			//reset the store load location.  Done this way for backwards
			//compatibility.  People use _jsonFileUrl (even though officially
			//private.
			if(this._jsonFileUrl !== this._ccUrl){
				kernel.deprecated(this.declaredClass + ": ",
					"To change the url, set the url property of the store," +
						" not _jsonFileUrl.  _jsonFileUrl support will be removed in 2.0");
				this._ccUrl = this._jsonFileUrl;
				this.url = this._jsonFileUrl;
			}else if(this.url !== this._ccUrl){
				this._jsonFileUrl = this.url;
				this._ccUrl = this.url;
			}

			//See if there was any forced reset of data.
			if(this.data != null && this._jsonData == null){
				this._jsonData = this.data;
				this.data = null;
			}

			if(this._jsonFileUrl){

				if(this._loadInProgress){
					this._queuedFetches.push({args: keywordArgs});
				}else{
					this._loadInProgress = true;
					var getArgs = {
						url: self._jsonFileUrl,
						handleAs: "json-comment-optional",
						preventCache: this.urlPreventCache,
						failOk: this.failOk
					};
					var getHandler = xhr.get(getArgs);
					getHandler.addCallback(function(data){
						var scope = keywordArgs.scope?keywordArgs.scope:kernel.global;
						try{
							self._getItemsFromLoadedData(data);
							self._loadFinished = true;
							self._loadInProgress = false;
							item = self._getItemByIdentity(keywordArgs.identity);
							if(keywordArgs.onItem){
								keywordArgs.onItem.call(scope, item);
							}
							self._handleQueuedFetches();
						}catch(error){
							self._loadInProgress = false;
							if(keywordArgs.onError){
								keywordArgs.onError.call(scope, error);
							}
						}
					});
					getHandler.addErrback(function(error){
						self._loadInProgress = false;
						if(keywordArgs.onError){
							var scope = keywordArgs.scope?keywordArgs.scope:kernel.global;
							keywordArgs.onError.call(scope, error);
						}
					});
				}

			}else if(this._jsonData){
				// Passed in data, no need to xhr.
				self._getItemsFromLoadedData(self._jsonData);
				self._jsonData = null;
				self._loadFinished = true;
				item = self._getItemByIdentity(keywordArgs.identity);
				if(keywordArgs.onItem){
					scope = keywordArgs.scope?keywordArgs.scope:kernel.global;
					keywordArgs.onItem.call(scope, item);
				}
			}
		}else{
			// Already loaded.  We can just look it up and call back.
			item = this._getItemByIdentity(keywordArgs.identity);
			if(keywordArgs.onItem){
				scope = keywordArgs.scope?keywordArgs.scope:kernel.global;
				keywordArgs.onItem.call(scope, item);
			}
		}
	},

	_getItemByIdentity: function(/* Object */ identity){
		// summary:
		//		Internal function to look an item up by its identity map.
		var item = null;
		if(this._itemsByIdentity){
			// If this map is defined, we need to just try to get it.  If it fails
			// the item does not exist.
			if(Object.hasOwnProperty.call(this._itemsByIdentity, identity)){
				item = this._itemsByIdentity[identity];
			}
		}else if (Object.hasOwnProperty.call(this._arrayOfAllItems, identity)){
			item = this._arrayOfAllItems[identity];
		}
		if(item === undefined){
			item = null;
		}
		return item; // Object
	},

	getIdentityAttributes: function(/* dojo/data/api/Item */ item){
		// summary:
		//		See dojo/data/api/Identity.getIdentityAttributes()

		var identifier = this._features['dojo.data.api.Identity'];
		if(identifier === Number){
			// If (identifier === Number) it means getIdentity() just returns
			// an integer item-number for each item.  The dojo/data/api/Identity
			// spec says we need to return null if the identity is not composed
			// of attributes
			return null; // null
		}else{
			return [identifier]; // Array
		}
	},

	_forceLoad: function(){
		// summary:
		//		Internal function to force a load of the store if it hasn't occurred yet.  This is required
		//		for specific functions to work properly.
		var self = this;
		//Do a check on the JsonFileUrl and crosscheck it.
		//If it doesn't match the cross-check, it needs to be updated
		//This allows for either url or _jsonFileUrl to he changed to
		//reset the store load location.  Done this way for backwards
		//compatibility.  People use _jsonFileUrl (even though officially
		//private.
		if(this._jsonFileUrl !== this._ccUrl){
			kernel.deprecated(this.declaredClass + ": ",
				"To change the url, set the url property of the store," +
					" not _jsonFileUrl.  _jsonFileUrl support will be removed in 2.0");
			this._ccUrl = this._jsonFileUrl;
			this.url = this._jsonFileUrl;
		}else if(this.url !== this._ccUrl){
			this._jsonFileUrl = this.url;
			this._ccUrl = this.url;
		}

		//See if there was any forced reset of data.
		if(this.data != null){
			this._jsonData = this.data;
			this.data = null;
		}

		if(this._jsonFileUrl){
			var getArgs = {
				url: this._jsonFileUrl,
				handleAs: "json-comment-optional",
				preventCache: this.urlPreventCache,
				failOk: this.failOk,
				sync: true
			};
			var getHandler = xhr.get(getArgs);
			getHandler.addCallback(function(data){
				try{
					//Check to be sure there wasn't another load going on concurrently
					//So we don't clobber data that comes in on it.  If there is a load going on
					//then do not save this data.  It will potentially clobber current data.
					//We mainly wanted to sync/wait here.
					//TODO:  Revisit the loading scheme of this store to improve multi-initial
					//request handling.
					if(self._loadInProgress !== true && !self._loadFinished){
						self._getItemsFromLoadedData(data);
						self._loadFinished = true;
					}else if(self._loadInProgress){
						//Okay, we hit an error state we can't recover from.  A forced load occurred
						//while an async load was occurring.  Since we cannot block at this point, the best
						//that can be managed is to throw an error.
						throw new Error(this.declaredClass + ":  Unable to perform a synchronous load, an async load is in progress.");
					}
				}catch(e){
					console.log(e);
					throw e;
				}
			});
			getHandler.addErrback(function(error){
				throw error;
			});
		}else if(this._jsonData){
			self._getItemsFromLoadedData(self._jsonData);
			self._jsonData = null;
			self._loadFinished = true;
		}
	}
});
//Mix in the simple fetch implementation to this class.
lang.extend(ItemFileReadStore,simpleFetch);

return ItemFileReadStore;

});

},
'dijit/form/DateTextBox':function(){
define("dijit/form/DateTextBox", [
	"dojo/_base/declare", // declare
	"../Calendar",
	"./_DateTimeTextBox"
], function(declare, Calendar, _DateTimeTextBox){

	// module:
	//		dijit/form/DateTextBox

	return declare("dijit.form.DateTextBox", _DateTimeTextBox, {
		// summary:
		//		A validating, serializable, range-bound date text box with a drop down calendar
		// example:
		// |	new DateTextBox({value: new Date(2009, 0, 20)})
		// example:
		// |	<input data-dojo-type='dijit/form/DateTextBox' value='2009-01-20'>

		baseClass: "dijitTextBox dijitComboBox dijitDateTextBox",
		popupClass: Calendar,
		_selector: "date",

		// value: Date
		//		The value of this widget as a JavaScript Date object, with only year/month/day specified.
		//		If specified in markup, use the format specified in `stamp.fromISOString`.
		//		set("value", ...) accepts either a Date object or a string.
		value: new Date("")	// value.toString()="NaN"
	});
});

},
'dijit/Calendar':function(){
define("dijit/Calendar", [
	"dojo/_base/array", // array.map
	"dojo/date",
	"dojo/date/locale",
	"dojo/_base/declare", // declare
	"dojo/dom-attr", // domAttr.get
	"dojo/dom-class", // domClass.add domClass.contains domClass.remove domClass.toggle
	"dojo/_base/event", // event.stop
	"dojo/_base/kernel", // kernel.deprecated
	"dojo/keys", // keys
	"dojo/_base/lang", // lang.hitch
	"dojo/sniff", // has("ie")
	"./CalendarLite",
	"./_Widget",
	"./_CssStateMixin",
	"./_TemplatedMixin",
	"./form/DropDownButton"
], function(array, date, local, declare, domAttr, domClass, event, kernel, keys, lang, has,
			CalendarLite, _Widget, _CssStateMixin, _TemplatedMixin, DropDownButton){

	// module:
	//		dijit/Calendar

	var Calendar = declare("dijit.Calendar",
		[CalendarLite, _Widget, _CssStateMixin], // _Widget for deprecated methods like setAttribute()
		{
		// summary:
		//		A simple GUI for choosing a date in the context of a monthly calendar.
		//
		// description:
		//		See CalendarLite for general description.   Calendar extends CalendarLite, adding:
		//
		//		- month drop down list
		//		- keyboard navigation
		//		- CSS classes for hover/mousepress on date, month, and year nodes
		//		- support of deprecated methods (will be removed in 2.0)

		// Set node classes for various mouse events, see dijit._CssStateMixin for more details
		cssStateNodes: {
			"decrementMonth": "dijitCalendarArrow",
			"incrementMonth": "dijitCalendarArrow",
			"previousYearLabelNode": "dijitCalendarPreviousYear",
			"nextYearLabelNode": "dijitCalendarNextYear"
		},

		setValue: function(/*Date*/ value){
			// summary:
			//		Deprecated.   Use set('value', ...) instead.
			// tags:
			//		deprecated
			kernel.deprecated("dijit.Calendar:setValue() is deprecated.  Use set('value', ...) instead.", "", "2.0");
			this.set('value', value);
		},

		_createMonthWidget: function(){
			// summary:
			//		Creates the drop down button that displays the current month and lets user pick a new one

			return new Calendar._MonthDropDownButton({
				id: this.id + "_mddb",
				tabIndex: -1,
				onMonthSelect: lang.hitch(this, "_onMonthSelect"),
				lang: this.lang,
				dateLocaleModule: this.dateLocaleModule
			}, this.monthNode);
		},

		postCreate: function(){
			this.inherited(arguments);

			// Events specific to Calendar, not used in CalendarLite
			this.connect(this.domNode, "onkeydown", "_onKeyDown");
			this.connect(this.dateRowsNode, "onmouseover", "_onDayMouseOver");
			this.connect(this.dateRowsNode, "onmouseout", "_onDayMouseOut");
			this.connect(this.dateRowsNode, "onmousedown", "_onDayMouseDown");
			this.connect(this.dateRowsNode, "onmouseup", "_onDayMouseUp");
		},

		_onMonthSelect: function(/*Number*/ newMonth){
			// summary:
			//		Handler for when user selects a month from the drop down list
			// tags:
			//		protected

			// move to selected month, bounding by the number of days in the month
			// (ex: jan 31 --> feb 28, not feb 31)
			var date  = new this.dateClassObj(this.currentFocus);
			date.setDate(1);
			date.setMonth(newMonth);
			var daysInMonth = this.dateModule.getDaysInMonth(date);
			var currentDate = this.currentFocus.getDate();
			date.setDate(Math.min(currentDate, daysInMonth));
			this._setCurrentFocusAttr(date);
		},

		_onDayMouseOver: function(/*Event*/ evt){
			// summary:
			//		Handler for mouse over events on days, sets hovered style
			// tags:
			//		protected

			// event can occur on <td> or the <span> inside the td,
			// set node to the <td>.
			var node =
				domClass.contains(evt.target, "dijitCalendarDateLabel") ?
				evt.target.parentNode :
				evt.target;

			if(node && (
				(node.dijitDateValue && !domClass.contains(node, "dijitCalendarDisabledDate"))
					|| node == this.previousYearLabelNode || node == this.nextYearLabelNode
				)){
				domClass.add(node, "dijitCalendarHoveredDate");
				this._currentNode = node;
			}
		},

		_onDayMouseOut: function(/*Event*/ evt){
			// summary:
			//		Handler for mouse out events on days, clears hovered style
			// tags:
			//		protected

			if(!this._currentNode){ return; }

			// if mouse out occurs moving from <td> to <span> inside <td>, ignore it
			if(evt.relatedTarget && evt.relatedTarget.parentNode == this._currentNode){ return; }
			var cls = "dijitCalendarHoveredDate";
			if(domClass.contains(this._currentNode, "dijitCalendarActiveDate")){
				cls += " dijitCalendarActiveDate";
			}
			domClass.remove(this._currentNode, cls);
			this._currentNode = null;
		},

		_onDayMouseDown: function(/*Event*/ evt){
			var node = evt.target.parentNode;
			if(node && node.dijitDateValue && !domClass.contains(node, "dijitCalendarDisabledDate")){
				domClass.add(node, "dijitCalendarActiveDate");
				this._currentNode = node;
			}
		},

		_onDayMouseUp: function(/*Event*/ evt){
			var node = evt.target.parentNode;
			if(node && node.dijitDateValue){
				domClass.remove(node, "dijitCalendarActiveDate");
			}
		},

		handleKey: function(/*Event*/ evt){
			// summary:
			//		Provides keyboard navigation of calendar.
			// description:
			//		Called from _onKeyDown() to handle keypress on a stand alone Calendar,
			//		and also from `dijit/form/_DateTimeTextBox` to pass a keydown event
			//		from the `dijit/form/DateTextBox` to be handled in this widget
			// returns:
			//		False if the key was recognized as a navigation key,
			//		to indicate that the event was handled by Calendar and shouldn't be propagated
			// tags:
			//		protected
			var increment = -1,
				interval,
				newValue = this.currentFocus;
			switch(evt.keyCode){
				case keys.RIGHT_ARROW:
					increment = 1;
					//fallthrough...
				case keys.LEFT_ARROW:
					interval = "day";
					if(!this.isLeftToRight()){ increment *= -1; }
					break;
				case keys.DOWN_ARROW:
					increment = 1;
					//fallthrough...
				case keys.UP_ARROW:
					interval = "week";
					break;
				case keys.PAGE_DOWN:
					increment = 1;
					//fallthrough...
				case keys.PAGE_UP:
					interval = evt.ctrlKey || evt.altKey ? "year" : "month";
					break;
				case keys.END:
					// go to the next month
					newValue = this.dateModule.add(newValue, "month", 1);
					// subtract a day from the result when we're done
					interval = "day";
					//fallthrough...
				case keys.HOME:
					newValue = new this.dateClassObj(newValue);
					newValue.setDate(1);
					break;
				case keys.ENTER:
				case keys.SPACE:
					this.set("value", this.currentFocus);
					break;
				default:
					return true;
			}

			if(interval){
				newValue = this.dateModule.add(newValue, interval, increment);
			}

			this._setCurrentFocusAttr(newValue);

			return false;
		},

		_onKeyDown: function(/*Event*/ evt){
			// summary:
			//		For handling keypress events on a stand alone calendar
			if(!this.handleKey(evt)){
				event.stop(evt);
			}
		},

		onValueSelected: function(/*Date*/ /*===== date =====*/){
			// summary:
			//		Deprecated.   Notification that a date cell was selected.  It may be the same as the previous value.
			// description:
			//		Formerly used by `dijit/form/_DateTimeTextBox` (and thus `dijit/form/DateTextBox`)
			//		to get notification when the user has clicked a date.  Now onExecute() (above) is used.
			// tags:
			//		protected
		},

		onChange: function(value){
			this.onValueSelected(value);	// remove in 2.0
		},

		getClassForDate: function(/*===== dateObject, locale =====*/){
			// summary:
			//		May be overridden to return CSS classes to associate with the date entry for the given dateObject,
			//		for example to indicate a holiday in specified locale.
			// dateObject: Date
			// locale: String?
			// tags:
			//		extension

/*=====
			return ""; // String
=====*/
		}
	});

	Calendar._MonthDropDownButton = declare("dijit.Calendar._MonthDropDownButton", DropDownButton, {
		// summary:
		//		DropDownButton for the current month.    Displays name of current month
		//		and a list of month names in the drop down

		onMonthSelect: function(){ },

		postCreate: function(){
			this.inherited(arguments);
			this.dropDown = new Calendar._MonthDropDown({
				id: this.id + "_mdd", //do not change this id because it is referenced in the template
				onChange: this.onMonthSelect
			});
		},
		_setMonthAttr: function(month){
			// summary:
			//		Set the current month to display as a label
			var monthNames = this.dateLocaleModule.getNames('months', 'wide', 'standAlone', this.lang, month);
			this.dropDown.set("months", monthNames);

			// Set name of current month and also fill in spacer element with all the month names
			// (invisible) so that the maximum width will affect layout.   But not on IE6 because then
			// the center <TH> overlaps the right <TH> (due to a browser bug).
			this.containerNode.innerHTML =
				(has("ie") == 6 ? "" : "<div class='dijitSpacer'>" + this.dropDown.domNode.innerHTML + "</div>") +
				"<div class='dijitCalendarMonthLabel dijitCalendarCurrentMonthLabel'>" +  monthNames[month.getMonth()] + "</div>";
		}
	});

	Calendar._MonthDropDown = declare("dijit.Calendar._MonthDropDown", [_Widget, _TemplatedMixin], {
		// summary:
		//		The list-of-months drop down from the MonthDropDownButton

		// months: String[]
		//		List of names of months, possibly w/some undefined entries for Hebrew leap months
		//		(ex: ["January", "February", undefined, "April", ...])
		months: [],

		templateString: "<div class='dijitCalendarMonthMenu dijitMenu' " +
			"data-dojo-attach-event='onclick:_onClick,onmouseover:_onMenuHover,onmouseout:_onMenuHover'></div>",

		_setMonthsAttr: function(/*String[]*/ months){
			this.domNode.innerHTML = array.map(months, function(month, idx){
					return month ? "<div class='dijitCalendarMonthLabel' month='" + idx +"'>" + month + "</div>" : "";
				}).join("");
		},

		_onClick: function(/*Event*/ evt){
			this.onChange(domAttr.get(evt.target, "month"));
		},

		onChange: function(/*Number*/ /*===== month =====*/){
			// summary:
			//		Callback when month is selected from drop down
		},

		_onMenuHover: function(evt){
			domClass.toggle(evt.target, "dijitCalendarMonthLabelHover", evt.type == "mouseover");
		}
	});

	return Calendar;
});

},
'dojo/date':function(){
define("dojo/date", ["./has", "./_base/lang"], function(has, lang){
// module:
//		dojo/date

var date = {
	// summary:
	//		Date manipulation utilities
};

date.getDaysInMonth = function(/*Date*/dateObject){
	// summary:
	//		Returns the number of days in the month used by dateObject
	var month = dateObject.getMonth();
	var days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
	if(month == 1 && date.isLeapYear(dateObject)){ return 29; } // Number
	return days[month]; // Number
};

date.isLeapYear = function(/*Date*/dateObject){
	// summary:
	//		Determines if the year of the dateObject is a leap year
	// description:
	//		Leap years are years with an additional day YYYY-02-29, where the
	//		year number is a multiple of four with the following exception: If
	//		a year is a multiple of 100, then it is only a leap year if it is
	//		also a multiple of 400. For example, 1900 was not a leap year, but
	//		2000 is one.

	var year = dateObject.getFullYear();
	return !(year%400) || (!(year%4) && !!(year%100)); // Boolean
};

// FIXME: This is not localized
date.getTimezoneName = function(/*Date*/dateObject){
	// summary:
	//		Get the user's time zone as provided by the browser
	// dateObject:
	//		Needed because the timezone may vary with time (daylight savings)
	// description:
	//		Try to get time zone info from toString or toLocaleString method of
	//		the Date object -- UTC offset is not a time zone.  See
	//		http://www.twinsun.com/tz/tz-link.htm Note: results may be
	//		inconsistent across browsers.

	var str = dateObject.toString(); // Start looking in toString
	var tz = ''; // The result -- return empty string if nothing found
	var match;

	// First look for something in parentheses -- fast lookup, no regex
	var pos = str.indexOf('(');
	if(pos > -1){
		tz = str.substring(++pos, str.indexOf(')'));
	}else{
		// If at first you don't succeed ...
		// If IE knows about the TZ, it appears before the year
		// Capital letters or slash before a 4-digit year
		// at the end of string
		var pat = /([A-Z\/]+) \d{4}$/;
		if((match = str.match(pat))){
			tz = match[1];
		}else{
		// Some browsers (e.g. Safari) glue the TZ on the end
		// of toLocaleString instead of putting it in toString
			str = dateObject.toLocaleString();
			// Capital letters or slash -- end of string,
			// after space
			pat = / ([A-Z\/]+)$/;
			if((match = str.match(pat))){
				tz = match[1];
			}
		}
	}

	// Make sure it doesn't somehow end up return AM or PM
	return (tz == 'AM' || tz == 'PM') ? '' : tz; // String
};

// Utility methods to do arithmetic calculations with Dates

date.compare = function(/*Date*/date1, /*Date?*/date2, /*String?*/portion){
	// summary:
	//		Compare two date objects by date, time, or both.
	// description:
	//		Returns 0 if equal, positive if a > b, else negative.
	// date1:
	//		Date object
	// date2:
	//		Date object.  If not specified, the current Date is used.
	// portion:
	//		A string indicating the "date" or "time" portion of a Date object.
	//		Compares both "date" and "time" by default.  One of the following:
	//		"date", "time", "datetime"

	// Extra step required in copy for IE - see #3112
	date1 = new Date(+date1);
	date2 = new Date(+(date2 || new Date()));

	if(portion == "date"){
		// Ignore times and compare dates.
		date1.setHours(0, 0, 0, 0);
		date2.setHours(0, 0, 0, 0);
	}else if(portion == "time"){
		// Ignore dates and compare times.
		date1.setFullYear(0, 0, 0);
		date2.setFullYear(0, 0, 0);
	}

	if(date1 > date2){ return 1; } // int
	if(date1 < date2){ return -1; } // int
	return 0; // int
};

date.add = function(/*Date*/date, /*String*/interval, /*int*/amount){
	// summary:
	//		Add to a Date in intervals of different size, from milliseconds to years
	// date: Date
	//		Date object to start with
	// interval:
	//		A string representing the interval.  One of the following:
	//		"year", "month", "day", "hour", "minute", "second",
	//		"millisecond", "quarter", "week", "weekday"
	// amount:
	//		How much to add to the date.

	var sum = new Date(+date); // convert to Number before copying to accomodate IE (#3112)
	var fixOvershoot = false;
	var property = "Date";

	switch(interval){
		case "day":
			break;
		case "weekday":
			//i18n FIXME: assumes Saturday/Sunday weekend, but this is not always true.  see dojo/cldr/supplemental

			// Divide the increment time span into weekspans plus leftover days
			// e.g., 8 days is one 5-day weekspan / and two leftover days
			// Can't have zero leftover days, so numbers divisible by 5 get
			// a days value of 5, and the remaining days make up the number of weeks
			var days, weeks;
			var mod = amount % 5;
			if(!mod){
				days = (amount > 0) ? 5 : -5;
				weeks = (amount > 0) ? ((amount-5)/5) : ((amount+5)/5);
			}else{
				days = mod;
				weeks = parseInt(amount/5);
			}
			// Get weekday value for orig date param
			var strt = date.getDay();
			// Orig date is Sat / positive incrementer
			// Jump over Sun
			var adj = 0;
			if(strt == 6 && amount > 0){
				adj = 1;
			}else if(strt == 0 && amount < 0){
			// Orig date is Sun / negative incrementer
			// Jump back over Sat
				adj = -1;
			}
			// Get weekday val for the new date
			var trgt = strt + days;
			// New date is on Sat or Sun
			if(trgt == 0 || trgt == 6){
				adj = (amount > 0) ? 2 : -2;
			}
			// Increment by number of weeks plus leftover days plus
			// weekend adjustments
			amount = (7 * weeks) + days + adj;
			break;
		case "year":
			property = "FullYear";
			// Keep increment/decrement from 2/29 out of March
			fixOvershoot = true;
			break;
		case "week":
			amount *= 7;
			break;
		case "quarter":
			// Naive quarter is just three months
			amount *= 3;
			// fallthrough...
		case "month":
			// Reset to last day of month if you overshoot
			fixOvershoot = true;
			property = "Month";
			break;
//		case "hour":
//		case "minute":
//		case "second":
//		case "millisecond":
		default:
			property = "UTC"+interval.charAt(0).toUpperCase() + interval.substring(1) + "s";
	}

	if(property){
		sum["set"+property](sum["get"+property]()+amount);
	}

	if(fixOvershoot && (sum.getDate() < date.getDate())){
		sum.setDate(0);
	}

	return sum; // Date
};

date.difference = function(/*Date*/date1, /*Date?*/date2, /*String?*/interval){
	// summary:
	//		Get the difference in a specific unit of time (e.g., number of
	//		months, weeks, days, etc.) between two dates, rounded to the
	//		nearest integer.
	// date1:
	//		Date object
	// date2:
	//		Date object.  If not specified, the current Date is used.
	// interval:
	//		A string representing the interval.  One of the following:
	//		"year", "month", "day", "hour", "minute", "second",
	//		"millisecond", "quarter", "week", "weekday"
	//
	//		Defaults to "day".

	date2 = date2 || new Date();
	interval = interval || "day";
	var yearDiff = date2.getFullYear() - date1.getFullYear();
	var delta = 1; // Integer return value

	switch(interval){
		case "quarter":
			var m1 = date1.getMonth();
			var m2 = date2.getMonth();
			// Figure out which quarter the months are in
			var q1 = Math.floor(m1/3) + 1;
			var q2 = Math.floor(m2/3) + 1;
			// Add quarters for any year difference between the dates
			q2 += (yearDiff * 4);
			delta = q2 - q1;
			break;
		case "weekday":
			var days = Math.round(date.difference(date1, date2, "day"));
			var weeks = parseInt(date.difference(date1, date2, "week"));
			var mod = days % 7;

			// Even number of weeks
			if(mod == 0){
				days = weeks*5;
			}else{
				// Weeks plus spare change (< 7 days)
				var adj = 0;
				var aDay = date1.getDay();
				var bDay = date2.getDay();

				weeks = parseInt(days/7);
				mod = days % 7;
				// Mark the date advanced by the number of
				// round weeks (may be zero)
				var dtMark = new Date(date1);
				dtMark.setDate(dtMark.getDate()+(weeks*7));
				var dayMark = dtMark.getDay();

				// Spare change days -- 6 or less
				if(days > 0){
					switch(true){
						// Range starts on Sat
						case aDay == 6:
							adj = -1;
							break;
						// Range starts on Sun
						case aDay == 0:
							adj = 0;
							break;
						// Range ends on Sat
						case bDay == 6:
							adj = -1;
							break;
						// Range ends on Sun
						case bDay == 0:
							adj = -2;
							break;
						// Range contains weekend
						case (dayMark + mod) > 5:
							adj = -2;
					}
				}else if(days < 0){
					switch(true){
						// Range starts on Sat
						case aDay == 6:
							adj = 0;
							break;
						// Range starts on Sun
						case aDay == 0:
							adj = 1;
							break;
						// Range ends on Sat
						case bDay == 6:
							adj = 2;
							break;
						// Range ends on Sun
						case bDay == 0:
							adj = 1;
							break;
						// Range contains weekend
						case (dayMark + mod) < 0:
							adj = 2;
					}
				}
				days += adj;
				days -= (weeks*2);
			}
			delta = days;
			break;
		case "year":
			delta = yearDiff;
			break;
		case "month":
			delta = (date2.getMonth() - date1.getMonth()) + (yearDiff * 12);
			break;
		case "week":
			// Truncate instead of rounding
			// Don't use Math.floor -- value may be negative
			delta = parseInt(date.difference(date1, date2, "day")/7);
			break;
		case "day":
			delta /= 24;
			// fallthrough
		case "hour":
			delta /= 60;
			// fallthrough
		case "minute":
			delta /= 60;
			// fallthrough
		case "second":
			delta /= 1000;
			// fallthrough
		case "millisecond":
			delta *= date2.getTime() - date1.getTime();
	}

	// Round for fractional values and DST leaps
	return Math.round(delta); // Number (integer)
};

// Don't use setObject() because it may overwrite dojo/date/stamp (if that has already been loaded)
 1  && lang.mixin(lang.getObject("dojo.date", true), date);

return date;
});

},
'dojo/date/locale':function(){
define("dojo/date/locale", [
	"../_base/lang",
	"../_base/array",
	"../date",
	/*===== "../_base/declare", =====*/
	"../cldr/supplemental",
	"../i18n",
	"../regexp",
	"../string",
	"../i18n!../cldr/nls/gregorian",
	"module"
], function(lang, array, date, /*===== declare, =====*/ supplemental, i18n, regexp, string, gregorian, module){

// module:
//		dojo/date/locale

var exports = {
	// summary:
	//		This modules defines dojo/date/locale, localization methods for Date.
};
lang.setObject(module.id.replace(/\//g, "."), exports);

// Localization methods for Date.   Honor local customs using locale-dependent dojo.cldr data.

// Load the bundles containing localization information for
// names and formats

//NOTE: Everything in this module assumes Gregorian calendars.
// Other calendars will be implemented in separate modules.

	// Format a pattern without literals
	function formatPattern(dateObject, bundle, options, pattern){
		return pattern.replace(/([a-z])\1*/ig, function(match){
			var s, pad,
				c = match.charAt(0),
				l = match.length,
				widthList = ["abbr", "wide", "narrow"];
			switch(c){
				case 'G':
					s = bundle[(l < 4) ? "eraAbbr" : "eraNames"][dateObject.getFullYear() < 0 ? 0 : 1];
					break;
				case 'y':
					s = dateObject.getFullYear();
					switch(l){
						case 1:
							break;
						case 2:
							if(!options.fullYear){
								s = String(s); s = s.substr(s.length - 2);
								break;
							}
							// fallthrough
						default:
							pad = true;
					}
					break;
				case 'Q':
				case 'q':
					s = Math.ceil((dateObject.getMonth()+1)/3);
//					switch(l){
//						case 1: case 2:
							pad = true;
//							break;
//						case 3: case 4: // unimplemented
//					}
					break;
				case 'M':
				case 'L':
					var m = dateObject.getMonth();
					if(l<3){
						s = m+1; pad = true;
					}else{
						var propM = [
							"months",
							c == 'L' ? "standAlone" : "format",
							widthList[l-3]
						].join("-");
						s = bundle[propM][m];
					}
					break;
				case 'w':
					var firstDay = 0;
					s = exports._getWeekOfYear(dateObject, firstDay); pad = true;
					break;
				case 'd':
					s = dateObject.getDate(); pad = true;
					break;
				case 'D':
					s = exports._getDayOfYear(dateObject); pad = true;
					break;
				case 'e':
				case 'c':
					var d = dateObject.getDay();
					if(l<2){
						s = (d - supplemental.getFirstDayOfWeek(options.locale) + 8) % 7
						break;
					}
					// fallthrough
				case 'E':
					d = dateObject.getDay();
					if(l<3){
						s = d+1; pad = true;
					}else{
						var propD = [
							"days",
							c == 'c' ? "standAlone" : "format",
							widthList[l-3]
						].join("-");
						s = bundle[propD][d];
					}
					break;
				case 'a':
					var timePeriod = dateObject.getHours() < 12 ? 'am' : 'pm';
					s = options[timePeriod] || bundle['dayPeriods-format-wide-' + timePeriod];
					break;
				case 'h':
				case 'H':
				case 'K':
				case 'k':
					var h = dateObject.getHours();
					// strange choices in the date format make it impossible to write this succinctly
					switch (c){
						case 'h': // 1-12
							s = (h % 12) || 12;
							break;
						case 'H': // 0-23
							s = h;
							break;
						case 'K': // 0-11
							s = (h % 12);
							break;
						case 'k': // 1-24
							s = h || 24;
							break;
					}
					pad = true;
					break;
				case 'm':
					s = dateObject.getMinutes(); pad = true;
					break;
				case 's':
					s = dateObject.getSeconds(); pad = true;
					break;
				case 'S':
					s = Math.round(dateObject.getMilliseconds() * Math.pow(10, l-3)); pad = true;
					break;
				case 'v': // FIXME: don't know what this is. seems to be same as z?
				case 'z':
					// We only have one timezone to offer; the one from the browser
					s = exports._getZone(dateObject, true, options);
					if(s){break;}
					l=4;
					// fallthrough... use GMT if tz not available
				case 'Z':
					var offset = exports._getZone(dateObject, false, options);
					var tz = [
						(offset<=0 ? "+" : "-"),
						string.pad(Math.floor(Math.abs(offset)/60), 2),
						string.pad(Math.abs(offset)% 60, 2)
					];
					if(l==4){
						tz.splice(0, 0, "GMT");
						tz.splice(3, 0, ":");
					}
					s = tz.join("");
					break;
//				case 'Y': case 'u': case 'W': case 'F': case 'g': case 'A':
//					console.log(match+" modifier unimplemented");
				default:
					throw new Error("dojo.date.locale.format: invalid pattern char: "+pattern);
			}
			if(pad){ s = string.pad(s, l); }
			return s;
		});
	}

/*=====
var __FormatOptions = exports.__FormatOptions = declare(null, {
	// selector: String
	//		choice of 'time','date' (default: date and time)
	// formatLength: String
	//		choice of long, short, medium or full (plus any custom additions).  Defaults to 'short'
	// datePattern:String
	//		override pattern with this string
	// timePattern:String
	//		override pattern with this string
	// am: String
	//		override strings for am in times
	// pm: String
	//		override strings for pm in times
	// locale: String
	//		override the locale used to determine formatting rules
	// fullYear: Boolean
	//		(format only) use 4 digit years whenever 2 digit years are called for
	// strict: Boolean
	//		(parse only) strict parsing, off by default
});
=====*/

exports._getZone = function(/*Date*/ dateObject, /*boolean*/ getName, /*__FormatOptions?*/ options){
	// summary:
	//		Returns the zone (or offset) for the given date and options.  This
	//		is broken out into a separate function so that it can be overridden
	//		by timezone-aware code.
	//
	// dateObject:
	//		the date and/or time being formatted.
	//
	// getName:
	//		Whether to return the timezone string (if true), or the offset (if false)
	//
	// options:
	//		The options being used for formatting
	if(getName){
		return date.getTimezoneName(dateObject);
	}else{
		return dateObject.getTimezoneOffset();
	}
};


exports.format = function(/*Date*/ dateObject, /*__FormatOptions?*/ options){
	// summary:
	//		Format a Date object as a String, using locale-specific settings.
	//
	// description:
	//		Create a string from a Date object using a known localized pattern.
	//		By default, this method formats both date and time from dateObject.
	//		Formatting patterns are chosen appropriate to the locale.  Different
	//		formatting lengths may be chosen, with "full" used by default.
	//		Custom patterns may be used or registered with translations using
	//		the dojo/date/locale.addCustomFormats() method.
	//		Formatting patterns are implemented using [the syntax described at
	//		unicode.org](http://www.unicode.org/reports/tr35/tr35-4.html#Date_Format_Patterns)
	//
	// dateObject:
	//		the date and/or time to be formatted.  If a time only is formatted,
	//		the values in the year, month, and day fields are irrelevant.  The
	//		opposite is true when formatting only dates.

	options = options || {};

	var locale = i18n.normalizeLocale(options.locale),
		formatLength = options.formatLength || 'short',
		bundle = exports._getGregorianBundle(locale),
		str = [],
		sauce = lang.hitch(this, formatPattern, dateObject, bundle, options);
	if(options.selector == "year"){
		return _processPattern(bundle["dateFormatItem-yyyy"] || "yyyy", sauce);
	}
	var pattern;
	if(options.selector != "date"){
		pattern = options.timePattern || bundle["timeFormat-"+formatLength];
		if(pattern){str.push(_processPattern(pattern, sauce));}
	}
	if(options.selector != "time"){
		pattern = options.datePattern || bundle["dateFormat-"+formatLength];
		if(pattern){str.push(_processPattern(pattern, sauce));}
	}

	return str.length == 1 ? str[0] : bundle["dateTimeFormat-"+formatLength].replace(/\'/g,'').replace(/\{(\d+)\}/g,
		function(match, key){ return str[key]; }); // String
};

exports.regexp = function(/*__FormatOptions?*/ options){
	// summary:
	//		Builds the regular needed to parse a localized date

	return exports._parseInfo(options).regexp; // String
};

exports._parseInfo = function(/*__FormatOptions?*/ options){
	options = options || {};
	var locale = i18n.normalizeLocale(options.locale),
		bundle = exports._getGregorianBundle(locale),
		formatLength = options.formatLength || 'short',
		datePattern = options.datePattern || bundle["dateFormat-" + formatLength],
		timePattern = options.timePattern || bundle["timeFormat-" + formatLength],
		pattern;
	if(options.selector == 'date'){
		pattern = datePattern;
	}else if(options.selector == 'time'){
		pattern = timePattern;
	}else{
		pattern = bundle["dateTimeFormat-"+formatLength].replace(/\{(\d+)\}/g,
			function(match, key){ return [timePattern, datePattern][key]; });
	}

	var tokens = [],
		re = _processPattern(pattern, lang.hitch(this, _buildDateTimeRE, tokens, bundle, options));
	return {regexp: re, tokens: tokens, bundle: bundle};
};

exports.parse = function(/*String*/ value, /*__FormatOptions?*/ options){
	// summary:
	//		Convert a properly formatted string to a primitive Date object,
	//		using locale-specific settings.
	//
	// description:
	//		Create a Date object from a string using a known localized pattern.
	//		By default, this method parses looking for both date and time in the string.
	//		Formatting patterns are chosen appropriate to the locale.  Different
	//		formatting lengths may be chosen, with "full" used by default.
	//		Custom patterns may be used or registered with translations using
	//		the dojo/date/locale.addCustomFormats() method.
	//
	//		Formatting patterns are implemented using [the syntax described at
	//		unicode.org](http://www.unicode.org/reports/tr35/tr35-4.html#Date_Format_Patterns)
	//		When two digit years are used, a century is chosen according to a sliding
	//		window of 80 years before and 20 years after present year, for both `yy` and `yyyy` patterns.
	//		year < 100CE requires strict mode.
	//
	// value:
	//		A string representation of a date

	// remove non-printing bidi control chars from input and pattern
	var controlChars = /[\u200E\u200F\u202A\u202E]/g,
		info = exports._parseInfo(options),
		tokens = info.tokens, bundle = info.bundle,
		re = new RegExp("^" + info.regexp.replace(controlChars, "") + "$",
			info.strict ? "" : "i"),
		match = re.exec(value && value.replace(controlChars, ""));

	if(!match){ return null; } // null

	var widthList = ['abbr', 'wide', 'narrow'],
		result = [1970,0,1,0,0,0,0], // will get converted to a Date at the end
		amPm = "",
		valid = array.every(match, function(v, i){
		if(!i){return true;}
		var token = tokens[i-1],
			l = token.length,
			c = token.charAt(0);
		switch(c){
			case 'y':
				if(l != 2 && options.strict){
					//interpret year literally, so '5' would be 5 A.D.
					result[0] = v;
				}else{
					if(v<100){
						v = Number(v);
						//choose century to apply, according to a sliding window
						//of 80 years before and 20 years after present year
						var year = '' + new Date().getFullYear(),
							century = year.substring(0, 2) * 100,
							cutoff = Math.min(Number(year.substring(2, 4)) + 20, 99);
						result[0] = (v < cutoff) ? century + v : century - 100 + v;
					}else{
						//we expected 2 digits and got more...
						if(options.strict){
							return false;
						}
						//interpret literally, so '150' would be 150 A.D.
						//also tolerate '1950', if 'yyyy' input passed to 'yy' format
						result[0] = v;
					}
				}
				break;
			case 'M':
			case 'L':
				if(l>2){
					var months = bundle['months-' +
							    (c == 'L' ? 'standAlone' : 'format') +
							    '-' + widthList[l-3]].concat();
					if(!options.strict){
						//Tolerate abbreviating period in month part
						//Case-insensitive comparison
						v = v.replace(".","").toLowerCase();
						months = array.map(months, function(s){ return s.replace(".","").toLowerCase(); } );
					}
					v = array.indexOf(months, v);
					if(v == -1){
//						console.log("dojo/date/locale.parse: Could not parse month name: '" + v + "'.");
						return false;
					}
				}else{
					v--;
				}
				result[1] = v;
				break;
			case 'E':
			case 'e':
			case 'c':
				var days = bundle['days-' +
						  (c == 'c' ? 'standAlone' : 'format') +
						  '-' + widthList[l-3]].concat();
				if(!options.strict){
					//Case-insensitive comparison
					v = v.toLowerCase();
					days = array.map(days, function(d){return d.toLowerCase();});
				}
				v = array.indexOf(days, v);
				if(v == -1){
//					console.log("dojo/date/locale.parse: Could not parse weekday name: '" + v + "'.");
					return false;
				}

				//TODO: not sure what to actually do with this input,
				//in terms of setting something on the Date obj...?
				//without more context, can't affect the actual date
				//TODO: just validate?
				break;
			case 'D':
				result[1] = 0;
				// fallthrough...
			case 'd':
				result[2] = v;
				break;
			case 'a': //am/pm
				var am = options.am || bundle['dayPeriods-format-wide-am'],
					pm = options.pm || bundle['dayPeriods-format-wide-pm'];
				if(!options.strict){
					var period = /\./g;
					v = v.replace(period,'').toLowerCase();
					am = am.replace(period,'').toLowerCase();
					pm = pm.replace(period,'').toLowerCase();
				}
				if(options.strict && v != am && v != pm){
//					console.log("dojo/date/locale.parse: Could not parse am/pm part.");
					return false;
				}

				// we might not have seen the hours field yet, so store the state and apply hour change later
				amPm = (v == pm) ? 'p' : (v == am) ? 'a' : '';
				break;
			case 'K': //hour (1-24)
				if(v == 24){ v = 0; }
				// fallthrough...
			case 'h': //hour (1-12)
			case 'H': //hour (0-23)
			case 'k': //hour (0-11)
				//TODO: strict bounds checking, padding
				if(v > 23){
//					console.log("dojo/date/locale.parse: Illegal hours value");
					return false;
				}

				//in the 12-hour case, adjusting for am/pm requires the 'a' part
				//which could come before or after the hour, so we will adjust later
				result[3] = v;
				break;
			case 'm': //minutes
				result[4] = v;
				break;
			case 's': //seconds
				result[5] = v;
				break;
			case 'S': //milliseconds
				result[6] = v;
//				break;
//			case 'w':
//TODO				var firstDay = 0;
//			default:
//TODO: throw?
//				console.log("dojo/date/locale.parse: unsupported pattern char=" + token.charAt(0));
		}
		return true;
	});

	var hours = +result[3];
	if(amPm === 'p' && hours < 12){
		result[3] = hours + 12; //e.g., 3pm -> 15
	}else if(amPm === 'a' && hours == 12){
		result[3] = 0; //12am -> 0
	}

	//TODO: implement a getWeekday() method in order to test
	//validity of input strings containing 'EEE' or 'EEEE'...

	var dateObject = new Date(result[0], result[1], result[2], result[3], result[4], result[5], result[6]); // Date
	if(options.strict){
		dateObject.setFullYear(result[0]);
	}

	// Check for overflow.  The Date() constructor normalizes things like April 32nd...
	//TODO: why isn't this done for times as well?
	var allTokens = tokens.join(""),
		dateToken = allTokens.indexOf('d') != -1,
		monthToken = allTokens.indexOf('M') != -1;

	if(!valid ||
		(monthToken && dateObject.getMonth() > result[1]) ||
		(dateToken && dateObject.getDate() > result[2])){
		return null;
	}

	// Check for underflow, due to DST shifts.  See #9366
	// This assumes a 1 hour dst shift correction at midnight
	// We could compare the timezone offset after the shift and add the difference instead.
	if((monthToken && dateObject.getMonth() < result[1]) ||
		(dateToken && dateObject.getDate() < result[2])){
		dateObject = date.add(dateObject, "hour", 1);
	}

	return dateObject; // Date
};

function _processPattern(pattern, applyPattern, applyLiteral, applyAll){
	//summary: Process a pattern with literals in it

	// Break up on single quotes, treat every other one as a literal, except '' which becomes '
	var identity = function(x){return x;};
	applyPattern = applyPattern || identity;
	applyLiteral = applyLiteral || identity;
	applyAll = applyAll || identity;

	//split on single quotes (which escape literals in date format strings)
	//but preserve escaped single quotes (e.g., o''clock)
	var chunks = pattern.match(/(''|[^'])+/g),
		literal = pattern.charAt(0) == "'";

	array.forEach(chunks, function(chunk, i){
		if(!chunk){
			chunks[i]='';
		}else{
			chunks[i]=(literal ? applyLiteral : applyPattern)(chunk.replace(/''/g, "'"));
			literal = !literal;
		}
	});
	return applyAll(chunks.join(''));
}

function _buildDateTimeRE(tokens, bundle, options, pattern){
	pattern = regexp.escapeString(pattern);
	if(!options.strict){ pattern = pattern.replace(" a", " ?a"); } // kludge to tolerate no space before am/pm
	return pattern.replace(/([a-z])\1*/ig, function(match){
		// Build a simple regexp.  Avoid captures, which would ruin the tokens list
		var s,
			c = match.charAt(0),
			l = match.length,
			p2 = '', p3 = '';
		if(options.strict){
			if(l > 1){ p2 = '0' + '{'+(l-1)+'}'; }
			if(l > 2){ p3 = '0' + '{'+(l-2)+'}'; }
		}else{
			p2 = '0?'; p3 = '0{0,2}';
		}
		switch(c){
			case 'y':
				s = '\\d{2,4}';
				break;
			case 'M':
			case 'L':
				s = (l>2) ? '\\S+?' : '1[0-2]|'+p2+'[1-9]';
				break;
			case 'D':
				s = '[12][0-9][0-9]|3[0-5][0-9]|36[0-6]|'+p2+'[1-9][0-9]|'+p3+'[1-9]';
				break;
			case 'd':
				s = '3[01]|[12]\\d|'+p2+'[1-9]';
				break;
			case 'w':
				s = '[1-4][0-9]|5[0-3]|'+p2+'[1-9]';
				break;
			case 'E':
			case 'e':
			case 'c':
				s = '\\S+';
				break;
			case 'h': //hour (1-12)
				s = '1[0-2]|'+p2+'[1-9]';
				break;
			case 'k': //hour (0-11)
				s = '1[01]|'+p2+'\\d';
				break;
			case 'H': //hour (0-23)
				s = '1\\d|2[0-3]|'+p2+'\\d';
				break;
			case 'K': //hour (1-24)
				s = '1\\d|2[0-4]|'+p2+'[1-9]';
				break;
			case 'm':
			case 's':
				s = '[0-5]\\d';
				break;
			case 'S':
				s = '\\d{'+l+'}';
				break;
			case 'a':
				var am = options.am || bundle['dayPeriods-format-wide-am'],
					pm = options.pm || bundle['dayPeriods-format-wide-pm'];
					s = am + '|' + pm;
				if(!options.strict){
					if(am != am.toLowerCase()){ s += '|' + am.toLowerCase(); }
					if(pm != pm.toLowerCase()){ s += '|' + pm.toLowerCase(); }
					if(s.indexOf('.') != -1){ s += '|' + s.replace(/\./g, ""); }
				}
				s = s.replace(/\./g, "\\.");
				break;
			default:
			// case 'v':
			// case 'z':
			// case 'Z':
				s = ".*";
//				console.log("parse of date format, pattern=" + pattern);
		}

		if(tokens){ tokens.push(match); }

		return "(" + s + ")"; // add capture
	}).replace(/[\xa0 ]/g, "[\\s\\xa0]"); // normalize whitespace.  Need explicit handling of \xa0 for IE.
}

var _customFormats = [];
exports.addCustomFormats = function(/*String*/ packageName, /*String*/ bundleName){
	// summary:
	//		Add a reference to a bundle containing localized custom formats to be
	//		used by date/time formatting and parsing routines.
	//
	// description:
	//		The user may add custom localized formats where the bundle has properties following the
	//		same naming convention used by dojo.cldr: `dateFormat-xxxx` / `timeFormat-xxxx`
	//		The pattern string should match the format used by the CLDR.
	//		See dojo/date/locale.format() for details.
	//		The resources must be loaded by dojo.requireLocalization() prior to use

	_customFormats.push({pkg:packageName,name:bundleName});
};

exports._getGregorianBundle = function(/*String*/ locale){
	var gregorian = {};
	array.forEach(_customFormats, function(desc){
		var bundle = i18n.getLocalization(desc.pkg, desc.name, locale);
		gregorian = lang.mixin(gregorian, bundle);
	}, this);
	return gregorian; /*Object*/
};

exports.addCustomFormats(module.id.replace(/\/date\/locale$/, ".cldr"),"gregorian");

exports.getNames = function(/*String*/ item, /*String*/ type, /*String?*/ context, /*String?*/ locale){
	// summary:
	//		Used to get localized strings from dojo.cldr for day or month names.
	//
	// item:
	//	'months' || 'days'
	// type:
	//	'wide' || 'abbr' || 'narrow' (e.g. "Monday", "Mon", or "M" respectively, in English)
	// context:
	//	'standAlone' || 'format' (default)
	// locale:
	//	override locale used to find the names

	var label,
		lookup = exports._getGregorianBundle(locale),
		props = [item, context, type];
	if(context == 'standAlone'){
		var key = props.join('-');
		label = lookup[key];
		// Fall back to 'format' flavor of name
		if(label[0] == 1){ label = undefined; } // kludge, in the absence of real aliasing support in dojo.cldr
	}
	props[1] = 'format';

	// return by copy so changes won't be made accidentally to the in-memory model
	return (label || lookup[props.join('-')]).concat(); /*Array*/
};

exports.isWeekend = function(/*Date?*/ dateObject, /*String?*/ locale){
	// summary:
	//	Determines if the date falls on a weekend, according to local custom.

	var weekend = supplemental.getWeekend(locale),
		day = (dateObject || new Date()).getDay();
	if(weekend.end < weekend.start){
		weekend.end += 7;
		if(day < weekend.start){ day += 7; }
	}
	return day >= weekend.start && day <= weekend.end; // Boolean
};

// These are used only by format and strftime.  Do they need to be public?  Which module should they go in?

exports._getDayOfYear = function(/*Date*/ dateObject){
	// summary:
	//		gets the day of the year as represented by dateObject
	return date.difference(new Date(dateObject.getFullYear(), 0, 1, dateObject.getHours()), dateObject) + 1; // Number
};

exports._getWeekOfYear = function(/*Date*/ dateObject, /*Number*/ firstDayOfWeek){
	if(arguments.length == 1){ firstDayOfWeek = 0; } // Sunday

	var firstDayOfYear = new Date(dateObject.getFullYear(), 0, 1).getDay(),
		adj = (firstDayOfYear - firstDayOfWeek + 7) % 7,
		week = Math.floor((exports._getDayOfYear(dateObject) + adj - 1) / 7);

	// if year starts on the specified day, start counting weeks at 1
	if(firstDayOfYear == firstDayOfWeek){ week++; }

	return week; // Number
};

return exports;
});

},
'dojo/cldr/supplemental':function(){
define("dojo/cldr/supplemental", ["../_base/lang", "../i18n"], function(lang, i18n){

// module:
//		dojo/cldr/supplemental


var supplemental = {
	// summary:
	//		TODOC
};
lang.setObject("dojo.cldr.supplemental", supplemental);

supplemental.getFirstDayOfWeek = function(/*String?*/locale){
	// summary:
	//		Returns a zero-based index for first day of the week
	// description:
	//		Returns a zero-based index for first day of the week, as used by the local (Gregorian) calendar.
	//		e.g. Sunday (returns 0), or Monday (returns 1)

	// from http://www.unicode.org/cldr/data/common/supplemental/supplementalData.xml:supplementalData/weekData/firstDay
	var firstDay = {/*default is 1=Monday*/
		bd:5,mv:5,
		ae:6,af:6,bh:6,dj:6,dz:6,eg:6,iq:6,ir:6,jo:6,kw:6,
		ly:6,ma:6,om:6,qa:6,sa:6,sd:6,sy:6,ye:6,
		ag:0,ar:0,as:0,au:0,br:0,bs:0,bt:0,bw:0,by:0,bz:0,ca:0,cn:0,
		co:0,dm:0,'do':0,et:0,gt:0,gu:0,hk:0,hn:0,id:0,ie:0,il:0,'in':0,
		jm:0,jp:0,ke:0,kh:0,kr:0,la:0,mh:0,mm:0,mo:0,mt:0,mx:0,mz:0,
		ni:0,np:0,nz:0,pa:0,pe:0,ph:0,pk:0,pr:0,py:0,sg:0,sv:0,th:0,
		tn:0,tt:0,tw:0,um:0,us:0,ve:0,vi:0,ws:0,za:0,zw:0
	};

	var country = supplemental._region(locale);
	var dow = firstDay[country];
	return (dow === undefined) ? 1 : dow; /*Number*/
};

supplemental._region = function(/*String?*/locale){
	locale = i18n.normalizeLocale(locale);
	var tags = locale.split('-');
	var region = tags[1];
	if(!region){
		// IE often gives language only (#2269)
		// Arbitrary mappings of language-only locales to a country:
		region = {de:"de", en:"us", es:"es", fi:"fi", fr:"fr", he:"il", hu:"hu", it:"it",
			ja:"jp", ko:"kr", nl:"nl", pt:"br", sv:"se", zh:"cn"}[tags[0]];
	}else if(region.length == 4){
		// The ISO 3166 country code is usually in the second position, unless a
		// 4-letter script is given. See http://www.ietf.org/rfc/rfc4646.txt
		region = tags[2];
	}
	return region;
};

supplemental.getWeekend = function(/*String?*/locale){
	// summary:
	//		Returns a hash containing the start and end days of the weekend
	// description:
	//		Returns a hash containing the start and end days of the weekend according to local custom using locale,
	//		or by default in the user's locale.
	//		e.g. {start:6, end:0}

	// from http://www.unicode.org/cldr/data/common/supplemental/supplementalData.xml:supplementalData/weekData/weekend{Start,End}
	var weekendStart = {/*default is 6=Saturday*/
			'in':0,
			af:4,dz:4,ir:4,om:4,sa:4,ye:4,
			ae:5,bh:5,eg:5,il:5,iq:5,jo:5,kw:5,ly:5,ma:5,qa:5,sd:5,sy:5,tn:5
		},

		weekendEnd = {/*default is 0=Sunday*/
			af:5,dz:5,ir:5,om:5,sa:5,ye:5,
			ae:6,bh:5,eg:6,il:6,iq:6,jo:6,kw:6,ly:6,ma:6,qa:6,sd:6,sy:6,tn:6
		},

		country = supplemental._region(locale),
		start = weekendStart[country],
		end = weekendEnd[country];

	if(start === undefined){start=6;}
	if(end === undefined){end=0;}
	return {start:start, end:end}; /*Object {start,end}*/
};

return supplemental;
});

},
'dijit/CalendarLite':function(){
require({cache:{
'url:dijit/templates/Calendar.html':"<table cellspacing=\"0\" cellpadding=\"0\" class=\"dijitCalendarContainer\" role=\"grid\" aria-labelledby=\"${id}_mddb ${id}_year\">\n\t<thead>\n\t\t<tr class=\"dijitReset dijitCalendarMonthContainer\" valign=\"top\">\n\t\t\t<th class='dijitReset dijitCalendarArrow' data-dojo-attach-point=\"decrementMonth\">\n\t\t\t\t<img src=\"${_blankGif}\" alt=\"\" class=\"dijitCalendarIncrementControl dijitCalendarDecrease\" role=\"presentation\"/>\n\t\t\t\t<span data-dojo-attach-point=\"decreaseArrowNode\" class=\"dijitA11ySideArrow\">-</span>\n\t\t\t</th>\n\t\t\t<th class='dijitReset' colspan=\"5\">\n\t\t\t\t<div data-dojo-attach-point=\"monthNode\">\n\t\t\t\t</div>\n\t\t\t</th>\n\t\t\t<th class='dijitReset dijitCalendarArrow' data-dojo-attach-point=\"incrementMonth\">\n\t\t\t\t<img src=\"${_blankGif}\" alt=\"\" class=\"dijitCalendarIncrementControl dijitCalendarIncrease\" role=\"presentation\"/>\n\t\t\t\t<span data-dojo-attach-point=\"increaseArrowNode\" class=\"dijitA11ySideArrow\">+</span>\n\t\t\t</th>\n\t\t</tr>\n\t\t<tr role=\"row\">\n\t\t\t${!dayCellsHtml}\n\t\t</tr>\n\t</thead>\n\t<tbody data-dojo-attach-point=\"dateRowsNode\" data-dojo-attach-event=\"onclick: _onDayClick\" class=\"dijitReset dijitCalendarBodyContainer\">\n\t\t\t${!dateRowsHtml}\n\t</tbody>\n\t<tfoot class=\"dijitReset dijitCalendarYearContainer\">\n\t\t<tr>\n\t\t\t<td class='dijitReset' valign=\"top\" colspan=\"7\" role=\"presentation\">\n\t\t\t\t<div class=\"dijitCalendarYearLabel\">\n\t\t\t\t\t<span data-dojo-attach-point=\"previousYearLabelNode\" class=\"dijitInline dijitCalendarPreviousYear\" role=\"button\"></span>\n\t\t\t\t\t<span data-dojo-attach-point=\"currentYearLabelNode\" class=\"dijitInline dijitCalendarSelectedYear\" role=\"button\" id=\"${id}_year\"></span>\n\t\t\t\t\t<span data-dojo-attach-point=\"nextYearLabelNode\" class=\"dijitInline dijitCalendarNextYear\" role=\"button\"></span>\n\t\t\t\t</div>\n\t\t\t</td>\n\t\t</tr>\n\t</tfoot>\n</table>\n"}});
define("dijit/CalendarLite", [
	"dojo/_base/array", // array.forEach array.map
	"dojo/_base/declare", // declare
	"dojo/cldr/supplemental", // cldrSupplemental.getFirstDayOfWeek
	"dojo/date", // date
	"dojo/date/locale",
	"dojo/date/stamp", // stamp.fromISOString
	"dojo/dom", // dom.setSelectable
	"dojo/dom-class", // domClass.contains
	"dojo/_base/event", // event.stop
	"dojo/_base/lang", // lang.getObject, lang.hitch
	"dojo/sniff", // has("ie") has("webkit")
	"dojo/string", // string.substitute
	"./_WidgetBase",
	"./_TemplatedMixin",
	"dojo/text!./templates/Calendar.html",
	"./hccss"	// not used directly, but sets CSS class on <body>
], function(array, declare, cldrSupplemental, date, locale, stamp, dom, domClass, event, lang, has, string,
			_WidgetBase, _TemplatedMixin, template){


	// module:
	//		dijit/CalendarLite

	var CalendarLite = declare("dijit.CalendarLite", [_WidgetBase, _TemplatedMixin], {
		// summary:
		//		Lightweight version of Calendar widget aimed towards mobile use
		//
		// description:
		//		A simple GUI for choosing a date in the context of a monthly calendar.
		//		This widget can't be used in a form because it doesn't serialize the date to an
		//		`<input>` field.  For a form element, use dijit/form/DateTextBox instead.
		//
		//		Note that the parser takes all dates attributes passed in the
		//		[RFC 3339 format](http://www.faqs.org/rfcs/rfc3339.html), e.g. `2005-06-30T08:05:00-07:00`
		//		so that they are serializable and locale-independent.
		//
		//		Also note that this widget isn't keyboard accessible; use dijit.Calendar for that
		// example:
		//	|	var calendar = new dijit.CalendarLite({}, dojo.byId("calendarNode"));
		//
		// example:
		//	|	<div data-dojo-type="dijit/CalendarLite"></div>

		// Template for main calendar
		templateString: template,

		// Template for cell for a day of the week (ex: M)
		dowTemplateString: '<th class="dijitReset dijitCalendarDayLabelTemplate" role="columnheader"><span class="dijitCalendarDayLabel">${d}</span></th>',

		// Templates for a single date (ex: 13), and for a row for a week (ex: 20 21 22 23 24 25 26)
		dateTemplateString: '<td class="dijitReset" role="gridcell" data-dojo-attach-point="dateCells"><span class="dijitCalendarDateLabel" data-dojo-attach-point="dateLabels"></span></td>',
		weekTemplateString: '<tr class="dijitReset dijitCalendarWeekTemplate" role="row">${d}${d}${d}${d}${d}${d}${d}</tr>',

		// value: Date
		//		The currently selected Date, initially set to invalid date to indicate no selection.
		value: new Date(""),
		// TODO: for 2.0 make this a string (ISO format) rather than a Date

		// datePackage: String
		//		JavaScript namespace to find calendar routines.	 If unspecified, uses Gregorian calendar routines
		//		at dojo/date and dojo/date/locale.
		datePackage: "",
		//		TODO: for 2.0, replace datePackage with dateModule and dateLocalModule attributes specifying MIDs,
		//		or alternately just get rid of this completely and tell user to use module ID remapping
		//		via require

		// dayWidth: String
		//		How to represent the days of the week in the calendar header. See locale
		dayWidth: "narrow",

		// tabIndex: String
		//		Order fields are traversed when user hits the tab key
		tabIndex: "0",

		// currentFocus: Date
		//		Date object containing the currently focused date, or the date which would be focused
		//		if the calendar itself was focused.   Also indicates which year and month to display,
		//		i.e. the current "page" the calendar is on.
		currentFocus: new Date(),

		baseClass:"dijitCalendar",

		_isValidDate: function(/*Date*/ value){
			// summary:
			//		Runs various tests on the value, checking that it's a valid date, rather
			//		than blank or NaN.
			// tags:
			//		private
			return value && !isNaN(value) && typeof value == "object" &&
				value.toString() != this.constructor.prototype.value.toString();
		},

		_getValueAttr: function(){
			// summary:
			//		Support get('value')

			// this.value is set to 1AM, but return midnight, local time for back-compat
			if(this.value && !isNaN(this.value)){
				var value = new this.dateClassObj(this.value);
				value.setHours(0, 0, 0, 0);

				// If daylight savings pushes midnight to the previous date, fix the Date
				// object to point at 1am so it will represent the correct day. See #9366
				if(value.getDate() < this.value.getDate()){
					value = this.dateModule.add(value, "hour", 1);
				}
				return value;
			}else{
				return null;
			}
		},

		_setValueAttr: function(/*Date|Number*/ value, /*Boolean*/ priorityChange){
			// summary:
			//		Support set("value", ...)
			// description:
			//		Set the current date and update the UI.  If the date is disabled, the value will
			//		not change, but the display will change to the corresponding month.
			// value:
			//		Either a Date or the number of seconds since 1970.
			// tags:
			//		protected
			if(typeof value == "string"){
				value = stamp.fromISOString(value);
			}
			value = this._patchDate(value);
			
			if(this._isValidDate(value) && !this.isDisabledDate(value, this.lang)){
				this._set("value", value);

				// Set focus cell to the new value.   Arguably this should only happen when there isn't a current
				// focus point.   This will also repopulate the grid to new month/year if necessary.
				this.set("currentFocus", value);

				// Mark the selected date
				this._markSelectedDates([value]);

				if(this._created && (priorityChange || typeof priorityChange == "undefined")){
					this.onChange(this.get('value'));
				}
			}else{
				// clear value, and mark all dates as unselected
				this._set("value", null);
				this._markSelectedDates([]);
			}
		},

		_patchDate: function(/*Date|Number*/ value){
			// summary:
			//		Convert Number into Date, or copy Date object.   Then, round to nearest day,
			//		setting to 1am to avoid issues when DST shift occurs at midnight, see #8521, #9366)
			if(value){
				value = new this.dateClassObj(value);
				value.setHours(1, 0, 0, 0);
			}
			return value;
		},

		_setText: function(node, text){
			// summary:
			//		This just sets the content of node to the specified text.
			//		Can't do "node.innerHTML=text" because of an IE bug w/tables, see #3434.
			// tags:
			//		private
			while(node.firstChild){
				node.removeChild(node.firstChild);
			}
			node.appendChild(node.ownerDocument.createTextNode(text));
		},

		_populateGrid: function(){
			// summary:
			//		Fills in the calendar grid with each day (1-31).
			//		Call this on creation, when moving to a new month.
			// tags:
			//		private

			var month = new this.dateClassObj(this.currentFocus);
			month.setDate(1);

			var firstDay = month.getDay(),
				daysInMonth = this.dateModule.getDaysInMonth(month),
				daysInPreviousMonth = this.dateModule.getDaysInMonth(this.dateModule.add(month, "month", -1)),
				today = new this.dateClassObj(),
				dayOffset = cldrSupplemental.getFirstDayOfWeek(this.lang);
			if(dayOffset > firstDay){ dayOffset -= 7; }

			// Mapping from date (as specified by number returned from Date.valueOf()) to corresponding <td>
			this._date2cell = {};

			// Iterate through dates in the calendar and fill in date numbers and style info
			array.forEach(this.dateCells, function(template, idx){
				var i = idx + dayOffset;
				var date = new this.dateClassObj(month),
					number, clazz = "dijitCalendar", adj = 0;

				if(i < firstDay){
					number = daysInPreviousMonth - firstDay + i + 1;
					adj = -1;
					clazz += "Previous";
				}else if(i >= (firstDay + daysInMonth)){
					number = i - firstDay - daysInMonth + 1;
					adj = 1;
					clazz += "Next";
				}else{
					number = i - firstDay + 1;
					clazz += "Current";
				}

				if(adj){
					date = this.dateModule.add(date, "month", adj);
				}
				date.setDate(number);

				if(!this.dateModule.compare(date, today, "date")){
					clazz = "dijitCalendarCurrentDate " + clazz;
				}

				if(this.isDisabledDate(date, this.lang)){
					clazz = "dijitCalendarDisabledDate " + clazz;
					template.setAttribute("aria-disabled", "true");
				}else{
					clazz = "dijitCalendarEnabledDate " + clazz;
					template.removeAttribute("aria-disabled");
					template.setAttribute("aria-selected", "false");
				}

				var clazz2 = this.getClassForDate(date, this.lang);
				if(clazz2){
					clazz = clazz2 + " " + clazz;
				}

				template.className = clazz + "Month dijitCalendarDateTemplate";

				// Each cell has an associated integer value representing it's date
				var dateVal = date.valueOf();
				this._date2cell[dateVal] = template;
				template.dijitDateValue = dateVal;

				// Set Date string (ex: "13").
				this._setText(this.dateLabels[idx], date.getDateLocalized ? date.getDateLocalized(this.lang) : date.getDate());
			}, this);
		},
		
		_populateControls: function(){
			// summary:
			//		Fill in localized month, and prev/current/next years
			// tags:
			//		protected

			var month = new this.dateClassObj(this.currentFocus);
			month.setDate(1);
			
			// set name of this month
			this.monthWidget.set("month", month);
			
			var y = month.getFullYear() - 1;
			var d = new this.dateClassObj();
			array.forEach(["previous", "current", "next"], function(name){
				d.setFullYear(y++);
				this._setText(this[name+"YearLabelNode"],
					this.dateLocaleModule.format(d, {selector:'year', locale:this.lang}));
			}, this);
		},

		goToToday: function(){
			// summary:
			//		Sets calendar's value to today's date
			this.set('value', new this.dateClassObj());
		},

		constructor: function(params /*===== , srcNodeRef =====*/){
			// summary:
			//		Create the widget.
			// params: Object|null
			//		Hash of initialization parameters for widget, including scalar values (like title, duration etc.)
			//		and functions, typically callbacks like onClick.
			//		The hash can contain any of the widget's properties, excluding read-only properties.
			// srcNodeRef: DOMNode|String?
			//		If a srcNodeRef (DOM node) is specified, replace srcNodeRef with my generated DOM tree

			this.dateModule = params.datePackage ? lang.getObject(params.datePackage, false) : date;
			this.dateClassObj = this.dateModule.Date || Date;
			this.dateLocaleModule = params.datePackage ? lang.getObject(params.datePackage+".locale", false) : locale;
		},

		_createMonthWidget: function(){
			// summary:
			//		Creates the drop down button that displays the current month and lets user pick a new one

			return CalendarLite._MonthWidget({
				id: this.id + "_mw",
				lang: this.lang,
				dateLocaleModule: this.dateLocaleModule
			}, this.monthNode);
		},

		buildRendering: function(){
			// Markup for days of the week (referenced from template)
			var d = this.dowTemplateString,
				dayNames = this.dateLocaleModule.getNames('days', this.dayWidth, 'standAlone', this.lang),
				dayOffset = cldrSupplemental.getFirstDayOfWeek(this.lang);
			this.dayCellsHtml = string.substitute([d,d,d,d,d,d,d].join(""), {d: ""}, function(){
				return dayNames[dayOffset++ % 7];
			});

			// Markup for dates of the month (referenced from template), but without numbers filled in
			var r = string.substitute(this.weekTemplateString, {d: this.dateTemplateString});
			this.dateRowsHtml = [r,r,r,r,r,r].join("");

			// Instantiate from template.
			// dateCells and dateLabels arrays filled when _Templated parses my template.
			this.dateCells = [];
			this.dateLabels = [];
			this.inherited(arguments);

			dom.setSelectable(this.domNode, false);

			var dateObj = new this.dateClassObj(this.currentFocus);

			this.monthWidget = this._createMonthWidget();

			this.set('currentFocus', dateObj, false);	// draw the grid to the month specified by currentFocus
		},

		postCreate: function(){
			this.inherited(arguments);
			this._connectControls();
		},

		_connectControls: function(){
			// summary:
			//		Set up connects for increment/decrement of months/years
			// tags:
			//		protected

			var connect = lang.hitch(this, function(nodeProp, part, amount){
				this.connect(this[nodeProp], "onclick", function(){
					this._setCurrentFocusAttr(this.dateModule.add(this.currentFocus, part, amount));
				});
			});
			
			connect("incrementMonth", "month", 1);
			connect("decrementMonth", "month", -1);
			connect("nextYearLabelNode", "year", 1);
			connect("previousYearLabelNode", "year", -1);
		},
		
		_setCurrentFocusAttr: function(/*Date*/ date, /*Boolean*/ forceFocus){
			// summary:
			//		If the calendar currently has focus, then focuses specified date,
			//		changing the currently displayed month/year if necessary.
			//		If the calendar doesn't have focus, updates currently
			//		displayed month/year, and sets the cell that will get focus
			//		when Calendar is focused.
			// forceFocus:
			//		If true, will focus() the cell even if calendar itself doesn't have focus

			var oldFocus = this.currentFocus,
				oldCell = this._getNodeByDate(oldFocus);
			date = this._patchDate(date);

			this._set("currentFocus", date);

			// If the focus is on a different month than the current calendar month, switch the displayed month.
			// Also will populate the grid initially, on Calendar creation.
			if(!this._date2cell || this.dateModule.difference(oldFocus, date, "month") != 0){
				this._populateGrid();
				this._populateControls();
				this._markSelectedDates([this.value]);
			}
			
			// set tabIndex=0 on new cell, and focus it (but only if Calendar itself is focused)
			var newCell = this._getNodeByDate(date);
			newCell.setAttribute("tabIndex", this.tabIndex);
			if(this.focused || forceFocus){
				newCell.focus();
			}

			// set tabIndex=-1 on old focusable cell
			if(oldCell && oldCell != newCell){
				if(has("webkit")){	// see #11064 about webkit bug
					oldCell.setAttribute("tabIndex", "-1");
				}else{
					oldCell.removeAttribute("tabIndex");
				}
			}
		},

		focus: function(){
			// summary:
			//		Focus the calendar by focusing one of the calendar cells
			this._setCurrentFocusAttr(this.currentFocus, true);
		},

		_onDayClick: function(/*Event*/ evt){
			// summary:
			//		Handler for day clicks, selects the date if appropriate
			// tags:
			//		protected
			event.stop(evt);
			for(var node = evt.target; node && !node.dijitDateValue; node = node.parentNode);
			if(node && !domClass.contains(node, "dijitCalendarDisabledDate")){
				this.set('value', node.dijitDateValue);
			}
		},

		_getNodeByDate : function(/*Date*/ value){
			// summary:
			//		Returns the cell corresponding to the date, or null if the date is not within the currently
			//		displayed month.
			value = this._patchDate(value);
			return value && this._date2cell ? this._date2cell[value.valueOf()] : null;
		},

		_markSelectedDates: function(/*Date[]*/ dates){
			// summary:
			//		Marks the specified cells as selected, and clears cells previously marked as selected.
			//		For CalendarLite at most one cell is selected at any point, but this allows an array
			//		for easy subclassing.

			// Function to mark a cell as selected or unselected
			function mark(/*Boolean*/ selected, /*DomNode*/ cell){
				domClass.toggle(cell, "dijitCalendarSelectedDate", selected);
				cell.setAttribute("aria-selected", selected ? "true" : "false");
			}

			// Clear previously selected cells.
			array.forEach(this._selectedCells || [], lang.partial(mark, false));

			// Mark newly selected cells.  Ignore dates outside the currently displayed month.
			this._selectedCells = array.filter(array.map(dates, this._getNodeByDate, this), function(n){ return n;});
			array.forEach(this._selectedCells, lang.partial(mark, true));
		},

		onChange: function(/*Date*/ /*===== date =====*/){
			// summary:
			//		Called only when the selected date has changed
		},

		isDisabledDate: function(/*===== dateObject, locale =====*/){
			// summary:
			//		May be overridden to disable certain dates in the calendar e.g. `isDisabledDate=dojo.date.locale.isWeekend`
			// dateObject: Date
			// locale: String?
			// tags:
			//		extension
/*=====
			return false; // Boolean
=====*/
		},

		getClassForDate: function(/*===== dateObject, locale =====*/){
			// summary:
			//		May be overridden to return CSS classes to associate with the date entry for the given dateObject,
			//		for example to indicate a holiday in specified locale.
			// dateObject: Date
			// locale: String?
			// tags:
			//		extension

/*=====
			return ""; // String
=====*/
		}
	});

	CalendarLite._MonthWidget = declare("dijit.CalendarLite._MonthWidget", _WidgetBase, {
		// summary:
		//		Displays name of current month padded to the width of the month
		//		w/the longest name, so that changing months doesn't change width.
		//
		//		Create as:
		// |	new Calendar._MonthWidget({
		// |			lang: ...,
		// |			dateLocaleModule: ...
		// |		})

		_setMonthAttr: function(month){
			// summary:
			//		Set the current month to display as a label
			var monthNames = this.dateLocaleModule.getNames('months', 'wide', 'standAlone', this.lang, month),
				spacer =
					(has("ie") == 6 ? "" :	"<div class='dijitSpacer'>" +
						array.map(monthNames, function(s){ return "<div>" + s + "</div>"; }).join("") + "</div>");

			// Set name of current month and also fill in spacer element with all the month names
			// (invisible) so that the maximum width will affect layout.   But not on IE6 because then
			// the center <TH> overlaps the right <TH> (due to a browser bug).
			this.domNode.innerHTML =
				spacer +
				"<div class='dijitCalendarMonthLabel dijitCalendarCurrentMonthLabel'>" +
				monthNames[month.getMonth()] + "</div>";
		}
	});

	return CalendarLite;
});

},
'url:dijit/templates/Calendar.html':"<table cellspacing=\"0\" cellpadding=\"0\" class=\"dijitCalendarContainer\" role=\"grid\" aria-labelledby=\"${id}_mddb ${id}_year\">\n\t<thead>\n\t\t<tr class=\"dijitReset dijitCalendarMonthContainer\" valign=\"top\">\n\t\t\t<th class='dijitReset dijitCalendarArrow' data-dojo-attach-point=\"decrementMonth\">\n\t\t\t\t<img src=\"${_blankGif}\" alt=\"\" class=\"dijitCalendarIncrementControl dijitCalendarDecrease\" role=\"presentation\"/>\n\t\t\t\t<span data-dojo-attach-point=\"decreaseArrowNode\" class=\"dijitA11ySideArrow\">-</span>\n\t\t\t</th>\n\t\t\t<th class='dijitReset' colspan=\"5\">\n\t\t\t\t<div data-dojo-attach-point=\"monthNode\">\n\t\t\t\t</div>\n\t\t\t</th>\n\t\t\t<th class='dijitReset dijitCalendarArrow' data-dojo-attach-point=\"incrementMonth\">\n\t\t\t\t<img src=\"${_blankGif}\" alt=\"\" class=\"dijitCalendarIncrementControl dijitCalendarIncrease\" role=\"presentation\"/>\n\t\t\t\t<span data-dojo-attach-point=\"increaseArrowNode\" class=\"dijitA11ySideArrow\">+</span>\n\t\t\t</th>\n\t\t</tr>\n\t\t<tr role=\"row\">\n\t\t\t${!dayCellsHtml}\n\t\t</tr>\n\t</thead>\n\t<tbody data-dojo-attach-point=\"dateRowsNode\" data-dojo-attach-event=\"onclick: _onDayClick\" class=\"dijitReset dijitCalendarBodyContainer\">\n\t\t\t${!dateRowsHtml}\n\t</tbody>\n\t<tfoot class=\"dijitReset dijitCalendarYearContainer\">\n\t\t<tr>\n\t\t\t<td class='dijitReset' valign=\"top\" colspan=\"7\" role=\"presentation\">\n\t\t\t\t<div class=\"dijitCalendarYearLabel\">\n\t\t\t\t\t<span data-dojo-attach-point=\"previousYearLabelNode\" class=\"dijitInline dijitCalendarPreviousYear\" role=\"button\"></span>\n\t\t\t\t\t<span data-dojo-attach-point=\"currentYearLabelNode\" class=\"dijitInline dijitCalendarSelectedYear\" role=\"button\" id=\"${id}_year\"></span>\n\t\t\t\t\t<span data-dojo-attach-point=\"nextYearLabelNode\" class=\"dijitInline dijitCalendarNextYear\" role=\"button\"></span>\n\t\t\t\t</div>\n\t\t\t</td>\n\t\t</tr>\n\t</tfoot>\n</table>\n",
'dijit/form/_DateTimeTextBox':function(){
require({cache:{
'url:dijit/form/templates/DropDownBox.html':"<div class=\"dijit dijitReset dijitInline dijitLeft\"\n\tid=\"widget_${id}\"\n\trole=\"combobox\"\n\t><div class='dijitReset dijitRight dijitButtonNode dijitArrowButton dijitDownArrowButton dijitArrowButtonContainer'\n\t\tdata-dojo-attach-point=\"_buttonNode, _popupStateNode\" role=\"presentation\"\n\t\t><input class=\"dijitReset dijitInputField dijitArrowButtonInner\" value=\"&#9660; \" type=\"text\" tabIndex=\"-1\" readonly=\"readonly\" role=\"button presentation\" aria-hidden=\"true\"\n\t\t\t${_buttonInputDisabled}\n\t/></div\n\t><div class='dijitReset dijitValidationContainer'\n\t\t><input class=\"dijitReset dijitInputField dijitValidationIcon dijitValidationInner\" value=\"&#935; \" type=\"text\" tabIndex=\"-1\" readonly=\"readonly\" role=\"presentation\"\n\t/></div\n\t><div class=\"dijitReset dijitInputField dijitInputContainer\"\n\t\t><input class='dijitReset dijitInputInner' ${!nameAttrSetting} type=\"text\" autocomplete=\"off\"\n\t\t\tdata-dojo-attach-point=\"textbox,focusNode\" role=\"textbox\" aria-haspopup=\"true\"\n\t/></div\n></div>\n"}});
define("dijit/form/_DateTimeTextBox", [
	"dojo/date", // date date.compare
	"dojo/date/locale", // locale.regexp
	"dojo/date/stamp", // stamp.fromISOString stamp.toISOString
	"dojo/_base/declare", // declare
	"dojo/_base/lang", // lang.getObject
	"./RangeBoundTextBox",
	"../_HasDropDown",
	"dojo/text!./templates/DropDownBox.html"
], function(date, locale, stamp, declare, lang, RangeBoundTextBox, _HasDropDown, template){

	// module:
	//		dijit/form/_DateTimeTextBox

	new Date("X"); // workaround for #11279, new Date("") == NaN

	var _DateTimeTextBox = declare("dijit.form._DateTimeTextBox", [RangeBoundTextBox, _HasDropDown], {
		// summary:
		//		Base class for validating, serializable, range-bound date or time text box.

		templateString: template,

		// hasDownArrow: [const] Boolean
		//		Set this textbox to display a down arrow button, to open the drop down list.
		hasDownArrow: true,

		// Set classes like dijitDownArrowButtonHover depending on mouse action over button node
		cssStateNodes: {
			"_buttonNode": "dijitDownArrowButton"
		},

		/*=====
		// constraints: _DateTimeTextBox.__Constraints
		//		Despite the name, this parameter specifies both constraints on the input
		//		(including starting/ending dates/times allowed) as well as
		//		formatting options like whether the date is displayed in long (ex: December 25, 2005)
		//		or short (ex: 12/25/2005) format.  See `dijit/form/_DateTimeTextBox.__Constraints` for details.
		constraints: {},
		======*/

		// Override ValidationTextBox.pattern.... we use a reg-ex generating function rather
		// than a straight regexp to deal with locale  (plus formatting options too?)
		pattern: locale.regexp,

		// datePackage: String
		//		JavaScript namespace to find calendar routines.	 If unspecified, uses Gregorian calendar routines
		//		at dojo/date and dojo/date/locale.
		datePackage: "",
		//		TODO: for 2.0, replace datePackage with dateModule and dateLocalModule attributes specifying MIDs,
		//		or alternately just get rid of this completely and tell user to use module ID remapping
		//		via require

		postMixInProperties: function(){
			this.inherited(arguments);
			this._set("type", "text"); // in case type="date"|"time" was specified which messes up parse/format
		},

		// Override _FormWidget.compare() to work for dates/times
		compare: function(/*Date*/ val1, /*Date*/ val2){
			var isInvalid1 = this._isInvalidDate(val1);
			var isInvalid2 = this._isInvalidDate(val2);
			return isInvalid1 ? (isInvalid2 ? 0 : -1) : (isInvalid2 ? 1 : date.compare(val1, val2, this._selector));
		},

		// flag to _HasDropDown to make drop down Calendar width == <input> width
		forceWidth: true,

		format: function(/*Date*/ value, /*locale.__FormatOptions*/ constraints){
			// summary:
			//		Formats the value as a Date, according to specified locale (second argument)
			// tags:
			//		protected
			if(!value){ return ''; }
			return this.dateLocaleModule.format(value, constraints);
		},

		"parse": function(/*String*/ value, /*locale.__FormatOptions*/ constraints){
			// summary:
			//		Parses as string as a Date, according to constraints
			// tags:
			//		protected

			return this.dateLocaleModule.parse(value, constraints) || (this._isEmpty(value) ? null : undefined);	 // Date
		},

		// Overrides ValidationTextBox.serialize() to serialize a date in canonical ISO format.
		serialize: function(/*anything*/ val, /*Object?*/ options){
			if(val.toGregorian){
				val = val.toGregorian();
			}
			return stamp.toISOString(val, options);
		},

		// dropDownDefaultValue: Date
		//		The default value to focus in the popupClass widget when the textbox value is empty.
		dropDownDefaultValue : new Date(),

		// value: Date
		//		The value of this widget as a JavaScript Date object.  Use get("value") / set("value", val) to manipulate.
		//		When passed to the parser in markup, must be specified according to `dojo/date/stamp.fromISOString()`
		value: new Date(""),	// value.toString()="NaN"

		_blankValue: null,	// used by filter() when the textbox is blank

		// popupClass: [protected extension] String
		//		Name of the popup widget class used to select a date/time.
		//		Subclasses should specify this.
		popupClass: "", // default is no popup = text only


		// _selector: [protected extension] String
		//		Specifies constraints.selector passed to dojo.date functions, should be either
		//		"date" or "time".
		//		Subclass must specify this.
		_selector: "",

		constructor: function(params /*===== , srcNodeRef =====*/){
			// summary:
			//		Create the widget.
			// params: Object|null
			//		Hash of initialization parameters for widget, including scalar values (like title, duration etc.)
			//		and functions, typically callbacks like onClick.
			//		The hash can contain any of the widget's properties, excluding read-only properties.
			// srcNodeRef: DOMNode|String?
			//		If a srcNodeRef (DOM node) is specified, replace srcNodeRef with my generated DOM tree

			this.dateModule = params.datePackage ? lang.getObject(params.datePackage, false) : date;
			this.dateClassObj = this.dateModule.Date || Date;
			this.dateLocaleModule = params.datePackage ? lang.getObject(params.datePackage+".locale", false) : locale;
			this._set('pattern', this.dateLocaleModule.regexp);
			this._invalidDate = this.constructor.prototype.value.toString();
		},

		buildRendering: function(){
			this.inherited(arguments);

			if(!this.hasDownArrow){
				this._buttonNode.style.display = "none";
			}

			// If hasDownArrow is false, we basically just want to treat the whole widget as the
			// button.
			if(!this.hasDownArrow){
				this._buttonNode = this.domNode;
				this.baseClass += " dijitComboBoxOpenOnClick";
			}
		},

		_setConstraintsAttr: function(/*Object*/ constraints){
			constraints.selector = this._selector;
			constraints.fullYear = true; // see #5465 - always format with 4-digit years
			var fromISO = stamp.fromISOString;
			if(typeof constraints.min == "string"){ constraints.min = fromISO(constraints.min); }
			if(typeof constraints.max == "string"){ constraints.max = fromISO(constraints.max); }
			this.inherited(arguments);
		},

		_isInvalidDate: function(/*Date*/ value){
			// summary:
			//		Runs various tests on the value, checking for invalid conditions
			// tags:
			//		private
			return !value || isNaN(value) || typeof value != "object" || value.toString() == this._invalidDate;
		},

		_setValueAttr: function(/*Date|String*/ value, /*Boolean?*/ priorityChange, /*String?*/ formattedValue){
			// summary:
			//		Sets the date on this textbox. Note: value can be a JavaScript Date literal or a string to be parsed.
			if(value !== undefined){
				if(typeof value == "string"){
					value = stamp.fromISOString(value);
				}
				if(this._isInvalidDate(value)){
					value = null;
				}
				if(value instanceof Date && !(this.dateClassObj instanceof Date)){
					value = new this.dateClassObj(value);
				}
			}
			this.inherited(arguments);
			if(this.value instanceof Date){
				this.filterString = "";
			}
			if(this.dropDown){
				this.dropDown.set('value', value, false);
			}
		},

		_set: function(attr, value){
			// Avoid spurious watch() notifications when value is changed to new Date object w/the same value
			if(attr == "value" && this.value instanceof Date && this.compare(value, this.value) == 0){
				return;
			}
			this.inherited(arguments);
		},

		_setDropDownDefaultValueAttr: function(/*Date*/ val){
			if(this._isInvalidDate(val)){
				// convert null setting into today's date, since there needs to be *some* default at all times.
				 val = new this.dateClassObj();
			}
			this.dropDownDefaultValue = val;
		},

		openDropDown: function(/*Function*/ callback){
			// rebuild drop down every time, so that constraints get copied (#6002)
			if(this.dropDown){
				this.dropDown.destroy();
			}
			var PopupProto = lang.isString(this.popupClass) ? lang.getObject(this.popupClass, false) : this.popupClass,
				textBox = this,
				value = this.get("value");
			this.dropDown = new PopupProto({
				onChange: function(value){
					// this will cause InlineEditBox and other handlers to do stuff so make sure it's last
					textBox.set('value', value, true);
				},
				id: this.id + "_popup",
				dir: textBox.dir,
				lang: textBox.lang,
				value: value,
				currentFocus: !this._isInvalidDate(value) ? value : this.dropDownDefaultValue,
				constraints: textBox.constraints,
				filterString: textBox.filterString, // for TimeTextBox, to filter times shown
				datePackage: textBox.params.datePackage,
				isDisabledDate: function(/*Date*/ date){
					// summary:
					//		disables dates outside of the min/max of the _DateTimeTextBox
					return !textBox.rangeCheck(date, textBox.constraints);
				}
			});

			this.inherited(arguments);
		},

		_getDisplayedValueAttr: function(){
			return this.textbox.value;
		},

		_setDisplayedValueAttr: function(/*String*/ value, /*Boolean?*/ priorityChange){
			this._setValueAttr(this.parse(value, this.constraints), priorityChange, value);
		}
	});


	/*=====
	 _DateTimeTextBox.__Constraints = declare([RangeBoundTextBox.__Constraints, locale.__FormatOptions], {
		 // summary:
		 //		Specifies both the rules on valid/invalid values (first/last date/time allowed),
		 //		and also formatting options for how the date/time is displayed.
		 // example:
		 //		To restrict to dates within 2004, displayed in a long format like "December 25, 2005":
		 //	|		{min:'2004-01-01',max:'2004-12-31', formatLength:'long'}
	 });
	 =====*/

	return _DateTimeTextBox;
});

},
'dijit/form/RangeBoundTextBox':function(){
define("dijit/form/RangeBoundTextBox", [
	"dojo/_base/declare", // declare
	"dojo/i18n", // i18n.getLocalization
	"./MappedTextBox"
], function(declare, i18n, MappedTextBox){

	// module:
	//		dijit/form/RangeBoundTextBox


	var RangeBoundTextBox = declare("dijit.form.RangeBoundTextBox", MappedTextBox, {
		// summary:
		//		Base class for textbox form widgets which defines a range of valid values.

		// rangeMessage: String
		//		The message to display if value is out-of-range
		rangeMessage: "",

		/*=====
		// constraints: RangeBoundTextBox.__Constraints
		constraints: {},
		======*/

		rangeCheck: function(/*Number*/ primitive, /*dijit/form/RangeBoundTextBox.__Constraints*/ constraints){
			// summary:
			//		Overridable function used to validate the range of the numeric input value.
			// tags:
			//		protected
			return	("min" in constraints? (this.compare(primitive,constraints.min) >= 0) : true) &&
				("max" in constraints? (this.compare(primitive,constraints.max) <= 0) : true); // Boolean
		},

		isInRange: function(/*Boolean*/ /*===== isFocused =====*/){
			// summary:
			//		Tests if the value is in the min/max range specified in constraints
			// tags:
			//		protected
			return this.rangeCheck(this.get('value'), this.constraints);
		},

		_isDefinitelyOutOfRange: function(){
			// summary:
			//		Returns true if the value is out of range and will remain
			//		out of range even if the user types more characters
			var val = this.get('value');
			if(val == null){ return false; } // not yet valid enough to compare to
			var outOfRange = false;
			if("min" in this.constraints){
				var min = this.constraints.min;
				outOfRange = this.compare(val, ((typeof min == "number") && min >= 0 && val != 0) ? 0 : min) < 0;
			}
			if(!outOfRange && ("max" in this.constraints)){
				var max = this.constraints.max;
				outOfRange = this.compare(val, ((typeof max != "number") || max > 0) ? max : 0) > 0;
			}
			return outOfRange;
		},

		_isValidSubset: function(){
			// summary:
			//		Overrides `dijit/form/ValidationTextBox._isValidSubset()`.
			//		Returns true if the input is syntactically valid, and either within
			//		range or could be made in range by more typing.
			return this.inherited(arguments) && !this._isDefinitelyOutOfRange();
		},

		isValid: function(/*Boolean*/ isFocused){
			// Overrides dijit/form/ValidationTextBox.isValid() to check that the value is also in range.
			return this.inherited(arguments) &&
				((this._isEmpty(this.textbox.value) && !this.required) || this.isInRange(isFocused)); // Boolean
		},

		getErrorMessage: function(/*Boolean*/ isFocused){
			// Overrides dijit/form/ValidationTextBox.getErrorMessage() to print "out of range" message if appropriate
			var v = this.get('value');
			if(v != null /* and !undefined */ && v !== '' && (typeof v != "number" || !isNaN(v)) && !this.isInRange(isFocused)){ // don't check isInRange w/o a real value
				return this.rangeMessage; // String
			}
			return this.inherited(arguments);
		},

		postMixInProperties: function(){
			this.inherited(arguments);
			if(!this.rangeMessage){
				this.messages = i18n.getLocalization("dijit.form", "validate", this.lang);
				this.rangeMessage = this.messages.rangeMessage;
			}
		},

		applyTextDir: function(/*===== element, text =====*/){
			// summary:
			//		The function overridden in the _BidiSupport module,
			//		originally used for setting element.dir according to this.textDir.
			//		In this case does nothing.
			// element: Object
			// text: String
			// tags:
			//		protected.
		}
	});
	/*=====
	RangeBoundTextBox.__Constraints = declare(null, {
		// min: Number
		//		Minimum signed value.  Default is -Infinity
		// max: Number
		//		Maximum signed value.  Default is +Infinity
	});
	=====*/
	return RangeBoundTextBox;
});

},
'dijit/form/MappedTextBox':function(){
define("dijit/form/MappedTextBox", [
	"dojo/_base/declare", // declare
	"dojo/dom-construct", // domConstruct.place
	"./ValidationTextBox"
], function(declare, domConstruct, ValidationTextBox){

	// module:
	//		dijit/form/MappedTextBox

	return declare("dijit.form.MappedTextBox", ValidationTextBox, {
		// summary:
		//		A dijit/form/ValidationTextBox subclass which provides a base class for widgets that have
		//		a visible formatted display value, and a serializable
		//		value in a hidden input field which is actually sent to the server.
		// description:
		//		The visible display may
		//		be locale-dependent and interactive.  The value sent to the server is stored in a hidden
		//		input field which uses the `name` attribute declared by the original widget.  That value sent
		//		to the server is defined by the dijit/form/MappedTextBox.serialize() method and is typically
		//		locale-neutral.
		// tags:
		//		protected

		postMixInProperties: function(){
			this.inherited(arguments);

			// we want the name attribute to go to the hidden <input>, not the displayed <input>,
			// so override _FormWidget.postMixInProperties() setting of nameAttrSetting
			this.nameAttrSetting = "";
		},

		// Override default behavior to assign name to focusNode
		_setNameAttr: null,

		serialize: function(val /*=====, options =====*/){
			// summary:
			//		Overridable function used to convert the get('value') result to a canonical
			//		(non-localized) string.  For example, will print dates in ISO format, and
			//		numbers the same way as they are represented in javascript.
			// val: anything
			// options: Object?
			// tags:
			//		protected extension
			return val.toString ? val.toString() : ""; // String
		},

		toString: function(){
			// summary:
			//		Returns widget as a printable string using the widget's value
			// tags:
			//		protected
			var val = this.filter(this.get('value')); // call filter in case value is nonstring and filter has been customized
			return val != null ? (typeof val == "string" ? val : this.serialize(val, this.constraints)) : ""; // String
		},

		validate: function(){
			// Overrides `dijit/form/TextBox.validate`
			this.valueNode.value = this.toString();
			return this.inherited(arguments);
		},

		buildRendering: function(){
			// Overrides `dijit/_TemplatedMixin/buildRendering`

			this.inherited(arguments);

			// Create a hidden <input> node with the serialized value used for submit
			// (as opposed to the displayed value).
			// Passing in name as markup rather than calling domConstruct.create() with an attrs argument
			// to make query(input[name=...]) work on IE. (see #8660)
			this.valueNode = domConstruct.place("<input type='hidden'" + (this.name ? ' name="' + this.name.replace(/"/g, "&quot;") + '"' : "") + "/>", this.textbox, "after");
		},

		reset: function(){
			// Overrides `dijit/form/ValidationTextBox.reset` to
			// reset the hidden textbox value to ''
			this.valueNode.value = '';
			this.inherited(arguments);
		}
	});
});

},
'dijit/form/TimeTextBox':function(){
define("dijit/form/TimeTextBox", [
	"dojo/_base/declare", // declare
	"dojo/keys", // keys.DOWN_ARROW keys.ENTER keys.ESCAPE keys.TAB keys.UP_ARROW
	"dojo/_base/lang", // lang.hitch
	"../_TimePicker",
	"./_DateTimeTextBox"
], function(declare, keys, lang, _TimePicker, _DateTimeTextBox){

	// module:
	//		dijit/form/TimeTextBox


	/*=====
	var __Constraints = declare([_DateTimeTextBox.__Constraints, _TimePicker.__Constraints], {
	});
	=====*/

	return declare("dijit.form.TimeTextBox", _DateTimeTextBox, {
		// summary:
		//		A validating, serializable, range-bound time text box with a drop down time picker

		baseClass: "dijitTextBox dijitComboBox dijitTimeTextBox",
		popupClass: _TimePicker,
		_selector: "time",

/*=====
		// constraints: __Constraints
		constraints:{},
=====*/

		// value: Date
		//		The value of this widget as a JavaScript Date object.  Note that the date portion implies time zone and daylight savings rules.
		//
		//		Example:
		// |	new dijit/form/TimeTextBox({value: stamp.fromISOString("T12:59:59", new Date())})
		//
		//		When passed to the parser in markup, must be specified according to locale-independent
		//		`stamp.fromISOString` format.
		//
		//		Example:
		// |	<input data-dojo-type='dijit/form/TimeTextBox' value='T12:34:00'>
		value: new Date(""),		// value.toString()="NaN"
		//FIXME: in markup, you have no control over daylight savings

		_onKey: function(evt){
			if(this.disabled || this.readOnly){ return; }
			this.inherited(arguments);

			// If the user has backspaced or typed some numbers, then filter the result list
			// by what they typed.  Maybe there's a better way to detect this, like _handleOnChange()?
			switch(evt.keyCode){
				case keys.ENTER:
				case keys.TAB:
				case keys.ESCAPE:
				case keys.DOWN_ARROW:
				case keys.UP_ARROW:
					// these keys have special meaning
					break;
				default:
					// defer() because the keystroke hasn't yet appeared in the <input>,
					// so the get('displayedValue') call below won't give the result we want.
					this.defer(function(){
						// set this.filterString to the filter to apply to the drop down list;
						// it will be used in openDropDown()
						var val = this.get('displayedValue');
						this.filterString = (val && !this.parse(val, this.constraints)) ? val.toLowerCase() : "";

						// close the drop down and reopen it, in order to filter the items shown in the list
						// and also since the drop down may need to be repositioned if the number of list items has changed
						// and it's being displayed above the <input>
						if(this._opened){
							this.closeDropDown();
						}
						this.openDropDown();
					});
			}
		}
	});
});

},
'dijit/_TimePicker':function(){
require({cache:{
'url:dijit/templates/TimePicker.html':"<div id=\"widget_${id}\" class=\"dijitMenu\"\n    ><div data-dojo-attach-point=\"upArrow\" class=\"dijitButtonNode dijitUpArrowButton\" data-dojo-attach-event=\"onmouseenter:_buttonMouse,onmouseleave:_buttonMouse\"\n\t\t><div class=\"dijitReset dijitInline dijitArrowButtonInner\" role=\"presentation\">&#160;</div\n\t\t><div class=\"dijitArrowButtonChar\">&#9650;</div></div\n    ><div data-dojo-attach-point=\"timeMenu,focusNode\" data-dojo-attach-event=\"onclick:_onOptionSelected,onmouseover,onmouseout\"></div\n    ><div data-dojo-attach-point=\"downArrow\" class=\"dijitButtonNode dijitDownArrowButton\" data-dojo-attach-event=\"onmouseenter:_buttonMouse,onmouseleave:_buttonMouse\"\n\t\t><div class=\"dijitReset dijitInline dijitArrowButtonInner\" role=\"presentation\">&#160;</div\n\t\t><div class=\"dijitArrowButtonChar\">&#9660;</div></div\n></div>\n"}});
define("dijit/_TimePicker", [
	"dojo/_base/array", // array.forEach
	"dojo/date", // date.compare
	"dojo/date/locale", // locale.format
	"dojo/date/stamp", // stamp.fromISOString stamp.toISOString
	"dojo/_base/declare", // declare
	"dojo/dom-class", // domClass.add domClass.contains domClass.toggle
	"dojo/dom-construct", // domConstruct.create
	"dojo/_base/event", // event.stop
	"dojo/_base/kernel", // deprecated
	"dojo/keys", // keys
	"dojo/_base/lang", // lang.mixin
	"dojo/sniff", // has(...)
	"dojo/query", // query
	"dojo/mouse", // mouse.wheel
	"./typematic",
	"./_Widget",
	"./_TemplatedMixin",
	"./form/_FormValueWidget",
	"dojo/text!./templates/TimePicker.html"
], function(array, ddate, locale, stamp, declare, domClass, domConstruct, event, kernel, keys, lang, has, query, mouse,
			typematic, _Widget, _TemplatedMixin, _FormValueWidget, template){

	// module:
	//		dijit/_TimePicker


	var TimePicker = declare("dijit._TimePicker", [_Widget, _TemplatedMixin], {
		// summary:
		//		A graphical time picker.
		//		This widget is used internally by other widgets and is not available
		//		as a standalone widget due to lack of accessibility support.

		templateString: template,

		// baseClass: [protected] String
		//		The root className to use for the various states of this widget
		baseClass: "dijitTimePicker",

		// clickableIncrement: String
		//		ISO-8601 string representing the amount by which
		//		every clickable element in the time picker increases.
		//		Set in local time, without a time zone.
		//		Example: `T00:15:00` creates 15 minute increments
		//		Must divide dijit/_TimePicker.visibleIncrement evenly
		clickableIncrement: "T00:15:00",

		// visibleIncrement: String
		//		ISO-8601 string representing the amount by which
		//		every element with a visible time in the time picker increases.
		//		Set in local time, without a time zone.
		//		Example: `T01:00:00` creates text in every 1 hour increment
		visibleIncrement: "T01:00:00",

		// visibleRange: String
		//		ISO-8601 string representing the range of this TimePicker.
		//		The TimePicker will only display times in this range.
		//		Example: `T05:00:00` displays 5 hours of options
		visibleRange: "T05:00:00",

		// value: String
		//		Date to display.
		//		Defaults to current time and date.
		//		Can be a Date object or an ISO-8601 string.
		//		If you specify the GMT time zone (`-01:00`),
		//		the time will be converted to the local time in the local time zone.
		//		Otherwise, the time is considered to be in the local time zone.
		//		If you specify the date and isDate is true, the date is used.
		//		Example: if your local time zone is `GMT -05:00`,
		//		`T10:00:00` becomes `T10:00:00-05:00` (considered to be local time),
		//		`T10:00:00-01:00` becomes `T06:00:00-05:00` (4 hour difference),
		//		`T10:00:00Z` becomes `T05:00:00-05:00` (5 hour difference between Zulu and local time)
		//		`yyyy-mm-ddThh:mm:ss` is the format to set the date and time
		//		Example: `2007-06-01T09:00:00`
		value: new Date(),

		_visibleIncrement:2,
		_clickableIncrement:1,
		_totalIncrements:10,

		// constraints: TimePicker.__Constraints
		//		Specifies valid range of times (start time, end time)
		constraints:{},

/*=====
		serialize: function(val, options){
			// summary:
			//		User overridable function used to convert the attr('value') result to a String
			// val: Date
			//		The current value
			// options: Object?
			// tags:
			//		protected
		},
=====*/
		serialize: stamp.toISOString,

/*=====
		// filterString: string
		//		The string to filter by
		filterString: "",
=====*/

		setValue: function(/*Date*/ value){
			// summary:
			//		Deprecated.  Used set('value') instead.
			// tags:
			//		deprecated
			kernel.deprecated("dijit._TimePicker:setValue() is deprecated.  Use set('value', ...) instead.", "", "2.0");
			this.set('value', value);
		},

		_setValueAttr: function(/*Date*/ date){
			// summary:
			//		Hook so set('value', ...) works.
			// description:
			//		Set the value of the TimePicker.
			//		Redraws the TimePicker around the new date.
			// tags:
			//		protected
			this._set("value", date);
			this._showText();
		},

		_setFilterStringAttr: function(val){
			// summary:
			//		Called by TimeTextBox to filter the values shown in my list
			this._set("filterString", val);
			this._showText();
		},

		isDisabledDate: function(/*===== dateObject, locale =====*/){
			// summary:
			//		May be overridden to disable certain dates in the TimePicker e.g. `isDisabledDate=locale.isWeekend`
			// dateObject: Date
			// locale: String?
			// type:
			//		extension
			return false; // Boolean
		},

		_getFilteredNodes: function(/*number*/ start, /*number*/ maxNum, /*Boolean*/ before, /*DOMnode*/ lastNode){
			// summary:
			//		Returns an array of nodes with the filter applied.  At most maxNum nodes
			//		will be returned - but fewer may be returned as well.  If the
			//		before parameter is set to true, then it will return the elements
			//		before the given index
			// tags:
			//		private
			var
				nodes = [],
				lastValue = lastNode ? lastNode.date : this._refDate,
				n,
				i = start,
				max = this._maxIncrement + Math.abs(i),
				chk = before ? -1 : 1,
				dec = before ? 1 : 0,
				inc = 1 - dec;
			do{
				i -= dec;
				n = this._createOption(i);
				if(n){
					if((before && n.date > lastValue) || (!before && n.date < lastValue)){
						break; // don't wrap
					}
					nodes[before ? "unshift" : "push"](n);
					lastValue = n.date;
				}
				i += inc;
			}while(nodes.length < maxNum && (i*chk) < max);
			return nodes;
		},

		_showText: function(){
			// summary:
			//		Displays the relevant choices in the drop down list
			// tags:
			//		private
			var fromIso = stamp.fromISOString;
			this.timeMenu.innerHTML = "";
			this._clickableIncrementDate=fromIso(this.clickableIncrement);
			this._visibleIncrementDate=fromIso(this.visibleIncrement);
			this._visibleRangeDate=fromIso(this.visibleRange);
			// get the value of the increments and the range in seconds (since 00:00:00) to find out how many divs to create
			var
				sinceMidnight = function(/*Date*/ date){
					return date.getHours() * 60 * 60 + date.getMinutes() * 60 + date.getSeconds();
				},
				clickableIncrementSeconds = sinceMidnight(this._clickableIncrementDate),
				visibleIncrementSeconds = sinceMidnight(this._visibleIncrementDate),
				visibleRangeSeconds = sinceMidnight(this._visibleRangeDate),
				// round reference date to previous visible increment
				time = (this.value || this.currentFocus).getTime();

			this._refDate = new Date(time - time % (clickableIncrementSeconds*1000));
			this._refDate.setFullYear(1970,0,1); // match parse defaults

			// assume clickable increment is the smallest unit
			this._clickableIncrement = 1;
			// divide the visible range by the clickable increment to get the number of divs to create
			// example: 10:00:00/00:15:00 -> display 40 divs
			this._totalIncrements = visibleRangeSeconds / clickableIncrementSeconds;
			// divide the visible increments by the clickable increments to get how often to display the time inline
			// example: 01:00:00/00:15:00 -> display the time every 4 divs
			this._visibleIncrement = visibleIncrementSeconds / clickableIncrementSeconds;
			// divide the number of seconds in a day by the clickable increment in seconds to get the
			// absolute max number of increments.
			this._maxIncrement = (60 * 60 * 24) / clickableIncrementSeconds;

			var
				// Find the nodes we should display based on our filter.
				// Limit to 10 nodes displayed as a half-hearted attempt to stop drop down from overlapping <input>.
				count = Math.min(this._totalIncrements, 10),
				after = this._getFilteredNodes(0, (count >> 1) + 1, false),
				moreAfter = [],
				estBeforeLength = count - after.length,
				before = this._getFilteredNodes(0, estBeforeLength, true, after[0]);
			if(before.length < estBeforeLength && after.length > 0){
				moreAfter = this._getFilteredNodes(after[after.length-1].idx + 1, estBeforeLength - before.length, false, after[after.length-1]);
			}
			array.forEach(before.concat(after, moreAfter), function(n){ this.timeMenu.appendChild(n); }, this);
			// never show empty due to a bad filter
			if(!before.length && !after.length && !moreAfter.length && this.filterString){
				this.filterString = '';
				this._showText();
			}
		},

		constructor: function(/*===== params, srcNodeRef =====*/){
			// summary:
			//		Create the widget.
			// params: Object|null
			//		Hash of initialization parameters for widget, including scalar values (like title, duration etc.)
			//		and functions, typically callbacks like onClick.
			//		The hash can contain any of the widget's properties, excluding read-only properties.
			// srcNodeRef: DOMNode|String?
			//		If a srcNodeRef (DOM node) is specified, replace srcNodeRef with my generated DOM tree

			this.constraints = {};
		},

		postMixInProperties: function(){
			this.inherited(arguments);
			this._setConstraintsAttr(this.constraints); // this needs to happen now (and later) due to codependency on _set*Attr calls
		},

		_setConstraintsAttr: function(/* Object */ constraints){
			// brings in visibleRange, increments, etc.
			lang.mixin(this, constraints);

			// locale needs the lang in the constraints as locale
			if(!constraints.locale){
				constraints.locale = this.lang;
			}
		},

		postCreate: function(){
			// assign typematic mouse listeners to the arrow buttons
			this.connect(this.timeMenu, mouse.wheel, "_mouseWheeled");
			this.own(
				typematic.addMouseListener(this.upArrow, this, "_onArrowUp", 33, 250),
				typematic.addMouseListener(this.downArrow, this, "_onArrowDown", 33, 250)
			);

			this.inherited(arguments);
		},

		_buttonMouse: function(/*Event*/ e){
			// summary:
			//		Handler for hover (and unhover) on up/down arrows
			// tags:
			//		private

			// in non-IE browser the "mouseenter" event will become "mouseover",
			// but in IE it's still "mouseenter"
			domClass.toggle(e.currentTarget, e.currentTarget == this.upArrow ? "dijitUpArrowHover" : "dijitDownArrowHover",
				e.type == "mouseenter" || e.type == "mouseover");
		},

		_createOption: function(/*Number*/ index){
			// summary:
			//		Creates a clickable time option
			// tags:
			//		private
			var date = new Date(this._refDate);
			var incrementDate = this._clickableIncrementDate;
			date.setTime(date.getTime()
				+ incrementDate.getHours() * index * 3600000
				+ incrementDate.getMinutes() * index * 60000
				+ incrementDate.getSeconds() * index * 1000);
			if(this.constraints.selector == "time"){
				date.setFullYear(1970,0,1); // make sure each time is for the same date
			}
			var dateString = locale.format(date, this.constraints);
			if(this.filterString && dateString.toLowerCase().indexOf(this.filterString) !== 0){
				// Doesn't match the filter - return null
				return null;
			}

			var div = this.ownerDocument.createElement("div");
			div.className = this.baseClass+"Item";
			div.date = date;
			div.idx = index;
			domConstruct.create('div',{
				"class": this.baseClass + "ItemInner",
				innerHTML: dateString
			}, div);

			if(index%this._visibleIncrement<1 && index%this._visibleIncrement>-1){
				domClass.add(div, this.baseClass+"Marker");
			}else if(!(index%this._clickableIncrement)){
				domClass.add(div, this.baseClass+"Tick");
			}

			if(this.isDisabledDate(date)){
				// set disabled
				domClass.add(div, this.baseClass+"ItemDisabled");
			}
			if(this.value && !ddate.compare(this.value, date, this.constraints.selector)){
				div.selected = true;
				domClass.add(div, this.baseClass+"ItemSelected");
				if(domClass.contains(div, this.baseClass+"Marker")){
					domClass.add(div, this.baseClass+"MarkerSelected");
				}else{
					domClass.add(div, this.baseClass+"TickSelected");
				}

				// Initially highlight the current value.   User can change highlight by up/down arrow keys
				// or mouse movement.
				this._highlightOption(div, true);
			}
			return div;
		},

		_onOptionSelected: function(/*Object*/ tgt){
			// summary:
			//		Called when user clicks an option in the drop down list
			// tags:
			//		private
			var tdate = tgt.target.date || tgt.target.parentNode.date;
			if(!tdate || this.isDisabledDate(tdate)){ return; }
			this._highlighted_option = null;
			this.set('value', tdate);
			this.onChange(tdate);
		},

		onChange: function(/*Date*/ /*===== time =====*/){
			// summary:
			//		Notification that a time was selected.  It may be the same as the previous value.
			// tags:
			//		public
		},

		_highlightOption: function(/*node*/ node, /*Boolean*/ highlight){
			// summary:
			//		Turns on/off highlight effect on a node based on mouse out/over event
			// tags:
			//		private
			if(!node){return;}
			if(highlight){
				if(this._highlighted_option){
					this._highlightOption(this._highlighted_option, false);
				}
				this._highlighted_option = node;
			}else if(this._highlighted_option !== node){
				return;
			}else{
				this._highlighted_option = null;
			}
			domClass.toggle(node, this.baseClass+"ItemHover", highlight);
			if(domClass.contains(node, this.baseClass+"Marker")){
				domClass.toggle(node, this.baseClass+"MarkerHover", highlight);
			}else{
				domClass.toggle(node, this.baseClass+"TickHover", highlight);
			}
		},

		onmouseover: function(/*Event*/ e){
			// summary:
			//		Handler for onmouseover event
			// tags:
			//		private
			this._keyboardSelected = null;
			var tgr = (e.target.parentNode === this.timeMenu) ? e.target : e.target.parentNode;
			// if we aren't targeting an item, then we return
			if(!domClass.contains(tgr, this.baseClass+"Item")){return;}
			this._highlightOption(tgr, true);
		},

		onmouseout: function(/*Event*/ e){
			// summary:
			//		Handler for onmouseout event
			// tags:
			//		private
			this._keyboardSelected = null;
			var tgr = (e.target.parentNode === this.timeMenu) ? e.target : e.target.parentNode;
			this._highlightOption(tgr, false);
		},

		_mouseWheeled: function(/*Event*/ e){
			// summary:
			//		Handle the mouse wheel events
			// tags:
			//		private
			this._keyboardSelected = null;
			event.stop(e);
			// we're not _measuring_ the scroll amount, just direction
			this[(e.wheelDelta>0 ? "_onArrowUp" : "_onArrowDown")](); // yes, we're making a new dom node every time you mousewheel, or click
		},

		_onArrowUp: function(count){
			// summary:
			//		Handler for up arrow key.
			// description:
			//		Removes the bottom time and add one to the top
			// tags:
			//		private
			if(count === -1){
				domClass.remove(this.upArrow, "dijitUpArrowActive");
				return;
			}else if(count === 0){
				domClass.add(this.upArrow, "dijitUpArrowActive");

			} // typematic end
			if(!this.timeMenu.childNodes.length){ return; }
			var index = this.timeMenu.childNodes[0].idx;
			var divs = this._getFilteredNodes(index, 1, true, this.timeMenu.childNodes[0]);
			if(divs.length){
				this.timeMenu.removeChild(this.timeMenu.childNodes[this.timeMenu.childNodes.length - 1]);
				this.timeMenu.insertBefore(divs[0], this.timeMenu.childNodes[0]);
			}
		},

		_onArrowDown: function(count){
			// summary:
			//		Handler for up arrow key.
			// description:
			//		Remove the top time and add one to the bottom
			// tags:
			//		private
			if(count === -1){
				domClass.remove(this.downArrow, "dijitDownArrowActive");
				return;
			}else if(count === 0){
				domClass.add(this.downArrow, "dijitDownArrowActive");
			} // typematic end
			if(!this.timeMenu.childNodes.length){ return; }
			var index = this.timeMenu.childNodes[this.timeMenu.childNodes.length - 1].idx + 1;
			var divs = this._getFilteredNodes(index, 1, false, this.timeMenu.childNodes[this.timeMenu.childNodes.length - 1]);
			if(divs.length){
				this.timeMenu.removeChild(this.timeMenu.childNodes[0]);
				this.timeMenu.appendChild(divs[0]);
			}
		},

		handleKey: function(/*Event*/ e){
			// summary:
			//		Called from `dijit/form/_DateTimeTextBox` to pass a keypress event
			//		from the `dijit/form/TimeTextBox` to be handled in this widget
			// tags:
			//		protected
			if(e.keyCode == keys.DOWN_ARROW || e.keyCode == keys.UP_ARROW){
				event.stop(e);
				// Figure out which option to highlight now and then highlight it
				if(this._highlighted_option && !this._highlighted_option.parentNode){
					this._highlighted_option = null;
				}
				var timeMenu = this.timeMenu,
					tgt = this._highlighted_option || query("." + this.baseClass + "ItemSelected", timeMenu)[0];
				if(!tgt){
					tgt = timeMenu.childNodes[0];
				}else if(timeMenu.childNodes.length){
					if(e.keyCode == keys.DOWN_ARROW && !tgt.nextSibling){
						this._onArrowDown();
					}else if(e.keyCode == keys.UP_ARROW && !tgt.previousSibling){
						this._onArrowUp();
					}
					if(e.keyCode == keys.DOWN_ARROW){
						tgt = tgt.nextSibling;
					}else{
						tgt = tgt.previousSibling;
					}
				}
				this._highlightOption(tgt, true);
				this._keyboardSelected = tgt;
				return false;
			}else if(e.keyCode == keys.ENTER || e.keyCode === keys.TAB){
				// mouse hover followed by TAB is NO selection
				if(!this._keyboardSelected && e.keyCode === keys.TAB){
					return true;	// true means don't call stopEvent()
				}

				// Accept the currently-highlighted option as the value
				if(this._highlighted_option){
					this._onOptionSelected({target: this._highlighted_option});
				}

				// Call stopEvent() for ENTER key so that form doesn't submit,
				// but not for TAB, so that TAB does switch focus
				return e.keyCode === keys.TAB;
			}
			return undefined;
		}
	});

	/*=====
	 TimePicker.__Constraints = declare(locale.__FormatOptions, {
		 // clickableIncrement: String
		 //		See `dijit/_TimePicker.clickableIncrement`
		 clickableIncrement: "T00:15:00",

		 // visibleIncrement: String
		 //		See `dijit/_TimePicker.visibleIncrement`
		 visibleIncrement: "T01:00:00",

		 // visibleRange: String
		 //		See `dijit/_TimePicker.visibleRange`
		 visibleRange: "T05:00:00"
	 });
	 =====*/

	return TimePicker;
});

},
'dijit/typematic':function(){
define("dijit/typematic", [
	"dojo/_base/array", // array.forEach
	"dojo/_base/connect", // connect.connect
	"dojo/_base/event", // event.stop
	"dojo/_base/kernel", // kernel.deprecated
	"dojo/_base/lang", // lang.mixin, lang.hitch
	"dojo/on",
	"dojo/sniff", // has("ie")
	"./main"		// setting dijit.typematic global
], function(array, connect, event, kernel, lang, on, has, dijit){

// module:
//		dijit/typematic

var typematic = (dijit.typematic = {
	// summary:
	//		These functions are used to repetitively call a user specified callback
	//		method when a specific key or mouse click over a specific DOM node is
	//		held down for a specific amount of time.
	//		Only 1 such event is allowed to occur on the browser page at 1 time.

	_fireEventAndReload: function(){
		this._timer = null;
		this._callback(++this._count, this._node, this._evt);

		// Schedule next event, timer is at most minDelay (default 10ms) to avoid
		// browser overload (particularly avoiding starving DOH robot so it never gets to send a mouseup)
		this._currentTimeout = Math.max(
			this._currentTimeout < 0 ? this._initialDelay :
				(this._subsequentDelay > 1 ? this._subsequentDelay : Math.round(this._currentTimeout * this._subsequentDelay)),
			this._minDelay);
		this._timer = setTimeout(lang.hitch(this, "_fireEventAndReload"), this._currentTimeout);
	},

	trigger: function(/*Event*/ evt, /*Object*/ _this, /*DOMNode*/ node, /*Function*/ callback, /*Object*/ obj, /*Number?*/ subsequentDelay, /*Number?*/ initialDelay, /*Number?*/ minDelay){
		// summary:
		//		Start a timed, repeating callback sequence.
		//		If already started, the function call is ignored.
		//		This method is not normally called by the user but can be
		//		when the normal listener code is insufficient.
		// evt:
		//		key or mouse event object to pass to the user callback
		// _this:
		//		pointer to the user's widget space.
		// node:
		//		the DOM node object to pass the the callback function
		// callback:
		//		function to call until the sequence is stopped called with 3 parameters:
		// count:
		//		integer representing number of repeated calls (0..n) with -1 indicating the iteration has stopped
		// node:
		//		the DOM node object passed in
		// evt:
		//		key or mouse event object
		// obj:
		//		user space object used to uniquely identify each typematic sequence
		// subsequentDelay:
		//		if > 1, the number of milliseconds until the 3->n events occur
		//		or else the fractional time multiplier for the next event's delay, default=0.9
		// initialDelay:
		//		the number of milliseconds until the 2nd event occurs, default=500ms
		// minDelay:
		//		the maximum delay in milliseconds for event to fire, default=10ms
		if(obj != this._obj){
			this.stop();
			this._initialDelay = initialDelay || 500;
			this._subsequentDelay = subsequentDelay || 0.90;
			this._minDelay = minDelay || 10;
			this._obj = obj;
			this._node = node;
			this._currentTimeout = -1;
			this._count = -1;
			this._callback = lang.hitch(_this, callback);
			this._evt = { faux: true };
			for(var attr in evt){
				if(attr != "layerX" && attr != "layerY"){ // prevent WebKit warnings
					var v = evt[attr];
					if(typeof v != "function" && typeof v != "undefined"){ this._evt[attr] = v }
				}
			}
			this._fireEventAndReload();
		}
	},

	stop: function(){
		// summary:
		//		Stop an ongoing timed, repeating callback sequence.
		if(this._timer){
			clearTimeout(this._timer);
			this._timer = null;
		}
		if(this._obj){
			this._callback(-1, this._node, this._evt);
			this._obj = null;
		}
	},

	addKeyListener: function(/*DOMNode*/ node, /*Object*/ keyObject, /*Object*/ _this, /*Function*/ callback, /*Number*/ subsequentDelay, /*Number*/ initialDelay, /*Number?*/ minDelay){
		// summary:
		//		Start listening for a specific typematic key.
		//		See also the trigger method for other parameters.
		// keyObject:
		//		an object defining the key to listen for:
		//
		//		- charOrCode: the printable character (string) or keyCode (number) to listen for.
		//		- keyCode: (deprecated - use charOrCode) the keyCode (number) to listen for (implies charCode = 0).
		//		- charCode: (deprecated - use charOrCode) the charCode (number) to listen for.
		//		- ctrlKey: desired ctrl key state to initiate the callback sequence:
		//			- pressed (true)
		//			- released (false)
		//			- either (unspecified)
		//		- altKey: same as ctrlKey but for the alt key
		//		- shiftKey: same as ctrlKey but for the shift key
		// returns:
		//		a connection handle

		if(keyObject.keyCode){
			keyObject.charOrCode = keyObject.keyCode;
			kernel.deprecated("keyCode attribute parameter for dijit.typematic.addKeyListener is deprecated. Use charOrCode instead.", "", "2.0");
		}else if(keyObject.charCode){
			keyObject.charOrCode = String.fromCharCode(keyObject.charCode);
			kernel.deprecated("charCode attribute parameter for dijit.typematic.addKeyListener is deprecated. Use charOrCode instead.", "", "2.0");
		}
		var handles = [
			on(node, connect._keypress, lang.hitch(this, function(evt){
				if(evt.charOrCode == keyObject.charOrCode &&
				(keyObject.ctrlKey === undefined || keyObject.ctrlKey == evt.ctrlKey) &&
				(keyObject.altKey === undefined || keyObject.altKey == evt.altKey) &&
				(keyObject.metaKey === undefined || keyObject.metaKey == (evt.metaKey || false)) && // IE doesn't even set metaKey
				(keyObject.shiftKey === undefined || keyObject.shiftKey == evt.shiftKey)){
					event.stop(evt);
					typematic.trigger(evt, _this, node, callback, keyObject, subsequentDelay, initialDelay, minDelay);
				}else if(typematic._obj == keyObject){
					typematic.stop();
				}
			})),
			on(node, "keyup", lang.hitch(this, function(){
				if(typematic._obj == keyObject){
					typematic.stop();
				}
			}))
		];
		return { remove: function(){ array.forEach(handles, function(h){ h.remove(); }); } };
	},

	addMouseListener: function(/*DOMNode*/ node, /*Object*/ _this, /*Function*/ callback, /*Number*/ subsequentDelay, /*Number*/ initialDelay, /*Number?*/ minDelay){
		// summary:
		//		Start listening for a typematic mouse click.
		//		See the trigger method for other parameters.
		// returns:
		//		a connection handle
		var handles = [
			on(node, "mousedown", lang.hitch(this, function(evt){
				evt.preventDefault();
				typematic.trigger(evt, _this, node, callback, node, subsequentDelay, initialDelay, minDelay);
			})),
			on(node, "mouseup", lang.hitch(this, function(evt){
				if(this._obj){
					evt.preventDefault();
				}
				typematic.stop();
			})),
			on(node, "mouseout", lang.hitch(this, function(evt){
				if(this._obj){
					evt.preventDefault();
				}
				typematic.stop();
			})),
			on(node, "dblclick", lang.hitch(this, function(evt){
				evt.preventDefault();
				if(has("ie") < 9){
					typematic.trigger(evt, _this, node, callback, node, subsequentDelay, initialDelay, minDelay);
					setTimeout(lang.hitch(this, typematic.stop), 50);
				}
			}))
		];
		return { remove: function(){ array.forEach(handles, function(h){ h.remove(); }); } };
	},

	addListener: function(/*Node*/ mouseNode, /*Node*/ keyNode, /*Object*/ keyObject, /*Object*/ _this, /*Function*/ callback, /*Number*/ subsequentDelay, /*Number*/ initialDelay, /*Number?*/ minDelay){
		// summary:
		//		Start listening for a specific typematic key and mouseclick.
		//		This is a thin wrapper to addKeyListener and addMouseListener.
		//		See the addMouseListener and addKeyListener methods for other parameters.
		// mouseNode:
		//		the DOM node object to listen on for mouse events.
		// keyNode:
		//		the DOM node object to listen on for key events.
		// returns:
		//		a connection handle
		var handles = [
			this.addKeyListener(keyNode, keyObject, _this, callback, subsequentDelay, initialDelay, minDelay),
			this.addMouseListener(mouseNode, _this, callback, subsequentDelay, initialDelay, minDelay)
		];
		return { remove: function(){ array.forEach(handles, function(h){ h.remove(); }); } };
	}
});

return typematic;

});

},
'url:dijit/templates/TimePicker.html':"<div id=\"widget_${id}\" class=\"dijitMenu\"\n    ><div data-dojo-attach-point=\"upArrow\" class=\"dijitButtonNode dijitUpArrowButton\" data-dojo-attach-event=\"onmouseenter:_buttonMouse,onmouseleave:_buttonMouse\"\n\t\t><div class=\"dijitReset dijitInline dijitArrowButtonInner\" role=\"presentation\">&#160;</div\n\t\t><div class=\"dijitArrowButtonChar\">&#9650;</div></div\n    ><div data-dojo-attach-point=\"timeMenu,focusNode\" data-dojo-attach-event=\"onclick:_onOptionSelected,onmouseover,onmouseout\"></div\n    ><div data-dojo-attach-point=\"downArrow\" class=\"dijitButtonNode dijitDownArrowButton\" data-dojo-attach-event=\"onmouseenter:_buttonMouse,onmouseleave:_buttonMouse\"\n\t\t><div class=\"dijitReset dijitInline dijitArrowButtonInner\" role=\"presentation\">&#160;</div\n\t\t><div class=\"dijitArrowButtonChar\">&#9660;</div></div\n></div>\n",
'dijit/form/NumberSpinner':function(){
define("dijit/form/NumberSpinner", [
	"dojo/_base/declare", // declare
	"dojo/_base/event", // event.stop
	"dojo/keys", // keys.END keys.HOME
	"./_Spinner",
	"./NumberTextBox"
], function(declare, event, keys, _Spinner, NumberTextBox){

// module:
//		dijit/form/NumberSpinner


return declare("dijit.form.NumberSpinner", [_Spinner, NumberTextBox.Mixin], {
	// summary:
	//		Extends NumberTextBox to add up/down arrows and pageup/pagedown for incremental change to the value
	//
	// description:
	//		A `dijit/form/NumberTextBox` extension to provide keyboard accessible value selection
	//		as well as icons for spinning direction. When using the keyboard, the typematic rules
	//		apply, meaning holding the key will gradually increase or decrease the value and
	//		accelerate.
	//
	// example:
	//	| new NumberSpinner({ constraints:{ max:300, min:100 }}, "someInput");

	baseClass: "dijitTextBox dijitSpinner dijitNumberTextBox",

	adjust: function(/*Object*/ val, /*Number*/ delta){
		// summary:
		//		Change Number val by the given amount
		// tags:
		//		protected

		var tc = this.constraints,
			v = isNaN(val),
			gotMax = !isNaN(tc.max),
			gotMin = !isNaN(tc.min)
		;
		if(v && delta != 0){ // blank or invalid value and they want to spin, so create defaults
			val = (delta > 0) ?
				gotMin ? tc.min : gotMax ? tc.max : 0 :
				gotMax ? this.constraints.max : gotMin ? tc.min : 0
			;
		}
		var newval = val + delta;
		if(v || isNaN(newval)){ return val; }
		if(gotMax && (newval > tc.max)){
			newval = tc.max;
		}
		if(gotMin && (newval < tc.min)){
			newval = tc.min;
		}
		return newval;
	},

	_onKeyPress: function(e){
		if((e.charOrCode == keys.HOME || e.charOrCode == keys.END) && !(e.ctrlKey || e.altKey || e.metaKey)
		&& typeof this.get('value') != 'undefined' /* gibberish, so HOME and END are default editing keys*/){
			var value = this.constraints[(e.charOrCode == keys.HOME ? "min" : "max")];
			if(typeof value == "number"){
				this._setValueAttr(value, false);
			}
			// eat home or end key whether we change the value or not
			event.stop(e);
		}
	}
});

});

},
'dijit/form/_Spinner':function(){
require({cache:{
'url:dijit/form/templates/Spinner.html':"<div class=\"dijit dijitReset dijitInline dijitLeft\"\n\tid=\"widget_${id}\" role=\"presentation\"\n\t><div class=\"dijitReset dijitButtonNode dijitSpinnerButtonContainer\"\n\t\t><input class=\"dijitReset dijitInputField dijitSpinnerButtonInner\" type=\"text\" tabIndex=\"-1\" readonly=\"readonly\" role=\"presentation\"\n\t\t/><div class=\"dijitReset dijitLeft dijitButtonNode dijitArrowButton dijitUpArrowButton\"\n\t\t\tdata-dojo-attach-point=\"upArrowNode\"\n\t\t\t><div class=\"dijitArrowButtonInner\"\n\t\t\t\t><input class=\"dijitReset dijitInputField\" value=\"&#9650; \" type=\"text\" tabIndex=\"-1\" readonly=\"readonly\" role=\"presentation\"\n\t\t\t\t\t${_buttonInputDisabled}\n\t\t\t/></div\n\t\t></div\n\t\t><div class=\"dijitReset dijitLeft dijitButtonNode dijitArrowButton dijitDownArrowButton\"\n\t\t\tdata-dojo-attach-point=\"downArrowNode\"\n\t\t\t><div class=\"dijitArrowButtonInner\"\n\t\t\t\t><input class=\"dijitReset dijitInputField\" value=\"&#9660; \" type=\"text\" tabIndex=\"-1\" readonly=\"readonly\" role=\"presentation\"\n\t\t\t\t\t${_buttonInputDisabled}\n\t\t\t/></div\n\t\t></div\n\t></div\n\t><div class='dijitReset dijitValidationContainer'\n\t\t><input class=\"dijitReset dijitInputField dijitValidationIcon dijitValidationInner\" value=\"&#935; \" type=\"text\" tabIndex=\"-1\" readonly=\"readonly\" role=\"presentation\"\n\t/></div\n\t><div class=\"dijitReset dijitInputField dijitInputContainer\"\n\t\t><input class='dijitReset dijitInputInner' data-dojo-attach-point=\"textbox,focusNode\" type=\"${type}\" data-dojo-attach-event=\"onkeypress:_onKeyPress\"\n\t\t\trole=\"spinbutton\" autocomplete=\"off\" ${!nameAttrSetting}\n\t/></div\n></div>\n"}});
define("dijit/form/_Spinner", [
	"dojo/_base/declare", // declare
	"dojo/_base/event", // event.stop
	"dojo/keys", // keys keys.DOWN_ARROW keys.PAGE_DOWN keys.PAGE_UP keys.UP_ARROW
	"dojo/_base/lang", // lang.hitch
	"dojo/sniff", // has("mozilla")
	"dojo/mouse", // mouse.wheel
	"../typematic",
	"./RangeBoundTextBox",
	"dojo/text!./templates/Spinner.html",
	"./_TextBoxMixin"	// selectInputText
], function(declare, event, keys, lang, has, mouse, typematic, RangeBoundTextBox, template, _TextBoxMixin){

	// module:
	//		dijit/form/_Spinner

	return declare("dijit.form._Spinner", RangeBoundTextBox, {
		// summary:
		//		Mixin for validation widgets with a spinner.
		// description:
		//		This class basically (conceptually) extends `dijit/form/ValidationTextBox`.
		//		It modifies the template to have up/down arrows, and provides related handling code.

		// defaultTimeout: Number
		//		Number of milliseconds before a held arrow key or up/down button becomes typematic
		defaultTimeout: 500,

		// minimumTimeout: Number
		//		minimum number of milliseconds that typematic event fires when held key or button is held
		minimumTimeout: 10,

		// timeoutChangeRate: Number
		//		Fraction of time used to change the typematic timer between events.
		//		1.0 means that each typematic event fires at defaultTimeout intervals.
		//		Less than 1.0 means that each typematic event fires at an increasing faster rate.
		timeoutChangeRate: 0.90,

		// smallDelta: Number
		//		Adjust the value by this much when spinning using the arrow keys/buttons
		smallDelta: 1,

		// largeDelta: Number
		//		Adjust the value by this much when spinning using the PgUp/Dn keys
		largeDelta: 10,

		templateString: template,

		baseClass: "dijitTextBox dijitSpinner",

		// Set classes like dijitUpArrowButtonHover or dijitDownArrowButtonActive depending on
		// mouse action over specified node
		cssStateNodes: {
			"upArrowNode": "dijitUpArrowButton",
			"downArrowNode": "dijitDownArrowButton"
		},

		adjust: function(val /*=====, delta =====*/){
			// summary:
			//		Overridable function used to adjust a primitive value(Number/Date/...) by the delta amount specified.
			//		The val is adjusted in a way that makes sense to the object type.
			// val: Object
			// delta: Number
			// tags:
			//		protected extension
			return val;
		},

		_arrowPressed: function(/*Node*/ nodePressed, /*Number*/ direction, /*Number*/ increment){
			// summary:
			//		Handler for arrow button or arrow key being pressed
			if(this.disabled || this.readOnly){ return; }
			this._setValueAttr(this.adjust(this.get('value'), direction*increment), false);
			_TextBoxMixin.selectInputText(this.textbox, this.textbox.value.length);
		},

		_arrowReleased: function(/*Node*/ /*===== node =====*/){
			// summary:
			//		Handler for arrow button or arrow key being released
			this._wheelTimer = null;
		},

		_typematicCallback: function(/*Number*/ count, /*DOMNode*/ node, /*Event*/ evt){
			var inc=this.smallDelta;
			if(node == this.textbox){
				var key = evt.charOrCode;
				inc = (key == keys.PAGE_UP || key == keys.PAGE_DOWN) ? this.largeDelta : this.smallDelta;
				node = (key == keys.UP_ARROW || key == keys.PAGE_UP) ? this.upArrowNode : this.downArrowNode;
			}
			if(count == -1){ this._arrowReleased(node); }
			else{ this._arrowPressed(node, (node == this.upArrowNode) ? 1 : -1, inc); }
		},

		_wheelTimer: null,
		_mouseWheeled: function(/*Event*/ evt){
			// summary:
			//		Mouse wheel listener where supported

			event.stop(evt);
			// FIXME: Safari bubbles

			// be nice to DOH and scroll as much as the event says to
			var wheelDelta = evt.wheelDelta / 120;
			if(Math.floor(wheelDelta) != wheelDelta){
				// If not an int multiple of 120, then its touchpad scrolling.
				// This can change very fast so just assume 1 wheel click to make it more manageable.
				wheelDelta = evt.wheelDelta > 0 ? 1 : -1;
			}
			var scrollAmount = evt.detail ? (evt.detail * -1) : wheelDelta;
			if(scrollAmount !== 0){
				var node = this[(scrollAmount > 0 ? "upArrowNode" : "downArrowNode" )];

				this._arrowPressed(node, scrollAmount, this.smallDelta);

				if(this._wheelTimer){
					this._wheelTimer.remove();
				}
				this._wheelTimer = this.defer(function(){ this._arrowReleased(node); }, 50);
			}
		},

		_setConstraintsAttr: function(/*Object*/ constraints){
			this.inherited(arguments);
			if(this.focusNode){ // not set when called from postMixInProperties
				if(this.constraints.min !== undefined){
					this.focusNode.setAttribute("aria-valuemin", this.constraints.min);
				}else{
					this.focusNode.removeAttribute("aria-valuemin");
				}
				if(this.constraints.max !== undefined){
					this.focusNode.setAttribute("aria-valuemax", this.constraints.max);
				}else{
					this.focusNode.removeAttribute("aria-valuemax");
				}
			}
		},

		_setValueAttr: function(/*Number*/ value, /*Boolean?*/ priorityChange){
			// summary:
			//		Hook so set('value', ...) works.

			this.focusNode.setAttribute("aria-valuenow", value);
			this.inherited(arguments);
		},

		postCreate: function(){
			this.inherited(arguments);

			// extra listeners
			this.connect(this.domNode, mouse.wheel, "_mouseWheeled");
			this.own(
				typematic.addListener(this.upArrowNode, this.textbox, {charOrCode:keys.UP_ARROW,ctrlKey:false,altKey:false,shiftKey:false,metaKey:false}, this, "_typematicCallback", this.timeoutChangeRate, this.defaultTimeout, this.minimumTimeout),
				typematic.addListener(this.downArrowNode, this.textbox, {charOrCode:keys.DOWN_ARROW,ctrlKey:false,altKey:false,shiftKey:false,metaKey:false}, this, "_typematicCallback", this.timeoutChangeRate, this.defaultTimeout, this.minimumTimeout),
				typematic.addListener(this.upArrowNode, this.textbox, {charOrCode:keys.PAGE_UP,ctrlKey:false,altKey:false,shiftKey:false,metaKey:false}, this, "_typematicCallback", this.timeoutChangeRate, this.defaultTimeout, this.minimumTimeout),
				typematic.addListener(this.downArrowNode, this.textbox, {charOrCode:keys.PAGE_DOWN,ctrlKey:false,altKey:false,shiftKey:false,metaKey:false}, this, "_typematicCallback", this.timeoutChangeRate, this.defaultTimeout, this.minimumTimeout)
			);
		}
	});
});

},
'url:dijit/form/templates/Spinner.html':"<div class=\"dijit dijitReset dijitInline dijitLeft\"\n\tid=\"widget_${id}\" role=\"presentation\"\n\t><div class=\"dijitReset dijitButtonNode dijitSpinnerButtonContainer\"\n\t\t><input class=\"dijitReset dijitInputField dijitSpinnerButtonInner\" type=\"text\" tabIndex=\"-1\" readonly=\"readonly\" role=\"presentation\"\n\t\t/><div class=\"dijitReset dijitLeft dijitButtonNode dijitArrowButton dijitUpArrowButton\"\n\t\t\tdata-dojo-attach-point=\"upArrowNode\"\n\t\t\t><div class=\"dijitArrowButtonInner\"\n\t\t\t\t><input class=\"dijitReset dijitInputField\" value=\"&#9650; \" type=\"text\" tabIndex=\"-1\" readonly=\"readonly\" role=\"presentation\"\n\t\t\t\t\t${_buttonInputDisabled}\n\t\t\t/></div\n\t\t></div\n\t\t><div class=\"dijitReset dijitLeft dijitButtonNode dijitArrowButton dijitDownArrowButton\"\n\t\t\tdata-dojo-attach-point=\"downArrowNode\"\n\t\t\t><div class=\"dijitArrowButtonInner\"\n\t\t\t\t><input class=\"dijitReset dijitInputField\" value=\"&#9660; \" type=\"text\" tabIndex=\"-1\" readonly=\"readonly\" role=\"presentation\"\n\t\t\t\t\t${_buttonInputDisabled}\n\t\t\t/></div\n\t\t></div\n\t></div\n\t><div class='dijitReset dijitValidationContainer'\n\t\t><input class=\"dijitReset dijitInputField dijitValidationIcon dijitValidationInner\" value=\"&#935; \" type=\"text\" tabIndex=\"-1\" readonly=\"readonly\" role=\"presentation\"\n\t/></div\n\t><div class=\"dijitReset dijitInputField dijitInputContainer\"\n\t\t><input class='dijitReset dijitInputInner' data-dojo-attach-point=\"textbox,focusNode\" type=\"${type}\" data-dojo-attach-event=\"onkeypress:_onKeyPress\"\n\t\t\trole=\"spinbutton\" autocomplete=\"off\" ${!nameAttrSetting}\n\t/></div\n></div>\n",
'dijit/form/NumberTextBox':function(){
define("dijit/form/NumberTextBox", [
	"dojo/_base/declare", // declare
	"dojo/_base/lang", // lang.hitch lang.mixin
	"dojo/number", // number._realNumberRegexp number.format number.parse number.regexp
	"./RangeBoundTextBox"
], function(declare, lang, number, RangeBoundTextBox){

	// module:
	//		dijit/form/NumberTextBox


	var NumberTextBoxMixin = declare("dijit.form.NumberTextBoxMixin", null, {
		// summary:
		//		A mixin for all number textboxes
		// tags:
		//		protected

		// Override ValidationTextBox.pattern.... we use a reg-ex generating function rather
		// than a straight regexp to deal with locale (plus formatting options too?)
		pattern: number.regexp,

		/*=====
		// constraints: NumberTextBox.__Constraints
		//		Despite the name, this parameter specifies both constraints on the input
		//		(including minimum/maximum allowed values) as well as
		//		formatting options like places (the number of digits to display after
		//		the decimal point).
		constraints: {},
		======*/

		// value: Number
		//		The value of this NumberTextBox as a Javascript Number (i.e., not a String).
		//		If the displayed value is blank, the value is NaN, and if the user types in
		//		an gibberish value (like "hello world"), the value is undefined
		//		(i.e. get('value') returns undefined).
		//
		//		Symmetrically, set('value', NaN) will clear the displayed value,
		//		whereas set('value', undefined) will have no effect.
		value: NaN,

		// editOptions: [protected] Object
		//		Properties to mix into constraints when the value is being edited.
		//		This is here because we edit the number in the format "12345", which is
		//		different than the display value (ex: "12,345")
		editOptions: { pattern: '#.######' },

		/*=====
		_formatter: function(value, options){
			// summary:
			//		_formatter() is called by format().  It's the base routine for formatting a number,
			//		as a string, for example converting 12345 into "12,345".
			// value: Number
			//		The number to be converted into a string.
			// options: number.__FormatOptions?
			//		Formatting options
			// tags:
			//		protected extension

			return "12345";		// String
		},
		 =====*/
		_formatter: number.format,

		postMixInProperties: function(){
			this.inherited(arguments);
			this._set("type", "text"); // in case type="number" was specified which messes up parse/format
		},

		_setConstraintsAttr: function(/*Object*/ constraints){
			var places = typeof constraints.places == "number"? constraints.places : 0;
			if(places){ places++; } // decimal rounding errors take away another digit of precision
			if(typeof constraints.max != "number"){
				constraints.max = 9 * Math.pow(10, 15-places);
			}
			if(typeof constraints.min != "number"){
				constraints.min = -9 * Math.pow(10, 15-places);
			}
			this.inherited(arguments, [ constraints ]);
			if(this.focusNode && this.focusNode.value && !isNaN(this.value)){
				this.set('value', this.value);
			}
		},

		_onFocus: function(){
			if(this.disabled){ return; }
			var val = this.get('value');
			if(typeof val == "number" && !isNaN(val)){
				var formattedValue = this.format(val, this.constraints);
				if(formattedValue !== undefined){
					this.textbox.value = formattedValue;
				}
			}
			this.inherited(arguments);
		},

		format: function(/*Number*/ value, /*number.__FormatOptions*/ constraints){
			// summary:
			//		Formats the value as a Number, according to constraints.
			// tags:
			//		protected

			var formattedValue = String(value);
			if(typeof value != "number"){ return formattedValue; }
			if(isNaN(value)){ return ""; }
			// check for exponential notation that dojo/number.format() chokes on
			if(!("rangeCheck" in this && this.rangeCheck(value, constraints)) && constraints.exponent !== false && /\de[-+]?\d/i.test(formattedValue)){
				return formattedValue;
			}
			if(this.editOptions && this.focused){
				constraints = lang.mixin({}, constraints, this.editOptions);
			}
			return this._formatter(value, constraints);
		},

		/*=====
		_parser: function(value, constraints){
			// summary:
			//		Parses the string value as a Number, according to constraints.
			// value: String
			//		String representing a number
			// constraints: number.__ParseOptions
			//		Formatting options
			// tags:
			//		protected

			return 123.45;		// Number
		},
		=====*/
		_parser: number.parse,

		parse: function(/*String*/ value, /*number.__FormatOptions*/ constraints){
			// summary:
			//		Replaceable function to convert a formatted string to a number value
			// tags:
			//		protected extension

			var v = this._parser(value, lang.mixin({}, constraints, (this.editOptions && this.focused) ? this.editOptions : {}));
			if(this.editOptions && this.focused && isNaN(v)){
				v = this._parser(value, constraints); // parse w/o editOptions: not technically needed but is nice for the user
			}
			return v;
		},

		_getDisplayedValueAttr: function(){
			var v = this.inherited(arguments);
			return isNaN(v) ? this.textbox.value : v;
		},

		filter: function(/*Number*/ value){
			// summary:
			//		This is called with both the display value (string), and the actual value (a number).
			//		When called with the actual value it does corrections so that '' etc. are represented as NaN.
			//		Otherwise it dispatches to the superclass's filter() method.
			//
			//		See `dijit/form/TextBox.filter()` for more details.
			return (value == null /* or undefined */ || value === '') ? NaN : this.inherited(arguments); // set('value', null||''||undefined) should fire onChange(NaN)
		},

		serialize: function(/*Number*/ value, /*Object?*/ options){
			// summary:
			//		Convert value (a Number) into a canonical string (ie, how the number literal is written in javascript/java/C/etc.)
			// tags:
			//		protected
			return (typeof value != "number" || isNaN(value)) ? '' : this.inherited(arguments);
		},

		_setBlurValue: function(){
			var val = lang.hitch(lang.mixin({}, this, { focused: true }), "get")('value'); // parse with editOptions
			this._setValueAttr(val, true);
		},

		_setValueAttr: function(/*Number*/ value, /*Boolean?*/ priorityChange, /*String?*/ formattedValue){
			// summary:
			//		Hook so set('value', ...) works.
			if(value !== undefined && formattedValue === undefined){
				formattedValue = String(value);
				if(typeof value == "number"){
					if(isNaN(value)){ formattedValue = '' }
					// check for exponential notation that number.format chokes on
					else if(("rangeCheck" in this && this.rangeCheck(value, this.constraints)) || this.constraints.exponent === false || !/\de[-+]?\d/i.test(formattedValue)){
						formattedValue = undefined; // lets format compute a real string value
					}
				}else if(!value){ // 0 processed in if branch above, ''|null|undefined flows through here
					formattedValue = '';
					value = NaN;
				}else{ // non-numeric values
					value = undefined;
				}
			}
			this.inherited(arguments, [value, priorityChange, formattedValue]);
		},

		_getValueAttr: function(){
			// summary:
			//		Hook so get('value') works.
			//		Returns Number, NaN for '', or undefined for unparseable text
			var v = this.inherited(arguments); // returns Number for all values accepted by parse() or NaN for all other displayed values

			// If the displayed value of the textbox is gibberish (ex: "hello world"), this.inherited() above
			// returns NaN; this if() branch converts the return value to undefined.
			// Returning undefined prevents user text from being overwritten when doing _setValueAttr(_getValueAttr()).
			// A blank displayed value is still returned as NaN.
			if(isNaN(v) && this.textbox.value !== ''){
				if(this.constraints.exponent !== false && /\de[-+]?\d/i.test(this.textbox.value) && (new RegExp("^"+number._realNumberRegexp(lang.mixin({}, this.constraints))+"$").test(this.textbox.value))){	// check for exponential notation that parse() rejected (erroneously?)
					var n = Number(this.textbox.value);
					return isNaN(n) ? undefined : n; // return exponential Number or undefined for random text (may not be possible to do with the above RegExp check)
				}else{
					return undefined; // gibberish
				}
			}else{
				return v; // Number or NaN for ''
			}
		},

		isValid: function(/*Boolean*/ isFocused){
			// Overrides dijit/form/RangeBoundTextBox.isValid() to check that the editing-mode value is valid since
			// it may not be formatted according to the regExp validation rules
			if(!this.focused || this._isEmpty(this.textbox.value)){
				return this.inherited(arguments);
			}else{
				var v = this.get('value');
				if(!isNaN(v) && this.rangeCheck(v, this.constraints)){
					if(this.constraints.exponent !== false && /\de[-+]?\d/i.test(this.textbox.value)){ // exponential, parse doesn't like it
						return true; // valid exponential number in range
					}else{
						return this.inherited(arguments);
					}
				}else{
					return false;
				}
			}
		}
	});

	var NumberTextBox = declare("dijit.form.NumberTextBox", [RangeBoundTextBox, NumberTextBoxMixin], {
		// summary:
		//		A TextBox for entering numbers, with formatting and range checking
		// description:
		//		NumberTextBox is a textbox for entering and displaying numbers, supporting
		//		the following main features:
		//
		//		1. Enforce minimum/maximum allowed values (as well as enforcing that the user types
		//			a number rather than a random string)
		//		2. NLS support (altering roles of comma and dot as "thousands-separator" and "decimal-point"
		//			depending on locale).
		//		3. Separate modes for editing the value and displaying it, specifically that
		//			the thousands separator character (typically comma) disappears when editing
		//			but reappears after the field is blurred.
		//		4. Formatting and constraints regarding the number of places (digits after the decimal point)
		//			allowed on input, and number of places displayed when blurred (see `constraints` parameter).

		baseClass: "dijitTextBox dijitNumberTextBox"
	});

	NumberTextBox.Mixin = NumberTextBoxMixin;	// for monkey patching

	/*=====
	 NumberTextBox.__Constraints = declare([RangeBoundTextBox.__Constraints, number.__FormatOptions, number.__ParseOptions], {
		 // summary:
		 //		Specifies both the rules on valid/invalid values (minimum, maximum,
		 //		number of required decimal places), and also formatting options for
		 //		displaying the value when the field is not focused.
		 // example:
		 //		Minimum/maximum:
		 //		To specify a field between 0 and 120:
		 //	|		{min:0,max:120}
		 //		To specify a field that must be an integer:
		 //	|		{fractional:false}
		 //		To specify a field where 0 to 3 decimal places are allowed on input:
		 //	|		{places:'0,3'}
	 });
	 =====*/

	return NumberTextBox;
});

},
'dijit/form/CurrencyTextBox':function(){
define("dijit/form/CurrencyTextBox", [
	"dojo/currency", // currency._mixInDefaults currency.format currency.parse currency.regexp
	"dojo/_base/declare", // declare
	"dojo/_base/lang", // lang.hitch
	"./NumberTextBox"
], function(currency, declare, lang, NumberTextBox){

	// module:
	//		dijit/form/CurrencyTextBox

	/*=====
	var __Constraints = declare([NumberTextBox.__Constraints, currency.__FormatOptions, currency.__ParseOptions], {
		// summary:
		//		Specifies both the rules on valid/invalid values (minimum, maximum,
		//		number of required decimal places), and also formatting options for
		//		displaying the value when the field is not focused (currency symbol,
		//		etc.)
		// description:
		//		Follows the pattern of `dijit/form/NumberTextBox.__Constraints`.
		//		In general developers won't need to set this parameter
		// example:
		//		To ensure that the user types in the cents (for example, 1.00 instead of just 1):
		//	|		{fractional:true}
	});
	=====*/

	return declare("dijit.form.CurrencyTextBox", NumberTextBox, {
		// summary:
		//		A validating currency textbox
		// description:
		//		CurrencyTextBox is similar to `dijit/form/NumberTextBox` but has a few
		//		extra features related to currency:
		//
		//		1. After specifying the currency type (american dollars, euros, etc.) it automatically
		//			sets parse/format options such as how many decimal places to show.
		//		2. The currency mark (dollar sign, euro mark, etc.) is displayed when the field is blurred
		//			but erased during editing, so that the user can just enter a plain number.

		// currency: [const] String
		//		the [ISO4217](http://en.wikipedia.org/wiki/ISO_4217) currency code, a three letter sequence like "USD"
		currency: "",

		/*=====
		// constraints: __Constraints
		//		Despite the name, this parameter specifies both constraints on the input
		//		(including minimum/maximum allowed values) as well as
		//		formatting options.
		constraints: {},
		======*/

		baseClass: "dijitTextBox dijitCurrencyTextBox",

		// Override pattern ValidationTextBox.pattern.... we use a reg-ex generating function rather
		// than a straight regexp to deal with locale (plus formatting options too?)
		pattern: function(constraints){
			// if focused, accept either currency data or NumberTextBox format
			return '(' + (this.focused ? this.inherited(arguments, [ lang.mixin({}, constraints, this.editOptions) ]) + '|' : '')
				+ currency.regexp(constraints) + ')';
		},

		// Override NumberTextBox._formatter to deal with currencies, ex: converts "123.45" to "$123.45"
		_formatter: currency.format,

		_parser: currency.parse,

		parse: function(/*String*/ value, /*Object*/ constraints){
			// summary:
			//		Parses string value as a Currency, according to the constraints object
			// tags:
			//		protected extension
			var v = this.inherited(arguments);
			if(isNaN(v) && /\d+/.test(value)){ // currency parse failed, but it could be because they are using NumberTextBox format so try its parse
				v = lang.hitch(lang.mixin({}, this, { _parser: NumberTextBox.prototype._parser }), "inherited")(arguments);
			}
			return v;
		},

		_setConstraintsAttr: function(/*Object*/ constraints){
			if(!constraints.currency && this.currency){
				constraints.currency = this.currency;
			}
			this.inherited(arguments, [ currency._mixInDefaults(lang.mixin(constraints, { exponent: false })) ]); // get places
		}
	});
});

},
'dojo/currency':function(){
define("dojo/currency", [
	"./_base/array",
	"./_base/lang",
	/*===== "./_base/declare", =====*/
	"./number",
	"./i18n", "./i18n!./cldr/nls/currency",
	"./cldr/monetary"
], function(darray, lang, /*===== declare, =====*/ dnumber, i18n, nlsCurrency, cldrMonetary){

// module:
//		dojo/currency

var currency = {
	// summary:
	//		localized formatting and parsing routines for currencies
	// description:
	//		extends dojo.number to provide culturally-appropriate formatting of values
	//		in various world currencies, including use of a currency symbol.  The currencies are specified
	//		by a three-letter international symbol in all uppercase, and support for the currencies is
	//		provided by the data in `dojo.cldr`.  The scripts generating dojo.cldr specify which
	//		currency support is included.  A fixed number of decimal places is determined based
	//		on the currency type and is not determined by the 'pattern' argument.  The fractional
	//		portion is optional, by default, and variable length decimals are not supported.
};
lang.setObject("dojo.currency", currency);

currency._mixInDefaults = function(options){
	options = options || {};
	options.type = "currency";

	// Get locale-dependent currency data, like the symbol
	var bundle = i18n.getLocalization("dojo.cldr", "currency", options.locale) || {};

	// Mixin locale-independent currency data, like # of places
	var iso = options.currency;
	var data = cldrMonetary.getData(iso);

	darray.forEach(["displayName","symbol","group","decimal"], function(prop){
		data[prop] = bundle[iso+"_"+prop];
	});

	data.fractional = [true, false];

	// Mixin with provided options
	return lang.mixin(data, options);
};

/*=====
currency.__FormatOptions = declare([dnumber.__FormatOptions], {
	// type: String?
	//		Should not be set.  Value is assumed to be "currency".
	// symbol: String?
	//		localized currency symbol. The default will be looked up in table of supported currencies in `dojo.cldr`
	//		A [ISO4217](http://en.wikipedia.org/wiki/ISO_4217) currency code will be used if not found.
	// currency: String?
	//		an [ISO4217](http://en.wikipedia.org/wiki/ISO_4217) currency code, a three letter sequence like "USD".
	//		For use with dojo.currency only.
	// places: Number?
	//		number of decimal places to show.  Default is defined based on which currency is used.
	type: "",
	symbol: "",
	currency: "",
	places: ""
});
=====*/

currency.format = function(/*Number*/ value, /*__FormatOptions?*/ options){
	// summary:
	//		Format a Number as a currency, using locale-specific settings
	//
	// description:
	//		Create a string from a Number using a known, localized pattern.
	//		[Formatting patterns](http://www.unicode.org/reports/tr35/#Number_Elements)
	//		appropriate to the locale are chosen from the [CLDR](http://unicode.org/cldr)
	//		as well as the appropriate symbols and delimiters and number of decimal places.
	//
	// value:
	//		the number to be formatted.

	return dnumber.format(value, currency._mixInDefaults(options));
};

currency.regexp = function(/*dnumber.__RegexpOptions?*/ options){
	//
	// summary:
	//		Builds the regular needed to parse a currency value
	//
	// description:
	//		Returns regular expression with positive and negative match, group and decimal separators
	//		Note: the options.places default, the number of decimal places to accept, is defined by the currency type.
	return dnumber.regexp(currency._mixInDefaults(options)); // String
};

/*=====
var __ParseOptions = currency.__ParseOptions = declare(dnumber.__ParseOptions, {
	// type: String?
	//		Should not be set.  Value is assumed to be currency.
	// currency: String?
	//		an [ISO4217](http://en.wikipedia.org/wiki/ISO_4217) currency code, a three letter sequence like "USD".
	//		For use with dojo.currency only.
	// symbol: String?
	//		localized currency symbol. The default will be looked up in table of supported currencies in `dojo.cldr`
	//		A [ISO4217](http://en.wikipedia.org/wiki/ISO_4217) currency code will be used if not found.
	// places: Number?
	//		fixed number of decimal places to accept.  The default is determined based on which currency is used.
	// fractional: Boolean|Array?
	//		Whether to include the fractional portion, where the number of decimal places are implied by the currency
	//		or explicit 'places' parameter.  The value [true,false] makes the fractional portion optional.
	//		By default for currencies, it the fractional portion is optional.
});
=====*/

currency.parse = function(/*String*/ expression, /*__ParseOptions?*/ options){
	//
	// summary:
	//		Convert a properly formatted currency string to a primitive Number,
	//		using locale-specific settings.
	// description:
	//		Create a Number from a string using a known, localized pattern.
	//		[Formatting patterns](http://www.unicode.org/reports/tr35/#Number_Format_Patterns)
	//		are chosen appropriate to the locale, as well as the appropriate symbols and delimiters
	//		and number of decimal places.
	// expression:
	//		A string representation of a currency value

	return dnumber.parse(expression, currency._mixInDefaults(options));
};

return currency;
});

},
'dojo/cldr/monetary':function(){
define("dojo/cldr/monetary", ["../_base/kernel", "../_base/lang"], function(dojo, lang){

// module:
//		dojo/cldr/monetary

var monetary = {
	// summary:
	//		TODOC
};
lang.setObject("dojo.cldr.monetary", monetary);

monetary.getData = function(/*String*/ code){
	// summary:
	//		A mapping of currency code to currency-specific formatting information. Returns a unique object with properties: places, round.
	// code:
	//		an [ISO 4217](http://en.wikipedia.org/wiki/ISO_4217) currency code

	// from http://www.unicode.org/cldr/data/common/supplemental/supplementalData.xml:supplementalData/currencyData/fractions

	var placesData = {
		ADP:0,AFN:0,ALL:0,AMD:0,BHD:3,BIF:0,BYR:0,CLF:0,CLP:0,
		COP:0,CRC:0,DJF:0,ESP:0,GNF:0,GYD:0,HUF:0,IDR:0,IQD:0,
		IRR:3,ISK:0,ITL:0,JOD:3,JPY:0,KMF:0,KPW:0,KRW:0,KWD:3,
		LAK:0,LBP:0,LUF:0,LYD:3,MGA:0,MGF:0,MMK:0,MNT:0,MRO:0,
		MUR:0,OMR:3,PKR:0,PYG:0,RSD:0,RWF:0,SLL:0,SOS:0,STD:0,
		SYP:0,TMM:0,TND:3,TRL:0,TZS:0,UGX:0,UZS:0,VND:0,VUV:0,
		XAF:0,XOF:0,XPF:0,YER:0,ZMK:0,ZWD:0
	};

	var roundingData = {CHF:5};

	var places = placesData[code], round = roundingData[code];
	if(typeof places == "undefined"){ places = 2; }
	if(typeof round == "undefined"){ round = 0; }

	return {places: places, round: round}; // Object
};

return monetary;
});

},
'dijit/form/HorizontalSlider':function(){
require({cache:{
'url:dijit/form/templates/HorizontalSlider.html':"<table class=\"dijit dijitReset dijitSlider dijitSliderH\" cellspacing=\"0\" cellpadding=\"0\" border=\"0\" rules=\"none\" data-dojo-attach-event=\"onkeypress:_onKeyPress,onkeyup:_onKeyUp\"\n\trole=\"presentation\"\n\t><tr class=\"dijitReset\"\n\t\t><td class=\"dijitReset\" colspan=\"2\"></td\n\t\t><td data-dojo-attach-point=\"topDecoration\" class=\"dijitReset dijitSliderDecoration dijitSliderDecorationT dijitSliderDecorationH\"></td\n\t\t><td class=\"dijitReset\" colspan=\"2\"></td\n\t></tr\n\t><tr class=\"dijitReset\"\n\t\t><td class=\"dijitReset dijitSliderButtonContainer dijitSliderButtonContainerH\"\n\t\t\t><div class=\"dijitSliderDecrementIconH\" style=\"display:none\" data-dojo-attach-point=\"decrementButton\"><span class=\"dijitSliderButtonInner\">-</span></div\n\t\t></td\n\t\t><td class=\"dijitReset\"\n\t\t\t><div class=\"dijitSliderBar dijitSliderBumper dijitSliderBumperH dijitSliderLeftBumper\" data-dojo-attach-event=\"press:_onClkDecBumper\"></div\n\t\t></td\n\t\t><td class=\"dijitReset\"\n\t\t\t><input data-dojo-attach-point=\"valueNode\" type=\"hidden\" ${!nameAttrSetting}\n\t\t\t/><div class=\"dijitReset dijitSliderBarContainerH\" role=\"presentation\" data-dojo-attach-point=\"sliderBarContainer\"\n\t\t\t\t><div role=\"presentation\" data-dojo-attach-point=\"progressBar\" class=\"dijitSliderBar dijitSliderBarH dijitSliderProgressBar dijitSliderProgressBarH\" data-dojo-attach-event=\"press:_onBarClick\"\n\t\t\t\t\t><div class=\"dijitSliderMoveable dijitSliderMoveableH\"\n\t\t\t\t\t\t><div data-dojo-attach-point=\"sliderHandle,focusNode\" class=\"dijitSliderImageHandle dijitSliderImageHandleH\" data-dojo-attach-event=\"press:_onHandleClick\" role=\"slider\"></div\n\t\t\t\t\t></div\n\t\t\t\t></div\n\t\t\t\t><div role=\"presentation\" data-dojo-attach-point=\"remainingBar\" class=\"dijitSliderBar dijitSliderBarH dijitSliderRemainingBar dijitSliderRemainingBarH\" data-dojo-attach-event=\"press:_onBarClick\"></div\n\t\t\t></div\n\t\t></td\n\t\t><td class=\"dijitReset\"\n\t\t\t><div class=\"dijitSliderBar dijitSliderBumper dijitSliderBumperH dijitSliderRightBumper\" data-dojo-attach-event=\"press:_onClkIncBumper\"></div\n\t\t></td\n\t\t><td class=\"dijitReset dijitSliderButtonContainer dijitSliderButtonContainerH\"\n\t\t\t><div class=\"dijitSliderIncrementIconH\" style=\"display:none\" data-dojo-attach-point=\"incrementButton\"><span class=\"dijitSliderButtonInner\">+</span></div\n\t\t></td\n\t></tr\n\t><tr class=\"dijitReset\"\n\t\t><td class=\"dijitReset\" colspan=\"2\"></td\n\t\t><td data-dojo-attach-point=\"containerNode,bottomDecoration\" class=\"dijitReset dijitSliderDecoration dijitSliderDecorationB dijitSliderDecorationH\"></td\n\t\t><td class=\"dijitReset\" colspan=\"2\"></td\n\t></tr\n></table>\n"}});
define("dijit/form/HorizontalSlider", [
	"dojo/_base/array", // array.forEach
	"dojo/_base/declare", // declare
	"dojo/dnd/move",
	"dojo/_base/event", // event.stop
	"dojo/_base/fx", // fx.animateProperty
	"dojo/dom-geometry", // domGeometry.position
	"dojo/dom-style", // domStyle.getComputedStyle
	"dojo/keys", // keys.DOWN_ARROW keys.END keys.HOME keys.LEFT_ARROW keys.PAGE_DOWN keys.PAGE_UP keys.RIGHT_ARROW keys.UP_ARROW
	"dojo/_base/lang", // lang.hitch
	"dojo/sniff", // has("ie") has("mozilla")
	"dojo/dnd/Moveable", // Moveable
	"dojo/dnd/Mover", // Mover Mover.prototype.destroy.apply
	"dojo/query", // query
	"dojo/mouse", // mouse.wheel
	"../registry", // registry.findWidgets
	"../focus",		// focus.focus()
	"../typematic",
	"./Button",
	"./_FormValueWidget",
	"../_Container",
	"dojo/text!./templates/HorizontalSlider.html"
], function(array, declare, move, event, fx, domGeometry, domStyle, keys, lang, has, Moveable, Mover, query, mouse,
			registry, focus, typematic, Button, _FormValueWidget, _Container, template){

// module:
//		dijit/form/HorizontalSlider


var _SliderMover = declare("dijit.form._SliderMover", Mover, {
	onMouseMove: function(e){
		var widget = this.widget;
		var abspos = widget._abspos;
		if(!abspos){
			abspos = widget._abspos = domGeometry.position(widget.sliderBarContainer, true);
			widget._setPixelValue_ = lang.hitch(widget, "_setPixelValue");
			widget._isReversed_ = widget._isReversed();
		}
		var pixelValue = e[widget._mousePixelCoord] - abspos[widget._startingPixelCoord];
		widget._setPixelValue_(widget._isReversed_ ? (abspos[widget._pixelCount]-pixelValue) : pixelValue, abspos[widget._pixelCount], false);
	},

	destroy: function(e){
		Mover.prototype.destroy.apply(this, arguments);
		var widget = this.widget;
		widget._abspos = null;
		widget._setValueAttr(widget.value, true);
	}
});

var HorizontalSlider = declare("dijit.form.HorizontalSlider", [_FormValueWidget, _Container], {
	// summary:
	//		A form widget that allows one to select a value with a horizontally draggable handle

	templateString: template,

	// Overrides FormValueWidget.value to indicate numeric value
	value: 0,

	// showButtons: [const] Boolean
	//		Show increment/decrement buttons at the ends of the slider?
	showButtons: true,

	// minimum: [const] Integer
	//		The minimum value the slider can be set to.
	minimum: 0,

	// maximum: [const] Integer
	//		The maximum value the slider can be set to.
	maximum: 100,

	// discreteValues: Integer
	//		If specified, indicates that the slider handle has only 'discreteValues' possible positions,
	//		and that after dragging the handle, it will snap to the nearest possible position.
	//		Thus, the slider has only 'discreteValues' possible values.
	//
	//		For example, if minimum=10, maxiumum=30, and discreteValues=3, then the slider handle has
	//		three possible positions, representing values 10, 20, or 30.
	//
	//		If discreteValues is not specified or if it's value is higher than the number of pixels
	//		in the slider bar, then the slider handle can be moved freely, and the slider's value will be
	//		computed/reported based on pixel position (in this case it will likely be fractional,
	//		such as 123.456789).
	discreteValues: Infinity,

	// pageIncrement: Integer
	//		If discreteValues is also specified, this indicates the amount of clicks (ie, snap positions)
	//		that the slider handle is moved via pageup/pagedown keys.
	//		If discreteValues is not specified, it indicates the number of pixels.
	pageIncrement: 2,

	// clickSelect: Boolean
	//		If clicking the slider bar changes the value or not
	clickSelect: true,

	// slideDuration: Number
	//		The time in ms to take to animate the slider handle from 0% to 100%,
	//		when clicking the slider bar to make the handle move.
	slideDuration: registry.defaultDuration,

	// Map widget attributes to DOMNode attributes.
	_setIdAttr: "",		// Override _FormWidget which sends id to focusNode

	baseClass: "dijitSlider",

	// Apply CSS classes to up/down arrows and handle per mouse state
	cssStateNodes: {
		incrementButton: "dijitSliderIncrementButton",
		decrementButton: "dijitSliderDecrementButton",
		focusNode: "dijitSliderThumb"
	},

	_mousePixelCoord: "pageX",
	_pixelCount: "w",
	_startingPixelCoord: "x",
	_handleOffsetCoord: "left",
	_progressPixelSize: "width",

	_onKeyUp: function(/*Event*/ e){
		if(this.disabled || this.readOnly || e.altKey || e.ctrlKey || e.metaKey){ return; }
		this._setValueAttr(this.value, true);
	},

	_onKeyPress: function(/*Event*/ e){
		if(this.disabled || this.readOnly || e.altKey || e.ctrlKey || e.metaKey){ return; }
		switch(e.charOrCode){
			case keys.HOME:
				this._setValueAttr(this.minimum, false);
				break;
			case keys.END:
				this._setValueAttr(this.maximum, false);
				break;
			// this._descending === false: if ascending vertical (min on top)
			// (this._descending || this.isLeftToRight()): if left-to-right horizontal or descending vertical
			case ((this._descending || this.isLeftToRight()) ? keys.RIGHT_ARROW : keys.LEFT_ARROW):
			case (this._descending === false ? keys.DOWN_ARROW : keys.UP_ARROW):
			case (this._descending === false ? keys.PAGE_DOWN : keys.PAGE_UP):
				this.increment(e);
				break;
			case ((this._descending || this.isLeftToRight()) ? keys.LEFT_ARROW : keys.RIGHT_ARROW):
			case (this._descending === false ? keys.UP_ARROW : keys.DOWN_ARROW):
			case (this._descending === false ? keys.PAGE_UP : keys.PAGE_DOWN):
				this.decrement(e);
				break;
			default:
				return;
		}
		event.stop(e);
	},

	_onHandleClick: function(e){
		if(this.disabled || this.readOnly){ return; }
		if(!has("ie")){
			// make sure you get focus when dragging the handle
			// (but don't do on IE because it causes a flicker on mouse up (due to blur then focus)
			focus.focus(this.sliderHandle);
		}
		event.stop(e);
	},

	_isReversed: function(){
		// summary:
		//		Returns true if direction is from right to left
		// tags:
		//		protected extension
		return !this.isLeftToRight();
	},

	_onBarClick: function(e){
		if(this.disabled || this.readOnly || !this.clickSelect){ return; }
		focus.focus(this.sliderHandle);
		event.stop(e);
		var abspos = domGeometry.position(this.sliderBarContainer, true);
		var pixelValue = e[this._mousePixelCoord] - abspos[this._startingPixelCoord];
		this._setPixelValue(this._isReversed() ? (abspos[this._pixelCount] - pixelValue) : pixelValue, abspos[this._pixelCount], true);
		this._movable.onMouseDown(e);
	},

	_setPixelValue: function(/*Number*/ pixelValue, /*Number*/ maxPixels, /*Boolean?*/ priorityChange){
		if(this.disabled || this.readOnly){ return; }
		var count = this.discreteValues;
		if(count <= 1 || count == Infinity){ count = maxPixels; }
		count--;
		var pixelsPerValue = maxPixels / count;
		var wholeIncrements = Math.round(pixelValue / pixelsPerValue);
		this._setValueAttr(Math.max(Math.min((this.maximum-this.minimum)*wholeIncrements/count + this.minimum, this.maximum), this.minimum), priorityChange);
	},

	_setValueAttr: function(/*Number*/ value, /*Boolean?*/ priorityChange){
		// summary:
		//		Hook so set('value', value) works.
		this._set("value", value);
		this.valueNode.value = value;
		this.focusNode.setAttribute("aria-valuenow", value);
		this.inherited(arguments);
		var percent = (value - this.minimum) / (this.maximum - this.minimum);
		var progressBar = (this._descending === false) ? this.remainingBar : this.progressBar;
		var remainingBar = (this._descending === false) ? this.progressBar : this.remainingBar;
		if(this._inProgressAnim && this._inProgressAnim.status != "stopped"){
			this._inProgressAnim.stop(true);
		}
		if(priorityChange && this.slideDuration > 0 && progressBar.style[this._progressPixelSize]){
			// animate the slider
			var _this = this;
			var props = {};
			var start = parseFloat(progressBar.style[this._progressPixelSize]);
			var duration = this.slideDuration * (percent-start/100);
			if(duration == 0){ return; }
			if(duration < 0){ duration = 0 - duration; }
			props[this._progressPixelSize] = { start: start, end: percent*100, units:"%" };
			this._inProgressAnim = fx.animateProperty({ node: progressBar, duration: duration,
				onAnimate: function(v){
					remainingBar.style[_this._progressPixelSize] = (100 - parseFloat(v[_this._progressPixelSize])) + "%";
				},
				onEnd: function(){
					delete _this._inProgressAnim;
				},
				properties: props
			});
			this._inProgressAnim.play();
		}else{
			progressBar.style[this._progressPixelSize] = (percent*100) + "%";
			remainingBar.style[this._progressPixelSize] = ((1-percent)*100) + "%";
		}
	},

	_bumpValue: function(signedChange, /*Boolean?*/ priorityChange){
		if(this.disabled || this.readOnly){ return; }
		var s = domStyle.getComputedStyle(this.sliderBarContainer);
		var c = domGeometry.getContentBox(this.sliderBarContainer, s);
		var count = this.discreteValues;
		if(count <= 1 || count == Infinity){ count = c[this._pixelCount]; }
		count--;
		var value = (this.value - this.minimum) * count / (this.maximum - this.minimum) + signedChange;
		if(value < 0){ value = 0; }
		if(value > count){ value = count; }
		value = value * (this.maximum - this.minimum) / count + this.minimum;
		this._setValueAttr(value, priorityChange);
	},

	_onClkBumper: function(val){
		if(this.disabled || this.readOnly || !this.clickSelect){ return; }
		this._setValueAttr(val, true);
	},

	_onClkIncBumper: function(){
		this._onClkBumper(this._descending === false ? this.minimum : this.maximum);
	},

	_onClkDecBumper: function(){
		this._onClkBumper(this._descending === false ? this.maximum : this.minimum);
	},

	decrement: function(/*Event*/ e){
		// summary:
		//		Decrement slider
		// tags:
		//		private
		this._bumpValue(e.charOrCode == keys.PAGE_DOWN ? -this.pageIncrement : -1);
	},

	increment: function(/*Event*/ e){
		// summary:
		//		Increment slider
		// tags:
		//		private
		this._bumpValue(e.charOrCode == keys.PAGE_UP ? this.pageIncrement : 1);
	},

	_mouseWheeled: function(/*Event*/ evt){
		// summary:
		//		Event handler for mousewheel where supported
		event.stop(evt);
		this._bumpValue(evt.wheelDelta < 0 ? -1 : 1, true); // negative scroll acts like a decrement
	},

	startup: function(){
		if(this._started){ return; }

		array.forEach(this.getChildren(), function(child){
			if(this[child.container] != this.containerNode){
				this[child.container].appendChild(child.domNode);
			}
		}, this);

		this.inherited(arguments);
	},

	_typematicCallback: function(/*Number*/ count, /*Object*/ button, /*Event*/ e){
		if(count == -1){
			this._setValueAttr(this.value, true);
		}else{
			this[(button == (this._descending? this.incrementButton : this.decrementButton)) ? "decrement" : "increment"](e);
		}
	},

	buildRendering: function(){
		this.inherited(arguments);
		if(this.showButtons){
			this.incrementButton.style.display="";
			this.decrementButton.style.display="";
		}

		// find any associated label element and add to slider focusnode.
		var label = query('label[for="'+this.id+'"]');
		if(label.length){
			if(!label[0].id){ label[0].id = this.id + "_label"; }
			this.focusNode.setAttribute("aria-labelledby", label[0].id);
		}

		this.focusNode.setAttribute("aria-valuemin", this.minimum);
		this.focusNode.setAttribute("aria-valuemax", this.maximum);
	},

	postCreate: function(){
		this.inherited(arguments);

		if(this.showButtons){
			this.own(
				typematic.addMouseListener(this.decrementButton, this, "_typematicCallback", 25, 500),
				typematic.addMouseListener(this.incrementButton, this, "_typematicCallback", 25, 500)
			);
		}
		this.connect(this.domNode, mouse.wheel, "_mouseWheeled");

		// define a custom constructor for a SliderMover that points back to me
		var mover = declare(_SliderMover, {
			widget: this
		});
		this._movable = new Moveable(this.sliderHandle, {mover: mover});

		this._layoutHackIE7();
	},

	destroy: function(){
		this._movable.destroy();
		if(this._inProgressAnim && this._inProgressAnim.status != "stopped"){
			this._inProgressAnim.stop(true);
		}
		this.inherited(arguments);
	}
});

HorizontalSlider._Mover = _SliderMover;	// for monkey patching

return HorizontalSlider;
});

},
'url:dijit/form/templates/HorizontalSlider.html':"<table class=\"dijit dijitReset dijitSlider dijitSliderH\" cellspacing=\"0\" cellpadding=\"0\" border=\"0\" rules=\"none\" data-dojo-attach-event=\"onkeypress:_onKeyPress,onkeyup:_onKeyUp\"\n\trole=\"presentation\"\n\t><tr class=\"dijitReset\"\n\t\t><td class=\"dijitReset\" colspan=\"2\"></td\n\t\t><td data-dojo-attach-point=\"topDecoration\" class=\"dijitReset dijitSliderDecoration dijitSliderDecorationT dijitSliderDecorationH\"></td\n\t\t><td class=\"dijitReset\" colspan=\"2\"></td\n\t></tr\n\t><tr class=\"dijitReset\"\n\t\t><td class=\"dijitReset dijitSliderButtonContainer dijitSliderButtonContainerH\"\n\t\t\t><div class=\"dijitSliderDecrementIconH\" style=\"display:none\" data-dojo-attach-point=\"decrementButton\"><span class=\"dijitSliderButtonInner\">-</span></div\n\t\t></td\n\t\t><td class=\"dijitReset\"\n\t\t\t><div class=\"dijitSliderBar dijitSliderBumper dijitSliderBumperH dijitSliderLeftBumper\" data-dojo-attach-event=\"press:_onClkDecBumper\"></div\n\t\t></td\n\t\t><td class=\"dijitReset\"\n\t\t\t><input data-dojo-attach-point=\"valueNode\" type=\"hidden\" ${!nameAttrSetting}\n\t\t\t/><div class=\"dijitReset dijitSliderBarContainerH\" role=\"presentation\" data-dojo-attach-point=\"sliderBarContainer\"\n\t\t\t\t><div role=\"presentation\" data-dojo-attach-point=\"progressBar\" class=\"dijitSliderBar dijitSliderBarH dijitSliderProgressBar dijitSliderProgressBarH\" data-dojo-attach-event=\"press:_onBarClick\"\n\t\t\t\t\t><div class=\"dijitSliderMoveable dijitSliderMoveableH\"\n\t\t\t\t\t\t><div data-dojo-attach-point=\"sliderHandle,focusNode\" class=\"dijitSliderImageHandle dijitSliderImageHandleH\" data-dojo-attach-event=\"press:_onHandleClick\" role=\"slider\"></div\n\t\t\t\t\t></div\n\t\t\t\t></div\n\t\t\t\t><div role=\"presentation\" data-dojo-attach-point=\"remainingBar\" class=\"dijitSliderBar dijitSliderBarH dijitSliderRemainingBar dijitSliderRemainingBarH\" data-dojo-attach-event=\"press:_onBarClick\"></div\n\t\t\t></div\n\t\t></td\n\t\t><td class=\"dijitReset\"\n\t\t\t><div class=\"dijitSliderBar dijitSliderBumper dijitSliderBumperH dijitSliderRightBumper\" data-dojo-attach-event=\"press:_onClkIncBumper\"></div\n\t\t></td\n\t\t><td class=\"dijitReset dijitSliderButtonContainer dijitSliderButtonContainerH\"\n\t\t\t><div class=\"dijitSliderIncrementIconH\" style=\"display:none\" data-dojo-attach-point=\"incrementButton\"><span class=\"dijitSliderButtonInner\">+</span></div\n\t\t></td\n\t></tr\n\t><tr class=\"dijitReset\"\n\t\t><td class=\"dijitReset\" colspan=\"2\"></td\n\t\t><td data-dojo-attach-point=\"containerNode,bottomDecoration\" class=\"dijitReset dijitSliderDecoration dijitSliderDecorationB dijitSliderDecorationH\"></td\n\t\t><td class=\"dijitReset\" colspan=\"2\"></td\n\t></tr\n></table>\n",
'dijit/Editor':function(){
define("dijit/Editor", [
	"require",
	"dojo/_base/array", // array.forEach
	"dojo/_base/declare", // declare
	"dojo/_base/Deferred", // Deferred
	"dojo/i18n", // i18n.getLocalization
	"dojo/dom-attr", // domAttr.set
	"dojo/dom-class", // domClass.add
	"dojo/dom-geometry",
	"dojo/dom-style", // domStyle.set, get
	"dojo/_base/event", // event.stop
	"dojo/keys", // keys.F1 keys.F15 keys.TAB
	"dojo/_base/lang", // lang.getObject lang.hitch
	"dojo/sniff", // has("ie") has("mac") has("webkit")
	"dojo/string", // string.substitute
	"dojo/topic", // topic.publish()
	"dojo/_base/window", // win.withGlobal
	"./_base/focus",	// dijit.getBookmark()
	"./_Container",
	"./Toolbar",
	"./ToolbarSeparator",
	"./layout/_LayoutWidget",
	"./form/ToggleButton",
	"./_editor/_Plugin",
	"./_editor/plugins/EnterKeyHandling",
	"./_editor/html",
	"./_editor/range",
	"./_editor/RichText",
	"./main",	// dijit._scopeName
	"dojo/i18n!./_editor/nls/commands"
], function(require, array, declare, Deferred, i18n, domAttr, domClass, domGeometry, domStyle,
			event, keys, lang, has, string, topic, win,
			focusBase, _Container, Toolbar, ToolbarSeparator, _LayoutWidget, ToggleButton,
			_Plugin, EnterKeyHandling, html, rangeapi, RichText, dijit){

	// module:
	//		dijit/Editor

	var Editor = declare("dijit.Editor", RichText, {
		// summary:
		//		A rich text Editing widget
		//
		// description:
		//		This widget provides basic WYSIWYG editing features, based on the browser's
		//		underlying rich text editing capability, accompanied by a toolbar (`dijit.Toolbar`).
		//		A plugin model is available to extend the editor's capabilities as well as the
		//		the options available in the toolbar.  Content generation may vary across
		//		browsers, and clipboard operations may have different results, to name
		//		a few limitations.  Note: this widget should not be used with the HTML
		//		&lt;TEXTAREA&gt; tag -- see dijit/_editor/RichText for details.

		// plugins: [const] Object[]
		//		A list of plugin names (as strings) or instances (as objects)
		//		for this widget.
		//
		//		When declared in markup, it might look like:
		//	|	plugins="['bold',{name:'dijit._editor.plugins.FontChoice', command:'fontName', generic:true}]"
		plugins: null,

		// extraPlugins: [const] Object[]
		//		A list of extra plugin names which will be appended to plugins array
		extraPlugins: null,

		constructor: function(/*===== params, srcNodeRef =====*/){
			// summary:
			//		Create the widget.
			// params: Object|null
			//		Initial settings for any of the attributes, except readonly attributes.
			// srcNodeRef: DOMNode
			//		The editor replaces the specified DOMNode.

			if(!lang.isArray(this.plugins)){
				this.plugins=["undo","redo","|","cut","copy","paste","|","bold","italic","underline","strikethrough","|",
				"insertOrderedList","insertUnorderedList","indent","outdent","|","justifyLeft","justifyRight","justifyCenter","justifyFull",
				EnterKeyHandling /*, "createLink"*/];
			}

			this._plugins=[];
			this._editInterval = this.editActionInterval * 1000;

			//IE will always lose focus when other element gets focus, while for FF and safari,
			//when no iframe is used, focus will be lost whenever another element gets focus.
			//For IE, we can connect to onBeforeDeactivate, which will be called right before
			//the focus is lost, so we can obtain the selected range. For other browsers,
			//no equivalent of onBeforeDeactivate, so we need to do two things to make sure
			//selection is properly saved before focus is lost: 1) when user clicks another
			//element in the page, in which case we listen to mousedown on the entire page and
			//see whether user clicks out of a focus editor, if so, save selection (focus will
			//only lost after onmousedown event is fired, so we can obtain correct caret pos.)
			//2) when user tabs away from the editor, which is handled in onKeyDown below.
			if(has("ie")){
				this.events.push("onBeforeDeactivate");
				this.events.push("onBeforeActivate");
			}
		},

		postMixInProperties: function(){
			// summary:
			//	Extension to make sure a deferred is in place before certain functions
			//	execute, like making sure all the plugins are properly inserted.

			// Set up a deferred so that the value isn't applied to the editor
			// until all the plugins load, needed to avoid timing condition
			// reported in #10537.
			this.setValueDeferred = new Deferred();
			this.inherited(arguments);
		},

		postCreate: function(){
			//for custom undo/redo, if enabled.
			this._steps=this._steps.slice(0);
			this._undoedSteps=this._undoedSteps.slice(0);

			if(lang.isArray(this.extraPlugins)){
				this.plugins=this.plugins.concat(this.extraPlugins);
			}

			this.inherited(arguments);

			this.commands = i18n.getLocalization("dijit._editor", "commands", this.lang);

			if(!this.toolbar){
				// if we haven't been assigned a toolbar, create one
				this.toolbar = new Toolbar({
					ownerDocument: this.ownerDocument,
					dir: this.dir,
					lang: this.lang
				});
				this.header.appendChild(this.toolbar.domNode);
			}

			array.forEach(this.plugins, this.addPlugin, this);

			// Okay, denote the value can now be set.
			this.setValueDeferred.resolve(true);

			domClass.add(this.iframe.parentNode, "dijitEditorIFrameContainer");
			domClass.add(this.iframe, "dijitEditorIFrame");
			domAttr.set(this.iframe, "allowTransparency", true);

			if(has("webkit")){
				// Disable selecting the entire editor by inadvertent double-clicks.
				// on buttons, title bar, etc.  Otherwise clicking too fast on
				// a button such as undo/redo selects the entire editor.
				domStyle.set(this.domNode, "KhtmlUserSelect", "none");
			}
			this.toolbar.startup();
			this.onNormalizedDisplayChanged(); //update toolbar button status
		},
		destroy: function(){
			array.forEach(this._plugins, function(p){
				if(p && p.destroy){
					p.destroy();
				}
			});
			this._plugins=[];
			this.toolbar.destroyRecursive();
			delete this.toolbar;
			this.inherited(arguments);
		},
		addPlugin: function(/*String||Object||Function*/ plugin, /*Integer?*/ index){
			// summary:
			//		takes a plugin name as a string or a plugin instance and
			//		adds it to the toolbar and associates it with this editor
			//		instance. The resulting plugin is added to the Editor's
			//		plugins array. If index is passed, it's placed in the plugins
			//		array at that index. No big magic, but a nice helper for
			//		passing in plugin names via markup.
			// plugin:
			//		String, args object, plugin instance, or plugin constructor
			// args:
			//		This object will be passed to the plugin constructor
			// index:
			//		Used when creating an instance from
			//		something already in this.plugins. Ensures that the new
			//		instance is assigned to this.plugins at that index.
			var args=lang.isString(plugin)?{name:plugin}:lang.isFunction(plugin)?{ctor:plugin}:plugin;
			if(!args.setEditor){
				var o={"args":args,"plugin":null,"editor":this};
				if(args.name){
					// search registry for a plugin factory matching args.name, if it's not there then
					// fallback to 1.0 API:
					// ask all loaded plugin modules to fill in o.plugin if they can (ie, if they implement args.name)
					// remove fallback for 2.0.
					if(_Plugin.registry[args.name]){
						o.plugin = _Plugin.registry[args.name](args);
					}else{
						topic.publish(dijit._scopeName + ".Editor.getPlugin", o);	// publish
					}
				}
				if(!o.plugin){
					try{
						// TODO: remove lang.getObject() call in 2.0
						var pc = args.ctor || lang.getObject(args.name) || require(args.name);
						if(pc){
							o.plugin = new pc(args);
						}
					}catch(e){
						throw new Error(this.id + ": cannot find plugin [" + args.name + "]");
					}
				}
				if(!o.plugin){
					throw new Error(this.id + ": cannot find plugin [" + args.name + "]");
				}
				plugin=o.plugin;
			}
			if(arguments.length > 1){
				this._plugins[index] = plugin;
			}else{
				this._plugins.push(plugin);
			}
			plugin.setEditor(this);
			if(lang.isFunction(plugin.setToolbar)){
				plugin.setToolbar(this.toolbar);
			}
		},

		//the following 2 functions are required to make the editor play nice under a layout widget, see #4070

		resize: function(size){
			// summary:
			//		Resize the editor to the specified size, see `dijit/layout/_LayoutWidget.resize()`
			if(size){
				// we've been given a height/width for the entire editor (toolbar + contents), calls layout()
				// to split the allocated size between the toolbar and the contents
				_LayoutWidget.prototype.resize.apply(this, arguments);
			}
			/*
			else{
				// do nothing, the editor is already laid out correctly.   The user has probably specified
				// the height parameter, which was used to set a size on the iframe
			}
			*/
		},
		layout: function(){
			// summary:
			//		Called from `dijit/layout/_LayoutWidget.resize()`.  This shouldn't be called directly
			// tags:
			//		protected

			// Converts the iframe (or rather the <div> surrounding it) to take all the available space
			// except what's needed for the header (toolbars) and footer (breadcrumbs, etc).
			// A class was added to the iframe container and some themes style it, so we have to
			// calc off the added margins and padding too. See tracker: #10662
			var areaHeight = (this._contentBox.h -
				(this.getHeaderHeight() + this.getFooterHeight() +
				 domGeometry.getPadBorderExtents(this.iframe.parentNode).h +
				 domGeometry.getMarginExtents(this.iframe.parentNode).h));
			this.editingArea.style.height = areaHeight + "px";
			if(this.iframe){
				this.iframe.style.height="100%";
			}
			this._layoutMode = true;
		},

		_onIEMouseDown: function(/*Event*/ e){
			// summary:
			//		IE only to prevent 2 clicks to focus
			// tags:
			//		private
			var outsideClientArea;
			// IE 8's componentFromPoint is broken, which is a shame since it
			// was smaller code, but oh well.  We have to do this brute force
			// to detect if the click was scroller or not.
			var b = this.document.body;
			var clientWidth = b.clientWidth;
			var clientHeight = b.clientHeight;
			var clientLeft = b.clientLeft;
			var offsetWidth = b.offsetWidth;
			var offsetHeight = b.offsetHeight;
			var offsetLeft = b.offsetLeft;

			//Check for vertical scroller click.
			if(/^rtl$/i.test(b.dir || "")){
				if(clientWidth < offsetWidth && e.x > clientWidth && e.x < offsetWidth){
					// Check the click was between width and offset width, if so, scroller
					outsideClientArea = true;
				}
			}else{
				// RTL mode, we have to go by the left offsets.
				if(e.x < clientLeft && e.x > offsetLeft){
					// Check the click was between width and offset width, if so, scroller
					outsideClientArea = true;
				}
			}
			if(!outsideClientArea){
				// Okay, might be horiz scroller, check that.
				if(clientHeight < offsetHeight && e.y > clientHeight && e.y < offsetHeight){
					// Horizontal scroller.
					outsideClientArea = true;
				}
			}
			if(!outsideClientArea){
				delete this._cursorToStart; // Remove the force to cursor to start position.
				delete this._savedSelection; // new mouse position overrides old selection
				if(e.target.tagName == "BODY"){
					this.defer("placeCursorAtEnd");
				}
				this.inherited(arguments);
			}
		},
		onBeforeActivate: function(){
			this._restoreSelection();
		},
		onBeforeDeactivate: function(e){
			// summary:
			//		Called on IE right before focus is lost.   Saves the selected range.
			// tags:
			//		private
			if(this.customUndo){
				this.endEditing(true);
			}
			//in IE, the selection will be lost when other elements get focus,
			//let's save focus before the editor is deactivated
			if(e.target.tagName != "BODY"){
				this._saveSelection();
			}
			//console.log('onBeforeDeactivate',this);
		},

		/* beginning of custom undo/redo support */

		// customUndo: Boolean
		//		Whether we shall use custom undo/redo support instead of the native
		//		browser support. By default, we now use custom undo.  It works better
		//		than native browser support and provides a consistent behavior across
		//		browsers with a minimal performance hit.  We already had the hit on
		//		the slowest browser, IE, anyway.
		customUndo: true,

		// editActionInterval: Integer
		//		When using customUndo, not every keystroke will be saved as a step.
		//		Instead typing (including delete) will be grouped together: after
		//		a user stops typing for editActionInterval seconds, a step will be
		//		saved; if a user resume typing within editActionInterval seconds,
		//		the timeout will be restarted. By default, editActionInterval is 3
		//		seconds.
		editActionInterval: 3,

		beginEditing: function(cmd){
			// summary:
			//		Called to note that the user has started typing alphanumeric characters, if it's not already noted.
			//		Deals with saving undo; see editActionInterval parameter.
			// tags:
			//		private
			if(!this._inEditing){
				this._inEditing=true;
				this._beginEditing(cmd);
			}
			if(this.editActionInterval>0){
				if(this._editTimer){
					this._editTimer.remove();
				}
				this._editTimer = this.defer("endEditing", this._editInterval);
			}
		},

		// TODO: declaring these in the prototype is meaningless, just create in the constructor/postCreate
		_steps:[],
		_undoedSteps:[],

		execCommand: function(cmd){
			// summary:
			//		Main handler for executing any commands to the editor, like paste, bold, etc.
			//		Called by plugins, but not meant to be called by end users.
			// tags:
			//		protected
			if(this.customUndo && (cmd == 'undo' || cmd == 'redo')){
				return this[cmd]();
			}else{
				if(this.customUndo){
					this.endEditing();
					this._beginEditing();
				}
				var r = this.inherited(arguments);
				if(this.customUndo){
					this._endEditing();
				}
				return r;
			}
		},

		_pasteImpl: function(){
			// summary:
			//		Over-ride of paste command control to make execCommand cleaner
			// tags:
			//		Protected
			return this._clipboardCommand("paste");
		},

		_cutImpl: function(){
			// summary:
			//		Over-ride of cut command control to make execCommand cleaner
			// tags:
			//		Protected
			return this._clipboardCommand("cut");
		},

		_copyImpl: function(){
			// summary:
			//		Over-ride of copy command control to make execCommand cleaner
			// tags:
			//		Protected
			return this._clipboardCommand("copy");
		},

		_clipboardCommand: function(cmd){
			// summary:
			//		Function to handle processing clipboard commands (or at least try to).
			// tags:
			//		Private
			var r;
			try{
				// Try to exec the superclass exec-command and see if it works.
				r = this.document.execCommand(cmd, false, null);
				if(has("webkit") && !r){ //see #4598: webkit does not guarantee clipboard support from js
					throw { code: 1011 }; // throw an object like Mozilla's error
				}
			}catch(e){
				//TODO: when else might we get an exception?  Do we need the Mozilla test below?
				if(e.code == 1011 /* Mozilla: service denied */ ||
					(e.code == 9 && has("opera") /* Opera not supported */)){
					// Warn user of platform limitation.  Cannot programmatically access clipboard. See ticket #4136
					var sub = string.substitute,
						accel = {cut:'X', copy:'C', paste:'V'};
					alert(sub(this.commands.systemShortcut,
						[this.commands[cmd], sub(this.commands[has("mac") ? 'appleKey' : 'ctrlKey'], [accel[cmd]])]));
				}
				r = false;
			}
			return r;
		},

		queryCommandEnabled: function(cmd){
			// summary:
			//		Returns true if specified editor command is enabled.
			//		Used by the plugins to know when to highlight/not highlight buttons.
			// tags:
			//		protected
			if(this.customUndo && (cmd == 'undo' || cmd == 'redo')){
				return cmd == 'undo' ? (this._steps.length > 1) : (this._undoedSteps.length > 0);
			}else{
				return this.inherited(arguments);
			}
		},
		_moveToBookmark: function(b){
			// summary:
			//		Selects the text specified in bookmark b
			// tags:
			//		private
			var bookmark = b.mark;
			var mark = b.mark;
			var col = b.isCollapsed;
			var r, sNode, eNode, sel;
			if(mark){
				if(has("ie") < 9){
					if(lang.isArray(mark)){
						//IE CONTROL, have to use the native bookmark.
						bookmark = [];
						array.forEach(mark,function(n){
							bookmark.push(rangeapi.getNode(n,this.editNode));
						},this);
						win.withGlobal(this.window,'moveToBookmark',focusBase,[{mark: bookmark, isCollapsed: col}]);
					}else{
						if(mark.startContainer && mark.endContainer){
							// Use the pseudo WC3 range API.  This works better for positions
							// than the IE native bookmark code.
							sel = rangeapi.getSelection(this.window);
							if(sel && sel.removeAllRanges){
								sel.removeAllRanges();
								r = rangeapi.create(this.window);
								sNode = rangeapi.getNode(mark.startContainer,this.editNode);
								eNode = rangeapi.getNode(mark.endContainer,this.editNode);
								if(sNode && eNode){
									// Okay, we believe we found the position, so add it into the selection
									// There are cases where it may not be found, particularly in undo/redo, when
									// IE changes the underlying DOM on us (wraps text in a <p> tag or similar.
									// So, in those cases, don't bother restoring selection.
									r.setStart(sNode,mark.startOffset);
									r.setEnd(eNode,mark.endOffset);
									sel.addRange(r);
								}
							}
						}
					}
				}else{//w3c range
					sel = rangeapi.getSelection(this.window);
					if(sel && sel.removeAllRanges){
						sel.removeAllRanges();
						r = rangeapi.create(this.window);
						sNode = rangeapi.getNode(mark.startContainer,this.editNode);
						eNode = rangeapi.getNode(mark.endContainer,this.editNode);
						if(sNode && eNode){
							// Okay, we believe we found the position, so add it into the selection
							// There are cases where it may not be found, particularly in undo/redo, when
							// formatting as been done and so on, so don't restore selection then.
							r.setStart(sNode,mark.startOffset);
							r.setEnd(eNode,mark.endOffset);
							sel.addRange(r);
						}
					}
				}
			}
		},
		_changeToStep: function(from, to){
			// summary:
			//		Reverts editor to "to" setting, from the undo stack.
			// tags:
			//		private
			this.setValue(to.text);
			var b=to.bookmark;
			if(!b){ return; }
			this._moveToBookmark(b);
		},
		undo: function(){
			// summary:
			//		Handler for editor undo (ex: ctrl-z) operation
			// tags:
			//		private
			var ret = false;
			if(!this._undoRedoActive){
				this._undoRedoActive = true;
				this.endEditing(true);
				var s=this._steps.pop();
				if(s && this._steps.length>0){
					this.focus();
					this._changeToStep(s,this._steps[this._steps.length-1]);
					this._undoedSteps.push(s);
					this.onDisplayChanged();
					delete this._undoRedoActive;
					ret = true;
				}
				delete this._undoRedoActive;
			}
			return ret;
		},
		redo: function(){
			// summary:
			//		Handler for editor redo (ex: ctrl-y) operation
			// tags:
			//		private
			var ret = false;
			if(!this._undoRedoActive){
				this._undoRedoActive = true;
				this.endEditing(true);
				var s=this._undoedSteps.pop();
				if(s && this._steps.length>0){
					this.focus();
					this._changeToStep(this._steps[this._steps.length-1],s);
					this._steps.push(s);
					this.onDisplayChanged();
					ret = true;
				}
				delete this._undoRedoActive;
			}
			return ret;
		},
		endEditing: function(ignore_caret){
			// summary:
			//		Called to note that the user has stopped typing alphanumeric characters, if it's not already noted.
			//		Deals with saving undo; see editActionInterval parameter.
			// tags:
			//		private
			if(this._editTimer){
				this._editTimer = this._editTimer.remove();
			}
			if(this._inEditing){
				this._endEditing(ignore_caret);
				this._inEditing=false;
			}
		},

		_getBookmark: function(){
			// summary:
			//		Get the currently selected text
			// tags:
			//		protected
			var b=win.withGlobal(this.window,focusBase.getBookmark);
			var tmp=[];
			if(b && b.mark){
				var mark = b.mark;
				if(has("ie") < 9){
					// Try to use the pseudo range API on IE for better accuracy.
					var sel = rangeapi.getSelection(this.window);
					if(!lang.isArray(mark)){
						if(sel){
							var range;
							if(sel.rangeCount){
								range = sel.getRangeAt(0);
							}
							if(range){
								b.mark = range.cloneRange();
							}else{
								b.mark = win.withGlobal(this.window,focusBase.getBookmark);
							}
						}
					}else{
						// Control ranges (img, table, etc), handle differently.
						array.forEach(b.mark,function(n){
							tmp.push(rangeapi.getIndex(n,this.editNode).o);
						},this);
						b.mark = tmp;
					}
				}
				try{
					if(b.mark && b.mark.startContainer){
						tmp=rangeapi.getIndex(b.mark.startContainer,this.editNode).o;
						b.mark={startContainer:tmp,
							startOffset:b.mark.startOffset,
							endContainer:b.mark.endContainer===b.mark.startContainer?tmp:rangeapi.getIndex(b.mark.endContainer,this.editNode).o,
							endOffset:b.mark.endOffset};
					}
				}catch(e){
					b.mark = null;
				}
			}
			return b;
		},
		_beginEditing: function(){
			// summary:
			//		Called when the user starts typing alphanumeric characters.
			//		Deals with saving undo; see editActionInterval parameter.
			// tags:
			//		private
			if(this._steps.length === 0){
				// You want to use the editor content without post filtering
				// to make sure selection restores right for the 'initial' state.
				// and undo is called.  So not using this.value, as it was 'processed'
				// and the line-up for selections may have been altered.
				this._steps.push({'text':html.getChildrenHtml(this.editNode),'bookmark':this._getBookmark()});
			}
		},
		_endEditing: function(){
			// summary:
			//		Called when the user stops typing alphanumeric characters.
			//		Deals with saving undo; see editActionInterval parameter.
			// tags:
			//		private
			
			// Avoid filtering to make sure selections restore.
			var v = html.getChildrenHtml(this.editNode);

			this._undoedSteps=[];//clear undoed steps
			this._steps.push({text: v, bookmark: this._getBookmark()});
		},
		onKeyDown: function(e){
			// summary:
			//		Handler for onkeydown event.
			// tags:
			//		private

			//We need to save selection if the user TAB away from this editor
			//no need to call _saveSelection for IE, as that will be taken care of in onBeforeDeactivate
			if(!has("ie") && !this.iframe && e.keyCode == keys.TAB && !this.tabIndent){
				this._saveSelection();
			}
			if(!this.customUndo){
				this.inherited(arguments);
				return;
			}
			var k = e.keyCode;
			if(e.ctrlKey && !e.altKey){//undo and redo only if the special right Alt + z/y are not pressed #5892
				if(k == 90 || k == 122){ //z
					event.stop(e);
					this.undo();
					return;
				}else if(k == 89 || k == 121){ //y
					event.stop(e);
					this.redo();
					return;
				}
			}
			this.inherited(arguments);

			switch(k){
					case keys.ENTER:
					case keys.BACKSPACE:
					case keys.DELETE:
						this.beginEditing();
						break;
					case 88: //x
					case 86: //v
						if(e.ctrlKey && !e.altKey && !e.metaKey){
							this.endEditing();//end current typing step if any
							if(e.keyCode == 88){
								this.beginEditing('cut');
							}else{
								this.beginEditing('paste');
							}
							//use timeout to trigger after the paste is complete
							this.defer("endEditing", 1);
							break;
						}
						//pass through
					default:
						if(!e.ctrlKey && !e.altKey && !e.metaKey && (e.keyCode<keys.F1 || e.keyCode>keys.F15)){
							this.beginEditing();
							break;
						}
						//pass through
					case keys.ALT:
						this.endEditing();
						break;
					case keys.UP_ARROW:
					case keys.DOWN_ARROW:
					case keys.LEFT_ARROW:
					case keys.RIGHT_ARROW:
					case keys.HOME:
					case keys.END:
					case keys.PAGE_UP:
					case keys.PAGE_DOWN:
						this.endEditing(true);
						break;
					//maybe ctrl+backspace/delete, so don't endEditing when ctrl is pressed
					case keys.CTRL:
					case keys.SHIFT:
					case keys.TAB:
						break;
				}
		},
		_onBlur: function(){
			// summary:
			//		Called from focus manager when focus has moved away from this editor
			// tags:
			//		protected

			//this._saveSelection();
			this.inherited(arguments);
			this.endEditing(true);
		},
		_saveSelection: function(){
			// summary:
			//		Save the currently selected text in _savedSelection attribute
			// tags:
			//		private
			try{
				this._savedSelection=this._getBookmark();
			}catch(e){ /* Squelch any errors that occur if selection save occurs due to being hidden simultaneously. */}
		},
		_restoreSelection: function(){
			// summary:
			//		Re-select the text specified in _savedSelection attribute;
			//		see _saveSelection().
			// tags:
			//		private
			if(this._savedSelection){
				// Clear off cursor to start, we're deliberately going to a selection.
				delete this._cursorToStart;
				// only restore the selection if the current range is collapsed
				// if not collapsed, then it means the editor does not lose
				// selection and there is no need to restore it
				if(win.withGlobal(this.window,'isCollapsed',focusBase)){
					this._moveToBookmark(this._savedSelection);
				}
				delete this._savedSelection;
			}
		},

		onClick: function(){
			// summary:
			//		Handler for when editor is clicked
			// tags:
			//		protected
			this.endEditing(true);
			this.inherited(arguments);
		},

		replaceValue: function(/*String*/ html){
			// summary:
			//		over-ride of replaceValue to support custom undo and stack maintenance.
			// tags:
			//		protected
			if(!this.customUndo){
				this.inherited(arguments);
			}else{
				if(this.isClosed){
					this.setValue(html);
				}else{
					this.beginEditing();
					if(!html){
						html = "&#160;";	// &nbsp;
					}
					this.setValue(html);
					this.endEditing();
				}
			}
		},

		_setDisabledAttr: function(/*Boolean*/ value){
			this.setValueDeferred.then(lang.hitch(this, function(){
				if((!this.disabled && value) || (!this._buttonEnabledPlugins && value)){
					// Disable editor: disable all enabled buttons and remember that list
					array.forEach(this._plugins, function(p){
						p.set("disabled", true);
					});
				}else if(this.disabled && !value){
					// Restore plugins to being active.
					array.forEach(this._plugins, function(p){
						p.set("disabled", false);
					});
				}
			}));
			this.inherited(arguments);
		},

		_setStateClass: function(){
			try{
				this.inherited(arguments);

				// Let theme set the editor's text color based on editor enabled/disabled state.
				// We need to jump through hoops because the main document (where the theme CSS is)
				// is separate from the iframe's document.
				if(this.document && this.document.body){
					domStyle.set(this.document.body, "color", domStyle.get(this.iframe, "color"));
				}
			}catch(e){ /* Squelch any errors caused by focus change if hidden during a state change */}
		}
	});

	// Register the "default plugins", ie, the built-in editor commands
	function simplePluginFactory(args){
		return new _Plugin({ command: args.name });
	}
	function togglePluginFactory(args){
		return new _Plugin({ buttonClass: ToggleButton, command: args.name });
	}
	lang.mixin(_Plugin.registry, {
		"undo": simplePluginFactory,
		"redo": simplePluginFactory,
		"cut": simplePluginFactory,
		"copy": simplePluginFactory,
		"paste": simplePluginFactory,
		"insertOrderedList": simplePluginFactory,
		"insertUnorderedList": simplePluginFactory,
		"indent": simplePluginFactory,
		"outdent": simplePluginFactory,
		"justifyCenter": simplePluginFactory,
		"justifyFull": simplePluginFactory,
		"justifyLeft": simplePluginFactory,
		"justifyRight": simplePluginFactory,
		"delete": simplePluginFactory,
		"selectAll": simplePluginFactory,
		"removeFormat": simplePluginFactory,
		"unlink": simplePluginFactory,
		"insertHorizontalRule": simplePluginFactory,

		"bold": togglePluginFactory,
		"italic": togglePluginFactory,
		"underline": togglePluginFactory,
		"strikethrough": togglePluginFactory,
		"subscript": togglePluginFactory,
		"superscript": togglePluginFactory,

		"|": function(){
			return new _Plugin({
				setEditor: function(editor){
					this.editor = editor;
					this.button = new ToolbarSeparator({ownerDocument: editor.ownerDocument});
				}
			});
		}
	});

	return Editor;
});

},
'dijit/_base/focus':function(){
define("dijit/_base/focus", [
	"dojo/_base/array", // array.forEach
	"dojo/dom", // dom.isDescendant
	"dojo/_base/lang", // lang.isArray
	"dojo/topic", // publish
	"dojo/_base/window", // win.doc win.doc.selection win.global win.global.getSelection win.withGlobal
	"../focus",
	"../main"	// for exporting symbols to dijit
], function(array, dom, lang, topic, win, focus, dijit){

	// module:
	//		dijit/_base/focus

	var exports = {
		// summary:
		//		Deprecated module to monitor currently focused node and stack of currently focused widgets.
		//		New code should access dijit/focus directly.

		// _curFocus: DomNode
		//		Currently focused item on screen
		_curFocus: null,

		// _prevFocus: DomNode
		//		Previously focused item on screen
		_prevFocus: null,

		isCollapsed: function(){
			// summary:
			//		Returns true if there is no text selected
			return dijit.getBookmark().isCollapsed;
		},

		getBookmark: function(){
			// summary:
			//		Retrieves a bookmark that can be used with moveToBookmark to return to the same range
			var bm, rg, tg, sel = win.doc.selection, cf = focus.curNode;

			if(win.global.getSelection){
				//W3C Range API for selections.
				sel = win.global.getSelection();
				if(sel){
					if(sel.isCollapsed){
						tg = cf? cf.tagName : "";
						if(tg){
							//Create a fake rangelike item to restore selections.
							tg = tg.toLowerCase();
							if(tg == "textarea" ||
									(tg == "input" && (!cf.type || cf.type.toLowerCase() == "text"))){
								sel = {
									start: cf.selectionStart,
									end: cf.selectionEnd,
									node: cf,
									pRange: true
								};
								return {isCollapsed: (sel.end <= sel.start), mark: sel}; //Object.
							}
						}
						bm = {isCollapsed:true};
						if(sel.rangeCount){
							bm.mark = sel.getRangeAt(0).cloneRange();
						}
					}else{
						rg = sel.getRangeAt(0);
						bm = {isCollapsed: false, mark: rg.cloneRange()};
					}
				}
			}else if(sel){
				// If the current focus was a input of some sort and no selection, don't bother saving
				// a native bookmark.  This is because it causes issues with dialog/page selection restore.
				// So, we need to create psuedo bookmarks to work with.
				tg = cf ? cf.tagName : "";
				tg = tg.toLowerCase();
				if(cf && tg && (tg == "button" || tg == "textarea" || tg == "input")){
					if(sel.type && sel.type.toLowerCase() == "none"){
						return {
							isCollapsed: true,
							mark: null
						}
					}else{
						rg = sel.createRange();
						return {
							isCollapsed: rg.text && rg.text.length?false:true,
							mark: {
								range: rg,
								pRange: true
							}
						};
					}
				}
				bm = {};

				//'IE' way for selections.
				try{
					// createRange() throws exception when dojo in iframe
					//and nothing selected, see #9632
					rg = sel.createRange();
					bm.isCollapsed = !(sel.type == 'Text' ? rg.htmlText.length : rg.length);
				}catch(e){
					bm.isCollapsed = true;
					return bm;
				}
				if(sel.type.toUpperCase() == 'CONTROL'){
					if(rg.length){
						bm.mark=[];
						var i=0,len=rg.length;
						while(i<len){
							bm.mark.push(rg.item(i++));
						}
					}else{
						bm.isCollapsed = true;
						bm.mark = null;
					}
				}else{
					bm.mark = rg.getBookmark();
				}
			}else{
				console.warn("No idea how to store the current selection for this browser!");
			}
			return bm; // Object
		},

		moveToBookmark: function(/*Object*/ bookmark){
			// summary:
			//		Moves current selection to a bookmark
			// bookmark:
			//		This should be a returned object from dijit.getBookmark()

			var _doc = win.doc,
				mark = bookmark.mark;
			if(mark){
				if(win.global.getSelection){
					//W3C Rangi API (FF, WebKit, Opera, etc)
					var sel = win.global.getSelection();
					if(sel && sel.removeAllRanges){
						if(mark.pRange){
							var n = mark.node;
							n.selectionStart = mark.start;
							n.selectionEnd = mark.end;
						}else{
							sel.removeAllRanges();
							sel.addRange(mark);
						}
					}else{
						console.warn("No idea how to restore selection for this browser!");
					}
				}else if(_doc.selection && mark){
					//'IE' way.
					var rg;
					if(mark.pRange){
						rg = mark.range;
					}else if(lang.isArray(mark)){
						rg = _doc.body.createControlRange();
						//rg.addElement does not have call/apply method, so can not call it directly
						//rg is not available in "range.addElement(item)", so can't use that either
						array.forEach(mark, function(n){
							rg.addElement(n);
						});
					}else{
						rg = _doc.body.createTextRange();
						rg.moveToBookmark(mark);
					}
					rg.select();
				}
			}
		},

		getFocus: function(/*Widget?*/ menu, /*Window?*/ openedForWindow){
			// summary:
			//		Called as getFocus(), this returns an Object showing the current focus
			//		and selected text.
			//
			//		Called as getFocus(widget), where widget is a (widget representing) a button
			//		that was just pressed, it returns where focus was before that button
			//		was pressed.   (Pressing the button may have either shifted focus to the button,
			//		or removed focus altogether.)   In this case the selected text is not returned,
			//		since it can't be accurately determined.
			//
			// menu: dijit/_WidgetBase|{domNode: DomNode} structure
			//		The button that was just pressed.  If focus has disappeared or moved
			//		to this button, returns the previous focus.  In this case the bookmark
			//		information is already lost, and null is returned.
			//
			// openedForWindow:
			//		iframe in which menu was opened
			//
			// returns:
			//		A handle to restore focus/selection, to be passed to `dijit.focus`
			var node = !focus.curNode || (menu && dom.isDescendant(focus.curNode, menu.domNode)) ? dijit._prevFocus : focus.curNode;
			return {
				node: node,
				bookmark: node && (node == focus.curNode) && win.withGlobal(openedForWindow || win.global, dijit.getBookmark),
				openedForWindow: openedForWindow
			}; // Object
		},

		// _activeStack: dijit/_WidgetBase[]
		//		List of currently active widgets (focused widget and it's ancestors)
		_activeStack: [],

		registerIframe: function(/*DomNode*/ iframe){
			// summary:
			//		Registers listeners on the specified iframe so that any click
			//		or focus event on that iframe (or anything in it) is reported
			//		as a focus/click event on the `<iframe>` itself.
			// description:
			//		Currently only used by editor.
			// returns:
			//		Handle to pass to unregisterIframe()
			return focus.registerIframe(iframe);
		},

		unregisterIframe: function(/*Object*/ handle){
			// summary:
			//		Unregisters listeners on the specified iframe created by registerIframe.
			//		After calling be sure to delete or null out the handle itself.
			// handle:
			//		Handle returned by registerIframe()

			handle && handle.remove();
		},

		registerWin: function(/*Window?*/targetWindow, /*DomNode?*/ effectiveNode){
			// summary:
			//		Registers listeners on the specified window (either the main
			//		window or an iframe's window) to detect when the user has clicked somewhere
			//		or focused somewhere.
			// description:
			//		Users should call registerIframe() instead of this method.
			// targetWindow:
			//		If specified this is the window associated with the iframe,
			//		i.e. iframe.contentWindow.
			// effectiveNode:
			//		If specified, report any focus events inside targetWindow as
			//		an event on effectiveNode, rather than on evt.target.
			// returns:
			//		Handle to pass to unregisterWin()

			return focus.registerWin(targetWindow, effectiveNode);
		},

		unregisterWin: function(/*Handle*/ handle){
			// summary:
			//		Unregisters listeners on the specified window (either the main
			//		window or an iframe's window) according to handle returned from registerWin().
			//		After calling be sure to delete or null out the handle itself.

			handle && handle.remove();
		}
	};

	// Override focus singleton's focus function so that dijit.focus()
	// has backwards compatible behavior of restoring selection (although
	// probably no one is using that).
	focus.focus = function(/*Object|DomNode */ handle){
		// summary:
		//		Sets the focused node and the selection according to argument.
		//		To set focus to an iframe's content, pass in the iframe itself.
		// handle:
		//		object returned by get(), or a DomNode

		if(!handle){ return; }

		var node = "node" in handle ? handle.node : handle,		// because handle is either DomNode or a composite object
			bookmark = handle.bookmark,
			openedForWindow = handle.openedForWindow,
			collapsed = bookmark ? bookmark.isCollapsed : false;

		// Set the focus
		// Note that for iframe's we need to use the <iframe> to follow the parentNode chain,
		// but we need to set focus to iframe.contentWindow
		if(node){
			var focusNode = (node.tagName.toLowerCase() == "iframe") ? node.contentWindow : node;
			if(focusNode && focusNode.focus){
				try{
					// Gecko throws sometimes if setting focus is impossible,
					// node not displayed or something like that
					focusNode.focus();
				}catch(e){/*quiet*/}
			}
			focus._onFocusNode(node);
		}

		// set the selection
		// do not need to restore if current selection is not empty
		// (use keyboard to select a menu item) or if previous selection was collapsed
		// as it may cause focus shift (Esp in IE).
		if(bookmark && win.withGlobal(openedForWindow || win.global, dijit.isCollapsed) && !collapsed){
			if(openedForWindow){
				openedForWindow.focus();
			}
			try{
				win.withGlobal(openedForWindow || win.global, dijit.moveToBookmark, null, [bookmark]);
			}catch(e2){
				/*squelch IE internal error, see http://trac.dojotoolkit.org/ticket/1984 */
			}
		}
	};

	// For back compatibility, monitor changes to focused node and active widget stack,
	// publishing events and copying changes from focus manager variables into dijit (top level) variables
	focus.watch("curNode", function(name, oldVal, newVal){
		dijit._curFocus = newVal;
		dijit._prevFocus = oldVal;
		if(newVal){
			topic.publish("focusNode", newVal);	// publish
		}
	});
	focus.watch("activeStack", function(name, oldVal, newVal){
		dijit._activeStack = newVal;
	});

	focus.on("widget-blur", function(widget, by){
		topic.publish("widgetBlur", widget, by);	// publish
	});
	focus.on("widget-focus", function(widget, by){
		topic.publish("widgetFocus", widget, by);	// publish
	});

	lang.mixin(dijit, exports);

	/*===== return exports; =====*/
	return dijit;	// for back compat :-(
});

},
'dijit/Toolbar':function(){
define("dijit/Toolbar", [
	"require",
	"dojo/_base/declare", // declare
	"dojo/has",
	"dojo/keys", // keys.LEFT_ARROW keys.RIGHT_ARROW
	"dojo/ready",
	"./_Widget",
	"./_KeyNavContainer",
	"./_TemplatedMixin"
], function(require, declare, has, keys, ready, _Widget, _KeyNavContainer, _TemplatedMixin){

	// module:
	//		dijit/Toolbar


	// Back compat w/1.6, remove for 2.0
	if(has("dijit-legacy-requires")){
		ready(0, function(){
			var requires = ["dijit/ToolbarSeparator"];
			require(requires);	// use indirection so modules not rolled into a build
		});
	}

	return declare("dijit.Toolbar", [_Widget, _TemplatedMixin, _KeyNavContainer], {
		// summary:
		//		A Toolbar widget, used to hold things like `dijit.Editor` buttons

		templateString:
			'<div class="dijit" role="toolbar" tabIndex="${tabIndex}" data-dojo-attach-point="containerNode">' +
			'</div>',

		baseClass: "dijitToolbar",

		postCreate: function(){
			this.inherited(arguments);

			this.connectKeyNavHandlers(
				this.isLeftToRight() ? [keys.LEFT_ARROW] : [keys.RIGHT_ARROW],
				this.isLeftToRight() ? [keys.RIGHT_ARROW] : [keys.LEFT_ARROW]
			);
		}
	});
});

},
'dijit/ToolbarSeparator':function(){
define("dijit/ToolbarSeparator", [
	"dojo/_base/declare", // declare
	"dojo/dom", // dom.setSelectable
	"./_Widget",
	"./_TemplatedMixin"
], function(declare, dom, _Widget, _TemplatedMixin){

	// module:
	//		dijit/ToolbarSeparator


	return declare("dijit.ToolbarSeparator", [_Widget, _TemplatedMixin], {
		// summary:
		//		A spacer between two `dijit.Toolbar` items

		templateString: '<div class="dijitToolbarSeparator dijitInline" role="presentation"></div>',

		buildRendering: function(){
			this.inherited(arguments);
			dom.setSelectable(this.domNode, false);
		},

		isFocusable: function(){
			// summary:
			//		This widget isn't focusable, so pass along that fact.
			// tags:
			//		protected
			return false;
		}
	});
});

},
'dijit/_editor/_Plugin':function(){
define("dijit/_editor/_Plugin", [
	"dojo/_base/connect", // connect.connect
	"dojo/_base/declare", // declare
	"dojo/_base/lang", // lang.mixin, lang.hitch
	"../form/Button"
], function(connect, declare, lang, Button){

// module:
//		dijit/_editor/_Plugin
// summary:
//		Base class for a "plugin" to the editor, which is usually
//		a single button on the Toolbar and some associated code


var _Plugin = declare("dijit._editor._Plugin", null, {
	// summary:
	//		Base class for a "plugin" to the editor, which is usually
	//		a single button on the Toolbar and some associated code

	constructor: function(args){
		// summary:
		//		Create the plugin.
		// args: Object?
		//		Initial settings for any of the attributes.

		this.params = args || {};
		lang.mixin(this, this.params);
		this._connects=[];
		this._attrPairNames = {};
	},

	// editor: [const] dijit.Editor
	//		Points to the parent editor
	editor: null,

	// iconClassPrefix: [const] String
	//		The CSS class name for the button node is formed from `iconClassPrefix` and `command`
	iconClassPrefix: "dijitEditorIcon",

	// button: dijit/_WidgetBase?
	//		Pointer to `dijit/form/Button` or other widget (ex: `dijit/form/FilteringSelect`)
	//		that is added to the toolbar to control this plugin.
	//		If not specified, will be created on initialization according to `buttonClass`
	button: null,

	// command: String
	//		String like "insertUnorderedList", "outdent", "justifyCenter", etc. that represents an editor command.
	//		Passed to editor.execCommand() if `useDefaultCommand` is true.
	command: "",

	// useDefaultCommand: Boolean
	//		If true, this plugin executes by calling Editor.execCommand() with the argument specified in `command`.
	useDefaultCommand: true,

	// buttonClass: Widget Class
	//		Class of widget (ex: dijit.form.Button or dijit/form/FilteringSelect)
	//		that is added to the toolbar to control this plugin.
	//		This is used to instantiate the button, unless `button` itself is specified directly.
	buttonClass: Button,

	// disabled: Boolean
	//		Flag to indicate if this plugin has been disabled and should do nothing
	//		helps control button state, among other things.  Set via the setter api.
	disabled: false,

	getLabel: function(/*String*/key){
		// summary:
		//		Returns the label to use for the button
		// tags:
		//		private
		return this.editor.commands[key];		// String
	},

	_initButton: function(){
		// summary:
		//		Initialize the button or other widget that will control this plugin.
		//		This code only works for plugins controlling built-in commands in the editor.
		// tags:
		//		protected extension
		if(this.command.length){
			var label = this.getLabel(this.command),
				editor = this.editor,
				className = this.iconClassPrefix+" "+this.iconClassPrefix + this.command.charAt(0).toUpperCase() + this.command.substr(1);
			if(!this.button){
				var props = lang.mixin({
					label: label,
					ownerDocument: editor.ownerDocument,
					dir: editor.dir,
					lang: editor.lang,
					showLabel: false,
					iconClass: className,
					dropDown: this.dropDown,
					tabIndex: "-1"
				}, this.params || {});
				this.button = new this.buttonClass(props);
			}
		}
		if(this.get("disabled") && this.button){
			this.button.set("disabled", this.get("disabled"));
		}
	},

	destroy: function(){
		// summary:
		//		Destroy this plugin

		var h;
		while(h = this._connects.pop()){ h.remove(); }
		if(this.dropDown){
			this.dropDown.destroyRecursive();
		}
	},

	connect: function(o, f, tf){
		// summary:
		//		Make a connect.connect() that is automatically disconnected when this plugin is destroyed.
		//		Similar to `dijit/_Widget.connect()`.
		// tags:
		//		protected
		this._connects.push(connect.connect(o, f, this, tf));
	},

	updateState: function(){
		// summary:
		//		Change state of the plugin to respond to events in the editor.
		// description:
		//		This is called on meaningful events in the editor, such as change of selection
		//		or caret position (but not simple typing of alphanumeric keys).   It gives the
		//		plugin a chance to update the CSS of its button.
		//
		//		For example, the "bold" plugin will highlight/unhighlight the bold button depending on whether the
		//		characters next to the caret are bold or not.
		//
		//		Only makes sense when `useDefaultCommand` is true, as it calls Editor.queryCommandEnabled(`command`).
		var e = this.editor,
			c = this.command,
			checked, enabled;
		if(!e || !e.isLoaded || !c.length){ return; }
		var disabled = this.get("disabled");
		if(this.button){
			try{
				enabled = !disabled && e.queryCommandEnabled(c);
				if(this.enabled !== enabled){
					this.enabled = enabled;
					this.button.set('disabled', !enabled);
				}
				if(enabled){
					if(typeof this.button.checked == 'boolean'){
						checked = e.queryCommandState(c);
						if(this.checked !== checked){
							this.checked = checked;
							this.button.set('checked', e.queryCommandState(c));
						}
					}
				}
			}catch(e){
				console.log(e); // FIXME: we shouldn't have debug statements in our code.  Log as an error?
			}
		}
	},

	setEditor: function(/*dijit/Editor*/ editor){
		// summary:
		//		Tell the plugin which Editor it is associated with.

		// TODO: refactor code to just pass editor to constructor.

		// FIXME: detach from previous editor!!
		this.editor = editor;

		// FIXME: prevent creating this if we don't need to (i.e., editor can't handle our command)
		this._initButton();

		// Processing for buttons that execute by calling editor.execCommand()
		if(this.button && this.useDefaultCommand){
			if(this.editor.queryCommandAvailable(this.command)){
				this.connect(this.button, "onClick",
					lang.hitch(this.editor, "execCommand", this.command, this.commandArg)
				);
			}else{
				// hide button because editor doesn't support command (due to browser limitations)
				this.button.domNode.style.display = "none";
			}
		}

		this.connect(this.editor, "onNormalizedDisplayChanged", "updateState");
	},

	setToolbar: function(/*dijit/Toolbar*/ toolbar){
		// summary:
		//		Tell the plugin to add it's controller widget (often a button)
		//		to the toolbar.  Does nothing if there is no controller widget.

		// TODO: refactor code to just pass toolbar to constructor.

		if(this.button){
			toolbar.addChild(this.button);
		}
		// console.debug("adding", this.button, "to:", toolbar);
	},

	set: function(/* attribute */ name, /* anything */ value){
		// summary:
		//		Set a property on a plugin
		// name:
		//		The property to set.
		// value:
		//		The value to set in the property.
		// description:
		//		Sets named properties on a plugin which may potentially be handled by a
		//		setter in the plugin.
		//		For example, if the plugin has a properties "foo"
		//		and "bar" and a method named "_setFooAttr", calling:
		//	|	plugin.set("foo", "Howdy!");
		//		would be equivalent to writing:
		//	|	plugin._setFooAttr("Howdy!");
		//		and:
		//	|	plugin.set("bar", 3);
		//		would be equivalent to writing:
		//	|	plugin.bar = 3;
		//
		//		set() may also be called with a hash of name/value pairs, ex:
		//	|	plugin.set({
		//	|		foo: "Howdy",
		//	|		bar: 3
		//	|	})
		//		This is equivalent to calling set(foo, "Howdy") and set(bar, 3)
		if(typeof name === "object"){
			for(var x in name){
				this.set(x, name[x]);
	}
			return this;
		}
		var names = this._getAttrNames(name);
		if(this[names.s]){
			// use the explicit setter
			var result = this[names.s].apply(this, Array.prototype.slice.call(arguments, 1));
		}else{
			this._set(name, value);
		}
		return result || this;
	},

	get: function(name){
		// summary:
		//		Get a property from a plugin.
		// name:
		//		The property to get.
		// description:
		//		Get a named property from a plugin. The property may
		//		potentially be retrieved via a getter method. If no getter is defined, this
		//		just retrieves the object's property.
		//		For example, if the plugin has a properties "foo"
		//		and "bar" and a method named "_getFooAttr", calling:
		//	|	plugin.get("foo");
		//		would be equivalent to writing:
		//	|	plugin._getFooAttr();
		//		and:
		//	|	plugin.get("bar");
		//		would be equivalent to writing:
		//	|	plugin.bar;
		var names = this._getAttrNames(name);
		return this[names.g] ? this[names.g]() : this[name];
	},

	_setDisabledAttr: function(disabled){
		// summary:
		//		Function to set the plugin state and call updateState to make sure the
		//		button is updated appropriately.
		this.disabled = disabled;
		this.updateState();
	},

	_getAttrNames: function(name){
		// summary:
		//		Helper function for get() and set().
		//		Caches attribute name values so we don't do the string ops every time.
		// tags:
		//		private

		var apn = this._attrPairNames;
		if(apn[name]){ return apn[name]; }
		var uc = name.charAt(0).toUpperCase() + name.substr(1);
		return (apn[name] = {
			s: "_set"+uc+"Attr",
			g: "_get"+uc+"Attr"
		});
	},

	_set: function(/*String*/ name, /*anything*/ value){
		// summary:
		//		Helper function to set new value for specified attribute
		this[name] = value;
	}
});

// Hash mapping plugin name to factory, used for registering plugins
_Plugin.registry = {};

return _Plugin;

});

},
'dijit/_editor/plugins/EnterKeyHandling':function(){
define("dijit/_editor/plugins/EnterKeyHandling", [
	"dojo/_base/declare", // declare
	"dojo/dom-construct", // domConstruct.destroy domConstruct.place
	"dojo/_base/event", // event.stop
	"dojo/keys", // keys.ENTER
	"dojo/_base/lang",
	"dojo/sniff", // has("ie") has("mozilla") has("webkit")
	"dojo/_base/window", // win.withGlobal
	"dojo/window", // winUtils.scrollIntoView
	"../_Plugin",
	"../RichText",
	"../range",
	"../../_base/focus"
], function(declare, domConstruct, event, keys, lang, has, win, winUtils, _Plugin, RichText, rangeapi, baseFocus){

// module:
//		dijit/_editor/plugins/EnterKeyHandling

return declare("dijit._editor.plugins.EnterKeyHandling", _Plugin, {
	// summary:
	//		This plugin tries to make all browsers behave consistently with regard to
	//		how ENTER behaves in the editor window.  It traps the ENTER key and alters
	//		the way DOM is constructed in certain cases to try to commonize the generated
	//		DOM and behaviors across browsers.
	//
	// description:
	//		This plugin has three modes:
	//
	//		- blockNodeForEnter=BR
	//		- blockNodeForEnter=DIV
	//		- blockNodeForEnter=P
	//
	//		In blockNodeForEnter=P, the ENTER key starts a new
	//		paragraph, and shift-ENTER starts a new line in the current paragraph.
	//		For example, the input:
	//
	//	|	first paragraph <shift-ENTER>
	//	|	second line of first paragraph <ENTER>
	//	|	second paragraph
	//
	//		will generate:
	//
	//	|	<p>
	//	|		first paragraph
	//	|		<br/>
	//	|		second line of first paragraph
	//	|	</p>
	//	|	<p>
	//	|		second paragraph
	//	|	</p>
	//
	//		In BR and DIV mode, the ENTER key conceptually goes to a new line in the
	//		current paragraph, and users conceptually create a new paragraph by pressing ENTER twice.
	//		For example, if the user enters text into an editor like this:
	//
	//	|		one <ENTER>
	//	|		two <ENTER>
	//	|		three <ENTER>
	//	|		<ENTER>
	//	|		four <ENTER>
	//	|		five <ENTER>
	//	|		six <ENTER>
	//
	//		It will appear on the screen as two 'paragraphs' of three lines each.  Markupwise, this generates:
	//
	//		BR:
	//	|		one<br/>
	//	|		two<br/>
	//	|		three<br/>
	//	|		<br/>
	//	|		four<br/>
	//	|		five<br/>
	//	|		six<br/>
	//
	//		DIV:
	//	|		<div>one</div>
	//	|		<div>two</div>
	//	|		<div>three</div>
	//	|		<div>&nbsp;</div>
	//	|		<div>four</div>
	//	|		<div>five</div>
	//	|		<div>six</div>

	// blockNodeForEnter: String
	//		This property decides the behavior of Enter key. It can be either P,
	//		DIV, BR, or empty (which means disable this feature). Anything else
	//		will trigger errors.  The default is 'BR'
	//
	//		See class description for more details.
	blockNodeForEnter: 'BR',

	constructor: function(args){
		if(args){
			if("blockNodeForEnter" in args){
				args.blockNodeForEnter = args.blockNodeForEnter.toUpperCase();
			}
			lang.mixin(this,args);
		}
	},

	setEditor: function(editor){
		// Overrides _Plugin.setEditor().
		if(this.editor === editor){ return; }
		this.editor = editor;
		if(this.blockNodeForEnter == 'BR'){
			// While Moz has a mode tht mostly works, it's still a little different,
			// So, try to just have a common mode and be consistent.  Which means
			// we need to enable customUndo, if not already enabled.
			this.editor.customUndo = true;
				editor.onLoadDeferred.then(lang.hitch(this,function(d){
					this.connect(editor.document, "onkeypress", function(e){
						if(e.charOrCode == keys.ENTER){
							// Just do it manually.  The handleEnterKey has a shift mode that
							// Always acts like <br>, so just use it.
							var ne = lang.mixin({},e);
							ne.shiftKey = true;
							if(!this.handleEnterKey(ne)){
								event.stop(e);
							}
						}
					});
					if(has("ie") >= 9){
						this.connect(editor.document, "onpaste", function(e){
							setTimeout(dojo.hitch(this, function(){
								// Use the old range/selection code to kick IE 9 into updating
								// its range by moving it back, then forward, one 'character'.
								var r = this.editor.document.selection.createRange();
								r.move('character',-1);
								r.select();
								r.move('character',1);
								r.select();
							}),0);
						});
					}
					return d;
				}));
		}else if(this.blockNodeForEnter){
			// add enter key handler
			// FIXME: need to port to the new event code!!
			var h = lang.hitch(this,this.handleEnterKey);
			editor.addKeyHandler(13, 0, 0, h); //enter
			editor.addKeyHandler(13, 0, 1, h); //shift+enter
			this.connect(this.editor,'onKeyPressed','onKeyPressed');
		}
	},
	onKeyPressed: function(){
		// summary:
		//		Handler for keypress events.
		// tags:
		//		private
		if(this._checkListLater){
			if(win.withGlobal(this.editor.window, 'isCollapsed', baseFocus)){
				var liparent = this.editor._sCall('getAncestorElement', ['LI']);
				if(!liparent){
					// circulate the undo detection code by calling RichText::execCommand directly
					RichText.prototype.execCommand.call(this.editor, 'formatblock',this.blockNodeForEnter);
					// set the innerHTML of the new block node
					var block = this.editor._sCall('getAncestorElement', [this.blockNodeForEnter]);
					if(block){
						block.innerHTML=this.bogusHtmlContent;
						if(has("ie") <= 9){
							// move to the start by moving backwards one char
							var r = this.editor.document.selection.createRange();
							r.move('character',-1);
							r.select();
						}
					}else{
						console.error('onKeyPressed: Cannot find the new block node'); // FIXME
					}
				}else{
					if(has("mozilla")){
						if(liparent.parentNode.parentNode.nodeName == 'LI'){
							liparent=liparent.parentNode.parentNode;
						}
					}
					var fc=liparent.firstChild;
					if(fc && fc.nodeType == 1 && (fc.nodeName == 'UL' || fc.nodeName == 'OL')){
						liparent.insertBefore(fc.ownerDocument.createTextNode('\xA0'),fc);
						var newrange = rangeapi.create(this.editor.window);
						newrange.setStart(liparent.firstChild,0);
						var selection = rangeapi.getSelection(this.editor.window, true);
						selection.removeAllRanges();
						selection.addRange(newrange);
					}
				}
			}
			this._checkListLater = false;
		}
		if(this._pressedEnterInBlock){
			// the new created is the original current P, so we have previousSibling below
			if(this._pressedEnterInBlock.previousSibling){
				this.removeTrailingBr(this._pressedEnterInBlock.previousSibling);
			}
			delete this._pressedEnterInBlock;
		}
	},

	// bogusHtmlContent: [private] String
	//		HTML to stick into a new empty block
	bogusHtmlContent: '&#160;',		// &nbsp;

	// blockNodes: [private] Regex
	//		Regex for testing if a given tag is a block level (display:block) tag
	blockNodes: /^(?:P|H1|H2|H3|H4|H5|H6|LI)$/,

	handleEnterKey: function(e){
		// summary:
		//		Handler for enter key events when blockNodeForEnter is DIV or P.
		// description:
		//		Manually handle enter key event to make the behavior consistent across
		//		all supported browsers. See class description for details.
		// tags:
		//		private

		var selection, range, newrange, startNode, endNode, brNode, doc=this.editor.document,br,rs,txt;
		if(e.shiftKey){		// shift+enter always generates <br>
			var parent = this.editor._sCall('getParentElement', []);
			var header = rangeapi.getAncestor(parent,this.blockNodes);
			if(header){
				if(header.tagName == 'LI'){
					return true; // let browser handle
				}
				selection = rangeapi.getSelection(this.editor.window);
				range = selection.getRangeAt(0);
				if(!range.collapsed){
					range.deleteContents();
					selection = rangeapi.getSelection(this.editor.window);
					range = selection.getRangeAt(0);
				}
				if(rangeapi.atBeginningOfContainer(header, range.startContainer, range.startOffset)){
						br=doc.createElement('br');
						newrange = rangeapi.create(this.editor.window);
						header.insertBefore(br,header.firstChild);
						newrange.setStartAfter(br);
						selection.removeAllRanges();
						selection.addRange(newrange);
				}else if(rangeapi.atEndOfContainer(header, range.startContainer, range.startOffset)){
					newrange = rangeapi.create(this.editor.window);
					br=doc.createElement('br');
						header.appendChild(br);
						header.appendChild(doc.createTextNode('\xA0'));
						newrange.setStart(header.lastChild,0);
					selection.removeAllRanges();
					selection.addRange(newrange);
				}else{
					rs = range.startContainer;
					if(rs && rs.nodeType == 3){
						// Text node, we have to split it.
						txt = rs.nodeValue;
						startNode = doc.createTextNode(txt.substring(0, range.startOffset));
						endNode = doc.createTextNode(txt.substring(range.startOffset));
						brNode = doc.createElement("br");

						if(endNode.nodeValue == "" && has("webkit")){
							endNode = doc.createTextNode('\xA0')
						}
						domConstruct.place(startNode, rs, "after");
						domConstruct.place(brNode, startNode, "after");
						domConstruct.place(endNode, brNode, "after");
						domConstruct.destroy(rs);
						newrange = rangeapi.create(this.editor.window);
						newrange.setStart(endNode,0);
						selection.removeAllRanges();
						selection.addRange(newrange);
						return false;
					}
					return true; // let browser handle
				}
			}else{
				selection = rangeapi.getSelection(this.editor.window);
				if(selection.rangeCount){
					range = selection.getRangeAt(0);
					if(range && range.startContainer){
						if(!range.collapsed){
							range.deleteContents();
							selection = rangeapi.getSelection(this.editor.window);
							range = selection.getRangeAt(0);
						}
						rs = range.startContainer;
						if(rs && rs.nodeType == 3){
							// Text node, we have to split it.
							var endEmpty = false;

							var offset = range.startOffset;
							if(rs.length < offset){
								//We are not splitting the right node, try to locate the correct one
								ret = this._adjustNodeAndOffset(rs, offset);
								rs = ret.node;
								offset = ret.offset;
							}
							txt = rs.nodeValue;

							startNode = doc.createTextNode(txt.substring(0, offset));
							endNode = doc.createTextNode(txt.substring(offset));
							brNode = doc.createElement("br");

							if(!endNode.length){
								endNode = doc.createTextNode('\xA0');
								endEmpty = true;
							}

							if(startNode.length){
								domConstruct.place(startNode, rs, "after");
							}else{
								startNode = rs;
							}
							domConstruct.place(brNode, startNode, "after");
							domConstruct.place(endNode, brNode, "after");
							domConstruct.destroy(rs);
							newrange = rangeapi.create(this.editor.window);
							newrange.setStart(endNode,0);
							newrange.setEnd(endNode, endNode.length);
							selection.removeAllRanges();
							selection.addRange(newrange);
							if(endEmpty && !has("webkit")){
								this.editor._sCall("remove", []);
							}else{
								this.editor._sCall("collapse", [true]);
							}
						}else{
							var targetNode;
							if(range.startOffset >= 0){
								targetNode = rs.childNodes[range.startOffset];
							}
							var brNode = doc.createElement("br");
							var endNode = doc.createTextNode('\xA0');
							if(!targetNode){
								rs.appendChild(brNode);
								rs.appendChild(endNode);
							}else{
								domConstruct.place(brNode, targetNode, "before");
								domConstruct.place(endNode, brNode, "after");
							}
							newrange = rangeapi.create(this.editor.window);
							newrange.setStart(endNode,0);
							newrange.setEnd(endNode, endNode.length);
							selection.removeAllRanges();
							selection.addRange(newrange);
							this.editor._sCall("collapse", [true]);
						}
					}
				}else{
					// don't change this: do not call this.execCommand, as that may have other logic in subclass
					RichText.prototype.execCommand.call(this.editor, 'inserthtml', '<br>');
				}
			}
			return false;
		}
		var _letBrowserHandle = true;

		// first remove selection
		selection = rangeapi.getSelection(this.editor.window);
		range = selection.getRangeAt(0);
		if(!range.collapsed){
			range.deleteContents();
			selection = rangeapi.getSelection(this.editor.window);
			range = selection.getRangeAt(0);
		}

		var block = rangeapi.getBlockAncestor(range.endContainer, null, this.editor.editNode);
		var blockNode = block.blockNode;

		// if this is under a LI or the parent of the blockNode is LI, just let browser to handle it
		if((this._checkListLater = (blockNode && (blockNode.nodeName == 'LI' || blockNode.parentNode.nodeName == 'LI')))){
			if(has("mozilla")){
				// press enter in middle of P may leave a trailing <br/>, let's remove it later
				this._pressedEnterInBlock = blockNode;
			}
			// if this li only contains spaces, set the content to empty so the browser will outdent this item
			if(/^(\s|&nbsp;|&#160;|\xA0|<span\b[^>]*\bclass=['"]Apple-style-span['"][^>]*>(\s|&nbsp;|&#160;|\xA0)<\/span>)?(<br>)?$/.test(blockNode.innerHTML)){
				// empty LI node
				blockNode.innerHTML = '';
				if(has("webkit")){ // WebKit tosses the range when innerHTML is reset
					newrange = rangeapi.create(this.editor.window);
					newrange.setStart(blockNode, 0);
					selection.removeAllRanges();
					selection.addRange(newrange);
				}
				this._checkListLater = false; // nothing to check since the browser handles outdent
			}
			return true;
		}

		// text node directly under body, let's wrap them in a node
		if(!block.blockNode || block.blockNode===this.editor.editNode){
			try{
				RichText.prototype.execCommand.call(this.editor, 'formatblock',this.blockNodeForEnter);
			}catch(e2){ /*squelch FF3 exception bug when editor content is a single BR*/ }
			// get the newly created block node
			// FIXME
			block = {blockNode: this.editor._sCall('getAncestorElement', [this.blockNodeForEnter]),
					blockContainer: this.editor.editNode};
			if(block.blockNode){
				if(block.blockNode != this.editor.editNode &&
					(!(block.blockNode.textContent || block.blockNode.innerHTML).replace(/^\s+|\s+$/g, "").length)){
					this.removeTrailingBr(block.blockNode);
					return false;
				}
			}else{	// we shouldn't be here if formatblock worked
				block.blockNode = this.editor.editNode;
			}
			selection = rangeapi.getSelection(this.editor.window);
			range = selection.getRangeAt(0);
		}

		var newblock = doc.createElement(this.blockNodeForEnter);
		newblock.innerHTML=this.bogusHtmlContent;
		this.removeTrailingBr(block.blockNode);
		var endOffset = range.endOffset;
		var node = range.endContainer;
		if(node.length < endOffset){
			//We are not checking the right node, try to locate the correct one
			var ret = this._adjustNodeAndOffset(node, endOffset);
			node = ret.node;
			endOffset = ret.offset;
		}
		if(rangeapi.atEndOfContainer(block.blockNode, node, endOffset)){
			if(block.blockNode === block.blockContainer){
				block.blockNode.appendChild(newblock);
			}else{
				domConstruct.place(newblock, block.blockNode, "after");
			}
			_letBrowserHandle = false;
			// lets move caret to the newly created block
			newrange = rangeapi.create(this.editor.window);
			newrange.setStart(newblock, 0);
			selection.removeAllRanges();
			selection.addRange(newrange);
			if(this.editor.height){
				winUtils.scrollIntoView(newblock);
			}
		}else if(rangeapi.atBeginningOfContainer(block.blockNode,
				range.startContainer, range.startOffset)){
			domConstruct.place(newblock, block.blockNode, block.blockNode === block.blockContainer ? "first" : "before");
			if(newblock.nextSibling && this.editor.height){
				// position input caret - mostly WebKit needs this
				newrange = rangeapi.create(this.editor.window);
				newrange.setStart(newblock.nextSibling, 0);
				selection.removeAllRanges();
				selection.addRange(newrange);
				// browser does not scroll the caret position into view, do it manually
				winUtils.scrollIntoView(newblock.nextSibling);
			}
			_letBrowserHandle = false;
		}else{ //press enter in the middle of P/DIV/Whatever/
			if(block.blockNode === block.blockContainer){
				block.blockNode.appendChild(newblock);
			}else{
				domConstruct.place(newblock, block.blockNode, "after");
			}
			_letBrowserHandle = false;

			// Clone any block level styles.
			if(block.blockNode.style){
				if(newblock.style){
					if(block.blockNode.style.cssText){
						newblock.style.cssText = block.blockNode.style.cssText;
					}
				}
			}

			// Okay, we probably have to split.
			rs = range.startContainer;
			var firstNodeMoved;
			if(rs && rs.nodeType == 3){
				// Text node, we have to split it.
				var nodeToMove, tNode;
				endOffset = range.endOffset;
				if(rs.length < endOffset){
					//We are not splitting the right node, try to locate the correct one
					ret = this._adjustNodeAndOffset(rs, endOffset);
					rs = ret.node;
					endOffset = ret.offset;
				}

				txt = rs.nodeValue;
				startNode = doc.createTextNode(txt.substring(0, endOffset));
				endNode = doc.createTextNode(txt.substring(endOffset, txt.length));

				// Place the split, then remove original nodes.
				domConstruct.place(startNode, rs, "before");
				domConstruct.place(endNode, rs, "after");
				domConstruct.destroy(rs);

				// Okay, we split the text.  Now we need to see if we're
				// parented to the block element we're splitting and if
				// not, we have to split all the way up.  Ugh.
				var parentC = startNode.parentNode;
				while(parentC !== block.blockNode){
					var tg = parentC.tagName;
					var newTg = doc.createElement(tg);
					// Clone over any 'style' data.
					if(parentC.style){
						if(newTg.style){
							if(parentC.style.cssText){
								newTg.style.cssText = parentC.style.cssText;
							}
						}
					}
					// If font also need to clone over any font data.
					if(parentC.tagName === "FONT"){
						if(parentC.color){
							newTg.color = parentC.color;
						}
						if(parentC.face){
							newTg.face = parentC.face;
						}
						if(parentC.size){  // this check was necessary on IE
							newTg.size = parentC.size;
						}
					}

					nodeToMove = endNode;
					while(nodeToMove){
						tNode = nodeToMove.nextSibling;
						newTg.appendChild(nodeToMove);
						nodeToMove = tNode;
					}
					domConstruct.place(newTg, parentC, "after");
					startNode = parentC;
					endNode = newTg;
					parentC = parentC.parentNode;
				}

				// Lastly, move the split out tags to the new block.
				// as they should now be split properly.
				nodeToMove = endNode;
				if(nodeToMove.nodeType == 1 || (nodeToMove.nodeType == 3 && nodeToMove.nodeValue)){
					// Non-blank text and non-text nodes need to clear out that blank space
					// before moving the contents.
					newblock.innerHTML = "";
				}
				firstNodeMoved = nodeToMove;
				while(nodeToMove){
					tNode = nodeToMove.nextSibling;
					newblock.appendChild(nodeToMove);
					nodeToMove = tNode;
				}
			}

			//lets move caret to the newly created block
			newrange = rangeapi.create(this.editor.window);
			var nodeForCursor;
			var innerMostFirstNodeMoved = firstNodeMoved;
			if(this.blockNodeForEnter !== 'BR'){
				while(innerMostFirstNodeMoved){
					nodeForCursor = innerMostFirstNodeMoved;
					tNode = innerMostFirstNodeMoved.firstChild;
					innerMostFirstNodeMoved = tNode;
				}
				if(nodeForCursor && nodeForCursor.parentNode){
					newblock = nodeForCursor.parentNode;
					newrange.setStart(newblock, 0);
					selection.removeAllRanges();
					selection.addRange(newrange);
					if(this.editor.height){
						winUtils.scrollIntoView(newblock);
					}
					if(has("mozilla")){
						// press enter in middle of P may leave a trailing <br/>, let's remove it later
						this._pressedEnterInBlock = block.blockNode;
					}
				}else{
					_letBrowserHandle = true;
				}
			}else{
				newrange.setStart(newblock, 0);
				selection.removeAllRanges();
				selection.addRange(newrange);
				if(this.editor.height){
					winUtils.scrollIntoView(newblock);
				}
				if(has("mozilla")){
					// press enter in middle of P may leave a trailing <br/>, let's remove it later
					this._pressedEnterInBlock = block.blockNode;
				}
			}
		}
		return _letBrowserHandle;
	},

	_adjustNodeAndOffset: function(/*DomNode*/node, /*Int*/offset){
		// summary:
		//		In the case there are multiple text nodes in a row the offset may not be within the node.  If the offset is larger than the node length, it will attempt to find
		//		the next text sibling until it locates the text node in which the offset refers to
		// node:
		//		The node to check.
		// offset:
		//		The position to find within the text node
		// tags:
		//		private.
		while(node.length < offset && node.nextSibling && node.nextSibling.nodeType==3){
			//Adjust the offset and node in the case of multiple text nodes in a row
			offset = offset - node.length;
			node = node.nextSibling;
		}
		return {"node": node, "offset": offset};
	},

	removeTrailingBr: function(container){
		// summary:
		//		If last child of container is a `<br>`, then remove it.
		// tags:
		//		private
		var para = /P|DIV|LI/i.test(container.tagName) ?
			container : this.editor._sCall("getParentOfType", [container,['P','DIV','LI']]);

		if(!para){ return; }
		if(para.lastChild){
			if((para.childNodes.length > 1 && para.lastChild.nodeType == 3 && /^[\s\xAD]*$/.test(para.lastChild.nodeValue)) ||
				para.lastChild.tagName=='BR'){

				domConstruct.destroy(para.lastChild);
			}
		}
		if(!para.childNodes.length){
			para.innerHTML=this.bogusHtmlContent;
		}
	}
});

});

},
'dijit/_editor/RichText':function(){
define("dijit/_editor/RichText", [
	"dojo/_base/array", // array.forEach array.indexOf array.some
	"dojo/_base/config", // config
	"dojo/_base/declare", // declare
	"dojo/_base/Deferred", // Deferred
	"dojo/dom", // dom.byId
	"dojo/dom-attr", // domAttr.set or get
	"dojo/dom-class", // domClass.add domClass.remove
	"dojo/dom-construct", // domConstruct.create domConstruct.destroy domConstruct.place
	"dojo/dom-geometry", // domGeometry.position
	"dojo/dom-style", // domStyle.getComputedStyle domStyle.set
	"dojo/_base/event", // event.stop
	"dojo/_base/kernel", // kernel.deprecated
	"dojo/keys", // keys.BACKSPACE keys.TAB
	"dojo/_base/lang", // lang.clone lang.hitch lang.isArray lang.isFunction lang.isString lang.trim
	"dojo/on", // on()
	"dojo/query", // query
	"dojo/ready", // ready
	"dojo/sniff", // has("ie") has("mozilla") has("opera") has("safari") has("webkit")
	"dojo/topic",	// topic.publish() (publish)
	"dojo/_base/unload", // unload
	"dojo/_base/url", // url
	"dojo/_base/window", // win.global
	"../_Widget",
	"../_CssStateMixin",
	"./selection",
	"./range",
	"./html",
	"../focus",
	"../main"	// dijit._scopeName
], function(array, config, declare, Deferred, dom, domAttr, domClass, domConstruct, domGeometry, domStyle,
	event, kernel, keys, lang, on, query, ready, has, topic, unload, _Url, win,
	_Widget, _CssStateMixin, selectionapi, rangeapi, htmlapi, focus, dijit){

// module:
//		dijit/_editor/RichText
// summary:
//		dijit/_editor/RichText is the core of dijit/Editor, which provides basic
//		WYSIWYG editing features.

// if you want to allow for rich text saving with back/forward actions, you must add a text area to your page with
// the id==dijit._scopeName + "._editor.RichText.value" (typically "dijit/_editor/RichText.value). For example,
// something like this will work:
//
//	<textarea id="dijit._editor.RichText.value" style="display:none;position:absolute;top:-100px;left:-100px;height:3px;width:3px;overflow:hidden;"></textarea>
//

var RichText = declare("dijit._editor.RichText", [_Widget, _CssStateMixin], {
	// summary:
	//		dijit/_editor/RichText is the core of dijit.Editor, which provides basic
	//		WYSIWYG editing features.
	//
	// description:
	//		dijit/_editor/RichText is the core of dijit.Editor, which provides basic
	//		WYSIWYG editing features. It also encapsulates the differences
	//		of different js engines for various browsers.  Do not use this widget
	//		with an HTML &lt;TEXTAREA&gt; tag, since the browser unescapes XML escape characters,
	//		like &lt;.  This can have unexpected behavior and lead to security issues
	//		such as scripting attacks.
	//
	// tags:
	//		private

	constructor: function(params /*===== , srcNodeRef =====*/){
		// summary:
		//		Create the widget.
		// params: Object|null
		//		Initial settings for any of the widget attributes, except readonly attributes.
		// srcNodeRef: DOMNode
		//		The widget replaces the specified DOMNode.

		// contentPreFilters: Function(String)[]
		//		Pre content filter function register array.
		//		these filters will be executed before the actual
		//		editing area gets the html content.
		this.contentPreFilters = [];

		// contentPostFilters: Function(String)[]
		//		post content filter function register array.
		//		These will be used on the resulting html
		//		from contentDomPostFilters. The resulting
		//		content is the final html (returned by getValue()).
		this.contentPostFilters = [];

		// contentDomPreFilters: Function(DomNode)[]
		//		Pre content dom filter function register array.
		//		These filters are applied after the result from
		//		contentPreFilters are set to the editing area.
		this.contentDomPreFilters = [];

		// contentDomPostFilters: Function(DomNode)[]
		//		Post content dom filter function register array.
		//		These filters are executed on the editing area dom.
		//		The result from these will be passed to contentPostFilters.
		this.contentDomPostFilters = [];

		// editingAreaStyleSheets: dojo._URL[]
		//		array to store all the stylesheets applied to the editing area
		this.editingAreaStyleSheets = [];

		// Make a copy of this.events before we start writing into it, otherwise we
		// will modify the prototype which leads to bad things on pages w/multiple editors
		this.events = [].concat(this.events);

		this._keyHandlers = {};

		if(params && lang.isString(params.value)){
			this.value = params.value;
		}

		this.onLoadDeferred = new Deferred();
	},

	baseClass: "dijitEditor",

	// inheritWidth: Boolean
	//		whether to inherit the parent's width or simply use 100%
	inheritWidth: false,

	// focusOnLoad: [deprecated] Boolean
	//		Focus into this widget when the page is loaded
	focusOnLoad: false,

	// name: String?
	//		Specifies the name of a (hidden) `<textarea>` node on the page that's used to save
	//		the editor content on page leave.   Used to restore editor contents after navigating
	//		to a new page and then hitting the back button.
	name: "",

	// styleSheets: [const] String
	//		semicolon (";") separated list of css files for the editing area
	styleSheets: "",

	// height: String
	//		Set height to fix the editor at a specific height, with scrolling.
	//		By default, this is 300px.  If you want to have the editor always
	//		resizes to accommodate the content, use AlwaysShowToolbar plugin
	//		and set height="".  If this editor is used within a layout widget,
	//		set height="100%".
	height: "300px",

	// minHeight: String
	//		The minimum height that the editor should have.
	minHeight: "1em",

	// isClosed: [private] Boolean
	isClosed: true,

	// isLoaded: [private] Boolean
	isLoaded: false,

	// _SEPARATOR: [private] String
	//		Used to concat contents from multiple editors into a single string,
	//		so they can be saved into a single `<textarea>` node.  See "name" attribute.
	_SEPARATOR: "@@**%%__RICHTEXTBOUNDRY__%%**@@",

	// _NAME_CONTENT_SEP: [private] String
	//		USed to separate name from content.  Just a colon isn't safe.
	_NAME_CONTENT_SEP: "@@**%%:%%**@@",

	// onLoadDeferred: [readonly] dojo.Deferred
	//		Deferred which is fired when the editor finishes loading.
	//		Call myEditor.onLoadDeferred.then(callback) it to be informed
	//		when the rich-text area initialization is finalized.
	onLoadDeferred: null,

	// isTabIndent: Boolean
	//		Make tab key and shift-tab indent and outdent rather than navigating.
	//		Caution: sing this makes web pages inaccessible to users unable to use a mouse.
	isTabIndent: false,

	// disableSpellCheck: [const] Boolean
	//		When true, disables the browser's native spell checking, if supported.
	//		Works only in Firefox.
	disableSpellCheck: false,

	postCreate: function(){
		if("textarea" === this.domNode.tagName.toLowerCase()){
			console.warn("RichText should not be used with the TEXTAREA tag.  See dijit._editor.RichText docs.");
		}

		// Push in the builtin filters now, making them the first executed, but not over-riding anything
		// users passed in.  See: #6062
		this.contentPreFilters = [lang.hitch(this, "_preFixUrlAttributes")].concat(this.contentPreFilters);
		if(has("mozilla")){
			this.contentPreFilters = [this._normalizeFontStyle].concat(this.contentPreFilters);
			this.contentPostFilters = [this._removeMozBogus].concat(this.contentPostFilters);
		}
		if(has("webkit")){
			// Try to clean up WebKit bogus artifacts.  The inserted classes
			// made by WebKit sometimes messes things up.
			this.contentPreFilters = [this._removeWebkitBogus].concat(this.contentPreFilters);
			this.contentPostFilters = [this._removeWebkitBogus].concat(this.contentPostFilters);
		}
		if(has("ie")){
			// IE generates <strong> and <em> but we want to normalize to <b> and <i>
			this.contentPostFilters = [this._normalizeFontStyle].concat(this.contentPostFilters);
			this.contentDomPostFilters = [lang.hitch(this, this._stripBreakerNodes)].concat(this.contentDomPostFilters);
		}
		this.inherited(arguments);

		topic.publish(dijit._scopeName + "._editor.RichText::init", this);
		this.open();
		this.setupDefaultShortcuts();
	},

	setupDefaultShortcuts: function(){
		// summary:
		//		Add some default key handlers
		// description:
		//		Overwrite this to setup your own handlers. The default
		//		implementation does not use Editor commands, but directly
		//		executes the builtin commands within the underlying browser
		//		support.
		// tags:
		//		protected
		var exec = lang.hitch(this, function(cmd, arg){
			return function(){
				return !this.execCommand(cmd,arg);
			};
		});

		var ctrlKeyHandlers = {
			b: exec("bold"),
			i: exec("italic"),
			u: exec("underline"),
			a: exec("selectall"),
			s: function(){ this.save(true); },
			m: function(){ this.isTabIndent = !this.isTabIndent; },

			"1": exec("formatblock", "h1"),
			"2": exec("formatblock", "h2"),
			"3": exec("formatblock", "h3"),
			"4": exec("formatblock", "h4"),

			"\\": exec("insertunorderedlist")
		};

		if(!has("ie")){
			ctrlKeyHandlers.Z = exec("redo"); //FIXME: undo?
		}

		var key;
		for(key in ctrlKeyHandlers){
			this.addKeyHandler(key, true, false, ctrlKeyHandlers[key]);
		}
	},

	// events: [private] String[]
	//		 events which should be connected to the underlying editing area
	events: ["onKeyPress", "onKeyDown", "onKeyUp"], // onClick handled specially

	// captureEvents: [deprecated] String[]
	//		 Events which should be connected to the underlying editing
	//		 area, events in this array will be addListener with
	//		 capture=true.
	// TODO: looking at the code I don't see any distinction between events and captureEvents,
	// so get rid of this for 2.0 if not sooner
	captureEvents: [],

	_editorCommandsLocalized: false,
	_localizeEditorCommands: function(){
		// summary:
		//		When IE is running in a non-English locale, the API actually changes,
		//		so that we have to say (for example) danraku instead of p (for paragraph).
		//		Handle that here.
		// tags:
		//		private
		if(RichText._editorCommandsLocalized){
			// Use the already generate cache of mappings.
			this._local2NativeFormatNames = RichText._local2NativeFormatNames;
			this._native2LocalFormatNames = RichText._native2LocalFormatNames;
			return;
		}
		RichText._editorCommandsLocalized = true;
		RichText._local2NativeFormatNames = {};
		RichText._native2LocalFormatNames = {};
		this._local2NativeFormatNames = RichText._local2NativeFormatNames;
		this._native2LocalFormatNames = RichText._native2LocalFormatNames;
		//in IE, names for blockformat is locale dependent, so we cache the values here

		//put p after div, so if IE returns Normal, we show it as paragraph
		//We can distinguish p and div if IE returns Normal, however, in order to detect that,
		//we have to call this.document.selection.createRange().parentElement() or such, which
		//could slow things down. Leave it as it is for now
		var formats = ['div', 'p', 'pre', 'h1', 'h2', 'h3', 'h4', 'h5', 'h6', 'ol', 'ul', 'address'];
		var localhtml = "", format, i=0;
		while((format=formats[i++])){
			//append a <br> after each element to separate the elements more reliably
			if(format.charAt(1) !== 'l'){
				localhtml += "<"+format+"><span>content</span></"+format+"><br/>";
			}else{
				localhtml += "<"+format+"><li>content</li></"+format+"><br/>";
			}
		}
		// queryCommandValue returns empty if we hide editNode, so move it out of screen temporary
		// Also, IE9 does weird stuff unless we do it inside the editor iframe.
		var style = { position: "absolute", top: "0px", zIndex: 10, opacity: 0.01 };
		var div = domConstruct.create('div', {style: style, innerHTML: localhtml});
		this.ownerDocumentBody.appendChild(div);

		// IE9 has a timing issue with doing this right after setting
		// the inner HTML, so put a delay in.
		var inject = lang.hitch(this, function(){
			var node = div.firstChild;
			while(node){
				try{
					this._sCall("selectElement", [node.firstChild]);
					var nativename = node.tagName.toLowerCase();
					this._local2NativeFormatNames[nativename] = document.queryCommandValue("formatblock");
					this._native2LocalFormatNames[this._local2NativeFormatNames[nativename]] = nativename;
					node = node.nextSibling.nextSibling;
					//console.log("Mapped: ", nativename, " to: ", this._local2NativeFormatNames[nativename]);
				}catch(e){ /*Sqelch the occasional IE9 error */ }
			}
			domConstruct.destroy(div);
		});
		this.defer(inject);
	},

	open: function(/*DomNode?*/ element){
		// summary:
		//		Transforms the node referenced in this.domNode into a rich text editing
		//		node.
		// description:
		//		Sets up the editing area asynchronously. This will result in
		//		the creation and replacement with an iframe.
		// tags:
		//		private

		if(!this.onLoadDeferred || this.onLoadDeferred.fired >= 0){
			this.onLoadDeferred = new Deferred();
		}

		if(!this.isClosed){ this.close(); }
		topic.publish(dijit._scopeName + "._editor.RichText::open", this);

		if(arguments.length === 1 && element.nodeName){ // else unchanged
			this.domNode = element;
		}

		var dn = this.domNode;

		// "html" will hold the innerHTML of the srcNodeRef and will be used to
		// initialize the editor.
		var html;

		if(lang.isString(this.value)){
			// Allow setting the editor content programmatically instead of
			// relying on the initial content being contained within the target
			// domNode.
			html = this.value;
			delete this.value;
			dn.innerHTML = "";
		}else if(dn.nodeName && dn.nodeName.toLowerCase() == "textarea"){
			// if we were created from a textarea, then we need to create a
			// new editing harness node.
			var ta = (this.textarea = dn);
			this.name = ta.name;
			html = ta.value;
			dn = this.domNode = this.ownerDocument.createElement("div");
			dn.setAttribute('widgetId', this.id);
			ta.removeAttribute('widgetId');
			dn.cssText = ta.cssText;
			dn.className += " " + ta.className;
			domConstruct.place(dn, ta, "before");
			var tmpFunc = lang.hitch(this, function(){
				//some browsers refuse to submit display=none textarea, so
				//move the textarea off screen instead
				domStyle.set(ta, {
					display: "block",
					position: "absolute",
					top: "-1000px"
				});

				if(has("ie")){ //nasty IE bug: abnormal formatting if overflow is not hidden
					var s = ta.style;
					this.__overflow = s.overflow;
					s.overflow = "hidden";
				}
			});
			if(has("ie")){
				this.defer(tmpFunc, 10);
			}else{
				tmpFunc();
			}

			if(ta.form){
				var resetValue = ta.value;
				this.reset = function(){
					var current = this.getValue();
					if(current !== resetValue){
						this.replaceValue(resetValue);
					}
				};
				on(ta.form, "submit", lang.hitch(this, function(){
					// Copy value to the <textarea> so it gets submitted along with form.
					// FIXME: should we be calling close() here instead?
					domAttr.set(ta, 'disabled', this.disabled); // don't submit the value if disabled
					ta.value = this.getValue();
				}));
			}
		}else{
			html = htmlapi.getChildrenHtml(dn);
			dn.innerHTML = "";
		}

		this.value = html;

		// If we're a list item we have to put in a blank line to force the
		// bullet to nicely align at the top of text
		if(dn.nodeName && dn.nodeName === "LI"){
			dn.innerHTML = " <br>";
		}

		// Construct the editor div structure.
		this.header = dn.ownerDocument.createElement("div");
		dn.appendChild(this.header);
		this.editingArea = dn.ownerDocument.createElement("div");
		dn.appendChild(this.editingArea);
		this.footer = dn.ownerDocument.createElement("div");
		dn.appendChild(this.footer);

		if(!this.name){
			this.name = this.id + "_AUTOGEN";
		}

		// User has pressed back/forward button so we lost the text in the editor, but it's saved
		// in a hidden <textarea> (which contains the data for all the editors on this page),
		// so get editor value from there
		if(this.name !== "" && (!config["useXDomain"] || config["allowXdRichTextSave"])){
			var saveTextarea = dom.byId(dijit._scopeName + "._editor.RichText.value");
			if(saveTextarea && saveTextarea.value !== ""){
				var datas = saveTextarea.value.split(this._SEPARATOR), i=0, dat;
				while((dat=datas[i++])){
					var data = dat.split(this._NAME_CONTENT_SEP);
					if(data[0] === this.name){
						html = data[1];
						datas = datas.splice(i, 1);
						saveTextarea.value = datas.join(this._SEPARATOR);
						break;
					}
				}
			}

			if(!RichText._globalSaveHandler){
				RichText._globalSaveHandler = {};
				unload.addOnUnload(function(){
					var id;
					for(id in RichText._globalSaveHandler){
						var f = RichText._globalSaveHandler[id];
						if(lang.isFunction(f)){
							f();
						}
					}
				});
			}
			RichText._globalSaveHandler[this.id] = lang.hitch(this, "_saveContent");
		}

		this.isClosed = false;

		var ifr = (this.editorObject = this.iframe = this.ownerDocument.createElement('iframe'));
		ifr.id = this.id+"_iframe";
		ifr.style.border = "none";
		ifr.style.width = "100%";
		if(this._layoutMode){
			// iframe should be 100% height, thus getting it's height from surrounding
			// <div> (which has the correct height set by Editor)
			ifr.style.height = "100%";
		}else{
			if(has("ie") >= 7){
				if(this.height){
					ifr.style.height = this.height;
				}
				if(this.minHeight){
					ifr.style.minHeight = this.minHeight;
				}
			}else{
				ifr.style.height = this.height ? this.height : this.minHeight;
			}
		}
		ifr.frameBorder = 0;
		ifr._loadFunc = lang.hitch( this, function(w){
			this.window = w;
			this.document = this.window.document;

			if(has("ie")){
				this._localizeEditorCommands();
			}

			// Do final setup and set initial contents of editor
			this.onLoad(html);
		});

		// Set the iframe's initial (blank) content.
		var src = this._getIframeDocTxt(),
			s = "javascript: '" + src.replace(/\\/g, "\\\\").replace(/'/g, "\\'") + "'";
		ifr.setAttribute('src', s);
		this.editingArea.appendChild(ifr);

		if(has("safari") <= 4){
			src = ifr.getAttribute("src");
			if(!src || src.indexOf("javascript") === -1){
				// Safari 4 and earlier sometimes act oddly
				// So we have to set it again.
				this.defer(function(){ ifr.setAttribute('src', s); });
			}
		}

		// TODO: this is a guess at the default line-height, kinda works
		if(dn.nodeName === "LI"){
			dn.lastChild.style.marginTop = "-1.2em";
		}

		domClass.add(this.domNode, this.baseClass);
	},

	//static cache variables shared among all instance of this class
	_local2NativeFormatNames: {},
	_native2LocalFormatNames: {},

	_getIframeDocTxt: function(){
		// summary:
		//		Generates the boilerplate text of the document inside the iframe (ie, `<html><head>...</head><body/></html>`).
		//		Editor content (if not blank) should be added afterwards.
		// tags:
		//		private
		var _cs = domStyle.getComputedStyle(this.domNode);

		// The contents inside of <body>.  The real contents are set later via a call to setValue().
		var html = "";
		var setBodyId = true;
		if(has("ie") || has("webkit") || (!this.height && !has("mozilla"))){
			// In auto-expand mode, need a wrapper div for AlwaysShowToolbar plugin to correctly
			// expand/contract the editor as the content changes.
			html = "<div id='dijitEditorBody'></div>";
			setBodyId = false;
		}else if(has("mozilla")){
			// workaround bug where can't select then delete text (until user types something
			// into the editor)... and/or issue where typing doesn't erase selected text
			this._cursorToStart = true;
			html = "&#160;";	// &nbsp;
		}

		var font = [ _cs.fontWeight, _cs.fontSize, _cs.fontFamily ].join(" ");

		// line height is tricky - applying a units value will mess things up.
		// if we can't get a non-units value, bail out.
		var lineHeight = _cs.lineHeight;
		if(lineHeight.indexOf("px") >= 0){
			lineHeight = parseFloat(lineHeight)/parseFloat(_cs.fontSize);
			// console.debug(lineHeight);
		}else if(lineHeight.indexOf("em")>=0){
			lineHeight = parseFloat(lineHeight);
		}else{
			// If we can't get a non-units value, just default
			// it to the CSS spec default of 'normal'.  Seems to
			// work better, esp on IE, than '1.0'
			lineHeight = "normal";
		}
		var userStyle = "";
		var self = this;
		this.style.replace(/(^|;)\s*(line-|font-?)[^;]+/ig, function(match){
			match = match.replace(/^;/ig,"") + ';';
			var s = match.split(":")[0];
			if(s){
				s = lang.trim(s);
				s = s.toLowerCase();
				var i;
				var sC = "";
				for(i = 0; i < s.length; i++){
					var c = s.charAt(i);
					switch(c){
						case "-":
							i++;
							c = s.charAt(i).toUpperCase();
						default:
							sC += c;
					}
				}
				domStyle.set(self.domNode, sC, "");
			}
			userStyle += match + ';';
		});


		// need to find any associated label element and update iframe document title
		var label=query('label[for="'+this.id+'"]');

		return [
			this.isLeftToRight() ? "<html>\n<head>\n" : "<html dir='rtl'>\n<head>\n",
			(has("mozilla") && label.length ? "<title>" + label[0].innerHTML + "</title>\n" : ""),
			"<meta http-equiv='Content-Type' content='text/html'>\n",
			"<style>\n",
			"\tbody,html {\n",
			"\t\tbackground:transparent;\n",
			"\t\tpadding: 1px 0 0 0;\n",
			"\t\tmargin: -1px 0 0 0;\n", // remove extraneous vertical scrollbar on safari and firefox

			// Set the html/body sizing.  Webkit always needs this, other browsers
			// only set it when height is defined (not auto-expanding), otherwise
			// scrollers do not appear.
			((has("webkit"))?"\t\twidth: 100%;\n":""),
			((has("webkit"))?"\t\theight: 100%;\n":""),
			"\t}\n",

			// TODO: left positioning will cause contents to disappear out of view
			//	   if it gets too wide for the visible area
			"\tbody{\n",
			"\t\ttop:0px;\n",
			"\t\tleft:0px;\n",
			"\t\tright:0px;\n",
			"\t\tfont:", font, ";\n",
				((this.height||has("opera")) ? "" : "\t\tposition: fixed;\n"),
			// FIXME: IE 6 won't understand min-height?
			"\t\tmin-height:", this.minHeight, ";\n",
			"\t\tline-height:", lineHeight,";\n",
			"\t}\n",
			"\tp{ margin: 1em 0; }\n",

			// Determine how scrollers should be applied.  In autoexpand mode (height = "") no scrollers on y at all.
			// But in fixed height mode we want both x/y scrollers.  Also, if it's using wrapping div and in auto-expand
			// (Mainly IE) we need to kill the y scroller on body and html.
			(!setBodyId && !this.height ? "\tbody,html {overflow-y: hidden;}\n" : ""),
			"\t#dijitEditorBody{overflow-x: auto; overflow-y:" + (this.height ? "auto;" : "hidden;") + " outline: 0px;}\n",
			"\tli > ul:-moz-first-node, li > ol:-moz-first-node{ padding-top: 1.2em; }\n",
			// Can't set min-height in IE9, it puts layout on li, which puts move/resize handles.
			(!has("ie") ? "\tli{ min-height:1.2em; }\n" : ""),
			"</style>\n",
			this._applyEditingAreaStyleSheets(),"\n",
			"</head>\n<body ",
			(setBodyId?"id='dijitEditorBody' ":""),

			// Onload handler fills in real editor content.
			// On IE9, sometimes onload is called twice, and the first time frameElement is null (test_FullScreen.html)
			"onload='frameElement && frameElement._loadFunc(window,document)' ",
			"style='"+userStyle+"'>", html, "</body>\n</html>"
		].join(""); // String
	},

	_applyEditingAreaStyleSheets: function(){
		// summary:
		//		apply the specified css files in styleSheets
		// tags:
		//		private
		var files = [];
		if(this.styleSheets){
			files = this.styleSheets.split(';');
			this.styleSheets = '';
		}

		//empty this.editingAreaStyleSheets here, as it will be filled in addStyleSheet
		files = files.concat(this.editingAreaStyleSheets);
		this.editingAreaStyleSheets = [];

		var text='', i=0, url;
		while((url=files[i++])){
			var abstring = (new _Url(win.global.location, url)).toString();
			this.editingAreaStyleSheets.push(abstring);
			text += '<link rel="stylesheet" type="text/css" href="'+abstring+'"/>';
		}
		return text;
	},

	addStyleSheet: function(/*dojo/_base/url*/ uri){
		// summary:
		//		add an external stylesheet for the editing area
		// uri:
		//		Url of the external css file
		var url=uri.toString();

		//if uri is relative, then convert it to absolute so that it can be resolved correctly in iframe
		if(url.charAt(0) === '.' || (url.charAt(0) !== '/' && !uri.host)){
			url = (new _Url(win.global.location, url)).toString();
		}

		if(array.indexOf(this.editingAreaStyleSheets, url) > -1){
//			console.debug("dijit/_editor/RichText.addStyleSheet(): Style sheet "+url+" is already applied");
			return;
		}

		this.editingAreaStyleSheets.push(url);
		this.onLoadDeferred.then(lang.hitch(this, function(){
			if(this.document.createStyleSheet){ //IE
				this.document.createStyleSheet(url);
			}else{ //other browser
				var head = this.document.getElementsByTagName("head")[0];
				var stylesheet = this.document.createElement("link");
				stylesheet.rel="stylesheet";
				stylesheet.type="text/css";
				stylesheet.href=url;
				head.appendChild(stylesheet);
			}
		}));
	},

	removeStyleSheet: function(/*dojo/_base/url*/ uri){
		// summary:
		//		remove an external stylesheet for the editing area
		var url=uri.toString();
		//if uri is relative, then convert it to absolute so that it can be resolved correctly in iframe
		if(url.charAt(0) === '.' || (url.charAt(0) !== '/' && !uri.host)){
			url = (new _Url(win.global.location, url)).toString();
		}
		var index = array.indexOf(this.editingAreaStyleSheets, url);
		if(index === -1){
//			console.debug("dijit/_editor/RichText.removeStyleSheet(): Style sheet "+url+" has not been applied");
			return;
		}
		delete this.editingAreaStyleSheets[index];
		query('link:[href="'+url+'"]', this.window.document).orphan();
	},

	// disabled: Boolean
	//		The editor is disabled; the text cannot be changed.
	disabled: false,

	_mozSettingProps: {'styleWithCSS':false},
	_setDisabledAttr: function(/*Boolean*/ value){
		value = !!value;
		this._set("disabled", value);
		if(!this.isLoaded){ return; } // this method requires init to be complete
		if(has("ie") || has("webkit") || has("opera")){
			var preventIEfocus = has("ie") && (this.isLoaded || !this.focusOnLoad);
			if(preventIEfocus){ this.editNode.unselectable = "on"; }
			this.editNode.contentEditable = !value;
			if(preventIEfocus){
				this.defer(function(){
					if(this.editNode){		// guard in case widget destroyed before timeout
						this.editNode.unselectable = "off";
					}
				});
			}
		}else{ //moz
			try{
				this.document.designMode=(value?'off':'on');
			}catch(e){ return; } // ! _disabledOK
			if(!value && this._mozSettingProps){
				var ps = this._mozSettingProps;
				var n;
				for(n in ps){
					if(ps.hasOwnProperty(n)){
						try{
							this.document.execCommand(n,false,ps[n]);
						}catch(e2){}
					}
				}
			}
//			this.document.execCommand('contentReadOnly', false, value);
//				if(value){
//					this.blur(); //to remove the blinking caret
//				}
		}
		this._disabledOK = true;
	},

/* Event handlers
 *****************/

	onLoad: function(/*String*/ html){
		// summary:
		//		Handler after the iframe finishes loading.
		// html: String
		//		Editor contents should be set to this value
		// tags:
		//		protected

		// TODO: rename this to _onLoad, make empty public onLoad() method, deprecate/make protected onLoadDeferred handler?

		if(!this.window.__registeredWindow){
			this.window.__registeredWindow = true;
			this._iframeRegHandle = focus.registerIframe(this.iframe);
		}
		if(!has("ie") && !has("webkit") && (this.height || has("mozilla"))){
			this.editNode=this.document.body;
		}else{
			// there's a wrapper div around the content, see _getIframeDocTxt().
			this.editNode=this.document.body.firstChild;
			var _this = this;
			if(has("ie")){ // #4996 IE wants to focus the BODY tag
				this.tabStop = domConstruct.create('div', { tabIndex: -1 }, this.editingArea);
				this.iframe.onfocus = function(){ _this.editNode.setActive(); };
			}
		}
		this.focusNode = this.editNode; // for InlineEditBox


		var events = this.events.concat(this.captureEvents);
		var ap = this.iframe ? this.document : this.editNode;
		array.forEach(events, function(item){
			this.connect(ap, item.toLowerCase(), item);
		}, this);

		this.connect(ap, "onmouseup", "onClick"); // mouseup in the margin does not generate an onclick event

		if(has("ie")){ // IE contentEditable
			this.connect(this.document, "onmousedown", "_onIEMouseDown"); // #4996 fix focus

			// give the node Layout on IE
			// TODO: this may no longer be needed, since we've reverted IE to using an iframe,
			// not contentEditable.   Removing it would also probably remove the need for creating
			// the extra <div> in _getIframeDocTxt()
			this.editNode.style.zoom = 1.0;
		}else{
			this.connect(this.document, "onmousedown", function(){
				// Clear the moveToStart focus, as mouse
				// down will set cursor point.  Required to properly
				// work with selection/position driven plugins and clicks in
				// the window. refs: #10678
				delete this._cursorToStart;
			});
		}

		if(has("webkit")){
			//WebKit sometimes doesn't fire right on selections, so the toolbar
			//doesn't update right.  Therefore, help it out a bit with an additional
			//listener.  A mouse up will typically indicate a display change, so fire this
			//and get the toolbar to adapt.  Reference: #9532
			this._webkitListener = this.connect(this.document, "onmouseup", "onDisplayChanged");
			this.connect(this.document, "onmousedown", function(e){
				var t = e.target;
				if(t && (t === this.document.body || t === this.document)){
					// Since WebKit uses the inner DIV, we need to check and set position.
					// See: #12024 as to why the change was made.
					this.defer("placeCursorAtEnd");
				}
			});
		}

		if(has("ie")){
			// Try to make sure 'hidden' elements aren't visible in edit mode (like browsers other than IE
			// do).  See #9103
			try{
				this.document.execCommand('RespectVisibilityInDesign', true, null);
			}catch(e){/* squelch */}
		}

		this.isLoaded = true;

		this.set('disabled', this.disabled); // initialize content to editable (or not)

		// Note that setValue() call will only work after isLoaded is set to true (above)

		// Set up a function to allow delaying the setValue until a callback is fired
		// This ensures extensions like dijit.Editor have a way to hold the value set
		// until plugins load (and do things like register filters).
		var setContent = lang.hitch(this, function(){
			this.setValue(html);
			if(this.onLoadDeferred){
				this.onLoadDeferred.resolve(true);
			}
			this.onDisplayChanged();
			if(this.focusOnLoad){
				// after the document loads, then set focus after updateInterval expires so that
				// onNormalizedDisplayChanged has run to avoid input caret issues
				ready(lang.hitch(this, "defer", "focus", this.updateInterval));
			}
			// Save off the initial content now
			this.value = this.getValue(true);
		});
		if(this.setValueDeferred){
			this.setValueDeferred.then(setContent);
		}else{
			setContent();
		}
	},

	onKeyDown: function(/* Event */ e){
		// summary:
		//		Handler for onkeydown event
		// tags:
		//		protected

		// we need this event at the moment to get the events from control keys
		// such as the backspace. It might be possible to add this to Dojo, so that
		// keyPress events can be emulated by the keyDown and keyUp detection.

		if(e.keyCode === keys.TAB && this.isTabIndent ){
			event.stop(e); //prevent tab from moving focus out of editor

			// FIXME: this is a poor-man's indent/outdent. It would be
			// better if it added 4 "&nbsp;" chars in an undoable way.
			// Unfortunately pasteHTML does not prove to be undoable
			if(this.queryCommandEnabled((e.shiftKey ? "outdent" : "indent"))){
				this.execCommand((e.shiftKey ? "outdent" : "indent"));
			}
		}
		if(has("ie")){
			if(e.keyCode == keys.TAB && !this.isTabIndent){
				if(e.shiftKey && !e.ctrlKey && !e.altKey){
					// focus the BODY so the browser will tab away from it instead
					this.iframe.focus();
				}else if(!e.shiftKey && !e.ctrlKey && !e.altKey){
					// focus the BODY so the browser will tab away from it instead
					this.tabStop.focus();
				}
			}else if(e.keyCode === keys.BACKSPACE && this.document.selection.type === "Control"){
				// IE has a bug where if a non-text object is selected in the editor,
				// hitting backspace would act as if the browser's back button was
				// clicked instead of deleting the object. see #1069
				event.stop(e);
				this.execCommand("delete");
			}else if((65 <= e.keyCode && e.keyCode <= 90) ||
				(e.keyCode>=37 && e.keyCode<=40) // FIXME: get this from connect() instead!
			){ //arrow keys
				e.charCode = e.keyCode;
				this.onKeyPress(e);
			}
		}
		if(has("ff")){
			if(e.keyCode === keys.PAGE_UP || e.keyCode === keys.PAGE_DOWN ){
				if(this.editNode.clientHeight >= this.editNode.scrollHeight){
					// Stop the event to prevent firefox from trapping the cursor when there is no scroll bar.
					e.preventDefault();
				}
			}
		}
		return true;
	},

	onKeyUp: function(/*===== e =====*/){
		// summary:
		//		Handler for onkeyup event
		// tags:
		//		callback
	},

	setDisabled: function(/*Boolean*/ disabled){
		// summary:
		//		Deprecated, use set('disabled', ...) instead.
		// tags:
		//		deprecated
		kernel.deprecated('dijit.Editor::setDisabled is deprecated','use dijit.Editor::attr("disabled",boolean) instead', 2.0);
		this.set('disabled',disabled);
	},
	_setValueAttr: function(/*String*/ value){
		// summary:
		//		Registers that attr("value", foo) should call setValue(foo)
		this.setValue(value);
	},
	_setDisableSpellCheckAttr: function(/*Boolean*/ disabled){
		if(this.document){
			domAttr.set(this.document.body, "spellcheck", !disabled);
		}else{
			// try again after the editor is finished loading
			this.onLoadDeferred.then(lang.hitch(this, function(){
				domAttr.set(this.document.body, "spellcheck", !disabled);
			}));
		}
		this._set("disableSpellCheck", disabled);
	},

	onKeyPress: function(e){
		// summary:
		//		Handle the various key events
		// tags:
		//		protected

		var c = (e.keyChar && e.keyChar.toLowerCase()) || e.keyCode,
			handlers = this._keyHandlers[c],
			args = arguments;
			
		if(handlers && !e.altKey){
			array.some(handlers, function(h){
				// treat meta- same as ctrl-, for benefit of mac users
				if(!(h.shift ^ e.shiftKey) && !(h.ctrl ^ (e.ctrlKey||e.metaKey))){ 
					if(!h.handler.apply(this, args)){
						e.preventDefault();
					}
					return true;
				}
			}, this);
		}

		// function call after the character has been inserted
		if(!this._onKeyHitch){
			this._onKeyHitch = lang.hitch(this, "onKeyPressed");
		}
		this.defer("_onKeyHitch", 1);
		return true;
	},

	addKeyHandler: function(/*String*/ key, /*Boolean*/ ctrl, /*Boolean*/ shift, /*Function*/ handler){
		// summary:
		//		Add a handler for a keyboard shortcut
		// description:
		//		The key argument should be in lowercase if it is a letter character
		// tags:
		//		protected
		if(!lang.isArray(this._keyHandlers[key])){
			this._keyHandlers[key] = [];
		}
		//TODO: would be nice to make this a hash instead of an array for quick lookups
		this._keyHandlers[key].push({
			shift: shift || false,
			ctrl: ctrl || false,
			handler: handler
		});
	},

	onKeyPressed: function(){
		// summary:
		//		Handler for after the user has pressed a key, and the display has been updated.
		//		(Runs on a timer so that it runs after the display is updated)
		// tags:
		//		private
		this.onDisplayChanged(/*e*/); // can't pass in e
	},

	onClick: function(/*Event*/ e){
		// summary:
		//		Handler for when the user clicks.
		// tags:
		//		private

		// console.info('onClick',this._tryDesignModeOn);
		this.onDisplayChanged(e);
	},

	_onIEMouseDown: function(){
		// summary:
		//		IE only to prevent 2 clicks to focus
		// tags:
		//		protected

		if(!this.focused && !this.disabled){
			this.focus();
		}
	},

	_onBlur: function(e){
		// summary:
		//		Called from focus manager when focus has moved away from this editor
		// tags:
		//		protected

		// console.info('_onBlur')

		this.inherited(arguments);

		var newValue = this.getValue(true);
		if(newValue !== this.value){
			this.onChange(newValue);
		}
		this._set("value", newValue);
	},

	_onFocus: function(/*Event*/ e){
		// summary:
		//		Called from focus manager when focus has moved into this editor
		// tags:
		//		protected

		// console.info('_onFocus')
		if(!this.disabled){
			if(!this._disabledOK){
				this.set('disabled', false);
			}
			this.inherited(arguments);
		}
	},

	// TODO: remove in 2.0
	blur: function(){
		// summary:
		//		Remove focus from this instance.
		// tags:
		//		deprecated
		if(!has("ie") && this.window.document.documentElement && this.window.document.documentElement.focus){
			this.window.document.documentElement.focus();
		}else if(this.ownerDocumentBody.focus){
			this.ownerDocumentBody.focus();
		}
	},

	focus: function(){
		// summary:
		//		Move focus to this editor
		if(!this.isLoaded){
			this.focusOnLoad = true;
			return;
		}
		if(this._cursorToStart){
			delete this._cursorToStart;
			if(this.editNode.childNodes){
				this.placeCursorAtStart(); // this calls focus() so return
				return;
			}
		}
		if(!has("ie")){
			focus.focus(this.iframe);
		}else if(this.editNode && this.editNode.focus){
			// editNode may be hidden in display:none div, lets just punt in this case
			//this.editNode.focus(); -> causes IE to scroll always (strict and quirks mode) to the top the Iframe
			// if we fire the event manually and let the browser handle the focusing, the latest
			// cursor position is focused like in FF
			this.iframe.fireEvent('onfocus', document.createEventObject()); // createEventObject only in IE
		//	}else{
		// TODO: should we throw here?
		// console.debug("Have no idea how to focus into the editor!");
		}
	},

	// _lastUpdate: 0,
	updateInterval: 200,
	_updateTimer: null,
	onDisplayChanged: function(/*Event*/ /*===== e =====*/){
		// summary:
		//		This event will be fired every time the display context
		//		changes and the result needs to be reflected in the UI.
		// description:
		//		If you don't want to have update too often,
		//		onNormalizedDisplayChanged should be used instead
		// tags:
		//		private

		// var _t=new Date();
		if(this._updateTimer){
			this._updateTimer.remove();
		}
		this._updateTimer = this.defer("onNormalizedDisplayChanged", this.updateInterval);

		// Technically this should trigger a call to watch("value", ...) registered handlers,
		// but getValue() is too slow to call on every keystroke so we don't.
	},
	onNormalizedDisplayChanged: function(){
		// summary:
		//		This event is fired every updateInterval ms or more
		// description:
		//		If something needs to happen immediately after a
		//		user change, please use onDisplayChanged instead.
		// tags:
		//		private
		delete this._updateTimer;
	},
	onChange: function(/*===== newContent =====*/){
		// summary:
		//		This is fired if and only if the editor loses focus and
		//		the content is changed.
	},
	_normalizeCommand: function(/*String*/ cmd, /*Anything?*/argument){
		// summary:
		//		Used as the advice function to map our
		//		normalized set of commands to those supported by the target
		//		browser.
		// tags:
		//		private

		var command = cmd.toLowerCase();
		if(command === "formatblock"){
			if(has("safari") && argument === undefined){ command = "heading"; }
		}else if(command === "hilitecolor" && !has("mozilla")){
			command = "backcolor";
		}

		return command;
	},

	_qcaCache: {},
	queryCommandAvailable: function(/*String*/ command){
		// summary:
		//		Tests whether a command is supported by the host. Clients
		//		SHOULD check whether a command is supported before attempting
		//		to use it, behaviour for unsupported commands is undefined.
		// command:
		//		The command to test for
		// tags:
		//		private

		// memoizing version. See _queryCommandAvailable for computing version
		var ca = this._qcaCache[command];
		if(ca !== undefined){ return ca; }
		return (this._qcaCache[command] = this._queryCommandAvailable(command));
	},

	_queryCommandAvailable: function(/*String*/ command){
		// summary:
		//		See queryCommandAvailable().
		// tags:
		//		private

		var ie = 1;
		var mozilla = 1 << 1;
		var webkit = 1 << 2;
		var opera = 1 << 3;

		function isSupportedBy(browsers){
			return {
				ie: Boolean(browsers & ie),
				mozilla: Boolean(browsers & mozilla),
				webkit: Boolean(browsers & webkit),
				opera: Boolean(browsers & opera)
			};
		}

		var supportedBy = null;

		switch(command.toLowerCase()){
			case "bold": case "italic": case "underline":
			case "subscript": case "superscript":
			case "fontname": case "fontsize":
			case "forecolor": case "hilitecolor":
			case "justifycenter": case "justifyfull": case "justifyleft":
			case "justifyright": case "delete": case "selectall": case "toggledir":
				supportedBy = isSupportedBy(mozilla | ie | webkit | opera);
				break;

			case "createlink": case "unlink": case "removeformat":
			case "inserthorizontalrule": case "insertimage":
			case "insertorderedlist": case "insertunorderedlist":
			case "indent": case "outdent": case "formatblock":
			case "inserthtml": case "undo": case "redo": case "strikethrough": case "tabindent":
				supportedBy = isSupportedBy(mozilla | ie | opera | webkit);
				break;

			case "blockdirltr": case "blockdirrtl":
			case "dirltr": case "dirrtl":
			case "inlinedirltr": case "inlinedirrtl":
				supportedBy = isSupportedBy(ie);
				break;
			case "cut": case "copy": case "paste":
				supportedBy = isSupportedBy( ie | mozilla | webkit | opera);
				break;

			case "inserttable":
				supportedBy = isSupportedBy(mozilla | ie);
				break;

			case "insertcell": case "insertcol": case "insertrow":
			case "deletecells": case "deletecols": case "deleterows":
			case "mergecells": case "splitcell":
				supportedBy = isSupportedBy(ie | mozilla);
				break;

			default: return false;
		}

		return (has("ie") && supportedBy.ie) ||
			(has("mozilla") && supportedBy.mozilla) ||
			(has("webkit") && supportedBy.webkit) ||
			(has("opera") && supportedBy.opera);	// Boolean return true if the command is supported, false otherwise
	},

	execCommand: function(/*String*/ command, argument){
		// summary:
		//		Executes a command in the Rich Text area
		// command:
		//		The command to execute
		// argument:
		//		An optional argument to the command
		// tags:
		//		protected
		var returnValue;

		//focus() is required for IE to work
		//In addition, focus() makes sure after the execution of
		//the command, the editor receives the focus as expected
		this.focus();

		command = this._normalizeCommand(command, argument);
		
		if(argument !== undefined){
			if(command === "heading"){
				throw new Error("unimplemented");
			}else if((command === "formatblock") && has("ie")){
				argument = '<'+argument+'>';
			}
		}

		//Check to see if we have any over-rides for commands, they will be functions on this
		//widget of the form _commandImpl.  If we don't, fall through to the basic native
		//exec command of the browser.
		var implFunc = "_" + command + "Impl";
		if(this[implFunc]){
			returnValue = this[implFunc](argument);
		}else{
			argument = arguments.length > 1 ? argument : null;
			if(argument || command !== "createlink"){
				returnValue = this.document.execCommand(command, false, argument);
			}
		}

		this.onDisplayChanged();
		return returnValue;
	},

	queryCommandEnabled: function(/*String*/ command){
		// summary:
		//		Check whether a command is enabled or not.
		// command:
		//		The command to execute
		// tags:
		//		protected
		if(this.disabled || !this._disabledOK){ return false; }

		command = this._normalizeCommand(command);

		//Check to see if we have any over-rides for commands, they will be functions on this
		//widget of the form _commandEnabledImpl.  If we don't, fall through to the basic native
		//command of the browser.
		var implFunc = "_" + command + "EnabledImpl";

		if(this[implFunc]){
			return  this[implFunc](command);
		}else{
			return this._browserQueryCommandEnabled(command);
		}
	},

	queryCommandState: function(command){
		// summary:
		//		Check the state of a given command and returns true or false.
		// tags:
		//		protected

		if(this.disabled || !this._disabledOK){ return false; }
		command = this._normalizeCommand(command);
		try{
			return this.document.queryCommandState(command);
		}catch(e){
			//Squelch, occurs if editor is hidden on FF 3 (and maybe others.)
			return false;
		}
	},

	queryCommandValue: function(command){
		// summary:
		//		Check the value of a given command. This matters most for
		//		custom selections and complex values like font value setting.
		// tags:
		//		protected

		if(this.disabled || !this._disabledOK){ return false; }
		var r;
		command = this._normalizeCommand(command);
		if(has("ie") && command === "formatblock"){
			r = this._native2LocalFormatNames[this.document.queryCommandValue(command)];
		}else if(has("mozilla") && command === "hilitecolor"){
			var oldValue;
			try{
				oldValue = this.document.queryCommandValue("styleWithCSS");
			}catch(e){
				oldValue = false;
			}
			this.document.execCommand("styleWithCSS", false, true);
			r = this.document.queryCommandValue(command);
			this.document.execCommand("styleWithCSS", false, oldValue);
		}else{
			r = this.document.queryCommandValue(command);
		}
		return r;
	},

	// Misc.

	_sCall: function(name, args){
		// summary:
		//		Run the named method of dijit/_editor/selection over the
		//		current editor instance's window, with the passed args.
		// tags:
		//		private
		return win.withGlobal(this.window, name, selectionapi, args);
	},

	// FIXME: this is a TON of code duplication. Why?

	placeCursorAtStart: function(){
		// summary:
		//		Place the cursor at the start of the editing area.
		// tags:
		//		private

		this.focus();

		//see comments in placeCursorAtEnd
		var isvalid=false;
		if(has("mozilla")){
			// TODO:  Is this branch even necessary?
			var first=this.editNode.firstChild;
			while(first){
				if(first.nodeType === 3){
					if(first.nodeValue.replace(/^\s+|\s+$/g, "").length>0){
						isvalid=true;
						this._sCall("selectElement", [ first ]);
						break;
					}
				}else if(first.nodeType === 1){
					isvalid=true;
					var tg = first.tagName ? first.tagName.toLowerCase() : "";
					// Collapse before childless tags.
					if(/br|input|img|base|meta|area|basefont|hr|link/.test(tg)){
						this._sCall("selectElement", [ first ]);
					}else{
						// Collapse inside tags with children.
						this._sCall("selectElementChildren", [ first ]);
					}
					break;
				}
				first = first.nextSibling;
			}
		}else{
			isvalid=true;
			this._sCall("selectElementChildren", [ this.editNode ]);
		}
		if(isvalid){
			this._sCall("collapse", [ true ]);
		}
	},

	placeCursorAtEnd: function(){
		// summary:
		//		Place the cursor at the end of the editing area.
		// tags:
		//		private

		this.focus();

		//In mozilla, if last child is not a text node, we have to use
		// selectElementChildren on this.editNode.lastChild otherwise the
		// cursor would be placed at the end of the closing tag of
		//this.editNode.lastChild
		var isvalid=false;
		if(has("mozilla")){
			var last=this.editNode.lastChild;
			while(last){
				if(last.nodeType === 3){
					if(last.nodeValue.replace(/^\s+|\s+$/g, "").length>0){
						isvalid=true;
						this._sCall("selectElement", [ last ]);
						break;
					}
				}else if(last.nodeType === 1){
					isvalid=true;
					this._sCall("selectElement", [ last.lastChild || last]);
					break;
				}
				last = last.previousSibling;
			}
		}else{
			isvalid=true;
			this._sCall("selectElementChildren", [ this.editNode ]);
		}
		if(isvalid){
			this._sCall("collapse", [ false ]);
		}
	},

	getValue: function(/*Boolean?*/ nonDestructive){
		// summary:
		//		Return the current content of the editing area (post filters
		//		are applied).  Users should call get('value') instead.
		// nonDestructive:
		//		defaults to false. Should the post-filtering be run over a copy
		//		of the live DOM? Most users should pass "true" here unless they
		//		*really* know that none of the installed filters are going to
		//		mess up the editing session.
		// tags:
		//		private
		if(this.textarea){
			if(this.isClosed || !this.isLoaded){
				return this.textarea.value;
			}
		}

		return this._postFilterContent(null, nonDestructive);
	},
	_getValueAttr: function(){
		// summary:
		//		Hook to make attr("value") work
		return this.getValue(true);
	},

	setValue: function(/*String*/ html){
		// summary:
		//		This function sets the content. No undo history is preserved.
		//		Users should use set('value', ...) instead.
		// tags:
		//		deprecated

		// TODO: remove this and getValue() for 2.0, and move code to _setValueAttr()

		if(!this.isLoaded){
			// try again after the editor is finished loading
			this.onLoadDeferred.then(lang.hitch(this, function(){
				this.setValue(html);
			}));
			return;
		}
		this._cursorToStart = true;
		if(this.textarea && (this.isClosed || !this.isLoaded)){
			this.textarea.value=html;
		}else{
			html = this._preFilterContent(html);
			var node = this.isClosed ? this.domNode : this.editNode;
			if(html && has("mozilla") && html.toLowerCase() === "<p></p>"){
				html = "<p>&#160;</p>";	// &nbsp;
			}

			// Use &nbsp; to avoid webkit problems where editor is disabled until the user clicks it
			if(!html && has("webkit")){
				html = "&#160;";	// &nbsp;
			}
			node.innerHTML = html;
			this._preDomFilterContent(node);
		}

		this.onDisplayChanged();
		this._set("value", this.getValue(true));
	},

	replaceValue: function(/*String*/ html){
		// summary:
		//		This function set the content while trying to maintain the undo stack
		//		(now only works fine with Moz, this is identical to setValue in all
		//		other browsers)
		// tags:
		//		protected

		if(this.isClosed){
			this.setValue(html);
		}else if(this.window && this.window.getSelection && !has("mozilla")){ // Safari
			// look ma! it's a totally f'd browser!
			this.setValue(html);
		}else if(this.window && this.window.getSelection){ // Moz
			html = this._preFilterContent(html);
			this.execCommand("selectall");
			if(!html){
				this._cursorToStart = true;
				html = "&#160;";	// &nbsp;
			}
			this.execCommand("inserthtml", html);
			this._preDomFilterContent(this.editNode);
		}else if(this.document && this.document.selection){//IE
			//In IE, when the first element is not a text node, say
			//an <a> tag, when replacing the content of the editing
			//area, the <a> tag will be around all the content
			//so for now, use setValue for IE too
			this.setValue(html);
		}

		this._set("value", this.getValue(true));
	},

	_preFilterContent: function(/*String*/ html){
		// summary:
		//		Filter the input before setting the content of the editing
		//		area. DOM pre-filtering may happen after this
		//		string-based filtering takes place but as of 1.2, this is not
		//		guaranteed for operations such as the inserthtml command.
		// tags:
		//		private

		var ec = html;
		array.forEach(this.contentPreFilters, function(ef){ if(ef){ ec = ef(ec); } });
		return ec;
	},
	_preDomFilterContent: function(/*DomNode*/ dom){
		// summary:
		//		filter the input's live DOM. All filter operations should be
		//		considered to be "live" and operating on the DOM that the user
		//		will be interacting with in their editing session.
		// tags:
		//		private
		dom = dom || this.editNode;
		array.forEach(this.contentDomPreFilters, function(ef){
			if(ef && lang.isFunction(ef)){
				ef(dom);
			}
		}, this);
	},

	_postFilterContent: function(
		/*DomNode|DomNode[]|String?*/ dom,
		/*Boolean?*/ nonDestructive){
		// summary:
		//		filter the output after getting the content of the editing area
		//
		// description:
		//		post-filtering allows plug-ins and users to specify any number
		//		of transforms over the editor's content, enabling many common
		//		use-cases such as transforming absolute to relative URLs (and
		//		vice-versa), ensuring conformance with a particular DTD, etc.
		//		The filters are registered in the contentDomPostFilters and
		//		contentPostFilters arrays. Each item in the
		//		contentDomPostFilters array is a function which takes a DOM
		//		Node or array of nodes as its only argument and returns the
		//		same. It is then passed down the chain for further filtering.
		//		The contentPostFilters array behaves the same way, except each
		//		member operates on strings. Together, the DOM and string-based
		//		filtering allow the full range of post-processing that should
		//		be necessaray to enable even the most agressive of post-editing
		//		conversions to take place.
		//
		//		If nonDestructive is set to "true", the nodes are cloned before
		//		filtering proceeds to avoid potentially destructive transforms
		//		to the content which may still needed to be edited further.
		//		Once DOM filtering has taken place, the serialized version of
		//		the DOM which is passed is run through each of the
		//		contentPostFilters functions.
		//
		// dom:
		//		a node, set of nodes, which to filter using each of the current
		//		members of the contentDomPostFilters and contentPostFilters arrays.
		//
		// nonDestructive:
		//		defaults to "false". If true, ensures that filtering happens on
		//		a clone of the passed-in content and not the actual node
		//		itself.
		//
		// tags:
		//		private

		var ec;
		if(!lang.isString(dom)){
			dom = dom || this.editNode;
			if(this.contentDomPostFilters.length){
				if(nonDestructive){
					dom = lang.clone(dom);
				}
				array.forEach(this.contentDomPostFilters, function(ef){
					dom = ef(dom);
				});
			}
			ec = htmlapi.getChildrenHtml(dom);
		}else{
			ec = dom;
		}

		if(!lang.trim(ec.replace(/^\xA0\xA0*/, '').replace(/\xA0\xA0*$/, '')).length){
			ec = "";
		}

		//	if(has("ie")){
		//		//removing appended <P>&nbsp;</P> for IE
		//		ec = ec.replace(/(?:<p>&nbsp;</p>[\n\r]*)+$/i,"");
		//	}
		array.forEach(this.contentPostFilters, function(ef){
			ec = ef(ec);
		});

		return ec;
	},

	_saveContent: function(){
		// summary:
		//		Saves the content in an onunload event if the editor has not been closed
		// tags:
		//		private

		var saveTextarea = dom.byId(dijit._scopeName + "._editor.RichText.value");
		if(saveTextarea){
			if(saveTextarea.value){
				saveTextarea.value += this._SEPARATOR;
			}
			saveTextarea.value += this.name + this._NAME_CONTENT_SEP + this.getValue(true);
		}
	},


	escapeXml: function(/*String*/ str, /*Boolean*/ noSingleQuotes){
		// summary:
		//		Adds escape sequences for special characters in XML.
		//		Optionally skips escapes for single quotes
		// tags:
		//		private

		str = str.replace(/&/gm, "&amp;").replace(/</gm, "&lt;").replace(/>/gm, "&gt;").replace(/"/gm, "&quot;");
		if(!noSingleQuotes){
			str = str.replace(/'/gm, "&#39;");
		}
		return str; // string
	},

	getNodeHtml: function(/* DomNode */ node){
		// summary:
		//		Deprecated.   Use dijit/_editor/html::_getNodeHtml() instead.
		// tags:
		//		deprecated
		kernel.deprecated('dijit.Editor::getNodeHtml is deprecated','use dijit/_editor/html::getNodeHtml instead', 2);
		return htmlapi.getNodeHtml(node); // String
	},

	getNodeChildrenHtml: function(/* DomNode */ dom){
		// summary:
		//		Deprecated.   Use dijit/_editor/html::getChildrenHtml() instead.
		// tags:
		//		deprecated
		kernel.deprecated('dijit.Editor::getNodeChildrenHtml is deprecated','use dijit/_editor/html::getChildrenHtml instead', 2);
		return htmlapi.getChildrenHtml(dom);
	},

	close: function(/*Boolean?*/ save){
		// summary:
		//		Kills the editor and optionally writes back the modified contents to the
		//		element from which it originated.
		// save:
		//		Whether or not to save the changes. If false, the changes are discarded.
		// tags:
		//		private

		if(this.isClosed){ return; }

		if(!arguments.length){ save = true; }
		if(save){
			this._set("value", this.getValue(true));
		}

		// line height is squashed for iframes
		// FIXME: why was this here? if(this.iframe){ this.domNode.style.lineHeight = null; }

		if(this.interval){ clearInterval(this.interval); }

		if(this._webkitListener){
			//Cleaup of WebKit fix: #9532
			this.disconnect(this._webkitListener);
			delete this._webkitListener;
		}

		// Guard against memory leaks on IE (see #9268)
		if(has("ie")){
			 this.iframe.onfocus = null;
		}
		this.iframe._loadFunc = null;

		if(this._iframeRegHandle){
			this._iframeRegHandle.remove();
			delete this._iframeRegHandle;
		}

		if(this.textarea){
			var s = this.textarea.style;
			s.position = "";
			s.left = s.top = "";
			if(has("ie")){
				s.overflow = this.__overflow;
				this.__overflow = null;
			}
			this.textarea.value = this.value;
			domConstruct.destroy(this.domNode);
			this.domNode = this.textarea;
		}else{
			// Note that this destroys the iframe
			this.domNode.innerHTML = this.value;
		}
		delete this.iframe;

		domClass.remove(this.domNode, this.baseClass);
		this.isClosed = true;
		this.isLoaded = false;

		delete this.editNode;
		delete this.focusNode;

		if(this.window && this.window._frameElement){
			this.window._frameElement = null;
		}

		this.window = null;
		this.document = null;
		this.editingArea = null;
		this.editorObject = null;
	},

	destroy: function(){
		if(!this.isClosed){ this.close(false); }
		if(this._updateTimer){
			this._updateTimer.remove();
		}
		this.inherited(arguments);
		if(RichText._globalSaveHandler){
			delete RichText._globalSaveHandler[this.id];
		}
	},

	_removeMozBogus: function(/* String */ html){
		// summary:
		//		Post filter to remove unwanted HTML attributes generated by mozilla
		// tags:
		//		private
		return html.replace(/\stype="_moz"/gi, '').replace(/\s_moz_dirty=""/gi, '').replace(/_moz_resizing="(true|false)"/gi,''); // String
	},
	_removeWebkitBogus: function(/* String */ html){
		// summary:
		//		Post filter to remove unwanted HTML attributes generated by webkit
		// tags:
		//		private
		html = html.replace(/\sclass="webkit-block-placeholder"/gi, '');
		html = html.replace(/\sclass="apple-style-span"/gi, '');
		// For some reason copy/paste sometime adds extra meta tags for charset on
		// webkit (chrome) on mac.They need to be removed.  See: #12007"
		html = html.replace(/<meta charset=\"utf-8\" \/>/gi, '');
		return html; // String
	},
	_normalizeFontStyle: function(/* String */ html){
		// summary:
		//		Convert 'strong' and 'em' to 'b' and 'i'.
		// description:
		//		Moz can not handle strong/em tags correctly, so to help
		//		mozilla and also to normalize output, convert them to 'b' and 'i'.
		//
		//		Note the IE generates 'strong' and 'em' rather than 'b' and 'i'
		// tags:
		//		private
		return html.replace(/<(\/)?strong([ \>])/gi, '<$1b$2')
			.replace(/<(\/)?em([ \>])/gi, '<$1i$2' ); // String
	},

	_preFixUrlAttributes: function(/* String */ html){
		// summary:
		//		Pre-filter to do fixing to href attributes on `<a>` and `<img>` tags
		// tags:
		//		private
		return html.replace(/(?:(<a(?=\s).*?\shref=)("|')(.*?)\2)|(?:(<a\s.*?href=)([^"'][^ >]+))/gi,
				'$1$4$2$3$5$2 _djrealurl=$2$3$5$2')
			.replace(/(?:(<img(?=\s).*?\ssrc=)("|')(.*?)\2)|(?:(<img\s.*?src=)([^"'][^ >]+))/gi,
				'$1$4$2$3$5$2 _djrealurl=$2$3$5$2'); // String
	},

	/*****************************************************************************
		The following functions implement HTML manipulation commands for various
		browser/contentEditable implementations.  The goal of them is to enforce
		standard behaviors of them.
	******************************************************************************/

	/*** queryCommandEnabled implementations ***/

	_browserQueryCommandEnabled: function(command){
		// summary:
		//		Implementation to call to the native queryCommandEnabled of the browser.
		// command:
		//		The command to check.
		// tags:
		//		protected
		if(!command) { return false; }
		var elem = has("ie") ? this.document.selection.createRange() : this.document;
		try{
			return elem.queryCommandEnabled(command);
		}catch(e){
			return false;
		}
	},

	_createlinkEnabledImpl: function(/*===== argument =====*/){
		// summary:
		//		This function implements the test for if the create link
		//		command should be enabled or not.
		// argument:
		//		arguments to the exec command, if any.
		// tags:
		//		protected
		var enabled = true;
		if(has("opera")){
			var sel = this.window.getSelection();
			if(sel.isCollapsed){
				enabled = true;
			}else{
				enabled = this.document.queryCommandEnabled("createlink");
			}
		}else{
			enabled = this._browserQueryCommandEnabled("createlink");
		}
		return enabled;
	},

	_unlinkEnabledImpl: function(/*===== argument =====*/){
		// summary:
		//		This function implements the test for if the unlink
		//		command should be enabled or not.
		// argument:
		//		arguments to the exec command, if any.
		// tags:
		//		protected
		var enabled = true;
		if(has("mozilla") || has("webkit")){
			enabled = this._sCall("hasAncestorElement", ["a"]);
		}else{
			enabled = this._browserQueryCommandEnabled("unlink");
		}
		return enabled;
	},

	_inserttableEnabledImpl: function(/*===== argument =====*/){
		// summary:
		//		This function implements the test for if the inserttable
		//		command should be enabled or not.
		// argument:
		//		arguments to the exec command, if any.
		// tags:
		//		protected
		var enabled = true;
		if(has("mozilla") || has("webkit")){
			enabled = true;
		}else{
			enabled = this._browserQueryCommandEnabled("inserttable");
		}
		return enabled;
	},

	_cutEnabledImpl: function(/*===== argument =====*/){
		// summary:
		//		This function implements the test for if the cut
		//		command should be enabled or not.
		// argument:
		//		arguments to the exec command, if any.
		// tags:
		//		protected
		var enabled = true;
		if(has("webkit")){
			// WebKit deems clipboard activity as a security threat and natively would return false
			var sel = this.window.getSelection();
			if(sel){ sel = sel.toString(); }
			enabled = !!sel;
		}else{
			enabled = this._browserQueryCommandEnabled("cut");
		}
		return enabled;
	},

	_copyEnabledImpl: function(/*===== argument =====*/){
		// summary:
		//		This function implements the test for if the copy
		//		command should be enabled or not.
		// argument:
		//		arguments to the exec command, if any.
		// tags:
		//		protected
		var enabled = true;
		if(has("webkit")){
			// WebKit deems clipboard activity as a security threat and natively would return false
			var sel = this.window.getSelection();
			if(sel){ sel = sel.toString(); }
			enabled = !!sel;
		}else{
			enabled = this._browserQueryCommandEnabled("copy");
		}
		return enabled;
	},

	_pasteEnabledImpl: function(/*===== argument =====*/){
		// summary:c
		//		This function implements the test for if the paste
		//		command should be enabled or not.
		// argument:
		//		arguments to the exec command, if any.
		// tags:
		//		protected
		var enabled = true;
		if(has("webkit")){
			return true;
		}else{
			enabled = this._browserQueryCommandEnabled("paste");
		}
		return enabled;
	},

	/*** execCommand implementations ***/

	_inserthorizontalruleImpl: function(argument){
		// summary:
		//		This function implements the insertion of HTML 'HR' tags.
		//		into a point on the page.  IE doesn't to it right, so
		//		we have to use an alternate form
		// argument:
		//		arguments to the exec command, if any.
		// tags:
		//		protected
		if(has("ie")){
			return this._inserthtmlImpl("<hr>");
		}
		return this.document.execCommand("inserthorizontalrule", false, argument);
	},

	_unlinkImpl: function(argument){
		// summary:
		//		This function implements the unlink of an 'a' tag.
		// argument:
		//		arguments to the exec command, if any.
		// tags:
		//		protected
		if((this.queryCommandEnabled("unlink")) && (has("mozilla") || has("webkit"))){
			var a = this._sCall("getAncestorElement", [ "a" ]);
			this._sCall("selectElement", [ a ]);
			return this.document.execCommand("unlink", false, null);
		}
		return this.document.execCommand("unlink", false, argument);
	},

	_hilitecolorImpl: function(argument){
		// summary:
		//		This function implements the hilitecolor command
		// argument:
		//		arguments to the exec command, if any.
		// tags:
		//		protected
		var returnValue;
		var isApplied = this._handleTextColorOrProperties("hilitecolor", argument);
		if(!isApplied){
			if(has("mozilla")){
				// mozilla doesn't support hilitecolor properly when useCSS is
				// set to false (bugzilla #279330)
				this.document.execCommand("styleWithCSS", false, true);
				console.log("Executing color command.");
				returnValue = this.document.execCommand("hilitecolor", false, argument);
				this.document.execCommand("styleWithCSS", false, false);
			}else{
				returnValue = this.document.execCommand("hilitecolor", false, argument);
			}
		}
		return returnValue;
	},

	_backcolorImpl: function(argument){
		// summary:
		//		This function implements the backcolor command
		// argument:
		//		arguments to the exec command, if any.
		// tags:
		//		protected
		if(has("ie")){
			// Tested under IE 6 XP2, no problem here, comment out
			// IE weirdly collapses ranges when we exec these commands, so prevent it
			//	var tr = this.document.selection.createRange();
			argument = argument ? argument : null;
		}
		var isApplied = this._handleTextColorOrProperties("backcolor", argument);
		if(!isApplied){
			isApplied = this.document.execCommand("backcolor", false, argument);
		}
		return isApplied;
	},

	_forecolorImpl: function(argument){
		// summary:
		//		This function implements the forecolor command
		// argument:
		//		arguments to the exec command, if any.
		// tags:
		//		protected
		if(has("ie")){
			// Tested under IE 6 XP2, no problem here, comment out
			// IE weirdly collapses ranges when we exec these commands, so prevent it
			//	var tr = this.document.selection.createRange();
			argument = argument? argument : null;
		}
		var isApplied = false;
		isApplied = this._handleTextColorOrProperties("forecolor", argument);
		if(!isApplied){
			isApplied = this.document.execCommand("forecolor", false, argument);
		}
		return isApplied;
	},

	_inserthtmlImpl: function(argument){
		// summary:
		//		This function implements the insertion of HTML content into
		//		a point on the page.
		// argument:
		//		The content to insert, if any.
		// tags:
		//		protected
		argument = this._preFilterContent(argument);
		var rv = true;
		if(has("ie")){
			var insertRange = this.document.selection.createRange();
			if(this.document.selection.type.toUpperCase() === 'CONTROL'){
				var n=insertRange.item(0);
				while(insertRange.length){
					insertRange.remove(insertRange.item(0));
				}
				n.outerHTML=argument;
			}else{
				insertRange.pasteHTML(argument);
			}
			insertRange.select();
			//insertRange.collapse(true);
		}else if(has("mozilla") && !argument.length){
			//mozilla can not inserthtml an empty html to delete current selection
			//so we delete the selection instead in this case
			this._sCall("remove"); // FIXME
		}else{
			rv = this.document.execCommand("inserthtml", false, argument);
		}
		return rv;
	},

	_boldImpl: function(argument){
		// summary:
		//		This function implements an over-ride of the bold command.
		// argument:
		//		Not used, operates by selection.
		// tags:
		//		protected
		var applied = false;
		if(has("ie")){
			this._adaptIESelection();		
			applied = this._adaptIEFormatAreaAndExec("bold");
		}
		if(!applied){
			applied = this.document.execCommand("bold", false, argument);
		}
		return applied;
	},

	_italicImpl: function(argument){
		// summary:
		//		This function implements an over-ride of the italic command.
		// argument:
		//		Not used, operates by selection.
		// tags:
		//		protected
		var applied = false;
		if(has("ie")){
			this._adaptIESelection();			
			applied = this._adaptIEFormatAreaAndExec("italic");
		}
		if(!applied){
			applied = this.document.execCommand("italic", false, argument);
		}
		return applied;
	},

	_underlineImpl: function(argument){
		// summary:
		//		This function implements an over-ride of the underline command.
		// argument:
		//		Not used, operates by selection.
		// tags:
		//		protected
		var applied = false;
		if(has("ie")){
			this._adaptIESelection();			
			applied = this._adaptIEFormatAreaAndExec("underline");
		}
		if(!applied){
			applied = this.document.execCommand("underline", false, argument);
		}
		return applied;
	},

	_strikethroughImpl: function(argument){
		// summary:
		//		This function implements an over-ride of the strikethrough command.
		// argument:
		//		Not used, operates by selection.
		// tags:
		//		protected
		var applied = false;
		if(has("ie")){
			this._adaptIESelection();			
			applied = this._adaptIEFormatAreaAndExec("strikethrough");
		}
		if(!applied){
			applied = this.document.execCommand("strikethrough", false, argument);
		}
		return applied;
	},

	_superscriptImpl: function(argument){
		// summary:
		//		This function implements an over-ride of the superscript command.
		// argument:
		//		Not used, operates by selection.
		// tags:
		//		protected
		var applied = false;
		if(has("ie")){
			this._adaptIESelection();			
			applied = this._adaptIEFormatAreaAndExec("superscript");
		}
		if(!applied){
			applied = this.document.execCommand("superscript", false, argument);
		}
		return applied;
	},

	_subscriptImpl: function(argument){
		// summary:
		//		This function implements an over-ride of the superscript command.
		// argument:
		//		Not used, operates by selection.
		// tags:
		//		protected
		var applied = false;
		if(has("ie")){
			this._adaptIESelection();			
			applied = this._adaptIEFormatAreaAndExec("subscript");
			
		}
		if(!applied){
			applied = this.document.execCommand("subscript", false, argument);
		}
		return applied;
	},
	
	_fontnameImpl: function(argument){
		// summary:
		//		This function implements the fontname command
		// argument:
		//		arguments to the exec command, if any.
		// tags:
		//		protected
		var isApplied;
		if(has("ie")){
			isApplied = this._handleTextColorOrProperties("fontname", argument);
		}
		if(!isApplied){
			isApplied = this.document.execCommand("fontname", false, argument);
		}
		return isApplied;
	},

	_fontsizeImpl: function(argument){
		// summary:
		//		This function implements the fontsize command
		// argument:
		//		arguments to the exec command, if any.
		// tags:
		//		protected
		var isApplied;
		if(has("ie")){
			isApplied = this._handleTextColorOrProperties("fontsize", argument);
		}
		if(!isApplied){
			isApplied = this.document.execCommand("fontsize", false, argument);
		}
		return isApplied;
	},
	
	_insertorderedlistImpl: function(argument){
		// summary:
		//		This function implements the insertorderedlist command
		// argument:
		//		arguments to the exec command, if any.
		// tags:
		//		protected
		var applied = false;
		if(has("ie")){
			applied = this._adaptIEList("insertorderedlist", argument);
		}
		if(!applied){
			applied = this.document.execCommand("insertorderedlist", false, argument);
		}
		return applied;
	},
	
	_insertunorderedlistImpl: function(argument){
		// summary:
		//		This function implements the insertunorderedlist command
		// argument:
		//		arguments to the exec command, if any.
		// tags:
		//		protected
		var applied = false;
		if(has("ie")){
			applied = this._adaptIEList("insertunorderedlist", argument);
		}
		if(!applied){
			applied = this.document.execCommand("insertunorderedlist", false, argument);
		}
		return applied;
	},
	
	getHeaderHeight: function(){
		// summary:
		//		A function for obtaining the height of the header node
		return this._getNodeChildrenHeight(this.header); // Number
	},

	getFooterHeight: function(){
		// summary:
		//		A function for obtaining the height of the footer node
		return this._getNodeChildrenHeight(this.footer); // Number
	},

	_getNodeChildrenHeight: function(node){
		// summary:
		//		An internal function for computing the cumulative height of all child nodes of 'node'
		// node:
		//		The node to process the children of;
		var h = 0;
		if(node && node.childNodes){
			// IE didn't compute it right when position was obtained on the node directly is some cases,
			// so we have to walk over all the children manually.
			var i;
			for(i = 0; i < node.childNodes.length; i++){
				var size = domGeometry.position(node.childNodes[i]);
				h += size.h;
			}
		}
		return h; // Number
	},

	_isNodeEmpty: function(node, startOffset){
		// summary:
		//		Function to test if a node is devoid of real content.
		// node:
		//		The node to check.
		// tags:
		//		private.
		if(node.nodeType === 1/*element*/){
			if(node.childNodes.length > 0){
				return this._isNodeEmpty(node.childNodes[0], startOffset);
	}
			return true;
		}else if(node.nodeType === 3/*text*/){
			return (node.nodeValue.substring(startOffset) === "");
		}
		return false;
	},

	_removeStartingRangeFromRange: function(node, range){
		// summary:
		//		Function to adjust selection range by removing the current
		//		start node.
		// node:
		//		The node to remove from the starting range.
		// range:
		//		The range to adapt.
		// tags:
		//		private
		if(node.nextSibling){
			range.setStart(node.nextSibling,0);
		}else{
			var parent = node.parentNode;
			while(parent && parent.nextSibling == null){
				//move up the tree until we find a parent that has another node, that node will be the next node
				parent = parent.parentNode;
			}
			if(parent){
				range.setStart(parent.nextSibling,0);
			}
		}
		return range;
	},

	_adaptIESelection: function(){
		// summary:
		//		Function to adapt the IE range by removing leading 'newlines'
		//		Needed to fix issue with bold/italics/underline not working if
		//		range included leading 'newlines'.
		//		In IE, if a user starts a selection at the very end of a line,
		//		then the native browser commands will fail to execute correctly.
		//		To work around the issue,  we can remove all empty nodes from
		//		the start of the range selection.
		var selection = rangeapi.getSelection(this.window);
		if(selection && selection.rangeCount && !selection.isCollapsed){
			var range = selection.getRangeAt(0);
			var firstNode = range.startContainer;
			var startOffset = range.startOffset;

			while(firstNode.nodeType === 3/*text*/ && startOffset >= firstNode.length && firstNode.nextSibling){
				//traverse the text nodes until we get to the one that is actually highlighted
				startOffset = startOffset - firstNode.length;
				firstNode = firstNode.nextSibling;
			}

			//Remove the starting ranges until the range does not start with an empty node.
			var lastNode=null;
			while(this._isNodeEmpty(firstNode, startOffset) && firstNode !== lastNode){
				lastNode =firstNode; //this will break the loop in case we can't find the next sibling
				range = this._removeStartingRangeFromRange(firstNode, range); //move the start container to the next node in the range
				firstNode = range.startContainer;
				startOffset = 0; //start at the beginning of the new starting range
			}
			selection.removeAllRanges();// this will work as long as users cannot select multiple ranges. I have not been able to do that in the editor.
			selection.addRange(range);
		}
	},
	
	_adaptIEFormatAreaAndExec: function(command){
		// summary:
		//		Function to handle IE's quirkiness regarding how it handles
		//		format commands on a word.  This involves a lit of node splitting
		//		and format cloning.
		// command:
		//		The format command, needed to check if the desired
		//		command is true or not.
		var selection = rangeapi.getSelection(this.window);
		var doc = this.document;
		var rs, ret, range, txt, startNode, endNode, breaker, sNode;
		if(command && selection && selection.isCollapsed){
			var isApplied = this.queryCommandValue(command);
			if(isApplied){
				
				// We have to split backwards until we hit the format
				var nNames = this._tagNamesForCommand(command);
				range = selection.getRangeAt(0);
				var fs = range.startContainer;
				if(fs.nodeType === 3){
					var offset = range.endOffset;
					if(fs.length < offset){
						//We are not looking from the right node, try to locate the correct one
						ret = this._adjustNodeAndOffset(rs, offset);
						fs = ret.node;
						offset = ret.offset;
					}
				}									
				var topNode;
				while(fs && fs !== this.editNode){
					// We have to walk back and see if this is still a format or not.
					// Hm, how do I do this?
					var tName = fs.tagName? fs.tagName.toLowerCase() : "";
					if(array.indexOf(nNames, tName) > -1){
						topNode = fs;
						break;
					}
					fs = fs.parentNode;
				}

				// Okay, we have a stopping place, time to split things apart.
				if(topNode){
					// Okay, we know how far we have to split backwards, so we have to split now.
					rs = range.startContainer;
					var newblock = doc.createElement(topNode.tagName);
					domConstruct.place(newblock, topNode, "after");
					if(rs && rs.nodeType === 3){
						// Text node, we have to split it.
						var nodeToMove, tNode;
						var endOffset = range.endOffset;
						if(rs.length < endOffset){
							//We are not splitting the right node, try to locate the correct one
							ret = this._adjustNodeAndOffset(rs, endOffset);
							rs = ret.node;
							endOffset = ret.offset;
						}
		
						txt = rs.nodeValue;
						startNode = doc.createTextNode(txt.substring(0, endOffset));
						var endText = txt.substring(endOffset, txt.length);
						if(endText){
							endNode = doc.createTextNode(endText);
						}
						// Place the split, then remove original nodes.
						domConstruct.place(startNode, rs, "before");
						if(endNode){
							breaker = doc.createElement("span");
							breaker.className = "ieFormatBreakerSpan";
							domConstruct.place(breaker, rs, "after");
							domConstruct.place(endNode, breaker, "after");
							endNode = breaker;
						}
						domConstruct.destroy(rs);
						
						// Okay, we split the text.  Now we need to see if we're
						// parented to the block element we're splitting and if
						// not, we have to split all the way up.  Ugh.
						var parentC = startNode.parentNode;
						var tagList = [];
						var tagData;
						while(parentC !== topNode){
							var tg = parentC.tagName;
							tagData = {tagName: tg};
							tagList.push(tagData);
														
							var newTg = doc.createElement(tg);
							// Clone over any 'style' data.
							if(parentC.style){
								if(newTg.style){
									if(parentC.style.cssText){
										newTg.style.cssText = parentC.style.cssText;
										tagData.cssText = parentC.style.cssText;
									}
								}
							}
							// If font also need to clone over any font data.
							if(parentC.tagName === "FONT"){
								if(parentC.color){
									newTg.color = parentC.color;
									tagData.color = parentC.color;
								}
								if(parentC.face){
									newTg.face = parentC.face;
									tagData.face = parentC.face;
								}
								if(parentC.size){  // this check was necessary on IE
									newTg.size = parentC.size;
									tagData.size = parentC.size;
								}
							}
							if(parentC.className){
								newTg.className = parentC.className;
								tagData.className = parentC.className;
							}
							
							// Now move end node and every sibling 
							// after it over into the new tag.
							if(endNode){
								nodeToMove = endNode;
								while(nodeToMove){
									tNode = nodeToMove.nextSibling;
									newTg.appendChild(nodeToMove);
									nodeToMove = tNode;
								}
							}
							if(newTg.tagName == parentC.tagName){
								breaker = doc.createElement("span");
								breaker.className = "ieFormatBreakerSpan";
								domConstruct.place(breaker, parentC, "after");
								domConstruct.place(newTg, breaker, "after");
							}else{
								domConstruct.place(newTg, parentC, "after");
							}
							startNode = parentC;
							endNode = newTg;
							parentC = parentC.parentNode;
						}

						// Lastly, move the split out all the split tags 
						// to the new block as they should now be split properly.
						if(endNode){
							nodeToMove = endNode;
							if(nodeToMove.nodeType === 1 || (nodeToMove.nodeType === 3 && nodeToMove.nodeValue)){
								// Non-blank text and non-text nodes need to clear out that blank space
								// before moving the contents.
								newblock.innerHTML = "";
							}
							while(nodeToMove){
								tNode = nodeToMove.nextSibling;
								newblock.appendChild(nodeToMove);
								nodeToMove = tNode;
							}
						}
						
						// We had intermediate tags, we have to now recreate them inbetween the split
						// and restore what styles, classnames, etc, we can.
						var newrange;
						if(tagList.length){
							tagData = tagList.pop();
							var newContTag = doc.createElement(tagData.tagName);
							if(tagData.cssText && newContTag.style){
								newContTag.style.cssText = tagData.cssText;
							}
							if(tagData.className){
								newContTag.className = tagData.className;
							}
							if(tagData.tagName === "FONT"){
								if(tagData.color){
									newContTag.color = tagData.color;
								}
								if(tagData.face){
									newContTag.face = tagData.face;
								}
								if(tagData.size){ 
									newContTag.size = tagData.size;
								}
							}								
							domConstruct.place(newContTag, newblock, "before");
							while(tagList.length){
								tagData = tagList.pop();
								var newTgNode = doc.createElement(tagData.tagName);
								if(tagData.cssText && newTgNode.style){
									newTgNode.style.cssText = tagData.cssText;
								}
								if(tagData.className){
									newTgNode.className = tagData.className;
								}
								if(tagData.tagName === "FONT"){
									if(tagData.color){
										newTgNode.color = tagData.color;
									}
									if(tagData.face){
										newTgNode.face = tagData.face;
									}
									if(tagData.size){ 
										newTgNode.size = tagData.size;
									}
								}	
								newContTag.appendChild(newTgNode);
								newContTag = newTgNode;
							}							
							
							// Okay, everything is theoretically split apart and removed from the content
							// so insert the dummy text to select, select it, then
							// clear to position cursor.
							sNode = doc.createTextNode(".");
							breaker.appendChild(sNode);
							newContTag.appendChild(sNode);
							newrange = rangeapi.create(this.window);
							newrange.setStart(sNode, 0);
							newrange.setEnd(sNode, sNode.length);
							selection.removeAllRanges();
							selection.addRange(newrange);
							this._sCall("collapse", [false]);
							sNode.parentNode.innerHTML = "";
						}else{
							// No extra tags, so we have to insert a breaker point and rely
							// on filters to remove it later.
							breaker = doc.createElement("span");
							breaker.className="ieFormatBreakerSpan";
							sNode = doc.createTextNode(".");
							breaker.appendChild(sNode);
							domConstruct.place(breaker, newblock, "before");
							newrange = rangeapi.create(this.window);
							newrange.setStart(sNode, 0);
							newrange.setEnd(sNode, sNode.length);
							selection.removeAllRanges();
							selection.addRange(newrange);
							this._sCall("collapse", [false]);
							sNode.parentNode.innerHTML = "";
						}
						if(!newblock.firstChild){
							// Empty, we don't need it.  Split was at end or similar
							// So, remove it.
							domConstruct.destroy(newblock);
						}					
						return true;
					}
				}
				return false;
			}else{
				range = selection.getRangeAt(0);
				rs = range.startContainer;
				if(rs && rs.nodeType === 3){
					// Text node, we have to split it.
					var offset = range.startOffset;
					if(rs.length < offset){
						//We are not splitting the right node, try to locate the correct one
						ret = this._adjustNodeAndOffset(rs, offset);
						rs = ret.node;
						offset = ret.offset;
					}
					txt = rs.nodeValue;
					startNode = doc.createTextNode(txt.substring(0, offset));
					var endText = txt.substring(offset);
					if(endText !== ""){
						endNode = doc.createTextNode(txt.substring(offset));
					}
					// Create a space, we'll select and bold it, so
					// the whole word doesn't get bolded
					breaker = doc.createElement("span");
					sNode = doc.createTextNode(".");
					breaker.appendChild(sNode);
					if(startNode.length){
						domConstruct.place(startNode, rs, "after");
					}else{
						startNode = rs;
					}
					domConstruct.place(breaker, startNode, "after");
					if(endNode){
						domConstruct.place(endNode, breaker, "after");
					}
					domConstruct.destroy(rs);
					var newrange = rangeapi.create(this.window);
					newrange.setStart(sNode, 0);
					newrange.setEnd(sNode, sNode.length);
					selection.removeAllRanges();
					selection.addRange(newrange);
					doc.execCommand(command);
					domConstruct.place(breaker.firstChild, breaker, "before");
					domConstruct.destroy(breaker);
					newrange.setStart(sNode, 0);
					newrange.setEnd(sNode, sNode.length);
					selection.removeAllRanges();
					selection.addRange(newrange);
					this._sCall("collapse", [false]);
					sNode.parentNode.innerHTML = "";
					return true;
				}
			}
		}else{
			return false;
		}
	},
	
	_adaptIEList: function(command /*===== , argument =====*/){
		// summary:
		//		This function handles normalizing the IE list behavior as 
		//		much as possible.
		// command:
		//		The list command to execute.
		// argument:
		//		Any additional argument.
		// tags:
		//		private
		var selection = rangeapi.getSelection(this.window);
		if(selection.isCollapsed){
			// In the case of no selection, lets commonize the behavior and
			// make sure that it indents if needed.
			if(selection.rangeCount && !this.queryCommandValue(command)){
				var range = selection.getRangeAt(0);
				var sc = range.startContainer;
				if(sc && sc.nodeType == 3){
					// text node.  Lets see if there is a node before it that isn't
					// some sort of breaker.
					if(!range.startOffset){
						// We're at the beginning of a text area.  It may have been br split
						// Who knows?  In any event, we must create the list manually
						// or IE may shove too much into the list element.  It seems to
						// grab content before the text node too if it's br split.
						// Why can't IE work like everyone else?

						// Create a space, we'll select and bold it, so
						// the whole word doesn't get bolded
						var lType = "ul";
						if(command === "insertorderedlist"){
							lType = "ol";
						}
						var list = this.document.createElement(lType);
						var li = domConstruct.create("li", null, list);
						domConstruct.place(list, sc, "before");
						// Move in the text node as part of the li.
						li.appendChild(sc);
						// We need a br after it or the enter key handler
						// sometimes throws errors.
						domConstruct.create("br", null, list, "after");
						// Okay, now lets move our cursor to the beginning.
						var newrange = rangeapi.create(this.window);
						newrange.setStart(sc, 0);
						newrange.setEnd(sc, sc.length);
						selection.removeAllRanges();
						selection.addRange(newrange);
						this._sCall("collapse", [true]);
						return true;
					}
				}
			}
		}
		return false;
	},
	
	_handleTextColorOrProperties: function(command, argument){
		// summary:
		//		This function handles appplying text color as best it is 
		//		able to do so when the selection is collapsed, making the
		//		behavior cross-browser consistent. It also handles the name
		//		and size for IE.
		// command:
		//		The command.
		// argument:
		//		Any additional arguments.
		// tags:
		//		private
		var selection = rangeapi.getSelection(this.window);
		var doc = this.document;
		var rs, ret, range, txt, startNode, endNode, breaker, sNode;
		argument = argument || null;
		if(command && selection && selection.isCollapsed){
			if(selection.rangeCount){
				range = selection.getRangeAt(0);
				rs = range.startContainer;
				if(rs && rs.nodeType === 3){
					// Text node, we have to split it.
					var offset = range.startOffset;
					if(rs.length < offset){
						//We are not splitting the right node, try to locate the correct one
						ret = this._adjustNodeAndOffset(rs, offset);
						rs = ret.node;
						offset = ret.offset;
					}
					txt = rs.nodeValue;
					startNode = doc.createTextNode(txt.substring(0, offset));
					var endText = txt.substring(offset);
					if(endText !== ""){
						endNode = doc.createTextNode(txt.substring(offset));
					}
					// Create a space, we'll select and bold it, so
					// the whole word doesn't get bolded
					breaker = doc.createElement("span");
					sNode = doc.createTextNode(".");
					breaker.appendChild(sNode);
					// Create a junk node to avoid it trying to style the breaker.
					// This will get destroyed later.
					var extraSpan = doc.createElement("span");
					breaker.appendChild(extraSpan);
					if(startNode.length){
						domConstruct.place(startNode, rs, "after");
					}else{
						startNode = rs;
					}
					domConstruct.place(breaker, startNode, "after");
					if(endNode){
						domConstruct.place(endNode, breaker, "after");
					}
					domConstruct.destroy(rs);
					var newrange = rangeapi.create(this.window);
					newrange.setStart(sNode, 0);
					newrange.setEnd(sNode, sNode.length);
					selection.removeAllRanges();
					selection.addRange(newrange);
					if(has("webkit")){
						// WebKit is frustrating with positioning the cursor.
						// It stinks to have a selected space, but there really
						// isn't much choice here.
						var style = "color";
						if(command === "hilitecolor" || command === "backcolor"){
							style = "backgroundColor";
						}
						domStyle.set(breaker, style, argument);
						this._sCall("remove", []);
						domConstruct.destroy(extraSpan);
						breaker.innerHTML = "&#160;";	// &nbsp;
						this._sCall("selectElement", [breaker]);
						this.focus();
					}else{
						this.execCommand(command, argument);
						domConstruct.place(breaker.firstChild, breaker, "before");
						domConstruct.destroy(breaker);
						newrange.setStart(sNode, 0);
						newrange.setEnd(sNode, sNode.length);
						selection.removeAllRanges();
						selection.addRange(newrange);
						this._sCall("collapse", [false]);
						sNode.parentNode.removeChild(sNode);
					}
					return true;
				}
			}				
		}
		return false;
	},
	
	_adjustNodeAndOffset: function(/*DomNode*/node, /*Int*/offset){
		// summary:
		//		In the case there are multiple text nodes in a row the offset may not be within the node.  
		//		If the offset is larger than the node length, it will attempt to find
		//		the next text sibling until it locates the text node in which the offset refers to
		// node:
		//		The node to check.
		// offset:
		//		The position to find within the text node
		// tags:
		//		private.
		while(node.length < offset && node.nextSibling && node.nextSibling.nodeType === 3){
			//Adjust the offset and node in the case of multiple text nodes in a row
			offset = offset - node.length;
			node = node.nextSibling;
		}
		return {"node": node, "offset": offset};
	},
	
	_tagNamesForCommand: function(command){
		// summary:
		//		Function to return the tab names that are associated
		//		with a particular style.
		// command: String
		//		The command to return tags for.
		// tags:
		//		private
		if(command === "bold"){
			return ["b", "strong"];
		}else if(command === "italic"){
			return ["i","em"];
		}else if(command === "strikethrough"){
			return ["s", "strike"];
		}else if(command === "superscript"){
			return ["sup"];
		}else if(command === "subscript"){
			return ["sub"];
		}else if(command === "underline"){
			return ["u"];
		}	
		return [];
	},

	_stripBreakerNodes: function(/*DOMNode*/ node){
		// summary:
		//		Function for stripping out the breaker spans inserted by the formatting command.
		//		Registered as a filter for IE, handles the breaker spans needed to fix up
		//		How bold/italic/etc, work when selection is collapsed (single cursor).
		if(!this.isLoaded){ return; } // this method requires init to be complete
		query(".ieFormatBreakerSpan", node).forEach(function(b){
			while(b.firstChild){
				domConstruct.place(b.firstChild, b, "before");
			}
			domConstruct.destroy(b);
		});
		return node;
	}
});

return RichText;

});

},
'dijit/_editor/selection':function(){
define("dijit/_editor/selection", [
	"dojo/dom", // dom.byId
	"dojo/_base/lang",
	"dojo/sniff", // has("ie") has("opera")
	"dojo/_base/window", // win.body win.doc win.doc.createElement win.doc.selection win.doc.selection.createRange win.doc.selection.type.toLowerCase win.global win.global.getSelection
	"../main"		// for exporting symbols to dijit._editor.selection (TODO: remove in 2.0)
], function(dom, lang, has, win, dijit){

// module:
//		dijit/_editor/selection
// summary:
//		Text selection API

// FIXME:
//		all of these methods branch internally for IE. This is probably
//		sub-optimal in terms of runtime performance. We should investigate the
//		size difference for differentiating at definition time.

var selection = {
	getType: function(){
		// summary:
		//		Get the selection type (like win.doc.select.type in IE).
		if(win.doc.getSelection){
			// W3C path
			var stype = "text";

			// Check if the actual selection is a CONTROL (IMG, TABLE, HR, etc...).
			var oSel;
			try{
				oSel = win.global.getSelection();
			}catch(e){ /*squelch*/ }

			if(oSel && oSel.rangeCount == 1){
				var oRange = oSel.getRangeAt(0);
				if(	(oRange.startContainer == oRange.endContainer) &&
					((oRange.endOffset - oRange.startOffset) == 1) &&
					(oRange.startContainer.nodeType != 3 /* text node*/)
				){
					stype = "control";
				}
			}
			return stype; //String
		}else{
			// IE6-8
			return win.doc.selection.type.toLowerCase();
		}
	},

	getSelectedText: function(){
		// summary:
		//		Return the text (no html tags) included in the current selection or null if no text is selected
		if(win.doc.getSelection){
			// W3C path
			var selection = win.global.getSelection();
			return selection ? selection.toString() : ""; //String
		}else{
			// IE6-8
			if(dijit._editor.selection.getType() == 'control'){
				return null;
			}
			return win.doc.selection.createRange().text;
		}
	},

	getSelectedHtml: function(){
		// summary:
		//		Return the html text of the current selection or null if unavailable
		if(win.doc.getSelection){
			// W3C path
			var selection = win.global.getSelection();
			if(selection && selection.rangeCount){
				var i;
				var html = "";
				for(i = 0; i < selection.rangeCount; i++){
					//Handle selections spanning ranges, such as Opera
					var frag = selection.getRangeAt(i).cloneContents();
					var div = win.doc.createElement("div");
					div.appendChild(frag);
					html += div.innerHTML;
				}
				return html; //String
			}
			return null;
		}else{
			// IE6-8
			if(dijit._editor.selection.getType() == 'control'){
				return null;
			}
			return win.doc.selection.createRange().htmlText;
		}
	},

	getSelectedElement: function(){
		// summary:
		//		Retrieves the selected element (if any), just in the case that
		//		a single element (object like and image or a table) is
		//		selected.
		if(dijit._editor.selection.getType() == "control"){
			if(win.doc.getSelection){
				// W3C path
				var selection = win.global.getSelection();
				return selection.anchorNode.childNodes[ selection.anchorOffset ];
			}else{
				// IE6-8
				var range = win.doc.selection.createRange();
				if(range && range.item){
					return win.doc.selection.createRange().item(0);
				}
			}
		}
		return null;
	},

	getParentElement: function(){
		// summary:
		//		Get the parent element of the current selection
		if(dijit._editor.selection.getType() == "control"){
			var p = this.getSelectedElement();
			if(p){ return p.parentNode; }
		}else{
			if(win.doc.getSelection){
				var selection = win.global.getSelection();
				if(selection){
					var node = selection.anchorNode;
					while(node && (node.nodeType != 1)){ // not an element
						node = node.parentNode;
					}
					return node;
				}
			}else{
				var r = win.doc.selection.createRange();
				r.collapse(true);
				return r.parentElement();
			}
		}
		return null;
	},

	hasAncestorElement: function(/*String*/ tagName /* ... */){
		// summary:
		//		Check whether current selection has a  parent element which is
		//		of type tagName (or one of the other specified tagName)
		// tagName: String
		//		The tag name to determine if it has an ancestor of.
		return this.getAncestorElement.apply(this, arguments) != null; //Boolean
	},

	getAncestorElement: function(/*String*/ tagName /* ... */){
		// summary:
		//		Return the parent element of the current selection which is of
		//		type tagName (or one of the other specified tagName)
		// tagName: String
		//		The tag name to determine if it has an ancestor of.
		var node = this.getSelectedElement() || this.getParentElement();
		return this.getParentOfType(node, arguments); //DOMNode
	},

	isTag: function(/*DomNode*/ node, /*String[]*/ tags){
		// summary:
		//		Function to determine if a node is one of an array of tags.
		// node:
		//		The node to inspect.
		// tags:
		//		An array of tag name strings to check to see if the node matches.
		if(node && node.tagName){
			var _nlc = node.tagName.toLowerCase();
			for(var i=0; i<tags.length; i++){
				var _tlc = String(tags[i]).toLowerCase();
				if(_nlc == _tlc){
					return _tlc; // String
				}
			}
		}
		return "";
	},

	getParentOfType: function(/*DomNode*/ node, /*String[]*/ tags){
		// summary:
		//		Function to locate a parent node that matches one of a set of tags
		// node:
		//		The node to inspect.
		// tags:
		//		An array of tag name strings to check to see if the node matches.
		while(node){
			if(this.isTag(node, tags).length){
				return node; // DOMNode
			}
			node = node.parentNode;
		}
		return null;
	},

	collapse: function(/*Boolean*/ beginning){
		// summary:
		//		Function to collapse (clear), the current selection
		// beginning: Boolean
		//		Indicates whether to collapse the cursor to the beginning of the selection or end.
		if(win.doc.getSelection){
			// W3C path
			var selection = win.global.getSelection();
			if(selection.removeAllRanges){ // Mozilla
				if(beginning){
					selection.collapseToStart();
				}else{
					selection.collapseToEnd();
				}
			}else{ // Safari
				// pulled from WebCore/ecma/kjs_window.cpp, line 2536
				selection.collapse(beginning);
			}
		}else{
			// IE6-8
			var range = win.doc.selection.createRange();
			range.collapse(beginning);
			range.select();
		}
	},

	remove: function(){
		// summary:
		//		Function to delete the currently selected content from the document.
		var sel = win.doc.selection;
		if(win.doc.getSelection){
			// W3C path
			sel = win.global.getSelection();
			sel.deleteFromDocument();
			return sel; //Selection
		}else{
			// IE6-8
			if(sel.type.toLowerCase() != "none"){
				sel.clear();
			}
			return sel; //Selection
		}
	},

	selectElementChildren: function(/*DomNode*/ element, /*Boolean?*/ nochangefocus){
		// summary:
		//		clear previous selection and select the content of the node
		//		(excluding the node itself)
		// element: DOMNode
		//		The element you wish to select the children content of.
		// nochangefocus: Boolean
		//		Indicates if the focus should change or not.
		var doc = win.doc;
		var range;
		element = dom.byId(element);
		if(win.doc.getSelection){
			// W3C
			var selection = win.global.getSelection();
			if(has("opera")){
				//Opera's selectAllChildren doesn't seem to work right
				//against <body> nodes and possibly others ... so
				//we use the W3C range API
				if(selection.rangeCount){
					range = selection.getRangeAt(0);
				}else{
					range = doc.createRange();
				}
				range.setStart(element, 0);
				range.setEnd(element,(element.nodeType == 3) ? element.length : element.childNodes.length);
				selection.addRange(range);
			}else{
				selection.selectAllChildren(element);
			}
		}else{
			// IE6-8
			range = element.ownerDocument.body.createTextRange();
			range.moveToElementText(element);
			if(!nochangefocus){
				try{
					range.select(); // IE throws an exception here if the widget is hidden.  See #5439
				}catch(e){ /* squelch */}
			}
		}
	},

	selectElement: function(/*DomNode*/ element, /*Boolean?*/ nochangefocus){
		// summary:
		//		clear previous selection and select element (including all its children)
		// element: DOMNode
		//		The element to select.
		// nochangefocus: Boolean
		//		Boolean indicating if the focus should be changed.  IE only.
		var range;
		element = dom.byId(element);	// TODO: remove for 2.0 or sooner, spec listed above doesn't allow for string
		var doc = element.ownerDocument;
		var global = win.global;	// TODO: use winUtils.get(doc)?
		if(doc.getSelection){
			// W3C path
			var selection = global.getSelection();
			range = doc.createRange();
			if(selection.removeAllRanges){ // Mozilla
				// FIXME: does this work on Safari?
				if(has("opera")){
					//Opera works if you use the current range on
					//the selection if present.
					if(selection.getRangeAt(0)){
						range = selection.getRangeAt(0);
					}
				}
				range.selectNode(element);
				selection.removeAllRanges();
				selection.addRange(range);
			}
		}else{
			// IE6-8
			try{
				var tg = element.tagName ? element.tagName.toLowerCase() : "";
				if(tg === "img" || tg === "table"){
					range = win.body(doc).createControlRange();
				}else{
					range = win.body(doc).createRange();
				}
				range.addElement(element);
				if(!nochangefocus){
					range.select();
				}
			}catch(e){
				this.selectElementChildren(element, nochangefocus);
			}
		}
	},

	inSelection: function(node){
		// summary:
		//		This function determines if 'node' is
		//		in the current selection.
		// tags:
		//		public
		if(node){
			var newRange;
			var doc = win.doc;
			var range;

			if(win.doc.getSelection){
				// WC3
				var sel = win.global.getSelection();
				if(sel && sel.rangeCount > 0){
					range = sel.getRangeAt(0);
				}
				if(range && range.compareBoundaryPoints && doc.createRange){
					try{
						newRange = doc.createRange();
						newRange.setStart(node, 0);
						if(range.compareBoundaryPoints(range.START_TO_END, newRange) === 1){
							return true;
						}
					}catch(e){ /* squelch */}
				}
			}else{
				// IE6-8, so we can't use the range object as the pseudo
				// range doesn't implement the boundary checking, we have to
				// use IE specific crud.
				range = doc.selection.createRange();
				try{
					newRange = node.ownerDocument.body.createControlRange();
					if(newRange){
						newRange.addElement(node);
					}
				}catch(e1){
					try{
						newRange = node.ownerDocument.body.createTextRange();
						newRange.moveToElementText(node);
					}catch(e2){/* squelch */}
				}
				if(range && newRange){
					// We can finally compare similar to W3C
					if(range.compareEndPoints("EndToStart", newRange) === 1){
						return true;
					}
				}
			}
		}
		return false; // Boolean
	}
};


lang.setObject("dijit._editor.selection", selection);

return selection;
});

},
'dijit/_editor/range':function(){
define("dijit/_editor/range", [
	"dojo/_base/array", // array.every
	"dojo/_base/declare", // declare
	"dojo/_base/lang", // lang.isArray
	"dojo/_base/window", // win.doc   TODO: remove in 2.0
	"../main"	// for exporting symbols to dijit, TODO: remove in 2.0
], function(array, declare, lang, win, dijit){

// module:
//		dijit/_editor/range
// summary:
//		W3C range API


dijit.range={};

dijit.range.getIndex = function(/*DomNode*/node, /*DomNode*/parent){
//	dojo.profile.start("dijit.range.getIndex");
	var ret = [], retR = [];
	var onode = node;

	var pnode, n;
	while(node != parent){
		var i = 0;
		pnode = node.parentNode;
		while((n = pnode.childNodes[i++])){
			if(n === node){
				--i;
				break;
			}
		}
		//if(i>=pnode.childNodes.length){
			//dojo.debug("Error finding index of a node in dijit/range.getIndex()");
		//}
		ret.unshift(i);
		retR.unshift(i - pnode.childNodes.length);
		node = pnode;
	}

	//normalized() can not be called so often to prevent
	//invalidating selection/range, so we have to detect
	//here that any text nodes in a row
	if(ret.length > 0 && onode.nodeType == 3){
		n = onode.previousSibling;
		while(n && n.nodeType == 3){
			ret[ret.length - 1]--;
			n = n.previousSibling;
		}
		n = onode.nextSibling;
		while(n && n.nodeType == 3){
			retR[retR.length - 1]++;
			n = n.nextSibling;
		}
	}
//	dojo.profile.end("dijit/range.getIndex()");
	return {o: ret, r:retR};
};

dijit.range.getNode = function(/*Array*/index, /*DomNode*/parent){
	if(!lang.isArray(index) || index.length == 0){
		return parent;
	}
	var node = parent;
//	if(!node)debugger
	array.every(index, function(i){
		if(i >= 0 && i < node.childNodes.length){
			node = node.childNodes[i];
		}else{
			node = null;
			//console.debug('Error: can not find node with index',index,'under parent node',parent );
			return false; //terminate array.every
		}
		return true; //carry on the every loop
	});

	return node;
};

dijit.range.getCommonAncestor = function(n1, n2, root){
	root = root || n1.ownerDocument.body;
	var getAncestors = function(n){
		var as = [];
		while(n){
			as.unshift(n);
			if(n !== root){
				n = n.parentNode;
			}else{
				break;
			}
		}
		return as;
	};
	var n1as = getAncestors(n1);
	var n2as = getAncestors(n2);

	var m = Math.min(n1as.length, n2as.length);
	var com = n1as[0]; //at least, one element should be in the array: the root (BODY by default)
	for(var i = 1; i < m; i++){
		if(n1as[i] === n2as[i]){
			com = n1as[i]
		}else{
			break;
		}
	}
	return com;
};

dijit.range.getAncestor = function(/*DomNode*/node, /*RegEx?*/regex, /*DomNode?*/root){
	root = root || node.ownerDocument.body;
	while(node && node !== root){
		var name = node.nodeName.toUpperCase();
		if(regex.test(name)){
			return node;
		}

		node = node.parentNode;
	}
	return null;
};

dijit.range.BlockTagNames = /^(?:P|DIV|H1|H2|H3|H4|H5|H6|ADDRESS|PRE|OL|UL|LI|DT|DE)$/;
dijit.range.getBlockAncestor = function(/*DomNode*/node, /*RegEx?*/regex, /*DomNode?*/root){
	root = root || node.ownerDocument.body;
	regex = regex || dijit.range.BlockTagNames;
	var block = null, blockContainer;
	while(node && node !== root){
		var name = node.nodeName.toUpperCase();
		if(!block && regex.test(name)){
			block = node;
		}
		if(!blockContainer && (/^(?:BODY|TD|TH|CAPTION)$/).test(name)){
			blockContainer = node;
		}

		node = node.parentNode;
	}
	return {blockNode:block, blockContainer:blockContainer || node.ownerDocument.body};
};

dijit.range.atBeginningOfContainer = function(/*DomNode*/container, /*DomNode*/node, /*Int*/offset){
	var atBeginning = false;
	var offsetAtBeginning = (offset == 0);
	if(!offsetAtBeginning && node.nodeType == 3){ //if this is a text node, check whether the left part is all space
		if(/^[\s\xA0]+$/.test(node.nodeValue.substr(0, offset))){
			offsetAtBeginning = true;
		}
	}
	if(offsetAtBeginning){
		var cnode = node;
		atBeginning = true;
		while(cnode && cnode !== container){
			if(cnode.previousSibling){
				atBeginning = false;
				break;
			}
			cnode = cnode.parentNode;
		}
	}
	return atBeginning;
};

dijit.range.atEndOfContainer = function(/*DomNode*/container, /*DomNode*/node, /*Int*/offset){
	var atEnd = false;
	var offsetAtEnd = (offset == (node.length || node.childNodes.length));
	if(!offsetAtEnd && node.nodeType == 3){ //if this is a text node, check whether the right part is all space
		if(/^[\s\xA0]+$/.test(node.nodeValue.substr(offset))){
			offsetAtEnd = true;
		}
	}
	if(offsetAtEnd){
		var cnode = node;
		atEnd = true;
		while(cnode && cnode !== container){
			if(cnode.nextSibling){
				atEnd = false;
				break;
			}
			cnode = cnode.parentNode;
		}
	}
	return atEnd;
};

dijit.range.adjacentNoneTextNode = function(startnode, next){
	var node = startnode;
	var len = (0 - startnode.length) || 0;
	var prop = next ? 'nextSibling' : 'previousSibling';
	while(node){
		if(node.nodeType != 3){
			break;
		}
		len += node.length;
		node = node[prop];
	}
	return [node,len];
};

dijit.range.create = function(/*Window?*/ win){	// TODO: for 2.0, replace optional window param w/mandatory window or document param
	win = win || window;
	if(win.getSelection){
		return win.document.createRange();
	}else{//IE
		return new dijit.range.W3CRange();
	}
};

dijit.range.getSelection = function(/*Window*/ window, /*Boolean?*/ ignoreUpdate){
	if(window.getSelection){
		return window.getSelection();
	}else{//IE
		var s = new dijit.range.ie.selection(window);
		if(!ignoreUpdate){
			s._getCurrentSelection();
		}
		return s;
	}
};

// TODO: convert to has() test?   But remember IE9 issues with quirks vs. standards in main frame vs. iframe.
if(!window.getSelection){
	dijit.range.ie = {
		cachedSelection: {},
		selection: function(window){
			this._ranges = [];
			this.addRange = function(r, /*boolean*/ internal){
				this._ranges.push(r);
				if(!internal){
					r._select();
				}
				this.rangeCount = this._ranges.length;
			};
			this.removeAllRanges = function(){
				//don't detach, the range may be used later
//				for(var i=0;i<this._ranges.length;i++){
//					this._ranges[i].detach();
//				}
				this._ranges = [];
				this.rangeCount = 0;
			};
			var _initCurrentRange = function(){
				var r = window.document.selection.createRange();
				var type = window.document.selection.type.toUpperCase();
				if(type == "CONTROL"){
					//TODO: multiple range selection(?)
					return new dijit.range.W3CRange(dijit.range.ie.decomposeControlRange(r));
				}else{
					return new dijit.range.W3CRange(dijit.range.ie.decomposeTextRange(r));
				}
			};
			this.getRangeAt = function(i){
				return this._ranges[i];
			};
			this._getCurrentSelection = function(){
				this.removeAllRanges();
				var r = _initCurrentRange();
				if(r){
					this.addRange(r, true);
					this.isCollapsed = r.collapsed;
				}else{
					this.isCollapsed = true;
				}
			};
		},
		decomposeControlRange: function(range){
			var firstnode = range.item(0), lastnode = range.item(range.length - 1);
			var startContainer = firstnode.parentNode, endContainer = lastnode.parentNode;
			var startOffset = dijit.range.getIndex(firstnode, startContainer).o[0];
			var endOffset = dijit.range.getIndex(lastnode, endContainer).o[0] + 1;
			return [startContainer, startOffset,endContainer, endOffset];
		},
		getEndPoint: function(range, end){
			var atmrange = range.duplicate();
			atmrange.collapse(!end);
			var cmpstr = 'EndTo' + (end ? 'End' : 'Start');
			var parentNode = atmrange.parentElement();

			var startnode, startOffset, lastNode;
			if(parentNode.childNodes.length > 0){
				array.every(parentNode.childNodes, function(node, i){
					var calOffset;
					if(node.nodeType != 3){
						atmrange.moveToElementText(node);

						if(atmrange.compareEndPoints(cmpstr, range) > 0){
							//startnode = node.previousSibling;
							if(lastNode && lastNode.nodeType == 3){
								//where shall we put the start? in the text node or after?
								startnode = lastNode;
								calOffset = true;
							}else{
								startnode = parentNode;
								startOffset = i;
								return false;
							}
						}else{
							if(i == parentNode.childNodes.length - 1){
								startnode = parentNode;
								startOffset = parentNode.childNodes.length;
								return false;
							}
						}
					}else{
						if(i == parentNode.childNodes.length - 1){//at the end of this node
							startnode = node;
							calOffset = true;
						}
					}
					//			try{
					if(calOffset && startnode){
						var prevnode = dijit.range.adjacentNoneTextNode(startnode)[0];
						if(prevnode){
							startnode = prevnode.nextSibling;
						}else{
							startnode = parentNode.firstChild; //firstChild must be a text node
						}
						var prevnodeobj = dijit.range.adjacentNoneTextNode(startnode);
						prevnode = prevnodeobj[0];
						var lenoffset = prevnodeobj[1];
						if(prevnode){
							atmrange.moveToElementText(prevnode);
							atmrange.collapse(false);
						}else{
							atmrange.moveToElementText(parentNode);
						}
						atmrange.setEndPoint(cmpstr, range);
						startOffset = atmrange.text.length - lenoffset;

						return false;
					}
					//			}catch(e){ debugger }
					lastNode = node;
					return true;
				});
			}else{
				startnode = parentNode;
				startOffset = 0;
			}

			//if at the end of startnode and we are dealing with start container, then
			//move the startnode to nextSibling if it is a text node
			//TODO: do this for end container?
			if(!end && startnode.nodeType == 1 && startOffset == startnode.childNodes.length){
				var nextnode = startnode.nextSibling;
				if(nextnode && nextnode.nodeType == 3){
					startnode = nextnode;
					startOffset = 0;
				}
			}
			return [startnode, startOffset];
		},
		setEndPoint: function(range, container, offset){
			//text node
			var atmrange = range.duplicate(), node, len;
			if(container.nodeType != 3){ //normal node
				if(offset > 0){
					node = container.childNodes[offset - 1];
					if(node){
						if(node.nodeType == 3){
							container = node;
							offset = node.length;
							//pass through
						}else{
							if(node.nextSibling && node.nextSibling.nodeType == 3){
								container = node.nextSibling;
								offset = 0;
								//pass through
							}else{
								atmrange.moveToElementText(node.nextSibling ? node : container);
								var parent = node.parentNode;
								var tempNode = parent.insertBefore(node.ownerDocument.createTextNode(' '), node.nextSibling);
								atmrange.collapse(false);
								parent.removeChild(tempNode);
							}
						}
					}
				}else{
					atmrange.moveToElementText(container);
					atmrange.collapse(true);
				}
			}
			if(container.nodeType == 3){
				var prevnodeobj = dijit.range.adjacentNoneTextNode(container);
				var prevnode = prevnodeobj[0];
				len = prevnodeobj[1];
				if(prevnode){
					atmrange.moveToElementText(prevnode);
					atmrange.collapse(false);
					//if contentEditable is not inherit, the above collapse won't make the end point
					//in the correctly position: it always has a -1 offset, so compensate it
					if(prevnode.contentEditable != 'inherit'){
						len++;
					}
				}else{
					atmrange.moveToElementText(container.parentNode);
					atmrange.collapse(true);

					// Correct internal cursor position
					// http://bugs.dojotoolkit.org/ticket/15578
					atmrange.move('character', 1);
					atmrange.move('character', -1);
				}

				offset += len;
				if(offset > 0){
					if(atmrange.move('character', offset) != offset){
						console.error('Error when moving!');
					}
				}
			}

			return atmrange;
		},
		decomposeTextRange: function(range){
			var tmpary = dijit.range.ie.getEndPoint(range);
			var startContainer = tmpary[0], startOffset = tmpary[1];
			var endContainer = tmpary[0], endOffset = tmpary[1];

			if(range.htmlText.length){
				if(range.htmlText == range.text){ //in the same text node
					endOffset = startOffset + range.text.length;
				}else{
					tmpary = dijit.range.ie.getEndPoint(range, true);
					endContainer = tmpary[0],endOffset = tmpary[1];
//					if(startContainer.tagName == "BODY"){
//						startContainer = startContainer.firstChild;
//					}
				}
			}
			return [startContainer, startOffset, endContainer, endOffset];
		},
		setRange: function(range, startContainer, startOffset, endContainer, endOffset, collapsed){
			var start = dijit.range.ie.setEndPoint(range, startContainer, startOffset);

			range.setEndPoint('StartToStart', start);
			if(!collapsed){
				var end = dijit.range.ie.setEndPoint(range, endContainer, endOffset);
			}
			range.setEndPoint('EndToEnd', end || start);

			return range;
		}
	};

declare("dijit.range.W3CRange",null, {
	constructor: function(){
		if(arguments.length>0){
			this.setStart(arguments[0][0],arguments[0][1]);
			this.setEnd(arguments[0][2],arguments[0][3]);
		}else{
			this.commonAncestorContainer = null;
			this.startContainer = null;
			this.startOffset = 0;
			this.endContainer = null;
			this.endOffset = 0;
			this.collapsed = true;
		}
	},
	_updateInternal: function(){
		if(this.startContainer !== this.endContainer){
			this.commonAncestorContainer = dijit.range.getCommonAncestor(this.startContainer, this.endContainer);
		}else{
			this.commonAncestorContainer = this.startContainer;
		}
		this.collapsed = (this.startContainer === this.endContainer) && (this.startOffset == this.endOffset);
	},
	setStart: function(node, offset){
		offset=parseInt(offset);
		if(this.startContainer === node && this.startOffset == offset){
			return;
		}
		delete this._cachedBookmark;

		this.startContainer = node;
		this.startOffset = offset;
		if(!this.endContainer){
			this.setEnd(node, offset);
		}else{
			this._updateInternal();
		}
	},
	setEnd: function(node, offset){
		offset=parseInt(offset);
		if(this.endContainer === node && this.endOffset == offset){
			return;
		}
		delete this._cachedBookmark;

		this.endContainer = node;
		this.endOffset = offset;
		if(!this.startContainer){
			this.setStart(node, offset);
		}else{
			this._updateInternal();
		}
	},
	setStartAfter: function(node, offset){
		this._setPoint('setStart', node, offset, 1);
	},
	setStartBefore: function(node, offset){
		this._setPoint('setStart', node, offset, 0);
	},
	setEndAfter: function(node, offset){
		this._setPoint('setEnd', node, offset, 1);
	},
	setEndBefore: function(node, offset){
		this._setPoint('setEnd', node, offset, 0);
	},
	_setPoint: function(what, node, offset, ext){
		var index = dijit.range.getIndex(node, node.parentNode).o;
		this[what](node.parentNode, index.pop()+ext);
	},
	_getIERange: function(){
		var r = (this._body || this.endContainer.ownerDocument.body).createTextRange();
		dijit.range.ie.setRange(r, this.startContainer, this.startOffset, this.endContainer, this.endOffset, this.collapsed);
		return r;
	},
	getBookmark: function(){
		this._getIERange();
		return this._cachedBookmark;
	},
	_select: function(){
		var r = this._getIERange();
		r.select();
	},
	deleteContents: function(){
		var s = this.startContainer, r = this._getIERange();
		if(s.nodeType === 3 && !this.startOffset){
			//if the range starts at the beginning of a
			//text node, move it to before the textnode
			//to make sure the range is still valid
			//after deleteContents() finishes
			this.setStartBefore(s);
		}
		r.pasteHTML('');
		this.endContainer = this.startContainer;
		this.endOffset = this.startOffset;
		this.collapsed = true;
	},
	cloneRange: function(){
		var r = new dijit.range.W3CRange([this.startContainer,this.startOffset,
			this.endContainer,this.endOffset]);
		r._body = this._body;
		return r;
	},
	detach: function(){
		this._body = null;
		this.commonAncestorContainer = null;
		this.startContainer = null;
		this.startOffset = 0;
		this.endContainer = null;
		this.endOffset = 0;
		this.collapsed = true;
}
});
} //if(!window.getSelection)


return dijit.range;
});

},
'dijit/_editor/html':function(){
define("dijit/_editor/html", [
	"dojo/_base/array",
	"dojo/_base/lang", // lang.setObject
	"dojo/sniff" // has("ie")
], function(array, lang, has){

// module:
//		dijit/_editor/html
// summary:
//		HTML serialization utility functions used by editor

var exports = {};
lang.setObject("dijit._editor.html", exports);

var escape = exports.escapeXml = function(/*String*/ str, /*Boolean?*/ noSingleQuotes){
	// summary:
	//		Adds escape sequences for special characters in XML: `&<>"'`.
	//		Optionally skips escapes for single quotes.
	str = str.replace(/&/gm, "&amp;").replace(/</gm, "&lt;").replace(/>/gm, "&gt;").replace(/"/gm, "&quot;");
	if(!noSingleQuotes){
		str = str.replace(/'/gm, "&#39;");
	}
	return str; // string
};


exports.getNodeHtml = function(/*DomNode*/ node){
	// summary:
	//		Return string representing HTML for node and it's children
	var output = [];
	exports.getNodeHtmlHelper(node, output);
	return output.join("");
};

exports.getNodeHtmlHelper = function(/*DomNode*/ node, /*String[]*/ output){
	// summary:
	//		Pushes array of strings into output[] which represent HTML for node and it's children
	switch(node.nodeType){
		case 1: // element node
			var lName = node.nodeName.toLowerCase();
			if(!lName || lName.charAt(0) == "/"){
				// IE does some strange things with malformed HTML input, like
				// treating a close tag </span> without an open tag <span>, as
				// a new tag with tagName of /span.  Corrupts output HTML, remove
				// them.  Other browsers don't prefix tags that way, so will
				// never show up.
				return "";
			}
			output.push('<', lName);

			// store the list of attributes and sort it to have the
			// attributes appear in the dictionary order
			var attrarray = [], attrhash = {};
			var attr;
			if(has("dom-attributes-explicit") || has("dom-attributes-specified-flag")){
				// IE8+ and all other browsers.
				var i = 0;
				while((attr = node.attributes[i++])){
					// ignore all attributes starting with _dj which are
					// internal temporary attributes used by the editor
					var n = attr.name;
					if(n.substr(0,3) !== '_dj' &&
						(!has("dom-attributes-specified-flag") || attr.specified) &&
						!(n in attrhash)){	// workaround repeated attributes bug in IE8 (LinkDialog test)
						var v = attr.value;
						if(n == 'src' || n == 'href'){
							if(node.getAttribute('_djrealurl')){
								v = node.getAttribute('_djrealurl');
							}
						}
						if(has("ie") === 8 && n === "style"){
							v = v.replace("HEIGHT:", "height:").replace("WIDTH:", "width:");
						}
						attrarray.push([n,v]);
						attrhash[n] = v;
					}
				}
			}else{
				// IE6-7 code path
				var clone = /^input$|^img$/i.test(node.nodeName) ? node : node.cloneNode(false);
				var s = clone.outerHTML;
				// Split up and manage the attrs via regexp
				// similar to prettyPrint attr logic.
				var rgxp_attrsMatch = /[\w-]+=("[^"]*"|'[^']*'|\S*)/gi
				var attrSplit = s.match(rgxp_attrsMatch);
				s = s.substr(0, s.indexOf('>'));
				array.forEach(attrSplit, function(attr){
					if(attr){
						var idx = attr.indexOf("=");
						if(idx > 0){
							var key = attr.substring(0,idx);
							if(key.substr(0,3) != '_dj'){
								if(key == 'src' || key == 'href'){
									if(node.getAttribute('_djrealurl')){
										attrarray.push([key,node.getAttribute('_djrealurl')]);
										return;
									}
								}
								var val, match;
								switch(key){
									case 'style':
										val = node.style.cssText.toLowerCase();
										break;
									case 'class':
										val = node.className;
										break;
									case 'width':
										if(lName === "img"){
											// This somehow gets lost on IE for IMG tags and the like
											// and we have to find it in outerHTML, known IE oddity.
											match=/width=(\S+)/i.exec(s);
											if(match){
												val = match[1];
											}
											break;
										}
									case 'height':
										if(lName === "img"){
											// This somehow gets lost on IE for IMG tags and the like
											// and we have to find it in outerHTML, known IE oddity.
											match=/height=(\S+)/i.exec(s);
											if(match){
												val = match[1];
											}
											break;
										}
									default:
										val = node.getAttribute(key);
								}
								if(val != null){
									attrarray.push([key, val.toString()]);
								}
							}
						}
					}
				}, this);
			}
			attrarray.sort(function(a,b){
				return a[0] < b[0] ? -1 : (a[0] == b[0] ? 0 : 1);
			});
			var j = 0;
			while((attr = attrarray[j++])){
				output.push(' ', attr[0], '="',
					(typeof attr[1] === "string" ? escape(attr[1], true) : attr[1]), '"');
			}
			switch(lName){
				case 'br':
				case 'hr':
				case 'img':
				case 'input':
				case 'base':
				case 'meta':
				case 'area':
				case 'basefont':
					// These should all be singly closed
					output.push(' />');
					break;
				case 'script':
					// Browsers handle script tags differently in how you get content,
					// but innerHTML always seems to work, so insert its content that way
					// Yes, it's bad to allow script tags in the editor code, but some people
					// seem to want to do it, so we need to at least return them right.
					// other plugins/filters can strip them.
					output.push('>', node.innerHTML, '</', lName, '>');
					break;
				default:
					output.push('>');
					if(node.hasChildNodes()){
						exports.getChildrenHtmlHelper(node, output);
					}
					output.push('</', lName, '>');
			}
			break;
		case 4: // cdata
		case 3: // text
			// FIXME:
			output.push(escape(node.nodeValue, true));
			break;
		case 8: // comment
			// FIXME:
			output.push('<!--', escape(node.nodeValue, true), '-->');
			break;
		default:
			output.push("<!-- Element not recognized - Type: ", node.nodeType, " Name: ", node.nodeName, "-->");
	}
};

exports.getChildrenHtml = function(/*DomNode*/ node){
	// summary:
	//		Returns the html content of a DomNode's children
	var output = [];
	exports.getChildrenHtmlHelper(node, output);
	return output.join("");
};

exports.getChildrenHtmlHelper = function(/*DomNode*/ dom, /*String[]*/ output){
	// summary:
	//		Pushes the html content of a DomNode's children into out[]

	if(!dom){ return; }
	var nodes = dom["childNodes"] || dom;

	// IE issue.
	// If we have an actual node we can check parent relationships on for IE,
	// We should check, as IE sometimes builds invalid DOMS.  If no parent, we can't check
	// And should just process it and hope for the best.
	var checkParent = !has("ie") || nodes !== dom;

	var node, i = 0;
	while((node = nodes[i++])){
		// IE is broken.  DOMs are supposed to be a tree.  But in the case of malformed HTML, IE generates a graph
		// meaning one node ends up with multiple references (multiple parents).  This is totally wrong and invalid, but
		// such is what it is.  We have to keep track and check for this because otherwise the source output HTML will have dups.
		// No other browser generates a graph.  Leave it to IE to break a fundamental DOM rule.  So, we check the parent if we can
		// If we can't, nothing more we can do other than walk it.
		if(!checkParent || node.parentNode == dom){
			exports.getNodeHtmlHelper(node, output);
		}
	}
};

return exports;
});

},
'*now':function(r){r(['dojo/i18n!*preload*JBrowse/View/TrackList/nls/Faceted*["ar","ca","cs","da","de","el","en-gb","en-us","es-es","fi-fi","fr-fr","he-il","hu","it-it","ja-jp","ko-kr","nl-nl","nb","pl","pt-br","pt-pt","ru","sk","sl","sv","th","tr","zh-tw","zh-cn","ROOT"]']);}
}});
define(
    "JBrowse/View/TrackList/Faceted", [
        'dojo/_base/declare',
        'dojo/_base/array',
        'dijit/layout/AccordionContainer',
        'dijit/layout/ContentPane',
        'JBrowse/Util',
        'dojox/grid/EnhancedGrid',
        'dojox/grid/enhanced/plugins/IndirectSelection'
    ],
    function ( declare, dArray, AccordionContainer, ContentPane, Util, EnhancedGrid ){

var dojof = Util.dojof;
return declare( 'JBrowse.View.TrackList.Faceted', null,
   /**
    * @lends JBrowse.View.TrackList.Faceted.prototype
    */
   {

   /**
     * Track selector with facets and text searching.
     * @constructs
     */
   constructor: function(args) {
       this.browser = args.browser;
       this.tracksActive = {};
       this.config = args;

       // construct the discriminator for whether we will display a
       // facet selector for this facet
       this._isSelectableFacet = this._coerceFilter(
               args.selectableFacetFilter
               // default facet filtering function
               || function( facetName, store ){
                   return (
                       // has an avg bucket size > 1
                       store.getFacetStats( facetName ).avgBucketSize > 1
                    &&
                       // and not an ident or label attribute
                       ! dojo.some( store.getLabelAttributes()
                                    .concat( store.getIdentityAttributes() ),
                                    function(l) {return l == facetName;}
                                  )
                   );
               }
           );

       // construct a similar discriminator for which columns will be displayed
       this.displayColumns = args.displayColumns;
       this._isDisplayableColumn = this._coerceFilter(
           args.displayColumnFilter || function(l) { return l.toLowerCase() != 'label'; }
       );

       // data store that fetches and filters our track metadata
       this.trackDataStore = args.trackMetaData;

       // subscribe to commands coming from the the controller
       this.browser.subscribe( '/jbrowse/v1/c/tracks/show',
                       dojo.hitch( this, 'setTracksActive' ));
       // subscribe to commands coming from the the controller
       this.browser.subscribe( '/jbrowse/v1/c/tracks/hide',
                       dojo.hitch( this, 'setTracksInactive' ));
       this.browser.subscribe( '/jbrowse/v1/c/tracks/delete',
                       dojo.hitch( this, 'setTracksInactive' ));

       this.renderInitial();

       // once its data is loaded and ready
       this.trackDataStore.onReady( this, function() {

           // render our controls and so forth
           this.renderSelectors();

           // connect events so that when a grid row is selected or
           // deselected (with the checkbox), publish a message
           // indicating that the user wants that track turned on or
           // off
           dojo.connect( this.dataGrid.selection, 'onSelected', this, function(index) {
                         this._ifNotSuppressed( 'selectionEvents', function() {
                             this._suppress( 'gridUpdate', function() {
                                 this.browser.publish( '/jbrowse/v1/v/tracks/show', [this.dataGrid.getItem( index ).conf] );
                             });
                         });

           });
           dojo.connect( this.dataGrid.selection, 'onDeselected', this, function(index) {
                         this._ifNotSuppressed( 'selectionEvents', function() {
                             this._suppress( 'gridUpdate', function() {
                                 this.browser.publish( '/jbrowse/v1/v/tracks/hide', [this.dataGrid.getItem( index ).conf] );
                             });
                         });
           });
       });

       this.trackDataStore.onReady( this, '_updateFacetCounts' ); // just once at start
       this.trackDataStore.onReady( this, '_updateMatchCount' ); // just once at start

       dojo.connect( this.trackDataStore, 'onFetchSuccess', this, '_updateGridSelections' );
       dojo.connect( this.trackDataStore, 'onFetchSuccess', this, '_updateMatchCount' );

    },

    /**
     * Coerces a string or array of strings into a function that,
     * given a string, returns true if the string matches one of the
     * given strings.  If passed a function, just returns that
     * function.
     * @private
     */
    _coerceFilter: function( filter ) {
        // if we have a non-function filter, coerce to an array,
        // then convert that array to a function
        if( typeof filter == 'string' )
            filter = [filter];
        if( dojo.isArray( filter ) ) {
            filter = function( store, facetName) {
                return dojo.some( filter, function(fn) {
                                      return facetName == fn;
                                  });
            };
        }
        return filter;
    },

    /**
     * Call the given callback if none of the given event suppression flags are set.
     * @private
     */
    _ifNotSuppressed: function( suppressFlags, callback ) {
        if( typeof suppressFlags == 'string')
            suppressFlags = [suppressFlags];
        if( !this.suppress)
            this.suppress = {};
        if( dojo.some( suppressFlags, function(f) {return this.suppress[f];}, this) )
            return undefined;
        return callback.call(this);
    },

    /**
     * Call the given callback while setting the given event suppression flags.
     * @private
     */
    _suppress: function( suppressFlags, callback ) {
        if( typeof suppressFlags == 'string')
            suppressFlags = [suppressFlags];
        if( !this.suppress)
            this.suppress = {};
        dojo.forEach( suppressFlags, function(f) {this.suppress[f] = true; }, this);
        var retval = callback.call( this );
        dojo.forEach( suppressFlags, function(f) {this.suppress[f] = false;}, this);
        return retval;
    },

    /**
     * Call a method of our object such that it cannot call itself
     * by way of event cycles.
     * @private
     */
    _suppressRecursion: function( methodName ) {
        var flag   = ['method_'+methodName];
        var method = this[methodName];
        return this._ifNotSuppressed( flag, function() { this._suppress( flag, method );});
    },

    renderInitial: function() {
        this.containerElem = dojo.create( 'div', {
            id: 'faceted_tracksel',
            style: {
                left: '-95%',
                width: '95%',
                zIndex: 500
            }
        },
        document.body );

        // make the tab that turns the selector on and off
        dojo.create('div',
                    {
                        className: 'faceted_tracksel_on_off tab',
                        innerHTML: '<img src="img/left_arrow.png"><div>Select<br>tracks</div>'
                    },
                    this.containerElem
                   );
        this.mainContainer = new dijit.layout.BorderContainer(
            { design: 'headline', gutters: false },
            dojo.create('div',{ className: 'mainContainer' }, this.containerElem)
        );


        this.topPane = new dijit.layout.ContentPane(
            { region: 'top',
              id: "faceted_tracksel_top",
              content: '<div class="title">Select Tracks</div> '
                       + '<div class="topLink" style="cursor: help"><a title="Track selector help">Help</a></div>'
            });
        dojo.query('div.topLink a[title="Track selector help"]',this.topPane.domNode)
            .forEach(function(helplink){
                var helpdialog = new dijit.Dialog({
                    "class": 'help_dialog',
                    refocus: false,
                    draggable: false,
                    title: 'Track Selection',
                    content: '<div class="main">'
                             + '<p>The JBrowse Faceted Track Selector makes it easy to search through'
                             + ' large numbers of available tracks to find exactly the ones you want.'
                             + ' You can incrementally filter the track display to narrow it down to'
                             + ' those your are interested in.  There are two types of filtering available,'
                             + ' which can be used together:'
                             + ' <b>filtering with data fields</b>, and free-form <b>filtering with text</b>.'
                             + '</p>'
                             + '  <dl><dt>Filtering with Data Fields</dt>'
                             + '  <dd>The left column of the display contains the available <b>data fields</b>.  Click on the data field name to expand it, and then select one or more values for that field.  This narrows the search to display only tracks that have one of those values for that field.  You can do this for any number of fields.<dd>'
                             + '  <dt>Filtering with Text</dt>'
                             + '  <dd>Type text in the "Contains text" box to filter for tracks whose data contains that text.  If you type multiple words, tracks are filtered such that they must contain all of those words, in any order.  Placing "quotation marks" around the text filters for tracks that contain that phrase exactly.  All text matching is case insensitive.</dd>'
                             + '  <dt>Activating Tracks</dt>'
                             + "  <dd>To activate and deactivate a track, click its check-box in the left-most column.  When the box contains a check mark, the track is activated.  You can also turn whole groups of tracks on and off using the check-box in the table heading.</dd>"
                             + "  </dl>"
                             + "</div>"
                 });
                dojo.connect( helplink, 'onclick', this, function(evt) {helpdialog.show(); return false;});
            },this);

        this.mainContainer.addChild( this.topPane );

        // make both buttons toggle this track selector
        dojo.query( '.faceted_tracksel_on_off' )
            .onclick( dojo.hitch( this, 'toggle' ));

        this.centerPane = new dijit.layout.BorderContainer({region: 'center', "class": 'gridPane', gutters: false});
        this.mainContainer.addChild( this.centerPane );
        var textFilterContainer = this.renderTextFilter();

        this.busyIndicator = dojo.create(
            'div', {
                innerHTML: '<img src="img/spinner.gif">',
                className: 'busy_indicator'
            }, this.containerElem );

        this.centerPane.addChild(
            new dijit.layout.ContentPane(
                { region: 'top',
                  "class": 'gridControls',
                  content: [
                      dojo.create( 'button', {
                                       className: 'faceted_tracksel_on_off',
                                       innerHTML: '<img src="img/left_arrow.png"> <div>Back to browser</div>',
                                       onclick: dojo.hitch( this, 'hide' )
                                   }
                                 ),
                      dojo.create( 'button', {
                                       className: 'clear_filters',
                                       innerHTML:'<img src="img/red_x.png">'
                                                 + '<div>Clear All Filters</div>',
                                       onclick: dojo.hitch( this, function(evt) {
                                           this._clearTextFilterControl();
                                           this._clearAllFacetControls();
                                           this._async( function() {
                                               this.updateQuery();
                                               this._updateFacetCounts();
                                           },this).call();
                                       })
                                   }
                                 ),
                      this.busyIndicator,
                      textFilterContainer,
                      dojo.create('div', { className: 'matching_record_count' })
                  ]
                }
            )
        );


    },
    renderSelectors: function() {

        // make our main components
        var facetContainer = this.renderFacetSelectors();
        // put them in their places in the overall layout of the track selector
        facetContainer.set('region','left');
        this.mainContainer.addChild( facetContainer );

        this.dataGrid = this.renderGrid();
        this.dataGrid.set('region','center');
        this.centerPane.addChild( this.dataGrid );

        this.mainContainer.startup();
    },

    /** do something in a timeout to avoid blocking the UI */
    _async: function( func, scope ) {
        var that = this;
        return function() {
            var args = arguments;
            var nativeScope = this;
            that._busy( true );
            window.setTimeout(
                function() {
                    func.apply( scope || nativeScope, args );
                    that._busy( false );
                },
                50
            );
        };
    },

    _busy: function( busy ) {
        this.busyCount = Math.max( 0, (this.busyCount || 0) + ( busy ? 1 : -1 ) );
        if( this.busyCount > 0 )
            dojo.addClass( this.containerElem, 'busy' );
        else
            dojo.removeClass( this.containerElem, 'busy' );
    },

    renderGrid: function() {

        var displayColumns = dojo.filter(
            this.displayColumns || this.trackDataStore.getFacetNames(),
            dojo.hitch(this, '_isDisplayableColumn')
        );
        var colWidth = 90/displayColumns.length;

        var grid = new EnhancedGrid({
               id: 'trackSelectGrid',
               store: this.trackDataStore,
               selectable: true,
               escapeHTMLInData: ('escapeHTMLInData' in this.config) ? this.config.escapeHTMLInData : false,
               noDataMessage: "No tracks match the filtering criteria.",
               structure: [
                   dojo.map(
                       displayColumns,
                       function(facetName) {
                           // rename name to key to avoid configuration confusion
                           facetName = {name: 'key'}[facetName.toLowerCase()] || facetName;
                           return {'name': this._facetDisplayName(facetName), 'field': facetName.toLowerCase(), 'width': colWidth+'%'};
                       },
                       this
                   )
               ],
               plugins: {
                   indirectSelection: {
                       headerSelector: true
                   }
               }
           }
        );

        this._monkeyPatchGrid( grid );
        return grid;
    },

    /**
     * Given a raw facet name, format it for user-facing display.
     * @private
     */
    _facetDisplayName: function( facetName ) {
        // make renameFacets if needed
        this.renameFacets = this.renameFacets || function(){
            var rename = this.config.renameFacets || {};
            rename.key = rename.key || 'Name';
            return rename;
        }.call(this);

        return this.renameFacets[facetName] || Util.ucFirst( facetName.replace('_',' ') );
    },

    /**
     * Apply several run-time patches to the dojox.grid.EnhancedGrid
     * code to fix bugs and customize the behavior in ways that aren't
     * quite possible using the regular Dojo APIs.
     * @private
     */
    _monkeyPatchGrid: function( grid ) {

        // 1. monkey-patch the grid's onRowClick handler to not do
        // anything.  without this, clicking on a row selects it, and
        // deselects everything else, which is quite undesirable.
        grid.onRowClick = function() {};

        // 2. monkey-patch the grid's range-selector to refuse to select
        // if the selection is too big
        var origSelectRange = grid.selection.selectRange;
        grid.selection.selectRange = function( inFrom, inTo ) {
            var selectionLimit = 30;
            if( inTo - inFrom > selectionLimit ) {
                alert( "Too many tracks selected, please select fewer than "+selectionLimit+" tracks." );
                return undefined;
            }
            return origSelectRange.apply( this, arguments );
        };
    },

    renderTextFilter: function( parent ) {
        // make the text input for text filtering
        this.textFilterLabel = dojo.create(
            'label',
            { className: 'textFilterControl',
              innerHTML: 'Contains text ',
              id: 'tracklist_textfilter',
              style: {position: 'relative'}
            },
            parent
        );
        this.textFilterInput = dojo.create(
            'input',
            { type: 'text',
              size: 40,
              disabled: true, // disabled until shown
              onkeypress: dojo.hitch( this, function(evt) {
                  // don't pay attention to modifier keys
                  if( evt.keyCode == dojo.keys.SHIFT || evt.keyCode == dojo.keys.CTRL || evt.keyCode == dojo.keys.ALT )
                      return;

                  // use a timeout to avoid updating the display too fast
                  if( this.textFilterTimeout )
                      window.clearTimeout( this.textFilterTimeout );
                  this.textFilterTimeout = window.setTimeout(
                      dojo.hitch( this, function() {
                                      // do a new search and update the display
                                      this._updateTextFilterControl();
                                      this._async( function() {
                                          this.updateQuery();
                                          this._updateFacetCounts();
                                          this.textFilterInput.focus();
                                      },this).call();
                                      this.textFilterInput.focus();
                                  }),
                      500
                  );
                  this._updateTextFilterControl();

                  evt.stopPropagation();
              })
            },
            this.textFilterLabel
        );
        // make a "clear" button for the text filtering input
        this.textFilterClearButton = dojo.create('img', {
            src: 'img/red_x.png',
            className: 'text_filter_clear',
            onclick: dojo.hitch( this, function() {
                this._clearTextFilterControl();
                this._async( function() {
                    this.updateQuery();
                    this._updateFacetCounts();
                },this).call();
            }),
            style: {
                position: 'absolute',
                right: '4px',
                top: '20%'
            }
        }, this.textFilterLabel );

        return this.textFilterLabel;
    },

   /**
    * Clear the text filter control input.
    * @private
    */
    _clearTextFilterControl: function() {
        this.textFilterInput.value = '';
        this._updateTextFilterControl();
    },
    /**
     * Update the display of the text filter control based on whether
     * it has any text in it.
     * @private
     */
    _updateTextFilterControl: function() {
        if( this.textFilterInput.value.length )
            dojo.addClass( this.textFilterLabel, 'selected' );
        else
            dojo.removeClass( this.textFilterLabel, 'selected' );

    },

    /**
     * Create selection boxes for each searchable facet.
     */
    renderFacetSelectors: function() {
        var container = new AccordionContainer({style: 'width: 200px'});

        var store = this.trackDataStore;
        this.facetSelectors = {};

        // render a facet selector for a pseudo-facet holding
        // attributes regarding the tracks the user has been working
        // with
        var usageFacet = this._renderFacetSelector(
            'My Tracks', ['Currently Active', 'Recently Used'] );
        usageFacet.set('class', 'myTracks' );
        container.addChild( usageFacet );

        // for the facets from the store, only render facet selectors
        // for ones that are not identity attributes, and have an
        // average bucket size greater than 1
        var selectableFacets =
            dojo.filter( this.config.selectableFacets || store.getFacetNames(),
                         function( facetName ) {
                             return this._isSelectableFacet( facetName, this.trackDataStore );
                         },
                         this
                       );

        dojo.forEach( selectableFacets, function(facetName) {

            // get the values of this facet
            var values = store.getFacetValues(facetName).sort();
            if( !values || !values.length )
                return;

            var facetPane = this._renderFacetSelector( facetName, values );
            container.addChild( facetPane );
        },this);

        return container;
    },

    /**
     * Make HTML elements for a single facet selector.
     * @private
     * @returns {dijit.layout.ContentPane}
     */
    _renderFacetSelector: function( /**String*/ facetName, /**Array[String]*/ values ) {

        var facetPane = new ContentPane(
            {
                title: '<span id="facet_title_' + facetName +'" '
                    + 'class="facetTitle">'
                    + this._facetDisplayName(facetName)
                    + ' <a class="clearFacet"><img src="img/red_x.png" /></a>'
                    + '</span>'
            });

        // make a selection control for the values of this facet
        var facetControl = dojo.create( 'table', {className: 'facetSelect'}, facetPane.containerNode );
        // populate selector's options
        this.facetSelectors[facetName] = dojo.map(
            values,
            function(val) {
                var that = this;
                var node = dojo.create(
                    'tr',
                    { className: 'facetValue',
                      innerHTML: '<td class="count"></td><td class="value">'+ val + '</td>',
                      onclick: function(evt) {
                          dojo.toggleClass(this, 'selected');
                          that._updateFacetControl( facetName );
                          that._async( function() {
                              that.updateQuery();
                              that._updateFacetCounts( facetName );
                          }).call();
                      }
                    },
                    facetControl
                );
                node.facetValue = val;
                return node;
            },
            this
        );

        return facetPane;
    },

    /**
     * Clear all the selections from all of the facet controls.
     * @private
     */
    _clearAllFacetControls: function() {
       dojo.forEach( dojof.keys( this.facetSelectors ), function( facetName ) {
           this._clearFacetControl( facetName );
       },this);
    },

    /**
     * Clear all the selections from the facet control with the given name.
     * @private
     */
    _clearFacetControl: function( facetName ) {
        dojo.forEach( this.facetSelectors[facetName] || [], function(selector) {
                          dojo.removeClass(selector,'selected');
                      },this);
        this._updateFacetControl( facetName );
    },

    /**
     * Incrementally update the facet counts as facet values are selected.
     * @private
     */
    _updateFacetCounts: function( /**String*/ skipFacetName ) {
        dojo.forEach( dojof.keys( this.facetSelectors ), function( facetName ) {
            if( facetName == 'My Tracks' )// || facetName == skipFacetName )
                return;
            var thisFacetCounts = this.trackDataStore.getFacetCounts( facetName );
            dojo.forEach( this.facetSelectors[facetName] || [], function( selectorNode ) {
                dojo.query('.count',selectorNode)
                    .forEach( function(countNode) {
                         var count = thisFacetCounts ? thisFacetCounts[ selectorNode.facetValue ] || 0 : 0;
                         countNode.innerHTML = Util.addCommas( count );
                         if( count )
                             dojo.removeClass( selectorNode, 'disabled');
                         else
                             dojo.addClass( selectorNode, 'disabled' );
                     },this);
                //dojo.removeClass(selector,'selected');
            },this);
            this._updateFacetControl( facetName );
       },this);
    },

    /**
     * Update the title bar of the given facet control to reflect
     * whether it has selected values in it.
     */
    _updateFacetControl: function( facetName ) {
        var titleContent = dojo.byId('facet_title_'+facetName);

        // if all our values are disabled, add 'disabled' to our
        // title's CSS classes
        if( dArray.every( this.facetSelectors[facetName] ||[], function(sel) {
                            return dojo.hasClass( sel, 'disabled' );
                        },this)
          ) {
                dojo.addClass( titleContent, 'disabled' );
        }

        // if we have some selected values, make a "clear" button, and
        // add 'selected' to our title's CSS classes
        if( dojo.some( this.facetSelectors[facetName] || [], function(sel) {
                return dojo.hasClass( sel, 'selected' );
            }, this ) ) {
                var clearFunc = dojo.hitch( this, function(evt) {
                    this._clearFacetControl( facetName );
                    this._async( function() {
                        this.updateQuery();
                        this._updateFacetCounts( facetName );
                    },this).call();
                    evt.stopPropagation();
                });
                dojo.addClass( titleContent.parentNode.parentNode, 'activeFacet' );
                dojo.query( '> a', titleContent )
                    .forEach(function(node) { node.onclick = clearFunc; },this)
                    .attr('title','clear selections');
        }
        // otherwise, no selected values
        else {
                dojo.removeClass( titleContent.parentNode.parentNode, 'activeFacet' );
                dojo.query( '> a', titleContent )
                    .onclick( function(){return false;})
                    .removeAttr('title');
        }
    },

    /**
     * Update the query we are using with the track metadata store
     * based on the values of the search form elements.
     */
    updateQuery: function() {
        this._suppressRecursion( '_updateQuery' );
    },
    _updateQuery: function() {
        var newQuery = {};

        var is_selected = function(node) {
            return dojo.hasClass(node,'selected');
        };

        // update from the My Tracks pseudofacet
        (function() {
             var mytracks_options = this.facetSelectors['My Tracks'];

             // index the optoins by name
             var byname = {};
             dojo.forEach( mytracks_options, function(opt){ byname[opt.facetValue] = opt;});

             // if filtering for active tracks, add the labels for the
             // currently selected tracks to the query
             if( is_selected( byname['Currently Active'] ) ) {
                 var activeTrackLabels = dojof.keys(this.tracksActive || {});
                 newQuery.label = Util.uniq(
                     (newQuery.label ||[])
                     .concat( activeTrackLabels )
                 );
             }

             // if filtering for recently used tracks, add the labels of recently used tracks
             if( is_selected( byname['Recently Used'])) {
                 var recentlyUsed = dojo.map(
                     this.browser.getRecentlyUsedTracks(),
                     function(t){
                         return t.label;
                     }
                 );

                 newQuery.label = Util.uniq(
                     (newQuery.label ||[])
                     .concat(recentlyUsed)
                 );
             }

             // finally, if something is selected in here, but we have
             // not come up with any track labels, then insert a dummy
             // track label value that will never match, because the
             // query engine ignores empty arrayrefs.
             if( ( ! newQuery.label || ! newQuery.label.length )
                 && dojo.some( mytracks_options, is_selected )
               ) {
                   newQuery.label = ['FAKE LABEL THAT IS HIGHLY UNLIKELY TO EVER MATCH ANYTHING'];
             }

        }).call(this);

        // update from the text filter
        if( this.textFilterInput.value.length ) {
            newQuery.text = this.textFilterInput.value;
        }

        // update from the data-based facet selectors
        dojo.forEach( this.trackDataStore.getFacetNames(), function(facetName) {
            var options = this.facetSelectors[facetName];
            if( !options ) return;

            var selectedFacets = dojo.map(
                dojo.filter( options, is_selected ),
                function(opt) {return opt.facetValue;}
            );
            if( selectedFacets.length )
                newQuery[facetName] = selectedFacets;
        },this);

        this.query = newQuery;
        this.dataGrid.setQuery( this.query );
        this._updateMatchCount();
    },

    /**
     * Update the match-count text in the grid controls bar based
     * on the last query that was run against the store.
     * @private
     */
    _updateMatchCount: function() {
        var count = this.dataGrid.store.getCount();
        dojo.query( '.matching_record_count', this.containerElem )
            .forEach( function(n) {
                          n.innerHTML =
                              Util.addCommas(count)
                              + ' '+( dojof.keys(this.query||{}).length ? 'matching ' : '' )
                              +'track' + ( count == 1 ? '' : 's' );
                      },
                      this
                    );
    },

    /**
     * Update the grid to have only rows checked that correspond to
     * tracks that are currently active.
     * @private
     */
    _updateGridSelections: function() {
        // keep selection events from firing while we mess with the
        // grid
        this._ifNotSuppressed('gridUpdate', function(){
            this._suppress('selectionEvents', function() {
                this.dataGrid.selection.deselectAll();

                // check the boxes that should be checked, based on our
                // internal memory of what tracks should be on.
                for( var i= 0; i < Math.min( this.dataGrid.get('rowCount'), this.dataGrid.get('rowsPerPage') ); i++ ) {
                    var item = this.dataGrid.getItem( i );
                    if( item ) {
                        var label = this.dataGrid.store.getIdentity( item );
                        if( this.tracksActive[label] )
                            this.dataGrid.rowSelectCell.toggleRow( i, true );
                    }
                }

            });
        });
    },

    /**
     * Given an array of track configs, update the track list to show
     * that they are turned on.
     */
    setTracksActive: function( /**Array[Object]*/ trackConfigs ) {
        dojo.forEach( trackConfigs, function(conf) {
            this.tracksActive[conf.label] = true;
        },this);
        this._updateGridSelections();
    },

    /**
     * Given an array of track configs, update the track list to show
     * that they are turned off.
     */
    setTracksInactive: function( /**Array[Object]*/ trackConfigs ) {
        dojo.forEach( trackConfigs, function(conf) {
            delete this.tracksActive[conf.label];
        },this);
        this._updateGridSelections();
    },

    /**
     * Make the track selector visible.
     */
    show: function() {
        window.setTimeout( dojo.hitch( this, function() {
            this.textFilterInput.disabled = false;
            this.textFilterInput.focus();
        }), 300);

        dojo.addClass( this.containerElem, 'active' );
        dojo.animateProperty({
            node: this.containerElem,
            properties: {
                left: { start: -95, end: 0, units: '%' }
            }
        }).play();

        this.shown = true;
    },

    /**
     * Make the track selector invisible.
     */
    hide: function() {

        dojo.removeClass( this.containerElem, 'active' );

        dojo.animateProperty({
            node: this.containerElem,
            properties: {
                left: { start: 0, end: -95, units: '%' }
            }
        }).play();

        this.textFilterInput.blur();
        this.textFilterInput.disabled = true;

        this.shown = false;
    },

    /**
     * Toggle whether the track selector is visible.
     */
    toggle: function() {
        this.shown ? this.hide() : this.show();
    }
});
});