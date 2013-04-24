define([], function(){
	// A basic binding to value, our base class
	function Binding(value){
		this.value= value;
		if(value){
			value._binding = this;
		}
	}
	Binding.prototype = {
		then: function(callback){
			// get the value of this binding, notifying the callback of changes
			var callbacks= this.callbacks;
			(callbacks || (this.callbacks = [])).push(callback);
			if(callbacks){
				if("value" in this){
					callback(this.value);
				}
			}else{
				var self = this;
				callbacks = this.callbacks;
				this.getValue(function(value){
					self.value = value;
					for(var i = 0; i < callbacks.length; i++){
						callbacks[i](value);
					}
				});
			}
		},
		getValue: function(callback){
			if(this.hasOwnProperty("value")){
				callback(this.value);
			}
		},
		get: function(key, callback){
			// use an existing child/property if it exists
			var value, child = this['_' + key] || 
				(this['_' + key] = this.hasOwnProperty("value") && typeof this.value == "object" ?
					(value = this.value[key]) && typeof value != "object" ? new PropertyBinding(this.value, key) :
						convertToBindable(value) : new Binding());
			if(callback){
				return child.then(callback);
			}
			return child;
		},
		put: function(value){
			if(this.source){
				this.source.put(value);
			}
			this.is(value);
		},
		is: function(value){
			if(value !== this.value){
				this.value = value;
				if(typeof this.value == "object"){
					for(var i in value){
						if(i.charAt(0) != '_'){
							this.get(i).is(value[i]);
						}
					}
				}
				if(this.callbacks){
					for(var i = 0; i < this.callbacks.length; i++){
						this.callbacks[i](value);
					}
				}
			}
		},
		keys: function(callback){
			if(this.source){
				this.source.keys(callback);
			}
			for(var i in this.value){
				if(i.charAt(0) != '_'){
					callback(i, this.get(i));
				}			}
		},
		to: function(source, property){
			source = convertToBindable(source);
			if(property){
				source = source.get(property);
			}
			var self = this;
			this.source = source;
			source.then(function(value){
				self.is(value);
			});
			for(var i in this){
				if(i.charAt(0) == '_'){
					i = i.slice(1);
					var child = self.get(i);
					var sourceChild = source.get(i);
					if(child != sourceChild){
						child.to(sourceChild);
					}
				}
			}
			source.keys(function(i, sourceChild){
				if(!(('_' + i) in self)){
					var child = self.get(i);
					if(child != sourceChild){
						child.to(sourceChild);
					}
				}
			});
			return this;
		}
	};
	// StatefulBinding is used for binding to Stateful objects, particularly Dijit widgets
	function StatefulBinding(stateful){
		this.stateful = stateful;
		stateful._binding = this;
	}
	StatefulBinding.prototype = new Binding({}); 
	StatefulBinding.prototype.to = function(source){
		Binding.prototype.to.apply(this, arguments);
		source = this.source;
		var stateful = this.stateful;
		source.then(function(value){
			stateful.set('value', value);
		});
		stateful.watch('value', function(property, oldValue, value){
			if(oldValue !== value){
				source.put(value);
			}
		});
		return this;
	};
	StatefulBinding.prototype.get = function(key, callback){
		return this['_' + key] || (this['_' + key] = new StatefulPropertyBinding(this.stateful, key));
	};
	StatefulBinding.prototype.keys = function(callback){
		for(var i in this.stateful){
			if(i.charAt(0) != '_'){
				callback(i, this.get(i));
			}
		}
	};

	function StatefulPropertyBinding(stateful, name){
		this.stateful = stateful;
		this.name = name;
	}
	StatefulPropertyBinding.prototype = new Binding;
	StatefulPropertyBinding.prototype.getValue = function(callback){
		// get the value of this property
		var stateful = this.stateful,
			name = this.name;
		// get the current value
		callback(stateful.get(name));
		// watch for changes
		stateful.watch(name, function(name, oldValue, newValue){
			if(oldValue !== newValue){
				callback(newValue);
			}
		});
	};
	StatefulPropertyBinding.prototype.is = function(value){
		// don't go through setters, it is bubbling up through the source
		this.stateful._changeAttrValue(this.name, value);
		//return Binding.prototype.is.call(this, value);
	}
	StatefulPropertyBinding.prototype.put = function(value){
		// put a value, go through setter
		this.stateful.set(this.name, value);
	}

	function ElementBinding(element, container){
		this.element= element;
		this.container = container;
		element._binding = this;
	}
	function ContainerBinding(element){
		return new ElementBinding(element, true);
	}
	ElementBinding.prototype = new Binding({});
	var checkable = {radio: 1, checkbox: 1};
	ElementBinding.prototype.is = function(value){
		var element = this.element;
		if(this.container || element.tagName == "FORM"){
			return Binding.prototype.is.call(this, value);
		};
		if("value" in element){
			if(element.type == "radio"){
				element.checked = element.value == value;
			}else if(element.type == "checkbox"){
				element.checked = value;
			}else{
				element.value = value || "";
			}
		}else{
			element.innerText = value || "";
		}
		this.oldValue = value;
	};
	var inputLike = {
		"INPUT":1,
		"SELECT":1,
		"TEXTAREA":1
	};
	ElementBinding.prototype.to = function(source){
		Binding.prototype.to.apply(this, arguments);
		source = this.source;
		var element = this.element;
		if(element.tagName == "FORM"){
			var binding = this;
			function findInputs(tag){
				var inputs = element.getElementsByTagName(tag);
				for(var i = 0; i < inputs.length; i++){
					var input = inputs[i];
					if(input.name){
						bind(input, binding.get(input.name));
					}
				}
			}
			findInputs("input");
			findInputs("select");
		}else if(element.tagName in inputLike){
			var value, binding = this,
				onchange = element.onchange = function(){
				if(element.type == "radio"){
					if(element.checked){
						value = element.value;
					}else{
						return; // not checked, don't do anything
					}
				}else{
					value = element.type == "checkbox" ? element.checked : element.value;
				}
				source.put(typeof binding.oldValue == "number" && !isNaN(value) ? +value : value);
			};
			if(element.getAttribute('data-continuous')){
				element.onkeyup = onchange;
			}
		}
		return this;
	}
	ElementBinding.prototype.then = function(callback){
		var element = this.element;
		if(this.container){
			return Binding.prototype.then.call(this, callback);
		}
		if("value" in element){
			callback(element.value);
			element.onchange = function(){
				callback(element.value);
			};
		}else{
			callback(element.innerText);
		}
	};
	function ArrayBinding(){
		Binding.apply(this, arguments);
	}
	ArrayBinding.prototype = new Binding({});
	ArrayBinding.prototype.getValue = function(callback){
		var currentValues = [],
			updates = 0;
			length = this.value.length;
		// watch all the items, and return a resulting array whenever an item is updated
		for(var i = 0; i < length; i++){
			(function(i, source){
				when(source, function(value){
					currentValues[i] = value;
					updates++;
					if(updates >= length){
						callback(currentValues);
					}
				});
			})(i, this.value[i]);
		}
	}
	function PropertyBinding(object, name){
		this.object = object;
		this.name = name;
	}
	PropertyBinding.prototype = new Binding;
	PropertyBinding.prototype.getValue = function(callback){
		callback(this.object[this.name]);
	};
	PropertyBinding.prototype.put = function(value){
		this.object[this.name] = value;
		this.is(value);
	}
	function FunctionBinding(func, reverseFunc){
		this.func = func;
		this.reverseFunc = reverseFunc;
	}
	FunctionBinding.prototype = {
		then: function(callback){
			if(callback){
				var func = this.func;
				return this.source.then(function(value){
					callback(value.slice ? func.apply(this, value) : func(value));
				});
			}
		},
		get: function(key){
			return this[key] || (this[key] = new Binding('None'));
		},
		put: function(value){
			this.source.put(this.reverseFunc(value));
		},
		is: function(){},
		to: function(source){
			source = bind.apply(this, arguments);
			this.source = source;
			return this;
		},
		keys: function(){}
	} 
	function convertToBindable(object){
		return object ?
			object._binding || 
				(object.get ? 
					object.is ?
				 		object :
					 	new StatefulBinding(object) :
					object.nodeType ?
						new ElementBinding(object) :
						typeof object == "function" ?
							new FunctionBinding(object) :
							object instanceof Array ?
								new ArrayBinding(object) :
				 				new Binding(object))
			: new Binding(object);
	}
	function bind(to){
		// first convert target object to the bindable interface
		to = convertToBindable(to);
	
		for(var i = 1; i < arguments.length; i++){
			var arg = arguments[i];
			if(typeof arg == "object" || typeof arg == "function"){
				arg = convertToBindable(arguments[i]);
				to.to(arg);
			}else{
				to = to.get(arg);
			}
		}
		return to;
	};
	var nativeWatch = Object.prototype.watch;
	function get(object, key, callback){
		if(key.call){
			// get(object ,callback) form
			if(object.get){
				if(object.get.binding){
					object.get(key);
				}else{
					key(object.get('value'));
				}
			}else{
				key(object);
			}
		}else{
			var value;
			if(object.get){
				if(object.get.binding){
					object.get(key, callback);
				}else{
					key(object.get('value'));
					if(object.watch === nativeWatch){
						object.watch('value', key);
					}
				}
			}else{
				return new PropertyBinding(object, key);
			}
		}
	};
	bind.get = get;
	bind.Element = ElementBinding;
	bind.Container = ContainerBinding;
	bind.Binding = Binding;


	function when(value, callback){
		if(value && value.then){
			return value.then(callback);
		}
		return callback(value);
	}
	bind.when = when;

	return bind;
});