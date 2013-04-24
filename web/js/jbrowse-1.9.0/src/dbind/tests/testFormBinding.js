define(['dbind/bind', 'dbind/Validator', 'put-selector/put'], function(bind, Validator, put){
	var get = bind.get;
	myObject = {quantity: 3, price: 5, discounted: true, color: "red", pattern: "solid"};
	return function(form){
		// TODO: put this in a model module
		var quantity = bind(
				new Validator({type:"number", maximum: 20, minimum: 10})).to 
					(get(myObject, 'quantity'));
		
		
		quantity.get("title").is("Quantity");
		
		function ValidationTextBox(){
			var mainElement = put('div');
			var binding = new bind.Container(mainElement);
			// the label
			bind(put(mainElement, 'label')).to(binding.get('title'));
			// the main value is bound to the input
			bind(put(mainElement, 'input[type=text]')).to(binding);
			// any errors go after it
			bind(put(mainElement, 'span.error-message')).to(binding.get('error'));
			return mainElement;
		}
		// create the form elements
		var quantityRow = put(form, 'div');
		var quantityTextBox = ValidationTextBox();
		put(quantityRow, quantityTextBox);
		bind(quantityTextBox).to(quantity);
		
		put(form, "div", "Price", "input[type=text]", {name: "price"});
		
		put(form, "div", "Discounted", "input[type=checkbox]", {name: "discounted"});
		
		var patternSelect = put(form, "div", "Pattern", "select[name=pattern]");
		put(patternSelect, "option[value=[striped]", "Striped");
		put(patternSelect, "option[value=[solid]", "Solid");
		
		var colorDiv = put(form, "div", "Color");
		put(colorDiv, "div", "Red", "input[type=radio][value=red]", {name: "color"});
		put(colorDiv, "div", "Green", "input[type=radio][value=green]", {name: "color"});
		put(colorDiv, "div", "Blue", "input[type=radio][value=blue]", {name: "color"});

		// now bind the form to the object, this will bind the elements to the properties of the object
		bind(form, myObject);
		
		bind(put(form, 'div label', 'Total Price: ', '< span'), bind(function(quantity, price, discounted){
			return "$" + quantity * price * (discounted ? 0.9 : 1);
		}).to([quantity, bind(myObject, "price"), bind(myObject, "discounted")]));
	}
});