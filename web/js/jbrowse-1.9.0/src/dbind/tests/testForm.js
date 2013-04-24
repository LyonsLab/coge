define(['dbind/bind', 'dbind/Validator', 'put-selector/put'], function(bind, Validator, put){
	var get = bind.get;
	myObject = {quantity: 3, price: 5, discounted: true, color: "red", pattern: "striped"};
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
			bind(put(mainElement, 'span.error-message')).to(binding, 'error');
			return mainElement;
		}
		// create the form elements
		var quantityRow = put(form, 'div');
		var quantityTextBox = ValidationTextBox();
		put(quantityRow, quantityTextBox);
		bind(quantityTextBox).to(quantity);
		
		bind(put(form, "div", "Price", "input[type=text]")).to(myObject, "price");
		
		bind(put(form, "div", "Discounted", "input[type=checkbox]")).to(myObject, "discounted");
		
		var patternSelect = put(form, "div", "Pattern", "select");
		put(patternSelect, "option[value=striped]", "Striped");
		put(patternSelect, "option[value=solid]", "Solid");
		bind(patternSelect).to(myObject, "pattern");
		
		var colorDiv = put(form, "div", "Color");
		var colorProperty = bind(myObject, "color");
		bind(put(colorDiv, "div", "Red", "input[type=radio][value=red]")).to(colorProperty);
		bind(put(colorDiv, "div", "Green", "input[type=radio][value=green]")).to(colorProperty);
		bind(put(colorDiv, "div", "Blue", "input[type=radio][value=blue]")).to(colorProperty);

		bind(put(form, 'div label', 'Total Price: ', '< span'), bind(function(quantity, price, discounted){
			return "$" + quantity * price * (discounted ? 0.9 : 1);
		}).to([quantity, bind(myObject, "price"), bind(myObject, "discounted")]));
	}
});