define(['dbind/bind', 'dbind/Validator', 'put-selector/put'], function(bind, Validator, put){
	var get = bind.get;
	myObject = {quantity: 3, price: 5};
	return function(form){
		// TODO: put this in a model module
		var quantity = bind(
				new Validator({type:"number", maximum: 20, minimum: 10})).to 
					(myObject, 'quantity');

		var quantityRow = put(form, 'div');

		quantity.get("title").is("Quantity");

		function ValidationTextBox(){
			var mainElement = put('div');
			var binding = new bind.Container(mainElement);
			// the label
			bind(put(mainElement, 'label')).to(binding, 'title');
			// the main value is bound to the input
			bind(put(mainElement, 'input[type=text]')).to(binding);
			// any errors go after it
			bind(put(mainElement, 'span.error-message')).to(binding, 'error');
			return mainElement;
		}
		var quantityTextBox = ValidationTextBox();
		put(quantityRow, quantityTextBox);
		bind(quantityTextBox).to(quantity);
		put(form, "div", "Price", "input[type=text]", {name: "price"});
		bind(form, myObject);
		
/*		bind(put(quantityRow, 'input[type=text][placeholder=number]'), quantity);
		bind(put(quantityRow, 'span.error-message'), quantity.get('error'));
		bind(put(form, 'div label'), quantity.get("title"))
		bind(put(form, 'div span'), quantity);*/
		bind(put(form, 'div label', 'Total Price: ', '< span'), bind(function(quantity, price){
			return "$" + quantity * price;
		}).to([quantity, bind(myObject, "price")]));
	}
});