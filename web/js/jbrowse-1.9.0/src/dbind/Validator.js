define(['json-schema/validate', './bind'], function(validate, bind){
	function ValidatorBinding(schema){
		this.schema = schema;
	}
	ValidatorBinding.prototype = new bind.Binding;
	ValidatorBinding.prototype.put = function(value){
		var results = validate(value, this.schema);
		if(results.valid){
			this.source.put(value);
			this.get("error").is("");
		}else{
			var errors = [];
			for(var i = 0; i < results.errors.length; i++){
				errors.push(results.errors[i].property + results.errors[i].message);
			}
			this.get("error").is(errors);
		}
	};
	return ValidatorBinding;
});