dbind is a functional reactive data binding package that provides straightforward binding of data to
components like form inputs, validation connectors, and more. The dbind framework 
is designed to help you create organized, well-structured, layered applications, facilitating
a clean separation between a data model with validation logic and presentation elements.
It is also intended to be compatible with Dojo and bindr, giving you the full capabilities of the
bindr reactive data binding language with Dojo and Dijit widgets. 

Note that not all of the features described here are implemented and/or tested, this
is a work in progress.

# Getting Started

The foundational module of dbind is `dbind/bind`, which returns a function that gives you
bindable objects. Usage is very simple, call bind with a target component and than indicate what
object or property you want to bind to:

	require(['dbind/bind'], function(bind){
		bind(anInputElement).to(myObject, "propertyName");
	});

And just like that we have a two-way binding. The value from the object's property will be provided to the
input. Any changes to the input will cause the object's property to be modified. We could even
bind another element to this property to easily see the binding in action.

	bind(myDiv).to(myObject, 'propertyName');

And the value of the property would be put in the div, and updated any time the input
was changed.

We can also create a property binding that can encapsulate a single property, and 
be directly bound to components:

	var myProperty = bind(myObject, 'propertyName');
	// now we can bind components to this property
	bind(anInputElement).to(myProperty);

We can bind object properties to inputs, where the property value is synchronized with
the input's value, and container elements like divs, where the property value is outputted
to the element's inner text.

In addition we can also bind an object to a form. In this case, dbind will search through
the form inputs and bind each one to the object's properties based on the input's "name" attribute.
If you have a form that you wish to create with HTML, this makes it very easy to bind an object to it 
without directly referencing each input: 

	bind(myForm).to(myObject);

## Dijit Components

With our bindings, we can easily use direct DOM elements or Dijit components interchangeably.
For example, we could bind a Dijit TextBox to myProperty as well:

	require(['dijit/form/TextBox', 'dbind/bind'], function(TextBox){
		var textBox = new TextBox({}, 'textbox');
		bind(textBox).to(myProperty);

## Transformations

We can also bind to functions to create a transformation for our binding. The function
will be called with a value to convert, allowing us to continuously apply transformation
to a source object. For example, we could create a functional transformation:

	function double(x){
		return x * 2;
	}
	var doubledValue = bind(double).to(sourceValue);

Now doubledValue will contain the value equal to twice the value of sourceValue. This
will remain true even as sourceValue varies in the future, doubleValue will continually
stay in sync.

We can also bind a transformation function to multiple source objects:

	function multiply(x, y){
		return x * y;
	}
	var productValue = bind(multiply).to([sourceValueA, sourceValueB]);

# Validation

With dbind, we can do more than bind data to elements, we can also bind simple
data objects to validation layers to compose more sophisticated data models, that can
then be bound to UI elements. To bind to a validator, first we create a validator, giving
it a validation definition (based on JSON Schema), and then we bind it to a property or
object: 

	require(['dbind/bind', 'dbind/Validator'], function(bind, Validator){
		// create a validator and bind it to a property of myObject 
		var myProperty = bind(new Validator({type:"number", maximum: 20, minimum: 10})).
			to(myObject, 'propertyName');
		// now we can bind the validated property to an element
		bind(anInputElement).to(myProperty); 

And now when a user enters a value that is not a number or doesn't fall in the given
range it will be rejected.

The validator also gives us access to the error message so the UI can properly display
information to the user on why the input is invalid:

	bind(errorMessageElement).to(myProperty, 'error');

Any time an error occurs in validation, the `errorMessageElement` will automatically
be updated with error message. This makes it easy to build coherent, manageable 
validated forms. The validation layer is distinct from the UI layer, and they can easily 
be wired together for responsive validated forms and UIs.

# dbind Interfaces

dbind relies on several interfaces for connecting components. You can interact with these objects
using the following API, or you can create your own implementations of the APIs. The bind() function 
returns bindable objects. Bindable objects have the follow method:  

* to(source, property?) - This binds this object to the provided source object. Any changes in the
source object will be propagated to the bindable target object.  If a property argument
is provided, the target object will be bound to the property of the source object. We speak 
of changes coming
from the source object as traveling *up* to the target. The target component may be
UI component that can support editing, sending requested changes from the user *down*
to the source.

* is(value) - This is called to change the value of the target object from a downstream
source.

The source object should be a reactive object. The bind() function generally returns
objects that are also reactive. A reactive object represents a value that may change
over time. The reactive object has a value at any given point in time, and may change
to different values over time. It has the following methods:

* then(callback) - This is called to get the value of the source object, both now and in the 
future. A function should be provided, and will be called with the current value of the
source object, and called again each time it changes in the future. It is worth noting
that a source object that is constant (value doesn't change) is the same as a promise,
and can be provided to consumers that expect a promise.
 
* put(value) - This is called to change the value of the source object from an upstream
target component. This may be rejected.
This method may be omitted if the source object can't be modified by upstream components.

Source objects may also be mappable; they can have properties. This is provided through
the following methods:

* get(property) - Returns a reactive object for the given property.

* get(property, callback) - Shorthand for get(property).then(callback).

* set(property, value) - Shorthand for get(property).put(value).

# Composition of Binding-Driven Components

TODO

